#!/usr/bin/env python3
import logging
import numpy as np
import copy

import crosstalk
import gates
import predistortion
import pulses
import qubits
import readout
import tomography

# Allow logging to Labber's instrument log
log = logging.getLogger('LabberDriver')

# TODO Reduce calc of CZ by finding all unique TwoQubitGates in seq and calc.
# TODO Make I(width=None) have the width of the longest gate in the step
# TODO Add checks so that not both t0 and dt are given
# TODO Two composite gates should be able to be parallel
# TODO check number of qubits in seq and in gate added to seq
# TODO Remove pulse from I gates


class GateOnQubit:
    def __init__(self, gate, qubit, pulse=None):
        self.gate = gate
        self.qubit = qubit
        self.pulse = pulse

        if pulse is None:
            self.duration = 0
        else:
            self.duration = pulse.total_duration()

    def __str__(self):
        return "Gate {} on qubit {}".format(self.gate, self.qubit)

    def __repr__(self):
        return self.__str__()


class Step:
    """Represent one step in a sequence.

    Parameters
    ----------
    n_qubit : int
        Number of qubits in the sequece.
    t0 : float
        Center of the sequence in seconds (the default is None).
    dt : float
        Spacing to previous pulse in seconds (the default is None). Use only
        either t0 or dt.
    align : str {'left', 'center', 'right'}
        The alignment of pulses if they have different lengths,
        (the default is 'center').

    Attributes
    ----------
    gates : list of :dict:
        The different gates in the step.

    """

    def __init__(self, t0=None, dt=None, align='center'):
        self.gates = []
        self.align = align
        self.t0 = t0
        self.dt = dt
        self.t_start = None
        self.t_end = None

    def add_gate(self, qubit, gate):
        """Add the given gate to the specified qubit(s).

        The number of gates must equal the number of qubits.

        If the number of qubits given are less than the number of qubits in the
        step, I gates are added to the other qubits.
        Parameters
        ----------
        qubit : int or list of int
            The qubit indices.
        gate : :obj:`BaseGate`
            The gate(s).

        """
        if gate.number_of_qubits() > 1 and not isinstance(qubit, list):
            raise ValueError(
                "Provide a list of qubits for gates with more than one qubit")

        if gate.number_of_qubits() > 1 and not gate.number_of_qubits() == len(
                qubit):
            raise ValueError(
                """Number of qubits in the gate must equal the number of qubit
                indices given""")

        if gate.number_of_qubits() == 1 and not isinstance(qubit, int):
            raise ValueError("Provide qubit as int for gates with one qubit")

        if isinstance(qubit, int):
            if self._qubit_in_step(qubit):
                raise ValueError("Qubit {} already in step.".format(qubit))
        else:
            for n in qubit:
                if self._qubit_in_step(n):
                    raise ValueError("Qubit {} already in step.".format(n))

        self.gates.append(GateOnQubit(gate, qubit))

    def time_shift(self, shift):
        """Shift the timings of the step.

        Parameters
        ----------
        shift : float
            The amount of shift to apply in seconds.

        """
        self.t_start += shift
        self.t0 += shift
        self.t_end += shift

    def _qubit_in_step(self, qubit):
        """Returns whatever the given qubit is in the step or not. """
        if not isinstance(qubit, int):
            raise ValueError("Qubit index should be int.")

        def _in(input_list, n):
            flat_list = []
            for sublist_or_el in input_list:
                if isinstance(sublist_or_el, list):
                    if _in(sublist_or_el, n) is True:
                        return True
                elif sublist_or_el == n:
                    return True
            return False

        return _in([x.qubit for x in self.gates], qubit)

    def __str__(self):
        return str(self.gates)

    def __repr__(self):
        return str(self.gates)


class Sequence:
    """A multi qubit seqence.

    Parameters
    ----------
    n_qubit : type
        The number of qubits in the sequence.

    Attributes
    ----------
    sequences : list of :obj:`Step`
        Holds the steps of the sequence.
    perform_process_tomography : bool
        Flag for performing process tomography.
    perform_state_tomography : bool
        Flag for performing state tomography.
    readout_delay : float
        Delay time between last pulse and readout, in seconds.
    n_qubit

    """

    def __init__(self, n_qubit):
        self.n_qubit = n_qubit
                
        # log.info('initiating empty seqence list')
        self.sequence_list = []

        # perform Heralding
        self.perform_heralding = False
        self.heralding_delay = 0.0
        
        self.perform_heralding_2 = False
        self.heralding_delay_2 = 0.0

        # process tomography
        self.perform_process_tomography = False
        self._process_tomography = tomography.ProcessTomography()

        # state tomography
        self.perform_state_tomography = False
        self._state_tomography = tomography.StateTomography()

        # readout
        self.readout_delay = 0.0

    # Public methods
    def generate_sequence(self, config):
        """Generate sequence by adding gates/pulses to waveforms.

        Parameters
        ----------
        config : dict
            Configuration as defined by Labber driver configuration window.

        """
        raise NotImplementedError()

    def get_sequence(self, config):
        """Compile sequence and return it.

        Parameters
        ----------
        config : dict
            Labber instrument configuration.

        Returns
        -------
        list of :obj:`Step`
            The compiled qubit sequence.

        """
        self.sequence_list = []
		
        if self.perform_initialization:
            self.add_gate([0],[gates.Zini], align='right')			
            if self.initialization_delay > 0:
                delay_initialization = gates.IdentityGate(width=self.initialization_delay)
                self.add_gate_to_all(delay_initialization, dt=0, align='right')
        
        if self.perform_heralding:
            self.add_gate_to_all(gates.HeraldingGate(), dt=0, align='right')			
            if self.heralding_delay > 0:
                delay_heralding = gates.IdentityGate(width=self.heralding_delay)
                self.add_gate_to_all(delay_heralding, dt=0, align='right')        
			
        if self.perform_process_tomography:
            self._process_tomography.add_pulses(self)

        self.generate_sequence(config)

        if self.perform_state_tomography:
            self._state_tomography.add_pulses(self)
            
        if self.perform_heralding_2:
            self.add_gate_to_all(gates.HeraldingGate2(), dt=0, align='left')	
            if self.heralding_delay_2 > 0:
                delay_heralding_2 = gates.IdentityGate(width=self.heralding_delay_2)
                self.add_gate_to_all(delay_heralding_2, dt=0)
                
        if self.readout_delay > 0:
            delay = gates.IdentityGate(width=self.readout_delay)
            self.add_gate_to_all(delay, dt=0)
        self.add_gate_to_all(gates.ReadoutGate(), dt=0, align='left')

        return self

    # Public methods for adding pulses and gates to the sequence.
    def add_single_pulse(self,
                         qubit,
                         pulse,
                         t0=None,
                         dt=None,
                         align_left=False):
        """Add single qubit pulse to specified qubit.

        This function still exist to not break existing
        funcationallity. You should really use the add_gate method.

        t0 or dt can be used to override the global pulse spacing.

        Parameters
        ----------
        qubit : int
            Qubit number, indexed from 0.

        pulse : :obj:`Pulse`
            Definition of pulse to add.

        t0 : float, optional
            Absolute pulse position.

        dt : float, optional
            Pulse spacing, referenced to the previous pulse.

        align_left: bool, optional
            If True, aligns the pulse to the left. Defaults to False.

        """
        gate = gates.CustomGate(pulse)
        if align_left is True:
            self.add_gate(qubit, gate, t0, dt, 'left')
        else:
            self.add_gate(qubit, gate, t0, dt, 'center')

    def add_single_gate(self, qubit, gate, t0=None, dt=None, align_left=False):
        """Add single gate to specified qubit sequence.

        Note, this function still exist is to not break existing
        funcationallity. You should really use the add_gate method.

        t0 or dt can be used to override the global pulse spacing.

        Parameters
        ----------
        qubit : int
            Qubit number, indexed from 0.

        gate : :obj:`Gate`
            Definition of gate to add.

        t0 : float, optional
            Absolute pulse position.

        dt : float, optional
            Pulse spacing, referenced to the previous pulse.

        align_left : boolean, optional
            If True, t0 is the start of the pulse, otherwise it is the center
            of the pulse. False is the default.

        """
        if align_left is True:
            self.add_gate(qubit, gate, t0, dt, 'left')
        else:
            self.add_gate(qubit, gate, t0, dt, 'center')

    def add_gate(self,
                 qubit,
                 gate,
                 t0=None,
                 dt=None,
                 align='center',
                 index=None):
        """Add a set of gates to the given qubit sequences.

        For the qubits with no specificied gate, an IdentityGate will be given.
        The length of the step is given by the longest pulse.

        Parameters
        ----------
        qubit : int or list of int
            The qubit(s) to add the gate(s) to.
        gate : :obj:`BaseGate` or list of :obj:`BaseGate`
            The gate(s) to add.
        t0 : float, optional
            Absolute gate position (the default is None).
        dt : float, optional
            Gate spacing, referenced to the previous pulse
            (the default is None).
        align : str, optional
            If two or more qubits have differnt pulse lengths, `align`
            specifies how those pulses should be aligned. 'Left' aligns the
            start, 'center' aligns the centers, and 'right' aligns the end,
            (the default is 'center').
        index : int, optional
            Where in the sequence to insert the new gate.

        """
        step = Step(t0=t0, dt=dt, align=align)
        if isinstance(gate, list):
            if not isinstance(qubit, list):
                raise ValueError(
                    """Provide qubit indices as a list when adding more than
                    one gate.""")
            if len(gate) != len(qubit):
                raise ValueError(
                    "Length of gate list must equal length of qubit list.")

            for q, g in zip(qubit, gate):
                step.add_gate(q, g)
        else:
            if gate.number_of_qubits() > 1:
                if not isinstance(qubit, list):
                    raise ValueError(
                        "Provide qubit list for gates with more than one qubit"
                    )
            else:
                if not isinstance(qubit, int):
                    raise ValueError(
                        "For single gates, give qubit as int (not list).")
            step.add_gate(qubit, gate)
            # log.info('adding gate {} to {}. 2qb gate: {}'.format(gate, qubit, isinstance(gate, gates.TwoQubitGate)))

        if index is None:
            self.sequence_list.append(step)
            # log.info('adding step to sequence list')
            # log.info('sequence len is {}'.format(len(self.sequence_list)))
        else:
            self.sequence_list.insert(index + 1, step)
            # log.info('inserting step in sequence list')
            # log.info('sequence len is {}'.format(len(self.sequence_list)))


    def add_gate_to_all(self, gate, t0=None, dt=None, align='center'):
        """Add a single gate to all qubits.

        Pulses are added at the end of the sequence, with the gate spacing set
        by either the spacing parameter or the aboslut position.
        """
        if isinstance(gate, list):
            raise ValueError("Only single gates allowed.")
        if isinstance(gate, (gates.BaseGate, gates.CompositeGate)):
            if gate.number_of_qubits() > 1:
                raise ValueError(
                    "Not clear how to add multi-qubit gates to all qubits.")

        qubit = list(range((self.n_qubit)))
        gate = [gate for n in range(self.n_qubit)]
        self.add_gate(qubit, gate, t0=t0, dt=dt, align=align)

    def add_gates(self, gates):
        """Add multiple gates to the qubit waveform.

        Pulses are added at the end of the sequence, with the gate spacing set
        by the spacing parameter.

        Examples
        --------
        Add three gates to a two-qubit sequence, first a positive pi-pulse
        around X to qubit 1, then a negative pi/2-pulse to qubit 2, finally
        simultaneous positive pi-pulses to qubits 1 and 2.

        >>> add_gates([[gates.Xp,  None    ],
                       [None,     gates.Y2m],
                       [gates.Xp,  gates.Xp]])

        Parameters
        ----------
        gates : list of list of :obj:`BaseGate`
            List of lists defining gates to add. The innermost list should
            have the same length as number of qubits in the sequence.

        """
        # make sure we have correct input
        if not isinstance(gates, (list, tuple)):
            raise Exception('The input must be a list of list with gates')
        if len(gates) == 0:
            return
        if not isinstance(gates[0], (list, tuple)):
            raise Exception('The input must be a list of list with gates')
        # add gates sequence to waveforms
        for gate in gates:
            # add gate to specific qubit waveform
            qubit = list(range(len(gate)))
            self.add_gate(qubit, gate)

    def set_parameters(self, config={}):
        """Set base parameters using config from from Labber driver.

        Parameters
        ----------
        config : dict
            Configuration as defined by Labber driver configuration window

        """
        # sequence parameters
        d = dict(
            Zero=0,
            One=1,
            Two=2,
            Three=3,
            Four=4,
            Five=5,
            Six=6,
            Seven=7,
            Eight=8,
            Nine=9)
        # If the number of qubits changed, we need to re-init
        if self.n_qubit != d[config.get('Number of qubits')]:
            self.__init__(d[config.get('Number of qubits')])

        # Readout
        self.readout_delay = config.get('Readout delay')

        #Heralding
        self.heralding_delay = config.get('Delay after heralding')
        self.heralding_delay_2 = config.get('Delay after heralding #2')

        #Initialization
        self.initialization_delay = config.get('Delay after initialization')

        # perform Initialization
        self.perform_initialization = \
            config.get('Activate Initialization', False)

        # perform heralding
        self.perform_heralding = \
            config.get('Activate Heralding', False)
        self.perform_heralding_2 = \
            config.get('Activate Heralding #2', False)

        # process tomography prepulses
        self.perform_process_tomography = \
            config.get('Generate process tomography prepulse', False)
        self._process_tomography.set_parameters(config)

        # state tomography
        self.perform_state_tomography = config.get(
            'Generate state tomography postpulse', False)
        self._state_tomography.set_parameters(config)


class SequenceToWaveforms:
    """Compile a multi qubit sequence into waveforms.

    Parameters
    ----------
    n_qubit : type
        The maximum number of qubits.

    Attributes
    ----------
    dt : float
        Pulse spacing, in seconds.
    local_xy : bool
        If False, collate all waveforms into one.
    simultaneous_pulses : bool
        If False, seperate all pulses in time.
    sample_rate : float
        AWG Sample rate.
    n_pts : float
        Number of points in the waveforms.
    first_delay : float
        Delay between start of waveform and start of the first pulse.
    trim_to_sequence : bool
        If True, adjust `n_points` to just fit the sequence.
    align_to_end : bool
        Align the whole sequence to the end of the waveforms.
        Only relevant if `trim_to_sequence` is False.
    sequences : list of :obj:`Step`
        The qubit sequences.
    qubits : list of :obj:`Qubit`
        Parameters of each qubit.
    wave_xy_delays : list of float
        Indiviudal delays for the XY waveforms.
    wave_z_delays : list of float
        Indiviudal delays for the Z waveforms.
    n_qubit

    """

    def __init__(self, n_qubit):
        self.n_qubit = n_qubit
        self.dt = 10E-9
        self.local_xy = True
        self.single_qubit_pulse_scaling = 'Amplitude'
        self.simultaneous_pulses = True

        # waveform parameter
        self.sample_rate = 1.2E9
        self.n_pts = 240E3
        self.first_delay = 100E-9
        self.trim_to_sequence = True
        self.align_to_end = False

        self.sequence_list = []
        self.qubits = [qubits.Qubit() for n in range(self.n_qubit)]

        # waveforms
        self._wave_xy = [
            np.zeros(0, dtype=np.complex) for n in range(self.n_qubit)
        ]
        # log.info('_wave_z initiated to 0s')
        self._wave_z = [np.zeros(0) for n in range(self.n_qubit)]
        self._wave_gate = [np.zeros(0) for n in range(self.n_qubit)]
        self._wave_gate_z = [np.zeros(0) for n in range(self.n_qubit)]

        # waveform delays
        self.wave_xy_delays = np.zeros(self.n_qubit)
        self.wave_z_delays = np.zeros(self.n_qubit)

        # define pulses
        self.pulses_1qb_xy = [None for n in range(self.n_qubit)]
        self.pulses_1qb_z = [None for n in range(self.n_qubit)]
        self.pulses_2qb = [None for n in range(self.n_qubit - 1)]
        self.pulses_readout = [None for n in range(self.n_qubit)]
        self.pulses_heralding = [None for n in range(self.n_qubit)]
        self.pulses_heralding_2 = [None for n in range(self.n_qubit)]

        # cross-talk
        self.compensate_crosstalk = False
        self._crosstalk = crosstalk.Crosstalk()

        # predistortion
        self.perform_predistortion = False
        self._predistortions = [
            predistortion.Predistortion(n) for n in range(self.n_qubit)
        ]
        self._predistortions_z = [
            predistortion.ExponentialPredistortion(n)
            for n in range(self.n_qubit)
        ]

        # gate switch waveform
        self.generate_gate_switch = False
        self.uniform_gate = False
        self.gate_delay = 0.0
        self.gate_overlap = 20E-9
        self.minimal_gate_time = 20E-9

        # filters
        self.use_gate_filter = False
        self.use_z_filter = False

        # readout trig settings
        self.readout_trig_generate = False
        self.heralding_trig_generate = False
        self.heralding_trig_generate_2 = False

        # readout wave object and settings
        self.readout = readout.Demodulation(self.n_qubit)
        self.readout_trig = np.array([], dtype=float)
        self.heralding_trig = np.array([], dtype=float)
        self.heralding_trig_2 = np.array([], dtype=float)
        self.readout_iq = np.array([], dtype=np.complex)
        self.heralding_iq = np.array([], dtype=np.complex)
        self.heralding_iq_2 = np.array([], dtype=np.complex)

    def get_waveforms(self, sequence):
        """Compile the given sequence into waveforms.

        Parameters
        ----------
        sequences : list of :obj:`Step`
            The qubit sequence to be compiled.

        Returns
        -------
        type
            Description of returned object.

        """
        self.sequence = sequence
        self.sequence_list = sequence.sequence_list
        # log.info('Start of get_waveforms. Len sequence list: {}'.format(len(self.sequence_list)))
        # log.info('Point 1: Sequence_list[3].gates = {}'.format(self.sequence_list[3].gates))

        if not self.simultaneous_pulses:
            self._seperate_gates()
        # log.info('Point 2: Sequence_list[6].gates = {}'.format(self.sequence_list[6].gates))

        self._explode_composite_gates()
        # log.info('Point 3: Sequence_list[6].gates = {}'.format(self.sequence_list[6].gates))

        if self.compose_singleqb_gate_virtualZ:
            self._compose_singleqb_gate_virtualZ_()

        self._add_pulses_and_durations()
        # log.info('Point 4: Sequence_list[6].gates = {}'.format(self.sequence_list[6].gates))

        self._add_timings()
        # log.info('Point 5: Sequence_list[6].gates = {}'.format(self.sequence_list[6].gates))

        self._init_waveforms()
        # log.info('Point 6: Sequence_list[6].gates = {}'.format(self.sequence_list[6].gates))


        if self.align_to_end:
            shift = self._round((self.n_pts - 2) / self.sample_rate -
                                self.sequence_list[-1].t_end)
            for step in self.sequence_list:
                step.time_shift(shift)

        self._perform_virtual_z()
        self._generate_waveforms()

        # collapse all xy pulses to one waveform if no local XY control or assigned them to the specified waveforms
        if not self.local_xy:
            if self.waveform_assignation == 'Multiplexed':
                # sum all waveforms to first one
                self._wave_xy[0] = np.sum(self._wave_xy[:self.n_qubit], 0)
                # clear other waveforms
                for n in range(1, self.n_qubit):
                    self._wave_xy[n][:] = 0.0
            elif self.waveform_assignation == 'Manual':
            
                self._wave_xy_temp = [np.zeros(0, dtype=np.complex) for n in range(self.n_qubit)]
                
                for n in range(self.n_qubit):
                
                    self._wave_xy_temp[n] = np.zeros(self.n_pts, dtype=np.complex)
                    
                # reAssign the waveform according to the dictionary
                for n in range(self.n_qubit):
                    self._wave_xy_temp[self.QubitWaveformDictionary[n]] = np.sum([self._wave_xy_temp[self.QubitWaveformDictionary[n]],self._wave_xy[n]], 0)
                    
                self._wave_xy = self._wave_xy_temp
        # log.info('before predistortion, _wave_z max is {}'.format(np.max(self._wave_z)))
        # if self.compensate_crosstalk:
        #     self._perform_crosstalk_compensation()
        if self.perform_predistortion:
            self._predistort_xy_waveforms()
        if self.perform_predistortion_z:
            self._predistort_z_waveforms()
        if self.readout_trig_generate:
            self._add_readout_trig()
        if self.heralding_trig_generate:
            self._add_heralding_trig()
        if self.heralding_trig_generate_2:
            self._add_heralding_trig(heralding_number=2)
        if self.generate_gate_switch:
            self._add_microwave_gate()
            self._add_microwave_gate_z()
        self._filter_output_waveforms()

        # Apply offsets
        self.readout_iq += self.readout_i_offset + 1j * self.readout_q_offset

        #log.info('self.readout_iq.shape is %s' % str(self.readout_iq.shape))
        if self.activate_cavity_leakage:
            t = np.linspace(0, (self.n_pts - 1) / self.sample_rate, self.n_pts)
            y = self.leakage_amplitude
            phase = self.leakage_phase
            # single-sideband mixing, get frequency
            omega = 2 * np.pi * self.leakage_frequency
            # apply SSBM transform
            data_i = (y.real * np.cos(omega * t - phase) +
                      -y.imag * np.cos(omega * t - phase + +np.pi / 2))
            data_q = (y.real * np.sin(omega * t - phase) +
                      -y.imag * np.sin(omega * t - phase + +np.pi / 2))
            self.readout_iq += data_i + 1j * data_q
        # create and return dictionary with waveforms
        waveforms = dict()
        waveforms['xy'] = self._wave_xy
        self._wave_z[0] += self.ovr_offset        
        waveforms['z'] = self._wave_z
        waveforms['gate'] = self._wave_gate
        waveforms['gate_z'] = self._wave_gate_z
        waveforms['readout_trig'] = np.sum([self.readout_trig, self.heralding_trig, self.heralding_trig_2],0)
        waveforms['readout_iq'] = np.sum([self.readout_iq, self.heralding_iq, self.heralding_iq_2], 0)
        
        # log.info('returning z waveforms in get_waveforms. Max is {}'.format(np.max(waveforms['z'])))
        return waveforms

    def _seperate_gates(self):
        new_sequences = []
        for step in self.sequence_list:
            if any(
                    isinstance(gate, (gates.ReadoutGate, gates.HeraldingGate, gates.HeraldingGate2, gates.IdentityGate))
                    for gate in step.gates):
                # Don't seperate I gates or readouts since we do
                # multiplexed readout
                new_sequences.append(step)
                continue
            for gate in step.gates:
                # log.info('In seperate gates, handling gate {}'.format(gate))
                if gate.gate is not None:
                    new_step = Step(
                        t0=step.t_start, dt=step.dt, align=step.align)
                    new_step.add_gate(gate.qubit, gate.gate)
                    new_sequences.append(new_step)
        # log.info('New sequence [6] is {}'.format(new_sequences[6].gates))
        self.sequence_list = new_sequences

    def _add_timings(self):
        t_start = 0
        for step in self.sequence_list:
            if step.dt is None and step.t0 is None:
                # Use global pulse spacing
                step.dt = self.dt
            # Find longest gate in sequence
            max_duration = np.max([x.duration for x in step.gates])
            if step.t0 is None:
                step.t_start = self._round(t_start + step.dt)
                step.t0 = self._round(step.t_start + max_duration / 2)
            else:
                step.t_start = self._round(step.t0 - max_duration / 2)
            step.t_end = self._round(step.t_start + max_duration)
            t_start = step.t_end  # Next step starts where this one ends
            # Avoid double spacing for steps with 0 duration
            if max_duration == 0:
                t_start = t_start - step.dt

        # Make sure that the sequence is sorted chronologically.
        # self.sequence_list.sort(key=lambda x: x.t_start) # TODO Fix this

        # Make sure that sequnce starts on first delay
        time_diff = self._round(self.first_delay -
                                self.sequence_list[0].t_start)
        for step in self.sequence_list:
            step.time_shift(time_diff)
        
    def _add_pulses_and_durations(self):
        for step in self.sequence_list:
            for gate in step.gates:
                if gate.pulse is None:
                    gate.pulse = self._get_pulse_for_gate(gate)                                  
                if gate.pulse is None:
                    gate.duration = 0
                else:
                    gate.duration = gate.pulse.total_duration()

    def _get_pulse_for_gate(self, gate):
        qubit = gate.qubit
        gate = gate.gate
        # Virtual Z is special since it has no length
        if isinstance(gate, gates.VirtualZGate):
            pulse = None
        # Get the corresponding pulse for other gates
        elif isinstance(gate, gates.SingleQubitXYRotation):
            # log.info('before the error:' + str(gate))
            # log.info(str(gate==gates.X2m))
            pulse = gate.get_adjusted_pulse(self.pulses_1qb_xy[qubit], \
                pulse_scaling = self.single_qubit_pulse_scaling, \
                scaling_factor_pi2 = self.single_qubit_pulse_scaling_factor[qubit], \
                detuning_pi2 = self.single_qubit_pulse_scaling_factor[qubit+2])
        elif isinstance(gate, gates.SingleQubitZRotation):
            pulse = gate.get_adjusted_pulse(self.pulses_1qb_z[qubit], \
                pulse_scaling = self.single_qubit_pulse_scaling)#, \
                #scaling_factor_pi2 = 0.5)
        elif isinstance(gate, gates.IdentityGate):
            if gate.width==None:
                pulse = None
            else:
                pulse = gate.get_adjusted_pulse(self.pulses_1qb_xy[qubit])
        elif isinstance(gate, gates.TwoQubitGate):
            pulse = gate.get_adjusted_pulse(self.pulses_2qb[qubit[0]])            
        elif isinstance(gate, gates.ReadoutGate):
            pulse = gate.get_adjusted_pulse(self.pulses_readout[qubit])
        elif isinstance(gate, gates.HeraldingGate):
            pulse = gate.get_adjusted_pulse(self.pulses_heralding[qubit])
        elif isinstance(gate, gates.HeraldingGate2):
            pulse = gate.get_adjusted_pulse(self.pulses_heralding_2[qubit])
        elif isinstance(gate, gates.CustomGate):
            pulse = gate.get_adjusted_pulse(gate.pulse)
        else:
            raise ValueError('Please provide a pulse for {}'.format(gate))

        return pulse

    def _predistort_xy_waveforms(self):
        """Pre-distort the waveforms."""
        # go through and predistort all xy waveforms

        n_wave = self.n_qubit if not(self.waveform == 'Multiplexed') else 1
        for n in range(n_wave):
            self._wave_xy[n] = self._predistortions[n].predistort(
                self._wave_xy[n])

    def _predistort_z_waveforms(self):
        # go through and predistort all waveforms
        for n in range(self.n_qubit):
            self._wave_z[n] = self._predistortions_z[n].predistort(
                self._wave_z[n])

    def _perform_crosstalk_compensation(self):
        """Compensate for Z-control crosstalk."""
        self._wave_z = self._crosstalk.compensate(self._wave_z)

    def _compose_singleqb_gate_virtualZ_(self):
        n = 0
        while n < len(self.sequence_list):
            step = self.sequence_list[n]
            #log.info(step)
            new_gate = []
            new_qubit = []
            for gate in step.gates:
                if isinstance(gate.gate, gates.SingleQubitXYRotation):
                    qubit_for_z_rotation, VZ_pi2, VZ_pi = self.compose_VirtualZ_dictionary[gate.qubit]
                    if abs(gate.gate.theta)==np.pi:
                        new_qubit.append(qubit_for_z_rotation-1)
                        new_gate.append(VZ_pi)
                    elif abs(gate.gate.theta)==np.pi/2:
                        new_qubit.append(qubit_for_z_rotation-1)
                        new_gate.append(VZ_pi2)
            for i in range(len(new_gate)):
                self.sequence.add_gate(new_qubit[i], new_gate[i], index = n + i)
            n += 1
            
    def _explode_composite_gates(self):
        # Loop through the sequence until all CompositeGates are removed
        # Note that there could be nested CompositeGates
        n = 0
        while n < len(self.sequence_list):
            step = self.sequence_list[n]
            i = 0
            while i < len(step.gates):
                gate = step.gates[i]
                if isinstance(gate.gate, gates.CompositeGate):
                    # # log.info('In exploded composite, handling composite gate {} at step {}'.format(gate, n))
                    for m, g in enumerate(gate.gate.sequence):
                        new_gate = [x.gate for x in g.gates]
                        # Single gates shouldn't be lists
                        if len(new_gate) == 1:
                            new_gate = new_gate[0]

                        # Translate gate qubit number to device qubit number
                        new_qubit = [x.qubit for x in g.gates]
                        for j, q in enumerate(new_qubit):
                            if isinstance(q, int):
                                if isinstance(gate.qubit, int):
                                    new_qubit[j] = gate.qubit
                                    continue
                                new_qubit[j] = gate.qubit[q]
                            else:
                                new_qubit[j] = []
                                for k in q:
                                    new_qubit[j].append(gate.qubit[k])

                        # Single qubit shouldn't be lists
                        if len(new_qubit) == 1:
                            new_qubit = new_qubit[0]
                        # # log.info('In explode composite; modifying {} by adding gate {} at index {}'.format(gate, new_gate, n+m))
                        self.sequence.add_gate(
                            new_qubit, new_gate, index=n + m)

                    del step.gates[i]
                    # # log.info('In composite gates, removing step {}', i)

                    continue
                i = i + 1
            n = n + 1

        # Remove any empty steps where the composite gates were
        i = 0
        while i < len(self.sequence_list):
            step = self.sequence_list[i]
            if len(step.gates) == 0:
                del self.sequence_list[i]
                # log.info('In composite gates, removing step {}', i)
                continue
            i = i + 1

        # for i, step in enumerate(self.sequence_list):
            # log.info('At end of explode, step {} is {}'.format(i, step.gates))

    def _perform_virtual_z(self):
        """Shifts the phase of pulses subsequent to virtual z gates."""
        for qubit in range(self.n_qubit):
            phase = 0
            for step in self.sequence_list:
                for gate in step.gates:
                    gate_obj = None
                    if qubit == gate.qubit or isinstance(gate.qubit, list):  # TODO Allow for 2 qb
                        gate_obj = gate.gate
                    if isinstance(gate_obj, gates.VirtualZGate):
                        phase -= gate_obj.theta
                        continue
                    if (isinstance(gate_obj, gates.SingleQubitXYRotation)
                            and phase != 0):
                        gate.gate = copy.copy(gate_obj)
                        gate.gate.phi += phase
                        # Need to recompute the pulse
                        gate.pulse = self._get_pulse_for_gate(gate)
                    #if (isinstance(gate_obj, gates.SingleQubitZRotation)
                    #        and phase != 0):
                    #    gate.gate = copy.copy(gate_obj)
                    #    gate.gate.phi += phase
                        # Need to recompute the pulse
                    #    gate.pulse = self._get_pulse_for_gate(gate)
                    if (isinstance(gate_obj, gates.TwoQubitGate) # TODO need to double check
                            and phase != 0):                        
                        if qubit == self.TargetQubit_VZ_2QB:    
                            gate.gate = copy.copy(gate_obj)                        
                            # gate.pulse.phase += phase
                            gate.gate.phi += phase
                            # log.info('gate.pulse.phase is ' + str(gate.pulse.phase))
                            # Need to recompute the pulse
                            gate.pulse = self._get_pulse_for_gate(gate)

    def _add_microwave_gate(self):
        """Create waveform for gating microwave switch."""
        n_wave = self.n_qubit if not(self.waveform_assignation == 'Multiplexed') else 1
        # go through all waveforms
        for n, wave in enumerate(self._wave_xy[:n_wave]):
            if self.uniform_gate:
                # the uniform gate is all ones
                gate = np.ones_like(wave)
                # if creating readout trig, turn off gate during readout
                if self.readout_trig_generate:
                    gate[-int((self.readout_trig_duration - self.gate_overlap -
                               self.gate_delay) * self.sample_rate):] = 0.0
            else:
                # non-uniform gate, find non-zero elements
                gate = np.array(np.abs(wave) > 0.0, dtype=float)
                # fix gate overlap
                n_overlap = int(np.round(self.gate_overlap * self.sample_rate))
                diff_gate = np.diff(gate)
                indx_up = np.nonzero(diff_gate > 0.0)[0]
                indx_down = np.nonzero(diff_gate < 0.0)[0]
                # add extra elements to left and right for overlap
                for indx in indx_up:
                    gate[max(0, indx - n_overlap):(indx + 1)] = 1.0
                for indx in indx_down:
                    gate[indx:(indx + n_overlap + 1)] = 1.0

                # fix gaps in gate shorter than min (look for 1>0)
                diff_gate = np.diff(gate)
                indx_up = np.nonzero(diff_gate > 0.0)[0]
                indx_down = np.nonzero(diff_gate < 0.0)[0]
                # ignore first transition if starting in zero
                if gate[0] == 0:
                    indx_up = indx_up[1:]
                n_down_up = min(len(indx_down), len(indx_up))
                len_down = indx_up[:n_down_up] - indx_down[:n_down_up]
                # find short gaps
                short_gaps = np.nonzero(
                    len_down < (self.minimal_gate_time * self.sample_rate))[0]
                for indx in short_gaps:
                    gate[indx_down[indx]:(1 + indx_up[indx])] = 1.0

                # shift gate in time
                n_shift = int(np.round(self.gate_delay * self.sample_rate))
                if n_shift < 0:
                    n_shift = abs(n_shift)
                    gate = np.r_[gate[n_shift:], np.zeros((n_shift, ))]
                elif n_shift > 0:
                    gate = np.r_[np.zeros((n_shift, )), gate[:(-n_shift)]]
            # make sure gate starts/ends in 0
            gate[0] = 0.0
            gate[-1] = 0.0
            # store results
            self._wave_gate[n] = gate
            
    def _add_microwave_gate_z(self):
        """Create waveform for gating microwave switch."""
        n_wave = self.n_qubit if not(self.waveform_assignation == 'Multiplexed') else 1
        # go through all waveforms
        for n, wave in enumerate(self._wave_z[:n_wave]):
            if self.uniform_gate:
                # the uniform gate is all ones
                gate = np.ones_like(wave)
                # if creating readout trig, turn off gate during readout
                if self.readout_trig_generate:
                    gate[-int((self.readout_trig_duration - self.gate_overlap -
                               self.gate_delay) * self.sample_rate):] = 0.0
            else:
                # non-uniform gate, find non-zero elements
                gate = np.array(np.abs(wave) > 0.0, dtype=float)
                # fix gate overlap
                n_overlap = int(np.round(self.gate_overlap * self.sample_rate))
                diff_gate = np.diff(gate)
                indx_up = np.nonzero(diff_gate > 0.0)[0]
                indx_down = np.nonzero(diff_gate < 0.0)[0]
                # add extra elements to left and right for overlap
                for indx in indx_up:
                    gate[max(0, indx - n_overlap):(indx + 1)] = 1.0
                for indx in indx_down:
                    gate[indx:(indx + n_overlap + 1)] = 1.0

                # fix gaps in gate shorter than min (look for 1>0)
                diff_gate = np.diff(gate)
                indx_up = np.nonzero(diff_gate > 0.0)[0]
                indx_down = np.nonzero(diff_gate < 0.0)[0]
                # ignore first transition if starting in zero
                if gate[0] == 0:
                    indx_up = indx_up[1:]
                n_down_up = min(len(indx_down), len(indx_up))
                len_down = indx_up[:n_down_up] - indx_down[:n_down_up]
                # find short gaps
                short_gaps = np.nonzero(
                    len_down < (self.minimal_gate_time * self.sample_rate))[0]
                for indx in short_gaps:
                    gate[indx_down[indx]:(1 + indx_up[indx])] = 1.0

                # shift gate in time
                n_shift = int(np.round(self.gate_delay * self.sample_rate))
                if n_shift < 0:
                    n_shift = abs(n_shift)
                    gate = np.r_[gate[n_shift:], np.zeros((n_shift, ))]
                elif n_shift > 0:
                    gate = np.r_[np.zeros((n_shift, )), gate[:(-n_shift)]]
            # make sure gate starts/ends in 0
            gate[0] = 0.0
            gate[-1] = 0.0
            # store results
            self._wave_gate_z[n] = gate

    def _filter_output_waveforms(self):
        """Filter output waveforms"""
        # start with gate
        if self.use_gate_filter and self.gate_filter_size > 1:
            # prepare filter
            window = self._get_filter_window(
                self.gate_filter_size, self.gate_filter,
                self.gate_filter_kaiser_beta)
            # apply filter to all output waveforms
            n_wave = self.n_qubit if not(self.waveform_assignation == 'Multiplexed') else 1
            for n in range(n_wave):
                self._wave_gate[n] = self._apply_window_filter(
                    self._wave_gate[n], window)
                self._wave_gate_z[n] = self._apply_window_filter(
                    self._wave_gate_z[n], window)
                # make sure gate starts/ends in 0
                self._wave_gate[n][0] = 0.0
                self._wave_gate[n][-1] = 0.0
                self._wave_gate_z[n][0] = 0.0
                self._wave_gate_z[n][-1] = 0.0

        # same for z waveforms
        if self.use_z_filter and self.z_filter_size > 1:
            # prepare filter
            window = self._get_filter_window(
                self.z_filter_size, self.z_filter, self.z_filter_kaiser_beta)
            # apply filter to all output waveforms
            for n in range(self.n_qubit):
                self._wave_z[n] = self._apply_window_filter(
                    self._wave_z[n], window)

    def _get_filter_window(self, size=11, window='Kaiser', kaiser_beta=14.0):
        """Get filter for waveform convolution"""
        if window == 'Rectangular':
            w = np.ones(size)
        elif window == 'Bartlett':
            # for filters that start/end in zeros, add 2 points and truncate
            w = np.bartlett(max(1, size+2))
            w = w[1:-1]
        elif window == 'Blackman':
            w = np.blackman(size + 2)
            w = w[1:-1]
        elif window == 'Hamming':
            w = np.hamming(size)
        elif window == 'Hanning':
            w = np.hanning(size + 2)
            w = w[1:-1]
        elif window == 'Kaiser':
            w = np.kaiser(size, kaiser_beta)
        else:
            raise('Unknown filter windows function %s.' % str(window))
        return w/w.sum()

    def _apply_window_filter(self, x, window):
        """Apply window filter to input waveform

        Parameters
        ----------
        x: np.array
            Input waveform.
        window: np.array
            Filter waveform.

        Returns
        -------
        np.array
            Filtered waveform.

        """
        # buffer waveform to avoid wrapping effects at boundaries
        n = len(window)
        s = np.r_[2*x[0] - x[n-1::-1], x, 2*x[-1] - x[-1:-n:-1]]
        # apply convolution
        y = np.convolve(s, window, mode='same')
        return y[n:-n+1]

    def _round(self, t, acc=1E-12):
        """Round the time `t` with a certain accuarcy `acc`.

        Parameters
        ----------
        t : float
            The time to be rounded.
        acc : float
            The accuarcy (the default is 1E-12).

        Returns
        -------
        float
            The rounded time.

        """
        return int(np.round(t / acc)) * acc

    def _add_readout_trig(self):
        """Create waveform for readout trigger."""
        trig = np.zeros_like(self.readout_iq)
        
        start = int(np.min(((np.abs(self.readout_iq) > 0.0).nonzero()[0][0] + self.readout_trig_delay * self.sample_rate,self.n_pts_readout)))
        end = int(
            np.min((start + self.readout_trig_duration * self.sample_rate,
                    self.n_pts_readout)))
        trig[start:end] = self.readout_trig_amplitude

        # make sure trig starts and ends in 0.
        trig[0] = 0.0
        trig[-1] = 0.0
        self.readout_trig = trig

    def _add_heralding_trig(self, heralding_number=1):
        """Create waveform for readout trigger."""
        trig = np.zeros_like(self.readout_iq)
        
        if heralding_number == 2:
            iq = self.heralding_iq_2
            delay = self.heralding_trig_delay_2
        else:
            iq = self.heralding_iq
            delay = self.heralding_trig_delay
        start = int(np.min(((np.abs(iq) > 0.0).nonzero()[0][0] + delay * self.sample_rate,self.n_pts_readout)))
        end = int(
            np.min((start + self.readout_trig_duration * self.sample_rate,
                    self.n_pts_readout)))
        trig[start:end] = self.readout_trig_amplitude

        # make sure trig starts and ends in 0.
        trig[0] = 0.0
        trig[-1] = 0.0
        if heralding_number == 2:
            self.heralding_trig_2 = trig
        else:
            self.heralding_trig = trig

    def _init_waveforms(self):
        """Initialize waveforms according to sequence settings."""
        # To keep the first pulse delay, use the smallest delay as reference.
        min_delay = np.min([
            self.wave_xy_delays[:self.n_qubit],
            self.wave_z_delays[:self.n_qubit]
        ])
        self.wave_xy_delays -= min_delay
        self.wave_z_delays -= min_delay
        max_delay = np.max([
            self.wave_xy_delays[:self.n_qubit],
            self.wave_z_delays[:self.n_qubit]
        ])

        # find the end of the sequence
        # only include readout in size estimate if all waveforms have same size
        if self.readout_match_main_size:
            end = np.max([s.t_end for s in self.sequence_list]) + max_delay
        else:
            end = np.max([s.t_end
                          for s in self.sequence_list[0:-1]]) + max_delay

        # create empty waveforms of the correct size
        if self.trim_to_sequence:
            self.n_pts = int(np.ceil(end * self.sample_rate)) + 1
            if self.n_pts % 2 == 1:
                # Odd n_pts give spectral leakage in FFT
                self.n_pts += 1
        for n in range(self.n_qubit):
            self._wave_xy[n] = np.zeros(self.n_pts, dtype=np.complex)
            # log.info('wave z {} initiated to 0'.format(n))
            self._wave_z[n] = np.zeros(self.n_pts, dtype=np.complex)
            self._wave_gate[n] = np.zeros(self.n_pts, dtype=float)
            self._wave_gate_z[n] = np.zeros(self.n_pts, dtype=float)

        # Waveform time vector
        self.t = np.arange(self.n_pts) / self.sample_rate

        # readout trig and i/q waveforms
        if self.readout_match_main_size:
            # same number of points for readout and main waveform
            self.n_pts_readout = self.n_pts
        else:
            # different number of points for readout and main waveform
            self.n_pts_readout = 1 + int(
                np.ceil(self.sample_rate * (self.sequence_list[-1].t_end -
                                            self.sequence_list[-1].t_start)))
            if self.n_pts_readout % 2 == 1:
                # Odd n_pts give spectral leakage in FFT
                self.n_pts_readout += 1

        self.readout_trig = np.zeros(self.n_pts_readout, dtype=float)
        self.heralding_trig = np.zeros(self.n_pts_readout, dtype=float)
        self.heralding_trig_2 = np.zeros(self.n_pts_readout, dtype=float)
        self.readout_iq = np.zeros(self.n_pts_readout, dtype=np.complex)
        self.heralding_iq = np.zeros(self.n_pts_readout, dtype=np.complex)
        self.heralding_iq_2 = np.zeros(self.n_pts_readout, dtype=np.complex)

    def _generate_waveforms(self):
        """Generate the waveforms corresponding to the sequence."""
        # log.info('generating waveform from sequence. Len is {}'.format(len(self.sequence_list)))
        for step in self.sequence_list:
            # log.info('Generating gates {}'.format(step.gates))
            for gate in step.gates:
                qubit = gate.qubit
                if isinstance(qubit, list):
                    qubit = qubit[0]
                gate_obj = gate.gate

                if isinstance(gate_obj,
                              (gates.IdentityGate, gates.VirtualZGate)):
                    continue
                elif isinstance(gate_obj, gates.SingleQubitZRotation):
                    waveform = self._wave_z[qubit]
                    delay = self.wave_z_delays[qubit]
                    if self.compensate_crosstalk:
                        crosstalk = self._crosstalk.compensation_matrix[:,
                                                                        qubit]
                elif isinstance(gate_obj, gates.TwoQubitGate):
                    # log.info('adding 2qb gate waveforms')
                    waveform = self._wave_z[qubit]
                    delay = self.wave_z_delays[qubit]
                    if self.compensate_crosstalk:
                        crosstalk = self._crosstalk.compensation_matrix[:,
                                                                        qubit]
                    
                elif isinstance(gate_obj, gates.SingleQubitXYRotation):
                    waveform = self._wave_xy[qubit]
                    delay = self.wave_xy_delays[qubit]
                    if self.compensate_classical_xtalk:
                        qubit2, theta2, phase2, frequency2 = self.compensate_classical_xtalk_dictionary[qubit]
                        waveform2 = self._wave_xy[qubit2]
                        gate2 = copy.copy(gate)
                        gate2.qubit = qubit2
                        #gate2.pulse = self._get_pulse_for_gate(gate2)
                        
                        gateII = gate2.gate
                        gate2.pulse = gateII.get_adjusted_pulse(self.pulses_1qb_xy[qubit2], \
                        pulse_scaling = self.single_qubit_pulse_scaling, \
                        scaling_factor_pi2 = self.single_qubit_pulse_scaling_factor[1-qubit2], \
                        detuning_pi2 = self.single_qubit_pulse_scaling_factor[1-qubit2+2])
                        
                        detuning_pi2 = self.single_qubit_pulse_scaling_factor[1-qubit2+2]
                        amp_ratio_scaling_pi2 = self.single_qubit_pulse_scaling_factor[1-qubit2+4]
                        
                        pulse2 = gate2.pulse
                        pulse2.width = self._get_pulse_for_gate(gate).width
                        pulse2.plateau = self._get_pulse_for_gate(gate).plateau
                        pulse2.drag_coefficient = self._get_pulse_for_gate(gate).drag_coefficient
                        pulse2.drag_detuning = self._get_pulse_for_gate(gate).drag_detuning
                        pulse2.amplitude *= theta2
                        pulse2.phase += phase2
                        if abs(gate_obj.theta)==np.pi/2:
                            pulse2.frequency = frequency2 + detuning_pi2
                            pulse2.amplitude *= amp_ratio_scaling_pi2
                        else:
                            pulse2.frequency = frequency2
                elif isinstance(gate_obj, gates.ReadoutGate):
                    waveform = self.readout_iq
                    delay = 0
                elif isinstance(gate_obj, gates.HeraldingGate):
                    waveform = self.heralding_iq
                    delay = 0
                elif isinstance(gate_obj, gates.HeraldingGate2):
                    waveform = self.heralding_iq_2
                    delay = 0
                else:
                    raise ValueError(
                        "Don't know which waveform to add {} to.".format(
                            gate_obj))

                # get the range of indices in use
                if (isinstance(gate_obj, gates.ReadoutGate)
                        and not self.readout_match_main_size):
                    # special case for readout if not matching main wave size
                    start = 0.0
                    end = self._round(step.t_end - step.t_start)
                else:
                    start = self._round(step.t_start + delay)
                    end = self._round(step.t_end + delay)

                if (self.compensate_crosstalk and
                    isinstance(gate_obj,
                               (gates.SingleQubitZRotation,
                                gates.TwoQubitGate))):
                    for q in range(self.n_qubit):
                        waveform = self._wave_z[q]
                        delay = self.wave_z_delays[q]
                        start = self._round(step.t_start + delay)
                        end = self._round(step.t_end + delay)
                        indices = np.arange(
                            max(np.floor(start * self.sample_rate), 0),
                            min(
                                np.ceil(end * self.sample_rate),
                                len(waveform)),
                            dtype=int)

                        # return directly if no indices
                        if len(indices) == 0:
                            continue

                        # calculate time values for the pulse indices
                        t = indices / self.sample_rate
                        max_duration = end - start
                        middle = end - max_duration / 2
                        if step.align == 'center':
                            t0 = middle
                        elif step.align == 'left':
                            t0 = middle - (max_duration - gate.duration) / 2
                        elif step.align == 'right':
                            t0 = middle + (max_duration - gate.duration) / 2

                        scaling_factor = float(crosstalk[q, 0])
                        if q != qubit:
                            scaling_factor = -scaling_factor
                        waveform[indices] += (
                            scaling_factor
                            * gate.pulse.calculate_waveform(t0, t))
                else:
                    # calculate the pulse waveform for the selected indices
                    indices = np.arange(
                        max(np.floor(start * self.sample_rate), 0),
                        min(np.ceil(end * self.sample_rate), len(waveform)),
                        dtype=int)

                    # return directly if no indices
                    if len(indices) == 0:
                        continue

                    # calculate time values for the pulse indices
                    t = indices / self.sample_rate
                    max_duration = end - start
                    middle = end - max_duration / 2
                    if step.align == 'center':
                        t0 = middle
                    elif step.align == 'left':
                        t0 = middle - (max_duration - gate.duration) / 2
                    elif step.align == 'right':
                        t0 = middle + (max_duration - gate.duration) / 2
                    if isinstance(gate_obj, gates.TwoQubitGate):
                        [phase_offset, phase_shift_frequency] = self.CZ_phase_correction
                        phase_correction = phase_offset + phase_shift_frequency * t0 * 2 * np.pi
                        gate.pulse.phase += phase_correction                      
                        
                        if self.CZ_leakage_control[0]:
                            waveform[indices] += gate.pulse.calculate_waveform(t0, t)                            
                            gate.pulse.amplitude = self._get_pulse_for_gate(gate).amplitude2
                            gate.pulse.frequency = self._get_pulse_for_gate(gate).frequency2
                            gate.pulse.drag_coefficient = self._get_pulse_for_gate(gate).drag_coefficient2
                            waveform[indices] += gate.pulse.calculate_waveform(t0, t)
                        else:
                            gate_duration = gate.duration
                            leakage_duration = gate_duration - gate.pulse.plateau + self.CZ_leakage_control[1]
                            gate.pulse.plateau -= leakage_duration
                            new_indices = np.arange(
                                max(np.floor((t0 - gate_duration / 2) * self.sample_rate), 0),
                                min(np.ceil((t0 + gate_duration / 2 - leakage_duration) * self.sample_rate), len(waveform)),
                                dtype=int)
                            new_t = new_indices / self.sample_rate
                            new_t0 = t0 - leakage_duration / 2
                            waveform[new_indices] += gate.pulse.calculate_waveform(new_t0, new_t)
                            gate.pulse.amplitude = self._get_pulse_for_gate(gate).amplitude2
                            gate.pulse.frequency = self._get_pulse_for_gate(gate).frequency2
                            gate.pulse.drag_coefficient = self._get_pulse_for_gate(gate).drag_coefficient2
                            waveform[indices] += gate.pulse.calculate_waveform(new_t0, new_t)
                    else:
                        waveform[indices] += gate.pulse.calculate_waveform(t0, t)
                    if self.compensate_classical_xtalk and isinstance(gate_obj, gates.SingleQubitXYRotation):
                        waveform2[indices] += pulse2.calculate_waveform(t0, t)
                    elif isinstance(gate_obj, gates.TwoQubitGate):
                        for XY_dict in self.CZ_XY_rotation_list:
                            qubit2 = XY_dict['index']
                            amplitude2 = XY_dict['amplitude']
                            phase2 = XY_dict['phase']
                            frequency2 = XY_dict['frequency']                            
                            phase_shift_frequency2 = XY_dict['phase shift frequency']
                            gate2 = copy.copy(GateOnQubit(gates.Xp, qubit2))
                            gate2.pulse = self._get_pulse_for_gate(gate2)
                            
                            gate3 = copy.copy(gate)
                            gate3.pulse = self._get_pulse_for_gate(gate3)
                            PhaseFromVirtualZ = gate3.pulse.phase
                            
                            pulse2 = gate2.pulse
                            pulse2.width = self._get_pulse_for_gate(gate).width
                            pulse2.plateau = self.CZ_leakage_control[1]
                            pulse2.amplitude = amplitude2
                            pulse2.phase += PhaseFromVirtualZ + phase2 + phase_shift_frequency2 * t0 * 2 * np.pi
                            pulse2.frequency = frequency2
                            # pulse2.shape = self._get_pulse_for_gate(gate).shape
                            pulse2.use_drag = self._get_pulse_for_gate(gate).use_drag
                            pulse2.drag_coefficient = self._get_pulse_for_gate(gate).drag_coefficient
                            pulse2.drag_detuning = self._get_pulse_for_gate(gate).drag_detuning
                            pulse2.truncation_range = self._get_pulse_for_gate(gate).truncation_range
                            pulse2.start_at_zero = self._get_pulse_for_gate(gate).start_at_zero

                            if self.CZ_leakage_control[0]:
                                self._wave_xy[qubit2][indices] += pulse2.calculate_waveform(t0, t)
                            else:
                                new_indices = np.arange(
                                    max(np.floor((t0 + gate_duration / 2 - leakage_duration) * self.sample_rate), 0),
                                    min(np.ceil((t0 + gate_duration / 2) * self.sample_rate), len(waveform)),
                                    dtype=int)
                                new_t = new_indices / self.sample_rate
                                new_t0 = t0 + gate_duration / 2 - leakage_duration / 2
                                self._wave_xy[qubit2][new_indices] += pulse2.calculate_waveform(new_t0, new_t)
                            

    def set_parameters(self, config={}):
        """Set base parameters using config from from Labber driver.

        Parameters
        ----------
        config : dict
            Configuration as defined by Labber driver configuration window

        """
        # sequence parameters
        d = dict(
            Zero=0,
            One=1,
            Two=2,
            Three=3,
            Four=4,
            Five=5,
            Six=6,
            Seven=7,
            Eight=8,
            Nine=9)

        # If the number of qubits changed, re-init to update pulses etc
        if self.n_qubit != d[config.get('Number of qubits')]:
            self.__init__(d[config.get('Number of qubits')])

        self.single_qubit_pulse_scaling = config.get('Single qubit pulse scaling')
        if config.get('Non trivial pulse scaling'):
            self.single_qubit_pulse_scaling_factor = []
            for n in range(self.n_qubit):
                m = n + 1
                self.single_qubit_pulse_scaling_factor.append(config.get('pi/2 scaling #%d' % m))
            for n in range(self.n_qubit):
                m = n + 1
                self.single_qubit_pulse_scaling_factor.append(config.get('pi/2 detuning #%d' % m))
            for n in range(self.n_qubit):
                m = n + 1
                self.single_qubit_pulse_scaling_factor.append(config.get('pi/2 leakage scaling #%d' % m))
        else:
            self.single_qubit_pulse_scaling_factor = [0.5 for _ in range(d[config.get('Number of qubits')])]
        self.dt = config.get('Pulse spacing')
        
        #Get the Waveform assignation
        self.waveform_assignation = config.get('Waveform assignation')
        
        #Define the Local XY boolean if the Waveform assignation is Local XY
        if self.waveform_assignation == 'Local XY':
            self.local_xy = True
        elif self.waveform_assignation == 'Manual':
            self.local_xy = False
            self.QubitWaveformDictionary = dict()
            #Dictionary with the associated waveforms
            for n in range(self.n_qubit):
                m = n + 1
                self.QubitWaveformDictionary[n] = int(config.get('Qubit %d waveform' % m))-1
        else:
            self.local_xy = False
            
        # default for simultaneous pulses is true, only option for benchmarking
        self.simultaneous_pulses = config.get('Simultaneous pulses', True)

        # waveform parameters
        self.sample_rate = config.get('Sample rate')
        self.n_pts = int(config.get('Number of points', 0))
        self.first_delay = config.get('First pulse delay')
        self.trim_to_sequence = config.get('Trim waveform to sequence')
        self.trim_start = config.get('Trim both start and end')
        self.align_to_end = config.get('Align pulses to end of waveform')

        # Compose gates
        self.compensate_classical_xtalk = config.get('Compensate Classical xtalks')
        self.compose_singleqb_gate_virtualZ = config.get('Compensate Stark Shift with Virtual Z')
        self.compose_VirtualZ_dictionary = dict()
        self.compensate_classical_xtalk_dictionary = dict()
        if self.compose_singleqb_gate_virtualZ:
            for n in range(self.n_qubit):
                m = n + 1
                self.compose_VirtualZ_dictionary[n] = \
                [d[config.get('Z rotated Qubit #{}'.format(m))],
                gates.VirtualZGate(config.get('Z rotation pi/2 #{}'.format(m))),
                gates.VirtualZGate(config.get('Z rotation pi #{}'.format(m)))]

        if self.compensate_classical_xtalk:
            for n in range(self.n_qubit):
                m = n + 1
                self.compensate_classical_xtalk_dictionary[n] = \
                [d[config.get('Classical rotated Qubit #{}'.format(m))]-1,
                config.get('Leakage amplitude #{}'.format(m)),
                config.get('Leakage phase #{}'.format(m)),
                config.get('Leakage frequency #{}'.format(m))]

        # qubit spectra
        for n in range(self.n_qubit):
            m = n + 1  # pulses are indexed from 1 in Labber
            qubit = qubits.Transmon(
                config.get('f01 max #{}'.format(m)),
                config.get('f01 min #{}'.format(m)),
                config.get('Ec #{}'.format(m)),
                config.get('Vperiod #{}'.format(m)),
                config.get('Voffset #{}'.format(m)),
                config.get('V0 #{}'.format(m)),
            )
            self.qubits[n] = qubit

        # single-qubit pulses XY
        for n, pulse in enumerate(self.pulses_1qb_xy):
            m = n + 1  # pulses are indexed from 1 in Labber
            pulse = (getattr(pulses, config.get('Pulse type'))(complex=True))
            # global parameters
            pulse.truncation_range = config.get('Truncation range')
            pulse.start_at_zero = config.get('Start at zero')
            pulse.use_drag = config.get('Use DRAG')
            # pulse shape
            if config.get('Uniform pulse shape'):
                pulse.width = config.get('Width')
                pulse.plateau = config.get('Plateau')
            else:
                pulse.width = config.get('Width #%d' % m)
                pulse.plateau = config.get('Plateau #%d' % m)

            if config.get('Uniform amplitude'):
                pulse.amplitude = config.get('Amplitude')
            else:
                pulse.amplitude = config.get('Amplitude #%d' % m)

            # pulse-specific parameters
            pulse.frequency = config.get('Frequency #%d' % m)
            pulse.drag_coefficient = config.get('DRAG scaling #%d' % m)
            pulse.drag_detuning = config.get('DRAG frequency detuning #%d' % m)

            self.pulses_1qb_xy[n] = pulse

        # single-qubit pulses Z
        for n, pulse in enumerate(self.pulses_1qb_z):
            # pulses are indexed from 1 in Labber
            m = n + 1
            # global parameters
            pulse = (getattr(pulses,
                             config.get('Pulse type, Z'))(complex=True))
            pulse.truncation_range = config.get('Truncation range, Z')
            pulse.start_at_zero = config.get('Start at zero, Z')            
            # pulse shape
            if config.get('Uniform pulse shape, Z'):
                pulse.width = config.get('Width, Z')
                pulse.plateau = config.get('Plateau, Z')
            else:
                pulse.width = config.get('Width #%d, Z' % m)
                pulse.plateau = config.get('Plateau #%d, Z' % m)

            if config.get('Uniform amplitude, Z'):
                pulse.amplitude = config.get('Amplitude, Z')
            else:
                pulse.amplitude = config.get('Amplitude #%d, Z' % m)

            pulse.frequency = config.get('Frequency #%d, Z' % m)

            self.pulses_1qb_z[n] = pulse

        # two-qubit pulses
        self.ovr_offset = config.get('Overall offset, 2QB #12')        
        for n, pulse in enumerate(self.pulses_2qb):
            # pulses are indexed from 1 in Labber
            s = ' #%d%d' % (n + 1, n + 2)
            # global parameters
            pulse = (getattr(pulses,
                             config.get('Pulse type, 2QB'))(complex=True))

            if config.get('Pulse type, 2QB') in ['CZ', 'NetZero']:
                pulse.F_Terms = d[config.get('Fourier terms, 2QB')]
                if config.get('Uniform 2QB pulses'):
                    pulse.width = config.get('Width, 2QB')
                    pulse.plateau = config.get('Plateau, 2QB')
                else:
                    pulse.width = config.get('Width, 2QB' + s)
                    pulse.plateau = config.get('Plateau, 2QB')

                # spectra
                if config.get('Assume linear dependence' + s, True):
                    pulse.qubit = None
                else:
                    pulse.qubit = self.qubits[n]

                # Get Fourier values
                if d[config.get('Fourier terms, 2QB')] == 4:
                    pulse.Lcoeff = np.array([
                        config.get('L1, 2QB' + s),
                        config.get('L2, 2QB' + s),
                        config.get('L3, 2QB' + s),
                        config.get('L4, 2QB' + s)
                    ])
                elif d[config.get('Fourier terms, 2QB')] == 3:
                    pulse.Lcoeff = np.array([
                        config.get('L1, 2QB' + s),
                        config.get('L2, 2QB' + s),
                        config.get('L3, 2QB' + s)
                    ])
                elif d[config.get('Fourier terms, 2QB')] == 2:
                    pulse.Lcoeff = np.array(
                        [config.get('L1, 2QB' + s),
                         config.get('L2, 2QB' + s)])
                elif d[config.get('Fourier terms, 2QB')] == 1:
                    pulse.Lcoeff = np.array([config.get('L1, 2QB' + s)])

                pulse.Coupling = config.get('Coupling, 2QB' + s)
                pulse.Offset = config.get('f11-f20 initial, 2QB' + s)
                pulse.amplitude = config.get('f11-f20 final, 2QB' + s)
                pulse.dfdV = config.get('df/dV, 2QB' + s)
                pulse.negative_amplitude = config.get('Negative amplitude' + s)

                pulse.calculate_cz_waveform()

            else:
                pulse.truncation_range = config.get('Truncation range, 2QB')
                pulse.start_at_zero = config.get('Start at zero, 2QB')
                # pulse shape
                if config.get('Uniform 2QB pulses'):
                    pulse.width = config.get('Width, 2QB')
                    pulse.plateau = config.get('Plateau, 2QB')
                else:
                    pulse.width = config.get('Width, 2QB' + s)
                    pulse.plateau = config.get('Plateau, 2QB' + s)
                # pulse-specific parameters
                pulse.amplitude = config.get('Amplitude, 2QB' + s)
                pulse.frequency = config.get('Frequency, 2QB' + s)
                pulse.use_drag = config.get('Use DRAG, 2QB')
                pulse.drag_coefficient = config.get('DRAG scaling, 2QB' + s)
                pulse.drag_detuning = config.get('DRAG frequency detuning, 2QB' + s)
                pulse.amplitude2 = config.get('Amplitude #2, 2QB #12')
                pulse.frequency2 = config.get('Frequency #2, 2QB #12')
                pulse.drag_coefficient2 = config.get('DRAG scaling #2, 2QB #12')
                

            gates.CZ.new_angles(
                config.get('QB1 Phi 2QB #12'), config.get('QB2 Phi 2QB #12'), config.get('QB1 Phi before and after 2QB #12'), config.get('QB2 Phi before and after 2QB #12'), VZ_factor=config.get('Virtual Z phase factor #12'),
                CPhase_pulse_number=config.get('CPhase pulse number #12'), phi_offset_step=config.get('Phi offset step #12'))
            self.TargetQubit_VZ_2QB = config.get('Virtual Z target qubit #12')            
            self.CZ_phase_correction = [config.get('Phase offset, 2QB #12'),
                                        config.get('Phase shift frequency, 2QB #12')]
            self.CZ_leakage_control = [config.get('Simultaneous leakage #12'),
                                        config.get('Leakage plateau #12')]
            self.CZ_XY_rotation_list = []
            for m in range(2):
                XY_dict = {
                'index': m,
                'amplitude': config.get('QB{} amplitude 2QB #12'.format(m + 1)),
                'phase': config.get('QB{} pulse phase 2QB #12'.format(m + 1)),
                'frequency': config.get('QB{} frequency 2QB #12'.format(m + 1)),
                'phase shift frequency': config.get('QB{} phase shift frequency 2QB #12'.format(m + 1)),
                }
                self.CZ_XY_rotation_list += [XY_dict]
                
            self.pulses_2qb[n] = pulse

        # predistortion
        self.perform_predistortion = config.get('Predistort waveforms', False)
        # update all predistorting objects
        for p in self._predistortions:
            p.set_parameters(config)

        # Z predistortion
        self.perform_predistortion_z = config.get('Predistort Z')
        for p in self._predistortions_z:
            p.set_parameters(config)

        # crosstalk
        self.compensate_crosstalk = config.get('Compensate cross-talk', False)
        self._crosstalk.set_parameters(config)

        # gate switch waveform
        self.generate_gate_switch = config.get('Generate gate')
        self.uniform_gate = config.get('Uniform gate')
        self.gate_delay = config.get('Gate delay')
        self.gate_overlap = config.get('Gate overlap')
        self.minimal_gate_time = config.get('Minimal gate time')

        # filters
        self.use_gate_filter = config.get('Filter gate waveforms', False)
        self.gate_filter = config.get('Gate filter', 'Kaiser')
        self.gate_filter_size = int(config.get('Gate - Filter size', 5))
        self.gate_filter_kaiser_beta = config.get(
            'Gate - Kaiser beta', 14.0)
        self.use_z_filter = config.get('Filter Z waveforms', False)
        self.z_filter = config.get('Z filter', 'Kaiser')
        self.z_filter_size = int(config.get('Z - Filter size', 5))
        self.z_filter_kaiser_beta = config.get(
            'Z - Kaiser beta', 14.0)

        # readout
        self.readout_match_main_size = config.get(
            'Match main sequence waveform size')
        self.readout_i_offset = config.get('Readout offset - I')
        self.readout_q_offset = config.get('Readout offset - Q')
        self.readout_trig_generate = config.get('Generate readout trig')
        self.readout_trig_amplitude = config.get('Readout trig amplitude')
        self.readout_trig_duration = config.get('Readout trig duration')
        self.readout_trig_delay = config.get('Readout trig delay')
        self.heralding_trig_generate = config.get('Generate heralding trigger')
        self.heralding_trig_delay = config.get('Heralding trigger delay')
        self.heralding_trig_generate_2 = config.get('Generate heralding trigger #2')
        self.heralding_trig_delay_2 = config.get('Heralding trigger delay #2')
        self.readout_predistort = config.get('Predistort readout waveform')
        self.readout.set_parameters(config)

        # cavity leakage
        self.activate_cavity_leakage = config.get('Activate cavity leakage')
        self.leakage_amplitude = config.get('Leakage amplitude')
        self.leakage_frequency = config.get('Leakage frequency')
        self.leakage_phase = config.get('Leakage phase')
        # get readout pulse parameters
        phases = 2 * np.pi * np.array([
            0.8847060, 0.2043214, 0.9426104, 0.6947334, 0.8752361, 0.2246747,
            0.6503154, 0.7305004, 0.1309068
        ])
        for n, pulse in enumerate(self.pulses_readout):
            # pulses are indexed from 1 in Labber
            m = n + 1
            pulse = (getattr(pulses,
                             config.get('Readout pulse type'))(complex=True))
            pulse.truncation_range = config.get('Readout truncation range')
            pulse.start_at_zero = config.get('Readout start at zero')
            pulse.iq_skew = config.get('Readout IQ skew') * np.pi / 180
            pulse.iq_ratio = config.get('Readout I/Q ratio')

            if config.get('Distribute readout phases'):
                pulse.phase = phases[n]
            else:
                pulse.phase = 0

            if config.get('Uniform readout pulse shape'):
                pulse.width = config.get('Readout width')
                pulse.plateau = config.get('Readout duration')
            else:
                pulse.width = config.get('Readout width #%d' % m)
                pulse.plateau = config.get('Readout duration #%d' % m)

            if config.get('Uniform readout amplitude') is True:
                pulse.amplitude = config.get('Readout amplitude')
            else:
                pulse.amplitude = config.get('Readout amplitude #%d' % (n + 1))

            pulse.frequency = config.get('Readout frequency #%d' % m)
            self.pulses_readout[n] = pulse
            
        #Heralding pulse
        for n, pulse in enumerate(self.pulses_heralding):
            # pulses are indexed from 1 in Labber
            m = n + 1
            pulse = (getattr(pulses,
                             config.get('Readout pulse type'))(complex=True))
            pulse.truncation_range = config.get('Readout truncation range')
            pulse.start_at_zero = config.get('Readout start at zero')
            pulse.iq_skew = config.get('Readout IQ skew') * np.pi / 180
            pulse.iq_ratio = config.get('Readout I/Q ratio')

            if config.get('Distribute readout phases'):
                pulse.phase = phases[n]
            else:
                pulse.phase = 0

            if config.get('Uniform readout pulse shape'):
                pulse.width = config.get('Readout width')
                pulse.plateau = config.get('Heralding duration factor')*config.get('Readout duration')
            else:
                pulse.width = config.get('Readout width #%d' % m)
                pulse.plateau = config.get('Heralding duration factor')*config.get('Readout duration #%d' % m)

            if config.get('Uniform readout amplitude') is True:
                pulse.amplitude = config.get('Heralding amplitude factor')*config.get('Readout amplitude')
            else:
                pulse.amplitude = config.get('Heralding amplitude factor')*config.get('Readout amplitude #%d' % (n + 1))

            pulse.frequency = pulse.frequency = config.get('Heralding frequency')
            self.pulses_heralding[n] = pulse
        
        #Heralding pulse 2
        for n, pulse in enumerate(self.pulses_heralding_2):
            # pulses are indexed from 1 in Labber
            m = n + 1
            pulse = (getattr(pulses,
                             config.get('Readout pulse type'))(complex=True))
            pulse.truncation_range = config.get('Readout truncation range')
            pulse.start_at_zero = config.get('Readout start at zero')
            pulse.iq_skew = config.get('Readout IQ skew') * np.pi / 180
            pulse.iq_ratio = config.get('Readout I/Q ratio')

            if config.get('Distribute readout phases'):
                pulse.phase = phases[n]
            else:
                pulse.phase = 0

            if config.get('Uniform readout pulse shape'):
                pulse.width = config.get('Readout width')
                pulse.plateau = config.get('Heralding duration factor #2')*config.get('Readout duration')
            else:
                pulse.width = config.get('Readout width #%d' % m)
                pulse.plateau = config.get('Heralding duration factor #2')*config.get('Readout duration #%d' % m)

            if config.get('Uniform readout amplitude') is True:
                pulse.amplitude = config.get('Heralding amplitude factor #2')*config.get('Readout amplitude')
            else:
                pulse.amplitude = config.get('Heralding amplitude factor #2')*config.get('Readout amplitude #%d' % (n + 1))

            pulse.frequency = config.get('Heralding frequency #2')
            self.pulses_heralding_2[n] = pulse
        # Delays
        self.wave_xy_delays = np.zeros(self.n_qubit)
        self.wave_z_delays = np.zeros(self.n_qubit)
        for n in range(self.n_qubit):
            m = n + 1
            self.wave_xy_delays[n] = config.get('Qubit %d XY Delay' % m)
            self.wave_z_delays[n] = config.get('Qubit %d Z Delay' % m)


if __name__ == '__main__':
    pass
