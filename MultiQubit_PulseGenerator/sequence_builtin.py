#!/usr/bin/env python3
# add logger, to allow logging to Labber's instrument log
import logging
import numpy as np

import gates
from sequence import Sequence

log = logging.getLogger('LabberDriver')


class Rabi(Sequence):
    """Sequence for driving Rabi oscillations in multiple qubits."""

    def generate_sequence(self, config):
        """Generate sequence by adding gates/pulses to waveforms."""
        # just add pi-pulses for the number of available qubits
        
        manipulated_qubits = config['Manipulated Qubits']
        
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
        if manipulated_qubits=='All':
            qubit = list(range((self.n_qubit)))
        else:
            qubit = [d[manipulated_qubits]-1]
        
        gate_pi = [gates.Xp for n in range(len(qubit))]
        
        self.add_gate(qubit,gate_pi, align='right')


class CPMG(Sequence):
    """Sequence for multi-qubit Ramsey/Echo/CMPG experiments."""

    def generate_sequence(self, config):
        """Generate sequence by adding gates/pulses to waveforms."""
        # get parameters
        n_pulse = int(config['# of pi pulses'])
        pi_to_q = config['Add pi pulses to Q']
        duration = config['Sequence duration']
        edge_to_edge = config['Edge-to-edge pulses']
        manipulated_qubits = config['Manipulated Qubits']
        
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
        if manipulated_qubits=='All':
            qubit = list(range((self.n_qubit)))
        else:
            qubit = [d[manipulated_qubits]-1]
        
        # select type of refocusing pi pulse
        gate_pi = [(gates.Yp if pi_to_q else gates.Xp) for n in range(len(qubit))]
        gate_X2p = [gates.X2p for n in range(len(qubit))]
        
        # always do T1 same way, regardless if edge-to-edge or center-center
        if n_pulse < 0:
            self.add_gate(qubit,gate_pi)
            self.add_gate_to_all(gates.IdentityGate(width=duration), dt=0)

        elif edge_to_edge:
            # edge-to-edge pulsing, set pulse separations
            self.add_gate(qubit,gate_X2p)
            # for ramsey, just add final pulse
            if n_pulse == 0:
                self.add_gate(qubit,gate_X2p, dt=duration)
            else:
                dt = duration / n_pulse
                # add first pi pulse after half duration
                self.add_gate(qubit,gate_pi, dt=dt/2)
                # add rest of pi pulses
                for i in range(n_pulse - 1):
                    self.add_gate(qubit,gate_pi, dt=dt)
                # add final pi/2 pulse
                self.add_gate(qubit,gate_X2p, dt=dt/2)

        else:
            # center-to-center spacing, set absolute pulse positions
            self.add_gate(qubit,gate_X2p, t0=0)
            # add pi pulses at right position
            for i in range(n_pulse):
                self.add_gate(qubit,gate_pi,
                                     t0=(i + 0.5) * (duration / n_pulse))
            # add final pi/2 pulse
            self.add_gate(qubit,gate_X2p, t0=duration)


class PulseTrain(Sequence):
    """Sequence for multi-qubit pulse trains, for pulse calibrations."""

    def generate_sequence(self, config):
        """Generate sequence by adding gates/pulses to waveforms."""
        # get parameters
        n_pulse = int(config['# of pulses'])
        alternate = config['Alternate pulse direction']
        manipulated_qubits = config['Manipulated Qubits']
        
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
        if manipulated_qubits=='All':
            qubit = list(range((self.n_qubit)))
        else:
            qubit = [d[manipulated_qubits]-1]


        if n_pulse == 0:
            self.add_gate_to_all(gates.I)
        for n in range(n_pulse):
            pulse_type = config['Pulse']
            if pulse_type == 'CPh':
                for i in range(self.n_qubit-1):
                    self.add_gate([i, i+1], gates.CPh)
            elif pulse_type =='CZ':
                for i in range(self.n_qubit-1):   
                    self.add_gate([i, i+1],gates.CZ)
            elif pulse_type == 'I':
                self.add_gate(qubit,[gates.I for n in range(len(qubit))])
            else:
                if alternate and (n % 2) == 1:
                    pulse_type = pulse_type.replace('p', 'm')
                gate = getattr(gates, pulse_type)
                gate_on_selected_qubits = [gate for n in range(len(qubit))]
                self.add_gate(qubit,gate_on_selected_qubits)

class AllXY(Sequence):
    """Sequence for multi-qubit pulse trains, for pulse calibrations."""

    def generate_sequence(self, config):
        """Generate sequence by adding gates/pulses to waveforms."""
        # get parameters
        pulse_index = config['AllXY pulse index']
        manipulated_qubits = config['Manipulated Qubits']
        interleave_gate = config['Interleave Gate']
        if interleave_gate:
            pulse_type = config['Interleaved pulse']
            interleaved_gate = getattr(gates, pulse_type)
            inverse_pulse_type = pulse_type.replace('p', 'u')
            inverse_pulse_type = inverse_pulse_type.replace('m', 'p')
            inverse_pulse_type = inverse_pulse_type.replace('u', 'm')
            inverse_interleaved_gate = getattr(gates, inverse_pulse_type)
            manipulated_qubits_for_interleaved_pulse = config['Manipulated Qubits for interleaved pulse']
        
        
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
        if manipulated_qubits=='All':
            qubit = list(range((self.n_qubit)))
        else:
            qubit = [d[manipulated_qubits]-1]
        if interleave_gate:
            if manipulated_qubits_for_interleaved_pulse=='All':
                qubit_interleave = list(range((self.n_qubit)))
            else:
                qubit_interleave = [d[manipulated_qubits_for_interleaved_pulse]-1]

        if pulse_index == 'I-I':
            gate1 = gates.I
            gate2 = gates.I
        elif pulse_index == 'Xp-Xp':
            gate1 = gates.Xp
            gate2 = gates.Xp
        elif pulse_index == 'Yp-Yp':
            gate1 = gates.Yp
            gate2 = gates.Yp
        elif pulse_index == 'Xp-Yp':
            gate1 = gates.Xp
            gate2 = gates.Yp
        elif pulse_index == 'Yp-Xp':
            gate1 = gates.Xp
            gate2 = gates.Yp
        elif pulse_index == 'X2p-I':
            gate1 = gates.X2p
            gate2 = gates.I
        elif pulse_index == 'Y2p-I':
            gate1 = gates.Y2p
            gate2 = gates.I
        elif pulse_index == 'X2p-Y2p':
            gate1 = gates.X2p
            gate2 = gates.Y2p
        elif pulse_index == 'Y2p-X2p':
            gate1 = gates.Y2p
            gate2 = gates.X2p
        elif pulse_index == 'X2p-Yp':
            gate1 = gates.X2p
            gate2 = gates.Yp
        elif pulse_index == 'Y2p-Xp':
            gate1 = gates.Y2p
            gate2 = gates.Xp
        elif pulse_index == 'Xp-Y2p':
            gate1 = gates.Xp
            gate2 = gates.Y2p
        elif pulse_index == 'Yp-X2p':
            gate1 = gates.Yp
            gate2 = gates.X2p
        elif pulse_index == 'X2p-Xp':
            gate1 = gates.X2p
            gate2 = gates.Xp
        elif pulse_index == 'Xp-X2p':
            gate1 = gates.Xp
            gate2 = gates.X2p
        elif pulse_index == 'Y2p-Yp':
            gate1 = gates.Y2p
            gate2 = gates.Yp
        elif pulse_index == 'Yp-Y2p':
            gate1 = gates.Yp
            gate2 = gates.Y2p
        elif pulse_index == 'Xp-I':
            gate1 = gates.Xp
            gate2 = gates.I
        elif pulse_index == 'Yp-I':
            gate1 = gates.Yp
            gate2 = gates.I
        elif pulse_index == 'X2p-X2p':
            gate1 = gates.X2p
            gate2 = gates.X2p
        elif pulse_index == 'Y2p-Y2p':
            gate1 = gates.Y2p
            gate2 = gates.Y2p
            
        gate1_on_selected_qubits = [gate1 for _ in range(len(qubit))]
        gate2_on_selected_qubits = [gate2 for _ in range(len(qubit))]
            
        self.add_gate(qubit, gate1_on_selected_qubits)
        if interleave_gate:
            self.add_gate(qubit_interleave, [interleaved_gate for _ in range(len(qubit_interleave))])
        self.add_gate(qubit, gate2_on_selected_qubits)
        if interleave_gate:
            self.add_gate(qubit_interleave, [inverse_interleaved_gate for _ in range(len(qubit_interleave))])


class SpinLocking(Sequence):
    """ Sequence for spin-locking experiment.

    """

    def generate_sequence(self, config):
        """Generate sequence by adding gates/pulses to waveforms."""

        pulse_amps = []
        for ii in range(9):
            pulse_amps.append(
                float(config['Drive pulse amplitude #' + str(ii + 1)]))
        pulse_duration = float(config['Drive pulse duration'])
        pulse_phase = float(config['Drive pulse phase']) / 180.0 * np.pi
        pulse_sequence = config['Pulse sequence']

        if pulse_sequence == 'SL-3':
            self.add_gate_to_all(gates.Y2p)
        if pulse_sequence == 'SL-5a':
            self.add_gate_to_all(gates.Y2m)
        if pulse_sequence == 'SL-5b':
            self.add_gate_to_all(gates.Y2p)

        if pulse_sequence != 'SL-3':
            self.add_gate_to_all(gates.Xp)

        rabi_gates = []
        for ii in range(self.n_qubit):
            rabi_gates.append(
                gates.RabiGate(pulse_amps[ii], pulse_duration, pulse_phase))
        self.add_gate(list(range(self.n_qubit)), rabi_gates)
        if pulse_sequence != 'SL-3':
            self.add_gate_to_all(gates.Xp)

        if pulse_sequence == 'SL-3':
            self.add_gate_to_all(gates.Y2p)
        if pulse_sequence == 'SL-5a':
            self.add_gate_to_all(gates.Y2m)
        if pulse_sequence == 'SL-5b':
            self.add_gate_to_all(gates.Y2p)

        return

if __name__ == '__main__':
    pass
