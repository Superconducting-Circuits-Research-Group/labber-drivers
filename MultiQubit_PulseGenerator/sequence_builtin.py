#!/usr/bin/env python3
# add logger, to allow logging to Labber's instrument log
import logging

import numpy as np

from gates import Gate, IdentityGate
from sequence import Sequence

log = logging.getLogger('LabberDriver')


class Rabi(Sequence):
    """Sequence for driving Rabi oscillations in multiple qubits."""

    def generate_sequence(self, config):
        """Generate sequence by adding gates/pulses to waveforms."""
        # just add pi-pulses for the number of available qubits
        self.add_gate_to_all(Gate.Xp, align='right')


class CPMG(Sequence):
    """Sequence for multi-qubit Ramsey/Echo/CMPG experiments."""

    def generate_sequence(self, config):
        """Generate sequence by adding gates/pulses to waveforms."""
        # get parameters
        n_pulse = int(config['# of pi pulses'])
        pi_to_q = config['Add pi pulses to Q']
        duration = config['Sequence duration']
        edge_to_edge = config['Edge-to-edge pulses']

        # select type of refocusing pi pulse
        gate_pi = Gate.Yp if pi_to_q else Gate.Xp

        if edge_to_edge:
            if n_pulse < 0:
                self.add_gate_to_all(gate_pi)
                self.add_gate_to_all(IdentityGate(width=0), dt=duration)
            else:
                self.add_gate_to_all(Gate.X2p)
                dt = duration/(n_pulse+1)
                for i in range(n_pulse):
                    self.add_gate_to_all(gate_pi, dt=dt)
                self.add_gate_to_all(Gate.X2p, dt=dt)
        else:
            if n_pulse < 0:
                self.add_gate_to_all(gate_pi, t0=0)
                self.add_gate_to_all(IdentityGate(width=0), t0=duration)
            else:
                self.add_gate_to_all(Gate.X2p, t0=0)
                for i in range(n_pulse):
                    self.add_gate_to_all(gate_pi, t0=duration/(n_pulse+1)*(i+1))
                self.add_gate_to_all(Gate.X2p, t0=duration)


class PulseTrain(Sequence):
    """Sequence for multi-qubit pulse trains, for pulse calibrations."""

    def generate_sequence(self, config):
        """Generate sequence by adding gates/pulses to waveforms."""
        # get parameters
        n_pulse = int(config['# of pulses'])
        alternate = config['Alternate pulse direction']

        if n_pulse == 0:
            self.add_gate_to_all(Gate.I)
        for n in range(n_pulse):
            pulse_type = config['Pulse']
            # check if alternate pulses
            if alternate and (n % 2) == 1:
                pulse_type = pulse_type.replace('p', 'm')
                gate = Gate.__getattr__(pulse_type)
            else:
                gate = Gate.__getattr__(pulse_type)
            self.add_gate_to_all(gate)


if __name__ == '__main__':
    pass
