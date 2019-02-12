#!/usr/bin/env python3
# add logger, to allow logging to Labber's instrument log
import logging

import numpy as np

from pulse import Pulse, PulseShape

from gates import Gate, BaseGate, IdentityGate, CustomGate
from sequence import Sequence


class CustomSequence(Sequence):
    """Generate sequence by adding gates/pulses to waveforms."""
    
    def generate_sequence(self, config):
        """Generate sequence by adding gates/pulses to waveforms."""
        pulse1 = Pulse(shape=PulseShape.SQUARE)
        pulse1.amplitude = config['Parameter #1']
        pulse1.width = float(config['Number of points']) / float(config['Sample rate'])
        pulse3 = Pulse(shape=PulseShape.SQUARE)
        pulse3.amplitude = config['Parameter #1']
        pulse3.width = config['Readout duration'] + float(config['Parameter #2'])
        self.add_gate_to_all(CustomGate(pulse1), t0=0)
        self.add_gate_to_all(Gate.Xp, dt=0)
        self.add_gate_to_all(CustomGate(pulse3), dt=0)

if __name__ == '__main__':
    pass
