#!/usr/bin/env python3
from copy import copy
import numpy as np
import logging
from sequence import Step
log = logging.getLogger('LabberDriver')

# TODO remove Step dep from CompositeGate


class BaseGate:
    """Base class for a qubit gate.

    """

    def get_adjusted_pulse(self, pulse):
        pulse = copy(pulse)
        return pulse

    def __repr__(self):
        return self.__str__()


class OneQubitGate(BaseGate):
    def number_of_qubits(self):
        return 1


class TwoQubitGate(BaseGate):
    def number_of_qubits(self):
        return 2


class SingleQubitXYRotation(OneQubitGate):
    """Single qubit rotations around the XY axes.

    Angles defined as in https://en.wikipedia.org/wiki/Bloch_sphere.

    Parameters
    ----------
    phi : float
        Rotation axis.
    theta : float
        Roation angle.

    """

    def __init__(self, phi, theta, name=None, plateau=None):
        self.phi = phi
        self.theta = theta
        self.name = name
        self.plateau = plateau
        # log.info(str(self.plateau))
        # log.info(str(self.name))

    def get_adjusted_pulse(self, pulse, pulse_scaling='Amplitude', scaling_factor_pi2=0.5):
        log.info('before the error:' + str(self.name))
        
        pulse = copy(pulse)
        pulse.phase = self.phi
        self.pulse_scaling = pulse_scaling
        if self.plateau is not None:
            pulse.plateau = self.plateau
            # log.info(str(pulse.plateau))
        if self.pulse_scaling == 'Amplitude':
        # pi pulse correspond to the full amplitude
            pulse.amplitude *= self.theta / np.pi
            if self.theta==np.pi / 2:
                pulse.amplitude *= 2*scaling_factor_pi2
        elif self.pulse_scaling == 'Duration':
            pulse.plateau *= self.theta / np.pi
            pulse.width *= self.theta / np.pi
            if self.theta==np.pi / 2:
                pulse.plateau *= 2 * scaling_factor_pi2
                pulse.width *= 2 * scaling_factor_pi2
        return pulse
        
    def new_scaling(self, new_pulse_scaling):
        """Update the pulse scaling

        """
        self.__init__(self.phi, self.theta, self.name, pulse_scaling=new_pulse_scaling)

    def __str__(self):
        if self.name is None:
            return "XYPhi={:+.6f}theta={:+.6f}".format(self.phi, self.theta)
        else:
            return self.name

    def __eq__(self, other):
        threshold = 1e-10
        if not isinstance(other, SingleQubitXYRotation):
            return False
        if np.abs(self.phi - other.phi) > threshold:
            return False
        if np.abs(self.theta - other.theta) > threshold:
            return False
        return True


class SingleQubitZRotation(OneQubitGate):
    """Single qubit rotation around the Z axis.

    Parameters
    ----------
    theta : float
        Roation angle.

    """

    def __init__(self, theta, name=None):
        self.theta = theta
        self.name = name

    def get_adjusted_pulse(self, pulse, pulse_scaling='Amplitude'):
        pulse = copy(pulse)
        pulse.phase = self.phi
        self.pulse_scaling = pulse_scaling
        if self.pulse_scaling == 'Amplitude':
        # pi pulse correspond to the full amplitude
            pulse.amplitude *= self.theta / np.pi
        elif self.pulse_scaling == 'Duration':
            pulse.plateau *= self.theta / np.pi
            pulse.width *= self.theta / np.pi
        return pulse

    def __str__(self):
        if self.name is None:
            return "Ztheta={:+.2f}".format(self.theta)
        else:
            return self.name

    def __eq__(self, other):
        threshold = 1e-10
        if not isinstance(other, SingleQubitZRotation):
            return False
        if np.abs(self.theta - other.theta) > threshold:
            return False
        return True


class IdentityGate(OneQubitGate):
    """Identity gate.

    Does nothing to the qubit. The width can be specififed to
    implement a delay in the sequence. If no width is given, the identity gate
    inherits the width of the given pulse.

    Parameters
    ----------
    width : float
        Width of the I gate in seconds,
        None uses the XY width (the default is None).

    """

    def __init__(self, width=None):
        self.width = width

    def get_adjusted_pulse(self, pulse):
        pulse = copy(pulse)
        pulse.amplitude = 0
        pulse.use_drag = False  # Avoids bug
        if self.width is not None:
            pulse.width = 0
            pulse.plateau = self.width
        return pulse

    def __str__(self):
        return "I"


class VirtualZGate(OneQubitGate):
    """Virtual Z Gate."""

    def __init__(self, theta, name=None):
        self.theta = theta
        self.name = name

    def __eq__(self, other):
        threshold = 1e-10
        if not isinstance(other, VirtualZGate):
            return False
        if np.abs(self.theta - other.theta) > threshold:
            return False
        return True

    def __str__(self):
        if self.name is None:
            return "VZtheta={:+.2f}".format(self.theta)
        else:
            return self.name


class CPHASE(TwoQubitGate):
    """ CPHASE gate. 
    
    Parameters
    ----------
    phi : float
        Rotation axis.
    VZ_factor : float
        Roation axis factor.
        
    """
    def __init__(self, phi, VZ_factor=0.5, phi_offset=0):
        self.phi = phi
        self.VZ_factor = VZ_factor
        self.phi_offset = phi_offset
        
    def get_adjusted_pulse(self, pulse):
        # log.info('self.phi is ' + str(self.phi))
        pulse = copy(pulse)
        pulse.phase = self.phi * self.VZ_factor + self.phi_offset
        return pulse

class ReadoutGate(OneQubitGate):
    """Readouts the qubit state."""
	
class HeraldingGate(OneQubitGate):
    """Heraldings the qubit state."""

class HeraldingGate2(OneQubitGate):
    """Heraldings the qubit state."""                                  

class CustomGate(BaseGate):
    """A gate using a given :obj:`Pulse`.

    Parameters
    ----------
    pulse : :obj:`Pulse`
        The corresponding pulse.

    """

    def __init__(self, pulse):
        self.pulse = pulse


class RabiGate(SingleQubitXYRotation):
    """Creates the Rabi gate used in the spin-locking sequence.

    Parameters
    ----------
    amplitude : Amplitude of the pulse
    plateau : The duration of the pulse.
    phase : Phase of the Rabi gate. 0 corresponds to rotation around X axis.
    """

    def __init__(self, amplitude, plateau, phase):
        self.amplitude = amplitude
        self.plateau = plateau
        self.phase = phase

    def get_adjusted_pulse(self, pulse):
        pulse = copy(pulse)
        pulse.amplitude = self.amplitude
        pulse.plateau = self.plateau
        pulse.phase = self.phase
        return pulse


class CompositeGate:
    """Multiple gates in one object.

    Parameters
    ----------
    n_qubit : int
        Number of qubits involved in the composite gate.

    Attributes
    ----------
    sequence : list of :Step:
        Holds the gates involved.

    """

    def __init__(self, n_qubit, name=None):
        self.n_qubit = n_qubit
        self.sequence = []
        self.name = name

    def add_gate(self, gate, qubit=None):
        """Add a set of gates to the given qubit.

        For the qubits with no specificied gate, an IdentityGate will be given.
        The length of the step is given by the longest pulse.

        Parameters
        ----------
        qubit : int or list of int
            The qubit(s) to add the gate(s) to.
        gate : :obj:`BaseGate` or list of :obj:`BaseGate`
            The gate(s) to add.
        """
        if qubit is None:
            if self.n_qubit == 1:
                qubit = 0
            else:
                qubit = [n for n in range(self.n_qubit)]

        step = Step()
        if isinstance(gate, list):
            if len(gate) == 1:
                raise ValueError(
                    "For single gates, don't provide gate as a list.")
            if not isinstance(qubit, list):
                raise ValueError(
                    """Please provide qubit indices as a list when adding more
                    than one gate.""")
            if len(gate) != len(qubit):
                raise ValueError(
                    "Length of gate list must equal length of qubit list.")

            for q, g in zip(qubit, gate):
                step.add_gate(q, g)

        else:
            if gate.number_of_qubits() > 1:
                if not isinstance(qubit, list):
                    raise ValueError(
                        """Please provide qubit list for gates with more than
                        one qubit.""")
            else:
                if not isinstance(qubit, int):
                    raise ValueError(
                        "For single gates, give qubit as int (not list).")
            step.add_gate(qubit, gate)

        self.sequence.append(step)

    def number_of_qubits(self):
        return self.n_qubit

    def __len__(self):
        return len(self.sequence)

    def __str__(self):
        if self.name is not None:
            return self.name
        else:
            super().__str__()

    def __repr__(self):
        return self.__str__()


class CPHASE_with_1qb_phases(CompositeGate):
    """CPHASE gate followed by single qubit Z rotations.

    Parameters
    ----------
    phi1 : float
        Z rotation angle for qubit 1.
    phi2 : float
        Z rotation angle for qubit 2.

    """

    def __init__(self, phi1, phi2, VZ_factor=0, CPhase_pulse_number=1, phi_offset_step=np.pi):
        super().__init__(n_qubit=2)
        self.phi1 = phi1
        self.phi2 = phi2
        self.VZ_factor = VZ_factor
        self.CPhase_pulse_number = max(int(CPhase_pulse_number), 1)
        self.phi_offset_step = phi_offset_step
        for i in range(self.CPhase_pulse_number):
            self.add_gate(CPHASE(0, VZ_factor=VZ_factor, phi_offset=phi_offset_step * i))
        self.add_gate([VirtualZGate(phi1), VirtualZGate(phi2)])

    def new_angles(self, phi1, phi2, VZ_factor=0, CPhase_pulse_number=1, phi_offset_step=np.pi):
        """Update the angles of the single qubit rotations and the VZ phase factor on its own. Update whether to use echo pulse.

        Parameters
        ----------
        phi1 : float
            Z rotation angle for qubit 1.
        phi2 : float
            Z rotation angle for qubit 2.
        VZ_factor : float
            VZ rotation angle factor.

        """
        self.__init__(phi1, phi2, VZ_factor=VZ_factor, CPhase_pulse_number=CPhase_pulse_number, phi_offset_step=phi_offset_step)
    def __str__(self):
        return "CZ"


I = IdentityGate(width=None)
I0 = IdentityGate(width=0)
Ilong = IdentityGate(width=75e-9)

# X gates
Xp = SingleQubitXYRotation(phi=0, theta=np.pi, name='Xp')
Xm = SingleQubitXYRotation(phi=0, theta=-np.pi, name='Xm')
X2p = SingleQubitXYRotation(phi=0, theta=np.pi / 2, name='X2p')
X2m = SingleQubitXYRotation(phi=0, theta=-np.pi / 2, name='X2m')
Xp2 = SingleQubitXYRotation(phi=0, theta = 2 * np.pi, name='Xp2')
Xm2 = SingleQubitXYRotation(phi=0, theta = - 2 * np.pi, name='Xp2')

# Y gates
Yp = SingleQubitXYRotation(phi=np.pi / 2, theta=np.pi, name='Yp')
Ym = SingleQubitXYRotation(phi=np.pi / 2, theta=-np.pi, name='Ym')
Y2m = SingleQubitXYRotation(phi=np.pi / 2, theta=-np.pi / 2, name='Y2m')
Y2p = SingleQubitXYRotation(phi=np.pi / 2, theta=np.pi / 2, name='Y2p')

# Z gates
Zp = SingleQubitZRotation(np.pi, name='Zp')
Z2p = SingleQubitZRotation(np.pi / 2, name='Z2p')
Zm = SingleQubitZRotation(-np.pi, name='Zm')
Z2m = SingleQubitZRotation(-np.pi / 2, name='Z2m')

# Virtual Z gates
VZp = VirtualZGate(np.pi, name='VZp')
VZ2p = VirtualZGate(np.pi / 2, name='VZ2p')
VZm = VirtualZGate(-np.pi, name='VZm')
VZ2m = VirtualZGate(-np.pi / 2, name='VZ2m')

# two-qubit gates
CPh = CPHASE(0)

# Composite gates
CZEcho = CompositeGate(n_qubit=2)
CZEcho.add_gate([X2p, I])
CZEcho.add_gate(CPh)
CZEcho.add_gate([Xp, Xp])
CZEcho.add_gate(CPh)
CZEcho.add_gate([X2p, Xp])

H = CompositeGate(n_qubit=1, name='H')
H.add_gate(VZp)
H.add_gate(Y2p)

CZ = CPHASE_with_1qb_phases(0, 0)  # Start with 0, 0 as the single qubit phase shifts.

CNOT = CompositeGate(n_qubit=2, name='CNOT')
CNOT.add_gate(H, 1)
CNOT.add_gate(CZ, [0, 1])
CNOT.add_gate(H, 1)

if __name__ == '__main__':
    pass
