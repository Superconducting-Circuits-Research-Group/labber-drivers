import numpy as np

import InstrumentDriver


class Driver(InstrumentDriver.InstrumentWorker):
    """This class implements a readout pulse signal generator."""

    def performGetValue(self, quant, options={}):
        """Perform the Get Value instrument operation."""
        if quant.name == 'Readout waveform':
            N = self.getValue('Number of points')
            delay = self.getValue('Readout pulse delay')
            length = self.getValue('Readout pulse length')
            zero = np.array([0])
            Ndelay = int((N - 2) * delay / (delay + length))
            Nlength = N - 2 - Ndelay
            return np.hstack([np.zeros(Ndelay + 1),
                              np.ones(Nlength),
                              np.zeros(1)])
        else:
            # just return the quantity value
            return quant.getValue()


if __name__ == '__main__':
    pass
