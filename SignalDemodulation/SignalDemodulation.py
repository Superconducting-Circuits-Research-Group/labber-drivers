import numpy as np

import InstrumentDriver


class Driver(InstrumentDriver.InstrumentWorker):
    """This class implements a demodulation driver."""

    def performOpen(self, options={}):
        """Perform the operation of opening the instrument connection."""
        pass

    def performClose(self, bError=False, options={}):
        """Perform the close instrument connection operation."""
        pass

    def performSetValue(self, quant, value, sweepRate=0.0, options={}):
        """Perform the Set Value instrument operation. This function
        should return the actual value set by the instrument."""
        return value

    def performGetValue(self, quant, options={}):
        """Perform the Get Value instrument operation."""
        if quant.name == 'Value':
            value = np.mean(self.getIQAmplitudes())
        elif quant.name == 'Value - Single shot':
            value = self.getIQAmplitudes()
        elif quant.name.startswith('Value #'):
            index = int(quant.name[-1])
            sDemodFreq = 'Mod. frequency #%d' % index
            dFreq = self.getValue(sDemodFreq)
            value = np.mean(self.getIQAmplitudes_MultiFreq(dFreq))
        else:
            # just return the quantity value
            value = quant.getValue()
        return value

    def getIQAmplitudes(self, dFreq=None):
        """Calculate complex signal from data and reference."""
        # get parameters
        if dFreq is None:
            dFreq = self.getValue('Modulation frequency')
        skipStart = self.getValue('Skip start')
        nSegment = int(self.getValue('Number of segments'))
        # get input data from dictionary, dictionary format:
        # {'y': value, 't0': t0, 'dt': dt}
        traceIn = self.getValue('Input data')
        if traceIn is None:
            return complex(0.0)
        vY = traceIn['y']
        nTotLength = vY.size
        if nTotLength % nSegment != 0:
            raise ValueError('Total record length should be exactly '
                             'a multiple of segments.')
        dt = traceIn['dt']
        # avoid exceptions if no time step is given
        if dt == 0:
            dt = 1.0
        skipIndex = int(round(skipStart / dt))
        length = 1 + int(round(self.getValue('Length') / dt))
        segmentLength = int(nTotLength / nSegment)
        length = min(length, segmentLength - skipIndex)
        if length <= 1:
            return complex(0.0)
        # define data to use, put in 2d array of segments
        vData = np.reshape(vY, (nSegment, segmentLength))
        # calculate cos/sin vectors, allow segmenting
        vTime = dt * (skipIndex + np.arange(length, dtype=float))
        vPhase = 2 * np.pi * vTime * dFreq
        vExp = np.exp(1j * vPhase)
        # calc I/Q
        signal = np.empty(nSegment, dtype=complex)
        np.dot(vData[:,skipIndex:skipIndex+length], vExp, signal)
        signal /= .5 * float(length - 1)
        bUseRef = bool(self.getValue('Use phase reference signal'))
        if bUseRef:
            traceRef = self.getValue('Reference data')
            # skip reference if trace length doesn't match
            if len(traceRef['y']) != len(vY):
                return signal
            vRef = np.reshape(traceRef['y'], (nSegment, segmentLength))
            ref = np.empty(nSegment, dtype=complex)
            np.dot(vData[:,skipIndex:skipIndex+length], vExp, ref)
            # subtract the reference angle
            signal -= np.exp(1j * np.angle(ref))
        return signal


if __name__ == '__main__':
    pass
