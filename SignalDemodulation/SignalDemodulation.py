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
        if self.isFirstCall(options):
            # clear record buffer
            self.data = {}
        if quant.name in ('Value', 'Value - Single shot'):
            if quant.name not in self.data:
                self.getIQData()
            return self.data[quant.name]
        elif quant.name.startswith('Value #'):
            if quant.name not in self.data:
                index = int(quant.name[-1])
                sDemodFreq = 'Mod. frequency #%d' % index
                dFreq = self.getValue(sDemodFreq)
                self.getIQData(index, dFreq)
            return self.data[quant.name]
        else:
            # just return the quantity value
            return quant.getValue()

    def getIQData(self, index=None, dFreq=None):
        """Calculate complex signal from data and reference."""
        # get parameters
        if index is None:
            sValue = 'Value'
            sSingleShot = 'Value - Single shot'
        else:
            sValue = 'Value #%d' % index
            sSingleShot = 'Value #%d - Single shot' % index
        if dFreq is None:
            dFreq = self.getValue('Modulation frequency')
        skipStart = self.getValue('Skip start')
        nSegment = int(self.getValue('Number of segments'))
        # get input data from dictionary, dictionary format:
        # {'y': value, 't0': t0, 'dt': dt}
        recordIn = self.getValue('Input data')
        if recordIn is None:
            self.data[sValue] = np.NaN
            self.data[sSingleShot] = np.NaN
            return
        vData = recordIn['y']
        nTotLength = vData.size
        if nTotLength % nSegment != 0:
            raise ValueError('Total record length should be exactly '
                             'a multiple of segments.')
        dt = recordIn['dt']
        # avoid exceptions if no time step is given
        if dt == 0:
            dt = 1.0
        skipIndex = int(round(skipStart / dt))
        length = 1 + int(round(self.getValue('Length') / dt))
        segmentLength = int(nTotLength / nSegment)
        length = min(length, segmentLength - skipIndex)
        if length <= 1:
            self.data[sValue] = np.NaN
            self.data[sSingleShot] = np.NaN
            return
        # define data to use, put in 2D array of segments
        vData.shape = (nSegment, segmentLength)
        # calculate cos/sin vectors, allow segmenting
        vTime = dt * (skipIndex + np.arange(length, dtype=float))
        vExp = np.exp(2.j * np.pi * vTime * dFreq)
        # calc I/Q
        signal = np.empty(nSegment, dtype=complex)
        np.dot(vData[:,skipIndex:skipIndex+length], vExp, signal)
        signal /= .5 * float(length - 1)
        bUseRef = bool(self.getValue('Use phase reference signal'))
        if bUseRef:
            vRef = self.getValue('Reference data')
            # skip reference if record length doesn't match
            vRef = vRef['y']
            if vRef.size != nTotLength:
                raise ValueError('The data and reference record '
                        'lengths do not match.')
            vRef.shape = (nSegment, segmentLength)
            ref = np.empty(nSegment, dtype=complex)
            np.dot(vRef[:,skipIndex:skipIndex+length], vExp, ref)
            # subtract the reference angle
            signal *= np.exp(-1j * np.angle(ref))
        self.data[sValue] = np.mean(signal)
        self.data[sSingleShot] = signal
        return


if __name__ == '__main__':
    pass
