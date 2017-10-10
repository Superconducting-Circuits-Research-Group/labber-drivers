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
        if quant.name in ('Value', 'Value - Single shot',
                          'Time average'):
            if quant.name not in self.data:
                self.getIQData()
            if quant.name == 'Time average':
                return quant.getTraceDict(self.data['Time average'],
                                          dt=self.data['dt'])
            else:
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
        sTimeAverage = 'Time average'
        if index is None:
            sValue = 'Value'
            sSingleShot = 'Value - Single shot'
        else:
            sValue = 'Value #%d' % index
            sSingleShot = 'Value #%d - Single shot' % index
        self.data[sTimeAverage] = np.NaN
        self.data[sSingleShot] = np.NaN
        self.data[sValue] = np.NaN
        
        # get parameters
        if dFreq is None:
            dFreq = self.getValue('Modulation frequency')
        skipStart = self.getValue('Skip start')
        nSegment = int(self.getValue('Number of segments'))
        # get input data from dictionary, dictionary format:
        # {'y': value, 't0': t0, 'dt': dt}
        recordIn = self.getValue('Input data')
        if recordIn is None:
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
        self.data['dt'] = dt
        
        skipIndex = int(round(skipStart / dt))
        length = 1 + int(round(self.getValue('Length') / dt))
        segmentLength = int(nTotLength / nSegment)
        length = min(length, segmentLength - skipIndex)
        if length < 1:
            return

        bRef = bool(self.getValue('Use phase reference signal'))

        # define data to use, put in 2D array of segments
        vData.shape = (nSegment, segmentLength)
        # calculate cos/sin vectors, allow segmenting
        if not hasattr(self, '_dt') or dt != self._dt or \
                not hasattr(self, '_skipIndex') or \
                skipIndex != self._skipIndex or \
                not hasattr(self, '_length') or \
                length != self._length or \
                not hasattr(self, '_dFreq') or \
                dFreq != self._dFreq:
            vTime = dt * (skipIndex + np.arange(length, dtype=np.float32))
            self._vExp = np.exp(2.j * np.pi * vTime * dFreq)
            self._dt = dt
            self._skipIndex = skipIndex
            self._length = length
            self._dFreq = dFreq
        vExp = self._vExp

        # calc I/Q
        if not hasattr(self, '_nSegment') or nSegment != self._nSegment:
            self._signal = np.empty(nSegment, dtype=np.complex64)
            if bRef:
                self._ref = np.empty(nSegment, dtype=np.complex64)
            self._nSegment = nSegment
        signal = self._signal
        np.dot(vData[:,skipIndex:skipIndex+length], vExp, signal)
        signal /= .5 * float(length)

        if bRef:
            vRef = self.getValue('Reference data')
            # skip reference if record length doesn't match
            vRef = vRef['y']
            if vRef.size != nTotLength:
                raise ValueError('The data and reference record '
                        'lengths do not match.')
            vRef.shape = (nSegment, segmentLength)
            ref = self._ref
            np.dot(vRef[:,skipIndex:skipIndex+length], vExp, ref)
            # subtract the reference angle
            expRef = np.exp(-1j * np.angle(ref))
            signal *= expRef

        self.data[sSingleShot] = signal
        self.data[sValue] = np.mean(signal)
        
        if bRef:
            if not hasattr(self, '_segmentLength') or \
                segmentLength != self._segmentLength:
                self._timeAverage = np.empty(segmentLength, dtype=np.complex64)
                self._segmentLength = segmentLength
            timeAverage = self._timeAverage
            np.dot(vData.T, expRef, timeAverage)
        else:
            timeAverage = np.sum(vData, 0)
        timeAverage /= nSegment
        self.data[sTimeAverage] = timeAverage


if __name__ == '__main__':
    pass
