import time
import numpy as np

import InstrumentDriver
import AlazarTech_SignalDemodulator_Wrapper as AlazarDig


class Driver(InstrumentDriver.InstrumentWorker):
    """This class implements the AlazarTech card driver."""

    def performOpen(self, options={}):
        """Perform the operation of opening the instrument connection."""
        # keep track of sampled records
        self.data = {}
        self.dt = 1.0
        # open connection
        boardId = int(self.comCfg.address)
        timeout = self.dComCfg['Timeout']
        self.dig = AlazarDig.AlazarTechDigitizer(systemId=1,
                boardId=boardId, timeout=timeout)
        self.dig.testLED()

    def performClose(self, bError=False, options={}):
        """Perform the close instrument connection operation."""
        # try to remove buffers
        self.dig.removeBuffersDMA()
        # remove digitizer object
        del self.dig

    def performSetValue(self, quant, value, sweepRate=0.0, options={}):
        """Perform the Set Value instrument operation. This function
        should return the actual value set by the instrument.
        """
        # start with setting current quant value
        quant.setValue(value)
        # don't do anything until all options are set, then set
        # complete config
        if self.isFinalCall(options):
            self.setConfiguration()
        return value

    def performGetValue(self, quant, options={}):
        """Perform the Get Value instrument operation."""
        # only implemented for records
        mode = self.getValue('Acquisition mode')
        if quant.name in ('Channel A - Records',
                          'Channel B - Records',
                          'Channel A - Average record',
                          'Channel B - Average record',
                          'Channel A - Demodulated values',
                          'Channel B - Demodulated values',
                          'Channel A - Average demodulated value',
                          'Channel B - Average demodulated value'):
            # special case for hardware looping
            if self.isHardwareLoop(options):
                self.getSignalHardwareLoop(quant, options)
                return self.getData(mode, quant)
            
            # check if first call, if so get new records
            if self.isFirstCall(options):
                # clear trace dictionary
                self.data = {}
                # read records
                if self.getValue('NPT AsyncDMA Enabled'):
                    self.getRecordsDMA(
                        hardware_trig=self.isHardwareTrig(options))
                else:
                    self.getRecordsSinglePort()
                    
            if mode.startswith('Referenced'):
                self._bRef = True
            else:
                self._bRef = False

            # cash demodulation parameters
            if mode.startswith('Individual') or \
                    mode.startswith('Referenced Individual'):
                self.cashDemodulationParametersForIndividual()

            # return correct data
            return self.getData(mode, quant)
        else:
            # just return the quantity value
            return quant.getValue()

    def performArm(self, quant_names, options={}):
        """Perform the instrument arm operation."""
        # get config
        nSample = int(self.getValue('Number of samples'))
        nRecord = int(self.getValue('Number of records'))
        recordsPerBuffer = int(self.getValue('Records per buffer'))
        maxBufferSize = int(self.getValue('Max buffer size'))
        nMaxBuffer = int(self.getValue('Max number of buffers'))
        mode = self.getValue('Acquisition mode')
        # configure and start acquisition
        if self.isHardwareLoop(options):
            # in hardware looping, number of records is set by
            # the hardware looping
            (seq_no, n_seq) = self.getHardwareLoopIndex(options)
            # need to re-configure the card since record size was not
            # known at config
            self.dig.readRecordsDMA(mode, nSample, n_seq,
                    bConfig=True, bArm=True, bMeasure=False,
                    maxBuffers=maxBuffers,
                    maxBufferSize=maxBufferSize)
        else:
            # if not hardware looping, just trigger the card, buffers
            # are already configured
            self.dig.readRecordsDMA(mode, nSample, nRecord,
                    bConfig=False, bArm=True, bMeasure=False,
                    maxBuffers=maxBuffers,
                    maxBufferSize=maxBufferSize)

    def _callbackProgress(self, progress):
        """Report progress to server, as text string."""
        s = 'Acquiring records (%.0f%%).' % (100 * progress)
        self.reportStatus(s)

    def getSignalHardwareLoop(self, quant, options):
        """Get data from round-robin type averaging."""
        (seq_no, n_seq) = self.getHardwareLoopIndex(options)
        # if first sequence call, get data
        if seq_no == 0 and self.isFirstCall(options):
            nSample = int(self.getValue('Number of samples'))
            maxBufferSize = int(self.getValue('Max buffer size'))
            maxBuffers = int(self.getValue('Max number of buffers'))
            mode = self.getValue('Acquisition mode')
            # show status before starting acquisition
            self.reportStatus('Digitizer - Waiting for signal')
            # get data
            self.data = self.dig.readRecordsDMA(mode, nSample, n_seq,
                    bConfig=False, bArm=False, bMeasure=True,
                    funcStop=self.isStopped,
                    funcProgress=self._callbackProgress,
                    firstTimeout=self.dComCfg['Timeout'] + 180.0,
                    maxBuffers=maxBuffers, maxBufferSize=maxBufferSize)

    def setConfiguration(self):
        """Set digitizer configuration based on driver settings."""
        sampleRateIndex = self.getValueIndex('Sample rate')
        # clock configuration
        SourceId = int(self.getCmdStringFromValue('Clock source'))
        if self.getValue('Clock source') == 'Internal':
            # internal
            SampleRateId = int(self.getCmdStringFromValue('Sample rate'), 0)
            lFreq = [1E3, 2E3, 5E3, 10E3, 20E3, 50E3, 100E3, 200E3,
                     500E3, 1E6, 2E6, 5E6, 10E6, 20E6, 50E6, 100E6,
                     200E6, 500E6, 1E9, 1.2E9, 1.5E9, 2E9, 2.4E9, 3E9,
                     3.6E9, 4E9]
            Decimation = 0
        elif self.getValue('Clock source') == '10 MHz Reference' and \
                self.getModel() in ('9373', '9360'):
            # 10 MHz ref, for 9373 - decimation is 1
            # for now don't allow DES mode; talk to Simon about best
            # implementation
            lFreq = [1E3, 2E3, 5E3, 10E3, 20E3, 50E3, 100E3, 200E3,
                     500E3, 1E6, 2E6, 5E6, 10E6, 20E6, 50E6, 100E6,
                     200E6, 500E6, 1E9, 1.2E9, 1.5E9, 2E9, 2.4E9, 3E9,
                     3.6E9, 4E9]
            SampleRateId = int(lFreq[sampleRateIndex])
            Decimation = 0
        else:
            # 10 MHz ref, use 1GHz rate + divider; divider must be 1, 2,
            # 4, or multiple of 10
            SampleRateId = int(1E9)
            lFreq = [1E3, 2E3, 5E3, 10E3, 20E3, 50E3, 100E3, 200E3,
                     500E3, 1E6, 2E6, 5E6, 10E6, 20E6, 50E6, 100E6,
                     250E6, 500E6, 1E9]
            Decimation = int(1E9 / lFreq[sampleRateIndex])
        try:
            self.dig.AlazarSetCaptureClock(SourceId, SampleRateId, 0,
                                           Decimation)
        except:
            self.log("'ApiPllNotLocked' has been thrown.")
            time.sleep(.1)
            self.dig.AlazarSetCaptureClock(SourceId, SampleRateId, 0,
                                           Decimation)
        # define time step from sample rate
        self.dt = 1.0 / lFreq[sampleRateIndex]
        # configure inputs
        chnls = {1: 'Channel A', 2: 'Channel B'}
        for n in (1, 2):
            # coupling and range
            if self.getModel() in ('9373', '9360'):
                # these options are not available for these models,
                # set to default
                Coupling = 2
                InputRange = 7
                Impedance = 2
            else:
                Coupling = int(self.getCmdStringFromValue('%s - '
                        'Coupling' % chnls[n]))
                InputRange = int(self.getCmdStringFromValue('%s - '
                        'Range' % chnls[n]))
                Impedance = int(self.getCmdStringFromValue('%s - '
                        'Impedance' % chnls[n]))
            # set coupling, input range, impedance
            self.dig.AlazarInputControl(n, Coupling, InputRange,
                                        Impedance)
            # bandwidth limit, only for model 9870
            if self.getModel() in ('9870',):
                self.dig.AlazarSetBWLimit(n, 0)

        # configure trigger
        Source = int(self.getCmdStringFromValue('Trigger source'))
        Slope = int(self.getCmdStringFromValue('Trigger slope'))
        Delay = self.getValue('Trigger delay')
        timeout = self.dComCfg['Timeout']
        # trigger level is relative to full range
        trigLevel = self.getValue('Trigger level')
        vAmp = np.array([4, 2, 1, 0.4, 0.2, 0.1, .04], dtype=float)
        if self.getValue('Trigger source') == 'Channel A':
            maxLevel = vAmp[self.getValueIndex('Channel A - Range')]
        elif self.getValue('Trigger source') == 'Channel B':
            maxLevel = vAmp[self.getValueIndex('Channel B - Range')]
        elif self.getValue('Trigger source') == 'External':
            maxLevel = 5.0
        elif self.getValue('Trigger source') == 'Immediate':
            maxLevel = 5.0
            # set timeout to very short with immediate triggering
            timeout = 0.001
        # convert relative level to U8
        if abs(trigLevel) > maxLevel:
            trigLevel = maxLevel * np.sign(trigLevel)
        Level = int(128 + 127 * trigLevel / maxLevel)
        # set configuration
        self.dig.AlazarSetTriggerOperation(0, 0, Source, Slope, Level)

        # configure external input, if in use
        if self.getValue('Trigger source') == 'External':
            Coupling = int(self.getCmdStringFromValue('Trigger coupling'))
            self.dig.AlazarSetExternalTrigger(Coupling)

        # set trigger delay and timeout
        Delay = int(self.getValue('Trigger delay') / self.dt)
        self.dig.AlazarSetTriggerDelay(Delay)
        self.dig.AlazarSetTriggerTimeOut(time=timeout)
        # configure memory buffers, only possible when using DMA read
        nPostSize = int(self.getValue('Number of samples'))
        nRecord = int(self.getValue('Number of records'))
        maxBufferSize = int(self.getValue('Max buffer size'))
        maxBuffers = int(self.getValue('Max number of buffers'))
        mode = self.getValue('Acquisition mode')
        # configure DMA read
        self.dig.readRecordsDMA(mode, nPostSize, nRecord,
                bConfig=True, bArm=False, bMeasure=False,
                maxBuffers=maxBuffers, maxBufferSize=maxBufferSize)

    def getRecordsDMA(self, hardware_trig=False):
        """Resample the data with DMA."""
        # get channels in use
        nPostSize = int(self.getValue('Number of samples'))
        nRecord = int(self.getValue('Number of records'))
        maxBufferSize = int(self.getValue('Max buffer size'))
        maxBuffers = int(self.getValue('Max number of buffers'))
        mode = self.getValue('Acquisition mode')
        # in hardware triggering mode, there is no noed to re-arm
        # the card
        bArm = not hardware_trig
        # get data
        self.data = self.dig.readRecordsDMA(mode,
                nPostSize, nRecord,
                bConfig=False, bArm=bArm, bMeasure=True,
                funcStop=self.isStopped,
                maxBuffers=maxBuffers, maxBufferSize=maxBufferSize)

    def getRecordsSinglePort(self):
        """Resample the data."""
        # get channels in use
        nPreSize = int(self.getValue('Pre-trigger samples'))
        nPostSize = int(self.getValue('Number of samples'))
        nRecord = int(self.getValue('Number of records'))

        model = self.getModel()
        if model == '9870':
            pAlgn = 64
            rAlgn = 256
            nPreAlgn = pAlgn * ((nPreSize + pAlgn - 1) // pAlgn)
            nPostAlgn = pAlgn * ((nPostSize + pAlgn - 1) // pAlgn)
            nTotal = rAlgn * ((nPreAlgn + nPostAlgn + rAlgn - 1) // rAlgn)
            nPostAlgn = nTotal - nPreAlgn
            start = nPreAlgn - nPreSize
            end = start + nPreSize + nPostSize
        else:
            raise NotImplementedError('Model ATS%s is not yet supported.'
                                      % model)

        self.dig.AlazarSetRecordSize(nPreAlgn, nPostAlgn)
        self.dig.AlazarSetRecordCount(nRecord)
        # start aquisition
        self.dig.AlazarStartCapture()
        nTry = self.dComCfg['Timeout'] / 0.05
        while nTry and self.dig.AlazarBusy() and not self.isStopped():
            # sleep for a while to save resources, then try again
            self.wait(0.050)
            nTry -= 1
        # check if timeout occurred
        if nTry <= 0:
            self.dig.AlazarAbortCapture()
            raise TimeoutError('Acquisition timed out.')
        # check if user stopped
        if self.isStopped():
            self.dig.AlazarAbortCapture()
            return

        # read data for channels in use
        records = self.dig.readRecordsSinglePort(1)
        self.data['Channel A'] = records[:,start:end]
        
        records = self.dig.readRecordsSinglePort(2)
        self.data['Channel B'] = records[:,start:end]

    def _raiseError(self, name, mode):
            raise NotImplementedError("Output data '%s' could not be "
                    "acquired in acquisition mode '%s'." % (name, mode))

    def getData(self, mode, quant):
        """Return data that corresponds to the seleted data acquisition
        mode."""
        name = quant.name
        if mode == 'Raw':
            if name in ('Channel A - Records',
                        'Channel B - Records'):
                flattened = self.data[name[:9]].ravel()
                return quant.getTraceDict(flattened, dt=self.dt)
            else:
                raise self._raiseError(name, mode)
                
        elif mode == 'Average Record':
            if name in ('Channel A - Average record',
                        'Channel B - Average record'):
                if name in self.data:
                    return quant.getTraceDict(self.data[name],
                                              dt=self.dt)
                raise self._raiseError(name, mode)
                
        elif mode == 'Individual Record Demodulation':
            if name in ('Channel A - Records',
                        'Channel B - Records'):
                flattened = self.data[name[:9]].ravel()
                return quant.getTraceDict(flattened, dt=self.dt)
            elif name in ('Channel A - Demodulated values',
                                'Channel B - Demodulated values'):
                if name not in self.data:
                    self.getDemodulatedValuesFromIndividual()
                return quant.getTraceDict(self.data[name],
                                          dt=1)
            elif name in ('Channel A - Average demodulated value',
                          'Channel B - Average demodulated value'):
                if name not in self.data:
                    t0 = time.clock()
                    self.getDemodulatedValuesFromIndividual()
                    self.log('Data processing %.6f s.' % (time.clock() - t0))
                    if hasattr(self, '_tcycle'):
                        self.log('Total time %.6f s.' % (time.clock() - self._tcycle))
                    self._tcycle = time.clock()
                return self.data[name]
            elif name in ('Channel A - Average record',
                                'Channel B - Average record'):
                if name not in self.data:
                    self.getAverageRecordFromIndiviual()
                return quant.getTraceDict(self.data[name],
                                          dt=self.dt)
            else:
                raise self._raiseError(name, mode)
                
        elif mode == 'Referenced Individual Record Demodulation':
            if name in ('Channel A - Records',
                        'Channel B - Records'):
                flattened = self.data[name[:9]].ravel()
                return quant.getTraceDict(flattened, dt=self.dt)
            elif name == 'Channel A - Demodulated values':
                if name not in self.data:
                    self.getDemodulatedValuesFromIndividual()
                return quant.getTraceDict(self.data[name], dt=1)
            elif name == 'Channel A - Average demodulated value':
                if name not in self.data:
                    t0 = time.clock()
                    self.getDemodulatedValuesFromIndividual()
                    self.log('Data processing %.6f s.' % (time.clock() - t0))
                return self.data[name]
            elif name == 'Channel A - Average record':
                if name not in self.data:
                    self.getAverageRecordFromIndiviual()
                return quant.getTraceDict(self.data[name],
                                          dt=self.dt)
            else:
                raise self._raiseError(name, mode)
                
        elif mode == 'Average Record Demodulation':
            if name in ('Channel A - Average record',
                        'Channel B - Average record'):
                if name in self.data:
                    return quant.getTraceDict(self.data[name],
                                              dt=self.dt)
                else:
                    self._raiseError(name, mode)
            elif name in ('Channel A - Average demodulated value',
                          'Channel B - Average demodulated value'):
                self.getDemodulatedValueFromAverage()
                return self.data[name]
            else:
                raise self._raiseError(name, mode)
                              
        elif mode == 'Referenced Average Record Demodulation':
            if name in ('Channel A - Average record',
                        'Channel B - Average record'):
                if name in self.data:
                    return quant.getTraceDict(self.data[name],
                                              dt=self.dt)
                else:
                    self._raiseError(name, mode)
            elif name == 'Channel A - Average demodulated value':
                t0 = time.clock()
                self.getDemodulatedValueFromAverage()
                self.log('Data processing %.6f s.' % (time.clock() - t0))
                if hasattr(self, '_tcycle'):
                    self.log('Total time %.6f s.' % (time.clock() - self._tcycle))
                self._tcycle = time.clock()
                return self.data[name]
            else:
                raise self._raiseError(name, mode)

    def cashDemodulationParametersForIndividual(self):
        """Cash parameters that are used for calculating output values."""
        # get parameters
        dFreq = self.getValue('Demodulation frequency')
        skipStart = self.getValue('Skip start')
        nSegment = int(self.getValue('Number of records'))
        mode = self.getValue('Acquisition mode')

        nTotLength = self.data['Channel A'].size
        dt = self.dt
        skip = int(round(skipStart / dt))
        length = 1 + int(round(self.getValue('Demodulation length') / dt))
        segmentLength = int(nTotLength / nSegment)
        length = min(length, segmentLength - skip)
        if length < 1:
            raise Exception('Nothing to demodulate: increase '
                    'the demodulation or record length.')

        if not hasattr(self, '_firstRun'):
            self._firstRun = True
        # calculate cos/sin vectors, allow segmenting
        if self._firstRun or dt != self._dt or \
                skip != self._skip or \
                length != self._length or \
                dFreq != self._dFreq:
            vTime = dt * (skip + np.arange(length, dtype=np.float32))
            self._vExp = np.exp(2.j * np.pi * vTime * dFreq).view('complex64')
            self._dt = dt
            self._skip = skip
            self._length = length
            self._dFreq = dFreq
        
        if self._firstRun or self._nSegment != nSegment:
            self._nSegment = nSegment
        
        if self._bRef:
            vChB = self.data['Channel B'][:,skip:skip+length]
            vDemodRefs = np.dot(vChB, self._vExp)
            self._vExpRef = np.exp(-1.j * np.angle(vDemodRefs))

        self._firstRun = False

    def getAverageRecordFromIndiviual(self):
        """Calculate time average from data and reference."""
        if self._bRef:
            timeAverage = np.dot(self.data['Channel A'].T, self._vExpRef)
            timeAverage /= self._nSegment
            self.data['Channel A - Average record'] = timeAverage
        else:
            for ch in ('Channel A', 'Channel B'):
                timeAverage = np.sum(self.data[ch], axis=0)
                timeAverage /= self._nSegment
                self.data['%s - Average record' % ch] = timeAverage

    def getDemodulatedValuesFromIndividual(self):
        """Calculate complex signal vector from data and reference."""
        skip = self._skip
        length = self._length
        
        if self._bRef:
            vChA = self.data['Channel A'][:,skip:skip+length]
            vDemodVals = np.dot(vChA, self._vExp)
            vDemodVals /= .5 * np.float32(length)
            vDemodVals *= self._vExpRef
            self.data['Channel A - Demodulated values'] = vDemodVals
            self.data['Channel A - Average demodulated value'] = np.mean(vDemodVals)
        else:
            for ch in ('Channel A', 'Channel B'):
                vCh = self.data[ch][:,skip:skip+length]
                vDemodVals = np.dot(vCh, self._vExp)
                vDemodVals /= .5 * np.float32(length)
                self.data['%s - Demodulated values' % ch] = vDemodVals
                self.data['%s - Average demodulated value' % ch] = np.mean(vDemodVals)
                
    def getDemodulatedValueFromAverage(self):
        """Calculate complex signal vector from data and reference."""
        # get parameters
        dFreq = self.getValue('Demodulation frequency')
        skipStart = self.getValue('Skip start')
        mode = self.getValue('Acquisition mode')

        dt = self.dt
        skip = int(round(skipStart / dt))
        length = 1 + int(round(self.getValue('Demodulation length') / dt))
        vTime = dt * (skip + np.arange(length, dtype=np.float32))
        vExp = np.exp(2.j * np.pi * vTime * dFreq).view('complex64')

        if self._bRef:
            vChA = self.data['Channel A - Average record'][skip:skip+length]
            vDemodVal = np.dot(vChA, vExp)
            vDemodVal /= .5 * np.float32(length)
            vChB = self.data['Channel B - Average record'][skip:skip+length]
            vDemodRef = np.dot(vChB, vExp)
            vDemodVal *= np.exp(-1.j * np.angle(vDemodRef))
            self.data['Channel A - Average demodulated value'] = vDemodVal
        else:
            for ch in ('Channel A', 'Channel B'):
                vCh = self.data['%s - Average record' % ch][skip:skip+length]
                vDemodVal = np.dot(vCh, vExp)
                vDemodVal /= .5 * np.float32(length)
                self.data['%s - Average demodulated value' % ch] = vDemodVal


if __name__ == '__main__':
    pass
