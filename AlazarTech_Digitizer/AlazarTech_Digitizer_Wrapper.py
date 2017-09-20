import os
import time
import numpy as np
import ctypes

from atsapi import *
from ctypes import (c_int, c_int32, c_uint8, c_uint16, c_uint32, 
                    c_float, c_char_p, c_void_p, c_long, byref)

# add logger, to allow logging to Labber's instrument log
import logging
log = logging.getLogger('LabberDriver')


class AlazarTechDigitizer():
    """Represents the AlazarTech digitizer, redefines the DLL functions
    in Python.
    """

    def __init__(self, systemId=1, boardId=1, timeout=10.0):
        """Defines a session ID, used to identify the instrument."""
        # range settings; default value of 400 mV for model 9373;
        # will be overwritten if model is 9870 and AlazarInputControl is
        # called
        self.dRange = {1: 0.4, 2: 0.4}
        self.buffers = []
        self.timeout = timeout
        # create a session ID
        func = getattr(ats, 'AlazarNumOfSystems')
        func.restype = U32
        func = getattr(ats, 'AlazarGetBoardBySystemID')
        func.restype = c_void_p
        handle = func(U32(systemId), U32(boardId))
        if handle is None:
            raise Exception('Device with system ID=%d and board ID=%d '
                            'could not be found.' % (systemId, boardId))
        self.handle = handle
        # get memory and bitsize
        (self.memorySize_samples, self.bitsPerSample) = \
            self.AlazarGetChannelInfo()

    def testLED(self):
        self.callFunc('AlazarSetLED', self.handle, U32(1))
        time.sleep(0.1)
        self.callFunc('AlazarSetLED', self.handle, U32(0))

    def callFunc(self, sFunc, *args, **kargs):
        """General function caller with restype=status,
        also checks for errors.
        """
        # get function from DLL
        func = getattr(ats, sFunc)
        func.restype = c_int
        # call function, raise error if needed
        status = func(*args)
        if 'bIgnoreError' in kargs:
            bIgnoreError = kargs['bIgnoreError']
        else:
            bIgnoreError = False
        if status is not None and status > 512 and not bIgnoreError:
            raise Exception(self.getError(status))

    def getError(self, status):
        """Convert the error in status to a string."""
        func = getattr(DLL, 'AlazarErrorToText')
        func.restype = c_char_p
        # const char* AlazarErrorToText(RETURN_CODE retCode)
        errorText = func(c_int(status))
        return str(errorText)

    def AlazarGetChannelInfo(self):
        """Get the on-board memory in samples per channel and sample
        size in bits per sample.
        """
        memorySize_samples = U32(0)
        bitsPerSample = U8(0)
        self.callFunc('AlazarGetChannelInfo', self.handle,
                      byref(memorySize_samples), byref(bitsPerSample))
        return (int(memorySize_samples.value), int(bitsPerSample.value))

    # RETURN_CODE AlazarSetCaptureClock(HANDLE h, U32 Source, U32 Rate,
    #       U32 Edge, U32 Decimation);
    def AlazarSetCaptureClock(self, SourceId, SampleRateId, EdgeId=0,
                              Decimation=0):
        self.callFunc('AlazarSetCaptureClock', self.handle, U32(SourceId),
                      U32(SampleRateId), U32(EdgeId), U32(Decimation))

    # RETURN_CODE AlazarInputControl(HANDLE h, U8 Channel,
    #       U32 Coupling, U32 InputRange, U32 Impedance);
    def AlazarInputControl(self, Channel, Coupling, InputRange, Impedance):
        # keep track of input range
        dConv = {12: 4.0, 11: 2.0, 10: 1.0, 7: 0.4,
                 6: 0.2, 5: 0.1, 2: 0.04}
        self.dRange[Channel] = dConv[InputRange]
        self.callFunc('AlazarInputControl', self.handle, U8(Channel),
                      U32(Coupling), U32(InputRange), U32(Impedance))

    # RETURN_CODE AlazarSetBWLimit(HANDLE h, U8 Channel, U32 enable);
    def AlazarSetBWLimit(self, Channel, enable):
        self.callFunc('AlazarSetBWLimit', self.handle, U32(Channel),
                      U32(enable))

    # RETURN_CODE AlazarSetTriggerOperation(HANDLE h, U32 TriggerOperation,
    #       U32 TriggerEngine1/*J,K*/, U32 Source1, U32 Slope1, U32 Level1,
    #       U32 TriggerEngine2/*J,K*/, U32 Source2, U32 Slope2, U32 Level2);
    def AlazarSetTriggerOperation(self, TriggerOperation=0,
            TriggerEngine1=0, Source1=0, Slope1=1, Level1=128,
            TriggerEngine2=1, Source2=3, Slope2=1, Level2=128):
        self.callFunc('AlazarSetTriggerOperation', self.handle,
            U32(TriggerOperation),
            U32(TriggerEngine1), U32(Source1), U32(Slope1), U32(Level1),
            U32(TriggerEngine2), U32(Source2), U32(Slope2), U32(Level2))

    # RETURN_CODE AlazarSetExternalTrigger(HANDLE h, U32 Coupling, U32 Range);
    def AlazarSetExternalTrigger(self, Coupling, Range=0):
        self.callFunc('AlazarSetExternalTrigger', self.handle,
                      U32(Coupling), U32(Range))

    # RETURN_CODE AlazarSetTriggerDelay( HANDLE h, U32 Delay);
    def AlazarSetTriggerDelay(self, Delay=0):
        self.callFunc('AlazarSetTriggerDelay', self.handle, U32(Delay))

    # RETURN_CODE AlazarSetTriggerTimeOut( HANDLE h, U32 to_ns);
    def AlazarSetTriggerTimeOut(self, time=0.0):
        tick = U32(int(time * 1E5))
        self.callFunc('AlazarSetTriggerTimeOut', self.handle, tick)

    # RETURN_CODE AlazarSetRecordSize(HANDLE h, U32 PreSize, U32 PostSize);
    def AlazarSetRecordSize(self, PreSize, PostSize):
        self.nPreSize = int(PreSize)
        self.nPostSize = int(PostSize)
        self.callFunc('AlazarSetRecordSize', self.handle, U32(PreSize),
                      U32(PostSize))

    # RETURN_CODE AlazarSetRecordCount(HANDLE h, U32 Count);
    def AlazarSetRecordCount(self, Count):
        self.nRecord = int(Count)
        self.callFunc('AlazarSetRecordCount', self.handle, U32(Count))

    # RETURN_CODE AlazarStartCapture(HANDLE h);
    def AlazarStartCapture(self):
        self.callFunc('AlazarStartCapture', self.handle)

    # RETURN_CODE AlazarAbortCapture(HANDLE h);
    def AlazarAbortCapture(self):
        self.callFunc('AlazarAbortCapture', self.handle)

    # U32 AlazarBusy(HANDLE h);
    def AlazarBusy(self):
        # get function from DLL
        func = getattr(ats, 'AlazarBusy')
        func.restype = U32
        # call function, return result
        return bool(func(self.handle))

    # U32 AlazarRead(HANDLE h, U32 Channel, void *buffer, int ElementSize,
    #       long Record, long TransferOffset, U32 TransferLength);
    def AlazarRead(self, Channel, buffer, ElementSize,
                   Record, TransferOffset, TransferLength):
        self.callFunc('AlazarRead', self.handle, U32(Channel), buffer,
            c_int(ElementSize), c_long(Record), c_long(TransferOffset),
            U32(TransferLength))

    def AlazarBeforeAsyncRead(self, channels, transferOffset,
                              samplesPerRecord, recordsPerBuffer,
                              recordsPerAcquisition, flags):
        """Prepares the board for an asynchronous acquisition."""
        self.callFunc('AlazarBeforeAsyncRead', self.handle,
                      channels, transferOffset, samplesPerRecord,
                      recordsPerBuffer, recordsPerAcquisition, flags)

    # RETURN_CODE AlazarAbortAsyncRead( HANDLE h);
    def AlazarAbortAsyncRead(self):
        """Cancels any asynchronous acquisition running on a board."""
        self.callFunc('AlazarAbortAsyncRead', self.handle)

    def AlazarPostAsyncBuffer(self, buffer, bufferLength):
        """Posts a DMA buffer to a board."""
        self.callFunc('AlazarPostAsyncBuffer', self.handle, buffer,
                      bufferLength)

    def AlazarWaitAsyncBufferComplete(self, buffer, timeout_ms):
        """Blocks until the board confirms that buffer is filled with
        data.
        """
        self.callFunc('AlazarWaitAsyncBufferComplete', self.handle,
                      buffer, timeout_ms)

    def readTracesDMA(self, bGetChA, bGetChB,
                      nSamples, nRecord,
                      bConfig=True, bArm=True, bMeasure=True,
                      funcStop=None, funcProgress=None,
                      timeout=None, firstTimeout=None,
                      maxBuffers=8, maxBufferSizeMB=64,
                      bGetAllTraces=False):
        """Reads traces in NPT AutoDMA mode, converts to float,
        demodulates/averages to a single trace.
        """
        t0 = time.clock()
        lT = []
        # use global timeout if not given
        timeout = self.timeout if timeout is None else timeout
        # first timeout can be different in case of slow initial arming
        firstTimeout = timeout if firstTimeout is None else firstTimeout

        # change alignment to be 64
        samplesPerRecord = 64 * ((nSamples + 63) // 64)

        # select the active channels
        channels = int(bGetChA) | int(bGetChB) << 1
        numberOfChannels = int(bGetChA) + int(bGetChB)

        data = {}
        if numberOfChannels == 0:
            return data

        # compute the number of bytes per record and per buffer
        bytesPerSample = (self.bitsPerSample + 7) // 8
        bytesPerRecord = bytesPerSample * samplesPerRecord
        maxMemSize = int(maxBufferSizeMB * 2**20)
        buffersPerAcquisition = 1
        recordsPerBuffer = nRecord
        consistencyFlag = False
        while not consistencyFlag:
            bytesPerBuffer = int(numberOfChannels * recordsPerBuffer *
                                 bytesPerRecord)
            if bytesPerBuffer > maxMemSize:
                for factor in range(1, int(np.sqrt(nRecord)) + 1):
                    if nRecord % factor == 0:
                        recordsPerBuffer = int(recordsPerBuffer / factor)
                    else:
                        break
            else:
                consistencyFlag = True
        if not consistencyFlag:
            raise MemoryError('The number of records and the number '
                'of samples in a single record is not consistent '
                'with the allowed maximum buffer size. Try to use '
                'values that could be factorized into small prime '
                'factors (values that are powers of 2 are most '
                'preferable).')
        buffersPerAcquisition = int(nRecord / recordsPerBuffer)
        # do not allocate more buffers than needed for all data
        bufferCount = int(min(2 * ((buffersPerAcquisition + 1) // 2),
                              maxBuffers))
        lT.append('Total buffers needed: %d' % buffersPerAcquisition)
        lT.append('Buffer count: %d' % bufferCount)
        lT.append('Buffer size [MB]: %f' % (float(bytesPerBuffer) / 2**20))
        lT.append('Records per buffer: %d' % recordsPerBuffer)

        # configure board, if wanted
        if bConfig:
            self.AlazarSetRecordSize(0, samplesPerRecord)
            self.AlazarSetRecordCount(nRecord)
            # allocate DMA buffers
            if bytesPerSample > 1:
                sample_type = ctypes.c_uint16
            else:
                sample_type = ctypes.c_uint8
            # clear old buffers
            self.removeBuffersDMA()
            # create new buffers
            self.buffers = []
            log.info('bufferCount %d' % bufferCount)
            for i in range(bufferCount):
                self.buffers.append(DMABuffer(sample_type,
                                              bytesPerBuffer))
        # arm and start capture, if wanted
        if bArm:
            # Configure the board to make a NPT AutoDMA acquisition
            self.AlazarBeforeAsyncRead(channels, 0,
                                       samplesPerRecord,
                                       recordsPerBuffer,
                                       nRecord,
                                       ADMA_EXTERNAL_STARTCAPTURE | ADMA_NPT)
            # post DMA buffers to board
            for buf in self.buffers:
                self.AlazarPostAsyncBuffer(buf.addr, buf.size_bytes)
            try:
                self.AlazarStartCapture()
            except BaseException:
                # make sure buffers release memory if failed
                self.removeBuffersDMA()
                raise

        # if not waiting for result, return here
        if not bMeasure:
            return

        lT.append('Post: %.1f ms' % (1000 * (time.clock() - t0)))
        
        # initialize data array
        avgTrace = np.zeros((numberOfChannels,
                             recordsPerBuffer,
                             samplesPerRecord), dtype=np.float32)
        if bGetAllTraces:
            bufferSize = recordsPerBuffer * numberOfChannels * samplesPerRecord
            samples = buffersPerAcquisition * bufferSize
            if bytesPerSample == 1:
                traces = np.empty(samples, dtype=np.uint8)
            elif bytesPerSample == 2:
                traces = np.empty(samples, dtype=np.uint16)
            else:
                traces = np.empty(samples, dtype=np.uint32)

            lT.append('Start: %.1f ms' % (1000 * (time.clock() - t0)))
            buffersCompleted = 0
            # range and zero for conversion to voltages
            codeZero = float(1 << (self.bitsPerSample - 1)) - .5
            codeRange = float(1 << (self.bitsPerSample - 1)) - .5
            # range and zero for each channel, combined with bit shifting
            rangeA = self.dRange[1] / codeRange
            rangeB = self.dRange[2] / codeRange

            timeout_ms = int(1000 * firstTimeout)
            log.info(str(lT))
            lT = []

        try:
            while (buffersCompleted < buffersPerAcquisition):
                # wait for the buffer at the head of the list of
                # available buffers to be filled by the board
                buf = self.buffers[buffersCompleted % len(self.buffers)]
                self.AlazarWaitAsyncBufferComplete(buf.addr,
                                                   timeout_ms=timeout_ms)

                # reset timeout time, can be different than first call
                timeout_ms = int(1000. * timeout)

                # break if stopped from outside
                if funcStop is not None and funcStop():
                    break
                # report progress
                if funcProgress is not None:
                    funcProgress(float(buffersCompleted) /
                                 float(buffersPerAcquisition))

                if bGetAllTraces:
                    bufferPosition = bufferSize * buffersCompleted
                    traces[bufferPosition:bufferPosition + bufferSize] = \
                        buf.buffer
                                 
                bufferRecords = buf.buffer.astype(np.float32, copy=False)
                bufferRecords.shape = (numberOfChannels,
                                       recordsPerBuffer,
                                       samplesPerRecord)
                bufferRecords -= codeZero
                if numberOfChannels == 2:
                    bufferRecords[0] *= rangeA
                    bufferRecords[1] *= rangeB
                elif bGetChA:
                    bufferRecords *= rangeA
                else:
                    bufferRecords *= rangeB
                avgTrace += bufferRecords

                buffersCompleted += 1
                
                # add the buffer to the end of the list of available
                # buffers
                self.AlazarPostAsyncBuffer(buf.addr, buf.size_bytes)
        finally:
            # release resources
            try:
                self.AlazarAbortAsyncRead()
            except BaseException:
                lT.append('Abort: %.1f ms' % (1000 * (time.clock() - t0)))
                pass

        # log timing information
        lT.append('Done: %.1f ms' % (1000 * (time.clock() - t0)))
        log.info(str(lT))
        
        if bGetAllTraces:
            tracesView = traces.view()
            tracesView.shape = (buffersPerAcquisition, numberOfChannels,
                                recordsPerBuffer, samplesPerRecord)
            tracesView = np.rollaxis(tracesView, 2, 1).reshape(nRecord,
                                numberOfChannels, samplesPerRecord)
            traces = tracesView.astype(np.float32, copy=False)

            traces -= codeZero
            if samplesPerRecord != nSamples:
                traces = traces[:,:,:nSamples]
            if numberOfChannels == 2:
                traces[:,0,:] *= rangeA
                traces[:,1,:] *= rangeB
                data['Channel A - Flattened Data'] = traces[:,0,:].flatten()
                data['Channel B - Flattened Data'] = traces[:,1,:].flatten()
            elif bGetChA:
                traces *= rangeA
                data['Channel A - Flattened Data'] = traces.flatten()
            else:
                traces *= rangeB
                data['Channel B - Flattened Data'] = traces.flatten()

        avgTrace /= float(buffersPerAcquisition)
        avgTraces = np.mean(avgTrace, axis=1)
        # return data - requested vector length, not restricted to
        # 64-bit multiple
        if samplesPerRecord != nSamples:
            avgTraces = avgTraces[:,:nSamples]
        if numberOfChannels == 2:
            data['Channel A - Averaged Data'] = avgTraces[0]
            data['Channel B - Averaged Data'] = avgTraces[1]
        elif bGetChA:
            data['Channel A - Averaged Data'] = avgTraces[0]
        else:
            data['Channel A - Averaged Data'] = avgTraces[1]
        return data

    def removeBuffersDMA(self):
        """Clear and remove DMA buffers to release memory."""
        # make sure buffers release memory
        for buf in self.buffers:
            buf.__exit__()
        # remove all
        self.buffers = []

    def readTracesSinglePort(self, Channel, bGetAllTraces):
        """Reads traces, converts to float, averages to a single trace."""
        # define sizes
        bytesPerSample = (self.bitsPerSample + 7) // 8
        samplesPerRecord = self.nPreSize + self.nPostSize
        dataBuffer = (c_uint8 * samplesPerRecord)()
        # define scale factors
        codeZero = float(1 << (self.bitsPerSample - 1)) - .5
        codeRange = float(1 << (self.bitsPerSample - 1)) - .5
        voltScale = self.dRange[Channel] / codeRange
        # initialize a scaled float vector
        avgTrace = np.zeros(samplesPerRecord, dtype=float)
        if bGetAllTraces:
            traces = np.empty((self.nRecord, samplesPerRecord), dtype=float)
        else:
            traces = None
        for n in range(self.nRecord):
            self.AlazarRead(Channel, dataBuffer, bytesPerSample, n + 1,
                            -self.nPreSize, samplesPerRecord)
            # convert and scale to float
            vBuffer = np.array(dataBuffer, dtype=float)
            vBuffer -= codeZero
            vBuffer *= voltScale
            if bGetAllTraces:
                traces[n] = vBuffer
            # add to output vector
            avgTrace += vBuffer
        # normalize
        avgTrace /= float(self.nRecord)
        return avgTrace, traces


if __name__ == '__main__':
    # test driver
    Digitizer = AlazarTechDigitizer()
