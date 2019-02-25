#!/usr/bin/env python

import InstrumentDriver
from VISA_Driver import VISA_Driver
import numpy as np


class Driver(VISA_Driver):
    """This class implements the Tektronix AWG driver."""
    def _output(self, status=True):
        sOutput = ''
        for n, bUpdate in enumerate(self.lInUse):
            if bUpdate:
                sOutput += (':OUTP%d:STAT %d;' % ((n + 1), int(status)))
        if sOutput != '':
            self.writeAndLog(sOutput)

    def _run(self):
        if self.bIsStopped:
            self.askAndLog('*OPC?')
            self.writeAndLog(':AWGC:RUN')
            self.bIsStopped = False
            
    def _stop(self):
        if not self.bIsStopped:
            self.writeAndLog(':AWGC:STOP')
            self.bIsStopped = True

    def _error(self):
        self.askAndLog('*OPC?')
        try:
            esr = int(self.askAndLog('*ESR?'))
        except:
            sBuffer = self.read()
            self.log('Extra data read from Tektronix AWG: %s'% sBuffer)
            esr = int(self.askAndLog('*ESR?'))
        if esr & 0b0001000:
            stb = int(self.askAndLog('*STB?'))
            # self.writeAndLog('*RST')
            raise InstrumentDriver.Error('Make sure that what '
                    'you are trying to do is compatible with '
                    'AWG specifications. The Standard Event '
                    'Status Register (SESR) value is %s '
                    'and the Status Byte Register (SBR) value '
                    'is %s.' % (bin(esr), bin(stb)))
                    
    def _clear_all(self):
        self.bWaveUpdated = False
        self.lOldU16 = [[np.array([], dtype=np.uint16)
                       for _ in range(self.nCh)]]
        self.lInUse = [False] * self.nCh
        self._clear()

    def _clear(self):
        self.nOldSeq = -1
        self.lStoredNames = [] # used to record the name of the waveforms 
        # that will be used to construct the sequence in fast sequence mode
        self.lStoredSubseq = []
        self.lValidChannel = []
        self.lWaveform = [] # used to record which waveforms have been
        # created on the AWG in fast sequence mode
        self.dSubseqList = {}
        self.iPulseIdx = 0

    def _getTrigChannel(self, options):
        """Helper function, get trig channel for instrument, or None if N/A"""
        trig_channel = options.get('trig_channel', None)
        return trig_channel

    def performOpen(self, options={}):
        """Perform the operation of opening the instrument connection"""
        # add compatibility with pre-Python 3 version of Labber
        if not hasattr(self, 'write_raw'):
            self.write_raw = self.write
        # add compatibility with pre-1.5.4 version of Labber
        if not hasattr(self, 'getTrigChannel'):
            self.getTrigChannel = self._getTrigChannel
        # start by calling the generic VISA open to make sure we have
        # a connection
        VISA_Driver.performOpen(self, options)
        # check for a strange bug by reading the status bit
        try:
            status = self.askAndLog('*STB?', bCheckError=False)
            status = int(status)
        except:
            # if conversion to int failed, re-read instrument buffer to
            # clear
            sBuffer = self.read()
            self.log('Extra data read from Tektronix AWG: %s' % sBuffer)
        # get model name and number of channels
        sModel = self.getModel()
        self.nCh = 4 if sModel in ('5004', '5014') else 2
        # turn off run mode
        self.bIsStopped = False
        self._stop()
        self.writeAndLog('WLIS:WAV:DEL ALL', bCheckError=False)
        self.writeAndLog('SLIS:SUBS:DEL ALL', bCheckError=False)
        self.bFastSeq = False
        self.bSubSeq = False
        self._clear_all()
        for n in range(self.nCh):
            self.createWaveformOnTek(n + 1, 0, bOnlyClear=True)

    def performClose(self, bError=False, options={}):
        """Perform the close instrument connection operation"""
        # close VISA connection
        VISA_Driver.performClose(self, bError, options)

    def initSetConfig(self):
        """This function is run before setting values in Set Config"""
        # turn off run mode
        self._stop()
        self._clear_all()
        for n in range(self.nCh):
            channel = n + 1
            self.setValue('Ch %d' % channel, [])
            self.setValue('Ch %d - Marker 1' % channel, [])
            self.setValue('Ch %d - Marker 2' % channel, [])
            self.createWaveformOnTek(channel, 0, bOnlyClear=True)

    def performSetValue(self, quant, value, sweepRate=0.0, options={}):
        """
        Perform the Set Value instrument operation. This method should
        return the actual value set by the instrument.
        """
        # keep track of if waveform is updated, to avoid sending
        # it many times
        if self.isFirstCall(options):
            for n in range(self.nCh):
                channel = n + 1
                self.setValue('Ch %d' % channel, [])
                self.setValue('Ch %d - Marker 1' % channel, [])
                self.setValue('Ch %d - Marker 2' % channel, [])
            self.bWaveUpdated = False
            # if sequence mode, make sure the buffer contains enough
            # waveforms
            if self.isHardwareLoop(options):
                (seq_no, n_seq) = self.getHardwareLoopIndex(options)
                mode = self.getValue('Sequence transfer mode')
                if mode == 'Subwaveforms':
                    self.bFastSeq = True
                    self.bSubSeq = False
                elif mode == 'Subsequences':
                    self.bFastSeq = True
                    self.bSubSeq = True
                else:
                    self.bFastSeq = False
                    self.bSubSeq = False
                # if first call, clear sequence and create buffer
                if seq_no == 0:
                    self._stop()
                    self.writeAndLog('WLIS:WAV:DEL ALL', bCheckError=False)
                    self.writeAndLog('SLIS:SUBS:DEL ALL', bCheckError=False)
                    self._clear()
                # if different sequence length, re-create buffer
                if seq_no == 0 and n_seq != len(self.lOldU16):
                    self.lOldU16 = [[np.array([], dtype=np.uint16)
                            for n1 in range(self.nCh)]
                            for n2 in range(n_seq)]
            elif self.isHardwareTrig(options):
                # if hardware triggered, always stop outputting before
                # setting
                self._stop()
                self.bFastSeq = False
                self.bSubSeq = False
            else:
                self.bFastSeq = False
                self.bSubSeq = False
            
        if quant.name in ('Ch 1', 'Ch 2', 'Ch 3', 'Ch 4',
                          'Ch 1 - Marker 1', 'Ch 1 - Marker 2',
                          'Ch 2 - Marker 1', 'Ch 2 - Marker 2',
                          'Ch 3 - Marker 1', 'Ch 3 - Marker 2',
                          'Ch 4 - Marker 1', 'Ch 4 - Marker 2'):
            # set value, then mark that waveform needs an update
            quant.setValue(value)
            self.bWaveUpdated = True
        elif quant.name in ('Run'):
            if value:
                # turn on channels again, to avoid issues when switch run mode
                sOutput = ''
                for n, bUpdate in enumerate(self.lInUse):
                    if bUpdate:
                        sOutput += (':OUTP%d:STAT 1;' % (n + 1))
                if sOutput != '':
                    self.writeAndLog(sOutput)
                self._run()
                self._error()
            else:
                # stop AWG
                self._stop()
        elif quant.name == 'Hardware loop forced stop':
            self._stop()
            value = True
        elif quant.name in ('Sequence transfer mode', 'Run Mode'):
            quant.setValue(value)
        else:
            # for all other cases, call VISA driver
            value = VISA_Driver.performSetValue(self, quant, value,
                    sweepRate, options=options)

        # if final call and wave is updated, send it to AWG
        if self.isFinalCall(options) and self.bWaveUpdated:
            (seq_no, n_seq) = self.getHardwareLoopIndex(options)
            mode = self.getValue('Run mode')
            if self.isHardwareLoop(options):
                if mode == 'Continuous':
                    self.setValue('Run mode', 'Sequence')
                seq = seq_no
                self.reportStatus('Sending waveform (%d out of %d)'
                        % (seq_no + 1, n_seq))
            else:
                if mode == 'Sequence':
                    self.setValue('Run mode', 'Continuous')
                seq = None
                if len(self.lOldU16) > 1:
                    self.lOldU16 = [[np.array([], dtype=np.uint16)
                            for n1 in range(self.nCh)]]
            
            # in trig mode, don't start AWG if trig channel will start
            # it later
            if (self.isHardwareTrig(options) and
                    self.getTrigChannel(options) == 'Run'):
                bStart = False
            else:
                bStart = True
            self.sendWaveformAndStartTek(seq=seq, n_seq=n_seq,
                    bStart=bStart)
        return value

    def performGetValue(self, quant, options={}):
        """Perform the Get Value instrument operation"""
        if quant.name in ('Ch 1', 'Ch 2', 'Ch 3', 'Ch 4',
                          'Ch 1 - Marker 1', 'Ch 1 - Marker 2',
                          'Ch 2 - Marker 1', 'Ch 2 - Marker 2',
                          'Ch 3 - Marker 1', 'Ch 3 - Marker 2',
                          'Ch 4 - Marker 1', 'Ch 4 - Marker 2',
                          'Sequence transfer mode'):
            # do nothing here
            value = quant.getValue()
        else:
            # for all other cases, call VISA driver
            value = VISA_Driver.performGetValue(self, quant, options)
        return value

    def sendWaveformAndStartTek(self, seq=None, n_seq=1, bStart=True):
        """Rescale and send waveform data to the Tek"""
        # get model name and number of channels
        self.nPrevData = 0
        bWaveUpdate = False
        if self.bFastSeq:
            if seq == 0:
                # first waveform in the sequence
                if len(self.lStoredNames) != n_seq:
                    # recreate lStoredNames
                    self.lStoredNames = [[]] * n_seq
                if len(self.lStoredSubseq) != n_seq:
                    self.lStoredSubseq = [[]] * n_seq
            self.decomposeWaveformAndSendFragments(seq)
        else:
            # go through all channels
            for n in range(self.nCh):
                # channels are numbered 1-4
                channel = n + 1
                vData = self.getValueArray('Ch %d' % channel)
                vMark1 = self.getValueArray('Ch %d - Marker 1' % channel)
                vMark2 = self.getValueArray('Ch %d - Marker 2' % channel)
                bWaveUpdate = self.sendWaveformToTek(channel, vData,
                        vMark1, vMark2, seq) or bWaveUpdate      
        # check if sequence mode
        if seq is not None:
            # if not final seq call, just return here
            if seq + 1 < n_seq:
                return
            # final call, check if sequence has changed
            if bWaveUpdate or n_seq != self.nOldSeq:
                # create sequence list, first clear to reset old values
                if self.bFastSeq:
                    # notice that the list of the names has already been
                    # updated
                    self.assembleSequence()
                else:
                    self.writeAndLog('SEQ:LENG 0')
                    self.writeAndLog('SEQ:LENG %d' % n_seq)
                    for n1 in range(n_seq):               
                        for n2, bUpdate in enumerate(self.lInUse):
                            if bUpdate:
                                name = 'wvfrm_ch%d_%05d' % (n2 + 1, n1 + 1)
                                self.writeAndLog('SEQ:ELEM%d:WAV%d "%s"'
                                        % (n1 + 1, n2 + 1, name))
                        # don't wait for trigger 
                        self.writeAndLog('SEQ:ELEM%d:TWA 0' % (n1 + 1))
                    # for last element, set jump to first
                    self.writeAndLog('SEQ:ELEM%d:GOTO:STAT 1' % n_seq)
                    self.writeAndLog('SEQ:ELEM%d:GOTO:IND 1' % n_seq)
                # save old sequence length
                self.nOldSeq = n_seq

            self.askAndLog('*OPC?')
            # turn on sequence mode
            self.writeAndLog(':AWGC:RMOD SEQ')
            # turn on channels in use
            self._output(True)
            return
        # turn on channels in use 
        self._output(True)
        # if not starting, make sure AWG is not running, then return
        if not bStart:
            self.askAndLog('*OPC?')
            iRunState = int(self.askAndLog(':AWGC:RST?'))
            nTry = 1000
            while nTry and iRunState and not self.isStopped():
                # sleep for while to save resources, then try again
                self.wait(0.05)
                # try again
                iRunState = int(self.askAndLog(':AWGC:RST?'))
                nTry -= 1
            return
        # check if AWG has been stopped, if not return here
        if not self.bIsStopped:
            # no waveforms updated, just turn on output, no need to wait
            # for run
            return
        # send command to turn on run mode to tek
        self._run()
        # wait for output to be turned on again
        iRunState = int(self.askAndLog(':AWGC:RST?'))
        nTry = 1000
        while nTry and iRunState == 0 and not self.isStopped():
            # sleep for while to save resources, then try again
            self.wait(0.05)
            # try again
            iRunState = int(self.askAndLog(':AWGC:RST?'))
            nTry -= 1
        # check if timeout occurred
        if nTry <= 0:
            # timeout
            raise InstrumentDriver.Error('Cannot turn on Run mode.')
        # turn on channels again, to avoid issues when turning on/off
        # run mode
        self._output(True)

    def scaleWaveformToU16(self, vData, dVpp, ch):
        """Scales the waveform and returns data in a string of U16"""
        # make sure waveform data is within the voltage range
        if np.any(vData > 1.003 * dVpp / 2.) or \
                np.any(vData < -1.003 * dVpp / 2.):
            raise InstrumentDriver.Error(
                    'Waveform for channel %d contains values that are ' 
                    'outside the channel voltage range.' % ch)
        # clip waveform and store in-place
        np.clip(vData, -dVpp / 2., dVpp / 2., vData)
        return np.array(16382 * (vData + dVpp / 2.) / dVpp, dtype=np.uint16)
        
    def createWaveformOnTek(self, channel, length, seq=None, bOnlyClear=False):
        """
        Remove old and create new waveform on the Tek. The waveform
        is named by the channel number.
        """
        if seq is None:
            name = 'wvfrm_ch%d' % channel
        else:
            name = 'wvfrm_ch%d_%05d' % (channel, seq+1)
        # first, turn off output
        self.writeAndLog(':OUTP%d:STAT 0' % channel)
        if bOnlyClear:
            # just clear this channel (from manual, only when Run Mode
            # is not Sequence?)
            self.writeAndLog(':SOUR%d:WAV ""' % channel)
        else:
            # remove old waveform, ignoring errors, then create new
            self.writeAndLog(':WLIS:WAV:DEL "%s";*CLS' % name,
                    bCheckError=False)
            self.writeAndLog(':WLIS:WAV:NEW "%s",%d,INT' % (name, length))
 
    def sendWaveformToTek(self, channel, vData, vMark1, vMark2, seq=None):
        """Send a waveform to Tek."""
        # check if sequence
        if seq is None:
            iSeq = 0
        else:
            iSeq = seq
        # channels are named 1-4
        n = channel - 1
        if len(vData) == 0:
            if len(vMark1) == 0 and len(vMark2) == 0:
                # if channel in use, turn off, clear, go to next channel
                if self.lInUse[n]:
                    self.createWaveformOnTek(channel, 0, seq, bOnlyClear=True)
                    self.lOldU16[iSeq][n] = np.array([], dtype=np.uint16)
                    self.lInUse[n] = False
                return False
            else:
                # no data, but markers exist, output zeros for data
                nMark = max(len(vMark1), len(vMark2))
                vData = np.zeros((nMark,), dtype=float)
        # make sure length of data is the same
        if (len(vMark1) and len(vData) != len(vMark1)) or \
           (len(vMark2) and len(vData) != len(vMark2)) or \
           (self.nPrevData and len(vData) != self.nPrevData):        
            raise InstrumentDriver.Error(
                    'All channels need to have the same number of '
                    'elements: Marker 1 length is %d, '
                              'Marker 2 length is %d, '
                              'Analog Waveform length is %d.'
                              'nPrevData is %d'
                              % (len(vMark1), len(vMark2),
                                 len(vData), self.nPrevData))
        self.nPrevData = len(vData)
        # channel in use, mark
        self.lInUse[n] = True
        # get range and scale to U16
        Vpp = self.getValue('Ch%d - Range' % channel)
        vU16 = self.scaleWaveformToU16(vData, Vpp, channel)
        # check for marker traces
        for m, marker in enumerate([vMark1, vMark2]):
            if len(marker) == len(vU16):
                # get marker trace
                vMU16 = np.array(marker != 0, dtype=np.uint16)
                # add marker trace to data trace, with bit shift
                vU16 += 2**(14 + m) * vMU16
        start, length = 0, len(vU16)
        # compare to previous trace
        if length != len(self.lOldU16[iSeq][n]):
            # stop AWG if still running
            if not self.bIsStopped:
                self._stop()
            # length has changed, delete old waveform and create new
            self.createWaveformOnTek(channel, length, seq)
        else:
            # same length, check for similarities
            vIndx = np.nonzero(vU16 != self.lOldU16[iSeq][n])[0]
            if len(vIndx) == 0:
                # nothing changed, don't update, go on to next
                return False
            # some elements changed, find start and length
            start = vIndx[0]
            length = vIndx[-1] - vIndx[0] + 1
        # stop AWG if still running
        self._stop()
        # create binary data as bytes with header
        sLen = b'%d' % (2 * length)
        sHead = b'#%d%s' % (len(sLen), sLen)
        # send to tek, start by turning off output
        if seq is None:
            # non-sequence mode, get name
            name = b'wvfrm_ch%d' % channel
        else:
            # sequence mode, get name
            name = b'wvfrm_ch%d_%05d' % (channel, seq+1)
        sSend = b':OUTP%d:STAT 0;' % channel
        sSend += b':WLIS:WAV:DATA "%s",%d,%d,%s' % (name, start, length,
                 sHead + vU16[start:start+length].tobytes())
        self.write_raw(sSend)
        if seq is None:
            # (re-)set waveform to channel
            self.writeAndLog(':SOUR%d:WAV "%s"' % (channel, name.decode()))
        # store new waveform for next call
        self.lOldU16[iSeq][n] = vU16
        return True

    def checkWaveformLengthBeforeSending(self, channel, vData,
            vMark1, vMark2, seq):
        """
        Check waveform before sending to Tek. Similar to the method
        self.sendWaveformToTek(), but not the same. This is for fast 
        sequence transfer.
        """
        # check if sequence
        if seq is None:
            iSeq = 0
        else:
            iSeq = seq
        # channels are named 1-4
        n = channel - 1
        if len(vData) == 0:
            if len(vMark1) == 0 and len(vMark2) == 0:
                # if channel in use, turn off, clear, go to next channel
                if self.lInUse[n]:
                    self.createWaveformOnTek(channel, 0, seq,
                            bOnlyClear=True)
                    self.lOldU16[iSeq][n] = np.array([], dtype=np.uint16)
                    self.lInUse[n] = False
                return (False, [])
            else:
                # no data, but markers exist, output zeros for data
                nMark = max(len(vMark1), len(vMark2))
                vData = np.zeros((nMark,), dtype=float)                
        else:
            nMark = len(vData)
            
        if len(vMark1) == 0:
            vMark1 = np.zeros((nMark,), dtype=float)
        if len(vMark2) == 0:
            vMark2 = np.zeros((nMark,), dtype=float)  
        mDataMark = np.array([vData, vMark1, vMark2]) 

        # make sure length of data is the same
        if (len(vMark1) and len(vData) != len(vMark1)) or \
           (len(vMark2) and len(vData) != len(vMark2)) or \
           (self.nPrevData and len(vData) != self.nPrevData):        
            raise InstrumentDriver.Error(
                    'All channels need to have the same number of '
                    'elements: Marker 1 length is %d, '
                              'Marker 2 length is %d, '
                              'Analog Waveform length is %d.'
                              %(len(vMark1), len(vMark2), len(vData)))
        self.nPrevData = len(vData)
        # channel in use, mark
        self.lInUse[n] = True
        return (True, mDataMark)

    def scaleWaveformAndMarkerToU16(self, vData, vMark1, vMark2, channel):
        # get range and scale to U16
        Vpp = self.getValue('Ch%d - Range' % channel)
        vU16 = self.scaleWaveformToU16(vData, Vpp, channel)
        # check for marker traces
        for m, marker in enumerate([vMark1, vMark2]):
            if len(marker) == len(vU16):
                # get marker trace
                vMU16 = np.array(marker != 0, dtype=np.uint16)
                # add marker trace to data trace, with bit shift
                vU16 += 2**(14 + m) * vMU16
        return vU16

    def _waveformName(self, mDataMark, channel, idx):
        unique = np.unique(mDataMark)
        length = mDataMark.shape[1]
        if unique.size == 1:
            # constant pulse
            if unique[0] == 0:
                # all zeros
                name = 'zeros_%d' % length
        else:
            uData = np.unique(mDataMark[0,:])
            uMark1 = np.unique(mDataMark[1,:])
            uMark2 = np.unique(mDataMark[2,:])
            if uData.size == 1 and uMark1.size == 1 and uMark2.size == 1:
                name = 'const_%d_%s_%s_%s' % (length,
                    str(np.float(uData[0])),
                    str(np.float(uMark1[0])),
                    str(np.float(uMark2[0])))
                name = name.replace('-', 'm').replace('.', 'p')
            else:
                name = 'pulse_%d_ch%d_%05d' % (length, channel, idx)
        return name
    
    def decomposeWaveformAndSendFragments(self, seq):
        """Decompose the waveform into fragments to assemble a sequence."""
        mDataMarkTot = 0 # initial value
        lValidChannel = []
        for n in range(self.nCh):
            # channels are numbered 1-4
            channel = n + 1
            vData = self.getValueArray('Ch %d' % channel)
            vMark1 = self.getValueArray('Ch %d - Marker 1' % channel)
            vMark2 = self.getValueArray('Ch %d - Marker 2' % channel)
            
            # similar process in self.sendWaveformToTek
            bValidChannel, mDataMark = \
                self.checkWaveformLengthBeforeSending(channel, vData,
                        vMark1, vMark2, seq)

            if bValidChannel:
                # this channel is used
                lValidChannel.append(channel)
                if isinstance(mDataMarkTot, int):
                    # if this is the first valid channel
                    mDataMarkTot = mDataMark
                else:
                    mDataMarkTot = np.concatenate((mDataMarkTot, mDataMark))

        # chunk waveform
        Len = mDataMarkTot.shape[1]
        if Len > 25000:
            FragmentSize = int(Len / 25)
        elif Len > 5000:
            FragmentSize = 1000
        else:
            raise InstrumentDriver.Error('The length of each sequence '
                    'is too short! Do not use Subwaveforms or '
                    'Subsequences transfer modes.')

        vNode = np.arange(0, Len, FragmentSize)
        if Len % FragmentSize:
            vNode = vNode[:-1]
        numFrag = len(vNode)
        vNode = np.concatenate((vNode, [Len]))
        lNames = [[]] * len(lValidChannel)
        for nC in range(len(lValidChannel)):
            channel = lValidChannel[nC]
            for nF in range(numFrag):
                start = vNode[nF]
                end = vNode[nF+1]
                # pick up from the 'start' element to 'end - 1' element
                thisPulse = mDataMarkTot[3*nC:3*nC+3,start:end]
                name = self._waveformName(thisPulse, channel, self.iPulseIdx)
                if name not in self.lWaveform:
                    # send waveform
                    self.createAndSendFragment(thisPulse, channel, name)
                    self.lWaveform = self.lWaveform + [name]
                    if name.startswith('pulse'):
                        self.iPulseIdx += 1
                # store name
                lNames[nC] = lNames[nC] + [name]

        # create subsequence
        mNames = np.array(lNames)
        lNames = [[i] for i in mNames[:,0]]
        ThisSubseq = mNames[:,0]

        SubseqName = 'subseq_%s' % str(hash(ThisSubseq.tostring()))
        # the name exists for this subsequence;
        # check whether the subsequences are the same.        
        if SubseqName in self.dSubseqList and \
                np.any(self.dSubseqList[SubseqName] != ThisSubseq):
            k = 1
            testname = '%s_%d' % (SubseqName, k)
            while testname in self.dSubseqList and \
                    np.any(self.dSubseqList[SubseqName] != ThisSubseq):
                k += 1
                testname = '%s_%d' % (SubseqName, k)
            SubseqName = testname
         
        if SubseqName not in self.dSubseqList:
            # does not exist; create then store it
            self.dSubseqList[SubseqName] = ThisSubseq
            if self.bSubSeq:
                self.writeAndLog(':SLIS:SUBS:DEL "%s";*CLS'
                        % SubseqName, bCheckError=False)
                self.writeAndLog(':SLIS:SUBS:NEW "%s",%d' % (SubseqName, 1))
                for k, channel in enumerate(lValidChannel):                    
                    self.writeAndLog(':SLIS:SUBS:ELEM%d:WAV%d "%s","%s";*WAI'
                            % (1, channel, SubseqName, ThisSubseq[k]))

        lSubseq = [SubseqName]
        for nF in range(numFrag - 1):
            name = mNames[0,nF+1]
            # check whether two subsequences are the same
            if np.all(mNames[:,nF] == mNames[:,nF+1]):
                LastName = lNames[0][-1]
                if LastName.startswith('R'):
                    # find repeat number
                    ind = 1
                    while LastName[ind].isnumeric():
                        ind += 1
                    LoopNum = int(LastName[1:ind]) + 1
                    lNames[0][-1] = 'R%d%s' % (LoopNum, name)
                else:
                    lNames[0][-1] = 'R2%s' % name
                    
                LastSubseq = lSubseq[-1]
                if LastSubseq.startswith('R'):
                    # find repeat number
                    ind = 1
                    while LastSubseq[ind].isnumeric():
                        ind += 1
                    LoopNum = int(LastSubseq[1:ind]) + 1
                    name = LastSubseq[ind:]
                    lSubseq[-1] = 'R%d%s' % (LoopNum, name)
                else:
                    lSubseq[-1] = 'R2%s' % LastSubseq
            else:
                ThisSubseq = mNames[:,nF+1]
                SubseqName = 'subseq_%s' % str(hash(ThisSubseq.tostring()))
                for n2, channel in enumerate(lValidChannel):
                    name = ThisSubseq[n2]
                    lNames[n2] = lNames[n2] + [name]
                if SubseqName in self.dSubseqList and \
                        np.any(self.dSubseqList[SubseqName] != ThisSubseq):
                    k = 1
                    testname = '%s_%d' % (SubseqName, k)
                    while testname in self.dSubseqList and \
                            np.any(self.dSubseqList[SubseqName] != ThisSubseq):
                        k += 1
                        testname = '%s_%d' % (SubseqName, k)
                    SubseqName = testname
                 
                if SubseqName not in self.dSubseqList:
                    # does not exist; create then store it
                    self.dSubseqList[SubseqName] = ThisSubseq
                    if self.bSubSeq:
                        self.writeAndLog(':SLIS:SUBS:DEL "%s";*CLS'
                            % SubseqName, bCheckError=False)
                        self.writeAndLog(':SLIS:SUBS:NEW "%s",%d' % (SubseqName, 1))
                        for k, channel in enumerate(lValidChannel):
                            self.writeAndLog(':SLIS:SUBS:ELEM%d:WAV%d "%s","%s";*WAI'
                                    % (1, channel, SubseqName, ThisSubseq[k]))

                lSubseq = lSubseq + [SubseqName]

        # transpose the list, so that the sublists are divided by 
        # different fragments (therefore different sublists are
        # different elements in the whole sequence)                    
        lNames = np.array(lNames).transpose().tolist()
        # compare to previous trace
        if lValidChannel != self.lValidChannel or \
                lNames != self.lStoredNames[seq] or \
                lSubseq != self.lStoredSubseq[seq]:
            self.lValidChannel = lValidChannel
            if len(lNames) > 1:
                if len(lNames[0]) != len(lNames[1]):
                    raise InstrumentDriver.Error('Some fragments are missing.')
            self.lStoredSubseq[seq] = lSubseq
            self.lStoredNames[seq] = lNames

    def createAndSendFragment(self, mDataMark, channel, name):
        """Determine the name of the pulse and send it to the AWG."""
        vData = mDataMark[0,:]
        vMark1 = mDataMark[1,:]
        vMark2 = mDataMark[2,:]
        
        self._stop()
        # create binary data as bytes with header
        vU16 = self.scaleWaveformAndMarkerToU16(vData,
                vMark1, vMark2, channel)

        start, length = 0, len(vU16)
        sLen = b'%d' % (2 * length)
        sHead = b'#%d%s' % (len(sLen), sLen)
        # send to tek, start by turning off output
        self.writeAndLog(':OUTP%d:STAT 0' % channel)
        # remove old waveform, ignoring errors, then create new
        self.writeAndLog(':WLIS:WAV:DEL "%s";*CLS' % name,
                bCheckError=False)
        self.askAndLog('*OPC?')
        self.writeAndLog(':WLIS:WAV:NEW "%s",%d,INT' % (name, length))
        sSend = b':WLIS:WAV:DATA "%s",%d,%d,%s;*WAI' % (name.encode(),
                start, length, sHead + vU16[start:start+length].tobytes())
        self.write_raw(sSend)

    def assembleSequence(self):
        """Assemble the fragment waveforms into a sequence."""
        self.writeAndLog('SEQ:LENG 0')
        lFlattenNames = [item for sublist in self.lStoredNames
                         for item in sublist]
        lFlattenSubseq = [item for sublist in self.lStoredSubseq
                         for item in sublist]
        if self.bSubSeq:
            seq_len = len(lFlattenSubseq)
        else:
            seq_len = len(lFlattenNames)
        # flatten the list, merge the lists for different seq together,
        # which means:
        #    for sublist in l:
        #        for item in sublist:
        #            lFlattenNames.append(item)
        self.reportStatus('Assembling sequence')
        self.writeAndLog('SEQ:LENG %d' % seq_len)
        for n1 in range(seq_len):
            if self.bSubSeq:
                name = lFlattenSubseq[n1]
                loopnum = 1
                if name.startswith('R'):
                    ind = 1
                    while name[ind].isnumeric():
                        ind += 1
                    loopnum = int(name[1:ind])
                    name = name[ind:]
                sSend = ':SEQ:ELEM%d:SUBS "%s";' % (n1 + 1, name)
                if loopnum > 1:
                    sSend += (':SEQ:ELEM%d:LOOP:COUNT %d;'
                            % (n1 + 1, loopnum))
            else:
                sSend = ''
                for n2, channel in enumerate(self.lValidChannel):
                    name = lFlattenNames[n1][n2]
                    loopnum = 1
                    if name.startswith('R'):
                        ind = 1
                        while name[ind].isnumeric():
                            ind += 1
                        loopnum = int(name[1:ind])
                        name = name[ind:]
                    sSend += (':SEQ:ELEM%d:WAV%d "%s";'
                            % (n1 + 1, channel, name))
                    if loopnum > 1 and n2 == 0:
                        sSend += (':SEQ:ELEM%d:LOOP:COUNT %d;'
                                % (n1 + 1, loopnum))
        
            sSend += ':SEQ:ELEM%d:TWA 0;*WAI' % (n1 + 1)
            self.write_raw(sSend.encode())

        # for last element, set jump to first
        self.writeAndLog('SEQ:ELEM%d:GOTO:STAT 1' % seq_len)
        self.writeAndLog('SEQ:ELEM%d:GOTO:IND 1' % seq_len)


if __name__ == '__main__':
    pass
