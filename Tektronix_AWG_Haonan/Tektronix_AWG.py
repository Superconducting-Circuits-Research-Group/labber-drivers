#!/usr/bin/env python

import InstrumentDriver
from VISA_Driver import VISA_Driver
import numpy as np

class Driver(VISA_Driver):
    """ This class implements the Tektronix AWG5014 driver"""
    
    def performOpen(self, options={}):
        """Perform the operation of opening the instrument connection"""
        # add compatibility with pre-python 3 version of Labber
        if not hasattr(self, 'write_raw'):
            self.write_raw = self.write
        # start by calling the generic VISA open to make sure we have a connection
        VISA_Driver.performOpen(self, options)
        # check for strange bug by reading the status bit
        try:
            status = self.askAndLog('*STB?', bCheckError=False)
            status = int(status)
        except:
            # if conversion to int failed, re-read instrument buffer to clear
            sBuffer = self.read()
            self.log('Extra data read from Tek: %s, %s' % (str(status), sBuffer))
        # get model name and number of channels
        sModel = self.getModel()
        self.nCh = 4 if sModel in ('5004', '5014') else 2
        # turn off run mode
        self.writeAndLog(':AWGC:STOP;')
        # init vectors with old values
        self.bWaveUpdated = False
        self.nOldSeq = -1
        self.nOldLength = 0
        self.lOldU16 = [[np.array([], dtype=np.uint16) \
                       for n1 in range(self.nCh)] for n2 in range(1)]
        self.lStoredNames = []
        self.lValidChannel = []
        self.lWaveform = []
        # clear old waveforms
        self.lInUse = [False]*self.nCh
        for n in range(self.nCh):
            self.createWaveformOnTek(n+1, 0, bOnlyClear=True)


    def performClose(self, bError=False, options={}):
        """Perform the close instrument connection operation"""
        # close VISA connection
        VISA_Driver.performClose(self, bError, options)


    def initSetConfig(self):
        """This function is run before setting values in Set Config"""
        # turn off run mode
        self.writeAndLog(':AWGC:STOP;')
        # init vectors with old values
        self.bWaveUpdated = False
        self.bLengthUpdate = True
        self.bFastSeq = False
        self.nOldLength = 0
        self.nOldSeq = -1
        self.lOldU16 = [[np.array([], dtype=np.uint16) \
                       for n1 in range(self.nCh)] for n2 in range(1)]
        # clear old waveforms
        self.lStoredNames = []
        self.lValidChannel = []
        self.lWaveform = []
        self.lInUse = [False]*self.nCh
        for n in range(self.nCh):
            channel = n+1
            self.setValue('Ch %d' % channel, [])
            self.setValue('Ch %d - Marker 1' % channel, [])
            self.setValue('Ch %d - Marker 2' % channel, [])
            self.createWaveformOnTek(channel, 0, bOnlyClear=True)


    def performSetValue(self, quant, value, sweepRate=0.0, options={}):
        """Perform the Set Value instrument operation. This function should
        return the actual value set by the instrument"""
        # keep track of if waveform is updated, to avoid sending it many times
        if self.isFirstCall(options):
            self.bWaveUpdated = False
            self.bLengthUpdate = True
            # if sequence mode, make sure the buffer contains enough waveforms
            if self.isHardwareLoop(options):
                (seq_no, n_seq) = self.getHardwareLoopIndex(options)
                self.bFastSeq = self.getValue('Fast Sequence Transfer')
                # if first call, clear sequence and create buffer
                if seq_no==0:
                    # variable for keepin track of sequence updating
                    self.writeAndLog(':AWGC:STOP;')
                    self.bSeqUpdate = False
                # if different sequence length, re-create buffer
                if seq_no==0 and n_seq != len(self.lOldU16):
                    self.lOldU16 = [[np.array([], dtype=np.uint16) \
                                   for n1 in range(self.nCh)] for n2 in range(n_seq)]
            elif self.isHardwareTrig(options):
                # if hardware triggered, always stop outputting before setting
                self.writeAndLog(':AWGC:STOP;')
        if quant.name in ('Ch 1', 'Ch 2', 'Ch 3', 'Ch 4',
                          'Ch 1 - Marker 1', 'Ch 1 - Marker 2',
                          'Ch 2 - Marker 1', 'Ch 2 - Marker 2',
                          'Ch 3 - Marker 1', 'Ch 3 - Marker 2',
                          'Ch 4 - Marker 1', 'Ch 4 - Marker 2'):
            # set value, then mark that waveform needs an update
            # self.log('***quant.name=%s, len(value["y"])=%d****' % (quant.name, len(value['y'])))
            # if the length of the waveform changes, clear all the waveforms and reset them
            
            if self.bLengthUpdate and len(value['y']) != self.nOldLength:
                for n in range(self.nCh):
                    channel = n+1
                    self.setValue('Ch %d' % channel, [])
                    self.setValue('Ch %d - Marker 1' % channel, [])
                    self.setValue('Ch %d - Marker 2' % channel, [])
                    # self.createWaveformOnTek(channel, 0, bOnlyClear=True)
                # self.log('======================all waveforms clear===========================')
                self.bLengthUpdate = False
            
            quant.setValue(value)
            self.bWaveUpdated = True
        elif quant.name in ('Run'):
            self.writeAndLog(':AWGC:RUN')
            # turn on channels again, to avoid issues when turning on/off run mode
            sOutput = ''
            for n, bUpdate in enumerate(self.lInUse):
                if bUpdate:
                    sOutput += (':OUTP%d:STAT 1;' % (n+1))
            if sOutput != '':
                self.writeAndLog(sOutput)
        elif quant.name == 'Fast Sequence Transfer':
            quant.setValue(value)
        else:
            # for all other cases, call VISA driver
            value = VISA_Driver.performSetValue(self, quant, value, sweepRate,
                                                options=options)
        # if final call and wave is updated, send it to AWG
        if self.isFinalCall(options) and self.bWaveUpdated:
            (seq_no, n_seq) = self.getHardwareLoopIndex(options)
            if self.isHardwareLoop(options):
                seq = seq_no
                self.reportStatus('Sending waveform (%d/%d)' % (seq_no+1, n_seq))
            else:
                seq = None
            bStart = not self.isHardwareTrig(options)
            self.sendWaveformAndStartTek(seq=seq, n_seq=n_seq, bStart=bStart)
        return value


    def performGetValue(self, quant, options={}):
        """Perform the Get Value instrument operation"""
        # check type of quantity
        if quant.name in ('Ch 1', 'Ch 2', 'Ch 3', 'Ch 4',
                          'Ch 1 - Marker 1', 'Ch 1 - Marker 2',
                          'Ch 2 - Marker 1', 'Ch 2 - Marker 2',
                          'Ch 3 - Marker 1', 'Ch 3 - Marker 2',
                          'Ch 4 - Marker 1', 'Ch 4 - Marker 2','Fast Sequence Transfer'):
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
        self.bIsStopped = False
        bWaveUpdate = False
        if self.bFastSeq:
            if not seq:
            # first waveform in the sequence
                if len(self.lStoredNames) != n_seq:
                # recreate lStoredNames 
                    self.lStoredNames = [[]] * n_seq
            self.decomposeWaveformAndSendFragments(seq)
        else:
            # go through all channels
            for n in range(self.nCh):
                # channels are numbered 1-4
                channel = n+1
                vData = self.getValueArray('Ch %d' % channel)
                vMark1 = self.getValueArray('Ch %d - Marker 1' % channel)
                vMark2 = self.getValueArray('Ch %d - Marker 2' % channel)
                bWaveUpdate = self.sendWaveformToTek(channel, vData, vMark1, vMark2, seq) or bWaveUpdate
                self.log('bWaveUpdate=%s' % bWaveUpdate)             
        # check if sequence mode
        if seq is not None:
            # if not final seq call, just return here
            self.bSeqUpdate = self.bSeqUpdate or bWaveUpdate
            if (seq+1) < n_seq:
                return
            # final call, check if sequence has changed
            self.log('bSeqUpdate=%s' % str(self.bSeqUpdate))
            if self.bSeqUpdate or n_seq != self.nOldSeq:
                # create sequence list, first clear to reset old values
                if self.bFastSeq:
                # notice that the list of the names has already been updated
                    self.assembleSequence()
                else:
                    self.writeAndLog('SEQ:LENG 0')
                    self.writeAndLog('SEQ:LENG %d' % n_seq)
                    for n1 in range(n_seq):               
                        for n2, bUpdate in enumerate(self.lInUse):
                            if bUpdate:
                                name = 'Labber_%d_%d' % (n2+1, n1+1)
                                self.writeAndLog('SEQ:ELEM%d:WAV%d "%s"' % \
                                                 (n1+1, n2+1, name))
                    # don't wait for trigger 
                        self.writeAndLog('SEQ:ELEM%d:TWA 0' % (n1+1))
                    # for last element, set jump to first
                    self.writeAndLog('SEQ:ELEM%d:GOTO:STAT 1' % n_seq)
                    self.writeAndLog('SEQ:ELEM%d:GOTO:IND 1' % n_seq)
                # save old sequence length
                self.nOldSeq = n_seq
            # turn on sequence mode
            self.writeAndLog(':AWGC:RMOD SEQ')
            # turn on channels in use 
            sOutput = ''
            for n, bUpdate in enumerate(self.lInUse):
                if bUpdate:
                    sOutput += (':OUTP%d:STAT 1;' % (n+1))
            if sOutput != '':
                self.writeAndLog(sOutput)
            return
        # turn on channels in use 
        sOutput = ''
        for n, bUpdate in enumerate(self.lInUse):
            if bUpdate:
                sOutput += (':OUTP%d:STAT 1;' % (n+1))
        if sOutput != '':
            self.writeAndLog(sOutput)
        # if not starting, make sure AWG is not running, then return
        if not bStart:
            iRunState = int(self.askAndLog(':AWGC:RST?'))
            nTry = 1000
            while nTry>0 and iRunState!=0 and not self.isStopped():
                # sleep for while to save resources, then try again
                self.wait(0.05)
                # try again
                iRunState = int(self.askAndLog(':AWGC:RST?'))
                nTry -= 1
            return
        # check if AWG has been stopped, if not return here
        if not self.bIsStopped:
            # no waveforms updated, just turn on output, no need to wait for run
            return
        # send command to turn on run mode to tek
        self.writeAndLog(':AWGC:RUN;')
        # wait for output to be turned on again
        iRunState = int(self.askAndLog(':AWGC:RST?'))
        nTry = 1000
        while nTry>0 and iRunState==0 and not self.isStopped():
            # sleep for while to save resources, then try again
            self.wait(0.05)
            # try again
            iRunState = int(self.askAndLog(':AWGC:RST?'))
            nTry -= 1
        # check if timeout occurred
        if nTry <= 0:
            # timeout
            raise InstrumentDriver.Error('Cannot turn on Run mode')
        # turn on channels again, to avoid issues when turning on/off run mode
        if sOutput != '':
            self.writeAndLog(sOutput)


    def scaleWaveformToU16(self, vData, dVpp, ch):
        """Scales the waveform and returns data in a string of U16"""
        # make sure waveform data is within the voltage range 
        if np.sum(vData > dVpp/2) or np.sum(vData < -dVpp/2):
            raise InstrumentDriver.Error(
                ('Waveform for channel %d contains values that are ' % ch) + 
                'outside the channel voltage range.')
        # clip waveform and store in-place
        np.clip(vData, -dVpp/2., dVpp/2., vData)
        vU16 = np.array(16382 * (vData + dVpp/2.)/dVpp, dtype=np.uint16)
        return vU16


    def createWaveformOnTek(self, channel, length, seq=None, bOnlyClear=False):
        """Remove old and create new waveform on the Tek. The waveform is named
        by the channel nunber"""
        if seq is None:
            name = 'Labber_%d' % channel
        else:
            name = 'Labber_%d_%d' % (channel, seq+1)
        # first, off output
        self.writeAndLog(':OUTP%d:STAT 0;' % channel)
        if bOnlyClear:
            # just clear this channel (only when Run Mode is not Sequence? From manual)
            self.writeAndLog(':SOUR%d:WAV ""' % channel)
        else:
            # remove old waveform, ignoring errors, then create new
            self.writeAndLog(':WLIS:WAV:DEL "%s"; *CLS' % name, bCheckError=False)
            self.writeAndLog(':WLIS:WAV:NEW "%s",%d,INT;' % (name, length))

                
    def sendWaveformToTek(self, channel, vData, vMark1, vMark2, seq=None):
        """Send waveform to Tek"""
        # check if sequence
        if seq is None:
            iSeq = 0
        else:
            iSeq = seq
        # channels are named 1-4
        n = channel-1
        if len(vData)==0:
            if len(vMark1)==0 and len(vMark2)==0:
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
        if (len(vMark1)>0 and len(vData)!=len(vMark1)) or \
           (len(vMark2)>0 and len(vData)!=len(vMark2)) or \
           (self.nPrevData>0 and len(vData)!=self.nPrevData):        
            raise InstrumentDriver.Error(\
                  'All channels need to have the same number of elements')
        self.nPrevData = len(vData)
        # channel in use, mark
        self.lInUse[n] = True
        # get range and scale to U16
        Vpp = self.getValue('Ch%d - Range' % channel)
        vU16 = self.scaleWaveformToU16(vData, Vpp, channel)
        # check for marker traces
        for m, marker in enumerate([vMark1, vMark2]):
            if len(marker)==len(vU16):
                # get marker trace
                vMU16 = np.array(marker != 0, dtype=np.uint16)
                # add marker trace to data trace, with bit shift
                vU16 += 2**(14+m) * vMU16
        start, length = 0, len(vU16)
        self.nOldLength = length
        # compare to previous trace
        self.log('---------channel=%d, length=%d, len=%d----------' % (channel, length, len(self.lOldU16[iSeq][n])))
        if length != len(self.lOldU16[iSeq][n]):
            # stop AWG if still running
            if not self.bIsStopped:
                self.writeAndLog(':AWGC:STOP;')
                self.bIsStopped = True
            # len has changed, del old waveform and create new
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
        if not self.bIsStopped:
            self.writeAndLog(':AWGC:STOP;')
            self.bIsStopped = True
        # create binary data as bytes with header
        sLen = b'%d' % (2*length)
        sHead = b'#%d%s' % (len(sLen), sLen)
        # send to tek, start by turning off output
        if seq is None:
            # non-sequence mode, get name
            name = b'Labber_%d' % channel
            sSend = b':OUTP%d:STAT 0;' % channel
            sSend += b':WLIS:WAV:DATA "%s",%d,%d,%s' % (name, start, length,
                     sHead + vU16[start:start+length].tobytes())
            self.write_raw(sSend)
            # (re-)set waveform to channel
            self.writeAndLog(':SOUR%d:WAV "%s"' % (channel, name.decode()))
            # self.log('-------waveform for channel%d has been sent--------' % channel)
        else:
            # sequence mode, get name
            name = b'Labber_%d_%d' % (channel, seq+1)
            sSend = b':OUTP%d:STAT 0;' % channel
            sSend += b':WLIS:WAV:DATA "%s",%d,%d,%s' % (name, start, length,
                     sHead + vU16[start:start+length].tobytes())
            # sSend += b'SEQ:ELEM%d:WAV%d' % (1, channel)
            self.write_raw(sSend)
        # store new waveform for next call
        self.lOldU16[iSeq][n] = vU16
        return True
    
    
    def checkWaveformLengthBeforeSending(self, channel, vData, vMark1, vMark2, seq):
        """Check waveform before sending to Tek. Similar to the method self.sendWaveformToTek(), but not the same. This is for fast sequence transfer"""
        # check if sequence
        if seq is None:
            iSeq = 0
        else:
            iSeq = seq
        # channels are named 1-4
        n = channel-1
        if len(vData)==0:
            if len(vMark1)==0 and len(vMark2)==0:
                # if channel in use, turn off, clear, go to next channel
                if self.lInUse[n]:
                    self.createWaveformOnTek(channel, 0, seq, bOnlyClear=True)
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
        if (len(vMark1)>0 and len(vData)!=len(vMark1)) or \
           (len(vMark2)>0 and len(vData)!=len(vMark2)) or \
           (self.nPrevData>0 and len(vData)!=self.nPrevData):        
            self.log('len(vData)=%d, len(vMark1)=%d, len(vMark2)=%d, self.nPrevData=%d' % (len(vData), len(vMark1), len(vMark2), self.nPrevData))
            raise InstrumentDriver.Error(\
                  'All channels need to have the same number of elements')
        self.nPrevData = len(vData)
        self.nOldLength = len(vData)
        # channel in use, mark
        self.lInUse[n] = True
        self.log('mDataMark = %s' % str(mDataMark))
        return (True, mDataMark)
           
    
    def scaleWaveformAndMarkerToU16(self, vData, vMark1, vMark2, channel):

        # get range and scale to U16
        Vpp = self.getValue('Ch%d - Range' % channel)
        vU16 = self.scaleWaveformToU16(vData, Vpp, channel)
        # check for marker traces
        for m, marker in enumerate([vMark1, vMark2]):
            if len(marker)==len(vU16):
                # get marker trace
                vMU16 = np.array(marker != 0, dtype=np.uint16)
                # add marker trace to data trace, with bit shift
                vU16 += 2**(14+m) * vMU16
        return vU16
    
    
    def decomposeWaveformAndSendFragments(self, seq):
        """"decompose the waveform into fragments (to assemble a sequence)"""
        self.log('decomposing the waveform for waveform %d' % seq)
        mDataMarkTot = 0 #initial value
        lValidChannel = []
        for n in range(self.nCh):
            # channels are numbered 1-4
            channel = n + 1
            vData = self.getValueArray('Ch %d' % channel)
            vMark1 = self.getValueArray('Ch %d - Marker 1' % channel)
            vMark2 = self.getValueArray('Ch %d - Marker 2' % channel)
            
            # similar process in self.sendWaveformToTek
            bValidChannel, mDataMark = self.checkWaveformLengthBeforeSending(channel, vData, vMark1, vMark2, seq)

            if bValidChannel:
                # this channel is used
                lValidChannel += [channel]
                if isinstance(mDataMarkTot, int):
                    # if this is the first valid channel
                    mDataMarkTot = mDataMark
                else:
                    # self.log('mDataMark = %s' % str(mDataMark))
                    # self.log('mDataMarkTot = %s' % str(mDataMarkTot))
                    mDataMarkTot = np.concatenate((mDataMarkTot, mDataMark))
             
                
        # find the zero points and non-zero points in a sequence. Because the 'repeat' function on AWG can not repeat different channels separately, all the channels should be considered together when trying to use 'repeat'.
        mNonzero = np.array(mDataMarkTot != 0, dtype = int)
        mDiff = np.diff(mNonzero)
        vDiff = np.array(np.any(mDiff != 0, axis = 0), dtype = int)
        vNode = np.nonzero(vDiff)[0]
        NumFrag = len(vNode) + 1 # number of fragments.
        vNode = np.concatenate(([0], vNode + 1, [len(vDiff) + 1]))
        lNames = [[]] * len(lValidChannel)
        for nC in range(len(lValidChannel)):
            channel = lValidChannel[nC]
            Vpp = self.getValue('Ch%d - Range' % channel)
            for nF in range(NumFrag):
                start = vNode[nF]
                end = vNode[nF + 1]
                # pick up from the 'start' element to 'end - 1' element
                ThisPulse = mDataMarkTot[3 * nC:3 * nC + 3, start:end]
                self.log('ThisPulse.shape = %s' % str(ThisPulse.shape))
                if np.all(np.diff(mDataMarkTot[:, start:end]) == 0):
                # constant pulse for all channels
                    DataV = ThisPulse[0, 0]
                    Mark1V = ThisPulse[1, 0]
                    Mark2V = ThisPulse[2, 0]
                    bWaveformsExist = True
                    
                    if seq:
                    #if not belong to the first waveform in a sequence, check whether the waveform already exists in the waveform list
                        if [DataV, Mark1V, Mark2V] == [0, 0, 0]:
                        # all zeros, check whether there are four waveforms in the waveform list ( length = 1, 10, 100, 1000 )
                            for i in range(4):
                                testname = 'zero_' + str(10 ** i)
                                if not testname in self.lWaveform:
                                    bWaveformsExist = False      
                        else:
                        # other constant pulses
                            for i in range(4):
                                testname = 'const_Vpp_' + str(Vpp) + str(DataV) + '_' + str(Mark1V) + '_' + str(Mark2V) + '_'  + str(10 ** i)
                                if not testname in self.lWaveform:
                                    bWaveformsExist = False                 
                    else:
                    # the first waveform, send fragments without checking
                        bWaveformsExist = False
                        
                    if [DataV, Mark1V, Mark2V] == [0, 0, 0]:
                    # all zeros
                        if not bWaveformsExist:
                        # send waveforms
                            for i in range(4):
                                mDataMark = np.zeros((3, 10 ** i))
                                self.createAndSendFragment(mDataMark, channel)
                        # store names        
                        ThisLength = len(ThisPulse[0, :])
                        sThisLength = str(ThisLength)
                        for i in range(min(4, len(sThisLength))):
                            if int(sThisLength[-( i + 1 )]): #nonzero
                                if i != 3:
                                    sLoopNum = sThisLength[-( i + 1 )]
                                else:
                                    sLoopNum = sThisLength[0:-i]
                                lNames[nC] = lNames[nC] + ['R' + sLoopNum + 'zero_' + str(10 ** ( i ))]    
                    else:
                    # other constant pulses
                        if not bWaveformsExist:
                        # send waveforms
                            for i in range(4):
                                mDataMark = np.zeros((3, 10 ** i))
                                mDataMark[0, :] = DataV
                                mDataMark[1, :] = Mark1V
                                mDataMark[2, :] = Mark2V
                                self.createAndSendFragment(mDataMark, channel)
                        # store names
                        ThisLength = len(ThisPulse[0, :])
                        sThisLength = str(ThisLength)
                        for i in range(min(4, len(sThisLength))):
                            if int(sThisLength[-( i + 1 )]): #nonzero
                                if i != 3:
                                    sLoopNum = sThisLength[-( i + 1 )]
                                else:
                                    sLoopNum = sThisLength[0:-i]
                                lNames[nC] = lNames[nC] + ['R' + sLoopNum + 'const_Vpp_' + str(Vpp) + str(DataV) + '_' + str(Mark1V) + '_' + str(Mark2V) + '_' + str(10 ** i)]  
                else:
                   # other complex shape pulses exist in at least one channel
                    lNames[nC] = lNames[nC] + [self.createAndSendFragment(ThisPulse, channel)]

        # compare to previous trace
        if lValidChannel != self.lValidChannel or lNames != self.lStoredNames[seq]:
            self.bSeqUpdate = True
            self.lValidChannel = lValidChannel
            lNames = np.array(lNames).transpose().tolist() # transpose the list, so that the sublists are divided by different fragments (therefore different sublists are different elements in the whole sequence)
            self.lStoredNames[seq] = lNames
            
        if not seq:
        # first waveform has been sent. get waveform list on the AWG
            self.getWaveformList()

    
    def createAndSendFragment(self, vFragment, channel):
        """ determine the name of the pulse and send it to the AWG """
        mDataMark = vFragment
        vData = mDataMark[0, :]
        vMark1 = mDataMark[1, :]
        vMark2 = mDataMark[2, :]
        Vpp = self.getValue('Ch%d - Range' % channel)
        if np.all(np.diff(mDataMark) == 0):
        # constant pulse
            if [vData[0], vMark1[0], vMark2[0]] == [0, 0, 0]:
            # all zero
                name = 'zero_' + str(len(vData))
            else:
                name = 'const_Vpp_' + str(Vpp) + str(vData[0]) + '_' + str(vMark1[0]) + '_' + str(vMark2[0]) + '_' + str(len(vData))
        else:
            time = self.askAndLog('SYST:TIME?', bCheckError=False)
            name = 'pulse' + time
            
        
        if not self.bIsStopped:
            self.writeAndLog(':AWGC:STOP;')
            self.bIsStopped = True
        

        # create binary data as bytes with header    
        vU16 = self.scaleWaveformAndMarkerToU16(vData, vMark1, vMark2, channel)
        start, length = 0, len(vU16)
        sLen = b'%d' % (2*length)
        sHead = b'#%d%s' % (len(sLen), sLen)
        # send to tek, start by turning off output
        self.writeAndLog(':OUTP%d:STAT 0;' % channel)    
        # remove old waveform, ignoring errors, then create new
        self.writeAndLog(':WLIS:WAV:DEL "%s"; *CLS' % name, bCheckError=False)
        self.writeAndLog(':WLIS:WAV:NEW "%s",%d,INT;' % (name, length))
        sname = name.encode()
        sSend = b':WLIS:WAV:DATA "%s",%d,%d,%s' % (sname, start, length,
                 sHead + vU16[start:start+length].tobytes())
        self.write_raw(sSend)
        return name
    
    
    def getWaveformList(self):
        """get the waveform list on the AWG"""
        size = int(self.askAndLog('WLIS:SIZE?', bCheckError=False))
        self.lWaveform = list(range(size))
        for i in range(size):
            self.lWaveform[i] = self.askAndLog('WLIS:NAME? ' + str(i), bCheckError=False)
            
            
    def assembleSequence(self):
        """assemble the fragment waveforms into a sequence"""
        self.writeAndLog('SEQ:LENG 0')
        lFlattenNames = [item for sublist in self.lStoredNames for item in sublist]
        #flatten the list, merge the lists for different seq together
        #which means:
        #    for sublist in l:
        #        for item in sublist:
        #            lFlattenNames.append(item)
        
        self.writeAndLog('SEQ:LENG %d' % len(lFlattenNames))
        
        for n1 in range(len(lFlattenNames)):               
            for n2, channel in enumerate(self.lValidChannel):
                name = lFlattenNames[n1][n2]
                loopnum = 1
                if name.startswith('R'):
                    ind = 1
                    while name[ind].isnumeric():
                        ind += 1
                    loopnum = int(name[1: ind])
                    name = name[ind:]
                self.writeAndLog('SEQ:ELEM%d:WAV%d "%s"' % \
                                 (n1 + 1, channel, name))
                if loopnum > 1 and n2 == 0:
                    self.writeAndLog('SEQ:ELEM%d:LOOP:COUNT %d' % (n1 + 1, loopnum))
        # don't wait for trigger 
            self.writeAndLog('SEQ:ELEM%d:TWA 0' % (n1+1))
        # for last element, set jump to first
        self.writeAndLog('SEQ:ELEM%d:GOTO:STAT 1' % len(lFlattenNames))
        self.writeAndLog('SEQ:ELEM%d:GOTO:IND 1' % len(lFlattenNames))
        self.log('Sequence has been established')
        # raise InstrumentDriver.Error()
    
    def debugPrint(self, arg):
        self.log('arg = %s' % str(arg))
    
if __name__ == '__main__':
    pass
