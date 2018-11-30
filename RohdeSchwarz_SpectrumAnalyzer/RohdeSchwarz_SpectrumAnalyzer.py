#!/usr/bin/env python

from VISA_Driver import VISA_Driver
import numpy as np


class Error(Exception):
    pass


class Driver(VISA_Driver):
    """This class implements the Rohde&Schwarz Network Analyzer driver."""

    def performOpen(self, options={}):
        """Perform the operation of opening the instrument connection."""
        # calling the generic VISA open to make sure we have a connection
        VISA_Driver.performOpen(self, options=options)

    def performSetValue(self, quant, value, sweepRate=0.0, options={}):
        """
        Perform the Set Value instrument operation. This function should
        return the actual value set by the instrument.
        """
        value = VISA_Driver.performSetValue(self, quant, value, sweepRate, options)
        return value

    def performGetValue(self, quant, options={}):
        """Perform the Get Value instrument operation."""
        # check type of quantity
        if quant.name in ['Signal', 'Peak power', 'Peak frequency']:
            if self.isFirstCall(options):
                # check if channel is on
                # if not in continous mode, trig from computer
                bWaitTrace = self.getValue('Wait for new trace')
                bAverage = self.getValue('Average')
                # wait for trace, either in averaging or normal mode
                if bWaitTrace:
                    if bAverage:
                        nAverage = self.getValue('# of averages')
                        self.writeAndLog(':ABOR;:INIT:CONT OFF;:SENS:AVER:COUN %d;:INIT:IMM;*OPC' % nAverage)
                    else:
                        self.writeAndLog(':ABOR;:INIT:CONT OFF;:SENS:AVER:COUN 1;:INIT:IMM;*OPC')
                    bDone = False
                    while (not bDone) and (not self.isStopped()):
                        stb = int(self.askAndLog('*ESR?'))
                        bDone = stb & 0b1
                        if not bDone:
                            self.wait(0.1)
                    # if stopped, don't get data
                    if self.isStopped():
                        self.writeAndLog('*CLS;:INIT:CONT ON;')
                        return []
                # get data as float32, convert to numpy array
                sData = self.ask(':FORM ASC;TRAC1? TRACE1')
                vData = np.array(sData.split(',')).astype(dtype=float)
                if bWaitTrace and not bAverage:
                    self.writeAndLog(':INIT:CONT ON;')
                startFreq = self.readValueFromOther('Start frequency')
                stopFreq = self.readValueFromOther('Stop frequency')
                self._vData = quant.getTraceDict(vData, x0=startFreq, x1=stopFreq)

            # create a trace dict
            if quant.name == 'Signal':
                value = self._vData
                self.log(value)
            elif quant.name == 'Peak power':
                value = np.max(self._vData['y'])
            elif quant.name == 'Peak frequency':
                idx = np.median(np.argmax(self._vData['y']))
                value = self._vData['t0'] + float(idx) * self._vData['dt']
        elif quant.name in ['Wait for new trace']:
            # do nothing, return local value
            value = quant.getValue()
        else:
            # for all other cases, call VISA driver
            value = VISA_Driver.performGetValue(self, quant, options)
        return value


if __name__ == '__main__':
    pass
