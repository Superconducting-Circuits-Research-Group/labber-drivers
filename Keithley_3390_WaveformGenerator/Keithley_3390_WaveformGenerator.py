import numpy as np

import InstrumentDriver
from VISA_Driver import VISA_Driver


class Driver(VISA_Driver):
    """This class implements the Keithley 3390 AWG."""

    def performOpen(self, options={}):
        """Perform the operation of opening the instrument connection."""
        # add compatibility with pre-python 3 version of Labber
        if not hasattr(self, 'write_raw'):
            self.write_raw = self.write
        # start by calling the generic VISA open to make sure we have
        # a connection
        VISA_Driver.performOpen(self, options)
        # clear value of waveform
        self.setValue('Arbitrary waveform', [])
        self._sent_waveform = np.array([])

    def performSetValue(self, quant, value, sweepRate=0.0, options={}):
        """Perform the Set Value instrument operation. This function
        should return the actual value set by the instrument."""
        # keep track of if waveform is updated, to avoid sending it many
        # times
        if self.isFirstCall(options):
            self.bWaveUpdated = False
        if quant.name == 'Arbitrary waveform':
            # set value, then mark that waveform needs an update
            quant.setValue(value)
            self.bWaveUpdated = True
        else:
            # there seem to be a descipancy on what voltage mean
            if quant.name == 'Voltage':
                value /= 2.
            # for all other cases, call VISA driver
            value = VISA_Driver.performSetValue(self, quant, value,
                    sweepRate, options=options)
            if quant.name == 'Voltage':
                value *= 2.
        # if final call and wave is updated, send it to AWG
        if self.isFinalCall(options) and self.bWaveUpdated:
            self.sendWaveform()
        return value

    def performGetValue(self, quant, options={}):
        """Perform the Get Value instrument operation."""
        value = VISA_Driver.performGetValue(self, quant, options)
        if quant.name == 'Voltage':
            value *= 2.
        return value

    def sendWaveform(self):
        """Rescale and send waveform data to the AWG."""
        # get data
        vData = self.getValueArray('Arbitrary waveform')
        # scale to U16
        vData = np.clip(vData, -1, 1)
        vI16 = np.array(4096 * vData, dtype=np.int16)
        if np.array_equal(self._sent_waveform, vI16):
            return
        self._sent_waveform = vI16
        length = len(vI16)
        # create data as bytes with header
        sLen = b'%d' % (2 * length)
        sHead = b':DATA:DAC VOLATILE, #%d%s' % (len(sLen), sLen)
        # write header + data
        trigger = self.readValueFromOther('Trigger source')
        self.sendValueToOther('Trigger source', 'Manual')
        self.write_raw(sHead + vI16.tobytes())
        # select volatile waveform
        self.write(':FUNC:USER VOLATILE')
        self.sendValueToOther('Trigger source', trigger)


if __name__ == '__main__':
    pass
