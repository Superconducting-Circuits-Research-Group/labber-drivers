from VISA_Driver import VISA_Driver


class Driver(VISA_Driver):
    """This class implements HP 8672A Signal Generator driver."""

    def performSetValue(self, quant, value, sweepRate=0.0, options={}):
        """Perform the Set Value instrument operation. This function
        should return the actual value set by the instrument."""
        # check quantity name
        if quant.name == 'Frequency':
            self.write('FR%011.0fHZ' % value)
        elif quant.name == 'Power':
            self.write('LE%+.1fDM' % value)
        else:
            # otherwise, call standard VISA case
            VISA_Driver.performSetValue(self, quant, value, sweepRate,
                    options)
        return value

    def performGetValue(self, quant, options={}):
        """Perform the Get Value instrument operation."""
        if quant.name == 'Frequency':
            value = self.ask('FR')
            value = float(value.strip().strip('FR').strip('HZ'))
        elif quant.name == 'Power':
            # value = self.as('LEOA;K')
            value = float(value.strip().strip('LE').strip('DM'))
        else:
            # otherwise, call standard VISA case
            value = VISA_Driver.performGetValue(self, quant, options)
        return value


if __name__ == '__main__':
	pass