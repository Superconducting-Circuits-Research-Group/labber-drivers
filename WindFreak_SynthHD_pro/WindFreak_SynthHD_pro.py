from VISA_Driver import VISA_Driver


class Driver(VISA_Driver):

    def performSetValue(self, quant, value, sweepRate=0.0, options={}):
        """Perform the Set Value instrument operation. This function
        should return the actual value set by the instrument."""
        # check quantity name
        if quant.name == 'Frequency':
            self.writeAndLog('f%f' % (1.e-6 * value))
        elif quant.name == 'Power':
            self.writeAndLog('W%f' % value)
        elif quant.name == 'Output':
            self.writeAndLog('r%d' % value)
        else:
            # otherwise, call standard VISA case
            VISA_Driver.performSetValue(self, quant, value, sweepRate,
                    options)
        return value

    def performGetValue(self, quant, options={}):
        """Perform the Get Value instrument operation."""
        if quant.name == 'Frequency':
            value = 1.e6 * float(self.ask('f?'))
        elif quant.name == 'Power':
            value = self.ask('W?')
        else:
            # otherwise, call standard VISA case
            value = VISA_Driver.performGetValue(self, quant, options)
        return value


if __name__ == '__main__':
    pass