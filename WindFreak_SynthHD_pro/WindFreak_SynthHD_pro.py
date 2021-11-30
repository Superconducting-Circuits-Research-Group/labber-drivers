from VISA_Driver import VISA_Driver

class Driver(VISA_Driver):

    def performSetValue(self, quant, value, sweepRate=0.0, options={}):
        """Perform the Set Value instrument operation. This function
        should return the actual value set by the instrument."""
        # check quantity name
        if quant.name == 'FrequencyA':
            self.writeAndLog('C0')
            self.writeAndLog('f%f' % (1.e-6 * value))
        elif quant.name == 'PowerA':
            self.writeAndLog('C0')
            self.writeAndLog('W%f' % value)
        elif quant.name == 'OutputA':
            self.writeAndLog('C0')
            self.writeAndLog('r%d' % value)
        elif quant.name == 'FrequencyB':
            self.writeAndLog('C1')
            self.writeAndLog('f%f' % (1.e-6 * value))
        elif quant.name == 'PowerB':
            self.writeAndLog('C1')
            self.writeAndLog('W%f' % value)
        elif quant.name == 'OutputB':
            self.writeAndLog('C1')
            self.writeAndLog('r%d' % value)
        else:
            # otherwise, call standard VISA case
            VISA_Driver.performSetValue(self, quant, value, sweepRate,
                    options)
        return value

    def performGetValue(self, quant, options={}):
        """Perform the Get Value instrument operation."""
        if quant.name == 'FrequencyA':
            self.writeAndLog('C0')
            value = 1.e6 * float(self.ask('f?'))
        elif quant.name == 'PowerA':
            self.writeAndLog('C0')
            value = self.ask('W?')
        elif quant.name == 'OutputA':
            self.writeAndLog('C0')
            value = self.ask('r?')
            value = int(value)
        elif quant.name == 'FrequencyB':
            self.writeAndLog('C1')
            value = 1.e6 * float(self.ask('f?'))
        elif quant.name == 'PowerB':
            self.writeAndLog('C1')
            value = self.ask('W?')
        elif quant.name == 'OutputB':
            self.writeAndLog('C1')
            value = self.ask('r?')
            value = int(value)
        else:
            # otherwise, call standard VISA case
            value = VISA_Driver.performGetValue(self, quant, options)
        return value


if __name__ == '__main__':
    pass