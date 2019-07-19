from VISA_Driver import VISA_Driver


class Driver(VISA_Driver):
    """This class implements HP 8672A Signal Generator driver."""

    def performSetValue(self, quant, value, sweepRate=0.0, options={}):
        """Perform the Set Value instrument operation. This function
        should return the actual value set by the instrument."""
        # check quantity name
        if quant.name == 'Frequency':
            self.write('P%08.fZ0K0L3M0N6O3' % (value / 1000.))
        else:
            # otherwise, call standard VISA case
            VISA_Driver.performSetValue(self, quant, value, sweepRate,
                    options)
        return value


if __name__ == '__main__':
	pass