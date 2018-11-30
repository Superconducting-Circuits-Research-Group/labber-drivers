#!/usr/bin/env python
import numpy as np

import InstrumentDriver


class Driver(InstrumentDriver.InstrumentWorker):
    """ This class implements a demodulation driver"""
    def performOpen(self, options={}):
        """Perform the operation of opening the instrument connection"""
        pass

    def performClose(self, bError=False, options={}):
        """Perform the close instrument connection operation"""
        pass

    def performSetValue(self, quant, value, sweepRate=0.0, options={}):
        """Perform the Set Value instrument operation. This function should
        return the actual value set by the instrument"""
        # do nothing here
        return value
        
    def pairwiseDifference(self, traces):
        return np.diff(traces)[:,0]
        
    def pairwiseRatio(self, traces):
        return traces[:,0] / traces[:,1]

    def performGetValue(self, quant, options={}):
        """Perform the Get Value instrument operation"""
        if quant.name in ('First value',
                          'Last value',
                          'Maximum',
                          'Minimum',
                          'Average',
                          'Median',
                          'Difference',
                          'Ratio',
                          'Pairwise difference',
                          'Pairwise ratio'):
            data = self.getValue('Input vector')
            traces = data['y']
            if quant.name == 'First value':
                value = traces[0]
            elif quant.name == 'Last value':
                value = traces[-1]
            elif quant.name == 'Maximum':
                value = np.max(traces)
            elif quant.name == 'Minimum':
                value = np.min(traces)
            elif quant.name == 'Average':
                value = np.mean(traces)
            elif quant.name == 'Median':
                value = np.median(traces)
            else:
                if not quant.name.startswith('Pairwise') and \
                        traces.size != 2:
                    traces = np.array([traces[0], traces[-1]])
                traces.shape = (-1, 2)
                if quant.name.lower().endswith('difference'):
                    value = self.pairwiseDifference(traces)
                else:
                    value = self.pairwiseRatio(traces)
                if quant.name.startswith('Pairwise'):
                    data['y'] = value
                    value = data
                else:
                    value = value[0]
        else:
            # just return the quantity value
            value = quant.getValue()
        return value

if __name__ == '__main__':
    pass
