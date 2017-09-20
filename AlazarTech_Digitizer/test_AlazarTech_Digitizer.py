"""
This test assumes that the ATS board is continiously triggered and
all parameters are properly set via the Labber GUI. 
"""
import time

import Labber

# Connect to the Labber server.
client = Labber.connectToServer('localhost')

# Get list of instruments.
instruments = client.getListOfInstrumentsString()
print('Instrument List:')
for instr in instruments:
    print(instr)

# Connect to the AlazarTech Digitizer.
ats = client.connectToInstrument(
    'AlazarTech Digitizer', {
        'name': 'ATS', 'interface': 'Other', 'address': 1})
#Start the AlazarTech Digitizer.
ats.startInstrument()

ats.setValue('Number of records', 1000)

t0 = time.time()
ats.setValue('NPT AsyncDMA Enabled', False)
try:
    for ch in ('Channel A - Averaged Data', 'Channel B - Averaged Data'):
        print('%s shape: %s' % (ch, ats.getValue(ch)['y'].shape))
except Exception as e:
    print(e)
print('Single Port execution time: %f sec.\n' % (time.time() - t0))

t0 = time.time()
ats.setValue('NPT AsyncDMA Enabled', True)
try:
    for ch in ('Channel A - Averaged Data', 'Channel B - Averaged Data'):
        print('%s shape: %s' % (ch, ats.getValue(ch)['y'].shape))
except Exception as e:
    print(e)
print('NPT AsyncDMA execution time: %f sec.' % (time.time() - t0))


# Close the connection.
client.close()
