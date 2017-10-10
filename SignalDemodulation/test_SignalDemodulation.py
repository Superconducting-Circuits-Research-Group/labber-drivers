import time
import numpy as np

length = 100000
R = 500
skipIndex = 100
vData = np.random.rand(R, skipIndex + length)
dFreq = 10e-9

t0 = time.clock()
vTime = np.arange(length, dtype=float)
vCos = np.cos(2 * np.pi * vTime * dFreq)
vSin = np.sin(2 * np.pi * vTime * dFreq)
# calc I/Q
dI = np.trapz(vCos * vData[:,skipIndex:skipIndex+length]) / float(length-1)
dQ = np.trapz(vSin * vData[:,skipIndex:skipIndex+length]) / float(length-1)
signal = dI + 1j*dQ
print(signal[:5])
print('Original: %f sec.\n' % (time.clock() - t0))

t0 = time.clock()
vTime = np.arange(length, dtype=float)
vPhase = 2 * np.pi * vTime * dFreq
vCos = np.cos(vPhase)
vSin = np.sin(vPhase)
# calc I/Q
dI = np.empty(R, dtype=float)
dQ = np.empty(R, dtype=float)
np.dot(vData[:,skipIndex:skipIndex+length], vCos, dI)
np.dot(vData[:,skipIndex:skipIndex+length], vSin, dQ)
signal = (dI + 1j * dQ) / float(length - 1)
print(signal[:5])
print('Method #1: %f sec.\n' % (time.clock() - t0))

t0 = time.clock()
vTime = np.arange(length, dtype=float)
vPhase = 2 * np.pi * vTime * dFreq
vExp = np.exp(1j * vPhase)
# calc I/Q
signal = np.empty(R, dtype=complex)
np.dot(vData[:,skipIndex:skipIndex+length], vExp, signal)
signal /= float(length - 1)
print(signal[:5])
print('Method #2: %f sec.\n' % (time.clock() - t0))