# Single Qubit Experiment Manual
## T1 measurement


# Common parameters
## Sequence
### Sequence
- Sequence: Rabi, CP/CPMG, Pulse train, C-phase Pulses, C-phase Echo, 1-QB Randomized Benchmarking, 2-QB Randomized Benchmarking, Custom. 
- Number of qubits: maximum is 9. For each qubit, there are 8 channels, 'Trace - I', 'Trace - Q', 'Trace - G', 'Trace - Z', 'Single-shot, QB - Real', 'Single-shot, QB - Imag', 'Single-shot, QB - Magnitude', 'Single-shot, QB - Phase' that can be sent to the AWG.
- Pulse period, 1-QB: Period for single-qubit gates.
- Pulse period, 2-QB: Period for two-qubit gates. Single-qubit and two-qubit gates are used in 'Pulse train', '1-QB Randomized Benchmarking' and '2-QB Randomized Benchmarking' modes. Note: the qubit gate here is different from the gate in 'Microwave gate switch'. The former exists in 'Trace - I', 'Trace - Q', or 'Trace - Z'. The latter exists in 'Trace - G'.
- Local XY control: If False, all control pulses are added to a single output waveform, i.e. the first qubit output. The outputs for other qubit will be 0.

### Rabi
- Uniform pulse amplitude for all qubits: If True, the amplitude of pulse #1 will be used for all qubits

## Waveform
### Waveform
- Sample rate: xxx samples per second
- Number of points: number of samples
- First pulse delays[s]: the time interval between the center of the pulse and the start of the output
- Trim waveform to sequence: send the waveform from the first point to the last non-zero point
- Trim both start and end: send the waveform from the first non-zero point to the last non-zero point
- Align pulses to end of waveform: If Ture, the pulse will be at the end of the waveform regardless of 'First pulse delays', unless it is very large or negative.

## 1-QB gates
### Pulse settings
- Pulse type: Square, Ramp, Gaussian
- Truncation range: Visible if 'Pulse type' is 'Gaussian'. Truncate the pulse. The duration of the Gaussian pulse is TruncRange * Width + Plateau.
- Use DRAG: certain scaling
- Uniform pulse shape: If True, set the width and plateau for all 'Pulse #%d'

### Pulse #1
- Amplitude: amplitude of the pulse #1
- Width: the width of the pulse #1
- Plateau: determines how long the maximum of the pulse will last
- Frequency: the frequency of pulse. The parameters mentioned earlier are describing the envelope of the pulse. 'Trace - I' is a cosine wave and 'Trace - Q' is a sine wave. The end of the pulse corresponds to zero phase point for each case.

## 2-QB gates
### 2-QB pulses
- Pulse type: Gaussion, Square, Ramp, CZ. For CZ pulse, the notation and calculations are based on the Paper "Fast adiabatic qubit gates using only sigma_z control" PRA 90, 022307 (2014)
- Uniform 2QB pulses: If True, set the width and plateau for all '2-QB pulse #%d%d'. There seems a bug when it is False while 'Pulse type' is 'CZ'. 

## Tomography
### State tomography
- Generate tomography pulse: This part still remains to be done. Don't use.

## Predistortion
### Predistortion
This driver appears to take a premade transfer function for each mixer to predistort I/Q waveforms for qubit XY control. The trick is then to load in the transfer function and perform the actual predistortion. We therefore need to find the correct transfer function and save this transfer function to file. These transfer functions need to be saved as Transfer function #1, etc.
No idea on how to write the transfer function...

## Cross-talk
### Cross-talk
Remians to be done. Don't use.
- Compensate cross-talk: if True, compensate for Z-control crosstalk
- Cross-talk(CT) matrix: remains to be done

## Readout
### Readout trig
- Generate readout trig: if True, generate readout trigger in 'Trace - Readout Trig'
- Readout trig amplitude: amplitude
- Readout trig duration: duration

### Readout
- Generate readout waveform: if True, generate readout pulse. The pulse will be in the 'Trace - Readout I' and 'Trace - Readout Q' channel
- Number of readout tones: readout is a combination of different tones listed below. For each tone you can set the frequency and the amplitude. 'Trace - I' is a summation of cosine waves and 'Trace - Q' is a summation of sine waves. The end of the pulse corresponds to zero phase point for each case.
- Uniform readout amplitude: if True, assign the same amplitude for all tones
- Readout duration: the width of the readout pulse
- Readout delay: the time interval between the readout pulse and the last pulse.
- Match main sequence waveform size: I think it is useless. (For the version downloaded from Github, there was a bug that if this option is False, the readout waveform will start from t=0.)
- Readout frequency: the frequency of that tone
- Readout amplitude: the amplitude of that tone

## Output
### Output
- Swap IQ: if True, swap IQ channels for all qubit pulses.

### Microwave gate switch
- Generate gate: if True, generate gate switch in 'Trace - G'
- Uniform gate: if True, gate pulses will cover almost the whole waveform for each output
- Gate delay: (for non-uniform gate) the time interval between the gate pulse and the output pulse
- Gate overlap: (for non-uniform gate) to increase the width of the gate pulse. If overlap is 0, the gate will have the same width as the pulse does.
- Minimal gate time: (for non-uniform gate) This should be "Minimal gap time"? If the gap time between two gates is less than this value, the gap will be filled up.

## Demodulation
We don't use this part.