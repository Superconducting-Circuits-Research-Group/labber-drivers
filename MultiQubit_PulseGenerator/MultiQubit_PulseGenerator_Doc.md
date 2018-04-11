# Multi Qubit Experiment Manual

## T1 measurement
1. Set sequence to be 'CP/CPMG'.
   Set 'Number of qubits' to be 1 (single qubit case), 'Local XY' to be True. 
   'Pulse period' won't be used.
   Set 'Number of pi pulses' to be -1 (This stands for 'T1 measurement').
   Set 'Sequence duration' to be 0. 
   'Add pi pulses to Q', 'Edge-to-edge pulses' won't be used.

2. Set 'Sample rate', 'Number of points', 'First pulse delay' properly.
If you want to fix the pi pulse and move the readout pulse for each
waveform in a sequence, turn off 'Align pulses to end of waveform'. If not, turn it on. In this case, 'First pulse delay' won't affect the sequence unless it is very small. Turn off 'Trim waveform to sequence'.
3. Set the parameters in '1-QB gates' properly. For detailed explanation, see 'Common parameters' below. Typical settings: 'Pulse type' is 'Square', 'Use DRAG' is off, 'Uniform pulse shape' is on.
4. Pulse #1 will be the pi pulse. Set the 'Width', 'Amplitude' according to the shape of the pi pulse. Set 'Plateau', 'Frequency' to be 0. Don't need to set the parameters for other pulses: in T1 measurement, only pulse #1 is used.
5. '2-QB gates', 'Tomography', 'Predistortion', 'Cross-talk' won't be used.
6. Turn on 'Generate readout trig'. Set 'Readout trig amplitude' to be 0.5V. Set 'Readout trig duration' to be 1e-6s. Turn on 'Generate readout waveform'. Set 'Number of readout tones' to be 1. Set 'Readout amplitude' = 0.5V, 'Readout duration' = 1e-5s. Set 'Readout delay' as the step parameter in 'Step sequence' list in 'Measurement Editor'.
7. 'Match main sequence waveform size' should be on. Set 'Readout frequency #1' = 0. Other frequencies are irrelevant. Turn off 'Distribute readout phases'. Set 'Readout I/Q ratio' = 1, 'Offset I' = 0, 'Offset Q' = 0, 'IQ skew' = 0. Turn off 'Predistort readout waveform'.
8. In 'Signal connections', set 'Trace - I1' to be the source for 'Ch1' of 'Tektronix AWG'. If 'Phase' for pulse #1 is non-zero, set 'Trace - Q' to be the source for 'Ch2'. Set 'Trace - Readout' to be the source for 'Ch3'. Set 'Trace - Readout trig' to be the source for 'Ch3 - Marker 1' as the trigger.
9. In 'Step sequence', set all the other parameters. Note: 
	- 'AlazarTech Signal Demodulator' 
		- 'Records per buffer' = the number of step points. 
		- 'Trigger level' = 450mV.
		- 'Sequence time step' = the step of the 'Readout delay'.
		- 'Acquisition mode' = 'Referenced Average Buffer Demodulation with AWG Hardware Loop' 
	- 'Tektronix AWG'
		- 'Run mode' = 'Sequence'.

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

### CPMG
- \Number of pi pulses: if -1, only one pi pulse (T1 experiment); 
if 0, only two pi/2 pulses (Ramsey experiment);
if 1, one pi pulse sandwiched by two pi/2 pulses (spin-echo experiment);
if more than one, number of pi pulses between the two pi/2 pulses.
The pi pulses are defined in section '1-QB gates'. In that section, different pulses represent pi pulses for different qubits. pi/2 pulses have half the amplitudes.
- Sequence duration: the time interval between two pi/2 pulses. If it is T1 experiment, increasing 'Sequence duration' moves readout pulse further from the pi pulse. If 'Sequence duration' is too small, the pi/2 pulses and the pi pulses will be added together.
- Add pi pulses to Q: the added pi pulse will be in 'Trace - Q' channel
- Edge-to-edge pulses: increase the time interval for pi pulses if it has non-zero plateau or truncation range (for Gaussian pulses).

## Waveform
### Waveform
- Sample rate: xxx samples per second
- Number of points: number of samples
- First pulse delays[s]: the time interval between the center of the pulse and the start of the output
- Trim waveform to sequence: send the waveform from the first point to the last non-zero point. But it doesn't trim the readout waveform.
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
- Match main sequence waveform size: if False, the readout waveform will be independent on the position of the other pulses. Otherwise the readout waveform will follow the readout trigger. If 'Trim waveform to sequence' is True but this option is False, the readout waveform isn't tirmmed.
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