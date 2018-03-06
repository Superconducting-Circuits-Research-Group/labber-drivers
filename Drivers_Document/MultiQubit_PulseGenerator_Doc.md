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

Below is for SingleQubit-PulseGenerator



### Sequence
- \# of pulses: number of pulses used in the sequence. The pulse can be a square pulse, a ramp pulse or a Gaussian pulse. The pulses can be in the same output or in different outputs. (except for the Rabi sequence)

## Pulses
### Pulse settings
- Pulse type: Square, Ramp, Gaussian
- Truncation range: only applies to Gaussian pulses. Truncate the pulse. The duration of the Gaussian pulse is TruncRange * Width + Plateau. See code: dTotTime = truncRange * dWidth + dPlateau
- Edge-to-edge pulses: after checked, increase 'Edge position' to make pulses more isolated from the other ones. The length of actual pulses is EdgePosition * Width + Plateau
- Use SSB mixing: certain transformation on pulses
- Use DRAG: certain scaling

### Pulse #1
- Amplitude: amplitude of the pulse #1
- Width: the width of the pulse #1
- Plateau: determines how long the maximum of the pulse will last
- Spacing: time interval between this pulse and the next pulse. Changing spacing will move the positions of pulse #2, #3, ...
- Phase: roughly speaking, the amplitude for trace-I is Amplitude * Cos(Phase) and the amplitude for trace-Q is Amplitude * Sin(Phase)
- Output: assign the pulse to certain output. The output with more than one pulses is a summation of these pulses.

## Readout
### state tomography
- Generate tomography pulse: generate tomography pulses if checked. The pulses will be pi/2 rotations 
- State index: state index is cycled. If StateIndex % 3 = 0, the pulse is empty (z measurement). If StateIndex % 3 = 1, the pulse is an X(or Y) rotation (y(or x) measurement). If StateIndex % 3 = 2, the pulse is an Y(or X) rotation (x(or y) measurement) 
- Tomography delay: the time interval between the tomography pulse and other pulses
- Definition, pi/2 pulse: use a pulse to define the tomography pi/2 pulse. The tomography pulse will have the same shape and output channel as the selected pulse do.

### Readout
- Generate readout: generate readout pulse if checked. The pulse will be in the Trace-Readout channel
- Readout delay: the time interval between the readout pulse and the last pulse.
- Readout amplitude: the amplitude of the readout pulse
- Readout duration: the width of the readout pulse
- Sample-and-hold readout: if checked, add a tail to the readout pulse.

## Output
### Output
- Swap IQ: if checked, swap IQ channels for all outputs.

### Pre-pulses
- Add pre-pulses: if checked, add some identical pulses before all the other pulses
- Number of pre-pulses: number of pre-pulses
- Pre-pulse period: period of pre-pulses
- Pre-pulse definition: use a pulse to define the pre-pulses. The pre-pulses will have the same shape and output channel as the selected pulse do.

### Gate
- Generate gate: if checked, add gate pulses in Trace-Gate channel. There will be gate pulses for all pulses in each output channel with "Uniform gate" option unchecked
- Uniform gate: if checked, gate pulses will cover almost the whole waveform for each output
- Gate delay: (for non-uniform gate) the time interval between the gate pulse and the output pulse
- Gate overlap: (for non-uniform gate) to increase the width of the gate pulse. If overlap is 0, the gate will have the same width as the pulse does.
- Minimal gate time: (for non-uniform gate) This should be "Minimal gap time"? If the gap time between two gates is less than this value, the gap will be filled up.