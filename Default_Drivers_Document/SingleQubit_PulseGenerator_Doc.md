# Single Qubit Experiment Manual
## T_1 measurement
1. Set sequence to be 'Rabi'. Set 'Sample rate', 'Number of points', 'First pulse delay' properly. 
2. If you want to fix the $\pi$ pulse and move the readout pulse in each circle, turn off 'Trim waveform to sequence'. If not, turn on 'Trim waveform to sequence' and then turn on 'Buffer start to restore size'.
3. Set 'Number of outputs' to be 'One'
4. Set the parameters in 'Pulse settings' properly. For detailed explanation, see 'Common parameters' below. Typical settings: 'Pulse type' = 'Square', 'Truncation range' = 2, turn off 'Edge-to-edge pulses', 'Use SSB mixing', 'Use DRAG'.
5. Pulse #1 will be the $\pi$ pulse. Set the 'Amplitude', 'Width', 'Phase' according to the shape of the $\pi$ pulse. Set 'Plateau', 'Spacing' to be 0. Don't need to set the parameters for other pulses (In T_1 measurement, typically only pulse #1 will be used).
6. Turn off 'Generate tomography pulse'. Turn on 'Generate readout' and set 'Readout amplitude', 'Readout duration' properly. Set 'Readout delay' as the step parameter in 'Step sequence' list in 'Measurement Editor'

# Common parameters
## Waveform
### Waveform
- Sequence: Rabi, CP/CPMG, Pulse train, Generic sequence. Selecting Rabi simply means setting \# of pulses to 1. 
- Sample rate: xxx samples per second
- Number of points: number of samples
- First pulse delays[s]: the time interval between the center of the pulse and the start of the output
- Trim waveform to sequence: send the waveform from the first non-zero point to the last non-zero point
- Number of outputs: maximum is Four. For each output, there are three channels, Trace-I, Trace-Q, Trace-Gate, that can be sent to the AWG.

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
- Generate tomography pulse: generate tomography pulses if checked. The pulses will be $\pi$/2 rotations 
- State index: state index is cycled. If StateIndex % 3 = 0, the pulse is empty (z measurement). If StateIndex % 3 = 1, the pulse is an X(or Y) rotation (y(or x) measurement). If StateIndex % 3 = 2, the pulse is an Y(or X) rotation (x(or y) measurement) 
- Tomography delay: the time interval between the tomography pulse and other pulses
- Definition, pi/2 pulse: use a pulse to define the tomography $\pi$/2 pulse. The tomography pulse will have the same shape and output channel as the selected pulse do.

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