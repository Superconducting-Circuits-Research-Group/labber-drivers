# Waveform
- Sequence: Rabi, CP/CPMG, Pulse train, Generic sequence
- Sample rate: xxx samples per second
- Number of points: number of samples in one output?
- First pulse delays[s]: the time interval between the center of the pulse and the start of the output
- Trim waveform to sequence: send the waveform from the first non-zero point to the last non-zero point
- Number of outputs: maximum is Four. For each output, there are three channels, Trace-I, Trace-Q, Trace-Gate, that can be sent to the AWG.

# Sequence
- &#35 of pulses: number of pulses. The pulse can be a square pulse, a ramp pulse or a Gaussian pulse. The pulses can be in the same output or in different outputs.

# Pulse settings
- Pulse type: Square, Ramp, Gaussian
- Truncation range: only applies to Gaussian pulses. Truncate the pulse. The duration of the Gaussian pulse is TruncRange * Width + Plateau. See code: dTotTime = truncRange*dWidth + dPlateau
- Edge-to-edge pulses: after checked, increase 'Edge position' to make pulses more isolated from the other ones. The length of actual pulses is EdgePosition * Width + Plateau
- Use SSB mixing: certain transformation on pulses
- Use DRAG: certain scaling

# Pulse #1
- Amplitude: