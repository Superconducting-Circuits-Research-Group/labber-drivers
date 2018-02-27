# Waveform
- Sequence: Rabi, CP/CPMG, Pulse train, Generic sequence
- Sample rate: xxx samples per second
- Number of points: number of samples
- First pulse delays[s]: the time interval between the center of the pulse and the start of the output
- Trim waveform to sequence: send the waveform from the first non-zero point to the last non-zero point
- Number of outputs: maximum is Four. For each output, there are three channels, Trace-I, Trace-Q, Trace-Gate, that can be sent to the AWG.

# Sequence
- \# of pulses: number of pulses used in the sequence. The pulse can be a square pulse, a ramp pulse or a Gaussian pulse. The pulses can be in the same output or in different outputs. (except for the Rabi sequence)

# Pulse settings
- Pulse type: Square, Ramp, Gaussian
- Truncation range: only applies to Gaussian pulses. Truncate the pulse. The duration of the Gaussian pulse is TruncRange * Width + Plateau. See code: dTotTime = truncRange * dWidth + dPlateau
- Edge-to-edge pulses: after checked, increase 'Edge position' to make pulses more isolated from the other ones. The length of actual pulses is EdgePosition * Width + Plateau
- Use SSB mixing: certain transformation on pulses
- Use DRAG: certain scaling

# Pulse #1
- Amplitude: amplitude of the pulse #1
- Width: the width of the pulse #1
- Plateau: determines how long the maximum of the pulse will last
- Spacing: distance between this pulse and the next pulse. Changing spacing will move the positions of pulse #2, #3, ...