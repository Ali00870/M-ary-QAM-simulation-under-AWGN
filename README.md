# M-QAM Simulation with RRC Pulse Shaping

This MATLAB script simulates M-ary Quadrature Amplitude Modulation (M-QAM) transmission over an AWGN channel using Root Raised Cosine (RRC) pulse shaping. It plots the simulated Symbol Error Rate (SER) and compares it to the theoretical SER over a range of Eb/N0 values.

## Overview

The simulation includes:
- Random M-QAM symbol generation
- RRC pulse shaping and matched filtering
- AWGN channel noise
- Symbol detection and SER calculation
- Comparison with theoretical SER

This code is general and supports any square M-QAM constellation (e.g., 16-QAM, 64-QAM) by specifying the modulation order `M`. It can be used to study performance degradation due to noise and to understand the role of pulse shaping in communication systems.

## File

- `simulate_MQAM.m` — Main simulation function

## Usage

In MATLAB, call the function as:

```matlab
[EbN0_dB, SER, theo_SER] = simulate_MQAM(M, d, show_plot);
