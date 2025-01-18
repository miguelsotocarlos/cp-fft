# Competitive Programming FFT
This is an implementation of fast convolutions using dicrete FFT over prime fields (also referred to as Number-Theoretic Transform) optimized for competitive programming, specifically for ICPC-style competitions with a printed reference material. This code is part of the reference material used by the InChaVoLa team in the ICPC World Finals.

## Usage
The code provides three functions
  - `conv_small` computes convolutions modulo `998244353`.
  - `conv_big` computes convolutions modulo a 62-bit prime.
  - `conv_sunzi` computes convolutions modulo an arbitrary 32-bit modulus by combining two 62-bit convolutions with the Sunzi theorem.

Both `conv_small` and `conv_big` can be adapted to different primes, provided that those primes have a suitable primitive root.

The global parameter `LG` determines the largest possible convolution size to be `2^LG`. Due to time or memory requirements, you may need to lower it from its default value of `23`.

## Optimizations
The three most important optimizations in this project are
  - The use of the Montgomery form allows much more efficient modular operations, especially for large primes of around 64 bits.
  - The use of incompletely-reduced modular values in intermediate computations reduces the number of reductions necessary.
  - The use of both decimation-in-time (for the inverse FFT) and decimation-in-frequency (for the forward FFT) allows for an implementation without a bit-reversal permutation.

## Project structure
The main code resides in `fft.cpp`. This is the part intended to be included in your reference material. A commented version of the main implementation can be found in `fft_commented.cpp`.

There are randomized tests against a naive quadratic implementation in `tests.cpp`, which can be compiled with `g++ tests.cpp -o tests` and ran with `./tests` (they take a while, around 7min). There are also submission tests to the SPOJ and codeforces online judges, which can be found in `submissions/`.
