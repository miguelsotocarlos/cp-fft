// Competitive Programming FFT
// Copyright (C) 2025 Carlos Miguel Soto <miguelsotocarlos@gmail.com>

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published by
// the Free Software Foundation, version 3.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.

// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include <bits/stdc++.h>
using namespace std;
#define forsn(i,s,n) for(int i = int(s);i<int(n);i++)
#define dforsn(i,s,n) for(int i = int(n)-1;i>=int(s);i--)
#define fore(i,s,n) forsn(i,s,n)
#define forn(i,n) forsn(i,0,n)
#define dforn(i,n) dforsn(i,0,n)
#define SZ(v) int(v.size())

// The maximum length of the resulting
// convolution vector is 2^LG
const int LG = 23;
using ll = long long;

// u must be uint32_t or uint64_t
// uu must be uint64_t or __uint128_t, and must have
// twice as many bits as u
// p must be a prime number such that 4*p fits in u
// root must be a number mod p such that the 2-adicity
// of its multiplicative order is maximal,
// for example a primitive root
template<class u, class uu, u p, u root>
struct FFT {
	u r[1+(2<<LG)];

	// m multiplies two numbers in montgomery form
	// that is, it returns ab/2^k mod p
	// where k = 32 or 64 depending on the size of u
	// the result is guaranteed to be in [0, 2p)
	// so it's not fully reduced
	constexpr u m(u a, u b) {
		uu k = uu(a)*b;
		// op represents a step in newton's method for calculating
		// the negative multiplicative inverse of p modulo 2^k
		// thus, (op(op(op(op(op(-p)))))) * p should be -1
		#define op(g) g*(g*p+2)
		k += u(k) * (op(op(op(op(op(-p)))))) * uu(p);
		#undef op
		return u(k>>(8*sizeof(u)));
	}

	constexpr u red(u k, u a) { return a-k*(a>=k); }

	FFT() {
		// Precalculates the necessary roots of unity
		// -p%p is the montgomery form for 1
		u k = r[2<<LG] = -p%p,
			// Here we use b directly instead of the montgomery
			// for of b. This is ok because the order of b has
			// sufficient 2-adicity if and only if the order of
			// b * 2^k has, provided k is even. And in our case
			// k can only be 32 or 64
			b=root,
			e = p>>LG;

		// First, the primitive pth-root of unity is raised
		// to the power of p/2^k to make it a primitive 2^k-th root
		for(; e; e/=2, b=m(b,b)) if(e%2) k=m(k, b);

		// the array r is divided into power of two blocks
		// [1]
		// [2, 3]
		// [4, 5, 6, 7]
		// ...
		// In each block, the ith element is w^-ij where j = 2^k/blocksize
		// (the negative sign in the exponent is not important)
		dforn(i, 2<<LG)
			// `k` stores the current power w^j
			r[i]=red(p, m(r[i+1], k)),
			// Each time we reach a power of two, k must be updated (squared)
			i&(i-1)?0:k=m(k,k);

		// r[2] != r[3] ensures the order of `root` had enough factors of two
		// to compute an fft this size. This fails if a value of LG is used
		// that is larger than the maximum supported for the prime
		assert(r[2] != r[3]);

		// This is just a sanity check
		assert(r[1] == r[2]);
	}

	// Convolve the vectors as and bs, using v as auxilliary space
	vector<ll> cv(const vector<ll> &as, const vector<ll> &bs, u *v) {
		int c=max(SZ(as)+SZ(bs)-1, 0), n=1;
		assert(c <= (1<<LG));
		u h=u(uu(-p)*-p%p), a=m(h, p/2+1), x, y;

		// This computes both the padded size and the
		// renormalization factor h = 2^-n
		while(n<c) n*=2, h=red(p, m(h, a));

		forn(i, n)
			v[i] = i<SZ(as) ? u(as[i]) : 0,
			v[i+n] = i<SZ(bs) ? u(bs[i]) : 0;

		// The forward fft is done using DIF and inverse fft
		// is done using DIT. This is to avoid having to apply
		// the bit-reversal permutation to the vector

		// Decimation-in-frequency FFT
		for(auto s:{v,v+n})
		dforsn(j, 2, n+1) for(int k=j&-j; k/=2;) forsn(i, j-k, j)
			x=s[i], y=s[i-k],
			s[i-k] = red(2*p, x+y),
			s[i] = m(2*p+y-x, r[3*k-j+i]);

		forn(i, n) v[i] = m(v[i], v[i+n]);

		// Decimation-in-time FFT
		forsn(j, 2, n+1) for(int k=1; !(k&j); k*=2) forsn(i, j-k, j)
			x = m(v[i], r[3*k+j-i]),
			y = red(2*p, v[i-k]),
			v[i-k]=x+y, v[i]=2*p+y-x;

		// Renormalize and reduce back to the [0, p) range
		forn(i, c) v[i] = red(p, m(v[i], h));

		return vector<ll>(v, v+c);
	}
};

// For modular convolutions modulo 998244353:
vector<ll> conv_small(const vector<ll> &as, const vector<ll> &bs) {
	static uint32_t v[2<<LG];
	static FFT<uint32_t, uint64_t, 998244353, 3> fft;
	return fft.cv(as, bs, v);
}

// For modular convolutions modulo a 62 bit prime:
vector<ll> conv_big(const vector<ll> &as, const vector<ll> &bs) {
	static uint64_t v[2<<LG];
	static FFT<uint64_t, __uint128_t, (1ull<<62)-(18ull<<32)+1, 3> fft;
	return fft.cv(as, bs, v);
}

// For modular convolutions modulo an arbitrary 32-bit modulus:
vector<ll> conv_sunzi(const vector<ll> &v1, const vector<ll> &v2, ll m) {
	const uint64_t inv = 2703402103339935109ull,
		mod1 = (1ull<<62)-(18ull<<32)+1,
		mod2 = (1ull<<62)-(76ull<<32)+1;
	static_assert(__uint128_t(inv)*mod2%mod1 == (-mod1)%mod1);
	static uint64_t v[2<<LG];
	static FFT<uint64_t, __uint128_t, mod1, 3> fft1;

	// 17 is not actually a primitive root, but its order
	// has enough two-adicity
	static FFT<uint64_t, __uint128_t, mod2, 17> fft2;

	auto as=fft1.cv(v1, v2, v), bs=fft2.cv(v1, v2, v);
	forn(i, SZ(as)) {
		// The Sunzi theorem is used to reconstruct
		// the value modulo m given the value modulo the two
		// large primes mod1 and mod2
		auto d = fft1.m(mod1+as[i]-bs[i], inv);
		d -= mod1*(d >= mod1); d %= m;
		as[i] = (bs[i] + mod2%m*d)%m;
	}
	return as;
}
