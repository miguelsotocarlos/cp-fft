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

#include "fft.cpp"
#include <random>

template<class F>
void test(mt19937 rng, int a, int b, ll mod, F f) {
	uniform_int_distribution<ll> dist(0, mod);
	vector<ll> as(a), bs(b);
	forn(i, a) as[i] = dist(rng);
	forn(i, b) bs[i] = dist(rng);
	vector<ll> conv(a+b-1, 0);
	forn(i, a) forn(j, b) {
		conv[i+j] = ll((conv[i+j] + __uint128_t(as[i])*bs[j]) % mod);
	}
	auto computed = f(as, bs);
	assert(conv == computed);
}

template<class F>
void test_all_sizes(mt19937 rng, ll mod, F f) {
	test(rng, 1, 1, mod, f);
	test(rng, 1, 2, mod, f);
	test(rng, 1, 10, mod, f);
	test(rng, 10, 1, mod, f);
	test(rng, 1000, 1, mod, f);
	test(rng, 1000, 1000, mod, f);
	test(rng, 10000, 10000, mod, f);
	test(rng, 50000, 1, mod, f);
	test(rng, 50000, 50000, mod, f);
	test(rng, 200000, 200000, mod, f);
}

int main() {
	mt19937 rng;
	test_all_sizes(rng, 998244353, conv_small);
	test_all_sizes(rng, (1ull<<62)-(18ull<<32)+1, conv_big);
	for(ll mod : {1ll<<31, ll(1e9), ll(1e9)+7}) {
		test_all_sizes(rng, mod, [&](vector<ll> as, vector<ll> bs) {
			return conv_sunzi(as, bs, mod);
		});
	}
}
