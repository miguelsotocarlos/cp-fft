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

// AC: https://codeforces.com/gym/104288/problem/G
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
const int LG = 22;
using ll = long long;

template<class u, class uu, u p, u root>
struct FFT {
	u r[1+(2<<LG)];
	constexpr u m(u a, u b) {
		uu k = uu(a)*b;
		#define op(g) g*(g*p+2)
		k += u(k) * (op(op(op(op(op(-p)))))) * uu(p);
		#undef op
		return u(k>>(8*sizeof(u)));
	}
	constexpr u red(u k, u a) { return a-k*(a>=k); }
	FFT() {
		u k = r[2<<LG] = -p%p, b=root, e = p>>LG;
		for(; e; e/=2, b=m(b,b)) if(e%2) k=m(k, b);
		dforn(i, 2<<LG) r[i]=red(p, m(r[i+1], k)), i&(i-1)?0:k=m(k,k);
		assert(r[2] != r[3]); assert(r[1] == r[2]);
	}
	vector<ll> cv(const vector<ll> &as, const vector<ll> &bs, u *v) {
		int c=max(SZ(as)+SZ(bs)-1, 0), n=1;
		assert(c <= (1<<LG));
		u h=u(uu(-p)*-p%p), a=m(h, p/2+1), x, y;
		while(n<c) n*=2, h=red(p, m(h, a));
		forn(i, n)
			v[i] = i<SZ(as) ? u(as[i]) : 0,
			v[i+n] = i<SZ(bs) ? u(bs[i]) : 0;
		for(auto s:{v,v+n})
		dforsn(j, 2, n+1) for(int k=j&-j; k/=2;) forsn(i, j-k, j)
			x=s[i], y=s[i-k],
			s[i-k] = red(2*p, x+y),
			s[i] = m(2*p+y-x, r[3*k-j+i]);
		forn(i, n) v[i] = m(v[i], v[i+n]);
		forsn(j, 2, n+1) for(int k=1; !(k&j); k*=2) forsn(i, j-k, j)
			x = m(v[i], r[3*k+j-i]),
			y = red(2*p, v[i-k]),
			v[i-k]=x+y, v[i]=2*p+y-x;
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

const int maxn = 1010;
const ll mod1 = 998244353;

int pattern[maxn][maxn];
int board[maxn][maxn];
int res[maxn][maxn];

int main(){
	cin.tie(0);
	cin.sync_with_stdio(0);
	int px, py; cin >> px >> py;
	forn(i, px) forn(j, py) cin >> pattern[i][j];
	int qx, qy; cin >> qx >> qy;
	forn(i, qx) forn(j, qy) cin >> board[i][j];
	forn(b, 7) {
		vector<ll> as(2*maxn*maxn, 0);
		vector<ll> bs(2*maxn*maxn, 0);
		int total = 0;
		forn(i, maxn) forn(j, maxn) {
			if(board[i][j] != 0) as[2*i*maxn + j] = (board[i][j] & (1<<b)) ? 1 : mod1-1;
			if(pattern[i][j] != 0) total++, bs[2*i*maxn + j] = (pattern[i][j] & (1<<b)) ? 1 : mod1-1;
		}
		reverse(bs.begin(), bs.end());
		auto conv = conv_small(as, bs);
		forn(x, qx - px + 1) forn(y, qy - py + 1) {
			if(total == conv[2*maxn * x + y + SZ(bs) - 1]) res[x][y]++;
		}
	}
	vector<pair<int, int>> poss;
	forn(i, maxn) forn(j, maxn) if(res[i][j] == 7) poss.push_back({i, j});
	cout << SZ(poss) << "\n";
	for(auto p : poss) {
		cout << p.first+1 << " " << p.second+1 << "\n";
	}
}
