# My ACM Templates

## 平面最近点对

```cpp
#include <bits/stdc++.h>
#define rep(i, l, r) for (int i = l; i <= r; ++i)
using namespace std;

const int N = 200005;
const double inf = 1e18;
int n;
class point {
   public:
    double x, y;
    point() {}
    point(double _x, double _y) : x(_x), y(_y) {}
    bool operator<(const point &rhs) const { return x < rhs.x; }
    double dis() { return sqrt(x * x + y * y); }
    point operator-(const point &rhs) const {
        return point(x - rhs.x, y - rhs.y);
    }
    void read() { scanf("%lf%lf", &x, &y); }
} a[N];

double dis(point a, point b) { return (a - b).dis(); }

double solve(int l, int r) {
    if (l == r) return inf;
    int mid = (l + r) / 2;
    double mid_x = a[mid].x;
    double dist = min(solve(l, mid), solve(mid + 1, r));
    deque<int> rq;
    int j = mid + 1;
    // search
    rep(i, l, mid) {
        if (a[i].x - mid_x <= -dist) continue;
        while (j <= r && a[j].y - a[i].y < dist) {
            if (a[j].x - mid_x < dist) rq.push_back(j);
            j++;
        }
        while (!rq.empty() && a[rq.front()].y - a[i].y < -dist) {
            rq.pop_front();
        }
        for (auto j : rq) {
            dist = min(dist, dis(a[i], a[j]));
        }
    }
    // merge
    static point b[N];
    int cur = l;
    for (int i = l, j = mid + 1; i <= mid || j <= r;) {
        if (j > r || i <= mid && a[i].y <= a[j].y) {
            b[cur++] = a[i++];
        } else {
            b[cur++] = a[j++];
        }
    }
    rep(i, l, r) a[i] = b[i];

    return dist;
}

int main() {
    scanf("%d", &n);
    rep(i, 1, n) { a[i].read(); }
    sort(a + 1, a + n + 1);
    printf("%.4lf\n", solve(1, n));
    return 0;
}
```

## 半平面交

```cpp
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>
#define rep(i, l, r) for (int i = l; i < r; ++i)
#define N 605
#define eps 1e-10
#define pi acos(-1.0)
using namespace std;
int n, m;
class point {
   public:
    double x, y;
    point(double _x = 0, double _y = 0) : x(_x), y(_y) {}
    point operator+(const point &rhs) const {
        return point(x + rhs.x, y + rhs.y);
    }
    point operator-(const point &rhs) const {
        return point(x - rhs.x, y - rhs.y);
    }
    point operator*(double k) const { return point(x * k, y * k); }
    void read() { scanf("%lf%lf", &x, &y); }
} p[N], a[N];
double cross(point p1, point p2) { return p1.x * p2.y - p2.x * p1.y; }
class line {
   public:
    point p, v;  // point+vector
    line() {}
    line(point _p, point _v) : p(_p), v(_v) {}
    bool operator<(const line &rhs) const {
        return atan2(v.y, v.x) < atan2(rhs.v.y, rhs.v.x);
    }
} l[N], q[N];
int cl = 0;                                  // count line
void newl(const point &a, const point &b) {  // new line ab
    l[++cl] = line(a, b - a);
}
bool left(const point &p, const line &l) {  // point left of vector
    return cross(l.v, p - l.p) > 0;
}
point pos(const line &a, const line &b) {  // intersection
    point x = a.p - b.p;
    double t = cross(b.v, x) / cross(a.v, b.v);  // t = -h_x/h_(a.v)
    return a.p + a.v * t;  // a.p+cross(a.p-b.p, b.v)/cross(a.v, b.v)*a.v
}
double halfplane() {
    sort(l + 1, l + cl + 1);
    int h = 1, t = 1;  // head, tail
    q[1] = l[1];
    rep(i, 2, cl + 1) {
        while (h < t && !left(p[t - 1], l[i])) --t;  // tail erased
        while (h < t && !left(p[h], l[i])) ++h;      // head erased
        if (fabs(cross(q[t].v, l[i].v)) < eps)
            q[t] = left(q[t].p, l[i]) ? q[t] : l[i];  // same direction
        else
            q[++t] = l[i];
        if (h < t) p[t - 1] = pos(q[t], q[t - 1]);  // intersection point
    }
    while (h < t && !left(p[t - 1], q[h])) --t;
    p[t] = pos(q[t], q[h]);  // last intersection
    if (t - h <= 1) return 0;
    double ans = 0;
    rep(i, h, t) ans += cross(p[i], p[i + 1]);
    return (ans + cross(p[t], p[h])) / 2;
}
int main() {
    scanf("%d", &n);
    while (n--) {
        scanf("%d", &m);
        rep(i, 0, m) a[i].read();  // counter-clockwise
        rep(i, 0, m) newl(a[i], a[(i + 1) % m]);
    }
    printf("%.3f", halfplane());
    return 0;
}
```

## 带花树

```cpp
#include <bits/stdc++.h>
#define rep(i, l, r) for (int i = l; i <= r; ++i)
using namespace std;
const int N = 1005;
const int M = 50005;
int n, m;
vector<int> E[N];

int tag[N];  //-1: none, 0: out, 1: in
int match[N], pre[N];
int vis[N];
int f[N];
void init() {
    rep(i, 1, n) {
        match[i] = 0;
        vis[i] = 0;
    }
}

int lca(int x, int y) {
    static int cnt = 0;
    ++cnt;  // time
    x = f[x], y = f[y];
    while (vis[x] != cnt) {
        if (x) {
            vis[x] = cnt;
            x = f[pre[match[x]]];
        }
        swap(x, y);
    }
    return x;  // double visited
}

queue<int> q;
void blossom(int x, int y, int fa) {  // from x to fa
    while (f[x] != fa) {
        pre[x] = y;  // bidirectional pre
        y = match[x];
        if (tag[y] == 1) {  // new outlet
            tag[y] = 0;
            q.push(y);
        }
        f[x] = f[y] = fa;
        x = pre[y];
    }
}

bool aug(int s) {
    rep(i, 1, n) { tag[i] = -1, f[i] = i; }
    q = queue<int>();
    q.push(s);
    tag[s] = 0;
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        for (auto v : E[u]) {
            if (tag[v] == -1) {
                pre[v] = u;
                tag[v] = 1;
                if (!match[v]) {  // found road
                    for (int from = u, to = v; to; from = pre[to]) {  // augment
                        match[to] = from;
                        swap(match[from], to);
                    }
                    return true;
                } else {  // matched
                    tag[match[v]] = 0;
                    q.push(match[v]);
                }
            } else if (tag[v] == 0 && f[u] != f[v]) {  // blossom
                int fa = lca(u, v);
                blossom(u, v, fa);
                blossom(v, u, fa);
            }
        }
    }
    return false;
}
int main() {
    scanf("%d%d", &n, &m);
    rep(i, 1, n) { E[i].clear(); }
    rep(i, 1, m) {
        int u, v;
        scanf("%d%d", &u, &v);
        E[u].push_back(v);
        E[v].push_back(u);
    }

    init();
    int ans = 0;
    rep(i, 1, n) {
        if (!match[i]) ans += aug(i);
    }
    printf("%d\n", ans);
    rep(i, 1, n) printf("%d ", match[i]);
    return 0;
}
```

## 类Euclid算法

求的是 $f_{a,b,c,n}(p, q)=\sum_{i=0}^n\sum_{j=1}^{[\frac{ai+b}{c}]}i^pj^q$

```cpp
#pragma GCC optimize(3)
#include <bits/stdc++.h>
#define rep(i, l, r) for (int i = l; i <= r; ++i)
using namespace std;
const int N = 125;
const int mod = 998244353;
int C[N][N];
int B[N];
bool vis[N * 3][N][N];
int dp[N * 3][N][N];
bool vs[N][N];
int ss[N][N];
int inv[N];
int A[N][N];
void init() {
   int lim = 120;
   C[0][0] = 1;
   rep(i, 1, lim) {
       C[i][0] = 1;
       rep(j, 1, i) { C[i][j] = (C[i - 1][j - 1] + C[i - 1][j]) % mod; }
   }

   inv[1] = 1;
   rep(i, 2, lim) { inv[i] = 1ll * (mod - mod / i) * inv[mod % i] % mod; }

   B[0] = 1;
   rep(i, 1, lim) {
       int tmp = 0;
       rep(j, 0, i - 1) { tmp = (tmp + 1ll * C[i + 1][j] * B[j]) % mod; }
       B[i] = (-1ll * tmp * inv[i + 1] % mod + mod) % mod;
       // printf("%d ", B[i]);
   }
   // printf("\n");

   rep(i, 0, lim) {
       rep(j, 0, i) {
           A[i][j] = 1ll * inv[i + 1] * C[i + 1][j] % mod * B[j] % mod;
           if (j & 1) A[i][j] = (mod - A[i][j]) % mod;
       }
   }
}

int S(int p, int n, int d = 0) {
   if (d && vs[d][p]) return ss[d][p];
   int ret = 0;
   rep(j, 0, p) { ret = (1ll * ret * n + A[p][j]) % mod; }
   ret = 1ll * ret * n % mod;
   // printf("S(%d,%d)=%d\n", p, n, ret);
   vs[d][p] = true;
   return ss[d][p] = ret;
}

int f(int a, int b, int c, int n, int p, int q, int d, int step) {
   if (vis[step][p][q]) return dp[step][p][q];
   int ret;
   if (b >= c) {
       int k = b / c, res = b - c * k;
       int s = 1ll * (S(p, n, d) + !p) * S(q, k) % mod;
       if (!a) return s;
       int s1 = 0;
       rep(r, 0, q) {
           int tmp = 1ll * C[q][r] * f(a, res, c, n, p, r, d, step + 1) % mod;
           s1 = (1ll * s1 * k + tmp) % mod;
       }
       ret = (s + s1) % mod;
   } else if (a >= c) {
       int k = a / c, res = a - c * k;
       int s = 0;
       rep(j, 0, q) {
           int tmp = 1ll * A[q][j] * S(p + q + 1 - j, n, d) % mod;
           s = (1ll * s * k + tmp) % mod;
       }
       s = 1ll * s * k % mod;
       int s1 = 0;
       rep(r, 0, q) {
           int tmp = 1ll * C[q][r] *
                     f(res, b, c, n, p + q - r, r, d, step + 1) % mod;
           s1 = (1ll * s1 * k + tmp) % mod;
       }
       ret = (s + s1) % mod;
   } else {
       int m = (1ll * a * n + b) / c;
       if (!m) return 0;

       int s = 1ll * S(p, n, d) * S(q, m) % mod;
       int s1 = 0;
       rep(r, 0, q) {
           s1 = (s1 + 1ll * C[q][r] *
                          f(c, c - b - 1, a, m - 1, r, p, d + 1, step + 1)) %
                mod;
       }
       ret = (s - s1 + mod) % mod;
       // printf("(%dx+%d)/%d %d (%d,%d): %d(%d)\n", a, b, c, n, p, q, ret, d);
   }
   vis[step][p][q] = true;

   return dp[step][p][q] = ret;
}

int main() {
   init();
   int a, b, c, p, q, n;
   scanf("%d%d%d%d%d%d", &a, &b, &c, &p, &q, &n);
   printf("%d", f(a, b, c, n, p, q, 1, 0));
   return 0;
}
```

## SAM

ms 存储的是只出现一次的串的长度

```cpp
#include<bits/stdc++.h>
#define rep(i, l, r) for(int i=l; i<=r; ++i)
#define N 1000006
using namespace std;
int n;
char s[N];
multiset<int> ms;
class SAM {    
    public:
        class state {
            public:
                int len; // length of longest string
                int link; // longest nonequivalent suffix
                int nxt[26]; // next states (default 0)
                bool vis;
        };
        state st[N*2]; // states
        int sz; // size
        int last; // last state

        void init() {
            last = 0;
            st[0].len = 0; // empty state
            st[0].link = -1;
            st[0].vis = false;
            sz = 1;
        }

        void extend(char c) {
            int cur = sz++; // new state
            st[cur].vis = false;
            st[cur].len = st[last].len + 1; // length
            int p = last;
            while(p != -1 && !st[p].nxt[c - 'a']) {
                st[p].nxt[c - 'a'] = cur;
                p = st[p].link;
            }
            if(p == -1) {
                st[cur].link = 0;
            } else {
                int q = st[p].nxt[c - 'a'];
                if(st[p].len + 1 == st[q].len) { // no other state
                    st[cur].link = q;
                    if(!st[q].vis){
                        st[q].vis=true;
                        ms.erase(ms.find(st[st[q].link].len+1));
                    }
                } else {
                    int clone = sz++;
                    st[clone] = st[q];
                    st[clone].len = st[p].len + 1;
                    st[clone].vis = true;
                    if(!st[q].vis) ms.erase(ms.find(st[st[q].link].len+1));
                    while(p != -1 && st[p].nxt[c - 'a'] == q){
                        st[p].nxt[c - 'a'] = clone; // transfer
                        p = st[p].link;
                    }
                    st[q].link = st[cur].link = clone;
                    if(!st[q].vis) ms.insert(st[st[q].link].len+1);
                }
            }
            last = cur;
            ms.insert(st[st[cur].link].len+1);
        }
};
SAM a;
int main() {
    scanf("%d%s", &n, s+1);
    a.init();
    rep(i, 1, n){
        a.extend(s[i]);
        printf("%d\n", *ms.begin());
    }
    return 0;
}
```

## 多项式齐次线性递推

```cpp
#include <bits/stdc++.h>
#define rep(i, l, r) for (int i = l; i <= r; ++i)
#define per(i, r, l) for (int i = r; i >= l; --i)

using namespace std;
const int mod = 998244353;
typedef vector<int> vi;
int add(int x, int y) { return (x + y) % mod; }
int sub(int x, int y) { return (x - y + mod) % mod; }
int mul(int x, int y) { return 1ll * x * y % mod; }
int pw(int x, int y) {
    int ret = 1;
    while (y) {
        if (y & 1) ret = mul(ret, x);
        x = mul(x, x);
        y >>= 1;
    }
    return ret;
}
int inv(int x) { return pw(x, mod - 2); }

vi rev;
void get_rev(int n) {
    static int lim = -1;
    if (n == lim) return;
    lim = n;
    rev.resize(n);
    int bit = 0;
    while ((1 << bit) < n) ++bit;
    rev[0] = 0;
    rep(i, 1, n - 1) { rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (bit - 1)); }
}

void DFT(vi &a, int n, int dir = 1) {  // n = 2^k
    get_rev(n);
    rep(i, 0, n - 1) {
        if (i < rev[i]) swap(a[i], a[rev[i]]);
    }
    for (int len = 1; (len << 1) <= n; len <<= 1) {
        int wn = pw(3, (mod - 1) / (len << 1));
        if (dir == -1) wn = inv(wn);
        for (int i = 0; i < n; i += len << 1) {
            int w = 1;
            rep(j, i, i + len - 1) {
                int tmp = mul(a[j + len], w);
                a[j + len] = sub(a[j], tmp);
                a[j] = add(a[j], tmp);
                w = mul(w, wn);
            }
        }
    }
    if (dir == -1) {
        int invn = inv(n);
        rep(i, 0, n - 1) { a[i] = mul(a[i], invn); }
    }
}

vi operator+(vi a, vi b) {
    if (a.size() < b.size()) a.resize(b.size());
    rep(i, 0, (int)b.size() - 1) { a[i] = add(a[i], b[i]); }
    return a;
}
vi operator-(vi a, vi b) {
    if (a.size() < b.size()) a.resize(b.size());
    rep(i, 0, (int)b.size() - 1) { a[i] = sub(a[i], b[i]); }
    return a;
}
vi operator*(vi a, vi b) {
    int n = a.size() + b.size() - 1;
    int lim = 1;
    while (lim < n) lim <<= 1;
    a.resize(lim);
    b.resize(lim);
    DFT(a, lim);
    DFT(b, lim);
    rep(i, 0, lim - 1) { a[i] = mul(a[i], b[i]); }
    DFT(a, lim, -1);
    return a;
}

vi inv(vi a) {
    int lim = 1;
    while (lim < a.size()) lim <<= 1;
    vi b(1, inv(a[0]));
    for (int len = 2; len <= lim; len <<= 1) {
        vi x(a);
        x.resize(len);
        x.resize(len << 1);
        b.resize(len << 1);
        DFT(x, len << 1);
        DFT(b, len << 1);
        rep(i, 0, (len << 1) - 1) {
            b[i] = (2ll - 1ll * x[i] * b[i] % mod + mod) * b[i] % mod;
        }
        DFT(b, len << 1, -1);
        b.resize(len);
    }
    return b;
}

vi operator/(vi a, vi b) {  // highest not zero
    if (a.size() < b.size()) return vi(1, 0);
    int l = a.size() - b.size() + 1;
    reverse(a.begin(), a.end());
    reverse(b.begin(), b.end());
    a.resize(l);
    b.resize(l);
    vi c = a * inv(b);
    c.resize(l);
    reverse(c.begin(), c.end());
    return c;
}

vi operator%(vi a, vi b) {
    vi r = a - a / b * b;
    r.resize(max((int)b.size() - 1, 1));
    return r;
}

vi dw(vi a) {
    int n = a.size();
    rep(i, 0, n - 2) a[i] = mul(a[i + 1], i + 1);
    a.resize(n - 1);
    return a;
}
vi up(vi a) {
    int n = a.size();
    a.resize(n + 1);
    per(i, n, 1) a[i] = mul(a[i - 1], inv(i));
    a[0] = 0;
    return a;
}
vi ln(vi a) { return up(dw(a) * inv(a)); }
vi exp(vi a) {
    int lim = 1;
    while (lim < a.size()) lim <<= 1;
    vi b(1, 1);
    for (int len = 2; len <= lim * 2; len <<= 1) {  // I don't know why lim*2
        b = b * (a + vi(1, 1) - ln(b));
        b.resize(len);
    }
    return b;
}
vi pw(vi a, int k, vi b) {
    vi ret(1, 1);
    while (k) {
        if (k & 1) ret = ret * a % b;
        a = a * a % b;
        k >>= 1;
    }
    // for (auto v : a) printf("%d ", v);
    return ret;
}
int linear(vi g, vi a, int n) {
    // for (auto v : g) printf("%d ", v);
    int k = g.size() - 1;
    vi t{0, 1};
    vi r = pw(t, n, g);
    int ret = 0;
    rep(i, 0, k - 1) { ret = (ret + mul(r[i], a[i])) % mod; }
    return ret;
}

void test_linear() {
    int n, k;
    scanf("%d%d", &n, &k);
    int x;
    vi f(1, 0);
    rep(i, 1, k) {
        scanf("%d", &x);
        f.push_back((x + mod) % mod);
    }
    vi a;
    rep(i, 0, k - 1) {
        scanf("%d", &x);
        a.push_back((x + mod) % mod);
    }
    vi g;
    per(i, k, 1) { g.push_back(sub(0, f[i])); }
    g.push_back(1);
    int ret = linear(g, a, n);
    printf("%d", ret);
}
int main() {
    test_linear();
    return 0;
}
```

## 三模 NTT

```cpp
#include <bits/stdc++.h>
#define rep(i, l, r) for (int i = l; i <= r; ++i)
#define per(i, r, l) for (int i = r; i >= l; --i)

using namespace std;

const int mods[3] = {998244353, 1004535809, 469762049};
typedef vector<int> vi;
int add(int x, int y, int mod) { return (x + y) % mod; }
int sub(int x, int y, int mod) { return (x - y + mod) % mod; }
int mul(int x, int y, int mod) { return 1ll * x * y % mod; }
int pw(int x, int y, int mod) {
    int ret = 1;
    while (y) {
        if (y & 1) ret = mul(ret, x, mod);
        x = mul(x, x, mod);
        y >>= 1;
    }
    return ret;
}
int inv(int x, int mod) { return pw(x, mod - 2, mod); }

vi rev;
void get_rev(int n) {
    static int lim = -1;
    if (n == lim) return;
    lim = n;
    rev.resize(n);
    int bit = 0;
    while ((1 << bit) < n) ++bit;
    rev[0] = 0;
    rep(i, 1, n - 1) { rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (bit - 1)); }
}

void DFT(vi &a, int n, int mod, int dir = 1) {  // n=2^k, k<=21
    get_rev(n);
    rep(i, 0, n - 1) {
        if (i < rev[i]) swap(a[i], a[rev[i]]);
    }
    for (int len = 1; (len << 1) <= n; len <<= 1) {
        int wn = pw(3, (mod - 1, mod) / (len << 1), mod);
        if (dir == -1) wn = inv(wn, mod);
        for (int i = 0; i < n; i += len << 1) {
            int w = 1;
            rep(j, i, i + len - 1) {
                int tmp = mul(a[j + len], w, mod);
                a[j + len] = sub(a[j], tmp, mod);
                a[j] = add(a[j], tmp, mod);
                w = mul(w, wn, mod);
            }
        }
    }
    if (dir == -1) {
        int invn = inv(n, mod);
        rep(i, 0, n - 1) { a[i] = mul(a[i], invn, mod); }
    }
}

vi add(vi a, vi b, int mod) {  // any
    if (a.size() < b.size()) a.resize(b.size());
    rep(i, 0, (int)b.size() - 1) { a[i] = add(a[i], b[i], mod); }
    return a;
}
vi sub(vi a, vi b, int mod) {  // any
    if (a.size() < b.size()) a.resize(b.size());
    rep(i, 0, (int)b.size() - 1) { a[i] = sub(a[i], b[i], mod); }
    return a;
}

vi mul(vi a, vi b, int mod) {  // in mods
    int n = a.size() + b.size() - 1;
    int lim = 1;
    while (lim < n) lim <<= 1;
    a.resize(lim);
    b.resize(lim);
    DFT(a, lim, mod);
    DFT(b, lim, mod);
    rep(i, 0, lim - 1) { a[i] = mul(a[i], b[i], mod); }
    DFT(a, lim, mod, -1);
    return a;
}

vi mul_(vi a, vi b, int mod) {  // any
    static int mod_ = 0;
    static __int128_t mul_mods = (__int128_t)mods[0] * mods[1] * mods[2], A[3];
    if (mod != mod_) {
        mod_ = mod;
        A[0] = (__int128_t)mods[1] * mods[2] *
               inv(1ll * mods[1] * mods[2] % mods[0], mods[0]);
        A[1] = (__int128_t)mods[2] * mods[0] *
               inv(1ll * mods[2] * mods[0] % mods[1], mods[1]);
        A[2] = (__int128_t)mods[0] * mods[1] *
               inv(1ll * mods[0] * mods[1] % mods[2], mods[2]);
    }
    vi c[3];
    for (int j = 0; j <= 2; ++j) c[j] = mul(a, b, mods[j]);
    int n = c[0].size();
    vi ret;
    for (int i = 0; i < n; ++i) {
        __int128_t now = 0;
        for (int j = 0; j <= 2; ++j) {
            // printf("%d ", c[j][i]);
            now = (now + A[j] * c[j][i]) % mul_mods;
        }
        now = now % mod;
        // printf("\n");
        ret.push_back(now);
    }
    return ret;
}

vi inv(vi a, int mod) {  // prime
    int lim = 1;
    while (lim < a.size()) lim <<= 1;
    vi b(1, inv(a[0], mod));
    for (int len = 2; len <= lim; len <<= 1) {
        vi x(a);
        x.resize(len);
        b = mul_(sub(vi(1, 2), mul_(x, b, mod), mod), b, mod);
        b.resize(len);
    }
    return b;
}

void test_mul() {
    int n, m, p;
    scanf("%d%d%d", &n, &m, &p);
    vi a, b;
    int x;
    rep(i, 0, n) {
        scanf("%d", &x);
        a.push_back(x);
    }
    rep(i, 0, m) {
        scanf("%d", &x);
        b.push_back(x);
    }
    vi c = mul_(a, b, p);
    rep(i, 0, n + m) printf("%d ", c[i]);
}

void test_inv() {
    int n, p = 1e9 + 7;
    scanf("%d", &n);
    vi a;
    int x;

    rep(i, 0, n - 1) {
        scanf("%d", &x);
        a.push_back(x);
    }
    vi b = inv(a, p);
    rep(i, 0, n - 1) printf("%d ", b[i]);
}

int main() {
    assert(mods[0] - 1 == 119 << 23);
    assert(mods[1] - 1 == 479 << 21);
    assert(mods[2] - 1 == 7 << 26);
    // test_mul();
    test_inv();
    return 0;
}
```

## ISAP

```cpp
#include <bits/stdc++.h>
#define rep(i, l, r) for (int i = l; i <= r; i++)
#define per(i, r, l) for (int i = r; i >= l; --i)
using namespace std;
typedef long long ll;
typedef vector<int> vi;
const int N = 55 + 5;
const int M = 205 + 10 + 55;
const int K = 10004;
const int inf = 1e9;
const int mod = 998244353;

int s, t;
int hd[N], last[N];
class edge {
   public:
    int to, r, nxt;
} e[M * 2];
void add(int u, int v, int r, int i) {
    e[i].to = v;
    e[i].r = r;
    e[i].nxt = hd[u];
    hd[u] = i;
}
int ind = 1;
void add2(int u, int v, int r) {
    add(u, v, r, ++ind);
    add(v, u, 0, ++ind);
}

ll dis[M];
int gap[M];
void bfs() {
    queue<int> q;
    rep(i, 1, t) dis[i] = inf, gap[i] = 0, last[i] = hd[i];
    gap[dis[t] = 1] = 1;
    q.push(t);
    while (!q.empty()) {
        int u = q.front();
        // printf("%d ",u);
        q.pop();
        for (int i = hd[u]; i; i = e[i].nxt) {
            if (!e[i ^ 1].r) continue;
            int v = e[i].to;
            if (dis[v] == inf) {
                gap[dis[v] = dis[u] + 1]++;
                q.push(v);
            }
        }
    }
}

ll dfs(int u, ll rest) {
    // printf("%d %lld\n", u, rest);
    if (u == t) return rest;
    int ret = 0;
    for (int &i = last[u]; i; i = e[i].nxt) {
        int v = e[i].to;
        if (e[i].r && dis[v] == dis[u] - 1) {
            int f = dfs(v, min(rest, 1ll * e[i].r));
            e[i].r -= f, e[i ^ 1].r += f;
            rest -= f, ret += f;
        }
        if (!rest) return ret;
    }
    if (--gap[dis[u]] == 0) dis[s] = t + 1;
    ++gap[++dis[u]];
    last[u] = hd[u];
    return ret;
}

ll flow() {
    bfs();
    ll ret = 0;
    while (dis[s] <= t + 1) {
        ret += dfs(s, inf);
    }
    return ret;
}
int n, m;
int V[N];
int num_cdn;
int cdn_node[6];
int cdn_son[N];
struct origin_edge {
    int u, v, c;
    void read() { scanf("%d%d%d", &u, &v, &c); }
    void add_edge(int x) {
        int new_u = (V[u] == 1) ? cdn_son[u] : u;
        add2(new_u, v, c / x);
    }
} oe[M];
bool check(int x) {
    // printf("%d:\n", x);
    rep(status, 0, (1 << num_cdn) - 1) {
        int goal = 0;
        rep(i, 1, t) { hd[i] = 0; }
        ind = 1;
        rep(i, 1, m) { oe[i].add_edge(x); }
        rep(i, 1, num_cdn) {
            int u = cdn_node[i];
            int v = cdn_son[u];
            add2(u, v, inf);
            if (status & (1 << i - 1)) {
                add2(u, t, 1);
                add2(s, v, inf);
                goal++;
            }
        }
        rep(u, 1, n) {
            if (V[u] == 2) {
                add2(u, t, 1);
                goal++;
            }
        }
        ll f = flow();
        // printf("%d %d %lld\n", status, goal, f);
        if (f == goal) return true;
    }
    return false;
}
void work() {
    scanf("%d%d", &n, &m);
    rep(i, 2, n) { scanf("%d", &V[i]); }
    rep(i, 1, m) { oe[i].read(); }

    num_cdn = 0;
    rep(i, 1, n) {
        if (V[i] == 1) {
            ++num_cdn;
            cdn_node[num_cdn] = i;
            cdn_son[i] = n + num_cdn;  // maybe useless
            // printf("!%d->%d\n", i, cdn_son[i]);
        } else {
            cdn_son[i] = -1;
        }
    }
    assert(num_cdn <= 5);
    s = 1;
    t = n + num_cdn + 1;
    int l = 0, r = 1e6;
    while (l < r) {
        int mid = (l + r + 1) / 2;
        if (check(mid))
            l = mid;
        else
            r = mid - 1;
    }
    printf("%d\n", r);
}
int main() {
    int T;
    scanf("%d", &T);
    while (T--) {
        work();
    }
    return 0;
}
```
