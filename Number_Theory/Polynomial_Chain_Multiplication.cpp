/*
 Polynomial Chain Multiplication
 takes input c[i]-> coefficients of polynomials
 computes (x-c0)(x-c0)(x-c0)..(x-c_{n-1})
 returns ans[] -> coefficients of resultant polynomial
 O(N*log^2(N))
 https://discuss.codechef.com/t/countway-editorial/13246/2
*/

#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace std;
using namespace __gnu_pbds;

#define ll long long
#define db long double
#define ii pair<int,int>
#define vi vector<ll>
#define fi first
#define se second
#define sz(a) (int)(a).size()
#define all(a) (a).begin(),(a).end()
#define pb push_back
#define mp make_pair
#define FN(i, n) for (int i = 0; i < (int)(n); ++i)
#define FEN(i,n) for (int i = 1;i <= (int)(n); ++i)
#define rep(i,a,b) for(int i=a;i<b;i++)
#define repv(i,a,b) for(int i=b-1;i>=a;i--)
#define SET(A, val) memset(A, val, sizeof(A))
typedef tree<int ,null_type,less<int>,rb_tree_tag,tree_order_statistics_node_update>ordered_set ;
// order_of_key (val): returns the no. of values less than val
// find_by_order (k): returns the kth largest element.(0-based)
#define TRACE
#ifdef TRACE
#define trace(...) __f(#__VA_ARGS__, __VA_ARGS__)
	template <typename Arg1>
	void __f(const char* name, Arg1&& arg1){
		cerr << name << " : " << arg1 << std::endl;
	}
	template <typename Arg1, typename... Args>
	void __f(const char* names, Arg1&& arg1, Args&&... args){
		const char* comma = strchr(names + 1, ','); cerr.write(names, comma - names) << " : " << arg1<<" | ";__f(comma+1, args...);
	}
#else
#define trace(...)
#endif

const db PI = acos((db)-1) ;
namespace fft
{
    template < class T>
    class comp
    {
        public:
            T real,img;
            comp (T a = (T)0, T b = (T)0) : real(a), img(b) {}

            comp conj() { return comp(this->real, -(this->img));  }
            comp operator = (const comp& a) { this->real = a.real, this->img = a.img; return *this; }
            comp operator + (const comp& b) { return comp(this->real + b.real, this->img + b.img); }
            comp operator - (const comp& b) { return comp(this->real - b.real, this->img - b.img); }
            comp operator * (const T& num) { return comp(this->real * num, this->img * num); }
            comp operator / (const T& num) { return comp(this->real / num, this->img / num); }
            comp operator * (const comp& b)
            {
                return comp(this->real * b.real - this->img * b.img, this->img * b.real + this->real * b.img);
            }
            comp operator / (const comp& b)
            {
                comp temp(b.real, -b.img) ; comp n = (*this) * temp;
                return n / (b.x * b.x + b.y * b.y);
            }
    };
    #define cd comp<double>
    vector<cd> w ;
    vi rev ;
    void revbits(int newlim)
    {
        static int lim = -1;
        int t, j;
        if(newlim == lim) return ;
        lim = newlim ; rev.resize(lim+1) ;
        int k = 0;
		while ((1 << k) < newlim) ++k;
		assert((1 << k) == newlim);
        FEN(i,lim)
        {
            j = rev[i-1];
            t = k-1;
            while(t >= 0 && ((j >> t) & 1) )  j ^= (1<<t) , --t ;
            if(t>=0) j ^= (1<<t) , --t ;
            rev[i] = j;
        }
    }

    void fft(vector<cd> &poly,int inv = false)
    {
        int len, l;
        revbits(sz(poly));
        if(inv) for(auto &x:poly) x = x.conj() ;

        FN(i,sz(poly)) if(rev[i] > i) swap(poly[i], poly[rev[i]]);
        cd u,v ;
        if(sz(w) < sz(poly)) w.resize(sz(poly)) ;
        for(len = 2, l = 1; len <= sz(poly); len +=len, l +=l)
        {
            if(w[l].real == 0 && w[l].img == 0)
            {
                double ang = PI / l;
                cd ww(cos(ang),sin(ang));
                if(l > 1)
                {
                    for(int j = 0; j < l; ++j)
                    {
                        if(j & 1) w[l + j] = w[(l + j) >> 1] * ww;
                        else w[l + j] = w[(l + j) >> 1];
                    }
                }
                else w[l] = cd(1.0,0.0);
            }

            for(int i = 0; i < sz(poly); i += len) FN(j,l)
            {
                u = poly[i + j], v = poly[i + j + l] * w[l + j] ;
                poly[i + j] = u + v , poly[i + j + l] = u - v ;
            }
        }
        if(inv) for(auto &x:poly) x = x / sz(poly) ;
    }
  vi multiply(vi &a,vi &b)
  {
      int bits = 1, sz1 = sz(a) + sz(b) , reqsz;
      while( (1 << bits) < sz1 ) ++bits;
      reqsz = (1 << bits);
      vector<cd> poly(reqsz) ;
      FN(i,sz(a)) poly[i].real = a[i] ;
      FN(i,sz(b)) poly[i].img = b[i] ;
      fft(poly);
      cd p , qtmp , q ;
      poly[0] = poly[0].real * poly[0].img ; // for i = 0
      FEN(i,reqsz>>1)
      {
          p = poly[i] + poly[reqsz - i].conj() , qtmp = poly[reqsz-i] - poly[i].conj();
          q.real = qtmp.img, q.img = qtmp.real ;
          poly[i] = (p * q) * 0.25;
          if(i) poly[reqsz-i] = poly[i].conj();
      }
      fft(poly,true);
      vi ans(sz1-1);
      FN(i,sz(ans)) ans[i] = (ll)(poly[i].real + 0.5) %13313;

      /*Uncomment for multiplication of two numbers
      int carry = 0;
      for (int i=0; i<(int)(ans.size()); ++i)
      {
        ans[i] += carry;
        carry = ans[i] / 10;
        ans[i] %= 10;
      }
      */
      return ans;
  }
}
int c[200005];
vi vals;
vi polynomial_chain_multiplication(int l,int r)
{
  if(l==r)
    {
      vi tmp(vals[l]+1,1);
      return tmp;
    }
  int mid=(l+r)>>1;
  vi left=polynomial_chain_multiplication(l,mid);
  vi right=polynomial_chain_multiplication(mid+1,r);
  vi ans=fft::multiply(left,right);
  return ans;
}
int main()
{
	std::ios::sync_with_stdio(false);
	cin.tie(NULL) ; cout.tie(NULL) ;
  int n,k,x;
  cin>>n>>k;
  rep(i,1,n+1)
    {
      cin>>x; c[x]++;
    }

  rep(i,1,200001)
    if(c[i])
      vals.pb(c[i]);
  vi ans=polynomial_chain_multiplication(0,sz(vals)-1);
  cout<<ans[k]<<endl;
  return 0;
}
