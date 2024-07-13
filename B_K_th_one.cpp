#include <bits/stdc++.h>
#define ll long long int
#define f first
#define vi vector<long long int> 
#define arrn int arr[n]
#define ld long double
using namespace std;
#define sortv sort(v.begin(),v.end());
#define allv v.begin(),v.end()
#define loop(p,k,n) for(ll p=k;p<n;p++)
#define revloop(i,n) for(ll i=n-1;i>=0;i--)
#define pb push_back
#define pq priority_queue
#define pll pair<ll,ll>
#define dispv for(int i=0;i<v.size();i++)cout<<v[i]<<" ";
#define vsz v.size()
#define ok {cout<<"YES"<<endl; }
#define sp {<<" "<< }
#define no {cout<<"NO"<<endl; }
typedef ll item;
struct segtree{
    item NEUTRAL_ELEMENT=0;
        ll size;
        vector<item> sums;
        item single(ll v){
            
            return v;
        }
        item merge(item a,item b)
        {
           return a+b;   
        }
        void init(ll n)
        {
            size=1;
            while(size<n)
                size*=2;
            sums.assign(2*size,0ll);
        }
        void set(ll i,ll v,ll x,ll lx,ll rx)
        {
            if(rx-lx==1)
            {
                sums[x]=single(v);
                return ;
            }
            int m=(lx+rx)/2;
            if(i<m)
                set(i,v,2*x+1,lx,m);
            else
                set(i,v,2*x+2,m,rx);
            sums[x]=merge(sums[2*x+1],sums[2*x+2]);
        }
        void set(ll i,ll v)
        {
            set(i,v,0,0,size);
        }
        item sum(ll l,ll r,ll x,ll lx,ll rx)
        {
            if(r<=lx||rx<=l)
                return NEUTRAL_ELEMENT;
            if(l<=lx&&r>=rx)
                return sums[x];
            int m=(lx+rx)/2;
            item sumleft=sum(l,r,2*x+1,lx,m);
            item sumright=sum(l,r,2*x+2,m,rx);
            return merge(sumleft,sumright);
        }
        item sum(ll l,ll r)
        {
            return sum(l,r,0,0,size);
        }
        void build(vector<ll>&a,ll x,ll lx,ll rx)
        {
            if(rx-lx==1)
            {
                if(lx<a.size())
                    sums[x]=single(a[lx]);
                return;
            }
            int m=(lx+rx)/2;
            build(a,2*x+1,lx,m);
            build(a,2*x+2,m,rx);
            sums[x]=merge(sums[2*x+1],sums[2*x+2]);
        }
        void build(vector<ll>&a)
        {
            build(a,0,0,size);
        }
        ll find(ll k,ll x,ll lx,ll rx)
        {
            if(rx-lx==1)return lx;
            ll sl=sums[2*x+1];
            ll m=(lx+rx)/2;
            if(k<sl)
            {
                return find(k,2*x+1,lx,m);
            }
            else
            return find(k-sl,2*x+2,m,rx);
        }
        ll find(ll k)
        {
            return find(k,0,0,size);
        }
    };
//Aliases
using lld= long double;
using ull= unsigned long long;
//Constants
const lld pi= 3.141592653589793238;
const ll mod=1e9+7;
ll binpow(ll a,ll b) {
    ll ans = 1;
    while(b > 0) {
        if((b & 1) == 1) ans *= a;
        a *= a;
        b = b >> 1;
    }
    return ans;
}
 
ll binpowmod(ll a,ll b) {
    ll ans = 1;
    while(b > 0) {
        if((b & 1) == 1) {
            ans *= a;
            ans %= mod;
        }
        a *= a;
        a %= mod;
        b = b >> 1;
    }
    return ans;
}
 
 
const ll MAX = 2e5 + 7;
vector<ll> fact(MAX);
 
ll mult(ll a, ll b) {
    return ((a * b) % mod);
}
 
ll inv(ll a) {
    return binpowmod(a, mod - 2);
}
 
ll divide(ll a, ll b) {
    return mult(a, inv(b));
}
 
ll nCr(ll n, ll r) {
    if(n < r) return 0;
    return divide(fact[n], mult(fact[r], fact[n - r]));
}
 
void preFactorial() {
    fact[0] = 1;
    for(ll i = 1; i < MAX; i++) fact[i] = mult(i, fact[i - 1]);
}
bool issorted(ll n, vector<ll> &v) {
    for(ll i = 1; i < n; i++) {
    if(v[i] < v[i - 1]) return false;
    }
    return true;
}
 
bool isvowel(char c) {
    if(c == 'a' || c == 'e' || c == 'i' || c == 'o' || c == 'u') return true;
    if(c == 'A' || c == 'E' || c == 'I' || c == 'O' || c == 'U') return true;
    return false;
}
 
ll sumv(ll n, vector<ll> &arr) {
    ll ans = 0;
    for(ll i = 0; i < n; i++) ans += arr[i];
    return ans;
}
 
ll countSetBits(ll n) {
    ll ans = 0;
    while(n) {
        ans++;
        n = n & (n - 1);
    }
    return ans;
}
 
vector<ll> primefactors(ll n) {
    vector<ll> factors;
    for(ll i = 2; i * i <= n; i++) {
        ll cnt = 0; 
        while(n % i == 0) {
            cnt++;
            n /= i;
            factors.push_back(i);
        }
    }
    if(n > 1) factors.push_back(n);
    return factors;
}
 
void factors(vector<ll> &v, int n)
{
    for (int i = 1; i * i <= n; i++)
    {
        if (n % i == 0)
        {
            if (i * i == n)
            {
                v.push_back(i);
            }
            else
            {
                v.push_back(i);
                v.push_back(n / i);
            }
        }
    }
}
 
bool ispalindrome(string s) {
    ll i = 0;
    ll j = s.size() - 1;
    while(i <= j) {
        if(s[i] != s[j]) return false;
        i++;
        j--;
    }
    return true;
}
// Mathematical functions
ll gcd(ll a, ll b){if (b == 0)return a;return gcd(b, a % b);} //__gcd 
ll lcm(ll a, ll b){return (a/gcd(a,b)*b);}
ll moduloMultiplication(ll a,ll b,ll mod){ll res = 0;a %= mod;while (b){if (b & 1)res = (res + a) % mod;b >>= 1;}return res;}
ll powermod(ll x, ll y, ll p){ll res = 1;x = x % p;if (x == 0) return 0;while (y > 0){if (y & 1)res = (res*x) % p;y = y>>1;x = (x*x) % p;}return res;}
//Check
bool isprime(ll n){if(n<=1)return false;if(n<=3)return true;if(n%2==0||n%3==0)return false;for(int i=5;i*i<=n;i=i+6)if(n%i==0||n%(i+2)==0)return false;return true;}
bool ispoweroftwo(ll n){if(n==0)return false;return (ceil(log2(n)) == floor(log2(n)));}
bool isperfectsquare(ll x){if (x >= 0) {ll sr = sqrt(x);return (sr * sr == x);}return false;}
 
//Bits
string dectobinary(ll n){string s="";ll i = 0;while (n > 0) {s =to_string(n % 2)+s;n = n / 2;i++;}return s;}
ll binarytodecimal(string n){string num = n;ll dec_value = 0;int base = 1;int len = num.length();for(int i = len - 1; i >= 0; i--){if (num[i] == '1')dec_value += base;base = base * 2;}return dec_value;}
 

// KMP_ALGOOOO BY MUUJIC
// returns the longest proper prefix array of pattern p
// where lps[i]=longest proper prefix which is also suffix of p[0...i]
vector<ll> build_lps(string p) {
  ll sz = p.size();vector<ll> lps;lps.assign(sz + 1, 0); ll j = 0;lps[0] = 0;
  for(ll i = 1; i < sz; i++) {while(j >= 0 && p[i] != p[j]) {if(j >= 1) j = lps[j - 1];else j = -1;}
    j++;lps[i] = j;}return lps;}
// returns matches in vector ans in 0-indexed
void kmp(vector<ll> lps, string s, string p,vector<ll>&ans) {ll psz = p.size(), sz = s.size(); ll j = 0;
  for(ll i = 0; i < sz; i++) { while(j >= 0 && p[j] != s[i]) if(j >= 1) j = lps[j - 1]; else j = -1;j++; if(j == psz) {   j = lps[j - 1];
      // pattern found in string s at position i-psz+1
     ans.push_back(i - psz + 1);}
    // after each loop we have j=longest common suffix of s[0..i] which is also prefix of p
}}
vector<ll> solve_kmp(string s,string p){  vector<ll>ans;vector<ll>lps = build_lps(p);kmp(lps, s, p,ans); return ans;}

//z_fxn_algo by muujic
vector<ll> z_function(string s) {ll n = (ll) s.length();vector<ll> z(n);
  for (ll i = 1, l = 0, r = 0; i < n; ++i) { if (i <= r) z[i] = min (r - i + 1, z[i - l]);
    while (i + z[i] < n && s[z[i]] == s[i + z[i]])
      ++z[i]; if (i + z[i] - 1 > r) l = i, r = i + z[i] - 1;}return z;}
//Speed
#define code ios_base::sync_with_stdio(false);
#define by cin.tie(NULL);
#define muujic cout.tie(NULL);
//Sorting
bool sorta(const pair<ll,ll> &a,const pair<ll,ll> &b){return (a.second < b.second);}
bool sortd(const pair<ll,ll> &a,const pair<ll,ll> &b){return (a.second > b.second);}
//auto pr = std::max_element(m.begin(), m.end(), [](const auto &x, const auto &y) {return x.second < y.second;});//
// ----------- Code Starts Here ----------- //
void luv_by_muujic()
{
  ll n,m;
        cin>>n>>m;
        segtree st;
        st.init(n);
        vector<ll> arr(n);
        loop(i,0,n)
            cin>>arr[i];
        st.build(arr);
       // cout<<st.sum(0,n).seg<<endl;
        ll op,x,y;
        while(m--)
        {
            cin>>op;
            if(op%2)
            {   ll i;cin>>i;
                arr[i]=1-arr[i];
                st.set(i,arr[i]); 
            }
            else  
            {ll k;cin>>k;
              cout<< st.find(k)<<endl;
            }
                
           
        }
}
int main()
{
code by muujic
// int test;cin>>test;
// while(test--){luv_by_muujic();}
 luv_by_muujic();
return 0;
}
