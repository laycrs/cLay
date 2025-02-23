// cLay ゲームにっき（仮）別館（仮） http://rsujskf.s602.xrea.com/

// This code is licensed under CC0.
// http://creativecommons.org/publicdomain/zero/1.0/

#include<bits/stdc++.h>
using namespace std;

#define REP(i,a,b) for(i=a;i<b;i++)
#define rep(i,n) REP(i,0,n)

#define mygc(c) (c)=getchar_unlocked()
#define mypc(c) putchar_unlocked(c)

#define ll long long
#define ull unsigned ll

void reader(int *x){int k,m=0;*x=0;for(;;){mygc(k);if(k=='-'){m=1;break;}if('0'<=k&&k<='9'){*x=k-'0';break;}}for(;;){mygc(k);if(k<'0'||k>'9')break;*x=(*x)*10+k-'0';}if(m)(*x)=-(*x);}
void reader(ll *x){int k,m=0;*x=0;for(;;){mygc(k);if(k=='-'){m=1;break;}if('0'<=k&&k<='9'){*x=k-'0';break;}}for(;;){mygc(k);if(k<'0'||k>'9')break;*x=(*x)*10+k-'0';}if(m)(*x)=-(*x);}
void reader(double *x){scanf("%lf",x);}
int reader(char c[]){int i,s=0;for(;;){mygc(i);if(i!=' '&&i!='\n'&&i!='\r'&&i!='\t'&&i!=EOF) break;}c[s++]=i;for(;;){mygc(i);if(i==' '||i=='\n'||i=='\r'||i=='\t'||i==EOF) break;c[s++]=i;}c[s]='\0';return s;}
template <class T, class S> void reader(T *x, S *y){reader(x);reader(y);}
template <class T, class S, class U> void reader(T *x, S *y, U *z){reader(x);reader(y);reader(z);}
template <class T, class S, class U, class V> void reader(T *x, S *y, U *z, V *w){reader(x);reader(y);reader(z);reader(w);}

void writer(int x, char c){int s=0,m=0;char f[10];if(x<0)m=1,x=-x;while(x)f[s++]=x%10,x/=10;if(!s)f[s++]=0;if(m)mypc('-');while(s--)mypc(f[s]+'0');mypc(c);}
void writer(ll x, char c){int s=0,m=0;char f[20];if(x<0)m=1,x=-x;while(x)f[s++]=x%10,x/=10;if(!s)f[s++]=0;if(m)mypc('-');while(s--)mypc(f[s]+'0');mypc(c);}
void writer(double x, char c){printf("%.15f",x);mypc(c);}
void writer(const char c[]){int i;for(i=0;c[i]!='\0';i++)mypc(c[i]);}
void writer(const char x[], char c){int i;for(i=0;x[i]!='\0';i++)mypc(x[i]);mypc(c);}
template<class T> void writerLn(T x){writer(x,'\n');}
template<class T, class S> void writerLn(T x, S y){writer(x,' ');writer(y,'\n');}
template<class T, class S, class U> void writerLn(T x, S y, U z){writer(x,' ');writer(y,' ');writer(z,'\n');}
template<class T> void writerArr(T x[], int n){int i;if(!n){mypc('\n');return;}rep(i,n-1)writer(x[i],' ');writer(x[n-1],'\n');}

char memarr[17000000]; void *mem = memarr;
#define MD 1000000007

ull modpow(ull a, ull b, ull m){
  ull r = 1;
  while(b){
    if(b&1) r = r*a%m;
    b>>=1;
    a = a*a%m;
  }
  return r;
}

ll primitiveRoot(ll p){
  ll i, r, m = p-1;
  if(p<=1) return -1;
  for(i=2;i*i<=p;i++) if(p%i==0) return -1;
  for(r=1;r<p;r++){
    for(i=1;i*i<=m;i++) if(m%i==0){
      if(i*i!=m && modpow(r,i,p)==1) break;
      if(i!=1 && modpow(r,m/i,p)==1) break;
    }
    if(i*i>m) return r;
  }
  return -1;
}


int user_code = 1;
set<string> g_flags;

struct insertfunctions{
  vector<string> name;
  std::set<string> doit, already;
  map<string,string> func, place, parent;
  map<string,vector<string> > need;
  map<string,vector<string> > del;

  void set(){
    {
      string n = "__int128_t";
      string c = "";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    {
      string n = "__uint128_t";
      string c = "";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "optimize";
      string c = "#pragma GCC optimize(\"Ofast\")\n#pragma GCC optimize(\"unroll-loops\")\n#pragma GCC optimize(\"inline\")\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }


    {
      string n = "stdc";
      string c = "#include<bits/stdc++.h>\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "sys_time";
      string c = "#include<sys/time.h>\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "chrono";
      string c = "#include<chrono>\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "namespace";
      string c = "using namespace std;\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "BoostMultiprecision";
      string c = "#include <boost/multiprecision/cpp_int.hpp>\nusing namespace boost::multiprecision;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }


    {
      string n = "define_MD";
      string c = "#define MD (1000000007U)\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "define_for_Mint";
      string c = "#define MINT_W (32U)\n#define MINT_R (294967268U)\n#define MINT_RR (582344008U)\n#define MINT_MDNINV (2226617417U)\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "define_MD_PRIMITIVE_ROOT";
      string c = "#define MD_PRIMITIVE_ROOT (5U)\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "define_PI";
      string c = "#define PI 3.14159265358979323846\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "cLtraits_identity";
      string c = "template<class T>\nstruct cLtraits_identity { using type = T; };\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "cLtraits_try_make_signed";
      string c = "template<class T>\nusing cLtraits_try_make_signed =\n  typename conditional<\n    is_integral<T>::value,\n    make_signed<T>,\n    cLtraits_identity<T>\n    >::type;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"cLtraits_identity");
      need[n] = d;
    }
    {
      string n = "cLtraits_try_make_unsigned";
      string c = "template<class T>\nusing cLtraits_try_make_unsigned =\n  typename conditional<\n    is_integral<T>::value,\n    make_unsigned<T>,\n    cLtraits_identity<T>\n    >::type;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"cLtraits_identity");
      need[n] = d;
    }
    {
      string n = "cLtraits_common_type";
      string c = "template <class S, class T> struct cLtraits_common_type {\n  using tS = typename cLtraits_try_make_signed<S>::type;\n  using tT = typename cLtraits_try_make_signed<T>::type;\n  using type = typename common_type<tS,tT>::type;\n};\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"cLtraits_identity");
      d.push_back((string)"cLtraits_try_make_signed");
      need[n] = d;
    }


    {
      string n = "workmemory";
      string c = "inplace_L void *wmem; char memarr[96000000];";
      string p = "first";
      vector<string> d;

      d.push_back((string)"workmemory_init");
      
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "workmemory_init";
      string c = "wmem = memarr;";
      string p = "main_first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "min_L";
      string c = "template<class S, class T>\ninline auto min_L(S a, T b)\n-> typename cLtraits_common_type<S,T>::type\n{\n  return (typename cLtraits_common_type<S,T>::type) a <= (typename cLtraits_common_type<S,T>::type) b ? a : b;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"cLtraits_common_type");
      need[n] = d;
    }
    {
      string n = "max_L";
      string c = "template<class S, class T>\ninline auto max_L(S a, T b)\n-> typename cLtraits_common_type<S,T>::type\n{\n  return (typename cLtraits_common_type<S,T>::type) a >= (typename cLtraits_common_type<S,T>::type) b ? a : b;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"cLtraits_common_type");
      need[n] = d;
    }


    {
      string n = "chmin";
      string c = "template<class S, class T> inline S chmin(S &a, T b){if(a>b)a=b;return a;}";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "chmax";
      string c = "template<class S, class T> inline S chmax(S &a, T b){if(a<b)a=b;return a;}";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "fDiv";
      string c = "template<class T, class S>\nT inline fDiv(T a, S b){\n  T m;\n  if(b < 0) a = -a, b = -b;\n  m = a % b;\n  if(m == 0) return a / b;\n  if(m < 0) m += b;\n  return (a - m) / b;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    
    {
      string n = "cDiv";
      string c = "template<class T, class S>\nT inline cDiv(T a, S b){\n  T m;\n  if(b < 0) a = -a, b = -b;\n  m = a % b;\n  if(m == 0) return a / b;\n  if(m < 0) m += b;\n  return (a + b - m) / b;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    
    {
      string n = "divup";
      string c = "template<class S, class T> inline S divup_L(S a, T b){ return (a+b-1)/b; }";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "moddw";
      string c = "template<class S, class T> inline S moddw_L(S a, const T b){\n  a %= b;\n  if(a < 0) a += b;\n  return a;\n}\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    

    {
      string n = "walloc1d";
      string c = "template<class T>\ninline void walloc1d(T **arr, int x, void **mem = &wmem){\n  static int skip[16] = {0, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1};\n  (*mem) = (void*)( ((char*)(*mem)) + skip[((ull)(*mem)) & 15] );\n  (*arr)=(T*)(*mem);\n  (*mem)=((*arr)+x);\n}\ntemplate<class T>\ninline void walloc1d(T **arr, int x1, int x2, void **mem = &wmem){\n  walloc1d(arr, x2-x1, mem);\n  (*arr) -= x1;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"workmemory");
      need[n] = d;
    }
    {
      string n = "walloc2d";
      string c = "template<class T>\ninline void walloc2d(T ***arr, int x, int y, void **mem = &wmem){\n  int i;\n  walloc1d(arr, x, mem);\n  rep(i,x) walloc1d(&((*arr)[i]), y, mem);\n}\ntemplate<class T>\ninline void walloc2d(T ***arr, int x1, int x2, int y1, int y2, void **mem = &wmem){\n  int i;\n  walloc1d(arr, x1, x2, mem);\n  rep(i,x1,x2) walloc1d(&((*arr)[i]), y1, y2, mem);\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
    }

    {
      string n = "malloc1d";
      string c = "template<class T>\nvoid malloc1d(T **arr, int x){\n  (*arr) = (T*)malloc(x*sizeof(T));\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "free1d";
      string c = "template<class T>\nvoid free1d(T *arr){\n  free(arr);\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "malloc2d";
      string c = "template<class T>\nvoid malloc2d(T ***arr, int x, int y){\n  int i;\n  (*arr) = (T**)malloc(x*sizeof(T*));\n  (*arr)[0] = (T*)malloc(x*y*sizeof(T));\n  REP(i,1,x)(*arr)[i]=(*arr)[i-1]+y;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "free2d";
      string c = "template<class T>\nvoid free2d(T **arr){\n  free(arr[0]);\n  free(arr);\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }


    {
      string n = "isPrime_head";
      string c = "#define ISPRIME_PRE_CALC_SIZE 1000000\nchar isPrime_prime_table[ISPRIME_PRE_CALC_SIZE];\ntemplate<class T> inline int isPrime(T n);\nvoid isPrime32_init(void);\nint isPrime32_sub(int b, unsigned n);\nint isPrime32(unsigned n);\nint isPrime64_sub(ll b, ull n);\nint isPrime64(ull n);\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "isPrime_init";
      string c = "{\n  isPrime32_init();\n}\n";
      string p = "main_first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "isPrime";
      string c = "template<class T>\ninline int isPrime(T n){\n  T i;\n  if(n<=1) return 0;\n  if(n <= (1ULL<<32) - 1) return isPrime32(n);\n  if(n <= (1ULL<<63) - 1 + (1ULL<<63)) return isPrime64(n);\n  if(n<=3) return 1;\n  if(n%2==0) return 0;\n  for(i=3;i*i<=n;i+=2) if(n%i==0) return 0;\n  return 1;\n}\nint isPrime32_sub(int b, unsigned n){\n  unsigned i, t = 0, u = n-1;\n  ull nw, nx;\n  while(!(u&1)) t++, u >>= 1;\n  nw = 1;\n  nx = b % n;\n  while(u){\n    if(u&1) nw = (nw * nx) % n;\n    nx = (nx * nx) % n;\n    u >>= 1;\n  }\n  rep(i,t){\n    nx = (nw * nw) % n;\n    if(nx == 1 && nw != 1 && nw != n-1) return 0;\n    nw = nx;\n  }\n  if(nw == 1) return 1;\n  return 0;\n}\nint isPrime32(unsigned n){\n  if(n < 1d5) return isPrime_prime_table[n];\n  if(n % 2 == 0) return 0;\n  if(!isPrime32_sub(2,n)) return 0;\n  if(n<=1000000){\n    if(!isPrime32_sub(3,n)) return 0;\n  } else {\n    if(!isPrime32_sub(7,n)) return 0;\n    if(!isPrime32_sub(61,n)) return 0;\n  }\n  return 1;\n}\nint isPrime64_sub(ll b, ull n){\n  ull i, t = 0, u = n-1;\n  __uint128_t nw, nx;\n  while(!(u&1)) t++, u >>= 1;\n  nw = 1;\n  nx = b % n;\n  while(u){\n    if(u&1) nw = (nw * nx) % n;\n    nx = (nx * nx) % n;\n    u >>= 1;\n  }\n  rep(i,t){\n    nx = (nw * nw) % n;\n    if(nx == 1 && nw != 1 && nw != n-1) return 0;\n    nw = nx;\n  }\n  if(nw == 1) return 1;\n  return 0;\n}\nint isPrime64(ull n){\n  if(n < 1d5) return isPrime_prime_table[n];\n  if(n < (1ULL<<32)) return isPrime32(n);\n  if(n % 2 == 0) return 0;\n  if(!isPrime64_sub(2,n)) return 0;\n  if(n <= 21652684502221ULL){\n    if(!isPrime64_sub(1215,n)) return 0;\n    if(!isPrime64_sub(34862,n)) return 0;\n    if(!isPrime64_sub(574237825,n)) return 0;\n  } else {\n    if(!isPrime64_sub(325,n)) return 0;\n    if(!isPrime64_sub(9375,n)) return 0;\n    if(!isPrime64_sub(28178,n)) return 0;\n    if(!isPrime64_sub(450775,n)) return 0;\n    if(!isPrime64_sub(9780504,n)) return 0;\n    if(!isPrime64_sub(1795265022,n)) return 0;\n  }\n  return 1;\n}\nvoid isPrime32_init(void){\n  int i, j, k;\n  k = Isqrt_f(ISPRIME_PRE_CALC_SIZE);\n  rep(i,2,ISPRIME_PRE_CALC_SIZE) isPrime_prime_table[i] = 1;\n  rep(i,2,k+1) if(isPrime_prime_table[i]) rep(j,i*i,ISPRIME_PRE_CALC_SIZE,i) isPrime_prime_table[j] = 0;\n}\n";
      string p = "last";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"isPrime_head");
      d.push_back((string)"isPrime_init");
      d.push_back((string)"Isqrt_f");
      need[n] = d;
    }
    {
      string n = "Factor_head";
      string c = "#define FACTOR_PRE_CALC_SIZE 1000000\nint factor_hasprime_table[FACTOR_PRE_CALC_SIZE];\ntemplate<class T, class R1, class R2> int Factor(T N, R1 fac[], R2 fs[], void *mem = wmem);\ntemplate<class T, class R1> int Factor(T N, R1 fac[], void *mem = wmem);\ntemplate<class T> int Factor(T N, void *mem = wmem);\nunsigned Factor32_rho(unsigned n);\ntemplate<class R1, class R2> int Factor32(unsigned N, R1 fac[], R2 fs[], void *mem = wmem);\null Factor64_rho(ull n);\ntemplate<class R1, class R2> int Factor64(ull N, R1 fac[], R2 fs[], void *mem = wmem);\nvoid Factor32_init(void);\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "Factor_init";
      string c = "{\n  Factor32_init();\n}\n";
      string p = "main_first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "Factor";
      string c = "template<class T, class R1, class R2>\nint Factor(T N, R1 fac[], R2 fs[], void *mem/* = wmem*/){\n  T i;\n  int sz = 0;\n  if(N <= 1) return sz;\n  if(N <= (1ULL<<32) - 1) return Factor32(N, fac, fs, mem);\n  if(N <= (1ULL<<63) - 1 + (1ULL<<63)) return Factor64(N, fac, fs, mem);\n  if(N%2==0){\n    fac[sz] = 2;\n    fs[sz] = 1;\n    N /= 2;\n    while(N%2==0){\n      N /= 2;\n      fs[sz]++;\n    }\n    sz++;\n  }\n  for(i=3;i*i<=N;i+=2) if(N%i==0){\n    fac[sz] = i;\n    fs[sz] = 1;\n    N /= i;\n    while(N%i==0){\n      N /= i;\n      fs[sz]++;\n    }\n    sz++;\n  }\n  if(N > 1){\n    fac[sz] = N;\n    fs[sz] = 1;\n    sz++;\n  }\n  return sz;\n}\ntemplate<class T, class R1>\nint Factor(T N, R1 fac[], void *mem/* = wmem*/){\n  int *fs;\n  walloc1d(&fs,128,&mem);\n  return Factor(N, fac, fs, mem);\n}\ntemplate<class T>\nint Factor(T N, void *mem/* = wmem*/){\n  T *fac;\n  int *fs;\n  walloc1d(&fac,128,&mem);\n  walloc1d(&fs,128,&mem);\n  return Factor(N, fac, fs, mem);\n}\nunsigned Factor32_rho(unsigned n){\n  static Rand rnd;\n  const int step = 16;\n  int i, s, nx, mx;\n  ull x, y, memo, c, m;\n  unsigned g;\n  ll lm;\n  lm = min(1ULL<<30, n - 1);\n  for(;;){\n    x = y = rnd.get(1LL, lm);\n    c = rnd.get(1LL, lm);\n    g = 1;\n    for(nx=1;g==1;nx<<=1){\n      x = y;\n      rep(i,nx) y = (y * y + c) % n;\n      for(s=0;s<nx&&g==1;s+=step){\n        m = 1;\n        memo = y;\n        mx = min(step, nx-s);\n        rep(i,mx){\n          y = (y * y + c) % n;\n          if(x >= y) m = (m * (x - y)) % n;\n          else       m = (m * (y - x)) % n;\n        }\n        g = gcd(n, m);\n        if(g != 1){\n          if(g != n) return g;\n          y = memo;\n          for(;;){\n            y = (y * y + c) % n;\n            if(x >= y) m = x - y;\n            else       m = y - x;\n            g = gcd(n, m);\n            if(g == n) break;\n            if(g != 1) return g;\n          }\n        }\n      }\n    }\n  }\n  return 0;\n}\ntemplate<class R1, class R2>\nint Factor32(unsigned N, R1 fac[], R2 fs[], void *mem/* = wmem*/){\n  int res = 0, sz = 0, i, k;\n  unsigned *val, *valtmp, pf, n;\n  if(N <= 1) return 0;\n  walloc1d(&val, 128, &mem);\n  walloc1d(&valtmp, 128, &mem);\n  while(N%2==0) val[res++] = 2, N /= 2;\n  while(N%3==0) val[res++] = 3, N /= 3;\n  while(N%5==0) val[res++] = 5, N /= 5;\n  if(N > 1){\n    valtmp[sz++] = N;\n  }\n  while(sz){\n    while(sz && isPrime32(valtmp[sz-1])){\n      val[res] = valtmp[sz-1];\n      res++;\n      sz--;\n    }\n    if(sz==0) break;\n    n = valtmp[sz-1];\n    if(n < FACTOR_PRE_CALC_SIZE){\n      while(n > 1){\n        val[res++] = factor_hasprime_table[n];\n        n /= factor_hasprime_table[n];\n      }\n      sz--;\n    } else {\n      pf = Factor32_rho(n);\n      valtmp[sz-1] = pf;\n      valtmp[sz] = n / pf;\n      sz++;\n    }\n  }\n  sortA(res, val, mem);\n  k = 0;\n  rep(i,res){\n    if(k && fac[k-1] == val[i]) fs[k-1]++, continue;\n    fac[k] = val[i];\n    fs[k] = 1;\n    k++;\n  }\n  res = k;\n  return res;\n}\null Factor64_rho(ull n){\n  static Rand rnd;\n  const int step = 16;\n  int i, s, nx, mx;\n  __uint128_t x, y, memo, c, m;\n  ull g;\n  ll lm;\n  lm = min(1ULL<<30, n - 1);\n  for(;;){\n    x = y = rnd.get(1LL, lm);\n    c = rnd.get(1LL, lm);\n    g = 1;\n    for(nx=1;g==1;nx<<=1){\n      x = y;\n      rep(i,nx) y = (y * y + c) % n;\n      for(s=0;s<nx&&g==1;s+=step){\n        m = 1;\n        memo = y;\n        mx = min(step, nx-s);\n        rep(i,mx){\n          y = (y * y + c) % n;\n          if(x >= y) m = (m * (x - y)) % n;\n          else       m = (m * (y - x)) % n;\n        }\n        g = gcd(n, m);\n        if(g != 1){\n          if(g != n) return g;\n          y = memo;\n          for(;;){\n            y = (y * y + c) % n;\n            if(x >= y) m = x - y;\n            else       m = y - x;\n            g = gcd(n, m);\n            if(g == n) break;\n            if(g != 1) return g;\n          }\n        }\n      }\n    }\n  }\n  return 0;\n}\ntemplate<class R1, class R2>\nint Factor64(ull N, R1 fac[], R2 fs[], void *mem/* = wmem*/){\n  int res = 0, sz = 0, i, k;\n  ull *val, *valtmp, pf, n;\n  if(N <= 1) return 0;\n  walloc1d(&val, 128, &mem);\n  walloc1d(&valtmp, 128, &mem);\n  while(N%2==0) val[res++] = 2, N /= 2;\n  while(N%3==0) val[res++] = 3, N /= 3;\n  while(N%5==0) val[res++] = 5, N /= 5;\n  if(N > 1){\n    valtmp[sz++] = N;\n  }\n  while(sz){\n    while(sz && isPrime64(valtmp[sz-1])){\n      val[res] = valtmp[sz-1];\n      res++;\n      sz--;\n    }\n    if(sz==0) break;\n    n = valtmp[sz-1];\n    if(n < FACTOR_PRE_CALC_SIZE){\n      while(n > 1){\n        val[res++] = factor_hasprime_table[n];\n        n /= factor_hasprime_table[n];\n      }\n      sz--;\n    } else if(n < (1ULL<<32)){\n      pf = Factor32_rho(n);\n      valtmp[sz-1] = pf;\n      valtmp[sz] = n / pf;\n      sz++;\n    } else {\n      pf = Factor64_rho(n);\n      valtmp[sz-1] = pf;\n      valtmp[sz] = n / pf;\n      sz++;\n    }\n  }\n  sortA(res, val, mem);\n  k = 0;\n  rep(i,res){\n    if(k && fac[k-1] == val[i]) fs[k-1]++, continue;\n    fac[k] = val[i];\n    fs[k] = 1;\n    k++;\n  }\n  res = k;\n  return res;\n}\nvoid Factor32_init(void){\n  int i, j, k;\n  k = Isqrt_f(FACTOR_PRE_CALC_SIZE);\n  rep(i,2,FACTOR_PRE_CALC_SIZE) factor_hasprime_table[i] = i;\n  rep(i,2,k+1) if(factor_hasprime_table[i]==i) rep(j,i*i,FACTOR_PRE_CALC_SIZE,i) factor_hasprime_table[j] = i;\n}\n";
      string p = "last";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Factor_head");
      d.push_back((string)"Factor_init");
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"sortA");
      d.push_back((string)"Rand");
      d.push_back((string)"isPrime");
      d.push_back((string)"Isqrt_f");
      d.push_back((string)"min_L");
      d.push_back((string)"gcd");
      need[n] = d;
    }
    {
      string n = "FactorM_head";
      string c = "template<class T, class R> int FactorM(T N, R fac[], void *mem = wmem);\ntemplate<class T> int FactorM(T N, void *mem = wmem);\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "FactorM";
      string c = "template<class T, class R> int FactorM(T N, R fac[], void *mem/* = wmem*/){\n  int i, j, msz, sz = 0;\n  T *mfac;\n  int *fs;\n  walloc1d(&mfac,128,&mem);\n  walloc1d(&fs,128,&mem);\n  msz = Factor(N, mfac, fs, mem);\n  rep(i,msz) rep(j,fs[i]) fac[sz++] = mfac[i];\n  return sz;\n}\ntemplate<class T> int FactorM(T N, void *mem/* = wmem*/){\n  int i, msz, res = 0;\n  T *mfac;\n  int *fs;\n  walloc1d(&mfac,128,&mem);\n  walloc1d(&fs,128,&mem);\n  msz = Factor(N, mfac, fs, mem);\n  rep(i,msz) res += fs[i];\n  return res;\n}\n";
      string p = "last";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"FactorM_head");
      d.push_back((string)"Factor");
      d.push_back((string)"workmemory");
      d.push_back((string)"walloc1d");
      need[n] = d;
    }
    {
      string n = "Divisor_head";
      string c = "template<class T, class R> int Divisor(T N, R res[], void *mem = wmem);\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "Divisor";
      string c = "template<class T, class R>\nint Divisor(T N, R res[], void *mem/* = wmem*/){\n  int i, j, k, s, sz = 0;\n  T *fc;\n  int *fs, fsz;\n  walloc1d(&fc, 128, &mem);\n  walloc1d(&fs, 128, &mem);\n  \n  fsz = Factor(N, fc, fs, mem);\n  res[sz++] = 1;\n  rep(i,fsz){\n    s = sz;\n    k = s * fs[i];\n    rep(j,k) res[sz++] = res[j] * fc[i];\n  }\n  sort(res, res+sz);\n  return sz;\n}\n";
      string p = "last";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Divisor_head");
      d.push_back((string)"workmemory");
      d.push_back((string)"walloc1d");
      d.push_back((string)"Factor");
      need[n] = d;
    }
    {
      string n = "DivisorSum_head";
      string c = "ll DivisorSum(ll n, void *mem = wmem);\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "DivisorSum";
      string c = "ll DivisorSum(ll n, void *mem/* = wmem*/){\n  int i;\n  ll res, t, s;\n  int fs, *fn;\n  ll *fc;\n  if(n<=0) return 0;\n  walloc1d(&fc, 128, &mem);\n  walloc1d(&fn, 128, &mem);\n  fs = Factor(n, fc, fn, mem);\n  res = 1;\n  rep(i,fs){\n    s = t = 1;\n    rep(j,fn[i]){\n      t *= fc[i];\n      s += t;\n    }\n    res *= s;\n  }\n  return res;\n}\n";
      string p = "last";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"DivisorSum_head");
      d.push_back((string)"workmemory");
      d.push_back((string)"walloc1d");
      d.push_back((string)"Factor");
      need[n] = d;
    }
    {
      string n = "Moebius_head";
      string c = "template<class T>\nint Moebius(T n, void *mem = wmem);\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "Moebius";
      string c = "template<class T>\nint Moebius(T n, void *mem /* = wmem*/){\n  int i;\n  T *fc;\n  int *fs, fsz;\n  walloc1d(&fc, 128, &mem);\n  walloc1d(&fs, 128, &mem);\n  fsz = Factor(n, fc, fs, mem);\n  rep(i,fsz) if(fs[i] > 1) return 0;\n  if(fsz%2) return -1;\n  return 1;\n}\n";
      string p = "last";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Moebius_head");
      d.push_back((string)"workmemory");
      d.push_back((string)"walloc1d");
      d.push_back((string)"Factor");
      need[n] = d;
    }
    {
      string n = "EulerPhi_head";
      string c = "template<class T>\nT EulerPhi(T n, void *mem = wmem);\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "EulerPhi";
      string c = "template<class T>\nT EulerPhi(T n, void *mem /* = wmem*/){\n  int i;\n  T *fc;\n  int *fs, fsz;\n  walloc1d(&fc, 128, &mem);\n  walloc1d(&fs, 128, &mem);\n  fsz = Factor(n, fc, fs, mem);\n  rep(i,fsz) n = n / fc[i] * (fc[i]-1);\n  return n;\n}\n";
      string p = "last";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"EulerPhi_head");
      d.push_back((string)"workmemory");
      d.push_back((string)"walloc1d");
      d.push_back((string)"Factor");
      need[n] = d;
    }

    {
      string n = "Ilog2_f";
      string c = "inline int Ilog2_f_L(const int n){\n  int res;\n  if(n <= 0) return -1;\n  res = sizeof(int) * 8 - __builtin_clz(n) - 1;\n  return res;\n}\ninline int Ilog2_f_L(const ll n){\n  int res;\n  if(n <= 0) return -1;\n  res = sizeof(ll) * 8 - __builtin_clzll(n) - 1;\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "Ilog2_c";
      string c = "inline int Ilog2_c_L(const int n){\n  int res;\n  if(n <= 0) return -1;\n  res = sizeof(int) * 8 - __builtin_clz(n) - 1;\n  if(n != (1<<res)) res++;\n  return res;\n}\ninline int Ilog2_c_L(const ll n){\n  int res;\n  if(n <= 0) return -1;\n  res = sizeof(ll) * 8 - __builtin_clzll(n) - 1;\n  if(n != (1LL<<res)) res++;\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "Ilog2_s";
      string c = "inline int Ilog2_s_L(const int n){\n  int res;\n  if(n <= 0) return -1;\n  res = __builtin_ffs(n) - 1;\n  if(n == (1<<res)) return res;\n  return -1;\n}\ninline int Ilog2_s_L(const ll n){\n  int res;\n  if(n <= 0) return -1;\n  res = __builtin_ffsll(n) - 1;\n  if(n == (1LL<<res)) return res;\n  return -1;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }


    {
      string n = "sortI";
      string c = "template<class T1>\nvoid sortI(int N, T1 a[], void *mem = wmem){\n  sort(a, a+N);\n}\ntemplate<class T1, class T2>\nvoid sortI(int N, T1 a[], T2 b[], void *mem = wmem){\n  int i;\n  pair<T1, T2> *arr;\n  walloc1d(&arr, N, &mem);\n  rep(i,N) arr[i].first = a[i], arr[i].second = b[i];\n  sort(arr, arr+N);\n  rep(i,N) a[i] = arr[i].first, b[i] = arr[i].second;\n}\ntemplate<class T1, class T2, class T3>\nvoid sortI(int N, T1 a[], T2 b[], T3 c[], void *mem = wmem){\n  int i;\n  pair<T1, pair<T2, T3> > *arr;\n  walloc1d(&arr, N, &mem);\n  rep(i,N) arr[i].first = a[i], arr[i].second.first = b[i], arr[i].second.second = c[i];\n  sort(arr, arr+N);\n  rep(i,N) a[i] = arr[i].first, b[i] = arr[i].second.first, c[i] = arr[i].second.second;\n}\ntemplate<class T1, class T2, class T3, class T4>\nvoid sortI(int N, T1 a[], T2 b[], T3 c[], T4 d[], void *mem = wmem){\n  int i;\n  pair<pair<T1, T2>, pair<T3, T4> > *arr;\n  walloc1d(&arr, N, &mem);\n  rep(i,N) arr[i].first.first = a[i], arr[i].first.second = b[i], arr[i].second.first = c[i], arr[i].second.second = d[i];\n  sort(arr, arr+N);\n  rep(i,N) a[i] = arr[i].first.first, b[i] = arr[i].first.second, c[i] = arr[i].second.first, d[i] = arr[i].second.second;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
    }

    {
      string n = "sortA";
      string c = "template<class T> inline int sort_helper_getbit(T A[]){ return -1; }\ntemplate<> inline int sort_helper_getbit(int A[]){ return sizeof(int)*8; }\ntemplate<> inline int sort_helper_getbit(unsigned A[]){ return sizeof(unsigned)*8; }\ntemplate<> inline int sort_helper_getbit(ll A[]){ return sizeof(ll)*8; }\ntemplate<> inline int sort_helper_getbit(ull A[]){ return sizeof(ull)*8; }\ntemplate<> inline int sort_helper_getbit(char A[]){ return sizeof(char)*8; }\ntemplate<class T>\nvoid sortA_1_int_L(int N, T A[], void *mem = wmem){\n  int i, j, k, b, s, ok;\n  ok = 1;\n  rep(i,1,N) if(A[i-1] > A[i]){\n    ok = 0;\n    break;\n  }\n  if(ok) return;\n  if(N < 128){\n    sort(A,A+N);\n    return;\n  }\n  b = sort_helper_getbit(A);\n  if(b==-1){\n    sort(A,A+N);\n    return;\n  }\n  T mn, mx;\n  mn = mx = A[0];\n  rep(i,1,N) mn <?= A[i];\n  rep(i,1,N) mx >?= A[i];\n  ok = 1;\n  if(mn < 0 && mx > 0 && (mn < -N || mx > N)) ok = 0;\n  if(ok && mx - mn > N) ok = 0;\n  if(ok){\n    int *tmp;\n    walloc1d(&tmp, mx-mn+1, &mem);\n    rep(i,mx-mn+1) tmp[i] = 0;\n    rep(i,N) tmp[A[i]-mn]++;\n    k = 0;\n    rep(i,mx-mn+1) while(tmp[i] > 0){\n      tmp[i]--;\n      A[k++] = i+mn;\n    }\n    return;\n  }\n  {\n    typename make_unsigned<T>::type *t[2], mask, cur, one = 1;\n    T tone = 1;\n    int *pos, nn = 0, ss;\n    s = Ilog2_f(N);\n    if(s > 8) s = (8 + (s-7)/2);\n    ss = 1;\n    for(;;){\n      if(ss >= b) break;\n      if( mx >= 0 && (tone << (ss-1)) < mx ) ss++, continue;\n      if( mn < 0 && -(tone << (ss-1)) >= mn ) ss++, continue;\n      break;\n    }\n    k = ss /+ s;\n    s = ss /+ k;\n    mask = 0;\n    rep(i,b) if(i < s*k) mask |= one << i;\n    t[0] = (typename make_unsigned<T>::type *) A;\n    walloc1d(&t[1], N, &mem);\n    walloc1d(&pos, (1<<s)+1, &mem);\n    rep(j,k){\n      cur = 0;\n      rep(i,b) if(s*j <= i && i < s*(j+1) && i < b) cur |= one << i;\n      rep(i,(1<<s)+1) pos[i] = 0;\n      rep(i,N) pos[((t[nn][i]&cur)>>(s*j))+1]++;\n      rep(i,(1<<s)) pos[i+1] += pos[i];\n      rep(i,N) t[nn^1][pos[(t[nn][i]&cur)>>(s*j)]++] = t[nn][i];\n      nn ^= 1;\n      mask ^= cur;\n    }\n    if(mn < 0 && mx >= 0){\n      k = 0;\n      rep(i,N) if(A[i] < 0) k++;\n      rep(i,k) t[nn^1][i] = t[nn][N-k+i];\n      rep(i,k,N) t[nn^1][i] = t[nn][i-k];\n      nn ^= 1;\n    }\n    if(nn==1){\n      rep(i,N) t[0][i] = t[1][i];\n    }\n    return;\n  }\n  sort(A, A+N);\n}\ntemplate<class T>\nvoid sortA_1_nonint_L(int N, T A[], void *mem = wmem){\n  sort(A,A+N);\n}\ntemplate<class T> void sortA_1_call_L(int N, T A[], void *mem = wmem){ sortA_1_nonint_L(N, A, mem); }\ntemplate<> void sortA_1_call_L(int N, int A[], void *mem){ sortA_1_int_L(N, A, mem); }\ntemplate<> void sortA_1_call_L(int N, unsigned A[], void *mem){ sortA_1_int_L(N, A, mem); }\ntemplate<> void sortA_1_call_L(int N, ll A[], void *mem){ sortA_1_int_L(N, A, mem); }\ntemplate<> void sortA_1_call_L(int N, ull A[], void *mem){ sortA_1_int_L(N, A, mem); }\ntemplate<> void sortA_1_call_L(int N, char A[], void *mem){ sortA_1_int_L(N, A, mem); }\ntemplate<class T1>\nvoid sortA(int N, T1 a[], void *mem = wmem){\n  sortA_1_call_L(N, a, mem);\n}\ntemplate<class T1, class T2>\nvoid sortA_2_int_L(int N, T1 A[], T2 B[], void *mem = wmem){\n  int b_a, b_b, s1, s2, so2;\n  T1 mn1, mx1;\n  T2 mn2, mx2;\n  typename cLtraits_try_make_unsigned<T1>::type r1;\n  typename cLtraits_try_make_unsigned<T2>::type r2;\n  so2 = 1;\n  rep(i,1,N){\n    if(A[i-1] > A[i] || (A[i-1]==A[i] && B[i-1] > B[i])){\n      so2 = 0;\n      break;\n    }\n  }\n  if(so2) return;\n  so2 = 1;\n  rep(i,1,N){\n    if(A[i-1] > A[i]){\n      so2 = 0;\n      break;\n    }\n  }\n  if(so2==1){\n    int k = 0;\n    rep(i,1,N) if(A[i] != A[i-1]){\n      sortA_1_call_L(i-k, B+k, mem);\n      k = i;\n    }\n    sortA_1_call_L(N-k, B+k, mem);\n    return;\n  }\n  if(N < 128){\n    sortI(N,A,B,mem);\n    return;\n  }\n  b_a = sort_helper_getbit(A);\n  b_b = sort_helper_getbit(B);\n  if(b_a == -1 || b_b == -1){\n    sortI(N,A,B,mem);\n    return;\n  }\n  mn1 = mx1 = A[0];\n  rep(i,1,N) mn1 <?= A[i];\n  rep(i,1,N) mx1 >?= A[i];\n  mn2 = mx2 = B[0];\n  rep(i,1,N) mn2 <?= B[i];\n  rep(i,1,N) mx2 >?= B[i];\n  if(mn1 < -ll_inf || mn2 < -ll_inf || mx1 > ll_inf || mx2 > ll_inf || mx1-mn1 > ll_inf || mx2-mn2 > ll_inf){\n    sortI(N,A,B,mem);\n    return;\n  }\n  r1 = (typename cLtraits_try_make_unsigned<T1>::type)(mx1) - (typename cLtraits_try_make_unsigned<T1>::type)(mn1);\n  r2 = (typename cLtraits_try_make_unsigned<T2>::type)(mx2) - (typename cLtraits_try_make_unsigned<T2>::type)(mn2);\n  if(r1 == 0){\n    sortA_1_call_L(N, B, mem);\n    return;\n  }\n  if(r2 == 0){\n    sortA_1_call_L(N, A, mem);\n    return;\n  }\n  if(r1 <= N){\n    so2 = 1;\n    rep(i,1,N) if(B[i-1] > B[i]){\n      so2 = 0;\n      break;\n    }\n    if(so2 == 1){\n      T1 *aa; T2 *bb;\n      int *pos, k;\n      walloc1d(&aa,N,&mem);\n      walloc1d(&bb,N,&mem);\n      walloc1d(&pos,r1+2,&mem);\n      rep(i,r1+2) pos[i] = 0;\n      rep(i,N) aa[i] = A[i];\n      rep(i,N) bb[i] = B[i];\n      rep(i,N) pos[(typename cLtraits_try_make_unsigned<T1>::type)((typename cLtraits_try_make_unsigned<T1>::type)aa[i]-(typename cLtraits_try_make_unsigned<T1>::type)mn1)+1]++;\n      rep(i,1,r1+2) pos[i] += pos[i-1];\n      rep(i,N){\n        k = pos[(typename cLtraits_try_make_unsigned<T1>::type)((typename cLtraits_try_make_unsigned<T1>::type)aa[i]-(typename cLtraits_try_make_unsigned<T1>::type)mn1)+0]++;\n        A[k] = aa[i];\n        B[k] = bb[i];\n      }\n      return;\n    }\n  }\n  s1 = s2 = 1;\n  while( s1 < 64 && r1 >= (1ULL<<s1) ) s1++;\n  while( s2 < 64 && r2 >= (1ULL<<s2) ) s2++;\n  if(s1 + s2 <= 32){\n    unsigned *tmp;\n    walloc1d(&tmp,N,&mem);\n    rep(i,N) tmp[i] = (((unsigned)((int)A[i]-(int)mn1)) << s2) | ((unsigned)((int)B[i]-(int)mn2));\n    sortA_1_call_L(N, tmp, mem);\n    rep(i,N){\n      A[i] = ((int)(tmp[i] >> s2)) + ((int)mn1);\n      B[i] = ((int)(tmp[i] & ((1U<<s2)-1))) + ((int)mn2);\n    }\n    return;\n  }\n  if(s1 + s2 <= 64){\n    ull *tmp;\n    walloc1d(&tmp,N,&mem);\n    rep(i,N) tmp[i] = (((ull)((ll)A[i]-(ll)mn1)) << s2) | ((ull)((ll)B[i]-(ll)mn2));\n    sortA_1_call_L(N, tmp, mem);\n    rep(i,N){\n      A[i] = ((ll)(tmp[i] >> s2)) + ((ll)mn1);\n      B[i] = ((ll)(tmp[i] & ((1ULL<<s2)-1))) + ((ll)mn2);\n    }\n    return;\n  }\n  sortI(N,A,B,mem);\n}\ntemplate<class T1, class T2>\nvoid sortA_2_nonint_L(int N, T1 A[], T2 B[], void *mem = wmem){\n  sortI(N,A,B,mem);\n}\ntemplate<class T1, class T2> void sortA_2_call_L(int N, T1 A[], T2 B[], void *mem = wmem){ sortA_2_nonint_L(N, A, B, mem); }\ntemplate<class T2> void sortA_2_call_L(int N, int A[], T2 B[], void *mem){ sortA_2_int_L(N, A, B, mem); }\ntemplate<class T2> void sortA_2_call_L(int N, unsigned A[], T2 B[], void *mem){ sortA_2_int_L(N, A, B, mem); }\ntemplate<class T2> void sortA_2_call_L(int N, ll A[], T2 B[], void *mem){ sortA_2_int_L(N, A, B, mem); }\ntemplate<class T2> void sortA_2_call_L(int N, ull A[], T2 B[], void *mem){ sortA_2_int_L(N, A, B, mem); }\ntemplate<class T2> void sortA_2_call_L(int N, char A[], T2 B[], void *mem){ sortA_2_int_L(N, A, B, mem); }\ntemplate<class T1, class T2>\nvoid sortA(int N, T1 a[], T2 b[], void *mem = wmem){\n  sortA_2_call_L(N, a, b, mem);\n}\ntemplate<class T1, class T2, class T3>\nvoid sortA(int N, T1 a[], T2 b[], T3 c[], void *mem = wmem){\n  int i;\n  pair<T1, pair<T2, T3> > *arr;\n  walloc1d(&arr, N, &mem);\n  rep(i,N) arr[i].first = a[i], arr[i].second.first = b[i], arr[i].second.second = c[i];\n  sort(arr, arr+N);\n  rep(i,N) a[i] = arr[i].first, b[i] = arr[i].second.first, c[i] = arr[i].second.second;\n}\ntemplate<class T1, class T2, class T3, class T4>\nvoid sortA(int N, T1 a[], T2 b[], T3 c[], T4 d[], void *mem = wmem){\n  int i;\n  pair<pair<T1, T2>, pair<T3, T4> > *arr;\n  walloc1d(&arr, N, &mem);\n  rep(i,N) arr[i].first.first = a[i], arr[i].first.second = b[i], arr[i].second.first = c[i], arr[i].second.second = d[i];\n  sort(arr, arr+N);\n  rep(i,N) a[i] = arr[i].first.first, b[i] = arr[i].first.second, c[i] = arr[i].second.first, d[i] = arr[i].second.second;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"chmin");
      d.push_back((string)"chmax");
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"divup");
      d.push_back((string)"Ilog2_f");
      d.push_back((string)"sortI");
      d.push_back((string)"cLtraits_try_make_unsigned");
      need[n] = d;
    }
    {
      string n = "rsortA";
      string c = "template<class T1>\nvoid rsortA(int N, T1 a[], void *mem = wmem){\n  sortA(N, a, mem);\n  reverse(a, a+N);\n}\ntemplate<class T1, class T2>\nvoid rsortA(int N, T1 a[], T2 b[], void *mem = wmem){\n  sortA(N, a, b, mem);\n  reverse(a, a+N);\n  reverse(b, b+N);\n}\ntemplate<class T1, class T2, class T3>\nvoid rsortA(int N, T1 a[], T2 b[], T3 c[], void *mem = wmem){\n  sortA(N, a, b, c, mem);\n  reverse(a, a+N);\n  reverse(b, b+N);\n  reverse(c, c+N);\n}\ntemplate<class T1, class T2, class T3, class T4>\nvoid rsortA(int N, T1 a[], T2 b[], T3 c[], T4 d[], void *mem = wmem){\n  sortA(N, a, b, c, d, mem);\n  reverse(a, a+N);\n  reverse(b, b+N);\n  reverse(c, c+N);\n  reverse(d, d+N);\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"sortA");
      need[n] = d;
    }

    {
      string n = "sortA_index";
      string c = "template<class T1, class T2>\nvoid sortA_index(int N, T1 a[], T2 b[], void *mem = wmem){\n  int i;\n  rep(i,N) b[i] = i;\n  sortA(N,a,b,mem);\n}\ntemplate<class T1, class T2, class T3>\nvoid sortA_index(int N, T1 a[], T2 b[], T3 c[], void *mem = wmem){\n  int i;\n  rep(i,N) c[i] = i;\n  sortA(N,a,b,c,mem);\n}\ntemplate<class T1, class T2, class T3, class T4>\nvoid sortA_index(int N, T1 a[], T2 b[], T3 c[], T4 d[], void *mem = wmem){\n  int i;\n  rep(i,N) d[i] = i;\n  sortA(N,a,b,c,d,mem);\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"sortA");
      need[n] = d;
    }

    {
      string n = "sortV";
      string c = "template<class T1>\nvoid sortV(vector<T1> &a, void *mem = wmem){\n  int i, n = a.size();\n  T1 *aa;\n  walloc1d(&aa, n, &mem);\n  rep(i,n) aa[i] = a[i];\n  sortA(n, aa, mem);\n  rep(i,n) a[i] = aa[i];\n}\ntemplate<class T1, class T2>\nvoid sortV(vector<T1> &a, vector<T2> &b, void *mem = wmem){\n  int i, n = a.size();\n  T1 *aa;\n  T2 *bb;\n  walloc1d(&aa, n, &mem);\n  walloc1d(&bb, n, &mem);\n  rep(i,n) aa[i] = a[i];\n  rep(i,n) bb[i] = b[i];\n  sortA(n, aa, bb, mem);\n  rep(i,n) a[i] = aa[i];\n  rep(i,n) b[i] = bb[i];\n}\ntemplate<class T1, class T2, class T3>\nvoid sortV(vector<T1> &a, vector<T2> &b, vector<T3> &c, void *mem = wmem){\n  int i, n = a.size();\n  T1 *aa;\n  T2 *bb;\n  T3 *cc;\n  walloc1d(&aa, n, &mem);\n  walloc1d(&bb, n, &mem);\n  walloc1d(&cc, n, &mem);\n  rep(i,n) aa[i] = a[i];\n  rep(i,n) bb[i] = b[i];\n  rep(i,n) cc[i] = c[i];\n  sortA(n, aa, bb, cc, mem);\n  rep(i,n) a[i] = aa[i];\n  rep(i,n) b[i] = bb[i];\n  rep(i,n) c[i] = cc[i];\n}\ntemplate<class T1, class T2, class T3, class T4>\nvoid sortV(vector<T1> &a, vector<T2> &b, vector<T3> &c, vector<T4> &d, void *mem = wmem){\n  int i, n = a.size();\n  T1 *aa;\n  T2 *bb;\n  T3 *cc;\n  T4 *dd;\n  walloc1d(&aa, n, &mem);\n  walloc1d(&bb, n, &mem);\n  walloc1d(&cc, n, &mem);\n  walloc1d(&dd, n, &mem);\n  rep(i,n) aa[i] = a[i];\n  rep(i,n) bb[i] = b[i];\n  rep(i,n) cc[i] = c[i];\n  rep(i,n) dd[i] = d[i];\n  sortA(n, aa, bb, cc, dd, mem);\n  rep(i,n) a[i] = aa[i];\n  rep(i,n) b[i] = bb[i];\n  rep(i,n) c[i] = cc[i];\n  rep(i,n) d[i] = dd[i];\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"sortA");
      need[n] = d;
    }
    {
      string n = "rsortV";
      string c = "template<class T1>\nvoid rsortV(vector<T1> &a, void *mem = wmem){\n  sortV(a, mem);\n  reverse(a.begin(), a.end());\n}\ntemplate<class T1, class T2>\nvoid rsortV(vector<T1> &a, vector<T2> &b, void *mem = wmem){\n  sortV(a, b, mem);\n  reverse(a.begin(), a.end());\n  reverse(b.begin(), b.end());\n}\ntemplate<class T1, class T2, class T3>\nvoid rsortV(vector<T1> &a, vector<T2> &b, vector<T3> &c, void *mem = wmem){\n  sortV(a, b, c, mem);\n  reverse(a.begin(), a.end());\n  reverse(b.begin(), b.end());\n  reverse(c.begin(), c.end());\n}\ntemplate<class T1, class T2, class T3, class T4>\nvoid rsortV(vector<T1> &a, vector<T2> &b, vector<T3> &c, vector<T4> &d, void *mem = wmem){\n  sortV(a, b, c, d, mem);\n  reverse(a.begin(), a.end());\n  reverse(b.begin(), b.end());\n  reverse(c.begin(), c.end());\n  reverse(d.begin(), d.end());\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"sortV");
      need[n] = d;
    }


    {
      string n = "sortF_unsigned";
      string c = "void sortF_L(int N, unsigned A[], void *mem = wmem){\n  int i, m, bt;\n  unsigned *arr, c;\n  int *sz;\n\n  if(N < 256){\n    sort(A, A+N);\n    return;\n  }\n\n  bt = sizeof(unsigned) * 8;\n  walloc1d(&arr, N, &mem);\n  walloc1d(&sz, N, &mem);\n\n  for(m=0;m<bt;m+=8){\n    rep(i,257) sz[i] = 0;\n    rep(i,N) sz[ 1+((A[i]>>m)&255u) ]++;\n    rep(i,1,257) sz[i] += sz[i-1];\n    rep(i,N){\n      c = ((A[i]>>m)&255u);\n      arr[sz[c]++] = A[i];\n    }\n    swap(A, arr);\n  }\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "sortF_int";
      string c = "void sortF_L(int N, int A[], void *mem = wmem){\n  int i, x, y, z;\n  int *arr;\n  unsigned *send;\n\n  if(N < 256){\n    sort(A, A+N);\n    return;\n  }\n\n  send = (unsigned*)A;\n  sortF_L(N, send, mem);\n  if(A[0] < 0 || A[N-1] >= 0) return;\n  \n  x = 0;\n  y = N;\n  while(x < y){\n    z = (x+y) / 2;\n    if(A[z] < 0) y = z; else x = z+1;\n  }\n\n  walloc1d(&arr, N, &mem);\n  z = 0;\n  rep(i,x,N) arr[z++] = A[i];\n  rep(i,x) arr[z++] = A[i];\n  rep(i,N) A[i] = arr[i];\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"sortF_unsigned");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "sortF_ull";
      string c = "void sortF_L(int N, ull A[], void *mem = wmem){\n  int i, m, bt;\n  ull *arr;\n  unsigned c;\n  int *sz;\n\n  if(N < 512){\n    sort(A, A+N);\n    return;\n  }\n\n  bt = sizeof(ull) * 8;\n  \n  walloc1d(&arr, N, &mem);\n  walloc1d(&sz, N, &mem);\n\n  for(m=0;m<bt;m+=8){\n    rep(i,257) sz[i] = 0;\n    rep(i,N) sz[ 1+((A[i]>>m)&255u) ]++;\n    rep(i,1,257) sz[i] += sz[i-1];\n    rep(i,N){\n      c = ((A[i]>>m)&255u);\n      arr[sz[c]++] = A[i];\n    }\n    swap(A, arr);\n  }\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "sortF_ll";
      string c = "void sortF_L(int N, ll A[], void *mem = wmem){\n  int i, x, y, z;\n  ll *arr;\n  ull *send;\n\n  if(N < 512){\n    sort(A, A+N);\n    return;\n  }\n\n  send = (ull*)A;\n  sortF_L(N, send, mem);\n  if(A[0] < 0 || A[N-1] >= 0) return;\n  \n  x = 0;\n  y = N;\n  while(x < y){\n    z = (x+y) / 2;\n    if(A[z] < 0) y = z; else x = z+1;\n  }\n\n  walloc1d(&arr, N, &mem);\n  z = 0;\n  rep(i,x,N) arr[z++] = A[i];\n  rep(i,x) arr[z++] = A[i];\n  rep(i,N) A[i] = arr[i];\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"sortF_ull");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "Kth0_size2";
      string c = "template<class T1, class T2>\ninline T1 Kth0_L(const T1 a, const T2 b){\n  if(a <= b) return a;\n  return b;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "Kth1_size2";
      string c = "template<class T1, class T2>\ninline T1 Kth1_L(const T1 a, const T2 b){\n  if(a >= b) return a;\n  return b;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "Kth0_size3";
      string c = "template<class T1, class T2, class T3>\ninline T1 Kth0_L(const T1 a, const T2 b, const T3 c){\n  if(a <= b){\n    if(a <= c) return a;\n    return c;\n  }\n  if(b <= c) return b;\n  return c;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "Kth1_size3";
      string c = "template<class T1, class T2, class T3>\ninline T1 Kth1_L(const T1 a, const T2 b, const T3 c){\n  if(a <= b){\n    if(b <= c) return b;\n    if(c <= a) return a;\n    return c;\n  }\n  if(a <= c) return a;\n  if(c <= b) return b;\n  return c;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "Kth2_size3";
      string c = "template<class T1, class T2, class T3>\ninline T1 Kth2_L(const T1 a, const T2 b, const T3 c){\n  if(a <= b){\n    if(b <= c) return c;\n    return b;\n  }\n  if(a <= c) return c;\n  return a;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "Kth0_size4";
      string c = "template<class T1, class T2, class T3, class T4>\ninline T1 Kth0_L(const T1 a, const T2 b, const T3 c, const T4 d){\n  if(a <= b){\n    if(a <= c){\n      if(a <= d) return a;\n      return d;\n    }\n    if(c <= d) return c;\n    return d;\n  }\n  if(b <= c){\n    if(b <= d) return b;\n    return d;\n  }\n  if(c <= d) return c;\n  return d;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "Kth1_size4";
      string c = "template<class T1, class T2, class T3, class T4>\ninline T1 Kth1_L(const T1 a, const T2 b, const T3 c, const T4 d){\n  if(a <= b){\n    if(c <= d){\n      if(a <= c){\n        if(b <= c) return b;\n        return c;\n      }\n      if(a <= d) return a;\n      return d;\n    }\n    if(a <= d){\n      if(b <= d) return b;\n      return d;\n    }\n    if(a <= c) return a;\n    return c;\n  }\n  if(c <= d){\n    if(b <= c){\n      if(a <= c) return a;\n      return c;\n    }\n    if(b <= d) return b;\n    return d;\n  }\n  if(b <= d){\n    if(a <= d) return a;\n    return d;\n  }\n  if(b <= c) return b;\n  return c;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "Kth2_size4";
      string c = "template<class T1, class T2, class T3, class T4>\ninline T1 Kth2_L(const T1 a, const T2 b, const T3 c, const T4 d){\n  if(a >= b){\n    if(c >= d){\n      if(a >= c){\n        if(b >= c) return b;\n        return c;\n      }\n      if(a >= d) return a;\n      return d;\n    }\n    if(a >= d){\n      if(b >= d) return b;\n      return d;\n    }\n    if(a >= c) return a;\n    return c;\n  }\n  if(c >= d){\n    if(b >= c){\n      if(a >= c) return a;\n      return c;\n    }\n    if(b >= d) return b;\n    return d;\n  }\n  if(b >= d){\n    if(a >= d) return a;\n    return d;\n  }\n  if(b >= c) return b;\n  return c;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "Kth3_size4";
      string c = "template<class T1, class T2, class T3, class T4>\ninline T1 Kth3_L(const T1 a, const T2 b, const T3 c, const T4 d){\n  if(a >= b){\n    if(a >= c){\n      if(a >= d) return a;\n      return d;\n    }\n    if(c >= d) return c;\n    return d;\n  }\n  if(b >= c){\n    if(b >= d) return b;\n    return d;\n  }\n  if(c >= d) return c;\n  return d;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "KthA";
      string c = "template<class T>\nT KthA_L(const int K, const int N, const T A[], void *mem = wmem){\n  int i;\n  T *a;\n  walloc1d(&a, N, &mem);\n  rep(i,N) a[i] = A[i];\n  nth_element(a, a+K, a+N);\n  return a[K];\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "next_mcomb";
      string c = "int next_mcomb(int len, int arr[], int lim){\n  int i;\n  rrep(i,len){\n    if(arr[i]+1 < lim){\n      arr[i]++;\n      i++;\n      while(i < len){\n        arr[i] = arr[i-1];\n        i++;\n      }\n      return 1;\n    }\n    arr[i] = 0;\n  }\n  rep(i,len) arr[i] = 0;\n  return 0;\n}\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "next_scomb";
      string c = "int next_scomb(int len, int arr[], int lim){\n  int i;\n  rrep(i,len){\n    if(arr[i] < lim+i-len){\n      arr[i]++;\n      i++;\n      while(i < len){\n        arr[i] = arr[i-1] + 1;\n        i++;\n      }\n      return 1;\n    }\n    arr[i] = 0;\n  }\n  rep(i,len) arr[i] = i;\n  return 0;\n}\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "next_marr";
      string c = "int next_marr(int len, int arr[], int lim){\n  int i;\n  rrep(i,len){\n    if(arr[i]+1 < lim) arr[i]++, return 1;\n    arr[i] = 0;\n  }\n  return 0;\n}\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "next_sarr";
      string c = "int next_sarr(int len, int arr[], int lim, void *mem = wmem){\n  int i, j;\n  char *use;\n  walloc1d(&use, lim, &mem);\n  rep(i,lim) use[i] = 0;\n  rep(i,len) use[arr[i]]++;\n  rrep(i,len){\n    use[arr[i]++] = 0;\n    while(arr[i] < lim && use[arr[i]]) arr[i]++;\n    if(arr[i] == lim) continue;\n\n    use[arr[i++]] = 1;\n    j = 0;\n    while(i < len){\n      while(use[j]) j++;\n      arr[i++] = j;\n      use[j] = 1;\n    }\n    return 1;\n  }\n  rep(i,len){\n    arr[i] = i;\n    use[i] = 1;\n  }\n  return 0;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "next_sarr_s";
      string c = "template<class T>\nint next_sarr_s(int len, int arr[], int lim, T use[]){\n  int i, j;\n  rrep(i,len){\n    use[arr[i]++] = 0;\n    while(arr[i] < lim && use[arr[i]]) arr[i]++;\n    if(arr[i] == lim) continue;\n\n    use[arr[i++]] = 1;\n    j = 0;\n    while(i < len){\n      while(use[j]) j++;\n      arr[i++] = j;\n      use[j] = 1;\n    }\n    return 1;\n  }\n  rep(i,len){\n    arr[i] = i;\n    use[i] = 1;\n  }\n  return 0;\n}\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "Timer";
      string c = "struct Timer{\n  double x;\n  double gettimeofday_sec(void){\n    timeval t;\n    gettimeofday(&t, 0);\n    return t.tv_sec + t.tv_usec * 1e-6;\n  }\n  void set(void){\n    x = gettimeofday_sec();\n  }\n  double get(void){\n    return gettimeofday_sec() - x;\n  }\n};\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"sys_time");
      need[n] = d;
    }
    {
      string n = "Timer2";
      string c = "struct Timer2{\n  std::chrono::time_point<std::chrono::high_resolution_clock> start_time;\n  void set(void) {\n    start_time = std::chrono::high_resolution_clock::now();\n  }\n  double get(void) {\n    auto end_time = std::chrono::high_resolution_clock::now();\n    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);\n    return duration.count() * 1e-6;\n  }\n};\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"chrono");
      need[n] = d;
    }


    {
      string n = "Rand";
      string c = "struct Rand{\n  unsigned x,y,z,w;\n\n  Rand(void){\n    x=123456789, y=362436069, z=521288629, w=(unsigned)time(NULL);\n  }\n  Rand(unsigned seed){\n    x=123456789, y=362436069, z=521288629, w=seed;\n  }\n  inline unsigned get(void){\n    unsigned t;\n    t = (x^(x<<11));\n    x=y; y=z; z=w;\n    w = (w^(w>>19))^(t^(t>>8));\n    return w;\n  }\n  inline double getUni(void){\n    return get()/4294967296.0;\n  }\n  inline int get(int a){\n    return (int)(a*getUni());\n  }\n  inline int get(int a, int b){\n    return a+(int)((b-a+1)*getUni());\n  }\n  inline ll get(ll a){\n    return(ll)(a*getUni());\n  }\n  inline ll get(ll a, ll b){\n    return a+(ll)((b-a+1)*getUni());\n  }\n  inline double get(double a, double b){\n    return a+(b-a)*getUni();\n  }\n  inline int getExp(int a){\n    return(int)(exp(getUni()*log(a+1.0))-1.0);\n  }\n  inline int getExp(int a, int b){\n    return a+(int)(exp(getUni()*log((b-a+1)+1.0))-1.0);\n  }\n};\n";
      string p = "first";
      vector<string> d;

      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "Modint";
      string c = "struct Modint{\n  unsigned val;\n  Modint(){val=0;}\n  Modint(int a){val = ord(a);}\n  Modint(unsigned a){val = ord(a);}\n  Modint(ll a){val = ord(a);}\n  Modint(ull a){val = ord(a);}\n  inline unsigned ord(unsigned a){\n    return a%MD;\n  }\n  inline unsigned ord(int a){\n    a %= (int)MD;\n    if(a < 0) a += MD;\n    return a;\n  }\n  inline unsigned ord(ull a){\n    return a%MD;\n  }\n  inline unsigned ord(ll a){\n    a %= (int)MD;\n    if(a < 0) a += MD;\n    return a;\n  }\n  inline unsigned get(){\n    return val;\n  }\n  inline Modint &operator++(){\n    val++;\n    if(val >= MD) val -= MD;\n    return *this;\n  }\n  inline Modint &operator--(){\n    if(val == 0) val = MD - 1;\n    else         --val;\n    return *this;\n  }\n  inline Modint operator++(int a){\n    Modint res(*this);\n    val++;\n    if(val >= MD) val -= MD;\n    return res;\n  }\n  inline Modint operator--(int a){\n    Modint res(*this);\n    if(val == 0) val = MD - 1;\n    else         --val;\n    return res;\n  }\n  inline Modint &operator+=(Modint a){\n    val += a.val;\n    if(val >= MD) val -= MD;\n    return *this;\n  }\n  inline Modint &operator-=(Modint a){\n    if(val < a.val) val = val + MD - a.val;\n    else val -= a.val;\n    return *this;\n  }\n  inline Modint &operator*=(Modint a){\n    val = ((ull)val*a.val)%MD;\n    return *this;\n  }\n  inline Modint &operator/=(Modint a){\n    return *this *= a.inverse();\n  }\n  inline Modint operator+(Modint a){ return Modint(*this)+=a; }\n  inline Modint operator-(Modint a){ return Modint(*this)-=a; }\n  inline Modint operator*(Modint a){ return Modint(*this)*=a; }\n  inline Modint operator/(Modint a){ return Modint(*this)/=a; }\n  inline Modint operator+(int a){ return Modint(*this)+=Modint(a); }\n  inline Modint operator-(int a){ return Modint(*this)-=Modint(a); }\n  inline Modint operator*(int a){ return Modint(*this)*=Modint(a); }\n  inline Modint operator/(int a){ return Modint(*this)/=Modint(a); }\n  inline Modint operator+(ll a){ return Modint(*this)+=Modint(a); }\n  inline Modint operator-(ll a){ return Modint(*this)-=Modint(a); }\n  inline Modint operator*(ll a){ return Modint(*this)*=Modint(a); }\n  inline Modint operator/(ll a){ return Modint(*this)/=Modint(a); }\n  inline Modint operator-(void){ Modint res; if(val) res.val=MD-val; else res.val=0; return res; }\n  \n  inline operator bool(void){\n    return val!=0;\n  }\n  inline operator int(void){\n    return get();\n  }\n  inline operator ll(void){\n    return get();\n  }\n  inline Modint inverse(){\n    int a = val, b = MD, u = 1, v = 0, t;\n    Modint res;\n    while(b){\n      t = a / b;\n      a -= t * b; swap(a, b);\n      u -= t * v; swap(u, v);\n    }\n    if(u < 0) u += MD;\n    res.val = u;\n    return res;\n  }\n  inline Modint pw(ull b){\n    Modint a(*this), res;\n    res.val = 1;\n    while(b){\n      if(b&1) res *= a;\n      b >>= 1;\n      a *= a;\n    }\n    return res;\n  }\n  inline bool operator==(int a){return ord(a)==val;}\n  inline bool operator!=(int a){return ord(a)!=val;}\n};\ninline Modint operator+(int a, Modint b){return Modint(a)+=b;}\ninline Modint operator-(int a, Modint b){return Modint(a)-=b;}\ninline Modint operator*(int a, Modint b){return Modint(a)*=b;}\ninline Modint operator/(int a, Modint b){return Modint(a)/=b;}\ninline Modint operator+(ll a, Modint b){return Modint(a)+=b;}\ninline Modint operator-(ll a, Modint b){return Modint(a)-=b;}\ninline Modint operator*(ll a, Modint b){return Modint(a)*=b;}\ninline Modint operator/(ll a, Modint b){return Modint(a)/=b;}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"define_MD");
      need[n] = d;
    }
    {
      string n = "modint";
      string c = "struct modint{\n  static unsigned md;\n  unsigned val;\n  modint(){val=0;}\n  modint(int a){val = ord(a);}\n  modint(unsigned a){val = ord(a);}\n  modint(ll a){val = ord(a);}\n  modint(ull a){val = ord(a);}\n  void setmod(unsigned m){\n    md = m;\n  }\n  unsigned ord(unsigned a){\n    return a%md;\n  }\n  unsigned ord(int a){\n    a %= (int)md;\n    if(a < 0) a += md;\n    return a;\n  }\n  unsigned ord(ull a){\n    return a%md;\n  }\n  unsigned ord(ll a){\n    a %= (int)md;\n    if(a < 0) a += md;\n    return a;\n  }\n  unsigned get(){\n    return val;\n  }\n  inline modint &operator++(){\n    val++;\n    if(val >= md) val -= md;\n    return *this;\n  }\n  inline modint &operator--(){\n    if(val == 0) val = md - 1;\n    else         --val;\n    return *this;\n  }\n  inline modint operator++(int a){\n    modint res(*this);\n    val++;\n    if(val >= md) val -= md;\n    return res;\n  }\n  inline modint operator--(int a){\n    modint res(*this);\n    if(val == 0) val = md - 1;\n    else         --val;\n    return res;\n  }\n  modint &operator+=(modint a){\n    val += a.val;\n    if(val >= md) val -= md;\n    return *this;\n  }\n  modint &operator-=(modint a){\n    if(val < a.val) val = val + md - a.val;\n    else val -= a.val;\n    return *this;\n  }\n  modint &operator*=(modint a){\n    val = ((ull)val*a.val)%md;\n    return *this;\n  }\n  modint &operator/=(modint a){\n    return *this *= a.inverse();\n  }\n  modint operator+(modint a){ return modint(*this)+=a; }\n  modint operator-(modint a){ return modint(*this)-=a; }\n  modint operator*(modint a){ return modint(*this)*=a; }\n  modint operator/(modint a){ return modint(*this)/=a; }\n  modint operator+(int a){ return modint(*this)+=modint(a); }\n  modint operator-(int a){ return modint(*this)-=modint(a); }\n  modint operator*(int a){ return modint(*this)*=modint(a); }\n  modint operator/(int a){ return modint(*this)/=modint(a); }\n  modint operator+(ll a){ return modint(*this)+=modint(a); }\n  modint operator-(ll a){ return modint(*this)-=modint(a); }\n  modint operator*(ll a){ return modint(*this)*=modint(a); }\n  modint operator/(ll a){ return modint(*this)/=modint(a); }\n  modint operator-(void){ modint res; if(val) res.val=md-val; else res.val=0; return res; }\n  \n  operator bool(void){\n    return val!=0;\n  }\n  operator int(void){\n    return get();\n  }\n  operator ll(void){\n    return get();\n  }\n  modint inverse(){\n    int a = val, b = md, u = 1, v = 0, t;\n    modint res;\n    while(b){\n      t = a / b;\n      a -= t * b; swap(a, b);\n      u -= t * v; swap(u, v);\n    }\n    if(u < 0) u += md;\n    res.val = u;\n    return res;\n  }\n  modint pw(ull b){\n    modint a(*this), res;\n    res.val = 1;\n    while(b){\n      if(b&1) res *= a;\n      b >>= 1;\n      a *= a;\n    }\n    return res;\n  }\n  bool operator==(int a){return ord(a)==val;}\n  bool operator!=(int a){return ord(a)!=val;}\n};\nunsigned modint::md;\nmodint operator+(int a, modint b){return modint(a)+=b;}\nmodint operator-(int a, modint b){return modint(a)-=b;}\nmodint operator*(int a, modint b){return modint(a)*=b;}\nmodint operator/(int a, modint b){return modint(a)/=b;}\nmodint operator+(ll a, modint b){return modint(a)+=b;}\nmodint operator-(ll a, modint b){return modint(a)-=b;}\nmodint operator*(ll a, modint b){return modint(a)*=b;}\nmodint operator/(ll a, modint b){return modint(a)/=b;}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"modint_init");
      d.push_back((string)"define_MD");
      need[n] = d;
    }
    {
      string n = "modint_init";
      string c = "{modint x; x.setmod(MD);}\n";
      string p = "main_first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "Mint";
      string c = "struct Mint{\n  unsigned val;\n  Mint(){val=0;}\n  Mint(int a){val = mulR(a);}\n  Mint(unsigned a){val = mulR(a);}\n  Mint(ll a){val = mulR(a);}\n  Mint(ull a){val = mulR(a);}\n  inline unsigned mulR(unsigned a){\n    return (ull)a*MINT_R%MD;\n  }\n  inline unsigned mulR(int a){\n    if(a < 0) a = a%((int)MD)+(int)MD;\n    return mulR((unsigned)a);\n  }\n  inline unsigned mulR(ull a){\n    return mulR((unsigned)(a%MD));\n  }\n  inline unsigned mulR(ll a){\n    a %= (int)MD;\n    if(a < 0) a += MD;\n    return mulR((unsigned)a);\n  }\n  inline unsigned reduce(unsigned T){\n    unsigned m = T * MINT_MDNINV;\n    unsigned t = (unsigned)((T + (ull)m*MD) >> MINT_W);\n    if(t >= MD) t -= MD;\n    return t;\n  }\n  inline unsigned reduce(ull T){\n    unsigned m = (unsigned)T * MINT_MDNINV;\n    unsigned t = (unsigned)((T + (ull)m*MD) >> MINT_W);\n    if(t >= MD) t -= MD;\n    return t;\n  }\n  inline unsigned get(){\n    return reduce(val);\n  }\n  inline Mint &operator++(){\n    (*this) += 1;\n    return *this;\n  }\n  inline Mint &operator--(){\n    (*this) -= 1;\n    return *this;\n  }\n  inline Mint operator++(int a){\n    Mint res(*this);\n    (*this) += 1;\n    return res;\n  }\n  inline Mint operator--(int a){\n    Mint res(*this);\n    (*this) -= 1;\n    return res;\n  }\n  inline Mint &operator+=(Mint a){\n    val += a.val;\n    if(val >= MD) val -= MD;\n    return *this;\n  }\n  inline Mint &operator-=(Mint a){\n    if(val < a.val) val = val + MD - a.val;\n    else val -= a.val;\n    return *this;\n  }\n  inline Mint &operator*=(Mint a){\n    val = reduce((ull)val*a.val);\n    return *this;\n  }\n  inline Mint &operator/=(Mint a){\n    return *this *= a.inverse();\n  }\n  inline Mint operator+(Mint a){ return Mint(*this)+=a; }\n  inline Mint operator-(Mint a){ return Mint(*this)-=a; }\n  inline Mint operator*(Mint a){ return Mint(*this)*=a; }\n  inline Mint operator/(Mint a){ return Mint(*this)/=a; }\n  inline Mint operator+(int a){ return Mint(*this)+=Mint(a); }\n  inline Mint operator-(int a){ return Mint(*this)-=Mint(a); }\n  inline Mint operator*(int a){ return Mint(*this)*=Mint(a); }\n  inline Mint operator/(int a){ return Mint(*this)/=Mint(a); }\n  inline Mint operator+(ll a){ return Mint(*this)+=Mint(a); }\n  inline Mint operator-(ll a){ return Mint(*this)-=Mint(a); }\n  inline Mint operator*(ll a){ return Mint(*this)*=Mint(a); }\n  inline Mint operator/(ll a){ return Mint(*this)/=Mint(a); }\n  inline Mint operator-(void){ Mint res; if(val) res.val=MD-val; else res.val=0; return res; }\n  \n  inline operator bool(void){\n    return val!=0;\n  }\n  inline operator int(void){\n    return get();\n  }\n  inline operator ll(void){\n    return get();\n  }\n  inline Mint inverse(){\n    int a = val, b = MD, u = 1, v = 0, t;\n    Mint res;\n    while(b){\n      t = a / b;\n      a -= t * b; swap(a, b);\n      u -= t * v; swap(u, v);\n    }\n    if(u < 0) u += MD;\n    res.val = (ull)u*MINT_RR % MD;\n    return res;\n  }\n  inline Mint pw(ull b){\n    Mint a(*this), res;\n    res.val = MINT_R;\n    while(b){\n      if(b&1) res *= a;\n      b >>= 1;\n      a *= a;\n    }\n    return res;\n  }\n  inline bool operator==(int a){return mulR(a)==val;}\n  inline bool operator!=(int a){return mulR(a)!=val;}\n};\ninline Mint operator+(int a, Mint b){return Mint(a)+=b;}\ninline Mint operator-(int a, Mint b){return Mint(a)-=b;}\ninline Mint operator*(int a, Mint b){return Mint(a)*=b;}\ninline Mint operator/(int a, Mint b){return Mint(a)/=b;}\ninline Mint operator+(ll a, Mint b){return Mint(a)+=b;}\ninline Mint operator-(ll a, Mint b){return Mint(a)-=b;}\ninline Mint operator*(ll a, Mint b){return Mint(a)*=b;}\ninline Mint operator/(ll a, Mint b){return Mint(a)/=b;}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"define_MD");
      d.push_back((string)"define_for_Mint");
      need[n] = d;
    }
    {
      string n = "mint";
      string c = "struct mint{\n  static unsigned md, W, R, Rinv, mdninv, RR;\n  unsigned val;\n  mint(){val=0;}\n  mint(int a){val = mulR(a);}\n  mint(unsigned a){val = mulR(a);}\n  mint(ll a){val = mulR(a);}\n  mint(ull a){val = mulR(a);}\n  int get_inv(ll a, int md){ll t=a,s=md,u=1,v=0,e;while(s){e=t/s;t-=e*s;u-=e*v;swap(t,s);swap(u,v);}if(u<0)u+=md;return u;}\n  void setmod(unsigned m){\n    int i;\n    unsigned t;\n    W = 32;\n    md = m;\n    R = (1ULL << W) % md;\n    RR = (ull)R*R % md;\n    switch(m){\n    case 104857601:\n      Rinv = 2560000;\n      mdninv = 104857599;\n      break;\n    case 998244353:\n      Rinv = 232013824;\n      mdninv = 998244351;\n      break;\n    case 1000000007:\n      Rinv = 518424770;\n      mdninv = 2226617417U;\n      break;\n    case 1000000009:\n      Rinv = 171601999;\n      mdninv = 737024967;\n      break;\n    case 1004535809:\n      Rinv = 234947584;\n      mdninv = 1004535807;\n      break;\n    case 1007681537:\n      Rinv = 236421376;\n      mdninv = 1007681535;\n      break;\n    case 1012924417:\n      Rinv = 238887936;\n      mdninv = 1012924415;\n      break;\n    case 1045430273:\n      Rinv = 254466304;\n      mdninv = 1045430271;\n      break;\n    case 1051721729:\n      Rinv = 257538304;\n      mdninv = 1051721727;\n      break;\n    default:\n      Rinv = get_inv(R, md);\n      mdninv = 0;\n      t = 0;\n      rep(i,(int)W){\n        if(t%2==0) t+=md, mdninv |= (1U<<i);\n        t /= 2;\n      }\n    }\n  }\n  unsigned mulR(unsigned a){\n    return (ull)a*R%md;\n  }\n  unsigned mulR(int a){\n    if(a < 0) a = a%((int)md)+(int)md;\n    return mulR((unsigned)a);\n  }\n  unsigned mulR(ull a){\n    return mulR((unsigned)(a%md));\n  }\n  unsigned mulR(ll a){\n    a %= (int)md;\n    if(a < 0) a += md;\n    return mulR((unsigned)a);\n  }\n  unsigned reduce(unsigned T){\n    unsigned m = T * mdninv;\n    unsigned t = (unsigned)((T + (ull)m*md) >> W);\n    if(t >= md) t -= md;\n    return t;\n  }\n  unsigned reduce(ull T){\n    unsigned m = (unsigned)T * mdninv;\n    unsigned t = (unsigned)((T + (ull)m*md) >> W);\n    if(t >= md) t -= md;\n    return t;\n  }\n  unsigned get(){\n    return reduce(val);\n  }\n  inline mint &operator++(){\n    (*this) += 1;\n    return *this;\n  }\n  inline mint &operator--(){\n    (*this) -= 1;\n    return *this;\n  }\n  inline mint operator++(int a){\n    mint res(*this);\n    (*this) += 1;\n    return res;\n  }\n  inline mint operator--(int a){\n    mint res(*this);\n    (*this) -= 1;\n    return res;\n  }\n  mint &operator+=(mint a){\n    val += a.val;\n    if(val >= md) val -= md;\n    return *this;\n  }\n  mint &operator-=(mint a){\n    if(val < a.val) val = val + md - a.val;\n    else val -= a.val;\n    return *this;\n  }\n  mint &operator*=(mint a){\n    val = reduce((ull)val*a.val);\n    return *this;\n  }\n  mint &operator/=(mint a){\n    return *this *= a.inverse();\n  }\n  mint operator+(mint a){ return mint(*this)+=a; }\n  mint operator-(mint a){ return mint(*this)-=a; }\n  mint operator*(mint a){ return mint(*this)*=a; }\n  mint operator/(mint a){ return mint(*this)/=a; }\n  mint operator+(int a){ return mint(*this)+=mint(a); }\n  mint operator-(int a){ return mint(*this)-=mint(a); }\n  mint operator*(int a){ return mint(*this)*=mint(a); }\n  mint operator/(int a){ return mint(*this)/=mint(a); }\n  mint operator+(ll a){ return mint(*this)+=mint(a); }\n  mint operator-(ll a){ return mint(*this)-=mint(a); }\n  mint operator*(ll a){ return mint(*this)*=mint(a); }\n  mint operator/(ll a){ return mint(*this)/=mint(a); }\n  mint operator-(void){ mint res; if(val) res.val=md-val; else res.val=0; return res; }\n  \n  operator bool(void){\n    return val!=0;\n  }\n  operator int(void){\n    return get();\n  }\n  operator ll(void){\n    return get();\n  }\n  mint inverse(){\n    int a = val, b = md, u = 1, v = 0, t;\n    mint res;\n    while(b){\n      t = a / b;\n      a -= t * b; swap(a, b);\n      u -= t * v; swap(u, v);\n    }\n    if(u < 0) u += md;\n    res.val = (ull)u*RR % md;\n    return res;\n  }\n  mint pw(ull b){\n    mint a(*this), res;\n    res.val = R;\n    while(b){\n      if(b&1) res *= a;\n      b >>= 1;\n      a *= a;\n    }\n    return res;\n  }\n  bool operator==(int a){return mulR(a)==val;}\n  bool operator!=(int a){return mulR(a)!=val;}\n};\nunsigned mint::md, mint::W, mint::R, mint::Rinv, mint::mdninv, mint::RR;\nmint operator+(int a, mint b){return mint(a)+=b;}\nmint operator-(int a, mint b){return mint(a)-=b;}\nmint operator*(int a, mint b){return mint(a)*=b;}\nmint operator/(int a, mint b){return mint(a)/=b;}\nmint operator+(ll a, mint b){return mint(a)+=b;}\nmint operator-(ll a, mint b){return mint(a)-=b;}\nmint operator*(ll a, mint b){return mint(a)*=b;}\nmint operator/(ll a, mint b){return mint(a)/=b;}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"mint_init");
      d.push_back((string)"define_MD");
      need[n] = d;
    }
    {
      string n = "mint_init";
      string c = "{mint x; x.setmod(MD);}\n";
      string p = "main_first";
      name.push_back(n); func[n] = c; place[n] = p;
    }


    {
      string n = "readerFile";
      string c = "inplace_L FILE *readerfp = stdin;\ninplace_L int readermode = 0;\nvoid readerFile(){\n  if(readermode) fclose(readerfp);\n  readerfp = stdin;\n  readermode = 0;\n}\nvoid readerFile(string filename, string mode = \"r\"){\n  if(readermode) fclose(readerfp);\n  readerfp = fopen(filename.c_str(), mode.c_str());\n  readermode = 1;\n}\nvoid readerFile(FILE *fp){\n  if(readermode) fclose(readerfp);\n  readerfp = fp;\n  readermode = 0;\n}\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "my_getchar_unlocked";
      string c = "inline int my_getchar_unlocked(){\n  static char buf[1048576];\n  static int s = 1048576, e = 1048576;\n  if(s == e && e == 1048576){\n    e = fread_unlocked(buf, 1, 1048576, stdin);\n    s = 0;\n  }\n  if(s == e) return EOF;\n  return buf[s++];\n}\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "reader_all";
      string c = "";
      string p = "first";
      vector<string> d;
      d.push_back((string)"reader_int");
      d.push_back((string)"reader_ll");
      d.push_back((string)"reader_unsigned");
      d.push_back((string)"reader_ull");
      d.push_back((string)"reader_int128");
      d.push_back((string)"reader_uint128");
      d.push_back((string)"reader_Modint");
      d.push_back((string)"reader_Mint");
      d.push_back((string)"reader_modint");
      d.push_back((string)"reader_mint");
      d.push_back((string)"reader_double");
      d.push_back((string)"reader_char");
      d.push_back((string)"reader_char_array");
      d.push_back((string)"reader_string");
      d.push_back((string)"reader_Point2d");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "reader_int";
      string c = "inline void rd(int &x){\n  int k, m=0;\n  x=0;\n  for(;;){\n    k = my_getchar_unlocked();\n    if(k=='-'){\n      m=1;\n      break;\n    }\n    if('0'<=k&&k<='9'){\n      x=k-'0';\n      break;\n    }\n  }\n  for(;;){\n    k = my_getchar_unlocked();\n    if(k<'0'||k>'9'){\n      break;\n    }\n    x=x*10+k-'0';\n  }\n  if(m){\n    x=-x;\n  }\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"my_getchar_unlocked");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "reader_unsigned";
      string c = "inline void rd(unsigned &x){\n  int k;\n  x=0;\n  for(;;){\n    k = my_getchar_unlocked();\n    if('0'<=k&&k<='9'){\n      x=k-'0';\n      break;\n    }\n  }\n  for(;;){\n    k = my_getchar_unlocked();\n    if(k<'0'||k>'9'){\n      break;\n    }\n    x=x*10+k-'0';\n  }\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"my_getchar_unlocked");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "reader_ll";
      string c = "inline void rd(ll &x){\n  int k, m=0;\n  x=0;\n  for(;;){\n    k = my_getchar_unlocked();\n    if(k=='-'){\n      m=1;\n      break;\n    }\n    if('0'<=k&&k<='9'){\n      x=k-'0';\n      break;\n    }\n  }\n  for(;;){\n    k = my_getchar_unlocked();\n    if(k<'0'||k>'9'){\n      break;\n    }\n    x=x*10+k-'0';\n  }\n  if(m){\n    x=-x;\n  }\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"my_getchar_unlocked");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "reader_ull";
      string c = "inline void rd(ull &x){\n  int k;\n  x=0;\n  for(;;){\n    k = my_getchar_unlocked();\n    if('0'<=k&&k<='9'){\n      x=k-'0';\n      break;\n    }\n  }\n  for(;;){\n    k = my_getchar_unlocked();\n    if(k<'0'||k>'9'){\n      break;\n    }\n    x=x*10+k-'0';\n  }\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"my_getchar_unlocked");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "reader_int128";
      string c = "inline void rd(__int128_t &x){\n  int k, m=0;\n  x=0;\n  for(;;){\n    k = my_getchar_unlocked();\n    if(k=='-'){\n      m=1;\n      break;\n    }\n    if('0'<=k&&k<='9'){\n      x=k-'0';\n      break;\n    }\n  }\n  for(;;){\n    k = my_getchar_unlocked();\n    if(k<'0'||k>'9'){\n      break;\n    }\n    x=x*10+k-'0';\n  }\n  if(m){\n    x=-x;\n  }\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"my_getchar_unlocked");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p; parent[n] = "__int128_t";
    }
    {
      string n = "reader_uint128";
      string c = "inline void rd(__uint128_t &x){\n  int k;\n  x=0;\n  for(;;){\n    k = my_getchar_unlocked();\n    if('0'<=k&&k<='9'){\n      x=k-'0';\n      break;\n    }\n  }\n  for(;;){\n    k = my_getchar_unlocked();\n    if(k<'0'||k>'9'){\n      break;\n    }\n    x=x*10+k-'0';\n  }\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"my_getchar_unlocked");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p; parent[n] = "__uint128_t";
    }
    {
      string n = "reader_Modint";
      string c = "inline void rd(Modint &x){int i; rd(i); x=i;}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"reader_int");
      d.push_back((string)"my_getchar_unlocked");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p; parent[n] = "Modint";
    }
    {
      string n = "reader_modint";
      string c = "inline void rd(modint &x){int i; rd(i); x=i;}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"reader_int");
      d.push_back((string)"my_getchar_unlocked");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p; parent[n] = "modint";
    }
    
    {
      string n = "reader_Mint";
      string c = "inline void rd(Mint &x){int i; rd(i); x=i;}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"reader_int");
      d.push_back((string)"my_getchar_unlocked");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p; parent[n] = "Mint";
    }
    
    {
      string n = "reader_mint";
      string c = "inline void rd(mint &x){int i; rd(i); x=i;}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"reader_int");
      d.push_back((string)"my_getchar_unlocked");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p; parent[n] = "mint";
    }
    
    {
      string n = "reader_double";
      string c = "inline void rd(double &x){\n  int k, m=0, p=0;\n  double r = 1;\n  x = 0;\n  for(;;){\n    k = my_getchar_unlocked();\n    if(k=='-') m = 1, break;\n    if(k=='.') p = 1, break;\n    if('0'<=k&&k<='9') x = k - '0', break;\n  }\n  for(;;){\n    k = my_getchar_unlocked();\n    if(k=='.') p = 1, continue;\n    if(k<'0'||k>'9') break;\n    if(p){\n      r *= 0.1;\n      x += r * (k - '0');\n    } else {\n      x = x * 10 + k - '0';\n    }\n  }\n  if(m) x = -x;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"my_getchar_unlocked");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "reader_char";
      string c = "inline void rd(char &c){\n  int i;\n    for(;;){\n    i = my_getchar_unlocked();\n    if(i!=' '&&i!='\\n'&&i!='\\r'&&i!='\\t'&&i!=EOF) break;\n  }\n  c = i;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"my_getchar_unlocked");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "reader_char_array";
      string c = "inline int rd(char c[]){\n  int i, sz = 0;\n  for(;;){\n    i = my_getchar_unlocked();\n    if(i!=' '&&i!='\\n'&&i!='\\r'&&i!='\\t'&&i!=EOF) break;\n  }\n  c[sz++] = i;\n  for(;;){\n    i = my_getchar_unlocked();\n    if(i==' '||i=='\\n'||i=='\\r'||i=='\\t'||i==EOF) break;\n    c[sz++] = i;\n  }\n  c[sz]='\\0';\n  return sz;\n}";
      string p = "first";
      vector<string> d;
      d.push_back((string)"my_getchar_unlocked");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "reader_string";
      string c = "inline void rd(string &x){\n  char *buf = (char *)wmem;\n  rd(buf);\n  x = buf;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"reader_char_array");
      d.push_back((string)"workmemory");
      d.push_back((string)"my_getchar_unlocked");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "rdLine";
      string c = "inline int rdLine_L(char c[]){\n  int i, sz = 0;\n  for(;;){\n    i = my_getchar_unlocked();\n    if(i=='\\r') continue;\n    if(i=='\\n') break;\n    if(i==EOF){\n      if(sz==0){\n        c[sz] = '\\0';\n        return -1;\n      }\n      break;\n    }\n    c[sz++] = i;\n  }\n  c[sz]='\\0';\n  return sz;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"my_getchar_unlocked");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "rd_int";
      string c = "inline int rd_int(void){int x; rd(x); return x;}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"reader_int");
      d.push_back((string)"my_getchar_unlocked");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "rd_ll";
      string c = "inline ll rd_ll(void){ll x; rd(x); return x;}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"reader_ll");
      d.push_back((string)"my_getchar_unlocked");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "rd_string";
      string c = "inline string rd_string(void){string x; rd(x); return x;}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"reader_string");
      d.push_back((string)"my_getchar_unlocked");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "writerFile";
      string c = "inplace_L FILE *writerfp = stdout;\ninplace_L int writermode = 0;\nvoid writerFile(){\n  if(writermode) fclose(writerfp);\n  writerfp = stdout;\n  writermode = 0;\n}\nvoid writerFile(string filename, string mode = \"w\"){\n  if(writermode) fclose(writerfp);\n  writerfp = fopen(filename.c_str(), mode.c_str());\n  writermode = 1;\n}\nvoid writerFile(FILE *fp){\n  if(writermode) fclose(writerfp);\n  writerfp = fp;\n  writermode = 0;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "my_putchar_unlocked";
      string c = "struct MY_WRITER {\n  char buf[1048576];\n  int s, e;\n  MY_WRITER(){\n    s = 0;\n    e = 1048576;\n  }\n  ~MY_WRITER(){\n    if(s) fwrite_unlocked(buf, 1, s, stdout);\n  }\n};\n\ninplace_L MY_WRITER MY_WRITER_VAR;\n\nvoid my_putchar_unlocked(int a){\n  if(MY_WRITER_VAR.s == MY_WRITER_VAR.e){\n    fwrite_unlocked(MY_WRITER_VAR.buf, 1, MY_WRITER_VAR.s, stdout);\n    MY_WRITER_VAR.s = 0;\n  }\n  MY_WRITER_VAR.buf[MY_WRITER_VAR.s++] = a;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    
    {
      string n = "writer_vector_head";
      string c = "template<class T>\ninline void wt_L(const vector<T> &x);\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "writer_vector";
      string c = "template<class T>\ninline void wt_L(const vector<T> &x){\n  int fg = 0;\n  for(auto a : x){\n    if(fg) my_putchar_unlocked(' ');\n    fg = 1;\n    wt_L(a);\n  }\n}\n";
      string p = "last";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"writer_vector_head");
      need[n] = d;
    }
    {
      string n = "writer_set_head";
      string c = "template<class T>\ninline void wt_L(const set<T> &x);\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "writer_set";
      string c = "template<class T>\ninline void wt_L(const set<T> &x){\n  int fg = 0;\n  for(auto a : x){\n    if(fg) my_putchar_unlocked(' ');\n    fg = 1;\n    wt_L(a);\n  }\n}\n";
      string p = "last";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"writer_set_head");
      need[n] = d;
    }
    {
      string n = "writer_multiset_head";
      string c = "template<class T>\ninline void wt_L(const multiset<T> &x);\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "writer_multiset";
      string c = "template<class T>\ninline void wt_L(const multiset<T> &x){\n  int fg = 0;\n  for(auto a : x){\n    if(fg) my_putchar_unlocked(' ');\n    fg = 1;\n    wt_L(a);\n  }\n}\n";
      string p = "last";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"writer_multiset_head");
      need[n] = d;
    }
    {
      string n = "writer_pair_head";
      string c = "template<class T1, class T2>\ninline void wt_L(const pair<T1,T2> x);\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "writer_pair";
      string c = "template<class T1, class T2>\ninline void wt_L(const pair<T1,T2> x){\n  wt_L(x.first);\n  my_putchar_unlocked(' ');\n  wt_L(x.second);\n}\n";
      string p = "last";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"writer_pair_head");
      need[n] = d;
    }


    {
      string n = "writer_all";
      string c = "";
      string p = "first";
      vector<string> d;
      d.push_back((string)"writer_char");
      d.push_back((string)"writer_int");
      d.push_back((string)"writer_unsigned");
      d.push_back((string)"writer_ll");
      d.push_back((string)"writer_ull");
      d.push_back((string)"writer_int128");
      d.push_back((string)"writer_uint128");
      d.push_back((string)"writer_Modint");
      d.push_back((string)"writer_Mint");
      d.push_back((string)"writer_modint");
      d.push_back((string)"writer_mint");
      d.push_back((string)"writer_double");
      d.push_back((string)"writer_char_array");
      d.push_back((string)"writer_string");
      d.push_back((string)"writer_Point2d");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "writer_char";
      string c = "inline void wt_L(const char a){\n  my_putchar_unlocked(a);\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"my_putchar_unlocked");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "writer_int";
      string c = "inline void wt_L(int x){\n  int s=0, m=0;\n  char f[10];\n  if(x<0) m=1, x=-x;\n  while(x) f[s++]=x%10, x/=10;\n  if(!s) f[s++]=0;\n  if(m) my_putchar_unlocked('-');\n  while(s--) my_putchar_unlocked(f[s]+'0');\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"my_putchar_unlocked");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "writer_int_withBase";
      string c = "inline void wt_L(int x, int b){\n  int s=0, m=0;\n  char f[35];\n  if(x<0) m=1, x=-x;\n  while(x) f[s++]=x%b, x/=b;\n  if(!s) f[s++]=0;\n  if(m) my_putchar_unlocked('-');\n  while(s--) my_putchar_unlocked(\"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz\"[f[s]]);\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"my_putchar_unlocked");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "writer_unsigned";
      string c = "inline void wt_L(unsigned x){\n  int s=0;\n  char f[10];\n  while(x) f[s++]=x%10, x/=10;\n  if(!s) f[s++]=0;\n  while(s--) my_putchar_unlocked(f[s]+'0');\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"my_putchar_unlocked");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "writer_ll";
      string c = "inline void wt_L(ll x){\n  int s=0, m=0;\n  char f[20];\n  if(x<0) m=1, x=-x;\n  while(x) f[s++]=x%10, x/=10;\n  if(!s) f[s++]=0;\n  if(m) my_putchar_unlocked('-');\n  while(s--) my_putchar_unlocked(f[s]+'0');\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"my_putchar_unlocked");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "writer_ll_withBase";
      string c = "inline void wt_L(ll x, int b){\n  int s=0, m=0;\n  char f[70];\n  if(x<0) m=1, x=-x;\n  while(x) f[s++]=x%b, x/=b;\n  if(!s) f[s++]=0;\n  if(m) my_putchar_unlocked('-');\n  while(s--) my_putchar_unlocked(\"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz\"[f[s]]);\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"my_putchar_unlocked");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "writer_ull";
      string c = "inline void wt_L(ull x){\n  int s=0;\n  char f[21];\n  while(x) f[s++]=x%10, x/=10;\n  if(!s) f[s++]=0;\n  while(s--) my_putchar_unlocked(f[s]+'0');\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"my_putchar_unlocked");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "writer_int128";
      string c = "inline void wt_L(__int128_t x){\n  int s=0, m=0;\n  char f[40];\n  if(x<0) m=1, x=-x;\n  while(x) f[s++]=x%10, x/=10;\n  if(!s) f[s++]=0;\n  if(m) my_putchar_unlocked('-');\n  while(s--) my_putchar_unlocked(f[s]+'0');\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"my_putchar_unlocked");
      d.push_back((string)"writer_int");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p; parent[n] = "__int128_t";
    }
    
    {
      string n = "writer_uint128";
      string c = "inline void wt_L(__uint128_t x){\n  int s=0;\n  char f[40];\n  while(x) f[s++]=x%10, x/=10;\n  if(!s) f[s++]=0;\n  while(s--) my_putchar_unlocked(f[s]+'0');\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"my_putchar_unlocked");
      d.push_back((string)"writer_int");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p; parent[n] = "__uint128_t";
    }
    
    {
      string n = "writer_Modint";
      string c = "inline void wt_L(Modint x){int i; i = (int)x; wt_L(i);}";
      string p = "first";
      vector<string> d;
      d.push_back((string)"my_putchar_unlocked");
      d.push_back((string)"writer_int");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p; parent[n] = "Modint";
    }
    
    {
      string n = "writer_modint";
      string c = "inline void wt_L(modint x){int i; i = (int)x; wt_L(i);}";
      string p = "first";
      vector<string> d;
      d.push_back((string)"my_putchar_unlocked");
      d.push_back((string)"writer_int");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p; parent[n] = "modint";
    }
    
    {
      string n = "writer_Mint";
      string c = "inline void wt_L(Mint x){int i; i = (int)x; wt_L(i);}";
      string p = "first";
      vector<string> d;
      d.push_back((string)"my_putchar_unlocked");
      d.push_back((string)"writer_int");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p; parent[n] = "Mint";
    }
    
    {
      string n = "writer_mint";
      string c = "inline void wt_L(mint x){int i; i = (int)x; wt_L(i);}";
      string p = "first";
      vector<string> d;
      d.push_back((string)"my_putchar_unlocked");
      d.push_back((string)"writer_int");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p; parent[n] = "mint";
    }
    
    {
      string n = "writer_double";
      string c = "inplace_L int WRITER_DOUBLE_DIGIT = 15;\n\ninline int writerDigit_double(){\n  return WRITER_DOUBLE_DIGIT;\n}\n\ninline void writerDigit_double(int d){\n  WRITER_DOUBLE_DIGIT = d;\n}\n\ninline void wt_L(double x){\n  const int d = WRITER_DOUBLE_DIGIT;\n  int k, r;\n  double v;\n  if(x!=x || (x==x+1 && x==2*x)){\n    my_putchar_unlocked('E');\n    my_putchar_unlocked('r');\n    my_putchar_unlocked('r');\n    return;\n  }\n  if(x < 0){\n    my_putchar_unlocked('-');\n    x = -x;\n  }\n  x += 0.5 * pow(0.1, d);\n  \n  r = 0;\n  v = 1;\n  while(x >= 10*v) v *= 10, r++;\n  while(r >= 0){\n    r--;\n    k = floor(x / v);\n    if(k >= 10) k = 9;\n    if(k <= -1) k = 0;\n    x -= k * v;\n    v *= 0.1;\n    my_putchar_unlocked(k + '0');\n  }\n\n  if(d > 0){\n    my_putchar_unlocked('.');\n    v = 1;\n    rep(r,d){\n      v *= 0.1;\n      k = floor(x / v);\n      if(k >= 10) k = 9;\n      if(k <= -1) k = 0;\n      x -= k * v;\n      my_putchar_unlocked(k + '0');\n    }\n  }\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"my_putchar_unlocked");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "writer_char_array";
      string c = "inline void wt_L(const char c[]){\n  int i=0;\n  for(i=0;c[i]!='\\0';i++) my_putchar_unlocked(c[i]);\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"my_putchar_unlocked");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "writer_string";
      string c = "inline void wt_L(const string &x){\n  int i=0;\n  for(i=0;x[i]!='\\0';i++){\n    my_putchar_unlocked(x[i]);\n  }\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"my_putchar_unlocked");
      d.push_back((string)"writer_char_array");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "Matrix";
      string c = "template<class T>\nstruct Matrix {\n  int r, c, mem;\n  T *dat;\n\n  Matrix(){r=c=mem = 0;}\n  Matrix(const int rr, const int cc){\n    if(rr == 0 || cc == 0){\n      r = c = 0;\n    } else {\n      r = rr;\n      c = cc;\n    }\n    mem = r * c;\n    if(mem > 0) dat = new T[mem];\n  }\n  Matrix(const Matrix<T> &a){\n    int i;\n    r = a.r;\n    c = a.c;\n    mem = r * c;\n    dat = new T[mem];\n    rep(i,mem) dat[i] = a.dat[i];\n  }\n  \n  ~Matrix(){\n    if(mem) delete [] dat;\n  }\n\n  void changeSize(const int rr, const int cc){\n    if(rr==0 || cc==0){\n      r = c = 0;\n    } else {\n      r = rr;\n      c = cc;\n    }\n    if(mem < r*c){\n      if(mem) delete [] dat;\n      mem = r*c;\n      dat = new T[mem];\n    }\n  }\n\n  Matrix<T>& operator=(const Matrix<T> &a){\n    int i, j;\n    r = a.r;\n    c = a.c;\n    j = r * c;\n    changeSize(r,c);\n    rep(i,j) dat[i] = a.dat[i];\n    return *this;\n  }\n\n  Matrix<T>& operator=(const int a){\n    int i, j;\n    j = r * c;\n    rep(i,j) dat[i] = 0;\n    j = min(r,c);\n    rep(i,j) dat[i*c+i] = a;\n    return *this;\n  }\n\n  Matrix<T>& operator+=(const Matrix<T> &a){\n    int i, j;\n    if(r==0 || r!=a.r || c!=a.c){\n      changeSize(0,0);\n      return *this;\n    }\n    j = r*c;\n    rep(i,j) dat[i] += a.dat[i];\n    return *this;\n  }\n  Matrix<T> operator+(const Matrix<T> &a){\n    return Matrix<T>(*this) += a;\n  }\n\n  Matrix<T>& operator-=(const Matrix<T> &a){\n    int i, j;\n    if(r==0 || r!=a.r || c!=a.c){\n      changeSize(0,0);\n      return *this;\n    }\n    j = r*c;\n    rep(i,j) dat[i] -= a.dat[i];\n    return *this;\n  }\n  Matrix<T> operator-(const Matrix<T> &a){\n    return Matrix<T>(*this) -= a;\n  }\n\n  Matrix<T>& operator*=(const Matrix<T> &a){\n    int i, j, k, x;\n    T *m;\n    if(r==0 || c!=a.r){\n      changeSize(0,0);\n      return *this;\n    }\n    m = (T*)wmem;\n    x = r * a.c;\n    rep(i,x) m[i] = 0;\n    rep(i,r) rep(k,c) rep(j,a.c) m[i*a.c+j] += dat[i*c+k] * a.dat[k*a.c+j];\n    changeSize(r, a.c);\n    rep(i,x) dat[i] = m[i];\n    return *this;\n  }\n  Matrix<T> operator*(const Matrix<T> &a){\n    return Matrix<T>(*this) *= a;\n  }\n\n  Matrix<T>& operator*=(const int a){\n    int i, j;\n    j = r * c;\n    rep(i,j) dat[i] *= a;\n    return *this;\n  }\n  Matrix<T>& operator*=(const ll a){\n    int i, j;\n    j = r * c;\n    rep(i,j) dat[i] *= a;\n    return *this;\n  }\n  Matrix<T>& operator*=(const double a){\n    int i, j;\n    j = r * c;\n    rep(i,j) dat[i] *= a;\n    return *this;\n  }\n\n\n  inline T* operator[](const int a){\n    return dat+a*c;\n  }\n};\ntemplate<class T> Matrix<T> operator*(const int a, const Matrix<T> &b){return Matrix<T>(b)*=a;}\ntemplate<class T> Matrix<T> operator*(const Matrix<T> &b, const int a){return Matrix<T>(b)*=a;}\ntemplate<class T> Matrix<T> operator*(const ll a, const Matrix<T> &b){return Matrix<T>(b)*=a;}\ntemplate<class T> Matrix<T> operator*(const Matrix<T> &b, const ll a){return Matrix<T>(b)*=a;}\ntemplate<class T> Matrix<T> operator*(const double a, const Matrix<T> &b){return Matrix<T>(b)*=a;}\ntemplate<class T> Matrix<T> operator*(const Matrix<T> &b, const double a){return Matrix<T>(b)*=a;}\n\n\ntemplate<class T, class S> inline Matrix<T> pow_L(Matrix<T> a, S b){\n  int i, j;\n  Matrix<T> res;\n  res.changeSize(a.r, a.c);\n  res = 1;\n  while(b){\n    if(b&1) res *= a;\n    b >>= 1;\n    a *= a;\n  }\n  return res;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"min_L");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "LIS_length";
      string c = "template<class T>\nint LIS_length(int n, T a[], void *mem = wmem){\n  int i, k, res;\n  T *arr;\n\n  if(n==0) return 0;\n  walloc1d(&arr, n, &mem);\n  arr[0] = a[0];\n  res = 1;\n  REP(i,1,n){\n    k = lower_bound(arr, arr+res, a[i]) - arr;\n    arr[k] = a[i];\n    if(res==k) res++;\n  }\n\n  return res;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "weaklyLIS_length";
      string c = "template<class T>\nint weaklyLIS_length(int n, T a[], void *mem = wmem){\n  int i, k, res;\n  T *arr;\n\n  if(n==0) return 0;\n  walloc1d(&arr, n, &mem);\n  arr[0] = a[0];\n  res = 1;\n  REP(i,1,n){\n    k = upper_bound(arr, arr+res, a[i]) - arr;\n    arr[k] = a[i];\n    if(res==k) res++;\n  }\n\n  return res;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "LIS_ends";
      string c = "template<class T, class S>\nint LIS_ends(int n, T a[], S res[], void *mem = wmem){\n  int i, k, sz;\n  T *arr;\n\n  if(n==0) return 0;\n  walloc1d(&arr, n, &mem);\n  arr[0] = a[0];\n  res[0] = 1;\n  sz = 1;\n  REP(i,1,n){\n    k = lower_bound(arr, arr+sz, a[i]) - arr;\n    arr[k] = a[i];\n    res[i] = k + 1;\n    if(sz==k) sz++;\n  }\n\n  return sz;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "weaklyLIS_ends";
      string c = "template<class T, class S>\nint weaklyLIS_ends(int n, T a[], S res[], void *mem = wmem){\n  int i, k, sz;\n  T *arr;\n\n  if(n==0) return 0;\n  walloc1d(&arr, n, &mem);\n  arr[0] = a[0];\n  res[0] = 1;\n  sz = 1;\n  REP(i,1,n){\n    k = upper_bound(arr, arr+sz, a[i]) - arr;\n    arr[k] = a[i];\n    res[i] = k + 1;\n    if(sz==k) sz++;\n  }\n\n  return sz;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "BIT_popcount";
      string c = "inline int BIT_popcount_L(const int x){\n  return __builtin_popcount(x);\n}\ninline int BIT_popcount_L(const ll x){\n  return __builtin_popcountll(x);\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    
    {
      string n = "BIT_ctz";
      string c = "inline int BIT_ctz_L(const int x){\n  return __builtin_ctz(x);\n}\ninline int BIT_ctz_L(const ll x){\n  return __builtin_ctzll(x);\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "BIT_parity";
      string c = "inline int BIT_parity(const int x){\n  return __builtin_parity(x);\n}\ninline int BIT_parity(const ll x){\n  return __builtin_parityll(x);\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "BIT_parity_pm";
      string c = "inline int BIT_parity_pm(const int x){\n  return 1 - 2 * __builtin_parity(x);\n}\ninline int BIT_parity_pm(const ll x){\n  return 1 - 2 * __builtin_parityll(x);\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }


    {
      string n = "Digit";
      string c = "template<class T> inline int Digit_L(T n){\n  int res = 0;\n  while(n) res++, n /= 10;\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    
    {
      string n = "Digit_base";
      string c = "template<class T, class S> inline int Digit_L(T n, S b){\n  int res = 0;\n  while(n) res++, n /= b;\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    
    {
      string n = "Digit_all";
      string c = "template<class T, class S> inline int Digit_L(T n, S res[]){\n  int sz = 0;\n  while(n) res[sz++] = n % 10, n /= 10;\n  return sz;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    
    {
      string n = "Digit_all_base";
      string c = "template<class T, class S, class U> inline int Digit_L(T n, S res[], U b){\n  int sz = 0;\n  while(n) res[sz++] = n % b, n /= b;\n  return sz;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    
    {
      string n = "DigitF_all";
      string c = "template<class T, class S> inline void DigitF_L(T n, int sz, S res[]){\n  int i;\n  rep(i,sz) res[i] = n % 10, n /= 10;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    
    {
      string n = "DigitF_all_base";
      string c = "template<class T, class S, class U> inline void DigitF_L(T n, int sz, S res[], U b){\n  int i;\n  rep(i,sz) res[i] = n % b, n /= b;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    
    {
      string n = "invDigit";
      string c = "template<class T> inline T invDigit_L(int sz, T d[]){\n  T res = 0;\n  int i;\n  rrep(i,sz) res = 10 * res + d[i];\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    
    {
      string n = "invDigit_base";
      string c = "template<class T, class S> inline T invDigit_L(int sz, T d[], S b){\n  T res = 0;\n  int i;\n  rrep(i,sz) res = b * res + d[i];\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    
    {
      string n = "sod";
      string c = "template<class T> inline int sod_L(T n){\n  int res = 0;\n  while(n) res += n%10, n /= 10;\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    
    {
      string n = "sod_base";
      string c = "template<class T, class S> inline S sod_L(T n, S b){\n  S res = 0;\n  while(n) res += n%b, n /= b;\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "prodDigits";
      string c = "template<class T>\nT prodDigits(T n){\n  T res = 1;\n  if(n==0) return 0;\n  while(n) res *= n % 10, n /= 10;\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "DigitHist";
      string c = "template<class T, class S> inline int DigitHist(T n, S res[]){\n  int i, len = 0;\n  rep(i,10) res[i] = 0;\n  while(n){\n    len++;\n    res[n%10]++;\n    n /= 10;\n  }\n  return len;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "invDigit_r";
      string c = "template<class T> inline T invDigit_r(int sz, T d[]){\n  T res = 0;\n  int i;\n  rep(i,sz) res = 10 * res + d[i];\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "STR2int";
      string c = "inline int STR2int(string s, const int b){\n  int i = 0, fg = 1;\n  int res = 0;\n  if(s[i] == '-') i++, fg = -fg;\n  while(i < s.size()){\n    res = b * res;\n    if('0' <= s[i] <= '9') res += s[i] - '0';\n    else if('A' <= s[i] <= 'Z') res += s[i] - 'A' + 10;\n    else if('a' <= s[i] <= 'z') res += s[i] - 'a' + 36;\n    i++;\n  }\n  return fg * res;\n}\ninline int STR2int(string s){\n  int i = 0, fg = 1;\n  int res = 0;\n  if(s[i] == '-') i++, fg = -fg;\n  while(i < s.size()){\n    res = 10 * res + s[i] - '0';\n    i++;\n  }\n  return fg * res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "STR2ll";
      string c = "inline ll STR2ll(string s, const int b){\n  int i = 0, fg = 1;\n  ll res = 0;\n  if(s[i] == '-') i++, fg = -fg;\n  while(i < s.size()){\n    res = b * res;\n    if('0' <= s[i] <= '9') res += s[i] - '0';\n    else if('A' <= s[i] <= 'Z') res += s[i] - 'A' + 10;\n    else if('a' <= s[i] <= 'z') res += s[i] - 'a' + 36;\n    i++;\n  }\n  return fg * res;\n}\ninline ll STR2ll(string s){\n  int i = 0, fg = 1;\n  ll res = 0;\n  if(s[i] == '-') i++, fg = -fg;\n  while(i < s.size()){\n    res = 10 * res + s[i] - '0';\n    i++;\n  }\n  return fg * res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "Prime";
      string c = "int Prime_L(int N, int res[], void *mem=wmem){\n  int i, a, b;\n  int sz = 1;\n  const int r = 23000;\n  bool *isprime;\n  int *sf, ss = 1;\n  walloc1d(&isprime, r, &mem);\n  walloc1d(&sf, r, &mem);\n  N /= 2;\n  res[0] = 2;\n  b = min(r, N);\n  rep(i,1,b) isprime[i] = 1;\n  rep(i,1,b) if(isprime[i]){\n    res[sz++] = 2i+1;\n    sf[ss] = 2i*(i+1);\n    if(sf[ss] < N){\n      while(sf[ss] < r) isprime[sf[ss]] = 0, sf[ss] += res[ss];\n      ss++;\n    }\n  }\n  for(a=r; a<N; a+=r){\n    b = min(a + r, N);\n    isprime -= r;\n    rep(i,a,b) isprime[i] = 1;\n    rep(i,1,ss){\n      while(sf[i] < b) isprime[sf[i]] = 0, sf[i] += res[i];\n    }\n    rep(i,a,b) if(isprime[i]) res[sz++] = 2i+1;\n  }\n  return sz;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"min_L");
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
    }
    {
      string n = "intervalSieve";
      string c = "template<class T>\nvoid intervalSieve(ll st, int len, T res[], int ps, int p[]){\n  int i;\n  ll k = 2-st;\n  rep(i,len) res[i] = 1;\n  rep(i,k) res[i] = 0;\n  rep(i,ps){\n    k = (ll)p[i]*p[i];\n    if(k >= st+len) break;\n    if(k < st) k = (st+p[i]-1) / p[i] * p[i];\n    while(k < st+len){\n      res[k-st] = 0;\n      k += p[i];\n    }\n  }\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }


    {
      string n = "pow2";
      string c = "template<class T> inline T pow2_L(T a){ return a*a; }";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "pow3";
      string c = "template<class T> inline T pow3_L(T a){ return a*a*a; }";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "pow4";
      string c = "template<class T> inline T pow4_L(T a){ return a*a*a*a; }";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "pow";
      string c = "template<class T, class S> inline T pow_L(T a, S b){\n  T res = 1;\n  res = 1;\n  for(;;){\n    if(b&1) res *= a;\n    b >>= 1;\n    if(b==0) break;\n    a *= a;\n  }\n  return res;\n}\ninline double pow_L(double a, double b){\n  return pow(a,b);\n}\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "powmod";
      string c = "ull powmod(ull a, ull b, ull m){\n  if(m==1) return 0;\n  ull r = 1;\n  while(b){\n    if(b&1) r = r * a % m;\n    b>>=1;\n    if(b) a = a * a % m;\n  }\n  return r;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "PowMod";
      string c = "template<class T, class P, class M>\nT PowMod(T a, P b, M m){\n  T r;\n  r = 1;\n  r %= m;\n  while(b > 0){\n    if(b % 2) r = r * a % m;\n    b /= 2;\n    if(b > 0) a = a * a % m;\n  }\n  return r;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }


    {
      string n = "primitiveRoot";
      string c = "ll primitiveRoot(ll p, void *mem = wmem){\n  int ys;\n  ll *y, r;\n  if(!isPrime(p)) return -1;\n  walloc1d(&y, 1d5, &mem);\n  ys = Divisor(p-1, y, mem);\n  rep(r,1,p){\n    rep(i,ys-1) if(powmod(r,y[i],p) == 1) break_continue;\n    return r;\n  }\n  return -1;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"workmemory");
      d.push_back((string)"walloc1d");
      d.push_back((string)"isPrime");
      d.push_back((string)"Divisor");
      d.push_back((string)"powmod");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "gcd";
      string c = "template<class T, class U> inline T GCD_L(T a, U b){T r; while(b) r=a, a=b, b=r%a; return a;}\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "lcm";
      string c = "template<class T, class U> inline T LCM_L(T a, U b){return a/GCD_L(a,b)*b;}\n";
      string p = "first";
      vector<string> d;

      d.push_back((string)"gcd");
      
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "multiset_popFirst";
      string c = "template<class T> inline T popFirst(multiset<T> &a){\n  T res = *(a.begin());\n  a.erase(a.begin());\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "multiset_getFirst";
      string c = "template<class T> inline T getFirst(multiset<T> &a){\n  return *(a.begin());\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "multiset_popLast";
      string c = "template<class T> inline T popLast(multiset<T> &a){\n  T res;\n  typename multiset<T>::iterator it;\n  it = a.end();\n  it--;\n  res = *it;\n  a.erase(it);\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "multiset_getLast";
      string c = "template<class T> inline T getLast(multiset<T> &a){\n  typename multiset<T>::iterator it;\n  it = a.end();\n  it--;\n  return *it;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "set_popFirst";
      string c = "template<class T> inline T popFirst(set<T> &a){\n  T res = *(a.begin());\n  a.erase(a.begin());\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "set_getFirst";
      string c = "template<class T> inline T getFirst(set<T> &a){\n  return *(a.begin());\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "set_popLast";
      string c = "template<class T> inline T popLast(set<T> &a){\n  T res;\n  typename set<T>::iterator it;\n  it = a.end();\n  it--;\n  res = *it;\n  a.erase(it);\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "set_getLast";
      string c = "template<class T> inline T getLast(set<T> &a){\n  typename set<T>::iterator it;\n  it = a.end();\n  it--;\n  return *it;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "multiset_getClosest";
      string c = "template<class T, class S>\nT getClosest(multiset<T> &s, S v){\n  typename multiset<T>::iterator it, jt;\n  it = s.lower_bound(v);\n  if(it==s.begin()) return *it;\n  if(it==s.end()) return *(--it);\n  jt = it;\n  it--;\n  if(v - (*it) <= (*jt) - v) return *it;\n  return *jt;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "multiset_getClosestL";
      string c = "template<class T, class S>\nT getClosestL(multiset<T> &s, S v){\n  typename multiset<T>::iterator it, jt;\n  it = s.lower_bound(v);\n  if(it==s.begin()) return *it;\n  if(it==s.end()) return *(--it);\n  jt = it;\n  it--;\n  if(v - (*it) < (*jt) - v) return *it;\n  return *jt;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "set_getClosest";
      string c = "template<class T, class S>\nT getClosest(set<T> &s, S v){\n  typename set<T>::iterator it, jt;\n  it = s.lower_bound(v);\n  if(it==s.begin()) return *it;\n  if(it==s.end()) return *(--it);\n  jt = it;\n  it--;\n  if(v - (*it) <= (*jt) - v) return *it;\n  return *jt;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "set_getClosestL";
      string c = "template<class T, class S>\nT getClosestL(set<T> &s, S v){\n  typename set<T>::iterator it, jt;\n  it = s.lower_bound(v);\n  if(it==s.begin()) return *it;\n  if(it==s.end()) return *(--it);\n  jt = it;\n  it--;\n  if(v - (*it) < (*jt) - v) return *it;\n  return *jt;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "arr_getClosest";
      string c = "template<class T, class S>\nT getClosest(int N, T A[], S val){\n  int i;\n  i = lower_bound(A, A+N, val) - A;\n  if(i==0) return A[0];\n  if(i==N) return A[N-1];\n  if(val - A[i-1] <= A[i] - val) return A[i-1];\n  return A[i];\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "arr_getClosestL";
      string c = "template<class T, class S>\nT getClosestL(int N, T A[], S val){\n  int i;\n  i = lower_bound(A, A+N, val) - A;\n  if(i==0) return A[0];\n  if(i==N) return A[N-1];\n  if(val - A[i-1] < A[i] - val) return A[i-1];\n  return A[i];\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }


    {
      string n = "combination_mint";
      string c = "struct combination_mint{\n  mint *fac, *ifac;\n\n  void init(int n, void **mem = &wmem){\n    int i;\n    \n    walloc1d(&fac, n, mem);\n    walloc1d(&ifac, n, mem);\n    \n    fac[0] = 1;\n    rep(i,1,n) fac[i] = fac[i-1] * i;\n    ifac[n-1] = 1 / fac[n-1];\n    for(i=n-2;i>=0;i--) ifac[i] = ifac[i+1] * (i+1);\n  }\n\n  mint C(int a, int b){\n    if(b < 0 || b > a) return 0;\n    return fac[a]*ifac[b]*ifac[a-b];\n  }\n\n  mint P(int a, int b){\n    if(b < 0 || b > a) return 0;\n    return fac[a]*ifac[a-b];\n  }\n\n  mint H(int a, int b){\n    if(a==0 && b==0) return 1;\n    if(a<=0 || b<0) return 0;\n    return C(a+b-1, b);\n  }\n};\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"mint");
      d.push_back((string)"workmemory");
      d.push_back((string)"walloc1d");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "reduceFraction";
      string c = "template<class T, class U> void reduceFraction(T&a, U&b){T g=GCD_L(a,b);a/=g;b/=g;}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"gcd");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "runLength";
      string c = "template<class T>\nint runLength(int N, T *arr, void *val_s = NULL, void *len_s = NULL){\n  int i, rN;\n  T *val = (T*) val_s;\n  int *len = (int*) len_s;\n  if(N==0) return 0;\n  if(val==NULL || len==NULL){\n    void *mem = wmem;\n    if(val==NULL) walloc1d(&val, N, &mem);\n    if(len==NULL) walloc1d(&len, N, &mem);\n  }\n  rN = 1;\n  val[0] = arr[0];\n  len[0] = 1;\n  rep(i,1,N){\n    if(val[rN-1] == arr[i]){\n      len[rN-1]++;\n    } else {\n      val[rN] = arr[i];\n      len[rN] = 1;\n      rN++;\n    }\n  }\n  return rN;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"workmemory");
      d.push_back((string)"walloc1d");
      need[n] = d;
    }


    {
      string n = "arrcmp";
      string c = "template<class S, class T> inline int arrcmp(int As, S A[], int Bs, T B[]){\n  int i;\n  for(i=0;;i++){\n    if(i==As==Bs) break;\n    if(i==As) return -1;\n    if(i==Bs) return 1;\n    if(A[i] < B[i]) return -1;\n    if(A[i] > B[i]) return 1;\n  }\n  return 0;\n}\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "arrErase";
      string c = "template<class S>\nS arrErase(int k, int &sz, S a[]){\n  int i;\n  S res;\n  res = a[k];\n  sz--;\n  rep(i,k,sz) a[i] = a[i+1];\n  return res;\n}\ntemplate<class S, class T>\nvoid arrErase(int k, int &sz, S a[], T b[]){\n  int i;\n  sz--;\n  rep(i,k,sz) a[i] = a[i+1];\n  rep(i,k,sz) b[i] = b[i+1];\n}\ntemplate<class S, class T, class U>\nvoid arrErase(int k, int &sz, S a[], T b[], U c[]){\n  int i;\n  sz--;\n  rep(i,k,sz) a[i] = a[i+1];\n  rep(i,k,sz) b[i] = b[i+1];\n  rep(i,k,sz) c[i] = c[i+1];\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }


    {
      string n = "arrEraseVal1";
      string c = "template<class S>\ninline int arrEraseVal_L(S val1, int &sz, S a[]){\n  int i, n = sz;\n  sz = 0;\n  rep(i,n) if(a[i]!=val1) a[sz++] = a[i];\n  return n - sz;\n}\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "arrEraseVal2";
      string c = "template<class S>\ninline int arrEraseVal_L(S val1, S val2, int &sz, S a[]){\n  int i, n = sz;\n  sz = 0;\n  rep(i,n) if(a[i]!=val1 && a[i]!=val2) a[sz++] = a[i];\n  return n - sz;\n}\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "arrEraseVal3";
      string c = "template<class S>\ninline int arrEraseVal_L(S val1, S val2, S val3, int &sz, S a[]){\n  int i, n = sz;\n  sz = 0;\n  rep(i,n) if(a[i]!=val1 && a[i]!=val2 && a[i]!=val3) a[sz++] = a[i];\n  return n - sz;\n}\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "arrEraseVal4";
      string c = "template<class S>\ninline int arrEraseVal_L(S val1, S val2, S val3, S val4, int &sz, S a[]){\n  int i, n = sz;\n  sz = 0;\n  rep(i,n) if(a[i]!=val1 && a[i]!=val2 && a[i]!=val3 && a[i]!=val4) a[sz++] = a[i];\n  return n - sz;\n}\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "vecEraseVal1";
      string c = "template<class S>\ninline int vecEraseVal_L(S val1, vector<S> &a){\n  int i, k = 0, n = a.size();\n  rep(i,n) if(a[i]!=val1) a[k++] = a[i];\n  i = n - k;\n  while(a.size() > k) a.pop_back();\n  return i;\n}\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "Mex";
      string c = "template<class T>\nint Mex_L(int n, T a[], int sorted = 0, void *mem = wmem){\n  int i, k = 0;\n  if(sorted==0){\n    T *arr;\n    walloc1d(&arr, n, &mem);\n    rep(i,n) arr[i] = a[i];\n    sort(arr, arr+n);\n    rep(i,n) if(arr[i]==k) k++;\n  } else {\n    rep(i,n) if(a[i]==k) k++;\n  }\n  return k;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"workmemory");
      d.push_back((string)"walloc1d");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "Mex1";
      string c = "template<class T>\ninline int Mex_L(T x){\n  int k = 0;\n  while(k==x) k++;\n  return k;\n}\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "Mex2";
      string c = "template<class T1, class T2>\ninline int Mex_L(T1 x, T2 y){\n  int k = 0;\n  while(k==x || k==y) k++;\n  return k;\n}\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "Mex3";
      string c = "template<class T1, class T2, class T3>\ninline int Mex_L(int x, int y, int z){\n  int k = 0;\n  while(k==x || k==y || k==z) k++;\n  return k;\n}\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "Distinct";
      string c = "template<class T>\nint Distinct(int N, T A[], int sorted=0, void *mem = wmem){\n  int i, k, res = 1;\n  T *a;\n  if(N==0) return 0;\n  if(sorted){\n    rep(i,1,N) if(A[i]!=A[i-1]) res++;\n  } else {\n    walloc1d(&a,N,&mem);\n    rep(i,N) a[i] = A[i];\n    sort(a,a+N);\n    rep(i,1,N) if(a[i]!=a[i-1]) res++;\n  }\n  return res;\n}\ntemplate<class T>\nint Distinct(vector<T> A, int sorted=0, void *mem = wmem){\n  int i, k, res = 1, N = A.size();\n  T *a;\n  if(N==0) return 0;\n  if(sorted){\n    rep(i,1,N) if(A[i]!=A[i-1]) res++;\n  } else {\n    walloc1d(&a,N,&mem);\n    rep(i,N) a[i] = A[i];\n    sort(a,a+N);\n    rep(i,1,N) if(a[i]!=a[i-1]) res++;\n  }\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"workmemory");
      d.push_back((string)"walloc1d");
      need[n] = d;
    }


    {
      string n = "DistinctE_2";
      string c = "template<class T1, class T2>\ninline int DistinctE_L(T1 a, T2 b){\n  if(a != b) return 2;\n  return 1;\n}\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "DistinctE_3";
      string c = "template<class T1, class T2, class T3>\ninline int DistinctE_L(T1 a, T2 b, T3 c){\n  if(a==b){\n    if(a==c) return 1;\n    return 2;\n  }\n  if(a==c || b==c) return 2;\n  return 3;\n}\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "DistinctE_4";
      string c = "template<class T1, class T2, class T3, class T4>\ninline int DistinctE_L(T1 a, T2 b, T3 c, T4 d){\n  if(a==b){\n    if(a==c){\n      if(a==d) return 1;\n      return 2;\n    }\n    if(a==d || c==d) return 2;\n    return 3;\n  }\n  if(a==c){\n    if(a==d || b==d) return 2;\n    return 3;\n  }\n  if(b==c){\n    if(a==d || b==d) return 2;\n    return 3;\n  }\n  if(a==d || b==d || c==d) return 3;\n  return 4;\n}\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "Count";
      string c = "template<class T, class S>\nint Count(const int N, const T A[], const S val){\n  int res = 0;\n  rep(i,N) if(A[i]==val) res++;\n  return res;\n}\ntemplate<class T, class S>\nint Count(const vector<T> &A, const S val){\n  int res = 0;\n  rep(i,A.size()) if(A[i]==val) res++;\n  return res;\n}\ntemplate<class S>\nint Count(const string &A, const S val){\n  int res = 0;\n  rep(i,A.size()) if(A[i]==val) res++;\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "arrCountVal";
      string c = "template<class T, class S>\nint arrCountVal(const int N, const T A[], const S val){\n  int res = 0;\n  rep(i,N) if(A[i]==val) res++;\n  return res;\n}\ntemplate<class T, class S>\nint arrCountVal(const vector<T> &A, const S val){\n  int res = 0;\n  rep(i,A.size()) if(A[i]==val) res++;\n  return res;\n}\ntemplate<class S>\nint arrCountVal(const string &A, const S val){\n  int res = 0;\n  rep(i,A.size()) if(A[i]==val) res++;\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "arrCountValSeqMax";
      string c = "template<class T, class S>\nint arrCountValSeqMax(const int N, const T A[], const S val){\n  int res = 0, cur = 0;\n  rep(i,N){\n    if(A[i]==val){\n      cur++;\n      if(res < cur) res = cur;\n    } else {\n      cur = 0;\n    }\n  }\n  return res;\n}\ntemplate<class S>\nint arrCountValSeqMax(const string A, const S val){\n  int res = 0, cur = 0;\n  rep(i,A.size()){\n    if(A[i]==val){\n      cur++;\n      if(res < cur) res = cur;\n    } else {\n      cur = 0;\n    }\n  }\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "PalindromeCost";
      string c = "template<class T>\nint PalindromeCost(const int N, const T A[]){\n  int i, res = 0;\n  rep(i,N/2) if(A[i] != A[N-1-i]) res++;\n  return res;\n}\ntemplate<class T>\nint PalindromeCost(const vector<T> &A){\n  const int N = A.size();\n  int i, res = 0;\n  rep(i,N/2) if(A[i] != A[N-1-i]) res++;\n  return res;\n}\nint PalindromeCost(const string &A){\n  const int N = A.size();\n  int i, res = 0;\n  rep(i,N/2) if(A[i] != A[N-1-i]) res++;\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }


    {
      string n = "arrInsert";
      string c = "template<class S>\ninline void arrInsert(const int k, int &sz, S a[], const S aval){\n  int i;\n  sz++;\n  for(i=sz-1;i>k;i--) a[i] = a[i-1];\n  a[k] = aval;\n}\n\ntemplate<class S, class T>\ninline void arrInsert(const int k, int &sz, S a[], const S aval, T b[], const T bval){\n  int i;\n  sz++;\n  for(i=sz-1;i>k;i--) a[i] = a[i-1];\n  for(i=sz-1;i>k;i--) b[i] = b[i-1];\n  a[k] = aval;\n  b[k] = bval;\n}\n\ntemplate<class S, class T, class U>\ninline void arrInsert(const int k, int &sz, S a[], const S aval, T b[], const T bval, U c[], const U cval){\n  int i;\n  sz++;\n  for(i=sz-1;i>k;i--) a[i] = a[i-1];\n  for(i=sz-1;i>k;i--) b[i] = b[i-1];\n  for(i=sz-1;i>k;i--) c[i] = c[i-1];\n  a[k] = aval;\n  b[k] = bval;\n  c[k] = cval;\n}\n\ntemplate<class S, class T, class U, class V>\ninline void arrInsert(const int k, int &sz, S a[], const S aval, T b[], const T bval, U c[], const U cval, V d[], const V dval){\n  int i;\n  sz++;\n  for(i=sz-1;i>k;i--) a[i] = a[i-1];\n  for(i=sz-1;i>k;i--) b[i] = b[i-1];\n  for(i=sz-1;i>k;i--) c[i] = c[i-1];\n  for(i=sz-1;i>k;i--) d[i] = d[i-1];\n  a[k] = aval;\n  b[k] = bval;\n  c[k] = cval;\n  d[k] = dval;\n}\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "floor_sum";
      string c = "ll floor_sum(ll n, ll m, ll a, ll b){\n  ll res = 0, s, t;\n  if(m < 0){\n    m = -m;\n    a = -a;\n    b = -b;\n  }\n\n  t = fDiv(a, m);\n  if(t){\n    res += (n-1) * n / 2 * t;\n    a -= t * m;\n  }\n  t = fDiv(b, m);\n  if(t){\n    res += n * t;\n    b -= t * m;\n  }\n  for(;;){\n    s = (a * n + b) / m;\n    if(s==0) break;\n    t = s * m - b;\n    res += (n - (t + a - 1) / a) * s;\n\n    n = s;\n    b = a - t % a;\n    if(b == a) b = 0;\n    swap(m,a);\n    \n    if(a >= m){\n      res += (n-1) * n / 2 * (a / m);\n      a %= m;\n    }\n    if(b >= m){\n      res += n * (b / m);\n      b %= m;\n    }\n  }\n  return res;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"fDiv");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "floor_sum2";
      string c = "ll floor_sum2(ll a, ll k){\n  ll i = 1, j, v, res = 0;\n  if(k > a) k = a;\n  while(i <= k){\n    v = a / i;\n    j = min(k, a / v) + 1;\n    res += v * (j-i);\n    i = j;\n  }\n  return res;\n}\n\nll floor_sum2(ll a){\n  ll i = 1, j, v, res = 0;\n  while(i <= a){\n    v = a / i;\n    j = a / v + 1;\n    res += v * (j-i);\n    i = j;\n  }\n  return res;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"min_L");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "Comb";
      string c = "template<class T>\nstruct Comb{\n  int mem_fact;\n  T *factri, *ifactri;\n  int mem_dfact;\n  T *dfactri;\n  int mem_pw2, mem_pw3, mem_pw10, mem_rep1;\n  T *pw2c, *pw3c, *pw10c, *rep1c;\n  int mem_ipw2, mem_ipw3, mem_ipw10;\n  T *ipw2c, *ipw3c, *ipw10c;\n  Comb(){\n    mem_fact = 0;\n    mem_dfact = 0;\n    mem_pw2 = mem_pw3 = mem_pw10 = mem_rep1 = 0;\n    mem_ipw2 = mem_ipw3 = mem_ipw10 = 0;\n  }\n  inline void expand_fact(int k){\n    int i;\n    if(k <= mem_fact) return;\n    k >?= 2 * mem_fact;\n    if(mem_fact == 0){\n      factri = (T*)malloc(k * sizeof(T));\n      ifactri = (T*)malloc(k * sizeof(T));\n      factri[0] = 1;\n      rep(i,1,k) factri[i] = i * factri[i-1];\n      ifactri[k-1] = 1 / factri[k-1];\n      rrep(i,k-1) ifactri[i] = (i+1) * ifactri[i+1];\n    } else {\n      factri = (T*)realloc(factri, k * sizeof(T));\n      ifactri = (T*)realloc(ifactri, k * sizeof(T));\n      rep(i,mem_fact,k) factri[i] = i * factri[i-1];\n      ifactri[k-1] = 1 / factri[k-1];\n      rrep(i,mem_fact,k-1) ifactri[i] = (i+1) * ifactri[i+1];\n    }\n    mem_fact = k;\n  }\n  inline T fac(int k){\n    if(mem_fact < k+1) expand_fact(k+1);\n    return factri[k];\n  }\n  inline T ifac(int k){\n    if(mem_fact < k+1) expand_fact(k+1);\n    return ifactri[k];\n  }\n  inline T C(int a, int b){\n    if(b < 0 || b > a) return 0;\n    if(mem_fact < a+1) expand_fact(a+1);\n    return factri[a] * ifactri[b] * ifactri[a-b];\n  }\n  inline T P(int a, int b){\n    if(b < 0 || b > a) return 0;\n    if(mem_fact < a+1) expand_fact(a+1);\n    return factri[a] * ifactri[a-b];\n  }\n  inline T H(int a, int b){\n    if(a==0 && b==0) return 1;\n    if(a <= 0 || b < 0) return 0;\n    if(mem_fact < a+b) expand_fact(a+b);\n    return C(a+b-1, b);\n  }\n  inline T Multinomial(int sz, int a[]){\n    int i, s = 0;\n    T res;\n    rep(i,sz) s += a[i];\n    if(mem_fact < s+1) expand_fact(s+1);\n    res = factri[s];\n    rep(i,sz) res *= ifactri[a[i]];\n    return res;\n  }\n  inline T Multinomial(int a){\n    return 1;\n  }\n  inline T Multinomial(int a, int b){\n    if(mem_fact < a+b+1) expand_fact(a+b+1);\n    return factri[a+b] * ifactri[a] * ifactri[b];\n  }\n  inline T Multinomial(int a, int b, int c){\n    if(mem_fact < a+b+c+1) expand_fact(a+b+c+1);\n    return factri[a+b+c] * ifactri[a] * ifactri[b] * ifactri[c];\n  }\n  inline T Multinomial(int a, int b, int c, int d){\n    if(mem_fact < a+b+c+d+1) expand_fact(a+b+c+d+1);\n    return factri[a+b+c+d] * ifactri[a] * ifactri[b] * ifactri[c] * ifactri[d];\n  }\n  inline T Catalan(int n){\n    if(n < 0) return 0;\n    if(mem_fact < 2*n+1) expand_fact(2*n+1);\n    return factri[2*n] * ifactri[n] * ifactri[n+1];\n  }\n  inline T Catalan(int n, int m, int k){\n    if(k <= 0) return C(n+m, n);\n    if(n < k || m < k) return 0;\n    return C(n+m, m) - C(n+m, k-1);\n  }\n  inline T Catalan_s(ll n, ll m, ll k){\n    if(k <= 0) return C_s(n+m, n);\n    if(n < k || m < k) return 0;\n    return C_s(n+m, m) - C_s(n+m, k-1);\n  }\n  inline T C_s(ll a, ll b){\n    ll i;\n    T res;\n    if(b < 0 || b > a) return 0;\n    if(b > a - b) b = a - b;\n    res = 1;\n    rep(i,b){\n      res *= a - i;\n      res /= i + 1;\n    }\n    return res;\n  }\n  inline T P_s(ll a, ll b){\n    ll i;\n    T res;\n    if(b < 0 || b > a) return 0;\n    res = 1;\n    rep(i,b) res *= a - i;\n    return res;\n  }\n  inline T H_s(ll a, ll b){\n    if(a==0 && b==0) return 1;\n    if(a <= 0 || b < 0) return 0;\n    return C_s(a+b-1, b);\n  }\n  inline T per_s(ll n, ll k){\n    T d;\n    int m;\n    if(n < 0 || k < 0) return 0;\n    if(n == k == 0) return 1;\n    if(n == 0 || k == 0) return 0;\n    if(k==1) return 1;\n    if(k==2){\n      d = n / 2;\n      return d;\n    }\n    if(k==3){\n      d = (n-1) / 6;\n      m = (n-1) % 6;\n      if(m==0) return 3 * d * d + d;\n      if(m==1) return 3 * d * d + 2 * d;\n      if(m==2) return 3 * d * d + 3 * d + 1;\n      if(m==3) return 3 * d * d + 4 * d + 1;\n      if(m==4) return 3 * d * d + 5 * d + 2;\n      if(m==5) return 3 * d * d + 6 * d + 3;\n    }\n    assert(0 && \"per_s should be k <= 3\");\n    return -1;\n  }\n  inline void expand_dfact(int k){\n    int i;\n    if(k <= mem_dfact) return;\n    k >?= 3;\n    k >?= 2 * mem_dfact;\n    if(mem_dfact==0){\n      dfactri = (T*)malloc(k * sizeof(T));\n      dfactri[0] = dfactri[1] = 1;\n      rep(i,2,k) dfactri[i] = i * dfactri[i-2];\n    } else {\n      dfactri = (T*)realloc(dfactri, k * sizeof(T));\n      rep(i,mem_dfact,k) dfactri[i] = i * dfactri[i-2];\n    }\n    mem_dfact = k;\n  }\n  inline void expand_pw2(int k){\n    int i;\n    if(k <= mem_pw2) return;\n    k >?= 2 * mem_pw2;\n    if(mem_pw2==0){\n      pw2c = (T*)malloc(k * sizeof(T));\n      pw2c[0] = 1;\n      rep(i,1,k) pw2c[i] = 2 * pw2c[i-1];\n    } else {\n      pw2c = (T*)realloc(pw2c, k * sizeof(T));\n      rep(i,mem_pw2,k) pw2c[i] = 2 * pw2c[i-1];\n    }\n    mem_pw2 = k;\n  }\n  inline void expand_ipw2(int k){\n    int i;\n    if(k <= mem_ipw2) return;\n    k >?= 2;\n    k >?= 2 * mem_ipw2;\n    if(mem_ipw2==0){\n      ipw2c = (T*)malloc(k * sizeof(T));\n      ipw2c[0] = 1;\n      ipw2c[1] = ipw2c[0] / 2;\n      rep(i,1,k) ipw2c[i] = ipw2c[1] * ipw2c[i-1];\n    } else {\n      ipw2c = (T*)realloc(ipw2c, k * sizeof(T));\n      rep(i,mem_ipw2,k) ipw2c[i] = ipw2c[1] * ipw2c[i-1];\n    }\n    mem_ipw2 = k;\n  }\n  inline void expand_pw3(int k){\n    int i;\n    if(k <= mem_pw3) return;\n    k >?= 2 * mem_pw3;\n    if(mem_pw3==0){\n      pw3c = (T*)malloc(k * sizeof(T));\n      pw3c[0] = 1;\n      rep(i,1,k) pw3c[i] = 3 * pw3c[i-1];\n    } else {\n      pw3c = (T*)realloc(pw3c, k * sizeof(T));\n      rep(i,mem_pw3,k) pw3c[i] = 3 * pw3c[i-1];\n    }\n    mem_pw3 = k;\n  }\n  inline void expand_ipw3(int k){\n    int i;\n    if(k <= mem_ipw3) return;\n    k >?= 2;\n    k >?= 2 * mem_ipw3;\n    if(mem_ipw3==0){\n      ipw3c = (T*)malloc(k * sizeof(T));\n      ipw3c[0] = 1;\n      ipw3c[1] = ipw3c[0] / 3;\n      rep(i,1,k) ipw3c[i] = ipw3c[1] * ipw3c[i-1];\n    } else {\n      ipw3c = (T*)realloc(ipw3c, k * sizeof(T));\n      rep(i,mem_ipw3,k) ipw3c[i] = ipw3c[1] * ipw3c[i-1];\n    }\n    mem_ipw3 = k;\n  }\n  inline void expand_pw10(int k){\n    int i;\n    if(k <= mem_pw10) return;\n    k >?= 2 * mem_pw10;\n    if(mem_pw10==0){\n      pw10c = (T*)malloc(k * sizeof(T));\n      pw10c[0] = 1;\n      rep(i,1,k) pw10c[i] = 10 * pw10c[i-1];\n    } else {\n      pw10c = (T*)realloc(pw10c, k * sizeof(T));\n      rep(i,mem_pw10,k) pw10c[i] = 10 * pw10c[i-1];\n    }\n    mem_pw10 = k;\n  }\n  inline void expand_ipw10(int k){\n    int i;\n    if(k <= mem_ipw10) return;\n    k >?= 2;\n    k >?= 2 * mem_ipw10;\n    if(mem_ipw10==0){\n      ipw10c = (T*)malloc(k * sizeof(T));\n      ipw10c[0] = 1;\n      ipw10c[1] = ipw10c[0] / 10;\n      rep(i,1,k) ipw10c[i] = ipw10c[1] * ipw10c[i-1];\n    } else {\n      ipw10c = (T*)realloc(ipw10c, k * sizeof(T));\n      rep(i,mem_ipw10,k) ipw10c[i] = ipw10c[1] * ipw10c[i-1];\n    }\n    mem_ipw10 = k;\n  }\n  inline void expand_rep1(int k){\n    int i;\n    if(k <= mem_rep1) return;\n    k >?= 2 * mem_rep1;\n    if(mem_rep1==0){\n      rep1c = (T*)malloc(k * sizeof(T));\n      rep1c[0] = 0;\n      rep(i,1,k) rep1c[i] = 10 * rep1c[i-1] + 1;\n    } else {\n      rep1c = (T*)realloc(rep1c, k * sizeof(T));\n      rep(i,mem_rep1,k) rep1c[i] = 10 * rep1c[i-1] + 1;\n    }\n    mem_rep1 = k;\n  }\n  inline T dfac(int k){\n    if(k >= 0){\n      if(mem_dfact < k+1) expand_dfact(k+1);\n      return dfactri[k];\n    }\n    if(k==-1) return 1;\n    k = - k - 2;\n    if(k % 4 == 1) return 1 / (-dfac(k));\n    return 1 / dfac(k);\n  }\n  inline T pw2(int k){\n    if(k >= 0){\n      if(mem_pw2 < k+1) expand_pw2(k+1);\n      return pw2c[k];\n    } else {\n      k = -k;\n      if(mem_ipw2 < k+1) expand_ipw2(k+1);\n      return ipw2c[k];\n    }\n  }\n  inline T pw3(int k){\n    if(k >= 0){\n      if(mem_pw3 < k+1) expand_pw3(k+1);\n      return pw3c[k];\n    } else {\n      k = -k;\n      if(mem_ipw3 < k+1) expand_ipw3(k+1);\n      return ipw3c[k];\n    }\n  }\n  inline T pw10(int k){\n    if(k >= 0){\n      if(mem_pw10 < k+1) expand_pw10(k+1);\n      return pw10c[k];\n    } else {\n      k = -k;\n      if(mem_ipw10 < k+1) expand_ipw10(k+1);\n      return ipw10c[k];\n    }\n  }\n  inline T repunit(int k){\n    if(mem_rep1 < k+1) expand_rep1(k+1);\n    return rep1c[k];\n  }\n};\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"chmax");
      d.push_back((string)"workmemory");
      d.push_back((string)"walloc1d");
      d.push_back((string)"Comb_Cs_Modint");
      d.push_back((string)"Comb_Cs_modint");
      d.push_back((string)"Comb_Cs_Mint");
      d.push_back((string)"Comb_Cs_mint");
      need[n] = d;
    }
    {
      string n = "Comb_Cs_Modint";
      string c = "template<>\ninline Modint Comb<Modint>::C_s(ll a, ll b){\n  ll i;\n  Modint res, d;\n  if(b < 0 || b > a) return 0;\n  if(b > a - b) b = a - b;\n  res = d = 1;\n  rep(i,b){\n    res *= a - i;\n    d *= i + 1;\n  }\n  return res / d;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Modint";
    }
    {
      string n = "Comb_Cs_modint";
      string c = "template<>\ninline modint Comb<modint>::C_s(ll a, ll b){\n  ll i;\n  modint res, d;\n  if(b < 0 || b > a) return 0;\n  if(b > a - b) b = a - b;\n  res = d = 1;\n  rep(i,b){\n    res *= a - i;\n    d *= i + 1;\n  }\n  return res / d;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "modint";
    }
    {
      string n = "Comb_Cs_Mint";
      string c = "template<>\ninline Mint Comb<Mint>::C_s(ll a, ll b){\n  ll i;\n  Mint res, d;\n  if(b < 0 || b > a) return 0;\n  if(b > a - b) b = a - b;\n  res = d = 1;\n  rep(i,b){\n    res *= a - i;\n    d *= i + 1;\n  }\n  return res / d;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Mint";
    }
    {
      string n = "Comb_Cs_mint";
      string c = "template<>\ninline mint Comb<mint>::C_s(ll a, ll b){\n  ll i;\n  mint res, d;\n  if(b < 0 || b > a) return 0;\n  if(b > a - b) b = a - b;\n  res = d = 1;\n  rep(i,b){\n    res *= a - i;\n    d *= i + 1;\n  }\n  return res / d;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "mint";
    }


    {
      string n = "arrRot";
      string c = "template<class T>\nvoid arrRot(int k, int N, T A[], T B[] = NULL, void *mem = wmem){\n  int fg = 0;\n  k %%= N;\n  if(B==NULL){\n    walloc1d(&B, N, &mem);\n    fg = 1;\n  }\n  rep(i,k,N) B[i-k] = A[i];\n  rep(i,k) B[N-k+i] = A[i];\n  if(fg) rep(i,N) A[i] = B[i];\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"workmemory");
      d.push_back((string)"walloc1d");
      d.push_back((string)"moddw");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "Polynomial";
      string c = "template<class T>\nstruct Polynomial {\n  int d, mem;\n  T *c;\n  Polynomial(){\n    mem = 1;\n    c = new T[mem];\n    d = 0;\n    c[0] = 0;\n  }\n  Polynomial(T a){\n    mem = 1;\n    c = new T[mem];\n    d = 0;\n    c[0] = a;\n  }\n  Polynomial(const Polynomial<T> &a){\n    d = a.d;\n    mem = d + 1;\n    c = new T[mem];\n    rep(i,d+1) c[i] = a.c[i];\n  }\n  ~Polynomial(){\n    delete [] c;\n  }\n  void expand(int z){\n    T *cc;\n    if(z <= mem) return;\n    mem = max(z, 2 * mem);\n    cc = new T[mem];\n    rep(i,d+1) cc[i] = c[i];\n    delete [] c;\n    c = cc;\n  }\n  inline void change(const int dg, const T cf){\n    expand(dg+1);\n    while(d < dg) c[++d] = 0;\n    c[dg] = cf;\n    while(d && c[d]==0) d--;\n  }\n  inline int deg(void){\n    return d;\n  }\n  inline T coef(const int k){\n    if(k > d) return 0;\n    return c[k];\n  }\n  Polynomial<T>& operator=(const T a){\n    d = 0;\n    expand(d + 1);\n    c[0] = a;\n    return *this;\n  }\n  Polynomial<T>& operator=(const Polynomial<T> &a){\n    int i;\n    d = a.d;\n    expand(d + 1);\n    rep(i,d+1) c[i] = a.c[i];\n    return *this;\n  }\n  Polynomial<T>& operator+=(const Polynomial<T> &a){\n    int i, k;\n    k = max(d, a.d);\n    expand(k+1);\n    while(d < k) c[++d] = 0;\n    rep(i,a.d+1) c[i] += a.c[i];\n    while(d && c[d]==0) d--;\n    return *this;\n  }\n  Polynomial<T> operator+(const Polynomial<T> &a){\n    return Polynomial<T>(*this) += a;\n  }\n  Polynomial<T>& operator-=(const Polynomial<T> &a){\n    int i, k;\n    k = max(d, a.d);\n    expand(k+1);\n    while(d < k) c[++d] = 0;\n    rep(i,a.d+1) c[i] -= a.c[i];\n    while(d && c[d]==0) d--;\n    return *this;\n  }\n  Polynomial<T> operator-(const Polynomial<T> &a){\n    return Polynomial<T>(*this) -= a;\n  }\n  Polynomial<T>& operator*=(const Polynomial<T> &a){\n    int i, j, k;\n    T *cc;\n    void *mem = wmem;\n    \n    k = d + a.d;\n    expand(k+1);\n    walloc1d(&cc, k+1, &mem);\n    rep(i,k+1) cc[i] = 0;\n    rep(i,d+1) rep(j,a.d+1) cc[i+j] += c[i] * a.c[j];\n    rep(i,k+1) c[i] = cc[i];\n    d = k;\n    while(d && c[d]==0) d--;\n    return *this;\n  }\n  Polynomial<T> operator*(const Polynomial<T> &a){\n    return Polynomial<T>(*this) *= a;\n  }\n  Polynomial<T>& operator/=(const Polynomial<T> &a){\n    int i, j, k;\n    T *cc, e;\n    void *mem = wmem;\n    walloc1d(&cc, d-a.d, &mem);\n    for(i=d; i>=a.d; i--){\n      cc[i-a.d] = e = c[i] / a.c[a.d];\n      rep(j, a.d+1) c[i-j] -= e * a.c[a.d-j];\n    }\n    d -= a.d;\n    rep(i,d+1) c[i] = cc[i];\n    return *this;\n  }\n  Polynomial<T> operator/(const Polynomial<T> &a){\n    return Polynomial<T>(*this) /= a;\n  }\n  Polynomial<T>& operator%=(const Polynomial<T> &a){\n    int i, j, k;\n    T *cc, e;\n    void *mem = wmem;\n    walloc1d(&cc, d-a.d, &mem);\n    for(i=d; i>=a.d; i--){\n      cc[i-a.d] = e = c[i] / a.c[a.d];\n      rep(j, a.d+1) c[i-j] -= e * a.c[a.d-j];\n    }\n    while(d && c[d]==0) d--;\n    return *this;\n  }\n  Polynomial<T> operator%(const Polynomial<T> &a){\n    return Polynomial<T>(*this) %= a;\n  }\n  Polynomial<T>& operator*=(const T &a){\n    rep(i,d+1) c[i] *= a;\n    while(d && c[d]==0) d--;\n    return *this;\n  }\n  Polynomial<T> operator*(const T &a){\n    return Polynomial<T>(*this) *= a;\n  }\n  Polynomial<T>& operator/=(const T &a){\n    rep(i,d+1) c[i] /= a;\n    while(d && c[d]==0) d--;\n    return *this;\n  }\n  Polynomial<T> operator/(const T &a){\n    return Polynomial<T>(*this) /= a;\n  }\n  inline T operator()(const T x){\n    int i;\n    T res;\n    res = 0;\n    for(i=d;i>=0;i--) res = res * x + c[i];\n    return res;\n  }\n};\ntemplate<class T> Polynomial<T> operator*(const T a, const Polynomial<T> &b){return Polynomial<T>(b)*=a;}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"chmax");
      d.push_back((string)"max_L");
      need[n] = d;
    }


    {
      string n = "SuffixArray";
      string c = "template<class T> void SuffixArray(T *s, int N, int K, int *SA, int *LCP = NULL, void *mem = wmem) {\n  int i, j, d, m, *s1;\n  int name, prev, pos;\n  char *t, *lms;\n  int *cnt, *cnt1, *cnt2;\n  \n  walloc1d(&t, N+1, &mem);\n  walloc1d(&lms, N+1, &mem);\n  walloc1d(&cnt, K+1, &mem);\n  walloc1d(&cnt1, K+1, &mem);\n  walloc1d(&cnt2, K+1, &mem);\n\n  N++;\n\n  s[N-1] = 0;\n\n  t[N-1] = 1;\n  t[N-2] = 0;\n  for(i=N-3;i>=0;i--){\n    if(s[i] < s[i+1] || (s[i]==s[i+1] && t[i+1])) t[i] = 1; else t[i] = 0;\n  }\n  lms[0] = 0;\n  REP(i,1,N){\n    if(t[i] && !t[i-1]) lms[i] = 1; else lms[i] = 0;\n  }\n\n  rep(i,K+1) cnt1[i] = 0;\n  rep(i,N) cnt1[s[i]]++;\n  j = 0;\n  rep(i,K+1){\n    j += cnt1[i];\n    cnt2[i] = j - cnt1[i];\n    cnt1[i] = j;\n  }\n\n  rep(i,K+1) cnt[i] = cnt1[i];\n  for(i=0; i<N; i++) SA[i] = -1;\n  for(i=1; i<N; i++) if(lms[i]) SA[--cnt[s[i]]]=i;\n  \n  rep(i,K+1) cnt[i] = cnt2[i];\n  rep(i,N){\n    j = SA[i]-1;\n    if(j>=0 && !t[j]) SA[cnt[s[j]]++] = j;\n  }\n\n  rep(i,K+1) cnt[i] = cnt1[i];\n  for(i=N-1;i>=0;i--){\n    j = SA[i] - 1;\n    if(j>=0 && t[j]) SA[--cnt[s[j]]] = j;\n  }\n\n  m = 0;\n  rep(i,N) if(lms[SA[i]]) SA[m++] = SA[i];\n  REP(i,m,N) SA[i] = -1;\n  \n  name=0;\n  prev=-1;\n  rep(i,m){\n    pos = SA[i];\n    rep(d,N){\n      if(prev==-1 || s[pos+d]!=s[prev+d] || t[pos+d]!=t[prev+d]){\n        name++;\n        prev=pos;\n        break;\n      } else if(d>0 && (lms[pos+d] || lms[prev+d])){\n        break;\n      }\n    }\n    pos /= 2;\n    SA[m+pos]=name-1;\n  }\n  for(i=N-1, j=N-1; i>=m; i--) if(SA[i]>=0) SA[j--]=SA[i];\n\n  s1 = SA+N-m;\n  if(name<m){\n    SuffixArray(s1, m-1, name-1, SA, NULL, mem);\n  } else {\n    for(i=0; i<m; i++) SA[s1[i]] = i;\n  }\n\n  rep(i,K+1) cnt[i] = cnt1[i];\n  \n  for(i=1, j=0; i<N; i++) if(lms[i]) s1[j++]=i;\n  for(i=0; i<m; i++) SA[i]=s1[SA[i]];\n  for(i=m; i<N; i++) SA[i]=-1;\n  for(i=m-1; i>=0; i--) {\n    j=SA[i]; SA[i]=-1;\n    SA[--cnt[s[j]]]=j;\n  }\n\n  rep(i,N){\n    j = SA[i]-1;\n    if(j>=0 && !t[j]) SA[cnt2[s[j]]++] = j;\n  }\n\n  for(i=N-1;i>=0;i--){\n    j = SA[i] - 1;\n    if(j>=0 && t[j]) SA[--cnt1[s[j]]] = j;\n  }\n\n  if(LCP != NULL){\n    cnt = (int*)t;\n    d = 0;\n    rep(i,N) cnt[SA[i]] = i;\n    rep(i,N){\n      if(cnt[i]){\n        for(j=SA[cnt[i]-1]; j+d<N-1&&i+d<N-1&&s[j+d]==s[i+d];d++);\n        LCP[cnt[i]]=d;\n      } else {\n        LCP[cnt[i]] = -1;\n      }\n      if(d>0) d--;\n    }\n  }\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "unionFind";
      string c = "struct unionFind{\n  int *d, N, M;\n  unionFind(){}\n  unionFind(const char mode, const int n, void **mem = &wmem){\n    assert(mode == 'm' || mode == 'w');\n    if(mode == 'm') malloc(n);\n    if(mode == 'w') walloc(n, mem);\n  }\n  unionFind(const char mode, const int n, const int fg, void **mem = &wmem){\n    assert(mode == 'm' || mode == 'w');\n    if(mode == 'm') malloc(n, fg);\n    if(mode == 'w') walloc(n, fg, mem);\n  }\n  inline void malloc(const int n){\n    d = (int*)std::malloc(n*sizeof(int));\n    M = n;\n  }\n  inline void malloc(const int n, const int fg){\n    d = (int*)std::malloc(n*sizeof(int));\n    M = n;\n    if(fg) init(n);\n  }\n  inline void free(void){\n    std::free(d);\n  }\n  inline void walloc(const int n, void **mem=&wmem){\n    walloc1d(&d, n, mem);\n    M = n;\n  }\n  inline void walloc(const int n, const int fg, void **mem=&wmem){\n    walloc1d(&d, n, mem);\n    M = n;\n    if(fg) init(n);\n  }\n  inline void init(const int n){\n    int i;\n    N = n;\n    rep(i,n) d[i] = -1;\n  }\n  inline void init(void){\n    init(M);\n  }\n  inline int get(int a){\n    int t = a, k;\n    while(d[t]>=0) t=d[t];\n    while(d[a]>=0) k=d[a], d[a]=t, a=k;\n    return a;\n  }\n  inline int connect(int a, int b){\n    if(d[a]>=0) a=get(a);\n    if(d[b]>=0) b=get(b);\n    if(a==b) return 0;\n    if(d[a] < d[b]) d[a] += d[b], d[b] = a;\n    else            d[b] += d[a], d[a] = b;\n    return 1;\n  }\n  inline int operator()(int a){\n    return get(a);\n  }\n  inline int operator()(int a, int b){\n    return connect(a,b);\n  }\n  inline int& operator[](const int a){\n    return d[a];\n  }\n  inline int size(int a){\n    a = get(a);\n    return -d[a];\n  }\n  inline int sizeList(int res[]){\n    int i, sz=0;\n    rep(i,N) if(d[i]<0) res[sz++] = -d[i];\n    return sz;\n  }\n  inline int comp(int res[], void *mem = wmem){\n    int i, sz=0;\n    int *cnt;\n    walloc1d(&cnt, N, &mem);\n    rep(i,N) cnt[i] = 0;\n    rep(i,N) cnt[get(i)] = 1;\n    rep(i,N) if(cnt[i]) cnt[i] = sz++;\n    rep(i,N) res[i] = cnt[get(i)];\n    return sz;\n  }\n};\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
    }
    {
      string n = "weightedUnionFind";
      string c = "template<class T>\nstruct weightedUnionFind{\n  int *d, N, M;\n  T *w;\n  weightedUnionFind(){}\n  weightedUnionFind(const char mode, const int n, void **mem = &wmem){\n    assert(mode == 'm' || mode == 'w');\n    if(mode == 'm') malloc(n);\n    if(mode == 'w') walloc(n, mem);\n  }\n  weightedUnionFind(const char mode, const int n, const int fg, void **mem = &wmem){\n    assert(mode == 'm' || mode == 'w');\n    if(mode == 'm') malloc(n, fg);\n    if(mode == 'w') walloc(n, fg, mem);\n  }\n  inline void malloc(const int n){\n    d = (int*)std::malloc(n*sizeof(int));\n    w = (T*)std::malloc(n*sizeof(T));\n    M = n;\n  }\n  inline void malloc(const int n, const int fg){\n    d = (int*)std::malloc(n*sizeof(int));\n    w = (T*)std::malloc(n*sizeof(T));\n    M = n;\n    if(fg) init(n);\n  }\n  inline void free(void){\n    std::free(d);\n  }\n  inline void walloc(const int n, void **mem=&wmem){\n    walloc1d(&d, n, mem);\n    walloc1d(&w, n, mem);\n    M = n;\n  }\n  inline void walloc(const int n, const int fg, void **mem=&wmem){\n    walloc1d(&d, n, mem);\n    walloc1d(&w, n, mem);\n    M = n;\n    if(fg) init(n);\n  }\n  inline void init(const int n){\n    int i;\n    N = n;\n    rep(i,n) d[i] = -1;\n    rep(i,n) w[i] = 0;\n  }\n  inline void init(void){\n    init(M);\n  }\n  inline int get(int a){\n    int t = a, k;\n    T s, m;\n    s = 0;\n    while(d[t]>=0) s+=w[t], t=d[t];\n    while(d[a]>=0) k=d[a], m=w[a], d[a]=t, w[a]=s, a=k, s-=m;\n    return a;\n  }\n  inline int connect(int a, int b, T c){\n    if(d[a]>=0) get(a), c += w[a], a = get(a);\n    if(d[b]>=0) get(b), c -= w[b], b = get(b);\n    if(a==b) return 0;\n    if(d[a] < d[b]) d[a] += d[b], d[b] = a, w[b] = c;\n    else            d[b] += d[a], d[a] = b, w[a] = -c;\n    return 1;\n  }\n  inline T diff(int a, int b){\n    get(a);\n    get(b);\n    return w[b]-w[a];\n  }\n  inline int operator()(int a){\n    return get(a);\n  }\n  inline int operator()(int a, int b, T c){\n    return connect(a,b,c);\n  }\n  inline int& operator[](const int a){\n    return d[a];\n  }\n  inline int size(int a){\n    a = get(a);\n    return -d[a];\n  }\n  inline int sizeList(int res[]){\n    int i, sz=0;\n    rep(i,N) if(d[i]<0) res[sz++] = -d[i];\n    return sz;\n  }\n  inline int comp(int res[], void *mem = wmem){\n    int i, sz=0;\n    int *cnt;\n    walloc1d(&cnt, N, &mem);\n    rep(i,N) cnt[i] = 0;\n    rep(i,N) cnt[get(i)] = 1;\n    rep(i,N) if(cnt[i]) cnt[i] = sz++;\n    rep(i,N) res[i] = cnt[get(i)];\n    return sz;\n  }\n};\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
    }
    {
      string n = "rollbackUnionFind";
      string c = "struct rollbackUnionFind{\n  int *d, N, M;\n  int *h1, *h2, sz, snapsz;\n  rollbackUnionFind(){}\n  rollbackUnionFind(const char mode, const int n, void **mem = &wmem){\n    assert(mode == 'm' || mode == 'w');\n    if(mode == 'm') malloc(n);\n    if(mode == 'w') walloc(n, mem);\n  }\n  rollbackUnionFind(const char mode, const int n, const int fg, void **mem = &wmem){\n    assert(mode == 'm' || mode == 'w');\n    if(mode == 'm') malloc(n, fg);\n    if(mode == 'w') walloc(n, fg, mem);\n  }\n  inline void malloc(const int n){\n    d = (int*)std::malloc(n*sizeof(int));\n    h1 = (int*)std::malloc(2*n*sizeof(int));\n    h2 = (int*)std::malloc(2*n*sizeof(int));\n    M = n;\n  }\n  inline void malloc(const int n, const int fg){\n    d = (int*)std::malloc(n*sizeof(int));\n    h1 = (int*)std::malloc(2*n*sizeof(int));\n    h2 = (int*)std::malloc(2*n*sizeof(int));\n    M = n;\n    if(fg) init(n);\n  }\n  inline void free(void){\n    std::free(d);\n    std::free(h1);\n    std::free(h2);\n  }\n  inline void walloc(const int n, void **mem=&wmem){\n    walloc1d(&d, n, mem);\n    walloc1d(&h1, 2*n, mem);\n    walloc1d(&h2, 2*n, mem);\n    M = n;\n  }\n  inline void walloc(const int n, const int fg, void **mem=&wmem){\n    walloc1d(&d, n, mem);\n    walloc1d(&h1, 2*n, mem);\n    walloc1d(&h2, 2*n, mem);\n    M = n;\n    if(fg) init(n);\n  }\n  inline void init(const int n){\n    int i;\n    N = n;\n    rep(i,n) d[i] = -1;\n    sz = 0;\n    snapsz = 0;\n  }\n  inline void init(void){\n    init(M);\n  }\n  inline int get(int a){\n    while(d[a]>=0) a=d[a];\n    return a;\n  }\n  inline int connect(int a, int b){\n    if(d[a]>=0) a=get(a);\n    if(d[b]>=0) b=get(b);\n    if(a==b) return 0;\n    h1[sz] = a; h2[sz] = d[a]; sz++;\n    h1[sz] = b; h2[sz] = d[b]; sz++;\n    if(d[a] < d[b]) d[a] += d[b], d[b] = a;\n    else            d[b] += d[a], d[a] = b;\n    return 1;\n  }\n  inline int operator()(int a){\n    return get(a);\n  }\n  inline int operator()(int a, int b){\n    return connect(a,b);\n  }\n  inline int& operator[](const int a){\n    return d[a];\n  }\n  inline int size(int a){\n    a = get(a);\n    return -d[a];\n  }\n  inline int sizeList(int res[]){\n    int i, sz=0;\n    rep(i,N) if(d[i]<0) res[sz++] = -d[i];\n    return sz;\n  }\n  inline int comp(int res[], void *mem = wmem){\n    int i, sz=0;\n    int *cnt;\n    walloc1d(&cnt, N, &mem);\n    rep(i,N) cnt[i] = 0;\n    rep(i,N) cnt[get(i)] = 1;\n    rep(i,N) if(cnt[i]) cnt[i] = sz++;\n    rep(i,N) res[i] = cnt[get(i)];\n    return sz;\n  }\n  inline int undo(void){\n    if(sz==0) return 0;\n    sz--; d[h1[sz]] = h2[sz];\n    sz--; d[h1[sz]] = h2[sz];\n    if(snapsz > sz) snapsz = sz;\n    return 1;\n  }\n  inline void snapshot(void){\n    snapsz = sz;\n  }\n  inline void rollback(void){\n    while(snapsz < sz){\n      sz--; d[h1[sz]] = h2[sz];\n      sz--; d[h1[sz]] = h2[sz];\n    }\n  }\n};\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
    }

    {
      string n = "fibonacci_mod";
      string c = "int fibonacci_mod_L(ull n, int md){\n  ull a=1, b=0, ma=0, mb=1, ta, tb;\n  while(n){\n    if(n%2){\n      ta = a*ma + b*mb;\n      tb = a*mb + b*(ma+mb);\n      a = ta % md;\n      b = tb % md;\n    }\n    ta = ma*ma + mb*mb;\n    tb = (ma*2 + mb)*mb;\n    ma = ta % md;\n    mb = tb % md;\n    n/=2;\n  }\n  return b;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "Fibonacci_mod";
      string c = "int Fibonacci_mod(ull n){\n  ull a=1, b=0, ma=0, mb=1, ta, tb;\n  while(n){\n    if(n%2){\n      ta = a*ma + b*mb;\n      tb = a*mb + b*(ma+mb);\n      a = ta % MD;\n      b = tb % MD;\n    }\n    ta = ma*ma + mb*mb;\n    tb = (ma*2 + mb)*mb;\n    ma = ta % MD;\n    mb = tb % MD;\n    n/=2;\n  }\n  return b;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"define_MD");
      need[n] = d;
    }
    {
      string n = "Fib_mod";
      string c = "ull Fib_mod_L_ma[64], Fib_mod_L_mb[64];\nvoid Fib_mod_init(){\n  int i;\n  Fib_mod_L_ma[0] = 0;\n  Fib_mod_L_mb[0] = 1;\n  rep(i,1,64){\n    Fib_mod_L_ma[i] = (Fib_mod_L_ma[i-1]*Fib_mod_L_ma[i-1] + Fib_mod_L_mb[i-1]*Fib_mod_L_mb[i-1])%MD;\n    Fib_mod_L_mb[i] = ((Fib_mod_L_ma[i-1]*2 + Fib_mod_L_mb[i-1]) * Fib_mod_L_mb[i-1])%MD;\n  }\n}\nint Fib_mod(ull n){\n  int i = 0;\n  ull a=1, b=0, ta, tb;\n  while(n){\n    if(n%2){\n      ta = a*Fib_mod_L_ma[i] + b*Fib_mod_L_mb[i];\n      tb = a*Fib_mod_L_mb[i] + b*(Fib_mod_L_ma[i]+Fib_mod_L_mb[i]);\n      a = ta % MD;\n      b = tb % MD;\n    }\n    i++;\n    n/=2;\n  }\n  return b;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"define_MD");
      d.push_back((string)"Fib_mod_init");
      need[n] = d;
    }
    {
      string n = "Fib_mod_init";
      string c = "{Fib_mod_init();}\n";
      string p = "main_first";
      name.push_back(n); func[n] = c; place[n] = p;
    }


    {
      string n = "get_inv_mod";
      string c = "int get_inv_mod(ll a, int md){\n  ll t = a, s = md, u = 1, v = 0, e;\n  while(s){\n    e = t/s;\n    t -= e*s;\n    u -= e*v;\n    swap(t,s);\n    swap(u,v);\n  }\n  if(u<0) u += md;\n  return u;\n}\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "extendedEuclid";
      string c = "void extendedEuclid(ll a, ll b, ll &x, ll &y, ll &z){\n  ll d=1, e=0, f, g=0, h=1, i, m, q, aa = a, bb = b;\n  if(a==0 && b==0){\n    x = y = z = 0;\n    return;\n  }\n  while(b){\n    q=a/b;\n    m=a%b;\n    f=d-q*e;\n    i=g-q*h;\n    a=b;\n    b=m;\n    d=e;\n    e=f;\n    g=h;\n    h=i;\n  }\n  if(a < 0){\n    d = -d;\n    g = -g;\n    a = -a;\n  }\n  if(a != 1){\n    aa /= a;\n    bb /= a;\n  }\n  if(d < 0 && bb > 0){\n    d += bb;\n    g -= aa;\n  }\n  if(d < 0 && bb < 0){\n    d -= bb;\n    g += aa;\n  }\n  x=d;\n  y=g;\n  z=a;\n}\n\nvoid extendedEuclid(ll a, ll b, ll &x, ll &y){\n  ll d=1, e=0, f, g=0, h=1, i, m, q, aa = a, bb = b;\n  if(a==0 && b==0){\n    x = y = 0;\n    return;\n  }\n  while(b){\n    q=a/b;\n    m=a%b;\n    f=d-q*e;\n    i=g-q*h;\n    a=b;\n    b=m;\n    d=e;\n    e=f;\n    g=h;\n    h=i;\n  }\n  if(a < 0){\n    d = -d;\n    g = -g;\n    a = -a;\n  }\n  if(a != 1){\n    aa /= a;\n    bb /= a;\n  }\n  if(d < 0 && bb > 0){\n    d += bb;\n    g -= aa;\n  }\n  if(d < 0 && bb < 0){\n    d -= bb;\n    g += aa;\n  }\n  x=d;\n  y=g;\n}\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "chineseRemainder";
      string c = "ll chineseRemainder(int sz, ll val[], ll md[]){\n  int i;\n  ll r = 0, m = 1, rn, mn, a, b, c, u, v;\n\n  rep(i,sz){\n    rn = val[i];\n    mn = md[i];\n    if(rn < 0 || rn >= mn) rn %= mn;\n    if(rn < 0) rn += mn;\n    if(m < mn){\n      swap(r, rn);\n      swap(m, mn);\n    }\n    if(m % mn == 0){\n      if(r % mn != rn) return -1;\n      continue;\n    }\n\n    extendedEuclid(m, mn, a, b, c);\n    u = mn / c;\n    if((rn-r)%c) return -1;\n    v = (rn - r) / c % u * a % u;\n    r += v * m;\n    m *= mn / c;\n    if(r < 0) r += m;\n  }\n\n  return r;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"extendedEuclid");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "Cmod2";
      string c = "template<class T>\ninline int Cmod2(T n, T k){\n  if((n&k)==k) return 1;\n  return 0;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "Cmod5";
      string c = "template<class T>\ninline int Cmod5(T n, T k){\n  static int c[5][5];\n  if(c[0][0]==0){\n    int i, j;\n    rep(i,5) c[i][0] = 1;\n    rep(i,1,5) rep(j,1,5) c[i][j] = (c[i-1][j-1] + c[i-1][j]) % 5;\n  }\n  if(k < 0 || k > n) return 0;\n  int res = 1;\n  while(n > 0 && res > 0){\n    res = (res * c[n%5][k%5]) % 5;\n    n /= 5;\n    k /= 5;\n  }\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "Cmod10";
      string c = "template<class T>\ninline int Cmod10(T n, T k){\n  int a, b;\n  if(k < 0 || k > n) return 0;\n  a = Cmod2(n,k);\n  b = Cmod5(n,k);\n  if(b%2 != a) b += 5;\n  return b;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Cmod2");
      d.push_back((string)"Cmod5");
      need[n] = d;
    }

    {
      string n = "Isqrt_f";
      string c = "inline ll Isqrt_f_L(const ll n){\n  ll r = sqrt(n);\n  r = max(r-2, 0);\n  while( (r+1)**2 <= n ) r++;\n  return r;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"max_L");
      d.push_back((string)"pow2");
      need[n] = d;
    }
    {
      string n = "Isqrt_c";
      string c = "inline ll Isqrt_c_L(const ll n){\n  ll r = sqrt(n);\n  r = max(r-2, 0);\n  while( r**2 < n ) r++;\n  return r;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"max_L");
      d.push_back((string)"pow2");
      need[n] = d;
    }
    {
      string n = "Isqrt_s";
      string c = "inline ll Isqrt_s_L(const ll n){\n  ll r = sqrt(n);\n  r = max(r-2, 0);\n  while( (r+1)**2 <= n ) r++;\n  if(r*r!=n) r = -1;\n  return r;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"max_L");
      d.push_back((string)"pow2");
      need[n] = d;
    }
    {
      string n = "Icbrt_f";
      string c = "inline ll Icbrt_f_L(const ll n){\n  ll r = pow(n, 1.0/3);\n  r = max(r-2, 0);\n  while( (r+1)**3 <= n ) r++;\n  return r;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"max_L");
      d.push_back((string)"pow3");
      need[n] = d;
    }
    {
      string n = "Icbrt_c";
      string c = "inline ll Icbrt_c_L(const ll n){\n  ll r = pow(n, 1.0/3);\n  r = max(r-2, 0);\n  while( r**3 < n ) r++;\n  return r;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"max_L");
      d.push_back((string)"pow3");
      need[n] = d;
    }
    {
      string n = "Icbrt_s";
      string c = "inline ll Icbrt_s_L(const ll n){\n  ll r = pow(n, 1.0/3);\n  r = max(r-2, 0);\n  while( (r+1)**3 <= n ) r++;\n  if(r**3!=n) r = -1;\n  return r;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"max_L");
      d.push_back((string)"pow3");
      need[n] = d;
    }
    {
      string n = "Iroot_f";
      string c = "inline ll Iroot_f_L(const ll N, const ll k){\n  ll res;\n  if(N<=0) return 0;\n  if(k==1) return N;\n  if(k==2) return Isqrt_f(N);\n  if(k==3) return Icbrt_f(N);\n  if(k>=64) return 1;\n  if(k>=40){\n    if((1LL<<k) <= N) return 2;\n    return 1;\n  }\n  res = pow(N, 1.0 / k);\n  if(res > 0) res--;\n  if( ((double)res+1)**k < 9.2233e18 && (res+1)**k <= N ) res++;\n  if( ((double)res+1)**k < 9.2233e18 && (res+1)**k <= N ) res++;\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Isqrt_f");
      d.push_back((string)"Icbrt_f");
      d.push_back((string)"pow");
      need[n] = d;
    }
    {
      string n = "Iroot_c";
      string c = "inline ll Iroot_c_L(const ll N, const ll k){\n  ll res;\n  if(N<=0) return 0;\n  if(N==1) return 1;\n  if(k==1) return N;\n  if(k==2) return Isqrt_c(N);\n  if(k==3) return Icbrt_c(N);\n  if(k>=64) return 1;\n  if(k>=40){\n    if((1LL<<k) < N) return 3;\n    return 2;\n  }\n  res = pow(N, 1.0 / k);\n  if( ((double)res)**k < 9.2233e18 && res**k < N ) res++;\n  if( ((double)res)**k < 9.2233e18 && res**k < N ) res++;\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Isqrt_c");
      d.push_back((string)"Icbrt_c");
      d.push_back((string)"pow");
      need[n] = d;
    }
    {
      string n = "Iroot_s";
      string c = "inline ll Iroot_s_L(const ll N, const ll k){\n  ll res;\n  if(N<=0) return 0;\n  if(N==1) return 1;\n  if(k==1) return N;\n  if(k==2) return Isqrt_s(N);\n  if(k==3) return Icbrt_s(N);\n  if(k>=64) return 1;\n  if(k>=40){\n    if((1LL<<k) == N) return 2;\n    return -1;\n  }\n  res = pow(N, 1.0 / k);\n  if( ((double)res)**k < 9.2233e18 && res**k < N ) res++;\n  if( ((double)res)**k > 9.2233e18 || res**k != N) return -1;\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Isqrt_s");
      d.push_back((string)"Icbrt_s");
      d.push_back((string)"pow");
      need[n] = d;
    }

    {
      string n = "Unique";
      string c = "template<class T>\nvoid Unique(int &N, T A[], int sorted=0, void *mwm = wmem){\n  int i, k;\n  if(!sorted) sortA(N, A);\n  k = 0;\n  rep(i,N) if(k==0 || A[k-1]!=A[i]) A[k++] = A[i];\n  N = k;\n}\ntemplate<class T, class S>\nvoid Unique(int &N, T A[], S B[], int sorted=0, void *mem = wmem){\n  int i, k = 0;\n  if(!sorted) sortA(N, A, B, mem);\n  rep(i,N){\n    if(!k || A[k-1]!=A[i]){\n      A[k] = A[i];\n      B[k] = B[i];\n      k++;\n    } else {\n      B[k-1] += B[i];\n    }\n  }\n  N=k;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"sortA");
      need[n] = d;
    }

    {
      string n = "coordcomp_1";
      string c = "template<class T>\nint coordcomp_L(int n, T arr[], void *mem = wmem){\n  int i, k = 0;\n  T *r;\n  int *ind;\n  walloc1d(&r, n, &mem);\n  walloc1d(&ind, n, &mem);\n  rep(i,n) r[i] = arr[i];\n  rep(i,n) ind[i] = i;\n  sortA(n, r, ind, mem);\n  rep(i,n){\n    if(i && r[i] != r[i-1]) k++;\n    arr[ind[i]] = k;\n  }\n  return k+1;\n}\ntemplate<class T, class S>\nint coordcomp_L(int n, T arr[], S res[], void *mem = wmem){\n  int i, k = 0;\n  T *r;\n  int *ind;\n  walloc1d(&r, n, &mem);\n  walloc1d(&ind, n, &mem);\n  rep(i,n) r[i] = arr[i];\n  rep(i,n) ind[i] = i;\n  sortA(n, r, ind, mem);\n  rep(i,n){\n    if(i && r[i] != r[i-1]) k++;\n    res[ind[i]] = k;\n  }\n  return k+1;\n}\ntemplate<class T, class S>\nint coordcomp_L(int n, const T arr[], S res[], void *mem = wmem){\n  int i, k = 0;\n  T *r;\n  int *ind;\n  walloc1d(&r, n, &mem);\n  walloc1d(&ind, n, &mem);\n  rep(i,n) r[i] = arr[i];\n  rep(i,n) ind[i] = i;\n  sortA(n, r, ind, mem);\n  rep(i,n){\n    if(i && r[i] != r[i-1]) k++;\n    res[ind[i]] = k;\n  }\n  return k+1;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"sortA");
      need[n] = d;
    }
    {
      string n = "coordcomp_2";
      string c = "template<class T>\nint coordcomp_L(int n1, T arr1[], int n2, T arr2[], void *mem = wmem){\n  int i, k = 0;\n  T *r;\n  int *ind;\n  walloc1d(&r, n1+n2, &mem);\n  walloc1d(&ind, n1+n2, &mem);\n  rep(i,n1) r[i] = arr1[i];\n  rep(i,n2) r[n1+i] = arr2[i];\n  rep(i,n1+n2) ind[i] = i;\n  sortA(n1+n2, r, ind, mem);\n  rep(i,n1+n2){\n    if(i && r[i] != r[i-1]) k++;\n    if(ind[i] < n1){\n      arr1[ind[i]] = k;\n    } else {\n      arr2[ind[i]-n1] = k;\n    }\n  }\n  return k+1;\n}\ntemplate<class T, class S1, class S2>\nint coordcomp_L(int n1, const T arr1[], int n2, const T arr2[], S1 res1[], S2 res2[], void *mem = wmem){\n  int i, k = 0;\n  T *r;\n  int *ind;\n  walloc1d(&r, n1+n2, &mem);\n  walloc1d(&ind, n1+n2, &mem);\n  rep(i,n1) r[i] = arr1[i];\n  rep(i,n2) r[n1+i] = arr2[i];\n  rep(i,n1+n2) ind[i] = i;\n  sortA(n1+n2, r, ind, mem);\n  rep(i,n1+n2){\n    if(i && r[i] != r[i-1]) k++;\n    if(ind[i] < n1){\n      res1[ind[i]] = k;\n    } else {\n      res2[ind[i]-n1] = k;\n    }\n  }\n  return k+1;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"sortA");
      need[n] = d;
    }
    {
      string n = "coordcomp_3";
      string c = "template<class T>\nint coordcomp_L(int n1, T arr1[], int n2, T arr2[], int n3, T arr3[], void *mem = wmem){\n  int i, k = 0, nn = n1 + n2 + n3;\n  T *r;\n  int *ind;\n  walloc1d(&r, nn, &mem);\n  walloc1d(&ind, nn, &mem);\n  rep(i,n1) r[i] = arr1[i];\n  rep(i,n2) r[n1+i] = arr2[i];\n  rep(i,n3) r[n1+n2+i] = arr3[i];\n  rep(i,nn) ind[i] = i;\n  sortA(nn, r, ind, mem);\n  rep(i,nn){\n    if(i && r[i] != r[i-1]) k++;\n    if(ind[i] < n1){\n      arr1[ind[i]] = k;\n    } else if(ind[i] < n1+n2){\n      arr2[ind[i]-n1] = k;\n    } else {\n      arr3[ind[i]-n1-n2] = k;\n    }\n  }\n  return k+1;\n}\ntemplate<class T, class S1, class S2, class S3>\nint coordcomp_L(int n1, const T arr1[], int n2, const T arr2[], int n3, const T arr3[], S1 res1[], S2 res2[], S3 res3[], void *mem = wmem){\n  int i, k = 0, nn = n1 + n2 + n3;\n  T *r;\n  int *ind;\n  walloc1d(&r, nn, &mem);\n  walloc1d(&ind, nn, &mem);\n  rep(i,n1) r[i] = arr1[i];\n  rep(i,n2) r[n1+i] = arr2[i];\n  rep(i,n3) r[n1+n2+i] = arr3[i];\n  rep(i,nn) ind[i] = i;\n  sortA(nn, r, ind, mem);\n  rep(i,nn){\n    if(i && r[i] != r[i-1]) k++;\n    if(ind[i] < n1){\n      res1[ind[i]] = k;\n    } else if(ind[i] < n1+n2){\n      res2[ind[i]-n1] = k;\n    } else {\n      res3[ind[i]-n1-n2] = k;\n    }\n  }\n  return k+1;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"sortA");
      need[n] = d;
    }


    {
      string n = "Heap";
      string c = "template<class T>\nstruct Heap {\n  int size;\n  T *val;\n\n  void malloc(const int N){\n    val = (T*) std::malloc(N*sizeof(T));\n    size = 0;\n  }\n\n  void walloc(const int N, void **mem = &wmem){\n    walloc1d(&val, N, mem);\n    size = 0;\n  }\n\n  void free(){\n    std::free(val);\n  }\n\n  void init(){\n    size = 0;\n  }\n\n  void up(){\n    int n = size - 1, m;\n    while(n){\n      m = (n-1) / 2;\n      if(val[m] <= val[n]) break;\n      swap(val[m], val[n]);\n      n = m;\n    }\n  }\n  \n  void down(){\n    int n = 0, m;\n    for(;;){\n      m=2n+1;\n      if(m>=size) break;\n      if(m+1<size && val[m] > val[m+1]) m++;\n      if(val[m] >= val[n]) break;\n      swap(val[m], val[n]);\n      n = m;\n    }\n  }\n\n  T top(){\n    return val[0];\n  }\n\n  T pop(){\n    T res = val[0];\n    size--;\n    if(size > 0) val[0] = val[size], down();\n    return res;\n  }\n\n  T push(const T x){\n    val[size++] = x;\n    up();\n    return x;\n  }\n};\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "Heap_max";
      string c = "template<class T>\nstruct Heap_max {\n  int size;\n  T *val;\n\n  void malloc(const int N){\n    val = (T*) std::malloc(N*sizeof(T));\n    size = 0;\n  }\n\n  void walloc(const int N, void **mem = &wmem){\n    walloc1d(&val, N, mem);\n    size = 0;\n  }\n\n  void free(){\n    std::free(val);\n  }\n\n  void init(){\n    size = 0;\n  }\n\n  void up(){\n    int n = size - 1, m;\n    while(n){\n      m = (n-1) / 2;\n      if(val[m] >= val[n]) break;\n      swap(val[m], val[n]);\n      n = m;\n    }\n  }\n  \n  void down(){\n    int n = 0, m;\n    for(;;){\n      m=2n+1;\n      if(m>=size) break;\n      if(m+1<size && val[m] < val[m+1]) m++;\n      if(val[m] <= val[n]) break;\n      swap(val[m], val[n]);\n      n=m;\n    }\n  }\n\n  T top(){\n    return val[0];\n  }\n\n  T pop(){\n    T res = val[0];\n    size--;\n    if(size > 0) val[0] = val[size], down();\n    return res;\n  }\n\n  T push(const T x){\n    val[size++] = x;\n    up();\n    return x;\n  }\n};\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }
    
    {
      string n = "LHeap";
      string c = "template <class T>\nstruct LHeap {\n  int *hp, *place, size;\n  T *val;\n\n  void malloc(int N){\n    hp = (int*)std::malloc(N*sizeof(int));\n    place=(int*)std::malloc(N*sizeof(int));\n    val=(T*)std::malloc(N*sizeof(T));\n  }\n  void malloc(int N, int ini){\n    hp = (int*)std::malloc(N*sizeof(int));\n    place=(int*)std::malloc(N*sizeof(int));\n    val=(T*)std::malloc(N*sizeof(T));\n    if(ini) init(N);\n  }\n  void walloc(int N, void **mem=&wmem){\n    walloc1d(&hp, N, mem);\n    walloc1d(&place, N, mem);\n    walloc1d(&val, N, mem);\n  }\n  void walloc(int N, int ini, void **mem=&wmem){\n    walloc1d(&hp, N, mem);\n    walloc1d(&place, N, mem);\n    walloc1d(&val, N, mem);\n    if(ini) init(N);\n  }\n  void free(){\n    std::free(hp);\n    std::free(place);\n    std::free(val);\n  }\n  void init(int N){\n    int i;\n    size=0;\n    rep(i,N) place[i]=-1;\n  }\n  void up(int n){\n    int m;\n    while(n){\n      m=(n-1)/2;\n      if(val[hp[m]]<=val[hp[n]])break;\n      swap(hp[m],hp[n]);\n      swap(place[hp[m]],place[hp[n]]);\n      n=m;\n    }\n  }\n  void down(int n){\n    int m;\n    for(;;){\n      m=2*n+1;\n      if(m>=size)break;\n      if(m+1<size&&val[hp[m]]>val[hp[m+1]])m++;\n      if(val[hp[m]]>=val[hp[n]])break;\n      swap(hp[m],hp[n]);\n      swap(place[hp[m]],place[hp[n]]);\n      n=m;\n    }\n  }\n  void change(int n, T v){\n    T f = val[n];\n    val[n] = v;\n    if(place[n]==-1){\n      place[n] = size;\n      hp[size++] = n;\n      up(place[n]);\n    } else {\n      if(f < v) down(place[n]);\n      else if(f > v) up(place[n]);\n    }\n  }\n  int pop(void){\n    int res = hp[0];\n    place[res] = -1;\n    size--;\n    if(size){\n      hp[0]=hp[size];\n      place[hp[0]]=0;\n      down(0);\n    }\n    return res;\n  }\n};\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "DijkstraHeap";
      string c = "template <class T>\nstruct DijkstraHeap {\n  int *hp, *place, size;\n  char *visited; T *val;\n  void malloc(int N){\n    hp = (int*)std::malloc(N*sizeof(int));\n    place = (int*)std::malloc(N*sizeof(int));\n    visited = (char*)std::malloc(N*sizeof(char));\n    val = (T*)std::malloc(N*sizeof(T));\n  }\n  void free(){\n    std::free(hp);\n    std::free(place);\n    std::free(visited);\n    std::free(val);\n  }\n  void walloc(int N, void **mem=&wmem){\n    walloc1d(&hp, N, mem);\n    walloc1d(&place, N, mem);\n    walloc1d(&visited, N, mem);\n    walloc1d(&val, N, mem);\n  }\n  void malloc(int N, int init_fg){\n    malloc(N);\n    if(init_fg) init(N);\n  }\n  void walloc(int N, int init_fg, void **mem=&wmem){\n    walloc(N,mem);\n    if(init_fg) init(N);\n  }\n  void init(int N){\n    int i;\n    size = 0;\n    rep(i,N) place[i]=-1;\n    rep(i,N) visited[i]=0;\n  }\n  void up(int n){\n    int m;\n    while(n){\n      m=(n-1)/2;\n      if(val[hp[m]]<=val[hp[n]])break;\n      swap(hp[m],hp[n]);\n      swap(place[hp[m]],place[hp[n]]);\n      n=m;\n    }\n  }\n  void down(int n){\n    int m;\n    for(;;){\n      m=2*n+1;\n      if(m>=size)break;\n      if(m+1<size&&val[hp[m]]>val[hp[m+1]])m++;\n      if(val[hp[m]]>=val[hp[n]]) break;\n      swap(hp[m],hp[n]);\n      swap(place[hp[m]],place[hp[n]]);\n      n=m;\n    }\n  }\n  void change(int n, T v){\n    if(visited[n]||(place[n]>=0&&val[n]<=v))return;\n    val[n]=v;\n    if(place[n]==-1)place[n]=size,hp[size++]=n,up(place[n]);\n    else up(place[n]);\n  }\n  int pop(void){\n    int res=hp[0];\n    place[res]=-1;\n    size--;\n    if(size)hp[0]=hp[size],place[hp[0]]=0,down(0);\n    visited[res]=1;\n    return res;\n  }\n};\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
    }

    {
      string n = "fenwick_header";
      string c = "template<class T>\nstruct fenwick{\n  int size, memory;\n  T *data;\n  void malloc(int mem);\n  void malloc(int mem, int fg);\n  void walloc(int mem, void **workMemory = &wmem);\n  void walloc(int mem, int fg, void **workMemory = &wmem);\n  void free(void);\n  void init(int N);\n  void add(int k, T val);\n  T get(int k);\n  T range(int a, int b);\n  int kth(T k);\n};\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "fenwick";
      string c = "template<class T> void fenwick<T>::malloc(int mem){\n  memory = mem;\n  data = (T*)std::malloc(sizeof(T)*mem);\n}\ntemplate<class T> void fenwick<T>::malloc(int mem, int fg){\n  memory = mem;\n  data = (T*)std::malloc(sizeof(T)*mem);\n  if(fg) init(mem);\n}\ntemplate<class T> void fenwick<T>::walloc(int mem, void **workMemory /* = &wmem*/){\n  memory = mem;\n  walloc1d(&data, mem, workMemory);\n}\ntemplate<class T> void fenwick<T>::walloc(int mem, int fg, void **workMemory /* = &wmem*/){\n  memory = mem;\n  walloc1d(&data, mem, workMemory);\n  if(fg) init(mem);\n}\ntemplate<class T> void fenwick<T>::free(void){\n  memory = 0;\n  free(data);\n}\ntemplate<class T> void fenwick<T>::init(int N){\n  size = N;\n  memset(data,0,sizeof(T)*N);\n}\ntemplate<class T> void fenwick<T>::add(int k, T val){\n  while(k < size) data[k] += val, k |= k+1;\n}\ntemplate<class T> T fenwick<T>::get(int k){\n  T res = 0;\n  while(k>=0) res += data[k], k = (k&(k+1))-1;\n  return res;\n}\ntemplate<class T> T fenwick<T>::range(int a, int b){\n  if(a < 0) a = 0;\n  if(b >= size) b = size - 1;\n  if(b < a) return 0;\n  return get(b) - get(a-1);\n}\ntemplate<class T> int fenwick<T>::kth(T k){\n  int i=0, j=size, c;\n  T v;\n  while(i<j){\n    c = (i+j)/2;\n    v = get(c);\n    if(v <= k) i=c+1; else j=c;\n  }\n  return i==size?-1:i;\n}\n";
      string p = "last";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"fenwick_header");
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
    }
    {
      string n = "fenwick_xor_header";
      string c = "template<class T>\nstruct fenwick_xor{\n  int size, memory;\n  T *data;\n  void malloc(int mem);\n  void malloc(int mem, int fg);\n  void walloc(int mem, void **workMemory = &wmem);\n  void walloc(int mem, int fg, void **workMemory = &wmem);\n  void free(void);\n  void init(int N);\n  void add(int k, T val);\n  T get(int k);\n  T range(int a, int b);\n};\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "fenwick_xor";
      string c = "template<class T> void fenwick_xor<T>::malloc(int mem){\n  memory = mem;\n  data = (T*)std::malloc(sizeof(T)*mem);\n}\ntemplate<class T> void fenwick_xor<T>::malloc(int mem, int fg){\n  memory = mem;\n  data = (T*)std::malloc(sizeof(T)*mem);\n  if(fg) init(mem);\n}\ntemplate<class T> void fenwick_xor<T>::walloc(int mem, void **workMemory /* = &wmem*/){\n  memory = mem;\n  walloc1d(&data, mem, workMemory);\n}\ntemplate<class T> void fenwick_xor<T>::walloc(int mem, int fg, void **workMemory /* = &wmem*/){\n  memory = mem;\n  walloc1d(&data, mem, workMemory);\n  if(fg) init(mem);\n}\ntemplate<class T> void fenwick_xor<T>::free(void){\n  memory = 0;\n  free(data);\n}\ntemplate<class T> void fenwick_xor<T>::init(int N){\n  size = N;\n  memset(data,0,sizeof(T)*N);\n}\ntemplate<class T> void fenwick_xor<T>::add(int k, T val){\n  while(k < size) data[k] ^= val, k |= k+1;\n}\ntemplate<class T> T fenwick_xor<T>::get(int k){\n  T res = 0;\n  while(k>=0) res ^= data[k], k = (k&(k+1))-1;\n  return res;\n}\ntemplate<class T> T fenwick_xor<T>::range(int a, int b){\n  if(a < 0) a = 0;\n  if(b >= size) b = size - 1;\n  if(b < a) return 0;\n  return get(b) ^ get(a-1);\n}\n";
      string p = "last";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"fenwick_xor_header");
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
    }

    {
      string n = "segtree";
      string c = "template<class T>\nstruct segtree{\n  int N, logN;\n  T *sum, *mn;\n  int *mnind;\n  T *fixval; char *fixed;\n  T *addval;\n  void malloc(int maxN, int once = 0){\n    int i;\n    for(i=1;i<maxN;i*=2);\n    sum = new T[2*i];\n    mn = new T[2*i];\n    mnind = new int[2*i];\n    fixval = new T[i];\n    addval = new T[i];\n    fixed = new char[i];\n    if(once) setN(maxN);\n  }\n  void walloc(int maxN, int once = 0, void **mem = &wmem){\n    int i;\n    for(i=1;i<maxN;i*=2);\n    walloc1d(&sum, 2i, mem);\n    walloc1d(&mn, 2i, mem);\n    walloc1d(&mnind, 2i, mem);\n    walloc1d(&fixval, i, mem);\n    walloc1d(&addval, i, mem);\n    walloc1d(&fixed, i, mem);\n    if(once) setN(maxN);\n  }\n  void free(void){\n    delete [] sum;\n    delete [] mn;\n    delete [] mnind;\n    delete [] fixval;\n    delete [] addval;\n    delete [] fixed;\n  }\n  T& operator[](int i){\n    return sum[N+i];\n  }\n \n  void setN(int n, int zerofill = 1, int dobuild = 1){\n    int i;\n    for(i=1,logN=0;i<n;i*=2,logN++);\n    N = i;\n    if(zerofill) rep(i,N) sum[N+i] = 0;\n    if(dobuild) build();\n  }\n \n  void build(void){\n    int i;\n    rep(i,N) mn[N+i] = sum[N+i], mnind[N+i] = i;\n    for(i=N-1;i;i--){\n      sum[i] = sum[2*i] + sum[2*i+1];\n      if(mn[2*i] <= mn[2*i+1]){\n        mn[i] = mn[2*i];\n        mnind[i] = mnind[2*i];\n      } else {\n        mn[i] = mn[2*i+1];\n        mnind[i] = mnind[2*i+1];\n      }\n    }\n    REP(i,1,N) fixed[i] = 0;\n    REP(i,1,N) addval[i] = 0;\n  }\n \n  inline void push_one(int a, int sz, int st){\n    if(fixed[a]){\n      if(sz > 1){\n        fixed[a*2] = fixed[a*2+1] = 1;\n        fixval[a*2] = fixval[a*2+1] = fixval[a];\n        sum[a*2] = sum[a*2+1] = sz * fixval[a];\n        mn[a*2] = mn[a*2+1] = fixval[a];\n        mnind[a*2] = st;\n        mnind[a*2+1] = st + sz;\n      } else {\n        sum[a*2] = sum[a*2+1] = sz * fixval[a];\n        mn[a*2] = mn[a*2+1] = fixval[a];\n        mnind[a*2] = st;\n        mnind[a*2+1] = st + sz;\n      }\n      fixed[a] = 0;\n      addval[a] = 0;\n      return;\n    }\n    if(addval[a] != 0){\n      if(sz > 1){\n        if(fixed[a*2]) fixval[a*2] += addval[a];\n        else           addval[a*2] += addval[a];\n        if(fixed[a*2+1]) fixval[a*2+1] += addval[a];\n        else             addval[a*2+1] += addval[a];\n        sum[a*2] += sz * addval[a];\n        sum[a*2+1] += sz * addval[a];\n        mn[a*2] += addval[a];\n        mn[a*2+1] += addval[a];\n      } else {\n        sum[a*2] += sz * addval[a];\n        sum[a*2+1] += sz * addval[a];\n        mn[a*2] += addval[a];\n        mn[a*2+1] += addval[a];\n      }\n      addval[a] = 0;\n      return;\n    }\n  }\n \n  inline void push(int a){\n    int i, aa = a - N, nd, sz, st;\n    for(i=logN;i;i--){\n      nd = a>>i;\n      sz = 1<<(i-1);\n      st = 2 * sz * (aa>>i);\n      push_one(nd, sz, st);\n    }\n  }\n \n  inline void build(int a){\n    int sz = 1, st = a - N;\n    while(a > 1){\n      if(a%2) st += sz;\n      a /= 2;\n      sz *= 2;\n      if(fixed[a]){\n        sum[a] = sz * fixval[a];\n        mn[a] = fixval[a];\n      } else {\n        sum[a] = sum[a*2] + sum[a*2+1];\n        if(mn[a*2] <= mn[a*2+1]){\n          mn[a] = mn[a*2];\n          mnind[a] = mnind[a*2];\n        } else {\n          mn[a] = mn[a*2+1];\n          mnind[a] = mnind[a*2+1];\n        }\n        if(addval[a] != 0){\n          mn[a] += addval[a];\n          sum[a] += sz * addval[a];\n        }\n      }\n    }\n  }\n \n  inline void change(int a, int b, T val){\n    int sz = 1, aa, bb, st_a = a, st_b = b;\n    if(a >= b) return;\n \n    aa = (a += N);\n    bb = (b += N);\n    push(a); push(b-1);\n \n    if(a%2){\n      sum[a] = mn[a] = val;\n      a++;\n      st_a += sz;\n    }\n    if(b%2){\n      b--;\n      st_b -= sz;\n      sum[b] = mn[b] = val;\n    }\n    a /= 2;\n    b /= 2;\n \n    while(a < b){\n      sz *= 2;\n      if(a%2){\n        fixed[a]=1, fixval[a]=val;\n        sum[a] = sz * val;\n        mn[a] = val;\n        mnind[a] = st_a;\n        a++;\n        st_a += sz;\n      }\n      if(b%2){\n        b--;\n        st_b -= sz;\n        fixed[b]=1, fixval[b]=val;\n        sum[b] = sz * val;\n        mn[b] = val;\n        mnind[b] = st_b;\n      }\n      a /= 2;\n      b /= 2;\n    }\n \n    build(aa);\n    build(bb-1);\n  }\n \n  inline void add(int a, int b, T val){\n    int sz = 1, aa, bb;\n    if(a >= b) return;\n \n    aa = (a += N);\n    bb = (b += N);\n    push(a); push(b-1);\n \n    if(a%2){\n      sum[a] += val;\n      mn[a] += val;\n      a++;\n    }\n    if(b%2){\n      b--;\n      sum[b] += val;\n      mn[b] += val;\n    }\n    a /= 2;\n    b /= 2;\n \n    while(a < b){\n      sz *= 2;\n      if(a%2){\n        if(fixed[a]) fixval[a] += val; else addval[a] += val;\n        sum[a] += sz * val;\n        mn[a] += val;\n        a++;\n      }\n      if(b%2){\n        b--;\n        if(fixed[b]) fixval[b] += val; else addval[b] += val;\n        sum[b] += sz * val;\n        mn[b] += val;\n      }\n      a /= 2;\n      b /= 2;\n    }\n \n    build(aa);\n    build(bb-1);\n  }\n \n  inline pair<T,int> getMin(int a, int b){\n    pair<T,int> res, tmp;\n    int fga = 0, fgb = 0;\n    \n    a += N;\n    b += N;\n    push(a); push(b-1);\n \n    while(a < b){\n      if(a%2){\n        if(fga){\n          res = min(res, make_pair(mn[a], mnind[a]));\n        } else {\n          res = make_pair(mn[a], mnind[a]);\n          fga = 1;\n        }\n        a++;\n      }\n      if(b%2){\n        b--;\n        if(fgb){\n          tmp = min(make_pair(mn[b], mnind[b]), tmp);\n        } else {\n          tmp = make_pair(mn[b], mnind[b]);\n          fgb = 1;\n        }\n      }\n      a /= 2;\n      b /= 2;\n    }\n    if(fga==1 && fgb==0) return res;\n    if(fga==0 && fgb==1) return tmp;\n    if(fga==1 && fgb==1){\n      res = min(res, tmp);\n      return res;\n    }\n    return res;\n  }\n  inline T getMinVal(int a, int b){\n    T res, tmp;\n    int fga = 0, fgb = 0;\n    \n    a += N;\n    b += N;\n    push(a); push(b-1);\n \n    while(a < b){\n      if(a%2){\n        if(fga){\n          res = min(res, mn[a]);\n        } else {\n          res = mn[a];\n          fga = 1;\n        }\n        a++;\n      }\n      if(b%2){\n        b--;\n        if(fgb){\n          tmp = min(mn[b], tmp);\n        } else {\n          tmp = mn[b];\n          fgb = 1;\n        }\n      }\n      a /= 2;\n      b /= 2;\n    }\n    if(fga==1 && fgb==0) return res;\n    if(fga==0 && fgb==1) return tmp;\n    if(fga==1 && fgb==1){\n      res = min(res, tmp);\n      return res;\n    }\n    return res;\n  }\n  inline int getMinInd(int a, int b){\n    return getMin(a,b).second;\n  }\n  inline T getSum(int a, int b){\n    T res, tmp;\n    int fga = 0, fgb = 0;\n    \n    a += N;\n    b += N;\n    push(a); push(b-1);\n \n    while(a < b){\n      if(a%2){\n        if(fga){\n          res = res + sum[a];\n        } else {\n          res = sum[a];\n          fga = 1;\n        }\n        a++;\n      }\n      if(b%2){\n        b--;\n        if(fgb){\n          tmp = sum[b] + tmp;\n        } else {\n          tmp = sum[b];\n          fgb = 1;\n        }\n      }\n      a /= 2;\n      b /= 2;\n    }\n    if(fga==1 && fgb==0) return res;\n    if(fga==0 && fgb==1) return tmp;\n    if(fga==1 && fgb==1){\n      res = res + tmp;\n      return res;\n    }\n    return res;\n  }\n};\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"min_L");
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
    }
    {
      string n = "segtree_Change_Sum";
      string c = "template<class T>\nstruct segtree_Change_Sum{\n  int N, logN;\n  T *sum;\n  T *fixval; char *fixed;\n  void malloc(int maxN, int once = 0){\n    int i;\n    for(i=1;i<maxN;i*=2);\n    sum = new T[2*i];\n    fixval = new T[i];\n    fixed = new char[i];\n    if(once) setN(maxN);\n  }\n  void walloc(int maxN, int once = 0, void **mem = &wmem){\n    int i;\n    for(i=1;i<maxN;i*=2);\n    walloc1d(&sum, 2i, mem);\n    walloc1d(&fixval, i, mem);\n    walloc1d(&fixed, i, mem);\n    if(once) setN(maxN);\n  }\n  void free(void){\n    delete [] sum;\n    delete [] fixval;\n    delete [] fixed;\n  }\n  T& operator[](int i){\n    return sum[N+i];\n  }\n \n  void setN(int n, int zerofill = 1, int dobuild = 1){\n    int i;\n    for(i=1,logN=0;i<n;i*=2,logN++);\n    N = i;\n    if(zerofill) rep(i,N) sum[N+i] = 0;\n    if(dobuild) build();\n  }\n \n  void build(void){\n    int i;\n    for(i=N-1;i;i--){\n      sum[i] = sum[2*i] + sum[2*i+1];\n    }\n    REP(i,1,N) fixed[i] = 0;\n  }\n \n  inline void push_one(int a, int sz, int st){\n    if(fixed[a]){\n      if(sz > 1){\n        fixed[a*2] = fixed[a*2+1] = 1;\n        fixval[a*2] = fixval[a*2+1] = fixval[a];\n        sum[a*2] = sum[a*2+1] = sz * fixval[a];\n      } else {\n        sum[a*2] = sum[a*2+1] = sz * fixval[a];\n      }\n      fixed[a] = 0;\n    }\n  }\n \n  inline void push(int a){\n    int i, aa = a - N, nd, sz, st;\n    for(i=logN;i;i--){\n      nd = a>>i;\n      sz = 1<<(i-1);\n      st = 2 * sz * (aa>>i);\n      push_one(nd, sz, st);\n    }\n  }\n \n  inline void build(int a){\n    int sz = 1, st = a - N;\n    while(a > 1){\n      if(a%2) st += sz;\n      a /= 2;\n      sz *= 2;\n      if(fixed[a]){\n        sum[a] = sz * fixval[a];\n      } else {\n        sum[a] = sum[a*2] + sum[a*2+1];\n      }\n    }\n  }\n \n  inline void change(int a, int b, T val){\n    int sz = 1, aa, bb, st_a = a, st_b = b;\n    if(a >= b) return;\n \n    aa = (a += N);\n    bb = (b += N);\n    push(a); push(b-1);\n \n    if(a%2){\n      sum[a] = val;\n      a++;\n      st_a += sz;\n    }\n    if(b%2){\n      b--;\n      st_b -= sz;\n      sum[b] = val;\n    }\n    a /= 2;\n    b /= 2;\n \n    while(a < b){\n      sz *= 2;\n      if(a%2){\n        fixed[a]=1, fixval[a]=val;\n        sum[a] = sz * val;\n        a++;\n        st_a += sz;\n      }\n      if(b%2){\n        b--;\n        st_b -= sz;\n        fixed[b]=1, fixval[b]=val;\n        sum[b] = sz * val;\n      }\n      a /= 2;\n      b /= 2;\n    }\n \n    build(aa);\n    build(bb-1);\n  }\n  inline T getSum(int a, int b){\n    T res, tmp;\n    int fga = 0, fgb = 0;\n    \n    a += N;\n    b += N;\n    push(a); push(b-1);\n \n    while(a < b){\n      if(a%2){\n        if(fga){\n          res = res + sum[a];\n        } else {\n          res = sum[a];\n          fga = 1;\n        }\n        a++;\n      }\n      if(b%2){\n        b--;\n        if(fgb){\n          tmp = sum[b] + tmp;\n        } else {\n          tmp = sum[b];\n          fgb = 1;\n        }\n      }\n      a /= 2;\n      b /= 2;\n    }\n    if(fga==1 && fgb==0) return res;\n    if(fga==0 && fgb==1) return tmp;\n    if(fga==1 && fgb==1){\n      res = res + tmp;\n      return res;\n    }\n    return res;\n  }\n};\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
    }


    {
      string n = "segtree_Point_Minval";
      string c = "template<class T>\nstruct segtree_Point_Minval{\n  int N, logN;\n  T *mn;\n\n  void malloc(int maxN, int once = 0){\n    int i;\n    for(i=1;i<maxN;i*=2);\n    mn = new T[2*i];\n    if(once) setN(maxN);\n  }\n\n  void walloc(int maxN, int once = 0, void **mem = &wmem){\n    int i;\n    for(i=1;i<maxN;i*=2);\n    walloc1d(&mn, 2i, mem);\n    if(once) setN(maxN);\n  }\n\n  void free(void){\n    delete [] mn;\n  }\n\n  T& operator[](int i){\n    return mn[N+i];\n  }\n \n  void setN(int n, int zerofill = 1, int dobuild = 1){\n    int i;\n    for(i=1,logN=0;i<n;i*=2,logN++);\n    N = i;\n    if(zerofill) rep(i,N) mn[N+i] = 0;\n    if(dobuild) build();\n  }\n \n  void build(void){\n    int i;\n    for(i=N-1;i;i--){\n      mn[i] = min(mn[2i], mn[2i+1]);\n    }\n  }\n \n  inline void build(int a){\n    while(a > 1){\n      a /= 2;\n      mn[a] = min(mn[2a], mn[2a+1]);\n    }\n  }\n \n  inline void change(int a, T val){\n    mn[a+N] = val;\n    build(a+N);\n  }\n \n  inline void add(int a, T val){\n    mn[a+N] += val;\n    build(a+N);\n  }\n \n  inline T getMinVal(int a, int b){\n    T res, tmp;\n    int fga = 0, fgb = 0;\n\n    a += N;\n    b += N;\n \n    while(a < b){\n      if(a%2){\n        if(fga){\n          res = min(res, mn[a]);\n        } else {\n          res = mn[a];\n          fga = 1;\n        }\n        a++;\n      }\n      if(b%2){\n        b--;\n        if(fgb){\n          tmp = min(mn[b], tmp);\n        } else {\n          tmp = mn[b];\n          fgb = 1;\n        }\n      }\n      a /= 2;\n      b /= 2;\n    }\n    if(fga==1 && fgb==0) return res;\n    if(fga==0 && fgb==1) return tmp;\n    if(fga==1 && fgb==1){\n      res = min(res, tmp);\n      return res;\n    }\n    return res;\n  }\n};\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"min_L");
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "segtree_Point_Maxval";
      string c = "template<class T>\nstruct segtree_Point_Maxval{\n  int N, logN;\n  T *mx;\n\n  void malloc(int maxN, int once = 0){\n    int i;\n    for(i=1;i<maxN;i*=2);\n    mx = new T[2*i];\n    if(once) setN(maxN);\n  }\n\n  void walloc(int maxN, int once = 0, void **mem = &wmem){\n    int i;\n    for(i=1;i<maxN;i*=2);\n    walloc1d(&mx, 2i, mem);\n    if(once) setN(maxN);\n  }\n\n  void free(void){\n    delete [] mx;\n  }\n\n  T& operator[](int i){\n    return mx[N+i];\n  }\n\n  void setN(int n, int zerofill = 1, int dobuild = 1){\n    int i;\n    for(i=1,logN=0;i<n;i*=2,logN++);\n    N = i;\n    if(zerofill) rep(i,N) mx[N+i] = 0;\n    if(dobuild) build();\n  }\n\n  void build(void){\n    int i;\n    for(i=N-1;i;i--){\n      mx[i] = max(mx[2i], mx[2i+1]);\n    }\n  }\n\n  inline void build(int a){\n    while(a > 1){\n      a /= 2;\n      mx[a] = max(mx[2a], mx[2a+1]);\n    }\n  }\n\n  inline void change(int a, T val){\n    mx[a+N] = val;\n    build(a+N);\n  }\n\n  inline void add(int a, T val){\n    mx[a+N] += val;\n    build(a+N);\n  }\n\n  inline T getMaxVal(int a, int b){\n    T res, tmp;\n    int fga = 0, fgb = 0;\n\n    a += N;\n    b += N;\n\n    while(a < b){\n      if(a%2){\n        if(fga){\n          res = max(res, mx[a]);\n        } else {\n          res = mx[a];\n          fga = 1;\n        }\n        a++;\n      }\n      if(b%2){\n        b--;\n        if(fgb){\n          tmp = max(mx[b], tmp);\n        } else {\n          tmp = mx[b];\n          fgb = 1;\n        }\n      }\n      a /= 2;\n      b /= 2;\n    }\n    if(fga==1 && fgb==0) return res;\n    if(fga==0 && fgb==1) return tmp;\n    if(fga==1 && fgb==1){\n      res = max(res, tmp);\n      return res;\n    }\n    return res;\n  }\n};\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"max_L");
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "segtree_Point_Min";
      string c = "template<class T>\nstruct segtree_Point_Min{\n  int N, logN;\n  T *mn;\n  int *mnind;\n\n  void malloc(int maxN, int once = 0){\n    int i;\n    for(i=1;i<maxN;i*=2);\n    mn = new T[2*i];\n    mnind = new int[2*i];\n    if(once) setN(maxN);\n  }\n\n  void walloc(int maxN, int once = 0, void **mem = &wmem){\n    int i;\n    for(i=1;i<maxN;i*=2);\n    walloc1d(&mn, 2i, mem);\n    walloc1d(&mnind, 2i, mem);\n    if(once) setN(maxN);\n  }\n\n  void free(void){\n    delete [] mn;\n    delete [] mnind;\n  }\n\n  T& operator[](int i){\n    return mn[N+i];\n  }\n \n  void setN(int n, int zerofill = 1, int dobuild = 1){\n    int i;\n    for(i=1,logN=0;i<n;i*=2,logN++);\n    N = i;\n    if(zerofill) rep(i,N) mn[N+i] = 0;\n    if(dobuild) build();\n  }\n \n  void build(void){\n    int i;\n    rep(i,N) mnind[N+i] = i;\n    for(i=N-1;i;i--){\n      if(mn[2*i] <= mn[2*i+1]){\n        mn[i] = mn[2*i];\n        mnind[i] = mnind[2*i];\n      } else {\n        mn[i] = mn[2*i+1];\n        mnind[i] = mnind[2*i+1];\n      }\n    }\n  }\n \n  inline void build(int a){\n    while(a > 1){\n      a /= 2;\n      if(mn[a*2] <= mn[a*2+1]){\n        mn[a] = mn[a*2];\n        mnind[a] = mnind[a*2];\n      } else {\n        mn[a] = mn[a*2+1];\n        mnind[a] = mnind[a*2+1];\n      }\n    }\n  }\n \n  inline void change(int a, T val){\n    mn[a+N] = val;\n    build(a+N);\n  }\n \n  inline void add(int a, T val){\n    mn[a+N] += val;\n    build(a+N);\n  }\n \n  inline pair<T,int> getMin(int a, int b){\n    pair<T,int> res, tmp;\n    int fga = 0, fgb = 0;\n    \n    a += N;\n    b += N;\n \n    while(a < b){\n      if(a%2){\n        if(fga){\n          res = min(res, make_pair(mn[a], mnind[a]));\n        } else {\n          res = make_pair(mn[a], mnind[a]);\n          fga = 1;\n        }\n        a++;\n      }\n      if(b%2){\n        b--;\n        if(fgb){\n          tmp = min(make_pair(mn[b], mnind[b]), tmp);\n        } else {\n          tmp = make_pair(mn[b], mnind[b]);\n          fgb = 1;\n        }\n      }\n      a /= 2;\n      b /= 2;\n    }\n    if(fga==1 && fgb==0) return res;\n    if(fga==0 && fgb==1) return tmp;\n    if(fga==1 && fgb==1){\n      res = min(res, tmp);\n      return res;\n    }\n    return res;\n  }\n\n  inline T getMinVal(int a, int b){\n    T res, tmp;\n    int fga = 0, fgb = 0;\n    \n    a += N;\n    b += N;\n \n    while(a < b){\n      if(a%2){\n        if(fga){\n          res = min(res, mn[a]);\n        } else {\n          res = mn[a];\n          fga = 1;\n        }\n        a++;\n      }\n      if(b%2){\n        b--;\n        if(fgb){\n          tmp = min(mn[b], tmp);\n        } else {\n          tmp = mn[b];\n          fgb = 1;\n        }\n      }\n      a /= 2;\n      b /= 2;\n    }\n    if(fga==1 && fgb==0) return res;\n    if(fga==0 && fgb==1) return tmp;\n    if(fga==1 && fgb==1){\n      res = min(res, tmp);\n      return res;\n    }\n    return res;\n  }\n  \n  inline int getMinInd(int a, int b){\n    return getMin(a,b).second;\n  }\n};\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"min_L");
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "segtree_Point_SumMin";
      string c = "template<class T>\nstruct segtree_Point_SumMin{\n  int N, logN;\n  T *sum, *mn;\n  int *mnind;\n\n  void malloc(int maxN, int once = 0){\n    int i;\n    for(i=1;i<maxN;i*=2);\n    \n    sum = new T[2*i];\n    mn = new T[2*i];\n    mnind = new int[2*i];\n\n    if(once) setN(maxN);\n  }\n\n  void walloc(int maxN, int once = 0, void **mem = &wmem){\n    int i;\n    for(i=1;i<maxN;i*=2);\n\n    walloc1d(&sum, 2i, mem);\n    walloc1d(&mn, 2i, mem);\n    walloc1d(&mnind, 2i, mem);\n\n    if(once) setN(maxN);\n  }\n\n  void free(void){\n    delete [] sum;\n    delete [] mn;\n    delete [] mnind;\n  }\n\n  T& operator[](int i){\n    return sum[N+i];\n  }\n \n  void setN(int n, int zerofill = 1, int dobuild = 1){\n    int i;\n    for(i=1,logN=0;i<n;i*=2,logN++);\n    N = i;\n    if(zerofill) rep(i,N) sum[N+i] = 0;\n    if(dobuild) build();\n  }\n \n  void build(void){\n    int i;\n    rep(i,N) mn[N+i] = sum[N+i], mnind[N+i] = i;\n    for(i=N-1;i;i--){\n      sum[i] = sum[2*i] + sum[2*i+1];\n      if(mn[2*i] <= mn[2*i+1]){\n        mn[i] = mn[2*i];\n        mnind[i] = mnind[2*i];\n      } else {\n        mn[i] = mn[2*i+1];\n        mnind[i] = mnind[2*i+1];\n      }\n    }\n  }\n \n  inline void build(int a){\n    while(a > 1){\n      a /= 2;\n      sum[a] = sum[a*2] + sum[a*2+1];\n      if(mn[a*2] <= mn[a*2+1]){\n        mn[a] = mn[a*2];\n        mnind[a] = mnind[a*2];\n      } else {\n        mn[a] = mn[a*2+1];\n        mnind[a] = mnind[a*2+1];\n      }\n    }\n  }\n \n  inline void change(int a, T val){\n    sum[a+N] = mn[a+N] = val;\n    build(a+N);\n  }\n \n  inline void add(int a, T val){\n    sum[a+N] = (mn[a+N] += val);\n    build(a+N);\n  }\n \n  inline pair<T,int> getMin(int a, int b){\n    pair<T,int> res, tmp;\n    int fga = 0, fgb = 0;\n    \n    a += N;\n    b += N;\n \n    while(a < b){\n      if(a%2){\n        if(fga){\n          res = min(res, make_pair(mn[a], mnind[a]));\n        } else {\n          res = make_pair(mn[a], mnind[a]);\n          fga = 1;\n        }\n        a++;\n      }\n      if(b%2){\n        b--;\n        if(fgb){\n          tmp = min(make_pair(mn[b], mnind[b]), tmp);\n        } else {\n          tmp = make_pair(mn[b], mnind[b]);\n          fgb = 1;\n        }\n      }\n      a /= 2;\n      b /= 2;\n    }\n    if(fga==1 && fgb==0) return res;\n    if(fga==0 && fgb==1) return tmp;\n    if(fga==1 && fgb==1){\n      res = min(res, tmp);\n      return res;\n    }\n    return res;\n  }\n\n  inline T getMinVal(int a, int b){\n    T res, tmp;\n    int fga = 0, fgb = 0;\n    \n    a += N;\n    b += N;\n \n    while(a < b){\n      if(a%2){\n        if(fga){\n          res = min(res, mn[a]);\n        } else {\n          res = mn[a];\n          fga = 1;\n        }\n        a++;\n      }\n      if(b%2){\n        b--;\n        if(fgb){\n          tmp = min(mn[b], tmp);\n        } else {\n          tmp = mn[b];\n          fgb = 1;\n        }\n      }\n      a /= 2;\n      b /= 2;\n    }\n    if(fga==1 && fgb==0) return res;\n    if(fga==0 && fgb==1) return tmp;\n    if(fga==1 && fgb==1){\n      res = min(res, tmp);\n      return res;\n    }\n    return res;\n  }\n  \n  inline int getMinInd(int a, int b){\n    return getMin(a,b).second;\n  }\n\n  inline T getSum(int a, int b){\n    T res, tmp;\n    int fga = 0, fgb = 0;\n    \n    a += N;\n    b += N;\n \n    while(a < b){\n      if(a%2){\n        if(fga){\n          res = res + sum[a];\n        } else {\n          res = sum[a];\n          fga = 1;\n        }\n        a++;\n      }\n      if(b%2){\n        b--;\n        if(fgb){\n          tmp = sum[b] + tmp;\n        } else {\n          tmp = sum[b];\n          fgb = 1;\n        }\n      }\n      a /= 2;\n      b /= 2;\n    }\n    if(fga==1 && fgb==0) return res;\n    if(fga==0 && fgb==1) return tmp;\n    if(fga==1 && fgb==1){\n      res = res + tmp;\n      return res;\n    }\n    return res;\n  }\n};\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"min_L");
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "segtree_Point_Prod";
      string c = "template<class T>\nstruct segtree_Point_Prod{\n  int N, logN;\n  T *mul;\n\n  void malloc(int maxN, int once = 0){\n    int i;\n    for(i=1;i<maxN;i*=2);\n    mul = new T[2*i];\n    if(once) setN(maxN);\n  }\n\n  void walloc(int maxN, int once = 0, void **mem = &wmem){\n    int i;\n    for(i=1;i<maxN;i*=2);\n    walloc1d(&mul, 2i, mem);\n    if(once) setN(maxN);\n  }\n\n  void free(void){\n    delete [] mul;\n  }\n\n  T& operator[](int i){\n    return mul[N+i];\n  }\n \n  void setN(int n, int zerofill = 1, int dobuild = 1){\n    int i;\n    for(i=1,logN=0;i<n;i*=2,logN++);\n    N = i;\n    if(zerofill) rep(i,N) mul[N+i] = 0;\n    if(dobuild) build();\n  }\n\n  void build(void){\n    int i;\n    for(i=N-1;i;i--){\n      mul[i] = mul[2*i] * mul[2*i+1];\n    }\n  }\n \n  inline void build(int a){\n    while(a > 1){\n      a /= 2;\n      mul[a] = mul[a*2] * mul[a*2+1];\n    }\n  }\n \n  inline void change(int a, T val){\n    mul[a+N] = val;\n    build(a+N);\n  }\n \n  inline void add(int a, T val){\n    mul[a+N] += val;\n    build(a+N);\n  }\n \n  inline T getProd(int a, int b){\n    T res, tmp;\n    int fga = 0, fgb = 0;\n    \n    a += N;\n    b += N;\n \n    while(a < b){\n      if(a%2){\n        if(fga){\n          res = res * mul[a];\n        } else {\n          res = mul[a];\n          fga = 1;\n        }\n        a++;\n      }\n      if(b%2){\n        b--;\n        if(fgb){\n          tmp = mul[b] * tmp;\n        } else {\n          tmp = mul[b];\n          fgb = 1;\n        }\n      }\n      a /= 2;\n      b /= 2;\n    }\n    if(fga==1 && fgb==0) return res;\n    if(fga==0 && fgb==1) return tmp;\n    if(fga==1 && fgb==1){\n      res = res * tmp;\n      return res;\n    }\n    return res;\n  }\n};\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "segtree_Point_Minval2";
      string c = "template<class T>\nstruct segtree_Point_Minval2{\n  int N, logN;\n  T *mn, *mn2;\n\n  void malloc(int maxN, int once = 0){\n    int i;\n    for(i=1;i<maxN;i*=2);\n    \n    mn = new T[2*i];\n    mn2 = new T[2*i];\n\n    if(once) setN(maxN);\n  }\n\n  void walloc(int maxN, int once = 0, void **mem = &wmem){\n    int i;\n    for(i=1;i<maxN;i*=2);\n\n    walloc1d(&mn, 2i, mem);\n    walloc1d(&mn2, 2i, mem);\n\n    if(once) setN(maxN);\n  }\n\n  void free(void){\n    delete [] mn;\n    delete [] mn2;\n  }\n\n  T& operator[](int i){\n    return mn[N+i];\n  }\n \n  void setN(int n, int zerofill = 1, int dobuild = 1){\n    int i;\n    for(i=1,logN=0;i<n;i*=2,logN++);\n    N = i;\n    if(zerofill) rep(i,N) mn[N+i] = 0;\n    if(dobuild) build();\n  }\n \n  void build(void){\n    int i;\n    rep(i,N) mn2[N+i] = numeric_limits<T>::max();\n    for(i=N-1;i;i--){\n      if(mn[2i] <= mn[2i+1]){\n        mn[i] = mn[2i];\n        mn2[i] = min(mn2[2i], mn[2i+1]);\n      } else {\n        mn[i] = mn[2i+1];\n        mn2[i] = min(mn[2i], mn2[2i+1]);\n      }\n    }\n  }\n \n  inline void build(int a){\n    while(a > 1){\n      a /= 2;\n      if(mn[2a] <= mn[2a+1]){\n        mn[a] = mn[2a];\n        mn2[a] = min(mn2[2a], mn[2a+1]);\n      } else {\n        mn[a] = mn[2a+1];\n        mn2[a] = min(mn[2a], mn2[2a+1]);\n      }\n    }\n  }\n \n  inline void change(int a, T val){\n    mn[a+N] = val;\n    build(a+N);\n  }\n \n  inline void add(int a, T val){\n    mn[a+N] += val;\n    build(a+N);\n  }\n \n  inline pair<T,T> getMinVal2(int a, int b){\n    T res1, res2;\n\n    a += N;\n    b += N;\n\n    res1 = res2 = numeric_limits<T>::max();\n    while(a < b){\n      if(a%2){\n        res2 <?= mn[a];\n        sortE(res1, res2);\n        res2 <?= mn2[a];\n        a++;\n      }\n      if(b%2){\n        b--;\n        res2 <?= mn[b];\n        sortE(res1, res2);\n        res2 <?= mn2[b];\n      }\n      a /= 2;\n      b /= 2;\n    }\n    return make_pair(res1, res2);\n  }\n\n  inline T getMinVal(int a, int b){\n    T res, tmp;\n    int fga = 0, fgb = 0;\n\n    a += N;\n    b += N;\n \n    while(a < b){\n      if(a%2){\n        if(fga){\n          res = min(res, mn[a]);\n        } else {\n          res = mn[a];\n          fga = 1;\n        }\n        a++;\n      }\n      if(b%2){\n        b--;\n        if(fgb){\n          tmp = min(mn[b], tmp);\n        } else {\n          tmp = mn[b];\n          fgb = 1;\n        }\n      }\n      a /= 2;\n      b /= 2;\n    }\n    if(fga==1 && fgb==0) return res;\n    if(fga==0 && fgb==1) return tmp;\n    if(fga==1 && fgb==1){\n      res = min(res, tmp);\n      return res;\n    }\n    return res;\n  }\n};\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"min_L");
      d.push_back((string)"chmin");
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "segtree_Point_Or";
      string c = "template<class T>\nstruct segtree_Point_Or{\n  int N, logN;\n  T *dat;\n\n  void malloc(int maxN, int once = 0){\n    int i;\n    for(i=1;i<maxN;i*=2);\n    dat = new T[2*i];\n    if(once) setN(maxN);\n  }\n\n  void walloc(int maxN, int once = 0, void **mem = &wmem){\n    int i;\n    for(i=1;i<maxN;i*=2);\n    walloc1d(&dat, 2i, mem);\n    if(once) setN(maxN);\n  }\n\n  void free(void){\n    delete [] dat;\n  }\n\n  T& operator[](int i){\n    return dat[N+i];\n  }\n \n  void setN(int n, int zerofill = 1, int dobuild = 1){\n    int i;\n    for(i=1,logN=0;i<n;i*=2,logN++);\n    N = i;\n    if(zerofill) rep(i,N) dat[N+i] = 0;\n    if(dobuild) build();\n  }\n \n  void build(void){\n    int i;\n    for(i=N-1;i;i--){\n      dat[i] = (dat[2i] | dat[2i+1]);\n    }\n  }\n \n  inline void build(int a){\n    while(a > 1){\n      a /= 2;\n      dat[a] = (dat[2a] | dat[2a+1]);\n    }\n  }\n \n  inline void change(int a, T val){\n    dat[a+N] = val;\n    build(a+N);\n  }\n \n  inline void add(int a, T val){\n    dat[a+N] += val;\n    build(a+N);\n  }\n \n  inline T getOr(int a, int b){\n    T res, tmp;\n    int fga = 0, fgb = 0;\n\n    a += N;\n    b += N;\n \n    while(a < b){\n      if(a%2){\n        if(fga){\n          res = (res | dat[a]);\n        } else {\n          res = dat[a];\n          fga = 1;\n        }\n        a++;\n      }\n      if(b%2){\n        b--;\n        if(fgb){\n          tmp = (dat[b] | tmp);\n        } else {\n          tmp = dat[b];\n          fgb = 1;\n        }\n      }\n      a /= 2;\n      b /= 2;\n    }\n    if(fga==1 && fgb==0) return res;\n    if(fga==0 && fgb==1) return tmp;\n    if(fga==1 && fgb==1){\n      res = (res | tmp);\n      return res;\n    }\n    return res;\n  }\n};\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "segtree_Point_And";
      string c = "template<class T>\nstruct segtree_Point_And{\n  int N, logN;\n  T *dat;\n\n  void malloc(int maxN, int once = 0){\n    int i;\n    for(i=1;i<maxN;i*=2);\n    dat = new T[2*i];\n    if(once) setN(maxN);\n  }\n\n  void walloc(int maxN, int once = 0, void **mem = &wmem){\n    int i;\n    for(i=1;i<maxN;i*=2);\n    walloc1d(&dat, 2i, mem);\n    if(once) setN(maxN);\n  }\n\n  void free(void){\n    delete [] dat;\n  }\n\n  T& operator[](int i){\n    return dat[N+i];\n  }\n \n  void setN(int n, int zerofill = 1, int dobuild = 1){\n    int i;\n    for(i=1,logN=0;i<n;i*=2,logN++);\n    N = i;\n    if(zerofill) rep(i,N) dat[N+i] = 0;\n    if(dobuild) build();\n  }\n \n  void build(void){\n    int i;\n    for(i=N-1;i;i--){\n      dat[i] = (dat[2i] & dat[2i+1]);\n    }\n  }\n \n  inline void build(int a){\n    while(a > 1){\n      a /= 2;\n      dat[a] = (dat[2a] & dat[2a+1]);\n    }\n  }\n \n  inline void change(int a, T val){\n    dat[a+N] = val;\n    build(a+N);\n  }\n \n  inline void add(int a, T val){\n    dat[a+N] += val;\n    build(a+N);\n  }\n \n  inline T getAnd(int a, int b){\n    T res, tmp;\n    int fga = 0, fgb = 0;\n\n    a += N;\n    b += N;\n \n    while(a < b){\n      if(a%2){\n        if(fga){\n          res = (res & dat[a]);\n        } else {\n          res = dat[a];\n          fga = 1;\n        }\n        a++;\n      }\n      if(b%2){\n        b--;\n        if(fgb){\n          tmp = (dat[b] & tmp);\n        } else {\n          tmp = dat[b];\n          fgb = 1;\n        }\n      }\n      a /= 2;\n      b /= 2;\n    }\n    if(fga==1 && fgb==0) return res;\n    if(fga==0 && fgb==1) return tmp;\n    if(fga==1 && fgb==1){\n      res = (res & tmp);\n      return res;\n    }\n    return res;\n  }\n};\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "segtree_Point_Xor";
      string c = "template<class T>\nstruct segtree_Point_Xor{\n  int N, logN;\n  T *dat;\n\n  void malloc(int maxN, int once = 0){\n    int i;\n    for(i=1;i<maxN;i*=2);\n    dat = new T[2*i];\n    if(once) setN(maxN);\n  }\n\n  void walloc(int maxN, int once = 0, void **mem = &wmem){\n    int i;\n    for(i=1;i<maxN;i*=2);\n    walloc1d(&dat, 2i, mem);\n    if(once) setN(maxN);\n  }\n\n  void free(void){\n    delete [] dat;\n  }\n\n  T& operator[](int i){\n    return dat[N+i];\n  }\n \n  void setN(int n, int zerofill = 1, int dobuild = 1){\n    int i;\n    for(i=1,logN=0;i<n;i*=2,logN++);\n    N = i;\n    if(zerofill) rep(i,N) dat[N+i] = 0;\n    if(dobuild) build();\n  }\n \n  void build(void){\n    int i;\n    for(i=N-1;i;i--){\n      dat[i] = (dat[2i] ^ dat[2i+1]);\n    }\n  }\n \n  inline void build(int a){\n    while(a > 1){\n      a /= 2;\n      dat[a] = (dat[2a] ^ dat[2a+1]);\n    }\n  }\n \n  inline void change(int a, T val){\n    dat[a+N] = val;\n    build(a+N);\n  }\n \n  inline void add(int a, T val){\n    dat[a+N] += val;\n    build(a+N);\n  }\n \n  inline T getXor(int a, int b){\n    T res, tmp;\n    int fga = 0, fgb = 0;\n\n    a += N;\n    b += N;\n \n    while(a < b){\n      if(a%2){\n        if(fga){\n          res = (res ^ dat[a]);\n        } else {\n          res = dat[a];\n          fga = 1;\n        }\n        a++;\n      }\n      if(b%2){\n        b--;\n        if(fgb){\n          tmp = (dat[b] ^ tmp);\n        } else {\n          tmp = dat[b];\n          fgb = 1;\n        }\n      }\n      a /= 2;\n      b /= 2;\n    }\n    if(fga==1 && fgb==0) return res;\n    if(fga==0 && fgb==1) return tmp;\n    if(fga==1 && fgb==1){\n      res = (res ^ tmp);\n      return res;\n    }\n    return res;\n  }\n};\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "segtree_Add_Minval";
      string c = "template<class T>\nstruct segtree_Add_Minval{\n  int N, logN;\n  T *mn;\n  T *addval;\n\n  void malloc(int maxN, int once = 0){\n    int i;\n    for(i=1;i<maxN;i*=2);\n\n    mn = new T[2*i];\n    addval = new T[i];\n\n    if(once) setN(maxN);\n  }\n\n  void walloc(int maxN, int once = 0, void **mem = &wmem){\n    int i;\n    for(i=1;i<maxN;i*=2);\n\n    walloc1d(&mn, 2i, mem);\n    walloc1d(&addval, i, mem);\n\n    if(once) setN(maxN);\n  }\n\n  void free(void){\n    delete [] mn;\n    delete [] addval;\n  }\n\n  T& operator[](int i){\n    return mn[N+i];\n  }\n \n  void setN(int n, int zerofill = 1, int dobuild = 1){\n    int i;\n    for(i=1,logN=0;i<n;i*=2,logN++);\n    N = i;\n    if(zerofill) rep(i,N) mn[N+i] = 0;\n    if(dobuild) build();\n  }\n\n  void build(void){\n    int i;\n    for(i=N-1;i;i--) mn[i] = min(mn[2*i], mn[2*i+1]);\n    REP(i,1,N) addval[i] = 0;\n  }\n\n  inline void push_one(int a, int sz, int st){\n    if(addval[a] != 0){\n      if(sz > 1){\n        addval[a*2] += addval[a];\n        addval[a*2+1] += addval[a];\n        mn[a*2] += addval[a];\n        mn[a*2+1] += addval[a];\n      } else {\n        mn[a*2] += addval[a];\n        mn[a*2+1] += addval[a];\n      }\n      addval[a] = 0;\n      return;\n    }\n  }\n\n  inline void push(int a){\n    int i, aa = a - N, nd, sz, st;\n    for(i=logN;i;i--){\n      nd = a>>i;\n      sz = 1<<(i-1);\n      st = 2 * sz * (aa>>i);\n      push_one(nd, sz, st);\n    }\n  }\n\n  inline void build(int a){\n    int sz = 1, st = a - N;\n    while(a > 1){\n      if(a%2) st += sz;\n      a /= 2;\n      sz *= 2;\n      mn[a] = min(mn[a*2], mn[a*2+1]);\n      if(addval[a] != 0) mn[a] += addval[a];\n    }\n  }\n\n  inline void add(int a, int b, T val){\n    int sz = 1, aa, bb;\n    if(a >= b) return;\n\n    aa = (a += N);\n    bb = (b += N);\n    push(a); push(b-1);\n\n    if(a%2){\n      mn[a] += val;\n      a++;\n    }\n    if(b%2){\n      b--;\n      mn[b] += val;\n    }\n    a /= 2;\n    b /= 2;\n\n    while(a < b){\n      sz *= 2;\n      if(a%2){\n        addval[a] += val;\n        mn[a] += val;\n        a++;\n      }\n      if(b%2){\n        b--;\n        addval[b] += val;\n        mn[b] += val;\n      }\n      a /= 2;\n      b /= 2;\n    }\n\n    build(aa);\n    build(bb-1);\n  }\n\n  inline T getMinVal(int a, int b){\n    T res, tmp;\n    int fga = 0, fgb = 0;\n    \n    a += N;\n    b += N;\n    push(a); push(b-1);\n\n    while(a < b){\n      if(a%2){\n        if(fga){\n          res = min(res, mn[a]);\n        } else {\n          res = mn[a];\n          fga = 1;\n        }\n        a++;\n      }\n      if(b%2){\n        b--;\n        if(fgb){\n          tmp = min(mn[b], tmp);\n        } else {\n          tmp = mn[b];\n          fgb = 1;\n        }\n      }\n      a /= 2;\n      b /= 2;\n    }\n    if(fga==1 && fgb==0) return res;\n    if(fga==0 && fgb==1) return tmp;\n    if(fga==1 && fgb==1){\n      res = min(res, tmp);\n      return res;\n    }\n    return res;\n  }\n};\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"min_L");
      d.push_back((string)"chmin");
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "segtree_ChangeAdd_Sum";
      string c = "template<class T>\nstruct segtree_ChangeAdd_Sum{\n  int N, logN;\n  T *sum;\n\n  T *fixval; char *fixed;\n  T *addval;\n\n  void malloc(int maxN, int once = 0){\n    int i;\n    for(i=1;i<maxN;i*=2);\n\n    sum = new T[2*i];\n    fixval = new T[i];\n    addval = new T[i];\n    fixed = new char[i];\n\n    if(once) setN(maxN);\n  }\n\n  void walloc(int maxN, int once = 0, void **mem = &wmem){\n    int i;\n    for(i=1;i<maxN;i*=2);\n\n    walloc1d(&sum, 2i, mem);\n    walloc1d(&fixval, i, mem);\n    walloc1d(&addval, i, mem);\n    walloc1d(&fixed, i, mem);\n\n    if(once) setN(maxN);\n  }\n\n  void free(void){\n    delete [] sum;\n    delete [] fixval;\n    delete [] addval;\n    delete [] fixed;\n  }\n\n  T& operator[](int i){\n    return sum[N+i];\n  }\n \n  void setN(int n, int zerofill = 1, int dobuild = 1){\n    int i;\n    for(i=1,logN=0;i<n;i*=2,logN++);\n    N = i;\n    if(zerofill) rep(i,N) sum[N+i] = 0;\n    if(dobuild) build();\n  }\n \n  void build(void){\n    int i;\n    for(i=N-1;i;i--){\n      sum[i] = sum[2*i] + sum[2*i+1];\n    }\n    REP(i,1,N) fixed[i] = 0;\n    REP(i,1,N) addval[i] = 0;\n  }\n \n  inline void push_one(int a, int sz, int st){\n    if(fixed[a]){\n      if(sz > 1){\n        fixed[a*2] = fixed[a*2+1] = 1;\n        fixval[a*2] = fixval[a*2+1] = fixval[a];\n        sum[a*2] = sum[a*2+1] = sz * fixval[a];\n      } else {\n        sum[a*2] = sum[a*2+1] = sz * fixval[a];\n      }\n      fixed[a] = 0;\n      addval[a] = 0;\n      return;\n    }\n    if(addval[a] != 0){\n      if(sz > 1){\n        if(fixed[a*2]) fixval[a*2] += addval[a];\n        else           addval[a*2] += addval[a];\n        if(fixed[a*2+1]) fixval[a*2+1] += addval[a];\n        else             addval[a*2+1] += addval[a];\n        sum[a*2] += sz * addval[a];\n        sum[a*2+1] += sz * addval[a];\n      } else {\n        sum[a*2] += sz * addval[a];\n        sum[a*2+1] += sz * addval[a];\n      }\n      addval[a] = 0;\n      return;\n    }\n  }\n \n  inline void push(int a){\n    int i, aa = a - N, nd, sz, st;\n    for(i=logN;i;i--){\n      nd = a>>i;\n      sz = 1<<(i-1);\n      st = 2 * sz * (aa>>i);\n      push_one(nd, sz, st);\n    }\n  }\n \n  inline void build(int a){\n    int sz = 1, st = a - N;\n    while(a > 1){\n      if(a%2) st += sz;\n      a /= 2;\n      sz *= 2;\n      if(fixed[a]){\n        sum[a] = sz * fixval[a];\n      } else {\n        sum[a] = sum[a*2] + sum[a*2+1];\n        if(addval[a] != 0){\n          sum[a] += sz * addval[a];\n        }\n      }\n    }\n  }\n \n  inline void change(int a, int b, T val){\n    int sz = 1, aa, bb, st_a = a, st_b = b;\n    if(a >= b) return;\n \n    aa = (a += N);\n    bb = (b += N);\n    push(a); push(b-1);\n \n    if(a%2){\n      sum[a] = val;\n      a++;\n      st_a += sz;\n    }\n    if(b%2){\n      b--;\n      st_b -= sz;\n      sum[b] = val;\n    }\n    a /= 2;\n    b /= 2;\n \n    while(a < b){\n      sz *= 2;\n      if(a%2){\n        fixed[a]=1, fixval[a]=val;\n        sum[a] = sz * val;\n        a++;\n        st_a += sz;\n      }\n      if(b%2){\n        b--;\n        st_b -= sz;\n        fixed[b]=1, fixval[b]=val;\n        sum[b] = sz * val;\n      }\n      a /= 2;\n      b /= 2;\n    }\n \n    build(aa);\n    build(bb-1);\n  }\n \n  inline void add(int a, int b, T val){\n    int sz = 1, aa, bb;\n    if(a >= b) return;\n \n    aa = (a += N);\n    bb = (b += N);\n    push(a); push(b-1);\n \n    if(a%2){\n      sum[a] += val;\n      a++;\n    }\n    if(b%2){\n      b--;\n      sum[b] += val;\n    }\n    a /= 2;\n    b /= 2;\n \n    while(a < b){\n      sz *= 2;\n      if(a%2){\n        if(fixed[a]) fixval[a] += val; else addval[a] += val;\n        sum[a] += sz * val;\n        a++;\n      }\n      if(b%2){\n        b--;\n        if(fixed[b]) fixval[b] += val; else addval[b] += val;\n        sum[b] += sz * val;\n      }\n      a /= 2;\n      b /= 2;\n    }\n \n    build(aa);\n    build(bb-1);\n  }\n \n  inline T getSum(int a, int b){\n    T res, tmp;\n    int fga = 0, fgb = 0;\n    \n    a += N;\n    b += N;\n    push(a); push(b-1);\n \n    while(a < b){\n      if(a%2){\n        if(fga){\n          res = res + sum[a];\n        } else {\n          res = sum[a];\n          fga = 1;\n        }\n        a++;\n      }\n      if(b%2){\n        b--;\n        if(fgb){\n          tmp = sum[b] + tmp;\n        } else {\n          tmp = sum[b];\n          fgb = 1;\n        }\n      }\n      a /= 2;\n      b /= 2;\n    }\n    if(fga==1 && fgb==0) return res;\n    if(fga==0 && fgb==1) return tmp;\n    if(fga==1 && fgb==1){\n      res = res + tmp;\n      return res;\n    }\n    return res;\n  }\n};\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "segtree_ChangeP1add_Sum";
      string c = "template<class T>\nstruct segtree_ChangeP1add_Sum{\n  int N, logN;\n  T *sum;\n\n  T *fixval; char *fixed;\n  T *addval, *p1addval;\n\n  void malloc(int maxN, int once = 0){\n    int i;\n    for(i=1;i<maxN;i*=2);\n\n    sum = new T[2*i];\n    fixval = new T[i];\n    addval = new T[i];\n    p1addval = new T[i];\n    fixed = new char[i];\n\n    if(once) setN(maxN);\n  }\n\n  void walloc(int maxN, int once = 0, void **mem = &wmem){\n    int i;\n    for(i=1;i<maxN;i*=2);\n\n    walloc1d(&sum, 2i, mem);\n    walloc1d(&fixval, i, mem);\n    walloc1d(&addval, i, mem);\n    walloc1d(&p1addval, i, mem);\n    walloc1d(&fixed, i, mem);\n\n    if(once) setN(maxN);\n  }\n\n  void free(void){\n    delete [] sum;\n    delete [] fixval;\n    delete [] addval;\n    delete [] p1addval;\n    delete [] fixed;\n  }\n\n  T& operator[](int i){\n    return sum[N+i];\n  }\n \n  void setN(int n, int zerofill = 1, int dobuild = 1){\n    int i;\n    for(i=1,logN=0;i<n;i*=2,logN++);\n    N = i;\n    if(zerofill) rep(i,N) sum[N+i] = 0;\n    if(dobuild) build();\n  }\n \n  void build(void){\n    int i;\n    for(i=N-1;i;i--){\n      sum[i] = sum[2*i] + sum[2*i+1];\n    }\n    REP(i,1,N) fixed[i] = 0;\n    REP(i,1,N) addval[i] = p1addval[i] = 0;\n  }\n \n  inline void push_one(int a, int sz, int st){\n    if(fixed[a]){\n      if(sz > 1){\n        fixed[a*2] = fixed[a*2+1] = 1;\n        fixval[a*2] = fixval[a*2+1] = fixval[a];\n        addval[a*2] = p1addval[a*2] = 0;\n        addval[a*2+1] = p1addval[a*2+1] = 0;\n        sum[a*2] = sum[a*2+1] = sz * fixval[a];\n      } else {\n        sum[a*2] = sum[a*2+1] = sz * fixval[a];\n      }\n      fixed[a] = 0;\n    }\n    if(addval[a] != 0 || p1addval[a] != 0){\n      if(sz > 1){\n        addval[a*2] += addval[a];\n        addval[a*2+1] += addval[a] + sz * p1addval[a];\n        p1addval[a*2] += p1addval[a];\n        p1addval[a*2+1] += p1addval[a];\n      }\n      sum[a*2] += sz * addval[a] + (ll) sz * (sz-1) / 2 * p1addval[a];\n      sum[a*2+1] += sz * (addval[a] + sz * p1addval[a]) + (ll) sz * (sz-1) / 2 * p1addval[a];\n      addval[a] = p1addval[a] = 0;\n    }\n  }\n\n  inline void push(int a){\n    int i, aa = a - N, nd, sz, st;\n    for(i=logN;i;i--){\n      nd = a>>i;\n      sz = 1<<(i-1);\n      st = 2 * sz * (aa>>i);\n      push_one(nd, sz, st);\n    }\n  }\n \n  inline void build(int a){\n    int sz = 1, st = a - N;\n    while(a > 1){\n      if(a%2) st += sz;\n      a /= 2;\n      sz *= 2;\n      if(fixed[a]){\n        sum[a] = sz * fixval[a];\n      } else {\n        sum[a] = sum[a*2] + sum[a*2+1];\n      }\n      if(addval[a] != 0 || p1addval[a] != 0){\n        sum[a] += sz * addval[a] + (ll) sz * (sz-1) / 2 * p1addval[a];\n      }\n    }\n  }\n\n  inline void change(int a, int b, T val){\n    int sz = 1, aa, bb, st_a = a, st_b = b;\n    if(a >= b) return;\n\n    aa = (a += N);\n    bb = (b += N);\n    push(a); push(b-1);\n\n    if(a%2){\n      sum[a] = val;\n      a++;\n      st_a += sz;\n    }\n    if(b%2){\n      b--;\n      st_b -= sz;\n      sum[b] = val;\n    }\n    a /= 2;\n    b /= 2;\n\n    while(a < b){\n      sz *= 2;\n      if(a%2){\n        fixed[a]=1, fixval[a]=val;\n        addval[a] = p1addval[a] = 0;\n        sum[a] = sz * val;\n        a++;\n        st_a += sz;\n      }\n      if(b%2){\n        b--;\n        st_b -= sz;\n        fixed[b]=1, fixval[b]=val;\n        addval[b] = p1addval[b] = 0;\n        sum[b] = sz * val;\n      }\n      a /= 2;\n      b /= 2;\n    }\n\n    build(aa);\n    build(bb-1);\n  }\n\n  inline void add(int a, int b, T val){\n    int sz = 1, aa, bb;\n    if(a >= b) return;\n\n    aa = (a += N);\n    bb = (b += N);\n    push(a); push(b-1);\n\n    if(a%2){\n      sum[a] += val;\n      a++;\n    }\n    if(b%2){\n      b--;\n      sum[b] += val;\n    }\n    a /= 2;\n    b /= 2;\n\n    while(a < b){\n      sz *= 2;\n      if(a%2){\n        addval[a] += val;\n        sum[a] += sz * val;\n        a++;\n      }\n      if(b%2){\n        b--;\n        addval[b] += val;\n        sum[b] += sz * val;\n      }\n      a /= 2;\n      b /= 2;\n    }\n\n    build(aa);\n    build(bb-1);\n  }\n\n  inline void p1add(int a, int b, T x1, T x0){\n    int sz = 1, aa, bb;\n    T xa, xb;\n    if(a >= b) return;\n\n    aa = (a += N);\n    bb = (b += N);\n    push(a); push(b-1);\n\n    xa = x0;\n    xb = x1 * (b-a) + x0;\n\n    if(a%2){\n      sum[a] += xa;\n      a++;\n      xa += sz * x1;\n    }\n    if(b%2){\n      b--;\n      xb -= sz * x1;\n      sum[b] += xb;\n    }\n    a /= 2;\n    b /= 2;\n\n    while(a < b){\n      sz *= 2;\n      if(a%2){\n        addval[a] += xa;\n        p1addval[a] += x1;\n        sum[a] += sz * xa + (ll) sz * (sz-1) / 2 * x1;\n        a++;\n        xa += sz * x1;\n      }\n      if(b%2){\n        b--;\n        xb -= sz * x1;\n        addval[b] += xb;\n        p1addval[b] += x1;\n        sum[b] += sz * xb + (ll) sz * (sz-1) / 2 * x1;\n      }\n      a /= 2;\n      b /= 2;\n    }\n\n    build(aa);\n    build(bb-1);\n  }\n\n  inline T getSum(int a, int b){\n    T res, tmp;\n    int fga = 0, fgb = 0;\n    \n    a += N;\n    b += N;\n    push(a); push(b-1);\n\n    while(a < b){\n      if(a%2){\n        if(fga){\n          res = res + sum[a];\n        } else {\n          res = sum[a];\n          fga = 1;\n        }\n        a++;\n      }\n      if(b%2){\n        b--;\n        if(fgb){\n          tmp = sum[b] + tmp;\n        } else {\n          tmp = sum[b];\n          fgb = 1;\n        }\n      }\n      a /= 2;\n      b /= 2;\n    }\n    if(fga==1 && fgb==0) return res;\n    if(fga==0 && fgb==1) return tmp;\n    if(fga==1 && fgb==1){\n      res = res + tmp;\n      return res;\n    }\n    return res;\n  }\n};\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "segtree_pg_header";
      string c = "template<class T>\nstruct segtree_pg{\n  int N, logN;\n  T *val;\n\n  void malloc(int maxN, int once = 0);\n  void walloc(int maxN, int once = 0, void **mem = &wmem);\n  void free(void);\n  T& operator[](int i);\n  void setN(int n, int zerofill = 1, int dobuild = 1);\n  void build(void);\n  inline void build(int a);\n  inline void change(int a, T v);\n  inline T get(int a, int b);\n};\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "segtree_pg";
      string c = "template<class T> void segtree_pg<T>::malloc(int maxN, int once /*= 0*/){\n  int i;\n  for(i=1;i<maxN;i*=2);\n  val = new T[2*i];\n  if(once) setN(maxN);\n}\n\ntemplate<class T> void segtree_pg<T>::walloc(int maxN, int once /*= 0*/, void **mem /*= &wmem*/){\n  int i;\n  for(i=1;i<maxN;i*=2);\n  walloc1d(&val, 2i, mem);\n  if(once) setN(maxN);\n}\n\ntemplate<class T> void segtree_pg<T>::free(void){\n  delete [] val;\n}\n\ntemplate<class T> T& segtree_pg<T>::operator[](int i){\n  return val[N+i];\n}\n\ntemplate<class T> void segtree_pg<T>::setN(int n, int zerofill /*= 1*/, int dobuild /*= 1*/){\n  int i;\n  for(i=1,logN=0;i<n;i*=2,logN++);\n  N = i;\n  // if(zerofill) rep(i,N) val[N+i] = 0;\n  if(dobuild) build();\n}\n\ntemplate<class T> void segtree_pg<T>::build(void){\n  for(int i=N-1;i;i--) segtree_pg_func(val[i], val[2i], val[2i+1]);\n}\n\ntemplate<class T> inline void segtree_pg<T>::build(int a){\n  while(a > 1){\n    a /= 2;\n    segtree_pg_func(val[a], val[2a], val[2a+1]);\n  }\n}\n\ntemplate<class T> inline void segtree_pg<T>::change(int a, T v){\n  val[a+N] = v;\n  build(a+N);\n}\n\ntemplate<class T> inline T segtree_pg<T>::get(int a, int b){\n  T res, tmp;\n  int fga = 0, fgb = 0;\n\n  a += N;\n  b += N;\n\n  while(a < b){\n    if(a%2){\n      if(fga){\n        segtree_pg_func(res, res, val[a]);\n      } else {\n        res = val[a];\n        fga = 1;\n      }\n      a++;\n    }\n    if(b%2){\n      b--;\n      if(fgb){\n        segtree_pg_func(tmp, val[b], tmp);\n      } else {\n        tmp = val[b];\n        fgb = 1;\n      }\n    }\n    a /= 2;\n    b /= 2;\n  }\n  if(fga==1 && fgb==0) return res;\n  if(fga==0 && fgb==1) return tmp;\n  if(fga==1 && fgb==1){\n    segtree_pg_func(res, res, tmp);\n    return res;\n  }\n  return res;\n}\n";
      string p = "last";
      vector<string> d;
      d.push_back((string)"segtree_pg_header");
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "segtree_ph_header";
      string c = "template<class T>\nstruct segtree_ph{\n  int N, logN;\n  T *val;\n\n  void malloc(int maxN, int once = 0);\n  void walloc(int maxN, int once = 0, void **mem = &wmem);\n  void free(void);\n  T& operator[](int i);\n  void setN(int n, int zerofill = 1, int dobuild = 1);\n  void build(void);\n  inline void build(int a);\n  inline void change(int a, T v);\n  inline T get(int a, int b);\n};\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "segtree_ph";
      string c = "template<class T> void segtree_ph<T>::malloc(int maxN, int once /*= 0*/){\n  int i;\n  for(i=1;i<maxN;i*=2);\n  val = new T[2*i];\n  if(once) setN(maxN);\n}\n\ntemplate<class T> void segtree_ph<T>::walloc(int maxN, int once /*= 0*/, void **mem /*= &wmem*/){\n  int i;\n  for(i=1;i<maxN;i*=2);\n  walloc1d(&val, 2i, mem);\n  if(once) setN(maxN);\n}\n\ntemplate<class T> void segtree_ph<T>::free(void){\n  delete [] val;\n}\n\ntemplate<class T> T& segtree_ph<T>::operator[](int i){\n  return val[N+i];\n}\n\ntemplate<class T> void segtree_ph<T>::setN(int n, int zerofill /*= 1*/, int dobuild /*= 1*/){\n  int i;\n  for(i=1,logN=0;i<n;i*=2,logN++);\n  N = i;\n  if(dobuild) build();\n}\n\ntemplate<class T> void segtree_ph<T>::build(void){\n  for(int i=N-1;i;i--) val[i] = segtree_ph_func(val[2i], val[2i+1]);\n}\n\ntemplate<class T> inline void segtree_ph<T>::build(int a){\n  while(a > 1){\n    a /= 2;\n    val[a] = segtree_ph_func(val[2a], val[2a+1]);\n  }\n}\n\ntemplate<class T> inline void segtree_ph<T>::change(int a, T v){\n  val[a+N] = v;\n  build(a+N);\n}\n\ntemplate<class T> inline T segtree_ph<T>::get(int a, int b){\n  T res, tmp;\n  int fga = 0, fgb = 0;\n\n  a += N;\n  b += N;\n\n  while(a < b){\n    if(a%2){\n      if(fga){\n        res = segtree_ph_func(res, val[a]);\n      } else {\n        res = val[a];\n        fga = 1;\n      }\n      a++;\n    }\n    if(b%2){\n      b--;\n      if(fgb){\n        tmp = segtree_ph_func(val[b], tmp);\n      } else {\n        tmp = val[b];\n        fgb = 1;\n      }\n    }\n    a /= 2;\n    b /= 2;\n  }\n  if(fga==1 && fgb==0) return res;\n  if(fga==0 && fgb==1) return tmp;\n  if(fga==1 && fgb==1) return segtree_ph_func(res, tmp);\n  return res;\n}\n";
      string p = "last";
      vector<string> d;
      d.push_back((string)"segtree_ph_header");
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "segtree_rg_header";
      string c = "template<class SVAL, class SFUN>\nstruct segtree_rg {\n  int N, logN;\n  SVAL *val;\n  SFUN *fun;\n\n  void malloc(int maxN, int once = 0);\n  void walloc(int maxN, int once = 0, void **mem = &wmem);\n  void free(void);\n  SVAL& operator[](int i);\n  void setN(int n, int zerofill = 1, int dobuild = 1);\n  void build(void);\n  inline void push_one(int a);\n  inline void push(int a);\n  inline void build(int a);\n  inline void change(int a, int b, SFUN f);\n  inline SVAL get(int a, int b);\n};\n\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "segtree_rg";
      string c = "template<class SVAL, class SFUN> void segtree_rg<SVAL, SFUN>::malloc(int maxN, int once /*= 0*/){\n  int i;\n  for(i=1;i<maxN;i*=2);\n  val = new SVAL[2*i];\n  fun = new SFUN[i];\n  if(once) setN(maxN);\n}\n\ntemplate<class SVAL, class SFUN> void segtree_rg<SVAL, SFUN>::walloc(int maxN, int once /*= 0*/, void **mem /*= &wmem*/){\n  int i;\n  for(i=1;i<maxN;i*=2);\n  walloc1d(&val, 2i, mem);\n  walloc1d(&fun, i, mem);\n  if(once) setN(maxN);\n}\n\ntemplate<class SVAL, class SFUN> void segtree_rg<SVAL, SFUN>::free(void){\n  delete [] val;\n  delete [] fun;\n}\n\ntemplate<class SVAL, class SFUN> SVAL& segtree_rg<SVAL, SFUN>::operator[](int i){\n  return val[N+i];\n}\n\ntemplate<class SVAL, class SFUN> void segtree_rg<SVAL, SFUN>::setN(int n, int zerofill /*= 1*/, int dobuild /*= 1*/){\n  int i;\n  for(i=1,logN=0;i<n;i*=2,logN++);\n  N = i;\n  if(dobuild) build();\n}\n\ntemplate<class SVAL, class SFUN> void segtree_rg<SVAL, SFUN>::build(void){\n  int i;\n  for(i=N-1;i;i--){\n    segtree_rg_func(val[i], val[2*i], val[2*i+1]);\n  }\n  REP(i,1,N) segtree_rg_id(fun[i]);\n}\n\ntemplate<class SVAL, class SFUN> inline void segtree_rg<SVAL, SFUN>::push_one(int a){\n  if(2*a < N) segtree_rg_func(fun[2*a], fun[a], fun[2*a]);\n  segtree_rg_func(val[2*a], fun[a], val[2*a]);\n\n  if(2*a+1 < N) segtree_rg_func(fun[2*a+1], fun[a], fun[2*a+1]);\n  segtree_rg_func(val[2*a+1], fun[a], val[2*a+1]);\n\n  segtree_rg_id(fun[a]);\n}\n\ntemplate<class SVAL, class SFUN> inline void segtree_rg<SVAL, SFUN>::push(int a){\n  int i;\n  for(i=logN;i;i--) push_one(a>>i);\n}\n\ntemplate<class SVAL, class SFUN> inline void segtree_rg<SVAL, SFUN>::build(int a){\n  while(a > 1){\n    a /= 2;\n    segtree_rg_func(val[a], val[2*a], val[2*a+1]);\n    segtree_rg_func(val[a], fun[a], val[a]);\n  }\n}\n\ntemplate<class SVAL, class SFUN> inline void segtree_rg<SVAL, SFUN>::change(int a, int b, SFUN f){\n  int aa, bb;\n  if(a >= b) return;\n\n  aa = (a += N);\n  bb = (b += N);\n  push(a); push(b-1);\n\n  if(a%2){\n    segtree_rg_func(val[a], f, val[a]);\n    a++;\n  }\n  if(b%2){\n    b--;\n    segtree_rg_func(val[b], f, val[b]);\n  }\n  a /= 2;\n  b /= 2;\n\n  while(a < b){\n    if(a%2){\n      segtree_rg_func(val[a], f, val[a]);\n      segtree_rg_func(fun[a], f, fun[a]);\n      a++;\n    }\n    if(b%2){\n      b--;\n      segtree_rg_func(val[b], f, val[b]);\n      segtree_rg_func(fun[b], f, fun[b]);\n    }\n    a /= 2;\n    b /= 2;\n  }\n\n  build(aa);\n  build(bb-1);\n}\n\ntemplate<class SVAL, class SFUN> inline SVAL segtree_rg<SVAL, SFUN>::get(int a, int b){\n  SVAL res, tmp;\n  int fga = 0, fgb = 0;\n\n  a += N;\n  b += N;\n  push(a); push(b-1);\n\n  while(a < b){\n    if(a%2){\n      if(fga){\n        segtree_rg_func(res, res, val[a]);\n      } else {\n        res = val[a];\n        fga = 1;\n      }\n      a++;\n    }\n    if(b%2){\n      b--;\n      if(fgb){\n        segtree_rg_func(tmp, val[b], tmp);\n      } else {\n        tmp = val[b];\n        fgb = 1;\n      }\n    }\n    a /= 2;\n    b /= 2;\n  }\n  if(fga==1 && fgb==0) return res;\n  if(fga==0 && fgb==1) return tmp;\n  if(fga==1 && fgb==1){\n    segtree_rg_func(res, res, tmp);\n    return res;\n  }\n  return res;\n}\n";
      string p = "last";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"segtree_rg_header");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "segtree_rh_header";
      string c = "template<class SVAL, class SFUN>\nstruct segtree_rh {\n  int N, logN, trueN;\n  SVAL *val;\n  SFUN *fun;\n  char *dofun;\n\n  void malloc(int maxN, int once = 0);\n  void walloc(int maxN, int once = 0, void **mem = &wmem);\n  void free(void);\n  SVAL& operator[](int i);\n  void setN(int n, int zerofill = 1, int dobuild = 1);\n  void build(void);\n  inline void push_one(int a);\n  inline void push(int a);\n  inline void build(int a);\n  inline void change(int a, int b, SFUN f);\n  inline SVAL get(int a, int b);\n  template<bool (*f)(SVAL)> int max_right(int a, int mx);\n  template<bool (*f)(SVAL)> int max_right(int a);\n  template<bool (*f)(SVAL)> int min_left(int b, int mn);\n  template<bool (*f)(SVAL)> int min_left(int b);\n};\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "segtree_rh";
      string c = "template<class SVAL, class SFUN> void segtree_rh<SVAL, SFUN>::malloc(int maxN, int once /*= 0*/){\n  int i;\n  for(i=1;i<maxN;i*=2);\n  val = new SVAL[2*i];\n  fun = new SFUN[i];\n  dofun = new char[i];\n  if(once) setN(maxN);\n}\n\ntemplate<class SVAL, class SFUN> void segtree_rh<SVAL, SFUN>::walloc(int maxN, int once /*= 0*/, void **mem /*= &wmem*/){\n  int i;\n  for(i=1;i<maxN;i*=2);\n  walloc1d(&val, 2i, mem);\n  walloc1d(&fun, i, mem);\n  walloc1d(&dofun, i, mem);\n  if(once) setN(maxN);\n}\n\ntemplate<class SVAL, class SFUN> void segtree_rh<SVAL, SFUN>::free(void){\n  delete [] val;\n  delete [] fun;\n}\n\ntemplate<class SVAL, class SFUN> SVAL& segtree_rh<SVAL, SFUN>::operator[](int i){\n  return val[N+i];\n}\n\ntemplate<class SVAL, class SFUN> void segtree_rh<SVAL, SFUN>::setN(int n, int zerofill /*= 1*/, int dobuild /*= 1*/){\n  int i;\n  for(i=1,logN=0;i<n;i*=2,logN++);\n  trueN = n;\n  N = i;\n  if(dobuild) build();\n}\n\ntemplate<class SVAL, class SFUN> void segtree_rh<SVAL, SFUN>::build(void){\n  int i;\n  for(i=N-1;i;i--) val[i] = segtree_rh_merge(val[2*i], val[2*i+1]);\n  REP(i,1,N) dofun[i] = 0;\n}\n\ntemplate<class SVAL, class SFUN> inline void segtree_rh<SVAL, SFUN>::push_one(int a){\n  if(dofun[a]){\n    if(2*a < N){\n      if(dofun[2*a]) fun[2*a] = segtree_rh_compose(fun[a], fun[2*a]);\n      else           fun[2*a] = fun[a], dofun[2*a] = 1;\n    }\n    val[2*a] = segtree_rh_apply(fun[a], val[2*a]);\n\n    if(2*a+1 < N){\n      if(dofun[2*a+1]) fun[2*a+1] = segtree_rh_compose(fun[a], fun[2*a+1]);\n      else             fun[2*a+1] = fun[a], dofun[2*a+1] = 1;\n    }\n    val[2*a+1] = segtree_rh_apply(fun[a], val[2*a+1]);\n\n    dofun[a] = 0;\n  }\n}\n\ntemplate<class SVAL, class SFUN> inline void segtree_rh<SVAL, SFUN>::push(int a){\n  int i;\n  for(i=logN;i;i--) push_one(a>>i);\n}\n\ntemplate<class SVAL, class SFUN> inline void segtree_rh<SVAL, SFUN>::build(int a){\n  while(a > 1){\n    a /= 2;\n    val[a] = segtree_rh_merge(val[2*a], val[2*a+1]);\n    if(dofun[a]) val[a] = segtree_rh_apply(fun[a], val[a]);\n  }\n}\n\ntemplate<class SVAL, class SFUN> inline void segtree_rh<SVAL, SFUN>::change(int a, int b, SFUN f){\n  int aa, bb;\n  if(a >= b) return;\n\n  aa = (a += N);\n  bb = (b += N);\n  push(a); push(b-1);\n\n  if(a%2){\n    val[a] = segtree_rh_apply(f, val[a]);\n    a++;\n  }\n  if(b%2){\n    b--;\n    val[b] = segtree_rh_apply(f, val[b]);\n  }\n  a /= 2;\n  b /= 2;\n\n  while(a < b){\n    if(a%2){\n      val[a] = segtree_rh_apply(f, val[a]);\n      if(dofun[a]) fun[a] = segtree_rh_compose(f, fun[a]);\n      else         fun[a] = f, dofun[a] = 1;\n      a++;\n    }\n    if(b%2){\n      b--;\n      val[b] = segtree_rh_apply(f, val[b]);\n      if(dofun[b]) fun[b] = segtree_rh_compose(f, fun[b]);\n      else         fun[b] = f, dofun[b] = 1;\n    }\n    a /= 2;\n    b /= 2;\n  }\n\n  build(aa);\n  build(bb-1);\n}\n\ntemplate<class SVAL, class SFUN> inline SVAL segtree_rh<SVAL, SFUN>::get(int a, int b){\n  SVAL res, tmp;\n  int fga = 0, fgb = 0;\n\n  a += N;\n  b += N;\n  push(a); push(b-1);\n\n  while(a < b){\n    if(a%2){\n      if(fga){\n        res = segtree_rh_merge(res, val[a]);\n      } else {\n        res = val[a];\n        fga = 1;\n      }\n      a++;\n    }\n    if(b%2){\n      b--;\n      if(fgb){\n        tmp = segtree_rh_merge(val[b], tmp);\n      } else {\n        tmp = val[b];\n        fgb = 1;\n      }\n    }\n    a /= 2;\n    b /= 2;\n  }\n  if(fga==1 && fgb==0) return res;\n  if(fga==0 && fgb==1) return tmp;\n  if(fga==1 && fgb==1) return segtree_rh_merge(res, tmp);\n  return res;\n}\n\ntemplate<class SVAL, class SFUN> template<bool (*f)(SVAL)> int segtree_rh<SVAL, SFUN>::max_right(int a, int mx){\n  int fg = 0;\n  int ta = a, sz = 1;\n  SVAL cur, tmp;\n\n  if(a>=mx) return mx;\n  a += N;\n  push(a);\n\n  for(;;){\n    while(a%2==0) a /= 2, sz *= 2;\n\n    if(ta + sz <= mx){\n      if(fg==0){\n        tmp = val[a];\n      } else {\n        tmp = segtree_rh_merge(cur, val[a]);\n      }\n    }\n    if(ta + sz > mx || !f(tmp)){\n      while(a < N){\n        push_one(a);\n        a *= 2;\n        sz /= 2;\n        if(ta + sz <= mx){\n          if(fg==0){\n            tmp = val[a];\n          } else {\n            tmp = segtree_rh_merge(cur, val[a]);\n          }\n        }\n        if(ta + sz <= mx && f(tmp)){\n          fg = 1;\n          cur = tmp;\n          a++;\n          ta += sz;\n        }\n      }\n      return a - N;\n    }\n    fg = 1;\n    cur = tmp;\n    if((a & (a+1)) == 0) break;\n    a++;\n    ta += sz;\n  }\n  return mx;\n}\n\ntemplate<class SVAL, class SFUN> template<bool (*f)(SVAL)> int segtree_rh<SVAL, SFUN>::max_right(int a){\n  return max_right<f>(a, trueN);\n}\n\n\ntemplate<class SVAL, class SFUN> template<bool (*f)(SVAL)> int segtree_rh<SVAL, SFUN>::min_left(int b, int mn){\n  int fg = 0;\n  int tb = b, sz = 1;\n  SVAL cur, tmp;\n\n  if(b <= mn) return mn;\n  b += N;\n  push(b-1);\n\n  for(;;){\n    while(b%2==0) b /= 2, sz *= 2;\n\n    if(tb - sz >= mn){\n      if(fg==0){\n        tmp = val[b-1];\n      } else {\n        tmp = segtree_rh_merge(val[b-1], cur);\n      }\n    }\n    if(tb - sz < mn || !f(tmp)){\n      while(b-1 < N){\n        push_one(b-1);\n        b *= 2;\n        sz /= 2;\n        if(tb - sz >= mn){\n          if(fg==0){\n            tmp = val[b-1];\n          } else {\n            tmp = segtree_rh_merge(val[b-1], cur);\n          }\n        }\n        if(tb - sz >= mn && f(tmp)){\n          fg = 1;\n          cur = tmp;\n          b--;\n          tb -= sz;\n        }\n      }\n      return b - N;\n    }\n    fg = 1;\n    cur = tmp;\n    b--;\n    tb -= sz;\n    if(tb <= mn) break;\n  }\n  return mn;\n}\n\ntemplate<class SVAL, class SFUN> template<bool (*f)(SVAL)> int segtree_rh<SVAL, SFUN>::min_left(int b){\n  return min_left<f>(b, 0);\n}\n";
      string p = "last";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"segtree_rh_header");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "static_segtree_Add_At";
      string c = "template<class T>\nstruct static_segtree_Add_At{\n  int N, mem, has_new;\n  T *dat, *ad;\n\n  void malloc(int maxN, int once = 0){\n    dat = new T[maxN];\n    ad = new T[maxN];\n    if(once) setN(maxN);\n  }\n\n  void walloc(int maxN, int once = 0, void **mem = &wmem){\n    walloc1d(&dat, maxN, mem);\n    walloc1d(&ad, maxN, mem);\n    if(once) setN(maxN);\n  }\n\n  void free(void){\n    delete [] dat;\n    delete [] ad;\n  }\n\n  T& operator[](int i){\n    return dat[i];\n  }\n\n  void setN(int n, int zerofill = 1, int dobuild = 1){\n    N = n;\n    if(zerofill) rep(i,N) dat[i] = 0;\n    if(dobuild) build();\n  }\n\n  void build(void){\n    int i;\n    rep(i,N) ad[i] = 0;\n    has_new = 0;\n  }\n\n  inline void apply(void){\n    if(has_new){\n      int i;\n      T s;\n      s = 0;\n      rep(i,N){\n        dat[i] += (s += ad[i]);\n        ad[i] = 0;\n      }\n    }\n    has_new = 0;\n  }\n\n  inline void add(int a, int b, T val){\n    has_new = 1;\n    ad[a] += val;\n    if(b < N) ad[b] -= val;\n  }\n\n  inline T getAt(int i){\n    apply();\n    return dat[i];\n  }\n};\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "wAdjEdge1";
      string c = "template<class S>\nvoid wAdjEdge_L(const int N, const int M, const int *A, const S *B, int **res_sz, S ***res_B, void **mem = &wmem){\n  int i, j, k;\n  walloc1d(res_sz, N, mem);\n  rep(i,N) (*res_sz)[i] = 0;\n  rep(i,M) (*res_sz)[A[i]]++;\n  walloc1d(res_B, N, mem);\n  rep(i,N) walloc1d(&((*res_B)[i]), (*res_sz)[i], mem);\n  rep(i,N) (*res_sz)[i] = 0;\n  rep(i,M) (*res_B)[A[i]][(*res_sz)[A[i]]++] = B[i];\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "wAdjEdge2";
      string c = "template<class S, class T>\nvoid wAdjEdge_L(const int N, const int M, const int *A, const S *B, const T *C, int **res_sz, S ***res_B, T ***res_C, void **mem = &wmem){\n  int i, j, k;\n  walloc1d(res_sz, N, mem);\n  rep(i,N) (*res_sz)[i] = 0;\n  rep(i,M) (*res_sz)[A[i]]++;\n  walloc1d(res_B, N, mem);\n  rep(i,N) walloc1d(&((*res_B)[i]), (*res_sz)[i], mem);\n  walloc1d(res_C, N, mem);\n  rep(i,N) walloc1d(&((*res_C)[i]), (*res_sz)[i], mem);\n  rep(i,N) (*res_sz)[i] = 0;\n  rep(i,M){\n    (*res_B)[A[i]][(*res_sz)[A[i]]] = B[i];\n    (*res_C)[A[i]][(*res_sz)[A[i]]] = C[i];\n    (*res_sz)[A[i]]++;\n  }\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "wAdjEdge3";
      string c = "template<class S, class T, class U>\nvoid wAdjEdge_L(const int N, const int M, const int *A, const S *B, const T *C, const U *D, int **res_sz, S ***res_B, T ***res_C, U ***res_D, void **mem = &wmem){\n  int i, j, k;\n  walloc1d(res_sz, N, mem);\n  rep(i,N) (*res_sz)[i] = 0;\n  rep(i,M) (*res_sz)[A[i]]++;\n  walloc1d(res_B, N, mem);\n  rep(i,N) walloc1d(&((*res_B)[i]), (*res_sz)[i], mem);\n  walloc1d(res_C, N, mem);\n  rep(i,N) walloc1d(&((*res_C)[i]), (*res_sz)[i], mem);\n  walloc1d(res_D, N, mem);\n  rep(i,N) walloc1d(&((*res_D)[i]), (*res_sz)[i], mem);\n  rep(i,N) (*res_sz)[i] = 0;\n  rep(i,M){\n    (*res_B)[A[i]][(*res_sz)[A[i]]] = B[i];\n    (*res_C)[A[i]][(*res_sz)[A[i]]] = C[i];\n    (*res_D)[A[i]][(*res_sz)[A[i]]] = D[i];\n    (*res_sz)[A[i]]++;\n  }\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "wAdjEdge4";
      string c = "template<class S, class T, class U, class V>\nvoid wAdjEdge_L(const int N, const int M, int *A, const S *B, const T *C, const U *D, const V *E, int **res_sz, S ***res_B, T ***res_C, U ***res_D, V ***res_E, void **mem = &wmem){\n  int i, j, k;\n  walloc1d(res_sz, N, mem);\n  rep(i,N) (*res_sz)[i] = 0;\n  rep(i,M) (*res_sz)[A[i]]++;\n  walloc1d(res_B, N, mem);\n  rep(i,M) walloc1d(&((*res_B)[i]), (*res_sz)[i], mem);\n  walloc1d(res_C, N, mem);\n  rep(i,M) walloc1d(&((*res_C)[i]), (*res_sz)[i], mem);\n  walloc1d(res_D, N, mem);\n  rep(i,M) walloc1d(&((*res_D)[i]), (*res_sz)[i], mem);\n  walloc1d(res_E, N, mem);\n  rep(i,M) walloc1d(&((*res_E)[i]), (*res_sz)[i], mem);\n  rep(i,N) (*res_sz)[i] = 0;\n  rep(i,M){\n    (*res_B)[A[i]][(*res_sz)[A[i]]] = B[i];\n    (*res_C)[A[i]][(*res_sz)[A[i]]] = C[i];\n    (*res_D)[A[i]][(*res_sz)[A[i]]] = D[i];\n    (*res_E)[A[i]][(*res_sz)[A[i]]] = E[i];\n    (*res_sz)[A[i]]++;\n  }\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "graph";
      string c = "struct graph{\n  int N, *es, **edge;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"graph_end");
      need[n] = d;
    }
    {
      string n = "graph_setEdge";
      string c = "  void setEdge(int N__, int M, int A[], int B[], void **mem = &wmem){\n    int i;\n    N = N__;\n    walloc1d(&es, N, mem);\n    walloc1d(&edge, N, mem);\n    rep(i,N) es[i] = 0;\n    rep(i,M) es[A[i]]++, es[B[i]]++;\n    rep(i,N) walloc1d(&edge[i], es[i], mem);\n    rep(i,N) es[i] = 0;\n    rep(i,M) edge[A[i]][es[A[i]]++] = B[i], edge[B[i]][es[B[i]]++] = A[i];\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
      parent[n] = "graph";
    }
    {
      string n = "graph_setDirectEdge";
      string c = "  void setDirectEdge(int N__, int M, int A[], int B[], void **mem = &wmem){\n    int i;\n    N = N__;\n    walloc1d(&es, N, mem);\n    walloc1d(&edge, N, mem);\n    walloc1d(&edge[0], M, mem);\n    rep(i,N) es[i] = 0;\n    rep(i,M) es[A[i]]++;\n    rep(i,N) walloc1d(&edge[i], es[i], mem);\n    rep(i,N) es[i] = 0;\n    rep(i,M) edge[A[i]][es[A[i]]++] = B[i];\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
      parent[n] = "graph";
    }
    {
      string n = "graph_setEdgeRootedTree";
      string c = "  void setEdgeRootedTree(int N__, int M, int A[], int B[], int root=0, int reorder=0, int cnv[] = NULL, void **mem = &wmem){\n    int i, j, k;\n    int *dist, *q, qs, qe, *ind;\n    void *tmem;\n    N = N__;\n    tmem = ((char*)(*mem)) + (sizeof(int) * N + 15) + (sizeof(int*) * N + 15) + (sizeof(int) * M + 15 * N);\n    walloc1d(&es, N, mem);\n    walloc1d(&edge, N, mem);\n    rep(i,N) es[i] = 0;\n    rep(i,M) es[A[i]]++, es[B[i]]++;\n    rep(i,N) walloc1d(&edge[i], es[i], &tmem);\n    rep(i,N) es[i] = 0;\n    rep(i,M) edge[A[i]][es[A[i]]++] = B[i], edge[B[i]][es[B[i]]++] = A[i];\n    walloc1d(&dist, N, &tmem);\n    walloc1d(&q, N, &tmem);\n    walloc1d(&ind, N, &tmem);\n    if(cnv==NULL) walloc1d(&cnv, N, &tmem);\n    rep(i,N) dist[i] = -1;\n    dist[root] = 0;\n    qs = qe = 0;\n    q[qe++] = root;\n    while(qs < qe){\n      i = q[qs++];\n      rep(j,es[i]){\n        k = edge[i][j];\n        if(dist[k]==-1) dist[k] = dist[i] + 1, q[qe++] = k;\n      }\n    }\n    if(reorder == 0){\n      rep(i,N) cnv[i] = i;\n      rep(i,N) ind[i] = i;\n    } else {\n      rep(i,N) cnv[i] = q[i];\n      rep(i,N) ind[cnv[i]] = i;\n    }\n    rep(i,N) es[i] = 0;\n    rep(i,M){\n      j = A[i];\n      k = B[i];\n      if(dist[j] > dist[k]) swap(j, k);\n      es[ind[j]]++;\n    }\n    rep(i,N) walloc1d(&edge[i], es[i], mem);\n    rep(i,N) es[i] = 0;\n    rep(i,M){\n      j = A[i];\n      k = B[i];\n      if(dist[j] > dist[k]) swap(j, k);\n      j = ind[j];\n      k = ind[k];\n      edge[j][es[j]++] = k;\n    }\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
      parent[n] = "graph";
    }
    {
      string n = "graph_reverse";
      string c = "  graph reverse(void **mem = &wmem){\n    int i, j, k;\n    graph g;\n    g.N = N;\n    walloc1d(&g.es, N, mem);\n    walloc1d(&g.edge, N, mem);\n    rep(i,N) g.es[i] = 0;\n    rep(i,N) rep(j,es[i]) g.es[edge[i][j]]++;\n    rep(i,N) walloc1d(&g.edge[i], g.es[i], mem);\n    rep(i,N) g.es[i] = 0;\n    rep(i,N) rep(j,es[i]){\n      k = edge[i][j];\n      g.edge[k][g.es[k]++] = i;\n    }\n    return g;\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
      parent[n] = "graph";
    }
    {
      string n = "graph_reduce";
      string c = "  graph reduce(int tn, int ind[], int self_e = 0, int dep_e = 0, void **mem = &wmem){\n    int i, j, k, M = 0;\n    int x, y;\n    graph g;\n    void *tmem;\n    pair<int,int> *A;\n    rep(i,N) M += es[i];\n    tmem = ((char*)(*mem)) + sizeof(int**) * N + sizeof(int*) * N + sizeof(int) * M + 16 * (N+2);\n    walloc1d(&A, M, &tmem);\n    M = 0;\n    rep(i,N){\n      x = ind[i];\n      if(x < 0) continue;\n      rep(j,es[i]){\n        y = ind[edge[i][j]];\n        if(y < 0) continue;\n        if(self_e==0 && x==y) continue;\n        A[M++] = make_pair(x, y);\n      }\n    }\n    if(dep_e==0){\n      sort(A, A+M);\n      k = 0;\n      rep(i,M){\n        if(k && A[k-1]==A[i]) continue;\n        A[k++] = A[i];\n      }\n      M = k;\n    }\n    g.N = tn;\n    walloc1d(&g.es, tn, mem);\n    walloc1d(&g.edge, tn, mem);\n    rep(i,tn) g.es[i] = 0;\n    rep(i,M) g.es[A[i].first]++;\n    rep(i,tn) walloc1d(&g.edge[i], g.es[i], mem);\n    rep(i,tn) g.es[i] = 0;\n    rep(i,M){\n      j = A[i].first;\n      k = A[i].second;\n      g.edge[j][g.es[j]++] = k;\n    }\n    return g;\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
      parent[n] = "graph";
    }
    {
      string n = "graph_getDist";
      string c = "  void getDist(int root, int res[], void *mem = wmem){\n    int i,j,k,*q,s,z;\n    walloc1d(&q, N, &mem);\n    rep(i,N)res[i]=-1;\n    res[root]=0;\n    s=0;\n    z=1;\n    q[0]=root;\n    while(z){\n      i=q[s++];\n      z--;\n      rep(j,es[i]){\n        k=edge[i][j];\n        if(res[k]>=0)continue;\n        res[k]=res[i]+1;\n        q[s+z++]=k;\n      }\n    }\n  }\n  int getDist(int a, int b, void *mem = wmem){\n    int i, j, k, *q, s, z, *d;\n    if(a==b) return 0;\n    walloc1d(&d, N, &mem);\n    walloc1d(&q, N, &mem);\n    rep(i,N) d[i] = -1;\n    d[a] = 0;\n    s = 0;\n    z = 1;\n    q[0] = a;\n    while(z){\n      i = q[s++];\n      z--;\n      rep(j,es[i]){\n        k = edge[i][j];\n        if(d[k] >= 0) continue;\n        d[k] = d[i] + 1;\n        if(k==b) return d[k];\n        q[s+z++] = k;\n      }\n    }\n    return -1;\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
      parent[n] = "graph";
    }
    {
      string n = "graph_TreeDiameter";
      string c = "  int TreeDiameter(int &a, int &b, void *mem = wmem){\n    int i, mx, *d;\n    walloc1d(&d, N, &mem);\n    getDist(0, d, mem);\n    mx = -1;\n    rep(i,N) if(mx < d[i]) mx = d[i], a = i;\n    getDist(a, d, mem);\n    mx = -1;\n    rep(i,N) if(mx < d[i]) mx = d[i], b = i;\n    return mx;\n  }\n  int TreeDiameter(void *mem = wmem){\n    int a, b;\n    return TreeDiameter(a, b, mem);\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"graph_getDist");
      need[n] = d;
      parent[n] = "graph";
    }
    {
      string n = "graph_getDistTree_WeightedNode_max";
      string c = "  template<class S1, class S2>\n  void getDistTree_WeightedNode_max(int root, S1 w[], S2 res[], void *mem = wmem){\n    int i, j, k, m;\n    int *q, qs = 0, qe = 1;\n    char *vis;\n    walloc1d(&q,N,&mem);\n    walloc1d(&vis,N,&mem);\n    rep(i,N) vis[i] = 0;\n    vis[root] = 1;\n    res[root] = w[root];\n    q[0] = root;\n    while(qs < qe){\n      i = q[qs++];\n      rep(j,es[i]){\n        k = edge[i][j];\n        if(vis[k]) continue;\n        res[k] = max(w[k], res[i]);\n        q[qe++] = k;\n        vis[k] = 1;\n      }\n    }\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"max_L");
      need[n] = d;
      parent[n] = "graph";
    }
    {
      string n = "graph_getDistPairMatrix";
      string c = "  template<class T>\n  void getDistPairMatrix(int k, int ind[], T **d, void *mem = wmem){\n    int *dist, i, j;\n    walloc1d(&dist, N, &mem);\n    if(k==0) d[0][0] = 0;\n    rep(i,k-1){\n      getDist(ind[i], dist, mem);\n      rep(j,i,k) d[i][j] = d[j][i] = dist[ind[j]];\n    }\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"graph_getDist");
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
      parent[n] = "graph";
    }
    {
      string n = "graph_SubTreeSize";
      string c = "  void SubTreeSize(int root, int res[], void *mem = wmem){\n    int i, j, k, m;\n    int *q, qs = 0, qe = 1;\n    walloc1d(&q,N,&mem);\n    rep(i,N) res[i] = -1;\n    res[root] = 0;\n    q[0] = root;\n    while(qs < qe){\n      i = q[qs++];\n      rep(j,es[i]){\n        k = edge[i][j];\n        if(res[k]==0) continue;\n        res[k] = 0;\n        q[qe++] = k;\n      }\n    }\n    rrep(m,N){\n      i = q[m];\n      res[i] = 1;\n      rep(j,es[i]){\n        k = edge[i][j];\n        res[i] += res[k];\n      }\n    }\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
      parent[n] = "graph";
    }
    {
      string n = "graph_SubTreeWeight";
      string c = "  template<class S>\n  void SubTreeWeight(int root, S weight[], S res[], void *mem = wmem){\n    int i, j, k, m;\n    int *q, qs = 0, qe = 1;\n    walloc1d(&q,N,&mem);\n    rep(i,N) res[i] = -1;\n    res[root] = 0;\n    q[0] = root;\n    while(qs < qe){\n      i = q[qs++];\n      rep(j,es[i]){\n        k = edge[i][j];\n        if(res[k]==0) continue;\n        res[k] = 0;\n        q[qe++] = k;\n      }\n    }\n    rrep(m,N){\n      i = q[m];\n      res[i] = weight[i];\n      rep(j,es[i]){\n        k = edge[i][j];\n        res[i] += res[k];\n      }\n    }\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
      parent[n] = "graph";
    }
    {
      string n = "graph_cntShortest";
      string c = "  template<class S> void cntShortest(int root, int dist[], S cnt[], void *mem = wmem){\n    int i,j,k,*q,s,z;\n    walloc1d(&q, N, &mem);\n    rep(i,N)dist[i]=-1;\n    rep(i,N)cnt[i]=0;\n    dist[root]=0;\n    cnt[root]=1;\n    s=0;\n    z=1;\n    q[0]=root;\n    while(z){\n      i=q[s++];\n      z--;\n      rep(j,es[i]){\n        k=edge[i][j];\n        if(dist[k]==-1) dist[k] = dist[i] + 1, cnt[k] = 0, q[s+z++] = k;\n        if(dist[k]==dist[i]+1) cnt[k] += cnt[i];\n      }\n    }\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
      parent[n] = "graph";
    }
    {
      string n = "graph_scc";
      string c = "  inline int sccDFS(int num[], int st, int mx){\n    int i,j;\n    num[st]=-2;\n    rep(i,es[st]) {\n      j=edge[st][i]; if(num[j]==-1) mx=sccDFS(num,j,mx);\n    }\n    num[st]=mx; return mx+1;\n  }\n  \n  int scc(int res[], void *mem = wmem){\n    int i, j, k, ret=0;\n    graph r;\n    int *st, st_size, *num, *nrv;\n    \n    r = reverse(&mem);\n    walloc1d(&st, N, &mem);\n    walloc1d(&num, N, &mem);\n    walloc1d(&nrv, N, &mem);\n    \n    rep(i,N) res[i] = num[i] = -1;\n    k = 0;\n    rep(i,N) if(num[i]==-1) k = sccDFS(num,i,k);\n    rep(i,N) nrv[num[i]] = i;\n    \n    for(k=N-1;k>=0;k--) {\n      i=nrv[k]; if(res[i]>=0)continue;\n      res[i]=ret; st_size=0; st[st_size++]=i;\n      while(st_size){\n        i=st[--st_size];\n        rep(j,r.es[i])\n          if(res[r.edge[i][j]]==-1) res[r.edge[i][j]]=ret, st[st_size++]=r.edge[i][j];\n      }\n      ret++;\n    }\n    \n    return ret;\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"graph_reverse");
      need[n] = d;
      parent[n] = "graph";
    }
    {
      string n = "graph_bcc";
      string c = "  inline void bccDFS(int v, int u, int *res, int *rt, int &rts, int *S, int &Ss, int *inS, int *num, int &tm){\n    int i, k;\n    \n    num[v] = ++tm;\n    S[Ss++] = v; inS[v] = 1;\n    rt[rts++] = v;\n    rep(i, es[v]){\n      int w = edge[v][i];\n      if(!num[w]){\n        bccDFS(w, v, res, rt, rts, S, Ss, inS, num, tm);\n      } else if(u != w && inS[w]){\n        while(num[rt[rts-1]] > num[w]) rts--;\n      }\n    }\n    \n    if(v == rt[rts-1]){\n      k = S[Ss-1];\n      for(;;){\n        int w = S[--Ss];\n        inS[w] = 0;\n        res[w] = k;\n        if(v==w) break;\n      }\n      rts--;\n    }\n  }\n  int bcc(int res[], void *mem=wmem){\n    int i, k;\n    int *rt, *S, *num, *inS;\n    pair<int,int> *arr;\n    int rts = 0, Ss = 0, tm = 0;\n    walloc1d(&num, N, &mem);\n    walloc1d(&rt, N, &mem);\n    walloc1d(&S, N, &mem);\n    walloc1d(&inS, N, &mem);\n    \n    memset(num, 0, sizeof(int)*N);\n    memset(inS, 0, sizeof(int)*N);\n    rep(i,N) if(!num[i]) bccDFS(i, N, res, rt, rts, S, Ss, inS, num, tm);\n    \n    arr = (pair<int,int>*)mem;\n    rep(i,N) arr[i].first = res[i], arr[i].second = i;\n    sort(arr, arr+N);\n    k = 0;\n    rep(i,N){\n      if(i && arr[i].first != arr[i-1].first) k++;\n      res[arr[i].second] = k;\n    }\n    return k+1;\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
      parent[n] = "graph";
    }
    {
      string n = "graph_articulation";
      string c = "  void articulationDFS(int n, int b, int &k, int od[], int lw[], int vs[], int &ress, int res[]){\n    int i, j, a = 0, c = 0;\n    vs[n] = 1;\n    lw[n] = od[n] = k++;\n    rep(i,es[n]){\n      j = edge[n][i];\n      if(j==b) continue;\n      \n      if(!vs[j]){\n        c++;\n        articulationDFS(j, n, k, od, lw, vs, ress, res);\n        lw[n] <?= lw[j];\n        if(b != -1 && od[n] <= lw[j]) a = 1;\n      } else {\n        lw[n] <?= od[j];\n      }\n    }\n    if(b == -1 && c >= 2) a = 1;\n    if(a) res[ress++] = n;\n  }\n  \n  int articulation(int res[], void *mem=wmem){\n    int i, k = 0, ress = 0;\n    int *od, *lw, *vs;\n    walloc1d(&od, N, &mem);\n    walloc1d(&lw, N, &mem);\n    walloc1d(&vs, N, &mem);\n    rep(i,N) vs[i] = 0;\n    rep(i,N) if(!vs[i]) articulationDFS(i, -1, k, od, lw, vs, ress, res);\n    return ress;\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"chmin");
      need[n] = d;
      parent[n] = "graph";
    }
    {
      string n = "graph_shortestPath";
      string c = "  int shortestPath(const int s, const int t, int res[], void *mem=wmem){\n    int i, j, k;\n    int *q, qs = 0, qe = 0, *b;\n    walloc1d(&b, N, &mem);\n    walloc1d(&q, N, &mem);\n    rep(i,N) b[i] = -1;\n    b[s] = -2;\n    q[qe++] = s;\n    while(qe > qs){\n      i = q[qs++];\n      rep(j,es[i]){\n        k = edge[i][j];\n        if(b[k]!=-1) continue;\n        b[k] = i;\n        q[qe++] = k;\n      }\n      if(b[t]!=-1) break;\n    }\n    if(b[t]==-1) return -1;\n    k = 0;\n    res[k] = i = t;\n    while(i != s) res[++k] = (i = b[i]);\n    std::reverse(res, res+k+1);\n    return k;\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
      parent[n] = "graph";
    }
    {
      string n = "graph_TopologicalSort";
      string c = "  int TopologicalSort(int res[], void *mem=wmem){\n    int i, j, k, rs;\n    int *deg, *q, qs = 0, qe = 0;\n    walloc1d(&deg, N, &mem);\n    walloc1d(&q, N, &mem);\n    rs = 0;\n    rep(i,N) deg[i] = 0;\n    rep(i,N) rep(j,es[i]) deg[edge[i][j]]++;\n    rep(i,N) if(deg[i]==0) q[qe++] = i;\n    while(qs < qe){\n      i = q[qs++];\n      res[rs++] = i;\n      rep(j,es[i]){\n        k = edge[i][j];\n        deg[k]--;\n        if(deg[k]==0) q[qe++] = k;\n      }\n    }\n    if(rs==N) return 1;\n    return 0;\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
      parent[n] = "graph";
    }
    {
      string n = "graph_longestPath_length";
      string c = "  int longestPath_length(void *mem = wmem){\n    int *dp, *arr, res = 0, i, j, ii, jj;\n    walloc1d(&dp, N, &mem);\n    walloc1d(&arr, N, &mem);\n    if(!TopologicalSort(arr, mem)) return -1;\n    rep(i,N) dp[i] = 0;\n    rep(ii,N){\n      i = arr[ii];\n      res >?= dp[i];\n      rep(jj,es[i]){\n        j = edge[i][jj];\n        dp[j] >?= dp[i] + 1;\n      }\n    }\n    return res;\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"chmax");
      d.push_back((string)"graph_TopologicalSort");
      need[n] = d;
      parent[n] = "graph";
    }
    {
      string n = "graph_Grundy";
      string c = "  int Grundy(int res[], void *mem=wmem){\n    int i, j, k, fg;\n    int *tp, *arr;\n    walloc1d(&tp, N, &mem);\n    fg = TopologicalSort(tp, mem);\n    if(fg==0) return 0;\n    walloc1d(&arr, N, &mem);\n    rep(i,N) arr[i] = -1;\n    rrep[tp](i,N){\n      rep[edge[i]](j,es[i]) arr[res[j]] = i;\n      for(k=0;;k++) if(arr[k] != i) break;\n      res[i] = k;\n    }\n    return 1;\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"graph_TopologicalSort");
      need[n] = d;
      parent[n] = "graph";
    }
    {
      string n = "graph_preorder";
      string c = "  int preorder(int res[], int root = 0, void *mem=wmem){\n    int i, j, k, sts, sz = 0;\n    ll *st;\n    char *vis;\n    walloc1d(&vis, N, &mem);\n    walloc1d(&st, N, &mem);\n    sts = 0;\n    st[sts++] = ((ll)root) << 32;\n    rep(i,N) vis[i] = 0;\n    vis[root] = 1;\n    while(sts){\n      i = st[--sts] >> 32;\n      j = st[sts] & 2147483647;\n      if(j==0) res[sz++] = i;\n      while(j < es[i]){\n        k = edge[i][j++];\n        if(vis[k]) continue;\n        if(j < es[i]) st[sts++] = (((ll)i) << 32) + j;\n        vis[k] = 1;\n        st[sts++] = ((ll)k) << 32;\n        break;\n      }\n    }\n    return sz;\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
      parent[n] = "graph";
    }
    {
      string n = "graph_anUndirectedCycle";
      string c = "  int anUndirectedCycle(int res[] = NULL, void *mem = wmem){\n    int i, j, k, m;\n    int *arr, *q, qs, qe, *bk;\n    if(res==NULL) walloc1d(&res, N+1, &mem);\n    rep(i,N) rep(j,es[i]) if(edge[i][j]==i){\n      res[0] = res[1] = i;\n      return 1;\n    }\n    walloc1d(&arr, N, &mem);\n    walloc1d(&q, N, &mem);\n    walloc1d(&bk, N, &mem);\n    rep(i,N) arr[i] = -1;\n    rep(i,N) rep(j,es[i]){\n      k = edge[i][j];\n      if(arr[k] == i){\n        res[0] = i;\n        res[1] = k;\n        res[2] = i;\n        return 2;\n      }\n      arr[k] = i;\n    }\n    rep(i,N) arr[i] = bk[i] = -1;\n    rep(m,N) if(arr[m]==-1){\n      qs = qe = 0;\n      q[qe++] = m;\n      arr[m] = 0;\n      while(qs < qe){\n        i = q[qs++];\n        rep(j,es[i]){\n          k = edge[i][j];\n          if(arr[k]==-1){\n            arr[k] = arr[i] + 1;\n            bk[k] = i;\n            q[qe++] = k;\n            continue;\n          }\n          if(arr[k] == arr[i] - 1) continue;\n          \n          qs = qe = 1;\n          res[0] = i;\n          q[0] = k;\n          while(i!=k){\n            if(arr[i] > arr[k]){\n              res[qs++] = (i = bk[i]);\n            } else {\n              q[qe++] = (k = bk[k]);\n            }\n          }\n          reverse(res, res+qs);\n          rep(i,qe) res[qs++] = q[i];\n          return qs - 1;\n        }\n      }\n    }\n    return -1;\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
      parent[n] = "graph";
    }
    {
      string n = "graph_shortestCycle";
      string c = "  int shortestCycle(int res[] = NULL, void *mem=wmem){\n    int i, j, k, r, sz, tp;\n    int *q, qs, qe, *arr, *bk;\n    if(res==NULL) walloc1d(&res, N+1, &mem);\n    rep(i,N) rep(j,es[i]){\n      k = edge[i][j];\n      if(k==i){\n        res[0] = res[1] = i;\n        return 1;\n      }\n    }\n    walloc1d(&arr, N, &mem);\n    walloc1d(&bk, N, &mem);\n    walloc1d(&q, N, &mem);\n    sz = int_inf;\n    rep(r,N){\n      rep(i,N) arr[i] = -1;\n      rep(i,N) bk[i] = -1;\n      arr[r] = 0;\n      qs = qe = 0;\n      q[qe++] = r;\n      tp = int_inf;\n      while(qs < qe && tp == int_inf){\n        i = q[qs++];\n        rep(j,es[i]){\n          k = edge[i][j];\n          if(k==r){\n            bk[k] = i;\n            tp = arr[i] + 1;\n            break;\n          }\n          if(arr[k]==-1){\n            bk[k] = i;\n            arr[k] = arr[i] + 1;\n            q[qe++] = k;\n            continue;\n          }\n        }\n      }\n      if(sz > tp){\n        sz = tp;\n        res[tp] = r;\n        k = r;\n        for(;;){\n          k = bk[k];\n          res[arr[k]] = k;\n          if(k==r) break;\n        }\n      }\n    }\n    \n    if(sz==int_inf) sz = -1;\n    return sz;\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
      parent[n] = "graph";
    }
    {
      string n = "graph_shortestUndirectedCycle_length";
      string c = "int shortestUndirectedCycle_length(void *mem=wmem){\n    int i, j, k, r, res;\n    int *arr, *q, qs, qe;\n    rep(i,N) rep(j,es[i]) if(edge[i][j]==i) return 1;\n    walloc1d(&arr, N, &mem);\n    rep(i,N) arr[i] = -1;\n    rep(i,N) rep(j,es[i]){\n      k = edge[i][j];\n      if(arr[k]==i) return 2;\n      arr[k] = i;\n    }\n    walloc1d(&q, N, &mem);\n    res = int_inf;\n    rep(r,N){\n      rep(i,N) arr[i] = -1;\n      arr[r] = 0;\n      qs = qe = 0;\n      q[qe++] = r;\n      while(qs < qe){\n        i = q[qs++];\n        rep(j,es[i]){\n          k = edge[i][j];\n          if(arr[k]==-1){\n            arr[k] = arr[i] + 1;\n            q[qe++] = k;\n            continue;\n          }\n          if(arr[k]==arr[i]) res <?= 2 arr[i] + 1;\n          if(arr[k]==arr[i]+1) res <?= 2 arr[k];\n        }\n      }\n    }\n    \n    if(res==int_inf) res = -1;\n    return res;\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"chmin");
      need[n] = d;
      parent[n] = "graph";
    }
    {
      string n = "graph_shortestUndirectedCycle_length_node";
      string c = "  int shortestUndirectedCycle_length(const int node, void *mem=wmem){\n    int i, j, k, res;\n    int *arr, *q, qs, qe;\n    const int r = node;\n    rep(j,es[r]) if(edge[r][j]==r) return 1;\n    walloc1d(&arr, N, &mem);\n    rep(i,N) arr[i] = 0;\n    rep(j,es[r]){\n      k = edge[r][j];\n      if(arr[k]) return 2;\n      arr[k] = 1;\n    }\n    walloc1d(&q, N, &mem);\n    res = int_inf;\n    rep(i,N) arr[i] = -1;\n    arr[r] = 0;\n    qs = qe = 0;\n    q[qe++] = r;\n    while(qs < qe){\n      i = q[qs++];\n      rep(j,es[i]){\n        k = edge[i][j];\n        if(arr[k]==-1){\n          arr[k] = arr[i] + 1;\n          q[qe++] = k;\n          continue;\n        }\n        if(arr[k]==arr[i]) res <?= 2 arr[i] + 1;\n        if(arr[k]==arr[i]+1) res <?= 2 arr[k];\n      }\n    }\n    \n    if(res==int_inf) res = -1;\n    return res;\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"chmin");
      need[n] = d;
      parent[n] = "graph";
    }
    {
      string n = "graph_maxIndependenceSet";
      string c = "  int maxIndependenceSet(int res[] = NULL, void *mem = wmem, int lim = -1, int skip = 0){\n    int i, j, k, m, x, y, fg;\n    int ress;\n    int *deg, *used, *ind, *rev, *val;\n    ll *pr;\n    graph g;\n    unionFind uf;\n    void *tmem;\n    if(res==NULL) walloc1d(&res, N, &mem);\n    if(N == 0) return 0;\n    if(N == 1) res[0] = 0, return 1;\n    if(N <= lim) return 0;\n    walloc1d(&deg, N, &mem);\n    walloc1d(&used, N, &mem);\n    walloc1d(&ind, N, &mem);\n    walloc1d(&rev, N, &mem);\n    walloc1d(&val, N, &mem);\n    walloc1d(&pr, N, &mem);\n    rep(i,N) deg[i] = es[i];\n    rep(i,N) used[i] = 0;\n    ress = 0;\n    if(!(skip&1)){\n      do{\n        fg = 0;\n        rep(i,N) if(!used[i] && deg[i] <= 1){\n          fg = 1;\n          res[ress++] = i;\n          used[i] = 1;\n          if(deg[i] == 0) continue;\n          rep(j,es[i]){\n            k = edge[i][j];\n            if(!used[k]) break;\n          }\n          deg[k]--;\n          used[k] = 1;\n          rep(j,es[k]){\n            m = edge[k][j];\n            if(!used[m]) deg[m]--;\n          }\n        }\n        \n        rep(i,N) if(!used[i] && deg[i]==2){\n          m = 0;\n          rep(j,es[i]){\n            k = edge[i][j];\n            if(used[k]) continue;\n            if(m==0) x = k, m++;\n            else     y = k, m++;\n          }\n          rep(j,es[x]) if(edge[x][j] == y) break;\n          if(j < es[x]){\n            fg = 1;\n            used[i] = used[x] = used[y] = 1;\n            res[ress++] = i;\n            rep(j,es[x]){\n              m = edge[x][j];\n              if(!used[m]) deg[m]--;\n            }\n            rep(j,es[y]){\n              m = edge[y][j];\n              if(!used[m]) deg[m]--;\n            }\n          }\n        }\n      }while(fg);\n    }\n    if(ress){\n      k = 0;\n      rep(i,N){\n        if(used[i]) ind[i] = -1, continue;\n        ind[i] = k;\n        rev[k] = i;\n        k++;\n      }\n      g = reduce(k, ind, 1, 1, &mem);\n      m = g.maxIndependenceSet(res+ress, mem, lim - ress, 1);\n      rep(i,m) res[ress+i] = rev[res[ress+i]];\n      ress += m;\n      sort(res, res+ress);\n      return ress;\n    }\n    if(N-2 <= lim) return 0;\n    if(lim >= 1){\n      rep(i,N) deg[i] = es[i];\n      sort(deg, deg+N);\n      int ss = 0;\n      rep(k,N) ss += es[k];\n      j = 0;\n      rep(i,N){\n        j += deg[i];\n        if(2*j > ss) break;\n      }\n      if(i <= lim) return 0;\n      \n      rep(i,N) if(deg[i] > N-i-1) break;\n      if(i <= lim) return 0;\n    }\n    \n    if(!(skip&2)){\n      k = 0;\n      uf.walloc(N, &mem);\n      uf.init(N);\n      rep(i,N) rep(j,es[i]) k += uf(i,edge[i][j]);\n      if(k < N-1){\n        rep(i,N) val[i] = i;\n        rep(i,N) deg[i] = 0;\n        rep(i,N) deg[uf(i)]++;\n        rsortA(N, deg, val);\n        y = N;\n        rrep(x,N) if(deg[x]){\n          k = 0;\n          rep(i,N){\n            if(uf(i)!=val[x]) ind[i] = -1, continue;\n            ind[i] = k;\n            rev[k] = i;\n            k++;\n          }\n          y -= k;\n          tmem = mem;\n          g = reduce(k, ind, 1, 1, &mem);\n          m = g.maxIndependenceSet(res+ress, mem, max(-1, lim-(y-2x)-ress), 3);\n          mem = tmem;\n          rep(i,m) res[ress+i] = rev[res[ress+i]];\n          ress += m;\n        }\n        sort(res,res+ress);\n        return ress;\n      }\n    }\n    k = articulation(ind, mem);\n    rep(i,N) pr[i] = 0;\n    rep(i,k) pr[ind[i]] += (1LL<<40);\n    rep(i,N) rep(j,es[i]) pr[i] += (1LL<<20) - es[edge[i][j]];\n    x = argmax(pr(N));\n    k = 0;\n    rep(i,N){\n      if(i == x) ind[i] = -1, continue;\n      ind[i] = k;\n      rev[k] = i;\n      k++;\n    }\n    tmem = mem;\n    g = reduce(k, ind, 1, 1, &mem);\n    ress = g.maxIndependenceSet(res, mem, lim);\n    mem = tmem;\n    rep(i,ress) res[i] = rev[res[i]];\n    \n    k = 0;\n    used[x] = 1;\n    rep(j,es[x]) used[edge[x][j]] = 1;\n    rep(i,N){\n      if(used[i]) ind[i] = -1, continue;\n      ind[i] = k;\n      rev[k] = i;\n      k++;\n    }\n    g = reduce(k, ind, 1, 1, &mem);\n    m = g.maxIndependenceSet(deg, mem, max(ress, lim)-1);\n    if(m+1 > ress){\n      rep(i,m) res[i] = rev[deg[i]];\n      res[m++] = x;\n      ress = m;\n    }\n    \n    sort(res,res+ress);\n    return ress;\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"graph_reduce");
      d.push_back((string)"graph_articulation");
      d.push_back((string)"max_L");
      d.push_back((string)"sortA");
      d.push_back((string)"rsortA");
      d.push_back((string)"unionFind");
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
      parent[n] = "graph";
    }
    {
      string n = "graph_countIndependenceSet";
      string c = "  ll countIndependenceSet(void *mem = wmem, int skip = 0){\n    int i, j, k, x;\n    ll res;\n    int *deg, *ind;\n    ll *pr;\n    graph g;\n    unionFind uf;\n    void *tmem;\n    if(N == 0) return 1;\n    if(N == 1) return 2;\n    walloc1d(&deg, N, &mem);\n    walloc1d(&ind, N, &mem);\n    walloc1d(&pr, N, &mem);\n    if(!skip){\n      k = 0;\n      uf.walloc(N, &mem);\n      uf.init(N);\n      rep(i,N) rep(j,es[i]) k += uf(i,edge[i][j]);\n      if(k < N-1){\n        res = 1;\n        rep(i,N) deg[i] = 0;\n        rep(i,N) deg[uf(i)]++;\n        rep(x,N) if(deg[x]){\n          k = 0;\n          rep(i,N){\n            if(uf(i)!=x) ind[i] = -1, continue;\n            ind[i] = k++;\n          }\n          tmem = mem;\n          g = reduce(k, ind, 1, 1, &mem);\n          res *= g.countIndependenceSet(mem, 1);\n          mem = tmem;\n        }\n        return res;\n      }\n    }\n    k = articulation(ind, mem);\n    rep(i,N) pr[i] = 0;\n    rep(i,k) pr[ind[i]] += (1LL<<40);\n    rep(i,N) rep(j,es[i]) pr[i] += (1LL<<20) - es[edge[i][j]];\n    x = argmax(pr(N));\n    res = 0;\n    \n    k = 0;\n    rep(i,N){\n      if(i == x) ind[i] = -1, continue;\n      ind[i] = k++;\n    }\n    tmem = mem;\n    g = reduce(k, ind, 1, 1, &mem);\n    res += g.countIndependenceSet(mem);\n    mem = tmem;\n    k = 0;\n    rep(i,N) deg[i] = 0;\n    deg[x] = 1;\n    rep(j,es[x]) deg[edge[x][j]] = 1;\n    rep(i,N){\n      if(deg[i]) ind[i] = -1, continue;\n      ind[i] = k++;\n    }\n    g = reduce(k, ind, 1, 1, &mem);\n    res += g.countIndependenceSet(mem);\n    return res;\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"graph_reduce");
      d.push_back((string)"graph_articulation");
      d.push_back((string)"max_L");
      d.push_back((string)"unionFind");
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
      parent[n] = "graph";
    }
    {
      string n = "graph_bipartite";
      string c = "  int bipartite(int res[] = NULL, void *mem = wmem){\n    int i, j, k, x;\n    int *st, sts;\n    if(res==NULL) walloc1d(&res, N, &mem);\n    walloc1d(&st, N, &mem);\n    rep(i,N) res[i] = -1;\n    rep(x,N) if(res[x]==-1){\n      res[x] = 0;\n      sts = 0;\n      st[sts++] = x;\n      while(sts){\n        i = st[--sts];\n        rep(j,es[i]){\n          k = edge[i][j];\n          if(res[k]==-1){\n            res[k] = 1 - res[i];\n            st[sts++] = k;\n          }\n          if(res[i] + res[k] != 1) return 0;\n        }\n      }\n    }\n    return 1;\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
      parent[n] = "graph";
    }
    {
      string n = "graph_Rerooting_head";
      string c = "template<class V> void Rerooting(V res[], void *mem = wmem);\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "graph";
    }
    {
      string n = "graph_Rerooting";
      string c = "template<class V>\nvoid graph::Rerooting(V res[], void *mem /*= wmem*/){\n  int i, j, k, ui, ii;\n  int *upind, *dwind, *vis, *q, qs, qe;\n  V **val, **pleft, **pright, tmp;\n  \n  walloc1d(&upind, N, &mem);\n  walloc1d(&dwind, N, &mem);\n  walloc1d(&vis, N, &mem);\n  walloc1d(&q, N, &mem);\n  walloc1d(&val, N, &mem);\n  rep(i,N) walloc1d(&val[i], es[i], &mem);\n  walloc1d(&pleft, N, &mem);\n  rep(i,N) walloc1d(&pleft[i], es[i]+1, &mem);\n  walloc1d(&pright, N, &mem);\n  rep(i,N) walloc1d(&pright[i], es[i]+1, &mem);\n  upind[0] = -1;\n  rep(i,N) vis[i] = 0;\n  qs = qe = 0;\n  q[qe++] = 0; vis[0] = 1;\n  while(qs < qe){\n    i = q[qs++];\n    rep(j,es[i]){\n      k = edge[i][j];\n      if(vis[k]){\n        upind[i] = j;\n        continue;\n      }\n      dwind[k] = j;\n      vis[k] = 1;\n      q[qe++] = k;\n    }\n  }\n  rrep(ii,N){\n    i = q[ii];\n    rep(j,es[i]) if(j != upind[i]){\n      k = edge[i][j];\n      ui = upind[k];\n      tmp = RerootingMerge(pleft[k][ui], pright[k][es[k]-1-ui]);\n      tmp = RerootingNode(tmp, k);\n      val[i][j] = RerootingEdge(tmp, i, k);\n    }\n    RerootingId(pleft[i][0]);\n    RerootingId(pright[i][0]);\n    rep(j,es[i]){\n      if(j == upind[i]) break;\n      pleft[i][j+1] = RerootingMerge(pleft[i][j], val[i][j]);\n    }\n    rep(j,es[i]){\n      if(es[i]-1-j == upind[i]) break;\n      pright[i][j+1] = RerootingMerge(pright[i][j], val[i][es[i]-1-j]);\n    }\n  }\n  rep(ii,N){\n    i = q[ii];\n    j = upind[i];\n    if(j != -1){\n      k = edge[i][j];\n      ui = dwind[i];\n      tmp = RerootingMerge(pleft[k][ui], pright[k][es[k]-1-ui]);\n      tmp = RerootingNode(tmp, k);\n      val[i][j] = RerootingEdge(tmp, i, k);\n      rep(j,upind[i],es[i]){\n        pleft[i][j+1] = RerootingMerge(pleft[i][j], val[i][j]);\n      }\n      rep(j,es[i]-upind[i]-1,es[i]){\n        pright[i][j+1] = RerootingMerge(pright[i][j], val[i][es[i]-1-j]);\n      }\n    }\n    res[i] = RerootingNode(pleft[i][es[i]], i);\n  }\n}\n";
      string p = "last";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"graph_Rerooting_head");
      need[n] = d;
      parent[n] = "graph";
    }
    {
      string n = "graph_end";
      string c = "};\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "graph";
    }
    {
      string n = "wgraph";
      string c = "template<class T>\nstruct wgraph{\n  int N, *es, **edge;\n  T **cost;\n  graph g;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"wgraph_end");
      d.push_back((string)"graph");
      need[n] = d;
    }
    {
      string n = "wgraph_setEdge";
      string c = "  void setEdge(int N__, int M, int A[], int B[], T C[], void **mem = &wmem){\n    int i;\n    N = N__;\n    walloc1d(&es, N, mem);\n    rep(i,N) es[i] = 0;\n    rep(i,M) es[A[i]]++, es[B[i]]++;\n    walloc1d(&edge, N, mem);\n    rep(i,N) walloc1d(&edge[i], es[i], mem);\n    walloc1d(&cost, N, mem);\n    rep(i,N) walloc1d(&cost[i], es[i], mem);\n    rep(i,N) es[i] = 0;\n    rep(i,M){\n      edge[A[i]][es[A[i]]] = B[i];\n      edge[B[i]][es[B[i]]] = A[i];\n      cost[A[i]][es[A[i]]++] = C[i];\n      cost[B[i]][es[B[i]]++] = C[i];\n    }\n    g.N = N;\n    g.es = es;\n    g.edge = edge;\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
      parent[n] = "wgraph";
    }
    {
      string n = "wgraph_setDirectEdge";
      string c = "  void setDirectEdge(int N__, int M, int A[], int B[], T C[], void **mem = &wmem){\n    int i;\n    N = N__;\n    walloc1d(&es, N, mem);\n    rep(i,N) es[i] = 0;\n    rep(i,M) es[A[i]]++;\n    walloc1d(&edge, N, mem);\n    rep(i,N) walloc1d(&edge[i], es[i], mem);\n    walloc1d(&cost, N, mem);\n    rep(i,N) walloc1d(&cost[i], es[i], mem);\n    rep(i,N) es[i] = 0;\n    rep(i,M){\n      edge[A[i]][es[A[i]]] = B[i];\n      cost[A[i]][es[A[i]]++] = C[i];\n    }\n    g.N = N;\n    g.es = es;\n    g.edge = edge;\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
      parent[n] = "wgraph";
    }
    {
      string n = "wgraph_setEdgeRootedTree";
      string c = "void setEdgeRootedTree(int N__, int M, int A[], int B[], T C[], int root=0, int reorder=0, int cnv[] = NULL, void **mem = &wmem){\n    int i, j, k;\n    int *dist, *q, qs, qe, *ind;\n    void *tmem;\n    N = N__;\n    tmem = ((char*)(*mem)) + (sizeof(int) * N + 15) + (sizeof(int*) * N + 15) + (sizeof(int) * M + 15 * N) + (sizeof(T*) * N + 15) + (sizeof(T) * M + 15 * N);\n    walloc1d(&es, N, mem);\n    walloc1d(&edge, N, mem);\n    walloc1d(&cost, N, mem);\n    rep(i,N) es[i] = 0;\n    rep(i,M) es[A[i]]++, es[B[i]]++;\n    rep(i,N) walloc1d(&edge[i], es[i], &tmem);\n    rep(i,N) es[i] = 0;\n    rep(i,M) edge[A[i]][es[A[i]]++] = B[i], edge[B[i]][es[B[i]]++] = A[i];\n    walloc1d(&dist, N, &tmem);\n    walloc1d(&q, N, &tmem);\n    walloc1d(&ind, N, &tmem);\n    if(cnv==NULL) walloc1d(&cnv, N, &tmem);\n    rep(i,N) dist[i] = -1;\n    dist[root] = 0;\n    qs = qe = 0;\n    q[qe++] = root;\n    while(qs < qe){\n      i = q[qs++];\n      rep(j,es[i]){\n        k = edge[i][j];\n        if(dist[k]==-1) dist[k] = dist[i] + 1, q[qe++] = k;\n      }\n    }\n    if(reorder == 0){\n      rep(i,N) cnv[i] = i;\n      rep(i,N) ind[i] = i;\n    } else {\n      rep(i,N) cnv[i] = q[i];\n      rep(i,N) ind[cnv[i]] = i;\n    }\n    rep(i,N) es[i] = 0;\n    rep(i,M){\n      j = A[i];\n      k = B[i];\n      if(dist[j] > dist[k]) swap(j, k);\n      es[ind[j]]++;\n    }\n    rep(i,N) walloc1d(&edge[i], es[i], mem);\n    rep(i,N) walloc1d(&cost[i], es[i], mem);\n    rep(i,N) es[i] = 0;\n    rep(i,M){\n      j = A[i];\n      k = B[i];\n      if(dist[j] > dist[k]) swap(j, k);\n      j = ind[j];\n      k = ind[k];\n      edge[j][es[j]] = k;\n      cost[j][es[j]] = C[i];\n      es[j]++;\n    }\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
      parent[n] = "wgraph";
    }
    {
      string n = "wgraph_getDist";
      string c = "  template<class S>\n  void getDist(int root, S res[], S unreachable = -1, void *mem = wmem){\n    int i, j;\n    DijkstraHeap<S> hp;\n    hp.walloc(N, &mem);\n    hp.init(N);\n    hp.change(root,0);\n    while(hp.size){\n      i = hp.pop();\n      rep(j,es[i]) hp.change(edge[i][j], hp.val[i]+cost[i][j]);\n    }\n    rep(i,N) res[i] = (hp.visited[i] ? hp.val[i] : unreachable);\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"DijkstraHeap");
      need[n] = d;
      parent[n] = "wgraph";
    }
    {
      string n = "wgraph_getDistT";
      string c = "  template<class S>\n  S getDistT(int a, int b, S unreachable = -1, void *mem = wmem){\n    int i, j;\n    DijkstraHeap<S> hp;\n    hp.walloc(N, &mem);\n    hp.init(N);\n    hp.change(a,0);\n    while(hp.size){\n      i = hp.pop();\n      if(i==b) return hp.val[i];\n      rep(j,es[i]) hp.change(edge[i][j], hp.val[i]+cost[i][j]);\n    }\n    return unreachable;\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"DijkstraHeap");
      need[n] = d;
      parent[n] = "wgraph";
    }
    {
      string n = "wgraph_getDistDense";
      string c = "  template<class S>\n  void getDistDense(int root, S res[], S unreachable = -1, void *mem = wmem){\n    int i, j, k;\n    char *vis, *done;\n    walloc1d(&vis, N, &mem);\n    walloc1d(&done, N, &mem);\n    rep(i,N) vis[i] = 0;\n    rep(i,N) done[i] = 0;\n    res[root] = 0;\n    vis[root] = 1;\n    for(;;){\n      i = -1;\n      rep(j,N) if(!done[j] && vis[j]){\n        if(i==-1) i = j, continue;\n        if(res[i] > res[j]) i = j;\n      }\n      if(i==-1) break;\n      done[i] = 1;\n      rep(j,es[i]){\n        k = edge[i][j];\n        if(vis[k]==0 || res[k] > res[i] + cost[i][j]){\n          vis[k] = 1;\n          res[k] = res[i] + cost[i][j];\n        }\n      }\n    }\n    rep(i,N) if(!vis[i]) res[i] = unreachable;\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
      parent[n] = "wgraph";
    }
    {
      string n = "wgraph_getDistForest";
      string c = "  template<class S>\n  void getDistForest(int root, S res[], S unreachable = -1, void *mem = wmem){\n    int i,j,k,*q,s,z;\n    char *r;\n    walloc1d(&q,N,&mem);\n    walloc1d(&r,N,&mem);\n    rep(i,N)r[i]=0;\n    res[root]=0; r[root]=1;\n    s=0;\n    z=1;\n    q[0]=root;\n    while(z){\n      i=q[s++];\n      z--;\n      rep(j,es[i]){\n        k=edge[i][j];\n        if(r[k])continue;\n        res[k]=res[i]+cost[i][j];\n        r[k]=1;\n        q[s+z++]=k;\n      }\n    }\n    rep(i,N)if(!r[i])res[i]=unreachable;\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
      parent[n] = "wgraph";
    }
    {
      string n = "wgraph_BellmanFord";
      string c = "  template<class S>\n  void BellmanFord(int root, S res[], S unreachable = -1, S minusInf = -2, int step = -1, void *mem = wmem){\n    int i, j, k, t;\n    int *inf, *q, qs, qe;\n    S *arr;\n    walloc1d(&q, N, &mem);\n    walloc1d(&inf, N, &mem);\n    walloc1d(&arr, N, &mem);\n    rep(i,N) inf[i] = 0;\n    rep(i,N) res[i] = arr[i] = std::numeric_limits<S>::max();\n    res[root] = arr[root] = 0;\n    t = step;\n    if(t==-1) t = N;\n    rep(t){\n      rep(i,N) if(res[i] != std::numeric_limits<S>::max()) rep(j,es[i]){\n        arr[edge[i][j]] <?= res[i] + cost[i][j];\n      }\n      rep(i,N) res[i] = arr[i];\n    }\n    if(step != -1){\n      rep(i,N) if(res[i]==std::numeric_limits<S>::max()) res[i] = unreachable;\n      return;\n    }\n    rep(i,N) if(res[i] != std::numeric_limits<S>::max()) rep(j,es[i]){\n      k = edge[i][j];\n      if(arr[k] > res[i] + cost[i][j]) inf[k] = 1;\n    }\n    qs = qe = 0;\n    rep(i,N) if(inf[i]) q[qe++] = i;\n    while(qs < qe){\n      i = q[qs++];\n      rep(j,es[i]){\n        k = edge[i][j];\n        if(inf[k]==0){\n          inf[k] = 1;\n          q[qe++] = k;\n        }\n      }\n    }\n    rep(i,N) if(res[i]==std::numeric_limits<S>::max()) res[i] = unreachable;\n    rep(i,N) if(inf[i]==1) res[i] = minusInf;\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"chmin");
      need[n] = d;
      parent[n] = "wgraph";
    }
    {
      string n = "wgraph_getDist01";
      string c = "  template<class S>\n  void getDist01(int root, S res[], void *mem = wmem){\n    int i, j, k, nd;\n    int *q, qs, qe;\n    char *don;\n    walloc1d(&q, 2*N, &mem);\n    walloc1d(&don, N, &mem);\n    rep(i,N) res[i] = -1;\n    rep(i,N) don[i] = 0;\n    qs = qe = N;\n    q[qe++] = root;\n    res[root] = 0;\n    while(qs < qe){\n      i = q[qs++];\n      if(don[i]) continue;\n      don[i] = 1;\n      rep(j,es[i]){\n        k = edge[i][j];\n        nd = res[i] + cost[i][j];\n        if(res[k]==-1 || res[k] > nd){\n          res[k] = nd;\n          if(cost[i][j]==0){\n            q[--qs] = k;\n          } else {\n            q[qe++] = k;\n          }\n        }\n      }\n    }\n  }\n  int getDist01(int a, int b, void *mem = wmem){\n    int i, j, k, nd;\n    int *res, *q, qs, qe;\n    char *don;\n    walloc1d(&res, N, &mem);\n    walloc1d(&q, 2*N, &mem);\n    walloc1d(&don, N, &mem);\n    rep(i,N) res[i] = -1;\n    rep(i,N) don[i] = 0;\n    qs = qe = N;\n    q[qe++] = a;\n    res[a] = 0;\n    while(qs < qe){\n      i = q[qs++];\n      if(i==b) return res[b];\n      if(don[i]) continue;\n      don[i] = 1;\n      rep(j,es[i]){\n        k = edge[i][j];\n        nd = res[i] + cost[i][j];\n        if(res[k]==-1 || res[k] > nd){\n          res[k] = nd;\n          if(cost[i][j]==0){\n            q[--qs] = k;\n          } else {\n            q[qe++] = k;\n          }\n        }\n      }\n    }\n    return -1;\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
      parent[n] = "wgraph";
    }
    {
      string n = "wgraph_MST_Prim_cost";
      string c = "  T MST_Prim_cost(void *mem = wmem){\n    int i, j;\n    T res = 0;\n    DijkstraHeap<T> hp;\n    hp.walloc(N, &mem);\n    hp.init(N);\n    hp.change(0,0);\n    while(hp.size){\n      i = hp.pop();\n      res += hp.val[i];\n      rep(j,es[i]) hp.change(edge[i][j], cost[i][j]);\n    }\n    return res;\n  }\n  int MST_Prim_cost(T &res, void *mem = wmem){\n    int i, j, cnt = 0;\n    res = 0;\n    DijkstraHeap<T> hp;\n    hp.walloc(N, &mem);\n    hp.init(N);\n    hp.change(0,0);\n    while(hp.size){\n      i = hp.pop();\n      res += hp.val[i];\n      cnt++;\n      rep(j,es[i]) hp.change(edge[i][j], cost[i][j]);\n    }\n    if(cnt==N) return 1;\n    return 0;\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"DijkstraHeap");
      need[n] = d;
      parent[n] = "wgraph";
    }
    {
      string n = "wgraph_Rerooting_head";
      string c = "template<class V> void Rerooting(V res[], void *mem = wmem);\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "wgraph";
    }
    {
      string n = "wgraph_Rerooting";
      string c = "template<class T> template<class V>\nvoid wgraph<T>::Rerooting(V res[], void *mem /*= wmem*/){\n  int i, j, k, ui, ii;\n  int *upind, *dwind, *vis, *q, qs, qe;\n  V **val, **pleft, **pright, tmp;\n  \n  walloc1d(&upind, N, &mem);\n  walloc1d(&dwind, N, &mem);\n  walloc1d(&vis, N, &mem);\n  walloc1d(&q, N, &mem);\n  walloc1d(&val, N, &mem);\n  rep(i,N) walloc1d(&val[i], es[i], &mem);\n  walloc1d(&pleft, N, &mem);\n  rep(i,N) walloc1d(&pleft[i], es[i]+1, &mem);\n  walloc1d(&pright, N, &mem);\n  rep(i,N) walloc1d(&pright[i], es[i]+1, &mem);\n  upind[0] = -1;\n  rep(i,N) vis[i] = 0;\n  qs = qe = 0;\n  q[qe++] = 0; vis[0] = 1;\n  while(qs < qe){\n    i = q[qs++];\n    rep(j,es[i]){\n      k = edge[i][j];\n      if(vis[k]){\n        upind[i] = j;\n        continue;\n      }\n      dwind[k] = j;\n      vis[k] = 1;\n      q[qe++] = k;\n    }\n  }\n  rrep(ii,N){\n    i = q[ii];\n    rep(j,es[i]) if(j != upind[i]){\n      k = edge[i][j];\n      ui = upind[k];\n      tmp = wRerootingMerge(pleft[k][ui], pright[k][es[k]-1-ui]);\n      tmp = wRerootingNode(tmp, k);\n      val[i][j] = wRerootingEdge(tmp, cost[i][j], i, k);\n    }\n    wRerootingId(pleft[i][0]);\n    wRerootingId(pright[i][0]);\n    rep(j,es[i]){\n      if(j == upind[i]) break;\n      pleft[i][j+1] = wRerootingMerge(pleft[i][j], val[i][j]);\n    }\n    rep(j,es[i]){\n      if(es[i]-1-j == upind[i]) break;\n      pright[i][j+1] = wRerootingMerge(pright[i][j], val[i][es[i]-1-j]);\n    }\n  }\n  rep(ii,N){\n    i = q[ii];\n    j = upind[i];\n    if(j != -1){\n      k = edge[i][j];\n      ui = dwind[i];\n      tmp = wRerootingMerge(pleft[k][ui], pright[k][es[k]-1-ui]);\n      tmp = wRerootingNode(tmp, k);\n      val[i][j] = wRerootingEdge(tmp, cost[i][j], i, k);\n      rep(j,upind[i],es[i]){\n        pleft[i][j+1] = wRerootingMerge(pleft[i][j], val[i][j]);\n      }\n      rep(j,es[i]-upind[i]-1,es[i]){\n        pright[i][j+1] = wRerootingMerge(pright[i][j], val[i][es[i]-1-j]);\n      }\n    }\n    res[i] = wRerootingNode(pleft[i][es[i]], i);\n  }\n}\n";
      string p = "last";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"wgraph_Rerooting_head");
      need[n] = d;
      parent[n] = "wgraph";
    }
    {
      string n = "wgraph_end";
      string c = "};\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "wgraph";
    }


    {
      string n = "HLD";
      string c = "struct HLD{\n  int N;\n  int *es, **edge;\n\n  int *group, *groupind;\n  int groupNum, *groupSize, **groupNode, *groupUpNode, *groupDepth;\n\n  void init(graph g, void **mem = &wmem){\n    init(g.N, g.es, g.edge, mem);\n  }\n  \n  void init(int N__, int *es__, int **edge__, void **mem = &wmem){\n    int i, j, k, x, y, mx;\n    int *q, q_st, q_ed, *sz;\n    char *vis;\n    void *tmpmem;\n\n    N = N__;\n    es = es__;\n    edge = edge__;\n\n    walloc1d(&group, N, mem);\n    walloc1d(&groupind, N, mem);\n\n    tmpmem = *mem;\n    walloc1d(&q, N, &tmpmem);\n    walloc1d(&sz, N, &tmpmem);\n    walloc1d(&vis, N, &tmpmem);\n    rep(i,N) vis[i] = 0;\n    q_st = 0; q_ed = 1;\n    q[0] = 0; vis[0] = 1;\n    while(q_st < q_ed){\n      i = q[q_st++];\n      rep(j,es[i]){\n        k = edge[i][j];\n        if(!vis[k]) vis[k] = 1, q[q_ed++] = k;\n      }\n    }\n    \n    rep(i,N) sz[i] = 0;\n    for(j=N-1;j>=0;j--){\n      i = q[j];\n      sz[i] = 1;\n      rep(k,es[i]) sz[i] += sz[edge[i][k]];\n    }\n\n    rep(i,N) group[i] = -1;\n\n    groupNum = 0;\n    rep(j,N){\n      i = q[j];\n      if(group[i]>=0) continue;\n\n      group[i] = groupNum++;\n      groupind[i] = 0;\n      for(;;){\n        mx = -1;\n        rep(k,es[i]){\n          if(group[edge[i][k]] != -1) continue;\n          if(mx==-1) mx = k;\n          else if(sz[edge[i][k]] > sz[edge[i][mx]]) mx = k;\n        }\n        if(mx==-1) break;\n        group[edge[i][mx]] = group[i];\n        groupind[edge[i][mx]] = groupind[i]+1;\n        i = edge[i][mx];\n      }\n    }\n\n    walloc1d(&groupSize, groupNum, mem);\n    walloc1d(&groupUpNode, groupNum, mem);\n    walloc1d(&groupDepth, groupNum, mem);\n\n    rep(i,groupNum) groupSize[i] = 0;\n    rep(i,N) groupSize[group[i]]++;\n    walloc1d(&groupNode, groupNum, mem);\n    rep(i,groupNum) walloc1d(&groupNode[i], groupSize[i], mem);\n    rep(i,N) groupNode[group[i]][groupind[i]] = i;\n\n    rep(i,groupNum) groupDepth[i] = -1;\n    groupUpNode[0] = -1;\n    groupDepth[0] = 0;\n    rep(x,groupNum) rep(y,groupSize[x]){\n      i = groupNode[x][y];\n      rep(j,es[i]){\n        k = edge[i][j];\n        if(x != group[k] && groupDepth[group[k]]==-1){\n          groupUpNode[group[k]] = i;\n          groupDepth[group[k]] = groupDepth[x] + 1;\n        }\n      }\n    }\n  }\n  \n  int lca(int x, int y){\n    int x1, y1, x2, y2;\n    x1 = group[x]; x2 = groupind[x];\n    y1 = group[y]; y2 = groupind[y];\n    while(groupDepth[x1] > groupDepth[y1]){\n      x = groupUpNode[x1];\n      x1 = group[x]; x2 = groupind[x];\n    }\n    while(groupDepth[x1] < groupDepth[y1]){\n      y = groupUpNode[y1];\n      y1 = group[y]; y2 = groupind[y];\n    }\n    while(x1 != y1){\n      x = groupUpNode[x1];\n      x1 = group[x]; x2 = groupind[x];\n      y = groupUpNode[y1];\n      y1 = group[y]; y2 = groupind[y];\n    }\n    \n    if(x2 <= y2) return x;\n    return y;\n  }\n\n  int depth(int x){\n    int x1, x2, res = 0;\n    x1 = group[x];\n    x2 = groupind[x];\n    while(groupUpNode[x1] != -1){\n      res += x2 + 1;\n      x = groupUpNode[x1];\n      x1 = group[x];\n      x2 = groupind[x];\n    }\n    return res + x2;\n  }\n\n  int dist(int x, int y){\n    int x1, y1, x2, y2, res = 0;\n    x1 = group[x]; x2 = groupind[x];\n    y1 = group[y]; y2 = groupind[y];\n    while(groupDepth[x1] > groupDepth[y1]){\n      res += x2 + 1;\n      x = groupUpNode[x1];\n      x1 = group[x]; x2 = groupind[x];\n    }\n    while(groupDepth[x1] < groupDepth[y1]){\n      res += y2 + 1;\n      y = groupUpNode[y1];\n      y1 = group[y]; y2 = groupind[y];\n    }\n    while(x1 != y1){\n      res += x2 + y2 + 2;\n      x = groupUpNode[x1];\n      x1 = group[x]; x2 = groupind[x];\n      y = groupUpNode[y1];\n      y1 = group[y]; y2 = groupind[y];\n    }\n\n    if(x2 <= y2) return res + y2 - x2;\n    return res + x2 - y2;\n  }\n\n  int up(int x){\n    int x1 = group[x];\n    int x2 = groupind[x];\n    if(x2==0) return groupUpNode[x1];\n    return groupNode[x1][x2-1];\n  }\n\n  int up(int x, int d){\n    int x1 = group[x];\n    int x2 = groupind[x];\n    while(d > x2){\n      if(groupUpNode[x1]==-1) return -1;\n      d -= x2 + 1;\n      x = groupUpNode[x1];\n      x1 = group[x];\n      x2 = groupind[x];\n    }\n    return groupNode[x1][x2-d];\n  }\n};\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"graph");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "HLD_fenwick";
      string c = "template<class T>\nstruct HLD_fenwick{\n  HLD *hld;\n  fenwick<T> *fen;\n\n  void init(HLD *hld__, void **mem = &wmem){\n    int i, j;\n\n    hld = hld__;\n    walloc1d(&fen, hld->groupNum, mem);\n\n    rep(i,hld->groupNum){\n      fen[i].walloc(hld->groupSize[i], mem);\n      fen[i].init(hld->groupSize[i]);\n    }\n  }\n\n  inline void add(int u, T val){\n    int ug, ui;\n    ug = hld->group[u];\n    ui = hld->groupind[u];\n    fen[ug].add(ui, val);\n  }\n\n  inline T get(int u, int v){\n    T res;\n    int ug, vg, ui, vi;\n    ug = hld->group[u];\n    vg = hld->group[v];\n\n    res = 0;\n    while(ug != vg){\n      if(hld->groupDepth[ug] < hld->groupDepth[vg]){\n        swap(u, v);\n        swap(ug, vg);\n      }\n      res += fen[ug].get(hld->groupind[u]);\n      u = hld->groupUpNode[ug];\n      ug = hld->group[u];\n    }\n    ui = hld->groupind[u];\n    vi = hld->groupind[v];\n    res += fen[ug].range(min(ui,vi), max(ui,vi));\n    return res;\n  }\n};\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"graph");
      d.push_back((string)"HLD");
      d.push_back((string)"min_L");
      d.push_back((string)"max_L");
      d.push_back((string)"fenwick");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "HLD_segtree";
      string c = "template<class T>\nstruct HLD_segtree{\n  HLD *hld;\n  segtree<T> *seg;\n\n  void init(HLD *hld__, T initval[], void **mem = &wmem){\n    int i, j;\n\n    hld = hld__;\n    walloc1d(&seg, hld->groupNum, mem);\n\n    rep(i,hld->groupNum){\n      seg[i].walloc(hld->groupSize[i], 0, mem);\n      seg[i].setN(hld->groupSize[i], 0, 0);\n      if(initval!=NULL) rep(j,hld->groupSize[i]) seg[i][j] = initval[ hld->groupNode[i][j] ];\n      else              rep(j,hld->groupSize[i]) seg[i][j] = 0;\n      seg[i].build();\n    }\n  }\n\n  inline void change(int u, int v, T val){\n    int ug, vg, ui, vi;\n    ug = hld->group[u];\n    vg = hld->group[v];\n    while(ug != vg){\n      if(hld->groupDepth[ug] < hld->groupDepth[vg]){\n        swap(u, v);\n        swap(ug, vg);\n      }\n      seg[ug].change(0, hld->groupind[u]+1, val);\n      u = hld->groupUpNode[ug];\n      ug = hld->group[u];\n    }\n    ui = hld->groupind[u];\n    vi = hld->groupind[v];\n    seg[ug].change(min(ui,vi), max(ui,vi)+1, val);\n  }\n\n  inline void add(int u, int v, T val){\n    int ug, vg, ui, vi;\n    ug = hld->group[u];\n    vg = hld->group[v];\n    while(ug != vg){\n      if(hld->groupDepth[ug] < hld->groupDepth[vg]){\n        swap(u, v);\n        swap(ug, vg);\n      }\n      seg[ug].add(0, hld->groupind[u]+1, val);\n      u = hld->groupUpNode[ug];\n      ug = hld->group[u];\n    }\n    ui = hld->groupind[u];\n    vi = hld->groupind[v];\n    seg[ug].add(min(ui,vi), max(ui,vi)+1, val);\n  }\n\n  inline pair<T,int> getMin(int u, int v){\n    pair<T,int> res, tmp;\n    int ug, vg, ui, vi;\n    ug = hld->group[u];\n    vg = hld->group[v];\n\n    res.first = numeric_limits<T>::max();\n    res.second = -1;\n    while(ug != vg){\n      if(hld->groupDepth[ug] < hld->groupDepth[vg]){\n        swap(u, v);\n        swap(ug, vg);\n      }\n      tmp = seg[ug].getMin(0, hld->groupind[u]+1);\n      tmp.second = hld->groupNode[ug][tmp.second];\n      res <?= tmp;\n      u = hld->groupUpNode[ug];\n      ug = hld->group[u];\n    }\n    ui = hld->groupind[u];\n    vi = hld->groupind[v];\n    tmp = seg[ug].getMin(min(ui,vi), max(ui,vi)+1);\n    tmp.second = hld->groupNode[ug][tmp.second];\n    res <?= tmp;\n    return res;\n  }\n\n  inline T getMinVal(int u, int v){\n    return getMin(u,v).first;\n  }\n\n  inline int getMinInd(int u, int v){\n    return getMin(u,v).second;\n  }\n\n  inline T getSum(int u, int v){\n    T res;\n    int ug, vg, ui, vi;\n    ug = hld->group[u];\n    vg = hld->group[v];\n\n    res = 0;\n    while(ug != vg){\n      if(hld->groupDepth[ug] < hld->groupDepth[vg]){\n        swap(u, v);\n        swap(ug, vg);\n      }\n      res += seg[ug].getSum(0, hld->groupind[u]+1);\n      u = hld->groupUpNode[ug];\n      ug = hld->group[u];\n    }\n    ui = hld->groupind[u];\n    vi = hld->groupind[v];\n    res += seg[ug].getSum(min(ui,vi), max(ui,vi)+1);\n    return res;\n  }\n\n  inline void change_edge(int u, int v, T val){\n    int x, z, d;\n    z = hld->lca(u, v);\n    d = hld->depth(z);\n    if(z != u){\n      x = hld->up(u, hld->depth(u) - d - 1);\n      change(u, x, val);\n    }\n    if(z != v){\n      x = hld->up(v, hld->depth(v) - d - 1);\n      change(v, x, val);\n    }\n  }\n\n  inline void add_edge(int u, int v, T val){\n    int x, z, d;\n    z = hld->lca(u, v);\n    d = hld->depth(z);\n    if(z != u){\n      x = hld->up(u, hld->depth(u) - d - 1);\n      add(u, x, val);\n    }\n    if(z != v){\n      x = hld->up(v, hld->depth(v) - d - 1);\n      add(v, x, val);\n    }\n  }\n\n  inline pair<T,int> getMin_edge(int u, int v){\n    int x, z, d;\n    pair<T,int> res, tmp;\n    res.first = numeric_limits<T>::max();\n    res.second = -1;\n    \n    z = hld->lca(u, v);\n    d = hld->depth(z);\n    if(z != u){\n      x = hld->up(u, hld->depth(u) - d - 1);\n      tmp = getMin(u, x);\n      res <?= tmp;\n    }\n    if(z != v){\n      x = hld->up(v, hld->depth(v) - d - 1);\n      tmp = getMin(v, x);\n      res <?= tmp;\n    }\n    return res;\n  }\n\n  inline T getMinVal_edge(int u, int v){\n    return getMin_edge(u,v).first;\n  }\n\n  inline int getMinInd_edge(int u, int v){\n    return getMin_edge(u,v).second;\n  }\n\n  inline T getSum_edge(int u, int v){\n    int x, z, d;\n    T res;\n    res = 0;\n    \n    z = hld->lca(u, v);\n    d = hld->depth(z);\n    if(z != u){\n      x = hld->up(u, hld->depth(u) - d - 1);\n      res += getSum(u, x);\n    }\n    if(z != v){\n      x = hld->up(v, hld->depth(v) - d - 1);\n      res += getSum(v, x);\n    }\n    return res;\n  }\n};\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"workmemory");
      d.push_back((string)"walloc1d");
      d.push_back((string)"graph");
      d.push_back((string)"HLD");
      d.push_back((string)"chmin");
      d.push_back((string)"min_L");
      d.push_back((string)"max_L");
      d.push_back((string)"segtree");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "maxflow";
      string c = "template<class T, class S>\nstruct maxflow{\n  int node, st, ed;\n  int *es, *emem, **edge, **rev, *level, *qq;\n  T **flow, eps;\n  void malloc(int N){\n    int i;\n    es = (int*)std::malloc(N*sizeof(int));\n    emem = (int*)std::malloc(N*sizeof(int));\n    level = (int*)std::malloc(N*sizeof(int));\n    qq = (int*)std::malloc(N*sizeof(int));\n    edge = (int**)std::malloc(N*sizeof(int*));\n    rev = (int**)std::malloc(N*sizeof(int*));\n    flow = (T**)std::malloc(N*sizeof(T*));\n    rep(i,N) emem[i] = 0, edge[i] = rev[i] = NULL, flow[i] = NULL;\n  }\n  void malloc(int N, int init_flag){\n    int i;\n    es = (int*)std::malloc(N*sizeof(int));\n    emem = (int*)std::malloc(N*sizeof(int));\n    level = (int*)std::malloc(N*sizeof(int));\n    qq = (int*)std::malloc(N*sizeof(int));\n    edge = (int**)std::malloc(N*sizeof(int*));\n    rev = (int**)std::malloc(N*sizeof(int*));\n    flow = (T**)std::malloc(N*sizeof(T*));\n    rep(i,N) emem[i] = 0, edge[i] = rev[i] = NULL, flow[i] = NULL;\n    if(init_flag) init(N);\n  }\n  void walloc(int N, void**mem = &wmem){\n    int i;\n    walloc1d(&es, N, mem);\n    walloc1d(&emem, N, mem);\n    walloc1d(&level, N, mem);\n    walloc1d(&qq, N, mem);\n    walloc1d(&edge, N, mem);\n    walloc1d(&rev, N, mem);\n    walloc1d(&flow, N, mem);\n    (*mem) = (flow + N);\n  }\n  void walloc(int N, int init_flag, void**mem = &wmem){\n    int i;\n    walloc1d(&es, N, mem);\n    walloc1d(&emem, N, mem);\n    walloc1d(&level, N, mem);\n    walloc1d(&qq, N, mem);\n    walloc1d(&edge, N, mem);\n    walloc1d(&rev, N, mem);\n    walloc1d(&flow, N, mem);\n    (*mem) = (flow + N);\n    if(init_flag) init(N);\n  }\n  void levelize(void){\n    int i, j, k, t;\n    int q_st = 0, q_ed = 1;\n    rep(i,node) level[i] = -1;\n    level[st] = 0;\n    qq[0] = st;\n    while(q_st != q_ed){\n      i = qq[q_st++];\n      t = level[i] + 1;\n      rep(j,es[i]) if(flow[i][j] > eps){\n        k = edge[i][j];\n        if(level[k]!=-1) continue;\n        level[k] = t;\n        qq[q_ed++] = k;\n        if(k==ed) return;\n      }\n    }\n  }\n  S pushflow(int i, S lim){\n    int j, k, ji;\n    S s, t, res = 0;\n    if(i==ed) return lim;\n    rep(j,es[i]) if(flow[i][j] > eps){\n      k = edge[i][j];\n      if(level[k] != level[i]+1) continue;\n      s = min(lim, (S)flow[i][j]);\n      t = pushflow(k, s); if(!t) continue;\n      res += t;\n      lim -= t;\n      ji = rev[i][j];\n      flow[i][j] -= t; flow[k][ji] += t;\n      if(!lim) break;\n    }\n    if(lim) level[i] = -1;\n    return res;\n  }\n  S solve(int st_, int ed_){\n    S res = 0;\n    st = st_; ed = ed_;\n    for(;;){\n      levelize();\n      if(level[ed] == -1) break;\n      res += pushflow(st, numeric_limits<S>::max());\n    }\n    return res;\n  }\n  void init(int N){\n    int i;\n    node = N;\n    rep(i,N) es[i] = 0;\n    eps = (T)1e-9;\n  }\n  void memoryExpand(int i, int sz){\n    if(sz <= emem[i]) return;\n    sz = max(sz, max(3, emem[i]*2));\n    emem[i]=sz;\n    edge[i] = (int*)realloc(edge[i], sz*sizeof(int));\n    rev[i] = (int*)realloc(rev[i], sz*sizeof(int));\n    flow[i] = (T*)realloc(flow[i], sz*sizeof(T));\n  }\n  void addEdge(int n1, int n2, T f1, T f2 = 0){\n    int s1 = es[n1]++, s2 = es[n2]++;\n    if(s1 >= emem[n1]) memoryExpand(n1, es[n1]);\n    if(s2 >= emem[n2]) memoryExpand(n2, es[n2]);\n    edge[n1][s1]=n2; edge[n2][s2]=n1;\n    flow[n1][s1]=f1; flow[n2][s2]=f2;\n    rev[n1][s1]=s2; rev[n2][s2]=s1;\n  }\n  \n  void addEdgeAdv(int n1, int n2, T f1, T f2 = 0){\n    int s1 = es[n1]++, s2 = es[n2]++;\n    edge[n1][s1]=n2; edge[n2][s2]=n1;\n    flow[n1][s1]=f1; flow[n2][s2]=f2;\n    rev[n1][s1]=s2; rev[n2][s2]=s1;\n  }\n  \n  void setGraph(int N, int M, int n1[], int n2[], T f1[], T f2[]){\n    int i;\n    node = N;\n    rep(i,N) es[i] = 0;\n    rep(i,M) es[n1[i]]++, es[n2[i]]++;\n    rep(i,N) memoryExpand(i, es[i]);\n    rep(i,N) es[i] = 0;\n    rep(i,M) addEdgeAdv(n1[i], n2[i], f1[i], f2[i]);\n    eps = (T)1e-9;\n  }\n  void setGraph_w(int N, int M, int n1[], int n2[], T f1[], T f2[], void **mem = wmem){\n    int i, j, k;\n    node = N;\n    rep(i,N) es[i] = emem[i] = 0;\n    rep(i,M) es[n1[i]]++, es[n2[i]]++;\n    \n    edge[0] = (int*)(*mem);\n    REP(i,1,N) edge[i] = edge[i-1] + es[i-1];\n    rev[0] = edge[N-1] + es[N-1];\n    REP(i,1,N) rev[i] = rev[i-1] + es[i-1];\n    flow[0] = (T*)(rev[N-1] + es[N-1]);\n    REP(i,1,N) flow[i] = flow[i-1] + es[i-1];\n    *mem = (void*)(flow[N-1] + es[N-1]);\n    \n    rep(i,N) es[i] = 0;\n    rep(i,M) addEdgeAdv(n1[i], n2[i], f1[i], f2[i]);\n    eps = (T)1e-9;\n  }\n};\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"min_L");
      d.push_back((string)"max_L");
      d.push_back((string)"workmemory");
      d.push_back((string)"walloc1d");
      need[n] = d;
    }

    {
      string n = "minCostFlow";
      string c = "template<class FT, class CT>\nstruct minCostFlow {\n  int node;\n  int *es, *emem, **edge, **rev;\n  FT **flow, f_eps;\n  CT **cost, *potential, c_eps;\n  \n  LHeap<CT> hp;\n  char *reached;\n  FT *cur_flow;\n  CT *cur_cost;\n  int *back_edge;\n  void malloc(int N){\n    int i;\n    es = (int*)std::malloc(N*sizeof(int));\n    emem = (int*)std::malloc(N*sizeof(int));\n    edge = (int**)std::malloc(N*sizeof(int*));\n    rev = (int**)std::malloc(N*sizeof(int*));\n    flow = (FT**)std::malloc(N*sizeof(FT*));\n    cost = (CT**)std::malloc(N*sizeof(CT*));\n    rep(i,N){\n      emem[i] = 0;\n      edge[i] = rev[i] = NULL;\n      flow[i] = NULL;\n      cost[i] = NULL;\n    }\n    \n    hp.malloc(N);\n    reached = (char*)std::malloc(N*sizeof(char));\n    cur_flow = (FT*)std::malloc(N*sizeof(FT));\n    cur_cost = (CT*)std::malloc(N*sizeof(CT));\n    potential = (CT*)std::malloc(N*sizeof(CT));\n    back_edge = (int*)std::malloc(N*sizeof(int));\n    node = N;\n    rep(i,N) es[i] = 0;\n    f_eps = (FT)1e-9;\n    c_eps = (CT)1e-9;\n  }\n  void init(int N){\n    int i;\n    node = N;\n    rep(i,N) es[i] = 0;\n    f_eps = (FT)1e-9;\n    c_eps = (CT)1e-9;\n  }\n  void memoryExpand(int i, int sz){\n    if(sz <= emem[i]) return;\n    sz = max(sz, 3, 2emem[i]);\n    emem[i] = sz;\n    edge[i] = (int*)realloc(edge[i], sz*sizeof(int));\n    rev[i] = (int*)realloc(rev[i], sz*sizeof(int));\n    flow[i] = (FT*)realloc(flow[i], sz*sizeof(FT));\n    cost[i] = (CT*)realloc(cost[i], sz*sizeof(CT));\n  }\n  void addEdge(int n1, int n2, FT f, CT c){\n    int s1 = es[n1]++;\n    int s2 = es[n2]++;\n    if(s1 >= emem[n1]) memoryExpand(n1, es[n1]);\n    if(s2 >= emem[n2]) memoryExpand(n2, es[n2]);\n    edge[n1][s1] = n2; edge[n2][s2] = n1;\n    rev[n1][s1]  = s2; rev[n2][s2]  = s1;\n    flow[n1][s1] = f; flow[n2][s2] = 0;\n    cost[n1][s1] = c; cost[n2][s2] = -c;\n  }\n  template<class FTS, class CTS>\n  void solve(int st, int ed, FTS &fres, CTS &cres, FT flim = -1, CT clim = 0){\n    int i, j, k, l;\n    FT f;\n    CT nc;\n    fres = 0;\n    cres = 0;\n    rep(i,node) potential[i] = 0;\n    for(;;){\n      if(flim >= -f_eps && flim <= f_eps) break;\n      hp.init(node);\n      rep(i,node) reached[i] = 0;\n      reached[st] = 1;\n      cur_cost[st] = 0;\n      l = 0;\n      hp.change(st, cur_cost[st]);\n      while(hp.size){\n        i = hp.pop();\n        rep(j, es[i]){\n          if(flow[i][j] <= f_eps) continue;\n          k = edge[i][j];\n          nc = cur_cost[i] + cost[i][j] + potential[i] - potential[k];\n          if(reached[k]==0 || cur_cost[k] > nc+c_eps){\n            reached[k] = 1;\n            cur_cost[k] = nc;\n            cur_flow[k] = flow[i][j];\n            if(i != st) cur_flow[k] <?= cur_flow[i];\n            back_edge[k] = rev[i][j];\n            hp.change(k, cur_cost[k]);\n          }\n        }\n      }\n      if(reached[ed]==0) break;\n      if(flim==-2 && cur_cost[ed] + potential[ed] >= clim) break;\n      f = cur_flow[ed];\n      if(flim >= -f_eps){\n        f <?= flim;\n        flim -= f;\n      }\n      if(f <= f_eps) break;\n      rep(i,node) if(reached[i]) potential[i] += cur_cost[i];\n      fres += f;\n      cres += f * potential[ed];\n      i = ed;\n      while(i != st){\n        j = back_edge[i];\n        k = edge[i][j];\n        flow[i][j] += f;\n        flow[k][rev[i][j]] -= f;\n        i = k;\n      }\n    }\n  }\n};\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"chmin");
      d.push_back((string)"min_L");
      d.push_back((string)"max_L");
      d.push_back((string)"LHeap");
      need[n] = d;
    }


    {
      string n = "KMP";
      string c = "template<class T>\nint KMP(T A[], int As, T B[], int Bs, int res[] = NULL, int *fail = (int*)wmem){\n  int i, k, cnt = 0;\n  k = fail[0] = -1;\n  rep(i,Bs){\n    while(k>=0 && B[k]!=B[i]) k = fail[k];\n    fail[i+1] = ++k;\n  }\n  if(res != NULL) rep(i,As) res[i] = 0;\n  k = 0;\n  rep(i,As){\n    while(k >= 0 && B[k] != A[i]) k = fail[k];\n    k++;\n    if(k == Bs){\n      cnt++;\n      if(res != NULL) res[i-Bs+1] = 1;\n      k = fail[k];\n    }\n  }\n  return cnt;\n}\n\ntemplate<class T, class S>\nint KMP(T A[], int As, T B[], int Bs, S res[], int *fail = (int*)wmem){\n  int i, k, cnt = 0;\n  k = fail[0] = -1;\n  rep(i,Bs){\n    while(k>=0 && B[k]!=B[i]) k = fail[k];\n    fail[i+1] = ++k;\n  }\n  if(res != NULL) rep(i,As) res[i] = 0;\n  k = 0;\n  rep(i,As){\n    while(k >= 0 && B[k] != A[i]) k = fail[k];\n    k++;\n    if(k == Bs){\n      cnt++;\n      if(res != NULL) res[i-Bs+1] = 1;\n      k = fail[k];\n    }\n  }\n  return cnt;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "isSubsequence";
      string c = "template<class T1, class T2>\nint isSubsequence(int As, const T1 A[], int Bs, const T2 B[]){\n  int i, j = 0;\n  if(Bs==0) return 1;\n  rep(i,As) if(A[i]==B[j]){\n    j++;\n    if(j==Bs) break;\n  }\n  return j == Bs;\n}\nint isSubsequence(string A, string B){\n  int i, j = 0;\n  if(B.size()==0) return 1;\n  rep(i,A.size()) if(A[i]==B[j]){\n    j++;\n    if(j==B.size()) break;\n  }\n  return j == B.size();\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "isSubsequence_r";
      string c = "int isSubsequence_r(string &A, string &B){\n  int i, j = 0;\n  if(B.size()==0) return 1;\n  rep(i,A.size()) if(A[i]==B[j]){\n    j++;\n    if(j==B.size()) break;\n  }\n  return j == B.size();\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "cntSubsequence";
      string c = "template<class R, class T>\nR cntSubsequence(int As, const T A[], int Bs, const T B[], void *mem = wmem){\n  int i, j, k;\n  int *aa, *bb, sz;\n  int *s, **arr;\n  R *dp;\n  if(Bs > As) return 0;\n  walloc1d(&aa, As, &mem);\n  walloc1d(&bb, Bs, &mem);\n  sz = coordcomp(As, A, Bs, B, aa, bb, mem);\n  walloc1d(&s, sz, &mem);\n  walloc1d(&arr, sz, &mem);\n  rep(i,sz) s[i] = 0;\n  rep(i,Bs) s[bb[i]]++;\n  rep(i,sz) if(s[i]) walloc1d(&arr[i], s[i], &mem);\n  rep(i,sz) s[i] = 0;\n  rrep(i,Bs) arr[bb[i]][s[bb[i]]++] = i;\n  walloc1d(&dp, Bs+1, &mem);\n  dp[0] = 1;\n  rep(i,1,Bs+1) dp[i] = 0;\n  rep(i,As){\n    k = aa[i];\n    rep(j,s[k]) dp[arr[k][j]+1] += dp[arr[k][j]];\n  }\n  return dp[Bs];\n}\ntemplate<class R, class T1, class T2>\nR cntSubsequence(int As, const T1 A[], int Bs, const T2 B[], void *mem = wmem){\n  int i;\n  typename cLtraits_common_type<T1,T2>::type *aa, *bb;\n  if(Bs > As) return 0;\n  walloc1d(&aa, As, &mem);\n  walloc1d(&bb, Bs, &mem);\n  rep(i,As) aa[i] = A[i];\n  rep(i,Bs) bb[i] = B[i];\n  return cntSubsequence<R>(As, aa, Bs, bb, mem);\n}\ntemplate<class R>\nR cntSubsequence(string A, string B, void *mem = wmem){\n  return cntSubsequence<R>((int)A.size(), A.c_str(), (int)B.size(), B.c_str(), mem);\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"workmemory");
      d.push_back((string)"walloc1d");
      d.push_back((string)"coordcomp_2");
      d.push_back((string)"cLtraits_common_type");
      need[n] = d;
    }


    {
      string n = "isSubstring_string";
      string c = "int isSubstring(string A, string B, void *mem = wmem){\n  int i = 0, k, *fail;\n  char *m;\n\n  if(B.size() > A.size()) return 0;\n  \n  walloc1d(&fail, B.size()+1, &mem);\n\n  k = fail[0] = -1;\n  rep(i,B.size()){\n    while(k>=0 && B[k] != B[i]) k = fail[k];\n    fail[i+1] = ++k;\n  }\n\n  k = 0;\n  rep(i,A.size()){\n    while(k >= 0 && B[k] != A[i]) k = fail[k];\n    if((++k) == B.size()) return 1;\n  }\n\n  return 0;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"workmemory");
      d.push_back((string)"walloc1d");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "isPalindrome";
      string c = "template<class T>\ninline int isPalindrome(const int N, const T A[]){\n  int i = 0, j = N-1;\n  while(i < j){\n    if(A[i] != A[j]) return 0;;\n    i++; j--;\n  }\n  return 1;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "longestSuffixPrefix";
      string c = "template<class T>\nint longestSuffixPrefix(int As, T A[], int Bs, T B[], void *mem = wmem){\n  int i, k, res;\n  int *fail;\n\n  if(As > Bs) A += As-Bs, As = Bs;\n  if(As < Bs) Bs = As;\n\n  walloc1d(&fail, Bs, &mem);\n  \n  k = fail[0] = -1;\n  rep(i,Bs){\n    while(k>=0 && B[k]!=B[i]) k = fail[k];\n    fail[i+1] = ++k;\n  }\n\n  res = 0;\n  rep(i,As){\n    while(res && A[i]!=B[res]) res = fail[res];\n    if(A[i]==B[res]) res++;\n  }\n  return res;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"min_L");
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "strReplace1";
      string c = "string strReplace_L(string str, string bef, string aft, void *mem = wmem){\n  int i = 0, k, *fail;\n  char *m;\n  string res;\n  \n  walloc1d(&fail, bef.size()+1, &mem);\n  walloc1d(&m, str.size(), &mem);\n  rep(i,str.size()) m[i] = 0;\n\n  k = fail[0] = -1;\n  rep(i,bef.size()){\n    while(k>=0 && bef[k] != bef[i]) k = fail[k];\n    fail[i+1] = ++k;\n  }\n\n  k = 0;\n  rep(i,str.size()){\n    while(k >= 0 && bef[k] != str[i]) k = fail[k];\n    k++;\n    if(k == bef.size()){\n      m[i-bef.size()+1] = 1;\n      k = fail[k];\n    }\n  }\n\n  i = 0;\n  while(i < str.size()){\n    if(m[i]){\n      res += aft;\n      i += bef.size();\n    } else {\n      res += str[i++];\n    }\n  }\n\n  return res;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "strReplace2";
      string c = "string strReplace_L(string str, vector<string> bef, vector<string> aft, void *mem = wmem){\n  int i = 0, k, q, *fail, *m;\n  string res;\n  \n  walloc1d(&fail, bef.size()+1, &mem);\n  walloc1d(&m, str.size(), &mem);\n  rep(i,str.size()) m[i] = -1;\n\n  rrep(q,bef.size()){\n    k = fail[0] = -1;\n    rep(i,bef[q].size()){\n      while(k>=0 && bef[q][k] != bef[q][i]) k = fail[k];\n      fail[i+1] = ++k;\n    }\n\n    k = 0;\n    rep(i,str.size()){\n      while(k >= 0 && bef[q][k] != str[i]) k = fail[k];\n      k++;\n      if(k == bef[q].size()){\n        m[i-bef[q].size()+1] = q;\n        k = fail[k];\n      }\n    }\n  }\n\n  i = 0;\n  while(i < str.size()){\n    if(m[i] >= 0){\n      res += aft[m[i]];\n      i += bef[m[i]].size();\n    } else {\n      res += str[i++];\n    }\n  }\n\n  return res;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "smallestSubsequenceLengthK";
      string c = "template<class T>\nvoid smallestSubsequenceLengthK(int N, T A[], int K, T res[], void *mem = wmem){\n  int i, d = N - K, s = 0;\n  T *st;\n  walloc1d(&st, N, &mem);\n  rep(i,N){\n    while(s > 0 && d > 0 && A[i] < st[s-1]) s--, d--;\n    st[s++] = A[i];\n  }\n  rep(i,K) res[i] = st[i];\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "Trie";
      string c = "struct Trie{\n  int node, alphabet;\n  int node_mem, alphabet_mem;\n  int **nx;\n\n  void init(const int k){\n    int i;\n    node = 1;\n    alphabet = k;\n    rep(i,alphabet) nx[0][i] = -1;\n  }\n\n  void init(void){\n    init(alphabet_mem);\n  }\n\n  void malloc(const int n, const int k){\n    malloc2d(&nx,n,k);\n    node_mem = n;\n    alphabet_mem = k;\n    init();\n  }\n\n  void walloc(const int n, const int k, void **mem = &wmem){\n    walloc2d(&nx,n,k,mem);\n    node_mem = n;\n    alphabet_mem = k;\n    init();\n  }\n\n  void free(void){\n    free2d(nx);\n  }\n\n  template<class T> int addWord(const T word[], const int len){\n    int i, j, k, now = 0;\n    rep(i,len){\n      if(nx[now][word[i]]==-1){\n        k = node++;\n        nx[now][word[i]] = k;\n        rep(j,alphabet) nx[k][j] = -1;\n      }\n      now = nx[now][word[i]];\n    }\n    return now;\n  }\n\n  template<class T> inline int addNext(const int n, const T c){\n    int j, k;\n    if(nx[n][c] != -1) return nx[n][c];\n    k = node++;\n    nx[n][c] = k;\n    rep(j,alphabet) nx[k][j] = -1;\n    return k;\n  }\n\n  template<class T> inline int next(const int n, const T c){\n    return nx[n][c];\n  }\n};\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"malloc1d");
      d.push_back((string)"malloc2d");
      d.push_back((string)"walloc1d");
      d.push_back((string)"walloc2d");
      d.push_back((string)"free1d");
      d.push_back((string)"free2d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "AhoCorasick";
      string c = "struct AhoCorasick{\n  int node, mem, alphabet;\n  int **nx, *failed;\n  int **ind, *indsz, *indmem;\n\n  void init(void){\n    int i;\n    node = 1;\n    rep(i,alphabet) nx[0][i] = -1;\n    failed[0] = 0;\n    indsz[0] = 0;\n  }\n\n  void malloc(const int n, const int k){\n    int i;\n    malloc2d(&nx,n,k);\n    malloc1d(&failed,n);\n    malloc1d(&ind,n);\n    malloc1d(&indsz,n);\n    malloc1d(&indmem,n);\n    node = n;\n    alphabet = k;\n    rep(i,n) indmem[i] = 0;\n    init();\n  }\n\n  void free(void){\n    free2d(nx);\n    free1d(failed);\n    free1d(ind);\n    free1d(indsz);\n    free1d(indmem);\n  }\n\n  inline void addEnd(const int n, const int id){\n    int s;\n    if(indsz[n]+1 > indmem[n]){\n      s = indmem[n] * 2 + 1;\n      if(indmem[n]==0) ind[n] = (int*) std::malloc(s * sizeof(int));\n      else             ind[n] = (int*) std::realloc(ind[n], s * sizeof(int));\n      indmem[n] = s;\n    }\n    ind[n][indsz[n]++] = id;\n  }\n\n  template<class T> int addWord(const T word[], const int len, int id){\n    int i, j, k, now = 0;\n    rep(i,len){\n      if(nx[now][word[i]]==-1){\n        k = node++;\n        nx[now][word[i]] = k;\n        rep(j,alphabet) nx[k][j] = -1;\n        indsz[k] = 0;\n      }\n      now = nx[now][word[i]];\n    }\n    addEnd(now, id);\n    return now;\n  }\n\n  void construct(void *mem = wmem){\n    int i, j, k, now;\n    int *q, qs, qe;\n\n    q = (int*) mem;\n    qs = qe = 0;\n\n    now = 0;\n    rep(k,alphabet) if(nx[now][k] != -1){\n      q[qe++] = nx[now][k];\n      failed[ nx[now][k] ] = now;\n    }\n\n    while(qs < qe){\n      now = q[qs++];\n      rep(k,alphabet) if(nx[now][k] != -1){\n        i = failed[now];\n        while(i){\n          if(nx[i][k] != -1) break;\n          i = failed[i];\n        }\n        if(nx[i][k] != -1) i = nx[i][k];\n        failed[ nx[now][k] ] = i;\n        rep(j,indsz[i]) addEnd(nx[now][k], ind[i][j]);\n        q[qe++] = nx[now][k];\n      }\n    }\n  }\n\n  template<class T> inline int next(const int n, const T c){\n    int i, now;\n    now = n;\n    if(nx[n][c]!=-1) return nx[n][c];\n    while(now && nx[now][c]==-1) now=failed[now];\n    if(nx[now][c]!=-1) now = nx[now][c];\n    return nx[n][c] = now;\n  }\n};\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"malloc1d");
      d.push_back((string)"malloc2d");
      d.push_back((string)"free1d");
      d.push_back((string)"free2d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "AhoCorasick_Sum";
      string c = "template<class S>\nstruct AhoCorasick_Sum{\n  int node, mem, alphabet;\n  int **nx, *failed;\n  S *sum;\n\n  void init(void){\n    int i;\n    node = 1;\n    rep(i,alphabet) nx[0][i] = -1;\n    failed[0] = 0;\n    sum[0] = 0;\n  }\n\n  void malloc(const int n, const int k){\n    int i;\n    malloc2d(&nx,n,k);\n    malloc1d(&failed,n);\n    malloc1d(&sum,n);\n    node = n;\n    alphabet = k;\n    init();\n  }\n\n  void free(void){\n    free2d(nx);\n    free1d(failed);\n    free1d(sum);\n  }\n\n  template<class T> int addWord(const T word[], const int len, S val){\n    int i, j, k, now = 0;\n    rep(i,len){\n      if(nx[now][word[i]]==-1){\n        k = node++;\n        nx[now][word[i]] = k;\n        rep(j,alphabet) nx[k][j] = -1;\n        sum[k] = 0;\n      }\n      now = nx[now][word[i]];\n    }\n    sum[now] += val;\n    return now;\n  }\n\n  void construct(void *mem = wmem){\n    int i, j, k, now;\n    int *q, qs, qe;\n\n    q = (int*) mem;\n    qs = qe = 0;\n\n    now = 0;\n    rep(k,alphabet) if(nx[now][k] != -1){\n      q[qe++] = nx[now][k];\n      failed[ nx[now][k] ] = now;\n    }\n\n    while(qs < qe){\n      now = q[qs++];\n      rep(k,alphabet) if(nx[now][k] != -1){\n        i = failed[now];\n        while(i){\n          if(nx[i][k] != -1) break;\n          i = failed[i];\n        }\n        if(nx[i][k] != -1) i = nx[i][k];\n        failed[ nx[now][k] ] = i;\n        sum[ nx[now][k] ] += sum[i];\n        q[qe++] = nx[now][k];\n      }\n    }\n  }\n\n  template<class T> inline int next(const int n, const T c){\n    int i, now;\n    now = n;\n    if(nx[n][c]!=-1) return nx[n][c];\n    while(now && nx[now][c]==-1) now=failed[now];\n    if(nx[now][c]!=-1) now = nx[now][c];\n    return nx[n][c] = now;\n  }\n};\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"malloc1d");
      d.push_back((string)"malloc2d");
      d.push_back((string)"free1d");
      d.push_back((string)"free2d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "twoMultisets";
      string c = "template<class T>\nstruct twoMultisets{\n  multiset<T> a, b;\n  T sa, sb;\n\n  twoMultisets(){\n    clear();\n  }\n\n  void clear(){\n    a.clear();\n    b.clear();\n    sa = sb = 0;\n  }\n\n  void insert(T x){\n    if(b.size() == 0 || x < *b.begin()){\n      a.insert(x);\n      sa += x;\n    } else {\n      b.insert(x);\n      sb += x;\n    }\n  }\n\n  int erase(T x){\n    typename multiset<T>::iterator it;\n    it = a.find(x);\n    if(it != a.end()){\n      a.erase(it);\n      sa -= x;\n      return 1;\n    }\n    it = b.find(x);\n    if(it != b.end()){\n      b.erase(it);\n      sb -= x;\n      return 1;\n    }\n    return 0;\n  }\n\n  int size(void){\n    return a.size() + b.size();\n  }\n\n  T allsum(void){\n    return sa + sb;\n  }\n\n  void assign(int K){\n    T x;\n    typename multiset<T>::iterator it;\n\n    while(a.size() < K){\n      x = *b.begin();\n      b.erase(b.begin());\n      a.insert(x);\n      sa += x;\n      sb -= x;\n    }\n\n    while(a.size() > K){\n      it = a.end();\n      it--;\n      x = *it;\n      a.erase(it);\n      b.insert(x);\n      sa -= x;\n      sb += x;\n    }\n  }\n\n  T Kth(int K){\n    assign(K);\n    return *b.begin();\n  }\n\n  T Ksum(int K){\n    assign(K);\n    return sa;\n  }\n\n  T rKth(int K){\n    return Kth(a.size() + b.size() - K - 1);\n  }\n\n  T rKsum(int K){\n    assign(a.size() + b.size() - K);\n    return sb;\n  }\n\n  T getMin(void){\n    if(a.size()) return *a.begin();\n    if(b.size()) return *b.begin();\n    return 0;\n  }\n\n  T getMin(T x){\n    if(a.size()) return *a.begin();\n    if(b.size()) return *b.begin();\n    return x;\n  }\n\n  T getMax(void){\n    if(b.size()) return *b.rbegin();\n    if(a.size()) return *a.rbegin();\n    return 0;\n  }\n\n  T getMax(T x){\n    if(b.size()) return *b.rbegin();\n    if(a.size()) return *a.rbegin();\n    return x;\n  }\n};\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "inversion_range";
      string c = "ll inversion_range(const int N, const int A[], const int mn, const int mx, void *mem=wmem){\n  int i, j, k;\n  ll res = 0;\n  fenwick<int> t;\n  t.walloc(mx-mn+1, &mem);\n  t.init(mx-mn+1);\n  for(i=N-1;i>=0;i--){\n    res += t.get(A[i]-mn-1);\n    t.add(A[i]-mn,1);\n  }\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"workmemory");
      d.push_back((string)"fenwick");
      need[n] = d;
    }
    {
      string n = "inversion";
      string c = "template<class T>\nll inversion(const int N, const T A[], void *mem=wmem){\n  int i, j, k, p;\n  int n1, n2;\n  T *x, *y;\n  ll res = 0;\n  walloc1d(&x, N, &mem);\n  walloc1d(&y, N, &mem);\n  rep(i,N) x[i] = A[i];\n  for(k=0;k<N;k+=4){\n    n1 = min(4, N-k);\n    for(j=n1;j;j--){\n      REP(i,1,j) if(x[k+i-1] > x[k+i]){\n        swap(x[k+i-1], x[k+i]);\n        res++;\n      }\n    }\n  }\n  p = 4;\n  while(p<N){\n    for(k=0;k<N;k+=2*p){\n      n1 = min(p,N-k);\n      n2 = min(p,N-k-n1);\n      \n      i = j = 0;\n      while(i<n1 && j<n2){\n        if(x[k+i] <= x[k+n1+j]){\n          y[k+i+j] = x[k+i];\n          i++;\n        } else {\n          y[k+i+j] = x[k+n1+j];\n          res += n1-i;\n          j++;\n        }\n      }\n      while(i<n1){\n        y[k+i+j] = x[k+i];\n        i++;\n      }\n      while(j<n2){\n        y[k+i+j] = x[k+n1+j];\n        j++;\n      }\n    }\n    \n    swap(x,y);\n    p *= 2;\n  }\n  return res;\n}\ntemplate<class T>\nll inversion(const vector<T> &A, void *mem = wmem){\n  const int N = A.size();\n  int i, j, k, p;\n  int n1, n2;\n  T *x, *y;\n  ll res = 0;\n  walloc1d(&x, N, &mem);\n  walloc1d(&y, N, &mem);\n  rep(i,N) x[i] = A[i];\n  for(k=0;k<N;k+=4){\n    n1 = min(4, N-k);\n    for(j=n1;j;j--){\n      REP(i,1,j) if(x[k+i-1] > x[k+i]){\n        swap(x[k+i-1], x[k+i]);\n        res++;\n      }\n    }\n  }\n  p = 4;\n  while(p<N){\n    for(k=0;k<N;k+=2*p){\n      n1 = min(p,N-k);\n      n2 = min(p,N-k-n1);\n      \n      i = j = 0;\n      while(i<n1 && j<n2){\n        if(x[k+i] <= x[k+n1+j]){\n          y[k+i+j] = x[k+i];\n          i++;\n        } else {\n          y[k+i+j] = x[k+n1+j];\n          res += n1-i;\n          j++;\n        }\n      }\n      while(i<n1){\n        y[k+i+j] = x[k+i];\n        i++;\n      }\n      while(j<n2){\n        y[k+i+j] = x[k+n1+j];\n        j++;\n      }\n    }\n    \n    swap(x,y);\n    p *= 2;\n  }\n  return res;\n}\nll inversion(const string &A, void *mem = wmem){\n  return inversion(A.size(), A.c_str(), mem);\n}\ntemplate<class T>\nll inversion(const int N, const T A[], const T B[], void *mem = wmem){\n  int i, k, sz, *aa, *bb, *hist1, *hist2, **ind;\n  walloc1d(&aa, N, &mem);\n  walloc1d(&bb, N, &mem);\n  sz = coordcomp(N,A,N,B,aa,bb,mem);\n  if(sz > N) return -1;\n  walloc1d(&hist1, sz, &mem);\n  walloc1d(&hist2, sz, &mem);\n  rep(i,sz) hist1[i] = 0;\n  rep(i,sz) hist2[i] = 0;\n  rep(i,N) hist1[aa[i]]++;\n  rep(i,N) hist2[bb[i]]++;\n  rep(i,sz) if(hist1[i] != hist2[i]) return -1;\n  walloc1d(&ind, sz, &mem);\n  rep(i,sz) walloc1d(&ind[i], hist2[i], &mem);\n  rep(i,sz) hist2[i] = 0;\n  rep(i,N){\n    k = bb[i];\n    ind[k][hist2[k]++] = i;\n  }\n  rrep(i,N){\n    k = aa[i];\n    aa[i] = ind[k][--hist1[k]];\n  }\n  return inversion_range(N,aa,0,N-1,mem);\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"min_L");
      d.push_back((string)"inversion_range");
      d.push_back((string)"coordcomp_2");
      need[n] = d;
    }


    {
      string n = "ZetaTransform";
      string c = "template<class T, class S>\nvoid ZetaTransform_L(int N, T A[], S res[]){\n  int i, j, x, m;\n  rep(j,N) res[j] = A[j];\n  for(i=1;i<N;i*=2){\n    rep(j,i,N,2*i){\n      m = min(N, j+i);\n      rep(x,j,m) res[x] += res[x^i];\n    }\n  }\n}\ntemplate<class T>\nvoid ZetaTransform(int N, T A[]){\n  int i, j, x, m;\n  for(i=1;i<N;i*=2){\n    rep(j,i,N,2*i){\n      m = min(N, j+i);\n      rep(x,j,m) A[x] += A[x^i];\n    }\n  }\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"min_L");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "ZetaTransform_min";
      string c = "template<class T, class S>\nvoid ZetaTransform_min_L(int N, T A[], S res[]){\n  int i, j, x, m;\n  rep(j,N) res[j] = A[j];\n  for(i=1;i<N;i*=2){\n    rep(j,i,N,2*i){\n      m = min(N, j+i);\n      rep(x,j,m) res[x] <?= res[x^i];\n    }\n  }\n}\ntemplate<class T>\nvoid ZetaTransform_min_L(int N, T A[]){\n  int i, j, x, m;\n  for(i=1;i<N;i*=2){\n    rep(j,i,N,2*i){\n      m = min(N, j+i);\n      rep(x,j,m) A[x] <?= A[x^i];\n    }\n  }\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"min_L");
      d.push_back((string)"chmin");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "ZetaTransform_min2";
      string c = "template<class T, class S>\nvoid ZetaTransform_min_L(int N, T A[], S r1[], S r2[]){\n  int i, j, k, x, m;\n  rep(j,N) r1[j] = A[j];\n  rep(j,N) r2[j] = numeric_limits<S>::max();\n  for(i=1;i<N;i*=2){\n    rep(j,i,N,2*i){\n      m = min(N, j+i);\n      rep(x,j,m){\n        k = (x^i);\n        r2[x] <?= r1[k];\n        sortE(r1[x], r2[x]);\n        r2[x] <?= r2[k];\n      }\n    }\n  }\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"min_L");
      d.push_back((string)"chmin");
      d.push_back((string)"sortE");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "ZetaTransform_min3";
      string c = "template<class T, class S>\nvoid ZetaTransform_min_L(int N, T A[], S r1[], S r2[], S r3[]){\n  int i, j, k, x, m;\n  rep(j,N) r1[j] = A[j];\n  rep(j,N) r2[j] = numeric_limits<S>::max();\n  rep(j,N) r3[j] = numeric_limits<S>::max();\n  for(i=1;i<N;i*=2){\n    rep(j,i,N,2*i){\n      m = min(N, j+i);\n      rep(x,j,m){\n        k = (x^i);\n        r3[x] <?= r1[k];\n        sortE(r2[x], r3[x]);\n        sortE(r1[x], r2[x]);\n        r3[x] <?= r2[k];\n        sortE(r2[x], r3[x]);\n        r3[x] <?= r3[k];\n      }\n    }\n  }\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"min_L");
      d.push_back((string)"chmin");
      d.push_back((string)"sortE");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "ZetaTransform_max";
      string c = "template<class T, class S>\nvoid ZetaTransform_max_L(int N, T A[], S res[]){\n  int i, j, x, m;\n  rep(j,N) res[j] = A[j];\n  for(i=1;i<N;i*=2){\n    rep(j,i,N,2*i){\n      m = min(N, j+i);\n      rep(x,j,m) res[x] >?= res[x^i];\n    }\n  }\n}\ntemplate<class T>\nvoid ZetaTransform_max_L(int N, T A[]){\n  int i, j, x, m;\n  for(i=1;i<N;i*=2){\n    rep(j,i,N,2*i){\n      m = min(N, j+i);\n      rep(x,j,m) A[x] >?= A[x^i];\n    }\n  }\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"min_L");
      d.push_back((string)"chmax");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "MoebiusTransform";
      string c = "template<class T, class S>\nvoid MoebiusTransform(int N, T A[], S res[]){\n  int i, j, x, m;\n  rep(j,N) res[j] = A[j];\n  for(i=1;i<N;i*=2){\n    rep(j,i,N,2*i){\n      m = min(N, j+i);\n      rep(x,j,m) res[x] -= res[x^i];\n    }\n  }\n}\ntemplate<class T>\nvoid MoebiusTransform(int N, T A[]){\n  int i, j, x, m;\n  for(i=1;i<N;i*=2){\n    rep(j,i,N,2*i){\n      m = min(N, j+i);\n      rep(x,j,m) A[x] -= A[x^i];\n    }\n  }\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"min_L");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "ZetaTransform2";
      string c = "template<class T, class S>\nvoid ZetaTransform2(int N, T A[], S res[]){\n  int i, j, x, m;\n  rep(j,N) res[j] = A[j];\n  for(i=1;i<N;i*=2){\n    rep(j,i,N,2*i){\n      m = min(N, j+i);\n      rep(x,j,m) res[x^i] += res[x];\n    }\n  }\n}\ntemplate<class T>\nvoid ZetaTransform2(int N, T A[]){\n  int i, j, x, m;\n  for(i=1;i<N;i*=2){\n    rep(j,i,N,2*i){\n      m = min(N, j+i);\n      rep(x,j,m) A[x^i] += A[x];\n    }\n  }\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"min_L");
      need[n] = d;
    }
    {
      string n = "MoebiusTransform2";
      string c = "template<class T, class S>\nvoid MoebiusTransform2(int N, T A[], S res[]){\n  int i, j, x, m;\n  rep(j,N) res[j] = A[j];\n  for(i=1;i<N;i*=2){\n    rep(j,i,N,2*i){\n      m = min(N, j+i);\n      rep(x,j,m) res[x^i] -= res[x];\n    }\n  }\n}\ntemplate<class T>\nvoid MoebiusTransform2(int N, T A[]){\n  int i, j, x, m;\n  for(i=1;i<N;i*=2){\n    rep(j,i,N,2*i){\n      m = min(N, j+i);\n      rep(x,j,m) A[x^i] -= A[x];\n    }\n  }\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"min_L");
      need[n] = d;
    }

    {
      string n = "HadamardTransform";
      string c = "template<class T>\nvoid HadamardTransform(int N, T A[]){\n  int i, j, k;\n  T x, y;\n  for(i=1;i<N;i*=2){\n    rep(j,0,N,2*i){\n      rep(k,i){\n        x = A[j+k];\n        y = A[j+k+i];\n        A[j+k] = x + y;\n        A[j+k+i] = x - y;\n      }\n    }\n  }\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "xorConvolution";
      string c = "template<class T1, class T2, class T3>\nvoid xorConvolution(int As, T1 A[], int Bs, T2 B[], int Rs, T3 R[], void *mem = wmem){\n  int i, n = 1, m;\n  T1 *aa;\n  T2 *bb;\n  T3 iv;\n  while(n < As) n *= 2;\n  while(n < Bs) n *= 2;\n  walloc1d(&aa, n, &mem);\n  walloc1d(&bb, n, &mem);\n  rep(i,As) aa[i] = A[i];\n  rep(i,As,n) aa[i] = 0;\n  rep(i,Bs) bb[i] = B[i];\n  rep(i,Bs,n) bb[i] = 0;\n  HadamardTransform(n,aa);\n  HadamardTransform(n,bb);\n  rep(i,n) aa[i] *= bb[i];\n  HadamardTransform(n,aa);\n  iv = 1;\n  iv /= n;\n  m = min(Rs, n);\n  rep(i,m) R[i] = aa[i] * iv;\n  rep(i,m,Rs) R[i] = 0;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"HadamardTransform");
      d.push_back((string)"min_L");
      need[n] = d;
    }
    {
      string n = "orConvolution";
      string c = "template<class T1, class T2, class T3>\nvoid orConvolution(int As, T1 A[], int Bs, T2 B[], int Rs, T3 R[], void *mem = wmem){\n  int i, n = 1, m;\n  T1 *aa;\n  T2 *bb;\n  while(n < As) n *= 2;\n  while(n < Bs) n *= 2;\n  walloc1d(&aa, n, &mem);\n  walloc1d(&bb, n, &mem);\n  rep(i,As) aa[i] = A[i];\n  rep(i,As,n) aa[i] = 0;\n  rep(i,Bs) bb[i] = B[i];\n  rep(i,Bs,n) bb[i] = 0;\n  ZetaTransform(n,aa);\n  ZetaTransform(n,bb);\n  rep(i,n) aa[i] *= bb[i];\n  MoebiusTransform(n,aa);\n  m = min(Rs, n);\n  rep(i,m) R[i] = aa[i];\n  rep(i,m,Rs) R[i] = 0;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"ZetaTransform");
      d.push_back((string)"MoebiusTransform");
      d.push_back((string)"min_L");
      need[n] = d;
    }
    {
      string n = "andConvolution";
      string c = "template<class T1, class T2, class T3>\nvoid andConvolution(int As, T1 A[], int Bs, T2 B[], int Rs, T3 R[], void *mem = wmem){\n  int i, n = 1, m;\n  T1 *aa;\n  T2 *bb;\n  while(n < As) n *= 2;\n  while(n < Bs) n *= 2;\n  walloc1d(&aa, n, &mem);\n  walloc1d(&bb, n, &mem);\n  rep(i,As) aa[i] = A[i];\n  rep(i,As,n) aa[i] = 0;\n  rep(i,Bs) bb[i] = B[i];\n  rep(i,Bs,n) bb[i] = 0;\n  ZetaTransform2(n,aa);\n  ZetaTransform2(n,bb);\n  rep(i,n) aa[i] *= bb[i];\n  MoebiusTransform2(n,aa);\n  m = min(Rs, n);\n  rep(i,m) R[i] = aa[i];\n  rep(i,m,Rs) R[i] = 0;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"ZetaTransform2");
      d.push_back((string)"MoebiusTransform2");
      d.push_back((string)"min_L");
      need[n] = d;
    }


    {
      string n = "slideMin";
      string c = "template<class T>\nvoid slideMin(int n, int k, T in[], T res[], void *mem = wmem){\n  int i, s = 0;\n  T *q;\n  int q_st = 0, q_size = 0;\n  walloc1d(&q, n);\n  rep(i,n){\n    while(q_size && q[q_st+q_size-1] > in[i]) q_size--;\n    q[q_st+q_size++] = in[i];\n    if(i>=k && in[i-k]==q[q_st]) q_st++, q_size--;\n    if(i>=k-1) res[s++] = q[q_st];\n  }\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"workmemory");
      d.push_back((string)"walloc1d");
      need[n] = d;
    }
    {
      string n = "slideMax";
      string c = "template<class T>\nvoid slideMax(int n, int k, T in[], T res[], void *mem = wmem){\n  int i, s = 0;\n  T *q;\n  int q_st = 0, q_size = 0;\n  walloc1d(&q, n);\n  rep(i,n){\n    while(q_size && q[q_st+q_size-1] < in[i]) q_size--;\n    q[q_st+q_size++] = in[i];\n    if(i>=k && in[i-k]==q[q_st]) q_st++, q_size--;\n    if(i>=k-1) res[s++] = q[q_st];\n  }\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"workmemory");
      d.push_back((string)"walloc1d");
      need[n] = d;
    }

    {
      string n = "isLeapYear";
      string c = "template<class T>\ninline int isLeapYear(const T y){\n  if(y%4) return 0;\n  if(y%100) return 1;\n  if(y%400) return 0;\n  return 1;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "numOfDaysInMonth1";
      string c = "inline int numOfDaysInMonth(const int m){\n  if(m==2) return 28;\n  if(m==4||m==6||m==9||m==11) return 30;\n  return 31;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"isLeapYear");
      need[n] = d;
    }
    {
      string n = "numOfDaysInMonth2";
      string c = "template<class T>\ninline int numOfDaysInMonth(const T y, const int m){\n  return numOfDaysInMonth(m) if[m==2, + isLeapYear(y)];\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"isLeapYear");
      d.push_back((string)"numOfDaysInMonth1");
      need[n] = d;
    }
    {
      string n = "dayOfWeek";
      string c = "template<class T>\ninline int dayOfWeek(T Y, int M, int D){\n  int i, j;\n  if(M <= 2) M += 12, Y--;\n  i = (Y / 100) % 4;\n  j = Y % 100;\n  return (D + 13*(M+1)/5 + j + j/4 + 5*i + i/4 + 5) % 7;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "dayOfWeekStr";
      string c = "const char* WeekStr[7] = {\"Monday\", \"Tuesday\", \"Wednesday\", \"Thursday\", \"Friday\", \"Saturday\", \"Sunday\"};\ninline const char* dayOfWeekStr(int w){\n  return WeekStr[w];\n}\ntemplate<class T>\ninline const char* dayOfWeekStr(T Y, int M, int D){\n  return WeekStr[dayOfWeek(Y, M, D)];\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"dayOfWeek");
      need[n] = d;
    }
    {
      string n = "prevDay";
      string c = "template<class T>\ninline void prevDay(T &y, int &m, int &d){\n  d--;\n  if(d==0){\n    m--;\n    if(m==0) y--, m = 12;\n    d = numOfDaysInMonth(y, m);\n  }\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"numOfDaysInMonth2");
      need[n] = d;
    }
    {
      string n = "nextDay";
      string c = "template<class T>\ninline void nextDay(T &y, int &m, int &d){\n  d++;\n  if(d > numOfDaysInMonth(y,m)){\n    d = 1;\n    m++;\n    if(m==13) m = 1, y++;\n  }\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"numOfDaysInMonth2");
      need[n] = d;
    }
    {
      string n = "dayIndex";
      string c = "inline ll dayIndex(ll y, int m, int d){\n  static int sm[13] = {0, 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};\n  ll p, res = 1;\n  p = (y-1) / 400;\n  res += p * 146097;\n  y -= p * 400;\n  res += (y-1) * 365 + (y-1) / 4 - (y-1) / 100 + sm[m] + d-1;\n  if(m >= 3 && isLeapYear(y)) res++;\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"isLeapYear");
      need[n] = d;
    }
    {
      string n = "dayFromIndex";
      string c = "template<class S, class T>\nvoid dayFromIndex(S ind, T &y, int &m, int &d){\n  int i, k;\n  static int mc[13] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};\n  ind--;\n  y = 1;\n  y += (ind / 146097) * 400;\n  ind %= 146097;\n  y += (ind / 36524) * 100;\n  ind %= 36524;\n  y += (ind / 1461) * 4;\n  ind %= 1461;\n  rep(i,3) if(ind >= 365) ind -= 365, y++;\n  rep(i,1,13){\n    k = mc[i];\n    if(i==2 && isLeapYear(y)) k++;\n    if(ind < k){\n      m = i;\n      d = ind + 1;\n      break;\n    } else {\n      ind -= k;\n    }\n  }\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"isLeapYear");
      need[n] = d;
    }


    {
      string n = "isVowel";
      string c = "inline int isVowel(const char c){\n  if(c=='a'||c=='i'||c=='u'||c=='e'||c=='o') return 1;\n  if(c=='A'||c=='I'||c=='U'||c=='E'||c=='O') return 1;\n  return 0;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "Arr1d";
      string c = "template<class T>\nstruct Arr1d{\n  int n, mem;\n  T *d;\n  T& operator[](int a){\n    return d[a];\n  }\n  void sort(){\n    reset();\n    std::sort(d, d+n);\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"Arr1d_end");
      d.push_back((string)"Arr1d_reset_head");
      d.push_back((string)"Arr1d_reset_foot");
      d.push_back((string)"Arr1d_constructor_head");
      d.push_back((string)"Arr1d_constructor_foot");
      d.push_back((string)"Arr1d_destructor_head");
      d.push_back((string)"Arr1d_destructor_foot");
      need[n] = d;
    }
    {
      string n = "Arr1d_getSum";
      string c = "  int set_cumulative_sum;\n  int cumulative_sum_mem;\n  T *cumulative_sum;\n  void setSum(void){\n    int i;\n    set_cumulative_sum = 1;\n    if(cumulative_sum_mem < n+1){\n      delete[] cumulative_sum;\n      cumulative_sum = new T[n+1];\n      cumulative_sum_mem = n+1;\n    }\n    cumulative_sum[0] = 0;\n    rep(i,n) cumulative_sum[i+1] = cumulative_sum[i] + d[i];\n  }\n  template<class T1, class T2>\n  T getSum(T1 i, T2 j){\n    if(i > j) return 0;\n    if(set_cumulative_sum==0) setSum();\n    return cumulative_sum[min(j+1, n)] - cumulative_sum[max(0, i)];\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Arr1d_getSum_reset");
      d.push_back((string)"Arr1d_getSum_constructor");
      d.push_back((string)"Arr1d_getSum_destructor");
      d.push_back((string)"min_L");
      d.push_back((string)"max_L");
      need[n] = d;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_ConstLenLeft";
      string c = "  int set_const_len_left;\n  int const_len_left_mem;\n  int *const_len_left;\n  void setConstLenLeft(void){\n    int i;\n    set_const_len_left = 1;\n    if(const_len_left_mem < n){\n      delete[] const_len_left;\n      const_len_left = new int[n];\n      const_len_left_mem = n;\n    }\n    rep(i,n) const_len_left[i] = 1;\n    rep(i,1,n) if(d[i]==d[i-1]) const_len_left[i] = const_len_left[i-1] + 1;\n  }\n  int ConstLenLeft(int st, T val){\n    if(!set_const_len_left) setConstLenLeft();\n    if(val != d[st]) return 0;\n    return const_len_left[st];\n  }\n  int ConstLenLeft(int st){\n    if(!set_const_len_left) setConstLenLeft();\n    return const_len_left[st];\n  }\n  int ConstLenLeftCyclic(int st, T val){\n    if(!set_const_len_left) setConstLenLeft();\n    st %= n;\n    if(st < 0) st += n;\n    if(val != d[st]) return 0;\n    if(const_len_left[st] != st+1 || d[st] != d[n-1]) return const_len_left[st];\n    if(const_len_left[n-1] == n) return int_inf;\n    return const_len_left[st] + const_len_left[n-1];\n  }\n  int ConstLenLeftCyclic(int st){\n    if(!set_const_len_left) setConstLenLeft();\n    st %= n;\n    if(st < 0) st += n;\n    if(const_len_left[st] != st+1 || d[st] != d[n-1]) return const_len_left[st];\n    if(const_len_left[n-1] == n) return int_inf;\n    return const_len_left[st] + const_len_left[n-1];\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Arr1d_ConstLenLeft_reset");
      d.push_back((string)"Arr1d_ConstLenLeft_constructor");
      d.push_back((string)"Arr1d_ConstLenLeft_destructor");
      need[n] = d;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_ConstLenRight";
      string c = "  int set_const_len_right;\n  int const_len_right_mem;\n  int *const_len_right;\n  void setConstLenRight(void){\n    int i;\n    set_const_len_right = 1;\n    if(const_len_right_mem < n){\n      delete[] const_len_right;\n      const_len_right = new int[n];\n      const_len_right_mem = n;\n    }\n    rep(i,n) const_len_right[i] = 1;\n    rrep(i,n-1) if(d[i]==d[i+1]) const_len_right[i] = const_len_right[i+1] + 1;\n  }\n  int ConstLenRight(int st, T val){\n    if(!set_const_len_right) setConstLenRight();\n    if(val != d[st]) return 0;\n    return const_len_right[st];\n  }\n  int ConstLenRight(int st){\n    if(!set_const_len_right) setConstLenRight();\n    return const_len_right[st];\n  }\n  int ConstLenRightCyclic(int st, T val){\n    if(!set_const_len_right) setConstLenRight();\n    if(val != d[st]) return 0;\n    st %= n;\n    if(st < 0) st += n;\n    if(const_len_right[st] != n-st || d[st] != d[0]) return const_len_right[st];\n    if(const_len_right[0] == n) return int_inf;\n    return const_len_right[st] + const_len_right[0];\n  }\n  int ConstLenRightCyclic(int st){\n    if(!set_const_len_right) setConstLenRight();\n    st %= n;\n    if(st < 0) st += n;\n    if(const_len_right[st] != n-st || d[st] != d[0]) return const_len_right[st];\n    if(const_len_right[0] == n) return int_inf;\n    return const_len_right[st] + const_len_right[0];\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Arr1d_ConstLenRight_reset");
      d.push_back((string)"Arr1d_ConstLenRight_constructor");
      d.push_back((string)"Arr1d_ConstLenRight_destructor");
      need[n] = d;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_dHist";
      string c = "  int set_dhist;\n  int dhist_mem;\n  int *dhist, *dhists;\n  T dhist_mn, dhist_mx;\n  void setDHist(void){\n    int i, len;\n    set_dhist = 1;\n    if(n==0) return;\n    dhist_mn = dhist_mx = d[0];\n    rep(i,1,n){\n      if(dhist_mn > d[i]) dhist_mn = d[i];\n      if(dhist_mx < d[i]) dhist_mx = d[i];\n    }\n    len = dhist_mx - dhist_mn + 1;\n    if(dhist_mem < len){\n      delete[] dhist;\n      dhist = new int[len];\n      delete[] dhists;\n      dhists = new int[len+1];\n      dhist_mem = len;\n    }\n    rep(i,len) dhist[i] = 0;\n    rep(i,n) dhist[d[i] - dhist_mn]++;\n    dhists[0] = 0;\n    rep(i,len) dhists[i+1] = dhists[i] + dhist[i];\n  }\n  int dHist(T x){\n    if(set_dhist==0) setDHist();\n    if(n == 0 || x < dhist_mn || x > dhist_mx) return 0;\n    return dhist[x - dhist_mn];\n  }\n  int dHist(T x, T y){\n    if(set_dhist==0) setDHist();\n    if(x < dhist_mn) x = dhist_mn;\n    if(y > dhist_mx) y = dhist_mx;\n    if(n == 0 || x > y) return 0;\n    return dhists[y-dhist_mn+1] - dhists[x-dhist_mn];\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Arr1d_dHist_reset");
      d.push_back((string)"Arr1d_dHist_constructor");
      d.push_back((string)"Arr1d_dHist_destructor");
      need[n] = d;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_sHist";
      string c = "  int set_shist;\n  int shist_mem;\n  T *shist;\n  void setSHist(void){\n    int i;\n    set_shist = 1;\n    if(shist_mem < n){\n      delete[] shist;\n      shist = new T[n];\n      shist_mem = n;\n    }\n    rep(i,n) shist[i] = d[i];\n    std::sort(shist, shist + n);\n  }\n  int sHist(T x){\n    if(set_shist==0) setSHist();\n    auto r = equal_range(shist, shist+n, x);\n    return r.second - r.first;\n  }\n  int sHist(T x, T y){\n    if(set_shist==0) setSHist();\n    return upper_bound(shist, shist+n, y) - lower_bound(shist, shist+n, x);\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Arr1d_sHist_reset");
      d.push_back((string)"Arr1d_sHist_constructor");
      d.push_back((string)"Arr1d_sHist_destructor");
      need[n] = d;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_PrevLE";
      string c = "  int set_prevLE;\n  int prevLE_mem;\n  int *prevLE;\n  void setPrevLE(void *mem = wmem){\n    int i;\n    int s = 0, *st;\n    set_prevLE = 1;\n    if(prevLE_mem < n){\n      delete[] prevLE;\n      prevLE = new int[n];\n      prevLE_mem = n;\n    }\n    walloc1d(&st, n, &mem);\n    rep(i,n){\n      while(s && d[st[s-1]] > d[i]) s--;\n      if(s==0){\n        prevLE[i] = -1;\n      } else {\n        prevLE[i] = st[s-1];\n      }\n      st[s++] = i;\n    }\n  }\n  int PrevLE(int i){\n    if(set_prevLE==0) setPrevLE();\n    return prevLE[i];\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Arr1d_PrevLE_reset");
      d.push_back((string)"Arr1d_PrevLE_constructor");
      d.push_back((string)"Arr1d_PrevLE_destructor");
      need[n] = d;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_PrevLT";
      string c = "  int set_prevLT;\n  int prevLT_mem;\n  int *prevLT;\n  void setPrevLT(void *mem = wmem){\n    int i;\n    int s = 0, *st;\n    set_prevLT = 1;\n    if(prevLT_mem < n){\n      delete[] prevLT;\n      prevLT = new int[n];\n      prevLT_mem = n;\n    }\n    walloc1d(&st, n, &mem);\n    rep(i,n){\n      while(s && d[st[s-1]] >= d[i]) s--;\n      if(s==0){\n        prevLT[i] = -1;\n      } else {\n        prevLT[i] = st[s-1];\n      }\n      st[s++] = i;\n    }\n  }\n  int PrevLT(int i){\n    if(set_prevLT==0) setPrevLT();\n    return prevLT[i];\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Arr1d_PrevLT_reset");
      d.push_back((string)"Arr1d_PrevLT_constructor");
      d.push_back((string)"Arr1d_PrevLT_destructor");
      need[n] = d;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_PrevGE";
      string c = "  int set_prevGE;\n  int prevGE_mem;\n  int *prevGE;\n  void setPrevGE(void *mem = wmem){\n    int i;\n    int s = 0, *st;\n    set_prevGE = 1;\n    if(prevGE_mem < n){\n      delete[] prevGE;\n      prevGE = new int[n];\n      prevGE_mem = n;\n    }\n    walloc1d(&st, n, &mem);\n    rep(i,n){\n      while(s && d[st[s-1]] < d[i]) s--;\n      if(s==0){\n        prevGE[i] = -1;\n      } else {\n        prevGE[i] = st[s-1];\n      }\n      st[s++] = i;\n    }\n  }\n  int PrevGE(int i){\n    if(set_prevGE==0) setPrevGE();\n    return prevGE[i];\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Arr1d_PrevGE_reset");
      d.push_back((string)"Arr1d_PrevGE_constructor");
      d.push_back((string)"Arr1d_PrevGE_destructor");
      need[n] = d;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_PrevGT";
      string c = "  int set_prevGT;\n  int prevGT_mem;\n  int *prevGT;\n  void setPrevGT(void *mem = wmem){\n    int i;\n    int s = 0, *st;\n    set_prevGT = 1;\n    if(prevGT_mem < n){\n      delete[] prevGT;\n      prevGT = new int[n];\n      prevGT_mem = n;\n    }\n    walloc1d(&st, n, &mem);\n    rep(i,n){\n      while(s && d[st[s-1]] <= d[i]) s--;\n      if(s==0){\n        prevGT[i] = -1;\n      } else {\n        prevGT[i] = st[s-1];\n      }\n      st[s++] = i;\n    }\n  }\n  int PrevGT(int i){\n    if(set_prevGT==0) setPrevGT();\n    return prevGT[i];\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Arr1d_PrevGT_reset");
      d.push_back((string)"Arr1d_PrevGT_constructor");
      d.push_back((string)"Arr1d_PrevGT_destructor");
      need[n] = d;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_NextLE";
      string c = "  int set_nextLE;\n  int nextLE_mem;\n  int *nextLE;\n  void setNextLE(void *mem = wmem){\n    int i;\n    int s = 0, *st;\n    set_nextLE = 1;\n    if(nextLE_mem < n){\n      delete[] nextLE;\n      nextLE = new int[n];\n      nextLE_mem = n;\n    }\n    walloc1d(&st, n, &mem);\n    rrep(i,n){\n      while(s && d[st[s-1]] > d[i]) s--;\n      if(s==0){\n        nextLE[i] = n;\n      } else {\n        nextLE[i] = st[s-1];\n      }\n      st[s++] = i;\n    }\n  }\n  int NextLE(int i){\n    if(set_nextLE==0) setNextLE();\n    return nextLE[i];\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Arr1d_NextLE_reset");
      d.push_back((string)"Arr1d_NextLE_constructor");
      d.push_back((string)"Arr1d_NextLE_destructor");
      need[n] = d;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_NextLT";
      string c = "  int set_nextLT;\n  int nextLT_mem;\n  int *nextLT;\n  void setNextLT(void *mem = wmem){\n    int i;\n    int s = 0, *st;\n    set_nextLT = 1;\n    if(nextLT_mem < n){\n      delete[] nextLT;\n      nextLT = new int[n];\n      nextLT_mem = n;\n    }\n    walloc1d(&st, n, &mem);\n    rrep(i,n){\n      while(s && d[st[s-1]] >= d[i]) s--;\n      if(s==0){\n        nextLT[i] = n;\n      } else {\n        nextLT[i] = st[s-1];\n      }\n      st[s++] = i;\n    }\n  }\n  int NextLT(int i){\n    if(set_nextLT==0) setNextLT();\n    return nextLT[i];\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Arr1d_NextLT_reset");
      d.push_back((string)"Arr1d_NextLT_constructor");
      d.push_back((string)"Arr1d_NextLT_destructor");
      need[n] = d;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_NextGE";
      string c = "  int set_nextGE;\n  int nextGE_mem;\n  int *nextGE;\n  void setNextGE(void *mem = wmem){\n    int i;\n    int s = 0, *st;\n    set_nextGE = 1;\n    if(nextGE_mem < n){\n      delete[] nextGE;\n      nextGE = new int[n];\n      nextGE_mem = n;\n    }\n    walloc1d(&st, n, &mem);\n    rrep(i,n){\n      while(s && d[st[s-1]] < d[i]) s--;\n      if(s==0){\n        nextGE[i] = n;\n      } else {\n        nextGE[i] = st[s-1];\n      }\n      st[s++] = i;\n    }\n  }\n  int NextGE(int i){\n    if(set_nextGE==0) setNextGE();\n    return nextGE[i];\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Arr1d_NextGE_reset");
      d.push_back((string)"Arr1d_NextGE_constructor");
      d.push_back((string)"Arr1d_NextGE_destructor");
      need[n] = d;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_NextGT";
      string c = "  int set_nextGT;\n  int nextGT_mem;\n  int *nextGT;\n  void setNextGT(void *mem = wmem){\n    int i;\n    int s = 0, *st;\n    set_nextGT = 1;\n    if(nextGT_mem < n){\n      delete[] nextGT;\n      nextGT = new int[n];\n      nextGT_mem = n;\n    }\n    walloc1d(&st, n, &mem);\n    rrep(i,n){\n      while(s && d[st[s-1]] <= d[i]) s--;\n      if(s==0){\n        nextGT[i] = n;\n      } else {\n        nextGT[i] = st[s-1];\n      }\n      st[s++] = i;\n    }\n  }\n  int NextGT(int i){\n    if(set_nextGT==0) setNextGT();\n    return nextGT[i];\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Arr1d_NextGT_reset");
      d.push_back((string)"Arr1d_NextGT_constructor");
      d.push_back((string)"Arr1d_NextGT_destructor");
      need[n] = d;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_reset_head";
      string c = "  void reset(){\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_getSum_reset";
      string c = "    set_cumulative_sum = 0;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_ConstLenLeft_reset";
      string c = "    set_const_len_left = 0;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_ConstLenRight_reset";
      string c = "    set_const_len_right = 0;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_dHist_reset";
      string c = "    set_dhist = 0;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_sHist_reset";
      string c = "    set_shist = 0;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_PrevLE_reset";
      string c = "    set_prevLE = 0;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_PrevLT_reset";
      string c = "    set_prevLT = 0;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_PrevGE_reset";
      string c = "    set_prevGE = 0;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_PrevGT_reset";
      string c = "    set_prevGT = 0;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_NextLE_reset";
      string c = "    set_nextLE = 0;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_NextLT_reset";
      string c = "    set_nextLT = 0;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_NextGE_reset";
      string c = "    set_nextGE = 0;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_NextGT_reset";
      string c = "    set_nextGT = 0;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_reset_foot";
      string c = "  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_constructor_head";
      string c = "  void constructor(){\n    n = mem = 0;\n    d = NULL;\n    \n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_getSum_constructor";
      string c = "    set_cumulative_sum = 0;\n    cumulative_sum_mem = 0;\n    cumulative_sum = NULL;\n    \n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_ConstLenLeft_constructor";
      string c = "    set_const_len_left = 0;\n    const_len_left_mem = 0;\n    const_len_left = NULL;\n    \n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_ConstLenRight_constructor";
      string c = "    set_const_len_right = 0;\n    const_len_right_mem = 0;\n    const_len_right = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_dHist_constructor";
      string c = "    set_dhist = 0;\n    dhist_mem = 0;\n    dhist = NULL;\n    dhists = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_sHist_constructor";
      string c = "    set_shist = 0;\n    shist_mem = 0;\n    shist = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_PrevLE_constructor";
      string c = "    set_prevLE = 0;\n    prevLE_mem = 0;\n    prevLE = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_PrevLT_constructor";
      string c = "    set_prevLT = 0;\n    prevLT_mem = 0;\n    prevLT = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_PrevGE_constructor";
      string c = "    set_prevGE = 0;\n    prevGE_mem = 0;\n    prevGE = NULL;\n    \n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_PrevGT_constructor";
      string c = "    set_prevGT = 0;\n    prevGT_mem = 0;\n    prevGT = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_NextLE_constructor";
      string c = "    set_nextLE = 0;\n    nextLE_mem = 0;\n    nextLE = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_NextLT_constructor";
      string c = "    set_nextLT = 0;\n    nextLT_mem = 0;\n    nextLT = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_NextGE_constructor";
      string c = "    set_nextGE = 0;\n    nextGE_mem = 0;\n    nextGE = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_NextGT_constructor";
      string c = "    set_nextGT = 0;\n    nextGT_mem = 0;\n    nextGT = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_constructor_foot";
      string c = "  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_destructor_head";
      string c = "  void destructor(){\n    delete[] d;\n    d = NULL;\n    mem = n = 0;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_getSum_destructor";
      string c = "    set_cumulative_sum = 0;\n    cumulative_sum_mem = 0;\n    delete[] cumulative_sum;\n    cumulative_sum = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_ConstLenLeft_destructor";
      string c = "    set_const_len_left = 0;\n    const_len_left_mem = 0;\n    delete[] const_len_left;\n    const_len_left = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_ConstLenRight_destructor";
      string c = "    set_const_len_right = 0;\n    const_len_right_mem = 0;\n    delete[] const_len_right;\n    const_len_right = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_dHist_destructor";
      string c = "    set_dhist = 0;\n    dhist_mem = 0;\n    delete[] dhist;\n    delete[] dhists;\n    dhist = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_sHist_destructor";
      string c = "    set_shist = 0;\n    shist_mem = 0;\n    delete[] shist;\n    shist = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_PrevLE_destructor";
      string c = "    set_prevLE = 0;\n    prevLE_mem = 0;\n    delete[] prevLE;\n    prevLE = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_PrevLT_destructor";
      string c = "    set_prevLT = 0;\n    prevLT_mem = 0;\n    delete[] prevLT;\n    prevLT = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_PrevGE_destructor";
      string c = "    set_prevGE = 0;\n    prevGE_mem = 0;\n    delete[] prevGE;\n    prevGE = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_PrevGT_destructor";
      string c = "    set_prevGT = 0;\n    prevGT_mem = 0;\n    delete[] prevGT;\n    prevGT = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_NextLE_destructor";
      string c = "    set_nextLE = 0;\n    nextLE_mem = 0;\n    delete[] nextLE;\n    nextLE = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_NextLT_destructor";
      string c = "    set_nextLT = 0;\n    nextLT_mem = 0;\n    delete[] nextLT;\n    nextLT = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_NextGE_destructor";
      string c = "    set_nextGE = 0;\n    nextGE_mem = 0;\n    delete[] nextGE;\n    nextGE = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_NextGT_destructor";
      string c = "    set_nextGT = 0;\n    nextGT_mem = 0;\n    delete[] nextGT;\n    nextGT = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_destructor_foot";
      string c = "  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }
    {
      string n = "Arr1d_end";
      string c = "  void constructor(int nn){\n    constructor();\n    malloc(nn);\n  }\n  void memory_expand(int nn){\n    if(mem < nn){\n      delete[] d;\n      d = new T[nn];\n      mem = nn;\n    }\n  }\n  void malloc(int nn){\n    reset();\n    memory_expand(nn);\n    n = nn;\n  }\n  void setN(int nn){\n    reset();\n    memory_expand(nn);\n    n = nn;\n  }\n  void setN(int nn, T val){\n    int i;\n    reset();\n    memory_expand(nn);\n    n = nn;\n    rep(i,n) d[i] = val;\n  }\n  template<class S> void set(vector<S> &a){\n    int i, nn = a.size();\n    setN(nn);\n    rep(i,nn) d[i] = a[i];\n  }\n  template<class S> void set_c(vector<S> a){\n    int i, nn = a.size();\n    setN(nn);\n    rep(i,nn) d[i] = a[i];\n  }\n  template<class S> void set(int nn, S a[]){\n    int i;\n    setN(nn);\n    rep(i,nn) d[i] = a[i];\n  }\n  void free(){\n    destructor();\n  }\n  Arr1d(){\n    constructor();\n  }\n  Arr1d(int nn){\n    constructor(nn);\n  }\n  ~Arr1d(){\n    destructor();\n  }\n};\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr1d";
    }

    {
      string n = "Arr2d";
      string c = "template<class T>\nstruct Arr2d{\n  int n1, n2, mem1, mem2;\n  T **d;\n  T* operator[](int a){\n    return d[a];\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Arr2d_end");
      d.push_back((string)"Arr2d_reset_head");
      d.push_back((string)"Arr2d_reset_foot");
      d.push_back((string)"Arr2d_constructor_head");
      d.push_back((string)"Arr2d_constructor_foot");
      d.push_back((string)"Arr2d_destructor_head");
      d.push_back((string)"Arr2d_destructor_foot");
      need[n] = d;
    }
    {
      string n = "Arr2d_getSum";
      string c = "  int set_cumulative_sum;\n  int cumulative_sum_mem1, cumulative_sum_mem2;\n  T **cumulative_sum;\n  void setSum(void){\n    int i, j;\n    set_cumulative_sum = 1;\n    if(cumulative_sum_mem1 < n1+1 || cumulative_sum_mem2 < n2+1){\n      rep(i,cumulative_sum_mem1) delete[] cumulative_sum[i];\n      delete[] cumulative_sum;\n      cumulative_sum_mem1 = n1+1;\n      cumulative_sum_mem2 = n2+1;\n      cumulative_sum = new T*[cumulative_sum_mem1];\n      rep(i,cumulative_sum_mem1) cumulative_sum[i] = new T[cumulative_sum_mem2];\n    }\n    rep(i,n1+1) cumulative_sum[i][0] = 0;\n    rep(i,n2+1) cumulative_sum[0][i] = 0;\n    rep(i,n1) rep(j,n2) cumulative_sum[i+1][j+1] = cumulative_sum[i+1][j] + cumulative_sum[i][j+1] - cumulative_sum[i][j] + d[i][j];\n  }\n  template<class T1, class T2, class T3, class T4>\n  T getSum(T1 r1, T2 c1, T3 r2, T4 c2){\n    if(!set_cumulative_sum) setSum();\n    if(r1 > r2 || c1 > c2) return 0;\n    r1 >?= 0;\n    c1 >?= 0;\n    r2 <?= n1-1;\n    c2 <?= n2-1;\n    return cumulative_sum[r2+1][c2+1] - cumulative_sum[r2+1][c1] - cumulative_sum[r1][c2+1] + cumulative_sum[r1][c1];\n  }\n  T getSumBorder(int r1, int c1, int r2, int c2){\n    T res;\n    if(!set_cumulative_sum) setSum();\n    res = cumulative_sum[r2+1][c2+1] - cumulative_sum[r2+1][c1] - cumulative_sum[r1][c2+1] + cumulative_sum[r1][c1];\n    if(r2 - r1 > 1 && c2 - c1 > 1) res -= cumulative_sum[r2][c2] - cumulative_sum[r2][c1+1] - cumulative_sum[r1+1][c2] + cumulative_sum[r1+1][c1+1];\n    return res;\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Arr2d_getSum_reset");
      d.push_back((string)"Arr2d_getSum_constructor");
      d.push_back((string)"Arr2d_getSum_destructor");
      d.push_back((string)"chmax");
      d.push_back((string)"chmin");
      need[n] = d;
      parent[n] = "Arr2d";
    }
    {
      string n = "Arr2d_getSum45";
      string c = "  int set_cumulative_sum45;\n  int cumulative_sum45_mem;\n  T **cumulative_sum45;\n  void setSum45(void){\n    int i, j;\n    set_cumulative_sum45 = 1;\n    if(cumulative_sum45_mem < n1+n2+1){\n      rep(i,cumulative_sum45_mem) delete[] cumulative_sum45[i];\n      delete[] cumulative_sum45;\n      cumulative_sum45_mem = n1+n2+1;\n      cumulative_sum45 = new T*[cumulative_sum45_mem];\n      rep(i,cumulative_sum45_mem) cumulative_sum45[i] = new T[cumulative_sum45_mem];\n    }\n    rep(i,n1+n2+1) rep(j,n1+n2+1) cumulative_sum45[i][j] = 0;\n    rep(i,n1) rep(j,n2) cumulative_sum45[n1-i+j][i+j+1] += d[i][j];\n    rep(i,n1+n2) rep(j,n1+n2) cumulative_sum45[i+1][j+1] += cumulative_sum45[i+1][j] + cumulative_sum45[i][j+1] - cumulative_sum45[i][j];\n  }\n  T getSum45(int r1, int c1, int r2, int c2){\n    int x1, x2, y1, y2;\n    if(!set_cumulative_sum45) setSum45();\n    x1 = n1 - 1 - r1 + c1;\n    y1 = r1 + c1;\n    x2 = n1 - 1 - r2 + c2;\n    y2 = r2 + c2;\n    if(x1 > x2) swap(x1, x2);\n    if(y1 > y2) swap(y1, y2);\n    return cumulative_sum45[x2+1][y2+1] - cumulative_sum45[x2+1][y1] - cumulative_sum45[x1][y2+1] + cumulative_sum45[x1][y1];\n  }\n  T getSum45Border(int r1, int c1, int r2, int c2){\n    int x1, x2, y1, y2;\n    T res;\n    if(!set_cumulative_sum45) setSum45();\n    x1 = n1 - 1 - r1 + c1;\n    y1 = r1 + c1;\n    x2 = n1 - 1 - r2 + c2;\n    y2 = r2 + c2;\n    if(x1 > x2) swap(x1, x2);\n    if(y1 > y2) swap(y1, y2);\n    res = cumulative_sum45[x2+1][y2+1] - cumulative_sum45[x2+1][y1] - cumulative_sum45[x1][y2+1] + cumulative_sum45[x1][y1];\n    if(x2 - x1 > 1 && y2 - y1 > 1) res -= cumulative_sum45[x2][y2] - cumulative_sum45[x2][y1+1] - cumulative_sum45[x1+1][y2] + cumulative_sum45[x1+1][y1+1];\n    return res;\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Arr2d_getSum45_reset");
      d.push_back((string)"Arr2d_getSum45_constructor");
      d.push_back((string)"Arr2d_getSum45_destructor");
      need[n] = d;
      parent[n] = "Arr2d";
    }
    {
      string n = "Arr2d_ConstLenLeft";
      string c = "  int set_const_len_left;\n  int const_len_left_mem1, const_len_left_mem2;\n  int **const_len_left;\n  void setConstLenLeft(void){\n    int i, j;\n    set_const_len_left = 1;\n    if(const_len_left_mem1 < n1 || const_len_left_mem2 < n2){\n      rep(i,const_len_left_mem1) delete[] const_len_left[i];\n      delete[] const_len_left;\n      const_len_left = new int*[n1];\n      rep(i,n1) const_len_left[i] = new int[n2];\n      const_len_left_mem1 = n1;\n      const_len_left_mem2 = n2;\n    }\n    rep(i,n1) rep(j,n2) const_len_left[i][j] = 1;\n    rep(i,n1) rep(j,1,n2) if(d[i][j]==d[i][j-1]) const_len_left[i][j] = const_len_left[i][j-1] + 1;\n  }\n  int ConstLenLeft(int i, int j, T val){\n    if(!set_const_len_left) setConstLenLeft();\n    if(val != d[i][j]) return 0;\n    return const_len_left[i][j];\n  }\n  int ConstLenLeft(int i, int j){\n    if(!set_const_len_left) setConstLenLeft();\n    return const_len_left[i][j];\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Arr2d_ConstLenLeft_reset");
      d.push_back((string)"Arr2d_ConstLenLeft_constructor");
      d.push_back((string)"Arr2d_ConstLenLeft_destructor");
      need[n] = d;
      parent[n] = "Arr2d";
    }
    {
      string n = "Arr2d_ConstLenRight";
      string c = "  int set_const_len_right;\n  int const_len_right_mem1, const_len_right_mem2;\n  int **const_len_right;\n  void setConstLenRight(void){\n    int i, j;\n    set_const_len_right = 1;\n    if(const_len_right_mem1 < n1 || const_len_right_mem2 < n2){\n      rep(i,const_len_right_mem1) delete[] const_len_right[i];\n      delete[] const_len_right;\n      const_len_right = new int*[n1];\n      rep(i,n1) const_len_right[i] = new int[n2];\n      const_len_right_mem1 = n1;\n      const_len_right_mem2 = n2;\n    }\n    rep(i,n1) rep(j,n2) const_len_right[i][j] = 1;\n    rep(i,n1) rrep(j,1,n2) if(d[i][j-1]==d[i][j]) const_len_right[i][j-1] = const_len_right[i][j] + 1;\n  }\n  int ConstLenRight(int i, int j, T val){\n    if(!set_const_len_right) setConstLenRight();\n    if(val != d[i][j]) return 0;\n    return const_len_right[i][j];\n  }\n  int ConstLenRight(int i, int j){\n    if(!set_const_len_right) setConstLenRight();\n    return const_len_right[i][j];\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Arr2d_ConstLenRight_reset");
      d.push_back((string)"Arr2d_ConstLenRight_constructor");
      d.push_back((string)"Arr2d_ConstLenRight_destructor");
      need[n] = d;
      parent[n] = "Arr2d";
    }
    {
      string n = "Arr2d_ConstLenUp";
      string c = "  int set_const_len_up;\n  int const_len_up_mem1, const_len_up_mem2;\n  int **const_len_up;\n  void setConstLenUp(void){\n    int i, j;\n    set_const_len_up = 1;\n    if(const_len_up_mem1 < n1 || const_len_up_mem2 < n2){\n      rep(i,const_len_up_mem1) delete[] const_len_up[i];\n      delete[] const_len_up;\n      const_len_up = new int*[n1];\n      rep(i,n1) const_len_up[i] = new int[n2];\n      const_len_up_mem1 = n1;\n      const_len_up_mem2 = n2;\n    }\n    rep(i,n1) rep(j,n2) const_len_up[i][j] = 1;\n    rep(i,1,n1) rep(j,n2) if(d[i][j]==d[i-1][j]) const_len_up[i][j] = const_len_up[i-1][j] + 1;\n  }\n  int ConstLenUp(int i, int j, T val){\n    if(!set_const_len_up) setConstLenUp();\n    if(val != d[i][j]) return 0;\n    return const_len_up[i][j];\n  }\n  int ConstLenUp(int i, int j){\n    if(!set_const_len_up) setConstLenUp();\n    return const_len_up[i][j];\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Arr2d_ConstLenUp_reset");
      d.push_back((string)"Arr2d_ConstLenUp_constructor");
      d.push_back((string)"Arr2d_ConstLenUp_destructor");
      need[n] = d;
      parent[n] = "Arr2d";
    }
    {
      string n = "Arr2d_ConstLenDown";
      string c = "  int set_const_len_down;\n  int const_len_down_mem1, const_len_down_mem2;\n  int **const_len_down;\n  void setConstLenDown(void){\n    int i, j;\n    set_const_len_down = 1;\n    if(const_len_down_mem1 < n1 || const_len_down_mem2 < n2){\n      rep(i,const_len_down_mem1) delete[] const_len_down[i];\n      delete[] const_len_down;\n      const_len_down = new int*[n1];\n      rep(i,n1) const_len_down[i] = new int[n2];\n      const_len_down_mem1 = n1;\n      const_len_down_mem2 = n2;\n    }\n    rep(i,n1) rep(j,n2) const_len_down[i][j] = 1;\n    rrep(i,1,n1) rrep(j,n2) if(d[i-1][j]==d[i][j]) const_len_down[i-1][j] = const_len_down[i][j] + 1;\n  }\n  int ConstLenDown(int i, int j, T val){\n    if(!set_const_len_down) setConstLenDown();\n    if(val != d[i][j]) return 0;\n    return const_len_down[i][j];\n  }\n  int ConstLenDown(int i, int j){\n    if(!set_const_len_down) setConstLenDown();\n    return const_len_down[i][j];\n  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Arr2d_ConstLenDown_reset");
      d.push_back((string)"Arr2d_ConstLenDown_constructor");
      d.push_back((string)"Arr2d_ConstLenDown_destructor");
      need[n] = d;
      parent[n] = "Arr2d";
    }
    {
      string n = "Arr2d_reset_head";
      string c = "  void reset(){\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr2d";
    }
    {
      string n = "Arr2d_getSum_reset";
      string c = "    set_cumulative_sum = 0;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr2d";
    }
    {
      string n = "Arr2d_getSum45_reset";
      string c = "    set_cumulative_sum45 = 0;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr2d";
    }
    {
      string n = "Arr2d_ConstLenLeft_reset";
      string c = "    set_const_len_left = 0;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr2d";
    }
    {
      string n = "Arr2d_ConstLenRight_reset";
      string c = "    set_const_len_right = 0;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr2d";
    }
    {
      string n = "Arr2d_ConstLenUp_reset";
      string c = "    set_const_len_up = 0;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr2d";
    }
    {
      string n = "Arr2d_ConstLenDown_reset";
      string c = "    set_const_len_down = 0;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr2d";
    }
    {
      string n = "Arr2d_reset_foot";
      string c = "  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr2d";
    }
    {
      string n = "Arr2d_constructor_head";
      string c = "  void constructor(){\n    n1 = n2 = mem1 = mem2 = 0;\n    d = NULL;\n    \n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr2d";
    }
    {
      string n = "Arr2d_getSum_constructor";
      string c = "    set_cumulative_sum = 0;\n    cumulative_sum_mem1 = cumulative_sum_mem2 = 0;\n    cumulative_sum = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr2d";
    }
    {
      string n = "Arr2d_getSum45_constructor";
      string c = "    set_cumulative_sum45 = 0;\n    cumulative_sum45_mem = 0;\n    cumulative_sum45 = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr2d";
    }
    {
      string n = "Arr2d_ConstLenLeft_constructor";
      string c = "    set_const_len_left = 0;\n    const_len_left_mem1 = const_len_left_mem2 = 0;\n    const_len_left = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr2d";
    }
    {
      string n = "Arr2d_ConstLenRight_constructor";
      string c = "    set_const_len_right = 0;\n    const_len_right_mem1 = const_len_right_mem2 = 0;\n    const_len_right = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr2d";
    }
    {
      string n = "Arr2d_ConstLenUp_constructor";
      string c = "    set_const_len_up = 0;\n    const_len_up_mem1 = const_len_up_mem2 = 0;\n    const_len_up = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr2d";
    }
    {
      string n = "Arr2d_ConstLenDown_constructor";
      string c = "    set_const_len_down = 0;\n    const_len_down_mem1 = const_len_down_mem2 = 0;\n    const_len_down = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr2d";
    }
    {
      string n = "Arr2d_constructor_foot";
      string c = "  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr2d";
    }
    {
      string n = "Arr2d_destructor_head";
      string c = "  void destructor(){\n    int i;\n    \n    if(d != NULL){\n      rep(i,mem1) delete[] d[i];\n      delete[] d;\n    }\n    d = NULL;\n    mem1 = mem2 = n1 = n2 = 0;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr2d";
    }
    {
      string n = "Arr2d_getSum_destructor";
      string c = "    set_cumulative_sum = 0;\n    if(cumulative_sum != NULL){\n      rep(i,cumulative_sum_mem1) delete[] cumulative_sum[i];\n      delete[] cumulative_sum;\n    }\n    cumulative_sum_mem1 = cumulative_sum_mem2 = 0;\n    cumulative_sum = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr2d";
    }
    {
      string n = "Arr2d_getSum45_destructor";
      string c = "    set_cumulative_sum45 = 0;\n    if(cumulative_sum45 != NULL){\n      rep(i,cumulative_sum45_mem) delete[] cumulative_sum45[i];\n      delete[] cumulative_sum45;\n    }\n    cumulative_sum45_mem = 0;\n    cumulative_sum45 = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr2d";
    }
    {
      string n = "Arr2d_ConstLenLeft_destructor";
      string c = "    set_const_len_left = 0;\n    if(const_len_left != NULL){\n      rep(i,const_len_left_mem1) delete[] const_len_left[i];\n      delete[] const_len_left;\n    }\n    const_len_left_mem1 = const_len_left_mem2 = 0;\n    const_len_left = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr2d";
    }
    {
      string n = "Arr2d_ConstLenRight_destructor";
      string c = "    set_const_len_right = 0;\n    if(const_len_right != NULL){\n      rep(i,const_len_right_mem1) delete[] const_len_right[i];\n      delete[] const_len_right;\n    }\n    const_len_right_mem1 = const_len_right_mem2 = 0;\n    const_len_right = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr2d";
    }
    {
      string n = "Arr2d_ConstLenUp_destructor";
      string c = "    set_const_len_up = 0;\n    if(const_len_up != NULL){\n      rep(i,const_len_up_mem1) delete[] const_len_up[i];\n      delete[] const_len_up;\n    }\n    const_len_up_mem1 = const_len_up_mem2 = 0;\n    const_len_up = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr2d";
    }
    {
      string n = "Arr2d_ConstLenDown_destructor";
      string c = "    set_const_len_down = 0;\n    if(const_len_down != NULL){\n      rep(i,const_len_down_mem1) delete[] const_len_down[i];\n      delete[] const_len_down;\n    }\n    const_len_down_mem1 = const_len_down_mem2 = 0;\n    const_len_down = NULL;\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr2d";
    }
    {
      string n = "Arr2d_destructor_foot";
      string c = "  }\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr2d";
    }
    {
      string n = "Arr2d_end";
      string c = "  void constructor(int nn1, int nn2){\n    constructor();\n    malloc(nn1, nn2);\n  }\n  void memory_expand(int nn1, int nn2){\n    int i;\n    if(mem1 < nn1 || mem2 < nn2){\n      if(d != NULL){\n        rep(i,mem1) delete[] d[i];\n        delete[] d;\n      }\n      d = new T*[nn1];\n      rep(i,nn1) d[i] = new T[nn2];\n      mem1 = nn1;\n      mem2 = nn2;\n    }\n  }\n  void malloc(int nn1, int nn2){\n    reset();\n    memory_expand(nn1, nn2);\n    n1 = nn1;\n    n2 = nn2;\n  }\n  void setN(int nn1, int nn2){\n    reset();\n    memory_expand(nn1, nn2);\n    n1 = nn1;\n    n2 = nn2;\n  }\n  void setN(int nn1, int nn2, T val){\n    int i, j;\n    reset();\n    memory_expand(nn1, nn2);\n    n1 = nn1;\n    n2 = nn2;\n    rep(i,n1) rep(j,n2) d[i][j] = val;\n  }\n  template<class S> void set(vector<vector<S>> &a){\n    int i, j, nn1 = a.size(), nn2 = a[0].size();\n    setN(nn1, nn2);\n    rep(i,nn1) rep(j,nn2) d[i][j] = a[i][j];\n  }\n  template<class S> void set_c(vector<vector<S>> a){\n    int i, j, nn1 = a.size(), nn2 = a[0].size();\n    setN(nn1, nn2);\n    rep(i,nn1) rep(j,nn2) d[i][j] = a[i][j];\n  }\n  template<class S> void set(int nn1, int nn2, S **a){\n    int i, j;\n    setN(nn1, nn2);\n    rep(i,nn1) rep(j,nn2) d[i][j] = a[i][j];\n  }\n  void free(){\n    destructor();\n  }\n  Arr2d(){\n    constructor();\n  }\n  Arr2d(int nn1, int nn2){\n    constructor(nn1, nn2);\n  }\n  ~Arr2d(){\n    destructor();\n  }\n};\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Arr2d";
    }


    {
      string n = "Grid1d";
      string c = "template<class T>\nstruct Grid1d {\n  int n;\n  T *d;\n\n  int set_s, set_d;\n  \n  T *d_s;\n  int *up, *dw, *lf, *rg;\n\n  void malloc(const int nn){\n    n = nn;\n    set_s = 0;\n    set_d = 0;\n    malloc1d(&d, n);\n  }\n\n  void free(void){\n    free1d(d);\n    if(set_s) free1d(d_s);\n    if(set_d) free1d(up), free1d(dw);\n  }\n\n  T& operator[](int a){\n    return d[a];\n  }\n\n  void setSum(void){\n    int i;\n    if(set_s == 0){\n      set_s = 1;\n      malloc1d(&d_s, n+1);\n    }\n    d_s[0] = 0;\n    rep(i,n) d_s[i+1] = d_s[i] + d[i];\n  }\n\n  void setDir(void){\n    int i;\n    if(set_d == 0){\n      set_d = 1;\n      malloc1d(&up, n);\n      malloc1d(&dw, n);\n      lf = dw;\n      rg = up;\n    }\n\n    lf[0] = 1;\n    rep(i,1,n) lf[i] = 1 if[d[i]==d[i-1], + lf[i-1]];\n\n    rg[n-1] = 1;\n    for(i=n-2;i>=0;i--) rg[i] = 1 if[d[i]==d[i+1], + rg[i+1]];\n  }\n\n  void setDirMatch(const T v){\n    int i;\n    if(set_d == 0){\n      set_d = 1;\n      malloc1d(&up, n);\n      malloc1d(&dw, n);\n      lf = dw;\n      rg = up;\n    }\n\n    lf[0] = if[d[0]==v, 1, 0];\n    rep(i,1,n) lf[i] = if[d[i]==v, 1 + lf[i-1], 0];\n\n    rg[n-1] = if[d[n-1]==v, 1, 0];\n    for(i=n-2;i>=0;i--) rg[i] = if[d[i]==v, 1 + rg[i+1], 0];\n  }\n\n  inline T getSum(const int a, const int b){\n    return d_s[b+1] - d_s[a];\n  }\n  \n};\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"malloc1d");
      d.push_back((string)"free1d");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "Grid2d";
      string c = "template<class T>\nstruct Grid2d {\n  int r, c;\n  T **d;\n\n  int set_s, set_d;\n  \n  T **d_s;\n  int **up, **dw, **lf, **rg;\n\n  void malloc(const int rr, const int cc){\n    r = rr;\n    c = cc;\n    set_s = 0;\n    set_d = 0;\n    malloc2d(&d, r, c);\n  }\n\n  void free(void){\n    free2d(d);\n    if(set_s) free2d(d_s);\n    if(set_d) free2d(up), free2d(dw), free2d(lf), free2d(rg);\n  }\n\n  T*operator[](int a){\n    return d[a];\n  }\n\n  void setSum(void){\n    int i, j;\n    if(set_s == 0){\n      set_s = 1;\n      malloc2d(&d_s, r+1, c+1);\n    }\n    rep(i,r+1) d_s[i][0] = 0;\n    rep(j,c+1) d_s[0][j] = 0;\n    rep(i,r) rep(j,c) d_s[i+1][j+1] = d_s[i][j+1] + d_s[i+1][j] - d_s[i][j] + d[i][j];\n  }\n\n  void setDir(void){\n    int i, j;\n    if(set_d == 0){\n      set_d = 1;\n      malloc2d(&up, r, c);\n      malloc2d(&dw, r, c);\n      malloc2d(&lf, r, c);\n      malloc2d(&rg, r, c);\n    }\n\n    rep(j,c) up[0][j] = 1;\n    rep(i,1,r) rep(j,c) up[i][j] = 1 if[d[i][j]==d[i-1][j], + up[i-1][j]];\n\n    rep(j,c) dw[r-1][j] = 1;\n    for(i=r-2;i>=0;i--) rep(j,c) dw[i][j] = 1 if[d[i][j]==d[i+1][j], + dw[i+1][j]];\n\n    rep(i,r){\n      lf[i][0] = 1;\n      rep(j,1,c) lf[i][j] = 1 if[d[i][j]==d[i][j-1], + lf[i][j-1]];\n    }\n\n    rep(i,r){\n      rg[i][c-1] = 1;\n      for(j=c-2;j>=0;j--) rg[i][j] = 1 if[d[i][j]==d[i][j+1], + rg[i][j+1]];\n    }\n  }\n\n  void setDirMatch(const T v){\n    int i, j;\n    if(set_d == 0){\n      set_d = 1;\n      malloc2d(&up, r, c);\n      malloc2d(&dw, r, c);\n      malloc2d(&lf, r, c);\n      malloc2d(&rg, r, c);\n    }\n\n    rep(j,c) up[0][j] = if[d[0][j]==v, 1, 0];\n    rep(i,1,r) rep(j,c) up[i][j] = if[d[i][j]==v, 1 + up[i-1][j], 0];\n\n    rep(j,c) dw[r-1][j] = if[d[r-1][j]==v, 1, 0];\n    for(i=r-2;i>=0;i--) rep(j,c) dw[i][j] = if[d[i][j]==v, 1 + dw[i+1][j], 0];\n\n    rep(i,r){\n      lf[i][0] = if[d[i][0]==v, 1, 0];\n      rep(j,1,c) lf[i][j] = if[d[i][j]==v, 1 + lf[i][j-1], 0];\n    }\n\n    rep(i,r){\n      rg[i][c-1] = if[d[i][c-1]==v, 1, 0];\n      for(j=c-2;j>=0;j--) rg[i][j] = if[d[i][j]==v, 1 + rg[i][j+1], 0];\n    }\n  }\n\n  inline T getSum(const int r1, const int c1, const int r2, const int c2){\n    return d_s[r2+1][c2+1] - d_s[r1][c2+1] - d_s[r2+1][c1] + d_s[r1][c1];\n  }\n\n  template<class S> inline void getDist4(int sr, int sc, S **res, void *mem = wmem){\n    int i, j, k;\n    DijkstraHeap<S> hp;\n    hp.walloc(r*c);\n    hp.init(r*c);\n    if(d[sr][sc] >= 0) hp.change(sr*c+sc, d[sr][sc]);\n    while(hp.size){\n      k = hp.pop();\n      i = k / c;\n      j = k % c;\n      if(i-1 >= 0 && d[i-1][j] >= 0) hp.change((i-1)*c+j, hp.val[k]+d[i-1][j]);\n      if(i+1 <  r && d[i+1][j] >= 0) hp.change((i+1)*c+j, hp.val[k]+d[i+1][j]);\n      if(j-1 >= 0 && d[i][j-1] >= 0) hp.change(i*c+(j-1), hp.val[k]+d[i][j-1]);\n      if(j+1 <  c && d[i][j+1] >= 0) hp.change(i*c+(j+1), hp.val[k]+d[i][j+1]);\n    }\n    rep(i,r) rep(j,c) res[i][j] = if[hp.visited[i*c+j], hp.val[i*c+j], -1];\n  }\n\n  template<class S> inline void getDist4_BFS(int sr, int sc, S **res, void *mem = wmem){\n    int i, j, k;\n    int *q, qs=0, qe=0;\n    walloc1d(&q,r*c,&mem);\n    rep(i,r) rep(j,c) res[i][j] = -1;\n    if(d[sr][sc] >= 0) res[sr][sc] = 1;\n    q[qe++] = sr*c+sc;\n    while(qs < qe){\n      k = q[qs++];\n      i = k / c;\n      j = k % c;\n      if(i-1 >= 0 && d[i-1][j] >= 0 && res[i-1][j]==-1) res[i-1][j] = res[i][j] + 1, q[qe++] = (i-1)*c + j;\n      if(i+1 <  r && d[i+1][j] >= 0 && res[i+1][j]==-1) res[i+1][j] = res[i][j] + 1, q[qe++] = (i+1)*c + j;\n      if(j-1 >= 0 && d[i][j-1] >= 0 && res[i][j-1]==-1) res[i][j-1] = res[i][j] + 1, q[qe++] = i*c + (j-1);\n      if(j+1 <  c && d[i][j+1] >= 0 && res[i][j+1]==-1) res[i][j+1] = res[i][j] + 1, q[qe++] = i*c + (j+1);\n    }\n  }\n};\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"malloc2d");
      d.push_back((string)"free2d");
      d.push_back((string)"DijkstraHeap");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "fft";
      string c = "struct fft_pnt{\n  double x, y;\n\n  fft_pnt(void){\n  }\n  fft_pnt(double a, double b){\n    x = a;\n    y = b;\n  }\n\n  void set(double a, double b){\n    x = a;\n    y = b;\n  }\n\n  fft_pnt& operator+=(fft_pnt a){ x+=a.x; y+=a.y; return *this; }\n  fft_pnt& operator-=(fft_pnt a){ x-=a.x; y-=a.y; return *this; }\n  fft_pnt& operator*=(fft_pnt a){ fft_pnt p = *this; x = p.x*a.x-p.y*a.y; y = p.x*a.y+p.y*a.x; return *this; }\n\n  fft_pnt operator+(fft_pnt a){ return fft_pnt(*this) += a; }\n  fft_pnt operator-(fft_pnt a){ return fft_pnt(*this) -= a; }\n  fft_pnt operator*(fft_pnt a){ return fft_pnt(*this) *= a; }\n};\n\nvoid fft(int n, fft_pnt x[], void *mem = wmem){\n  int i, j;\n  int n1, n2, n3, step = 1;\n  double theta = 2*PI / n, tmp;\n  fft_pnt w1, w2, w3, a, b, c, d, aa, bb, cc, dd, *y = (fft_pnt*)mem;\n\n  while(n > 2){\n    n1 = n / 4;\n    n2 = n1 + n1;\n    n3 = n1 + n2;\n    rep(i,n1){\n      w1 = fft_pnt(cos(i*theta),-sin(i*theta));\n      w2 = w1*w1;\n      w3 = w1*w2;\n      rep(j,step){\n        a = x[j+step*i];\n        b = x[j+step*(i+n1)];\n        c = x[j+step*(i+n2)];\n        d = x[j+step*(i+n3)];\n        aa = a + c;\n        bb = a - c;\n        cc = b + d;\n        dd = b - d;\n        tmp = dd.y; dd.y = dd.x; dd.x = -tmp;\n        y[j+step*(4*i  )] = aa + cc;\n        y[j+step*(4*i+1)] = w1*(bb - dd);\n        y[j+step*(4*i+2)] = w2*(aa - cc);\n        y[j+step*(4*i+3)] = w3*(bb + dd);\n      }\n    }\n    n /= 4;\n    step *= 4;\n    theta *= 4;\n    swap(x,y);\n  }\n\n  if(n==2){\n    rep(i,step){\n      y[i] = x[i] + x[i+step];\n      y[i+step] = x[i] - x[i+step];\n    }\n    n /= 2;\n    step *= 2;\n    theta *= 2;\n    swap(x,y);\n  }\n  \n  rep(i,step) y[i] = x[i];\n}\n\nvoid fftinv(int n, fft_pnt x[], void *mem = wmem){\n  int i, j;\n  int n1, n2, n3, step = 1;\n  double theta = 2*PI / n, tmp;\n  fft_pnt w1, w2, w3, a, b, c, d, aa, bb, cc, dd, *y = (fft_pnt*)mem;\n\n  while(n > 2){\n    n1 = n / 4;\n    n2 = n1 + n1;\n    n3 = n1 + n2;\n    rep(i,n1){\n      w1 = fft_pnt(cos(i*theta),sin(i*theta));\n      w2 = w1*w1;\n      w3 = w1*w2;\n      rep(j,step){\n        a = x[j+step*i];\n        b = x[j+step*(i+n1)];\n        c = x[j+step*(i+n2)];\n        d = x[j+step*(i+n3)];\n        aa = a + c;\n        bb = a - c;\n        cc = b + d;\n        dd = b - d;\n        tmp = dd.y; dd.y = dd.x; dd.x = -tmp;\n        y[j+step*(4*i  )] = aa + cc;\n        y[j+step*(4*i+1)] = w1*(bb + dd);\n        y[j+step*(4*i+2)] = w2*(aa - cc);\n        y[j+step*(4*i+3)] = w3*(bb - dd);\n      }\n    }\n    n /= 4;\n    step *= 4;\n    theta *= 4;\n    swap(x,y);\n  }\n\n  if(n==2){\n    rep(i,step){\n      y[i] = x[i] + x[i+step];\n      y[i+step] = x[i] - x[i+step];\n    }\n    n /= 2;\n    step *= 2;\n    theta *= 2;\n    swap(x,y);\n  }\n  \n  rep(i,step) y[i] = x[i];\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"define_PI");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "convolution";
      string c = "void convolution_L(double A[], int As, double B[], int Bs, double res[], int Rs, void *mem = wmem){\n  int i, n, n2;\n  double mul;\n  fft_pnt *a, *b;\n\n  n = max(As+Bs, Rs);\n  for(n2=1;n2<n;n2*=2);\n  walloc1d(&a, n2, &mem);\n  walloc1d(&b, n2, &mem);\n\n  rep(i,As) a[i].set(A[i], 0);\n  REP(i,As,n2) a[i].set(0,0);\n  rep(i,Bs) b[i].set(B[i], 0);\n  REP(i,Bs,n2) b[i].set(0,0);\n\n  fft(n2, a, mem);\n  fft(n2, b, mem);\n  rep(i,n2) a[i] *= b[i];\n  fftinv(n2, a, mem);\n  mul = 1.0 / n2;\n  rep(i,Rs) res[i] = a[i].x * mul;\n}\n\nvoid convolution_L(double A[], int As, double res[], int Rs, void *mem = wmem){\n  int i, n, n2;\n  double mul;\n  fft_pnt *a;\n\n  n = max(As+As, Rs);\n  for(n2=1;n2<n;n2*=2);\n  walloc1d(&a, n2, &mem);\n\n  rep(i,As) a[i].set(A[i], 0);\n  REP(i,As,n2) a[i].set(0,0);\n\n  fft(n2, a, mem);\n  rep(i,n2) a[i] *= a[i];\n  fftinv(n2, a, mem);\n  mul = 1.0 / n2;\n  rep(i,Rs) res[i] = a[i].x * mul;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"max_L");
      d.push_back((string)"fft");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "fft-mint";
      string c = "void fft(int n, mint x[], mint root = MD_PRIMITIVE_ROOT, void *mem = wmem){\n  int i, j;\n  int n1, n2, n3, step = 1;\n  mint w1, w2, w3, a, b, c, d, aa, bb, cc, dd, tmp, *y;\n  walloc1d(&y, n, &mem);\n\n  tmp = root.pw((mint::md-1)/4*3);\n  root = root.pw((mint::md-1)/n);\n\n  while(n > 2){\n    n1 = n / 4;\n    n2 = n1 + n1;\n    n3 = n1 + n2;\n    w1.val = mint::R;\n    rep(i,n1){\n      w2 = w1*w1;\n      w3 = w1*w2;\n      rep(j,step){\n        a = x[j+step*i];\n        b = x[j+step*(i+n1)];\n        c = x[j+step*(i+n2)];\n        d = x[j+step*(i+n3)];\n        aa = a + c;\n        bb = a - c;\n        cc = b + d;\n        dd = (b - d) * tmp;\n        y[j+step*(4*i  )] = aa + cc;\n        y[j+step*(4*i+1)] = w1*(bb - dd);\n        y[j+step*(4*i+2)] = w2*(aa - cc);\n        y[j+step*(4*i+3)] = w3*(bb + dd);\n      }\n      w1 *= root;\n    }\n    n /= 4;\n    step *= 4;\n    root *= root;\n    root *= root;\n    swap(x,y);\n  }\n\n  if(n==2){\n    rep(i,step){\n      y[i] = x[i] + x[i+step];\n      y[i+step] = x[i] - x[i+step];\n    }\n    n /= 2;\n    step *= 2;\n    root *= root;\n    swap(x,y);\n  }\n  \n  rep(i,step) y[i] = x[i];\n}\n\nvoid fftinv(int n, mint x[], mint root = MD_PRIMITIVE_ROOT, void *mem = wmem){\n  int i, j;\n  int n1, n2, n3, step = 1;\n  mint w1, w2, w3, a, b, c, d, aa, bb, cc, dd, tmp, *y;\n  walloc1d(&y, n, &mem);\n\n  root = root.inverse();\n  tmp = root.pw((mint::md-1)/4);\n  root = root.pw((mint::md-1)/n);\n\n  while(n > 2){\n    n1 = n / 4;\n    n2 = n1 + n1;\n    n3 = n1 + n2;\n    w1.val = mint::R;\n    rep(i,n1){\n      w2 = w1*w1;\n      w3 = w1*w2;\n      rep(j,step){\n        a = x[j+step*i];\n        b = x[j+step*(i+n1)];\n        c = x[j+step*(i+n2)];\n        d = x[j+step*(i+n3)];\n        aa = a + c;\n        bb = a - c;\n        cc = b + d;\n        dd = (b - d) * tmp;\n        y[j+step*(4*i  )] = aa + cc;\n        y[j+step*(4*i+1)] = w1*(bb + dd);\n        y[j+step*(4*i+2)] = w2*(aa - cc);\n        y[j+step*(4*i+3)] = w3*(bb - dd);\n      }\n      w1 *= root;\n    }\n    n /= 4;\n    step *= 4;\n    root *= root;\n    root *= root;\n    swap(x,y);\n  }\n\n  if(n==2){\n    rep(i,step){\n      y[i] = x[i] + x[i+step];\n      y[i+step] = x[i] - x[i+step];\n    }\n    n /= 2;\n    step *= 2;\n    root *= root;\n    swap(x,y);\n  }\n  \n  rep(i,step) y[i] = x[i];\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"define_MD_PRIMITIVE_ROOT");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p; parent[n] = "mint";
    }

    {
      string n = "convolution-mint";
      string c = "void convolution_L(mint A[], int As, mint B[], int Bs, mint res[], int Rs,  mint root = MD_PRIMITIVE_ROOT, void *mem = wmem){\n  int i, n, n2;\n  mint *a, *b, r;\n\n  n = max(As+Bs, Rs);\n  for(n2=1;n2<n;n2*=2);\n  walloc1d(&a, n2, &mem);\n  walloc1d(&b, n2, &mem);\n\n  rep(i,As) a[i] = A[i];\n  REP(i,As,n2) a[i].val = 0;\n  rep(i,Bs) b[i] = B[i];\n  REP(i,Bs,n2) b[i].val = 0;\n\n  fft(n2, a, root, mem);\n  fft(n2, b, root, mem);\n  rep(i,n2) a[i] *= b[i];\n  fftinv(n2, a, root, mem);\n  r = mint(n2).inverse();\n  rep(i,Rs) res[i] = a[i] * r;\n}\n\nvoid convolution_L(mint A[], int As, mint res[], int Rs, mint root = MD_PRIMITIVE_ROOT, void *mem = wmem){\n  int i, n, n2;\n  mint *a, r;\n\n  n = max(2*As, Rs);\n  for(n2=1;n2<n;n2*=2);\n  walloc1d(&a, n2, &mem);\n\n  rep(i,As) a[i] = A[i];\n  REP(i,As,n2) a[i].val = 0;\n\n  fft(n2, a, root, mem);\n  rep(i,n2) a[i] *= a[i];\n  fftinv(n2, a, root, mem);\n  r = mint(n2).inverse();\n  rep(i,Rs) res[i] = a[i]*r;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"fft-mint");
      d.push_back((string)"max_L");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p; parent[n] = "mint";
    }

    {
      string n = "fft-Mint";
      string c = "void fft(int n, Mint x[], Mint root = MD_PRIMITIVE_ROOT, void *mem = wmem){\n  int i, j;\n  int n1, n2, n3, step = 1;\n  Mint w1, w2, w3, a, b, c, d, aa, bb, cc, dd, tmp, *y;\n  walloc1d(&y, n, &mem);\n\n  tmp = root.pw((MD-1)/4*3);\n  root = root.pw((MD-1)/n);\n\n  while(n > 2){\n    n1 = n / 4;\n    n2 = n1 + n1;\n    n3 = n1 + n2;\n    w1.val = MINT_R;\n    rep(i,n1){\n      w2 = w1*w1;\n      w3 = w1*w2;\n      rep(j,step){\n        a = x[j+step*i];\n        b = x[j+step*(i+n1)];\n        c = x[j+step*(i+n2)];\n        d = x[j+step*(i+n3)];\n        aa = a + c;\n        bb = a - c;\n        cc = b + d;\n        dd = (b - d) * tmp;\n        y[j+step*(4*i  )] = aa + cc;\n        y[j+step*(4*i+1)] = w1*(bb - dd);\n        y[j+step*(4*i+2)] = w2*(aa - cc);\n        y[j+step*(4*i+3)] = w3*(bb + dd);\n      }\n      w1 *= root;\n    }\n    n /= 4;\n    step *= 4;\n    root *= root;\n    root *= root;\n    swap(x,y);\n  }\n\n  if(n==2){\n    rep(i,step){\n      y[i] = x[i] + x[i+step];\n      y[i+step] = x[i] - x[i+step];\n    }\n    n /= 2;\n    step *= 2;\n    root *= root;\n    swap(x,y);\n  }\n  \n  rep(i,step) y[i] = x[i];\n}\n\nvoid fftinv(int n, Mint x[], Mint root = MD_PRIMITIVE_ROOT, void *mem = wmem){\n  int i, j;\n  int n1, n2, n3, step = 1;\n  Mint w1, w2, w3, a, b, c, d, aa, bb, cc, dd, tmp, *y;\n  walloc1d(&y, n, &mem);\n\n  root = root.inverse();\n  tmp = root.pw((MD-1)/4);\n  root = root.pw((MD-1)/n);\n\n  while(n > 2){\n    n1 = n / 4;\n    n2 = n1 + n1;\n    n3 = n1 + n2;\n    w1.val = MINT_R;\n    rep(i,n1){\n      w2 = w1*w1;\n      w3 = w1*w2;\n      rep(j,step){\n        a = x[j+step*i];\n        b = x[j+step*(i+n1)];\n        c = x[j+step*(i+n2)];\n        d = x[j+step*(i+n3)];\n        aa = a + c;\n        bb = a - c;\n        cc = b + d;\n        dd = (b - d) * tmp;\n        y[j+step*(4*i  )] = aa + cc;\n        y[j+step*(4*i+1)] = w1*(bb + dd);\n        y[j+step*(4*i+2)] = w2*(aa - cc);\n        y[j+step*(4*i+3)] = w3*(bb - dd);\n      }\n      w1 *= root;\n    }\n    n /= 4;\n    step *= 4;\n    root *= root;\n    root *= root;\n    swap(x,y);\n  }\n\n  if(n==2){\n    rep(i,step){\n      y[i] = x[i] + x[i+step];\n      y[i+step] = x[i] - x[i+step];\n    }\n    n /= 2;\n    step *= 2;\n    root *= root;\n    swap(x,y);\n  }\n  \n  rep(i,step) y[i] = x[i];\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"define_MD_PRIMITIVE_ROOT");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p; parent[n] = "Mint";
    }

    {
      string n = "convolution-Mint";
      string c = "void convolution_L(Mint A[], int As, Mint B[], int Bs, Mint res[], int Rs,  Mint root = MD_PRIMITIVE_ROOT, void *mem = wmem){\n  int i, n, n2;\n  Mint *a, *b, r;\n\n  n = max(As+Bs, Rs);\n  for(n2=1;n2<n;n2*=2);\n  walloc1d(&a, n2, &mem);\n  walloc1d(&b, n2, &mem);\n\n  rep(i,As) a[i] = A[i];\n  REP(i,As,n2) a[i].val = 0;\n  rep(i,Bs) b[i] = B[i];\n  REP(i,Bs,n2) b[i].val = 0;\n\n  fft(n2, a, root, mem);\n  fft(n2, b, root, mem);\n  rep(i,n2) a[i] *= b[i];\n  fftinv(n2, a, root, mem);\n  r = Mint(n2).inverse();\n  rep(i,Rs) res[i] = a[i] * r;\n}\n\nvoid convolution_L(Mint A[], int As, Mint res[], int Rs, Mint root = MD_PRIMITIVE_ROOT, void *mem = wmem){\n  int i, n, n2;\n  Mint *a, r;\n\n  n = max(2*As, Rs);\n  for(n2=1;n2<n;n2*=2);\n  walloc1d(&a, n2, &mem);\n\n  rep(i,As) a[i] = A[i];\n  REP(i,As,n2) a[i].val = 0;\n\n  fft(n2, a, root, mem);\n  rep(i,n2) a[i] *= a[i];\n  fftinv(n2, a, root, mem);\n  r = Mint(n2).inverse();\n  rep(i,Rs) res[i] = a[i]*r;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"fft-Mint");
      d.push_back((string)"max_L");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p; parent[n] = "Mint";
    }

    {
      string n = "fft-modint";
      string c = "void fft(int n, modint x[], modint root = MD_PRIMITIVE_ROOT, void *mem = wmem){\n  int i, j;\n  int n1, n2, n3, step = 1;\n  modint w1, w2, w3, a, b, c, d, aa, bb, cc, dd, tmp, *y;\n  walloc1d(&y, n, &mem);\n\n  tmp = root.pw((modint::md-1)/4*3);\n  root = root.pw((modint::md-1)/n);\n\n  while(n > 2){\n    n1 = n / 4;\n    n2 = n1 + n1;\n    n3 = n1 + n2;\n    w1.val = 1;\n    rep(i,n1){\n      w2 = w1*w1;\n      w3 = w1*w2;\n      rep(j,step){\n        a = x[j+step*i];\n        b = x[j+step*(i+n1)];\n        c = x[j+step*(i+n2)];\n        d = x[j+step*(i+n3)];\n        aa = a + c;\n        bb = a - c;\n        cc = b + d;\n        dd = (b - d) * tmp;\n        y[j+step*(4*i  )] = aa + cc;\n        y[j+step*(4*i+1)] = w1*(bb - dd);\n        y[j+step*(4*i+2)] = w2*(aa - cc);\n        y[j+step*(4*i+3)] = w3*(bb + dd);\n      }\n      w1 *= root;\n    }\n    n /= 4;\n    step *= 4;\n    root *= root;\n    root *= root;\n    swap(x,y);\n  }\n\n  if(n==2){\n    rep(i,step){\n      y[i] = x[i] + x[i+step];\n      y[i+step] = x[i] - x[i+step];\n    }\n    n /= 2;\n    step *= 2;\n    root *= root;\n    swap(x,y);\n  }\n  \n  rep(i,step) y[i] = x[i];\n}\n\nvoid fftinv(int n, modint x[], modint root = MD_PRIMITIVE_ROOT, void *mem = wmem){\n  int i, j;\n  int n1, n2, n3, step = 1;\n  modint w1, w2, w3, a, b, c, d, aa, bb, cc, dd, tmp, *y;\n  walloc1d(&y, n, &mem);\n\n  root = root.inverse();\n  tmp = root.pw((modint::md-1)/4);\n  root = root.pw((modint::md-1)/n);\n\n  while(n > 2){\n    n1 = n / 4;\n    n2 = n1 + n1;\n    n3 = n1 + n2;\n    w1.val = 1;\n    rep(i,n1){\n      w2 = w1*w1;\n      w3 = w1*w2;\n      rep(j,step){\n        a = x[j+step*i];\n        b = x[j+step*(i+n1)];\n        c = x[j+step*(i+n2)];\n        d = x[j+step*(i+n3)];\n        aa = a + c;\n        bb = a - c;\n        cc = b + d;\n        dd = (b - d) * tmp;\n        y[j+step*(4*i  )] = aa + cc;\n        y[j+step*(4*i+1)] = w1*(bb + dd);\n        y[j+step*(4*i+2)] = w2*(aa - cc);\n        y[j+step*(4*i+3)] = w3*(bb - dd);\n      }\n      w1 *= root;\n    }\n    n /= 4;\n    step *= 4;\n    root *= root;\n    root *= root;\n    swap(x,y);\n  }\n\n  if(n==2){\n    rep(i,step){\n      y[i] = x[i] + x[i+step];\n      y[i+step] = x[i] - x[i+step];\n    }\n    n /= 2;\n    step *= 2;\n    root *= root;\n    swap(x,y);\n  }\n  \n  rep(i,step) y[i] = x[i];\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"define_MD_PRIMITIVE_ROOT");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p; parent[n] = "modint";
    }

    {
      string n = "convolution-modint";
      string c = "void convolution_L(modint A[], int As, modint B[], int Bs, modint res[], int Rs,  modint root = MD_PRIMITIVE_ROOT, void *mem = wmem){\n  int i, n, n2;\n  modint *a, *b, r;\n\n  n = max(As+Bs, Rs);\n  for(n2=1;n2<n;n2*=2);\n  walloc1d(&a, n2, &mem);\n  walloc1d(&b, n2, &mem);\n\n  rep(i,As) a[i] = A[i];\n  REP(i,As,n2) a[i].val = 0;\n  rep(i,Bs) b[i] = B[i];\n  REP(i,Bs,n2) b[i].val = 0;\n\n  fft(n2, a, root, mem);\n  fft(n2, b, root, mem);\n  rep(i,n2) a[i] *= b[i];\n  fftinv(n2, a, root, mem);\n  r = modint(n2).inverse();\n  rep(i,Rs) res[i] = a[i] * r;\n}\n\nvoid convolution_L(modint A[], int As, modint res[], int Rs, modint root = MD_PRIMITIVE_ROOT, void *mem = wmem){\n  int i, n, n2;\n  modint *a, r;\n\n  n = max(2*As, Rs);\n  for(n2=1;n2<n;n2*=2);\n  walloc1d(&a, n2, &mem);\n\n  rep(i,As) a[i] = A[i];\n  REP(i,As,n2) a[i].val = 0;\n\n  fft(n2, a, root, mem);\n  rep(i,n2) a[i] *= a[i];\n  fftinv(n2, a, root, mem);\n  r = modint(n2).inverse();\n  rep(i,Rs) res[i] = a[i]*r;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"fft-modint");
      d.push_back((string)"max_L");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p; parent[n] = "modint";
    }

    {
      string n = "fft-Modint";
      string c = "void fft(int n, Modint x[], Modint root = MD_PRIMITIVE_ROOT, void *mem = wmem){\n  int i, j;\n  int n1, n2, n3, step = 1;\n  Modint w1, w2, w3, a, b, c, d, aa, bb, cc, dd, tmp, *y;\n  walloc1d(&y, n, &mem);\n\n  tmp = root.pw((MD-1)/4*3);\n  root = root.pw((MD-1)/n);\n\n  while(n > 2){\n    n1 = n / 4;\n    n2 = n1 + n1;\n    n3 = n1 + n2;\n    w1.val = 1;\n    rep(i,n1){\n      w2 = w1*w1;\n      w3 = w1*w2;\n      rep(j,step){\n        a = x[j+step*i];\n        b = x[j+step*(i+n1)];\n        c = x[j+step*(i+n2)];\n        d = x[j+step*(i+n3)];\n        aa = a + c;\n        bb = a - c;\n        cc = b + d;\n        dd = (b - d) * tmp;\n        y[j+step*(4*i  )] = aa + cc;\n        y[j+step*(4*i+1)] = w1*(bb - dd);\n        y[j+step*(4*i+2)] = w2*(aa - cc);\n        y[j+step*(4*i+3)] = w3*(bb + dd);\n      }\n      w1 *= root;\n    }\n    n /= 4;\n    step *= 4;\n    root *= root;\n    root *= root;\n    swap(x,y);\n  }\n\n  if(n==2){\n    rep(i,step){\n      y[i] = x[i] + x[i+step];\n      y[i+step] = x[i] - x[i+step];\n    }\n    n /= 2;\n    step *= 2;\n    root *= root;\n    swap(x,y);\n  }\n  \n  rep(i,step) y[i] = x[i];\n}\n\nvoid fftinv(int n, Modint x[], Modint root = MD_PRIMITIVE_ROOT, void *mem = wmem){\n  int i, j;\n  int n1, n2, n3, step = 1;\n  Modint w1, w2, w3, a, b, c, d, aa, bb, cc, dd, tmp, *y;\n  walloc1d(&y, n, &mem);\n\n  root = root.inverse();\n  tmp = root.pw((MD-1)/4);\n  root = root.pw((MD-1)/n);\n\n  while(n > 2){\n    n1 = n / 4;\n    n2 = n1 + n1;\n    n3 = n1 + n2;\n    w1.val = 1;\n    rep(i,n1){\n      w2 = w1*w1;\n      w3 = w1*w2;\n      rep(j,step){\n        a = x[j+step*i];\n        b = x[j+step*(i+n1)];\n        c = x[j+step*(i+n2)];\n        d = x[j+step*(i+n3)];\n        aa = a + c;\n        bb = a - c;\n        cc = b + d;\n        dd = (b - d) * tmp;\n        y[j+step*(4*i  )] = aa + cc;\n        y[j+step*(4*i+1)] = w1*(bb + dd);\n        y[j+step*(4*i+2)] = w2*(aa - cc);\n        y[j+step*(4*i+3)] = w3*(bb - dd);\n      }\n      w1 *= root;\n    }\n    n /= 4;\n    step *= 4;\n    root *= root;\n    root *= root;\n    swap(x,y);\n  }\n\n  if(n==2){\n    rep(i,step){\n      y[i] = x[i] + x[i+step];\n      y[i+step] = x[i] - x[i+step];\n    }\n    n /= 2;\n    step *= 2;\n    root *= root;\n    swap(x,y);\n  }\n  \n  rep(i,step) y[i] = x[i];\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"define_MD_PRIMITIVE_ROOT");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p; parent[n] = "Modint";
    }

    {
      string n = "convolution-Modint";
      string c = "void convolution_L(Modint A[], int As, Modint B[], int Bs, Modint res[], int Rs,  Modint root = MD_PRIMITIVE_ROOT, void *mem = wmem){\n  int i, n, n2;\n  Modint *a, *b, r;\n\n  n = max(As+Bs, Rs);\n  for(n2=1;n2<n;n2*=2);\n  walloc1d(&a, n2, &mem);\n  walloc1d(&b, n2, &mem);\n\n  rep(i,As) a[i] = A[i];\n  REP(i,As,n2) a[i].val = 0;\n  rep(i,Bs) b[i] = B[i];\n  REP(i,Bs,n2) b[i].val = 0;\n\n  fft(n2, a, root, mem);\n  fft(n2, b, root, mem);\n  rep(i,n2) a[i] *= b[i];\n  fftinv(n2, a, root, mem);\n  r = Modint(n2).inverse();\n  rep(i,Rs) res[i] = a[i] * r;\n}\n\nvoid convolution_L(Modint A[], int As, Modint res[], int Rs, Modint root = MD_PRIMITIVE_ROOT, void *mem = wmem){\n  int i, n, n2;\n  Modint *a, r;\n\n  n = max(2*As, Rs);\n  for(n2=1;n2<n;n2*=2);\n  walloc1d(&a, n2, &mem);\n\n  rep(i,As) a[i] = A[i];\n  REP(i,As,n2) a[i].val = 0;\n\n  fft(n2, a, root, mem);\n  rep(i,n2) a[i] *= a[i];\n  fftinv(n2, a, root, mem);\n  r = Modint(n2).inverse();\n  rep(i,Rs) res[i] = a[i]*r;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"fft-Modint");
      d.push_back((string)"max_L");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p; parent[n] = "Modint";
    }

    {
      string n = "Hungarian";
      string c = "template<class T>\nT Hungarian(T **mat, int n, int m, int match[] = NULL, void *mem = wmem){\n  int i, a, b, c, r, z;\n  int *toright;\n  int *toleft;\n  T *ofsleft;\n  T *ofsright;\n  int *left, *right;\n  int *trace, *ptr;\n  T d, t, res = 0;\n\n  walloc1d(&toright, n, &mem);\n  walloc1d(&toleft, m, &mem);\n  walloc1d(&ofsleft, n, &mem);\n  walloc1d(&ofsright, m, &mem);\n\n  walloc1d(&left, n, &mem);\n  walloc1d(&right, m, &mem);\n  walloc1d(&trace, m, &mem);\n  walloc1d(&ptr, m, &mem);\n\n  rep(i,n) toright[i] = -1, ofsleft[i] = 0;\n  rep(i,m) toleft[i] = -1, ofsright[i] = 0;\n\n  rep(r,n){\n    rep(i,n) left[i] = 0;\n    rep(i,m) right[i] = 0;\n    rep(i,m) trace[i] = -1, ptr[i] = r;\n    left[r] = 1;\n\n    for(;;){\n      d = std::numeric_limits<T>::max();\n      rep(i,m) if(!right[i]){\n        t = mat[ptr[i]][i] + ofsleft[ptr[i]] + ofsright[i];\n        if(d > t) d = t, b = i;\n      }\n\n      res += d;\n      rep(i,n) if(left[i]) ofsleft[i] -= d;\n      rep(i,m) if(right[i]) ofsright[i] += d;\n\n      trace[b] = ptr[b];\n      c = toleft[b];\n      if(c < 0){\n        while(b>=0){\n          a = trace[b];\n          z = toright[a];\n          toleft[b] = a;\n          toright[a] = b;\n          b = z;\n        }\n        break;\n      }\n      right[b] = left[c] = 1;\n      rep(i,m) if(mat[c][i] + ofsleft[c] + ofsright[i] < mat[ptr[i]][i] + ofsleft[ptr[i]] + ofsright[i]) ptr[i] = c;\n    }\n  }\n\n  if(match!=NULL) rep(i,n) match[i] = toright[i];\n  return res;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "polationVal";
      string c = "template<class T>\nT polationVal(int n, T x[], T y[], T t){\n  int i, j;\n  T res, tmp;\n  rep(i,n) if(x[i]==t) return y[i];\n  res = 0;\n  rep(i,n){\n    tmp = 1;\n    rep(j,n) if(i!=j) tmp *= x[i] - x[j];\n    tmp = 1 / tmp;\n    rep(j,n) if(i!=j) tmp *= t - x[j];\n    res += y[i] * tmp;\n  }\n  return res;\n}\ntemplate<class T>\nT polationVal(int n, T y[], T t, void *mem = wmem){\n  int i, j;\n  T res, m, *ifac, *lf, *rg;\n  if(n==1) return y[0];\n  rep(i,n) if(t==i) return y[i];\n  walloc1d(&ifac, n+1);\n  walloc1d(&lf, n+1);\n  walloc1d(&rg, n+1);\n  ifac[0] = ifac[1] = m = 1;\n  rep(i,2,n+1) m *= i;\n  m = ifac[n] = 1 / m;\n  rrep(i,2,n) ifac[i] = ifac[i+1] * (i+1);\n  lf[0] = 1;\n  rep(i,n) lf[i+1] = lf[i] * (t-i);\n  rg[0] = 1;\n  rep(i,n) rg[i+1] = rg[i] * (t-(n-1-i));\n  if(n%2==0) m = -m;\n  res = 0;\n  rep(i,n){\n    m *= i - n;\n    res -= y[i] * lf[i] * rg[n-1-i] * m * ifac[i];\n  }\n  \n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"workmemory");
      d.push_back((string)"walloc1d");
      need[n] = d;
    }


    {
      string n = "polationPoly";
      string c = "template<class T>\nPolynomial<T> polationPoly_L(int n, T x[], T y[]){\n  int i, j;\n  T c;\n  Polynomial<T> res, tmp, t1;\n\n  tmp.change(0, 1);\n  t1.change(1, 1);\n  rep(i,n){\n    t1.change(0, -x[i]);\n    tmp *= t1;\n  }\n\n  rep(i,n){\n    c = 1;\n    rep(j,n) if(j!=i) c *= (x[i] - x[j]);\n    c = y[i] / c;\n\n    t1.change(0, -x[i]);\n    res += c * tmp / t1;\n  }\n\n  return res;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"Polynomial");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "Explode";
      string c = "vector<string> Explode_L(const string &str, const string &d){\n  int s = 0, i = 0, j;\n  vector<string> res;\n  while(i + d.size() - 1 < str.size()){\n    rep(j,d.size()) if(str[i+j] != d[j]) break;\n    if(j != d.size()) i++, continue;\n    res.push_back(str.substr(s, i-s));\n    s = (i += d.size());\n  }\n  res.push_back(str.substr(s));\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "Implode";
      string c = "string Implode_L(const vector<string> &v, const string &d){\n  int i;\n  string res;\n  if(v.size()==0) return res;\n  res += v[0];\n  rep(i,1,v.size()){\n    res += d;\n    res += v[i];\n  }\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "knightDistance";
      string c = "template<class T>\ninline T knightDistance(T x, T y){\n  T res;\n  if(x < 0) x = -x;\n  if(y < 0) y = -y;\n  if(x+y==1) return 3;\n  if(x==y==2) return 4;\n  res = max( x/+2, y/+2, (x+y)/+3 );\n  if(res%2 != (x+y)%2) res++;\n  return res;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"max_L");
      d.push_back((string)"divup");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "InnerProd_array";
      string c = "template<class S>\nS InnerProd(int n, S a[]){\n  S res = 0;\n  rep(i,n) res += a[i];\n  return res;\n}\ntemplate<class S, class T>\nS InnerProd(int n, S a[], T b[]){\n  S res = 0;\n  rep(i,n) res += a[i] * b[i];\n  return res;\n}\ntemplate<class S, class T, class U>\nS InnerProd(int n, S a[], T b[], U c[]){\n  S res = 0;\n  rep(i,n) res += a[i] * b[i] * c[i];\n  return res;\n}\ntemplate<class S, class T, class U, class V>\nS InnerProd(int n, S a[], T b[], U c[], V d[]){\n  S res = 0;\n  rep(i,n) res += a[i] * b[i] * c[i] * d[i];\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }


    {
      string n = "crossProd";
      string c = "ll crossProd_L(ll x1, ll y1, ll x2, ll y2){\n  return x1 * y2 - x2 * y1;\n}\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "LineIntersection_size";
      string c = "int LineIntersection_size_L(ll x1, ll y1, ll x2, ll y2, ll x3, ll y3, ll x4, ll y4, int p1 = 1, int p2 = 1, int p3 = 1, int p4 = 1){\n  ll a, b, dx1, dy1, dx2, dy2;\n\n  if(p1 == p3 == 2){\n    dx1 = x2 - x1;\n    dy1 = y2 - y1;\n    if(dx1 < 0) dx1 = -dx1, dy1 = -dy1;\n    dx2 = x4 - x3;\n    dy2 = y4 - y3;\n    if(dx2 < 0) dx2 = -dx2, dy2 = -dy2;\n    a = gcd(dx1, abs(dy1));\n    dx1 /= a;\n    dy1 /= a;\n    a = gcd(dx2, abs(dy2));\n    dx2 /= a;\n    dy2 /= a;\n    if(dx1!=dx2 || dy1 != dy2) return 1;\n\n    a = crossProd(x2-x1, y2-y1, x3-x1, y3-y1);\n    if(a==0) return 2;\n    return 0;\n  }\n\n  if(p3 == 2){\n    swap(x1, x3); swap(y1, y3); swap(p1, p3);\n    swap(x2, x4); swap(y2, y4); swap(p2, p4);\n  }\n\n  a = crossProd(x2-x1, y2-y1, x3-x1, y3-y1);\n  b = crossProd(x2-x1, y2-y1, x4-x1, y4-y1);\n\n  if(a==b==0){\n    if(p1==2) return 2;\n    if(x1==x2==x3==x4) x1 = y1, x2 = y2, x3 = y3, x4 = y4;\n    if(x1 > x2) swap(x1, x2), swap(p1, p2);\n    if(x3 > x4) swap(x3, x4), swap(p3, p4);\n    if(x2 == x3 && p2 == p3 == 1) return 1;\n    if(x4 == x1 && p4 == p1 == 1) return 1;\n    if(x2 < x3 || x4 < x1) return 0;\n    return 2;\n  }\n\n  if(a > b) a = -a, b = -b;\n  if(a > 0 || (p3==0 && a >= 0)) return 0;\n  if(b < 0 || (p4==0 && b <= 0)) return 0;\n  if(p1==2) return 1;\n\n  a = crossProd(x4-x3, y4-y3, x1-x3, y1-y3);\n  b = crossProd(x4-x3, y4-y3, x2-x3, y2-y3);\n  if(a > b) a = -a, b = -b;\n  if(a > 0 || (p1==0 && a >= 0)) return 0;\n  if(b < 0 || (p2==0 && b <= 0)) return 0;\n\n  return 1;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"gcd");
      d.push_back((string)"crossProd");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "dimcomp2";
      string c = "struct dimcomp2 {\n  int B;\n  dimcomp2(){}\n  dimcomp2(int b){\n    B = b;\n  }\n  dimcomp2(int a, int b){\n    B = b;\n  }\n  inline void set(int b){\n    B = b;\n  }\n  inline void set(int a, int b){\n    B = b;\n  }\n  inline int mask(int a, int b){\n    return a * B + b;\n  }\n  inline int operator()(int a, int b){\n    return a * B + b;\n  }\n  inline void para(int mask, int &a, int &b){\n    a = mask / B;\n    b = mask % B;\n  }\n  inline void operator()(int mask, int &a, int &b){\n    a = mask / B;\n    b = mask % B;\n  }\n};\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "dimcomp3";
      string c = "struct dimcomp3 {\n  int B, C;\n  dimcomp3(){};\n  dimcomp3(int b, int c){\n    B = b; C = c;\n  }\n  dimcomp3(int a, int b, int c){\n    B = b; C = c;\n  }\n  inline void set(int b, int c){\n    B = b; C = c;\n  }\n  inline void set(int a, int b, int c){\n    B = b; C = c;\n  }\n  inline int mask(int a, int b, int c){\n    return (a * B + b) * C + c;\n  }\n  inline int operator()(int a, int b, int c){\n    return (a * B + b) * C + c;\n  }\n  inline void para(int mask, int &a, int &b, int &c){\n    a = mask / (B*C);\n    b = mask % (B*C) / C;\n    c = mask % C;\n  }\n  inline void operator()(int mask, int &a, int &b, int &c){\n    a = mask / (B*C);\n    b = mask % (B*C) / C;\n    c = mask % C;\n  }\n};\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "dimcomp4";
      string c = "struct dimcomp4 {\n  int B, C, D;\n  dimcomp4(){}\n  dimcomp4(int b, int c, int d){\n    B = b; C = c; D = d;\n  }\n  dimcomp4(int a, int b, int c, int d){\n    B = b; C = c; D = d;\n  }\n  inline void set(int b, int c, int d){\n    B = b; C = c; D = d;\n  }\n  inline void set(int a, int b, int c, int d){\n    B = b; C = c; D = d;\n  }\n  inline int operator()(int a, int b, int c, int d){\n    return ((a * B + b) * C + c) * D + d;\n  }\n  inline int mask(int a, int b, int c, int d){\n    return ((a * B + b) * C + c) * D + d;\n  }\n  inline void operator()(int mask, int &a, int &b, int &c, int &d){\n    a = mask / (B*C*D);\n    b = mask % (B*C*D) / (C*D);\n    c = mask % (C*D) / D;\n    d = mask % D;\n  }\n  inline void para(int mask, int &a, int &b, int &c, int &d){\n    a = mask / (B*C*D);\n    b = mask % (B*C*D) / (C*D);\n    c = mask % (C*D) / D;\n    d = mask % D;\n  }\n};\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "HammingDistance";
      string c = "template<class S, class T>\nint HammingDistance(int N, S A[], T B[]){\n  int i, res = 0;\n  rep(i,N) if(A[i] != B[i]) res++;\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "editDistance";
      string c = "template<class T>\nint editDistance(int As, T A[], int Bs, T B[], void *mem = wmem){\n  int i, j, k, *dp;\n\n  walloc1d(&dp, (As+1)*(Bs+1), &mem);\n\n  rep(i,As+1) dp[i*(Bs+1)] = i;\n  rep(i,Bs+1) dp[i] = i;\n  \n  REP(i,1,As+1) REP(j,1,Bs+1){\n    k = min(dp[i*(Bs+1)+j-1], dp[(i-1)*(Bs+1)+j]) + 1;\n    if(A[i-1]==B[j-1]) k = min(k, dp[(i-1)*(Bs+1)+(j-1)]);\n    else               k = min(k, dp[(i-1)*(Bs+1)+(j-1)]+1);\n    dp[i*(Bs+1)+j] = k;\n  }\n\n  return dp[As*(Bs+1)+Bs];\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"min_L");
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "counterSumIsLT";
      string c = "template<class T>\nll counterSumIsLT(int As, T A[], int Bs, T B[], T val){\n  int i = 0, j = Bs;\n  ll res = 0;\n\n  while(i < As){\n    while(j && A[i] + B[j-1] >= val) j--;\n    if(!j) break;\n    while(i<As && A[i] + B[j-1] < val) i++, res += j;\n  }\n\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "counterProdIsLT";
      string c = "template<class T>\nll counterProdIsLT(int As, T A[], int Bs, T B[], T val){\n  int i, j;\n  int ma = 0, za = 0, pa = 0;\n  int mb = 0, zb = 0, pb = 0;\n  ll res = 0;\n\n  i = 0;\n  while(i < As && A[i] < 0) i++;\n  ma = i;\n  while(i < As && A[i] == 0) i++;\n  za = i - ma;\n  pa = As - i;\n\n  i = 0;\n  while(i < Bs && B[i] < 0) i++;\n  mb = i;\n  while(i < Bs && B[i] == 0) i++;\n  zb = i - mb;\n  pb = Bs - i;\n\n  if(val < 0){\n    j = 0;\n    rep(i,pa){\n      while(j < mb && A[As-pa+i] * B[j] < val) j++;\n      res += j;\n    }\n    j = 0;\n    rep(i,pb){\n      while(j < ma && B[Bs-pb+i] * A[j] < val) j++;\n      res += j;\n    }\n  } else if(val == 0){\n    res = (ll) ma * pb + (ll) pa * mb;\n  } else {\n    res = (ll) As * Bs - (ll) pa * pb - (ll) ma * mb;\n    j = pb;\n    rep(i,pa){\n      while(j && A[As-pa+i] * B[Bs-pb+j-1] >= val) j--;\n      if(j==0) break;\n      res += j;\n    }\n    j = mb;\n    rep(i,ma){\n      while(j && A[ma-1-i] * B[mb-j] >= val) j--;\n      if(j==0) break;\n      res += j;\n    }\n  }\n\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "counterD2SumIsLT";
      string c = "template<class T>\nll counterD2SumIsLT(int As, T A[], T val){\n  int i = 0, j = As - 1;\n  ll res = 0;\n\n  while(i < j){\n    while(i < j && A[i] + A[j] >= val) j--;\n    if(i >= j) break;\n    while(i < j && A[i] + A[j] < val) res += j - i, i++;\n  }\n\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "counterM2SumIsLT";
      string c = "template<class T>\nll counterM2SumIsLT(int As, T A[], T val){\n  int i = 0, j = As - 1;\n  ll res = 0;\n\n  while(i <= j){\n    while(i <= j && A[i] + A[j] >= val) j--;\n    if(i > j) break;\n    while(i <= j && A[i] + A[j] < val) res += j - i + 1, i++;\n  }\n\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "counterD2ProdIsLT";
      string c = "template<class T>\nll counterD2ProdIsLT(int As, T A[], T val){\n  int i, j;\n  int ma = 0, za = 0, pa = 0;\n  ll res = 0;\n\n  i = 0;\n  while(i < As && A[i] < 0) i++;\n  ma = i;\n  while(i < As && A[i] == 0) i++;\n  za = i - ma;\n  pa = As - i;\n\n  if(val < 0){\n    j = 0;\n    rep(i,pa){\n      while(j < ma && A[As-pa+i] * A[j] < val) j++;\n      res += j;\n    }\n  } else if(val == 0){\n    res = (ll) ma * pa;\n  } else {\n    res = (ll) As * (As - 1) / 2 - (ll) pa * (pa - 1) / 2 - (ll) ma * (ma - 1) / 2;\n    j = pa;\n    rep(i,pa){\n      while(j - 1 > i && A[As-pa+i] * A[As-pa+j-1] >= val) j--;\n      if(j - 1 <= i) break;\n      res += j - i - 1;\n    }\n    j = ma;\n    rep(i,ma){\n      while(j - 1 > i && A[ma-1-i] * A[ma-j] >= val) j--;\n      if(j - 1 <= i) break;\n      res += j - i - 1;\n    }\n  }\n\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "counterM2ProdIsLT";
      string c = "template<class T>\nll counterM2ProdIsLT(int As, T A[], T val){\n  int i, j;\n  int ma = 0, za = 0, pa = 0;\n  ll res = 0;\n\n  i = 0;\n  while(i < As && A[i] < 0) i++;\n  ma = i;\n  while(i < As && A[i] == 0) i++;\n  za = i - ma;\n  pa = As - i;\n\n  if(val < 0){\n    j = 0;\n    rep(i,pa){\n      while(j < ma && A[As-pa+i] * A[j] < val) j++;\n      res += j;\n    }\n  } else if(val == 0){\n    res = (ll) ma * pa;\n  } else {\n    res = (ll) As * (As + 1) / 2 - (ll) pa * (pa + 1) / 2 - (ll) ma * (ma + 1) / 2;\n    j = pa;\n    rep(i,pa){\n      while(j - 1 >= i && A[As-pa+i] * A[As-pa+j-1] >= val) j--;\n      if(j - 1 < i) break;\n      res += j - i;\n    }\n    j = ma;\n    rep(i,ma){\n      while(j - 1 >= i && A[ma-1-i] * A[ma-j] >= val) j--;\n      if(j - 1 < i) break;\n      res += j - i;\n    }\n  }\n\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "rangeTree2d";
      string c = "template<class S, class T1, class T2>\nstruct rangeTree2d{\n  int N, N2;\n  int *sz;\n  S *tot;\n  fenwick<S> *w;\n  T1 **d1, *ddd1;\n  T2 *d2;\n  inline void build(int nn, T1 dd1[], T2 dd2[], S ww[] = NULL, void **mem = &wmem){\n    int i, j, i1, i2, k1, k2, s, s1, s2;\n    S *www;\n    int *ind;\n    N = nn;\n    for(N2=1;N2<N;N2*=2);\n    walloc1d(&sz,2*N2,mem);\n    walloc1d(&tot,2*N2,mem);\n    walloc1d(&w,2*N2,mem);\n    walloc1d(&d1,2*N2,mem);\n    walloc1d(&d2,nn,mem);\n    malloc1d(&www,nn);\n    walloc1d(&ddd1,nn);\n    walloc1d(&ind,nn);\n    rep(i,N) ddd1[i] = dd1[i];\n    rep(i,N) d2[i] = dd2[i];\n    if(ww==NULL){\n      rep(i,N) www[i] = 1;\n      sortA(N,d2,ddd1);\n    } else {\n      rep(i,N) ind[i] = i;\n      sortA(N,d2,ddd1,ind);\n      rep(i,N) www[i] = ww[ind[i]];\n    }\n    rep(i,N){\n      sz[N2+i] = 1;\n      walloc1d(&d1[N2+i], 1, mem);\n      d1[N2+i][0] = ddd1[i];\n      w[N2+i].walloc(1, mem);\n      w[N2+i].init(1);\n      w[N2+i].add(0,www[i]);\n      tot[N2+i] = www[i];\n    }\n    rep(i,N,N2){\n      sz[N2+i] = 0;\n      tot[N2+i] = 0;\n    }\n    rrep(i,1,N2){\n      i1 = 2i;\n      i2 = 2i + 1;\n      s1 = sz[i1];\n      s2 = sz[i2];\n      sz[i] = s1 + s2;\n      s = k1 = k2 = 0;\n      walloc1d(&d1[i], sz[i], mem);\n      w[i].walloc(sz[i], mem);\n      w[i].init(sz[i]);\n      while(k1 < s1 || k2 < s2){\n        if(k2==s2){\n          d1[i][s] = d1[i1][k1];\n          w[i].add(s,w[i1].range(k1,k1));\n          s++; k1++; continue;\n        }\n        if(k1==s1){\n          d1[i][s] = d1[i2][k2];\n          w[i].add(s,w[i2].range(k2,k2));\n          s++; k2++; continue;\n        }\n        if(d1[i1][k1] < d1[i2][k2]){\n          d1[i][s] = d1[i1][k1];\n          w[i].add(s,w[i1].range(k1,k1));\n          s++; k1++; continue;\n        }else{\n          d1[i][s] = d1[i2][k2];\n          w[i].add(s,w[i2].range(k2,k2));\n          s++; k2++; continue;\n        }\n      }\n    }\n    free1d(www);\n  }\n  inline void add(T1 x, T2 y, S v){\n    int a, b, z;\n    a = lower_bound(d2, d2+N, y) - d2;\n    b = upper_bound(d2, d2+N, y) - d2;\n    z = lower_bound(ddd1+a, ddd1+b, x) - ddd1 + N2;\n    while(z){\n      a = lower_bound(d1[z], d1[z]+sz[z], x) - d1[z];\n      w[z].add(a, v);\n      z /= 2;\n    }\n  }\n  inline S query(T1 x1, T1 x2, T2 y1, T2 y2){\n    S res = 0;\n    int a, b, z1, z2;\n    a = lower_bound(d2, d2+N, y1) - d2 + N2;\n    b = lower_bound(d2, d2+N, y2) - d2 + N2;\n    while(a < b){\n      if(a%2){\n        z1 = lower_bound(d1[a], d1[a]+sz[a], x1) - d1[a];\n        z2 = lower_bound(d1[a], d1[a]+sz[a], x2) - d1[a];\n        if(z1 < z2) res += w[a].range(z1,z2-1);\n        a++;\n      }\n      if(b%2){\n        b--;\n        z1 = lower_bound(d1[b], d1[b]+sz[b], x1) - d1[b];\n        z2 = lower_bound(d1[b], d1[b]+sz[b], x2) - d1[b];\n        if(z1 < z2) res += w[b].range(z1,z2-1);\n      }\n      a /= 2;\n      b /= 2;\n    }\n    return res;\n  }\n};\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"fenwick");
      d.push_back((string)"walloc1d");
      d.push_back((string)"malloc1d");
      d.push_back((string)"free1d");
      d.push_back((string)"sortA");
      d.push_back((string)"workmemory");
      need[n] = d;
    }
    {
      string n = "rangeTree2d_nw";
      string c = "template<class T1, class T2>\nstruct rangeTree2d_nw{\n  int N, N2;\n  int *sz;\n  T1 **d1;\n  T2 *d2;\n  inline void build(int nn, T1 dd1[], T2 dd2[], void **mem = &wmem){\n    int i, j, i1, i2, k1, k2, s, s1, s2;\n    T1 *ddd1;\n    N = nn;\n    for(N2=1;N2<N;N2*=2);\n    walloc1d(&sz,2*N2,mem);\n    walloc1d(&d1,2*N2,mem);\n    walloc1d(&d2,nn,mem);\n    malloc1d(&ddd1,nn);\n    rep(i,N) ddd1[i] = dd1[i];\n    rep(i,N) d2[i] = dd2[i];\n    sortA(N,d2,ddd1);\n    rep(i,N){\n      sz[N2+i] = 1;\n      walloc1d(&d1[N2+i], 1, mem);\n      d1[N2+i][0] = ddd1[i];\n    }\n    rep(i,N,N2){\n      sz[N2+i] = 0;\n    }\n    rrep(i,1,N2){\n      i1 = 2i;\n      i2 = 2i + 1;\n      s1 = sz[i1];\n      s2 = sz[i2];\n      sz[i] = s1 + s2;\n      s = k1 = k2 = 0;\n      walloc1d(&d1[i], sz[i], mem);\n      while(k1 < s1 || k2 < s2){\n        if(k2==s2){\n          d1[i][s] = d1[i1][k1];\n          s++; k1++; continue;\n        }\n        if(k1==s1){\n          d1[i][s] = d1[i2][k2];\n          s++; k2++; continue;\n        }\n        if(d1[i1][k1] < d1[i2][k2]){\n          d1[i][s] = d1[i1][k1];\n          s++; k1++; continue;\n        }else{\n          d1[i][s] = d1[i2][k2];\n          s++; k2++; continue;\n        }\n      }\n    }\n    free1d(ddd1);\n  }\n  inline int query(T1 x1, T1 x2, T2 y1, T2 y2){\n    int res = 0;\n    int a, b, z1, z2;\n    a = lower_bound(d2, d2+N, y1) - d2 + N2;\n    b = lower_bound(d2, d2+N, y2) - d2 + N2;\n    while(a < b){\n      if(a%2){\n        z1 = lower_bound(d1[a], d1[a]+sz[a], x1) - d1[a];\n        z2 = lower_bound(d1[a], d1[a]+sz[a], x2) - d1[a];\n        if(z1 < z2) res += z2-z1;\n        a++;\n      }\n      if(b%2){\n        b--;\n        z1 = lower_bound(d1[b], d1[b]+sz[b], x1) - d1[b];\n        z2 = lower_bound(d1[b], d1[b]+sz[b], x2) - d1[b];\n        if(z1 < z2) res += z2-z1;\n      }\n      a /= 2;\n      b /= 2;\n    }\n    return res;\n  }\n};\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"malloc1d");
      d.push_back((string)"free1d");
      d.push_back((string)"sortA");
      d.push_back((string)"workmemory");
      need[n] = d;
    }
    {
      string n = "rangeTree2d_pf_header";
      string c = "template<class T>\nstruct segtree_ph_inRT{\n  int N, logN;\n  T *val;\n  void malloc(int maxN, int once = 0);\n  void walloc(int maxN, void **mem = &wmem);\n  void walloc(int maxN, int once = 0, void **mem = &wmem);\n  void free(void);\n  T& operator[](int i);\n  void setN(int n, int zerofill = 1, int dobuild = 1);\n  void build(void);\n  inline void build(int a);\n  inline void change(int a, T v);\n  inline void add(int a, T v);\n  inline T get(int a, int b);\n};\ntemplate<class S, class T1, class T2>\nstruct rangeTree2d_pf{\n  int N, N2;\n  int *sz;\n  S *tot, defval;\n  segtree_ph_inRT<S> *w;\n  T1 **d1, *ddd1;\n  T2 *d2;\n  inline void build(int nn, T1 dd1[], T2 dd2[], S ww[] = NULL, void **mem = &wmem);\n  inline void change(T1 x, T2 y, S v);\n  inline void add(T1 x, T2 y, S v);\n  inline void setDefault(const S val){\n    defval = val;\n  }\n  inline S query(T1 x1, T1 x2, T2 y1, T2 y2);\n};\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "rangeTree2d_pf";
      string c = "template<class T> void segtree_ph_inRT<T>::malloc(int maxN, int once /*= 0*/){\n  int i;\n  for(i=1;i<maxN;i*=2);\n  val = new T[2*i];\n  if(once) setN(maxN);\n}\ntemplate<class T> void segtree_ph_inRT<T>::walloc(int maxN, void **mem /*= &wmem*/){\n  int i;\n  for(i=1;i<maxN;i*=2);\n  walloc1d(&val, 2i, mem);\n}\ntemplate<class T> void segtree_ph_inRT<T>::walloc(int maxN, int once /*= 0*/, void **mem /*= &wmem*/){\n  int i;\n  for(i=1;i<maxN;i*=2);\n  walloc1d(&val, 2i, mem);\n  if(once) setN(maxN);\n}\ntemplate<class T> void segtree_ph_inRT<T>::free(void){\n  delete [] val;\n}\ntemplate<class T> T& segtree_ph_inRT<T>::operator[](int i){\n  return val[N+i];\n}\ntemplate<class T> void segtree_ph_inRT<T>::setN(int n, int zerofill /*= 1*/, int dobuild /*= 1*/){\n  int i;\n  for(i=1,logN=0;i<n;i*=2,logN++);\n  N = i;\n  if(dobuild) build();\n}\ntemplate<class T> void segtree_ph_inRT<T>::build(void){\n  for(int i=N-1;i;i--) val[i] = rangeTree2d_pf_func(val[2i], val[2i+1]);\n}\ntemplate<class T> inline void segtree_ph_inRT<T>::build(int a){\n  while(a > 1){\n    a /= 2;\n    val[a] = rangeTree2d_pf_func(val[2a], val[2a+1]);\n  }\n}\ntemplate<class T> inline void segtree_ph_inRT<T>::change(int a, T v){\n  val[a+N] = v;\n  build(a+N);\n}\ntemplate<class T> inline void segtree_ph_inRT<T>::add(int a, T v){\n  val[a+N] += v;\n  build(a+N);\n}\ntemplate<class T> inline T segtree_ph_inRT<T>::get(int a, int b){\n  T res, tmp;\n  int fga = 0, fgb = 0;\n  a += N;\n  b += N;\n  while(a < b){\n    if(a%2){\n      if(fga){\n        res = rangeTree2d_pf_func(res, val[a]);\n      } else {\n        res = val[a];\n        fga = 1;\n      }\n      a++;\n    }\n    if(b%2){\n      b--;\n      if(fgb){\n        tmp = rangeTree2d_pf_func(val[b], tmp);\n      } else {\n        tmp = val[b];\n        fgb = 1;\n      }\n    }\n    a /= 2;\n    b /= 2;\n  }\n  if(fga==1 && fgb==0) return res;\n  if(fga==0 && fgb==1) return tmp;\n  if(fga==1 && fgb==1) return rangeTree2d_pf_func(res, tmp);\n  return res;\n}\ntemplate<class S, class T1, class T2> inline void rangeTree2d_pf<S,T1,T2>::build(int nn, T1 dd1[], T2 dd2[], S ww[]/* = NULL*/, void **mem /* = &wmem*/){\n  int i, j, i1, i2, k1, k2, s, s1, s2;\n  S *www;\n  int *ind;\n  N = nn;\n  for(N2=1;N2<N;N2*=2);\n  walloc1d(&sz,2*N2,mem);\n  walloc1d(&tot,2*N2,mem);\n  walloc1d(&w,2*N2,mem);\n  walloc1d(&d1,2*N2,mem);\n  walloc1d(&d2,nn,mem);\n  malloc1d(&www,nn);\n  walloc1d(&ddd1,nn);\n  walloc1d(&ind,nn);\n  rep(i,N) ddd1[i] = dd1[i];\n  rep(i,N) d2[i] = dd2[i];\n  if(ww==NULL){\n    sortA(N,d2,ddd1);\n  } else {\n    rep(i,N) ind[i] = i;\n    sortA(N,d2,ddd1,ind);\n    rep(i,N) www[i] = ww[ind[i]];\n  }\n  rep(i,N){\n    sz[N2+i] = 1;\n    walloc1d(&d1[N2+i], 1, mem);\n    d1[N2+i][0] = ddd1[i];\n    w[N2+i].walloc(1, mem);\n    w[N2+i].setN(1);\n    w[N2+i][0] = www[i];\n    w[N2+i].build();\n    tot[N2+i] = www[i];\n  }\n  rep(i,N,N2){\n    sz[N2+i] = 0;\n    tot[N2+i] = 0;\n  }\n  rrep(i,1,N2){\n    i1 = 2i;\n    i2 = 2i + 1;\n    s1 = sz[i1];\n    s2 = sz[i2];\n    sz[i] = s1 + s2;\n    s = k1 = k2 = 0;\n    walloc1d(&d1[i], sz[i], mem);\n    w[i].walloc(sz[i], mem);\n    w[i].setN(sz[i]);\n    while(k1 < s1 || k2 < s2){\n      if(k2==s2){\n        d1[i][s] = d1[i1][k1];\n        w[i][s] = w[i1].get(k1,k1+1);\n        s++; k1++; continue;\n      }\n      if(k1==s1){\n        d1[i][s] = d1[i2][k2];\n        w[i][s] = w[i2].get(k2,k2+1);\n        s++; k2++; continue;\n      }\n      if(d1[i1][k1] < d1[i2][k2]){\n        d1[i][s] = d1[i1][k1];\n        w[i][s] = w[i1].get(k1,k1+1);\n        s++; k1++; continue;\n      }else{\n        d1[i][s] = d1[i2][k2];\n        w[i][s] = w[i2].get(k2,k2+1);\n        s++; k2++; continue;\n      }\n    }\n    w[i].build();\n  }\n  free1d(www);\n}\ntemplate<class S, class T1, class T2> inline void rangeTree2d_pf<S,T1,T2>::change(T1 x, T2 y, S v){\n  int a, b, z;\n  a = lower_bound(d2, d2+N, y) - d2;\n  b = upper_bound(d2, d2+N, y) - d2;\n  z = lower_bound(ddd1+a, ddd1+b, x) - ddd1 + N2;\n  while(z){\n    a = lower_bound(d1[z], d1[z]+sz[z], x) - d1[z];\n    w[z].change(a, v);\n    z /= 2;\n  }\n}\ntemplate<class S, class T1, class T2> inline void rangeTree2d_pf<S,T1,T2>::add(T1 x, T2 y, S v){\n  int a, b, z;\n  a = lower_bound(d2, d2+N, y) - d2;\n  b = upper_bound(d2, d2+N, y) - d2;\n  z = lower_bound(ddd1+a, ddd1+b, x) - ddd1 + N2;\n  while(z){\n    a = lower_bound(d1[z], d1[z]+sz[z], x) - d1[z];\n    w[z].add(a, v);\n    z /= 2;\n  }\n}\ntemplate<class S, class T1, class T2> inline S rangeTree2d_pf<S,T1,T2>::query(T1 x1, T1 x2, T2 y1, T2 y2){\n  S res;\n  int a, b, z1, z2, fg = 0;\n  a = lower_bound(d2, d2+N, y1) - d2 + N2;\n  b = lower_bound(d2, d2+N, y2) - d2 + N2;\n  while(a < b){\n    if(a%2){\n      z1 = lower_bound(d1[a], d1[a]+sz[a], x1) - d1[a];\n      z2 = lower_bound(d1[a], d1[a]+sz[a], x2) - d1[a];\n      if(z1 < z2){\n        if(fg == 0){\n          fg = 1;\n          res = w[a].get(z1,z2);\n        } else {\n          res = rangeTree2d_pf_func(res, w[a].get(z1,z2));\n        }\n      }\n      a++;\n    }\n    if(b%2){\n      b--;\n      z1 = lower_bound(d1[b], d1[b]+sz[b], x1) - d1[b];\n      z2 = lower_bound(d1[b], d1[b]+sz[b], x2) - d1[b];\n      if(z1 < z2){\n        if(fg == 0){\n          fg = 1;\n          res = w[b].get(z1,z2);\n        } else {\n          res = rangeTree2d_pf_func(res, w[b].get(z1,z2));\n        }\n      }\n    }\n    a /= 2;\n    b /= 2;\n  }\n  if(!fg) return defval;\n  return res;\n}\n";
      string p = "last";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"rangeTree2d_pf_header");
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"malloc1d");
      d.push_back((string)"free1d");
      d.push_back((string)"sortA");
      need[n] = d;
    }


    {
      string n = "Permutation";
      string c = "struct Permutation {\n  int n, mem;\n  int *dat;\n\n  Permutation(){n = mem = 0;}\n  Permutation(const int nn){\n    n = mem = nn;\n    if(mem > 0) dat = new int[mem];\n  }\n  Permutation(const Permutation &a){\n    int i;\n    mem = n = a.n;\n    dat = new int[mem];\n    rep(i,mem) dat[i] = a.dat[i];\n  }\n  \n  ~Permutation(){\n    if(mem) delete [] dat;\n  }\n\n  void changeSize(const int nn){\n    n = nn;\n    if(mem < n){\n      if(mem) delete [] dat;\n      mem = n;\n      dat = new int[mem];\n    }\n  }\n\n  Permutation& operator=(const Permutation &a){\n    int i;\n    changeSize(a.n);\n    n = a.n;\n    rep(i,n) dat[i] = a.dat[i];\n    return *this;\n  }\n\n  Permutation& operator=(const int a){\n    int i;\n    rep(i,n) dat[i] = i;\n    return *this;\n  }\n\n  Permutation& operator*=(const Permutation &a){\n    int i, *m;\n    void *mv = wmem;\n    \n    if(n==0 || n!=a.n){\n      changeSize(0);\n      return *this;\n    }\n    walloc1d(&m, n, &mv);\n    rep(i,n) m[i] = dat[a.dat[i]];\n    rep(i,n) dat[i] = m[i];\n    return *this;\n  }\n  Permutation operator*(const Permutation &a){\n    return Permutation(*this) *= a;\n  }\n\n  bool operator==(const Permutation &a){\n    int i;\n    if(n != a.n) return false;\n    rep(i,n) if(dat[i] != a.dat[i]) return false;\n    return true;\n  }\n\n  template<class T>\n  void apply(T A[]){\n    int i;\n    T *B; void *mv = wmem;\n    walloc1d(&B, n, &mv);\n    rep(i,n) B[dat[i]] = A[i];\n    rep(i,n) A[i] = B[i];\n  }\n\n  template<class T>\n  void apply(T A[], T B[]){\n    int i;\n    rep(i,n) B[dat[i]] = A[i];\n  }\n\n  int cycle_len(int res[] = NULL){\n    int i, j, k, sz = 0;\n    int *vis;\n    void *mv = wmem;\n\n    if(res==NULL) walloc1d(&res, n, &mv);\n    walloc1d(&vis, n, &mv);\n\n    rep(i,n) vis[i] = 0;\n    rep(i,n) if(!vis[i]){\n      k = 0;\n      j = i;\n      while(vis[j]==0){\n        vis[j] = 1;\n        j = dat[j];\n        k++;\n      }\n      res[sz++] = k;\n    }\n    return sz;\n  }\n\n  void cycle_len_EachElement(int res[]){\n    int i, j, k, sz = 0;\n    int *vis;\n    void *mv = wmem;\n\n    walloc1d(&vis, n, &mv);\n\n    rep(i,n) vis[i] = 0;\n    rep(i,n) if(!vis[i]){\n      k = 0;\n      j = i;\n      while(vis[j]==0){\n        vis[j] = 1;\n        j = dat[j];\n        k++;\n      }\n      j = i;\n      while(vis[j]==1){\n        res[j] = k;\n        vis[j] = 2;\n        j = dat[j];\n      }\n    }\n  }\n\n  template<class T>\n  inline T getIndex(void *mem = wmem){\n    int i;\n    fenwick<int> t;\n    T res, *fac;\n\n    walloc1d(&fac,n,&mem);\n    fac[0] = 1;\n    rep(i,1,n) fac[i] = i * fac[i-1];\n\n    t.walloc(n,&mem);\n    t.init(n);\n    rep(i,n) t.add(i,1);\n\n    res = 0;\n    rep(i,n){\n      t.add(dat[i], -1);\n      res += fac[n-1-i] * t.get(dat[i]-1);\n    }\n    return res;\n  }\n\n  inline int& operator[](const int a){\n    return dat[a];\n  }\n};\n\ntemplate<class S> inline Permutation pow_L(Permutation a, S b){\n  Permutation res;\n  res.changeSize(a.n);\n  res = 1;\n  while(b){\n    if(b&1) res *= a;\n    b >>= 1;\n    a *= a;\n  }\n  return res;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"fenwick");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "IntMap";
      string c = "struct IntMap {\n  int n, mem, logn;\n  int *dat;\n\n  int **nx, calc_nx;\n  int numCycle, *cycleLen, **cycle;\n  int *int2c, *int2cInd;\n\n  IntMap(){\n    n = mem = logn = 0;\n    calc_nx = 0;\n    numCycle = 0;\n  }\n\n  IntMap(int nn){\n    n = mem = nn;\n    for(logn=0; (1<<logn)<n; logn++);\n    logn++;\n    if(mem > 0) dat = new int[mem];\n    calc_nx = 0;\n    numCycle = 0;\n  }\n\n  ~IntMap(){\n    if(mem) delete[] dat;\n    if(calc_nx){\n      rep(i,logn) delete[] nx[i];\n      delete[] nx;\n    }\n    if(numCycle > 0){\n      delete[] cycleLen;\n      rep(i,numCycle) delete[] cycle[i];\n      delete[] cycle;\n      delete[] int2c;\n      delete[] int2cInd;\n    }\n  }\n\n  void changeSize(const int nn){\n    int old_logn = logn;\n    n = nn;\n    for(logn=0; (1<<logn)<n; logn++);\n    logn++;\n\n    if(mem < n){\n      if(mem) delete [] dat;\n      mem = n;\n      dat = new int [mem];\n    }\n\n    if(calc_nx){\n      rep(i,old_logn) delete[] nx[i];\n      delete[] nx;\n      calc_nx = 0;\n    }\n\n    if(numCycle > 0){\n      delete[] cycleLen;\n      rep(i,numCycle) delete[] cycle[i];\n      delete[] cycle;\n      delete[] int2c;\n      delete[] int2cInd;\n      numCycle = 0;\n    }\n  }\n\n  int calcCycle(void){\n    int i, j, k;\n\n    if(numCycle){\n      delete[] cycleLen;\n      rep(i,numCycle) delete[] cycle[i];\n      delete[] cycle;\n      delete[] int2c;\n      delete[] int2cInd;\n      numCycle = 0;\n    }\n\n    int2c = new int[n];\n    int2cInd = new int[n];\n    rep(i,n) int2c[i] = -2;\n    rep(i,n) int2cInd[i] = -1;\n\n    numCycle = 0;\n    rep(i,n) if(int2c[i] == -2){\n      j = i;\n      for(;;){\n        if(int2c[j] != -2) break;\n        int2c[j] = -3;\n        j = dat[j];\n      }\n      if(int2c[j] == -3){\n        k = 0;\n        for(;;){\n          if(int2c[j] != -3) break;\n          int2c[j] = numCycle;\n          int2cInd[j] = k++;\n          j = dat[j];\n        }\n        numCycle++;\n      }\n      j = i;\n      for(;;){\n        if(int2c[j] != -3) break;\n        int2c[j] = -1;\n        j = dat[j];\n      }\n    }\n\n    cycleLen = new int[numCycle];\n    rep(i,numCycle) cycleLen[i] = 0;\n    rep(i,n) if(int2c[i] >= 0) cycleLen[int2c[i]]++;\n    cycle = new int*[numCycle];\n    rep(i,numCycle) cycle[i] = new int[cycleLen[i]];\n    rep(i,n) if(int2c[i] >= 0) cycle[int2c[i]][int2cInd[i]] = i;\n\n    return numCycle;\n  }\n\n  void calcNext(int recalc = 1){\n    if(recalc || numCycle==0) calcCycle();\n    if(calc_nx == 0){\n      calc_nx = 1;\n      nx = new int*[logn];\n      rep(i,logn) nx[i] = new int[n];\n    }\n    rep(i,n) nx[0][i] = dat[i];\n    rep(k,1,logn) rep(i,n) nx[k][i] = nx[k-1][nx[k-1][i]];\n  }\n\n  template<class T>\n  int getNext(int x, T s){\n    if(calc_nx==0) calcNext();\n    if(s >= (1<<(logn-1))){\n      x = nx[logn-1][x];\n      s -= (1<<(logn-1));\n      return cycle[int2c[x]][(int2cInd[x]+s)%cycleLen[int2c[x]]];\n    }\n    rep(i,logn) if(s&1<<i) x = nx[i][x];\n    return x;\n  }\n\n  inline int& operator[](const int a){\n    return dat[a];\n  }\n};\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "TSP_cycle";
      string c = "template<class T>\nT TSP_cycle(int n, T **dist, void *mem = wmem){\n  const T m = numeric_limits<T>::max();\n  int i, j, ii, jj, mask;\n  int as, bs, *a, *b;\n  T **dp, res = m;\n  if(n==1) return 0;\n  n--;\n  walloc2d(&dp, n, 1<<n, &mem);\n  walloc1d(&a, n, &mem);\n  walloc1d(&b, n, &mem);\n  rep(i,n) rep(j,1<<n) dp[i][j] = m;\n  rep(i,n) dp[i][1<<i] = dist[n][i];\n  rep(mask,1<<n){\n    as = bs = 0;\n    rep(i,n) if(dp[i][mask] < m) a[as++] = i;\n    rep(i,n) if(!(mask&1<<i)) b[bs++] = i;\n    rep(ii,as){\n      i = a[ii];\n      rep(jj,bs){\n        j = b[jj];\n        dp[j][mask|(1<<j)] <?= dp[i][mask] + dist[i][j];\n      }\n    }\n  }\n  rep(i,n) if(dp[i][(1<<n)-1] < m) res <?= dp[i][(1<<n)-1] + dist[i][n];\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"walloc2d");
      d.push_back((string)"workmemory");
      d.push_back((string)"chmin");
      need[n] = d;
    }
    {
      string n = "TSP_path";
      string c = "template<class T>\nT TSP_path(int n, T **dist, void *mem = wmem){\n  const T m = numeric_limits<T>::max();\n  int i, j, ii, jj, mask;\n  int as, bs, *a, *b;\n  T **dp, res = m;\n  if(n==1) return 0;\n  walloc2d(&dp, n, 1<<n, &mem);\n  walloc1d(&a, n, &mem);\n  walloc1d(&b, n, &mem);\n  rep(i,n) rep(j,1<<n) dp[i][j] = m;\n  rep(i,n) dp[i][1<<i] = 0;\n  rep(mask,1<<n){\n    as = bs = 0;\n    rep(i,n) if(dp[i][mask] < m) a[as++] = i;\n    rep(i,n) if(!(mask&1<<i)) b[bs++] = i;\n    rep(ii,as){\n      i = a[ii];\n      rep(jj,bs){\n        j = b[jj];\n        dp[j][mask|(1<<j)] <?= dp[i][mask] + dist[i][j];\n      }\n    }\n  }\n  rep(i,n) if(dp[i][(1<<n)-1] < m) res <?= dp[i][(1<<n)-1];\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"walloc2d");
      d.push_back((string)"workmemory");
      d.push_back((string)"chmin");
      need[n] = d;
    }
    {
      string n = "TSP_path_s";
      string c = "template<class T>\nT TSP_path_s(int n, T **dist, int s = 0, void *mem = wmem){\n  const T m = numeric_limits<T>::max();\n  int i, j, ii, jj, mask;\n  int as, bs, *a, *b;\n  T **dp, **d, res = m;\n  if(n==1) return 0;\n  walloc2d(&d, n, n, &mem);\n  rep(i,n) rep(j,n){\n    if(i < s) ii = i;\n    if(j < s) jj = j;\n    if(i == s) ii = n - 1;\n    if(j == s) jj = n - 1;\n    if(i > s) ii = i - 1;\n    if(j > s) jj = j - 1;\n    d[ii][jj] = dist[i][j];\n  }\n  n--;\n  walloc2d(&dp, n, 1<<n, &mem);\n  walloc1d(&a, n, &mem);\n  walloc1d(&b, n, &mem);\n  rep(i,n) rep(j,1<<n) dp[i][j] = m;\n  rep(i,n) dp[i][1<<i] = d[n][i];\n  rep(mask,1<<n){\n    as = bs = 0;\n    rep(i,n) if(dp[i][mask] < m) a[as++] = i;\n    rep(i,n) if(!(mask&1<<i)) b[bs++] = i;\n    rep(ii,as){\n      i = a[ii];\n      rep(jj,bs){\n        j = b[jj];\n        dp[j][mask|(1<<j)] <?= dp[i][mask] + d[i][j];\n      }\n    }\n  }\n  rep(i,n) if(dp[i][(1<<n)-1] < m) res <?= dp[i][(1<<n)-1];\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"walloc2d");
      d.push_back((string)"workmemory");
      d.push_back((string)"chmin");
      need[n] = d;
    }


    {
      string n = "rollingHash_init";
      string c = "{\n  rollingHashInit();\n}\n";
      string p = "main_first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "rollingHash";
      string c = "#define ROLLING_HASH_MOD (2305843009213693951ULL)\n#define ROLLING_HASH_PRIMITIVE_ROOT (3)\n#define ROLLING_HASH_MAX_MEMORY (2000000)\n\nint ROLLING_HASH_MEM;\null ROLLING_HASH_BASE, ROLLING_HASH_IBASE, *ROLLING_HASH_PW = NULL, *ROLLING_HASH_IPW = NULL;\n\ninline ull rollingHash61_mul(ull a, ull b){\n  __uint128_t r = (__uint128_t) a * b;\n  a = (r >> 61) + (r & ROLLING_HASH_MOD);\n  if(a >= ROLLING_HASH_MOD) a -= ROLLING_HASH_MOD;\n  return a;\n}\n\ninline ull rollingHash61_pow(ull a, ull b){\n  ull r = 1;\n  for(;;){\n    if(b&1){\n      r = rollingHash61_mul(r, a);\n    }\n    if(b==0) break;\n    b >>= 1;\n    a = rollingHash61_mul(a, a);\n  }\n  return r;\n}\n\nvoid rollingHashInit(){\n  int i;\n  Rand rnd;\n  ull x;\n  rep(i,20) rnd.get(2);\n  do{\n    x = rnd.get(1.0, (double)(ROLLING_HASH_MOD-2));\n  }while(gcd(x,ROLLING_HASH_MOD-1) != 1);\n  ROLLING_HASH_BASE = rollingHash61_pow(ROLLING_HASH_PRIMITIVE_ROOT, x);\n  ROLLING_HASH_IBASE = rollingHash61_pow(ROLLING_HASH_BASE, ROLLING_HASH_MOD - 2);\n}\n\nvoid rollingHash_expand(int k){\n  int i;\n  if(ROLLING_HASH_MEM >= k) return;\n  ROLLING_HASH_MEM = max(2 * ROLLING_HASH_MEM, k);\n  assert(ROLLING_HASH_MEM <= 2 * ROLLING_HASH_MAX_MEMORY);\n  ROLLING_HASH_PW = (ull*) realloc(ROLLING_HASH_PW, ROLLING_HASH_MEM * sizeof(ull));\n  ROLLING_HASH_IPW = (ull*) realloc(ROLLING_HASH_IPW, ROLLING_HASH_MEM * sizeof(ull));\n  ROLLING_HASH_PW[0] = 1;\n  rep(i,1,ROLLING_HASH_MEM) ROLLING_HASH_PW[i] = rollingHash61_mul(ROLLING_HASH_PW[i-1], ROLLING_HASH_BASE);\n  ROLLING_HASH_IPW[0] = 1;\n  rep(i,1,ROLLING_HASH_MEM) ROLLING_HASH_IPW[i] = rollingHash61_mul(ROLLING_HASH_IPW[i-1], ROLLING_HASH_IBASE);\n}\n\n\nstruct rollingHash{\n  ll len;\n  ull hs;\n\n  template<class T>\n  void set(int N, T A[]){\n    int i;\n    ll tmp;\n    hs = 0;\n    len = N;\n    rollingHash_expand(N);\n    rep(i,N){\n      tmp = A[i] % ((ll)ROLLING_HASH_MOD);\n      if(tmp < 0) tmp += ROLLING_HASH_MOD;\n      hs += rollingHash61_mul(tmp, ROLLING_HASH_PW[i]);\n      if(hs >= ROLLING_HASH_MOD) hs -= ROLLING_HASH_MOD;\n    }\n  }\n\n  template<class S, class T>\n  void change(ll ind, S bef, T aft){\n    ll tmp1, tmp2;\n    tmp1 = bef % ((ll)ROLLING_HASH_MOD);\n    tmp2 = aft % ((ll)ROLLING_HASH_MOD);\n    tmp1 = tmp2 - tmp1;\n    if(tmp1 < 0) tmp1 += ROLLING_HASH_MOD;\n    if(tmp1 < 0) tmp1 += ROLLING_HASH_MOD;\n    if(tmp1 >= ROLLING_HASH_MOD) tmp1 -= ROLLING_HASH_MOD;\n    if(ind+1 <= ROLLING_HASH_MAX_MEMORY || ind+1 >= ROLLING_HASH_MEM){\n      rollingHash_expand(ind+1);\n      hs += rollingHash61_mul(tmp1, ROLLING_HASH_PW[ind]);\n    } else {\n      hs += rollingHash61_mul(tmp1, rollingHash61_pow(ROLLING_HASH_BASE, ind));\n    }\n    if(hs >= ROLLING_HASH_MOD) hs -= ROLLING_HASH_MOD;\n  }\n\n  void push_front(rollingHash a){\n    if(a.len + 1 <= ROLLING_HASH_MAX_MEMORY || a.len + 1 >= ROLLING_HASH_MEM){\n      rollingHash_expand(a.len + 1);\n      hs = rollingHash61_mul(hs, ROLLING_HASH_PW[a.len]);\n    } else {\n      hs = rollingHash61_mul(hs, rollingHash61_pow(ROLLING_HASH_BASE, a.len));\n    }\n    hs += a.hs;\n    if(hs >= ROLLING_HASH_MOD) hs -= ROLLING_HASH_MOD;\n    len += a.len;\n  }\n\n  void push_back(rollingHash a){\n    if(len + 1 <= ROLLING_HASH_MAX_MEMORY || len + 1 >= ROLLING_HASH_MEM){\n      rollingHash_expand(len + 1);\n      hs += rollingHash61_mul(a.hs, ROLLING_HASH_PW[len]);\n    } else {\n      hs += rollingHash61_mul(a.hs, rollingHash61_pow(ROLLING_HASH_BASE, len));\n    }\n    if(hs >= ROLLING_HASH_MOD) hs -= ROLLING_HASH_MOD;\n    len += a.len;\n  }\n\n  void pop_front(rollingHash a){\n    if(hs >= a.hs){\n      hs -= a.hs;\n    } else {\n      hs = hs + ROLLING_HASH_MOD - a.hs;\n    }\n    if(a.len + 1 <= ROLLING_HASH_MAX_MEMORY || a.len + 1 >= ROLLING_HASH_MEM){\n      rollingHash_expand(a.len + 1);\n      hs = rollingHash61_mul(hs, ROLLING_HASH_IPW[a.len]);\n    } else {\n      hs = rollingHash61_mul(hs, rollingHash61_pow(ROLLING_HASH_IBASE, a.len));\n    }\n    len -= a.len;\n  }\n\n  void pop_back(rollingHash a){\n    ull tmp;\n    if(len + 1 <= ROLLING_HASH_MAX_MEMORY || len + 1 >= ROLLING_HASH_MEM){\n      rollingHash_expand(len + 1);\n      tmp = rollingHash61_mul(a.hs, ROLLING_HASH_PW[len]);\n    } else {\n      tmp = rollingHash61_mul(a.hs, rollingHash61_pow(ROLLING_HASH_BASE, len));\n    }\n    if(hs >= tmp){\n      hs -= tmp;\n    } else {\n      hs = hs + ROLLING_HASH_MOD - tmp;\n    }\n    len -= a.len;\n  }\n\n  bool operator==(const rollingHash a){\n    return len == a.len && hs == a.hs;\n  }\n  bool operator!=(const rollingHash a){\n    return len != a.len || hs != a.hs;\n  }\n};\n\n\ntemplate<class T>\nrollingHash calcRollingHash(int N, T A[]){\n  rollingHash res;\n  res.set(N, A);\n  return res;\n}\n\n\nstruct rollingHashSubarrays{\n  ull *hs;\n  int mem, len;\n\n  void set(){\n    hs = NULL;\n    mem = len = 0;\n  }\n  void free(){\n    if(mem) delete[] hs;\n  }\n  void expand(int k){\n    if(mem >= k) return;\n    free();\n    mem = max(2*mem, k);\n    hs = new ull[mem];\n  }\n\n  template<class T>\n  void set(int N, T A[]){\n    int i;\n    ll tmp;\n    if(N <= 0) return;\n    rollingHash_expand(N);\n    expand(N);\n    len = N;\n\n    tmp = A[0] % ((ll)ROLLING_HASH_MOD);\n    if(tmp < 0) tmp += ROLLING_HASH_MOD;\n    hs[0] = tmp;\n\n    rep(i,1,N){\n      tmp = A[i] % ((ll)ROLLING_HASH_MOD);\n      if(tmp < 0) tmp += ROLLING_HASH_MOD;\n      hs[i] = hs[i-1] + rollingHash61_mul(tmp, ROLLING_HASH_PW[i]);\n      if(hs[i] >= ROLLING_HASH_MOD) hs[i] -= ROLLING_HASH_MOD;\n    }\n  }\n\n  rollingHash get_len(int s, int len){\n    ull x;\n    rollingHash res;\n    res.len = len;\n    rollingHash_expand(s+1);\n    if(s == 0){\n      res.hs = hs[len-1];\n    } else {\n      if(hs[s+len-1] >= hs[s-1]){\n        res.hs = hs[s+len-1] - hs[s-1];\n      } else {\n        res.hs = hs[s+len-1] + ROLLING_HASH_MOD - hs[s-1];\n      }\n      res.hs = rollingHash61_mul(res.hs, ROLLING_HASH_IPW[s]);\n    }\n    return res;\n  }\n\n  rollingHash get(int a, int b){\n    return get_len(a, b - a + 1);\n  }\n\n  rollingHashSubarrays(){\n    set();\n  }\n  ~rollingHashSubarrays(){\n    free();\n  }\n};\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"max_L");
      d.push_back((string)"gcd");
      d.push_back((string)"Rand");
      d.push_back((string)"rollingHash_init");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "HashMap_init";
      string c = "{\n  int i, j, k;\n  Rand rnd;\n  rep(i,20) rnd.get(2);\n  rep(i,4){\n    rep(j,32){\n      k = rnd.get(1,62);\n      HashMap_ullP_L[i] |= (1ULL << k);\n    }\n    HashMap_ullP_L[i] |= (1ULL << 0);\n    HashMap_ullP_L[i] |= (1ULL << 63);\n  }\n}\n";
      string p = "main_first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Rand");
      need[n] = d;
    }
    {
      string n = "HashMap";
      string c = "ull HashMap_ullP_L[4];\ntemplate<class KEY, class VAL>\nstruct HashMap{\n  char *used;\n  KEY *key;\n  VAL *val;\n  int mem, n, mask;\n  int init_flag; VAL init_val;\n  HashMap(){\n    mem = 0;\n    init_flag = 0;\n  }\n  ~HashMap(){\n    free();\n  }\n  void expand(int nn){\n    if(mem >= nn) return;\n    if(mem) free();\n    mem = nn;\n    used = new char[nn];\n    key = new KEY[nn];\n    val = new VAL[nn];\n  }\n  void free(){\n    if(mem){\n      mem = 0;\n      delete[] used;\n      delete[] key;\n      delete[] val;\n    }\n  }\n  void init(int nn){\n    n = 1;\n    nn = nn + (nn + 1) / 2;\n    while(n < nn) n *= 2;\n    mask = n - 1;\n    expand(n);\n    rep(i,n) used[i] = 0;\n    init_flag = 0;\n  }\n  void init(int nn, VAL ini){\n    n = 1;\n    nn = nn + (nn + 1) / 2;\n    while(n < nn) n *= 2;\n    mask = n - 1;\n    expand(n);\n    rep(i,n) used[i] = 0;\n    init_flag = 1;\n    init_val = ini;\n  }\n  inline int getHash(const int a){\n    ull d = a;\n    d = (((d * HashMap_ullP_L[0]) >> 32) * HashMap_ullP_L[1]) & mask;\n    return d;\n  }\n  inline int getHash(const unsigned a){\n    ull d = a;\n    d = (((d * HashMap_ullP_L[0]) >> 32) * HashMap_ullP_L[1]) & mask;\n    return d;\n  }\n  inline int getHash(const ll a){\n    ull d = a;\n    d = (((((d * HashMap_ullP_L[0]) >> 32) * HashMap_ullP_L[1]) >> 32) * HashMap_ullP_L[2]) & mask;\n    return d;\n  }\n  inline int getHash(const ull a){\n    ull d = a;\n    d = (((((d * HashMap_ullP_L[0]) >> 32) * HashMap_ullP_L[1]) >> 32) * HashMap_ullP_L[2]) & mask;\n    return d;\n  }\n  inline int getHash(const pair<int,int> a){\n    ull d = (((ull)a.first) << 32) + ((ull)a.second);\n    d = (((((d * HashMap_ullP_L[0]) >> 32) * HashMap_ullP_L[1]) >> 32) * HashMap_ullP_L[2]) & mask;\n    return d;\n  }\n  inline VAL& operator[](const KEY a){\n    int k = getHash(a);\n    for(;;){\n      if(used[k]==1 && key[k]==a) break;\n      if(used[k]==0){\n        used[k] = 1, key[k] = a;\n        if(init_flag) val[k] = init_val;\n        break;\n      }\n      k = (k+1) & mask;\n    }\n    return val[k];\n  }\n  inline bool exist(const KEY a){\n    int k = getHash(a);\n    for(;;){\n      if(used[k]==1 && key[k]==a) return true;\n      if(used[k]==0) break;\n      k = (k+1) & mask;\n    }\n    return false;\n  }\n  template<class S>\n  inline bool exist(const KEY a, S &res){\n    int k = getHash(a);\n    for(;;){\n      if(used[k]==1 && key[k]==a) res = val[k], return true;\n      if(used[k]==0) break;\n      k = (k+1) & mask;\n    }\n    return false;\n  }\n};\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"HashMap_init");
      need[n] = d;
    }

    {
      string n = "subsetSum";
      string c = "template<class T, class S>\nint subsetSum(int n, T a[], S res[]){\n  int i, k, sz = 1;\n  res[0] = 0;\n  rep(k,n){\n    rep(i,sz) res[sz+i] = res[i] + a[k];\n    sz *= 2;\n  }\n  return sz;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "subsetSumS";
      string c = "template<class T, class S>\nint subsetSumS(int n, T a[], S res[], void *mem = wmem){\n  int i, k, sz = 1, bf, x, y;\n  S *arr;\n  walloc1d(&arr, 1, &mem);\n  res[0] = 0;\n  rep(k,n){\n    bf = sz;\n    sz = x = y = 0;\n    rep(i,bf) arr[i] = res[i];\n    rep(i,bf) arr[bf+i] = res[i] + a[k];\n    while(x < bf && y < bf){\n      if(arr[x] < arr[bf+y]){\n        res[sz++] = arr[x];\n        x++;\n      } else {\n        res[sz++] = arr[bf+y];\n        y++;\n      }\n    }\n    while(x < bf){\n      res[sz++] = arr[x];\n      x++;\n    }\n    while(y < bf){\n      res[sz++] = arr[bf+y];\n      y++;\n    }\n  }\n  return sz;\n}\ntemplate<class T, class S, class U>\nint subsetSumS(int n, T a[], S res[], U lim, void *mem = wmem){\n  int i, k, sz = 1, bf, x, y;\n  S *arr;\n  walloc1d(&arr, 1, &mem);\n  res[0] = 0;\n  rep(k,n){\n    bf = sz;\n    sz = x = y = 0;\n    rep(i,bf) arr[i] = res[i];\n    rep(i,bf) arr[bf+i] = res[i] + a[k];\n    while(x < bf && y < bf){\n      if(arr[x] < arr[bf+y]){\n        res[sz++] = arr[x];\n        x++;\n      } else {\n        res[sz++] = arr[bf+y];\n        y++;\n      }\n    }\n    while(x < bf){\n      res[sz++] = arr[x];\n      x++;\n    }\n    while(y < bf){\n      res[sz++] = arr[bf+y];\n      y++;\n    }\n    while(sz && res[sz-1] > lim) sz--;\n  }\n  return sz;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"workmemory");
      d.push_back((string)"walloc1d");
      need[n] = d;
    }
    {
      string n = "subsetSumSD";
      string c = "template<class T, class S>\nint subsetSumSD(int n, T a[], S res[], void *mem = wmem){\n  int i, k, sz = 1, bf, x, y;\n  S *arr;\n  walloc1d(&arr, 1, &mem);\n  res[0] = 0;\n  rep(k,n){\n    bf = sz;\n    sz = x = y = 0;\n    rep(i,bf) arr[i] = res[i];\n    rep(i,bf) arr[bf+i] = res[i] + a[k];\n    while(x < bf && y < bf){\n      if(arr[x] < arr[bf+y]){\n        if(sz==0 || res[sz-1]!=arr[x]) res[sz++] = arr[x];\n        x++;\n      } else {\n        if(sz==0 || res[sz-1]!=arr[bf+y]) res[sz++] = arr[bf+y];\n        y++;\n      }\n    }\n    while(x < bf){\n      if(sz==0 || res[sz-1]!=arr[x]) res[sz++] = arr[x];\n      x++;\n    }\n    while(y < bf){\n      if(sz==0 || res[sz-1]!=arr[bf+y]) res[sz++] = arr[bf+y];\n      y++;\n    }\n  }\n  return sz;\n}\ntemplate<class T, class S, class U>\nint subsetSumSD(int n, T a[], S res[], U lim, void *mem = wmem){\n  int i, k, sz = 1, bf, x, y;\n  S *arr;\n  walloc1d(&arr, 1, &mem);\n  res[0] = 0;\n  rep(k,n){\n    bf = sz;\n    sz = x = y = 0;\n    rep(i,bf) arr[i] = res[i];\n    rep(i,bf) arr[bf+i] = res[i] + a[k];\n    while(x < bf && y < bf){\n      if(arr[x] < arr[bf+y]){\n        if(sz==0 || res[sz-1]!=arr[x]) res[sz++] = arr[x];\n        x++;\n      } else {\n        if(sz==0 || res[sz-1]!=arr[bf+y]) res[sz++] = arr[bf+y];\n        y++;\n      }\n    }\n    while(x < bf){\n      if(sz==0 || res[sz-1]!=arr[x]) res[sz++] = arr[x];\n      x++;\n    }\n    while(y < bf){\n      if(sz==0 || res[sz-1]!=arr[bf+y]) res[sz++] = arr[bf+y];\n      y++;\n    }\n    while(sz && res[sz-1] > lim) sz--;\n  }\n  return sz;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"workmemory");
      d.push_back((string)"walloc1d");
      need[n] = d;
    }

    {
      string n = "maxSubsetDP";
      string c = "template<class T, class S>\nvoid maxSubsetDP(int N, T cost[], S res[], void *mem = wmem){\n  int i, j, k, r;\n  S **dp;\n  walloc2d(&dp, N, 1<<N, &mem);\n  rep(j,1,1<<N) dp[0][j] = cost[j];\n  rep(i,1,N) rep(j,1,1<<N) dp[i][j] = numeric_limits<S>::max();\n  rep(i,1,N) rep(j,1,1<<N){\n    r = ((1<<N)-1) ^ j;\n    for(k=r;k;k=(k-1)&r) dp[i][j|k] <?= max(dp[i-1][j], cost[k]);\n  }\n  rep(i,N) res[i] = dp[i][(1<<N)-1];\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"chmin");
      d.push_back((string)"walloc2d");
      d.push_back((string)"workmemory");
      d.push_back((string)"max_L");
      need[n] = d;
    }


    {
      string n = "graph_minColor";
      string c = "int graph_minColor(int N, int **mat, void *mem = wmem){\n  static int fg = 0, p1;\n  int i, j, k, g, logN, res, mx;\n  int *edge, *I, *hist, **hv;\n  char *bt, *bfg;\n\n  if(fg==0){\n    Rand rnd;\n    rep(i,20) rnd.get(2);\n    p1 = rnd.get(1d9-1d8, 1d9);\n    while(!isPrime(p1)) p1++;\n    fg = 1;\n  }\n\n  if(N<=1) return N;\n  rep(i,N){\n    mx = 0;\n    rep(j,N) if(i!=j && mat[i][j]) mx++;\n    if(mx == 0 || mx == N-1){\n      walloc2d(&hv, N-1, N-1, &mem);\n      rep(j,N-1) rep(k,N-1) hv[j][k] = mat[if[j>=i, j+1, j]][if[k>=i, k+1, k]];\n      return graph_minColor(N-1, hv, mem) + if[mx==0, 0, 1];\n    }\n  }\n\n  logN = 1;\n  while(N >> logN) logN++;\n\n  walloc1d(&edge, N, &mem);\n  rep(i,N) edge[i] = 0;\n  rep(i,N) rep(j,N) if(i!=j && !mat[i][j]) edge[i] |= 1<<j;\n\n  walloc1d(&I, 1<<N, &mem);\n  walloc1d(&bt, 1<<N, &mem);\n  walloc1d(&bfg, 1<<N, &mem);\n\n  rep(i,N) bt[1<<i] = i;\n  bfg[0] = 1;\n  rep(i,1,1<<N) bfg[i] = 1 - bfg[i & (i-1)];\n\n  I[0] = 1;\n  rep(i,1,1<<N){\n    k = bt[i&(-i)];\n    I[i] = I[i^(1<<k)] + I[i&edge[k]];\n  }\n  \n  k = I[(1<<N) - 1];\n  walloc1d(&hist, k, &mem);\n  rep(i,k) hist[i] = 0;\n  rep(i,1<<N){\n    if(bfg[i]) hist[I[i]-1]++;\n    else       hist[I[i]-1]--;\n  }\n\n  mx = 0;\n  rep(i,1,k) if(hist[i]) mx++;\n  walloc2d(&hv, logN, mx, &mem);\n  mx = 0;\n  rep(i,1,k) if(hist[i]){\n    hist[mx] = hist[i];\n    hv[0][mx] = i;\n    mx++;\n  }\n  rep(k,1,logN) rep(i,mx) hv[k][i] = (ll)hv[k-1][i] * hv[k-1][i] % p1;\n\n  res = 1;\n  rrep(k,logN){\n    g = 0;\n    rep(i,mx) g = (g + (ll)hist[i] * hv[k][i]) % p1;\n    if(!g){\n      res += (1<<k);\n      rep(i,mx) hist[i] = (ll)hist[i] * hv[k][i] % p1;\n    }\n  }\n  return res;\n}\n";
      string p = "first";
      vector<string> d;
      d.push_back((string)"workmemory");
      d.push_back((string)"walloc1d");
      d.push_back((string)"walloc2d");
      d.push_back((string)"Rand");
      d.push_back((string)"isPrime");
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "xorMin";
      string c = "int xorMin(int X, int N, int A[], int bt = 31){\n  int i, k;\n  rrep(i,bt){\n    if(!(A[N-1] & 1 << i)){\n      continue;\n    }\n    if(A[0] & 1 << i){\n      X ^= (1 << i);\n      continue;\n    }\n    k = lower_bound(A, A+N, A[N-1] & (~((1 << i) - 1))) - A;\n    if(X & 1 << i){\n      X ^= (1 << i);\n      A += k;\n      N -= k;\n    } else {\n      N = k;\n    }\n  }\n  return X;\n}\n\nll xorMin(ll X, int N, ll A[], int bt = 63){\n  int i, k;\n  rrep(i,bt){\n    if(!(A[N-1] & 1LL << i)){\n      continue;\n    }\n    if(A[0] & 1LL << i){\n      X ^= (1LL << i);\n      continue;\n    }\n    k = lower_bound(A, A+N, A[N-1] & (~((1LL << i) - 1))) - A;\n    if(X & 1LL << i){\n      X ^= (1LL << i);\n      A += k;\n      N -= k;\n    } else {\n      N = k;\n    }\n  }\n  return X;\n}\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "xorMax";
      string c = "int xorMax(int X, int N, int A[], int bt = 31){\n  int i, k;\n  rrep(i,bt){\n    if(!(A[N-1] & 1 << i)){\n      continue;\n    }\n    if(A[0] & 1 << i){\n      X ^= (1 << i);\n      continue;\n    }\n    k = lower_bound(A, A+N, A[N-1] & (~((1 << i) - 1))) - A;\n    if(!(X & 1 << i)){\n      X ^= (1 << i);\n      A += k;\n      N -= k;\n    } else {\n      N = k;\n    }\n  }\n  return X;\n}\n\nll xorMax(ll X, int N, ll A[], int bt = 63){\n  int i, k;\n  rrep(i,bt){\n    if(!(A[N-1] & 1LL << i)){\n      continue;\n    }\n    if(A[0] & 1LL << i){\n      X ^= (1LL << i);\n      continue;\n    }\n    k = lower_bound(A, A+N, A[N-1] & (~((1LL << i) - 1))) - A;\n    if(!(X & 1LL << i)){\n      X ^= (1LL << i);\n      A += k;\n      N -= k;\n    } else {\n      N = k;\n    }\n  }\n  return X;\n}\n";
      string p = "first";
      vector<string> d;
      name.push_back(n); func[n] = c; need[n] = d; place[n] = p;
    }

    {
      string n = "maxRectArea";
      string c = "template<class T>\nT maxRectArea(int N, T H[], void *mem = wmem){\n  int i, sz = 0;\n  T res, *st_h;\n  int *st_w;\n  if(N == 0) return 0;\n  walloc1d(&st_h, N, &mem);\n  walloc1d(&st_w, N, &mem);\n  res = 0;\n  rep(i,N){\n    if(sz == 0 || st_h[sz-1] <= H[i]){\n      st_h[sz] = H[i];\n      st_w[sz] = i;\n      sz++;\n    } else {\n      while(sz && st_h[sz-1] > H[i]){\n        res >?= st_h[sz-1] * (i - st_w[sz-1]);\n        sz--;\n      }\n      st_h[sz] = H[i];\n      sz++;\n    }\n  }\n  while(sz){\n    res >?= st_h[sz-1] * (N - st_w[sz-1]);\n    sz--;\n  }\n  return res;\n}\ntemplate<class T, class S>\nauto maxRectArea(int N, T H[], S W[], void *mem = wmem)\n-> decltype(H[0]*W[0])\n{\n  int i, sz = 0;\n  decltype(H[0]*W[0]) res, ws, *st_h, *st_w;\n  if(N == 0) return 0;\n  walloc1d(&st_h, N, &mem);\n  walloc1d(&st_w, N, &mem);\n  res = ws = 0;\n  rep(i,N){\n    if(sz == 0 || st_h[sz-1] <= H[i]){\n      st_h[sz] = H[i];\n      st_w[sz] = ws;\n      sz++;\n    } else {\n      while(sz && st_h[sz-1] > H[i]){\n        res >?= st_h[sz-1] * (ws - st_w[sz-1]);\n        sz--;\n      }\n      st_h[sz] = H[i];\n      sz++;\n    }\n    ws += W[i];\n  }\n  while(sz){\n    res >?= st_h[sz-1] * (ws - st_w[sz-1]);\n    sz--;\n  }\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"chmax");
      need[n] = d;
    }

    {
      string n = "isValidBracket1";
      string c = "int isValidBracket1(int N, char S[]){\n  int i, k = 0;\n  if(N%2) return 0;\n  rep(i,N){\n    if(S[i] == '(') k++;\n    else if(S[i] == ')') k--;\n    else return 0;\n    if(k < 0) return 0;\n  }\n  if(k==0) return 1;\n  return 0;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "isValidBracket2";
      string c = "int isValidBracket2(int N, char S[], void *mem = wmem){\n  int i, sz = 0;\n  char *s;\n  if(N%2) return 0;\n  walloc1d(&s, N, &mem);\n  rep(i,N){\n    if(S[i] == ')'){\n      if(sz == 0 || s[sz-1] != '(') return 0;\n      sz--;\n      continue;\n    }\n    if(S[i] == ']'){\n      if(sz == 0 || s[sz-1] != '[') return 0;\n      sz--;\n      continue;\n    }\n    if(S[i] == '}'){\n      if(sz == 0 || s[sz-1] != '{') return 0;\n      sz--;\n      continue;\n    }\n    if(S[i] == '>'){\n      if(sz == 0 || s[sz-1] != '<') return 0;\n      sz--;\n      continue;\n    }\n    s[sz++] = S[i];\n  }\n  if(sz==0) return 1;\n  return 0;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
    }

    {
      string n = "swapV";
      string c = "template<class T, class T1, class T2>\ninline void swapV(T &a, T1 x, T2 y){\n  if(a == x) a = y;\n  else if(a == y) a = x;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "cntSubarrayFreq";
      string c = "template<class T, class S>\nll cntSubarrayFreq(int N, T A[], S lim[], void *mem = wmem){\n  int i, j, *f;\n  T mx, mn;\n  ll res = 0;\n  if(N == 0) return 0;\n  rep(i,N) if(lim[A[i]] < 0) return 0;\n  mn = mx = A[0];\n  rep(i,1,N){\n    mn <?= A[i];\n    mx >?= A[i];\n  }\n  walloc1d(&f, mn, mx+1, &mem);\n  rep(i,N) f[A[i]] = lim[A[i]];\n  i = 0;\n  rep(j,N){\n    f[A[j]]--;\n    while(f[A[j]] < 0){\n      f[A[i]]++;\n      i++;\n    }\n    res += j - i + 1;\n  }\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"chmin");
      d.push_back((string)"chmax");
      need[n] = d;
    }
    {
      string n = "cntSubarrayDistinct";
      string c = "template<class T, class S, class S2>\nll cntSubarrayDistinct(int N, T A[], S cost[], S2 lim, void *mem = wmem){\n  int i, j, *f;\n  T mx, mn;\n  S2 cur = 0;\n  ll res = 0;\n  if(N == 0 || lim < 0) return 0;\n  mn = mx = A[0];\n  rep(i,1,N){\n    mn <?= A[i];\n    mx >?= A[i];\n  }\n  walloc1d(&f, mn, mx+1, &mem);\n  rep(i,N) f[A[i]] = 0;\n  i = 0;\n  rep(j,N){\n    if(f[A[j]] == 0) cur += cost[A[j]];\n    f[A[j]]++;\n    while(cur > lim){\n      f[A[i]]--;\n      if(f[A[i]] == 0) cur -= cost[A[i]];\n      i++;\n    }\n    res += j - i + 1;\n  }\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"chmin");
      d.push_back((string)"chmax");
      need[n] = d;
    }
    {
      string n = "maxSubarrayDistinct";
      string c = "template<class T, class S, class S2>\nint maxSubarrayDistinct(int N, T A[], S cost[], S2 lim, void *mem = wmem){\n  int i, j, *f;\n  T mx, mn;\n  S2 cur = 0;\n  int res = 0;\n  if(N == 0 || lim < 0) return 0;\n  mn = mx = A[0];\n  rep(i,1,N){\n    mn <?= A[i];\n    mx >?= A[i];\n  }\n  walloc1d(&f, mn, mx+1, &mem);\n  rep(i,N) f[A[i]] = 0;\n  i = 0;\n  rep(j,N){\n    if(f[A[j]] == 0) cur += cost[A[j]];\n    f[A[j]]++;\n    while(cur > lim){\n      f[A[i]]--;\n      if(f[A[i]] == 0) cur -= cost[A[i]];\n      i++;\n    }\n    res >?= j - i + 1;\n  }\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"chmin");
      d.push_back((string)"chmax");
      need[n] = d;
    }

    {
      string n = "vec2arr";
      string c = "template<class T, class S>\ninline int vec2arr(vector<T> &v, S arr[]){\n  int i, N = v.size();\n  rep(i,N) arr[i] = v[i];\n  return N;\n}\ntemplate<class T, class S1, class S2>\ninline int vec2arr(vector<vector<T>> &v, S1 arr1[], S2 arr2[]){\n  int i, N = v.size();\n  rep(i,N){\n    arr1[i] = v[i][0];\n    arr2[i] = v[i][1];\n  }\n  return N;\n}\ntemplate<class T, class S1, class S2, class S3>\ninline int vec2arr(vector<vector<T>> &v, S1 arr1[], S2 arr2[], S3 arr3[]){\n  int i, N = v.size();\n  rep(i,N){\n    arr1[i] = v[i][0];\n    arr2[i] = v[i][1];\n    arr3[i] = v[i][2];\n  }\n  return N;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "Rot90";
      string c = "template<class T>\nvector<vector<T>> Rot90(vector<vector<T>> a){\n  int i, j, x, y;\n  x = a.size();\n  y = a[0].size();\n  vector<vector<T>> b(y, vector<T>(x));\n  rep(i,x) rep(j,y) b[j][x-1-i] = a[i][j];\n  return b;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "Point2dEPS";
      string c = "template<class T>\nstruct Point2dEPS {\n  static T eps;\n  void set(T a){eps = a;}\n};\ntemplate<class T> T Point2dEPS<T>::eps;\nvoid Point2d_init(void){\n  Point2dEPS<float> x; x.set(1e-4);\n  Point2dEPS<double> y; y.set(1e-10);\n  Point2dEPS<long double> z; z.set(1e-10);\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Point2dEPS_init");
      need[n] = d;
    }
    {
      string n = "Point2dEPS_init";
      string c = "{Point2d_init();}\n";
      string p = "main_first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "Point2d";
      string c = "template<class T, class S, class F>\nstruct Point2d{\n  T x, y;\n  Point2d(){x = y = 0;}\n  Point2d(T a){x = a; y = 0;}\n  Point2d(T a, T b){x = a; y = b;}\n  void set(T a, T b){ x = a; y = b; }\n  Point2d<T,S,F> &operator+=(Point2d<T,S,F> a){\n    x += a.x;\n    y += a.y;\n    return *this;\n  }\n  Point2d<T,S,F> &operator-=(Point2d<T,S,F> a){\n    x -= a.x;\n    y -= a.y;\n    return *this;\n  }\n  Point2d<T,S,F> operator+(Point2d<T,S,F> a){ return Point2d<T,S,F>(*this)+=a; }\n  Point2d<T,S,F> operator-(Point2d<T,S,F> a){ return Point2d<T,S,F>(*this)-=a; }\n  inline F dist(void){\n    F tx, ty;\n    tx = x;\n    ty = y;\n    return sqrt(tx*tx + ty*ty);\n  }\n  inline F dist(Point2d<T,S,F> a){\n    F tx, ty;\n    tx = ((F)x) - ((F)a.x);\n    ty = ((F)y) - ((F)a.y);\n    return sqrt(tx*tx + ty*ty);\n  }\n  inline S dist2(void){\n    S tx, ty;\n    tx = x;\n    ty = y;\n    return tx*tx + ty*ty;\n  }\n  inline S dist2(Point2d<T,S,F> a){\n    S tx, ty;\n    tx = ((S)x) - ((S)a.x);\n    ty = ((S)y) - ((S)a.y);\n    return tx*tx + ty*ty;\n  }\n  inline F arg(void){\n    F res;\n    if(x==0 && y==0) return 0;\n    res = atan2(y, x);\n    if(res <= -PI + Point2dEPS<F>::eps) res += 2*PI;\n    return res;\n  }\n  inline F arg(Point2d<T,S,F> a){\n    F res;\n    res = arg() - a.arg();\n    if(res <= -PI + Point2dEPS<F>::eps) res += 2*PI;\n    if(res > PI + Point2dEPS<F>::eps) res -= 2*PI;\n    return res;\n  }\n};\ntemplate<class T, class S, class F> inline Point2d<T,S,F> operator*(T a, Point2d<T,S,F> b){return Point2d<T,S,F>(a*b.x, a*b.y);}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"define_PI");
      d.push_back((string)"Point2dEPS");
      need[n] = d;
    }
    {
      string n = "reader_Point2d";
      string c = "template<class T, class S, class F>\nvoid rd(Point2d<T,S,F> &a){\n  reader_ignore_error(a.x, a.y);\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Point2d";
    }
    {
      string n = "writer_Point2d";
      string c = "template<class T, class S, class F>\nvoid wt_L(Point2d<T,S,F> a){\n  wt_L(a.x);\n  wt_L(' ');\n  wt_L(a.y);\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Point2d";
    }
    {
      string n = "Point2d_InnerProd";
      string c = "template<class T, class S, class F>\nS InnerProd(Point2d<T,S,F> a, Point2d<T,S,F> b){\n  return ((S)a.x) * ((S)b.x) + ((S)a.y) * ((S)b.y);\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Point2d";
    }
    {
      string n = "Point2d_CrossProd";
      string c = "template<class T, class S, class F>\nS CrossProd(Point2d<T,S,F> a, Point2d<T,S,F> b){\n  return ((S)a.x) * ((S)b.y) - ((S)b.x) * ((S)a.y);\n}\ntemplate<class T, class S, class F>\nS CrossProd(Point2d<T,S,F> c, Point2d<T,S,F> a, Point2d<T,S,F> b){\n  S x1, x2, y1, y2;\n  x1 = ((S)a.x) - ((S)c.x);\n  y1 = ((S)a.y) - ((S)c.y);\n  x2 = ((S)b.x) - ((S)c.x);\n  y2 = ((S)b.y) - ((S)c.y);\n  return x1 * y2 - x2 * y1;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Point2d";
    }
    {
      string n = "Point2d_xysortA";
      string c = "template<class T, class S, class F>\nvoid xysortA(int N, Point2d<T,S,F> A[]){\n  sort(\n    A,\n    A+N,\n    [](auto &a, auto &b){\n      if(a.x < b.x - Point2dEPS<T>::eps) return true;\n      if(b.x < a.x - Point2dEPS<T>::eps) return false;\n      if(a.y < b.y - Point2dEPS<T>::eps) return true;\n      return false;\n    }\n  );\n}\ntemplate<class T, class S, class F, class D>\nvoid xysortA(int N, Point2d<T,S,F> A[], D ind[], void *mem = wmem){\n  int i;\n  pair<Point2d<T,S,F>,D> *arr;\n  arr = walloc1d(N,&arr,&mem);\n  rep(i,N) arr[i] = make_pair(A[i],ind[i]);\n  sort(\n    arr,\n    arr+N,\n    [](auto &a, auto &b){\n      if(a.first.x < b.first.x - Point2dEPS<T>::eps) return true;\n      if(b.first.x < a.first.x - Point2dEPS<T>::eps) return false;\n      if(a.first.y < b.first.y - Point2dEPS<T>::eps) return true;\n      return false;\n    }\n  );\n  rep(i,N) A[i] = arr[i].first, ind[i] = arr[i].second;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"workmemory");
      d.push_back((string)"walloc1d");
      need[n] = d;
      parent[n] = "Point2d";
    }
    {
      string n = "Point2d_argsortA";
      string c = "template<class T, class S, class F>\nvoid argsortA(int N, Point2d<T,S,F> A[], void *mem = wmem){\n  int i, k, c, sz[4] = {}, *sp;\n  Point2d<T,S,F> *arr[4];\n  walloc1d(&sp, N, &mem);\n  rep(i,N){\n    if(A[i].x >= -Point2dEPS<T>::eps){\n      if(A[i].y >= -Point2dEPS<T>::eps){\n        sp[i] = 2;\n      } else {\n        sp[i] = 1;\n      }\n    } else {\n      if(A[i].y >= -Point2dEPS<T>::eps){\n        sp[i] = 3;\n      } else {\n        sp[i] = 0;\n      }\n    }\n    sz[sp[i]]++;\n  }\n  rep(k,4) walloc1d(&arr[k], sz[k], &mem);\n  rep(k,4) sz[k] = 0;\n  rep(i,N){\n    k = sp[i];\n    arr[k][sz[k]++] = A[i];\n  }\n  rep(k,4) if(sz[k]){\n    sort(\n      arr[k],\n      arr[k]+sz[k],\n      [](auto &a, auto &b){\n        return CrossProd(b, a) < -Point2dEPS<S>::eps;\n      }\n    );\n  }\n  c = 0;\n  rep(k,4) rep(i,sz[k]) A[c++] = arr[k][i];\n}\ntemplate<class T, class S, class F, class D>\nvoid argsortA(int N, Point2d<T,S,F> A[], D ind[], void *mem = wmem){\n  int i, k, c, sz[4] = {}, *sp;\n  pair<Point2d<T,S,F>,D> *arr[4];\n  walloc1d(&sp, N, &mem);\n  rep(i,N){\n    if(A[i].x >= -Point2dEPS<T>::eps){\n      if(A[i].y >= -Point2dEPS<T>::eps){\n        sp[i] = 2;\n      } else {\n        sp[i] = 1;\n      }\n    } else {\n      if(A[i].y >= -Point2dEPS<T>::eps){\n        sp[i] = 3;\n      } else {\n        sp[i] = 0;\n      }\n    }\n    sz[sp[i]]++;\n  }\n  rep(k,4) walloc1d(&arr[k], sz[k], &mem);\n  rep(k,4) sz[k] = 0;\n  rep(i,N){\n    k = sp[i];\n    arr[k][sz[k]++] = make_pair(A[i], ind[i]);\n  }\n  rep(k,4) if(sz[k]){\n    sort(\n      arr[k],\n      arr[k]+sz[k],\n      [](auto &a, auto &b){\n        return CrossProd(b.first, a.first) < -Point2dEPS<S>::eps;\n      }\n    );\n  }\n  c = 0;\n  rep(k,4) rep(i,sz[k]){\n    A[c] = arr[k][i].first;\n    ind[c] = arr[k][i].second;\n    c++;\n  }\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"workmemory");
      d.push_back((string)"walloc1d");
      d.push_back((string)"Point2d_CrossProd");
      need[n] = d;
      parent[n] = "Point2d";
    }
    {
      string n = "Point2d_CCW";
      string c = "template<class T, class S, class F>\nint CCW(Point2d<T,S,F> b, Point2d<T,S,F> c){\n  S p;\n  p = CrossProd(b, c);\n  if(p < -Point2dEPS<S>::eps) return -1;\n  if(p > Point2dEPS<S>::eps) return 1;\n  p = InnerProd(b, c);\n  if(p < -Point2dEPS<S>::eps) return 2;\n  p = c.dist2() - b.dist2();\n  if(p > Point2dEPS<S>::eps) return -2;\n  return 0;\n}\ntemplate<class T, class S, class F>\nint CCW(Point2d<T,S,F> a, Point2d<T,S,F> b, Point2d<T,S,F> c){\n  S p;\n  b -= a;\n  c -= a;\n  p = CrossProd(b, c);\n  if(p < -Point2dEPS<S>::eps) return -1;\n  if(p > Point2dEPS<S>::eps) return 1;\n  p = InnerProd(b, c);\n  if(p < -Point2dEPS<S>::eps) return 2;\n  p = c.dist2() - b.dist2();\n  if(p > Point2dEPS<S>::eps) return -2;\n  return 0;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Point2d_InnerProd");
      d.push_back((string)"Point2d_CrossProd");
      need[n] = d;
      parent[n] = "Point2d";
    }
    {
      string n = "Point2d_ConvexHull";
      string c = "template<class T, class S, class F>\nint ConvexHull_sorted(int N, Point2d<T,S,F> A[], Point2d<T,S,F> res[]){\n  int i, sz = 0, t;\n  if(N <= 2){\n    rep(i,N) res[i] = A[i];\n    res[N] = A[0];\n    return N;\n  }\n  for(i=0; i<N; res[sz++] = A[i++]){\n    while(sz >= 2 && CCW(res[sz-2], res[sz-1], A[i]) <= 0) sz--;\n  }\n  t = sz;\n  for(i=N-2; i>=0; res[sz++] = A[i--]){\n    while(sz > t && CCW(res[sz-2], res[sz-1], A[i]) <= 0) sz--;\n  }\n  return sz - 1;\n}\ntemplate<class T, class S, class F>\nint ConvexHull(int N, Point2d<T,S,F> A[], Point2d<T,S,F> res[], void *mem = wmem){\n  int i;\n  Point2d<T,S,F> *arr;\n  walloc1d(&arr, N, &mem);\n  rep(i,N) arr[i] = A[i];\n  xysortA(N, arr);\n  return ConvexHull_sorted(N, arr, res);\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Point2d_CCW");
      d.push_back((string)"Point2d_xysortA");
      d.push_back((string)"workmemory");
      d.push_back((string)"walloc1d");
      need[n] = d;
      parent[n] = "Point2d";
    }
    {
      string n = "Point2d_PolygonArea2";
      string c = "template<class T, class S, class F>\nS PolygonArea2(int N, Point2d<T,S,F> A[]){\n  int i;\n  S res = 0;\n  if(N <= 2) return res;\n  rep(i,1,N) res += CrossProd(A[i-1], A[i]);\n  res += CrossProd(A[N-1], A[0]);\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Point2d_CrossProd");
      need[n] = d;
      parent[n] = "Point2d";
    }

    {
      string n = "minFactorList";
      string c = "template<class T>\nvoid minFactorList(int N, T res[]){\n  int i, j;\n  if(N <= 2){\n    rep(i,N) res[i] = i;\n    return;\n  }\n  res[0] = 0;\n  res[1] = 1;\n  for(i=2;i<N;i+=2) res[i] = 2;\n  for(i=3;i<N;i+=2) res[i] = i;\n  for(i=3;i*i<=N;i+=2) if(res[i]==i){\n    for(j=i*i;j<=N;j+=i) if(res[j]==j) res[j] = i;\n  }\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "maxFactorList";
      string c = "template<class T>\nvoid maxFactorList(int N, T res[]){\n  int i, j;\n  if(N <= 2){\n    rep(i,N) res[i] = i;\n    return;\n  }\n  res[0] = 0;\n  res[1] = 1;\n  for(i=2;i<N;i+=2) res[i] = 2;\n  for(i=3;i<N;i+=2) res[i] = i;\n  for(i=3;i*i<=N;i+=2) if(res[i]==i){\n    for(j=i*i;j<=N;j+=i) if(res[j]==j) res[j] = i;\n  }\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "FactorList";
      string c = "template<class T>\nvoid FactorList(int N, T res[], void *mem = wmem){\n  int i, j, k, *f;\n  walloc1d(&f, N, &mem);\n  maxFactorList(N, f);\n  rep(i,N) res[i] = 0;\n  rep(i,2,N){\n    k = i;\n    while(k > 1){\n      j = f[k];\n      k /= j;\n      while(f[k]==j) k /= j;\n      res[i]++;\n    }\n  }\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"maxFactorList");
      need[n] = d;
    }
    {
      string n = "FactorMList";
      string c = "template<class T>\nvoid FactorMList(int N, T res[], void *mem = wmem){\n  int i, k, *f;\n  walloc1d(&f, N, &mem);\n  maxFactorList(N, f);\n  rep(i,N) res[i] = 0;\n  rep(i,2,N){\n    k = i;\n    while(k > 1){\n      k /= f[k];\n      res[i]++;\n    }\n  }\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"maxFactorList");
      need[n] = d;
    }
    {
      string n = "EulerPhiList";
      string c = "template<class T>\nvoid EulerPhiList(int N, T res[], void *mem = wmem){\n  int i, j, k, *f;\n  walloc1d(&f, N, &mem);\n  maxFactorList(N, f);\n  rep(i,N) res[i] = i;\n  rep(i,2,N){\n    res[i] = i;\n    k = i;\n    while(k > 1){\n      j = f[k];\n      k /= j;\n      while(f[k]==j) k /= j;\n      res[i] = res[i] / j * (j-1);\n    }\n  }\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"maxFactorList");
      need[n] = d;
    }

    {
      string n = "LinearEquationMod2";
      string c = "template<class T1, class T2>\nint LinearEquationMod2(int R, int C, T1 **A, T2 *b, void *mem = wmem){\n  const int LBIT = 64;\n  const int bc = C / LBIT + 1;\n  int i, j, r, c, bcs;\n  ull **mat;\n  walloc2d(&mat, R, bc, &mem);\n  rep(i,R) rep(j,bc) mat[i][j] = 0;\n  rep(i,R) rep(j,C) if(A[i][j]) mat[i][j/LBIT] |= (1ULL<<(j%LBIT));\n  rep(i,R) if(b[i]) mat[i][C/LBIT] |= (1ULL<<(C%LBIT));\n  r = 0;\n  rep(c,C){\n    rep(i,r,R) if(mat[i][c/LBIT] & (1ULL << (c%LBIT))) break;\n    if(i==R) continue;\n    bcs = c / LBIT;\n    if(i != r){\n      rep(j,bcs,bc) swap(mat[r][j], mat[i][j]);\n    }\n    rep(i,r+1,R) if(mat[i][c/LBIT] & (1ULL << (c%LBIT))){\n      rep(j,bcs,bc) mat[i][j] ^= mat[r][j];\n    }\n    r++;\n  }\n  rep(i,r,R) if(mat[i][bc-1]) return -1;\n  return C - r;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"workmemory");
      d.push_back((string)"walloc2d");
      need[n] = d;
    }
    {
      string n = "LinearEquation";
      string c = "double LinearEquation_EPS = 1e-9;\ntemplate<class T>\nint LinearEquation_for_real(int R, int C, T **A, T *b, T *x, void *mem = wmem){\n  int i, j, p, rr = 0, cc;\n  T tmp, tmp2;\n  T **aa, *bb;\n  walloc2d(&aa, R, C, &mem);\n  walloc1d(&bb, R, &mem);\n  rep(i,R) rep(j,C) aa[i][j] = A[i][j];\n  rep(i,R) bb[i] = b[i];\n  rep(cc,C) if(rr < R){\n    tmp = -1;\n    rep(p,rr,R){\n      tmp2 = abs(aa[p][cc]);\n      if(tmp < tmp2) tmp = tmp2, i = p;\n    }\n    if(tmp < LinearEquation_EPS) continue;\n    if(i != rr){\n      rep(j,cc,C) swap(aa[rr][j], aa[i][j]);\n      swap(bb[rr], bb[i]);\n    }\n    rep(i,rr+1,R) if(aa[i][cc]){\n      tmp = aa[i][cc] / aa[rr][cc];\n      rep(j,cc,C) aa[i][j] -= tmp * aa[rr][j];\n      bb[i] -= tmp * bb[rr];\n    }\n    rr++;\n  }\n  rep(i,rr,R) if(abs(bb[i]) >= LinearEquation_EPS) return -1;\n  rep(i,C) x[i] = 0;\n  rrep(i,rr){\n    rep(j,C) if(aa[i][j] >= LinearEquation_EPS) break;\n    x[j] = bb[i];\n    rep(k,j+1,C) x[j] -= x[k] * aa[i][k];\n    x[j] /= aa[i][j];\n  }\n  return C - rr;\n}\ntemplate<class T>\nint LinearEquation_for_modint(int R, int C, T **A, T *b, T *x, void *mem = wmem){\n  int i, j, rr = 0, cc;\n  T tmp;\n  T **aa, *bb;\n  walloc2d(&aa, R, C, &mem);\n  walloc1d(&bb, R, &mem);\n  rep(i,R) rep(j,C) aa[i][j] = A[i][j];\n  rep(i,R) bb[i] = b[i];\n  rep(cc,C) if(rr < R){\n    rep(i,rr,R) if(aa[i][cc]) break;\n    \n    if(i == R) continue;\n    if(i != rr){\n      rep(j,cc,C) swap(aa[rr][j], aa[i][j]);\n      swap(bb[rr], bb[i]);\n    }\n    rep(i,rr+1,R) if(aa[i][cc]){\n      tmp = aa[i][cc] / aa[rr][cc];\n      rep(j,cc,C) aa[i][j] -= tmp * aa[rr][j];\n      bb[i] -= tmp * bb[rr];\n    }\n    rr++;\n  }\n  rep(i,rr,R) if(bb[i]) return -1;\n  rep(i,C) x[i] = 0;\n  rrep(i,rr){\n    rep(j,C) if(aa[i][j]) break;\n    x[j] = bb[i];\n    rep(k,j+1,C) x[j] -= x[k] * aa[i][k];\n    x[j] /= aa[i][j];\n  }\n  return C - rr;\n}\ntemplate<class T>\nint LinearEquation(int R, int C, T **A, T *b, T *x, void *mem = wmem){\n  return LinearEquation_for_real(R, C, A, b, x, mem);\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"workmemory");
      d.push_back((string)"walloc2d");
      d.push_back((string)"LinearEquation_Modint");
      d.push_back((string)"LinearEquation_Mint");
      d.push_back((string)"LinearEquation_modint");
      d.push_back((string)"LinearEquation_mint");
      need[n] = d;
    }
    {
      string n = "LinearEquation_Modint";
      string c = "int LinearEquation(int R, int C, Modint **A, Modint *b, Modint *x, void *mem = wmem){\n  return LinearEquation_for_modint(R, C, A, b, x, mem);\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Modint";
    }
    {
      string n = "LinearEquation_Mint";
      string c = "int LinearEquation(int R, int C, Mint **A, Mint *b, Mint *x, void *mem = wmem){\n  return LinearEquation_for_modint(R, C, A, b, x, mem);\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "Mint";
    }
    {
      string n = "LinearEquation_modint";
      string c = "int LinearEquation(int R, int C, modint **A, modint *b, modint *x, void *mem = wmem){\n  return LinearEquation_for_modint(R, C, A, b, x, mem);\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "modint";
    }
    {
      string n = "LinearEquation_mint";
      string c = "int LinearEquation(int R, int C, mint **A, mint *b, mint *x, void *mem = wmem){\n  return LinearEquation_for_modint(R, C, A, b, x, mem);\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      parent[n] = "mint";
    }

    {
      string n = "LexicographicGE";
      string c = "template<class T1, class T2, class T3>\nint LexicographicGE(int N, T1 A[], int sz, T2 available_num[], T3 res[], void *mem = wmem){\n  int i, k, *nx;\n  T1 *av;\n  walloc1d(&av, sz, &mem);\n  rep(i,sz) av[i] = available_num[i];\n  rep(i,1,sz) if(av[i-1] > av[i]){\n    sort(av, av+sz);\n    break;\n  }\n  walloc1d(&nx, N, &mem);\n  rep(i,N){\n    nx[i] = lower_bound(av, av+sz, A[i]) - av;\n    if(nx[i] == sz){\n      while(i >= 0 && nx[i] >= sz-1) i--;\n      if(i < 0){\n        rep(k,N) res[k] = av[0];\n        return 0;\n      }\n      rep(k,i) res[k] = A[k];\n      res[i] = av[nx[i]+1];\n      rep(k,i+1,N) res[k] = av[0];\n      return 1;\n    }\n    if(A[i] != av[nx[i]]){\n      rep(k,i) res[k] = A[k];\n      res[i] = av[nx[i]];\n      rep(k,i+1,N) res[k] = av[0];\n      return 1;\n    }\n  }\n  rep(i,N) res[i] = A[i];\n  return 1;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
    }

    {
      string n = "cReader_ll";
      string c = "ll cReader_ll(ll mn, ll mx, char nx){\n  int i, fg = 0, m = 1, f = -1;\n  ll res = 0; double tmp = 0;\n  for(;;){\n    i = my_getchar_unlocked();\n    if(fg==0 && i=='-'){\n      fg++;\n      m = -1;\n    } else if('0' <= i <= '9'){\n      fg++;\n      if(f == -1) f = i - '0';\n      res = 10 * res + i - '0';\n      tmp = 10 * tmp + i - '0';\n      assert(tmp < 1e20);\n    } else {\n      break;\n    }\n  }\n  assert(tmp / 2 <= res);\n  assert((m==1 && fg >= 1) || (m==-1 && fg >= 2));\n  assert(mn <= m * res <= mx);\n  assert(!(res == 0 && m == -1));\n  assert(!(res != 0 && f == 0));\n  assert(!(res == 0 && fg >= 2));\n  assert(i == nx);\n  return m * res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"my_getchar_unlocked");
      need[n] = d;
    }
    {
      string n = "cReader_eof";
      string c = "void cReader_eof(){\n  int i;\n  i = my_getchar_unlocked();\n  assert(i == EOF);\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"my_getchar_unlocked");
      need[n] = d;
    }

    {
      string n = "cntPrime";
      string c = "ll cntPrime(ll n, void *mem = wmem){\n  int i, j, k, m, sn, ssn, c;\n  char *isp;\n  ll *s1, *s2, x, tmp;\n  const double nn = n;\n  if(n <= 1) return 0;\n  if(n == 2) return 1;\n  c = 0;\n  sn = Isqrt_f(n);\n  ssn = Isqrt_f(sn);\n  walloc1d(&s1, sn+1, &mem);\n  walloc1d(&s2, sn+1, &mem);\n  walloc1d(&isp, sn+1, &mem);\n  s1[0] = 0;\n  rep(i,1,sn+1) s1[i] = i - 1;\n  rep(i,1,sn+1) s2[i] = n/i - 1;\n  rep(i,2,sn+1) isp[i] = 1;\n  rep(i,2,ssn+1) if(isp[i]){\n    for(j=i*i;j<=sn;j+=i) isp[j] = 0;\n    rep(j,1,sn+1){\n      x = (ll) i * j;\n      if(x > sn) s2[j] -= s1[(int)(nn / x)] - c;\n      else       s2[j] -= s2[x] - c;\n    }\n    for(k = (sn+1)/i; k >= i; k--){\n      m = min(sn+1, i*(k+1));\n      tmp = s1[k] - c;\n      rep(j,i*k,m) s1[j] -= tmp;\n    }\n    c++;\n  }\n  rep(i,ssn+1,sn+1) if(isp[i]){\n    rep(j,1,sn+1){\n      x = (ll) i * j;\n      if(x * i > n) break;\n      if(x > sn) s2[j] -= s1[(int)(nn / x)] - c;\n      else       s2[j] -= s2[x] - c;\n    }\n    c++;\n  }\n  return s2[1];\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"workmemory");
      d.push_back((string)"walloc1d");
      d.push_back((string)"min_L");
      d.push_back((string)"Isqrt_f");
      need[n] = d;
    }
    {
      string n = "sumPrime";
      string c = "template<class T>\nT sumPrime(ll n, void *mem = wmem){\n  int i, j, k, m, sn, ssn;\n  char *isp;\n  ll x;\n  T *s1, *s2, c, tmp;\n  const double nn = n;\n  if(n <= 1) return 0;\n  if(n == 2) return 2;\n  c = 0;\n  sn = Isqrt_f(n);\n  ssn = Isqrt_f(sn);\n  walloc1d(&s1, sn+1, &mem);\n  walloc1d(&s2, sn+1, &mem);\n  walloc1d(&isp, sn+1, &mem);\n  s1[0] = 0;\n  rep(i,1,sn+1) s1[i] = ((T)i) * ((T)(i+1)) / 2 - 1;\n  rep(i,1,sn+1) s2[i] = ((T)(n/i)) * ((T)(n/i+1)) / 2 - 1;\n  rep(i,2,sn+1) isp[i] = 1;\n  rep(i,2,ssn+1) if(isp[i]){\n    for(j=i*i;j<=sn;j+=i) isp[j] = 0;\n    rep(j,1,sn+1){\n      x = (ll) i * j;\n      if(x > sn) s2[j] -= i * (s1[(int)(nn / x)] - c);\n      else       s2[j] -= i * (s2[x] - c);\n    }\n    for(k = (sn+1)/i; k >= i; k--){\n      m = min(sn+1, i*(k+1));\n      tmp = i * (s1[k] - c);\n      rep(j,i*k,m) s1[j] -= tmp;\n    }\n    c += i;\n  }\n  rep(i,ssn+1,sn+1) if(isp[i]){\n    rep(j,1,sn+1){\n      x = (ll) i * j;\n      if(x * i > n) break;\n      if(x > sn) s2[j] -= i * (s1[(int)(nn / x)] - c);\n      else       s2[j] -= i * (s2[x] - c);\n    }\n    c += i;\n  }\n  return s2[1];\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"workmemory");
      d.push_back((string)"walloc1d");
      d.push_back((string)"min_L");
      d.push_back((string)"Isqrt_f");
      need[n] = d;
    }

    {
      string n = "Determinant_Modint";
      string c = "Modint Determinant(int n, Modint **mat, void *mem = wmem){\n  int i, j, c;\n  Modint **m;\n  Modint res = 1, tmp;\n  walloc2d(&m, n, n, &mem);\n  rep(i,n) rep(j,n) m[i][j] = mat[i][j];\n  rep(c,n){\n    rep(i,c,n) if(m[i][c]) break;\n    if(i==n) return 0;\n    if(i!=c){\n      rep(j,c,n) swap(m[i][j], m[c][j]);\n      res = -res;\n    }\n    rep(i,c+1,n) if(m[i][c]){\n      tmp = m[i][c] / m[c][c];\n      rep(j,c+1,n) m[i][j] -= tmp * m[c][j];\n    }\n    res *= m[c][c];\n  }\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"walloc2d");
      d.push_back((string)"workmemory");
      need[n] = d;
      parent[n] = "Modint";
    }
    {
      string n = "Determinant_modint";
      string c = "modint Determinant(int n, modint **mat, void *mem = wmem){\n  int i, j, c;\n  modint **m;\n  modint res = 1, tmp;\n  walloc2d(&m, n, n, &mem);\n  rep(i,n) rep(j,n) m[i][j] = mat[i][j];\n  rep(c,n){\n    rep(i,c,n) if(m[i][c]) break;\n    if(i==n) return 0;\n    if(i!=c){\n      rep(j,c,n) swap(m[i][j], m[c][j]);\n      res = -res;\n    }\n    rep(i,c+1,n) if(m[i][c]){\n      tmp = m[i][c] / m[c][c];\n      rep(j,c+1,n) m[i][j] -= tmp * m[c][j];\n    }\n    res *= m[c][c];\n  }\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"walloc2d");
      d.push_back((string)"workmemory");
      need[n] = d;
      parent[n] = "modint";
    }
    {
      string n = "Determinant_Mint";
      string c = "Mint Determinant(int n, Mint **mat, void *mem = wmem){\n  int i, j, c;\n  Mint **m;\n  Mint res = 1, tmp;\n  walloc2d(&m, n, n, &mem);\n  rep(i,n) rep(j,n) m[i][j] = mat[i][j];\n  rep(c,n){\n    rep(i,c,n) if(m[i][c]) break;\n    if(i==n) return 0;\n    if(i!=c){\n      rep(j,c,n) swap(m[i][j], m[c][j]);\n      res = -res;\n    }\n    rep(i,c+1,n) if(m[i][c]){\n      tmp = m[i][c] / m[c][c];\n      rep(j,c+1,n) m[i][j] -= tmp * m[c][j];\n    }\n    res *= m[c][c];\n  }\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"walloc2d");
      d.push_back((string)"workmemory");
      need[n] = d;
      parent[n] = "Mint";
    }
    {
      string n = "Determinant_mint";
      string c = "mint Determinant(int n, mint **mat, void *mem = wmem){\n  int i, j, c;\n  mint **m;\n  mint res = 1, tmp;\n  walloc2d(&m, n, n, &mem);\n  rep(i,n) rep(j,n) m[i][j] = mat[i][j];\n  rep(c,n){\n    rep(i,c,n) if(m[i][c]) break;\n    if(i==n) return 0;\n    if(i!=c){\n      rep(j,c,n) swap(m[i][j], m[c][j]);\n      res = -res;\n    }\n    rep(i,c+1,n) if(m[i][c]){\n      tmp = m[i][c] / m[c][c];\n      rep(j,c+1,n) m[i][j] -= tmp * m[c][j];\n    }\n    res *= m[c][c];\n  }\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"walloc2d");
      d.push_back((string)"workmemory");
      need[n] = d;
      parent[n] = "mint";
    }

    {
      string n = "arrMerge";
      string c = "template<class T1, class T2, class T3>\nint arrMerge(int As, T1 A[], int Bs, T2 B[], T3 res[]){\n  int aa = 0, bb = 0, sz = 0;\n  while(aa < As && bb < Bs){\n    if(A[aa] <= B[bb]){\n      res[sz++] = A[aa++];\n    } else {\n      res[sz++] = B[bb++];\n    }\n  }\n  while(aa < As) res[sz++] = A[aa++];\n  while(bb < Bs) res[sz++] = B[bb++];\n  return sz;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "arrMergeD";
      string c = "template<class T1, class T2, class T3>\nint arrMergeD(int As, T1 A[], int Bs, T2 B[], T3 res[]){\n  int aa = 0, bb = 0, sz = 0;\n  while(aa < As && bb < Bs){\n    if(A[aa] == B[bb]){\n      res[sz++] = A[aa++];\n      bb++;\n    } else if(A[aa] <= B[bb]){\n      res[sz++] = A[aa++];\n    } else {\n      res[sz++] = B[bb++];\n    }\n  }\n  while(aa < As) res[sz++] = A[aa++];\n  while(bb < Bs) res[sz++] = B[bb++];\n  return sz;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "opt01SubsetSum_brute";
      string c = "template<class T>\nT opt01SubsetSum_brute(int N, T A[], T t, T notfound = -1){\n  T res, tmp;\n  int ok = 0, mask, i;\n  rep(mask,1<<N){\n    tmp = 0;\n    rep(i,N) if(mask & (1<<i)) tmp += A[i];\n    if(tmp <= t){\n      if(!ok) ok = 1, res = tmp;\n      else if(res < tmp) res = tmp;\n    }\n  }\n  if(!ok) return notfound;\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "opt01SubsetSum_mim";
      string c = "template<class T>\nT opt01SubsetSum_mim(int N, T A[], T t, T notfound = -1, void *mem = wmem){\n  int i, j, n1, n2, s1, s2;\n  T *a, *a1, *a2, res, sm, ad;\n  walloc1d(&a, N, &mem);\n  rep(i,N) a[i] = A[i];\n  ad = 0;\n  rep(i,N) if(a[i] < 0) a[i] = -a[i], ad += a[i];\n  t += ad;\n  if(t < 0) return notfound;\n  sort(a, a+N);\n  while(N && a[N-1] > t) N--;\n  if(N==0){\n    if(t >= 0) return -ad;\n    return notfound;\n  }\n  sm = 0;\n  rep(i,N){\n    if(a[i] > 0) sm += a[i];\n    if(sm == t) return t - ad;\n  }\n  if(sm < t) return sm - ad;\n  n1 = N / 2;\n  n2 = N - n1;\n  walloc1d(&a1, 1<<n1, &mem);\n  walloc1d(&a2, 1<<n2, &mem);\n  s1 = subsetSumSD(n1, a, a1, t, mem);\n  s2 = subsetSumSD(n2, a+n1, a2, t, mem);\n  res = 0;\n  j = s2 - 1;\n  rep(i,s1){\n    while(j >= 0 && a1[i] + a2[j] > t) j--;\n    if(j < 0) break;\n    res >?= a1[i] + a2[j];\n  }\n  return res - ad;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"subsetSumSD");
      d.push_back((string)"chmax");
      need[n] = d;
    }
    {
      string n = "opt01SubsetSum_sdp";
      string c = "template<class T>\nT opt01SubsetSum_sdp(int N, T A[], T t, T notfound = -1, void *mem = wmem){\n  int i, k;\n  T *a, sm, ad, g;\n  char *arr;\n  walloc1d(&a, N, &mem);\n  rep(i,N) a[i] = A[i];\n  ad = 0;\n  rep(i,N) if(a[i] < 0) a[i] = -a[i], ad += a[i];\n  t += ad;\n  if(t < 0) return notfound;\n  sort(a, a+N);\n  while(N && a[N-1] > t) N--;\n  if(N==0){\n    if(t >= 0) return -ad;\n    return notfound;\n  }\n  walloc1d(&arr, t+1, &mem);\n  while(N > 0 && a[N-1] > t) N--;\n  sm = 0;\n  rep(i,N){\n    sm += a[i];\n    if(sm == t) return t - ad;\n  }\n  if(sm < t) return sm - ad;\n  g = a[0];\n  rep(i,1,N) g = gcd(g, a[i]);\n  if(g > 1){\n    rep(i,N) a[i] /= g;\n    t = fDiv(t, g);\n  }\n  arr[0] = 1;\n  rep(i,1,t+1) arr[i] = 0;\n  sm = 0;\n  rep(k,N){\n    sm = min(sm + a[k], t);\n    rrep(i,a[k],sm+1) arr[i] |= arr[i-a[k]];\n  }\n  rrep(i,t+1) if(arr[i]) return i*g - ad;\n  return notfound;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"min_L");
      d.push_back((string)"gcd");
      d.push_back((string)"fDiv");
      need[n] = d;
    }
    {
      string n = "opt01SubsetSum";
      string c = "template<class T>\nT opt01SubsetSum(int N, T A[], T t, T notfound = -1, int mem_lim = sizeof(memarr)/2, void *mem = wmem){\n  int i; T g;\n  double x;\n  double min_time;\n  double time_brute, time_mim, time_sdp;\n  double memo_brute, memo_mim, memo_sdp;\n  min_time = double_inf;\n  time_brute = pow(2.0, N) * N;\n  memo_brute = 1;\n  if(memo_brute > mem_lim) time_brute = double_inf;\n  min_time <?= time_brute;\n  time_mim = 2 * pow(2.0, (N+1)/2);\n  memo_mim = time_mim * sizeof(T);\n  if(memo_mim > mem_lim) time_mim = double_inf;\n  min_time <?= time_mim;\n  g = 0;\n  rep(i,N) g = gcd(g, abs(A[i]));\n  x = t;\n  rep(i,N) if(A[i] < 0) x -= A[i];\n  time_sdp = (double) x * N;\n  if(g > 1) time_sdp /= g;\n  memo_sdp = time_sdp * sizeof(T);\n  if(memo_sdp > mem_lim) time_sdp = double_inf;\n  min_time <?= time_sdp;\n  if(min_time == time_brute) return opt01SubsetSum_brute(N, A, t, notfound);\n  if(min_time == time_mim) return opt01SubsetSum_mim(N, A, t, notfound, mem);\n  if(min_time == time_sdp) return opt01SubsetSum_sdp(N, A, t, notfound, mem);\n  return notfound;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"opt01SubsetSum_brute");
      d.push_back((string)"opt01SubsetSum_mim");
      d.push_back((string)"opt01SubsetSum_sdp");
      d.push_back((string)"chmin");
      d.push_back((string)"gcd");
      need[n] = d;
    }
    {
      string n = "opt01SubsetSumF_brute";
      string c = "template<class T>\nT opt01SubsetSumF_brute(int N, int F, T A[], T t, T notfound = -1){\n  T res, tmp;\n  int ind[N];\n  int ok = 0, mask, i;\n  if(F < 0 || F > N) return notfound;\n  rep(i,N) ind[i] = 0;\n  rep(i,F) ind[N-1-i] = 1;\n  do{\n    tmp = 0;\n    rep(i,N) if(ind[i]) tmp += A[i];\n    if(tmp <= t){\n      if(!ok) ok = 1, res = tmp;\n      else if(res < tmp) res = tmp;\n    }\n  }while(next_permutation(ind,ind+N));\n  if(!ok) return notfound;\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "opt01SubsetSumF_mim";
      string c = "template<class T>\nT opt01SubsetSumF_mim(int N, int F, T A[], T t, T notfound = -1, void *mem = wmem){\n  int i, j, k, n1, n2, nn, *s1, *s2, *mis;\n  T *a, **a1, **a2, res, sm, ad, x, y;\n  T *arr1, *arr2; int arrsz1, arrsz2, st1, st2;\n  ll **c;\n  if(F < 0 || F > N) return notfound;\n  walloc1d(&a, N, &mem);\n  walloc1d(&mis, N, &mem);\n  rep(i,N) a[i] = A[i];\n  sort(a, a+N);\n  sm = 0;\n  rep(i,F) sm += a[i];\n  if(sm > t) return notfound;\n  sm = 0;\n  rep(i,F) sm += a[N-1-i];\n  if(sm <= t) return sm;\n  rep(i,N) a[i] = A[i];\n  rep(i,N) mis[i] = 0;\n  ad = 0;\n  rep(i,N) if(a[i] < 0){\n    a[i] = -a[i];\n    ad += a[i];\n    mis[i] = 1;\n  }\n  t += ad;\n  if(t < 0) return notfound;\n  sortA(N, a, mis, mem);\n  while(N && a[N-1] > t){\n    F -= mis[N-1];\n    N--;\n  }\n  n1 = N / 2;\n  n2 = N - n1;\n  nn = max(n1, n2);\n  walloc2d(&c, nn+1, nn+1, &mem);\n  rep(i,nn+1) c[i][0] = 1;\n  rep(j,1,nn+1) c[0][j] = 0;\n  rep(i,1,nn+1) rep(j,1,nn+1) c[i][j] = c[i-1][j-1] + c[i-1][j];\n  walloc1d(&s1, n1+1, &mem);\n  walloc1d(&a1, n1+1, &mem);\n  rep(i,n1+1) walloc1d(&a1[i], c[n1][i], &mem);\n  walloc1d(&s2, n2+1, &mem);\n  walloc1d(&a2, n2+1, &mem);\n  rep(i,n2+1) walloc1d(&a2[i], c[n2][i], &mem);\n  k = c[nn][nn/2];\n  walloc1d(&arr1, k, &mem);\n  walloc1d(&arr2, k, &mem);\n  x = 0;\n  rep(i,n1) x += mis[i];\n  y = x;\n  rep(i,n1+1) s1[i] = 0;\n  a1[x][s1[x]++] = 0;\n  rep(k,n1){\n    if(mis[k]==0){\n      rrep(j,x,y+1){\n        arrsz1 = s1[j];\n        arrsz2 = s1[j+1];\n        rep(i,arrsz1) arr1[i] = a1[j][i] + a[k];\n        rep(i,arrsz2) arr2[i] = a1[j+1][i];\n        s1[j+1] = arrMergeD(arrsz1, arr1, arrsz2, arr2, a1[j+1]);\n      }\n      y++;\n    } else {\n      rep(j,x,y+1){\n        arrsz1 = s1[j];\n        arrsz2 = s1[j-1];\n        rep(i,arrsz1) arr1[i] = a1[j][i] + a[k];\n        rep(i,arrsz2) arr2[i] = a1[j-1][i];\n        s1[j-1] = arrMergeD(arrsz1, arr1, arrsz2, arr2, a1[j-1]);\n      }\n      x--;\n    }\n  }\n  x = 0;\n  rep(i,n2) x += mis[n1+i];\n  y = x;\n  rep(i,n2+1) s2[i] = 0;\n  a2[x][s2[x]++] = 0;\n  rep(k,n2){\n    if(mis[n1+k]==0){\n      rrep(j,x,y+1){\n        arrsz1 = s2[j];\n        arrsz2 = s2[j+1];\n        rep(i,arrsz1) arr1[i] = a2[j][i] + a[n1+k];\n        rep(i,arrsz2) arr2[i] = a2[j+1][i];\n        s2[j+1] = arrMergeD(arrsz1, arr1, arrsz2, arr2, a2[j+1]);\n      }\n      y++;\n    } else {\n      rep(j,x,y+1){\n        arrsz1 = s2[j];\n        arrsz2 = s2[j-1];\n        rep(i,arrsz1) arr1[i] = a2[j][i] + a[n1+k];\n        rep(i,arrsz2) arr2[i] = a2[j-1][i];\n        s2[j-1] = arrMergeD(arrsz1, arr1, arrsz2, arr2, a2[j-1]);\n      }\n      x--;\n    }\n  }\n  res = -int_inf;\n  rep(st1,n1+1){\n    st2 = F - st1;\n    if(st2 < 0 || st2 > n2) continue;\n    j = s2[st2] - 1;\n    rep(i,s1[st1]){\n      while(j >= 0 && a1[st1][i] + a2[st2][j] > t) j--;\n      if(j < 0) break;\n      res >?= a1[st1][i] + a2[st2][j];\n    }\n  }\n  return res - ad;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"sortA");
      d.push_back((string)"walloc2d");
      d.push_back((string)"arrMergeD");
      d.push_back((string)"chmax");
      d.push_back((string)"max_L");
      need[n] = d;
    }
    {
      string n = "opt01SubsetSumF_sdp";
      string c = "template<class T>\nT opt01SubsetSumF_sdp(int N, int F, T A[], T t, T notfound = -1, void *mem = wmem){\n  int i, j, k, x, y, *mis;\n  T *a, sm, ad, *mx, g;\n  char **arr;\n  if(F < 0 || F > N) return notfound;\n  walloc1d(&a, N, &mem);\n  walloc1d(&mis, N, &mem);\n  rep(i,N) a[i] = A[i];\n  sort(a, a+N);\n  sm = 0;\n  rep(i,F) sm += a[i];\n  if(sm > t) return notfound;\n  sm = 0;\n  rep(i,F) sm += a[N-1-i];\n  if(sm <= t) return sm;\n  rep(i,N) a[i] = A[i];\n  rep(i,N) mis[i] = 0;\n  ad = 0;\n  rep(i,N) if(a[i] < 0){\n    a[i] = -a[i];\n    ad += a[i];\n    mis[i] = 1;\n  }\n  t += ad;\n  if(t < 0) return notfound;\n  sortA(N, a, mis, mem);\n  while(N && a[N-1] > t){\n    F -= mis[N-1];\n    N--;\n  }\n  g = 0;\n  rep(i,N) g = gcd(g, a[i]);\n  if(g > 1){\n    rep(i,N) a[i] /= g;\n    t = fDiv(t, g);\n  }\n  walloc2d(&arr, N+1, t+1, &mem);\n  walloc1d(&mx, N+1, &mem);\n  x = 0;\n  rep(i,N) x += mis[i];\n  y = x;\n  rep(j,N+1) rep(i,t+1) arr[j][i] = 0;\n  rep(j,N+1) mx[j] = -1;\n  arr[x][0] = 1;\n  mx[x] = 0;\n  rep(k,N){\n    if(mis[k]==0){\n      rrep(j,x,y+1){\n        mx[j+1] = min(max(mx[j+1], mx[j] + a[k]), t);\n        rep(i,a[k],mx[j+1]+1) arr[j+1][i] |= arr[j][i-a[k]];\n      }\n      y++;\n    } else {\n      rep(j,x,y+1){\n        mx[j-1] = min(max(mx[j-1], mx[j] + a[k]), t);\n        rep(i,a[k],mx[j-1]+1) arr[j-1][i] |= arr[j][i-a[k]];\n      }\n      x--;\n    }\n  }\n  rrep(i,t+1) if(arr[F][i]) return i*g - ad;\n  return notfound;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"sortA");
      d.push_back((string)"min_L");
      d.push_back((string)"max_L");
      d.push_back((string)"walloc2d");
      d.push_back((string)"gcd");
      d.push_back((string)"fDiv");
      need[n] = d;
    }
    {
      string n = "opt01SubsetSumF";
      string c = "template<class T>\nT opt01SubsetSumF(int N, int F, T A[], T t, T notfound = -1, int mem_lim = sizeof(memarr)/2, void *mem = wmem){\n  int i, k; T g;\n  double x;\n  double min_time;\n  double time_brute, time_mim, time_sdp;\n  double memo_brute, memo_mim, memo_sdp;\n  if(F < 0 || F > N) return notfound;\n  min_time = double_inf;\n  time_brute = N;\n  k = min(F, N-F);\n  rep(i,k) time_brute = time_brute * (N-i) / (i+1);\n  memo_brute = N * sizeof(int);\n  if(memo_brute > mem_lim) time_brute = double_inf;\n  min_time <?= time_brute;\n  time_mim = 2 * pow(2.0, (N+1)/2);\n  memo_mim = time_mim * sizeof(T);\n  if(memo_mim > mem_lim) time_mim = double_inf;\n  min_time <?= time_mim;\n  g = 0;\n  rep(i,N) g = gcd(g, abs(A[i]));\n  x = t;\n  rep(i,N) if(A[i] < 0) x -= A[i];\n  time_sdp = (double) x * N * (F+1);\n  if(g > 1) time_sdp /= g;\n  memo_sdp = time_sdp * sizeof(T);\n  if(memo_sdp > mem_lim) time_sdp = double_inf;\n  min_time <?= time_sdp;\n  if(min_time == time_brute) return opt01SubsetSumF_brute(N, F, A, t, notfound);\n  if(min_time == time_mim) return opt01SubsetSumF_mim(N, F, A, t, notfound, mem);\n  if(min_time == time_sdp) return opt01SubsetSumF_sdp(N, F, A, t, notfound, mem);\n  return notfound;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      d.push_back((string)"opt01SubsetSumF_brute");
      d.push_back((string)"opt01SubsetSumF_mim");
      d.push_back((string)"opt01SubsetSumF_sdp");
      d.push_back((string)"gcd");
      d.push_back((string)"chmin");
      need[n] = d;
    }

    {
      string n = "cntArrayNecessaryElement";
      string c = "template<class T>\nvoid cntArrayNecessaryElement(int N, int K, T **res){\n  int i, j;\n  rep(i,N+1) rep(j,K+1) res[i][j] = 0;\n  res[0][0] = 1;\n  rep(i,N) rep(j,K+1){\n    if(j < K) res[i+1][j] += (K - j) * res[i][j];\n    if(j > 0) res[i+1][j] += j * res[i][j-1];\n  }\n}\ntemplate<class T>\nvoid cntArrayNecessaryElement(int N, int K, T **res, T **cnt1, T **cnt2, void *mem = wmem){\n  int i, j; T *inv;\n  walloc1d(&inv, K+1, &mem);\n  rep(i,1,K+1) inv[i] = T(1) / T(i);\n  rep(i,N+1) rep(j,K+1) res[i][j] = 0;\n  rep(i,N+1) rep(j,K+1) cnt1[i][j] = 0;\n  rep(i,N+1) rep(j,K+1) cnt2[i][j] = 0;\n  res[0][0] = 1;\n  rep(i,N) rep(j,K+1){\n    if(j < K){\n      res[i+1][j] += (K - j) * res[i][j];\n      cnt1[i+1][j] += (K - j) * cnt1[i][j];\n    }\n    if(j > 0){\n      res[i+1][j] += j * res[i][j-1];\n      cnt1[i+1][j] += j * (cnt1[i][j-1] + res[i][j-1] + cnt2[i][j-1] * inv[K+1-j]);\n    }\n    cnt2[i+1][j] = res[i+1][j] * (i+1) - cnt1[i+1][j];\n  }\n  rep(i,N) rep(j,1,K+1) cnt1[i+1][j] *= inv[j];\n  rep(i,N) rep(j,K) cnt2[i+1][j] *= inv[K-j];\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc1d");
      d.push_back((string)"workmemory");
      need[n] = d;
    }
    {
      string n = "cntArrayNecessaryElement_walloc";
      string c = "template<class T>\nvoid cntArrayNecessaryElement_walloc(int N, int K, T ***res, void **mem = &wmem){\n  walloc2d(res, N+1, K+1, mem);\n  cntArrayNecessaryElement(N, K, *res);\n}\ntemplate<class T>\nvoid cntArrayNecessaryElement_walloc(int N, int K, T ***res, T ***cnt1, T ***cnt2, void **mem = &wmem){\n  walloc2d(res, N+1, K+1, mem);\n  walloc2d(cnt1, N+1, K+1, mem);\n  walloc2d(cnt2, N+1, K+1, mem);\n  cntArrayNecessaryElement(N, K, *res, *cnt1, *cnt2, *mem);\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"walloc2d");
      d.push_back((string)"cntArrayNecessaryElement");
      need[n] = d;
    }

    {
      string n = "arr2bit_int";
      string c = "template<class T>\nint arr2bit_int(const int N, const T A[], const T base = 0){\n  int i, res = 0;\n  rep(i,N) res |= (1 << (A[i] - base));\n  return res;\n}\nint arr2bit_int(const string S, const char base = 0){\n  int i, res = 0;\n  rep(i,(int)S.size()) res |= (1 << (S[i] - base));\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "arr2bit_ll";
      string c = "template<class T>\nll arr2bit_ll(const int N, const T A[], const T base = 0){\n  int i; ll res = 0;\n  rep(i,N) res |= (1LL << (A[i] - base));\n  return res;\n}\nll arr2bit_ll(const string S, const char base = 0){\n  int i; ll res = 0;\n  rep(i,(int)S.size()) res |= (1LL << (S[i] - base));\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "toLower";
      string c = "string toLower(const string S){\n  int i;\n  string res = S;\n  rep(i,res.size()) res[i] = tolower(res[i]);\n  return res;\n}\nstring toUpper(const string S){\n  int i;\n  string res = S;\n  rep(i,res.size()) res[i] = toupper(res[i]);\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "isSorted";
      string c = "template<class T>\nint isSorted(int N, const T A[]){\n  int i;\n  rep(i,1,N) if(A[i-1] > A[i]) return 0;\n  return 1;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "Slice";
      string c = "template<class T, class S>\nint Slice(const int N, const T A[], const ll bg, const ll ed, const ll step, S res[]){\n  ll i, k, s;\n  int len = 0;\n  if(N == 0) return 0;\n  s = step %% N;\n  assert(step != 0);\n  if(step > 0){\n    for(i = bg, k = bg %% N; i < ed; i += step, k += s){\n      if(k >= N) k -= N;\n      res[len++] = A[k];\n    }\n  } else {\n    for(i = bg, k = bg %% N; i > ed; i += step, k += s){\n      if(k >= N) k -= N;\n      res[len++] = A[k];\n    }\n  }\n  return len;\n}\ntemplate<class T>\nvector<T> Slice(const vector<T> &A, const ll bg, const ll ed, const ll step){\n  const int N = A.size();\n  ll i, k, s;\n  vector<T> res;\n  if(N == 0) return res;\n  s = step %% N;\n  assert(step != 0);\n  if(step > 0){\n    for(i = bg, k = bg %% N; i < ed; i += step, k += s){\n      if(k >= N) k -= N;\n      res.push_back(A[k]);\n    }\n  } else {\n    for(i = bg, k = bg %% N; i > ed; i += step, k += s){\n      if(k >= N) k -= N;\n      res.push_back(A[k]);\n    }\n  }\n  return res;\n}\nstring Slice(const string &A, const ll bg, const ll ed, const ll step){\n  const int N = A.size();\n  ll i, k, s;\n  string res;\n  if(N == 0) return res;\n  s = step %% N;\n  assert(step != 0);\n  if(step > 0){\n    for(i = bg, k = bg %% N; i < ed; i += step, k += s){\n      if(k >= N) k -= N;\n      res += A[k];\n    }\n  } else {\n    for(i = bg, k = bg %% N; i > ed; i += step, k += s){\n      if(k >= N) k -= N;\n      res += A[k];\n    }\n  }\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"moddw");
      need[n] = d;
    }

    {
      string n = "RoundUp";
      string c = "template<class T, class S>\ninline T RoundUp(T a, S b){\n  T m;\n  if(b < 0) b = -b;\n  if(b <= 1) return a;\n  m = a % b;\n  if(m == 0) return a;\n  if(m < 0) m += b;\n  return ((a + b - m) / b) * b;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "RoundDown";
      string c = "template<class T, class S>\ninline T RoundDown(T a, S b){\n  T m;\n  if(b < 0) b = -b;\n  if(b <= 1) return a;\n  m = a % b;\n  if(m == 0) return a;\n  if(m < 0) m += b;\n  return ((a - m) / b) * b;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "startWith";
      string c = "template<class T1, class T2>\nint startWith(const int As, const T1 A[], const int Bs, const T2 B[]){\n  int i;\n  if(As < Bs) return 0;\n  rep(i,Bs) if(A[i] != B[i]) return 0;\n  return 1;\n}\ntemplate<class T1, class T2>\nint startWith(const vector<T1> &A, const vector<T2> &B){\n  int i;\n  if(A.size() < B.size()) return 0;\n  rep(i,B.size()) if(A[i] != B[i]) return 0;\n  return 1;\n}\nint startWith(const string &A, const string &B){\n  int i;\n  if(A.size() < B.size()) return 0;\n  rep(i,B.size()) if(A[i] != B[i]) return 0;\n  return 1;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }
    {
      string n = "endWith";
      string c = "template<class T1, class T2>\nint endWith(const int As, const T1 A[], const int Bs, const T2 B[]){\n  int i;\n  if(As < Bs) return 0;\n  const int k = As - Bs;\n  rep(i,Bs) if(A[k+i] != B[i]) return 0;\n  return 1;\n}\ntemplate<class T1, class T2>\nint endWith(const vector<T1> &A, const vector<T2> &B){\n  int i;\n  if(A.size() < B.size()) return 0;\n  const int k = A.size() - B.size();\n  rep(i,B.size()) if(A[k+i] != B[i]) return 0;\n  return 1;\n}\nint endWith(const string &A, const string &B){\n  int i;\n  if(A.size() < B.size()) return 0;\n  const int k = A.size() - B.size();\n  rep(i,B.size()) if(A[k+i] != B[i]) return 0;\n  return 1;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "MergeTech";
      string c = "template<class T>\nvoid MergeTech(set<T> &A, set<T> &B){\n  if(A.size() < B.size()) swap(A,B);\n  for(auto x : B) A.insert(x);\n  B.clear();\n}\ntemplate<class T>\nvoid MergeTech(multiset<T> &A, multiset<T> &B){\n  if(A.size() < B.size()) swap(A,B);\n  for(auto x : B) A.insert(x);\n  B.clear();\n}\ntemplate<class T>\nvoid MergeTech(vector<T> &A, vector<T> &B){\n  if(A.size() < B.size()) swap(A,B);\n  for(auto x : B) A.push_back(x);\n  B.clear();\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "WildEQ";
      string c = "int WildEQ(const string &A, const string &B, const char wc){\n  if(A.size() != B.size()) return 0;\n  const int N = A.size();\n  for(int i = 0; i < N; i++){\n    if(A[i] != wc && B[i] != wc && A[i] != B[i]) return 0;\n  }\n  return 1;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
    }

    {
      string n = "kthPalindromicNumber64";
      string c = "__int128_t kthPalindromicNumber64(ll k){\n  int i, d = 0, d2;\n  ll cnt = 1;\n  __int128_t res = 0;\n  int dig[20];\n  while(k >= cnt){\n    k -= cnt;\n    d++;\n    if(d==1) cnt = 9;\n    else if(d%2) cnt *= 10;\n  }\n  d2 = d /+ 2;\n  rep(i,d2) dig[i] = k%10, k/=10;\n  dig[d2-1]++;\n  rrep(i,d2) res = 10 * res + dig[i];\n  rep(i,d%2,d2) res = 10 * res + dig[i];\n  return res;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"__int128_t");
      d.push_back((string)"divup");
      need[n] = d;
    }
    {
      string n = "reverseNumber";
      string c = "ll reverseNumber(ll n){\n  int i, ds, d[20];\n  ds = Digit(n, d);\n  n = 0;\n  rep(i,ds) n = 10 * n + d[i];\n  return n;\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Digit_all");
      need[n] = d;
    }
    {
      string n = "isPalindromicNumber";
      string c = "int isPalindromicNumber(int n){\n  int ds, d[10];\n  ds = Digit(n, d);\n  return isPalindrome(ds, d);\n}\nint isPalindromicNumber(ll n){\n  int ds, d[20];\n  ds = Digit(n, d);\n  return isPalindrome(ds, d);\n}\n";
      string p = "first";
      name.push_back(n); func[n] = c; place[n] = p;
      vector<string> d;
      d.push_back((string)"Digit_all");
      d.push_back((string)"isPalindrome");
      need[n] = d;
    }


    if(0){
      string n = "";
      string c = "";
      string p = "first";
      vector<string> d;
      vector<string> L;
      
      name.push_back(n);
      func[n] = c;

      need[n] = d;
      place[n] = p;
      del[n] = L;
    }

    doit.insert((string)"optimize");
    doit.insert((string)"stdc");
    doit.insert((string)"namespace");
  }


  string get_insert_string(string p){
    int i, fg;
    map<string,vector<string> >::iterator it;
    vector<string> vs;
    string res, tmp;

    for(;;){
      fg = 0;
      for(it=need.begin(); it!=need.end(); it++){
        tmp = it->first;
        vs = it->second;
        if(doit.count(tmp)==0 || g_flags.count("no-insert-"+tmp)) continue;
        if(parent.count(tmp) && doit.count(parent[tmp])==0) continue;
        rep(i,vs.size()){
          if(doit.count(vs[i])==0){
            fg++;
            doit.insert(vs[i]);
          }
        }
      }
      for(it=del.begin(); it!=del.end(); it++){
        tmp = it->first;
        vs = it->second;
        if(doit.count(tmp)==0 || g_flags.count("no-insert-"+tmp)) continue;
        if(parent.count(tmp) && doit.count(parent[tmp])==0) continue;
        rep(i,vs.size()){
          if(already.count(vs[i])==0){
            fg++;
            already.insert(vs[i]);
          }
        }
      }
      if(!fg) break;
    }

    rep(i,name.size()){
      if(parent.count(name[i]) && doit.count(parent[name[i]])==0) continue;
      if(doit.count(name[i]) && g_flags.count("no-insert-"+name[i])==0 && already.count(name[i])==0 && place[name[i]] == p){
        res += func[name[i]];
      }
    }

    return res;
  }
};

insertfunctions ifun;
set<string> used_label;


struct code{
  code *up;

  vector<string> str, strtype;
  vector<int> nxt;
  vector<code*> nxtlst;

  string name, type, loop_type, continue_label, break_label;
  std::set<string> vartype, tvartype, tmp_vartype;
  map<string,string> localvar, globalvar, argvar;

  int insert_count;

  void set_init(){
    vartype.insert((string)"void");
    vartype.insert((string)"char");
    vartype.insert((string)"signed char");
    vartype.insert((string)"unsigned char");
    vartype.insert((string)"int");
    vartype.insert((string)"signed int");
    vartype.insert((string)"unsigned int");
    vartype.insert((string)"signed");
    vartype.insert((string)"unsigned");
    vartype.insert((string)"uint");
    vartype.insert((string)"long");
    vartype.insert((string)"signed long");
    vartype.insert((string)"unsigned long");
    vartype.insert((string)"ulong");
    vartype.insert((string)"long long");
    vartype.insert((string)"signed long long");
    vartype.insert((string)"unsigned long long");
    vartype.insert((string)"VI");
    vartype.insert((string)"VII");
    vartype.insert((string)"VVI");
    vartype.insert((string)"VIII");
    vartype.insert((string)"VLL");
    vartype.insert((string)"VVLL");
    vartype.insert((string)"VVVLL");
    vartype.insert((string)"VVVI");
    vartype.insert((string)"VS");
    vartype.insert((string)"ll");
    vartype.insert((string)"ull");
    vartype.insert((string)"size_t");
    vartype.insert((string)"int8_t");
    vartype.insert((string)"int16_t");
    vartype.insert((string)"int32_t");
    vartype.insert((string)"int64_t");
    vartype.insert((string)"uint8_t");
    vartype.insert((string)"uint16_t");
    vartype.insert((string)"uint32_t");
    vartype.insert((string)"uint64_t");
    vartype.insert((string)"__int8_t");
    vartype.insert((string)"__int16_t");
    vartype.insert((string)"__int32_t");
    vartype.insert((string)"__int64_t");
    vartype.insert((string)"__int128_t");
    vartype.insert((string)"__uint8_t");
    vartype.insert((string)"__uint16_t");
    vartype.insert((string)"__uint32_t");
    vartype.insert((string)"__uint64_t");
    vartype.insert((string)"__uint128_t");
    vartype.insert((string)"float");
    vartype.insert((string)"double");
    vartype.insert((string)"long double");
    vartype.insert((string)"bool");
    vartype.insert((string)"string");
    vartype.insert((string)"FILE");
    vartype.insert((string)"Timer");
    vartype.insert((string)"Rand");
    vartype.insert((string)"fft_pnt");
    vartype.insert((string)"combination_mint");
    vartype.insert((string)"HLD");
    vartype.insert((string)"Trie");
    vartype.insert((string)"AhoCorasick");
    vartype.insert((string)"fft_pnt");
    vartype.insert((string)"Permutation");
    vartype.insert((string)"IntMap");
    vartype.insert((string)"dimcomp2");
    vartype.insert((string)"dimcomp3");
    vartype.insert((string)"dimcomp4");
    vartype.insert((string)"dimcomp5");
    vartype.insert((string)"MY_WRITER");
    vartype.insert((string)"rollingHash");
    vartype.insert((string)"rollingHashSubarrays");
    vartype.insert((string)"auto");
    vartype.insert((string)"cpp_int");
    vartype.insert((string)"Timer");
    vartype.insert((string)"Timer2");
    vartype.insert((string)"Modint");
    vartype.insert((string)"modint");
    vartype.insert((string)"Mint");
    vartype.insert((string)"mint");
    vartype.insert((string)"unionFind");
    vartype.insert((string)"graph");


    tvartype.insert((string)"pair");
    tvartype.insert((string)"tuple");
    tvartype.insert((string)"vector");
    tvartype.insert((string)"stack");
    tvartype.insert((string)"queue");
    tvartype.insert((string)"deque");
    tvartype.insert((string)"priority_queue");
    tvartype.insert((string)"set");
    tvartype.insert((string)"multiset");
    tvartype.insert((string)"map");
    tvartype.insert((string)"bitset");
    tvartype.insert((string)"unordered_set");
    tvartype.insert((string)"unordered_map");
    tvartype.insert((string)"Heap");
    tvartype.insert((string)"Heap_max");
    tvartype.insert((string)"LHeap");
    tvartype.insert((string)"twoMultisets");
    tvartype.insert((string)"segtree_Point_Minval");
    tvartype.insert((string)"segtree_Point_Maxval");
    tvartype.insert((string)"segtree_Point_Min");
    tvartype.insert((string)"segtree_Point_SumMin");
    tvartype.insert((string)"segtree_Point_Prod");
    tvartype.insert((string)"segtree_Point_Minval2");
    tvartype.insert((string)"segtree_Point_Or");
    tvartype.insert((string)"segtree_Point_And");
    tvartype.insert((string)"segtree_Point_Xor");
    tvartype.insert((string)"segtree_Add_Minval");
    tvartype.insert((string)"segtree_ChangeAdd_Sum");
    tvartype.insert((string)"segtree_ChangeP1add_Sum");
    tvartype.insert((string)"segtree_pg");
    tvartype.insert((string)"segtree_ph");
    tvartype.insert((string)"segtree_rg");
    tvartype.insert((string)"segtree_rh");
    tvartype.insert((string)"static_segtree_Add_At");
    tvartype.insert((string)"HLD_fenwick");
    tvartype.insert((string)"HLD_segtree");
    tvartype.insert((string)"Matrix");
    tvartype.insert((string)"AhoCorasick_Sum");
    tvartype.insert((string)"Grid1d");
    tvartype.insert((string)"Grid2d");
    tvartype.insert((string)"Comb");
    tvartype.insert((string)"Polynomial");
    tvartype.insert((string)"weightedUnionFind");
    tvartype.insert((string)"rollbackUnionFind");
    tvartype.insert((string)"DijkstraHeap");
    tvartype.insert((string)"fenwick");
    tvartype.insert((string)"fenwick_xor");
    tvartype.insert((string)"segtree");
    tvartype.insert((string)"segtree_Change_Sum");
    tvartype.insert((string)"wgraph");
    tvartype.insert((string)"maxflow");
    tvartype.insert((string)"minCostFlow");
    tvartype.insert((string)"Arr1d");
    tvartype.insert((string)"Arr2d");
    tvartype.insert((string)"rangeTree2d");
    tvartype.insert((string)"rangeTree2d_nw");
    tvartype.insert((string)"rangeTree2d_pf");
    tvartype.insert((string)"HashMap");
    tvartype.insert((string)"Point2d");

  }


  void merge(std::set<string> &a, std::set<string> &b){
    std::set<string>::iterator it;
    for(it=b.begin();it!=b.end();it++){
      a.insert(*it);
    }
  }
  void merge(map<string,string> &a, map<string,string> &b){
    map<string,string>::iterator it;
    for(it=b.begin();it!=b.end();it++){
      a[it->first] = it->second;
    }
  }
  void set_next_types(std::set<string> &s_vartype, set<string> &s_tvartype, map<string,string> &s_localvar, map<string,string> &s_globalvar, map<string,string> &s_argvar){
    merge(vartype, s_vartype);
    merge(tvartype, s_tvartype);
    merge(globalvar, s_globalvar);
    merge(globalvar, s_localvar);
    merge(globalvar, s_argvar);
  }

  void setUpnode(code *cc){
    up = cc;
    set_next_types(cc->vartype, cc->tvartype, cc->localvar, cc->globalvar, cc->argvar);
  }


  int is_empty_block(void){
    int i, j;
    int res = 1;

    if(localvar.size()) res = 0;
    if(nxtlst.size()) res = 0;
    rep(i,(int)str.size()) if(str[i]!=";") res = 0;

    return res;
  }

  code* get_root(void){
    code *r = this;
    while(r->up != NULL) r = r->up;
    return r;
  }

  pair<code*, string> find_struct(string in){
    int i;
    string tmp;
    pair<code*, string> res;

    res.first = NULL;
    res.second = "";

    rep(i,(int)str.size()){
      if(strtype[i] == "block-struct"){
        tmp = str[i];
        if(tmp == "struct " + in){
          res.first = nxtlst[nxt[i]];
          res.second = str[i];
          return res;
        }
      }

      if(strtype[i] == "block-class"){
        tmp = str[i];
        if(tmp == "class " + in){
          res.first = nxtlst[nxt[i]];
          res.second = str[i];
          return res;
        }
      }

      if(strtype[i].substr(0,5) == "block"){
        res = nxtlst[nxt[i]]->find_struct(in);
        if(res.first != NULL) return res;
      }
    }

    return res;
  }

  void add_vartype(string in){
    int i;
    vartype.insert(in);
    rep(i,nxtlst.size()) nxtlst[i]->add_vartype(in);
  }

  void add_tvartype(string in){
    int i;
    tvartype.insert(in);
    rep(i,nxtlst.size()) nxtlst[i]->add_tvartype(in);
  }
  
  void add_localvar(string in1, string in2){
    int i;
    code *cc = this;

    while(cc->type == "block-inserted") cc = cc->up;
    cc->localvar[in1] = in2;
    rep(i,cc->nxtlst.size()) cc->nxtlst[i]->add_globalvar(in1, in2);
  }
  
  void add_globalvar(string in1, string in2){
    int i;
    globalvar[in1] = in2;
    rep(i,nxtlst.size()) nxtlst[i]->add_globalvar(in1, in2);
  }
  
  void add_argvar(string in1, string in2){
    int i;
    argvar[in1] = in2;
    rep(i,nxtlst.size()) nxtlst[i]->add_globalvar(in1, in2);
  }


  void ftrim(string &in){
    while(in.size() && isspace(in[0])) in = in.substr(1);
  }
  void etrim(string &in){
    while(in.size() && isspace(in[in.size()-1])) in = in.substr(0, in.size()-1);
  }
  void trim(string &in){
    ftrim(in);
    etrim(in);
  }
  void alltrim(string &in){
    int i;
    string t;
    rep(i,in.size()) if(!isspace(in[i])) t += in[i];
    in = t;
  }

  void ftrim_until(string &in, char s){
    while(in.size() && in[0] != s) in = in.substr(1);
  }
  void etrim_until(string &in, char t){
    while(in.size() && in[in.size()-1] != t) in = in.substr(0, in.size()-1);
  }
  void trim_until(string &in, char s, char t){
    ftrim_until(in,s);
    etrim_until(in,t);
  }

  int MD_get_inv(ll a, int md){ll t=a,s=md,u=1,v=0,e;while(s){e=t/s;t-=e*s;u-=e*v;swap(t,s);swap(u,v);}if(u<0)u+=md;return u;}


  void code_replace(string &in){
    int i, j, k, dot, dotp, d, dp, ex;
    string str;
    char buf[30];
    ull val;
    
    trim(in);
    replaceAll_ns_t(in, "ll", "long long");
    replaceAll_ns_t(in, "ull", "unsigned long long");
    replaceAll_ns_t(in, "int_inf", "1073709056");
    replaceAll_ns_t(in, "ll_inf", "4611686016279904256LL");
    replaceAll_ns_t(in, "double_inf", "1e150");
    replaceAll_ns_t(in, "VI", "vector<int>");
    replaceAll_ns_t(in, "VII", "vector<vector<int>>");
    replaceAll_ns_t(in, "VVI", "vector<vector<int>>");
    replaceAll_ns_t(in, "VIII", "vector<vector<vector<int>>>");
    replaceAll_ns_t(in, "VVVI", "vector<vector<vector<int>>>");
    replaceAll_ns_t(in, "VLL", "vector<long long>");
    replaceAll_ns_t(in, "VVLL", "vector<vector<long long>>");
    replaceAll_ns_t(in, "VVVLL", "vector<vector<vector<long long>>>");
    replaceAll_ns_t(in, "VS", "vector<string>");

    if(strpos_ns_t(in, (string)"scanf") >= 0 || strpos_ns_t(in, (string)"cin") >= 0){
      g_flags.insert((string)"no-fread");
    }
    if(strpos_ns_t(in, (string)"printf") >= 0 || strpos_ns_t(in, (string)"cout") >= 0 || strpos_ns_t(in, (string)"puts") >= 0){
      g_flags.insert((string)"no-fwrite");
    }
    
    if(strpos_ns_t(in, (string)"__int128_t") >= 0) ifun.doit.insert((string)"__int128_t");
    if(strpos_ns_t(in, (string)"__uint128_t") >= 0) ifun.doit.insert((string)"__uint128_t");
    if(strpos_ns_t(in, (string)"MD") >= 0) ifun.doit.insert((string)"define_MD");
    if(strpos_ns_t(in, (string)"PI") >= 0) ifun.doit.insert((string)"define_PI");
    if(strpos_ns_t(in, (string)"Timer") >= 0) ifun.doit.insert((string)"Timer");
    if(strpos_ns_t(in, (string)"Rand") >= 0) ifun.doit.insert((string)"Rand");
    if(strpos_ns_t(in, (string)"combination_mint") >= 0) ifun.doit.insert((string)"combination_mint");
    if(strpos_ns_t(in, (string)"Heap") >= 0) ifun.doit.insert((string)"Heap");
    if(strpos_ns_t(in, (string)"Heap_max") >= 0) ifun.doit.insert((string)"Heap_max");
    if(strpos_ns_t(in, (string)"LHeap") >= 0) ifun.doit.insert((string)"LHeap");
    if(strpos_ns_t(in, (string)"twoMultisets") >= 0) ifun.doit.insert((string)"twoMultisets");
    if(strpos_ns_t(in, (string)"segtree_Point_Minval") >= 0) ifun.doit.insert((string)"segtree_Point_Minval");
    if(strpos_ns_t(in, (string)"segtree_Point_Maxval") >= 0) ifun.doit.insert((string)"segtree_Point_Maxval");
    if(strpos_ns_t(in, (string)"segtree_Point_Min") >= 0) ifun.doit.insert((string)"segtree_Point_Min");
    if(strpos_ns_t(in, (string)"segtree_Point_SumMin") >= 0) ifun.doit.insert((string)"segtree_Point_SumMin");
    if(strpos_ns_t(in, (string)"segtree_Point_Prod") >= 0) ifun.doit.insert((string)"segtree_Point_Prod");
    if(strpos_ns_t(in, (string)"segtree_Point_Minval2") >= 0) ifun.doit.insert((string)"segtree_Point_Minval2");
    if(strpos_ns_t(in, (string)"segtree_Point_Or") >= 0) ifun.doit.insert((string)"segtree_Point_Or");
    if(strpos_ns_t(in, (string)"segtree_Point_And") >= 0) ifun.doit.insert((string)"segtree_Point_And");
    if(strpos_ns_t(in, (string)"segtree_Point_Xor") >= 0) ifun.doit.insert((string)"segtree_Point_Xor");
    if(strpos_ns_t(in, (string)"segtree_Add_Minval") >= 0) ifun.doit.insert((string)"segtree_Add_Minval");
    if(strpos_ns_t(in, (string)"segtree_ChangeAdd_Sum") >= 0) ifun.doit.insert((string)"segtree_ChangeAdd_Sum");
    if(strpos_ns_t(in, (string)"segtree_ChangeP1add_Sum") >= 0) ifun.doit.insert((string)"segtree_ChangeP1add_Sum");
    if(strpos_ns_t(in, (string)"segtree_pg") >= 0) ifun.doit.insert((string)"segtree_pg");
    if(strpos_ns_t(in, (string)"segtree_ph") >= 0) ifun.doit.insert((string)"segtree_ph");
    if(strpos_ns_t(in, (string)"segtree_rg") >= 0) ifun.doit.insert((string)"segtree_rg");
    if(strpos_ns_t(in, (string)"segtree_rh") >= 0) ifun.doit.insert((string)"segtree_rh");
    if(strpos_ns_t(in, (string)"static_segtree_Add_At") >= 0) ifun.doit.insert((string)"static_segtree_Add_At");
    if(strpos_ns_t(in, (string)"HLD") >= 0) ifun.doit.insert((string)"HLD");
    if(strpos_ns_t(in, (string)"HLD_fenwick") >= 0) ifun.doit.insert((string)"HLD_fenwick");
    if(strpos_ns_t(in, (string)"HLD_segtree") >= 0) ifun.doit.insert((string)"HLD_segtree");
    if(strpos_ns_t(in, (string)"Matrix") >= 0) ifun.doit.insert((string)"Matrix");
    if(strpos_ns_t(in, (string)"Permutation") >= 0) ifun.doit.insert((string)"Permutation");
    if(strpos_ns_t(in, (string)"IntMap") >= 0) ifun.doit.insert((string)"IntMap");
    if(strpos_ns_t(in, (string)"Trie") >= 0) ifun.doit.insert((string)"Trie");
    if(strpos_ns_t(in, (string)"AhoCorasick") >= 0) ifun.doit.insert((string)"AhoCorasick");
    if(strpos_ns_t(in, (string)"AhoCorasick_Sum") >= 0) ifun.doit.insert((string)"AhoCorasick_Sum");
    if(strpos_ns_t(in, (string)"Grid1d") >= 0) ifun.doit.insert((string)"Grid1d");
    if(strpos_ns_t(in, (string)"Grid2d") >= 0) ifun.doit.insert((string)"Grid2d");
    if(strpos_ns_t(in, (string)"dimcomp2") >= 0) ifun.doit.insert((string)"dimcomp2");
    if(strpos_ns_t(in, (string)"dimcomp3") >= 0) ifun.doit.insert((string)"dimcomp3");
    if(strpos_ns_t(in, (string)"dimcomp4") >= 0) ifun.doit.insert((string)"dimcomp4");
    if(strpos_ns_t(in, (string)"dimcomp5") >= 0) ifun.doit.insert((string)"dimcomp5");
    if(strpos_ns_t(in, (string)"rollingHash") >= 0) ifun.doit.insert((string)"rollingHash");
    if(strpos_ns_t(in, (string)"rollingHashSubarrays") >= 0) ifun.doit.insert((string)"rollingHash");
    if(strpos_ns_t(in, (string)"cLtraits_identity") >= 0) ifun.doit.insert((string)"cLtraits_identity");
    if(strpos_ns_t(in, (string)"cLtraits_try_make_signed") >= 0) ifun.doit.insert((string)"cLtraits_try_make_signed");
    if(strpos_ns_t(in, (string)"cLtraits_try_make_unsigned") >= 0) ifun.doit.insert((string)"cLtraits_try_make_unsigned");
    if(strpos_ns_t(in, (string)"cLtraits_common_type") >= 0) ifun.doit.insert((string)"cLtraits_common_type");
    if(strpos_ns_t(in, (string)"Timer") >= 0) ifun.doit.insert((string)"Timer");
    if(strpos_ns_t(in, (string)"Timer2") >= 0) ifun.doit.insert((string)"Timer2");
    if(strpos_ns_t(in, (string)"Modint") >= 0) ifun.doit.insert((string)"Modint");
    if(strpos_ns_t(in, (string)"modint") >= 0) ifun.doit.insert((string)"modint");
    if(strpos_ns_t(in, (string)"Mint") >= 0) ifun.doit.insert((string)"Mint");
    if(strpos_ns_t(in, (string)"mint") >= 0) ifun.doit.insert((string)"mint");
    if(strpos_ns_t(in, (string)"Comb") >= 0) ifun.doit.insert((string)"Comb");
    if(strpos_ns_t(in, (string)"Polynomial") >= 0) ifun.doit.insert((string)"Polynomial");
    if(strpos_ns_t(in, (string)"unionFind") >= 0) ifun.doit.insert((string)"unionFind");
    if(strpos_ns_t(in, (string)"weightedUnionFind") >= 0) ifun.doit.insert((string)"weightedUnionFind");
    if(strpos_ns_t(in, (string)"rollbackUnionFind") >= 0) ifun.doit.insert((string)"rollbackUnionFind");
    if(strpos_ns_t(in, (string)"DijkstraHeap") >= 0) ifun.doit.insert((string)"DijkstraHeap");
    if(strpos_ns_t(in, (string)"fenwick") >= 0) ifun.doit.insert((string)"fenwick");
    if(strpos_ns_t(in, (string)"fenwick_xor") >= 0) ifun.doit.insert((string)"fenwick_xor");
    if(strpos_ns_t(in, (string)"segtree") >= 0) ifun.doit.insert((string)"segtree");
    if(strpos_ns_t(in, (string)"segtree_Change_Sum") >= 0) ifun.doit.insert((string)"segtree_Change_Sum");
    if(strpos_ns_t(in, (string)"graph") >= 0) ifun.doit.insert((string)"graph");
    if(strpos_ns_t(in, (string)"wgraph") >= 0) ifun.doit.insert((string)"wgraph");
    if(strpos_ns_t(in, (string)"maxflow") >= 0) ifun.doit.insert((string)"maxflow");
    if(strpos_ns_t(in, (string)"minCostFlow") >= 0) ifun.doit.insert((string)"minCostFlow");
    if(strpos_ns_t(in, (string)"Arr1d") >= 0) ifun.doit.insert((string)"Arr1d");
    if(strpos_ns_t(in, (string)"Arr2d") >= 0) ifun.doit.insert((string)"Arr2d");
    if(strpos_ns_t(in, (string)"rangeTree2d") >= 0) ifun.doit.insert((string)"rangeTree2d");
    if(strpos_ns_t(in, (string)"rangeTree2d_nw") >= 0) ifun.doit.insert((string)"rangeTree2d_nw");
    if(strpos_ns_t(in, (string)"rangeTree2d_pf") >= 0) ifun.doit.insert((string)"rangeTree2d_pf");
    if(strpos_ns_t(in, (string)"HashMap") >= 0) ifun.doit.insert((string)"HashMap");
    if(strpos_ns_t(in, (string)"subsetSum") >= 0) ifun.doit.insert((string)"subsetSum");
    if(strpos_ns_t(in, (string)"subsetSumS") >= 0) ifun.doit.insert((string)"subsetSumS");
    if(strpos_ns_t(in, (string)"subsetSumSD") >= 0) ifun.doit.insert((string)"subsetSumSD");
    if(strpos_ns_t(in, (string)"Point2d") >= 0) ifun.doit.insert((string)"Point2d");
    if(strpos_ns_t(in, (string)"opt01SubsetSum_brute") >= 0) ifun.doit.insert((string)"opt01SubsetSum_brute");
    if(strpos_ns_t(in, (string)"opt01SubsetSum_mim") >= 0) ifun.doit.insert((string)"opt01SubsetSum_mim");
    if(strpos_ns_t(in, (string)"opt01SubsetSum_sdp") >= 0) ifun.doit.insert((string)"opt01SubsetSum_sdp");
    if(strpos_ns_t(in, (string)"opt01SubsetSum") >= 0) ifun.doit.insert((string)"opt01SubsetSum");
    if(strpos_ns_t(in, (string)"opt01SubsetSumF_brute") >= 0) ifun.doit.insert((string)"opt01SubsetSumF_brute");
    if(strpos_ns_t(in, (string)"opt01SubsetSumF_mim") >= 0) ifun.doit.insert((string)"opt01SubsetSumF_mim");
    if(strpos_ns_t(in, (string)"opt01SubsetSumF_sdp") >= 0) ifun.doit.insert((string)"opt01SubsetSumF_sdp");
    if(strpos_ns_t(in, (string)"opt01SubsetSumF") >= 0) ifun.doit.insert((string)"opt01SubsetSumF");

    if(strpos_ns_t(in, (string)"cpp_int") >= 0) ifun.doit.insert((string)"BoostMultiprecision");


    for(;;){
      int k4 = 0, k5 = 0;
      int fg = 0;
      rep(i,in.size()){
        if(k5 == 0 && in[i]=='\'') k4 ^= 1;
        if(k4 == 0 && in[i]=='"') k5 ^= 1;
        if(k4 || k5) continue;
        if(!isdigit(in[i])) continue;
        if(i && isalnum(in[i-1])) continue;

        dot = d = ex = 0;
        REP(j,i,in.size()){
          if(in[j]=='d'){ d++; dp=j; if(in[j+1]=='-' || in[j+1]=='+') j++; continue; }
          if(in[j]=='.'){ dot++; dotp=j; continue; }
          if(!isdigit(in[j])) break;
        }
        if(j < in.size() && isalpha(in[j])) continue;
        if(dot > 1 || d != 1) continue;
        if(dot && dotp > dp) continue;
        if(j == dp+1) continue;

        str = in.substr(i, dp-i);
        ex = atoi(in.c_str() + dp + 1);

        if(ex > 30) ex = 30;
        if(ex < -30) ex = -30;
        while(ex > 0){
          rep(k,str.size()) if(str[k]=='.') break;
          if(k==str.size()){
            str += '0';
          } else if(k==str.size()-1) {
            str[k] = '0';
          } else {
            str = str.substr(0,k) + str.substr(k+1,1) + "." +  str.substr(k+2);
          }
          ex--;
        }
        while(ex < 0){
          rep(k,str.size()) if(str[k]=='.') break;
          if(k!=str.size()) str = str.substr(0, k);
          if(str.size()) str = str.substr(0, str.size()-1);
          ex++;
        }
        
        val = 0;
        rep(k,str.size()){
          if(str[k]=='.') break;
          val = val * 10 + (str[k] - '0');
        }
        
        if(val < (1ULL << 31)){
          sprintf(buf, "%llu", val);
        } else if(val < (1ULL << 63)) {
          sprintf(buf, "%lluLL", val);
        } else {
          sprintf(buf, "%lluULL", val);
        }
        str = buf;
        in = in.substr(0,i) + str + in.substr(j);
        
        fg = 1;
        break;
      }
      if(!fg) break;
    }
  }

  vector<string> split_p(string in, char c){
    int i;
    int k1 = 0, k2 = 0, k3 = 0, k4 = 0, k5 = 0;
    vector<string> res;
    string tmp;

    rep(i,in.size()){
      if(k4==0 && k5==0){
        if(in[i]=='(') k1++;
        if(in[i]==')') k1--;
        if(in[i]=='[') k2++;
        if(in[i]==']') k2--;
        if(in[i]=='{') k3++;
        if(in[i]=='}') k3--;
      }
      if(k5 == 0 && in[i]=='\'') k4 ^= 1;
      if(k4 == 0 && in[i]=='"') k5 ^= 1;
      if(k1==0 && k2==0 && k3==0 && k4==0 && k5==0 && in[i]==c){
        res.push_back(tmp);
        tmp = "";
      } else {
        tmp += in[i];
      }
    }

    res.push_back(tmp);
    return res;
  }

  vector<string> split_p2(string in, char c){
    int i;
    int k1 = 0, k2 = 0, k3 = 0, k4 = 0, k5 = 0, k6 = 0;
    vector<string> res;
    string tmp;

    rep(i,in.size()){
      if(k4==0 && k5==0){
        if(in[i]=='(') k1++;
        if(in[i]==')') k1--;
        if(in[i]=='[') k2++;
        if(in[i]==']') k2--;
        if(in[i]=='{') k3++;
        if(in[i]=='}') k3--;
      }
      if(k2==0 && in[i]=='<' && in.substr(i,3) != "<?=") k6++;
      if(k2==0 && in[i]=='>' && in.substr(i,3) != ">?=") k6--;
      if(k5 == 0 && in[i]=='\'') k4 ^= 1;
      if(k4 == 0 && in[i]=='"') k5 ^= 1;
      if(k1==0 && k2==0 && k3==0 && k4==0 && k5==0 && k6==0 && in[i]==c){
        res.push_back(tmp);
        tmp = "";
      } else {
        tmp += in[i];
      }
    }

    res.push_back(tmp);
    return res;
  }

  vector<string> split_p3(string in, char c){
    int i, j, k;
    int k1 = 0, k2 = 0, k3 = 0, k4 = 0, k5 = 0;
    vector<string> res;
    string tmp;

    rep(i,in.size()){
      if(i==0 || (!isalnum(in[i-1]) && in[i-1] != '_')){
        j = -1;
        REP(k,i,in.size()){
          if(isValidVarType(in.substr(i,k-i), in[k])) j = k;
        }
        if(j >= 0){
          tmp += in.substr(i,j-i);
          i = j - 1;
          continue;
        }
      }
      if(k4==0 && k5==0){
        if(in[i]=='(') k1++;
        if(in[i]==')') k1--;
        if(in[i]=='[') k2++;
        if(in[i]==']') k2--;
        if(in[i]=='{') k3++;
        if(in[i]=='}') k3--;
      }
      if(k5 == 0 && in[i]=='\'') k4 ^= 1;
      if(k4 == 0 && in[i]=='"') k5 ^= 1;
      if(k1==0 && k2==0 && k3==0 && k4==0 && k5==0 && in[i]==c){
        res.push_back(tmp);
        tmp = "";
      } else {
        tmp += in[i];
      }
    }

    res.push_back(tmp);
    return res;
  }

  int getBlockLength(string in){
    int i = 0;
    int k1 = 0, k2 = 0, k3 = 0, k4 = 0, k5 = 0;
    rep(i,in.size()){
      if(k4==0 && k5==0){
        if(in[i]=='(') k1++;
        if(in[i]==')') k1--;
        if(in[i]=='[') k2++;
        if(in[i]==']') k2--;
        if(in[i]=='{') k3++;
        if(in[i]=='}'){
          k3--;
          if(k3==0) break;
        }
      }
      if(k5==0 && in[i]=='\'') k4 ^= 1;
      if(k4==0 && in[i]=='"') k5 ^= 1;
      if(k1==0 && k2==0 && k3==0 && k4==0 && k5==0 && in[i]==';') break;
    }
    return i+1;
  }

  int isValidVarType(string in, char nxt){
    int i, j, cnt;

    if(in.size()==0) return 0;
    if(in[in.size()-1] != '>' && isalnum(nxt)) return 0;
    if( (isalnum(in[in.size()-1]) || in[in.size()-1]=='_') && (isalnum(nxt) || nxt=='_') ) return 0;

    for(;;){
      if(in.substr(0,1) == "~"){
        in = in.substr(1);
        ftrim(in);
        continue;
      }
      
      if(in.substr(0,6) == "const "){
        in = in.substr(6);
        ftrim(in);
        continue;
      }
      if(in.substr(0,7) == "static "){
        in = in.substr(7);
        ftrim(in);
        continue;
      }
      if(in.substr(0,9) == "register "){
        in = in.substr(9);
        ftrim(in);
        continue;
      }
      if(in.substr(0,9) == "register "){
        in = in.substr(9);
        ftrim(in);
        continue;
      }
      if(in.substr(0,7) == "extern "){
        in = in.substr(7);
        ftrim(in);
        continue;
      }
      if(in.substr(0,9) == "typename "){
        in = in.substr(9);
        ftrim(in);
        continue;
      }
      if(in.substr(0,10) == "inplace_L "){
        in = in.substr(10);
        ftrim(in);
        continue;
      }
      break;
    }

    //for(string hoge : tmp_vartype) fprintf(stderr, "tmpvar [%s] [%s]\n", hoge.c_str(), in.c_str());

    if(vartype.count(in) || tmp_vartype.count(in)) return 1;
    if(tvartype.count(in)) return 1;

    //fprintf(stderr, "dame\n");

    if(in.substr(0,8) == "decltype"){
      REP(j,8,in.size()) if(in[j]=='(') break;
      if(j < in.size()){
        j = pairBracket(in,j);
        if(j != -1){
          j++;
          while(j < in.size()){
            if(!isspace(in[j])) break;
            j++;
          }
          if(j == in.size()) return 1;
        }
      }
    }

    rep(i,in.size()) if(!isalnum(in[i]) && in[i]!='.' && in[i]!='_') break;
    if(i < in.size() && in[i]=='<'){
      string tmp;
      j = pairBracketIq(in, i);
      if(j >= 0) tmp = in.substr(j+1);
      trim(tmp);
      if(tmp.size()){
        while(tmp.substr(0,6)=="::type") tmp = tmp.substr(6);
        if(tmp.size() == 0) return 1;
      }
    }

    rep(i,in.size()){
      if(tvartype.count( in.substr(0,i) )){
        cnt = 0;
        REP(j,i,in.size()){
          if(cnt==0 && isalnum(in[j])) break;
          if(in[j]=='<') cnt++;
          if(in[j]=='>') cnt--;
          if(cnt==0 && in[j]=='>'){
            j++;
            if(j==in.size()) return 1;
            if(in.substr(j) == "::iterator") return 1;
            break;
          }
        }
      }
    }

    return 0;
  }

  string isValidVarType_string(string in, char nxt){ /* valid なら type の string を返す */
    int i, j, cnt;

    if(in.size()==0) return "";
    if(in[in.size()-1] != '>' && isalnum(nxt)) return "";
    if( (isalnum(in[in.size()-1]) || in[in.size()-1]=='_') && (isalnum(nxt) || nxt=='_') ) return "";

    for(;;){
      if(in.substr(0,1) == "~"){
        in = in.substr(1);
        ftrim(in);
        continue;
      }
      
      if(in.substr(0,6) == "const "){
        in = in.substr(6);
        ftrim(in);
        continue;
      }
      if(in.substr(0,7) == "static "){
        in = in.substr(7);
        ftrim(in);
        continue;
      }
      if(in.substr(0,9) == "register "){
        in = in.substr(9);
        ftrim(in);
        continue;
      }
      if(in.substr(0,9) == "register "){
        in = in.substr(9);
        ftrim(in);
        continue;
      }
      if(in.substr(0,7) == "extern "){
        in = in.substr(7);
        ftrim(in);
        continue;
      }
      if(in.substr(0,9) == "typename "){
        in = in.substr(9);
        ftrim(in);
        continue;
      }
      if(in.substr(0,10) == "inplace_L "){
        in = in.substr(10);
        ftrim(in);
        continue;
      }
      break;
    }
    return in;
  }

  pair<string, char> nextToken(string &in){
    int i, dig = 0;
    string res1 = "";
    char res2 = '\0';

    i = 0;
    while(i < in.size() && isspace(in[i])) i++;
    if(i==in.size()) return make_pair(res1,res2);

    if(isdigit(in[0]) || (in[0]=='.' && isdigit(in[1]))) dig = 1;
    while(i < in.size() && (isalnum(in[i]) || in[i]=='.' || (!dig && in[i]=='_'))) res1 += in[i++];
//    while(i < in.size() && (isalnum(in[i]) || (dig && in[i]=='.') || (!dig && in[i]=='_'))) res1 += in[i++];
    while(i < in.size() && isspace(in[i])) i++;
    res2 = in[i];
    return make_pair(res1,res2);
  }

  pair<string,string> var_definition(string in){
    int i, j, tt;
    string type, name, pre, suf, eq;
    vector<string> tmp;

    rep(i,in.size()+1) if(isValidVarType(in.substr(0,i), in[i])) tt = i;

    type = in.substr(0, tt);
    
    in = in.substr(tt);
    trim(in);

    if(in[in.size()-1] == ';'){
      in = in.substr(0, in.size()-1);
      trim(in);
    }

    //fprintf(stderr, "[%s]\n", in.c_str());
    tmp = split_p(in, '=');
    while(tmp.size() > 2){
      string a, b;
      a = tmp[tmp.size() - 2];
      b = tmp[tmp.size() - 1];
      tmp.pop_back();
      tmp.pop_back();
      tmp.push_back(a + "=" + b);
    }
    if(tmp.size()==2){
      eq = tmp[1];
      trim(eq);
      eq = equation_main(eq);
      in = tmp[0];
    }

    rep(i,in.size()) if(isalpha(in[i]) || in[i]=='_' || in[i]==':') break;
    REP(j,i,in.size()) if(!(isalnum(in[j]) || in[j]=='_' || in[j]==':')) break;
    name = in.substr(i, j-i);
    pre = in.substr(0, i);
    suf = in.substr(j);
    trim(name);
    trim(pre);
    trim(suf);
    
    return make_pair(name, type+","+pre+","+suf+","+eq);
  }

  string getElementalyVarType(string in){
    int i, k;
    vector<string> vs;

    k = strpos(in, ".");
    if(k >= 0){
      string tp, res;
      code *cc;
      pair<code*,string> rs;
      
      tp = getElementalyVarType( in.substr(0,k) );
      in = in.substr(k+1);
      trim(in);

      //fprintf(stderr, "--- %s - %s\n", tp.c_str(), in.c_str());

      if(tp == "graph"){
        if(in=="es" || in=="edge" || in=="N") return "int";
        return "error";
      }
      if(tp == "maxflow"){
        // todo
      }

      cc = get_root();
      rs = cc->find_struct(tp);
      if(rs.first == NULL) return "error";
      res = rs.first->getElementalyVarType(in);
      return res;
    }

    while(in.size() && !isalpha(in[0])) in = in.substr(1);
    if(in.size()==0) return (string)"error";
    rep(i,in.size()) if(!isalnum(in[i]) && in[i]!='_') break;
    in = in.substr(0, i);

    if(in == "MD") return "int";
    if(in == "PI") return "double";

    if(localvar.count(in)){
      in = localvar[in];
    } else if(argvar.count(in)){
      in = argvar[in];
    } else if(globalvar.count(in)){
      in = globalvar[in];
    } else {
      return (string)"error";
    }

    vs = split_p(in, ',');
    in = vs[0];

    for(;;){
      if(in.substr(0,6) == "const "){
        in = in.substr(6);
        ftrim(in);
        continue;
      }
      if(in.substr(0,7) == "static "){
        in = in.substr(7);
        ftrim(in);
        continue;
      }
      if(in.substr(0,8) == "typename "){
        in = in.substr(8);
        ftrim(in);
        continue;
      }
      if(in.substr(0,10) == "inplace_L "){
        in = in.substr(10);
        ftrim(in);
        continue;
      }
      break;
    }

    for(;;){
      int fg = 0;
      if(in.substr(0,5) == "Arr1d"){
        trim_until(in, '<', '>');
        in = in.substr(1, in.size()-2);
        trim(in);
        fg = 1;
      }
      if(in.substr(0,5) == "Arr2d"){
        trim_until(in, '<', '>');
        in = in.substr(1, in.size()-2);
        trim(in);
        fg = 1;
      }
      if(in.substr(0,6) == "Grid1d"){
        trim_until(in, '<', '>');
        in = in.substr(1, in.size()-2);
        trim(in);
        fg = 1;
      }
      if(in.substr(0,6) == "Grid2d"){
        trim_until(in, '<', '>');
        in = in.substr(1, in.size()-2);
        trim(in);
        fg = 1;
      }
      if(in.substr(0,6) == "vector"){
        trim_until(in, '<', '>');
        in = in.substr(1, in.size()-2);
        trim(in);
        fg = 1;
      }
      if(in.substr(0,6) == "Matrix"){
        trim_until(in, '<', '>');
        in = in.substr(1, in.size()-2);
        trim(in);
        fg = 1;
      }
      if(in == "Permutation"){
        in = "int";
        fg = 1;
      }
      if(in == "IntMap"){
        in = "int";
        fg = 1;
      }
      if(!fg) break;
    }
    return in;
  }

  string getVarType(string in){
    int i, k;
    vector<string> vs;

    k = strpos(in, ".");
    if(k >= 0){
      string tp, res;
      code *cc;
      pair<code*,string> rs;
      
      tp = getElementalyVarType( in.substr(0,k) );
      in = in.substr(k+1);
      trim(in);

      //fprintf(stderr, "--- %s - %s\n", tp.c_str(), in.c_str());

      if(tp == "graph"){
        if(in=="es" || in=="edge" || in=="N") return "int";
        return "error";
      }
      if(tp == "maxflow"){
        // todo
      }

      cc = get_root();
      rs = cc->find_struct(tp);
      if(rs.first == NULL) return "error";
      res = rs.first->getElementalyVarType(in);
      return res;
    }

    while(in.size() && !isalpha(in[0])) in = in.substr(1);
    if(in.size()==0) return (string)"error";
    rep(i,in.size()) if(!isalnum(in[i]) && in[i]!='_') break;
    in = in.substr(0, i);

    if(in == "MD") return "int";
    if(in == "PI") return "double";

    if(localvar.count(in)){
      in = localvar[in];
    } else if(argvar.count(in)){
      in = argvar[in];
    } else if(globalvar.count(in)){
      in = globalvar[in];
    } else {
      return (string)"error";
    }

    vs = split_p(in, ',');
    in = vs[0];

    for(;;){
      if(in.substr(0,6) == "const "){
        in = in.substr(6);
        ftrim(in);
        continue;
      }
      if(in.substr(0,7) == "static "){
        in = in.substr(7);
        ftrim(in);
        continue;
      }
      if(in.substr(0,8) == "typename "){
        in = in.substr(8);
        ftrim(in);
        continue;
      }
      if(in.substr(0,10) == "inplace_L "){
        in = in.substr(10);
        ftrim(in);
        continue;
      }
      break;
    }

    return in;
  }

  string getEquationType(string in, int *fg = NULL, int *char_fg = NULL){
    int i;
    map<string, int> flags;
    pair<string,char> stchar;
    vector<string> do_typename;
    string tp, str;

    for(;;){
      trim(in);
      if(in.size()==0) break;

      if(in[0] == '['){
        i = pairBracket(in, 0);
        if(i == -1) return (string)"error";
        in = in.substr(i+1);
        continue;
      }

      if(in[0] == '"'){
        flags["char"] = 1;
        if(char_fg != NULL) char_fg++;
        for(i=1;;i++){
          if(in[i]=='\\'){ i++; continue; }
          if(in[i]=='"') break;
        }
        in = in.substr(i+1);
        continue;
      }

      if(in[0] == '\''){
        flags["char"] = 1;
        if(char_fg != NULL) char_fg++;
        for(i=1;;i++) if(in[i]=='\'') break;
        in = in.substr(i+1);
        continue;
      }

      stchar = nextToken(in);
      str = stchar.first;
      //fprintf(stderr,"[%s - %c]\n", stchar.first.c_str(), stchar.second);
      if((stchar.second == '?' || stchar.second == ':') && fg != NULL) (*fg)++;
      in = in.substr(str.size());
      if(str.size()==0) in = in.substr(1);

      if(str.size()){
        if(isdigit(str[0]) || (str[0]=='.' && isdigit(str[1]))) {
          if(strpos(str, (string)".") >= 0 || strpos(str,(string)"e") >= 0 || strpos(str,(string)"E") >= 0){
            if(str[str.size()-1] == 'f' || str[str.size()-1] == 'F') tp = "float";
            else                                                     tp = "double";
          } else {
            if(str.size() >= 3 && (str.substr(str.size()-3) == "ull" ||  str.substr(str.size()-3) == "ULL"))    tp = "unsigned long long";
            else if(str.size() >= 2 && (str.substr(str.size()-2) == "ll" ||  str.substr(str.size()-2) == "LL")) tp = "long long";
            else if(str.size() >= 1 && (str.substr(str.size()-1) == "u" ||  str.substr(str.size()-1) == "U"))   tp = "unsigned";
            else                                                                                                tp = "int";
          }
        } else if(isalpha(str[0])) {
          tp = getElementalyVarType(str);
        }
        //fprintf(stderr, "%s\n",tp.c_str());

        flags[tp] = 1;
        if(tp == "char" && char_fg != NULL) (*char_fg)++;
        if(tp == "error" && fg!=NULL) (*fg)++;
        if(tp.substr(0,8) == "cLtraits" && fg!=NULL) (*fg)++;
        if(tp.substr(0,16) == "remove_reference" && fg!=NULL) (*fg)++;
        if(tp.substr(0,8) == "decltype" && fg!=NULL) (*fg)++;
      }
    }

    if(flags["string"]) return "string";
    if(flags["double"]) return "double";
    if(flags["float"]) return "float";
    if(flags["Modint"]) return "Modint";
    if(flags["Mint"]) return "Mint";
    if(flags["modint"]) return "modint";
    if(flags["mint"]) return "mint";
    if(flags["__uint128_t"]) return "__uint128_t";
    if(flags["__int128_t"]) return "__int128_t";
    if(flags["unsigned long long"]) return "unsigned long long";
    if(flags["ull"]) return "unsigned long long";
    if(flags["long long"]) return "long long";
    if(flags["ll"]) return "long long";
    if(flags["uint64_t"]) return "uint64_t";
    if(flags["__uint64_t"]) return "__uint64_t";
    if(flags["int64_t"]) return "int64_t";
    if(flags["__int64_t"]) return "__int64_t";
    if(flags["unsigned"]) return "unsigned";
    if(flags["unsigned int"]) return "unsigned";
    if(flags["uint"]) return "unsigned";
    if(flags["int"]) return "int";
    if(flags["signed int"]) return "int";
    if(flags["signed"]) return "int";
    if(flags["uint32_t"]) return "uint32_t";
    if(flags["__uint32_t"]) return "__uint32_t";
    if(flags["int32_t"]) return "int32_t";
    if(flags["__int32_t"]) return "__int32_t";
    if(flags["uint16_t"]) return "uint16_t";
    if(flags["__uint16_t"]) return "__uint16_t";
    if(flags["int16_t"]) return "int16_t";
    if(flags["__int16_t"]) return "__int16_t";
    if(flags["char"]) return "char";
    if(flags["uint8_t"]) return "uint8_t";
    if(flags["__uint8_t"]) return "__uint8_t";
    if(flags["int8_t"]) return "int8_t";
    if(flags["__int8_t"]) return "__int8_t";
    if(flags["void"]) return "void";
    if(fg != NULL) (*fg)++;
    if(char_fg != NULL) (*char_fg)++;
    return (string)"error";
  }

  string getUnusedVarName(void){
    int i, k, p, r;
    string f = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_";
    string res;

    for(k=8;;k++) rep(p,100){
      res = "";

      r = rand()%52;
      res += f[r];
      REP(i,1,k){
        r = rand()%f.size();
        res += f[r];
      }

      if(localvar.count(res)) continue;
      if(argvar.count(res)) continue;
      if(globalvar.count(res)) continue;
      return res;
    }

    return res;
  }

  void insert(string &in, int p){
    int pure;
    for(;;){
      ftrim(in);
      if(in.size()==0) break;

      insert_count++;
      if(in[0]=='{') pure = 1; else  pure = 0;
      
      code *cc = new code;
      cc->setUpnode(this);
      if(pure){
        cc->set(in, (string)"block");
      } else {
        cc->set(in, (string)"block-inserted");
      }
      nxt.insert(nxt.begin()+p, nxtlst.size());
      nxtlst.push_back(cc);
      str.insert(str.begin()+p, (string)"");
      if(pure){
        strtype.insert(strtype.begin()+p, (string)"block");
      } else {
        strtype.insert(strtype.begin()+p, (string)"block-inserted");
      }

      p++;
    }
  }

  int strpos(string &A, string B, int st = 0){
    int i, As, Bs;
    As = A.size();
    Bs = B.size();
    REP(i,st,As-Bs+1) if(A.substr(i,Bs) == B) return i;
    return -1;
  }

  int strpos_ns(string &A, string B, int st = 0){ // "", '' の中は検索対象外

    int i, k4, k5, As, Bs;
    As = A.size();
    Bs = B.size();
    k4 = k5 = 0;
    REP(i,st,As-Bs+1){
      if( (k4 || k5) && A[i] == '\\' ){i++; continue;}
      if(k5 == 0 && A[i] == '"') k4 ^= 1;
      if(k4 == 0 && A[i] == '\'') k5 ^= 1;
      if(k4==0 && k5==0 && A.substr(i,Bs) == B) return i;
    }
    return -1;
  }

  int strrpos_ns(string &A, string B, int st = -1){
    int i, k4, k5, As, Bs;
    As = A.size();
    Bs = B.size();
    k4 = k5 = 0;

    if(st == -1) st = As-Bs+1;
    for(i=st;i>=0;i--){
      if( (k4 || k5) && A[i] == '\\' ){i++; continue;}
      if(k5 == 0 && A[i] == '"') k4 ^= 1;
      if(k4 == 0 && A[i] == '\'') k5 ^= 1;
      if(k4==0 && k5==0 && A.substr(i,Bs) == B) return i;
    }
    return -1;
  }

  int strpos_t(string &A, string B, int st = 0){
    int i, As, Bs;
    As = A.size();
    Bs = B.size();
    REP(i,st,As-Bs+1) if(A.substr(i,Bs) == B){
      if(i && (isalnum(A[i-1]) || A[i-1]=='_')) continue;
      if(i+Bs<As && (isalnum(A[i+Bs]) || A[i+Bs]=='_')) continue;
      return i;
    }
    return -1;
  }

  int strpos_ns_t(string &A, string B, int st = 0){
    int i, k4, k5, As, Bs;
    As = A.size();
    Bs = B.size();
    k4 = k5 = 0;
    REP(i,st,As-Bs+1){
      if( (k4 || k5) && A[i] == '\\' ){i++; continue;}
      if(k5 == 0 && A[i] == '"') k4 ^= 1;
      if(k4 == 0 && A[i] == '\'') k5 ^= 1;
      if(k4==0 && k5==0 && A.substr(i,Bs) == B){
        if(i && (isalnum(A[i-1]) || A[i-1]=='_')) continue;
        if(i+Bs<As && (isalnum(A[i+Bs]) || A[i+Bs]=='_')) continue;
        return i;
      }
    }
    return -1;
  }

  int replaceAll(string &A, string B, string C){
    int i, res = 0;
    for(;;){
      i = strpos(A, B);
      if(i==-1) break;
      A = A.substr(0, i) + C + A.substr(i+B.size());
    }
    return res;
  }

  int replaceAll_ns(string &A, string B, string C){
    int i, res = 0;
    for(;;){
      i = strpos_ns(A, B);
      if(i==-1) break;
      A = A.substr(0, i) + C + A.substr(i+B.size());
    }
    return res;
  }

  int replaceAll_t(string &A, string B, string C){
    int i, res = 0;
    for(;;){
      i = strpos_t(A, B);
      if(i==-1) break;
      A = A.substr(0, i) + C + A.substr(i+B.size());
    }
    return res;
  }

  int replaceAll_ns_t(string &A, string B, string C){
    int i, res = 0;
    for(;;){
      i = strpos_ns_t(A, B);
      if(i==-1) break;
      A = A.substr(0, i) + C + A.substr(i+B.size());
    }
    return res;
  }

  int replaceAll_ns_t_with_leading_digit(string &A, string &B, string &C){
    int i, j, k, ls, k4, k5, As, Bs;
    As = A.size();
    Bs = B.size();
    k4 = k5 = 0;
    rep(i,As-Bs+1){
      if( (k4 || k5) && A[i] == '\\' ){i++; continue;}
      if(k5 == 0 && A[i] == '"') k4 ^= 1;
      if(k4 == 0 && A[i] == '\'') k5 ^= 1;
      if(k4==0 && k5==0 && A.substr(i,Bs) == B){
        //if(i && (isalnum(A[i-1]) || A[i-1]=='_')) continue;
        if(i+Bs<As && (isalnum(A[i+Bs]) || A[i+Bs]=='_')) continue;
        j = k = 0; ls = 1;
        while(i-j-1 >= 0 && (isdigit(A[i-j-1]) || isspace(A[i-j-1]))){
          ls = 0;
          if(isdigit(A[i-j-1])) k++, ls = 1;
          j++;
        }
        //fprintf(stderr, "[%s] [%s] [%s] : %d %d %d %d\n", A.c_str(), B.c_str(), C.c_str(), i, j, k, ls);
        if(ls != 0 && (i-j && (isalnum(A[i-j-1]) || A[i-j-1]=='_'))) continue;
        if(k) A = A.substr(0,i) + "*" + C + A.substr(i+B.size());
        else  A = A.substr(0,i) + C + A.substr(i+B.size());
        return replaceAll_ns_t_with_leading_digit(A, B, C) + 1;
      }
    }
    return 0;
  }

  vector<string> rd_wt_array(string in){ // 1+((a,b)(N)+2) returns "1+(", "(a,b)", "N", "+2)", hoge returns "hoge"
    int i, j, k, bf = -1, bf_ind = -1, fg, lasf = 0;
    string chk;
    vector<string> res;

    rep(i,in.size()){
      fg = 0;
      if(bf == ')'){
        string type_ka;
        fg = 1;
        j = pairBracket(in, bf_ind);
        type_ka = in.substr(j+1, bf_ind-j-1);
        trim(type_ka);
        if(isValidVarType(type_ka, ' ')) fg = 0;
      }
      if(bf == ']') fg = lasf;
      if(isalnum(bf) && in[i] == '('){
        if(localvar.count(chk) || globalvar.count(chk) || argvar.count(chk)) fg = 1;
        //fprintf(stderr, "rd_wt [%s] [%d]\n", chk.c_str(), fg);
        if(getElementalyVarType(chk) != "error") fg = 1;
      }
      if(isalnum(bf) && in[i] == '[' && lasf == 0){
        if(localvar.count(chk) || globalvar.count(chk) || argvar.count(chk)) lasf = 1;
        //fprintf(stderr, "rd_wt [%s] [%d]\n", chk.c_str(), fg);
        if(getElementalyVarType(chk) != "error") lasf = 1;
      }
      if(fg && in[i] == '('){
        rep(k,4) res.push_back((string)"");
        j = pairBracket(in, i);
        res[2] = in.substr(i+1, j-i-1);
        res[3] = in.substr(j+1);
        in = in.substr(0, i);

        trim(in);
        i = in.size() - 1;
        while(in[i] == ']'){
          i = pairBracket(in, i) - 1;
          while(isspace(in[i])) i--;
        }
        if(in[i]==')'){
          i = pairBracket(in, i);
          res[0] = in.substr(0, i);
          res[1] = in.substr(i);
        } else {
          while(i >= 0 && (isalnum(in[i]) || in[i]=='_' || in[i]=='@' || in[i]=='.')) i--;
          i++;
          res[0] = in.substr(0, i);
          res[1] = in.substr(i);
        }

        rep(k,4) trim(res[k]);
        return res;
      }
      if(!isspace(in[i])){
        bf = in[i];
        bf_ind = i;
        if(isalnum(bf) || bf=='.' || bf=='_') chk += bf;
        else                                  chk = "";
      }
    }

    res.push_back(in);
    return res;
  }

  int checkBracketsCoressponding(string &in){
    int i;
    int k4, k5;
    stack<int> s;

    k4 = k5 = 0;
    rep(i,in.size()){
      if(k4==0 && k5==0){
        if(in[i] == '(') s.push(1);
        if(in[i] == '[') s.push(2);
        if(in[i] == '{') s.push(3);
        if(in[i] == ')'){
          if(s.size()==0 || s.top()!=1) return 0;
          s.pop();
        }
        if(in[i] == ']'){
          if(s.size()==0 || s.top()!=2) return 0;
          s.pop();
        }
        if(in[i] == '}'){
          if(s.size()==0 || s.top()!=3) return 0;
          s.pop();
        }
      }
      if(k5==0 && in[i] == '\'') k4^=1;
      if(k4==0 && in[i] == '"') k5^=1;
      if( (k4||k5) && in[i] == '\\') {i++; continue;}
    }
    if(k4 || k5 || s.size()) return 0;
    return 1;
  }

  // 最初と最後は<>だが，内部で()とかの中の<>はカッコでなく不等号とみなすもの
  int pairBracketIq2(string &in, int p){
    int i;
    int k1, k2, k3, k4, k5, k6;
    string st;

    k1 = k2 = k3 = k4 = k5 = k6 = 0;
    if(in[p]=='(' || in[p]=='[' || in[p]=='{' || in[p] == '<'){
      for(i=p;i<in.size();i++){
        if(k4==0 && k5==0){
          if(in[i] == '(') st += in[i];
          if(in[i] == ')') st += in[i];
          if(in[i] == '[') st += in[i];
          if(in[i] == ']') st += in[i];
          if(in[i] == '{') st += in[i];
          if(in[i] == '}') st += in[i];
          if((st.size() == 0 || st[st.size()-1] == '>') && in[i] == '<') st += in[i];
          if((st.size() == 0 || st[st.size()-1] == '<') && in[i] == '>') st += in[i];
          if(st.size() >= 2 && st.substr(st.size()-2) == "()") st = st.substr(0, st.size()-2);
          if(st.size() >= 2 && st.substr(st.size()-2) == "[]") st = st.substr(0, st.size()-2);
          if(st.size() >= 2 && st.substr(st.size()-2) == "{}") st = st.substr(0, st.size()-2);
          if(st.size() >= 2 && st.substr(st.size()-2) == "<>") st = st.substr(0, st.size()-2);
        }
        if(k5==0 && in[i] == '\'') k4^=1;
        if(k4==0 && in[i] == '"') k5^=1;
        if( (k4||k5) && in[i] == '\\') {i++; continue;}
        if(st.size()==0) return i;
      }
    } else if(in[p]==')' || in[p]==']' || in[p]=='}' || in[p] == '>'){
      for(i=p;i>=0;i--){
        if(k4==0 && k5==0){
          if(in[i] == '(') st += in[i];
          if(in[i] == ')') st += in[i];
          if(in[i] == '[') st += in[i];
          if(in[i] == ']') st += in[i];
          if(in[i] == '{') st += in[i];
          if(in[i] == '}') st += in[i];
          if((st.size() == 0 || st[st.size()-1] == '>') && in[i] == '<') st += in[i];
          if((st.size() == 0 || st[st.size()-1] == '<') && in[i] == '>') st += in[i];
          if(st.size() >= 2 && st.substr(st.size()-2) == ")(") st = st.substr(0, st.size()-2);
          if(st.size() >= 2 && st.substr(st.size()-2) == "][") st = st.substr(0, st.size()-2);
          if(st.size() >= 2 && st.substr(st.size()-2) == "}{") st = st.substr(0, st.size()-2);
          if(st.size() >= 2 && st.substr(st.size()-2) == "><") st = st.substr(0, st.size()-2);
        }
        if(k5==0 && in[i] == '\'') k4^=1;
        if(k4==0 && in[i] == '"') k5^=1;
        if(st.size()==0) return i;
      }
    }
    return -1;
  }

  int pairBracketIq(string &in, int p){
    int i;
    int k1, k2, k3, k4, k5, k6;
    string st;

    k1 = k2 = k3 = k4 = k5 = k6 = 0;
    if(in[p]=='(' || in[p]=='[' || in[p]=='{' || in[p] == '<'){
      for(i=p;i<in.size();i++){
        if(k4==0 && k5==0){
          if(in[i] == '(') st += in[i];
          if(in[i] == ')') st += in[i];
          if(in[i] == '[') st += in[i];
          if(in[i] == ']') st += in[i];
          if(in[i] == '{') st += in[i];
          if(in[i] == '}') st += in[i];
          if(in[i] == '<') st += in[i];
          if(in[i] == '>') st += in[i];
          if(st.size() >= 2 && st.substr(st.size()-2) == "()") st = st.substr(0, st.size()-2);
          if(st.size() >= 2 && st.substr(st.size()-2) == "[]") st = st.substr(0, st.size()-2);
          if(st.size() >= 2 && st.substr(st.size()-2) == "{}") st = st.substr(0, st.size()-2);
          if(st.size() >= 2 && st.substr(st.size()-2) == "<>") st = st.substr(0, st.size()-2);
        }
        if(k5==0 && in[i] == '\'') k4^=1;
        if(k4==0 && in[i] == '"') k5^=1;
        if( (k4||k5) && in[i] == '\\') {i++; continue;}
        if(st.size()==0) return i;
      }
    } else if(in[p]==')' || in[p]==']' || in[p]=='}' || in[p] == '>'){
      for(i=p;i>=0;i--){
        if(k4==0 && k5==0){
          if(in[i] == '(') st += in[i];
          if(in[i] == ')') st += in[i];
          if(in[i] == '[') st += in[i];
          if(in[i] == ']') st += in[i];
          if(in[i] == '{') st += in[i];
          if(in[i] == '}') st += in[i];
          if(in[i] == '<') st += in[i];
          if(in[i] == '>') st += in[i];
          if(st.size() >= 2 && st.substr(st.size()-2) == ")(") st = st.substr(0, st.size()-2);
          if(st.size() >= 2 && st.substr(st.size()-2) == "][") st = st.substr(0, st.size()-2);
          if(st.size() >= 2 && st.substr(st.size()-2) == "}{") st = st.substr(0, st.size()-2);
          if(st.size() >= 2 && st.substr(st.size()-2) == "><") st = st.substr(0, st.size()-2);
        }
        if(k5==0 && in[i] == '\'') k4^=1;
        if(k4==0 && in[i] == '"') k5^=1;
        if(st.size()==0) return i;
      }
    }
    return -1;
  }

  int pairBracket(string &in, int p){
    int i;
    int k1, k2, k3, k4, k5;

    k1 = k2 = k3 = k4 = k5 = 0;
    if(in[p]=='(' || in[p]=='[' || in[p]=='{'){
      for(i=p;i<in.size();i++){
        if(k4==0 && k5==0){
          if(in[i] == '(') k1++;
          if(in[i] == ')') k1--;
          if(in[i] == '[') k2++;
          if(in[i] == ']') k2--;
          if(in[i] == '{') k3++;
          if(in[i] == '}') k3--;
        }
        if(k5==0 && in[i] == '\'') k4^=1;
        if(k4==0 && in[i] == '"') k5^=1;
        if( (k4||k5) && in[i] == '\\') {i++; continue;}
        if(k1==0 && k2==0 && k3==0 && k4==0 && k5==0) return i;
      }
    } else if(in[p]==')' || in[p]==']' || in[p]=='}'){
      for(i=p;i>=0;i--){
        if(k4==0 && k5==0){
          if(in[i] == '(') k1++;
          if(in[i] == ')') k1--;
          if(in[i] == '[') k2++;
          if(in[i] == ']') k2--;
          if(in[i] == '{') k3++;
          if(in[i] == '}') k3--;
        }
        if(k5==0 && in[i] == '\'') k4^=1;
        if(k4==0 && in[i] == '"') k5^=1;
        if(k1==0 && k2==0 && k3==0 && k4==0 && k5==0) return i;
      }
    }
    return -1;
  }

  /* in = bb + argmin[hoge](piyo) + aa; */
  /* format = argmin[]() */
  /* return {"bb + ", "hoge", "piyo", " + aa;"} */
  /* if not found -> return {} */
  vector<string> findFunction(string &in, string format){
    int i, j, k, m, ok;
    int k4, k5;
    vector<string> fmt;
    vector<string> res;

    for(;;){
      string tmp;
      
      trim(format);
      if(format.size()==0) break;

      if(format.substr(0,2) == "()" || format.substr(0,2) == "[]" || format.substr(0,2) == "<>"){
        fmt.push_back(format.substr(0,2));
        format = format.substr(2);
        continue;
      }

      tmp = "";
      while(format.size() && (isalnum(format[0]) || format[0]=='_')){
        tmp += format[0];
        format = format.substr(1);
      }
      if(tmp == ""){
        tmp += format[0];
        format = format.substr(1);
      }
      fmt.push_back(tmp);
    }

    k4 = k5 = 0;
    rep(i,in.size()){
      if(k5 == 0 && in[i] == '"'){ k4 ^= 1; continue; }
      if(k4 == 0 && in[i] == '\''){ k5 ^= 1; continue; }
      if(k4 || k5){
        if(in[i] == '\\') i++;
        continue;
      }

      if(i && (isalnum(fmt[0][0]) || fmt[0][0]=='_') && (isalnum(in[i-1]) || in[i-1]=='_')) continue;

      res.clear();
      res.push_back( in.substr(0, i) );
      j = i;
      ok = 1;
      rep(k,fmt.size()){
        while(j < in.size() && isspace(in[j])) j++;
        if(j >= in.size()){ ok = 0; break; }
        
        if(fmt[k] == "()"){
          if(in[j]!='('){ ok = 0; break; }
          m = pairBracket(in, j);
          res.push_back( in.substr(j+1, m-j-1) );
          j = m + 1;
          while(j < in.size() && isspace(in[j])) j++;
          continue;
        }
        if(fmt[k] == "[]"){
          if(in[j]!='['){ ok = 0; break; }
          m = pairBracket(in, j);
          res.push_back( in.substr(j+1, m-j-1) );
          j = m + 1;
          while(j < in.size() && isspace(in[j])) j++;
          continue;
        }
        if(fmt[k] == "<>"){
          if(in[j]!='<'){ ok = 0; break; }
          m = pairBracketIq(in, j);
          res.push_back( in.substr(j+1, m-j-1) );
          j = m + 1;
          while(j < in.size() && isspace(in[j])) j++;
          continue;
        }
        
        if(in.substr(j, fmt[k].size()) == fmt[k]){
          j += fmt[k].size();
        } else {
          ok = 0;
          break;
        }
      }
      if(ok){
        res.push_back( in.substr(j) );
        return res;
      }
    }

    res.clear();
    return res;
  }


  int getExprLength_firstind(string tmp, int k){
    int i;
    
    tmp = tmp.substr(k);
    k = 0;

    for(;;){
      if(k==tmp.size()) break;
      if(isspace(tmp[k])){ k++; continue; }
      if(isOperator(tmp[k])){ k++; continue; }
      
      if(tmp[k]=='('){
        i = pairBracket(tmp, k);
        if(isValidVarType(tmp.substr(k+1, i-k-1), ' ')){
          k = i + 1;
          continue;
        }
      }
      
      break;
    }

    while(k < tmp.size() && (isalnum(tmp[k]) || tmp[k]=='.' || tmp[k]=='_')) k++;
    
    for(;;){
      while(k < tmp.size() && isspace(tmp[k])) k++;
      if(k==tmp.size()) break;
      
      if(tmp[k]=='('){
        k = pairBracket(tmp, k) + 1;
        tmp = tmp.substr(0, k);
        break;
      }

      if(tmp[k]=='['){
        k = pairBracket(tmp, k) + 1;
        continue;
      }

      tmp = tmp.substr(0, k);
      break;
    }

    return tmp.size();
  }

  int getExprLength_lastind(string tmp, int k){
    int i;

    tmp = tmp.substr(0, k+1);
    for(;;){
      if(k >= 6 && isspace(tmp[k]) && tmp.substr(k-6, 6) == "return") break;
      if(k >= 0 && isspace(tmp[k])){ k--; continue; }
      if(k >= 0 && (isalnum(tmp[k]) || tmp[k]=='.' || tmp[k]=='_')){ k--; continue; }
      if(tmp[k]==']' || tmp[k]==')'){
        k = pairBracket(tmp, k) - 1;
        if(k >= 5 && tmp.substr(k-5,6) == "return") break;
        continue;
      }
      break;
    }

    tmp = tmp.substr(k+1);
    return tmp.size();
  }


  string sentence_hatena_minmax_operator(string tmp){
    int i, j;
    int k1, k2, k3, k4, k5;
    int p, p1, p2, p3, p4, p5;
    string mode, arg1, arg2, bef, aft, exaft;

    for(;;){
      p1 = strpos_ns(tmp, (string)">?=");
      p2 = strpos_ns(tmp, (string)"<?=");
      p3 = strpos_ns(tmp, (string)"**=");
      p4 = strpos_ns(tmp, (string)"%%=");
      p5 = strpos_ns(tmp, (string)"/+=");
      if(p1<0 && p2<0 && p3<0 && p4<0 && p5<0) break;

      if(p1 > p2 && p1 > p3 && p1 > p4 && p1 > p5)      p = p1, mode = "chmax";
      else if(p2 > p1 && p2 > p3 && p2 > p3 && p2 > p5) p = p2, mode = "chmin";
      else if(p3 > p1 && p3 > p2 && p3 > p4 && p3 > p5) p = p3, mode = "pow_L";
      else if(p4 > p1 && p4 > p2 && p4 > p3 && p4 > p5) p = p4, mode = "moddw";
      else                                              p = p5, mode = "divup";
      ifun.doit.insert(mode);

      exaft = "";
      k1 = k2 = k3 = k4 = k5 = 0;
      for(i=p+3;i<tmp.size();i++){
        if((k4||k5) && tmp[i] == '\\'){ i++; continue; }
        if(k4==0 && k5==0){
          if(tmp[i]=='(') k1++;
          if(tmp[i]==')') k1--;
          if(tmp[i]=='[') k2++;
          if(tmp[i]==']') k2--;
          if(tmp[i]=='{') k3++;
          if(tmp[i]=='}') k3--;
        }
        if(k1 < 0 || k2 < 0 || k3 < 0){
          exaft = tmp.substr(i);
          tmp = tmp.substr(0,i);
          break;
        }
        if(k1 == 0 && k2 == 0 & k3 == 0 && k4 == 0 && k5 == 0 && (tmp[i] == ';' || tmp[i] == ',')){
          exaft = tmp.substr(i);
          tmp = tmp.substr(0,i);
          break;
        }
      }
      //fprintf(stderr, "1 [%s] [%s]\n", tmp.c_str(), exaft.c_str());
      
      k1 = k2 = k3 = k4 = k5 = 0;
      for(i=p-1;;i--){
        if(i >= 1 && k4 && tmp[i-1]=='/' && tmp[i]=='\''){ i-=2; continue; }
        if(i >= 0 && k5 == 0 && tmp[i]=='\'') k4 ^= 1;
        if(i >= 0 && k4 == 0 && tmp[i]=='"')  k5 ^= 1;
        if(i >= 0 && k4==0 && k5==0){
          if(tmp[i]==')') k1++;
          if(tmp[i]=='(') k1--;
          if(tmp[i]=='}') k2++;
          if(tmp[i]=='{') k2--;
          if(tmp[i]==']') k3++;
          if(tmp[i]=='[') k3--;
        }

        j = 0;
        if(k1 < 0 || k2 < 0 || k3 < 0) j = pairBracket(tmp, i);

        if(j==-1 || i==-1 || ((tmp[i] == '=' || tmp[i] == ',') && k1==0 && k2==0 && k3==0)){
          aft = "";
          bef = tmp.substr(0, i+1);
          arg1 = tmp.substr(i+1,p-(i+1));
          arg2 = tmp.substr(p+3);
          arg2 = arg2.substr(0, arg2.size());
        } else if(k1 < 0 || k2 < 0 || k3 < 0){
          bef = tmp.substr(0, i+1);
          arg1 = tmp.substr(i+1, p-(i+1));
          arg2 = tmp.substr(p+3, j-(p+3));
          aft = tmp.substr(j);
          aft = aft.substr(0, aft.size());
        } else {
          continue;
        }

        trim(bef);
        trim(aft);
        trim(arg1);
        trim(arg2);
        //fprintf(stderr, "2 [%s] [%s] [%s] [%s] [%s]\n", bef.c_str(), arg1.c_str(), arg2.c_str(), aft.c_str(), exaft.c_str());
        if(mode=="pow_L"){
          tmp = bef + "(" + arg1 + " = pow_L(" + arg1 + "," + arg2 + "))" + aft + exaft;
          ifun.doit.insert((string)"pow");
        } else if(mode=="moddw"){
          tmp = bef + "(" + arg1 + " = moddw_L(" + arg1 + "," + arg2 + "))" + aft + exaft;
          ifun.doit.insert((string)"moddw");
        } else if(mode=="divup"){
          tmp = bef + "(" + arg1 + " = divup_L(" + arg1 + "," + arg2 + "))" + aft + exaft;
          ifun.doit.insert((string)"divup");
        } else {
          tmp = bef + mode + "(" + arg1 + ", " + arg2 + ")" + aft + exaft;
        }
        break;
      }
    }

    return tmp;
  }

  string sentence_inequation(string tmp){
    int i, j, k, ls = 0, ok;
    int k1 = 0, k2 = 0, k3 = 0, k4 = 0, k5 = 0;
    string send, recv, tmp2, soec;
    vector<string> vs, vs2;

    rep(i,tmp.size()){

      k = -1;
      REP(j,1,tmp.size()-i){
        if(j>100) break;
        if(isValidVarType(tmp.substr(i,j), tmp[i+j])) k = j;
      }
      if(k > 0){
        i += k-1;
        continue;
      }

      if(tmp.substr(i,14) == "numeric_limits"){
        REP(j,14,tmp.size()-i) if(!isspace(tmp[i+j])) break;
        if(j < tmp.size()-i && tmp[i+j]=='<'){
          k = pairBracketIq2(tmp, i+j);
          //fprintf(stderr, "%d %d %d : %d\n", i, j, k, tmp.size());
          if(!(i < k)) fprintf(stderr, "error (Specify the type explicitly)\n");
          assert(i < k);
          i = k;
          continue;
        }
      }

      if(tmp.substr(i,16) == "remove_reference"){
        REP(j,16,tmp.size()-i) if(!isspace(tmp[i+j])) break;
        if(j < tmp.size()-i && tmp[i+j]=='<'){
          k = pairBracketIq2(tmp, i+j);
          //fprintf(stderr, "%d %d %d : %d\n", i, j, k, tmp.size());
          if(!(i < k)) fprintf(stderr, "error (Specify the type explicitly)\n");
          assert(i < k);
          i = k;
          continue;
        }
      }

      if(tmp[i]=='<'){
        int k6 = 0;
        REP(j,i,tmp.size()){
          if(tmp[j] == '<') k6++;
          if(tmp[j] == '>') k6--;
          if(k6==0) break;
        }
        if(j<tmp.size()){
          j++;
          tmp2 = tmp.substr(i+1, j-i-2);
          trim(tmp2);
          vs2 = split_p(tmp2, ',');
          ok = 0;
          rep(k,vs2.size()){
            trim(vs2[k]);
            if(isValidVarType(vs2[k], ' ')) ok++;
          }
          if(ok){
            i = j;
            continue;
          }
        }
      }
      
      if( (k4 || k5) && tmp[i]=='\\' ){ i++; continue; }
      if(k4==0 && k5==0){
        if(tmp[i]=='(') k1++;
        if(tmp[i]==')') k1--;
        if(tmp[i]=='[') k2++;
        if(tmp[i]==']') k2--;
        if(tmp[i]=='{') k3++;
        if(tmp[i]=='}') k3--;
      }
      if(k5 == 0 && tmp[i]=='\'') k4 ^= 1;
      if(k4 == 0 && tmp[i]=='"')  k5 ^= 1;
    
      if(k1 == 0 && k2 == 0 & k3 == 0 && k4 == 0 && k5 == 0){
        if(tmp.substr(i,2) == "->"){ i++; continue; }
        if(tmp.substr(i,2) == ">>"){ i++; continue; }
        if(tmp.substr(i,2) == "<<"){ i++; continue; }
        
        if(tmp.substr(i,2) == ">=" || tmp.substr(i,2) == "<=" || tmp.substr(i,2) == "==" || tmp.substr(i,2) == "!=" || tmp.substr(i,2) == "||" || tmp.substr(i,2) == "&&"){
          if(ls < i){
            vs.push_back( tmp.substr(ls, i-ls) );
            ls = i;
          }
          i++;
          vs.push_back( tmp.substr(ls, i+1-ls) );
          ls = i+1;
          continue;
        }
        if(tmp[i]=='<' || tmp[i]=='>' || tmp[i]==';' || tmp[i]==':' || tmp[i]=='?' || tmp[i]=='|' || tmp[i]=='&'){
          if(ls < i){
            vs.push_back( tmp.substr(ls, i-ls) );
            ls = i;
          }
          vs.push_back( tmp.substr(ls, i+1-ls) );
          ls = i+1;
          continue;
        }
      }
    }
    i = tmp.size();
    if(ls < i){
      vs.push_back( tmp.substr(ls, i-ls) );
      ls = i;
    }

    tmp = "";
    rep(i,vs.size()){
      if(vs[i]=="<" || vs[i]==">" || vs[i]=="<=" || vs[i]==">=" || vs[i]=="==" || vs[i]=="!="){
        if(i-2>=0 && (vs[i-2]=="<" || vs[i-2]==">" || vs[i-2]=="<=" || vs[i-2]==">=" || vs[i-2]=="==" || vs[i-2]=="!=")){
          if(!( vs[i-2]=="<" && vs[i]==">") && !(i-4 >= 0 && i+2 < vs.size() && vs[i-4] == "<" && vs[i-2] == ">" && vs[i] == "<" && vs[i+2] == ">")){
            tmp += " && " + vs[i-1];
          }
        }
      }
      tmp += vs[i];
    }
    
    
    rep(i,tmp.size()){
      if( (k4 || k5) && tmp[i]=='\\' ){ i++; continue; }
      if(k4==0 && k5==0){
        if(tmp[i]=='(') k1++;
        if(tmp[i]==')') k1--;
        if(tmp[i]=='[') k2++;
        if(tmp[i]==']') k2--;
        if(tmp[i]=='{') k3++;
        if(tmp[i]=='}') k3--;
      }
      if(k5 == 0 && tmp[i]=='\'') k4 ^= 1;
      if(k4 == 0 && tmp[i]=='"')  k5 ^= 1;

      if(k1==1 && k2==0 && k3==0 && k4==0 && k5==0 && tmp[i]=='('){
        j = pairBracket(tmp, i);
        send = tmp.substr(i+1, j-i-1);
        recv = sentence_inequation(send);
        tmp = tmp.substr(0, i+1) + recv + tmp.substr(j);
      }
      if(k1==0 && k2==1 && k3==0 && k4==0 && k5==0 && tmp[i]=='['){
        j = pairBracket(tmp, i);
        send = tmp.substr(i+1, j-i-1);
        recv = sentence_inequation(send);
        tmp = tmp.substr(0, i+1) + recv + tmp.substr(j);
      }
      if(k1==0 && k2==0 && k3==1 && k4==0 && k5==0 && tmp[i]=='{'){
        j = pairBracket(tmp, i);
        send = tmp.substr(i+1, j-i-1);
        recv = sentence_inequation(send);
        tmp = tmp.substr(0, i+1) + recv + tmp.substr(j);
      }
    }

    return tmp;
  }
  
  string sentence_times_operator(string tmp){
    int i, j, k, ok, t, doit = 0;
    int k1, k2;
    string now;
    pair<string,char> stchar;

    for(;;){
      int fg = 0;
      k1 = k2 = 0;
      now = "";
      rep(i,tmp.size()){
        ok = 1;

        if((k1 || k2) && tmp[i]=='\\'){ i++; continue; }
        if(k2==0 && tmp[i]=='"') k1 ^= 1;
        if(k1==0 && tmp[i]=='\'') k2 ^= 1;
        if(k1 || k2) continue;
        
        if(isalnum(tmp[i]) || tmp[i]=='_' || tmp[i]=='.') now += tmp[i];
        else                                              now = "";

        if(now[0]=='.' && !isdigit(now[1])) ok = 0;
        if(now == "0x" || now == "0X" || now == "0b" || now == "0B") ok = 0;
        if(now.substr(0,2) == "0x" || now.substr(0,2) == "0x"){
          REP(t,2,now.size()) if(!( ('0'<=now[t] && now[t]<='9') || ('a'<=now[t] && now[t]<='f') || ('A'<=now[t] && now[t]<='F') )) ok = 0;
        }
        if(now.substr(0,2) == "0b" || now.substr(0,2) == "0b"){
          REP(t,2,now.size()) if(!( now[t]=='0' || now[t]=='1' )) ok = 0;
        }

        if(!isdigit(now[0]) && now[0]!='.') ok = 0;
        if(now.size()==1 && now[0]=='0' && (tmp[i+1]=='b' || tmp[i+1]=='B') && (tmp[i+2]=='0' || tmp[i+2]=='1')) ok = 0;
        if(now.size()==1 && now[0]=='0' && (tmp[i+1]=='x' || tmp[i+1]=='X') && (('0'<=tmp[i+2] && tmp[i+2]<='9') || ('a'<=tmp[i+2] && tmp[i+2]<='f') || ('A'<=tmp[i+2] && tmp[i+2]<='F'))) ok = 0;
        if((tmp[i+1]=='d' || tmp[i+1]=='e' || tmp[i+1]=='D' || tmp[i+1]=='E') && (isdigit(tmp[i+2])||tmp[i+2]=='-')) ok = 0;
        if(tmp.substr(i+1,1)=="u" && (!isalnum(tmp[i+2]) && tmp[i+2]!='_')) ok = 0;
        if(tmp.substr(i+1,1)=="U" && (!isalnum(tmp[i+2]) && tmp[i+2]!='_')) ok = 0;
        if(tmp.substr(i+1,1)=="l" && (!isalnum(tmp[i+2]) && tmp[i+2]!='_')) ok = 0;
        if(tmp.substr(i+1,1)=="L" && (!isalnum(tmp[i+2]) && tmp[i+2]!='_')) ok = 0;
        if(tmp.substr(i+1,1)=="f" && (!isalnum(tmp[i+2]) && tmp[i+2]!='_')) ok = 0;
        if(tmp.substr(i+1,1)=="F" && (!isalnum(tmp[i+2]) && tmp[i+2]!='_')) ok = 0;
        if(tmp.substr(i+1,2)=="ll" && (!isalnum(tmp[i+3]) && tmp[i+3]!='_')) ok = 0;
        if(tmp.substr(i+1,2)=="LL" && (!isalnum(tmp[i+3]) && tmp[i+3]!='_')) ok = 0;
        if(tmp.substr(i+1,3)=="ull" && (!isalnum(tmp[i+4]) && tmp[i+4]!='_')) ok = 0;
        if(tmp.substr(i+1,3)=="ULL" && (!isalnum(tmp[i+4]) && tmp[i+4]!='_')) ok = 0;

        j = i + 1;
        while(j < tmp.size() && isspace(tmp[j])) j++;
        if(j==i+1 && (isdigit(tmp[j]) || tmp[j]=='.')) ok = 0;
        if(!(isalnum(tmp[j]) || tmp[j]=='_' || tmp[j]=='.' || tmp[j]=='(')) ok = 0;

        if(ok){
          tmp = tmp.substr(0,i+1) + "*" + tmp.substr(i+1);
//          int len = getExprLength_firstind(tmp, j);
//          tmp = tmp.substr(0,i+1-now.size()) + "(" + tmp.substr(i+1-now.size(),now.size()) + "*" + tmp.substr(i+1,len) + ")" + tmp.substr(i+1+len);
          fg = 1;
          doit++;
          break;
        }
      }
      if(!fg) break;
    }

    if(doit) code_replace(tmp);
    return tmp;
  }
  
  string sentence_pow_operator(string tmp){
    int i, j, k, st = 0;
    int bi, bj, ei, ej;
    int len_i, len_j;
    string bef, t1, t2, aft;
    
    while(strpos_ns(tmp, (string)"**", st) >= 0){
      k = strpos_ns(tmp, (string)"**", st);
      st = k+1;
      
      i = k-1;
      while(i >= 0 && isspace(tmp[i])) i--;
      if(i < 0 || isOperator(tmp[i])) continue;
      for(j=i;j>=0&&j>=i-100;j--) if(isValidVarType(tmp.substr(j,i-j+1),tmp[i+1])) break;
      if(j>=0 && j>=i-100) continue;

      j = k+2;
      while(j < tmp.size() && isspace(tmp[j])) j++;
      if(j == tmp.size()) continue;

      len_i = getExprLength_lastind(tmp, i);
      len_j = getExprLength_firstind(tmp, j);

      bi = i - len_i + 1;
      ei = i;

      bj = j;
      ej = j + len_j - 1;

      bef = tmp.substr(0, bi);
      t1 = tmp.substr(bi, ei-bi+1);
      t2 = tmp.substr(bj, ej-bj+1);
      aft = tmp.substr(ej+1);

      trim(t1);
      trim(t2);

      if(t2 == "2"){
        tmp = bef + "(pow2_L(" + t1 + "))" + aft;
        ifun.doit.insert((string)"pow2");
      } else if(t2 == "3"){
        tmp = bef + "(pow3_L(" + t1 + "))" + aft;
        ifun.doit.insert((string)"pow3");
      } else if(t2 == "4"){
        tmp = bef + "(pow4_L(" + t1 + "))" + aft;
        ifun.doit.insert((string)"pow4");
      } else {
        tmp = bef + "(pow_L(" + t1 + "," + t2 + "))" + aft;
        ifun.doit.insert((string)"pow");
      }
      st = 0;
    }
    
    return tmp;
  }

  string sentence_div_operator(string tmp){
    int i, j, k, st = 0, loop;
    int bi, bj, ei, ej;
    int len_i, len_j;
    string bef, t1, t2, aft, op;

    rep(loop,2){
      if(loop==0) op = "/+";
      if(loop==1) op = "%%";
      
      while(strpos_ns(tmp, op, st) >= 0){
        k = strpos_ns(tmp, op, st);
        st = k+1;
        
        i = k-1;
        while(i >= 0 && isspace(tmp[i])) i--;
        if(i < 0 || isOperator(tmp[i])) continue;
        
        j = k+2;
        while(j < tmp.size() && isspace(tmp[j])) j++;
        if(j == tmp.size()) continue;
        
        len_i = getExprLength_lastind(tmp, i);
        len_j = getExprLength_firstind(tmp, j);
        
        bi = i - len_i + 1;
        ei = i;
        
        bj = j;
        ej = j + len_j - 1;
        
        bef = tmp.substr(0, bi);
        t1 = tmp.substr(bi, ei-bi+1);
        t2 = tmp.substr(bj, ej-bj+1);
        aft = tmp.substr(ej+1);
        
        trim(t1);
        trim(t2);

        if(loop==0){
          tmp = bef + "(divup_L(" + t1 + "," + t2 + "))" + aft;
          ifun.doit.insert((string)"divup");
        }
        if(loop==1){
          tmp = bef + "(moddw_L(" + t1 + "," + t2 + "))" + aft;
          ifun.doit.insert((string)"moddw");
        }
        st = 0;
      }
    }
    
    return tmp;
  }


  // tmp = hoge -> return (*((hoge*)NULL))
  string get_type_with_decltype_basic_type(string tmp){
    return "(*((" + tmp + "*)NULL))";
  }

  string get_type_with_decltype(string tmp){
    vector<string> vs;
    string str;

    vs = findFunction(tmp, "sum[][]()");
    if(vs.size()==0) vs = findFunction(tmp, "mul[][]()");
    if(vs.size()==0) vs = findFunction(tmp, "gcd[][]()");
    if(vs.size()==0) vs = findFunction(tmp, "lcm[][]()");
    if(vs.size()==0) vs = findFunction(tmp, "min[][]()");
    if(vs.size()==0) vs = findFunction(tmp, "max[][]()");
    if(vs.size()==0) vs = findFunction(tmp, "XOR[][]()");
    if(vs.size()==0) vs = findFunction(tmp, "OR[][]()");
    if(vs.size()==0) vs = findFunction(tmp, "AND[][]()");
    if(vs.size() != 0) return get_type_with_decltype(vs[0] + get_type_with_decltype_basic_type(vs[1]) + vs[4]);

    if(vs.size()==0) vs = findFunction(tmp, "argmin[][]()");
    if(vs.size()==0) vs = findFunction(tmp, "argmax[][]()");
    if(vs.size()==0) vs = findFunction(tmp, "argminL[][]()");
    if(vs.size()==0) vs = findFunction(tmp, "argmaxL[][]()");
    if(vs.size() != 0) return get_type_with_decltype(vs[0] + get_type_with_decltype_basic_type("int") + vs[4]);

    if(vs.size()==0) vs = findFunction(tmp, "argmin[]()");
    if(vs.size()==0) vs = findFunction(tmp, "argmax[]()");
    if(vs.size()==0) vs = findFunction(tmp, "argminL[]()");
    if(vs.size()==0) vs = findFunction(tmp, "argmaxL[]()");
    if(vs.size() != 0) return get_type_with_decltype(vs[0] + get_type_with_decltype_basic_type("int") + vs[3]);

    vs = findFunction(tmp, "sum[]()");
    if(vs.size()==0) vs = findFunction(tmp, "mul[]()");
    if(vs.size()==0) vs = findFunction(tmp, "gcd[]()");
    if(vs.size()==0) vs = findFunction(tmp, "lcm[]()");
    if(vs.size()==0) vs = findFunction(tmp, "min[]()");
    if(vs.size()==0) vs = findFunction(tmp, "max[]()");
    if(vs.size()==0) vs = findFunction(tmp, "XOR[]()");
    if(vs.size()==0) vs = findFunction(tmp, "OR[]()");
    if(vs.size()==0) vs = findFunction(tmp, "AND[]()");
    if(vs.size() != 0){
      int i, j, k;
      string var = vs[1];
      string eq = vs[2];
      vector<string> tmpvs = split_p(var, ',');
      vector<string> tmpvs2, loopvar, loopvar_type;
      
      rep(i,tmpvs.size()){
        string cur = tmpvs[i];
        string cur_name, cur_L, cur_R;
        vector<string> cur_names;
        tmpvs2 = split_p(cur, '=');
        if(tmpvs2.size()==2){
          int pos = strpos_ns(tmpvs2[1], "---");
          assert(pos >= 0 && "min|max|sum|mul[][]()");
          cur_name = tmpvs2[0];
          cur_L = tmpvs2[1].substr(0,pos);
          cur_R = tmpvs2[1].substr(pos+3);
        } else {
          cur_name = cur;
          cur_L = tmpvs[++i];
          cur_R = tmpvs[++i];
          cur_R = "(" + cur_R + ")-1";
        }
        trim(cur_name);
        trim(cur_L);
        trim(cur_R);

        if(cur_name[0]=='('){
          cur_name = cur_name.substr(1, cur_name.size()-2);
          cur_names = split_p(cur_name, ',');
          rep(i,cur_names.size()) trim(cur_names[i]);
        } else {
          cur_names.push_back(cur_name);
        }

        rep(j,cur_names.size()){
          cur_name = cur_names[j];
          for(k=cur_name.size()-1;k>=0;k--) if(cur_name[k]==' ') break;
          tmpvs2.clear();
          if(k < 0){
            tmpvs2.push_back("int");
            tmpvs2.push_back(cur_name);
          } else {
            tmpvs2.push_back(cur_name.substr(0,k));
            tmpvs2.push_back(cur_name.substr(k));
            trim(tmpvs2[0]);
            trim(tmpvs2[1]);
          }
          loopvar.push_back(tmpvs2[1]);
          loopvar_type.push_back(tmpvs2[0]);
        }
      }
      rep(i,loopvar.size()) replaceAll_ns_t(eq, loopvar[i], get_type_with_decltype_basic_type(loopvar_type[i]));

      return get_type_with_decltype(vs[0] + eq + vs[3]);
    }

    vs = findFunction(tmp, "sum()");
    if(vs.size()==0) vs = findFunction(tmp, "SUM()");
    if(vs.size()==0) vs = findFunction(tmp, "mul()");
    if(vs.size()==0) vs = findFunction(tmp, "MUL()");
    if(vs.size()==0) vs = findFunction(tmp, "gcd()");
    if(vs.size()==0) vs = findFunction(tmp, "GCD()");
    if(vs.size()==0) vs = findFunction(tmp, "lcm()");
    if(vs.size()==0) vs = findFunction(tmp, "LCM()");
    if(vs.size()==0) vs = findFunction(tmp, "min()");
    if(vs.size()==0) vs = findFunction(tmp, "MIN()");
    if(vs.size()==0) vs = findFunction(tmp, "max()");
    if(vs.size()==0) vs = findFunction(tmp, "MAX()");
    if(vs.size()==0) vs = findFunction(tmp, "XOR()");
    if(vs.size()==0) vs = findFunction(tmp, "OR()");
    if(vs.size()==0) vs = findFunction(tmp, "AND()");
    if(vs.size() != 0){
      string eq;
      vector<string> arr, barr;
      arr = split_p(vs[1], ',');
      trim(arr[0]);
      barr = rd_wt_array(arr[0]);
      if(barr.size()==4){
        eq = get_type_with_decltype_basic_type(getEquationType(barr[1]));
      } else {
        eq = arr[0];
      }
      return get_type_with_decltype(vs[0] + eq + vs[2]);
    }
    
    return tmp;
  }

  string sentence_minmax_function_one(string tmp, int blocked = 1){
    int i, j, k, mode;
    string tp, cd, var, eqn, ieq, resvar, flagvar, indvar, candvar;
    vector<string> chkvs, vs, tmpvs, tmpvs2, tmpvs3;
    vector<string> loopvar_tmpname, loopvar, loopvar_type, loopL, loopR;
    blocked = 0;

    for(;;){
      vs.clear();

      if(vs.size()==0){ vs = findFunction(tmp, "argmin()"); tp = "argmin"; }
      if(vs.size()==0){ vs = findFunction(tmp, "argmax()"); tp = "argmax"; }
      if(vs.size()==0){ vs = findFunction(tmp, "argminL()"); tp = "argminL"; }
      if(vs.size()==0){ vs = findFunction(tmp, "argmaxL()"); tp = "argmaxL"; }
      if(vs.size()==3){
        var = getUnusedVarName();
        trim(vs[1]);
        j = vs[1].size() - 1;
        assert(vs[1][j]==')');
        i = pairBracket(vs[1], j);
        vs[1] = tp+"["+var+"=0---("+vs[1].substr(i+1,j-i-1)+")-1]("+vs[1].substr(0,i)+"["+var+"])";
        tmp = vs[0] + vs[1] + vs[2];
        break;
      }


      chkvs = findFunction(tmp, "sum[][]()");
      if(chkvs.size() && (vs.size()==0 || vs[0].size() > chkvs[0].size())){
        vs = chkvs;
        ieq = "+", mode = 0;
      }

      chkvs = findFunction(tmp, "sum[]()");
      if(chkvs.size() && (vs.size()==0 || vs[0].size() > chkvs[0].size())){
        vs = chkvs;
        ieq = "+", mode = 0;
        vs.push_back(""); vs[4] = vs[3]; vs[3] = vs[2]; vs[2] = vs[1]; vs[1] = "auto";
      }

      chkvs = findFunction(tmp, "mul[][]()");
      if(chkvs.size() && (vs.size()==0 || vs[0].size() > chkvs[0].size())){
        vs = chkvs;
        ieq = "*", mode = 0;
      }

      chkvs = findFunction(tmp, "mul[]()");
      if(chkvs.size() && (vs.size()==0 || vs[0].size() > chkvs[0].size())){
        vs = chkvs;
        ieq = "*", mode = 0;
        vs.push_back(""); vs[4] = vs[3]; vs[3] = vs[2]; vs[2] = vs[1]; vs[1] = "auto";
      }

      chkvs = findFunction(tmp, "gcd[][]()");
      if(chkvs.size() && (vs.size()==0 || vs[0].size() > chkvs[0].size())){
        vs = chkvs;
        ieq = "gcd", mode = 0;
      }

      chkvs = findFunction(tmp, "gcd[]()");
      if(chkvs.size() && (vs.size()==0 || vs[0].size() > chkvs[0].size())){
        vs = chkvs;
        ieq = "gcd", mode = 0;
        vs.push_back(""); vs[4] = vs[3]; vs[3] = vs[2]; vs[2] = vs[1]; vs[1] = "auto";
      }

      chkvs = findFunction(tmp, "lcm[][]()");
      if(chkvs.size() && (vs.size()==0 || vs[0].size() > chkvs[0].size())){
        vs = chkvs;
        ieq = "lcm", mode = 0;
      }

      chkvs = findFunction(tmp, "lcm[]()");
      if(chkvs.size() && (vs.size()==0 || vs[0].size() > chkvs[0].size())){
        vs = chkvs;
        ieq = "lcm", mode = 0;
        vs.push_back(""); vs[4] = vs[3]; vs[3] = vs[2]; vs[2] = vs[1]; vs[1] = "auto";
      }

      chkvs = findFunction(tmp, "min[][]()");
      if(chkvs.size() && (vs.size()==0 || vs[0].size() > chkvs[0].size())){
        vs = chkvs;
        ieq = ">", mode = 0;
      }

      chkvs = findFunction(tmp, "min[]()");
      if(chkvs.size() && (vs.size()==0 || vs[0].size() > chkvs[0].size())){
        vs = chkvs;
        ieq = ">", mode = 0;
        vs.push_back(""); vs[4] = vs[3]; vs[3] = vs[2]; vs[2] = vs[1]; vs[1] = "auto";
      }

      chkvs = findFunction(tmp, "max[][]()");
      if(chkvs.size() && (vs.size()==0 || vs[0].size() > chkvs[0].size())){
        vs = chkvs;
        ieq = "<", mode = 0;
      }

      chkvs = findFunction(tmp, "max[]()");
      if(chkvs.size() && (vs.size()==0 || vs[0].size() > chkvs[0].size())){
        vs = chkvs;
        ieq = "<", mode = 0;
        vs.push_back(""); vs[4] = vs[3]; vs[3] = vs[2]; vs[2] = vs[1]; vs[1] = "auto";
      }

      chkvs = findFunction(tmp, "argmin[][]()");
      if(chkvs.size() && (vs.size()==0 || vs[0].size() > chkvs[0].size())){
        vs = chkvs;
        ieq = ">", mode = 1;
      }

      chkvs = findFunction(tmp, "argmin[]()");
      if(chkvs.size() && (vs.size()==0 || vs[0].size() > chkvs[0].size())){
        vs = chkvs;
        ieq = ">", mode = 1;
        vs.push_back(""); vs[4] = vs[3]; vs[3] = vs[2]; vs[2] = vs[1]; vs[1] = "auto";
      }

      chkvs = findFunction(tmp, "argmax[][]()");
      if(chkvs.size() && (vs.size()==0 || vs[0].size() > chkvs[0].size())){
        vs = chkvs;
        ieq = "<", mode = 1;
      }

      chkvs = findFunction(tmp, "argmax[]()");
      if(chkvs.size() && (vs.size()==0 || vs[0].size() > chkvs[0].size())){
        vs = chkvs;
        ieq = "<", mode = 1;
        vs.push_back(""); vs[4] = vs[3]; vs[3] = vs[2]; vs[2] = vs[1]; vs[1] = "auto";
      }

      chkvs = findFunction(tmp, "argminL[][]()");
      if(chkvs.size() && (vs.size()==0 || vs[0].size() > chkvs[0].size())){
        vs = chkvs;
        ieq = ">=", mode = 1;
      }

      chkvs = findFunction(tmp, "argminL[]()");
      if(chkvs.size() && (vs.size()==0 || vs[0].size() > chkvs[0].size())){
        vs = chkvs;
        ieq = ">=", mode = 1;
        vs.push_back(""); vs[4] = vs[3]; vs[3] = vs[2]; vs[2] = vs[1]; vs[1] = "auto";
      }

      chkvs = findFunction(tmp, "argmaxL[][]()");
      if(chkvs.size() && (vs.size()==0 || vs[0].size() > chkvs[0].size())){
        vs = chkvs;
        ieq = "<=", mode = 1;
      }

      chkvs = findFunction(tmp, "argmaxL[]()");
      if(chkvs.size() && (vs.size()==0 || vs[0].size() > chkvs[0].size())){
        vs = chkvs;
        ieq = "<=", mode = 1;
        vs.push_back(""); vs[4] = vs[3]; vs[3] = vs[2]; vs[2] = vs[1]; vs[1] = "auto";
      }

      chkvs = findFunction(tmp, "XOR[][]()");
      if(chkvs.size() && (vs.size()==0 || vs[0].size() > chkvs[0].size())){
        vs = chkvs;
        ieq = "^", mode = 0;
      }

      chkvs = findFunction(tmp, "XOR[]()");
      if(chkvs.size() && (vs.size()==0 || vs[0].size() > chkvs[0].size())){
        vs = chkvs;
        ieq = "^", mode = 0;
        vs.push_back(""); vs[4] = vs[3]; vs[3] = vs[2]; vs[2] = vs[1]; vs[1] = "auto";
      }

      chkvs = findFunction(tmp, "OR[][]()");
      if(chkvs.size() && (vs.size()==0 || vs[0].size() > chkvs[0].size())){
        vs = chkvs;
        ieq = "|", mode = 0;
      }

      chkvs = findFunction(tmp, "OR[]()");
      if(chkvs.size() && (vs.size()==0 || vs[0].size() > chkvs[0].size())){
        vs = chkvs;
        ieq = "|", mode = 0;
        vs.push_back(""); vs[4] = vs[3]; vs[3] = vs[2]; vs[2] = vs[1]; vs[1] = "auto";
      }

      chkvs = findFunction(tmp, "AND[][]()");
      if(chkvs.size() && (vs.size()==0 || vs[0].size() > chkvs[0].size())){
        vs = chkvs;
        ieq = "&", mode = 0;
      }

      chkvs = findFunction(tmp, "AND[]()");
      if(chkvs.size() && (vs.size()==0 || vs[0].size() > chkvs[0].size())){
        vs = chkvs;
        ieq = "&", mode = 0;
        vs.push_back(""); vs[4] = vs[3]; vs[3] = vs[2]; vs[2] = vs[1]; vs[1] = "auto";
      }

      if(vs.size()==0) break;
      // vs[0] 
      // vs[1] type
      // vs[2] range
      // vs[3] eq
      // vs[4]
      //fprintf(stderr, "[%s] [%s] [%s] [%s] [%s]\n", vs[0].c_str(), vs[1].c_str(), vs[2].c_str(), vs[3].c_str(), vs[4].c_str());
      vs[0] = vs[0] + "(";
      vs[4] = ")" + vs[4];

      string count_condition, init_value;
      for(;;){
        int fg = 0;
        for(i=vs[2].size()-1; i>=0; i--){
          if(vs[2][i] == '@'){
            count_condition = vs[2].substr(i+1);
            vs[2] = vs[2].substr(0,i);
            fg = 1;
            break;
          }
          if(vs[2][i] == '$'){
            init_value = vs[2].substr(i+1);
            vs[2] = vs[2].substr(0,i);
            fg = 1;
            break;
          }
        }
        if(!fg) break;
      }
      trim(vs[2]);
      trim(count_condition);
      trim(init_value);

      tmpvs = split_p(vs[2], ',');
      rep(i,tmpvs.size()) trim(tmpvs[i]);
      rep(i,tmpvs.size()){
        string cur = tmpvs[i];
        string cur_name, cur_L, cur_R;
        vector<string> cur_names;
        tmpvs2 = split_p(cur, '=');
        if(tmpvs2.size()==2){
          int pos = strpos_ns(tmpvs2[1], "---");
          assert(pos >= 0 && "min|max|sum|mul[][]()");
          cur_name = tmpvs2[0];
          cur_L = tmpvs2[1].substr(0,pos);
          cur_R = tmpvs2[1].substr(pos+3);
        } else {
          cur_name = cur;
          cur_L = tmpvs[++i];
          cur_R = tmpvs[++i];
          cur_R = "(" + cur_R + ")-1";
        }
        trim(cur_name);
        trim(cur_L);
        trim(cur_R);

        if(cur_name[0]=='('){
          cur_name = cur_name.substr(1, cur_name.size()-2);
          cur_names = split_p(cur_name, ',');
          rep(i,cur_names.size()) trim(cur_names[i]);
        } else {
          cur_names.push_back(cur_name);
        }

        rep(j,cur_names.size()){
          cur_name = cur_names[j];
          for(k=cur_name.size()-1;k>=0;k--) if(cur_name[k]==' ') break;
          tmpvs2.clear();
          if(k < 0){
            tmpvs2.push_back("int");
            tmpvs2.push_back(cur_name);
          } else {
            tmpvs2.push_back(cur_name.substr(0,k));
            tmpvs2.push_back(cur_name.substr(k));
            trim(tmpvs2[0]);
            trim(tmpvs2[1]);
          }
          loopvar_tmpname.push_back( getUnusedVarName() );
          loopvar.push_back(tmpvs2[1]);
          loopvar_type.push_back(tmpvs2[0]);
          loopL.push_back(cur_L);
          loopR.push_back(cur_R);
        }
      }
      
      resvar = getUnusedVarName();
      flagvar = getUnusedVarName();
      indvar = getUnusedVarName();
      candvar = getUnusedVarName();

      eqn = vs[3];
      rep(i,loopvar.size()) replaceAll_ns_t(eqn, loopvar[i], loopvar_tmpname[i]);
      rep(i,loopvar.size()) replaceAll_ns_t(count_condition, loopvar[i], loopvar_tmpname[i]);

      if(vs[1] == "auto") vs[1] = "remove_reference<decltype(" + get_type_with_decltype(eqn) + ")>::type";

      if(init_value == ""){
        if(mode || ieq == "+" || ieq == "gcd" || ieq == "^" || ieq == "|" || ieq == "&"){
          init_value = "0";
        } else if(ieq == "*" || ieq == "lcm"){
          init_value = "1";
        } else if(ieq == ">"){
          init_value = "numeric_limits<" + vs[1] + ">::max()";
        } else {
          init_value = "numeric_limits<" + vs[1] + ">::lowest()";
        }
      }

      if(blocked) cd += "{";
      rep(i,loopvar.size()) cd += loopvar_type[i] + " " + loopvar_tmpname[i] + ";";
      cd += vs[1] + " " + resvar + ";";
      cd += "int " + flagvar + " = 0;";
      if(mode) cd += "int " + indvar + ";";

      cd += "if(";
      {
        rep(i,loopvar.size()){
          if(i) cd += "||";
          cd += "(" + loopL[i] + ") > (" + loopR[i] + ")";
        }
      }
      cd += "){";
      {
        cd += resvar + " = " + init_value + ";";
      }
      cd += "} else {";
      {
        rep(i,loopvar.size()) cd += "for(" + loopvar_tmpname[i] + " = " + loopL[i] + "; " + loopvar_tmpname[i] + " <= " + loopR[i] + "; " + loopvar_tmpname[i] + "++)";
        if(count_condition != "") cd += "if(" + count_condition + ")";
        cd += "{";
        
        cd += "if(" + flagvar + " == 0){";
        cd += resvar + " = " + eqn + ";";
        if(mode) cd += indvar + " = " + loopvar_tmpname[0] + ";";
        cd += flagvar + " = 1;";
        cd += "continue;";
        cd += "}";

        if(ieq == "+" || ieq == "*" || ieq == "^" || ieq == "|" || ieq == "&"){
          cd += resvar + " " + ieq + "= " + eqn + ";";
        } else if(ieq == "gcd" || ieq == "lcm"){
          cd += resvar + " = " + ieq + "(" + resvar + ", " + eqn + ");";
        } else {
          cd += "const auto " + candvar + " = " + eqn + ";";
          cd += "if(" + resvar + " " + ieq + " " + candvar + "){";
          cd += resvar + " = " + candvar + ";";
          if(mode) cd += indvar + " = " + loopvar_tmpname[0] + ";";
          cd += "}";
        }

        cd += "}";
        
        cd += "if(" + flagvar + "==0){";
        cd += resvar + " = " + init_value + ";";
        cd += "}";
      }
      cd += "}";

      if(blocked){
        fprintf(stderr,"hoge6-1\n");
        if(mode==0) cd += vs[0] + resvar + vs[4];
        if(mode==1) cd += vs[0] + indvar + vs[4];
        
        cd += "}";
        insert(cd, str.size());
        return "";
      } else {
        insert(cd, str.size());
        if(mode==0) tmp = vs[0] + resvar + vs[4];
        if(mode==1) tmp = vs[0] + indvar + vs[4];
      }
      break;
    }
    return tmp;
  }


  string sentence_sortE(string tmp){
    int i, j, n;
    string stc, op;
    vector<string> vs, e;

    if(vs.size() == 0){
      vs = findFunction(tmp, "sortE()");
      if(vs.size() > 0) op = ">";
    }
    if(vs.size() == 0){
      vs = findFunction(tmp, "rsortE()");
      if(vs.size() > 0) op = "<";
    }
    if(vs.size()==0) return tmp;

    e = split_p(vs[1], ',');
    n = e.size();
    rep(i,n) trim(e[i]);
    for(j=n-1;j;j--){
      rep(i,j){
        stc += "if(" + e[i] + " " + op + " " + e[i+1] + ") swap(" + e[i] + ", " + e[i+1] + ");";
      }
    }

    stc = vs[0] + stc + vs[2];
    insert(stc, str.size());
    return "";
  }

  string sentence_reader(string tmp){
    int i, j, k, m, ignore_error = 0;
    int inc;
    pair<string, char> stchar;
    vector<string> vtmp, vvtmp;
    string at, stmp;

    stchar = nextToken(tmp);
    if( (stchar.first == "rd" || stchar.first == "reader" || stchar.first == "reader_ignore_error") && stchar.second == '(' ){
      ignore_error = 0;
      if(stchar.first == "reader_ignore_error") ignore_error = 1;
      trim_until(tmp, '(', ')');
      tmp = tmp.substr(1,tmp.size()-2);
      vtmp = split_p(tmp, ',');
      rep(k,vtmp.size()){
        trim(vtmp[k]);

        string stc, trail;
        string strtmp;
        vector<string> vstrtmp;

        strtmp = vtmp[k];
        trim(strtmp);
        vstrtmp = rd_wt_array(strtmp);
        if(vstrtmp.size() > 1){
          //fprintf(stderr, "rdwt [%s] [%s] [%s] [%s]\n", vstrtmp[0].c_str(), vstrtmp[1].c_str(), vstrtmp[2].c_str(), vstrtmp[3].c_str());

          vector<string> vn, looptm;
          string var;

          looptm = split_p(vstrtmp[2], ',');
          rep(m,looptm.size()) vn.push_back(getUnusedVarName());

          vector<string> vars;
          if(vstrtmp[1][0]=='('){
            vstrtmp[1] = vstrtmp[1].substr(1, vstrtmp[1].size()-2);
            vars = split_p(vstrtmp[1], ',');
          } else {
            vars.push_back(vstrtmp[1]);
          }
          var = "";
          rep(j,vars.size()){
            if(j) var += ", ";
            trim(vars[j]);
            trail = "";
            while(vars[j].size() >= 2 && (vars[j].substr(vars[j].size()-2)=="++" || vars[j].substr(vars[j].size()-2)=="--")){
              trail = trail + vars[j].substr(vars[j].size()-2);
              vars[j] = vars[j].substr(0, vars[j].size()-2);
              trim(vars[j]);
            }
            vvtmp = split_p(vars[j], '@');
            stmp = "";
            rep(i,vvtmp.size()){
              if(i) stmp += (string)"@";
              stmp += vvtmp[i];
              rep(m,looptm.size()) stmp += (string)"[" + vn[m] + "]";
              stmp += trail;
            }
            //fprintf(stderr, "%s\n", stmp.c_str());
            var += vstrtmp[0] + stmp + vstrtmp[3];

            string vtype = getVarType(vvtmp[0]);
            if(vtype.substr(0,5) == "Arr1d"){
              string hoge_tmp = vvtmp[0] + ".malloc(" + looptm[0] + ");";
              insert(hoge_tmp, str.size());
            }
            if(vtype.substr(0,5) == "Arr2d"){
              string hoge_tmp = vvtmp[0] + ".malloc(" + looptm[0] + "," + looptm[1] + ");";
              insert(hoge_tmp, str.size());
            }
          }
          stc = (string)"{";
          rep(m,looptm.size()) stc += (string)"int " + vn[m] + ";";
          rep(m,looptm.size()) stc += (string)" rep(" + vn[m] + ", " + looptm[m] + ")";
          stc += (string)"rd(" + var + ");";
          stc += (string)"}";
          insert(stc, str.size());
        } else {
          vvtmp = split_p(vstrtmp[0], '@');
          at = "";
          inc = 0;
          if(vvtmp.size()>=2){
            vstrtmp[0] = vvtmp[0];
            at = vvtmp[1];
            trim(vstrtmp[0]);
            trim(at);
          }
          for(;;){
            if(vstrtmp[0].length() >= 2 && vstrtmp[0].substr(vstrtmp[0].length()-2) == "++"){
              inc++;
              vstrtmp[0] = vstrtmp[0].substr(0, vstrtmp[0].length()-2);
              trim(vstrtmp[0]);
              continue;
            }
            if(vstrtmp[0].length() >= 2 && vstrtmp[0].substr(vstrtmp[0].length()-2) == "--"){
              inc--;
              vstrtmp[0] = vstrtmp[0].substr(0, vstrtmp[0].length()-2);
              trim(vstrtmp[0]);
              continue;
            }
            break;
          }

          int fg = 0;
          string etype = getElementalyVarType(vstrtmp[0]);
          getEquationType(vstrtmp[0], &fg);
          //fprintf(stderr, "type estimation: [%s] [%s] [fg=%d]\n",etype.c_str(), vstrtmp[0].c_str(), fg);

          if(fg) ifun.doit.insert((string)"reader_all");
          if(etype=="int" || etype=="signed" || etype=="signed int"){
            ifun.doit.insert((string)"reader_int");
          } else if(etype=="long long" || etype=="signed long long" || etype=="ll"){
            ifun.doit.insert((string)"reader_ll");
          } else if(etype=="unsigned" || etype=="unsigned int" || etype=="uint"){
            ifun.doit.insert((string)"reader_unsigned");
          } else if(etype=="unsigned long long" || etype=="ull"){
            ifun.doit.insert((string)"reader_ull");
          } else if(etype=="int8_t" || etype=="int16_t" || etype=="int32_t" || etype=="int64_t" || etype=="__int8_t" || etype=="__int16_t" || etype=="__int32_t" || etype=="__int64_t" || etype=="__int128_t"){
            ifun.doit.insert((string)"reader_int");
            ifun.doit.insert((string)"reader_ll");
            ifun.doit.insert((string)"reader_int128");
          } else if(etype=="uint8_t" || etype=="uint16_t" || etype=="uint32_t" || etype=="uint64_t" || etype=="__uint8_t" || etype=="__uint16_t" || etype=="__uint32_t" || etype=="__uint64_t" || etype=="__uint128_t"){
            ifun.doit.insert((string)"reader_unsigned");
            ifun.doit.insert((string)"reader_ull");
            ifun.doit.insert((string)"reader_uint128");
          } else if(etype=="Modint"){
            ifun.doit.insert((string)"reader_Modint");
          } else if(etype=="Mint"){
            ifun.doit.insert((string)"reader_Mint");
          } else if(etype=="modint"){
            ifun.doit.insert((string)"reader_modint");
          } else if(etype=="mint"){
            ifun.doit.insert((string)"reader_mint");
          } else if(etype=="double"){
            ifun.doit.insert((string)"reader_double");
          } else if(etype=="char"){
            ifun.doit.insert((string)"reader_char");
            ifun.doit.insert((string)"reader_char_array");
          } else if(etype=="string"){
            ifun.doit.insert((string)"reader_string");
          } else {
            if(ignore_error || etype.substr(0,8) == "Point2d<"){
            } else {
              fprintf(stderr, "warning: unknown type [%s] for rd (reader) : %s\n", etype.c_str(), vtmp[k].c_str());
            }
            //assert(0);
          }
          
          if(at!="") stc = at + " = rd(" + vstrtmp[0] + ");";
          else       stc = (string)"rd(" + vstrtmp[0] + ");";
          if(inc){
            char tmpbuf[10];
            sprintf(tmpbuf, "%d", inc);
            stc += vstrtmp[0] + " += (" + tmpbuf + ");";
          }
          
          str.push_back(stc);
          nxt.push_back(-1);
          strtype.push_back((string)"sentence");
        }
      }
      return "";
    }
    
    return tmp;
  }

  string sentence_writer(string tmp){
    int i, j, k, extend = 1;
    string mode;
    pair<string, char> stchar;
    vector<string> vtmp;

    stchar = nextToken(tmp);
    if( (stchar.first == "wt"   || stchar.first == "writer"  ) && (stchar.second == '[' || stchar.second == '(') ) mode = "wt";
    if( (stchar.first == "wtSp" || stchar.first == "writerSp") && (stchar.second == '[' || stchar.second == '(') ) mode = "wtSp";
    if( (stchar.first == "wtLn" || stchar.first == "writerLn") && (stchar.second == '[' || stchar.second == '(') ) mode = "wtLn";
    if( (stchar.first == "wtN"  || stchar.first == "writerN" ) && (stchar.second == '[' || stchar.second == '(') ) mode = "wtN";
    if( (stchar.first == "wtF"  || stchar.first == "writerF" ) && (stchar.second == '[' || stchar.second == '(') ) mode = "wtF";

    if( (stchar.first == "Wt"   || stchar.first == "Writer"  ) && (stchar.second == '[' || stchar.second == '(') ) mode = "wt", extend = 0;
    if( (stchar.first == "WtSp" || stchar.first == "WriterSp") && (stchar.second == '[' || stchar.second == '(') ) mode = "wtSp", extend = 0;
    if( (stchar.first == "WtLn" || stchar.first == "WriterLn") && (stchar.second == '[' || stchar.second == '(') ) mode = "wtLn", extend = 0;
    if( (stchar.first == "WtN"  || stchar.first == "WriterN" ) && (stchar.second == '[' || stchar.second == '(') ) mode = "wtN", extend = 0;
    if( (stchar.first == "WtF"  || stchar.first == "WriterF" ) && (stchar.second == '[' || stchar.second == '(') ) mode = "wtF", extend = 0;

    if(mode=="") return tmp;

    string optarg = "", basestr = "";
    if(stchar.second == '['){
      int ss, tt;
      string stmp; vector<string> vtmp1, vtmp2;
      pair<string,char> ctmp;
      
      ss = strpos(tmp, "[");
      tt = pairBracket(tmp, ss);
      optarg = tmp.substr(ss, tt-ss+1);
      tmp = tmp.substr(tt+1);
      
      stmp = optarg.substr(1, optarg.size()-2);
      vtmp1 = split_p(stmp, ',');
      rep(k,vtmp1.size()){
        trim(vtmp1[k]);
        ctmp = nextToken(vtmp1[k]);
        if(ctmp.first == "B" && ctmp.second == '='){
          vtmp2 = split_p(vtmp1[k], '=');
          trim(vtmp2[1]);
          basestr = "," + vtmp2[1];
        }
      }
    }

    
    trim_until(tmp, '(', ')');
    tmp = tmp.substr(1,tmp.size()-2);


    if(mode == "wtF"){
      string wstr, estr, cstr;
      
      trim_until(tmp, '"', '"');
      tmp = tmp.substr(1, tmp.size() - 2);

      i = k = 0;
      rep(i,tmp.size()){
        if(tmp[i]=='{'){
          if(wstr.size()) cstr += (string)"wtN(\"" + wstr + "\");\n";
          wstr = "";
          k = 1;
          continue;
        }
        if(tmp[i]=='}'){
          if(estr.size()) cstr += (string)"wtN(" + estr + ");\n";
          estr = "";
          k = 0;
          continue;
        }
        if(tmp[i]=='\\' && (tmp[i+1]=='{' || tmp[i+1]=='}')) i++;
        if(k==0) wstr += tmp[i];
        else     estr += tmp[i];
      }
      if(wstr.size()) cstr += (string)"wtN(\"" + wstr + "\");\n";
      if(estr.size()) cstr += (string)"wtN(" + estr + ");\n";
      //fprintf(stderr, "[[%s]]\n", cstr.c_str());
      insert(cstr, str.size());
      
    } else {
      
      vtmp = split_p(tmp, ',');
      rep(k,vtmp.size()){
        trim(vtmp[k]);
        int fg = 0, char_fg = 0;
        string etype = getEquationType(vtmp[k], &fg, &char_fg);
        string stc;

        //fprintf(stderr, "wt [%s] [%s] [%d] [%d]\n", etype.c_str(), vtmp[k].c_str(), fg, char_fg);
        if(fg) ifun.doit.insert((string)"writer_all");
        if(char_fg) ifun.doit.insert((string)"writer_char_array");

        ifun.doit.insert((string)"writer_char");

        ifun.doit.insert((string)"writer_vector");
        ifun.doit.insert((string)"writer_set");
        ifun.doit.insert((string)"writer_multiset");
        ifun.doit.insert((string)"writer_pair");

        if(etype=="int" || etype=="signed" || etype=="signed int"){
          if(basestr=="") ifun.doit.insert((string)"writer_int");
          else            ifun.doit.insert((string)"writer_int_withBase");
        } else if(etype=="long long" || etype=="signed long long" || etype=="ll"){
          if(basestr=="") ifun.doit.insert((string)"writer_ll");
          else            ifun.doit.insert((string)"writer_ll_withBase");
        } else if(etype=="unsigned" || etype=="unsigned int" || etype=="uint"){
          ifun.doit.insert((string)"writer_unsigned");
        } else if(etype=="unsigned long long" || etype=="ull"){
          ifun.doit.insert((string)"writer_ull");
        } else if(etype=="int8_t" || etype=="int16_t" || etype=="int32_t" || etype=="int64_t" || etype=="__int8_t" || etype=="__int16_t" || etype=="__int32_t" || etype=="__int64_t" || etype=="__int128_t"){
          ifun.doit.insert((string)"writer_int");
          ifun.doit.insert((string)"writer_ll");
          ifun.doit.insert((string)"writer_int128");
        } else if(etype=="uint8_t" || etype=="uint16_t" || etype=="uint32_t" || etype=="uint64_t" || etype=="__uint8_t" || etype=="__uint16_t" || etype=="__uint32_t" || etype=="__uint64_t" || etype=="__uint128_t"){
          ifun.doit.insert((string)"writer_unsigned");
          ifun.doit.insert((string)"writer_ull");
          ifun.doit.insert((string)"writer_uint128");
        } else if(etype=="Modint"){
          ifun.doit.insert((string)"writer_Modint");
        } else if(etype=="Mint"){
          ifun.doit.insert((string)"writer_Mint");
        } else if(etype=="modint"){
          ifun.doit.insert((string)"writer_modint");
        } else if(etype=="mint"){
          ifun.doit.insert((string)"writer_mint");
        } else if(etype=="double"){
          ifun.doit.insert((string)"writer_double");
        } else if(etype=="char"){
          ifun.doit.insert((string)"writer_char_array");
        } else if(etype=="string"){
          ifun.doit.insert((string)"writer_string");
        } else {
          //fprintf(stderr, "unknown type [%s] for wt (writer) : %s\n", etype.c_str(), vtmp[k].c_str());
          //assert(0);
        }
        
        string strtmp;
        vector<string> vstrtmp;
        
        strtmp = vtmp[k];
        trim(strtmp);
        vstrtmp = rd_wt_array(strtmp);
        if(extend == 0){
          vstrtmp.clear();
          vstrtmp.push_back(strtmp);
        }
        if(vstrtmp.size() > 1){
          int m;
          vector<string> vn, looptm;
          string var;

          looptm = split_p(vstrtmp[2], ',');
          rep(m,looptm.size()) vn.push_back(getUnusedVarName());
          
          if(vstrtmp[1][0]=='('){
            vector<string> vars;
            vstrtmp[1] = vstrtmp[1].substr(1, vstrtmp[1].size()-2);
            vars = split_p(vstrtmp[1], ',');
            var = "";
            rep(j,vars.size()){
              if(j) var += ", ";
              trim(vars[j]);
              var += vstrtmp[0] + vars[j];
              rep(m,looptm.size()) var += (string)"[" + vn[m] + "]";
              var += vstrtmp[3];
            }
          } else {
            var = vstrtmp[0] + vstrtmp[1];
            rep(m,looptm.size()) var += (string)"[" + vn[m] + "]";
            var += vstrtmp[3];
          }

          stc = (string)"{";
          rep(m,looptm.size()) stc += (string)"int " + vn[m] + ";";
          rep(m,looptm.size()-1) stc += (string)" rep(" + vn[m] + ", " + looptm[m] + ")";
          if(looptm.size() >= 2) stc += "{";
          m = looptm.size() - 1;
          if(mode == "wt"){
            if(k==vtmp.size()-1){
              stc += (string)"if("+looptm[m]+"==0) wt_L('\\n');\n  else {\n    rep("+vn[m]+","+looptm[m]+"-1) wtSp" + optarg + "("+var+");\n    wt" + optarg + "("+var+");\n  }";
            } else {
              stc += (string)"rep(" + vn[m] + "," + looptm[m] + ") wtSp" + optarg + "(" + var + ");";
            }
          } else {
            stc += (string)"rep(" + vn[m] + "," + looptm[m] + ") " + mode + optarg + "(" + var + ");";
          }
          if(looptm.size() >= 2) stc += "}";
          stc += (string)"}";
          //fprintf(stderr, "----\n%s\n----\n",stc.c_str());
          insert(stc, str.size());
        } else {
          if(etype=="int" || etype=="long long") vstrtmp[0] += basestr;
          stc = (string)"wt_L(" + vstrtmp[0] + ");";
          
          str.push_back(stc);
          nxt.push_back(-1);
          strtype.push_back((string)"sentence");
          
          stc = "";
          if(mode == "wtLn" || (mode=="wt" && k==vtmp.size()-1)){
            stc = (string)"wt_L('\\n');";
          }
          if(mode == "wtSp" || (mode=="wt" && k!=vtmp.size()-1)){
            stc = (string)"wt_L(' ');";
          }
          if(stc != ""){
            str.push_back(stc);
            nxt.push_back(-1);
            strtype.push_back((string)"sentence");
          }
        }
      }
    }
    return "";
  }

  string sentence_dot_loop(string tmp){
    int k, pl, pb, pa, fg = 0;
    int k1, k2, k3, k4, k5;
    string dots = "..";
    string vn, stc, bg, ed, fbg, fed;

    if(strpos_ns(tmp, (string)"..") == -1) return tmp;
    
    while(strpos_ns(tmp, dots) >= 0) dots += ".";
    dots = dots.substr(1);
    
    vn = getUnusedVarName();
    stc = "";
    stc += (string)"{ int " + vn + ";";
    
    for(;;){
      pl = strpos_ns(tmp, dots);
      if(pl<0) break;
      
      k1 = k2 = k3 = k4 = k5 = 0;
      for(k=pl-1;k>=0;k--){
        if(k4==0 && k5==0){
          if(tmp[k] == '(') k1++;
          if(tmp[k] == ')') k1--;
          if(tmp[k] == '[') k2++;
          if(tmp[k] == ']') k2--;
          if(tmp[k] == '{') k3++;
          if(tmp[k] == '}') k3--;
        }
        if(k5==0 && tmp[k] == '\'') k4^=1;
        if(k4==0 && tmp[k] == '"') k5^=1;
        if(k1 > 0 || k2 > 0 || k3 > 0) break;
      }
      pb = k + 1;
      
      k1 = k2 = k3 = k4 = k5 = 0;
      REP(k,pl,tmp.size()){
        if(k4==0 && k5==0){
          if(tmp[k] == '(') k1++;
          if(tmp[k] == ')') k1--;
          if(tmp[k] == '[') k2++;
          if(tmp[k] == ']') k2--;
          if(tmp[k] == '{') k3++;
          if(tmp[k] == '}') k3--;
        }
        if(k5==0 && tmp[k] == '\'') k4^=1;
        if(k4==0 && tmp[k] == '"') k5^=1;
        if(k1 < 0 || k2 < 0 || k3 < 0) break;
      }
      pa = k;
      
      bg = tmp.substr(pb, pl-pb);
      ed = tmp.substr(pl+dots.size(), pa-pl-dots.size());
      
      if(fg==0){
        fg = 1;
        fbg = bg;
        fed = ed;
      }
      if(bg==fbg) tmp = tmp.substr(0, pb) + vn + tmp.substr(pa);
      else        tmp = tmp.substr(0, pb) + vn + " - (" + fbg + ") + (" + bg + ")" + tmp.substr(pa);
    }
    
    stc += "rep(" + vn + ", " + fbg + ", (" + fed + ") + 1) ";
    stc += tmp;
    
    stc += "}";
    insert(stc, str.size());
    return "";
  }



  string sentence_gcdlcm_sub(vector<string> vs, string op){
    while(vs.size() > 1){
      if(op=="+" || op=="*" || op=="^" || op=="|" || op=="&"){
        vs[0] = "(" + vs[0] + ")" + op + "(" + vs[1] +")";
      } else {
        vs[0] = op + "(" + vs[0] + ", " + vs[1] + ")";
      }
      vs.erase(vs.begin()+1);
    }
    return vs[0];
  }
  
  string sentence_gcdlcm(string tmp, int blocked = 1){
    int i, j, k;
    string stc, func, func_n1, func_n2;
    string vn, loop_vn, vtype;
    vector<string> vs, vsnx, arr, barr, varr, tmpvs, vsub;
    blocked = 0;

    for(;;){
      vs.clear();
      vs.push_back("");
      vs.push_back(tmp);
      vs.push_back("");
      func = func_n1 = func_n2 = "";
      for(;;){
        if(g_flags.count("no-gcd()")==0){
          vsnx = findFunction(vs[1],"gcd()");
          if(vsnx.size()==3){
            vs[0] = vs[0] + func_n1 + vsnx[0];
            vs[1] = vsnx[1];
            vs[2] = func_n2 + vsnx[2] + vs[2];
            func_n1 = "gcd(";
            func_n2 = ")";
            func = "GCD_L", ifun.doit.insert((string)"gcd");
            continue;
          }
        }
        if(g_flags.count("no-GCD()")==0){
          vsnx = findFunction(vs[1],"GCD()");
          if(vsnx.size()==3){
            vs[0] = vs[0] + func_n1 + vsnx[0];
            vs[1] = vsnx[1];
            vs[2] = func_n2 + vsnx[2] + vs[2];
            func_n1 = "GCD(";
            func_n2 = ")";
            func = "GCD_L", ifun.doit.insert((string)"gcd");
            continue;
          }
        }
        if(g_flags.count("no-lcm()")==0){
          vsnx = findFunction(vs[1],"lcm()");
          if(vsnx.size()==3){
            vs[0] = vs[0] + func_n1 + vsnx[0];
            vs[1] = vsnx[1];
            vs[2] = func_n2 + vsnx[2] + vs[2];
            func_n1 = "lcm(";
            func_n2 = ")";
            func = "LCM_L", ifun.doit.insert((string)"lcm");
            continue;
          }
        }
        if(g_flags.count("no-LCM()")==0){
          vsnx = findFunction(vs[1],"LCM()");
          if(vsnx.size()==3){
            vs[0] = vs[0] + func_n1 + vsnx[0];
            vs[1] = vsnx[1];
            vs[2] = func_n2 + vsnx[2] + vs[2];
            func_n1 = "LCM(";
            func_n2 = ")";
            func = "LCM_L", ifun.doit.insert((string)"lcm");
            continue;
          }
        }
        if(g_flags.count("no-min()")==0){
          vsnx = findFunction(vs[1],"min()");
          if(vsnx.size()==3 && vsnx[1].size()==0) vsnx.clear();
          if(vsnx.size()){
            string chk = vsnx[0];
            trim(chk);
            int len = chk.size();
            if(len >= 5 && chk.substr(len-5)=="std::") vsnx.clear();
          }
          if(vsnx.size()==3){
            vs[0] = vs[0] + func_n1 + vsnx[0];
            vs[1] = vsnx[1];
            vs[2] = func_n2 + vsnx[2] + vs[2];
            func_n1 = "min(";
            func_n2 = ")";
            func = "min_L", ifun.doit.insert((string)"min_L");
            continue;
          }
        }
        if(g_flags.count("no-MIN()")==0){
          vsnx = findFunction(vs[1],"MIN()");
          if(vsnx.size()==3){
            vs[0] = vs[0] + func_n1 + vsnx[0];
            vs[1] = vsnx[1];
            vs[2] = func_n2 + vsnx[2] + vs[2];
            func_n1 = "MIN(";
            func_n2 = ")";
            func = "min_L", ifun.doit.insert((string)"min_L");
            continue;
          }
        }
        if(g_flags.count("no-max()")==0){
          vsnx = findFunction(vs[1],"max()");
          if(vsnx.size()==3 && vsnx[1].size()==0) vsnx.clear();
          if(vsnx.size()){
            string chk = vsnx[0];
            trim(chk);
            int len = chk.size();
            if(len >= 5 && chk.substr(len-5)=="std::") vsnx.clear();
          }
          if(vsnx.size()==3){
            vs[0] = vs[0] + func_n1 + vsnx[0];
            vs[1] = vsnx[1];
            vs[2] = func_n2 + vsnx[2] + vs[2];
            func_n1 = "max(";
            func_n2 = ")";
            func = "max_L", ifun.doit.insert((string)"max_L");
            continue;
          }
        }
        if(g_flags.count("no-MAX()")==0){
          vsnx = findFunction(vs[1],"MAX()");
          if(vsnx.size()==3){
            vs[0] = vs[0] + func_n1 + vsnx[0];
            vs[1] = vsnx[1];
            vs[2] = func_n2 + vsnx[2] + vs[2];
            func_n1 = "MAX(";
            func_n2 = ")";
            func = "max_L", ifun.doit.insert((string)"max_L");
            continue;
          }
        }
        if(g_flags.count("no-sum()")==0){
          vsnx = findFunction(vs[1],"sum()");
          if(vsnx.size()==3){
            vs[0] = vs[0] + func_n1 + vsnx[0];
            vs[1] = vsnx[1];
            vs[2] = func_n2 + vsnx[2] + vs[2];
            func_n1 = "sum(";
            func_n2 = ")";
            func = "+";
            continue;
          }
        }
        if(g_flags.count("no-SUM()")==0){
          vsnx = findFunction(vs[1],"SUM()");
          if(vsnx.size()==3){
            vs[0] = vs[0] + func_n1 + vsnx[0];
            vs[1] = vsnx[1];
            vs[2] = func_n2 + vsnx[2] + vs[2];
            func_n1 = "SUM(";
            func_n2 = ")";
            func = "+";
            continue;
          }
        }
        if(g_flags.count("no-mul()")==0){
          vsnx = findFunction(vs[1],"mul()");
          if(vsnx.size()==3){
            vs[0] = vs[0] + func_n1 + vsnx[0];
            vs[1] = vsnx[1];
            vs[2] = func_n2 + vsnx[2] + vs[2];
            func_n1 = "mul(";
            func_n2 = ")";
            func = "*";
            continue;
          }
        }
        if(g_flags.count("no-MUL()")==0){
          vsnx = findFunction(vs[1],"MUL()");
          if(vsnx.size()==3){
            vs[0] = vs[0] + func_n1 + vsnx[0];
            vs[1] = vsnx[1];
            vs[2] = func_n2 + vsnx[2] + vs[2];
            func_n1 = "MUL(";
            func_n2 = ")";
            func = "*";
            continue;
          }
        }
        break;
      }
      if(func=="") break;


      stc = "";
      loop_vn = getUnusedVarName();
      arr = split_p(vs[1], ',');

      rep(k,arr.size()){
        trim(arr[k]);
        barr = rd_wt_array(arr[k]);
        if(barr.size()==4){
          vn = getUnusedVarName();
          vtype = getEquationType(barr[1]);

          stc += vtype + " " + vn + ";";

          if(barr[1][0]=='('){
            barr[1] = barr[1].substr(1, barr[1].size()-2);
            varr = split_p(barr[1], ',');
          } else {
            varr.clear();
            varr.push_back(barr[1]);
          }

          stc += "if(" + barr[2] + "==0){";
          if(func=="*" || func=="LCM_L"){
            stc += vn + " = 1;";
          } else if(vtype == "string"){
            stc += vn + ";";
          } else {
            stc += vn + " = 0;";
          }
          stc += "} else {";

          vsub.clear();
          rep(i,varr.size()) vsub.push_back(varr[i] + "[0]");
          stc += vn + " = " + sentence_gcdlcm_sub(vsub,func) + ";";

          stc += "rep(" + loop_vn + ", 1, " + barr[2] + "){";

          vsub.clear();
          rep(i,varr.size()) vsub.push_back(varr[i] + "[" + loop_vn + "]");
          if(func=="+" || func=="*"){
            stc += vn + " " + func + "= " + sentence_gcdlcm_sub(vsub,func) + ";";
          } else {
            stc += vn + " = " + func + "(" + vn + ", " + sentence_gcdlcm_sub(vsub,func) + ");";
          }
          stc += "}";
          stc += "}";

          arr[k] = vn;
        }
      }

      if(stc != ""){
        stc = "int " + loop_vn + ";" + stc;
        tmp = sentence_gcdlcm_sub(arr, func);
        if(blocked){
          stc = "{" + stc;
          stc += vs[0] + tmp + vs[2];
          stc += "}";
          insert(stc, str.size());
          return "";
        } else {
          insert(stc, str.size());
          tmp = vs[0] + tmp + vs[2];
        }
      } else {
        tmp = sentence_gcdlcm_sub(arr, func);
        tmp = vs[0] + tmp + vs[2];
      }
    }

    return tmp;
  }


  string sentence_gcdlcm_one(string tmp, int blocked = 1){
    int i, j, k;
    string stc, func, func_n1, func_n2, mode;
    string vn, loop_vn, vtype, eq;
    vector<string> vs, vsnx, arr, barr, varr, tmpvs, vsub;
    blocked = 0;

    for(;;){
      get_first_function(tmp, vs, func, mode);
      assert(mode == "gcdlcm");
      if(func == "gcd()" || func == "GCD()") func_n1 = "gcd(", func_n2 = ")", func = "GCD_L", ifun.doit.insert("gcd");
      if(func == "lcm()" || func == "LCM()") func_n1 = "lcm(", func_n2 = ")", func = "LCM_L", ifun.doit.insert("lcm");
      if(func == "min()" || func == "MIN()") func_n1 = "min(", func_n2 = ")", func = "min_L", ifun.doit.insert("min_L");
      if(func == "max()" || func == "MAX()") func_n1 = "max(", func_n2 = ")", func = "max_L", ifun.doit.insert("max_L");
      if(func == "sum()" || func == "SUM()") func_n1 = "sum(", func_n2 = ")", func = "+";
      if(func == "mul()" || func == "MUL()") func_n1 = "mul(", func_n2 = ")", func = "*";
      if(func == "XOR()") func_n1 = "XOR(", func_n2 = ")", func = "^";
      if(func == "OR()")  func_n1 = "OR(", func_n2 = ")", func = "|";
      if(func == "AND()") func_n1 = "AND(", func_n2 = ")", func = "&";
      if(func=="") break;

      stc = "";
      loop_vn = getUnusedVarName();
      arr = split_p(vs[1], ',');

      eq = "";
      rep(k,arr.size()){
        trim(arr[k]);
        if(k) eq += "+";
        barr = rd_wt_array(arr[k]);
        if(barr.size()==4) eq += get_type_with_decltype_basic_type(getEquationType(barr[1]));
        else               eq += "("+arr[k]+")";
      }
      vtype = "cLtraits_try_make_signed<remove_reference<decltype("+eq+")>::type>::type";

      rep(k,arr.size()){
        barr = rd_wt_array(arr[k]);
        if(barr.size()==4 && getEquationType(barr[1]) == "string") vtype = "string";
        if(barr.size()!=4 && getEquationType(arr[k]) == "string") vtype = "string";
      }
      if(vtype != "string") ifun.doit.insert("cLtraits_try_make_signed");

      rep(k,arr.size()){
        trim(arr[k]);
        barr = rd_wt_array(arr[k]);
        if(barr.size()==4){
          trim(barr[0]);
          trim(barr[3]);
          if(barr[0].size() || barr[3].size()) barr.clear();
        }
        /*if(barr.size()==4){
          vector<string> tmpvs;
          string chkstr, funstr, modestr;
          string s1 = barr[0], s2 = barr[3];
          int len;
          for(;;){
            chkstr = s1 + s2;
            get_first_function(chkstr, tmpvs, funstr, modestr);
            if(funstr == "not found") break;
            if(tmpvs[0].size() >= s1.size()) break;
            if(tmpvs[2].size() < s2.size()){ barr.clear(); break; }
            len = tmpvs[2].size() - s2.size();
            s1 = s1.substr(s1.size() - len);
          }
        }*/
        if(barr.size()==4){
          vn = getUnusedVarName();
          //vtype = getEquationType(barr[1]);

          stc += vtype + " " + vn + ";";

          if(barr[1][0]=='('){
            barr[1] = barr[1].substr(1, barr[1].size()-2);
            varr = split_p(barr[1], ',');
          } else {
            varr.clear();
            varr.push_back(barr[1]);
          }

          stc += "if(" + barr[2] + "==0){";
          if(func=="*" || func=="LCM_L"){
            stc += vn + " = 1;";
          } else if(vtype == "string"){
            stc += vn + ";";
          } else {
            stc += vn + " = 0;";
          }
          stc += "} else {";

          vsub.clear();
          rep(i,varr.size()) vsub.push_back(varr[i] + "[0]");
          stc += vn + " = " + sentence_gcdlcm_sub(vsub,func) + ";";

          stc += "rep(" + loop_vn + ", 1, " + barr[2] + "){";

          vsub.clear();
          rep(i,varr.size()) vsub.push_back(varr[i] + "[" + loop_vn + "]");
          if(func=="+" || func=="*" || func=="^" || func=="|" || func=="&"){
            stc += vn + " " + func + "= " + sentence_gcdlcm_sub(vsub,func) + ";";
          } else {
            stc += vn + " = " + func + "(" + vn + ", " + sentence_gcdlcm_sub(vsub,func) + ");";
          }
          stc += "}";
          stc += "}";

          arr[k] = vn;
        }
      }

      if(stc != ""){
        stc = "int " + loop_vn + ";" + stc;
        tmp = sentence_gcdlcm_sub(arr, func);
        if(blocked){
          stc = "{" + stc;
          stc += vs[0] + tmp + vs[2];
          stc += "}";
          insert(stc, str.size());
          return "";
        } else {
          insert(stc, str.size());
          tmp = vs[0] + tmp + vs[2];
        }
      } else {
        tmp = sentence_gcdlcm_sub(arr, func);
        tmp = vs[0] + tmp + vs[2];
      }
      break;
    }

    return tmp;
  }


  string sentence_twopointers(string tmp){
    string f_insert, f_erase, f_check, str_bef, str_aft;
    string varname_n, varname_i, varname_j, fname_insert, fname_erase, fname_check;
    string stc;
    vector<string> vs, varg;

    vs = findFunction(tmp, "TwoPointers()[][][]");
    if(vs.size() != 6) return tmp;
    str_bef = vs[0]; trim(str_bef);
    f_insert = vs[2]; trim(f_insert);
    f_erase = vs[3]; trim(f_erase);
    f_check = vs[4]; trim(f_check);
    str_aft = vs[5]; trim(str_aft);
    if(str_bef != "") fprintf(stderr, "twopointers: has some prefix\n");
    if(str_aft != ";") fprintf(stderr, "twopointers: has some surfix\n");

    varname_n = getUnusedVarName();
    varname_i = getUnusedVarName();
    varname_j = getUnusedVarName();
    fname_insert = getUnusedVarName();
    fname_erase = getUnusedVarName();
    fname_check = getUnusedVarName();

    varg = split_p2(vs[1], ','); // n, res, ind
    assert(varg.size() == 3 && "twopointers must has 3 args in ()");
    trim(varg[0]);
    trim(varg[1]);
    trim(varg[2]);

    stc += "auto " + fname_insert + " = [&](int " + varg[2] + "){\n";
    stc += f_insert;
    stc += "};\n";
    stc += "auto " + fname_erase + " = [&](int " + varg[2] + "){\n";
    stc += f_erase;
    stc += "};\n";
    stc += "auto " + fname_check + " = [&](void){\n";
    stc += f_check;
    stc += "};\n";

    stc += "if(true){\n";
    
    stc += "int " + varname_i + " = 0;\n";
    stc += "int " + varname_j + " = 0;\n";
    stc += "const int " + varname_n + " = " + varg[0] + ";\n";

    stc += "for(;;){\n";

    stc += "  if(" + fname_check + "()){\n";
    stc += "    " + varg[1] + "[" + varname_i + "] = " + varname_j + ";\n";
    stc += "    if(" + varname_j + " == " + varname_n + ") break;\n";
    stc += "    " + fname_insert + "(" + varname_j + "++);\n";
    stc += "  } else {\n";
    stc += "    " + fname_erase + "(" + varname_i + "++);\n";
    stc += "    " + varg[1] + "[" + varname_i + "] = " + varg[1] + "[" + varname_i + "-1];\n";
    stc += "  }\n";
    
    stc += "}\n"; // for(;;)
    stc += "while(" + varname_i + "+1 < " + varname_n + ") " + varg[1] + "[++" + varname_i + "] = " + varname_j + ";";
    stc += "}\n"; // if(true)

    insert(stc, str.size());
    return "";
  }


  string sentence_bsearch(string tmp){
    int i, j, k;
    int min_fg, max_fg, mode_int, setL, setU;
    string stc, varnameL, varnameU, varnameX, block, cond;
    vector<string> vs, var, vs_tmp;

    int AE_fg, RE_fg; string AE, RE;

    for(;;){
      vs.clear();
      min_fg = max_fg = 0;

      vs_tmp = findFunction(tmp, "bsearch_min[]()");
      if(vs_tmp.size()==4 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        min_fg=1;
        max_fg=0;
        vs.push_back("");
        vs[4] = vs[3];
        vs[3] = vs[2];
        vs[2] = "";
      }

      vs_tmp = findFunction(tmp, "bsearch_max[]()");
      if(vs_tmp.size()==4 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        max_fg=1;
        min_fg=0;
        vs.push_back("");
        vs[4] = vs[3];
        vs[3] = vs[2];
        vs[2] = "";
      }

      vs_tmp = findFunction(tmp, "bsearch_min[][]()");
      if(vs_tmp.size()==5 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        min_fg=1;
        max_fg=0;
      }

      vs_tmp = findFunction(tmp, "bsearch_max[][]()");
      if(vs_tmp.size()==5 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        max_fg=1;
        min_fg=0;
      }

      if(vs.size()==0) break;

      //fprintf(stderr, "bsearch\n");

      var = split_p2(vs[1], ',');
      rep(i,var.size()) trim(var[i]);

      AE_fg = RE_fg = 0;
      rep(i,var.size()){
        vector<string> arg = split_p2(var[i], '=');
        rep(j,arg.size()) trim(arg[j]);
        if(arg.size() == 2){
          if(arg[0] == "E"){
            AE_fg = RE_fg = 1;
            AE = RE = arg[1];
          } else if(arg[0] == "AE"){
            AE_fg = 1;
            AE = arg[1];
          } else if(arg[0] == "RE"){
            RE_fg = 1;
            RE = arg[1];
          } else {
            assert(0 && "invalid argment in bsearch_min|max");
          }
          var.erase( var.begin() + i );
          i--;
        }
      }
      
      varnameL = getUnusedVarName();
      varnameU = getUnusedVarName();
      varnameX = getUnusedVarName();

      setL = setU = 0;
      if(var.size() >= 3 && var[2] != "") setL = 1;
      if(var.size() >= 4 && var[3] != "") setU = 1;
      block = vs[2];
      cond = vs[3];
      //fprintf(stderr, "bsearch replace [%s] [%s] [%s]\n", block.c_str(), var[1].c_str(), varnameX.c_str());
      replaceAll_ns_t_with_leading_digit(block, var[1], varnameX);
      replaceAll_ns_t_with_leading_digit(cond, var[1], varnameX);

      stc = var[0] + " " + varnameL + ", " + varnameU + ", " + varnameX + ";";
      stc += varnameL + " = " + var[2] + ";";
      stc += varnameU + " = " + var[3] + ";";

      if(var[0]=="double" || var[0]=="float") mode_int = 0; else mode_int = 1;

      if(mode_int){
        stc += "while(" + varnameL + " < " + varnameU + "){\n";
        stc += "if(("+varnameL+" + "+varnameU+")%2==0){\n";
        stc += "  "+varnameX+" = ("+varnameL+" + "+varnameU+") / 2;\n";
        stc += "} else {\n";
        if(min_fg) stc += "  "+varnameX+" = ("+varnameL+" + "+varnameU+" - 1) / 2;\n";
        if(max_fg) stc += "  "+varnameX+" = ("+varnameL+" + "+varnameU+" + 1) / 2;\n";
        stc += "}\n";
        stc += block;
        if(min_fg) stc += "if("+cond+") "+varnameU+" = "+varnameX+"; else "+varnameL+" = "+varnameX+" + 1;\n";
        if(max_fg) stc += "if("+cond+") "+varnameL+" = "+varnameX+"; else "+varnameU+" = "+varnameX+" - 1;\n";
        stc += "}\n";
      } else {
        stc += "for(;;){\n";
        stc += varnameX+" = ("+varnameL+" + "+varnameU+") / 2;\n";
        if(AE_fg) stc += "if("+varnameU+" - "+varnameL+" < "+AE+") break;\n";
        if(RE_fg) stc += "if("+varnameL+" > 0 && "+varnameU+" - "+varnameL+" < "+varnameL+" * "+RE+") break;\n";
        if(RE_fg) stc += "if("+varnameU+" < 0 && "+varnameU+" - "+varnameL+" < (-"+varnameU+") * "+RE+") break;\n";
        stc += "if("+varnameX+" == "+varnameL+" || "+varnameX+" == "+varnameU+") break;\n";
        stc += block;
        if(min_fg) stc += "if("+cond+") "+varnameU+" = "+varnameX+"; else "+varnameL+" = "+varnameX+";\n";
        if(max_fg) stc += "if("+cond+") "+varnameL+" = "+varnameX+"; else "+varnameU+" = "+varnameX+";\n";
        stc += "}\n";
      }

      insert(stc, str.size());
      if(mode_int){
        tmp = stc + vs[0] + varnameU + vs[4];
      } else {
        tmp = stc + vs[0] + "(("+varnameL+" + "+varnameU+") / 2)" + vs[4];
      }
    }

    return tmp;
  }


  string sentence_bsearch_one(string tmp){
    int i, j, k;
    int min_fg, max_fg, mode_int, setL, setU;
    string stc, varnameL, varnameU, varnameX, block, cond;
    vector<string> vs, var, vs_tmp;

    int AE_fg, RE_fg; string AE, RE;

    for(;;){
      vs.clear();
      min_fg = max_fg = 0;

      vs_tmp = findFunction(tmp, "bsearch_min[]()");
      if(vs_tmp.size()==4 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        min_fg=1;
        max_fg=0;
        vs.push_back("");
        vs[4] = vs[3];
        vs[3] = vs[2];
        vs[2] = "";
      }

      vs_tmp = findFunction(tmp, "bsearch_max[]()");
      if(vs_tmp.size()==4 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        max_fg=1;
        min_fg=0;
        vs.push_back("");
        vs[4] = vs[3];
        vs[3] = vs[2];
        vs[2] = "";
      }

      vs_tmp = findFunction(tmp, "bsearch_min[][]()");
      if(vs_tmp.size()==5 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        min_fg=1;
        max_fg=0;
      }

      vs_tmp = findFunction(tmp, "bsearch_max[][]()");
      if(vs_tmp.size()==5 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        max_fg=1;
        min_fg=0;
      }

      if(vs.size()==0) break;

      //fprintf(stderr, "bsearch\n");

      var = split_p2(vs[1], ',');
      rep(i,var.size()) trim(var[i]);

      AE_fg = RE_fg = 0;
      rep(i,var.size()){
        vector<string> arg = split_p2(var[i], '=');
        rep(j,arg.size()) trim(arg[j]);
        if(arg.size() == 2){
          if(arg[0] == "E"){
            AE_fg = RE_fg = 1;
            AE = RE = arg[1];
          } else if(arg[0] == "AE"){
            AE_fg = 1;
            AE = arg[1];
          } else if(arg[0] == "RE"){
            RE_fg = 1;
            RE = arg[1];
          } else {
            assert(0 && "invalid argment in bsearch_min|max");
          }
          var.erase( var.begin() + i );
          i--;
        }
      }
      
      varnameL = getUnusedVarName();
      varnameU = getUnusedVarName();
      varnameX = getUnusedVarName();

      setL = setU = 0;
      if(var.size() >= 3 && var[2] != "") setL = 1;
      if(var.size() >= 4 && var[3] != "") setU = 1;
      block = vs[2];
      cond = vs[3];
      //fprintf(stderr, "bsearch replace [%s] [%s] [%s]\n", block.c_str(), var[1].c_str(), varnameX.c_str());
      replaceAll_ns_t_with_leading_digit(block, var[1], varnameX);
      replaceAll_ns_t_with_leading_digit(cond, var[1], varnameX);

      stc = var[0] + " " + varnameL + ", " + varnameU + ", " + varnameX + ";";
      stc += varnameL + " = " + var[2] + ";";
      stc += varnameU + " = " + var[3] + ";";

      if(var[0]=="double" || var[0]=="float") mode_int = 0; else mode_int = 1;

      if(mode_int){
        stc += "while(" + varnameL + " < " + varnameU + "){\n";
        stc += "if(("+varnameL+" + "+varnameU+")%2==0){\n";
        stc += "  "+varnameX+" = ("+varnameL+" + "+varnameU+") / 2;\n";
        stc += "} else {\n";
        if(min_fg) stc += "  "+varnameX+" = ("+varnameL+" + "+varnameU+" - 1) / 2;\n";
        if(max_fg) stc += "  "+varnameX+" = ("+varnameL+" + "+varnameU+" + 1) / 2;\n";
        stc += "}\n";
        stc += block;
        if(min_fg) stc += "if("+cond+") "+varnameU+" = "+varnameX+"; else "+varnameL+" = "+varnameX+" + 1;\n";
        if(max_fg) stc += "if("+cond+") "+varnameL+" = "+varnameX+"; else "+varnameU+" = "+varnameX+" - 1;\n";
        stc += "}\n";
      } else {
        stc += "for(;;){\n";
        stc += varnameX+" = ("+varnameL+" + "+varnameU+") / 2;\n";
        if(AE_fg) stc += "if("+varnameU+" - "+varnameL+" < "+AE+") break;\n";
        if(RE_fg) stc += "if("+varnameL+" > 0 && "+varnameU+" - "+varnameL+" < "+varnameL+" * "+RE+") break;\n";
        if(RE_fg) stc += "if("+varnameU+" < 0 && "+varnameU+" - "+varnameL+" < (-"+varnameU+") * "+RE+") break;\n";
        stc += "if("+varnameX+" == "+varnameL+" || "+varnameX+" == "+varnameU+") break;\n";
        stc += block;
        if(min_fg) stc += "if("+cond+") "+varnameU+" = "+varnameX+"; else "+varnameL+" = "+varnameX+";\n";
        if(max_fg) stc += "if("+cond+") "+varnameL+" = "+varnameX+"; else "+varnameU+" = "+varnameX+";\n";
        stc += "}\n";
      }

      insert(stc, str.size());
      if(mode_int){
        tmp = stc + vs[0] + varnameU + vs[4];
      } else {
        tmp = stc + vs[0] + "(("+varnameL+" + "+varnameU+") / 2)" + vs[4];
      }
      break;
    }

    return tmp;
  }


  string sentence_tsearch(string tmp){
    int i, j, k;
    int min_fg, max_fg, arg_fg, las_fg, mode_int, setL, setU;
    string stc, varnameL, varnameU, varnameM, varnameX, varnameLV, varnameUV, varnameMV, varnameXV, block, cond;
    string modestr;
    vector<string> vs, var, vs_tmp;

    int AE_fg, RE_fg; string AE, RE, TYPE;

    for(;;){
      vs.clear();
      min_fg = max_fg = 0;

      vs_tmp = findFunction(tmp, "tsearch_min[]()");
      if(vs_tmp.size()==4 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        min_fg = 1;
        max_fg = 0;
        arg_fg = 0;
        las_fg = 0;
        modestr = "<=";
        vs.push_back("");
        vs[4] = vs[3];
        vs[3] = vs[2];
        vs[2] = "";
      }

      vs_tmp = findFunction(tmp, "tsearch_argmin[]()");
      if(vs_tmp.size()==4 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        min_fg = 1;
        max_fg = 0;
        arg_fg = 1;
        las_fg = 0;
        modestr = "<=";
        vs.push_back("");
        vs[4] = vs[3];
        vs[3] = vs[2];
        vs[2] = "";
      }

      vs_tmp = findFunction(tmp, "tsearch_argminL[]()");
      if(vs_tmp.size()==4 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        min_fg = 1;
        max_fg = 0;
        arg_fg = 1;
        las_fg = 1;
        modestr = "<";
        vs.push_back("");
        vs[4] = vs[3];
        vs[3] = vs[2];
        vs[2] = "";
      }

      vs_tmp = findFunction(tmp, "tsearch_max[]()");
      if(vs_tmp.size()==4 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        max_fg = 1;
        min_fg = 0;
        arg_fg = 0;
        las_fg = 0;
        modestr = ">=";
        vs.push_back("");
        vs[4] = vs[3];
        vs[3] = vs[2];
        vs[2] = "";
      }

      vs_tmp = findFunction(tmp, "tsearch_argmax[]()");
      if(vs_tmp.size()==4 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        max_fg = 1;
        min_fg = 0;
        arg_fg = 1;
        las_fg = 0;
        modestr = ">=";
        vs.push_back("");
        vs[4] = vs[3];
        vs[3] = vs[2];
        vs[2] = "";
      }

      vs_tmp = findFunction(tmp, "tsearch_argmaxL[]()");
      if(vs_tmp.size()==4 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        max_fg = 1;
        min_fg = 0;
        arg_fg = 1;
        las_fg = 1;
        modestr = ">";
        vs.push_back("");
        vs[4] = vs[3];
        vs[3] = vs[2];
        vs[2] = "";
      }

      vs_tmp = findFunction(tmp, "tsearch_min[][]()");
      if(vs_tmp.size()==5 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        min_fg = 1;
        max_fg = 0;
        arg_fg = 0;
        las_fg = 0;
        modestr = "<=";
      }

      vs_tmp = findFunction(tmp, "tsearch_argmin[][]()");
      if(vs_tmp.size()==5 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        min_fg = 1;
        max_fg = 0;
        arg_fg = 1;
        las_fg = 0;
        modestr = "<=";
      }

      vs_tmp = findFunction(tmp, "tsearch_argminL[][]()");
      if(vs_tmp.size()==5 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        min_fg = 1;
        max_fg = 0;
        arg_fg = 1;
        las_fg = 1;
        modestr = "<";
      }

      vs_tmp = findFunction(tmp, "tsearch_max[][]()");
      if(vs_tmp.size()==5 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        max_fg = 1;
        min_fg = 0;
        arg_fg = 0;
        las_fg = 0;
        modestr = ">=";
      }
      
      vs_tmp = findFunction(tmp, "tsearch_argmax[][]()");
      if(vs_tmp.size()==5 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        max_fg = 1;
        min_fg = 0;
        arg_fg = 1;
        las_fg = 0;
        modestr = ">=";
      }
      
      vs_tmp = findFunction(tmp, "tsearch_argmaxL[][]()");
      if(vs_tmp.size()==5 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        max_fg = 1;
        min_fg = 0;
        arg_fg = 1;
        las_fg = 1;
        modestr = ">";
      }

      if(vs.size()==0) break;

      //fprintf(stderr, "tsearch\n");
      // vs[0] : 
      // vs[1] : 設定 [int,x,0,100,hoge]
      // vs[2] : ブロック
      // vs[3] : 式
      // vs[4] : 

      var = split_p2(vs[1], ',');
      rep(i,var.size()) trim(var[i]);

      AE_fg = RE_fg = 0;
      TYPE = "";
      rep(i,var.size()){
        vector<string> arg = split_p2(var[i], '=');
        rep(j,arg.size()) trim(arg[j]);
        if(arg.size() == 2){
          if(arg[0] == "E"){
            AE_fg = RE_fg = 1;
            AE = RE = arg[1];
          } else if(arg[0] == "AE"){
            AE_fg = 1;
            AE = arg[1];
          } else if(arg[0] == "RE"){
            RE_fg = 1;
            RE = arg[1];
          } else if(arg[0] == "T"){
            TYPE = arg[1];
          } else {
            assert(0 && "invalid argment in tsearch_[arg]min|max");
          }
          var.erase( var.begin() + i );
          i--;
        }
      }

      varnameL = getUnusedVarName();
      varnameU = getUnusedVarName();
      varnameM = getUnusedVarName();
      varnameX = getUnusedVarName();
      varnameLV = getUnusedVarName();
      varnameUV = getUnusedVarName();
      varnameMV = getUnusedVarName();
      varnameXV = getUnusedVarName();

      setL = setU = 0;
      if(var.size() >= 3 && var[2] != "") setL = 1;
      if(var.size() >= 4 && var[3] != "") setU = 1;
      block = vs[2];
      cond = vs[3];
      //fprintf(stderr, "bsearch replace [%s] [%s] [%s]\n", block.c_str(), var[1].c_str(), varnameX.c_str());
      replaceAll_ns_t_with_leading_digit(block, var[1], varnameX);
      replaceAll_ns_t_with_leading_digit(cond, var[1], varnameX);

      if(TYPE=="") TYPE = "remove_reference<decltype(" + cond + ")>::type";
      stc = "";

      stc += var[0] + " " + varnameL + ", " + varnameU + ", " + varnameM + ", " + varnameX + ";";
      stc += varnameL + " = " + var[2] + ";";
      stc += varnameU + " = " + var[3] + ";";
      stc += TYPE + " " + varnameLV + ";";
      stc += TYPE + " " + varnameUV + ";";
      stc += TYPE + " " + varnameMV + ";";
      stc += TYPE + " " + varnameXV + ";";

      if(var[0]=="double" || var[0]=="float") mode_int = 0; else mode_int = 1;

      stc += varnameX + " = " + varnameL + ";";
      stc += "{";
      stc += block;
      stc += varnameLV + " = " + cond + ";";
      stc += "}";
      stc += varnameX + " = " + varnameU + ";";
      stc += "{";
      stc += block;
      stc += varnameUV + " = " + cond + ";";
      stc += "}";
      stc += varnameX + " = " + varnameM + " = (" + varnameL + " + " + varnameU + ") / 2;";
      stc += "{";
      stc += block;
      stc += varnameMV + " = " + cond + ";";
      stc += "}";
      stc += "for(;;){";
      stc += "  if(" + varnameU + " - " + varnameM + " > " + varnameM + " - " + varnameL + "){";
      stc += varnameX + " = (" + varnameM + " + " + varnameU + ") / 2;";
      stc += "    if(" + varnameX + " == " + varnameM + " || " + varnameX + " == " + varnameU + ") break;";
      stc += "  } else {";
      stc += varnameX + " = (" + varnameL + " + " + varnameM + ") / 2;";
      stc += "    if(" + varnameX + " == " + varnameL + " || " + varnameX + " == " + varnameM + ") break;";
      stc += "  }";
      stc += "  {";
      stc += block;
      stc += varnameXV + " = " + cond + ";";
      stc += "  }";
      stc += "  if(" + varnameU + " - " + varnameM + " > " + varnameM + " - " + varnameL + "){";
      stc += "    if(" + varnameMV + " " +  modestr + " " + varnameXV + "){";
      stc += varnameU + " = " + varnameX + ";";
      stc += varnameUV + " = " + varnameXV + ";";
      stc += "    } else {";
      stc += varnameL + " = " + varnameM + ";";
      stc += varnameLV + " = " + varnameMV + ";";
      stc += varnameM + " = " + varnameX + ";";
      stc += varnameMV + " = " + varnameXV + ";";
      stc += "    }";
      stc += "  } else {";
      stc += "    if(" + varnameXV + " " + modestr + " " + varnameMV + "){";
      stc += varnameU + " = " + varnameM + ";";
      stc += varnameUV + " = " + varnameMV + ";";
      stc += varnameM + " = " + varnameX + ";";
      stc += varnameMV + " = " + varnameXV + ";";
      stc += "    } else {";
      stc += varnameL + " = " + varnameX + ";";
      stc += varnameLV + " = " + varnameXV + ";";
      stc += "    }";
      stc += "  }";
      stc += "}";

      stc += varnameX + " = " + varnameU + ";";
      stc += varnameXV + " = " + varnameUV + ";";
      stc += "if(" + varnameMV + " " + modestr + " " + varnameXV + "){";
      stc += varnameX + " = " + varnameM + ";";
      stc += varnameXV + " = " + varnameMV + ";";
      stc += "}";
      stc += "if(" + varnameLV + " " + modestr + " " + varnameXV + "){";
      stc += varnameX + " = " + varnameL + ";";
      stc += varnameXV + " = " + varnameLV + ";";
      stc += "}";

      insert(stc, str.size());
      if(arg_fg){
        tmp = stc + vs[0] + varnameX + vs[4];
      } else {
        tmp = stc + vs[0] + varnameXV + vs[4];
      }

      stc += "}";
    }

    return tmp;
  }


  string sentence_tsearch_one(string tmp){
    int i, j, k;
    int min_fg, max_fg, arg_fg, las_fg, mode_int, setL, setU;
    string stc, varnameL, varnameU, varnameM, varnameX, varnameLV, varnameUV, varnameMV, varnameXV, block, cond;
    string modestr;
    vector<string> vs, var, vs_tmp;

    int AE_fg, RE_fg; string AE, RE, TYPE;

    for(;;){
      vs.clear();
      min_fg = max_fg = 0;

      vs_tmp = findFunction(tmp, "tsearch_min[]()");
      if(vs_tmp.size()==4 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        min_fg = 1;
        max_fg = 0;
        arg_fg = 0;
        las_fg = 0;
        modestr = "<=";
        vs.push_back("");
        vs[4] = vs[3];
        vs[3] = vs[2];
        vs[2] = "";
      }

      vs_tmp = findFunction(tmp, "tsearch_argmin[]()");
      if(vs_tmp.size()==4 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        min_fg = 1;
        max_fg = 0;
        arg_fg = 1;
        las_fg = 0;
        modestr = "<=";
        vs.push_back("");
        vs[4] = vs[3];
        vs[3] = vs[2];
        vs[2] = "";
      }

      vs_tmp = findFunction(tmp, "tsearch_argminL[]()");
      if(vs_tmp.size()==4 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        min_fg = 1;
        max_fg = 0;
        arg_fg = 1;
        las_fg = 1;
        modestr = "<";
        vs.push_back("");
        vs[4] = vs[3];
        vs[3] = vs[2];
        vs[2] = "";
      }

      vs_tmp = findFunction(tmp, "tsearch_max[]()");
      if(vs_tmp.size()==4 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        max_fg = 1;
        min_fg = 0;
        arg_fg = 0;
        las_fg = 0;
        modestr = ">=";
        vs.push_back("");
        vs[4] = vs[3];
        vs[3] = vs[2];
        vs[2] = "";
      }

      vs_tmp = findFunction(tmp, "tsearch_argmax[]()");
      if(vs_tmp.size()==4 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        max_fg = 1;
        min_fg = 0;
        arg_fg = 1;
        las_fg = 0;
        modestr = ">=";
        vs.push_back("");
        vs[4] = vs[3];
        vs[3] = vs[2];
        vs[2] = "";
      }

      vs_tmp = findFunction(tmp, "tsearch_argmaxL[]()");
      if(vs_tmp.size()==4 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        max_fg = 1;
        min_fg = 0;
        arg_fg = 1;
        las_fg = 1;
        modestr = ">";
        vs.push_back("");
        vs[4] = vs[3];
        vs[3] = vs[2];
        vs[2] = "";
      }

      vs_tmp = findFunction(tmp, "tsearch_min[][]()");
      if(vs_tmp.size()==5 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        min_fg = 1;
        max_fg = 0;
        arg_fg = 0;
        las_fg = 0;
        modestr = "<=";
      }

      vs_tmp = findFunction(tmp, "tsearch_argmin[][]()");
      if(vs_tmp.size()==5 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        min_fg = 1;
        max_fg = 0;
        arg_fg = 1;
        las_fg = 0;
        modestr = "<=";
      }

      vs_tmp = findFunction(tmp, "tsearch_argminL[][]()");
      if(vs_tmp.size()==5 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        min_fg = 1;
        max_fg = 0;
        arg_fg = 1;
        las_fg = 1;
        modestr = "<";
      }

      vs_tmp = findFunction(tmp, "tsearch_max[][]()");
      if(vs_tmp.size()==5 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        max_fg = 1;
        min_fg = 0;
        arg_fg = 0;
        las_fg = 0;
        modestr = ">=";
      }
      
      vs_tmp = findFunction(tmp, "tsearch_argmax[][]()");
      if(vs_tmp.size()==5 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        max_fg = 1;
        min_fg = 0;
        arg_fg = 1;
        las_fg = 0;
        modestr = ">=";
      }
      
      vs_tmp = findFunction(tmp, "tsearch_argmaxL[][]()");
      if(vs_tmp.size()==5 && (vs.size()==0 || vs[0].size() > vs_tmp[0].size())){
        vs = vs_tmp;
        max_fg = 1;
        min_fg = 0;
        arg_fg = 1;
        las_fg = 1;
        modestr = ">";
      }

      if(vs.size()==0) break;

      //fprintf(stderr, "tsearch\n");
      // vs[0] : 
      // vs[1] : 設定 [int,x,0,100,hoge]
      // vs[2] : ブロック
      // vs[3] : 式
      // vs[4] : 

      var = split_p2(vs[1], ',');
      rep(i,var.size()) trim(var[i]);

      AE_fg = RE_fg = 0;
      TYPE = "";
      rep(i,var.size()){
        vector<string> arg = split_p2(var[i], '=');
        rep(j,arg.size()) trim(arg[j]);
        if(arg.size() == 2){
          if(arg[0] == "E"){
            AE_fg = RE_fg = 1;
            AE = RE = arg[1];
          } else if(arg[0] == "AE"){
            AE_fg = 1;
            AE = arg[1];
          } else if(arg[0] == "RE"){
            RE_fg = 1;
            RE = arg[1];
          } else if(arg[0] == "T"){
            TYPE = arg[1];
          } else {
            assert(0 && "invalid argment in tsearch_[arg]min|max");
          }
          var.erase( var.begin() + i );
          i--;
        }
      }

      varnameL = getUnusedVarName();
      varnameU = getUnusedVarName();
      varnameM = getUnusedVarName();
      varnameX = getUnusedVarName();
      varnameLV = getUnusedVarName();
      varnameUV = getUnusedVarName();
      varnameMV = getUnusedVarName();
      varnameXV = getUnusedVarName();

      setL = setU = 0;
      if(var.size() >= 3 && var[2] != "") setL = 1;
      if(var.size() >= 4 && var[3] != "") setU = 1;
      block = vs[2];
      cond = vs[3];
      //fprintf(stderr, "bsearch replace [%s] [%s] [%s]\n", block.c_str(), var[1].c_str(), varnameX.c_str());
      replaceAll_ns_t_with_leading_digit(block, var[1], varnameX);
      replaceAll_ns_t_with_leading_digit(cond, var[1], varnameX);

      if(TYPE=="") TYPE = "remove_reference<decltype(" + cond + ")>::type";
      stc = "";

      stc += var[0] + " " + varnameL + ", " + varnameU + ", " + varnameM + ", " + varnameX + ";";
      stc += varnameL + " = " + var[2] + ";";
      stc += varnameU + " = " + var[3] + ";";
      stc += TYPE + " " + varnameLV + ";";
      stc += TYPE + " " + varnameUV + ";";
      stc += TYPE + " " + varnameMV + ";";
      stc += TYPE + " " + varnameXV + ";";

      if(var[0]=="double" || var[0]=="float") mode_int = 0; else mode_int = 1;

      stc += varnameX + " = " + varnameL + ";";
      stc += "{";
      stc += block;
      stc += varnameLV + " = " + cond + ";";
      stc += "}";
      stc += varnameX + " = " + varnameU + ";";
      stc += "{";
      stc += block;
      stc += varnameUV + " = " + cond + ";";
      stc += "}";
      stc += varnameX + " = " + varnameM + " = (" + varnameL + " + " + varnameU + ") / 2;";
      stc += "{";
      stc += block;
      stc += varnameMV + " = " + cond + ";";
      stc += "}";
      stc += "for(;;){";
      stc += "  if(" + varnameU + " - " + varnameM + " > " + varnameM + " - " + varnameL + "){";
      stc += varnameX + " = (" + varnameM + " + " + varnameU + ") / 2;";
      stc += "    if(" + varnameX + " == " + varnameM + " || " + varnameX + " == " + varnameU + ") break;";
      stc += "  } else {";
      stc += varnameX + " = (" + varnameL + " + " + varnameM + ") / 2;";
      stc += "    if(" + varnameX + " == " + varnameL + " || " + varnameX + " == " + varnameM + ") break;";
      stc += "  }";
      stc += "  {";
      stc += block;
      stc += varnameXV + " = " + cond + ";";
      stc += "  }";
      stc += "  if(" + varnameU + " - " + varnameM + " > " + varnameM + " - " + varnameL + "){";
      stc += "    if(" + varnameMV + " " +  modestr + " " + varnameXV + "){";
      stc += varnameU + " = " + varnameX + ";";
      stc += varnameUV + " = " + varnameXV + ";";
      stc += "    } else {";
      stc += varnameL + " = " + varnameM + ";";
      stc += varnameLV + " = " + varnameMV + ";";
      stc += varnameM + " = " + varnameX + ";";
      stc += varnameMV + " = " + varnameXV + ";";
      stc += "    }";
      stc += "  } else {";
      stc += "    if(" + varnameXV + " " + modestr + " " + varnameMV + "){";
      stc += varnameU + " = " + varnameM + ";";
      stc += varnameUV + " = " + varnameMV + ";";
      stc += varnameM + " = " + varnameX + ";";
      stc += varnameMV + " = " + varnameXV + ";";
      stc += "    } else {";
      stc += varnameL + " = " + varnameX + ";";
      stc += varnameLV + " = " + varnameXV + ";";
      stc += "    }";
      stc += "  }";
      stc += "}";

      stc += varnameX + " = " + varnameU + ";";
      stc += varnameXV + " = " + varnameUV + ";";
      stc += "if(" + varnameMV + " " + modestr + " " + varnameXV + "){";
      stc += varnameX + " = " + varnameM + ";";
      stc += varnameXV + " = " + varnameMV + ";";
      stc += "}";
      stc += "if(" + varnameLV + " " + modestr + " " + varnameXV + "){";
      stc += varnameX + " = " + varnameL + ";";
      stc += varnameXV + " = " + varnameLV + ";";
      stc += "}";

      insert(stc, str.size());
      if(arg_fg){
        tmp = stc + vs[0] + varnameX + vs[4];
      } else {
        tmp = stc + vs[0] + varnameXV + vs[4];
      }

      stc += "}";
      break;
    }

    return tmp;
  }


  string sentence_if(string tmp){
    int i;
    string stc;
    vector<string> vs, var;

    vs = findFunction(tmp, "if[]");
    if(vs.size()==0) return tmp;

    var = split_p(vs[1], ',');
    if(var.size()%2==0) var.push_back((string)"");
    rep(i,var.size()) trim(var[i]);
    // rep(i,var.size())fprintf(stderr, "--- %d %s\n", i, var[i].c_str());
    
    stc = "";
    for(i=0;i+1<var.size();i+=2){
      if(i) stc += "else ";
      stc += "if(" + var[i] + "){";
      stc += vs[0] + var[i+1] + vs[2];
      stc += "}";
    }
    stc += "else{";
    stc += vs[0] + var[i] + vs[2];
    stc += "}";

    insert(stc, str.size());
    return "";
  }



  string sentence_arraylike_operations(string tmp){
    int i, j, k, m, n, x, y;
    string stc, op, lef, rigstr, op1, op2, tmpstr;
    char lc;
    vector<string> vl, vr, vname;
    vector<vector<string>> rig;
    if(g_flags.count((string)"no-arraylike-operations")) return tmp;

    // ++(a,b,c)--++;
    rep(i,tmp.size()) if(tmp[i] == '(') break;
    if(i == tmp.size()) return tmp;
    j = pairBracket(tmp, i);
    x = y = 0;
    rep(k,tmp.size()) if(k < i || j < k) if(!isspace(tmp[k]) && tmp[k] != ';'){
      if(tmp[k] == '+' || tmp[k] == '-') x++;
      else                               y++;
    }
    if(x > 0 && y == 0){
      lef = tmp.substr(i+1,j-i-1);
      vl = split_p(lef, ',');
      if(vl.size() <= 1) return tmp;
      op1 = tmp.substr(0, i);
      op2 = tmp.substr(j+1);
      rep(i,vl.size()) stc += op1 + vl[i] + op2;
      insert(stc, str.size());
      return "";
    }

    // (a,b) += (c,d) + e * (f,g)
    if(tmp[0]!='(') return tmp;
    i = 0;
    j = pairBracket(tmp, i);

    lef = tmp.substr(i+1, j-i-1);
    vl = split_p(lef, ',');
    if(vl.size() <= 1) return tmp;

    k = j+1;
    while(k < tmp.size() && isspace(tmp[k])) k++;
    if(k==tmp.size()) return tmp;

    m = k;
    while(m < tmp.size() && (tmp[m]=='=' || tmp[m]=='+' || tmp[m]=='-' || tmp[m]=='*' || tmp[m]=='/' || tmp[m]=='>' || tmp[m]=='<' || tmp[m]=='|' || tmp[m]=='&' || tmp[m]=='^' || tmp[m]=='%' || tmp[m]=='?')) m++;
    if(m==tmp.size()) return tmp;

    op = tmp.substr(k,m-k);
    rep(x,op.size()) if(op[x]=='=') break;
    if(x==op.size()) return tmp;

    rigstr = tmp.substr(m, tmp.size()-m-1);
    for(;;){
      trim(rigstr);
      if(rigstr.size()==0) break;

      if(rigstr[0] == '('){
        i = 0;
        j = pairBracket(rigstr, i);
        tmpstr = rigstr.substr(i+1, j-i-1);
        vr = split_p(tmpstr, ',');
        if(vr.size() > 1){
          rep(k,vr.size()) vr[k] = "(" + vr[k] + ")";
          rig.push_back(vr);
          rigstr = rigstr.substr(j+1);
          continue;
        }
      }

      lc = ' ';
      rep(i,rigstr.size()){
        if(rigstr[i] == '(' && (lc == '+' || lc == '-' || lc == '*' || lc == '/' || lc == '%' || lc == '&' || lc == '|' || lc == '^' || lc == '<' || lc == '>' || lc == '=' || lc == '(')){
          vr.clear();
          vr.push_back(rigstr.substr(0, i));
          rig.push_back(vr);
          rigstr = rigstr.substr(i);
          i = -1;
          break;
        }
        if(!isspace(rigstr[i])) lc = rigstr[i];
      }
      if(i == -1) continue;

      vr.clear();
      vr.push_back(rigstr);
      rig.push_back(vr);
      rigstr = "";
    }

    x = 1;
    rep(i,rig.size()) x = max(x, (int)rig[i].size());
    rep(i,x) vname.push_back(getUnusedVarName());

    //rep(i,rig.size()) if(rig[i].size() < x) rep(j,rig[i].size()){
    //  string nm = getUnusedVarName();
    //  stc += "auto " + nm + " = (" + rig[i][j] + ");";
    //  rig[i][j] = nm;
    //}

    rep(i,x){
      stc += "auto " + vname[i] + " = (";
      rep(j,rig.size()) stc += rig[j][i%rig[j].size()];
      stc += ");";
    }
    rep(i,vl.size()) stc += vl[i] + op + vname[i%x] + ";";
    insert(stc, str.size());

    return "";
  }

  string sentence_arraylike_sentence(string tmp){
    int i, k, mx;
    string stc, sep, logstr = tmp;
    vector<string> gt, vs;
    vector<vector<string>> varr;
    if(g_flags.count((string)"no-arraylike-sentence")) return tmp;

    REP(i,-1,10){
      sep = "@";
      if(i >= 0) sep += to_string(i);
      gt = findFunction(tmp, sep+"[]");
      if(gt.size() == 0) continue;

      vs.push_back(gt[0]);
      varr.push_back(split_p(gt[1],','));
      tmp = gt[2];
      for(;;){
        gt = findFunction(tmp, sep+"[]");
        if(gt.size()==0) break;
        vs.push_back(gt[0]);
        varr.push_back(split_p(gt[1],','));
        tmp = gt[2];
      }
      vs.push_back(tmp);

      mx = 0;
      rep(i,varr.size()) if(mx < varr[i].size()) mx = varr[i].size();

      rep(k,mx){
        rep(i,varr.size()){
          stc += vs[i];
          stc += varr[i][k % varr[i].size()];
        }
        stc += vs[i];
      }
      
      insert(stc, str.size());
      return "";
    }
    return tmp;
  }

  string sentence_otherfunctions(string tmp){
    vector<string> vs, vtmp;
    static int readerfile = 0;
    static int writerfile = 0;

    for(;;){
      vs = findFunction(tmp, "readerFile()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"readerFile");
        ifun.already.insert((string)"my_getchar_unlocked");
        if(!readerfile){
          int i;
          readerfile = 1;
          rep(i, ifun.name.size()) if(ifun.name[i].substr(0,7)=="reader_" || ifun.name[i]=="rdLine"){
            string str = ifun.func[ifun.name[i]];
            for(;;){
              vtmp = findFunction(str, "my_getchar_unlocked()");
              if(vtmp.size()!=3) break;
              str = vtmp[0] + "getc(readerfp)" + vtmp[2];
            }
            for(;;){
              vtmp = findFunction(str, "getchar_unlocked()");
              if(vtmp.size()!=3) break;
              str = vtmp[0] + "getc(readerfp)" + vtmp[2];
            }
            for(;;){
              vtmp = findFunction(str, "my_getchar()");
              if(vtmp.size()!=3) break;
              str = vtmp[0] + "getc(readerfp)" + vtmp[2];
            }
            for(;;){
              vtmp = findFunction(str, "getchar()");
              if(vtmp.size()!=3) break;
              str = vtmp[0] + "getc(readerfp)" + vtmp[2];
            }
            ifun.func[ifun.name[i]] = str;
          }
        }
      }

      vs = findFunction(tmp, "writerFile()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"writerFile");
        ifun.already.insert((string)"my_putchar_unlocked");
        if(!writerfile){
          int i;
          writerfile = 1;
          rep(i, ifun.name.size()) if(ifun.name[i].substr(0,7)=="writer_"){
            string str = ifun.func[ifun.name[i]];
            for(;;){
              vtmp = findFunction(str, "my_putchar_unlocked()");
              if(vtmp.size()!=3) break;
              str = vtmp[0] + "putc( " + vtmp[1] + ", writerfp)" + vtmp[2];
            }
            for(;;){
              vtmp = findFunction(str, "putchar_unlocked()");
              if(vtmp.size()!=3) break;
              str = vtmp[0] + "putc( " + vtmp[1] + ", writerfp)" + vtmp[2];
            }
            for(;;){
              vtmp = findFunction(str, "my_putchar()");
              if(vtmp.size()!=3) break;
              str = vtmp[0] + "putc( " + vtmp[1] + ", writerfp)" + vtmp[2];
            }
            for(;;){
              vtmp = findFunction(str, "putchar()");
              if(vtmp.size()!=3) break;
              str = vtmp[0] + "putc( " + vtmp[1] + ", writerfp)" + vtmp[2];
            }
            ifun.func[ifun.name[i]] = str;
          }
        }
      }

      vs = findFunction(tmp, "b[]()");
      if(vs.size()==4){
        int i, j;
        string s;
        vector<string> v1, v2;
        v1 = split_p(vs[1], ',');
        v2 = split_p(vs[2], ',');

        rep(i,v2.size()){
          if(i){
            j = i;
            if(j >= v1.size()) j = v1.size() - 1;
            s = "(" + s + "*(" + v1[j] + "))";
          }

          if(i){
            s = "(" + s + "+(" + v2[i] + "))";
          } else {
            s = "(" + v2[i] + ")";
          }
        }

        tmp = vs[0] + s + vs[3];
        continue;
      }

      vs = findFunction(tmp, "rdLine()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"rdLine");
        tmp = vs[0] + "rdLine_L(" + vs[1] + ")" + vs[2];
        continue;
      }

      vs = findFunction(tmp, "sortF()");
      if(vs.size() == 3){ // 手抜き
        vtmp = split_p(vs[1], ',');
        ifun.doit.insert((string)"sortF_int");
        ifun.doit.insert((string)"sortF_ll");
        tmp = vs[0] + "sortF_L(" + vs[1] + ")" + vs[2];
        continue;
      }
      
      vs = findFunction(tmp, "fib_mod()");
      if(vs.size()!=3) vs = findFunction(tmp, "fibonacci_mod()");
      if(vs.size() == 3){
        vtmp = split_p(vs[1], ',');
        if(vtmp.size()==1){
          vs[1] = vs[1] + ", MD";
          ifun.doit.insert((string)"define_MD");
        }
        ifun.doit.insert((string)"fibonacci_mod");
        tmp = vs[0] + "fibonacci_mod_L(" + vs[1] + ")" + vs[2];
        continue;
      }
      
      vs = findFunction(tmp, "coordcomp()");
      if(vs.size() != 3) vs = findFunction(tmp, "coord_comp()");
      if(vs.size() == 3){
        trim(vs[0]);
        if(vs[0].size() == 0 || vs[0][vs[0].size()-1] != '.'){
          vtmp = split_p(vs[1], ',');
          if(2 <= vtmp.size() && vtmp.size() <= 4) ifun.doit.insert((string)"coordcomp_1");
          if(4 <= vtmp.size() && vtmp.size() <= 7) ifun.doit.insert((string)"coordcomp_2");
          if(6 <= vtmp.size() && vtmp.size() <= 10) ifun.doit.insert((string)"coordcomp_3");
          tmp = vs[0] + "coordcomp_L(" + vs[1] + ")" + vs[2];
          continue;
        }
      }

      vs = findFunction(tmp, "Digit()");
      if(vs.size() == 3){
        vtmp = split_p(vs[1], ',');
        if(vtmp.size() == 1) ifun.doit.insert((string)"Digit");
        if(vtmp.size() == 2) ifun.doit.insert((string)"Digit_base");
        if(vtmp.size() == 2) ifun.doit.insert((string)"Digit_all");
        if(vtmp.size() == 3) ifun.doit.insert((string)"Digit_all_base");
        tmp = vs[0] + "Digit_L(" + vs[1] + ")" + vs[2];
        continue;
      }

      vs = findFunction(tmp, "DigitF()");
      if(vs.size() == 3){
        vtmp = split_p(vs[1], ',');
        if(vtmp.size() == 3) ifun.doit.insert((string)"DigitF_all");
        if(vtmp.size() == 4) ifun.doit.insert((string)"DigitF_all_base");
        tmp = vs[0] + "DigitF_L(" + vs[1] + ")" + vs[2];
        continue;
      }

      vs = findFunction(tmp, "invDigit()");
      if(vs.size() == 3){
        vtmp = split_p(vs[1], ',');
        if(vtmp.size() == 2) ifun.doit.insert((string)"invDigit");
        if(vtmp.size() == 3) ifun.doit.insert((string)"invDigit_base");
        tmp = vs[0] + "invDigit_L(" + vs[1] + ")" + vs[2];
        continue;
      }

      vs = findFunction(tmp, "sod()");
      if(vs.size() == 3){
        vtmp = split_p(vs[1], ',');
        if(vtmp.size() == 1) ifun.doit.insert((string)"sod");
        if(vtmp.size() == 2) ifun.doit.insert((string)"sod_base");
        tmp = vs[0] + "sod_L(" + vs[1] + ")" + vs[2];
        continue;
      }

      vs = findFunction(tmp, "wAdjEdge()");
      if(vs.size() == 3){
        vtmp = split_p(vs[1], ',');
        if(vtmp.size() == 6 || vtmp.size() == 7) ifun.doit.insert((string)"wAdjEdge1");
        else if(vtmp.size() == 8 || vtmp.size() == 9) ifun.doit.insert((string)"wAdjEdge2");
        else if(vtmp.size() == 10 || vtmp.size() == 11) ifun.doit.insert((string)"wAdjEdge3");
        else if(vtmp.size() == 12 || vtmp.size() == 13) ifun.doit.insert((string)"wAdjEdge4");
        else{
          fprintf(stderr, "#arg of wAdjEdge is invalid\n");
          assert(0);
        }
        tmp = vs[0] + "wAdjEdge_L(" + vs[1] + ")" + vs[2];
        continue;
      }

      vs = findFunction(tmp, "convolution()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"convolution");
        ifun.doit.insert((string)"convolution-mint");
        ifun.doit.insert((string)"convolution-Mint");
        ifun.doit.insert((string)"convolution-modint");
        ifun.doit.insert((string)"convolution-Modint");
        tmp = vs[0] + "convolution_L(" + vs[1] + ")" + vs[2];
        continue;
      }

      vs = findFunction(tmp, "ZetaTransform()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"ZetaTransform");
        tmp = vs[0] + "ZetaTransform_L(" + vs[1] + ")" + vs[2];
        continue;
      }
      
      vs = findFunction(tmp, "ZetaTransform_min()");
      if(vs.size() == 3){
        vtmp = split_p(vs[1], ',');
        if(vtmp.size() == 2 || vtmp.size() == 3) ifun.doit.insert((string)"ZetaTransform_min");
        if(vtmp.size() == 4) ifun.doit.insert((string)"ZetaTransform_min2");
        if(vtmp.size() == 5) ifun.doit.insert((string)"ZetaTransform_min3");
        tmp = vs[0] + "ZetaTransform_min_L(" + vs[1] + ")" + vs[2];
        continue;
      }
      
      vs = findFunction(tmp, "ZetaTransform_max()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"ZetaTransform_max");
        tmp = vs[0] + "ZetaTransform_max_L(" + vs[1] + ")" + vs[2];
        continue;
      }

      vs = findFunction(tmp, "polationPoly()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"polationPoly");
        tmp = vs[0] + "polationPoly_L(" + vs[1] + ")" + vs[2];
        continue;
      }
      
      vs = findFunction(tmp, "Explode()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"Explode");
        tmp = vs[0] + "Explode_L(" + vs[1] + ")" + vs[2];
        continue;
      }
      
      vs = findFunction(tmp, "Implode()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"Implode");
        tmp = vs[0] + "Implode_L(" + vs[1] + ")" + vs[2];
        continue;
      }
      
      vs = findFunction(tmp, "crossProd()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"crossProd");
        tmp = vs[0] + "crossProd_L(" + vs[1] + ")" + vs[2];
        continue;
      }
      
      vs = findFunction(tmp, "LineIntersection_size()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"LineIntersection_size");
        tmp = vs[0] + "LineIntersection_size_L(" + vs[1] + ")" + vs[2];
        continue;
      }
      
      vs = findFunction(tmp, "arrEraseVal()");
      if(vs.size() == 3){
        vtmp = split_p(vs[1], ',');
        if(vtmp.size() == 3) ifun.doit.insert((string)"arrEraseVal1");
        if(vtmp.size() == 4) ifun.doit.insert((string)"arrEraseVal2");
        if(vtmp.size() == 5) ifun.doit.insert((string)"arrEraseVal3");
        if(vtmp.size() == 6) ifun.doit.insert((string)"arrEraseVal4");
        tmp = vs[0] + "arrEraseVal_L(" + vs[1] + ")" + vs[2];
        continue;
      }
      
      vs = findFunction(tmp, "vecEraseVal()");
      if(vs.size() == 3){
        vtmp = split_p(vs[1], ',');
        if(vtmp.size() == 2) ifun.doit.insert((string)"vecEraseVal1");
        tmp = vs[0] + "vecEraseVal_L(" + vs[1] + ")" + vs[2];
        continue;
      }
      
      vs = findFunction(tmp, "DistinctE()");
      if(vs.size() == 3){
        vtmp = split_p(vs[1], ',');
        if(vtmp.size() == 2) ifun.doit.insert((string)"DistinctE_2");
        if(vtmp.size() == 3) ifun.doit.insert((string)"DistinctE_3");
        if(vtmp.size() == 4) ifun.doit.insert((string)"DistinctE_4");
        tmp = vs[0] + "DistinctE_L(" + vs[1] + ")" + vs[2];
        continue;
      }
      
      vs = findFunction(tmp, "Mex()");
      if(vs.size() == 3){
        vtmp = split_p(vs[1], ',');
        if(vtmp.size() == 1) ifun.doit.insert((string)"Mex1");
        if(vtmp.size() == 2) ifun.doit.insert((string)"Mex2"), ifun.doit.insert((string)"Mex");
        if(vtmp.size() == 3) ifun.doit.insert((string)"Mex3"), ifun.doit.insert((string)"Mex");
        if(vtmp.size() == 4) ifun.doit.insert((string)"Mex");
        tmp = vs[0] + "Mex_L(" + vs[1] + ")" + vs[2];
        continue;
      }
      
      vs = findFunction(tmp, "Kth0()");
      if(vs.size() == 3){
        vtmp = split_p(vs[1], ',');
        if(vtmp.size() == 2) ifun.doit.insert((string)"Kth0_size2");
        if(vtmp.size() == 3) ifun.doit.insert((string)"Kth0_size3");
        if(vtmp.size() == 4) ifun.doit.insert((string)"Kth0_size4");
        tmp = vs[0] + "Kth0_L(" + vs[1] + ")" + vs[2];
        continue;
      }
      
      vs = findFunction(tmp, "Kth1()");
      if(vs.size() == 3){
        vtmp = split_p(vs[1], ',');
        if(vtmp.size() == 2) ifun.doit.insert((string)"Kth1_size2");
        if(vtmp.size() == 3) ifun.doit.insert((string)"Kth1_size3");
        if(vtmp.size() == 4) ifun.doit.insert((string)"Kth1_size4");
        tmp = vs[0] + "Kth1_L(" + vs[1] + ")" + vs[2];
        continue;
      }
      
      vs = findFunction(tmp, "Kth2()");
      if(vs.size() == 3){
        vtmp = split_p(vs[1], ',');
        if(vtmp.size() == 3) ifun.doit.insert((string)"Kth2_size3");
        if(vtmp.size() == 4) ifun.doit.insert((string)"Kth2_size4");
        tmp = vs[0] + "Kth2_L(" + vs[1] + ")" + vs[2];
        continue;
      }
      
      vs = findFunction(tmp, "Kth3()");
      if(vs.size() == 3){
        vtmp = split_p(vs[1], ',');
        if(vtmp.size() == 4) ifun.doit.insert((string)"Kth3_size4");
        tmp = vs[0] + "Kth3_L(" + vs[1] + ")" + vs[2];
        continue;
      }
      
      vs = findFunction(tmp, "KthA()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"KthA");
        tmp = vs[0] + "KthA_L(" + vs[1] + ")" + vs[2];
        continue;
      }
      
      vs = findFunction(tmp, "strReplace()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"strReplace1");
        ifun.doit.insert((string)"strReplace2");
        tmp = vs[0] + "strReplace_L(" + vs[1] + ")" + vs[2];
        continue;
      }
      
      vs = findFunction(tmp, "BIT_popcount()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"BIT_popcount");
        tmp = vs[0] + "BIT_popcount_L(" + vs[1] + ")" + vs[2];
        continue;
      }
      
      vs = findFunction(tmp, "BIT_ctz()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"BIT_ctz");
        tmp = vs[0] + "BIT_ctz_L(" + vs[1] + ")" + vs[2];
        continue;
      }
      
      vs = findFunction(tmp, "BIT_Ith()");
      if(vs.size() == 3){
        vtmp = split_p(vs[1], ',');
        if(vtmp.size()==2){
          trim(vtmp[0]);
          trim(vtmp[1]);
          tmp = vs[0] + "(((" + vtmp[0] + ") >> (" + vtmp[1] + "))&1)" + vs[2];
        } else {
          trim(vs[1]);
          tmp = vs[0] + "(1LL<<(" + vs[1] + "))" + vs[2];
        }
        continue;
      }

      vs = findFunction(tmp, "BIT_ith()");
      if(vs.size() == 3){
        vtmp = split_p(vs[1], ',');
        if(vtmp.size()==2){
          trim(vtmp[0]);
          trim(vtmp[1]);
          tmp = vs[0] + "((" + vtmp[0] + ") & BIT_ith(" + vtmp[1] + "))" + vs[2];
        } else {
          trim(vs[1]);
          if(isdigit(vs[1][0]) && atoi(vs[1].c_str()) >= 31) tmp = vs[0] + "(1LL<<(" + vs[1] + "))" + vs[2];
          else                                               tmp = vs[0] + "(1<<(" + vs[1] + "))" + vs[2];
        }
        continue;
      }

      vs = findFunction(tmp, "BIT_lowest()");
      if(vs.size() == 3){
        tmp = vs[0] + "(-(" + vs[1] + ") & (" + vs[1] + "))" + vs[2];
        continue;
      }

      vs = findFunction(tmp, "BIT_nonlowest()");
      if(vs.size() == 3){
        tmp = vs[0] + "((" + vs[1] + ") & ((" + vs[1] + ")-1))" + vs[2];
        continue;
      }

      vs = findFunction(tmp, "Ilog2_f()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"Ilog2_f");
        tmp = vs[0] + "Ilog2_f_L(" + vs[1] + ")" + vs[2];
        continue;
      }

      vs = findFunction(tmp, "Ilog2_c()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"Ilog2_c");
        tmp = vs[0] + "Ilog2_c_L(" + vs[1] + ")" + vs[2];
        continue;
      }

      vs = findFunction(tmp, "Ilog2_s()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"Ilog2_s");
        tmp = vs[0] + "Ilog2_s_L(" + vs[1] + ")" + vs[2];
        continue;
      }

      vs = findFunction(tmp, "Prime()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"Prime");
        tmp = vs[0] + "Prime_L(" + vs[1] + ")" + vs[2];
        continue;
      }

      vs = findFunction(tmp, "Isqrt_f()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"Isqrt_f");
        tmp = vs[0] + "Isqrt_f_L(" + vs[1] + ")" + vs[2];
        continue;
      }

      vs = findFunction(tmp, "Isqrt_c()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"Isqrt_c");
        tmp = vs[0] + "Isqrt_c_L(" + vs[1] + ")" + vs[2];
        continue;
      }

      vs = findFunction(tmp, "Isqrt_s()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"Isqrt_s");
        tmp = vs[0] + "Isqrt_s_L(" + vs[1] + ")" + vs[2];
        continue;
      }

      vs = findFunction(tmp, "Icbrt_f()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"Icbrt_f");
        tmp = vs[0] + "Icbrt_f_L(" + vs[1] + ")" + vs[2];
        continue;
      }

      vs = findFunction(tmp, "Icbrt_c()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"Icbrt_c");
        tmp = vs[0] + "Icbrt_c_L(" + vs[1] + ")" + vs[2];
        continue;
      }

      vs = findFunction(tmp, "Icbrt_s()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"Icbrt_s");
        tmp = vs[0] + "Icbrt_s_L(" + vs[1] + ")" + vs[2];
        continue;
      }

      vs = findFunction(tmp, "Iroot_f()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"Iroot_f");
        tmp = vs[0] + "Iroot_f_L(" + vs[1] + ")" + vs[2];
        continue;
      }

      vs = findFunction(tmp, "Iroot_c()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"Iroot_c");
        tmp = vs[0] + "Iroot_c_L(" + vs[1] + ")" + vs[2];
        continue;
      }

      vs = findFunction(tmp, "Iroot_s()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"Iroot_s");
        tmp = vs[0] + "Iroot_s_L(" + vs[1] + ")" + vs[2];
        continue;
      }


      vs = findFunction(tmp, "Isqrt()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"Isqrt_f");
        tmp = vs[0] + "Isqrt_f_L(" + vs[1] + ")" + vs[2];
        continue;
      }

      vs = findFunction(tmp, "Icbrt()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"Icbrt_f");
        tmp = vs[0] + "Icbrt_f_L(" + vs[1] + ")" + vs[2];
        continue;
      }

      vs = findFunction(tmp, "Iroot()");
      if(vs.size() == 3){
        ifun.doit.insert((string)"Iroot_f");
        tmp = vs[0] + "Iroot_f_L(" + vs[1] + ")" + vs[2];
        continue;
      }



      vs = findFunction(tmp, "primitiveRoot()");
      if(vs.size() == 3) ifun.doit.insert((string)"primitiveRoot");
      
      vs = findFunction(tmp, "knightDistance()");
      if(vs.size() == 3) ifun.doit.insert((string)"knightDistance");
      
      vs = findFunction(tmp, "LIS_length()");
      if(vs.size() == 3) ifun.doit.insert((string)"LIS_length");
      
      vs = findFunction(tmp, "weaklyLIS_length()");
      if(vs.size() == 3) ifun.doit.insert((string)"weaklyLIS_length");
      
      vs = findFunction(tmp, "LIS_ends()");
      if(vs.size() == 3) ifun.doit.insert((string)"LIS_ends");
      
      vs = findFunction(tmp, "weaklyLIS_ends()");
      if(vs.size() == 3) ifun.doit.insert((string)"weaklyLIS_ends");
      
      vs = findFunction(tmp, "longestSuffixPrefix()");
      if(vs.size() == 3) ifun.doit.insert((string)"longestSuffixPrefix");
      
      vs = findFunction(tmp, "SuffixArray()");
      if(vs.size() == 3) ifun.doit.insert((string)"SuffixArray");
      
      vs = findFunction(tmp, "smallestSubsequenceLengthK()");
      if(vs.size() == 3) ifun.doit.insert((string)"smallestSubsequenceLengthK");
      
      vs = findFunction(tmp, "isSubstring()");
      if(vs.size() == 3) ifun.doit.insert((string)"isSubstring_string");
      
      vs = findFunction(tmp, "isPalindrome()");
      if(vs.size() == 3) ifun.doit.insert((string)"isPalindrome");
      
      vs = findFunction(tmp, "rd_int()");
      if(vs.size() == 3) ifun.doit.insert((string)"rd_int");
      
      vs = findFunction(tmp, "rd_ll()");
      if(vs.size() == 3) ifun.doit.insert((string)"rd_ll");
      
      vs = findFunction(tmp, "rd_string()");
      if(vs.size() == 3) ifun.doit.insert((string)"rd_string");
      
      vs = findFunction(tmp, "fDiv()");
      if(vs.size() == 3) ifun.doit.insert((string)"fDiv");
      
      vs = findFunction(tmp, "cDiv()");
      if(vs.size() == 3) ifun.doit.insert((string)"cDiv");
      
      vs = findFunction(tmp, "floor_sum()");
      if(vs.size() == 3) ifun.doit.insert((string)"floor_sum");
      
      vs = findFunction(tmp, "floor_sum2()");
      if(vs.size() == 3) ifun.doit.insert((string)"floor_sum2");
      
      vs = findFunction(tmp, "graph_minColor()");
      if(vs.size() == 3) ifun.doit.insert((string)"graph_minColor");
      
      vs = findFunction(tmp, "xorMin()");
      if(vs.size() == 3) ifun.doit.insert((string)"xorMin");
      vs = findFunction(tmp, "xorMax()");
      if(vs.size() == 3) ifun.doit.insert((string)"xorMax");
      
      vs = findFunction(tmp, "calcRollingHash()");
      if(vs.size() == 3) ifun.doit.insert((string)"rollingHash");
      
      vs = findFunction(tmp, "fft()");
      if(vs.size() == 3) ifun.doit.insert((string)"fft");
      if(vs.size() == 3) ifun.doit.insert((string)"fft-mint");
      if(vs.size() == 3) ifun.doit.insert((string)"fft-Mint");
      if(vs.size() == 3) ifun.doit.insert((string)"fft-modint");
      if(vs.size() == 3) ifun.doit.insert((string)"fft-Modint");
      
      vs = findFunction(tmp, "get_inv_mod()");
      if(vs.size() == 3) ifun.doit.insert((string)"get_inv_mod");
      
      vs = findFunction(tmp, "extendedEuclid()");
      if(vs.size() == 3) ifun.doit.insert((string)"extendedEuclid");
      
      vs = findFunction(tmp, "chineseRemainder()");
      if(vs.size() == 3) ifun.doit.insert((string)"chineseRemainder");
      
      vs = findFunction(tmp, "writerDigit_double()");
      if(vs.size() == 3) ifun.doit.insert((string)"writer_double");
      
      vs = findFunction(tmp, "arrcmp()");
      if(vs.size() == 3) ifun.doit.insert((string)"arrcmp");
      
      vs = findFunction(tmp, "arrRot()");
      if(vs.size() == 3) ifun.doit.insert((string)"arrRot");
      
      vs = findFunction(tmp, "arrInsert()");
      if(vs.size() == 3) ifun.doit.insert((string)"arrInsert");
      
      vs = findFunction(tmp, "reduceFraction()");
      if(vs.size() == 3) ifun.doit.insert((string)"reduceFraction");

      vs = findFunction(tmp, "MoebiusTransform()");
      if(vs.size() == 3) ifun.doit.insert((string)"MoebiusTransform");

      vs = findFunction(tmp, "isVowel()");
      if(vs.size() == 3) ifun.doit.insert((string)"isVowel");

      vs = findFunction(tmp, "KMP()");
      if(vs.size() == 3) ifun.doit.insert((string)"KMP");

      vs = findFunction(tmp, "Hungarian()");
      if(vs.size() == 3) ifun.doit.insert((string)"Hungarian");

      vs = findFunction(tmp, "counterSumIsLT()");
      if(vs.size() == 3) ifun.doit.insert((string)"counterSumIsLT");

      vs = findFunction(tmp, "counterProdIsLT()");
      if(vs.size() == 3) ifun.doit.insert((string)"counterProdIsLT");

      vs = findFunction(tmp, "counterD2SumIsLT()");
      if(vs.size() == 3) ifun.doit.insert((string)"counterD2SumIsLT");

      vs = findFunction(tmp, "counterM2SumIsLT()");
      if(vs.size() == 3) ifun.doit.insert((string)"counterM2SumIsLT");

      vs = findFunction(tmp, "counterD2ProdIsLT()");
      if(vs.size() == 3) ifun.doit.insert((string)"counterD2ProdIsLT");

      vs = findFunction(tmp, "counterM2ProdIsLT()");
      if(vs.size() == 3) ifun.doit.insert((string)"counterM2ProdIsLT");

      vs = findFunction(tmp, "HammingDistance()");
      if(vs.size() == 3) ifun.doit.insert((string)"HammingDistance");

      vs = findFunction(tmp, "editDistance()");
      if(vs.size() == 3) ifun.doit.insert((string)"editDistance");

      vs = findFunction(tmp, "Cmod2()");
      if(vs.size() == 3) ifun.doit.insert((string)"Cmod2");

      vs = findFunction(tmp, "next_mcomb()");
      if(vs.size() == 3) ifun.doit.insert((string)"next_mcomb");

      vs = findFunction(tmp, "next_scomb()");
      if(vs.size() == 3) ifun.doit.insert((string)"next_scomb");

      vs = findFunction(tmp, "next_marr()");
      if(vs.size() == 3) ifun.doit.insert((string)"next_marr");

      vs = findFunction(tmp, "next_sarr()");
      if(vs.size() == 3) ifun.doit.insert((string)"next_sarr");

      vs = findFunction(tmp, "next_sarr_s()");
      if(vs.size() == 3) ifun.doit.insert((string)"next_sarr_s");

      vs = findFunction(tmp, "walloc1d()");
      if(vs.size() == 3) ifun.doit.insert((string)"walloc1d");

      vs = findFunction(tmp, "walloc2d()");
      if(vs.size() == 3) ifun.doit.insert((string)"walloc2d");

      vs = findFunction(tmp, "malloc1d()");
      if(vs.size() == 3) ifun.doit.insert((string)"malloc1d");

      vs = findFunction(tmp, "free1d()");
      if(vs.size() == 3) ifun.doit.insert((string)"free1d");

      vs = findFunction(tmp, "malloc2d()");
      if(vs.size() == 3) ifun.doit.insert((string)"malloc2d");

      vs = findFunction(tmp, "free2d()");
      if(vs.size() == 3) ifun.doit.insert((string)"free2d");

      vs = findFunction(tmp, "isPrime()");
      if(vs.size() == 3) ifun.doit.insert((string)"isPrime");

      vs = findFunction(tmp, "Factor()");
      if(vs.size() == 3) ifun.doit.insert((string)"Factor");

      vs = findFunction(tmp, "FactorM()");
      if(vs.size() == 3) ifun.doit.insert((string)"FactorM");

      vs = findFunction(tmp, "Divisor()");
      if(vs.size() == 3) ifun.doit.insert((string)"Divisor");

      vs = findFunction(tmp, "DivisorSum()");
      if(vs.size() == 3) ifun.doit.insert((string)"DivisorSum");

      vs = findFunction(tmp, "Moebius()");
      if(vs.size() == 3) ifun.doit.insert((string)"Moebius");

      vs = findFunction(tmp, "EulerPhi()");
      if(vs.size() == 3) ifun.doit.insert((string)"EulerPhi");

      vs = findFunction(tmp, "sortI()");
      if(vs.size() == 3) ifun.doit.insert((string)"sortI");

      vs = findFunction(tmp, "sortA()");
      if(vs.size() == 3) ifun.doit.insert((string)"sortA");

      vs = findFunction(tmp, "rsortA()");
      if(vs.size() == 3) ifun.doit.insert((string)"rsortA");

      vs = findFunction(tmp, "sortA_index()");
      if(vs.size() == 3) ifun.doit.insert((string)"sortA_index");

      vs = findFunction(tmp, "sortV()");
      if(vs.size() == 3) ifun.doit.insert((string)"sortV");

      vs = findFunction(tmp, "rsortV()");
      if(vs.size() == 3) ifun.doit.insert((string)"rsortV");

      vs = findFunction(tmp, "BIT_parity()");
      if(vs.size() == 3) ifun.doit.insert((string)"BIT_parity");

      vs = findFunction(tmp, "BIT_parity_pm()");
      if(vs.size() == 3) ifun.doit.insert((string)"BIT_parity_pm");

      vs = findFunction(tmp, "prodDigits()");
      if(vs.size() == 3) ifun.doit.insert((string)"prodDigits");

      vs = findFunction(tmp, "DigitHist()");
      if(vs.size() == 3) ifun.doit.insert((string)"DigitHist");

      vs = findFunction(tmp, "invDigit_r()");
      if(vs.size() == 3) ifun.doit.insert((string)"invDigit_r");

      vs = findFunction(tmp, "STR2int()");
      if(vs.size() == 3) ifun.doit.insert((string)"STR2int");

      vs = findFunction(tmp, "STR2ll()");
      if(vs.size() == 3) ifun.doit.insert((string)"STR2ll");

      vs = findFunction(tmp, "intervalSieve()");
      if(vs.size() == 3) ifun.doit.insert((string)"intervalSieve");

      vs = findFunction(tmp, "powmod()");
      if(vs.size() == 3) ifun.doit.insert((string)"powmod");

      vs = findFunction(tmp, "PowMod()");
      if(vs.size() == 3) ifun.doit.insert((string)"PowMod");

      vs = findFunction(tmp, "runLength()");
      if(vs.size() == 3) ifun.doit.insert((string)"runLength");

      vs = findFunction(tmp, "arrErase()");
      if(vs.size() == 3) ifun.doit.insert((string)"arrErase");

      vs = findFunction(tmp, "Distinct()");
      if(vs.size() == 3) ifun.doit.insert((string)"Distinct");

      vs = findFunction(tmp, "Count()");
      if(vs.size() == 3) ifun.doit.insert((string)"Count");

      vs = findFunction(tmp, "arrCountVal()");
      if(vs.size() == 3) ifun.doit.insert((string)"arrCountVal");

      vs = findFunction(tmp, "arrCountValSeqMax()");
      if(vs.size() == 3) ifun.doit.insert((string)"arrCountValSeqMax");

      vs = findFunction(tmp, "PalindromeCost()");
      if(vs.size() == 3) ifun.doit.insert((string)"PalindromeCost");

      vs = findFunction(tmp, "Fibonacci_mod()");
      if(vs.size() == 3) ifun.doit.insert((string)"Fibonacci_mod");

      vs = findFunction(tmp, "Fib_mod()");
      if(vs.size() == 3) ifun.doit.insert((string)"Fib_mod");

      vs = findFunction(tmp, "Cmod2()");
      if(vs.size() == 3) ifun.doit.insert((string)"Cmod2");

      vs = findFunction(tmp, "Cmod5()");
      if(vs.size() == 3) ifun.doit.insert((string)"Cmod5");

      vs = findFunction(tmp, "Cmod10()");
      if(vs.size() == 3) ifun.doit.insert((string)"Cmod10");

      vs = findFunction(tmp, "Unique()");
      if(vs.size() == 3) ifun.doit.insert((string)"Unique");

      vs = findFunction(tmp, "isSubsequence()");
      if(vs.size() == 3) ifun.doit.insert((string)"isSubsequence");

      vs = findFunction(tmp, "isSubsequence_r()");
      if(vs.size() == 3) ifun.doit.insert((string)"isSubsequence_r");

      vs = findFunction(tmp, "inversion_range()");
      if(vs.size() == 3) ifun.doit.insert((string)"inversion_range");

      vs = findFunction(tmp, "inversion()");
      if(vs.size() == 3) ifun.doit.insert((string)"inversion");

      vs = findFunction(tmp, "ZetaTransform2()");
      if(vs.size() == 3) ifun.doit.insert((string)"ZetaTransform2");

      vs = findFunction(tmp, "MoebiusTransform2()");
      if(vs.size() == 3) ifun.doit.insert((string)"MoebiusTransform2");

      vs = findFunction(tmp, "HadamardTransform()");
      if(vs.size() == 3) ifun.doit.insert((string)"HadamardTransform");

      vs = findFunction(tmp, "xorConvolution()");
      if(vs.size() == 3) ifun.doit.insert((string)"xorConvolution");

      vs = findFunction(tmp, "orConvolution()");
      if(vs.size() == 3) ifun.doit.insert((string)"orConvolution");

      vs = findFunction(tmp, "andConvolution()");
      if(vs.size() == 3) ifun.doit.insert((string)"andConvolution");

      vs = findFunction(tmp, "slideMin()");
      if(vs.size() == 3) ifun.doit.insert((string)"slideMin");

      vs = findFunction(tmp, "slideMax()");
      if(vs.size() == 3) ifun.doit.insert((string)"slideMax");

      vs = findFunction(tmp, "isLeapYear()");
      if(vs.size() == 3) ifun.doit.insert((string)"isLeapYear");

      vs = findFunction(tmp, "dayOfWeek()");
      if(vs.size() == 3) ifun.doit.insert((string)"dayOfWeek");

      vs = findFunction(tmp, "dayOfWeekStr()");
      if(vs.size() == 3) ifun.doit.insert((string)"dayOfWeekStr");

      vs = findFunction(tmp, "prevDay()");
      if(vs.size() == 3) ifun.doit.insert((string)"prevDay");

      vs = findFunction(tmp, "nextDay()");
      if(vs.size() == 3) ifun.doit.insert((string)"nextDay");

      vs = findFunction(tmp, "dayIndex()");
      if(vs.size() == 3) ifun.doit.insert((string)"dayIndex");

      vs = findFunction(tmp, "dayFromIndex()");
      if(vs.size() == 3) ifun.doit.insert((string)"dayFromIndex");

      vs = findFunction(tmp, "polationVal()");
      if(vs.size() == 3) ifun.doit.insert((string)"polationVal");

      vs = findFunction(tmp, "TSP_cycle()");
      if(vs.size() == 3) ifun.doit.insert((string)"TSP_cycle");

      vs = findFunction(tmp, "TSP_path()");
      if(vs.size() == 3) ifun.doit.insert((string)"TSP_path");

      vs = findFunction(tmp, "TSP_path_s()");
      if(vs.size() == 3) ifun.doit.insert((string)"TSP_path_s");

      vs = findFunction(tmp, "maxSubsetDP()");
      if(vs.size() == 3) ifun.doit.insert((string)"maxSubsetDP");

      vs = findFunction(tmp, "maxRectArea()");
      if(vs.size() == 3) ifun.doit.insert((string)"maxRectArea");

      vs = findFunction(tmp, "isValidBracket1()");
      if(vs.size() == 3) ifun.doit.insert((string)"isValidBracket1");

      vs = findFunction(tmp, "isValidBracket2()");
      if(vs.size() == 3) ifun.doit.insert((string)"isValidBracket2");

      vs = findFunction(tmp, "swapV()");
      if(vs.size() == 3) ifun.doit.insert((string)"swapV");

      vs = findFunction(tmp, "cntSubarrayFreq()");
      if(vs.size() == 3) ifun.doit.insert((string)"cntSubarrayFreq");

      vs = findFunction(tmp, "cntSubarrayDistinct()");
      if(vs.size() == 3) ifun.doit.insert((string)"cntSubarrayDistinct");

      vs = findFunction(tmp, "maxSubarrayDistinct()");
      if(vs.size() == 3) ifun.doit.insert((string)"maxSubarrayDistinct");

      vs = findFunction(tmp, "vec2arr()");
      if(vs.size() == 3) ifun.doit.insert((string)"vec2arr");

      vs = findFunction(tmp, "Rot90()");
      if(vs.size() == 3) ifun.doit.insert((string)"Rot90");

      vs = findFunction(tmp, "minFactorList()");
      if(vs.size() == 3) ifun.doit.insert((string)"minFactorList");

      vs = findFunction(tmp, "maxFactorList()");
      if(vs.size() == 3) ifun.doit.insert((string)"maxFactorList");

      vs = findFunction(tmp, "FactorList()");
      if(vs.size() == 3) ifun.doit.insert((string)"FactorList");

      vs = findFunction(tmp, "FactorMList()");
      if(vs.size() == 3) ifun.doit.insert((string)"FactorMList");

      vs = findFunction(tmp, "EulerPhiList()");
      if(vs.size() == 3) ifun.doit.insert((string)"EulerPhiList");

      vs = findFunction(tmp, "LinearEquationMod2()");
      if(vs.size() == 3) ifun.doit.insert((string)"LinearEquationMod2");

      vs = findFunction(tmp, "LinearEquation()");
      if(vs.size() == 3) ifun.doit.insert((string)"LinearEquation");

      vs = findFunction(tmp, "LexicographicGE()");
      if(vs.size() == 3) ifun.doit.insert((string)"LexicographicGE");

      vs = findFunction(tmp, "cReader_ll()");
      if(vs.size() == 3) ifun.doit.insert((string)"cReader_ll");

      vs = findFunction(tmp, "cReader_eof()");
      if(vs.size() == 3) ifun.doit.insert((string)"cReader_eof");

      vs = findFunction(tmp, "cntPrime()");
      if(vs.size() == 3) ifun.doit.insert((string)"cntPrime");

      vs = findFunction(tmp, "arrMerge()");
      if(vs.size() == 3) ifun.doit.insert((string)"arrMerge");

      vs = findFunction(tmp, "arrMergeD()");
      if(vs.size() == 3) ifun.doit.insert((string)"arrMergeD");

      vs = findFunction(tmp, "cntArrayNecessaryElement()");
      if(vs.size() == 3) ifun.doit.insert((string)"cntArrayNecessaryElement");

      vs = findFunction(tmp, "cntArrayNecessaryElement_walloc()");
      if(vs.size() == 3) ifun.doit.insert((string)"cntArrayNecessaryElement_walloc");

      vs = findFunction(tmp, "arr2bit_int()");
      if(vs.size() == 3) ifun.doit.insert((string)"arr2bit_int");

      vs = findFunction(tmp, "arr2bit_ll()");
      if(vs.size() == 3) ifun.doit.insert((string)"arr2bit_ll");

      vs = findFunction(tmp, "toLower()");
      if(vs.size() == 3) ifun.doit.insert((string)"toLower");

      vs = findFunction(tmp, "isSorted()");
      if(vs.size() == 3) ifun.doit.insert((string)"isSorted");

      vs = findFunction(tmp, "Slice()");
      if(vs.size() == 3) ifun.doit.insert((string)"Slice");

      vs = findFunction(tmp, "RoundUp()");
      if(vs.size() == 3) ifun.doit.insert((string)"RoundUp");

      vs = findFunction(tmp, "RoundDown()");
      if(vs.size() == 3) ifun.doit.insert((string)"RoundDown");

      vs = findFunction(tmp, "startWith()");
      if(vs.size() == 3) ifun.doit.insert((string)"startWith");

      vs = findFunction(tmp, "endWith()");
      if(vs.size() == 3) ifun.doit.insert((string)"endWith");

      vs = findFunction(tmp, "MergeTech()");
      if(vs.size() == 3) ifun.doit.insert((string)"MergeTech");

      vs = findFunction(tmp, "WildEQ()");
      if(vs.size() == 3) ifun.doit.insert((string)"WildEQ");

      vs = findFunction(tmp, "kthPalindromicNumber64()");
      if(vs.size() == 3) ifun.doit.insert((string)"kthPalindromicNumber64");

      vs = findFunction(tmp, "reverseNumber()");
      if(vs.size() == 3) ifun.doit.insert((string)"reverseNumber");

      vs = findFunction(tmp, "isPalindromicNumber()");
      if(vs.size() == 3) ifun.doit.insert((string)"isPalindromicNumber");


      vs = findFunction(tmp, "popFirst()");
      if(vs.size() == 3) ifun.doit.insert((string)"multiset_popFirst");

      vs = findFunction(tmp, "getFirst()");
      if(vs.size() == 3) ifun.doit.insert((string)"multiset_getFirst");

      vs = findFunction(tmp, "popLast()");
      if(vs.size() == 3) ifun.doit.insert((string)"multiset_popLast");

      vs = findFunction(tmp, "getLast()");
      if(vs.size() == 3) ifun.doit.insert((string)"multiset_getLast");

      vs = findFunction(tmp, "popFirst()");
      if(vs.size() == 3) ifun.doit.insert((string)"set_popFirst");

      vs = findFunction(tmp, "getFirst()");
      if(vs.size() == 3) ifun.doit.insert((string)"set_getFirst");

      vs = findFunction(tmp, "popLast()");
      if(vs.size() == 3) ifun.doit.insert((string)"set_popLast");

      vs = findFunction(tmp, "getLast()");
      if(vs.size() == 3) ifun.doit.insert((string)"set_getLast");

      vs = findFunction(tmp, "getClosest()");
      if(vs.size() == 3) ifun.doit.insert((string)"multiset_getClosest");

      vs = findFunction(tmp, "getClosestL()");
      if(vs.size() == 3) ifun.doit.insert((string)"multiset_getClosestL");

      vs = findFunction(tmp, "getClosest()");
      if(vs.size() == 3) ifun.doit.insert((string)"set_getClosest");

      vs = findFunction(tmp, "getClosestL()");
      if(vs.size() == 3) ifun.doit.insert((string)"set_getClosestL");

      vs = findFunction(tmp, "getClosest()");
      if(vs.size() == 3) ifun.doit.insert((string)"arr_getClosest");

      vs = findFunction(tmp, "getClosestL()");
      if(vs.size() == 3) ifun.doit.insert((string)"arr_getClosestL");

      vs = findFunction(tmp, "setEdge()");
      if(vs.size() == 3) ifun.doit.insert((string)"graph_setEdge");

      vs = findFunction(tmp, "setDirectEdge()");
      if(vs.size() == 3) ifun.doit.insert((string)"graph_setDirectEdge");

      vs = findFunction(tmp, "setEdgeRootedTree()");
      if(vs.size() == 3) ifun.doit.insert((string)"graph_setEdgeRootedTree");

      vs = findFunction(tmp, "reverse()");
      if(vs.size() == 3) ifun.doit.insert((string)"graph_reverse");

      vs = findFunction(tmp, "reduce()");
      if(vs.size() == 3) ifun.doit.insert((string)"graph_reduce");

      vs = findFunction(tmp, "getDist()");
      if(vs.size() == 3) ifun.doit.insert((string)"graph_getDist");

      vs = findFunction(tmp, "TreeDiameter()");
      if(vs.size() == 3) ifun.doit.insert((string)"graph_TreeDiameter");

      vs = findFunction(tmp, "getDistTree_WeightedNode_max()");
      if(vs.size() == 3) ifun.doit.insert((string)"graph_getDistTree_WeightedNode_max");

      vs = findFunction(tmp, "getDistPairMatrix()");
      if(vs.size() == 3) ifun.doit.insert((string)"graph_getDistPairMatrix");

      vs = findFunction(tmp, "SubTreeSize()");
      if(vs.size() == 3) ifun.doit.insert((string)"graph_SubTreeSize");

      vs = findFunction(tmp, "SubTreeWeight()");
      if(vs.size() == 3) ifun.doit.insert((string)"graph_SubTreeWeight");

      vs = findFunction(tmp, "cntShortest()");
      if(vs.size() == 3) ifun.doit.insert((string)"graph_cntShortest");

      vs = findFunction(tmp, "scc()");
      if(vs.size() == 3) ifun.doit.insert((string)"graph_scc");

      vs = findFunction(tmp, "bcc()");
      if(vs.size() == 3) ifun.doit.insert((string)"graph_bcc");

      vs = findFunction(tmp, "articulation()");
      if(vs.size() == 3) ifun.doit.insert((string)"graph_articulation");

      vs = findFunction(tmp, "shortestPath()");
      if(vs.size() == 3) ifun.doit.insert((string)"graph_shortestPath");

      vs = findFunction(tmp, "TopologicalSort()");
      if(vs.size() == 3) ifun.doit.insert((string)"graph_TopologicalSort");

      vs = findFunction(tmp, "longestPath_length()");
      if(vs.size() == 3) ifun.doit.insert((string)"graph_longestPath_length");

      vs = findFunction(tmp, "Grundy()");
      if(vs.size() == 3) ifun.doit.insert((string)"graph_Grundy");

      vs = findFunction(tmp, "preorder()");
      if(vs.size() == 3) ifun.doit.insert((string)"graph_preorder");

      vs = findFunction(tmp, "anUndirectedCycle()");
      if(vs.size() == 3) ifun.doit.insert((string)"graph_anUndirectedCycle");

      vs = findFunction(tmp, "shortestCycle()");
      if(vs.size() == 3) ifun.doit.insert((string)"graph_shortestCycle");

      vs = findFunction(tmp, "shortestUndirectedCycle_length()");
      if(vs.size() == 3) ifun.doit.insert((string)"graph_shortestUndirectedCycle_length");

      vs = findFunction(tmp, "shortestUndirectedCycle_length()");
      if(vs.size() == 3) ifun.doit.insert((string)"graph_shortestUndirectedCycle_length_node");

      vs = findFunction(tmp, "maxIndependenceSet()");
      if(vs.size() == 3) ifun.doit.insert((string)"graph_maxIndependenceSet");

      vs = findFunction(tmp, "countIndependenceSet()");
      if(vs.size() == 3) ifun.doit.insert((string)"graph_countIndependenceSet");

      vs = findFunction(tmp, "bipartite()");
      if(vs.size() == 3) ifun.doit.insert((string)"graph_bipartite");

      vs = findFunction(tmp, "Rerooting()");
      if(vs.size() == 3) ifun.doit.insert((string)"graph_Rerooting");

      vs = findFunction(tmp, "setEdge()");
      if(vs.size() == 3) ifun.doit.insert((string)"wgraph_setEdge");

      vs = findFunction(tmp, "setDirectEdge()");
      if(vs.size() == 3) ifun.doit.insert((string)"wgraph_setDirectEdge");

      vs = findFunction(tmp, "setEdgeRootedTree()");
      if(vs.size() == 3) ifun.doit.insert((string)"wgraph_setEdgeRootedTree");

      vs = findFunction(tmp, "getDist()");
      if(vs.size() == 3) ifun.doit.insert((string)"wgraph_getDist");

      vs = findFunction(tmp, "getDistDense()");
      if(vs.size() == 3) ifun.doit.insert((string)"wgraph_getDistDense");

      vs = findFunction(tmp, "getDistForest()");
      if(vs.size() == 3) ifun.doit.insert((string)"wgraph_getDistForest");

      vs = findFunction(tmp, "BellmanFord()");
      if(vs.size() == 3) ifun.doit.insert((string)"wgraph_BellmanFord");

      vs = findFunction(tmp, "getDist01()");
      if(vs.size() == 3) ifun.doit.insert((string)"wgraph_getDist01");

      vs = findFunction(tmp, "MST_Prim_cost()");
      if(vs.size() == 3) ifun.doit.insert((string)"wgraph_MST_Prim_cost");

      vs = findFunction(tmp, "Rerooting()");
      if(vs.size() == 3) ifun.doit.insert((string)"wgraph_Rerooting");

      vs = findFunction(tmp, "numOfDaysInMonth()");
      if(vs.size() == 3) ifun.doit.insert((string)"numOfDaysInMonth1");

      vs = findFunction(tmp, "numOfDaysInMonth()");
      if(vs.size() == 3) ifun.doit.insert((string)"numOfDaysInMonth2");

      vs = findFunction(tmp, "getSum()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr1d_getSum");

      vs = findFunction(tmp, "setSum()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr1d_getSum");

      vs = findFunction(tmp, "setConstLenLeft()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr1d_ConstLenLeft");

      vs = findFunction(tmp, "ConstLenLeft()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr1d_ConstLenLeft");

      vs = findFunction(tmp, "ConstLenLeftCyclic()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr1d_ConstLenLeft");

      vs = findFunction(tmp, "setConstLenRight()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr1d_ConstLenRight");

      vs = findFunction(tmp, "ConstLenRight()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr1d_ConstLenRight");

      vs = findFunction(tmp, "ConstLenRightCyclic()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr1d_ConstLenRight");

      vs = findFunction(tmp, "setDHist()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr1d_dHist");

      vs = findFunction(tmp, "dHist()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr1d_dHist");

      vs = findFunction(tmp, "setSHist()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr1d_sHist");

      vs = findFunction(tmp, "sHist()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr1d_sHist");

      vs = findFunction(tmp, "setPrevLE()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr1d_PrevLE");

      vs = findFunction(tmp, "PrevLE()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr1d_PrevLE");

      vs = findFunction(tmp, "setPrevLT()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr1d_PrevLT");

      vs = findFunction(tmp, "PrevLT()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr1d_PrevLT");

      vs = findFunction(tmp, "setPrevGE()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr1d_PrevGE");

      vs = findFunction(tmp, "PrevGE()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr1d_PrevGE");

      vs = findFunction(tmp, "setPrevGT()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr1d_PrevGT");

      vs = findFunction(tmp, "PrevGT()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr1d_PrevGT");

      vs = findFunction(tmp, "setNextLE()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr1d_NextLE");

      vs = findFunction(tmp, "NextLE()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr1d_NextLE");

      vs = findFunction(tmp, "setNextLT()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr1d_NextLT");

      vs = findFunction(tmp, "NextLT()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr1d_NextLT");

      vs = findFunction(tmp, "setNextGE()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr1d_NextGE");

      vs = findFunction(tmp, "NextGE()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr1d_NextGE");

      vs = findFunction(tmp, "setNextGT()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr1d_NextGT");

      vs = findFunction(tmp, "NextGT()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr1d_NextGT");

      vs = findFunction(tmp, "getSum()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr2d_getSum");

      vs = findFunction(tmp, "getSumBorder()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr2d_getSum");

      vs = findFunction(tmp, "getSum45()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr2d_getSum45");

      vs = findFunction(tmp, "getSum45Border()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr2d_getSum45");

      vs = findFunction(tmp, "setConstLenLeft()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr2d_ConstLenLeft");

      vs = findFunction(tmp, "ConstLenLeft()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr2d_ConstLenLeft");

      vs = findFunction(tmp, "setConstLenRight()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr2d_ConstLenRight");

      vs = findFunction(tmp, "ConstLenRight()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr2d_ConstLenRight");

      vs = findFunction(tmp, "setConstLenUp()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr2d_ConstLenUp");

      vs = findFunction(tmp, "ConstLenUp()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr2d_ConstLenUp");

      vs = findFunction(tmp, "setConstLenDown()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr2d_ConstLenDown");

      vs = findFunction(tmp, "ConstLenDown()");
      if(vs.size() == 3) ifun.doit.insert((string)"Arr2d_ConstLenDown");

      vs = findFunction(tmp, "InnerProd()");
      if(vs.size() == 3) ifun.doit.insert((string)"InnerProd_array");

      vs = findFunction(tmp, "InnerProd()");
      if(vs.size() == 3) ifun.doit.insert((string)"Point2d_InnerProd");

      vs = findFunction(tmp, "CrossProd()");
      if(vs.size() == 3) ifun.doit.insert((string)"Point2d_CrossProd");

      vs = findFunction(tmp, "xysortA()");
      if(vs.size() == 3) ifun.doit.insert((string)"Point2d_xysortA");

      vs = findFunction(tmp, "argsortA()");
      if(vs.size() == 3) ifun.doit.insert((string)"Point2d_argsortA");

      vs = findFunction(tmp, "CCW()");
      if(vs.size() == 3) ifun.doit.insert((string)"Point2d_CCW");

      vs = findFunction(tmp, "ConvexHull_sorted()");
      if(vs.size() == 3) ifun.doit.insert((string)"Point2d_ConvexHull");

      vs = findFunction(tmp, "ConvexHull()");
      if(vs.size() == 3) ifun.doit.insert((string)"Point2d_ConvexHull");

      vs = findFunction(tmp, "PolygonArea2()");
      if(vs.size() == 3) ifun.doit.insert((string)"Point2d_PolygonArea2");

      vs = findFunction(tmp, "Determinant()");
      if(vs.size() == 3) ifun.doit.insert((string)"Determinant_Modint");

      vs = findFunction(tmp, "Determinant()");
      if(vs.size() == 3) ifun.doit.insert((string)"Determinant_modint");

      vs = findFunction(tmp, "Determinant()");
      if(vs.size() == 3) ifun.doit.insert((string)"Determinant_Mint");

      vs = findFunction(tmp, "Determinant()");
      if(vs.size() == 3) ifun.doit.insert((string)"Determinant_mint");


      vs = findFunction(tmp, "cntSubsequence<>()");
      if(vs.size() == 4) ifun.doit.insert((string)"cntSubsequence");

      vs = findFunction(tmp, "sumPrime<>()");
      if(vs.size() == 4) ifun.doit.insert((string)"sumPrime");




      vs = findFunction(tmp, "getDistT()");
      if(vs.size() == 3) ifun.doit.insert((string)"wgraph_getDistT");
      vs = findFunction(tmp, "getDistT<>()");
      if(vs.size() == 4) ifun.doit.insert((string)"wgraph_getDistT");
      
      break;
    }
    return tmp;
  }

  int isOperator(char c){
    if(c=='+' || c=='-' || c=='*' || c=='/' || c=='%'|| c=='=') return 1;
    if(c=='&' || c=='|' || c=='^' || c=='~') return 1;
    return 0;
  }


  /*
    vs = findFunction(tmp, "bsearch_min[]()"); とかで最初に見つかったもの
    rvs = vs;
    funstr = "bsearch_min[]()" or "not found"
    mode = "bsearch" or ""
  */
  void get_first_function(string tmp, vector<string> &rvs, string &funstr, string &modestr, string lim = ""){
    int len = 1000000000;
    vector<string> vs;

    funstr = "not found";
    modestr = "";

    if(lim == "" || lim == "bsearch"){
      vs = findFunction(tmp, "bsearch_min[]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "bsearch_min[]()", modestr = "bsearch";
      vs = findFunction(tmp, "bsearch_max[]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "bsearch_max[]()", modestr = "bsearch";
      vs = findFunction(tmp, "bsearch_min[][]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "bsearch_min[][]()", modestr = "bsearch";
      vs = findFunction(tmp, "bsearch_max[][]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "bsearch_max[][]()", modestr = "bsearch";
    }

    if(lim == "" || lim == "tsearch"){
      vs = findFunction(tmp, "tsearch_min[]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "tsearch_min[]()", modestr = "tsearch";
      vs = findFunction(tmp, "tsearch_argmin[]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "tsearch_argmin[]()", modestr = "tsearch";
      vs = findFunction(tmp, "tsearch_argminL[]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "tsearch_argminL[]()", modestr = "tsearch";
      vs = findFunction(tmp, "tsearch_max[]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "tsearch_max[]()", modestr = "tsearch";
      vs = findFunction(tmp, "tsearch_argmax[]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "tsearch_argmax[]()", modestr = "tsearch";
      vs = findFunction(tmp, "tsearch_argmaxL[]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "tsearch_argmaxL[]()", modestr = "tsearch";
      vs = findFunction(tmp, "tsearch_min[][]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "tsearch_min[][]()", modestr = "tsearch";
      vs = findFunction(tmp, "tsearch_argmin[][]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "tsearch_argmin[][]()", modestr = "tsearch";
      vs = findFunction(tmp, "tsearch_argminL[][]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "tsearch_argminL[][]()", modestr = "tsearch";
      vs = findFunction(tmp, "tsearch_max[][]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "tsearch_max[][]()", modestr = "tsearch";
      vs = findFunction(tmp, "tsearch_argmax[][]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "tsearch_argmax[][]()", modestr = "tsearch";
      vs = findFunction(tmp, "tsearch_argmaxL[][]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "tsearch_argmaxL[][]()", modestr = "tsearch";
    }

    if(lim == "" || lim == "minmax_function"){
      vs = findFunction(tmp, "argmin()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "argmin()", modestr = "minmax_function";
      vs = findFunction(tmp, "argmax()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "argmax()", modestr = "minmax_function";
      vs = findFunction(tmp, "argminL()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "argminL()", modestr = "minmax_function";
      vs = findFunction(tmp, "argmaxL()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "argmaxL()", modestr = "minmax_function";
      vs = findFunction(tmp, "sum[][]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "sum[][]()", modestr = "minmax_function";
      vs = findFunction(tmp, "sum[]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "sum[]()", modestr = "minmax_function";
      vs = findFunction(tmp, "mul[][]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "mul[][]()", modestr = "minmax_function";
      vs = findFunction(tmp, "mul[]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "mul[]()", modestr = "minmax_function";
      vs = findFunction(tmp, "gcd[][]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "gcd[][]()", modestr = "minmax_function";
      vs = findFunction(tmp, "gcd[]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "gcd[]()", modestr = "minmax_function";
      vs = findFunction(tmp, "lcm[][]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "lcm[][]()", modestr = "minmax_function";
      vs = findFunction(tmp, "lcm[]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "lcm[]()", modestr = "minmax_function";
      vs = findFunction(tmp, "min[][]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "min[][]()", modestr = "minmax_function";
      vs = findFunction(tmp, "min[]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "min[]()", modestr = "minmax_function";
      vs = findFunction(tmp, "max[][]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "max[][]()", modestr = "minmax_function";
      vs = findFunction(tmp, "max[]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "max[]()", modestr = "minmax_function";
      vs = findFunction(tmp, "argmin[][]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "argmin[][]()", modestr = "minmax_function";
      vs = findFunction(tmp, "argmin[]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "argmin[]()", modestr = "minmax_function";
      vs = findFunction(tmp, "argmax[][]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "argmax[][]()", modestr = "minmax_function";
      vs = findFunction(tmp, "argmax[]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "argmax[]()", modestr = "minmax_function";
      vs = findFunction(tmp, "argminL[][]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "argminL[][]()", modestr = "minmax_function";
      vs = findFunction(tmp, "argminL[]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "argminL[]()", modestr = "minmax_function";
      vs = findFunction(tmp, "argmaxL[][]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "argmaxL[][]()", modestr = "minmax_function";
      vs = findFunction(tmp, "argmaxL[]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "argmaxL[]()", modestr = "minmax_function";
      vs = findFunction(tmp, "XOR[][]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "XOR[][]()", modestr = "minmax_function";
      vs = findFunction(tmp, "XOR[]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "XOR[]()", modestr = "minmax_function";
      vs = findFunction(tmp, "OR[][]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "OR[][]()", modestr = "minmax_function";
      vs = findFunction(tmp, "OR[]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "OR[]()", modestr = "minmax_function";
      vs = findFunction(tmp, "AND[][]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "AND[][]()", modestr = "minmax_function";
      vs = findFunction(tmp, "AND[]()");
      if(vs.size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "AND[]()", modestr = "minmax_function";
    }

    if(lim == "" || lim == "gcdlcm"){
      if(g_flags.count("no-gcd()")==0){
        vs = findFunction(tmp, "gcd()");
        if(vs.size()) trim(vs[1]);
        if(vs.size() && vs[1].size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "gcd()", modestr = "gcdlcm";
      }
      if(g_flags.count("no-GCD()")==0){
        vs = findFunction(tmp, "GCD()");
        if(vs.size()) trim(vs[1]);
        if(vs.size() && vs[1].size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "GCD()", modestr = "gcdlcm";
      }
      if(g_flags.count("no-lcm()")==0){
        vs = findFunction(tmp, "lcm()");
        if(vs.size()) trim(vs[1]);
        if(vs.size() && vs[1].size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "lcm()", modestr = "gcdlcm";
      }
      if(g_flags.count("no-LCM()")==0){
        vs = findFunction(tmp, "LCM()");
        if(vs.size()) trim(vs[1]);
        if(vs.size() && vs[1].size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "LCM()", modestr = "gcdlcm";
      }
      if(g_flags.count("no-min()")==0){
        vs = findFunction(tmp, "min()");
        if(vs.size()) trim(vs[1]);
        if(vs.size()){
          string chk = vs[0];
          trim(chk);
          int len = chk.size();
          if(len >= 5 && chk.substr(len-5)=="std::") vs.clear();
        }
        if(vs.size() && vs[1].size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "min()", modestr = "gcdlcm";
      }
      if(g_flags.count("no-MIN()")==0){
        vs = findFunction(tmp, "MIN()");
        if(vs.size()) trim(vs[1]);
        if(vs.size() && vs[1].size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "MIN()", modestr = "gcdlcm";
      }
      if(g_flags.count("no-max()")==0){
        vs = findFunction(tmp, "max()");
        if(vs.size()){
          string chk = vs[0];
          trim(chk);
          int len = chk.size();
          if(len >= 5 && chk.substr(len-5,5)=="std::") vs.clear();
        }
        if(vs.size()) trim(vs[1]);
        if(vs.size() && vs[1].size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "max()", modestr = "gcdlcm";
      }
      if(g_flags.count("no-MAX()")==0){
        vs = findFunction(tmp, "MAX()");
        if(vs.size()) trim(vs[1]);
        if(vs.size() && vs[1].size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "MAX()", modestr = "gcdlcm";
      }
      if(g_flags.count("no-sum()")==0){
        vs = findFunction(tmp, "sum()");
        if(vs.size()) trim(vs[1]);
        if(vs.size() && vs[1].size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "sum()", modestr = "gcdlcm";
      }
      if(g_flags.count("no-SUM()")==0){
        vs = findFunction(tmp, "SUM()");
        if(vs.size()) trim(vs[1]);
        if(vs.size() && vs[1].size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "SUM()", modestr = "gcdlcm";
      }
      if(g_flags.count("no-mul()")==0){
        vs = findFunction(tmp, "mul()");
        if(vs.size()) trim(vs[1]);
        if(vs.size() && vs[1].size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "mul()", modestr = "gcdlcm";
      }
      if(g_flags.count("no-MUL()")==0){
        vs = findFunction(tmp, "MUL()");
        if(vs.size()) trim(vs[1]);
        if(vs.size() && vs[1].size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "MUL()", modestr = "gcdlcm";
      }
      if(g_flags.count("no-XOR()")==0){
        vs = findFunction(tmp, "XOR()");
        if(vs.size()) trim(vs[1]);
        if(vs.size() && vs[1].size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "XOR()", modestr = "gcdlcm";
      }
      if(g_flags.count("no-OR()")==0){
        vs = findFunction(tmp, "OR()");
        if(vs.size()) trim(vs[1]);
        if(vs.size() && vs[1].size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "OR()", modestr = "gcdlcm";
      }
      if(g_flags.count("no-AND()")==0){
        vs = findFunction(tmp, "AND()");
        if(vs.size()) trim(vs[1]);
        if(vs.size() && vs[1].size() && len > vs[0].size()) len = vs[0].size(), rvs = vs, funstr = "AND()", modestr = "gcdlcm";
      }
    }
  }

  
  
  void sentence_main(string tmpstr, string tmp, int tt, int &fg_return){
    pair<string, char> stchar;
    vector<string> vtmp;
    
    code_replace(tmp);
    
    //fprintf(stderr, "sentence main [tt %d] [%s] [%s]\n", tt, tmpstr.c_str(), tmp.c_str());
    
    stchar = nextToken(tmp);
    if(stchar.first == "return") fg_return = 1; else fg_return = 0;

    {
      string chk = "break";
      int iii, chknum = 2;
      while(chk.size() < stchar.first.size()){
        if(stchar.first == chk + "_break"){
          code *codec = this;
          iii = chknum;
          for(;;){
            if(codec->loop_type != "") iii--;
            if(iii==0) break;
            codec = codec->up;
          }
          tmp = "goto " + codec->break_label + ";";
          used_label.insert(codec->break_label);
          break;
        }
        if(stchar.first == chk + "_continue"){
          code *codec = this;
          iii = chknum;
          for(;;){
            if(codec->loop_type != "") iii--;
            if(iii==0) break;
            codec = codec->up;
          }
          tmp = "goto " + codec->continue_label + ";";
          used_label.insert(codec->continue_label);
          break;
        }
        chk += "_break";
        chknum++;
      }
    }

    tmp = sentence_twopointers(tmp); // TwoPointers()[][][];
    if(tmp=="") return;

    for(;;){
      string mode, funstr;
      vector<string> vs;

      get_first_function(tmp, vs, funstr, mode);

      if(mode == "") break;
      if(mode == "bsearch") tmp = sentence_bsearch_one(tmp);
      if(mode == "tsearch") tmp = sentence_tsearch_one(tmp);
      if(mode == "minmax_function") tmp = sentence_minmax_function_one(tmp);
      if(mode == "gcdlcm") tmp = sentence_gcdlcm_one(tmp);
      if(tmp=="") return;
    }

    /*
    tmp = sentence_bsearch(tmp); // bsearch_min[][](), bsearch_max[][]() など
    if(tmp=="") return;

    tmp = sentence_tsearch(tmp); // tsearch_min[][](), tsearch_max[][]() など
    if(tmp=="") return;

    tmp = sentence_minmax_function(tmp); // min[](), max[](), argmin[](), argmax[]()
    if(tmp=="") return;

    tmp = sentence_gcdlcm(tmp); // gcd(), lcm(), min(), max(), argmin(), argmax() など
    if(tmp=="") return;
    */

    tmp = sentence_if(tmp); // if[]
    if(tmp=="") return;

    tmp = sentence_dot_loop(tmp); // A[0..N-1] = 0;
    if(tmp=="") return;

    tmp = sentence_arraylike_operations(tmp); // (a,b) += (c,d); (a,b,c)++;など
    if(tmp=="") return;

    tmp = sentence_arraylike_sentence(tmp); // @[a,b].set(N);など
    if(tmp=="") return;

    tmp = sentence_otherfunctions(tmp); // runLength(), reduceFraction()
    if(tmp=="") return;

    tmp = sentence_inequation(tmp); // &&の省略
    if(tmp=="") return;

    tmp = sentence_hatena_minmax_operator(tmp); // >?=, <?=, **=
    if(tmp=="") return;

    tmp = sentence_times_operator(tmp); // *の省略
    if(tmp=="") return;

    tmp = sentence_pow_operator(tmp); // **
    if(tmp=="") return;

    tmp = sentence_div_operator(tmp); // /+
    if(tmp=="") return;

    tmp = sentence_sortE(tmp); // sortE
    if(tmp=="") return;

    tmp = sentence_reader(tmp); // rd()
    if(tmp=="") return;
    
    tmp = sentence_writer(tmp); // wt(), wtLn(), wtSp(), wtN()
    if(tmp=="") return;

    if(tmp.substr(0,10) == "inplace_L "){
      tmpstr = "inplace_L " + tmpstr;
      tmp = tmp.substr(10);
    }
    //fprintf(stderr, "sentence main last [%s] [%s]\n", tmpstr.c_str(), tmp.c_str());

    tmp = tmpstr + tmp;
    while(tmp.size() > 0){
      int i, j, k, fg = 0;
      char bef = ';';
      string stc;

      rep(i,tmp.size()-1){
        if(tmp[i] == '[' && tmp[i+1] == ']' && (bef==';' || bef==',' || bef=='=' || bef=='(' || bef==')' || bef=='{' || bef=='}')){
          REP(j,i,tmp.size()) if(tmp[j] == '{') break;
          if(j == tmp.size()) break;
          k = getBlockLength(tmp.substr(j));

          tt = -1;
          stc = tmp.substr(0,j);
          if(stc.substr(0,10) == "inplace_L ") stc = stc.substr(10);
          str.push_back(stc);
          nxt.push_back(-1);
          strtype.push_back((string)"sentence");

          stc = tmp.substr(j,k);
          insert(stc, str.size());

          tmp = tmp.substr(j+k);
          fg = 1;
          break;
        }
        if(!isspace(tmp[i])) bef = tmp[i];
      }
      if(fg) continue;
      
      str.push_back(tmp);
      nxt.push_back(-1);
      if(tt==-1){
        strtype.push_back((string)"sentence");
      } else {
        strtype.push_back((string)"sentence-var-def");
      }
      tmp = "";
    }
  }

  string equation_main(string tmp){
    string chk;
    tmp = sentence_times_operator(tmp);
    tmp = sentence_inequation(tmp);
    tmp = sentence_hatena_minmax_operator(tmp);
    tmp = sentence_pow_operator(tmp);
    tmp = sentence_div_operator(tmp);

    for(;;){
      string mode, funstr;
      vector<string> vs;

      get_first_function(tmp, vs, funstr, mode);

      if(mode == "") break;
      if(mode == "bsearch") tmp = sentence_bsearch_one(tmp);
      if(mode == "tsearch") tmp = sentence_tsearch_one(tmp);
      if(mode == "minmax_function") tmp = sentence_minmax_function_one(tmp);
      if(mode == "gcdlcm") tmp = sentence_gcdlcm_one(tmp);
      if(tmp=="") break;
    }
    return tmp;
  }

  // type_str: int
  // vardef_str @a--[n];
  // -> return : vec({"def", "int a[n];"}, {"exec", "rd((a--)(n));"})
  vector<pair<string,string>> set_vardef_doit(string type_str, string vardef_str){
    static map<pair<int,int>,string> dimstr;
    int i, j, k, sx, sy, adm, cdm, cfg, mm;
    string dstr; vector<string> dchk;
    int read_flag = 0;
    int regged = 0; // auto [a,b] = ...
    vector<string> type_vs, vardef_vs;
    vector<string> var_names;
    vector<pair<string,string>> res;

    //fprintf(stderr, "set_vardef_doit: [%s] [%s]\n", type_str.c_str(), vardef_str.c_str());

    trim(vardef_str);

    // [] no hokan st
    adm = 0; // [] no dim
    rep(sx, vardef_str.size()){
      if(vardef_str[sx] == '='){ adm = 0; break; }
      if(vardef_str[sx] == '['){
        sy = pairBracket(vardef_str, sx);
        sx = sy;
        adm++;
      }
    }
    dchk = findFunction(vardef_str, "operator[]");
    if(dchk.size()) adm = 0;

    if(adm){
      cdm = 0;
      if(adm){
        cdm = 0;
        rep(sx,vardef_str.size()){
          if(vardef_str[sx] == '='){ adm = 0; break; }
          if(vardef_str[sx] == '['){
            sy = pairBracket(vardef_str, sx);
            dstr = vardef_str.substr(sx+1, sy-sx-1);
            trim(dstr);
            if(dstr.size()){
              dimstr[make_pair(adm,cdm)] = dstr;
            } else if(dimstr.count(make_pair(adm,cdm))) {
              vardef_str = vardef_str.substr(0,sx+1) + dimstr[make_pair(adm,cdm)] + vardef_str.substr(sy);
            }
            sx = pairBracket(vardef_str, sx);
            cdm++;
          }
        }
      }
    }
    // [] no hokan ed

    sx = -1;
    if(type_str == "auto"){
      sx = 0;
      while(sx < vardef_str.size() && (isspace(vardef_str[sx]) || vardef_str[sx] == '&')) sx++;
    }
    if(type_str == "auto" && sx >= 0 && sx < vardef_str.size() && vardef_str[sx] == '['){
      sy = pairBracket(vardef_str, sx);
      if(sy >= 0){
        vector<string> vs;
        vs = split_p(vardef_str.substr(sx+1, sy-sx-1), ',');
        regged = 1;
        for(string s : vs){
          s = "inplace_L auto " + s + ";";
          pair<string,string> pss = var_definition(s);
          string str_var = pss.first;
          string str_type = pss.second;
          //fprintf(stderr,"auto -> [%s] [%s]\n", str_var.c_str(), str_type.c_str());
          add_localvar(str_var, str_type);
        }
      }
    }

    if(vardef_str[0] == '@'){
      read_flag = 1;
      vardef_str = vardef_str.substr(1);
    }

    cfg = 1;
    rep(sx, vardef_str.size()){
      if(read_flag == 0) break;
      if(vardef_str[sx] == '=') break;

      if(sx){
        if(isalnum(vardef_str[sx-1]) || vardef_str[sx-1]=='_') cfg = 0;
        else if(!isspace(vardef_str[sx-1])) cfg = 1;
      }

      if(cfg && vardef_str[sx] == '('){
        //fprintf(stderr, "cfg && vardef_str (\n");
        int ok = 1, tlen, i;
        vector<string> vs;
        sy = pairBracket(vardef_str, sx);
        vs = split_p(vardef_str.substr(sx+1, sy-sx-1), ',');
        for(string s : vs){
          string a, b;
          tlen = 0;
          rep(i,s.size()) if(isValidVarType(s.substr(0,i), s[i])) tlen = i;
          a = s.substr(0,tlen);
          b = s.substr(tlen);
          rep(i,b.size()){
            if(b[i]=='['){ i = pairBracket(b, i); continue; }
            if(isspace(b[i])) continue;
            if(isalnum(b[i]) || b[i]=='_') continue;
            if(b.substr(i,2)=="++" || b.substr(i,2)=="--"){i++; continue;}
            ok = 0;
          }
        }
        if(ok){
          for(string s : vs){
            string a, b;
            tlen = 0;
            rep(i,s.size()) if(isValidVarType(s.substr(0,i), s[i])) tlen = i;
            a = s.substr(0,tlen);
            b = s.substr(tlen);
            if(a.size() > 0) type_vs.push_back(a);
            else             type_vs.push_back(type_str);
            vardef_vs.push_back(vardef_str.substr(0,sx) + " " + b + vardef_str.substr(sy+1));
          }
          break;
        }
      }
    }
    if(vardef_vs.size()==0){
      type_vs.push_back(type_str);
      vardef_vs.push_back(vardef_str);
    }

    rep(mm, vardef_vs.size()){
      string reg, ad;

      //if(type_vs[i].substr(0,10) != "inplace_L") reg += "inplace_L ";
      reg += type_vs[mm];

      ad = vardef_vs[mm];
      //fprintf(stderr, "[1] [ad] %s\n", ad.c_str());
      if(read_flag){
        rep(i,ad.size()){
          if(ad[i]=='=') break;
          if(ad[i]=='['){
            i = pairBracket(ad, i);
            continue;
          }
          if(ad.substr(i,2) == "++" || ad.substr(i,2) == "--"){
            ad.erase(ad.begin()+i);
            ad.erase(ad.begin()+i);
            i--;
          }
        }
      }
      if(isalnum(ad[0]) || ad[0]=='_') reg += " ";
      reg += ad;

      if(reg.substr(0,10) != "inplace_L ") reg = "inplace_L " + reg;
      code_replace(reg);

      {
        int before_count = insert_count;
        pair<string,string> pss = var_definition(reg);
        string str_var = pss.first;
        string str_type = pss.second;

        if(insert_count != before_count){ // int i = min[k=0---4](a[k]); みたいなの
          vector<string> vvs = split_p(str_type, ',');
          reg = vvs[0] + " " + str_var + " = " + vvs[3] + ";";
          vvs[3] = "";
          str_type = vvs[0] + "," + vvs[1] + "," + vvs[2] + "," + vvs[3];
          if(str_type.substr(0,10) != "inplace_L ") str_type = "inplace_L " + str_type;
          res.push_back({"def", reg});
        } else {
          res.push_back({"def", reg});
        }
        
        int is_lambda = 0; // ラムダ式っぽいものはlocalvarに登録しない
        {
          int k;
          vector<string> vvs = split_p(str_type, ',');
          string rh = vvs[3];
          for(;;){
            trim(rh);
            if(rh[0] != '[') break;
            k = pairBracket(rh, 0);
            if(k <= 0) break;
            rh = rh.substr(k+1);

            trim(rh);
            if(rh[0] != '(') break;
            k = pairBracket(rh, 0);
            if(k <= 0) break;
            rh = rh.substr(k+1);

            // if(rh[0] != '{') break;

            is_lambda = 1;
            break;
          }
        }
        
        if(!is_lambda && !regged) add_localvar(str_var, str_type);
      }
    }

    if(read_flag){
      vector<string> dim_str;
      vardef_str = vardef_vs[0];
      //fprintf(stderr, "[2] [vardef_str] %s\n", vardef_str.c_str());
      for(string s : vardef_vs){
        string incdec;
        rep(i,s.size()){
          if(s[i] == '=') break;
          if(s[i] == '['){ i = pairBracket(s, i); continue; }
          if(s.substr(i,2) == "++" || s.substr(i,2) == "--"){
            incdec += s.substr(i,2);
            s.erase(s.begin()+i); s.erase(s.begin()+i);
            i--;
          }
        }
        while(s.size() && !(isalpha(s[0]) || s[0]=='_' || s.substr(0,2)=="++" || s.substr(0,2)=="--")) s = s.substr(1);
        rep(i,s.size()){
          if(isalnum(s[i])) continue;
          if(s[i]=='_') continue;
          if(s.substr(i,2)=="++"){ i++; continue; }
          if(s.substr(i,2)=="--"){ i++; continue; }
          s = s.substr(0,i);
          break;
        }
        var_names.push_back(s+incdec);
      }
      rep(sx,vardef_str.size()){
        if(vardef_str[sx] == '=') break;
        if(vardef_str[sx] == '[' || vardef_str[sx] == '('){
          sy = pairBracket(vardef_str, sx);
          dim_str.push_back( vardef_str.substr(sx+1, sy-sx-1) );
          sx = sy;
        }
      }

      if(dim_str.size()==0){
        string a = "exec", b;
        b = "rd(";
        rep(i,var_names.size()){
          if(i) b += ",";
          b += var_names[i];
        }
        b += ");";
        res.push_back({a,b});
      } else {
        string a = "exec", b;
        b = "rd((";
        rep(i,var_names.size()){
          if(i) b += ",";
          b += var_names[i];
        }
        b += ")(";
        rep(i,dim_str.size()){
          if(i) b += ",";
          b += dim_str[i];
        }
        b += "));";
        res.push_back({a,b});
      }
    }

    //fprintf(stderr,"return size %d\n", (int)res.size());
    //rep(i,res.size()) fprintf(stderr, "-- [%s] [%s]\n", res[i].first.c_str(), res[i].second.c_str());
    
    return res;
  }

  void set(string &in, string tp, string additional = ""){
    int i, j, k, st, cnt, tt, fg;
    int k1, k2, k3, k4, k5, ks, do_flag = 0;
    int end_type = 0, extra_sentence = 0;
    int fg_main = 0, fg_return = 0;
    string tmp, tmpstr, send_additional;
    string str_type, str_var;
    vector<string> str_arg;
    pair<string, char> stchar;
    vector<string> vtmp;

    type = tp;
    set_init();

    ftrim(in);
    if(tp != "program" && in[0]=='{') end_type = 1, in = in.substr(1);
    if(tp == "program") end_type = 2;
    if(additional != ""){
      in = additional + in;
      if(end_type == 0) extra_sentence = 1;
    }

    //fprintf(stderr, "name [%s] type [%s] loop_type [%s] continue_label[%s] break_label[%s]\n", name.c_str(), type.c_str(), loop_type.c_str(), continue_label.c_str(), break_label.c_str());

    for(;;){
      ftrim(in);
      tmp_vartype.clear();
      //fprintf(stderr, "--- %s\n", in.substr(0,10).c_str());

      if(in.substr(0,12) == "//no-insert-"){
        for(i=12;;i++) if(isspace(in[i]) || in[i]=='\0') break;
        g_flags.insert(in.substr(2,i-2));
        in = in.substr(i);
      }
      if(in.substr(0,13) == "//no-unlocked" && (isspace(in[13]) || in[13]=='\0')){
        g_flags.insert((string)"no-unlocked");
        in = in.substr(13);
        continue;
      }
      if(in.substr(0,10) == "//no-fread" && (isspace(in[10]) || in[10]=='\0')){
        g_flags.insert((string)"no-fread");
        in = in.substr(10);
        continue;
      }
      if(in.substr(0,11) == "//no-fwrite" && (isspace(in[11]) || in[11]=='\0')){
        g_flags.insert((string)"no-fwrite");
        in = in.substr(11);
        continue;
      }
      if(in.substr(0,13) == "//interactive" && (isspace(in[13]) || in[13]=='\0')){
        g_flags.insert((string)"no-unlocked");
        g_flags.insert((string)"interactive");
        in = in.substr(13);
        continue;
      }
      if(in.substr(0,10) == "//no-gcd()" && (isspace(in[10]) || in[10]=='\0')){
        g_flags.insert((string)"no-gcd()");
        in = in.substr(10);
        continue;
      }
      if(in.substr(0,10) == "//no-GCD()" && (isspace(in[10]) || in[10]=='\0')){
        g_flags.insert((string)"no-GCD()");
        in = in.substr(10);
        continue;
      }
      if(in.substr(0,10) == "//no-lcm()" && (isspace(in[10]) || in[10]=='\0')){
        g_flags.insert((string)"no-lcm()");
        in = in.substr(10);
        continue;
      }
      if(in.substr(0,10) == "//no-LCM()" && (isspace(in[10]) || in[10]=='\0')){
        g_flags.insert((string)"no-LCM()");
        in = in.substr(10);
        continue;
      }
      if(in.substr(0,10) == "//no-min()" && (isspace(in[10]) || in[10]=='\0')){
        g_flags.insert((string)"no-min()");
        in = in.substr(10);
        continue;
      }
      if(in.substr(0,10) == "//no-MIN()" && (isspace(in[10]) || in[10]=='\0')){
        g_flags.insert((string)"no-MIN()");
        in = in.substr(10);
        continue;
      }
      if(in.substr(0,10) == "//no-max()" && (isspace(in[10]) || in[10]=='\0')){
        g_flags.insert((string)"no-max()");
        in = in.substr(10);
        continue;
      }
      if(in.substr(0,10) == "//no-MAX()" && (isspace(in[10]) || in[10]=='\0')){
        g_flags.insert((string)"no-MAX()");
        in = in.substr(10);
        continue;
      }
      if(in.substr(0,10) == "//no-sum()" && (isspace(in[10]) || in[10]=='\0')){
        g_flags.insert((string)"no-sum()");
        in = in.substr(10);
        continue;
      }
      if(in.substr(0,10) == "//no-SUM()" && (isspace(in[10]) || in[10]=='\0')){
        g_flags.insert((string)"no-SUM()");
        in = in.substr(10);
        continue;
      }
      if(in.substr(0,10) == "//no-mul()" && (isspace(in[10]) || in[10]=='\0')){
        g_flags.insert((string)"no-mul()");
        in = in.substr(10);
        continue;
      }
      if(in.substr(0,10) == "//no-MUL()" && (isspace(in[10]) || in[10]=='\0')){
        g_flags.insert((string)"no-MUL()");
        in = in.substr(10);
        continue;
      }
      if(in.substr(0,25) == "//no-arraylike-operations" && (isspace(in[25]) || in[25]=='\0')){
        g_flags.insert((string)"no-arraylike-operations");
        in = in.substr(25);
        continue;
      }

      if(in.substr(0,2) == "//"){
        string tmp;
        rep(i,in.size()){
          if(in[i] == '\n'){ i++; break; }
          if(!isspace(in[i])) tmp += in[i];
        }
        in = in.substr(i);

        if(tmp.substr(0,17) == "//working_memory="){
          ll val;
          char buf[1000];
          tmp = tmp.substr(17);
          if(tmp[tmp.size()-1]=='B' || tmp[tmp.size()-1]=='b') tmp = tmp.substr(0,tmp.size()-1);
          val = atoi(tmp.c_str());
          if(tmp[tmp.size()-1]=='K' || tmp[tmp.size()-1]=='k') val *= 1024;
          if(tmp[tmp.size()-1]=='M' || tmp[tmp.size()-1]=='m') val *= 1024 * 1024;
          if(tmp[tmp.size()-1]=='G' || tmp[tmp.size()-1]=='g') val *= 1024 * 1024 * 1024;
          sprintf(buf, "inplace_L void *wmem; char memarr[%lld];", val);
          ifun.func["workmemory"] = buf;
        }
        continue;
      }
      
      if(in.substr(0,2) == "//"){
        rep(i,in.size()) if(in[i] == '\n'){ i++; break; }
        in = in.substr(i);
        continue;
      }
      if(in.substr(0,2) == "/*"){
        REP(i,1,in.size()) if(in.substr(i-1,2) == "*/"){ i++; break; }
        in = in.substr(i);
        continue;
      }

      if(end_type==-1) break;
      if(end_type==0){
        if(extra_sentence > 0) extra_sentence--;
        else                   end_type = -1;
      }

      if(in.size()==0 && end_type==2) break;
      if(in.size() && in[0]=='}' && end_type==1){
        in = in.substr(1);
        if(name == "int main()" && fg_return == 0){
          str.push_back((string)"return 0;");
          nxt.push_back(-1);
          strtype.push_back((string)"sentence");
        }
        break;
      }

      if(in.substr(0,8) == "#include"){
        cnt = 0;
        for(i=0;;i++){
          if(in[i]=='"') cnt++;
          if(cnt==2 || in[i]=='>'){ i++; break; }
        }
        tmp = in.substr(0, i);
        in = in.substr(i);
        str.push_back(tmp);
        nxt.push_back(-1);
        strtype.push_back((string)"sentence-include");
        continue;
      }

      if(in.substr(0,7) == "#define"){
        int fg_puts = 1;
        string in_tmp = in.substr(7);
        stchar = nextToken(in_tmp);
        if(user_code && stchar.first == "MD"){
          unsigned md, W, R, mdninv, RR, t;
          string sss, alls;
          char buf[100];
          for(i=0;;i++) if(in_tmp[i]=='\n' || in_tmp[i]=='\r') break;
          sss = in_tmp.substr(0, i);
          trim(sss);
          sss = sss.substr(2);
          trim(sss);
          sscanf(sss.c_str(), "%u", &md);
          sprintf(buf, "#define MD (%uU)\n", md);
          ifun.func["define_MD"] = (string)buf;

          W = 32;
          R = (1ULL << W) % md;
          RR = (ull)R*R % md;

          mdninv = 0;
          t = 0;
          rep(i,(int)W){
            if(t%2==0) t+=md, mdninv |= (1U<<i);
            t /= 2;
          }

          sprintf(buf, "#define MINT_W (%uU)\n", W);
          alls += buf;
          sprintf(buf, "#define MINT_R (%uU)\n", R);
          alls += buf;
          sprintf(buf, "#define MINT_RR (%uU)\n", RR);
          alls += buf;
          sprintf(buf, "#define MINT_MDNINV (%uU)\n", mdninv);
          alls += buf;
          ifun.func["define_for_Mint"] = alls;
          fg_puts = 0;

          sprintf(buf, "#define MD_PRIMITIVE_ROOT (%lldU)\n", primitiveRoot(md));
          ifun.func["define_MD_PRIMITIVE_ROOT"] = buf;
        }
        if(stchar.first == "PI") ifun.already.insert((string)"define_PI");
        for(i=0;;i++){
          if(in[i]=='\n' || in[i]=='\r') break;
        }
        tmp = in.substr(0, i);
        in = in.substr(i);
        if(fg_puts){
          str.push_back(tmp);
          nxt.push_back(-1);
          strtype.push_back((string)"sentence-define");
        }
        continue;
      }

      if(in.substr(0,6) == "#undef"){
        string in_tmp = in.substr(6);
        stchar = nextToken(in_tmp);
        for(i=0;;i++){
          if(in[i]=='\n' || in[i]=='\r') break;
        }
        tmp = in.substr(0, i);
        in = in.substr(i);
        str.push_back(tmp);
        nxt.push_back(-1);
        strtype.push_back((string)"sentence-undef");
        continue;
      }

      if(in.substr(0,7) == "#pragma"){
        string in_tmp = in.substr(7);
        stchar = nextToken(in_tmp);
        for(i=0;;i++){
          if(in[i]=='\n' || in[i]=='\r') break;
        }
        tmp = in.substr(0, i);
        in = in.substr(i);
        str.push_back(tmp);
        nxt.push_back(-1);
        strtype.push_back((string)"sentence-pragma");
        continue;
      }

      if(in.substr(0,6) == "#ifdef"){
        string in_tmp = in.substr(6);
        stchar = nextToken(in_tmp);
        for(i=0;;i++){
          if(in[i]=='\n' || in[i]=='\r') break;
        }
        tmp = in.substr(0, i);
        in = in.substr(i);
        str.push_back(tmp);
        nxt.push_back(-1);
        strtype.push_back((string)"sentence-ifdef");
        continue;
      }

      if(in.substr(0,6) == "#endif"){
        string in_tmp = in.substr(6);
        stchar = nextToken(in_tmp);
        for(i=0;;i++){
          if(in[i]=='\n' || in[i]=='\r') break;
        }
        tmp = in.substr(0, i);
        in = in.substr(i);
        str.push_back(tmp);
        nxt.push_back(-1);
        strtype.push_back((string)"sentence-endif");
        continue;
      }

      if(in.substr(0,5) == "#else"){
        string in_tmp = in.substr(5);
        stchar = nextToken(in_tmp);
        for(i=0;;i++){
          if(in[i]=='\n' || in[i]=='\r') break;
        }
        tmp = in.substr(0, i);
        in = in.substr(i);
        str.push_back(tmp);
        nxt.push_back(-1);
        strtype.push_back((string)"sentence-else");
        continue;
      }

      stchar = nextToken(in);
      if(stchar.second == ':'){
        for(i=0;;i++) if(in[i]==':'){ i++; break; }
        if(in[i]!=':'){
          tmp = in.substr(0, i);
          in = in.substr(i);
          str.push_back(tmp);
          nxt.push_back(-1);
          strtype.push_back((string)"sentence-label");
          continue;
        }
      }
      
      if(stchar.first == "case" || stchar.first == "default"){
        for(i=0;;i++) if(in[i]==':'){ i++; break; }
        tmp = in.substr(0, i);
        in = in.substr(i);
        str.push_back(tmp);
        nxt.push_back(-1);
        strtype.push_back((string)"sentence-case-label");
        continue;
      }
      
      tmpstr = "";
      for(;;){
        stchar = nextToken(in);
//        printf("[%s] [%c]\n",stchar.first.c_str(),stchar.second);
        if(stchar.first == "template"){
          cnt = 0;
          for(i=0;;i++){
            if(in[i]=='<') cnt++;
            if(in[i]=='>') cnt--;
            if(cnt==0 && in[i]=='>'){ i++; break; }
          }
          tmpstr += in.substr(0, i) + " ";
          in = in.substr(i);
          ftrim(in);
          continue;
        }

        if(stchar.first == "inline"){
          tmpstr += "inline ";
          in = in.substr(6);
          ftrim(in);
          continue;
        }
        
        break;
      }

      //fprintf(stderr,"1: %s\n",tmpstr.c_str());
      if(tmpstr.size()){
        string tt = tmpstr;
        vector<string> vtt;
        for(;;){
          int fg = -1;
          rep(i,tt.size()) if(tt.substr(i,8)=="template") fg = i;
          if(fg==-1) break;
          vtt.push_back(tt.substr(fg));
          tt = tt.substr(0,fg);
        }
        for(string tt : vtt){
          trim_until(tt, '<', '>');
          if(tt.size() >= 2){
            tt = tt.substr(1, tt.size()-2);
            vtmp = split_p(tt, ',');
            rep(i,vtmp.size()){
              trim(vtmp[i]);
              if(vtmp[i].substr(0,6)=="class "){
                vtmp[i] = vtmp[i].substr(5);
                trim(vtmp[i]);
                tmp_vartype.insert(vtmp[i]);
              }
              if(vtmp[i].substr(0,9)=="typename "){
                vtmp[i] = vtmp[i].substr(8);
                trim(vtmp[i]);
                tmp_vartype.insert(vtmp[i]);
              }
            }
          }
        }
      }

      fg = 0;
      tt = -1;
      {
        int lim = min(100, (int)in.size()+1);
        REP(i,1,lim) if(isValidVarType(in.substr(0,i), in[i])) tt = i;
      }

      stchar = nextToken(in);

      if(tt >= 0 || in.substr(0,8)=="operator"){
        for(i=0;;i++){
          int j;
          string tmp = in.substr(0,i);
          if(in[i] == '='){
            if(strpos_ns_t(tmp, "operator") == -1) break;
          }
          if(in[i] == ';') break;
          j = strpos_ns_t(tmp, "operator");
          if(j >= 0){
            j += 8;
            while(j < in.size() && isspace(in[j])) j++;
          }
          if(j == i && in[i] == '(' && in[i+1]==')'){
            i++; continue;
          }
          if(in[i] == '('){
            j = pairBracket(in, i);
            j++; while(j < in.size() && isspace(in[j])) j++;
            if(in[j]!='{' && in[j]!=':' && !(in[j]=='-' && in[j+1]=='>')) break;
            fg++; stchar.first = "function"; break;
          }
        }
      }
      //fprintf(stderr, "[function check] %d [%s] [%s]\n", tt, in.substr(0,20).c_str(), stchar.first.c_str());

      if(stchar.first == "if" && stchar.second =='(') fg = 1;
      if(stchar.first == "for" && stchar.second =='(') fg = 1;
      if(stchar.first == "rep" && (stchar.second =='(' || stchar.second == '[')) fg = 1;
      if(stchar.first == "REP" && (stchar.second =='(' || stchar.second == '[')) fg = 1;
      if(stchar.first == "rrep" && (stchar.second =='(' || stchar.second == '[')) fg = 1;
      if(stchar.first == "RREP" && (stchar.second =='(' || stchar.second == '[')) fg = 1;
      if(stchar.first == "rep_perm" && stchar.second == '(') fg = 1;
      if(stchar.first == "rep_mcomb" && stchar.second == '(') fg = 1;
      if(stchar.first == "rep_scomb" && stchar.second == '(') fg = 1;
      if(stchar.first == "rep_marr" && stchar.second == '(') fg = 1;
      if(stchar.first == "rep_sarr" && stchar.second == '(') fg = 1;
      if(stchar.first == "rep_dist" && stchar.second == '(') fg = 1;
      if(stchar.first == "while" && stchar.second =='(') fg = 1;
      if(stchar.first == "switch" && stchar.second =='(') fg = 1;
      if(stchar.first == "do") fg = 4;
      if(stchar.first == "struct") fg = 3;
      if(stchar.first == "class") fg = 3;
      if(stchar.first == "else") fg = 2;

      if(fg==2){
        string str_t = in.substr(4,100);
        pair<string,char> stch_t = nextToken(str_t);
        if(stch_t.first == "if" && stch_t.second == '('){
          fg = 1;
          stchar.first = "else if";
        }
      }

      if(fg || in[0]=='{'){
        fg_return = 0;

        if(fg==1 && stchar.first == "function") {
          rep(i,in.size()) if(in[i]=='{' || in[i] == ';') break;
        } else if(fg==1){
          cnt = 0;
          k4 = k5 = 0;
          rep(i,in.size()){
            string tmp = in.substr(0,i);
            if(tt>=0 && i>=tt) tmp=tmp.substr(tt);
            alltrim(tmp);
            while(tmp.size() && (tmp[0]=='*' || tmp[0]=='&')) tmp = tmp.substr(1);

            if(k4==0 && k5==0){
              if(in[i] == '(') cnt++;
              if(in[i] == ')') cnt--;
              if(in[i] == ')' && cnt == 0 && tmp != "operator("){
                i++; break;
              }
            }
            if(k5 == 0 && in[i]=='\'') k4 ^= 1;
            if(k4 == 0 && in[i]=='"') k5 ^= 1;
          }
        } else if(fg==2 || fg==4) {
          i = stchar.first.size();
        } else if(fg==3) {
          rep(i,in.size()) if(in[i]=='{' || in[i] == ';') break;
          tmp = in.substr(0, i);
          trim(tmp);
          vtmp = split_p(tmp, ' ');
          if(tmpstr.size() == 0) add_vartype(vtmp[vtmp.size()-1]);
          else                   add_tvartype(vtmp[vtmp.size()-1]);
        } else {
          i = 0;
        }
        tmp = in.substr(0, i);
        in = in.substr(i);

        //fprintf(stderr, "[[ %s /// %s /// fg = %d]]\n", tmp.c_str(), in.substr(0,10).c_str(), fg);
        //exit(0);

        code_replace(tmp);
        tmp = equation_main(tmp);
        
        tt = -1;
        REP(i,1,tmp.size()) if(isValidVarType(tmp.substr(0,i), tmp[i])) tt = i;

        if(type=="program" && tmp=="" && fg_main==0) tmp = "int main()", fg_main = 1;

        if(tmp.substr(0,8) == "rep_perm" || tmp.substr(0,9) == "rep_mcomb" || tmp.substr(0,9) == "rep_scomb" || tmp.substr(0,8) == "rep_marr" || tmp.substr(0,8) == "rep_sarr"){
          int len;
          string exec1, exec2, vn;
          string rep_mode;

          if(tmp.substr(0,8) == "rep_perm") rep_mode = "perm";
          if(tmp.substr(0,9) == "rep_mcomb") rep_mode = "mcomb";
          if(tmp.substr(0,9) == "rep_scomb") rep_mode = "scomb";
          if(tmp.substr(0,8) == "rep_marr") rep_mode = "marr";
          if(tmp.substr(0,8) == "rep_sarr") rep_mode = "sarr";
          
          while(tmp[0]!='(') tmp = tmp.substr(1);
          tmp = tmp.substr(1, tmp.size()-2);
          str_arg = split_p(tmp, ',');
          rep(i,str_arg.size()) trim(str_arg[i]);

          if(rep_mode == "perm") assert(rep_mode == "perm" && str_arg.size() == 2);
          if(rep_mode == "mcomb") assert(rep_mode == "mcomb" && str_arg.size() == 3);
          if(rep_mode == "scomb") assert(rep_mode == "scomb" && str_arg.size() == 3);
          if(rep_mode == "marr") assert(rep_mode == "marr" && str_arg.size() == 3);
          if(rep_mode == "sarr") assert(rep_mode == "sarr" && str_arg.size() == 3);

          if(rep_mode == "mcomb") ifun.doit.insert((string)"next_mcomb");
          if(rep_mode == "scomb") ifun.doit.insert((string)"next_scomb");
          if(rep_mode == "marr") ifun.doit.insert((string)"next_marr");
          if(rep_mode == "sarr") ifun.doit.insert((string)"next_sarr");

          if(rep_mode == "scomb" || rep_mode == "sarr") exec1 += "if(" + str_arg[2] + " >= " + str_arg[1] + "){";

          vn = getUnusedVarName();
          exec1 += "int " + vn + ";";
          if(rep_mode == "perm" || rep_mode == "scomb" || rep_mode == "sarr"){
            exec1 += "rep(" + vn + ", " + str_arg[1] + ") " + str_arg[0] + "[" + vn + "] = " + vn + ";";
          } else {
            exec1 += "rep(" + vn + ", " + str_arg[1] + ") " + str_arg[0] + "[" + vn + "] = 0;";
          }
          exec1 += "do{";

          exec2 += "}";
          if(rep_mode == "perm") exec2 += "while(next_permutation("+str_arg[0]+","+str_arg[0]+"+"+str_arg[1]+"));";
          if(rep_mode == "mcomb") exec2 += "while(next_mcomb("+str_arg[1]+","+str_arg[0]+","+str_arg[2]+"));";
          if(rep_mode == "scomb") exec2 += "while(next_scomb("+str_arg[1]+","+str_arg[0]+","+str_arg[2]+"));";
          if(rep_mode == "marr") exec2 += "while(next_marr("+str_arg[1]+","+str_arg[0]+","+str_arg[2]+"));";
          if(rep_mode == "sarr") exec2 += "while(next_sarr("+str_arg[1]+","+str_arg[0]+","+str_arg[2]+"));";
          
          if(rep_mode == "scomb" || rep_mode == "sarr") exec2 += "}";

          len = getBlockLength(in);
          in = "{" + exec1 + in.substr(0,len) + exec2 + "}" + in.substr(len);
          if(end_type == -1) end_type = 0;
          continue;
        }

        if(tmp.substr(0,8) == "rep_dist"){
          int j, k, len, mask, mx, mtmp, ddmx, zero_fg = 0;
          vector<string> rep_var, rep_center, rep_dvar;
          vector< vector<int> > vec_d; vector<int> chk_d;
          char cnt_buf[10], buf[10];
          string exec1, exec2, vn;
          string dist_str, dist_mode = "m"; int dist_mx = 1;
          
          while(tmp[0]!='(') tmp = tmp.substr(1);
          tmp = tmp.substr(1, tmp.size()-2);
          str_arg = split_p(tmp, ',');
          rep(i,str_arg.size()) trim(str_arg[i]);

          if(str_arg.size() % 2){
            dist_str = str_arg[str_arg.size()-1];
            str_arg.pop_back();
          }
          rep(k,str_arg.size()/2){
            rep_var.push_back(str_arg[k]);
            rep_center.push_back(str_arg[k + str_arg.size()/2]);
          }

          rep(k,rep_var.size()){
            if(localvar.count(rep_var[k])==0 && argvar.count(rep_var[k])==0 && globalvar.count(rep_var[k])==0){
              int fg = 1, idx;
              rep(idx,rep_var[k].size()) if(!isalnum(rep_var[k][idx]) && rep_var[k][idx]!='_') fg = 0;
              if(fg){
                string vardef;
                pair<string,string> vardefs;
                vardef = "int " + rep_var[k] + ";";
                vardefs = var_definition(vardef);
                add_localvar(vardefs.first, vardefs.second);
              }
            }
          }

          if(dist_str.size()){
            string x1, x2;
            rep(k,dist_str.size()){
              if(dist_str[k] == 'z'){
                zero_fg = 1;
                continue;
              }
              if(isdigit(dist_str[k])){
                x1 += dist_str[k];
                continue;
              }
              x2 += dist_str[k];
            }
            if(x1.size()) dist_mx = atoi(x1.c_str());
            if(x2.size()) dist_mode = x2;
          }
          if(dist_mode == "c") dist_mode = "t";
          if(dist_mode == "e"){
            dist_mx = dist_mx * dist_mx;
            dist_mode = "ee";
          }
          assert(dist_mode == "m" || dist_mode == "t" || dist_mode == "ee");

          ddmx = dist_mx;
          if(dist_mode == "ee"){
            ddmx = 0;
            while(ddmx * ddmx < dist_mx) ddmx++;
          }

          mx = 1;
          rep(k,rep_var.size()) mx *= 2 * ddmx + 1;
          rep(k,rep_var.size()) vec_d.push_back(vector<int>());
          rep(mask,mx){
            chk_d.resize(rep_var.size());
            mtmp = mask;
            rep(k,rep_var.size()){
              chk_d[rep_var.size()-1-k] = mtmp%(2*ddmx+1) - ddmx;
              mtmp /= (2*ddmx+1);
            }
            mtmp = 0;
            rep(k,rep_var.size()){
              if(dist_mode == "m") mtmp += abs(chk_d[k]);
              if(dist_mode == "t") mtmp = max(mtmp, abs(chk_d[k]));
              if(dist_mode == "ee") mtmp += chk_d[k] * chk_d[k];
            }
            if(mtmp == 0 && zero_fg == 0) continue;
            if(mtmp <= dist_mx){
              rep(k,rep_var.size()) vec_d[k].push_back(chk_d[k]);
            }
          }

          sprintf(cnt_buf, "%d", (int)vec_d[0].size());
          rep(k,rep_var.size()){
            rep_dvar.push_back(getUnusedVarName());
            exec1 += "static int " + rep_dvar[k] + "[" + cnt_buf + "] = {";
            rep(j,vec_d[k].size()){
              sprintf(buf, "%d", vec_d[k][j]);
              if(j) exec1 += ", ";
              exec1 += buf;
            }
            exec1 += "};\n";
          }
          
          vn = getUnusedVarName();
          exec1 += "int " + vn + ";";
          exec1 += "rep(" + vn + ", " + cnt_buf + "){\n";

          rep(k,rep_var.size()){
            exec1 += rep_var[k] + " = (" + rep_center[k] + ") + " + rep_dvar[k] + "[" + vn + "];\n";
          }

          exec2 += "}";

          len = getBlockLength(in);
          in = "{" + exec1 + in.substr(0,len) + exec2 + "}" + in.substr(len);
          if(end_type == -1) end_type = 0;
          continue;
        }

        if(tmp.substr(0,3) == "rep" || tmp.substr(0,3) == "REP" || tmp.substr(0,4) == "rrep" || tmp.substr(0,4) == "RREP"){
          int REP_fg = 0, rrep_fg = 0;
          string arr_name, arr_key, arr_str;

          if(tmp.substr(0,3) == "REP") REP_fg = 1;
          if(tmp.substr(0,4) == "rrep") rrep_fg = 1;
          if(tmp.substr(0,4) == "RREP") REP_fg = rrep_fg = 1;

          while(tmp[0]!='(' && tmp[0]!='[') tmp = tmp.substr(1);
          if(tmp[0]=='['){
            int t_ind = pairBracket(tmp, 0);
            vector<string> t_str;
            arr_str = tmp.substr(1, t_ind-1);
            tmp = tmp.substr(t_ind);
            trim(arr_str);
            trim(tmp);
            t_str = split_p(arr_str, ',');
            if(t_str.size() >= 1) arr_name = t_str[0];
            if(t_str.size() >= 2) arr_key = t_str[1];
          }
          while(tmp[0]!='(') tmp = tmp.substr(1);
          while(tmp[tmp.size()-1]!=')') tmp = tmp.substr(0, tmp.size()-1);
          tmp = tmp.substr(1, tmp.size()-2);
          str_arg = split_p(tmp, ',');
          rep(i,str_arg.size()) trim(str_arg[i]);

          if(str_arg.size()==1){
            str_arg.push_back(str_arg[0]);
            str_arg[0] = getUnusedVarName();
          }
          if(str_arg.size()==2){
            str_arg.push_back(str_arg[1]);
            str_arg[1] = "0";
          }

          if(arr_name != ""){
            if(arr_key == "") arr_key = getUnusedVarName();
            send_additional = "auto &" + str_arg[0] + " = " + arr_name + "[" + arr_key + "]" + ";";
            str_arg[0] = arr_key;
          }
          
          if(localvar.count(str_arg[0])==0 && argvar.count(str_arg[0])==0 && globalvar.count(str_arg[0])==0){
            int fg = 1, idx;
            rep(idx,str_arg[0].size()) if(!isalnum(str_arg[0][idx]) && str_arg[0][idx]!='_') fg = 0;
            if(fg){
              string vardef;
              pair<string,string> vardefs;
              vardef = "int " + str_arg[0] + ";";
              vardefs = var_definition(vardef);
              add_localvar(vardefs.first, vardefs.second);
            }
          }
          if(REP_fg){
            int agn;
            string doit, tvn;
            tvn = getUnusedVarName();
            
            if(rrep_fg==0) agn = 2;
            else           agn = 1;
            
            doit = "inplace_L int " + tvn + " = " + str_arg[agn] + ";";
            insert(doit, str.size());
            str_arg[agn] = tvn;
          }
          if(REP_fg && str_arg.size() >= 4){
            int agn = 3;
            string doit, tvn;
            tvn = getUnusedVarName();
            doit = "inplace_L int " + tvn + " = " + str_arg[agn] + ";";
            insert(doit, str.size());
            str_arg[agn] = tvn;
          }

          if(str_arg.size()==3 && rrep_fg==0){
            tmp = (string)"for(" + str_arg[0] + "=(" + str_arg[1] + ");" + str_arg[0] + "<(" + str_arg[2] + ");" + str_arg[0] + "++)";
          } else if(str_arg.size()==3 && rrep_fg==1){
            tmp = (string)"for(" + str_arg[0] + "=(" + str_arg[2] + ")-1;" + str_arg[0] + ">=(" + str_arg[1] + ");" + str_arg[0] + "--)";
          } else if(str_arg.size()==4 && rrep_fg==0){
            tmp = (string)"for(" + str_arg[0] + "=(" + str_arg[1] + ");" + str_arg[0] + "<(" + str_arg[2] + ");" + str_arg[0] + "+=(" + str_arg[3] + "))";
          } else if(str_arg.size()==4 && rrep_fg==1){
            tmp = (string)"for(" + str_arg[0] + "=(" + str_arg[2] + ")-1;" + str_arg[0] + ">=(" + str_arg[1] + ");" + str_arg[0] + "-=(" + str_arg[3] + "))";
          } else {
            assert(0);
          }

          //fprintf(stderr, "rep -> [%s]\n", tmp.c_str());
          
          stchar.first = "for";
        }

        code *cc = new code;
        cc->setUpnode(this);

        if(tmpstr.size()){
          string tt = tmpstr;
          vector<string> vtt;
          for(;;){
            int fg = -1;
            rep(i,tt.size()) if(tt.substr(i,8)=="template") fg = i;
            if(fg==-1) break;
            vtt.push_back(tt.substr(fg));
            tt = tt.substr(0,fg);
          }
          for(string tt : vtt){
            trim_until(tt, '<', '>');
            if(tt.size() >= 2){
              tt = tt.substr(1, tt.size()-2);
              vtmp = split_p(tt, ',');
              rep(i,vtmp.size()){
                trim(vtmp[i]);
                if(vtmp[i].substr(0,6)=="class "){
                  vtmp[i] = vtmp[i].substr(5);
                  trim(vtmp[i]);
                  cc->vartype.insert(vtmp[i]);
                }
                if(vtmp[i].substr(0,9)=="typename "){
                  vtmp[i] = vtmp[i].substr(8);
                  trim(vtmp[i]);
                  cc->vartype.insert(vtmp[i]);
                }
              }
            }
          }
        }
        if(stchar.first == "function"){
          int j;
          string ttt = tmp, tttt;
          string memostr, funname, restype;
          vector<string> var_def, var_type, var_name;
          pair<string,string> ss;

          //fprintf(stderr,"function [%d %d] %s\n",ttt.size(),tt,ttt.c_str());
          if(tt == -1){
            restype = "";
          } else {
            restype = isValidVarType_string(tmp.substr(0,tt), tmp[tt]);
          }

          //fprintf(stderr,"+-+-- %s\n",ttt.c_str());
          trim_until(ttt, '(', ')');
          assert(ttt.size()>=2);
          j = pairBracket(ttt, 0);
          ttt = ttt.substr(0,j+1);
          //fprintf(stderr,"+-+-- %s\n",ttt.c_str());

          tttt = ttt;
          alltrim(tttt);
          if(tttt.size() > 3 && tttt.substr(0,3)=="()("){
            int tind = pairBracket(ttt, 0);
            ttt = ttt.substr(tind+1);
          }
          
          ttt = ttt.substr(1, ttt.size()-2);

          rep(i,ttt.size()) if(ttt[i]==':'){
            REP(j,i+1,ttt.size()){
              if(!isspace(ttt[j])) break;
            }
            if(ttt.substr(j,7) == "Memoize"){
              memostr = ttt.substr(i);
              rep(i,ttt.size()) if(ttt.substr(i,memostr.size()) == memostr){
                ttt = ttt.substr(0,i) + ttt.substr(i+memostr.size());
                break;
              }
              rep(i,tmp.size()) if(tmp.substr(i,memostr.size()) == memostr){
                tmp = tmp.substr(0,i) + tmp.substr(i+memostr.size());
                break;
              }
              while(memostr.size() && memostr[0] != 'M') memostr = memostr.substr(1);
              break;
            }
          }
          rep(j,tmp.size()) if(tmp[j]=='(') break;
          j--;
          while(isspace(tmp[j])) j--;
          i = j;
          while(i >= 0 && (isalnum(tmp[i]) || tmp[i]=='_')) i--;
          funname = tmp.substr(i+1, j-i);

          vtmp = split_p3(ttt, ',');
          rep(i,vtmp.size()){
            trim(vtmp[i]);
            if(vtmp[i].size()==0 || vtmp[i]=="void") continue;
            ss = cc->var_definition(vtmp[i]);
            cc->argvar[ss.first] = ss.second;
            rep(j,ss.second.size()) if(ss.second[j]==',') break;
            var_type.push_back(ss.second.substr(0,j));
            var_name.push_back(ss.first);
            var_def.push_back(vtmp[i]);
          }

          if(memostr.size()){
            vector<int> range1, range2;
            string varnames, varnames_memo, vardefs;
            string str_add;

            rep(i,(int)memostr.size()) if(memostr[i]=='[') break;
            if(i < (int)memostr.size()){
              vector<string> memovs, memovs2;
              trim_until(memostr, '[', ']');
              memostr = memostr.substr(1, memostr.size()-2);
              memovs = split_p3(memostr, ',');
              rep(i,(int)memovs.size()){
                memovs2 = split_p3(memovs[i], ':');
                rep(j,(int)memovs2.size()) trim(memovs2[j]);
                if(memovs2.size()==1){
                  range1.push_back(0);
                  range2.push_back(atoi(memovs2[0].c_str())-1);
                } else {
                  range1.push_back(atoi(memovs2[0].c_str()));
                  range2.push_back(atoi(memovs2[1].c_str()));
                }
              }
            }

            assert(range1.size() == 0 || range1.size() == var_type.size());

            rep(i,(int)var_type.size()){
              if(i) varnames += ",";
              varnames += var_name[i];
              if(i) vardefs += ",";
              vardefs += var_def[i];

            }
            if(range1.size()==0){
              varnames_memo = varnames;
              if(var_type.size() != 1) varnames_memo = "{" + varnames_memo + "}";
            } else {
              rep(i,(int)var_type.size()){
                varnames_memo += "[" + var_name[i];
                if(range1[i]) varnames_memo += " - (" + to_string(range1[i]) + ")";
                varnames_memo += "]";
              }
            }

            rep(i,tmp.size()) if(tmp[i]=='(') break;
            i--;
            while(isspace(tmp[i])) i--;
            tmp = tmpstr + tmp.substr(0,i+1) + "_cL_func" + tmp.substr(i+1);

            if(range1.size()){
              str_add += restype + " " + funname + "_cL_memo_val";
              rep(i,range1.size()) str_add += "[" + to_string(range2[i] - range1[i] + 1) + "]";
              str_add += ";\n";
              str_add += "char " + funname + "_cL_memo_exist";
              rep(i,range1.size()) str_add += "[" + to_string(range2[i] - range1[i] + 1) + "]";
              str_add += ";\n";
            } else {
              str_add += "map<";
              if(var_type.size()==1){
                str_add += var_type[0];
              } else {
                str_add += "tuple<";
                rep(i,(int)var_type.size()){
                  if(i) str_add += ",";
                  str_add += var_type[i];
                }
                str_add += ">";
              }
              str_add += "," + restype + "> " + funname + "_cL_memo_map;\n";
            }

            str_add += tmp + ";\n";

            str_add += restype + " " + funname + "(" + vardefs + "){\n";
            if(range1.size()){
              str_add += "if(!" + funname + "_cL_memo_exist" + varnames_memo + "){\n";
              str_add += funname + "_cL_memo_val" + varnames_memo + " = " + funname + "_cL_func(" + varnames + ");\n";
              str_add += funname + "_cL_memo_exist" + varnames_memo + " = 1;\n";
              str_add += "}\n";
              str_add += "return " + funname + "_cL_memo_val" + varnames_memo + ";\n";
            } else {
              str_add += "auto it = " + funname + "_cL_memo_map.find(" + varnames_memo + ");\n";
              str_add += "if(it != " + funname + "_cL_memo_map.end()) return it->second;\n";
              str_add += "return " + funname + "_cL_memo_map[" + varnames_memo + "] = " + funname + "_cL_func(" + varnames + ");\n";
            }
            str_add += "}\n";

            str_add += "void " + funname + "_clear(){\n";
            if(range1.size()){
              vector<string> vnnn;
              rep(i,range1.size()) vnnn.push_back(getUnusedVarName());
              str_add += "int ";
              rep(i,range1.size()){
                if(i) str_add += ",";
                str_add += vnnn[i];
              }
              str_add += ";\n";
              rep(i,range1.size()){
                str_add += "rep(" + vnnn[i] + "," + to_string(range2[i] - range1[i] + 1) + ")";
              }
              str_add += funname + "_cL_memo_exist";
              rep(i,range1.size()){
                str_add += "[" + vnnn[i] + "]";
              }
              str_add += " = 0;\n";
            } else {
              str_add += funname + "_cL_memo_map.clear();\n";
            }
            str_add += "}\n";

            in = str_add + tmp + in;
            continue;
            //fprintf(stderr, "---\n%s\n---\n",str_add.c_str());
          }

          //fprintf(stderr, "function\n");
          //fprintf(stderr, "ttt [%s]\n", ttt.c_str());
          //fprintf(stderr, "tmp [%s]\n", tmp.c_str());
          //fprintf(stderr, "memostr [%s]\n", memostr.c_str());
          //fprintf(stderr, "funname [%s]\n", funname.c_str());
          //fprintf(stderr, "restype [%s]\n", restype.c_str());
          //rep(i,var_type.size()) fprintf(stderr, "var_type [%s] var_name [%s] var_def[%s]\n", var_type[i].c_str(), var_name[i].c_str(), var_def[i].c_str());

        }
        if(stchar.first == "for"){
        }
        
        if(stchar.first == "for"){
          //fprintf(stderr,"for-1 [%s]\n",tmp.c_str());
          string a, b; vector<string> c;
          int sx, sy = -1;
          rep(sx,tmp.size()) if(tmp[sx]=='(') break;
          if(sx < tmp.size()) sy = pairBracket(tmp, sx);
          if(sy >= 0){
            a = tmp.substr(sx+1, sy-sx-1);
            //fprintf(stderr,"for-2 [%s]\n",a.c_str());
            rep(sx,a.size()) if(a[sx]==';' || a[sx]==':' || a[sx]=='"'){ a = a.substr(0,sx); break; }
            a += ";";
            //fprintf(stderr,"for-3 [%s]\n",a.c_str());
            for(;;){
              trim(a);
              if(a.substr(0,7)=="static "){a = a.substr(7); continue;}
              if(a.substr(0,6)=="const "){a = a.substr(6); continue;}
              if(a.substr(0,6)=="final "){a = a.substr(6); continue;}
              break;
            }
            sy = -1;
            rep(sx,a.size()) if(isValidVarType(a.substr(0,sx), a[sx])) sy = sx;
            if(sy >= 0){
              b = a.substr(0,sy);
              c = split_p(a.substr(sy), ',');
              //fprintf(stderr,"for-var [%s] [%s]\n", a.substr(0,sy).c_str(), a.substr(sy).c_str());
              for(string cc : c) set_vardef_doit(b, cc);
              for(auto it = localvar.begin(); it != localvar.end(); it++){
                if(cc->localvar.count(it->first)==0) cc->localvar[it->first] = it->second;
              }
            }
          }
        }

        tmp = sentence_otherfunctions(tmp);

        if(stchar.first == "for" || (do_flag==0 && stchar.first == "while") || stchar.first == "do"){
          cc->loop_type = stchar.first;
          cc->continue_label = getUnusedVarName();
          cc->break_label = getUnusedVarName();
        }
        if(stchar.first == "do") do_flag = 1;
        if(stchar.first == "while") do_flag = 0;
        cc->name = tmpstr + tmp;
        cc->set(in, (string)"block", send_additional); send_additional = "";
        nxt.push_back(nxtlst.size());
        nxtlst.push_back(cc);
        str.push_back(tmpstr + tmp);
        strtype.push_back((string)"block-"+stchar.first);
        if(cc->loop_type != ""){
          str.push_back(cc->break_label+":;");
          nxt.push_back(-1);
          strtype.push_back((string)"auto-label");
        }
        continue;
      }

      k1 = k2 = k3 = k4 = k5 = ks = 0;
      st = 0;
      if(tt != -1) st = tt;
      
      for(i=st;;i++){
        if(i==in.size()){
          fprintf(stderr, "something wrong? : around %s\n", in.substr(st,100).c_str());
          exit(1);
        }
        
        if(k4==0 && k5==0){
          if(in[i] == '(') k1++;
          if(in[i] == ')') k1--;
          if(in[i] == '[') k2++;
          if(in[i] == ']') k2--;
          if(in[i] == '{') k3++;
          if(in[i] == '}') k3--;
        }
        if(k5==0 && in[i] == '\'') k4^=1;
        if(k4==0 && in[i] == '"') k5^=1;
        if( (k4||k5) && in[i] == '\\') {i++; continue;}

        if(k1==0 && k2==0 && k3==0 && k4==0 && k5==0){
          string insubstr = in.substr(0,i);
          int stpos = 0, hatena_flag = 0, using_typedef_flag = 0;
          while(stpos < i){
            stpos = strpos_ns(insubstr, "?", stpos);
            if(stpos == -1) break;
            if(insubstr.substr(stpos-1,3) == ">?=" || insubstr.substr(stpos-1,3) == "<?="){
              stpos++;
              continue;
            }
            hatena_flag = 1;
            break;
          }
          if(in.substr(st,5) == "using" && isspace(in[st+5])){
            int pp;
            string chkstr;
            vector<string> chk;
            REP(pp,st,i){
              if(isspace(in[pp])){
                if(chkstr.size()) chk.push_back(chkstr);
                chkstr = "";
              } else {
                chkstr += in[pp];
              }
            }
            if(chk.size() >= 3 && chk[2] == "="){
              using_typedef_flag = 1;
              if(tmpstr.substr(0,8) == "template"){
                tvartype.insert(chk[1]);
              } else {
                vartype.insert(chk[1]);
              }
            }
          }
          if(i==0 || (!isalnum(in[i-1]) && in[i-1]!='_')){
            rep(k,50){
              if(i+k > in.size()) break;
              if(isValidVarType(in.substr(i,k), in[i+k])) ks = max(ks, i+k);
            }
          }
          if((hatena_flag == 0 && using_typedef_flag == 0 && ks <= i && in[i] == ',') || in[i] == ';'){
            static map<pair<int,int>,string> dimstr;
            string additional = "", advarname;
            int sx, sy, adm, cdm; string dstr; vector<string> dchk;
            tmp = in.substr(st, i+1-st);
            /*if(tt>=0)*/ tmp[tmp.size()-1] = ';';

            if(tt >= 0){
              int exi;
              string type_str, vardef_str;
              vector<pair<string,string>> ex;

              type_str = in.substr(0, tt);
              vardef_str = tmp;
              ex = set_vardef_doit(type_str, vardef_str);

              rep(exi, ex.size()){
                string type, sen;
                type = ex[exi].first;
                sen = ex[exi].second;
                //fprintf(stderr,"-send [%s] [%s] [%d]\n", tmpstr.c_str(), tmp.c_str(), tt);
                if(type=="exec") sentence_main("", sen, -1, fg_return);
                if(type=="def") sentence_main(tmpstr, sen, tt, fg_return);
              }
            } else {
              sentence_main(tmpstr, tmp, tt, fg_return);
            }

            if(in[i] == ','){
              st = i+1;
              continue;
            } else {
              i++;
              break;
            }
          }
        }
      }
      in = in.substr(i);
      continue;
    }
    if(loop_type != ""){
      str.push_back(continue_label+":;");
      nxt.push_back(-1);
      strtype.push_back((string)"auto-label");
    }
  }

  void debug_writer(int tab = 0){
    int i, j;
    
    std::set<string>::iterator it1;
    map<string,string>::iterator it2;
    printf("vartype  :");
    for(it1=vartype.begin();it1!=vartype.end();it1++) printf(" (%s)", it1->c_str()); puts("");
    printf("tvartype :");
    for(it1=tvartype.begin();it1!=tvartype.end();it1++) printf(" (%s)", it1->c_str()); puts("");
    printf("localvar :");
    for(it2=localvar.begin();it2!=localvar.end();it2++) printf(" (%s->%s)", it2->first.c_str(), it2->second.c_str()); puts("");
    printf("globalvar:");
    for(it2=globalvar.begin();it2!=globalvar.end();it2++) printf(" (%s->%s)", it2->first.c_str(), it2->second.c_str()); puts("");
    printf("argvar   :");
    for(it2=argvar.begin();it2!=argvar.end();it2++) printf(" (%s->%s)", it2->first.c_str(), it2->second.c_str()); puts("");
    
    rep(i,str.size()){
      rep(j,tab) printf(" ");
      printf("%4d: %s   [%s]\n",i,str[i].c_str(),strtype[i].c_str());
      if(nxt[i]>=0) nxtlst[nxt[i]]->debug_writer(tab+4);
    }
  }

  void output_vardef(int mode, int tab){
    int i, f;
    std::set<string> tp;
    std::set<string>::iterator it;
    vector<string> tmp;
    map<string,string>::iterator mt;

    for(mt=localvar.begin(); mt!=localvar.end(); mt++){
      //fprintf(stderr, "%s\n", mt->second.c_str());
      tmp = split_p2(mt->second, ',');
      if(tmp.size()!=4){
        fprintf(stderr,"tmp: size = %d\n", (int)tmp.size());
        rep(i,tmp.size()) fprintf(stderr,"tmp[%d] = %s\n",i,tmp[i].c_str());
      }
      assert(tmp.size()==4);
      if(tmp[0].substr(0,10) == "inplace_L ") continue;
      tp.insert(tmp[0]);
    }
    for(it=tp.begin(); it!=tp.end(); it++){
      rep(i,tab) printf(" ");
      f = 0;
      for(mt=localvar.begin(); mt!=localvar.end(); mt++){
        tmp = split_p2(mt->second, ',');
        assert(tmp.size()==4);
        if(*it == tmp[0]){
          if(f) printf(", ");
          else  printf("%s ", it->c_str());
          f++;
          printf("%s%s%s", tmp[1].c_str(), mt->first.c_str(), tmp[2].c_str());
          if(tmp[3].size() != 0) printf("=%s",tmp[3].c_str());
        }
      }
      printf(";\n");
    }
  }

  void output(int mode, int tab = 0){
    int i, j;
    int var = 0;
    string s;
    
    rep(i,str.size()+1){

      if(i==str.size()){
        if(var==0) output_vardef(mode, tab);
        var = 1;
        break;
      }
      
      if(strtype[i]=="sentence-var-def"){
        if(str[i].substr(0,10) != "inplace_L ") continue;
        str[i] = str[i].substr(10);
      }

      if(strtype[i]=="auto-label"){
        s = str[i];
        s = s.substr(0, s.size()-2);
        if(used_label.count(s)==0) continue;
      }

      if(var==0 && (type != "program" || strtype[i] != "block-inserted" || i==str.size()-1)){
        var = 1;
        output_vardef(mode, tab);
      }
      
      if(strtype[i]=="block-inserted"){
        nxtlst[nxt[i]]->output(mode, tab);
        continue;
      }
      
      rep(j,tab) printf(" ");
      s = str[i];
      if(g_flags.count((string)"no-unlocked")) replaceAll_ns_t(s, "my_putchar_unlocked", "my_putchar");
      if(g_flags.count((string)"no-unlocked")) replaceAll_ns_t(s, "my_getchar_unlocked", "my_getchar");
      if(g_flags.count((string)"no-unlocked")) replaceAll_ns_t(s, "putchar_unlocked", "putchar");
      if(g_flags.count((string)"no-unlocked")) replaceAll_ns_t(s, "getchar_unlocked", "getchar");
      if(g_flags.count((string)"no-unlocked")) replaceAll_ns_t(s, "fread_unlocked", "fread");
      if(g_flags.count((string)"no-unlocked")) replaceAll_ns_t(s, "fwrite_unlocked", "fwrite");
//      if(g_flags.count((string)"no-fread"))    replaceAll_ns_t(s, "my_getchar_unlocked", "getchar_unlocked");
//      if(g_flags.count((string)"no-fread"))    replaceAll_ns_t(s, "my_getchar", "getchar");
//      if(g_flags.count((string)"no-fwrite"))   replaceAll_ns_t(s, "my_putchar_unlocked", "putchar_unlocked");
//      if(g_flags.count((string)"no-fwrite"))   replaceAll_ns_t(s, "my_putchar", "putchar");
      printf("%s",s.c_str());

      if(strtype[i]=="block-while"){
        if(nxtlst[nxt[i]]->is_empty_block()){
          printf(";\n");
          continue;
        }
      }
      
      if(strtype[i].substr(0,5)=="block"){
        printf("{\n");
      } else {
        printf("\n");
      }
      if(nxt[i]>=0) nxtlst[nxt[i]]->output(mode, tab+(mode==0?2:0));
      if(strtype[i].substr(0,5)=="block"){
        rep(j,tab) printf(" ");
        printf("}\n");
      }
    }
  }


};


int main(int argc, char **argv){
  int i, k, f;
  char buf[10001];
  string str, str_store;

  if(argc >= 2) str = argv[1];
  if(str == "-lib"){

    str = "";
    for(;;){
      k = fread(buf, 1, 10000, stdin);
      buf[k] = '\0';
      str += buf;
      if(k < 10000) break;
    }

    putchar('"');
    rep(i,(int)str.size()){
      if(str.substr(i,6)=="-----\n"){
        putchar('"');
        putchar('\n');
        putchar('"');
        i += 5;
        continue;
      }
      if(str.substr(i,7)=="-----\r\n"){
        putchar('"');
        putchar('\n');
        putchar('"');
        i += 6;
        continue;
      }
      if(str[i]=='"'){
        putchar('\\');
        putchar('"');
      } else if(str[i]=='\\'){
        putchar('\\');
        putchar('\\');
      } else if(str[i]=='\n'){
        putchar('\\');
        putchar('n');
      } else if(str[i]=='\r'){
      } else {
        putchar(str[i]);
      }
    }
    putchar('"');
    
    return 0;
  }

  if(str == "-libshow"){
    int i;
    
    ifun.set();

    rep(i,ifun.name.size()){
      puts("-------------------------------");
      printf("%s\n------\n", ifun.name[i].c_str());
      printf("%s\n\n\n", ifun.func[ifun.name[i]].c_str());
    }

    return 0;
  }

  for(;;){
    k = fread(buf, 1, 10000, stdin);
    buf[k] = '\0';
    str += buf;
    if(k < 10000) break;
  }
  str_store = str;

  { // delete comments
    vector<string> functional_comments;
    functional_comments.push_back("//no-insert-");
    functional_comments.push_back("//no-unlocked");
    functional_comments.push_back("//no-fread");
    functional_comments.push_back("//no-fwrite");
    functional_comments.push_back("//interactive");
    functional_comments.push_back("//no-gcd()");
    functional_comments.push_back("//no-GCD()");
    functional_comments.push_back("//no-lcm()");
    functional_comments.push_back("//no-LCM()");
    functional_comments.push_back("//no-min()");
    functional_comments.push_back("//no-MIN()");
    functional_comments.push_back("//no-max()");
    functional_comments.push_back("//no-MAX()");
    functional_comments.push_back("//no-sum()");
    functional_comments.push_back("//no-SUM()");
    functional_comments.push_back("//no-mul()");
    functional_comments.push_back("//no-MUL()");
    functional_comments.push_back("//no-arraylike-operations");
    functional_comments.push_back("//working_memory");

    int i = 0, j, k1 = 0, k2 = 0, c1 = 0, c2 = 0;
    str = "";
    while(i < str_store.size()){
      if(c1){
        if(str_store.substr(i,2) == "*/"){
          c1 ^= 1;
          i += 2;
          continue;
        }
        i++;
        continue;
      }
      if(c2){
        if(str_store[i] == '\n') c2 ^= 1;
        else i++;
        continue;
      }
      if(str_store[i] == '\\'){
        str += str_store[i++];
        str += str_store[i++];
        continue;
      }
      if(k2 == 0 && str_store[i] == '\''){
        k1 ^= 1;
        str += str_store[i++];
        continue;
      }
      if(k1 == 0 && str_store[i] == '"'){
        k2 ^= 1;
        str += str_store[i++];
        continue;
      }
      if(k1 == 0 && k2 == 0){
        if(str_store.substr(i,2) == "/*"){
          c1 ^= 1;
          i += 2;
          continue;
        }
        if(str_store.substr(i,2) == "//"){
          rep(j,functional_comments.size()){
            if(str_store.substr(i,functional_comments[j].size()) == functional_comments[j]) break;
          }
          if(j == functional_comments.size()){
            c2 ^= 1;
            i += 2;
            continue;
          }
        }
      }
      str += str_store[i++];
    }
  }
  str += "\n";

  // find main
  {
    int i = 0, j, k;
    int bf = ';', bfns = ';';
    int k1 = 0, k2 = 0;
    int found = 0;

    while(i < str.size() && !found){
      if(str[i] == '\\'){
        i += 2;
        continue;
      }
      if(k2 == 0 && str[i] == '\''){
        k1 ^= 1;
        i++;
        continue;
      }
      if(k1 == 0 && str[i] == '"'){
        k2 ^= 1;
        i++;
        continue;
      }
      if(k1 || k2){
        i++;
        continue;
      }

      if(str[i] == '#'){
        bf = bfns = ';';
        while(i < str.size() && str[i] != '\n') i++;
        continue;
      }

      if(bf == ';' || bf == '}' || isspace(bf)){
        if((bfns == ';' || bfns == '}') && str[i] == '{'){
          found = 1;
        }
        if(str.substr(i,4) == "main"){
          j = i + 4;
          while(j < (int)str.size() && isspace(str[j])) j++;
          if(j < (int)str.size() && str[j] == '(') found = 1;
        }
      }

      if(!isspace(str[i])) bfns = str[i];
      bf = str[i];
      i++;
    }
    if(!found){
      str = "{" + str + "}\n";
    }
  }

  ifun.set();
  string tmp;

  code c;
  c.up = NULL;
  c.set(str, (string)"program");
  user_code = 0;

  if(g_flags.count((string)"no-fread")){
    ifun.func["my_getchar_unlocked"] = "inline int my_getchar_unlocked(){ return getchar_unlocked(); }";
  }
  if(g_flags.count((string)"no-fwrite")){
    ifun.func["my_putchar_unlocked"] = "inline void my_putchar_unlocked(const int k){ putchar_unlocked(k); }";
  }
  if(g_flags.count((string)"interactive")){
    ifun.func["my_getchar_unlocked"] = "inline int my_getchar_unlocked(){ return getchar_unlocked(); }";
    ifun.func["my_putchar_unlocked"] = "inline void my_putchar_unlocked(const int k){ putchar_unlocked(k); if(k=='\\n') fflush(stdout); }";
  }

  tmp = ifun.get_insert_string((string)"first");
  c.insert(tmp, 0);
  
  tmp = ifun.get_insert_string((string)"main_first");
  rep(i,c.str.size()){
    if(c.str[i].substr(0,9) == "int main(" || c.str[i].substr(0,10) == "void main("){
      c.nxtlst[c.nxt[i]]->insert(tmp,0);
    }
  }
  
  tmp = ifun.get_insert_string((string)"last");
  c.insert(tmp, c.str.size());
  
//  c.debug_writer(); return 0;
  c.output(0);
  printf("// cLay version 20250224-1\n");


  str = str_store;
  f = 0;
  printf("\n// --- original code ---\n");
  rep(i,(int)str.size()){
    if(!f) printf("// "), f = 1;
    if(str[i]=='\r') continue;
    putchar(str[i]);
    if(str[i]=='\n') f = 0;
  }

  return 0;
}



/*
メモ：

【機能追加っぽいの】

【バグ修正っぽいの】

【その他っぽいの】

【ドキュメントに追記したもの（元々使えるもの）】

【ドキュメントの修正】

【追加したい】

det

多倍長??
LList?

graph 枝数
graph 連結成分
graph isTree
graph isForest
graph 重心分解


*/


