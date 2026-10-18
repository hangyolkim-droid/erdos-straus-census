#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static uint64_t primes[42000]; static int np=0;
static void sieve(uint64_t lim){char*c=calloc(lim+1,1);for(uint64_t i=2;i<=lim;i++)if(!c[i]){primes[np++]=i;for(uint64_t j=i*i;j<=lim;j+=i)c[j]=1;}free(c);}
static inline uint64_t gcd64(uint64_t a,uint64_t b){while(b){uint64_t t=a%b;a=b;b=t;}return a;}
static uint64_t mulmod(uint64_t a,uint64_t b,uint64_t m){return (__uint128_t)a*b%m;}
static uint64_t powmod(uint64_t a,uint64_t e,uint64_t m){uint64_t r=1;a%=m;while(e){if(e&1)r=mulmod(r,a,m);a=mulmod(a,a,m);e>>=1;}return r;}
static int isp(uint64_t n){
  if(n<2)return 0;
  for(int i=0;i<12&&(uint64_t)primes[i]*(uint64_t)primes[i]<=n;i++) if(n%primes[i]==0) return n==primes[i];
  uint64_t d=n-1;int s=0;while(!(d&1)){d>>=1;s++;}
  uint64_t b[]={2,3,5,7,11,13,17,19,23,29,31,37};
  for(int i=0;i<12;i++){
    uint64_t x=powmod(b[i],d,n); if(x==1||x==n-1)continue;
    int ok=0; for(int r=1;r<s;r++){x=mulmod(x,x,n); if(x==n-1){ok=1;break;}}
    if(!ok) return 0;
  }
  return 1;
}
static long ndiv=0;
static int works(uint64_t n,uint32_t delta){
  uint64_t a=(n+delta)>>2;
  uint64_t na=n*a;
  uint64_t g=gcd64(delta,na); uint32_t p=delta/g;
  if(p==1) return 1;
  uint64_t q=na/g;
  unsigned char reach[128]; memset(reach,0,p); reach[1%p]=1;
  // a/g is fully factored; n is prime
  uint64_t x=a/g;
  for(int i=0;i<np && primes[i]*primes[i]<=x;i++){
    if(x%primes[i]==0){
      int e=0; while(x%primes[i]==0){x/=primes[i];e++; ndiv++;}
      uint64_t pk=1; unsigned char nr[128];
      memset(nr,0,p);
      for(int k=0;k<=2*e;k++){
        if(pk<p) for(uint32_t r=0;r<p;r++) if(reach[r]) nr[(r*pk)%p]=1;
        pk=(pk*(primes[i]%p))%p;
      }
      memcpy(reach,nr,p);
    }
    ndiv++;
  }
  if(x>1){ // x prime, exponent 2
    uint64_t b=x%p;
    unsigned char nr[128]; memset(nr,0,p);
    for(uint32_t r=0;r<p;r++) if(reach[r]){ nr[(r*b)%p]=1; nr[(r*b%p*b)%p]=1; }
    memcpy(reach,nr,p);
  }
  // multiply by n^j, j=0,1,2
  uint64_t nm=n%p, n2=nm*nm%p;
  unsigned char fin[128]; memset(fin,0,p);
  for(uint32_t r=0;r<p;r++) if(reach[r]){ fin[r]=1; fin[(r*nm)%p]=1; fin[(r*n2)%p]=1; }
  return fin[(p-(q%p))%p];
}
int main(int argc,char**argv){
  int N = argc>1?atoi(argv[1]):20000;
  sieve(500000);
  struct timespec t0,t1; clock_gettime(CLOCK_MONOTONIC,&t0);
  uint64_t found=0, tested=0; int maxd=0; uint64_t maxn=0; long tot=0;
  uint64_t seed=88172645463325252ULL;
  while(tested<(uint64_t)N){
    seed^=seed<<13; seed^=seed>>7; seed^=seed<<17;
    uint64_t n=200000000000ULL + (seed%800000000000ULL);
    n|=1; n=(n/4)*4+1;
    if(!isp(n)) continue;
    tested++;
    uint32_t d=0;
    for(uint32_t delta=3;delta<=107;delta+=4){ tot++; if(works(n,delta)){d=delta;break;} }
    if(!d) found++;
    if((int)d>maxd){maxd=d;maxn=n;}
  }
  clock_gettime(CLOCK_MONOTONIC,&t1);
  double s=(t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)/1e9;
  printf("primes tested: %llu in %.3f s\n",(unsigned long long)tested,s);
  printf("per prime: %.2f us | delta-tests/prime: %.2f | trial divisions/prime: %.1f\n",1e6*s/tested,(double)tot/tested,(double)ndiv/tested);
  printf("max delta %d at n=%llu ; needs >107: %llu\n",maxd,(unsigned long long)maxn,(unsigned long long)found);
  printf("extrapolate: 1.4e7 primes (decade 1e10-1e11) = %.0f s single-core; 10x for next decade\n", 1.4e7*s/tested);
  return 0;
}
