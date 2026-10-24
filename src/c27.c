#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <time.h>

#define KS (1u<<20)
static uint64_t primes[80000]; static int np=0;
static uint64_t inv9240[80000];
static uint64_t LO,HI;
static uint32_t RES[]={1,169,289,361,529,841,961,1369,1681,1849,2041,2209,2521,2641,2689,2809,
                       3361,3481,3529,3721,4321,4489,5041,5161,5329,5569,6169,6241,6889,
                       7561,7681,7921,8089,8761,1201,6001};
static int NRES=36;
static volatile int next_res=0;
static pthread_mutex_t mtx=PTHREAD_MUTEX_INITIALIZER;

static inline uint64_t gcd64(uint64_t a,uint64_t b){while(b){uint64_t t=a%b;a=b;b=t;}return a;}
static uint64_t mulmod(uint64_t a,uint64_t b,uint64_t m){return (__uint128_t)a*b%m;}
static uint64_t powmod(uint64_t a,uint64_t e,uint64_t m){uint64_t r=1;a%=m;while(e){if(e&1)r=mulmod(r,a,m);a=mulmod(a,a,m);e>>=1;}return r;}

static void sieve_primes(uint64_t lim){
  char*c=calloc(lim+1,1);
  for(uint64_t i=2;i<=lim;i++) if(!c[i]){
    if(np>=80000){ fprintf(stderr,"primes[] overflow at %llu (cap 80000)\n",(unsigned long long)i); exit(1); }
    primes[np++]=i; for(uint64_t j=i*i;j<=lim;j+=i) c[j]=1;
  }
  free(c);
  for(int i=0;i<np;i++) if(primes[i]>11) inv9240[i]=powmod(9240%primes[i],primes[i]-2,primes[i]);
}

static int works(uint64_t n,uint32_t delta){
  uint64_t a=(n+delta)>>2, na=n*a;
  uint64_t g=gcd64(delta,na); uint32_t p=delta/g;
  if(p==1) return 1;
  uint64_t q=na/g;
  uint8_t reach[128], nr[128];
  memset(reach,0,p); reach[1%p]=1;
  uint64_t x=a/g;
  for(int i=0;i<np && primes[i]*primes[i]<=x;i++){
    if(x%primes[i]==0){
      int e=0; while(x%primes[i]==0){x/=primes[i];e++;}
      uint64_t pw=1; memset(nr,0,p);
      for(int k=0;k<=2*e;k++){
        for(uint32_t r=0;r<p;r++) if(reach[r]) nr[(r*pw)%p]=1;
        pw=(pw*(primes[i]%p))%p;
      }
      memcpy(reach,nr,p);
    }
  }
  if(x>1){
    uint64_t b=x%p; memset(nr,0,p);
    for(uint32_t r=0;r<p;r++) if(reach[r]){ nr[r]=1; nr[(r*b)%p]=1; nr[(r*b%p*b)%p]=1; }
    memcpy(reach,nr,p);
  }
  uint64_t nm=n%p, n2=nm*nm%p; memset(nr,0,p);
  for(uint32_t r=0;r<p;r++) if(reach[r]){ nr[r]=1; nr[(r*nm)%p]=1; nr[(r*n2)%p]=1; }
  return nr[(p-(q%p))%p];
}

typedef struct { int idx; uint64_t cnt; uint32_t maxd; uint64_t maxn; uint64_t over; } Res;
static Res res[64];

static void* worker(void*arg){
  uint8_t *buf=malloc(KS);
  for(;;){
    pthread_mutex_lock(&mtx); int idx=next_res++; pthread_mutex_unlock(&mtx);
    if(idx>=NRES) break;
    uint64_t r=RES[idx];
    uint64_t kA0=(LO>r)?(LO-r+9239)/9240:0;
    uint64_t kB=(HI-1-r)/9240;
    uint64_t cnt=0; uint32_t maxd=0; uint64_t maxn=0, over=0;
    for(uint64_t kA=kA0;kA<=kB;kA+=KS){
      uint64_t len=kB-kA+1; if(len>KS) len=KS;
      memset(buf,1,len);
      for(int i=0;i<np;i++){
        uint64_t p=primes[i];
        if(p<=11) continue;
        if(p*p>HI) break;
        uint64_t k0=(p-(r*inv9240[i])%p)%p;
        uint64_t j=(k0+p-(kA%p))%p;
        for(;j<len;j+=p) buf[j]=0;
      }
      for(uint64_t j=0;j<len;j++){
        if(!buf[j]) continue;
        uint64_t n=r+9240*(kA+j);
        cnt++;
        uint32_t d=0;
        for(uint32_t delta=3;delta<=107;delta+=4) if(works(n,delta)){d=delta;break;}
        if(!d){
          over++;
          for(uint32_t delta=111;delta<=2000;delta+=4) if(works(n,delta)){d=delta;break;}
          if(d>maxd){maxd=d;maxn=n;}
        } else if(d>maxd){maxd=d;maxn=n;}
      }
    }
    res[idx].idx=idx; res[idx].cnt=cnt; res[idx].maxd=maxd; res[idx].maxn=maxn; res[idx].over=over;
  }
  return NULL;
}

int main(int argc,char**argv){
  LO=argc>1?strtoull(argv[1],0,10):1000000000ULL;
  HI=argc>2?strtoull(argv[2],0,10):10000000000ULL;
  struct timespec t0,t1; clock_gettime(CLOCK_MONOTONIC,&t0);
  sieve_primes(1000000);
  /* self-test */
  if(!works(8803369,107)){ printf("SELFTEST FAIL: 8803369 not solved at 107\n"); return 1; }
  for(uint32_t d=3;d<107;d+=4) if(works(8803369,d)){ printf("SELFTEST FAIL: 8803369 solved at %u\n",d); return 1; }
  printf("selftest ok (8803369 -> delta 107, none below)\n"); fflush(stdout);
  pthread_t th[4];
  for(int i=0;i<4;i++) pthread_create(&th[i],NULL,worker,NULL);
  for(int i=0;i<4;i++) pthread_join(th[i],NULL);
  uint64_t cnt=0,over=0; uint32_t maxd=0; uint64_t maxn=0;
  for(int i=0;i<NRES;i++){ cnt+=res[i].cnt; over+=res[i].over;
    if(res[i].maxd>maxd){maxd=res[i].maxd;maxn=res[i].maxn;} }
  clock_gettime(CLOCK_MONOTONIC,&t1);
  double s=(t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)/1e9;
  printf("range [%llu,%llu)  primes tested=%llu  in %.1f s\n",(unsigned long long)LO,(unsigned long long)HI,(unsigned long long)cnt,s);
  printf("MAX delta=%u at n=%llu\n",maxd,(unsigned long long)maxn);
  printf("primes needing delta>107: %llu\n",(unsigned long long)over);
  return 0;
}
