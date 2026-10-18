#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#define BLK 30000000ULL
#define SEG 10000000ULL
#define CAP 20000U

static uint8_t *b3c;
static uint32_t *spfN, *spfA;
static uint64_t SEGL,SEGH,ALo,AHi;
static uint8_t *R1,*R2;
static uint64_t base[16384]; static int nb=0;
static uint64_t N,NA,KMAX;

static inline void bset(uint8_t*b,uint64_t i){ b[i>>3] |= (uint8_t)(1u<<(i&7)); }
static inline int  bget(const uint8_t*b,uint64_t i){ return (b[i>>3]>>(i&7))&1; }
static inline uint64_t gcd64(uint64_t a,uint64_t b){while(b){uint64_t t=a%b;a=b;b=t;}return a;}

static void sieve_base(uint64_t lim){
  char *c=calloc(lim+1,1);
  for(uint64_t i=2;i<=lim;i++) if(!c[i]){ base[nb++]=i; for(uint64_t j=i*i;j<=lim;j+=i) c[j]=1; }
  free(c);
}
static void build_spf(uint32_t*spf,uint64_t lo,uint64_t hi){
  memset(spf,0,(hi-lo)*sizeof(uint32_t));
  for(int i=0;i<nb;i++){
    uint64_t p=base[i]; if(p*p>=hi) break;
    uint64_t s=(lo+p-1)/p*p; if(s<p*p) s=p*p;
    for(uint64_t x=s;x<hi;x+=p) if(spf[x-lo]==0) spf[x-lo]=(uint32_t)p;
  }
}
static int bad3_of(uint64_t x,uint64_t lo,uint64_t hi,const uint32_t*spf){
  while(x>1){
    uint64_t pp;
    if(x>=lo&&x<hi){ uint32_t s=spf[x-lo]; pp=s?(uint64_t)s:x; }
    else { pp=0; for(int i=0;i<nb;i++){ if(base[i]*base[i]>x){pp=x;break;} if(x%base[i]==0){pp=base[i];break;} } if(!pp)pp=x; }
    if(pp==3||pp%3==2) return 1;
    while(x%pp==0) x/=pp;
  }
  return 0;
}
static int factor_into(uint64_t x,uint64_t lo,uint64_t hi,const uint32_t*spf,uint64_t*pr,int*ex,int k){
  while(x>1){
    uint64_t pp;
    if(x>=lo && x<hi){ uint32_t sp=spf[x-lo]; pp = sp ? (uint64_t)sp : x; }
    else { pp=0; for(int i=0;i<nb;i++){ if(base[i]*base[i]>x){pp=x;break;} if(x%base[i]==0){pp=base[i];break;} } if(!pp)pp=x; }
    int e=0; while(x%pp==0){x/=pp;e++;}
    int f=0; for(int i=0;i<k;i++) if(pr[i]==pp){ex[i]+=e;f=1;break;}
    if(!f){ pr[k]=pp; ex[k]=e; k++; }
  }
  return k;
}
static int works(uint64_t n,uint32_t delta){
  uint64_t a=(n+delta)/4, na=n*a;
  uint64_t g=gcd64(delta,na); uint32_t p=delta/g;
  if(p==1) return 1;
  uint64_t q=na/g;
  uint64_t pr[64]; int ex[64]; int k=0;
  k=factor_into(n,SEGL,SEGH,spfN,pr,ex,k);
  k=factor_into(a,ALo,AHi,spfA,pr,ex,k);
  uint64_t x=g;
  for(int i=0;i<nb && base[i]*base[i]<=x;i++){
    uint64_t pp=base[i];
    if(x%pp==0){ int e=0; while(x%pp==0){x/=pp;e++;} for(int j=0;j<k;j++) if(pr[j]==pp){ex[j]-=e;break;} }
  }
  if(x>1) for(int j=0;j<k;j++) if(pr[j]==x){ex[j]-=1;break;}
  memset(R1,0,p); R1[1%p]=1;
  for(int i=0;i<k;i++){
    if(ex[i]<=0) continue;
    uint64_t pw[160]; int nw=0, lim=2*ex[i]; if(lim>159)lim=159;
    uint64_t cur=1%p;
    for(int e=0;e<=lim;e++){ pw[nw++]=cur; cur=(cur*(pr[i]%p))%p; }
    memset(R2,0,p);
    for(uint32_t r=0;r<p;r++) if(R1[r]) for(int w=0;w<nw;w++) R2[(r*pw[w])%p]=1;
    memcpy(R1,R2,p);
  }
  return R1[(p-(q%p))%p];
}

int main(int argc,char**argv){
  N = argc>1?strtoull(argv[1],0,10):10000000000ULL;
  uint32_t cap = argc>2?(uint32_t)strtoull(argv[2],0,10):CAP;
  clock_t t0=clock();
  NA=N/4+8; KMAX=N/4;
  fprintf(stderr,"N=%llu NA=%llu cap=%u\n",(unsigned long long)N,(unsigned long long)NA,cap);

  sieve_base(100000);
  fprintf(stderr,"base primes=%d\n",nb);

  b3c=calloc(KMAX/8+1,1);
  if(!b3c){fprintf(stderr,"alloc failed\n");return 1;}

  uint8_t *blk=malloc(BLK);
  for(uint64_t lo=2;lo<N;lo+=BLK){
    uint64_t hi=lo+BLK; if(hi>N) hi=N;
    memset(blk,1,hi-lo);
    for(int i=0;i<nb;i++){
      uint64_t p=base[i]; if(p*p>=hi) break;
      uint64_t s=(lo+p-1)/p*p; if(s<p*p) s=p*p;
      for(uint64_t x=s;x<hi;x+=p) blk[x-lo]=0;
    }
    for(uint64_t x=lo;x<hi;x++){
      if(!blk[x-lo]) continue;
      uint64_t p=x;
      if(p==2) continue;
      if(p==3 || p%3==2){
        uint64_t m0 = (p%4==1)?1:3;
        uint64_t n0 = p*m0;
        if(n0<4*KMAX){ for(uint64_t k=(n0-1)/4;k<KMAX;k+=p) bset(b3c,k); }
      }
    }
  }
  free(blk);
  fprintf(stderr,"marking done %.1fs\n",(double)(clock()-t0)/CLOCKS_PER_SEC);

  spfN=malloc(SEG*4); spfA=malloc((SEG/4+cap+8)*4);
  R1=malloc(cap+2); R2=malloc(cap+2);
  if(!spfN||!spfA||!R1||!R2){fprintf(stderr,"alloc2 failed\n");return 1;}

  uint64_t surv=0,unres=0,segs=0; uint32_t gmax=0; uint64_t gmaxn=0;
  uint32_t dmax[12]={0}; uint64_t dmaxn[12]={0}, cnt[12]={0};

  for(uint64_t L=2;L<N;L+=SEG){
    uint64_t R=L+SEG; if(R>N) R=N;
    SEGL=L; SEGH=R; ALo=L/4; AHi=(R+cap)/4+4;
    build_spf(spfN,L,R); build_spf(spfA,ALo,AHi);
    uint64_t n0=L; if(n0%4!=1) n0 += (1+4-(n0%4))%4;
    for(uint64_t n=n0;n<R;n+=4){
      if(n<5) continue;
      if(bget(b3c,(n-1)/4)) continue;
      uint64_t a3=(n+3)/4;
      if(bad3_of(a3,ALo,AHi,spfA)) continue;
      surv++;
      uint32_t d=0;
      for(uint32_t delta=7;delta<=cap;delta+=4) if(works(n,delta)){d=delta;break;}
      if(!d){unres++;continue;}
      int dec=0; uint64_t t=n; while(t>=10){t/=10;dec++;}
      if(dec>11)dec=11;
      cnt[dec]++;
      if(d>dmax[dec]){dmax[dec]=d;dmaxn[dec]=n;}
      if(d>gmax){gmax=d;gmaxn=n;}
    }
    segs++;
    if(segs%100==0) fprintf(stderr,"  seg %llu surv=%llu gmax=%u %.1fs\n",(unsigned long long)segs,(unsigned long long)surv,gmax,(double)(clock()-t0)/CLOCKS_PER_SEC);
  }
  printf("N=%llu survivors=%llu unresolved=%llu\n",(unsigned long long)N,(unsigned long long)surv,(unsigned long long)unres);
  printf("GLOBAL MAX delta=%u at n=%llu\n",gmax,(unsigned long long)gmaxn);
  for(int i=1;i<=11;i++) if(cnt[i]) printf("decade 1e%d : max delta=%u at n=%llu (count %llu)\n",i,dmax[i],(unsigned long long)dmaxn[i],(unsigned long long)cnt[i]);
  fprintf(stderr,"total %.1fs\n",(double)(clock()-t0)/CLOCKS_PER_SEC);
  return 0;
}
