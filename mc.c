/* mc.c -- Metropolis dynamics of the Ising array of a sign matrix H of order N.
   Energy E = -(1/sqrt N) sum_jk H_jk a_j b_k = -Omega/sqrt N, and we report m = Omega/N^(3/2) = -E/N.
   We keep the integer local fields A = Hb and B = H^t a, so that a flip costs O(N).
   Usage: mc matrixfile seed sweeps T_1 T_2 ... T_k [-start startfile]
   For each temperature in the given order (cooling or heating schedule), runs the given number of sweeps
   (1 sweep = 2N attempted flips), and prints T, the mean of m over the second half of the sweeps, its variance
   times N (for the specific heat), and the largest m seen.  At the end, quenches at T = 0 to a local minimum
   and prints the final m.  The start is random, unless a file with the 2N signs of (a,b) is given. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
static unsigned long long s[2];
static unsigned long long nxt(void){unsigned long long s1=s[0],s0=s[1];s[0]=s0;s1^=s1<<23;s[1]=s1^s0^(s1>>17)^(s0>>26);return s[1]+s0;}
static double unif(void){return (nxt()>>11)*(1.0/9007199254740992.0);}
int main(int argc,char**argv){
  if(argc<5){fprintf(stderr,"usage: mc matrixfile seed sweeps T1 T2 ... [-start file]\n");return 1;}
  FILE*f=fopen(argv[1],"r"); int N=0,x,cap=1<<20,n=0; int*buf=malloc(cap*sizeof(int));
  while(fscanf(f,"%d",&x)==1){ if(n==cap){cap*=2; buf=realloc(buf,cap*sizeof(int));} buf[n++]=x; } fclose(f);
  N=(int)llround(sqrt((double)n)); if(N*N!=n){fprintf(stderr,"not square\n");return 1;}
  int*H=buf; s[0]=0x9E3779B97F4A7C15ULL*(atoll(argv[2])+1); s[1]=0xBF58476D1CE4E5B9ULL^atoll(argv[2]); for(int i=0;i<20;i++)nxt();
  long sweeps=atol(argv[3]); int nT=0; double T[1000]; char*start=NULL;
  for(int i=4;i<argc;i++){ if(!strcmp(argv[i],"-start")){start=argv[++i];} else T[nT++]=atof(argv[i]); }
  int*a=malloc(N*sizeof(int)),*b=malloc(N*sizeof(int)); long*A=calloc(N,sizeof(long)),*B=calloc(N,sizeof(long));
  if(start){FILE*g=fopen(start,"r"); for(int j=0;j<N;j++) if(fscanf(g,"%d",&a[j])!=1) return 1; for(int k=0;k<N;k++) if(fscanf(g,"%d",&b[k])!=1) return 1; fclose(g);}
  else {for(int j=0;j<N;j++) a[j]=(nxt()&1)?1:-1; for(int k=0;k<N;k++) b[k]=(nxt()&1)?1:-1;}
  for(int j=0;j<N;j++) for(int k=0;k<N;k++){A[j]+=H[j*N+k]*b[k]; B[k]+=H[j*N+k]*a[j];}
  long Om=0; for(int j=0;j<N;j++) Om+=a[j]*A[j];
  double sq=sqrt((double)N), n32=N*sq; double*acc=malloc((2*N+1)*sizeof(double));
  for(int t=0;t<nT;t++){
    double beta=1.0/T[t]; for(int y=-N;y<=N;y++) acc[y+N]=(y<=0)?1.0:exp(-2.0*beta*y/sq);
    double sm=0,sm2=0; long cnt=0, best=Om;
    for(long sw=0;sw<sweeps;sw++){
      for(int it=0;it<2*N;it++){
        int r=nxt()%(2*N);
        if(r<N){int j=r; long y=a[j]*A[j]; if(y<=0||unif()<acc[y+N]){ a[j]=-a[j]; Om-=2*y; int*h=H+j*N; int aj=a[j]; for(int k=0;k<N;k++) B[k]+=2*aj*h[k]; }}
        else {int k=r-N; long y=b[k]*B[k]; if(y<=0||unif()<acc[y+N]){ b[k]=-b[k]; Om-=2*y; int bk=b[k]; for(int j=0;j<N;j++) A[j]+=2*bk*H[j*N+k]; }}
        if(Om>best) best=Om;
      }
      if(sw>=sweeps/2){double m=Om/n32; sm+=m; sm2+=m*m; cnt++;}
    }
    double mean=sm/cnt, var=sm2/cnt-mean*mean;
    printf("%.4f %.6f %.6f %.6f\n",T[t],mean,var*N,best/n32); fflush(stdout);
  }
  int changed=1; while(changed){changed=0;
    for(int j=0;j<N;j++){long y=a[j]*A[j]; if(y<0){a[j]=-a[j];Om-=2*y;for(int k=0;k<N;k++)B[k]+=2*a[j]*H[j*N+k];changed=1;}}
    for(int k=0;k<N;k++){long y=b[k]*B[k]; if(y<0){b[k]=-b[k];Om-=2*y;for(int j=0;j<N;j++)A[j]+=2*b[k]*H[j*N+k];changed=1;}}}
  printf("quench %.6f\n",Om/n32);
  return 0;
}
