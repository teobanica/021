/* dos.c -- for a Hadamard matrix H of order N <= 44, enumerates the 2^(N-1) sign vectors b with b_1 = 1 in Gray
   code order, keeps R = Hb up to date (O(N) per step), and counts the vectors b according to the multiset of the
   numbers |R_j|, encoded by the histogram h_v = #{j : |R_j| = 2v + (N mod 2)}.  The density of states of the array
   only depends on these counts, see glowlaw.py.  Usage: dos matrixfile > out.  Prints one line per multiset:
   count h_0 h_1 ... h_{N/2}. */
#include <stdio.h>
#include <stdlib.h>
typedef unsigned __int128 u128;
#define TS (1<<22)
static u128 keys[TS]; static long long cnts[TS]; static char used[TS];
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r"); int H[64*64],n=0,x; while(fscanf(f,"%d",&x)==1) H[n++]=x; fclose(f);
  int N=0; while(N*N<n) N++; int V=N/2+1;
  u128 pw[70]; pw[0]=1; for(int v=1;v<=V;v++) pw[v]=pw[v-1]*(u128)(N+1);
  int R[64],b[64]; for(int k=0;k<N;k++) b[k]=1;
  for(int j=0;j<N;j++){R[j]=0; for(int k=0;k<N;k++) R[j]+=H[j*N+k];}
  u128 key=0; for(int j=0;j<N;j++) key+=pw[abs(R[j])/2];
  long long total=1LL<<(N-1);
  for(long long i=0;i<total;i++){
    if(i){ int k=__builtin_ctzll(i)+1; b[k]=-b[k]; int d=2*b[k];
      for(int j=0;j<N;j++){ int o=R[j]; R[j]+=d*H[j*N+k]; key+=pw[abs(R[j])/2]; key-=pw[abs(o)/2]; } }
    unsigned long long hsh=(unsigned long long)(key^(key>>64))*0x9E3779B97F4A7C15ULL; unsigned long p=hsh>>(64-22);
    while(used[p]&&keys[p]!=key) p=(p+1)&(TS-1);
    if(!used[p]){used[p]=1;keys[p]=key;} cnts[p]++;
    if(((i+1)&((1LL<<26)-1))==0) fprintf(stderr,"%s %lld/%lld vectors done\n",argv[1],i+1,total);
  }
  for(long p=0;p<TS;p++) if(used[p]){ printf("%lld",cnts[p]); u128 k=keys[p]; for(int v=0;v<V;v++){ printf(" %d",(int)(k%(u128)(N+1))); k/=(u128)(N+1);} printf("\n"); }
  return 0;
}
