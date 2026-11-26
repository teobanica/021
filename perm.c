/* perm.c -- the permanent of a (+1,-1)-matrix H of order N <= 32, read from a file, by the formula of Glynn,
   per(H) = 2^(1-N) sum over b in {-1,1}^N with b_1 = 1 of (prod_k b_k) (prod_j (Hb)_j), the vectors b being
   enumerated in Gray code order, with Hb updated in O(N).  The sum is computed in 128-bit arithmetic modulo 2^128,
   which is exact since its true value 2^(N-1) per(H) is much smaller than 2^127 for Hadamard matrices of order 24.
   Usage: ./perm matrix.txt, prints per(H). */
#include <stdio.h>
#include <stdlib.h>
typedef unsigned __int128 u128;
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r"); int H[32][32],N=0,x,c=0; long R[32];
  int vals[1024]; while(c<1024&&fscanf(f,"%d",&x)==1) vals[c++]=x; fclose(f);
  for(N=1;N*N<c;N++);
  for(int i=0;i<N;i++)for(int j=0;j<N;j++)H[i][j]=vals[i*N+j];
  for(int j=0;j<N;j++){R[j]=0;for(int k=0;k<N;k++)R[j]+=H[j][k];}   /* b = all ones */
  int b[32]; for(int k=0;k<N;k++)b[k]=1; int sg=1; u128 S=0;
  unsigned long long M=1ULL<<(N-1);
  for(unsigned long long t=0;;t++){
    u128 p=(u128)(__int128)sg; for(int j=0;j<N;j++) p*=(u128)(__int128)R[j];
    S+=p;
    if(t+1==M)break;
    int k=__builtin_ctzll(t+1)+1;               /* flip coordinate k (1..N-1), keeping b_1 = 1 */
    for(int j=0;j<N;j++) R[j]-=2*b[k]*H[j][k];
    b[k]=-b[k]; sg=-sg;
  }
  __int128 v=(__int128)S; v>>=(N-1);              /* exact division by 2^(N-1) */
  long long r=(long long)v; printf("%lld\n",r); return 0;
}
