// probSAT for k-SAT. usage: psk file.cnf seed maxsec cb [init.txt]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
static unsigned long long rs; static inline unsigned rnd(){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return (unsigned)(rs>>11);}
int main(int argc,char**argv){
 FILE*f=fopen(argv[1],"r"); rs=atoll(argv[2])*2654435761ULL+12345; double maxs=atof(argv[3]),cb=atof(argv[4]);
 int nv,nc; char buf[512]; while(fgets(buf,512,f)) if(buf[0]=='p'){sscanf(buf,"p cnf %d %d",&nv,&nc);break;}
 int *cs=malloc(4*(nc+1)),*lits=malloc(4*nc*8),n=0;
 for(int i=0;i<nc;i++){cs[i]=n;int l;while(fscanf(f,"%d",&l)==1&&l)lits[n++]=l;} cs[nc]=n;
 int *occn=calloc(2*nv+2,4); for(int i=0;i<n;i++)occn[lits[i]+nv]++;
 int *os=malloc(4*(2*nv+3)); os[0]=0; for(int j=0;j<=2*nv;j++)os[j+1]=os[j]+occn[j];
 int *occ=malloc(4*n),*fl=calloc(2*nv+2,4); for(int i=0;i<nc;i++)for(int j=cs[i];j<cs[i+1];j++){int x=lits[j]+nv;occ[os[x]+fl[x]++]=i;}
 char*a=malloc(nv+1); for(int v=1;v<=nv;v++)a[v]=rnd()&1;
 if(argc>5){FILE*g=fopen(argv[5],"r");if(g){for(int v=1;v<=nv;v++){int c=fgetc(g);if(c=='0'||c=='1')a[v]=c-'0';}fclose(g);}}
 int *nt=calloc(nc,4),*cr=calloc(nc,4),*br=calloc(nv+1,4),*un=malloc(4*nc),*pos=malloc(4*nc),nun=0;
 #define T(l) ((l)>0?a[l]:!a[-(l)])
 for(int i=0;i<nc;i++){for(int j=cs[i];j<cs[i+1];j++)if(T(lits[j])){nt[i]++;cr[i]+=abs(lits[j]);} if(nt[i]==1)br[cr[i]]++; if(!nt[i]){pos[i]=nun;un[nun++]=i;}}
 double pw[256]; for(int k=0;k<256;k++)pw[k]=pow(cb,-k);
 clock_t t0=clock(); long long it=0; int best=nun;
 while(nun){
  if((++it&0xfffff)==0){ if((clock()-t0)>maxs*CLOCKS_PER_SEC)break; }
  int c=un[rnd()%nun]; double w[16],s=0; int k=cs[c+1]-cs[c];
  for(int j=0;j<k;j++){int b=br[abs(lits[cs[c]+j])];w[j]=pw[b>255?255:b];s+=w[j];}
  double r=s*((rnd()&0xFFFFFF)/16777216.0); int j=0; while(j<k-1&&r>=w[j]){r-=w[j];j++;}
  int v=abs(lits[cs[c]+j]); int tl= a[v]? -v: v;
  for(int q=os[tl+nv];q<os[tl+nv+1];q++){int d=occ[q]; if(nt[d]==0){int p=pos[d];un[p]=un[--nun];pos[un[p]]=p;} else if(nt[d]==1)br[cr[d]]--; nt[d]++;cr[d]+=v; if(nt[d]==1)br[cr[d]]++;}
  for(int q=os[-tl+nv];q<os[-tl+nv+1];q++){int d=occ[q]; if(nt[d]==1)br[cr[d]]--; nt[d]--;cr[d]-=v; if(nt[d]==0){pos[d]=nun;un[nun++]=d;} else if(nt[d]==1)br[cr[d]]++;}
  a[v]=!a[v]; if(nun<best)best=nun;
 }
 double el=(double)(clock()-t0)/CLOCKS_PER_SEC;
 fprintf(stderr,"c flips %lld time %.2f best %d final %d\n",it,el,best,nun);
 if(!nun){printf("s SATISFIABLE\nv");for(int v=1;v<=nv;v++)printf(" %d",a[v]?v:-v);printf(" 0\n");}
 return nun?1:0;
}
