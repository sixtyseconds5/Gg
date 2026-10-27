// clause-weighted local search (PAWS-like, clause-focused). usage: wls file.cnf seed maxsec [init.txt] [noise] [smooth]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
static unsigned long long rs; static inline unsigned rnd(){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return (unsigned)(rs>>11);}
int main(int argc,char**argv){
 FILE*f=fopen(argv[1],"r"); rs=atoll(argv[2])*2654435761ULL+12345; double maxs=atof(argv[3]);
 double noise=argc>5?atof(argv[5]):0.005, sp=argc>6?atof(argv[6]):0.01;
 unsigned noiseT=(unsigned)(noise*16777216), spT=(unsigned)(sp*16777216);
 int nv,nc; char buf[512]; while(fgets(buf,512,f)) if(buf[0]=='p'){sscanf(buf,"p cnf %d %d",&nv,&nc);break;}
 int *cs=malloc(4*(nc+1)),*lits=malloc(4*nc*8),n=0;
 for(int i=0;i<nc;i++){cs[i]=n;int l;while(fscanf(f,"%d",&l)==1&&l)lits[n++]=l;} cs[nc]=n;
 int *occn=calloc(2*nv+2,4); for(int i=0;i<n;i++)occn[lits[i]+nv]++;
 int *os=malloc(4*(2*nv+3)); os[0]=0; for(int j=0;j<=2*nv;j++)os[j+1]=os[j]+occn[j];
 int *occ=malloc(4*n),*fl=calloc(2*nv+2,4); for(int i=0;i<nc;i++)for(int j=cs[i];j<cs[i+1];j++){int x=lits[j]+nv;occ[os[x]+fl[x]++]=i;}
 char*a=malloc(nv+1),*ba=malloc(nv+1); for(int v=1;v<=nv;v++)a[v]=rnd()&1;
 if(argc>4&&strcmp(argv[4],"-")){FILE*g=fopen(argv[4],"r");if(g){for(int v=1;v<=nv;v++){int c=fgetc(g);if(c=='0'||c=='1')a[v]=c-'0';}fclose(g);}}
 int *nt=calloc(nc,4),*cr=calloc(nc,4),*w=malloc(4*nc),*un=malloc(4*nc),*pos=malloc(4*nc),nun=0;
 long long *wb=calloc(nv+1,8),*wm=calloc(nv+1,8);
 #define T(l) ((l)>0?a[l]:!a[-(l)])
 for(int i=0;i<nc;i++){w[i]=1;for(int j=cs[i];j<cs[i+1];j++)if(T(lits[j])){nt[i]++;cr[i]+=abs(lits[j]);}
  if(nt[i]==1)wb[cr[i]]+=1; if(!nt[i]){pos[i]=nun;un[nun++]=i;for(int j=cs[i];j<cs[i+1];j++)wm[abs(lits[j])]+=1;}}
 int *big=malloc(4*nc),nbig=0,*bpos=malloc(4*nc); for(int i=0;i<nc;i++)bpos[i]=-1;
 clock_t t0=clock(); long long it=0; int best=nun; memcpy(ba,a,nv+1);
 while(nun){
  if((++it&0xfffff)==0){ if((clock()-t0)>maxs*CLOCKS_PER_SEC)break; }
  int c=un[rnd()%nun]; int k=cs[c+1]-cs[c]; int v=0; long long bs=-(1LL<<62);
  for(int j=0;j<k;j++){int u=abs(lits[cs[c]+j]); long long s=wm[u]-wb[u]; if(s>bs||(s==bs&&(rnd()&1))){bs=s;v=u;}}
  if(bs<=0){
   if((rnd()&0xFFFFFF)<noiseT){ v=abs(lits[cs[c]+rnd()%k]); }
   else {
    // bump weights of unsat clauses
    for(int q=0;q<nun;q++){int d=un[q]; w[d]++; if(bpos[d]<0){bpos[d]=nbig;big[nbig++]=d;} for(int j=cs[d];j<cs[d+1];j++)wm[abs(lits[j])]++;}
    if((rnd()&0xFFFFFF)<spT){ // smooth: decrement all weights>1
     for(int q=0;q<nbig;){int d=big[q]; w[d]--; if(nt[d]==1)wb[cr[d]]--; else if(nt[d]==0)for(int j=cs[d];j<cs[d+1];j++)wm[abs(lits[j])]--;
      if(w[d]==1){bpos[d]=-1; big[q]=big[--nbig]; if(q<nbig)bpos[big[q]]=q;} else q++;}
    }
    continue;
   }
  }
  int tl= a[v]? -v: v;
  for(int q=os[tl+nv];q<os[tl+nv+1];q++){int d=occ[q],wd=w[d];
   if(nt[d]==0){int p=pos[d];un[p]=un[--nun];pos[un[p]]=p; for(int j=cs[d];j<cs[d+1];j++)wm[abs(lits[j])]-=wd; wb[v]+=wd;}
   else if(nt[d]==1) wb[cr[d]]-=wd;
   nt[d]++;cr[d]+=v;}
  for(int q=os[-tl+nv];q<os[-tl+nv+1];q++){int d=occ[q],wd=w[d];
   nt[d]--;cr[d]-=v;
   if(nt[d]==0){wb[v]-=wd; pos[d]=nun;un[nun++]=d; for(int j=cs[d];j<cs[d+1];j++)wm[abs(lits[j])]+=wd;}
   else if(nt[d]==1) wb[cr[d]]+=wd;}
  a[v]=!a[v]; if(nun<best){best=nun;memcpy(ba,a,nv+1);}
 }
 double el=(double)(clock()-t0)/CLOCKS_PER_SEC;
 fprintf(stderr,"c steps %lld time %.2f best %d final %d\n",it,el,best,nun);
 FILE*o=fopen("wls_best.txt","w"); for(int v=1;v<=nv;v++)fputc('0'+ba[v],o); fclose(o);
 if(!nun){printf("s SATISFIABLE\nv");for(int v=1;v<=nv;v++)printf(" %d",a[v]?v:-v);printf(" 0\n");}
 return nun?1:0;
}
