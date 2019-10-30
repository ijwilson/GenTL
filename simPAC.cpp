#include <numeric>                         // for accumulate
#include <deque>
#include "simPAC.H"
#include "tnt/tnt.h"

/** simulate using a PAC likelihood.  Note that this is for a recombination rate
 * of rho between equally spaced markers and for theta as defined in 
 * Li and Stephens
 */
HaplotypeArray simPAC(const int nloci, const int n,  const double rho, rng &r)
{
  HaplotypeArray a(n,nloci);

  PACmutate pm(n,r);
  PACrecombine pr(n,rho,r);

  // the first chromosome
  for (int i=0;i<nloci;i++) a[0][i]=zero;
// {
//     if (r()<0.5) a[0][i]=one;
//     else a[0][i]=zero;
//   }
  // the second
  if (pm(1)) a[1][0]=flip(a[0][0]);
  else a[1][0]=a[0][0];
   
  for (int k=1;k<nloci;k++) {
    if (pm(1)) a[1][k]=flip(a[0][k]);
    else a[1][k]=a[0][k];
  }
  // and the rest
  for (int j=2;j<n;j++) {
    int copychromosome=r.rint(j-1); 
    if (pm(j)) a[j][0]=flip(a[copychromosome][0]);
    else a[j][0]=a[copychromosome][0];
   
    for (int k=1;k<nloci;k++) {
      if (pr(j)) copychromosome=r.rint(j-1);
      if (pm(j)) a[j][k]=flip(a[copychromosome][k]);
      else a[j][k]=a[copychromosome][k];
    }  
  }
  return a;
}

/** generate from my hacked together PAC                             
 *  use a random ordering of samples                                 */
HaplotypeArray simPACsubdiv(const int nloci, const vector<int> &n
			    , const double rho, const double f, rng &r)
{
  int total=std::accumulate(n.begin(),n.end(),0);
  HaplotypeArray a(total,nloci);
 
  // use a random ordering (not sure if needed)
  vector<int> order=r.integer_permutations(total,total);
  // what are the populations for the samples
  vector<int> where(total);
  int count=0;
  for (int i=0;i<n.size();i++) {
    for (int j=0;j<n[i];j++) where[count++]=i;
  }
  // and where are they in HaplotypeArray a
  vector<int> start=n;
  for (int i=1;i<n.size;i++) start[i]+=start[i-1];

  PACmutate pm(total,r);
  PACrecombine pr(total,rho,r);
  PACline pl(f,n.size(),r);
  CountLines cl(int npops);

  std::deque<int> general;
  std::vector<std::deque<int> > local(n.size());
  
  // the first chromosome is in the general population
  for (int i=0;i<nloci;i++) {
    if (r()<0.5) a[order[0]][i]=one;
    else a[order[0]][i]=zero;
  }
  general.push_pack(order[0]);
  pl.add(where[order[0]]);
  // the second chromosome - is it a copy
  // or from the "general" population
  int nxt=order[1]; 
  pair<bool,int> which=pl(where[nxt]);
  if (which->first) { // general population
    copychromosome=general[wh.second];
  } else { //copy from population
    copychromosome=local[where[nxt]][wh.second];
  }
  


  // the second - copy with mutations from the first
  if (pm(1)) a[nxt][0]=flip(a[use][0]);
  else a[nxt][0]=a[use][0];
   
  for (int k=1;k<nloci;k++) {
    if (pm(1)) a[nxt][k]=flip(a[use][k]);
    else a[nxt][k]=a[use][k];
  }
  // and the rest 
  for (int j=0;j<n;j++) {
    int copychromosome=pl(); 
    if (pm(j)) a[j][0]=flip(a[copychromosome][0]);
    else a[j][0]=a[copychromosome][0];
   
    for (int k=1;k<nloci;k++) {
      if (pr(j)) copychromosome=r.rint(j-1);
      if (pm(j)) a[j][k]=flip(a[copychromosome][k]);
      else a[j][k]=a[copychromosome][k];
    }  
  }
  }

  int count=2;
  for (int pop=0;pop<n.size();i++) {
    for (int j=0;j<n;j++) {
    int copychromosome=pl(); 
    if (pm(j)) a[j][0]=flip(a[copychromosome][0]);
    else a[j][0]=a[copychromosome][0];
   
    for (int k=1;k<nloci;k++) {
      if (pr(j)) copychromosome=r.rint(j-1);
      if (pm(j)) a[j][k]=flip(a[copychromosome][k]);
      else a[j][k]=a[copychromosome][k];
    }  
  }
  }
  return a;
}



/** generate from my hacked together PAC                             
 *  use a random ordering of samples                                 */
HaplotypeArray simPACsubdiv(const int nloci, const vector<int> &n
			    , const double rho, const double f, rng &r)
{
  int total=std::accumulate(n.begin(),n.end(),0);
  HaplotypeArray a(total,nloci);
  TNT::Array2D<int> copyfrom(a.dim1(),a.dim2());
 
  // use a random ordering (not sure if needed)
  vector<int> order=r.integer_permutations(total,total);
  // what are the populations for the samples
  vector<int> where(total);
  int count=0;
  for (int i=0;i<n.size();i++) {
    for (int j=0;j<n[i];j++) where[count++]=i;
  }
  // and where are they in HaplotypeArray a
  vector<int> start=n;
  for (int i=1;i<n.size;i++) start[i]+=start[i-1];

  PACmutate pm(total,r);
  PACrecombine pr(total,rho,r);
 
  // the first chromosome is in the general population
  for (int i=0;i<nloci;i++) {
    if (r()<0.5) a[0][i]=one;
    else a[0][i]=zero;
    copyfrom[0][i]=where[order[0]]; 
  }
  // the second chromosome - is it a copy
  // or from the "general" population
  int nxtpop=where[order[1]];
  int nlocal=(copyfrom[0][0]==nxtpop); // just for line 2

  


  


  pair<bool,int> which=pl(where[nxt]);
  if (which->first) { // general population
    copychromosome=general[wh.second];
  } else { //copy from population
    copychromosome=local[where[nxt]][wh.second];
  }
  


  // the second - copy with mutations from the first
  if (pm(1)) a[nxt][0]=flip(a[use][0]);
  else a[nxt][0]=a[use][0];
   
  for (int k=1;k<nloci;k++) {
    if (pm(1)) a[nxt][k]=flip(a[use][k]);
    else a[nxt][k]=a[use][k];
  }
  // and the rest 
  for (int j=0;j<n;j++) {
    int copychromosome=pl(); 
    if (pm(j)) a[j][0]=flip(a[copychromosome][0]);
    else a[j][0]=a[copychromosome][0];
   
    for (int k=1;k<nloci;k++) {
      if (pr(j)) copychromosome=r.rint(j-1);
      if (pm(j)) a[j][k]=flip(a[copychromosome][k]);
      else a[j][k]=a[copychromosome][k];
    }  
  }
  }

  int count=2;
  for (int pop=0;pop<n.size();i++) {
    for (int j=0;j<n;j++) {
    int copychromosome=pl(); 
    if (pm(j)) a[j][0]=flip(a[copychromosome][0]);
    else a[j][0]=a[copychromosome][0];
   
    for (int k=1;k<nloci;k++) {
      if (pr(j)) copychromosome=r.rint(j-1);
      if (pm(j)) a[j][k]=flip(a[copychromosome][k]);
      else a[j][k]=a[copychromosome][k];
    }  
  }
  }
  return a;
}
