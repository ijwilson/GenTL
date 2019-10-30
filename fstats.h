#ifndef FSTATS_H_
#define FSTATS_H_

#include "tnt/tnt.h"
#include "utilcpp/tntutils.H"
#include "utilcpp/gsl_minimiser.H"
#include <map>
double log_dmulti_dirichlet(const int *x,const  double *alpha, int n);

class MND {  // multinomial dirichlet utility.
  /** The initial class assumes that the first population is the continent
   *  and the rest are islands (for simplicity, and because this is the model I
   *  am using for diseases)
   */
public:
  /** constructor -- we assume that we have x[allele][population]  
   * and create another view of the data  */
  MND(const TNT::Array2D<int> &x)
    :alleles(x.dim1()),npops(x.dim2()),
     view(npops,alleles),p(alleles),alpha(alleles) {
    double sum=0.0;
    for (int i=0;i<alleles;i++) {
      sum += double(x[i][0]);
      for (int j=0;j<npops;j++) view[j][i]=x[i][j];
    }
    for (int i=0;i<alleles;i++) p[i] = double(view[0][i])/sum;
  };
   
  double operator()(double F)  {
    double ret=0.0;
    double theta=1./F-1.;
    for (int i=0;i<alleles;i++) alpha[i]=p[i]*theta;
    for (int i=1;i<npops;i++) 
      ret += log_dmulti_dirichlet(view[i],alpha,alleles);
    return ret;
  }
private:
  int alleles,npops;
  TNT::Array2D<int> view;
  TNT::Array1D<double> p;
  TNT::Array1D<double> alpha;
};


double theta1(int **data,int n_pops,int alleles);
void pair_fst_onelocus(double **fst,int **data, int alleles,  int npops);
double gst_islandmainlandonelocus(int **data, int alleles);
double fst_islandmainlandonelocus(int **data, int alleles);


template <typename T, typename L>
TNT::Array1D<double> MNDirichletF(const TNT::Array2D<T> &data,const TNT::Array1D<L> &location, int &mono) 
{
  TNT::Array1D<double> f(data.dim2());
  int n=location.dim();  
  int npops;
  mono=0;
  TNT::Array1D<int> convloc=factor(location,npops);
  for (int site=0;site<data.dim2();site++) {
    std::map<T,int> mapping;
    int count=0;
    for (int i=0;i<n;i++) {
      if (mapping.find(data[i][site])==mapping.end()) {
	mapping[data[i][site]]=count++;
      }
    }
    size_t alleles=mapping.size();
    if (alleles==1) { // no variation
      f[site]=0.0;
      mono+=1;
    } else {
      //std::cout << "alleles: " << alleles << " npops:" << npops<<  std::endl;
      TNT::Array2D<int> x(alleles,npops,0);
      for (int i=0;i<n;i++) x[mapping[data[i][site]]][convloc[i]] +=1;
      for (size_t i=0;i<alleles;i++) {
	if (x[i][0]==0) { //not variable in base
	  x[i][0]=1; // make it variable- rough ans ready!!!!!
	}
      }
      //   std::cout << x;
      MND mnd(x);
      double xx=0.05,fxx,l,fl,u,fu;
      if (bracket_maximum0(mnd,xx, fxx, l,fl,u,fu,0.001, 0.999)) {
	//   std::cerr << "xx = " << xx << " "<< fxx <<" " 
	//	<< l << " " << fl << " "   
	//	<< u << " " << fu << std::endl; 
	//  std::cerr << x ;
	
	gsl_minimiser<MND> mndmax(mnd,xx,-fxx,l,-fl,u,-fu,false);
	f[site]=mndmax.minimise();
      } else if (fabs(xx-0.001)<1E-8) {
	f[site]=0.0;
      } else if (xx>=0.999) {
	f[site]=0.999;
      }
    }
    // std::cout << " " << f[site] << std::endl;
  }
  return f;
}

template<typename T, typename L>
TNT::Array1D<double> estimateF_WC(const TNT::Array2D<T> &data,const TNT::Array1D<L> &location)
{
  TNT::Array1D<double> f(data.dim2());
  int n=location.dim();  
  int npops;
  TNT::Array1D<int> convloc=factor(location,npops);
  
  for (int site=0;site<data.dim2();site++) {
    std::map<T,int> mapping;
    int count=0;
    for (int i=0;i<n;i++) {
      if (mapping.find(data[i][site])==mapping.end()) {
	mapping[data[i][site]]=count++;
      }
    }
  
    assert(mapping.size()==count);
   
    size_t alleles=mapping.size();
				  
    TNT::Array2D<int> x(alleles,npops,0);
    for (int i=0;i<n;i++) x[mapping[data[i][site]]][convloc[i]] +=1;
    f[site] = theta1(x,npops,alleles);
  }
  return f;
}



template<typename T>
double theta(T *data,const TNT::Array1D<int> &location)
{
  int n=data.dim1();  
  int npops= *std::max_element(&location[0], &location[n]);
  
 
  std::map<T,int> mapping;
  int count=0;
  for (int i=0;i<n;i++) {
    if (mapping.find(data[i])==mapping.end()) {
      mapping[data[i]]=count++;
    }
  }
  assert(mapping.size()==count);
  size_t alleles=mapping.size();
  
  TNT::Array2D<int> x(alleles,npops);
  for (int i=0;i<n;i++) x[mapping[data[i]]][location[i]] +=1;
  return theta1(x,npops,alleles);
}
#endif
