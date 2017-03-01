/** @file     */
#ifndef HUDSON_H__
#define HUDSON_H__

#include <iostream>

/// Wrapper for Husdon simulations
class Hudson {
 struct hapdata {
    int n;
    int n00;
    int n01;
    int n10;
    int n11;
    int n0q;
    int n1q;
    int nq0;
    int nq1;
    int nqq;
    double  d12;
    int ad ;
  };
public:
  /** Ctor and Dtor                   */
  Hudson(int nfiles, const char **fname);
  Hudson(const char *fname);
  ~Hudson();
  /** Load some data into the class   */
  void setdata(int **mat, double *pos,size_t nsam, size_t nsites);
  std::ostream& print(std::ostream &o); 
  double operator()(double rho) {
    return lnlikemshap(rho,0.0, 1.0, nsitesp,recrates,prob,ppoly);
  }
private:
  /** First the functions             */
  double lnlikemshap( double c,double conv,double conlen, int nsites
		      , double *recrates, double *****prob, double **ppoly  );
  double ****getprobmat(const char *fname, int *pnsam, int *pnsites, double **precrates
			, double **pppoly );
  double  getprobcond(int n, int n1, int n2, int n11, double r, int nsites, double *recrates, 
		      double ****prob, double *ppoly );
  double getprobcondu(int n, int n1, int n2, int n11, double r, int nsites, double *recrates, 
		       double ****prob, double *ppoly  );
  double getprobmscond(hapdata hv,  double c, int nsites, double *recrates, double ****prob, double *ppoly );
  double  extrap(  double p1, double p2, double r1, double r2, double r );
  void setup(int nfiles, const char **fname);

  hapdata *hv;
  int numfiles,*n;

  double conv, conlen;
  double *****prob, *recrates , **ppoly ;
  int nsitesp;

  size_t maxpair,npairs;
};


#endif
