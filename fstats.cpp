#include <cstdlib>
#include <cstdio>
#include <map>
#include "tnt/tnt.h"
#include "utilcpp/utilityfunctionals.h"

lfactorl lfactrl;

/************************************************************************/
void avevar(double *data,int n,double *ave, double *var)
{
	int j;
	double s,ep;
	
	for (*ave=0.0,j=0;j<n;j++) *ave += data[j];
	*ave /= n;
	*var=ep=0.0;
	for (j=0;j<n;j++)	{
		s=data[j] - *ave;
		ep +=s;
		*var += s*s;
	}
	*var = (*var-ep*ep/n)/(n-1);
}
/************************************************************************/
void n_c_calc(int *data, int n,double *ave, double *n_c)

{
	double s2=0.0;
	*ave=0.0;

	for (int j=0;j<n;j++) {
	  *ave += (double)data[j];
	  s2+=(double)data[j]*(double)data[j];
	}
	*n_c = (*ave-s2/(*ave))/((double)n-1.0);
	*ave /= (double)n;
}


/****************************************************************************/
/* note that I have transposed the original as my data have alleles in rows
 * and populations in columns                 

The data look like

*/
double theta1(int **data,int n_pops,int alleles)
/*Theta for one locus*/
{
  int allele,pop,i,j;
  double nbar,nc,pbar,ss,atot,dtot;
	
  TNT::Array1D<int> n(n_pops,0);

  for (i=0;i<n_pops;i++) {
    for (j=0;j<alleles;j++) {
      n[i] += data[j][i];
    }
  }
  n_c_calc(n,n_pops,&nbar,&nc);
	
  TNT::Array2D<double> f(alleles,n_pops);
  TNT::Array1D<double> a(alleles);
  TNT::Array1D<double> d(alleles);
  for (allele=0;allele<alleles;allele++){
    for (pop=0;pop<n_pops;pop++)
      f[allele][pop] = (double)data[allele][pop]/(double)n[pop];
	  
    avevar(f[allele],n_pops,&pbar,&ss);	
    a[allele] = (nbar/nc)*(ss - 1.0/(nbar-1.0)*
			   (pbar*(1.0-pbar)-(n_pops-1.0)/((double)n_pops)*ss));
    d[allele] = nbar/(nbar-1.0)*
      (pbar*(1.0-pbar)-((double)n_pops-1.0)/(double)n_pops * ss);
  }
	
  atot=0.0;dtot=0.0;
    
  for (allele=0;allele<alleles;allele++){
    atot+= a[allele];
    dtot+= d[allele];
  }
  return atot/(atot+dtot);
}


/****************************************************************************/
/* note that I have transposed the original as my data have alleles in rows
 * and populations in columns                 

The data look like

mainland[allele0] mainland[allele1] ...
island[allele0] ...

*/


/************************************************************************/
/*  note data is 

data[pop1][allele1] data[pop2][allele1] data[pop3][allele1] ...
data[pop1][allele2] ..........                                          */

void pair_fst_onelocus(double **fst,int **data, int alleles,  int npops)
{
  int i,j,k;
  
  TNT::Array2D<int> d(alleles,2);
 
  for (i=0;i<npops;i++) {
    for (j=i+1;j<npops;j++) {
       for (k=0;k<alleles;k++) {
	 d[k][0]=data[k][i];
	 d[k][1]=data[k][j];
      }
      fst[i][j]=theta1(d,2,alleles);
    }
  }

  return;
}
/************************************************************************/
/***************************************************************************/
double log_D(const double *b, int n)
{
  double sum=0.0,temp=0.0;
  
  for (int i=0;i<n;i++) {
    sum += b[i];
    temp += lgamma(b[i]);
  }
  temp -= lgamma(sum);
  return temp;
}
/**********************************************************************/
/*    The log density for the multinomial-Dirichlet distribution for 
      vector of counts x, where x[i] is the number of observations i 
	  with parameters alpha, a vector of length k                     */
/**********************************************************************/ 
double log_dmulti_dirichlet(const int *x, const double *alpha, int n)
{
  int i,tot=0;
  double temp;
  
  double *ax = new double[n];
  temp = -log_D(alpha,n);
  for (i=0;i<n;i++) {
    tot+=x[i];
    ax[i]= alpha[i]+(double)x[i];
    temp -= lfactrl(x[i]);
  }
  temp += lfactrl(tot)+log_D(ax,n);
  delete [] ax;
  return temp;
}
/***************************************************************************/
double lddirichlet(double *x, double a, int n)
{
  /* expected frequencies   */
  double *alpha = new double[n];
  
  double lx=0.0;
  for (int i=0;i<n;i++) {
    alpha[i]=a;
    lx += (a-1.0)*log(x[i]);
  }
  lx -= log_D(alpha,n);
  delete[] alpha;
  return lx;
}

/************************************************************************/
double gst_islandmainlandonelocus(int **data, int alleles)
{
	
  TNT::Array1D<int> n(2,0);
  
  for (int j=0;j<alleles;j++) {
    n[0] += data[j][0];
    n[1] += data[j][1];
  }
 
  

  TNT::Array1D<double> ht(alleles);
  TNT::Array1D<double> dst(alleles);

  for (int l=0;l<alleles;l++) {
    double pmainland =double(data[l][0])/double(n[0]);
    ht[l]= 2.*pmainland*(1.-pmainland);
    double p =(double)data[l][1]/(double)n[1];
    dst[l] = 2.*p*(1.-p);
  }   
  double F=0.0;
  for (int i=0;i<alleles;i++) 
    F += 1.- dst[i]/ht[i];
    
  return F/double(alleles);
}
/************************************************************************/
double fst_islandmainlandonelocus(int **data, int alleles)
{
	
  TNT::Array1D<int> n(2,0);
  


  for (int j=0;j<alleles;j++) {
    n[0] += data[j][0];
    n[1] += data[j][1];
  }

  double sumsq=0.0; 

  for (int l=0;l<alleles;l++) {
    double pmainland =double(data[l][0])/double(n[0]);
    double p =(double)data[l][1]/(double)n[1];
    double pmean=(pmainland+p)/2.0;
      
      sumsq += (p-pmainland)*(p-pmainland)/(pmean*(1.-pmean));
  }   
    
  return sumsq/double(alleles-1);
}
