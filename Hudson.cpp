#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "Hudson.H"

std::ostream& Hudson::print(std::ostream &o) 
  {
    o << npairs << " pairs of sites " << std::endl;
    for (size_t i=0;i<npairs;i++) 
      o << hv[i].d12 << " " << hv[i].n00 << " " << hv[i].n01 << " " 
	<< hv[i].n10 << " " << hv[i].n11 << std::endl;
    return o;
  }


/** setup the probability matrices                                       */
Hudson::Hudson(int nfiles, const char **fname){
  setup(nfiles,fname);
}
Hudson::Hudson(const char *fname){
  setup(1,&fname);
}

void Hudson::setup(int nfiles, const char **fname)
{
  numfiles=nfiles;maxpair=200000;
  prob = (double *****)malloc( (size_t)numfiles*sizeof( double ****) ) ;
  ppoly = (double **)malloc( (size_t)numfiles*sizeof( double * ) ) ;
  n = (int *)malloc( (size_t)numfiles*sizeof( int) ) ;
  for( int i=0; i<numfiles; i++) {
    prob[i] = getprobmat(fname[i], &(n[i]),  &nsitesp, &recrates , &(ppoly[i]) ) ; 
    if( (i>0) && ( n[i] > n[i-1] ) ) { 
      std::cerr << "likefiles in wrong order.\n"; 
      exit(0);
    }
  }
  /* allocate memory for the hv files                                     */
  hv  = (struct hapdata *) malloc( maxpair*sizeof( struct hapdata ) ) ; 
}


/** The destructor   */
    Hudson::~Hudson() {
      free(hv);
      free(n);
      free(ppoly);
      free(prob);
      // this is not a proper destructor -- need to get rid of prob mat
}
  /** convert array type data to hv data                        */
void Hudson::setdata(int **mat, double *pos,size_t nsam, size_t nsites) 
  {
    size_t s1,s2;
    npairs=0;
    for( s1=0; s1 < nsites-1; s1++){
      for( s2 = s1+1; s2 < nsites; s2++){
	if( (npairs+1)  >= maxpair ){
	  maxpair = npairs + 10 ;
	  hv = (struct hapdata *) realloc( hv , maxpair*sizeof( struct hapdata) ) ;
	}
	hv[npairs].n00 = hv[npairs].n01 =  hv[npairs].n10 = hv[npairs].n11 = 0;
	hv[npairs].n0q =  hv[npairs].n1q =  hv[npairs].nq0 =  0;
	hv[npairs].nq1 =  hv[npairs].nqq = 0;
	hv[npairs].n=nsam;
	for(size_t i=0; i<nsam; i++){
	  if( (mat[i][s1] == 0 )&&( mat[i][s2] == 0 ) ) hv[npairs].n00++ ; 
	  else if( (mat[i][s1] == 0 )&&( mat[i][s2] == 1 ) ) hv[npairs].n01++ ; 
	  else if( (mat[i][s1] == 1 )&&( mat[i][s2] == 0 ) ) hv[npairs].n10++ ; 
	  else if( (mat[i][s1] == 1 )&&( mat[i][s2] == 1 ) ) hv[npairs].n11++ ; 
	  else if( (mat[i][s1] == 0 )&&( mat[i][s2] == 3 ) ) hv[npairs].n0q++ ; 
	  else if( (mat[i][s1] == 1 )&&( mat[i][s2] == 3 ) ) hv[npairs].n1q++ ; 
	  else if( (mat[i][s1] == 3 )&&( mat[i][s2] == 0 ) ) hv[npairs].nq0++ ;
	  else if( (mat[i][s1] == 3 )&&( mat[i][s2] == 1 ) ) hv[npairs].nq1++ ; 
	  else if( (mat[i][s1] == 3 )&&( mat[i][s2] == 3 ) ) hv[npairs].nqq++ ; 
	}
	hv[npairs].d12=fabs(pos[s2]-pos[s1]);
	hv[npairs].ad=0;  // missing ancestors for now
	npairs++;
      }
    }
  }


  


#define MIN(x, y) ( (x)<(y) ? (x) : (y) )
#define MAX(x, y) ( (x)>(y) ? (x) : (y) )     


double Hudson::getprobcond(int n, int n1, int n2, int n11, double r, int nsites, double *recrates, 
	      double ****prob, double *ppoly )
  {
    int m1, m2,  k  ;
     
    m1 = MAX( n1, n2);
    m2 = MIN( n1, n2) ;

    for( k=0; k<nsites; k++){
      if( recrates[k] == r ) return(  prob[k][m1][m2][n11]/ppoly[k]) ;
      if( recrates[k] >r )  return( extrap(prob[k-1][m1][m2][n11]/ppoly[k-1],
					   prob[k][m1][m2][n11]/ppoly[k] , recrates[k-1],recrates[k],r));
    }
    return(  prob[nsites-1][m1][m2][n11]/ppoly[nsites-1]) ;

  }

  double
  lnval( double p2)
  {
    if( p2 > 0.0 ) return( log(p2 ) ) ;
    else return( -9999. ) ;
  }


double  Hudson::extrap(  double p1, double p2, double r1, double r2, double r )
  {
    double linextrap( double, double, double, double, double), lnval( double) ;

    if( (r1 == 0 ) && ( (p1 == 0) || (p2 == 0.0) ) ) {
      return(  linextrap(p1, p2, r1, r2, r)  ) ;
    }
    else if( r1 == 0.0 )  return( exp( linextrap( log(p1), log(p2), r1, r2, r )) ) ;
    else if( (p1 == 0.0) || ( p2 == 0.0 ) ) return(  linextrap( p1, p2, log(r1),log(r2), log(r)));  
    else return( exp(  linextrap( log(p1), log(p2), log(r1), log(r2), log(r) )) ) ;
  }


  double 
  linextrap(  double p1, double p2, double r1, double r2, double r )
  {
    return(  p1 + (p2-p1)*(r-r1)/(r2-r1) ) ;  

  }

  double  Hudson::getprobcondu(int n, int n1, int n2, int n11, double r, int nsites, double *recrates, 
	       double ****prob, double *ppoly  )
  {
    int m1, m2,  k  ;
    double nadprob( int, int, int, int, int, double ****);
       
    m1 = MAX( n1, n2);
    m2 = MIN( n1, n2) ;

    for( k=0; k<nsites; k++){
      if(recrates[k] == r)return( nadprob(n, k, m1, m2, n11, prob)/ppoly[k ] ) ;
      if(recrates[k] >r) return(extrap(nadprob(n, k-1, m1, m2, n11, prob)/ppoly[k -1] ,
				       nadprob(n,k,m1,m2,n11,prob)/ppoly[k] , recrates[k-1],
				       recrates[k],r));
    }
    return( nadprob(n,nsites-1,m1,m2,n11,prob)/ppoly[nsites-1] ) ;
  }

  double
  nadprob(int n,  int k, int m1, int m2, int n11, double ****prob )
  {
    int n00, n01, n10 ;
    double sum  ;

    n10 = m1 - n11 ;
    n01 = m2 - n11 ;
    n00 = n - n10 - n01 - n11 ;
    if( (n00 == n01) && ( n01 == n10 ) && ( n10 == n11 ) ) return( prob[k][m1][m2
									       ][n11] ) ;
    else if( (n00 == n01) && ( n10 == n11) )
      return(  prob[k][m1][m2][n11] + prob[k][MAX(m2,n-m1)][MIN(m2,n-m1)][n01] );
    else if( ( n00 == n10 ) && ( n01 == n11 ) )
      return( prob[k][m1][m2][n11] + prob[k][MAX(m1,n-m2)][MIN(n-m2,m1)][n10] );
    else if( ( n00 == n11 ) && ( n01 == n10 ) )
      return(  prob[k][m1][m2][n11] + prob[k][n-m2][n-m1][n00] ) ;
    else {
      sum = prob[k][m1][m2][n11] ;
      sum += prob[k][ MAX( m2, n-m1) ][ MIN( m2, n-m1)][n01] ;
      sum += prob[k][ MAX( m1, n-m2) ][ MIN( m1, n-m2)][n10] ;
      sum += prob[k][ n-m2][ n-m1][n00] ;
      return(  sum) ;
    }
  }   

double  Hudson::getprobmscond( Hudson::hapdata hv,  double c, int nsites, double *recrates, double ****prob, double *ppoly )
  {
    int n00, n01, n10, n11, n00_0q, n01_0q, n10_1q, n11_1q, n00_q0, n10_q0,
      n01_q1, n11_q1, n00_qq, n10_qq, n01_qq, n11_qq, n1dot, ndot1 ;
    double sum = 0.0, pr, coef , x ;
    float bico(int, int), factln( int) ;
    int hvtotal, n ;
    double co1, co2, co3, co4, r ;
   

    /* printf(" getprob:  %d %d %d %d %d %d %lf %d\n", hv.n, hv.n00, hv.n01, hv.n10, hv.n11, hv.n0q, hv.d12, hv.ad);
     */
    r =  c ;
    n = hv.n ;
    n00 = hv.n00 ;
    n01 = hv.n01 ;
    n10 = hv.n10 ;
    n11 = hv.n11 ;
    hvtotal = n00 + n01 + n10 + n11 ;
    x = factln( hvtotal ) - factln( n) - factln( hv.n00) - factln( hv.n01) -
      factln( hv.n10 ) - factln( hv.n11) ;
    for( n00_0q = 0 ; n00_0q <=  hv.n0q; n00_0q++) {
      n01_0q = hv.n0q - n00_0q ;
      co1 = bico( hv.n0q, n00_0q) ;
      for( n10_1q = 0 ; n10_1q <=  hv.n1q; n10_1q++) {
	n11_1q = hv.n1q - n10_1q ;
	co2 = bico( hv.n1q, n10_1q);
	for( n00_q0 = 0 ; n00_q0 <=  hv.nq0; n00_q0++) {
	  n10_q0 = hv.nq0 - n00_q0 ;
	  co3 = bico( hv.nq0, n00_q0);
	  for( n01_q1 = 0 ; n01_q1 <=  hv.nq1; n01_q1++) {
	    n11_q1 = hv.nq1 - n01_q1 ;
	    co4 = bico( hv.nq1, n01_q1);
	    for( n00_qq = 0 ; n00_qq <= hv.nqq ; n00_qq++ ){
	      for( n01_qq = 0; n01_qq <= hv.nqq - n00_qq ; n01_qq++){
		for( n10_qq = 0; n10_qq <= hv.nqq - n00_qq - n01_qq ; n10_qq++){
		  n11_qq = hv.nqq - n00_qq - n01_qq - n10_qq ;
		  n00 = hv.n00 +  n00_0q +  n00_q0 + n00_qq ;
		  n01 =  hv.n01 + n01_0q + n01_q1 + n01_qq ;
		  n10 = hv.n10 + n10_1q + n10_q0 + n10_qq ;
		  n11 = hv.n11 + n11_1q + n11_q1 + n11_qq ;
		  n1dot = n10 + n11 ;
		  ndot1 = n01 + n11 ;

		  if( hv.ad == 1 )
		    pr = getprobcond( n, n1dot, ndot1, n11, r, nsites, recrates, prob, ppoly) ;
		  else
		    pr = getprobcondu( n, n1dot, ndot1, n11, r, nsites, recrates, prob, ppoly) ;
	
		  coef = x + factln( n00) + factln( n01) + factln(n10) + factln( n11) +
		    factln( hv.nqq) - factln( n00_qq) - factln( n01_qq) - factln( n10_qq) -
		    factln( n11_qq) ;
		  coef = exp( coef) *co1 *co2 * co2 * co3 * co4 ; 
		  sum += pr*coef  ;
		} 
	      } 
	    } 
	  } 
	} 
      } 
    }
    return( sum  ) ;
  }
 

  double factln(int n)
  {
    // double gammln(double xx);
    static double a[101];

    if (n < 0){
      std::cerr << " Negative factorial in routine factln\n";
      exit(0);
    }
    if (n <= 1) return 0.0;
    if (n <= 100) return a[n] ? a[n] : (a[n]=lgamma(n+1.0));
    else return lgamma(n+1.0);
  }

  double bico(int n,int k)
  {
    return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
  }


//   double gammln(double xx)
//   {
//     double x,y,tmp,ser;
//     static double cof[6]={76.18009172947146,-86.50532032941677,
// 			  24.01409824083091,-1.231739572450155,
// 			  0.1208650973866179e-2,-0.5395239384953e-5};
//     y=x=xx;
//     tmp=x+5.5;
//     tmp -= (x+0.5)*log(tmp);
//     ser=1.000000000190015;
//     for (int j=0;j<=5;j++) ser += cof[j]/++y;
//     return -tmp+log(2.5066282746310005*ser/x);
//   }



double Hudson::lnlikemshap( double c,double conv,double conlen, int nsites, double *recrates, double *****prob, double **ppoly  )
  {
    double sum=0.0, r ;
    size_t i,j;
   
    for( i =0 ; i< npairs; i++){
      r = c*( hv[i].d12  + 2.*conv*conlen*(1.0 - exp( -hv[i].d12/conlen) ))  ; 
      for( j=numfiles-1;  n[j] != hv[i].n ; j--);
      sum += lnval(  getprobmscond(hv[i], r, nsites,recrates, prob[j], ppoly[j] )) ;
    }

    return( sum ) ;
  }



double ****Hudson::getprobmat(const char *fname, int *pnsam, int *pnsites, double **precrates, double **pppoly )
{

  int  i ; 
  FILE *pf;
  int    j, k, m  ;
  double ****prob;
  int nsam, nsites, npop, *config, pop;
  double alphag;
  double  y , *pp ;
  int d1, d2 , ng ;

  pf = fopen( fname, "r");
  if( pf == NULL ) { fprintf(stderr,"No such hxxrho file.\n"); exit(1); } 

  fscanf(pf, " %d", &npop);
  config = (int *)malloc((unsigned)((npop+1)*sizeof(int)));
  for( pop=0  ; pop<npop; pop++) {
    fscanf(pf, " %d", config+pop ) ;
    if( config[pop] < 0 ) break;
  }
  fscanf(pf," %d", &nsites) ;
  *pnsites = nsites ;
  *precrates = (double *)malloc( (unsigned)nsites*sizeof(double) ) ;
  pp  = (double *)malloc( (unsigned)nsites*sizeof(double) ) ;
  *pppoly = pp ;
  for( i=0; i<nsites; i++){
    fscanf(pf," %lf", *precrates + i ) ;
  }
  fscanf(pf," %lf", &alphag ) ;
  if( alphag != 0.0 ) fscanf(pf, " %*lf" ) ;
  if( npop > 1 )  fscanf(pf," %*lf" ) ;
  fscanf(pf, " %*d" ) ;
  for(i=0;i<3;i++){
    fscanf(pf," %*d");
  }

  while( pop<npop ){ config[pop++]=0 ;  }
  nsam = 0 ;
  for(i=0;i<npop; i++) nsam += config[i] ;
  *pnsam = nsam ;

  prob = (double ****)malloc( (unsigned)nsites*sizeof( double *** ) ) ;
  for( i=0; i<nsites; i++){
    prob[i] = (double ***)malloc( (unsigned)nsam*sizeof( double **) ) ;
    for( j=0; j<nsam; j++) {
      prob[i][j] = (double **)malloc( (unsigned)(j+1)*sizeof( double *) ) ;
      for( k=0; k<=j ; k++) {

	prob[i][j][k] = (double *)malloc( (unsigned)( MIN(j,k) - MAX( 0, (j+k-nsam) )+1 )*sizeof( double) ) ;
	if( prob[i][j][k] == 0 ) printf( " whoaa!\n" ) ;
	prob[i][j][k] -= MAX( 0, j+k-nsam ) ;
      }
    }
  }

  for( m=0; m<nsites; m++) pp[m] = 0.0 ;

  for( i=1; i<nsam; i++)
    for( j=1; j<= i ; j++){
      fscanf(pf," freq: %d %d", &d1, &d2) ;
      for( m=0; m<nsites; m++){
	ng =  fscanf(pf," rho: %lf", &y) ;
	for( k= MAX(0,(i+j-nsam))  ; k<=  MIN(i,j)  ; k++) {
	  fscanf(pf," %le", &( prob[m][i][j][k]) ) ;
	  if( i != j ) pp[m] += 2.0*prob[m][i][j][k] ;
	  else pp[m] += prob[m][i][j][k] ;
	}

      }
    }
  fclose(pf) ;

  return( prob ) ;
}



