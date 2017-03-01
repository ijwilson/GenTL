#include "pact.H"

/***************************************************************************
 *
 *  utility function to get data in the correct format for SampleSNP
 *
 *****************************************************************************/
void getequalspacedpositionsanddistances(const TNT::Array1D<double> &distances
					 , TNT::Array1D<int> &pos
					 ,TNT::Array1D<double> &gap
					 , TNT::Array1D<double> &x)
{
  int n=pos.dim();
  double pmin=distances[0];
  double pmax=distances[distances.dim()-1]; 
  double trygap=(pmax-pmin)/(n+1.);
  int index=0;
  for (int i=0;i<n;i++) {
    double position=pmin+(i+1.)*trygap;
    for (;;) {
      assert(distances[index]<position);
      if (distances[index+1] >= position) break;
      index++;
      assert(index<distances.dim()-1);
    }
    if (fabs(distances[index]-position)<1E-8) {
      position -= 1E-6; // so not equal
    }
    pos[i]=index;
    gap[i]=position-distances[index];
    assert(gap[i]>0.0);
    x[i]=position;
  }
}
/****************************************************************************/
