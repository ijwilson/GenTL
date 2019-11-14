#include "mutmodel.h"

using namespace GenTL;

base randbase(rng &r) {
  double p=lr();
  if (p<0.25) return A;
  else if (p<0.5) return C;
  else if (p<0.75)  return G;
  return T; 
}
/**
 * mutate from above to below in time t for a Jukes-Cantor model
 */
// void JCmutation::mutate(std::vector<bases> &below, const std::vector<bases> above, 
// 		double time)
// {
//   assert(below.size()==above.size());
//   assert(below.size()==local_mu.size());

//   for (size_t i=0;i<below.size();i++) {
//     if (lr()<exp(-2.*time*local_mu[i])) below[i]=above[i];
//     else below[i]=randbase(lr);
//   }
//   return;
// }
