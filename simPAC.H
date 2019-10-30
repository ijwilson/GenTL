#ifndef SIMPAC_H__
#define SIMPAC_H__

#include <vector>

#include "utilcpp/haplotype.H"
#include "utilcpp/gsl_rand.H"
#include "utilcpp/GenEnums.H"


HaplotypeArray simPAC(const int nloci, const int n
		      , const double rho,rng &r);


HaplotypeArray simPACsubdiv(const int nloci, const vector<int> &n
			    , const double rho, const double f, rng &r);

/** class for lookup tables for mutation */
class PACmutate {
public:
  PACmutate(int mx,  rng &rnd): r(rnd),prob(mx) {
    double thetahat=0.0;
    for (int i=1;i<mx;i++) thetahat+= 1./double(i);
    thetahat=1./thetahat;
    for (int i=1;i<mx;i++) {
      prob[i]=0.5*thetahat/(thetahat+i);
    }
  }
  bool operator()(int k) {
    if (r() < prob[k]) return true;
    return false;
  }
  
private:
  rng &r;
  vector<double> prob;
};
/** class for lookup tables for recombination - the simplest version */
class PACrecombine {
public:
  PACrecombine(int mx, double rho, rng &rnd): r(rnd),prob(mx) {
    for (int i=1;i<mx;i++) {
      prob[i] = 1.-exp(-rho/double(i));
    }
  }
  bool operator()(int k) {
    if (r() < prob[k]) return true;
    return false;
  }
  
private:
  rng &r;
  vector<double> prob;
};

/** class for selecting the correct line to copy with subdivision   */
class PACline {
public:
  PACline(double lf, int npops,rng &rn) :    f(lf),r(rn),count(npops,0) {
  }
  add(int pop) {
    count[pop]+=1;
  }
  int choose(int pop) {
    vector<double> prob(npops);
    prob(
  }
  
    
private:
  rng &r;
  vector<int> count;
  double f;
  };
};

/** CountLines - a class to count the lines at all loci for a sample   */
class CountLines {
public:
  CountLines(int npops): x_(npops) {
  }
  /** Add a line to population pop                                                */
  void add(int pop, int line) {
    x_[pop].push_back(line);
  }
  /** Add a line to the general group                                             */
  void addgeneral(int line) {
    general_.push_back(line);
  }
  /** Add a recombination *before* locus from the pop to the General group        */
  void recombToGeneral(int locus, int pop) {
  }
  /** Add a recombination *before* locus from the General group to population pop */
  void recombFromGeneral(int locus, int pop) {
  }
  /** Count the lines at a locus in population pop                                */
  int operator(int locus, int pop) {
  }
  /** Count the lines at a locus in the general population                        */
  int operator(int locus) {
  }
private:
  /** x_ and general hold the positions of the lines for the local 
   * and general classes of haplotypes for the first locus                        */
  std::vector<std::deque<int> > x_;        // counts for the first position
  std::deque<int> general_;
  /** this vector details changes to the numbers of lines caused by recombinations*/
  std::vector<std::list<int > changes_;  // changes to the number of lines
  CountLines();
};
#endif
