#ifndef FLOWGRAM_H__
#define FLOWGRAM_H__

#include <vector>
#include <string>
#include <iosfwd>

#include "fasta.H"

/** A class to hold flowgram information                                          
 * This class holds the flowgram information for a single sequence     */
class flowGram {
 public:
  /** Constructors and destructors  */
  flowGram(const std::vector<std::string> &lines);
  ~flowGram(){};
  /** print information from a flowGram    */
  std::ostream &print(std::ostream &o,const std::string &gap=",") const;
  size_t length() const {
    return intensity.size();
  }
// The data
  std::vector<char> Bases;
  std::vector<double> intensity;  // the intensity of the flowGram peak
  std::string name,run;
  int x,y,region;
};


/** A custom class to create a predicate function for the 454  
 * That allows us to filter flowGrams                                 */
class flowGramFilter {
public:
  flowGramFilter(double minrun, int minlength):
    minRun(minrun),minLength(minlength) {};
  bool operator()(const flowGram &x) {
    if (x.length()<static_cast<size_t>(minLength)) return true;
    int count=0;
    for (size_t ii=0;ii<x.length();ii++) {
      if (x.intensity[ii]<minRun) ++count;
      else count=0;
      if (count==3) return true;
    }
    return false;
  };

private:
  double minRun;
  int minLength;
};


/** A custom class to create a predicate function for the 454  
 * That allows us to filter flowGrams                                 */
class flowGramMaxFilter {
public:
  flowGramMaxFilter(double maxintensity):
    maxIntensity(maxintensity){};

  bool operator()(const flowGram &x) {
    if (*std::max_element(
			  x.intensity.begin()
			  ,x.intensity.end())
	>maxIntensity
	) return true;
    return false;
  }
private:
  double maxIntensity;
};

std::vector<std::vector<flowGram> > SplitbyPrimers(const std::vector<flowGram> xx, const std::vector<std::string> &primers
					       , const std::vector<double> &zp,
						   const std::vector<double> &pp, double minBayes=100);
#endif
