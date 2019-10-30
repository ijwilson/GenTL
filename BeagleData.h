#ifndef  BEAGKEDATA_H__
#define  BEAGKEDATA_H__

#include "tnt/tnt.h"
#include <ios>
#include <fstream>
#include <vector>

/** A class to deal with the Beagle program 
 *
 * It also needs to be able to deal with different types
 * of data so that we can get the positions of the SNPs and
 * the case control variables too.
 * I have added extra information so that we can add continuous 
 * random variables rather than just discrete presence/absence of 
 * trait data
 */

class CCData {
public:
  CCData(const std::string &filename);
  void addTraits(const std::string &filename);
  void addQuantitativeTraits(const std::string &filename);
  void addPositions(const std::string &filename);
  void addRegions(const std::string &filename);
  /** Get a vector with the indices of cases          */
  std::vector<int> GetCases(int val, int trait=0) const;
public:
  TNT::Array2D<int> haplotypes;
  std::vector<std::string> MarkerNames;
  TNT::Array2D<int> DiseaseTraits;
  TNT::Array2D<double> QuantitativeTraits;
  std::vector<std::string> DiseaseTraitNames;
  std::vector<std::string> QuantitativeTraitNames;
  std::vector<int> position;
  std::vector<int> region;
  size_t nsamples;
};



#endif
