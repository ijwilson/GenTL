/**   @file  */
#ifndef READDATA_H__
#define READDATA_H__

#include <iostream>
#include <string>
#include "tnt/tnt.h"
#include <vector>

using std::vector;


/** A structure to keep track of data - 
 * should mirror the structures seen within R to read and write 
 * sima xml files                                                 
 */

struct genodata {
public:
  genodata(std::istream &in);
  
  vector<int> location;
  vector<int> haplogroup;
  vector<int> pos;
  
  TNT::Array2D<int> d;
};



bool readdata(std::istream &in,TNT::Array2D<int> &a,TNT::Array1D<double> &d,const std::string &mode="data");
bool readdatalocation(std::istream &in,TNT::Array2D<int> &a,TNT::Array1D<double> &d, TNT::Array1D<int> &location
		      ,const std::string &mode="data");

bool readdatahaplength(std::istream &in,TNT::Array2D<int> &a,TNT::Array1D<double> &d
		       ,TNT::Array2D<int> &mn, TNT::Array2D<int> &mx, const std::string &mode);





#endif
