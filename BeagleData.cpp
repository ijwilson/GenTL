#include "BeagleData.H"
#include <iostream>
#include <iterator>
#include <stdexcept>

typedef std::istream_iterator<int> ifstream_iit;
typedef std::istream_iterator<double> ifstream_dit;

CCData::CCData(const std::string &filename) {
  std::vector<std::vector<int> > xx;   // SNPS
  std::ifstream in(filename.c_str());
  if (!in) {
    std::cerr << "Error, cannot open input file " << filename << " in CCData exiting\n";
    exit(EXIT_FAILURE);
  }
  // first ignore any starting lines that begin with #
  std::string line;
  for (;;) {
    getline(in,line);
    if (line[0]!='#') break;
  }
  // read the first line
  for (;;) {
    std::istringstream iss(line);
    char c;
    std::string str;
    iss >> c >> str;
    if (c=='M') { // a marker
      MarkerNames.push_back(str);
      std::vector<int> a;
      std::copy(ifstream_iit(iss),ifstream_iit(),std::back_inserter(a));
      xx.push_back(a);
    }
    getline(in,line);
    if (in.fail()) break;
  }
  if (xx.size()==0) {  
    std::ostringstream oss;
    oss << "have not read any marker data from file " << filename << "\nexiting ....\n";
    throw std::runtime_error(oss.str());
  }
    
  nsamples=xx[0].size();
  haplotypes=TNT::Array2D<int>(nsamples,xx.size());
  for (size_t i=0;i<xx.size();i++) {
    if (xx[i].size()!=nsamples) {
      std::ostringstream oss;
      oss << "incorrect number of columns for row "<< i 
	  << "\nHave " << xx[i].size() << " and should have " << nsamples
	  << "\nin  CCdata constructor\n" << std::endl;
      throw std::runtime_error(oss.str());
    }
    for (size_t j=0;j<nsamples;j++) haplotypes[j][i]=xx[i][j];
  }
}

void CCData::addQuantitativeTraits(const std::string &filename) {
  std::vector<std::vector<double> > yy;   // Traits
  std::ifstream in(filename.c_str());
  if (!in) {
    std::cerr << "Error, cannot open qtl file " << filename << " exiting\n";
    exit(EXIT_FAILURE);
  }
  // first ignore any starting lines that begin with #
  std::string line;
  for (;;) {
    getline(in,line);
    if (line[0]!='#') break;
    //   std::cerr << "skipping line\n";
  }
  // read the first line
  for (;;) {
    std::istringstream iss(line);
    char c;
    std::string str;
    iss >> c >> str;
    if (in.fail()) break;
    if (c=='Q') { // a marker
      QuantitativeTraitNames.push_back(str);
      std::vector<double> a;
      std::copy(ifstream_dit(iss),ifstream_dit(),std::back_inserter(a));
      yy.push_back(a);
    }
    getline(in,line);
    if (line=="") break;
  }
  if (yy.size()>0) {
    if (haplotypes.dim1()==0) nsamples=yy[0].size(); 
    QuantitativeTraits=TNT::Array2D<double>(yy.size(),nsamples);
    for (size_t i=0;i<yy.size();i++) {
      if (yy[i].size()!=nsamples) {
        std::ostringstream oss;
        oss << "incorrect number of columns for row "<< i  << " in addQuantitativeTraits"
	      << "\nhave " << yy[i].size() << " and should have " << nsamples
	      << " in AddTraits\n" << std::endl;
	  throw std::runtime_error(oss.str());
      }
      for (size_t j=0;j<nsamples;j++) QuantitativeTraits[i][j]=yy[i][j];
    }
  }
}
void CCData::addTraits(const std::string &filename) {
  std::vector<std::vector<int> > yy;   // Traits
  std::ifstream in(filename.c_str());
  if (!in) {
    std::cerr << "Error, cannot open input file " << filename << " exiting\n";
    exit(EXIT_FAILURE);
  }
  // first ignore any starting lines that begin with #
  std::string line;
  for (;;) {
    getline(in,line);
    if (line[0]!='#') break;
    //   std::cerr << "skipping line\n";
  }
  // read the first line
  for (;;) {
    std::istringstream iss(line);
    char c;
    std::string str;
    iss >> c >> str;
    if (in.fail()) break;
    if (c=='A') { // a marker
      DiseaseTraitNames.push_back(str);
      std::vector<int> a;
      std::copy(ifstream_iit(iss),ifstream_iit(),std::back_inserter(a));
      yy.push_back(a);
    }
    getline(in,line);
  }
  if (yy.size()>0) {
    if (haplotypes.dim1()==0) nsamples=yy[0].size(); 
    DiseaseTraits=TNT::Array2D<int>(yy.size(),nsamples);
    for (size_t i=0;i<yy.size();i++) {
      if (yy[i].size()!=nsamples) {
	  std::ostringstream oss;
	  oss << "incorrect number of columns for row "<< i 
	      << "\nhave " << yy[i].size() << " and should have " << nsamples
	      << " in AddTraits\n" << std::endl;
	  throw std::runtime_error(oss.str());
      }
      for (size_t j=0;j<nsamples;j++) DiseaseTraits[i][j]=yy[i][j];
    }
  }
}

void CCData::addPositions(const std::string &filename) {
  position.resize(MarkerNames.size());
  std::ifstream in(filename.c_str());
  if (!in) {
    std::cerr << "Error, cannot open position file " << filename << " exiting\n";
    exit(EXIT_FAILURE);
  }
  // first ignore any starting lines that begin with #
  std::string line;
  for (;;) {
    getline(in,line);
    if (line[0]!='#') break;
    //   std::cerr << "skipping line\n";
  }
  int currpos=0;
  for (;;) {
    std::istringstream iss(line);
    //    std::cerr << line << std::endl;
    std::string name;
    char al1,al2;
    int pos;
    iss >> name >> pos >> al1 >> al2;
    
    if (in.fail()) break;
    // have to account for dropped markers
    if (name == MarkerNames[currpos]) position[currpos++]=pos;
    if (currpos==static_cast<int>(position.size())) break;
    getline(in,line);
  }

  if (currpos!=static_cast<int>(position.size())) {
    std::ostringstream oss;
    oss << "Error, not enough values have been read from positions\n";
    oss << "Needed " << position.size() << " and got " << currpos << std::endl;
    throw std::runtime_error(oss.str());
  }
}
/** Extract the cases and return as a vector of labels             */
std::vector<int> CCData::GetCases(int val,int trait) const
{
  std::vector<int> aa;
  if (trait>DiseaseTraits.dim1()) {
    std::cerr << "error, this trait not available\n";
  }
  for (int i=0;i<DiseaseTraits.dim2();i++) 
    if (DiseaseTraits[trait][i]==val) aa.push_back(i);
  
  return aa;
}



void CCData::addRegions(const std::string &filename) {
  region.resize(0);
  std::ifstream in(filename.c_str());
  if (!in) {
    std::cerr << "Error, cannot open region file " << filename << " exiting\n";
    exit(EXIT_FAILURE);
  }
  // get the comment in the first line
  std::string line;
  getline(in,line);
  //std::cerr << line;
 
  for (;;) {
    int pos;
    in >> pos;
    if (in.fail()) break;
    region.push_back(pos);
 
  }
  //std::cerr << region.size() << " " << nsamples << std::endl;
  if (region.size()!=nsamples)  {
    std::cerr << "Error, length incorrect in regions file, have length " << region.size() << " and require " << nsamples << " exiting \n";
    exit(EXIT_FAILURE);
  }
  // need some check
}
