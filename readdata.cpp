#include "readdata.H"
#include "newio.H"
#include "myxml.H"

#include <stdexcept>
#include <deque>

/** Note that the new sima file format has optional SNP and 
    locations that are surrounded by 
    <location> </location> 
    <SNP> .. </SNP> 
    if they are present
*/ 
bool readdata(std::istream &in,TNT::Array2D<int> &a,TNT::Array1D<double> &d
	      ,const std::string &mode)
{
  static bool first=true;
  static int n;
   if (mode=="ms"&&first==true) { // need to know how many samples
    std::string chstr;
    in >> chstr >> n;
    if (chstr != "ms") {
      std::cerr << "Expected ms input for this mode" << std::endl;
      exit(EXIT_FAILURE);
    }
   }
   if (first) first=false;
   if (mode=="sima") {   
     if (!skipto(in,"//")) return false;
     std::string rs;
     getline(in,rs);
     //string rs=in.;
     //double simrho;
     //in >> rs;assert(rs=="rho");
     //skipto(in,"=");
     //in >> simrho;
     //   std::cerr << "read: " << rs << std::endl;
     in >> a;
     in >> d;
     
    } else if (mode=="data") {
      if (!skipto(in,"//")) return false;
      std::string rs;
      getline(in,rs);
      in >> a;
      in >> d;
    } else {
      if (!skipto(in,"//")) return false;
      int segsites;
      std::string chstr;
      in >> chstr >> segsites;
      if (chstr != "segsites:") {
	std::cerr << "Expected <segsites:> in ms input for this mode, got" 
		  << chstr << std::endl;
	exit(EXIT_FAILURE);
      }
      in >> chstr;
      if (chstr!="positions:") {
	std::cerr << "Expected <positions:> in ms input for this mode, got" 
		  << chstr << std::endl;
	exit(EXIT_FAILURE);
      }
      d=TNT::Array1D<double>(segsites);
      for (int i=0;i<segsites;i++) in >> d[i];
      skipblank(in);
      a=TNT::Array2D<int>(n,segsites);
      for (int i=0;i<n;i++) {
	for (int j=0;j<segsites;j++) {
	  a[i][j]=(in.get()=='1');
	}
	if (in.get()!='\n') 
	  std::cerr << "error with input file" << std::endl;
      }
   }
   if (d.dim()!=a.dim2()) {
     std::ostringstream oss;
     oss << "Positions dimension (" << d.dim() << ") does not match number of data columns" 
	 << " (" << a.dim2() << ") in readdatalocation";
     throw std::range_error(oss.str().c_str());
   }
   return true;
}

bool readdatahaplength(std::istream &in,TNT::Array2D<int> &a,TNT::Array1D<double> &d
		       ,TNT::Array2D<int> &mn, TNT::Array2D<int> &mx, const std::string &mode)
{
  if (mode=="ms") { // need to know how many samples
     std::cerr << "need sima output for readdatahaplength" << std::endl;
     exit(EXIT_FAILURE);  
   }
  if (mode=="sima") {   
     if (!skipto(in,"//")) return false;
     std::string rs;
     double simrho;
     in >> rs;assert(rs=="rho");
     skipto(in,"=");
     in >> simrho;
     in >> a;
     in >> d;
     in >> mn;
     in >> mx;
    } else if (mode=="data") {
      if (!skipto(in,"//")) return false;
      std::string rs;
      getline(in,rs);
      in >> a;
      in >> d;
      in >> mn;
      in >> mx;
  } 
  
  if (d.dim()!=a.dim2()) {
    std::ostringstream oss;
    oss << "Positions dimension (" << d.dim() << ") does not match number of data columns" 
	<< " (" << a.dim2() << ") in readdatalocation";
    throw std::range_error(oss.str().c_str());
  }
  return true;
}

bool readdatalocation(std::istream &in,TNT::Array2D<int> &a,TNT::Array1D<double> &d, TNT::Array1D<int> &location,const std::string &mode)
{
  static bool first=true;
  static int n,npops;
  static std::deque<int> popsize;
  if (mode=="ms"&&first==true) { // need to know how many samples
    std::string chstr;
    in >> chstr >> n;
    if (chstr != "ms") 
      throw std::domain_error("Expected ms input for this mode");
    if (!skipto(in,"-I")) 
      throw std::domain_error("this is not output from a split population model");
    in >> npops;
    int siZ;
    for (int i=0;i<npops;i++) {
      in >> siZ;
      popsize.push_back(siZ);
    }
   }
   if (first) first=false;
   if (mode=="sima") {   
     if (!skipto(in,"//")) return false;
     std::string rs;
     double simrho;
     in >> rs;assert(rs=="rho");
     skipto(in,"=");
     in >> simrho;
     in >> a;
     in >> d;
     in >> location;
   } else if (mode=="data") {
     if (!skipto(in,"//")) return false;
      std::string rs;
      getline(in,rs);
      // std::cerr << "read " << rs;
      in >> a;
      // std::cerr << "read a";
      in >> d;
      in >> location;
   } else {
      if (!skipto(in,"//")) return false;
      int segsites;
      std::string chstr;
      in >> chstr >> segsites;
      if (chstr != "segsites:") {
	std::cerr << "Expected <segsites:> in ms input for this mode, got" 
		  << chstr << std::endl;
	exit(EXIT_FAILURE);
      }
      in >> chstr;
      if (chstr!="positions:") {
	std::cerr << "Expected <positions:> in ms input for this mode, got" 
		  << chstr << std::endl;
	exit(EXIT_FAILURE);
      }
      d=TNT::Array1D<double>(segsites);
      for (int i=0;i<segsites;i++) in >> d[i];
      skipblank(in);
      a=TNT::Array2D<int>(n,segsites);
      for (int i=0;i<n;i++) {
	for (int j=0;j<segsites;j++) {
	  a[i][j]=(in.get()=='1');
	}
	if (in.get()!='\n') 
	  std::cerr << "error with input file" << std::endl;
      }
      int pos=0;
      location=TNT::Array1D<int>(n);
      for (int i=0;i<npops;i++) {
	for (int j=0;j<popsize[i];j++) location[pos++]=i;
      }
   }
   if (location.dim()!=a.dim1()) {
     std::ostringstream oss;
     oss << "Locations dimension (" << location.dim() << ") does not match number of data rows" 
	 << " (" << a.dim1() << ") in readdatalocation";
     throw std::range_error(oss.str().c_str());
   }
   if (d.dim()!=a.dim2()) {
     std::ostringstream oss;
     oss << "Positions dimension (" << d.dim() << ") does not match number of data columns" 
	 << " (" << a.dim2() << ") in readdatalocation";
     throw std::range_error(oss.str().c_str());
   }

   return true;
}


// genodata::genodata(std::istream &in) 
// {
//   std::vector<int> itmp;
//   int dim[2];
  
//   //  readBetweenTags<2>(filename,"haplotypes",dim,ivec);
//   //readBetweenTags
  
  

// }
