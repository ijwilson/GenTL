#include "SffInfo.h"
#include "newio.h"

#include <iostream>
#include <fstream>
#include <sstream>

std::string RunSffinfo(const char *sfffile, sffoption opt1, bool notrim)
{
  std::ostringstream oss;
  std::string tmpfile = tempfilename("/tmp",20);
   
  oss  << "sffinfo ";
  switch (opt1) {
  case qual:
    oss << "-qual ";break;
  case Seq:
    oss << "-seq ";break;
  case flow:
    oss << "-flow ";break;
  default:
    std::cerr << "This option not supported yet in RunSffinfo";
  }

  if (notrim) oss << "-notrim ";

  oss << sfffile << " | cat >  " << tmpfile;
  std::cerr << "Calling system command\n" << oss.str().c_str() << std::endl;
  int warning = system(oss.str().c_str());
  if (warning!=0) 
    std::cerr << "possible error in RunSffinfo\n";

  return tmpfile;
}

std::string RunSffextract(const char *sfffile, sffoption opt1, bool notrim)
{
  std::ostringstream oss;
  std::string tmpfile = tempfilename("/tmp",20);
   
  oss  << "sff_extract.py ";

  if (notrim) oss << "-notrim ";

  oss << sfffile << " | cat >  " << tmpfile;
  std::cerr << "Calling system command\n" << oss.str().c_str() << std::endl;
  int warning = system(oss.str().c_str());
  if (warning!=0) 
    std::cerr << "possible error in RunSffinfo\n";

  return tmpfile;
}



std::string RunPyrobayes(const char *sfffile)
{
  std::ostringstream oss;
  std::string tmpfile = tempfilename("/tmp",20);
  oss << "PyroBayes -i " << sfffile << " -o  " << tmpfile << " > /tmp/Pyrotmpfile";
  int warning = system(oss.str().c_str());
  if (warning!=0) 
    std::cerr << "possible error in RunSffinfo\n";

  return tmpfile;

} 
