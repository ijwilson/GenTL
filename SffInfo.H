#ifndef SFFINFO_H__
#define SFFINFO_H__

#include <string>

enum sffoption {qual,Seq,flow,accno};

std::string RunSffinfo(const char *sfffile, sffoption opt1, bool notrim=true);
std::string RunPyrobayes(const char *sfffile);


#endif
