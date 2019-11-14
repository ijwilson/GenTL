#ifndef PYROSEQ_H__
#define PYROSEQ_H__
#include <string>
#include "fasta.h"
#include "SffInfo.h"
/** A class to call the Pyrobayes software and read in the results from an sff file */
class PyroSeq {
public:
  PyroSeq(const char *sfffile)  {
    std::string filename =  RunPyrobayes(sfffile);
    std::string fastafile = filename+".fasta";
    std::string qualfile=filename+".fasta.qual";
    std::ifstream inf(fastafile.c_str());
    std::ifstream inq(qualfile.c_str());
    while (!inf.eof()) seq_.push_back(fasta<char>(inf));
    while (!inq.eof()) qual_.push_back(fasta<int>(inq));
    std::cerr << "read " << seq_.size() << " sequences " << " and " 
	      << qual_.size() << " quality scores\n";
    makeMap();
  }
  
  std::pair<fasta<char>,fasta<int> > seq(const std::string &name)  {
    if (seqmap.find(name)==seqmap.end()) {
      std::cerr << "unable to find that sequence - expect a crash\n\n";
      exit(EXIT_FAILURE);
    }
    int num = seqmap[name];
     return std::pair<fasta<char>, fasta<int> >(seq_[num],qual_[num]); 
  }

  std::pair<fasta<char>,fasta<int> > seq(int ii)  {
     return std::pair<fasta<char>, fasta<int> >(seq_[ii],qual_[ii]); 
  }

  std::string sequence(int ii) {
    return seq_[ii].tostring();
  }
  const std::vector<fasta<char> > &sequences() const {
    return seq_;
  }

  void makeMap() {
    for (size_t ii=0;ii<seq_.size();ii++) {
      if (seq_[ii].name()!=qual_[ii].name() ) {
	std::cerr << "problem, names in sequences and quality scores do not match at position " << ii 
		  << "\nhave " << seq_[ii].name() << " and "<< qual_[ii].name() << std::endl;
	exit(EXIT_FAILURE);
      }
      seqmap[seq_[ii].name()]=ii;
    }
  }
private:
  std::map<std::string,int> seqmap;
  std::vector<fasta<char> > seq_;  // the Roche 454 sequences
  std::vector<fasta<int> > qual_;  // the Roche 454 quality scores
};
#endif
