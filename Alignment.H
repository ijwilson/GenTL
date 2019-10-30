#ifndef ALIGNMENT_H__
#define ALIGNMENT_H__

#include <map>
#include <algorithm>
#include "fasta.H"
#include "flowSeq.H"

class GaplessAlignment {
	public:
	GaplessAlignment(const std::vector<fasta<char> > &refseqs);
	std::pair<int,int> operator()(const flowSeq &f);
	
	const flowSeq &operator()(int ii) const {
		return refs[ii];
	}
	const vector<int> refCopies(size_t num, size_t start, size_t length) {
      return refs[num].subcopies(start,length);
	}
private:
	std::vector<flowSeq> refs;	
};

#endif
