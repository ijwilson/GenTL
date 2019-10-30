#ifndef SEQUENCE_COLLLECTION_H__
#define SEQUENCE_COLLLECTION_H__

#include "flowSeq.H"
#include <map>
#include <vector>

class FlowSeqCollection{
	typedef std::map<CFS, std::vector<int> > fsMap;
	typedef std::map<CFS, std::vector<int> >::const_iterator fsMap_citor;
public:
	FlowSeqCollection(const std::vector<CFS> &fs) {
		for (size_t ii=0;ii<fs.size();ii++) 
			aa[fs[ii]].push_back(ii);
	};
	
	  
	size_t nseq() const {
			return aa.size();
	};
	
	std::pair<size_t,CFS > maxfreq() const {
		fsMap_citor i=aa.begin();
		fsMap_citor mxi;
		size_t mx=1;
		while (i!=aa.end()) {
			if (i->second.size()>mx) {
				mx =  i->second.size();
				mxi=i;
			}
			i++;	
		}
		return std::pair<size_t,CFS >(mx,mxi->first);
	};
	
	std::pair<size_t,std::vector<int> > whmax() const {
		fsMap_citor i=aa.begin();
		fsMap_citor mxi;
		size_t mx=1;
		while (i!=aa.end()) {
			if (i->second.size()>mx) {
				mx =  i->second.size();
				mxi=i;
			}
			i++;	
		}
		return std::pair<size_t,std::vector<int> >(mx,mxi->second);
	};
	
	std::vector<size_t> counts() const {
		std::vector<size_t> xx;
		xx.reserve(aa.size());
		fsMap_citor i=aa.begin();
		while (i!=aa.end()) {
			xx.push_back(i->second.size());
			i++;
		}
		std::sort(xx.begin(),xx.end());
		return xx;
	}

private:
	 fsMap aa; 
};


#endif
