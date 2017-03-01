#ifndef CompressFlowSeq_H__
#define CompressFlowSeq_H__

#include <string>
#include <vector>
#include <map>
#include <cassert>
#include "flowSeq.H"
#include "newio.H"

/** The assumed order for the flowsequences                          */
const char order[4]={'G','T','A','C'};

const unsigned short mask16 =    15;
const unsigned short mask240 =   240;
const unsigned short mask3804 =  3840;
const unsigned short mask61440 = 61440;     	

template <typename I>
int val(I v , int w) {
	switch(w) {
		case 0:
		return v&mask16;
		case 1:
		return (v&mask240)<<4;
		case 2:
		return (v&mask3804) << 8;
		case 3:
		return (v&mask61440) << 12;
		default:
		std::cerr << "This position not known\n";
		return 0;
	}
}


/** A class to hold the sequences that correspond to the short integers
 * hopefully to speed up the calculation of the sequences from flowsequences  */
class CFSmap{
public:
   CFSmap() {
  }
  std::string &operator()(int v) {
    if (cmap.find(v)==cmap.end()) {  	
    	int curr=v;
      std::string s; 
      for (int i=0;i<4;i++) {
       int copies = curr&mask16 ; 
       for (int j=0;j<copies;j++) s = order[3-i]+s;
       curr = curr >> 4;
      }
      cmap[v] = s; 
    } 
    return cmap[v]; 
  }
private:
  std::map<unsigned short,std::string> cmap;
};
extern CFSmap cfsmap;
/** Can we get a compressed flowSeq so that it is easier to check for
 * alignments?  */
class CFS {
public:
/** construct from a standard flowsequence                           
 * Note that this chops off any partly constructed FlowSequence 
 * I need to make sure that this matches up with the calling routine   */
  CFS(const flowSeq &fs) {
    for (size_t ii=0;ii<fs.length()/4;ii++) { 
     ushort x = std::min(fs[4*ii+3],16) + (std::min(fs[4*ii+2],16) << 4)
     +  (std::min(fs[4*ii+1],16) << 8) + (std::min(fs[4*ii],16) << 12);
      data_.push_back(x);
    }
  }
  bool operator<(const CFS &RHS) const  {
  	return this->data_<RHS.data_;
 	};
 	/** return the number of copies at position ii for uncompressed sequence  */
	int operator[](int ii) const {
		return val(data_[ii/4],ii%4);
	}
	
	void assign(int ii, int d) { // not checked
		assert(d>0&&d<16);
		ushort x=data_[ii/4];
		switch(ii%4) {
			case(0):
			x = (x&(!mask16)) + d;
			break;
			case(1):
			x = (x&(!mask240)) + (d>>4);
			break;
			case(2):
			x=(x&!mask3804)+(d>>8);
			break;
			case(3):
			x=(x&!mask61440)+(d>>12);
			break;
		}
		data_[ii/4] = x;
	}
  /** Return the flow sequence as a string                           */
  std::string tostring() {
    std::string aa;
    for (size_t i=0;i<data_.size();i++) aa.append(cfsmap(data_[i]));
    return aa;
  }
  const std::vector<unsigned short> &operator()() const {
  	return data_;
  }

  private:
  std::vector<unsigned short> data_;
};

// need to specify a distance
int CFSdistance(int a,int b) 
{
	return abs((a%10-b%10))
	+abs(a/10-b/10)
	+abs(a/100-b/100)
	+abs(a/1000-b/1000);
}


#endif
