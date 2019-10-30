
#ifndef FLOWSEQ_H__
#define FLOWSEQ_H__
#include <string>
#include <algorithm>
#include <vector>
#include <utility>  // for pair

using std::ostream;
using std::string;
using std::vector;

/** A class to represent a flow-wise version of a sequence   */
class flowSeq {
public:
  /** constructors for flowseq                                      */
  flowSeq(const std::string &seq,
	  bool TrimFloatingStart, 
	  bool TrimFloatingEnd,const string &name="", 
	  const string &Floworder="TACG");
  
  flowSeq(const std::string &seq,const string &name="",const string &Floworder="TACG");
  flowSeq(size_t len,const string &Floworder="TACG");
  flowSeq(){};
  /** write the flowsequence to a stream                                */
  ostream &print(ostream &o,const string &sep=",") const;
  /** convert a flowseq back into a sequence    */
  string sequence() const;
  /** What is the probability of a flowgram given this sequence 
   * the flowgram starting at the correct point              
   * assume means of 0,1,2,... and variances of sigma_i.  Note that 
   * I assume a lognormal distribution for the zeros and normals for the ones 
   **/
  double lprob(const vector<double> &flowgram, 
	       const vector<double> &zp, const vector<double> &pp
	       ,const std::vector<int> &tag) const;
  double lprobsimple(const std::vector<double> &flowgram, const std::vector<double> &zp
		    , const std::vector<double> &pp) const;
    int operator[](size_t ii) const {
    	return copies[ii];
  	}
	/** helper function to return the number of copies at position ii */
	int &operator[](size_t ii)  {
    	return copies[ii];
  	}
	/** What is the length of a sequence                              */
	size_t length() const {
    	return copies.size();
  	}
  	const std::vector<int> &operator()() const {
  		return copies;
	}
	/** an overloaded equality operator                               */
  	bool operator==(const flowSeq &rhs) const {
  		return copies==rhs.copies;	
  	}
  /** extract subsequence copies   */
  std::vector<int> subcopies(size_t start,size_t len) {
    size_t mx=std::min(start+len,length());
    return std::vector<int>(copies.begin()+start,copies.begin()+mx);
  }
  /** Is flowSeq a a subset of this flowSeq                           
   * If it is then this returns the first and last positions, otherwise return first<0  
   * This should allow us to categorise sequences by which ones partially match        */
  std::pair<int,int> subflowgram(const flowSeq &a) const {
    int start=0;
    size_t alen=a.copies.size();
    size_t left=copies.size();
    if (alen>left) return std::pair<int,int>(-1,-1);
    std::vector<int>::const_iterator ii=copies.begin();
    std::vector<int>::const_iterator jj=a.copies.begin();

    for (;;) {
      if (jj==a.copies.end()) {
        return std::pair<int,int>(start,std::distance(copies.begin(),ii));  // reached the end so a match
      } else if (*ii != *jj) {                // a mismatch here - start again
        jj=a.copies.begin();
        alen=a.copies.size();
        left--;ii++;
        if (alen>left) return std::pair<int,int>(-1,-1);  // can we fit it in
        start=std::distance(copies.begin(),ii);
      } else {
        alen--;left--;ii++;jj++;
      }
    }
  }

  void checkseq(string message,int maxCopies) const;
private:
  std::vector<int> copies;
  string order;
  string name;
};

std::ostream &operator<<(ostream &out, const flowSeq &a);

#endif
