#ifndef TRINODE_H__
#define TRINODE_H__

#include <vector>
#include <iterator>
#include <iosfwd>
#include <stack>
#include "tnt/tnt.h"

#include "testStats.h"


/** Note that both of these should be sorted                               */
template<typename T>
int CountIntersection(const std::vector<T> &a, const std::vector<T> &b)
{
  typename std::vector<T>::const_iterator al=a.end(),bl=b.end(),ii=a.begin(),jj=b.begin();
  int count=0;
  while (ii!=al && jj!=bl) {
    if (*ii<*jj) ++ii;
    else if (*jj<*ii) ++jj;
    else {
      count++;
      ii++;
      jj++;
    }
  }
  return count;
} 

using std::pair;
using std::stack;
/** Node, left centre right iterator    */
template <class T> class NLCRIterator;
/** A class for representing haplotype trees                               
 *  The trees are represented as binary nodes with the addition of a 
 * vector of labels for each node, and the position of the SNP that causes 
 * the particular split.  Up is included, but not used yet, so that it
 * is easier to calculate distances between nodes.
 *
 */
class trinode {
public:
  trnode(const std::vector<int> &a,trinode *father=0):
    left(0),centre(0),right(0),up(father),position(-1),labels(a) {}
  /** Print the labels at a node                                           */
  std::ostream & printlabels(std::ostream &o) {
    std::copy(labels.begin(),labels.end(),std::ostream_iterator<int>(o," "));
    return o;
  }  
  
  /** A function to turn a tree into something that looks like an ape tree object    
   * count is the number that should be assigned to this node                    
   * the function returns the number of the next available node            */
  int recurseApeSplit(std::vector<std::pair<int,int> > &Edge, std::vector<std::vector<int> > &labs
		      ,int count)
  {
    assert(!isleaf());
    int newcount=count+1;
    if (left->isleaf()) {
      labs.push_back(left->labels);
      Edge.push_back(pair<int,int>(count,labs.size()));
    } else {
      Edge.push_back(std::pair<int,int>(count,newcount));
      newcount=left->recurseApeSplit(Edge,labs,newcount);
    } 
    if (right->isleaf()) {
      labs.push_back(right->labels);
      Edge.push_back(std::pair<int,int>(count,labs.size()));
    } else {
      Edge.push_back(std::pair<int,int>(count,newcount));
      newcount=right->recurseApeSplit(Edge,labs,newcount);
    }
    return newcount;
  }
  /** is the node a leaf?    by definition a node mustr have at least two descendents   */
  bool isleaf() const {
    return (descendents() < 2);
  }
  /** how many descendents does a node have   */
  bool descendents() const {
    return (left==0)+(right==0)+(centre==0);
  }
   std::pair<int,int> GetLength() {
     if (labels.size()==1) return std::pair<int,int>(-1,-1);
   }

  /** returns a string which gives the set of mutations 
   * that occur to get here on the tree    
   */

  std::vector<std::pair<int,int> > ToGetHere() const {
    if (up==0) return std::vector<std::pair<int,int> >(0);
    
    std::vector<std::pair<int,int> > res = up->ToGetHere();
    if (up->left==this) res.push_back(std::pair<int,int>(up->position,0));
    else  res.push_back(std::pair<int,int>(up->position,1));
    if (isleaf())  res.push_back(std::pair<int,int>(0,0));
    return res;
  }
   
   int GetLength(int centre,const TNT::Array2D<int> &haps) 
   {
   		if (labels.size()==1) return -1;
   		int nSNP=haps.dim2();
        int Left,Right;
        size_t jj;
        // find the maximum value to the left that is shared
        for (Left=centre;Left>=0;Left--) {
          int MatchSNP=haps[labels[0]][Left];
          for (jj=1;jj<labels.size();jj++) {
            if (haps[labels[jj]][Left]!=MatchSNP) break;
          }
          // we we have reached the end of labels and there was no split
          // then continue - otherwise break
          if (jj!=labels.size()) break;
        }
        // Left should hold the maximum split to the left
        for (  Right=centre;Right<nSNP;Right++) {
          int MatchSNP=haps[labels[0]][Right];
          for (jj=1;jj<labels.size();jj++) {
            if (haps[labels[jj]][Right]!=MatchSNP) break;
          }
          if (jj!=labels.size()) break;
        }
        // Right should hold the maximum split to the right
        // or nSNP if it has reached the end
        return Right-Left;	
   }

 
  std::vector<double> RecurseTestStatistics(const std::vector<int> &Cases, const double p,const int maxk, int CC[2], double (*CalcStat)(int,int, double)) {

    std::vector<double> newMaxValues(maxk,0.0);
    
    if (isleaf()) {
      CC[0]=CountIntersection(Cases,labels);
      CC[1]=static_cast<int>(labels.size());
      if (CC[1]>1) newMaxValues[0]=CalcStat(CC[0],CC[1],p);
    } else {
      int LeftCC[2],RightCC[2];
      std::vector<double> LeftStat=left->RecurseTestStatistics(Cases,p,maxk,LeftCC,CalcStat);
      std::vector<double> RightStat=right->RecurseTestStatistics(Cases,p,maxk,RightCC,CalcStat);
      CC[0]=LeftCC[0]+RightCC[0];
      CC[1]=LeftCC[1]+RightCC[1];
      // First calculate the maximum node below this one (including this one
      // maximum of this node and the maxima to the below left and right
      newMaxValues[0]=std::max(CalcStat(CC[0],CC[1],p),std::max(LeftStat[0],RightStat[0]));
      // now work out the maximum disjoint values for k = 2 up to its max value
      for (int d=1;d<maxk;d++) {  // This the the maximum CountStat values
        double maxVal=std::max(LeftStat[d],RightStat[d]);    
        int sumIndices=d-1;
        // what is the maximum made up of 0 to the left and sumIndicies to the right and the opposite
        for (int i=0;i<=sumIndices;i++) // i+1 to the left and CountStat-i-1 to the right
          if (LeftStat[i]+RightStat[sumIndices-i]>maxVal) maxVal = LeftStat[i]+RightStat[sumIndices-i];
        newMaxValues[d] = maxVal;
      }
    }
    return newMaxValues;  
  }

  std::vector<double> RecurseZStatA(const std::vector<double> &xx, double &sumx,const int maxk,double meanvar[2]) {

    std::vector<double> newMaxValues(maxk,0.0);
    
    if (isleaf()) {
      /** need to calculate test statistics for a leaf node - get the mean  */
      sumx=0.0;
      for (size_t ii=0;ii<labels.size();ii++) sumx += xx[labels[ii]];
  
      double n=static_cast<double>(labels.size());
      newMaxValues[0]= fabs((sumx/n - meanvar[0])/(sqrt(meanvar[1]/n)));
    } else {
      double Leftsumx,Rightsumx;
      std::vector<double> LeftStat=left->RecurseZStatA(xx,Leftsumx,maxk,meanvar);
      std::vector<double> RightStat=right->RecurseZStatA(xx,Rightsumx,maxk,meanvar);
      sumx=Leftsumx+Rightsumx;
      double n=static_cast<double>(labels.size());
      // First calculate the maximum node below this one (including this one
      // maximum of this node and the maxima to the below left and right
      newMaxValues[0]=std::max(fabs((sumx/n - meanvar[0])/(sqrt(meanvar[1]/n))),std::max(LeftStat[0],RightStat[0]));
      // now work out the maximum disjoint values for k = 2 up to its max value
      for (int d=1;d<maxk;d++) {  // This the the maximum CountStat values
        double maxVal=std::max(LeftStat[d],RightStat[d]);    
        int sumIndices=d-1;
        // what is the maximum made up of 0 to the left and sumIndicies to the right and the opposite
        for (int i=0;i<=sumIndices;i++) // i+1 to the left and CountStat-i-1 to the right
          if (LeftStat[i]+RightStat[sumIndices-i]>maxVal) maxVal = LeftStat[i]+RightStat[sumIndices-i];
        newMaxValues[d] = maxVal;
      }
    }
    return newMaxValues;  
  }
  std::vector<double> RecurseZStatP(const std::vector<double> &xx, double &sumx,const int maxk,double meanvar[2]) {
    std::vector<double> newMaxValues(maxk,0.0);
    if (isleaf()) {
      /** need to calculate test statistics for a leaf node - get the mean  */
      sumx=0.0;
      for (size_t ii=0;ii<labels.size();ii++) sumx += xx[labels[ii]];
  
      double n=static_cast<double>(labels.size());
      newMaxValues[0]= (sumx/n - meanvar[0])/(sqrt(meanvar[1]/n));
    } else {
      double Leftsumx,Rightsumx;
      std::vector<double> LeftStat=left->RecurseZStatP(xx,Leftsumx,maxk,meanvar);
      std::vector<double> RightStat=right->RecurseZStatP(xx,Rightsumx,maxk,meanvar);
      sumx=Leftsumx+Rightsumx;
      double n=static_cast<double>(labels.size());
      // First calculate the maximum node below this one (including this one
      // maximum of this node and the maxima to the below left and right
      newMaxValues[0]=std::max((sumx/n - meanvar[0])/(sqrt(meanvar[1]/n)),std::max(LeftStat[0],RightStat[0]));
      // now work out the maximum disjoint values for k = 2 up to its max value
      for (int d=1;d<maxk;d++) {  // This the the maximum CountStat values
        double maxVal=std::max(LeftStat[d],RightStat[d]);    
        int sumIndices=d-1;
        // what is the maximum made up of 0 to the left and sumIndicies to the right and the opposite
        for (int i=0;i<=sumIndices;i++) // i+1 to the left and CountStat-i-1 to the right
          if (LeftStat[i]+RightStat[sumIndices-i]>maxVal) maxVal = LeftStat[i]+RightStat[sumIndices-i];
        newMaxValues[d] = maxVal;
      }
    }
    return newMaxValues;  
  }

  std::vector<double> RecurseZStatN(const std::vector<double> &xx, double &sumx,const int maxk,double meanvar[2]) {

    std::vector<double> newMaxValues(maxk,0.0);
    
    if (isleaf()) {
      /** need to calculate test statistics for a leaf node - get the mean  */
      sumx=0.0;
      for (size_t ii=0;ii<labels.size();ii++) sumx += xx[labels[ii]];
  
      double n=static_cast<double>(labels.size());
      newMaxValues[0]= -(sumx/n - meanvar[0])/(sqrt(meanvar[1]/n));
    } else {
      double Leftsumx,Rightsumx;
      std::vector<double> LeftStat=left->RecurseZStatN(xx,Leftsumx,maxk,meanvar);
      std::vector<double> RightStat=right->RecurseZStatN(xx,Rightsumx,maxk,meanvar);
      sumx=Leftsumx+Rightsumx;
      double n=static_cast<double>(labels.size());
      // First calculate the maximum node below this one (including this one
      // maximum of this node and the maxima to the below left and right
      newMaxValues[0]=std::max(-(sumx/n - meanvar[0])/(sqrt(meanvar[1]/n)),std::max(LeftStat[0],RightStat[0]));
      // now work out the maximum disjoint values for k = 2 up to its max value
      for (int d=1;d<maxk;d++) {  // This the the maximum CountStat values
        double maxVal=std::max(LeftStat[d],RightStat[d]);    
        int sumIndices=d-1;
        // what is the maximum made up of 0 to the left and sumIndicies to the right and the opposite
        for (int i=0;i<=sumIndices;i++) // i+1 to the left and CountStat-i-1 to the right
          if (LeftStat[i]+RightStat[sumIndices-i]>maxVal) maxVal = LeftStat[i]+RightStat[sumIndices-i];
        newMaxValues[d] = maxVal;
      }
    }
    return newMaxValues;  
  }

  std::vector<double> RecurseZStat2(const std::vector<double> &xx, double &sumx,const int maxk,double meanvar[2]) {

    std::vector<double> newMaxValues(maxk,0.0);
    
    if (isleaf()) {
      /** need to calculate test statistics for a leaf node - get the mean  */
      sumx=0.0;
      for (size_t ii=0;ii<labels.size();ii++) sumx += xx[labels[ii]];
  
      double n=static_cast<double>(labels.size());
      newMaxValues[0]= n*(sumx/n - meanvar[0])*(sumx/n - meanvar[0])/meanvar[1];
    } else {
      double Leftsumx,Rightsumx;
      std::vector<double> LeftStat=left->RecurseZStat2(xx,Leftsumx,maxk,meanvar);
      std::vector<double> RightStat=right->RecurseZStat2(xx,Rightsumx,maxk,meanvar);
      sumx=Leftsumx+Rightsumx;
      double n=static_cast<double>(labels.size());
      // First calculate the maximum node below this one (including this one
      // maximum of this node and the maxima to the below left and right
      newMaxValues[0]=std::max(n*(sumx/n - meanvar[0])*(sumx/n - meanvar[0])/meanvar[1],std::max(LeftStat[0],RightStat[0]));
      // now work out the maximum disjoint values for k = 2 up to its max value
      for (int d=1;d<maxk;d++) {  // This the the maximum CountStat values
        double maxVal=std::max(LeftStat[d],RightStat[d]);    
        int sumIndices=d-1;
        // what is the maximum made up of 0 to the left and sumIndicies to the right and the opposite
        for (int i=0;i<=sumIndices;i++) // i+1 to the left and CountStat-i-1 to the right
          if (LeftStat[i]+RightStat[sumIndices-i]>maxVal) maxVal = LeftStat[i]+RightStat[sumIndices-i];
        newMaxValues[d] = maxVal;
      }
    }
    return newMaxValues;  
  }
  std::vector<double> RecurseContinuousTestStatistics(const std::vector<double> &xx,const double meanvar[2],const int maxk,double sms[2], double (*CalcStat)(double *,const double *, int)) {

    std::vector<double> newMaxValues(maxk,0.0);
    
    if (isleaf()) {
      /** need to calculate test statistics for a leaf node - get the mean  */
      sms[0]=sms[1]=0.0;
      for (size_t ii=0;ii<labels.size();ii++) {
        sms[0] += xx[labels[ii]];
        sms[1] += xx[labels[ii]]*xx[labels[ii]];
      }
      int n=static_cast<int>(labels.size());
      newMaxValues[0]=CalcStat(sms,meanvar,n);
    } else {
      double Leftsms[2],Rightsms[2];
      std::vector<double> LeftStat=left->RecurseContinuousTestStatistics(xx,meanvar,maxk,Leftsms,CalcStat);
      std::vector<double> RightStat=right->RecurseContinuousTestStatistics(xx,meanvar,maxk,Rightsms,CalcStat);
      sms[0] = Leftsms[0]+Rightsms[0];
      sms[1] = Leftsms[1]+Rightsms[1];
      int n=static_cast<int>(labels.size());
      // First calculate the maximum node below this one (including this one
      // maximum of this node and the maxima to the below left and right
      newMaxValues[0]=std::max(CalcStat(sms,meanvar,n),std::max(LeftStat[0],RightStat[0]));
      // now work out the maximum disjoint values for k = 2 up to its max value
      for (int d=1;d<maxk;d++) {  // This the the maximum CountStat values
        double maxVal=std::max(LeftStat[d],RightStat[d]);    
        int sumIndices=d-1;
        // what is the maximum made up of 0 to the left and sumIndicies to the right and the opposite
        for (int i=0;i<=sumIndices;i++) // i+1 to the left and CountStat-i-1 to the right
          if (LeftStat[i]+RightStat[sumIndices-i]>maxVal) maxVal = LeftStat[i]+RightStat[sumIndices-i];
        newMaxValues[d] = maxVal;
      }
    }
    return newMaxValues;  
  }


  /** The sum of minimum distances to the nearest case, nearest control and nearest different 
   *   The returned value consists of D0 - sum of cases, D1, sum of controls, d2 - sum of differences
   * d3 - distance to nearest case to left, d4 distance to nearest control to left
   * d5 - distance to nearest case to right, d6 - distance to nearest control to right         
   */
  std::vector<double> TreeDistanceStatistic(const std::vector<int> &Cases, int CC[2], int d[2]);

  std::vector<std::pair<double, std::vector< binode *> > >  RecurseTestNodes(const std::vector<int> &Cases, const double p,const int maxk, int CC[2], double (*CalcStat)(int,int, double)) {
    
    std::vector<std::pair<double, std::vector< binode *> > > newMaxValues(maxk);
    for (int i=0;i<maxk;i++) {
      newMaxValues[i].first=0.0;
      newMaxValues[i].second.resize(i+1,0);
    }
    
    if (isleaf()) {
      CC[0]=CountIntersection(Cases,labels);
      CC[1]=static_cast<int>(labels.size());
      if (CC[1]>1) newMaxValues[0].first=CalcStat(CC[0],CC[1],p);
      newMaxValues[0].second[0]=this;
    } else {
      int LeftCC[2],RightCC[2];
      std::vector<std::pair<double, std::vector< binode *> > > LeftStat=left->RecurseTestNodes(Cases,p,maxk,LeftCC,CalcStat);
      std::vector<std::pair<double, std::vector< binode *> > >  RightStat=right->RecurseTestNodes(Cases,p,maxk,RightCC,CalcStat);
      CC[0]=LeftCC[0]+RightCC[0];
      CC[1]=LeftCC[1]+RightCC[1];
      // maximum of this node and the maxima to the below left and right
      if (LeftStat[0].first>RightStat[0].first) {
        newMaxValues[0].first = LeftStat[0].first;
        newMaxValues[0].second[0] =  LeftStat[0].second[0];
      } else {
        newMaxValues[0].first = RightStat[0].first;
        newMaxValues[0].second[0] =  RightStat[0].second[0];
      }
      if (CalcStat(CC[0],CC[1],p) > newMaxValues[0].first) {
        newMaxValues[0].first = CalcStat(CC[0],CC[1],p);
        newMaxValues[0].second[0] = this;
      } 
      
      // now work out the maximum disjoint values for k up to its max value
      for (int d=1;d<maxk;d++) {  // This the the maximum CountStat+1 values
        // what is the maximum made up of 1 to the left and CountStat to the right
        int UseVal;
        double maxVal;
        if (LeftStat[d].first>RightStat[d].first) {
           maxVal = LeftStat[d].first;
           UseVal=-1;
        } else {
          maxVal = RightStat[d].first;
          UseVal=d;
        }
        int sumIndices=d-1;
        for (int i=0;i<=sumIndices;i++) {// i+1 to the left and CountStat-i-1 to the right
          if (LeftStat[i].first+RightStat[sumIndices-i].first>maxVal) {
            maxVal = LeftStat[i].first+RightStat[sumIndices-i].first;
            UseVal = i;
          }
        }
        newMaxValues[d].first = maxVal;
        if (UseVal<0) newMaxValues[d].second = LeftStat[d].second;
        else if (UseVal==d) newMaxValues[d].second = RightStat[d].second;
        else {
          std::vector<binode *>::iterator next = std::copy(LeftStat[UseVal].second.begin(),LeftStat[UseVal].second.end(),newMaxValues[d].second.begin());
          std::copy(RightStat[sumIndices-UseVal].second.begin(),RightStat[sumIndices-UseVal].second.end(),next);
        }
      }
    }
    return newMaxValues;  
  }



  std::vector<double> SevonTestStatisticsFaster(const std::vector<int> &Cases, const double p,const int maxk, int CC[2]) {

    std::vector<double> newMaxValues(maxk,0.0);
    
    if (isleaf()) {
      CC[0]=CountIntersection(Cases,labels);
      CC[1]=static_cast<int>(labels.size());
      newMaxValues[0]=(static_cast<double>(CC[0])-CC[1]*p)/sqrt(CC[1]*p*(1.-p));    
    } else {
      int LeftCC[2],RightCC[2];
      std::vector<double> LeftStat=left->SevonTestStatisticsFaster(Cases,p,maxk,LeftCC);
      std::vector<double> RightStat=right->SevonTestStatisticsFaster(Cases,p,maxk,RightCC);
      CC[0]=LeftCC[0]+RightCC[0];
      CC[1]=LeftCC[1]+RightCC[1];
      double z=(static_cast<double>(CC[0])-CC[1]*p)/sqrt(CC[1]*p*(1.-p));    
      // first get the maximum value from the left, the right and the current value
      // Note that this is the only place the current value can be used
      newMaxValues[0]=std::max(LeftStat[0],RightStat[0]);
      newMaxValues[0]=std::max(newMaxValues[0],z);
      // now work out the maximum disjoint values for k up to its max value
      for (int CountStat=1;CountStat<maxk;CountStat++) {  // This the the maximum CountStat+1 values
        // what is the maximum made up of 1 to the left and CountStat to the right
        double maxVal=LeftStat[0]+RightStat[CountStat-1];
        for (int i=1;i<CountStat;i++) // i+1 to the left and CountStat-i-1 to the right
          if (LeftStat[i]+RightStat[CountStat-i-1]>maxVal) maxVal = LeftStat[i]+RightStat[CountStat-i-1];
        newMaxValues[CountStat] = maxVal;
        // not sure if this is right but allow it to be the same as the one below if needed
        if (newMaxValues[CountStat]<newMaxValues[CountStat-1]) newMaxValues[CountStat]=newMaxValues[CountStat-1];
      }
    }
    return newMaxValues;  
  }



  std::vector<std::pair<int,double> > 
  maxLengths(int centre, int k, int maxn, const TNT::Array2D<int> &data) {
    std::vector<std::pair<int,double> > Left;
    if (labels.size()<static_cast<size_t>(k)) return Left;
    if (isleaf()) {
      Left.push_back(std::pair<int,double>(labels.size(),GetLength(centre,data)));
      return Left;
    } 
    std::pair<int,double> here(labels.size(),GetLength(centre,data));
    //Left.second=left->maxLengths(centre,k,maxn,data);
    // Left.first=left->labels.size();
    //int Right = right->maxLengths(centre,k,maxn,data);
    
    //std::vector<std::pair<int,double> > Right=right->maxLengths(centre,k,maxn,data);
    return Left; 
  }
  void RecurseUpDistances() {
    assert(!isleaf());
    left->distance[0][0] = std::min(distance[0][2],distance[0][0])+1;
    left->distance[1][0] = std::min(distance[1][2],distance[1][0])+1;    
    if (!left->isleaf())   left->RecurseUpDistances();
    
    right->distance[0][0] = std::min(distance[0][0],distance[0][1])+1;
    right->distance[1][0] = std::min(distance[1][0],distance[1][1])+1;    
    if (!right->isleaf())   right->RecurseUpDistances();  
  }

  
public:
  /** The data                                                                        */
  trinode *left;
  trinode *centre;
  trinode *right;
  trinode *up;
  int position;     // position of the SNP that causes this split
  std::vector<int> labels;
  TNT::Array2D<int> distance;   /** Distances for tree distance calculations    */
  int localCC[2];
};



/** A class to step through a tree in Node - LEFT - RIGHT order 
 * (This is known as preorder tree traversal in CS literature)
 * which is the same order that the ape phylo class seems to use
 * so this should enable us to link statistics calculated 
 * on these internal nodes to tree graphs produced in the R ape class
 * This returns a 0 when it gets to the last node
 */
template <class T> 
class NLRIterator {
public:
  /** Construct from the root                                  */
  NLRIterator(T *root) {
    current = root;
  }
  /** Step through one unit                                    */
  T *operator++() {
    next();
    return current;
  }
  /** Move the iterator to the next leaf and return it         */
  T *nextLeaf() {
    for(;;) {
      next();
      if (current==0||current->isleaf()) return current;
    }
  }
  /** Move the iterator to the next internal node and return   */
  T *nextInternal() {
    for(;;) {
      next();
      if (current==0||!current->isleaf()) return current;
    }
  }
  /** Is the iterator at the end?                              */
  bool isend() const {
    return (current==0);
  }
  /** Return the object that the iterator is pointing to       */
  T *operator *() const {	 
    return current;
  }
private:
  /** Get the next node                                 */
  void next() {
    if (current->isleaf()) {
      if (visited.size()==0) {
	current=0;
      } else {
	current=visited.top()->right;
	visited.pop();
      }
    } else {
      visited.push(current);
      current=current->left;
    }
  }
  /** A pointer to the current T                              */
  T* current;
  stack<T *> visited;  /** which nodes have already been visited */
};

/** A class to step through a tree in LEFT - RIGHT - Node order 
 * (commonly known by computer scientists as postorder tree traversal)
 */
template <class T> 
class LRNIterator {
public:
  /** Construct from the root                                  */
  LCRNIterator(T *start) {
  	// need to step from the start until we get get to the first
  	// leaf, keeping track of where we are.  Note that this does not
  	// have to use the whole tree, and can just look below a single node
    current=start;
  	while (!current->isleaf()) {
      visited.push(std::pair<T*,bool>(current,false));
      current=current->left;	
  	}
  }
  /** Step through one unit                                    */
  T *operator++() {
    next();
    return current;
  }
  /** Move the iterator to the next leaf and return it         */
  T *nextLeaf() {
    for(;;) {
      next();
      if (current==0||current->isleaf()) return current;
    }
  }
  /** Move the iterator to the next internal node and return   */
  T *nextInternal() {
    for(;;) {
      next();
      if (current==0||!current->isleaf()) return current;
    }
  }
  /** Is the iterator at the end?                              */
  bool isend() const {
    return (current==0);
  }
  /** Return the object that the iterator is pointing to       */
  T *operator *() const {	 
    return current;
  }
private:
  /** Go to the next node                                 */
  void next() {
    if (visited.size()==0) {
      current=0;
  	} else if (visited.top().second==true) {
      current=visited.top().first;
      visited.pop();
   	} else {
      visited.top().second=true;
      current=visited.top().first->right;
      while (!current->isleaf()) {
        visited.push(std::pair<T*,bool>(current,false));
        current=current->left;
      }
  	}
  }
  /** A pointer to the current T                              */
  T* current;
  stack<std::pair<T *, bool> > visited;  /** which nodes have already been visited */
};


#endif
