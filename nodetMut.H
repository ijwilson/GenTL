/** @file  */ 
#ifndef __NODE__H__T_IJW_MUT
#define __NODE__H__T_IJW_MUT

#include <stack>
#include <vector>
#include <iostream>
#include <stdexcept>
#include "gsl_rand.H"
#include "gsl_distributions.H"
#include "coalescent.H"

/** This is aa different version of the nodet.H header file that has been altered to 
 * allow for generic mutation functors - so that we just need to pass a class that
 * knows how to mutate itself and we can simulate using the machinery that is 
 * already there.  At the moment this is just for the non-arg classes but I 
 * should be able to extend it for the other classes.  I shall include 
 * an example mutable class to show how it can be done.

 * any mutate class needs to be able to 
 * copy itself to a new location 
 * A copy constructor
 * mutate using/not using  a vector of mutation parameters
 * print itself
 * I have got rid of all the specialisations  - they are no longer needed
 */




namespace GenTL {

  /** \brief template node class.  
   *
   * A class for building simple trees.  the class for extracting trees
   * from rnodes is built directly from this class
   */

  template <class MUTABLE>
  class node {
    typedef node* nodeptr; /// pointer to a node
  public:
    node()                                        /// default constructor
      :left(0),right(0),time(0.0),data_(){ };                            
    node(int lab,MUTABLE data)                          /// leaf constructor
      :left(0),right(0),time(0.0),label(lab),data_(data){}; 
    node(int lab)                          /// leaf constructor
      :left(0),right(0),time(0.0),label(lab){}; 
    node(MUTABLE data)                     /// leaf constructor
      :left(0),right(0),time(0.0),data_(data){ };                   
    node (double t, nodeptr l, nodeptr r)
      :left(l),right(r),time(t) {
      left->up=this;
      right->up=this;
    };
    node (double t, nodeptr l, nodeptr r,int lab,MUTABLE &a)
      :left(l),right(r),time(t),label(lab),data_(a){
      left->up=this;
      right->up=this;
    };
    node (double t, nodeptr l, nodeptr r,MUTABLE &data)
      :left(l),right(r),time(t),data_(data){
      left->up=this;
      right->up=this;
    };
    
    void assign(double t, nodeptr l, nodeptr r,MUTABLE &a){
      left=l;
      right=r;
      time=t; // assumes data already correct
      left->up=this;
      right->up=this;
      data_=a;
    }
    /** Recursively find the position of the mutation              */
    void recursemutate(const std::vector<double> &theta, rng &r) {
      left->data_ = data_.mutate(theta,time-left->time,r);
      right->data_ = data_.mutate(theta,time-right->time,r);
      if (!left->isleaf()) {
        left->recursemutate(theta,r);
      }
      if (!right->isleaf()) {
        right->recursemutate(theta,r);
      }  
    }


    void recursiveprint(std::ostream &o,bool quote=true, bool internalnodes=false,bool labels=false) const; 
    void recursiveprintnodata(std::ostream &o,bool printtimes=true,bool internalnodes=true) const;
    ///recursively print the tree in Newick format.
    int classifyancestralnode(const std::vector<int> &location);
    int count_leaves() const {
      if (left==0) return 1;
      else return left->count_leaves()+right->count_leaves();
    }
    /** get a vector of labels below a node       */
    void leaf_labels(std::vector<int> &leaves) {
       if (left==0) {
	 leaves.push_back(label);
       } else {
	 left->leaf_labels(leaves);
	 right->leaf_labels(leaves);
       }
       return;
    }
    /** get the list of leaves below a node       */
    void get_leaves(std::vector<node<MUTABLE> *> &leaves) {
       if (left==0) {
	 leaves.push_back(this);
       } else {
	 left->get_leaves(leaves);
	 right->get_leaves(leaves);
       }
       return;
    }
    /** change everything at position pos to value   */
    void change_below(int pos, char value);


    /** get the list of internal node below (and including) a node       */
    void get_internal(std::vector<node<MUTABLE> *> &a) {
       if (left!=0) {
	 a.push_back(this);
	 left->get_internal(a);
	 right->get_internal(a);
       } 
       return;
    }
    /** calculate the length of tree below a node */
    double length() const {
      if (left==0) return 0.0;
      double t=2.*time-left->time-right->time;
      return t+right->length()+left->length();
    }
    /** is this node a leaf ?                                            */
    bool isleaf() const {
      return (left==0);
    } 
    /** Find a node with this label                  */
    node<MUTABLE> *find_node(int lab) {
      if (isleaf()) {
	if (label==lab) return this;
	return 0;
      }
      node <MUTABLE> *ret=left->find_node(lab);
      if (ret==0) ret=right->find_node(lab);
      return ret;
    };
    void write_leaves(int pos);

   
    void nodefreq(int *nf, int npops);

    void print(std::ostream &o){ o << data_;}
    nodeptr left;  /// pointer to "left" node 
    nodeptr right; /// pointer to "right" node
    nodeptr up;  /// pointer to node above 
    std::vector<int> countleft;
    std::vector<int> countright;
    double time;   /// coalescence time
    int label;     /// label for the node
    const MUTABLE &data() const {return data_;}; ///extract the data (but not change it).
    MUTABLE &data() { return data_;}; //extract the data
  
  protected:
    MUTABLE data_;
  };

  /** forward declaration of iterator                       */
  template <class MUTABLE> class iterator; 
  /** tree class                                            */
  template <class MUTABLE> 
  class tree {
  public:
    /** default constructor                        */
    tree() {
    };
    /** constructor for n haplotypes               */
    tree(int n):sample_(n),ancestors_(n-1) {
      root = &(ancestors_[n-2]);
    };
    /** constructor for a coalescent               */
    tree(coalescent *c):
      sample_(c->nleft())
      ,ancestors_(c->nleft()-1) {
      root=&(ancestors_[c->nleft()-2]);
      if (c->nleft()>1)
        simulate(c);
    }; 
    /** constructor for a branching process               */
    tree(int n, ctsdistribution &nexteventprior, rng &r):
      sample_(n),ancestors_(n-1) {
      root=&(ancestors_[n-2]);
      bpsimulate(n,nexteventprior,r);
    };

    /** declaration for the simulation             */		
    void simulate(coalescent *c);
    /** declaration for the branching process simulation             */
    void bpsimulate(int n, ctsdistribution &nexteventprior, rng &r);
    /** iterate through the tree                   */
    iterator<MUTABLE> begin() {
      iterator<MUTABLE> i(*this);
      return i;
    };
    /** remakes a tree from the root */
    node<MUTABLE> *root;
 
    void getnodefreq(std::vector<int> &location, int npops);
    /** Count to the left and right from without locations    */
    void getnodefreq();

 
    node<MUTABLE> *readnode(std::istream &in, std::vector<node<MUTABLE> > &anc
		      , std::vector<node<MUTABLE> > &sample);
    
    iterator<MUTABLE> end();

    node<MUTABLE> *first() {
      return &sample_[0];
    }

    const node<MUTABLE> *ancestorptr(int i) const {
      return &ancestors_[i];
    }

    node<MUTABLE> *ancestorptr(int i)  {
      return &ancestors_[i];
    }

    const node<MUTABLE> &ancestor(int i) const {
      return &ancestors_[i];
    }

    node<MUTABLE> sample(int i) const {
      return sample_[i];
    }

    const MUTABLE &data(int i) const {
      return sample_[i].data();
    }

    MUTABLE &data(int i) {
      return sample_[i].data();
    }

    double anctime(int i) const {
      return ancestors_[i].time;
    }

    /** How many leaves dows the tree have ?                                  */
    size_t nleaves() const {
      return sample_.size();
    }
    /** mutate - template  specialisation for different types of node         */

    void mutate(int nloci, const std::vector<double> &theta,rng &r,bool randomStart=false) {
      root->data().initial(nloci,r,randomStart);
      root->recursemutate(theta,r);
    }

    /** mutate - mutate, but start from "above the root"                      */
    void mutate(int nloci, const std::vector<double> &theta,rng &r, double toptime,bool randomStart=false) {
      MUTABLE m;
      m.initial(nloci,r,randomStart);
      root->data()=m.mutate(theta,toptime-root->time,r);
      if (nleaves()>1)
        root->recursemutate(theta,r);
    }



    /** friend to the input operator                  */
    template<class V>
    friend  std::istream &operator>>(std::istream &in,  tree<V>  &a);
  protected:
    std::vector<node<MUTABLE> > sample_;
    std::vector<node<MUTABLE> > ancestors_;
  private:
    tree(const tree &a);
  };
  /** iterator class
   *
   */
  template <class MUTABLE> class iterator {
  public:
    iterator(const tree<MUTABLE> &t) {
      add(t.root);
    };
    node<MUTABLE> *operator++();

    node<MUTABLE> *operator *() {	 
      return where_.top();
    };
    bool operator != (const node<MUTABLE> *tmp) {
      if (where_.empty())  return (tmp!=0);
      else return (where_.top()!=tmp);
    };
               
  private:
    void add(node<MUTABLE> *here){
      for (;;) {
        if (here==0) break;
        where_.push(here);
        here=here->left;
      } 
    };
    /** we keep a stack of the  nodes.  MUTABLEhe current position
     * of the iterator is given by where_.top() */
    std::stack<node<MUTABLE> *> where_;
  };

  /**
   * output a tree on ostream 
   */
  template <class MUTABLE>
  std::ostream &operator<<(std::ostream &o,const  node <MUTABLE> &root) {
    root.recursiveprint(o);
    return o;
  }
  /**
   * output a tree on ostream 
   */
  template <class MUTABLE>
  std::ostream &operator<<(std::ostream &o,const  tree<MUTABLE> &t) {
    t.root->recursiveprint(o,true,false);
    o << ";";
    return o;
  }

  /** recursive read nodes         */
  template <class MUTABLE>
  node<MUTABLE> *readnode(std::istream &in, std::vector<node<MUTABLE> > &anc, std::vector<node<MUTABLE> > &sample
		    , bool quote=true, bool internalnodes=false,bool uselabel=true)
  {
    int ch=skipblank(in);
    if (ch=='(') {
      node<MUTABLE> *l,*r;
      l=readnode<MUTABLE>(in,anc,sample,quote,internalnodes,uselabel);
   
      checkreadcharacter(in,':',"after left node");
      
      double ltime,rtime;
      in >> ltime;
      
      checkreadcharacter(in,',',"comma after time");  
            
      r=readnode<MUTABLE>(in,anc,sample,quote,internalnodes,uselabel);
   
      checkreadcharacter(in,':',"reading right node");
      in >> rtime;

      if (fabs(ltime+l->time-rtime- r->time)>1E-5) 
	throw std::domain_error("error with times - I expected a ... tree");     
      
      checkreadcharacter(in,')',"error after reading right node");
      int lab(0);
      MUTABLE d(0);
      
      if (internalnodes) {
	ch=skipblank(in);
	if (quote) {
	  if (ch!='\'') {
	    std::cerr << "warning - no information at internal node, trying to continue";
	    std::cerr << " putting back " << char(ch) << std::endl;
	    in.putback(ch);
	  } else {
	    if (uselabel) in >> lab;
	    in >> d;
	    checkreadcharacter(in,'\'',"error - expected quote");  
	  }
	} else {
	  in.putback(ch);
	  if (uselabel) in >> lab;
	  in >> d;
	}
      }
      anc.push_back(node<MUTABLE>(ltime+l->time,l,r,lab,d)); 
      return &anc[anc.size()-1];
    } else {
      int lab(0);
      MUTABLE d(0);
      if (quote) {
	if (ch!='\'') {
	  std::cerr << "warning - no information at node, trying to continue";
	  in.putback(ch);
	} else {
	  if (uselabel) in >> lab;
	  in >> d;
	  checkreadcharacter(in,'\'',"error - expected quote");  
	}
      } else {
	in.putback(ch);
	if (uselabel) in >> lab; 
	in >> d;
      }
      sample.push_back(node<MUTABLE>(lab,d));
      return &sample[sample.size()-1];
    }
  }
  /**
   * input a tree from stream
   */
  template <class MUTABLE>
  std::istream &operator>>(std::istream &in,  tree <MUTABLE> &a) {
    a.sample_.resize(0);
    a.sample_.reserve(1000);
    a.ancestors_.reserve(1000);
    a.ancestors_.resize(0);
    a.root=readnode<MUTABLE>(in,a.ancestors_,a.sample_,true,true);
    return in;
  }

  template <class MUTABLE>
  node<MUTABLE> *iterator<MUTABLE>::operator++() {
    if (where_.empty()) return 0;
    node<MUTABLE> *curr=where_.top(); 
    where_.pop();
    if (curr->right !=0) {
      add(curr->right);
    }
    return curr;
  };
  template <class MUTABLE>
  void node<MUTABLE>::recursiveprint(std::ostream &o, bool quote, bool internalnodes,bool labels) const 
  {
    if (left!=0) {
      o << "(";
      left->recursiveprint(o,quote,internalnodes,labels);
      o << ":" << time-left->time << ",";
      right->recursiveprint(o,quote,internalnodes,labels);
      o << ":" << time-right->time << ")";
      if (internalnodes) {
	if (quote) o << "'";
	if (labels) o << label << " " ;
      o <<  data_;
	if (quote) o << "'";
      }
    } else {
      if (quote) o << "'";
      if (labels) o << label << " " ;
      o << data_;
      if (quote) o << "'";
    }
  }

 template <class MUTABLE>
   void node<MUTABLE>::recursiveprintnodata(std::ostream &o,bool printtime,bool internalnodes) const  {
    if (left!=0) {
      o << "(";
      left->recursiveprintnodata(o,printtime,internalnodes);
      if (printtime) o << ":" << time-left->time;
      o << ",";
      right->recursiveprintnodata(o,printtime,internalnodes);
      if (printtime) o << ":" << time-right->time;
      o << ")";
      if (internalnodes) o << label;
    } else {
      o << label;
    }
  }
  /** What sort of nodes is this node ancestral to?
   * types of location i only ( returns i) 
   * or to a mixture of different nodes - returns -1        */
  template <class MUTABLE>
  int node<MUTABLE>::classifyancestralnode(const std::vector<int> &location)
  {
    int l,r;
    if (left!=0) {
      int l= left->classifyancestralnode(location);
      if (l==-1) return -1;
      int r=right->classifyancestralnode(location);
      if (l==r) return l;
      else return -1;
    }
    else {
      return location[label];
    }
  }
  /** A utility function to do the simulations for all constructors                             */
  template <class MUTABLE>
  void tree<MUTABLE>::simulate(coalescent *c)
  {
    double t=0.0;
    std::vector<node<MUTABLE> *> lines(c->nleft());
    for (int i=0;i<c->nleft();i++) lines[i]=&sample_[i];
    int pos=0;    
    for (;;) {
      double NextEvent=c->next(t);
      if (NextEvent>1E5) {
	std::cerr << "Problem with coalescent: " << NextEvent << std::endl;
	exit(0);
      }
      t += NextEvent;
     
      std::pair<int,int> chromo = c->coalchromo();
      if (chromo.first!=chromo.second) { // a coalescence 
        ancestors_[pos].left=lines[chromo.first];
        ancestors_[pos].right=lines[chromo.second];
        ancestors_[pos].time=t;
        lines[chromo.second]=&ancestors_[pos];
        lines[chromo.first]=lines[c->nleft()-1];
        pos++;
        c->UpdateCoalescence(chromo,1);
        if (c->nleft()==1) break;
      }
    }
    return;
  }
  template <class MUTABLE> 
  void tree<MUTABLE>::getnodefreq(std::vector<int> &location, int npops)
  {
    for (int i=0;i<ancestors_.size();i++) {
      ancestors_[i].countleft.resize(npops,0);
      ancestors_[i].countright.resize(npops,0);
      if (ancestors_[i].left->left==0) { // a leaf
	ancestors_[i].countleft[location[&sample_[0]-ancestors_[i].left]] =1;
      } else {
	for (int j=0;j<npops;j++) ancestors_[i].countright[j]=
	  ancestors_[i].right->countleft[j]+ancestors_[i].right->countright[j];
      }
      if (ancestors_[i].right->left==0) { // a leaf
	ancestors_[i].countright[location[&sample_[0]-ancestors_[i].right]] =1;
      } else {
	for (int j=0;j<npops;j++) ancestors_[i].countleft[j]=
	  ancestors_[i].left->countleft[j]+ancestors_[i].left->countright[j];
      }
    }
  }
 template <class MUTABLE> 
  void tree<MUTABLE>::getnodefreq()
  {
    for (size_t i=0;i<ancestors_.size();i++) {
      ancestors_[i].countleft.resize(1,0);
      ancestors_[i].countright.resize(1,0);
      if (ancestors_[i].left->left==0) { // a leaf
	ancestors_[i].countleft[0] =1;
      } else {
	ancestors_[i].countright[0]=
	  ancestors_[i].right->countleft[0]+ancestors_[i].right->countright[0];
      }
      if (ancestors_[i].right->left==0) { // a leaf
	ancestors_[i].countright[0]=1;
      } else {
	 ancestors_[i].countleft[0]=
	  ancestors_[i].left->countleft[0]+ancestors_[i].left->countright[0];
      }
    }
  }
  template <typename MUTABLE>
    void recursively_destroy_node(node <MUTABLE> *here) {
    if (!here->isleaf()) {
      recursively_destroy_node(here->left);
      recursively_destroy_node(here->right);
    }
    delete here;
  }

 template <typename MUTABLE>
   int rotate_tree_labels(node <MUTABLE> *here) {
   if (here->isleaf()) return here->label;
   int lft=rotate_tree_labels(here->left);
   int rgt=rotate_tree_labels(here->right);
   if (lft>rgt) {
     node<MUTABLE> *tmp=here->left;
     here->left=here->right;
     here->right=tmp;
     return rgt;
   }
   return lft;
 }



}
#endif
