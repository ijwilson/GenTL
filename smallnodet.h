#ifndef SMALLNODE_T_H__
#define SMALLNODE_T_H__

#include <iostream>
#include <vector>

namespace GenTL {

  /** \brief template node class.  
   *
   * A class for buildins simple trees.  the class for extracting trees
   * from rnodes is built directly from this class
   */

  template <class T>
  class node {
    //  friend class tree;
    //    friend class iterator;
    typedef node* nodeptr; /// pointer to a node
  public:
    node():left(0),right(0),time(0.0),data_(){ };                             /// default constructor
    node(int lab,T data):left(0),right(0),time(0.0),label(lab),data_(data){}; /// leaf constructor
    node(T data):left(0),right(0),time(0.0),data_(data){ };                   /// leaf constructor
	//    node(double t):left(0),right(0),time(t){};                      /// leaf constructor 
    node (double t, nodeptr l, nodeptr r)
      :left(l),right(r),time(t) {
      left->up=this;
      right->up=this;
    };
    node (double t, nodeptr l, nodeptr r,int lab,T &a)
      :left(l),right(r),time(t),label(lab),data_(a){
      left->up=this;
      right->up=this;
    };
    node (double t, nodeptr l, nodeptr r,T &data)
      :left(l),right(r),time(t),data_(data){
      left->up=this;
      right->up=this;
    };
    
    void assign(double t, nodeptr l, nodeptr r,T &a){
      left=l;
      right=r;
      time=t; // assumes data already correct
      left->up=this;
      right->up=this;
      data_=a;
    }
    void recursiveprint(std::ostream &o,bool quote=true
			, bool internalnodes=false,bool label=true) const; 
    void recursiveprintnodata(std::ostream &o,bool printtimes=true,bool internalnodes=true) const;
    ///recursively print the tree in Newick format.
    int classifyancestralnode(const std::vector<int> &location);
    int count_leaves() const {
      if (left==0) return 1;
      else return left->count_leaves()+right->count_leaves();
    }
    
    /** calculate the length of tree below a node */
    double length() const {
      if (left==0) return 0.0;
      double t=2.*time-left->time-right->time;
      return t+right->length()+left->length();
    }
    /** is this node a leaf ?                                                   */
    bool isleaf() const {
      return (left==0);
    } 
    /** Find a node with this label                  */
    node<T> *find_node(int lab) {
      if (isleaf()) {
	if (label==lab) return this;
	return 0;
      }
      node <T> *ret=left->find_node(lab);
      if (ret==0) ret=right->find_node(lab);
      return ret;
    };

    void rotate() {
      nodeptr tmp=left;
      left=right;
      right=tmp;
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
    void get_leaves(std::vector<node<T> *> &leaves) {
       if (left==0) {
	 leaves.push_back(this);
       } else {
	 left->get_leaves(leaves);
	 right->get_leaves(leaves);
       }
       return;
    } 
/** get the list of internal node below (and including) a node       */
    void get_internal(std::vector<node<T> *> &a) {
       if (left!=0) {
	 a.push_back(this);
	 left->get_internal(a);
	 right->get_internal(a);
       } 
       return;
    }
    void write_leaves(int pos);

    double recursemutate(double mtime ,const int pos);

    void nodefreq(int *nf, int npops);

    void print(std::ostream &o){ o << data_;}
    nodeptr left;  /// pointer to "left" node 
    nodeptr right; /// pointer to "right" node
    nodeptr up;  /// pointer to node above 
    std::vector<int> countleft;
    std::vector<int> countright;
    double time;   /// coalescence time
    int label;     /// label for the node
    const T &data() const {return data_;}; ///extract the data (but not change it).
    T &data() { return data_;}; //extract the data
  
    //  friend node<T> readnode(std::istream &in);
  protected:
    T data_;
  };

  template <typename T>
    void recursively_destroy_node(node <T> *here) {
    if (!here->isleaf()) {
      recursively_destroy_node(here->left);
      recursively_destroy_node(here->right);
    }
    delete here;
  }

  template <class T>
    void node<T>::recursiveprint(std::ostream &o, bool quote, bool internalnodes, bool label) const 
  {
    if (left!=0) {
      o << "(";
      left->recursiveprint(o,quote,internalnodes,label);
      o << ":" << time-left->time << ",";
      right->recursiveprint(o,quote,internalnodes,label);
      o << ":" << time-right->time << ")";
      if (internalnodes) {
	if (quote) o << "'";
	if (label) 
	  o << label << " " << data_;
	else 
	  o << data_;
	if (quote) o << "'";
      }
    } else {
      if (quote) o << "'";
      if (label) 
	o << label << " " << data_;
      else o << data_;
      if (quote) o << "'";
    }
  }

 template <class T>
   void node<T>::recursiveprintnodata(std::ostream &o,bool printtime,bool internalnodes) const  {
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
}
#endif
