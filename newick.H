/**                        @file               */
#ifndef NEWICK_H
#define NEWICK_H
/***********************************************************************/
/*  charnode is a structure that is used to read from Newick files     */
class charnode {
 public:
  charnode():d1(0),d2(0),time(0.0) {
  }
  /// delete this tree recursively
  void destroy_chartree();
  charnode *d1;   //< descendent 1
  charnode *d2;   //< descendent2
  char *val;      //< the string at a node
  int len;        //< the length of the string
  double time;    //< the node time
};

class node;
class tree;

#include <iosfwd>
#include <vector>
#include <string>

using std::string;
using std::vector;

void write_Newick(node *root, node *sample, const char *filename, std::ostream &o, int npop,
		  int ninf, int nstr, bool label);//, const std::vector<int> &whmodel);

std::ostream &write_Newickshape(node *root, node *sample, std::ostream &o);
charnode *readcharnodeutil(std::istream &in, int *count);

int getposition(const char *info);
double getproportion(const char *info);
int getlocation(char *info);
void gettreeinfo(int &ninf,int &nstr, char *info,int len);
node *convertcharnodesample(node *sample, charnode *any, 
			    int nstr, int ninf, int posanc);

void write_Newick_label(node *root, node *sample, const char *filename, std::ostream &o
			,  const vector<string> &labels);

#endif
