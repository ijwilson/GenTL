/********************************************************************
 *   @file
 *   Time-stamp: <2011-07-19 16:47:52 ijw>
 *
 ********************************************************************/
#include <iostream>
#include <fstream>
#include <vector>

#include "util.H"
#include "newick.H"
#include "poptree.H"
#include "time.H"
#include "node.H"
#include "lhood.H"
#include "muts.H"

#include "newio.H"


/******************************************************************/
void charnode::destroy_chartree()
{
  if (d1!=0) d1->destroy_chartree();
  if (d2!=0) d2->destroy_chartree();
  delete val;
  delete this;
}

charnode *readcharnodeutil(std::istream &in, int *count)
/* reads in a tree from a file and the mutation and other
   information contained in a Newick file to enable the program to
   be restarted */
{
  charnode *here=new charnode;
  int ch = skipspace(in);

  if (ch == '(') {
    here->d1=readcharnodeutil(in,count);
    ch = skipspace(in);
    if (ch != ':') {
      std::cerr << "error expected a colon 1: got " << ch << std::endl;
      exit(EXIT_FAILURE);
    }
    double ttime;
    in >> ttime;
    here->time = ttime + here->d1->time;
    ch=skipspace(in);
    if (ch != ',')
      std::cerr << "error expected a comma got " << ch << std::endl;
    here->d2= readcharnodeutil(in,count);
    ch=skipspace(in);
    if (ch != ':') {
      std::cerr<< "error expected a colon 2: got " << ch << std::endl;
      exit(0);
    }
    in >> ttime;
    here->time = ttime + here->d2->time;
    ch=skipspace(in);
  } else {
    in.unget();
    *count +=1;
  }
  here->val=readfromquotes(in,&(here->len));
  return here;
}
/*********************************************/
int getposition(const char *info)
{
  int tmp;
  sscanf(info,"%d:",&tmp);
  return tmp;
}
/*********************************************/
double getproportion(const char *info)
{
  double tmp;
  int i;
  for (i=0;;i++) {
    if (info[i]=='~') break;
  }
  sscanf(info+i+1,"%lg ",&tmp);
  return tmp;
}
/*********************************************/
int getlocation(char *info)
{
	int tmp,i;
	for (i=0;;i++) if (info[i]=='<') break;
	sscanf(info+i+1,"%d>",&tmp);
	return tmp;
}
/**********************************************/
int  *getgenotype(int ninf,int nstr, char *info)
{
  int count=0,i,*gen;

  gen= new int[1+ninf+nstr];
  for (i=0;;i++) if (info[i]=='>') break;
  count=i;
  for (i=count+1;;i++) if (info[i]!=' ') break;
  count=i;
  for (i=0;i<ninf;i++) gen[i+1] = info[count+i]-48;
  if (ninf && info[count+i]!='~')
    error("should be a ~ in getgenotype");
  count+=ninf+1;
  for (i=1;i<nstr;i++) {
    sscanf(info+count,"%d-",gen+ninf+i);
    for (;;)
      if (info[count++]=='-') break;
  }
  sscanf(info+count,"%d-",gen+ninf+nstr);
  return gen;
}
/**********************************************/
void gettreeinfo(int &ninf,int &nstr, char *info,int len)
{
  int count=0,i;
  ninf=nstr=0;
  for (i=0;;i++) {
    if (info[i]=='>') break;
  }
  i++;
  for (;;i++) {
    if (info[i]=='~') break;
    else if (isalpha(info[i])||isdigit(info[i])) count++;
  }
  ninf=count;count=0;
  i++;
  for (;;) {
    for (;i<len;i++) {
      if (!isspace(info[i])) {
	count+=1;
	i++;
	break;
      }
    }
    for (;i<len;i++) if (info[i]=='-') {
      i++;
      break;
    }
    if (i==len) break;
  }
  nstr=count;
}
/**********************************************/
#ifndef NOINF
node *convertcharnodesample(node *sample, charnode *any,
						   int nstr, int ninf, int posanc)
{
	int p;

	p=getposition(any->val);
	sample[p].infgeno=getgenotype(ninf,nstr,any->val);
	sample[p].STRgeno=sample[p].infgeno+ninf;
	sample[p].location=getlocation(any->val);
	sample[p].time=any->time;
	if (posanc>0)
	sample[p].ancestor=sample+posanc;
	else sample[p].ancestor=NULL;

	if (any->d1==0) {
		sample[p].desc_left=NULL;
		sample[p].desc_right=NULL;
	} else {
		sample[p].desc_left=
			convertcharnodesample(sample,any->d1,nstr,ninf,p);
		sample[p].desc_right=
			convertcharnodesample(sample,any->d2,nstr,ninf,p);
	}
	return sample+p;
}
#endif

/***********************************************************************/
void writegeno(std::ostream &o,int *geno,int n,const std::vector<int> &type)
{
   int sep=0;
  for (int i=1;i<=n;i++) {
//    switch(type[i]) {
//    case KALLELES:
//    case SS:
      if (sep==1) o << "-";
      o << geno[i];break;
      sep=1;
//    default:
//      sep=0;
//      switch(geno[i]) {
//      case A:o << "A";break;
//      case G:o << "G";break;
//      case T:o << "T";break;
//      case C:o << "C";break;
//      default:
//	std::cerr << "error in print base?, what is " << char(geno[i]) << std::endl;
//	exit(EXIT_FAILURE);
//      }
 //   }
  }
}
/***********************************************************************/
void writenode(std::ostream &o,node *any,
	       int npop, int ninf, int nstr, bool label,node *samp)
{

  o << "'";
  if (label) {
    o << any-samp << ": ";
  }
  if (npop>1)  o << "<" << any->location << "> ";

  if (ninf>0) {
    for (int i=1;i<=ninf;i++)
      o << any->infgeno[i];
    o << "~";
  }

  for (int i = 1; i < nstr; i++)
    o << any->STRgeno[i] << "-";
  if (nstr>0) o << any->STRgeno[nstr];

  o << "'";
}
/***********************************************************************/
void writelabelutil(node *anynode, std::ostream &o, node *sample
		    , const std::vector<std::string> &labels)
{
  if (anynode->desc_left != NULL) {
    o << "(";
    writelabelutil(anynode->desc_left, o,sample,labels);
    o << ":" << std::setw(10) << std::fixed << anynode->time - anynode->desc_left->time << ",";
    writelabelutil(anynode->desc_right, o,sample,labels);
    o << ":" <<  std::setw(10) << std::fixed<<  anynode->time - anynode->desc_right->time << ")";
  } else {
    assert(anynode-sample-1>=0&&anynode-sample-1<static_cast<int>(labels.size()));
    o << "'" << labels.at(anynode-sample-1) << "'";
  }
}
/***********************************************************************/
void write_Newick_label(node *root, node *sample, const char *filename, std::ostream &o
			, const std::vector<std::string> &labels)
{
  if (filename==NULL) {
    writelabelutil(root, o, sample,labels);
    o << " ;" << std::endl;
  } else {
    std::ofstream of(filename);
    writelabelutil(root, of, sample,labels);
    of << " ;";
    if (filename)of.close();
  }
  return;
}
/***********************************************************************/
void writeutil(node *anynode,std::ostream &o, int npop,
	       int ninf, int nstr, bool label, node *sample)
{
  if (anynode->desc_left != NULL) {
    o << "(";
    writeutil(anynode->desc_left, o,npop,ninf,nstr,label,sample);
    o << ":" <<  std::setw(10) << std::fixed<< anynode->time - anynode->desc_left->time << ",";
    writeutil(anynode->desc_right, o,npop,ninf,nstr,label,sample);
    o << ":" <<  std::setw(10) << std::fixed<< anynode->time - anynode->desc_right->time << ")";
  }
  writenode(o,anynode,npop,ninf,nstr,label,sample);
}
/***********************************************************************/
void write_Newick(node *root, node *sample, const char *filename, std::ostream &o, int npop,
		  int ninf, int nstr, bool label)  //, const std::vector<int> &whmodel)
{
  std::ofstream of;
  if (filename!=NULL) of.open(filename);
  std::ostream &output=(filename==NULL)?o:of;

  writeutil(root,output,npop,ninf,nstr,label,sample);
  output << " ;";
  if (filename) of.close();
  return;
}
/***********************************************************************/
void writeshapenode(std::ostream &o,node *any,node *samp)
{
  if (any->desc_left==0) o << "'" << any-samp << "'";
}
/***********************************************************************/
void writeshape(node *anynode, std::ostream &o, node *sample)
{
  if (anynode->desc_left != 0) {
    o << "(";
    writeshape(anynode->desc_left, o,sample);
    o << "," ;
    writeshape(anynode->desc_right, o,sample);
    o << ")";
  }
  writeshapenode(o,anynode,sample);
}
/***********************************************************************/
 std::ostream &write_Newickshape(node *root, node *sample, std::ostream &o)
{
  writeshape(root, o,sample);
  o << " ;\n";
  return o;
}
/***********************************************************************/
