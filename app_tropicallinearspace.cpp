#include "vektor.h"
#include "printer.h"
#include "parser.h"
#include "gfanapplication.h"
#include "minors.h"
#include "field_rationals.h"
#include <vector>
#include <utility>
#include <sstream>
#include <algorithm>
#include "regularsubdivision.h"
#include "lp.h"

using namespace std;




static int64 fak(int a)
{
  if(a==0)return 1;
  return a*fak(a-1);
}

static int choose(int n, int d)
{
  return fak(n)/fak(n-d)/fak(d);
}



static int isInEdge(int v, pair<int, int> edge)
  {
    if(edge.first==v)return edge.second;
    if(edge.second==v)return edge.first;
    return -1;
  }
static int isLeaf(int v, list<pair<int,int> > edges) //returns parent
  {
    int c=0;
    int parent=0;
    for(list<pair<int,int> >::const_iterator i=edges.begin();i!=edges.end();i++)
      {
	if(isInEdge(v,*i)>=0)
	  {
	    c++;
	    parent=isInEdge(v,*i);
	  }
      }
    if(c<=1)return parent;
    return -1;
  }
static int isCherry(int v, list<pair<int,int> > edges)
  {
    int parent=isLeaf(v,edges);
    if(parent<0)return -1;
    int count=0;
    int other=-1;
    for(list<pair<int,int> >::const_iterator i=edges.begin();i!=edges.end();i++)
      {
	if(isInEdge(parent,*i)>=0)
	  {
	    int child=isInEdge(parent,*i);
	    //fprintf(stderr,"in edge v:%i p:%i c:%i childleaf?%i\n",v,parent,child,isLeaf(child,edges));
	    if((v!=child) && (isLeaf(child,edges)>=0))
	      {

		count++;
		other=child;
	      }
	  }
      }
    if(count==1)return other;
    return -1;
  }

list<string> treeStrings(int d, int n, IntegerVector const &v)
{
    list<string> ret;

    IntegerVector w(choose(n-1,d-1));
    for(int k=0;k<n;k++)
      {
	IntegerVectorList hyperSimplex;
	vector<int> I;
	IntegerVector p(n-1);
	for(int i=0;i<d;i++)I.push_back(1);
	for(int i=0;i<n-d;i++)I.push_back(0);
	int i=0,wi=0;
	do
	  {
	    if(I[k])
	      {
		w[wi++]=v[i];//CHANGING SIGN

		for(int j=0;j<k;j++)p[j]=I[j];
		for(int j=k+1;j<n;j++)p[j-1]=I[j];
		hyperSimplex.push_back(p);
	      }
	    i++;
	  }
	while(prev_permutation(I.begin(),I.end()));
	/*	    AsciiPrinter(Stderr).printVector(w);
	  AsciiPrinter(Stderr).printNewLine();
	*/
	IntegerMatrix hyperSimplex2=rowsToIntegerMatrix(hyperSimplex);
	set<set<int> > subd=regularSubdivision(hyperSimplex2,w);
	//  printSetSetInt(Stderr,subd);

	vector<char> vertexLabels(subd.size()+n-1);
	vector<IntegerVectorList> vertexCell(subd.size()+n-1);
	int J=0;
	for(int j=0;j<n;j++)
	  {
	    if(j!=k)
	      {
		vertexLabels[J++]=j+'1';//CHANGE THIS
		for(IntegerVectorList::const_iterator K=hyperSimplex.begin();K!=hyperSimplex.end();K++)
		  if((*K)[J-1])vertexCell[J-1].push_back(*K);
	      }
	  }
	for(set<set<int> >::const_iterator K=subd.begin();K!=subd.end();K++)
	  {
	    vertexLabels[J++]='I';
	    for(set<int>::const_iterator L=K->begin();L!=K->end();L++)
	      vertexCell[J-1].push_back(hyperSimplex2[*L]);
	  }
	/*
	  for(int i=0;i<vertexLabels.size();i++)
	  {
	  fprintf(Stderr,"\"%c\"\n",vertexLabels[i]);
	  AsciiPrinter(Stderr).printVectorList(vertexCell[i]);
	  }
	*/
	list<pair<int,int> > edges;

	for(int b=0;b<vertexLabels.size();b++)
	  for(int a=0;a<b;a++)
	    {
	      IntegerVectorList intersection;


	      /*		  AsciiPrinter(Stderr).printVectorList(vertexCell[a]);
		AsciiPrinter(Stderr).printVectorList(vertexCell[b]);
	      */
	      for(IntegerVectorList::const_iterator A=vertexCell[a].begin();A!=vertexCell[a].end();A++)
		for(IntegerVectorList::const_iterator B=vertexCell[b].begin();B!=vertexCell[b].end();B++)
		  if((*A-*B).isZero())intersection.push_back(*A);

	      //		  AsciiPrinter(Stderr).printVectorList(intersection);


	      //fprintf(Stderr,"RANK : %i\n",rankOfMatrix(intersection));
	      if(rankOfMatrix(intersection)==n-2)
		{
		  edges.push_back(pair<int,int>(a,b));
		}
	    }
	/*
	  for(list<pair<int,int> >::const_iterator j=edges.begin();j!=edges.end();j++)
	  {
	  fprintf(Stderr,"%i %i\n",j->first,j->second);
	  }
	*/
	set<int> allLeafs;
	set<int> allCherries;
	int numberOfCherries=0;
	for(int i=0;i<vertexLabels.size();i++)
	  {
	    //fprintf(Stderr,"%i isLeaf %i isCherry: %i\n",i,isLeaf(i,edges),isCherry(i,edges));
	    if(isCherry(i,edges)>=0)
	      {
		allCherries.insert(i);
		numberOfCherries++;
	      }
	    if(isLeaf(i,edges)>=0)allLeafs.insert(i);
	  }
	stringstream s;
	if(n==7 && numberOfCherries==6)
	  {
	    set<int> left=allLeafs;
	    //fprintf(Stdout,"S(");
	    s << "S(";
	    while(left.size())
	      {
		int a=*left.begin();
		int b=isCherry(a,edges);
		left.erase(a);
		left.erase(b);
		//fprintf(Stdout,"%c%c",vertexLabels[a],vertexLabels[b]);
		s << vertexLabels[a] << vertexLabels[b];
		if(left.size())
		  //fprintf(Stdout,",");
		  s << ",";
	      }
	    //fprintf(Stdout,")\n");
	    s << ")";
	  }
	else if(n==7 && numberOfCherries==4)
	  {
	    set<int> leftCherries=allCherries;
	    set<int> left=allLeafs;
	    int a=*leftCherries.begin();
	    int b=isCherry(a,edges);
	    leftCherries.erase(a);
	    leftCherries.erase(b);
	    int e=*leftCherries.begin();
	    int f=isCherry(e,edges);

	    int ab=isLeaf(a,edges);
	    set<int> allConnections;
	    for(list<pair<int,int> >::const_iterator i=edges.begin();i!=edges.end();i++)
	      if(isInEdge(ab,*i)>=0)allConnections.insert(isInEdge(ab,*i));
	    allConnections.erase(a);
	    allConnections.erase(b);
	    int C=*allConnections.begin();
	    int c=-1;
	    for(list<pair<int,int> >::const_iterator i=edges.begin();i!=edges.end();i++)
	      if(isInEdge(C,*i)>=0)
		if(isLeaf(isInEdge(C,*i),edges)>=0)c=isInEdge(C,*i);

	    left.erase(a);
	    left.erase(b);
	    left.erase(c);
	    left.erase(e);
	    left.erase(f);
	    int d=*left.begin();
	    //fprintf(Stdout,"C(%i%i,%i%i,%i%i)\n",a,b,c,d,e,f);
	    //		fprintf(Stdout,"C(%c%c,%c%c,%c%c)\n",vertexLabels[a],vertexLabels[b],vertexLabels[c],vertexLabels[d],vertexLabels[e],vertexLabels[f]);
	    s<<"C("<< vertexLabels[a] << vertexLabels[b] << "," << vertexLabels[c] << vertexLabels[d] <<","<<vertexLabels[e] << vertexLabels[f]<<")";
	  }
	else if(n==6 && numberOfCherries==4)
	  {
	    int a=*allCherries.begin();
	    allCherries.erase(a);
	    int b=isCherry(a,edges);
	    allCherries.erase(b);
	    int d=*allCherries.begin();
	    allCherries.erase(d);
	    int e=isCherry(d,edges);
	    allCherries.erase(e);
	    int c=0+1+2+3+4-a-b-d-e;
	    //fprintf(Stdout,"C(%c%c,%c,%c%c)\n",vertexLabels[a],vertexLabels[b],vertexLabels[c],vertexLabels[d],vertexLabels[e]);
	    s<<"C("<< vertexLabels[a] << vertexLabels[b] << "," << vertexLabels[c] <<"," <<vertexLabels[d] << vertexLabels[e] <<")";
	  }
	else
	  {
	    fprintf(Stderr,"I don't know how to print this tree\n");
	    assert(0);
	  }
	ret.push_back(s.str());
      }

    return ret;
  }



static int permutationIndex(vector<int> v)
{
  int ret=0;

  for(int i=0;i<v.size();i++)fprintf(stderr,"%i,",v[i]);

  while(prev_permutation(v.begin(),v.end()))ret++;

  fprintf(Stderr,"=%i\n",ret);
  return ret;
}


class TropicalLinearSpaceApplication : public GFanApplication
{
  IntegerOption dOption;
  IntegerOption nOption;
  SimpleOption boundaryTreesOption;
public:
  const char *helpText()
  {
    return "This program generates tropical equations for a tropical linear space in the Speyer sense given the tropical Pluecker coordinates as input.\n";
  }
  TropicalLinearSpaceApplication():
    dOption("-d","Specify d.",1),
    nOption("-n","Specify n.",1),
    boundaryTreesOption("--trees","list the boundary trees (assumes d=3)")
  {
    registerOptions();
  }
  const char *name()
  {
    return "_tropicallinearspace";
  }


  int main()
  {
    FileParser P(Stdin);

    int d=dOption.getValue();
    int n=nOption.getValue();

    assert(d<=n);

    IntegerVector v=P.parseIntegerVector();

    if(boundaryTreesOption.getValue())
      {
	list<string> s=treeStrings(d,n,v);
	for(list<string>::const_iterator i=s.begin();i!=s.end();i++)
	  fprintf(stderr,"%s\n",i->c_str());
      }
    else
      {
	vector<string> names=vectorVariableNames("x",n);
	vector<string> names2(names.size()+1);
	for(int i=0;i<names.size();i++)names2[i+1]=names[i];
	names2[0]="t";
	PolynomialRing R(Q,names2);
	PolynomialSet ret(R);

	if(d<n)
	  {
	    vector<int> I;
	    for(int i=0;i<d+1;i++)I.push_back(1);
	    for(int i=0;i<n-d-1;i++)I.push_back(0);
	    do
	      {
		Polynomial p(R);
		int sign=-1;
		for(int jr=0;jr<n;jr++)
		  if(I[jr])
		    {
		      int tExponent=1000;
		      //		  sign*=-1;
		      vector<int> J=I;
		      J[jr]=1-J[jr];
		      tExponent+=sign*v[permutationIndex(J)];

		      IntegerVector exponentVector(n+1);
		      exponentVector[jr+1]=1;
		      exponentVector[0]=tExponent;
		      p+=Term(R.getField().zHomomorphism(1),Monomial(R,exponentVector));
		}
		ret.push_back(p);
	      }
	    while(prev_permutation(I.begin(),I.end()));
	  }
	AsciiPrinter(Stdout).printPolynomialRing(R);
	AsciiPrinter(Stdout).printPolynomialSet(ret);
      }
    return 0;
  }
};

static TropicalLinearSpaceApplication theApplication;
