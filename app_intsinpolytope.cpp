#include "vektor.h"
#include "printer.h"
#include "parser.h"
#include "gfanapplication.h"
#include "intsinpolytope.h"
#include "lattice.h"
#include "determinant.h"

class IntsInPolytopeApplication : public GFanApplication
{
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  IntsInPolytopeApplication()
  {
    registerOptions();
  }
  const char *name()
  {
    return "_intsinpolytope";
  }

  // a is a squarematrix whose rows span a simplicial cone
  bool isInSimplicialCone(IntegerMatrix const &a, IntegerVector const &v)
  {

    for(int i=0;i<a.getHeight();i++)
      {
	IntegerMatrix b=a;
	b[i]=v;
	IntegerVectorList A=a.getRows();
	IntegerVectorList B=b.getRows();
	if((determinantSign(A)*determinantSign(B))<0)return false;
      }
    return true;
  }


  set<set<int> > partsOfZBases(set<set<int> > const &s, IntegerMatrix const &M)
  {
    set<set<int> > setsToCheck;

    for(set<set<int> >::const_iterator i=s.begin();i!=s.end();i++)
      {
	for(int j=0;j<M.getHeight();j++)
	  {
	    set<int> temp=*i;
	    temp.insert(j);
	    if(temp.size()!=i->size())
	      {
		bool isOK=true;
		/*		if(temp.size()==5)
		  {
		    IntegerMatrix MM(5,5);
		    int I=0;
		    for(set<int>::const_iterator k=temp.begin();k!=temp.end();k++,I++)
		      {
			MM[I]=M[*k];
		      }
		    IntegerVector v(5);
		    v[0]=2;//?????????????????????????????????????????????????
		    v[1]=2;
		    v[2]=2;
		    v[3]=17;
		    v[4]=8;
		    if(!isInSimplicialCone(MM,v))
		      isOK=false;
		      }*/

		if(isOK)
		  for(set<int>::const_iterator k=temp.begin();k!=temp.end();k++)
		    {
		      set<int> temp2=temp;
		      temp2.erase(*k);
		      if(s.count(temp2)==0)
			{
			  isOK=false;
			  break;
			}
		    }
		if(isOK)
		  {
		    static int c;
		    c++;
		    if(!(c&4095))fprintf(stderr,"tocheck:%i\n",c);

		    /*		    if(temp.size()==5)
		      {
			IntegerMatrix MM(5,5);
			int I=0;
			for(set<int>::const_iterator k=temp.begin();k!=temp.end();k++,I++)
			  {
			    MM[I]=M[*k];
			  }

			fprintf(stderr,"Candidate:\n");
			AsciiPrinter(Stderr).printVectorList(MM.getRows());
			}*/
		    setsToCheck.insert(temp);
		  }
	      }
	  }
      }

    fprintf(stderr,"Sets to test: %i\n",setsToCheck.size());

    set<set<int> > ret;
    for(set<set<int> >::const_iterator i=setsToCheck.begin();i!=setsToCheck.end();i++)
      {
	static int c;
	c++;
	if(!(c&4095))fprintf(stderr,"%i\n",c);
	IntegerVectorList l;
	for(set<int>::const_iterator j=i->begin();j!=i->end();j++)
	  l.push_back(M[*j]);
	if(isPartOfAZBasis(l))ret.insert(*i);
      }

    fprintf(stderr,"Produced sets: %i\n",ret.size());

    return ret;
  }

  int main()
  {
    FileParser p(Stdin);

    IntegerVectorList ivl=p.parseIntegerVectorList();
    IntegerMatrix A=rowsToIntegerMatrix(ivl);//.transposed();

    IntegerVector rightHandSide=p.parseIntegerVector();
        IntegerVector v=p.parseIntegerVector();

    AsciiPrinter P(stdout);

    fprintf(Stdout,"Lattice Kernel:\n");
    IntegerVectorList l=intsInPolytopeGivenIneqAndPt(A,rightHandSide,v);
    P.printVectorList(l);

    //    return 0;//!!!!!!!!!!!!!!!!!!!!!!!!

    //    P.printVectorList(intsInPolytopeGivenIneq(A,rightHandSide));






    {

      IntegerVectorList l2;
      for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
	{
	  IntegerVector temp(i->size()+1);
	  temp[0]=1;
	  for(int j=0;j<i->size();j++)temp[j+1]=(*i)[j];
	  l2.push_back(temp);
	  /*	  IntegerVectorList k;
	  k.push_back(*i);
	  if(isPartOfAZBasis(k))
	    P.printVector(*i);
	  */
	}

      IntegerVectorList l3;
      for(IntegerVectorList::const_iterator i=l2.begin();i!=l2.end();i++)
	//if(((*i)[1]<=2)&&((*i)[2]<=2)&&((*i)[3]<=17)&&((*i)[4]<=8))
	  l3.push_back(*i);
      //	if(((*i)[0]<=2)&&((*i)[1]<=2)&&((*i)[2]<=17)&&((*i)[3]<=8))l3.push_back(*i);

      IntegerMatrix M=rowsToIntegerMatrix(l3);

      P.printVectorList(M.getRows());
      fprintf(stdout,"size:%i\n",M.getHeight());

      set<set<int> > s;
      s.insert(set<int>());
      for(int i=0;i<5;i++)
	{
	  s=partsOfZBases(s,M);


      for(set<set<int> >::const_iterator i=s.begin();i!=s.end();i++)
	{
	  for(set<int>::const_iterator j=i->begin();j!=i->end();j++)
	    {
	      //	      fprintf(stderr,"%i ",*j);
	    }
	  //	  fprintf(stderr,"\n");
	}

      //	  fprintf(stderr,"\n\n\n\n");

	}
    }



    return 0;
  }
  const char *helpText()
  {
    return "This program computes the integer points in a polytope.\n";
  }
};

static IntsInPolytopeApplication theApplication;
