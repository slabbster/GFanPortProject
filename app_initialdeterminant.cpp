#include  <iostream>
#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "buchberger.h"
#include "wallideal.h"
#include "termorder.h"
#include "gfanapplication.h"
#include "determinantpoly.h"
#include "log.h"

using namespace std;


class InitialDeterminantApplication : public GFanApplication
{
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  SimpleOption optionTakeDerivatives;
  const char *helpText()
  {
    return "This program \n";
  }
  InitialDeterminantApplication():
    optionTakeDerivatives("-j","Take derivatives of input vector.")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_initialdeterminant";
  }


  int main()
  {
    FileParser P(Stdin);

    PolynomialSet g=P.parsePolynomialSetWithRing();

    if(optionTakeDerivatives.getValue())
    {
    	PolynomialSet G=P.parsePolynomialSet(g.getRing());
        G.markAndScale(LexicographicTermOrder());
        PolynomialSet::iterator I=g.begin();
        for(PolynomialSet::const_iterator i=G.begin();i!=G.end();i++,I++)
        {
        	*I=Term(i->getMarked().c,Monomial(g.getRing(),-(i->getMarked().m.exponent)))**I;
        }
        debug<<g;
    }
    IntegerVector w=P.parseIntegerVector();
    int degree=P.parseInt();



//    PolyMatrix m(g,w,optionTakeDerivatives.getValue());
//    m.print();
/*
    int n=g.size();
    for(int i=0;i<n;i++)
    {
    	for(int j=0;j<n;j++)
    		cerr<<!m.data[i][j].isZero;
    	cerr<<endl;
    }
*/

/*    list<int> rows=interval(g.size(),!optionTakeDerivatives.getValue());
    list<int> columns=interval(g.size(),!optionTakeDerivatives.getValue());
    cerr<<rows.size();
    cerr<<columns.size();*/
//    pout<<m.determinantForm(rows,columns,degree,0)<<"\n";
    pout<<g.getRing()<<"{\n";
//    pout<<m.determinantForm(rows,columns,degree,0)/*.numberOfTerms()*/
//    pout<<m.initialFormOfDeterminant()<<"\n}\n";
    pout<<initialFormOfDeterminant(g,w,optionTakeDerivatives.getValue())<<"\n}\n";
    return 0;
  }
};

static InitialDeterminantApplication theApplication;







#if 0


#include  <iostream>
#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "buchberger.h"
#include "wallideal.h"
#include "termorder.h"
#include "gfanapplication.h"
#include "tropical2.h"
#include "tropicaldeterminant.h"
#include "log.h"

using namespace std;



class PolyMatrix
{
public:
	class PolyMatrixEntry
	{
	public:
		bool isZero;
		Polynomial p;
		int maxDegree,minDegree;
		vector<Polynomial> forms;
		PolyMatrixEntry(Polynomial const &p_, IntegerVector const &w):
			p(p_),
			isZero(p_.isZero())
			{
			if(p.isZero())
			{
				maxDegree=-1;
				minDegree=0;
				return;
			}
			maxDegree=p_.degree(w);
			minDegree=-p_.degree(-w);


	//		debug<<maxDegree<<minDegree<<"\n";
			for(int i=0;i<=maxDegree;i++)forms.push_back(Polynomial(p.getRing()));
			for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
				{
					forms[dot(i->first.exponent,w)]+=Term(i->second,i->first);
				}
		}
	};

	list<int>::iterator indexOfSparsestRow(list<int> &rows, list<int> const &columns)
	{
		list<int>::iterator ret=rows.end();
		int bestNumberOfZeros=-1;
		for(list<int>::iterator i=rows.begin();i!=rows.end();i++)
		{
			int numberOfZeros=0;
			for(list<int>::const_iterator j=columns.begin();j!=columns.end();j++)numberOfZeros+=data[*i][*j].isZero;
			if(numberOfZeros>bestNumberOfZeros)
			{
				bestNumberOfZeros=numberOfZeros;
				ret=i;
			}
		}
//		cerr<<":"<<bestNumberOfZeros<<endl;
		return ret;
	}
	vector<vector<PolyMatrixEntry> > data;
	PolynomialRing theRing;
	//constructor for Jacobi matrix
	PolyMatrix(PolynomialSet const &generators, IntegerVector const &w):
		theRing(generators.getRing())
	{
		int n=generators.size();
		PolynomialSet::const_iterator I=generators.begin();
		for(int i=0;i<n;i++,I++)
	{
			vector<PolyMatrixEntry> row;
			for(int j=0;j<n;j++)
			{
	//			debug<<"P"<<*I<<"\n";
				Polynomial f=I->derivative(j);
	//			debug<<f<<"\n";
				row.push_back(PolyMatrixEntry(f,w));
			}
			data.push_back(row);
		}
	}
	void print()
	{
		int n=data.size();
		for(int i=0;i<n;i++)
		{
			for(int j=0;j<n;j++)
			{
				debug<<data[i][j].p<<",";
			}
			debug<<"\n";
		}
	}
	void p(list<int> l)
	{
		cerr<<"{";
		for(list<int>::const_iterator i=l.begin();i!=l.end();i++)
			cerr<<*i<<",";
		cerr<<"}\n";
	}
	Polynomial determinantForm(list<int> rows, list<int> columns, int degree, int level)
	{
		Polynomial ret(theRing);
		if(rows.size()==1)return data[rows.front()][columns.front()].p;
		list<int>::iterator chosenRowIterator=indexOfSparsestRow(rows, columns);//rows.begin();
		int chosenrow=*chosenRowIterator;
		{
			list<int>::iterator temp=chosenRowIterator;temp++;
			rows.erase(chosenRowIterator);
			//no need to update rows afterwards since it is stored on the stack
		}
		if(level==0)cerr<<"-"<<chosenrow<<endl;
		if(level==1)cerr<<"++"<<chosenrow<<endl;
		if(level==2)cerr<<"***"<<chosenrow<<endl;
		if(level==3)cerr<<"%%%%%%%%"<<chosenrow<<endl;
		int sign=1;
		for(list<int>::iterator i=columns.begin();i!=columns.end();i++)
		{
			int chosencol=*i;
			if(!data[chosenrow][chosencol].isZero)
			{
				list<int>::iterator temp=i;temp++;
				columns.erase(i);
				i=temp;
				ret+=data[chosenrow][chosencol].p*determinantForm(rows,columns,degree,level+1)*theRing.getField().zHomomorphism(sign);
				columns.insert(i,chosencol);
				i--;
			}
			sign*=-1;
		}
		rows.push_front(chosenrow);
		return ret;
	}
};


class InitialDeterminantApplication : public GFanApplication
{
public:
  const char *helpText()
  {
    return "This program \n";
  }
  InitialDeterminantApplication()
  {
    registerOptions();
  }

  char *name()
  {
    return "_initialdeterminant";
  }

  list<int> interval(int n)
  {
	  list<int> ret;
	  for(int i=0;i<n;i++)ret.push_back(i);
	  return ret;
  }

  int main()
  {
    FileParser P(Stdin);

    PolynomialSet g=P.parsePolynomialSetWithRing();
    IntegerVector w=P.parseIntegerVector();
    int degree=P.parseInt();

    PolyMatrix m(g,w);
//    m.print();

    int n=g.size();
    for(int i=0;i<n;i++)
    {
    	for(int j=0;j<n;j++)
    		cerr<<!m.data[i][j].isZero;
    	cerr<<endl;
    }


    list<int> rows=interval(g.size());
    list<int> columns=interval(g.size());
    pout<<m.determinantForm(rows,columns,0,0)<<"\n";

    return 0;
  }
};

static InitialDeterminantApplication theApplication;
#endif
