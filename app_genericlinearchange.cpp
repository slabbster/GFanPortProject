#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "buchberger.h"
#include "wallideal.h"
#include "termorder.h"
#include "gfanapplication.h"
#include "tropical2.h"

class GenericLinearChangeApplication : public GFanApplication
{
public:
  const char *helpText()
  {
    return "This program takes a list of polynomials and performs a generic linear change of coordinates by introducing nxn new variables.\n";
  }
  GenericLinearChangeApplication()
  {
    registerOptions();
  }

  const char *name()
  {
    return "_genericlinearchange";
  }

  Polynomial substitueMonomial(Monomial const &m, vector<Polynomial> const &substitutionTable, PolynomialRing const &R2)
  {
	  IntegerVector v=m.exponent;
	  Polynomial ret=R2.one();
	  for(int i=0;i<v.size();i++)
		  for(int j=0;j<v[i];j++)
			  ret*=substitutionTable[i];
	  return ret;
  }
  Polynomial substitute(Polynomial const &f, vector<Polynomial> const &substitutionTable, PolynomialRing const &R2)
  {
	  Polynomial ret(R2);
	  for(TermMap::const_iterator i=f.terms.begin();i!=f.terms.end();i++)
	  {
		  ret+=R2.polynomialFromField(i->second)*substitueMonomial(i->first,substitutionTable,R2);
	  }
	  return ret;
  }

  int main()
  {
    FileParser P(Stdin);

    PolynomialSet g=P.parsePolynomialSetWithRing();

    PolynomialRing R1=g.getRing();
    int n=g.getRing().getNumberOfVariables();
    vector<string> names=matrixVariableNames("g",n,n);
    for(int i=0;i<n;i++)names.push_back(R1.getVariableName(i));
    PolynomialRing R2(R1.getField(),names);
    vector<Polynomial> substitutionTable;
    for(int i=0;i<n;i++)
    {
    	Polynomial p(R2);
    	for(int j=0;j<n;j++)
    	{
    		p+=R2.ithVariable(n*n+j)*R2.ithVariable(i+j*n);
    	}
    	substitutionTable.push_back(p);
    }

    PolynomialSet g2(R2);
    for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
    	g2.push_back(substitute(*i,substitutionTable,R2));

    AsciiPrinter(Stdout).printPolynomialRing(R2);
    AsciiPrinter(Stdout).printPolynomialSet(g2);

    IntegerVectorList symmetries;
    {
    	IntegerVector v(n+n*n);

    for(int i=0;i<n;i++)
    {
    	v[n*n+i]=(n*n+((i+1)%n));
    	for(int j=0;j<n;j++)
			v[j*n+i]=(((j+1)%n)*n+((i)%n));
    }
    symmetries.push_back(v);
    }
    	{
    	IntegerVector v(n+n*n);

    	for(int i=0;i<n;i++)
    	{
    		int offset=0;
    		if(i==0)offset=+1;
    		if(i==1)offset=-1;

    		v[n*n+i]=(n*n+((i+offset)%n));
    	for(int j=0;j<n;j++)
			v[i*n+j]=(((i+offset)%n)*n+((j)%n));
    	}
        symmetries.push_back(v);
    }
    AsciiPrinter(Stdout)<<symmetries;

    return 0;
  }
};

static GenericLinearChangeApplication theApplication;
