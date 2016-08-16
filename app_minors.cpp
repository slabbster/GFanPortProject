#include "vektor.h"
#include "printer.h"
#include "parser.h"
#include "gfanapplication.h"
#include "minors.h"
#include "field_rationals.h"
#include <iostream>

class MinorsApplication : public GFanApplication
{
  IntegerOption rOption;
  IntegerOption dOption;
  IntegerOption nOption;
  SimpleOption M2Option;
  SimpleOption optionWithNames;
  SimpleOption dressianOption;
  SimpleOption pluckerSymmetriesOption;
  SimpleOption grassmannNormalizOption;
  SimpleOption symmetryOption;
  SimpleOption parametrizeOption;
  //  SimpleOption threeTermOption;
public:
  const char *helpText()
  {
    return "This program will generate the r*r minors of a d*n matrix of indeterminates.\n";
  }
  MinorsApplication():
    rOption("-r","Specify r.",1),
    dOption("-d","Specify d.",1),
    nOption("-n","Specify n.",1),
    M2Option("-M2","Use Macaulay2 conventions for order of variables."),
    dressianOption("--dressian","Produce tropical defining the Dressian(3,n) instead. (The signs may not be correct, that is the equations may not be Pluecker relations.)"),
    optionWithNames("--names","Assign names to the minors."),
    pluckerSymmetriesOption("--pluckersymmetries","Do nothing but produce symmetry generators for the Pluecker ideal."),
    grassmannNormalizOption("--grassmannnormalize","Produce polynomials describing the tropical polynomial map taking a plucker vector to a plucker vector where the leaf edges have length one."),
    symmetryOption("--symmetry","Produces a list of generators for the group of symmetries keeping the set of minors fixed. (Only without --names)."),
    parametrizeOption("--parametrize","Parametrize the set of d times n matrices of Barvinok rank less than or equal to r-1 by a list of tropical polynomials.")
/*,
															threeTermOption("--threeTerm","Do nothing but output the three term Plucker relations.")*/
  {
    registerOptions();
    grassmannNormalizOption.hide();
  }
  const char *name()
  {
    return "_minors";
  }

  IntegerVectorList symmetries(int r, int d, int n)
  {
    IntegerVectorList ret;
    {
      IntegerVector v1(d*n);
      IntegerVector v2(d*n);
      IntegerVector v3(d*n);
      IntegerVector v4(d*n);
      for(unsigned int i=0;i<d;i++)
	for(unsigned int j=0;j<n;j++)
	  {
	    v1[i*n+j]=((i+1)%d)*n+j;
	    v2[i*n+j]=(i-(i==1)+(i==0))*n+j;
	    v3[i*n+j]=i*n+(j+1)%n;
	    v4[i*n+j]=i*n+(j-(j==1)+(j==0));
	  }
      ret.push_back(v1);
      ret.push_back(v2);
      ret.push_back(v3);
      ret.push_back(v4);
    }
    return ret;
  }


static int lookup(vector<int> const &v, int i)
{
  for(int j=0;j<v.size();j++)
    {
      if(v[j])
	{
	  if(i==0)
	    return j;
	  i--;
	}
    }
  return 0;
}

  PolynomialSet normalizingMap(int d, int r, int n, bool M2)
  {
    assert(d==2);
    assert(r==2);

    vector<string> pnames=subsetVariableNames("p",n,d,M2);
    vector<string> qnames=subsetVariableNames("q",n,d,M2);
    vector<string> names(pnames.size()*2);
    for(int i=0;i<pnames.size();i++)names[i]=pnames[i];
    for(int i=0;i<qnames.size();i++)names[i+pnames.size()]=qnames[i];

    PolynomialRing R(Q,names);

    PolynomialSet ret(R);

    IntegerVector w(names.size());
    for(int i=0;i<pnames.size();i++)w[i]=1;

    vector<int> I;
    for(int i=0;i<n-d;i++)I.push_back(0);
    for(int i=0;i<d;i++)I.push_back(1);
    do
      {
	int k=lookup(I,0);
	int l=lookup(I,1);
	set<int> setkl;
	setkl.insert(k);
	setkl.insert(l);
	Polynomial f(R);
	Polynomial g(R);
	vector<int> J;
	for(int i=0;i<n-d;i++)J.push_back(0);
	for(int i=0;i<d;i++)J.push_back(1);
	do
	  {
	    int i=lookup(J,0);
	    int j=lookup(J,1);

	    set<int> setij;
	    setij.insert(i);
	    setij.insert(j);
	    set<int> setik;
	    setik.insert(i);
	    setik.insert(k);
	    set<int> setjk;
	    setjk.insert(j);
	    setjk.insert(k);
	    set<int> setil;
	    setil.insert(i);
	    setil.insert(l);
	    set<int> setjl;
	    setjl.insert(j);
	    setjl.insert(l);
	    if((J[l]==0))
	      {
		IntegerVector v=w;
		v[subsetToVariableIndex(setil,n,d,M2)]+=-1;
		v[subsetToVariableIndex(setjl,n,d,M2)]+=-1;
		v[subsetToVariableIndex(setij,n,d,M2)]+=1;
		v[subsetToVariableIndex(setkl,n,d,M2)]+=1;
		f+=Polynomial(Term(R.getField().zHomomorphism(1),Monomial(R,v)));
	      }
	    if((J[k]==0))
	      {
		IntegerVector v=w;
		v[subsetToVariableIndex(setik,n,d,M2)]+=-1;
		v[subsetToVariableIndex(setjk,n,d,M2)]+=-1;
		v[subsetToVariableIndex(setij,n,d,M2)]+=1;
		v[subsetToVariableIndex(setkl,n,d,M2)]+=1;
		g+=Polynomial(Term(R.getField().zHomomorphism(1),Monomial(R,v)));
	      }
	  }
	while(next_permutation(J.begin(),J.end()));

	IntegerVector v(names.size());
	v[pnames.size()+subsetToVariableIndex(setkl,n,d,M2)]=1;

	ret.push_back(Polynomial(Term(R.getField().zHomomorphism(-1),Monomial(R,v)))+f*g);
      }
    while(next_permutation(I.begin(),I.end()));

    return ret;
  }

  int main()
  {
    FileParser P(Stdin);

    int d=dOption.getValue();
    int n=nOption.getValue();
    int r=rOption.getValue();
    bool M2=M2Option.getValue();

    assert(r<=d);
    assert(r<=n);

    PolynomialRing R(Q,matrixVariableNames("m",d,n));

    if(parametrizeOption.getValue())
      {
	vector<string> A=matrixVariableNames("a",d,r-1);
	vector<string> B=matrixVariableNames("b",r-1,n);
	A.insert(A.end(),B.begin(),B.end());
	PolynomialRing R2(Q,A);
	PolynomialSet s(R2);
	for(int i=0;i<d;i++)
	  for(int j=0;j<n;j++)
	    {
	      Polynomial p(R2);
	      for(int k=0;k<r-1;k++)
		{
		  IntegerVector e(A.size());
		  e[(r-1)*i+k]=1;
		  e[(n)*k+j+(d*(r-1))]=1;
		  p+=Term(R2.getField().zHomomorphism(1),Monomial(R2,e));
		}
	      s.push_back(p);
	    }
	AsciiPrinter(Stdout).printPolynomialRing(s.getRing());
	AsciiPrinter(Stdout).printPolynomialSet(s);
      }
    else
      if(grassmannNormalizOption.getValue())
      {
	PolynomialSet s = normalizingMap(d,r,n,M2);
	AsciiPrinter(Stdout).printPolynomialRing(s.getRing());
	AsciiPrinter(Stdout).printPolynomialSet(s);
      }
    else if(dressianOption.getValue())
      {
	//	int d=3;
	PolynomialRing R(Q,subsetVariableNames("p",n,d,M2));
	PolynomialSet g(R);
	vector<int> I;
	if(n-d-2>=0)if(d-2>=0)
	for(int i=0;i<n-d-2;i++)I.push_back(0);
	for(int i=0;i<4;i++)I.push_back(1);
	for(int i=0;i<d-2;i++)I.push_back(2);

	do
	  {
	    vector<int> ijkl;
	    set<int> S;
	    int m=32873;
	    for(int i=0;i<I.size();i++)
	      {
		if(I[i]==1)
		  ijkl.push_back(i);
		if(I[i]==2)
		  {
		    S.insert(i);
		    m=i;
		  }
	      }
	    IntegerVector pA(R.getNumberOfVariables());
	    IntegerVector pB(R.getNumberOfVariables());
	    IntegerVector pC(R.getNumberOfVariables());
	    {
	      int i=ijkl[0];
	      int j=ijkl[1];
	      int k=ijkl[2];
	      int l=ijkl[3];
	      fprintf(Stderr,"(%i,%i,%i,%i)%i\n",i,j,k,l,m);
	      {
		set<int> s=S;s.insert(i);s.insert(j);//s.insert(m);
		pA[subsetToVariableIndex(s,n,d,M2)]=1;
	      }
	      {
		set<int> s=S;s.insert(k);s.insert(l);//s.insert(m);
		pA[subsetToVariableIndex(s,n,d,M2)]=1;
	      }
	      {
		set<int> s=S;s.insert(i);s.insert(k);//s.insert(m);
		pB[subsetToVariableIndex(s,n,d,M2)]=1;
	      }
	      {
		set<int> s=S;s.insert(j);s.insert(l);//s.insert(m);
		pB[subsetToVariableIndex(s,n,d,M2)]=1;
	      }
	      {
		set<int> s=S;s.insert(i);s.insert(l);//s.insert(m);
		pC[subsetToVariableIndex(s,n,d,M2)]=1;
	      }
	      {
		set<int> s=S;s.insert(j);s.insert(k);//s.insert(m);
		pC[subsetToVariableIndex(s,n,d,M2)]=1;
	      }

	      Polynomial p=
		Polynomial(Term(R.getField().zHomomorphism(1),Monomial(R,pA)))
		-Polynomial(Term(R.getField().zHomomorphism(1),Monomial(R,pB)))
		+Polynomial(Term(R.getField().zHomomorphism(1),Monomial(R,pC)));
	      g.push_back(p);
	    }
	  }
	while(next_permutation(I.begin(),I.end()));
	AsciiPrinter(Stdout).printPolynomialRing(R);
	AsciiPrinter(Stdout).printNewLine();
	AsciiPrinter(stdout).printPolynomialSet(g);
      }
    else if(pluckerSymmetriesOption.getValue())
      {
	int N=subsetVariableNames("p",n,d,M2).size();
	IntegerVectorList permutations1;
	{
	  IntegerVector p1(n);
	  IntegerVector p2(n);
	  for(int i=0;i<n;i++)p1[i]=(i+1)%n;
	  for(int i=0;i<n;i++)p2[i]=i;
	  if(p2.size()>2)
	    {
	      p2[0]=1;
	      p2[1]=0;
	    }
	  permutations1.push_back(p1);
	  permutations1.push_back(p2);
	}
	IntegerVectorList permutations2;
	IntegerVectorList signs2;
	for(IntegerVectorList::const_iterator k=permutations1.begin();k!=permutations1.end();k++)
	  {
	    IntegerVector p(N);
	    IntegerVector signs(N);
	    vector<int> I;
	    for(int i=0;i<n-r;i++)I.push_back(0);
	    for(int i=0;i<r;i++)I.push_back(1);
	    do
	      {
		set<int> indexSet1;
		set<int> indexSet2;
		list<int> indexSet2List;
		for(int i=0;i<n;i++)if(I[i])
		  {
		    indexSet1.insert(i);
		    indexSet2.insert((*k)[i]);
		    indexSet2List.push_back((*k)[i]);
		  }
		int permSign=1;
		for(list<int>::const_iterator i=indexSet2List.begin();i!=indexSet2List.end();i++)
		  for(list<int>::const_iterator j=indexSet2List.begin();j!=i;j++)
		    if(*i<*j)permSign*=-1;



		p[subsetToVariableIndex(indexSet1,n,r,M2)]=subsetToVariableIndex(indexSet2,n,r,M2);
		signs[subsetToVariableIndex(indexSet2,n,r,M2)]=permSign;
	      }
	    while(next_permutation(I.begin(),I.end()));
	    signs2.push_back(signs);
	    permutations2.push_back(p);
	  }
	AsciiPrinter(Stdout).printVectorList(permutations2);
	AsciiPrinter(Stdout).printVectorList(signs2);
      }
    /*    else if(threeTermOption.getValue())
      {
	for(int i=0;i<n-d-2-4;i++)I.push_back(0);
	for(int i=0;i<d-2;i++)I.push_back(1);
	for(int i=0;i<4;i++)I.push_back(2);
	do
	  {
	    set<int> S;
	    vector<int> ijkl;
	    list<int> indexSet2List;
	    for(int i=0;i<n;i++)
	      {
		if(I[i]==1)S.insert(i);
		if(I[i]==1)ijkl.insert(i);
	      }
	  }
	while(next_permutation(I.begin(),I.end()));
	}*/
    else
      {

	if(optionWithNames.getValue())
	  {
	    assert(d==r);
	    vector<string> a=matrixVariableNames("m",d,n);
	    vector<string> b=subsetVariableNames("p",n,d,M2);
	    for(vector<string>::const_iterator i=b.begin();i!=b.end();i++)a.push_back(*i);
	    R=PolynomialRing(Q,a);
	  }

	PolynomialSet p=minors(R,r,d,n,optionWithNames.getValue(),M2);
	AsciiPrinter(Stdout).printPolynomialRing(p.getRing());
	AsciiPrinter(Stdout).printNewLine();
	AsciiPrinter(Stdout).printPolynomialSet(p);
      }

    if(symmetryOption.getValue())
      {
	AsciiPrinter(Stdout)<<symmetries(r,d,n);
      }
    return 0;
  }
};

static MinorsApplication theApplication;
