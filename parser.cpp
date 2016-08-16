#include "parser.h"
#include <stdio.h>
#include <assert.h>
#include "printer.h"
#include "field_rationals.h"
#include "field_zmodpz.h"
#include <vector>
#include <string>
#include "log.h"

using namespace std;


static string variableIndexToString(int i)
{
  //  fprintf(Stderr,"I:%i\n",i);
  static char s[4];
  if(i>=0 && i<26)sprintf(s,"%c",i+'a');
  else if(i>=26 && i<52)sprintf(s,"%c",i+'A'-26);
  else assert(0);
  return string(s);
}


//--------------------------------------------------
// Parser
//--------------------------------------------------

void Parser::parserError(const char *expected, char c)
{
  fprintf(Stderr,"PARSER ERROR: Expected: \"%s\"    Found: \"%c\"\n",expected,c);
}


//--------------------------------------------------
// CharacterBasedParser
//--------------------------------------------------



PolynomialRing CharacterBasedParser::azAZ(Field const &f, int n)
{
  assert(n<=52);
  vector<string> na(n);
  for(int i=0;i<n;i++)
    {
      na[i]=variableIndexToString(i);
    }
  return PolynomialRing(f,na);
}

int CharacterBasedParser::nextNonBlank()
{
 int c=getChar();

 while(c!=0 && (c=='\t' || c==' ' || c=='\n' || c==13 || c==10))c=getChar();

 return c;
}


int CharacterBasedParser::nextNonBlankDoNotGet()
{
  int k=nextNonBlank();
  ungetChar(k);
  return k;
}



bool CharacterBasedParser::isLeftBracket(int c)
{
  return c=='(' || c=='[' || c=='{';
}


bool CharacterBasedParser::isRightBracket(int c)
{
  return c==')' || c==']' || c=='}';
}


int CharacterBasedParser::parseChar()
{
  return getChar();
}


FieldElement CharacterBasedParser::parseFieldElementFromInteger(Field const &f)
{
  FieldElement ret=f.zHomomorphism(0);
  int c=getChar();

  while(c!=0 && (c=='\t' || c==' ' || c=='\n' || c==10 || c==13))c=getChar();
  if(c=='-')return f.zHomomorphism(-1)*parseFieldElementFromInteger(f);
  while(c!=0)
    {
      if(c>='0' && c<='9')
	{
	  ret=ret*f.zHomomorphism(10)+f.zHomomorphism(c-'0');
	  c=getChar();
	}
      else
	break;
    }
  ungetChar(c);

  return ret;
}

int CharacterBasedParser::parseInt()
{
  int ret=0;
  int c=getChar();

  int numberOfDigits=0;

  while(c!=0 && (c=='\t' || c==' ' || c=='\n' || c==10 || c==13))c=getChar();
  if(c=='-')return -parseInt();
  if(c=='(')//allow parentheses for Laurent exponents
    {
      ungetChar(c);
      IntegerVector v=parseIntegerVector();
      assert(v.size()==1);
      return v[0];
    }
  while(c!=0)
    {
      if(c>='0' && c<='9')
	{
	  numberOfDigits++;
	  assert(numberOfDigits<10);
	  ret=ret*10+(c-'0');
	  c=getChar();
	}
      else
	break;
    }
  ungetChar(c);

  return ret;
}


double CharacterBasedParser::parseFloat()
{
  double ret=0;
  static char s[64];
  int i=0;

  int c=nextNonBlank();
  while(i<64 && ((c>='0' && c<='9') || c=='.' || c=='-') )
    {
      s[i++]=c;
      c=getChar();
    }
  assert(i<64);
  s[i]=0;

  ungetChar(c);

  sscanf(s,"%lf",&ret); //does this work? Are we reading floats or doubles?
  return ret;
}


ComplexNumber CharacterBasedParser::parseComplexNumber()
{
  FloatVector v=parseFloatVector();
  assert(v.size()==2);

  return ComplexNumber(v[0],v[1]);
}


FieldElement CharacterBasedParser::parseFieldElement(Field const &f)
{
  //  return FieldElement(parseInt());
  FieldElement a=parseFieldElementFromInteger(f);
  int c=nextNonBlank();
  FieldElement b=f.zHomomorphism(1);
  if(c=='/')
    b=parseFieldElementFromInteger(f);
  else
    ungetChar(c);

  //  FieldElement A=Field::staticZHomomorphism(a);
  //  FieldElement B=Field::staticZHomomorphism(b);

  a*=b.inverse();

  return a;
}


bool CharacterBasedParser::isVariable(int c)
{
  return (c>='a' && c<='z') || (c>='A' && c<='Z');
}


bool CharacterBasedParser::isDigit(int c)
{
  return (c>='0' && c<='9');
}


int CharacterBasedParser::variableIndex(int c)
{
  if(c>='a' && c<='z')return c-'a';
  if(c>='A' && c<='Z')return c-'A'+26;
  assert(0);
  return 0;
}

static bool isVariableCharacter(int c)
{
  return ((c>='a' && c<='z')||(c>='A' && c<='Z')||(c>='0' && c<='9')||(c=='_' || c==']' || c=='['));
}

Monomial CharacterBasedParser::parseMonomial(PolynomialRing const &r)
{
  //  fprintf(Stderr,"Atest\n");
  int dimension=r.getNumberOfVariables();
  //  fprintf(Stderr,"Dim%i\n",dimension);
  /* Parses a monomial without spaces and without sign */

  //  Monomial m(IntegerVector((dimension==-1)?0:dimension));
  Monomial m(r,IntegerVector(dimension));
  IntegerVector &ret=m.exponent;

  int c=nextNonBlank();
  //  int c=getChar();

  if(c=='1')return m;

  while(c!=0 && (c=='\t' || c==' ' || c=='\n' || c==10 || c==13))c=getChar();
  while(c!=0)
    {
      string name;
      while(isVariableCharacter(c))
	{
	  name+=c;
	  if(r.variableIndex(name)!=-1)break;
	  c=getChar();
	}

      if(r.variableIndex(name)!=-1)
	{
	  //fprintf(Stderr,"NAME:%s\n",name.c_str());
	  int index=r.variableIndex(name);
	  int exponent=1;
	  if(dimension==-1)
	    ret.grow(index+1);
	  else
	    assert(index<dimension);
	  //	  c=getChar();
	  c=nextNonBlank();
	  if(isDigit(c)) // This makes parsing of Singular polynomials possible
	    {
	      ungetChar(c);
	      c='^';
	    }
	  if(c=='^')
	    {
	      exponent=parseInt();
	      //c=getChar();
	      c=nextNonBlank();
	    }
          //          printf("Index:%i, Exponent:%i\n",index,exponent);
	  ret[index]+=exponent;
	}
      else if(name.size())
	{
	  fprintf(Stderr,"Unknown variable:%s\n",name.c_str());
	  assert(0);
	}
      else if(c=='*')
	{
	  //	  c=getChar();
	  c=nextNonBlank();
	}
      else break;
    }
  ungetChar(c);


  //  fprintf(Stderr,"test\n");

  return m;
}


Term CharacterBasedParser::parseTerm(PolynomialRing const &r)
{
  //  fprintf(Stderr,"testC\n");
  //  fprintf(Stderr,"parseTerm - Begin\n");
  int c=nextNonBlank();
  if(c=='-')
    {
      c=nextNonBlank();
      ungetChar(c);
      FieldElement k(r.getField());
      k=r.getField().zHomomorphism(1);
  //      k=k.one();
      if(isDigit(c))
	{
	  k=parseFieldElement(r.getField());
	}
      Monomial m=parseMonomial(r);
      //  fprintf(Stderr,"parseTerm - End\n");
      return Term(-k,m);
    }
  else
    {
      ungetChar(c);
      FieldElement k(r.getField());
      //    fprintf(Stderr,"A\n");
      k=r.getField().zHomomorphism(1);
  //      k=k.one();
  //fprintf(Stderr,"B\n");
      if(isDigit(c))
	{
	  k=parseFieldElement(r.getField());
	}
      Monomial m=parseMonomial(r);
      //   fprintf(Stderr,"parseTerm - End\n");
      return Term(k,m);
    }

  FieldElement k=parseFieldElement(r.getField());
  Monomial m=parseMonomial(r);

  // fprintf(Stderr,"parseTerm :");
  // AsciiPrinter(Stderr).printFieldElement(c);
  // AsciiPrinter(Stderr).printMonomial(m);

  //  fprintf(Stderr,"parseTerm - End\n");
  return Term(k,m);
}


Polynomial CharacterBasedParser::parsePolynomial(PolynomialRing const &r)
{
  //fprintf(Stderr,"%i\n",r.getNumberOfVariables());
  //  fprintf(Stderr,"parsePolynomial - Begin\n");
  bool first=true;
  Polynomial p(r);
  Monomial firstMonomial(r);

  // fprintf(Stderr,"parsePolynomialBegin\n");

  while(true)
    {
      //  fprintf(Stderr,"test1\n");
      int c=nextNonBlank();

      //  fprintf(Stderr,"test2\n");
      ungetChar(c);

      if(isDigit(c)|| c=='+' || c=='-' || isVariable(c))
        {
	  //  fprintf(Stderr,"test3\n");
          if(c=='+')nextNonBlank();
	  Term t=parseTerm(r);
	  if(!t.c.isZero())
	  {
	    if(first)
	      {
		firstMonomial=t.m;
		first=false;
	      }
	    //  fprintf(Stderr,"test4\n");
	    Polynomial q=Polynomial(t);
	    // AsciiPrinter(Stderr).printPolynomial(q);
	    //if(q.getNumberOfVariables()>p.getNumberOfVariables())p.changeNumberOfVariables(q.getNumberOfVariables());
	    //if(p.getNumberOfVariables()>q.getNumberOfVariables())q.changeNumberOfVariables(p.getNumberOfVariables());
	    p+=q;
	    //  fprintf(Stderr,"test5\n");
	  }
        }
      else
        break;
    }
  // fprintf(Stderr,"parsePolynomialEnd\n");
  // AsciiPrinter(Stderr).printPolynomial(p);
  // fprintf(Stderr,"parsePolynomialEnd\n");

  if(!first)
    {
      firstMonomial.exponent.resize(p.getNumberOfVariables());
      p.mark(firstMonomial);
      //AsciiPrinter(Stderr).printVector(firstMonomial.exponent);
    }

  //fprintf(Stderr,"parsePolynomial - End\n");
  return p;
}


PolynomialSet CharacterBasedParser::parsePolynomialSet(PolynomialRing const &r)
{
  //  fprintf(Stderr,"%i\n",r.getNumberOfVariables());
  PolynomialSet ret(r);
  int c=nextNonBlank();

  if(!isLeftBracket(c))
    {
      parserError("left bracket",c);
      assert(0);
    }

  c=nextNonBlank();
  if(isRightBracket(c))return ret;
  ungetChar(c);
  do
    {
      ret.push_back(parsePolynomial(r));
    }
  while((c=nextNonBlank())==',');
  assert(isRightBracket(c));

  return ret;
}


PolynomialSetList CharacterBasedParser::parsePolynomialSetList(PolynomialRing const &r)
{
  PolynomialSetList ret;
  int c=nextNonBlank();

  assert(isLeftBracket(c));
  do
    {
      ret.push_back(parsePolynomialSet(r));
    }
  while((c=nextNonBlank())==',');
  assert(isRightBracket(c));

  return ret;
}


IntegerVector CharacterBasedParser::parseIntegerVector()
{
  list<int> temp;
  int c=nextNonBlank();

  assert(isLeftBracket(c));
  do
    temp.push_back(parseInt());
  while((c=nextNonBlank())==',');
  assert(isRightBracket(c));

  IntegerVector ret(temp.size());

  {
    int j=0;
    for(list<int>::iterator i=temp.begin();i!=temp.end();i++)
      ret[j++]=*i;
  }

  return ret;
}


FloatVector CharacterBasedParser::parseFloatVector()
{
  list<double> temp;
  int c=nextNonBlank();

  assert(isLeftBracket(c));
  do
    temp.push_back(parseFloat());
  while((c=nextNonBlank())==',');
  assert(isRightBracket(c));

  FloatVector ret(temp.size());

  {
    int j=0;
    for(list<double>::iterator i=temp.begin();i!=temp.end();i++)
      ret[j++]=*i;
  }
  return ret;
}


ComplexVector CharacterBasedParser::parseComplexVector()
{
  list<ComplexNumber> temp;
  int c=nextNonBlank();

  assert(isLeftBracket(c));
  do
    temp.push_back(parseComplexNumber());
  while((c=nextNonBlank())==',');
  assert(isRightBracket(c));

  ComplexVector ret(temp.size());

  {
    int j=0;
    for(list<ComplexNumber>::iterator i=temp.begin();i!=temp.end();i++)
      ret[j++]=*i;
  }
  return ret;
}


IntegerVectorList CharacterBasedParser::parseIntegerVectorList()
{
  IntegerVectorList ret;
  int c=nextNonBlank();

  assert(isLeftBracket(c));
  if((c=nextNonBlank())=='}')return ret;
  ungetChar(c);

  do
    ret.push_back(parseIntegerVector());
  while((c=nextNonBlank())==',');
  assert(isRightBracket(c));

  return ret;
}


IntegerVectorList CharacterBasedParser::parseIntegerVectorList4ti2()
{
  IntegerVectorList ret;

  int m=parseInt();
  int n=parseInt();
  for(int i=0;i<m;i++)
    {
      IntegerVector v(n);
      for(int j=0;j<n;j++)
	{
	  v[j]=parseInt();
	}
      ret.push_back(v);
    }
  return ret;
}


Field CharacterBasedParser::parseField()
{
  int c=nextNonBlank();

  if(c=='Q')
    {
      return Q;
    }

  if(c=='Z')
    {
      c=nextNonBlank();
      assert(c=='/');
      int p=parseInt();
      c=nextNonBlank();
      assert(c=='Z');

      return FieldZModPZ(p);
    }

  parserError("field",c);
  assert(0);
}


string CharacterBasedParser::parseVariableName()
{
  string name;
  int c=nextNonBlank();
  while((c!=',') && (c!=']') && (c!=' ') && (c!='\t') && (c!='\n') && (c!=10) && (c!=13))
    {
      name+=c;
      c=getChar();
    }
  ungetChar(c);
  return name;
}


vector<string> CharacterBasedParser::parseVariableList()
{
  vector<string> ret;

  int c=nextNonBlank();
  assert(c=='[');
  c=nextNonBlank();
  if(c==']') return ret;
  ungetChar(c);
  do
    {
      ret.push_back(parseVariableName());
      c=nextNonBlank();
    }
  while((c)==',');
  assert(c==']');
  return ret;
}


PolynomialRing CharacterBasedParser::parsePolynomialRing()
{
  Field f=parseField();
  vector<string> variables=parseVariableList();

  return PolynomialRing(f,variables);
}


Polynomial CharacterBasedParser::parsePolynomialWithRing()
{
  PolynomialRing theRing=parsePolynomialRing();
  return parsePolynomial(theRing);
}


PolynomialSet CharacterBasedParser::parsePolynomialSetWithRing()
{
  //  PolynomialSet ret(r);
  int c=nextNonBlankDoNotGet();

  if(!isLeftBracket(c))
    {
      PolynomialRing theRing=parsePolynomialRing();
      return parsePolynomialSet(theRing);
    }

  PolynomialRing r=azAZ(Q);
  // fprintf(Stderr,"%i\n",r.getNumberOfVariables());
  log3  AsciiPrinter(Stderr).printPolynomialRing(r);
  //fprintf(Stderr,"a\n");
  PolynomialSet temp=parsePolynomialSet(r);
  //fprintf(Stderr,"a\n");
  int n=temp.numberOfVariablesInUseConsecutive();
  //fprintf(Stderr,"a\n");
  PolynomialRing r2=azAZ(Q,n);
  //fprintf(Stderr,"a\n");

  return temp.embeddedInto(r2);
}


PolynomialSetList CharacterBasedParser::parsePolynomialSetListWithRing()
{
  //  PolynomialSet ret(r);
  int c=nextNonBlankDoNotGet();

  if(!isLeftBracket(c))
    {
      PolynomialRing theRing=parsePolynomialRing();
      return parsePolynomialSetList(theRing);
    }

  PolynomialRing r=azAZ(Q);
  PolynomialSetList temp=parsePolynomialSetList(r);
  assert(temp.size());
  int n=temp.begin()->numberOfVariablesInUseConsecutive();//!!!!!!!
  PolynomialRing r2=azAZ(Q,n);

  PolynomialSetList ret;
  for(PolynomialSetList::const_iterator i=temp.begin();i!=temp.end();i++)
    ret.push_back(i->embeddedInto(r2));

  return ret;
}

//--------------------------------------------------
// FileParser
//--------------------------------------------------

int FileParser::getChar()
{
  int c=getc(f);
  //  fprintf(Stderr,"getChar c=\'%c\',%i\n",c,c);
  if(c==EOF)return 0;

  return c;
}


void FileParser::ungetChar(int c)
{
  ungetc(c,f);
}


FileParser::FileParser(FILE *f)
{
  this->f=f;
}


//--------------------------------------------------
// StringParser
//--------------------------------------------------

int StringParser::getChar()
{
  if(hasUngotten)
    {
      hasUngotten=false;
      return ungotten;
    }
  return s[index++];
}

void StringParser::ungetChar(int c)
{
  assert(hasUngotten==false);
  hasUngotten=true;
  ungotten=c;
}

StringParser::StringParser(const char *s):
  index(0),
  hasUngotten(false)
{
  this->s=s;
}
