#include "gfanapplication.h"
#include <stdio.h>
#include <iostream>
#include "field.h"
#include "field_rationals.h"
#include "field_zmodpz.h"
#include "division.h"
#include "printer.h"
#include "symmetry.h"
#include "log.h"
#include "linalg.h"
#include "polymakefile.h"

GFanApplication::FieldOption::FieldOption():
  IntegerOption("--mod","Set the field to Z / pZ. The value must be a prime number p between 2 and 32749.\n")
{
}


void GFanApplication::FieldOption::onOptionsParsed()
{
  if(getValue())
    {
      fprintf(Stderr,"CANNOT CHANGE FIELD\n");
      assert(0);
      //      Field::setField(new FieldZModPZ(getValue()));
      //      Field::setField(Field::find("Zmod2Z"));
    }
}


GFanApplication::LogLevelOption::LogLevelOption():
  IntegerOption("--log","Set the logging level.\n")
{
  hide();
}


GFanApplication::XmlOption::XmlOption():
  SimpleOption("--xml","Switch to XML polyhedral data format output.\n")
  {
    hide();
  }


GFanApplication::StdinFileOption::StdinFileOption():
	StringOption("--stdin","Set a file to be read as stdin.\n","")
	{
	hide();
	}

void GFanApplication::StdinFileOption::onOptionsParsed()
{
	if(strlen(getValue())>0)freopen(getValue(),"r",stdin);
}

GFanApplication::StdoutFileOption::StdoutFileOption():
	StringOption("--stdout","Set a file to be read as stdin.\n","")
	{
	hide();
	}

void GFanApplication::StdoutFileOption::onOptionsParsed()
{
	if(strlen(getValue())>0)freopen(getValue(),"w",stdout);
}



void GFanApplication::LogLevelOption::onOptionsParsed()
{
  if(getValue()==9)
    {
      Printer::setAssertOnPrinting(true);
    }
  else
    if(getValue())
      {
	setLogLevel(getValue());
      }
}


void GFanApplication::XmlOption::onOptionsParsed()
{
  if(getValue())PolymakeFile::forceXml();
}


void GFanApplication::assertSymmetriesMatch(IntegerVectorList const &g, class PolynomialSet &markedGroebnerBasis, FieldMatrix const *torusActions, bool asPolynomials)
{
  PolynomialSet markedGroebnerBasis2=markedGroebnerBasis;
  markedGroebnerBasis2.sort_();
  int n=markedGroebnerBasis2.numberOfVariablesInRing();
  int I=0;
  for(IntegerVectorList::const_iterator i=g.begin();i!=g.end();i++)
    {
      if(i->size()!=n)
	{
	  fprintf(Stderr,"PERMUTATION ");
	  AsciiPrinter(Stderr).printVector(*i);
	  fprintf(Stderr," HAS WRONG LENGTH.\n");
	  assert(0);
	}
      if(!SymmetryGroup::isPermutation(*i))
	{
	  fprintf(Stderr,"VECTOR ");
	  AsciiPrinter(Stderr).printVector(*i);
	  fprintf(Stderr," DOES NOT REPRESENT A PERMUTATION.\n");
	  assert(0);
	}
      //log0 AsciiPrinter(Stderr).printPolynomialSet(markedGroebnerBasis);
      PolynomialSet b2=SymmetryGroup::permutePolynomialSet(markedGroebnerBasis2,*i);
      //log0 AsciiPrinter(Stderr).printPolynomialSet(b2);
      if(torusActions)b2=b2.torusAct((*torusActions)[I]);
      //log0 AsciiPrinter(Stderr).printPolynomialSet(b2);

      if(!asPolynomials)
	{
	  if(!areIdealsEqual(markedGroebnerBasis2,b2))
	    {
	      fprintf(Stderr,"PERMUTATION ");
	      AsciiPrinter(Stderr).printVector(*i);
	      fprintf(Stderr," DOES NOT KEEP THE IDEAL FIXED.\n");
	      log1 fprintf(Stderr,"Original ideal:\n");
	      log1 AsciiPrinter(Stderr).printPolynomialSet(markedGroebnerBasis2);
	      log1 fprintf(Stderr,"Permuted  ideal:\n");
	      log1 AsciiPrinter(Stderr).printPolynomialSet(b2);

	      log1
		{
		  for(PolynomialSet::const_iterator i=b2.begin();i!=b2.end();i++)
		    {
		      Polynomial remainder=division(*i,markedGroebnerBasis2,LexicographicTermOrder());
		      log1 AsciiPrinter(Stderr).printString("Remainder: ");
		      log1 AsciiPrinter(Stderr).printPolynomial(remainder);
		      log1 AsciiPrinter(Stderr).printNewLine();
		      if(!remainder.isZero())
			{
			  log1 AsciiPrinter(Stderr).printString("Polynomial not in ideal: ");
			  log1 AsciiPrinter(Stderr).printPolynomial(*i);
			  log1 AsciiPrinter(Stderr).printNewLine();
			  break;
			}
		    }

		}

	      assert(0);
	    }
	}
      else
	{
	  b2.markAndScale(LexicographicTermOrder());
	  markedGroebnerBasis2.markAndScale(LexicographicTermOrder());
	  b2.sort_();
	  if(!(markedGroebnerBasis2==b2))
	    {
	      fprintf(Stderr,"PERMUTATION ");
	      AsciiPrinter(Stderr).printVector(*i);
	      fprintf(Stderr," DOES NOT KEEP POLYNOMIAL SET FIXED.\n");

	      log1
		{
		  AsciiPrinter(Stderr).printPolynomialSet(markedGroebnerBasis2);
		  AsciiPrinter(Stderr).printPolynomialSet(b2);

		  PolynomialSet::const_iterator j=markedGroebnerBasis2.begin();
		  for(PolynomialSet::const_iterator i=b2.begin();i!=b2.end();i++)
		    {
		      if(!(*j-*i).isZero())
			{
			  fprintf(Stderr,"Polynomials are different:\n");
			  AsciiPrinter(Stderr).printPolynomial(*i);
			  fprintf(Stderr,"\n");
			  AsciiPrinter(Stderr).printPolynomial(*j);

			  fprintf(Stderr,"\nDifference:\n");
			  AsciiPrinter(Stderr).printPolynomial(*i-*j);

			  fprintf(Stderr,"\n");
			  break;
			}
		      j++;
		    }
		}

	      assert(0);
	    }
	}
      I++;
    }
}

void GFanApplication::onExit()
{
  log1
    {
      fprintf(Stderr,"Number of living FieldImplementation objects:%i\n",FieldImplementation::getNumberOfLivingFieldImplementations());
      fprintf(Stderr,"Number of living FieldElementImplementation objects:%i\n",FieldElementImplementation::getNumberOfLivingFieldElementImplementations());
    }
}

#include "gmp.h"
#include "versioninfo.h"
#include "lp.h"
#include "groebnerengine.h"

class VersionApplication : public Application
{
  bool includeInDefaultInstallation()
  {
    return true;
  }
  const char *helpText()
  {
    return "This program writes out version information of the Gfan installation.\n";
  }
  const char *name()
  {
    return "_version";
  }
  int main()
  {
    cout<<"Gfan version:\n"<<GFAN_RELEASEDIR<<endl<<endl;
    cout<<"Forked from source tree on:\n"<<GFAN_FORKTIME<<endl<<endl;
    cout<<"Linked libraries:"<<endl;
    cout<<"GMP "<<gmp_version<<endl;
    cout<<"Cddlib       YES"<<endl;
    cout<<"SoPlex       "<<((LpSolver::find("SoPlexCddGmp"))?"YES":" NO")<<endl;
    cout<<"Singular     "<<((GroebnerEngine::find("singular"))?"YES":" NO")<<endl;
    return 0;
  }
};

VersionApplication theVersionApplication;

class ListApplication : public Application
{
  SimpleOption hiddenOption;
  bool includeInDefaultInstallation()
  {
    return true;
  }
public:
  ListApplication():
    hiddenOption("--hidden","Show hidden commands which are not officially supported.")
  {
    registerOptions();
  }
  const char *helpText()
  {
    return "This program lists all subcommands of the Gfan installation.\n";
  }
  const char *name()
  {
    return "_list";
  }
  int main()
  {
    cout<<"Available subcommands:"<<endl;
    std::list<Application*> alist2=getSortedApplicationList();
    for(std::list<Application*>::const_iterator p=alist2.begin();p!=alist2.end();p++)
      {
        if(((*p)->includeInDefaultInstallation())^hiddenOption.getValue())
          cout << "gfan"<<(*p)->name()<<"\n";
      }

    cout<<"\nTo get info on a command run the command with option --help.\nFor example run \"gfan _basis --help\".\n";

    return 0;
  }
};

ListApplication theListApplication;
