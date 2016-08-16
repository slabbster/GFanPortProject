#include "vektor.h"
#include "printer.h"
#include "parser.h"
#include "gfanapplication.h"
#include "minkowskisum.h"
#include "newtonpolytope.h"
#include "buchberger.h"
#include "wallideal.h"
#include "lp.h"
#include "tropical.h"
#include "division.h"
#include "bergman.h"
#include "tropical2.h"
#include "dimension.h"
#include "timer.h"
#include "log.h"
#include "linalg.h"
#include "tropicaltraverse.h"
#include "traverser_tropical.h"
#include "symmetrictraversal.h"
#include "traverser_stableintersection.h"

class SymmetriesApplication : public GFanApplication
{
  SimpleOption optionSymmetry;
  SimpleOption optionTorusSymmetry;
public:
  bool includeInDefaultInstallation()
  {
    return true;
  }
  const char *helpText()
  {
    return "This program computes the symmetries of a polynomial ideal. The program is slow, so think before using it. Use --symmetry to give hints about which subgroup of the symmetry group could be useful. The program checks each element of the specified subgroup to see if it preserves the ideal.\n";
  }
  SymmetriesApplication():
    optionSymmetry("--symmetry","Specify subgroup to be searched for permutations keeping the ideal fixed."),
    optionTorusSymmetry("--symsigns","Specify for each generator of the group specified wiht --symmetry an element of ${-1,+1}^n$ which by its multiplication on the variables together with the permutation is expected to keep the ideal fixed.")
  {
    registerOptions();
  }
  const char *name()
  {
    return "_symmetries";
  }
  int main()
  {
    FileParser P(Stdin);

    AsciiPrinter p(Stdout);
    PolynomialSet gb=P.parsePolynomialSetWithRing();
    int n=gb.numberOfVariablesInRing();

    WeightReverseLexicographicTermOrder T(IntegerVector::allOnes(n));
    buchberger(&gb,T);

    assert(n>=2);

    SymmetryGroup signSymmetries(2*n);
    if(optionSymmetry.getValue())
      {
    	IntegerVectorList generators=P.parseIntegerVectorList();

    	IntegerMatrix torusAction(generators.size(),n);
        if(optionTorusSymmetry.getValue())
        {
        	torusAction=rowsToIntegerMatrix(P.parseIntegerVectorList());
        }
        else
        {
        	for(int i=0;i<torusAction.getHeight();i++)
            	for(int j=0;j<torusAction.getWidth();j++)
            		torusAction[i][j]=1;
        }
        IntegerVectorList doubleGenerators;
        int I=0;
        for(IntegerVectorList::const_iterator i=generators.begin();i!=generators.end();i++,I++)
        	doubleGenerators.push_back(SymmetryGroup::combinePermutationAndSignChanges(*i,torusAction[I]));
    	signSymmetries.computeClosure(doubleGenerators);
      }
    else
    {
    	if(optionTorusSymmetry.getValue())
    	{
    		debug<<"Option --symsigns can only be used together with --symmetry\n";
    		assert(0);
    	}
    	IntegerVectorList doubleGenerators;
    	IntegerVector cycle(n);
    	for(int i=0;i<n;i++)cycle[i]=((i+1)%n);
    	IntegerVector transposition=SymmetryGroup::identity(n);
    	transposition[0]=1;
    	transposition[1]=0;
    	doubleGenerators.push_back(SymmetryGroup::combinePermutationAndSignChanges(cycle,IntegerVector::allOnes(n)));
    	doubleGenerators.push_back(SymmetryGroup::combinePermutationAndSignChanges(transposition,IntegerVector::allOnes(n)));
    	signSymmetries.computeClosure(doubleGenerators);
    }

    IntegerVectorList doubleGenerators2;

    for(SymmetryGroup::ElementContainer::const_iterator i=signSymmetries.elements.begin();i!=signSymmetries.elements.end();i++)
    {
    	IntegerVector permutation,signChanges;
    	SymmetryGroup::extractPermuationAndSignChanges(*i,permutation,signChanges);

    	log1 debug<<"Checking "<<permutation<<" "<<signChanges<<"\n";
        PolynomialSet b2=SymmetryGroup::permutePolynomialSet(gb,permutation);
        //log0 AsciiPrinter(Stderr).printPolynomialSet(b2);

        b2=b2.torusAct(integerVectorToFieldVector(signChanges,Q));
        if(areIdealsEqual(gb,b2))
        {
        	doubleGenerators2.push_back(*i);
        	log1 debug<<"OK";
        }
    }

    IntegerVectorList permutations;
    IntegerVectorList signChangess;
    for(IntegerVectorList::const_iterator i=doubleGenerators2.begin();i!=doubleGenerators2.end();i++)
    {
    	IntegerVector permutation,signChanges;
    	SymmetryGroup::extractPermuationAndSignChanges(*i,permutation,signChanges);
    	permutations.push_back(permutation);
    	signChangess.push_back(signChanges);
    }

    pout<<permutations<<signChangess;

    return 0;
  }
};

static SymmetriesApplication theApplication;
