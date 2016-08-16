#include "singularconversion.h"

#include "polynomialring.h"
#include "polynomial.h"
#include "field_rationals.h"
#include "buchberger.h"
#include "division.h"
#include "printer.h"
#include "log.h"
#include <iostream>


void singularBuchberger(PolynomialSet *g, TermOrder const &termOrder)
{
  ring R=singularRing(g->getRing());

  ideal i=singularPolynomialSet(*g);

  //  ideal j=kStd(i,NULL,testHomog,NULL);
//  test|=(Sy_bit(OPT_REDSB)|Sy_bit(OPT_REDTAIL)|Sy_bit(OPT_INTSTRATEGY));
  test|=(Sy_bit(OPT_REDSB)|Sy_bit(OPT_REDTAIL));

  ideal j=kStd(i,NULL,testHomog,NULL);
  //idShow(j);

  idDelete(&i);
  PolynomialSet ret=fromSingularIdeal(g->getRing(),j);
  idDelete(&j);


  freeSingularRing(R);

  *g=ret;
  //  return ret;
}

/**************************************************************************************************************************/

#include "groebnerengine.h"

static ring mySingularRingDegRevLex(PolynomialRing const &r, IntegerVector const &weight)
{
  //  FieldRationalsImplementation *q=dynamic_cast<FieldRationalsImplementation>r.getField().implementingObject;

  ring ret=(ring)omAlloc0(sizeof(sip_sring));

  if(r.getField().isRationals())
    {
      ret->ch=0;
    }
  else
    {
      assert(0);
    }


  ret->N=r.getNumberOfVariables();

  ret->names=(char**) omAlloc(ret->N*sizeof(char*));
  for(int i=0;i<ret->N;i++)
    ret->names[i]=omStrDup(r.getVariableName(i).c_str());

  ret->order=(int*)omAlloc0(3*sizeof(int));
  ret->block0=(int*)omAlloc0(3*sizeof(int));
  ret->block1=(int*)omAlloc0(3*sizeof(int));

  ret->order[0]=ringorder_wp;//degree revlex
  // ret->order[0]=ringorder_Wp;//degree lex
  ret->block0[0]=1;ret->block1[0]=ret->N;
  ret->order[1]=ringorder_C;

  ret->wvhdl=(int**)omAlloc0(3*sizeof(int*));
  ret->wvhdl[0]=(int*)omAlloc(ret->N*sizeof(int));
  {
	  IntegerVector weight2=weight+(1-weight.min())*IntegerVector::allOnes(weight.size());
	  //debug<<"WEIGHT FOR SINGULAR"<<weight2<<"\n";
  for(int i=0;i<ret->N;i++)
    ret->wvhdl[0][i]=weight2[i];//FIXMEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
  }
  rComplete(ret);
 // ret->options|=(Sy_bit(OPT_REDSB)|Sy_bit(OPT_REDTAIL)|Sy_bit(OPT_INTSTRATEGY)|Sy_bit(OPT_REDTHROUGH));
//  cerr<<"ringoptions"<<int(ret->options)<<"\n";

 // cerr<<"test"<<int(test)<<"\n";
  rChangeCurrRing(ret);
 // cerr<<"test"<<int(test)<<"\n";

  return ret;
}


class GroebnerEngineSingular : public GroebnerEngine
{
  virtual PolynomialSet groebnerBasis(bool &success, PolynomialSet const &idealGenerators, TermOrder const &termOrder, bool autoreduce)
  {
    PolynomialSet ret(idealGenerators.getRing());
    if(!ret.getRing().getField().isRationals())
      {
	success=false;
	return ret;
      }
    if(!idealGenerators.isHomogeneous())
    {
    	success=false;
    	return ret;
    }
    if(dynamic_cast<const WeightReverseLexicographicTermOrder*> (&termOrder))
      {
const WeightReverseLexicographicTermOrder *T=dynamic_cast<const WeightReverseLexicographicTermOrder*> (&termOrder);

ring R=mySingularRingDegRevLex(idealGenerators.getRing(),T->getWeight());
	ideal i=singularPolynomialSet(idealGenerators);
        test|=(Sy_bit(OPT_REDSB)|Sy_bit(OPT_REDTAIL)|Sy_bit(OPT_INTSTRATEGY));
        test|=(Sy_bit(OPT_REDTHROUGH));

	log2 cerr<<"calling singular\n";
	//  debug<<"test"<<int(test)<<"\n";
	ideal j=kStd(i,NULL,testHomog,NULL);
	log2 cerr<<"returning from singular\n";

	idDelete(&i);
	ret=fromSingularIdeal(ret.getRing(),j);
	idDelete(&j);
	freeSingularRing(R);

	ret.markAndScale(termOrder);
        minimize(&ret);//<---- This is needed for test case 1000 in the suite. Of course it would be nicer if Singular could compute the reduced GB itself.
        autoReduce_(&ret,termOrder);//------------------------------------------REMOVE AND TRUST SINGULAR
	if(!isReduced(ret))
		{
			/*
			 * This test fails on the example
			 * gfan _tropicaltraverse --symmetry --log2 <hannah.cone
			 * The monomial d^3*e*f^3*g^8*h^4*i^2 appears as a leading term, but also in a tail.
			 */

		debug<<idealGenerators;
		debug<<ret;
		  autoReduce_(&ret,termOrder);
			debug<<ret;
			assert(0);
		}
	assert(isMinimal(ret));
	assert(isReduced(ret));
	minimize(&ret);
	if(autoreduce)
	  autoReduce_(&ret,termOrder);

	success=true;
	return ret;
      }
    success=false;
    return ret;
  }
  virtual PolynomialSet autoReduce(bool &success, PolynomialSet const &idealGenerators)
  {
    PolynomialSet ret(idealGenerators.getRing());
    if(!ret.getRing().getField().isRationals())
      {
	success=false;
	return ret;
      }
    if(!idealGenerators.isHomogeneous())
    {
    	success=false;
    	return ret;
    }

    IntegerVector weight=termorderWeight(idealGenerators);

    ring R=mySingularRingDegRevLex(idealGenerators.getRing(),weight);
    ideal i=singularPolynomialSet(idealGenerators);
    test|=(Sy_bit(OPT_REDSB)|Sy_bit(OPT_REDTAIL)|Sy_bit(OPT_INTSTRATEGY));

    log2 cerr<<"calling singular\n";
    ideal j=kStd(i,NULL,testHomog,NULL);
	//    ideal j=kInterRed(i);
    log2 cerr<<"returning from singular\n";

    idDelete(&i);
    ret=fromSingularIdeal(ret.getRing(),j);
    idDelete(&j);
    freeSingularRing(R);

    ret.markAndScale(WeightTermOrder(weight));
    assert(isReduced(ret));

    success=true;
    return ret;
  }
  virtual const char* name()
  {
    return "singular";
  }
};

static GroebnerEngineSingular groebnerEngineSingular;
