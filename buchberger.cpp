#include "buchberger.h"

#include <set>
#include <algorithm>
#include <iostream>

#include "division.h"
#include "printer.h"
#include "timer.h"
#include "parser.h"
#include "log.h"

static Timer buchbergerTimer("Buchberger",10);

Polynomial sPolynomial(Polynomial a, Polynomial b)
{
  bool comments=false;
  if(comments)
    {
      AsciiPrinter(Stderr).printString("S(");
      AsciiPrinter(Stderr).printPolynomial(a);
      AsciiPrinter(Stderr).printString(",");
      AsciiPrinter(Stderr).printPolynomial(b);
      AsciiPrinter(Stderr).printString(")=");
    }

  //marked coefficient of a and b must be one

  IntegerVector ina=a.getMarked().m.exponent;
  IntegerVector inb=b.getMarked().m.exponent;

  IntegerVector L=max(ina,inb);

  if(comments)
    AsciiPrinter(Stderr).printVector(L);


  FieldElement const &f=a.getMarked().c;


  Polynomial A=a;A*=Term(f.one(),Monomial(a.getRing(),L-ina));
  Polynomial B=b;B*=Term(f.one(),Monomial(b.getRing(),L-inb));

  if(comments)
    {
      AsciiPrinter(Stderr).printPolynomial(A-B);
      AsciiPrinter(Stderr).printString("\n");
    }

  return A-B;
}


// Simple Buchberger

void buchberger/*Simple*/(PolynomialSet *g, TermOrder const &termOrder, bool allowSaturation)
{
  PolynomialRing theRing=g->getRing();
  //  log2 fprintf(Stderr,"ENTERING buchberger\n");
  TimerScope ts(&buchbergerTimer);
  PolynomialSet sPolynomials(theRing);

  for(PolynomialSet::const_iterator i=g->begin();i!=g->end();i++)
    if(!i->isZero())sPolynomials.push_back(*i); // It is safe and useful to ignore the 0 polynomial

  if(allowSaturation)sPolynomials.saturate();
  sPolynomials.markAndScale(termOrder);

  *g=PolynomialSet(theRing);

  while(!sPolynomials.empty())
    {
      Polynomial p=*sPolynomials.begin();
      sPolynomials.pop_front();

      p=division(p,*g,termOrder);
      if(!p.isZero())
        {
          if(allowSaturation)p.saturate();
          p.mark(termOrder);
          p.scaleMarkedCoefficientToOne();
	  bool isMonomial=p.isMonomial();
          for(PolynomialSet::const_iterator i=g->begin();i!=g->end();i++)
	    if((!isMonomial) || (!i->isMonomial())) // 2 % speed up!
            {
              if(!relativelyPrime(i->getMarked().m.exponent,p.getMarked().m.exponent))
                {
                  Polynomial s=sPolynomial(*i,p);
                  s.mark(termOrder); // with respect to some termorder
                  s.scaleMarkedCoefficientToOne();
                  sPolynomials.push_back(s);
                }
            }
          g->push_back(p);
	  {
	    static int t;
	    t++;
	    //	    if((t&31)==0)fprintf(Stderr," gsize %i  spolys:%i\n",g->size(),sPolynomials.size());
	  }
        }
    }
  //log2  fprintf(Stderr," buchberger minimize\n");
  minimize(g);
  //log2 fprintf(Stderr," buchberger autoreduce\n");
  autoReduce(g,termOrder);
  //log2 fprintf(Stderr,"LEAVING buchberger\n\n");
}



// Buchberger with chain criterion

struct ChainPair
{
  PolynomialSet::const_iterator a,b;
  int A,B;
  IntegerVector lcm;
  ChainPair(PolynomialSet::const_iterator const &a_,PolynomialSet::const_iterator const &b_,int A_,int B_):
    a(a_),
    b(b_),
    A(A_),
    B(B_),
    lcm(max(a_->getMarked().m.exponent,b_->getMarked().m.exponent))
  {
  }
  bool operator<(const ChainPair & b)const
  {
    if(b.lcm.sum()<lcm.sum())return false;
    if(lcm.sum()<b.lcm.sum())return true;
    if(b.lcm<lcm)return true;
    if(lcm<b.lcm)return false;
    if(A<b.A)return true;
    if(A>b.A)return false;
    if(B<b.B)return true;
    if(B>b.B)return false;
    assert(0);
  }
};

typedef set<ChainPair> ChainPairList;

static bool canBeRemoved(ChainPairList const &P, IntegerVector const &lcm, PolynomialSet::const_iterator i, PolynomialSet::const_iterator l)
{
  //  return false;
  for(ChainPairList::const_iterator t=P.begin();t!=P.end();t++)
    {
      if(t->a==i && t->b!=l && t->b->getMarked().m.exponent.divides(lcm)/* ||
	 t->b==i && t->a!=l && t->a->getMarked().m.exponent.divides(lcm) ||
	 t->a==l && t->b!=i && t->b->getMarked().m.exponent.divides(lcm) ||
	 t->b==l && t->a!=i && t->a->getMarked().m.exponent.divides(lcm)*/)return true;
    }
  return false;
}

void printPairs(ChainPairList const &P)
{
  return;
  for(ChainPairList::const_iterator t=P.begin();t!=P.end();t++)
    {
      cerr<<"("<<t->A<<","<<t->B<<")[";
      AsciiPrinter(Stderr)<<t->a->getMarked()<<","<<t->b->getMarked()<<"]<"<<t->lcm <<">";
    }
  cerr<<endl;
}

void buchbergerChain(PolynomialSet *g, TermOrder const &termOrder)
{
  PolynomialRing theRing=g->getRing();
  TimerScope ts(&buchbergerTimer);

  {
    PolynomialSet g2(theRing);
    for(PolynomialSet::const_iterator i=g->begin();i!=g->end();i++)
      if(!i->isZero())g2.push_back(*i); // It is safe and useful to ignore the 0 polynomial

    *g=g2;
  }
  g->markAndScale(termOrder);

  ChainPairList P;//use better data structure for this
  ChainPairList Pchecked;//use better data structure for this
  int I=0;
  for(PolynomialSet::const_iterator i=g->begin();i!=g->end();i++,I++)
    {
      int J=0;
      for(PolynomialSet::const_iterator j=g->begin();j!=i;j++,J++)
	{
	  P.insert(ChainPair(j,i,J,I));//here
	}
    }

  while(!P.empty())
    {
      //   cerr<<"P";printPairs(P);cerr<<"Pchecked";printPairs(Pchecked);


      PolynomialSet::const_iterator i=P.begin()->a;
      PolynomialSet::const_iterator l=P.begin()->b;
      int I=P.begin()->A;
      int L=P.begin()->B;

      if(relativelyPrime(i->getMarked().m.exponent,l->getMarked().m.exponent) || (i->isMonomial() && l->isMonomial()))
	{
	  // Pchecked.push_back(P.front());
	  // P.pop_front();
	  Pchecked.insert(*P.begin());
	  P.erase(P.begin());
	}
      else
	{
	  IntegerVector lcm=max(i->getMarked().m.exponent,l->getMarked().m.exponent);

	  if(canBeRemoved(P,lcm,i,l) || canBeRemoved(Pchecked,lcm,i,l))
	    {
	      //cerr<<"removin"<<endl;
	      //P.pop_front(); // This might remove elements from P with a trivial S-poly
	      P.erase(P.begin()); // This might remove elements from P with a trivial S-poly
	    }
	  else
	    {
	      Polynomial p=sPolynomial(*i,*l);

	      p.mark(termOrder);
	      p=division(p,*g,termOrder);

	      // Pchecked.push_back(P.front());
	      // P.pop_front();
	      Pchecked.insert(*P.begin());
	      P.erase(P.begin());

	      if(!p.isZero())
		{
		  p.mark(termOrder);
		  p.scaleMarkedCoefficientToOne();
		  int K=g->size();
		  g->push_back(p);
		  PolynomialSet::const_iterator k=g->end();k--;
		  int I=0;
		  for(PolynomialSet::const_iterator i=g->begin();i!=k;i++,I++)
		    {
		      P.insert(ChainPair(i,k,I,K));//here
		    }
		}
	    }
	}
    }
  //  AsciiPrinter(Stderr)<<*g;
  minimize(g);
  autoReduce(g,termOrder);
}


// Sugar Cube Buchberger

static const TermOrder * SPairTermOrder;
class SPair
{
public:
  PolynomialSet::const_iterator i,j;
  //int indexi;
  int indexj;
  IntegerVector v;
  int sugar;
  SPair(PolynomialSet::const_iterator i_, PolynomialSet::const_iterator j_, int indexj_):
    i(i_),
    j(j_),
    //indexi(indexi_),
    indexj(indexj_)
  {
    //    assert(indexi<indexj);
    int sugar1=j->getSugar()-j->getMarked().m.exponent.sum();
    int sugar2=i->getSugar()-i->getMarked().m.exponent.sum();
    if(sugar1>sugar2)
      sugar=sugar1;
    else
      sugar=sugar2;
    sugar+=max(i->getMarked().m.exponent,j->getMarked().m.exponent).sum();
    //  fprintf(Stderr,"%i\n",sugar);
    //    sugar=(max(i->getMarked().m.exponent,j->getMarked().m.exponent)).sum();
    v=max(i->getMarked().m.exponent,j->getMarked().m.exponent);
  }
  Polynomial sPolynomial_()const{return sPolynomial(*i,*j);};
  bool relativelyPrime_()const{return relativelyPrime(i->getMarked().m.exponent,j->getMarked().m.exponent);}
  bool operator<(const class SPair &a)const;
  void print(Printer &p)const
  {
    p.printString("S");
    //p.printInteger(indexi);
    p.printString(",");
    p.printInteger(indexj);
    p.printString("(");
    p.printPolynomial(*i);
    p.printString(",");
    p.printPolynomial(*j);
    p.printString(")=");
    p.printPolynomial(sPolynomial_());
    p.printString(" ");
    p.printVector(v);
    p.printString(" sugar: ");
    p.printInteger(sugar);
    p.printString("\n");
  };
};


bool SPair::operator<(const class SPair &a)const //partial order
{
  //  if(sugar<a.sugar)return true;
  //  if(sugar>a.sugar)return false;
  //  int d=v.sum()-a.v.sum();
  //  if(d<0)return false;
  //  if(d>0)return true;

  if((*SPairTermOrder)(v,a.v))return true;
  if((*SPairTermOrder)(a.v,v))return false;

      if(indexj<a.indexj)return true;
      if(indexj>a.indexj)return false;
      //if(indexj>a.indexj)return true;
      //if(indexj<a.indexj)return false;
      // fprintf(Stderr,"Tie breaking\n");

      AsciiPrinter P(Stderr);
      print(P);
      a.print(P);
      assert(0);
}


struct SPairLCMCompare
{
  bool operator()(const SPair &a, const SPair &b)
  {
    return (*SPairTermOrder)(a.v,b.v);
    /*    if(LexicographicTermOrder()(a.v,b.v))return true;
    if(LexicographicTermOrder()(b.v,a.v))return false;
    if(a.indexj<b.indexj)return true;
    if(a.indexj>b.indexj)return false;
    return false;
    */
  }
};

typedef set<SPair> SPairSet;
typedef list<SPair> SPairList;

void printSPairSet(SPairSet const &l)
{
  AsciiPrinter p(Stderr);
  for(SPairSet::const_iterator i=l.begin();i!=l.end();i++)i->print(p);
}
void printSPairList(SPairList const &l)
{
  AsciiPrinter p(Stderr);
  for(SPairList::const_iterator i=l.begin();i!=l.end();i++)i->print(p);
}


void updatePairs(SPairSet &sPolynomials, PolynomialSet *g, PolynomialSet *redundantOnes, PolynomialSet::iterator const &j, int indexj, int truncationDegree, IntegerVector const *grading)
{
  if(truncationDegree)
    if(dot(*grading,j->getMarked().m.exponent)>truncationDegree)return;
  //divisionPolynomials.push_front(*j);
  // fprintf(Stderr,"1\n");
  // 1 : Delete pairs in "sPolynomials"
  bool skip=false;
  //if(0)
  for(SPairSet::iterator ij=sPolynomials.begin();ij!=sPolynomials.end();)
    {
      skip=false;
      if(j->getMarked().m.exponent.divides(ij->v))
	if(ij->v!=max(j->getMarked().m.exponent,ij->i->getMarked().m.exponent))
	  if(ij->v!=max(j->getMarked().m.exponent,ij->j->getMarked().m.exponent))
	      {
		SPairSet::iterator temp=ij;
		ij++;skip=true;
		sPolynomials.erase(temp);
	      }
      if(!skip)ij++;
    }

  // fprintf(Stderr,"1b\n");
  //   : Generate new pairs in "newOnes" and move redundant polynomials from "g" into "redundantOnes"
  SPairList newOnes;
  {
    int indexi=0;
    for(PolynomialSet::iterator i=g->begin();i!=j;)
      {
	SPair s=SPair(i,j,indexj);
	if(truncationDegree==0 || dot(*grading,s.v)<=truncationDegree)
	  newOnes.push_back(s);
	if(j->getMarked().m.exponent.divides(i->getMarked().m.exponent))  //redundancy removal
	  {
	    PolynomialSet::iterator temp=i;
	    i++;
	    redundantOnes->splice(redundantOnes->begin(),*g,temp);
	  }
	else
	  {
	    i++;
	  }
	/*	if(j->getMarked().m.exponent.divides(i->getMarked().m.exponent))  //redundancy removal
	  {
	    PolynomialSet::iterator temp=i;
	    i++;
	    redundantOnes->splice(redundantOnes->begin(),*g,temp);
	    newOnes.push_back(SPair(redundantOnes->begin(),j,indexj));
	  }
	else
	  {
	    newOnes.push_back(SPair(i,j,indexj));
	    i++;
	    }*/
	indexi++;
      }
  }

  // fprintf(Stderr,"2\n");
  // 2 : Move one element from "newOnes" to "sPolynomials" for each equivalence class,
  //     but none if tau is "divisible" or representable as "a disjoint union"

  if(!newOnes.empty())
    {
      newOnes.sort(SPairLCMCompare());
      IntegerVectorList divisionList;
      SPairList::const_iterator representative=newOnes.begin();
      bool ignore=representative->relativelyPrime_();
      divisionList.push_back(representative->v);
      for(SPairList::iterator i=newOnes.begin();i!=newOnes.end();i++)
	if(representative->v==i->v)
	  {
	    if(!ignore && i->relativelyPrime_())ignore=true;
	  }
	else
	  {
	    if(!ignore)sPolynomials.insert(*representative);
	    //else divisionList.pop_back();//stupid?
	    representative=i;
	    ignore=representative->relativelyPrime_();
	    for(IntegerVectorList::const_iterator k=divisionList.begin();k!=divisionList.end();k++)
	      if(k->divides(i->v))
		{
		  ignore=true;
		  break;
		}
	    divisionList.push_back(i->v);
	  }
      if(!ignore)sPolynomials.insert(*representative);
    }

  /*
  // fprintf(Stderr,"2\n");
  // 2 : Remove elements from newOnes if tau of element is divisible
  for(SPairList::iterator s=newOnes.begin();s!=newOnes.end();)
    {
      bool remove=false;
      for(SPairList::const_iterator t=newOnes.begin();t!=newOnes.end();t++)
	if(t->v!=s->v)
	  if(t->v.divides(s->v))
	    {
	      remove=true;
	      break;
	    }
      if(remove)
	{
	  SPairList::iterator temp=s;
	  s++;
	  newOnes.erase(temp);
	}
      else
	s++;
    }

  // fprintf(Stderr,"3+4\n");
  // 3+4 : Add "newOnes" to "sPolynomials", but only one for each tau
  newOnes.sort(SPairLCMCompare());
  SPairList::const_iterator best=newOnes.begin();
  for(SPairList::const_iterator i=newOnes.begin();i!=newOnes.end();i++)
    if(best->v==i->v)
      {
	if(i->relativelyPrime_())
	  best=i;
      }
    else
      {
	if(!best->relativelyPrime_())sPolynomials.insert(*best);
	best=i;
      }
  if(best!=newOnes.end() && !best->relativelyPrime_())sPolynomials.insert(*best);
  */
}


void buchberger2(PolynomialSet *g, TermOrder const &termOrder)
{
  PolynomialRing theRing=g->getRing();
  assert(0); // We do not use the Gebauer-Moeller implementation since it is not complete
  bool printComments=false;

  int truncationDegree=0;
  IntegerVector grading;

  grading=StringParser("(1,12223,12224,36674,61119,85569)").parseIntegerVector();
  truncationDegree=89643482;


  //  PolynomialSet divisionPolynomials;
  TimerScope ts(&buchbergerTimer);

  SPairTermOrder=&termOrder;
  g->markAndScale(termOrder);
  g->sort(PolynomialCompareMarkedTermsReverse(termOrder));
  g->computeInitialSugar();
  PolynomialSet redundantOnes(theRing);

  SPairSet sPolynomials;
  int indexj=0;
  for(PolynomialSet::iterator j=g->begin();j!=g->end();j++)
    {
      if(j!=g->begin())
	updatePairs(sPolynomials,g,&redundantOnes,j,indexj,truncationDegree,&grading);
      indexj++;
    }

  //  printSPairSet(sPolynomials);

  int numberOfCriticalPairsConsidered=0;
  int numberOfUsefulCriticalPairs=0;

  while(!sPolynomials.empty())
    {
      Polynomial p=sPolynomial(*sPolynomials.begin()->i,*sPolynomials.begin()->j);
      p.mark(termOrder);
      p.scaleMarkedCoefficientToOne();
      sPolynomials.erase(sPolynomials.begin());

      numberOfCriticalPairsConsidered++;

      //      if(printComments)
	{
	  static int t;
	  t++;
	  if((t&31)==0)fprintf(Stderr,"         spolys:%i\n",sPolynomials.size());
	}

      {
	static int t;
	if(++t==10)
	  {
	    //assert(divisionPolynomials.size()==g->size());
	    //    divisionPolynomials.sort(PolynomialCompareMarkedTerms(termOrder));
	    t=0;
	    //	    fprintf(Stderr,"test\n");
	  }
      }
      //p=division(p,divisionPolynomials,termOrder);
      p=division(p,*g,termOrder);
      if(!p.isZero())
        {
	  p.mark(termOrder);
          p.scaleMarkedCoefficientToOne();
	  g->push_back(p);
	  numberOfUsefulCriticalPairs++;
	  //	  if(printComments)
	    {
	      static int t;
	      if(((++t)&=31)==0)
		fprintf(Stderr,"gsize:%i spolys:%i\n",g->size()+1,sPolynomials.size());
	    }
	  PolynomialSet::iterator j=g->end();j--;
	  updatePairs(sPolynomials,g,&redundantOnes,j,indexj,truncationDegree,&grading);
	  indexj++;
        }
    }
  minimize(g);
  autoReduce(g,termOrder);

  //  if(printComments)
    fprintf(Stderr,"Number of critical pairs considered %i Number of useful critical pairs %i",numberOfCriticalPairsConsidered,numberOfUsefulCriticalPairs);
}


struct polynomialOrder
{
  TermOrder const &order;
  polynomialOrder(TermOrder const &order_):order(order_){}

  bool operator()(const Polynomial &a, const Polynomial &b)
  {
    return order(a.getMarked().m.exponent, b.getMarked().m.exponent);
  }
};

void minimize(PolynomialSet *g)
{//CHECK THAT THIS ROUTINE WORKS IF TWO GENERATORS HAVE THE SAME INITIAL TERM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PolynomialRing theRing=g->getRing();
  PolynomialSet ret(theRing);

  g->sort(polynomialOrder(LexicographicTermOrder()));

  for(PolynomialSet::const_iterator i=g->begin();i!=g->end();i++)
    {
      bool someDivides=false;
      for(PolynomialSet::const_iterator j=ret.begin();j!=ret.end();j++)
	//      for(PolynomialSet::const_iterator j=g->begin();j!=i;j++) //changed Feb 2009
        {
          if(j->getMarked().m.exponent.divides(i->getMarked().m.exponent))
            {
              someDivides=true;
            }
        }
      if(!someDivides)
        ret.push_back(*i);
    }

  *g=ret;
}


void autoReduce(PolynomialSet *g, TermOrder const &termOrder)
{
	/**
	 * TODO: there should be two options : supplying a termorder, or not supplying a termorder. In the latter case this routine should decide if it wants to compute one.
	 */
//  WeightTermOrder termOrder2(termorderWeight(*g));//REMOVE ME ?? JAN 2009


  for(PolynomialSet::iterator i=g->begin();i!=g->end();i++)
    {
      Polynomial temp(*i);
      PolynomialSet::iterator tempIterator=i;
      tempIterator++;
      g->erase(i);
      Monomial monomial=temp.getMarked().m;
      g->insert(tempIterator,division(temp,*g,termOrder));
      //g->insert(tempIterator,smartDivision(temp,*g,termOrder));
      tempIterator--;
      i=tempIterator;
      i->mark(monomial);
    }
}


bool isMarkedGroebnerBasis(PolynomialSet const &g)
{
  int counter=0;
  for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
    {
    log2 fprintf(Stderr,"%i ",counter++);
    for(PolynomialSet::const_iterator j=i;j!=g.end();j++)
      if(!relativelyPrime(i->getMarked().m.exponent,j->getMarked().m.exponent))
	{
	  Polynomial s=sPolynomial(*i,*j);
	  if(!division(s,g,LexicographicTermOrder()).isZero())
	    {
	      log3{AsciiPrinter(Stderr)<<"Spoly("<<*i<<","<<*j<<")="<<sPolynomial(*i,*j)<<" gives remainder "<< division(s,g,LexicographicTermOrder()) <<"\n";}
	      return false;
	    }
	}
    }
  return true;
}


bool isMinimal(PolynomialSet const &g)
{
	PolynomialSet temp=g.markedTermIdeal();
	minimize(&temp);
	return temp.size()==g.size();
}


bool isReduced(PolynomialSet const &g)
{
	if(!isMinimal(g))return false;
	PolynomialSet temp=g.markedTermIdeal();
	for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
		for(TermMap::const_iterator j=i->terms.begin();j!=i->terms.end();j++)
			if(!(j->first.exponent==i->getMarked().m.exponent))
				for(PolynomialSet::const_iterator k=temp.begin();k!=temp.end();k++)
					if(k->getMarked().m.exponent.divides(j->first.exponent))return false;
	return true;
}
