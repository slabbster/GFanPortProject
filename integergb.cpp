/*
 * integergb.cpp
 *
 *  Created on: Dec 14, 2010
 *      Author: anders
 */

#include "integergb.h"
#include "polynomial.h"
#include "field_rationals.h"
#include "printer.h"
#include "polyhedralcone.h"
#include "wallideal.h"
#include "tropical2.h"

/*
 * Implemented according to [Becker, Weispfenning] Chapter 10.1.
 */

Polynomial dDivision(Polynomial p, PolynomialSet const &l, TermOrder const &termOrder)
{
  PolynomialRing theRing=p.getRing();
  Polynomial r(p.getRing());

  while(!p.isZero())
    {
      p.mark(termOrder);
      Term initial=p.getMarked();
      PolynomialSet::const_iterator i;
      PolynomialSet::iterator j;
      for(i=l.begin();i!=l.end();i++)
        {
          if(i->getMarked().m.exponent.divides(initial.m.exponent))
            if((initial.c*(i->getMarked().c.inverse())).isInteger())break;
        }

      {
        if(i!=l.end())
          {
            Term s(-initial.c*i->getMarked().c.inverse(),Monomial(p.getRing(),initial.m.exponent-i->getMarked().m.exponent));
            p.madd(s,*i);
          }
        else
          {
            p-=initial;
            r+=initial;
          }
      }
    }
  return r;
}


Polynomial eDivision(Polynomial p, PolynomialSet const &l, TermOrder const &termOrder)
{
  PolynomialRing theRing=p.getRing();
  Polynomial r(p.getRing());

//  debug<<"eDivision input "<< p<<l<<termOrder;

  while(!p.isZero())
    {
      p.mark(termOrder);
//      debug<<"P "<<p<<"\n";
      Term initial=p.getMarked();
      PolynomialSet::const_iterator i;
      PolynomialSet::iterator j;
      for(i=l.begin();i!=l.end();i++)
        {
          if(i->getMarked().m.exponent.divides(initial.m.exponent))
            {
//              debug<<"CHECKING:"<<initial<<"\n";
              FieldElement q=integerDivision(initial.c,i->getMarked().c);
              if(!q.isZero())break;
            }
        }

      {
        if(i!=l.end())
          {
//            debug<<"dividing by  "<<*i<<"\n";
            Term s(-integerDivision(initial.c,i->getMarked().c),Monomial(p.getRing(),initial.m.exponent-i->getMarked().m.exponent));
            p.madd(s,*i);
          }
        else
          {
            p-=initial;
            r+=initial;
          }
      }
    }
  return r;
}

static Polynomial sgpol(Polynomial const &g1, Polynomial const &g2, bool s)
{
  PolynomialRing R=g1.getRing();
  FieldElement a1=g1.getMarked().c;
  FieldElement a2=g2.getMarked().c;
  IntegerVector t1=g1.getMarked().m.exponent;
  IntegerVector t2=g2.getMarked().m.exponent;
  FieldElement c1(Q);
  FieldElement c2(Q);
  FieldElement g=gcd(a1,a2,c1,c2);
  FieldElement b1=a2*g.inverse();
  FieldElement b2=a1*g.inverse();
  IntegerVector s1=max(t1,t2)-t1;
  IntegerVector s2=max(t1,t2)-t2;

/*  debug
  <<"a1 "<<a1<<"\n"
  <<"a2 "<<a2<<"\n"
//  <<"t1 "<<t1<<"\n"
//  <<"t2 "<<t2<<"\n"
  <<"c1 "<<c1<<"\n"
  <<"c2 "<<c2<<"\n"
  <<"g  "<<g<<"\n"
  <<"b1 "<<b1<<"\n"
  <<"b2 "<<b2<<"\n"
  <<"s1 "<<s1<<"\n"
  <<"s2 "<<s2<<"\n";
*/
  if(s)return Term(b1,Monomial(R,s1))*g1-Term(b2,Monomial(R,s2))*g2;
  return Term(c1,Monomial(R,s1))*g1+Term(c2,Monomial(R,s2))*g2;
}


Polynomial spol(Polynomial const &g1, Polynomial const &g2)
{
  return sgpol(g1,g2,true);
}


Polynomial gpol(Polynomial const &g1, Polynomial const &g2)
{
  return sgpol(g1,g2,false);
}

void zMinimize(PolynomialSet &F)
{//can be done in nlogn time, but that requirese sorting according to coefficients.
  for(PolynomialSet::iterator i=F.begin();i!=F.end();)
    {
      bool doDelete=false;
      for(PolynomialSet::const_iterator j=F.begin();j!=F.end();j++)
        if(i!=j)
          if(j->getMarked().m.exponent.divides(i->getMarked().m.exponent))
            if((j->getMarked().c.inverse()*i->getMarked().c).isInteger()){doDelete=true;break;}
      if(doDelete)
        {
          PolynomialSet::iterator i2=i;
          i++;
          F.erase(i2);
        }
      else
        i++;
    }
}


void zAutoReduce(PolynomialSet *g, TermOrder const &termOrder)
{
  for(PolynomialSet::iterator i=g->begin();i!=g->end();i++)
    {
//      debug<<"1\n";
      Polynomial temp(*i);
      PolynomialSet::iterator tempIterator=i;
      tempIterator++;
      g->erase(i);
      Monomial monomial=temp.getMarked().m;
      if(temp.getMarked().c.sign()>0)
        g->insert(tempIterator,eDivision(temp,*g,termOrder));
      else
        g->insert(tempIterator,eDivision(temp.getRing().zero()-temp,*g,termOrder));
      tempIterator--;
      i=tempIterator;
      if(i->isZero())
        {
          assert(0);
        }
      else
        i->mark(monomial);
    }
}



typedef pair<Polynomial,Polynomial> Pair;
void zBuchberger(PolynomialSet &F, TermOrder const &T)
{
  list<Pair> B;
  PolynomialSet G(F.getRing());
  for(PolynomialSet::const_iterator i=F.begin();i!=F.end();i++)if(!i->isZero())G.push_back(*i);
  G.mark(T);
  for(PolynomialSet::const_iterator i=G.begin();i!=G.end();i++)
    for(PolynomialSet::const_iterator j=G.begin();j!=i;j++)
      B.push_back(Pair(*i,*j));
  list<Pair> D;
  list<Pair> C=B;

  while(!B.empty())
    {
//      debug<<"Looping\n";
//      debug<<G;

      while(!C.empty())
        {
          Pair f=C.front();
          C.pop_front();
          PolynomialSet::const_iterator gi=G.begin();
          f.first.mark(T);
          f.second.mark(T);
          IntegerVector lcm=max(f.first.getMarked().m.exponent,f.second.getMarked().m.exponent);
          FieldElement HCf1=f.first.getMarked().c;
          FieldElement HCf2=f.second.getMarked().c;
          for(;gi!=G.end();gi++)
            {
              if(gi->getMarked().m.exponent.divides(lcm)
                  && (gi->getMarked().c.inverse()*HCf1).isInteger()
                  && (gi->getMarked().c.inverse()*HCf2).isInteger())break;
            }
          if(gi==G.end())
            {
  //            debug<<f.first.getMarked()<<"  "<<f.second.getMarked()<<"\n";

              Polynomial h=gpol(f.first,f.second);
              Polynomial h0=dDivision(h,G,T);
              if(h0.isZero())
                {
  //                debug<<"remainder is zero!\n";
                }
              else
                {
                  h0.mark(T);
                  for(PolynomialSet::const_iterator gi=G.begin();gi!=G.end();gi++)
                    D.push_back(Pair(*gi,h0));
                  G.push_back(h0);
                }
            }
        }
      Pair f=B.front();
      B.pop_front();
      Polynomial h=spol(f.first,f.second);
      Polynomial h0=dDivision(h,G,T);
      h0.mark(T);
      if(!h0.isZero())
        {
          for(PolynomialSet::const_iterator gi=G.begin();gi!=G.end();gi++)
            D.push_back(Pair(*gi,h0));
          G.push_back(h0);
        }
      C=D;
      B.splice(B.end(),D);
    }
  F=G;
  zMinimize(F);
//debug<<F<<"---------------------------------------------------------------------\n";
  zAutoReduce(&F,T);
}




void IntegerGroebnerFanTraverser::updatePolyhedralCone()
{
  IntegerVectorList inequalities=fastNormals(wallInequalities(groebnerBasis));
  IntegerVectorList empty;
  theCone=PolyhedralCone(inequalities,empty,n);
  theCone.canonicalize();
 // debug<<theCone;
 /// assert(0);
}

IntegerGroebnerFanTraverser::IntegerGroebnerFanTraverser(PolynomialSet const &generators):
  n(generators.getRing().getNumberOfVariables()),
  groebnerBasis(generators)
{
//  LexicographicTermOrder tieBreaker;
//  zBuchberger(groebnerBasis,tieBreaker);
//debug<<"WE ASSUME THAT WE ALREADY HAVE A REDUCED GB\n";

  updatePolyhedralCone();
}

void IntegerGroebnerFanTraverser::changeCone(IntegerVector const &ridgeVector, IntegerVector const &rayVector)
{
  IntegerVectorList tOld;
  tOld.push_back(ridgeVector);
  tOld.push_back(-rayVector);
  MatrixTermOrder TOld(tOld);

  IntegerVectorList t;
  t.push_back(ridgeVector);
  t.push_back(rayVector);
  MatrixTermOrder T(t);

//debug<<"FLIP!\n";
//debug<<"OLDGB\n"<<groebnerBasis;
PolynomialSet g=initialFormsAssumeMarked(groebnerBasis,ridgeVector);
//debug<<"OLD INITIAL GB\n"<<g;
zBuchberger(g,T);
//debug<<"NEW INITIAL GB\n"<<g;



PolynomialSet liftedBasis(groebnerBasis.getRing());
for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
  {
    Polynomial lt=*i;//Polynomial(i->getMarked());
    Polynomial r=eDivision(lt, groebnerBasis,TOld);
    Polynomial rIn=eDivision(lt,initialFormsAssumeMarked(groebnerBasis,ridgeVector),TOld);
    liftedBasis.push_back(lt-r);

//   debug<<"--------------------\n";
//    debug<<"lt:"<<lt<<"\n";
//    debug<<"r:"<<r<<"\n";
//    debug<<"rIn:"<<rIn<<"\n";
//    debug<<"--------------------\n";
  }
liftedBasis.mark(T);
//debug<<"to be autoreduced:"<<liftedBasis;
//debug<<"autoreducing......\n";
liftedBasis.mark(T);
zAutoReduce(&liftedBasis,T);
PolynomialSet newInitialForms=initialForms(liftedBasis,ridgeVector);
//debug<<"LIFTED"<<liftedBasis;
//debug<<"NEWINITIALFORMS"<<newInitialForms;

/*zBuchberger(groebnerBasis,T);

{
  groebnerBasis.sort_();
  liftedBasis.sort_();
  bool areEqual=(groebnerBasis.size()==liftedBasis.size());
  if(areEqual)
    {
      PolynomialSet::const_iterator i=groebnerBasis.begin();
      for(PolynomialSet::const_iterator j=liftedBasis.begin();j!=liftedBasis.end();j++,i++)
        areEqual&=(*i-*j).isZero();
    }
  if(!areEqual)
    {
      debug<<groebnerBasis<<liftedBasis;
      assert(0);
    }
}*/

groebnerBasis=liftedBasis;
groebnerBasis.mark(T);


//debug<<"NEW BASIS"<<groebnerBasis;

  updatePolyhedralCone();
}

IntegerVectorList IntegerGroebnerFanTraverser::link(IntegerVector const &ridgeVector)
{
  IntegerVectorList ret;
  IntegerVector v=theCone.link(ridgeVector).getUniquePoint();
  ret.push_back(v);

  PolyhedralCone temp=intersection(PolyhedralCone::positiveOrthant(n),theCone.faceContaining(ridgeVector));
  IntegerVector temp2=temp.getRelativeInteriorPoint();
  if(temp2.min()>0)
    {
      ret.push_back(-v);
    }
//  debug<<"LINK"<<ret;
  return ret;
}


PolyhedralCone & IntegerGroebnerFanTraverser::refToPolyhedralCone()
{
  return theCone;
}
