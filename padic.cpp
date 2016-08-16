/*
 * padic.cpp
 *
 *  Created on: Nov 30, 2010
 *      Author: anders
 */

#define DEBUG if(0)

#include "padic.h"
#include "polynomial.h"
#include <algorithm>
#include "field_zmodpz.h"
#include "printer.h"
#include "wallideal.h"

using namespace std;

static int64 smartDot(IntegerVector const &v, IntegerVector const &omega)
{
  int64 ret=0;
  int n=v.size();
  for(int i=0;i<n;i++)ret+=v[i]*omega[i+1];
  return ret;
}

PolynomialRing residuePolynomialRing(PolynomialRing const &R, int prime)
{
  return PolynomialRing(FieldZModPZ(prime),R.getVariableNames());
}

static int64 omegaDegree(TermMap::const_iterator i, int prime, IntegerVector const &omega)
{
//  DEBUG debug<<"-----\n"<<i->second<<"  -----   "<<i->second.pAdicValuation(prime)<<"\n";

//  return omega[0]*(i->second.pAdicValuation(prime))+dotLong(i->first.exponent,omega.subvector(1,omega.size()));
  return omega[0]*(i->second.pAdicValuation(prime))+smartDot(i->first.exponent,omega);
}

static int64 omegaDegree(Term const t, int prime, IntegerVector const &omega)
{
//  return omega[0]*t.c.pAdicValuation(prime)+dotLong(t.m.exponent,omega.subvector(1,omega.size()));
  return omega[0]*t.c.pAdicValuation(prime)+smartDot(t.m.exponent,omega);
}

Term iteratorToTerm(TermMap::const_iterator i)
{
  return Term(i->second,i->first);
}

Polynomial pAdicInitialForm(PolynomialRing const &ZModPZRing, Polynomial const &f, int prime, IntegerVector const &omega)
{
  if(f.isZero())return ZModPZRing.zero();

  Polynomial ret(ZModPZRing);
  int64 d=omegaDegree(f.terms.begin(),prime,omega);

//  debug<<"INPUT"<<omega<<f;

  for(TermMap::const_iterator i=f.terms.begin();i!=f.terms.end();i++)
    {
      if(omegaDegree(i,prime,omega)==d)
        {
        ret+=Term(i->second.pAdicRemainder(ZModPZRing.getField()),Monomial(ZModPZRing,i->first.exponent));
        }
      if(omegaDegree(i,prime,omega)<d)
        {
          ret=Term(i->second.pAdicRemainder(ZModPZRing.getField()),Monomial(ZModPZRing,i->first.exponent));
          d=omegaDegree(i,prime,omega);
          }
    }
//  debug<<"OUTPUT"<<ret;

  return ret;
}

PolynomialSet pAdicInitialForms(PolynomialRing const &ZModPZRing, PolynomialSet const &g, int prime, IntegerVector const &omega)
{
  PolynomialSet ret(g.getRing());

  for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
    ret.push_back(pAdicInitialForm(ZModPZRing,*i,prime,omega));

  return ret;
}


static bool GReater(Term const &a, Term const &b, int prime, IntegerVector const &omega, TermOrder const &tieBreaker)
{
  int aD=omegaDegree(a,prime,omega);
  int bD=omegaDegree(b,prime,omega);
  if (aD>bD) return true;
  if (bD>aD) return false;
  return tieBreaker(b.m.exponent,a.m.exponent);
}


static Term ST(Polynomial const &f, int prime, IntegerVector const &omega, TermOrder const &tieBreaker)
{
  assert(!f.isZero());

  Term ret=iteratorToTerm(f.terms.begin());
  for(TermMap::const_iterator i=f.terms.begin();i!=f.terms.end();i++)
    {
      if(GReater(ret,iteratorToTerm(i),prime,omega,tieBreaker))
        {
          ret=iteratorToTerm(i);
        }
    }
  return ret;
}

static Monomial SM(Polynomial const &f, int prime, IntegerVector const &omega, TermOrder const &tieBreaker)
{
  assert(!f.isZero());
  return ST(f,prime,omega,tieBreaker).m;
}

Polynomial pAdicInitialTerm(PolynomialRing const &ZModPZRing, Polynomial const &f, int prime, IntegerVector const &omega, TermOrder const &tieBreaker)
{
  assert(!f.isZero());

  Term T=ST(f,prime,omega,tieBreaker);

  return Polynomial(Term(T.c.pAdicRemainder(ZModPZRing.getField()),Monomial(ZModPZRing,T.m.exponent)));
}

PolynomialSet pAdicInitialTerms(PolynomialRing const &ZModPZRing, PolynomialSet const &g, int prime, IntegerVector const &omega, TermOrder const &tieBreaker)
{
  PolynomialSet ret(g.getRing());

  for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
    ret.push_back(pAdicInitialTerm(ZModPZRing,*i,prime,omega,tieBreaker));

  return ret;
}

/***
 * Omega is considered to be perturbed.
 */
IntegerVectorList normalPolyhedralInequalities(Polynomial const &f, int prime, IntegerVector const &omega, TermOrder const &tieBreaker)
{//first coordinate is special
  IntegerVectorList ret;
  if(!f.isZero())
    {
      Term T=ST(f,prime,omega,tieBreaker);
      int Tval=T.c.pAdicValuation(prime);
      IntegerVector Texp=T.m.exponent;
      for(TermMap::const_iterator i=f.terms.begin();i!=f.terms.end();i++)
        {
          int Cval=i->second.pAdicValuation(prime);
          IntegerVector Cexp=i->first.exponent;
          IntegerVector temp(1);
          temp[0]=Cval-Tval;
          ret.push_back(concatenation(temp,Cexp-Texp));
        }
    }
  return ret;
}


IntegerVectorList normalPolyhedralInequalities(PolynomialSet const &g, int prime, IntegerVector const &omega, TermOrder const &tieBreaker)
{
  IntegerVectorList ret;
  for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
    {
      IntegerVectorList ret2=normalPolyhedralInequalities(*i,prime,omega,tieBreaker);
      ret.splice(ret.end(),ret2);
    }
  return ret;
}
static int E(Polynomial const &f, Polynomial const &g)
{
  int ret=0;
  for(TermMap::const_iterator i=f.terms.begin();i!=f.terms.end();i++)
    {
      if(!g.terms.count(i->first))ret++;
    }
  return ret;
}

//Based on Diane's draft:
Polynomial longDivision(Polynomial f, PolynomialSet const &G2, int prime, IntegerVector const &omega, TermOrder const &tieBreaker, PolynomialSet &H2, Polynomial &u2)
{
//  debug<<"DIVISION INPUT"<<f<<G2<<omega<<tieBreaker;

  {
    assert(f.isHomogeneous());
    assert(G2.isHomogeneous());
    for(PolynomialSet::const_iterator i=G2.begin();i!=G2.end();i++)assert(!i->isZero());
  }
  FieldZModPZ residueField(prime);
  PolynomialRing S=G2.getRing();
  FieldElement pFieldElement=S.getField().zHomomorphismImplementation(prime);
  int s=G2.size();
  vector<Polynomial> G;
  for(PolynomialSet::const_iterator k=G2.begin();k!=G2.end();k++)G.push_back(*k);
  vector<Polynomial> H(s,Polynomial(S));
  Polynomial r(S);
  Polynomial u=S.one();
  vector<vector<Polynomial> > HCache;
  vector<Polynomial> rCache;
  vector<Polynomial> uCache;
  Polynomial q=f;
  vector<Polynomial> T=G;

  while(!q.isZero())
    {
      DEBUG debug<<"looping\n";
      bool doesSomeSMgDivide=false;
      for(int ii=0;ii<T.size();ii++)
        {
//          debug<<"T[ii]"<<T[ii]<<"\n";
 //         debug<<"q"<<q<<"\n";
          if(SM(T[ii],prime,omega,tieBreaker).exponent.divides(SM(q,prime,omega,tieBreaker).exponent))
              {
                doesSomeSMgDivide=true;
                break;
              }
        }
      if(!doesSomeSMgDivide)
          {
            DEBUG debug<<"a-----\n";
            r+=ST(q,prime,omega,tieBreaker);
            q-=ST(q,prime,omega,tieBreaker);
            DEBUG debug<<"-----a\n";
          }
      else
        {
          DEBUG debug<<"b----\n";
          DEBUG debug<<"i----\n";
          int bestE=0;
          int bestIndex=-1;
          for(int i=0;i<T.size();i++)
            {
              if(SM(T[i],prime,omega,tieBreaker).exponent.divides(SM(q,prime,omega,tieBreaker).exponent))
                  {
                    if((bestIndex==-1) || (E(T[i],q)<bestE))
                        {
                          bestE=E(T[i],q);
                          bestIndex=i;
                        }
                  }
            }
          assert(bestIndex!=-1);
          Polynomial g=T[bestIndex];
          DEBUG debug<<"----i\n";
          //ii
          DEBUG debug<<"ii----\n";
          if(bestE>0)
            {
              T.push_back(q);
              uCache.push_back(u);
              HCache.push_back(H);
              rCache.push_back(r);
            }
          DEBUG debug<<"----ii\n";
          //iii
          DEBUG debug<<"iii----\n";
          Term cuxu=ST(q,prime,omega,tieBreaker)/ST(g,prime,omega,tieBreaker);
          DEBUG debug<<"g="<< g << "   ST   "<< ST(g,prime,omega,tieBreaker)<<"\n";
          DEBUG debug<<"qj="<< q << "   ST   "<< ST(q,prime,omega,tieBreaker)<<"\n";
          Polynomial a=cuxu/SM(cuxu,prime,omega,tieBreaker);
          Polynomial b=S.one();
          int aval=cuxu.c.pAdicValuation(prime);
          while(aval<0)
            {
              b*=pFieldElement;
              a*=pFieldElement;
              aval++;
            }
          Polynomial p=b*q-a*SM(cuxu,prime,omega,tieBreaker)*g;
          if(p.isZero()){DEBUG debug<<"p IS ZERO!\n";}
          else
            {
              DEBUG debug<<"p="<< p << "   ST   "<< ST(p,prime,omega,tieBreaker)<<"\n";
              DEBUG debug<<"bq="<< b*q << "   ST   "<< ST(b*q,prime,omega,tieBreaker)<<"\n";
              assert(GReater(ST(p,prime,omega,tieBreaker),ST(b*q,prime,omega,tieBreaker),prime,omega,tieBreaker));
            }
          DEBUG debug<<"----iii\n";
          //iv
          DEBUG debug<<"iv----\n";
          for(int iprime=0;iprime<s;iprime++)
            if((G[iprime]-g).isZero())
              {
                q=p;
                u=b*u;
                for(int i=0;i<H.size();i++)
                  if(i==iprime)
                    H[i]=b*H[i]+a*SM(cuxu,prime,omega,tieBreaker);
                  else
                    H[i]=b*H[i];
                r=b*r;
                break;
              }
          DEBUG debug<<"----iv\n";
          //v
          DEBUG debug<<"v----\n";
          for(int m=0;m<HCache.size();m++)
              if((g-T[s+m]).isZero())
                {
                  q=p;
                  u=b*u-a*uCache[m];
                  for(int i=0;i<s;i++)
                    H[i]=b*H[i]-a*HCache[m][i];
                  r=b*r-a*rCache[m];
                  break;
                }
          DEBUG debug<<"----v\n";
          DEBUG debug<<"----b\n";
        }
    }
  DEBUG debug<<"Loop end\n";

  H2=PolynomialSet(S);
  for(vector<Polynomial>::const_iterator k=H.begin();k!=H.end();k++)H2.push_back(*k);
  u2=u;

  //Check correctness of result
  {
    Polynomial sum=u*f-r;
    PolynomialSet::const_iterator i=H2.begin();
    for(PolynomialSet::const_iterator j=G2.begin();j!=G2.end();i++,j++)sum-=(*i)*(*j);
    assert(sum.isZero());
    {for(PolynomialSet::const_iterator i=H2.begin();i!=H2.end();i++)assert(i->isHomogeneous());}
    assert(r.isHomogeneous());
    if(!(u*f).isZero())
      {
        if(!r.isZero())
          assert(!GReater(ST(u*f,prime,omega,tieBreaker),ST(r,prime,omega,tieBreaker),prime,omega,tieBreaker));
        i=H2.begin();
        for(PolynomialSet::const_iterator j=G2.begin();j!=G2.end();i++,j++)
          {
            if(!((*i)*(*j)).isZero())
              assert(!GReater(ST(u*f,prime,omega,tieBreaker),ST((*i)*(*j),prime,omega,tieBreaker),prime,omega,tieBreaker));
          }
      }
    DEBUG debug<<"1\n";
    for(TermMap::const_iterator k=r.terms.begin();k!=r.terms.end();k++)
      for(PolynomialSet::const_iterator j=G2.begin();j!=G2.end();j++)
        if(!j->isZero())
          assert(!(SM(*j,prime,omega,tieBreaker).exponent.divides(k->first.exponent)));
  }

  return r;
}

Polynomial SPolynomial(Polynomial const &a, Polynomial const &b, int prime, IntegerVector const &omega, TermOrder const &tieBreaker)
{
  Term lcm=Term(a.getRing().getField().zHomomorphism(1),Monomial(a.getRing(),max(
              ST(a,prime,omega,tieBreaker).m.exponent,
              ST(b,prime,omega,tieBreaker).m.exponent)));

  return (lcm/ST(a,prime,omega,tieBreaker))*a-(lcm/ST(b,prime,omega,tieBreaker))*b;
}

void pAdicAutoReduce(PolynomialSet &g, int prime, IntegerVector const &omega, TermOrder const &tieBreaker)
{
  PolynomialSet temp(g.getRing());
  for(PolynomialSet::iterator i=g.begin();i!=g.end();i++)
    {
      PolynomialSet::iterator ipp=i;ipp++;
      temp.splice(temp.begin(),g,i,ipp);
      i=ipp;

      DEBUG debug<<"TEMP"<<temp;
      DEBUG debug<<"G"<<g;


      PolynomialSet H2(g.getRing());
      Polynomial u2(g.getRing());

      DEBUG debug<<"REDUCING:\n"<<*temp.begin()<<"\n";
      *temp.begin()=longDivision(*temp.begin(),g,prime,omega,tieBreaker,H2,u2);
      DEBUG debug<<"REMAINDER:\n"<<*temp.begin()<<"\n";
//      debug<<"U:\n"<<u2<<"\n";
      if(!(temp.begin()->isZero()))
        {
          (*temp.begin())*=ST(*temp.begin(),prime,omega,tieBreaker).c.inverse();
          g.splice(i,temp,temp.begin(),temp.end());
        }
      else
        {
          temp.pop_back();
        }
//      debug<<g;
      i--;
    }
}


void pAdicMark(Polynomial &g, int prime, IntegerVector const &omega, TermOrder const &tieBreaker)
{
  g.mark(SM(g,prime,omega,tieBreaker));
}

void pAdicMark(PolynomialSet &g, int prime, IntegerVector const &omega, TermOrder const &tieBreaker)
{
  for(PolynomialSet::iterator i=g.begin();i!=g.end();i++)
    i->mark(SM(*i,prime,omega,tieBreaker));
}

void pAdicBuchberger(PolynomialSet &g, int prime, IntegerVector const &omega, TermOrder const &tieBreaker)
{
  for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
    for(PolynomialSet::const_iterator j=g.begin();j!=i;j++)
      {
//        debug<<"checking s poly\n"<<*i<<"("<<ST(*i,prime,omega,tieBreaker)<<")"<<"\n"<<*j<<"("<<ST(*j,prime,omega,tieBreaker)<<")"<<"\n";
        Polynomial f=SPolynomial(*i,*j,prime,omega,tieBreaker);
//        debug<<"S="<<f<<"\n";
// Enable the relatively prime criterion by uncommenting this following:
        /*        if(relativelyPrime(ST(*i,prime,omega,tieBreaker).m.exponent,ST(*j,prime,omega,tieBreaker).m.exponent))
          {
            debug<<"SKIPPING\n";
          }
        else
*/          {

        PolynomialSet H(g.getRing());
        Polynomial u(g.getRing());
        f=longDivision(f,g,prime,omega,tieBreaker,H,u);
        if(!f.isZero())
          {
//            debug<<"Adding:"<<f<<"\n";
            g.push_back(f);
          }
        }
      }
  pAdicAutoReduce(g,prime,omega,tieBreaker);
  pAdicMark(g,prime,omega,tieBreaker);
}

void PAdicGroebnerFanTraverser::updatePolyhedralCone(IntegerVector const &omega, TermOrder const &tieBreaker)
{
  IntegerVectorList inequalities=normalPolyhedralInequalities(groebnerBasis,prime,omega,tieBreaker);
  inequalities.push_back(IntegerVector::standardVector(n+1,0));
  IntegerVectorList empty;

  theCone=PolyhedralCone(inequalities,empty,n+1);
  theCone.canonicalize();
//  debug<<theCone;
}

PAdicGroebnerFanTraverser::PAdicGroebnerFanTraverser(PolynomialSet const &generators, int prime_):
  n(generators.getRing().getNumberOfVariables()),
  groebnerBasis(generators),
  prime(prime_)
{
  IntegerVector omega(n+1);omega[0]=1;
  LexicographicTermOrder tieBreaker;
  pAdicBuchberger(groebnerBasis,prime,omega,tieBreaker);
  updatePolyhedralCone(omega,tieBreaker);
}

static void print(PolynomialSet const &p, int prime, IntegerVector const &omega, TermOrder const &tieBreaker)
{
  for(PolynomialSet::const_iterator i=p.begin();i!=p.end();i++)
    {
      debug<<"("<<ST(*i,prime,omega,tieBreaker)<<")"<<*i<<"\n";
    }

}

void PAdicGroebnerFanTraverser::changeCone(IntegerVector const &ridgeVector, IntegerVector const &rayVector)
{
  IntegerVector omega=ridgeVector;
  WeightTermOrder tieBreaker(rayVector.subvector(1,rayVector.size()));
  WeightTermOrder tieBreaker2(-rayVector.subvector(1,rayVector.size()));

//  debug<<"OLD:\n";
//  print(groebnerBasis,prime,omega,tieBreaker2);

  {
    PolynomialRing R2=residuePolynomialRing(groebnerBasis.getRing(),prime);
    PolynomialSet wallIdeal=pAdicInitialForms(R2,groebnerBasis,prime,omega);
    WeightTermOrder T1(rayVector.subvector(1,rayVector.size()));
//    WeightTermOrder T1(rayVector.subvector(1,rayVector.size()));
    wallIdeal.markAndScale(T1);
//    debug<<"Wall ideal:\n"<<wallIdeal<<"\n";

    PolynomialSet wall2=flip(wallIdeal,rayVector.subvector(1,rayVector.size()));
//    debug<<"Flipped basis of wall ideal:\n"<<wall2<<"\n";
//    WeightTermOrder T2(rayVector.subvector(1,rayVector.size()));

    PolynomialSet liftedBasis(groebnerBasis.getRing());
    for(PolynomialSet::const_iterator i=wall2.begin();i!=wall2.end();i++)
      {
        Polynomial inOldRing=groebnerBasis.getRing().zero();
        for(TermMap::const_iterator j=i->terms.begin();j!=i->terms.end();j++)
          {
            inOldRing+=Term(groebnerBasis.getRing().getField().zHomomorphism(j->second.integerRepresentative()),Monomial(groebnerBasis.getRing(),j->first.exponent));
          }
        PolynomialSet H(groebnerBasis.getRing());
        Polynomial u=groebnerBasis.getRing().zero();
        Polynomial r=longDivision(inOldRing, groebnerBasis,prime,omega,tieBreaker2,H,u);
        liftedBasis.push_back(u*inOldRing-r);
      }
    pAdicAutoReduce(liftedBasis,prime,omega,tieBreaker);
    PolynomialSet newInitialForms=pAdicInitialForms(R2,liftedBasis,prime,omega);
//    debug<<"LIFTED"<<liftedBasis;
//    debug<<"NEWINITIALFORMS"<<newInitialForms;

//  }
//  pAdicBuchberger(groebnerBasis,prime,omega,tieBreaker);
    groebnerBasis=liftedBasis;
  }



//  debug<<"NEW:\n";
//  print(groebnerBasis,prime,omega,tieBreaker);

  updatePolyhedralCone(omega,tieBreaker);
}

IntegerVectorList PAdicGroebnerFanTraverser::link(IntegerVector const &ridgeVector)
{
  IntegerVectorList ret;
/*  IntegerVectorList temp;temp.push_back(IntegerVector::standardVector(n+1,0));
  IntegerVectorList empty;
  PolyhedralCone tempCone=intersection(theCone.link(ridgeVector),PolyhedralCone(empty,temp,n+1));
  tempCone.canonicalize();

  IntegerVector v=tempCone.getUniquePoint();
*/
  IntegerVector v=theCone.link(ridgeVector).getUniquePoint();
//  debug<<"VVVVVVVVVVVVVVVVVVVV"<<v<<ridgeVector<<"\n";
//  debug<<theCone;
//  debug<<theCone.extremeRays();

  ret.push_back(v);

  if(ridgeVector[0]>0)
    {
//      debug<<ridgeVector;
//      assert(ridgeVector[0]==1);
      ret.push_back(-v);
    }
//  debug<<"LINK AT"<<ridgeVector<<"\n"<<ret<<"\n";

  return ret;
}


PolyhedralCone & PAdicGroebnerFanTraverser::refToPolyhedralCone()
{
  return theCone;
}
