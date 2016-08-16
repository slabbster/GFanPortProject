#include "determinantpoly.h"

#include  <iostream>
#include "printer.h"
#include "buchberger.h"
#include "wallideal.h"
#include "termorder.h"
#include "tropical2.h"
#include "tropicaldeterminant.h"
#include "log.h"

using namespace std;

static list<int> interval(int n, bool sqrt)
{
        list<int> ret;
        for(int i=0;((sqrt)?i*i:i)<n;i++)ret.push_back(i);
        return ret;
}


class PolyMatrix
{
public:
        class PolyMatrixEntry
        {
        public:
                bool isZero;
                Polynomial p;
                int maxDegree,mminDegree;
                vector<Polynomial> forms;
                PolyMatrixEntry(Polynomial const &p_, IntegerVector const &w):
                        p(p_),
                        isZero(p_.isZero())
                        {
                        if(p.isZero())
                        {
                                maxDegree=td_minusInfinity;
                                mminDegree=td_minusInfinity;
                                return;
                        }
                        maxDegree=p_.degree(w);
                        mminDegree=p_.degree(-w);


        //              debug<<maxDegree<<minDegree<<"\n";
                        for(int i=0;i<=maxDegree+mminDegree;i++)forms.push_back(Polynomial(p.getRing()));
                        for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
                                {
                                        forms[dot(i->first.exponent,w)+mminDegree]+=Term(i->second,i->first);
                                }
                }
        };

        list<int>::iterator indexOfSparsestRow(list<int> &rows, list<int> const &columns, int &chosenRowIteratorIndex)
        {
                list<int>::iterator ret=rows.end();
                int bestNumberOfZeros=-1;
                int I=0;
                for(list<int>::iterator i=rows.begin();i!=rows.end();i++,I++)
                {
                        int numberOfZeros=0;
                        for(list<int>::const_iterator j=columns.begin();j!=columns.end();j++)numberOfZeros+=data[*i][*j].isZero;
                        if(numberOfZeros>bestNumberOfZeros)
                        {
                                bestNumberOfZeros=numberOfZeros;
                                ret=i;
                                chosenRowIteratorIndex=I;
                        }
                }
//              cerr<<":"<<bestNumberOfZeros<<endl;
                return ret;
        }
        vector<vector<PolyMatrixEntry> > data;
        PolynomialRing theRing;
        //constructor for Jacobi matrix
        PolyMatrix(PolynomialSet const &generators, IntegerVector const &w, bool takeDerivatives, int numberOfDVariables=-1):
                theRing(generators.getRing())
        {
          if(takeDerivatives)
            {
              int n=generators.size();
              if(numberOfDVariables==-1)numberOfDVariables=n;
              PolynomialSet::const_iterator I=generators.begin();
              for(int i=0;i<n;i++)
                {
                  vector<PolyMatrixEntry> row;
                  for(int j=0;j<numberOfDVariables;j++)
                    row.push_back(PolyMatrixEntry(I->derivative(j),w));
                  data.push_back(row);
                  I++;
                }
            }
          else
            {
              int n=generators.size();
                for(int i=0;i<n;i++)if(i*i>=n){n=i;break;}assert(n*n==generators.size());
                PolynomialSet::const_iterator I=generators.begin();
                for(int i=0;i<n;i++)
                  {
                        vector<PolyMatrixEntry> row;
                        for(int j=0;j<n;j++)
                        {
                                row.push_back(PolyMatrixEntry(*I,w));
                                I++;
                        }
                        data.push_back(row);
                  }
            }
        }
        void print()
        {
                int n=data.size();
                for(int i=0;i<n;i++)
                {
                        for(int j=0;j<data[i].size();j++)
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
        IntegerMatrix subTropicalMatrix(bool max, list<int> const &rows, list<int> const &columns)
        {
                IntegerMatrix ret(rows.size(),columns.size());
                list<int>::const_iterator I=rows.begin();
                for(int i=0;i<ret.getHeight();i++,I++)
                {
                        list<int>::const_iterator J=columns.begin();
                        for(int j=0;j<ret.getWidth();j++,J++)
                                ret[i][j]=(max)?data[*I][*J].maxDegree:data[*I][*J].mminDegree;
                }

                return ret;
        }

        Polynomial determinantForm(list<int> rows, list<int> columns, int degree, int level)
        {
                Polynomial ret(theRing);

                IntegerMatrix matMax=subTropicalMatrix(true,rows,columns);
                IntegerMatrix matMMin=subTropicalMatrix(false,rows,columns);

/*                if(level==0)
                {
                        debug<<matMax.getRows();
                        debug<<matMMin.getRows();

                        int maxDet=tropicalDeterminant(matMax);
                        int mminDet=tropicalDeterminant(matMMin);

                        debug<<"Max:"<<maxDet<<" Min:"<<-mminDet<<"\n";
                }
*/
//              debug<<matMax.getRows();
//              debug<<matMMin.getRows();

                int maxDet=tropicalDeterminant(matMax);
                int mminDet=tropicalDeterminant(matMMin);

//              assert(degree==0);

//              maxDet=0;//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//              mminDet=0;//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//              debug<<maxDet<<" "<<mminDet<<"\n";
                if(degree>maxDet)return Polynomial(theRing);
                if(-degree>mminDet)return Polynomial(theRing);
                if(rows.size()==1)
                  {
                    if(data[rows.front()][columns.front()].isZero)return Polynomial(theRing);
                   // cerr<<degree<<data[rows.front()][columns.front()].mminDegree<<endl;
                    return data[rows.front()][columns.front()].forms[degree+data[rows.front()][columns.front()].mminDegree];
                  }
                int chosenRowIteratorIndex=0;
                list<int>::iterator chosenRowIterator=indexOfSparsestRow(rows, columns, chosenRowIteratorIndex);//rows.begin();
                int chosenrow=*chosenRowIterator;
                {
                        list<int>::iterator temp=chosenRowIterator;temp++;
                        rows.erase(chosenRowIterator);
                        //no need to update rows afterwards since it is stored on the stack
                }
/*                if(level==0)cerr<<"-"<<chosenrow<<endl;
                if(level==1)cerr<<"++"<<chosenrow<<endl;
                if(level==2)cerr<<"***"<<chosenrow<<endl;
                if(level==3)cerr<<"%%%%%%%%"<<chosenrow<<endl;
  */              int sign=1;
                sign*=(1-2*(chosenRowIteratorIndex&1));
                for(list<int>::iterator i=columns.begin();i!=columns.end();i++)
                {
                        int chosencol=*i;
                        if(!data[chosenrow][chosencol].isZero)
                        {
                                list<int>::iterator temp=i;temp++;
                                columns.erase(i);
                                i=temp;
//                              int D=degree+1;
//                              if(data[chosenrow][chosencol].forms.size()<D)D=data[chosenrow][chosencol].forms.size();
//                              for(int k=0;k<D;k++)
                                for(int k=0;k<=data[chosenrow][chosencol].maxDegree+data[chosenrow][chosencol].mminDegree;k++)
                                {
                                  //cerr<<chosenrow<<chosencol<<endl;
                                  //cerr<<degree<<-k<<data[chosenrow][chosencol].mminDegree<<endl;
                                  assert(k<data[chosenrow][chosencol].forms.size());
                                        if(!data[chosenrow][chosencol].forms[k].isZero())
                                        ret+=data[chosenrow][chosencol].forms[k]*determinantForm(rows,columns,degree-k+data[chosenrow][chosencol].mminDegree,level+1)*theRing.getField().zHomomorphism(sign);
                                }
                                columns.insert(i,chosencol);
                                i--;
                        }
                        sign*=-1;
                }
                rows.push_front(chosenrow);
                return ret;
        }
        Polynomial initialFormOfDeterminant(bool onlyComputeInitial=true)
        {
          list<int> rows=interval(data.size(),0);
          list<int> columns=interval(data.size(),0);

          IntegerMatrix matMax=subTropicalMatrix(true,rows,columns);
          IntegerMatrix matMMin=subTropicalMatrix(false,rows,columns);

          int maxDet=tropicalDeterminant(matMax);
          int mminDet=tropicalDeterminant(matMMin);

          cerr<<"Max"<<maxDet<<"Min"<<-mminDet<<endl;

          Polynomial ret(theRing);
          for(int d=maxDet;d>=-mminDet;d--)
            {
              log1 debug<<"Computing degree "<<d<<"\n";
              ret+=determinantForm(rows,columns,d,0);
              if(onlyComputeInitial&&!ret.isZero())break;
            }
          return ret;
        }
        Polynomial determinant()
        {
          return initialFormOfDeterminant(false);
        }
};


Polynomial initialFormOfDeterminant(PolynomialSet const &g, IntegerVector const &w, bool takeDerivatives)
{
  PolyMatrix m(g,w,takeDerivatives);

  return m.initialFormOfDeterminant();
}

list<int> toIndexed(list<int> const &l)
{
  list<int> ret;
  int a=0;
  for(list<int>::const_iterator i=l.begin();i!=l.end();i++,a++)if(*i)ret.push_back(a);
  return ret;
}

PolynomialSet jacobiMinors(PolynomialSet const &g, int codim)
{
  IntegerVector w(g.getRing().getNumberOfVariables());
  PolyMatrix m(g,w,true,w.size());

  //log1 m.print();
  PolynomialSet ret(g.getRing());

  if(codim<=g.size())if(codim<=w.size())
    {
  list<int> rows;
  for(int i=0;i<g.size()-codim;i++)rows.push_back(0);
  for(int i=0;i<codim;i++)rows.push_back(1);
  do
  {
    list<int> cols;
    for(int i=0;i<w.size()-codim;i++)cols.push_back(0);
    for(int i=0;i<codim;i++)cols.push_back(1);
    do
      {
 //       cerr<<rows.size()<<cols.size()<<endl;
        ret.push_back(m.determinantForm(toIndexed(rows),toIndexed(cols),0,0));
      }
    while(next_permutation(cols.begin(),cols.end()));
  }
  while(next_permutation(rows.begin(),rows.end()));
    }
  return ret;
}
