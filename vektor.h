#ifndef VEKTOR_H_INCLUDED
#define VEKTOR_H_INCLUDED

#include <vector>
#include <list>
#include <assert.h>
#include <algorithm>
#include <complex>
#include <stdio.h>

using namespace std;

typedef signed long long int int64;

void outOfRange(int i, int n);

/*
 * TODO: in the process of making gfan a library to be used by POLYMAKE etc.
 * we need to:
 * - make an integer class with overflow checking. At every addition the difference
 *   between bit 31 and 32 is or'ed to a register, which can be checked later.
 * - group the current asserts into categories. Some will assert in the usual way.
 *   Others will throw exceptions. In the library version all throw exceptions.
 */

template <class typ> class Vektor{
public:
  vector<typ> v;
  int support;
public:
  //--------------
  // Constructors
  //--------------
  Vektor(const Vektor &a):v(a.v),support(a.support){
  }
  Vektor(int n):v(n){
    assert(n>=0);
    support=0;
    for(int i=0;i<n;i++)v[i]=typ(0.0);
  };
  Vektor(const typ *p, int n):v(n){
    assert(n>=0);
    support=0;for(int i=0;i<n;i++)v[i]=p[i];
  };
  ~Vektor(){
  };
  Vektor(){
    support=0;
  };
  static Vektor standardVector(int n, int i)
    {
      Vektor v(n);
      v[i]=typ(1.0);
      return v;
    }
  static Vektor allOnes(int n)
    {
      Vektor v(n);
      for(int i=0;i<n;i++)
	v[i]=typ(1.0);
      return v;
    }
  //-------------
  // Conversions
  //-------------
  template<class T> Vektor(const Vektor<T>& c)
	  :v(c.size())
	{
	  for(int i=0;i<size();i++)v[i]=typ(c[i]);
	}

  //--------
  // Access
  //--------
  typ& operator[](int n)
    {
      if(!(n>=0 && n<v.size()))outOfRange(n,v.size());
      return (v[n]);
    }
  const typ& operator[](int n)const{assert(n>=0 && n<v.size());return (v[n]);}
  const typ& UNCHECKEDACCESS(int n)const{return (v[n]);}

  //-------------
  // STL related
  //-------------
  unsigned int size()const{return v.size();};
  void resize(int n){v.resize(n,0);};
  void grow(int i){if(size()<i)resize(i);}
  void push_back(typ a)
  {
    v.push_back(a);
  }
  void sort()
    {
      std::sort(v.begin(),v.end());
    }
  bool nextPermutation()
    {
      return std::next_permutation(v.begin(),v.end());
    }
  bool operator<(const Vektor & b)const
    {
      if(size()<b.size())return true;
      if(size()>b.size())return false;
      for(int i=0;i<size();i++)
	{
	  if(v[i]<b[i])return true;
	  if(b[i]<v[i])return false;
	}
      return false;
    }

  //-----------------
  // Arithmetic fast
  //-----------------
  typ sum()const{typ f=0;for(int i=0;i<size();i++)f+=v[i];return f;};
  Vektor& operator+=(const Vektor& q){assert(size()==q.size());for(int i=0;i<size();i++)v[i]+=q.v[i];return *this;}
  Vektor& operator-=(const Vektor& q){assert(size()==q.size());for(int i=0;i<size();i++)v[i]-=q.v[i];return *this;}
  inline friend typ dot(const Vektor& p, const Vektor& q){assert(p.size()==q.size());typ s=0;for(int i=0;i<p.size();i++)s+=p[i]*q[i];return s;}
 // inline friend int64 dotLong(const Vektor& p, const Vektor& q){assert(p.size()==q.size());int64 s=0;for(int i=0;i<p.size();i++)s+=(int64)p[i]*(int64)q[i];return s;}
  inline friend int64 dotLong(const Vektor& p, const Vektor& q)
  {assert(p.v.size()==q.v.size());
  int64 s=0;
  typename vector<typ>::const_iterator j=q.v.begin();
  for(typename std::vector<typ>::const_iterator i=p.v.begin();i!=p.v.end();i++,j++)s+=((int64)*i)*((int64)*j);
  return s;}
  bool operator==(const Vektor & q)const{if(size()!=q.size())return false;for(int i=0;i<size();i++)if(v[i]!=q[i])return false;return true;}
  bool operator!=(const Vektor & q)const {return !(operator==(q));}
  bool isZero() const
    {
      int n=v.size();
      for(int i=0;i<n;i++)if(v[i]!=0)return 0;
      return 1;
    }
  bool isPositive() const
    {
      int n=v.size();
      for(int i=0;i<n;i++)if(v[i]<=0)return 0;
      return 1;
    }
  bool isNonNegative() const
    {
      int n=v.size();
      for(int i=0;i<n;i++)if(v[i]<0)return 0;
      return 1;
    }
  int max()const
  {
    int ret=-0x7fffffff; //not completely correct, but kind of works for 64bit
    for(int i=0;i<v.size();i++)if(ret<v[i])ret=v[i];
    return ret;
  }
  int min()const
  {
    int ret=0x7fffffff;
    for(int i=0;i<v.size();i++)if(ret>v[i])ret=v[i];
    return ret;
  }
  friend bool dependent(const Vektor& p, const Vektor& q)
    {
	  /*
      typ pp=dot(p,p);
      typ qq=dot(q,q);
      typ pq=dot(p,q);
      return pq*pq==pp*qq;
*/
	  int n=p.size();
	  assert(n==q.size());
	  int i;
	  for(i=0;i<n;i++)
	  {
		  if(p.v[i])break;
	  }
	  if(i==n)return true;
	  if(q.v[i]==0)return q.isZero();
	  int64 a=p.v[i];
	  int64 b=q.v[i];
	  for(int j=0;j<n;j++)
		  if(a*q.v[j]!=b*p.v[j])return false;
	  return true;
    }

  //-----------------
  // Arithmetic slow
  //-----------------
  inline friend Vektor operator-(const Vektor& q){return -1*q;};
  inline friend Vektor operator*(typ s, const Vektor& q){Vektor p=q;for(int i=0;i<q.size();i++)p[i]*=s;return p;}
  inline friend Vektor operator/(const Vektor& q, typ s){Vektor p=q;for(int i=0;i<q.size();i++)p[i]/=s;return p;}
  inline friend Vektor operator*(const Vektor& p, const Vektor& q){assert(p.size()==q.size());Vektor p1=p;for(int i=0;i<p.size();i++)p1.v[i]*=q.v[i];return p1;}
//  inline friend Vektor operator+(const Vektor& p, const Vektor& q){assert(p.size()==q.size());Vektor p1=p;for(int i=0;i<p.size();i++)p1[i]+=q[i];return p1;}
  inline friend Vektor operator+(const Vektor& p, const Vektor& q){if(p.size()!=q.size()){fprintf(stderr,"%i %i\n",p.size(),q.size());assert(p.size()==q.size());};Vektor p1=p;for(int i=0;i<p.size();i++)p1[i]+=q[i];return p1;}
  inline friend Vektor operator-(const Vektor& p, const Vektor& q){assert(p.size()==q.size());Vektor p1=p;for(int i=0;i<p.size();i++)p1[i]-=q[i];return p1;}
  friend Vektor max(const Vektor& p, const Vektor& q){assert(p.size()==q.size());Vektor p1=p;for(int i=0;i<p.size();i++)if(p1[i]<q[i])p1[i]=q[i];return p1;}
  friend Vektor min(const Vektor& p, const Vektor& q){assert(p.size()==q.size());Vektor p1=p;for(int i=0;i<p.size();i++)if(p1[i]>q[i])p1[i]=q[i];return p1;}

  //------------------
  // Monomial related
  //------------------
  int divides(const Vektor& q) const
    {
      assert(size()==q.size());
      int n=v.size();
      for(int i=0;i<n;i++)
        {
          if(v[i]>0)if(q.v[i]<v[i])return 0;
        }
      return 1;
    }
  inline friend bool relativelyPrime(const Vektor& p, const Vektor& q)
    {
      assert(p.size()==q.size());
      int n=p.size();
      for(int t=0;t<n;t++)if((p[t]>0)&&(q[t]>0)) return false;
      return true;
    }
  Vektor supportVector()const
    {
      Vektor r(v.size());
      for(int i=0;i<size();i++)
	r[i]=(v[i]!=0);
      return r;
    }

  //------------------------------
  // Subvectors and concatenation
  //------------------------------
  Vektor subvector(int begin, int end)const
    {
      assert(begin>=0);
      assert(end<=size());
      assert(end>=begin);
      Vektor ret(end-begin);
      for(int i=0;i<end-begin;i++)
	ret[i]=v[begin+i];
      return ret;
    }
  Vektor subvector(list<int> const &chosenIndices)const
  {
    Vektor ret(chosenIndices.size());
    int I=0;
    for(list<int>::const_iterator i=chosenIndices.begin();i!=chosenIndices.end();i++,I++)ret[I]=v[*i];
    return ret;
  }
  friend Vektor concatenation(Vektor const &a, Vektor const &b)
  {
    Vektor ret(a.size()+b.size());
    for(int i=0;i<a.size();i++)ret[i]=a[i];
    for(int i=0;i<b.size();i++)ret[i+a.size()]=b[i];
    return ret;
  }

  //----------------------
  // Zero/nonZero entries
  //----------------------
  int indexOfLargestNonzeroEntry()const
  {
    int ret=-1;
    for(int i=0;i<v.size();i++)
      {
	if(v[i])ret=i;
      }
    return ret;
  }
  Vektor supportIndices()const
  {
    Vektor ret(0);
    for(int i=0;i<v.size();i++)
      if(v[i]!=0)ret.push_back(i);
    return ret;
  }
  Vektor supportAsZeroOneVector()const
  {
    Vektor ret(v.size());
    for(int i=0;i<v.size();i++)ret[i]=bool(v[i]);
    return ret;
  }
  void calcsupport(void)
    {
      support=0;
      for(int i=0;i<v.size();i++)support=(support<<1)|(((v[i]>0)==true)&1);
    }
};

typedef complex<double> ComplexNumber;

typedef Vektor<ComplexNumber> ComplexVector;
typedef Vektor<double> FloatVector;
typedef Vektor<int> IntegerVector;
typedef list<IntegerVector> IntegerVectorList;

IntegerVectorList transposeIntegerVectorList(IntegerVectorList const &l);
IntegerVectorList multiplyIntegerVectorList(IntegerVectorList const &A, IntegerVectorList const &B);
IntegerVectorList subvectorsOfIntegerVectorList(IntegerVectorList const &l, list<int> const &chosen);
int gcdOfVector(IntegerVector const &v);
void normalizedLowLevel(IntegerVector const &v, IntegerVector &dest);
IntegerVector normalized(IntegerVector const &v);

/**
 * Removes duplicates and reorders list.
 */
void removeDuplicates(IntegerVectorList &l);

#endif



int gcdGFAN(int r, int s);

