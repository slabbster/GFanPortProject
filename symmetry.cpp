#include "symmetry.h"

#include "printer.h"
#include "log.h"
#include <iostream>

#include <map>
using namespace std;

class Trie
{
	class TrieNode
	{
		typedef map<int,class TrieNode> Map;
		Map m;
	public:
		TrieNode()
		{

		}
		TrieNode(IntegerVector const &v, int i)
		{
			if(i<v.size())
			m[v[i]]=TrieNode(v,i+1);
		}
		int stabilizerSize(IntegerVector const &v, int i)const
		{
		  int ret=0;
		  if(i==v.size())return 1;
                  for(Map::const_iterator j=m.begin();j!=m.end();j++)
                    {
                      if(v[i]==v[j->first])
                        ret+=j->second.stabilizerSize(v,i+1);
                    }
                  return ret;
		}
		void search(IntegerVector const &v, IntegerVector  &building, IntegerVector &tempPerm, IntegerVector &ret, IntegerVector &optimal, int i, bool &isImproving)const
		{
			if(i==v.size()){ret=tempPerm;optimal=building;isImproving=false;return;}
			if(isImproving)
				building[i]=-0x7fffffff;
			else
				building[i]=optimal[i];
			for(Map::const_iterator j=m.begin();j!=m.end();j++)
				if(v[j->first]>building[i])
				{
					isImproving=true;
					building[i]=v[j->first];
				}
				for(Map::const_iterator j=m.begin();j!=m.end();j++)
				if(v[j->first]==building[i])
				{
					tempPerm[i]=j->first;
					j->second.search(v,building,tempPerm,ret,optimal,i+1,isImproving);
				}
		}
		void searchStabalizer(IntegerVector const &v, IntegerVector  &building, IntegerVector &tempPerm, IntegerVector &ret, IntegerVector &optimal, int i, bool &isImproving, IntegerVector const &toBeFixed)const
		{
			if(i==v.size())
				if(!(SymmetryGroup::compose(tempPerm,v)<optimal))
					{
					ret=tempPerm;
					optimal=SymmetryGroup::compose(tempPerm,v);
					return;
					}
				for(Map::const_iterator j=m.begin();j!=m.end();j++)
					if(toBeFixed[i]==toBeFixed[j->first])
				{
					tempPerm[i]=j->first;
					j->second.searchStabalizer(v,building,tempPerm,ret,optimal,i+1,isImproving,toBeFixed);
				}
		}
/* this code contains mistakes		void searchStabalizer(IntegerVector const &v, IntegerVector  &building, IntegerVector &tempPerm, IntegerVector &ret, IntegerVector &optimal, int i, bool &isImproving, IntegerVector const &toBeFixed)const
		{
			if(i==v.size()){ret=tempPerm;optimal=building;isImproving=false;debug<<"DEEP";return;}
			if(isImproving)
				building[i]=-0x7fffffff;
			else
				building[i]=optimal[i];
			for(Map::const_iterator j=m.begin();j!=m.end();j++)
				if(toBeFixed[i]==toBeFixed[j->first])
				if(v[j->first]>building[i])
				{
					isImproving=true;
					building[i]=v[j->first];
				}
				for(Map::const_iterator j=m.begin();j!=m.end();j++)
					if(toBeFixed[i]==toBeFixed[j->first])
				if(v[j->first]==building[i])
				{
					debug.printInteger(i);debug<<":";
					debug.printInteger(j->first);debug<<" ";
					tempPerm[i]=j->first;
					j->second.searchStabalizer(v,building,tempPerm,ret,optimal,i+1,isImproving,toBeFixed);
				}
		}*/
	//	void doubleSearch();
		void insert(IntegerVector const &v, int i)
		{
			if(i==v.size())return;
			if(m.count(v[i]))
				m[v[i]].insert(v,i+1);
			else
				m[v[i]]=		TrieNode(v,i+1);
		}
		void print(int i, int n)const
		{
			if(i==n)return;
			for(Map::const_iterator j=m.begin();j!=m.end();j++)
			{
				{for(int j=0;j<2*i;j++)debug<<" ";}
				debug.printInteger(j->first);
				debug<<"\n";
				j->second.print(i+1,n);
			}
			}
		int size(int i,int n)const
		{
			if(i==n)return 1;
			int ret=0;
			for(Map::const_iterator j=m.begin();j!=m.end();j++)
				ret+=j->second.size(i+1,n);
			return ret;
		}
	};
public:
	TrieNode theTree;
	int n;
	Trie(int n_):
		n(n_),
		theTree(SymmetryGroup::identity(n_),0)
	{
	}
	int size()const
	{
		return theTree.size(0,n);
	}
	void insert(IntegerVector const &v)
	{
		theTree.insert(v,0);
//		debug<<v;
//		theTree.print(0,v.size());

//		debug<<"---------------------------------------------\n";
	}
	/**
	 * returns the sigma from the set with sigma(v) maximal in the lexicographic ordering.
	 */
	IntegerVector search(IntegerVector const &v)
	{
		IntegerVector tempPerm(v.size());
		IntegerVector ret(v.size());
		IntegerVector building(v.size());
		IntegerVector optimal=v;//the identity is always in the trie
		bool isImproving=true;
		theTree.search(v,building,tempPerm,ret,optimal,0,isImproving);
		return ret;
	}
	IntegerVector searchStabalizer(IntegerVector const &v, IntegerVector const &toBeFixed)
	{
		IntegerVector tempPerm=SymmetryGroup::identity(v.size());
		IntegerVector ret(v.size());
		IntegerVector building(v.size());
		IntegerVector optimal=v;//the identity is always in the trie
		bool isImproving=true;
		theTree.searchStabalizer(v,building,tempPerm,ret,optimal,0,isImproving,toBeFixed);
		return ret;
	}
        int stabilizerSize(IntegerVector const &v)const
        {
          return theTree.stabilizerSize(v,0);
        }
};







IntegerVector SymmetryGroup::identity(int n)
{
  IntegerVector v(n);
  for(int i=0;i<n;i++)v[i]=i;

  return v;
}


IntegerVector SymmetryGroup::inverse(IntegerVector const &a)
{
  return composeInverse(a,identity(a.size()));
}


SymmetryGroup::SymmetryGroup(int n):
  byteTable(0),
  trie(0)
{
  elements.insert(identity(n));
}


int SymmetryGroup::sizeOfBaseSet()const
{
  assert(!elements.empty());
  return elements.begin()->size();
}

void SymmetryGroup::computeClosure(IntegerVector const &v) //does this work??
{
  ElementContainer newOnes;

  newOnes.insert(v);

  while(!newOnes.empty())
    {
      static int i;
      i++;
      if((i&127)==0)fprintf(Stderr,"%i\n",i);


      IntegerVector v=*newOnes.begin();
      for(ElementContainer::const_iterator i=elements.begin();i!=elements.end();i++)
	{
	  {
	    IntegerVector n(compose(*i,v));
	    if(0==elements.count(n))
	      newOnes.insert(n);
	  }
	  {
	    IntegerVector n(compose(v,*i));
	    if(0==elements.count(n))
	      newOnes.insert(n);
	  }
	}
      newOnes.erase(v);
      elements.insert(v);
    }
}


void SymmetryGroup::computeClosure(IntegerVectorList const &l)
{
  //  for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
  //  computeClosure(*i);

  bool growing=true;
  while(growing)
    {
      growing=false;
      for(ElementContainer::const_iterator i=elements.begin();i!=elements.end();i++)
	{
	  for(IntegerVectorList::const_iterator j=l.begin();j!=l.end();j++)
	    {
	      {
		IntegerVector n(compose(*i,*j));
		growing|=(0==elements.count(n));
		elements.insert(n);
	      }
	      {
		IntegerVector n(compose(*i,*j));
		growing|=(0==elements.count(n));
		elements.insert(n);
	      }
	    }
	}
    }
}

IntegerVectorList SymmetryGroup::getUniqueGenerators()const
{
  int n=sizeOfBaseSet();
  IntegerVectorList ret;

restart:
  SymmetryGroup temp(n);
  temp.computeClosure(ret);
  ElementContainer::const_iterator j=temp.elements.begin();
  for(ElementContainer::const_iterator i=elements.begin();i!=elements.end();i++,j++)
    {
      if(*j!=*i)
        {
          ret.push_back(*i);
          goto restart;
        }
    }
  return ret;
}

void SymmetryGroup::print(FILE *f)
{
  AsciiPrinter P(f);
  P.printString("Printing SymmetryGroup\n");
  IntegerVectorList l;
  for(ElementContainer::const_iterator i=elements.begin();i!=elements.end();i++)
    {
      //      P.printVector(*i);
      //      P.printNewLine();
      l.push_back(*i);
    }
  P.printVectorList(l);
  fprintf(f,"Group order=%i\n",elements.size());
  P.printString("Done printing SymmetryGroup.\n");
}


IntegerVector SymmetryGroup::compose(IntegerVector const &perm, IntegerVector const &b)
{
  IntegerVector v(perm);
  assert(perm.size()==b.size());
  for(int i=0;i<perm.size();i++)v[i]=b[perm[i]];
  return v;
}


IntegerVector SymmetryGroup::composeInverse(IntegerVector const &a, IntegerVector const &b)
{
  IntegerVector v(a);
  assert(a.size()==b.size());
  for(int i=0;i<a.size();i++)v[a[i]]=b[i];
  return v;
}


IntegerVectorList SymmetryGroup::permuteIntegerVectorList(IntegerVectorList const &l, IntegerVector const &v)
{
  IntegerVectorList ret;
  for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
    ret.push_back(compose(v,*i));

  return ret;
}


Polynomial SymmetryGroup::permutePolynomial(Polynomial const &p, IntegerVector const &v)
{
  Polynomial q(p.getRing());

  for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
    {
      q+=Term(i->second,Monomial(p.getRing(),compose(v,i->first.exponent)));
    }

  q.mark(Monomial(p.getRing(),compose(v,p.getMarked().m.exponent)));

  return q;
}


PolynomialSet SymmetryGroup::permutePolynomialSet(PolynomialSet const &s, IntegerVector const &v)
{
  PolynomialRing theRing=s.getRing();
  PolynomialSet ret(theRing);
  for(PolynomialSet::const_iterator i=s.begin();i!=s.end();i++)
    {
      ret.push_back(permutePolynomial(*i,v));
    }

  return ret;
}


Polynomial SymmetryGroup::computeUniqueRepresentative(Polynomial p)
{
  Polynomial best=p;

  for(ElementContainer::const_iterator i=elements.begin();i!=elements.end();i++)
    {
      Polynomial q=permutePolynomial(p,*i);
      if(PolynomialCompare()(best,q))best=q;
    }
  return best;
}


IntegerVector SymmetryGroup::orbitRepresentative(IntegerVector const &v, IntegerVector *usedPermutation)const
{
	if(trie){
		  if(usedPermutation)
			  {
			  *usedPermutation=trie->search(v);
				return compose(*usedPermutation,v);
			  }
		return compose(trie->search(v),v);
	}
  IntegerVector ret=v;
  ElementContainer::const_iterator usedPerm;
  for(ElementContainer::const_iterator i=elements.begin();i!=elements.end();i++)
    {
      IntegerVector q=compose(*i,v);
      if(! (q<ret))//negation to make sure that usedPerm is set
	{
	  usedPerm=i;
	  ret=q;
	}
    }

  if(usedPermutation)*usedPermutation=*usedPerm;

  if(trie)
  {
//	  debug<<"Input"<<v<<"\n";
//	  debug<<"Bruteforce"<<ret<<"\n";
	  IntegerVector triePerm=trie->search(v);
//	  debug<<"Trie"<<compose(triePerm,v)<<"\n";
	  assert((compose(triePerm,v)-ret).isZero());
  }

  return ret;
}

IntegerVector SymmetryGroup::orbitRepresentativeFixing(IntegerVector const &v, IntegerVector const &fixed)const
{
	if(trie){
		return compose(trie->searchStabalizer(v,fixed),v);
	}
  IntegerVector ret=v;

  for(ElementContainer::const_iterator i=elements.begin();i!=elements.end();i++)
    if(compose(*i,fixed)==fixed)
      {
	IntegerVector q=compose(*i,v);
	if(ret<q)ret=q;
      }
	if(trie){
		IntegerVector temp=compose(trie->searchStabalizer(v,fixed),v);
//		debug<<"Input"<<v;
//		debug<<"Brute"<<ret;
//		debug<<"Quick"<<temp;
		assert((temp-ret).isZero());
//		return compose(trie->searchStabalizer(v,fixed),v);
	}
  return ret;
}

bool SymmetryGroup::isPermutation(IntegerVector const &a)
{
  int n=a.size();
  IntegerVector temp(n);
  for(int i=0;i<n;i++)temp[i]=-1;
  for(int i=0;i<n;i++)
    {
      if(a[i]<0 || a[i]>=n)return false;
      temp[i]=i;
    }
  for(int i=0;i<n;i++)if(temp[i]<0)return false;
  return true;
}


IntegerVector SymmetryGroup::combinePermutationAndSignChanges(IntegerVector const &permutation, IntegerVector const &signChanges)
{
	assert(isPermutation(permutation));
	int n=permutation.size();
	assert(n==signChanges.size());
	IntegerVector ret(2*n);
	for(int i=0;i<n;i++)
		if(signChanges[i]==1)
		{
			ret[i]=permutation[i];
			ret[i+n]=n+permutation[i];
		}
		else
		{
			ret[i]=n+permutation[i];
			ret[i+n]=permutation[i];
		}
	return ret;
}

void SymmetryGroup::extractPermuationAndSignChanges(IntegerVector const &v, IntegerVector &permutation, IntegerVector &signChanges)
{
	int n=v.size()/2;
	permutation=IntegerVector(n);
	signChanges=IntegerVector(n);

	for(int i=0;i<n;i++)
	{
		permutation[i]=v[i]%n;
		if(v[i]<n)
			signChanges[i]=1;
		else
			signChanges[i]=-1;
	}
}


int SymmetryGroup::orbitSize(IntegerVector const &stable)const
{
  int groupSize=elements.size();

  int n=stable.size();
  int numFixed=0;

  if(trie)
    {
      numFixed=trie->stabilizerSize(stable);
    }
  else
    {
      for(SymmetryGroup::ElementContainer::const_iterator j=elements.begin();j!=elements.end();j++)
        {
          bool doesFix=true;

          for(int i=0;i<n;i++)
            if(stable[i]!=stable[(*j)[i]])
              {
                doesFix=false;
                break;
              }
          if(doesFix)numFixed++;
        }
    }
  return groupSize/numFixed;
}


IntegerVector SymmetryGroup::fundamentalDomainInequality(IntegerVector const &perm)
{
  for(int i=0;i<perm.size();i++)
    if(perm[i]!=i)
      return IntegerVector::standardVector(perm.size(),i)-IntegerVector::standardVector(perm.size(),perm[i]);
      //      return -IntegerVector::standardVector(perm.size(),i)-IntegerVector::standardVector(perm.size(),perm[i]);//USE this for 4x4 of 5x5
  return IntegerVector(perm.size());
}


IntegerVectorList SymmetryGroup::fundamentalDomainInequalities()const
{
  set<IntegerVector> ret2;
  for(ElementContainer::const_iterator i=elements.begin();i!=elements.end();i++)
    ret2.insert(fundamentalDomainInequality(*i));

  IntegerVectorList ret;
  for(set<IntegerVector>::const_iterator i=ret2.begin();i!=ret2.end();i++)
    if(!i->isZero())ret.push_back(*i);

  return ret;
}


void SymmetryGroup::createByteTable()
{
  assert(!byteTable);
  int n=sizeOfBaseSet();
  byteTableHeight=elements.size();
  byteTable=(unsigned char*)malloc(n*byteTableHeight*sizeof(unsigned char));
  assert(byteTable);
  int j=0;
  for(ElementContainer::const_iterator i=elements.begin();i!=elements.end();i++,j++)
    {
      for(int k=0;k<n;k++)
	byteTable[j*n+k]=(*i)[k];
    }
}

void SymmetryGroup::createTrie()
{
	log1 debug<<"Creating symmetry trie.\n";
	trie=new Trie(sizeOfBaseSet());
	for(ElementContainer::const_iterator i=elements.begin();i!=elements.end();i++)
		trie->insert(*i);
	log2 debug<<"Number of elements";log2 debug.printInteger(trie->size());log2 debug<<"\n";
	log1 debug<<"Done creating symmetry trie.\n";
}


unsigned char * SymmetryGroup::getByteTable()const
{
  return byteTable;
}


int  SymmetryGroup::getByteTableHeight()const
{
  return byteTableHeight;
}


bool SymmetryGroup::isTrivial()const
{
	ElementContainer::const_iterator i=elements.begin();
	assert(i!=elements.end());
	i++;
	return i==elements.end();
}

static int mergeSortRek(IntegerVector &v, int begin, int end, IntegerVector &temp)
{
  if(end-begin<2)return 0;
  int med=(begin+end)>>1;
  int nswaps=mergeSortRek(v,begin,med,temp);
  nswaps+=mergeSortRek(v,med,end,temp);

  {
    int Astart=begin;
    int Alength=med-begin;
    int Bstart=med;
    int Blength=end-med;
    int nextFree=begin;
    while(nextFree!=end)
      {
//        debug<<"Astart:"<<Astart<<"Alength:"<<Alength<<"Bstart:"<<Bstart<<"Blength:"<<Blength<<"nextFree:"<<nextFree<<"\n";
        if(Blength==0 || (Alength!=0 && v[Astart]<v[Bstart]))
          {
            temp[nextFree++]=v[Astart++];
            Alength--;
          }
        else
          {
            temp[nextFree++]=v[Bstart++];
            nswaps+=Alength;
            Blength--;
          }
      }
    for(int i=begin;i!=end;i++)v[i]=temp[i];
  }
//debug<<"return\n";
  return nswaps;
}

int mergeSort(IntegerVector &v)
{
  IntegerVector temp(v.size());
  return mergeSortRek(v,0,v.size(),temp);
}

