#include "symmetriccomplex.h"

#include <sstream>
#include "polyhedralcone.h"
#include "printer.h"
#include "lp.h"
#include "linalg.h"
#include "determinant.h"
#include "log.h"
#include <iostream>

// {ab+aa+cc+dd,ba+ca+db}

SymmetricComplex::Cone::Cone(set<int> const &indices_, int dimension_, int multiplicity_, bool sortWithSymmetry, SymmetricComplex const &complex):
  dimension(dimension_),
  multiplicity(multiplicity_),
  isKnownToBeNonMaximalFlag(false)
{
  indices=vector<int>(indices_.size());
  int j=0;
  for(set<int>::const_iterator i=indices_.begin();i!=indices_.end();i++,j++)
    indices[j]=*i;

  IntegerMatrix const &vertices=complex.getVertices();
  IntegerVector sum(vertices.getWidth());
  for(int i=0;i<indices.size();i++)
    sum+=vertices[indices[i]];


  if(sortWithSymmetry)
    {
      sortKey=complex.sym.orbitRepresentative(sum,&sortKeyPermutation);
//      sortKey=sum;

      /*
	LexicographicTermOrder myOrder;

      for(SymmetryGroup::ElementContainer::const_iterator k=complex.sym.elements.begin();k!=complex.sym.elements.end();k++)
	if(myOrder(SymmetryGroup::compose(*k,sum),sortKey))
	  sortKey=SymmetryGroup::compose(*k,sum);
      */


      /*
      int n=sum.size();
      for(SymmetryGroup::ElementContainer::const_iterator k=complex.sym.elements.begin();k!=complex.sym.elements.end();k++)
	{
	  bool isBetter=true;
	  for(int i=0;i<n;i++)
	    {
	      if(sum[(*k)[i]]>sortKey[i]){isBetter=false;break;}
	      if(sum[(*k)[i]]<sortKey[i])break;
	    }
	  if(isBetter)
	    {
	      sortKey=SymmetryGroup::compose(*k,sum);
	    }
	}
*/

    }
  else
    {
      sortKey=sum;
    }
}


int SymmetricComplex::indexOfVertex(IntegerVector const &v)const
{
  map<IntegerVector,int>::const_iterator it=indexMap.find(v);
  assert(it!=indexMap.end());
  return it->second;
}


void SymmetricComplex::Cone::remap(SymmetricComplex &complex)
{
  IntegerMatrix const &vertices=complex.getVertices();
  IntegerVector sum(vertices.getWidth());
  for(int i=0;i<indices.size();i++)
    sum+=vertices[indices[i]];

  int n=sum.size();
/*  IntegerVector bestPermutation;
  for(SymmetryGroup::ElementContainer::const_iterator k=complex.sym.elements.begin();k!=complex.sym.elements.end();k++)
    {
      if(SymmetryGroup::compose(*k,sum)==sortKey)
	bestPermutation=*k;
    }
    */
  IntegerVector const &bestPermutation=sortKeyPermutation;

  assert(bestPermutation.size()==n);

  vector<int> indicesNew(indices.size());
  int I=0;
  for(vector<int>::const_iterator i=indices.begin();i!=indices.end();i++,I++)
    {
      IntegerVector ny=SymmetryGroup::compose(bestPermutation,complex.vertices[*i]);
      map<IntegerVector,int>::const_iterator it=complex.indexMap.find(ny);
      assert(it!=complex.indexMap.end());
      indicesNew[I]=it->second;
    }
  indices=indicesNew;
}


set<int> SymmetricComplex::Cone::indexSet()const
{
  set<int> ret;
  for(vector<int>::const_iterator i=indices.begin();i!=indices.end();i++)
    ret.insert(*i);

  return ret;
}

bool SymmetricComplex::Cone::isSubsetOf(Cone const &c)const
{
  int next=0;
  for(int i=0;i<indices.size();i++)
    {
      while(1)
	{
	  if(next>=c.indices.size())return false;
	  if(indices[i]==c.indices[next])break;
	  next++;
	}
    }
  return true;
  /*
  set<int> b=c.indexSet();

  for(vector<int>::const_iterator i=indices.begin();i!=indices.end();i++)
    if(b.count(*i)==0)return false;
  return true;
  */
}


SymmetricComplex::Cone SymmetricComplex::Cone::permuted(IntegerVector const &permutation, SymmetricComplex const &complex, bool withSymmetry)const
{
  /*  Cone ret;
  ret.dimension=dimension;
  ret.multiplicity=multiplicity;*/

  set<int> r;
  for(vector<int>::const_iterator i=indices.begin();i!=indices.end();i++)
    {
      IntegerVector ny=SymmetryGroup::compose(permutation,complex.vertices[*i]);
      map<IntegerVector,int>::const_iterator it=complex.indexMap.find(ny);
      if(it==complex.indexMap.end())
	{
	  AsciiPrinter(Stderr).printVector(complex.vertices[*i]);
	  AsciiPrinter(Stderr).printVector(ny);

	  assert(0);
	}
      r.insert(it->second);
    }


  return Cone(r,dimension,multiplicity,withSymmetry,complex);
}


/*
void SymmetricComplex::Cone::computeRelativeInteriorPoint(SymmetricComplex const &complex)
{
  IntegerMatrix const &vertices=complex.getVertices();
  IntegerVector sum(vertices.getWidth());
  for(const_iterator i=begin();i!=end();i++)
    sum+=vertices[*i];
  relativeInteriorPoint=sum;
  sum.sort();
  summary=sum;
}
*/

 /*void SymmetricComplex::Cone::computeSmallestRepresentative(SymmetricComplex const &complex)
{
  if(relativeInteriorPoint.size()==0)computeRelativeInteriorPoint(complex);

  LexicographicTermOrder myOrder;

  smallestRepresentative=relativeInteriorPoint;

  for(SymmetryGroup::ElementContainer::const_iterator k=complex.sym.elements.begin();k!=complex.sym.elements.end();k++)
    if(myOrder(SymmetryGroup::compose(*k,relativeInteriorPoint),smallestRepresentative))
      smallestRepresentative=SymmetryGroup::compose(*k,relativeInteriorPoint);
}
 */

bool SymmetricComplex::Cone::operator<(Cone const & b)const
{
  /*  if(ignoreSymmetry)return ((set<int>)*this)<(b);*/
  return sortKey<b.sortKey;
}


bool SymmetricComplex::Cone::isSimplicial(int linealityDim)const
{
  return (indices.size()+linealityDim)==dimension;
}


IntegerVectorList SymmetricComplex::Cone::orthogonalComplement(SymmetricComplex &complex)const
{
	IntegerVectorList l;
	for(int i=0;i<indices.size();i++)
		l.push_back(complex.vertices[indices[i]]);

	FieldMatrix m=integerMatrixToFieldMatrix(rowsToIntegerMatrix(l,complex.n),Q);
	return fieldMatrixToIntegerMatrixPrimitive(m.reduceAndComputeKernel()).getRows();
}


SymmetricComplex::SymmetricComplex(int n_, IntegerVectorList const &v, SymmetryGroup const &sym_):
  n(n_),
  sym(sym_),
  dimension(-1)
{
  vertices=rowsToIntegerMatrix(v,n);
  for(int i=0;i<vertices.getHeight();i++)indexMap[vertices[i]]=i;
}


bool SymmetricComplex::contains(Cone const &c)const
{
  Cone temp=c;//#1
  /*  temp.computeRelativeInteriorPoint(*this);

  temp.computeSmallestRepresentative(*this);
  */
 return cones.find(temp)!=cones.end();///////////////////!!!!!!!!!!!!!!!!!!!!!!!


 /*
  set<IntegerVector> possibleMatches;
  for(ConeContainer::const_iterator i=cones.begin();i!=cones.end();i++)
    {
 */
   /*      AsciiPrinter(Stderr).printVector(i->summary);
      AsciiPrinter(Stderr).printVector(temp.summary);
      fprintf(stderr,"\n");*/
 /*      if(i->dimension==temp.dimension)
      if(i->summary==temp.summary)
	{
	  possibleMatches.insert(i->relativeInteriorPoint);
	}
    }
  for(SymmetryGroup::ElementContainer::const_iterator k=sym.elements.begin();k!=sym.elements.end();k++)
    {
      if(possibleMatches.find(SymmetryGroup::compose(*k,temp.relativeInteriorPoint))!=possibleMatches.end())
	return true;
    }
 */
  /*  for(SymmetryGroup::ElementContainer::const_iterator k=sym.elements.begin();k!=sym.elements.end();k++)
    {
      Cone c2=c.permuted(*k,*this);
      for(list<Cone>::const_iterator i=cones.begin();i!=cones.end();i++)
	{
	  if(c2==*i)return true;
	}
	}*/

 //  return false;
}


void SymmetricComplex::insert(Cone const &c)
{
	if(c.dimension>dimension)dimension=c.dimension;
  //  Cone temp=c;
  /*  temp.computeRelativeInteriorPoint(*this);
  temp.computeSmallestRepresentative(*this);
  */
  if(!contains(c))//#2
    {
      cones.insert(c);
      //cones.push_back(temp);
      //      fprintf(Stderr,"INSERTING\n");
      //      cones.back().computeRelativeInteriorPoint(*this);
    }
  else
    {
      //      if(c.isKnownToBeNonMaximal())cones.find(c)->setKnownToBeNonMaximal();
      if(c.isKnownToBeNonMaximal()){cones.erase(c);cones.insert(c);}// mark as non-maximal
    }
}


int SymmetricComplex::getMaxDim()const
{/*
  int ret=-1;
  for(ConeContainer::const_iterator i=cones.begin();i!=cones.end();i++)
    {
      if(i->dimension>ret)ret=i->dimension;
    }
  return ret;*/
	return dimension;
}


int SymmetricComplex::getMinDim()const
{
  int ret=100000;
  for(ConeContainer::const_iterator i=cones.begin();i!=cones.end();i++)
    {
      if(i->dimension<ret)ret=i->dimension;
    }
  return ret;
}


bool SymmetricComplex::isMaximal(Cone const &c)const
{
  if(c.isKnownToBeNonMaximal())return false;
  if(c.dimension==dimension)return true;
  for(SymmetryGroup::ElementContainer::const_iterator k=sym.elements.begin();k!=sym.elements.end();k++)
    {
      Cone c2=c.permuted(*k,*this,false);
      for(ConeContainer::const_iterator i=cones.begin();i!=cones.end();i++)
	{
	  if(i->dimension>c.dimension)
	    if(c2.isSubsetOf(*i) && !i->isSubsetOf(c2))return false;
	}
    }
  return true;
}


IntegerVector SymmetricComplex::dimensionsAtInfinity()const
{
  /* Using a double description like method this routine computes the
     dimension of the intersection of each cone in the complex with
     the plane x_0=0 */

  IntegerVector ret(cones.size());

  int I=0;
  for(ConeContainer::const_iterator i=cones.begin();i!=cones.end();i++,I++)
    {
      IntegerVectorList raysAtInfinity;
      for(vector<int>::const_iterator j=i->indices.begin();j!=i->indices.end();j++)
	{
	  if(vertices[*j][0]==0)raysAtInfinity.push_back(vertices[*j]);
	  for(vector<int>::const_iterator k=j;k!=i->indices.end();k++)
	    if(vertices[*j][0]*vertices[*k][0]<0)
	      raysAtInfinity.push_back(((vertices[*j][0]>0)?1:-1)*(vertices[*j][0])*vertices[*k]+
				       ((vertices[*k][0]>0)?1:-1)*(vertices[*k][0])*vertices[*j]);
	}
      ret[I]=rankOfMatrix(raysAtInfinity);
    }
  return ret;
}


string SymmetricComplex::toString(int dimLow, int dimHigh, bool onlyMaximal, bool group, ostream *multiplicities, bool compressed, bool tPlaneSort, bool xml)const
{
  stringstream ret;

  if(!onlyMaximal)
    {
      if(xml)ret<<"<m>\n";
      if(xml)if(multiplicities)*multiplicities<<"<m>\n";
    }
  IntegerVector additionalSortKeys(cones.size());
  if(tPlaneSort)additionalSortKeys=dimensionsAtInfinity();
  int lowKey=additionalSortKeys.min();
  int highKey=additionalSortKeys.max();

  if(xml)if(onlyMaximal)
    {
      if(compressed)
        ret<<"<m>\n";
      else
        ret<<"<m cols=\""<<this->vertices.getHeight()<<"\">\n";
      if(multiplicities)*multiplicities<<"<m>\n";
    }
  for(int d=dimLow;d<=dimHigh;d++)
    {
      log1 cerr << "Processing dimension "<<d<<".\n";
      int numberOfOrbitsOutput=0;
      int numberOfOrbitsOfThisDimension=0;
      bool newDimension=true;
      for(int key=lowKey;key<=highKey;key++)
	{
	  if(xml)
	    {
	      if(!onlyMaximal)
	        {
	          if(compressed)
	            ret<<"<m>\n";
	          else
	            ret<<"<m cols=\""<<this->vertices.getHeight()<<"\">\n";
	          if(multiplicities)*multiplicities<<"<m>\n";
	        }
	    }
	  int I=0;
	  for(ConeContainer::const_iterator i=cones.begin();i!=cones.end();i++,I++)
	    if(additionalSortKeys[I]==key)
		  if(i->dimension==d)
		    {
		      numberOfOrbitsOfThisDimension++;
	      if(!onlyMaximal || isMaximal(*i))
		{
		  numberOfOrbitsOutput++;
		  bool isMax=isMaximal(*i);
		  bool newOrbit=true;
		  set<set<int> > temp;
		    for(SymmetryGroup::ElementContainer::const_iterator k=sym.elements.begin();k!=sym.elements.end();k++)
		      {
			Cone temp1=i->permuted(*k,*this,false);
			temp.insert(temp1.indexSet());
			if(compressed)break;
		    }
		  for(set<set<int> >::const_iterator j=temp.begin();j!=temp.end();j++)
		    {
		      if(!xml)ret << "{";else ret<<"<v>";
		      for(set<int>::const_iterator a=j->begin();a!=j->end();a++)
			{
			  if(a!=j->begin())ret<<" ";
			  ret << *a;
			}
                      if(!xml)ret << "}";else ret<<"</v>";
                      if(!xml)
                        {
                           if(group)if(newOrbit)ret << "\t# New orbit";
                           if(newDimension)ret << "\t# Dimension "<<d;
                        }
		      ret <<endl;
		      if(isMax)if(multiplicities)
			{
			  if(xml)*multiplicities<<"<v>";
			  *multiplicities << i->multiplicity;
                          if(xml)*multiplicities<<"</v>";
			  if(!xml)
			    {
			      if(group)if(newOrbit)*multiplicities << "\t# New orbit";
			      if(newDimension)*multiplicities << "\t# Dimension "<<d;
			    }
			  *multiplicities << endl;
			}
		      newOrbit=false;
		      newDimension=false;
		    }
	      }
		    }
	}
      if(xml)
        {
          if(!onlyMaximal)
            {
              ret<<"</m>\n";
              if(multiplicities)*multiplicities<<"</m>\n";
            }
        }

      log1 cerr<<"Number of orbits of this dimension: " << numberOfOrbitsOfThisDimension << endl;
      log1 cerr<<"Number of orbits output: " << numberOfOrbitsOutput << endl;
    }
  if(xml)if(onlyMaximal)
    {
      ret<<"</m>\n";
      if(multiplicities)*multiplicities<<"</m>\n";
    }

  if(!onlyMaximal)
    {
      if(xml)ret<<"</m>\n";
      if(xml)if(multiplicities)*multiplicities<<"</m>\n";
    }

  return ret.str();
}


IntegerVector SymmetricComplex::fvector(bool boundedPart)const
{
  int min=getMinDim();
  IntegerVector ret(getMaxDim()-min+1);

  for(ConeContainer::const_iterator i=cones.begin();i!=cones.end();i++)
    {
      /*      set<Cone> temp;
      for(SymmetryGroup::ElementContainer::const_iterator k=sym.elements.begin();k!=sym.elements.end();k++)
      temp.insert(i->permuted(*k,*this));*/
      /*      set<IntegerVector> temp;
      for(SymmetryGroup::ElementContainer::const_iterator k=sym.elements.begin();k!=sym.elements.end();k++)
	temp.insert(SymmetryGroup::compose(*k,i->sortKey));
      */
      bool doAdd=!boundedPart;
      if(boundedPart)
	{
	  bool isBounded=true;
	  for(vector<int>::const_iterator j=i->indices.begin();j!=i->indices.end();j++)
	    if(vertices[*j][0]==0)isBounded=false;
	  doAdd=isBounded;
	}
      if(doAdd)
	ret[i->dimension-min]+=sym.orbitSize(i->sortKey);
    }
  return ret;
}


bool SymmetricComplex::isPure()const
{
  int dim=-1;
  for(ConeContainer::const_iterator i=cones.begin();i!=cones.end();i++)
    {
      log2{static int a;if(!((a++)&63))fprintf(Stderr,"%i\n",a);}//log0
    if(isMaximal(*i))
      {
	int dim2=i->dimension;
	if(dim==-1)dim=dim2;
	if(dim!=dim2)return false;
      }
    }
  return true;
}


bool SymmetricComplex::isSimplicial()const
{
  int linealityDim=getMinDim();
  for(ConeContainer::const_iterator i=cones.begin();i!=cones.end();i++)
    if(!i->isSimplicial(linealityDim))
      return false;
  return true;
}


void SymmetricComplex::remap()
{
  for(ConeContainer::iterator i=cones.begin();i!=cones.end();i++)
    {
      Cone const&j=*i;
      Cone &j2=const_cast<Cone&>(j);//DANGER: cast away const. This does not change the sort key in the container, so should be OK.
      j2.remap(*this);
    }
}


int SymmetricComplex::numberOfConesOfDimension(int d)const
{
	assert(sym.isTrivial());

	int ret=0;
	for(ConeContainer::const_iterator i=cones.begin();i!=cones.end();i++)
		if(d==i->dimension)
		{
				ret++;
		}
	return ret;
}


int SymmetricComplex::dimensionIndex(Cone const &c)
{
	assert(sym.isTrivial());

	int ret=0;
	for(ConeContainer::const_iterator i=cones.begin();i!=cones.end();i++)
		if(c.dimension==i->dimension)
		{
			if(!(c<*i)&&!(*i<c))
				return ret;
			else
				ret++;
		}
	return ret;
}


void SymmetricComplex::boundary(Cone const &c, vector<int> &indices_, vector<int> &signs)
{
	indices_=vector<int>();
	signs=vector<int>();
	int d=c.dimension;


	IntegerVectorList l;
	for(int i=0;i<c.indices.size();i++)
		l.push_back(vertices[c.indices[i]]);
	IntegerVectorList facetNormals=PolyhedralCone(l,IntegerVectorList(),n).extremeRays();
	IntegerVectorList complementBasis=c.orthogonalComplement(*this);
	for(IntegerVectorList::const_iterator i=facetNormals.begin();i!=facetNormals.end();i++)
	{
		IntegerVectorList complementBasis1=complementBasis;
		complementBasis1.push_back(*i);
		FieldMatrix m=integerMatrixToFieldMatrix(rowsToIntegerMatrix(complementBasis1,n),Q);
		IntegerVectorList completion=fieldMatrixToIntegerMatrixPrimitive(m.reduceAndComputeKernel()).getRows();
		for(IntegerVectorList::const_iterator j=completion.begin();j!=completion.end();j++)complementBasis1.push_back(*j);
		int sign=determinantSign(complementBasis1);



		set<int> indices;
		for(vector<int>::const_iterator j=c.indices.begin();j!=c.indices.end();j++)if(dotLong(vertices[*j],*i)==0)indices.insert(*j);
		Cone facet(indices,d-1,1,true,*this);
		IntegerVectorList complementBasis2=facet.orthogonalComplement(*this);
		for(IntegerVectorList::const_iterator j=completion.begin();j!=completion.end();j++)complementBasis2.push_back(*j);
		indices_.push_back(dimensionIndex(facet));
		signs.push_back(sign*determinantSign(complementBasis2));
	}
}


IntegerMatrix SymmetricComplex::boundaryMap(int d)
{
	assert(sym.isTrivial());

	IntegerMatrix ret(numberOfConesOfDimension(d-1),numberOfConesOfDimension(d));

	for(ConeContainer::const_iterator i=cones.begin();i!=cones.end();i++)
		if(d==i->dimension)
		{
			int I=dimensionIndex(*i);
			vector<int> indices;
			vector<int> signs;
			boundary(*i,indices,signs);
			for(int j=0;j<indices.size();j++)
			{
				ret[indices[j]][I]+=signs[j];
			}
		}
	return ret;
}

