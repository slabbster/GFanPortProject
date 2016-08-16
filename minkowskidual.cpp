#include "minkowskidual.h"
#include "field_rationals.h"
#include "linalg.h"



class MyCone{
public:
  struct Pair{
    int rayIndex;
    vector<IntegerVector> const &rays;
    IntegerVector const &permutation;
    Pair(vector<IntegerVector> const &rays_, int rayIndex_, IntegerVector const & permutation_):
      rays(rays_),
      rayIndex(rayIndex_),
      permutation(permutation_)
    {
    }
    bool operator<(Pair const &b)const
    {
      return SymmetryGroup::compose(permutation,rays[rayIndex])<SymmetryGroup::compose(b.permutation,b.rays[b.rayIndex]);
    }
  };
  set<Pair> rays;    
  void print()
  {
    for(set<Pair>::const_iterator i=rays.begin();i!=rays.end();i++)
      {
	fprintf(Stderr,"--:%i ",i->rayIndex);
	AsciiPrinter(Stderr).printVector((i->permutation));
	fprintf(Stderr,"\n");
      }
  }
};

 static list<SymmetricComplex::Cone> computeFacets(SymmetricComplex::Cone const &theCone, IntegerMatrix const &rays, IntegerVectorList const &facetCandidates, SymmetricComplex const &theComplex/*, int linealityDim*/)
{
  set<SymmetricComplex::Cone> ret;

  for(IntegerVectorList::const_iterator i=facetCandidates.begin();i!=facetCandidates.end();i++)
    {
      set<int> indices;

      bool notAll=false;
      for(vector<int>::const_iterator j=theCone.indices.begin();j!=theCone.indices.end();j++)
	if(dotLong(rays[*j],*i)==0)
	  indices.insert(*j);
	else
	  notAll=true;

      SymmetricComplex::Cone temp(indices,theCone.dimension-1,0,false,theComplex);
      /*      temp.multiplicity=0;
      temp.dimension=theCone.dimension-1;
      temp.setIgnoreSymmetry(true);
      */
      if(notAll)ret.insert(temp);

    }
  //  fprintf(Stderr,"HEJ!!!!\n");

  list<SymmetricComplex::Cone> ret2;
  for(set<SymmetricComplex::Cone>::const_iterator i=ret.begin();i!=ret.end();i++)
    {
      bool isMaximal=true;

      /*      if(i->indices.size()+linealityDim<i->dimension)//#3
	isMaximal=false;
	else*/
	for(set<SymmetricComplex::Cone>::const_iterator j=ret.begin();j!=ret.end();j++)
	  if(i!=j && i->isSubsetOf(*j))
	    {
	      isMaximal=false;
	      break;
	    }
      if(isMaximal)
	{
	  SymmetricComplex::Cone temp(i->indexSet(),i->dimension,i->multiplicity,true,theComplex);
	  //	  temp.setIgnoreSymmetry(false);
	  ret2.push_back(temp);
	}
    }
  return ret2;
}


SymmetricComplex dualMinkowskiMixed(PolynomialSet const &g, SymmetryGroup const &sym, PolyhedralFan const &cones)
{
  /* The set of rays of the dual fan is the vertices of the Minkowski
     sum after centering around zero. We only generate them up to symmetry. */
  int n=cones.getAmbientDimension();
  int nMaxOrbits=0;
  for(PolyhedralFan::coneIterator i=cones.conesBegin();i!=cones.conesEnd();i++)nMaxOrbits++;
  vector<IntegerVector> rays(nMaxOrbits);
  FieldVector translationVector(Q,n);
  for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
    translationVector+=Q.zHomomorphism(i->numberOfTerms()).inverse()*integerVectorToFieldVector(i->exponentsSum(),Q);
  {
    int j=0;
    for(PolyhedralFan::coneIterator i=cones.conesBegin();i!=cones.conesEnd();i++,j++)
      {
	FieldVector temp=integerVectorToFieldVector(i->getRelativeInteriorPoint(),Q);
	temp+=(Q.zHomomorphism(-1)*translationVector);
	rays[j]=temp.primitive();
      }
  }

  fprintf(Stderr,"The new rays:\n");
  for(int i=0;i<nMaxOrbits;i++)
    {
      fprintf(Stderr,"%i:",i);
      AsciiPrinter(Stderr).printVector(rays[i]);
      fprintf(Stderr,"\n");      
    }

  /* The new cones of the dual fan comes from the rays of the origanal
     Minkowski sum's normal fan. For every primary ray we add the set of dual rays whose primary cone contains the primary ray.

     That is, add the vertices to the facets of the minkowskisum.
   */
  IntegerVectorList primaryRays=cones.getRaysInPrintingOrder(&sym,true);

  vector<MyCone> newCones(primaryRays.size());
  IntegerVectorList::const_iterator K=primaryRays.begin();
  for(int k=0;k<newCones.size();k++,K++)
    {
      for(SymmetryGroup::ElementContainer::const_iterator l=sym.elements.begin();l!=sym.elements.end();l++)
	{
	  IntegerVector v=SymmetryGroup::compose(*l,*K);
	  int j=0;
	  for(PolyhedralFan::coneIterator i=cones.conesBegin();i!=cones.conesEnd();i++,j++)
	    {
	      if(i->contains(v))
		{
		  newCones[k].rays.insert(MyCone::Pair(rays,j,*l));
		}
	    }
	}
    }

  fprintf(Stderr,"Facets with their vertices:\n");
  for(int k=0;k<newCones.size();k++)
    {
      fprintf(Stderr,"Printing facet:%i\n\n",k);
      //      newCones[k].print();
    }
  


  /* To make things easier we store for each primary cone the indices of rays it contains (among all rays)*/
  IntegerVectorList allPrimaryRays=cones.getRaysInPrintingOrder(&sym,false);
  vector<set<int> > allRaysIndices(nMaxOrbits);

  {
    int j=0;
    for(PolyhedralFan::coneIterator i=cones.conesBegin();i!=cones.conesEnd();i++,j++)
      {
	int k=0;
	for(IntegerVectorList::const_iterator l=allPrimaryRays.begin();l!=allPrimaryRays.end();l++,k++)
	  if(i->contains(*l))allRaysIndices[j].insert(k);
      }
  }

  fprintf(Stderr,"Facets with their vertices (as indices):\n");
  for(int k=0;k<nMaxOrbits;k++)
    {
      
      fprintf(Stderr,"Cone %i\n",k);
      for(set<int>::const_iterator j=allRaysIndices[k].begin();j!=allRaysIndices[k].end();j++)
	fprintf(Stderr,"%i ",*j);
      fprintf(Stderr,"\n\n");
      //      newCones[k].print();
    }

  
  /* We prepare the complex to return. */

  SymmetricComplex ret(n,allPrimaryRays,sym);

  /* For every dual cone we need to compute facet normal. This is done
     by cddlib. Afterward we use these facet normals to add the faces
     of the dual cone to our collection of cones (pre-complex) - each
     time checking that the corresponding face is mixed. */
  
  fprintf(Stderr,"Adding mixed faces for cones:\n");
  for(int k=0;k<newCones.size();k++)
    {
      IntegerVectorList a;
      for(set<MyCone::Pair>::const_iterator i=newCones[k].rays.begin();i!=newCones[k].rays.end();i++)
	{
	  a.push_back(SymmetryGroup::compose(i->permutation,rays[i->rayIndex]));	  
	}

      AsciiPrinter(Stderr).printVectorList(a);

      IntegerVectorList empty;
      PolyhedralCone C(a,empty,n);
      IntegerVectorList facetCandidates=C.extremeRays();


      fprintf(Stderr,"Cone %i\nFacetNormals:\n",k);
      AsciiPrinter(Stderr).printVectorList(facetCandidates);
      /*      list<SymmetricComplex::Cone> clist;
      {
	set<int> indices;
	for(int j=0;j<rays.getHeight();j++)if(cone.contains(rays[j]))indices.insert(j);
	SymmetricComplex::Cone temp(indices,cone.dimension(),cone.getMultiplicity(),true,c);
	clist.push_back(temp);
      }

      while(!clist.empty())
	{
	  SymmetricComplex::Cone &theCone=clist.front();
	  
	  if()//Check mixed face property
	  if(!ret.contains(theCone))
	    {
	      c.insert(theCone);
	      
	      list<SymmetricComplex::Cone> facets=computeFacets(theCone,rays,facetCandidates,c);
	      clist.splice(clist.end(),facets);
	    }
	  clist.pop_front();
	}
      */      
    }



  /* For every cone in the final "pre-complex" we translate it into a cone in the primary */

  return ret;
}
