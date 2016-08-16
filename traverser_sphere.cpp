#include "traverser_sphere.h"
#include "determinant.h"

#include "printer.h"

//void updatePolyhedralCone();
SphereTraverser::SphereTraverser(vector<PolyhedralCone> const &cones_, map<IntegerVector,list<int> > const &adjacency_, IntegerVector const &startCone, IntegerVector const &normal):
  adjacency(adjacency_),
  cones(cones_)
{
  currentConeIndex=0;
  for(;currentConeIndex<cones.size();currentConeIndex++)
    {
      if(cones[currentConeIndex].contains(startCone))
        break;
    }
  assert(currentConeIndex!=cones.size());

  currentNormal=normal;
}

void SphereTraverser::changeCone(IntegerVector const &ridgeVector, IntegerVector const &rayVector)
{
  list<int> const &candidates=adjacency[ridgeVector];
  for(list<int>::const_iterator i=candidates.begin();i!=candidates.end();i++)
    {
        PolyhedralCone l=cones[*i].link(ridgeVector);
        if(l.contains(rayVector))
          {
            IntegerVectorList M;
            M.push_back(cones[currentConeIndex].getRelativeInteriorPoint());
            M.push_back(currentNormal);
            currentConeIndex=*i;
            IntegerVectorList eq=cones[*i].getEquations();
            assert(eq.size()==1);

            MatrixTermOrder T(M);
            if(!T(*eq.begin(),currentNormal-currentNormal))//HERE
            //            if(dotLong(currentNormal,*eq.begin())>0)
              {currentNormal=*eq.begin();}
            else
              {currentNormal=-*eq.begin();}

            break;
          }
//        assert(0);
    }
}

IntegerVectorList SphereTraverser::link(IntegerVector const &ridgeVector)
{
  list<int> const &candidates=adjacency[ridgeVector];

  PolyhedralCone ridge=cones[currentConeIndex].faceContaining(ridgeVector);
  ridge.canonicalize();
//  IntegerVectorList partialBasis=ridge.getEquations();
  IntegerVectorList partialBasis=ridge.link(ridge.getRelativeInteriorPoint()).dualCone().getEquations();

  IntegerVector myRay=cones[currentConeIndex].link(ridgeVector).getUniquePoint();

  partialBasis.push_back(currentNormal);
  partialBasis.push_back(myRay);
  int sign=determinantSign(partialBasis);
  partialBasis.pop_back();
  partialBasis.pop_back();

  IntegerVectorList c2;

  for(list<int>::const_iterator i=candidates.begin();i!=candidates.end();i++)
    {
      debug << "TESTET\n";
      PolyhedralCone temp=cones[*i].link(ridgeVector);
      IntegerVector v=temp.getUniquePoint();
      if(dotLong(v,currentNormal)>=0)c2.push_back(v);
    }
  assert(c2.size()>0);
  IntegerVectorList ret2;
  ret2.push_back(myRay);
  IntegerVector ret;
  for(IntegerVectorList::const_iterator i=c2.begin();i!=c2.end();i++)
    if(*i!=myRay){
      ret=*i;
      break;
    }

  if(ret.size()==0)
    {
      PolyhedralFan f(theCone.ambientDimension());
      f.insert(theCone);
      f.printWithIndices(&pout);
      pout<<currentNormal;
    }
  for(IntegerVectorList::const_iterator i=c2.begin();i!=c2.end();i++)
    {
      if(*i==myRay) continue;
      partialBasis.push_back(ret);
      partialBasis.push_back(*i);
      int s=determinantSign(partialBasis);
      partialBasis.pop_back();
      partialBasis.pop_back();
      if(sign==s)
        {
          ret=*i;
        }
    }
  ret2.push_back(ret);

  debug<<ret2;

  return ret2;
}


PolyhedralCone  & SphereTraverser::refToPolyhedralCone()
{
  theCone=cones[currentConeIndex];
  return theCone;//cones[currentCone];
}
