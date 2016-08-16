#include "fieldlp.h"

FieldVector FieldLP::subvector(FieldVector const &v, set<int> const &b)
{
  FieldVector ret(v.getField(),b.size());
  int I=0;
  for(set<int>::const_iterator i=b.begin();i!=b.end();i++,I++)
    ret[I]=v[*i];

  return ret;
}


FieldMatrix FieldLP::submatrix(FieldMatrix const &m, set<int> const &b)
{
  FieldMatrix ret(m.getField(),b.size(),m.getWidth());
  int I=0;
  for(set<int>::const_iterator i=b.begin();i!=b.end();i++,I++)
    ret[I]=m[*i];

  return ret;
}

FieldVector FieldLP::edgeDirection(int i)const
{
  set<int> basis2=basis;
  assert(basis2.count(i));
  basis2.erase(i);
  FieldMatrix A2=submatrix(A,basis2);
  FieldMatrix ker=A2.reduceAndComputeKernel();
  assert(ker.getHeight()==1);
  FieldVector e=ker[0];
  if(dot(A[i],e).sign()<0)return e;
  return theField.zHomomorphism(-1)*e;
}

bool FieldLP::isImprovingDirection(int i)const
{
  return dot(edgeDirection(i),w).sign()>0;
}

FieldElement FieldLP::improvement(int i, int &newBasisMember)const
{
  FieldVector e=edgeDirection(i);
  FieldElement ew=dot(e,w);
  FieldElement ret(theField);
  bool first=true;
  newBasisMember=-1;
  for(int s=0;s<A.getHeight();s++)//This is _the_ line where an anticycling rule is implemented
  //  for(int s=A.getHeight()-1;s>=0;s--)//This is _the_ line where an anticycling rule is implemented
    {
      FieldElement denominator=dot(A[s],e);
      if(denominator.sign()>0)
	{
	  FieldElement imp=(b[s]-dot(A[s],x))*denominator.inverse()*ew;
	  if(first){ret=imp;newBasisMember=s;first=false;}
	  if((imp-ret).sign()<=0){newBasisMember=s;ret=imp;}
	}
    }
  return ret;
}

int FieldLP::step()
{
  FieldElement impBest(theField);
  int bestIndex=-1;
  int bestNewElement=-1;
  for(set<int>::const_iterator i=basis.begin();i!=basis.end();i++)
    {
      //  AsciiPrinter P(Stderr);
      //  edgeDirection(*i).print(P);
    if(isImprovingDirection(*i))
      {
	//fprintf(Stderr,"IMPROVING DIRECTION\n");
	int newBasisElement;
	FieldElement imp=improvement(*i,newBasisElement);
	//P.printFieldElement(imp);
	//P.printNewLine();
	if(newBasisElement==-1)return -1; //UNBOUNDED
	if((imp-impBest).sign()>=0)//needed for perturbation
	  {
	    impBest=imp;
	    bestIndex=*i;
	    bestNewElement=newBasisElement;
	  }
      }
    }
  if(bestIndex==-1)
    return 0; //OPTIMAL

  set<int> basis2=basis;
  //fprintf(Stderr,"REMOVING %i INSERTING %i\n",bestIndex,bestNewElement);
  basis2.erase(bestIndex);
  basis2.insert(bestNewElement);
  setCurrentBasis(basis2);
  return 1;
}

FieldLP::FieldLP(FieldMatrix const &A_, FieldVector const &b_):
  A(A_),
  b(b_),
  theField(A_.getField()),
  x(A_.getField(),0),
  w(A.getField(),A_.getWidth())
{
  //  if(A.getHeight()!=b.size())fprintf(Stderr,"%i %i",A.getHeight(),b.size());
  assert(A.getHeight()==b.size());
}


void FieldLP::setObjectiveFunction(FieldVector const &w_)
{
  w=w_;
  assert(A.getWidth()==w.size());
}


FieldVector FieldLP::basisToPoint(set<int> const &basis)const
{
  AsciiPrinter P(Stdout);
  //  FieldMatrix AT=A.transposed();
  //  AT.printMatrix(P);
  FieldMatrix A2(A.getField(),basis.size(),A.getWidth());
  FieldVector b2(A.getField(),basis.size());
  int I=0;
  for(set<int>::const_iterator i=basis.begin();i!=basis.end();i++,I++)
    {
      A2[I]=A[*i];
      b2[I]=b[*i];
    }
  //A2.printMatrix(P);
  //b2.print(P);
  FieldMatrix s=A2.solver();
  //fprintf(Stderr,"%i %i\n",s.getWidth(),b2.size()+A.getHeight());
  //s.printMatrix(P);
  //concatenation(b2,FieldVector(b2.getField(),A.getWidth())).print(P);
  FieldVector r=s.canonicalize(concatenation(b2,FieldVector(b2.getField(),A.getWidth())));
  //r.print(P);
  
  assert(r.subvector(0,basis.size()).isZero());

  return r.subvector(basis.size(),r.size());
}


void FieldLP::setCurrentBasis(set<int> const &basis_)
{
  basis=basis_;
  x=basisToPoint(basis);
}


void FieldLP::print(Printer &P)const
{
  P.printString("LPprogram\n");
  P.printString("A:\n");
  A.printMatrix(P);
  P.printString("w:\n");
  w.print(P);
  P.printString("b:\n");
  b.print(P);
  P.printString("current basis:\n");
  for(set<int>::const_iterator i=basis.begin();i!=basis.end();i++)P.printInteger(*i);  
  P.printString("x:\n");
  x.print(P);  
}

FieldMatrix FieldLP::computeLinealitySpace()
{
  FieldMatrix temp=A;
  return temp.reduceAndComputeKernel();
}


FieldLP FieldLP::buildLPForFeasibilityCheck()
{
  assert(computeLinealitySpace().getHeight()==0);
  set<int> S;
  {
    FieldMatrix temp=A.flipped().transposed();//important that inequalities which do not appear in basis set are
    temp.reduce();
    int i=-1;
    int j=-1;
    while(temp.nextPivot(i,j))
      {
	S.insert(A.getHeight()-1-j);//perturbed furthest
	//fprintf(Stderr,"INSERTING %i\n",A.getHeight()-1-j);
      }
  }
  FieldMatrix AS=submatrix(A,S);
  FieldMatrix ASSolver=AS.solver();
  FieldVector xs2=ASSolver.canonicalize(concatenation(subvector(b,S),FieldVector(b.getField(),A.getWidth()))).subvector(S.size(),S.size()+A.getWidth());



  set<int> inequalitiesWithSlack;
  for(int i=0;i<A.getHeight();i++)
    {
      if(!S.count(i))
	if((dot(A[i],xs2)-b[i]).sign()>0)
	  {
	    inequalitiesWithSlack.insert(i);
	    //fprintf(Stderr,"adding slack variable\n");
	  }
    }
  FieldMatrix newA=combineOnTop(
				combineLeftRight(
						 A,
						 FieldMatrix(
							     theField,
							     A.getHeight(),
							     inequalitiesWithSlack.size()
							     )
						 ),
				FieldMatrix(
					    theField,
					    inequalitiesWithSlack.size(),
					    A.getWidth()+inequalitiesWithSlack.size()
					    )
				);
  FieldVector newW(theField,inequalitiesWithSlack.size()+A.getWidth());
  int I=0;
  for(set<int>::const_iterator i=inequalitiesWithSlack.begin();i!=inequalitiesWithSlack.end();i++,I++)
    {
      newA[*i][I+A.getWidth()]=theField.zHomomorphism(-1);
      //      int sign=dot(A[*i],xs2).sign();
      newA[I+A.getHeight()][I+A.getWidth()]=theField.zHomomorphism(-1);//(sign>0)?-1:1);
      newW[I+A.getWidth()]=theField.zHomomorphism(-1);//(sign>0)?-1:1);
    }
  set<int> newBasis;

  //  for(int i=0;i<A.getHeight();i++)newBasis.insert(i);
  for(set<int>::const_iterator i=S.begin();i!=S.end();i++)newBasis.insert(*i);
  for(set<int>::const_iterator i=inequalitiesWithSlack.begin();i!=inequalitiesWithSlack.end();i++)newBasis.insert(*i);

  FieldLP ret(newA,concatenation(b,FieldVector(theField,inequalitiesWithSlack.size())));
  ret.setObjectiveFunction(newW);
  ret.setCurrentBasis(newBasis);

  return ret;
}


bool FieldLP::findFeasibleBasis()
{
  FieldLP A2=buildLPForFeasibilityCheck();

      AsciiPrinter Q(Stdout);
  int status;
  do
    {
      //A2.print(Q);
      status=A2.step();
    }
  while(status==1);
  fprintf(Stderr,status?"LP is unbounded.\n":"Optimal solution found.\n");
  //A2.print(Q);

  set<int> newBasis;
  for(set<int>::const_iterator i=A2.basis.begin();i!=A2.basis.end();i++)
    {
      if(*i<A.getHeight())
	newBasis.insert(*i);
    }

  if(newBasis.size()!=A.getWidth())return false;
  setCurrentBasis(newBasis);
  return true;
}


FieldLP FieldLP::withNoLineality()
{
  set<int> S;
  {
    FieldMatrix temp=A;
    temp.reduce();
    int i=-1;
    int j=-1;
    while(temp.nextPivot(i,j))S.insert(j);
  }
  //  FieldMatrix(theField,A.getHeight(),S.size());
  FieldLP ret(submatrix(A.transposed(),S).transposed(),b);

  //fprintf(Stderr,"New height=%i New width=%i\n",ret.A.getHeight(),ret.A.getWidth());

  //  AsciiPrinter P(Stderr);
  //ret.print(P);

  if(w.size())ret.setObjectiveFunction(subvector(w,S));

  return ret;
}
