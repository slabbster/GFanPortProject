#include "lll.h"

/*double lensq(vektor_f *a)
{
  vektor_f A1=(*a* *a);
  return A1.sum();
}
*/
int down(double f)
{
  if(f>=0)return int(f);
  return int(f-1);
}

void calcmy(IntegerMatrix const &b, FloatMatrix &my, Vektor<double> &B)
{
  FloatMatrix bs=integerToFloatMatrix(b);

  for(int k=0;k<b.getHeight();k++)
    {
      bs[k]=b[k];
      for(int j=0;j<k;j++)
	{
	  if(B[j]==0)
	    my[k][j]=0;
	  else
	    my[k][j]=dot(Vektor<double>(b[k]),bs[j])/B[j];

	  bs[k]=bs[k]-my[k][j]*bs[j];
	}
      B[k]=dot(bs[k],bs[k]);
    }
}

/*void calclambda(basis_i b, basis_f &lambda)
{
  for(int k=0;k<b.size();k++)
    {
      for(int j=0;j<=k;j++)
	{
	  double u=dot(b[k],b[j]);

	  for(int i=0;i<j;i++)
	    {
	      if(i-1==-1)
		{
		  u=lambda[i][i]*u-lambda[k][i]*lambda[j][i];
		}
	      else
		{
		  if(lambda[i-1][i-1]!=0)
		    {
		      u=lambda[i][i]*u-lambda[k][i]*lambda[j][i];
		      u/=lambda[i-1][i-1];
		    }
		}
	    }
	  lambda[k][j]=u;
	}
    }
}*/

IntegerMatrix mlll(IntegerMatrix &b, IntegerMatrix *inverseTransposedM)
{
  int k=1;
  int n=b.getHeight();// number of generators

  FloatMatrix my(n,n);//n*n

  Vektor<double> B(n);
  IntegerMatrix M=IntegerMatrix::identity(n);//n*n;
  IntegerMatrix MInverseTransposed=M;

  calcmy(b,my,B);
  while(k<n)
    {
      //size reduction
      for(int l=k-1;l>=0;l--)
	{
	  //	  calclambda(b,lambda);
	  int q=down(my[k][l]+0.5);
	  if(q)
	    {
	      b[k]=b[k]-q*b[l];
	      M[k]=M[k]-q*M[l];
	      MInverseTransposed[l]+=q*MInverseTransposed[k];
	      calcmy(b,my,B);
	    }
	}
      //      calcmy(b,my,B);
      //      calclambda(b,lambda);
      //test Lovasz' condition
      if(B[k]<(0.75-my[k][k-1]*my[k][k-1])*B[k-1])
	{
	  Vektor<double> temp=b[k];b[k]=b[k-1];b[k-1]=temp;
	  Vektor<double> temp1=M[k];M[k]=M[k-1];M[k-1]=temp1;
	  temp1=MInverseTransposed[k];MInverseTransposed[k]=MInverseTransposed[k-1];MInverseTransposed[k-1]=temp1;
      calcmy(b,my,B);
	  k--;
	  if(k<1)k=1;
	}
      else k++;
    }

  if(inverseTransposedM)*inverseTransposedM=MInverseTransposed;
  return M;
}


/*int rank(basis_i B)
{
  int n=B.size();
  basis_i M=mlll(B);
  int kerdim=0;
  while((kerdim<B.size())&&(B[kerdim].iszero()))kerdim++;

  return n-kerdim;
}
*/


IntegerMatrix latticeKernelOfTransposed(IntegerMatrix const &B)
{
  IntegerMatrix B2=B;
  IntegerMatrix M=mlll(B2);
  int kerdim=0;
  while((kerdim<B.getHeight())&&(B2[kerdim].isZero()))
    kerdim++;

  IntegerMatrix ret(kerdim,B.getHeight());

  for(int i=0;i<kerdim;i++)
    ret[i]=M[i];

  // add asserts to check that lll reduction works correctly

  return ret;
}
