#include "tropicaldeterminant.h"

#include <algorithm>

#include "printer.h"
#include "log.h"

using namespace std;

static int sub(int a, int b)
{//a-b
	assert(b!=td_minusInfinity);
	if(a==td_minusInfinity)return td_minusInfinity;
	return a-b;
}
static int max(int a, int b)
{
	if(a==td_minusInfinity)return b;
	if(b==td_minusInfinity)return a;
	if(a>b)return a;
	return b;
}

/*
 * Computes the tropical max determinant using the Hungarian method.
 * The implementation follows notes by Golin on Bipartite matchings.
*/
int tropicalDeterminant(IntegerMatrix const &m_)
{
	IntegerMatrix m=m_;
	int d=m_.getHeight();
	int n=m_.getWidth();

	// make a square matrix by introducing zeros
	if(d<=n)
	{
		m=m_;
		m.append(IntegerMatrix(n-d,n));
		d=n;
	}
	else
	{
		m=m_.transposed();
		m.append(IntegerMatrix(d-n,d));
		n=d;
	}

	int total=(n*n-n)/2;

	IntegerVector labelingRows(d);
	IntegerVector labelingColumns(n);
	for(int j=0;j<n;j++)
	{
		int colmax=td_minusInfinity;
		for(int i=0;i<d;i++)colmax=max(m[i][j],colmax);
		labelingColumns[j]=colmax;
		if(colmax==td_minusInfinity)return td_minusInfinity;
	}

	// choose trivial matching M
	IntegerVector M(n);
	for(int i=0;i<n;i++)M[i]=-1;
	IntegerVector Minv(d);
	for(int i=0;i<d;i++)Minv[i]=-1;


	while(M.sum()!=total)
	{
		log3 debug<<"M"<<M<<"\n";
		log3 debug<<"Minv"<<Minv<<"\n";
		log3 debug<<"Lrows"<<labelingRows<<"\n";
		log3 debug<<"Lcolumns"<<labelingColumns<<"\n";

		int chosenRow=0;
		while(chosenRow<d){if(Minv[chosenRow]==-1)break;chosenRow++;}
		log3 debug<<"ChosenRow"<<chosenRow<<"\n";
		IntegerVector S(d);//-1=not present, n=root
		for(int k=0;k<d;k++)S[k]=-1;
		S[chosenRow]=d;
		IntegerVector T(n);//-1=not present
		for(int k=0;k<n;k++)T[k]=-1;

		l:
		log3 debug<<"looping\n";
		log3 debug<<"T"<<T<<"\n";
		log3 debug<<"S"<<S<<"\n";
		// compute NLS
		IntegerVector NLS(n);
		for(int k=0;k<n;k++)NLS[k]=-1;
		for(int i=0;i<d;i++)
		{
			if(S[i]!=-1)
			for(int j=0;j<n;j++)
			{
				if(NLS[j]==-1)if(labelingRows[i]+labelingColumns[j]==m[i][j])NLS[j]=i;
			}
		}

		log3 debug<<"NLS"<<NLS<<"\n";
		// find element in difference
		int chosenColumn=0;
		while(chosenColumn<n)
		{
			if(NLS[chosenColumn]!=-1 && T[chosenColumn]==-1)
				break;
			chosenColumn++;
		}
		log3 debug<<"ChosenCol"<<chosenColumn<<"\n";
		if(chosenColumn==n)
		{
			// if no element in difference then	update labels
			log3 debug<<"Updating labels\n";
			int minusal=td_minusInfinity;
			for(int i=0;i<d;i++)if(S[i]!=-1)
				for(int j=0;j<n;j++)if(T[j]==-1)
					minusal=max(minusal,sub(sub(m[i][j],labelingRows[i]),labelingColumns[j]));
			log3 debug<<"M"<<M<<"\n";
			log3 debug<<"Minv"<<Minv<<"\n";
			log3 debug<<"Lrows"<<labelingRows<<"\n";
			log3 debug<<"Lcolumns"<<labelingColumns<<"\n";
			log3 debug<<"T"<<T<<"\n";
			log3 debug<<"S"<<S<<"\n";

			log3 debug<<"Minus al"<<minusal<<"\n";
			if(minusal==td_minusInfinity)return td_minusInfinity;//is this right?
			for(int i=0;i<d;i++)if(S[i]!=-1)labelingRows[i]+=minusal;
			for(int i=0;i<n;i++)if(T[i]!=-1)labelingColumns[i]-=minusal;

			log3 debug<<"Lrows"<<labelingRows<<"\n";
			log3 debug<<"Lcolumns"<<labelingColumns<<"\n";

			goto l;
		}
		else
		{
			if(M[chosenColumn]==-1)
			{
				log3 debug<<"Updating using alternating path\n";
				log3 debug<<"T"<<T<<"\n";
				log3 debug<<"S"<<S<<"\n";
				log3 debug<<"chosenrow"<<chosenRow<<"\n";
				// we have found alternating path... update

				T[chosenColumn]=NLS[chosenColumn];
				//				while(n!=S[NLS[chosenColumn]])
					while(n!=S[T[chosenColumn]])
				{
					int c=chosenColumn;
//					int b=NLS[c];
					int b=T[c];
					int a=S[b];
					//add path b->c
					M[c]=b;
					Minv[b]=c;
					chosenColumn=a;
				}
				M[chosenColumn]=chosenRow;
				Minv[chosenRow]=chosenColumn;
				log3 debug<<"Done updating\n";
			}
			else
			{
				int z=M[chosenColumn];
				T[chosenColumn]=NLS[chosenColumn];//Minv[z];//or chosen row?
				S[z]=chosenColumn;
				goto l;
			}
		}
	}

	log3 debug<<"M"<<M<<"\n";
	int ret=0;
	for(int i=0;i<n;i++)
		ret+=m[M[i]][i];
	return ret;
}

int tropicalDeterminantSlow(IntegerMatrix const &m_)
{
	IntegerMatrix m=m_;
	if(m.getHeight()>m.getWidth())m=m.transposed();
	int d=m.getHeight();
	int n=m.getWidth();
	int ret=td_minusInfinity;

	list<int> perm;
	for(int i=0;i<d;i++)perm.push_back(i);
	for(int i=d;i<n;i++)perm.push_back(d);

	do
	{
		int prod=0;
		int i=0;
		for(list<int>::const_iterator I=perm.begin();I!=perm.end();i++,I++)if(*I!=d)
		{if(m[*I][i]==td_minusInfinity){prod=td_minusInfinity;break;}prod+=m[*I][i];}
		ret=max(ret,prod);
	}while(next_permutation(perm.begin(),perm.end()));
	return ret;
}


void tropicalDeterminantTest()
{
	for(int i=0;i<100000;i++)
	{
		debug<<i<<"\n";
		int d=((unsigned int)rand())%6;
		int n=((unsigned int)rand())%6;
		IntegerMatrix m(d,n);
		for(int y=0;y<d;y++)
			for(int x=0;x<n;x++)
			{
				if(rand()&1)
					m[y][x]=td_minusInfinity;
				else
					m[y][x]=rand()%16;
			}
		if(tropicalDeterminantSlow(m)!=tropicalDeterminant(m))
		{
			debug<<m.getRows();
			debug<<tropicalDeterminantSlow(m)<<"\n";
			debug<<tropicalDeterminant(m)<<"\n";
			assert(0);
		}
	}
}
