#include "parser.h"
#include "printer.h"
#include "lp.h"
#include "gfanapplication.h"
#include "polyhedralcone.h"
#include "polyhedralfan.h"
#include "polymakefile.h"
#include "symmetriccomplex.h"
#include "lll.h"
#include "linalg.h"
#include "log.h"

/*
 * About the homology computation:
 * The boundary map matrices are computed by SymmetricComplex, but since these are very large
 * this program does not store the complete matrices in the -o mode. Instead sparse matrices are used.
 * For each homology group we need to consider two boundary map matrices B2 and B1. We have to compute the kernel
 * of B1 modulo the image of B2. One way to do this is by first computing a lattice basis of the
 * kernel of B1 and then write the columns of B2 in new coordinates. In these new coordinates
 * we check how much of the space is spanned by computing an integer row echelon form of the transposed matrix
 * of generators. The quotient can now be read of from the pivots, while the number of non-pivot rows
 * are the Betti numbers.
 * In the optimized -o version, we instead start by investigating the second matrix B1 which we need to
 * compute the kernel of. It often turns out that some coordinates of a vector in the kernel are implicit
 * in the sense that given the others we can determine what they are. This happens every time we have a
 * non-zero entry in the matrix. We can then erase that column and row of the matrix if we update the
 * other columns of B1 accordingly. In an attempt to avoid coefficient growth we only pick an entry which
 * is +/-1 when removing its rows and columns. If we were not worried about coefficient growth we could
 * repeat this process until B1 has the full space as a kernel - but we avoid this in the current
 * implementation. We update B2 by removing the rows of B2 which correspond to the removed columns of B1.
 * We now row reduce B1 transposed - or at least we do this partly by identifying possible pivots with
 *
 */

class SparseMatrix
{
	/* We choose a cache friendly vector for storing a column:
	 * even entries refer to a row index, while odd entries are values.
	 * The pairs are kept sorted.
	 */
	vector<IntegerVector> columns;
	set<int> remainingRows;
	set<int> remainingColumns;
public:
	SparseMatrix(int height, int width):
		columns(width)
		{
				for(int i=0;i<height;i++)remainingRows.insert(i);
				for(int j=0;j<width;j++)remainingColumns.insert(j);
		}
	void assignColumn(int col, vector<int> const &indices, vector<int> const &values)
	{
		IntegerVector temp(indices.size()*2);
		set<pair<int,int> > s;
		for(int k=0;k<indices.size();k++)
		{
			s.insert(pair<int,int>(indices[k],values[k]));
		}

		int k=0;
		for(set<pair<int,int> >::const_iterator i=s.begin();i!=s.end();i++,k++)
		{
			temp[k*2]=i->first;
			temp[k*2+1]=i->second;
		}
		columns[col]=temp;
	}
	set<int> getRemainingColumns()
	{
		return remainingColumns;
	}
	IntegerMatrix remainingToIntegerMatrix()
	{
		int height=remainingRows.size();
		int width=remainingColumns.size();
		IntegerMatrix ret(height,width);

		int J=0;
		for(set<int>::const_iterator j=remainingColumns.begin();j!=remainingColumns.end();j++,J++)
		{
			int a=0;
			int I=0;
			for(set<int>::const_iterator i=remainingRows.begin();i!=remainingRows.end();i++,I++)
			{
				if(a>=columns[*j].size())break;
				if(columns[*j][a]==*i)
					{
						ret[I][J]=columns[*j][a+1];
						a+=2;
					}
			}
		}

		return ret;
	}
	/* Returns -1 if entry does not exist. Otherwise index into column data (multiplied by 2 already)*/
	int lookUpIndex(int i, int j)
	{
		for(int k=0;k<columns[j].size();k+=2)
			if(columns[j][k]==i)return k;
		return -1;
/*		int hi=columns[j].size();
		int lo=0;
		while(hi>lo)
		{
			int med=((hi>>1)+(lo>>1))&(-2);
			if(columns[j][med]<i)
			{

			}
		}
*/
		}
	int lookUp(int i, int j)
	{
		int k=lookUpIndex(i,j);
		if(k!=-1)
			return columns[j][k+1];
		return 0;
	}
	void eraseRow(int i)
	{
		for(set<int>::const_iterator j=remainingColumns.begin();j!=remainingColumns.end();j++)
		{
			int k=lookUpIndex(i,*j);
			if(k!=-1)
			{
				columns[*j]=concatenation(columns[*j].subvector(0,k),columns[*j].subvector(k+2,columns[*j].size()));
			}
		}
		remainingRows.erase(i);
	}
	void eraseColumn(int col)
	{
		assert(col<columns.size());
		assert(col>=0);
		remainingColumns.erase(col);
		columns[col]=IntegerVector();
	}
	void multiplyAddColumn(int source, int dest, int multiplier)
	{
		if(multiplier==0)return;
		IntegerVector temp;
		int si=0;
		int di=0;
		while(si<columns[source].size() || di<columns[dest].size())
		{
			if(si==columns[source].size())
			{
				temp.grow(temp.size()+2);
				temp[temp.size()-2]=columns[dest][di];
				temp[temp.size()-1]=columns[dest][di+1];
				di+=2;
			}
			else if(di==columns[dest].size())
			{
				temp.grow(temp.size()+2);
				temp[temp.size()-2]=columns[source][si];
				temp[temp.size()-1]=multiplier*columns[source][si+1];
				si+=2;
			}
			else
			{
				if(columns[source][si]<columns[dest][di])
				{
					temp.grow(temp.size()+2);
					temp[temp.size()-2]=columns[source][si];
					temp[temp.size()-1]=multiplier*columns[source][si+1];
					si+=2;
				}
				else if(columns[source][si]>columns[dest][di])
				{
					temp.grow(temp.size()+2);
					temp[temp.size()-2]=columns[dest][di];
					temp[temp.size()-1]=columns[dest][di+1];
					di+=2;
				}
				else
				{
					int sum=columns[dest][di+1]+multiplier*columns[source][si+1];
					if(sum<-32000 || sum>32000)
					{
						cerr<<"Overflow in homology computation"<<endl;
						assert(0);
					}
					if(sum!=0)
					{
						temp.grow(temp.size()+2);
						temp[temp.size()-2]=columns[dest][di];
						temp[temp.size()-1]=sum;
					}
					di+=2;
					si+=2;
				}
			}
		}
		columns[dest]=temp;
	}
	void debug()
	{
		AsciiPrinter(Stderr)<<remainingToIntegerMatrix().getRows();
		for(set<int>::const_iterator i=remainingRows.begin();i!=remainingRows.end();i++)cerr<<*i;cerr<<endl;
		for(set<int>::const_iterator i=remainingColumns.begin();i!=remainingColumns.end();i++)cerr<<*i;cerr<<endl;
		for(vector<IntegerVector>::const_iterator i=columns.begin();i!=columns.end();i++)AsciiPrinter(Stderr)<<*i;
	}
	void killEntry(int i, int j)
	{
		int v=lookUp(i,j);
		assert(v==-1 || v==1);
		for(set<int>::const_iterator a=remainingColumns.begin();a!=remainingColumns.end();a++)
			if(*a!=j)
			{
				multiplyAddColumn(j, *a, -v*lookUp(i,*a));
			}
		eraseColumn(j);
		eraseRow(i);
		log3 debug();
	}
	bool getEntryToKill(int &reti, int &retj)
	{
		for(set<int>::const_iterator i=remainingColumns.begin();i!=remainingColumns.end();i++)
		{
			for(int j=0;j<columns[*i].size();j+=2)
				if(columns[*i][j+1]==-1 || columns[*i][j+1]==1)
				{
					retj=*i;
					reti=columns[*i][j];
					return true;
				}
		}
		return false;
	}
};


class FanHomologyApplication : public GFanApplication
{
  StringOption inputOption;
  SimpleOption noOptimizeOption;
public:
  const char *helpText()
  {
    return "This program takes a polyhedral fan and computes its reduced homology groups. "
          "Of course the support of a fan is contractible, so what is really computed is the "
    "reduced homology groups of the support of the fan after quotienting out with the lineality space and intersecting with a sphere. "
    "Notice that taking the quotient with the lineality space results in an inverted suspension which just results in a shift of the reduced homology groups.\n";
  }
  FanHomologyApplication():
    inputOption("-i","Specify the name of the input file.","polymake.out"),
    noOptimizeOption("--no-optimize","Disable preprocessing of boundary maps before doing lattice computations.")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_fanhomology";
  }

  void printHomologyGroup(Printer &P, FieldMatrix const &newB)
  {
		bool first=true;
		int i=-1;
		int j=-1;
		int lasti=-1;
		while(newB.nextPivot(i,j))
		{
			int a=toInteger(newB[i][j]);
			if(a<0)a=-a;
			if(a!=1)
				{
					if(!first)
					{
						P<<"x";
					}
					P<<"(Z/";
					P.printInteger(a);
					P<<"Z)";
					first=false;
				}
			lasti=i;
		}
		int quotientRank=newB.nonPivotColumns().size();
//		int quotientRank=newB.getHeight()-lasti-1;
		if(quotientRank>0)
		{
			if(!first)P<<"x";
			P<<"Z^";
			P.printInteger(quotientRank);
			first=false;
		}
		if(first)P<<"{0}";
		P<<"\n";
	}
  int main()
  {
    PolyhedralFan f=PolyhedralFan::readFan(inputOption.getValue());

    AsciiPrinter P(Stdout);
    AsciiPrinter PP(Stderr);

    SymmetricComplex s=f.toSymmetricComplex(0);

    log1 PP<<s.fvector();

    int linealityDim=f.dimensionOfLinealitySpace();
    int d=f.getMaxDimension();

	log1 PP<<"\n------------------------\n";
    for(int i=d;i>linealityDim;i--)
    {
    	log1 PP<<"H_";
    	log1 PP.printInteger(i-linealityDim-1);
    	log1 PP<<"()\n";

    	IntegerMatrix b1,b2;
    	if(noOptimizeOption.getValue())
    	{
    		b1=s.boundaryMap(i);
    		b2=s.boundaryMap(i+1);

    		log1{
    			PP<<"Boundary map ";
    		PP.printInteger(i-linealityDim);
    		PP<<"\n"<<b2.getRows();
    		PP<<"Boundary map ";
    		PP.printInteger(i-linealityDim-1);
    		PP<<"\n"<<b1.getRows();
    		//P<<m.getRows();
    		}
        	FieldMatrix kernel=integerMatrixToFieldMatrix(latticeKernelOfTransposed(b1.transposed()),Q);
        	kernel.reduce(false,true);
        	kernel.removeZeroRows();
        	log1 PP<<"Kernel\n"<<kernel;
        	log1 PP.printInteger(kernel.getHeight());
        	log1 PP.printInteger(kernel.getWidth());
        	FieldMatrix image=integerMatrixToFieldMatrix(b2,Q).transposed();
        	image.reduce(false,true);
        	image.removeZeroRows();
        	log1 PP<<"Image\n"<<image;
        	log1 PP.printInteger(image.getHeight());
        	log1 PP.printInteger(image.getWidth());

        	FieldMatrix combined=combineLeftRight(kernel.transposed(),image.transposed());
        	combined.reduce(false,true);
        	combined.REformToRREform(false);
        	combined.removeZeroRows();
        	FieldMatrix newB=combined.submatrix(0,combined.getHeight(),kernel.getHeight(),combined.getWidth());
        	newB=newB.transposed();

        	newB.reduce(false,true);
        	log1 PP<<"In new basis\n:"<<newB;
        	P<<"H_";
        	P.printInteger(i-linealityDim-1);
        	P<<"()=";
        	printHomologyGroup(P,newB);
    	}
    	else
    	{
//    		assert(s.sym.isTrivial());
    		int d=i;
    		int width=s.numberOfConesOfDimension(d);
    		int height=s.numberOfConesOfDimension(d-1);
    		SparseMatrix B1( height, width);

    		{
    			for(SymmetricComplex::ConeContainer::const_iterator i=s.cones.begin();i!=s.cones.end();i++)

    				if(d==i->dimension)
    				{
    					int I=s.dimensionIndex(*i);
    					vector<int> indices;
    					vector<int> signs;
    					s.boundary(*i,indices,signs);
    					B1.assignColumn(I, indices, signs);
    				}
    		}
			log3 AsciiPrinter(Stderr)<<B1.remainingToIntegerMatrix().getRows();
    		{
    			{
    				int i,j;

    			while(B1.getEntryToKill(i,j))
    				{
						log3 AsciiPrinter(Stderr)<<B1.remainingToIntegerMatrix().getRows();
						B1.killEntry(i,j);
    				}
    			}
    			if(1)//kill even more columns in B1. At the end we are left with an all zero matrix.
    			{
	    			map<int,int> remap;
	    			int L=0;
    	    		set<int> remainingColumns=B1.getRemainingColumns();
	    			for(set<int>::const_iterator l=remainingColumns.begin();l!=remainingColumns.end();l++,L++)remap[L]=*l;
	       			FieldMatrix BRat=integerMatrixToFieldMatrix(B1.remainingToIntegerMatrix(),Q);
	       			BRat.reduce(false,false);//no need for Z arithmetics
	       			{
	       				int i=-1;
	       				int j=-1;
	       				while(BRat.nextPivot(i,j))
	       				{
	       					B1.eraseColumn(remap[i]);
	       				}
	       			}
	       			//b1=IntegerMatrix(0,B1.getRemainingColumns().size());//We don't have to make the number of rows match
    			}
    			else
    			{
    			//	b1=B1.remainingToIntegerMatrix();
    			}
    			{
    	    		int d=i+1;
    	    		int width=s.numberOfConesOfDimension(d);
    	    		set<int> remainingColumns=B1.getRemainingColumns();
    	    		int height=remainingColumns.size();
    	    		if(0)
    	    		{
    	    			b2=IntegerMatrix( height, width);
    	    			for(SymmetricComplex::ConeContainer::const_iterator i=s.cones.begin();i!=s.cones.end();i++)
    	    				if(d==i->dimension)
    	    				{
    	    					int I=s.dimensionIndex(*i);
    	    					vector<int> indices;
    	    					vector<int> signs;
    	    					s.boundary(*i,indices,signs);
    	    					for(int k=0;k<indices.size();k++)
    	    						if(remainingColumns.count(indices[k]))
    	    						{
    	    							int newIndex=0;
    	    							for(set<int>::const_iterator l=remainingColumns.begin();l!=remainingColumns.end();l++,newIndex++)
    	    								if(*l==indices[k])break;
											b2[newIndex][I]=signs[k];
    	    						}
    	    				}
    	    		}
    	    		else
    	    		{
    	    			SparseMatrix B2(height,width);
    	    			for(SymmetricComplex::ConeContainer::const_iterator i=s.cones.begin();i!=s.cones.end();i++)
    	    				if(d==i->dimension)
    	    				{
    	    					int I=s.dimensionIndex(*i);
    	    					vector<int> indices;
    	    					vector<int> signs;
    	    					s.boundary(*i,indices,signs);
    	    					vector<int> indices2;
    	    					vector<int> signs2;
    	    					for(int k=0;k<indices.size();k++)
    	    						if(remainingColumns.count(indices[k]))
    	    						{
    	    							int newIndex=0;
    	    							for(set<int>::const_iterator l=remainingColumns.begin();l!=remainingColumns.end();l++,newIndex++)
    	    								if(*l==indices[k])break;
										indices2.push_back(newIndex);
										signs2.push_back(signs[k]);
    	    						}
    	    					B2.assignColumn(I,indices2,signs2);
    	    				}
    	    			map<int,int> remap;
    	    			int L=0;
    	    			for(set<int>::const_iterator l=remainingColumns.begin();l!=remainingColumns.end();l++,L++)remap[L]=*l;
    	    			{
    	    				int i,j;

    	    				while(B2.getEntryToKill(i,j))
    	    				{
    	    					log3 AsciiPrinter(Stderr)<<B2.remainingToIntegerMatrix().getRows();
    	    					B2.killEntry(i,j);
    	    					//B1.eraseColumn(remap[i]);
    	    				}
    	    			}
    	    			b2=B2.remainingToIntegerMatrix();
    	    			//b1=B1.remainingToIntegerMatrix();
    	    		}
    			}
    		}
    		log1{
    			PP<<"Boundary map restricted to \"non-basis\" coordinates";
//    		PP.printInteger(i-linealityDim);
    		PP<<"\n"<<b2.getRows();
    		//P<<m.getRows();
    		}
        	FieldMatrix image=integerMatrixToFieldMatrix(b2,Q).transposed();
        	image.reduce(false,true);
        	image.removeZeroRows();
        	log1 PP<<"Image\n"<<image;
        	log1 PP.printInteger(image.getHeight());
        	log1 PP.printInteger(image.getWidth());

        	P<<"H_";
        	P.printInteger(i-linealityDim-1);
        	P<<"()=";
        	printHomologyGroup(P,image);
    	}
    	log1 PP<<"\n------------------------\n";
    }

    return 0;
  }
};

static FanHomologyApplication theApplication;
