#include <iostream>
#include <stdlib.h>
#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "division.h"
#include "buchberger.h"
#include "wallideal.h"
#include "lp.h"
#include "reversesearch.h"
#include "termorder.h"
#include "ep_standard.h"
#include "ep_xfig.h"
#include "polyhedralcone.h"
#include "gfanapplication.h"
#include "saturation.h"
#include "field_rationals.h"
#include "field_zmodpz.h"
#include "field_rationalfunctions.h"
#include "symmetry.h"
#include "linalg.h"
#include "fieldlp.h"

class TestApplication : public GFanApplication
{
  StringOption testSuiteFolderOption;
  StringOption executableOption;
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This runs the test suite and checks against the stored result. If no result exists, it is generated.\n";
  }
  TestApplication():
	testSuiteFolderOption("--suite","Specify the folder which contains the test suite.","testsuite"),
	executableOption("--gfan","Specify name of gfan executable to test.","./gfan")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_test";
  }

  void lpRationalFunctionTest()
  {
    //int n=3;
	//     IntegerVectorList L=StringParser("{(1,0,2),(1,2,0)}").parseIntegerVectorList();
    //  IntegerVectorList L=StringParser("{(1,2,0)}").parseIntegerVectorList();

    //     int n=4;
    //  IntegerVectorList L=StringParser("{(1,0,2,3),(1,2,3,0)}").parseIntegerVectorList();

    /*     int n=16;
     IntegerVectorList L=StringParser("{(1,0,2,3,5,4,6,7,9,8,10,11,13,12,14,15),"
				      "(3,0,1,2,7,4,5,6,11,8,9,10,15,12,13,14),"
				      "(0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15)}").parseIntegerVectorList();
    */
    int n=9;
     IntegerVectorList L=StringParser("{(1,0,2,4,3,5,7,6,8),"
				      "(1,2,0,4,5,3,7,8,6),"
				      "(3,4,5,0,1,2,6,7,8),"
				      "(3,4,5,6,7,8,0,1,2)"
				      "}").parseIntegerVectorList();
     /*
    int n=9;
     IntegerVectorList L=StringParser("{(1,0,2,4,3,5,7,6,8),"
				      "(1,2,0,4,5,3,7,8,6),"
				      "(0,3,6,1,4,7,2,5,8)}").parseIntegerVectorList();
     */

    SymmetryGroup s(n);
    s.computeClosure(L);

    FieldRationalFunctions F(Q,"t");
    FieldMatrix M(F,s.elements.size(),n);
    int I=0;
    for(SymmetryGroup::ElementContainer::const_iterator i=s.elements.begin();i!=s.elements.end();i++,I++)
      {
	for(int j=0;j<n;j++)
	  {
	    M[I][j]=M[I][j]+F.exponent(j);
	    M[I][j]=M[I][j]-F.exponent(((SymmetryGroup::inverse(*i))[j]));
	  }
      }
    AsciiPrinter P(Stderr);
    //    M.printMatrix(P);


    FieldMatrix M2(Q,M.getHeight(),M.getWidth());
    for(int i=0;i<M.getHeight();i++)
      for(int j=0;j<M.getWidth();j++)
	M2[i][j]=F.substitute(M[i][j],Q.zHomomorphism(10000).inverse());

        Field &myField=F;
	{
      M=M2;
      myField=Q;//HERE
      }
    IntegerVectorList extreme;

    M.printMatrix(P);


    I=0;
    FieldElement minusOne=myField.zHomomorphism(-1);//HERE
    for(SymmetryGroup::ElementContainer::const_iterator i=s.elements.begin();i!=s.elements.end();i++,I++)
      //    for(int i=0;i<s.elements.size();i++)
      {
	fprintf(Stderr,"IIII:%i\n",I);
	M[I]=minusOne*M[I];
	FieldVector b(myField,s.elements.size());//HERE
	b[I]=minusOne;
	FieldLP lp(M,b);
	FieldLP lp2=lp.withNoLineality();

	AsciiPrinter P(Stderr);
	//lp2.print(P);

	//	cerr <<"Result" << lp2.findFeasibleBasis()<<endl;

	if(lp2.findFeasibleBasis())
	  extreme.push_back(*i);

	M[I]=minusOne*M[I];
      }
    fprintf(Stdout,"Extreme permutations (%i):\n",extreme.size());
    P.printVectorList(extreme);
  }


  class TestCase
  {
  public:
	  string folder;
	  TestCase(string const &folder_):
		  folder(folder_)
		  {

		  }
	  void fail()
	  {
		  cerr<<"Test failed:"<<folder<<endl;
		  assert(0);
	  }
	  void compare(string a, string b)
	  {
		  FILE *A=fopen(a.c_str(),"r");
		  FILE *B=fopen(b.c_str(),"r");
		  assert(A);
		  assert(B);
		  while((!feof(A))&&(!feof(B)))
		  {
			  if(fgetc(A)!=fgetc(B))fail();
		  }
		  if(feof(A)!=feof(B))fail();
	  }
	  bool fileExists(string name)
	  {
		  FILE *f=fopen(name.c_str(),"r");
		  if(f)fclose(f);
		  return f;
	  }
	  /*
	   * Returns true if test was successful.
	   * Returns false if test was performed for the first time.
	   * Asserts if test fails
	   */
	  bool perform(const char *exe)
	  {
		  string fileName=folder+"/command";
		  FILE *f=fopen(fileName.c_str(),"r");
		  if(!f)
		  {
			  cerr<<"Could not open file:\""<<fileName<<"\""<<endl;
			  assert(f);
		  }
		  char command[4096];
		  char *temp=fgets(command,4095,f);
		  fclose(f);
		  assert(temp);
		  for(int i=0;i<4096 && command[i];i++)if(command[i]=='\n'){command[i]=0;}
		  char command2[4096];
		  string input=folder+"/input";
		  sprintf(command2,command,exe,exe,exe,exe);
		  bool outputExists=fileExists(folder+"/output");
		  string outputName=folder+"/output";
		  if(outputExists)
		  {
			  outputName=outputName+"New";
		  }
		  {
			  string t="rm "+outputName;
			  system(t.c_str());
		  }
		  string command3="cat <"+input+"|"+string(command2)+">"+outputName;
		  cerr<<"Running command:\""<<command3<<"\""<<endl;
		  system(command3.c_str());
		  if(outputExists)compare(folder+"/output",folder+"/outputNew");
		  return outputExists;
	  }
  };

  list<string> subFolderNames()
  {
#define tempName "GfAnTeMpTeStS"
	  char command[256];
	  system("rm "tempName);
	  sprintf(command,"ls %s>" tempName ,testSuiteFolderOption.getValue());
	  system(command);

	  list<string> ret;
	  FILE *f=fopen(tempName,"r");
	  assert(f);
	  char name[256];
	  while(fgets(name,255,f))
	  {
		  for(int i=0;i<255 && name[i];i++)if(name[i]=='\n'){name[i]=0;}
		  if(name[0]>='0' && name[0]<='9')ret.push_back(string(testSuiteFolderOption.getValue())+"/"+string(name));
	  }
	  fclose(f);
	  return ret;
  }

  int main()
  {
//    lpRationalFunctionTest();
//    testRationalFunctionField();

	  list<string> testFolders=subFolderNames();
	  list<TestCase> testList;
	  for(list<string>::const_iterator i=testFolders.begin();i!=testFolders.end();i++)
	  {
		  testList.push_back(TestCase(*i));
	  }

	  cout<<"Number of tests to perform "<<testList.size()<<endl;

	  int good=0;
	  int bad=0;
	  for(list<TestCase>::iterator i=testList.begin();i!=testList.end();i++)
		  if(i->perform(executableOption.getValue()))
			  good++;
		  else
			  bad++;
	  cout<<"Number of succesful tests "<<good<<endl;
	  cout<<"Number of initialized tests "<<bad<<endl;

	  return 0;
  }
};

static TestApplication theApplication;
