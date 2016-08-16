#include "primarydecomposition.h"
#include "python2.6/Python.h"

extern "C" void minimal_associated_primes(PolynomialSetList *l, PolynomialRing const *r, PolynomialSet const *g);

#include <iostream>

class SageSingularPrimDec : public PrimaryDecompositionEngine
{
 public:
	 SageSingularPrimDec() : PrimaryDecompositionEngine()
	 {
		 cerr<<"constructing SageSingularPrimDec"<<endl;
		 if (!Py_IsInitialized()) {
			 Py_Initialize();
	   char* argv = {"gfan_python"};
	         PySys_SetArgv(1, &argv); // Sage needs argc, argv set for some reason

	         //Py_SetPythonHome("/home/burcin/sage/sage-4.3.2/local");
	         PyRun_SimpleString("import sys");
	         //PyRun_SimpleString("print sys.path");
	         PyRun_SimpleString("sys.path.append('/home/burcin/sage/sage-4.3.2/local/lib/python2.6/lib-dynload/')");
	         //PyRun_SimpleString("print sys.path");
	         cout<<"initialized python"<<endl;

	         PyRun_SimpleString("import sage.all");

	         PyObject* sage_module = PyImport_ImportModule("sage_link");
	         }
	 }
	 ~SageSingularPrimDec()
	 {
		 Py_Finalize();
			 }
  virtual const char *name()
  {
	  return "sagesingular";
  }
  virtual PolynomialSetList minimalAssociatedPrimes(PolynomialSet const &idealGenerators)
  {
//	  ring r=singularRing(idealGenerators.getRing());
//	  ideal singularPolynomialSet(idealGenerators);

	  cerr<<"CALLING SINGULAR THROUGH SAGE"<<endl;
	  PolynomialSetList ret;


	  PolynomialRing r=idealGenerators.getRing();
	  minimal_associated_primes(&ret,&r,&idealGenerators);
	  cerr<<"RETURNiNG FROM SAGE"<<endl;

	  return ret;
  }
};

static SageSingularPrimDec sageSingularPrimDec;
