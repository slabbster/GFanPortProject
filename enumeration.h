#ifndef ENUMERATION_H_INCLUDED
#define ENUMERATION_H_INCLUDED

#include "polynomial.h"

class EnumerationTarget
{
 public:
  virtual void beginEnumeration(const PolynomialSet &groebnerBasis)=0;
  virtual void endEnumeration()=0;
  virtual bool basis(const PolynomialSet &groebnerBasis)=0; /* return false to break enumaration */
};


class EnumerationTargetCollector: public EnumerationTarget
{
 public:
  PolynomialSetList theList;
  void beginEnumeration(const PolynomialSet& g);
  void endEnumeration();
  bool basis(const PolynomialSet &groebnerBasis);
  PolynomialSetList getList();
};


class EnumerationPrinter: public EnumerationTarget
{
};


class EnumerationFilePrinter: public EnumerationPrinter
{
  FILE *initialisedFile;
 protected:
  string filename;
  FILE *file;
 public:
  EnumerationFilePrinter();
  ~EnumerationFilePrinter();

  void open(std::string filename);
  void open(FILE *file);
  void close();

  virtual void onOpened()=0;
  virtual void onClose()=0;
  virtual void onClosed()=0;

  virtual string extension();
};


class EnumerationAlgorithm
{
 protected:
  EnumerationTarget *target;
  int progressCounter;
  virtual void printProgress(int step=1);
  void targetBeginEnumeration(const PolynomialSet &groebnerBasis){if(target)target->beginEnumeration(groebnerBasis);}
  void targetEndEnumeration(){if(target)target->endEnumeration();}
  bool targetBasis(const PolynomialSet &groebnerBasis){bool ret=true;if(target)ret=target->basis(groebnerBasis);printProgress();return ret;}
 public:
  EnumerationAlgorithm(){target=0;progressCounter=0;}
  void setEnumerationTarget(EnumerationTarget *target){this->target=target;}
  virtual void enumerate(const PolynomialSet &groebnerBasis){}
};

// The following is used to glue the old code to the new code
#include "symmetrictraversal.h"

class TargetGlue:public SymmetricTarget
{
	EnumerationTarget &target;
public:
	TargetGlue(EnumerationTarget &target_);
	bool process(ConeTraverser &traverser);
};

#endif
