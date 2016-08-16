#ifndef MACAULAY2_H_INCLUDED
#define MACAULAY2_H_INCLUDED

#include "printer.h"
#include "vektor.h"
#include "decompose.h"

class Pipe
{
  int fdM2Input[2];
  int fdM2Output[2];
protected:
  FILE *pipeInput;
  FILE *pipeOutput;
public:
  Pipe();
  ~Pipe();
};


class Macaulay2Pipe : public Pipe
{
  AsciiPrinter printer;
public:
  Macaulay2Pipe();
  ~Macaulay2Pipe();

  void skipStartOfLine();
  char *readLine();
  int readInt();
  bool readBool();
  void setPolynomialRing(int numberOfVariables);
  void setPolynomialRing(gbasis const &g);

  int getPdimCokerGensMonomial(gbasis const &monomialIdeal);
  bool isHomogeneousIdeal(gbasis &ideal);
  StandardPairList standardPairs(gbasis &monomialIdeal);
};

#endif
