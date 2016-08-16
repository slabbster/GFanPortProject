/*
 * app_matrixproduct.cpp
 *
 *  Created on: Jan 4, 2011
 *      Author: anders
 */

#include "parser.h"
#include "printer.h"
#include "gfanapplication.h"
#include "matrix.h"

class MatrixProductApplication : public GFanApplication
{
public:
  SimpleOption optionTropical;
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This program computes the product of two matrices.\n";
  }
  MatrixProductApplication():
    optionTropical("--tropical","Do the computation in the max-plus semi-ring.")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_matrixproduct";
  }

  int main()
  {
    FileParser P(Stdin);

    IntegerMatrix A=rowsToIntegerMatrix(P.parseIntegerVectorList());
    IntegerMatrix B=rowsToIntegerMatrix(P.parseIntegerVectorList());

    pout<<((optionTropical.getValue())?tropicalProduct(A,B):A*B).getRows();

    return 0;
  }
};

static MatrixProductApplication theApplication;
