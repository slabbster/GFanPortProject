#ifndef GFANAPPLICATION_H_INCLUDED
#define GFANAPPLICATION_H_INCLUDED

#include "application.h"
#include "vektor.h"

/** @mainpage The Gfan Source Code
This is the documentation of the Gfan header
 files. Please note that many cpp files do not have header files. In
 particular many Application files do not have header files and the
 best way to find documentation of these is through their help text
 which can be viewed when running the associated Gfan program with the
 option --help.

For the best overview of the Gfan classes go to the hierarchical view using the Classes/Class Hierarchy link above.
*/

class GFanApplication : public Application
{
 protected:
  class FieldOption: public IntegerOption
    {
    public:
      FieldOption();
      void onOptionsParsed();
    };
  class LogLevelOption: public IntegerOption
    {
    public:
      LogLevelOption();
      void onOptionsParsed();
    };
  class StdinFileOption: public StringOption
  {
  public:
	  StdinFileOption();
	  void onOptionsParsed();
  };
  class StdoutFileOption: public StringOption
  {
  public:
	  StdoutFileOption();
	  void onOptionsParsed();
  };
  class XmlOption: public SimpleOption
  {
  public:
    XmlOption();
    void onOptionsParsed();
  };
  void assertSymmetriesMatch(IntegerVectorList const &g, class PolynomialSet &markedGroebnerBasis, class FieldMatrix const *torusActions=0, bool asPolynomials=false);
  virtual void onExit();
  //  FieldOption theFieldOption;
  LogLevelOption theLogLevelOption;
  StdinFileOption theStdinFileOption;
  StdoutFileOption theStdoutFileOption;
  XmlOption theXmlOption;
};

#endif

