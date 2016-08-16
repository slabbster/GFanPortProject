#ifndef POLYMAKEFILE_H_INCLUDED
#define POLYMAKEFILE_H_INCLUDED

#include <list>
#include <vector>
#include <string>
#include <iostream>

#include "matrix.h"

using namespace std;

class PolymakeProperty
{
 public:
  std::string value;
  std::string name;
  bool embedded;
  PolymakeProperty(const string &name_, const string &value_, bool embedded_=true);
};

class PolymakeFile
{
  static bool xmlForced;
  string application,type;
  string fileName;
  list<PolymakeProperty> properties;
  list<PolymakeProperty>::iterator findProperty(const char *p);
  void writeProperty(const char *p, const string &data, bool embedded=false);
  bool isXml;
 public:
   bool isXmlFormat()const{return isXml;}
   /**
    * This function takes the Gfan property names and maps them to the corresponding names in
    * polymake type SymmetricFan (if the file is an xml file). This preserves compatibility with old
    * Gfan files. Maybe Gfan should switch completely to polymake names in the future.
    */
   const char *mapToPolymakeNames(const char *s);
/**
 * Calling this function will set the xmlForced flag, with the consequence that every subsequently
 * created PolymakeFile is output in XML format no matter with what parameters it is created.
 */
  static void forceXml(){xmlForced=true;}
  void open(const char *fileName_);
  void create(const char *fileName_, const char *application_, const char *type_, bool isXml_=false);
  void writeStream(ostream &file);
  void close();
  bool hasProperty(const char *p, bool doAssert=false);

  // The following could be part of a subclass to avoid dependencies on gfan
  int readCardinalProperty(const char *p);
  void writeCardinalProperty(const char *p, int n);

  bool readBooleanProperty(const char *p);
  void writeBooleanProperty(const char *p, bool n);

  IntegerMatrix readMatrixProperty(const char *p, int height, int width);
  void writeMatrixProperty(const char *p, const IntegerMatrix &m, bool indexed=false, const vector<string> *comments=0);

  IntegerMatrix readArrayArrayIntProperty(const char *p, int width);
  void writeArrayArrayIntProperty(const char *p, const IntegerMatrix &m);

  vector<list<int> > readMatrixIncidenceProperty(const char *p);
  void writeIncidenceMatrixProperty(const char *p, const vector<list<int> > &m, int baseSetSize);

  IntegerVector readCardinalVectorProperty(const char *p);
  void writeCardinalVectorProperty(const char *p, IntegerVector const &v);

  void writeStringProperty(const char *p, const string &s);
};

#endif
