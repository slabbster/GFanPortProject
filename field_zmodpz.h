#ifndef FIELD_ZMODPZ_H_INCLUDED
#define FIELD_ZMODPZ_H_INCLUDED

#include "field.h"

class FieldZModPZImplementation : public FieldImplementation
{
  FieldElementImplementation *zHomomorphismImplementation(int n);
  FieldElement zHomomorphism(int n);

  std::string toString()const;
  const char *name();
  const int p;
 public:
  FieldZModPZImplementation(int p_);
  int getP()const;
};

class FieldZModPZ : public Field
{
 public:
  FieldZModPZ(int p_);
  int getP();
};

Field *theZMod2ZField();



#endif
