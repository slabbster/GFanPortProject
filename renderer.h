#ifndef RENDERER_H_INCLUDED
#define RENDERER_H_INCLUDED

#include "polynomial.h"
#include <stdio.h>

class StandardMonomialRenderer
{
  FILE *f;
  int numberOfDrawingsPerLine;
  int maxEntry;
  int boxSize;
  int position;
  int getOffsetX();
  int getOffsetY();
  bool isInInitialIdeal(const PolynomialSet &s,int x, int y, int z);
  int trace(const PolynomialSet &s, int x, int y, int z, int d);
  void putPoint(float x, float y, int rotation, float size);
  void drawQuad(int x, int y, int z, int rotation, int color, int intensity=25);
 public:
  StandardMonomialRenderer(FILE *f);
  void render(const PolynomialSet &p);
  void setMaxEntry(int maxEntry);
  void setBoxSize(int boxSize);
  void setNumberOfDrawingsPerLine(int numberOfDrawingsPerLine);
};

#endif
