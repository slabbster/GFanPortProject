#ifndef TERMORDER_H_INCLUDED
#define TERMORDER_H_INCLUDED

class TermOrder;

#include "vektor.h"
class Printer;
//#include "printer.h"

// All term orders must be represented by a matrix. This is needed in the generic Groebner walk.

class TermOrder
{
 public:
  virtual int rowDot(int row, const IntegerVector &v)const=0;
  virtual bool operator()(const IntegerVector &a, const IntegerVector &b, int scaleA=1, int scaleB=1, int perturbationDegree=-1)const=0;
  //  virtual bool operator()(IntegerVector const &a, IntegerVector const &b)const=0;
  virtual void print(Printer &p)const;
  void printMatrix(Printer &p, int dim)const;
};


class LexicographicTermOrder : public TermOrder
{
  int largest;
 public:
  LexicographicTermOrder(int largest=0);
  int rowDot(int row, const IntegerVector &v)const;
  bool operator()(const IntegerVector &a, const IntegerVector &b, int scaleA=1, int scaleB=1, int perturbationDegree=-1)const;
  //  bool operator()(const IntegerVector &a, const IntegerVector &b)const;
  void print(Printer &p)const;
};


class ReverseLexicographicTermOrder : public TermOrder
{
  int largest;
  int index(int row, const IntegerVector &a)const;
 public:
  ReverseLexicographicTermOrder(int largest=0);
  virtual int rowDot(int row, const IntegerVector &v)const;
  virtual bool operator()(const IntegerVector &a, const IntegerVector &b, int scaleA=1, int scaleB=1, int perturbationDegree=-1)const;
  //  bool operator()(const IntegerVector &a, const IntegerVector &b)const;
  void print(Printer &p)const;
};

class StandardGradedLexicographicTermOrder : public TermOrder
{
  int largest;
 public:
  StandardGradedLexicographicTermOrder(int largest=0);
  //StandardGradedLexicographicTermOrder(const IntegerVector &order);
  virtual int rowDot(int row, const IntegerVector &v)const;
  virtual bool operator()(const IntegerVector &a, const IntegerVector &b, int scaleA=1, int scaleB=1, int perturbationDegree=-1)const;
  //  bool operator()(const IntegerVector &a, const IntegerVector &b)const;
  void print(Printer &p)const;
};


/*class StandardGradedReverseLexicographicTermOrder : public TermOrder
{
  int largest;
 public:
  StandardGradedReverseLexicographicTermOrder(int largest=0);
  virtual int rowDot(int row, const IntegerVector &v)const;
  virtual bool operator()(const IntegerVector &a, const IntegerVector &b, int scaleA=1, int scaleB=1, int perturbationDegree=-1)const;
  void print(Printer &p)const;
};
*/

class WeightTermOrder : public TermOrder  // tie broken with lexicographic
{
  IntegerVector weight;
  //  TermOrder &tieBreaker;  // how do I implement tie breaking. pointer or ref... what happens when you copy a term order?
 public:
  WeightTermOrder(const IntegerVector weight_):weight(weight_){}
  // WeightTermOrder(const IntegerVector weight_, TermOrder ):weight(weight_),tieBreaker(0){}
  virtual int rowDot(int row, const IntegerVector &v)const;
  virtual bool operator()(const IntegerVector &a, const IntegerVector &b, int scaleA=1, int scaleB=1, int perturbationDegree=-1)const;
  //  bool operator()(const IntegerVector &a, const IntegerVector &b)const;
  void print(Printer &p)const;
};


class WeightReverseLexicographicTermOrder : public TermOrder  // tie broken with reverse lexicographic
{
  IntegerVector weight;
 public:
  IntegerVector getWeight()const;
  WeightReverseLexicographicTermOrder(const IntegerVector weight_):weight(weight_){}
  virtual int rowDot(int row, const IntegerVector &v)const;
  virtual bool operator()(const IntegerVector &a, const IntegerVector &b, int scaleA=1, int scaleB=1, int perturbationDegree=-1)const;
  //  bool operator()(const IntegerVector &a, const IntegerVector &b)const;
  void print(Printer &p)const;
};


class MatrixTermOrder : public TermOrder  // tie broken with reverse lexicographic
{
  IntegerVectorList weights;
 public:
  MatrixTermOrder(const IntegerVectorList weights_):weights(weights_){}
  virtual int rowDot(int row, const IntegerVector &v)const;
  virtual bool operator()(const IntegerVector &a, const IntegerVector &b, int scaleA=1, int scaleB=1, int perturbationDegree=-1)const;
  //  bool operator()(const IntegerVector &a, const IntegerVector &b)const;
  void print(Printer &p)const;
};


/*class EliminationTermOrder
{
 public:
   };
*/

#endif
