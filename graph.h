#ifndef GRAPH_H_INCLUDED
#define GRAPH_H_INCLUDED

#include <set>
#include <string>

using namespace std;

class Graph{
  int n;
  bool isDirected;
  class Edge{
  public:
    int a,b;
  Edge(int a_, int b_):
    a(a_),
      b(b_)
      {
      }
    bool operator<(const Edge &e)const
    {
      return (a<e.a)||((a==e.a)&&(b<e.b));
    }
  };
  set<Edge> edges;
 public:
 Graph(int numberOfVertices, bool isDirected_=false):
  n(numberOfVertices),
    isDirected(isDirected_)
  {
  }
  void addEdge(int a, int b)
    {
      edges.insert(Edge(a,b));
    }
  int diameter(int *timesAttained=0)const;
  string toString()const;
};
#endif
