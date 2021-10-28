//acknowledgement: cs294-73 colella@eecs.berkeley.edu
#ifndef _NODE_H_
#define _NODE_H_
#define DIM 2
#include <cstdio>
#include <cmath>
#include <cassert>
#include <vector>
using namespace std;
class Node
{
 public:
  Node();
  Node(float a_position[DIM],
       const int&  a_interiorNodeID, 
       const bool& a_isInterior);
  void  getPosition(float a_position[DIM]) const;
  const int& getInteriorNodeID() const;
  const bool& isInterior() const;

 private:
  float m_position[DIM];
  bool  m_isInterior;
  int m_interiorNodeID;
};
#endif
