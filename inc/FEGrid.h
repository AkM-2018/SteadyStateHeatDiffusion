//acknowledgement: cs294-73 colella@eecs.berkeley.edu
#ifndef _FEGRID_H_
#define _FEGRID_H_

#include <cstdio>
#include <cmath>
#include <cassert>
#include <cstring>
#include <vector>
#include <string>
#include "Node.h"
#include "Element.h"

using namespace std;

class FEGrid
{
public:
  FEGrid();
  FEGrid(const std::string &nodeFile, const std::string &a_elementFileName);
  void gradient(float a_gradient[DIM],
                const int &a_eltNumber,
                const int &a_localNodeNumber) const;
  float elementArea(const int &a_eltNumber) const;

  const Node &getNode(const int &a_eltNumber,
                      const int &a_localNodeNumber) const;
  int getNumElts() const;
  int getNumNodes() const;
  int getNumInteriorNodes() const;
  const Element &element(int i) const;
  const Node &node(int i) const;

private:
  vector<Node> m_nodes;
  vector<Element> m_elements;
  int m_numInteriorNodes;
  double length;
  double width;
};

#endif
