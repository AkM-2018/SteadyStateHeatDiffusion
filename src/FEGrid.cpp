//acknowledgement: cs294-73 colella@eecs.berkeley.edu
#include <cstdio>
#include <iostream>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
#include <fstream>
#include "Node.h"
#include "Element.h"
#include "FEGrid.h"

FEGrid::FEGrid() : m_numInteriorNodes(0)
{
}

FEGrid::FEGrid(const std::string &a_nodeFileName, const std::string &a_elementFileName)
{
  ifstream nodes(a_nodeFileName.c_str());
  int ncount;
  nodes >> ncount;
  length = 0.6;
  width = 0.4;
  m_nodes.resize(ncount);
  m_numInteriorNodes = 0;
  for (int i = 0; i < ncount; i++)
  {
    int vertex;
    float x[DIM];
    nodes >> vertex >> x[0] >> x[1];
    vertex--;
    if (x[0] == length)
    {
      m_nodes[vertex] = Node(x, -1, false);
    }
    else if (x[1] == width)
    {
      m_nodes[vertex] = Node(x, -1, false);
    }
    else
    {
      m_nodes[vertex] = Node(x, m_numInteriorNodes, true);
      m_numInteriorNodes++;
    }
  }

  ifstream elements(a_elementFileName.c_str());
  int ncell;
  elements >> ncell;
  int vert[VERTICES];
  m_elements.resize(ncell);
  for (int i = 0; i < ncell; i++)
  {
    int cellID;
    elements >> cellID >> vert[0] >> vert[1] >> vert[2];
    vert[0]--;
    vert[1]--;
    vert[2]--;
    cellID--;
    m_elements[cellID] = Element(vert);
  }
}

//returns the gradient of shape functions N1 N2 N3.
void FEGrid::gradient(float a_gradient[DIM],
                      const int &a_eltNumber,
                      const int &a_nodeNumber) const
{
  const Element &e = m_elements[a_eltNumber];
  const Node &n = m_nodes[e[a_nodeNumber]];
  assert(n.isInterior()); //make sure that you are computing it for interior node only.
  float xbase[DIM];
  n.getPosition(xbase);
  float dx[VERTICES - 1][DIM];
  for (int ivert = 0; ivert < VERTICES - 1; ivert++)
  {
    int otherNodeNumber = e[(a_nodeNumber + ivert + 1) % VERTICES];
    m_nodes[otherNodeNumber].getPosition(dx[ivert]);
    for (int idir = 0; idir < DIM; idir++)
    {
      dx[ivert][idir] -= xbase[idir];
    }
  }
  // for(int i=0;i<2;i++)
  // {
  //   for(int j=0;j<DIM;j++)
  //   {
  //     std::cout<<dx[i][j]<<" ";
  //   }
  //   std::cout<<endl;
  // }
  // WARNING: the following calculation is correct for triangles in 2D *only*.
  float det = dx[0][0] * dx[1][1] - dx[1][0] * dx[0][1];

  a_gradient[0] = (-(dx[1][1] - dx[0][1]) / det);
  a_gradient[1] = (-(dx[1][0] - dx[0][0]) / det);
}

float FEGrid::elementArea(const int &a_eltNumber) const
{
  const Element &e = m_elements[a_eltNumber];
  const Node &n = m_nodes[e[0]];
  float xbase[DIM];
  n.getPosition(xbase);
  float dx[VERTICES - 1][DIM];
  for (int ivert = 1; ivert < VERTICES; ivert++)
  {
    int otherNodeNumber = e[ivert];
    m_nodes[otherNodeNumber].getPosition(dx[ivert - 1]);
    for (int idir = 0; idir < DIM; idir++)
    {
      dx[ivert - 1][idir] -= xbase[idir];
    }
  }

  // WARNING: the following calculation is correct for triangles in 2D *only*.
  float area = fabs(dx[0][0] * dx[1][1] - dx[1][0] * dx[0][1]) / 2;
  return area;
}

const Node &FEGrid::getNode(const int &a_eltNumber, const int &a_localNodeNumber) const
{
  return m_nodes[m_elements[a_eltNumber][a_localNodeNumber]];
}

int FEGrid::getNumElts() const
{
  return m_elements.size();
}

int FEGrid::getNumNodes() const
{
  return m_nodes.size();
}

int FEGrid::getNumInteriorNodes() const
{
  return m_numInteriorNodes;
}

const Element &FEGrid::element(int i) const
{
  return m_elements[i];
}
const Node &FEGrid::node(int i) const
{
  return m_nodes[i];
}
