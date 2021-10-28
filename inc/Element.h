//acknowledgement: cs294-73 colella@eecs.berkeley.edu
#ifndef _ELEMENT_H_
#define _ELEMENT_H_
#define VERTICES 3

class Element
{
public: 
  Element();
  Element(int a_vertices[VERTICES]);
  const int& operator[](const int& a_localNodeNumber) const;
  void vertices(int a_vertices[VERTICES]) const;

private:
  int m_vertices[VERTICES];
};
#endif
