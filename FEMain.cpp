//acknowledgement: cs294-73 colella@eecs.berkeley.edu
#include "FEGrid.h"

#include <iostream>


using namespace std;

int main(int argc, char** argv)
{
  if(argc != 2)
    {
      cout << "this program takes one argument that is the common name prefix of .node and .elem files ";
      cout << "file prefix"<<endl;
      return 1;
    }
  string prefix(argv[1]);
  string nodeFile = prefix+".node";
  string eleFile  = prefix+".elem";

  FEGrid grid(nodeFile, eleFile);
  //insert code here to form the global stiffness matrix and the N.C. vector. Also, compute temperature at each node of the grid based on the stiffness matrix and NC vector. 
  

  
  return 0;
  
}
