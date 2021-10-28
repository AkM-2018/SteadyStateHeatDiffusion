#acknowledgement: cs294-73 colella@eecs.berkeley.edu
CFLAGS = -g -Wall
#CFLAGS = -O3
CXX = g++

all: pa4
 
pa4: FEMain.o FEGrid.o Element.o Node.o 
	  $(CXX) FEMain.o FEGrid.o Element.o Node.o 

FEMain.o: FEMain.cpp FEGrid.H 
	$(CXX) $(CFLAGS) -c -o FEMain.o FEMain.cpp

Node.o: Node.cpp Node.H
	$(CXX) $(CFLAGS) -c -o Node.o Node.cpp

Element.o: Element.H Element.cpp 
	$(CXX) $(CFLAGS) -c -o Element.o Element.cpp

FEGrid.o: FEGrid.H FEGrid.cpp Element.H Node.H $(VW)/VisitWriter.H
	$(CXX) -I$(VW) $(CFLAGS) -c -o FEGrid.o FEGrid.cpp
