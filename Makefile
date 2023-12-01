#acknowledgement: cs294-73 colella@eecs.berkeley.edu
CFLAGS = -g -Wall
#CFLAGS = -O3
CXX = g++

all: pa4
 
pa4: FEMain.o FEGrid.o Element.o Node.o 
	  $(CXX) -I inc/ FEMain.o FEGrid.o Element.o Node.o 

FEMain.o:
	$(CXX) -I inc/ $(CFLAGS) -c -o FEMain.o src/FEMain.cpp

Node.o:
	$(CXX) -I inc/ $(CFLAGS) -c -o Node.o src/Node.cpp

Element.o:
	$(CXX) -I inc/ $(CFLAGS) -c -o Element.o src/Element.cpp

FEGrid.o:
	$(CXX) -I inc/ $(CFLAGS) -c -o FEGrid.o src/FEGrid.cpp

team:
	@echo "Akhil Manoj"
	@echo "Shivam Chaturvedi"

clean:
	rm -f *.o a.out