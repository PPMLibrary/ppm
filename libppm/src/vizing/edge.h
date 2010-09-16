#ifndef EDGE_H
#define EDGE_H
#include <iostream>
using namespace std;

typedef int Color;

class Edge {

    friend  ostream& operator<<(ostream &output, const Edge &E);

public:

    static const Edge zero;
    
    Edge(){startVx=0; endVx=0; color=0;}
    Edge(int);
    Edge(int, int);
    Edge(int, int, int);
    const Edge& operator=(const Edge&);
  
    int startVx;
    int endVx;
    Color color;
};

const Edge Edge::zero = 0;  //if color = 0, the edge is uncolored.

/*---------------------implementation-------------------------*/

Edge::Edge(int x)
{
    startVx = x;
    endVx = 0;
    color = 0;
}

Edge::Edge(int x, int y)
{
    startVx  = x;
    endVx = y;
    color = 0;
}

Edge::Edge(int x, int y, int c)
{
    startVx = x;
    endVx = y;
    color = c;
}

const Edge&
Edge::operator=(const Edge& old) 
{
    if (&old == this) {} 
	else {
	    startVx = old.startVx;
	    endVx = old.endVx;
		color = old.color;
    }
    return *this;
}

ostream& operator<<(ostream &output, const Edge &E)
{
    output << E.startVx << "   "<< E.endVx 
	   << "  "<< E.color<< endl;
    return output;
}

#endif
