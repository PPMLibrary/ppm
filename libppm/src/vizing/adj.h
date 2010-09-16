#ifndef ADJ_H
#define ADJ_H
#include "list.h"
#include "Array.h"
#include "edge.h"
#include <iostream>
using namespace std;

class Adj {

    friend  ostream& operator<<(ostream &output, const Adj &graph);
    friend class Graph;
    friend class ColorGraph;

public:
    Adj(int);
    Adj(int, int); 
    Adj(int v, Array<Edge> A);
    //copy constructor;
    Adj(const Adj &);
    ~Adj();// { delete [] adj; } 
    int getVertices() { return vertices; }
    int getEdges() { return edges; }
    const Adj & operator=(const Adj &);
    List getList(int u); //return the list of edges joining u;
    int maxDegree();
    List **getadj() { return adj; }
    
private:
    typedef List* PtrToList;
    int vertices;    //# of vertices;
    int edges;       //# of edges;
    PtrToList *adj;  //adjacent list; 
};

/*---------------------implementation-------------------------*/

Adj::Adj(int v)
{
    vertices = v;
    edges = 0;
  
    adj = new PtrToList[vertices];
    for(int i=0; i<vertices; i++)
    adj[i] = new List;
}

Adj::~Adj()
{
    int i;
    for(i=0; i<vertices; i++)
        delete adj[i];
  
    delete [] adj;
}

Adj::Adj(int v, int e)
{
    int i;
    vertices = v;
    edges = e;

    adj = new PtrToList[vertices];

    for(i=0; i<vertices; i++)
        adj[i] = new List;


    for(i=0; i<e; i++) {
        int x, y;
      
        cin >> x >> y;
        assert( x < y );
        adj[x-1]->Prepend(y);
        adj[y-1]->Prepend(x);
    }
}

Adj::Adj(int v, Array<Edge> A)
{
    int i;
    vertices = v;
    edges = A.getSize();
    adj = new PtrToList[vertices];

    for(i=0; i<vertices; i++)
        adj[i] = new List;

    for(i=0; i<edges; i++) {
        int x, y;
        x = A[i].startVx;
        y = A[i].endVx;
      
        adj[x-1]->Append(y);
        adj[y-1]->Append(x);
    }
}

Adj::Adj(const Adj &old)
{
    adj = 0;
    *this = old;
}

const Adj&
Adj::operator=(const Adj& old)
{
    if (&old == this) {}
        //Self assignment; do nothing;
    else {
        delete [] adj;
        vertices = old.vertices;
        edges = old.edges;
        adj = new PtrToList[vertices];
        for (int i = 0; i < vertices; i++)  
		    adj[i] = old.adj[i];
    }
    return *this;
}
                                                       
List
Adj::getList(int u)
{
    List L;
    L = *(adj[u-1]);
    return L;
}

//return maximum of degrees at each vertex
int
Adj::maxDegree()
{
    int max = 0;

    for (int i = 0; i < vertices; i++)
        if ( max < (adj[i] -> length()))
	        max = adj[i] -> length();
    
    return max;
}
    
ostream& operator<<(ostream &output, const Adj &graph)
{
    output << "Number of Vertices  " << graph.vertices << endl;
    output << "Number of Edges  " << graph.edges << endl;
    output << "Edges of The Graph:"<< endl;

    for (int i=0; i<graph.vertices; i++) {
	    ListIterator I;
	    for (I.start(*(graph.adj[i])); !I.done(); ++I)
	        if ( (i+1) < I())
			    output << i+1 <<' '<< I() << endl; 
    }
    
    return output;
}

#endif
