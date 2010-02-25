#ifndef GRAPH_H
#define GRAPH_H
#include "Boolean.h"
#include "adj.h"
#include "edge.h"
using namespace std;

class Graph {

    friend class ColorGraph;
    friend ostream& operator<<(ostream &output, const Graph &graph);
    enum NodeColor {white, gray, black};
    
public:
    Graph(int v);
    Graph(int v, int e);
    Graph(int v, Array<Edge> A);
    const Graph & operator=(const Graph &);
    ~Graph() { delete theAdj; delete color; }
    Graph(const Graph &old);
    int getVertices() { return theAdj->getVertices(); }
    int getEdges() { return theAdj->getEdges(); }
    int maxDegree() { return theAdj->maxDegree(); }
    Boolean isConnected(int, int);
    void BFS(int);
    void whiten();
    void AddOneEdge(Edge E);
    Adj *gettheAdj() { return theAdj; }
   
private:
    Adj *theAdj;
    NodeColor *color;

};

/*---------------------implementation-------------------------*/

Graph::Graph(int v)
{
    theAdj = new Adj(v);
    color = new NodeColor[v];

    for ( int i = 0; i < v; i++ )
        color[i] = white;
}
  
Graph::Graph(int v, int e) 
{
    theAdj = new Adj(v, e);
    color = new NodeColor[v];

    for ( int i = 0; i < v; i++ )
        color[i] = white;
}

Graph::Graph(int v, Array<Edge> A)
{
    theAdj = new Adj(v, A);
    color = new NodeColor[v];

    for (int i = 0; i < v; i++)
        color[i] = white;
}

const Graph & 
Graph::operator=(const Graph & old) 
{
    if (&old == this) {}
    //Self assignment, do nothing;
    else {
        delete theAdj;
        int v = old.theAdj->vertices;
        theAdj = new Adj(v);
        *theAdj = *old.theAdj;
   
        delete [] color;
        color = new NodeColor[v];

        for (int i=0; i<v; i++)
            color[i] = old.color[i];
    }
    return *this;
}

Graph:: Graph(const Graph &old)
{
    theAdj = 0;
    color = 0;
    *this = old;
}

void
Graph::BFS(int v)
{
    List Q, L;

    color[v-1] = gray;
    Q.Append(v);

    while (!Q.isEmpty()) {
        int u = Q.nth(1);
        L = theAdj -> getList(u);
        ListIterator I;
        for (I.start(L); !I.done(); ++I)
            if (color[I()-1] == white) {
	            color[I()-1] = gray;
	            Q.Append(I());
            }
        Q.deleteNth(1);
        color[u-1] = black;
    }
} 

void
Graph::whiten()
{
    int v = theAdj->vertices;
    for (int i = 0; i < v; i++)
	    color[i] = white;
}

Boolean
Graph::isConnected(int v, int w)
{
    whiten();
    BFS(v);
    if (color[w-1] == black)
        return TRUE;
    else
        return FALSE;
}

void 
Graph::AddOneEdge(Edge)
{
}

ostream& operator<<(ostream &output, const Graph &graph)
{
    output << *graph.theAdj;
    return output;
}

#endif
