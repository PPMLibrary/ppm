#ifndef COLOR_H
#define COLOR_H
#include "graph.h"
#include <assert.h>
using namespace std;

typedef Color* ColorPtr;
  
//if color = 0, the edge is uncolored (default value);
//color starts from 1 to max.

class ColorGraph {
    friend  ostream& operator<<(ostream &output, const ColorGraph &graph);

public:
    ColorGraph(int, Array<Edge>);
    ColorGraph(int v, int e);
    
    Boolean isColor(int v, Color s);
    Color missingColor(int vertex);
    Edge locateEdge(int v, Color s);
    Boolean isIn(int edVx, Array<int> endVerices, int range, int& position);
    void AddOneEdge(Array<Edge> A, int k);
    void recolor(int, Array<int>, Array<Color>, int,  int);
    Array<Edge> searchEdges(Array<Edge> A, Color s, Color t);
    Graph* subGraph(Array<Edge> A, Color s, Color t);
    void switchColor(Graph *subPtr, Color s, Color t, int vertex);
    void getResult(int* res);

private:
    Graph *G;
    //record the colors of each edge, for efficient search;
    ColorPtr *colorMatrix;

    int max; //maximum degree+1 (maximum colors allowed);
    int count; //count # of edges already colored;
};

/*---------------------implementation-------------------------*/

ColorGraph::ColorGraph(int v, Array<Edge> A)
{
    G = new Graph(v, A);
    int size = A.getSize();
    count = 0;
    max = G -> maxDegree() +1;

    //color initialization;

    colorMatrix = new ColorPtr[v];
    for (int i = 0; i < v; i++) {
        colorMatrix[i] = new Color[v];
        for (int j = 0; j < v; j++)
	        colorMatrix[i][j] = 0;
    }

    //start coloring processing
    while ( count < size ) {
        if ( count < max ) {
            int start = A[count].startVx;
            int end = A[count].endVx;
            colorMatrix[start-1][end-1] = ++count;
            colorMatrix[end-1][start-1] = count;
        } else {
        //coloring one more edge;
            AddOneEdge(A, count++);
        }
    }
}

ColorGraph::ColorGraph(int v, int e)
{
    G = new Graph(v);
    colorMatrix = new ColorPtr[v];

    for (int i = 0; i < v; i++) {
        colorMatrix[i] = new Color[v];
        for (int j = 0; j < v; j++)
            colorMatrix[i][j] = 0;
    }
    count = 0;
}

Color
ColorGraph::missingColor(int vertex)
{
    Color tmp;
    for (Color color = 1; color <= max; color++) 
        if (!isColor(vertex, color)) {
            tmp = color;
            color = max+1;
    }
    return tmp;
}

//Is "s" a color of some edge at v
Boolean
ColorGraph::isColor(int v, Color s)
{
    int size = G->getVertices();
   
    for (int i=0; i<size; i++) {
        assert( colorMatrix[v-1][i] == colorMatrix[i][v-1]);
        if ( colorMatrix[v-1][i] == s)
	        return TRUE;
    }
    return FALSE;
 }

//search for an edge with vertex v and color s
//assume that such an edge exists
Edge
ColorGraph::locateEdge(int v, Color s)
{
    Edge E;
    assert (isColor(v, s));
    int vertices = G->getVertices();
    for (int i = 0; i < vertices; i++) {
        if (colorMatrix[v-1][i] == s) {
            E.startVx = v;
            E.endVx = i+1;
            E.color = s;
            i=vertices;
        }
	}	
   return E;
}

Boolean 
ColorGraph::isIn(int edVx, Array<int> endVertices, int range, int& position) 
{
    for (int i = 0; i < range; i++) {
        if (edVx == endVertices[i]) {
		    position = i;
            return TRUE;
		}
    }
    return FALSE;
}    

void
ColorGraph::AddOneEdge(Array<Edge> A, int k)
{
    Edge E = A[k];
    int x = E.startVx;
    Array<Edge> xEdges(max-1);
    Array<int> endVertices(max-1);
  
    int v = G -> getVertices();
    Array<Color> missingCol(v);
    int xEdgeCount = 0;
    Color t;
    Color s = missingColor(x);
    missingCol.insert(x-1, s);
    int tmp = E.endVx;
    Boolean flag = TRUE;
    int position;
  
	while (!isIn(tmp, endVertices, xEdgeCount, position) && flag) {
        endVertices.insert(xEdgeCount, tmp);
        xEdges.insert(xEdgeCount, E);
        xEdgeCount++;
        t = missingColor(tmp);
        missingCol.insert(tmp-1, t);
    
        if (isColor(x, t)) {
		    E = locateEdge(x, t);
            tmp = E.endVx;
        } else {
            flag = FALSE;
        }
    }

    if (flag) {
        //recolor the edges for i < position;
        int y = endVertices[position];
        assert(y==tmp);
        int z = endVertices[xEdgeCount-1];

        recolor(x,  endVertices, missingCol, 0, position);

        Graph* subPtr=subGraph(A, s, t);

        if (!(subPtr->isConnected(x, y))) {
            switchColor(subPtr, s, t, y);
            colorMatrix[x-1][y-1] = s;
            colorMatrix[y-1][x-1] = s;
        } else {
            assert (!(subPtr->isConnected(x, z)));
            recolor(x, endVertices, missingCol, position, xEdgeCount-1);
            switchColor(subPtr, s, t, z);
            colorMatrix[x-1][z-1] = s;
            colorMatrix[z-1][x-1] = s;
		}
        delete subPtr;
    } else {
        for (int i = xEdgeCount-2; i >= 0; i--) {
            int j = endVertices[i]-1;
			colorMatrix[x-1][j] = missingCol[j];
            colorMatrix[j][x-1] = missingCol[j];
        }
        colorMatrix[x-1][tmp-1] = t;
        colorMatrix[tmp-1][x-1] = t;
    }
}

void
ColorGraph::recolor(int x, Array<int> endVertices, Array<Color> missingCol,
           int start,  int position) 
{
    int end;
    assert (start <= position);
    for (int i = start; i < position; i++) {
        end = endVertices[i];
        colorMatrix[x-1][end-1] = missingCol[end-1];
        colorMatrix[end-1][x-1] = missingCol[end-1];
    }
    end = endVertices[position];
    colorMatrix[x-1][end-1] = 0;
    colorMatrix[end-1][x-1] = 0;
}

Array<Edge>
ColorGraph::searchEdges(Array<Edge> A, Color s, Color t)
{
    Array<Edge> first(count);
    int i;
    int tmp = 0;
    for (i = 0; i < count; i++) {
        int start = A[i].startVx;
        int end = A[i].endVx;

        if ( colorMatrix[start-1][end-1] == s || colorMatrix[start-1][end-1] == t){
            Edge choose(start, end, colorMatrix[start-1][end-1]);
            first.insert(tmp++, choose);
        }
    }

    Array<Edge> final(tmp);
    for (i = 0; i < tmp; i++)
        final.insert(i, first[i]);
    
	return final;
}
     
Graph* 
ColorGraph::subGraph(Array<Edge> A, Color s, Color t)
{
     int v = G->getVertices();
     Graph *graphPtr = new Graph(v, searchEdges(A, s, t));
     return graphPtr;
}

void
ColorGraph::switchColor(Graph *subPtr, Color s, Color t, int vertex)
{
    int v = subPtr -> getVertices();

    for (int i = 0; i < v; i++) {
        List L = *((subPtr -> theAdj)->adj[i]);
        ListIterator I;
        if ( subPtr -> isConnected(i+1, vertex)) {
            for (I.start(L); !I.done(); ++I) {
                if ( (i+1) < I()) {
	                if (colorMatrix[i][I()-1] == s) {
					     colorMatrix[i][I()-1] = t;
		                 colorMatrix[I()-1][i] = t;
	                } else {
		                 assert ( colorMatrix[i][I()-1] == t );
		                 colorMatrix[i][I()-1] = s;
		                 colorMatrix[I()-1][i] = s;
					}
                }
            }
        }
    }
}

void ColorGraph::getResult(int *res) 
{
    // assume res is allocated to proper size 3*e !!!
    int v = G->getVertices();
    int e = G->getEdges();
    int k = 0;

    for (int i = 0; i < v; i++) {
         ListIterator I;
         List L;
         L = *(((G->gettheAdj())->getadj())[i]);
         for (I.start(L); !I.done(); ++I)
             if ( (i+1) < I()) {
		         res[k] = i+1;
	             res[k+1] = I();
		         res[k+2] = colorMatrix[i][I()-1]; 
		         k+=3;
             }
    }
    delete colorMatrix;
    delete G;
    return;
}

ostream& operator<<(ostream &output, const ColorGraph &graph)
{
    output << "Edges  " << "  Color " << endl;
    int v = graph.G->getVertices();
    for (int i = 0; i < v; i++) {
         ListIterator I;
         List L;
         L = *(((graph.G->gettheAdj())->getadj())[i]);
         for (I.start(L); !I.done(); ++I)
              if ( (i+1) < I())
                  output << i+1 <<' '<< I() <<"       "
		          << graph.colorMatrix[i][I()-1]<< endl; 
    }
    return output;
}

#endif
