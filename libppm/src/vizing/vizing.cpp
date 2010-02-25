// This is a wrapper to the C++ code of approximate edge coloring from:
// The Stony Brook Algorithm Repository 
// http://www.cs.sunysb.edu/~algorith/implement/stony/distrib/Vizing/

#include "color.h"
#include "graph.h"
using namespace std;

extern "C" {
#ifdef __XLF
  int vizing_coloring(int *nvertices, int *nedges, int *edges, int *coloring);
#else
  int vizing_coloring_(int *nvertices, int *nedges, int *edges, int *coloring);
#endif
}

// when receiving data from fortran, everything is a pointer.
// also pay attention to the row/column major in arrays
#ifdef __XLF
int vizing_coloring(int *nvertices, int *nedges, int *edges, int *coloring)
#else
int vizing_coloring_(int *nvertices, int *nedges, int *edges, int *coloring)
#endif
{
    int v, e, i;
    v = *nvertices;
    e = *nedges;
    Array<Edge> A(e);

    for (i = 0; i < e; i++) {
        Edge tmp(edges[2*i], edges[2*i+1]);
        A.insert(i, tmp);
    }

    // do the magic ...
    ColorGraph G(v, A);
    // coloring receives the results in the form: (p1,p2,c),(),(),...
    G.getResult(coloring);

    return 0;
}
