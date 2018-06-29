#include <stdio.h>
#include <stdlib.h>

#include "graph.h"
#include "sorting.h"
#include "util.h"

double inc_deg_key(struct Graph *g, int v) { return g->degree[v]; }
double dec_deg_key(struct Graph *g, int v) { return -g->degree[v]; }
void order_vertices(int *vv, struct Graph *g, int vtx_ordering) {
    for (int i=0; i<g->n; i++)
        vv[i] = i;

    switch(vtx_ordering) {
    case  0: break;  // no sorting
    case  1: INSERTION_SORT_VV(inc_deg_key) break;
    case -1: INSERTION_SORT_VV(dec_deg_key) break;
    case 10: INSERTION_SORT_VV(dec_deg_key) break;
    default: fail("Unrecognised vertex order");
    }
}


