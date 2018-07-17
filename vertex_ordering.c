#include <stdio.h>
#include <stdlib.h>

#include "bitset.h"
#include "graph.h"
#include "sorting.h"
#include "util.h"

double inc_deg_key(struct Graph *g, int v) { return g->degree[v]; }
double dec_deg_key(struct Graph *g, int v) { return -g->degree[v]; }

void calc_weighted_degs(struct Graph *g, struct Weight *residual_weighted_deg) {
    for (int i=0; i<g->n; i++) {
        residual_weighted_deg[i] = (struct Weight) {};
        for (int j=0; j<g->n; j++) {
            if (!test_bit(g->bit_complement_nd[i], j) && i != j) {
                residual_weighted_deg[i] = weight_sum(residual_weighted_deg[i], g->weight[j]);
            }
        }
    }
}

void carraghan_pardalos_order(int *vv, struct Graph *g, bool reverse) {
    struct Weight *residual_weighted_deg = malloc(g->n * sizeof(*residual_weighted_deg));

    calc_weighted_degs(g, residual_weighted_deg);

    for (int i=0; i<g->n; i++) {
        // find vertex with lowest residual_weighted_deg
        int best_v_pos = i;
        struct Weight best_wt_deg = residual_weighted_deg[vv[i]];
        for (int j=i+1; j<g->n; j++) {
            int v = vv[j];
            if (weight_lt(residual_weighted_deg[v], best_wt_deg)) {
                best_wt_deg = residual_weighted_deg[v];
                best_v_pos = j;
            }
        }
        int v = vv[best_v_pos];
        vv[best_v_pos] = vv[i];
        vv[i] = v;

        for (int j=i+1; j<g->n; j++) {
            int w = vv[j];
            if (!test_bit(g->bit_complement_nd[v], w) && v != w)
                residual_weighted_deg[w] = weight_difference(residual_weighted_deg[w], g->weight[v]);
        }
    }

    if (reverse) {
        for (int i=0; i<g->n/2; i++) {
            int tmp = vv[i];
            vv[i] = vv[g->n - 1 - i];
            vv[g->n - 1 - i] = tmp;

        }
    }
    free(residual_weighted_deg);
}

void order_vertices(int *vv, struct Graph *g, int vtx_ordering) {
    for (int i=0; i<g->n; i++)
        vv[i] = i;

    switch(vtx_ordering) {
    case  0: break;  // no sorting
    case  1: INSERTION_SORT_VV(inc_deg_key) break;
    case -1: INSERTION_SORT_VV(dec_deg_key) break;
    case  5: carraghan_pardalos_order(vv, g, false); break;
    case -5: carraghan_pardalos_order(vv, g, true); break;
    case 10: INSERTION_SORT_VV(dec_deg_key) break;
    default: fail("Unrecognised vertex order");
    }
}


