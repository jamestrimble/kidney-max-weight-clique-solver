#define _GNU_SOURCE

#include "graph.h"
#include "util.h"
#include "bitset.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

struct Weight default_weight()
{
    return (struct Weight) {{0, 0, 0, 0, 1}};
}

void print_weight(struct Weight const wt)
{
    printf("%ld", wt.weight[0]);
    for (int i=1; i<WEIGHT_SIZE; i++) {
        printf(",%ld", wt.weight[i]);
    }
}

void add_edge(struct Graph *g, int v, int w) {
    unset_bit(g->bit_complement_nd[v], w);
    unset_bit(g->bit_complement_nd[w], v);
}

void calculate_all_degrees(struct Graph *g) {
    for (int v=0; v<g->n; v++) {
        g->degree[v] = g->n - 1 - bitset_popcount(g->bit_complement_nd[v], (g->n+BITS_PER_WORD-1)/BITS_PER_WORD);
    }
}

// Checks if a set of vertices induces a clique
bool check_clique(struct Graph* g, struct VtxList* clq) {
    struct Weight total_wt = {};
    for (int i=0; i<clq->size; i++)
        total_wt = weight_sum(total_wt, g->weight[clq->vv[i]]);
    if (weight_eq(total_wt, clq->total_wt))
        return true;

    for (int i=0; i<clq->size-1; i++)
        for (int j=i+1; j<clq->size; j++)
            if (test_bit(g->bit_complement_nd[clq->vv[i]], clq->vv[j]))
                return false;
    return true;
}

struct Graph *new_graph(int n)
{
    struct Graph *g = calloc(1, sizeof(*g));
    g->n = n;
    g->degree = calloc(n, sizeof(*g->degree));
    g->weight = calloc(n, sizeof(*g->weight));
    g->bit_complement_nd = calloc(n, sizeof(*g->bit_complement_nd));
    for (int i=0; i<n; i++) {
        g->bit_complement_nd[i] = calloc((n+BITS_PER_WORD-1)/BITS_PER_WORD, sizeof *g->bit_complement_nd[i]);
    }

    // initialise bit_complement_nd
    unsigned long long *bit_complement_nd_template = calloc((n+BITS_PER_WORD-1)/BITS_PER_WORD, sizeof *bit_complement_nd_template);
    int word = 0;
    int remaining = n;
    while (remaining >= 64) {
        bit_complement_nd_template[word++] = ~0ull;
        remaining -= 64;
    }
    for (int i=0; i<remaining; i++) {
        bit_complement_nd_template[word] |= (1ull << i);
    }
    for (int i=0; i<n; i++) {
        for (int j=0; j<(n+BITS_PER_WORD-1)/BITS_PER_WORD; j++) {
            g->bit_complement_nd[i][j] = bit_complement_nd_template[j];
        }
        unset_bit(g->bit_complement_nd[i], i);
    }
    free(bit_complement_nd_template);

    return g;
}

void free_graph(struct Graph *g)
{
    for (int i=0; i<g->n; i++) {
        free(g->bit_complement_nd[i]);
    }
    free(g->degree);
    free(g->weight);
    free(g->bit_complement_nd);
    free(g);
}

struct Graph *induced_subgraph(struct Graph *g, int *vv, int vv_len) {
    struct Graph* subg = new_graph(vv_len);
    for (int i=0; i<subg->n; i++)
        for (int j=0; j<subg->n; j++)
            if (!test_bit(g->bit_complement_nd[vv[i]], vv[j]))
                unset_bit(subg->bit_complement_nd[i], j);

    for (int i=0; i<subg->n; i++)
        subg->weight[i] = g->weight[vv[i]];
    return subg;
}

struct Graph *readGraph(char* filename) {
    FILE* f;
    
    if ((f=fopen(filename, "r"))==NULL)
        fail("Cannot open file");

    char* line = NULL;
    size_t nchar = 0;

    int nvertices = 0;
    int medges = 0;
    int v, w;
    int edges_read = 0;
    long wt;

    struct Graph *g = NULL;

    struct Weight weight;
    int offset;
    int read_char_count;

    while (getline(&line, &nchar, f) != -1) {
        if (nchar > 0) {
            switch (line[0]) {
            case 'p':
                if (sscanf(line, "p edge %d %d", &nvertices, &medges)!=2)
                    fail("Error reading a line beginning with p.\n");
                printf("%d vertices\n", nvertices);
                printf("%d edges\n", medges);
                g = new_graph(nvertices);
                for (int i=0; i<nvertices; i++)
                    g->weight[i] = default_weight();
                break;
            case 'e':
                if (sscanf(line, "e %d %d", &v, &w)!=2)
                    fail("Error reading a line beginning with e.\n");
                add_edge(g, v-1, w-1);
                edges_read++;
                break;
            case 'n':
                if (sscanf(line, "n %d%n", &v, &read_char_count)!=1)
                    fail("Error reading a line beginning with n.\n");
                offset = read_char_count;
                for (int i=0; i<WEIGHT_SIZE; i++) {
                    if (sscanf(line + offset, "%ld%n", &wt, &read_char_count)!=1)
                        fail("Error reading a line beginning with n.\n");
                    offset += read_char_count;
                    weight.weight[i] = wt;
                }
                g->weight[v-1] = weight;
                break;
            }
        }
    }

    if (medges>0 && edges_read != medges) fail("Unexpected number of edges.");

    fclose(f);
    return g;
}

void init_VtxList(struct VtxList *l, int capacity)
{
    l->vv = malloc(capacity * sizeof *l->vv);
    l->size = 0;
    l->total_wt = (struct Weight) {};
}

void destroy_VtxList(struct VtxList *l)
{
    free(l->vv);
}

void init_UnweightedVtxList(struct UnweightedVtxList *l, int capacity)
{
    l->vv = malloc(capacity * sizeof *l->vv);
    l->size = 0;
}

void destroy_UnweightedVtxList(struct UnweightedVtxList *l)
{
    free(l->vv);
}

void vtxlist_push_vtx(struct Graph *g, struct VtxList *L, int v)
{
    L->vv[L->size++] = v;
    L->total_wt = weight_sum(L->total_wt, g->weight[v]);
}

void vtxlist_pop_vtx(struct Graph *g, struct VtxList *L)
{
    L->size--;
    L->total_wt = weight_difference(L->total_wt, g->weight[L->vv[L->size]]);
}

void copy_VtxList(struct VtxList *src, struct VtxList *dest)
{
    dest->size = src->size;
    dest->total_wt = src->total_wt;
    for (int i=0; i<src->size; i++)
        dest->vv[i] = src->vv[i];
}
