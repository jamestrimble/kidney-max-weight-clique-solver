#ifndef GRAPH_H
#define GRAPH_H

#include <limits.h>
#include <stdbool.h>

#define BYTES_PER_WORD sizeof(unsigned long long)
#define BITS_PER_WORD (CHAR_BIT * BYTES_PER_WORD)

struct Weight
{
    long weight;
};

bool weight_lt(struct Weight const wt0, struct Weight const wt1);

bool weight_gt(struct Weight const wt0, struct Weight const wt1);

bool weight_eq(struct Weight const wt0, struct Weight const wt1);

bool weight_gt_zero(struct Weight const wt0);

struct Weight weight_sum(struct Weight const wt0, struct Weight const wt1);

struct Weight weight_difference(struct Weight const wt0, struct Weight const wt1);

struct Weight default_weight();

void print_weight(struct Weight const wt);

struct Graph {
    int n;
    int *degree;
    struct Weight *weight;
    bool **adjmat;
    unsigned long long **bit_complement_nd;
};

struct VtxList {
    struct Weight total_wt;
    int size;
    int *vv;
};

struct UnweightedVtxList {
    int size;
    int *vv;
};

void init_VtxList(struct VtxList *l, int capacity);
void destroy_VtxList(struct VtxList *l);
void init_UnweightedVtxList(struct UnweightedVtxList *l, int capacity);
void destroy_UnweightedVtxList(struct UnweightedVtxList *l);

void add_edge(struct Graph *g, int v, int w);

void calculate_all_degrees(struct Graph *g);

// Checks if a set of vertices induces a clique
bool check_clique(struct Graph* g, struct VtxList* clq);

void populate_bit_complement_nd(struct Graph *g);

struct Graph *induced_subgraph(struct Graph *g, int *vv, int vv_len);

struct Graph *new_graph(int n);

void free_graph(struct Graph *g);

struct Graph *readGraph(char* filename);

void copy_VtxList(struct VtxList *src, struct VtxList *dest);

void vtxlist_push_vtx(struct Graph *g, struct VtxList *L, int v);

void vtxlist_pop_vtx(struct Graph *g, struct VtxList *L);

#endif
