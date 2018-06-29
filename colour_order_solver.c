#define _GNU_SOURCE
#define _POSIX_SOURCE

#include "c_program_timing.h"
#include "graph.h"
#include "sorting.h"
#include "bitset.h"
#include "vertex_ordering.h"
#include "util.h"
#include "colour_order_solver.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void colouring_bound(struct Graph *g,
        unsigned long long * P_bitset,
        unsigned long long * branch_vv_bitset,
        struct Weight target,
        int numwords)
{
    unsigned long long *to_colour = malloc((g->n+BITS_PER_WORD-1)/BITS_PER_WORD * sizeof *to_colour);
    unsigned long long *candidates = malloc((g->n+BITS_PER_WORD-1)/BITS_PER_WORD * sizeof *candidates);

    copy_bitset(P_bitset, to_colour, numwords);

    int v;
    struct Weight bound = {};

    int *col_class = malloc(g->n * sizeof *col_class);
    struct Weight *residual_wt = malloc(g->n * sizeof *residual_wt);
    for (int i=0; i<g->n; i++)
        residual_wt[i] = g->weight[i];

    int k = 14;
    int K = 14;
    int v_options[14];

    int pc = bitset_popcount(to_colour, numwords);

    while ((v=first_set_bit(to_colour, numwords))!=-1) {
        if (pc >= K) {
            for (int i=0; i<k; i++) {
                int w = first_set_bit(to_colour, numwords);
                v_options[i] = w;
                unset_bit(to_colour, w);
            }
            for (int i=0; i<k; i++) {
                int w = v_options[i];
                set_bit(to_colour, w);
            }
            double rnd = (double)rand() / ((double)RAND_MAX+1);
            if (rnd == 0)
                rnd = 0.00000001;
            long rand_num = (long)k - 1 - (unsigned)(-6 * log(rnd));
            if (rand_num < 0)
                rand_num = 0;
//            printf("%ld\n", rand_num);
            v = v_options[rand_num];
        }
        copy_bitset(to_colour, candidates, numwords);
        struct Weight class_min_wt = residual_wt[v];
        struct Weight class_max_wt = residual_wt[v];
//        unset_bit(to_colour, v);
        int col_class_size = 1;
        col_class[0] = v;
        bitset_intersect_with(candidates, g->bit_complement_nd[v], numwords);
        while ((v=first_set_bit(candidates, numwords))!=-1) {
            if (weight_lt(residual_wt[v], class_min_wt))
                class_min_wt = residual_wt[v];
            if (weight_gt(residual_wt[v], class_max_wt))
                class_max_wt = residual_wt[v];
//            unset_bit(to_colour, v);
            col_class[col_class_size++] = v;
            bitset_intersect_with(candidates, g->bit_complement_nd[v], numwords);
        }
        if (!weight_gt(weight_sum(bound, class_max_wt), target)) {
            bound = weight_sum(bound, class_min_wt);
            for (int i=0; i<col_class_size; i++) {
                int w = col_class[i];
                residual_wt[w] = weight_difference(residual_wt[w], class_min_wt);
                if (weight_eq_zero(residual_wt[w])) {
                    unset_bit(to_colour, w);
                    --pc;
                }
            }
        } else {
            for (int i=0; i<numwords; i++) {
                unsigned long long word = to_colour[i];
                while (word) {
                    int bit = __builtin_ctzll(word);
                    word ^= (1ull << bit);
                    int v = i * BITS_PER_WORD + bit;
                    if (weight_gt(weight_sum(bound, residual_wt[v]), target)) {
                        set_bit(branch_vv_bitset, v);
                        unset_bit(to_colour, v);
                        --pc;
                    }
                }
            }
        }
    }
    free(residual_wt);
    free(col_class);

    free(to_colour);
    free(candidates);
}

void expand(struct Graph *g, struct VtxList *C, unsigned long long *P_bitset,
        struct VtxList *incumbent, int level, long *expand_call_count,
        bool quiet, int numwords)
{
    (*expand_call_count)++;
    if (*expand_call_count % 1000 == 0)
        check_for_timeout();
    if (is_timeout_flag_set()) return;

    if (!quiet && weight_gt(C->total_wt, incumbent->total_wt)) {
        copy_VtxList(C, incumbent);
        long elapsed_msec = get_elapsed_time_msec();
        printf("New incumbent: weight ");
        print_weight(incumbent->total_wt);
        printf(" at time %ld ms after %ld expand calls\n", elapsed_msec, *expand_call_count);
    }

    struct Weight target = weight_difference(incumbent->total_wt, C->total_wt);

    unsigned long long *branch_vv_bitset = malloc(numwords * sizeof *branch_vv_bitset);
    unsigned long long *bvvb = malloc(numwords * sizeof *branch_vv_bitset);

    int top = 1;
    bool can_backtrack = false;
    for (int i=0; i<top; i++) {
        for (int i=0; i<numwords; i++)
            bvvb[i] = 0;
        colouring_bound(g, P_bitset, bvvb, target, numwords);
        if (0 == bitset_popcount(bvvb, numwords)) {
            can_backtrack = true;
            break;
        }
        if (i == 0 || bitset_popcount(bvvb, numwords) < bitset_popcount(branch_vv_bitset, numwords)) {
            copy_bitset(bvvb, branch_vv_bitset, numwords);
            if (C->size != 0)
                top = bitset_popcount(bvvb, numwords) * 2;
        }
    }

    if (!can_backtrack) {
        bitset_intersect_with_complement(P_bitset, branch_vv_bitset, numwords);

        unsigned long long *new_P_bitset = malloc(numwords * sizeof *new_P_bitset);

        int v;
        while ((v=first_set_bit(branch_vv_bitset, numwords))!=-1) {
            unset_bit(branch_vv_bitset, v);

            copy_bitset(P_bitset, new_P_bitset, numwords);
            bitset_intersect_with_complement(new_P_bitset, g->bit_complement_nd[v], numwords);

            vtxlist_push_vtx(g, C, v);
            expand(g, C, new_P_bitset, incumbent, level+1, expand_call_count, quiet, numwords);
            set_bit(P_bitset, v);
            vtxlist_pop_vtx(g, C);
        }
        free(new_P_bitset);
    }

    free(bvvb);
    free(branch_vv_bitset);
}

void mc(struct Graph* g, long *expand_call_count,
        bool quiet, int vtx_ordering, struct VtxList *incumbent)
{
    srand(time(NULL));

    calculate_all_degrees(g);

    int *vv = malloc(g->n * sizeof *vv);
    order_vertices(vv, g, vtx_ordering);

    struct Graph *ordered_graph = induced_subgraph(g, vv, g->n);
    populate_bit_complement_nd(ordered_graph);

    int numwords = (ordered_graph->n + BITS_PER_WORD - 1) / BITS_PER_WORD;

    unsigned long long *P_bitset = calloc(numwords, sizeof *P_bitset);
    for (int v=0; v<ordered_graph->n; v++)
        set_bit(P_bitset, v);

    struct VtxList C;
    init_VtxList(&C, ordered_graph->n);
    expand(ordered_graph, &C, P_bitset, incumbent, 0, expand_call_count, quiet, numwords);
    destroy_VtxList(&C);
    free(P_bitset);

    // Use vertex indices from original graph
    for (int i=0; i<incumbent->size; i++)
        incumbent->vv[i] = vv[incumbent->vv[i]];

    free_graph(ordered_graph);
    free(vv);
}
