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

#define NUM_RANDOMS 9973
#define NUM_RANDOM_VALUES 20

int randoms[NUM_RANDOMS];

void make_randoms()
{
    for (int i=0; i<NUM_RANDOMS; i++) {
//        int val = NUM_RANDOM_VALUES - 1;
//        while (rand() > RAND_MAX * 0.2 && val > 0)
//            --val;
        int val = rand() % NUM_RANDOM_VALUES;
        randoms[i] = val;
    }
}

int kth_set_bit(int k, unsigned long long const * const bitset, int numwords)
{
    int num_bits_seen = 0;
    for (int word_index = 0; word_index < numwords; ++word_index) {
        unsigned long long word = bitset[word_index];
        int popcount = __builtin_popcountll(word);
        if (num_bits_seen + popcount >= k) {
            for (;;) {
                int bit = __builtin_ctzll(word);
                ++num_bits_seen;
                if (num_bits_seen == k) {
                    return word_index * BITS_PER_WORD + bit;
                }
                word ^= (1ull << bit);
            }
        }
        num_bits_seen += popcount;
    }
    return -1;
}

void remove_vertices_heavier_than_max_permitted(unsigned long long *to_colour,
        int *pc, unsigned long long *branch_vv_bitset,
        struct Weight *residual_wt, struct Weight max_permitted_weight, int numwords)
{
    for (int i=0; i<numwords; i++) {
        unsigned long long word = to_colour[i];
        while (word) {
            int bit = __builtin_ctzll(word);
            word ^= (1ull << bit);
            int v = i * BITS_PER_WORD + bit;
            if (weight_gt(residual_wt[v], max_permitted_weight)) {
                set_bit(branch_vv_bitset, v);
                unset_bit(to_colour, v);
                --(*pc);
            }
        }
    }
}

void colouring_bound(struct Graph *g,
        unsigned long long * P_bitset,
        unsigned long long * branch_vv_bitset,
        struct Weight target,
        int numwords)
{
    unsigned long long *to_colour = malloc(numwords * sizeof *to_colour);
    unsigned long long *candidates = malloc(numwords * sizeof *candidates);
    unsigned long long *col_class_bitset = malloc(numwords * sizeof *candidates);
    unsigned long long *prev_col_class_bitset = malloc(numwords * sizeof *candidates);
    for (int i=0; i<numwords; i++)
        prev_col_class_bitset[i] = 0;

    copy_bitset(P_bitset, to_colour, numwords);

    int u = -1;
    int v;
    struct Weight bound = {};

    int *col_class = malloc(g->n * sizeof *col_class);
    struct Weight *residual_wt = malloc(g->n * sizeof *residual_wt);
    for (int i=0; i<g->n; i++)
        residual_wt[i] = g->weight[i];

    int pc = bitset_popcount(to_colour, numwords);

    int rand_index = rand() % NUM_RANDOMS;
    int rand_step = 1 + (rand() % (NUM_RANDOMS - 1));

    while (0 != pc) {
        struct Weight max_permitted_weight = weight_difference(target, bound);
        if (u == -1 || !test_bit(to_colour, u)) {
            if (pc >= NUM_RANDOM_VALUES) {
                int rand_num = randoms[rand_index];
                rand_index = (rand_index + rand_step) % NUM_RANDOMS;
                u = kth_set_bit(rand_num + 1, to_colour, numwords);
            } else {
                u=first_set_bit(to_colour, numwords);
            }
        }
        for (int i=0; i<numwords; i++)
            col_class_bitset[i] = 0;
        copy_bitset(to_colour, candidates, numwords);
        if (weight_gt(residual_wt[u], max_permitted_weight)) {
            remove_vertices_heavier_than_max_permitted(to_colour, &pc, branch_vv_bitset,
                    residual_wt, max_permitted_weight, numwords);
            goto next_colour_class;
        }
        struct Weight class_min_wt = residual_wt[u];
        int col_class_size = 1;
        col_class[0] = u;
        set_bit(col_class_bitset, u);
        bitset_intersect_with(candidates, g->bit_complement_nd[u], numwords);
        v = 0;

        while ((v=first_set_bit_from_word(candidates, v/BITS_PER_WORD, numwords))!=-1) {
            if (weight_lt(residual_wt[v], class_min_wt)) {
                class_min_wt = residual_wt[v];
            }
            if (weight_gt(residual_wt[v], max_permitted_weight)) {
                remove_vertices_heavier_than_max_permitted(to_colour, &pc, branch_vv_bitset,
                        residual_wt, max_permitted_weight, numwords);
                goto next_colour_class;
            }
            col_class[col_class_size++] = v;
            set_bit(col_class_bitset, v);

            bitset_intersect_with_from_word(candidates, g->bit_complement_nd[v], v/BITS_PER_WORD, numwords);
            if (union_is_subset_of(col_class_bitset, candidates, prev_col_class_bitset, numwords)) {
                bitset_union_with(col_class_bitset, candidates, numwords);
                for (int i=v/BITS_PER_WORD; i<numwords; i++) {
                    unsigned long long word = candidates[i];
                    while (word) {
                        int bit = __builtin_ctzll(word);
                        word ^= (1ull << bit);
                        int v = i * BITS_PER_WORD + bit;
                        if (weight_lt(residual_wt[v], class_min_wt)) {
                            class_min_wt = residual_wt[v];
                        }
                        if (weight_gt(residual_wt[v], max_permitted_weight)) {
                            remove_vertices_heavier_than_max_permitted(to_colour, &pc, branch_vv_bitset,
                                    residual_wt, max_permitted_weight, numwords);
                            goto next_colour_class;
                        }
                        col_class[col_class_size++] = v;
                    }
                }
                break;
            }
        }
        bound = weight_sum(bound, class_min_wt);
        for (int i=0; i<numwords; i++)
            prev_col_class_bitset[i] = 0;
        for (int i=0; i<col_class_size; i++) {
            int w = col_class[i];
            set_bit(prev_col_class_bitset, w);
            weight_subtract(&residual_wt[w], class_min_wt);
            if (weight_eq_zero(residual_wt[w])) {
                unset_bit(to_colour, w);
                --pc;
            }
        }
        copy_bitset(col_class_bitset, prev_col_class_bitset, numwords);
next_colour_class:
        ;
    }
    free(residual_wt);
    free(col_class);

    free(to_colour);
    free(candidates);
    free(col_class_bitset);
    free(prev_col_class_bitset);
}

void expand(struct Graph *g, struct VtxList *C, unsigned long long *P_bitset,
        struct VtxList *incumbent, int level,
        long *expand_call_count,
        long *colouring_count,
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
        for (int j=0; j<numwords; j++)
            bvvb[j] = 0;
        colouring_bound(g, P_bitset, bvvb, target, numwords);
        ++(*colouring_count);
        int pc = bitset_popcount(bvvb, numwords);
        if (0 == pc) {
            can_backtrack = true;
            break;
        }
        if (i == 0 || pc < bitset_popcount(branch_vv_bitset, numwords)) {
            copy_bitset(bvvb, branch_vv_bitset, numwords);
            if (!weight_eq_zero(incumbent->total_wt)) {
                top = pc * 5;
            }
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
            expand(g, C, new_P_bitset, incumbent, level+1, expand_call_count, colouring_count, quiet, numwords);
            set_bit(P_bitset, v);
            vtxlist_pop_vtx(g, C);
        }
        free(new_P_bitset);
    }

    free(bvvb);
    free(branch_vv_bitset);
}

void mc(struct Graph* g, long *expand_call_count, long *colouring_count,
        bool quiet, int vtx_ordering, struct VtxList *incumbent)
{
//    srand(time(NULL));
    make_randoms();

    calculate_all_degrees(g);

    int *vv = malloc(g->n * sizeof *vv);
    order_vertices(vv, g, vtx_ordering);

    struct Graph *ordered_graph = induced_subgraph(g, vv, g->n);
    populate_bit_complement_nd(ordered_graph);

    int numwords = (ordered_graph->n + BITS_PER_WORD - 1) / BITS_PER_WORD;

    struct Weight *saved_weights = malloc(ordered_graph->n * sizeof *saved_weights);
    for (int i=0; i<ordered_graph->n; i++) {
        saved_weights[i] = ordered_graph->weight[i];
        for (int j=1; j<WEIGHT_SIZE; j++) {
            ordered_graph->weight[i].weight[j] = 0;
        }
    }

    for (int i=0; i<WEIGHT_SIZE; i++) {
        if (i != 0) {
            for (int j=0; j<ordered_graph->n; j++) {
                ordered_graph->weight[j].weight[i] = saved_weights[j].weight[i];
            }
        }
        unsigned long long *P_bitset = calloc(numwords, sizeof *P_bitset);
        for (int v=0; v<ordered_graph->n; v++)
            set_bit(P_bitset, v);

        struct VtxList C;
        init_VtxList(&C, ordered_graph->n);
        expand(ordered_graph, &C, P_bitset, incumbent, 0, expand_call_count, colouring_count, quiet, numwords);
        destroy_VtxList(&C);
        free(P_bitset);
    }

    // Use vertex indices from original graph
    for (int i=0; i<incumbent->size; i++)
        incumbent->vv[i] = vv[incumbent->vv[i]];

    free_graph(ordered_graph);
    free(vv);
}
