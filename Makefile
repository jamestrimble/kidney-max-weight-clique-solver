HEADERS = graph.h c_program_timing.h sorting.h bitset.h vertex_ordering.h util.h
SHARED_C = graph.c c_program_timing.c vertex_ordering.c util.c bitset.c 

all: colour_order

colour_order: colour_order.c colour_order_solver.c colour_order_solver.h $(SHARED_C) $(HEADERS)
	gcc -O3 -march=native -Wall -std=c11 -o colour_order colour_order.c colour_order_solver.c $(SHARED_C) -lm

clean:
	rm colour_order
