HEADERS = graph.h c_program_timing.h sorting.h bitset.h vertex_ordering.h util.h
SHARED_C = graph.c c_program_timing.c vertex_ordering.c util.c bitset.c 

all: kidney_mwc

kidney_mwc: kidney_mwc.c kidney_mwc_solver.c kidney_mwc_solver.h $(SHARED_C) $(HEADERS)
	gcc -O3 -march=native -Wall -std=c11 -Wno-unused-function -o kidney_mwc kidney_mwc.c kidney_mwc_solver.c $(SHARED_C) -lm

clean:
	rm kidney_mwc
