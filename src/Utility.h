//
// Created by ssunah on 11/30/16.
//

#ifndef FSC_UTILITY_H
#define FSC_UTILITY_H

#include "Graph.h"
#include <string.h>
#include <sys/time.h>

class Utility {
private:
    Utility() {}

public:
    // Refer to paper: An O(m) Algorithm for Cores Decomposition of Networks.
    // The difference is that the vertices in our implementation are ranged from 0 to n-1.
    static void GetKCore(const Graph* graph, int* core_table);
    static timeval GetTime();
    static unsigned long TimeDiffInMicroseconds(const timeval start, const timeval end);
    static double TimeDiffInSeconds(const timeval start, const timeval end);
    static unsigned long UpperPowerOfTwo(unsigned long v);
    static void old_cheap(int* col_ptrs, int* col_ids, int* match, int* row_match, int n, int m);
    static void match_bfs(int* col_ptrs, int* col_ids, int* match, int* row_match, int* visited,
                          int* queue, int* previous, int n, int m);

};

#endif //FSM_UTILITY_H
