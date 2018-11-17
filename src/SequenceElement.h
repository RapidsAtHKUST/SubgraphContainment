//
// Created by ssunah on 11/29/16.
//

#ifndef FSM_MATCHINGELEMENT_H
#define FSM_MATCHINGELEMENT_H

#include "Config.h"
#include <stddef.h>
#include <utility>

using namespace std;


class SequenceElement {
public:
    int id_;                    // The id of query vertex.
    int pivot_;                 // The pivot of query vertex.
    int bn_count_;              // The number of this vertex's neighbors that have been matched.
    int* bn_;                   // The neighbors of this vertex, its value is the vertex id.
public:
    SequenceElement() {
        bn_ = new int[MAX_QUERY_GRAPH_VERTICES_COUNT];
    }

    ~SequenceElement() {
        delete[] bn_;
        bn_ = NULL;
    }
};

#endif //FSM_MATCHINGELEMENT_H
