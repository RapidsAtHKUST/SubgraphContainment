//
// Created by ssunah on 5/31/17.
//

#ifndef SRC_GRAPHREADER_H
#define SRC_GRAPHREADER_H

#include "Graph.h"

class GraphReader {
public:
    static vector<Graph*>* LoadGraphFormatA(const string& file, unordered_map<string, int>*& label_map,
                                            unordered_map<string, int>*& graph_map, unordered_map<int, int>*& label_histogram);

    static vector<Graph*>* LoadGraphFormatB(const string& file, unordered_map<string, int>*& label_map,
                                            unordered_map<string, int>*& graph_map, unordered_map<int, int>*& label_histogram);

    static vector<Graph*>* LoadGraphFormatC(const string& file, unordered_map<string, int>*& label_map,
                                            unordered_map<string, int>*& graph_map, unordered_map<int, int>*& label_histogram);
};


#endif //SRC_GRAPHREADER_H
