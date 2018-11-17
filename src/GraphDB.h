//
// Created by ssunah on 5/31/17.
//

#ifndef FSC_GRAPHDB_H
#define FSC_GRAPHDB_H

#include <unordered_map>
#include <string>
#include "Types.h"
#include "Graph.h"
using namespace std;

class GraphDB {
private:
    unordered_map<string, int>* label_map_;                  // Map label value to label id, label id starts from 0.

    unordered_map<string, int>* graph_map_;                  // Map graph name to graph id, graph id starts from 0.

    unordered_map<int, int>* data_graph_label_histogram_;    // The number of distinct labels in data graphs.

    unordered_map<int, int>* query_graph_label_histogram_;   // The number of distinct labels in query graphs.

    vector<Graph*>* graphs_;                                 // The data graphs stored in db.

public:
    GraphDB();
    ~GraphDB();

    void BuildDB(const string& file, GraphFormat format);
    vector<Graph*>* GetDataGraphs() {
        return graphs_;
    }
    unordered_map<string, int>* GetLabelMap() {
        return label_map_;
    }
    unordered_map<string, int>* GetGraphMap() {
        return graph_map_;
    }
    vector<Graph*>* LoadQueryGraphs(const string& file, GraphFormat format);

    const int LabelsBinSize(const int label) const {
        return (*data_graph_label_histogram_)[label];
    }
};


#endif //SRC_GRAPHDB_H
