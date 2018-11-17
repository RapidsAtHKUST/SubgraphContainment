//
// Created by ssunah on 5/31/17.
//

#include "GraphDB.h"
#include "GraphReader.h"

#include <iostream>

GraphDB::GraphDB() {
    label_map_ = NULL;
    graph_map_ = NULL;
    data_graph_label_histogram_ = NULL;
    query_graph_label_histogram_ = NULL;
    graphs_ = NULL;
}

GraphDB::~GraphDB() {
    delete label_map_;
    delete graph_map_;
    delete data_graph_label_histogram_;
    delete query_graph_label_histogram_;

    for (int i = 0; i < graphs_->size(); ++i) {
        delete (*graphs_)[i];
    }

    delete graphs_;
}

void GraphDB::BuildDB(const string &file, GraphFormat format) {
    switch (format) {
        case GraphFormat::FA:
            graphs_ = GraphReader::LoadGraphFormatA(file, label_map_, graph_map_, data_graph_label_histogram_);
            break;
        case GraphFormat::FB:
            graphs_ = GraphReader::LoadGraphFormatB(file, label_map_, graph_map_, data_graph_label_histogram_);
            break;
        case GraphFormat::FC:
            graphs_ = GraphReader::LoadGraphFormatC(file, label_map_, graph_map_, data_graph_label_histogram_);
            break;
        default:
            break;
    }

    // Build inverse label index and nlf filter.
    for (int i = 0; i < graphs_->size(); ++i) {
        (*graphs_)[i]->BuildReverseLabelIndex(label_map_->size());
        (*graphs_)[i]->ComputeMaxDegree();
#ifdef ENABLE_NEIGHBOR_LABEL_FREQUENCY_FILTER
        (*graphs_)[i]->BuildNLF();
#endif
    }

}

vector<Graph*>* GraphDB::LoadQueryGraphs(const string &file, GraphFormat format) {
    vector<Graph*>* query_graphs = NULL;
    unordered_map<string, int>* graph_map = NULL;
    switch (format) {
        case GraphFormat::FA:
            query_graphs = GraphReader::LoadGraphFormatA(file, label_map_, graph_map, query_graph_label_histogram_);
            break;
        case GraphFormat::FB:
            query_graphs = GraphReader::LoadGraphFormatB(file, label_map_, graph_map, query_graph_label_histogram_);
            break;
        case GraphFormat::FC:
            query_graphs = GraphReader::LoadGraphFormatC(file, label_map_, graph_map, query_graph_label_histogram_);
            break;
        default:
            break;
    }

    for (int i = 0; i < query_graphs->size(); ++i) {
        (*query_graphs)[i]->ComputeLabelFrequency();
        (*query_graphs)[i]->ComputeMaxDegree();

#ifdef ENABLE_NEIGHBOR_LABEL_FREQUENCY_FILTER
        (*query_graphs)[i]->BuildNLF();
#endif

        (*query_graphs)[i]->SortEdgeByLabel();
        (*query_graphs)[i]->BuildPrunedGraph();
        (*query_graphs)[i]->BuildDegreeOneVertices();
    }

    delete graph_map;
    return query_graphs;
}