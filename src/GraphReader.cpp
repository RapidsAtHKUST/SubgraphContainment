//
// Created by ssunah on 5/31/17.
//

#include "GraphReader.h"
#include <fstream>
using namespace std;

vector<Graph *> *GraphReader::LoadGraphFormatA(const string &file_path, unordered_map<string, int> *&label_map,
                                                unordered_map<string, int>*& graph_map, unordered_map<int, int>*& label_histogram) {
    ifstream infile(file_path.c_str());

    if (infile) {
        printf("Open graph file %s successfully.\n", file_path.c_str());
    }
    else {
        printf("Cannot open graph file %s.\n", file_path.c_str());
        exit(-1);
    }

    if (label_map == NULL) {
        label_map = new unordered_map<string, int>();
    }

    graph_map = new unordered_map<string, int>();
    label_histogram = new unordered_map<int, int>();

    Graph* graph = NULL;
    vector<Graph*>* graph_list = new vector<Graph*>();
    vector<pair<int, int> > edges;
    vector<int> degree;
    int type = 1;
    int vertex_count = 0;
    int edge_count = 0;
    while(!infile.eof()) {
        switch (type) {
            case 1: {
                string graph_name;
                infile >> graph_name;
                if (graph_name.empty())
                    break;

                (*graph_map)[graph_name] = graph_list->size();
                graph = new Graph();
                graph_list->push_back(graph);
                vertex_count = 0;
                edge_count = 0;
                edges.clear();
                degree.clear();
                type = 2;
                break;
            }
            case 2: {
                int vertex_num;
                infile >> vertex_num;
                degree.resize(vertex_num, 0);
                graph->SetVerticesCount(vertex_num);
                type = 3;
                break;
            }
            case 3: {
                string label;
                infile >> label;
                int label_id;
                unordered_map<string, int>::iterator iter = label_map->find(label);
                if (iter == label_map->end()) {
                    label_id = label_map->size();
                    (*label_map)[label] = label_id;
                }
                else {
                    label_id = iter->second;
                }

                unordered_map<int, int>::iterator label_histogram_iter = label_histogram->find(label_id);
                if (label_histogram_iter == label_histogram->end()) {
                    (*label_histogram)[label_id] = 1;
                }
                else {
                    (*label_histogram)[label_id] += 1;
                }

                graph->SetLabel(vertex_count, label_id);
                vertex_count += 1;
                if (vertex_count >= graph->VerticesCount())
                    type = 4;
                break;
            }
            case 4: {
                int edge_num;
                infile >> edge_num;
                graph->SetEdgesCount(edge_num);
                edges.reserve(edge_num);
                type = 5;
                break;
            }
            case 5: {
                int begin, end;
                infile >> begin >> end;
                degree[begin] += 1;
                degree[end] += 1;
                edges.push_back(make_pair(begin, end));

                edge_count += 1;
                if (edge_count >= graph->EdgesCount()) {
                    int offset = 0;

                    for (int i = 0; i < graph->VerticesCount(); ++i) {
                        graph->SetOffset(i, offset);
                        offset += degree[i];
                    }

                    graph->SetOffset(graph->VerticesCount(), offset);
                    fill(degree.begin(), degree.end(), 0);

                    for (int i = 0; i < edges.size(); ++i) {
                        begin = edges[i].first;
                        end = edges[i].second;

                        graph->SetNeighbor(graph->GetOffset(begin) + degree[begin], end);
                        graph->SetNeighbor(graph->GetOffset(end) + degree[end], begin);

                        degree[begin] += 1;
                        degree[end] += 1;
                    }

                    graph->SortEdge();
                    type = 1;
                }

                break;
            }
            default: {
                exit(-1);
            }
        }
    }

    return graph_list;
}

vector<Graph *> *GraphReader::LoadGraphFormatB(const string &file_path, unordered_map<string, int> *&label_map,
                                               unordered_map<string, int>*& graph_map, unordered_map<int, int>*& label_histogram) {
    ifstream infile(file_path.c_str());

    if (infile) {
        printf("Open graph file %s successfully.\n", file_path.c_str());
    }
    else {
        printf("Cannot open graph file %s.\n", file_path.c_str());
        exit(-1);
    }

    if (label_map == NULL) {
        label_map = new unordered_map<string, int>();
    }

    graph_map = new unordered_map<string, int>();
    label_histogram = new unordered_map<int, int>();

    Graph* graph = NULL;
    vector<Graph*>* graph_list = new vector<Graph*>();
    vector<pair<int, int> > edges;
    vector<int> vertices;
    while(!infile.eof()) {
        char type;
        infile >> type;
        if (infile.eof())
            break;
        switch (type) {
            case 't': {
                if (graph != NULL) {
                    graph->SetEdgesCount(edges.size());
                    int offset = 0;
                    for (int i = 0; i < graph->VerticesCount(); ++i) {
                        graph->SetOffset(i, offset);
                        offset += vertices[i];
                    }

                    graph->SetOffset(graph->VerticesCount(), offset);
                    fill(vertices.begin(), vertices.end(), 0);

                    for (int i = 0; i < edges.size(); ++i) {
                        int begin = edges[i].first;
                        int end = edges[i].second;

                        graph->SetNeighbor(graph->GetOffset(begin) + vertices[begin], end);
                        graph->SetNeighbor(graph->GetOffset(end) + vertices[end], begin);

                        vertices[begin] += 1;
                        vertices[end] += 1;
                    }

                    graph->SortEdge();
                }

                char symbol;
                string graph_name;
                infile >> symbol >> graph_name;
                if (graph_name.empty())
                    break;

                (*graph_map)[graph_name] = graph_list->size();
                graph = new Graph();
                graph_list->push_back(graph);
                edges.clear();
                vertices.clear();
                break;
            }
            case 'v': {
                int id;
                int label_id;
                string label;
                infile >> id >> label;

                unordered_map<string, int>::iterator iter = label_map->find(label);
                if (iter == label_map->end()) {
                    label_id = label_map->size();
                    (*label_map)[label] = label_id;
                }
                else {
                    label_id = iter->second;
                }

                unordered_map<int, int>::iterator label_histogram_iter = label_histogram->find(label_id);
                if (label_histogram_iter == label_histogram->end()) {
                    (*label_histogram)[label_id] = 1;
                }
                else {
                    (*label_histogram)[label_id] += 1;
                }

                vertices.push_back(label_id);
                break;
            }
            case 'e': {
                if (edges.size() == 0) {
                    graph->SetVerticesCount(vertices.size());
                    for (int i = 0; i < vertices.size(); ++i) {
                        graph->SetLabel(i, vertices[i]);
                    }

                    fill(vertices.begin(), vertices.end(), 0);
                }
                int begin_id;
                int end_id;
                int label;
                infile >> begin_id >> end_id >> label;

                vertices[begin_id] += 1;
                vertices[end_id] += 1;

                edges.push_back(make_pair(begin_id, end_id));

                break;
            }
            default: {
                exit(-1);
            }
        }
    }

    if (graph != NULL) {
        graph->SetEdgesCount(edges.size());
        int offset = 0;

        for (int i = 0; i < graph->VerticesCount(); ++i) {
            graph->SetOffset(i, offset);
            offset += vertices[i];
        }

        graph->SetOffset(graph->VerticesCount(), offset);
        fill(vertices.begin(), vertices.end(), 0);

        for (int i = 0; i < edges.size(); ++i) {
            int begin = edges[i].first;
            int end = edges[i].second;

            graph->SetNeighbor(graph->GetOffset(begin) + vertices[begin], end);
            graph->SetNeighbor(graph->GetOffset(end) + vertices[end], begin);

            vertices[begin] += 1;
            vertices[end] += 1;
        }

        graph->SortEdge();
    }

    return graph_list;
}

vector<Graph *> *GraphReader::LoadGraphFormatC(const string &file_path, unordered_map<string, int> *&label_map,
                                               unordered_map<string, int>*& graph_map, unordered_map<int, int>*& label_histogram) {
    ifstream infile(file_path.c_str());

    if (infile) {
        printf("Open graph file %s successfully.\n", file_path.c_str());
    }
    else {
        printf("Cannot open graph file %s.\n", file_path.c_str());
        exit(-1);
    }

    if (label_map == NULL) {
        label_map = new unordered_map<string, int>();
    }

    graph_map = new unordered_map<string, int>();
    label_histogram = new unordered_map<int, int>();

    Graph* graph = NULL;
    vector<Graph*>* graph_list = new vector<Graph*>();
    vector<int> vertices;
    int offset = 0;
    while(!infile.eof()) {
        char type;
        infile >> type;
        switch (type) {
            case 't': {
                string graph_name;
                int vertices_num;
                int edges_num;
                infile >> graph_name >> vertices_num >> edges_num;
                if (graph_name.empty())
                    break;

                if (graph != NULL) {
                    graph->SetOffset(graph->VerticesCount(), offset);
                }
                (*graph_map)[graph_name] = graph_list->size();
                graph = new Graph();
                graph->SetVerticesCount(vertices_num);
                graph->SetEdgesCount(edges_num);
                graph_list->push_back(graph);
                offset = 0;
                vertices.resize(vertices_num);
                fill(vertices.begin(), vertices.end(), 0);
                break;
            }
            case 'v': {
                int id;
                int label_id;
                int degree;
                string label;
                infile >> id >> label >> degree;

                unordered_map<string, int>::iterator iter = label_map->find(label);
                if (iter == label_map->end()) {
                    label_id = label_map->size();
                    (*label_map)[label] = label_id;
                }
                else {
                    label_id = iter->second;
                }

                unordered_map<int, int>::iterator label_histogram_iter = label_histogram->find(label_id);
                if (label_histogram_iter == label_histogram->end()) {
                    (*label_histogram)[label_id] = 1;
                }
                else {
                    (*label_histogram)[label_id] += 1;
                }

                graph->SetLabel(id, label_id);
                graph->SetOffset(id, offset);
                offset += degree;
                break;
            }
            case 'e': {
                int begin_id;
                int end_id;

                infile >> begin_id >> end_id;

                graph->SetNeighbor(graph->GetOffset(begin_id) + vertices[begin_id], end_id);
                graph->SetNeighbor(graph->GetOffset(end_id) + vertices[end_id], begin_id);

                vertices[begin_id] += 1;
                vertices[end_id] += 1;

                break;
            }
            default: {
                exit(-1);
            }
        }
    }

    if (graph != NULL) {
        graph->SetOffset(graph->VerticesCount(), offset);
    }

    for (int i = 0; i < graph_list->size(); ++i){
        (*graph_list)[i]->SortEdge();
    }

    return graph_list;
}
