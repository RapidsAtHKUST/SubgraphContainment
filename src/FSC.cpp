//
// Created by ssunah on 6/13/17.
//

#include "FSC.h"
#include "Utility.h"
#include <iostream>
#include <queue>
using namespace std;

FSC::FSC(GraphDB* graphDB) {
    graphDB_ = graphDB;
    InitializeDataGraphResource(graphDB_->GetDataGraphs());
}

FSC::~FSC() {
    ClearDataGraphResource();
}

void FSC::Query(const Graph *query_graph, int *query_results, int &results_count) {

    InitializeQueryGraphResource(query_graph);

    filtering_time_ = 0;
    verification_time_ = 0;
    ordering_time_ = 0;
    enumeration_time_ = 0;
    save_time_ = 0;

    for (int i = 0; i < (int)data_graphs_->size(); ++i) {
        Reset();

        const Graph *data_graph = data_graphs_->at(i);

        timeval simulation_start = Utility::GetTime();

        bool is_valid = Simulation(query_graph, data_graph);

        timeval simulation_end = Utility::GetTime();
        double temp_filtering_time = Utility::TimeDiffInMicroseconds(simulation_start, simulation_end);
        filtering_time_ += temp_filtering_time;

        if (is_valid) {
            candidate_graph_num_ += 1;

            timeval verification_start = Utility::GetTime();

            is_valid = Verification(query_graph, data_graph);

            timeval verification_end = Utility::GetTime();
            double temp_verification_time = Utility::TimeDiffInMicroseconds(verification_start, verification_end);

            verification_time_ += temp_verification_time;
            ordering_time_ += temp_ordering_time_;
            enumeration_time_ += temp_verification_time - temp_ordering_time_;
            if (is_valid) {
                query_results[valid_graph_num_++] = i;
            } else {
                save_time_ += temp_verification_time;
            }
        }

    }
    cout << endl;

    GetMemoryCost();
    double filtering_precision = (valid_graph_num_) / (double)candidate_graph_num_;

    processing_time_ = filtering_time_ + verification_time_;
    cout << "Candidate Graphs Count :" << candidate_graph_num_ << " ." << endl;
    cout << "Answer Graphs Count :" << valid_graph_num_ << " ." << endl;
    cout << "Filtering Time : " << filtering_time_ << " us." << endl;
    cout << "Ordering Time : " << ordering_time_ << " us." << endl;
    cout << "Enumeration Time : " << enumeration_time_ << " us." << endl;
    cout << "Verification Time : " << verification_time_ << " us." << endl;
    cout << "Per Verification Time : " << verification_time_ / candidate_graph_num_  << " us." << endl;
    cout << "Processing Time : " << processing_time_ << " us." << endl;
    cout << "Verification Time due to false prediction : " << save_time_ << " us." << endl;
    cout << "Filtering Precision : " << filtering_precision << " ." << endl;

    total_candidate_num_ += candidate_graph_num_;
    total_answer_num_ += valid_graph_num_;
    total_time_ += processing_time_;
    total_enumeration_time_ += enumeration_time_;
    total_ordering_time_ += ordering_time_;
    total_save_time_ += save_time_;
    total_filtering_time_ += filtering_time_;
    total_verification_time_ += verification_time_;
    total_per_verification_time_ += verification_time_ / candidate_graph_num_;
    total_filtering_precision += filtering_precision;

    results_count = valid_graph_num_;
    ClearQueryGraphResource();
}

void FSC::InitializeDataGraphResource(const vector<Graph *> *data_graphs) {
    data_graphs_ = data_graphs;
    total_label_num_ = graphDB_->GetLabelMap()->size();
    total_time_ = 0;
    total_filtering_precision = 0;
    total_filtering_time_ = 0;
    total_verification_time_ = 0;
    total_per_verification_time_ = 0;
    total_ordering_time_ = 0;
    total_save_time_ = 0;
    total_candidate_num_ = 0;
    total_answer_num_ = 0;
    peak_memory_cost_ = 0;
    max_data_vertices_num_ = -1;
    max_data_vertices_with_same_label_ = -1;
    max_data_edges_num_ = -1;
    max_degree_ = -1;
    for (int i = 0; i < (int)data_graphs->size(); ++i) {
        if ((*data_graphs)[i]->MaxVerticesNumWithSameLabel() > max_data_vertices_with_same_label_) {
            max_data_vertices_with_same_label_ = (*data_graphs)[i]->MaxVerticesNumWithSameLabel();
        }

        if ((*data_graphs)[i]->VerticesCount() > max_data_vertices_num_) {
            max_data_vertices_num_ = (*data_graphs)[i]->VerticesCount();
        }

        if ((*data_graphs)[i]->EdgesCount() > max_data_edges_num_) {
            max_data_edges_num_ = (*data_graphs)[i]->EdgesCount();
        }

        if ((*data_graphs)[i]->MaxDegree() > max_degree_) {
            max_degree_ = (*data_graphs)[i]->MaxDegree();
        }
    }

    max_data_edges_num_ = max_data_edges_num_ * 2;

    data_vertices_flag_ = new int[max_data_vertices_num_];
    updated_data_vertices_flag_ = new int[max_data_vertices_num_];
    label_frequency_ = new int[total_label_num_];

    data_vertex_distinct_neighbor_count_ = new int[max_data_vertices_num_];
    data_vertex_visited_update_ = new int[max_data_vertices_num_];
}

void FSC::InitializeQueryGraphResource(const Graph *query_graph) {
    candidate_graph_num_ = 0;
    valid_graph_num_ = 0;

    simulation_order_count_ = query_graph->VerticesCount();

    query_vertices_flags_ = new bool[simulation_order_count_];
    candidates_count_ = new int[simulation_order_count_];
    candidates_ = new int*[simulation_order_count_];
    for (int i = 0; i < simulation_order_count_; ++i) {
        candidates_[i] = new int[max_data_vertices_with_same_label_];
    }

    query_vertices_cores_ = new int[simulation_order_count_];
    backward_value_ = new int[simulation_order_count_];
    memset(query_vertices_cores_, 0, sizeof(int) * simulation_order_count_);
    Utility::GetKCore(query_graph, query_vertices_cores_);

    core_length_ = 0;
    for (int i = 0; i < simulation_order_count_; ++i) {
        if (query_vertices_cores_[i] > 1) {
            core_length_ += 1;
        }
    }

    sequence_ = new SequenceElement[simulation_order_count_];

    idx_ = new int[simulation_order_count_];
    matching_ = new int[simulation_order_count_];
    embedding_ = new int[simulation_order_count_];
    data_vertex_visited_ = new bool[max_data_vertices_num_];

    candidates_index_count_ = new int[simulation_order_count_];
    candidates_index_ = new const int*[simulation_order_count_];
    candidates_index_matching_ = new int[simulation_order_count_];

    weight_ = new double[simulation_order_count_];
    pivot_ = new int[simulation_order_count_];
    directed_weighted_graph_ = new double[simulation_order_count_ * simulation_order_count_];

    target_ = new int[simulation_order_count_];
    result_ = new int[simulation_order_count_];

    bigraph_offset_ = new int*[simulation_order_count_];
    bigraph_edges_ = new int*[simulation_order_count_];

    for (int i = 0; i < simulation_order_count_; ++i) {
        bigraph_offset_[i] = new int[max_data_vertices_with_same_label_ + 1];
        bigraph_edges_[i] = new int[max_data_edges_num_];
    }

    // Data structures used for CFL
    tree_node_ = new TreeNode[simulation_order_count_];
    bfs_order_ = new int[simulation_order_count_];
    bfs_order_index_ = new int[simulation_order_count_];
    level_offset_ = new int[simulation_order_count_ + 1];
    path_cost_ = new size_t*[simulation_order_count_];

    for (int i = 0; i < simulation_order_count_; ++i) {
        path_cost_[i] = new size_t[max_data_vertices_with_same_label_];
    }

    ResetQueryGraphResource();
}

void FSC::ClearQueryGraphResource() {
    delete[] query_vertices_flags_;
    delete[] candidates_count_;
    for (int i = 0; i < simulation_order_count_; ++i) {
        delete[] candidates_[i];
        delete[] bigraph_offset_[i];
        delete[] bigraph_edges_[i];
        delete[] path_cost_[i];
    }

    delete[] path_cost_;
    delete[] bigraph_offset_;
    delete[] bigraph_edges_;
    delete[] candidates_;

    delete[] sequence_;

    delete[] query_vertices_cores_;
    delete[] backward_value_;

    delete[] idx_;
    delete[] matching_;
    delete[] embedding_;
    delete[] data_vertex_visited_;

    delete[] candidates_index_;
    delete[] candidates_index_count_;
    delete[] candidates_index_matching_;

    delete[] tree_node_;
    delete[] bfs_order_;
    delete[] bfs_order_index_;
    delete[] level_offset_;

    delete[] weight_;
    delete[] pivot_;
    delete[] directed_weighted_graph_;

    delete[] target_;
    target_ = NULL;
    delete[] result_;
    result_ = NULL;
}

void FSC::ClearDataGraphResource() {
    delete[] updated_data_vertices_flag_;
    delete[] data_vertices_flag_;
    delete[] label_frequency_;

    delete[] data_vertex_distinct_neighbor_count_;
    delete[] data_vertex_visited_update_;
}

void FSC::Reset() {
    updated_data_vertices_count_ = 0;
    core_paths_.clear();
    tree_paths_.clear();

    memset(directed_weighted_graph_, 0, sizeof(double) * simulation_order_count_ * simulation_order_count_);
    memset(query_vertices_flags_, false, sizeof(bool) * simulation_order_count_);
    memset(data_vertices_flag_, 0, sizeof(int) * max_data_vertices_num_);
    memset(updated_data_vertices_flag_, 0, sizeof(int) * max_data_vertices_num_);
    memset(candidates_count_, 0, sizeof(int) * simulation_order_count_);
    memset(label_frequency_, 0, sizeof(int) * total_label_num_);
    memset(data_vertex_visited_, 0, sizeof(bool) * max_data_vertices_num_);
    memset(data_vertex_distinct_neighbor_count_, 0, sizeof(int) * max_data_vertices_num_);
    memset(data_vertex_visited_update_, 0, sizeof(int) * max_data_vertices_num_);
}

void FSC::ResetQueryGraphResource() {
    memset(query_vertices_flags_, false, sizeof(bool) * simulation_order_count_);
    memset(backward_value_, 0, sizeof(int) * simulation_order_count_);
}

bool FSC::GraphMetadataFilter(const Graph *query_graph, const Graph *data_graph) {
    if (query_graph->VerticesCount() > data_graph->VerticesCount() ||
        query_graph->EdgesCount() > data_graph->EdgesCount() ||
        query_graph->MaxDegree() > data_graph->MaxDegree() ||
        query_graph->LabelsCount() > data_graph->LabelsCount()) {
        return false;
    }

    const pair<int, int>* query_graph_label_frequency = query_graph->LabelFrequency();
    for (int i = 0; i < query_graph->LabelsCount(); ++i) {
        if (query_graph_label_frequency[i].second > data_graph->VerticesCount(query_graph_label_frequency[i].first)) {
            return false;
        }
    }

    return true;
}

bool FSC::Simulation(const Graph *query_graph, const Graph *data_graph) {
    // We first filter data graphs based on the static information that can be obtained in advance.
    if (!GraphMetadataFilter(query_graph, data_graph))
        return false;

    // Generate BFS tree.
    GenerateBFSTree(query_graph, data_graph);

    // Get the candidates of root.
    if (!GenerateRootCandidate(data_graph))
        return false;

    // Top-down generation level by level.
    for (int level = 1; level < level_count_; ++level) {
        // Forward generation
        for (int i = level_offset_[level]; i < level_offset_[level + 1]; ++i) {
            const int cur_vertex = bfs_order_[i];
            TreeNode& node = tree_node_[cur_vertex];

            int count = 0;
            // Loop each backward neighbors of current query vertex.
            for (int j = 0; j < node.bn_count_; ++j) {
                int query_vertex_neighbor = node.bn_[j];

                // Loop each candidate of backward neighbors
                for (int k = 0; k < candidates_count_[query_vertex_neighbor]; ++k) {
                    int backward_neighbor_candidate = candidates_[query_vertex_neighbor][k];

                    if (backward_neighbor_candidate == -1)
                        continue;

                    int neighbor_count;
                    const int* neighbors = data_graph->Neighbors(backward_neighbor_candidate, neighbor_count);

                    // Loop each neighbors of candidate
                    for (int l = 0; l < neighbor_count; ++l) {
                        int neighbor = neighbors[l];

                        if (data_vertices_flag_[neighbor] == count &&
                            data_graph->Label(neighbor) == node.label_ &&
                            data_graph->Degree(neighbor) >= node.degree_) {

                            data_vertices_flag_[neighbor] += 1;

                            if (count == 0)
                                updated_data_vertices_flag_[updated_data_vertices_count_++] = neighbor;
                        }
                    }
                }

                count += 1;
            }

            bool has_candidate = false;
            // Generate candidate for current query vertex.
            for (int j = 0; j < updated_data_vertices_count_; ++j) {
                int candidate = updated_data_vertices_flag_[j];
                if (data_vertices_flag_[candidate] == count) {
#ifdef ENABLE_NEIGHBOR_LABEL_FREQUENCY_FILTER
                    bool is_valid = true;
                    const unordered_map<int, int>* data_vertex_nlf = data_graph->NLF(candidate);

                    if (data_vertex_nlf->size() >= (size_t)node.nlf_label_count_) {
                        for (int k = 0; k < node.nlf_label_count_; ++k) {
                            pair<int, int> label_count = node.nlf_[k];
                            auto data_vertex_iter = data_vertex_nlf->find(label_count.first);
                            if (data_vertex_iter == data_vertex_nlf->end() ||
                                data_vertex_iter->second < label_count.second) {
                                is_valid = false;
                                break;
                            }
                        }
                    }
                    else {
                        is_valid = false;
                    }

                    if (is_valid) {
                        has_candidate = true;
                        candidates_[node.id_][candidates_count_[node.id_]++] = candidate;
                    }
#endif
#ifndef ENABLE_NEIGHBOR_LABEL_FREQUENCY_FILTER
                    has_candidate = true;
                    candidates_[node.id_][candidates_count_[node.id_]++] = candidate;
#endif
                }
            }

            // Clear the updated flag count.
            for (int j = 0; j < updated_data_vertices_count_; ++j) {
                int candidate = updated_data_vertices_flag_[j];
                data_vertices_flag_[candidate] = 0;
            }
            updated_data_vertices_count_ = 0;

            if (!has_candidate) {
                return false;
            }
        }

        // Backward prune.
        for (int i = level_offset_[level + 1] - 1; i >= level_offset_[level]; --i) {
            const int cur_vertex = bfs_order_[i];
            TreeNode& node = tree_node_[cur_vertex];

            int count = 0;
            // Loop each forward neighbors of current query vertex.
            for (int j = 0; j < node.fn_count_; ++j) {
                int query_vertex_neighbor = node.fn_[j];

                // Loop each candidate of forward neighbors
                for (int k = 0; k < candidates_count_[query_vertex_neighbor]; ++k) {
                    int fn_candidate_neighbor = candidates_[query_vertex_neighbor][k];

                    if (fn_candidate_neighbor == -1)
                        continue;

                    int neighbor_count;
                    const int* neighbors = data_graph->Neighbors(fn_candidate_neighbor, neighbor_count);

                    // Loop each neighbors of candidate
                    for (int l = 0; l < neighbor_count; ++l) {
                        int neighbor = neighbors[l];

                        if (data_vertices_flag_[neighbor] == count &&
                            data_graph->Label(neighbor) == node.label_ &&
                            data_graph->Degree(neighbor) >= node.degree_) {

                            data_vertices_flag_[neighbor] += 1;

                            if (count == 0)
                                updated_data_vertices_flag_[updated_data_vertices_count_++] = neighbor;
                        }
                    }
                }

                count += 1;
            }

            bool has_candidate = false;
            for (int j = 0; j < candidates_count_[node.id_]; ++j) {
                if (candidates_[node.id_][j] == -1)
                    continue;

                if (data_vertices_flag_[candidates_[node.id_][j]] != count) {
                    candidates_[node.id_][j] = -1;
                }
                else {
                    has_candidate = true;
                }
            }

            // Clear the updated flag count.
            for (int j = 0; j < updated_data_vertices_count_; ++j) {
                int candidate = updated_data_vertices_flag_[j];
                data_vertices_flag_[candidate] = 0;
            }
            updated_data_vertices_count_ = 0;

            if (!has_candidate) {
                return false;
            }
        }
    }

    // Bottom-up refinement.
    for (int level = level_count_ - 2; level >= 0; --level) {

        for (int i = level_offset_[level]; i < level_offset_[level + 1]; ++i) {
            const int cur_vertex = bfs_order_[i];
            TreeNode& node = tree_node_[cur_vertex];

            int count = 0;
            // Loop each children of current query vertex.
            for (int j = 0; j < node.under_level_count_; ++j) {
                int query_vertex_neighbor = node.under_level_[j];

                // Loop each candidate of children neighbors
                for (int k = 0; k < candidates_count_[query_vertex_neighbor]; ++k) {
                    int fn_candidate_neighbor = candidates_[query_vertex_neighbor][k];

                    if (fn_candidate_neighbor == -1)
                        continue;

                    int neighbor_count;
                    const int* neighbors = data_graph->Neighbors(fn_candidate_neighbor, neighbor_count);

                    // Loop each neighbors of candidate
                    for (int l = 0; l < neighbor_count; ++l) {
                        int neighbor = neighbors[l];

                        if (data_vertices_flag_[neighbor] == count &&
                            data_graph->Label(neighbor) == node.label_ &&
                            data_graph->Degree(neighbor) >= node.degree_) {

                            data_vertices_flag_[neighbor] += 1;

                            if (count == 0)
                                updated_data_vertices_flag_[updated_data_vertices_count_++] = neighbor;
                        }
                    }
                }

                count += 1;
            }

            bool has_candidate = false;
            for (int j = 0; j < candidates_count_[node.id_]; ++j) {
                if (candidates_[node.id_][j] == -1)
                    continue;

                if (data_vertices_flag_[candidates_[node.id_][j]] != count) {
                    candidates_[node.id_][j] = -1;
                }
                else {
                    has_candidate = true;
                }
            }

            // Clear the updated flag count.
            for (int j = 0; j < updated_data_vertices_count_; ++j) {
                int candidate = updated_data_vertices_flag_[j];
                data_vertices_flag_[candidate] = 0;
            }
            updated_data_vertices_count_ = 0;

            if (!has_candidate) {
                return false;
            }
        }
    }

    // Compact each candidate.
    for (int i = 0; i < simulation_order_count_; ++i) {
        int cur_vertex = bfs_order_[i];
        TreeNode& node = tree_node_[cur_vertex];
        int next_position = 0;
        for (int j = 0; j < candidates_count_[node.id_]; ++j) {
            int candidate = candidates_[node.id_][j];
            if (candidate != -1) {
                if (j != next_position) {
                    candidates_[node.id_][next_position] = candidate;
                }
                next_position += 1;

                if (data_vertices_flag_[candidate] == 0) {
                    label_frequency_[node.label_] += 1;
                    data_vertices_flag_[candidate] = 1;
                    updated_data_vertices_flag_[updated_data_vertices_count_++] = candidate;
                }
            }
        }

        if (next_position == 0) {
            for (int j = 0; j < updated_data_vertices_count_; ++j) {
                int candidate = updated_data_vertices_flag_[j];
                data_vertices_flag_[candidate] = 0;
            }
            updated_data_vertices_count_ = 0;
            return false;
        }

        candidates_count_[node.id_] = next_position;
    }

    const pair<int, int>* query_graph_label_frequency = query_graph->LabelFrequency();
    for (int i = 0; i < query_graph->LabelsCount(); ++i) {
        if (label_frequency_[query_graph_label_frequency->first] < query_graph_label_frequency->second) {
            for (int j = 0; j < updated_data_vertices_count_; ++j) {
                int candidate = updated_data_vertices_flag_[j];
                data_vertices_flag_[candidate] = 0;
            }
            updated_data_vertices_count_ = 0;
            return false;
        }
    }

    for (int j = 0; j < updated_data_vertices_count_; ++j) {
        int candidate = updated_data_vertices_flag_[j];
        data_vertices_flag_[candidate] = 0;
    }
    updated_data_vertices_count_ = 0;

    return true;
}

bool FSC::Verification(const Graph* query_graph, const Graph *data_graph) {
    timeval ordering_start = Utility::GetTime();
    ComputeDirectedWeightedGraph(query_graph, data_graph);
    GenerateVerificationOrder(query_graph);
    GenerateBigraphIndex(data_graph);
    timeval ordering_end = Utility::GetTime();

    temp_ordering_time_ = Utility::TimeDiffInMicroseconds(ordering_start, ordering_end);

    int cur_level = 0;
    candidates_index_count_[0] = candidates_count_[verification_start_vertex_];

    for (int i = 0; i < candidates_index_count_[0]; ++i) {
        int u = sequence_[0].id_;
        int v = Candidate(u, i);

        candidates_index_matching_[u] = i;
        matching_[u] = v;
        embedding_[0] = v;
        data_vertex_visited_[v] = true;

        cur_level += 1;
        idx_[cur_level] = 0;
        int pivot = sequence_[cur_level].pivot_;
        int pivot_index = candidates_index_matching_[pivot];
        candidates_index_[cur_level] = CandidateIndex(sequence_[cur_level].id_, pivot_index,
                                                      candidates_index_count_[cur_level]);

        while (true) {
            while (idx_[cur_level] < candidates_index_count_[cur_level]) {
                int candidate_idx = candidates_index_[cur_level][idx_[cur_level]];
                u = sequence_[cur_level].id_;
                v = Candidate(u, candidate_idx);

                idx_[cur_level] += 1;

                if (IsValid(cur_level, v, data_graph)) {
                    candidates_index_matching_[u] = candidate_idx;
                    matching_[u] = v;
                    embedding_[cur_level] = v;
                    data_vertex_visited_[v] = true;
                    if (cur_level + 1 == simulation_order_count_) {
                        for (int j = 0; j < simulation_order_count_; ++j) {
                            data_vertex_visited_[embedding_[j]] = false;
                        }
                        return true;
                    } else {
                        cur_level += 1;
                        idx_[cur_level] = 0;
                        pivot = sequence_[cur_level].pivot_;
                        pivot_index = candidates_index_matching_[pivot];
                        candidates_index_[cur_level] = CandidateIndex(sequence_[cur_level].id_, pivot_index,
                                                                      candidates_index_count_[cur_level]);
                    }
                }
            }

            cur_level -= 1;
            data_vertex_visited_[embedding_[cur_level]] = false;

            if (cur_level <= 0) {
                break;
            }
        }
    }
    return false;
}

void FSC::UpdateSequenceElement(const Graph* query_graph, const int pivot, const int u, SequenceElement &element) {
    element.id_ = u;
    element.pivot_ = pivot;
    element.bn_count_ = 0;

    int count;
    const int* neighbors = query_graph->Neighbors(u, count);

    for (int i = 0; i < count; ++i) {
        int n_u = neighbors[i];
        if (query_vertices_flags_[n_u]) {
            element.bn_[element.bn_count_++] = n_u;
        }
    }
}

void FSC::GenerateNLF(const Graph* query_graph) {
    for (int i = 0; i < query_graph->VerticesCount(); ++i) {
        TreeNode& node = tree_node_[i];
        node.id_ = i;
        node.label_ = query_graph->Label(i);
        node.degree_ = query_graph->Degree(i);

        int count;
        const int* neighbors = query_graph->Neighbors(i, count);
        node.nlf_label_count_ = 0;

        for (int j = 0; j < count; ++j) {
            int label = query_graph->Label(neighbors[j]);
            bool new_label = true;

            for (int k = 0; k < node.nlf_label_count_; ++k) {
                if (label == node.nlf_[k].first) {
                    new_label = false;
                    node.nlf_[k].second += 1;
                    break;
                }
            }

            if (new_label) {
                node.nlf_[node.nlf_label_count_++] = make_pair(label, 1);
            }
        }
    }
}

void FSC::GetMemoryCost() {
    size_t total_memory_cost = (max_data_vertices_num_ * 4) + (max_data_vertices_with_same_label_ * simulation_order_count_) + (max_data_vertices_num_ * simulation_order_count_) + (max_degree_ * simulation_order_count_ * 2);
    double memory_cost_mb = (total_memory_cost * 4) / (1024.0 * 1024.0);
    if (peak_memory_cost_ < memory_cost_mb)
        peak_memory_cost_ = memory_cost_mb;
    printf("Memory Cost: %.4f MB.\n", memory_cost_mb);
}

bool FSC::IsValid(const int level, const int v, const Graph* data_graph) {
    SequenceElement &element = sequence_[level];

    if (!data_vertex_visited_[v]) {
        for (int j = 0; j < element.bn_count_; ++j) {
            int neighbor = element.bn_[j];
            int temp = matching_[neighbor];
            if (!data_graph->IsEdge(v, temp)) {
                return false;
            }
        }

        return true;
    }

    return false;
}

void FSC::GenerateBFSTree(const Graph *query_graph, const Graph *data_graph) {
    // Select the root vertex
    SelectRootVertex(query_graph, data_graph);

    // BFS start from the root vertex.
    queue<int> bfs_queue;
    int count = 0;

    bfs_queue.push(root_vertex_);
    query_vertices_flags_[root_vertex_] = true;

    tree_node_[root_vertex_].level_ = 0;

    while(!bfs_queue.empty()) {
        const int cur_vertex = bfs_queue.front();
        bfs_queue.pop();
        bfs_order_index_[cur_vertex] = count;
        bfs_order_[count++] = cur_vertex;

        tree_node_[cur_vertex].children_count_ = 0;

        int n_count;
        const int* n = query_graph->Neighbors(cur_vertex, n_count);
        for (int i = 0; i < n_count; ++i) {
            const int neighbor = n[i];

            if (!query_vertices_flags_[neighbor]) {
                bfs_queue.push(neighbor);
                query_vertices_flags_[neighbor] = true;

                tree_node_[neighbor].parent_ = cur_vertex;
                tree_node_[neighbor].level_ = tree_node_[cur_vertex] .level_ + 1;

                tree_node_[cur_vertex].children_[tree_node_[cur_vertex].children_count_++] = neighbor;
            }
        }
    }

    level_count_ = -1;
    for (int i = 0; i < query_graph->VerticesCount(); ++i) {
        int cur_vertex = bfs_order_[i];
        tree_node_[cur_vertex].under_level_count_ = 0;
        tree_node_[cur_vertex].bn_count_ = 0;
        tree_node_[cur_vertex].fn_count_ = 0;

        if (tree_node_[cur_vertex].level_ != level_count_) {
            level_count_ += 1;
            level_offset_[level_count_] = 0;
        }

        level_offset_[level_count_] += 1;

        int n_count;
        const int* n = query_graph->Neighbors(cur_vertex, n_count);
        for (int j = 0; j < n_count; ++j) {
            int neighbor = n[j];

            if (tree_node_[cur_vertex].level_ == tree_node_[neighbor].level_) {
                if (bfs_order_index_[neighbor] < bfs_order_index_[cur_vertex]) {
                    tree_node_[cur_vertex].bn_[tree_node_[cur_vertex].bn_count_++] = neighbor;
                }
                else {
                    tree_node_[cur_vertex].fn_[tree_node_[cur_vertex].fn_count_++] = neighbor;
                }
            }
            else if (tree_node_[cur_vertex].level_ > tree_node_[neighbor].level_) {
                tree_node_[cur_vertex].bn_[tree_node_[cur_vertex].bn_count_++] = neighbor;
            }
            else {
                tree_node_[cur_vertex].under_level_[tree_node_[cur_vertex].under_level_count_++] = neighbor;
            }
        }
    }

    level_count_ += 1;

    int prev_value = 0;
    for (int i = 1; i <= level_count_; ++i) {
        int temp = level_offset_[i];
        level_offset_[i] = level_offset_[i - 1] + prev_value;
        prev_value = temp;
    }
    level_offset_[0] = 0;

    GenerateNLF(query_graph);
}

void FSC::SelectRootVertex(const Graph *query_graph, const Graph *data_graph) {
    double min_score = numeric_limits<double>::max();
    for (int i = 0; i < simulation_order_count_; ++i) {
        if ((query_vertices_cores_[i] > 1 || core_length_ == 0) && query_graph->Degree(i) >= 1) {
            const int label = query_graph->Label(i);
            const int degree = query_graph->Degree(i);
            double cur_score = data_graph->VerticesCount(label) / (double)degree;
             if (cur_score < min_score) {
                 min_score = cur_score;
                 root_vertex_ = i;
             }
        }
    }
}

bool FSC::GenerateRootCandidate(const Graph *data_graph) {
    TreeNode& node = tree_node_[root_vertex_];
    int candidate_count;
    const int * candidates = data_graph->VerticesWithLabel(node.label_, candidate_count);

    for (int i = 0; i < candidate_count; ++i) {
        int v = candidates[i];

        if (data_graph->Degree(v) >= node.degree_) {
#ifdef ENABLE_NEIGHBOR_LABEL_FREQUENCY_FILTER
            bool is_valid = true;
            const unordered_map<int, int>* data_vertex_nlf = data_graph->NLF(v);

            if (data_vertex_nlf->size() >= (size_t)node.nlf_label_count_) {
                for (int k = 0; k < node.nlf_label_count_; ++k) {
                    pair<int, int> label_count = node.nlf_[k];
                    auto data_vertex_iter = data_vertex_nlf->find(label_count.first);
                    if (data_vertex_iter == data_vertex_nlf->end() ||
                        data_vertex_iter->second < label_count.second) {
                        is_valid = false;
                        break;
                    }
                }
            }
            else {
                is_valid = false;
            }

            if (is_valid) {
                candidates_[node.id_][candidates_count_[node.id_]++] = v;
            }
#endif
#ifndef ENABLE_NEIGHBOR_LABEL_FREQUENCY_FILTER
            candidates_[node.id_][candidates_count_[node.id_]++] = v;
#endif
        }
    }

    return candidates_count_[node.id_] != 0;
}

void FSC::GenerateBigraphIndex(const Graph *data_graph) {
    for (int i = 1; i < simulation_order_count_; ++i) {
        // First compact the vertex.
        SequenceElement& element = sequence_[i];

        // Update the flag array.
        for (int j = 0; j < candidates_count_[element.id_]; ++j) {
            int query_vertex_candidate = candidates_[element.id_][j];
            data_vertices_flag_[query_vertex_candidate] = j + 1;
            updated_data_vertices_flag_[updated_data_vertices_count_++] = query_vertex_candidate;
        }

        int pivot = element.pivot_;

        // Generate the CSR graph of bigraph.
        int edge_count = 0;

        for (int j = 0; j < candidates_count_[pivot]; ++j) {
            int backward_neighbor_candidate = candidates_[pivot][j];
            int neighbor_count;
            const int *neighbors = data_graph->Neighbors(backward_neighbor_candidate, neighbor_count);

            bigraph_offset_[element.id_][j] = edge_count;
            // Loop each neighbors of candidate
            for (int l = 0; l < neighbor_count; ++l) {
                int neighbor = neighbors[l];
                if (data_vertices_flag_[neighbor] != 0) {
                    bigraph_edges_[element.id_][edge_count++] = data_vertices_flag_[neighbor] - 1;
                }
            }
        }
        bigraph_offset_[element.id_][candidates_count_[pivot]] = edge_count;

        // Clear flag_ array
        for (int j = 0; j < updated_data_vertices_count_; ++j) {
            int candidate = updated_data_vertices_flag_[j];
            data_vertices_flag_[candidate] = 0;
        }

        updated_data_vertices_count_ = 0;
    }
}


void FSC::CorePath(const int depth, const int cur_vertex, vector<pair<int, size_t >> &path) {
    bool core_leaf = true;
    path[depth].first = cur_vertex;
    TreeNode& cur_node = tree_node_[cur_vertex];

    for (int i = 0; i < cur_node.children_count_; ++i) {
        const int child = cur_node.children_[i];
        if (query_vertices_cores_[child] > 1) {
            CorePath(depth + 1, child, path);
            core_leaf = false;
        }
    }

    if (core_leaf) {
        // Add core path
        const int path_length = depth + 1;
        core_paths_.emplace_back(path_length);
        copy(path.begin(), path.begin() + path_length, core_paths_.back().begin());
    }
}

void FSC::TreePath(const int depth, const int cur_vertex, vector<pair<int, size_t>> &path) {
    bool tree_leaf = query_vertices_cores_[cur_vertex] == 1;
    path[depth].first = cur_vertex;
    TreeNode& cur_node = tree_node_[cur_vertex];

    for (int i = 0; i < cur_node.children_count_; ++i) {
        const int child = cur_node.children_[i];
        if (tree_node_[child].degree_ > 1) {
            TreePath(depth + 1, child, path);
            tree_leaf = false;
        }
    }

    if (tree_leaf) {
        const int tree_root = path[0].first;
        auto paths = tree_paths_.find(tree_root);
        const int path_length = depth + 1;

        if (paths != tree_paths_.end()) {
            paths->second.emplace_back(path_length);
            copy(path.begin(), path.begin() + path_length, paths->second.back().begin());
        }
        else {
            vector<vector<pair<int, size_t >>> tree_paths;
            tree_paths.emplace_back(path_length);
            copy(path.begin(), path.begin() + path_length, tree_paths.back().begin());
            tree_paths_.emplace(tree_root, tree_paths);
        }
    }
}

void FSC::GenerateVerificationOrder(const Graph *query_graph) {
    memset(query_vertices_flags_, 0, simulation_order_count_ * sizeof(bool));
    memset(backward_value_, 0, simulation_order_count_ * sizeof(int));
    for (int i = 0; i < simulation_order_count_; ++i)
        weight_[i] = numeric_limits<double>::max();
    SelectStartVertex(query_graph);
    // start_vertex_ = 10;
    UpdateSequenceElement(query_graph, 0, verification_start_vertex_, sequence_[0]);
    UpdateBackwardValue(query_graph, verification_start_vertex_);
    query_vertices_flags_[verification_start_vertex_] = true;

    for (int i = 1; i < simulation_order_count_; ++i) {
        int pivot;
        int next_vertex;
        SelectNextVertex(query_graph, i, next_vertex, pivot);
        UpdateSequenceElement(query_graph, pivot, next_vertex, sequence_[i]);
        UpdateBackwardValue(query_graph, next_vertex);
        query_vertices_flags_[next_vertex] = true;
    }
}

void FSC::SelectNextVertex(const Graph *query_graph, const int idx, int &next_query_vertex, int &pivot) {
    /** Select the next matching vertex.
      * A. Core structure.
      *  1. backward_value
      *  2. core number
      *  3. weight
      *  4. label id (TODO).
      * B. Tree structure.
      *  1. backward_value
      *  2. weight
      *  3. degree
      *  4. label id (TODO).
      * C. Leaf structure.
      *  1. backward_value
      *  2. weight
      *  3. label id (TODO).
      * */

    target_size_ = 0;
    result_size_ = 0;

    int max_backward_value = -1;
    for (int i = 0; i < query_graph->VerticesCount(); ++i) {
        if (!query_vertices_flags_[i] && backward_value_[i] > 0) {
            if (idx < core_length_ && query_vertices_cores_[i] < 2)
                continue;

            if (idx < query_graph->PrunedVerticesCount() && query_graph->Degree(i) == 1)
                continue;

            if (backward_value_[i] > max_backward_value) {
                max_backward_value = backward_value_[i];
                result_size_ = 0;
                result_[result_size_++] = i;
            } else if (backward_value_[i] == max_backward_value) {
                result_[result_size_++] = i;
            }
        }
    }

    if (result_size_ == 1) {
        next_query_vertex = result_[0];
        pivot = pivot_[next_query_vertex];
        // cout << next_query_vertex << ' ' << selected_backward_value << ' ' << max_backward_value << endl;
        return;
    }

    if (max_backward_value == 1) {
        int* temp_pointer = target_;
        target_ = result_;
        result_ = temp_pointer;

        target_size_ = result_size_;
        result_size_ = 0;

        int max_core_number = -1;
        for (int i = 0; i < target_size_; ++i) {
            int u = target_[i];

            if (query_vertices_cores_[u] > max_core_number) {
                max_core_number = query_vertices_cores_[u];
                result_size_ = 0;
                result_[result_size_++] = u;
            } else if (query_vertices_cores_[u] == max_core_number) {
                result_[result_size_++] = u;
            }
        }

        if (result_size_ == 1) {
            next_query_vertex = result_[0];
            pivot = pivot_[next_query_vertex];
            // cout << next_query_vertex << ' ' << selected_backward_value << ' ' << max_backward_value << endl;
            return;
        }

        temp_pointer = target_;
        target_ = result_;
        result_ = temp_pointer;

        target_size_ = result_size_;
        result_size_ = 0;

        double min_weight = numeric_limits<double>::max();

        for (int i = 0; i < target_size_; ++i) {
            int u = target_[i];

            if (weight_[u] < min_weight) {
                min_weight = weight_[u];
                result_size_ = 0;
                result_[result_size_++] = u;
            } else if (weight_[u] == min_weight) {
                result_[result_size_++] = u;
            }
        }

        if (result_size_ == 1) {
            next_query_vertex = result_[0];
            pivot = pivot_[next_query_vertex];
            // cout << next_query_vertex << ' ' << min_weight << endl;
            return;
        }

        temp_pointer = target_;
        target_ = result_;
        result_ = temp_pointer;

        target_size_ = result_size_;
        result_size_ = 0;

        int max_degree = -1;
        for (int i = 0; i < target_size_; ++i) {
            int u = target_[i];

            if (query_graph->Degree(u) > max_degree) {
                max_degree = query_graph->Degree(u);
                result_size_ = 0;
                result_[result_size_++] = u;
            } else if (query_graph->Degree(u) == max_degree) {
                result_[result_size_++] = u;
            }
        }

        next_query_vertex = result_[0];
        pivot = pivot_[next_query_vertex];
    }
    else {
        int * temp_pointer = target_;
        target_ = result_;
        result_ = temp_pointer;

        target_size_ = result_size_;
        result_size_ = 0;

        double min_weight = numeric_limits<double>::max();

        for (int i = 0; i < target_size_; ++i) {
            int u = target_[i];
            if (weight_[u] < min_weight) {
                min_weight = weight_[u];
                result_size_ = 0;
                result_[result_size_++] = u;
            } else if (weight_[u] == min_weight) {
                result_[result_size_++] = u;
            }
        }

        if (result_size_ == 1) {
            next_query_vertex = result_[0];
            pivot = pivot_[next_query_vertex];
            // cout << next_query_vertex << ' ' << min_weight << endl;
            return;
        }

        temp_pointer = target_;
        target_ = result_;
        result_ = temp_pointer;

        target_size_ = result_size_;
        result_size_ = 0;

        int max_core_number = -1;
        for (int i = 0; i < target_size_; ++i) {
            int u = target_[i];

            if (query_vertices_cores_[u] > max_core_number) {
                max_core_number = query_vertices_cores_[u];
                result_size_ = 0;
                result_[result_size_++] = u;
            } else if (query_vertices_cores_[u] == max_core_number) {
                result_[result_size_++] = u;
            }
        }

        if (result_size_ == 1) {
            next_query_vertex = result_[0];
            pivot = pivot_[next_query_vertex];
            // cout << next_query_vertex << ' ' << selected_backward_value << ' ' << max_backward_value << endl;
            return;
        }

        temp_pointer = target_;
        target_ = result_;
        result_ = temp_pointer;

        target_size_ = result_size_;
        result_size_ = 0;

        int max_degree = -1;
        for (int i = 0; i < target_size_; ++i) {
            int u = target_[i];

            if (query_graph->Degree(u) > max_degree) {
                max_degree = query_graph->Degree(u);
                result_size_ = 0;
                result_[result_size_++] = u;

            } else if (query_graph->Degree(u) == max_degree) {
                result_[result_size_++] = u;
            }
        }

        next_query_vertex = result_[0];
        pivot = pivot_[next_query_vertex];
    }
}

void FSC::UpdateBackwardValue(const Graph *query_graph, const int u) {
    int count;
    const int* neighbors = query_graph->Neighbors(u, count);

    for (int i = 0; i < count; ++i) {
        int v = neighbors[i];
        if (!query_vertices_flags_[v]) {
            backward_value_[v] += 1;

            double new_weight = directed_weighted_graph_[u * simulation_order_count_ + v];
            if (new_weight < weight_[v]) {
                weight_[v] = new_weight;
                pivot_[v] = u;
            }
        }
    }
}

void FSC::SelectStartVertex(const Graph *query_graph) {
    /** Select the start vertex of the matching order. The priority is as follows:
    * A. If |V(q)| > 2, do not consider the leaf vertex.
    *  1. core value;
    *  2. the number of candidates;
    *  3. degree;
    *  4. label id.
    * B. If |V(q)| == 2.
    *  1. the number of candidates;
    *  2. label id
    */

    target_size_ = 0;
    result_size_ = 0;

    if (simulation_order_count_ > 2) {

        int max_core_number = -1;
        for (int i = 0; i < simulation_order_count_; ++i) {
            // Filter the leaf vertices.
            if (query_graph->Degree(i) > 1 ) {
                if (query_vertices_cores_[i] > max_core_number) {
                    max_core_number = query_vertices_cores_[i];
                    result_size_ = 0;
                    result_[result_size_++] = i;
                } else if (query_vertices_cores_[i] == max_core_number) {
                    result_[result_size_++] = i;
                }
            }
        }

        if (result_size_ == 1) {
            verification_start_vertex_ = result_[0];
            return;
        }

        int* temp_pointer = target_;
        target_ = result_;
        result_ = temp_pointer;

        target_size_ = result_size_;
        result_size_ = 0;

        int min_candidate_count = numeric_limits<int>::max();
        for (int i = 0; i < target_size_; ++i) {
            int u = target_[i];

            if (candidates_count_[u] < min_candidate_count) {
                min_candidate_count = candidates_count_[u];
                result_size_ = 0;
                result_[result_size_++] = u;
            } else if (candidates_count_[u] == min_candidate_count) {
                result_[result_size_++] = u;
            }
        }

        if (result_size_ == 1) {
            verification_start_vertex_ = result_[0];
            return;
        }

        temp_pointer = target_;
        target_ = result_;
        result_ = temp_pointer;

        target_size_ = result_size_;
        result_size_ = 0;

        int max_degree = -1;
        for (int i = 0; i < target_size_; ++i) {
            int u = target_[i];

            if (query_graph->Degree(u) > max_degree) {
                max_degree = query_graph->Degree(u);
                result_size_ = 0;
                result_[result_size_++] = u;
            } else if (query_graph->Degree(u) == max_degree) {
                result_[result_size_++] = u;
            }
        }

        verification_start_vertex_ = result_[0];
    }
    else {
        if (candidates_count_[0] > candidates_count_[1])
            verification_start_vertex_ = 1;
        else
            verification_start_vertex_ = 0;

    }
}

void FSC::ComputeDirectedWeightedGraph(const Graph *query_graph, const Graph* data_graph) {
    for (int i = 0; i < simulation_order_count_; ++i) {
        int query_vertex = i;

        int count;
        const int* neighbors = query_graph->PrunedNeighbors(query_vertex, count);

        bool generated = false;
        for (int j = 0; j < count; ++j) {
            int neighbor = neighbors[j];
            int edge_count = 0;
            if (directed_weighted_graph_[query_vertex * simulation_order_count_ + neighbor] == 0) {

                if (!generated) {
                    generated = true;

                    for (int k = 0; k < candidates_count_[query_vertex]; ++k) {
                        int candidate = candidates_[query_vertex][k];
                        data_vertices_flag_[candidate] = 1;
                        updated_data_vertices_flag_[updated_data_vertices_count_++] = candidate;
                    }
                }

                for (int k = 0; k < candidates_count_[neighbor]; ++k) {
                    int candidate = candidates_[neighbor][k];

                    int data_vertex_neighbor_count;
                    const int *data_vertex_neighbors = data_graph->Neighbors(candidate, data_vertex_neighbor_count);

                    // Loop each neighbors of candidate
                    for (int l = 0; l < data_vertex_neighbor_count; ++l) {
                        int data_vertex_neighbor = data_vertex_neighbors[l];
                        if (data_vertices_flag_[data_vertex_neighbor] != 0) {
                            edge_count += 1;
                        }
                    }
                }

                directed_weighted_graph_[query_vertex * simulation_order_count_ + neighbor] = edge_count / (double)candidates_count_[query_vertex];
                directed_weighted_graph_[neighbor * simulation_order_count_ + query_vertex] = edge_count / (double)candidates_count_[neighbor];
            }
        }

        for (int j = 0; j < updated_data_vertices_count_; ++j) {
            int candidate = updated_data_vertices_flag_[j];
            data_vertices_flag_[candidate] = 0;
        }
        updated_data_vertices_count_ = 0;
    }
}