#include "Graph.h"
#include <iostream>

void Graph::BuildReverseLabelIndex(const int label_num) {
    reverse_label_offset_ = new int[label_num + 1];
    reverse_label_vertices_ = new int[vertices_num_];
    labels_num_ = 0;
    memset(reverse_label_offset_, 0, sizeof(int) * (label_num + 1));

    // Histogram
    for (int i = 0; i < vertices_num_; ++i) {
        int label = labels_[i];
        if (reverse_label_offset_[label] == 0)
            labels_num_ += 1;
        reverse_label_offset_[label] += 1;
    }

    // Prefix Sum
    int temp = reverse_label_offset_[0];
    reverse_label_offset_[0] = 0;
    for (int i = 1; i <= label_num; ++i) {
        int cur = reverse_label_offset_[i];
        reverse_label_offset_[i] = reverse_label_offset_[i - 1] + temp;
        temp = cur;
    }

    // Set index
    for (int i = 0; i < vertices_num_; ++i) {
        int label = labels_[i];
        reverse_label_vertices_[reverse_label_offset_[label]++] = i;
    }

    // Recalculate the offset
    for (int i = label_num; i >= 1; --i) {
        reverse_label_offset_[i] = reverse_label_offset_[i - 1];
    }

    reverse_label_offset_[0] = 0;

    max_vertices_num_with_same_label_ = -1;
    for (int i = 0; i < label_num; ++i) {
        if (reverse_label_offset_[i + 1] - reverse_label_offset_[i] > max_vertices_num_with_same_label_)
            max_vertices_num_with_same_label_ = reverse_label_offset_[i + 1] - reverse_label_offset_[i];
    }
}

void Graph::ComputeLabelFrequency() {
    map<int, int> label_frequency;

    for (int i = 0; i < vertices_num_; ++i) {
        int label = labels_[i];

        if (label_frequency.find(label) == label_frequency.end())
            label_frequency[label] = 0;
        label_frequency[label] += 1;
    }

    label_frequency_ = new pair<int, int>[label_frequency.size()];
    labels_num_ = label_frequency.size();

    int i = 0;
    for (map<int, int>::iterator iter = label_frequency.begin(); iter != label_frequency.end(); iter++) {
        label_frequency_[i++] = make_pair(iter->first, iter->second);
    }
}

void Graph::SortEdge() {
    for (int i = 0; i < vertices_num_; ++i) {
        sort(edges_ + edges_offset_[i], edges_ + edges_offset_[i + 1]);
    }
}

void Graph::SortEdgeByLabel() {
    edges_sort_by_label_ = new int[edges_num_ * 2];
    copy(edges_, edges_ + edges_num_ * 2, edges_sort_by_label_);

    for (int i = 0; i < vertices_num_; ++i) {
        sort(edges_sort_by_label_ + edges_offset_[i], edges_sort_by_label_ + edges_offset_[i + 1],
            [this](const int& x, const int& y) -> bool {
                if (labels_[x] < labels_[y])
                {
                    return true;
                }
                else if (labels_[x] == labels_[y]) {
                    if (x < y) {
                        return true;
                    }
                }
                return false;
            });
    }
}

void Graph::BuildNLF() {
    nlf_ = new unordered_map<int, int>[vertices_num_];
    for (int i = 0; i < vertices_num_; ++i) {
        int count;
        const int* neighbors = Neighbors(i, count);
        for (int j = 0; j < count; ++j) {
            int neighbor = neighbors[j];
            int label = Label(neighbor);

            if (nlf_[i].find(label) == nlf_[i].end()) {
                nlf_[i][label] = 1;
            }
            else {
                nlf_[i][label] += 1;
            }
        }
    }
}

void Graph::ComputeMaxDegree() {
    max_degree_ = -1;
    for (int i = 0; i < vertices_num_; ++i) {
        if (Degree(i) > max_degree_)
            max_degree_ = Degree(i);
    }
}

void Graph::BuildPrunedGraph() {
    pruned_edges_offset_ = new int[vertices_num_ + 1];
    pruned_edges_ = new int[edges_num_ * 2];

    pruned_edges_offset_[0] = 0;

    pruned_vertices_num_ = 0;

    int pruned_edges_num = 0;

    map<int, int> label_frequency;

    for (int i = 0; i < vertices_num_; ++i) {
        int degree = Degree(i);
        if (degree > 1) {
            int count;
            const int* neighbors = SortedNeighbors(i, count);

            for (int j = 0; j < count; ++j) {
                int neighbor = neighbors[j];
                int neighbor_degree = Degree(neighbor);
                if (neighbor_degree > 1) {
                    pruned_edges_[pruned_edges_num++] = neighbor;
                }
            }

            pruned_vertices_num_ += 1;

            int label = labels_[i];

            if (label_frequency.find(label) == label_frequency.end())
                label_frequency[label] = 0;
            label_frequency[label] += 1;
        }

        pruned_edges_offset_[i + 1] = pruned_edges_num;
    }

    pruned_label_frequency_ = new pair<int, int>[label_frequency.size()];
    pruned_labels_num_ = label_frequency.size();

    int i = 0;
    for (map<int, int>::iterator iter = label_frequency.begin(); iter != label_frequency.end(); iter++) {
        pruned_label_frequency_[i++] = make_pair(iter->first, iter->second);
    }
}

void Graph::BuildDegreeOneVertices() {
    degree_one_vertices_num_ = 0;
    degree_one_labels_num_ = 0;
    degree_one_vertices_ = new int[vertices_num_];
    degree_one_offset_ = new int[vertices_num_ + 1];
    degree_one_labels_count_ = new int[vertices_num_];
    memset(degree_one_labels_count_, 0, sizeof(int) * vertices_num_);
    degree_one_parent_labels_count_ = new int[vertices_num_];
    memset(degree_one_parent_labels_count_, 0, sizeof(int) * vertices_num_);

    map<int, int> label_frequency;
    for (int i = 0; i < vertices_num_; ++i) {
        if (Degree(i) == 1) {
            degree_one_vertices_[degree_one_vertices_num_++] = i;
            int label = labels_[i];

            if (label_frequency.find(label) == label_frequency.end())
                label_frequency[label] = 0;

            label_frequency[label] += 1;
        }
        else {
            int count;
            const int* neighbors = Neighbors(i, count);

            map<int, int> parent_label_frequency;
            for (int j = 0; j < count; ++j) {
                int neighbor = neighbors[j];
                if (Degree(neighbor) == 1) {
                    int label = labels_[neighbor];
                    if (parent_label_frequency.find(label) == parent_label_frequency.end())
                        parent_label_frequency[label] = 0;
                    parent_label_frequency[label] += 1;
                }
            }

            for (int j = 0; j < count; ++j) {
                int neighbor = neighbors[j];

                if (Degree(neighbor) == 1) {
                    int label = labels_[neighbor];
                    degree_one_parent_labels_count_[neighbor] = parent_label_frequency[label];
                }
            }
        }
    }

    for (int i = 0; i < degree_one_vertices_num_; ++i) {
        int u = degree_one_vertices_[i];
        int label = labels_[u];
        degree_one_labels_count_[u] = label_frequency[label];
    }

    sort(degree_one_vertices_, degree_one_vertices_ + degree_one_vertices_num_, [this](const int& x, const int& y) -> bool {
        if (this->degree_one_labels_count_[x] < this->degree_one_labels_count_[y]) {
            return true;
        }
        else if (this->degree_one_labels_count_[x] == this->degree_one_labels_count_[y]) {
            if (this->labels_[x] < this->labels_[y]) {
                return true;
            }
            else if (this->labels_[x] == this->labels_[y]) {
                if (this->degree_one_parent_labels_count_[x] < this->degree_one_parent_labels_count_[y]) {
                    return true;
                }
                else if (this->degree_one_parent_labels_count_[x] == this->degree_one_parent_labels_count_[y]){
                    int x_parent = edges_[edges_offset_[x]];
                    int y_parent = edges_[edges_offset_[y]];

                    if (x_parent < y_parent) {
                        return true;
                    }
                }
            }
        }

        return false;
    });

    degree_one_labels_num_ = 0;
    int cur_label = -1;
    for (int i = 0; i < degree_one_vertices_num_; ++i) {
        int u = degree_one_vertices_[i];
        int label = labels_[u];

        if (cur_label != label) {
            cur_label = label;
            degree_one_offset_[degree_one_labels_num_++] = i;
        }
    }

    degree_one_offset_[degree_one_labels_num_] = degree_one_vertices_num_;
}