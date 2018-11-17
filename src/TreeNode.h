//
// Created by ssunah on 5/12/18.
//

#ifndef FSC_TREENODE_H
#define FSC_TREENODE_H

#include "Config.h"
class TreeNode {
public:
    int id_;
    int parent_;
    int label_;
    int degree_;
    int level_;
    int under_level_count_;
    int children_count_;
    int bn_count_;
    int fn_count_;
    int nlf_label_count_;
    int* under_level_;
    int* children_;
    int* bn_;
    int* fn_;
    pair<int, int>* nlf_;

public:
    TreeNode() {
        under_level_ = new int[MAX_QUERY_GRAPH_VERTICES_COUNT];
        bn_ = new int[MAX_QUERY_GRAPH_VERTICES_COUNT];
        fn_ = new int[MAX_QUERY_GRAPH_VERTICES_COUNT];
        nlf_ = new pair<int, int>[MAX_QUERY_GRAPH_VERTICES_COUNT];
        children_ = new int[MAX_QUERY_GRAPH_VERTICES_COUNT];
    }

    ~TreeNode() {
        delete[] under_level_;
        delete[] bn_;
        delete[] fn_;
        delete[] nlf_;
        delete[] children_;
    }
};

#endif //FSC_TREENODE_H
