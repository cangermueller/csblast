// Copyright 2009, Andreas Biegert

#ifndef CS_GUIDE_TREE_H_
#define CS_GUIDE_TREE_H_

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/graph_utility.hpp>

using boost::adjacency_list;
using boost::vecS;
using boost::bidirectionalS;
using boost::graph_traits;
using boost::property_map;
using boost::vertex_index_t;
using boost::vertex_index;

namespace cs {

// POD class representing a tree node.
struct TreeNodeProperties {
    // Constructs an empty node
    TreeNodeProperties(size_t d = 0, size_t n = 1, double t = 0.0, double h = 0.0)
            : depth(d), size(n), dist(t), height(h) {}

    size_t depth;   // depth of leafs is 0
    size_t size;    // number of seqs in this cluster
    double dist;    // distance time of child nodes
    double height;  // decimal height of this node in tree
};

// POD class representing a tree edge.
struct TreeEdgeProperties {
    // Constructs an edge with given length
    TreeEdgeProperties(double l = 0.0) : length(l) {}

    double length;  // edge length
};

typedef adjacency_list<vecS, vecS, bidirectionalS, TreeNodeProperties, TreeEdgeProperties> GuideTree;
typedef graph_traits<GuideTree> GuideTreeTraits;
typedef GuideTreeTraits::vertex_descriptor TreeNode;
typedef property_map<GuideTree, vertex_index_t>::type TreeNodeIndex;
typedef GuideTreeTraits::adjacency_iterator TreeAdjacencyIter;

// Constructs a guide tree by UPGMA using given distance matrix
GuideTree BuildGuideTree(Matrix<float> dist, bool verbose = true);

// Constructs simply a leftist guide tree, no need for a distance matrix
GuideTree BuildLeftistGuideTree(size_t num_seqs, double branch_len = 0.0857, bool verbose = true);

}  // namespace cs

#endif  // CS_GUIDE_TREE_H_
