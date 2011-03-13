// Copyright 2010, Andreas Biegert

#include "cs.h"
#include "guide_tree.h"

namespace cs {

// Constructs a guide tree by UPGMA using given distance matrix
GuideTree BuildGuideTree(Matrix<float> dist, bool verbose) {
    const size_t num_seqs = dist.nrows();
    GuideTree tree;
    Vector<TreeNode> nodes(num_seqs);  // nodes[i] is node of i-th cluster
    Bitset disabled(num_seqs);     // indicates disabled indices

    if (verbose) puts("Building UPGMA guide tree ...");

    // Add leaf nodes to guide tree
    for (size_t k = 0; k < num_seqs; ++k)
        nodes[k] = add_vertex(TreeNodeProperties(), tree);

    // Using dist[k][l] find the closest 2 clusters, amalgamate, and update
    // dist[i][j]'s. Repeatedly do this, until there is only 1 cluster left.
    size_t num_clusters = num_seqs;
    while (num_clusters > 1) {
        // Find clusters k and l with maximal distance
        size_t k = 0, l = 0;
        float dmin = FLT_MAX;
        for (size_t i = 0; i < num_seqs - 1; ++i) {
            for (size_t j = i + 1; j < num_seqs; ++j) {
                if (!disabled.test(i) && !disabled.test(j) && dist[i][j] < dmin) {
                    dmin = dist[i][j];
                    k = i;
                    l = j;
                }
            }
        }

        // Create new node p with child nodes c[k] and c[l]
        size_t depth = 1 + MAX(tree[nodes[k]].depth, tree[nodes[l]].depth);
        double height = dmin / 2.0;
        size_t sk = tree[nodes[k]].size;
        size_t sl = tree[nodes[l]].size;
        double bk = height - tree[nodes[k]].height;
        double bl = height - tree[nodes[l]].height;
        TreeNode p = add_vertex(TreeNodeProperties(depth, sk + sl, dmin, height), tree);
        add_edge(p, nodes[k], TreeEdgeProperties(bk), tree);
        add_edge(p, nodes[l], TreeEdgeProperties(bl), tree);
        nodes[k] = p;
        disabled.set(l);
        num_clusters--;

        // Update distance matrix
        float size = sk + sl;
        for (size_t j = 0; j < num_seqs; ++j) {
            if (j != l)
                dist[k][j] = dist[j][k] = (sk * dist[k][j] + sl * dist[l][j]) / size;
        }
        for (size_t i = 0; i < num_seqs; ++i)
            dist[i][l] = dist[l][i] = -1;

        LOG(ERROR) << strprintf("merging nodes %zu and %zu into node %zu (distance = %.4f  height = %.4f  length[%zu] = %.4f  length[%zu] = %.4f)",
                                k, l, p, dmin, height, k, bk, l, bl);
        if (verbose) printf("  merging nodes %zu and %zu (distance = %.4f)\n", k, l, dmin);
    }

    return tree;
}


GuideTree BuildLeftistGuideTree(size_t num_seqs, double branch_len, bool verbose) {
    GuideTree tree;
    Vector<TreeNode> nodes(num_seqs);  // nodes[i] is node of i-th cluster
    Bitset disabled(num_seqs);         // indicates disabled indices

    if (verbose) puts("Building leftist guide tree ...");

    // Add leaf nodes to guide tree
    for (size_t k = 0; k < num_seqs; ++k)
        nodes[k] = add_vertex(TreeNodeProperties(), tree);

    TreeNode p, l = 0;
    for (size_t r = 1; r < num_seqs; ++r) {
        p = add_vertex(TreeNodeProperties(0, 0, 0, 0), tree);
        add_edge(p, l, TreeEdgeProperties(branch_len), tree);
        add_edge(p, r, TreeEdgeProperties(branch_len * r), tree);
        l = p;
    }

    // // Add node 8
    // p = add_vertex(TreeNodeProperties(0, 0, 0, 0), tree);
    // add_edge(p, 0, TreeEdgeProperties(branch_len), tree);
    // add_edge(p, 1, TreeEdgeProperties(branch_len), tree);

    // // Add node 9
    // p = add_vertex(TreeNodeProperties(0, 0, 0, 0), tree);
    // add_edge(p, 2, TreeEdgeProperties(branch_len), tree);
    // add_edge(p, 8, TreeEdgeProperties(branch_len), tree);

    // // Add node 10
    // p = add_vertex(TreeNodeProperties(0, 0, 0, 0), tree);
    // add_edge(p, 4, TreeEdgeProperties(branch_len), tree);
    // add_edge(p, 5, TreeEdgeProperties(branch_len), tree);

    // // Add node 11
    // p = add_vertex(TreeNodeProperties(0, 0, 0, 0), tree);
    // add_edge(p, 6, TreeEdgeProperties(branch_len), tree);
    // add_edge(p, 7, TreeEdgeProperties(branch_len), tree);

    // // Add node 12
    // p = add_vertex(TreeNodeProperties(0, 0, 0, 0), tree);
    // add_edge(p, 8, TreeEdgeProperties(branch_len), tree);
    // add_edge(p, 9, TreeEdgeProperties(branch_len), tree);

    // // Add node 13
    // p = add_vertex(TreeNodeProperties(0, 0, 0, 0), tree);
    // add_edge(p, 10, TreeEdgeProperties(branch_len), tree);
    // add_edge(p, 11, TreeEdgeProperties(branch_len), tree);

    // // Add node 14
    // p = add_vertex(TreeNodeProperties(0, 0, 0, 0), tree);
    // add_edge(p, 12, TreeEdgeProperties(branch_len), tree);
    // add_edge(p, 13, TreeEdgeProperties(branch_len), tree);

    return tree;
}


}  // namespace cs
