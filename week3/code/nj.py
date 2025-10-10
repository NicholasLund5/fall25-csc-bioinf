# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

__author__ = "Patrick Kunzmann"
__all__ = ["neighbor_joining"]

import numpy as np
from tree import Tree, TreeNode


MAX_FLOAT = np.finfo(np.float32).max


def neighbor_joining(distances):
    distances = np.asarray(distances)
    
    if distances.shape[0] != distances.shape[1] or not np.allclose(distances.T, distances):
        raise ValueError("Distance matrix must be symmetric")
    if np.isnan(distances).any():
        raise ValueError("Distance matrix contains NaN values")
    if (distances >= MAX_FLOAT).any():
        raise ValueError("Distance matrix contains infinity")
    if distances.shape[0] < 4:
        raise ValueError("At least 4 nodes are required")
    if (distances < 0).any():
        raise ValueError("Distances must be positive")

    nodes = np.array([TreeNode(index=i) for i in range(distances.shape[0])])
    is_clustered = np.full(distances.shape[0], False, dtype=np.uint8)
    n_rem_nodes = len(distances) - np.count_nonzero(is_clustered)
    divergence = np.zeros(distances.shape[0], dtype=np.float32)
    corr_distances = np.zeros((distances.shape[0],) * 2, dtype=np.float32)
    distances_mat = distances.astype(np.float32, copy=True)

    while True:
        # Calculate divergence
        for i in range(distances_mat.shape[0]):
            if is_clustered[i]:
                continue
            dist_sum = 0
            for k in range(distances_mat.shape[0]):
                if is_clustered[k]:
                    continue
                dist_sum += distances_mat[i, k]
            divergence[i] = dist_sum
        
        # Calculate corrected distance matrix
        for i in range(distances_mat.shape[0]):
            if is_clustered[i]:
                continue
            for j in range(i):
                if is_clustered[j]:
                    continue
                corr_distances[i, j] = (n_rem_nodes - 2) * distances_mat[i, j] - divergence[i] - divergence[j]

        # Find minimum corrected distance
        dist_min = MAX_FLOAT
        i_min = -1
        j_min = -1
        for i in range(corr_distances.shape[0]):
            if is_clustered[i]:
                continue
            for j in range(i):
                if is_clustered[j]:
                    continue
                dist = corr_distances[i, j]
                if dist < dist_min:
                    dist_min = dist
                    i_min = i
                    j_min = j
        
        if i_min == -1 or j_min == -1:
            break
        
        node_dist_i = 0.5 * (distances_mat[i_min, j_min] + 1/(n_rem_nodes-2) * (divergence[i_min] - divergence[j_min]))
        node_dist_j = 0.5 * (distances_mat[i_min, j_min] + 1/(n_rem_nodes-2) * (divergence[j_min] - divergence[i_min]))
        
        if n_rem_nodes > 3:
            nodes[i_min] = TreeNode((nodes[i_min], nodes[j_min]), (node_dist_i, node_dist_j))
            nodes[j_min] = None
            is_clustered[j_min] = True
        else:
            is_clustered[i_min] = True
            is_clustered[j_min] = True
            k = np.where(~is_clustered.astype(bool))[0][0]
            node_dist_k = 0.5 * (distances_mat[i_min, k] + distances_mat[j_min, k] - distances_mat[i_min, j_min])
            root = TreeNode((nodes[i_min], nodes[j_min], nodes[k]), (node_dist_i, node_dist_j, node_dist_k))
            return Tree(root)
        
        # Update distance matrix
        for k in range(distances_mat.shape[0]):
            if not is_clustered[k] and k != i_min:
                dist = 0.5 * (distances_mat[i_min, k] + distances_mat[j_min, k] - distances_mat[i_min, j_min])
                distances_mat[i_min, k] = dist
                distances_mat[k, i_min] = dist

        n_rem_nodes = len(distances) - np.count_nonzero(is_clustered)