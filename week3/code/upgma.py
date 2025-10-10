# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

__author__ = "Patrick Kunzmann"
__all__ = ["upgma"]

import numpy as np
from tree import Tree, TreeNode


MAX_FLOAT = np.finfo(np.float32).max


def upgma(distances):
    distances = np.asarray(distances)
    
    if distances.shape[0] != distances.shape[1] \
        or not np.allclose(distances.T, distances):
            raise ValueError("Distance matrix must be symmetric")
    if np.isnan(distances).any():
        raise ValueError("Distance matrix contains NaN values")
    if (distances >= MAX_FLOAT).any():
        raise ValueError("Distance matrix contains infinity")
    if (distances < 0).any():
        raise ValueError("Distances must be positive")

    nodes = np.array(
        [TreeNode(index=i) for i in range(distances.shape[0])]
    )
    is_clustered = np.full(distances.shape[0], False, dtype=np.uint8)
    cluster_size = np.ones(distances.shape[0], dtype=np.uint32)
    node_heights = np.zeros(distances.shape[0], dtype=np.float32)

    distances_mat = distances.astype(np.float32, copy=True)
    
    while True:
        # Find minimum distance
        dist_min = MAX_FLOAT
        i_min = -1
        j_min = -1
        for i in range(distances_mat.shape[0]):
            if is_clustered[i]:
                continue
            for j in range(i):
                if is_clustered[j]:
                    continue
                dist = distances_mat[i, j]
                if dist < dist_min:
                    dist_min = dist
                    i_min = i
                    j_min = j
        
        if i_min == -1 or j_min == -1:
            break
        
        height = dist_min / 2
        nodes[i_min] = TreeNode(
            (nodes[i_min], nodes[j_min]),
            (height - node_heights[i_min], height - node_heights[j_min])
        )
        node_heights[i_min] = height
        nodes[j_min] = None
        is_clustered[j_min] = True
        
        # Calculate arithmetic mean distances
        for k in range(distances_mat.shape[0]):
            if not is_clustered[k] and k != i_min:
                mean = (
                    (
                        distances_mat[i_min, k] * cluster_size[i_min]
                        + distances_mat[j_min, k] * cluster_size[j_min]
                    ) / (cluster_size[i_min] + cluster_size[j_min])
                )
                distances_mat[i_min, k] = mean
                distances_mat[k, i_min] = mean
        
        cluster_size[i_min] = cluster_size[i_min] + cluster_size[j_min]

    return Tree(nodes[len(nodes) - 1])