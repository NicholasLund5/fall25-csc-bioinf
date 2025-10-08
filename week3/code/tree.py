# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

__author__ = "Patrick Kunzmann, Tom David Müller"
__all__ = ["Tree", "TreeNode", "TreeError"]

import copy
import numpy as np


class InvalidFileError(Exception):
    pass


class Copyable:
    def copy(self):
        return self.__copy_create__()


class Tree(Copyable):
    def __init__(self, root):
        if root is None:
            raise TypeError("Root node cannot be None")
        root.as_root()
        self._root = root
        
        leaves_unsorted = self._root.get_leaves()
        leaf_count = len(leaves_unsorted)
        indices = np.array([leaf.index for leaf in leaves_unsorted])
        self._leaves = [None] * leaf_count
        
        for i in range(len(indices)):
            index = indices[i]
            if index >= leaf_count or index < 0:
                raise TreeError("The tree's indices are out of range")
            self._leaves[index] = leaves_unsorted[i]
    
    def __copy_create__(self):
        return Tree(self._root.copy())
    
    @property
    def root(self):
        return self._root
    
    @property
    def leaves(self):
        return copy.copy(self._leaves)

    def get_distance(self, index1, index2, topological=False):
        return self._leaves[index1].distance_to(
            self._leaves[index2], topological
        )
    
    @staticmethod
    def from_newick(newick, labels=None):
        newick = newick.strip()
        if len(newick) == 0:
            raise InvalidFileError("Newick string is empty")
        if newick[-1] == ";":
            newick = newick[:-1]
        root, distance = TreeNode.from_newick(newick, labels)
        return Tree(root)

    def __str__(self):
        return self.to_newick()
    
    def __len__(self):
        return len(self._leaves)
    
    def __eq__(self, item):
        if not isinstance(item, Tree):
            return False
        return self._root == item._root
    
    def __hash__(self):
        return hash(self._root)


class TreeNode:
    def __init__(self, children=None, distances=None, index=None):
        self._is_root = False
        self._distance = 0
        self._parent = None
        
        if index is None:
            if children is None or distances is None:
                raise TypeError(
                    "Either reference index (for terminal node) or "
                    "child nodes including the distance "
                    "(for intermediate node) must be set"
                )
            for item in children:
                if not isinstance(item, TreeNode):
                    raise TypeError(
                        f"Expected 'TreeNode', but got '{type(item).__name__}'"
                    )
            for item in distances:
                if not isinstance(item, (float, int, np.floating, np.integer)):
                    raise TypeError(
                        f"Expected 'float' or 'int', "
                        f"but got '{type(item).__name__}'"
                    )
            if len(children) == 0:
                raise TreeError(
                    "Intermediate nodes must at least contain one child node"
                )
            if len(children) != len(distances):
                raise ValueError(
                    "The number of children must equal the number of distances"
                )
            for i in range(len(children)):
                for j in range(len(children)):
                    if i != j and children[i] is children[j]:
                        raise TreeError(
                            "Two child nodes cannot be the same object"
                        )
            self._index = -1
            self._children = tuple(children)
            for child, distance in zip(children, distances):
                child._set_parent(self, distance)
        elif index < 0:
            raise ValueError("Index cannot be negative")
        else:
            if children is not None or distances is not None:
                raise TypeError(
                    "Reference index and child nodes are mutually exclusive"
                )
            self._index = index
            self._children = None
    
    def _set_parent(self, parent, distance):
        if parent is None:
            raise TypeError("Parent node cannot be None")
        if self._parent is not None or self._is_root:
            raise TreeError("Node already has a parent")
        self._parent = parent
        self._distance = distance
    
    def copy(self):
        if self.is_leaf():
            return TreeNode(index=self._index)
        else:
            distances = [child.distance for child in self._children]
            children_clones = [child.copy() for child in self._children]
            return TreeNode(children_clones, distances)

    @property
    def index(self):
        return None if self._index == -1 else self._index
    
    @property
    def children(self):
        return self._children
    
    @property
    def parent(self):
        return self._parent
    
    @property
    def distance(self):
        return None if self._parent is None else self._distance

    def is_leaf(self):
        return False if self._index == -1 else True
    
    def is_root(self):
        return bool(self._is_root)
    
    def as_root(self):
        if self._parent is not None:
            raise TreeError("Node has parent, cannot be a root node")
        self._is_root = True
    
    def distance_to(self, node, topological=False):
        distance = 0
        lca = self.lowest_common_ancestor(node)
        if lca is None:
            raise TreeError("The nodes do not have a common ancestor")
        
        current_node = self
        while current_node is not lca:
            if topological:
                distance += 1
            else:
                distance += current_node._distance
            current_node = current_node._parent
        
        current_node = node
        while current_node is not lca:
            if topological:
                distance += 1
            else:
                distance += current_node._distance
            current_node = current_node._parent
        
        return distance
    
    def lowest_common_ancestor(self, node):
        lca = None
        self_path = _create_path_to_root(self)
        other_path = _create_path_to_root(node)
        
        for i in range(-1, -min(len(self_path), len(other_path)) - 1, -1):
            if self_path[i] is other_path[i]:
                lca = self_path[i]
            else:
                break
        return lca

    def get_leaves(self):
        leaf_list = []
        _get_leaves(self, leaf_list)
        return leaf_list
    
    @staticmethod
    def from_newick(newick, labels=None):
        newick = "".join(newick.split())

        subnewick_start_i = -1
        subnewick_stop_i = -1
        
        for i in range(len(newick)):
            char = newick[i]
            if char == "(":
                subnewick_start_i = i
                break
            if char == ")":
                raise InvalidFileError("Bracket closed before it was opened")
        
        for i in reversed(range(len(newick))):
            char = newick[i]
            if char == ")":
                subnewick_stop_i = i + 1
                break
            if char == "(":
                raise InvalidFileError("Bracket was opened but not closed")
        
        if subnewick_start_i == -1 and subnewick_stop_i == -1:
            label_and_distance = newick
            try:
                label, distance = label_and_distance.split(":")
                distance = float(distance)
            except ValueError:
                distance = 0
                label = label_and_distance
            index = int(label) if labels is None else labels.index(label)
            return TreeNode(index=index), distance
        
        else:
            if subnewick_stop_i == len(newick):
                label = None
                distance = 0
            else:
                label_and_distance = newick[subnewick_stop_i:]
                try:
                    label, distance = label_and_distance.split(":")
                    distance = float(distance)
                except ValueError:
                    distance = 0
                    label = label_and_distance
                distance = float(distance)
            
            subnewick = newick[subnewick_start_i + 1 : subnewick_stop_i - 1]
            if len(subnewick) == 0:
                raise InvalidFileError(
                    "Intermediate node must at least have one child"
                )
            
            comma_pos = []
            level = 0
            for i, char in enumerate(subnewick):
                if char == "(":
                    level += 1
                elif char == ")":
                    level -= 1
                elif char == ",":
                    if level == 0:
                        comma_pos.append(i)
                if level < 0:
                    raise InvalidFileError(
                        "Bracket closed before it was opened"
                    )
        
            children = []
            distances = []
            for i, pos in enumerate(comma_pos):
                if i == 0:
                    child, dist = TreeNode.from_newick(
                        subnewick[:pos], labels=labels
                    )
                else:
                    prev_pos = comma_pos[i - 1]
                    child, dist = TreeNode.from_newick(
                        subnewick[prev_pos + 1 : pos], labels=labels
                    )
                children.append(child)
                distances.append(dist)
            
            if len(comma_pos) != 0:
                child, dist = TreeNode.from_newick(
                    subnewick[comma_pos[-1] + 1:], labels=labels
                )
            else:
                child, dist = TreeNode.from_newick(
                    subnewick, labels=labels
                )
            children.append(child)
            distances.append(dist)
            return TreeNode(children, distances), distance

    def __str__(self):
        return self.to_newick()
    
    def __eq__(self, item):
        if not isinstance(item, TreeNode):
            return False
        node = item
        if self._distance != node._distance:
            return False
        if self._index != -1:
            if self._index != node._index:
                return False
        else:
            if frozenset(self._children) != frozenset(node._children):
                return False
        return True
    
    def __hash__(self):
        children_set = frozenset(self._children) \
                       if self._children is not None else None
        return hash((self._index, children_set, self._distance))


def _get_leaves(node, leaf_list):
    if node._index == -1:
        for child in node._children:
            _get_leaves(child, leaf_list)
    else:
        leaf_list.append(node)


def _create_path_to_root(node):
    path = []
    current_node = node
    while current_node is not None:
        path.append(current_node)
        current_node = current_node._parent
    return path


class TreeError(Exception):
    pass