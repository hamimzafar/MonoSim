import abc

class Tree(object):
    """
    Represents a tree
    """
    def __init__(self, name, node_props = {}):
        """
        Creates a tree

        Arguments
        name - an identifier for the root of the tree
        data - any information that one wants in the form of a dictionary
        """
        self.node_props = node_props
        self.name = name
        self.children = set([])

    def get_name(self):
        """
        :return: teh name of the root node of teh tree
        """
        return self.name

    def add_children(self, children):
        """
        :param child: the subtree that we are about to append to the original tree
        :return: updates the tree to have a new child
        """
        for child in children:
            self.children.add(child)

    def remove_children(self, delete):
        """
        :param remove: the names of the children that the user wants to remove
        :return: updates the tree to have the children of the specified node removed
        """
        for child in delete:
            self.children.remove(child)

    def get_children(self):
        """
        :return: the children of the node that you are currently at
        """
        return list(self.children)

    def set_node_prop(self, key, val):
        """
        :param key: the key of the property that the user is specifying
        :param val: the value of the property that the user is specifying
        :return: the altered node properties
        """
        self.node_props[key] = val

    def get_node_prop(self, key):
        """
        :param key: the key of the property the user is trying to reach
        :return: the value of what is stored in that key of the node properties
        """
        return self.node_props[key]

    def to_string(self, owner):
        "Contract from super."
        if owner.get_children() != []:
            new = str(owner.get_name()) + '('
            for child in owner.get_children():
                new += self.to_string(child)
                new += ')'
        else:
            new = str(owner.get_name())
        return new

    def __str__(self):
        """
        :return: what will be printed out when one prints the tree
        """
        return self.to_string(self)

class FullBiTree(object):
    """
    Represents a full binary tree.
    """

    def __init__(self, name, left_tree = None, right_tree = None):
        """
        Creates a full binary tree.

        This constructor must be called with exactly one or three parameters.
        That is, a name alone or a name and both a left and right child.

        Arguments:
        name - an identifier for the root node of the tree.
        left_tree - the FullBiTree left substree if the tree's root has children. (optional)
        right_tree - the FullBiTree left substree if the tree's root has children. (optional)
        """

        self.__name = name
        self.__node_props = {}
        if left_tree == None and right_tree == None:
            self.__set_state(TreeNodeStateLeaf())
        elif left_tree != None and right_tree != None:
            self.__set_state(TreeNodeStateInternal(left_tree, right_tree))
        else:
            raise Exception("A full binary tree must have either 0 or 2 children")

    def get_name(self):
        """
        Gets the name of the root node of the tree.

        Returns:
        The name of the root node.
        """
        return self.__name

    def get_left_child(self):
        """
        Gets the left subtree of the tree's root if it has children or generates an exception if the root has no children.

        Returns:
        The left subtree of the tree.
        """
        return self.__get_state().get_left_child()

    def get_all_props(self):
        return self.__node_props

    def get_right_child(self):
        """
        Gets the right subtree of the tree's root if it has children or generates an exception if the root has no children.

        Returns:
        The left subtree of the tree.
        """
        return self.__get_state().get_right_child()

    def set_children(self, left_tree, right_tree):
        """
        Updates the tree's root to contain new children.

        Arguments:
        left_tree - the new left subtree for the tree.
        right_tree - the new right subtree for the tree.
        """
        self.__set_state(TreeNodeStateInternal(left_tree, right_tree))

    def remove_children(self):
        """
        Updates the tree's root to contain no children.

        Arguments:
        left_tree - the new left subtree for the tree.
        right_tree - the new right subtree for the tree.
        """
        self.__set_state(TreeNodeStateLeaf())

    def is_leaf(self):
        """
        Tests whether the tree's root has no children.

        Returns:
        True if the tree is only a single node, else false.
        """
        return self.__get_state().is_leaf()

    def __set_state(self, new_state):
        """
        Sets the internal node/leaf node state for the node.

        Arguments:
        new_state - the new node state.
        """
        self.__node_state = new_state

    def __get_state(self):
        """
        Gets the internal node/leaf node state for the node.

        Returns:
        The current node state.
        """
        return self.__node_state

    def __str__(self):
        " Contract from super. "
        return self.__get_state().to_string(self)

    def get_node_property(self, key):
        """
        Accesses a user specified property of the tree's root.

        Arguments:
        key - the property of the desired key value pair.

        Returns:
        The value of the given key for the tree's root.
        """
        return self.__node_props[key]

    def set_node_property(self, key, value):
        """
        Defines a user specified property of the tree's root.

        Arguments:
        key - the key of the desired property.
        value - the value of the desired property.
        """
        self.__node_props[key] = value

    def get_left_edge_property(self, key):
        """
        Accesses a user specified property of the tree's left subtree edge.
        Throws exception if the tree has no left subtree.

        Arguments:
        key - the property of the desired key value pair.

        Returns:
        The value of the given key for the tree's left subtree edge.
        """
        return self.__get_state().get_left_edge_property(key)

    def set_left_edge_property(self, key, value):
        """
        Defines a user specified property of the tree's left subtree edge.
        Throws exception if the tree has no left subtree.

        Arguments:
        key - the key of the desired property.
        value - the value of the desired property.
        """
        self.__get_state().set_left_edge_property(key, value)

    def get_right_edge_property(self, key):
        """
        Accesses a user specified property of the tree's right subtree edge.
        Throws exception if the tree has no left subtree.

        Arguments:
        key - the property of the desired key value pair.

        Returns:
        The value of the given key for the tree's right subtree edge.
        """
        return self.__get_state().get_right_edge_property(key)

    def set_right_edge_property(self, key, value):
        """
        Defines a user specified property of the tree's right subtree edge.
        Throws exception if the tree has no left subtree.

        Arguments:
        key - the key of the desired property.
        value - the value of the desired property.
        """
        self.__get_state().set_right_edge_property(key, value)

class TreeNodeState(object):
    """
    Abstract class for defining all operations for a node state.
    """

    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def is_leaf(self):
        """
        Tests whether the node state represents a leaf.

        Returns:
        True if the node state represents a leaf, else false.
        """
        pass

    @abc.abstractmethod
    def to_string(self, owner):
        """
        Returns a prefix string representation of the whole tree rooted by the node state.

        Returns:
        A prefix string representation of the tree.
        """
        pass

    @abc.abstractmethod
    def get_left_child(self):
        """
        Returns the left child of this node if in the internal state, or generate exeption if in leaf state.

        Returns:
        The left subtree.
        """
        pass

    @abc.abstractmethod
    def get_right_child(self):
        """
        Returns the right child of this node if in the internal state, or generate exeption if in leaf state.

        Returns:
        The right subtree.
        """
        pass

    @abc.abstractmethod
    def get_left_edge_property(self, key):
        """
        Accesses a user specified property of the node state's left subtree edge.
        Throws exception if the tree has no left subtree.

        Arguments:
        key - the property of the desired key value pair.

        Returns:
        The value of the given key for the tree's left subtree edge.
        """
        pass

    @abc.abstractmethod
    def set_left_edge_property(self, key, value):
        """
        Accesses a user specified property of the node state's left subtree edge.
        Throws exception if the node state has no left subtree.

        Arguments:
        key - the property of the desired key value pair.

        Returns:
        The value of the given key for the tree's right subtree edge.
        """
        pass

    @abc.abstractmethod
    def get_right_edge_property(self, key):
        """
        Accesses a user specified property of the node state's right subtree edge.
        Throws exception if the tree has no right subtree.

        Arguments:
        key - the property of the desired key value pair.

        Returns:
        The value of the given key for the tree's right subtree edge.
        """
        pass

    @abc.abstractmethod
    def set_right_edge_property(self, key, value):
        """
        Accesses a user specified property of the node state's right subtree edge.
        Throws exception if the node state has no left subtree.

        Arguments:
        key - the property of the desired key value pair.

        Returns:
        The value of the given key for the tree's right subtree edge.
        """
        pass

class TreeNodeStateLeaf(TreeNodeState):
    """
    TreeNodeState representing a leaf.
    """

    def is_leaf(self):
        "Contract from super."
        return True

    def to_string(self, owner):
        "Contract from super."
        return str(owner.get_name())

    def get_left_child(self):
        "Contract from super."
        raise Exception("A leaf does not have a left child.")

    def get_right_child(self):
        "Contract from super."
        raise Exception("A leaf does not have a right child.")

    def get_left_edge_property(self, key):
        "Contract from super."
        raise Exception("A leaf does not have a left edge.")

    def set_left_edge_property(self, key, value):
        "Contract from super."
        raise Exception("A leaf does not have a left edge.")

    def get_right_edge_property(self, key):
        "Contract from super."
        raise Exception("A leaf does not have a right edge.")

    def set_right_edge_property(self, key, value):
        "Contract from super."
        raise Exception("A leaf does not have a right edge.")

class TreeNodeStateInternal(TreeNodeState):
    """
    TreeNodeState for an internal node.
    """

    def __init__(self, left_tree, right_tree):
        """
        Creates a new TreeNodeState instance.

        Arguments:
        left_tree - The FullBiTree left subtree of this node.
        right_tree - The FullBiTree right subtree of this node.
        """
        self.__left_tree = left_tree
        self.__right_tree = right_tree
        self.__left_edge_props = {}
        self.__right_edge_props = {}

    def is_leaf(self):
        "Contract from super."
        return False

    def get_left_child(self):
        "Contract from super."
        return self.__left_tree

    def get_right_child(self):
        "Contract from super."
        return self.__right_tree

    def get_left_edge_property(self, key):
        "Contract from super."
        return self.__left_edge_props[key]

    def set_left_edge_property(self, key, value):
        "Contract from super."
        self.__left_edge_props[key] = value

    def get_right_edge_property(self, key):
        "Contract from super."
        return self.__right_edge_props[key]

    def set_right_edge_property(self, key, value):
        "Contract from super."
        self.__right_edge_props[key] = value

    def to_string(self, owner):
        "Contract from super."
        return str(owner.get_name()) + '(' + str(self.get_left_child()) + ', ' + str(self.get_right_child()) + ')'

"""
Counts the number of leaves in a tree in dict form.
"""
def count_leaves_full_bitree_dict(tree):
    num_leaves = 0
    for adj_set in tree.values():
        adj_set_size = len(adj_set)
        if adj_set_size == 0 or adj_set_size == 1:
            num_leaves+=1

    return num_leaves

"""
Counts the number of leaves in a FullBiTree
"""
def count_leaves_fullbitree(tree):
    if tree.is_leaf():
        return 1
    else:
        return count_leaves_full_bitree_dict(tree.get_right_child()) + \
               count_leaves_full_bitree_dict(tree.get_left_child())


"""
Computes the hight of a FullBiTree.

Arguments:
tree - a full binary tree in FullBiTree form.

Returns:
The height of the tree
"""
def tree_height(tree):
    if tree.is_leaf():
        return 0
    else:
        left_height = tree_height(tree.get_left_child())
        right_height = tree_height(tree.get_right_child())
        if left_height > right_height:
            return left_height + 1
        else:
            return right_height + 1

def infix_string(tree):
    """
    Computes the infix order string of a tree.

    Arguments:
    tree - a full binary tree in FullBiTree form.

    Returns:
    An infix string of the tree.
    """
    if tree.is_leaf():
        return tree.get_name()
    else:
        return infix_string(tree.get_left_child()) + tree.get_name() + infix_string(tree.get_right_child())


def find_paths(tree):
    """
    Computes a string for each path in the givn tree starting at the root and terminating at a leaf.

    Arguments:
    tree - a FullBiTree

    Returns:
    A set of strings encoding the order of nodes in each path from the root to all leaves.
    """
    found_paths = set()
    find_paths_help(tree, "", found_paths)
    return found_paths

def find_paths_help(tree, path_so_far, found_paths):
    """
    Computes a string for each path in the given sub-tree starting at the sub-tree root and terminating at a leaf.
    Stores completed paths in the given set

    Arguments:
    tree - a FullBiTree
    path_so_far - a string encoding the path seen so far from the global tree root to this sub-tree root.
    found_paths - all complete paths seen so far from the tree root to a leaf.
    """
    if tree.is_leaf():
        path = path_so_far + tree.get_name()
        found_paths.add(path)
    else:
        find_paths_help(tree.get_left_child(), path_so_far + tree.get_name(), found_paths)
        find_paths_help(tree.get_right_child(), path_so_far + tree.get_name(), found_paths)

def is_valid_bst(tree):
    """
    Tests to see if the given tree has the binary search property.

    Arguments:
    tree - a FullBiTree where the value of each node is an integer stored as the node's name.

    Returns:
    True if the tree has the binary search property, else false.
    """
    infix_list = list()
    is_valid_bst_help(tree, infix_list)

    prev_element = infix_list[0]

    for element in infix_list:
        if element < prev_element:
            return False
        prev_element = element

    return True

def is_valid_bst_help(tree, infix_list):
    if tree.is_leaf():
        infix_list.append(tree.get_name())
    else:
        is_valid_bst_help(tree.get_left_child(), infix_list)
        infix_list.append(tree.get_name())
        is_valid_bst_help(tree.get_right_child(), infix_list)
