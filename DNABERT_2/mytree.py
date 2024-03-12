class TreeNode():
    def __init__(self, value) -> None:
        self.value = value
        self.children = []


def print_tree(node, level=0):
    if node.value is not None:
        try:
            print('  ' * level + node.value)
        except Exception as e:
            print(e)
            print(f"{node.value=}, {node.children=}")
            return
    for child in node.children:
        print_tree(child, level + 1)


def add_node(root, value, depth):
    current = root
    
    for _ in range(depth - 1):
        if current.children:
            current = current.children[-1]
        else:
            new_child = TreeNode(None)
            current.children.append()
            current = new_child
    
    new_node = TreeNode(value)
    current.children.append(new_node)


def cal_tree_distance(root, node1, node2):
    lca = find_lca(root, node1, node2)
    distance1 = find_distance_to_ancestor(lca, node1)
    distance2 = find_distance_to_ancestor(lca, node2)
    return distance1 + distance2

def find_lca(root, node1, node2):
    if root is None:
        return None
    
    if not root.value is None:
        if root.value.lower() == node1.value.lower() or root.value.lower() == node2.value.lower():
            return root
    
    lca_node = None
    for child in root.children:
        lca = find_lca(child, node1, node2)
        if lca:
            if lca_node:
                return root
            lca_node = lca
    
    return lca_node

def find_distance_to_ancestor(ancestor, node):
    if ancestor is None:
        return -1
    
    if ancestor.value.lower() == node.value.lower():
        return 0
    
    for child in ancestor.children:
        distance = find_distance_to_ancestor(child, node)
        if distance != -1:
            return distance + 1
    
    return -1

def find_deepest_node(root, name_begin):
    if not root:
        return None

    deepest_node = None
    max_depth = -1

    def dfs(node, depth):
        nonlocal deepest_node, max_depth

        if (not node.value is None) and depth > max_depth:
            if (node.value.startswith(name_begin) and depth == 0) or depth > 0:
                max_depth = depth
                deepest_node = node

        for child in node.children:
            dfs(child, depth + 1)

    dfs(root, 0)
    return deepest_node, max_depth