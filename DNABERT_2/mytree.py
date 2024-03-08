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