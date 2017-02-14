class MaskObjNodeTree:
    '''
    This node_tree object is made to contain the nodes which contain masks
    produced by fil_finder.
    Each tree contains a starting node and a string of nodes which define the
    tree.

    The first node contains the aggregated mask and is merged with each a
    subsequent addition
    '''

    def __init__(self, node_obj):
        self.root_node = node_obj
        self.root_v_slice = node_obj.v_slice_index[0]

        self.node_list = [node_obj]
        self.length = 1

        self.has_ended = False

    def addNode(self, new_node):
        self.root_node.mergeNode(new_node)
        self.node_list.append(new_node)
        self.length += 1
        return self.length

    def getNode(self, node_number):
        return self.node_list[node_number]

    def getLastNode(self):
        return self.getNode(self.length - 1)

    def getTreeMask(self):
        return self.root_node.mask

    '''
    def removeLastNode(self):
    def visitAllNodes(self):
    '''
