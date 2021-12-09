from make_prg.recursion_tree import RecursiveTreeNode, SingleClusterNode
import networkx as nx
import matplotlib.pyplot as plt
from pathlib import Path


class RecursiveTreeDrawer:
    def __init__(self, root: RecursiveTreeNode):
        self._graph: nx.DiGraph = self._build_recursive_tree_graph(root)

    @staticmethod
    def _build_recursive_tree_graph(root: RecursiveTreeNode) -> nx.DiGraph:
        graph = nx.DiGraph()
        RecursiveTreeDrawer._preorder_visit(root, graph)
        return graph

    @staticmethod
    def _preorder_visit(node: RecursiveTreeNode, graph: nx.DiGraph):
        RecursiveTreeDrawer._visit(node, graph)
        for child in node.children:
            RecursiveTreeDrawer._preorder_visit(child, graph)

    @staticmethod
    def _visit(node: RecursiveTreeNode, graph: nx.DiGraph):
        node_attributes = {}
        node_is_SingleClusterNode = isinstance(node, SingleClusterNode)

        if node.is_leaf():
            assert node_is_SingleClusterNode, "Error, a leaf is not a SingleClusterNode"
            prg_as_list = []
            node._get_prg(prg_as_list)
            node_attributes["label"] = "".join(prg_as_list)
        else:
            node_attributes["label"] = str(node.id)

        node_attributes["color"] = "blue" if node_is_SingleClusterNode else "red"
        graph.add_node(node, **node_attributes)

        if node.parent is not None:
            graph.add_edge(node.parent, node)

    def output_graph(self, filename: Path):
        plt.figure(figsize=(20, 10))
        a_graph = nx.drawing.nx_agraph.to_agraph(self._graph)
        a_graph.layout(prog="dot")
        a_graph.draw(filename)
