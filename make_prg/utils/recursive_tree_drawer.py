from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx

from make_prg.recursion_tree import (
    LeafNode,
    MultiClusterNode,
    MultiIntervalNode,
    RecursiveTreeNode,
)


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

        if node.is_leaf():
            prg_as_list = []
            node.preorder_traversal_to_build_prg(prg_as_list, do_indexing=False)
            node_attributes["label"] = "".join(prg_as_list)
        else:
            node_attributes["label"] = str(node.node_id)

        node_attributes["color"] = RecursiveTreeDrawer._get_node_colour(node)
        graph.add_node(node, **node_attributes)

        if node.parent is not None:
            graph.add_edge(node.parent, node)

    @staticmethod
    def _get_node_colour(node: RecursiveTreeNode) -> str:
        if isinstance(node, MultiClusterNode):
            return "red"
        elif isinstance(node, MultiIntervalNode):
            return "blue"
        elif isinstance(node, LeafNode):
            return "green"
        else:
            assert False, f"Unknown node class when drawing: {node.__class__}"

    def output_graph(self, filename: Path):
        plt.figure(figsize=(20, 10))
        a_graph = nx.drawing.nx_agraph.to_agraph(self._graph)
        a_graph.layout(prog="dot")
        a_graph.draw(filename)
