
class Xml():

    @staticmethod
    def get_child_by_tag_name(node, tag_name):
        print(f"sl: node.childNodes: {node.childNodes}")
        for child in node.childNodes:
            if node.nodeType == "ELEMENT_NODE":
                if child.tagName == tag_name:
                    return child
