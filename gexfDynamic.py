class DynamicGraph:
    def __init__(self, nodes, edges, nodesizes):
        '''Initializes an object of type DynamicGraph'''
        '''@arg nodes: dict with key=nodeId and value=lifetime list'''
        '''@arg edges: dict with key=edgeId and value=lifetime list'''
        '''@arg nodesizes: dict with key=nodeId and value=list of values'''
        self.nodes = nodes
        self.edges = edges
        self.nodesizes = nodesizes


    def makeAttributesDeclaration(self):
        '''Makes the part which declares the attributes'''
        return '<attributes class ="node" mode ="dynamic">\n' \
               ' <attribute id="0" title ="size" type ="float"/> ' \
               '\n</attributes>'

    def makeNodeSizeAttributes(self,id):
        '''Makes the part where the size attribute of a node are declared'''
        '''@arg id: id of the node'''
        sizes = self.nodesizes[id]
        res = ""
        for i in range(len(sizes)): res += '<attvalue for="{attfor}" value="{val}" start="{start}" end="{end}"/>\n'.format(attfor=0,val=sizes[i],start=(i-1),end=i)
        return ('<attvalues>\n' + res + '</attvalues>')

    def makeNodeSpells(self,id):
        '''Makes the part where the lifespans of a node are declared'''
        '''@arg id: id of the node'''
        lifetime = self.nodes[id]
        res = ""
        for i in range(len(lifetime)):
            if i%2 == 0:
                if len(lifetime)-i>=2:
                    res += '<spell start="{start}" end="{end}" />\n'.format(start=lifetime[i], end=lifetime[i+1])
                else:
                    res += '<spell start="{start}" />\n'.format(start=lifetime[i])
        return ('<spells>\n' + res + '</spells>')

    def makeEdgeSpells(self,id):
        '''Makes the part where the lifespans of an edge are declared'''
        '''@arg id: id of the node'''
        lifetime = self.edges[id]
        res = ""
        for i in range(len(lifetime)):
            if i%2 == 0:
                if len(lifetime)-i>=2:
                    res += '<spell start="{start}" end="{end}" />\n'.format(start=lifetime[i], end=lifetime[i+1])
                else:
                    res += '<spell start="{start}" />\n'.format(start=lifetime[i])
        return ('<spells>\n' + res + '</spells>')

    def makeNode(self, id):
        '''Makes a whole description of a node'''
        '''@arg id: id of the node'''
        return '<node id="{id}">\n'.format(id=id) + self.makeNodeSizeAttributes(id) + '\n' + self.makeNodeSpells(id) + '\n</node>'

    def makeEdge(self, id):
        '''Makes a whole description of an edge'''
        '''@arg id: id of the node'''
        return '<edge id="{id}" source="{source}" target="{target}">\n'.format(id=id, target=int(int(id)/len(self.nodes)), source=int(id)%len(self.nodes)) + self.makeEdgeSpells(id) + '\n</edge>'

    def makeNodes(self):
        '''Makes the nodes part'''
        res = ""
        for key in self.nodes:
            res += self.makeNode(key) + "\n"
        return '<nodes>\n' + res + '</nodes>'

    def makeEdges(self):
        '''Makes the edges part'''
        res = ""
        for key in self.edges:
            res += self.makeEdge(key) + "\n"
        return '<edges>\n' + res + '</edges>'

    def makeGraph(self):
        '''Makes the whole graph description including the header'''
        return '<?xml' \
        'version = "1.0"' \
        ' encoding = "utf-8"?> <gexf' \
        ' version = "1.2"' \
        ' xmlns = "http://www.gexf.net/1.2draft"' \
        ' xmlns:xsi = "http://www.w3.org/2001/XMLSchema-instance"' \
        ' xsi:schemaLocation = "http://www.w3.org/2001/XMLSchema-instance">\n<graph' \
        ' defaultedgetype = "directed"' \
        ' mode = "dynamic"' \
        ' name = "Dynamischer Graph"' \
        ' timeformat = "long" >\n' \
        + self.makeAttributesDeclaration() + '\n'\
        + self.makeNodes() + '\n' \
        + self.makeEdges() + '\n'  \
        + '</graph>\n</gexf>'


# Test: Uncomment to test
#nodes = {'0':[-1,5,6,8,12,14],'1':[-1],'2':[3,9,15]}
#edges = {'1':[-1,3,6,8,12,13],'7':[-1]}
#nodesizes = {'0':[4.5,7.3,4,6,23,6,97,3,5,7.6],'1':[4,2,6,2,8,9,5,4],'2':[2,42,56,3,345,5,10,5.4]}
#d = DynamicGraph(nodes, edges, nodesizes)
# print(d.makeNodeSizeAttributes('0'))
# print(d.makeNodeSizeAttributes('1'))
# print(d.makeNodeSizeAttributes('2'))
# print(d.makeNodeSpells('0'))
# print(d.makeNodeSpells('1'))
# print(d.makeNodeSpells('2'))
#print(d.makeGraph())