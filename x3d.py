# Desenvolvido por: Luciano Soares <lpsoares@insper.edu.br>
# Disciplina: Computação Gráfica
# Data: 31 de Agosto de 2020

# XML
import xml.etree.ElementTree as ET

# Outras
import re
import math

def clean(child):
    _, _, child.tag = child.tag.rpartition('}') # remove os namespaces

class X3D:

    current_color = [1.0, 1.0, 1.0] # controle de cor instantânea
    preview = None # artibuto que aponta para o sistema de preview
    render = {} # dicionario dos métodos de renderização

    def __init__(self, filename):
        self.x3d = ET.parse(filename)
        self.root = self.x3d.getroot()

    def set_preview(self, preview):
        X3D.preview = preview

    def set_resolution(self, width, height):
        self.width = width
        self.height = height

    def parse(self):
        """ parse começando da raiz do X3D. """
        for child in self.root:
            clean(child) # remove namespace
            if child.tag == "Scene":
                self.scene = Scene(child)

class Scene:
    def __init__(self, node):
        self.children = []
        for child in node:
            clean(child) # remove namespace
            if child.tag == "Transform":
                self.children.append(Transform(child))
            elif child.tag == "Viewpoint":
                self.children.append(Viewpoint(child))

# Core component

class X3DNode:
    pass

class X3DChildNode(X3DNode):
    def __init__(self):
        super().__init__()

class X3DBindableNode(X3DChildNode):
    def __init__(self):
        super().__init__()

# Grouping component

class X3DGroupingNode(X3DChildNode):
    def __init__(self):
        super().__init__()
        self.children = []

class Transform(X3DGroupingNode):
    def __init__(self, node):
        super().__init__()
        self.rotation = [0, 0, 1, 0]
        self.scale = [1, 1, 1]
        self.translation = [0, 0, 0]

        if 'rotation' in node.attrib:
            rotation_str = re.split(r'[,\s]\s*',node.attrib['rotation'])
            self.rotation = [ float(value) for value in rotation_str]

        if 'scale' in node.attrib:
            scale_str = re.split(r'[,\s]\s*',node.attrib['scale'])
            self.scale = [ float(value) for value in scale_str]

        if 'translation' in node.attrib:
            translation_str = re.split(r'[,\s]\s*',node.attrib['translation'])
            self.translation = [ float(value) for value in translation_str]

        # Render
        if "Transform" in X3D.render:
            X3D.render["Transform"](translation=self.translation, scale=self.scale, rotation=self.rotation)

        for child in node:
            clean(child) # remove namespace
            if child.tag == "Shape":
                self.children.append(Shape(child))
            elif child.tag == "Transform":
                self.children.append(Transform(child))

# Shape component

class X3DShapeNode(X3DChildNode):
    def __init__(self):
        super().__init__()
        self.geometry = None
        self.appearance = None

class X3DAppearanceNode(X3DNode):
    def __init__(self):
        super().__init__()

class X3DMaterialNode():
    def __init__(self):
        super().__init__()

class Material(X3DMaterialNode):
    def __init__(self, node):
        super().__init__()
        self.diffuseColor = [0.8, 0.8, 0.8]
        if 'diffuseColor' in node.attrib:
            diffuseColor_str = re.split(r'[,\s]\s*',node.attrib['diffuseColor'])
            self.diffuseColor = [ float(color) for color in diffuseColor_str]
        X3D.current_color = self.diffuseColor

class Appearance(X3DAppearanceNode):
    def __init__(self, node):
        super().__init__()
        self.material = None
        for child in node:
            clean(child) # remove namespace
            if child.tag == "Material":
                self.material = Material(child)

class Shape(X3DShapeNode):
    def __init__(self, node):
        super().__init__()
        for child in node:
            clean(child) # remove namespace
            if child.tag == "Appearance":
                self.appearance = Appearance(child)
        for child in node:
            clean(child) # remove namespace
            if child.tag == "Polypoint2D":
                self.geometry = Polypoint2D(child)
            elif child.tag == "Polyline2D":
                self.geometry = Polyline2D(child)
            elif child.tag == "TriangleSet2D":
                self.geometry = TriangleSet2D(child)
            elif child.tag == "TriangleSet":
                self.geometry = TriangleSet(child)


# Rendering component

class X3DGeometryNode(X3DNode):
    def __init__(self):
        super().__init__()

class X3DComposedGeometryNode(X3DGeometryNode):
    def __init__(self):
        super().__init__()

class X3DGeometricPropertyNode(X3DNode):
    def __init__(self):
        super().__init__()

class X3DCoordinateNode(X3DGeometricPropertyNode):
    def __init__(self):
        super().__init__()

class Coordinate(X3DCoordinateNode):
    def __init__(self, node):
        super().__init__()
        point_str = re.split(r'[,\s]\s*',node.attrib['point'].strip())
        self.point = [ float(p) for p in point_str]

# Geometry2D component

class Polypoint2D(X3DGeometryNode):
    def __init__(self, node):
        super().__init__()
        point_str = re.split(r'[,\s]\s*',node.attrib['point'].strip())
        self.point = [ float(p) for p in point_str]

        # Preview
        if X3D.preview:
            polypoint2D = []
            for i in range(0, len(self.point), 2):
                polypoint2D.append([self.point[i], self.point[i+1]])
            X3D.preview._pontos.append({'color': X3D.current_color, 'points': polypoint2D})

        # Render
        if "Polypoint2D" in X3D.render:
            X3D.render["Polypoint2D"](point=self.point, color=X3D.current_color)

class Polyline2D(X3DGeometryNode):
    def __init__(self, node):
        super().__init__()
        lineSegments_str = re.split(r'[,\s]\s*',node.attrib['lineSegments'].strip())
        self.lineSegments = [ float(point) for point in lineSegments_str]

        # Preview
        if X3D.preview:
            polyline2D = []
            for i in range(0, len(self.lineSegments), 2):
                polyline2D.append([self.lineSegments[i], self.lineSegments[i+1]])
            X3D.preview._linhas.append({'color': X3D.current_color, 'lines': polyline2D})

        # Render
        if "Polyline2D" in X3D.render:
            X3D.render["Polyline2D"](lineSegments=self.lineSegments, color=X3D.current_color)

class TriangleSet2D(X3DGeometryNode):
    def __init__(self, node):
        super().__init__()
        vertices_str = re.split(r'[,\s]\s*',node.attrib['vertices'].strip())
        self.vertices = [ float(point) for point in vertices_str]

        # Preview
        if X3D.preview:
            for i in range(0, len(self.vertices), 6):
                X3D.preview._poligonos.append({'color': X3D.current_color,
                                                    'vertices': [[self.vertices[i  ], self.vertices[i+1]],
                                                                 [self.vertices[i+2], self.vertices[i+3]],
                                                                 [self.vertices[i+4], self.vertices[i+5]]]})

        # Render
        if "TriangleSet2D" in X3D.render:
            X3D.render["TriangleSet2D"](vertices=self.vertices, color=X3D.current_color)

class TriangleSet(X3DComposedGeometryNode):
    def __init__(self, node):
        super().__init__()
        self.coord = None
        for child in node:
            clean(child) # remove namespace
            if child.tag == "Coordinate":
                self.coord = Coordinate(child)

        # Preview
        # Implemente se desejar
        
        # Render
        if "TriangleSet" in X3D.render:
            X3D.render["TriangleSet"](point=self.coord.point, color=X3D.current_color)

# Navigation component

class X3DViewpointNode(X3DBindableNode):
    def __init__(self):
        super().__init__()

class Viewpoint(X3DViewpointNode):
    def __init__(self, node):
        super().__init__()
        self.position = [0, 0, 10]         # Valores padrão
        self.orientation = [0, 0, 1, 0]    # Valores padrão
        self.fieldOfView = math.pi/4       # Valores padrão
        if 'position' in node.attrib:
            position_str = re.split(r'[,\s]\s*',node.attrib['position'].strip())
            self.position = [ float(point) for point in position_str]

        if 'orientation' in node.attrib:
            orientation_str = re.split(r'[,\s]\s*',node.attrib['orientation'].strip())
            self.orientation = [ float(point) for point in orientation_str]

        if 'fieldOfView' in node.attrib:
            self.fieldOfView = float(node.attrib['fieldOfView'].strip())

        # Render
        if "Viewpoint" in X3D.render:
            X3D.render["Viewpoint"](position=self.position, orientation=self.orientation, fieldOfView=self.fieldOfView)
