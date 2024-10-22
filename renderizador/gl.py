#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# pylint: disable=invalid-name

"""
Biblioteca Gráfica / Graphics Library.

Desenvolvido por: Marcelo Rabello Barranco
Disciplina: Computação Gráfica
Data: 12/08/2024
"""

import time         # Para operações com tempo
import gpu          # Simula os recursos de uma GPU
import math         # Funções matemáticas
import numpy as np  # Biblioteca do Numpy

# Sendo bem sincero eu tentei, fiz de tudo pra fazer funcionar, mas nao foi, ent escrevi esse comentario
# Eu fiz, refiz, fiz e refiz e nada funcionou, ent nao vou nem dar commit da merda q eu fiz, pq ta tudo errado

class GL:
    """Classe que representa a biblioteca gráfica (Graphics Library)."""

    width = 800   # largura da tela
    height = 600  # altura da tela
    near = 0.01   # plano de corte próximo
    far = 1000    # plano de corte distante
    transformation_stack = [] # pilha de transformações


    @staticmethod
    def setup(width, height, near=0.01, far=1000):
        """Definr parametros para câmera de razão de aspecto, plano próximo e distante."""
        GL.width = width
        GL.height = height
        GL.near = near
        GL.far = far
        GL.view_matrix = np.identity(4)
        GL.perspective_matrix = np.identity(4)

    @staticmethod
    def polypoint2D(point, colors):
        """Função usada para renderizar Polypoint2D."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry2D.html#Polypoint2D
        # Nessa função você receberá pontos no parâmetro point, esses pontos são uma lista
        # de pontos x, y sempre na ordem. Assim point[0] é o valor da coordenada x do
        # primeiro ponto, point[1] o valor y do primeiro ponto. Já point[2] é a
        # coordenada x do segundo ponto e assim por diante. Assuma a quantidade de pontos
        # pelo tamanho da lista e assuma que sempre vira uma quantidade par de valores.
        # O parâmetro colors é um dicionário com os tipos cores possíveis, para o Polypoint2D
        # você pode assumir inicialmente o desenho dos pontos com a cor emissiva (emissiveColor).

        color = [int(255 * colors["emissiveColor"][0]),
                int(255 * colors["emissiveColor"][1]),
                int(255 * colors["emissiveColor"][2])]

        x = int(point[0])
        y = int(point[1])

        if 0 <= x < GL.width and 0 <= y < GL.height:
            gpu.GPU.draw_pixel([x, y], gpu.GPU.RGB8, color)
        
    @staticmethod
    def polyline2D(lineSegments, colors):
        """Função usada para renderizar Polyline2D."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry2D.html#Polyline2D
        # Nessa função você receberá os pontos de uma linha no parâmetro lineSegments, esses
        # pontos são uma lista de pontos x, y sempre na ordem. Assim point[0] é o valor da
        # coordenada x do primeiro ponto, point[1] o valor y do primeiro ponto. Já point[2] é
        # a coordenada x do segundo ponto e assim por diante. Assuma a quantidade de pontos
        # pelo tamanho da lista. A quantidade mínima de pontos são 2 (4 valores), porém a
        # função pode receber mais pontos para desenhar vários segmentos. Assuma que sempre
        # vira uma quantidade par de valores.
        # O parâmetro colors é um dicionário com os tipos cores possíveis, para o Polyline2D
        # você pode assumir inicialmente o desenho das linhas com a cor emissiva (emissiveColor).

        x_pos_1 = lineSegments[0]
        y_pos_1 = lineSegments[1]
        x_pos_2 = lineSegments[2]
        y_pos_2 = lineSegments[3]

        dx = x_pos_2 - x_pos_1
        dy = y_pos_2 - y_pos_1

        steps = abs(dx) if abs(dx) > abs(dy) else abs(dy)

        x_inc = dx / steps
        y_inc = dy / steps

        color = [int(255 * colors["emissiveColor"][0]),
                int(255 * colors["emissiveColor"][1]),
                int(255 * colors["emissiveColor"][2])]

        for i in range(int(steps)):
            if 0 <= int(x_pos_1) < GL.width and 0 <= int(y_pos_1) < GL.height:
                gpu.GPU.draw_pixel([int(x_pos_1), int(y_pos_1)], gpu.GPU.RGB8, color)
            x_pos_1 += x_inc
            y_pos_1 += y_inc

    @staticmethod
    def circle2D(radius, colors):
        """Função usada para renderizar Circle2D."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry2D.html#Circle2D
        # Nessa função você receberá um valor de raio e deverá desenhar o contorno de
        # um círculo.
        # O parâmetro colors é um dicionário com os tipos cores possíveis, para o Circle2D
        # você pode assumir o desenho das linhas com a cor emissiva (emissiveColor).

        print("Circle2D : radius = {0}".format(radius)) # imprime no terminal
        print("Circle2D : colors = {0}".format(colors)) # imprime no terminal as cores
        
        # Exemplo:
        pos_x = GL.width//2
        pos_y = GL.height//2
        gpu.GPU.draw_pixel([pos_x, pos_y], gpu.GPU.RGB8, [255, 0, 255])  # altera pixel (u, v, tipo, r, g, b)
        # cuidado com as cores, o X3D especifica de (0,1) e o Framebuffer de (0,255)


    @staticmethod
    def triangleSet2D(vertices, colors):
        """Função usada para renderizar TriangleSet2D."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry2D.html#TriangleSet2D
        # Nessa função você receberá os vertices de um triângulo no parâmetro vertices,
        # esses pontos são uma lista de pontos x, y sempre na ordem. Assim point[0] é o
        # valor da coordenada x do primeiro ponto, point[1] o valor y do primeiro ponto.
        # Já point[2] é a coordenada x do segundo ponto e assim por diante. Assuma que a
        # quantidade de pontos é sempre multiplo de 3, ou seja, 6 valores ou 12 valores, etc.
        # O parâmetro colors é um dicionário com os tipos cores possíveis, para o TriangleSet2D
        # você pode assumir inicialmente o desenho das linhas com a cor emissiva (emissiveColor).

        def draw_line(x1, y1, x2, y2, color):
            dx = x2 - x1
            dy = y2 - y1

            steps = abs(dx) if abs(dx) > abs(dy) else abs(dy)

            x_inc = dx / steps
            y_inc = dy / steps

            for i in range(int(steps)):
                if 0 <= int(x1) < GL.width and 0 <= int(y1) < GL.height:
                    gpu.GPU.draw_pixel([int(x1), int(y1)], gpu.GPU.RGB8, color)
                x1 += x_inc
                y1 += y_inc
        
        color = [int(255 * colors["emissiveColor"][0]),
                int(255 * colors["emissiveColor"][1]),
                int(255 * colors["emissiveColor"][2])]
        
        for i in range(0, len(vertices), 6):
            draw_line(vertices[i], vertices[i+1], vertices[i+2], vertices[i+3], color)
            draw_line(vertices[i+2], vertices[i+3], vertices[i+4], vertices[i+5], color)
            draw_line(vertices[i+4], vertices[i+5], vertices[i], vertices[i+1], color)

        def L(x, y, x0, y0, x1, y1):
            return (y1 - y0) * x - (x1 - x0) * y + y0 * (x1 - x0) - x0 * (y1 - y0)
        
        def is_inside(x, y, x0, y0, x1, y1, x2, y2):
            return L(x, y, x0, y0, x1, y1) >= 0 and L(x, y, x1, y1, x2, y2) >= 0 and L(x, y, x2, y2, x0, y0) >= 0
        
        for i in range(0, len(vertices), 6):
            x0 = vertices[i]
            y0 = vertices[i+1]
            x1 = vertices[i+2]
            y1 = vertices[i+3]
            x2 = vertices[i+4]
            y2 = vertices[i+5]

            for x in range(GL.width):
                for y in range(GL.height):
                    if is_inside(x, y, x0, y0, x1, y1, x2, y2):
                        gpu.GPU.draw_pixel([x, y], gpu.GPU.RGB8, color)


    @staticmethod
    def triangleSet(point, colors):
        """Função usada para renderizar TriangleSet."""
        color = [int(255 * colors["emissiveColor"][0]),
                 int(255 * colors["emissiveColor"][1]),
                 int(255 * colors["emissiveColor"][2])]

        for i in range(0, len(point), 9):
            p1 = np.array([[point[i]], [point[i+1]], [point[i+2]], [1]])
            p2 = np.array([[point[i+3]], [point[i+4]], [point[i+5]], [1]])
            p3 = np.array([[point[i+6]], [point[i+7]], [point[i+8]], [1]])

            # Aplicar transformações
            transform_matrix = GL.transformation_stack[-1]
            p1_world = transform_matrix @ p1
            p2_world = transform_matrix @ p2
            p3_world = transform_matrix @ p3

            # Aplicar visualização
            p1_view = GL.view_matrix @ p1_world
            p2_view = GL.view_matrix @ p2_world
            p3_view = GL.view_matrix @ p3_world

            # Aplicar projeção
            p1_proj = GL.perspective_matrix @ p1_view
            p2_proj = GL.perspective_matrix @ p2_view
            p3_proj = GL.perspective_matrix @ p3_view

            # Normalizar
            p1_ndc = p1_proj / p1_proj[3][0]
            p2_ndc = p2_proj / p2_proj[3][0]
            p3_ndc = p3_proj / p3_proj[3][0]

            # Converter para coordenadas de tela
            x0 = int((p1_ndc[0][0] + 1) * (GL.width / 2))
            y0 = int((1 - p1_ndc[1][0]) * (GL.height / 2))
            x1 = int((p2_ndc[0][0] + 1) * (GL.width / 2))
            y1 = int((1 - p2_ndc[1][0]) * (GL.height / 2))
            x2 = int((p3_ndc[0][0] + 1) * (GL.width / 2))
            y2 = int((1 - p3_ndc[1][0]) * (GL.height / 2))

            # Desenhar o triângulo
            GL.fill_triangle(x0, y0, x1, y1, x2, y2, color)

    @staticmethod
    def viewpoint(position, orientation, fieldOfView):
        """Função usada para renderizar (na verdade coletar os dados) de Viewpoint."""
        # Na função de viewpoint você receberá a posição, orientação e campo de visão da
        # câmera virtual. Use esses dados para poder calcular e criar a matriz de projeção
        # perspectiva para poder aplicar nos pontos dos objetos geométricos.

        translation_matrix = np.array([
            [1,0,0,-position[0]],
            [0,1,0,-position[1]],
            [0,0,1,-position[2]],
            [0,0,0,1]
        ])

        x, y, z, w = orientation[0], orientation[1], orientation[2], orientation[3]/2
        mag = math.sqrt(x*x + y*y + z*z)
        x, y, z = x/mag, y/mag, z/mag
        escalar = math.cos(w)
        v_x = x*math.sin(w)
        v_y = y*math.sin(w)
        v_z = z*math.sin(w)

        rotation_matrix = np.array([
            [1 - 2*(v_y**2 + v_z**2), 2*(v_x*v_y - v_z*escalar), 2*(v_x*v_z + v_y*escalar), 0],
            [2*(v_x*v_y + v_z*escalar), 1 - 2*(v_x**2 + v_z**2), 2*(v_y*v_z - v_x*escalar), 0],
            [2*(v_x*v_z - v_y*escalar), 2*(v_y*v_z + v_x*escalar), 1 - 2*(v_x**2 + v_y**2), 0],
            [0, 0, 0, 1]
        ])
        inv_rotation_matrix = np.linalg.inv(rotation_matrix)

        view_matrix = inv_rotation_matrix @ translation_matrix

        fovy = 2*math.atan(math.tan(fieldOfView/2)*(GL.height/math.sqrt(GL.height**2 + GL.width**2)))
        top = GL.near*math.tan(fovy)
        right = top*(GL.width/GL.height)
        perspective_matrix = np.array([
            [GL.near/right, 0, 0, 0],
            [0, GL.near/top, 0, 0],
            [0, 0, -(GL.far + GL.near)/(GL.far - GL.near), -2*GL.far*GL.near/(GL.far - GL.near)],
            [0, 0, -1, 0]
        ])
        
        GL.view_matrix = view_matrix
        GL.perspective_matrix = perspective_matrix

    @staticmethod
    def transform_in(translation, scale, rotation):
        """Função usada para renderizar (na verdade coletar os dados) de Transform."""
        # A função transform_in será chamada quando se entrar em um nó X3D do tipo Transform
        # do grafo de cena. Os valores passados são a escala em um vetor [x, y, z]
        # indicando a escala em cada direção, a translação [x, y, z] nas respectivas
        # coordenadas e finalmente a rotação por [x, y, z, t] sendo definida pela rotação
        # do objeto ao redor do eixo x, y, z por t radianos, seguindo a regra da mão direita.
        # Quando se entrar em um nó transform se deverá salvar a matriz de transformação dos
        # modelos do mundo para depois potencialmente usar em outras chamadas. 
        # Quando começar a usar Transforms dentre de outros Transforms, mais a frente no curso
        # Você precisará usar alguma estrutura de dados pilha para organizar as matrizes.

        translation_matrix = np.array([
            [1, 0, 0, translation[0]],
            [0, 1, 0, translation[1]],
            [0, 0, 1, translation[2]],
            [0, 0, 0, 1]
        ])

        scale_matrix = np.array([
            [scale[0], 0, 0, 0],
            [0, scale[1], 0, 0],
            [0, 0, scale[2], 0],
            [0, 0, 0, 1]
        ])
        
        x, y, z, w = rotation[0], rotation[1], rotation[2], rotation[3]/2
        mag = math.sqrt(x*x + y*y + z*z)
        x, y, z = x/mag, y/mag, z/mag
        escalar = math.cos(w)
        v_x = x*math.sin(w)
        v_y = y*math.sin(w)
        v_z = z*math.sin(w)

        rotation_matrix = np.array([
            [1 - 2*(v_y**2 + v_z**2), 2*(v_x*v_y - v_z*escalar), 2*(v_x*v_z + v_y*escalar), 0],
            [2*(v_x*v_y + v_z*escalar), 1 - 2*(v_x**2 + v_z**2), 2*(v_y*v_z - v_x*escalar), 0],
            [2*(v_x*v_z - v_y*escalar), 2*(v_y*v_z + v_x*escalar), 1 - 2*(v_x**2 + v_y**2), 0],
            [0, 0, 0, 1]
        ])
        
        GL.transformation_stack.append(translation_matrix @ rotation_matrix @ scale_matrix)

    @staticmethod
    def transform_out():
        """Função usada para renderizar (na verdade coletar os dados) de Transform."""
        # A função transform_out será chamada quando se sair em um nó X3D do tipo Transform do
        # grafo de cena. Não são passados valores, porém quando se sai de um nó transform se
        # deverá recuperar a matriz de transformação dos modelos do mundo da estrutura de
        # pilha implementada.

        if GL.transformation_stack:
            GL.transformation_stack.pop()

    @staticmethod
    def triangleStripSet(point, stripCount, colors):
        """Função usada para renderizar TriangleStripSet."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/rendering.html#TriangleStripSet
        # A função triangleStripSet é usada para desenhar tiras de triângulos interconectados,
        # você receberá as coordenadas dos pontos no parâmetro point, esses pontos são uma
        # lista de pontos x, y, e z sempre na ordem. Assim point[0] é o valor da coordenada x
        # do primeiro ponto, point[1] o valor y do primeiro ponto, point[2] o valor z da
        # coordenada z do primeiro ponto. Já point[3] é a coordenada x do segundo ponto e assim
        # por diante. No TriangleStripSet a quantidade de vértices a serem usados é informado
        # em uma lista chamada stripCount (perceba que é uma lista). Ligue os vértices na ordem,
        # primeiro triângulo será com os vértices 0, 1 e 2, depois serão os vértices 1, 2 e 3,
        # depois 2, 3 e 4, e assim por diante. Cuidado com a orientação dos vértices, ou seja,
        # todos no sentido horário ou todos no sentido anti-horário, conforme especificado.

        index = 0
        vertices = []

        for count in stripCount:
            num_vertices = count
            if num_vertices < 3:
                index += num_vertices * 3
                continue

            strip_points = point[index : index + num_vertices * 3]
            index += num_vertices * 3

            v0 = strip_points[0:3]
            v1 = strip_points[3:6]

            for i in range(2, num_vertices):
                v2 = strip_points[i*3 : i*3 + 3]

                if (i % 2 == 0):
                    vertices.extend(v0 + v1 + v2)
                else:
                    vertices.extend(v0 + v2 + v1)

                v0 = v1
                v1 = v2

        GL.triangleSet(vertices, colors)

    @staticmethod
    def indexedTriangleStripSet(point, index, colors):
        """Função usada para renderizar IndexedTriangleStripSet."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/rendering.html#IndexedTriangleStripSet
        # A função indexedTriangleStripSet é usada para desenhar tiras de triângulos
        # interconectados, você receberá as coordenadas dos pontos no parâmetro point, esses
        # pontos são uma lista de pontos x, y, e z sempre na ordem. Assim point[0] é o valor
        # da coordenada x do primeiro ponto, point[1] o valor y do primeiro ponto, point[2]
        # o valor z da coordenada z do primeiro ponto. Já point[3] é a coordenada x do
        # segundo ponto e assim por diante. No IndexedTriangleStripSet uma lista informando
        # como conectar os vértices é informada em index, o valor -1 indica que a lista
        # acabou. A ordem de conexão será de 3 em 3 pulando um índice. Por exemplo: o
        # primeiro triângulo será com os vértices 0, 1 e 2, depois serão os vértices 1, 2 e 3,
        # depois 2, 3 e 4, e assim por diante. Cuidado com a orientação dos vértices, ou seja,
        # todos no sentido horário ou todos no sentido anti-horário, conforme especificado.

        vertices = []
        strip_indices = []
        current_strip = []

        for i in index:
            if i == -1:
                if len(current_strip) >= 3:
                    strip_indices.append(current_strip)
                current_strip = []
            else:
                current_strip.append(i)

        if len(current_strip) >= 3:
            strip_indices.append(current_strip)

        for strip in strip_indices:
            num_vertices = len(strip)
            if num_vertices < 3:
                continue

            v0 = strip[0]
            v1 = strip[1]

            for j in range(2, num_vertices):
                v2 = strip[j]

                if (j % 2 == 0):
                    vertices.extend(point[v0*3 : v0*3+3] + point[v1*3 : v1*3+3] + point[v2*3 : v2*3+3])
                else:
                    vertices.extend(point[v0*3 : v0*3+3] + point[v2*3 : v2*3+3] + point[v1*3 : v1*3+3])

                v0 = v1
                v1 = v2

        # Chamar a função triangleSet para renderizar os triângulos
        GL.triangleSet(vertices, colors)

    @staticmethod
    def indexedFaceSet(coord, coordIndex, colorPerVertex, color, colorIndex,
                       texCoord, texCoordIndex, colors, current_texture):
        """Função usada para renderizar IndexedFaceSet."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#IndexedFaceSet
        # A função indexedFaceSet é usada para desenhar malhas de triângulos. Ela funciona de
        # forma muito simular a IndexedTriangleStripSet porém com mais recursos.
        # Você receberá as coordenadas dos pontos no parâmetro cord, esses
        # pontos são uma lista de pontos x, y, e z sempre na ordem. Assim coord[0] é o valor
        # da coordenada x do primeiro ponto, coord[1] o valor y do primeiro ponto, coord[2]
        # o valor z da coordenada z do primeiro ponto. Já coord[3] é a coordenada x do
        # segundo ponto e assim por diante. No IndexedFaceSet uma lista de vértices é informada
        # em coordIndex, o valor -1 indica que a lista acabou.
        # A ordem de conexão não possui uma ordem oficial, mas em geral se o primeiro ponto com os dois
        # seguintes e depois este mesmo primeiro ponto com o terçeiro e quarto ponto. Por exemplo: numa
        # sequencia 0, 1, 2, 3, 4, -1 o primeiro triângulo será com os vértices 0, 1 e 2, depois serão
        # os vértices 0, 2 e 3, e depois 0, 3 e 4, e assim por diante, até chegar no final da lista.
        # Adicionalmente essa implementação do IndexedFace aceita cores por vértices, assim
        # se a flag colorPerVertex estiver habilitada, os vértices também possuirão cores
        # que servem para definir a cor interna dos poligonos, para isso faça um cálculo
        # baricêntrico de que cor deverá ter aquela posição. Da mesma forma se pode definir uma
        # textura para o poligono, para isso, use as coordenadas de textura e depois aplique a
        # cor da textura conforme a posição do mapeamento. Dentro da classe GPU já está
        # implementadado um método para a leitura de imagens.

        vertices = []
        i = 0
        face = []

        while i < len(coordIndex):
            index = coordIndex[i]
            if index == -1:
                if len(face) >= 3:
                    v0 = face[0]
                    for j in range(1, len(face)-1):
                        v1 = face[j]
                        v2 = face[j+1]
                        # Adicionar os vértices do triângulo à lista
                        vertices.extend(coord[v0*3 : v0*3+3] +
                                        coord[v1*3 : v1*3+3] +
                                        coord[v2*3 : v2*3+3])
                face = []
            else:
                face.append(index)
            i += 1

        if len(face) >= 3:
            v0 = face[0]
            for j in range(1, len(face)-1):
                v1 = face[j]
                v2 = face[j+1]
                vertices.extend(coord[v0*3 : v0*3+3] +
                                coord[v1*3 : v1*3+3] +
                                coord[v2*3 : v2*3+3])

        GL.triangleSet(vertices, colors)

    @staticmethod
    def fill_triangle(x0, y0, x1, y1, x2, y2, color):
        """Preenche um triângulo na tela usando o algoritmo de varredura."""
        # Ordenar os pontos pelo y
        if y0 > y1:
            x0, y0, x1, y1 = x1, y1, x0, y0
        if y0 > y2:
            x0, y0, x2, y2 = x2, y2, x0, y0
        if y1 > y2:
            x1, y1, x2, y2 = x2, y2, x1, y1

        # Função auxiliar para desenhar linhas
        def draw_line(xa, ya, xb, yb):
            dx = xb - xa
            dy = yb - ya
            steps = max(abs(dx), abs(dy))
            if steps == 0:
                gpu.GPU.draw_pixel([int(xa), int(ya)], gpu.GPU.RGB8, color)
                return
            x_inc = dx / steps
            y_inc = dy / steps
            x = xa
            y = ya
            for _ in range(int(steps)+1):
                gpu.GPU.draw_pixel([int(x), int(y)], gpu.GPU.RGB8, color)
                x += x_inc
                y += y_inc

        # Preencher o triângulo usando varredura de linhas
        total_height = y2 - y0
        if total_height == 0:
            return
        for y in range(int(y0), int(y2)+1):
            segment_height = y2 - y1 if y > y1 else y1 - y0
            alpha = (y - y0) / total_height if total_height != 0 else 0
            beta = (y - y0) / (y1 - y0) if y < y1 and y1 - y0 != 0 else (y - y1) / (y2 - y1) if y1 != y2 and y >= y1 else 0

            xa = x0 + (x2 - x0) * alpha
            xb = x0 + (x1 - x0) * beta if y < y1 else x1 + (x2 - x1) * beta

            if xa > xb:
                xa, xb = xb, xa

            for x in range(int(xa), int(xb)+1):
                gpu.GPU.draw_pixel([x, y], gpu.GPU.RGB8, color)


    @staticmethod
    def box(size, colors):
        """Função usada para renderizar Boxes."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#Box
        # A função box é usada para desenhar paralelepípedos na cena. O Box é centrada no
        # (0, 0, 0) no sistema de coordenadas local e alinhado com os eixos de coordenadas
        # locais. O argumento size especifica as extensões da caixa ao longo dos eixos X, Y
        # e Z, respectivamente, e cada valor do tamanho deve ser maior que zero. Para desenha
        # essa caixa você vai provavelmente querer tesselar ela em triângulos, para isso
        # encontre os vértices e defina os triângulos.
        hx = size[0] / 2
        hy = size[1] / 2
        hz = size[2] / 2

        # Definir os 8 vértices do paralelepípedo
        vertices = [
            [-hx, -hy, -hz],  # v0
            [hx, -hy, -hz],   # v1
            [hx, hy, -hz],    # v2
            [-hx, hy, -hz],   # v3
            [-hx, -hy, hz],   # v4
            [hx, -hy, hz],    # v5
            [hx, hy, hz],     # v6
            [-hx, hy, hz]     # v7
        ]

        # Definir as faces laterais com a ordem dos vértices ajustada
        faces = [
            [4, 5, 6, 7],  # Frente
            [0, 1, 2, 3],  # Traseira
            [0, 4, 7, 3],  # Esquerda
            [1, 5, 6, 2]   # Direita
        ]

        # Coletar os triângulos
        triangles = []

        for face in faces:
            v0 = vertices[face[0]]
            v1 = vertices[face[1]]
            v2 = vertices[face[2]]
            v3 = vertices[face[3]]

            # Primeiro triângulo
            triangles.extend(v0 + v1 + v2)
            # Segundo triângulo
            triangles.extend(v0 + v2 + v3)

        GL.triangleSet(triangles, colors)

    @staticmethod
    def sphere(radius, colors):
        """Função usada para renderizar Esferas."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#Sphere
        # A função sphere é usada para desenhar esferas na cena. O esfera é centrada no
        # (0, 0, 0) no sistema de coordenadas local. O argumento radius especifica o
        # raio da esfera que está sendo criada. Para desenha essa esfera você vai
        # precisar tesselar ela em triângulos, para isso encontre os vértices e defina
        # os triângulos.

        latitude_bands = 50
        longitude_bands = 50

        vertices = []

        for lat_number in range(latitude_bands + 1):
            theta = lat_number * math.pi / latitude_bands
            sin_theta = math.sin(theta)
            cos_theta = math.cos(theta)

            for long_number in range(longitude_bands + 1):
                phi = long_number * 2 * math.pi / longitude_bands
                sin_phi = math.sin(phi)
                cos_phi = math.cos(phi)

                x = radius * sin_theta * cos_phi
                y = radius * cos_theta
                z = radius * sin_theta * sin_phi

                vertices.append([x, y, z])

        # Coletar triângulos
        indices = []

        for lat_number in range(latitude_bands):
            for long_number in range(longitude_bands):
                first = (lat_number * (longitude_bands + 1)) + long_number
                second = first + longitude_bands + 1

                indices.extend([
                    first, second, first + 1,
                    second, second + 1, first + 1
                ])

        # Converter vértices e índices em lista de coordenadas para triangleSet
        triangles = []
        for i in range(0, len(indices), 3):
            v0 = vertices[indices[i]]
            v1 = vertices[indices[i + 1]]
            v2 = vertices[indices[i + 2]]
            triangles.extend(v0 + v1 + v2)

        GL.triangleSet(triangles, colors)

    @staticmethod
    def cone(bottomRadius, height, colors):
        """Função usada para renderizar Cones."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#Cone
        # A função cone é usada para desenhar cones na cena. O cone é centrado no
        # (0, 0, 0) no sistema de coordenadas local. O argumento bottomRadius especifica o
        # raio da base do cone e o argumento height especifica a altura do cone.
        # O cone é alinhado com o eixo Y local. O cone é fechado por padrão na base.
        # Para desenha esse cone você vai precisar tesselar ele em triângulos, para isso
        # encontre os vértices e defina os triângulos.

        N = 20  # Número de segmentos
        angle_step = 2 * math.pi / N

        # Gerar pontos da base
        base_points = []
        for i in range(N):
            angle = i * angle_step
            x = bottomRadius * math.cos(angle)
            z = bottomRadius * math.sin(angle)
            y = 0  # Base em y=0
            base_points.append([x, y, z])

        # Ponto do ápice
        apex = [0, height, 0]

        # Coletar triângulos apenas para as faces laterais
        triangles = []

        for i in range(N):
            v0 = base_points[i]
            v1 = apex
            v2 = base_points[(i+1) % N]
            triangles.extend(v0 + v1 + v2)

        GL.triangleSet(triangles, colors)

    @staticmethod
    def cylinder(radius, height, colors):
        """Função usada para renderizar Cilindros."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#Cylinder
        # A função cylinder é usada para desenhar cilindros na cena. O cilindro é centrado no
        # (0, 0, 0) no sistema de coordenadas local. O argumento radius especifica o
        # raio da base do cilindro e o argumento height especifica a altura do cilindro.
        # O cilindro é alinhado com o eixo Y local. O cilindro é fechado por padrão em ambas as extremidades.
        # Para desenha esse cilindro você vai precisar tesselar ele em triângulos, para isso
        # encontre os vértices e defina os triângulos.

        N = 1000  # Número de segmentos
        angle_step = 2 * math.pi / N

        # Gerar pontos da base e do topo
        base_points = []
        top_points = []

        for i in range(N + 1):
            angle = i * angle_step
            x = radius * math.cos(angle)
            z = radius * math.sin(angle)
            base_points.append([x, 0, z])
            top_points.append([x, height, z])

        # Coletar triângulos para a superfície lateral
        triangles = []

        for i in range(N):
            # Vértices do quadrilátero lateral
            p0 = base_points[i]
            p1 = top_points[i]
            p2 = base_points[i + 1]
            p3 = top_points[i + 1]

            # Primeiro triângulo
            triangles.extend(p0 + p2 + p1)
            # Segundo triângulo
            triangles.extend(p1 + p2 + p3)

        GL.triangleSet(triangles, colors)

    @staticmethod
    def navigationInfo(headlight):
        """Características físicas do avatar do visualizador e do modelo de visualização."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/navigation.html#NavigationInfo
        # O campo do headlight especifica se um navegador deve acender um luz direcional que
        # sempre aponta na direção que o usuário está olhando. Definir este campo como TRUE
        # faz com que o visualizador forneça sempre uma luz do ponto de vista do usuário.
        # A luz headlight deve ser direcional, ter intensidade = 1, cor = (1 1 1),
        # ambientIntensity = 0,0 e direção = (0 0 −1).

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("NavigationInfo : headlight = {0}".format(headlight)) # imprime no terminal

    @staticmethod
    def directionalLight(ambientIntensity, color, intensity, direction):
        """Luz direcional ou paralela."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/lighting.html#DirectionalLight
        # Define uma fonte de luz direcional que ilumina ao longo de raios paralelos
        # em um determinado vetor tridimensional. Possui os campos básicos ambientIntensity,
        # cor, intensidade. O campo de direção especifica o vetor de direção da iluminação
        # que emana da fonte de luz no sistema de coordenadas local. A luz é emitida ao
        # longo de raios paralelos de uma distância infinita.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("DirectionalLight : ambientIntensity = {0}".format(ambientIntensity))
        print("DirectionalLight : color = {0}".format(color)) # imprime no terminal
        print("DirectionalLight : intensity = {0}".format(intensity)) # imprime no terminal
        print("DirectionalLight : direction = {0}".format(direction)) # imprime no terminal

    @staticmethod
    def pointLight(ambientIntensity, color, intensity, location):
        """Luz pontual."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/lighting.html#PointLight
        # Fonte de luz pontual em um local 3D no sistema de coordenadas local. Uma fonte
        # de luz pontual emite luz igualmente em todas as direções; ou seja, é omnidirecional.
        # Possui os campos básicos ambientIntensity, cor, intensidade. Um nó PointLight ilumina
        # a geometria em um raio de sua localização. O campo do raio deve ser maior ou igual a
        # zero. A iluminação do nó PointLight diminui com a distância especificada.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("PointLight : ambientIntensity = {0}".format(ambientIntensity))
        print("PointLight : color = {0}".format(color)) # imprime no terminal
        print("PointLight : intensity = {0}".format(intensity)) # imprime no terminal
        print("PointLight : location = {0}".format(location)) # imprime no terminal

    @staticmethod
    def fog(visibilityRange, color):
        """Névoa."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/environmentalEffects.html#Fog
        # O nó Fog fornece uma maneira de simular efeitos atmosféricos combinando objetos
        # com a cor especificada pelo campo de cores com base nas distâncias dos
        # vários objetos ao visualizador. A visibilidadeRange especifica a distância no
        # sistema de coordenadas local na qual os objetos são totalmente obscurecidos
        # pela névoa. Os objetos localizados fora de visibilityRange do visualizador são
        # desenhados com uma cor de cor constante. Objetos muito próximos do visualizador
        # são muito pouco misturados com a cor do nevoeiro.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Fog : color = {0}".format(color)) # imprime no terminal
        print("Fog : visibilityRange = {0}".format(visibilityRange))

    @staticmethod
    def timeSensor(cycleInterval, loop):
        """Gera eventos conforme o tempo passa."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/time.html#TimeSensor
        # Os nós TimeSensor podem ser usados para muitas finalidades, incluindo:
        # Condução de simulações e animações contínuas; Controlar atividades periódicas;
        # iniciar eventos de ocorrência única, como um despertador;
        # Se, no final de um ciclo, o valor do loop for FALSE, a execução é encerrada.
        # Por outro lado, se o loop for TRUE no final de um ciclo, um nó dependente do
        # tempo continua a execução no próximo ciclo. O ciclo de um nó TimeSensor dura
        # cycleInterval segundos. O valor de cycleInterval deve ser maior que zero.

        # Deve retornar a fração de tempo passada em fraction_changed

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("TimeSensor : cycleInterval = {0}".format(cycleInterval)) # imprime no terminal
        print("TimeSensor : loop = {0}".format(loop))

        # Esse método já está implementado para os alunos como exemplo
        epoch = time.time()  # time in seconds since the epoch as a floating point number.
        fraction_changed = (epoch % cycleInterval) / cycleInterval

        return fraction_changed

    @staticmethod
    def splinePositionInterpolator(set_fraction, key, keyValue, closed):
        """Interpola não linearmente entre uma lista de vetores 3D."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/interpolators.html#SplinePositionInterpolator
        # Interpola não linearmente entre uma lista de vetores 3D. O campo keyValue possui
        # uma lista com os valores a serem interpolados, key possui uma lista respectiva de chaves
        # dos valores em keyValue, a fração a ser interpolada vem de set_fraction que varia de
        # zeroa a um. O campo keyValue deve conter exatamente tantos vetores 3D quanto os
        # quadros-chave no key. O campo closed especifica se o interpolador deve tratar a malha
        # como fechada, com uma transições da última chave para a primeira chave. Se os keyValues
        # na primeira e na última chave não forem idênticos, o campo closed será ignorado.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("SplinePositionInterpolator : set_fraction = {0}".format(set_fraction))
        print("SplinePositionInterpolator : key = {0}".format(key)) # imprime no terminal
        print("SplinePositionInterpolator : keyValue = {0}".format(keyValue))
        print("SplinePositionInterpolator : closed = {0}".format(closed))

        # Abaixo está só um exemplo de como os dados podem ser calculados e transferidos
        value_changed = [0.0, 0.0, 0.0]
        
        return value_changed

    @staticmethod
    def orientationInterpolator(set_fraction, key, keyValue):
        """Interpola entre uma lista de valores de rotação especificos."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/interpolators.html#OrientationInterpolator
        # Interpola rotações são absolutas no espaço do objeto e, portanto, não são cumulativas.
        # Uma orientação representa a posição final de um objeto após a aplicação de uma rotação.
        # Um OrientationInterpolator interpola entre duas orientações calculando o caminho mais
        # curto na esfera unitária entre as duas orientações. A interpolação é linear em
        # comprimento de arco ao longo deste caminho. Os resultados são indefinidos se as duas
        # orientações forem diagonalmente opostas. O campo keyValue possui uma lista com os
        # valores a serem interpolados, key possui uma lista respectiva de chaves
        # dos valores em keyValue, a fração a ser interpolada vem de set_fraction que varia de
        # zeroa a um. O campo keyValue deve conter exatamente tantas rotações 3D quanto os
        # quadros-chave no key.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("OrientationInterpolator : set_fraction = {0}".format(set_fraction))
        print("OrientationInterpolator : key = {0}".format(key)) # imprime no terminal
        print("OrientationInterpolator : keyValue = {0}".format(keyValue))

        # Abaixo está só um exemplo de como os dados podem ser calculados e transferidos
        value_changed = [0, 0, 1, 0]

        return value_changed

    # Para o futuro (Não para versão atual do projeto.)
    def vertex_shader(self, shader):
        """Para no futuro implementar um vertex shader."""

    def fragment_shader(self, shader):
        """Para no futuro implementar um fragment shader."""
