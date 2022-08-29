# The Shoelace Algorithm - www.101computing.net/the-shoelace-algorithm


def polygonArea(vertices):
    # A function to apply the Shoelace algorithm
    numberOfVertices = len(vertices)
    sum1 = 0
    sum2 = 0

    for i in range(0, numberOfVertices - 1):
        sum1 = sum1 + vertices[i][0] * vertices[i + 1][1]
        sum2 = sum2 + vertices[i][1] * vertices[i + 1][0]

    # Add xn.y1
    sum1 = sum1 + vertices[numberOfVertices - 1][0] * vertices[0][1]
    # Add x1.yn
    sum2 = sum2 + vertices[0][0] * vertices[numberOfVertices - 1][1]

    area = abs(sum1 - sum2) / 2
    return area


# Vertices (x,y) Coordinates
A = [2, 7]
B = [10, 1]
C = [8, 6]
D = [11, 7]
E = [7, 10]
# Define a polygon as being a list of vertices, (on anticlockwise order)
polygon = [A, B, C, D, E]

area = polygonArea(polygon)
print("Polygon Vertices:")
print(polygon)
print("")
print("Area = " + str(area) + "cm2")
