import matplotlib.pyplot as plt

def read_coords(filename):
    coords = []
    with open(filename, 'r') as f:
        for line in f:
            x, y = map(int, line.strip().split())
            coords.append((x, y))
    return coords

def read_edges(filename):
    edges = []
    with open(filename, 'r') as f:
        for line in f:
            u, v = map(int, line.strip().split())
            edges.append((u, v))
    return edges

def plot_edges(coords, edges):
    xs, ys = zip(*coords)
    plt.figure(figsize=(10, 10))
    plt.scatter(xs, ys, s=0.5, color='black')

    for u, v in edges:
        x1, y1 = coords[u]
        x2, y2 = coords[v]
        plt.plot([x1, x2], [y1, y2], color='blue', linewidth=0.4)

    plt.title("Kruskal MST Tour (Approximate TSP)")
    plt.axis("equal")
    plt.show()

coords = read_coords("monalisa_coords_result.txt")
edges = read_edges("monalisa_result.txt")
plot_edges(coords, edges)