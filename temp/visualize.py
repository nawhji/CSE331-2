import matplotlib.pyplot as plt
import math

def rotate_coords(coords, angle_deg):
    angle = math.radians(angle_deg)
    cos_a, sin_a = math.cos(angle), math.sin(angle)
    return [(x*cos_a - y*sin_a, x*sin_a + y*cos_a) for x, y in coords]

def read_tsp_coords(filename):
    coords = []
    with open(filename, 'r') as f:
        found = False
        for line in f:
            if "NODE_COORD_SECTION" in line:
                found = True
                continue
            if not found:
                continue
            if "EOF" in line or line.strip() == "":
                break
            parts = line.strip().split()
            x, y = float(parts[1]), float(parts[2])
            coords.append((x, y))
    coords = rotate_coords(coords, 270)
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
    plt.figure(figsize=(20, 20))
    plt.scatter(xs, ys, c='black', s=5)  # 점 크기 줄임

    for u, v in edges:
        x1, y1 = coords[u]
        x2, y2 = coords[v]
        plt.plot([x1, x2], [y1, y2], 'r-', linewidth=0.8)

    plt.gca().invert_yaxis()
    plt.axis('equal')
    plt.axis('off')  # 축도 안 보이게
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    coords = read_tsp_coords("kz9976.tsp")
    edges = read_edges("edges_kz.txt")
    plot_edges(coords, edges)
