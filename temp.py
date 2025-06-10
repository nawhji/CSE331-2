import itertools
import math

# 정점 좌표 (1-based → 0-based index로 저장)
coords = [
    (6, 7),
    (19, 2),
    (14, 1),
    (10, 11),
    (7, 5),
    (6, 1),
    (18, 0),
    (10, 11),
    (10, 11),
    (3, 16)
]

# 거리 행렬 계산
def euclidean(p1, p2):
    return math.hypot(p1[0] - p2[0], p1[1] - p2[1])

n = len(coords)
dist = [[0] * n for _ in range(n)]
for i in range(n):
    for j in range(n):
        dist[i][j] = euclidean(coords[i], coords[j])

# 브루트포스 TSP
min_cost = float('inf')
best_path = []

for perm in itertools.permutations(range(1, n)):  # 0번 정점 고정
    path = [0] + list(perm)
    cost = sum(dist[path[i]][path[i + 1]] for i in range(n - 1))
    cost += dist[path[-1]][path[0]]  # 돌아오기
    if cost < min_cost:
        min_cost = cost
        best_path = path

# 출력
print(f"Best path: {[x+1 for x in best_path]} → {best_path[0]+1}")
print(f"Cost: {min_cost:.6f}")
