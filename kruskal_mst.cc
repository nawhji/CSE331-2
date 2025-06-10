#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
using namespace std;

struct Edge {
    int u, v;
    double w;
    bool operator<(const Edge& other) const {
        return w < other.w;
    }
};

struct Path {
    vector<int> nodes;
};

void connect_remaining_paths(unordered_map<int, Path>& node_to_path,
                             unordered_map<int, int>& node_to_path_id,
                             double** matrix,
                             int& next_path_id,
                             double& total_cost) {
    vector<Edge> extra_edges;

    vector<int> path_ids;
    for (auto& [id, _] : node_to_path) path_ids.push_back(id);

    for (size_t i = 0; i < path_ids.size(); ++i) {
        for (size_t j = i + 1; j < path_ids.size(); ++j) {
            int id1 = path_ids[i], id2 = path_ids[j];
            Path& p1 = node_to_path[id1];
            Path& p2 = node_to_path[id2];

            vector<pair<int, int>> pairs = {
                {p1.nodes.front(), p2.nodes.front()},
                {p1.nodes.front(), p2.nodes.back()},
                {p1.nodes.back(),  p2.nodes.front()},
                {p1.nodes.back(),  p2.nodes.back()}
            };

            for (auto [u, v] : pairs) {
                extra_edges.push_back({u, v, matrix[u][v]});
            }
        }
    }

    sort(extra_edges.begin(), extra_edges.end());

    for (auto& e : extra_edges) {
        int u = e.u, v = e.v;
        int pu = node_to_path_id[u];
        int pv = node_to_path_id[v];
        if (pu == pv) continue;

        Path& path_u = node_to_path[pu];
        Path& path_v = node_to_path[pv];

        bool u_head = (path_u.nodes.front() == u);
        bool u_tail = (path_u.nodes.back() == u);
        bool v_head = (path_v.nodes.front() == v);
        bool v_tail = (path_v.nodes.back() == v);

        vector<int> new_nodes;
        if (u_tail && v_head) {
            new_nodes = path_u.nodes;
            new_nodes.insert(new_nodes.end(), path_v.nodes.begin(), path_v.nodes.end());
        } else if (v_tail && u_head) {
            new_nodes = path_v.nodes;
            new_nodes.insert(new_nodes.end(), path_u.nodes.begin(), path_u.nodes.end());
        } else if (u_tail && v_tail) {
            new_nodes = path_u.nodes;
            reverse(path_v.nodes.begin(), path_v.nodes.end());
            new_nodes.insert(new_nodes.end(), path_v.nodes.begin(), path_v.nodes.end());
        } else if (u_head && v_head) {
            reverse(path_u.nodes.begin(), path_u.nodes.end());
            new_nodes = path_u.nodes;
            new_nodes.insert(new_nodes.end(), path_v.nodes.begin(), path_v.nodes.end());
        } else {
            continue;
        }

        total_cost += e.w;
        int new_id = next_path_id++;
        node_to_path[new_id] = Path{new_nodes};

        for (int node : new_nodes)
            node_to_path_id[node] = new_id;

        node_to_path.erase(pu);
        node_to_path.erase(pv);

        if (node_to_path.size() == 1) {
            auto it = node_to_path.begin();
            Path& final_path = it->second;
            int first = final_path.nodes.front();
            int last = final_path.nodes.back();
            total_cost += matrix[last][first];
        }
    }
}

double compute_distance_limit(double** matrix, int dim, double percentile) {
    vector<double> distances;

    for (int i = 0; i < dim; ++i) {
        for (int j = i + 1; j < dim; ++j) {
            distances.push_back(matrix[i][j]);
        }
    }

    sort(distances.begin(), distances.end());

    int index = static_cast<int>(percentile * distances.size());
    return distances[index];
}

double kruskal_mst(double** matrix, int dim) {
    vector<Edge> edges;

    double dist_limit = compute_distance_limit(matrix, dim, 0.05);

    for (int i = 0; i < dim; ++i) {
        for (int j = i + 1; j < dim; ++j) {
            if (matrix[i][j] <= dist_limit)
                edges.push_back({i, j, matrix[i][j]});
        }
    }

    sort(edges.begin(), edges.end());

    unordered_map<int, Path> node_to_path;
    unordered_map<int, int> node_to_path_id;
    int next_path_id = 0;

    for (int i = 0; i < dim; ++i) {
        Path p;
        p.nodes.push_back(i);
        node_to_path[next_path_id] = p;
        node_to_path_id[i] = next_path_id;
        next_path_id++;
    }

    double total_cost = 0.0;

    for (const auto& e : edges) {
        int u = e.u, v = e.v;
        int path_u = node_to_path_id[u];
        int path_v = node_to_path_id[v];

        if (path_u == path_v) continue;

        Path& pu = node_to_path[path_u];
        Path& pv = node_to_path[path_v];

        bool u_head = (pu.nodes.front() == u);
        bool u_tail = (pu.nodes.back() == u);
        bool v_head = (pv.nodes.front() == v);
        bool v_tail = (pv.nodes.back() == v);

        if (!(u_head || u_tail) || !(v_head || v_tail)) continue;

        vector<int> new_nodes;
        if (u_tail && v_head) {
            new_nodes = pu.nodes;
            new_nodes.insert(new_nodes.end(), pv.nodes.begin(), pv.nodes.end());
        } else if (v_tail && u_head) {
            new_nodes = pv.nodes;
            new_nodes.insert(new_nodes.end(), pu.nodes.begin(), pu.nodes.end());
        } else if (u_tail && v_tail) {
            new_nodes = pu.nodes;
            reverse(pv.nodes.begin(), pv.nodes.end());
            new_nodes.insert(new_nodes.end(), pv.nodes.begin(), pv.nodes.end());
        } else if (u_head && v_head) {
            reverse(pu.nodes.begin(), pu.nodes.end());
            new_nodes = pu.nodes;
            new_nodes.insert(new_nodes.end(), pv.nodes.begin(), pv.nodes.end());
        }

        Path new_path;
        new_path.nodes = new_nodes;
        int new_id = next_path_id++;
        node_to_path[new_id] = new_path;

        for (int node : new_nodes) {
            node_to_path_id[node] = new_id;
        }

        node_to_path.erase(path_u);
        node_to_path.erase(path_v);

        total_cost += e.w;

        int single_count = 0;
        for (auto& [id, path] : node_to_path) {
            if (path.nodes.size() == 1) single_count++;
        }

        if (node_to_path.size() > 1) {
            connect_remaining_paths(node_to_path, node_to_path_id, matrix, next_path_id, total_cost);
            break;
        }
    }

    return total_cost;
}