#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <queue>
#include <fstream>
using namespace std;

namespace kruskal_ondemand {

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

double euclidean_dist(const pair<int, int>& a, const pair<int, int>& b) {
    int dx = a.first - b.first;
    int dy = a.second - b.second;
    return sqrt(dx * dx + dy * dy);
}

double compute_distance_limit(const vector<pair<int, int>>& coords, double percentile = 0.1) {
    size_t n = coords.size();
    uint64_t total_pairs = n * (n - 1) / 2;
    uint64_t MAX_DISTANCES = total_pairs * percentile;  // 1%

    vector<double> distances;
    distances.reserve(MAX_DISTANCES + 10000);  // 여유

    uint64_t count = 0;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            double d = euclidean_dist(coords[i], coords[j]);
            distances.push_back(d);

            // 초과 시 상위 10% 제거
            if (distances.size() > MAX_DISTANCES) {
                size_t keep_size = MAX_DISTANCES * 90 / 100;
                nth_element(distances.begin(), distances.begin() + keep_size, distances.end());
                double cutoff = distances[keep_size];

                vector<double> filtered;
                filtered.reserve(keep_size);
                for (double x : distances)
                    if (x <= cutoff) filtered.push_back(x);

                distances.swap(filtered);
            }

            if (++count % 10000000 == 0)
                cout << "[INFO] Processed " << count << " pairs, current vector size: " << distances.size() << endl;
        }
    }

    return *max_element(distances.begin(), distances.end());
}


double compute_distance_cutoff(const vector<pair<int, int>>& coords, double percentile = 0.1) {
    using namespace std;
    size_t n = coords.size();
    uint64_t total_pairs = n * (n - 1) / 2;
    uint64_t max_keep = total_pairs * percentile;

    // 최대 힙 유지 (작은 값만 남기기 위해 큰 값은 빼버림)
    priority_queue<double> max_heap;

    uint64_t count = 0;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            double d = euclidean_dist(coords[i], coords[j]);
            if (max_heap.size() < max_keep)
                max_heap.push(d);
            else if (d < max_heap.top()) {
                max_heap.pop();
                max_heap.push(d);
            }

            if (++count % 10000000 == 0)
                cout << "[INFO] Processed " << count << " pairs" << endl;
        }
    }

    return max_heap.top();  // top이 가장 큰 값 = cutoff
}

vector<Edge> collect_edges_below_cutoff(const vector<pair<int, int>>& coords, double cutoff) {
    vector<Edge> result;
    size_t n = coords.size();
    for (size_t i = 0; i < n; ++i)
        for (size_t j = i + 1; j < n; ++j) {
            double d = euclidean_dist(coords[i], coords[j]);
            if (d <= cutoff)
                result.push_back({(int)i, (int)j, d});
        }
    return result;
}


void connect_remaining_paths(unordered_map<int, Path>& node_to_path,
                             unordered_map<int, int>& node_to_path_id,
                             const vector<pair<int, int>>& coords,
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
                double w = euclidean_dist(coords[u], coords[v]);
                extra_edges.push_back({u, v, w});
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
        } else continue;

        total_cost += e.w;
        int new_id = next_path_id++;
        node_to_path[new_id] = Path{new_nodes};
        for (int node : new_nodes)
            node_to_path_id[node] = new_id;

        node_to_path.erase(pu);
        node_to_path.erase(pv);

        if (node_to_path.size() == 1) {
            int first = new_nodes.front();
            int last = new_nodes.back();
            total_cost += euclidean_dist(coords[first], coords[last]);
        }
    }
}

void save_result(unordered_map<int, Path>& node_to_path, const vector<pair<int, int>>& coords) {
    ofstream fout("monalisa_result.txt");
    for (const auto& [id, path] : node_to_path) {
        const vector<int>& nodes = path.nodes;
        for (size_t i = 0; i + 1 < nodes.size(); ++i) {
            fout << nodes[i] << " " << nodes[i+1] << "\n";
        }
        if (!nodes.empty())
            fout << nodes.back() << " " << nodes.front() << "\n";
    }
    fout.close();

    ofstream fcoord("monalisa_coords_result.txt");
    for (const auto& [x,y] : coords)
        fcoord << x << " " << y << "\n";
    fcoord.close();
}

double kruskal_mst_ondemand(const vector<pair<int, int>>& coords) {
    int dim = coords.size();

    double cutoff = compute_distance_cutoff(coords, 0.1);
    auto edges = collect_edges_below_cutoff(coords, cutoff);

    cout << "[INFO] Distance cutoff: " << cutoff << endl;
    cout << "[INFO] Number of edges after filtering: " << edges.size() << endl;

    sort(edges.begin(), edges.end());

    unordered_map<int, Path> node_to_path;
    unordered_map<int, int> node_to_path_id;
    int next_path_id = 0;

    for (int i = 0; i < dim; ++i) {
        node_to_path[next_path_id] = Path{{i}};
        node_to_path_id[i] = next_path_id++;
    }

    double total_cost = 0.0;
    for (auto& e : edges) {
        int u = e.u, v = e.v;
        int path_u = node_to_path_id[u];
        int path_v = node_to_path_id[v];
        if (path_u == path_v) continue;

        Path& pu = node_to_path[path_u];
        Path& pv = node_to_path[path_v];
        bool u_head = pu.nodes.front() == u, u_tail = pu.nodes.back() == u;
        bool v_head = pv.nodes.front() == v, v_tail = pv.nodes.back() == v;

        if (!(u_head || u_tail) || !(v_head || v_tail)) continue;

        vector<int> new_nodes;
        if (u_tail && v_head)
            new_nodes = pu.nodes, new_nodes.insert(new_nodes.end(), pv.nodes.begin(), pv.nodes.end());
        else if (v_tail && u_head)
            new_nodes = pv.nodes, new_nodes.insert(new_nodes.end(), pu.nodes.begin(), pu.nodes.end());
        else if (u_tail && v_tail) {
            new_nodes = pu.nodes;
            reverse(pv.nodes.begin(), pv.nodes.end());
            new_nodes.insert(new_nodes.end(), pv.nodes.begin(), pv.nodes.end());
        } else if (u_head && v_head) {
            reverse(pu.nodes.begin(), pu.nodes.end());
            new_nodes = pu.nodes;
            new_nodes.insert(new_nodes.end(), pv.nodes.begin(), pv.nodes.end());
        }

        node_to_path[next_path_id] = Path{new_nodes};
        for (int node : new_nodes)
            node_to_path_id[node] = next_path_id;
        next_path_id++;

        node_to_path.erase(path_u);
        node_to_path.erase(path_v);

        total_cost += e.w;

        int single_count = 0;
        for (auto& [_, p] : node_to_path)
            if (p.nodes.size() == 1) single_count++;

        if (single_count == 0) {
            cout << "[INFO] Connecting " << node_to_path.size() << " remaining paths..." << endl;
            connect_remaining_paths(node_to_path, node_to_path_id, coords, next_path_id, total_cost);
            break;
        }
    }

    save_result(node_to_path, coords);
    return total_cost;
}

}