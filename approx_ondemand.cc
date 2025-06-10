#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <limits>
#include <vector>

namespace approx_2_ondemand {

    struct Edge {
        int u, v;
        double w;
        bool operator<(Edge const& o) const { return w < o.w; }
    };

    struct Vertex {
        int u;
        Vertex* sib;
        Vertex* par;
        Vertex* child;
        Edge*  edge_to_par;
        Vertex(int _u)
          : u(_u), sib(nullptr), par(nullptr), child(nullptr), edge_to_par(nullptr) {}
    };

    inline double dist(int x1,int y1,int x2,int y2) {
        double dx = x1 - x2, dy = y1 - y2;
        return std::sqrt(dx*dx + dy*dy);
    }

    Vertex* prim_mst(int coords[][2], int dim) {
        const double INF = std::numeric_limits<double>::infinity();
        std::vector<double> min_cost(dim, INF);
        std::vector<int>    from(dim,    -1);
        std::vector<bool>   used(dim, false);
        std::vector<Vertex*> nodes(dim, nullptr);

        min_cost[0] = 0;
        for (int i = 0; i < dim; ++i)
            nodes[i] = new Vertex(i);

        for (int k = 0; k < dim; ++k) {
            int u = -1;
            for (int i = 0; i < dim; ++i)
                if (!used[i] && (u < 0 || min_cost[i] < min_cost[u]))
                    u = i;
            used[u] = true;
            
            if (k % 1000 == 0) {
                printf("Phase %d: currently added vertex %d\n", k, u);
            }

            if (from[u] != -1) {
                Vertex* child  = nodes[u];
                Vertex* parent = nodes[from[u]];
                child->par = parent;

                Edge* e = new Edge{ from[u], u,
                    dist(coords[from[u]][0], coords[from[u]][1],
                         coords[u][0],           coords[u][1]) };
                child->edge_to_par = e;

                if (!parent->child) parent->child = child;
                else {
                    Vertex* t = parent->child;
                    while (t->sib) t = t->sib;
                    t->sib = child;
                }
            }

            for (int v = 0; v < dim; ++v) {
                if (!used[v]) {
                    double d = dist(
                        coords[u][0], coords[u][1],
                        coords[v][0], coords[v][1]
                    );
                    if (d < min_cost[v]) {
                        min_cost[v] = d;
                        from[v]      = u;
                    }
                }
            }
        }

        return nodes[0];  // 0번 정점이 루트
    }

    // preorder 순회해서 list에 담는다
    int preorder(Vertex** list, Vertex* cur, int cnt=0) {
        list[cnt++] = cur;
        if (cur->child) cnt = preorder(list, cur->child, cnt);
        if (cur->sib)   cnt = preorder(list, cur->sib,   cnt);
        return cnt;
    }

    // tour cost + 파일 저장
    double mst_based_approx_2(int coords[][2], int dim) {
        // 1) MST 만들기
        Vertex* root = prim_mst(coords, dim);

        // 2) preorder 트리 순회
        Vertex** order = new Vertex*[dim];
        int cnt = preorder(order, root);

        // 3) tour cost 계산
        double cost = 0;
        for (int i = 0; i < cnt - 1; ++i) {
            int u = order[i]->u;
            int v = order[i+1]->u;
            cost += dist(
                coords[u][0], coords[u][1],
                coords[v][0], coords[v][1]
            );
        }
        cost += dist(
            coords[order[cnt-1]->u][0], coords[order[cnt-1]->u][1],
            coords[order[0]    ->u][0], coords[order[0]    ->u][1]
        );

        // std::ofstream fo("mst-based.txt");
        // for (int i = 0; i < cnt-1; ++i)
        //     fo << order[i]->u << ' ' << order[i+1]->u << '\n';
        // fo << order[cnt-1]->u << ' ' << order[0]->u << '\n';

        delete[] order;

        return cost;
    }

} // namespace approx_2