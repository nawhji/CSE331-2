#include <iostream>
#include <fstream>

namespace approx_2 {
        
    struct Edge {
        int u, v; // vertex
        double w; // cost
        bool operator<(const Edge& other) const {
            return w < other.w;
        }
    };

    struct Vertex {
        int u; // vertex num;
        Vertex* sib; // next next next
        Vertex* par;
        Vertex* child;
        Edge* edge_to_par;
    };

    const double MAX_COST = 1e9;

    Vertex* prim_mst(double **matrix, int dim);
    int preorder_tree_traversal(Vertex** list, Vertex* cur, int count);

    void save_edges_as_txt(const char* filename, Vertex** preorder, int count);

    double mst_based_approx_2(double **matrix, int dim) {
        Vertex* root_node = prim_mst(matrix, dim);
        Vertex** preorder = new Vertex*[dim];
        int count = preorder_tree_traversal(preorder, root_node, 0);

        double cost = 0;
        for (int i = 0; i < dim - 1; i++) {
            int u = preorder[i]->u;
            int v = preorder[i+1]->u;
            cost += matrix[u][v];
        }
        cost += matrix[preorder[0]->u][preorder[dim-1]->u];

        save_edges_as_txt("mst-based.txt", preorder, count);
        
        return cost;
    }

    int preorder_tree_traversal(Vertex** list, Vertex* cur, int count) {
        list[count++] = cur;
        if (cur->child)
            count = preorder_tree_traversal(list, cur->child, count);
        if (cur->sib)
            count = preorder_tree_traversal(list, cur->sib, count);
        return count;
    }

    Vertex* prim_mst(double **matrix, int dim) {
        int mst_vertex_num = 0;
        Vertex* loaded_vertex[dim];
        double min_cost[dim];
        int from[dim];
        bool visited[dim];

        for (int i = 0; i < dim; i++) {
            loaded_vertex[i] = nullptr;
            min_cost[i] = MAX_COST;
            from[i] = -1;
            visited[i] = false;
        }

        min_cost[0] = 0;
        int cur = 0;

        Vertex* vertex_obj[dim];
        for (int i = 0; i < dim; i++) vertex_obj[i] = nullptr;

        while (mst_vertex_num < dim) {
            double best = MAX_COST;
            for (int i = 0; i < dim; i++) {
                if (!visited[i] && min_cost[i] < best) {
                    best = min_cost[i];
                    cur = i;
                }
            }

            visited[cur] = true;
            Vertex* v = new Vertex;
            v->u = cur;
            v->sib = nullptr;
            v->par = nullptr;
            v->child = nullptr;
            v->edge_to_par = nullptr;
            vertex_obj[cur] = v;

            if (from[cur] != -1) {
                v->par = vertex_obj[from[cur]];

                Edge* e = new Edge;
                e->u = from[cur];
                e->v = cur;
                e->w = matrix[from[cur]][cur];
                v->edge_to_par = e;

                if (v->par->child == nullptr)
                    v->par->child = v;
                else {
                    Vertex* temp = v->par->child;
                    while (temp->sib != nullptr)
                        temp = temp->sib;
                    temp->sib = v;
                }
            }

            for (int j = 0; j < dim; j++) {
                if (!visited[j] && matrix[cur][j] < min_cost[j]) {
                    min_cost[j] = matrix[cur][j];
                    from[j] = cur;
                }
            }

            ++mst_vertex_num;
        }

        return vertex_obj[0]; // root는 항상 0번 정점
    }

    void save_edges_as_txt(const char* filename, Vertex** preorder, int count) {
        std::ofstream fout(filename);
        if (!fout.is_open()) {
            // std::cerr << "Failed to open file for writing: " << filename << std::endl;
            return;
        }

        for (int i = 0; i < count - 1; ++i) {
            fout << preorder[i]->u << " " << preorder[i + 1]->u << "\n";
        }
        fout << preorder[count - 1]->u << " " << preorder[0]->u << "\n"; // 마지막 → 처음으로
        fout.close();
    }
}