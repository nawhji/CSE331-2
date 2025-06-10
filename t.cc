#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <chrono>
#include <utility>
#include "held-carp.cc"
#include "kruskal_mst.cc"
#include "kruskal_ondemand.cc"
#include "mst-based-2-approx.cc"
#include "approx_ondemand.cc"

using namespace std;
using namespace std::chrono;

int mona_lisa_my_algo();
int mona_lisa_2_approx();

int main() {
    string file_name;
    cout << "test file name: ";
    cin >> file_name;

    string algo_name;
    cout << "test algorithm: (ex. held-carp, kruskal_mst, 2-approx) ";
    cin >> algo_name;

    if (file_name == "mona-lisa" && algo_name == "kruskal_mst") {
        mona_lisa_my_algo();
        return 1;
    } else if (file_name == "mona-lisa" && algo_name == "2-approx") {
        mona_lisa_2_approx();
    }

    file_name = "testcases/" + file_name + ".tsp";
    ifstream file(file_name);
    

    int dim;
    int** array_2d;

    string line;
    bool found = false;
    while (getline(file, line)) {
        if (found) {
            if (line.find("EOF") != string::npos) break;
            stringstream ss(line);
            string num;
            int *nums = new int[3];
            int i = 0;

            while (ss >> num) {
                nums[i] = stoi(num);
                i++;
            }
            array_2d[nums[0] - 1][0] = nums[1];
            array_2d[nums[0] - 1][1] = nums[2];
        }

        if (line.find("DIMENSION") != string::npos) {
            size_t pos = line.find("DIMENSION: ");
            string num_str = line.substr(pos + 12);
            cout << "num_str: " << num_str << endl;
            dim = stoi(num_str);
            cout << "after stoi: " << dim << endl;

            array_2d = new int*[dim];
            for (int i = 0; i < dim; i++)
                array_2d[i] = new int[2];
        }

        if (line.find("NODE_COORD_SECTION") != string::npos) {
            found = true;
        }
    } file.close();

    double **matrix;
    matrix = new double*[dim];
    for (int i = 0; i < dim; i++) {
        matrix[i] = new double[dim];
    }

    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            if (i == j) {
                matrix[i][j] = 0;
                continue;
            }
            matrix[i][j] = sqrt(pow(array_2d[i][0] - array_2d[j][0], 2) + pow(array_2d[i][1] - array_2d[j][1], 2));
        }
    }

    if (algo_name == "held-carp") {
        auto start = high_resolution_clock::now();
        cout << "result of held-carp: " << held_carp(matrix, dim) << endl;
        auto end = high_resolution_clock::now();
        cout << "Elapsed time: " << duration<double>(end - start).count() << endl;
    }
    else if (algo_name == "kruskal_mst") {
        auto start = high_resolution_clock::now();
        cout << "result of kruskal_mst: " << kruskal_mst(matrix, dim) << endl;
        auto end = high_resolution_clock::now();
        cout << "Elapsed time: " << duration<double>(end - start).count() << endl;
    }
    else {//if (algo_name == "2-approx") {
        auto start = high_resolution_clock::now();
        cout << "result of mst-approx-2: " << approx_2::mst_based_approx_2(matrix, dim) << endl;
        auto end = high_resolution_clock::now();
        cout << "Elapsed time: " << duration<double>(end - start).count() << endl;
    }

    for (int i = 0; i < dim; i++) {
        delete[] array_2d[i];
        delete[] matrix[i];
    }
    delete array_2d;
    delete matrix;

    return 0;
}

int mona_lisa_my_algo() {
    ifstream file("testcases/mona-lisa100K.tsp");
    vector<pair<int, int>> coords;
    string line;
    bool found = false;
    while (getline(file, line)) {
        if (found) {
            if (line.find("EOF") != string::npos) break;
            stringstream ss(line);
            int idx, x, y;
            ss >> idx >> x >> y;
            coords.push_back({x, y});
        }
        if (line.find("NODE_COORD_SECTION") != string::npos) found = true;
    }
    file.close();
    
    auto start = high_resolution_clock::now();
    double result = kruskal_ondemand::kruskal_mst_ondemand(coords);
    auto end = high_resolution_clock::now();
    cout << "result: " << result << endl;
    cout << "Elapsed time: " << duration<double>(end - start).count() << endl;
    return 0;
}

int mona_lisa_2_approx() {
    static int coords[100000][2];
    int dim = 0;

    ifstream file("testcases/mona-lisa100K.tsp");
    string line;
    bool found = false;

    while (getline(file, line)) {
        if (found) {
            if (line.find("EOF") != string::npos) break;
            stringstream ss(line);
            int idx, x, y;
            ss >> idx >> x >> y;
            coords[dim][0] = x;
            coords[dim][1] = y;
            ++dim;
        }
        if (!found && line.find("NODE_COORD_SECTION") != string::npos) {
            found = true;
        }
    }
    file.close();

    clock_t start = clock();
    double result = approx_2_ondemand::mst_based_approx_2(coords, dim);
    clock_t end = clock();

    printf("result: %.6f\n", result);
    printf("Elapsed time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    return 0;

}