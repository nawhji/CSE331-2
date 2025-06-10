// starting point
// 출발점에서 모든 점까지(1번째 방문) -> 초기화

// 그 다음부터 최소 경로만 남기기

double double_max = 9007199254740992.0;
double min(double x, double y) { return (x <= y) ? x : y; }

double held_carp(double **matrix, int dim) {
    double** dp = new double*[1 << dim];
    for (int i = 0; i < (1 << dim); i++) {
        dp[i] = new double[dim];
        for (int j = 0; j < dim; j++) {
            dp[i][j] = double_max;
        }
    }

    // dp -> [비트 개수][도시 개수]
    // 도시 0부터 n-1 까지. 0은 출발점

    /* initialize */
    dp[1 << 0][0] = 0;
    for (int i = 1; i < dim; i++) {
        dp[(1 << 0) | (1 << i)][i] = matrix[0][i];
    }

    /* 도시 i개(2개부터 시작) 선택 -> n-1 P i * 도시 개수 */

    for (int S = 1; S < (1 << dim); S++) {
        for (int i = 1; i < dim; i++) {
            if (!(S & (1 << i))) continue;
            for (int j = 1; j < dim; j++) {
                if (i == j || !(S & (1 << j))) continue;

                dp[S][i] = min(dp[S][i], dp[S ^ (1 << i)][j] + matrix[j][i]);
            }
        }
    }   

    int end = (1 << dim) - 1;
    double result = double_max;
    // end = end ^ (1 << 0);

    for (int i = 1; i < dim; i++) {
        // dp[end][0] = min(dp[end ^ (1 << 0)][0], dp[end ^ (1 << 0)][i] + matrix[i][0]);
        result = min(result, dp[end][i] + matrix[i][0]);
    }

    return result;
}
