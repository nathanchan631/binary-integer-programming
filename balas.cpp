#include <iostream>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <climits>
#include <array>
#include <chrono>
#include <vector>

const int n_vars = 50;

std::array<std::array<double, 3>, n_vars> A;
//  = {{
//     {.3, 0., .3},
//     {.2, 0., .3},
//     {.3, .1, .2},
//     {.1, .3, .15},
//     {.2, .4, .1},
//     {0., .5, 0.},
// }};

const double threshold = 6.5;

// may need to unroll remaining_scores and scores
// but maybe not because of spatial locality
void nx(int x, int ck, int objective, std::array<int, n_vars>& path,
        std::array<double, 3>& scores, std::array<double, 3>& remaining_scores, 
        int& min_objective, double& solution_sum, std::array<int, n_vars>& solution,
        std::array<int, n_vars>& candidates, 
        const std::array<std::array<int, n_vars>, n_vars>& worse_keyframes) {
            
    if (ck == 1) {
        // Balas's second test
        if (objective > min_objective) {
            return;
        }

        else if (objective == min_objective) {
            if (    scores[0] > threshold &&
                    scores[1] > threshold &&
                    scores[2] > threshold &&
                    scores[0] + scores[1] + scores[2] > solution_sum ) {
                solution_sum = scores[0] + scores[1] + scores[2];
                min_objective = objective;
                solution = path;
            }
            return;
        }

        int violated = 0;
        for (int i = 0; i < 3; ++i) {
            const double amount_off = threshold - scores[i];

            // Balas's first test
            if (amount_off > remaining_scores[i]) {
                return;
            }

            else if (amount_off > 0) {
                violated = 1;

                // Glover's first modification
                int possible = 0;
                const double goal_ratio = amount_off / (min_objective - objective);
                for (int j = x; j < n_vars; ++j) {
                    if (A[j][i] >= goal_ratio) {
                        possible = 1;
                        break;
                    }
                }

                if (!possible) {
                    return;
                }
            }
        }

        // If the current solution is feasible, update the best result
        if (violated == 0) {
            if (objective < min_objective) {
                min_objective = objective;
                solution_sum = scores[0] + scores[1] + scores[2];
                solution = path;
            }
            return;
        }
    }

    // Out of bounds
    if (x >= n_vars) {
        return;
    }

    // One branch
    if (candidates[x] == 1) {
        path[x] = 1;
        for (int i = 0; i < 3; ++i) {
            scores[i] += A[x][i];
            remaining_scores[i] -= A[x][i];
        }

        nx(x + 1, 1, objective+1, path, scores, remaining_scores,
           min_objective, solution_sum, solution, candidates, worse_keyframes);

        // undo changes
        path[x] = 0;
        for (int i = 0; i < 3; ++i) {
            scores[i] -= A[x][i];
            remaining_scores[i] += A[x][i];
        }
    }

    // Zero branch

    // Custom modification
    // If a variable is worse than the one we're skipping here
    // then there is no point in exploring it
    std::array<int, n_vars> new_candidates = candidates;
    for (int i = x + 1; i < n_vars; ++i) {
        if (worse_keyframes[x][i] == 1) {
            new_candidates[i] = 0;
        }
    }

    // Don't need to check anything if path are unchanged
    nx(x + 1, 0, objective, path, scores, remaining_scores,
       min_objective, solution_sum, solution, new_candidates, worse_keyframes);
}


int main() {
    // Random path between 0 and 0.5
    std::srand(std::time(nullptr));
    for (int i = 0; i < n_vars; ++i) {
        for (int j = 0; j < 3; ++j) {
            A[i][j] = (std::rand() % 50) / 100.0;
        }
    }

    // fill remaining_scores
    std::array<double, 3> remaining_scores {};
    for (int i = 0; i < n_vars; ++i) {
        for (int j = 0; j < 3; ++j) {
            remaining_scores[j] += A[i][j];
        }
    }

    // get the order to sort by
    std::array<int, 3> sortOrder = {0, 1, 2};
    std::sort(sortOrder.begin(), sortOrder.end(), [&](int a, int b) {
        return remaining_scores[a] < remaining_scores[b];
    });

    // sort A based on which constraints are lowest
    std::sort(A.begin(), A.end(),
        [&](const std::array<double, 3>& a, const std::array<double, 3>& b) {
            for (int i = 0; i < 3; ++i) {
                int col = sortOrder[i];
                if (a[col] != b[col]) {
                    return a[col] > b[col];
                }
            }
            return false;
        });

    // fill worse keyframes array
    std::array<std::array<int, n_vars>, n_vars> worse_keyframes {};
    for (int i = 0; i < n_vars; ++i) {
        for (int j = i + 1; j < n_vars; ++j) {
            std::array<double, 3>& frame1 = A[i];
            std::array<double, 3>& frame2 = A[j];
            bool is_worse = true;

            for (int k = 0; k < 3; ++k) {
                if (frame1[k] < frame2[k]) {
                    is_worse = false;
                    break;
                }
            }
            if (is_worse) {
                worse_keyframes[i][j] = 1;
            }
        }
    }

    std::array<int, n_vars> candidates;
    candidates.fill(1);

    std::array<double, 3> scores {};
    std::array<int, n_vars> path {};
    std::array<int, n_vars> solution {};
    int min_objective = INT_MAX;
    double solution_sum = 0;


    auto start = std::chrono::steady_clock::now();
    
    // Assume all zeros is not a solution
    nx(0, 0, 0, path, scores, remaining_scores, min_objective,
       solution_sum, solution, candidates, worse_keyframes);
    
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    std::cout << "time: " << std::chrono::duration<double, std::milli>(diff).count() << std::endl;

    
    // min_objective = INT_MAX;
    // auto start2 = std::chrono::steady_clock::now();
    
    // // Assume all zeros is not a solution
    // nx_2(solution, path, min_objective,
    //    remaining_scores, scores, worse_keyframes);
    
    // auto end2 = std::chrono::steady_clock::now();
    // auto diff2 = end2 - start2;
    // std::cout << "time: " << std::chrono::duration<double, std::milli>(diff2).count() << std::endl;

    std::cout << "min_objective = " << min_objective << std::endl;
    std::cout << "solution sum = " << solution_sum << std::endl;
    std::cout << "solution = ";
    for (int i = 0; i < n_vars; ++i) {
        std::cout << solution[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}
