/*
implementing balas's algorithm: https://ieeexplore.ieee.org/document/9803720
collection of improvements from Peterson: https://www.jstor.org/stable/2628090?seq=8
balas's paper: https://www.jstor.org/stable/167850?seq=12
glover + zionts improvements: https://homes.di.unimi.it/righini/Didattica/ComplementiRicercaOperativa/MaterialeCRO/Glover%20Zionts%201965%20-%20Note%20on%20Balas%20algorithm.pdf

Also tried (but slower):
- Surrogate constraint (Glover)
- Testing if an entry is greater than (leeway=remaining-amount_needed) (Fleischmann 2nd Modification)
- Gomory Cut + reduction (Lemke-Spielberg) (this one doesn't fit well within our framework)
- Peterson's Modification R-1
- testing if max (k=min_objective-objective) elements add up to amount_needed

Many tests do not apply since all entries of A are nonnegative, ie
- Fleischmann's 1st modification
- Brauer's modification
- Tests GZ and P1 from Hrouda
- and parts of other tests, such as Balas's P-Test No. 3 and Peterson's Modification R-1

Also looked at Balas's constraint on the objective function from the filter problem,
but unlikely to help since the algorithm converges quickly, and spends most of its time
exploring dead ends.

Tested with 50 random keyframes, expecting to choose 20 keyframes.
*/

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <climits>
#include <array>
#include <vector>
#include <chrono>
#include <random>


void solve(const std::vector<std::pair<int, std::vector<double>>>& A,
        const int n_vars, const double threshold,
        int x, int objective, std::vector<int>& path,
        std::array<double, 3>& scores, std::array<double, 3>& remaining_scores, 
        int& min_objective, double& solution_sum, std::vector<int>& solution,
        std::vector<int>& excluded, 
        const std::vector<std::vector<int>>& worse_keyframes) {

    // out of bounds
    if (x >= n_vars) {
        return;
    }

    // cuts
    for (int i = 0; i < 3; ++i) {
        const double amount_needed = threshold - scores[i];

        // Balas step 3
        if (amount_needed > remaining_scores[i]) {
            return;
        }

        // Glover's 2nd modification
        else if (amount_needed > 0.0) {
            int possible = 0;
            const int objective_diff = min_objective - objective;
            for (int j = x; j < n_vars; ++j) {
                if (A[j].second[i] * objective_diff >= amount_needed && excluded[j] == 0) {
                    possible = 1;
                    break;
                }
            }

            if (!possible) {
                return;
            }
        }
    }

    // branching
    if (excluded[x] == 0) {

        // "1" branch

        // update values
        objective += 1;
        path[x] = 1;
        for (int i = 0; i < 3; ++i) {
            scores[i] += A[x].second[i];
            remaining_scores[i] -= A[x].second[i];
        }
        
        // check feasibility and update the best result
        if (scores[0] > threshold && scores[1] > threshold && scores[2] > threshold) {
            if (objective < min_objective || (objective == min_objective &&
                                              scores[0] + scores[1] + scores[2] > solution_sum)) {
                min_objective = objective;
                solution_sum = scores[0] + scores[1] + scores[2];
                solution = path;
            }
        }

        // if not a solution, but an improved solution is possible,
        // then explore the branch further
        else if (objective < min_objective) {
            solve(A, n_vars, threshold, x+1, objective, path, scores, remaining_scores,
                  min_objective, solution_sum, solution, excluded, worse_keyframes);
        }

        // revert changes
        objective -= 1;
        path[x] = 0;
        for (int i = 0; i < 3; ++i) {
            scores[i] -= A[x].second[i];
        }


        // "0" branch

        // Custom modification: if a variable is worse than the one we're skipping,
        // then there is no point in exploring it

        // Note: we only do this when the current keyframe has not already been excluded,
        // because if the current keyframe is excluded, then all keyframes worse than it
        // would have already been excluded too by the same keyframe that excluded it
        
        std::vector<int> removed;
        for (int i = x + 1; i < n_vars; ++i) {
            if (worse_keyframes[x][i] == 1 && excluded[i] == 0) {

                excluded[i] = 1;
                for (int j = 0; j < 3; ++j) {
                    remaining_scores[j] -= A[i].second[j];
                }
                
                removed.push_back(i);
            }
        }

        solve(A, n_vars, threshold, x+1, objective, path, scores, remaining_scores,
            min_objective, solution_sum, solution, excluded, worse_keyframes);

        // revert changes
        for (std::vector<int>::size_type i = 0; i < removed.size(); ++i) {
            int index = removed[i];
            
            excluded[index] = 0;
            for (int j = 0; j < 3; ++j) {
                remaining_scores[j] += A[index].second[j];
            }
        }

        for (int i = 0; i < 3; ++i) {
            remaining_scores[i] += A[x].second[i];
        }
    }

    else {
        // "0" branch
        solve(A, n_vars, threshold, x+1, objective, path, scores, remaining_scores,
               min_objective, solution_sum, solution, excluded, worse_keyframes);
    }
}


void start(const std::vector<std::vector<double>>& A,
           const double threshold, const double min_score_sum,
           std::vector<int>& selected_indices) {

    // package A into pairs of index, keyframe scores
    // to keep track of indices after filtering and sorting
    std::vector<std::pair<int, std::vector<double>>> A_indexed;
    for (std::vector<int>::size_type i = 0; i < A.size(); ++i)
        A_indexed.push_back({i, A[i]});
    
    // Iterate over A_indexed and remove elements that do not meet the sum condition
    A_indexed.erase(std::remove_if(A_indexed.begin(), A_indexed.end(),
        [min_score_sum](const std::pair<int, std::vector<double>>& keyframe) {
            double rowSum = 0;
            for (double value : keyframe.second)
                rowSum += value;

            return rowSum < min_score_sum;
        }
    ), A_indexed.end());

    const int n_vars = A_indexed.size();

    // fill remaining_scores
    std::array<double, 3> remaining_scores {};
    for (int i = 0; i < n_vars; ++i)
        for (int j = 0; j < 3; ++j)
            remaining_scores[j] += A_indexed[i].second[j];

    // get the order to sort by
    std::array<int, 3> sortOrder = {0, 1, 2};
    std::sort(sortOrder.begin(), sortOrder.end(), [&](int a, int b) {
        return remaining_scores[a] < remaining_scores[b];
    });

    // sort A based on which constraints are lowest (Peterson's modification R-2)
    std::sort(A_indexed.begin(), A_indexed.end(),
        [&](const std::pair<int, std::vector<double>>& a,
            const std::pair<int, std::vector<double>>& b) {
            for (int i = 0; i < 3; ++i) {
                int col = sortOrder[i];
                if (a.second[col] != b.second[col]) {
                    return a.second[col] > b.second[col];
                }
            }
            return false;
        }
    );

    // fill worse keyframes array
    // Note: could make multiple of these for each combination of unsatisftied
    // constraints, but it's probably not needed for our case since keyframes with
    // a higher value in one dimension will likely have higher values in all dimensions
    std::vector<std::vector<int>> worse_keyframes(n_vars, std::vector<int>(n_vars));
    for (int i = 0; i < n_vars; ++i) {
        for (int j = i + 1; j < n_vars; ++j) {
            const std::vector<double>& frame1 = A_indexed[i].second;
            const std::vector<double>& frame2 = A_indexed[j].second;
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

    int min_objective = INT_MAX;
    double solution_sum = 0;
    std::array<double, 3> scores {};
    std::vector<int> path(n_vars);
    std::vector<int> solution(n_vars);
    std::vector<int> excluded(n_vars);

    auto start = std::chrono::steady_clock::now();
    
    // Assumes all zeros is not a solution
    solve(A_indexed, n_vars, threshold, 0, 0, path, scores, remaining_scores,
          min_objective, solution_sum, solution, excluded, worse_keyframes);
    
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;

    std::cout << "time: " << std::chrono::duration<double, std::milli>(diff).count() << std::endl;
    std::cout << "min_objective = " << min_objective << std::endl;
    std::cout << "solution sum = " << solution_sum << std::endl;

    // min_objective = INT_MAX;
    // scores = {0.0, 0.0, 0.0};
    // std::fill(path.begin(), path.end(), 0);
    // std::fill(excluded.begin(), excluded.end(), 0);

    // auto start2 = std::chrono::steady_clock::now();
    
    // // Assumes all zeros is not a solution
    // solve2(A_indexed, n_vars, threshold, 0, 0, path, scores, remaining_scores,
    //       min_objective, solution_sum, solution, excluded, worse_keyframes);
    
    // auto end2 = std::chrono::steady_clock::now();
    // auto diff2 = end2 - start2;

    // std::cout << "time: " << std::chrono::duration<double, std::milli>(diff2).count() << std::endl;
    // std::cout << "min_objective = " << min_objective << std::endl;
    // std::cout << "solution sum = " << solution_sum << std::endl;

    selected_indices.clear();
    for (std::vector<int>::size_type i = 0; i < solution.size(); ++i)
        if (solution[i] == 1)
            selected_indices.push_back(A_indexed[i].first);
}


int main() {
    const int n_vars = 50;
    std::vector<std::vector<double>> A;
    //  = {{
        // {.2, 0., .3},
        // {.1, .4, 0.},
        // {.3, .1, .2},
        // {.1, .3, .1},
        // {0., .5, 0.1},
        // {.3, 0., .3}
    //     {0.18,0.24,0.55},
    //     {0.2,0.25,0.54},
    //     {0.18,0.25,0.49},
    //     {0.16,0.22,0.42},
    //     {0.14,0.19,0.31},
    //     {0.13,0.15,0.21},
    //     {0.1,0.14,0.13},
    //     {0.1,0.13,0.07},
    //     {0.09,0.12,0.04},
    //     {0.08,0.11,0.02},
    //     {0.08,0.1,0.02},
    //     {0.07,0.07,0.01},
    //     {0.05,0.06,0.01},
    //     {0.03,0.03,0.01},
    //     {0.01,0.01,0.01},
    //     {0,0,0},
    //     {0.01,0.01,0},
    //     {0,0,0},
    //     {0.01,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0.01,0.01,0},
    //     {0.01,0.01,0},
    //     {0,0,0},
    //     {0.01,0.01,0},
    //     {0.01,0.01,0},
    //     {0.01,0.01,0},
    //     {0.01,0.01,0},
    //     {0.02,0.03,0.01},
    //     {0.03,0.03,0.01},
    //     {0.02,0.03,0},
    //     {0.02,0.03,0},
    //     {0.02,0.02,0},
    //     {0.04,0.05,0.01},
    //     {0.03,0.05,0.01},
    //     {0.02,0.03,0.01},
    //     {0.03,0.03,0.01},
    //     {0.03,0.04,0.01},
    //     {0.05,0.08,0.02},
    //     {0.09,0.12,0.06},
    //     {0.09,0.13,0.07},
    //     {0.1,0.14,0.09},
    //     {0.11,0.15,0.12},
    //     {0.12,0.16,0.14},
    //     {0.12,0.16,0.17},
    //     {0.12,0.16,0.19},
    //     {0.12,0.17,0.23},
    //     {0.13,0.17,0.3},
    //     {0.13,0.17,0.34},
    //     {0.14,0.17,0.35},
    //     {0.14,0.19,0.35},
    //     {0.16,0.21,0.39},
    //     {0.17,0.21,0.41},
    //     {0.16,0.21,0.47},
    //     {0.16,0.21,0.48},
    //     {0.18,0.23,0.46},
    //     {0.19,0.25,0.51},
    //     {0.18,0.24,0.49},
    //     {0.17,0.21,0.49},
    // }};

    const double threshold = 6.5;
    const double min_score_sum = 0.025;
    std::vector<int> selected_indices;

    // Initialize random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 0.5);

    // Generate random values and populate the matrix
    for (int i = 0; i < n_vars; ++i) {
        std::vector<double> row;
        for (int j = 0; j < 3; ++j) {
            row.push_back(dis(gen));
        }
        A.push_back(row);
    }

    start(A, threshold, min_score_sum, selected_indices);

    std::cout << "solution = ";
    for (std::vector<int>::size_type i = 0; i < selected_indices.size(); ++i)
        std::cout << selected_indices[i] << " ";
    std::cout << std::endl;

    return 0;
}