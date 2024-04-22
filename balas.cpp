#include <iostream>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <climits>
#include <array>

const int n_vars = 50;

std::array<std::array<double, 3>, n_vars> A;
const double threshold = 6.5;

// may need to unroll possibleremain and feasiblevec
// but maybe not because of spatial locality
void nx(int x, int ck, std::array<int, n_vars>& candidates,
        std::array<int, n_vars>& solution, std::array<int, n_vars>& values,
        int goalval, int& mingoalval, std::array<double, 3>& possibleremain,
        std::array<double, 3>& feasiblevec,
        const std::array<std::array<int, n_vars>, n_vars>& worse_keyframes) {
            
    if (ck == 1) {
        // Balas's second test
        if (goalval >= mingoalval) {
            return;
        }

        int violated = 0;
        for (int i = 0; i < 3; ++i) {
            double amount_off = threshold - feasiblevec[i];

            // Balas's first test
            if (amount_off > possibleremain[i]) {
                return;
            }

            else if (amount_off > 0) {
                violated = 1;

                // Modification G-1
                int possible = 0;
                double goal_ratio = amount_off / (mingoalval - goalval);
                for (int j = x + 1; j < n_vars; ++j) {
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
            // critical section
            if (goalval < mingoalval) {
                mingoalval = goalval;
                solution = values;
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
        values[x] = 1;
        for (int i = 0; i < 3; ++i) {
            feasiblevec[i] += A[x][i];
            possibleremain[i] -= A[x][i];
        }

        nx(x + 1, 1, candidates, solution, values, goalval+1, mingoalval,
           possibleremain, feasiblevec, worse_keyframes);

        // undo changes
        values[x] = 0;
        for (int i = 0; i < 3; ++i) {
            feasiblevec[i] -= A[x][i];
            possibleremain[i] += A[x][i];
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

    // Don't need to check anything if values are unchanged
    nx(x + 1, 0, new_candidates, solution, values, goalval, mingoalval,
       possibleremain, feasiblevec, worse_keyframes);
}

int main() {
    // Random values between 0 and 0.5
    std::srand(std::time(nullptr));
    for (int i = 0; i < n_vars; ++i) {
        for (int j = 0; j < 3; ++j) {
            A[i][j] = (std::rand() % 50) / 100.0;
        }
    }

    // fill possibleremain
    std::array<double, 3> possibleremain {};
    for (int i = 0; i < n_vars; ++i) {
        for (int j = 0; j < 3; ++j) {
            possibleremain[j] += A[i][j];
        }
    }

    // get the order to sort by
    std::array<int, 3> sortOrder = {0, 1, 2};
    std::sort(sortOrder.begin(), sortOrder.end(), [&](int a, int b) {
        return possibleremain[a] < possibleremain[b];
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

    std::array<double, 3> feasiblevec {};
    std::array<int, n_vars> values {};
    std::array<int, n_vars> solution {};
    int mingoalval = INT_MAX;

    // Assume all zeros is not a solution
    nx(0, 0, candidates, solution, values, 0, mingoalval,
       possibleremain, feasiblevec, worse_keyframes);

    std::cout << "mingoalval = " << mingoalval << std::endl;
    std::cout << "solution = ";
    for (int i = 0; i < n_vars; ++i) {
        std::cout << solution[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}
