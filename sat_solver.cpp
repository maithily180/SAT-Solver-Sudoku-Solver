#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <array>
#include <numeric>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <future>
#include <memory>
#include <fstream>
#include <sstream>     // ADD THIS for std::ostringstream
#include <sys/stat.h>  // ADD THIS for mkdir
#include <sys/types.h> // ADD THIS for mkdir


// --- Timer Utility ---
class Timer {
public:
    Timer() : start_time(std::chrono::high_resolution_clock::now()) {}

    void reset() {
        start_time = std::chrono::high_resolution_clock::now();
    }

    double elapsed_s() const {
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end_time - start_time;
        return diff.count();
    }

private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
};

// --- Dynamic Grid Structure ---
struct Grid {
    int size;
    int box_size;
    std::vector<std::vector<int>> cells;
    
    Grid(int n = 9) : size(n) {
        box_size = static_cast<int>(std::sqrt(n));
        if (box_size * box_size != n) {
            std::cerr << "Error: Grid size must be a perfect square (4, 9, 16, 25, etc.)" << std::endl;
            exit(1);
        }
        cells.resize(n, std::vector<int>(n, 0));
    }
    
    std::vector<int>& operator[](int r) {
        return cells[r];
    }
    
    const std::vector<int>& operator[](int r) const {
        return cells[r];
    }
};


void print_grid(const Grid& grid) {
    int n = grid.size;
    int box = grid.box_size;
    
    for (int i = 0; i < n + box + 1; ++i) std::cout << "-";
    std::cout << std::endl;
    
    for (int r = 0; r < n; ++r) {
        std::cout << "|";
        for (int c = 0; c < n; ++c) {
            if (grid[r][c] == 0) {
                std::cout << " .";
            } else {
                std::cout << std::setw(2) << grid[r][c];
            }
            
            if ((c + 1) % box == 0) {
                std::cout << " |";
            }
        }
        std::cout << std::endl;
        
        if ((r + 1) % box == 0 && r + 1 < n) {
            for (int i = 0; i < n + box + 1; ++i) std::cout << "-";
            std::cout << std::endl;
        }
    }
    
    for (int i = 0; i < n + box + 1; ++i) std::cout << "-";
    std::cout << std::endl;
}


Grid load_puzzle_from_string(const std::string& s, int grid_size = 9) {
    Grid grid(grid_size);
    int idx = 0;
    
    for (int r = 0; r < grid_size; ++r) {
        for (int c = 0; c < grid_size; ++c) {
            if (idx >= (int)s.length()) {
                std::cerr << "Error: Puzzle string too short" << std::endl;
                exit(1);
            }
            
            char ch = s[idx++];
            if (ch >= '0' && ch <= '9') {
                grid[r][c] = ch - '0';
            } else if (ch >= 'A' && ch <= 'G') {
                grid[r][c] = ch - 'A' + 10;
            }
        }
    }
    
    return grid;
}

Grid load_puzzle_from_file(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        exit(1);
    }
    
    int grid_size;
    if (!(file >> grid_size)) {
        std::cerr << "Error: Cannot read grid size from " << filename << std::endl;
        exit(1);
    }
    
    int box_size = static_cast<int>(std::sqrt(grid_size));
    if (box_size * box_size != grid_size) {
        std::cerr << "Error: Grid size " << grid_size 
                  << " is not a perfect square (must be 4, 9, 16, 25, etc.)" << std::endl;
        exit(1);
    }
    
    Grid grid(grid_size);
    
    for (int r = 0; r < grid_size; ++r) {
        for (int c = 0; c < grid_size; ++c) {
            if (!(file >> grid[r][c])) {
                std::cerr << "Error: Cannot read value at position (" 
                          << r << "," << c << ") in " << filename << std::endl;
                exit(1);
            }
            
            if (grid[r][c] < 0 || grid[r][c] > grid_size) {
                std::cerr << "Error: Invalid value " << grid[r][c] 
                          << " at position (" << r << "," << c 
                          << ") - must be 0-" << grid_size << std::endl;
                exit(1);
            }
        }
    }
    
    file.close();
    
    std::cout << "Successfully loaded " << grid_size << "x" << grid_size 
              << " puzzle from " << filename << std::endl;
    
    return grid;
}


bool verify_sudoku_solution(const Grid& grid) {
    int n = grid.size;
    int box = grid.box_size;
    
    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < n; ++c) {
            if (grid[r][c] < 1 || grid[r][c] > n) {
                std::cout << "  ✗ Invalid value at (" << r << "," << c 
                          << "): " << grid[r][c] << std::endl;
                return false;
            }
        }
    }
    
    for (int r = 0; r < n; ++r) {
        std::set<int> seen;
        for (int c = 0; c < n; ++c) {
            if (seen.count(grid[r][c])) {
                std::cout << "  ✗ Duplicate " << grid[r][c] 
                          << " in row " << r << std::endl;
                return false;
            }
            seen.insert(grid[r][c]);
        }
    }
    
    for (int c = 0; c < n; ++c) {
        std::set<int> seen;
        for (int r = 0; r < n; ++r) {
            if (seen.count(grid[r][c])) {
                std::cout << "  ✗ Duplicate " << grid[r][c] 
                          << " in column " << c << std::endl;
                return false;
            }
            seen.insert(grid[r][c]);
        }
    }
    
    for (int box_r = 0; box_r < box; ++box_r) {
        for (int box_c = 0; box_c < box; ++box_c) {
            std::set<int> seen;
            for (int r = box_r * box; r < (box_r + 1) * box; ++r) {
                for (int c = box_c * box; c < (box_c + 1) * box; ++c) {
                    if (seen.count(grid[r][c])) {
                        std::cout << "  ✗ Duplicate " << grid[r][c] 
                                  << " in box (" << box_r << "," << box_c << ")" << std::endl;
                        return false;
                    }
                    seen.insert(grid[r][c]);
                }
            }
        }
    }
    
    return true;
}

// --- CSV Export Functions ---
void init_csv_files() {
    // Create results directory using C++17 filesystem
    #include <sys/stat.h>
    #include <sys/types.h>
    
    // Create directory (mkdir returns 0 on success, -1 if already exists)
    #ifdef _WIN32
        mkdir("results");
    #else
        mkdir("results", 0755);
    #endif
    
    // Initialize backtracking CSV
    std::ofstream bt_csv("results/backtracking_results.csv");
    bt_csv << "Puzzle,Clues,Configuration,Time_s,Decisions,Backtracks\n";
    bt_csv.close();
    
    // Initialize SAT CSV
    std::ofstream sat_csv("results/sat_results.csv");
    sat_csv << "Puzzle,Clues,Configuration,Time_s,Decisions,UnitProps,Backtracks,Backjumps,Timeout\n";    // ADD Backjumps
    sat_csv.close();
}


void export_backtracking_result(const std::string& puzzle_name, int clues,
                                 const std::string& config, double time,
                                 long long decisions, long long backtracks) {
    std::ofstream csv("results/backtracking_results.csv", std::ios::app);
    csv << puzzle_name << ","
        << clues << ","
        << config << ","
        << std::fixed << std::setprecision(6) << time << ","
        << decisions << ","
        << backtracks << "\n";
    csv.close();
}

void export_sat_result(const std::string& puzzle_name, int clues,
                       const std::string& config, double time,
                       long long decisions, long long unit_props,
                       long long backtracks, long long backjumps, bool timeout) {    // ADD backjumps
    std::ofstream csv("results/sat_results.csv", std::ios::app);
    csv << puzzle_name << ","
        << clues << ","
        << config << ",";
    
    if (timeout) {
    csv << "TIMEOUT,0,0,0,0,YES\n";    // ADD 0 for backjumps
    } else {
    csv << std::fixed << std::setprecision(6) << time << ","
        << decisions << ","
        << unit_props << ","
        << backtracks << ","
        << backjumps << ",NO\n";    // ADD backjumps
    }

    csv.close();
}

int count_clues(const Grid& grid) {
    int count = 0;
    for (int r = 0; r < grid.size; ++r) {
        for (int c = 0; c < grid.size; ++c) {
            if (grid[r][c] != 0) count++;
        }
    }
    return count;
}


// [BACKTRACKING SOLVER CLASS - keeping your existing one]
class BacktrackingSolver {
public:
    Grid grid;
    Grid solution;
    int size;
    int box_size;
    
    struct {
        long long decisions = 0;
        long long backtracks = 0;
    } stats;
    
    BacktrackingSolver(const Grid& puzzle) 
        : grid(puzzle), solution(puzzle), size(puzzle.size), box_size(puzzle.box_size) {}
    
    bool solve(bool use_mrv, bool use_fc) {
        stats = {0, 0};
        
        std::map<std::pair<int,int>, std::set<int>> domains;
        initialize_domains(domains);
        
        return solve_recursive(use_mrv, use_fc, domains);
    }
    
private:
    void initialize_domains(std::map<std::pair<int,int>, std::set<int>>& domains) {
        for (int r = 0; r < size; ++r) {
            for (int c = 0; c < size; ++c) {
                if (grid[r][c] == 0) {
                    for (int val = 1; val <= size; ++val) {
                        domains[{r, c}].insert(val);
                    }
                }
            }
        }
        
        for (int r = 0; r < size; ++r) {
            for (int c = 0; c < size; ++c) {
                if (grid[r][c] != 0) {
                    int val = grid[r][c];
                    std::vector<std::pair<int,int>> dummy_changes;
                    propagate_constraints(r, c, val, domains, dummy_changes);
                }
            }
        }
    }
    
    bool find_empty_cell(int& r, int& c) {
        for (r = 0; r < size; ++r) {
            for (c = 0; c < size; ++c) {
                if (grid[r][c] == 0) return true;
            }
        }
        return false;
    }
    
    bool find_mrv_cell(int& r_best, int& c_best, 
                       const std::map<std::pair<int,int>, std::set<int>>& domains, 
                       bool use_fc) {
        int min_remaining = size + 1;
        bool found = false;
        
        for (int r = 0; r < size; ++r) {
            for (int c = 0; c < size; ++c) {
                if (grid[r][c] == 0) {
                    int remaining = 0;
                    
                    if (use_fc) {
                        auto it = domains.find({r, c});
                        if (it != domains.end()) {
                            remaining = it->second.size();
                        }
                    } else {
                        for (int val = 1; val <= size; ++val) {
                            if (is_valid(r, c, val)) {
                                remaining++;
                            }
                        }
                    }
                    
                    if (remaining == 0) {
                        r_best = r;
                        c_best = c;
                        return true;
                    }
                    
                    if (remaining < min_remaining) {
                        min_remaining = remaining;
                        r_best = r;
                        c_best = c;
                        found = true;
                    }
                }
            }
        }
        return found;
    }
    
    bool propagate_constraints(int r, int c, int val,
                               std::map<std::pair<int,int>, std::set<int>>& domains,
                               std::vector<std::pair<int,int>>& changes) {
        for (int col = 0; col < size; ++col) {
            if (col != c && domains[{r, col}].count(val)) {
                domains[{r, col}].erase(val);
                changes.push_back({r, col});
                if (domains[{r, col}].empty() && grid[r][col] == 0) {
                    return false;
                }
            }
        }
        
        for (int row = 0; row < size; ++row) {
            if (row != r && domains[{row, c}].count(val)) {
                domains[{row, c}].erase(val);
                changes.push_back({row, c});
                if (domains[{row, c}].empty() && grid[row][c] == 0) {
                    return false;
                }
            }
        }
        
        int box_start_r = (r / box_size) * box_size;
        int box_start_c = (c / box_size) * box_size;
        for (int i = 0; i < box_size; ++i) {
            for (int j = 0; j < box_size; ++j) {
                int row = box_start_r + i;
                int col = box_start_c + j;
                if ((row != r || col != c) && domains[{row, col}].count(val)) {
                    domains[{row, col}].erase(val);
                    changes.push_back({row, col});
                    if (domains[{row, col}].empty() && grid[row][col] == 0) {
                        return false;
                    }
                }
            }
        }
        
        return true;
    }
    
    bool is_valid(int r, int c, int val) {
        for (int col = 0; col < size; ++col) {
            if (grid[r][col] == val) return false;
        }
        
        for (int row = 0; row < size; ++row) {
            if (grid[row][c] == val) return false;
        }
        
        int box_start_r = (r / box_size) * box_size;
        int box_start_c = (c / box_size) * box_size;
        for (int i = 0; i < box_size; ++i) {
            for (int j = 0; j < box_size; ++j) {
                if (grid[box_start_r + i][box_start_c + j] == val) {
                    return false;
                }
            }
        }
        
        return true;
    }
    
    bool solve_recursive(bool use_mrv, bool use_fc,
                        std::map<std::pair<int,int>, std::set<int>> domains) {
        int r, c;
        
        if (use_mrv) {
            if (!find_mrv_cell(r, c, domains, use_fc)) {
                solution = grid;
                return true;
            }
        } else {
            if (!find_empty_cell(r, c)) {
                solution = grid;
                return true;
            }
        }
        
        std::vector<int> possible_values;
        if (use_fc) {
            auto it = domains.find({r, c});
            if (it != domains.end()) {
                for (int val : it->second) {
                    possible_values.push_back(val);
                }
            }
        } else {
            for (int val = 1; val <= size; ++val) {
                if (is_valid(r, c, val)) {
                    possible_values.push_back(val);
                }
            }
        }
        
        if (possible_values.empty() && grid[r][c] == 0) {
            stats.backtracks++;
            return false;
        }
        
        for (int val : possible_values) {
            stats.decisions++;
            grid[r][c] = val;
            
            if (use_fc) {
                std::map<std::pair<int,int>, std::set<int>> domains_copy = domains;
                std::vector<std::pair<int,int>> changes;
                
                if (propagate_constraints(r, c, val, domains_copy, changes)) {
                    if (solve_recursive(use_mrv, use_fc, domains_copy)) {
                        return true;
                    }
                }
            } else {
                if (solve_recursive(use_mrv, use_fc, domains)) {
                    return true;
                }
            }
            
            stats.backtracks++;
            grid[r][c] = 0;
        }
        
        return false;
    }
};

class SudokuSATSolver {
public:
    int size;
    int box_size;
    int num_vars;
    std::vector<std::vector<int>> clauses;
    std::vector<int> assignment;
    std::vector<double> vsids_scores;
    
    struct {
    long long decisions = 0;
    long long unit_props = 0;
    long long backtracks = 0;
    long long backjumps = 0;        // ADD THIS
    } stats;

    // ADD THESE NEW MEMBERS:
    int current_level;
    std::vector<int> var_level;         // Decision level for each variable
    std::vector<int> decision_stack;    // Variables in decision order
    std::vector<int> trail;             // Trail of all assignments in order
    std::vector<int> trail_lim;         // Separator indices for each decision level

    SudokuSATSolver(const Grid& puzzle) 
        : size(puzzle.size), box_size(puzzle.box_size) {
        num_vars = size * size * size;
        assignment.resize(num_vars, -1);
        vsids_scores.resize(num_vars, 0.0);
        var_level.resize(num_vars, -1);
        generate_cnf(puzzle);
    }

    SudokuSATSolver(int num_vars,
                const std::vector<std::vector<int>>& cnf_clauses)
    : num_vars(num_vars),
      assignment(num_vars, -1),
      clauses(cnf_clauses),
      vsids_scores(num_vars, 0.0)
    {
    var_level.resize(num_vars, -1);
    current_level = 0;
    // Any other needed custom initialization
    }

    
    int var_index(int r, int c, int val) const {
        return r * size * size + c * size + (val - 1);
    }
    
    void generate_cnf(const Grid& puzzle) {
        clauses.clear();
        
        // R1: Each cell has at least one value
        for (int r = 0; r < size; ++r) {
            for (int c = 0; c < size; ++c) {
                std::vector<int> clause;
                for (int val = 1; val <= size; ++val) {
                    clause.push_back(var_index(r, c, val) + 1);
                }
                clauses.push_back(clause);
            }
        }
        
        // R2: Each cell has at most one value
        for (int r = 0; r < size; ++r) {
            for (int c = 0; c < size; ++c) {
                for (int v1 = 1; v1 <= size; ++v1) {
                    for (int v2 = v1 + 1; v2 <= size; ++v2) {
                        clauses.push_back({
                            -(var_index(r, c, v1) + 1),
                            -(var_index(r, c, v2) + 1)
                        });
                    }
                }
            }
        }
        
        // R3: Each row has each value exactly once
        for (int r = 0; r < size; ++r) {
            for (int val = 1; val <= size; ++val) {
                std::vector<int> clause;
                for (int c = 0; c < size; ++c) {
                    clause.push_back(var_index(r, c, val) + 1);
                }
                clauses.push_back(clause);
                
                // At most once
                for (int c1 = 0; c1 < size; ++c1) {
                    for (int c2 = c1 + 1; c2 < size; ++c2) {
                        clauses.push_back({
                            -(var_index(r, c1, val) + 1),
                            -(var_index(r, c2, val) + 1)
                        });
                    }
                }
            }
        }
        
        // R4: Each column has each value exactly once
        for (int c = 0; c < size; ++c) {
            for (int val = 1; val <= size; ++val) {
                std::vector<int> clause;
                for (int r = 0; r < size; ++r) {
                    clause.push_back(var_index(r, c, val) + 1);
                }
                clauses.push_back(clause);
                
                // At most once
                for (int r1 = 0; r1 < size; ++r1) {
                    for (int r2 = r1 + 1; r2 < size; ++r2) {
                        clauses.push_back({
                            -(var_index(r1, c, val) + 1),
                            -(var_index(r2, c, val) + 1)
                        });
                    }
                }
            }
        }
        
        // R5: Each box has each value exactly once
        for (int box_r = 0; box_r < box_size; ++box_r) {
            for (int box_c = 0; box_c < box_size; ++box_c) {
                for (int val = 1; val <= size; ++val) {
                    std::vector<int> clause;
                    std::vector<std::pair<int,int>> cells;
                    
                    for (int i = 0; i < box_size; ++i) {
                        for (int j = 0; j < box_size; ++j) {
                            int r = box_r * box_size + i;
                            int c = box_c * box_size + j;
                            clause.push_back(var_index(r, c, val) + 1);
                            cells.push_back({r, c});
                        }
                    }
                    clauses.push_back(clause);
                    
                    // At most once
                    for (size_t i = 0; i < cells.size(); ++i) {
                        for (size_t j = i + 1; j < cells.size(); ++j) {
                            clauses.push_back({
                                -(var_index(cells[i].first, cells[i].second, val) + 1),
                                -(var_index(cells[j].first, cells[j].second, val) + 1)
                            });
                        }
                    }
                }
            }
        }
        
        // R6: Pre-filled cells (initial clues)
        for (int r = 0; r < size; ++r) {
            for (int c = 0; c < size; ++c) {
                if (puzzle[r][c] != 0) {
                    clauses.push_back({var_index(r, c, puzzle[r][c]) + 1});
                }
            }
        }
    }
    
    bool solve(bool do_unit_prop, bool do_pure_literal, bool do_vsids) {
        stats = {0, 0, 0};
        assignment.assign(num_vars, -1);
        vsids_scores.assign(num_vars, 0.0);
        var_level.assign(num_vars, -1);          
        decision_stack.clear();                  
        trail.clear();                           // Clear trail
        trail_lim.clear();                       // Clear trail limits
        current_level = 0;
        
        if (do_vsids) {
            for (const auto& clause : clauses) {
                for (int lit : clause) {
                    vsids_scores[std::abs(lit) - 1] += 1.0;
                }
            }
        }
        
        return dpll(do_unit_prop, do_pure_literal, do_vsids);
    }
    
    Grid get_solution_grid() {
        Grid solution(size);
        for (int var = 0; var < num_vars; ++var) {
            if (assignment[var] == 1) {
                int r = var / (size * size);
                int c = (var / size) % size;
                int val = (var % size) + 1;
                solution[r][c] = val;
            }
        }
        return solution;
    }

private:
    bool is_clause_satisfied(const std::vector<int>& clause) const {
        for (int lit : clause) {
            int var = std::abs(lit) - 1;
            int required_val = (lit > 0) ? 1 : 0;
            if (assignment[var] == required_val) {
                return true;
            }
        }
        return false;
    }
    
    bool is_clause_conflicting(const std::vector<int>& clause) const {
        for (int lit : clause) {
            int var = std::abs(lit) - 1;
            if (assignment[var] == -1) {
                return false; // Has unassigned literal
            }
            int required_val = (lit > 0) ? 1 : 0;
            if (assignment[var] == required_val) {
                return false; // At least one literal is satisfied
            }
        }
        return true; // All literals are false = conflict
    }
    
    bool unit_propagate() {
        bool changed = true;
        while (changed) {
            changed = false;
            
            for (const auto& clause : clauses) {
                if (is_clause_satisfied(clause)) {
                    continue;
                }
                
                // Check for conflict
                if (is_clause_conflicting(clause)) {
                    return false; // Conflict found
                }
                
                // Check for unit clause
                int unassigned_lit = 0;
                int unassigned_count = 0;
                
                for (int lit : clause) {
                    int var = std::abs(lit) - 1;
                    if (assignment[var] == -1) {
                        unassigned_lit = lit;
                        unassigned_count++;
                    }
                }
                
                if (unassigned_count == 1) {
                    // Unit clause - must assign this literal to true
                    int var = std::abs(unassigned_lit) - 1;
                    assignment[var] = (unassigned_lit > 0) ? 1 : 0;
                    var_level[var] = current_level; 
                    trail.push_back(var);  // Add to trail
                    stats.unit_props++;
                    changed = true;
                }
            }
        }
        return true;
    }

    // Check for conflicts without propagating (used when unit propagation is disabled)
    bool check_for_conflicts() {
        for (const auto& clause : clauses) {
            if (is_clause_satisfied(clause)) {
                continue;
            }
            
            // Check if clause is completely falsified
            if (is_clause_conflicting(clause)) {
                return false; // Conflict found
            }
        }
        return true; // No conflicts
    }
    
    void pure_literal_assign() {
        std::map<int, int> polarity; // -1=both, 0=neg only, 1=pos only
        
        for (const auto& clause : clauses) {
            if (is_clause_satisfied(clause)) {
                continue;
            }
            
            for (int lit : clause) {
                int var = std::abs(lit) - 1;
                if (assignment[var] != -1) continue;
                
                int pol = (lit > 0) ? 1 : 0;
                
                if (polarity.find(var) == polarity.end()) {
                    polarity[var] = pol;
                } else if (polarity[var] != pol) {
                    polarity[var] = -1; // Both polarities seen
                }
            }
        }
        
        // FIX: Assign variables AND update their decision levels
        // Also add to trail for proper backtracking
        for (const auto& [var, pol] : polarity) {
            if (pol != -1 && assignment[var] == -1) {
                assignment[var] = pol;
                var_level[var] = current_level;
                trail.push_back(var);  // Add to trail
            }
        }
    }
    
    int pick_branch_var(bool do_vsids) {
        if (do_vsids) {
            double max_score = -1.0;
            int best_var = -1;
            
            for (int var = 0; var < num_vars; ++var) {
                if (assignment[var] == -1 && vsids_scores[var] > max_score) {
                    max_score = vsids_scores[var];
                    best_var = var;
                }
            }
            
            if (best_var != -1) return best_var;
        }
        
        for (int var = 0; var < num_vars; ++var) {
            if (assignment[var] == -1) return var;
        }
        
        return -1;
    }
    
    bool dpll(bool do_unit_prop, bool do_pure_literal, bool do_vsids) {
    // Unit propagation or conflict checking
    if (do_unit_prop) {
        if (!unit_propagate()) {
            // Conflict detected during unit propagation
            stats.backtracks++;
            return false;
        }
    } else {
        // When unit propagation is disabled, still check for conflicts
        if (!check_for_conflicts()) {
            stats.backtracks++;
            return false;
        }
    }
    
    // Pure literal elimination
    if (do_pure_literal) {
        pure_literal_assign();
    }
    
    // Check if all variables are assigned
    int var = pick_branch_var(do_vsids);
    if (var == -1) {
        // All assigned and no conflicts = SAT
        return true;
    }
    
    // Make a decision (branch)
    stats.decisions++;
    current_level++;
    decision_stack.push_back(var);
    
    // FIXED: Use trail-based backtracking
    // Remember where this decision level starts in the trail
    int trail_marker = trail.size();
    trail_lim.push_back(trail_marker);
    
    // Try assigning true
    assignment[var] = 1;
    var_level[var] = current_level;
    trail.push_back(var);  // Add decision to trail
    
    if (dpll(do_unit_prop, do_pure_literal, do_vsids)) {
        return true;
    }
    
    // FIXED: Undo all assignments made at this level using trail
    while (trail.size() > trail_marker) {
        int v = trail.back();
        trail.pop_back();
        assignment[v] = -1;
        var_level[v] = -1;
    }
    
    // Try assigning false
    assignment[var] = 0;
    var_level[var] = current_level;
    trail.push_back(var);  // Add to trail
    
    if (dpll(do_unit_prop, do_pure_literal, do_vsids)) {
        return true;
    }
    
    // Both failed - backtrack
    // Undo all assignments made at this level
    while (trail.size() > trail_marker) {
        int v = trail.back();
        trail.pop_back();
        assignment[v] = -1;
        var_level[v] = -1;
    }
    
    current_level--;
    decision_stack.pop_back();
    trail_lim.pop_back();
    
    // FIXED VSIDS: Bump score on conflicts (simplified - still not full VSIDS)
    if (do_vsids) {
        vsids_scores[var] += 1.0;
        // Decay all scores periodically
        if (stats.decisions % 100 == 0) {
            for (int v = 0; v < num_vars; ++v) {
                vsids_scores[v] *= 0.95;
            }
        }
    }
    
    return false;
}

    
    // NEW FUNCTION: Analyze conflict and find backjump level
    int analyze_conflict() {
        int backjump_level = 0;
        
        // Find the highest decision level involved in the conflict
        // (excluding the current level)
        for (const auto& clause : clauses) {
            if (is_clause_conflicting(clause)) {
                // Look at all variables in conflicting clauses
                for (int lit : clause) {
                    int var = std::abs(lit) - 1;
                    if (var >= 0 && var < num_vars && var_level[var] >= 0) {
                        // Find the second-highest level (to jump to)
                        if (var_level[var] < current_level && var_level[var] > backjump_level) {
                            backjump_level = var_level[var];
                        }
                    }
                }
            }
        }
        
        return backjump_level;
    }
    
    // NEW FUNCTION: Backtrack to a specific level
    void backtrack_to_level(int target_level) {
        // Unassign all variables assigned after target_level
        for (int var = 0; var < num_vars; ++var) {
            if (var_level[var] > target_level) {
                assignment[var] = -1;
                var_level[var] = -1;
            }
        }
        
        // Clean up decision stack
        while (!decision_stack.empty() && 
               (decision_stack.back() >= num_vars || 
                var_level[decision_stack.back()] > target_level)) {
            decision_stack.pop_back();
        }
        
        current_level = target_level;
    }

};


// ---- N-Queens CNF ----
struct NQueensCNF {
    int N;
    std::vector<std::vector<int>> clauses;
    NQueensCNF(int n) : N(n) { generate(); }
    int var(int i, int j) { return i*N + j + 1; }
    void generate() {
        for (int i = 0; i < N; ++i) {
            std::vector<int> clause;
            for (int j = 0; j < N; ++j) clause.push_back(var(i, j));
            clauses.push_back(clause); // At least one in row
        }
        for (int i = 0; i < N; ++i)
            for (int j1 = 0; j1 < N; ++j1)
                for (int j2 = j1+1; j2 < N; ++j2)
                    clauses.push_back({-var(i, j1), -var(i, j2)}); // At most one in row
        for (int j = 0; j < N; ++j)
            for (int i1 = 0; i1 < N; ++i1)
                for (int i2 = i1+1; i2 < N; ++i2)
                    clauses.push_back({-var(i1, j), -var(i2, j)}); // At most one in col
        for (int i1 = 0; i1 < N; ++i1)
            for (int j1 = 0; j1 < N; ++j1)
                for (int i2 = i1+1; i2 < N; ++i2) {
                    int j2a = j1 + (i2 - i1); // major diag
                    int j2b = j1 - (i2 - i1); // minor diag
                    if (j2a < N)
                        clauses.push_back({-var(i1, j1), -var(i2, j2a)});
                    if (j2b >= 0)
                        clauses.push_back({-var(i1, j1), -var(i2, j2b)});
                }
    }
};

// ---- Graph Coloring CNF ----
struct GraphColoringCNF {
    int V, K;
    std::vector<std::pair<int, int>> edges;
    std::vector<std::vector<int>> clauses;
    GraphColoringCNF(int V_, int K_, const std::vector<std::pair<int,int>>& E)
        : V(V_), K(K_), edges(E) { generate(); }
    int var(int v, int c) { return v*K + c + 1; }
    void generate() {
        for (int v = 0; v < V; ++v) {
            std::vector<int> clause;
            for (int c = 0; c < K; ++c) clause.push_back(var(v, c));
            clauses.push_back(clause); // at least one color
        }
        for (int v = 0; v < V; ++v)
            for (int c1 = 0; c1 < K; ++c1)
                for (int c2 = c1+1; c2 < K; ++c2)
                    clauses.push_back({-var(v, c1), -var(v, c2)}); // at most one color
        for (auto [u, v] : edges)
            for (int c = 0; c < K; ++c)
                clauses.push_back({-var(u, c), -var(v, c)}); // neighbors diff color
    }
};

std::vector<std::pair<int,int>> read_edges_from_file(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Failed to open edge file: " + filename);
    int V, E;
    file >> V >> E;
    std::vector<std::pair<int,int>> edges;
    for (int i = 0; i < E; i++) {
        int u, v;
        file >> u >> v;
        edges.emplace_back(u, v);
    }
    return edges;
}

// --- Main ---
int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage:\n";
        std::cerr << "  " << argv[0] << " <puzzle_file.txt>\n";
        std::cerr << "  " << argv[0] << " --nqueens [N]\n";
        std::cerr << "  " << argv[0] << " --gcolor <K> <edge_file>\n";
        return 1;
    }

    std::string arg1 = argv[1];

    if (arg1 == "--nqueens") {
        int N = (argc >= 3) ? std::stoi(argv[2]) : 8;
        NQueensCNF cnf(N);
        SudokuSATSolver solver(N*N, cnf.clauses);
        bool solved = solver.solve(true, true, true);
        if (solved) {
            std::cout << "N-Queens solution for N=" << N << ":\n";
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    std::cout << (solver.assignment[i*N + j] == 1 ? "Q " : ". ");
                }
                std::cout << "\n";
            }
        } else {
            std::cout << "No solution found.\n";
        }
        return 0;
    }

    if (arg1 == "--gcolor") {
        if (argc < 4) {
            std::cerr << "Usage: " << argv[0] << " --gcolor <K> <edge_file>\n";
            return 1;
        }
        int K = std::stoi(argv[2]);
        std::string edge_file = argv[3];
        std::vector<std::pair<int,int>> edges;
        try {
            edges = read_edges_from_file(edge_file);
        } catch (const std::exception &e) {
            std::cerr << e.what() << std::endl;
            return 1;
        }
        int max_vertex = -1;
        for (const auto& e : edges)
            max_vertex = std::max({max_vertex, e.first, e.second});
        int V = max_vertex + 1;
        GraphColoringCNF cnf(V, K, edges);
        SudokuSATSolver solver(V*K, cnf.clauses);
        bool solved = solver.solve(true, true, true);
        if (solved) {
            std::cout << "Graph Coloring solution:\n";
            for (int v = 0; v < V; ++v) {
                std::cout << "Vertex " << v << ": ";
                for (int c = 0; c < K; ++c) {
                    if (solver.assignment[v*K + c] == 1)
                        std::cout << (c + 1) << " ";
                }
                std::cout << "\n";
            }
        } else {
            std::cout << "No coloring possible.\n";
        }
        return 0;
    }
  // --- Original Sudoku Solver Mode ---
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <puzzle_file.txt>" << std::endl;
        std::cerr << "\nExample:" << std::endl;
        std::cerr << "  " << argv[0] << " test_puzzles/easy.txt" << std::endl;
        std::cerr << "\nFile format:" << std::endl;
        std::cerr << "  First line: grid size (e.g., 9)" << std::endl;
        std::cerr << "  Following lines: grid values (0 for empty)" << std::endl;
        return 1;
    }

    // Initialize CSVs and load puzzle
    init_csv_files();
    Grid puzzle = load_puzzle_from_file(argv[1]);
   
    // Extract puzzle name from filename
    std::string puzzle_name = argv[1];
    size_t last_slash = puzzle_name.find_last_of("/\\");
    if (last_slash != std::string::npos) {
        puzzle_name = puzzle_name.substr(last_slash + 1);
    }
    size_t last_dot = puzzle_name.find_last_of(".");
    if (last_dot != std::string::npos) {
        puzzle_name = puzzle_name.substr(0, last_dot);
    }
    
    int clues = count_clues(puzzle);

    // Open output file
    std::string output_filename = "results/" + puzzle_name + "_output.txt";
    std::ofstream outfile(output_filename);
    
    // Redirect cout to both console and file using a custom streambuf
    // For simplicity, we'll just write to file and console separately
    
   // Only print to file (for detailed logs)
    auto print_to_file_only = [&](const std::string& text) {
        outfile << text;
    };

    // Print only to terminal (for final concise output)
    auto print_to_terminal_only = [&](const std::string& text) {
        std::cout << text;
    };

    
    std::ostringstream header;
    header << std::string(70, '=') << "\n"
           << "       SUDOKU SOLVER BENCHMARK (C++)\n"
           << std::string(70, '=') << "\n"
           << "\nPuzzle: " << puzzle_name << "\n"
           << "Size: " << puzzle.size << "x" << puzzle.size << "\n"
           << "Clues: " << clues << " / " << (puzzle.size * puzzle.size) << "\n"
           << "\nInitial Puzzle:\n";
    print_to_file_only(header.str());
    
    // Print grid to both
    std::ostringstream grid_str;
    int n = puzzle.size;
    int box = puzzle.box_size;
    for (int i = 0; i < n + box + 1; ++i) grid_str << "-";
    grid_str << "\n";
    for (int r = 0; r < n; ++r) {
        grid_str << "|";
        for (int c = 0; c < n; ++c) {
            if (puzzle[r][c] == 0) {
                grid_str << " .";
            } else {
                grid_str << std::setw(2) << puzzle[r][c];
            }
            if ((c + 1) % box == 0) {
                grid_str << " |";
            }
        }
        grid_str << "\n";
        if ((r + 1) % box == 0 && r + 1 < n) {
            for (int i = 0; i < n + box + 1; ++i) grid_str << "-";
            grid_str << "\n";
        }
    }
    for (int i = 0; i < n + box + 1; ++i) grid_str << "-";
    grid_str << "\n";
    print_to_file_only(grid_str.str());

    // BACKTRACKING
    std::ostringstream bt_header;
    bt_header << "\n" << std::string(70, '=') << "\n"
              << "              BACKTRACKING SOLVER\n"
              << std::string(70, '=') << "\n";
    print_to_file_only(bt_header.str());

    std::vector<std::tuple<std::string, bool, bool>> bt_configs = {
        {"MRV + Forward Checking", true, true},
        {"MRV only", true, false},
        {"Forward Checking only", false, true},
        {"No heuristics", false, false},
    };

    struct BtResult { std::string name; double time; long long dec, back; bool success; };
    std::vector<BtResult> bt_results;

    for (const auto& [name, use_mrv, use_fc] : bt_configs) {
        BacktrackingSolver solver(puzzle);
        
        std::ostringstream config_output;
        config_output << "\n[" << name << "]\n";
        print_to_file_only(config_output.str());
        
        Timer t;
        bool solved = solver.solve(use_mrv, use_fc);
        double runtime = t.elapsed_s();
        
        if (solved) {
            std::ostringstream result;
            result << "  ✓ Solved in " << std::fixed << std::setprecision(6) << runtime << "s\n"
                   << "    Decisions: " << solver.stats.decisions 
                   << ", Backtracks: " << solver.stats.backtracks << "\n";
            print_to_file_only(result.str());
            
            bt_results.push_back({name, runtime, solver.stats.decisions, solver.stats.backtracks, true});
            export_backtracking_result(puzzle_name, clues, name, runtime,
                                       solver.stats.decisions, solver.stats.backtracks);
        } else {
            print_to_file_only("  ✗ Failed to solve\n");
            bt_results.push_back({name, runtime, 0, 0, false});
        }
    }
    
    // Print solution
    BacktrackingSolver temp_solver(puzzle);
    temp_solver.solve(true, true);
    
    std::ostringstream solution_output;
    solution_output << "\n  Solved Grid (Backtracking):\n";
    print_to_file_only(solution_output.str());
    
    // Print solved grid
    std::ostringstream solved_grid_str;
    for (int i = 0; i < n + box + 1; ++i) solved_grid_str << "-";
    solved_grid_str << "\n";
    for (int r = 0; r < n; ++r) {
        solved_grid_str << "|";
        for (int c = 0; c < n; ++c) {
            solved_grid_str << std::setw(2) << temp_solver.solution[r][c];
            if ((c + 1) % box == 0) {
                solved_grid_str << " |";
            }
        }
        solved_grid_str << "\n";
        if ((r + 1) % box == 0 && r + 1 < n) {
            for (int i = 0; i < n + box + 1; ++i) solved_grid_str << "-";
            solved_grid_str << "\n";
        }
    }
    for (int i = 0; i < n + box + 1; ++i) solved_grid_str << "-";
    solved_grid_str << "\n";
    print_to_file_only(solved_grid_str.str());

    if (verify_sudoku_solution(temp_solver.solution)) {
        print_to_file_only("  ✓ Backtracking solution verified as CORRECT!\n\n");
    } else {
        print_to_file_only("  ✗ Backtracking solution INVALID!\n\n");
    }

    // SAT SOLVER
    std::ostringstream sat_header;
    sat_header << "\n" << std::string(70, '=') << "\n"
               << "               SAT SOLVER (DPLL)\n"
               << std::string(70, '=') << "\n";
    print_to_file_only(sat_header.str());
    
    SudokuSATSolver initial_solver(puzzle);
    std::ostringstream cnf_info;
    cnf_info << "Generating CNF clauses...\n"
             << "Generated " << initial_solver.clauses.size() << " clauses for " 
             << initial_solver.num_vars << " variables.\n\n";
    print_to_file_only(cnf_info.str());
    
    struct SATConfig {
        std::string name;
        bool unit_prop;
        bool pure_literal;
        bool vsids;
    };
    
    std::vector<SATConfig> sat_configs = {
        {"All heuristics", true, true, true},
        {"Unit Prop + VSIDS", true, false, true},
        {"Unit Prop + Pure Literal", true, true, false},
        {"Unit Propagation only", true, false, false},
        {"VSIDS only", false, false, true},
        {"Pure Literal only", false, true, false},
        {"No heuristics", false, false, false},
    };

    struct SatResult { 
        std::string name; 
        double time; 
        long long dec, props, back, jumps; 
        bool success; 
        bool timed_out;
    };
    std::vector<SatResult> sat_results;
    
    Grid sat_solved_grid(puzzle.size);
    bool sat_solution_found = false;

    for (const auto& config : sat_configs) {
        std::ostringstream config_output;
        config_output << "\n[" << config.name << "]\n";
        print_to_file_only(config_output.str());
        
        Timer t;
        double runtime = 0.0;
        bool solved = false;
        bool timed_out = false;
        
        auto sat_solver_ptr = std::make_shared<SudokuSATSolver>(puzzle);
        auto stats_copy = sat_solver_ptr->stats;

        auto solve_task = [sat_solver_ptr, unit=config.unit_prop, pure=config.pure_literal, vsids=config.vsids]() -> bool {
            return sat_solver_ptr->solve(unit, pure, vsids);
        };

        auto future = std::async(std::launch::async, solve_task);
        std::chrono::seconds timeout_duration(30);
        std::future_status status = future.wait_for(timeout_duration);

        if (status == std::future_status::timeout) {
            print_to_file_only("  ✗ TIMEOUT (exceeded 30 seconds)\n");
            timed_out = true;
            runtime = 30.0;
            static std::vector<std::future<bool>> abandoned_futures;
            abandoned_futures.push_back(std::move(future));
        } else if (status == std::future_status::ready) {
            runtime = t.elapsed_s();
            try {
                solved = future.get();
                stats_copy = sat_solver_ptr->stats;
                if (solved) {
                    std::ostringstream result;
                    result << "  ✓ Solved in " << std::fixed << std::setprecision(6) << runtime << "s\n"
                           << "    Decisions: " << stats_copy.decisions 
                           << ", Unit Props: " << stats_copy.unit_props 
                           << ", Backtracks: " << stats_copy.backtracks
                           << ", Backjumps: " << stats_copy.backjumps << std::endl;
                    print_to_file_only(result.str());
                    
                    if (!sat_solution_found) {
                        sat_solved_grid = sat_solver_ptr->get_solution_grid();
                        sat_solution_found = true;
                    }
                } else {
                    print_to_file_only("  ✗ Failed to solve (Unsatisfiable)\n");
                }
            } catch (const std::exception& e) {
                std::ostringstream error;
                error << "  ✗ Error: " << e.what() << "\n";
                print_to_file_only(error.str());
                solved = false;
            }
        }

        export_sat_result(puzzle_name, clues, config.name, runtime,
                         solved ? stats_copy.decisions : 0,
                         solved ? stats_copy.unit_props : 0,
                         solved ? stats_copy.backtracks : 0,
                         solved ? stats_copy.backjumps : 0, 
                         timed_out);

        if (timed_out) {
            sat_results.push_back({config.name, runtime, 0, 0, 0, 0, false, true});    // ADD 0 for jumps
        } else if (solved) {
            sat_results.push_back({config.name, runtime, stats_copy.decisions, stats_copy.unit_props, stats_copy.backtracks, stats_copy.backjumps, true, false});    // ADD backjumps
        } else {
            sat_results.push_back({config.name, runtime, stats_copy.decisions, stats_copy.unit_props, stats_copy.backtracks, stats_copy.backjumps, false, false});    // ADD backjumps
        }

    }
    
    if (sat_solution_found) {
        print_to_file_only("\n  Solved Grid (SAT Solver):\n");
        
        std::ostringstream sat_grid_str;
        for (int i = 0; i < n + box + 1; ++i) sat_grid_str << "-";
        sat_grid_str << "\n";
        for (int r = 0; r < n; ++r) {
            sat_grid_str << "|";
            for (int c = 0; c < n; ++c) {
                sat_grid_str << std::setw(2) << sat_solved_grid[r][c];
                if ((c + 1) % box == 0) {
                    sat_grid_str << " |";
                }
            }
            sat_grid_str << "\n";
            if ((r + 1) % box == 0 && r + 1 < n) {
                for (int i = 0; i < n + box + 1; ++i) sat_grid_str << "-";
                sat_grid_str << "\n";
            }
        }
        for (int i = 0; i < n + box + 1; ++i) sat_grid_str << "-";
        sat_grid_str << "\n";
        print_to_file_only(sat_grid_str.str());

        if (verify_sudoku_solution(sat_solved_grid)) {
            print_to_file_only("  ✓ SAT solution verified as CORRECT!\n\n");
        } else {
            print_to_file_only("  ✗ SAT solution INVALID!\n\n");
        }
    }

    // SUMMARY
    std::ostringstream summary;
    summary << "\n" << std::string(70, '=') << "\n"
            << "                    SUMMARY\n"
            << std::string(70, '=') << "\n"
            << std::left << "\nBacktracking Solver:\n"
            << std::setw(30) << "Configuration" << std::setw(15) << "Time (s)" 
            << std::setw(15) << "Decisions" << "Backtracks\n"
            << std::string(70, '-') << "\n";
    
    for(const auto& res : bt_results) {
        summary << std::setw(30) << res.name 
                << std::setw(15) << std::fixed << std::setprecision(6) << res.time 
                << std::setw(15) << res.dec << res.back << "\n";
    }
    
    summary << "\nSAT Solver:\n"
        << std::setw(30) << "Configuration" << std::setw(15) << "Time (s)" 
        << std::setw(15) << "Decisions" << std::setw(15) << "Unit Props" 
        << std::setw(12) << "Backtracks" << "Backjumps\n"
        << std::string(90, '-') << "\n";    // Changed 70 to 90
    
    for(const auto& res : sat_results) {
    if (res.timed_out) {
        summary << std::setw(30) << res.name 
                << std::setw(15) << "> 30.000000"
                << std::setw(15) << "N/A" << std::setw(15) << "N/A" 
                << std::setw(12) << "N/A" << "N/A\n";    // ADDED one more N/A
    } else if (res.success) {
        summary << std::setw(30) << res.name 
                << std::setw(15) << std::fixed << std::setprecision(6) << res.time 
                << std::setw(15) << res.dec << std::setw(15) << res.props 
                << std::setw(12) << res.back << res.jumps << "\n";    // ADD res.jumps
    } else {
         summary << std::setw(30) << res.name 
                << std::setw(15) << std::fixed << std::setprecision(6) << res.time 
                << std::setw(15) << res.dec << std::setw(15) << res.props 
                << std::setw(12) << res.back << res.jumps << " (FAILED)\n";    // ADD res.jumps
    }
    }

    
    summary << "\n" << std::string(70, '=') << "\n"
            << "✓ RESULTS SAVED\n"
            << std::string(70, '=') << "\n"
            << "  Output:       " << output_filename << "\n"
            << "  Backtracking: results/backtracking_results.csv\n"
            << "  SAT Solver:   results/sat_results.csv\n"
            << std::string(70, '=') << "\n";
    
    print_to_file_only(summary.str());

SudokuSATSolver sat_solver(puzzle);
bool sat_solved = sat_solver.solve(true, true, true);
if (sat_solved) {
    std::cout << "SAT solver found a solution.\n";
} else {
    std::cout << "No solution found by SAT solver.\n";
}

 
outfile.close();
    
    return 0;
}
