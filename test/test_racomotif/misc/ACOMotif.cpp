#include "ACO.h"
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

using namespace std;

using Motif = pair<int, string>;

void ReadFile(vector<string>& input) {
	string input_dir = "./data/CRP.fasta"; // sample data file
	ifstream is;
	is.open(input_dir);
	string str;
	while (getline(is, str)) {
		for (char& c : str) {
			if ('a' <= c && c <= 'z') {
				c = c - 'a' + 'A';
			}
		}
		bool is_valid = !str.empty();
		for (char c : str) {
			is_valid &= (c == 'A' || c == 'C' || c == 'G' || c == 'T');
			if (!is_valid) break;
		}
		if (is_valid) {
			input.emplace_back(str);	
		}
	}
	is.close();
}

int main(int argc, char* argv[]) {
	// Set random seed
	srand(time(NULL));
	
	// TODO: Handling flags and arguments
	int num_ants = 10;
	int num_loops = 100;
	double rho = 0.02;
	int l = 22;
	int d = 15;
	double alpha = 1.0;
	double beta = 1.0;
	int delta = 15;

	// Set starting time
	cout << "Starting...\n";
	double start_time = clock();

	// Read file
	vector<string> input;
	ReadFile(input);

	// Set output stream
	string output_dir = "./result/group1_cpp.txt";
	ofstream os;
	os.open(output_dir);

	// Initialize heuristic information
	cout << "Initializing heuristic information...\n";
	vector<double> eta_vertex;
	vector<vector<double>> eta_edge;
	tie(eta_vertex, eta_edge) = InitializeHeuristicInformation(input);

	// Initialize pheromone
	vector<double> pheromone_vertex(4, 1);
	vector<vector<double>> pheromone_edge(l - 1, vector<double>(16, 1));

	vector<Motif> all_motif_candidates;

	// Expand solution candidates set
	for (int loop = 0; loop < num_loops; loop++) { // A loop of ant colony
		vector<Motif> colony_motifs;

		for (int ant = 0; ant < num_ants; ant++) { // Each ant constructs a solution
			auto motif = ConstructSolution(eta_vertex, eta_edge, pheromone_vertex, pheromone_edge, l, alpha, beta);
			auto score = GetHammingScore(input, motif);
			colony_motifs.emplace_back(score, motif);
		}

		if (colony_motifs.empty()) { // if no solution is obtained, construct it again
			loop--;
			continue;
		}

		// Filter motifs that have the best (minimum) score
		auto best_colony_motifs = FilterBestMotifs(colony_motifs, 0);

		// Apply local search on each motif candidate
		for (const Motif& p : best_colony_motifs) {
			auto ls_result = LocalSearchACO(input, p);
			auto best_motif_overall = ls_result.first;
			auto best_motifs_ls = ls_result.second;

			all_motif_candidates.emplace_back(best_motif_overall);

			for (const auto& motif : best_motifs_ls) {
				all_motif_candidates.emplace_back(motif);
			}

			UpdatePheromoneVertex(pheromone_vertex, best_motif_overall.second, rho);
			UpdatePheromoneEdge(pheromone_edge, best_motif_overall.second, rho);
		}

		cout << "Ending iteration " << loop + 1 << '\n';
	}
	cout << "Ending.\n";

	sort(all_motif_candidates.begin(), all_motif_candidates.end());

	vector<vector<Motif>> best_motifs;
	
	// Filter delta best solutions	
	int ptr = 0;
	for (int t = 0; t < delta; t++) {
		best_motifs.emplace_back(vector<Motif>());
		bool flag = false;

		for (int i = ptr; i < (int) all_motif_candidates.size(); i++) {
			const auto& motif = all_motif_candidates[ptr];
			const auto& cur_motif = all_motif_candidates[i];
			
			if (cur_motif.first == motif.first) { // If current motif has the same score with those in the group, add it to the group
				bool is_duplicate = false;
				for (const auto& included_motif : best_motifs.back()) {
					if (cur_motif.second == included_motif.second) {
						is_duplicate = true;
						break;
					}
				}

				if (!is_duplicate) {
					best_motifs.back().emplace_back(cur_motif);
				}

			} else { // Else, it will be the start of next group
				ptr = i;
				break;
			}

			if (i == (int) all_motif_candidates.size() - 1) {
				flag = true;
			}
		}

		if (flag) {
			break;
		}
	}
	assert((int) best_motifs.size() <= delta);

	// Check d condition and print solution
	int num_solutions = 0;
	for (int t = 0; t < (int) best_motifs.size(); t++) {
		for (const auto& motif : best_motifs[t]) {
			auto distances = FindDistances(input, motif.second);
			
			bool is_valid = true;
			for (int dist : distances) {
				if (dist > d) {
					is_valid = false;
					break;
				}
			}

			if (is_valid) {
				auto result = FindMSA(input, motif.second);
				vector<int> a(input.size());
				assert(a.size() == input.size());
				for (int i = 0; i < (int) input.size(); i++) {
					a[i] = result[i][rand() % result[i].size()];
				}

				os << "Solution #" << ++num_solutions << '\n';	
				// Print motif
				os << "\nMotif: " << motif.second << '\n';

				// Print positions of instances
				os << "Position: ";
				for (int position : a) {
					os << position + 1 << ' ';
				}
				os << '\n';

				// Print Concensus Score 
				os << "Concensus Score (CSc) = " << (long long) input.size() * l - motif.first << '\n';

				// Print Hamming Score
				os << "Hamming Score (HSc) = " << motif.first << '\n';
				os << "------------------------------\n";
			}
		}
	}

	// Print time of execution
	double end_time = clock();
	double time_of_execution = (end_time - start_time) / CLOCKS_PER_SEC;
	os << "Time of execution - " << fixed << setprecision(3) << time_of_execution << " second(s)\n";
	cout << time_of_execution << '\n';
	os.close();

	return 0;
}