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
	string input_dir = "./data/E2F200.fa"; // sample data file
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
	double rho = 0.01;
	int l = 11;
	int d = 11;
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

		// Choose the potential solutions with best Hamming scores
		auto best_colony_motifs = FilterBestMotifs(colony_motifs, 0);

		// Apply local search on each motif candidate
		for (const Motif& p : best_colony_motifs) {
			auto ls_result = LocalSearchRACO(input, p);
			int best_score = numeric_limits<int>::max();
			for (const auto& motif : ls_result) {
				all_motif_candidates.emplace_back(motif);
				best_score = min(best_score, motif.first);
			}
			for (const auto& motif : ls_result) {
				if (motif.first == best_score) {
					UpdatePheromoneVertex(pheromone_vertex, motif.second, rho);
					UpdatePheromoneEdge(pheromone_edge, motif.second, rho);
				}
			}
		}

		cout << "Ending iteration " << loop + 1 << '\n';
	}
	cout << "Ending.\n";
	cout << "No. of candidates before relax = " << all_motif_candidates.size() << '\n';

	vector<tuple<double, Motif, vector<int>>> best_relaxed_solutions;

	for (const Motif& motif : all_motif_candidates) {
		auto relaxed_solutions = Relax(input, motif, eta_vertex, eta_edge);
		assert(!relaxed_solutions.empty());
		auto pos = relaxed_solutions[0].second;
		best_relaxed_solutions.emplace_back(-relaxed_solutions[0].first, FindConsensus(input, l, pos), pos);
	}

	sort(best_relaxed_solutions.begin(), best_relaxed_solutions.end());


	// Print the stats of true msa
	// vector<int> v = {61,55,76,63,50,7,42,39,9,14,61,41,48,71,17,53,84,78};
	// vector<int> v = {95,95,95,95,95,95,95,95,95,95,95,95,95,95,95,95,95,95,95,95,95,95,95,95,95};
	vector<int> v = {95,95,95,126,95,95,95,95,95,95,95,95,95,95,95,46,95,95,95,95,95,97,95,95,27};
	for (int& x : v) x--;
	{
		os << "TRUE SOLUTION\n";
		auto motif = FindConsensus(input, l, v);
		os << "\nMotif: " << motif.second << '\n';

		// Print positions of instances
		os << "Position: ";
		for (int p : v) {
			os << p + 1 << ' ';
		}
		os << '\n';

		// Print Hamming Score
		os << "Hamming Score (HSc) = " << motif.first << '\n';

		// Print Concensus Score 
		os << "Concensus Score (CSc) = " << (long long) input.size() * l - motif.first << '\n';

		// Print IC Score 1
		os << "IC Score 1 (IC1) = " << fixed << setprecision(3) << GetICScore(input, eta_vertex, l, v) << '\n';

		// Print IC Score 2
		os << "IC Score 2 (IC2) = " << fixed << setprecision(3) << GetICScore2(input, eta_edge, l, v) << '\n';

		os << "------------------------------\n\n";
	}


	for (const auto& relaxed_solution : best_relaxed_solutions) {
		double score;
		Motif motif;
		vector<int> pos;
		tie(score, motif, pos) = relaxed_solution;

		static int num_sols = 0;
		os << "Solution #" << ++num_sols << '\n';	
		// Print motif
		os << "\nMotif: " << motif.second << '\n';

		// Print positions of instances
		os << "Position: ";
		for (int p : pos) {
			os << p + 1 << ' ';
		}
		os << '\n';
		assert(input.size() * 1LL * l - motif.first == -score);

		// Print Hamming Score
		os << "Hamming Score (HSc) = " << motif.first << '\n';

		// Print Concensus Score 
		os << "Concensus Score (CSc) = " << (long long) input.size() * l - motif.first << '\n';

		// Print IC Score 1
		os << "IC Score 1 (IC1) = " << fixed << setprecision(3) << GetICScore(input, eta_vertex, l, pos) << '\n';

		// Print IC Score 2
		os << "IC Score 2 (IC2) = " << fixed << setprecision(3) << GetICScore2(input, eta_edge, l, pos) << '\n';

		os << "------------------------------\n";
	}

	// Print time of execution
	double end_time = clock();
	double time_of_execution = (end_time - start_time) / CLOCKS_PER_SEC;
	os << "Time of execution - " << fixed << setprecision(3) << time_of_execution << " second(s)\n";
	cout << time_of_execution << '\n';
	os.close();

	return 0;
}