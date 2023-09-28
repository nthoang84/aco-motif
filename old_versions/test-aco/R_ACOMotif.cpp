/*
	Compile:
		g++ relax.cpp -O2 -std=c++11 -o relax
		g++ ACOMotif.cpp -O2 -std=c++11 -o ACOMotif
	
	Run:
		R_ACOMotif <dataset name> <no. of ants> <no. of colonies> <rho> <l> <d> <order to sort the results>

	Order to sort by: 
		- 0: alphabetical order of the positions
		- 1: alphabetical order of motif string
		- 2: descending order of Consensus Score
		- 3: ascending order of Hamming Score
		- 4: descending order of IC Score 1
		- 5: descending order of IC SCore 2
*/
#include "ACO.h"
#include <bits/stdc++.h>

using namespace std;
using Motif = pair<double, string>;

double GetPearsonCorrelation(const vector<vector<int>>& va, const vector<vector<int>>& vb) {
	int n = va.size() * va[0].size();
	vector<int> a;
	vector<int> b;
	for (const auto& v : va) for (int x : v) a.push_back(x);
	for (const auto& v : vb) for (int x : v) b.push_back(x);
	assert(a.size() == n && b.size() == n);
	int sum_a = accumulate(a.begin(), a.end(), 0);
	int sum_b = accumulate(b.begin(), b.end(), 0);
	double avg_a = sum_a * 1.0 / n;
	double avg_b = sum_b * 1.0 / n;
	double nume = 0;
	for (int i = 0; i < n; i++) {
		nume += (a[i] - avg_a) * (b[i] - avg_b);
	}
	double denom = 1;
	for (int t = 0; t < 2; t++) {
		const auto& v = (t == 0 ? a : b);
		const double& avg_v = (t == 0 ? avg_a : avg_b);
		double sum = 0;
		for (int i = 0; i < n; i++) {
			sum += pow(v[i] - avg_v, 2);
		}
		denom *= sqrt(sum);
	}
	return nume / denom;
}


int main(int argc, char* argv[]) {
	// Set random seed
	srand(time(NULL));
	
	// Set parameters
	string file_name, file_format;
	char* input_dir = argv[1];
	for (int i = 0, flag = 0; i < strlen(input_dir); i++) {
		if (input_dir[i] == '.') flag = 1;
		string& s = (flag ? file_format : file_name);
		s += input_dir[i];
	}
	int num_ants = atoi(argv[2]);
	int num_loops = atoi(argv[3]);
	double rho = atof(argv[4]);
	int l = atoi(argv[5]);
	int d = atoi(argv[6]);
	int sort_type = atoi(argv[7]);
	double alpha = 1.0;
	double beta = 1.0;
	int delta = 10;

	// Set starting time
	cerr << "Starting...\n";
	double start_time = clock();

	// Read file
	vector<string> input;
	vector<string> input_name;
	ReadFile(input, file_name, file_format, input_name);
	FilterDust(input, input_name, file_name);


	// Initialize heuristic information
	cerr << "Initializing heuristic information...\n";
	vector<double> eta_vertex;
	vector<vector<double>> eta_edge;
	tie(eta_vertex, eta_edge) = InitializeHeuristicInformation(input);

	// Initialize pheromone
	vector<double> pheromone_vertex(4, 1);
	vector<vector<double>> pheromone_edge(l - 1, vector<double>(16, 1));

	// vector<vector<int>> true_msa = {
	// 	{95},
	// 	{95},
	// 	{95},
	// 	{63, 126},
	// 	{95},
	// 	{95},
	// 	{95},
	// 	{95},
	// 	{95},
	// 	{95},
	// 	{95},
	// 	{95},
	// 	{95},
	// 	{95},
	// 	{95},
	// 	{46,144},
	// 	{95},
	// 	{95},
	// 	{95},
	// 	{95},
	// 	{95},
	// 	{97},
	// 	{95},
	// 	{95},
	// 	{27}
	// };

	// vector<vector<int>> true_msa = {
	// 	{17,61},
	// 	{17,55},
	// 	{76},
	// 	{63},
	// 	{50},
	// 	{7, 60},
	// 	{42},
	// 	{39},
	// 	{9,80},
	// 	{14},
	// 	{61},
	// 	{41},
	// 	{48},
	// 	{71},
	// 	{17},
	// 	{53},
	// 	{1,84},
	// 	{78}
	// };

	vector<vector<int>> true_msa = {
		{95},
		{95},
		{95},
		{95},
		{95},
		{95},
		{95},
		{95},
		{95},
		{95},
		{95},
		{95},
		{95},
		{95},
		{95},
		{95},
		{95},
		{95},
		{95},
		{95},
		{95},
		{95},
		{95},
		{95},
		{95}
	};

	assert(true_msa.size() == input.size());

	for (auto& pos : true_msa) for (int& p : pos) p--;
	vector<vector<int>> true_freq(l, vector<int>(4));
	for (int i = 0; i < true_msa.size(); i++) {
		int start = true_msa[i][rand() % true_msa[i].size()];
		for (int j = start; j < start + l; j++) {
			true_freq[j - start][Encode(input[i][j])]++;
		}
	}

	vector<Motif> all_motif_candidates;

	// Expand solution candidates set
	cerr << "Constructing solutions...\n";
	for (int loop = 0; loop < num_loops; loop++) { // A loop of ant colony
		vector<pair<Motif, vector<vector<int>>>> colony_motifs;

		for (int ant = 0; ant < num_ants; ant++) { // Each ant constructs a solution
			auto motif = ConstructSolution(eta_vertex, eta_edge, pheromone_vertex, pheromone_edge, l, alpha, beta);
			auto pos = FindMSA(input, motif);
			vector<int> pp(pos.size());
			for (int i = 0; i < pos.size(); i++) {
				pp[i] = pos[i][rand() % pos[i].size()];
			}
			auto score = GetICScore(input, eta_vertex, l, pp);
			
			colony_motifs.emplace_back(make_pair(score, motif), pos);
		}

		if (colony_motifs.empty()) { // if no solution is obtained, construct it again
			loop--;
			continue;
		}

		for (const auto& p : colony_motifs) {
			all_motif_candidates.emplace_back(p.first);
		}

		// Choose the potential solutions with best Hamming scores
		// auto best_colony_motifs = FilterBestMotifs(colony_motifs, 0);
		sort(colony_motifs.begin(), colony_motifs.end(), [&](pair<Motif, vector<vector<int>>> x, pair<Motif, vector<vector<int>>> y) {
			return x.first.first > y.first.first;
		});
		vector<pair<Motif, vector<vector<int>>>> best_colony_motifs = {colony_motifs[0]};

		// Apply local search on each motif candidate
		for (const auto& p : best_colony_motifs) {
			Motif motif1 = p.first;
			vector<vector<int>> pos1 = p.second;

			vector<int> normalized_pos;
			for (int i = 0; i < (int) pos1.size(); i++) {
				int selected_pos = pos1[i][rand() % pos1[i].size()];
				normalized_pos.emplace_back(selected_pos);
			}

			auto motif2 = FindConsensus(input, l, normalized_pos);

			auto ls_1 = LocalSearchRACO(input, motif1, eta_vertex);
			auto all_ls_1 = ls_1.first;
			auto best_ls_1 = ls_1.second;

			all_motif_candidates.emplace_back(best_ls_1);
			for (const auto& p : all_ls_1) {
				all_motif_candidates.emplace_back(p);
			}
			
			if (motif2.second != best_ls_1.second) {
				auto ls_2 = LocalSearchRACO(input, motif2, eta_vertex);
				auto all_ls_2 = ls_2.first;
				all_motif_candidates.emplace_back(motif2);
				for (const auto& p : all_ls_2) {
					all_motif_candidates.emplace_back(p);
				}
			}

			auto optimal_motif = (motif2.first <= best_ls_1.first ? motif2 : best_ls_1);
			UpdatePheromoneVertex(pheromone_vertex, optimal_motif.second, rho);
			UpdatePheromoneEdge(pheromone_edge, optimal_motif.second, rho);


			auto optimal_msa = FindMSA(input, optimal_motif.second);
			vector<vector<int>> freq(l, vector<int>(4));
			vector<int> pp;
			for (int i = 0; i < optimal_msa.size(); i++) {
				int start = optimal_msa[i][rand() % optimal_msa[i].size()];
				pp.push_back(start);
				for (int j = start; j < start + l; j++) {
					freq[j - start][Encode(input[i][j])]++;
				}
			}
			cerr << GetICScore(input, eta_vertex, l, pp) << '\n';

			cerr << fixed << setprecision(4) << "Loop #" << loop + 1 << ": " << GetPearsonCorrelation(freq, true_freq) << '\n';
		}

		cerr << "Ending iteration " << loop + 1 << '\n';
	}

	sort(all_motif_candidates.begin(), all_motif_candidates.end());
	cerr << "Construction ended.\n";

	return 0;

	/*
	vector<pair<double, string>> check_best_motifs;
	for (const auto& p : all_motif_candidates) {
		string motif = p.second;
		auto msa = FindMSA(input, motif);
		vector<int> pos(input.size());
		for (int i = 0; i < msa.size(); i++) {
			pos[i] = msa[i][rand() % msa[i].size()];
		}
		check_best_motifs.emplace_back(-GetICScore(input, eta_vertex, l, pos), motif);
	}
	sort(check_best_motifs.begin(), check_best_motifs.end());

	cerr << "TRUE MOTIF:\n";
	for (string s : true_consensus) {
		assert(s.size() == l);
		cerr << GetICScore(input, eta_vertex, l, pp) << ' ' << s << '\n';
	}
	cerr << '\n';

	ofstream oss("check.txt");
	for (auto p : check_best_motifs) {
		oss << fixed << setprecision(3) << -p.first << ' ' << p.second << ' ';
		// oss << "Diff: ";
		// for (auto true_motif : true_consensus) {
		// 	oss << GetHammingDistance(true_motif, p.second) << ' ';
		// }
		oss << '\n';
	}
	oss.close();
	*/

	cerr << "Number of candidates before filtering: " << all_motif_candidates.size() << '\n';
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

	// Create new output file or clear existing output file
	string output_dir = file_name + ".txt";
	ofstream os;
	os.open(output_dir, ofstream::out | ofstream::trunc);
	os.close();

	// Relax
	cerr << "Relaxing...\n";
	for (int t = 0; t < (int) best_motifs.size(); t++) {
		for (const auto& motif : best_motifs[t]) {
			system(("./relax data/" + file_name + file_format + " " + motif.second + " 1 1 >> " + file_name + ".txt").c_str());
		}
	}
	cerr << "Relax ended.\n";

	// Read motifs after relax
	cerr << "Retrieving solutions after relax...\n";
	string motifs_dir = file_name + ".txt";
	ifstream is;
	is.open(motifs_dir);

	using Result = tuple<vector<int>, string, int, int, double, double, string>; 
	vector<Result> results; 
	vector<int> pos;
	string motif, origin;
	int csc, hsc;
	string str;
	while (getline(is, str)) {
		if (str[0] == '-') {
			results.emplace_back(pos, motif, csc, hsc, GetICScore(input, eta_vertex, l, pos), GetICScore2(input, eta_edge, l, pos), origin);
			continue;
		}
		istringstream iss(str);
		if (str.find("Positions") != string::npos) {
			pos.clear();
			iss >> str;
			int p;
			while (iss >> p) pos.push_back(p - 1);
			assert(pos.size() == input.size());
			continue;
		}
		if (str.find("Motif") != string::npos) {
			iss >> str;
			iss >> motif;
			continue;
		}
		if (str.find("CSc") != string::npos) {
			iss >> str;
			iss >> csc;
			continue;
		}
		if (str.find("HSc") != string::npos) {
			iss >> str;
			iss >> hsc;
			continue;
		}
		if (str.find("Origin") != string::npos) {
			iss >> str;
			iss >> origin;
			continue;
		}
	}
	is.close();
	cerr << "Retrieving ended.\n";

	// Print results sorted by the pre-defined order
	cerr << "Printing solutions...\n";
	os.open(output_dir);
	sort(results.begin(), results.end(), [&](Result x, Result y) {
		if (sort_type == 0) return get<0>(x) < get<0>(y);
		if (sort_type == 1) return get<1>(x) < get<1>(y);
		if (sort_type == 2) return get<2>(x) > get<2>(y);
		if (sort_type == 3) return get<3>(x) < get<3>(y);
		if (sort_type == 4) return get<4>(x) > get<4>(y);
		if (sort_type == 5) return get<5>(x) > get<5>(y);
	});
	int num_solutions = 0;
	for (const auto& result : results) {
		os << "Solution #" << ++num_solutions << '\n';
		os << "Positions: ";
		for (int x : get<0>(result)) os << x + 1 << ' ';
		os << '\n';
		os << "Motif: " << get<1>(result) << '\n';
		os << "CSc: " << get<2>(result) << '\n';
		os << "HSc: " << get<3>(result) << '\n';
		os << "ICSc1: " << get<4>(result) << '\n';
		os << "ICSc2: " << get<5>(result) << '\n';
		os << "Origin: " << get<6>(result) << '\n';
		os << "------------------------------\n";
	}
	cerr << "Printing ended.\n";

	// Print time of execution
	double end_time = clock();
	double time_of_execution = (end_time - start_time) / CLOCKS_PER_SEC;
	cerr << "Time of execution: " << fixed << setprecision(3) << time_of_execution << " second(s)\n";
	os << "Time of execution: " << fixed << setprecision(3) << time_of_execution << " second(s)\n";
	os.close();
	return 0;
}