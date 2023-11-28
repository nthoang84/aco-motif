/*
	Compile:
		g++ relax.cpp -O2 -std=c++11 -o relax
		g++ acomotif_hm.cpp -O2 -std=c++11 -o acomotif_hm
	
	Run:
		acomotif_hm <dataset_name> <num_ants> <num_loops> <rho> <l> <d>
*/
#include "aco.h"
#include <bits/stdc++.h>

using namespace std;
using motif_t = pair<double, string>;

int main(int argc, char* argv[]) {
	srand(time(NULL));
	
	string dataset_name = argv[1];
	int num_ants = atoi(argv[2]);
	int num_loops = atoi(argv[3]);
	double rho = atof(argv[4]);
	int l = atoi(argv[5]);
	int d = atoi(argv[6]);
	double alpha = 1.0;
	double beta = 1.0;
	int delta = 10;

	cerr << "Starting...\n";
	double start_time = clock();

	vector<string> dataset;
	vector<string> sequence_names;
	read_file(dataset, sequence_names, dataset_name);
	filter_dust(dataset, sequence_names, dataset_name);

	cerr << "Initializing heuristic information...\n";
	vector<double> eta_vertex;
	vector<vector<double>> eta_edge;
	tie(eta_vertex, eta_edge) = get_heuristic_info(dataset);

	vector<double> pheromone_vertex(4, 1);
	vector<vector<double>> pheromone_edge(l - 1, vector<double>(16, 1));

	vector<motif_t> all_motif_candidates;
	cerr << "Constructing solutions...\n";
	for (int loop = 0; loop < num_loops; loop++) {
		vector<pair<motif_t, vector<int>>> colony_motifs;
		for (int ant = 0; ant < num_ants; ant++) {
			auto motif = ConstructSolution(eta_vertex, eta_edge, pheromone_vertex, pheromone_edge, l, alpha, beta);
			auto pos = FindMSAFixed(dataset, motif);
			auto score = GetICScore(dataset, eta_vertex, l, pos);
			colony_motifs.emplace_back(make_pair(score, motif), pos);
		}
		if (colony_motifs.empty()) {
			loop--;
			continue;
		}
		for (const auto& p : colony_motifs) {
			all_motif_candidates.emplace_back(p.first);
		}
		
		double best_score = 0;
		vector<pair<motif_t, vector<int>>> best_colony_motifs;
		for (const auto& p : colony_motifs) {
			if (p.first.first > best_score) {
				best_score = p.first.first;
				best_colony_motifs.clear();
				best_colony_motifs.emplace_back(p);
			} else if (p.first.first == best_score) {
				best_colony_motifs.emplace_back(p);
			}
		}
		for (const auto& p : best_colony_motifs) {
			motif_t motif1 = p.first;
			vector<int> pos1 = p.second;
			auto motif2 = FindConsensus(dataset, l, pos1);
			auto ls_1 = LocalSearchRACO(dataset, motif1, eta_vertex);
			auto all_ls_1 = ls_1.first;
			auto best_ls_1 = ls_1.second;
			all_motif_candidates.emplace_back(best_ls_1);
			for (const auto& p : all_ls_1) {
				all_motif_candidates.emplace_back(p);
			}
			
			if (motif2.second != best_ls_1.second) {
				auto ls_2 = LocalSearchRACO(dataset, motif2, eta_vertex);
				auto all_ls_2 = ls_2.first;
				all_motif_candidates.emplace_back(motif2);
				for (const auto& p : all_ls_2) {
					all_motif_candidates.emplace_back(p);
				}
			}
			auto optimal_motif = (motif2.first >= best_ls_1.first ? motif2 : best_ls_1);
			UpdatePheromoneVertex(pheromone_vertex, optimal_motif.second, rho);
			UpdatePheromoneEdge(pheromone_edge, optimal_motif.second, rho);
			cerr << GetICScore(dataset, eta_vertex, l, FindMSAFixed(dataset, optimal_motif.second)) << '\n';
		}
		cerr << "Ending iteration " << loop + 1 << '\n';
	}
	sort(all_motif_candidates.begin(), all_motif_candidates.end(), greater<Motif>());
	cerr << "Construction ended.\n";
	cerr << "Number of candidates before filtering: " << all_motif_candidates.size() << '\n';
	string output_dir = file_name + ".txt";
	ofstream os;
	os.open(output_dir, ofstream::out | ofstream::trunc);
	os.close();
	return 0;	
}