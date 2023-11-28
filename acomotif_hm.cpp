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
}