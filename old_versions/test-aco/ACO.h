#include <mpfr.h>
#include <bits/stdc++.h>

using namespace std;
using Motif = pair<double, string>;

const string alphabet = "ACGT";

// Map char values to int values
inline int Encode(char c) {
	return (c == 'A' ? 0 : (c == 'C' ? 1 : (c == 'G' ? 2 : 3)));
}

// Read data from dataset file
void ReadFile(vector<string>& dataset, const string& dataset_name, const string& dataset_format, vector<string>& seq_names) {
	string dataset_dir = "./data/" + dataset_name + dataset_format;
	ifstream is;
	is.open(dataset_dir);
	string line;
	while (getline(is, line)) {
		if (line[0] == '>') {
			dataset.emplace_back("");
			istringstream iss(line);
			seq_names.emplace_back("");
			iss >> seq_names.back();
			assert(seq_names.back()[0] == '>');
			seq_names.back() = seq_names.back().substr(1);
		} else if (!line.empty()) {
			for (char& c : line) {
				if ('a' <= c && c <= 'z') {
					c = c - 'a' + 'A';
				}
			}
			dataset.back() += line;
		}
	}
	is.close();
}

// Filter impossible segments in dataset, using Dust program
void FilterDust(vector<string>& dataset, const vector<string>& seq_names, const string& dataset_name) {
	ifstream is;
	is.open("data/" + dataset_name + ".dust");
	string line;
	while (getline(is, line)) {
		istringstream iss(line);
		string name;
		int start, finish;
		iss >> name >> start >> finish;
		for (int i = 0; i < (int) dataset.size(); i++) {
			if (seq_names[i] == name) {
				for (int j = start; j <= finish; j++) {
					dataset[i][j] = '$';
				}
				break;
			}
		}
	}
	is.close();
}

// Get Hamming distance of 2 strings with identical length (can be further improved?)
int GetHammingDistance(const string& s, const string& t) {
	assert(s.size() == t.size());
	int dist = 0;
	for (int i = 0; i < (int) s.size(); i++) {
		dist += (s[i] != t[i]);
	}
	return dist;
}

// Get Hamming score for a given motif
int GetHammingScore(const vector<string>& dataset, const string& motif) {
	int l = (int) motif.size();
	int score = 0;
	for (const string& seq : dataset) {
		int i = 0;
		int min_hamming_dist = numeric_limits<int>::max();
		while (i + l - 1 < (int) seq.size()) {
			string m = seq.substr(i, l);
			int hamming_dist = GetHammingDistance(m, motif);
			min_hamming_dist = min(min_hamming_dist, hamming_dist);
			i++;
		}
		score += min_hamming_dist; 
	}
	return score;
}

// Get consensus score
int GetConsensusScore(const vector<string>& dataset, int l, vector<int> pos) {
	vector<vector<int>> pfm(l, vector<int>(4));
	int score = 0;
	for (int i = 0; i < (int) dataset.size(); i++) {
		for (int j = 0; j < l; j++) {
			pfm[j][Encode(dataset[i][pos[i] + j])]++;
		}
	}
	for (int i = 0; i < l; i++) {
		score += *max_element(pfm[i].begin(), pfm[i].end());
	}
	return score;
}

// Get consensus from a given msa and Hamming score of it 
Motif FindConsensus(const vector<string>& dataset, int l, vector<int> pos) {
	vector<string> msa;
	vector<vector<int>> pfm(l, vector<int>(4));
	for (int i = 0; i < (int) dataset.size(); i++) {
		msa.emplace_back(dataset[i].substr(pos[i], l));
		for (int j = 0; j < l; j++) {
			pfm[j][Encode(msa[i][j])]++;
		}
	}
	string consensus;
	for (int i = 0; i < l; i++) {
		int mx = *max_element(pfm[i].begin(), pfm[i].end());
		for (int j = 0; j < 4; j++) {
			if (pfm[i][j] == mx) {
				consensus += alphabet[j];
				break;
			}
		}
	}
	int score = 0;
	for (const string& instance : msa) {
		int dist = GetHammingDistance(instance, consensus);
		score += dist;
	}
	return make_pair(score, consensus);
} 

// Get IC score of a given msa
double GetICScore(const vector<string>& dataset, const vector<double>& eta_vertex, int l, const vector<int>& pos) {
	vector<vector<int>> pfm(l, vector<int>(4));
	for (int i = 0; i < dataset.size(); i++) {
		for (int j = pos[i]; j < pos[i] + l; j++) {
			pfm[j - pos[i]][Encode(dataset[i][j])]++;
		}
	}
	double score = 0;
	for (int i = 0; i < l; i++) {
		for (int j = 0; j < 4; j++) {
			double ppm = (double) (pfm[i][j] + 1) / (dataset.size() + 4);
			double pwm = (double) ppm / eta_vertex[j];
			pwm = log2(pwm);
			pwm *= ppm;
			score += pwm;
			// score += (pfm[i][j] + 1) * pwm;
		}
	}
	return score;
}

// Get second order Markov IC score of a given set of instances positions
double GetICScore2(const vector<string>& dataset, const vector<vector<double>>& eta_edge, int l, const vector<int>& pos) {
	int N = (int) dataset.size();
	vector<vector<int>> cnt(l - 1, vector<int>(16));
	for (int i = 0; i < (int) dataset.size(); i++) {
		string line = dataset[i].substr(pos[i], l);
		for (int j = 0; j < l - 1; j++) {
			cnt[j][4 * Encode(line[j]) + Encode(line[j + 1])]++;
		}
	}

	double score = 0;
	for (int i = 0; i < l - 1; i++) {
		for (int j = 0; j < 16; j++) {
			double freq = (double) cnt[i][j] / N;
			if (freq != 0) {
				score += freq * log2(freq / eta_edge[j / 4][j % 4]);
			}
		}
	}
	return score;
}

// Get Complexity score of a given msa
double GetComplexityScore(vector<string>& dataset, int l, const vector<int>& pos) {
	mpfr_t l_fact;
	mpfr_init(l_fact);
	mpfr_fac_ui(l_fact, l, MPFR_RNDN);
	// mpfr_printf("%d! = %.0RF\n", l, l_fact);

	vector<vector<int>> pfm(l, vector<int>(4));
	for (int i = 0; i < dataset.size(); i++) {
		for (int j = pos[i]; j < pos[i] + l; j++) {
			pfm[j - pos[i]][Encode(dataset[i][j])]++;
		}
	}

	mpfr_t log4;
	mpfr_init(log4);
	mpfr_log_ui(log4, 4, MPFR_RNDN);
	// mpfr_printf("ln(4) = %.5RF\n", log4);
	
	double cs = 0;
	for (int i = 0; i < l; i++) {
		mpfr_t denom;
		mpfr_init(denom);
		mpfr_set_d(denom, 1, MPFR_RNDN);
		// mpfr_printf("Column: ");
		for (int j = 0; j < 4; j++) {
			mpfr_t fact;
			mpfr_init(fact);
			mpfr_fac_ui(fact, pfm[i][j], MPFR_RNDN);
			mpfr_mul(denom, denom, fact, MPFR_RNDN);
			// mpfr_printf("%d! ", pfm[i][j]);
		}
		// mpfr_printf("denom = %.0RF ", denom);
		mpfr_t result;
		mpfr_init(result);
		mpfr_div(result, l_fact, denom, MPFR_RNDN);
		// mpfr_printf("l!/denom = %.10RF ", result);
		mpfr_log(result, result, MPFR_RNDN);
		// mpfr_printf("ln(l!/denom) = %.10RF ", result);
		mpfr_div(result, result, log4, MPFR_RNDN);
		// mpfr_printf("log4(..) = %.10RF ", result);
		double db = mpfr_get_d(result, MPFR_RNDN);
		// mpfr_printf("db = %.10f\n", db);
		cs += db;
	}

	return cs;
}

// Get fitness score
double GetFitness(vector<string>& dataset, const vector<double>& eta_vertex, int l, const vector<int>& pos) {
	double v = 0.8;
	double IC = GetICScore(dataset, eta_vertex, l, pos);
	double CS = GetComplexityScore(dataset, l, pos);
	return v * IC + (1 - v) * CS;
}

// Update pheromone trail for vertices in the first column
void UpdatePheromoneVertex(vector<double>& pheromone_vertex, const string& motif, double rho) {
	int l = (int) motif.size();
	double tau_max = 1.0;
	double tau_min = 1.0 / pow(4, l);

	for (int i = 0; i < 4; i++) {
		if (motif[0] == alphabet[i]) {
			pheromone_vertex[i] = (1 - rho) * pheromone_vertex[i] + rho * tau_max;
		} else {
			pheromone_vertex[i] = (1 - rho) * pheromone_vertex[i] + rho * tau_min;
		}
	}
}

// Update pheromone trail for edges
void UpdatePheromoneEdge(vector<vector<double>>& pheromone_edge, const string& motif, double rho) {
	int l = (int) motif.size();
	double tau_max = 1;
	double tau_min = 1.0 / pow(4, l);
	
	for (int i = 0; i < l - 1; i++) {
		int u = Encode(motif[i]);
		int v = Encode(motif[i + 1]);
		int edge_index = 4 * u + v;
		for (int j = 0; j < 16; j++) {
			pheromone_edge[i][j] = (1 - rho) * pheromone_edge[i][j] + rho * (j == edge_index ? tau_max : tau_min);
		}
	}
}

// Calculate frequency of each nucleotide and nucleotide couple in dataset data set
pair<vector<double>, vector<vector<double>>> InitializeHeuristicInformation(const vector<string>& dataset) {
	int N = (int) dataset.size();
	int total_nucleotide = 0;
	int total_nucleotide_pair = 0;
	vector<int> eta_vertex(4);
	vector<vector<int>> eta_edge(4, vector<int>(4));

	for (string s : dataset) {
		int len = (int) s.size();
		total_nucleotide += len;
		total_nucleotide_pair += len - 1;

		for (char c : s) {
			eta_vertex[Encode(c)]++;
		}

		for (int i = 0; i < len - 1; i++) {
			eta_edge[Encode(s[i])][Encode(s[i + 1])]++;
		}
	}

	vector<double> r_eta_vertex(4);
	vector<vector<double>> r_eta_edge(4, vector<double>(4));
	for (int i = 0; i < 4; i++) {
		r_eta_vertex[i] = (double) eta_vertex[i] / total_nucleotide;
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			r_eta_edge[i][j] = (double) eta_edge[i][j] / total_nucleotide_pair;
		}
	}

	return make_pair(r_eta_vertex, r_eta_edge);
}

// Construct a solution for each ant
string ConstructSolution(const vector<double>& eta_vertex, const vector<vector<double>>& eta_edge, const vector<double>& pheromone_vertex, const vector<vector<double>>& pheromone_edge, int l, double alpha, double beta) {
	string motif;

	// Pick the first vertex
	vector<double> pvertex(4);
	double sum_pvertex = 0;
	for (int i = 0; i < 4; i++) {
		pvertex[i] = pow(pheromone_vertex[i], alpha) * pow(eta_vertex[i], beta);
		sum_pvertex += pvertex[i];
	}
	for (int i = 0; i < 4; i++) {
		pvertex[i] /= sum_pvertex;
	}

	double random = (double) rand() / RAND_MAX;
	sum_pvertex = 0;
	bool has_picked = false;
	for (int i = 0; i < 4; i++) {
		if (sum_pvertex <= random && random <= sum_pvertex + pvertex[i]) {
			motif += alphabet[i];
			has_picked = true;
			break;
		}
		sum_pvertex += pvertex[i];
	}
	if (!has_picked) {
		motif += 'T';
	}

	// Pick the next vertices
	for (int i = 1; i < l; i++) {
		vector<double> pedge(4);
		double sum_pedge = 0;
		for (int j = 0; j < 4; j++) {
			pedge[j] = pow(pheromone_edge[i - 1][4 * Encode(motif[i - 1]) + j], alpha) * 
					   pow(eta_edge[Encode(motif[i - 1])][j], beta);
			sum_pedge += pedge[j];
		}
		for (int j = 0; j < 4; j++) {
			pedge[j] /= sum_pedge;
		}

		random = (double) rand() / RAND_MAX;
		sum_pedge = 0;
		has_picked = false;
		for (int j = 0; j < 4; j++) {
			if (sum_pedge <= random && random <= sum_pedge + pedge[j]) {
				motif += alphabet[j];
				has_picked = true;
				break;
			}
			sum_pedge += pedge[j];
		}
		if (!has_picked) {
			motif += 'T';
		}
	}

	if ((int) motif.size() != l) {
		cout << "Error in constructing solution!\n" << l << ' ' << motif.size() << '\n';
	}	
	assert((int) motif.size() == l);

	return motif;
}

// Find the best solutions among the ant colony solutions that have score not greater than allowed_excess from the minimum
vector<Motif> FilterBestMotifs(const vector<Motif>& colony_motifs, int allowed_excess) {
	vector<Motif> best_colony_motifs;
	int best_score = numeric_limits<int>::max();

	for (const auto& p : colony_motifs) {
		int score = p.first;
		string motif = p.second;

		if (score < best_score) {
			best_score = score;
			best_colony_motifs.clear();
			best_colony_motifs.emplace_back(score, motif);
		} else if (score == best_score) {
			bool is_duplicate = false;
			for (const auto& s : best_colony_motifs) {
				if (s.second == motif) {
					is_duplicate = true;
					break;
				}
			}
			if (!is_duplicate) {
				best_colony_motifs.emplace_back(score, motif);
			}
		}
	}

	for (int excess = 1; excess <= allowed_excess; excess++) {
		vector<Motif> added_candidates;
		for (const auto& p : colony_motifs) {
			int score = p.first;
			string motif = p.second;
			if (score == best_score + excess) {
				bool is_duplicate = false;
				for (const auto & s : added_candidates) {
					if (s.second == motif) {
						is_duplicate = true;
						break;
					}
				}
				if (!is_duplicate) {
					added_candidates.emplace_back(score, motif);
				}
			}
			for (const auto& s : added_candidates) {
				best_colony_motifs.emplace_back(s);
			}
		}
	}
	return best_colony_motifs;
}

// Find the best solution among the ant colony solutions that have scores not exceed allowed_excess from the minimum, with positions
vector<pair<Motif, vector<vector<int>>>> FilterBestMotifs(const vector<pair<Motif, vector<vector<int>>>>& colony_motifs, int allowed_excess) {
	vector<pair<Motif, vector<vector<int>>>> best_colony_motifs;
	int best_score = numeric_limits<int>::max();

	for (const auto& p : colony_motifs) {
		int score = p.first.first;
		string motif = p.first.second;

		if (score < best_score) {
			best_score = score;
			best_colony_motifs.clear();
			best_colony_motifs.emplace_back(make_pair(score, motif), p.second);
		} else if (score == best_score) {
			bool is_duplicate = false;
			for (const auto& s : best_colony_motifs) {
				if (s.first.second == motif) {
					is_duplicate = true;
					break;
				}
			}
			if (!is_duplicate) {
				best_colony_motifs.emplace_back(make_pair(score, motif), p.second);
			}
		}
	}

	for (int excess = 1; excess <= allowed_excess; excess++) {
		vector<pair<Motif, vector<vector<int>>>> added_candidates;
		for (const auto& p : colony_motifs) {
			int score = p.first.first;
			string motif = p.first.second;
			if (score == best_score + excess) {
				bool is_duplicate = false;
				for (const auto & s : added_candidates) {
					if (s.first.second == motif) {
						is_duplicate = true;
						break;
					}
				}
				if (!is_duplicate) {
					added_candidates.emplace_back(make_pair(score, motif), p.second);
				}
			}
			for (const auto& s : added_candidates) {
				best_colony_motifs.emplace_back(s);
			}
		}
	}
	return best_colony_motifs;
}

// Find msa for a given motif (multiple positions on one line allowed)
vector<vector<int>> FindMSA(const vector<string>& dataset, const string& motif) {
	int l = (int) motif.size();
	vector<vector<int>> pos;
	int score = 0;
	for (const auto& seq : dataset) {
		int i = 0;
		pos.emplace_back(vector<int>());
		int min_hamming_dist = numeric_limits<int>::max();
		while (i + l - 1 < (int) seq.size()) {
			string m = seq.substr(i, l);
			int hamming_dist = GetHammingDistance(m, motif);
			if (hamming_dist < min_hamming_dist) {
				min_hamming_dist = hamming_dist;
				pos.back().clear();
				pos.back().emplace_back(i);
			} else if (hamming_dist == min_hamming_dist) {
				pos.back().emplace_back(i);
			}
			i++;
		}
		score += min_hamming_dist;
	}
	return pos;
}

// --------------------- ACO Motif -----------------------------

// Find instances in the dataset data set with distance not l/3 greater than the minimum
vector<vector<string>> FindInstances(const vector<string>& dataset, const Motif& motif) {
	int l = (int) motif.second.size();
	vector<vector<string>> instances(dataset.size());
	int allowed_excess = l / 3; // NEW APPROACH - NEED REVISION!

	for (int line_index = 0; line_index < (int) dataset.size(); line_index++) {
		const auto& line = dataset[line_index];
		int i = 0;
		int min_hamming_dist = numeric_limits<int>::max();

		while (i + l - 1 < (int) line.size()) {
			string m = line.substr(i, l);
			int hamming_dist = GetHammingDistance(m, motif.second);
			if (hamming_dist < min_hamming_dist) {
				min_hamming_dist = hamming_dist;
				instances[line_index].clear();
				instances[line_index].emplace_back(m);
			} else if (hamming_dist == min_hamming_dist) {
				instances[line_index].emplace_back(m);
			}
			i++;
		}

		for (int excess = 1; excess <= allowed_excess; excess++) {
			i = 0;
			while (i + l - 1 < (int) line.size()) {
				string m = line.substr(i, l);
				int hamming_dist = GetHammingDistance(m, motif.second);
				if (hamming_dist == min_hamming_dist + excess) {
					instances[line_index].emplace_back(m);
				}
				i++;
			}
		}
	}

	return instances;
}

// Get score for a given motif with given list of instances (used in local search)
int GetScoreFromInstances(const vector<vector<string>>& instances, const string& motif) {
	int l = (int) motif.size();
	int score = 0;

	for (int line_index = 0; line_index < (int) instances.size(); line_index++) {
		int min_hamming_dist = numeric_limits<int>::max();
		for (const auto& instance : instances[line_index]) {
			int hamming_dist = GetHammingDistance(instance, motif);
			min_hamming_dist = min(min_hamming_dist, hamming_dist);
		}
		score += min_hamming_dist;
	}

	return score;
}

// Apply local search on a given motif
pair<Motif, vector<Motif>> LocalSearchACO(const vector<string>& dataset, const Motif& original_motif) {
	int l = (int) original_motif.second.size();
	Motif prev_motif(numeric_limits<int>::max(), "");
	Motif cur_motif = original_motif;
	vector<Motif> best_motifs_ls = {original_motif};

	auto instances = FindInstances(dataset, original_motif); // Find "special" instances - NEED REVISION!

	while (cur_motif.first < prev_motif.first && cur_motif.second != prev_motif.second) { // While the objective function can still be improved
		prev_motif = cur_motif;
		string motif = cur_motif.second;

		for (int i = 0; i < l; i++) {
			vector<Motif> modified_motifs;
			for (char c : alphabet) {
				if (motif[i] != c) {
					string modified_motif = motif;
					modified_motif[i] = c;
					modified_motifs.emplace_back(GetScoreFromInstances(instances, modified_motif), modified_motif);
				}
			}

			sort(modified_motifs.begin(), modified_motifs.end());

			for (const auto& p : modified_motifs) {
				best_motifs_ls.emplace_back(p);
			}

			if (modified_motifs[0].first < cur_motif.first) {
				cur_motif = modified_motifs[0];
				best_motifs_ls.emplace_back(modified_motifs[0]);
			}
		}
	}

	return make_pair(cur_motif, best_motifs_ls);
}

// Find minimum distance on each line in dataset data set
vector<int> FindDistances(const vector<string>& dataset, const string& motif) {
	int l = (int) motif.size();
	vector<int> distances(dataset.size(), numeric_limits<int>::max());
	for (int line_index = 0; line_index < (int) dataset.size(); line_index++) {
		const auto& line = dataset[line_index];
		int i = 0;
		while (i + l - 1 < (int) line.size()) {
			string m = line.substr(i, l);
			int hamming_dist = GetHammingDistance(m, motif);
			if (hamming_dist < distances[line_index]) {
				distances[line_index] = hamming_dist;
			}
			i++;
		}
	}
	return distances;
}

// ----------------------- Relax-ACO Motif -----------------------------

// Apply local search on a given motif
pair<vector<Motif>, Motif> LocalSearchRACO(const vector<string>& dataset, Motif motif, const vector<double>& eta_vertex) {
	int l = motif.second.size();
	vector<Motif> motifs_ls;

	int current_best_score = motif.first;
	string current_best_motif = motif.second;

	int past_best_score = 0;
	string past_best_motif = current_best_motif;

	string original_motif = motif.second;

	bool has_started = false;

	while (!has_started || (current_best_score > past_best_score && current_best_motif != past_best_motif)) {
		has_started = true;
		past_best_score = current_best_score;
		past_best_motif = current_best_motif;
		
		for (int i = 0; i < l; i++) {
			vector<Motif> generated_motifs;

			for (char c : alphabet) {
				if (original_motif[i] != c) {
					char t = original_motif[i];
					original_motif[i] = c;
					auto msa = FindMSA(dataset, original_motif);	
					vector<int> pos(msa.size());
					for (int i = 0; i < msa.size(); i++) {
						pos[i] = msa[i][rand() % msa[i].size()];
					}
					generated_motifs.emplace_back(GetICScore(dataset, eta_vertex, l, pos), original_motif);
					original_motif[i] = t;
				}
			}

			sort(generated_motifs.begin(), generated_motifs.end(), greater<Motif>());
			for (const auto& p : generated_motifs) {
				motifs_ls.emplace_back(p);
			}
			
			assert(!generated_motifs.empty());
			auto best_motif = generated_motifs[0];
			if (best_motif.first > current_best_score) {
				current_best_score = best_motif.first;
				current_best_motif = best_motif.second;
				original_motif = current_best_motif;
			}
		}
	}

	auto best_motif_ls = make_pair(current_best_score, current_best_motif);
	assert(original_motif == current_best_motif);
	return make_pair(motifs_ls, best_motif_ls);
}