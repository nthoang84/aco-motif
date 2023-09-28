#include <mpfr.h>
#include <bits/stdc++.h>

using namespace std;
using motif_t = pair<double, string>;

const string alphabet = "ACGT";

inline int encode(char c) {
	return (c == 'A' ? 0 : (c == 'C' ? 1 : (c == 'G' ? 2 : 3)));
}

void read_file(vector<string>& dataset, vector<string>& sequence_names, const string& dataset_name) {
	string dataset_path = "./data/" + dataset_name + "/dataset.fa";
	ifstream is;
	is.open(dataset_path);
	string sequence;
	while (getline(is, sequence)) {
		if (sequence[0] == '>') {
			dataset.emplace_back("");
			sequence_names.emplace_back("");

			istringstream iss(sequence);
			iss >> sequence_names.back();
			assert(sequence_names.back()[0] == '>');
			sequence_names.back() = sequence_names.back().substr(1);
		
		} else if (!sequence.empty()) {
			for (char& c : sequence) {
				if ('a' <= c && c <= 'z') {
					c = c - 'a' + 'A';
				}
			}
			dataset.back() += sequence;
		}
	}
	is.close();
}

void filter_dust(vector<string>& dataset, const vector<string>& sequence_names, const string& dataset_name) {
	string dust_file_path = "./data/" + dataset_name + "/result.dust";
	ifstream is;
	is.open(dust_file_path);
	string line;
	while (getline(is, line)) {
		istringstream iss(line);
		string name;
		int start, finish;
		iss >> name >> start >> finish;
		for (int i = 0; i < (int) dataset.size(); i++) {
			if (sequence_names[i] == name) {
				for (int j = start; j <= finish; j++) {
					dataset[i][j] = 'N';
				}
				break;
			}
		}
	}
	is.close();
}

int get_hamming_distance(const string& s, const string& t) {
	assert(s.size() == t.size());
	int distance = 0;
	for (int i = 0; i < (int) s.size(); i++) {
		distance += (s[i] != t[i]);
	}
	return distance;
}

int get_hamming_score(const vector<string>& dataset, const string& motif) {
	int l = motif.size();
	int score = 0;
	for (const string& sequence : dataset) {
		int min_distance = numeric_limits<int>::max();
		for (int i = 0; i + l - 1 < sequence.size(); i++) {
			string s = sequence.substr(i, l);
			int distance = get_hamming_distance(s, motif);
			min_distance = min(min_distance, distance);
		}
		score += min_distance; 
	}
	return score;
}

int get_consensus_score(const vector<string>& dataset, int l, const vector<int> positions) {
	vector<vector<int>> frequency(l, vector<int>(4));
	int score = 0;
	for (int i = 0; i < (int) dataset.size(); i++) {
		for (int j = 0; j < l; j++) {
			frequency[j][encode(dataset[i][positions[i] + j])]++;
		}
	}
	for (int i = 0; i < l; i++) {
		score += *max_element(frequency[i].begin(), frequency[i].end());
	}
	return score;
}

motif_t get_consensus_string(const vector<string>& dataset, int l, const vector<int>& positions) {
	vector<string> msa;
	vector<vector<int>> frequency(l, vector<int>(4));
	for (int i = 0; i < (int) dataset.size(); i++) {
		msa.emplace_back(dataset[i].substr(positions[i], l));
		for (int j = 0; j < l; j++) {
			frequency[j][encode(msa[i][j])]++;
		}
	}
	string consensus_string;
	for (int i = 0; i < l; i++) {
		int max_frequency = *max_element(frequency[i].begin(), frequency[i].end());
		for (int j = 0; j < 4; j++) {
			if (frequency[i][j] == max_frequency) {
				consensus_string += alphabet[j];
				break;
			}
		}
	}
	int score = 0;
	for (const string& alignment : msa) {
		int distance = get_hamming_distance(alignment, consensus_string);
		score += distance;
	}
	return make_pair(score, consensus_string);
} 

double get_information_content_score(const vector<string>& dataset, const vector<double>& eta_vertex, int l, const vector<int>& positions) {
	vector<vector<int>> frequency(l, vector<int>(4));
	for (int i = 0; i < dataset.size(); i++) {
		for (int j = positions[i]; j < positions[i] + l; j++) {
			frequency[j - positions[i]][encode(dataset[i][j])]++;
		}
	}
	double score = 0;
	for (int i = 0; i < l; i++) {
		for (int j = 0; j < 4; j++) {
			double relative_frequency = (double) (frequency[i][j] + 1) / (dataset.size() + 4);
			score += relative_frequency * log2(relative_frequency / eta_vertex[j]);
		}
	}
	return score;
}

double get_information_content_score_2(const vector<string>& dataset, const vector<vector<double>>& eta_edge, int l, const vector<int>& positions) {
	vector<vector<int>> frequency(l - 1, vector<int>(16));
	for (int i = 0; i < (int) dataset.size(); i++) {
		string subsequence = dataset[i].substr(positions[i], l);
		for (int j = 0; j < l - 1; j++) {
			frequency[j][4 * encode(subsequence[j]) + encode(subsequence[j + 1])]++;
		}
	}
	double score = 0;
	for (int i = 0; i < l - 1; i++) {
		for (int j = 0; j < 16; j++) {
			double relative_frequency = (double) (frequency[i][j] + 1) / (dataset.size() + 16);
			score += relative_frequency * log2(relative_frequency / eta_edge[j / 4][j % 4]);
		}
	}
	return score;
}

double get_complexity_score(const vector<string>& dataset, int l, const vector<int>& positions) {
	mpfr_t l_factorial;
	mpfr_init(l_factorial);
	mpfr_fac_ui(l_factorial, l, MPFR_RNDN);

	vector<vector<int>> frequency(l, vector<int>(4));
	for (int i = 0; i < dataset.size(); i++) {
		for (int j = positions[i]; j < positions[i] + l; j++) {
			frequency[j - positions[i]][encode(dataset[i][j])]++;
		}
	}

	mpfr_t ln4;
	mpfr_init(ln4);
	mpfr_log_ui(ln4, 4, MPFR_RNDN);
	
	double score = 0;
	for (int i = 0; i < l; i++) {
		mpfr_t product;
		mpfr_init(product);
		mpfr_set_d(product, 1, MPFR_RNDN);
		for (int j = 0; j < 4; j++) {
			mpfr_t factorial;
			mpfr_init(factorial);
			mpfr_fac_ui(factorial, frequency[i][j], MPFR_RNDN);
			mpfr_mul(product, product, factorial, MPFR_RNDN);
		}
		mpfr_t added_term;
		mpfr_init(added_term);
		mpfr_div(added_term, l_factorial, product, MPFR_RNDN);
		mpfr_log(added_term, added_term, MPFR_RNDN);
		mpfr_div(added_term, added_term, ln4, MPFR_RNDN);
		score += mpfr_get_d(added_term, MPFR_RNDN);
	}

	return score;
}

double get_fitness_score(const vector<string>& dataset, const vector<double>& eta_vertex, int l, const vector<int>& positions) {
	double v = 0.8;
	double IC = get_information_content_score(dataset, eta_vertex, l, positions);
	double CS = get_complexity_score(dataset, l, positions);
	return v * IC + (1 - v) * CS;
}

double get_pearson_correlation(const vector<vector<int>>& va, const vector<vector<int>>& vb) {
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
	double numerator = 0;
	for (int i = 0; i < n; i++) {
		numerator += (a[i] - avg_a) * (b[i] - avg_b);
	}
	double denominator = 1;
	for (int t = 0; t < 2; t++) {
		const auto& v = (t == 0 ? a : b);
		const double& avg_v = (t == 0 ? avg_a : avg_b);
		double sum = 0;
		for (int i = 0; i < n; i++) {
			sum += pow(v[i] - avg_v, 2);
		}
		denominator *= sqrt(sum);
	}
	return numerator / denominator;
}

void update_pheromone_vertex(vector<double>& pheromone_vertex, const string& motif, double rho) {
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

void update_pheromone_edge(vector<vector<double>>& pheromone_edge, const string& motif, double rho) {
	int l = (int) motif.size();
	double tau_max = 1;
	double tau_min = 1.0 / pow(4, l);
	
	for (int i = 0; i < l - 1; i++) {
		int u = encode(motif[i]);
		int v = encode(motif[i + 1]);
		int edge_index = 4 * u + v;
		for (int j = 0; j < 16; j++) {
			pheromone_edge[i][j] = (1 - rho) * pheromone_edge[i][j] + rho * (j == edge_index ? tau_max : tau_min);
		}
	}
}

pair<vector<double>, vector<vector<double>>> get_heuristic_info(const vector<string>& dataset) {
	int total_nucleotides = 0;
	int total_nucleotides_pairs = 0;
	vector<int> eta_vertex(4);
	vector<vector<int>> eta_edge(4, vector<int>(4));

	for (const string& sequence : dataset) {
		int len = (int) sequence.size();
		total_nucleotides += len;
		total_nucleotides_pairs += len - 1;

		for (char c : sequence) {
			eta_vertex[encode(c)]++;
		}

		for (int i = 0; i < len - 1; i++) {
			eta_edge[encode(sequence[i])][encode(sequence[i + 1])]++;
		}
	}

	vector<double> r_eta_vertex(4);
	vector<vector<double>> r_eta_edge(4, vector<double>(4));
	for (int i = 0; i < 4; i++) {
		r_eta_vertex[i] = (double) eta_vertex[i] / total_nucleotides;
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			r_eta_edge[i][j] = (double) eta_edge[i][j] / total_nucleotides_pairs;
		}
	}

	return make_pair(r_eta_vertex, r_eta_edge);
}

string construct_solution(const vector<double>& eta_vertex, const vector<vector<double>>& eta_edge, const vector<double>& pheromone_vertex, const vector<vector<double>>& pheromone_edge, int l, double alpha, double beta) {
	string solution;
	double random;
	bool has_picked;

	vector<double> pvertex(4);
	double sum_pvertex = 0;
	for (int i = 0; i < 4; i++) {
		pvertex[i] = pow(pheromone_vertex[i], alpha) * pow(eta_vertex[i], beta);
		sum_pvertex += pvertex[i];
	}
	for (int i = 0; i < 4; i++) {
		pvertex[i] /= sum_pvertex;
	}

	random = (double) rand() / RAND_MAX;
	sum_pvertex = 0;
	has_picked = false;
	for (int i = 0; i < 4; i++) {
		if (sum_pvertex <= random && random <= sum_pvertex + pvertex[i]) {
			solution += alphabet[i];
			has_picked = true;
			break;
		}
		sum_pvertex += pvertex[i];
	}
	if (!has_picked) {
		solution += 'T';
	}

	for (int i = 1; i < l; i++) {
		vector<double> pedge(4);
		double sum_pedge = 0;
		for (int j = 0; j < 4; j++) {
			pedge[j] = pow(pheromone_edge[i - 1][4 * encode(solution[i - 1]) + j], alpha) * 
					   pow(eta_edge[encode(solution[i - 1])][j], beta);
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
				solution += alphabet[j];
				has_picked = true;
				break;
			}
			sum_pedge += pedge[j];
		}
		if (!has_picked) {
			solution += 'T';
		}
	}

	if ((int) solution.size() != l) {
		cout << "Error in constructing solution!\n" << l << ' ' << solution.size() << '\n';
	}	
	assert((int) solution.size() == l);

	return solution;
}

vector<pair<motif_t, vector<vector<int>>>> get_best_solutions_by_hamming_score(const vector<pair<motif_t, vector<vector<int>>>>& colony_motifs, int allowed_excess) {
	vector<pair<motif_t, vector<vector<int>>>> best_solutions;
	int best_score = numeric_limits<int>::max();

	for (const auto& p : colony_motifs) {
		int score = p.first.first;
		string motif = p.first.second;

		if (score < best_score) {
			best_score = score;
			best_solutions.clear();
			best_solutions.emplace_back(make_pair(score, motif), p.second);

		} else if (score == best_score) {
			bool is_duplicate = false;
			for (const auto& s : best_solutions) {
				if (s.first.second == motif) {
					is_duplicate = true;
					break;
				}
			}
			if (!is_duplicate) {
				best_solutions.emplace_back(make_pair(score, motif), p.second);
			}
		}
	}

	for (int excess = 1; excess <= allowed_excess; excess++) {
		vector<pair<motif_t, vector<vector<int>>>> added_candidates;
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
		}
		for (const auto& s : added_candidates) {
			best_solutions.emplace_back(s);
		}
	}
	return best_solutions;
}

vector<vector<int>> get_all_positions(const vector<string>& dataset, const string& motif) {
	int l = (int) motif.size();
	vector<vector<int>> positions;
	for (const auto& sequence : dataset) {
		positions.emplace_back(vector<int>());
		int min_distance = numeric_limits<int>::max();
		for (int i = 0; i + l - 1 < (int) sequence.size(); i++) {
			string s = sequence.substr(i, l);
			int distance = get_hamming_distance(s, motif);
			if (distance < min_distance) {
				min_distance = distance;
				positions.back().clear();
				positions.back().emplace_back(i);
			} else if (distance == min_distance) {
				positions.back().emplace_back(i);
			}
		}
	}
	return positions;
}

vector<int> get_random_positions(const vector<string>& dataset, const string& motif) {
	auto all_positions = get_all_positions(dataset, motif);
	vector<int> positions(all_positions.size());
	for (int i = 0; i < (int) positions.size(); i++) {
		positions[i] = all_positions[i][rand() % all_positions[i].size()];
	}
	return positions;
}

vector<int> get_first_positions(const vector<string>& dataset, const string& motif) {
	auto all_positions = get_all_positions(dataset, motif);
	vector<int> positions(all_positions.size());
	for (int i = 0; i < (int) positions.size(); i++) {
		positions[i] = all_positions[i][0];
	}
	return positions;
}

pair<vector<motif_t>, motif_t> local_search_by_hamming_score(const vector<string>& dataset, motif_t motif) {
	int l = motif.second.size();
	vector<motif_t> motifs_ls;

	int current_best_score = motif.first;
	string current_best_motif = motif.second;

	int past_best_score = 0;
	string past_best_motif = current_best_motif;

	string original_motif = motif.second;

	bool has_started = false;

	while (!has_started || (current_best_score < past_best_score && current_best_motif != past_best_motif)) {
		has_started = true;
		past_best_score = current_best_score;
		past_best_motif = current_best_motif;
		
		for (int i = 0; i < l; i++) {
			vector<motif_t> generated_motifs;

			for (char c : alphabet) {
				if (original_motif[i] != c) {
					char t = original_motif[i];
					original_motif[i] = c;
					generated_motifs.emplace_back(get_hamming_score(dataset, original_motif), original_motif);
					original_motif[i] = t;
				}
			}

			sort(generated_motifs.begin(), generated_motifs.end());
			for (const auto& p : generated_motifs) {
				motifs_ls.emplace_back(p);
			}
			
			assert(!generated_motifs.empty());
			auto best_motif = generated_motifs[0];
			if (best_motif.first < current_best_score) {
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