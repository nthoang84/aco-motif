#include <bits/stdc++.h>

using namespace std;

const string alphabet = "ACGT";

vector<string> input;
string motif;
vector<vector<int>> position;
vector<pair<int, vector<int>>> results;

// Map char values to int values
inline int Encode(char c) {
	return (c == 'A' ? 0 : (c == 'C' ? 1 : (c == 'G' ? 2 : 3)));
}

// Get data from arguments
void GetData(char* input_dir, char* con) {
	ifstream is;
	is.open(input_dir);
	string str;
	while (getline(is, str)) {
		if (str[0] == '>') {
			input.emplace_back("");
		} else if (!str.empty()) {
			for (char& c : str) {
				if ('a' <= c && c <= 'z') {
					c = c - 'a' + 'A';
				}
			}
			input.back() += str;
		}
	}
	is.close();
	string tmp2(con);
	for (int i = 0; i < tmp2.size(); ++i) {
		motif.push_back(tmp2[i]);
	}
}

// Get Hamming distance of two strings with equal length
int GetHammingDistance(string a, int starta, string b, int startb, int len) {
	int dist = 0;
	for (int i = 0; i < len; ++i){
		if (a[starta + i] != b[startb + i]) {
			++dist;
		}
	}
	return dist;
}

// Get Consensus score of a given msa
int GetConsensusScore(vector<int> start, int l, int N) {
	vector<vector<int>> freq(l, vector<int>(4));
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < l; ++j) {
			freq[j][Encode(input[i][j + start[i]])]++;
		}
	}
	int score = 0;
	for (int i = 0; i < l; ++i) {
		score += *max_element(freq[i].begin(), freq[i].end());
	}
	return score;
}

// Get consensus motif of a given msa
string GetConsensus(vector<int> start, int l, int N) {
	vector<vector<int>> freq(l, vector<int>(4));
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < l; ++j) {
			freq[j][Encode(input[i][j + start[i]])]++;
		}
	}
	string consensus;
	for (int i = 0; i < l; ++i) {
		int mx = *max_element(freq[i].begin(), freq[i].end());
		for (int j = 0; j < 4; j++) {
			if (freq[i][j] == mx) {
				consensus += alphabet[j];
				break;
			}
		}
	}
	return consensus;
}

// Find msa for a given motif (multiple positions on one line allowed)
void FindMSA(int allowed_excess) {
	int l = motif.size();
	for (int i = 0; i < input.size(); ++i) {
		position.push_back(vector<int>());
		int lastPos = input[i].size() - l;
		int score = INT_MAX;
		for (int j = 0; j <= lastPos; ++j) {
			int hamming = GetHammingDistance(motif, 0, input[i], j, motif.size());
			if (hamming < score) {
				score = hamming;
				vector<int> v;
				for (int k = 0; k < position[i].size(); ++k) {
					if (GetHammingDistance(motif, 0, input[i], position[i][k], motif.size()) <= score + allowed_excess) {
						v.push_back(position[i][k]);
					}
				}
				position[i].clear();
				position[i] = v;
				position[i].push_back(j);
				score = hamming;
			} else if (hamming <= score + allowed_excess) {
				position[i].push_back(j);
			}
		}
	}
}

// Get Hamming score of a given motif
int GetHammingScore(const string& motif, int N) {
	int l = (int) motif.size();
	int score = 0;
	for (int i = 0; i < N; i++) {
		int min_dist = numeric_limits<int>::max();
		for (int j = 0; j + l - 1 < (int) input[i].size(); j++) {
			int dist = GetHammingDistance(motif, 0, input[i], j, l);
			min_dist = min(min_dist, dist);
		}
		score += min_dist;
	}
	return score;
}

// Get Hamming Score of the consensus motif built from a given msa
int GetHammingScoreMSA(const vector<int>& pos) {
	int l = motif.size();
	int N = input.size();
    vector<vector<int>> cnt(l, vector<int>(4));
	for (int i = 0; i < N; ++i) {
		int p = position[i][pos[i]];
		for (int j = p; j < p + l; ++j) {
            cnt[j - p][Encode(input[i][j])]++;
		}
	}
	int score = 0;
	for (int i = 0; i < l; i++) {
		score += N - *max_element(cnt[i].begin(), cnt[i].end());
	}
	return score;
}

// Insert a solution to the solutions set
void Insert(int results_size, int score, const vector<int>& pos) {
	assert(score >= 0);
	int N = input.size();
	int i = results_size - 1;

	while (i > 0 && results[i].first > score) {
		results[i] = results[i - 1];
		i--;
	}

	if (i == 0) {
		if (score < results[0].first) {
			results[0].first = score;
			for (int j = 0; j < N; j++) {
				results[0].second[j] = position[j][pos[j]];
			}
		} else if (results_size > 1) {
			results[0].first = score;
			for (int j = 0; j < N; j++) {
				results[0].second[j] = position[j][pos[j]];
			}
		}
	} else if (i < results_size - 1) {
		results[i + 1].first = score;
		for (int j = 0; j < N; j++) {
			results[i + 1].second[j] = position[j][pos[j]];
		}
	}
}

int main(int argc, char* argv[]) {
	GetData(argv[1], argv[2]);

	int N = input.size();
	int l = motif.size();

    FindMSA(atoi(argv[3]));
	int init_sum = 0;
    for(int i = 0; i < position.size(); ++i) {
		init_sum += position[i].size();
    }

	auto input_init = input;
	double epsilon = 0.5;
	
	while (true) {
		if (init_sum < 2 * position.size()) {
			epsilon = 0;
		}

		vector<int> index(N);
		iota(index.begin(), index.end(), 0);

		for (int i = 0; i < N; ++i) {
			for (int j = i + 1; j < N; ++j) {
				if (position[i].size() > position[j].size()) {
					swap(position[i], position[j]);
					swap(index[i], index[j]);
				}
			}
		}

		int iter = 0;
		for (int i = 0; i < N; ++i) {
			input[i] = input_init[index[i]];
			if (position[i].size() == 1) {
				iter = i;
			}
		}

		for (int i = iter + 1; i < N; ++i) {
			vector<int> h(position[i].size());
			int g = numeric_limits<int>::max();
			vector<int> selected_pos;

			for (int k = 0; k < position[i].size(); ++k) {
				h[k] = 0;
				for (int j = 0; j < i; ++j) {
					for (int t = 0; t < (int) position[j].size(); t++) {
						if (position[j].size() == 1) {
							h[k] += GetHammingDistance(input[j], position[j][t], input[i], position[i][k], l);
						}
					}
				}
				double amount = epsilon;
				if (epsilon < 1) amount = epsilon * g;

				if (h[k] < g + amount) {
					if (h[k] < g){
						g = h[k];
						int sz = selected_pos.size();
						selected_pos.clear();
						for (int m = 0; m < sz; ++m) {
							if (h[m] < g + amount) {
								selected_pos.push_back(m);	
							}
						}
					}
					selected_pos.push_back(k);
				}
			}

			for (int m = 0; m < selected_pos.size(); ++m) {
				selected_pos[m] = position[i][selected_pos[m]];
			}
			position[i].clear();
			position[i] = selected_pos;
		}
	
		epsilon *= 0.98;

		auto pos_temp = position;
		auto input_temp = input;
		for (int i = 0; i < N ; ++i) {
			position[index[i]] = pos_temp[i];
			input[index[i]] = input_temp[i];
		}

		int cur_sum = 0;
		for (int i = 0; i < N; ++i) {
			cur_sum += position[i].size();
	    }

		if (cur_sum < 40) {
			if (cur_sum * 2 < init_sum || cur_sum == N) {
				break;
			}
		}
	}

	assert(input == input_init);

	int results_size = atoi(argv[4]);
	results.assign(results_size, make_pair(numeric_limits<int>::max(), vector<int>(N)));

	vector<int> pos(N);
	int score = GetHammingScoreMSA(pos);
	Insert(results_size, score, pos);

	int iter = N - 1;
	while (true) {
		while (iter >= 0 && pos[iter] >= (position[iter].size() - 1)) iter--;
		if (iter < 0) break;
		pos[iter]++;
		for (int i = iter + 1; i < N; ++i) {
			pos[i] = 0;
		}
		iter = N - 1;	
		score = GetHammingScoreMSA(pos);
		Insert(results_size, score, pos);
	}

	cout << "Positions: ";
	vector<int> result_msa;
	for (int i = 0; i < results_size; ++i) {
		for (int j = 0; j < N; ++j) {
			result_msa.push_back(results[i].second[j]);
			cout << result_msa.back() + 1 << ' ';
		}
		cout << '\n';
	}

	cout << "Motif: " << GetConsensus(result_msa, l, N) << '\n';
    int csc = GetConsensusScore(result_msa, l, N);
    cout << "CSc: " << csc << '\n';
    cout << "HSc: "<< l * N - csc << '\n';
    cout << "Origin: " << motif << '\n';
    cout << "------------------------------\n";

	return 0;
}
