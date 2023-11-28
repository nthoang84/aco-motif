#include<bits/stdc++.h>

using namespace std;

vector<vector<char>> sequence;
vector<char> consensus;
vector<vector<int>> position;
vector<int> topHamming;
vector<vector<int>> topPosition;

void readFile(char* seq, char* con) {
	ifstream in;
	in.open(seq);
	string tmp;
	while (!in.eof()) {
		vector<char> vec;
		getline(in, tmp);
		getline(in, tmp);
		for(int i = 0; i < tmp.length() - 1; ++i) {
			vec.push_back(toupper(tmp[i]));
		}
		sequence.push_back(vec);
		getline(in, tmp);
	}
	sequence[sequence.size() - 1].push_back(toupper(tmp[tmp.size() - 1]));
	in.close(); 
	string tmp2(con);
	for (int i = 0; i < tmp2.size(); ++i) {
		consensus.push_back(tmp2[i]);
	}
}

int hammingDistance(vector<char> a, int starta, vector<char> b, int startb, int len) {
	int tmp = 0;
	for (int i = 0; i < len; ++i){
		if ( a[starta + i] != b[startb + i]) {
			++tmp;
		}
	}
	return tmp;
}

int consensusScore(vector<vector<char>> data, vector<int> start, int len, int n) { // len: motif size, n: input size
	int fre[len][4];	
	for (int i = 0; i < len; ++i) {
		for (int j = 0; j < 4; ++j) {
			fre[i][j] = 0;
		}
	}
	for (int i = 0; i < n; ++i){
		for (int j = 0; j < len; ++j) {
			switch(data[i][j + start[i]]) {
				case 'A':
					fre[j][0]++;
					break;
				case 'T':
					fre[j][1]++;
					break;
				case 'C':
					fre[j][2]++;
					break;
				case 'G':
					fre[j][3]++;
			}
		}
	}

	int score = 0;
	for (int i = 0; i < len; ++i) {
		int max = fre[i][0];
		for (int j = 1; j < 4; ++j) {
			if (fre[i][j] > max){
				max = fre[i][j];
			}
		}
		score += max;
	}
	return score;
}

void determinePosition(int thres) {
	int l = consensus.size();
	for (int i = 0; i < sequence.size(); ++i) {
		position.push_back(vector<int>());
		int lastPos = sequence[i].size() - l;
		int score = INT_MAX;
		for (int j = 0; j <= lastPos; ++j) {
			int hamming = hammingDistance(consensus, 0, sequence[i], j, consensus.size());
			if (hamming < score) {
				score = hamming;
				vector<int> v;
				for (int k = 0; k < position[i].size(); ++k) {
					if (hammingDistance(consensus, 0, sequence[i], position[i][k], consensus.size()) <= score + thres) {
						v.push_back(position[i][k]);
					}
				}
				position[i].clear();
				position[i] = v;
				position[i].push_back(j);
				score = hamming;
			} else if (hamming <= score + thres) {
				position[i].push_back(j);
			}
		}
	}
}

int consensusCalculate(int* pos) {
	int l = consensus.size();
	int select = sequence.size();

    vector<vector<int>> tmp(l, vector<int>(4));
	for (int i = 0; i < sequence.size(); ++i) {
		int a = position[i][pos[i]];
		for (int j = a; j < a + l; ++j) {
			if (sequence[i][j] == 'A') {
				++tmp[j - a][0];
			} else if (sequence[i][j] == 'T') {
                ++tmp[j - a][1];
            } else if (sequence[i][j] == 'C') {
                ++tmp[j - a][2];
            } else if (sequence[i][j] == 'G') {
                ++tmp[j - a][3];
            }
		}
	}

	int hamming = 0;
	for (int i = 0; i < l; ++i) {
		if (tmp[i][0] >= tmp[i][1]) {
			if (tmp[i][0] >= tmp[i][2]) {
				if (tmp[i][0] >= tmp[i][3]) {
					hamming += select - tmp[i][0];
				} else {
					hamming += select - tmp[i][3];
				}
			} else {
				if (tmp[i][2] >= tmp[i][3]) {
					hamming += select - tmp[i][2];
				} else {
					hamming += select - tmp[i][3];
				}
			}
		} else {
			if (tmp[i][1] >= tmp[i][2]) {
	            if (tmp[i][1] >= tmp[i][3]) {
					hamming += select - tmp[i][1];
	            } else {
					hamming += select - tmp[i][3];
	            }
        	} else {
	            if (tmp[i][2] >= tmp[i][3]){
					hamming += select - tmp[i][2];
	            } else {
					hamming += select - tmp[i][3];
	            }
        	}
		}
	}
	return hamming;
}

void top(int num, int score, int* pos) {
	assert(score >= 0);
	int i = num - 1;
	int len = sequence.size();	
	while (topHamming[i] > score && i > 0) {
		topHamming[i] = topHamming[i - 1];
		for (int j = 0; j < len; ++j) {
			topPosition[i][j] = topPosition[i - 1][j];
		}
		--i;
	}

	if (i == 0) {
		if (score < topHamming[0]) {
			topHamming[0] = score;
			for (int j = 0; j < len; ++j) {
				topPosition[0][j] = position[j][pos[j]];
			}
		} else {
			if (num > 1) {
				topHamming[1] = score;
				for (int j = 0; j < len; ++j) {	
					topPosition[1][j] = position[j][pos[j]];
				}
			}
		}

	} else if (i < num - 1){
		topHamming[i + 1] = score;
		for (int j = 0; j < len; ++j) {
			topPosition[i + 1][j] = position[j][pos[j]];
		}
	}
}


int main(int argc, char* argv[]) {
	clock_t start = clock();
	readFile(argv[1], argv[2]);

	int l = sequence.size();
	int len = consensus.size();
	int candidate[10][l];

    determinePosition(atoi(argv[3]));
	int posNum = 0;
    for(int i = 0; i < position.size(); ++i) {
		posNum += position[i].size();
    }

	auto sequenceTemp = sequence;
	double epsilon = 0.5;
	
	while (true) {
		if (posNum < 2 * position.size()) {
			epsilon = 0;
		}

		vector<int> index(l);
		for(int i = 0; i < l; ++i) {
			index[i] = i;
		}

		for(int i = 0; i < l; ++i) {
			for (int j = i + 1; j < l; ++j) {
				if (position[i].size() > position[j].size()) {
					swap(position[i], position[j]);
					swap(index[i], index[j]);
				}
			}
		}

		int firstOccurent = 0;
		for(int i = 0; i < l; ++i) {
			sequence[i] = sequenceTemp[index[i]];
			if (position[i].size() == 1) {
				firstOccurent = i;
			}
		}

		for (int i = firstOccurent + 1; i < l; ++i) {
			vector<int> totalHamming(position[i].size());
			int g = numeric_limits<int>::max();
			vector<int> in;

			for (int k = 0; k < position[i].size(); ++k) {
				totalHamming[k] = 0;
				for (int j = 0; j < i; ++j) {
					if (position[j].size() == 1) {
						totalHamming[k] += hammingDistance(sequence[j], position[j][0], sequence[i], position[i][k], len);
					}
				}
				double amount = epsilon;
				if (epsilon < 1) amount = epsilon * g;

				if (totalHamming[k] < g + amount) {
					if (totalHamming[k] < g){
						g = totalHamming[k];
						int inSize = in.size();
						in.clear();
						for (int m = 0; m < inSize; ++m) {
							if (totalHamming[m] < g + amount) {
								in.push_back(m);	
							}
						}
					}
					in.push_back(k);
				}
			}

			for (int m = 0; m < in.size(); ++m) {
				in[m] = position[i][in[m]];
			}
			position[i].clear();
			position[i] = in;
		}
	
		epsilon *= 0.98;

		sequenceTemp = sequence;

		vector<vector<int>> tempPos = position;
		for (int i = 0; i < l ; ++i) {
			position[index[i]] = tempPos[i];
		}

		int posNumTmp = 0;
		for (int i = 0; i < position.size(); ++i) {
			posNumTmp += position[i].size();
	    }

		if (posNumTmp < 40) {
			if (posNumTmp * 2 < posNum || posNumTmp == position.size()) {
				break;
			}
		}
	}

	int topNumber = atoi(argv[4]);
	topHamming.assign(topNumber, numeric_limits<int>::max());
	topPosition.assign(topNumber, vector<int>(sequence.size()));

	int pos[l];
	for (int i = 0; i < l; ++i){
		pos[i] = 0;
	}

	int hamming = consensusCalculate(pos);
	top(topNumber, hamming, pos);

	int mark = l - 1;
	while (true) {
		while (mark >= 0 && pos[mark] >= (position[mark].size() - 1)) mark--;
		if (mark < 0) break;
		pos[mark]++;
		for (int i = mark + 1; i < l; ++i) {
			pos[i] = 0;
		}
		mark = l - 1;
		
		hamming = consensusCalculate(pos);
		top(topNumber, hamming, pos);
	}

	vector<int> posA;
	for (int i = 0; i < topNumber; ++i) {
		for (int j = 0; j < l; ++j) {
			posA.push_back(topPosition[i][j]);
			cout << topPosition[i][j] + 1 << ' ';
		}
		cout << '\n';
	}

	cout << "Complete's Motif: " << argv[2] << "\t";
    int csc = consensusScore(sequence, posA, len, l);
    cout << "Consensus Score: " << csc;
    cout << " Hamming Score: "<< len * l - csc << "\n";

	return 0;
}
