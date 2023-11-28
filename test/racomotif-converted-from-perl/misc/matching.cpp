/*
 ==> Description: Program list all the starting positions of a consensus in a set of input sequences
 ==> Input: 
	+ a set of input sequences stored in file argv[1]
	+ a consensus sequence stored in file argv[2]
 ==> Output: 
	+ list of all the starting position stored in file argv[3]
 ==> Usage:
	./matching <sequence file name> <consensus file name> <allowed threshold>
*/

#include<iostream>
#include<fstream>
#include<cstring>
#include<vector>
#include<cstdlib>
#include<climits>
#include<ctime>

using namespace std;

vector<vector<char> > sequence;
vector<char> consensus;
vector<vector<int> > position;
int ** topPosition;
int * topHamming;
/*
Description: read sequences and consensus
Input: two file names
Output: a two dimensional array storing sequences and a char array storing consensus
*/
void readFile(char * seq, char * con)
{
	ifstream in;
	// read sequence
	in.open(seq);
	string tmp;
	while(!in.eof()){
		vector<char> vec;
		getline(in, tmp);
		getline(in, tmp);
		for(int i = 0; i < tmp.length() - 1; ++i){
			vec.push_back(toupper(tmp[i]));
		}
		sequence.push_back(vec);
		getline(in, tmp);
	}
	sequence[sequence.size() - 1].push_back(toupper(tmp[tmp.size() - 1]));
	in.close(); 
	// read consensus
//	in.open(con);
//	getline(in, tmp);
	string tmp2(con);
	for (int i = 0; i < tmp2.size(); ++i){
		consensus.push_back(tmp2[i]);
	}
//	in.close();
}

/*
calculates Hamming distance score between two subsequence length l
input: two sequences equal in size
output: hamming distance between them
*/
int hammingDistance(vector<char> a, int starta, vector<char> b, int startb, int len)
{
	int tmp = 0;
	for (int i = 0; i < len; ++i){
		if ( a[starta + i] != b[startb + i]){
			++tmp;
		}
	}
	return tmp;
}

int consensusScore(vector<vector<char> >data, vector<int>start, int len, int n){
	int tmp = 0;
	int fre[len][4];	
	for (int i = 0; i < len; ++i){
		for (int j = 0; j < 4; ++j){
			fre[i][j] = 0;
		}
	}
	int x =0;
	for (int j = 0; j < len; ++j){
               switch(data[0][j]){
                       case 'A':
                               fre[j][0] += x;
                               break;
                       case 'T':
                               fre[j][1] += x;
                               break;
                       case 'C':
                               fre[j][2] += x;
                               break;
                       case 'G':
                               fre[j][3] += x;
               }
        }


	for (int i = 0; i < n; ++i){
//		cout << "Start " << start[i] << endl;
		for (int j = 0; j < len; ++j){
			switch(data[i][j + start[i]]){
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
	for (int i = 0; i < len; ++i){
		int max = fre[i][0];
		for (int j = 1; j < 4; ++j){
			if (fre[i][j] > max){
				max = fre[i][j];
			}
		}
//		cout << max << endl;
		tmp += max;
	}
//	cout << "Consensus: " << tmp << endl;
	return tmp;
}

/*
list all the posistion of the consensus with a threshold in hamming distance
*/
void determinePosition(int thres){
	int lastPos;
	int l = consensus.size();
	int score;
	int hamming;
	// scan all the sequences in the set
	for (int i = 0; i < sequence.size(); ++i){
		vector<int> n;
		position.push_back(n);
//		cout << "Working with sequence: " << i << ":" << sequence[i].size() << endl;
/*		for (int j = 0; j < sequence[i].size(); ++j){
                        cout << sequence[i][j];
                }
		cout << endl << "start finding" << endl;*/
		lastPos = sequence[i].size() - l;
		score = INT_MAX;
		// with each sequence we consider all the possible subsequence with length equal consensus sequence
		for (int j = 0; j <= lastPos; ++j){
			// calculate hamming distance for that subsequence
			hamming = hammingDistance(consensus, 0, sequence[i], j, consensus.size());
//			cout << "Hamming Score: "<< j << ":" << hamming << endl;
			// if hamming is smaller than current score then 
			// if hamming is smaller than current score minus permissive threshold then clear all the list positioned i
			// if not erase the position not in range
//			cout << "Current Score: " << score << endl;
			if (hamming < score){
//				cout << j << " Hamming < score" << endl;
//				cout << hamming << " < " << score << endl;
				score = hamming;
//				cout << "OK" << endl;
				vector<int> tmp;
				for (int k = 0; k < position[i].size(); ++k){
//					cout << position[i][k] << " " << hammingDistance(consensus, sequence[i], position[i][k]) << " ";
					if (hammingDistance(consensus, 0, sequence[i], position[i][k], consensus.size()) <= score + thres){
						tmp.push_back(position[i][k]);
					}
				}
//				cout << endl << endl;
				position[i].clear();
				position[i] = tmp;
				position[i].push_back(j);
				score = hamming;
			} else if (hamming <= score + thres){
//				cout << "Hamming < score + thres" << endl;
				position[i].push_back(j);
			}
		}
//		for (int k = 0; k < position[i].size(); ++k){
//			cout << position[i][k] << " ";
//		}
//		cout << endl;
	}
}

/*
list all the posistion of the consensus with a threshold in total hamming distance
*/
void determinePositionBasedOnTotalHamming(int thres){
	int lastPos;
	int l = consensus.size();
	int totalHamming;
	// scan all the sequences in the set
	for (int i = 0; i < sequence.size(); ++i){
		vector<int> n;
//		cout << "Working with sequence: " << i << ":" << sequence[i].size() << endl;
/*		for (int j = 0; j < sequence[i].size(); ++j){
                        cout << sequence[i][j];
                }
		cout << endl << "start finding" << endl;*/
		lastPos = sequence[i].size() - l;
		int totalHamming[lastPos];
		// with each sequence we consider all the possible subsequence with length equal consensus sequence
		for (int j = 0; j <= lastPos; ++j){
			//cout << "Total Hamming Distance " << j << endl;
			// calculate total hamming distance for that subsequence
			totalHamming[j] = 2*2*hammingDistance(consensus, 0, sequence[i], j, l);
			for (int u = 0; u < i; ++u){
				//cout << u << endl;
				totalHamming[j] += hammingDistance(sequence[u], position[u][0], sequence[i], j, l);
			}
		}
		int minScore = INT_MAX - thres;
		vector<int> minPos;
//		cout << "Sequence " << i << endl;
		for (int j = 0; j <= lastPos; ++j){
//			cout << j << ":" << totalHamming[j] << " ";
			if (totalHamming[j] <= minScore + thres){
				if (totalHamming[j] < minScore){
					minScore = totalHamming[j];
					vector<int> tmp;
					for (int u = 0; u < minPos.size(); ++u){
						if (totalHamming[minPos[u]] < minScore + thres){
							tmp.push_back(minPos[u]);
						}
					}
//					cout << "Create new vector" << endl;
					minPos = tmp;
					minPos.push_back(j);
				} else {
					minPos.push_back(j);
				}
			}
		}
//		cout << "Choose " << minPos.size() << endl;
//		cout << "Min Hamming: " << minScore << endl;
		position.push_back(minPos);
	}
}


/*
list all the posistion of the consensus with a threshold in total hamming distance
*/
void determinePositionBasedOnConsensus(int thres){
	int lastPos;
	int l = consensus.size();
	int totalHamming;
	// scan all the sequences in the set
	vector<vector<char> > data;
	data.push_back(consensus);
	vector<int> start;
	start.push_back(0);
	for (int i = 0; i < sequence.size(); ++i){
		data.push_back(sequence[i]);
//		vector<int> n;
//		n.push_back(0);
//		cout << "Working with sequence: " << i << ":" << sequence[i].size() << endl;
/*		for (int j = 0; j < sequence[i].size(); ++j){
                        cout << sequence[i][j];
                }
		cout << endl << "start finding" << endl;*/
		lastPos = sequence[i].size() - l;
		int totalHamming[lastPos];
		// with each sequence we consider all the possible subsequence with length equal consensus sequence
//		position.push_back(n);
		start.push_back(0);
		for (int j = 0; j <= lastPos; ++j){
//			cout << j << endl;
			start[i + 1] = j;
			// calculate total hamming distance for that subsequence
//			totalHamming[j] = 2*hammingDistance(consensus, 0, sequence[i], j, l);
			//cout << "Total Hamming Distance " << j << endl;
//			cout << start[i] << endl;
			totalHamming[j] = l * (i + 2) - consensusScore(data, start, l, i + 2);
		}
		int minScore = INT_MAX - thres;
		vector<int> minPos;
//		cout << "Sequence " << i << endl;
		for (int j = 0; j <= lastPos; ++j){
			//cout << j << ":" << totalHamming[j] << " ";
			if (totalHamming[j] <= minScore + thres){
				if (totalHamming[j] < minScore){
					minScore = totalHamming[j];
					vector<int> tmp;
					for (int u = 0; u < minPos.size(); ++u){
						if (totalHamming[minPos[u]] < minScore + thres){
							tmp.push_back(minPos[u]);
						}
					}
			//		cout << "Create new vector" << endl;
					minPos = tmp;
					minPos.push_back(j);
				} else {
					minPos.push_back(j);
				}
			}
		}
//		cout << "Choose " << minPos.size() << endl;
//		cout << "Min Hamming: " << minScore << endl;
//		position.pop_back();
		position.push_back(minPos);
		start[i + 1] = minPos[0];
	}
}



/*
calculates consensus of the found positions
*/
int consensusCalculate(int * pos){
	int l = consensus.size();
	int select = sequence.size();
	int tmp[l][4];
	for (int i = 0; i < l; ++i){
                tmp[i][0] = tmp[i][1] = tmp[i][2] = tmp[i][3] = 0;
        }

	for (int i = 0; i < sequence.size(); ++i){
		int a = position[i][pos[i]];
//		cout << i << ":" << a << endl;
		for (int j = a; j < a + l; ++j){
//			cout << j << " ";
			if (sequence[i][j] == 'A'){
				++tmp[j - a][0];
			} else if (sequence[i][j] == 'T'){
                                ++tmp[j - a][1];
                        } else if (sequence[i][j] == 'C'){
                                ++tmp[j - a][2];
                        } else if (sequence[i][j] == 'G'){
                                ++tmp[j - a][3];
                        }
		}
	}
	vector<char> candidate;
	int hamming = 0;
	for (int i = 0; i < l; ++i){
//		cout << tmp[i][0] << " " << tmp[i][1] << " " << tmp[i][2] << " " << tmp[i][3] << endl;
		if (tmp[i][0] >= tmp[i][1]){
			if (tmp[i][0] >= tmp[i][2]){
				if (tmp[i][0] >= tmp[i][3]){
					candidate.push_back('A');
					hamming += select - tmp[i][0];
//					cout << 'A' << endl;
				} else {
					candidate.push_back('G');
					hamming += select - tmp[i][3];
//					cout << 'G' << endl;
				}
			} else {
				if (tmp[i][2] >= tmp[i][3]){
					candidate.push_back('C');
					hamming += select - tmp[i][2];
//					cout << 'C' << endl;
				} else {
					candidate.push_back('G');
					hamming += select - tmp[i][3];
//					cout << 'G' << endl;
				}
			}
		} else {
			if (tmp[i][1] >= tmp[i][2]){
                                if (tmp[i][1] >= tmp[i][3]){
                                        candidate.push_back('T');
					hamming += select - tmp[i][1];
//					cout << 'T' << endl;
                                } else {
                                        candidate.push_back('G');
					hamming += select - tmp[i][3];
//					cout << 'G' << endl;
                                }
                        } else {
                                if (tmp[i][2] >= tmp[i][3]){
                                        candidate.push_back('C');
					hamming += select - tmp[i][2];
//					cout << 'C' << endl;
                                } else {
                                        candidate.push_back('G');
					hamming += select - tmp[i][3];
//					cout << 'G' << endl;
                                }
                        }

		}
	}
//	cout << hamming << endl;
//	if (consensus == candidate){
/*		for(int i = 0; i < sequence.size(); ++i){
        	        cout << pos[i] << " ";
	        }
		cout << endl;

		for(int i = 0; i < candidate.size(); ++i){
        	        cout << candidate[i];
	        }

		cout << "Bingo" << endl;*/
		return hamming;
//	}
//	return -1;
/*	for(int i = 0; i < sequence.size(); ++i){
		cout << pos[i] << " ";
	}
	cout << endl;

	for(int i = 0; i < candidate.size(); ++i){
                cout << candidate[i];
        }

	cout << endl;

//	cout << hammingDistance(consensus, candidate, 0) << endl;
	if (hammingDistance(consensus, candidate, 0) == 0){
		cout << "Bingo" << endl;
	}
	*/
}

/*
top best positions found
*/
void top(int num, int score, int * pos)
{
	
	int i = num - 1;
	int len = sequence.size();	
	while (topHamming[i] > score && i > 0){
		topHamming[i] = topHamming[i - 1];
		for (int j = 0; j < len; ++j){
			topPosition[i][j] = topPosition[i - 1][j];
		}
		--i;
	}
	//cout << i << endl;
	if (i == 0){
		//cout << endl << "Base on Hamming" << endl;
		if (score < topHamming[0]){
			topHamming[0] = score;
			for (int j = 0; j < len; ++j){
				topPosition[0][j] = position[j][pos[j]];
			}
		} else {
			if (num > 1){
			topHamming[1] = score;
				for (int j = 0; j < len; ++j){	
					topPosition[1][j] = position[j][pos[j]];
				}
			}
		}
	} else if (i < num - 1){
		topHamming[i + 1] = score;
		for (int j = 0; j < len; ++j){
			topPosition[i + 1][j] = position[j][pos[j]];
		}
	}
	
}


int main(int argc, char * argv[])
{
//	cout << "Finding Position..." << argv[1] << " " << argv[2] << " " << argv[3] << " " << argv[4] << endl;
	clock_t start = clock();
	readFile(argv[1], argv[2]);
//	cout << argv[2] << ":" << strlen(argv[2]) << endl;
//	for (int i = 0; i < strlen(argv[2]); ++i){
//		cout << argv[2][i] << endl;
//		consensus.push_back(argv[2][i]);
//		cout << argv[2][i] << endl;
//	}
//	cout << argv[2] << ":" << consensus.size() << endl;
//	for (int i = 0; i < strlen(argv[2]); ++i){
//		cout << consensus[i] << " ";
//	}
//	consensus = argv[2];
	int l = sequence.size();
	int len = consensus.size();
	int candidate[10][l];
//	cout << "Sequence Size: " << sequence.size() << endl;
/*	for(int i = 0; i < sequence.size(); ++i){
		cout << "Seq " << i << endl;
		for (int j = 0; j < sequence[i].size(); ++j){
			cout << sequence[i][j];
		}
                cout << endl;
        }
	
	cout << "Consensus: ";
	for (int i = 0; i < consensus.size(); ++i){
		cout << consensus[i];
	}
	cout << endl;
*/

        //cout << endl << "Base on Nothing" << endl;
        position.clear(); 
        determinePosition(atoi(argv[3]));
	int posNum = 0;
        for(int i = 0; i < position.size(); ++i){
		posNum += position[i].size();
/*	        for (int j = 0; j < position[i].size(); ++j){
                        cout << position[i][j] << " ";
                }
                cout << endl;*/
        }
//	cout << "Position Number: " << posNum << endl;

	vector<vector<int> > posit = position;
	int iter = 0;
	vector<vector<char> > sequenceTemp = sequence;
	float frac = 0.5;
	do{
		if (posNum < 2*position.size()){
			frac = 0;
		}
		//posit = position;
		int index[l];
		vector<int> tmp;
		for(int i = 0; i < l; ++i){
			index[i] = i;
		}
		for(int i = 0; i < l; ++i){
			for (int j = i + 1; j < l; ++j){
				if (position[i].size() > position[j].size()){
					tmp = position[i];
					position[i] = position[j];
					position[j] = tmp;
					int t = index[i];
					index[i] = index[j];
					index[j] = t;
				}
			}
		}
		int firstOccurent = 0;
		for(int i = 0; i < l; ++i){
			sequence[i] = sequenceTemp[index[i]];
			if (position[i].size() == 1){
				firstOccurent = i;
			}
/*		        for (int j = 0; j < position[i].size(); ++j){
		                cout << position[i][j] << " ";
		        }
		        cout << endl;*/

		}
		//cout << 10 * frac << endl;
		for (int i = firstOccurent + 1; i < l; ++i){
			int posSize = position[i].size();
			int totalHamming[posSize];
			int min = INT_MAX;
			vector<int> in;
			for (int k = 0; k < posSize; ++k){
				totalHamming[k] = 0;
				for (int j = 0; j < i; ++j){
					if (position[j].size() == 1){
						totalHamming[k] += hammingDistance(sequence[j], position[j][0], sequence[i], position[i][k], len);
					}
				}
				float amount = frac;
				if (frac < 1){
					amount = frac*min;
				}
				//cout << "Fraction: " << amount << endl;
				if (totalHamming[k] < min + amount){
					if (totalHamming[k] < min){
						min = totalHamming[k];
						int inSize = in.size();
						in.clear();
						for (int m = 0; m < inSize; ++m){
							if (totalHamming[m] < min + amount){
								in.push_back(m);	
							}
						}
					}
					in.push_back(k);
				}
				//cout << k << ":" << totalHamming[k] << endl;
			}
			/*cout << "Hamming List: " << endl;
			for (int m = 0; m < in.size(); ++m){
				cout << position[i][in[m]] << ":" << totalHamming[in[m]] << " ";
			}
			cout << endl;
			if (in.size() > 2){
				for (int m = 0; m < 2; ++m){
					for (int n = 0; n < in.size(); ++n){
						if (totalHamming[in[m]] > totalHamming[in[n]]){
							min = in[m];
							in[m] = in[n];
							in[n] = min;
						}
					}
				}
				in.resize(2);
			}*/
		
			for (int m = 0; m < in.size(); ++m){
				in[m] = position[i][in[m]];
			}
			position[i].clear();
			position[i] = in;
			/*for (int k = 0; k < l; ++k){
				for (int j = 0; j < position[k].size(); ++j){
			                cout << position[k][j] << " ";
		       		}
			        cout << endl;
			}*/

		}
	
		vector<vector<int> > tempPos = position;
		sequenceTemp = sequence;
		for(int i = 0; i < l ; ++i){
	//		cout << index[i] << " ";
			position[index[i]] = tempPos[i];
	//		sequence[index[i]] = sequenceTemp[i];
		}
		tempPos.clear();
/*		cout << "Two matrix: " << endl;
		for(int i = 0; i < posit.size(); ++i){
	                for (int j = 0; j < posit[i].size(); ++j){
	                        cout << posit[i][j] << " ";
	                }
	                cout << endl;
	        }*/
//		cout << endl;
		int posNumTmp = 0;
		for(int i = 0; i < position.size(); ++i){
			posNumTmp += position[i].size();
/*	                for (int j = 0; j < position[i].size(); ++j){
	                        cout << position[i][j] << " ";
	                }
	                cout << endl;*/
	        }
//		cout << "Temporary Position Number: " << posNumTmp << endl;
		if (posNumTmp < 40){
			if ((posNumTmp * 2 < posNum) || posNumTmp == position.size()){
				break;
			}
		}
		frac *= 0.98;
//		position = posit;
		//cin.get();
		//iter++;
	} while (true);
	//} while (position != posit);

//	cout << endl << "Result: " << endl;

	
//	cout << endl << "Base on Hamming" << endl;
//	position.clear();

/*	determinePositionBasedOnTotalHamming(atoi(argv[3]));
	for(int i = 0; i < position.size(); ++i){
                for (int j = 0; j < position[i].size(); ++j){
                        cout << position[i][j] << " ";
                }
                cout << endl;
        }
*/
/*	cout << endl << "Base on Consensus: " << endl;
        position.clear();
        determinePositionBasedOnConsensus(atoi(argv[3]));
        for(int i = 0; i < position.size(); ++i){
                for (int j = 0; j < position[i].size(); ++j){
                        cout << position[i][j] << " ";
                }
                cout << endl;
        }*/
	// store the top of the positions
	int topNumber = atoi(argv[4]);
	topHamming = new int[topNumber];
	topPosition = new int *[topNumber];
	for (int i = 0; i < topNumber; ++i){
		topHamming[i] = INT_MAX;
		topPosition[i] = new int[sequence.size()];
	}
	int pos[l];
//	cout << "Hamming 16: " << hammingDistance(consensus, sequence[16], 83) << endl;
//	cout << "Hamming 1: " << hammingDistance(consensus, sequence[1], 54) << endl;
//	cout << "Hamming 3: " << hammingDistance(consensus, sequence[3], 62) << endl;
//	consensusCalculate(pos);
	for (int i = 0; i < l; ++i){
		pos[i] = 0;
	}
	int hamming = consensusCalculate(pos);
	if (hamming >= 0){
		top(topNumber, hamming, pos);
	}
	int mark = l - 1;
	//cout << endl << "Base on Hamming" << endl;
	//cout << position[mark].size() << endl;
	while (1){
		while(mark >= 0 && pos[mark] >= (position[mark].size() - 1)){
			--mark;
		}
		if (mark < 0)
			break;
		++pos[mark];
		for (int i = mark + 1; i < l; ++i){
			pos[i] = 0;
		}
		mark = l - 1;
		
		hamming = consensusCalculate(pos);
		//cout << endl << "Base on Hamming" << endl;
		if (hamming >= 0){
			top(topNumber, hamming, pos);
		}
	}

	//cout << "Result: " << topNumber << ":" << endl;
	vector<int> posA;
//	cout<<"Position: ";
	for (int i = 0; i < topNumber; ++i){
//		cout << i << ":" << topHamming[i] << ":";
		for (int j = 0; j < l; ++j){
			posA.push_back(topPosition[i][j]);
			cout << topPosition[i][j]+1 << " ";
		}
		cout << endl;
	}
	//cout << endl;

	cout<<"Hoang's Motif: "<<argv[2]<<"\t";
        int x = consensusScore(sequence, posA, len, l);
        cout<<"Consensus Score: "<<x;
        cout<<" Hamming Distance: "<<len*l-x<<"\n";


	sequence = sequenceTemp;
//	cout << "Execution Time: " << clock() - start << endl;

	sequence.clear();
	sequenceTemp.clear();
	consensus.clear();
	position.clear();
	delete(topHamming);
	delete [] (topPosition);

	return 0;
}
