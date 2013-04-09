#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <functional>
#include "NGS++.h"

using namespace std;
using namespace NGS;

const int BIN_SIZE = 200;
const int MIN_SUM = 30;

void addRegion(vector<int>& raw_data_plus, vector<int>& raw_data_minus, uToken& token);
void binChr(const std::string& chr_name, const std::string& chr_size,  vector<int>& raw_data_plus,  vector<int>& raw_data_minus, float averageLength);
void resizeDataStructure(vector<int>& raw_data_plus, vector<int>& raw_data_minus, vector<string>& chr_size, int currentIndex);

int main(int argc, char* argv[]) {
	if (argc != 2) {
		cout << "Usage: DensityStrand filename.sam" << endl;
		cout << "Note: sam file must be sorted!" << endl;
	}
	else {
		// Fetch header infos
		uParser parser(argv[1], "SAM");	
		vector<string> chr_names = parser.getHeaderParamVector(header_param::CHR);
		vector<string> chr_size = parser.getHeaderParamVector(header_param::CHR_SIZE);

		// Create data structures
		vector<int> raw_data_plus;
		vector<int> raw_data_minus;

		// Parse the file
		int currentIndex = 0;
		long lengthSum = 0;
		int numberOfElem = 0;
		resizeDataStructure(raw_data_plus, raw_data_minus, chr_size, currentIndex);
		while (!parser.eof()) {
			uToken token = parser.getNextEntry();
			if (token.getParam(token_param::CHR) != chr_names[currentIndex]) { // Note: file must be sorted!
				if (numberOfElem > 0) {
					float averageLength = (float)(lengthSum) / (float)(numberOfElem);
					binChr(chr_names[currentIndex], chr_size[currentIndex], raw_data_plus, raw_data_minus, averageLength);
					currentIndex++;
					resizeDataStructure(raw_data_plus, raw_data_minus, chr_size, currentIndex);
					lengthSum = 0;
					numberOfElem = 0; 
				} 
			} 
			addRegion(raw_data_plus, raw_data_minus, token); stringstream ss2;
			ss2 << token.getParam(token_param::END_POS) << "\t" << token.getParam(token_param::START_POS);
			long end = 0;
			long start = 0;
			ss2 >> end >> start;
			lengthSum += end - start;
			numberOfElem++;
		}
		float averageLength = (float)(lengthSum) / (float)(numberOfElem);
		binChr(chr_names[currentIndex], chr_size[currentIndex], raw_data_plus, raw_data_minus, averageLength);
	}
	return 0;
}

void addRegion(vector<int>& raw_data_plus, vector<int>& raw_data_minus, uToken& token) {
	stringstream ss;
	ss << token.getParam(token_param::START_POS) << "\t" << token.getParam(token_param::END_POS);
	ss << "\t" << token.getParam(token_param::STRAND);
	int begin;
	int end;
	string strand;
	ss >> begin >> end >> strand;
	for (int i = begin; i < end; i++) {
		if (strand == "+") {
			raw_data_plus[(int)(i)]++;
		}
		else if (strand == "-") {
			raw_data_minus[(int)(i)]++;
		}
	}
}

void binChr(const std::string& chr_name, const std::string& chr_size,  vector<int>& raw_data_plus,  vector<int>& raw_data_minus, float averageLength) {
	// Open file to save
	stringstream ss; 
	ss << chr_name << ".txt";
	ofstream out;
	out.open(ss.str().c_str(), ios_base::out);
	// Parse each bins
	stringstream ss2;
	ss2 << chr_size;
	int size;
	ss2 >> size;
	for (int i = 0; i < size - BIN_SIZE; i++) {
		int sum_plus = 0;
		int sum_minus = 0;
		for (int j = i; j < i + BIN_SIZE; j++) {
			sum_plus += raw_data_plus[j];
			sum_minus += raw_data_minus[j];
		}
		if (((sum_plus + sum_minus) / averageLength) >= MIN_SUM) {
			out << (float)(sum_plus)/averageLength << "\t" << (float)(sum_minus)/averageLength << endl;
		}
	}
	out.close();
}

void resizeDataStructure(vector<int>& raw_data_plus, vector<int>& raw_data_minus, vector<string>& chr_size, int currentIndex) {
	// Reset vectors
	raw_data_plus.clear();
	raw_data_minus.clear();
	// Get size
	stringstream ss;
	ss << chr_size[currentIndex];
	int size;
	ss >> size;
	// Resize
	raw_data_plus.resize(size);
	raw_data_minus.resize(size);
}
