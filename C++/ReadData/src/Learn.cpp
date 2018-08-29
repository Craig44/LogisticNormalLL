//============================================================================
// Name        : Learn.cpp
// Author      : C.Marsh
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <map>
#include <stdlib.h>
#include <fstream>
#include <sstream>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

MatrixXd readMatrix(const string& filename);

int main() {
	cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!

	//MatrixXd obs_mat = readMatrix("observed_data.txt");
	MatrixXd temp_exp_mat = readMatrix("error_data.txt");
	MatrixXd exp_mat(temp_exp_mat.rows(),temp_exp_mat.cols());
	exp_mat = temp_exp_mat;

	double log_sum = exp_mat.array().log().sum();
	cout << "rows = " << exp_mat.rows() << " cols = " << exp_mat.cols() << endl;


	int N_bins = 20;
	MatrixXd covar(N_bins,N_bins);
	covar.setRandom();
	int ages = 2;
	cout << covar.block(ages,ages,ages,ages) << endl;
	MatrixXd subcovar(ages,ages);
	subcovar.setZero();
	covar.block(ages,ages,ages,ages) = subcovar;

	cout << covar.block(ages,ages,ages,ages) << endl;

	cout << "finished" << endl;

	return 0;
}


MatrixXd readMatrix(const string& filename) {
	ifstream infile;
	infile.open(filename);
	unsigned nrows = 0;
	unsigned ncols = 1;
	if (infile.is_open()) {
		  //cout << "start reading the file" << endl;
		  string current_line;

		  while(getline(infile, current_line)) {
			  if (nrows == 0) {
				  // count cols
				  size_t pos = 0;
				  string token;
				  string delimiter = " ";
				  while ((pos = current_line.find(delimiter)) != std::string::npos) {
					  token = current_line.substr(0, pos);
					  std::cout << token << std::endl;
					  current_line.erase(0, pos + delimiter.length());
					  ++ncols;
				  }
			  }
			  ++nrows;
		  }
		  infile.close();
	} else {
		  cout << "could not open the file" << endl;
	}
	cout << "matrix had n-rows = " << nrows << " and cols = " << ncols;

	/// Now create a matrix and save the informations
	MatrixXd return_mat(nrows,ncols);
	infile.open(filename);
	string current_line;
	unsigned row = 0;
	unsigned col = 0;

	while(getline(infile, current_line)) {
		  stringstream stream(current_line);
		  col = 0;
		  while(! stream.eof()) {
			stream >> return_mat(row,col);
			++col;
		  }
		  ++row;
	}
	//cout << "mat" << endl;
	return return_mat;
}
