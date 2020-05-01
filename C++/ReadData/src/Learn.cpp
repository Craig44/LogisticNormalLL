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
#include <list>
#include <random>
//#include <numeric>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

MatrixXd readMatrix(const string& filename);

template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}


class WorldCell {
public:
  // Methods
  WorldCell() = default;
  virtual                     ~WorldCell() = default;
  // iterator accessors for agents
  vector<unsigned>::iterator  agent_begin() { return agents_.begin();}
  vector<unsigned>::iterator  agent_end() { return agents_.end();}
  unsigned                    agent_size() {return agents_.size();};
  // iterator accessors for tagged agents
  vector<unsigned>::iterator  tag_agent_begin() { return tag_agents_.begin();}
  vector<unsigned>::iterator  tag_agent_end() { return tag_agents_.end();}
  unsigned                    tag_agent_size() {return tag_agents_.size();};

  // iterator accessors for ALL agents
  vector<unsigned>::iterator  begin() { return tag_agents_.begin();}
  vector<unsigned>::iterator  end() { return tag_agents_.end();}
  unsigned                    size() {return tag_agents_.size();};

  // iterator accessors for ALL agents
  vector<unsigned>::iterator  begin() { return vector_1_.begin();}
  vector<unsigned>::iterator  end() { return vector_2_.end();}

  for(unsigned i = 0; i < vector_1_.size(); ++i) {
	  Do some function on vector_1_[i]
  }



  // Overload the ++ operator to deal with traversing the end element
  vector<unsigned>::iterator()++;
protected:
  vector<unsigned> 			vector_1_;
  vector<unsigned> 			vector_2_;
};



int main() {
	cout << "\n!!!Hello World!!!" << endl; // prints !!!Hello World!!!

	unsigned max_attempt = 100;
	unsigned attempt = 0;
	unsigned count_values_less_than_10 = 0;
	unsigned values_we_want = 1000;
	vector<unsigned> iterators_stored;
	random_device rd;  //Will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    uniform_real_distribution<> dis(0.0, 100.0);
    double val = 0;
	while (count_values_less_than_10 <= values_we_want) {
		attempt++;
		cout << attempt << " ";
		if (attempt > max_attempt)
			break;
		val = dis(gen);
		if (val >= 10) {
			continue;
		}
		iterators_stored.push_back(attempt);
		count_values_less_than_10++;
	}

	cout << "\n\n\n\nvalues stored = " << iterators_stored.size() << endl;

	vector<int> ages;
	vector<int> ages1;
	for(int age = 0; age < 21; ++age) {
		if (age > 0)
			ages1.push_back(age);
		ages.push_back(age);
	}

	int spread = 20 - 0 + 1;
	int spread1 = 20 - 1 + 1;
	cout << "spread = " << spread<< " size = " << ages.size() << endl;
	cout << "spread1 = " << spread1 << " size = " << ages1.size() << endl;


	return 0;

	vector<vector<double>> biomass;
	vector<vector<double>> order_mat;
	int rows =4;
	int cols = 3;
	vector<double>	biomass_vec(rows * cols, 0.0);
	vector<double>	effort_vec(rows * cols, 0.0);
	vector<double>	ordered_error(rows * cols, 0.0);

	biomass.resize(rows);
	order_mat.resize(rows);
	for (int i = 0; i < rows; ++i) {
		biomass[i].resize(cols,0.0);
		order_mat[i].resize(cols,0);
		for (int j = 0; j < cols; ++j) {
			biomass[i][j]= (rand() % 1000);
			cout << biomass[i][j] << " ";
			biomass_vec[i*cols + j] = biomass[i][j];
			effort_vec[i*cols + j]= (rand() % 1000);
		}
		cout << "\n";
	}
	cout << "\n";



	cout << "\n";

//	for(auto val : biomass_vec)
//		cout << val << " ";
//
	//cout << "\n\neffort vec ordered = " << endl;
	sort(effort_vec.begin(),effort_vec.end());


	cout << "biomass\n";
		for(auto val : biomass_vec)
			cout << val << " ";
	cout << "\n";

	cout << "effort\n";
	for(auto val : effort_vec)
		cout << val << " ";
	cout << "\n";

	//sort(temp_vec.begin(),temp_vec.end());
	//for(auto val : temp_vec)
	//	cout << val << " ";
	//cout << "\n";


	vector<size_t> biomass_indices = sort_indexes(biomass_vec);

	cout << "sorted biomass\n";
	for(auto val : biomass_indices)
		cout << biomass_vec[val] << " ";
	cout << "\n";
	cout << "biomass index\n";
	for(auto val : biomass_indices)
		cout << val << " ";
	cout << "\n";
	unsigned counter = 0;

	for(unsigned i = 0; i < ordered_error.size(); ++i) {
		ordered_error[biomass_indices[i]] = effort_vec[i];
	}
	cout << "Effort where it should be\n";
	for(auto val : ordered_error)
		cout << val << " ";
	cout << "\n";


	counter = 0;
	for (auto i: biomass_indices) {
	  cout << biomass_vec[counter] << " " << effort_vec[i] << endl;
	  ++counter;
	}

	cout << "\nsorted \n";
	counter = 0;
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			cout << biomass_vec[i*cols + j] << " " << biomass_vec[biomass_indices[counter]]  << " " << effort_vec[biomass_indices[counter]] <<"  " << biomass_indices[counter] << " ";
			++counter;
		}
		cout << "\n";
	}
	cout << "\n";
	return 0;
	// Don't care pass here


	int A = 15;
   vector<double> pred_N(A,0.0);
   vector<double> numbers_at_age(A,0.0);
   numbers_at_age = {22464300, 20277814, 18304142, 16522571, 14914403, 13462760, 12152409, 10969595,  9901907,  8938139 , 8068176 , 7282887,6574032 , 5934171 ,55034435};
   for (int y = 0; y < 30; ++y) {
	   pred_N[0] = 22464300;
	   double M = 0.1024;
		// Ageing and Z
		for (int age = 1; age < A; ++age) {
		  //cout << "age = " << age << endl;
		  pred_N[age] = numbers_at_age[age - 1] * exp(-M);
		}
		// deal with plus group
		cout << "plus group current numbers from previous age = " << pred_N[A - 1] << " current numbers in plus group = " << numbers_at_age[A - 1];
		pred_N[A - 1] += numbers_at_age[A - 1] * exp(-M);
		cout << " after = " << pred_N[A - 1] << endl;
		numbers_at_age = pred_N;
   }




	std::list<unsigned> practice_list;
	for (unsigned i = 0; i < 30; ++i) {
		practice_list.push_back(i);
	}

	cout << "size = " << practice_list.size() << "\n";
	counter = 0;
	auto first_iter = practice_list.begin();
	for (auto iter = practice_list.begin(); iter != practice_list.end(); ++counter) {
		if (counter == 18) {
			practice_list.insert(practice_list.begin(),200);
		}
		if (counter > 15) {
			practice_list.erase(iter);
		}

		++iter;
	}
	cout << "look at values now " << endl;
	for (auto iter = practice_list.begin(); iter != practice_list.end(); ++iter) {
		cout << *iter <<  " ";
	}
	cout << "\n";
	//MatrixXd obs_mat = readMatrix("observed_data.txt");
	MatrixXd temp_exp_mat = readMatrix("error_data.txt");
	MatrixXd exp_mat(temp_exp_mat.rows(),temp_exp_mat.cols());
	exp_mat = temp_exp_mat;

	double log_sum = exp_mat.array().log().sum();
	cout << "rows = " << exp_mat.rows() << " cols = " << exp_mat.cols() << endl;

/*
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

	return 0;*/
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
