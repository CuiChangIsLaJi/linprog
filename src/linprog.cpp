#include "../include/linprog.h"
#include <eigen3/Eigen/Dense>
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <algorithm>

using namespace std;

void Solver::parse_input(string filename){
	fstream input_file;
	int n_eq, n_ineq, n_vars;
	try{
		input_file.open(filename.c_str(),fstream::in);
		if(!input_file) throw invalid_argument("Input file does not exist.");
	}
	catch(invalid_argument err){
		cout << err.what() << endl;
		exit(0);
	}
	input_file >> n_eq >> n_ineq >> n_vars;
	n_variables = n_vars;
	n_constraints = n_eq + n_ineq;
	_n_variables = n_vars + n_ineq;
	A.resize(n_constraints, _n_variables);
	b.resize(A.rows());
	for (int i = 0; i < A.rows(); i++){
		for (int j = 0; j < A.cols(); j++){
			input_file >> A(i,j);
		}
		input_file >> b(i);
	}
        x.resize(A.cols());
        c.resize(A.cols());
	for (int k = 0; k < A.rows(); k++){
		input_file >> c(k);
	}
	_M.resize(A.rows(), A.rows());
	_omega.resize(A.rows());
	_b_tilde.resize(A.rows());
	_y.resize(A.rows());
	_xM.resize(A.rows());
	_cM.resize(A.rows());
	_base_indices.resize(A.rows());
	_is_base_variable.resize(A.cols());
	_discriminants.resize(A.cols());
	_M = A.rightCols(A.rows());
	_cM = c.tail(A.rows());
	for (int i = A.rows() - 1; i > -1; i--){
		_base_indices[i] = i - A.rows() + A.cols();
	}
	cout << endl;
	for (int j = A.cols() - 1; j > -1; j--){
		if (j >= A.cols() - A.rows())	_is_base_variable[j] = true;
		else _is_base_variable[j] = false;
	}
	input_file.close();
}

bool Solver::check_b(){
	bool check = true;
	for (int i = 0; i < b.size(); i++){
		if (b(i) < 0) check = false;
	}
	return check;
}

void Solver::solve(string output_file){
	int i = 0, k;
	fstream file;
	file.open(output_file.c_str(), fstream::out);
	while(true){
		k = _step();
		file << "Iteration No. " << i + 1 << endl;
		file << "x =\n" << x << endl;
		i++;
		if (k == -1)	break;
	}
	file << "Final result:\n" << x << endl;
}

void Solver::report(){
	cout << "Result of simplex method computation:" << endl;
	double sum = 0;
	for (int i = 0; i < n_variables; i++){
		cout << "x[" << i + 1 << "] = " << x(i) << endl;
		sum += c(i) * x(i);
	}
	cout << "Minimum value = " << sum << endl;
}

void Solver::_update_x(){
	x.fill(0);
        for (int i = 0; i < A.rows(); i++){
		x(_base_indices[i]) = _xM(i);
	}
}

void Solver::_update_base_components(int j_old, int j_new){
	int j = _base_indices[j_old];
	_base_indices[j_old] = j_new;
	_cM(j_old) = c(j_new);
	for (int i = 0; i < n_constraints; i++){
		_M(i,j_old) = A(i,j_new);
	}
	_is_base_variable[j] = false;
	_is_base_variable[j_new] = true;
}

int Solver::_step(){
	int k, r;
	_b_tilde = _M.inverse() * b;
	_xM = _b_tilde;
	_update_x();
	z = _cM.dot(_xM);
	_omega = _M.inverse().transpose() * _cM;
	for (int j = 0; j < A.cols(); j++){
		if (_is_base_variable[j])	_discriminants(j) = 0;
		else _discriminants(j) = _omega.dot(A.col(j)) - c(j);
	}
	k = _search_max_discriminant();
	if (k != -1){
		_y = _M.inverse() * A.col(k);
		r = _assign_xk();
		if (r == -1){
			cout << "Finite solution does not exist." << endl;
			exit(0);
		}
		else{
			x(k) = _b_tilde(r) / _y(r);
			_update_base_components(r, k);
		}
	}
	return k;
}

int Solver::_search_max_discriminant(){
	int k = -1;
	double max = 0;
	for (int j = 0; j < A.cols(); j++){
                if (!_is_base_variable[j]){
			if (_discriminants(j) > max){
				max = _discriminants(j);
				k = j;
			}
		}
        }
	return k;
}

int Solver::_assign_xk(){
	int r = -1;
	vector<int> candidates;
	double min;
	for (int i = 0; i < A.rows(); i++)
		if (_y(i) > 0)	candidates.push_back(i);
	min = _b_tilde(candidates[0]) / _y(candidates[0]);
	for (int j = 0; j < candidates.size(); j++){
		if (_b_tilde(candidates[j]) / _y(candidates[j]) <= min){
			min = _b_tilde(candidates[j]) / _y(candidates[j]);
			r = candidates[j];
		}
	}
	return r;
}
