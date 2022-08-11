#pragma once
#include <eigen3/Eigen/Dense>
#include <string>
#include <fstream>
#include <vector>

using namespace Eigen;
using namespace std;

class Solver{
public:
	Matrix<double, Dynamic, Dynamic> A;
        VectorXd x, b, c;
	double z;
	int n_constraints, n_variables;
	void parse_input(string);
	bool check_b();
	void solve(string);
	void report();
private:
	Matrix<double, Dynamic, Dynamic> _M;
	VectorXd _omega, _b_tilde, _y, _xM, _cM, _discriminants;
	vector<int> _base_indices;
	vector<bool> _is_base_variable;
	int _n_variables;
	void _update_x();
	void _update_base_components(int, int);
	int _step();
	int _search_max_discriminant();
	int _assign_xk();
};


