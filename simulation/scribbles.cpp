#include <iostream>
#include <numeric>
#include <vector>
#include "Eigen/Dense"
#include "sph_sim.h"

using namespace Eigen;
using namespace std;

double myfunction (double a, double b) {
	return a + b;
}
int main()
{
	VectorXd v1(5);
	v1 << 1,2,3,4,5;
	int N = 0;
	N += v1.sum();
	cout << "N: " << N << endl << endl;

	MatrixXd m1(10,1);
	m1 << 1,2,44,4,5,6,7,8,9,0;
	MatrixXi I(5,1);
	I << 0,2,4,6,8;
	double max_m1 = I.unaryExpr([&](int x) {
		return m1(x);
	}).maxCoeff();
	cout << "max_m1: "<< max_m1 << endl << endl;



}
