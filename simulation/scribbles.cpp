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


}
