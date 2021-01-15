#include <iostream>
#include "sph_sim.h"
#include "matrix.h"

using namespace std;

int main() {

	cout << "test" << endl;
	
	// Testing stuff
	
	matrix m(2,2);
	m.print();
	m = 3 + m; m.print(); // 3
	m = m + 3; m.print(); // 6
	m = m - 2; m.print(); // 4
	m = 8 - m; m.print(); // 4
	m = m * 3; m.print(); // 12
	m = 2 * m; m.print(); // 24
	m = m / 3; m.print(); // 8
	m = 16 / m; m.print(); // 2
	m = m + m*2; m.print(); // 6
	m = pow(m,2); m.print(); // 36
	m = m - 72/m; m.print(); // 34
	m = m/m; m.print();	// 1
	matrix n(m.dim()); n.print();
	matrix o = (n<3); o.print();
	o = o*2; o.print();
	m = ( ( 1.0-3.0/2.0*pow(o,2.0) + 3.0/4.0*pow(o,3.0) )); m.print();
	m(3,4); cout << endl; m.print();
	m(5) = 4; m(16) = 9; cout << endl; m.print();
	cout << m.max() << endl;
	matrix mr = m(range(0,7)); cout << endl; mr.print();
	cout << endl;
	o(3,4) = 1; cout << endl; o.print();
	m += 1; cout << endl; m.print();
	matrix mo = (m == o); cout << endl; mo.print();
	matrix fmo = find(mo); cout << endl; fmo.print();

	matrix ma = {{1,2},{3},{4,5,6,7}}; cout << endl; ma.print();

	matrix midx = {{0, 1, 2}, {5, 8, 11}};
	matrix ma_midx = ma(midx); cout << endl; ma_midx.print();
	matrix fma = find(ma); cout << endl; fma.print();
	matrix ma_fma = ma(fma); cout << endl; ma_fma.print();
	double ma_fma_max = ma(fma).max(); cout << endl << "ma_fma.max() = " << ma_fma_max << endl;

	matrix mm1 = {{1,2,3},{4,5,6}};
	matrix mm2 = {{1,2},{3,4},{5,6}};
	matrix mm1_mm2 = mm1*mm2; cout << endl; mm1_mm2.print();

	matrix cv = {{0},{1},{2},{3},{4},{5}}; cout << endl; cv.print();
	matrix rv = {{0,1,2,3,4,5}}; cout << endl; rv.print();
	cout << endl; cv(range(0,4)).print();
	cout << endl; (cv(range(1,4)) * ones(1,4)).print();
	cout << endl; (rv(range(1,4)) * cv(range(1,4))).print();
	cout << endl; (cv(range(1,4)) * ones(1,4)).T().print();

	cout << endl; max({{1},{2},{3},{4}},{{1,3,2,5}}).print();

	matrix mapp = {{1,2},{3,4},{5,6}};
	cout << endl; mapp.append_right({{1,7},{2,8},{3,9}}); mapp.print();
	cout << endl; mapp.append_down({{1,2,3,4},{5,6,7,8}}); mapp.print();

	cout << endl; mlogical_and({{1,2},{3,4},{5,6}},{{1,1},{0,1},{1,0}}).print();
	cout << endl; mlogical_and({{1,2},{3,4},{5,6}},{{1,0}}).print();
	cout << endl; mlogical_and({{1},{3},{5}},{{1,0,1}}).print();

	matrix mba = {{1,2},{3,4},{5,6}};
	matrix mba_f_i = find({{0,1},{1,0},{1,1}});
	cout << endl; mba.bulk_assign(mba_f_i, 9); mba.print();

	cout << endl; mba.assign_row(2,{1,2}); mba.print();
	
	/* // Testing access
	matrix m2(2,2);
	m2 += 4;
	cout << m2(1,1) << endl;
	m2.print();
	m2(1,3) = 5;
	m2.print();
	*/

	return 0;
}