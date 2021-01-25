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
	MatrixXd m(2,2);
	m(0,0) = 3;
	m(1,0) = 0;
	m(0,1) = -1;
	m(1,1) = m(0,0) + m(0,1);
	cout << m << endl << endl;
	cout << m(2) << endl << endl;
	auto mfind = find(m);
	cout << mfind << endl << endl;
	cout << index(m,mfind) << endl << endl;

	MatrixXi n = MatrixXi::Random(4,4);
	cout << n << endl << endl;
	cout << n(all,seq(1,3)) << endl << endl;

	MatrixXd a3(3,3);
	a3 << 1,2,3,
		  4,5,6,
		  7,8,9;
	MatrixXd b3(3,3);
	b3 << 1,1,9,
		  4,5,5,
		  7,8,9;
	MatrixXd a53 = (b3.array() == 5).cast<double>();
	cout << a53 << endl << endl;
	MatrixXd ab3 = (a3.array() == b3.array()).cast<double>();
	cout << ab3 << endl << endl;
	MatrixXd a54ab3 = (a53.array() != 0 && ab3.array() != 0).cast<double>();
	cout << a54ab3 << endl << endl;
	MatrixXd maxa3b3 = ( a3.array() > b3.array() ).select(a3,b3);
	cout << maxa3b3 << endl << endl;

	MatrixXd matA(2, 2);
	matA << 1, 2, 3, 4;
	cout << matA << endl << endl;
	// append right
	matA = append_right(matA,MatrixXd::Ones(matA.rows(),2));
	cout << matA << endl << endl;

	// append down
	matA = append_down(matA,MatrixXd::Ones(2,matA.cols()).array() + 3);
	cout << matA << endl << endl;

	/*a3 << 0,2,3,
		  0,0,6,
		  0,0,0;
	b3 << 1,2,3,
		  4,5,6,
		  7,8,9;*/
	MatrixXd andab = (a3.array() != 0 && b3.array() != 0).cast<double>();
	cout << "andab" << endl << andab << endl << endl;

	ArrayXi fandab = find(andab);
	cout << fandab << endl << endl;

	cout << assign_d_by_index(a3,fandab,9) << endl << endl;

	MatrixXd l6(1,6);
	l6 << 5,6,7,8,9,0;
	cout << l6 << endl << endl;

	int N = 5;
	MatrixXd r(3,N);
	r << vseq(0,N-1), MatrixXd::Zero(2,N);
	cout << r << endl << endl;

	cout << "madd" << endl;
	MatrixXd ma = MatrixXd::Ones(2,2);
	MatrixXd madd = ma.array() + 3;
	cout << madd << endl << endl;
	MatrixXd app(2,1);
	app << 5,9;
	madd = append_right(madd,app);
	cout << madd << endl << endl;
	MatrixXd madd_upsidedown(madd.rows(),madd.cols());
	madd_upsidedown << madd(1,all), madd(0,all).array()+2;
	cout << madd_upsidedown << endl << endl;
	//madd_upsidedown(seq(0,1),all) = madd_upsidedown(seq(0,1),all).array() + MatrixXd::Random(2,madd_upsidedown.cols()).array()/2*1e-8;

	MatrixXd cv5(5,1);
	VectorXd seqcv5 = cv5(seq(0,3),0);
	cout << cv5 << endl << cv5.rows() << ", " << cv5.cols() << endl << endl;

	MatrixXd rstoo(2,3);
	rstoo << 1,3,2, 9,8,7;
	cout << rstoo << endl << endl;
	MatrixXd rstoo_d = rstoo.array().pow(2).colwise().sum().sqrt();
	cout << rstoo_d << endl << endl;
	//cout << index(rstoo(0,all),rstoo_d) << endl << endl;
	//MatrixXd d = madd_upsidedown.array().pow(2).colwise().sum().sqrt();
	//cout << d << endl << endl;
	MatrixXd unrstoo = rstoo.unaryExpr([](double d) {
		return d * 2;
	});
	cout << unrstoo << endl << endl;
	//MatrixXd nullrstoo = MatrixXd::NullaryExpr([&]() { return 5; });
	MatrixXd nullrstoo = MatrixXd::NullaryExpr(rstoo.rows(), rstoo.cols(),
		[&rstoo](Index i) { return rstoo(i) - i; });

	cout << nullrstoo << endl << endl;

	MatrixXi rstoo_d_I = sort(rstoo_d);
	cout << rstoo_d_I << endl << endl;
	for(int i = 0; i < rstoo.rows(); ++i) {
		rstoo(i,all) = index(rstoo(i,all), rstoo_d_I);
	}
	cout << rstoo << endl << endl;

	for(auto row : rstoo.rowwise()) {
		row(1) = 1;
		row(2) = 2;
	}
	cout << rstoo<< endl << endl;
	cout << rstoo(all,0).mean() << endl << endl;

	cout << (rstoo.array().pow(2) + rstoo.array().pow(2)).sqrt() << endl << endl;
	cout << rstoo(seq(0,1),seq(0,1)) + rstoo(seq(0,1),seq(0,1)).transpose() << endl <<  endl;

	/*MatrixXd mx = MatrixXd::Ones(3,3) *1;
	MatrixXd my = MatrixXd::Ones(3,3) *2;
	MatrixXd mz = MatrixXd::Ones(3,3) *3;
	
	vector<double> vi(3*3*3);
	iota(vi.begin(),vi.end(),1);
	Eigen::TensorMap<Eigen::Tensor<double,3>> txyz(vi.data(), 3, 3, 3 );
	//Tensor<double,3> txyz(3,3,3);
	//txyz.setZero();
	cout << txyz << endl << endl;

	
	Eigen::array<int, 3> startIdx = {0,0,1};	// start a 0,0,0
	Eigen::array<int, 3> extent = {3,3,1};		// take 3x3x1 slice starting at (0,0,0)
	cout << txyz.slice(startIdx,extent) << endl << endl;
	cout << txyz.chip(1,2) << endl << endl;*/

	MatrixXd mdiv = MatrixXd::Ones(3,3) * 6;
	MatrixXd mdivisor = MatrixXd::Ones(3,3) * 2;
	mdivisor(1,1) = 0;
	MatrixXd mquot = mdiv.array()/mdivisor.array();
	cout << mquot << endl << endl;
	MatrixXd mquot_i = mquot.array().isInf().cast<double>();
	cout << mquot_i << endl << endl;
	//assign_d_by_index(mquot, find(mquot_i), 0);
	mquot = (mquot.array().isInf()).select(0, mquot);	// the easier way
	cout << mquot << endl << endl;


	/*cout << "andab" << endl << andab << endl << endl;
	SparseMatrix<double> sm_andab = andab.sparseView();
	cout << sm_andab << endl << endl;
	cout << sm_andab(seqN(0,2),all) << endl << endl;*/

	DiagonalMatrix<double, Dynamic> dm(5);
	dm.diagonal() = VectorXd::Ones(5);
	//MatrixXd dm_xd = dm;
	cout << (MatrixXd)dm << endl << endl;
	cout << MatrixXd::Ones(5,5)*5 - (MatrixXd)dm << endl << endl;

	MatrixXd a1 = MatrixXd::Ones(3,3)*2;
	MatrixXd a2 = MatrixXd::Ones(3,3)*4;
	cout << a1.array()*a2.array() << endl << endl;

	MatrixXd cv1 = MatrixXd::Ones(5,1)*3;
	MatrixXd rv1 = MatrixXd::Ones(5,1)*2;
	cout << (cv1.array()+1).matrix() * rv1.transpose() << endl << endl;

	MatrixXd mbase(5,5);
	mbase << 1,2,3,4,5, 6,7,8,9,10, 11,12,13,14,15, 16,17,18,19,20, 21,22,23,24,25;
	cout << mbase << endl << endl;
	// variation 1: using full 2D array as mask
	MatrixXd mask = (mbase.array() > 12).cast<double>();
	cout << mask << endl << endl; // note: didn't use find!
	MatrixXd mbasemasked = MatrixXd::NullaryExpr(mask.rows(), mask.cols(), [&](Index i) { 
		// apply function if Mask is 1
		if(mask(i)) {
			return myfunction(mbase(i), 1);
		}
		// Mask is 0
		else {
			return mbase(i);
		}
	});
	cout << mbasemasked << endl << endl;
	// variation 2: using find(), i.e., a vector of indices to access/modify
	MatrixXi mmf = find(mask);
	MatrixXd mbase2 = mbase;
	mmf = mmf.unaryExpr([&](int x) {
		mbase2(x) = myfunction(mbase2(x), 1);
		return x;
	});
	cout << mbase2 << endl << endl;
	
	MatrixXd mcol(1,5);
	mcol << 1,2,3,4,5;
	MatrixXi mfcol = find( (mcol.array() > 3 || mcol.array() == 2).cast<double>() );
	cout << mfcol << endl << endl;
	cout << mbase2 << endl << endl;
	// reshaped reshapes matrix into a col vector
	mbase2(mfcol.reshaped(),1) = mbase2(mfcol.reshaped(),1).array() * 10;
	cout << mbase2 << endl << endl;

	MatrixXd matan2 = MatrixXd::NullaryExpr(mbase.rows(), mbase.cols(), [&](Index i) { 
		return atan2(mbase(i),1);
	});
	cout << matan2 << endl << endl;
	cout << exp(mcol.array()) << endl << endl;

	MatrixXd meven = mbase.unaryExpr([&](double x) {
		if( (int)x % 2 == 0) {
			return 1.;
		} else {
			return 0.;
		}
	});
	cout << meven << endl << endl;
	MatrixXi mevenmask = find(meven);
	cout << mevenmask << endl << endl;
	//cout << MatrixXd::NullaryExpr([&](meven))
	cout << sqrt( pow(mevenmask.unaryExpr([&](int x) {return mbase(x);}).mean(),2) * mevenmask.size() ) / 2.0 << endl << endl;
	MatrixXd massign = vseq(1, mevenmask.rows());
	MatrixXd mdest = MatrixXd::Zero(mbase.rows(), mbase.cols());
	mevenmask = MatrixXi::NullaryExpr(mevenmask.rows(), mevenmask.cols(), [&](Index i) {
		mdest(mevenmask(i)) = massign(i);
		return mevenmask(i);
	});
	cout << mdest << endl << endl;
	// max of first 12 Index position coefficients of mbase
	cout << MatrixXd::NullaryExpr(12, 1, [&](Index i) {
		return mbase(i);
	}).maxCoeff() << endl;

	//sph_sim sph;

}
