#ifndef MATRIX_H
#define MATRIX_H

#include <algorithm>
#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

struct range {
	int begin, end;
	range(int b, int e) : begin(b), end(e) {};
};

// A 2D matrix implemented using vectors
// Supports element-wise operations
class matrix {
private:
	// 2D vector, visually [row][col]
	vector<vector<double>> mat;

	// Number of rows and columns
	int nrows, ncols;

public:
	// Constructor
	matrix() {
		mat.resize(1);
		mat[0].resize(1, 0);
		nrows = ncols = 1;
	}
	matrix(const int& n, const int& m, double init_value=0) {
		mat.resize(n);
		for(int i = 0; i < n; ++i) {
			mat[i].resize(m, init_value);
		}
		nrows = n;
		ncols = m;
	}
	matrix(pair<int,int> s) : matrix(s.first, s.second) {}
	matrix(const initializer_list<initializer_list<double>>& l) {
		int size_x = 0;
		for(const auto it: l) {
			mat.push_back(it);
			if(it.size() > size_x) {
				size_x = it.size();
			}
		}
		// "Square out" the matrix if it is uneven
		for(size_t i = 0; i < mat.size(); ++i) {
			if(mat[i].size() < size_x) {
				mat[i].resize(size_x, 0);
			}
		}
		nrows = mat.size();
		ncols = size_x;
	}

	// Appends a matrix to the right side of this matrix
	// Number of rows must match
	matrix& append_right(const matrix& m) {
		// Check that number of rows are equal
		if(mat.size() != m.dim().first) {
			throw "Invalid dimensions";
		} else {
			int offset = ncols;
			int p = nrows;
			int q = m.dim().second;
			for(size_t i = 0; i < p; ++i) {
				mat[i].resize(q+offset);
				for(size_t j = 0; j < q; ++j) {
					mat[i][j+offset] = m(i,j);
				}
			}
			ncols = q+offset;
		}

		return *this;
	};
	// Appends a matrix to below this matrix
	// Number of columns must match
	matrix& append_down(const matrix& m) {
		// Check that number of columns are equal
		if(mat[0].size() != m.dim().second) {
			throw "Invalid dimensions";
		} else {
			int offset = nrows;
			int p = m.dim().first;
			int q = ncols;
			mat.resize(p+offset);
			for(size_t i = 0; i < p; ++i) {
				mat[i+offset].resize(q);
				for(size_t j = 0; j < q; ++j) {
					mat[i+offset][j] = m(i,j);
				}
			}
			ncols = p+offset;
		}

		return *this;
	};


	// Overload operator ()
	double& operator()(const int& i, const int& j) {
		if(i+1 > nrows) {
			mat.resize(i+1);
			nrows = i+1;
		}
		if(j+1 > mat[i].size()) {
			for(size_t u = 0; u < mat.size(); ++u) {
				mat[u].resize(j+1, 0);
			}
			ncols = j+1;
		}
		return mat.at(i).at(j);
	}
	// const ver.
	double operator()(const int& i, const int& j) const {
		return mat.at(i).at(j);
	}
	// Elements are accessed (0,0), (1,0), (2,0), ... (0,1), (1,1), (2,1), ... etc
	double& operator()(const int& n) {
		int i = n%nrows;
		int j = n/nrows;
		return mat.at(i).at(j);
	}
	// const ver.
	double operator()(const int& n) const {
		int i = n%nrows;
		int j = n/nrows;
		return mat.at(i).at(j);
	}
	// When called with a range, return a 1D slice of those indices
	// If the matrix being sliced is a column vector (ie, order nx1), then returned matrix is also a column oriented
	//	Otherwise, it is row oriented.
	// range is inclusive
	matrix operator()(const range& r) {
		matrix v;
		int idx = 0;
		// Column vector
		if(ncols == 1) {
			for(size_t i = r.begin; i <= r.end; ++i) {
				v(idx,0) = (*this)(i);
				idx++;
			}
		}
		// Return as row vector 
		else {
			for(size_t i = r.begin; i <= r.end; ++i) {
				v(0,idx) = (*this)(i);
				idx++;
			}
		}
		
		return v;
	}
	// When indexed with another matrix, the returned matrix is the same order as the indexing matrix
	// For output matrix c, indexing matrix b, and indexed matrix a (*this):
	//	c(i,j) = a(b(i,j))
	matrix operator()(const matrix& m) {
		matrix v(m.dim());
		for(size_t i = 0; i < m.dim().first; ++i) {
			for(size_t j = 0; j < m.dim().second; ++j) {
				v(i,j) = (*this)(m(i,j));
			}
		}
		
		return v;
	}

	// Assign a single double value to multiple indices at once
	// I is a row or column vector containing indices to assign s to
	void bulk_assign(const matrix& I, const double& s) {
		for(size_t i = 0; i < I.num_elements(); ++i) {
			(*this)(I(i)) = s;
		}
	}

	// Assignment to a specified row
	void assign_row(const unsigned int& row, const initializer_list<double>& l) {
		int idx = 0;
		for(const double it : l) {
			(*this)(row,idx) = it;
			idx++;
		}
	}

	// Assignment to a specified column
	void assign_col(const unsigned int& col, const initializer_list<double>& l) {
		throw "Not yet implemented";
	}

	
	// Overload operator +=
	// matrix + matrix
	matrix& operator+=(const matrix& b) {
		for(size_t i = 0; i < nrows; ++i) {
			for(size_t j = 0; j < ncols; ++j) {
				mat[i][j] += b.mat[i][j];
			}
		}
        return *this;
    }

	// Overload operator +=
	// matrix + double
	matrix& operator+=(double b) {
		for(size_t i = 0; i < nrows; ++i) {
			for(size_t j = 0; j < ncols; ++j) {
				mat[i][j] += b;
			}
		}
		return *this;
	}

	// Overload operator -=
	// matrix - matrix
	matrix& operator-=(const matrix& b) {
		for(size_t i = 0; i < nrows; ++i) {
			for(size_t j = 0; j < ncols; ++j) {
				mat[i][j] -= b.mat[i][j];
			}
		}
        return *this;
    }

	// Overload operator -=
	// matrix - double
	matrix& operator-=(double b) {
		for(size_t i = 0; i < nrows; ++i) {
			for(size_t j = 0; j < ncols; ++j) {
				mat[i][j] -= b;
			}
		}
		return *this;
	}

	// Overload operator *=
	// matrix * matrix
	// element-wise multiplication
	/*
	matrix& operator*=(const matrix& b) {
		for(size_t i = 0; i < mat.size(); ++i) {
			for(size_t j = 0; j < mat[i].size(); ++j) {
				mat[i][j] *= b.mat[i][j];
			}
		}
        return *this;
    }*/
	// Proper matrix multiplication
	matrix operator*(const matrix& b) {
		// Check dimensions
		if(nrows != b.dim().second || ncols != b.dim().first) {
			throw "Invalid dimensions";
		} else {
			// c(i,j) = a(i,0)*b(0,j) + a(i,1)*b(1,j) + ... + a(i,n)*b(m,j)
			matrix v(pair<int,int>(nrows,b.dim().second));

			for(size_t i = 0; i < v.dim().first; ++i) {
				for(size_t j = 0; j < v.dim().second; ++j) {
					for(size_t k = 0; k < ncols; ++k) {
						v(i,j) += mat[i][k] * b(k,j);
					}
				}
			}
			return v;
		}
    }
	// Overload operator *=
	// matrix * double
	matrix& operator*=(double b) {
		for(size_t i = 0; i < nrows; ++i) {
			for(size_t j = 0; j < ncols; ++j) {
				mat[i][j] *= b;
			}
		}
		return *this;
	}

	// Overload operator /=
	// matrix / matrix
	// element-wise right-division
	matrix& operator/=(const matrix& b) {
		for(size_t i = 0; i < nrows; ++i) {
			for(size_t j = 0; j < ncols; ++j) {
				mat[i][j] /= b.mat[i][j];
			}
		}
        return *this;
    }
	// matrix / double
	matrix operator/=(double b) {
		for(size_t i = 0; i < nrows; ++i) {
			for(size_t j = 0; j < ncols; ++j) {
				mat[i][j] /= b;
			}
		}
		return *this;
	}

	// Utility functions
	// The dimensions of the matrix, nxm
	pair<int,int> dim() {
		pair<int,int> s;
		s.first = nrows;
		s.second = ncols;
		return s;
	}
	pair<int,int> dim() const {
		pair<int,int> s;
		s.first = nrows;
		s.second = ncols;
		return s;
	}

	// Returns the largest value in the matrix
	double max() const {
		double largest_value = (*this)(0);
		for(size_t i = 0; i < nrows; ++i) {
			double it = *max_element(mat[i].begin(), this->mat[i].end());
			if(it > largest_value) {
				largest_value = it;
			}
		}

		return largest_value;
	}

	// Returns the total number of elements in the matrix, n*m
	int num_elements() {
		return nrows * ncols;
	}
	// const ver.
	int num_elements() const {
		return nrows * ncols;
	}

	void print() {
		for(size_t i = 0; i < mat.size(); ++i) {
			for(size_t j = 0; j < mat[i].size(); ++j) {
				cout << mat[i][j] << ' ';
			}
			cout << endl;
		}
	}

	// Returns the transpose of the matrix
	matrix T() {
		matrix v(ncols, nrows);
		for(size_t i = 0; i < nrows; ++i) {
			for(size_t j = 0; j < ncols; ++j) {
				v(j,i) = mat[i][j];
			}
		}
		return v;
	}
	
};

// ======================================================
// ================== END MATRIX CLASS ==================
// ======================================================

// Overload operator +
inline matrix operator+(matrix lhs, const matrix& rhs) {
	lhs += rhs;
	return lhs;
}
inline matrix operator+(matrix lhs, double s) {
	lhs += s;
	return lhs;
}
inline matrix operator+(double s, matrix rhs) {
	rhs += s;
	return rhs;
}

// Overload operator -
inline matrix operator-(matrix lhs, const matrix& rhs) {
	lhs -= rhs;
	return lhs;
}
inline matrix operator-(matrix lhs, double s) {
	lhs -= s;
	return lhs;
}
inline matrix operator-(double s, matrix rhs) {
	rhs *= -1;
	rhs += s;
	return rhs;
}



matrix pow(matrix lhs, double s) {
	pair<int,int> lhs_size = lhs.dim();
	for(size_t i = 0; i < lhs_size.first; ++i) {
		for(size_t j = 0; j < lhs_size.second; ++j) {
			lhs(i,j) = pow(lhs(i,j), s);
		}
	}
	return lhs;
}

// Overload operator /
inline matrix operator/(matrix lhs, const matrix& rhs) {
	lhs /= rhs;
	return lhs;
}
inline matrix operator/(matrix lhs, double s) {
	lhs /= s;
	return lhs;
}
inline matrix operator/(double s, matrix rhs) {
	rhs = pow(rhs,-1.0);
	rhs *= s;
	return rhs;
}

// Overload operator *
inline matrix operator*(matrix lhs, double s) {
	lhs *= s;
	return lhs;
}
inline matrix operator*(double s, matrix rhs) {
	rhs *= s;
	return rhs;
}


// Overload operator ==
// matrix == matrix
// element-wise equality
// Note, double comparisons may cause problems?
inline matrix operator==(const matrix& lhs, const matrix& rhs) {
	pair<int,int> lhs_size = lhs.dim();
	matrix result(lhs_size);

	for(size_t i = 0; i < lhs_size.first; ++i) {
		for(size_t j = 0; j < lhs_size.second; ++j) {
			result(i,j) = (lhs(i,j) == rhs(i,j));
		}
	}
	return result;
}
// matrix == double
// element-wise equality
inline matrix operator==(const matrix& lhs, const double& s) {
	pair<int,int> lhs_size = lhs.dim();
	matrix result(lhs_size);

	for(size_t i = 0; i < lhs_size.first; ++i) {
		for(size_t j = 0; j < lhs_size.second; ++j) {
			result(i,j) = (lhs(i,j) == s);
		}
	}
	return result;
}

// Overload operator <
// matrix < double
inline matrix operator<(const matrix& lhs, double s) {
	pair<int,int> lhs_size = lhs.dim();
	matrix result(lhs_size);

	for(size_t i = 0; i < lhs_size.first; ++i) {
		for(size_t j = 0; j < lhs_size.second; ++j) {
			result(i,j) = (lhs(i,j) < s);
		}
	}
	return result;
}
// Overload operator <=
// matrix <= double
inline matrix operator<=(const matrix& lhs, double s) {
	pair<int,int> lhs_size = lhs.dim();
	matrix result(lhs_size);

	for(size_t i = 0; i < lhs_size.first; ++i) {
		for(size_t j = 0; j < lhs_size.second; ++j) {
			result(i,j) = (lhs(i,j) <= s);
		}
	}
	return result;
}
// Overload operator >=
// matrix >= double
inline matrix operator>=(const matrix& lhs, double s) {
	pair<int,int> lhs_size = lhs.dim();
	matrix result(lhs_size);

	for(size_t i = 0; i < lhs_size.first; ++i) {
		for(size_t j = 0; j < lhs_size.second; ++j) {
			result(i,j) = (lhs(i,j) >= s);
		}
	}
	return result;
}

// Other functions
// Returns vector of indices of nonzero elements of a matrix,
// as a row vector if m is a row vector, as a column vector otherwise.
// Note: does double comparison cause problems?
matrix find(const matrix& m) {
	matrix v;
	int idx = 0;
	int total_elem = m.num_elements();
	bool is_row_vector = (m.dim().first == 1);

	for(size_t i = 0; i < total_elem; ++i) {
		if(m(i) != 0) {
			if(is_row_vector) {
				v(0,idx) = i;
			} else {
				v(idx,0) = i;
			}
			idx++;
		}
		
	}
	return v;
}

// Logical And on two matrices
// Return shape is that of the 2D matrix, otherwise if the inputs are row and col vectors respectively (or vice versa),
// then returns a rectangularized matrix
matrix mlogical_and(const matrix& m1, const matrix& m2) {
	// A points to the "taller" matrix
	const matrix *A, *B;
	if(m1.dim().first >= m2.dim().first) {
		A = &m1;
		B = &m2;
	} else {
		A = &m2;
		B = &m1;
	}
	matrix v(A->dim());
	// Take A by rows
	for(size_t ai = 0; ai < A->dim().first; ++ai) {
		// Take B by rows too if it matches with A
		if(B->dim().first == A->dim().first) {
			// Compare columns
			// Element-wise comparison if number of columns match (ie: A and B are the same shape)
			if(B->dim().second == A->dim().second) {
				for(size_t aj = 0; aj < A->dim().second; ++aj) {
					v(ai,aj) = ( (*A)(ai,aj) != 0 ) && ( (*B)(ai,aj) != 0 );
				}
			}
			// A is a column vector
			else if(A->dim().second == 1) {
				for(size_t bj = 0; bj < B->dim().second; ++bj) {
					v(ai,bj) = ( (*A)(ai,0) != 0 ) && ( (*B)(ai,bj) != 0 );
				}
			}
			// B is a column vector
			else if(B->dim().second == 1) {
				for(size_t aj = 0; aj < A->dim().second; ++aj) {
					v(ai,aj) = ( (*A)(ai,aj) != 0 ) && ( (*B)(ai,0) != 0 );
				}
			}
			else {
				throw "Invalid dimensions";
			}
		}
		// Otherwise compare against the single row of B (ie: B is a row vector)
		else if(B->dim().first == 1) {
			// Iterate through elements in row of A if there are any, also check that the number of columns match
			if( (A->dim().second > 1) && (A->dim().second == B->dim().second) ) {
				// Compare columns
				for(size_t aj = 0; aj < A->dim().second; ++aj) {
					v(ai,aj) = ( ( (*A)(ai,aj) != 0 ) && ( (*B)(0,aj) != 0 ) );
				}
			}
			// A is a column vector (ie: this is a row vector & column vector operation. Return a rectangled matrix)
			else if(A->dim().second == 1) {
				for(size_t bj = 0; bj < B->dim().second; ++bj) {
					v(ai,bj) = ( ( (*A)(ai,0) != 0 ) && ( (*B)(0,bj) != 0 ) );
				}
			}
			else {
				throw "Invalid dimensions";
			}
		}
		else {
			throw "Invalid dimensions";
		}
	}
	return v;
}

// Returns the largest value in the matrix
double max(const matrix& m) {
	return m.max();
}
// Max function on two matrices
// Return shape is that of the 2D matrix, otherwise if the inputs are row and col vectors respectively (or vice versa),
// then returns a rectangularized matrix
matrix max(const matrix& m1, const matrix& m2) {
	// A points to the "taller" matrix
	const matrix *A, *B;
	if(m1.dim().first >= m2.dim().first) {
		A = &m1;
		B = &m2;
	} else {
		A = &m2;
		B = &m1;
	}
	matrix v(A->dim());
	// Take A by rows
	for(size_t ai = 0; ai < A->dim().first; ++ai) {
		// Take B by rows too if it matches with A
		if(B->dim().first == A->dim().first) {
			// Compare columns
			// Element-wise comparison if number of columns match (ie: A and B are the same shape)
			if(B->dim().second == A->dim().second) {
				for(size_t aj = 0; aj < A->dim().second; ++aj) {
					throw "Not yet implemented";
				}
			}
			// A is a column vector
			else if(A->dim().second == 1) {
				for(size_t bj = 0; bj < B->dim().second; ++bj) {
					throw "Not yet implemented";
				}
			}
			// B is a column vector
			else if(B->dim().second == 1) {
				for(size_t aj = 0; aj < A->dim().second; ++aj) {
					throw "Not yet implemented";
				}
			}
			else {
				throw "Invalid dimensions";
			}
		}
		// Otherwise compare against the single row of B (ie: B is a row vector)
		else if(B->dim().first == 1) {
			// Iterate through elements in row of A if there are any, also check that the number of columns match
			if( (A->dim().second > 1) && (A->dim().second == B->dim().second) ) {
				// Compare columns
				for(size_t aj = 0; aj < A->dim().second; ++aj) {
					throw "Not yet implemented";
				}
			}
			// A is a column vector (ie: this is a row vector & column vector operation. Return a rectangled matrix)
			else if(A->dim().second == 1) {
				for(size_t bj = 0; bj < B->dim().second; ++bj) {
					v(ai,bj) = ( (*A)(ai) > (*B)(bj) ) ? (*A)(ai) : (*B)(bj);
				}
			}
			else {
				throw "Invalid dimensions";
			}
		}
		else {
			throw "Invalid dimensions";
		}
	}
	return v;
}

// Return a nxn matrix initialized with 1
matrix ones(const int& n) {
	return matrix(n,n,1);
}
// Return a nxm matrix initialized with 1
matrix ones(const int& n, const int& m) {
	return matrix(n,m,1);
}

// TODO:
// add dimension checks on matrix-matrix operations

#endif