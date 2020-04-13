#pragma once
#include <iostream>
#include <exception>
#include <cmath> 
#include <cstdlib> 

using namespace std;

class Vector {
public:
	Vector();
	Vector(size_t);
	size_t Size() const;
	double Mod() const; 
    double& operator[](size_t);
    const double& operator[](size_t) const;
    
    friend ostream& operator<<(ostream&, const Vector&); 
    friend double operator*(const Vector&, const Vector&);
    friend Vector operator+(const Vector&, const Vector&);
    friend Vector operator-(const Vector&, const Vector&);
    friend istream& operator>>(istream& is, Vector &v); 
    friend Vector operator*(const Vector&, int);
    friend Vector operator*(int, const Vector&);
	~Vector();
private:
	double* vect;
	size_t size;
}; 



/////////////////////////////////////////////////////////////////
double& Vector::operator[](size_t index) {
    return vect[index];
}
const double& Vector::operator[](size_t index) const{
    return vect[index];
}
Vector::Vector() : vect(new double[0]), size(0) {};
Vector::Vector(size_t s) {
    size = s;
	vect = new double[size];
}

size_t Vector::Size() const{
	return size;
} 

double Vector::Mod() const{
	double sum = 0;
	for(size_t i = 0; i < size; ++i)
		sum += pow(vect[i], 2);
	return pow(sum, 0.5);
}

Vector::~Vector() {
    if (vect != nullptr)
	    delete[] vect;
}

ostream& operator<<(ostream& os, const Vector& v) {
    if (v.size != 0) {
	    for (size_t i = 0; i < v.size; ++i) 
		    os << v[i] << " ";
	    os << endl;
    }
	return os;
}

double operator*(const Vector& lhs, const Vector& rhs) {
	double res = 0;
	for (size_t i = 0; i < lhs.Size(); ++i) {
		res += lhs[i] * rhs[i];
	}
    return res; 
} 
Vector operator*(const Vector& v, int val) {
    Vector res(v.size);
    for(size_t i = 0; i < v.size; ++i) 
        res[i] = val * v[i];
    return res;
} 
Vector operator*(int val, const Vector& v) {
    Vector res(v.size);
    for(size_t i = 0; i < v.size; ++i) 
        res[i] = val * v[i];
    return res;
}
Vector operator+(const Vector& lhs, const Vector& rhs) {
	Vector res(lhs.size);
	for(size_t i = 0; i < lhs.size; ++i) {
		res[i] = lhs[i] + rhs[i];
	}
	return res;
} 
Vector operator-(const Vector& lhs, const Vector& rhs) {
	Vector res(lhs.size);
	for(size_t i = 0; i < lhs.size; ++i) {
		res[i] = lhs[i] - rhs[i];
	}
	return res;
}
istream& operator>>(istream& is, Vector &v) {
    if(v.size != 0) {
    cout << "Input the vector: " << endl;
	    for (size_t i = 0; i < v.Size(); ++i) 
		    is >> v[i];
    }
    else {
        cout << "Input a size of vector: " << endl;
        int s;
        is >> s;
        v.size = s;
        cout << "Input the vector: " << endl;
        for (size_t i = 0; i < v.size; ++i) 
		    is >> v.vect[i];
    }
	return is;	
}


////////////////////////////////////////////////////////////////////////////////
class Complex {
public:
	Complex();
	Complex(double, double);
	double Mod() const;

    friend Complex operator+(const Complex&, const Complex&);
    friend Complex operator-(const Complex&, const Complex&);
    friend Complex operator*(const Complex&, const Complex&); 
    friend Complex operator/(const Complex&, const Complex&);
    friend ostream& operator<<(ostream&, const Complex&);
    friend istream& operator>>(istream& is, Complex &c);
	~Complex();
private: 
	double im;
	double re;
};  
//////////////////////////////
Complex::Complex() : re(0.0), im(0.0) {};
Complex::Complex(double re_val, double im_val) : re(re_val), im(im_val) {}
double Complex::Mod() const {
	return pow(pow(re, 2) + pow(im, 2), 0.5);
} 
Complex::~Complex() {};

Complex operator+(const Complex& c1, const Complex& c2) {
	Complex res;
	res.re = c1.re + c2.re;
	res.im = c1.im + c2.im;
	return res;
}
Complex operator-(const Complex& c1, const Complex& c2) {
	Complex res;
	res.re = c1.re - c2.re;
	res.im = c1.im - c2.im;
	return res;
}
Complex operator*(const Complex& c1, const Complex& c2) {
	Complex res;
	res.re = c1.re * c2.re - c1.im * c2.im; 
	res.im = c1.im * c2.re + c1.re * c2.im;
    return res;
}

Complex operator/(const Complex& c1, const Complex& c2) {
    Complex res;
    res.re = (c1.re * c2.re + c1.im * c2.im) / pow(c2.Mod(), 2); 
    res.im = (c1.im * c2.re - c1.re * c2.im) / pow(c2.Mod(), 2);
    return res; 
}
ostream& operator<<(ostream& os, const Complex& c) {
	if (c.im >= 0) {
		os << c.re << " + " << c.im << " * i" << endl;
	}
	else {
		os << c.re << " " << c.im << " * i" << endl;
	}
	return os;
} 
istream& operator>>(istream& is, Complex& c) {
	cout << "Input re and im: " << endl;
	is >> c.re >> c.im;
	return is;
};


////////////////////////////////////////////////////////////
class Matrix {
public: 
    Matrix();
    Matrix(size_t);
    Matrix(const Matrix&);
    double Determinant() const;
    Matrix Transpose();
    Matrix Reverse();
    Matrix MinorMatr(size_t, size_t) const;
    double AlgAddition(size_t, size_t) const;
    void Resize(const size_t&);
    friend ostream& operator<<(ostream&, const Matrix&); 
    friend Matrix operator*(const Matrix&, const Matrix&);
    friend Matrix operator+(const Matrix&, const Matrix&); 
	friend Matrix operator*(double, const Matrix&);
    friend Matrix operator-(const Matrix&, const Matrix&);
    friend istream& operator>>(istream& is, Matrix &v); 
    ~Matrix(); 

	Matrix& operator = (const Matrix&);
	double Minor(const size_t& i, const size_t& j) const { 
		return this->MinorMatr(i, j).Determinant(); }
private: 
    size_t size;
	double** matr;
};

/////////////////////////////////////////////////////////////////////
double Matrix::Determinant() const {
		if (size == 1) {
			return (*this).matr[0][0];
		}
		double result = 0;
		for (size_t i = 0; i < size; ++i) {
			double minor = (*this).matr[0][i] * this->Minor(0, i);
			if (i % 2 == 1) {
				minor = minor * (-1);
			}
			result = result + minor;
		}
		return result;
	}

Matrix& Matrix::operator = (const Matrix& m) {
	if (this == &m) { return *this; }
	//size = m.size;
	delete[] matr;
    size = m.size;
    matr = new double*[size];
    for (size_t i = 0; i < size; ++i) {
	    (*this).matr[i] = new double[size];
    }
	for (size_t i = 0; i < size; ++i) {
		(*this).matr[i] = m.matr[i];
	}
	return *this;
}
Matrix::Matrix() : size(0), matr(nullptr) {};
void Matrix::Resize(const size_t& n) {
	delete[] matr;
	size = n;
	matr = new double*[size];
	for (size_t i = 0; i < size; ++i) {
		(*this).matr[i] = new double[size];
	}
}
Matrix::Matrix(size_t s) : size(s) {
	matr = new double*[size];  
	for(size_t i = 0; i < size; ++i) 
    	matr[i]= new double[size];
	for(int i = 0; i < size; ++i) {
		for(int j = 0; j < size; ++j)
			matr[i][j] = 0;
	}
}

Matrix::Matrix(const Matrix& m) {
    size = m.size;
    matr = new double*[size];
    for (size_t i = 0; i < size; ++i)
        matr[i] = new double[size];

    for(size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < size; ++j)
            matr[i][j] = m.matr[i][j];
    }
}
Matrix::~Matrix() {
    for (size_t i = 0; i < size; ++i) {
		    delete[] matr[i];
    }
	delete[] matr;
}
Matrix Matrix::Transpose() {
    Matrix res(size);
    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < size; ++j)
            res.matr[i][j] = matr[j][i];
    }
    return res;
}
Matrix Matrix::MinorMatr(size_t _i, size_t _j) const {
    Matrix result(size - 1);
	for (size_t i = 0; i < size - 1; ++i) {
		for (size_t j = 0; j < size - 1; ++j) {
			if (i < _i) {
				if (j < _j) {
					result.matr[i][j] = (*this).matr[i][j];
				}
				else {
					result.matr[i][j] = (*this).matr[i][j + 1];
				}
			}
			else {
				if (j < _j) {
					result.matr[i][j] = (*this).matr[i + 1][j];
				}
				else {
					result.matr[i][j] = (*this).matr[i + 1][j + 1];
				}
			}
		}
	}
	return result;
}

double Matrix::AlgAddition(size_t _i, size_t _j) const{
    return pow(-1, _i + _j) * (*this).Minor(_i, _j);
} 
Matrix Matrix::Reverse() {
		Matrix result(size); 
		if (this->Determinant() == 0) {
			cout << "Error: Determinant = 0";
			return Matrix(0);
		}
		for (size_t i = 0; i < size; ++i) {
			for (size_t j = 0; j < size; ++j) {
				result.matr[i][j] = this->AlgAddition(i, j) / this->Determinant();
			}
		}
		return result.Transpose();
}
///////////////////////////////////////////////////////////// 
Matrix operator*(double val, const Matrix& m) {
	Matrix res(m.size);
	for(int i = 0; i < m.size; ++i) {
		for(int j = 0; j < m.size; ++j)
			res.matr[i][j] = val * m.matr[i][j];
	}
	return res;
}
istream& operator>>(istream& is, Matrix &m) {
    if(m.size != 0) {
        cout << "Input the Matrix: " << endl;
	    for (size_t i = 0; i < m.size; ++i) {
	        for (size_t j = 0; j < m.size; ++j)
	            is >> m.matr[i][j];
	    } 
    }
    else {
        cout << "Input a size of Matrix: " << endl;
        int s;
        is >> s;
        //m.size = s;
        m.Resize(s);
        cout << "Input the Matrix: " << endl;
	    for (size_t i = 0; i < m.size; ++i) {
	        for (size_t j = 0; j < m.size; ++j)
	            is >> m.matr[i][j];
	    } 
    }
	return is;
}

ostream& operator<<(ostream& os, const Matrix& m) {
	for (size_t i = 0; i < m.size; ++i) {
		for (size_t j = 0; j < m.size; ++j) {
			os << m.matr[i][j] << " ";
		}
		os << endl;
	}
	return os;
}
Matrix operator+(const Matrix& m1, const Matrix& m2) {
	if (m1.size != m2.size) {
	    cout << "Different sizes of Matrix" << endl;
	}
		//throw invalid_argument("Different sizes of Matrix");
	Matrix res(m1.size);
		for (size_t i = 0; i < m1.size; ++i) {
			for (size_t j = 0; j < m1.size; ++j) {
				res.matr[i][j] = m1.matr[i][j] + m2.matr[i][j]; //can i do it???
		    }		
	    }
	return res;
}
Matrix operator-(const Matrix& m1, const Matrix& m2) {
	if (m1.size != m2.size) {
	    cout << "Different sizes of Matrix" << endl;
	    return Matrix(0);
	}
		//throw invalid_argument("Different sizes of Matrix");
	Matrix res(m1.size);
		for (size_t i = 0; i < m1.size; ++i) {
			for (size_t j = 0; j < m1.size; ++j) {
				res.matr[i][j] = m1.matr[i][j] - m2.matr[i][j]; 
		    }		
	    }
	return res;
}
Matrix operator*(const Matrix& m1, const Matrix& m2) {
	if (m1.size != m2.size) {
	    cout << "Different sizes of Matrix" << endl;
	}
	Matrix res(m1.size);

		for (size_t i = 0; i < m1.size; ++i) {
			for (size_t j = 0; j < m1.size; ++j) { 

				for (size_t m = 0; m < m1.size; ++m)
					res.matr[i][j] += m1.matr[i][m] * m2.matr[m][j]; 
		    }		
	    }
	return res;
}
