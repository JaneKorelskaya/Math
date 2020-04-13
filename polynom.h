#include <iostream>
#include <vector>
#include <cmath>
#include "rational.h"
using namespace std;

template <typename T> 
T Pow(T x, int i) {
	T result = 1;
	for(int j = 0; j < i; ++j) 
		result = result * x;
	return result;
}

template <typename T>
class Polynom {
private:
	vector<T> odds;
public:
	Polynom() {
		
		odds.push_back(0);
	}
	Polynom(T a) {
		odds.push_back(a);
	}
	Polynom(T a, int n) {
		for (int i = 0; i < n + 1; ++i) {
			if (i == n) 
				odds.push_back(a);
			else 
				odds.push_back(0);
		} 
	}
	Polynom(vector<T>& val) {
		for (int i = 0; i < val.size(); ++i) 
			odds.push_back(val[i]); 
	}
	T Value(T x) const {
        	T result = static_cast<T>(0);
		//cout << result << endl;
        	for(int i = 0; i < odds.size(); ++i) {
			//T tmp = result;
                	result = result + odds.at(i) * Pow(x, i);
			//cout << result << endl;
		}
        	return result;
	}
	Polynom<T> Derivative() const {
		Polynom<T> result;
		if (odds.size() == 1)
			return result;
		result.odds.resize(odds.size() - 1);
		for (size_t i = 0; i < odds.size() - 1; ++i) 
			result.odds[i] = odds.at(i + 1) * (i + 1);
		return result;
	}
	
	Polynom<T> Integral() const {
		Polynom<T> result;
		result.odds.resize(odds.size() + 1);
		for(size_t i = 0; i < odds.size(); ++i) {
			result.odds[i + 1] = odds.at(i) / (i + 1);
		}
		return result;
	}
     	template<typename T4> friend Polynom<T4> operator+(const Polynom<T4>& p1, const Polynom<T4>& p2);
	template<typename T3> friend Polynom<T3> operator-(const Polynom<T3>& p1, const Polynom<T3>& p2);
        template<typename T2> friend Polynom<T2> operator*(const Polynom<T2>& p1, const Polynom<T2>& p2);
        template<typename T1> friend ostream& operator<<(ostream& os, const Polynom<T1>& p);
	friend ostream& operator<<(ostream& os, const Polynom<Rational>& p);
};

template<typename T>
Polynom<T> operator+(const Polynom<T>& p1, const Polynom<T>& p2) {
	Polynom<T> result;
	size_t size = 0;
	if(p1.odds.size() >= p2.odds.size())
		size = p1.odds.size();
	else 
		size = p2.odds.size();
	result.odds.resize(size);
	for (size_t i = 0; i < size; ++i) {
		if(i >= p1.odds.size()) {
			result.odds[i] = p2.odds.at(i);
		}
		else if (i >= p2.odds.size()) {
			result.odds[i] = p1.odds.at(i);
		}
		else 
			result.odds[i] = p1.odds.at(i) + p2.odds.at(i);
	}
	return result;
}

template<typename T>
Polynom<T> operator-(const Polynom<T>& p1, const Polynom<T>& p2) {
        Polynom<T> result;
        size_t size = 0;
        if(p1.odds.size() >= p2.odds.size())
                size = p1.odds.size();
        else
                size = p2.odds.size();
        result.odds.resize(size);
        for (size_t i = 0; i < size; ++i) {
                if(i >= p1.odds.size()) {
                        result.odds[i] = static_cast<T>(0) - p2.odds.at(i);
                }
                else if (i >= p2.odds.size()) {
                        result.odds[i] = p1.odds.at(i);
                }
                else
                        result.odds[i] = p1.odds.at(i) - p2.odds.at(i);
        }
        return result;
}

template<typename T>
ostream& operator<<(ostream& os, const Polynom<T>& p) {
	if (p.odds.size() == 1) {
		os << p.odds.at(0);
		return os;
	}
        for(size_t i = p.odds.size(); i > 0; --i) {
		if (p.odds.at(i - 1) == 0)
			continue;
                if (i == p.odds.size())
                        os << p.odds.at(i - 1) << "*x^" << i - 1;
		else if(i == 1) {
			if (p.odds.at(i - 1) > 0)
				os << " + " << p.odds.at(i - 1);
			else
				os << " - " << abs(p.odds.at(i - 1));
		}
                else if(i == 2) {
			if (p.odds.at(i - 1) > 0)
                        	os << " + " << p.odds.at(i - 1) << "*x";
			else 
				os << " - " << abs(p.odds.at(i - 1)) << "*x";
		}
                else {
			if (p.odds.at(i - 1) > 0)
                        	os << " + " << p.odds.at(i - 1) << "*x^" << i - 1;
			else 
				os << " - " << abs(p.odds.at(i - 1)) << "*x^" << i - 1;
		}
        }
        return os;
}
ostream& operator<<(ostream& os, const Polynom<Rational>& p) {
        if (p.odds.size() == 1) {
                os << p.odds.at(0);
                return os;
        }
        for(size_t i = p.odds.size(); i > 0; --i) {
                if (p.odds.at(i - 1).Numerator() == 0)
                        continue;
                if (i == p.odds.size())
                        os << p.odds.at(i - 1) << "*x^" << i - 1;
                else if(i == 1) {
                                os << " + " << p.odds.at(i - 1);
                }
                else if(i == 2) {
                                os << " + " << p.odds.at(i - 1) << "*x";
                }
                else {
                                os << " + " << p.odds.at(i - 1) << "*x^" << i - 1;
                }
        }
        return os;
}


template<typename T>
Polynom<T> operator*(const Polynom<T>& p1, const Polynom<T>& p2) {
	Polynom<T> result;
	result.odds.resize(p1.odds.size() + p2.odds.size() - 1);
	for (size_t i = 0; i < p1.odds.size(); ++i) {
		for (size_t j = 0; j < p2.odds.size(); ++j) {
			result.odds[i + j] = result.odds[i + j] +  p1.odds.at(i) * p2.odds.at(j);
		} 
	}
	return result;
}

	



