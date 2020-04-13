#pragma once
#include <iostream>
#include <cmath>
int NokSearch(int x, int y);
int NodSearch(int x, int y);
iusing namespace std;
class Rational {
public:
    Rational();
    Rational(int);
    Rational(int numerator, int denominator);
    int Numerator() const;

    int Denominator() const;

private:
    int p;
    int q;
};

Rational operator-(const Rational& lhs, const Rational& rhs);

bool operator==(const Rational& lhs, const Rational& rhs);
Rational operator*(const Rational& lhs, const Rational& rhs);

Rational operator/(const Rational& lhs, const Rational& rhs);
Rational operator+(const Rational& lhs, const Rational& rhs);
ostream& operator << (ostream& os, const Rational& r);
//Rational operator+=(const Rational& lhs, const Rational& rhs);
//Rational operator*=(const Rational& lhs, const Rational& rhs);
Rational operator*(const Rational& lhs, const int& val);
Rational operator/(const Rational& lhs, const int& val);

