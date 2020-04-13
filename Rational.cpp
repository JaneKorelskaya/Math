#include "rational.h"
using namespace std;

int NodSearch(int x, int y) {
    int a = x;
    int b = y;
    while (a > 0 && b > 0) {
        if (a > b)
            a %= b;
        else
            b %= a;
    }
    return a + b;
}

int NokSearch(int x, int y) {
    return x * y / NodSearch(x, y);
}
Rational::Rational(int n) : p(n), q(1) {}
Rational::Rational() {
        p = 0;
        q = 1;
    }

Rational::Rational(int numerator, int denominator) {
        int counter = 0, p1 = 0, q1 = 0;
        if (denominator == 0)
                throw invalid_argument("Divided by zero");
        else if (numerator == 0) {
            p = numerator;
            q = 1;
        }
        else {
            if ((numerator > 0 && denominator < 0) || (numerator < 0 && denominator > 0))
                counter++;

            numerator = abs(numerator);
            denominator = abs(denominator);
            p1 = numerator;
            q1 = denominator;

            int nod = NodSearch(numerator, denominator);
            p = p1 / nod;
            q = q1 / nod;

            if (counter != 0)
                p = -p;
        }
    }

    int Rational::Numerator() const {
        return p;
    }

    int Rational::Denominator() const {
        return q;
    }

Rational operator+(const Rational& lhs, const Rational& rhs) {
    int x = lhs.Numerator() * rhs.Denominator() + rhs.Numerator() * lhs.Denominator();
    int y = lhs.Denominator() * rhs.Denominator();
    return Rational(x, y);
}

Rational operator-(const Rational& lhs, const Rational& rhs) {
    int x = lhs.Numerator() * rhs.Denominator() - rhs.Numerator() * lhs.Denominator();
    int y = lhs.Denominator() * rhs.Denominator();
    return Rational(x, y);
}

bool operator==(const Rational& lhs, const Rational& rhs) {
    if (lhs.Numerator() == rhs.Numerator() && lhs.Denominator() == rhs.Denominator())
        return true;
    else
        return false;
}
Rational operator*(const Rational& lhs, const Rational& rhs) {
    int y = lhs.Denominator() * rhs.Denominator();
    int x = lhs.Numerator() * rhs.Numerator();
    return Rational(x, y);
}

Rational operator/(const Rational& lhs, const Rational& rhs) {
    int y = lhs.Denominator() * rhs.Numerator();
    int x = lhs.Numerator() * rhs.Denominator();
    return Rational(x, y);
}

ostream& operator << (ostream& os, const Rational& r) {
        os << r.Numerator() << "/" << r.Denominator();
        return os;
}
//Rational operator+=(const Rational& lhs, const Rational& rhs) {
//	return lhs + rhs;
//}
Rational operator*(const Rational& lhs, const int& val) {
	return Rational(lhs.Numerator() * val, lhs.Denominator());
} 
Rational operator/(const Rational& lhs, const int& val) {
        return Rational(lhs.Numerator(), lhs.Denominator() * val);
}
//Rational operator*=(const Rational& lhs, const Rational& rhs) {
//	return lhs * rhs;
//}

