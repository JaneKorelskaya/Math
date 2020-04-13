#include "polynom.h"
#include "rational.h"

using namespace std;

void TestPolynom() {
	cout << "_____________TEST CONSTRUCT_______________" << endl << endl;
	cout << "______________EMPTY POLYNOM_______________" << endl;
	Polynom<double> p1;
	cout << "0  ==  " << p1 << endl;
	cout << "_____________POLYNOM X^0__________________" << endl;
	Polynom<double> p2(-2.5);
	cout << "-2.5  ==  " << p2 << endl;
        cout << "_____________POLYNOM A*X^N________________" << endl;
        Polynom<double> p3(3.7, 5);
        cout << "3.7*x^5  ==  " << p3 << endl;
        cout << "______POLYNOM A*X^N + B*X^N-1 + ....______" << endl;
	vector<double> v = {3.4, -2.6, 1.2};
        Polynom<double> p4(v);
        cout << "1.2*x^2 - 2.6*x + 3.4  ==  " << p4 << endl;
        vector<double> v1 = {-4.1, 0, 2.2, 3.3};
        Polynom<double> p5(v1);
        cout << "3.3*x^3 + 2.2*x^2 - 4.1  ==  " << p5 << endl << endl;

	cout << "_____________TEST OPERATORS_______________" << endl << endl;
	cout << "0 + (- 2.5)  ==  " << p1 + p2 << endl;
        cout << "0 - (- 2.5)  ==  " << p1 - p2 << endl;
	cout << "3.7*x^5 + (1.2*x^2 - 2.6*x + 3.4)  ==  " << p3 + p4 << endl;
	cout << "3.3*x^3 + 2.2*x^2 - 4.1 - (1.2*x^2 - 2.6*x + 3.4)  ==  " 
		<< p5 - p4 << endl;
	cout << "(3.3*x^3 + 2.2*x^2 - 4.1) * (1.2*x^2 - 2.6*x + 3.4)  ==  " 
		<< p5 * p4 << endl << endl; 
        cout << "_____________TEST METHODS_______________" << endl << endl;
	cout << "Value (3.3*x^3 + 2.2*x^2 - 4.1) in x = 2  ==  " <<
		p5.Value(2.0) << endl;
	cout << "Integral (1.2*x^2 - 2.6*x + 3.4)  ==  " << 
		p4.Integral() << endl;
	cout << "Derivative (3.3*x^3 + 2.2*x^2 - 4.1)  ==  " << p5.Derivative() << endl;
	cout << "Derivative (3.7*x^5)  ==  " << p3.Derivative() << endl << endl;
	cout << "_____________TEST RATIONAL_______________" << endl << endl;	
	vector<Rational> v3 = {Rational(1, 3), Rational(2, 5), Rational(5, 9)};
        vector<Rational> v4 = {Rational(-7, 9), Rational(-11, 17), Rational(6, 10)};
	Polynom<Rational> p6(v3);
        Polynom<Rational> p7(v4);
	cout << "p1 : 5/9*x^2 + 2/5*x + 1/3  ==  " << p6 << endl;
	cout << "p2 : 6/10*x^2 - 11/17*x - 7/9  ==  " << p7 << endl;
	cout << "p1 + p2  ==  " << p6 + p7 << endl;
	cout << "p1 - p2  ==  " << p6 - p7 << endl;
	cout << "p1 * p2  ==  " << p6 * p7 << endl;
	cout << "p1 Value in x = 1/2  ==  " << p6.Value(Rational(1, 2)) << endl; 
	cout << "Integral p1  ==  " << p6.Integral() << endl;
	cout << "Derivative p2  ==  " << p7.Derivative() << endl;
		
}

int main() {
try{
	TestPolynom();
}
catch (invalid_argument& val) {
        cout << val.what() << endl; }
cout << "OK" << endl;
	return 0;
}
