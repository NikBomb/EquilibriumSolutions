#pragma once 

#include <vector>
#include <math.h>
#include <stdio.h>

class Gas {
protected:
	double Tb;
	double b;
	double Tc;
	double omega;
	double pc;

public:
	 double kEq(double p, double T) {
		/*Temperature R and pressure psia*/
		return 1/p * pow(10, a(p) + c(p) * F(T));
	}
	 double a(double p ) {
		 return 1.2 + 4.5E-4 * p + 15E-8 * p * p;
	 }

	 double c(double p) {
		 return 0.890 - 1.7E-4 * p - 3.5E-8 * p * p;
	 }

	 double F(double T) {
		 return b* (1 / Tb - 1 / T);
	 }

	 double getPc() {
		 return pc;
	 }

	 double getOmega() {
		 return omega;
	 }

	 double getTc() {
		 return Tc;
	 }
};

class C3 : public Gas {

public:
	C3() {
		b = 1888;
		Tb = 420;
	}

};

class nC4 : public Gas {

public:
	nC4() {
		b = 2293;
		Tb = 491;
	}
};

class nC5 : public Gas {

public:
	nC5() {
		b = 2750;
		Tb = 557;
	}

	
};

class C1 : public Gas {
public:
	C1() {
		b = 300;
		Tb = 96.8;
	}
};

class IsoC4 : public Gas {
public:
	IsoC4() {
		Tc = 274.46 + 459.67;
		pc = 527.9;
		omega = 0.1852;
	}
};