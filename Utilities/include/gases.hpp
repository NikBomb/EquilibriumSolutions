#pragma once 

#include <vector>
#include <math.h>
#include <stdio.h>
#include <unordered_map>
#include <typeindex>
#include <typeinfo>
#include <utility>
#include <functional>


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

	 virtual double getID() = 0;
};

class C3 : public Gas {

public:
	C3() {
		b = 1888;
		Tb = 420;
		pc = 616;
		omega = 0.1522;
		Tc = 206.06 + 459.67;
	}

	double getID() override {
		return 1;
	};
};

class nC4 : public Gas {

public:
	nC4() {
		b = 2293;
		Tb = 491;
		Tc = 765.3;
		pc = 550.6;
		omega = 0.1995;
	}
	
	double getID() override {
		return 2; 
	};


};

class nC5 : public Gas {

public:
	nC5() {
		b = 2750;
		Tb = 557;
		pc = 488.6;
		omega = 0.2514;
		Tc = 385.8 + 459.67;
	}
	double getID() override {
		return 3;
	};
	
};

class C1 : public Gas {
public:
	C1() {
		b = 300;
		Tb = 96.8;
		Tc = 343.0;
		pc = 666.4;
		omega = 0.0104;
	}
	double getID() override {
		return 4;
	};
};

class IsoC4 : public Gas {
public:
	IsoC4() {
		Tc = 274.46 + 459.67;
		pc = 527.9;
		omega = 0.1852;
		b = 2037;
		Tb = 471;
	}
	double getID() override {
		return 5;
	};
};

class nC10 : public Gas {

public:
	nC10() {
		b = 3828;
		Tb = 805;
		Tc = 1111.7;
		pc = 305.2;
		omega = 0.4898;
	}

	double getID() override {
		return 6;
	};
};


struct GasMixtures {

	struct pair_hash
	{
		template <class T1, class T2>
		std::size_t operator() (const std::pair<T1, T2>& pair) const
		{
			return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
		}
	};
	typedef std::pair< std::type_index, std::type_index> pair;
	const std::unordered_map<pair, double, pair_hash> interactions{
		{ {std::type_index(typeid(C1)),std::type_index(typeid(nC4))}, 0.02 },
		{ {std::type_index(typeid(C1)),std::type_index(typeid(C3))}, 0.001 },
		{ {std::type_index(typeid(C1)),std::type_index(typeid(nC5))}, 0.02 },
		{ {std::type_index(typeid(nC5)),std::type_index(typeid(C1))}, 0.02 },
		{ {std::type_index(typeid(C3)),std::type_index(typeid(C1))}, 0.00 },
		{ {std::type_index(typeid(C3)),std::type_index(typeid(C3))}, 0.00 },
		{ {std::type_index(typeid(C3)),std::type_index(typeid(nC5))}, 0.01 },
		{ {std::type_index(typeid(nC5)),std::type_index(typeid(C3))}, 0.00 },
		{ {std::type_index(typeid(nC5)),std::type_index(typeid(nC5))}, 0.00 },
		{ {std::type_index(typeid(nC4)),std::type_index(typeid(C1))}, 0.02 },
		{ {std::type_index(typeid(C1)),std::type_index(typeid(nC10))}, 0.04 },
		{ {std::type_index(typeid(nC10)),std::type_index(typeid(C1))}, 0.04 },
		{ {std::type_index(typeid(nC4)),std::type_index(typeid(nC10))}, 0.00 },
		{ {std::type_index(typeid(nC10)),std::type_index(typeid(nC4))}, 0.00 },
		{ {std::type_index(typeid(C1)),std::type_index(typeid(C1))}, 0.00 },
		{ {std::type_index(typeid(nC4)),std::type_index(typeid(nC4))}, 0.00 },
		{ {std::type_index(typeid(nC10)),std::type_index(typeid(nC10))}, 0.00 },
		{ {std::type_index(typeid(IsoC4)),std::type_index(typeid(IsoC4))}, 0.00 },

	};

public:
	std::vector<Gas*> gases;
	std::vector<double> comp;
	
	
	GasMixtures(std::vector<Gas*>& g, std::vector<double>& c ) {
		
		gases = g;
		comp = c;
	};

	
	double getBinaryInteraction( Gas* g1, Gas* g2) {
		auto p = std::make_pair(std::type_index(typeid(*g1)), std::type_index(typeid(*g2)));
		return interactions.at(p);
	}
};