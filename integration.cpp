#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES
#include "math.h"
#include <vector>
int NUMBER_OF_INTERVALS = 10000;
struct AddX
{
	float x = 0;
	AddX(float x):  x(x)
	{
	}
	float operator()(const float& y) const
	{
		return x + y;
	}
};

struct SinNX
{
	int n;
	SinNX(int n):n(n){}
	float operator()(const float& y) const
	{
		return sin(n*y);
	}
};
struct CosNX
{
	int n;
	CosNX(int n):n(n){}
	float operator()(const float& y) const
	{
		return cos(n*y);
	}
};

struct EAdPlusX
{
	float plus = 0;
	float mul = 1;
	EAdPlusX(float plus, float mul) : plus(plus), mul(mul){}

	float operator()(const float& x) const
	{
		return std::pow(M_E, (mul*x + plus));
	}

};

template <typename Function>
class MathematicalFunction {
public:
	MathematicalFunction(Function f, std::pair<float, float> d) : function(f), domain(d){ }
	float operator()(float x) const
	{
		if (x >= domain.first && x <= domain.second)
		{
			return function(x);
		}
		else
		{
			std::cout << "X not within bounds of the function";
			return 0.0;
		}

	}
	std::vector<float> operator()(float start, float end, int n) const
	{
		float interval = (end - start)/n;
		std::vector<float> ret;
		for(float x = start; x <= end; x += interval)
		{
			ret.push_back(function(x));
		}
		return ret;
	}

private:
	Function function;
	std::pair<float, float> domain;

};
template <typename Function>
float calculateIntegral(const float& start, const float& end,  const int& n, const Function& function)
{
	float interval = (end - start)/(n);
	float integral = 0;
	for (int i = 0; i < n; ++i)
	{
		float x_i = start + i * interval;
		float x_next = x_i + interval;
		integral += interval * (function(x_i) + function(x_next)) / 2;
	}
	return integral;
}

template <typename Function1, typename Function2>
auto operator*(const MathematicalFunction<Function1>& f1, const MathematicalFunction<Function2>& f2)
{
	auto newFunction = [f1, f2](float x) {
        return f1(x) * f2(x);
    };

    return MathematicalFunction<decltype(newFunction)>(newFunction, std::make_pair(0.0,10.0));
}
template <typename Function1, typename Function2>
auto operator+(const MathematicalFunction<Function1>& f1, const MathematicalFunction<Function2>& f2)
{
	auto newFunction = [f1, f2](float x) {
        return f1(x) + f2(x);
    };
    std::pair<float, float> newDomain = {
        std::max(f1.domain.first, f2.domain.first),
        std::min(f1.domain.second, f2.domain.second)
    };
    return MathematicalFunction<decltype(newFunction)>(newFunction, newDomain);
}


template <typename Function>
float calculateAZero(MathematicalFunction<Function> f)
{
	return (1/M_PI)*calculateIntegral(0,2*M_PI, NUMBER_OF_INTERVALS, f);
}
template <typename Function>
float calculateAN(MathematicalFunction<Function> f, int n)
{
	MathematicalFunction<CosNX> cosnx(CosNX(n), std::make_pair(0.0, 2*M_PI));
	return (1/(2*M_PI))*calculateIntegral(0,2*M_PI, NUMBER_OF_INTERVALS, f*cosnx);
}
template <typename Function>
float calculateBN(MathematicalFunction<Function> f, int n)
{
	MathematicalFunction<SinNX> sinnx(SinNX(n), std::make_pair(0.0, 2*M_PI));
	return (1/(2*M_PI))*calculateIntegral(0,2*M_PI, NUMBER_OF_INTERVALS, f*sinnx);
}

int main()
{
	MathematicalFunction<AddX> add(AddX(2.0), std::make_pair(0.0,100.0));
	MathematicalFunction<EAdPlusX> e_function(EAdPlusX(2.0,3.0), std::make_pair(0.0,100.0));
	MathematicalFunction<SinNX> sine(SinNX(1), std::make_pair(0.0, 10.0));
	
	auto f = add*sine;
	while(true)
	{
		int coef;
		std::cout << "Enter a number for the n-th Fourier coefficient" << std::endl;
		std::cin >> coef; 
		std::cout << calculateAN(add, coef) << std::endl;
		std::cout << calculateBN(add, coef) << std::endl;
	}
}



