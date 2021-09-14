#include<iostream>
#include<algorithm>
#include<cmath>
#define USE MATH DEFINES
#include <chrono>
using namespace std::chrono;
auto start = high_resolution_clock::now();

using namespace std;
/// Cumulative distribution function
double normal_CDF(const double& x)
{
	return erfc(-x / sqrt(2.)) / 2;
}
/// Generation of number Pi : ) 
constexpr double pi = 3.14159265358979323846;

/// Probability density function
double normal_PDF(const double& x)
{
	return (1.0 / (pow(2 * pi, 0.5))) * exp(-0.5*x*x);
}

/// Calculate d1 separately
double d1_func(double& S,double& K,double& T,double& r,double& sigma)
{
	double d1 = (log(S / K) + (r + sigma * sigma / 2) * T) / (sigma * sqrt(T));
	return d1;
}
/// Calculate d2 separately
double d2_func(double& S, double& K, double& T, double& r, double& sigma)
{
	double d2 = (d1_func(S, K, T, r, sigma) - sigma * sqrt(T));
	return d2; 
}
/// Black Scholes for Call Option
double BS_call(double& S, double& K, double& T, double& r, double& sigma, double& d1,double& d2)
{
	double c = (S * normal_CDF(d1)) - K * exp(- r * T) * normal_CDF(d2);
	return c;
}
/// Black Scholes for Put Option
double BS_put(double& S, double& K, double& T, double& r, double& sigma, double& d1,double& d2)
{
	double p = K * exp(- r * (T))-S + (S * normal_CDF(d1)) - K * exp(-r * (T)) * normal_CDF(d2);
	///double p = K * exp(-r * (T)) - S + S * normal_CDF(d1) - K * exp(-1 * r * (T)) * normal_CDF(d2);
	return p;
}
/// Calculate Greeks for sensitivity - DELTA
double delta_call(double& S, double& K, double& T, double& r, double& sigma,double& d1)
{
	double delta_c = normal_CDF(d1);
	return delta_c;
}
/// Calculate Greeks for sensitivity - GAMMA
double gamma_call(double& S, double& K, double& T, double& r, double& sigma, double& d1)
{
	double gamma_c = (normal_PDF(d1)) / (S * sigma * sqrt(T));
	return gamma_c;
}
/// Calculate Greeks for sensitivity - THETA
double theta_call(double& S, double& K, double& T, double& r, double& sigma, double& d1, double& d2)
{
	double theta_c = ((-0.5 * S * normal_PDF(d1)) * sigma / sqrt (T) - r * K * exp(-r * T) * normal_CDF(d2))/365;
	return theta_c;
}
/// Calculate Greeks for sensitivity - VEGA
double vega_call(double& S, double& K, double& T, double& r, double& sigma, double& d1, double& d2)
{
	double vega_c = 0.01 * (S * normal_PDF(d1)) * sqrt(T);
	return vega_c;
}
/// Calculate Greeks for sensitivity - RHO 
double rho_call(double& S, double& K, double& T, double& r, double& sigma, double& d2)
{
	double rho_c = 0.01 * (K * (T) * exp(-r * T) * (normal_CDF(d2)));
	return rho_c;
}

int main(int argc, char ** argv) 
{
	double S = 18.75;
	double K = 29.00;
	double T = 0.547945;
	double r = 0.003;
	double sigma = 0.24;

	double d1 = d1_func(S, K, T, r, sigma);
	double d2 = d2_func(S, K, T, r, sigma);

	double c_bsm = BS_call(S, K, T, r, sigma, d1, d2);
	double p_bsm = BS_put(S, K, T, r, sigma, d1, d2);

	double delta_c = delta_call(S, K, T, r, sigma, d1);
	double gamma_c = gamma_call(S, K, T, r, sigma, d1);
	double theta_c = theta_call(S, K, T, r, sigma, d1, d2);
	double vega_c = vega_call(S, K, T, r, sigma, d1, d2);
	double rho_c = rho_call(S, K, T, r, sigma, d2);

	cout << "d1: " << d1 << endl;
	cout << "d2: " << d2 << endl << endl;

	cout << "Black Scholes - Call Option: " << c_bsm << endl;
	cout << "Black Scholes - Put Option: " << p_bsm << endl << endl;

	cout << "Delta call: " << delta_c << endl;
	cout << "Gamma call: " << gamma_c << endl;
	cout << "Theta call: " << theta_c << endl;
	cout << "Vega call: " << vega_c << endl;
	cout << "Rho call: " << rho_c << endl << endl;

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	cout << "Time taken by function: " << duration.count() << " microseconds" << endl;
	return 0;
}