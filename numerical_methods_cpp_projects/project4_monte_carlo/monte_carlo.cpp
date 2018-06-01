#include <cstdlib>
#include <cmath>
#include <iostream>
#include <ctime>

using namespace std;

// interpreter values related to the running of the script
#define LOWER_END 0.13
#define HIGHER_END 0.87
#define PI 3.1415926535
#define REPEATS 100
#define INTEGRAL 0.525835

// function declarations
double CrudeMonteCarloIntegration(int N);
double Function(double x);
double Crude_MC_Algorithm(int runs);

int main() {
	srand((unsigned)time(0)); // random numbers seeded by the current time

	// function called to evaluate crude Monte Carlo algorithm at N = 10, 100, and 1000
	double value10 = Crude_MC_Algorithm(10);
	double value100 = Crude_MC_Algorithm(100);
	double value1000 = Crude_MC_Algorithm(1000);

	// analytically evaluated integral stated at end of program
	cout << "\n" << INTEGRAL;

	// wait for user to end program
	int a;
	cin >> a;
	return 0;
}

double CrudeMonteCarloIntegration(int N) {
	double I_est = 0;
	for (int i = 0; i < N; i++)
	{
		// random number generated between lower and higher limits (0.13 - 0.87)
		double x = (HIGHER_END - LOWER_END) * (static_cast <float> (rand()) / (static_cast <float> (RAND_MAX)));
		x += LOWER_END;

		// cumulative f(x) values kept in variable I_est
		I_est += Function(x);
	}

	// divide by N and multiplied by (higher limit - lower limit) to find the average I_est
	I_est /= N;
	I_est *= (HIGHER_END - LOWER_END);
	return I_est;
}

/// Function which returns f(x) for supplied x where f(x) = 0.9*sin(PI*x)
double Function(double x) {
	return 0.9*sin(PI*x);
}

/// Crude Monte Carlo algorithm including 100 repeats, to achieve a better average 
double Crude_MC_Algorithm(int runs) {
	double summedIntegrals = 0;
	double summedErrors = 0;
	for (int i = 0; i < REPEATS; i++)
	{
		// Integrals and errors are cumulatively summed up in this for loop
		double calculatedIntegral = CrudeMonteCarloIntegration(runs);
		summedIntegrals += calculatedIntegral;
		summedErrors += abs(calculatedIntegral - INTEGRAL);
	}
	// the cumulative integrals and errors are divided by the number of repeats to attain averages
	double area = summedIntegrals / REPEATS;
	double error = summedErrors / REPEATS;

	cout << "\nFor this many runs: " << runs << ", ran 100 times, the area found is:";
	cout << "\n" << area;
	cout << " +/- " << error;
	return area;
}
