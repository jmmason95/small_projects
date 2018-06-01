/*
The Arrhenius equation is used to fit a set of data:

k=A*exp(Ea/RT)

where k is the rate coefficient (in litres per mol per second, L mol-1 s-1),
R is the ideal gas constant (in kiloJoules per Mole per Kelvin, kJ mol-1 K-1), 
Ea is the activation energy (in kiloJoules, kJ), 
T is the temperature (in Kelvin, K),
and A is a constant (in litres per mol per second, L mol-1 s-1).

This is done through looped calculations of chi squared whilst shifting parameters Ea
and A between loops. This is done from multiple different starting points in order to
find a minimum chi squared value from a set of low values.

STARTING POINTS
Multiple starting points are chosen through a system which calculates a grid of 
points based on an inital Ea and A which are given by the identifiers INITIAL_Ea and
INITIAL_A. 

The number of starting points depends on GRID_LENGTH which on an axis of 
Ea against A defines the number of points are on each side of a square. At this time, 
this value must be odd and the centre of the square will be the inital values.

The initial values are calculated by adding/subtracting (1/GRID_LENGTH)*param*N to the 
parameter where N is 0, 1, -1, 2, -2 up tp GRID_LENGTH, and matching each Ea with 
each A value.

VARIABLES
vectors:	datak				used used to store hardcoded values for k
			dataT				used to store hardcoded values for T
doubles:	Ea 					stores current Ea value
			A					stores current A value
			chiSq				stores current chiSq value
			minEa				stores minimum Ea value so far
			minA				stores minimum A value so far
			minChiSq			stores minimum chiSq value so far
ints:		minIndex			stores the calculation index of the minimum chiSq value
			startingPointsN		stores number of starting points used
			counter				stores index of ongoing calculation

FUNCTIONS
The following functions are used:
CalculateChiSq(double stdev) -	Chi Squared is calculated using values standard deviation,
								the data vectors datak and dataT and the current values
								of Ea and A which are simply put into the chi squared 
								distribution function

Calculatek(double A,			Calculates k using the arhenius equation and the values 
double Ea, double T)			R, Ea, A and T

CalculateStandardDev			Returns the standard deviation of a vector list
(vector<double> v)				of values.

MinimiseParam(double stdevk,	This takes the standard deviation which is calulated
double *param)					at the beginning of the program and a pointer towards
								either &A or &Ea. The anderson sign indicates we
								pass a reference as the argument rather than the actual
								variable, so we can pass either A or Ea into the 
								function.
MAIN


*/


#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <iostream>

// Initial values of variables
#define INITIAL_Ea 150				// kJ mol-1
#define INITIAL_A 9E11				// L mol-1 s-1

// Running parameters
#define MAX_LOOP_RUNS 75
#define GRID_LENGTH 11	// set as odd number which is more than 1 (squared value of this is number of starting points)

// Constants
#define R 0.00831446				// kJ mol-1 K-1

using namespace std;

// Declaring methods
double CalculateChiSq(double stdev);
double Calculatek(double A, double Ea, double T);
double CalculateStandardDev(vector<double> v);
double MinimiseParam(double stdevk, double *param); // 1 for A, 2 for Ea
//double CalculateChiSqDiff(double stdevk);

// Variables accessible in all of main - better than passing through functions 

vector<double> dataT = { 700.f,730.f,760.f,790.f,810.f,840.f,910.f,1000.f };		// needs c++11 or later version for this declaration
vector<double> datak = { 0.011f,0.035f,0.105f,0.343f,0.789f,2.17f,20.f,145.f };

double minChiSq = 1001;	// minimum values of parameters and associated calculations:
double minEa;			// these variables hold the lowest shi and assosiated calculation number,
double minA;			// Ea and A.
int minIndex;

int startingPointsN = GRID_LENGTH * GRID_LENGTH; 

int counter = 0;
double chiSq = 0;
double Ea = INITIAL_Ea;
double A = INITIAL_A;

int main() {

	cout << "STANDARD DEVIATION\n";
	static double stdevk = CalculateStandardDev(datak);
	cout << "k:\t" << stdevk << "\n";

	int squareLength = GRID_LENGTH;
	double startingPointSteps = 1. / squareLength;
	cout << "step percentage of starting value added/subtracted to get starting points for minimisation: " << startingPointSteps << "\n\n";

	for (int i = 0; i < squareLength; i++)
	{
		for (int j = 0; j < squareLength; j++)
		{
			counter++;

			Ea = INITIAL_Ea + ((i - squareLength/2)*INITIAL_Ea*startingPointSteps);
			A = INITIAL_A + ((j - squareLength/2)*INITIAL_A*startingPointSteps);

			cout << "##########################- " << counter << " -##########################\n";
			cout << "STARTING POINT:\t" << Ea << "\t" << A << "\n";

			for (size_t k = 0; k < 5; k++)
			{
				cout << "----------------------CONVERGING Ea----------------------\n";
				MinimiseParam(stdevk, &Ea);

				cout << "----------------------CONVERGING A-----------------------\n";
				chiSq = MinimiseParam(stdevk, &A);
			}

			printf("RESULTS\nchi:\t%6.6f \nEa:\t%6.4f\tkJ mol-1 \nA:\t%6.3e\tL mol-1 s-1 \n",chiSq,Ea,A);

			if (chiSq < minChiSq) {
				minIndex = counter;
				minChiSq = chiSq;
				minEa = Ea;
				minA = A;
			}
		}
	}

	if (minChiSq > 1000) cout << "\nFAILURE TO FIND CONVERGANCE, CHECK INPUT PARAMETERS";

	else printf("\nMinimum chi found in calc no. %i:\nchi:\t%6.6f \nEa:\t%6.3f \nA:\t%6.3e \n", minIndex, minChiSq, minEa, minA);

	string b;
	cin >> b;
}

// methods (or functions)

/// calculates and returns chiSquared using vectors within scope and passing parameter stdevk. This prevents calculation of stdevk every time chiSquared is calculated 
double CalculateChiSq(double stdevk) {

	double chiSqrd = 0;

	for (size_t i = 0; i < datak.size(); i++)
	{
		chiSqrd += (1.f / pow(stdevk, 2))*(pow(datak[i] - Calculatek(A, Ea, dataT[i]), 2)); 
	}

	return chiSqrd;
}

/// use to calculate expected k values
double Calculatek(double A, double Ea, double T) {
	return A*exp(-Ea / (R*T));
}

/// calculate standard deviation of a list of numbers in a vector with type double
double CalculateStandardDev(vector<double> v) {

	double sum = 0;
	double mean = 0;
	double sumMeanDiffSqrd = 0;
	double variance;
	int N = v.size();

	for (int i = 0; i < v.size(); i++)
	{
		sum += v[i];
	}
	mean = sum / (int)v.size();

	for (int i = 0; i < v.size(); i++)
	{
		sumMeanDiffSqrd += pow((v[i] - mean), 2);
	}
	variance = sumMeanDiffSqrd / (v.size() - 1);

	return sqrt(variance);
}

/// return a minimised ChiSquared, minimising the specific parameter which is entered in the second argument as a reference for pointer *param
double MinimiseParam(double stdevk, double *param) { 

	int counter = 0;
	double chiSq1 = 0;				// values holding chiSquared 
	double chiSq2 = 0;
	double deltaMultiplier = 0.01;	// value holding the multiplier of delta param
	
	double deltaParam = deltaMultiplier * *param;	// initialise deltaParam as a percentage of the param (based on deltaMultiplier)

	printf("N  \tParam \t\tdChiSq \t\tdParam \t\tchiSq/8\n");

	while (counter < MAX_LOOP_RUNS)
	{
		chiSq1 = CalculateChiSq(stdevk);
		*param += deltaParam;
		chiSq2 = CalculateChiSq(stdevk);
		printf("%3i\t%6.4e\t%+7.6e\t%+7.4e\t%+6.4e\n", counter, *param, -(chiSq1 - chiSq2), deltaParam, chiSq2/8);
		
		if (chiSq2 - chiSq1 > 0) {
			*param -= deltaParam;					// return to original param value (because the new value gave worse chiSq)

			if (counter > 0) deltaParam *= 0.05f;	// if the first step increases chiSq, counter = 0, so the step doesn't drop to 10%. There may be many steps in the opposite direction to achieve minimisation.

			*param += deltaParam;
			chiSq1 = CalculateChiSq(stdevk);
			*param -= 2 * deltaParam;
			chiSq2 = CalculateChiSq(stdevk);
			*param += deltaParam;					// This finds chiSq for values just above and below 

			if (chiSq2 < chiSq1) deltaParam *= -1;
		}
		if (abs(chiSq2 - chiSq1) == 0.000) break;	// once this is true it has converged beyond further necessary calculation
		counter++;
	}

	return chiSq2;
}

/// Calculates differential of Chi Squares with respect to parameter which is chosen in the second argument
double CalculateChiSqDiff(double stdevk, int param) { // param: 0 for Ea, 1 for k
	double sum = 0;
	if (param == 0) {
		for (int i = 0; i < datak.size(); i++)
		{
			sum += dataT[i] - Ea - Ea * datak[i];
		}
	}
	else if (param == 1) {
		for (int i = 0; i < datak.size(); i++)
		{
			//sum += dataT[i] - k - k * datak[i];
		}
	}
	
	return ((-2 / pow(stdevk, 2))*sum);
}