/*
This program generates the energies and intensities for rovibrational transitions in hydrogen chloride and deuterium
chloride. To do so, it uses the equations below. Rotational energy (E_rot) is:

(1)		E_rot = B_v*J (J + 1)

where J is the rotational energy level, and B_v is the rotation constant dependant on the vibrational quantum
number. To calculate B_v, we use equation

(2)		B_v = B_e - a_e (v + 0.5)

Where B_e is the hypothetical value of the rotation constant at the bottom of the potential well, a_e is the vibration-
rotation interaction constant, and v is the vibrational quantum number. Vibrational energy (E_vib) is:

(3)		E_vib = w_e (v + 0.5) - w_e_x_e (v + 0.5)^2

where v is the vibrational quantum number, w_e is the harmonic vibration wavenumber, and w_e_x_e is the anharmonic
vibration wavenumber. Total energy (E_tot) for a rovibrational level is therefore:

(4)		E_tot = E_vib + E_rot

Transition intensities are directly related to the population of the rotational energy level which is excited in a
transition. Populations are calculated relative to the ground level populations, and the values are sued arbritary units
for intensity:

(5)		N_J / N_0 = (2J + 1) * exp(-hcB_vJ(J+1)/(kT))

where N_J is a rotational level intensity, N_0 is the ground state intensity, h is Plancks constant, c is the speed of light,
k is the Boltzmann constant, and T is the temperature the population is measured at.

Hence, these equations are implemented here to gather the intensities and energies of different rovibrational levels for diatomic molecules
HCl and DCl. The method EvaluateDiatomicClSpectra(...) takes all the constants associated with a diatomic molecule, and calculates the required amount
of energy levels to then take differences to get the transitional energies of P and R branches. It uses a nested for loop to gather the energies and
store them in 2D array energyLevels[][]. It then uses a for loop to gather the intensities the transitions which are stored in arrays intesityR[] and
intensityP[]. In this same for loop, the energy level differences relating to transitional energies are also taken, and stored in arrays PositionR[]
and PositionP[].

Afterward, method PrintSpectra(...) is used to display the spectra in the console window. This prints the titles for the columns, then the actual energy
and intensity values, followed by a string of dashes which are proportional to intensity based on SPECTRA_RESOLUTION. It first prints the values for the
R branch, then the P branch.
*/

#include <math.h>
#include <iostream>

// Declaration of preprocessor identifiers/replacements
// enhances readability and allows value errors to be easily corrected

// PHYSICAL VALUES					   UNITS				CONSTANT
#define kB 1.38064852E-23			// J K -1				Boltzmann constant
#define PLANCKS_C 6.62607004E-34	// m 2 kg s -1			Plancks constant
#define SPEED_OF_LIGHT 2.99792E8	// m s -1				Speed of light

#define H_B_e 10.59341				// cm -1				Hypothetical rotation constant for hydrogen (H)
#define H_ALPHA_e 0.307				// cm -1				Vibration-rotation interaction constant for H
#define H_OMEGA_e 2990.946			// cm -1				Harmonic vibration number for H
#define H_OMEGA_e_x_e 52.8186		// cm -1				Anharmonic vibration number for H

#define D_B_e 5.44879				// cm -1				Hypothetical rotation constant for deuterium (D)
#define D_ALPHA_e 0.113291			// cm -1				Vibration-rotation interaction constant for D
#define D_OMEGA_e 2145.16			// cm -1				Harmonic vibration number for D
#define D_OMEGA_e_x_e 27.1825		// cm -1				Anharmonic vibration number for D

#define TEMPERATURE 298				// Kelvin				Temperature at which rovibrational spectrum is measured

// PROGRAM PARAMETERS
#define NUMBER_OF_V_LEVELS 2		// These two parameters specify how many rovibrational 
#define NUMBER_OF_J_LEVELS 11		// energy levels are calculated.
#define P_R_AMOUNT 10				// This specifys how many P and R branches energies are calculated
#define SPECTRA_RESOLUTION 16		// Specifies how many dashes (-) are represented by an intensity of 1 in the visual representation given for intensities

using namespace std;

// Declaration of methods
double CalculateTotalEnergy(double Evib, double Erot);
double CalculateVibEnergy(double v, double omega_e, double omega_e_x_e);
double CalculateRotEnergy(double B, double J, double alpha, double v);
double CalculateVibIntesity(double j, double B, double t, double v, double alpha);
double CalculateBvib(double Be, double alpha, double v);
double GetEnergy(double v, double omega_e, double omega_e_x_e, double alpha, double B, double J);
void PrintArray2(double a[NUMBER_OF_V_LEVELS][NUMBER_OF_J_LEVELS]);
void PrintSpectra(double positionP[], double intensityP[], double positionR[], double intensityR[], int resolution);
void EvaluateDiatomicClSpectra(double omega_e, double omega_e_x_e, double alpha, double b_e, double temperature);

// Declaration of global variables (access from anywhere in the script)
double positionR[P_R_AMOUNT];		// Holds transitional energies for R branch
double intensityR[P_R_AMOUNT];		// Holds intensities for R branch transitions

double positionP[P_R_AMOUNT];		// Holds transitional energies for P branch
double intensityP[P_R_AMOUNT];		// Holds intensities for P branch transitions

double energyLevels[NUMBER_OF_V_LEVELS][NUMBER_OF_J_LEVELS];	// Holds energies calculated for v = 1, 2, ... and J = 1, 2, 3 ... in 2D array

int main() {
	printf("__HCl spectra__\n\n");
	EvaluateDiatomicClSpectra(H_OMEGA_e, H_OMEGA_e_x_e, H_ALPHA_e, H_B_e, TEMPERATURE);

	printf("__DCl spectra__\n\n");
	EvaluateDiatomicClSpectra(D_OMEGA_e, D_OMEGA_e_x_e, D_ALPHA_e, D_B_e, TEMPERATURE);

	int a;
	cin >> a;
	return 0;
}

/// Implementation of equation 1
double CalculateRotEnergy(double B, double J, double alpha, double v) {
	return (CalculateBvib(B, alpha, v)*J*(J + 1));
}

/// Implementation of equation 2
double CalculateBvib(double Be, double alpha, double v) {
	return (Be - alpha*(v + 0.5));
}

/// Implementation of equation 3
double CalculateVibEnergy(double v, double omega_e, double omega_e_x_e) {
	return ((v + 0.5)*omega_e - (v + 0.5) * (v + 0.5) * omega_e_x_e);
}

/// Implementation of equation 4
double CalculateTotalEnergy(double Evib, double Erot) {
	return (Evib + Erot);
}

/// Implementation of equation 5
double CalculateVibIntesity(double J, double B, double T, double v, double alpha) {
	return ((2. * J + 1.)*exp((-J*(J + 1.)*(CalculateBvib(B, alpha, v)*1E2*PLANCKS_C*SPEED_OF_LIGHT)) / (kB*T)));
}

/// Function used to combine physical constants and return energy of rovibrational level
double GetEnergy(double v, double omega_e, double omega_e_x_e, double alpha, double B, double J) {
	return CalculateTotalEnergy(CalculateVibEnergy(v, omega_e, omega_e_x_e), CalculateRotEnergy(B, J, alpha, v));
}

/// Prints energies of different rovibrational levels from a 2-dimensional array. Used for testing/debugging.
void PrintArray2(double a[NUMBER_OF_V_LEVELS][NUMBER_OF_J_LEVELS]) {
	for (int v = 0; v < NUMBER_OF_V_LEVELS; v++)
	{
		for (int j = 0; j < NUMBER_OF_J_LEVELS; j++)
		{
			printf("v:%d j:%d   \t%f cm-1\n", v + 1, j + 1, a[v][j]);
		}
	}
	printf("\n");
}

/// Prints the spectra given the arrays of transitional ernegies and intensities
void PrintSpectra(double positionP[], double intensityP[], double positionR[], double intensityR[], int resolution) {

	printf("      Position /cm-1, Intensity\n");

	for (int i = (P_R_AMOUNT - 1); i >= 0; i--)
	{
		printf("R%-3d: %-14.2f, %-9.2f ", i + 1, positionR[i], intensityR[i]);
		for (int counts = 0; counts < (int)(intensityR[i] * (float)resolution); counts++)
		{
			printf("-");
		}
		printf("\n");
	}

	printf("\n");

	for (int i = 0; i < P_R_AMOUNT; i++)
	{

		printf("P%-3d: %-14.2f, %-9.2f ", i + 1, positionP[i], intensityP[i]);
		for (int j = 0; j < (int)(intensityP[i] * resolution); j++)
		{
			printf("-");
		}
		printf("\n");
	}
	printf("\n");
}

/// Given the constants for the diatomic molecule, this calculates energies and intensities of transitions, and stores
/// these in the global arrays
void EvaluateDiatomicClSpectra(double omega_e, double omega_e_x_e, double alpha, double b_e, double temperature) {
	for (int v = 0; v < NUMBER_OF_V_LEVELS; v++)
	{
		for (int j = 0; j < NUMBER_OF_J_LEVELS; j++)
		{
			energyLevels[v][j] = GetEnergy(v, omega_e, omega_e_x_e, alpha, b_e, j);
		}
	}

	for (int i = 0; i < P_R_AMOUNT; i++)
	{
		intensityP[i] = CalculateVibIntesity(i + 2, b_e, temperature, 0, alpha);
		intensityR[i] = CalculateVibIntesity(i + 1, b_e, temperature, 0, alpha);

		positionP[i] = -energyLevels[0][i + 1] + energyLevels[1][i];
		positionR[i] = -energyLevels[0][i] + energyLevels[1][i + 1];
	}

	PrintSpectra(positionP, intensityP, positionR, intensityR, SPECTRA_RESOLUTION);
}