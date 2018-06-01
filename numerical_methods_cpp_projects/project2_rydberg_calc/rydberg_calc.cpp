#include <math.h>
#include <iostream>
#include <string>

// constant values
#define ELEMENTARY_Ch 1.602E-19		// C
#define PERM_OFS 8.854E-12			// F m-1
#define PLANCKS_C 6.626E-34			// J s
#define LIGHT_V 2.998E8				// m s-1
#define ELECTRON_M 9.109E-31		// kg
#define M_AMU 1.6605E-27			// kg
#define R_INFINITY 10973700.		// m-1

// method declarations
double Calculate_gj(int j, int s, int l);
double Calculate_R(double mass, int equation);
double Calculate_Reduced_Mass(double mass);
double Calculate_Reduced_Mass(double mass1, double mass2);
double Calculate_Rydberg(int p, int n);
double Calculate_Transition_Wavenumber(double rydberg, int n_2, int n_1, int p);

void Predict_Transition_Series(double rydberg, int p, int n_1, int i);

using namespace std;

int main() {
	bool quit = false;	// bool which is used to break from while loops
	char a = 0;			// char a and b used to check if user wants to break from while loops
	char b = 0;	
	int p = 0;			// p and n integers store proton and neutron numbers
	int n = 0;

	int j = 0;			// lande-g factor variables
	int s = 0;
	int l = 0;

	printf("POSITRONIUM\nCalculating positronium's three lowest energy lines using calculated Rydberg constant\n");
	double pRyd0 = Calculate_R(ELECTRON_M, 0);		// electron mass passed as the mass of the nucleus - reduced mass is calculated from two electrons for R
	Predict_Transition_Series(pRyd0, 1, 1, 3);

	while (true) {
		printf("Would you like to calculate Rydberg & Balmers series (option 1), or a Lande g-factor (option 2)?\n[1/2]: ");
		cin >> b;		// user gives option number. ASCII code 49 and 50 correspond to integer values 1 and 2
		if (b == 49) {
			while (true) {
				unsigned p = 0; // reset proton and neutron values to 0
				unsigned n = 0;
				
				printf("Please enter the number of protons and neutrons of the hydrogenic atom you wish to calculate Rydbergs Constant for\nPROTONS:  ");

				// taking user input for the number of protons and neutrons in the hydrogenic atom.
				// while loops are used to prevent none integer values being given for these variables
				while (!(cin >> p)) {
					cin.clear();			// cin bool value is cleared so it can be reevaluated 
					cin.ignore(100,'\n');	// cin is flushed (up to 100 characters) so that new input is fresh
					printf("PROTONS:  ");
				}

				printf("NEUTRONS: ");
				cin.clear();
				cin.ignore(1);


				while (!(cin >> n)) {
					cin.clear();
					cin.ignore(100, '\n');
					printf("NEUTRONS: ");
				}
				cin.clear();
				cin.ignore(1);

				// Rydberg and Balmer series calculations commence with functions.
				double rydberg = Calculate_Rydberg(p, n);
				Predict_Transition_Series(rydberg, p, 2, 5);

				// option to quit calculations of rydberg/balmer - if y or n are not inputed, the question is repeated.
				while (a != 121 && a != 110) {
					printf("quit this calculation? [y/n]: ");
					cin >> a;
					if (a == 121) quit = true;
					else if (a == 110) quit = false;
				}
				a = 0;				// a is reset so when the while loop comes back around, it asks again.
				if (quit) break;	// if user responds y, break from this calculation loop.
			}
		}
		else if (b == 50) {
			while(true){
				printf("\nPlease give integer values of j, s and l:\n");

				// takes user inputs for the values of j, s and l
				cout << "j:";
				cin >> j;
				cout << "s:";
				cin >> s;
				cout << "l:";
				cin >> l;

				// calculate lande-g factor using inputs passed to functions:
				double gj = Calculate_gj(j, s, l);
				printf("\ng_j in this case is: %.2f\n", gj);
				cout << "\n----------------------------\n\n";

				// option to quit calculation of lande g-factor 
				while (a != 121 && a != 110) {
					printf("quit this calculation? [y/n]: ");
					cin >> a;
					if (a == 121) quit = true;
					else if (a == 110) quit = false;
				}
				a = 0;
				if (quit) break;
			}
			
		}

		// option to quit program - if y or n are not inputed, the question is repeated.
		while (a != 121 && a != 110) {
				printf("\nquit program [y/n]: ");
				cin >> a;
				if (a == 121) quit = true;
				else if (a == 110) quit = false;
			}
		a = 0;
		if (quit) break;
	}
	return 0;
}

/// Rydberg constant calculated using this function. Mass is passed through regarding the proton and neutron nucleus,
/// and an int variable is used to choose which equation is used to calculate R_calc.
double Calculate_R(double mass, int equation) {
	double rm = Calculate_Reduced_Mass(mass);
	if (equation == 0) {
		double numerator0 = pow(ELEMENTARY_Ch, 4) * rm;
		double denominator0 = 8 * PERM_OFS*PERM_OFS*pow(PLANCKS_C, 3)*LIGHT_V;
		return (numerator0 / denominator0);
	}
	else if (equation == 1) {
		return (rm*R_INFINITY / ELECTRON_M);
	}
	else return 0;
}

/// Reduced mass calculated in this program requires only one mass as it always involves one electron, however extra functionality is left intact
double Calculate_Reduced_Mass(double mass1, double mass2) {
	return (mass1 * mass2 / (mass1 + mass2));
}

/// Reduced mass is calculated for hydrogenic atom when given nucleus mass in kg
double Calculate_Reduced_Mass(double mass) {
	return Calculate_Reduced_Mass(mass, ELECTRON_M);
}

/// Transition wavenumber is calculated for given n_2 and n_1 (energy levels), using Rydberg and Z numbers 
/// passed into the function
double Calculate_Transition_Wavenumber(double rydberg, int n_2, int n_1, int p) {
	double brackets = 1. / (n_2*n_2) - 1. / (n_1*n_1);
	return (-rydberg*p*p*brackets);
}

/// Rydberg constant is calculated for hydrogenic atom given the number of neutrons and protons. Neutron
/// count has very little effect on the outcome of the calculation as the reduced mass is dictated by the 
/// mass of the electron due to it being orders of magnitude lighter.
double Calculate_Rydberg(int p, int n) {
	double mass = ((p + n)*M_AMU);
	double answer1 = Calculate_R(mass, 0);
	double answer2 = Calculate_R(mass, 1);
	return answer1;
}

/// Transition series is calculated using a given rydberg constant, number of protons, final n value, and number of lines.
/// The energy lines are calculated up to a passed index i. These are printed in order into the console.
void Predict_Transition_Series(double rydberg, int p, int n_1, int i) {
	if (n_1 == 1) printf("\nUsing this Rydberg Constant %.1f cm-1, we find the Lymann series for the first %d lowest energy lines:\n\n", rydberg / 100, i);
	else if (n_1 == 2) printf("\nUsing this Rydberg Constant %.1f cm-1, we find the Balmer series for the first %d lowest energy lines:\n\n", rydberg / 100, i);
	else printf("\nUsing this Rydberg Constant %.1f cm-1, we find the following series for the first %d lowest energy lines:\n\n",rydberg / 100, i);
	
	for (int n = n_1 + 1; n < i + n_1 + 1; n++)
	{
		double wavenumber = Calculate_Transition_Wavenumber(rydberg, n, n_1, p);
		printf("n_%i ---> n_%i: %.1f nm\n", n, n_1, (1E9 / wavenumber));
	}
	cout << "\n----------------------------\n\n";
}

/// Lande g-factor is calculated given variables j, s and l.
double Calculate_gj(int j, int s, int l) {
	if (j == 0) {
		printf("\nj:0 results in unsolvable equation.\n");
		return 0;
	}
	double numerator = 3.*j*(j + 1) + s*(s + 1.) - l*(l + 1.);
	double denominator = 2.*j*(j + 1.);
	return numerator / denominator;
}
