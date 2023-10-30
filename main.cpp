#include <iostream>
#include <stdio.h>	// printf and fprintf live here.
#include <stdlib.h> // For dynamic memory allocation
#include <cstring>  // For memset
#include <omp.h>	// OpenMP, for multicore parallelization.

#include "cluster_dynamics.hpp"

//#define CSV		// Useful if we want to change the output format from text
					// to pure CSV.


#ifndef CONCENTRATION_BOUNDS
#define CONCENTRATION_BOUNDS 100 
#endif

#ifndef NUMT
#define	NUMT	2
#endif;

// result arrays of interstitials and vacancies

// Dynamic memory allocation and command-line arguments can help
// optimize the simulation's memory usage. - Sean H.
//double interstitials[CONCENTRATION_BOUNDS];
//double vacancies[CONCENTRATION_BOUNDS];

int concentration_bounds;   // Setting an overall global variable to hold the concentration_bounds
                            // so that we aren't forces to pass the variable to every single function. - Sean H.

int num_threads;			// Set the default number of OpenMP threads to use in the for-loop section.

// number of clusters of N interstitials (in) per unit volume
/* Pokor et al. 2004, 2a
                  (1)     (2)                    (3)              (4)
    dCi(n) / dt = Gi(n) + a[i,n+1] * Ci(n + 1) - b[i,n] * Ci(n) + c[i,n-1] * Ci(n-1)
*/
// The interstitial and vacancies arrays that were allocated in the main() function
// can be passed by reference. In this case, we only need the interstitials. - Sean H.
double i_clusters(int in, double* inters, NuclearReactor& reactor, Material& material)
{
    // boundary
    if (in > concentration_bounds || in < 1) return 0;

    // if n + 1 has not yet been calculated
    if (inters[in + 1] < 0) {
        // Recurse down to the boundry base-case, see above.
        inters[in + 1] = i_clusters(in + 1, inters, reactor, material);
    }

    // Using the dynamically allocated array to store the interstitial values,
    // calculate the rates, and calculate the new quantity of interstitials. - Sean H.
    return
        // (1)
        reactor.i_defect_production(in) +
        // (2)
        iemission_vabsorption_np1(in + 1) * inters[in + 1] -
        // (3)
        iemission_vabsorption_n(in) * inters[in] +
        // (4)
        iemission_vabsorption_nm1(in - 1) * inters[in - 1];
    /*
    if (interstitials[in + 1] < 0)
    {
        // recurse to the boundary
        interstitials[in + 1] = i_clusters(in + 1, reactor, material);
    }
    */

   /*
    return
        // (1)
        reactor.i_defect_production(in) +
        // (2)
        iemission_vabsorption_np1(in + 1) * interstitials[in + 1] -
        // (3)
        iemission_vabsorption_n(in) * interstitials[in] +
        // (4)
        iemission_vabsorption_nm1(in - 1) * interstitials[in - 1];
    */
}

// number of clusters of N vacancies (vn) per unit volume
/* Pokor et al. 2004, 2a
                  (1)     (2)                    (3)              (4)
    dCv(n) / dt = Gv(n) + a[v,n+1] * Cv(n + 1) - b[v,n] * Cv(n) + c[v,n-1] * Cv(n-1)
*/
// The interstitial and vacancies arrays that were allocated in the main() function
// can be passed by reference. In this case, we only need the vacancies. - Sean H.
double v_clusters(int vn, double* vacans, NuclearReactor& reactor, Material& material)
{
    // boundary
    if (vn > concentration_bounds || vn < 1) return 0;

    // If n + 1 has not yet been calculated
    if (vacans[vn + 1] < 0) {
        // Recurse to the boundry and the base-case.
        vacans[vn + 1] = v_clusters(vn + 1, vacans, reactor, material);
    }

    return
        // (1)
        reactor.v_defect_production(vn) +
        // (2)
        iemission_vabsorption_np1(vn + 1) * vacans[vn + 1] -
        // (3)
        iemission_vabsorption_n(vn) * vacans[vn] +
        // (4)
        iemission_vabsorption_nm1(vn - 1) * vacans[vn - 1];


    // if n + 1 has not yet been calculated
    /*
    if (vacancies[vn + 1] < 0)
    {
        // recurse to the boundary
        vacancies[vn + 1] = v_clusters(vn + 1, reactor, material);
    }
    */

   /*
    return
        // (1)
        reactor.v_defect_production(vn) +
        // (2)
        iemission_vabsorption_np1(vn + 1) * vacancies[vn + 1] -
        // (3)
        iemission_vabsorption_n(vn) * vacancies[vn] +
        // (4)
        iemission_vabsorption_nm1(vn - 1) * vacancies[vn - 1];
    */

}

int main(int argc, char* argv[])
{
    // Pass the values of the CONCENTRATION_BOUNDS and NUMT via the command-line.
    // If no value was passed, we can simply rely on the previously
    // provided global CONCENTRATION_BOUNDS. If a value was provided, we can
    // update the global integer of concentration_bounds. - Sean H.
    if (argc >= 2) {
		num_threads = atoi(argv[1]);
    } else if (argv >= 3) {
        concentration_bounds = atoi(argv[2]);
    } else {
		num_threads = NUMT;
		concentration_bounds = CONCENTRATION_BOUNDS;
	}

#ifdef _OPENMP
		fprintf( stderr, "OpenMP is supported -- version = %d\n", _OPENMP );
#else
        fprintf( stderr, "No OpenMP support!\n" );
        return 1;
#endif

	omp_set_num_threads( num_threads );    // set the number of threads to use in parallelizing the for-loop:`

    // initialize result arrays ----------------------------------------

    // All of the values of the results arrays should be set to -1.f, except for the very
    // first elements, which must be set to 0.f
	// Equivalent to:
	// double* interstitials = (double*)malloc(concentration_bounds * sizeof(double));
    double* interstitials = (double*)calloc(concentration_bounds, sizeof(double));
	// double* vacancies = (double*)malloc(concentration_bounds * sizeof(double));
    double* vacancies = (double*)calloc(concentration_bounds, sizeof(double));

    // malloc() and calloc() return a value of NULL if the memory allocation failed. We need to
    // test for that.
    if (interstitials == NULL) {
        fprintf(stderr, "An error occurred when allocating memory for the interstitial array.\n");
        return 2;
    } else if (vacancies == NULL) {
        fprintf(stderr, "An error occurred when allocating memory for the vacancies array.\n");
        return 3;
    } else
        fprintf(stdout, "$d Bytes of memory was successfully allocated for both the interstitial and vacancy arrays.", concentration_bounds * sizeof(double));


    // If the arrays were successfully allocated memory, let the first characters of the interstitials
	// and vacancies arrays be set to 0.f.
    interstitials[0] = 0.f;
    vacancies[0] = 0.f;

    // Using memset to initialize the remaining values of each array to be -1.f
    // IMPORTANT: We need to offset the arrays by the size of the double datatype to
    // ensure we only overwrite the indices from 1 to concentration_bounds - 1.
    // As such, we are only setting (concentration_bounds - 1) * sizeof(double)
    // bytes of memory. - Sean H.
    memset(interstitials + sizeof(double), -1.f, (concentration_bounds - 1) * sizeof(double));
    memset(vacancies + sizeof(double), -1.f, (concentration_bounds - 1) * sizeof(double));

    /*
    for (int i = 1; i < concentration_bounds; ++i)
    {
        interstitials[i] = -1.f;
        vacancies[i] = -1.f;
    }
    */

    // -----------------------------------------------------------------

    // calculate interstitial / vacancy concentrations -----------------
    // Candidate for parallelization on the CPU through functional decomposition.
	// Need to research further to see if we can parallize these functions on a GPU. - Sean H.
	#pragma omp parallel for default(none) shared(concentration_bounds, interstitials, vacancies)
    for (int n = 1; n < concentration_bounds; ++n)
    {
        interstitials[n] = i_clusters(n, interstitials, OSIRIS, SA304);
        vacancies[n] = v_clusters(n, vacancies, OSIRIS, SA304);
    }
    // -----------------------------------------------------------------

    // print results ---------------------------------------------------
#ifdef	CSV
    for (int n = 1; n < concentration_bounds; ++n)
    {
        fprintf(stdout, "%d, %8.1f, %8.1f", n, interstitials[n], vacancies[n]);
    }
#else
    fprintf(stdout, "Cluster Size\t-\tInterstitials\t-\tVacancies\n\n");
    for (int n = 1; n < concentration_bounds; ++n)
    {
        fprintf(stdout, "%d\t\t\t%8.1f\t\t%8.1f\n\n", n, interstitials[n], vacancies[n]);
    }
#endif
    // -----------------------------------------------------------------

    // Free up the dynamically allocated arrays. THIS IS ESSENTIAL TO AVOID MEMORY LEAKS! - Sean H.
    free(interstitials);
    free(vacancies);

    return 0;
}