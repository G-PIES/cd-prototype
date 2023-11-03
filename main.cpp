#include <iostream>
#include <cstring>
#include <cmath>
#include <array>

#include "cluster_dynamics.hpp"

// --------------------------------------------------------------------------------------------
//  GLOBALS
/*  Nuclear Reactor
    (1) Species Name
    (2) Flux (cm^2 / s)
    (3) Temperature (Kelvin)
    (4) Recombination in Cascades
    (5) Interstitials in Cascades
    (6) Vancacies in Cascades
*/
NuclearReactor OSIRIS =
{ 
    "OSIRIS", // (1)
    2.9e-7, // (2)
    330.f + CELCIUS_KELVIN_CONV,  // (3)
    .3f, // (4)
    // (5)
    .5f, // bi
    .2f, // tri
    .06f, // quad
    // (6)
    .06f, // bi
    .03f, // tri
    .02f  // quad
};

/*  Material
    (1) Species Name
    (2) Migration Energy (eV)
    (3) Diffusion Coefficients (cm^2 / s)
    (4) Formation Energy (eV)
    (5) Binding Energy (eV)
    (6) Recombination Radius (cm)
    (7) Interstitial Loop Bias
    (8) Interstitial Dislocation Bias
    (9) Vacancy Loop Bias
    (10) Vacancy Dislocation Bias
    (11) Initial Dislocation Density (cm^2)
    (12) Grain Size (cm)
*/
Material SA304 = { 
    "SA304", // (1)
    // (2)
    .45f, // i
    1.35f, // v
    // (3)
    1e-3, // i
    .6f, // v
    // (4)
    4.1f, // i
    1.7f, // v
    // (5)
    .6f,  // i
    .5f,  // v
    .7e-7, // (6)
    63.f, // (7)
    // (8)
    .8f, 
    1.1f, // param
    33, // (9)
    // (10)
    .65f, 
    1.f, // param
    10e10 * M_CM_CONV, // (11)
    4e-3 // (12)
};

int concentration_boundary;
int simulation_time;
int delta_time;

std::array<double, CONCENTRATION_BOUNDARY> interstitials;
std::array<double, CONCENTRATION_BOUNDARY> vacancies;
std::array<double, CONCENTRATION_BOUNDARY> interstitials_temp;
std::array<double, CONCENTRATION_BOUNDARY> vacancies_temp;

NuclearReactor* reactor = &OSIRIS;
Material* material = &SA304;
// --------------------------------------------------------------------------------------------


int main(int argc, char* argv[])
{
    concentration_boundary = CONCENTRATION_BOUNDARY;
    simulation_time = SIMULATION_TIME;
    delta_time = DELTA_TIME;

    interstitials.fill(0.f);
    vacancies.fill(0.f);
    interstitials_temp.fill(0.f);
    vacancies_temp.fill(0.f);

    // --------------------------------------------------------------------------------------------
    // main simulation loop
    for (int t = 0; t < simulation_time; t += delta_time)
    {
        #if VPRINT
        fprintf(stdout, "\n--------------------------------------------------------------------------------------- t = %d\n", t);
        #endif

        // calculate interstitial / vacancy concentrations for this time slice
        for (int n = 1; n < concentration_boundary; ++n)
        {
            interstitials_temp[n] = i_clusters(n);
            vacancies_temp[n] = v_clusters(n);
        }

        interstitials = interstitials_temp;
        vacancies = vacancies_temp;
    }
    // --------------------------------------------------------------------------------------------


    // --------------------------------------------------------------------------------------------
    // print results
    fprintf(stdout, "Cluster Size\t-\tInterstitials\t-\tVacancies\n\n");
    for (int n = 1; n < concentration_boundary; ++n)
    {
        fprintf(stdout, "%d\t\t\t%8.10f\t\t%8.10f\n\n", n, interstitials[n], vacancies[n]);
    }
    // --------------------------------------------------------------------------------------------


    return 0;
}