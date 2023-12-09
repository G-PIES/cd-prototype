#include <iostream>
#include <cstring>
#include <cmath>
#include <array>


#include "material.hpp"
#include "cluster_dynamics.hpp"
#include "nuclear_reactor.hpp"

#ifndef VPRINT
#define VPRINT false
#endif

#ifndef VBREAK
#define VBREAK false
#endif

#ifndef CSV
#define CSV false
#endif

double concentration_boundary;
double simulation_time;
double delta_time;
double sample_interval; // How often (in seconds) to record the state

void print_start_message()
{
    fprintf(stderr, "\nSimulation Started: ");
    fprintf(stderr, "delta_time: %g, ", delta_time);
    fprintf(stderr, "simulation_time: %g, ", simulation_time);
    fprintf(stderr, "concentration_boundary: %llu\n", (unsigned long long)concentration_boundary);
}

void print_state(ClusterDynamicsState& state)
{
    if (!state.valid) fprintf(stdout, "\nINVALID SIM @ Time=%g", state.time);
    else fprintf(stdout, "\nTime=%g", state.time);

    if (state.interstitials.size() != concentration_boundary || state.vacancies.size() != concentration_boundary)
    {
        fprintf(stderr, "\nError: Output data is incorrect size.\n");
        return;
    }

    fprintf(stdout, "\nCluster Size\t\t-\t\tInterstitials\t\t-\t\tVacancies\n\n");
    for (uint64_t n = 1; n < concentration_boundary; ++n)
    {
        fprintf(stdout, "%llu\t\t\t\t\t%13g\t\t\t  %15g\n\n", (unsigned long long)n, state.interstitials[n], state.vacancies[n]);
    }
}

void print_csv(ClusterDynamicsState& state)
{
    for (uint64_t n = 1; n < concentration_boundary; ++n)
    {
        fprintf(stdout, "%g,%llu,%g,%g\n", state.time, (unsigned long long) n, state.interstitials[n], state.vacancies[n]);
    }
}

int main(int argc, char* argv[])
{
    NuclearReactor reactor = nuclear_reactors::OSIRIS();
    Material material = materials::SA304();

    // Default values
    concentration_boundary = 10;
    simulation_time = 1.0;
    delta_time = 1e-5;
    sample_interval = delta_time;

    // Override default values with CLI arguments
    switch (argc)
    {
        case 4:
            concentration_boundary = strtod(argv[3], NULL);
        case 3:
            simulation_time = strtod(argv[2], NULL);
        case 2:
            delta_time = strtod(argv[1], NULL);
        default:
            break;
    }

    ClusterDynamics cd(concentration_boundary, reactor, material);

    print_start_message();

    #if CSV
    fprintf(stdout, "Time (s),Cluster Size,Interstitials / cm^3,Vacancies / cm^3\n");
    #endif

    ClusterDynamicsState state;

    // --------------------------------------------------------------------------------------------
    // main simulation loop
    for (double t = 0; t < simulation_time; t = state.time)
    {
        // run simulation for this time slice
        state = cd.run(delta_time, sample_interval);

        #if VPRINT 
            print_state(state);
        #elif CSV
            print_csv(state);
        #endif

        if (!state.valid) 
        {
            break;
        }

        #if VBREAK
        fgetc(stdin);
        #endif
    }
    // --------------------------------------------------------------------------------------------


    // --------------------------------------------------------------------------------------------
    // print results
    #if !VPRINT && !CSV
    print_state(state);
    #endif
    // --------------------------------------------------------------------------------------------

    return 0;
}
