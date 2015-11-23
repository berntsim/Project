#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <algorithm>
#include <complex>
#include <limits>

#include "header.h"
#include "source.cpp"
using namespace std;
int main(){
//---------------------------------Parameters-----------------------------------
    int nbr_particles = 1000;         //Number of particles we want in the syst.
    int system_length = 100;         //Dimension of the array used for the syst.
    float r_p = 1.0;                //Radius of the particles.
    float L_min = 1.0;              //Minimum step corresponding to rand. walk.
    double pi = 3.1415926535897932384626433832; // Declaring pi.
    bool hit = false;
    int i_max = std::numeric_limits<int>::max();
    int counter = 0;
    int i = 0;
    int counter2 = 0;
    int len = system_length - 1;
    double D = 0.000001;
    double dt = 0.01;
    double t = 0;
    int s = 3;
    int iteration_length = 10;

    nbr_particles++;
    std::vector<std::vector<complex<double>>> clusters(nbr_particles,
                                  std::vector<complex<double>>(1));
    std::vector<std::complex<double>> particles;
    std::vector<std::vector<int>> num_grid_C(system_length,
                                           std::vector<int>(system_length));
    std::vector<std::vector<int>> num_grid_N(system_length,
                                           std::vector<int>(system_length));
//------------------------------Seeding MT RNG----------------------------------
    std::mt19937::result_type seed = time(0);
    //std::mt19937::result_type seed = 1447952689;
    writeSeed(seed);
    auto rand_dir = std::bind(std::uniform_real_distribution<float>(0,2*pi),
				   std::mt19937(seed));
    auto rand_seed = std::bind(std::uniform_int_distribution<int>(0,i_max),
				   std::mt19937(seed));
//-------------------------------Running prog.----------------------------------
    initGrid(system_length, rand_seed(), nbr_particles, r_p,
             clusters, num_grid_C, num_grid_N);
    std::ofstream out_stream;
    out_stream.open("cluster_size_distribution.txt");
    for (int i = 0; i < iteration_length; ++i){
        writeConfigColor(clusters, r_p, system_length, i);    
        plotConfig(system_length, std::to_string(i));
        out_stream << t << " " << clusterSize(clusters, s, system_length)
                   << std::endl;
        std::cout << std::endl;
        std::cout << "Iteration number " << i << std::endl;
        std::cout << " t = " << t << std::endl;
        std::cout << "number of clusters = " << clusters.size()-1 << std::endl;
        walkOnGrid(clusters, num_grid_C, num_grid_N, rand_seed(), L_min, pi, r_p,
                   counter, D, dt, nbr_particles);
        t += dt;
    }
    std::cout << "clusters.size() = " << clusters.size()-1 << std::endl;
    clustersTime(clusters, system_length);
    writeConfigColor(clusters, r_p, system_length, iteration_length);    
    plotConfig(system_length, std::to_string(iteration_length));
//-------testing section--------
//    std::cout << "dx(10) = " << findStepLength(D, dt, 10, nbr_particles) << std::endl;
    return 0;
}

