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
    int nbr_particles = 15;         //Number of particles we want in the syst.
    int system_length = 10;         //Dimension of the array used for the syst.
    float r_p = 1.0;                //Radius of the particles.
    float L_min = 0.5;              //Minimum step corresponding to rand. walk.
    double pi = 3.1415926535897932384626433832; // Declaring pi.
    bool hit = false;
    int i_max = std::numeric_limits<int>::max();
    int counter = 0;
    int counter2 = 0;
    int len = system_length - 1;

    nbr_particles++;
    std::vector<std::vector<complex<double>>> clusters(nbr_particles,
                                  std::vector<complex<double>>(1));
    std::vector<std::complex<double>> particles;
    std::vector<std::vector<int>> num_grid_C(system_length,
                                           std::vector<int>(system_length));
    std::vector<std::vector<int>> num_grid_N(system_length,
                                           std::vector<int>(system_length));
//------------------------------Seeding MT RNG----------------------------------
    //std::mt19937::result_type seed = time(0);
    std::mt19937::result_type seed = 1447081009;
    writeSeed(seed);
    auto rand_dir = std::bind(std::uniform_real_distribution<float>(0,2*pi),
				   std::mt19937(seed));
    auto rand_seed = std::bind(std::uniform_int_distribution<int>(0,i_max),
				   std::mt19937(seed));
//-------------------------------Running prog.----------------------------------
    initGrid(system_length, rand_seed(), nbr_particles, r_p,
             clusters, num_grid_C, num_grid_N);
//    printNumGrid(num_grid_C);
//    printClusters(clusters);
    writeConfigTest(clusters, r_p, system_length, 0);
    writeConfig(clusters, r_p);
    writeConfigColor(clusters, r_p, system_length, 0);
    plotConfig(system_length, std::to_string(0)); 
    std::cout << "clusters.size() = " << clusters.size()-1 << std::endl;
    for (int i = 1; i < 4; ++i){
        std::cout << std::endl;
        std::cout << "Iteration number " << i << std::endl;
        std::cout << "number of clusters = " << clusters.size()-1 << std::endl;
        walkOnGrid(clusters, num_grid_C, num_grid_N, rand_seed(), L_min, pi, r_p,
                   counter);
        writeConfigColor(clusters, r_p, system_length, i);    
        plotConfig(system_length, std::to_string(i));
    }
    std::cout << "clusters.size() = " << clusters.size()-1 << std::endl;
    writeConfig1(clusters, r_p, system_length);
    plotConfigTest(system_length, 0, 1);
//-------testing section--------
    return 0;
}

