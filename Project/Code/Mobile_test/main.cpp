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
    int nbr_particles = 500;         //Number of particles we want in the syst.
    int system_length = 200;         //Dimension of the array used for the syst.
    float r_p = 0.5;                //Radius of the particles.
    float L_min = 1.0;              //Minimum step corresponding to rand. walk.
    double pi = 3.1415926535897932384626433832; // Declaring pi.
    float d_f;                      //Fractal dimension (to be calculated).
    bool hit = false;
    int i_max = std::numeric_limits<int>::max();
    int counter1 = 0;
    int counter2 = 0;

    nbr_particles++;
    std::vector<std::vector<complex<double>>> on_lattice_cluster(
            system_length, std::vector<complex<double>>(system_length));
    std::vector<std::vector<complex<double>>> clusters(nbr_particles,
                                  std::vector<complex<double>>(1));
    std::vector<std::complex<double>> particles;
    std::vector<std::vector<int>> num_grid(system_length,
                                           std::vector<int>(system_length));

//------------------------------Seeding MT RNG----------------------------------
    std::mt19937::result_type seed = time(0);
    auto rand_dir = std::bind(std::uniform_real_distribution<float>(0,2*pi),
				   std::mt19937(seed));
    auto rand_seed = std::bind(std::uniform_int_distribution<int>(0,i_max),
				   std::mt19937(seed));

//-------------------------------Running prog.----------------------------------
    initGrid(on_lattice_cluster, system_length, rand_seed(), nbr_particles,
             clusters, num_grid);
    joinInit(clusters, num_grid);
    writeConfig(clusters, r_p);
    std::cout << "starting the walk" << std::endl;
    for (int i = 0; i < 100; ++i){
        std::cout << i << std::endl;
        walkOnGrid(on_lattice_cluster, clusters, num_grid, rand_seed(), L_min, pi, r_p);
    }
    std::cout << "Walk is performed" << std::endl;
    writeConfig1(clusters, r_p, system_length);
    plotConfig(1,system_length);
    return 0;
}

