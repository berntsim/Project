#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <algorithm>
#include <complex>
#include <limits>

#include "header.h"
#include "off_grid.cpp"
#include "on_grid.cpp"
#include "print.cpp"
using namespace std;
int main(){
//---------------------------------Parameters-----------------------------------
    int nbr_particles = 5;         //Number of particles we want in the syst.
    int system_length = 5;         //Dimension of the array used for the syst.
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

//------------------------------Program not using Grid--------------------------
//    initClusters(clusters, seed, nbr_particles, system_length);
//    writeConfig(clusters,r_p);

//    for (unsigned int i = 0; i < 10; ++i){
//        walk(clusters, rand_seed(), L_min, pi, system_length, counter1,
//                   counter2);
//        std::cout << i << std::endl;
//    }
//    printClusters(clusters);
//    writeConfig1(clusters, r_p, system_length);
//    plotConfig(2, system_length);
//    std::cout << "counter 1 = " << counter1 << std::endl;
//    std::cout << "counter 2 = " << counter2 << std::endl;
//------------------------------Program using Grid------------------------------
    initGrid(on_lattice_cluster, system_length, seed, nbr_particles, clusters,
             num_grid);         
    printGrid(on_lattice_cluster);
    printTest(num_grid);
    std::cout << std::endl;

    printClusters(clusters);
    //printPartPos(on_lattice_cluster);
//    printParticles(particles);
    writeConfigGrid(on_lattice_cluster, r_p);
    walkOnGrid(on_lattice_cluster, clusters, num_grid, rand_seed(), L_min, pi,
               r_p);
    printGrid(on_lattice_cluster);
    printTest(num_grid);
    printClusters(clusters);
//    printGrid(on_lattice_cluster);
//    printPartPos(on_lattice_cluster);
    writeConfigGrid1(on_lattice_cluster, r_p);
    plotConfig(1, system_length);
    return 0;
}
