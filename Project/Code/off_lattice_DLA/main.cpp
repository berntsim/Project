#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <algorithm>

#include "header.h"
#include "source.cpp"
using namespace std;

int main(){
//---------------------------------Parameters-----------------------------------
    int nbr_particles = 1000000;    // Number of particles wanted in cluster.
    float l_0 = 0.00000001;     // Typical length scale of the system.
    float t_0 = 0.001;          // Typical time scale for the system.
    int D_max = 50;              // Parameter used to optimize algorithm.
    int system_length = 13000;     // Size of grids used to describe the syst.
    float r_p = 1.0;            // Radius of particles. Should be alterable.
    float L_min = 1.0;
    double pi = 3.1415926535897932384626433832; // Declaring pi.
    unsigned int counter = 0;
    float r_c = 1.0;


    std::vector<std::vector<float>> ex_pos(nbr_particles,std::vector<float>(2));
    std::vector<float> x_pos;
    std::vector<float> y_pos;
    std::vector<std::vector<int>> on_lattice_cluster(system_length,
                                  std::vector<int>(system_length));
    std::vector<std::vector<int>> on_lattice_distance(system_length,
                                  std::vector<int>(system_length));
    std::vector<std::vector<int>> vicinity((2*D_max+1),
                                            std::vector<int>(2*D_max+1));
    float walker_x_pos = 0;     // Variable used for walking particles.
    float walker_y_pos = 0;     // Used for walking particles.
    float step_L = L_min;           // Step size (recalculated at later stage).
    float step_dir = pi;         // Direction of step.


//------------------------------Seeding MT RNG----------------------------------
    std::mt19937::result_type seed = time(0);
    auto rand_dir = std::bind(std::uniform_real_distribution<float>(0,2*pi),
				   std::mt19937(seed));

//------------------------------Actual program----------------------------------
    initEx_pos(x_pos, y_pos,system_length);
    initOnLatticeCluster(on_lattice_cluster,system_length);
    initOnLatticeDistance(on_lattice_distance, system_length, D_max);
    initVicinity(vicinity,D_max);
    neighbourDist(on_lattice_distance, D_max, system_length/2, system_length/2,
                  system_length);
    runAll(step_L, step_L, seed, walker_x_pos, walker_y_pos, pi, r_p,
           on_lattice_distance, on_lattice_cluster, L_min, D_max, x_pos, y_pos,
           system_length, nbr_particles, counter, r_c);   
    writeGrid(x_pos, y_pos, r_p);
    plotGnuplot(system_length);
    std::cout << "Fractal dimension = " << findFractalDim(x_pos, y_pos, r_p,
                                                          r_c) << std::endl;
//------------------Test Print Section------------------------------------------
    //printDistances(x_pos, y_pos, r_p);
//    for (unsigned int j = 0; j < 1000; ++j){
//        x_pos.push_back(j);
//        y_pos.push_back(j);
//    }
//    for (unsigned int i = 0; i < 10000000; ++i){
//        startWalker(walker_x_pos, walker_y_pos, D_max, x_pos, y_pos, rand_dir(),
//                    pi, r_p, r_c);
//        createStep(step_L, rand_dir(), walker_x_pos, walker_y_pos, pi,
//                             r_p, on_lattice_distance, on_lattice_cluster,
//                             L_min, D_max, x_pos, y_pos, system_length, r_c);
//        killParticle(walker_x_pos, walker_y_pos, x_pos, y_pos, r_p, r_c);
//        counter++;
//    }
//    std::cout << "counter = " << counter << std::endl;
    return 0;
}
