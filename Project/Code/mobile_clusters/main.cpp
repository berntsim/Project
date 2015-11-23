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
    int nbr_particles = 25000;         //Number of particles we want in the syst.
    int system_length = 400;         //Dimension of the array used for the syst.
    float r_p = 0.5;                //Radius of the particles.
    float L_min = 1.0;              //Minimum step corresponding to rand. walk.
    double pi = 3.1415926535897932384626433832; // Declaring pi.
    float d_f;                      //Fractal dimension (to be calculated).
    bool hit = false;
    int i_max = std::numeric_limits<int>::max();
    int counter = 0;
    int counter2 = 0;
    int len = system_length - 1;
    int s = 10;
    double dt = 0;
    double t = 0;

    nbr_particles++;
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
    initGrid(system_length, rand_seed(), nbr_particles,
             clusters, num_grid);
    joinInit(clusters, num_grid);
//    writeConfigTest(clusters, r_p, system_length, 0);
//    plot(system_length, std::to_string(0));
    std::cout << "starting the walk" << std::endl;
    //int iteration_len = 100;
    std::ofstream out_stream;
	out_stream.open("cluster_size_distribution.txt");
    //for (int i = 0; i < iteration_len; ++i){
    while ( t < 0.1){    
        std::cout << "t = " << t << std::endl;
        out_stream << t << " " << clusterSize(clusters, s, system_length)
                   << std::endl;
        t += double(s)/double(nbr_particles);
        std::cout << double(s)/double(nbr_particles) << std::endl;
        writeConfigColor(clusters, r_p, system_length, t);    
        plotConfig(system_length, std::to_string(t));
        walkOnGrid(clusters, num_grid, rand_seed(), L_min, pi, r_p, counter);
    }
    clustersTime(clusters, system_length);
    std::cout << "counter = " << counter << std::endl;
    std::cout << "#clusters = " << clusters.size()-1 << std::endl;
//    walkOnGrid(clusters, num_grid, rand_seed(), L_min, pi, r_p, counter);
//    writeConfigColor(clusters, r_p, system_length, iteration_len+1);
//    plotConfig(system_length, std::to_string(iteration_len+1));
//    plot(system_length, std::to_string(iteration_len+1));

    
//-------testing section--------
//    double nextX, nextY;
//    std::vector<double> moved_to_X, moved_to_Y, deleted_X, deleted_Y;
//    for (int i = 0; i < 10000; ++i){
//        int orientation = checkOrientation(1, 1, len);
        //makeStepOrg(num_grid, 1, 1, orientation, 1, len, L_min, nextX, nextY);
        //makeStep(1, 1, orientation, 1, len, L_min, nextX, nextY);
        //checkMovedTo(num_grid, deleted_X, deleted_Y, 1, 1, nextX, nextY);
        //updateNum_grid(num_grid, 1, clusters);
        //updateNum_gridTest(num_grid, 1, 1, 2, 1, 1);
//    }

    return 0;
}

