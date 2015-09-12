#include <iostream>
#include <vector>
#include <fstream>

#include "header.h"
#include "source.cpp"
using namespace std;

int main(){
//---------------------------------Parameters-----------------------------------
    int nbr_particles = 100;
    float l_0 = 0.00000001; 
    float t_0 = 0.001;
    int D_max = 10;
    float system_length = 100; 

    std::vector<std::vector<float>> ex_pos(nbr_particles,std::vector<float>(2));
    std::vector<float> x_pos;
    std::vector<float> y_pos;
    std::vector<std::vector<int>> on_lattice_cluster(system_length,
                                  std::vector<int>(system_length));
    std::vector<std::vector<int>> on_lattice_distance(system_length,
                                  std::vector<int>(system_length));
    std::vector<std::vector<int>> vicinity((2*D_max+1),
                                            std::vector<int>(2*D_max+1));
//------------------------------Actual program----------------------------------
    initEx_pos(x_pos, y_pos,system_length);
    initOnLatticeCluster(on_lattice_cluster,system_length);
    
//------------------Test Print Section------------------------------------------
    std::cout << x_pos[0] << std::endl;
    std::cout << y_pos[0] << std::endl;

    return 0;
}
