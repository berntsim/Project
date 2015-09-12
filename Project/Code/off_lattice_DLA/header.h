#include <iostream>
#include <vector>
#include <fstream>

void initEx_pos(std::vector<float> &x_pos, std::vector<float> &y_pos,
                int sys_len);

void initOnLatticeCluster(std::vector<std::vector<int>> &on_lattice_cluster, 
                                int system_length);

void initOnLatticeDistance(std::vector<std::vector<int>> &on_lattice_distance,
                           int sys_len, int D_max);

void neighbourDist(std::vector<std::vector<int>> &on_lattice_distance,
                   int D_max, int row, int col, int sys_len);

void testDist(std::vector<std::vector<int>> &on_lattice_distance,
                   int D_max, int row, int col, int sys_len);







void printGrid(std::vector<std::vector<int>> &grid, int sys_len);



