#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <algorithm>


//------------------------------Using Grid------------------------------
void initGrid(std::vector<std::vector<std::complex<double>>> &grid,
              int system_length, std::mt19937::result_type seed, 
              int nbr_particles,
              std::vector<std::vector<std::complex<double>>> &clusters,
              std::vector<std::vector<int>> &num_grid);

void updateGrid(std::vector<std::vector<std::complex<double>>> clusters,
                std::vector<std::vector<std::complex<double>>> &grid);

void walkOnGrid(std::vector<std::vector<std::complex<double>>> &grid,
                std::vector<std::vector<std::complex<double>>> &clusters,
                std::vector<std::vector<int>> &num_grid,
                std::mt19937::result_type seed, float L_step, double pi,
                float r_p);

int joinClustersGrid(std::vector<std::vector<std::complex<double>>> &clusters,
                     std::vector<std::vector<int>> &num_grid, int A, int B,
                     int &reduced);

int checkHit(std::vector<std::vector<std::complex<double>>> &clusters,
              std::vector<std::vector<std::complex<double>>> grid,
              std::vector<std::vector<int>> &num_grid,
              int row, int col, float r_p, int &reduced);

//------------------------------Not using Grid----------------------------------

void initClusters(std::vector<std::vector<std::complex<double>>> &clusters,
                  std::mt19937::result_type seed, int nbr_particles,
                  int system_length);

void findCM(std::vector<std::complex<double>> pos,
            std::complex<double> &ret_var);

void joinClusters(std::vector<std::vector<std::complex<double>>> &clusters,
                  int A, int B);

bool checkHit(std::vector<std::vector<std::complex<double>>> &clusters,
              int cluster_number, int particle_number, int &A, int &B,
              int &counter2);

void walk(std::vector<std::vector<std::complex<double>>> &clusters,
                std::mt19937::result_type seed, float L_step, double pi,
                int system_length, int &counter1, int &counter2);

//---------------------------------Print file-----------------------------------
void printClusters(std::vector<std::vector<std::complex<double>>> &clusters);

void writeConfig(std::vector<std::vector<std::complex<double>>> grid,
                 float r_p);

void plotConfig(int iteration);

void writeConfig1(std::vector<std::vector<std::complex<double>>> clusters,
                 float r_p);
void printGrid(std::vector<std::vector<std::complex<double>>> const grid);
