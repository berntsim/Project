#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <algorithm>
#include <complex>
#include <limits>

void initGrid(int system_length, std::mt19937::result_type seed,
              int nbr_particles, float r_p,
              std::vector<std::vector<std::complex<double>>> &clusters,
              std::vector<std::vector<int>> &num_grid_C,
              std::vector<std::vector<int>> &num_grid_N);
//initalizes the domain of interest. This function makes sure that nbr_particles
//are placed on a system_length by system_length domain. The particles have
//radius r_p, and they may not overlap or tang

void walkOnGrid(std::vector<std::vector<std::complex<double>>> &clusters,
                std::vector<std::vector<int>> &num_grid,
                std::mt19937::result_type seed, float L_step, double pi,
                float r_p, int &counter);

void findShortestDist(int len, double &x_c, double &y_c, double &X, double &Y,
                      float r_p);

void checkDestination(std::vector<std::vector<std::complex<double>>> clusters,
                      std::vector<std::vector<int>> num_grid_C, int len,
                      std::vector<std::vector<int>> num_grid_N, int row,
                      int col, std::vector<std::complex<double>> &to_check,
                      float r_p);

void checkMovedTo(std::vector<std::vector<int>> num_grid, 
                  std::vector<double> &moved_to_X,
                  std::vector<double> &moved_to_Y,
                  double X, double Y, double nextX,
                  double nextY);

void takeStep(std::vector<std::vector<std::complex<double>>> clusters,
              std::vector<std::vector<int>> num_grid, double X, double Y,
              int step_dir, int len, float L_step);

void updateClusters(std::vector<std::vector<std::complex<double>>> &clusters,
                    int i, int j, double nextX, double nextY);

void updateNum_grid_C(std::vector<std::vector<int>> &num_grid_C,
                      std::vector<std::complex<int>> mt,
                      std::complex<double> from, std::complex<double> to);

void updateNum_grid_N(std::vector<std::vector<int>> num_grid_N,
                      std::vector<std::vector<int>> &ngn_tmp,
                      int rt, int ct, int rf, int cf);

void updateNum_gridWalk(std::vector<std::vector<int>> &num_grid,
                        std::complex<double> f,
                        std::complex<double> t,
                        std::vector<std::complex<double>> mt);

void makeStep(double X, double Y, double step_dir, int len, float L_step,
              double &nextX, double &nextY);

float LHit(float step_L, double step_dir, double x_c, double y_c, double x_p,
           double y_p, float r_p);

void joinClusters(std::vector<std::vector<std::complex<double>>> &clusters,
                  std::vector<std::vector<int>> &num_grid, int A, int B,
                  std::vector<std::vector<int>> &num_grid_N);
    //Joins the elements of clusters[A] and clusters[B], and deletes the highest
    //indexed afterwards. The reduced parameter is there to keep track of how
    //many elements have been erased. 

void updateNum_grid(std::vector<std::vector<int>> &num_grid, int start,
                    std::vector<std::vector<std::complex<double>>> clusters);








void printClusters(std::vector<std::vector<std::complex<double>>> &clusters);

void writeConfig(std::vector<std::vector<std::complex<double>>> clusters,
                 float r_p);

void writeConfig1(std::vector<std::vector<std::complex<double>>> clusters,
                 float r_p, int system_length);

void writeConfigColor(std::vector<std::vector<std::complex<double>>> clusters,
                 float r_p, int system_length, int iteration);

void plotConfig(int iteration, int system_length);

void plot(int system_length);

void printNumGrid(std::vector<std::vector<int>> grid);

void writeSeed(std::mt19937::result_type seed);

void printSingleDist(double X, double Y, double x_c, double y_c,
                     double nextX, double nextY,        
                     std::vector<std::vector<int>> num_grid_C,
                     std::vector<std::vector<int>> num_grid_N,
                     std::vector<std::vector<std::complex<double>>> clusters);

void printColInfo(double X, double Y, double x_c, double y_c,
                  std::vector<std::vector<std::complex<double>>> clusters);
