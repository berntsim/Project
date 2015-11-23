#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <algorithm>
#include <complex>
#include <limits>

void initGrid(int system_length, std::mt19937::result_type seed,
              int nbr_particles,
              std::vector<std::vector<std::complex<double>>> &clusters,
              std::vector<std::vector<int>> &num_grid);


void walkOnGrid(std::vector<std::vector<std::complex<double>>> &clusters,
                std::vector<std::vector<int>> &num_grid,
                std::mt19937::result_type seed, float L_step, double pi,
                float r_p, int &counter);

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

void updateNum_gridTest(std::vector<std::vector<int>> &num_grid, int X,
                        int Y, double nextX, double nextY);

void updateNum_gridWalk(std::vector<std::vector<int>> &num_grid,
                        std::complex<double> f,
                        std::complex<double> t,
                        std::vector<std::complex<double>> mt);


int checkOrientation(double X, double Y, int len);
    //Checks where on the grid the particle is. This takes into account that it
    //may be on the perimeter, either on the sides or in the corners.



bool isFree(std::vector<std::vector<int>> num_grid, double nextX, double nextY,
            double X, double Y);
    //This function checks if the cell X,Y is empty or not.



void makeStepOrg(std::vector<std::vector<int>> num_grid, double X,
              double Y, int orientation, int step_dir, int len, float L_step,
              double &nextX, double &nextY);
    //This function changes the values of nextX&Y to be that of the walker. I.e.
    //it will only change if a valid step is possible.



int checkHit(std::vector<std::vector<int>> num_grid, double X, double Y,
              int orientation, int len);
    //Check if a point has any neighbours. If the cell has a non-zero neighbour,
    //the function returns the cluster ID of the neighbour, so that one may join
    //them together. Currently only checks nearest neighbours.


void makeStep(double X, double Y, int orientation, int step_dir, int len,
              float L_step, double &nextX, double &nextY);


void joinClusters(std::vector<std::vector<std::complex<double>>> &clusters,
                  std::vector<std::vector<int>> &num_grid, int A, int B);
    //Joins the elements of clusters[A] and clusters[B], and deletes the highest
    //indexed afterwards. The reduced parameter is there to keep track of how
    //many elements have been erased. 



void joinInit(std::vector<std::vector<std::complex<double>>> &clusters,
              std::vector<std::vector<int>> &num_grid);
    //To be used after initiating the arrays to check if some particles start
    //out as clusters, and join these together.


void updateNum_grid(std::vector<std::vector<int>> &num_grid, int start,
                    std::vector<std::vector<std::complex<double>>> clusters);

double clusterSize(std::vector<std::vector<std::complex<double>>> clusters,
                 int s, int system_length);

void clustersTime(std::vector<std::vector<std::complex<double>>> clusters,
                  int system_length);







void printClusters(std::vector<std::vector<std::complex<double>>> &clusters);

void printGrid(std::vector<std::vector<std::complex<double>>> grid);

void writeConfig(std::vector<std::vector<std::complex<double>>> clusters,
                 float r_p);

void writeConfig1(std::vector<std::vector<std::complex<double>>> clusters,
                 float r_p, int system_length);

void writeConfigColor(std::vector<std::vector<std::complex<double>>> clusters,
                 float r_p, int system_length, double iteration);

void plotConfig(int iteration, int system_length);

void plot(int system_length);

void printNumGrid(std::vector<std::vector<int>> grid);
