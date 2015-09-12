#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

void randTest(std::mt19937::result_type seed);

vector<vector<int>> makeArray(int arr_len);
//This function creates an array of dimension arrLen x arrLen with zeros,
//except for at the very centre. 

void printGrid(std::vector<std::vector<int>> &lattice, int arr_len);
//This function just prints the array.

void startWalkingParticle(std::vector<std::vector<int>> &lattice, int arr_len,
			  int perim, int &x_pos, int &y_pos,
			  int st);
// This function will place out a particle on the perimeter of our grid, which
// will perform a random walk.

bool validStep(int &arr_len, int x_pos, int y_pos, int step);
// This checks if the step taken by the particle will bring it off grid.

bool nextToCluster(std::vector<std::vector<int>> &lattice,int x_pos, int y_pos,
		   int arr_len);
// This function checks if the walking particle has arrived at cluster.

bool takeStep(std::vector<std::vector<int>> &lattice, int arr_len, int &x_pos,
		int &y_pos, int step, int perim, int st);
// This function performs the step. If the particle walks off the grid, this
// function returns false, otherwise true.

void performWalk(std::vector<std::vector<int>> &lattice, int arr_len,
		 int x_pos, int y_pos, int nbr_particles, int perim);
// This function puts all the other functions together, and performs the walk.

void walk(std::vector<std::vector<int>> &lattice, int arr_len,
		int &x_pos,int &y_pos,int perim,std::mt19937::result_type seed,
		int nbr_particles);
void writeGrid(std::vector<std::vector<int>> &lattice, int arr_len);
// This function writes the grid to a file, so that it's easy to plot in gnuplot.

void testWriteGrid(std::vector<std::vector<int>> &lattice, int arr_len);

void plotGnuplot(int &arr_size);

float findRadius(std::vector<std::vector<int>> &lattice, int arr_len);

float findFractalDim(std::vector<std::vector<int>> &lattice, int arr_len,
		     int nbr_part);
void plotJonas(int &arr_len, int part_cluster);

void writeTmp(std::vector<std::vector<int>> &lattice, int arr_len);
