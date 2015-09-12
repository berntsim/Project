#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include "header.h"
#include "source.cpp"
#include "parameters.cpp"

using namespace std;

int main(){
//----------------Parameter section. Declare your parameters here:------------
	int arr_size = 51; //must be an odd number
	int nbr_particles = 150;
	
	int perimeter = 4*(arr_size-1)-1;
	vector<vector<int>> grid;
	int x_pos, y_pos;
	bool dummy;
//-------------------------Parameter section ends-----------------------------
//-------------------------Seeding for the MT RNG-----------------------------

	std::mt19937::result_type seed = time(0);
	grid = makeArray(arr_size);
	mt19937 mt_rand(time(0));
	auto rand_step = std::bind(std::uniform_int_distribution<int>(0,3),
				   std::mt19937(seed));
	auto rand_int_perim = std::bind(std::uniform_int_distribution<int>
					(0,perimeter), std::mt19937(seed));
//-------------------------Here the actual program is run---------------------

	int start = rand_int_perim();

	startWalkingParticle(grid, arr_size, perimeter, x_pos, y_pos, start);	
	walk(grid, arr_size, x_pos, y_pos, perimeter, seed, nbr_particles);

	plotGnuplot(arr_size);
	cout << findRadius(grid, arr_size) << endl;
	cout << findFractalDim(grid, arr_size, nbr_particles) << endl;
	return 0;
}
