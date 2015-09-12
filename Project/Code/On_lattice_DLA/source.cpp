#include "header.h"
#include <vector>
#include <random>
#include <fstream>

using namespace std;

std::vector<std::vector<int>> makeArray(int arr_len){
	std::vector<std::vector<int>> lattice(arr_len,
	std::vector<int>(arr_len));
	for (unsigned int i = 0; i < arr_len; ++i){
		for (unsigned int j = 0; j < arr_len; ++j){
			if (( i == arr_len/2) && (j == arr_len/2)){
				lattice[i][j] = 2;
			}
			else {
				lattice[i][j] = 0;
			}
		}
	}
	return lattice;
}

void printGrid(std::vector<std::vector<int>> &lattice, int arr_len){
	for (unsigned int i = 0; i < arr_len; ++i){
		for (unsigned int j = 0; j < arr_len; ++j){
			cout << lattice[i][j];
		}
		cout << endl;
	}
	cout << endl;
}

void startWalkingParticle(std::vector<std::vector<int>> &lattice, int arr_len,
			 int perim, int &x_pos, int &y_pos,
			 int st){
	for (unsigned int i = 0; i < arr_len; ++i){
		for (unsigned int j = 0; j < arr_len; ++j){
			if ( lattice[i][j] == 2 ){
				lattice[i][j] = 2;
			}
			else {
				lattice[i][j] = 0;
			}
		}
	}

	int start_pos = st;
	if (start_pos < arr_len){
		lattice[0][start_pos] = 1;
		x_pos = start_pos;
		y_pos = 0;
	}
	else if ((start_pos >= arr_len) && (start_pos < 2*arr_len-1)){ 
		lattice[start_pos-(arr_len-1)][arr_len-1] = 1;
		x_pos = arr_len-1;
		y_pos = start_pos-(arr_len-1);
	}
	else if ((start_pos >= 2*arr_len-1) && (start_pos < 3*arr_len-2)){
		lattice[arr_len-1][3*arr_len-3-start_pos] = 1;
		x_pos = 3*arr_len-3-start_pos;
		y_pos = arr_len-1;
	}
	else if ((start_pos >= 3*arr_len-2) && (start_pos <= perim)){
		lattice[4*arr_len-4-start_pos][0] = 1;
		x_pos = 0;
		y_pos = 4*arr_len-4-start_pos;
	}
	else {
		cout << "Error in startWalkingParticle!" << endl;
	}

}

bool validStep(int &arr_len, int x_pos, int y_pos, int step){
	if ((y_pos == arr_len-1) && (step == 0)){
		return false;
	}
	else if ((y_pos == 0) && (step == 2)){
		return false;
	}
	else if ((x_pos == 0) && (step == 3)){
		return false;
	}
	else if ((x_pos == arr_len-1) && (step == 1)){
		return false;
	}
	else {
		return true;
	}
}

bool nextToCluster(std::vector<std::vector<int>> &lattice, int x_pos,int y_pos,
		   int arr_len){
	if ((x_pos != arr_len-1) && (y_pos != arr_len-1) && (x_pos != 0) && 
	    (y_pos != 0)){
		if ((lattice[y_pos][x_pos+1]==2) || (lattice[y_pos][x_pos-1]==2)
	   	|| (lattice[y_pos-1][x_pos]==2) || (lattice[y_pos+1][x_pos]==2)){
			return true;
		}
		else {
			return false;
		}
	}
	else {
		return false;
	}
}

bool takeStep(std::vector<std::vector<int>> &lattice, int arr_len, int &x_pos,
		int &y_pos, int step, int perim, int st){
	lattice[y_pos][x_pos]=0;
	if ((validStep(arr_len, x_pos, y_pos, step)) && (step == 0)){
		y_pos++;
		lattice[y_pos][x_pos]=1;
		return true;
	}
	else if ((validStep(arr_len, x_pos, y_pos, step)) && (step == 1)){
		x_pos++;
		lattice[y_pos][x_pos]=1;
		return true;
	}
	else if ((validStep(arr_len, x_pos, y_pos, step)) && (step == 2)){
		y_pos--;
		lattice[y_pos][x_pos]=1;
		return true;
	}
	else if ((validStep(arr_len, x_pos, y_pos, step)) && (step == 3)){
		x_pos--;
		lattice[y_pos][x_pos]=1;
		return true;
	}
	else {
		lattice[y_pos][x_pos]=0;
		startWalkingParticle(lattice, arr_len, perim, x_pos, y_pos,st);
		return false;
	}
}

void walk(std::vector<std::vector<int>> &lattice, int arr_len,
		 int &x_pos, int &y_pos, int perim,
		 std::mt19937::result_type seed, int nbr_particles){
	auto rand_step = std::bind(std::uniform_int_distribution<int>(0,3),
				   std::mt19937(seed));		
	auto rand_int_perim = std::bind(std::uniform_int_distribution<int>
					(0,perim), std::mt19937(seed));	
	int start = rand_int_perim();
	int step = rand_step();
	bool on_grid = true;
	bool attached = false;
	int particles_in_cluster = 1;

	while (particles_in_cluster < nbr_particles){
		startWalkingParticle(lattice, arr_len, perim, x_pos, y_pos, 
				     rand_int_perim());
		attached = false;
		while (attached == false){
			
			step = rand_step();
			on_grid = takeStep(lattice, arr_len, x_pos, y_pos, step, 
				  	   perim, rand_int_perim());
			if (nextToCluster(lattice, x_pos, y_pos, arr_len)){
				particles_in_cluster++;
				lattice[y_pos][x_pos] = 2;
				attached = true;
			}
		}
	writeTmp(lattice, arr_len);
	plotJonas(arr_len, particles_in_cluster);
	cout << particles_in_cluster << endl;
	}
	testWriteGrid(lattice, arr_len);
}

void writeGrid(std::vector<std::vector<int>> &lattice, int arr_len){
	ofstream out_stream;
	out_stream.open("data/test.txt");
	for (unsigned int i = 0; i < arr_len; ++i){
		for (unsigned int j = 0; j < arr_len; ++j){
			out_stream << i << ","  << j<< ","  << lattice[i][j] << endl;
			if (j == arr_len-1){
				out_stream << endl;
			}
		}
	}
	out_stream.close( );
}

void writeTmp(std::vector<std::vector<int>> &lattice, int arr_len){
	ofstream out_stream;
	out_stream.open("tmp.txt");
	for (unsigned int i = 0; i < arr_len; ++i){
		for (unsigned int j = 0; j < arr_len; ++j){
			if (lattice[i][j] == 2){
				out_stream << i << "	" << j << endl;
			}
		}
	}
	out_stream.close( );
}


void testWriteGrid(std::vector<std::vector<int>> &lattice, int arr_len){
	ofstream out_stream;
	out_stream.open("data/test.txt");
	for (unsigned int i = 0; i < arr_len; ++i){
		for (unsigned int j = 0; j < arr_len; ++j){
			if (lattice[i][j] == 2){
				out_stream << i << "	" << j << endl;
			}
		}
	}
	out_stream.close( );
}

void plotGnuplot(int &arr_len){
	ofstream out_stream;
	out_stream.open("gnuplotter.gnu");
	out_stream << "set terminal png size 600,500 enhanced font \"Helvetica,12\" "
		   << endl;
	out_stream << "set output \"figures/test.png\" " << endl;
	out_stream << "set xlabel \"x-dim\" " << endl;
	out_stream << "set ylabel \"y-dim\"" << endl;
	out_stream << "set xrange [0:" << arr_len << "]" << endl;
	out_stream << "set yrange [0:" << arr_len << "]" << endl;
	out_stream << "plot \"data/test.txt\" with points pointtype 16" << endl;
	system("gnuplot gnuplotter.gnu");
}

float findRadius(std::vector<std::vector<int>> &lattice, int arr_len){
	float rad= 0;
	int seed_x = arr_len/2;
	int seed_y = seed_x;
	for (int i = 0; i < arr_len; ++i){
		for (int j = 0; j < arr_len; ++j){
			if ((sqrt(pow(i-seed_x,2)+pow(j-seed_y,2)) > rad) &&
			    (lattice[i][j] == 2)){
				rad = sqrt(pow(i-seed_x,2)+pow(j-seed_y,2));
			}
		}
	}
	return rad;
}

float findFractalDim(std::vector<std::vector<int>> &lattice, int arr_len, int nbr_part){
	float r = findRadius(lattice, arr_len);
	int M = nbr_part;
	return log(M)/log(r);
}

void plotJonas(int &arr_len, int part_cluster){
	ofstream out_stream;
	out_stream.open("gnuplotter.gnu");
	out_stream << "set terminal png size 1200,1000 enhanced font \"Helvetica,12\" "
		   << endl;
	out_stream << "set output \"animations/test" << part_cluster << ".png\" " 
		   << endl;
	out_stream << "set xlabel \"x-dim\" " << endl;
	out_stream << "set ylabel \"y-dim\"" << endl;
	out_stream << "set xrange [0:" << arr_len << "]" << endl;
	out_stream << "set yrange [0:" << arr_len << "]" << endl;
	out_stream << "plot \"tmp.txt\" with points pointtype 16" << endl;
	system("gnuplot gnuplotter.gnu");
}
