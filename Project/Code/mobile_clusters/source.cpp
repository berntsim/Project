#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <algorithm>
#include <math.h>
#include <string>

#include "header.h"

void initGrid(int system_length, std::mt19937::result_type seed,
              int nbr_particles,
              std::vector<std::vector<std::complex<double>>> &clusters,
              std::vector<std::vector<int>> &num_grid){
    int counter = 1;
    double tmp1;
    double tmp2;
    for (unsigned int i = 0; i < num_grid.size(); ++i){
        for (unsigned int j = 0; j < num_grid.size(); ++j){
            num_grid[i][j]=0;
        }
    }
    //auto coord = std::bind(std::uniform_real_distribution<float>
    //                       (0, double(grid.size()-1)), std::mt19937(seed));
    auto coord = std::bind(std::uniform_int_distribution<int>
                           (1, num_grid.size()-1), std::mt19937(seed));
    while (counter < nbr_particles){
        tmp1 = coord();
        tmp2 = coord();
        if (num_grid[floor(tmp1)][floor(tmp2)] == 0){
            clusters[counter][0] = {tmp1, tmp2};
            num_grid[floor(tmp1)][floor(tmp2)] = counter;
            ++counter;
        }
    }
}

void walkOnGrid(std::vector<std::vector<std::complex<double>>> &clusters,
                std::vector<std::vector<int>> &num_grid,
                std::mt19937::result_type seed, float L_step, double pi,
                float r_p, int &counter){
    int step_dir, orientation;
    int len = num_grid.size()-1;
    double X, Y, nextX, nextY;
    int reduced = 0;
    int hit = 0;
    int has_hit_number = 0;
    bool has_hit = false;
    bool can_move = false;
    auto rand_dir_int = std::bind(std::uniform_real_distribution<float>(0,4),
	                std::mt19937(seed));
    for (int i = 1; i < clusters.size(); ++i){
        can_move = true;
        has_hit = false;
        step_dir = rand_dir_int();
        std::vector<std::complex<double>> moved_to, moved_from;
        for (int j = 0; j < clusters[i].size(); ++j){
            X = std::real(clusters[i][j]);
            Y = std::imag(clusters[i][j]);
            orientation = checkOrientation(X, Y, len);
            makeStep(X, Y, orientation, step_dir, len, L_step,
                      nextX, nextY);
            if ((num_grid[nextX][nextY] != 0) && 
                (num_grid[X][Y] != num_grid[nextX][nextY])){
                //if the particle cannot go to the new position, it remains at
                //at rest. This relies on propperly updating num_grid.
                nextX = X;
                nextY = Y;
            }
            if ((X != nextX) || (Y != nextY)){
                //if the particle can move, we also have to make sure that the
                //entire cluster can move. Thus, we have to store where it was
                //and where it was going, and then once we know that the move
                //is possible perform it.
                moved_to.push_back({nextX, nextY});
                moved_from.push_back({X,Y});
                }
            else {
                //There is a particle which cannot move -> We have a hit,
                //and the cluster cannot move.
                can_move = false;
            }
        }
        if (can_move){
            for (int j = 0; j < moved_from.size(); ++j){
                updateClusters(clusters, i, j, real(moved_to[j]),
                               imag(moved_to[j]));
                //updateNum_gridTest(num_grid, moved_from_X[j], moved_from_Y[j],
                //                   moved_to_X[j], moved_to_Y[j]);
                updateNum_gridWalk(num_grid, moved_from[j], moved_to[j],
                                   moved_to);
                orientation = checkOrientation(real(moved_to[j]),
                                               imag(moved_to[j]), len);
                hit = checkHit(num_grid, real(moved_to[j]), imag(moved_to[j]),
                               orientation, len);
                if (hit != 0){
                    has_hit = true;
                    has_hit_number = hit;
                }
            }
            if (has_hit != 0){
                joinClusters(clusters, num_grid, i, has_hit_number);
                updateNum_grid(num_grid, std::max(i, has_hit_number), clusters);
            }
            //for (int k = 0; k < deleted_X.size(); ++k){
            //    num_grid[deleted_X[k]][deleted_Y[k]] = i;
            //}
            
        }
        else if (hit != 0){
            //std::cout << "Hit registered" << std::endl;
            counter++;
        } 
    }
}


void checkMovedTo(std::vector<std::vector<int>> num_grid, 
                  std::vector<double> &moved_to_X,
                  std::vector<double> &moved_to_Y,
                  double X, double Y, double nextX,
                  double nextY){
    if (num_grid[floor(X)][floor(Y)] == num_grid[floor(nextX)][floor(nextY)]){
        moved_to_X.push_back(nextX);
        moved_to_Y.push_back(nextY);
    } 
}

void updateClusters(std::vector<std::vector<std::complex<double>>> &clusters,
                    int i, int j, double nextX, double nextY){
    clusters[i][j] = {nextX, nextY};
}

void updateNum_gridTest(std::vector<std::vector<int>> &num_grid, int X,
                        int Y, double nextX, double nextY){
    num_grid[floor(nextX)][floor(nextY)] = num_grid[X][Y];
    num_grid[X][Y] = 0;
}

void updateNum_gridWalk(std::vector<std::vector<int>> &num_grid,
                        std::complex<double> f,
                        std::complex<double> t,
                        std::vector<std::complex<double>> mt){
    if (std::find(mt.begin(), mt.end(), f) != mt.end()){
        num_grid[floor(real(t))][floor(imag(t))] = 
                num_grid[floor(real(f))][floor(imag(f))];    
    }
    else {
        num_grid[floor(real(t))][floor(imag(t))] = 
                num_grid[floor(real(f))][floor(imag(f))];
        num_grid[floor(real(f))][floor(imag(f))] = 0;
    }
}


int checkOrientation(double X, double Y, int len){
    if ((Y == len) && ((X != 0) && (X != len))){
        return 0;
    }
    else if ((X == len) && ((Y != 0) && (Y != len))){
        return 1;
    }
    else if ((Y == 0) && ((X != 0) && (X != len))){
        return 2;
    }
    else if ((X == 0) && ((Y != 0) && (Y != len))){
        return 3;
    }
    else if ((X == len) && (Y == len)){
        return 5;
    }
    else if ((X == len) && (Y == 0)){
        return 6;
    }
    else if ((X == 0) && (Y == 0)){
        return 7;
    }
    else if ((X == 0) && (Y == len)){
        return 8;
    }
    else {
        return 4;
    }
}
//onPerim returns the orientation of the walkers position. See notes in textbook

bool isFree(std::vector<std::vector<int>> num_grid, double nextX, double nextY,
            double X, double Y){
    int A = num_grid[X][Y];
    if ((num_grid[nextX][nextY] == 0) || (num_grid[nextX][nextY] == A)){
        //We return free if it is in the same cluster, since the particle there
        //will move silimarly.
        return true;
    }
    else {
        return false;
    }

}
//Checks if the cell at X,Y is occupied or not.

void makeStep(double X, double Y, int orientation, int step_dir, int len,
              float L_step, double &nextX, double &nextY){ 
    if (orientation == 0){
        // we must now chech if we can move to the cell in Step_dir direction.
        if (step_dir == 0){
            nextX = X;
            nextY = Y - len;
        }
        else if (step_dir == 1){
            nextX = X + L_step;
            nextY = Y;
        }
        else if (step_dir == 2){
            nextX = X;
            nextY = Y - L_step;
        }
        else if (step_dir == 3){
            nextX = X - L_step;
            nextY = Y;
        }
    }
    else if (orientation == 1){
        if (step_dir == 0){
            nextX = X;
            nextY = Y + L_step;
        }
        else if (step_dir == 1){
            nextX = X - len;
            nextY = Y;
        }
        else if (step_dir == 2){
            nextX = X;
            nextY = Y - L_step;
        }
        else if (step_dir == 3){
            nextX = X - L_step;
            nextY = Y;
        }
    }
    else if (orientation == 2){
        if (step_dir == 0){
            nextX = X;
            nextY = Y + L_step;
        }
        else if (step_dir == 1){
            nextX = X + L_step;
            nextY = Y;
        }
        else if (step_dir == 2){
            nextX = X;
            nextY = Y + len;
        }
        else if (step_dir == 3){
            nextX = X - L_step;
            nextY = Y;
        }
    }
    else if (orientation == 3){
        if (step_dir == 0){
            nextX = X;
            nextY = Y + L_step;
        }
        else if (step_dir == 1){
            nextX = X + L_step;
            nextY = Y;
        }
        else if (step_dir == 2){
            nextX = X;
            nextY = Y - L_step;
        }
        else if (step_dir == 3){
            nextX = X + len;
            nextY = Y;
        }
    }
    else if (orientation == 4){
        if (step_dir == 0){
            nextX = X;
            nextY = Y + L_step;
        }
        else if (step_dir == 1){
            nextX = X + L_step;
            nextY = Y;
        }
        else if (step_dir == 2){
            nextX = X;
            nextY = Y - L_step;
        }
        else if (step_dir == 3){
            nextX = X - L_step;
            nextY = Y;
        }
    }
    else if (orientation == 5){
        if (step_dir == 0){
            nextX = X;
            nextY = Y - len;
        }
        else if (step_dir == 1){
            nextX = X - len;
            nextY = Y;
        }
        else if (step_dir == 2){
            nextX = X;
            nextY = Y - L_step;
        }
        else if (step_dir == 3){
            nextX = X - L_step;
            nextY = Y;
        }
    }
    else if (orientation == 6){
        if (step_dir == 0){
            nextX = X;
            nextY = Y + L_step;
        }
        else if (step_dir == 1){
            nextX = X - len;
            nextY = Y;
        }
        else if (step_dir == 2){
            nextX = X;
            nextY = Y + len;
        }
        else if (step_dir == 3){
            nextX = X - L_step;
            nextY = Y;
        }
    }
    else if (orientation == 7){
        if (step_dir == 0){
            nextX = X;
            nextY = Y + L_step;
        }
        else if (step_dir == 1){
            nextX = X + L_step;
            nextY = Y;
        }
        else if (step_dir == 2){
            nextX = X;
            nextY = Y + len;
        }
        else if (step_dir == 3){
            nextX = X + len;
            nextY = Y;
        }
    }
    else if (orientation == 8){
        if (step_dir == 0){
            nextX = X;
            nextY = Y - len;
        }
        else if (step_dir == 1){
            nextX = X + L_step;
            nextY = Y;
        }
        else if (step_dir == 2){
            nextX = X;
            nextY = Y - L_step;
        }
        else if (step_dir == 3){
            nextX = X + len;
            nextY = Y;
        }
    }
}
//Performs step if possible, and updates the nextX&Y with where the particle 
//should be after a step. This is then processed to see if the particle CAN make
//the step or not.

int checkHit(std::vector<std::vector<int>> num_grid, double X, double Y,
              int orientation, int len){
    int A = num_grid[X][Y];
    if (orientation == 0){
        if ((num_grid[X+1][Y] != 0) && (num_grid[X+1][Y] != A)){
            return num_grid[X+1][Y];
        }
        else if ((num_grid[X-1][Y] != 0) && (num_grid[X-1][Y] != A)){
            return num_grid[X-1][Y];
        }
        else if ((num_grid[X][0] != 0) && (num_grid[X][0] != A)){
            return num_grid[X][0];
        }
        else if ((num_grid[X][Y-1] != 0) && (num_grid[X][Y-1] != A)){
            return num_grid[X][Y-1];
        }
        else {
            return 0;
        }
    }
    else if (orientation == 1){
        if ((num_grid[0][Y] != 0) && (num_grid[0][Y] != A)){
            return num_grid[0][Y];
        }
        else if ((num_grid[X-1][Y] != 0) && (num_grid[X-1][Y] != A)){
            return num_grid[X-1][Y];
        }
        else if ((num_grid[X][Y+1] != 0) && (num_grid[X][Y+1] != A)){
            return num_grid[X][Y+1];
        }
        else if ((num_grid[X][Y-1] != 0) && (num_grid[X][Y-1] != A)){
            return num_grid[X][Y-1];
        }
        else {
            return 0;
        }
    } 
    else if (orientation == 2){
        if ((num_grid[X+1][Y] != 0) && (num_grid[X+1][Y] != A)){
            return num_grid[X+1][Y];
        }
        else if ((num_grid[X-1][Y] != 0) && (num_grid[X-1][Y] != A)){
            return num_grid[X-1][Y];
        }
        else if ((num_grid[X][Y+1] != 0) && (num_grid[X][Y+1] != A)){
            return num_grid[X][Y+1];
        }
        else if ((num_grid[X][len] != 0) && (num_grid[X][len] != A)){
            return num_grid[X][len];
        }
        else {
            return 0;
        }
    }
    else if (orientation == 3){
        if ((num_grid[X+1][Y] != 0) && (num_grid[X+1][Y] != A)){
            return num_grid[X+1][Y];
        }
        else if ((num_grid[len][Y] != 0) && (num_grid[len][Y] != A)){
            return num_grid[len][Y];
        }
        else if ((num_grid[X][Y+1] != 0) && (num_grid[X][Y+1] != A)){
            return num_grid[X][Y+1];
        }
        else if ((num_grid[X][Y-1] != 0) && (num_grid[X][Y-1] != A)){
            return num_grid[X][Y-1];
        }
        else {
            return 0;
        }
    }
    else if (orientation == 4){
        if ((num_grid[X+1][Y] != 0) && (num_grid[X+1][Y] != A)){
            return num_grid[X+1][Y];
        }
        else if ((num_grid[X-1][Y] != 0) && (num_grid[X-1][Y] != A)){
            return num_grid[X-1][Y];
        }
        else if ((num_grid[X][Y+1] != 0) && (num_grid[X][Y+1] != A)){
            return num_grid[X][Y+1];
        }
        else if ((num_grid[X][Y-1] != 0) && (num_grid[X][Y-1] != A)){
            return num_grid[X][Y-1];
        }
        else {
            return 0;
        }
    }
    else if (orientation == 5){
        if ((num_grid[0][Y] != 0) && (num_grid[0][Y] != A)){
            return num_grid[0][Y];
        }
        else if ((num_grid[X-1][Y] != 0) && (num_grid[X-1][Y] != A)){
            return num_grid[X-1][Y];
        }
        else if ((num_grid[X][0] != 0) && (num_grid[X][0] != A)){
            return num_grid[X][0];
        }
        else if ((num_grid[X][Y-1] != 0) && (num_grid[X][Y-1] != A)){
            return num_grid[X][Y-1];
        }
        else {
            return 0;
        }
    }
    else if (orientation == 6){
        if ((num_grid[0][Y] != 0) && (num_grid[0][Y] != A)){
            return num_grid[0][Y];
        }
        else if ((num_grid[X-1][Y] != 0) && (num_grid[X-1][Y] != A)){
            return num_grid[X-1][Y];
        }
        else if ((num_grid[X][Y+1] != 0) && (num_grid[X][Y+1] != A)){
            return num_grid[X][Y+1];
        }
        else if ((num_grid[X][len] != 0) && (num_grid[X][len] != A)){
            return num_grid[X][len];
        }
        else {
            return 0;
        }
    }
    else if (orientation == 7){
        if ((num_grid[X+1][Y] != 0) && (num_grid[X+1][Y] != A)){
            return num_grid[X+1][Y];
        }
        else if ((num_grid[len][Y] != 0) && (num_grid[len][Y] != A)){
            return num_grid[len][Y];
        }
        else if ((num_grid[X][Y+1] != 0) && (num_grid[X][Y+1] != A)){
            return num_grid[X][Y+1];
        }
        else if ((num_grid[X][len] != 0) && (num_grid[X][len] != A)){
            return num_grid[X][len];
        }
        else {
            return 0;
        }
    }
    else if (orientation == 8){
        if ((num_grid[X+1][Y] != 0) && (num_grid[X+1][Y] != A)){
            return num_grid[X+1][Y];
        }
        else if ((num_grid[len][Y] != 0) && (num_grid[len][Y] != A)){
            return num_grid[len][Y];
        }
        else if ((num_grid[X][0] != 0) && (num_grid[X][0] != A)){
            return num_grid[X][0];
        }
        else if ((num_grid[X][Y-1] != 0) && (num_grid[X][Y-1] != A)){
            return num_grid[X][Y-1];
        }
        else {
            return 0;
        }
    }
    else {
        std::cout << "Exception in checkHitInit!" << std::endl;
        return 0;
    }
}
//Checks if a cell has non-zero neighbours. Currently only nearest neighbours.

void joinClusters(std::vector<std::vector<std::complex<double>>> &clusters,
                  std::vector<std::vector<int>> &num_grid, int A, int B){
    int C = std::min(A,B);
    int D = std::max(A,B);
    for (unsigned int i = 0; i < clusters[D].size(); ++i){
        clusters[C].push_back(clusters[D][i]);
        num_grid[floor(real(clusters[D][i]))][floor(imag(clusters[D][i]))] = C;
    }
    clusters.erase(clusters.begin() + D);
}
//Joins the clusters labeled A & B together, so that the lowest ID'd cluster
//gains all the other particles, and the highest is erased.

void joinInit(std::vector<std::vector<std::complex<double>>> &clusters,
              std::vector<std::vector<int>> &num_grid){
    double X,Y;
    int len = num_grid.size()-1;
    int A, B, orientation;
    B = 0;
    for (int i = 1; i < clusters.size(); ++i){
        for (int j = 0; j < clusters[i].size(); ++j){
            X = real(clusters[i][j]);
            Y = imag(clusters[i][j]);
            A = num_grid[floor(X)][floor(Y)];
            orientation = checkOrientation(X, Y, len);
            B = checkHit(num_grid, X, Y, orientation, len);
            if (B != 0){
                joinClusters(clusters, num_grid, A, B);
                updateNum_grid(num_grid, std::max(A, B), clusters);
            }
        }
    }
}
//Checks all particles to begin with, and joins them together in clusters if
//they are nearest neighbours.

void updateNum_grid(std::vector<std::vector<int>> &num_grid, int start,
                    std::vector<std::vector<std::complex<double>>> clusters){
    for (int i = start; i < clusters.size(); ++i){
        for (int j = 0; j < clusters[i].size(); ++j){
            num_grid[floor(real(clusters[i][j]))][floor(imag(clusters[i][j]))]--;
        }
    }
}
//Updates the num_grid after two clusters join, so all later clusters will
//decrease their ID by 1.

double clusterSize(std::vector<std::vector<std::complex<double>>> clusters,
                 int s, int system_length){
    int N = 0;
    for (int i = 0; i < clusters.size(); ++i){
        if (clusters[i].size() == s){
            N++;
        }
    }
    double n = double(N)/std::pow(double(system_length), 2.0);
    return n;
}

void clustersTime(std::vector<std::vector<std::complex<double>>> clusters,
                  int system_length){
    std::ofstream out_stream;
	out_stream.open("cluster_size_distribution(time).txt");
    int max = 0;
    for (int i = 0; i < clusters.size(); ++i){
        std::cout << clusters[i].size() << std::endl;
        if (max < clusters[i].size()){
            max = clusters[i].size();
        }
    }
    int N;
    for (int s = 0; s < max+1; ++s){
        N = 0;
        for (int i = 0; i < clusters.size(); ++i){
            if (clusters[i].size() == s){
                N++;
            }
        }
        out_stream << s << " " << double(N)/std::pow(double(system_length), 2.0)
                   << std::endl;
    }
    out_stream.close( );
}











void printClusters(std::vector<std::vector<std::complex<double>>> &clusters){
    for (unsigned int i = 0; i < clusters.size(); ++i){
        for (unsigned int j = 0; j < clusters[i].size(); ++j){
            std::cout << clusters[i][j] << std::endl;
        }
    std::cout << std::endl;
    }
}

void printGrid(std::vector<std::vector<std::complex<double>>> grid){
// Like the title suggests, this function prints whatever grid you would like to
// the screen. 
    for (int i = grid.size()-1; i >= 0; --i){
        for (int j = 0; j < grid.size(); ++j){
            if ((real(grid[j][i]) != 0) || (imag(grid[j][i]) != 0)){
                std::cout << 1;
            }
            else {
                std::cout << 0;
            }
        }
    std::cout << std::endl;
    }
    std::cout << std::endl;
}

void writeConfig(std::vector<std::vector<std::complex<double>>> clusters,
                 float r_p){
    std::ofstream out_stream;
    out_stream.open("data/before.txt");
    for (unsigned int i = 1; i < clusters.size(); ++i){
        for (unsigned int j = 0; j < clusters[i].size(); ++j){
            out_stream << real(clusters[i][j]) << " "
                       << imag(clusters[i][j]) << " " << r_p << std::endl;
        }
    }
    out_stream.close( );
}

void writeConfig1(std::vector<std::vector<std::complex<double>>> clusters,
                 float r_p, int system_length){
    std::ofstream out_stream;
    out_stream.open("data/after.txt");
    for (unsigned int i = 1; i < clusters.size(); ++i){
        for (unsigned int j = 0; j < clusters[i].size(); ++j){
            out_stream << std::fmod(real(clusters[i][j]), system_length) << " "
                       << std::fmod(imag(clusters[i][j]), system_length) << " "
                       << r_p << std::endl;
        }
    }
    out_stream.close( );
}

void writeConfigTest(std::vector<std::vector<std::complex<double>>> clusters,
                 float r_p, int system_length, double iteration){
    std::ofstream out_stream;
    std::string name = "data/" + std::to_string(iteration) + ".txt";
    out_stream.open(name);
    for (unsigned int i = 1; i < clusters.size(); ++i){
        for (unsigned int j = 0; j < clusters[i].size(); ++j){
            out_stream << std::fmod(real(clusters[i][j]), system_length) << " "
                       << std::fmod(imag(clusters[i][j]), system_length) << " "
                       << r_p << std::endl;
        }
    }
    out_stream.close( );
}

void writeConfigColor(std::vector<std::vector<std::complex<double>>> clusters,
                 float r_p, int system_length, double iteration){
    std::ofstream out_stream;
    std::string name = "data/00" + std::to_string(iteration) + ".txt";
    out_stream.open(name);
    for (unsigned int i = 1; i < clusters.size(); ++i){
        for (unsigned int j = 0; j < clusters[i].size(); ++j){
            out_stream << std::fmod(real(clusters[i][j]), system_length) << " "
                       << std::fmod(imag(clusters[i][j]), system_length) << " "
                       << r_p << " " << clusters[i].size() << std::endl;
        }
    }
    out_stream.close( );
}

void plotConfig(int system_length, std::string n){
    std::ofstream out_stream;
	out_stream.open("gnuplotter.gnu");
    out_stream << "set terminal png size 1200,1000"
                  " enhanced font \"Helvetica,12\" " << std::endl;
	out_stream << "set output \"fig/00" << n << ".png\" "
               << std::endl;
    out_stream << "set xrange [0:" << system_length-1 << "]" << std::endl;
    out_stream << "set yrange [0:" << system_length-1 << "]" << std::endl;
    out_stream << "set style fill transparent solid 1.0 noborder" << std::endl;
    out_stream << "unset key" << std::endl;
	out_stream << "plot \"data/00" << n << ".txt\" u 1:2:3:4 w circles palette"
               << std::endl;
	system("gnuplot gnuplotter.gnu");
}

void plot(int system_length, std::string n){
    std::ofstream out_stream;
	out_stream.open("gnuplotter.gnu");
    out_stream << "set terminal png size 1200,1000"
                  " enhanced font \"Helvetica,12\" " << std::endl;
	out_stream << "set output \"fig/" << n << ".png\" "
               << std::endl;
    out_stream << "set xrange [0:" << system_length-1 << "]" << std::endl;
    out_stream << "set yrange [0:" << system_length-1 << "]" << std::endl;
    out_stream << "unset key" << std::endl;
    out_stream << "set style fill transparent solid 1.0 noborder" << std::endl;
	out_stream << "plot \"data/" << n <<".txt\" with circles fc rgb \"navy\" "
               << std::endl;
	system("gnuplot gnuplotter.gnu");
    std::cout << "plot is OK!" << std::endl;
}

void plotConfigTest(int system_length, int i, int j){
    std::ofstream out_stream;
	out_stream.open("gnuplotter.gnu");
    out_stream << "set terminal png size 1200,1000"
                  " enhanced font \"Helvetica,12\" " << std::endl;
	out_stream << "set output \"fig/" << i << "-" << j << ".png\" "
               << std::endl;
    out_stream << "set xrange [0:" << system_length-1 << "]" << std::endl;
    out_stream << "set yrange [0:" << system_length-1 << "]" << std::endl;
    out_stream << "set style fill transparent solid 0.5 noborder" << std::endl;
//    out_stream << "unset key; unset tics; unset border" << std::endl;
	out_stream << "plot \"data/" << i << ".txt\" with circles fc rgb \"red\", \
                   \"data/" << j << ".txt\" with circles fc rgb \"navy\" "
               << std::endl;
	system("gnuplot gnuplotter.gnu");
    std::cout << "plotConfig is OK!" << std::endl;
}


void printNumGrid(std::vector<std::vector<int>> grid){
    for (int i = grid.size()-1; i >= 0; --i){
        for (int j = 0; j < grid.size(); ++j){
                std::cout << grid[j][i] << " ";
            }
        std::cout << std::endl; 
    }
    std::cout << std::endl;
}


void makeStepOrg(std::vector<std::vector<int>> num_grid, double X,
              double Y, int orientation, int step_dir, int len, float L_step,
              double &nextX, double &nextY){ 
    if (orientation == 0){
        // we must now chech if we can move to the cell in Step_dir direction.
        if (step_dir == 0){
            nextX = X;
            nextY = Y - len;
        }
        else if (step_dir == 1){
            nextX = X + L_step;
            nextY = Y;
        }
        else if (step_dir == 2){
            nextX = X;
            nextY = Y - L_step;
        }
        else if (step_dir == 3){
            nextX = X - L_step;
            nextY = Y;
        }
    }
    else if (orientation == 1){
        if (step_dir == 0){
            nextX = X;
            nextY = Y + L_step;
        }
        else if (step_dir == 1){
            nextX = X - len;
            nextY = Y;
        }
        else if (step_dir == 2){
            nextX = X;
            nextY = Y - L_step;
        }
        else if (step_dir == 3){
            nextX = X - L_step;
            nextY = Y;
        }
    }
    else if (orientation == 2){
        if (step_dir == 0){
            nextX = X;
            nextY = Y + L_step;
        }
        else if (step_dir == 1){
            nextX = X + L_step;
            nextY = Y;
        }
        else if (step_dir == 2){
            nextX = X;
            nextY = Y + len;
        }
        else if (step_dir == 3){
            nextX = X - L_step;
            nextY = Y;
        }
    }
    else if (orientation == 3){
        if (step_dir == 0){
            nextX = X;
            nextY = Y + L_step;
        }
        else if (step_dir == 1){
            nextX = X + L_step;
            nextY = Y;
        }
        else if (step_dir == 2){
            nextX = X;
            nextY = Y - L_step;
        }
        else if (step_dir == 3){
            nextX = X + len;
            nextY = Y;
        }
    }
    else if (orientation == 4){
        if (step_dir == 0){
            nextX = X;
            nextY = Y + L_step;
        }
        else if (step_dir == 1){
            nextX = X + L_step;
            nextY = Y;
        }
        else if (step_dir == 2){
            nextX = X;
            nextY = Y - L_step;
        }
        else if (step_dir == 3){
            nextX = X - L_step;
            nextY = Y;
        }
    }
    else if (orientation == 5){
        if (step_dir == 0){
            nextX = X;
            nextY = Y - len;
        }
        else if (step_dir == 1){
            nextX = X - len;
            nextY = Y;
        }
        else if (step_dir == 2){
            nextX = X;
            nextY = Y - L_step;
        }
        else if (step_dir == 3){
            nextX = X - L_step;
            nextY = Y;
        }
    }
    else if (orientation == 6){
        if (step_dir == 0){
            nextX = X;
            nextY = Y + L_step;
        }
        else if (step_dir == 1){
            nextX = X - len;
            nextY = Y;
        }
        else if (step_dir == 2){
            nextX = X;
            nextY = Y + len;
        }
        else if (step_dir == 3){
            nextX = X - L_step;
            nextY = Y;
        }
    }
    else if (orientation == 7){
        if (step_dir == 0){
            nextX = X;
            nextY = Y + L_step;
        }
        else if (step_dir == 1){
            nextX = X + L_step;
            nextY = Y;
        }
        else if (step_dir == 2){
            nextX = X;
            nextY = Y + len;
        }
        else if (step_dir == 3){
            nextX = X + len;
            nextY = Y;
        }
    }
    else if (orientation == 8){
        if (step_dir == 0){
            nextX = X;
            nextY = Y - len;
        }
        else if (step_dir == 1){
            nextX = X + L_step;
            nextY = Y;
        }
        else if (step_dir == 2){
            nextX = X;
            nextY = Y - L_step;
        }
        else if (step_dir == 3){
            nextX = X + len;
            nextY = Y;
        }
    }
    if (isFree(num_grid, nextX, nextY, X, Y) == false){
        nextX = X;
        nextY = Y;
    }
}
