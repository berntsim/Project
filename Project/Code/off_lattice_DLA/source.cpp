#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
//#include <cmath>
#include <random>
#include <algorithm>

#include "header.h"

void initEx_pos(std::vector<float> &x_pos, std::vector<float> &y_pos,
                int sys_len){
    x_pos.push_back(float(sys_len)/2.0);
    y_pos.push_back(float(sys_len)/2.0);
}


void initOnLatticeCluster(std::vector<std::vector<int>> &on_lattice_cluster, 
                                int system_length){
    for (unsigned int i = 0; i < system_length; ++i){
        for (unsigned int j = 0; j < system_length; ++j){
            if ((i == system_length/2) && (j == system_length/2)){
                on_lattice_cluster[i][j] = 1;
            }
            else{
                on_lattice_cluster[i][j] = 0;
            }
        }
    }
}

void initOnLatticeDistance(std::vector<std::vector<int>> &on_lattice_distance,
                           int sys_len, int D_max){
    for (unsigned int i = 0; i < sys_len; ++i){
        for (unsigned int j = 0; j < sys_len; ++j){
            if ((i == sys_len/2) && (j == sys_len/2)){
                on_lattice_distance[i][j] = 0;
            }
            else {
                on_lattice_distance[i][j] = D_max;
            }
        }
    }
}


void neighbourDist(std::vector<std::vector<int>> &on_lattice_distance,
                   int D_max, int row, int col, int sys_len){
    if ((row - D_max < 0) && (col - D_max)){ 
        int min_row = std::min(row,D_max);
        int min_col = std::min(col,D_max);
        for (int d = 0; d <= D_max; ++d){
            for (int n = 0; n <= D_max; ++n){
                if(on_lattice_distance[row][col+d] > d){ 
                    on_lattice_distance[row][col+d] = d;
                }
                if (on_lattice_distance[row+n][col] > n){
                    on_lattice_distance[row+n][col] = n;
                }
                if ((floor(sqrt(float(n)*float(n)+float(d)*float(d)))<D_max) &&
                    (on_lattice_distance[row+n][col+d] > 
                     floor(sqrt(float(n)*float(n)+float(d)*float(d))))){
                    on_lattice_distance[row+n][col+d] = 
                    floor(sqrt(float(n)*float(n)+float(d)*float(d)));
                }
            }
        }
        for (int d = 0; d <= D_max; ++d){
            for (int n = 0; n <= min_row; ++n){
                if(on_lattice_distance[row][col+d] > d){ 
                    on_lattice_distance[row][col+d] = d;
                }
                if ((floor(sqrt(float(n)*float(n)+float(d)*float(d)))<D_max) &&
                    (on_lattice_distance[row-n][col+d] > 
                     floor(sqrt(float(n)*float(n)+float(d)*float(d))))){
                    on_lattice_distance[row-n][col+d] =
                    floor(sqrt(float(n)*float(n)+float(d)*float(d)));
                }
            }
        }
        for (int d = 0; d <= min_col; ++d){
            for (int n = 0; n <= D_max; ++n){
                if (on_lattice_distance[row+n][col] > n){
                    on_lattice_distance[row+n][col] = n;
                }
                if ((floor(sqrt(float(n)*float(n)+float(d)*float(d)))<D_max) &&
                    (on_lattice_distance[row+n][col-d] > 
                     floor(sqrt(float(n)*float(n)+float(d)*float(d))))) {
                    on_lattice_distance[row+n][col-d] =
                    floor(sqrt(float(n)*float(n)+float(d)*float(d)));
                }
            }
        }
        for (int d = 0; d <= min_col; ++d){
            for (int n = 0; n <= min_row; ++n){
                if ((floor(sqrt(float(n)*float(n)+float(d)*float(d)))<D_max) &&
                    (on_lattice_distance[row-n][col-d] > 
                     floor(sqrt(float(n)*float(n)+float(d)*float(d))))) {
                    on_lattice_distance[row-n][col-d] =
                    floor(sqrt(float(n)*float(n)+float(d)*float(d)));
                }
            }
        }
    }
    else if ((row + D_max > sys_len) && (col + D_max > sys_len)){
    int min_col = sys_len - col;
    int min_row = sys_len - row;
        for (int d = 0; d < D_max; ++d) {
            for (int n = 0; n < D_max; ++n){
                if (on_lattice_distance[row][col-d] > d){
                    on_lattice_distance[row][col-d] = d;
                }
                if(on_lattice_distance[row-n][col] > n) {
                    on_lattice_distance[row-n][col] = n;
                }
                if ((floor(sqrt(float(n)*float(n)+float(d)*float(d)))<D_max) &&
                    (on_lattice_distance[row-n][col-d] > 
                     floor(sqrt(float(n)*float(n)+float(d)*float(d))))){
                    on_lattice_distance[row-n][col-d] =
                    floor(sqrt(float(n)*float(n)+float(d)*float(d)));
                }
            }
        }
        for (int d = 0; d < D_max; ++d){
            for (int n = 0; n < min_row; ++n){
                if (on_lattice_distance[row][col-d] > d){
                    on_lattice_distance[row][col-d] = d;
                }
                if ((floor(sqrt(float(n)*float(n)+float(d)*float(d)))<D_max) &&
                    (on_lattice_distance[row+n][col-d] > 
                     floor(sqrt(float(n)*float(n)+float(d)*float(d))))) {
                    on_lattice_distance[row+n][col-d] =
                    floor(sqrt(float(n)*float(n)+float(d)*float(d)));
                }
            }
        }
        for (int d = 0; d < min_col; ++d){
            for (int n = 0; n < D_max; ++n){
                if(on_lattice_distance[row-n][col] > n) {
                    on_lattice_distance[row-n][col] = n;
                }
                if ((floor(sqrt(float(n)*float(n)+float(d)*float(d)))<D_max) &&
                    (on_lattice_distance[row-n][col+d] > 
                     floor(sqrt(float(n)*float(n)+float(d)*float(d))))){
                    on_lattice_distance[row-n][col+d] =
                    floor(sqrt(float(n)*float(n)+float(d)*float(d)));
                }
            }
        }
        for (int d = 0; d < min_col; ++d){
            for (int n = 0; n < min_row; ++n) {
                if ((floor(sqrt(float(n)*float(n)+float(d)*float(d)))<D_max) &&
                    (on_lattice_distance[row+n][col+d] > 
                     floor(sqrt(float(n)*float(n)+float(d)*float(d))))){
                    on_lattice_distance[row+n][col+d] =
                    floor(sqrt(float(n)*float(n)+float(d)*float(d)));
                }
            }
        }
    }
    else if ((row + D_max > sys_len) && (col - D_max < 0)){
    int min_col = std::min(D_max,col);
    int min_row = sys_len-row;
    std::cout << "min_row = " << min_row << std::endl;
    std::cout << "min_col = " << min_col << std::endl;
    for (int d = 0; d < D_max; ++d){
        for (int n = 0; n < D_max; ++n){
            if(on_lattice_distance[row][col+d] > d){ 
                    on_lattice_distance[row][col+d] = d;
            }
            if(on_lattice_distance[row-n][col] > n) {
                    on_lattice_distance[row-n][col] = n;
            }
            if ((floor(sqrt(float(n)*float(n)+float(d)*float(d)))<D_max) &&
                    (on_lattice_distance[row-n][col+d] > 
                     floor(sqrt(float(n)*float(n)+float(d)*float(d))))){
                    on_lattice_distance[row-n][col+d] =
                    floor(sqrt(float(n)*float(n)+float(d)*float(d)));
            }
        }
    }
    for (int d = 0; d < D_max; ++d) {
        for (int n = 0; n < min_row; ++n) {
            if(on_lattice_distance[row+n][col] > n) {
                    on_lattice_distance[row+n][col] = n;
            }
            if ((floor(sqrt(float(n)*float(n)+float(d)*float(d)))<D_max) &&
                    (on_lattice_distance[row+n][col+d] > 
                     floor(sqrt(float(n)*float(n)+float(d)*float(d))))){
                    on_lattice_distance[row+n][col+d] = 
                    floor(sqrt(float(n)*float(n)+float(d)*float(d)));
            }
        }
    }
    for (int d = 0; d <= min_col; ++d){
        for (int n = 0; n < min_row; ++n){
            if (on_lattice_distance[row][col-d] > d){
                    on_lattice_distance[row][col-d] = d;
            }
            if (on_lattice_distance[row+n][col] > n){
                    on_lattice_distance[row+n][col] = n;
            }
            if ((floor(sqrt(float(n)*float(n)+float(d)*float(d)))<D_max) &&
                    (on_lattice_distance[row+n][col-d] > 
                     floor(sqrt(float(n)*float(n)+float(d)*float(d))))){
                    on_lattice_distance[row+n][col-d] =
                    floor(sqrt(float(n)*float(n)+float(d)*float(d)));
            }
        }
    }
    for (int d = 0; d <= min_col; ++d){
        for (int n = 0; n < D_max;++n){
            if (on_lattice_distance[row][col-d] > d){
                    on_lattice_distance[row][col-d] = d;
            }
            if(on_lattice_distance[row-n][col] > n) {
                    on_lattice_distance[row-n][col] = n;
            }
            if ((floor(sqrt(float(n)*float(n)+float(d)*float(d)))<D_max) &&
                    (on_lattice_distance[row-n][col-d] > 
                     floor(sqrt(float(n)*float(n)+float(d)*float(d))))){
                    on_lattice_distance[row-n][col-d] =
                    floor(sqrt(float(n)*float(n)+float(d)*float(d)));
            }
        }
    }
    }
    else if ((row-D_max < 0) && (col + D_max > sys_len)){
        int min_row = std::min(row,D_max);
        int min_col = sys_len - col;
        for (int d = 0; d < D_max; ++d){
            for (int n = 0; n < D_max; ++n){
                if (on_lattice_distance[row][col-d] > d){
                    on_lattice_distance[row][col-d] = d;
                }
                if (on_lattice_distance[row+n][col] > n){
                    on_lattice_distance[row+n][col] = n;
                }
                if ((floor(sqrt(float(n)*float(n)+float(d)*float(d)))<D_max) &&
                    (on_lattice_distance[row+n][col-d] > 
                     floor(sqrt(float(n)*float(n)+float(d)*float(d))))) {
                    on_lattice_distance[row+n][col-d] =
                    floor(sqrt(float(n)*float(n)+float(d)*float(d)));
                }
            }
        }
        for (int d = 0; d < min_col; ++d){
            for (int n = 0; n < D_max; ++n){
                if(on_lattice_distance[row][col+d] > d){ 
                    on_lattice_distance[row][col+d] = d;
                }
                if (on_lattice_distance[row+n][col] > n){
                    on_lattice_distance[row+n][col] = n;
                }
                if ((floor(sqrt(float(n)*float(n)+float(d)*float(d)))<D_max) &&
                    (on_lattice_distance[row+n][col+d] > 
                     floor(sqrt(float(n)*float(n)+float(d)*float(d))))){
                    on_lattice_distance[row+n][col+d] = 
                    floor(sqrt(float(n)*float(n)+float(d)*float(d)));
                }
            }
        }
    }
    else{
    int min_col = D_max;
    int min_row = D_max;
        for ( int d = 0; d < min_col; ++d){
            for ( int n = 0; n < min_row; ++n){
                if(on_lattice_distance[row][col+d] > d){ 
                    on_lattice_distance[row][col+d] = d;
                }
                else if (on_lattice_distance[row][col-d] > d){
                    on_lattice_distance[row][col-d] = d;
                }
                else if (on_lattice_distance[row+n][col] > n){
                    on_lattice_distance[row+n][col] = n;
                }
                else if(on_lattice_distance[row-n][col] > n) {
                    on_lattice_distance[row-n][col] = n;
                }
                if ((floor(sqrt(float(n)*float(n)+float(d)*float(d)))<D_max) &&
                    (on_lattice_distance[row+n][col+d] > 
                     floor(sqrt(float(n)*float(n)+float(d)*float(d))))){
                    on_lattice_distance[row+n][col+d] = 
                    floor(sqrt(float(n)*float(n)+float(d)*float(d)));
                }
                if ((floor(sqrt(float(n)*float(n)+float(d)*float(d)))<D_max) &&
                    (on_lattice_distance[row-n][col-d] > 
                     floor(sqrt(float(n)*float(n)+float(d)*float(d))))){
                    on_lattice_distance[row-n][col-d] =
                    floor(sqrt(float(n)*float(n)+float(d)*float(d)));
                }
                if ((floor(sqrt(float(n)*float(n)+float(d)*float(d)))<D_max) &&
                    (on_lattice_distance[row+n][col-d] > 
                     floor(sqrt(float(n)*float(n)+float(d)*float(d))))) {
                    on_lattice_distance[row+n][col-d] =
                    floor(sqrt(float(n)*float(n)+float(d)*float(d)));
                }
                if ((floor(sqrt(float(n)*float(n)+float(d)*float(d)))<D_max) &&
                    (on_lattice_distance[row-n][col+d] > 
                     floor(sqrt(float(n)*float(n)+float(d)*float(d))))){
                    on_lattice_distance[row-n][col+d] =
                    floor(sqrt(float(n)*float(n)+float(d)*float(d)));
                }
            }
        }
    }
}

void printGrid(std::vector<std::vector<int>> &grid){
    for (int i = grid.size()-1; i >= 0; --i){
        for (int j = 0; j < grid.size(); ++j){
            std::cout << grid[i][j];
        }
    std::cout << std::endl;
    }
    std::cout << std::endl;
}

void initVicinity(std::vector<std::vector<int>> &vicinity, int D_max){
    for (unsigned int i = 0; i < (2*D_max + 1)/2+1; ++i){
        for (unsigned int j = 0; j < (2*D_max + 1)/2+1; ++j){
            vicinity[(2*D_max+1)/2+i][(2*D_max+1)/2+j] = floor(sqrt(j*j + i*i));
            vicinity[(2*D_max+1)/2-i][(2*D_max+1)/2+j] = floor(sqrt(j*j + i*i));
            vicinity[(2*D_max+1)/2+i][(2*D_max+1)/2-j] = floor(sqrt(j*j + i*i));
            vicinity[(2*D_max+1)/2-i][(2*D_max+1)/2-j] = floor(sqrt(j*j + i*i));
        }
    }
}

float findRadius(std::vector<float> &x_pos, std::vector<float> &y_pos,
                 float r_p){
    float rad  = 0;
    for (unsigned int i = 0; i < x_pos.size(); ++i){
        if (rad < sqrt(pow(x_pos[i] - x_pos[0],2) + pow(y_pos[i] - y_pos[0],2))
            +r_p){
            rad = sqrt(pow(x_pos[i] - x_pos[0],2) + pow(y_pos[i] - y_pos[0],2))
                  +r_p; 
        }
    }
    return rad;
}

void startWalker(float &walker_x_pos, float &walker_y_pos, int D_max,
                 std::vector<float> &x_pos, std::vector<float> &y_pos,
                 float theta, double pi, float r_p, float &r_c){
    float circ_rad = r_c + D_max;
    walker_x_pos = x_pos[0] + circ_rad*cos(theta);
    walker_y_pos = y_pos[0] + circ_rad*sin(theta);
//    walker_x_pos = circ_rad*cos(theta);
//    walker_y_pos = circ_rad*sin(theta);
}

float updateCluster(std::vector<std::vector<int>> &on_lattice_cluster,
                   std::vector<std::vector<int>> &on_lattice_distance,
                   std::vector<float> &x_pos, std::vector<float> &y_pos,
                   float walker_x_pos, float walker_y_pos, int D_max,
                   int sys_len, float &r_c){
    x_pos.push_back(walker_x_pos);
    y_pos.push_back(walker_y_pos);
    on_lattice_cluster[floor(walker_y_pos)][floor(walker_x_pos)] = x_pos.size();
    neighbourDist(on_lattice_distance, D_max, floor(walker_y_pos),
                  floor(walker_x_pos), sys_len);
    float tmp = std::sqrt(std::pow(walker_x_pos - x_pos[0],2) + 
                          std::pow(walker_y_pos - y_pos[0],2));
    if (tmp > r_c){
        r_c = tmp;
    } 
    return r_c;
}

float findD_wc(float walker_x_pos, float walker_y_pos, int D_max,
              std::vector<std::vector<int>> &on_lattice_distance){
    int len = on_lattice_distance.size();
    if ((walker_x_pos < len-1) && (walker_x_pos > 0) && (walker_y_pos < len-1)
        && (walker_y_pos > 0)) {
        return on_lattice_distance[round(walker_y_pos)][round(walker_x_pos)];
    }
    else {
        return D_max;
    }
}

float LHit(float step_L, float step_dir, float x_c, float y_c, float x_p,
               float y_p, float r_p){
    float d_p = 2*r_p;
    float a = 1.0;
    float b = 2*(std::cos(step_dir)*(x_p-x_c) + std::sin(step_dir)*(y_p-y_c));
    float c = (x_c-x_p)*(x_c-x_p) + (y_c-y_p)*(y_c-y_p) - d_p*d_p;
    float res[2];
    bool sol = true;
    if ((b*b - 4*a*c) < 0){
        sol = false;
        return step_L;
    }
    else {
        res[0] = (-b + sqrt(b*b - 4*a*c))/(2*a);
        res[1] = (-b - sqrt(b*b - 4*a*c))/(2*a);
        if(std::abs(res[0]) < std::abs(res[1])){
            if(res[0] < 0){
                return step_L;
            }
            else if(std::abs(res[0]) < step_L){
                return res[0];
            }
            else {
                return step_L;
            }
        }
        else{
            if(res[1] < 0) {
                return step_L;
            }
            else if(std::abs(res[1]) < step_L){
                return res[1];
            }
            else {
                return step_L;
            }
        }
    }
}

bool killParticle(float walker_x_pos, float walker_y_pos,
                  std::vector<float> &x_pos, std::vector<float> &y_pos,
                  float r_p, float &r_c){
    if (sqrt(pow(walker_x_pos-x_pos[0],2) + pow(walker_y_pos-y_pos[0],2)) >= 
             3.0*(r_c+1)){ // 3 should be 5 i paper
    return true;
    }
    else {
        return false;
    }
}

bool createStep(float &step_L, float step_dir,float &walker_x_pos,
                float &walker_y_pos, double pi, float r_p,
                std::vector<std::vector<int>> &on_lattice_distance,
                std::vector<std::vector<int>> &on_lattice_cluster,
                float L_min, int D_max, std::vector<float> &x_pos,
                std::vector<float> &y_pos, int sys_len, float &r_c){
    int d_wc = findD_wc(walker_x_pos, walker_y_pos, D_max, on_lattice_distance);
    float len = 2*r_p + L_min + 1;
    std::vector<float> min_L;
    float step_org = step_L;
    std::vector<int> labels;
    if (d_wc <= len) {
        for (int i=round(walker_y_pos)-len; i<=round(walker_y_pos)+len;++i){
            for (int j=round(walker_x_pos)-len; j<=round(walker_x_pos)+len;++j){
                if (on_lattice_cluster[i][j] != 0){
                    labels.push_back(on_lattice_cluster[i][j]);

                }
            }
        }
        for (unsigned int i = 0; i < labels.size(); ++i){
            if (LHit(step_L, step_dir, x_pos[labels[i]-1], y_pos[labels[i]-1],
                     walker_x_pos, walker_y_pos, r_p) < step_org){
                min_L.push_back(LHit(step_L, step_dir, x_pos[labels[i]-1],
                         y_pos[labels[i]-1], walker_x_pos, walker_y_pos, r_p));
            }
        }
        for (unsigned int i = 0; i < min_L.size(); ++i){
            if (step_L >= min_L[i]){
                step_L = min_L[i];
            }
        }

        if (step_L < step_org){
            walker_x_pos += step_L*std::cos(step_dir);
            walker_y_pos += step_L*std::sin(step_dir);
            r_c = updateCluster(on_lattice_cluster, on_lattice_distance,
                          x_pos, y_pos, walker_x_pos,walker_y_pos, D_max,
                          sys_len, r_c);
            return true;
        }
        else {
            walker_x_pos += step_L*std::cos(step_dir);
            walker_y_pos += step_L*std::sin(step_dir);
            return false;
        }
    }
    else if ((len < d_wc) && (d_wc < D_max)){
        step_L = d_wc - len;
        walker_x_pos += step_L*std::cos(step_dir);
        walker_y_pos += step_L*std::sin(step_dir);
        return false;
    }
    else {
        //for (unsigned int i = 0; i < x_pos.size(); ++i){
        //}
        step_L = std::max(float(sqrt(pow(walker_x_pos-x_pos[0],2) + 
                          pow(walker_y_pos-y_pos[0],2)+r_p) -
                          r_c - len),
                          float(D_max - len));

        walker_x_pos += step_L*std::cos(step_dir);
        walker_y_pos += step_L*std::sin(step_dir);
        return false;
    }
}

void runAll(float &step_L, float &step_dir, std::mt19937::result_type seed,
                float &walker_x_pos, float &walker_y_pos, double pi, float r_p,
                std::vector<std::vector<int>> &on_lattice_distance,
                std::vector<std::vector<int>> &on_lattice_cluster,
                float L_min, int D_max, std::vector<float> &x_pos,
                std::vector<float> &y_pos, int sys_len, int nbr_particles,
                unsigned int counter, float &r_c){

    auto rand_dir = std::bind(std::uniform_real_distribution<float>(0,2*pi),
				   std::mt19937(seed));
    bool hit = false;
    bool kill = false;
    while (x_pos.size() < nbr_particles){
        startWalker(walker_x_pos, walker_y_pos, D_max, x_pos, y_pos, rand_dir(),
                    pi, r_p, r_c);
        bool hit = false;
        bool kill = false;
        while ((hit == false) && (kill == false)){
            hit = createStep(step_L, rand_dir(), walker_x_pos, walker_y_pos, pi,
                             r_p, on_lattice_distance, on_lattice_cluster,
                             L_min, D_max, x_pos, y_pos, sys_len, r_c);
            kill = killParticle(walker_x_pos, walker_y_pos, x_pos, y_pos, r_p,
                                r_c);
            ++counter;
        }
    }
    std::cout << "counter = " << counter << std::endl;
}





void writeGrid(std::vector<float> &x_pos, std::vector<float> &y_pos, float r_p){
    std::ofstream out_stream;
    out_stream.open("data/dla.txt");
    for (unsigned int i = 0; i < x_pos.size(); ++i){
        out_stream << x_pos[i] << " " << y_pos[i] << " " << r_p << std::endl;
    }
    out_stream.close( );
}

void plotGnuplot(int &arr_len){
    std::ofstream out_stream;
	out_stream.open("gnuplotter.gnu");
    out_stream << "set terminal png size 1200,1000"
                  " enhanced font \"Helvetica,12\" " << std::endl;
	out_stream << "set output \"figures/dla.png\" " << std::endl;
	out_stream << "set xlabel \"x-dim\" " << std::endl;
	out_stream << "set ylabel \"y-dim\"" << std::endl;
	out_stream << "set xrange [0:" << arr_len << "]" << std::endl;
	out_stream << "set yrange [0:" << arr_len << "]" << std::endl;
	out_stream << "plot \"data/dla.txt\" with circles fc rgb \"navy\" "
                  "fill solid " << std::endl;
	system("gnuplot gnuplotter.gnu");
}


float findFractalDim(std::vector<float> x_pos, std::vector<float> y_pos,
                     float r_p, float &r_c){
    return log(x_pos.size())/log(r_c);
}

void printDistances(std::vector<float> &x_pos, std::vector<float> &y_pos,
                    float r_p){
    float d = 10;
    std::vector<float> shorter;
    for (unsigned int i = 0; i < x_pos.size()-1; ++i){
        for (unsigned int j = 0; j < x_pos.size()-1; ++j){
            d = sqrt(pow(x_pos[i+1]-x_pos[j],2) + pow(y_pos[i+1]-y_pos[j],2));
            if ((d < 2*r_p) && (d != 0)){
                shorter.push_back(d);
            }
        }
    }
    for (unsigned int i = 0; i < shorter.size(); ++i){
        std::cout << shorter[i] << std::endl;
    }
    std::cout << "len Shorter = " << shorter.size() << std::endl;
}

