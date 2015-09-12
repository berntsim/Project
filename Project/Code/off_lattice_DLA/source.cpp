#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

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


void testDist(std::vector<std::vector<int>> &on_lattice_distance,
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

void printGrid(std::vector<std::vector<int>> &grid, int sys_len){
    for (unsigned int i = 0; i < sys_len; ++i){
        for (unsigned int j = 0; j < sys_len; ++j){
            std::cout << grid[i][j];
        }
    std::cout << std::endl;
    }
    std::cout << std::endl;
}
