#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <algorithm>
#include <math.h>

#include "header.h"

void initGrid(std::vector<std::vector<std::complex<double>>> &grid,
              int system_length, std::mt19937::result_type seed,
              int nbr_particles,
              std::vector<std::vector<std::complex<double>>> &clusters,
              std::vector<std::vector<int>> &num_grid){
    int counter = 1;
    double tmp1;
    double tmp2;
    for (unsigned int i = 0; i < grid.size(); ++i){
        for (unsigned int j = 0; j < grid.size(); ++j){
            grid[i][j] = {0,0};
            num_grid[i][j]=0;
        }
    }
    //auto coord = std::bind(std::uniform_real_distribution<float>
    //                       (0, double(grid.size()-1)), std::mt19937(seed));
    auto coord = std::bind(std::uniform_int_distribution<int>
                           (0, grid.size()-1), std::mt19937(seed));
    clusters[1][0] = {2,1};
    clusters[5][0] = {4,3};
    clusters[3][0] = {2,4};
    clusters[4][0] = {1,3};
    clusters[2][0] = {3,4};
    grid[2][1] = {2,1};
    grid[4][3] = {4,3};
    grid[2][4] = {2,4};
    grid[1][3] = {1,3};
    grid[3][4] = {3,4};
    num_grid[2][1] = 1;
    num_grid[4][3] = 5;
    num_grid[2][4] = 3;
    num_grid[1][3] = 4;
    num_grid[3][4] = 2;
//    while (counter < nbr_particles){
//        tmp1 = coord();
//        tmp2 = coord();
//        if ((real(grid[round(tmp1)][round(tmp2)]) == 0) && 
//            (imag(grid[round(tmp1)][round(tmp2)])) == 0){
//            clusters[counter][0] = {tmp1, tmp2};
//            grid[round(tmp1)][round(tmp2)] = {tmp1, tmp2};
//            num_grid[round(tmp1)][round(tmp2)] = counter;
//            ++counter;
//        }
//    }
}

void updateGrid(std::vector<std::vector<std::complex<double>>> clusters,
                std::vector<std::vector<std::complex<double>>> &grid){
    int counter = 0;
    int tmp1;
    int tmp2;
    for (unsigned int i = 0; i < grid.size(); ++i){
        for (unsigned int j = 0; j < grid.size(); ++j){
            grid[i][j] = 0;
        }
    }
    for (unsigned int i = 0; i < clusters.size(); ++i){
        for (unsigned int j = 1; j < clusters[i].size(); ++j){
            grid[int(floor(real(clusters[i][j])))]
                [int(floor(imag(clusters[i][j])))]
            = {clusters[i][j]};
        }
    }
}

void walkOnGrid(std::vector<std::vector<std::complex<double>>> &grid,
                std::vector<std::vector<std::complex<double>>> &clusters,
                std::vector<std::vector<int>> &num_grid,
                std::mt19937::result_type seed, float L_step, double pi,
                float r_p){
    int step_dir;
    int len = grid.size()-1;
    int X,Y;
    int reduced = 0;
    std::vector<std::vector<std::complex<double>>> tmp = clusters;
    auto rand_dir_int = std::bind(std::uniform_real_distribution<float>(0,4),
	                std::mt19937(seed));
    for (int i = 1; i < clusters.size(); ++i){
        //printGrid(grid);
        //printClusters(tmp);
//        step_dir = rand_dir_int();
        //std::cout << std::endl;
        step_dir = 1;
        std::cout << "len = " << clusters.size() << std::endl;
        std::cout << "We are in " << real(clusters[i][0]) << "," 
                  << imag(clusters[i][0]) << std::endl;
            std::cout << "We move in " << step_dir << " direction" << std::endl;
        for (int j = 0; j < clusters[i].size(); ++j){
            if (step_dir == 0){
                X = real(clusters[i][j]);
                Y = imag(clusters[i][j]);
                if (Y == len){
                    if ((real(grid[X][0]) == 0) && (imag(grid[X][0]) == 0)){
                        grid[X][0] = {real(clusters[i][j]),
                                      imag(clusters[i][j]) -(len)};
                        grid[X][Y] = {0,0};
                        clusters[i][j] = {real(clusters[i][j]),
                                         imag(clusters[i][j]) - (len)};
                        tmp[i][j] = {real(tmp[i-reduced][j]),
                                         imag(tmp[i-reduced][j]) - (len)};
                        num_grid[X][0] = num_grid[X][Y];
                        num_grid[X][Y]=0;
                        num_grid[X][0] = checkHit(tmp, grid, num_grid,
                                                  X, 0, r_p, reduced);
                    }
                    else {
                        std::cout << "Already something there" << std::endl;
                        //Hit!
                    }
                }
                else {
                    if ((real(grid[X][Y+1]) == 0) && (imag(grid[X][Y+1]) == 0)){
                        grid[X][Y+1] = {real(clusters[i][j]),
                                        imag(clusters[i][j]) + L_step};
                        grid[X][Y] = {0,0};
                        clusters[i][j] = {real(clusters[i][j]),
                                          imag(clusters[i][j]) + L_step};
                        tmp[i][j] = {real(tmp[i-reduced][j]),
                                          imag(tmp[i-reduced][j]) + L_step};
                        num_grid[X][Y+1] = num_grid[X][Y];
                        num_grid[X][Y]=0;
                        num_grid[X][Y+1] = checkHit(tmp, grid, num_grid,
                                                    X, Y+1, r_p, reduced);
                    }
                    else {
                        std::cout << "Already something there" << std::endl;
                        //Hit!
                    }
                }
            }
            else if (step_dir == 1){
                X = real(clusters[i][j]);
                Y = imag(clusters[i][j]);
                if (X == len){
                    if ((real(grid[0][Y]) == 0) && (imag(grid[0][Y]) == 0)){
                        grid[0][Y] = {real(clusters[i][j]) - (len),
                                      imag(clusters[i][j])};
                        grid[X][Y] = {0,0};
                        clusters[i][j] = {real(clusters[i][j]) - (len),
                                          imag(clusters[i][j])};
                        std::cout << "len(tmp) = " << tmp.size() << std::endl;
                        std::cout << "len(clusters) = " << clusters.size() << std::endl;
                        std::cout << "len = " << len << std::endl;
                        std::cout << "X = " << X << std::endl;
                        std::cout << "Y = " << Y << std::endl;
                        std::cout << " i = " << i << std::endl; 
                        std::cout << " j = " << j << std::endl; 
                        std::cout << "reduced = " << reduced << std::endl;
                        std::cout << real(tmp[i-reduced][j]) << std::endl;
                        std::cout << imag(tmp[i-reduced][j]) << std::endl;
                        tmp[i-reduced][j] = {real(tmp[i-reduced][j]) - len,
                                          imag(tmp[i-reduced][j])};
                        std::cout << "hei" << std::endl;
                        num_grid[0][Y] = num_grid[X][Y];
                        num_grid[X][Y]=0;
                        num_grid[0][Y] = checkHit(tmp, grid, num_grid,
                                                  0, Y, r_p, reduced);
                    }
                    else {
                        std::cout << "Already something there" << std::endl;
                        // Hit!
                    }
                }
                else {
                    if ((real(grid[X+1][Y]) == 0) && (imag(grid[X+1][Y]) == 0)){
                        grid[X+1][Y] = {real(clusters[i][j]) + L_step,
                                        imag(clusters[i][j])};
                        grid[X][Y] = {0,0};
                        clusters[i][j] = {real(clusters[i][j]) + L_step,
                                          imag(clusters[i][j])};
                        tmp[i][j] = {real(tmp[i-reduced][j]) + L_step,
                                          imag(tmp[i-reduced][j])};
                        num_grid[X+1][Y] = num_grid[X][Y];
                        num_grid[X][Y]=0;
                        num_grid[X+1][Y] = checkHit(tmp, grid, num_grid,
                                                    X+1, Y, r_p, reduced);
                    }
                    else {
                        std::cout << "Already something there" << std::endl;
                        //Hit!
                    }
                }
            }
            else if (step_dir == 2){
                X = real(clusters[i][j]);
                Y = imag(clusters[i][j]);
                if (Y == 0) {
                    if ((real(grid[X][len]) == 0) && (imag(grid[X][len]) == 0)){
                        grid[X][len] = {real(clusters[i][j]),
                                        imag(clusters[i][j]) + len};
                        grid[X][Y] = {0,0};
                        clusters[i][j] = {real(clusters[i][j]),
                                          imag(clusters[i][j]) + len};
                        tmp[i][j] = {real(tmp[i-reduced][j]),
                                          imag(tmp[i-reduced][j]) + len};
                        num_grid[X][len]=num_grid[X][Y];
                        num_grid[X][Y]=0;
                        num_grid[X][len] = checkHit(tmp, grid, num_grid,
                                                    X, len, r_p, reduced);
                    }
                    else {
                        std::cout << "Already something there" << std::endl;
                        //Hit!
                    }
                }
                else {
                    if ((real(grid[X][Y-1]) == 0) && (imag(grid[X][Y-1]) == 0)){
                        grid[X][Y-1] = {real(clusters[i][j]),
                                        imag(clusters[i][j]) - L_step};
                        grid[X][Y] = {0,0};
                        clusters[i][j] = {real(clusters[i][j]),
                                          imag(clusters[i][j]) - L_step};
                        tmp[i][j] = {real(tmp[i-reduced][j]),
                                          imag(tmp[i-reduced][j]) - L_step};
                        num_grid[X][Y-1]=num_grid[X][Y];
                        num_grid[X][Y]=0;
                        num_grid[X][Y-1] = checkHit(tmp, grid, num_grid,
                                                    X, Y-1, r_p, reduced);
                    }
                    else {
                        std::cout << "Already something there" << std::endl;
                        //Hit!
                    }
                }
            }
            else if (step_dir == 3){
                X = real(clusters[i][j]);
                Y = imag(clusters[i][j]);
                if (X == 0) {
                    if ((real(grid[len][Y]) == 0) && (imag(grid[len][Y]) == 0)) {
                        grid[len][Y] = {real(clusters[i][j]) + len,
                                        imag(clusters[i][j])};
                        grid[X][Y] = {0,0};
                        clusters[i][j] = {real(clusters[i][j]) + len,
                                          imag(clusters[i][j])};
                        tmp[i][j] = {real(tmp[i-reduced][j]) + len,
                                          imag(tmp[i-reduced][j])};
                        num_grid[len][Y] = num_grid[X][Y];
                        num_grid[X][Y]=0;
                        num_grid[len][Y] = checkHit(tmp, grid, num_grid,
                                                    len, Y, r_p, reduced);
                    }
                    else {
                        std::cout << "Already something there" << std::endl;
                        //Hit!
                    }
                }
                else {
                    if ((real(grid[X-1][Y]) == 0) && (imag(grid[X-1][Y]) == 0)){
                        grid[X-1][Y] = {real(clusters[i][j]) - L_step,
                                        imag(clusters[i][j])};
                        grid[X][Y] = {0,0};
                        clusters[i][j] = {real(clusters[i][j]) - L_step,
                                          imag(clusters[i][j])};
                        tmp[i][j] = {real(tmp[i-reduced][j]) - L_step,
                                          imag(tmp[i-reduced][j])};
                        num_grid[X-1][Y]=num_grid[X][Y];
                        num_grid[X][Y]=0;
                        num_grid[X-1][Y] = checkHit(tmp, grid, num_grid,
                                                    X-1, Y, r_p, reduced);
                    }
                    else {
                        std::cout << "Already something there" << std::endl;
                        //Hit!
                    }
                }
            }
        }
    }
    clusters = tmp;
}

int joinClustersGrid(std::vector<std::vector<std::complex<double>>> &clusters,
                     std::vector<std::vector<int>> &num_grid, int A, int B,
                     int &reduced){
    std::cout << "A,B = " << A << "," << B << std::endl;
    int ret_var = 9;
    int C = std::min(A,B);
    int D = std::max(A,B);
    for (unsigned int i = 0; i < clusters[D].size(); ++i){
        clusters[C].push_back(clusters[D][i]);
    }
    clusters.erase(clusters.begin() + D);
    reduced++; 
    ret_var = C;
    return ret_var;
}

int checkHit(std::vector<std::vector<std::complex<double>>> &clusters,
              std::vector<std::vector<std::complex<double>>> grid,
              std::vector<std::vector<int>> &num_grid,
              int row, int col, float r_p, int &reduced){
    double r = grid.size();
    int ret_var;
    int len = grid.size()-1;
    if (row == len){
        if ((num_grid[0][col] != 0) && (num_grid[0][col] < num_grid[row][col])){
            for (int i = 0; i < clusters[num_grid[0][col]].size(); ++i){
                r = std::sqrt(std::pow(real(grid[row][col]) -
                                   real(clusters[num_grid[0][col]][i]),2) +
                              std::pow(imag(grid[row][col]) -
                                   imag(clusters[num_grid[0][col]][i]),2));
//                if (r == 2*r_p){
                    ret_var = joinClustersGrid(clusters, num_grid,
                                               num_grid[row][col],
                                               num_grid[0][col], reduced);
                    std::cout << "We are in 1, num_grid[row][col] = "
                              << num_grid[row][col];
                    return ret_var;
//              }
            } 
        }
        else if ((num_grid[row-1][col] != 0) &&
                 (num_grid[row-1][col] < num_grid[row][col])){
            for (int i = 0; i < clusters[num_grid[row-1][col]].size(); ++i){
                r = std::sqrt(std::pow(real(grid[row][col]) -
                                   real(clusters[num_grid[row-1][col]][i]),2) +
                              std::pow(imag(grid[row][col]) -
                                   imag(clusters[num_grid[row-1][col]][i]),2));
//                if (r == 2*r_p){
                    ret_var = joinClustersGrid(clusters, num_grid,
                                               num_grid[row][col],
                                               num_grid[row-1][col], reduced);
                    std::cout << "We are in 2, num_grid[row][col] = "
                              << num_grid[row][col];
                    return ret_var;
//               }
            }
        }
        else if ((num_grid[row][col+1] != 0) &&
                 (num_grid[row][col+1] < num_grid[row][col])){
            for (int i = 0; i < clusters[num_grid[row][col+1]].size(); ++i){
                r = std::sqrt(std::pow(real(grid[row][col]) -
                                   real(clusters[num_grid[row][col+1]][i]),2) +
                              std::pow(imag(grid[row][col]) -
                                   imag(clusters[num_grid[row][col+1]][i]),2));
//                if (r == 2*r_p){
                    ret_var = joinClustersGrid(clusters, num_grid,
                                               num_grid[row][col],
                                               num_grid[row][col+1], reduced);
                    std::cout << "We are in 3, num_grid[row][col] = "
                              << num_grid[row][col] << std::endl;
                    return ret_var;
//                }
            }
        }
        else if ((num_grid[row][col-1] != 0) &&
                 (num_grid[row][col-1] < num_grid[row][col])){
            for (int i = 0; i < clusters[num_grid[row][col-1]].size(); ++i){
                r = std::sqrt(std::pow(real(grid[row][col]) -
                                   real(clusters[num_grid[row][col-1]][i]),2) +
                              std::pow(imag(grid[row][col]) -
                                   imag(clusters[num_grid[row][col-1]][i]),2));
//                if (r == 2*r_p){
                    ret_var = joinClustersGrid(clusters, num_grid,
                                               num_grid[row][col],
                                               num_grid[row][col-1], reduced);
                    std::cout << "We are in 4, num_grid[row][col] = "
                              << num_grid[row][col] << std::endl;
                    return ret_var;
//                }
            }
        }
    }
    else if (row == 0){
        if ((num_grid[row+1][col] != 0) &&
            (num_grid[row+1][col] < num_grid[row][col])){
            for (int i = 0; i < clusters[num_grid[row+1][col]].size(); ++i){
                r = std::sqrt(std::pow(real(grid[row][col]) -
                                   real(clusters[num_grid[row+1][col]][i]),2) +
                              std::pow(imag(grid[row][col]) -
                                   imag(clusters[num_grid[row+1][col]][i]),2));
//                if (r == 2*r_p){
                    ret_var = joinClustersGrid(clusters, num_grid,
                                               num_grid[row][col],
                                                num_grid[row+1][col], reduced);
                    std::cout << "We are in 5, num_grid[row][col] = "
                              << num_grid[row][col] << std::endl;
                    return ret_var;
//                }
            } 
        }
        else if ((num_grid[len][col] != 0) &&
                 (num_grid[len][col] < num_grid[row][col])){
            for (int i = 0; i < clusters[num_grid[len][col]].size(); ++i){
                r = std::sqrt(std::pow(real(grid[row][col]) -
                                   real(clusters[num_grid[len][col]][i]),2) +
                              std::pow(imag(grid[row][col]) -
                                   imag(clusters[num_grid[len][col]][i]),2));
//                if (r == 2*r_p){
                    ret_var = joinClustersGrid(clusters, num_grid,
                                               num_grid[row][col],
                                               num_grid[len][col], reduced);
                    std::cout << "We are in 6, num_grid[row][col] = "
                              << num_grid[row][col] << std::endl;
                    return ret_var;
//               }
            }
        }
        else if ((num_grid[row][col+1] != 0) &&
                (num_grid[row][col+1] < num_grid[row][col])){
            for (int i = 0; i < clusters[num_grid[row][col+1]].size(); ++i){
                r = std::sqrt(std::pow(real(grid[row][col]) -
                                   real(clusters[num_grid[row][col+1]][i]),2) +
                              std::pow(imag(grid[row][col]) -
                                   imag(clusters[num_grid[row][col+1]][i]),2));
//                if (r == 2*r_p){
                    ret_var = joinClustersGrid(clusters, num_grid,
                                               num_grid[row][col],
                                               num_grid[row][col+1], reduced);
                    std::cout << "We are in 7, num_grid[row][col] = "
                              << num_grid[row][col] << std::endl;
                    return ret_var;
//                }
            }
        }
        else if ((num_grid[row][col-1] != 0) &&
                 (num_grid[row][col-1] < num_grid[row][col])){
            for (int i = 0; i < clusters[num_grid[row][col-1]].size(); ++i){
                r = std::sqrt(std::pow(real(grid[row][col]) -
                                   real(clusters[num_grid[row][col-1]][i]),2) +
                              std::pow(imag(grid[row][col]) -
                                   imag(clusters[num_grid[row][col-1]][i]),2));
//                if (r == 2*r_p){
                    ret_var = joinClustersGrid(clusters, num_grid,
                                               num_grid[row][col],
                                               num_grid[row][col-1], reduced);
                    std::cout << "We are in 8, num_grid[row][col] = "
                              << num_grid[row][col] << std::endl;
                    return ret_var;
//                }
            }
        }
    }
    else if (col == len){
        if ((num_grid[row+1][col] != 0) &&
            (num_grid[row+1][col] < num_grid[row][col])){
            for (int i = 0; i < clusters[num_grid[row+1][col]].size(); ++i){
                r = std::sqrt(std::pow(real(grid[row][col]) -
                                   real(clusters[num_grid[row+1][col]][i]),2) +
                              std::pow(imag(grid[row][col]) -
                                   imag(clusters[num_grid[row+1][col]][i]),2));
//                if (r == 2*r_p){
                    ret_var = joinClustersGrid(clusters, num_grid,
                                               num_grid[row][col],
                                               num_grid[row+1][col], reduced);
                    std::cout << "We are in 9, num_grid[row][col] = "
                              << num_grid[row][col] << std::endl;
                    return ret_var;
//                }
            } 
        }
        else if ((num_grid[row-1][col] != 0) &&
                 (num_grid[row-1][col] < num_grid[row][col])){
            for (int i = 0; i < clusters[num_grid[row-1][col]].size(); ++i){
                r = std::sqrt(std::pow(real(grid[row][col]) -
                                   real(clusters[num_grid[row-1][col]][i]),2) +
                              std::pow(imag(grid[row][col]) -
                                   imag(clusters[num_grid[row-1][col]][i]),2));
//                if (r == 2*r_p){
                    ret_var = joinClustersGrid(clusters, num_grid,
                                               num_grid[row][col],
                                               num_grid[row-1][col], reduced);
                    std::cout << "We are in 10, num_grid[row][col] = "
                              << num_grid[row][col] << std::endl;
                    return ret_var;
//               }
            }
        }
        else if ((num_grid[row][0] != 0) &&
                 (num_grid[row][0] < num_grid[row][col])){
            for (int i = 0; i < clusters[num_grid[row][0]].size(); ++i){
                r = std::sqrt(std::pow(real(grid[row][col]) -
                                   real(clusters[num_grid[row][0]][i]),2) +
                              std::pow(imag(grid[row][col]) -
                                   imag(clusters[num_grid[row][0]][i]),2));
//                if (r == 2*r_p){
                    ret_var = joinClustersGrid(clusters, num_grid,
                                               num_grid[row][col],
                                               num_grid[row][0], reduced);
                    std::cout << "We are in 11, num_grid[row][col] = "
                              << num_grid[row][col] << std::endl;
                    return ret_var;
//                }
            }
        }
        else if ((num_grid[row][col-1] != 0) &&
                (num_grid[row][col-1] < num_grid[row][col])){
            for (int i = 0; i < clusters[num_grid[row][col-1]].size(); ++i){
                r = std::sqrt(std::pow(real(grid[row][col]) -
                                   real(clusters[num_grid[row][col-1]][i]),2) +
                              std::pow(imag(grid[row][col]) -
                                   imag(clusters[num_grid[row][col-1]][i]),2));
//                if (r == 2*r_p){
                    ret_var = joinClustersGrid(clusters, num_grid,
                                               num_grid[row][col],
                                               num_grid[row][col-1], reduced);
                    std::cout << "We are in 12, num_grid[row][col] = "
                              << num_grid[row][col] << std::endl;
                    return ret_var;
//                }
            }
        }
    }
    else if (col == 0){
        if ((num_grid[row+1][col] != 0) &&
            (num_grid[row+1][col] < num_grid[row][col])){
            for (int i = 0; i < clusters[num_grid[row+1][col]].size(); ++i){
                r = std::sqrt(std::pow(real(grid[row][col]) -
                                   real(clusters[num_grid[row+1][col]][i]),2) +
                              std::pow(imag(grid[row][col]) -
                                   imag(clusters[num_grid[row+1][col]][i]),2));
//                if (r == 2*r_p){
                    ret_var = joinClustersGrid(clusters, num_grid,
                                               num_grid[row][col],
                                               num_grid[row+1][col], reduced);
                    std::cout << "We are in 13, num_grid[row][col] = "
                              << num_grid[row][col] << std::endl;
                    return ret_var;
//                }
            } 
        }
        else if ((num_grid[row-1][col] != 0) &&
                 (num_grid[row-1][col] < num_grid[row][col])){
            for (int i = 0; i < clusters[num_grid[row-1][col]].size(); ++i){
                r = std::sqrt(std::pow(real(grid[row][col]) -
                                   real(clusters[num_grid[row-1][col]][i]),2) +
                              std::pow(imag(grid[row][col]) -
                                   imag(clusters[num_grid[row-1][col]][i]),2));
//                if (r == 2*r_p){
                    ret_var = joinClustersGrid(clusters, num_grid,
                                               num_grid[row][col],
                                               num_grid[row-1][col], reduced);
                    std::cout << "We are in 14, num_grid[row][col] = "
                              << num_grid[row][col] << std::endl;
                    return ret_var;
//               }
            }
        }
        else if ((num_grid[row][col+1] != 0) &&
                 (num_grid[row][col+1] < num_grid[row][col])){
            for (int i = 0; i < clusters[num_grid[row][col+1]].size(); ++i){
                r = std::sqrt(std::pow(real(grid[row][col]) -
                                   real(clusters[num_grid[row][col+1]][i]),2) +
                              std::pow(imag(grid[row][col]) -
                                   imag(clusters[num_grid[row][col+1]][i]),2));
//                if (r == 2*r_p){
                    ret_var = joinClustersGrid(clusters, num_grid,
                                               num_grid[row][col],
                                               num_grid[row][col+1], reduced);
                    std::cout << "We are in 15, num_grid[row][col] = "
                              << num_grid[row][col] << std::endl;
                    return ret_var;
//                }
            }
        }
        else if ((num_grid[row][len] != 0) &&
                (num_grid[row][len] < num_grid[row][col])){
            for (int i = 0; i < clusters[num_grid[row][len]].size(); ++i){
                r = std::sqrt(std::pow(real(grid[row][col]) -
                                   real(clusters[num_grid[row][len]][i]),2) +
                              std::pow(imag(grid[row][col]) -
                                   imag(clusters[num_grid[row][len]][i]),2));
//                if (r == 2*r_p){
                    ret_var = joinClustersGrid(clusters, num_grid,
                                               num_grid[row][col],
                                               num_grid[row][len], reduced);
                    std::cout << "We are in 16, num_grid[row][col] = "
                              << num_grid[row][col] << std::endl;
                    return ret_var;
//                }
            }
        }
    }
    else {
        if ((num_grid[row+1][col] != 0) &&
            (num_grid[row+1][col]<num_grid[row][col])){
            for (int i = 0; i < clusters[num_grid[row+1][col]].size(); ++i){
                r = std::sqrt(std::pow(real(grid[row][col]) -
                                   real(clusters[num_grid[row+1][col]][i]),2) +
                              std::pow(imag(grid[row][col]) -
                                   imag(clusters[num_grid[row+1][col]][i]),2));
//                if (r == 2*r_p){
                    ret_var = joinClustersGrid(clusters, num_grid,
                                               num_grid[row][col],
                                               num_grid[row+1][col], reduced);
                    std::cout << "We are in 17, num_grid[row][col] = "
                              << num_grid[row][col] << std::endl;
                    return ret_var;
//                }
            } 
        }
        else if ((num_grid[row-1][col] != 0) &&
                 (num_grid[row-1][col]) < (num_grid[row][col])){
            for (int i = 0; i < clusters[num_grid[row-1][col]].size(); ++i){
                r = std::sqrt(std::pow(real(grid[row][col]) -
                                   real(clusters[num_grid[row-1][col]][i]),2) +
                              std::pow(imag(grid[row][col]) -
                                   imag(clusters[num_grid[row-1][col]][i]),2));
//                if (r == 2*r_p){
                    ret_var = joinClustersGrid(clusters, num_grid,
                                               num_grid[row][col],
                                               num_grid[row-1][col], reduced);
                    std::cout << "We are in 18, num_grid[row][col] = "
                              << num_grid[row][col] << std::endl;
                    return ret_var;
//                }
            }
        }
        else if ((num_grid[row][col+1] != 0) && 
                 (num_grid[row][col+1] < num_grid[row][col])){
            for (int i = 0; i < clusters[num_grid[row][col+1]].size(); ++i){
                r = std::sqrt(std::pow(real(grid[row][col]) -
                                   real(clusters[num_grid[row][col+1]][i]),2) +
                              std::pow(imag(grid[row][col]) -
                                   imag(clusters[num_grid[row][col+1]][i]),2));
//                if (r == 2*r_p){
                    ret_var = joinClustersGrid(clusters, num_grid,
                                               num_grid[row][col],
                                               num_grid[row][col+1], reduced);
                    std::cout << "We are in 19, num_grid[row][col] = "
                              << num_grid[row][col] << std::endl;
                    return ret_var;
//                }
            }
        }
        else if ((num_grid[row][col-1] != 0) &&
                (num_grid[row][col-1]) < num_grid[row][col]){
            for (int i = 0; i < clusters[num_grid[row][col-1]].size(); ++i){
                r = std::sqrt(std::pow((real(grid[row][col]) -
                                   real(clusters[num_grid[row][col-1]][i])),2) +
                              std::pow((imag(grid[row][col]) -
                                   imag(clusters[num_grid[row][col-1]][i])),2));
//                if (r == 2*r_p){
                    ret_var = joinClustersGrid(clusters, num_grid,
                                               num_grid[row][col],
                                               num_grid[row][col-1], reduced);
                    std::cout << "We are in 20, num_grid[row][col] = "
                              << num_grid[row][col] << std::endl;
                    return ret_var;
//                }
            }
        }
        else if (num_grid[row-1][col-1] != 0){
        //do nothing atm.
        return num_grid[row][col];
        }
        else if (num_grid[row+1][col-1] != 0){
        //do nothing atm.
        return num_grid[row][col];
        }
        else if (num_grid[row+1][col+1] != 0){
        //do nothing atm.
        }
        else if (num_grid[row-1][col+1] != 0){
        //do nothing atm.
        }
    }
    return num_grid[row][col];
}


bool checkHitGrid(std::vector<std::vector<std::complex<double>>> &grid, int row,
              int col){
    if ((int(real(grid[row+2][col])) != 0) || (int(imag(grid[row+2][col])))){
        return true;
    }
    else if ((int(real(grid[row-2][col])) != 0) ||
             (int(imag(grid[row-2][col])))){
        return true;
    }
    else if ((int(real(grid[row][col+2])) != 0) ||
             (int(imag(grid[row][col+2])))){
        return true;
    }
    else if ((int(real(grid[row][col-2])) != 0) ||
             (int(imag(grid[row][col-2])))){
        return true;
    }
    else{
        return false;
    }
}
