#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <algorithm>
#include <math.h>
#include <string>

#include "header.h"

void initGrid(int system_length, std::mt19937::result_type seed,
              int nbr_particles, float r_p,
              std::vector<std::vector<std::complex<double>>> &clusters,
              std::vector<std::vector<int>> &num_grid_C,
              std::vector<std::vector<int>> &num_grid_N){
    int counter = 1;
    double tmp1, tmp2, tmp3, tmp4;
    bool occupied;
    float r;
    double x_c, y_c;
    for (unsigned int i = 0; i < num_grid_C.size(); ++i){
        for (unsigned int j = 0; j < num_grid_C.size(); ++j){
            num_grid_C[i][j] = 0;
            num_grid_N[i][j] = 0;
        }
    }
    auto coord = std::bind(std::uniform_real_distribution<double>
                           (0, double(num_grid_C.size())),
                           std::mt19937(seed));
    while (counter < nbr_particles){
        tmp1 = coord();
        tmp2 = coord();
        occupied = false;
        for (int i = 1; i < counter; ++i){
            tmp3 = tmp1;
            tmp4 = tmp2;
            x_c = real(clusters[i][0]);
            y_c = imag(clusters[i][0]);
            findShortestDist(num_grid_C.size(), x_c, y_c, tmp3, tmp4, r_p);
            r = std::sqrt(pow(x_c - tmp3, 2) + pow(y_c - tmp4, 2));
            if (r < 2.0*r_p){
                occupied = true;
            }
        }
        if (occupied == false){
            clusters[counter][0] = {tmp1, tmp2};
            num_grid_C[floor(tmp1)][floor(tmp2)] = counter;
            ++counter;
        }
    }
}

void walkOnGrid(std::vector<std::vector<std::complex<double>>> &clusters,
                std::vector<std::vector<int>> &num_grid_C,
                std::vector<std::vector<int>> &num_grid_N,
                std::mt19937::result_type seed, float L_step, double pi,
                float r_p, int &counter){
    int row, col;
    int len = num_grid_C.size();
    double X, Y, nextX, nextY, step_dir, L_min, L_org, x_c, y_c, x_1, y_1;
    int has_hit_number = 0;
    bool has_hit = false;
    bool can_move = false;
    //auto rand_dir_int = std::bind(std::uniform_real_distribution<float>(0,4),
	//                std::mt19937(seed));
    auto rand_dir = std::bind(std::uniform_real_distribution<double>(0,2*pi),
				   std::mt19937(seed));
    L_org = L_step;
    for (int i = 1; i < clusters.size(); ++i){
        L_step = L_org;
        can_move = true;
        L_min = L_step;
        has_hit_number = 0;
        step_dir = rand_dir();
        std::vector<std::complex<double>> moved_from;
        std::vector<std::complex<double>> to_check;
        std::vector<std::complex<int>> moved_from_labels;
        std::vector<std::complex<int>> moved_to_labels;
        for (int j = 0; j < clusters[i].size(); ++j){
            X = std::real(clusters[i][j]);
            Y = std::imag(clusters[i][j]);
            moved_from.push_back({X,Y});
            moved_from_labels.push_back({floor(X),floor(Y)});
            makeStep(X, Y, step_dir, len, L_step,
                     nextX, nextY);
            row = floor(nextX);
            col = floor(nextY);
            checkDestination(clusters, num_grid_C, len, num_grid_N,
                             row, col, to_check, r_p);
            for (int k = 0; k < to_check.size(); ++k){
                x_c = real(to_check[k]);
                y_c = imag(to_check[k]);
                x_1 = X;
                y_1 = Y;
                findShortestDist(len, x_c, y_c, x_1, y_1, r_p);
                if (num_grid_C[X][Y] !=
                    num_grid_C[floor(real(to_check[k]))]
                              [floor(imag(to_check[k]))]){
                    L_step = LHit(L_step, step_dir, x_c, y_c , X, Y, r_p);
                }
                if (L_step < L_min){
                    L_min = L_step;
                    has_hit_number = num_grid_C[floor(real(to_check[k]))]
                                               [floor(imag(to_check[k]))];
                }
            }
        }
        std::vector<std::vector<int>> ngn_tmp = num_grid_N;        
        for (int j = 0; j < moved_from.size(); ++j){
            makeStep(real(moved_from[j]), imag(moved_from[j]), step_dir, len,
                     L_min, nextX, nextY);
            moved_to_labels.push_back({floor(nextX),floor(nextY)});
            updateClusters(clusters, i, j, nextX, nextY);
            updateNum_grid_C(num_grid_C, moved_to_labels, moved_from[j],
                            {nextX, nextY});
            updateNum_grid_N(num_grid_N, ngn_tmp, real(moved_to_labels[j]),
                             imag(moved_to_labels[j]), 
                             real(moved_from_labels[j]),
                             imag(moved_from_labels[j]));
        }
        num_grid_N = ngn_tmp;
        if (has_hit_number != 0){
            joinClusters(clusters, num_grid_C, i, has_hit_number, num_grid_N);
            updateNum_grid(num_grid_C, std::max(i, has_hit_number), clusters);
        }
    }
}

void findShortestDist(int len, double &x_c, double &y_c, double &X, double &Y,
                      float r_p){
    if (std::abs(X - x_c) > (len-1 - ceil(2*r_p))){
        if (x_c > X){
            X = X + len;
        }
        else {
            x_c = x_c + len;
        }
    }
    if (std::abs(Y - y_c) > (len-1 - ceil(2*r_p))){
        if (y_c > Y){
            Y = Y + len;
        }
        else {
            y_c = y_c + len;
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

void checkDestination(std::vector<std::vector<std::complex<double>>> clusters,
                      std::vector<std::vector<int>> num_grid_C, int len,
                      std::vector<std::vector<int>> num_grid_N, int row,
                      int col, std::vector<std::complex<double>> &to_check,
                      float r_p){
    int A, B;
    if ((row - ceil(2*r_p) <= 0) && (col - ceil(2*r_p) <= 0)){
//        std::cout << "6" << std::endl;
        A = row + ceil(2*r_p)+1;
        B = col + ceil(2*r_p)+1;
        for (int i = 0; i < A; ++i){
            for (int j = 0; j < B; ++j){
                if (num_grid_C[i][j] != 0){
                    to_check.push_back(clusters[num_grid_C[i][j]]
                                               [num_grid_N[i][j]]);
                }
            }
        }
        for (int i = row + (len-1) - ceil(2*r_p); i < len; ++i){
            for (int j = col + (len-1) - ceil(2*r_p); j < len; ++j){
                if (num_grid_C[i][j] != 0){
                    to_check.push_back(clusters[num_grid_C[i][j]]
                                               [num_grid_N[i][j]]);
                }                
            }
        }
        A = row + ceil(2*r_p)+1;
        for (int i = 0; i < A; ++i){
            for (int j = col + (len-1) - ceil(2*r_p); j < len; ++j){
                if (num_grid_C[i][j] != 0){
                    to_check.push_back(clusters[num_grid_C[i][j]]
                                               [num_grid_N[i][j]]);
                } 
            }
        }
        B = col + ceil(2*r_p)+1;
        for (int i = row + (len-1) - ceil(2*r_p); i < len; ++i){
            for (int j = 0; j < B; ++j){
                if (num_grid_C[i][j] != 0){
                    to_check.push_back(clusters[num_grid_C[i][j]]
                                               [num_grid_N[i][j]]);
                } 
            }
        }
    }
    else if ((row + ceil(2*r_p) > len-1) && (col + ceil(2*r_p) > len-1)){
//        std::cout << "4" << std::endl;
        for (int i = row - ceil(2*r_p); i < len; ++i){
            for (int j = col - ceil(2*r_p); j < len; ++j){
                if (num_grid_C[i][j] != 0){
                    to_check.push_back(clusters[num_grid_C[i][j]]
                                               [num_grid_N[i][j]]);
                }
            }
        }
        A = ceil(2*r_p) + row - (len-1);
        B = ceil(2*r_p) + col - (len-1);
        for (int i = 0; i < A; ++i){
            for (int j = 0; j < B; ++j){
                if (num_grid_C[i][j] != 0){
                    to_check.push_back(clusters[num_grid_C[i][j]]
                                               [num_grid_N[i][j]]);
                }
            }
        }
        B = ceil(2*r_p) + col - (len-1);
        for (int i = row - ceil(2*r_p); i < len; ++i){
           for (int j = 0; j < B; ++j){
                if (num_grid_C[i][j] != 0){
                    to_check.push_back(clusters[num_grid_C[i][j]]
                                               [num_grid_N[i][j]]);
                }
            }
        }
        A = ceil(2*r_p) + row - (len-1);
        for (int i = 0; i < A; ++i){
            for (int j = col - ceil(2*r_p); j < len; ++j){
                if (num_grid_C[i][j] != 0){
                    to_check.push_back(clusters[num_grid_C[i][j]]
                                               [num_grid_N[i][j]]);
                }
            }
        }
    }
    else if ((row -ceil(2*r_p) <= 0) && (col + ceil(2*r_p) > len-1)){
//        std::cout << "7" << std::endl;
        A = row + ceil(2*r_p)+1;
        for (int i = 0; i < A; ++i){
            for (int j = col - ceil(2*r_p); j < len; ++j){
                if (num_grid_C[i][j] != 0){
                    to_check.push_back(clusters[num_grid_C[i][j]]
                                               [num_grid_N[i][j]]);
                }
            }
        }
        B = col - (len-1) + ceil(2*r_p);
        for (int i = len-1; i > len - ceil(2*r_p); --i){
            for (int j = 0; j < B; ++j){
                if (num_grid_C[i][j] != 0){
                    to_check.push_back(clusters[num_grid_C[i][j]]
                                               [num_grid_N[i][j]]);
                }
            }
        }
        for (int i = 0; i < row + ceil(2*r_p)+1; ++i){
            for (int j = 0; j < B; ++j){
                if (num_grid_C[i][j] != 0){
                    to_check.push_back(clusters[num_grid_C[i][j]]
                                               [num_grid_N[i][j]]);
                }
            }
        }
        for (int i = len-1; i > len - ceil(2*r_p); --i){
            for (int j = col - ceil(2*r_p); j < len; ++j) {
                if (num_grid_C[i][j] != 0){
                    to_check.push_back(clusters[num_grid_C[i][j]]
                                               [num_grid_N[i][j]]);
                }
            }
        }
    }
    else if ((row + ceil(2*r_p) > len-1) && (col - ceil(2*r_p) <= 0)){
//        std::cout << "5" << std::endl;
        B = col + ceil(2*r_p)+1;
        for (int i = row - ceil(2*r_p); i < len; ++i){
            for (int j = 0; j < B; ++j){
                if (num_grid_C[i][j] != 0){
                    to_check.push_back(clusters[num_grid_C[i][j]]
                                               [num_grid_N[i][j]]);
                }
            }
        }
        for (int i = row - ceil(2*r_p); i < len; ++i){
            for (int j = col + (len-1) - ceil(2*r_p); j < len; ++j){
                if (num_grid_C[i][j] != 0){
                    to_check.push_back(clusters[num_grid_C[i][j]]
                                               [num_grid_N[i][j]]);
                }
            }
        }
        A = ceil(2*r_p) + row - (len-1);
        for (int i = 0; i < A; ++i){
            for (int j = 0; j < col + ceil(2*r_p)+1; ++j){
                if (num_grid_C[i][j] != 0){
                    to_check.push_back(clusters[num_grid_C[i][j]]
                                               [num_grid_N[i][j]]);
                }
            }
        }
        A = ceil(2*r_p) + row - (len-1);
        for (int i = 0; i < A; ++i){
            for (int j = col + (len-1) - ceil(2*r_p); j < len; ++j){
                if (num_grid_C[i][j] != 0){
                    to_check.push_back(clusters[num_grid_C[i][j]]
                                               [num_grid_N[i][j]]);
                }
            }
        }
    }
    else if ((col + ceil(2*r_p) > len-1) && ((row + ceil(2*r_p)) <= len-1 &&
                                             (row - ceil(2*r_p)) >= 0)){
//        std::cout << "0" << std::endl;
        for (int i = row - ceil(2*r_p); i < row + ceil(2*r_p)+1; ++i){
            for (int j = col - ceil(2*r_p); j < len; ++j){
                if (num_grid_C[i][j] != 0){
                    to_check.push_back(clusters[num_grid_C[i][j]]
                                               [num_grid_N[i][j]]);
                }
            }
        }
        B = ceil(2*r_p) + col - (len-1);
        for (int i = row - ceil(2*r_p); i < row + ceil(2*r_p)+1; ++i){
            for (int j = 0; j < B; ++j){
                if (num_grid_C[i][j] != 0){
                    to_check.push_back(clusters[num_grid_C[i][j]]
                                               [num_grid_N[i][j]]);
                }
            }
        }
    }
    else if (((row + ceil(2*r_p)) > len-1) && ((col + ceil(2*r_p) <= len-1) &&
                                               (col + ceil(2*r_p) >= 0))){
//        std::cout << "1" << std::endl;
        for (int i = row - ceil(2*r_p); i < len; ++i){
            for (int j = col - ceil(2*r_p); j < col + ceil(2*r_p)+1; ++j){
                if (num_grid_C[i][j] != 0){
                    to_check.push_back(clusters[num_grid_C[i][j]]
                                               [num_grid_N[i][j]]);
                }
            }
        }
        A = ceil(2*r_p) + row - (len-1);
        for (int i = 0; i < A; ++i){
            for (int j = col - ceil(2*r_p); j < col + ceil(2*r_p)+1; ++j){
                if (num_grid_C[i][j] != 0){
                    to_check.push_back(clusters[num_grid_C[i][j]]
                                               [num_grid_N[i][j]]);
                }
            }
        }
    }
    else if ((col - ceil(2*r_p) <= 0) && ((row + ceil(2*r_p)) <= len-1 &&
                                         (row - ceil(2*r_p)) >= 0)){
//        std::cout << "2" << std::endl;
        for (int i = row - ceil(2*r_p); i < row + ceil(2*r_p)+1; ++i){
            for (int j = 0; j < col + ceil(2*r_p)+1; ++j){
                if (num_grid_C[i][j] != 0){
                    to_check.push_back(clusters[num_grid_C[i][j]]
                                               [num_grid_N[i][j]]);
                }
            }
        }
        for (int i = row - ceil(2*r_p); i < row + ceil(2*r_p)+1; ++i){
            for (int j = col + (len-1) - ceil(2*r_p); j < len; ++j){
                if (num_grid_C[i][j] != 0){
                    to_check.push_back(clusters[num_grid_C[i][j]]
                                               [num_grid_N[i][j]]);
                }
            }
        }
    }
    else if ((row - ceil(2*r_p) <= 0) && ((col - ceil(2*r_p) >= 0) && 
                                         (col + ceil(2*r_p) <= len-1))){
//        std::cout << "3" << std::endl;
        A = row + ceil(2*r_p)+1;
        for (int i = 0; i < A; ++i){
            for (int j = col - ceil(2*r_p); j < col + ceil(2*r_p)+1; ++j){
                if (num_grid_C[i][j] != 0){
                    to_check.push_back(clusters[num_grid_C[i][j]]
                                               [num_grid_N[i][j]]);
                }
            }
        }
        for (int i = row + (len-1) - ceil(2*r_p); i < len; ++i){
            for (int j = col - ceil(2*r_p); j < col + ceil(2*r_p)+1; ++j){
                if (num_grid_C[i][j] != 0){
                    to_check.push_back(clusters[num_grid_C[i][j]]
                                               [num_grid_N[i][j]]);
                }
            }
        }
    }
    else {
//        std::cout << "8" << std::endl;
        for (int i = row - ceil(2*r_p); i < row + ceil(2*r_p) + 1; ++i){
            for (int j = col - ceil(2*r_p); j < col + ceil(2*r_p) + 1; ++j){
                if (num_grid_C[i][j] != 0){
                    to_check.push_back(clusters[num_grid_C[i][j]]
                                               [num_grid_N[i][j]]);
                }
            }
        }
    }
}

void updateNum_grid_C(std::vector<std::vector<int>> &num_grid_C,
                      std::vector<std::complex<int>> mt,
                      std::complex<double> from, std::complex<double> to){
    //mt is the vector containing all labels moved to.
    std::complex<int> check;
    check = {floor(real(from)), floor(imag(from))};
    if (std::find(mt.begin(), mt.end(), check) != mt.end()){
        num_grid_C[floor(real(to))][floor(imag(to))] = 
            num_grid_C[floor(real(from))][floor(imag(from))];
    }
    else {
        num_grid_C[floor(real(to))][floor(imag(to))] = 
            num_grid_C[floor(real(from))][floor(imag(from))];
        num_grid_C[floor(real(from))][floor(imag(from))] = 0;
    }
}

void updateNum_grid_N(std::vector<std::vector<int>> num_grid_N,
                      std::vector<std::vector<int>> &ngn_tmp,
                      int rt, int ct, int rf, int cf){
    ngn_tmp[rf][cf] = 0;
    ngn_tmp[rt][ct] = num_grid_N[rf][cf];
}

void updateNum_gridWalk(std::vector<std::vector<int>> &num_grid,
                        std::complex<double> f,
                        std::complex<double> t,
                        std::vector<std::complex<double>> mt){
    //mt is the vector containing all positions moved to.
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

void makeStep(double X, double Y, double step_dir, int len, float L_step,
              double &nextX, double &nextY){
    nextX = X + (L_step * std::cos(step_dir));
    nextY = Y + (L_step * std::sin(step_dir));
    if (nextX < 0){
        nextX = std::abs(nextX + len);
    }
    else if (nextX > len){
        nextX = std::fmod(nextX, len);
    }
    if (nextY < 0){
        nextY = std::abs(nextY + len);
    }
    else if (nextY > len){
        nextY = std::fmod(nextY, len);
    }
}

float LHit(float step_L, double step_dir, double x_c, double y_c, double x_p,
           double y_p, float r_p){
    double d_p = 2*r_p;
    double a = 1.0;
    double b = 2*(std::cos(step_dir)*(x_p-x_c) + std::sin(step_dir)*(y_p-y_c));
    double c = (x_c-x_p)*(x_c-x_p) + (y_c-y_p)*(y_c-y_p) - d_p*d_p;
    double res[2];
    bool sol = true;
    //std::cout << "x_c = " << x_c << std::endl;
    //std::cout << "y_c = " << y_c << std::endl;
    //std::cout << "x_p = " << x_p << std::endl;
    //std::cout << "y_p = " << y_p << std::endl;
    //std::cout << "a = " << a << std::endl;
    //std::cout << "b = " << b << std::endl;
    //std::cout << "c = " << c << std::endl;
    //std::cout << "4*a*c = " << 4*a*c << std::endl;
    //std::cout << std::endl;
    if ((b*b - 4*a*c) < 0){
        sol = false;
        return step_L;
    }
    else {
        res[0] = (-b + sqrt(b*b - 4.0*a*c))/(2.0*a);
        res[1] = (-b - sqrt(b*b - 4.0*a*c))/(2.0*a);
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

void joinClusters(std::vector<std::vector<std::complex<double>>> &clusters,
                  std::vector<std::vector<int>> &num_grid_C, int A, int B,
                  std::vector<std::vector<int>> &num_grid_N){
    int C = std::min(A,B);
    int D = std::max(A,B);
    int row, col;
    int org_size = clusters[C].size(); 
    for (unsigned int i = 0; i < clusters[D].size(); ++i){
        row = floor(real(clusters[D][i]));
        col = floor(imag(clusters[D][i]));
        clusters[C].push_back(clusters[D][i]);
        num_grid_C[row][col] = C;
        num_grid_N[row][col] = org_size + i;
    }
    clusters.erase(clusters.begin() + D);
}
//Joins the clusters labeled A & B together, so that the lowest ID'd cluster
//gains all the other particles, and the highest is erased.

void updateNum_grid(std::vector<std::vector<int>> &num_grid, int start,
                    std::vector<std::vector<std::complex<double>>> clusters){
    for (int i = start; i < clusters.size(); ++i){
        for (int j = 0; j < clusters[i].size(); ++j){
            if (num_grid[floor(real(clusters[i][j]))]
                        [floor(imag(clusters[i][j]))] != 0){
            num_grid[floor(real(clusters[i][j]))][floor(imag(clusters[i][j]))]--;
            }
        }
    }
}
//Updates the num_grid after two clusters join, so all later clusters will
//decrease their ID by 1.

void printDist(std::vector<std::vector<std::complex<double>>> clusters, int len,
               float r_p){
    double d, d1;
    std::vector<std::complex<double>> tmp;
    double X, Y, x_c, y_c;
    for (int i = 1; i < clusters.size(); ++i){
        for (int j = 0; j < clusters[i].size(); ++j){
            tmp.push_back(clusters[i][j]);
        }
    }
    for (int i = 0; i < tmp.size(); ++i){
        for (int j = i; j < tmp.size(); ++j){
            X = real(tmp[i]);
            Y = imag(tmp[i]);
            x_c = real(tmp[j]);
            y_c = imag(tmp[j]);
            d1 = std::sqrt(pow(X - x_c, 2) +
                          pow(Y - y_c, 2));
            findShortestDist(len, x_c, y_c, X, Y, r_p);
            d = std::sqrt(pow(X - x_c, 2) +
                          pow(Y - y_c, 2));
            std::cout << "d" << i+1 << "," << j+1 << " = " << d << "    ";
            //std::cout << tmp[i] << " -> " << tmp[j] << std::endl;
            std::cout << "npbc" << i+1 << "," << j+1 << " = " << d1 << "    "
                      << std::endl;
        }
    }
}







void printClusters(std::vector<std::vector<std::complex<double>>> &clusters){
    for (unsigned int i = 0; i < clusters.size(); ++i){
        for (unsigned int j = 0; j < clusters[i].size(); ++j){
            std::cout << clusters[i][j] << std::endl;
        }
    std::cout << std::endl;
    }
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
                 float r_p, int system_length, int iteration){
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
                 float r_p, int system_length, int iteration){
    std::ofstream out_stream;
    std::string name = "data/00" + std::to_string(iteration) + ".txt";
    out_stream.open(name);
    for (unsigned int i = 1; i < clusters.size(); ++i){
        for (unsigned int j = 0; j < clusters[i].size(); ++j){
            out_stream << std::fmod(real(clusters[i][j]), system_length) << " "
                       << std::fmod(imag(clusters[i][j]), system_length) << " "
                       << r_p << " " << i << std::endl;
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
    out_stream << "set xrange [0:" << system_length << "]" << std::endl;
    out_stream << "set yrange [0:" << system_length << "]" << std::endl;
    out_stream << "set grid xtics lt 0 lw 1 lc rgb \"#000000\" " << std::endl;
    out_stream << "set grid ytics lt 0 lw 1 lc rgb \"#000000\" " << std::endl;
    out_stream << "set style fill transparent solid 1.0 noborder" << std::endl;
    out_stream << "unset key" << std::endl;
    out_stream << "set xtics 0, 1, " << system_length << std::endl;
    out_stream << "set ytics 0, 1, " << system_length << std::endl;
	out_stream << "plot \"data/00" << n << ".txt\" u 1:2:3:4 w circles palette,\
                    \"data/00" << n << ".txt\" " << std::endl;
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
    out_stream << "set terminal png size 1200,1200"
                  " enhanced font \"Helvetica,12\" " << std::endl;
	out_stream << "set output \"fig/" << i << "-" << j << ".png\" "
               << std::endl;
    out_stream << "set xrange [0:" << system_length << "]" << std::endl;
    out_stream << "set yrange [0:" << system_length << "]" << std::endl;
    out_stream << "set style fill transparent solid 0.5 noborder" << std::endl;
//    out_stream << "unset key; unset tics; unset border" << std::endl;
	out_stream << "plot \"data/before.txt\" with circles fc rgb \"red\", \
                   \"data/after.txt\" with circles fc rgb \"navy\" "
               << std::endl;
	system("gnuplot gnuplotter.gnu");
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

void writeSeed(std::mt19937::result_type seed){
    std::ofstream out_stream;
    out_stream.open("seed.txt");
    out_stream << seed << std::endl;
    out_stream.close( );

}

void printSingleDist(double X, double Y, double x_c, double y_c,
                     double nextX, double nextY,
                     std::vector<std::vector<int>> num_grid_C,
                     std::vector<std::vector<int>> num_grid_N,
                     std::vector<std::vector<std::complex<double>>> clusters){
    int A, B;
    int counter = 1;
    int counter2 = 1;
    for (int i = 1; i < clusters.size(); ++i){
        for (int j = 0; j < clusters[i].size(); ++j){
            if ((real(clusters[i][j]) == X) && (imag(clusters[i][j]) == Y)){
                A = counter;
            }
        counter++;
        }
    }
    for (int i = 1; i < clusters.size(); ++i){
        for (int j = 0; j < clusters[i].size(); ++j){
            if ((real(clusters[i][j]) == x_c) && (imag(clusters[i][j]) == y_c)){
                B = counter2;
            }
        counter2++;
        }
    }
    double d1 = std::sqrt(pow(nextX - x_c, 2) +
                          pow(nextY - y_c, 2));
    std::cout << "d" << A << "," << B << " = " << d1
              << std::endl;
}

void printColInfo(double X, double Y, double x_c, double y_c,
                  std::vector<std::vector<std::complex<double>>> clusters){
    int A, B;
    int counter = 1;
    int counter2 = 1;
    for (int i = 1; i < clusters.size(); ++i){
        for (int j = 0; j < clusters[i].size(); ++j){
            if ((real(clusters[i][j]) == X) && (imag(clusters[i][j]) == Y)){
                A = counter;
            }
        counter++;
        }
    }
    for (int i = 1; i < clusters.size(); ++i){
        for (int j = 0; j < clusters[i].size(); ++j){
            if ((real(clusters[i][j]) == x_c) && (imag(clusters[i][j]) == y_c)){
                B = counter2;
            }
        counter2++;
        }
    }
    std::cout << "particle " << A << " collided with " << B << std::endl;
}
