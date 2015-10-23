#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <algorithm>
#include <math.h>

#include "header.h"

void initClusters(std::vector<std::vector<std::complex<double>>> &clusters,
                  std::mt19937::result_type seed, int nbr_particles,
                  int system_length){
    auto rand_pos = std::bind(std::uniform_int_distribution<int>
                              (0, system_length), std::mt19937(seed));
    int counter = 0;
    int tmp1;
    int tmp2;
    while(counter < nbr_particles){
        tmp1 = rand_pos();
        tmp2 = rand_pos();
        //std::cout << "Particle number " << counter << " is placed at (" << tmp1
        //          << "," << tmp2 << ")" << std::endl;
        clusters[counter].push_back({double(tmp1), double(tmp2)});
        ++counter;
    }
}

void findCM(std::vector<std::complex<double>> pos,
            std::complex<double> &ret_var){
    float x_sum = 0;
    float y_sum = 0;
    int j = pos.size();
    for (unsigned int i = 1; i < j; ++i){
        x_sum += std::real(pos[i]);
        y_sum += std::imag(pos[i]);
    }
    ret_var = {x_sum/(j-1),y_sum/(j-1)};
}

void joinClusters(std::vector<std::vector<std::complex<double>>> &clusters,
                 int A, int B){
    std::complex<double> tmp;
    if (clusters[A].size() >= clusters[B].size()){
        for (unsigned int i = 1; i < clusters[B].size(); ++i){
            clusters[A].push_back(clusters[B][i]);
        }
        clusters.erase(clusters.begin() + B);
        findCM(clusters[A],tmp);
        clusters[A][0] = tmp;
    }
    else {
        for(unsigned int i = 1; i < clusters[A].size(); ++i){
            clusters[B].push_back(clusters[A][i]);
        }
        clusters.erase(clusters.begin() + A);
        findCM(clusters[B],tmp);
        clusters[B][0] = tmp;
    }
}

bool checkHit(std::vector<std::vector<std::complex<double>>> &clusters,
              int cluster_number, int particle_number, int &A, int &B,
              int &counter2){
    for (unsigned int i = 0; i < clusters.size(); ++i){
        for (unsigned int j = 1; j < clusters[i].size(); ++j){
            ++counter2;
            if ((real(clusters[cluster_number][particle_number]) == 
                    real(clusters[i][j]) + 2) &&
                (imag(clusters[cluster_number][particle_number]) ==
                 imag(clusters[i][j]))){
                A = fmin(i, cluster_number);
                B = fmax(i, cluster_number);
                clusters[i][0] = {real(clusters[i][j]),imag(clusters[i][j])};
                return true;
            }
            else if ((real(clusters[cluster_number][particle_number]) == 
                      real(clusters[i][j]) - 2) &&
                    (imag(clusters[cluster_number][particle_number]) ==
                     imag(clusters[i][j]))){
                A = fmin(i, cluster_number);
                B = fmax(i, cluster_number);
                clusters[i][0] = {real(clusters[i][j]),imag(clusters[i][j])};
                return true;
            }
            else if ((real(clusters[cluster_number][particle_number]) == 
                    real(clusters[i][j])) &&
                (imag(clusters[cluster_number][particle_number]) ==
                 imag(clusters[i][j]) + 2)){
                A = fmin(i, cluster_number);
                B = fmax(i, cluster_number);
                clusters[i][0] = {real(clusters[i][j]),imag(clusters[i][j])};
                return true;
            }
            else if ((real(clusters[cluster_number][particle_number]) == 
                    real(clusters[i][j])) &&
                (imag(clusters[cluster_number][particle_number]) ==
                 imag(clusters[i][j]) - 2)){
                A = fmin(i, cluster_number);
                B = fmax(i, cluster_number);
                clusters[i][0] = {real(clusters[i][j]),imag(clusters[i][j])};
                return true;
            }
        }
    }
    return false;
}


void walk(std::vector<std::vector<std::complex<double>>> &clusters,
                std::mt19937::result_type seed, float L_step, double pi,
                int system_length, int &counter1, int &counter2){

    std::complex<double> tmp1;
    int step_dir;
    int A;
    int B;
    std::vector<int> As;
    std::vector<int> Bs;
    bool hit = false;
    std::vector<std::vector<std::complex<double>>> tmp2 = clusters;
    auto rand_dir_int = std::bind(std::uniform_real_distribution<float>(0,3),
	                std::mt19937(seed));
    for (unsigned int i = 0; i < clusters.size(); ++i){
        step_dir = rand_dir_int();
        for (unsigned int j = 1; j < clusters[i].size(); ++j){
            ++counter1;
            if (step_dir == 0){
                tmp2[i][j]={fmod(real(clusters[i][j]), (system_length)),
                            fmod(imag(clusters[i][j]), (system_length)) + 1};
            }
            else if (step_dir == 1){
                tmp2[i][j]={fmod(real(clusters[i][j]), system_length) + 1,
                            fmod(imag(clusters[i][j]), system_length)}; 
            }
            else if (step_dir == 2){
                if (fmod(imag(clusters[i][j]), system_length) - 1 == -1){
                    tmp2[i][j]={fmod(real(clusters[i][j]), system_length),
                                double(system_length)};
                }
                else{
                    tmp2[i][j]={fmod(real(clusters[i][j]), system_length),
                             fmod(imag(clusters[i][j]), system_length) - 1};
                }
            }
            else if (step_dir == 3){
                if (fmod(real(clusters[i][j]), system_length) - 1 == -1){
                    tmp2[i][j]={double(system_length),
                                fmod(imag(clusters[i][j]), system_length)};
                }
                else{
                    tmp2[i][j]={fmod(real(clusters[i][j]), system_length) -1,
                                fmod(imag(clusters[i][j]), system_length)}; 
                }
            }
            else {
                std::cout << "\n" << "EXCEPTION!" << "\n" << std::endl;
            }
        }
    }
    std::vector<int>::iterator it;
    for (unsigned int i = 0; i < tmp2.size(); ++i){
        for (unsigned int j = 1; j < tmp2[i].size(); ++j){
            hit = checkHit(tmp2, i, j, A, B, counter2);
            if (hit == true){
                it = std::find(As.begin(), As.end(), A);
                if (it == As.end() && A != B){
                    As.push_back(A);
                    Bs.push_back(B);
                }
            }
        }
    }
    for (int i = 0; i < As.size(); ++i){
        joinClusters(tmp2, As[i]-i, Bs[i]-i);
            tmp2[As[i]-i][0] = {0,0};
        }
    clusters = tmp2;
}
