#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <algorithm>

void initEx_pos(std::vector<float> &x_pos, std::vector<float> &y_pos,
                int sys_len);

void initOnLatticeCluster(std::vector<std::vector<int>> &on_lattice_cluster, 
                                int system_length);

void initOnLatticeDistance(std::vector<std::vector<int>> &on_lattice_distance,
                           int sys_len, int D_max);


void neighbourDist(std::vector<std::vector<int>> &on_lattice_distance,
                   int D_max, int row, int col, int sys_len);

void printGrid(std::vector<std::vector<int>> &grid);

void initVicinity(std::vector<std::vector<int>> &vicinity);

float findRadius(std::vector<float> &x_pos, std::vector<float> &y_pos,
                 float r_p);

void startWalker(float &walker_x_pos, float &walker_y_pos, int D_max,
                 std::vector<float> &x_pos, std::vector<float> &y_pos,
                 float theta, double pi, float r_p, float &r_c);

void updateCluster(std::vector<std::vector<int>> &on_lattice_cluster,
                   std::vector<std::vector<int>> &on_lattice_distance,
                   std::vector<float> &x_pos, std::vector<float> &y_pos,
                   float walker_x_pos, float walker_y_pos, int D_max,
                   int sys_len, float &r_c, float r_p);

float findD_wc(float walker_x_pos, float walker_y_pos, int D_max,
              std::vector<std::vector<int>> &on_lattice_distance);

float LHit(float step_L, float step_dir, float x_c, float y_c, float x_p,
              float y_p, float r_p);

bool isHit();

bool killParticle(float walker_x_pos, float walker_y_pos,
                  std::vector<float> &x_pos, std::vector<float> &y_pos,
                  float r_p, float &r_c);

bool createStep(float &step_L, float step_dir, float &walker_x_pos,
                float &walker_y_pos, double pi, float r_p,
                std::vector<std::vector<int>> &on_lattice_distance,
                std::vector<std::vector<int>> &on_lattice_cluster,
                float L_min, int D_max, std::vector<float> &x_pos,
                std::vector<float> &y_pos, float &r_c);

void runAll(float &step_L, float &step_dir, std::mt19937::result_type seed,
                float &walker_x_pos, float &walker_y_pos, double pi, float r_p,
                std::vector<std::vector<int>> &on_lattice_distance,
                std::vector<std::vector<int>> &on_lattice_cluster,
                float L_min, int D_max, std::vector<float> &x_pos,
                std::vector<float> &y_pos, int sys_len, int nbr_particles,
                unsigned int counter, float &r_c);

void writeGrid(std::vector<float> &x_pos, std::vector<float> &y_pos, float r_p);

void plotGnuplot(int &arr_len, float d_f, float nbr_particles);

void findCM(std::vector<float> x_pos, std::vector<float> y_pos, float &x_CM,
            float &y_CM, unsigned int j);

float findRG(std::vector<float> x_pos, std::vector<float> y_pos, float x_CM,
             float y_CM, unsigned int j);

void fit(std::vector<float> x_pos, std::vector<float> y_pos,
         std::vector<float> &R_g_vec, std::vector<float> &nbr_particles_vec,
         float &x_CM, float &y_CM, int nbr_particles);

void writeFractalDim(std::vector<float> R_g_vec,
                     std::vector<float> nbr_particles_vec);

void plotLogLog(float d_f);

float slope(std::vector<float> &x, std::vector<float> &y);
