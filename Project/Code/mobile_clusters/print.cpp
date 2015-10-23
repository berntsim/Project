#include <iostream>
#include <fstream>
#include <vector>

#include "header.h"
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

void printGridTest(std::vector<std::vector<std::complex<double>>> const grid){
// Like the title suggests, this function prints whatever grid you would like to
// the screen. 
    for (int i = grid.size(); i >= 0; --i){
        for (int j = 0; j < grid.size(); ++j){
            if (grid[j][i] != double(0)){
                std::cout << grid[i][j];
            }
            else {
                std::cout << 0;
            }
        }
    std::cout << std::endl;
    }
    std::cout << std::endl;
}


void printPartPos(std::vector<std::vector<std::complex<double>>> const grid){
    for (unsigned int i = 0; i < grid.size(); ++i){
        for (unsigned int j = 0; j < grid.size(); ++j){
            if (grid[i][j] != double(0)){
                std::cout << real(grid[i][j]) << ", "
                          << imag(grid[i][j]) << std::endl;
            }
        }
    }
}

void writeConfigGrid(std::vector<std::vector<std::complex<double>>> grid,
                 float r_p){
    std::ofstream out_stream;
    int count = 0;
    out_stream.open("data/config.txt");
    for (unsigned int i = 0; i < grid.size(); ++i){
        for (unsigned int j = 0; j < grid.size(); ++j){
            if ((int(real(grid[i][j])) != 0) || (int(imag(grid[i][j]) != 0))){
                out_stream << real(grid[i][j])  << " "
                           << imag(grid[i][j]) << " " << r_p <<  std::endl;
                count++;
            }
        }
    }
    out_stream.close( );
}

void writeConfig(std::vector<std::vector<std::complex<double>>> clusters,
                 float r_p){
    std::ofstream out_stream;
    out_stream.open("data/config.txt");
    for (unsigned int i = 0; i < clusters.size(); ++i){
        for (unsigned int j = 1; j < clusters[i].size(); ++j){
            out_stream << real(clusters[i][j]) << " "
                       << imag(clusters[i][j]) << " " << r_p << std::endl;
        }
    }
    out_stream.close( );
}

void plotConfig(int iteration, int system_length){
    std::ofstream out_stream;
	out_stream.open("gnuplotter.gnu");
    out_stream << "set terminal png size 1200,1000"
                  " enhanced font \"Helvetica,12\" " << std::endl;
	out_stream << "set output \"fig/" << iteration << ".png\" "
               << std::endl;
    out_stream << "set xrange [0:" << system_length-1 << "]" << std::endl;
    out_stream << "set yrange [0:" << system_length-1 << "]" << std::endl;
    out_stream << "set style fill transparent solid 1.0 noborder" << std::endl;
//    out_stream << "unset key; unset tics; unset border" << std::endl;
	out_stream << "plot \"data/config.txt\" with circles fc rgb \"red\", \
                   \"data/config1.txt\" with circles fc rgb \"navy\" "
               << std::endl;
	system("gnuplot gnuplotter.gnu");
    std::cout << "plotConfig is OK!" << std::endl;
}

void writeConfigGrid1(std::vector<std::vector<std::complex<double>>> grid,
                 float r_p){
    std::ofstream out_stream;
    int count = 0;
    out_stream.open("data/config1.txt");
    for (unsigned int i = 0; i < grid.size(); ++i){
        for (unsigned int j = 0; j < grid.size(); ++j){
            if ((int(real(grid[i][j])) != 0) || (int(imag(grid[i][j]) != 0))){
                out_stream << real(grid[i][j])  << " "
                           << imag(grid[i][j]) << " " << r_p <<  std::endl;
                count++;
            }
        }
    }
    out_stream.close( );
    std::cout << "writeConfig is OK!" << count <<  std::endl;
}

void writeConfig1(std::vector<std::vector<std::complex<double>>> clusters,
                 float r_p, int system_length){
    std::ofstream out_stream;
    out_stream.open("data/config1.txt");
    for (unsigned int i = 0; i < clusters.size(); ++i){
        for (unsigned int j = 1; j < clusters[i].size(); ++j){
            out_stream << std::fmod(real(clusters[i][j]), system_length) << " "
                       << std::fmod(imag(clusters[i][j]), system_length) << " "
                       << r_p << std::endl;
        }
    }
    out_stream.close( );
}

void printParticles(std::vector<std::complex<double>> const part){
    for (unsigned int i = 0; i < part.size(); ++i){
        std::cout << part[i] << std::endl;
    }
}

void speedTest(){
    std::complex<double> a;
    std::complex<double> b;
    int min = 1;
    int max = 3;
    std::complex<double> retVal1 = {0,0};
    int retVal2 = 0;
    retVal2 = fmax(min,max);
}

void printTest(std::vector<std::vector<int>> grid){
    for (int i = grid.size()-1; i >= 0; --i){
        for (int j = 0; j < grid.size(); ++j){
                std::cout << grid[j][i];
            }
        std::cout << std::endl; 
    }
    std::cout << std::endl;
}
