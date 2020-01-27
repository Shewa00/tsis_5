#include <iostream>

#include <iomanip>

#include <random>

#include <vector>

#include <cmath>

#include "fstream"

const double a = 0.25, P = 0.95, E = 0.01,

x_min = 0, x_max = M_PI;

const int K = 100, L = 10;

std::random_device Device;

std::mt19937 Mt(Device());

double Function(double x) {
    
    return sin(x) + 0.5;
    
}

std::vector<double> AlphaCount(int r, int M) {
    
    std::vector<double> Alpha(r);
    
    std::uniform_real_distribution<double> Rand(0, 1);
    
    Alpha[M] = Rand(Mt);
    
    for (int m = 1; m < M; m++) {
        
        double Sum = 0;
        
        for (int s = m + 1; s < r - m - 1; s++)
            
            Sum += Alpha[s];
        
        std::uniform_real_distribution<double> Rand2(0, 1 - Sum);
        
        Alpha[r - m - 1] = 0.5 * Rand2(Mt);
        
        Alpha[m] = Alpha[r - m - 1];
        
    }
    
    double Sum = 0;
    
    for (int s = 1; s < r - 1; s++)
        
        Sum += Alpha[s];
    
    Alpha[r - 1] = 0.5 * (1 - Sum);
    
    Alpha[0] = Alpha[r - 1];
    
    return Alpha;
    
}

double AverageCount(std::vector<double> &Y, std::vector<double> &Alpha, int k, int M) {
    
    double Return = 0;
    
    for (int j = k - M; j <= k + M; j++) {
        
        if (j < 1 || j > K + 1)
            
            continue;
        
        Return += Y.at(j - 1) * Alpha.at(j + M - k);
        
    }
    
    return Return;
    
}

std::vector<double> ResultCount(std::vector<double> &Y, std::vector<double> &Alpha, int M) {
    
    std::vector<double> Return;
    
    for (int k = 0; k <= K; k++) {
        
        double Temp = AverageCount(Y, Alpha, k, M);
        
        Return.emplace_back(Temp);
        
    }
    
    return Return;
    
}

double OmegaCount(std::vector<double> &Y) {
    
    double Sum = 0;
    
    for (int k = 1; k < K + 1; k++)
        
        Sum += std::abs(Y[k] - Y[k - 1]);
    
    return Sum;
    
}

double DeltaCount(std::vector<double> &Y, std::vector<double> &_Y) {
    
    double Sum = 0;
    
    for (int k = 0; k < K + 1; k++)
        
        Sum += std::abs(Y[k] - _Y[k]) / K;
    
    return Sum;
    
}

double JCount(double Lmd, double Omega, double Delta) {
    
    return (Lmd * Omega + (1 - Lmd) * Delta * 10000.) / 10000.;
    
}

double Dist(double Omega, double Delta) {
    return std::max(std::abs(Omega), std::abs(Delta));
}

void PrintLine(int r) {
    
    int Amount = 83;
    
    if (r == 5)
        
        Amount +=18;
    
    for (int i = 0; i < Amount; i++)
        
        std::cout << "-";
    
    std::cout << std::endl;
    
}

void Filter(int r) {
    
    std::vector<double> X, Y, _Y;
    
    for (int i = 0; i <= K; i++) {
        
        X.emplace_back(x_min + i * (x_max - x_min) / K);
        
        Y.emplace_back(Function(X[i]));
        
        std::uniform_real_distribution<double> Rand(-a, a);
        
        double Shift = Rand(Mt);
        
        _Y.emplace_back(Y[i] + Shift);
        
    }
    
    std::vector<double> Lmd;
    
    for (int i = 0; i <= L; i++) {
        
        Lmd.emplace_back(double(i) / L);
        
    }
    
    int N = log(1 - P) / log(1 - (E / x_max - x_min)),
    
    M = (r - 1) / 2;
    
    std::vector<std::vector<double>> Alpha;
    
    Alpha.emplace_back(AlphaCount(r, M));
    
    std::vector<double> Y_min = ResultCount(_Y, Alpha[0], M);
    
    size_t Alpha_min = 0,
    
    Lmd_min = 0;
    
    double Omega_min = OmegaCount(Y_min),
    
    Delta_min = DeltaCount(Y_min, _Y),
    
    J_min = JCount(Lmd[Lmd_min], Omega_min, Delta_min);
    
    std::vector<double> omegas, deltas, distances;
    
    std::vector<std::vector<double>> Alpha_optimal;
    
    for (size_t i = 0; i < Lmd.size(); i++) {
        
        for (int j = 0; j < N; j++) {
            
            Alpha.emplace_back(AlphaCount(r, M));
            
            std::vector<double> Y_current = ResultCount(_Y, Alpha[Alpha.size() - 1], M);
            
            double Omega = OmegaCount(Y_current),
            
            Delta = DeltaCount(Y_current, _Y);
            
            if (Dist(JCount(Lmd.at(i), Omega, Delta), 0) < Dist(J_min, 0)) {
                
                Lmd_min = i;
                
                Alpha_min = Alpha.size() - 1;
                
                J_min = JCount(Lmd[i], Omega, Delta);
                
                Omega_min = Omega;
                
                Delta_min = Delta;
                
                Y_min = Y_current;
                
            }
            
        }
        
        Alpha_optimal.emplace_back(Alpha[Alpha_min]);
        
        omegas.emplace_back(Omega_min);
        
        deltas.emplace_back(Delta_min);
        
        distances.emplace_back(Dist(J_min, 0));
        
    }
    
    std::cout << "\n r = " << r << "\n";
    
    PrintLine(r);
    
    std::cout << "|    Lmd     |      W     |      D     |  Distance  |";
    
    if (r == 5)
        
        std::cout << "                       a                       |\n";
    
    else
        
        std::cout << "               a             |\n";
    
    PrintLine(r);
    
    std::ofstream out;
    out.open("../dots.txt");
    
    
    for (int j = 0; j < 11; j++) {
       out << "("<< omegas[j] << ";" << deltas[j] << ") , ";
        std::cout << "| " << std::setprecision(4) << std::setw(10) << Lmd[j];
        
        std::cout << " | " << std::setprecision(6) << std::setw(10) << omegas[j];
        
        std::cout << " | " << std::setprecision(4) << std::setw(10) << deltas[j];
        
        std::cout << " | " << std::setprecision(4) << std::setw(10) << distances[j];
        
        std::cout << " | " ;
        
        for (const auto &A : Alpha_optimal[j])
            
            std::cout << std::setprecision(4) << std::setw(8) << A << " ";
        
        std::cout << " |\n";
        
    }
    
    PrintLine(r);
    
    std::cout << "Results:\n";
    
    std::cout << "J: " << J_min << std::endl << "Lmd: " << Lmd[Lmd_min] << std::endl
    
    << "W: " << Omega_min << std::endl << "D: " << Delta_min << std::endl;
    
}

int main() {
    
    int r = 3;
    
    Filter(r);
    
    r = 5;
    
    Filter(r);
    
    return 0;
    
}
