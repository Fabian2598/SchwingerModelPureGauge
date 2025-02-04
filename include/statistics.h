#ifndef STATISTICS_H_INCLUDED
#define STATISTICS_H_INCLUDED
//inline functions have to be defined in the header
//templates also have to be defined in the header. There is a solution to this problem, but I would
//have to create another file and include it here.
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>

//mean of a vector
template <typename T>
double mean(std::vector<T> x){ 
    double prom = 0;
    for (T i : x) {
        prom += i*1.0;
    }   
    prom = prom / x.size();
    return prom;
}


//absolute value
template <typename T>
inline T absVal(T z){
    if (z < 0){
        return -z;
    }
    else{
        return z;
    }
}

//random double number in the inteval [a,b] a = min, b = max
inline double rand_range(double a, double b){
    double cociente = ((double) rand() / (RAND_MAX));
    double x = (b-a) * cociente + a;
    return x;
}

//----------Jackknife---------//
std::vector<double> samples_mean(std::vector<double> dat, int bin); 
double Jackknife_error(std::vector<double> dat, int bin); 
double Jackknife(std::vector<double> dat, std::vector<int> bins); 

//---------------Linspace (similar to python)----------------------//
template <typename T>
std::vector<double> linspace(T min, T max, int n) {
    std::vector<double> linspace;
    double h = (1.0*max - 1.0*min) / (n - 1);
    for (int i = 0; i < n; ++i) {
        linspace.insert(linspace.begin() + i, min*1.0 + i * h); 
    }
    return linspace;
}
//---------------Logspace (similar to python)----------------------//
template <typename T>
std::vector<int> logspace(T min, T max, int n) {
    std::vector<int> logspace(n);
    double h = (max*1.0 - min*1.0) / (n - 1);
    for (int i = 0; i < n; ++i) {
        logspace[i] = (int) pow(10.0, min*1.0 + i * h); 
    }
    return logspace;
}

//n modulus m 
inline int modulo(int n, int m) {
    if (n < 0) {
        return (n + m) % m;
    }
    else {
        return n % m;
    }
}
//n modulus m 
inline double fmodulo(double n, double m) {
    if (n < 0) {
        return fmod(n + m, m);
    }
    else {
        return fmod(n,m);
    }
}
#endif
