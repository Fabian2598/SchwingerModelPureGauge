#include <time.h> 
#include <ctime>
#include <fstream>
#include <string>
#include "gauge_conf.h"

//For beta=1 the target value is  Ep = 0.444099 dEp = 0.00147966

int main() {
    srand(time(0));
    int Ntherm, Nmeas, Nsteps, Nbeta;
    double beta_min, beta_max;
    int algorithm;
    double trajectory_length;
    int MD_steps;
    std::string ALGORITHM[2] = { "Metropolis","HMC" };
    //---Input data---//
    std::cout << "-----------------------" << std::endl;
    std::cout << "|Pure Gauge U(1) theory|" << std::endl;
    std::cout << "-----------------------" << std::endl;
    std::cout << "Ns " << Ns << " Nt " << Nt << std::endl;
    std::cout << "Algorithm (0 Metropolis, 1 HMC): ";
    std::cin >> algorithm;
    std::cout << "----" + ALGORITHM[algorithm] + "----" << std::endl;
    if (algorithm == 1) {
        std::cout << "Molecular dynamics steps: ";
        std::cin >> MD_steps;
        std::cout << "Trajectory length: ";
        std::cin >> trajectory_length;
    }
    std::cout << "beta min: ";
    std::cin >> beta_min;
    std::cout << "beta max: ";
    std::cin >> beta_max;
    std::cout << "Number of betas: ";
    std::cin >> Nbeta;
    std::cout << "Thermalization: ";
    std::cin >> Ntherm;
    std::cout << "Measurements: ";
    std::cin >> Nmeas;
    std::cout << "Step (sweeps between measurements): ";
    std::cin >> Nsteps;
    std::cout << " " << std::endl;
    std::vector<double> Betas(Nbeta);
    GaugeConf Configuration = GaugeConf(Ns, Nt);

    if (Nbeta == 1) {
        Betas = { beta_min };
    }
    else {
        Betas = linspace(beta_min, beta_max, Nbeta);
    }
    char NameData[500], Data_str[500];
    sprintf(NameData, "2D_U1_Ns%d_Nt%d_Meas%d.txt", Ns, Nt, Nmeas);

    std::ofstream Datfile;
    Datfile.open(NameData);
	Configuration.Coordinates(); //Compute vectorized coordinates
    for (double beta : Betas) {
        clock_t begin = clock();
        std::cout << "beta = " << beta << std::endl;
        Configuration.setBeta(beta);
        if (algorithm == 0) {
            Configuration.Metropolis(Ntherm, Nmeas, Nsteps);
        }
        else if (algorithm == 1) {
            Configuration.HMC(MD_steps, trajectory_length, Ntherm, Nmeas, Nsteps);
        }
        
        sprintf(Data_str, "%-30.17g%-30.17g%-30.17g\n", beta, Configuration.getEp(), Configuration.getdEp());
        std::cout << "Ep = " << Configuration.getEp() << " dEp = " << Configuration.getdEp() << std::endl;
        Datfile << Data_str;
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        std::cout << "Time = " << elapsed_secs << " s" << std::endl;
        std::cout << "------------------------------" << std::endl;

    }
    Datfile.close();

	return 0;
}



/*
std::cout << "TESTING" << std::endl;
GaugeConf Configuration = GaugeConf(Ns, Nt);
Configuration.setBeta(1.0);
Configuration.Coordinates(); //Compute vectorized coordinates
Configuration.initialization();
SaveConf(Configuration.getConf(),"testing_conf.txt");

Configuration.Compute_Plaquette01();
int COORDINATE;
std::cout << "HMC Sp " << Configuration.MeasureSp_HMC() << std::endl;
std::cout << "Metropolis  Sp"  << Configuration.MeasureSp() << std::endl;
std::cout << std::endl;
int x = 1, t = 1, mu = 0;
std::cout << "Staple at x, t, mu: " << x << " "  << t << " "  << mu << Configuration.Staple(x, t, mu) << std::endl;
mu = 1;
std::cout << "Staple at x, t, mu: " << x <<   " " << t << " " << mu << Configuration.Staple(x, t, mu) << std::endl;
std::cout << std::endl;
mu = 0;
COORDINATE = x * Ns + t;
std::cout << "Conf at mu=0" << Configuration.getConf()[COORDINATE][mu] << std::endl;
std::cout << std::endl;
mu = 1;
std::cout << "Conf at mu=1" << Configuration.getConf()[COORDINATE][mu] << std::endl;
std::cout << "-------------" << std::endl;
Configuration.Leapfrog(5, 1.0);
//std::cout << "-------------" << std::endl;


return 0;
*/