#include "gauge_conf.h"
#include <iomanip>

void GaugeConf::Compute_Plaquette01() {
    for (int x = 0; x < Ns; x++) {
        for (int t = 0; t < Nt; t++) {
            int Coord0 = Coords[x][t], Coord1 = Coords[modulo(x + 1, Ns)][t], Coord2 = Coords[x][modulo(t + 1, Nt)];
            Plaquette01[Coords[x][t]] = Conf[Coord0][0] * Conf[Coord1][1] * std::conj(Conf[Coord2][0]) * std::conj(Conf[Coord0][1]);
            Plaquette01_prime[Coords[x][t]] = Conf_copy[Coord0][0] * Conf_copy[Coord1][1] * std::conj(Conf_copy[Coord2][0]) * std::conj(Conf_copy[Coord0][1]);
        }
    }
}

void GaugeConf::Force(const std::vector<std::vector<std::complex<double>>>& U) {
    StapleHMC(U); //Computes staples
	for (int x = 0; x < Ns; x++) {
		for (int t = 0; t < Nt; t++) {
			int i = Coords[x][t];
			for (int mu = 0; mu < 2; mu++) {
				Forces[i][mu] = -beta * std::imag(U[i][mu] * std::conj(staples[i][mu]));
			}
		}
	}
}

double GaugeConf::DeltaH() {
    Compute_Plaquette01();
    double deltaH = 0.0;
    double dS = 0.0;
    for (int x = 0; x < Ns; x++) {
        for (int t = 0; t < Nt; t++) {
            dS += beta * std::real(Plaquette01[Coords[x][t]] - Plaquette01_prime[Coords[x][t]]);
            for (int mu = 0; mu < 2; mu++) {
                deltaH += 0.5 * (PConf_copy[Coords[x][t]][mu] * PConf_copy[Coords[x][t]][mu] - PConf[Coords[x][t]][mu] * PConf[Coords[x][t]][mu]);
            }  
        }
    }
    deltaH += dS;
    return deltaH;
}

double GaugeConf::MeasureSp_HMC() {
    //Plaquettes have to be computed at the HMC update
    double Sp = 0.0;
    for (int x = 0; x < Ns; x++) {
        for (int t = 0; t < Nt; t++) {
            Sp += std::real(Plaquette01[Coords[x][t]]);
        }
    }
    return Sp;

}


//Generates new configuration [U,Pi]
void GaugeConf::Leapfrog(const int& MD_steps, const double& trajectory_length){
        double StepSize = trajectory_length / (MD_steps * 1.0);
        PConf_copy = PConf;
        Conf_copy = Conf;
        std::complex<double> inumber(0.0, 1.0); //imaginary number

        //Conf_copy = Conf*exp(0.5i * StepSize * PConf_copy)
        for (int x = 0; x < Ns; x++) {
            for (int t = 0; t < Nt; t++) {
                int i = Coords[x][t];
                for (int mu = 0; mu < 2; mu++) {
                    Conf_copy[i][mu] = Conf_copy[i][mu] * exp(0.5 * inumber * StepSize * PConf_copy[i][mu]);
                }
            }
        }
		Force(Conf_copy);
        for (int step = 1; step < MD_steps - 1; step++) {
            //PConf_copy += StepSize*force
            //Conf_copy *= exp(i * StepSize * PConf_copy)
            for (int x = 0; x < Ns; x++) {
                for (int t = 0; t < Nt; t++) {
                    int i = Coords[x][t];
                    for (int mu = 0; mu < 2; mu++) {
                        PConf_copy[i][mu] += StepSize *  Forces[i][mu];
                        Conf_copy[i][mu] *= exp(inumber * StepSize * PConf_copy[i][mu]);        
                    }
                }
            }
            Force(Conf_copy);
        }
        //PConf_copy += StepSize*force
        //Conf_copy = Conf*exp(0.5i * StepSize* PConf_copy)
        for (int x = 0; x < Ns; x++) {
            for (int t = 0; t < Nt; t++) {
                int i = Coords[x][t];
                for (int mu = 0; mu < 2; mu++) {
                    PConf_copy[i][mu] += StepSize * Forces[i][mu];
                    Conf_copy[i][mu] *= exp(0.5 * inumber * StepSize * PConf_copy[i][mu]);
                    
                }
            }
        }
       // SaveConf(Conf_copy,  "testing_confCopy.txt");

}

void GaugeConf::HMC_Update(const int& MD_steps, const double& trajectory_length){
    PConf = RandomMomentum();//random momentum conf sampled from a normal distribution
    Leapfrog(MD_steps, trajectory_length); //Evolve [Pi] and [U] (updates PConf_copy and Conf_copy)
    double deltaH = DeltaH(); //deltaH = Hamiltonian[U'][Pi'] - [U][Pi]
    double r = rand_range(0, 1);
    if (r<= exp(-deltaH)){
        //Accept the new configuration
        Conf = Conf_copy;
    }
    //Else configuration is not modified.
}

//HMC algorithm
void GaugeConf::HMC(const int& MD_steps, 
	const double& trajectory_length, const int& Ntherm, 
	const int& Nmeas, const int& Nsteps){

    std::vector<double> SpVector(Nmeas);
    initialization();
    for(int i = 0; i < Ntherm; i++) {HMC_Update(MD_steps,trajectory_length);} //Thermalization
    for(int i = 0; i < Nmeas; i++) {
        HMC_Update(MD_steps,trajectory_length);
        SpVector[i] = MeasureSp_HMC(); 
        for(int j = 0; j < Nsteps; j++) {HMC_Update(MD_steps,trajectory_length);} //Decorrelation
    }
    Ep = mean(SpVector) / (Ntot * 1.0); dEp = Jackknife_error(SpVector, 20) / (Ntot * 1.0);
} 
