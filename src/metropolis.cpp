#include "gauge_conf.h"

//Change in the action
double GaugeConf::DeltaS(const std::complex<double>& U_prime, const std::complex<double>& U, const std::complex<double>& Sigma) {
	return -beta * std::real( (U_prime - U) * std::conj(Sigma) );
}

//Updates one site
void GaugeConf::LocalUpdate(const int& x, const int& t, const int& mu) {
	std::complex<double> U_prime = RandomU1();
	//Sigma is the staple
	double	r = ((double) rand() / (RAND_MAX)); //rand_range(0, 1); 
	double deltaS = DeltaS(U_prime, Conf[Coords[x][t]][mu], Staple(x, t, mu)); 	//Staple defined on gauge_conf.cpp
	double p = exp(-deltaS);
	if (r <= p) {
		Conf[Coords[x][t]][mu] = U_prime;
	}
}

//Updates the whole lattice
void GaugeConf::Sweep() {
	for (int x = 0; x < Ns; x++) {
		for (int t = 0; t < Nt; t++) {
			for (int mu = 0; mu < 2; mu++) {
				LocalUpdate(x, t, mu);
			}
		}
	}
}


double GaugeConf::MeasureSp() {
	//Plaquette U_mv(x) = U_m(x) U_v(x+m) U*_m(x+v) U*_v(x)
	//Sp = 1/Ntot Sum_x Sum(m<v) U_mv(x) (In two dimensions this is just 1/Ntot Sum_x U_01(x)
	std::complex<double> Plaquette;
	double Sp = 0.0;
	for (int x = 0; x < Ns; x++) {
		for (int t = 0; t < Nt; t++) {
			int Coord0, Coord1, Coord2;
			Coord0 = Coords[x][t];
			Coord1 = Coords[modulo(x + 1, Ns)][t]; Coord2 = Coords[x][modulo(t + 1, Nt)];
			Plaquette = Conf[Coord0][0] * Conf[Coord1][1] * std::conj(Conf[Coord2][0]) * std::conj(Conf[Coord0][1]);
			Sp += std::real(Plaquette);
		}
	}
	return Sp;
}

	
void GaugeConf::Metropolis(const int& Ntherm, const int& Nmeas,const int& Nstep) {
	std::vector<double> SpVector(Nmeas);
	initialization();
	for (int i = 0; i < Ntherm; i++) { Sweep(); } //Thermalization

	for (int i = 0; i < Nmeas; i++) {
		Sweep();
		SpVector[i] = MeasureSp(); 
		//Decorrelation sweeps
		for (int j = 0; j < Nstep; j++) {Sweep();}
	}
	Ep = mean(SpVector) / (Ntot * 1.0); dEp = Jackknife_error(SpVector, 20) / (Ntot * 1.0);
}

	

