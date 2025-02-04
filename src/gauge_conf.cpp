#include "gauge_conf.h"

//Random U1 variable
std::complex<double> RandomU1() {
	double cociente = ((double) rand() / (RAND_MAX));
    double theta = 2.0*pi * cociente;
	std::complex<double> z(cos(theta), sin(theta));
	return z;
}


//Initialize a random conf
void GaugeConf::initialization() {
	for (int i = 0; i < Ntot; i++) {
		for (int mu = 0; mu < 2; mu++) {
			Conf[i][mu] = RandomU1(); //Conf[Ns x Nt][mu in {0,1}]
		}
	}
}


void GaugeConf::printPlaquette() {
	for (int i = 0; i < Ntot; i++) {
		std::cout << Plaquette01[i] << " ";
	}
}


void GaugeConf::Coordinates() {
	for (int x = 0; x < Ns; x++) {
		for (int t = 0; t < Nt; t++) {
			Coords[x][t] = Coordinate(x, t);
		}
	}
}

//Compute staple at coordinate (x,t) in the mu-direction
std::complex<double> GaugeConf::Staple(const int& x, const int& t, const int& mu) {
	//x and t --> lattice sites. x: rows, t: columns
	// WARNING: Some references define the staple as the conjugate of this:
	//U_v(x) U_m(x+v) U*_v(x+m) + U*_v(x-v) U_m(x-v) U_v(x+m-v)
	int x1 = modulo(x + 1, Ns);
	int x_1 = modulo(x - 1, Ns);
	int t1 = modulo(t + 1, Nt);
	int t_1 = modulo(t - 1, Nt);

	if (mu == 0) {
		const std::complex<double>& conf1 = Conf[Coords[x][t]][1];
		const std::complex<double>& conf2 = Conf[Coords[x][t1]][0];
		const std::complex<double>& conf3 = Conf[Coords[x1][t]][1];
		const std::complex<double>& conf4 = Conf[Coords[x][t_1]][1];
		const std::complex<double>& conf5 = Conf[Coords[x][t_1]][0];
		const std::complex<double>& conf6 = Conf[Coords[x1][t_1]][1];
		return  conf1 * conf2 * std::conj(conf3) +
			std::conj(conf4) * conf5 * conf6;
	}
	else {
		const std::complex<double>& conf1 = Conf[Coords[x][t]][0];
		const std::complex<double>& conf2 = Conf[Coords[x1][t]][1];
		const std::complex<double>& conf3 = Conf[Coords[x][t1]][0];
		const std::complex<double>& conf4 = Conf[Coords[x_1][t]][0];
		const std::complex<double>& conf5 = Conf[Coords[x_1][t]][1];
		const std::complex<double>& conf6 = Conf[Coords[x_1][t1]][0];
		return conf1 * conf2 * std::conj(conf3) +
			std::conj(conf4) * conf5 * conf6;
	}
}

//Compute staple at coordinate (x,t) in the mu-direction
std::complex<double> GaugeConf::StapleHMC(const std::vector<std::vector<std::complex<double>>>& U,const int& x, const int& t, const int& mu) {
	//x and t --> lattice sites. x: rows, t: columns
	// WARNING: Some references define the staple as the conjugate of this:
	//U_v(x) U_m(x+v) U*_v(x+m) + U*_v(x-v) U_m(x-v) U_v(x+m-v)
	int x1 = modulo(x + 1, Ns);
	int x_1 = modulo(x - 1, Ns);
	int t1 = modulo(t + 1, Nt);
	int t_1 = modulo(t - 1, Nt);

	if (mu == 0) {
		const std::complex<double>& conf1 = U[Coords[x][t]][1];
		const std::complex<double>& conf2 = U[Coords[x][t1]][0];
		const std::complex<double>& conf3 = U[Coords[x1][t]][1];
		const std::complex<double>& conf4 = U[Coords[x][t_1]][1];
		const std::complex<double>& conf5 = U[Coords[x][t_1]][0];
		const std::complex<double>& conf6 = U[Coords[x1][t_1]][1];
		return  conf1 * conf2 * std::conj(conf3) +
			std::conj(conf4) * conf5 * conf6;
	}
	else {
		const std::complex<double>& conf1 = U[Coords[x][t]][0];
		const std::complex<double>& conf2 = U[Coords[x1][t]][1];
		const std::complex<double>& conf3 = U[Coords[x][t1]][0];
		const std::complex<double>& conf4 = U[Coords[x_1][t]][0];
		const std::complex<double>& conf5 = U[Coords[x_1][t]][1];
		const std::complex<double>& conf6 = U[Coords[x_1][t1]][0];
		return conf1 * conf2 * std::conj(conf3) +
			std::conj(conf4) * conf5 * conf6;
	}
}

//Print conf
std::ostream& operator<<(std::ostream& out, const GaugeConf& GConf) {
	std::vector<std::vector<std::complex<double>>> Conf = GConf.getConf();
	for (int i = 0; i < GConf.getNs(); i++) {
		for (int j = 0; j < GConf.getNt(); j++) {
			int N0 = GConf.Coordinate(i, j);
			out << "[" << Conf[N0][0] << ", " << Conf[N0][1] << "] "; //Conf[Ns x Nt][mu in {0,1}]
		}
		out << "\n";
	}
	return out;
}

