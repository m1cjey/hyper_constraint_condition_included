#ifndef HYPERELASTIC
#define HYPERELASTIC

class hyperelastic
{
public:

	unsigned ID;
	int NEI[1000];
	int N;
	double pnd0;
	int flag_wall;

	double dec_ET;
	double lambda;
	double lam_w;
	double half_p[DIMENSION];
	double stress[DIMENSION][DIMENSION];
	double differential_p[DIMENSION];
	double p[DIMENSION];
	double ang_p[DIMENSION];
	double Ai[DIMENSION][DIMENSION];
	double inverse_Ai[DIMENSION][DIMENSION];
	double t_inverse_Ai[DIMENSION][DIMENSION];
	double t_inverse_Fi[DIMENSION][DIMENSION];
	double J;
	double pnd;
	double Fi[DIMENSION][DIMENSION];
	double vis_force[DIMENSION];
	double vec_norm[DIMENSION];
	double Fw[DIMENSION];
	hyperelastic();
};

class hyperelastic2
{
public:
	unsigned ID;
	double wiin;
	double DgDq[DIMENSION];
	double DgDq_w[DIMENSION];
	double aiin[DIMENSION];
	double n0ij[DIMENSION];
	double spl_f;

	hyperelastic2();
};

#endif