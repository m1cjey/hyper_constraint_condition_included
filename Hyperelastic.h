#ifndef HYPERELASTIC
#define HYPERELASTIC

class hyperelastic
{
public:

	unsigned ID;
	int NEI0[200];
	int N0;
	double pnd;
	
	double lambda;
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
	double Fi[DIMENSION][DIMENSION];
	double vis_force[DIMENSION];
	double vec_norm[DIMENSION];
	
	hyperelastic();
};

class hyperelastic2
{
public:
	unsigned ID;
	double wiin;
	double DgDq[DIMENSION];
	double aiin[DIMENSION];
	double n0ij[DIMENSION];
	double spl_f;
	double DwDq[DIMENSION];

	hyperelastic2();
};

#endif