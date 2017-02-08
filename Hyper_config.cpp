#include "stdafx.h"

hyperelastic::hyperelastic()
{
	for(int i=0;i<60;i++)
	{
		NEI[i]=0;
	}
	N=0;
	pnd0=0;
	W=0;
	lambda=1;
	J=0;
	pnd=0;

	lam=1;
	mu=1;

	old_lam=1;
	old_mu=1;

	
	W_n=0;
	fw=0;
	W0=0;
	Ef=0.;
	E0=0.;

	for(int D=0;D<DIMENSION;D++)
	{
		half_p[D]=0;
		differential_p[D]=0;
		p[D]=0;	
		ang_p[D]=0;
		vis_force[D]=0;
		p0[D]=0;
		q_n[D]=0;
		p_n[D]=0;
		ph_n[D]=0;
		for(int D2=0;D2<DIMENSION;D2++)
		{
			stress[D][D2]=0;
			Ai[D][D2]=0;
			inverse_Ai[D][D2]=0;
			t_inverse_Ai[D][D2]=0;
			t_inverse_Fi[D][D2]=0;
			Fi[D][D2]=0;

			//pi[D][D2]=0;
			//pi_n[D][D2]=0;
			stress_n[D][D2]=0;
		}
	}
}

hyperelastic2::hyperelastic2()
{
	wiin=0;
	spl_f=0;
	for(int D=0;D<DIMENSION;D++)
	{
		DgDq[D]=0;
		aiin[D]=0;
		n0ij[D]=0;

		DgDq_n[D]=0;
		DpiDq[D]=0.;
	}
}

