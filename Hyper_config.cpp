#include "stdafx.h"

hyperelastic::hyperelastic()
{
	for(int i=0;i<600;i++)
	{
		NEI[i]=0;
	}
	N=0;
	pnd0=0;
	flag_wall=0;
	lambda=1;
	J=0;
	pnd=0;
	for(int D=0;D<DIMENSION;D++)
	{
		vec_norm[D]=0;
		half_p[D]=0;
		differential_p[D]=0;
		p[D]=0;	
		ang_p[D]=0;
		vis_force[D]=0;
		for(int D2=0;D2<DIMENSION;D2++)
		{
			stress[D][D2]=0;
			Ai[D][D2]=0;
			inverse_Ai[D][D2]=0;
			t_inverse_Ai[D][D2]=0;
			t_inverse_Fi[D][D2]=0;
			Fi[D][D2]=0;
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
	}
}

