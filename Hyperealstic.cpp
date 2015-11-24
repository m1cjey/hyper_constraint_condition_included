#include "stdafx.h"		

void calc_hyper(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int rigid_number,int t)
{	
	int h_num=HYPER.size();
	int r_num=rigid_number;
	int flag_w=OFF;
	double Dt=CON.get_dt();
	cout<<"Hypercalculation starts."<<endl;

	//calc_gravity(CON,HYPER,h_num);
	
	if(t==1)
	{
		for(int i=0;i<h_num;i++)	for(int D=0;D<DIMENSION;D++)	PART[i].q0[D]=0;
		for(int i=0;i<h_num;i++)	for(int D=0;D<DIMENSION;D++)	PART[i].q0[D]=PART[i].r[D];
		calc_constant(CON,PART,HYPER,HYPER1,r_num);
		calc_stress(CON,HYPER,r_num);
	}
	
	if(t==1 || t%CON.get_interval()==0)
	{
		output_hyper_data(PART,HYPER,HYPER1,r_num,t);
		momentum_movie_AVS(CON,t,PART,HYPER,r_num);
		output_energy(CON,PART,HYPER,r_num,t);
	}

	/*
	if(t!=1)//for tension test
	{
		double Dt=CON.get_dt();
		for(int i=0;i<h_num;i++)	HYPER[i].p[A_Z]+=10*Dt;
	}*/

	newton_raphson(CON,PART,HYPER,HYPER1,r_num,flag_w,t);

	calc_half_p(CON,PART,HYPER,HYPER1,0,r_num,flag_w,t);

	calc_F(CON,PART,HYPER,HYPER1,r_num,flag_w,t);

	calc_stress(CON,HYPER,r_num);
	
	calc_differential_p(CON,HYPER,HYPER1,r_num);

	renew_lambda(CON,HYPER,HYPER1,r_num,flag_w,t);

	calc_half_p(CON,PART,HYPER,HYPER1,1,r_num,flag_w,t);

	cout<<"Hypercalculation ends."<<endl;

}

void calc_constant(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int rigid_number)
{
	cout<<"�����l�v�Z";

	double le=CON.get_distancebp();
	double r=CON.get_h_dis();
	double Dt=CON.get_dt();
	double V=get_volume(&CON);	//�l�����K�v��������Ȃ�
	double mi=V*CON.get_hyper_density();
	int h_num=HYPER.size();
	int r_num=rigid_number;
	int model=CON.get_model_number();

	////�����^����
	//�Ȃ��˂���
	/*if(model==21)
	{		
		int t=30,b=2;
		double max=0,min=0;
		for(int i=0;i<h_num;i++)
		{
			if(max<PART[i].q0[A_Z])	max=PART[i].q0[A_Z];
			if(min>PART[i].q0[A_Z])	min=PART[i].q0[A_Z];
		}
		double H=max-min+le;
		cout<<H;
		//double H=1.8;
		for(int i=0;i<h_num;i++)	
		{
			double Z=PART[i].q0[A_Z];
			double Y=PART[i].q0[A_Y];
			double X=PART[i].q0[A_X];
			double part_p=(Z/H)*2;
			HYPER[i].p[A_X]=mi*(-t*part_p*part_p*part_p*Y+b*(3*part_p*part_p-1));
			HYPER[i].p[A_Y]=mi*t*part_p*part_p*part_p*X;
			HYPER[i].p[A_Z]=0;
		}
	}*/

	//�Ȃ�
	/*if(model==21)
	{
		int b=2;
		double max=0,min=0;
		for(int i=0;i<h_num;i++)
		{
			if(max<PART[i].q0[A_Z])	max=PART[i].q0[A_Z];
			if(min>PART[i].q0[A_Z])	min=PART[i].q0[A_Z];
		}
		double H=max-min+le;
		cout<<H;
		//double H=1.8;
		for(int i=0;i<h_num;i++)	
		{
			double Z=PART[i].q0[A_Z];
			double part_p=(Z/H)*2;
			HYPER[i].p[A_X]=mi*b*(3*part_p*part_p-1);
		}
	}//*/

	//�˂���
	/*if(model==21)
	{
		int t=30;
		double max=0,min=0;
		for(int i=0;i<h_num;i++)
		{
			if(max<PART[i].q0[A_Z])	max=PART[i].q0[A_Z];
			if(min>PART[i].q0[A_Z])	min=PART[i].q0[A_Z];
		}
		double H=max-min+le;
		cout<<H;
		//double H=1.8;
		for(int i=0;i<h_num;i++)	
		{
			double Z=PART[i].q0[A_Z];
			double Y=PART[i].q0[A_Y];
			double X=PART[i].q0[A_X];
			double part_p=(Z/H)*2;
			HYPER[i].p[A_X]=-mi*t*part_p*part_p*part_p*Y;
			HYPER[i].p[A_Y]=mi*t*part_p*part_p*part_p*X;
			HYPER[i].p[A_Z]=0;
		}
	}//*/


	//��]
	if(model==22)
	{
		for(int i=0;i<h_num-r_num;i++)
		{
			PART[i].p[A_X]=mi*0.4*(PART[i].q0[A_Z]-PART[i].q0[A_Y]);
			PART[i].p[A_Y]=mi*0.4*(PART[i].q0[A_X]-PART[i].q0[A_Z]);
			PART[i].p[A_Z]=mi*0.4*(PART[i].q0[A_Y]-PART[i].q0[A_X]);
		}
	}

	////�p�^���ʌv�Z
	for(int i=0;i<h_num;i++)
	{
		HYPER[i].ang_p[A_X]=PART[i].r[A_Y]*HYPER[i].p[A_Z]-PART[i].r[A_Z]*HYPER[i].p[A_Y];
		HYPER[i].ang_p[A_Y]=PART[i].r[A_Z]*HYPER[i].p[A_X]-PART[i].r[A_X]*HYPER[i].p[A_Z];
		HYPER[i].ang_p[A_Z]=PART[i].r[A_X]*HYPER[i].p[A_Y]-PART[i].r[A_Y]*HYPER[i].p[A_X];
	}


	////�ߖT���q�̋L����aiin,wiin,Ai�̌v�Z
	for(int i=0;i<h_num-r_num;i++)
	{
		int N0=0;
		double dis=0;
		double wiin=0;
		double aiin[DIMENSION];
		for(int j=0;j<h_num-r_num;j++)
		{
			wiin=0;
			aiin[A_X]=PART[j].q0[A_X]-PART[i].q0[A_X];	aiin[A_Y]=PART[j].q0[A_Y]-PART[i].q0[A_Y];	aiin[A_Z]=PART[j].q0[A_Z]-PART[i].q0[A_Z];
			
			HYPER1[i*h_num+j].aiin[A_X]=aiin[A_X];
			HYPER1[i*h_num+j].aiin[A_Y]=aiin[A_Y];
			HYPER1[i*h_num+j].aiin[A_Z]=aiin[A_Z];

			dis=sqrt(aiin[A_X]*aiin[A_X]+aiin[A_Y]*aiin[A_Y]+aiin[A_Z]*aiin[A_Z]);
			if(dis<r && j!=i)
			{	
				wiin=kernel4(r,dis);
				HYPER[i].NEI0[N0]=j;
				N0++;
			}
			else	wiin=0;

			HYPER[i].pnd+=wiin;
			HYPER1[i*h_num+j].wiin=wiin;
		}
		HYPER[i].N0=N0;
	}

	for(int k=h_num-r_num;k<h_num;k++)
	{
		int N0=0;
		double dis=0;
		double wiin=0;
		double aiin[DIMENSION];
		for(int l=h_num-r_num;l<h_num;l++)
		{
			wiin=0;
			aiin[A_X]=PART[l].q0[A_X]-PART[k].q0[A_X];	aiin[A_Y]=PART[l].q0[A_Y]-PART[k].q0[A_Y];	aiin[A_Z]=PART[l].q0[A_Z]-PART[k].q0[A_Z];
			
			HYPER1[k*h_num+l].aiin[A_X]=aiin[A_X];
			HYPER1[k*h_num+l].aiin[A_Y]=aiin[A_Y];
			HYPER1[k*h_num+l].aiin[A_Z]=aiin[A_Z];

			dis=sqrt(aiin[A_X]*aiin[A_X]+aiin[A_Y]*aiin[A_Y]+aiin[A_Z]*aiin[A_Z]);
			if(dis<r && l!=k)
			{	
				wiin=kernel4(r,dis);
				HYPER[k].NEI0[N0]=l;
				N0++;
			}
			else	wiin=0;

			HYPER[k].pnd+=wiin;
			HYPER1[k*h_num+l].wiin=wiin;
		}
		HYPER[k].N0=N0;
	}
	
	
	////Ai, Fi�֘A�̌v�Z
	double **p_Ai=new double *[DIMENSION];
	double **p_Fi=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)
	{
		p_Ai[D]=new double [DIMENSION];
		p_Fi[D]=new double[DIMENSION];
	}
	for(int i=0;i<h_num;i++)
	{
		//Ai�̌v�Z
		double Ai[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

		int Ni=HYPER[i].N0;
		for(int j=0;j<Ni;j++)	
		{
			int k=HYPER[i].NEI0[j];
			double w=HYPER1[i*h_num+k].wiin;		
			double a[DIMENSION]={HYPER1[i*h_num+k].aiin[A_X],	HYPER1[i*h_num+k].aiin[A_Y], HYPER1[i*h_num+k].aiin[A_Z]};
			
			Ai[0][0]+=w*a[0]*a[0];		Ai[0][1]+=w*a[0]*a[1];		Ai[0][2]+=w*a[0]*a[2];
			Ai[1][0]+=w*a[1]*a[0];		Ai[1][1]+=w*a[1]*a[1];		Ai[1][2]+=w*a[1]*a[2];
			Ai[2][0]+=w*a[2]*a[0];		Ai[2][1]+=w*a[2]*a[1];		Ai[2][2]+=w*a[2]*a[2];
		}
		HYPER[i].Ai[0][0]=Ai[0][0];	HYPER[i].Ai[0][1]=Ai[0][1];	HYPER[i].Ai[0][2]=Ai[0][2];
		HYPER[i].Ai[1][0]=Ai[1][0];	HYPER[i].Ai[1][1]=Ai[1][1];	HYPER[i].Ai[1][2]=Ai[1][2];
		HYPER[i].Ai[2][0]=Ai[2][0];	HYPER[i].Ai[2][1]=Ai[2][1];	HYPER[i].Ai[2][2]=Ai[2][2];
		
		//inverse_Ai,t_inverse_Ai�̌v�Z
		p_Ai[0][0]=Ai[0][0];	p_Ai[0][1]=Ai[0][1];	p_Ai[0][2]=Ai[0][2];
		p_Ai[1][0]=Ai[1][0];	p_Ai[1][1]=Ai[1][1];	p_Ai[1][2]=Ai[1][2];
		p_Ai[2][0]=Ai[2][0];	p_Ai[2][1]=Ai[2][1];	p_Ai[2][2]=Ai[2][2];
		inverse(p_Ai,DIMENSION);
		HYPER[i].inverse_Ai[0][0]=p_Ai[0][0];		HYPER[i].inverse_Ai[0][1]=p_Ai[0][1];		HYPER[i].inverse_Ai[0][2]=p_Ai[0][2];
		HYPER[i].inverse_Ai[1][0]=p_Ai[1][0];		HYPER[i].inverse_Ai[1][1]=p_Ai[1][1];		HYPER[i].inverse_Ai[1][2]=p_Ai[1][2];
		HYPER[i].inverse_Ai[2][0]=p_Ai[2][0];		HYPER[i].inverse_Ai[2][1]=p_Ai[2][1];		HYPER[i].inverse_Ai[2][2]=p_Ai[2][2];

		HYPER[i].t_inverse_Ai[0][0]=p_Ai[0][0];		HYPER[i].t_inverse_Ai[0][1]=p_Ai[1][0];		HYPER[i].t_inverse_Ai[0][2]=p_Ai[2][0];
		HYPER[i].t_inverse_Ai[1][0]=p_Ai[0][1];		HYPER[i].t_inverse_Ai[1][1]=p_Ai[1][1];		HYPER[i].t_inverse_Ai[1][2]=p_Ai[2][1];
		HYPER[i].t_inverse_Ai[2][0]=p_Ai[0][2];		HYPER[i].t_inverse_Ai[2][1]=p_Ai[1][2];		HYPER[i].t_inverse_Ai[2][2]=p_Ai[2][2];
		
		//Fi�̌v�Z
		p_Fi[0][0]=Ai[0][0]*p_Ai[0][0]+Ai[0][1]*p_Ai[1][0]+Ai[0][2]*p_Ai[2][0];	p_Fi[0][1]=Ai[0][0]*p_Ai[0][1]+Ai[0][1]*p_Ai[1][1]+Ai[0][2]*p_Ai[2][1];	p_Fi[0][2]=Ai[0][0]*p_Ai[0][2]+Ai[0][1]*p_Ai[1][2]+Ai[0][2]*p_Ai[2][2];
		p_Fi[1][0]=Ai[1][0]*p_Ai[0][0]+Ai[1][1]*p_Ai[1][0]+Ai[1][2]*p_Ai[2][0];	p_Fi[1][1]=Ai[1][0]*p_Ai[0][1]+Ai[1][1]*p_Ai[1][1]+Ai[1][2]*p_Ai[2][1];	p_Fi[1][2]=Ai[1][0]*p_Ai[0][2]+Ai[1][1]*p_Ai[1][2]+Ai[1][2]*p_Ai[2][2];
		p_Fi[2][0]=Ai[2][0]*p_Ai[0][0]+Ai[2][1]*p_Ai[1][0]+Ai[2][2]*p_Ai[2][0];	p_Fi[2][1]=Ai[2][0]*p_Ai[0][1]+Ai[2][1]*p_Ai[1][1]+Ai[2][2]*p_Ai[2][1];	p_Fi[2][2]=Ai[2][0]*p_Ai[0][2]+Ai[2][1]*p_Ai[1][2]+Ai[2][2]*p_Ai[2][2];		

		HYPER[i].Fi[0][0]=p_Fi[0][0];	HYPER[i].Fi[0][1]=p_Fi[0][1];	HYPER[i].Fi[0][2]=p_Fi[0][2];	
		HYPER[i].Fi[1][0]=p_Fi[1][0];	HYPER[i].Fi[1][1]=p_Fi[1][1];	HYPER[i].Fi[1][2]=p_Fi[1][2];	
		HYPER[i].Fi[2][0]=p_Fi[2][0];	HYPER[i].Fi[2][1]=p_Fi[2][1];	HYPER[i].Fi[2][2]=p_Fi[2][2];	
		//J�̌v�Z
		double J=calc_det3(p_Fi);
		HYPER[i].J=J;
//		cout<<"HYPER["<<i<<"].J="<<HYPER[i].J<<endl;

		//t_inverse_Fi�̌v�Z
		inverse(p_Fi,DIMENSION);
		HYPER[i].t_inverse_Fi[0][0]=p_Fi[0][0];		HYPER[i].t_inverse_Fi[0][1]=p_Fi[1][0];		HYPER[i].t_inverse_Fi[0][2]=p_Fi[2][0];
		HYPER[i].t_inverse_Fi[1][0]=p_Fi[0][1];		HYPER[i].t_inverse_Fi[1][1]=p_Fi[1][1];		HYPER[i].t_inverse_Fi[1][2]=p_Fi[2][1];
		HYPER[i].t_inverse_Fi[2][0]=p_Fi[0][2];		HYPER[i].t_inverse_Fi[2][1]=p_Fi[1][2];		HYPER[i].t_inverse_Fi[2][2]=p_Fi[2][2];
	}
	
	for(int D=0;D<DIMENSION;D++)
	{
		delete[]	p_Ai[D];
		delete[]	p_Fi[D];
	}
	delete[]	p_Ai;
	delete[]	p_Fi;


	////n0ij�̌v�Z
	for(int i=0;i<h_num;i++)
	{
		int Ni=HYPER[i].N0;
		double p_n0ij[DIMENSION]={0,0,0};
		for(int j=0;j<Ni;j++)
		{
			int k=HYPER[i].NEI0[j];			
			p_n0ij[A_X]+=HYPER1[k*h_num+i].wiin*HYPER1[i*h_num+k].aiin[A_X];
			p_n0ij[A_Y]+=HYPER1[k*h_num+i].wiin*HYPER1[i*h_num+k].aiin[A_Y];
			p_n0ij[A_Z]+=HYPER1[k*h_num+i].wiin*HYPER1[i*h_num+k].aiin[A_Z];

			HYPER1[i*h_num+k].n0ij[A_X]=V*HYPER1[k*h_num+i].wiin*(HYPER[k].t_inverse_Ai[A_X][0]*HYPER1[i*h_num+k].aiin[0]+HYPER[k].t_inverse_Ai[A_X][1]*HYPER1[i*h_num+k].aiin[1]+HYPER[k].t_inverse_Ai[A_X][2]*HYPER1[i*h_num+k].aiin[2]);
			HYPER1[i*h_num+k].n0ij[A_Y]=V*HYPER1[k*h_num+i].wiin*(HYPER[k].t_inverse_Ai[A_Y][0]*HYPER1[i*h_num+k].aiin[0]+HYPER[k].t_inverse_Ai[A_Y][1]*HYPER1[i*h_num+k].aiin[1]+HYPER[k].t_inverse_Ai[A_Y][2]*HYPER1[i*h_num+k].aiin[2]);
			HYPER1[i*h_num+k].n0ij[A_Z]=V*HYPER1[k*h_num+i].wiin*(HYPER[k].t_inverse_Ai[A_Z][0]*HYPER1[i*h_num+k].aiin[0]+HYPER[k].t_inverse_Ai[A_Z][1]*HYPER1[i*h_num+k].aiin[1]+HYPER[k].t_inverse_Ai[A_Z][2]*HYPER1[i*h_num+k].aiin[2]);
		}
		HYPER1[i*h_num+i].n0ij[A_X]=V*(HYPER[i].t_inverse_Ai[A_X][0]*p_n0ij[0]+HYPER[i].t_inverse_Ai[A_X][1]*p_n0ij[1]+HYPER[i].t_inverse_Ai[A_X][2]*p_n0ij[2]);
		HYPER1[i*h_num+i].n0ij[A_Y]=V*(HYPER[i].t_inverse_Ai[A_Y][0]*p_n0ij[0]+HYPER[i].t_inverse_Ai[A_Y][1]*p_n0ij[1]+HYPER[i].t_inverse_Ai[A_Y][2]*p_n0ij[2]);
		HYPER1[i*h_num+i].n0ij[A_Z]=V*(HYPER[i].t_inverse_Ai[A_Z][0]*p_n0ij[0]+HYPER[i].t_inverse_Ai[A_Z][1]*p_n0ij[1]+HYPER[i].t_inverse_Ai[A_Z][2]*p_n0ij[2]);
	}

	////DgDq�̌v�Z
	for(int i=0;i<h_num-r_num;i++)
	{
		int Ni=HYPER[i].N0;
		for(int j=0;j<Ni;j++)
		{
			int k=HYPER[i].NEI0[j];
			HYPER1[k*h_num+i].DgDq[A_X]=HYPER[k].J*(HYPER[k].t_inverse_Fi[A_X][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].t_inverse_Fi[A_X][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].t_inverse_Fi[A_X][2]*HYPER1[i*h_num+k].n0ij[2]);
			HYPER1[k*h_num+i].DgDq[A_Y]=HYPER[k].J*(HYPER[k].t_inverse_Fi[A_Y][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].t_inverse_Fi[A_Y][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].t_inverse_Fi[A_Y][2]*HYPER1[i*h_num+k].n0ij[2]);
			HYPER1[k*h_num+i].DgDq[A_Z]=HYPER[k].J*(HYPER[k].t_inverse_Fi[A_Z][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].t_inverse_Fi[A_Z][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].t_inverse_Fi[A_Z][2]*HYPER1[i*h_num+k].n0ij[2]);		
		}
		HYPER1[i*h_num+i].DgDq[A_X]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_X][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_X][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_X][2]*HYPER1[i*h_num+i].n0ij[2]);
		HYPER1[i*h_num+i].DgDq[A_Y]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_Y][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_Y][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_Y][2]*HYPER1[i*h_num+i].n0ij[2]);
		HYPER1[i*h_num+i].DgDq[A_Z]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_Z][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_Z][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_Z][2]*HYPER1[i*h_num+i].n0ij[2]);
		
		double DwDq[DIMENSION]={0,0,0};
		double n0=HYPER[i].pnd;
		for(int k=h_num-r_num;k<h_num;k++)
		{
			double qk[DIMENSION]={PART[k].r[A_X],PART[k].r[A_Y],PART[k].r[A_Z]};
			double qi[DIMENSION]={PART[i].r[A_X],PART[i].r[A_Y],PART[i].r[A_Z]};
			double qik[DIMENSION]={qk[A_X]-qi[A_X],	qk[A_Y]-qi[A_Y],	qk[A_Z]-qi[A_Z]};
			double dis=sqrt(qik[A_X]*qik[A_X]+qik[A_Y]*qik[A_Y]+qik[A_Z]*qik[A_Z]);
			if(dis<r)
			{
				if(qik[A_X]<0)
				{
					DwDq[A_X]+=r/(qi[A_X]*qi[A_X]);
					HYPER1[i*h_num+k].DgDq[A_X]=-1/n0*r/(qk[A_X]*qk[A_X]);
				}
				if(qik[A_X]>0)
				{
					DwDq[A_X]+=-r/(qi[A_X]*qi[A_X]);
					HYPER1[i*h_num+k].DgDq[A_X]=1/n0*r/(qk[A_X]*qk[A_X]);
				}
				if(qik[A_Y]<0)
				{
					DwDq[A_Y]+=r/(qi[A_Y]*qi[A_Y]);
					HYPER1[i*h_num+k].DgDq[A_Y]=-1/n0*r/(qk[A_Y]*qk[A_Y]);
				}
				if(qik[A_Y]>0)
				{
					DwDq[A_Y]+=-r/(qi[A_Y]*qi[A_Y]);
					HYPER1[i*h_num+k].DgDq[A_Y]=1/n0*r/(qk[A_Y]*qk[A_Y]);
				}
				if(qik[A_Z]<0)
				{
					DwDq[A_Z]+=r/(qi[A_Z]*qi[A_Z]);
					HYPER1[i*h_num+k].DgDq[A_Z]=-1/n0*r/(qk[A_Z]*qk[A_Z]);
				}
				if(qik[A_Z]>0)
				{
					DwDq[A_Z]+=-r/(qi[A_Z]*qi[A_Z]);					
					HYPER1[i*h_num+k].DgDq[A_Z]=1/n0*r/(qk[A_Z]*qk[A_Z]);
				}
			}
		}
		HYPER1[i*h_num+i].DgDq[A_X]+=1/n0*DwDq[A_X];
		HYPER1[i*h_num+i].DgDq[A_X]+=1/n0*DwDq[A_Y];
		HYPER1[i*h_num+i].DgDq[A_X]+=1/n0*DwDq[A_Z];
	}

	for(int l=h_num-r_num;l<h_num;l++)
	{
		int Nl=HYPER[l].N0;
		for(int ln=0;ln<Nl;ln++)
		{
			int lnn=HYPER[l].NEI0[ln];
			HYPER1[lnn*h_num+l].DgDq[A_X]=HYPER[lnn].J*(HYPER[lnn].t_inverse_Fi[A_X][0]*HYPER1[l*h_num+lnn].n0ij[0]+HYPER[lnn].t_inverse_Fi[A_X][1]*HYPER1[l*h_num+lnn].n0ij[1]+HYPER[lnn].t_inverse_Fi[A_X][2]*HYPER1[l*h_num+lnn].n0ij[2]);
			HYPER1[lnn*h_num+l].DgDq[A_Y]=HYPER[lnn].J*(HYPER[lnn].t_inverse_Fi[A_Y][0]*HYPER1[l*h_num+lnn].n0ij[0]+HYPER[lnn].t_inverse_Fi[A_Y][1]*HYPER1[l*h_num+lnn].n0ij[1]+HYPER[lnn].t_inverse_Fi[A_Y][2]*HYPER1[l*h_num+lnn].n0ij[2]);
			HYPER1[lnn*h_num+l].DgDq[A_Z]=HYPER[lnn].J*(HYPER[lnn].t_inverse_Fi[A_Z][0]*HYPER1[l*h_num+lnn].n0ij[0]+HYPER[lnn].t_inverse_Fi[A_Z][1]*HYPER1[l*h_num+lnn].n0ij[1]+HYPER[lnn].t_inverse_Fi[A_Z][2]*HYPER1[l*h_num+lnn].n0ij[2]);		
		}
		HYPER1[l*h_num+l].DgDq[A_X]=HYPER[l].J*(HYPER[l].t_inverse_Fi[A_X][0]*HYPER1[l*h_num+l].n0ij[0]+HYPER[l].t_inverse_Fi[A_X][1]*HYPER1[l*h_num+l].n0ij[1]+HYPER[l].t_inverse_Fi[A_X][2]*HYPER1[l*h_num+l].n0ij[2]);
		HYPER1[l*h_num+l].DgDq[A_Y]=HYPER[l].J*(HYPER[l].t_inverse_Fi[A_Y][0]*HYPER1[l*h_num+l].n0ij[0]+HYPER[l].t_inverse_Fi[A_Y][1]*HYPER1[l*h_num+l].n0ij[1]+HYPER[l].t_inverse_Fi[A_Y][2]*HYPER1[l*h_num+l].n0ij[2]);
		HYPER1[l*h_num+l].DgDq[A_Z]=HYPER[l].J*(HYPER[l].t_inverse_Fi[A_Z][0]*HYPER1[l*h_num+l].n0ij[0]+HYPER[l].t_inverse_Fi[A_Z][1]*HYPER1[l*h_num+l].n0ij[1]+HYPER[l].t_inverse_Fi[A_Z][2]*HYPER1[l*h_num+l].n0ij[2]);
		
		double n0=HYPER[l].pnd;
		double DwDq[DIMENSION]={0,0,0};
		for(int i=0;i<h_num-r_num;i++)
		{
			double qi[DIMENSION]={PART[i].r[A_X],PART[i].r[A_Y],PART[i].r[A_Z]};
			double ql[DIMENSION]={PART[l].r[A_X],PART[l].r[A_Y],PART[l].r[A_Z]};
			double qli[DIMENSION]={qi[A_X]-ql[A_X],	qi[A_Y]-ql[A_Y],	qi[A_Z]-ql[A_Z]};
			double dis=sqrt(qli[A_X]*qli[A_X]+qli[A_Y]*qli[A_Y]+qli[A_Z]*qli[A_Z]);
			if(dis<r)
			{
				if(qli[A_X]<0)
				{
					DwDq[A_X]+=r/(ql[A_X]*ql[A_X]);
					HYPER1[l*h_num+i].DgDq[A_X]=-1/n0*r/(qi[A_X]*qi[A_X]);
				}
				if(qli[A_X]>0)
				{
					DwDq[A_X]+=-r/(ql[A_X]*ql[A_X]);
					HYPER1[l*h_num+i].DgDq[A_X]=1/n0*r/(qi[A_X]*qi[A_X]);
				}
				if(qli[A_Y]<0)
				{
					DwDq[A_Y]+=r/(ql[A_Y]*ql[A_Y]);
					HYPER1[l*h_num+i].DgDq[A_Y]=-1/n0*r/(qi[A_Y]*qi[A_Y]);
				}
				if(qli[A_Y]>0)
				{
					DwDq[A_Y]+=-r/(ql[A_Y]*ql[A_Y]);
					HYPER1[l*h_num+i].DgDq[A_Y]=1/n0*r/(qi[A_Y]*qi[A_Y]);
				}
				if(qli[A_Z]<0)
				{
					DwDq[A_Z]+=r/(ql[A_Z]*ql[A_Z]);
					HYPER1[l*h_num+i].DgDq[A_Z]=-1/n0*r/(qi[A_Z]*qi[A_Z]);
				}
				if(qli[A_Z]>0)
				{
					DwDq[A_Z]+=r/(ql[A_Z]*ql[A_Z]);
					HYPER1[l*h_num+i].DgDq[A_Z]=-1/n0*r/(qi[A_Z]*qi[A_Z]);
				}
			}
		}

		HYPER1[l*h_num+l].DgDq[A_X]+=1/n0*DwDq[A_X];
		HYPER1[l*h_num+l].DgDq[A_Y]+=1/n0*DwDq[A_Y];
		HYPER1[l*h_num+l].DgDq[A_Z]+=1/n0*DwDq[A_Z];
	}
	cout<<"----------OK"<<endl;
 }


/////�j���[�g�����t�\���@ 
void newton_raphson(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int rigid_number,int flag_w,int t)
{
	/////fx(N*N�s��̊e����)��1�����z��Ŋi�[�A(i,j)�����Ȃ�[j*N+i]�ŎQ��
	/////DfDx(N*N�s��̊e����)��1�����z��Ŋi�[�A(i,j)�����Ȃ�[j*N+i]�ŎQ��

	int calc_type=1;//�j���[�g���t�\���̔������@ 0:�Δ������̋t�s������̂܂܋��߂�@1:���`�������𗘗p

	//pn=2;//test,�Ƃ肠����2���łƂ��邩�ǂ����m�F 
	//////////////////�@f1(x1,x2) = x1^2 + x2^2 -5 = 0 f2(x1,x2) = x1^2/9+ x2^2 -1 = 0  http://homepage1.nifty.com/gfk/excel_newton_ren.htm

	int h_num=HYPER.size();
	int r_num=rigid_number;
	double *fx=new double [h_num];//�֐��l�B
	double *DfDx=new double [h_num*h_num];//�֐��̕Δ����l�B
	double *fx_Nw0=new double [h_num-r_num];
	double *DfDx_Nw0=new double [(h_num-r_num)*(h_num-r_num)];

	double *XX=new double [h_num];//���݂̉��B	
	double *XX_old=new double [h_num];//1�X�e�b�v�O�̉��B
	double ep=1e-5;//��������
	double E=1;//���݂̌덷
/*	double start=0;
	double end=0;
	double newton_t=0;*/
	int count=0;//������
	double d;
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	double Dt=CON.get_interval();
	double sum=0;
	double E_old=0;
	int dec_flag=ON;

	for(int i=0;i<h_num-r_num;i++)
	{
		XX[i]=0;
		XX_old[i]=0;
		fx[i]=0;
		fx_Nw0[i]=0;
		for(int j=0;j<h_num-r_num;j++)
		{
			DfDx[i*h_num+j]=0;
			DfDx_Nw0[i*(h_num-r_num)+j]=0;
		}
		for(int j=h_num-r_num;j<h_num;j++)	DfDx[i*h_num+j]=0;
	}
	for(int i=h_num-r_num;i<h_num;i++)
	{
		XX[i]=1;
		XX_old[i]=0;
		fx[i]=0;
		for(int j=0;j<h_num;j++)	DfDx[i*h_num+j]=0;
	}


	//	for(int i=0; i<N; i++) XX[i]=1;///�����l��^����B�Ƃ肠����1��
	cout<<"NR�@�J�n";
//	start=clock();
	while(E>ep)
	{
		count++;
		for(int i=0; i<h_num; i++)	XX_old[i]=XX[i];	//�����L��

//		if(count==1)	for(int i=0;i<h_num;i++)	for(int j=0;j<h_num;j++)	for(int D=0;D<DIMENSION;D++)	HYPER1[i*h_num+j].newton_DgDq[D]=HYPER1[i*N+j].DgDq[D];

		calc_newton_function(CON,PART,HYPER,HYPER1,XX,fx,DfDx,r_num,count,flag_w,t);


/*		//���݂̊֐��l�����߂�
		if(count==1) cout<<fx[0]<<" "<<fx[1]<<endl;
		//���݂̕Δ����l�����߂�
		//calc_DfDx(XX)////���݂̕Δ����l�����߂�B���e���̂Ȃ�΁Acalc_DgDq()�Ȃǂŋ��߂���͂�
		DfDx[0*N+0]=2*XX[0];
		DfDx[0*N+1]=2*XX[1];
		DfDx[1*N+0]=2*XX[0]/9;
		DfDx[1*N+1]=2*XX[1];
		if(count==1) cout<<DfDx[0]<<" "<<DfDx[1]<<" "<<DfDx[2]<<" "<<DfDx[3]<<endl;*/

		///�l�̍X�V
		if(calc_type==0)//�t�s��𗘗p �t�s�񂪋��܂肳������Α����͂�
		{
			calc_inverse_matrix_for_NR(h_num,DfDx);

			for(int i=0; i<h_num; i++) 
			{
				d=0; //�ω���
				for(int j=0; j<h_num; j++)	d+=DfDx[i*h_num+j]*fx[j];
				XX[i]-=d;
			}
		}
		else if(calc_type==1)//�t�s���p���Ȃ��A���肷��͂������A�x���Ȃ�͂�
		{
/*			if(flag_w==OFF)
			{
				for(int i=0;i<h_num-r_num;i++)
				{
					fx_Nw0[i]=fx[i];
					for(int j=0;j<h_num-r_num;j++)	DfDx_Nw0[i*(h_num-r_num)+j]=DfDx[i*h_num+j];
				}
				gauss(DfDx_Nw0,fx_Nw0,h_num-r_num);				
				for(int i=0;i<h_num-r_num;i++)	XX[i]-=fx_Nw0[i];//0.5*mi/(Dt*Dt)*V*fx[i];
			}
			else*/
			{
				gauss(DfDx,fx,h_num);
				for(int i=0;i<h_num;i++)	XX[i]-=fx[i];//0.5*mi/(Dt*Dt)*V*fx[i];
			}
		}

		//�덷�̕]��
		E_old=E;
		E=0;
		sum=0;
		for(int i=0; i<h_num; i++)
		{
			E+=fabs(XX[i]-XX_old[i]);
			sum+=fabs(XX[i]);
		}
		E/=sum;


		if(count==1 || count%200==0)
		{		
			/*
			cout<<"XX_old["<<i<<"]-d=X["<<i<<"]	";
			cout<<XX_old[i]<<" - ";
			cout<<d<<" = ";
			cout<<XX[i]<<endl;*/

			cout<<"������	"<<count<<" E="<<E<<endl;
//			output_newton_data2(E,XX,h_num,count,t);

		}
		if(count>CON.get_nr())	break;
		else if(dec_flag==ON)	if(E_old-E<0)	break;	
	}
//	end=clock();
//	newton_t=(end-start)/CLOCKS_PER_SEC;

	cout<<"��������";

	for(int i=0;i<h_num;i++)
	{
		HYPER[i].lambda=XX[i];
		//cout<<"lambda"<<i<<"	"<<HYPER[i].lambda<<endl;
	}

//	for(int i=0;i<N;i++)	cout<<"lambda["<<i<<"]="<<HYPER[i].lambda<<endl;
	delete[]	fx;
	delete[]	DfDx;
	delete[]	XX;
	delete[]	XX_old;
	delete[]	DfDx_Nw0;
	delete[]	fx_Nw0;		


	cout<<"---------- OK"<<endl;
}


void calc_newton_function(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> HYPER,vector<hyperelastic2> HYPER1,double *lambda,double *fx,double *DfDx,int rigid_number,int count,int flag_w,int t)
{
	flag_w=OFF;
	int h_num=HYPER.size();
	int r_num=rigid_number;
	double Dt=CON.get_dt();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	double mr=V*CON.get_rigid_density();
	double r=CON.get_h_dis();

	double *n_rx=new double[h_num];
	double *n_ry=new double[h_num];
	double *n_rz=new double[h_num];

	double **n_DgDq_x=new double *[h_num];
	double **n_DgDq_y=new double *[h_num];
	double **n_DgDq_z=new double *[h_num];

	for(int i=0;i<h_num;i++)
	{
		n_rx[i]=0;
		n_ry[i]=0;
		n_rz[i]=0;
		n_DgDq_x[i]=new double [h_num];
		n_DgDq_y[i]=new double [h_num];
		n_DgDq_z[i]=new double [h_num];
	}
	for(int i=0;i<h_num;i++)
	{
		for(int j=0;j<h_num;j++)
		{
			n_DgDq_x[i][j]=0;
			n_DgDq_y[i][j]=0;
			n_DgDq_z[i][j]=0;
		}
	}

	double **p_Fi=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	p_Fi[D]=new double [DIMENSION];

	////�ʒu���W�̍X�V	
	for(int i=0;i<h_num-r_num;i++)
	{
		double p_half_p[DIMENSION]={0,0,0};

		//�������q�e��
		int Ni=HYPER[i].N0;
		for(int j=0;j<Ni;j++)
		{	
			int jn=HYPER[i].NEI0[j];

			p_half_p[A_X]+=(HYPER[jn].stress[0][0]-lambda[jn])*HYPER1[jn*h_num+i].DgDq[0]+HYPER[jn].stress[0][1]*HYPER1[jn*h_num+i].DgDq[1]+HYPER[jn].stress[0][2]*HYPER1[jn*h_num+i].DgDq[2];
			p_half_p[A_Y]+=HYPER[jn].stress[1][0]*HYPER1[jn*h_num+i].DgDq[0]+(HYPER[jn].stress[1][1]-lambda[jn])*HYPER1[jn*h_num+i].DgDq[1]+HYPER[jn].stress[1][2]*HYPER1[jn*h_num+i].DgDq[2];
			p_half_p[A_Z]+=HYPER[jn].stress[2][0]*HYPER1[jn*h_num+i].DgDq[0]+HYPER[jn].stress[2][1]*HYPER1[jn*h_num+i].DgDq[1]+(HYPER[jn].stress[2][2]-lambda[jn])*HYPER1[jn*h_num+i].DgDq[2];
		}
		p_half_p[A_X]+=(HYPER[i].stress[0][0]-lambda[i])*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[0][1]*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[0][2]*HYPER1[i*h_num+i].DgDq[2];
		p_half_p[A_Y]+=HYPER[i].stress[1][0]*HYPER1[i*h_num+i].DgDq[0]+(HYPER[i].stress[1][1]-lambda[i])*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[1][2]*HYPER1[i*h_num+i].DgDq[2];
		p_half_p[A_Z]+=HYPER[i].stress[2][0]*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[2][1]*HYPER1[i*h_num+i].DgDq[1]+(HYPER[i].stress[2][2]-lambda[i])*HYPER1[i*h_num+i].DgDq[2];
		
		//�O�����q�e��
		for(int l=h_num-r_num;l<h_num;l++)
		{
			p_half_p[A_X]-=lambda[l]*HYPER1[l*h_num+i].DgDq[A_X];
			p_half_p[A_Y]-=lambda[l]*HYPER1[l*h_num+i].DgDq[A_Y];
			p_half_p[A_Z]-=lambda[l]*HYPER1[l*h_num+i].DgDq[A_Z];
		}

		//�d�͉e��
		p_half_p[A_Z]-=mi*9.8;

		//�ʒu���W�̌v�Z
		n_rx[i]=PART[i].r[A_X]+Dt*(HYPER[i].p[A_X]+Dt*0.5*p_half_p[A_X])/mi;
		n_ry[i]=PART[i].r[A_Y]+Dt*(HYPER[i].p[A_Y]+Dt*0.5*p_half_p[A_Y])/mi;
		n_rz[i]=PART[i].r[A_Z]+Dt*(HYPER[i].p[A_Z]+Dt*0.5*p_half_p[A_Z])/mi;
	}

	for(int k=h_num-r_num;k<h_num;k++)
	{

		double p_half_p[DIMENSION]={0,0,0};
		//�������q�e��
		for(int j=0;j<h_num;j++)
		{
			p_half_p[A_X]-=lambda[j]*HYPER1[j*h_num+k].DgDq[A_X];
			p_half_p[A_Y]-=lambda[j]*HYPER1[j*h_num+k].DgDq[A_Y];
			p_half_p[A_Z]-=lambda[j]*HYPER1[j*h_num+k].DgDq[A_Z];
		}
		//�O�����q�e��
		for(int l=h_num-r_num;l<h_num;l++)
		{
			p_half_p[A_X]-=lambda[l]*HYPER1[l*h_num+k].DgDq[A_X];
			p_half_p[A_Y]-=lambda[l]*HYPER1[l*h_num+k].DgDq[A_Y];
			p_half_p[A_Z]-=lambda[l]*HYPER1[l*h_num+k].DgDq[A_Z];
		}
		//�d�͉e��
		//p_half_p[A_Z]-=mr*9.8;
		n_rx[k]=PART[k].r[A_X]+Dt*(HYPER[k].p[A_X]+Dt*0.5*p_half_p[A_X])/mr;
		n_ry[k]=PART[k].r[A_Y]+Dt*(HYPER[k].p[A_Y]+Dt*0.5*p_half_p[A_Y])/mr;
		n_rz[k]=PART[k].r[A_Z]+Dt*(HYPER[k].p[A_Z]+Dt*0.5*p_half_p[A_Z])/mr;
	}


	////DgDq��fx�̍X�V
	for(int i=0;i<h_num-r_num;i++)
	{
		//Fi�̌v�Z
		int Ni=HYPER[i].N0;
		double fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};	

		for(int in=0;in<Ni;in++)
		{
			int inn=HYPER[i].NEI0[in];
			double w=HYPER1[i*h_num+inn].wiin;

			fi[0][0]+=w*(n_rx[inn]-n_rx[i])*HYPER1[i*h_num+inn].aiin[A_X];
			fi[0][1]+=w*(n_rx[inn]-n_rx[i])*HYPER1[i*h_num+inn].aiin[A_Y];
			fi[0][2]+=w*(n_rx[inn]-n_rx[i])*HYPER1[i*h_num+inn].aiin[A_Z];
			fi[1][0]+=w*(n_ry[inn]-n_ry[i])*HYPER1[i*h_num+inn].aiin[A_X];
			fi[1][1]+=w*(n_ry[inn]-n_ry[i])*HYPER1[i*h_num+inn].aiin[A_Y];
			fi[1][2]+=w*(n_ry[inn]-n_ry[i])*HYPER1[i*h_num+inn].aiin[A_Z];
			fi[2][0]+=w*(n_rz[inn]-n_rz[i])*HYPER1[i*h_num+inn].aiin[A_X];
			fi[2][1]+=w*(n_rz[inn]-n_rz[i])*HYPER1[i*h_num+inn].aiin[A_Y];
			fi[2][2]+=w*(n_rz[inn]-n_rz[i])*HYPER1[i*h_num+inn].aiin[A_Z];
		}

		p_Fi[0][0]=fi[0][0]*HYPER[i].inverse_Ai[0][0]+fi[0][1]*HYPER[i].inverse_Ai[1][0]+fi[0][2]*HYPER[i].inverse_Ai[2][0];
		p_Fi[0][1]=fi[0][0]*HYPER[i].inverse_Ai[0][1]+fi[0][1]*HYPER[i].inverse_Ai[1][1]+fi[0][2]*HYPER[i].inverse_Ai[2][1];
		p_Fi[0][2]=fi[0][0]*HYPER[i].inverse_Ai[0][2]+fi[0][1]*HYPER[i].inverse_Ai[1][2]+fi[0][2]*HYPER[i].inverse_Ai[2][2];
		p_Fi[1][0]=fi[1][0]*HYPER[i].inverse_Ai[0][0]+fi[1][1]*HYPER[i].inverse_Ai[1][0]+fi[1][2]*HYPER[i].inverse_Ai[2][0];	
		p_Fi[1][1]=fi[1][0]*HYPER[i].inverse_Ai[0][1]+fi[1][1]*HYPER[i].inverse_Ai[1][1]+fi[1][2]*HYPER[i].inverse_Ai[2][1];
		p_Fi[1][2]=fi[1][0]*HYPER[i].inverse_Ai[0][2]+fi[1][1]*HYPER[i].inverse_Ai[1][2]+fi[1][2]*HYPER[i].inverse_Ai[2][2];
		p_Fi[2][0]=fi[2][0]*HYPER[i].inverse_Ai[0][0]+fi[2][1]*HYPER[i].inverse_Ai[1][0]+fi[2][2]*HYPER[i].inverse_Ai[2][0];
		p_Fi[2][1]=fi[2][0]*HYPER[i].inverse_Ai[0][1]+fi[2][1]*HYPER[i].inverse_Ai[1][1]+fi[2][2]*HYPER[i].inverse_Ai[2][1];
		p_Fi[2][2]=fi[2][0]*HYPER[i].inverse_Ai[0][2]+fi[2][1]*HYPER[i].inverse_Ai[1][2]+fi[2][2]*HYPER[i].inverse_Ai[2][2];

		//J�̌v�Z
		double J=calc_det3(p_Fi);

		//fx�̌v�Z
		fx[i]=V*(1-J);//1-J;//
		//t_inverse_Fi�̌v�Z
		inverse(p_Fi,DIMENSION);

		//DgDq�̌v�Z
		for(int j=0;j<Ni;j++)
		{
			int k=HYPER[i].NEI0[j];
			n_DgDq_x[i][k]=J*(p_Fi[0][0]*HYPER1[k*h_num+i].n0ij[0]+p_Fi[1][0]*HYPER1[k*h_num+i].n0ij[1]+p_Fi[2][0]*HYPER1[k*h_num+i].n0ij[2]);
			n_DgDq_y[i][k]=J*(p_Fi[0][1]*HYPER1[k*h_num+i].n0ij[0]+p_Fi[1][1]*HYPER1[k*h_num+i].n0ij[1]+p_Fi[2][1]*HYPER1[k*h_num+i].n0ij[2]);
			n_DgDq_z[i][k]=J*(p_Fi[0][2]*HYPER1[k*h_num+i].n0ij[0]+p_Fi[1][2]*HYPER1[k*h_num+i].n0ij[1]+p_Fi[2][2]*HYPER1[k*h_num+i].n0ij[2]);
		}
		n_DgDq_x[i][i]=J*(p_Fi[0][0]*HYPER1[i*h_num+i].n0ij[0]+p_Fi[1][0]*HYPER1[i*h_num+i].n0ij[1]+p_Fi[2][0]*HYPER1[i*h_num+i].n0ij[2]);
		n_DgDq_y[i][i]=J*(p_Fi[0][1]*HYPER1[i*h_num+i].n0ij[0]+p_Fi[1][1]*HYPER1[i*h_num+i].n0ij[1]+p_Fi[2][1]*HYPER1[i*h_num+i].n0ij[2]);
		n_DgDq_z[i][i]=J*(p_Fi[0][2]*HYPER1[i*h_num+i].n0ij[0]+p_Fi[1][2]*HYPER1[i*h_num+i].n0ij[1]+p_Fi[2][2]*HYPER1[i*h_num+i].n0ij[2]);
		
		double DwDq[DIMENSION]={0,0,0};		
		double n0=HYPER[i].pnd;
		for(int k=h_num-r_num;k<h_num;k++)
		{			
			double qik[DIMENSION]={n_rx[k]-n_rx[i],	n_ry[k]-n_ry[i],	n_rz[k]-n_rz[i]};
			double dis=sqrt(qik[A_X]*qik[A_X]+qik[A_Y]*qik[A_Y]+qik[A_Z]*qik[A_Z]);
			if(dis<r)
			{
				flag_w=ON;
				if(qik[A_X]<0)
				{
					DwDq[A_X]+=r/(n_rx[i]*n_rx[i]);
					n_DgDq_x[i][k]=-1/n0*r/(n_rx[k]*n_rx[k]);
				}
				if(qik[A_X]>0)
				{
					DwDq[A_X]+=-r/(n_rx[i]*n_rx[i]);
					n_DgDq_x[i][k]=1/n0*r/(n_rx[k]*n_rx[k]);
				}
				if(qik[A_Y]<0)
				{
					DwDq[A_Y]+=r/(n_ry[i]*n_ry[i]);
					n_DgDq_y[i][k]=-1/n0*r/(n_ry[k]*n_ry[k]);
				}
				if(qik[A_Y]>0)
				{
					DwDq[A_Y]+=-r/(n_ry[i]*n_ry[i]);
					n_DgDq_y[i][k]=1/n0*r/(n_ry[k]*n_ry[k]);
				}
				if(qik[A_Z]<0)
				{
					DwDq[A_Z]+=r/(n_rz[i]*n_rz[i]);
					n_DgDq_z[i][k]=-1/n0*r/(n_rz[k]*n_rz[k]);
				}
				if(qik[A_Z]>0)
				{
					DwDq[A_Z]+=-r/(n_rz[i]*n_rz[i]);					
					n_DgDq_z[i][k]=1/n0*r/(n_rz[k]*n_rz[k]);
				}
			}
		}
		n_DgDq_x[i][i]+=1/n0*DwDq[A_X];
		n_DgDq_y[i][i]+=1/n0*DwDq[A_Y];
		n_DgDq_z[i][i]+=1/n0*DwDq[A_Z];
	}

//	if(flag_w==ON)
	{
		for(int k=h_num-r_num;k<h_num;k++)
		{
			//Fi�̌v�Z
			int Nk=HYPER[k].N0;
			double fk[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};	

			for(int kn=0;kn<Nk;kn++)
			{
				int knn=HYPER[k].NEI0[kn];
				double w=HYPER1[k*h_num+knn].wiin;

				fk[0][0]+=w*(n_rx[knn]-n_rx[k])*HYPER1[k*h_num+knn].aiin[A_X];
				fk[0][1]+=w*(n_rx[knn]-n_rx[k])*HYPER1[k*h_num+knn].aiin[A_Y];
				fk[0][2]+=w*(n_rx[knn]-n_rx[k])*HYPER1[k*h_num+knn].aiin[A_Z];
				fk[1][0]+=w*(n_ry[knn]-n_ry[k])*HYPER1[k*h_num+knn].aiin[A_X];
				fk[1][1]+=w*(n_ry[knn]-n_ry[k])*HYPER1[k*h_num+knn].aiin[A_Y];
				fk[1][2]+=w*(n_ry[knn]-n_ry[k])*HYPER1[k*h_num+knn].aiin[A_Z];
				fk[2][0]+=w*(n_rz[knn]-n_rz[k])*HYPER1[k*h_num+knn].aiin[A_X];
				fk[2][1]+=w*(n_rz[knn]-n_rz[k])*HYPER1[k*h_num+knn].aiin[A_Y];
				fk[2][2]+=w*(n_rz[knn]-n_rz[k])*HYPER1[k*h_num+knn].aiin[A_Z];
			}

			p_Fi[0][0]=fk[0][0]*HYPER[k].inverse_Ai[0][0]+fk[0][1]*HYPER[k].inverse_Ai[1][0]+fk[0][2]*HYPER[k].inverse_Ai[2][0];
			p_Fi[0][1]=fk[0][0]*HYPER[k].inverse_Ai[0][1]+fk[0][1]*HYPER[k].inverse_Ai[1][1]+fk[0][2]*HYPER[k].inverse_Ai[2][1];
			p_Fi[0][2]=fk[0][0]*HYPER[k].inverse_Ai[0][2]+fk[0][1]*HYPER[k].inverse_Ai[1][2]+fk[0][2]*HYPER[k].inverse_Ai[2][2];
			p_Fi[1][0]=fk[1][0]*HYPER[k].inverse_Ai[0][0]+fk[1][1]*HYPER[k].inverse_Ai[1][0]+fk[1][2]*HYPER[k].inverse_Ai[2][0];	
			p_Fi[1][1]=fk[1][0]*HYPER[k].inverse_Ai[0][1]+fk[1][1]*HYPER[k].inverse_Ai[1][1]+fk[1][2]*HYPER[k].inverse_Ai[2][1];
			p_Fi[1][2]=fk[1][0]*HYPER[k].inverse_Ai[0][2]+fk[1][1]*HYPER[k].inverse_Ai[1][2]+fk[1][2]*HYPER[k].inverse_Ai[2][2];
			p_Fi[2][0]=fk[2][0]*HYPER[k].inverse_Ai[0][0]+fk[2][1]*HYPER[k].inverse_Ai[1][0]+fk[2][2]*HYPER[k].inverse_Ai[2][0];
			p_Fi[2][1]=fk[2][0]*HYPER[k].inverse_Ai[0][1]+fk[2][1]*HYPER[k].inverse_Ai[1][1]+fk[2][2]*HYPER[k].inverse_Ai[2][1];
			p_Fi[2][2]=fk[2][0]*HYPER[k].inverse_Ai[0][2]+fk[2][1]*HYPER[k].inverse_Ai[1][2]+fk[2][2]*HYPER[k].inverse_Ai[2][2];
		
			//J�̌v�Z
			double J=calc_det3(p_Fi);

			//fx�̌v�Z
			fx[k]=V*(1-J);//1-J;//

			//t_inverse_Fi�̌v�Z
			inverse(p_Fi,DIMENSION);


			//DgDq�̌v�Z
			for(int kn=0;kn<Nk;kn++)
			{
				int knn=HYPER[k].NEI0[kn];
				n_DgDq_x[k][knn]=J*(p_Fi[0][0]*HYPER1[knn*h_num+k].n0ij[0]+p_Fi[1][0]*HYPER1[knn*h_num+k].n0ij[1]+p_Fi[2][0]*HYPER1[knn*h_num+k].n0ij[2]);
				n_DgDq_y[k][knn]=J*(p_Fi[0][1]*HYPER1[knn*h_num+k].n0ij[0]+p_Fi[1][1]*HYPER1[knn*h_num+k].n0ij[1]+p_Fi[2][1]*HYPER1[knn*h_num+k].n0ij[2]);
				n_DgDq_z[k][knn]=J*(p_Fi[0][2]*HYPER1[knn*h_num+k].n0ij[0]+p_Fi[1][2]*HYPER1[knn*h_num+k].n0ij[1]+p_Fi[2][2]*HYPER1[knn*h_num+k].n0ij[2]);
			}
			n_DgDq_x[k][k]=J*(p_Fi[0][0]*HYPER1[k*h_num+k].n0ij[0]+p_Fi[1][0]*HYPER1[k*h_num+k].n0ij[1]+p_Fi[2][0]*HYPER1[k*h_num+k].n0ij[2]);
			n_DgDq_y[k][k]=J*(p_Fi[0][1]*HYPER1[k*h_num+k].n0ij[0]+p_Fi[1][1]*HYPER1[k*h_num+k].n0ij[1]+p_Fi[2][1]*HYPER1[k*h_num+k].n0ij[2]);
			n_DgDq_z[k][k]=J*(p_Fi[0][2]*HYPER1[k*h_num+k].n0ij[0]+p_Fi[1][2]*HYPER1[k*h_num+k].n0ij[1]+p_Fi[2][2]*HYPER1[k*h_num+k].n0ij[2]);

			double n0=HYPER[k].pnd;
			double DwDq[DIMENSION]={0,0,0};
			for(int i=0;i<h_num-r_num;i++)
			{
				double qki[DIMENSION]={n_rx[i]-n_rx[k],	n_ry[i]-n_ry[k],	n_rz[i]-n_rz[k]};
				double dis=sqrt(qki[A_X]*qki[A_X]+qki[A_Y]*qki[A_Y]+qki[A_Z]*qki[A_Z]);
				if(dis<r)
				{
					if(qki[A_X]<0)
					{
						DwDq[A_X]+=r/(n_rx[k]*n_rx[k]);
						n_DgDq_x[k][i]=-1/n0*r/(n_rx[i]*n_rx[i]);
					}
					if(qki[A_X]>0)
					{
						DwDq[A_X]+=-r/(n_rx[k]*n_rx[k]);
						n_DgDq_x[k][i]=1/n0*r/(n_rx[i]*n_rx[i]);
					}
					if(qki[A_Y]<0)
					{
						DwDq[A_Y]+=r/(n_ry[k]*n_ry[k]);
						n_DgDq_y[k][i]=-1/n0*r/(n_ry[i]*n_ry[i]);
					}
					if(qki[A_Y]>0)
					{
						DwDq[A_Y]+=-r/(n_ry[k]*n_ry[k]);
						n_DgDq_y[k][i]=1/n0*r/(n_ry[i]*n_ry[i]);
					}
					if(qki[A_Z]<0)
					{
						DwDq[A_Z]+=r/(n_rz[k]*n_rz[k]);
						n_DgDq_z[k][i]=-1/n0*r/(n_rz[i]*n_rz[i]);
					}
					if(qki[A_Z]>0)
					{
						DwDq[A_Z]+=-r/(n_rz[k]*n_rz[k]);
						n_DgDq_z[k][i]=1/n0*r/(n_rz[i]*n_rz[i]);
					}
				}
			}
			n_DgDq_x[k][k]+=1/n0*DwDq[A_X];
			n_DgDq_y[k][k]+=1/n0*DwDq[A_Y];
			n_DgDq_z[k][k]+=1/n0*DwDq[A_Z];
		}
	}
	
	////DfDx�̍X�V
/*	if(flag_w==OFF)
	{
		for(int i=0;i<h_num-r_num;i++)
		{
			for(int j=0;j<h_num-r_num;j++)
			{
				double DFDlam=0;
				for(int t=0;t<h_num-r_num;t++)	DFDlam+=n_DgDq_x[i][t]*HYPER1[j*h_num+t].DgDq[A_X]+n_DgDq_y[i][t]*HYPER1[j*h_num+t].DgDq[A_Y]+n_DgDq_z[i][t]*HYPER1[j*h_num+t].DgDq[A_Z];
				DfDx[i*h_num+j]=-Dt*Dt*0.5*1/mi*DFDlam;//-DFDlambda;//	
			}
		}			
	}
	else*/
	{
		for(int i=0;i<h_num-r_num;i++)
		{
			for(int j=0;j<h_num-r_num;j++)
			{
				double DFDlam_hhh=0;
				for(int t=0;t<h_num-r_num;t++)	DFDlam_hhh+=n_DgDq_x[i][t]*HYPER1[j*h_num+t].DgDq[A_X]+n_DgDq_y[i][t]*HYPER1[j*h_num+t].DgDq[A_Y]+n_DgDq_z[i][t]*HYPER1[j*h_num+t].DgDq[A_Z];
				double DFDlam_hhr=0;
				for(int s=h_num-r_num;s<h_num;s++)	DFDlam_hhr+=n_DgDq_x[i][s]*HYPER1[j*h_num+s].DgDq[A_X]+n_DgDq_y[i][s]*HYPER1[j*h_num+s].DgDq[A_Y]+n_DgDq_z[i][s]*HYPER1[j*h_num+s].DgDq[A_Z];
				DfDx[i*h_num+j]=-Dt*Dt*0.5*(1/mi*DFDlam_hhh+1/mr*DFDlam_hhr);//-DFDlambda;//			
			}
			for(int l=h_num-r_num;l<h_num;l++)
			{
				double DFDlam_hrh=0;
				for(int t=0;t<h_num-r_num;t++)	DFDlam_hrh+=n_DgDq_x[i][t]*HYPER1[l*h_num+t].DgDq[A_X]+n_DgDq_y[i][t]*HYPER1[l*h_num+t].DgDq[A_Y]+n_DgDq_z[i][t]*HYPER1[l*h_num+t].DgDq[A_Z];
				double DFDlam_hrr=0;
				for(int s=h_num-r_num;s<h_num;s++)	DFDlam_hrr+=n_DgDq_x[i][s]*HYPER1[l*h_num+s].DgDq[A_X]+n_DgDq_y[i][s]*HYPER1[l*h_num+s].DgDq[A_Y]+n_DgDq_z[i][s]*HYPER1[l*h_num+s].DgDq[A_Z];
				DfDx[i*h_num+l]=-Dt*Dt*0.5*(1/mi*DFDlam_hrh+1/mr*DFDlam_hrr);//-DFDlambda;//			
			}
		}
	
		for(int k=h_num-r_num;k<h_num;k++)
		{
			for(int j=0;j<h_num-r_num;j++)
			{
				double DFDlam_rhh=0;
				for(int t=0;t<h_num-r_num;t++)	DFDlam_rhh+=n_DgDq_x[k][t]*HYPER1[j*h_num+t].DgDq[A_X]+n_DgDq_y[k][t]*HYPER1[j*h_num+t].DgDq[A_Y]+n_DgDq_z[k][t]*HYPER1[j*h_num+t].DgDq[A_Z];
				double DFDlam_rhr=0;
				for(int s=h_num-r_num;s<h_num;s++)	DFDlam_rhr+=n_DgDq_x[k][s]*HYPER1[j*h_num+s].DgDq[A_X]+n_DgDq_y[k][s]*HYPER1[j*h_num+s].DgDq[A_Y]+n_DgDq_z[k][s]*HYPER1[j*h_num+s].DgDq[A_Z];			
				DfDx[k*h_num+j]=-Dt*Dt*0.5*(1/mi*DFDlam_rhh+1/mr*DFDlam_rhr);//-DFDlambda;//	
			}
			for(int l=h_num-r_num;l<h_num;l++)
			{
				double DFDlam_rrh=0;
				for(int t=0;t<h_num-r_num;t++)	DFDlam_rrh+=n_DgDq_x[k][t]*HYPER1[l*h_num+t].DgDq[A_X]+n_DgDq_y[k][t]*HYPER1[l*h_num+t].DgDq[A_Y]+n_DgDq_z[k][t]*HYPER1[l*h_num+t].DgDq[A_Z];
				double DFDlam_rrr=0;
				for(int s=h_num-r_num;s<h_num;s++)	DFDlam_rrr+=n_DgDq_x[k][s]*HYPER1[l*h_num+s].DgDq[A_X]+n_DgDq_y[k][s]*HYPER1[l*h_num+s].DgDq[A_Y]+n_DgDq_z[k][s]*HYPER1[l*h_num+s].DgDq[A_Z];
				DfDx[k*h_num+l]=-Dt*Dt*0.5*(1/mi*DFDlam_rrh+1/mr*DFDlam_rrr);//-DFDlambda;//			
			}
		}
	}
	

	////�o��
//	if(count%200==0 && count>CON.get_nr()/2)
	if(t==1||t%CON.get_interval()==0)	if(count%200==0||count==1)	output_newton_data1(fx,DfDx,n_rx,n_ry,n_rz,h_num,count,t);

	for(int D=0;D<DIMENSION;D++)	delete[]	p_Fi[D];
	delete[]	p_Fi;

	for(int i=0;i<h_num;i++)
	{
		delete[]	n_DgDq_x[i];
		delete[]	n_DgDq_y[i];
		delete[]	n_DgDq_z[i];
	}
	delete[]	n_DgDq_x;
	delete[]	n_DgDq_y;
	delete[]	n_DgDq_z;

	delete[]	n_rx;
	delete[]	n_ry;
	delete[]	n_rz;
}

void calc_half_p(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,bool repetation,int rigid_number,int flag_w,int t)
{
	if(repetation==0)	cout<<"���̉^���ʁ��ʒu���W�v�Z";
	else	cout<<"�^���ʌv�Z";

	int h_num=HYPER.size();
	int r_num=rigid_number;
	double Dt=CON.get_dt();
	double le=CON.get_distancebp();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	double mr=V*CON.get_rigid_density();

	for(int i=0;i<h_num;i++)
	{
		double p_half_p[DIMENSION]={0,0,0};
		int Ni=HYPER[i].N0;
		//�������q�̉e��
		for(int j=0;j<Ni;j++)
		{		
			int jn=HYPER[i].NEI0[j];
			p_half_p[A_X]+=(HYPER[jn].stress[0][0]-HYPER[jn].lambda)*HYPER1[jn*h_num+i].DgDq[0]+HYPER[jn].stress[0][1]*HYPER1[jn*h_num+i].DgDq[1]+HYPER[jn].stress[0][2]*HYPER1[jn*h_num+i].DgDq[2];
			p_half_p[A_Y]+=HYPER[jn].stress[1][0]*HYPER1[jn*h_num+i].DgDq[0]+(HYPER[jn].stress[1][1]-HYPER[jn].lambda)*HYPER1[jn*h_num+i].DgDq[1]+HYPER[jn].stress[1][2]*HYPER1[jn*h_num+i].DgDq[2];
			p_half_p[A_Z]+=HYPER[jn].stress[2][0]*HYPER1[jn*h_num+i].DgDq[0]+HYPER[jn].stress[2][1]*HYPER1[jn*h_num+i].DgDq[1]+(HYPER[jn].stress[2][2]-HYPER[jn].lambda)*HYPER1[jn*h_num+i].DgDq[2];
		}
		p_half_p[A_X]+=(HYPER[i].stress[0][0]-HYPER[i].lambda)*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[0][1]*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[0][2]*HYPER1[i*h_num+i].DgDq[2];
		p_half_p[A_Y]+=HYPER[i].stress[1][0]*HYPER1[i*h_num+i].DgDq[0]+(HYPER[i].stress[1][1]-HYPER[i].lambda)*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[1][2]*HYPER1[i*h_num+i].DgDq[2];
		p_half_p[A_Z]+=HYPER[i].stress[2][0]*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[2][1]*HYPER1[i*h_num+i].DgDq[1]+(HYPER[i].stress[2][2]-HYPER[i].lambda)*HYPER1[i*h_num+i].DgDq[2];
		
		//�O�����q�̉e��
	//	if(flag_w==ON)
		{
			for(int l=h_num-r_num;l<h_num;l++)
			{
				p_half_p[A_X]-=HYPER[l].lambda*HYPER1[l*h_num+i].DgDq[A_X];
				p_half_p[A_Y]-=HYPER[l].lambda*HYPER1[l*h_num+i].DgDq[A_Y];
				p_half_p[A_Z]-=HYPER[l].lambda*HYPER1[l*h_num+i].DgDq[A_Z];
			}
		}
		
		//�d�͂̉e��
		p_half_p[A_Z]-=mi*9.8;
		
		//�ʒu���W�̌v�Z
		if(repetation==0)
		{
			HYPER[i].half_p[A_X]=HYPER[i].p[A_X]+Dt*0.5*p_half_p[A_X];
			HYPER[i].half_p[A_Y]=HYPER[i].p[A_Y]+Dt*0.5*p_half_p[A_Y];
			HYPER[i].half_p[A_Z]=HYPER[i].p[A_Z]+Dt*0.5*p_half_p[A_Z];
			
			PART[i].r[A_X]+=Dt*(HYPER[i].p[A_X]+Dt*0.5*p_half_p[A_X])/mi;
			PART[i].r[A_Y]+=Dt*(HYPER[i].p[A_Y]+Dt*0.5*p_half_p[A_Y])/mi;
			PART[i].r[A_Z]+=Dt*(HYPER[i].p[A_Z]+Dt*0.5*p_half_p[A_Z])/mi;
		}
		else
		{
			//�^���ʂ̍X�V
			HYPER[i].p[A_X]=HYPER[i].half_p[A_X]+Dt*0.5*p_half_p[A_X];
			HYPER[i].p[A_Y]=HYPER[i].half_p[A_Y]+Dt*0.5*p_half_p[A_Y];
			HYPER[i].p[A_Z]=HYPER[i].half_p[A_Z]+Dt*0.5*p_half_p[A_Z];
			//���x�̍X�V
			PART[i].u[A_X]=HYPER[i].half_p[A_X]/mi;
			PART[i].u[A_Y]=HYPER[i].half_p[A_Y]/mi;
			PART[i].u[A_Z]=HYPER[i].half_p[A_Z]/mi;
			//�p�^���ʂ̍X�V
			HYPER[i].ang_p[A_X]=PART[i].r[A_Y]*HYPER[i].p[A_Z]-PART[i].r[A_Z]*HYPER[i].p[A_Y];
			HYPER[i].ang_p[A_Y]=PART[i].r[A_Z]*HYPER[i].p[A_X]-PART[i].r[A_X]*HYPER[i].p[A_Z];
			HYPER[i].ang_p[A_Z]=PART[i].r[A_X]*HYPER[i].p[A_Y]-PART[i].r[A_Y]*HYPER[i].p[A_X];
		}
	}

//	if(flag_w==ON)
	{

		for(int k=h_num-r_num;k<h_num;k++)
		{
			double p_half_p[DIMENSION]={0,0,0};

			//�������q�̉e��
			for(int j=0;j<h_num-r_num;j++)
			{	
				p_half_p[A_X]-=HYPER[j].lambda*HYPER1[j*h_num+k].DgDq[A_X];
				p_half_p[A_Y]-=HYPER[j].lambda*HYPER1[j*h_num+k].DgDq[A_Y];
				p_half_p[A_Z]-=HYPER[j].lambda*HYPER1[j*h_num+k].DgDq[A_Z];
			}

			//�O�����q�̉e��
			for(int l=h_num-r_num;l<h_num;l++)
			{
				p_half_p[A_X]-=HYPER[l].lambda*HYPER1[l*h_num+k].DgDq[A_X];
				p_half_p[A_Y]-=HYPER[l].lambda*HYPER1[l*h_num+k].DgDq[A_Y];
				p_half_p[A_Z]-=HYPER[l].lambda*HYPER1[l*h_num+k].DgDq[A_Z];
			}

			//�d�͂̉e��
			//p_half_p[A_Z]-=9.8*mr;

			if(repetation==0)
			{
				HYPER[k].half_p[A_X]=HYPER[k].p[A_X]+Dt*0.5*p_half_p[A_X];
				HYPER[k].half_p[A_Y]=HYPER[k].p[A_Y]+Dt*0.5*p_half_p[A_Y];
				HYPER[k].half_p[A_Z]=HYPER[k].p[A_Z]+Dt*0.5*p_half_p[A_Z];
			
				PART[k].r[A_X]+=Dt*HYPER[k].half_p[A_X]/mr;
				PART[k].r[A_Y]+=Dt*HYPER[k].half_p[A_Y]/mr;
				PART[k].r[A_Z]+=Dt*HYPER[k].half_p[A_Z]/mr;
			}
			else
			{
				//�^���ʂ̍X�V
				HYPER[k].p[A_X]=HYPER[k].half_p[A_X]+Dt*0.5*p_half_p[A_X];
				HYPER[k].p[A_Y]=HYPER[k].half_p[A_Y]+Dt*0.5*p_half_p[A_Y];
				HYPER[k].p[A_Z]=HYPER[k].half_p[A_Z]+Dt*0.5*p_half_p[A_Z];
				//���x�̍X�V
				PART[k].u[A_X]=HYPER[k].half_p[A_X]/mr;
				PART[k].u[A_Y]=HYPER[k].half_p[A_Y]/mr;
				PART[k].u[A_Z]=HYPER[k].half_p[A_Z]/mr;
				//�p�^���ʂ̍X�V
				HYPER[k].ang_p[A_X]=PART[k].r[A_Y]*HYPER[k].p[A_Z]-PART[k].r[A_Z]*HYPER[k].p[A_Y];
				HYPER[k].ang_p[A_Y]=PART[k].r[A_Z]*HYPER[k].p[A_X]-PART[k].r[A_X]*HYPER[k].p[A_Z];
				HYPER[k].ang_p[A_Z]=PART[k].r[A_X]*HYPER[k].p[A_Y]-PART[k].r[A_Y]*HYPER[k].p[A_X];
			}
		}
	}
	
/*	if(repetation==0)
	{
		for(int i=0;i<h_num;i++)
		{
			cout<<"half_p["<<i<<"]=";
			for(int D=0;D<DIMENSION;D++)	cout<<HYPER[i].half_p[D]<<" ";
			cout<<endl;
		}
		cout<<endl;
		*/
	cout<<"----------OK"<<endl;
}

void calc_F(mpsconfig &CON, vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int rigid_number,int flag_w,int t)
{
	cout<<"Fi�v�Z";
	////Fi�̍X�V
	int h_num=HYPER.size();
	int r_num=rigid_number;
	double r=CON.get_h_dis();

	double **p_Fi=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	p_Fi[D]=new double[DIMENSION];

	for(int i=0;i<h_num-r_num;i++)
	{
		double fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
		//Fi�̌v�Z

		int Ni=HYPER[i].N0;	
		for(int in=0;in<Ni;in++)
		{
			int inn=HYPER[i].NEI0[in];
			double w=HYPER1[i*h_num+inn].wiin;
			double a[DIMENSION]={HYPER1[i*h_num+inn].aiin[A_X],	HYPER1[i*h_num+inn].aiin[A_Y],	HYPER1[i*h_num+inn].aiin[A_Z]};
			
			fi[0][0]+=w*(PART[inn].r[A_X]-PART[i].r[A_X])*a[A_X];	fi[0][1]+=w*(PART[inn].r[A_X]-PART[i].r[A_X])*a[A_Y];	fi[0][2]+=w*(PART[inn].r[A_X]-PART[i].r[A_X])*a[A_Z];
			fi[1][0]+=w*(PART[inn].r[A_Y]-PART[i].r[A_Y])*a[A_X];	fi[1][1]+=w*(PART[inn].r[A_Y]-PART[i].r[A_Y])*a[A_Y];	fi[1][2]+=w*(PART[inn].r[A_Y]-PART[i].r[A_Y])*a[A_Z];
			fi[2][0]+=w*(PART[inn].r[A_Z]-PART[i].r[A_Z])*a[A_X];	fi[2][1]+=w*(PART[inn].r[A_Z]-PART[i].r[A_Z])*a[A_Y];	fi[2][2]+=w*(PART[inn].r[A_Z]-PART[i].r[A_Z])*a[A_Z];
		}

		p_Fi[0][0]=fi[0][0]*HYPER[i].inverse_Ai[0][0]+fi[0][1]*HYPER[i].inverse_Ai[1][0]+fi[0][2]*HYPER[i].inverse_Ai[2][0];
		p_Fi[0][1]=fi[0][0]*HYPER[i].inverse_Ai[0][1]+fi[0][1]*HYPER[i].inverse_Ai[1][1]+fi[0][2]*HYPER[i].inverse_Ai[2][1];
		p_Fi[0][2]=fi[0][0]*HYPER[i].inverse_Ai[0][2]+fi[0][1]*HYPER[i].inverse_Ai[1][2]+fi[0][2]*HYPER[i].inverse_Ai[2][2];
		p_Fi[1][0]=fi[1][0]*HYPER[i].inverse_Ai[0][0]+fi[1][1]*HYPER[i].inverse_Ai[1][0]+fi[1][2]*HYPER[i].inverse_Ai[2][0];
		p_Fi[1][1]=fi[1][0]*HYPER[i].inverse_Ai[0][1]+fi[1][1]*HYPER[i].inverse_Ai[1][1]+fi[1][2]*HYPER[i].inverse_Ai[2][1];
		p_Fi[1][2]=fi[1][0]*HYPER[i].inverse_Ai[0][2]+fi[1][1]*HYPER[i].inverse_Ai[1][2]+fi[1][2]*HYPER[i].inverse_Ai[2][2];
		p_Fi[2][0]=fi[2][0]*HYPER[i].inverse_Ai[0][0]+fi[2][1]*HYPER[i].inverse_Ai[1][0]+fi[2][2]*HYPER[i].inverse_Ai[2][0];
		p_Fi[2][1]=fi[2][0]*HYPER[i].inverse_Ai[0][1]+fi[2][1]*HYPER[i].inverse_Ai[1][1]+fi[2][2]*HYPER[i].inverse_Ai[2][1];
		p_Fi[2][2]=fi[2][0]*HYPER[i].inverse_Ai[0][2]+fi[2][1]*HYPER[i].inverse_Ai[1][2]+fi[2][2]*HYPER[i].inverse_Ai[2][2];
		
		HYPER[i].Fi[0][0]=p_Fi[0][0];	HYPER[i].Fi[0][1]=p_Fi[0][1];	HYPER[i].Fi[0][2]=p_Fi[0][2];	
		HYPER[i].Fi[1][0]=p_Fi[1][0];	HYPER[i].Fi[1][1]=p_Fi[1][1];	HYPER[i].Fi[1][2]=p_Fi[1][2];	
		HYPER[i].Fi[2][0]=p_Fi[2][0];	HYPER[i].Fi[2][1]=p_Fi[2][1];	HYPER[i].Fi[2][2]=p_Fi[2][2];	
		//J�̌v�Z
		double J=calc_det3(p_Fi);
	//	for(int i=0;i<h_num;i++)	cout<<"J["<<i<<"]="<<J<<endl;
		HYPER[i].J=J;

		//t_inverse_Fi�̌v�Z
		inverse(p_Fi,DIMENSION);
		HYPER[i].t_inverse_Fi[0][0]=p_Fi[0][0];	HYPER[i].t_inverse_Fi[0][1]=p_Fi[1][0];	HYPER[i].t_inverse_Fi[0][2]=p_Fi[2][0];
		HYPER[i].t_inverse_Fi[1][0]=p_Fi[0][1];	HYPER[i].t_inverse_Fi[1][1]=p_Fi[1][1];	HYPER[i].t_inverse_Fi[1][2]=p_Fi[2][1];
		HYPER[i].t_inverse_Fi[2][0]=p_Fi[0][2];	HYPER[i].t_inverse_Fi[2][1]=p_Fi[1][2];	HYPER[i].t_inverse_Fi[2][2]=p_Fi[2][2];
	}
	for(int k=h_num-r_num;k<h_num;k++)
	{
		double fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
		//Fi�̌v�Z

		int Nk=HYPER[k].N0;	
		for(int kn=0;kn<Nk;kn++)
		{
			int knn=HYPER[k].NEI0[kn];
			double w=HYPER1[k*h_num+knn].wiin;
			double a[DIMENSION]={HYPER1[k*h_num+knn].aiin[A_X],	HYPER1[k*h_num+knn].aiin[A_Y],	HYPER1[k*h_num+knn].aiin[A_Z]};
			
			fi[0][0]+=w*(PART[knn].r[A_X]-PART[k].r[A_X])*a[A_X];	fi[0][1]+=w*(PART[knn].r[A_X]-PART[k].r[A_X])*a[A_Y];	fi[0][2]+=w*(PART[knn].r[A_X]-PART[k].r[A_X])*a[A_Z];
			fi[1][0]+=w*(PART[knn].r[A_Y]-PART[k].r[A_Y])*a[A_X];	fi[1][1]+=w*(PART[knn].r[A_Y]-PART[k].r[A_Y])*a[A_Y];	fi[1][2]+=w*(PART[knn].r[A_Y]-PART[k].r[A_Y])*a[A_Z];
			fi[2][0]+=w*(PART[knn].r[A_Z]-PART[k].r[A_Z])*a[A_X];	fi[2][1]+=w*(PART[knn].r[A_Z]-PART[k].r[A_Z])*a[A_Y];	fi[2][2]+=w*(PART[knn].r[A_Z]-PART[k].r[A_Z])*a[A_Z];
		}

		p_Fi[0][0]=fi[0][0]*HYPER[k].inverse_Ai[0][0]+fi[0][1]*HYPER[k].inverse_Ai[1][0]+fi[0][2]*HYPER[k].inverse_Ai[2][0];
		p_Fi[0][1]=fi[0][0]*HYPER[k].inverse_Ai[0][1]+fi[0][1]*HYPER[k].inverse_Ai[1][1]+fi[0][2]*HYPER[k].inverse_Ai[2][1];
		p_Fi[0][2]=fi[0][0]*HYPER[k].inverse_Ai[0][2]+fi[0][1]*HYPER[k].inverse_Ai[1][2]+fi[0][2]*HYPER[k].inverse_Ai[2][2];
		p_Fi[1][0]=fi[1][0]*HYPER[k].inverse_Ai[0][0]+fi[1][1]*HYPER[k].inverse_Ai[1][0]+fi[1][2]*HYPER[k].inverse_Ai[2][0];
		p_Fi[1][1]=fi[1][0]*HYPER[k].inverse_Ai[0][1]+fi[1][1]*HYPER[k].inverse_Ai[1][1]+fi[1][2]*HYPER[k].inverse_Ai[2][1];
		p_Fi[1][2]=fi[1][0]*HYPER[k].inverse_Ai[0][2]+fi[1][1]*HYPER[k].inverse_Ai[1][2]+fi[1][2]*HYPER[k].inverse_Ai[2][2];
		p_Fi[2][0]=fi[2][0]*HYPER[k].inverse_Ai[0][0]+fi[2][1]*HYPER[k].inverse_Ai[1][0]+fi[2][2]*HYPER[k].inverse_Ai[2][0];
		p_Fi[2][1]=fi[2][0]*HYPER[k].inverse_Ai[0][1]+fi[2][1]*HYPER[k].inverse_Ai[1][1]+fi[2][2]*HYPER[k].inverse_Ai[2][1];
		p_Fi[2][2]=fi[2][0]*HYPER[k].inverse_Ai[0][2]+fi[2][1]*HYPER[k].inverse_Ai[1][2]+fi[2][2]*HYPER[k].inverse_Ai[2][2];
		
		HYPER[k].Fi[0][0]=p_Fi[0][0];	HYPER[k].Fi[0][1]=p_Fi[0][1];	HYPER[k].Fi[0][2]=p_Fi[0][2];	
		HYPER[k].Fi[1][0]=p_Fi[1][0];	HYPER[k].Fi[1][1]=p_Fi[1][1];	HYPER[k].Fi[1][2]=p_Fi[1][2];	
		HYPER[k].Fi[2][0]=p_Fi[2][0];	HYPER[k].Fi[2][1]=p_Fi[2][1];	HYPER[k].Fi[2][2]=p_Fi[2][2];	
		//J�̌v�Z
		double J=calc_det3(p_Fi);
	//	for(int i=0;i<h_num;i++)	cout<<"J["<<i<<"]="<<J<<endl;
		HYPER[k].J=J;

		//t_inverse_Fi�̌v�Z
		inverse(p_Fi,DIMENSION);
		HYPER[k].t_inverse_Fi[0][0]=p_Fi[0][0];	HYPER[k].t_inverse_Fi[0][1]=p_Fi[1][0];	HYPER[k].t_inverse_Fi[0][2]=p_Fi[2][0];
		HYPER[k].t_inverse_Fi[1][0]=p_Fi[0][1];	HYPER[k].t_inverse_Fi[1][1]=p_Fi[1][1];	HYPER[k].t_inverse_Fi[1][2]=p_Fi[2][1];
		HYPER[k].t_inverse_Fi[2][0]=p_Fi[0][2];	HYPER[k].t_inverse_Fi[2][1]=p_Fi[1][2];	HYPER[k].t_inverse_Fi[2][2]=p_Fi[2][2];
	}
	
	for(int D=0;D<DIMENSION;D++)	delete[]	p_Fi[D];
	delete[]	p_Fi;

	//calculation of DgDq
	for(int i=0;i<h_num-r_num;i++)
	{
		int Ni=HYPER[i].N0;
		for(int j=0;j<Ni;j++)
		{
			int k=HYPER[i].NEI0[j];
			HYPER1[k*h_num+i].DgDq[A_X]=HYPER[k].J*(HYPER[k].t_inverse_Fi[A_X][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].t_inverse_Fi[A_X][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].t_inverse_Fi[A_X][2]*HYPER1[i*h_num+k].n0ij[2]);
			HYPER1[k*h_num+i].DgDq[A_Y]=HYPER[k].J*(HYPER[k].t_inverse_Fi[A_Y][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].t_inverse_Fi[A_Y][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].t_inverse_Fi[A_Y][2]*HYPER1[i*h_num+k].n0ij[2]);
			HYPER1[k*h_num+i].DgDq[A_Z]=HYPER[k].J*(HYPER[k].t_inverse_Fi[A_Z][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].t_inverse_Fi[A_Z][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].t_inverse_Fi[A_Z][2]*HYPER1[i*h_num+k].n0ij[2]);		
		}
		HYPER1[i*h_num+i].DgDq[A_X]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_X][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_X][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_X][2]*HYPER1[i*h_num+i].n0ij[2]);
		HYPER1[i*h_num+i].DgDq[A_Y]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_Y][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_Y][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_Y][2]*HYPER1[i*h_num+i].n0ij[2]);
		HYPER1[i*h_num+i].DgDq[A_Z]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_Z][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_Z][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_Z][2]*HYPER1[i*h_num+i].n0ij[2]);
		
	//	if(flag_w==ON)
		{

			double DwDq[DIMENSION]={0,0,0};		
			double n0=HYPER[i].pnd;
			for(int k=h_num-r_num;k<h_num;k++)
			{
				double qk[DIMENSION]={PART[k].r[A_X],PART[k].r[A_Y],PART[k].r[A_Z]};
				double qi[DIMENSION]={PART[i].r[A_X],PART[i].r[A_Y],PART[i].r[A_Z]};
				double qik[DIMENSION]={qk[A_X]-qi[A_X],	qk[A_Y]-qi[A_Y],	qk[A_Z]-qi[A_Z]};
				double dis=sqrt(qik[A_X]*qik[A_X]+qik[A_Y]*qik[A_Y]+qik[A_Z]*qik[A_Z]);
				if(dis<r)
				{
					if(qik[A_X]<0)
					{
						DwDq[A_X]+=r/(qi[A_X]*qi[A_X]);
						HYPER1[i*h_num+k].DgDq[A_X]=-1/n0*r/(qk[A_X]*qk[A_X]);
					}
					if(qik[A_X]>0)
					{
						DwDq[A_X]+=-r/(qi[A_X]*qi[A_X]);
						HYPER1[i*h_num+k].DgDq[A_X]=1/n0*r/(qk[A_X]*qk[A_X]);
					}
					if(qik[A_Y]<0)
					{
						DwDq[A_Y]+=r/(qi[A_Y]*qi[A_Y]);
						HYPER1[i*h_num+k].DgDq[A_Y]=-1/n0*r/(qk[A_Y]*qk[A_Y]);
					}
					if(qik[A_Y]>0)
					{
						DwDq[A_Y]+=-r/(qi[A_Y]*qi[A_Y]);
						HYPER1[i*h_num+k].DgDq[A_Y]=1/n0*r/(qk[A_Y]*qk[A_Y]);
					}
					if(qik[A_Z]<0)
					{
						DwDq[A_Z]+=r/(qi[A_Z]*qi[A_Z]);
						HYPER1[i*h_num+k].DgDq[A_Z]=-1/n0*r/(qk[A_Z]*qk[A_Z]);
					}
					if(qik[A_Z]>0)
					{
						DwDq[A_Z]+=-r/(qi[A_Z]*qi[A_Z]);					
						HYPER1[i*h_num+k].DgDq[A_Z]=1/n0*r/(qk[A_Z]*qk[A_Z]);
					}
				}
			}
			HYPER1[i*h_num+i].DgDq[A_X]+=1/n0*DwDq[A_X];
			HYPER1[i*h_num+i].DgDq[A_X]+=1/n0*DwDq[A_Y];
			HYPER1[i*h_num+i].DgDq[A_X]+=1/n0*DwDq[A_Z];
		}

		for(int l=h_num-r_num;l<h_num;l++)
		{
			int Nl=HYPER[l].N0;
			for(int ln=0;ln<Nl;ln++)
			{
				int lnn=HYPER[l].NEI0[ln];
				HYPER1[lnn*h_num+l].DgDq[A_X]=HYPER[lnn].J*(HYPER[lnn].t_inverse_Fi[A_X][0]*HYPER1[l*h_num+lnn].n0ij[0]+HYPER[lnn].t_inverse_Fi[A_X][1]*HYPER1[l*h_num+lnn].n0ij[1]+HYPER[lnn].t_inverse_Fi[A_X][2]*HYPER1[l*h_num+lnn].n0ij[2]);
				HYPER1[lnn*h_num+l].DgDq[A_Y]=HYPER[lnn].J*(HYPER[lnn].t_inverse_Fi[A_Y][0]*HYPER1[l*h_num+lnn].n0ij[0]+HYPER[lnn].t_inverse_Fi[A_Y][1]*HYPER1[l*h_num+lnn].n0ij[1]+HYPER[lnn].t_inverse_Fi[A_Y][2]*HYPER1[l*h_num+lnn].n0ij[2]);
				HYPER1[lnn*h_num+l].DgDq[A_Z]=HYPER[lnn].J*(HYPER[lnn].t_inverse_Fi[A_Z][0]*HYPER1[l*h_num+lnn].n0ij[0]+HYPER[lnn].t_inverse_Fi[A_Z][1]*HYPER1[l*h_num+lnn].n0ij[1]+HYPER[lnn].t_inverse_Fi[A_Z][2]*HYPER1[l*h_num+lnn].n0ij[2]);		
			}
			HYPER1[l*h_num+l].DgDq[A_X]=HYPER[l].J*(HYPER[l].t_inverse_Fi[A_X][0]*HYPER1[l*h_num+l].n0ij[0]+HYPER[l].t_inverse_Fi[A_X][1]*HYPER1[l*h_num+l].n0ij[1]+HYPER[l].t_inverse_Fi[A_X][2]*HYPER1[l*h_num+l].n0ij[2]);
			HYPER1[l*h_num+l].DgDq[A_Y]=HYPER[l].J*(HYPER[l].t_inverse_Fi[A_Y][0]*HYPER1[l*h_num+l].n0ij[0]+HYPER[l].t_inverse_Fi[A_Y][1]*HYPER1[l*h_num+l].n0ij[1]+HYPER[l].t_inverse_Fi[A_Y][2]*HYPER1[l*h_num+l].n0ij[2]);
			HYPER1[l*h_num+l].DgDq[A_Z]=HYPER[l].J*(HYPER[l].t_inverse_Fi[A_Z][0]*HYPER1[l*h_num+l].n0ij[0]+HYPER[l].t_inverse_Fi[A_Z][1]*HYPER1[l*h_num+l].n0ij[1]+HYPER[l].t_inverse_Fi[A_Z][2]*HYPER1[l*h_num+l].n0ij[2]);
		
			double n0=HYPER[l].pnd;
			double DwDq[DIMENSION]={0,0,0};
			for(int i=0;i<h_num-r_num;i++)
			{
				double qi[DIMENSION]={PART[i].r[A_X],PART[i].r[A_Y],PART[i].r[A_Z]};
				double ql[DIMENSION]={PART[l].r[A_X],PART[l].r[A_Y],PART[l].r[A_Z]};
				double qli[DIMENSION]={qi[A_X]-ql[A_X],	qi[A_Y]-ql[A_Y],	qi[A_Z]-ql[A_Z]};
				double dis=sqrt(qli[A_X]*qli[A_X]+qli[A_Y]*qli[A_Y]+qli[A_Z]*qli[A_Z]);
				if(dis<r)
				{
					if(qli[A_X]<0)
					{
						DwDq[A_X]+=r/(ql[A_X]*ql[A_X]);
						HYPER1[l*h_num+i].DgDq[A_X]=-1/n0*r/(qi[A_X]*qi[A_X]);
					}
					if(qli[A_X]>0)
					{
						DwDq[A_X]+=-r/(ql[A_X]*ql[A_X]);
						HYPER1[l*h_num+i].DgDq[A_X]=1/n0*r/(qi[A_X]*qi[A_X]);
					}
					if(qli[A_Y]<0)
					{
						DwDq[A_Y]+=r/(ql[A_Y]*ql[A_Y]);
						HYPER1[l*h_num+i].DgDq[A_Y]=-1/n0*r/(qi[A_Y]*qi[A_Y]);
					}
					if(qli[A_Y]>0)
					{
						DwDq[A_Y]+=-r/(ql[A_Y]*ql[A_Y]);
						HYPER1[l*h_num+i].DgDq[A_Y]=1/n0*r/(qi[A_Y]*qi[A_Y]);
					}
					if(qli[A_Z]<0)
					{
						DwDq[A_Z]+=r/(ql[A_Z]*ql[A_Z]);
						HYPER1[l*h_num+i].DgDq[A_Z]=-1/n0*r/(qi[A_Z]*qi[A_Z]);
					}
					if(qli[A_Z]>0)
					{
						DwDq[A_Z]+=-r/(ql[A_Z]*ql[A_Z]);
						HYPER1[l*h_num+i].DgDq[A_Z]=1/n0*r/(qi[A_Z]*qi[A_Z]);
					}
				}
			}
			HYPER1[l*h_num+l].DgDq[A_X]+=1/n0*DwDq[A_X];
			HYPER1[l*h_num+l].DgDq[A_Y]+=1/n0*DwDq[A_Y];
			HYPER1[l*h_num+l].DgDq[A_Z]+=1/n0*DwDq[A_Z];
		}
	}

	cout<<"----------OK"<<endl;

}

void calc_stress(mpsconfig &CON,vector<hyperelastic> &HYPER,int rigid_number)
{
	cout<<"���͌v�Z";
	int d_num=HYPER.size()-rigid_number;

	double **d_Fi=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	d_Fi[D]=new double [DIMENSION];

	double c10=CON.get_c10();
	double c01=CON.get_c01();
	
	double b[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double bb[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

	for(int j=0;j<d_num;j++)
	{	
		double J=HYPER[j].J;
		if(J<0){
			d_Fi[0][0]=-1/pow(-J,1/3)*HYPER[j].Fi[0][0];	d_Fi[0][1]=-1/pow(-J,1/3)*HYPER[j].Fi[0][1];	d_Fi[0][2]=-1/pow(-J,1/3)*HYPER[j].Fi[0][2];
			d_Fi[1][0]=-1/pow(-J,1/3)*HYPER[j].Fi[1][0];	d_Fi[1][1]=-1/pow(-J,1/3)*HYPER[j].Fi[1][1];	d_Fi[1][2]=-1/pow(-J,1/3)*HYPER[j].Fi[1][2];
			d_Fi[2][0]=-1/pow(-J,1/3)*HYPER[j].Fi[2][0];	d_Fi[2][1]=-1/pow(-J,1/3)*HYPER[j].Fi[2][1];	d_Fi[2][2]=-1/pow(-J,1/3)*HYPER[j].Fi[2][2];
		}
		else
		{
			d_Fi[0][0]=1/pow(J,1/3)*HYPER[j].Fi[0][0];	d_Fi[0][1]=1/pow(J,1/3)*HYPER[j].Fi[0][1];	d_Fi[0][2]=1/pow(J,1/3)*HYPER[j].Fi[0][2];
			d_Fi[1][0]=1/pow(J,1/3)*HYPER[j].Fi[1][0];	d_Fi[1][1]=1/pow(J,1/3)*HYPER[j].Fi[1][1];	d_Fi[1][2]=1/pow(J,1/3)*HYPER[j].Fi[1][2];
			d_Fi[2][0]=1/pow(J,1/3)*HYPER[j].Fi[2][0];	d_Fi[2][1]=1/pow(J,1/3)*HYPER[j].Fi[2][1];	d_Fi[2][2]=1/pow(J,1/3)*HYPER[j].Fi[2][2];
		}

		b[0][0]=d_Fi[0][0]*d_Fi[0][0]+d_Fi[0][1]*d_Fi[0][1]+d_Fi[0][2]*d_Fi[0][2];
		b[0][1]=d_Fi[0][0]*d_Fi[1][0]+d_Fi[0][1]*d_Fi[1][1]+d_Fi[0][2]*d_Fi[1][2];
		b[0][2]=d_Fi[0][0]*d_Fi[2][0]+d_Fi[0][1]*d_Fi[2][1]+d_Fi[0][2]*d_Fi[2][2];
		b[1][0]=d_Fi[1][0]*d_Fi[0][0]+d_Fi[1][1]*d_Fi[0][1]+d_Fi[1][2]*d_Fi[0][2];
		b[1][1]=d_Fi[1][0]*d_Fi[1][0]+d_Fi[1][1]*d_Fi[1][1]+d_Fi[1][2]*d_Fi[1][2];
		b[1][2]=d_Fi[1][0]*d_Fi[2][0]+d_Fi[1][1]*d_Fi[2][1]+d_Fi[1][2]*d_Fi[2][2];
		b[2][0]=d_Fi[2][0]*d_Fi[0][0]+d_Fi[2][1]*d_Fi[0][1]+d_Fi[2][2]*d_Fi[0][2];
		b[2][1]=d_Fi[2][0]*d_Fi[1][0]+d_Fi[2][1]*d_Fi[1][1]+d_Fi[2][2]*d_Fi[1][2];
		b[2][2]=d_Fi[2][0]*d_Fi[2][0]+d_Fi[2][1]*d_Fi[2][1]+d_Fi[2][2]*d_Fi[2][2];

		bb[0][0]=b[0][0]*b[0][0]+b[0][1]*b[1][0]+b[0][2]*b[2][0];
		bb[0][1]=b[0][0]*b[0][1]+b[0][1]*b[1][1]+b[0][2]*b[2][1];
		bb[0][2]=b[0][0]*b[0][2]+b[0][1]*b[1][2]+b[0][2]*b[2][2];
		bb[1][0]=b[1][0]*b[0][0]+b[1][1]*b[1][0]+b[1][2]*b[2][0];
		bb[1][1]=b[1][0]*b[0][1]+b[1][1]*b[1][1]+b[1][2]*b[2][1];
		bb[1][2]=b[1][0]*b[0][2]+b[1][1]*b[1][2]+b[1][2]*b[2][2];
		bb[2][0]=b[2][0]*b[0][0]+b[2][1]*b[1][0]+b[2][2]*b[2][0];
		bb[2][1]=b[2][0]*b[0][1]+b[2][1]*b[1][1]+b[2][2]*b[2][1];
		bb[2][2]=b[2][0]*b[0][2]+b[2][1]*b[1][2]+b[2][2]*b[2][2];

		double trace_b=b[0][0]+b[1][1]+b[2][2];
		double trace_bb=bb[0][0]+bb[1][1]+bb[2][2];
	
		HYPER[j].stress[0][0]=2/J*((c10+c01*trace_b)*(b[0][0]-1.0/3.0*trace_b)-c01*(bb[0][0]-1.0/3.0*trace_bb));
		HYPER[j].stress[0][1]=2/J*((c10+c01*trace_b)*b[0][1]-c01*bb[0][1]);
		HYPER[j].stress[0][2]=2/J*((c10+c01*trace_b)*b[0][2]-c01*bb[0][2]);
		HYPER[j].stress[1][0]=2/J*((c10+c01*trace_b)*b[1][0]-c01*bb[1][0]);
		HYPER[j].stress[1][1]=2/J*((c10+c01*trace_b)*(b[1][1]-1.0/3.0*trace_b)-c01*(bb[1][1]-1.0/3.0*trace_bb));
		HYPER[j].stress[1][2]=2/J*((c10+c01*trace_b)*b[1][2]-c01*bb[1][2]);
		HYPER[j].stress[2][0]=2/J*((c10+c01*trace_b)*b[2][0]-c01*bb[2][0]);
		HYPER[j].stress[2][1]=2/J*((c10+c01*trace_b)*b[2][1]-c01*bb[2][1]);
		HYPER[j].stress[2][2]=2/J*((c10+c01*trace_b)*(b[2][2]-1.0/3.0*trace_b)-c01*(bb[2][2]-1.0/3.0*trace_bb));
		
		/*cout<<"stress="<<endl;
		cout<<HYPER[j].stress[0][0]<<"	"<<HYPER[j].stress[0][1]<<"	"<<HYPER[j].stress[0][2]<<endl;
		cout<<HYPER[j].stress[1][0]<<"	"<<HYPER[j].stress[1][1]<<"	"<<HYPER[j].stress[1][2]<<endl;
		cout<<HYPER[j].stress[2][0]<<"	"<<HYPER[j].stress[2][1]<<"	"<<HYPER[j].stress[2][2]<<endl;*/
	}

	for(int D=0;D<DIMENSION;D++)	delete[]	d_Fi[D];
	delete[]	d_Fi;

	cout<<"----------OK"<<endl;
}

void calc_differential_p(mpsconfig &CON,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int rigid_number)
{
	cout<<"�^���ʔ����l�v�Z";

	int h_num=HYPER.size();
	int r_num=rigid_number;
	double Dt=CON.get_dt();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	for(int i=0;i<h_num-r_num;i++)
	{
		double p_differential_p[DIMENSION]={0,0,0};
		int Ni=HYPER[i].N0;
		for(int j=0;j<Ni;j++)
		{				
			int k=HYPER[i].NEI0[j];
			p_differential_p[A_X]+=HYPER[k].stress[A_X][0]*HYPER1[k*h_num+i].DgDq[0]+HYPER[k].stress[A_X][1]*HYPER1[k*h_num+i].DgDq[1]+HYPER[k].stress[A_X][2]*HYPER1[k*h_num+i].DgDq[2];
			p_differential_p[A_Y]+=HYPER[k].stress[A_Y][0]*HYPER1[k*h_num+i].DgDq[0]+HYPER[k].stress[A_Y][1]*HYPER1[k*h_num+i].DgDq[1]+HYPER[k].stress[A_Y][2]*HYPER1[k*h_num+i].DgDq[2];	
			p_differential_p[A_Z]+=HYPER[k].stress[A_Z][0]*HYPER1[k*h_num+i].DgDq[0]+HYPER[k].stress[A_Z][1]*HYPER1[k*h_num+i].DgDq[1]+HYPER[k].stress[A_Z][2]*HYPER1[k*h_num+i].DgDq[2];	
		}
		HYPER[i].differential_p[A_X]=HYPER[i].half_p[A_X]+Dt*0.5*(p_differential_p[A_X]+HYPER[i].stress[A_X][0]*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[A_X][1]*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[A_X][2]*HYPER1[i*h_num+i].DgDq[2]);
		HYPER[i].differential_p[A_X]=HYPER[i].half_p[A_X]+Dt*0.5*(p_differential_p[A_X]+HYPER[i].stress[A_Y][0]*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[A_Y][1]*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[A_Y][2]*HYPER1[i*h_num+i].DgDq[2]);	
		HYPER[i].differential_p[A_Z]=HYPER[i].half_p[A_Z]+Dt*0.5*(p_differential_p[A_Z]+HYPER[i].stress[A_Z][0]*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[A_Z][1]*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[A_Z][2]*HYPER1[i*h_num+i].DgDq[2])-9.8*mi;//p_differential_p[A_Z];////
		
		//�d�͉e��
		HYPER[i].differential_p[A_Z]-=Dt*0.5*9.8*mi;
	}
	cout<<"----------OK"<<endl;

}

void renew_lambda(mpsconfig &CON,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int rigid_number,int flag_w,int t)
{
	cout<<"Lambda�v�Z";

	int h_num=HYPER.size();
	int r_num=rigid_number;
	double le=CON.get_distancebp();
	double Dt=CON.get_dt();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	double mr=V*CON.get_rigid_density();


	double *N_L=new double[h_num*h_num];
	double *N_R=new double[h_num];
	double *N_L_Nw0=new double[(h_num-r_num)*(h_num-r_num)];
	double *N_R_Nw0=new double[h_num-r_num];
	for(int i=0;i<h_num;i++)
	{
		N_R[i]=0;
		for(int j=0;j<h_num;j++)	N_L[j*h_num+i]=0;
	}

/*	if(flag_w==OFF)
	{
		for(int i=0;i<h_num-r_num;i++)
		{
			double N_Rh=0;
			for(int j=0;j<h_num-r_num;j++)
			{
				double N_Lh=0;
				for(int t=0;t<h_num-r_num;t++)
				{
					N_Lh+=1/mi*(HYPER1[i*h_num+t].DgDq[A_X]*HYPER1[j*h_num+t].DgDq[A_X]+HYPER1[i*h_num+t].DgDq[A_Y]*HYPER1[j*h_num+t].DgDq[A_Y]+HYPER1[i*h_num+t].DgDq[A_Z]*HYPER1[j*h_num+t].DgDq[A_Z]);
				}
				N_L_Nw0[i*(h_num-r_num)+j]=Dt*0.5*N_Lh;
				N_Rh+=1/mi*(HYPER1[i*h_num+j].DgDq[A_X]*HYPER[j].differential_p[A_X]+HYPER1[i*h_num+j].DgDq[A_Y]*HYPER[j].differential_p[A_Y]+HYPER1[i*h_num+j].DgDq[A_Z]*HYPER[j].differential_p[A_Z]);
			}
			N_R_Nw0[i]=N_Rh;
		}
		gauss(N_L_Nw0,N_R_Nw0,h_num-r_num);
	}
	else*/
	{
		for(int i=0;i<h_num-r_num;i++)
		{
			double N_Rh=0;
			for(int j=0;j<h_num-r_num;j++)
			{
				double N_Lh=0;
				for(int t=0;t<h_num-r_num;t++)
				{
					N_Lh+=1/mi*(HYPER1[i*h_num+t].DgDq[A_X]*HYPER1[j*h_num+t].DgDq[A_X]+HYPER1[i*h_num+t].DgDq[A_Y]*HYPER1[j*h_num+t].DgDq[A_Y]+HYPER1[i*h_num+t].DgDq[A_Z]*HYPER1[j*h_num+t].DgDq[A_Z]);
				}
				for(int s=h_num-r_num;s<h_num;s++)
				{
					N_Lh+=1/mr*(HYPER1[i*h_num+s].DgDq[A_X]*HYPER1[j*h_num+s].DgDq[A_X]+HYPER1[i*h_num+s].DgDq[A_Y]*HYPER1[j*h_num+s].DgDq[A_Y]+HYPER1[i*h_num+s].DgDq[A_Z]*HYPER1[j*h_num+s].DgDq[A_Z]);
				}
				N_L[i*h_num+j]=Dt*0.5*N_Lh;

				N_Rh+=1/mi*(HYPER1[i*h_num+j].DgDq[A_X]*HYPER[j].differential_p[A_X]+HYPER1[i*h_num+j].DgDq[A_Y]*HYPER[j].differential_p[A_Y]+HYPER1[i*h_num+j].DgDq[A_Z]*HYPER[j].differential_p[A_Z]);
			}
			for(int l=h_num-r_num;l<h_num;l++)
			{
				double N_Lr=0;
				for(int t=0;t<h_num-r_num;t++)
				{
					N_Lr+=1/mi*(HYPER1[i*h_num+t].DgDq[A_X]*HYPER1[l*h_num+t].DgDq[A_X]+HYPER1[i*h_num+t].DgDq[A_Y]*HYPER1[l*h_num+t].DgDq[A_Y]+HYPER1[i*h_num+t].DgDq[A_Z]*HYPER1[l*h_num+t].DgDq[A_Z]);
				}
				for(int s=h_num-r_num;s<h_num;s++)
				{
					N_Lr+=1/mr*(HYPER1[i*h_num+s].DgDq[A_X]*HYPER1[l*h_num+s].DgDq[A_X]+HYPER1[i*h_num+s].DgDq[A_Y]*HYPER1[l*h_num+s].DgDq[A_Y]+HYPER1[i*h_num+s].DgDq[A_Z]*HYPER1[l*h_num+s].DgDq[A_Z]);
				}
				N_L[i*h_num+l]=Dt*0.5*N_Lr;

				N_Rh+=1/mr*(HYPER1[i*h_num+l].DgDq[A_X]*HYPER[l].half_p[A_X]+HYPER1[i*h_num+l].DgDq[A_Y]*HYPER[l].half_p[A_Y]+HYPER1[i*h_num+l].DgDq[A_Z]*HYPER[l].half_p[A_Z]);
			}
			N_R[i]=N_Rh;
		}
	
		for(int k=h_num-r_num;k<h_num;k++)
		{
			double N_Rr=0;
			for(int j=0;j<h_num-r_num;j++)
			{
				double N_Lh=0;
				for(int t=0;t<h_num-r_num;t++)
				{
					N_Lh+=1/mi*(HYPER1[k*h_num+t].DgDq[A_X]*HYPER1[j*h_num+t].DgDq[A_X]+HYPER1[k*h_num+t].DgDq[A_Y]*HYPER1[j*h_num+t].DgDq[A_Y]+HYPER1[k*h_num+t].DgDq[A_Z]*HYPER1[j*h_num+t].DgDq[A_Z]);
				}
				for(int s=h_num-r_num;s<h_num;s++)
				{
					N_Lh+=1/mr*(HYPER1[k*h_num+s].DgDq[A_X]*HYPER1[j*h_num+s].DgDq[A_X]+HYPER1[k*h_num+s].DgDq[A_Y]*HYPER1[j*h_num+s].DgDq[A_Y]+HYPER1[k*h_num+s].DgDq[A_Z]*HYPER1[j*h_num+s].DgDq[A_Z]);
				}
				N_L[k*h_num+j]=Dt*0.5*N_Lh;

				N_Rr+=1/mi*(HYPER1[k*h_num+j].DgDq[A_X]*HYPER[j].differential_p[A_X]+HYPER1[k*h_num+j].DgDq[A_Y]*HYPER[j].differential_p[A_Y]+HYPER1[k*h_num+j].DgDq[A_Z]*HYPER[j].differential_p[A_Z]);
			}
			for(int l=h_num-r_num;l<h_num;l++)
			{
				double N_Lr=0;
				for(int t=0;t<h_num-r_num;t++)
				{
					N_Lr+=1/mi*(HYPER1[k*h_num+t].DgDq[A_X]*HYPER1[l*h_num+t].DgDq[A_X]+HYPER1[k*h_num+t].DgDq[A_Y]*HYPER1[l*h_num+t].DgDq[A_Y]+HYPER1[k*h_num+t].DgDq[A_Z]*HYPER1[l*h_num+t].DgDq[A_Z]);
				}
				for(int s=h_num-r_num;s<h_num;s++)
				{
					N_Lr+=1/mr*(HYPER1[k*h_num+s].DgDq[A_X]*HYPER1[l*h_num+s].DgDq[A_X]+HYPER1[k*h_num+s].DgDq[A_Y]*HYPER1[l*h_num+s].DgDq[A_Y]+HYPER1[k*h_num+s].DgDq[A_Z]*HYPER1[l*h_num+s].DgDq[A_Z]);
				}
				N_L[k*h_num+l]=Dt*0.5*N_Lr;

				N_Rr+=1/mr*(HYPER1[k*h_num+l].DgDq[A_X]*HYPER[l].half_p[A_X]+HYPER1[k*h_num+l].DgDq[A_Y]*HYPER[l].half_p[A_Y]+HYPER1[k*h_num+l].DgDq[A_Z]*HYPER[l].half_p[A_Z]);
			}
			N_R[k]=N_Rr;
		}
	

		stringstream ssl;
		stringstream ssr;
		ssl<<"./NL"<<t<<".csv";
		ssr<<"./NR"<<t<<".csv";
		ofstream nl("ssl.str()");
		ofstream nr("ssr.str()");

		for(int i=0;i<h_num;i++)
		{
			nr<<i<<","<<N_R[i]<<endl;
			for(int j=0;j<h_num;j++)	nl<<N_L[i*h_num+j]<<",";
			nl<<endl;
		}
		//lambda�����߂�
		gauss(N_L,N_R,h_num);
	}

	for(int i=0;i<h_num;i++)	HYPER[i].lambda=N_R[i];
	//for(int i=0;i<h_num;i++)	cout<<"lambda"<<i<<"	"<<HYPER[i].lambda<<endl;

	delete [] N_L;
	delete [] N_R;

	cout<<"----------OK"<<endl;
}

//det�����߂�֐��@�����ȗ��̂��ߎ��M��
double calc_det(double **M,int N)
{
	double det=0;
	double det_plus;
	double det_minus;
	int	a=0;

	for(int j=0;j<N;j++)
	{
		det_plus=1;
		for(int i=0;i<N;i++)	det_plus*=M[i%N][(j+i)%N];
		det+=det_plus;

		det_minus=1;
		for(int i=0;i<N;i++)
		{
			if(j-i<0)
			{
				a=j-i;
				a+=N;
				det_minus*=M[i%N][a%N];
			}
			else
			{
				det_minus*=M[i%N][(j-i)%N];
			}
		}
		det-=det_minus;
	}
	return det;
}

double calc_det3(double **M)
{
	double det=0;
	det+=M[0][0]*M[1][1]*M[2][2]+M[0][1]*M[1][2]*M[2][0]+M[0][2]*M[1][0]*M[2][1];
	det-=M[0][2]*M[1][1]*M[2][0]+M[0][1]*M[1][0]*M[2][2]+M[0][0]*M[1][2]*M[2][1];
	return det;
}


//�t�s������߂�֐� function.h�̈��͌v�Z�p�֐��ɗގ��������O�̊֐����������̂ł��݂킯
void calc_inverse_matrix_for_NR(int N, double *a)
{
	//N=���m��
	double buf=0;

	double *inv_a=new  double[N*N];					//�t�s��i�[
	
	//�P�ʍs������
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<N;j++)
		{
			if(i==j) inv_a[j*N+i]=1;
			else		inv_a[j*N+i]=0;
		}
	}

	//�|���o���@
	for(int i=0;i<N;i++)
	{
		if(a[i*N+i]<DBL_EPSILON) cout<<"a[i][i]=0"<<endl;
		buf=1/a[i*N+i];
		for(int j=0;j<N;j++)
		{
			a[j*N+i]*=buf;
			inv_a[j*N+i]*=buf;
		}
		for(int j=0;j<N;j++)
		{
			if(i!=j)
			{
				buf=a[i*N+j];//a[j][i];
				for(int k=0;k<N;k++)
				{
					a[k*N+j]-=a[k*N+i]*buf;//a[j][k]-=a[i][k]*buf;
					inv_a[k*N+j]-=inv_a[k*N+i]*buf;
				}
			}
		}
	}

	for(int i=0;i<N;i++)	for(int j=0;j<N;j++) a[i*N+j]=inv_a[i*N+j];

	delete [] inv_a;
}



void inverse(double **a,int N)
{
	double d=0;
	double *col=new double [N];
	double **y=new double *[N];
	for(int i=0;i<N;i++)	y[i]=new double [N];
	int *index=new int [N];

	//int p_revision=0,m_revision=0;
	int fig=0,x=0;
	int	flag=OFF;

	for(int D=0;D<N;D++)
	{
		col[D]=0;
		index[D]=0;
		for(int D2=0;D2<N;D2++)	y[D][D2]=0;
	}

	//for(int D=0;D<N;D++)	for(int D2=0;D2<N;D2++)	if(fabs(a[D][D2])<DBL_EPSILON)	a[D][D2]=0;

	for(int D=0;D<N;D++)
	{
		for(int D2=0;D2<N;D2++)
		{
			if(fabs(a[D][D2])>pow(10.0,5.0))
			{
				flag=ON;
				x=log10(fabs(a[D][D2]));
				if(x>fig)	fig=x;
			}
		}
	}	
	
	if(flag==ON)	for(int D=0;D<N;D++) for(int D2=0;D2<N;D2++) a[D][D2]/=pow(10.0,fig);

	ludcmp(a,N,index,&d);
	for(int j=0;j<N;j++)
	{
		for(int i=0;i<N;i++)	col[i]=0.0;
		col[j]=1.0;
		lubksb(a,N,index,col);
		for(int i=0;i<N;i++)	y[i][j]=col[i];
	}

	for(int i=0;i<N;i++)	for(int j=0;j<N;j++)	a[i][j]=y[i][j];
	if(flag==ON)	for(int i=0;i<N;i++)	for(int j=0;j<N;j++)	a[i][j]/=pow(10.0,fig);

	for(int i=0;i<N;i++)	delete[]	y[i];
	delete[]	y;
	delete[]	col;
	delete[]	index;
}




void ludcmp(double **a,int N,int *index,double *d)
{
	int imax=0;
	double max=0,temp=0,sum=0,dum=0;
	*d=1.0;
	vector<double> buf(N);

	for(int i=0;i<N;i++)	buf[i]=0;

	for(int i=0;i<N;i++)
	{
		max=0.0;
		for(int j=0;j<N;j++)
		{
			temp=fabs(a[i][j]);
			if(temp>max)	max=temp;
		}
		if(max==0.0)
		{
			cout<<"���ٍs��ł���"<<endl;
		}
		buf[i]=1.0/max;
	}

	for(int j=0;j<N;j++)
	{
		for(int i=0;i<j;i++)
		{
			sum=a[i][j];
			for(int k=0;k<i;k++)	sum-=a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		max=0.0;
		for(int i=j;i<N;i++)
		{
			sum=a[i][j];
			for(int k=0;k<j;k++)	sum-=a[i][k]*a[k][j];
			a[i][j]=sum;
			dum=buf[i]*fabs(sum);
			if(dum>=max)
			{
				max=dum;
				imax=i;
			}
		}

		if(j!=imax)
		{
			for(int k=0;k<N;k++)
			{
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d=-(*d);
			buf[imax]=buf[j];
		}
		index[j]=imax;
		if(fabs(a[j][j])<DBL_EPSILON)
		{
			cout<<"error:a["<<j<<"]["<<j<<"]=0"<<endl;
			a[j][j]=0;
		}
		if(j!=N-1)
		{
			dum=1.0/(a[j][j]);
			for(int i=j+1;i<N;i++)a[i][j]*=dum;
		}
	}
	buf.clear();
}

void lubksb(double **a,int N,int *index,double b[])
{
	int ip=0,ii=0,count=0;
	double sum=0;

	for(int i=0;i<N;i++)
	{
		ip=index[i];
		sum=b[ip];
		b[ip]=b[i];
		if(count)	for(int j=ii;j<i;j++)	sum-=a[i][j]*b[j];
		else if(sum)
		{
			count++;
			ii=i;
		}
		b[i]=sum;
	}

	for(int i=N-1;i>=0;i--)
	{
		sum=b[i];
		for(int j=i+1;j<N;j++)	sum-=a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}



//���͂ȂǗ��q�̎������R���^�[�}�ŕ\������֐�
void momentum_movie_AVS(mpsconfig &CON,int t,vector<mpselastic> PART,vector<hyperelastic> HYPER,int rigid_number)
{
	//�Q�l�ɂ��Ă��鏑����microAVS�̃w���v�ł��Ȃ��̃f�[�^�́H���u��\���i�q�^�f�[�^�i�A�X�L�[�j�̏����v
	int h_num=HYPER.size();
	double TIME=CON.get_step()*CON.get_dt();
	double le=CON.get_distancebp();
	int STEP=CON.get_step()/CON.get_interval()+1;		//�o�͂��鑍�X�e�b�v��
	int step;

	if(CON.get_interval()==1)	step=t/CON.get_interval();
	else step=t/CON.get_interval()+1;
	

	if(t==1) 
	{
		ofstream fp("momentum.inp", ios::trunc);			
		fp<<STEP<<endl;//���X�e�b�v��
		fp<<"data_geom"<<endl;
		fp.close();
	}

	//main�t�@�C����������
	ofstream fp("momentum.inp",ios :: app);
	fp<<"step"<<step<<" TIME="<<TIME<<endl;

	//fp<<step<<endl;
	
	//fp<<"data_geom"<<endl;
	//fp<<"step1"<<endl;
	//fp<<"step"<<t/CON->get_interval()+1<<" TIME="<<TIME<<endl;
	fp<<h_num<<" "<<h_num<<endl;	//�ߓ_���Ɨv�f���o��
	
	//�ߓ_�ԍ��Ƃ��̍��W�̏o�� 
	for(int i=0;i<h_num;i++) fp<<i<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
	
	//�v�f�ԍ��Ɨv�f�`��̎�ށA�����ėv�f���\������ߓ_�ԍ��o��
	for(int i=0;i<h_num;i++)	fp<<i<<"  0 pt "<<i<<endl;

	//fp<<"2 3"<<endl;//�ߓ_�̏��ʂ�2�ŁA�v�f�̏��ʂ�3�Ƃ������ƁB
	fp<<"7 0"<<endl;//�ߓ_�̏��ʂ�8�ŁA�v�f�̏��ʂ�0�Ƃ������ƁB
	fp<<"7 1 1 1 1 1 1 1"<<endl;	//���̍s�̏ڍׂ̓w���v���Q��
	//fp<<"8 1 1 1 1 1 1 1 1"<<endl;	//���̍s�̏ڍׂ̓w���v���Q��
	fp<<"p_x,"<<endl;
	fp<<"p_y,"<<endl;
	fp<<"p_z,"<<endl;
	fp<<"lambda,"<<endl;
	fp<<"ap_x,"<<endl;
	fp<<"ap_y,"<<endl;
	fp<<"ap_z,"<<endl;
	//fp<<"P,N/m^2"<<endl;
	//fp<<"value1,??"<<endl;

	//�e�ߓ_�̏��l����
	for(int i=0;i<h_num;i++)
	{
		fp<<i<<" "<<HYPER[i].p[A_X]<<" "<<HYPER[i].p[A_Y]<<" "<<HYPER[i].p[A_Z]<<" "<<HYPER[i].lambda<<" "<<HYPER[i].ang_p[A_X]<<" "<<HYPER[i].ang_p[A_Y]<<" "<<HYPER[i].ang_p[A_Z]<<endl;
		//fp<<i<<" "<<NODE[i].depth<<" "<<NODE[i].L/le<<" "<<NODE[i].potential<<" "<<NODE[i].Fs<<" "<<endl;
	}
	fp.close();
}


void contact_judge_hyper(mpsconfig CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,int t)
{
	cout<<"�ǉe���v�Z_���z���f��";

	//�A���S���Y��
	// 0. i���ӂ̗��q�����x�����������ꍇ�C�e�����a���ɂ��闱�q��T�����C�ȉ����s��
	// 1. �u�ڐG�̉\�������闱�q�v�i(PART[j].PND>PART[j].PND0)���^�H�j�𒲂ׂ�
	// 2. ���͂�u��
	// 3. �����z�u�̗��q�Əd�����Ȃ��悤�ɐڐG�̉\�������闱�q�Ƃ̊Ԃŗ͂��v�Z����
	int dim=3;
	double r=CON.get_h_dis();
	double le=CON.get_distancebp();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	double Dt=CON.get_dt();
	int p_num=PART.size();
	int h_num=HYPER.size();
	int w_num=p_num-h_num;

	double **w=new double *[h_num];
	double **dis=new double *[h_num];
	double **NEI_w=new double*[h_num];

	for(int i=0;i<h_num;i++)
	{
		w[i]=new double [w_num];
		dis[i]=new double [w_num];
		NEI_w[i]=new double[w_num];
	}

	if(t!=1)	for(int i=h_num;i<w_num+h_num;i++)
	{
		PART[i].r[0]+=PART[i].u[0]*Dt;
		PART[i].r[1]+=PART[i].u[1]*Dt;
		PART[i].r[2]+=PART[i].u[2]*Dt;
	}
		
//		if(PART[i].type==ELASTIC || PART[i].type==MAGELAST)
		for(int i=0;i<h_num;i++)
		{
			int N_w=0;
			double pnd=0;

			for(int j=0;j<w_num;j++)
			{
				int k=j+h_num;
			//�����z�u�̗��q�Ƃ͕��ʂɈ��͌��z���v�Z�i���͌v�Z�j
			//else if(PART[i].contact==true)
				double dis_temp0=0;
				double dis_temp1=0;
				double dis_temp=0;
				for(int D=0;D<DIMENSION;D++)
				{
					dis_temp0=PART[k].r[D]-PART[i].r[D];
					dis_temp1+=dis_temp0*dis_temp0;
				}
				dis_temp=sqrt(dis_temp1);
				if(dis_temp<r)
				{
					dis[i][N_w]=dis_temp;
					w[i][N_w]=kernel(r,dis_temp);
					pnd+=w[i][N_w];
					NEI_w[i][N_w]=k;
					N_w++;				
				}
				//���݈ʒu�ł̎��ӗ��q�����擾			
			}
			if(N_w>0)
			{
				for(int j=0;j<h_num;j++)
				{
					if(j!=i)
					{					
						double qiin[DIMENSION];
						for(int D=0;D<DIMENSION;D++)	qiin[D]=PART[j].r[D]-PART[i].r[D];
						double dis0=sqrt(qiin[A_X]*qiin[A_X]+qiin[A_Y]*qiin[A_Y]+qiin[A_Z]*qiin[A_Z]);
						double wiin=kernel(r,dis0);
						pnd+=wiin;
					}
				}
				HYPER[i].pnd=pnd;

				double gra_accel_i[DIMENSION];
				for(int D=0;D<DIMENSION;D++)	gra_accel_i[D]=0;

				for(int j=0;j<N_w;j++)
				{
					int k=NEI_w[i][j];
					double rij[DIMENSION];
					double accel_i[DIMENSION];
					for(int D=0;D<DIMENSION;D++)
					{
						rij[D]=0;
						accel_i[D]=0;

						rij[D]=PART[k].r[D]-PART[i].r[D];
						accel_i[D]*=rij[D]*w[i][j]/(dis[i][j]*dis[i][j]);
						gra_accel_i[D]+=accel_i[D];
					}
				}
				for(int D=0;D<DIMENSION;D++)	HYPER[i].p[D]-=Dt*mi*3/pnd*gra_accel_i[D];
			}
		}

	for(int i=0;i<h_num;i++)
	{
		delete[] w[i];
		delete[] dis[i];
		delete[] NEI_w[i];
	}
	delete[] w;
	delete[] dis;
	delete[]NEI_w;
	cout<<"----------OK"<<endl;
}

void contact_judge_hyper2(mpsconfig CON, vector<mpselastic> &PART, vector<hyperelastic> &HYPER, int hyper_number, int t)
{
	//�A���S���Y��
	// 0. i���ӂ̗��q�����x�����������ꍇ�C�e�����a���ɂ��闱�q��T�����C�ȉ����s��
	// 1. �u�ڐG�̉\�������闱�q�v�i(PART[j].PND>PART[j].PND0)���^�H�j�𒲂ׂ�
	// 2. ���͂�u��
	// 3. �����z�u�̗��q�Əd�����Ȃ��悤�ɐڐG�̉\�������闱�q�Ƃ̊Ԃŗ͂��v�Z����

	cout<<"�ǉe���v�Z_�����֐�";

	int h_num=hyper_number;
	int p_num=PART.size();
	double r=CON.get_h_dis();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	double Dt=CON.get_dt();
	
	double y_min=(atan(-10.0)+PI*0.5)/PI;
	double y_max=(atan(10.0)+PI*0.5)/PI;

	//�ǂ����ʂ̏ꍇ
	/*
	if(t==1)	calculation_vec_norm(PART,HYPER,h_num,p_num,t);
	else if(t!=1||PART[h_num].u[A_Z]!=0)
	{
		for(int i=h_num;i<p_num;i++)	for(int D=0;D<DIMENSION;D++)	PART[i].r[D]+=PART[i].u[D]*Dt;
		calculation_vec_norm(PART,HYPER,h_num,p_num,t);
	}*/

	stringstream s;
	s<<"./Wall/dp"<<t<<".dat";
	ofstream dp(s.str());

	ofstream s_pnd("pnd.csv", ios::app);
	if(t==1)
	{
		s_pnd<<"t"<<",";
		for(int i=0;i<h_num;i++)	s_pnd<<i<<",";
	}
	s_pnd<<t<<",";

	double vec_norm[DIMENSION]={0,0,1};
	//for(int i=0;i<h_num;i++)	for(int D=0;D<DIMENSION;D++)	HYPER[h_num].vec_norm[D]=vec_norm[D];

	for(int i=0;i<h_num;i++)
	{
		double pnd=0;
		double max_dp[DIMENSION]={-2*HYPER[i].p[0]*vec_norm[0],	-2*HYPER[i].p[1]*vec_norm[1],	-2*HYPER[i].p[2]*vec_norm[2]};
		double dp_i[DIMENSION]={0,	0,	0};

		for(int j=h_num;j<p_num;j++)
		{
			double vec_dis[DIMENSION]={PART[j].r[0]-PART[i].r[0],	PART[j].r[1]-PART[i].r[1],	PART[j].r[2]-PART[i].r[2]};
			double dis=0;
			double	in_vec=vec_dis[0]*vec_dis[0]+vec_dis[1]*vec_dis[1]+vec_dis[2]*vec_dis[2];
			dis=sqrt(in_vec);
			if(dis<r)
			{
				double w=kernel(r,dis);
				pnd+=w;
				
				//���������̋����v�Z
				double h=vec_norm[0]*vec_dis[0]+vec_norm[1]*vec_dis[1]+vec_norm[2]*vec_dis[2];
				double avs_h=fabs(h);
				double x=10-avs_h/r*20;
				double y=((atan(x)+PI*0.5)/PI-y_min)/(y_max-y_min);
				double dp_j[DIMENSION]={max_dp[0]*y,	max_dp[1]*y,	max_dp[2]*y};
				
				dp_i[A_X]+=dp_j[A_X]*w;
				dp_i[A_Y]+=dp_j[A_Y]*w;
				dp_i[A_Z]+=dp_j[A_Z]*w;
			}			
		}
		if(pnd!=0 && HYPER[i].p[A_Z]<=0)
		{
			HYPER[i].p[A_X]+=dp_i[A_X]/pnd;
			HYPER[i].p[A_Y]+=dp_i[A_Y]/pnd;
			HYPER[i].p[A_Z]+=dp_i[A_Z]/pnd;
		}

		s_pnd<<pnd<<",";
		dp<<i<<","<<dp_i[A_X]<<","<<dp_i[A_Y]<<","<<dp_i[A_Z]<<","<<endl;
	}
	s_pnd<<endl;
	s_pnd.close();
	dp.close();
	cout<<"----------OK"<<endl;
}

void calculation_vec_norm(vector<mpselastic> PART, vector<hyperelastic> &HYPER,int rigid_number,int t)
{
	int p_num=HYPER.size();
	int h_num=p_num-rigid_number;

	//�@���x�N�g���v�Z	�ǂ����ʂł��邱�ƑO��
	int id_norm[3];
	int count_max=0;
	double maxZ=0;
	double maxY=0;
	double maxX=0;
	for(int i=h_num;i<p_num;i++)
	{			
		if(PART[i].surface==ON)
		{
			if(maxZ<PART[i].r[A_Z])
			{
				maxZ=PART[i].r[A_Z];
				maxY=PART[i].r[A_Y];
				maxX=PART[i].r[A_X];
			}
			else if(maxZ=PART[i].r[A_Z])
			{
				if(maxY!=PART[i].r[A_Y]&&maxX!=PART[i].r[A_X])
				{
					maxZ=PART[i].r[A_Z];
					maxY=PART[i].r[A_Y];
					maxX=PART[i].r[A_X];
					id_norm[count_max]=i;
					cout<<"���q�ԍ�"<<id_norm[count_max]<<"Z���W"<<maxZ<<endl;
					count_max++;
				}
			}
		}
		if(count_max==3)	break;
	}

	double *r=new double [3*3];
	double *vec=new double [DIMENSION];
	for(int D=0;D<DIMENSION;D++)
	{
		vec[D]=0;
		for(int D2=0;D2<DIMENSION;D2++)	r[D*DIMENSION+D2]=0;
	}

	double r_n[DIMENSION];
	for(int D=0;D<DIMENSION;D++)
	{					
		r[0*DIMENSION+D]=PART[id_norm[0]].r[D];
		r[1*DIMENSION+D]=PART[id_norm[1]].r[D];
		r[2*DIMENSION+D]=PART[id_norm[2]].r[D];
		r_n[D]=(r[0*DIMENSION+D]+r[1*DIMENSION+D]+r[2*DIMENSION+D])/3;
		cout<<r_n[D];
	}
	cout<<endl;
	for(int nn=0;nn<3;nn++)	for(int D=0;D<DIMENSION;D++)	r[nn*DIMENSION+D]-=r_n[D];
	for(int nn=0;nn<3;nn++)	for(int D=0;D<DIMENSION;D++)	cout<<r[nn*DIMENSION+D]<<endl;;

	gauss(r,vec,DIMENSION);
	cout<<"x"<<vec[A_X]<<"y"<<vec[A_Y]<<"z"<<vec[A_Z]<<endl;
	double in_vec=vec[A_X]*vec[A_X]+vec[A_Y]*vec[A_Y]+vec[A_Z]*vec[A_Z];
	double ab_vec=sqrt(in_vec);
	for(int D=0;D<DIMENSION;D++)	HYPER[h_num].vec_norm[D]=vec[D]/ab_vec;
	cout<<"x"<<HYPER[h_num].vec_norm[A_X]<<"y"<<HYPER[h_num].vec_norm[A_Y]<<"z"<<HYPER[h_num].vec_norm[A_Z]<<endl;
	
	delete[]	r;
	delete[]	vec;
}


void output_hyper_data(vector<mpselastic> PART,vector<hyperelastic> HYPER,vector<hyperelastic2> HYPER1,int rigid_number,int t)
{
	int h_num=HYPER.size();
	int r_num=rigid_number;

	if(t==1)
	{
		////�v�Z�����e�萔�̏o��
		ofstream ai("Ai.csv");
		ofstream inai("inverse_Ai.csv");
		ofstream aiin("aiin.csv");
		ofstream n0("noij.csv");
		ofstream wiin("wiin.csv");

		aiin<<"i"<<","<<"j"<<","<<"x"<<","<<"y"<<","<<"z"<<endl;
		n0<<"i"<<","<<"j"<<","<<"x"<<","<<"y"<<","<<"z"<<endl;

		for(int j=0;j<h_num;j++)	wiin<<j<<",";
		wiin<<endl;

		for(int i=0;i<h_num;i++)
		{
			inai<<i<<",";
			inai<<HYPER[i].inverse_Ai[0][A_X]<<","<<HYPER[i].inverse_Ai[0][A_Y]<<","<<HYPER[i].inverse_Ai[0][A_Z]<<","<<endl;		
			inai<<","<<HYPER[i].inverse_Ai[1][A_X]<<","<<HYPER[i].inverse_Ai[1][A_Y]<<","<<HYPER[i].inverse_Ai[1][A_Z]<<","<<endl;	
			inai<<","<<HYPER[i].inverse_Ai[2][A_X]<<","<<HYPER[i].inverse_Ai[2][A_Y]<<","<<HYPER[i].inverse_Ai[2][A_Z]<<","<<endl;	

			ai<<i<<",";
			ai<<HYPER[i].Ai[0][A_X]<<","<<HYPER[i].Ai[0][A_Y]<<","<<HYPER[i].Ai[0][A_Z]<<","<<endl;
			ai<<","<<HYPER[i].Ai[1][A_X]<<","<<HYPER[i].Ai[1][A_Y]<<","<<HYPER[i].Ai[1][A_Z]<<","<<endl;
			ai<<","<<HYPER[i].Ai[2][A_X]<<","<<HYPER[i].Ai[2][A_Y]<<","<<HYPER[i].Ai[2][A_Z]<<","<<endl;
			
			for(int j=0;j<h_num;j++)
			{
				aiin<<i<<","<<j<<","<<HYPER1[i*h_num+j].aiin[A_X]<<","<<HYPER1[i*h_num+j].aiin[A_Y]<<","<<HYPER1[i*h_num+j].aiin[A_Z]<<","<<endl;
				n0<<i<<","<<j<<","<<HYPER1[i*h_num+j].n0ij[A_X]<<","<<HYPER1[i*h_num+j].n0ij[A_Y]<<","<<HYPER1[i*h_num+j].n0ij[A_Z]<<","<<endl;
				wiin<<i<<","<<HYPER1[i*h_num+j].wiin<<",";
			}
			wiin<<endl;
		}
		ai.close();
		inai.close();
		n0.close();
		wiin.close();
	}

	stringstream ss_dgdq;
	ss_dgdq<<"./DgDq/DgDq"<<t<<".csv";
	ofstream dg(ss_dgdq.str());
	ofstream d_p("d_P.csv", ios::app);
	ofstream h_p("h_P.csv", ios::app);
	ofstream lam("lambda.csv", ios::app);
	ofstream stress("stress.csv", ios::app);
	ofstream p("P.csv", ios::app);
	ofstream J("J.csv", ios::app);
	ofstream ti_Fi("ti_Fi.csv", ios::app);
	ofstream Fi("Fi.csv", ios::app);

	if(t==1)
	{
		p<<"t"<<",";		
		d_p<<"t"<<",";
		h_p<<"t"<<",";
		lam<<"t"<<",";
		stress<<"t"<<",";
		ti_Fi<<"t"<<",";
		Fi<<"t"<<",";
		J<<"t"<<",";

		for(int i=0;i<h_num-r_num;i++)
		{
			p<<i<<","<<","<<",";
			d_p<<i<<","<<","<<",";
			h_p<<i<<","<<","<<",";
			lam<<i<<",";
			stress<<i<<","<<","<<",";
			ti_Fi<<i<<","<<","<<",";
			Fi<<i<<","<<","<<",";
			J<<i<<",";
		}
		for(int i=h_num-r_num;i<h_num;i++)
		{
			p<<i<<","<<","<<",";
			d_p<<i<<","<<","<<",";
			h_p<<i<<","<<","<<",";
			lam<<i<<",";
		}

		p<<endl;
		d_p<<endl;
		h_p<<endl;
		lam<<endl;
		stress<<endl;
		ti_Fi<<endl;
		Fi<<endl;
		J<<endl;
	}

	p<<t<<",";
	d_p<<t<<",";
	h_p<<t<<",";
	lam<<t<<",";
	stress<<t<<",";
	ti_Fi<<t<<",";
	Fi<<t<<",";
	J<<t<<",";
	for(int i=0;i<h_num-r_num;i++)
	{
		
		p<<HYPER[i].p[A_X]<<","<<HYPER[i].p[A_Y]<<","<<HYPER[i].p[A_Z]<<",";
		d_p<<HYPER[i].differential_p[A_X]<<","<<HYPER[i].differential_p[A_Y]<<","<<HYPER[i].differential_p[A_Z]<<","<<endl;
		h_p<<HYPER[i].half_p[A_X]<<","<<HYPER[i].half_p[A_Y]<<","<<HYPER[i].half_p[A_Z]<<","<<endl;
		lam<<HYPER[i].lambda<<",";

		stress<<HYPER[i].stress[A_X][A_X]<<","<<HYPER[i].stress[A_X][A_Y]<<","<<HYPER[i].stress[A_X][A_Z]<<",";
		ti_Fi<<HYPER[i].t_inverse_Fi[A_X][A_X]<<","<<HYPER[i].t_inverse_Fi[A_X][A_Y]<<","<<HYPER[i].t_inverse_Fi[A_X][A_Z]<<",";	
		Fi<<HYPER[i].Fi[0][0]<<","<<HYPER[i].Fi[0][1]<<","<<HYPER[i].Fi[0][2]<<",";
		J<<HYPER[i].J<<",";
		for(int j=0;j<h_num-r_num;j++)	dg<<i<<","<<j<<","<<HYPER1[j*h_num+i].DgDq[A_X]<<","<<HYPER1[j*h_num+i].DgDq[A_Y]<<","<<HYPER1[j*h_num+i].DgDq[A_Z]<<","<<endl;
	}
	for(int i=h_num-r_num;i<h_num;i++)
	{
		p<<HYPER[i].p[A_X]<<","<<HYPER[i].p[A_Y]<<","<<HYPER[i].p[A_Z]<<",";
		d_p<<HYPER[i].differential_p[A_X]<<","<<HYPER[i].differential_p[A_Y]<<","<<HYPER[i].differential_p[A_Z]<<","<<endl;
		h_p<<HYPER[i].half_p[A_X]<<","<<HYPER[i].half_p[A_Y]<<","<<HYPER[i].half_p[A_Z]<<","<<endl;
		lam<<HYPER[i].lambda<<",";
	}
	stress<<endl<<",";
	ti_Fi<<endl<<",";
	Fi<<endl<<",";
	for(int i=0;i<h_num-r_num;i++)
	{
		stress<<HYPER[i].stress[A_Y][A_X]<<","<<HYPER[i].stress[A_Y][A_Y]<<","<<HYPER[i].stress[A_Y][A_Z]<<",";
		ti_Fi<<HYPER[i].t_inverse_Fi[A_Y][A_X]<<","<<HYPER[i].t_inverse_Fi[A_Y][A_Y]<<","<<HYPER[i].t_inverse_Fi[A_Y][A_Z]<<",";
		Fi<<HYPER[i].Fi[1][0]<<","<<HYPER[i].Fi[1][1]<<","<<HYPER[i].Fi[1][2]<<",";
	}
	stress<<endl<<",";
	ti_Fi<<endl<<",";
	Fi<<endl<<",";
	for(int i=0;i<h_num-r_num;i++)
	{
		stress<<HYPER[i].stress[A_Z][A_X]<<","<<HYPER[i].stress[A_Z][A_Y]<<","<<HYPER[i].stress[A_Z][A_Z]<<",";
		ti_Fi<<HYPER[i].t_inverse_Fi[A_Z][A_X]<<","<<HYPER[i].t_inverse_Fi[A_Z][A_Y]<<","<<HYPER[i].t_inverse_Fi[A_Z][A_Z]<<",";
		Fi<<HYPER[i].Fi[2][0]<<","<<HYPER[i].Fi[2][1]<<","<<HYPER[i].Fi[2][2]<<",";
	}

	p<<endl;
	d_p<<endl;
	h_p<<endl;
	lam<<endl;
	stress<<endl;
	ti_Fi<<endl;
	Fi<<endl;
	J<<endl;

	dg.close();
	stress.close();
	lam.close();
	ti_Fi.close();
	Fi.close();
	p.close();
	d_p.close();
	h_p.close();
	J.close();
}

void output_newton_data1(double *fx, double *DfDx, double *n_rx, double *n_ry, double *n_rz,int hyper_number,int count,int t)
{
	int h_num=hyper_number;
	stringstream ss_r;
	ss_r<<"./Newton_raphson/position"<<t<<".csv";
	stringstream ss_Df;
	ss_Df<<"./Newton_raphson/DfDx "<<t<<".csv";
	stringstream ss_fx;
	ss_fx<<"./Newton_raphson/fx"<<t<<".csv";
		
	if(count==1)
	{
		ofstream init0(ss_r.str(), ios::trunc);
		ofstream init1(ss_Df.str(), ios::trunc);
		ofstream init3(ss_fx.str(), ios::trunc);
	
		init0.close();
		init1.close();
		init3.close();
	}

	ofstream r(ss_r.str(), ios::app);
	ofstream sfx(ss_fx.str(), ios::app);
	ofstream Df(ss_Df.str(), ios::app);

	if(count==1)
	{
		r<<"count"<<",";
		sfx<<"count"<<",";
		for(int i=0;i<h_num;i++)
		{
			r<<i<<","<<","<<",";
			sfx<<i<<",";
		}
		r<<endl;
		sfx<<endl;
	}
	r<<count<<",";
	sfx<<count<<",";

	Df<<"count"<<","<<count<<endl;
	for(int i=0;i<h_num;i++)
	{
		r<<n_rx[i]<<","<<n_ry[i]<<","<<n_rz[i]<<","<<endl;
		sfx<<fx[i]<<",";
		for(int j=0;j<h_num;j++) Df<<DfDx[i*h_num+j]<<",";
		Df<<endl;
	}
	sfx<<endl;

	r.close();
	sfx.close();

	Df.close();
}

void output_newton_data2(double E, double *XX,int count,int hyper_number, int t)
{
	int h_num=hyper_number;

	stringstream ss_E;
	ss_E<<"./Newton_raphson/E"<<t<<".dat";
	
	stringstream ss_lam;
	ss_lam<<"./Newton_raphson/lambda"<<t<<".dat";
		
	if(count==1)
	{
		ofstream init0(ss_E.str(), ios::trunc);
		ofstream init1(ss_lam.str(), ios::trunc);
	
		init0.close();
		init1.close();
	}

	ofstream e(ss_E.str(), ios::app);
	ofstream lam(ss_lam.str(), ios::app);

	if(count==1)
	{
		e<<"������"<<","<<"E"<<endl;
		lam<<"������"<<","<<"lambda"<<endl;
		for(int i=0;i<h_num;i++)	lam<<","<<i;
		lam<<endl;
	}
	
	e<<count<<","<<E<<endl;

	lam<<count;
	for(int i=0;i<h_num;i++)	lam<<","<<XX[i];
	lam<<endl;

	e.close();
	lam.close();
}

void output_energy(mpsconfig CON, vector<mpselastic> PART, vector<hyperelastic> HYPER, int rigid_number, int t)
{
	cout<<"�e���|�e���V�����v�Z";
	int h_num=HYPER.size();
	int r_num=rigid_number;
	double d_Fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double c10=CON.get_c10();
	double c01=CON.get_c01();
	vector<double>	W;
	for(int i=0;i<h_num-r_num;i++)	W.push_back(0);
	double dC[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double dC2[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

	for(int i=0;i<h_num-r_num;i++)
	{
		double J=HYPER[i].J;
		if(J<0){
			d_Fi[0][0]=-1/pow(-J,1/3)*HYPER[i].Fi[0][0];	d_Fi[0][1]=-1/pow(-J,1/3)*HYPER[i].Fi[0][1];	d_Fi[0][2]=-1/pow(-J,1/3)*HYPER[i].Fi[0][2];
			d_Fi[1][0]=-1/pow(-J,1/3)*HYPER[i].Fi[1][0];	d_Fi[1][1]=-1/pow(-J,1/3)*HYPER[i].Fi[1][1];	d_Fi[1][2]=-1/pow(-J,1/3)*HYPER[i].Fi[1][2];
			d_Fi[2][0]=-1/pow(-J,1/3)*HYPER[i].Fi[2][0];	d_Fi[2][1]=-1/pow(-J,1/3)*HYPER[i].Fi[2][1];	d_Fi[2][2]=-1/pow(-J,1/3)*HYPER[i].Fi[2][2];
		}
		else
		{
			d_Fi[0][0]=1/pow(J,1/3)*HYPER[i].Fi[0][0];	d_Fi[0][1]=1/pow(J,1/3)*HYPER[i].Fi[0][1];	d_Fi[0][2]=1/pow(J,1/3)*HYPER[i].Fi[0][2];
			d_Fi[1][0]=1/pow(J,1/3)*HYPER[i].Fi[1][0];	d_Fi[1][1]=1/pow(J,1/3)*HYPER[i].Fi[1][1];	d_Fi[1][2]=1/pow(J,1/3)*HYPER[i].Fi[1][2];
			d_Fi[2][0]=1/pow(J,1/3)*HYPER[i].Fi[2][0];	d_Fi[2][1]=1/pow(J,1/3)*HYPER[i].Fi[2][1];	d_Fi[2][2]=1/pow(J,1/3)*HYPER[i].Fi[2][2];
		}

		dC[0][0]=d_Fi[0][0]*d_Fi[0][0]+d_Fi[1][0]*d_Fi[1][0]+d_Fi[2][0]*d_Fi[2][0];
		dC[0][1]=d_Fi[0][0]*d_Fi[0][1]+d_Fi[1][0]*d_Fi[1][1]+d_Fi[2][0]*d_Fi[2][1];
		dC[0][2]=d_Fi[0][0]*d_Fi[0][2]+d_Fi[1][0]*d_Fi[1][2]+d_Fi[2][0]*d_Fi[2][2];
		dC[1][0]=d_Fi[0][1]*d_Fi[0][0]+d_Fi[1][1]*d_Fi[1][0]+d_Fi[2][1]*d_Fi[2][0];
		dC[1][1]=d_Fi[0][1]*d_Fi[0][1]+d_Fi[1][1]*d_Fi[1][1]+d_Fi[2][1]*d_Fi[2][1];
		dC[1][2]=d_Fi[0][1]*d_Fi[0][2]+d_Fi[1][1]*d_Fi[1][2]+d_Fi[2][1]*d_Fi[2][2];
		dC[2][0]=d_Fi[0][2]*d_Fi[0][0]+d_Fi[1][2]*d_Fi[1][0]+d_Fi[2][2]*d_Fi[2][0];
		dC[2][1]=d_Fi[0][2]*d_Fi[0][1]+d_Fi[1][2]*d_Fi[1][1]+d_Fi[2][2]*d_Fi[2][1];
		dC[2][2]=d_Fi[0][2]*d_Fi[0][2]+d_Fi[1][2]*d_Fi[1][2]+d_Fi[2][2]*d_Fi[2][2];

		dC2[0][0]=dC[A_X][0]*dC[0][A_X]+dC[A_X][1]*dC[1][A_X]+dC[A_X][2]*dC[2][A_X];
		dC2[0][1]=dC[A_X][0]*dC[0][A_Y]+dC[A_X][1]*dC[1][A_Y]+dC[A_X][2]*dC[2][A_Y];
		dC2[0][2]=dC[A_X][0]*dC[0][A_Z]+dC[A_X][1]*dC[1][A_Z]+dC[A_X][2]*dC[2][A_Z];
		dC2[1][0]=dC[A_Y][0]*dC[0][A_X]+dC[A_Y][1]*dC[1][A_X]+dC[A_Y][2]*dC[2][A_X];
		dC2[1][1]=dC[A_Y][0]*dC[0][A_Y]+dC[A_Y][1]*dC[1][A_Y]+dC[A_Y][2]*dC[2][A_Y];
		dC2[1][2]=dC[A_Y][0]*dC[0][A_Z]+dC[A_Y][1]*dC[1][A_Z]+dC[A_Y][2]*dC[2][A_Z];
		dC2[2][0]=dC[A_Z][0]*dC[0][A_X]+dC[A_Z][1]*dC[1][A_X]+dC[A_Z][2]*dC[2][A_X];
		dC2[2][1]=dC[A_Z][0]*dC[0][A_Y]+dC[A_Z][1]*dC[1][A_Y]+dC[A_Z][2]*dC[2][A_Y];
		dC2[2][2]=dC[A_Z][0]*dC[0][A_Z]+dC[A_Z][1]*dC[1][A_Z]+dC[A_Z][2]*dC[2][A_Z];

		double trace_dC=dC[0][0]+dC[1][1]+dC[2][2];
		double trace_dC2=dC2[0][0]+dC2[1][1]+dC2[2][2];

		double Ic=trace_dC;
		double IIc=0.50*(trace_dC*trace_dC-trace_dC2);
		W[i]=c10*(Ic-3)+c01*(IIc-3);
	}
	cout<<"----------OK"<<endl;

	ofstream e("E.csv", ios::app);
	ofstream e_T("E_T.csv", ios::app);
	ofstream e_g("E_g.csv", ios::app);
	ofstream e_W("E_W.csv", ios::app);
	ofstream e_lam("E_lam.csv", ios::app);

	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();

	if(t==1)
	{
		e<<"t"<<",";
		e_T<<"t"<<",";
		e_g<<"t"<<",";
		e_W<<"t"<<",";
		e_lam<<"t"<<",";
		for(int i=0;i<h_num-r_num;i++)
		{
			e<<i<<",";
			e_T<<i<<",";
			e_g<<i<<",";
			e_W<<i<<",";
			e_lam<<i<<",";
		}
		for(int i=h_num-r_num;i<h_num;i++)
		{
			e<<i<<"r"<<",";
			e_T<<i<<"r"<<",";
			e_g<<i<<"r"<<",";
		}

		e<<endl;
		e_T<<endl;
		e_g<<endl;
		e_W<<endl;
		e_lam<<endl;
	}

	e<<t<<",";
	e_T<<t<<",";
	e_g<<t<<",";
	e_W<<t<<",";
	e_lam<<t<<",";
	double sum_e_T=0;
	double sum_e_g=0;
	double sum_e_lam=0;
	double sum_e_W=0;
	double sum_e=0;
	double vv=0;
	double energy=0;

	for(int i=0;i<h_num-r_num;i++)
	{
		vv=HYPER[i].p[0]*HYPER[i].p[0]+HYPER[i].p[1]*HYPER[i].p[1]+HYPER[i].p[2]*HYPER[i].p[2];
		energy=0.5/mi*vv+W[i]*V+HYPER[i].lambda*(1-HYPER[i].J)*V;
		//energy=0.5/mi*vv+mi*9.8*PART[i].r[A_Z]+W[i]*V+HYPER[i].lambda*(1-HYPER[i].J)*V;
		sum_e_T+=0.5/mi*vv;
		sum_e_g+=mi*9.8*PART[i].r[A_Z];
		sum_e_lam+=HYPER[i].lambda*(1-HYPER[i].J)*V;
		sum_e_W+=W[i]*V;
		sum_e+=energy;
		e<<energy<<",";
		e_T<<0.5/mi*vv<<",";
		e_g<<mi*9.8*PART[i].r[A_Z]<<",";
		e_W<<W[i]*V<<",";
		e_lam<<HYPER[i].lambda*(1-HYPER[i].J)*V<<",";
	}
	for(int i=h_num-r_num;i<h_num;i++)
	{
		vv=HYPER[i].p[0]*HYPER[i].p[0]+HYPER[i].p[1]*HYPER[i].p[1]+HYPER[i].p[2]*HYPER[i].p[2];
		energy=0.5/mi*vv+HYPER[i].lambda*(1-HYPER[i].J)*V;
		//energy=0.5/mi*vv+mi*9.8*PART[i].r[A_Z]+W[i]*V+HYPER[i].lambda*(1-HYPER[i].J)*V;
		sum_e_T+=0.5/mi*vv;
		sum_e_g+=mi*9.8*PART[i].r[A_Z];
		sum_e+=energy;
		e<<energy<<",";
		e_T<<0.5/mi*vv<<",";
		e_g<<mi*9.8*PART[i].r[A_Z]<<",";
	}

	e<<sum_e<<endl;
	e_T<<sum_e_T<<endl;
	e_g<<sum_e_g<<endl;
	e_W<<sum_e_W<<endl;
	e_lam<<sum_e_lam<<endl;

	e.close();
	e_T.close();
	e_g.close();
	e_W.close();
	e_lam.close();
}

void transpose(double **M, double **N)
{
	N[0][0]=M[0][0];
	N[1][1]=M[1][1];
	N[2][2]=M[2][2];
	N[0][1]=M[1][0];
	N[0][2]=M[2][0];
	N[1][0]=M[0][1];
	N[1][2]=M[2][1];
	N[2][0]=M[0][2];
	N[2][1]=M[1][2];
}

hyperelastic::hyperelastic()
{
	for(int i=0;i<200;i++)
	{
		NEI0[i]=0;
	}
	N0=0;
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



