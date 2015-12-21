#include "stdafx.h"		

void calc_half_p(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,bool repetation,double **F);
void renew_lambda(mpsconfig &CON,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int t);
void calc_differential_p(mpsconfig &CON,vector<mpselastic>PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,double **F);
void calc_transposed_inverse_matrix(double **M,bool transport,bool inversion);
double calc_det(double **M,int N);
double calc_det3(double **M);
void calc_stress(mpsconfig &CON,vector<hyperelastic> &HYPER);
void calc_constant(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1);
void calc_inverse_matrix_for_NR(int N, double *a);
void newton_raphson(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int t,double **F);
void calc_F(vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1);
void calc_newton_function(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> HYPER,vector<hyperelastic2> HYPER1,double *lambda,double *fx,double *DfDx,int hyper_number,int count,int t,double **F);
void inverse(double **a,int N);
void ludcmp(double **a,int N,int *index,double *d);
void lubksb(double **a,int N,int *index,double b[]);
void momentum_movie_AVS(mpsconfig &CON,int t,vector<mpselastic> PART,vector<hyperelastic> HYPER,double **F);
void contact_judge_hyper(mpsconfig &CON,vector<mpselastic> &PART, vector<hyperelastic> &HYPER,int t);
void contact_judge_hyper2(mpsconfig CON, vector<mpselastic> &PART, vector<hyperelastic> &HYPER, int hyper_number, int t);
void output_hyper_data(vector<mpselastic> PART,vector<hyperelastic> HYPER,vector<hyperelastic2> HYPER1, int t);
void transpose(double **M,double **N);
void output_newton_data1(double *fx, double *DfDx, double *n_rx, double *n_ry, double *n_rz, int hyper_number, int count, int t);
void output_newton_data2(double E, double *XX, int hyper_number, int count, int t);
void calc_gravity(mpsconfig CON,vector<hyperelastic> &HYPER,int hyper_number);
void calculation_vec_norm(vector<mpselastic> PART, vector<hyperelastic> &HYPER, int hyper_number,int particle_number,int t);
void output_energy(mpsconfig CON, vector<mpselastic> PART, vector<hyperelastic> HYPER,int t);
void contact_judge(mpsconfig &CON, vector<mpselastic> PART,vector<hyperelastic> &HYPER,double max_h,int t);


void calc_hyper(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t,double **F)
{	
	int h_num=HYPER.size();
	int p_num=PART.size();
	double Dt=CON.get_dt();
	cout<<"Hypercalculation starts."<<endl;
	double max_h=PART[p_num-h_num].r[A_Z];
	//calc_gravity(CON,HYPER,h_num);
	
//	contact_judge(CON,PART,HYPER,max_h,t);

	if(t==1)
	{
		for(int i=0;i<h_num;i++)	for(int D=0;D<DIMENSION;D++)	PART[i].q0[D]=0;
		for(int i=0;i<h_num;i++)	for(int D=0;D<DIMENSION;D++)	PART[i].q0[D]=PART[i].r[D];
		for(int k=p_num-h_num+1;k<p_num;k++)	if(max_h<PART[k].r[A_Z])	max_h=PART[k].r[A_Z];
		calc_constant(CON,PART,HYPER,HYPER1);
		calc_stress(CON,HYPER);
	}
	
	if(CON.get_flag_vis()==ON)	calc_vis_f(CON,PART,HYPER,HYPER1,0,t);

	if(t==1 || t%CON.get_interval()==0)
	{
		output_hyper_data(PART,HYPER,HYPER1,t);
		momentum_movie_AVS(CON,t,PART,HYPER,F);
		output_energy(CON,PART,HYPER,t);
	}


	newton_raphson(CON,PART,HYPER,HYPER1,t,F);

	calc_half_p(CON,PART,HYPER,HYPER1,0,F);

	calc_F(PART,HYPER,HYPER1);

	calc_stress(CON,HYPER);
	
	calc_differential_p(CON,PART,HYPER,HYPER1,F);

	renew_lambda(CON,HYPER,HYPER1,t);


	calc_half_p(CON,PART,HYPER,HYPER1,1,F);

	cout<<"Hypercalculation ends."<<endl;
}

void calc_constant(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1)
{
	cout<<"�����l�v�Z";

	double le=CON.get_distancebp();
	double r=CON.get_h_dis();
	double Dt=CON.get_dt();
	double V=get_volume(&CON);	//�l�����K�v��������Ȃ�
	double mi=V*CON.get_hyper_density();
	int h_num=HYPER.size();
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
	if(model==21)
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
		for(int i=0;i<h_num;i++)
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
	for(int i=0;i<h_num;i++)
	{
		int N=0;
		double dis=0;
		double wiin=0;
		double aiin[DIMENSION];
		for(int j=0;j<h_num;j++)
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
				HYPER[i].NEI[N]=j;
				N++;
			}
			else	wiin=0;

			HYPER[i].pnd+=wiin;
			HYPER1[i*h_num+j].wiin=wiin;
		}
		HYPER[i].N=N;
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
		Ai[0][0]=0;	Ai[0][1]=0;	Ai[0][2]=0;
		Ai[1][0]=0;	Ai[1][1]=0;	Ai[1][2]=0;
		Ai[2][0]=0;	Ai[2][1]=0;	Ai[2][2]=0;

		int Ni=HYPER[i].N;
		double a[DIMENSION];
		for(int j=0;j<Ni;j++)	
		{
			int k=HYPER[i].NEI[j];
			double w=HYPER1[i*h_num+k].wiin;		
			double a[DIMENSION]={HYPER1[i*h_num+k].aiin[A_X],	a[A_Y]=HYPER1[i*h_num+k].aiin[A_Y], a[A_Z]=HYPER1[i*h_num+k].aiin[A_Z]};
			
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
		int Ni=HYPER[i].N;

		double p_n0ij[DIMENSION]={0,0,0};
		for(int j=0;j<Ni;j++)
		{
			int k=HYPER[i].NEI[j];			
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
	for(int i=0;i<h_num;i++)
	{
		int Ni=HYPER[i].N;
		for(int j=0;j<Ni;j++)
		{
			int k=HYPER[i].NEI[j];
			HYPER1[k*h_num+i].DgDq[A_X]=HYPER[k].J*(HYPER[k].t_inverse_Fi[A_X][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].t_inverse_Fi[A_X][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].t_inverse_Fi[A_X][2]*HYPER1[i*h_num+k].n0ij[2]);
			HYPER1[k*h_num+i].DgDq[A_Y]=HYPER[k].J*(HYPER[k].t_inverse_Fi[A_Y][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].t_inverse_Fi[A_Y][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].t_inverse_Fi[A_Y][2]*HYPER1[i*h_num+k].n0ij[2]);
			HYPER1[k*h_num+i].DgDq[A_Z]=HYPER[k].J*(HYPER[k].t_inverse_Fi[A_Z][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].t_inverse_Fi[A_Z][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].t_inverse_Fi[A_Z][2]*HYPER1[i*h_num+k].n0ij[2]);		
		}
		HYPER1[i*h_num+i].DgDq[A_X]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_X][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_X][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_X][2]*HYPER1[i*h_num+i].n0ij[2]);
		HYPER1[i*h_num+i].DgDq[A_Y]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_Y][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_Y][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_Y][2]*HYPER1[i*h_num+i].n0ij[2]);
		HYPER1[i*h_num+i].DgDq[A_Z]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_Z][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_Z][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_Z][2]*HYPER1[i*h_num+i].n0ij[2]);
	}	
	cout<<"----------OK"<<endl;
 }


/////�j���[�g�����t�\���@ 
void newton_raphson(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int t,double **F)
{
	/////fx(N*N�s��̊e����)��1�����z��Ŋi�[�A(i,j)�����Ȃ�[j*N+i]�ŎQ��
	/////DfDx(N*N�s��̊e����)��1�����z��Ŋi�[�A(i,j)�����Ȃ�[j*N+i]�ŎQ��

	int calc_type=1;//�j���[�g���t�\���̔������@ 0:�Δ������̋t�s������̂܂܋��߂�@1:���`�������𗘗p

	//pn=2;//test,�Ƃ肠����2���łƂ��邩�ǂ����m�F 
	//////////////////�@f1(x1,x2) = x1^2 + x2^2 -5 = 0 f2(x1,x2) = x1^2/9+ x2^2 -1 = 0  http://homepage1.nifty.com/gfk/excel_newton_ren.htm

	int h_num=HYPER.size();
	double *fx=new double [h_num];//�֐��l�B
	double *DfDx=new double [h_num*h_num];//�֐��̕Δ����l�B
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

	for(int i=0;i<h_num;i++)
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

		calc_newton_function(CON,PART,HYPER,HYPER1,XX,fx,DfDx,h_num,count,t,F);


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
			gauss(DfDx,fx,h_num);
			for(int i=0;i<h_num;i++)	XX[i]-=fx[i];/*0.5*mi/(Dt*Dt)*V*fx[i];
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
			//output_newton_data2(E,XX,h_num,count,t);

		}
		if(count>CON.get_nr())	break;
		else if(dec_flag==ON)	if(E_old-E<0)	break;	
	}
//	end=clock();
//	newton_t=(end-start)/CLOCKS_PER_SEC;

	cout<<"��������";

	for(int i=0;i<h_num;i++) HYPER[i].lambda=XX[i];

//	for(int i=0;i<N;i++)	cout<<"lambda["<<i<<"]="<<HYPER[i].lambda<<endl;
	delete[]	fx;
	delete[]	DfDx;
	delete[]	XX;
	delete[]	XX_old;
	cout<<"---------- OK"<<endl;
}


void calc_newton_function(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> HYPER,vector<hyperelastic2> HYPER1,double *lambda,double *fx,double *DfDx,int hyper_number,int count,int t,double **F)
{
	clock_t t3=clock();

	int h_num=hyper_number;
	int flag_vis=CON.get_flag_vis();
	bool flag_FEM=CON.get_FEM_flag();
	int flag_G=CON.get_flag_G();
	int model_num=CON.get_model_number();
	double Dt=CON.get_dt();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	double density=CON.get_hyper_density();
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
	for(int i=0;i<h_num;i++)
	{
		if(model_num==30||model_num==23)
		{
			if(PART[i].q0[A_Z]!=0)
			{
				//half_p�̌v�Z
				double p_half_p[DIMENSION]={0,0,0};

				int Ni=HYPER[i].N;
				for(int j=0;j<Ni;j++)
				{	
					int k=HYPER[i].NEI[j];
					p_half_p[A_X]+=(HYPER[k].stress[0][0]-lambda[k])*HYPER1[k*h_num+i].DgDq[0]+HYPER[k].stress[0][1]*HYPER1[k*h_num+i].DgDq[1]+HYPER[k].stress[0][2]*HYPER1[k*h_num+i].DgDq[2];
					p_half_p[A_Y]+=HYPER[k].stress[1][0]*HYPER1[k*h_num+i].DgDq[0]+(HYPER[k].stress[1][1]-lambda[k])*HYPER1[k*h_num+i].DgDq[1]+HYPER[k].stress[1][2]*HYPER1[k*h_num+i].DgDq[2];
					p_half_p[A_Z]+=HYPER[k].stress[2][0]*HYPER1[k*h_num+i].DgDq[0]+HYPER[k].stress[2][1]*HYPER1[k*h_num+i].DgDq[1]+(HYPER[k].stress[2][2]-lambda[k])*HYPER1[k*h_num+i].DgDq[2];
				}//j�Ɋւ���for���̏I���
				p_half_p[A_X]+=(HYPER[i].stress[0][0]-lambda[i])*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[0][1]*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[0][2]*HYPER1[i*h_num+i].DgDq[2];
				p_half_p[A_Y]+=HYPER[i].stress[1][0]*HYPER1[i*h_num+i].DgDq[0]+(HYPER[i].stress[1][1]-lambda[i])*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[1][2]*HYPER1[i*h_num+i].DgDq[2];
				p_half_p[A_Z]+=HYPER[i].stress[2][0]*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[2][1]*HYPER1[i*h_num+i].DgDq[1]+(HYPER[i].stress[2][2]-lambda[i])*HYPER1[i*h_num+i].DgDq[2];
		
				//�d�͂̉e��
				if(flag_G==ON)	p_half_p[A_Z]-=9.8*mi;
				//�S�����̉e��
				if(flag_vis==ON)
				{
					p_half_p[A_X]+=HYPER[i].vis_force[A_X];
					p_half_p[A_Y]+=HYPER[i].vis_force[A_Y];
					p_half_p[A_Z]+=HYPER[i].vis_force[A_Z];
				}
				//����̍l��
				if(flag_FEM==ON)
				{
					p_half_p[A_X]+=F[A_X][i]*mi;//density;
					p_half_p[A_Y]+=F[A_Y][i]*mi;//density;
					p_half_p[A_Z]+=F[A_Z][i]*mi;//density;
				}
				//�ʒu���W�̌v�Z
				n_rx[i]=PART[i].r[A_X]+Dt*(HYPER[i].p[A_X]+Dt*0.5*p_half_p[A_X])/mi;
				n_ry[i]=PART[i].r[A_Y]+Dt*(HYPER[i].p[A_Y]+Dt*0.5*p_half_p[A_Y])/mi;
				n_rz[i]=PART[i].r[A_Z]+Dt*(HYPER[i].p[A_Z]+Dt*0.5*p_half_p[A_Z])/mi;
			}
			else
			{
				n_rx[i]=PART[i].q0[A_X];
				n_ry[i]=PART[i].q0[A_Y];
				n_rz[i]=PART[i].q0[A_Z];
			}
		}
		else
		{
			//half_p�̌v�Z
			double p_half_p[DIMENSION]={0,0,0};

			int Ni=HYPER[i].N;
			for(int j=0;j<Ni;j++)
			{	
				int k=HYPER[i].NEI[j];
				p_half_p[A_X]+=(HYPER[k].stress[0][0]-lambda[k])*HYPER1[k*h_num+i].DgDq[0]+HYPER[k].stress[0][1]*HYPER1[k*h_num+i].DgDq[1]+HYPER[k].stress[0][2]*HYPER1[k*h_num+i].DgDq[2];
				p_half_p[A_Y]+=HYPER[k].stress[1][0]*HYPER1[k*h_num+i].DgDq[0]+(HYPER[k].stress[1][1]-lambda[k])*HYPER1[k*h_num+i].DgDq[1]+HYPER[k].stress[1][2]*HYPER1[k*h_num+i].DgDq[2];
				p_half_p[A_Z]+=HYPER[k].stress[2][0]*HYPER1[k*h_num+i].DgDq[0]+HYPER[k].stress[2][1]*HYPER1[k*h_num+i].DgDq[1]+(HYPER[k].stress[2][2]-lambda[k])*HYPER1[k*h_num+i].DgDq[2];
			}//j�Ɋւ���for���̏I���
			p_half_p[A_X]+=(HYPER[i].stress[0][0]-lambda[i])*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[0][1]*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[0][2]*HYPER1[i*h_num+i].DgDq[2];
			p_half_p[A_Y]+=HYPER[i].stress[1][0]*HYPER1[i*h_num+i].DgDq[0]+(HYPER[i].stress[1][1]-lambda[i])*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[1][2]*HYPER1[i*h_num+i].DgDq[2];
			p_half_p[A_Z]+=HYPER[i].stress[2][0]*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[2][1]*HYPER1[i*h_num+i].DgDq[1]+(HYPER[i].stress[2][2]-lambda[i])*HYPER1[i*h_num+i].DgDq[2];
		
			//�d�͂̉e��
			if(flag_G==ON)	p_half_p[A_Z]-=9.8*mi;
			//�S�����̉e��
			if(flag_vis==ON)
			{
				p_half_p[A_X]+=HYPER[i].vis_force[A_X];
				p_half_p[A_Y]+=HYPER[i].vis_force[A_Y];
				p_half_p[A_Z]+=HYPER[i].vis_force[A_Z];
			}
			//����̍l��
			if(flag_FEM==ON)
			{
				p_half_p[A_X]+=F[A_X][i]*mi;//density;
				p_half_p[A_Y]+=F[A_Y][i]*mi;//density;
				p_half_p[A_Z]+=F[A_Z][i]*mi;//density;
			}
			//�ʒu���W�̌v�Z
			n_rx[i]=PART[i].r[A_X]+Dt*(HYPER[i].p[A_X]+Dt*0.5*p_half_p[A_X])/mi;
			n_ry[i]=PART[i].r[A_Y]+Dt*(HYPER[i].p[A_Y]+Dt*0.5*p_half_p[A_Y])/mi;
			n_rz[i]=PART[i].r[A_Z]+Dt*(HYPER[i].p[A_Z]+Dt*0.5*p_half_p[A_Z])/mi;
		}
	}


	////DgDq��fx�̍X�V
	for(int i=0;i<h_num;i++)
	{
		//Fi�̌v�Z
		int Ni=HYPER[i].N;
		double fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};	

		for(int in=0;in<Ni;in++)
		{
			int inn=HYPER[i].NEI[in];
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
			int k=HYPER[i].NEI[j];
			n_DgDq_x[i][k]=J*(p_Fi[0][0]*HYPER1[k*h_num+i].n0ij[0]+p_Fi[1][0]*HYPER1[k*h_num+i].n0ij[1]+p_Fi[2][0]*HYPER1[k*h_num+i].n0ij[2]);
			n_DgDq_y[i][k]=J*(p_Fi[0][1]*HYPER1[k*h_num+i].n0ij[0]+p_Fi[1][1]*HYPER1[k*h_num+i].n0ij[1]+p_Fi[2][1]*HYPER1[k*h_num+i].n0ij[2]);
			n_DgDq_z[i][k]=J*(p_Fi[0][2]*HYPER1[k*h_num+i].n0ij[0]+p_Fi[1][2]*HYPER1[k*h_num+i].n0ij[1]+p_Fi[2][2]*HYPER1[k*h_num+i].n0ij[2]);
		}
		n_DgDq_x[i][i]=J*(p_Fi[0][0]*HYPER1[i*h_num+i].n0ij[0]+p_Fi[1][0]*HYPER1[i*h_num+i].n0ij[1]+p_Fi[2][0]*HYPER1[i*h_num+i].n0ij[2]);
		n_DgDq_y[i][i]=J*(p_Fi[0][1]*HYPER1[i*h_num+i].n0ij[0]+p_Fi[1][1]*HYPER1[i*h_num+i].n0ij[1]+p_Fi[2][1]*HYPER1[i*h_num+i].n0ij[2]);
		n_DgDq_z[i][i]=J*(p_Fi[0][2]*HYPER1[i*h_num+i].n0ij[0]+p_Fi[1][2]*HYPER1[i*h_num+i].n0ij[1]+p_Fi[2][2]*HYPER1[i*h_num+i].n0ij[2]);
	}

	////DfDx�̍X�V
/*	for(int i=0;i<h_num;i++)
	{
		for(int j=0;j<h_num;j++)
		{
			double DFDlambda=0;
			for(int k=0;k<h_num;k++)	DFDlambda+=n_DgDq_x[i][k]*HYPER1[j*h_num+k].DgDq[A_X]+n_DgDq_y[i][k]*HYPER1[j*h_num+k].DgDq[A_Y]+n_DgDq_z[i][k]*HYPER1[j*h_num+k].DgDq[A_Z];
			if(DFDlambda!=0)	DfDx[i*h_num+j]=-Dt*Dt*0.5/mi*DFDlambda;//-DFDlambda;//
		}
	}//*/
	
	
	for(int k=0;k<h_num;k++)
	{
		double DFDlambda=0;
		int Nk=HYPER[k].N;
		for(int i=0;i<Nk;i++)
		{
			int in=HYPER[k].NEI[i];
			for(int j=0;j<Nk;j++)
			{
				int jn=HYPER[k].NEI[j];
				DfDx[in*h_num+jn]-=Dt*Dt*0.5/mi*(n_DgDq_x[in][k]*HYPER1[jn*h_num+k].DgDq[A_X]+n_DgDq_y[in][k]*HYPER1[jn*h_num+k].DgDq[A_Y]+n_DgDq_z[in][k]*HYPER1[jn*h_num+k].DgDq[A_Z]);
			}
			DfDx[k*h_num+in]-=Dt*Dt*0.5/mi*(n_DgDq_x[k][k]*HYPER1[in*h_num+k].DgDq[A_X]+n_DgDq_y[k][k]*HYPER1[in*h_num+k].DgDq[A_Y]+n_DgDq_z[k][k]*HYPER1[in*h_num+k].DgDq[A_Z]);
			DfDx[in*h_num+k]-=Dt*Dt*0.5/mi*(n_DgDq_x[in][k]*HYPER1[k*h_num+k].DgDq[A_X]+n_DgDq_y[in][k]*HYPER1[k*h_num+k].DgDq[A_Y]+n_DgDq_z[in][k]*HYPER1[k*h_num+k].DgDq[A_Z]);
		}
		DfDx[k*h_num+k]-=Dt*Dt*0.5/mi*(n_DgDq_x[k][k]*HYPER1[k*h_num+k].DgDq[A_X]+n_DgDq_y[k][k]*HYPER1[k*h_num+k].DgDq[A_Y]+n_DgDq_z[k][k]*HYPER1[k*h_num+k].DgDq[A_Z]);		
	}//*/

	////�o��
//	if(count%200==0 && count>CON.get_nr()/2)
//	if(t==1||t%CON.get_interval()==0)	if(count%200==0||count==1)	output_newton_data1(fx,DfDx,n_rx,n_ry,n_rz,h_num,count,t);
//	if(t==1||t%(10*CON.get_interval())==0)	if(count%200==0||count==1)	output_newton_data1(fx,DfDx,n_rx,n_ry,n_rz,h_num,count,t);

	ofstream t_loge("time_log_newton.dat", ios::app);
	clock_t t4=clock();
	t_loge<<"step="<<t<<", count="<<count<<", time="<<1000*(t4-t3)/CLOCKS_PER_SEC<<"[e-3sec]"<<endl;
	t_loge.close();


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

void calc_half_p(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,bool repetation,double **F)
{
	if(repetation==0)	cout<<"���̉^���ʁ��ʒu���W�v�Z";
	else	cout<<"�^���ʌv�Z";

	int h_num=HYPER.size();
	double Dt=CON.get_dt();
	double le=CON.get_distancebp();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	int flag_vis=CON.get_flag_vis();
	bool flag_FEM=CON.get_FEM_flag();
	int flag_G=CON.get_flag_G();
	double density=CON.get_hyper_density();
	int model_num=CON.get_model_number();

	for(int i=0;i<h_num;i++)
	{
		double p_half_p[DIMENSION]={0,0,0};
		int Ni=HYPER[i].N;
		for(int j=0;j<Ni;j++)
		{		
			int k=HYPER[i].NEI[j];
			p_half_p[A_X]+=(HYPER[k].stress[0][0]-HYPER[k].lambda)*HYPER1[k*h_num+i].DgDq[0]+HYPER[k].stress[0][1]*HYPER1[k*h_num+i].DgDq[1]+HYPER[k].stress[0][2]*HYPER1[k*h_num+i].DgDq[2];
			p_half_p[A_Y]+=HYPER[k].stress[1][0]*HYPER1[k*h_num+i].DgDq[0]+(HYPER[k].stress[1][1]-HYPER[k].lambda)*HYPER1[k*h_num+i].DgDq[1]+HYPER[k].stress[1][2]*HYPER1[k*h_num+i].DgDq[2];
			p_half_p[A_Z]+=HYPER[k].stress[2][0]*HYPER1[k*h_num+i].DgDq[0]+HYPER[k].stress[2][1]*HYPER1[k*h_num+i].DgDq[1]+(HYPER[k].stress[2][2]-HYPER[k].lambda)*HYPER1[k*h_num+i].DgDq[2];
		}//j�Ɋւ���for���̏I���	
		p_half_p[A_X]+=(HYPER[i].stress[0][0]-HYPER[i].lambda)*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[0][1]*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[0][2]*HYPER1[i*h_num+i].DgDq[2];
		p_half_p[A_Y]+=HYPER[i].stress[1][0]*HYPER1[i*h_num+i].DgDq[0]+(HYPER[i].stress[1][1]-HYPER[i].lambda)*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[1][2]*HYPER1[i*h_num+i].DgDq[2];
		p_half_p[A_Z]+=HYPER[i].stress[2][0]*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[2][1]*HYPER1[i*h_num+i].DgDq[1]+(HYPER[i].stress[2][2]-HYPER[i].lambda)*HYPER1[i*h_num+i].DgDq[2];

		//�d�͍�
		if(flag_G==ON)	p_half_p[A_Z]-=9.8*mi;
		//�S����
		if(flag_vis==ON)
		{
			p_half_p[A_X]+=HYPER[i].vis_force[A_X];
			p_half_p[A_Y]+=HYPER[i].vis_force[A_Y];
			p_half_p[A_Z]+=HYPER[i].vis_force[A_Z];
		}
		//���͍�
		if(flag_FEM==ON)
		{
			p_half_p[A_X]+=F[A_X][i]*mi;//density;
			p_half_p[A_Y]+=F[A_Y][i]*mi;//density;
			p_half_p[A_Z]+=F[A_Z][i]*mi;//density;
		/*	if(i==0){
			cout<<" 1="<<mi;
			cout<<" 2="<<F[A_X][i]*mi;
			cout<<" 3="<<F[A_X][i]/density<<endl;
			cout<<" 4="<<HYPER[i].vis_force[A_X];
			cout<<" 5="<<9.8*mi<<endl;}
		}*/
		}
		if(repetation==0)
		{

			//half_p�̍X�V
			HYPER[i].half_p[A_X]=HYPER[i].p[A_X]+Dt*0.5*p_half_p[A_X];
			HYPER[i].half_p[A_Y]=HYPER[i].p[A_Y]+Dt*0.5*p_half_p[A_Y];
			HYPER[i].half_p[A_Z]=HYPER[i].p[A_Z]+Dt*0.5*p_half_p[A_Z];//
			//�ʒu���W�̍X�V
			if(model_num==30||model_num==23)
			{
				if(PART[i].q0[A_Z]!=0)
				{
					PART[i].r[A_X]+=Dt*HYPER[i].half_p[A_X]/mi;
					PART[i].r[A_Y]+=Dt*HYPER[i].half_p[A_Y]/mi;
					PART[i].r[A_Z]+=Dt*HYPER[i].half_p[A_Z]/mi;					
				}
				else
				{
					PART[i].r[A_X]=PART[i].q0[A_X];
					PART[i].r[A_Y]=PART[i].q0[A_Y];
					PART[i].r[A_Z]=PART[i].q0[A_Z];
				}
			}
			else
			{
				PART[i].r[A_X]+=Dt*HYPER[i].half_p[A_X]/mi;
				PART[i].r[A_Y]+=Dt*HYPER[i].half_p[A_Y]/mi;
				PART[i].r[A_Z]+=Dt*HYPER[i].half_p[A_Z]/mi;
			}
		}
		else
		{
			//�^���ʂ̍X�V
			HYPER[i].p[A_X]=HYPER[i].half_p[A_X]+Dt*0.5*p_half_p[A_X];
			HYPER[i].p[A_Y]=HYPER[i].half_p[A_Y]+Dt*0.5*p_half_p[A_Y];
			HYPER[i].p[A_Z]=HYPER[i].half_p[A_Z]+Dt*0.5*p_half_p[A_Z];////
			//���x�̍X�V
			PART[i].u[A_X]=HYPER[i].half_p[A_X]/mi;
			PART[i].u[A_Y]=HYPER[i].half_p[A_Y]/mi;
			PART[i].u[A_Z]=HYPER[i].half_p[A_Z]/mi;
			//�p�^���ʂ̍X�V
			HYPER[i].ang_p[A_X]=PART[i].r[A_Y]*HYPER[i].p[A_Z]-PART[i].r[A_Z]*HYPER[i].p[A_Y];
			HYPER[i].ang_p[A_Y]=PART[i].r[A_Z]*HYPER[i].p[A_X]-PART[i].r[A_X]*HYPER[i].p[A_Z];
			HYPER[i].ang_p[A_Z]=PART[i].r[A_X]*HYPER[i].p[A_Y]-PART[i].r[A_Y]*HYPER[i].p[A_X];
		}
	}//i�Ɋւ���for���̏I���
	
	cout<<"----------OK"<<endl;
}

void calc_F(vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1)
{
	cout<<"Fi�v�Z";
	////Fi�̍X�V
	int h_num=HYPER.size();

	double **p_Fi=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	p_Fi[D]=new double[DIMENSION];

	for(int i=0;i<h_num;i++)
	{
		double fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
		//Fi�̌v�Z

		int Ni=HYPER[i].N;	

		for(int in=0;in<Ni;in++)
		{
			int inn=HYPER[i].NEI[in];
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
	
	for(int D=0;D<DIMENSION;D++)	delete[]	p_Fi[D];
	delete[]	p_Fi;

	//calculation of DgDq
	for(int i=0;i<h_num;i++)
	{
		int Ni=HYPER[i].N;
		for(int j=0;j<Ni;j++)
		{			
			int k=HYPER[i].NEI[j];
			HYPER1[k*h_num+i].DgDq[A_X]=HYPER[k].J*(HYPER[k].t_inverse_Fi[A_X][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].t_inverse_Fi[A_X][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].t_inverse_Fi[A_X][2]*HYPER1[i*h_num+k].n0ij[2]);
			HYPER1[k*h_num+i].DgDq[A_Y]=HYPER[k].J*(HYPER[k].t_inverse_Fi[A_Y][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].t_inverse_Fi[A_Y][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].t_inverse_Fi[A_Y][2]*HYPER1[i*h_num+k].n0ij[2]);
			HYPER1[k*h_num+i].DgDq[A_Z]=HYPER[k].J*(HYPER[k].t_inverse_Fi[A_Z][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].t_inverse_Fi[A_Z][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].t_inverse_Fi[A_Z][2]*HYPER1[i*h_num+k].n0ij[2]);
			//cout<<"i"<<i<<"j"<<k<<"	"<<HYPER1[k*h_num+i].DgDq[A_X]<<","<<HYPER1[k*h_num+i].DgDq[A_Y]<<","<<HYPER1[k*h_num+i].DgDq[A_Z]<<endl;
		}
		HYPER1[i*h_num+i].DgDq[A_X]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_X][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_X][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_X][2]*HYPER1[i*h_num+i].n0ij[2]);
		HYPER1[i*h_num+i].DgDq[A_Y]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_Y][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_Y][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_Y][2]*HYPER1[i*h_num+i].n0ij[2]);
		HYPER1[i*h_num+i].DgDq[A_Z]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_Z][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_Z][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_Z][2]*HYPER1[i*h_num+i].n0ij[2]);
		//cout<<"i"<<i<<"j"<<i<<"	"<<HYPER1[i*h_num+i].DgDq[A_X]<<","<<HYPER1[i*h_num+i].DgDq[A_Y]<<","<<HYPER1[i*h_num+i].DgDq[A_Z]<<endl;
	}
	cout<<"----------OK"<<endl;

}

void calc_stress(mpsconfig &CON,vector<hyperelastic> &HYPER)
{
	cout<<"���͌v�Z";
	int h_num=HYPER.size();

	double **d_Fi=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	d_Fi[D]=new double [DIMENSION];

	double c10=CON.get_c10();
	double c01=CON.get_c01();
	
	double b[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double bb[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

	for(int j=0;j<h_num;j++)
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

	for(int D=0;D>DIMENSION;D++)	delete[]	d_Fi[D];
	delete[]	d_Fi;

	cout<<"----------OK"<<endl;
}

void calc_differential_p(mpsconfig &CON,vector<mpselastic>PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,double **F)
{
	cout<<"�^���ʔ����l�v�Z";

	int h_num=HYPER.size();
	int flag_vis=CON.get_flag_vis();
	int flag_G=CON.get_flag_G();
	int flag_FEM=CON.get_FEM_flag();
	double Dt=CON.get_dt();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	double density=CON.get_hyper_density();
	for(int i=0;i<h_num;i++)
	{
		double p_differential_p[DIMENSION]={0,0,0};
		int Ni=HYPER[i].N;
		for(int j=0;j<Ni;j++)
		{				
			int k=HYPER[i].NEI[j];
			p_differential_p[A_X]+=HYPER[k].stress[A_X][0]*HYPER1[k*h_num+i].DgDq[0]+HYPER[k].stress[A_X][1]*HYPER1[k*h_num+i].DgDq[1]+HYPER[k].stress[A_X][2]*HYPER1[k*h_num+i].DgDq[2];
			p_differential_p[A_Y]+=HYPER[k].stress[A_Y][0]*HYPER1[k*h_num+i].DgDq[0]+HYPER[k].stress[A_Y][1]*HYPER1[k*h_num+i].DgDq[1]+HYPER[k].stress[A_Y][2]*HYPER1[k*h_num+i].DgDq[2];	
			p_differential_p[A_Z]+=HYPER[k].stress[A_Z][0]*HYPER1[k*h_num+i].DgDq[0]+HYPER[k].stress[A_Z][1]*HYPER1[k*h_num+i].DgDq[1]+HYPER[k].stress[A_Z][2]*HYPER1[k*h_num+i].DgDq[2];	
		}
		HYPER[i].differential_p[A_X]=HYPER[i].half_p[A_X]+Dt*0.5*(p_differential_p[A_X]+HYPER[i].stress[A_X][0]*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[A_X][1]*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[A_X][2]*HYPER1[i*h_num+i].DgDq[2]);
		HYPER[i].differential_p[A_Y]=HYPER[i].half_p[A_X]+Dt*0.5*(p_differential_p[A_X]+HYPER[i].stress[A_Y][0]*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[A_Y][1]*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[A_Y][2]*HYPER1[i*h_num+i].DgDq[2]);	
		HYPER[i].differential_p[A_Z]=HYPER[i].half_p[A_Z]+Dt*0.5*(p_differential_p[A_Z]+HYPER[i].stress[A_Z][0]*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[A_Z][1]*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[A_Z][2]*HYPER1[i*h_num+i].DgDq[2]);////p_differential_p[A_Z];////
		
		//�d�͉e��
		if(flag_G==ON)	HYPER[i].differential_p[A_Z]-=Dt*0.5*9.8*mi;
		if(flag_vis==ON)
		{
			HYPER[i].differential_p[A_X]+=Dt*0.5*HYPER[i].vis_force[A_X];
			HYPER[i].differential_p[A_Y]+=Dt*0.5*HYPER[i].vis_force[A_Y];
			HYPER[i].differential_p[A_Z]+=Dt*0.5*HYPER[i].vis_force[A_Z];
		}
		if(flag_FEM==ON)
		{
			HYPER[i].differential_p[A_X]+=Dt*0.5*F[A_X][i]*mi;//density;
			HYPER[i].differential_p[A_Y]+=Dt*0.5*F[A_Y][i]*mi;//density;
			HYPER[i].differential_p[A_Z]+=Dt*0.5*F[A_Z][i]*mi;//density;
		}
	}
	cout<<"----------OK"<<endl;

}

void renew_lambda(mpsconfig &CON,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int t)
{
	clock_t t3=clock();

	cout<<"Lambda�v�Z";

	int h_num=HYPER.size();
	double le=CON.get_distancebp();
	double Dt=CON.get_dt();
	double V=get_volume(&CON);
	double mk=V*CON.get_hyper_density();

	double *N_Left=new double[h_num*h_num];
	double *N_Right=new double[h_num];
	for(int i=0;i<h_num;i++)
	{
		N_Right[i]=0;
		for(int j=0;j<h_num;j++)	N_Left[j*h_num+i]=0;
	}

	
/*	for(int i=0;i<h_num;i++)
	{
		for(int j=0;j<h_num;j++)
		{
			double N_left=0;
			for(int k=0;k<h_num;k++)
			{
				N_left+=HYPER1[i*h_num+k].DgDq[0]*HYPER1[j*h_num+k].DgDq[0]+HYPER1[i*h_num+k].DgDq[1]*HYPER1[j*h_num+k].DgDq[1]+HYPER1[i*h_num+k].DgDq[2]*HYPER1[j*h_num+k].DgDq[2];
				//cout<<"i"<<i<<"j"<<j<<"k"<<k<<"	"<<HYPER1[i*h_num+k].DgDq[0]<<"	"<<HYPER1[j*h_num+k].DgDq[0]<<"	"<<HYPER1[i*h_num+k].DgDq[1]<<"	"<<HYPER1[j*h_num+k].DgDq[1]<<"	"<<HYPER1[i*h_num+k].DgDq[2]<<"	"<<HYPER1[j*h_num+k].DgDq[2]<<endl;
				//cout<<"i"<<i<<"j"<<j<<"k"<<k<<"	"<<HYPER1[i*h_num+k].DgDq[0]*HYPER1[j*h_num+k].DgDq[0]+HYPER1[i*h_num+k].DgDq[1]*HYPER1[j*h_num+k].DgDq[1]+HYPER1[i*h_num+k].DgDq[2]*HYPER1[j*h_num+k].DgDq[2]<<endl;
			}
			N_Left[i*h_num+j]=Dt*0.5*N_left;//Dt/2/mk*N_left;
		}//j�Ɋւ���for���̏I���
	}//i�Ɋւ���for���̏I���
	for(int i=0;i<h_num;i++)
	{
		double N_right=0;
		for(int j=0;j<h_num;j++)
		{
			N_right+=HYPER[j].differential_p[0]*HYPER1[i*h_num+j].DgDq[0]+HYPER[j].differential_p[1]*HYPER1[i*h_num+j].DgDq[1]+HYPER[j].differential_p[2]*HYPER1[i*h_num+j].DgDq[2];
			//cout<<"i"<<i<<"j"<<j<<"	"<<HYPER[j].differential_p[0]*HYPER1[i*h_num+j].DgDq[0]+HYPER[j].differential_p[1]*HYPER1[i*h_num+j].DgDq[1]+HYPER[j].differential_p[2]*HYPER1[i*h_num+j].DgDq[2]<<endl;
		}
		N_Right[i]=N_right;//1/mk*N_right;	
	}//*/
	
	for(int k=0;k<h_num;k++)
	{
		int Nk=HYPER[k].N;
		double N_right=0;
		for(int i=0;i<Nk;i++)
		{
			int in=HYPER[k].NEI[i];
			for(int j=0;j<Nk;j++)
			{
				int jn=HYPER[k].NEI[j];
				N_Left[in*h_num+jn]+=Dt*0.5*(HYPER1[in*h_num+k].DgDq[0]*HYPER1[jn*h_num+k].DgDq[0]+HYPER1[in*h_num+k].DgDq[1]*HYPER1[jn*h_num+k].DgDq[1]+HYPER1[in*h_num+k].DgDq[2]*HYPER1[jn*h_num+k].DgDq[2]);
			}
			N_Left[in*h_num+k]+=Dt*0.5*(HYPER1[in*h_num+k].DgDq[0]*HYPER1[k*h_num+k].DgDq[0]+HYPER1[in*h_num+k].DgDq[1]*HYPER1[k*h_num+k].DgDq[1]+HYPER1[in*h_num+k].DgDq[2]*HYPER1[k*h_num+k].DgDq[2]);
			N_Left[k*h_num+in]+=Dt*0.5*(HYPER1[k*h_num+k].DgDq[0]*HYPER1[in*h_num+k].DgDq[0]+HYPER1[k*h_num+k].DgDq[1]*HYPER1[in*h_num+k].DgDq[1]+HYPER1[k*h_num+k].DgDq[2]*HYPER1[in*h_num+k].DgDq[2]);
			N_right+=HYPER[in].differential_p[0]*HYPER1[k*h_num+in].DgDq[0]+HYPER[in].differential_p[1]*HYPER1[k*h_num+in].DgDq[1]+HYPER[in].differential_p[2]*HYPER1[k*h_num+in].DgDq[2];
		}//j�Ɋւ���for���̏I���
		N_Right[k]=N_right+HYPER[k].differential_p[0]*HYPER1[k*h_num+k].DgDq[0]+HYPER[k].differential_p[1]*HYPER1[k*h_num+k].DgDq[1]+HYPER[k].differential_p[2]*HYPER1[k*h_num+k].DgDq[2];//1/mk*N_right;	
		N_Left[k*h_num+k]+=Dt*0.5*(HYPER1[k*h_num+k].DgDq[0]*HYPER1[k*h_num+k].DgDq[0]+HYPER1[k*h_num+k].DgDq[1]*HYPER1[k*h_num+k].DgDq[1]+HYPER1[k*h_num+k].DgDq[2]*HYPER1[k*h_num+k].DgDq[2]);
	}//i�Ɋւ���for���̏I���*/

/*
	for(int i=0;i<h_num;i++)	cout<<"N_Right["<<i<<"]="<<N_Right[i]<<endl;
	cout<<endl;
	for(int i=0;i<h_num;i++)
	{
		for(int j=0;j<h_num;j++)	cout<<"N_Left["<<i<<"]["<<j<<"]="<<N_Left[i*h_num+j]<<endl;
		cout<<endl;
	}*/

	//lambda�����߂�

	/*ofstream fl("N_Left.csv");
	ofstream fr("N_Reft.csv");

	for(int i=0;i<h_num;i++)
	{
		for(int j=0;j<h_num;j++)
		{
			fl<<N_Left[i*h_num+j]<<",";
		}
		fl<<endl;
		fr<<N_Right[i]<<endl;
	}
	fl.close();
	fr.close();//*/
		
	gauss(N_Left,N_Right,h_num);

	for(int i=0;i<h_num;i++)	HYPER[i].lambda=N_Right[i];
//	for(int i=0;i<h_num;i++)	cout<<"lambda["<<i<<"]="<<HYPER[i].lambda<<endl;

	ofstream t_loge("time_log_gauss.dat", ios::app);
	clock_t t4=clock();
	t_loge<<"step="<<t<<", time="<<1000*(t4-t3)/CLOCKS_PER_SEC<<"[e-3sec]"<<endl;
	t_loge.close();

	delete [] N_Left;
	delete [] N_Right;

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

	for(int D=0;D<N;D++)	for(int D2=0;D2<N;D2++)	if(fabs(a[D][D2])<DBL_EPSILON)	a[D][D2]=0;

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
	int count=0;
	double max=0,temp=0,sum=0,dum=0;
	*d=1.0;
	vector<float> buf(N);

	for(int i=0;i<N;i++)	buf[i]=0;

	for(int i=0;i<N;i++)
	{
		max=0.0;
		for(int j=0;j<N;j++)
		{
			temp=fabs(a[i][j]);
			if(temp>max)	max=temp;
		}
		if(max==0.0&&count==0)
		{
			count++;
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
void momentum_movie_AVS(mpsconfig &CON,int t,vector<mpselastic> PART,vector<hyperelastic> HYPER,double **F)
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
	fp<<"10 0"<<endl;//�ߓ_�̏��ʂ�8�ŁA�v�f�̏��ʂ�0�Ƃ������ƁB
	fp<<"10 1 1 1 1 1 1 1 1 1 1"<<endl;	//���̍s�̏ڍׂ̓w���v���Q��
	//fp<<"8 1 1 1 1 1 1 1 1"<<endl;	//���̍s�̏ڍׂ̓w���v���Q��
	fp<<"p_x,"<<endl;
	fp<<"p_y,"<<endl;
	fp<<"p_z,"<<endl;
	fp<<"lambda,"<<endl;
	fp<<"F_x,"<<endl;
	fp<<"F_y,"<<endl;
	fp<<"F_z,"<<endl;
	fp<<"v_x,"<<endl;
	fp<<"v_y,"<<endl;
	fp<<"v_z,"<<endl;
	//fp<<"P,N/m^2"<<endl;
	//fp<<"value1,??"<<endl;

	//�e�ߓ_�̏��l����
	for(int i=0;i<h_num;i++)
	{
		fp<<i<<" "<<HYPER[i].p[A_X]<<" "<<HYPER[i].p[A_Y]<<" "<<HYPER[i].p[A_Z]<<" "<<HYPER[i].lambda<<" "<<F[A_X][i]<<" "<<F[A_Y][i]<<" "<<F[A_Z][i]<<" "<<HYPER[i].vis_force[A_X]<<" "<<HYPER[i].vis_force[A_Y]<<" "<<HYPER[i].vis_force[A_Z]<<endl;
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

void calculation_vec_norm(vector<mpselastic> PART, vector<hyperelastic> &HYPER,int hyper_number,int particle_number,int t)
{
	int h_num=hyper_number;
	int p_num=particle_number;

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

void contact_judge(mpsconfig &CON, vector<mpselastic> PART,vector<hyperelastic> &HYPER,double max_h,int t)
{
	int h_num=HYPER.size();
	double r=CON.get_h_dis();
	double Dt=CON.get_interval();

	for(int i=0;i<h_num;i++)
	{
		double qik_z=PART[i].r[A_Z]-max_h;
		if(qik_z<r)	HYPER[i].p[A_Z]-=HYPER[i].p[A_Z]/qik_z*Dt;
	}	
}


void output_hyper_data(vector<mpselastic> PART,vector<hyperelastic> HYPER,vector<hyperelastic2> HYPER1,int t)
{
	int h_num=HYPER.size();

	/*if(t==1)
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
	}*/

/*	stringstream ss_dgdq;
	ss_dgdq<<"./DgDq/DgDq"<<t<<".csv";
	ofstream dg(ss_dgdq.str());
	ofstream d_p("d_P.csv", ios::app);
	ofstream h_p("h_P.csv", ios::app);*/
	ofstream lam("lambda.csv", ios::app);
/*	ofstream stress("stress.csv", ios::app);*/
	ofstream p("P.csv", ios::app);
/*	ofstream J("J.csv", ios::app);
	ofstream ti_Fi("ti_Fi.csv", ios::app);
	ofstream Fi("Fi.csv", ios::app);*/

	if(t==1)
	{
		p<<"t"<<",";		
/*		d_p<<"t"<<",";
		h_p<<"t"<<",";*/
		lam<<"t"<<",";
/*		stress<<"t"<<",";
		ti_Fi<<"t"<<",";
		Fi<<"t"<<",";
		J<<"t"<<",";*/

		for(int i=0;i<h_num;i++)
		{
			p<<i<<","<<","<<",";
			lam<<i<<",";
/*			d_p<<i<<","<<","<<",";
			h_p<<i<<","<<","<<",";
			stress<<i<<","<<","<<",";
			ti_Fi<<i<<","<<","<<",";
			Fi<<i<<","<<","<<",";
			J<<i<<",";*/
		}
		p<<endl;
		lam<<endl;
/*		d_p<<endl;
		h_p<<endl;
		stress<<endl;
		ti_Fi<<endl;
		Fi<<endl;
		J<<endl;*/
	}

	p<<t<<",";
	lam<<t<<",";
/*	d_p<<t<<",";
	h_p<<t<<",";
	stress<<t<<",";
	ti_Fi<<t<<",";
	Fi<<t<<",";
	J<<t<<",";*/
	for(int i=0;i<h_num;i++)
	{
		
		p<<HYPER[i].p[A_X]<<","<<HYPER[i].p[A_Y]<<","<<HYPER[i].p[A_Z]<<",";
		lam<<HYPER[i].lambda<<",";
/*		d_p<<HYPER[i].differential_p[A_X]<<","<<HYPER[i].differential_p[A_Y]<<","<<HYPER[i].differential_p[A_Z]<<","<<endl;
		h_p<<HYPER[i].half_p[A_X]<<","<<HYPER[i].half_p[A_Y]<<","<<HYPER[i].half_p[A_Z]<<","<<endl;

		stress<<HYPER[i].stress[A_X][A_X]<<","<<HYPER[i].stress[A_X][A_Y]<<","<<HYPER[i].stress[A_X][A_Z]<<",";
		ti_Fi<<HYPER[i].t_inverse_Fi[A_X][A_X]<<","<<HYPER[i].t_inverse_Fi[A_X][A_Y]<<","<<HYPER[i].t_inverse_Fi[A_X][A_Z]<<",";	
		Fi<<HYPER[i].Fi[0][0]<<","<<HYPER[i].Fi[0][1]<<","<<HYPER[i].Fi[0][2]<<",";
		J<<HYPER[i].J<<",";
		for(int j=0;j<h_num;j++)	dg<<i<<","<<j<<","<<HYPER1[j*h_num+i].DgDq[A_X]<<","<<HYPER1[j*h_num+i].DgDq[A_Y]<<","<<HYPER1[j*h_num+i].DgDq[A_Z]<<","<<endl;
		*/
	}
/*	stress<<endl<<",";
	ti_Fi<<endl<<",";
	Fi<<endl<<",";
	for(int i=0;i<h_num;i++)
	{
		stress<<HYPER[i].stress[A_Y][A_X]<<","<<HYPER[i].stress[A_Y][A_Y]<<","<<HYPER[i].stress[A_Y][A_Z]<<",";
		ti_Fi<<HYPER[i].t_inverse_Fi[A_Y][A_X]<<","<<HYPER[i].t_inverse_Fi[A_Y][A_Y]<<","<<HYPER[i].t_inverse_Fi[A_Y][A_Z]<<",";
		Fi<<HYPER[i].Fi[1][0]<<","<<HYPER[i].Fi[1][1]<<","<<HYPER[i].Fi[1][2]<<",";
	}
	stress<<endl<<",";
	ti_Fi<<endl<<",";
	Fi<<endl<<",";
	for(int i=0;i<h_num;i++)
	{
		stress<<HYPER[i].stress[A_Z][A_X]<<","<<HYPER[i].stress[A_Z][A_Y]<<","<<HYPER[i].stress[A_Z][A_Z]<<",";
		ti_Fi<<HYPER[i].t_inverse_Fi[A_Z][A_X]<<","<<HYPER[i].t_inverse_Fi[A_Z][A_Y]<<","<<HYPER[i].t_inverse_Fi[A_Z][A_Z]<<",";
		Fi<<HYPER[i].Fi[2][0]<<","<<HYPER[i].Fi[2][1]<<","<<HYPER[i].Fi[2][2]<<",";
	}*/

	p<<endl;
	lam<<endl;
/*	d_p<<endl;
	h_p<<endl;
	stress<<endl;
	ti_Fi<<endl;
	Fi<<endl;
	J<<endl;

	dg.close();
	stress.close();
	ti_Fi.close();
	Fi.close();
	d_p.close();
	h_p.close();
	J.close();*/
	p.close();
	lam.close();
}

void output_newton_data1(double *fx, double *DfDx, double *n_rx, double *n_ry, double *n_rz,int hyper_number,int count, int t)
{
	int h_num=hyper_number;

	stringstream ss_r;
	ss_r<<"./Newton_raphson/position"<<t<<".xlsx";
	stringstream ss_Df;
	ss_Df<<"./Newton_raphson/DfDx "<<t<<".xlsx";
	stringstream ss_fx;
	ss_fx<<"./Newton_raphson/fx"<<t<<".xlsx";
		
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
		r<<n_rx[i]<<","<<n_ry[i]<<","<<n_rz[i]<<",";
		sfx<<fx[i]<<",";
		for(int j=0;j<h_num;j++) Df<<DfDx[i*h_num+j]<<",";
		Df<<endl;

	}
	r<<endl;
	sfx<<endl;

	r.close();
	sfx.close();

	Df.close();
}

void output_newton_data2(double E, double *XX, int hyper_number, int count, int t)
{
	int h_num=hyper_number;
	stringstream ss_E;
	ss_E<<"./Newton_raphson/E"<<t<<".xlsx";
	
	stringstream ss_lam;
	ss_lam<<"./Newton_raphson/lambda"<<t<<".xlsx";
		
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

void output_energy(mpsconfig CON, vector<mpselastic> PART, vector<hyperelastic> HYPER,int t)
{
	cout<<"�e���|�e���V�����v�Z";
	int h_num=HYPER.size();
	int p_num=PART.size();
	double d_Fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double c10=CON.get_c10();
	double c01=CON.get_c01();
	vector<double>	W;
	for(int i=0;i<h_num;i++)	W.push_back(0);
	double dC[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double dC2[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

	for(int i=0;i<h_num;i++)
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
		for(int i=0;i<h_num;i++)
		{
			e<<i<<",";
			e_T<<i<<",";
			e_g<<i<<",";
			e_W<<i<<",";
			e_lam<<i<<",";
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

	for(int i=0;i<h_num;i++)
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

void calc_gravity(mpsconfig CON,vector<hyperelastic> &HYPER,int hyper_number)
{
	int h_num=hyper_number;
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	//�d�͍~��
	double Dt=CON.get_dt();
	for(int i=0;i<h_num;i++)	HYPER[i].p[A_Z]-=9.8*Dt;
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
		NEI[i]=0;
	}
	N=0;
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

