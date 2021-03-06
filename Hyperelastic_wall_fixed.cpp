#include "stdafx.h"	

//諸定義
//void hyper_initial(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int *Nw);
void calc_half_p(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,bool repetation,double **F,int Nw);
void renew_lambda(mpsconfig &CON,vector<mpselastic>PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int t);
void calc_differential_p(mpsconfig &CON,vector<mpselastic>PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,double **F);
void calc_stress(mpsconfig &CON,vector<hyperelastic> &HYPER);
//void calc_pi(mpsconfig &CON,vector<hyperelastic> &HYPER);
void calc_constant(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1);
void newton_raphson(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t,double **F,int Nw);
void calc_F(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1);
void calc_newton_function(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,double *lambda,double *XX_old,double *fx,double *DfDx,int hyper_number,int count,int t,double **F,double *hp_x,double *hp_y,double *hp_z,int Nw);
void previous_strage(vector<mpselastic>PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int h_num);
void calc_W(mpsconfig &CON,vector<hyperelastic> &HYPER,int h_num);

void momentum_movie_AVS(mpsconfig &CON,int t,vector<mpselastic> PART,vector<hyperelastic> HYPER,double **F);
void output_hyper_data(vector<mpselastic> PART,vector<hyperelastic> HYPER,vector<hyperelastic2> HYPER1,int t,double V);
void output_newton_data(vector<mpselastic> PART,vector<hyperelastic> HYPER,vector<hyperelastic2> HYPER1,double *fx,double *DfDx,int hyper_number,int count, int t);
void output_newton_data2(double E, double *XX, double *fx, int hyper_number, int count, int t);
void output_energy(mpsconfig CON, vector<mpselastic> PART, vector<hyperelastic> HYPER,int t);

double calc_det(double **M,int N);
void calc_inverse_matrix_for_NR(int N, double *a);
void BiCGStab2_method(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X);
void iccg2(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X);
void CG3D(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X);
void GaussSeidelvh(double *A, int pn, double *b,double ep);

void hyper_initial(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int *Nw)
{
	int h_num=HYPER.size();
	double Dt=CON.get_dt();
	double mh=CON.get_hyper_density()*get_volume(&CON);
	double ms=CON.get_silicone_density()*get_volume(&CON);
	double V=get_volume(&CON);
	int Nw_p=0;

	for(int i=0;i<h_num;i++)
	{
		for(int D=0;D<DIMENSION;D++)
		{
			PART[i].q0[D]=0.;
			PART[i].q0[D]=PART[i].r[D];
		}
		if(PART[i].q0[A_Z]>-1.e-20&&PART[i].q0[A_Z]<1.e-20)
		{
			HYPER[i].fw=1;
			Nw_p++;
		}
	}
	*Nw=Nw_p;
	//cout<<"Nw="<<Nw_p<<endl;
	//ofstream fm("h_min.csv");
	//for(int i=0;i<h_num;i++)
	//{
	//	fm<<i<<","<<HYPER[i].fw<<","<<PART[i].toFEM<<endl;
	//}
	//fm.close();


	calc_constant(CON,PART,HYPER,HYPER1);
	calc_stress(CON,HYPER);
	calc_W(CON,HYPER,h_num);
	for(int i=0;i<h_num;i++)
	{
		HYPER[i].W0=HYPER[i].W;
		for(int D=0;D<DIMENSION;D++)	HYPER[i].p0[D]=HYPER[i].p[D];
	}
}

void calc_constant(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1)
{
	//	cout<<"初期値計算";

	double le=CON.get_distancebp();
	double r=CON.get_h_dis();
	double Dt=CON.get_dt();
	double V=get_volume(&CON);	//考慮が必要かもしれない
	double mh=V*CON.get_hyper_density();
	double ms=V*CON.get_silicone_density();
	int h_num=HYPER.size();
	int model=CON.get_model_number();
	//	cout<<"V"<<V<<endl;

	ofstream fc("constant.dat");
	fc<<"le"<<","<<le<<endl;
	fc<<"r"<<","<<r<<endl;
	fc<<"Dt"<<","<<Dt<<endl;

	//垂直降下
	//for(int i=0;i<h_num;i++)	HYPER[i].p[A_Z]=-1.*mi;

	//曲げねじり
	/*if(model==21)
	{		
	int t=30,b=2;
	double min=PART[0].q0[A_Z];
	double max=PART[0].q0[A_Z];

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
	}//*/

	//////曲げ
	//if(model==21)
	//{
	//	int b=10;
	//	double max=0,min=0;

	//	for(int i=0;i<h_num;i++)
	//	{
	//		if(max<PART[i].q0[A_Z])	max=PART[i].q0[A_Z];
	//		if(min>PART[i].q0[A_Z])	min=PART[i].q0[A_Z];
	//	}
	//	double H=max-min+le;
	//	cout<<H;
	//	//double H=1.8;
	//	for(int i=0;i<h_num;i++)	
	//	{
	//		double Z=PART[i].q0[A_Z];
	//		double part_p=(Z/H)*2;
	//		if(PART[i].toFEM==1)
	//		{
	//			HYPER[i].p[A_X]=mh*b*(3*part_p*part_p-1);
	//		}
	//		else
	//		{
	//			HYPER[i].p[A_X]=ms*b*(3*part_p*part_p-1);
	//		}

	//	}
	//}//*/

	//////ねじり
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


	//回転
	//if(model==22)
	//{
	//	for(int i=0;i<h_num;i++)
	//	{
	//		PART[i].p[A_X]=mi*0.4*(PART[i].q0[A_Z]-PART[i].q0[A_Y]);
	//		PART[i].p[A_Y]=mi*0.4*(PART[i].q0[A_X]-PART[i].q0[A_Z]);
	//		PART[i].p[A_Z]=mi*0.4*(PART[i].q0[A_Y]-PART[i].q0[A_X]);
	//	}
	//}

	////角運動量計算
	for(int i=0;i<h_num;i++)
	{
		HYPER[i].ang_p[A_X]=PART[i].r[A_Y]*HYPER[i].p[A_Z]-PART[i].r[A_Z]*HYPER[i].p[A_Y];
		HYPER[i].ang_p[A_Y]=PART[i].r[A_Z]*HYPER[i].p[A_X]-PART[i].r[A_X]*HYPER[i].p[A_Z];
		HYPER[i].ang_p[A_Z]=PART[i].r[A_X]*HYPER[i].p[A_Y]-PART[i].r[A_Y]*HYPER[i].p[A_X];
	}


	////近傍粒子の記憶とaiin,wiin,Aiの計算
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
				HYPER[i].pnd+=wiin;
				N++;
			}
			else	wiin=0;

			HYPER1[i*h_num+j].wiin=wiin;
		}
		HYPER[i].N=N;
		HYPER[i].pnd0=N;
	}


	////Ai, Fi関連の計算
	double **p_Ai=new double *[DIMENSION];
	double **p_Fi=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)
	{
		p_Ai[D]=new double [DIMENSION];
		p_Fi[D]=new double[DIMENSION];
	}
	for(int i=0;i<h_num;i++)
	{
		//Aiの計算
		double Ai[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
		Ai[0][0]=0;	Ai[0][1]=0;	Ai[0][2]=0;
		Ai[1][0]=0;	Ai[1][1]=0;	Ai[1][2]=0;
		Ai[2][0]=0;	Ai[2][1]=0;	Ai[2][2]=0;

		int Ni=HYPER[i].N;
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

		//inverse_Ai,t_inverse_Aiの計算
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

		//Fiの計算
		p_Fi[0][0]=Ai[0][0]*p_Ai[0][0]+Ai[0][1]*p_Ai[1][0]+Ai[0][2]*p_Ai[2][0];	p_Fi[0][1]=Ai[0][0]*p_Ai[0][1]+Ai[0][1]*p_Ai[1][1]+Ai[0][2]*p_Ai[2][1];	p_Fi[0][2]=Ai[0][0]*p_Ai[0][2]+Ai[0][1]*p_Ai[1][2]+Ai[0][2]*p_Ai[2][2];
		p_Fi[1][0]=Ai[1][0]*p_Ai[0][0]+Ai[1][1]*p_Ai[1][0]+Ai[1][2]*p_Ai[2][0];	p_Fi[1][1]=Ai[1][0]*p_Ai[0][1]+Ai[1][1]*p_Ai[1][1]+Ai[1][2]*p_Ai[2][1];	p_Fi[1][2]=Ai[1][0]*p_Ai[0][2]+Ai[1][1]*p_Ai[1][2]+Ai[1][2]*p_Ai[2][2];
		p_Fi[2][0]=Ai[2][0]*p_Ai[0][0]+Ai[2][1]*p_Ai[1][0]+Ai[2][2]*p_Ai[2][0];	p_Fi[2][1]=Ai[2][0]*p_Ai[0][1]+Ai[2][1]*p_Ai[1][1]+Ai[2][2]*p_Ai[2][1];	p_Fi[2][2]=Ai[2][0]*p_Ai[0][2]+Ai[2][1]*p_Ai[1][2]+Ai[2][2]*p_Ai[2][2];		

		HYPER[i].Fi[0][0]=p_Fi[0][0];	HYPER[i].Fi[0][1]=p_Fi[0][1];	HYPER[i].Fi[0][2]=p_Fi[0][2];	
		HYPER[i].Fi[1][0]=p_Fi[1][0];	HYPER[i].Fi[1][1]=p_Fi[1][1];	HYPER[i].Fi[1][2]=p_Fi[1][2];	
		HYPER[i].Fi[2][0]=p_Fi[2][0];	HYPER[i].Fi[2][1]=p_Fi[2][1];	HYPER[i].Fi[2][2]=p_Fi[2][2];	
		//Jの計算
		double J=calc_det3(p_Fi);
		HYPER[i].J=J;
		//		cout<<"J["<<i<<"]="<<HYPER[i].J<<endl;

		//t_inverse_Fiの計算
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


	////n0ijの計算
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

	////DgDqの計算
	for(int i=0;i<h_num;i++)
	{
		int Ni=HYPER[i].N;
		for(int j=0;j<Ni;j++)
		{
			int k=HYPER[i].NEI[j];
			HYPER1[k*h_num+i].DgDq[A_X]=HYPER[k].J*(HYPER[k].t_inverse_Fi[A_X][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].t_inverse_Fi[A_X][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].t_inverse_Fi[A_X][2]*HYPER1[i*h_num+k].n0ij[2]);
			HYPER1[k*h_num+i].DgDq[A_Y]=HYPER[k].J*(HYPER[k].t_inverse_Fi[A_Y][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].t_inverse_Fi[A_Y][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].t_inverse_Fi[A_Y][2]*HYPER1[i*h_num+k].n0ij[2]);
			HYPER1[k*h_num+i].DgDq[A_Z]=HYPER[k].J*(HYPER[k].t_inverse_Fi[A_Z][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].t_inverse_Fi[A_Z][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].t_inverse_Fi[A_Z][2]*HYPER1[i*h_num+k].n0ij[2]);		
			//			cout<<"DgDq["<<k<<","<<i<<"]="<<HYPER1[k*h_num+i].DgDq[A_X]<<","<<HYPER1[k*h_num+i].DgDq[A_Y]<<","<<HYPER1[k*h_num+i].DgDq[A_Z]<<endl;
		}
		HYPER1[i*h_num+i].DgDq[A_X]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_X][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_X][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_X][2]*HYPER1[i*h_num+i].n0ij[2]);
		HYPER1[i*h_num+i].DgDq[A_Y]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_Y][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_Y][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_Y][2]*HYPER1[i*h_num+i].n0ij[2]);
		HYPER1[i*h_num+i].DgDq[A_Z]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_Z][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_Z][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_Z][2]*HYPER1[i*h_num+i].n0ij[2]);
		//		cout<<"DgDq["<<i<<","<<i<<"]="<<HYPER1[i*h_num+i].DgDq[A_X]<<","<<HYPER1[i*h_num+i].DgDq[A_Y]<<","<<HYPER1[i*h_num+i].DgDq[A_Z]<<endl;
	}	
	//	cout<<"----------OK"<<endl;


}



////主計算
void calc_hyper(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t,double **F,int Nw)
{	
	ofstream time("time_log.dat",ios::app);
	clock_t	start_t=clock();
	time<<start_t*CLOCKS_PER_SEC<<"	";

	int h_num=HYPER.size();
	double Dt=CON.get_dt();
	double mh=CON.get_hyper_density()*get_volume(&CON);
	double ms=CON.get_silicone_density()*get_volume(&CON);
	double V=get_volume(&CON);

	cout<<"h_num="<<h_num<<endl;
	cout<<"Hypercalculation starts."<<endl;


	///////////外力計算
	if(CON.get_flag_vis()==ON)	calc_vis_f(CON,PART,HYPER,t);


	///////////データ出力
	if(t==1 || t%CON.get_interval()==0)
	{
		output_hyper_data(PART,HYPER,HYPER1,t,V);
		momentum_movie_AVS(CON,t,PART,HYPER,F);
		output_energy(CON,PART,HYPER,t);
	}

	///////////前ステップ計算結果保管
	previous_strage(PART,HYPER,HYPER1,h_num);


	/////////////////通常計算プロセス
	newton_raphson(CON,PART,HYPER,HYPER1,t,F,Nw);
	calc_half_p(CON,PART,HYPER,HYPER1,0,F,Nw);
	calc_F(CON,PART,HYPER,HYPER1);
	calc_stress(CON,HYPER);
	//calc_pi(CON,HYPER);
	calc_differential_p(CON,PART,HYPER,HYPER1,F);
	renew_lambda(CON,PART,HYPER,HYPER1,t);
	calc_half_p(CON,PART,HYPER,HYPER1,1,F,Nw);


	/////////////壁計算プロセス
	calc_W(CON,HYPER,h_num);
	//////calc_HYPER_QP_gh(CON,PART,HYPER,HYPER1,t,F);

	//calc_HYPER_QP_g(CON,PART,HYPER,HYPER1,t,F,Nw);

	//vector<double> NEIw;
	//double nG[DIMENSION]={0,0,1};
	//double aG[DIMENSION]={0,0,0};
	//double hi=0.;
	//for(int i=0;i<h_num;i++)
	//{
	//	hi=-1.*PART[i].r[A_Z];
	//	if(hi>0)
	//	{
	//		NEIw.push_back(i);
	//		HYPER[i].fw=1;
	//	}
	//	else
	//	{
	//		HYPER[i].fw=0;
	//	}
	//}
	//int Nw=NEIw.size();
	//if(Nw>0)
	//{
	//	cout<<"接触 ";
	//	for(int i=0;i<Nw;i++)	cout<<NEIw[i]<<", ";
	//	cout<<endl;
	//	calc_HYPER_QP_gh(CON,PART,HYPER,HYPER1,t,F,NEIw);
	//}
	//NEIw.clear();
	//calc_HYPER_QP_gh(CON,PART,HYPER,HYPER1,t,F,NEIw);
	cout<<"Hypercalculation ends."<<endl;

	clock_t end_t=clock();
	time<<end_t*CLOCKS_PER_SEC<<"	";
	time.close();
}


void newton_raphson(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t,double **F,int Nw)
{
	/////fx(N*N行列の各成分)は1次元配列で格納、(i,j)成分なら[j*N+i]で参照
	/////DfDx(N*N行列の各成分)は1次元配列で格納、(i,j)成分なら[j*N+i]で参照

#ifdef _OPENMP
#pragma omp parallel for
	omp_set_num_threads(8);
#endif
	{
		int calc_type=1;//ニュートラフソンの反復方法 0:偏微分項の逆行列をそのまま求める　1:線形方程式を利用

		//OPENMPがオンであれば書き込まれる。threadsはpcで扱える最大スレッド数
		//最大スレッド数で計算し続けると過負荷でCPUが非常に熱くなるため、並列化数を指定している
		//なお、最大スレッドが12のときは8〜10ぐらいが目安

		//pn=2;//test,とりあえず2元でとけるかどうか確認 
		//////////////////　f1(x1,x2) = x1^2 + x2^2 -5 = 0 f2(x1,x2) = x1^2/9+ x2^2 -1 = 0  http://homepage1.nifty.com/gfk/excel_newton_ren.htm

		int h_num=HYPER.size();
		double *fx=new double [h_num];//関数値。
		double *DfDx=new double [h_num*h_num];//関数の偏微分値。
		double *XX=new double [h_num];//現在の解。	
		double *XX_old=new double [h_num];//1ステップ前の解。
		double *hp_x=new double [h_num];
		double *hp_y=new double [h_num];
		double *hp_z=new double [h_num];
		double ep=1e-5;//収束判定
		double ep_gs=1e-30;//収束判定
		double E=1;//現在の誤差
		int count=0;//反復回数
		double d;
		double V=get_volume(&CON);
		double mh=V*CON.get_hyper_density();
		double ms=V*CON.get_silicone_density();
		double Dt=CON.get_interval();
		double sum=0;
		double E_old=0;
		int dec_flag=ON;

#ifdef _OPENMP
#pragma omp for
#endif
		for(int i=0;i<h_num;i++)
		{
			XX[i]=1;
			XX_old[i]=0;
			fx[i]=0;
			hp_x[i]=0.;
			hp_y[i]=0.;
			hp_z[i]=0.;
			for(int j=0;j<h_num;j++)	DfDx[i*h_num+j]=0;
		}

		//hp計算
		double G=9.8;
		double p_half_p[DIMENSION]={0,0,0};

#ifdef _OPENMP
#pragma omp for
#endif
		for(int i=0;i<h_num;i++)
		{
			//half_pの計算
			p_half_p[A_X]=0.;
			p_half_p[A_Y]=0.;
			p_half_p[A_Z]=0.;

			int Ni=HYPER[i].N;
			for(int j=0;j<Ni;j++)
			{	
				int k=HYPER[i].NEI[j];
				p_half_p[A_X]+=HYPER[k].stress_n[0][0]*HYPER1[k*h_num+i].DgDq_n[0]+HYPER[k].stress_n[0][1]*HYPER1[k*h_num+i].DgDq_n[1]+HYPER[k].stress_n[0][2]*HYPER1[k*h_num+i].DgDq_n[2];
				p_half_p[A_Y]+=HYPER[k].stress_n[1][0]*HYPER1[k*h_num+i].DgDq_n[0]+HYPER[k].stress_n[1][1]*HYPER1[k*h_num+i].DgDq_n[1]+HYPER[k].stress_n[1][2]*HYPER1[k*h_num+i].DgDq_n[2];
				p_half_p[A_Z]+=HYPER[k].stress_n[2][0]*HYPER1[k*h_num+i].DgDq_n[0]+HYPER[k].stress_n[2][1]*HYPER1[k*h_num+i].DgDq_n[1]+HYPER[k].stress_n[2][2]*HYPER1[k*h_num+i].DgDq_n[2];
			}//jに関するfor文の終わり
			p_half_p[A_X]+=HYPER[i].stress_n[0][0]*HYPER1[i*h_num+i].DgDq_n[0]+HYPER[i].stress_n[0][1]*HYPER1[i*h_num+i].DgDq_n[1]+HYPER[i].stress_n[0][2]*HYPER1[i*h_num+i].DgDq_n[2];
			p_half_p[A_Y]+=HYPER[i].stress_n[1][0]*HYPER1[i*h_num+i].DgDq_n[0]+HYPER[i].stress_n[1][1]*HYPER1[i*h_num+i].DgDq_n[1]+HYPER[i].stress_n[1][2]*HYPER1[i*h_num+i].DgDq_n[2];
			p_half_p[A_Z]+=HYPER[i].stress_n[2][0]*HYPER1[i*h_num+i].DgDq_n[0]+HYPER[i].stress_n[2][1]*HYPER1[i*h_num+i].DgDq_n[1]+HYPER[i].stress_n[2][2]*HYPER1[i*h_num+i].DgDq_n[2];

			//重力の影響
			p_half_p[A_Z]-=G*mh;
			//粘性項の影響
			p_half_p[A_X]+=HYPER[i].vis_force[A_X];
			p_half_p[A_Y]+=HYPER[i].vis_force[A_Y];
			p_half_p[A_Z]+=HYPER[i].vis_force[A_Z];
			//磁場の考慮
			//p_half_p[A_X]+=F[A_X][i]*V;//density;
			//p_half_p[A_Y]+=F[A_Y][i]*V;//density;
			//p_half_p[A_Z]+=F[A_Z][i]*V;//density;

			//仮の運動量計算
			hp_x[i]=HYPER[i].p_n[A_X]+Dt*0.5*p_half_p[A_X];
			hp_y[i]=HYPER[i].p_n[A_Y]+Dt*0.5*p_half_p[A_Y];
			hp_z[i]=HYPER[i].p_n[A_Z]+Dt*0.5*p_half_p[A_Z];
		}


		ofstream time("time_log.dat",ios::app);
		clock_t start_t=clock();
		time<<start_t*CLOCKS_PER_SEC<<"	";

		//	for(int i=0; i<N; i++) XX[i]=1;///初期値を与える。とりあえず1で
		cout<<"NR法開始";
		while(E>ep)
		{
			count++;

			//#ifdef _OPENMP
			//#pragma omp for
			//#endif


			//		if(count==1)	for(int i=0;i<h_num;i++)	for(int j=0;j<h_num;j++)	for(int D=0;D<DIMENSION;D++)	HYPER1[i*h_num+j].newton_DgDq[D]=HYPER1[i*N+j].DgDq[D];

			calc_newton_function(CON,PART,HYPER,HYPER1,XX,XX_old,fx,DfDx,h_num,count,t,F,hp_x,hp_y,hp_z,Nw);


			/*		//現在の関数値を求める
			if(count==1) cout<<fx[0]<<" "<<fx[1]<<endl;
			//現在の偏微分値を求める
			//calc_DfDx(XX)////現在の偏微分値を求める。超弾性体ならば、calc_DgDq()などで求められるはず
			DfDx[0*N+0]=2*XX[0];
			DfDx[0*N+1]=2*XX[1];
			DfDx[1*N+0]=2*XX[0]/9;
			DfDx[1*N+1]=2*XX[1];
			if(count==1) cout<<DfDx[0]<<" "<<DfDx[1]<<" "<<DfDx[2]<<" "<<DfDx[3]<<endl;*/

			/////値の更新
			//if(calc_type==0)//逆行列を利用 逆行列が求まりさえすれば速いはず
			//{
			//	calc_inverse_matrix_for_NR(h_num,DfDx);

			//	for(int i=0; i<h_num; i++) 
			//	{
			//		d=0; //変化量
			//		for(int j=0; j<h_num; j++)	d+=DfDx[i*h_num+j]*fx[j];
			//		XX[i]-=d;
			//	}
			//}
			//else if(calc_type==1)//逆行列を用いない、安定するはずだが、遅くなるはず
			{	
				/*
				int *b_ind=new int [h_num*h_num];
				double *b_val=new double [h_num*h_num];
				int *ptr=new int [h_num+1];
				int all_ind_num=0;

				ptr[0]=0;
				for(int i=0;i<h_num;i++)
				{
				ptr[i+1]=0;
				for(int j=0;j<h_num;j++)
				{
				b_ind[i*h_num+j]=0;
				b_val[i*h_num+j]=0;
				}
				}


				for(int i=0;i<h_num;i++)
				{
				int ind_num=0;
				for(int j=0;j<h_num;j++)
				{
				if(DfDx[i*h_num+j]!=0)
				{
				b_ind[all_ind_num+ind_num]=j;	
				b_val[all_ind_num+ind_num]=DfDx[i*h_num+j];
				ind_num++;
				ptr[i+1]=ind_num+ptr[i];
				}			
				}
				//cout<<"ptr"<<i+1<<" "<<ptr[i+1]<<endl;
				all_ind_num+=ind_num;			
				}

				//cout<<all_ind_num<<endl;

				int *ind=new int [all_ind_num];
				double *val=new double [all_ind_num];
				for(int i=0;i<all_ind_num;i++)
				{
				ind[i]=b_ind[i];
				val[i]=b_val[i];
				}
				delete[] b_ind;
				delete[] b_val;

				//for(int i=0;i<h_num;i++)	cout<<"ptr"<<i<<"	"<<ptr[i]<<endl;
				//CG3D(&CON,val,ind,ptr,h_num,fx,all_ind_num,XX);
				///iccg2(&CON,val,ind,ptr,h_num,fx,all_ind_num,XX);
				BiCGStab2_method(&CON,val,ind,ptr,h_num,fx,all_ind_num,XX);
				//for(int i=0;i<h_num;i++)	cout<<"lambda"<<i<<"="<<XX[i]<<endl;

				delete[] ind;
				delete[] val;
				delete[] ptr;*/



				gauss(DfDx,fx,h_num);
				//GaussSeidelvh(DfDx,h_num,fx,ep_gs);
				sum=0;

#ifdef _OPENMP
#pragma omp for
#endif
				for(int i=0;i<h_num;i++)
				{
					XX[i]-=fx[i];//*0.5*mi/(Dt*Dt)*V*fx[i];//*/
					sum+=fabs(fx[i]);
				}

			}

			//誤差の評価
			E_old=E;
			E=sum;

			//sum=0;
			//for(int i=0; i<h_num; i++)
			//{
			//	sum+=fabs(XX[i]-XX_old[i]);
			//	//sum+=fabs(XX[i]);	//絶対誤差で評価
			//}
			//E=sum;


			if(count==1 || count%100000==0)
			{		
				/*
				cout<<"XX_old["<<i<<"]-d=X["<<i<<"]	";
				cout<<XX_old[i]<<" - ";
				cout<<d<<" = ";
				cout<<XX[i]<<endl;*/

				cout<<"E"<<count<<"="<<E<<endl;
				//output_newton_data2(E,XX,fx,h_num,count,t);

			}
			//if(count>CON.get_nr())	break;
			//if(count>CON.get_nr()&&E<1.e-1)	break;
			else if(dec_flag==ON)	if(E_old-E<0)	break;	
		}
		//cout<<"反復回数	"<<count<<" E="<<E<<endl;
		ofstream fs("Newton_E.csv", ios::app);
		fs<<count<<","<<E<<endl;
		fs.close();

		//	end=clock();
		//	newton_t=(end-start)/CLOCKS_PER_SEC;

		cout<<"反復完了";

#ifdef _OPENMP
#pragma omp for
#endif
		for(int i=0;i<h_num;i++) HYPER[i].lambda=XX[i];

		//	for(int i=0;i<h_num;i++)	cout<<"lambda["<<i<<"]="<<HYPER[i].lambda<<endl;
		delete[]	fx;
		delete[]	DfDx;
		delete[]	XX;
		delete[]	XX_old;
		delete[]	hp_x;
		delete[]	hp_y;
		delete[]	hp_z;

		cout<<"---------- OK"<<endl;
	}
}
void calc_newton_function(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,double *lambda,double *XX_old,double *fx,double *DfDx,int hyper_number,int count,int t,double **F,double *hp_x,double *hp_y,double *hp_z, int Nw)
{


#ifdef _OPENMP
#pragma omp parallel for
	printf("OpenMP : On, threads = %d\n", omp_get_max_threads());
	omp_set_num_threads(8);
#endif
	{
		clock_t t3=clock();

		int h_num=hyper_number;
		int flag_vis=CON.get_flag_vis();
		bool flag_FEM=CON.get_FEM_flag();
		int flag_G=CON.get_flag_G();
		int model_num=CON.get_model_number();
		double Dt=CON.get_dt();
		double V=get_volume(&CON);
		double mh=V*CON.get_hyper_density();
		double ms=V*CON.get_silicone_density();
		double density_h=CON.get_hyper_density();
		double density_s=CON.get_silicone_density();
		double mi=V*CON.get_hyper_density();
		double G=9.8;

		double **p_Fi=new double *[DIMENSION];
		for(int D=0;D<DIMENSION;D++)	p_Fi[D]=new double [DIMENSION];

		////位置座標の更新	
		double p_half_p[DIMENSION]={0,0,0};


#ifdef _OPENMP
#pragma omp for
#endif
		for(int i=0; i<Nw; i++)
		{
			fx[i]=0.;
			XX_old[i]=lambda[i];	//解を記憶

			//half_pの計算
			p_half_p[A_X]=0.;
			p_half_p[A_Y]=0.;
			p_half_p[A_Z]=0.;

			int Ni=HYPER[i].N;
			for(int j=0;j<Ni;j++)
			{	
				int k=HYPER[i].NEI[j];
				p_half_p[A_X]+=-lambda[k]*HYPER1[k*h_num+i].DgDq_n[0];
				p_half_p[A_Y]+=-lambda[k]*HYPER1[k*h_num+i].DgDq_n[1];
				p_half_p[A_Z]+=-lambda[k]*HYPER1[k*h_num+i].DgDq_n[2];
			}//jに関するfor文の終わり
			p_half_p[A_X]+=-lambda[i]*HYPER1[i*h_num+i].DgDq_n[0];
			p_half_p[A_Y]+=-lambda[i]*HYPER1[i*h_num+i].DgDq_n[1];
			p_half_p[A_Z]+=-lambda[i]*HYPER1[i*h_num+i].DgDq_n[2];

			//位置座標の計算
			PART[i].r[A_X]=HYPER[i].q_n[A_X]+Dt*(hp_x[i]+Dt*0.5*p_half_p[A_X])/mh;
			PART[i].r[A_Y]=HYPER[i].q_n[A_Y]+Dt*(hp_y[i]+Dt*0.5*p_half_p[A_Y])/mh;
			PART[i].r[A_Z]=HYPER[i].q_n[A_Z]+Dt*(hp_z[i]+Dt*0.5*p_half_p[A_Z])/mh;
			//fx[i]=(PART[i].r[A_X]-PART[i].q0[A_X])*(PART[i].r[A_X]-PART[i].q0[A_X])+(PART[i].r[A_Y]-PART[i].q0[A_Y])*(PART[i].r[A_Y]-PART[i].q0[A_Y])+(PART[i].r[A_Z]-PART[i].q0[A_Z])*(PART[i].r[A_Z]-PART[i].q0[A_Z]);
			fx[i]=PART[i].r[A_X]-PART[i].q0[A_X]+PART[i].r[A_Y]-PART[i].q0[A_Y]+PART[i].r[A_Z]-PART[i].q0[A_Z];
			//PART[i].r[A_X]=PART[i].q0[A_X];
			//PART[i].r[A_Y]=PART[i].q0[A_Y];
			PART[i].r[A_Z]=PART[i].q0[A_Z];

			//cout<<"fx"<<fx[i]<<endl;
		}

#ifdef _OPENMP
#pragma omp for
#endif
		for(int i=Nw;i<h_num;i++)
		{
			fx[i]=0.;
			XX_old[i]=lambda[i];	//解を記憶

			//half_pの計算
			p_half_p[A_X]=0.;
			p_half_p[A_Y]=0.;
			p_half_p[A_Z]=0.;

			int Ni=HYPER[i].N;
			for(int j=0;j<Ni;j++)
			{	
				int k=HYPER[i].NEI[j];
				p_half_p[A_X]+=-lambda[k]*HYPER1[k*h_num+i].DgDq_n[0];
				p_half_p[A_Y]+=-lambda[k]*HYPER1[k*h_num+i].DgDq_n[1];
				p_half_p[A_Z]+=-lambda[k]*HYPER1[k*h_num+i].DgDq_n[2];
			}//jに関するfor文の終わり
			p_half_p[A_X]+=-lambda[i]*HYPER1[i*h_num+i].DgDq_n[0];
			p_half_p[A_Y]+=-lambda[i]*HYPER1[i*h_num+i].DgDq_n[1];
			p_half_p[A_Z]+=-lambda[i]*HYPER1[i*h_num+i].DgDq_n[2];

			//位置座標の計算
			PART[i].r[A_X]=HYPER[i].q_n[A_X]+Dt*(hp_x[i]+Dt*0.5*p_half_p[A_X])/mh;
			PART[i].r[A_Y]=HYPER[i].q_n[A_Y]+Dt*(hp_y[i]+Dt*0.5*p_half_p[A_Y])/mh;
			PART[i].r[A_Z]=HYPER[i].q_n[A_Z]+Dt*(hp_z[i]+Dt*0.5*p_half_p[A_Z])/mh;
		}



		////DgDqとfxの更新
		int p_num=PART.size();
		double r=CON.get_h_dis();
		//
		//#ifdef _OPENMP
		//#pragma omp for
		//#endif
		//		for(int i=0;i<Nw;i++)
		//		{
		//			//Fiの計算
		//			int Ni=HYPER[i].N;
		//			double fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};	
		//
		//			double part_dgdq_x=0;
		//			double part_dgdq_y=0;
		//			double part_dgdq_z=0;
		//			double part_g=0;
		//			double pnd=0;
		//
		//			//		cout<<"Ni="<<Ni<<endl;
		//			for(int in=0;in<Ni;in++)
		//			{
		//				int inn=HYPER[i].NEI[in];
		//				double w=HYPER1[i*h_num+inn].wiin;
		//
		//				//			cout<<"i["<<i<<"]j["<<inn<<"]="<<n_rx[inn]-n_rx[i]<<", "<<n_ry[inn]-n_ry[i]<<", "<<n_rz[inn]-n_rz[i]<<endl;
		//				//			cout<<"rz_inn["<<inn<<"]="<<n_rz[inn]<<", rz_i["<<i<<"]="<<n_rz[i]<<endl;
		//				fi[0][0]+=w*(PART[inn].r[A_X]-PART[i].r[A_X])*HYPER1[i*h_num+inn].aiin[A_X];
		//				fi[0][1]+=w*(PART[inn].r[A_X]-PART[i].r[A_X])*HYPER1[i*h_num+inn].aiin[A_Y];
		//				fi[0][2]+=w*(PART[inn].r[A_X]-PART[i].r[A_X])*HYPER1[i*h_num+inn].aiin[A_Z];
		//				fi[1][0]+=w*(PART[inn].r[A_Y]-PART[i].r[A_Y])*HYPER1[i*h_num+inn].aiin[A_X];
		//				fi[1][1]+=w*(PART[inn].r[A_Y]-PART[i].r[A_Y])*HYPER1[i*h_num+inn].aiin[A_Y];
		//				fi[1][2]+=w*(PART[inn].r[A_Y]-PART[i].r[A_Y])*HYPER1[i*h_num+inn].aiin[A_Z];
		//				fi[2][0]+=w*(PART[inn].r[A_Z]-PART[i].r[A_Z])*HYPER1[i*h_num+inn].aiin[A_X];
		//				fi[2][1]+=w*(PART[inn].r[A_Z]-PART[i].r[A_Z])*HYPER1[i*h_num+inn].aiin[A_Y];
		//				fi[2][2]+=w*(PART[inn].r[A_Z]-PART[i].r[A_Z])*HYPER1[i*h_num+inn].aiin[A_Z];
		//
		//				/*			double dis=sqrt((n_rx[inn]-n_rx[i])*(n_rx[inn]-n_rx[i])+(n_ry[inn]-n_ry[i])*(n_ry[inn]-n_ry[i])+(n_rz[inn]-n_rz[i])*(n_rz[inn]-n_rz[i]));
		//				pnd+=kernel4(r,dis);*/
		//				/*			if(n_rz[inn]<0)
		//				{
		//				//double dis=sqrt((n_rx[inn]-n_rx[i])*(n_rx[inn]-n_rx[i])+(n_ry[inn]-n_ry[i])*(n_ry[inn]-n_ry[i])+(n_rz[inn]-n_rz[i])*(n_rz[inn]-n_rz[i]));
		//				double dis=sqrt(HYPER1[i*h_num+inn].aiin[A_X]*HYPER1[i*h_num+inn].aiin[A_X]+HYPER1[i*h_num+inn].aiin[A_Y]*HYPER1[i*h_num+inn].aiin[A_Y]+HYPER1[i*h_num+inn].aiin[A_Z]*HYPER1[i*h_num+inn].aiin[A_Z]);
		//
		//				double d_wij=-r/dis/dis;
		//				//				double wij=kernel4(r,dis);
		//				part_g+=w/HYPER[i].pnd;
		//
		//				n_DgDq_x[i][inn]+=V/HYPER[i].pnd*d_wij*HYPER1[i*h_num+inn].aiin[A_X]/dis;
		//				n_DgDq_y[i][inn]+=V/HYPER[i].pnd*d_wij*HYPER1[i*h_num+inn].aiin[A_Y]/dis;
		//				n_DgDq_z[i][inn]+=V/HYPER[i].pnd*d_wij*HYPER1[i*h_num+inn].aiin[A_Z]/dis;
		//
		//				n_DgDq_x[i][i]+=-V/HYPER[i].pnd*d_wij*HYPER1[i*h_num+inn].aiin[A_X]/dis;
		//				n_DgDq_y[i][i]+=-V/HYPER[i].pnd*d_wij*HYPER1[i*h_num+inn].aiin[A_Y]/dis;
		//				n_DgDq_z[i][i]+=-V/HYPER[i].pnd*d_wij*HYPER1[i*h_num+inn].aiin[A_Z]/dis;
		//				}*/
		//			}
		//			p_Fi[0][0]=fi[0][0]*HYPER[i].inverse_Ai[0][0]+fi[0][1]*HYPER[i].inverse_Ai[1][0]+fi[0][2]*HYPER[i].inverse_Ai[2][0];
		//			p_Fi[0][1]=fi[0][0]*HYPER[i].inverse_Ai[0][1]+fi[0][1]*HYPER[i].inverse_Ai[1][1]+fi[0][2]*HYPER[i].inverse_Ai[2][1];
		//			p_Fi[0][2]=fi[0][0]*HYPER[i].inverse_Ai[0][2]+fi[0][1]*HYPER[i].inverse_Ai[1][2]+fi[0][2]*HYPER[i].inverse_Ai[2][2];
		//			p_Fi[1][0]=fi[1][0]*HYPER[i].inverse_Ai[0][0]+fi[1][1]*HYPER[i].inverse_Ai[1][0]+fi[1][2]*HYPER[i].inverse_Ai[2][0];	
		//			p_Fi[1][1]=fi[1][0]*HYPER[i].inverse_Ai[0][1]+fi[1][1]*HYPER[i].inverse_Ai[1][1]+fi[1][2]*HYPER[i].inverse_Ai[2][1];
		//			p_Fi[1][2]=fi[1][0]*HYPER[i].inverse_Ai[0][2]+fi[1][1]*HYPER[i].inverse_Ai[1][2]+fi[1][2]*HYPER[i].inverse_Ai[2][2];
		//			p_Fi[2][0]=fi[2][0]*HYPER[i].inverse_Ai[0][0]+fi[2][1]*HYPER[i].inverse_Ai[1][0]+fi[2][2]*HYPER[i].inverse_Ai[2][0];
		//			p_Fi[2][1]=fi[2][0]*HYPER[i].inverse_Ai[0][1]+fi[2][1]*HYPER[i].inverse_Ai[1][1]+fi[2][2]*HYPER[i].inverse_Ai[2][1];
		//			p_Fi[2][2]=fi[2][0]*HYPER[i].inverse_Ai[0][2]+fi[2][1]*HYPER[i].inverse_Ai[1][2]+fi[2][2]*HYPER[i].inverse_Ai[2][2];
		//
		//			//		cout<<"F"<<i<<"="<<p_Fi[0][0]<<", "<<p_Fi[0][1]<<", "<<p_Fi[0][2]<<endl;
		//			//		cout<<p_Fi[1][0]<<", "<<p_Fi[1][1]<<", "<<p_Fi[1][2]<<endl;
		//			//		cout<<p_Fi[2][0]<<", "<<p_Fi[2][1]<<", "<<p_Fi[2][2]<<endl;
		//
		//			//Jの計算
		//			double Ji=calc_det3(p_Fi);
		//			HYPER[i].J=Ji;
		//			//		cout<<"J"<<i<<"="<<J<<endl;
		//
		//			/*		double part_dgdq_x=0;
		//			double part_dgdq_y=0;
		//			double part_dgdq_z=0;
		//			double part_g=0;
		//			for(int k=h_num;k<p_num;k++)
		//			{
		//			double dis=sqrt((PART[k].r[A_X]-n_rx[i])*(PART[k].r[A_X]-n_rx[i])+(PART[k].r[A_Y]-n_ry[i])*(PART[k].r[A_Y]-n_ry[i])+(PART[k].r[A_Z]-n_rz[i])*(PART[k].r[A_Z]-n_rz[i]));
		//			if(dis<r)
		//			{
		//			double wij=kernel4(r,dis);
		//			double d_wij=-r/dis/dis;
		//			part_g+=wij/HYPER[i].pnd;
		//			part_dgdq_x+=-V/HYPER[i].pnd*d_wij*(PART[k].r[A_X]-n_rx[i])/dis;
		//			part_dgdq_y+=-V/HYPER[i].pnd*d_wij*(PART[k].r[A_Y]-n_ry[i])/dis;
		//			part_dgdq_z+=-V/HYPER[i].pnd*d_wij*(PART[k].r[A_Z]-n_rz[i])/dis;
		//			}
		//			}//*/
		//			//fxの計算
		//			fx[i]=V*(1-Ji);//1-J;//
		//			//		cout<<"fx"<<i<<"="<<V*(1-J)<<endl;
		//
		//			//t_inverse_Fiの計算
		//			inverse(p_Fi,DIMENSION);
		//			HYPER[i].t_inverse_Fi[0][A_X]=p_Fi[A_X][0];	HYPER[i].t_inverse_Fi[0][A_Y]=p_Fi[A_Y][0];	HYPER[i].t_inverse_Fi[0][A_Z]=p_Fi[A_Z][0];
		//			HYPER[i].t_inverse_Fi[1][A_X]=p_Fi[A_X][1];	HYPER[i].t_inverse_Fi[1][A_Y]=p_Fi[A_Y][1];	HYPER[i].t_inverse_Fi[1][A_Z]=p_Fi[A_Z][1];
		//			HYPER[i].t_inverse_Fi[2][A_X]=p_Fi[A_X][2];	HYPER[i].t_inverse_Fi[2][A_Y]=p_Fi[A_Y][2];	HYPER[i].t_inverse_Fi[2][A_Z]=p_Fi[A_Z][2];
		//
		//			////DgDqの計算
		//			//for(int j=0;j<Ni;j++)
		//			//{
		//			//	int k=HYPER[i].NEI[j];
		//			//	HYPER1[i*h_num+k].DgDq[A_X]=J*(p_Fi[0][0]*HYPER1[k*h_num+i].n0ij[0]+p_Fi[1][0]*HYPER1[k*h_num+i].n0ij[1]+p_Fi[2][0]*HYPER1[k*h_num+i].n0ij[2]);
		//			//	HYPER1[i*h_num+k].DgDq[A_Y]=J*(p_Fi[0][1]*HYPER1[k*h_num+i].n0ij[0]+p_Fi[1][1]*HYPER1[k*h_num+i].n0ij[1]+p_Fi[2][1]*HYPER1[k*h_num+i].n0ij[2]);
		//			//	HYPER1[i*h_num+k].DgDq[A_Z]=J*(p_Fi[0][2]*HYPER1[k*h_num+i].n0ij[0]+p_Fi[1][2]*HYPER1[k*h_num+i].n0ij[1]+p_Fi[2][2]*HYPER1[k*h_num+i].n0ij[2]);
		//			//	//		cout<<"n_DgDq["<<i<<"*p_num+"<<k<<"]={"<<n_DgDq_x[i][k]<<","<<n_DgDq_y[i][k]<<","<<n_DgDq_z[i][k]<<"}"<<endl;
		//			//}
		//			//HYPER1[i*h_num+i].DgDq[A_X]=J*(p_Fi[0][0]*HYPER1[i*h_num+i].n0ij[0]+p_Fi[1][0]*HYPER1[i*h_num+i].n0ij[1]+p_Fi[2][0]*HYPER1[i*h_num+i].n0ij[2]);
		//			//HYPER1[i*h_num+i].DgDq[A_Y]=J*(p_Fi[0][1]*HYPER1[i*h_num+i].n0ij[0]+p_Fi[1][1]*HYPER1[i*h_num+i].n0ij[1]+p_Fi[2][1]*HYPER1[i*h_num+i].n0ij[2]);
		//			//HYPER1[i*h_num+i].DgDq[A_Z]=J*(p_Fi[0][2]*HYPER1[i*h_num+i].n0ij[0]+p_Fi[1][2]*HYPER1[i*h_num+i].n0ij[1]+p_Fi[2][2]*HYPER1[i*h_num+i].n0ij[2]);
		//
		//			//		cout<<"n_DgDq["<<i<<"*p_num+"<<i<<"]={"<<n_DgDq_x[i][i]<<","<<n_DgDq_y[i][i]<<","<<n_DgDq_z[i][i]<<"}"<<endl;
		//			//		cout<<"fx["<<i<<"]="<<fx[i]<<endl;
		//		}

#ifdef _OPENMP
#pragma omp for
#endif
		for(int i=0;i<h_num;i++)
		{
			//Fiの計算
			int Ni=HYPER[i].N;
			double fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};	

			double part_dgdq_x=0;
			double part_dgdq_y=0;
			double part_dgdq_z=0;
			double part_g=0;
			double pnd=0;

			//		cout<<"Ni="<<Ni<<endl;
			for(int in=0;in<Ni;in++)
			{
				int inn=HYPER[i].NEI[in];
				double w=HYPER1[i*h_num+inn].wiin;

				//			cout<<"i["<<i<<"]j["<<inn<<"]="<<n_rx[inn]-n_rx[i]<<", "<<n_ry[inn]-n_ry[i]<<", "<<n_rz[inn]-n_rz[i]<<endl;
				//			cout<<"rz_inn["<<inn<<"]="<<n_rz[inn]<<", rz_i["<<i<<"]="<<n_rz[i]<<endl;
				fi[0][0]+=w*(PART[inn].r[A_X]-PART[i].r[A_X])*HYPER1[i*h_num+inn].aiin[A_X];
				fi[0][1]+=w*(PART[inn].r[A_X]-PART[i].r[A_X])*HYPER1[i*h_num+inn].aiin[A_Y];
				fi[0][2]+=w*(PART[inn].r[A_X]-PART[i].r[A_X])*HYPER1[i*h_num+inn].aiin[A_Z];
				fi[1][0]+=w*(PART[inn].r[A_Y]-PART[i].r[A_Y])*HYPER1[i*h_num+inn].aiin[A_X];
				fi[1][1]+=w*(PART[inn].r[A_Y]-PART[i].r[A_Y])*HYPER1[i*h_num+inn].aiin[A_Y];
				fi[1][2]+=w*(PART[inn].r[A_Y]-PART[i].r[A_Y])*HYPER1[i*h_num+inn].aiin[A_Z];
				fi[2][0]+=w*(PART[inn].r[A_Z]-PART[i].r[A_Z])*HYPER1[i*h_num+inn].aiin[A_X];
				fi[2][1]+=w*(PART[inn].r[A_Z]-PART[i].r[A_Z])*HYPER1[i*h_num+inn].aiin[A_Y];
				fi[2][2]+=w*(PART[inn].r[A_Z]-PART[i].r[A_Z])*HYPER1[i*h_num+inn].aiin[A_Z];
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

			//Jの計算
			double Ji=calc_det3(p_Fi);
			HYPER[i].J=Ji;
			//fxの計算
			fx[i]+=V*(1-Ji);//1-J;//
			//		cout<<"fx"<<i<<"="<<V*(1-J)<<endl;

			//t_inverse_Fiの計算
			inverse(p_Fi,DIMENSION);
			HYPER[i].t_inverse_Fi[0][A_X]=p_Fi[A_X][0];	HYPER[i].t_inverse_Fi[0][A_Y]=p_Fi[A_Y][0];	HYPER[i].t_inverse_Fi[0][A_Z]=p_Fi[A_Z][0];
			HYPER[i].t_inverse_Fi[1][A_X]=p_Fi[A_X][1];	HYPER[i].t_inverse_Fi[1][A_Y]=p_Fi[A_Y][1];	HYPER[i].t_inverse_Fi[1][A_Z]=p_Fi[A_Z][1];
			HYPER[i].t_inverse_Fi[2][A_X]=p_Fi[A_X][2];	HYPER[i].t_inverse_Fi[2][A_Y]=p_Fi[A_Y][2];	HYPER[i].t_inverse_Fi[2][A_Z]=p_Fi[A_Z][2];
		}

		////DfDxの更新
		double DgDq_kk[DIMENSION]={0,0,0};
		double DgDq_ik[DIMENSION]={0,0,0};
		double DgDq_jk[DIMENSION]={0,0,0};

		double Dgkk[DIMENSION]={0,0,0};
		double Dgik[DIMENSION]={0,0,0};

		double Jk=0., Jin=0.;
		double Fk[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
		double Fin[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

#ifdef _OPENMP
#pragma omp for
#endif
		for(int k=0;k<h_num;k++)
		{			
			int Nk=HYPER[k].N;
			DgDq_kk[A_X]=HYPER1[k*h_num+k].DgDq_n[A_X];
			DgDq_kk[A_Y]=HYPER1[k*h_num+k].DgDq_n[A_Y];
			DgDq_kk[A_Z]=HYPER1[k*h_num+k].DgDq_n[A_Z];

			Jk=HYPER[k].J;
			Fk[0][A_X]=HYPER[k].t_inverse_Fi[0][A_X];	Fk[0][A_Y]=HYPER[k].t_inverse_Fi[0][A_Y];	Fk[0][A_Z]=HYPER[k].t_inverse_Fi[0][A_Z];
			Fk[1][A_X]=HYPER[k].t_inverse_Fi[1][A_X];	Fk[1][A_Y]=HYPER[k].t_inverse_Fi[1][A_Y];	Fk[1][A_Z]=HYPER[k].t_inverse_Fi[1][A_Z];
			Fk[2][A_X]=HYPER[k].t_inverse_Fi[2][A_X];	Fk[2][A_Y]=HYPER[k].t_inverse_Fi[2][A_Y];	Fk[2][A_Z]=HYPER[k].t_inverse_Fi[2][A_Z];

			Dgkk[A_X]=Jk*(Fk[0][0]*HYPER1[k*h_num+k].n0ij[0]+Fk[0][1]*HYPER1[k*h_num+k].n0ij[1]+Fk[0][2]*HYPER1[k*h_num+k].n0ij[2]);
			Dgkk[A_Y]=Jk*(Fk[1][0]*HYPER1[k*h_num+k].n0ij[0]+Fk[1][1]*HYPER1[k*h_num+k].n0ij[1]+Fk[1][2]*HYPER1[k*h_num+k].n0ij[2]);
			Dgkk[A_Z]=Jk*(Fk[2][0]*HYPER1[k*h_num+k].n0ij[0]+Fk[2][1]*HYPER1[k*h_num+k].n0ij[1]+Fk[2][2]*HYPER1[k*h_num+k].n0ij[2]);


			for(int i=0;i<Nk;i++)
			{
				int in=HYPER[k].NEI[i];
				//if(in>=Nw)
				{
					Jin=HYPER[in].J;

					DgDq_ik[A_X]=HYPER1[in*h_num+k].DgDq_n[A_X];
					DgDq_ik[A_Y]=HYPER1[in*h_num+k].DgDq_n[A_Y];
					DgDq_ik[A_Z]=HYPER1[in*h_num+k].DgDq_n[A_Z];

					Fin[0][A_X]=HYPER[in].t_inverse_Fi[0][A_X];	Fin[0][A_Y]=HYPER[in].t_inverse_Fi[0][A_Y];	Fin[0][A_Z]=HYPER[in].t_inverse_Fi[0][A_Z];
					Fin[1][A_X]=HYPER[in].t_inverse_Fi[1][A_X];	Fin[1][A_Y]=HYPER[in].t_inverse_Fi[1][A_Y];	Fin[1][A_Z]=HYPER[in].t_inverse_Fi[1][A_Z];
					Fin[2][A_X]=HYPER[in].t_inverse_Fi[2][A_X];	Fin[2][A_Y]=HYPER[in].t_inverse_Fi[2][A_Y];	Fin[2][A_Z]=HYPER[in].t_inverse_Fi[2][A_Z];

					Dgik[A_X]=Jin*(Fin[0][0]*HYPER1[k*h_num+in].n0ij[0]+Fin[0][1]*HYPER1[k*h_num+in].n0ij[1]+Fin[0][2]*HYPER1[k*h_num+in].n0ij[2]);
					Dgik[A_Y]=Jin*(Fin[1][0]*HYPER1[k*h_num+in].n0ij[0]+Fin[1][1]*HYPER1[k*h_num+in].n0ij[1]+Fin[1][2]*HYPER1[k*h_num+in].n0ij[2]);
					Dgik[A_Z]=Jin*(Fin[2][0]*HYPER1[k*h_num+in].n0ij[0]+Fin[2][1]*HYPER1[k*h_num+in].n0ij[1]+Fin[2][2]*HYPER1[k*h_num+in].n0ij[2]);

					for(int j=0;j<Nk;j++)
					{
						int jn=HYPER[k].NEI[j];

						DgDq_jk[A_X]=HYPER1[jn*h_num+k].DgDq_n[A_X];
						DgDq_jk[A_Y]=HYPER1[jn*h_num+k].DgDq_n[A_Y];
						DgDq_jk[A_Z]=HYPER1[jn*h_num+k].DgDq_n[A_Z];

						DfDx[in*h_num+jn]-=Dt*Dt*0.5/mi*(Dgik[A_X]*DgDq_jk[A_X]+Dgik[A_Y]*DgDq_jk[A_Y]+Dgik[A_Z]*DgDq_jk[A_Z]);
					}
					DfDx[in*h_num+k]-=Dt*Dt*0.5/mi*(Dgik[A_X]*DgDq_kk[A_X]+Dgik[A_Y]*DgDq_kk[A_Y]+Dgik[A_Z]*DgDq_kk[A_Z]);
				}
				DfDx[k*h_num+in]-=Dt*Dt*0.5/mi*(Dgkk[A_X]*DgDq_ik[A_X]+Dgkk[A_Y]*DgDq_ik[A_Y]+Dgkk[A_Z]*DgDq_ik[A_Z]);
			}//*/
			DfDx[k*h_num+k]-=Dt*Dt*0.5/mi*(Dgkk[A_X]*DgDq_kk[A_X]+Dgkk[A_Y]*Dgkk[A_Y]+Dgkk[A_Z]*Dgkk[A_Z]);		
		}
		//#ifdef _OPENMP
		//#pragma omp for
		//#endif
		//		for(int k=Nw;k<h_num;k++)
		//		{			
		//			int Nk=HYPER[k].N;
		//			DgDq_kk[A_X]=HYPER1[k*h_num+k].DgDq_n[A_X];
		//			DgDq_kk[A_Y]=HYPER1[k*h_num+k].DgDq_n[A_Y];
		//			DgDq_kk[A_Z]=HYPER1[k*h_num+k].DgDq_n[A_Z];
		//
		//			Jk=HYPER[k].J;
		//			Fk[0][A_X]=HYPER[k].t_inverse_Fi[0][A_X];	Fk[0][A_Y]=HYPER[k].t_inverse_Fi[0][A_Y];	Fk[0][A_Z]=HYPER[k].t_inverse_Fi[0][A_Z];
		//			Fk[1][A_X]=HYPER[k].t_inverse_Fi[1][A_X];	Fk[1][A_Y]=HYPER[k].t_inverse_Fi[1][A_Y];	Fk[1][A_Z]=HYPER[k].t_inverse_Fi[1][A_Z];
		//			Fk[2][A_X]=HYPER[k].t_inverse_Fi[2][A_X];	Fk[2][A_Y]=HYPER[k].t_inverse_Fi[2][A_Y];	Fk[2][A_Z]=HYPER[k].t_inverse_Fi[2][A_Z];
		//
		//			Dgkk[A_X]=Jk*(Fk[0][0]*HYPER1[k*h_num+k].n0ij[0]+Fk[0][1]*HYPER1[k*h_num+k].n0ij[1]+Fk[0][2]*HYPER1[k*h_num+k].n0ij[2]);
		//			Dgkk[A_Y]=Jk*(Fk[1][0]*HYPER1[k*h_num+k].n0ij[0]+Fk[1][1]*HYPER1[k*h_num+k].n0ij[1]+Fk[1][2]*HYPER1[k*h_num+k].n0ij[2]);
		//			Dgkk[A_Z]=Jk*(Fk[2][0]*HYPER1[k*h_num+k].n0ij[0]+Fk[2][1]*HYPER1[k*h_num+k].n0ij[1]+Fk[2][2]*HYPER1[k*h_num+k].n0ij[2]);
		//
		//
		//			for(int i=0;i<Nk;i++)
		//			{
		//				int in=HYPER[k].NEI[i];
		//				if(in>=Nw)
		//				{
		//					Jin=HYPER[in].J;
		//
		//					DgDq_ik[A_X]=HYPER1[in*h_num+k].DgDq_n[A_X];
		//					DgDq_ik[A_Y]=HYPER1[in*h_num+k].DgDq_n[A_Y];
		//					DgDq_ik[A_Z]=HYPER1[in*h_num+k].DgDq_n[A_Z];
		//
		//					Fin[0][A_X]=HYPER[in].t_inverse_Fi[0][A_X];	Fin[0][A_Y]=HYPER[in].t_inverse_Fi[0][A_Y];	Fin[0][A_Z]=HYPER[in].t_inverse_Fi[0][A_Z];
		//					Fin[1][A_X]=HYPER[in].t_inverse_Fi[1][A_X];	Fin[1][A_Y]=HYPER[in].t_inverse_Fi[1][A_Y];	Fin[1][A_Z]=HYPER[in].t_inverse_Fi[1][A_Z];
		//					Fin[2][A_X]=HYPER[in].t_inverse_Fi[2][A_X];	Fin[2][A_Y]=HYPER[in].t_inverse_Fi[2][A_Y];	Fin[2][A_Z]=HYPER[in].t_inverse_Fi[2][A_Z];
		//
		//					Dgik[A_X]=Jin*(Fin[0][0]*HYPER1[k*h_num+in].n0ij[0]+Fin[0][1]*HYPER1[k*h_num+in].n0ij[1]+Fin[0][2]*HYPER1[k*h_num+in].n0ij[2]);
		//					Dgik[A_Y]=Jin*(Fin[1][0]*HYPER1[k*h_num+in].n0ij[0]+Fin[1][1]*HYPER1[k*h_num+in].n0ij[1]+Fin[1][2]*HYPER1[k*h_num+in].n0ij[2]);
		//					Dgik[A_Z]=Jin*(Fin[2][0]*HYPER1[k*h_num+in].n0ij[0]+Fin[2][1]*HYPER1[k*h_num+in].n0ij[1]+Fin[2][2]*HYPER1[k*h_num+in].n0ij[2]);
		//
		//					for(int j=0;j<Nk;j++)
		//					{
		//						int jn=HYPER[k].NEI[j];
		//
		//						DgDq_jk[A_X]=HYPER1[jn*h_num+k].DgDq_n[A_X];
		//						DgDq_jk[A_Y]=HYPER1[jn*h_num+k].DgDq_n[A_Y];
		//						DgDq_jk[A_Z]=HYPER1[jn*h_num+k].DgDq_n[A_Z];
		//
		//						DfDx[in*h_num+jn]-=Dt*Dt*0.5/mi*(Dgik[A_X]*DgDq_jk[A_X]+Dgik[A_Y]*DgDq_jk[A_Y]+Dgik[A_Z]*DgDq_jk[A_Z]);
		//					}
		//					DfDx[in*h_num+k]-=Dt*Dt*0.5/mi*(Dgik[A_X]*DgDq_kk[A_X]+Dgik[A_Y]*DgDq_kk[A_Y]+Dgik[A_Z]*DgDq_kk[A_Z]);
		//				}
		//
		//				DfDx[k*h_num+in]-=Dt*Dt*0.5/mi*(Dgkk[A_X]*DgDq_ik[A_X]+Dgkk[A_Y]*DgDq_ik[A_Y]+Dgkk[A_Z]*DgDq_ik[A_Z]);
		//			}
		//			DfDx[k*h_num+k]-=Dt*Dt*0.5/mi*(Dgkk[A_X]*DgDq_kk[A_X]+Dgkk[A_Y]*Dgkk[A_Y]+Dgkk[A_Z]*Dgkk[A_Z]);		
		//		}//*/
		//
		double DgDq_ji[DIMENSION]={0,0,0};
		double DgDq_ii[DIMENSION]={0,0,0};

#ifdef _OPENMP
#pragma omp for
#endif
		for(int i=0;i<Nw;i++)
		{
			int Ni=HYPER[i].N;

			DgDq_ii[A_X]=HYPER1[i*h_num+i].DgDq_n[A_X];
			DgDq_ii[A_Y]=HYPER1[i*h_num+i].DgDq_n[A_Y];
			DgDq_ii[A_Z]=HYPER1[i*h_num+i].DgDq_n[A_Z];

			for(int j=0;j<Ni;j++)
			{
				int jn=HYPER[i].NEI[j];
				DgDq_ji[A_X]=HYPER1[jn*h_num+i].DgDq_n[A_X];
				DgDq_ji[A_Y]=HYPER1[jn*h_num+i].DgDq_n[A_Y];
				DgDq_ji[A_Z]=HYPER1[jn*h_num+i].DgDq_n[A_Z];

				//DfDx[i*h_num+jn]+=-Dt*Dt/mi*(DgDq_ji[A_X]*(PART[i].r[A_X]-PART[i].q0[A_X])+DgDq_ji[A_Y]*(PART[i].r[A_Y]-PART[i].q0[A_Y])+DgDq_ji[A_Z]*(PART[i].r[A_Z]-PART[i].q0[A_Z]));
				DfDx[i*h_num+jn]+=-0.5*Dt*Dt/mi*(DgDq_ji[A_X]+DgDq_ji[A_Y]+DgDq_ji[A_Z]);
				//cout<<"DfDx"<<jn<<"_"<<i<<"="<<DfDx[jn*h_num+i]<<endl;
			}//*/
			DfDx[i*h_num+i]+=-0.5*Dt*Dt/mi*(DgDq_ii[A_X]+DgDq_ii[A_Y]+DgDq_ii[A_Z]);
			//DfDx[i*h_num+i]+=-Dt*Dt/mi*(DgDq_ii[A_X]*(PART[i].r[A_X]-PART[i].q0[A_X])+DgDq_ii[A_Y]*(PART[i].r[A_Y]-PART[i].q0[A_Y])+DgDq_ii[A_Z]*(PART[i].r[A_Z]-PART[i].q0[A_Z]));
			//cout<<"DfDx"<<i<<"_"<<i<<"="<<DfDx[i*h_num+i]<<endl;
		}
		//double *part_fx=new double [h_num];
		//double *w_DfDx=new double [h_num*h_num];
		//double *w_fx=new double [h_num];
		//int *Nw=new int [h_num];

		//for(int i=0;i<h_num;i++)
		//{
		//for(int j=0;j<h_num;j++)	w_DfDx[i*h_num+j]=0;
		//w_fx[i]=0;
		//Nw[i]=0;
		//}

		//int number=0;
		//for(int i=0;i<h_num;i++)
		//{
		//part_fx[i]=0;
		//if(n_rz[i]<0)
		//{
		//int Ni=HYPER[i].N;
		//DfDx[i*h_num+i]=1;
		//for(int j=0;j<Ni;j++)
		//{
		//Nw[number]=i;
		//int jn=HYPER[i].NEI[j];
		//w_DfDx[i*h_num+jn]=Dt/2*HYPER1[jn*h_num+i].DgDq[A_Z];
		//part_fx[i]+=Dt/2*(HYPER[jn].stress[A_Z][A_X]*HYPER1[jn*h_num+i].DgDq[A_X]+HYPER[jn].stress[A_Z][A_Y]*HYPER1[jn*h_num+i].DgDq[A_Y]+HYPER[jn].stress[A_Z][A_Z]*HYPER1[jn*h_num+i].DgDq[A_Z]);
		//DfDx[i*h_num+jn]=0;
		//}
		//w_fx[i]=HYPER[i].p[A_Z]+part_fx[i];
		//number++;
		//}
		////		w_DfDx[i*h_num+i]=1;
		//}
		//if(number>0)
		//{

		//output_newton_data3(w_fx,w_DfDx,n_rx,n_ry,n_rz,t_fx,h_num,count,t);		
		//gauss(w_DfDx,w_fx,h_num);

		//for(int i=0;i<number;i++)
		//{
		//int in=Nw[i];
		//fx[in]=w_fx[in];
		//}
		//}
		//delete[]	Nw;
		//delete[]	w_DfDx;
		//delete[]	w_fx;
		//delete[]	part_fx;//*/

		//if(count==1 || count%10000==0)	output_newton_data(PART,HYPER,HYPER1,fx,DfDx,h_num,count,t);




		////出力
		//	if(count%200==0 && count>CON.get_nr()/2)
		//	if(t==1||t%CON.get_interval()==0)	if(count%200==0||count==1)	
		//	if(t==1||t%(10*CON.get_interval())==0)	if(count%200==0||count==1)	output_newton_data1(fx,DfDx,n_rx,n_ry,n_rz,h_num,count,t);

		//ofstream t_loge("time_log_newton.dat", ios::app);
		clock_t t4=clock();
		//t_loge<<"step="<<t<<", count="<<count<<", time="<<1000*(t4-t3)/CLOCKS_PER_SEC<<"[e-3sec]"<<endl;
		//t_loge.close();


		for(int D=0;D<DIMENSION;D++)	delete[]	p_Fi[D];
		delete[]	p_Fi;
	}
}

void calc_half_p(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,bool repetation,double **F,int Nw)
{
	//	if(repetation==0)	cout<<"仮の運動量＆位置座標計算";
	//	else	cout<<"運動量計算";

	int h_num=HYPER.size();
	double Dt=CON.get_dt();
	double le=CON.get_distancebp();
	double V=get_volume(&CON);
	double mh=V*CON.get_hyper_density();
	double ms=V*CON.get_silicone_density();
	double G=9.8;
	int flag_vis=CON.get_flag_vis();
	bool flag_FEM=CON.get_FEM_flag();
	int flag_G=CON.get_flag_G();
	double density=CON.get_hyper_density();
	int model_num=CON.get_model_number();
	double mi=V*CON.get_hyper_density();

#ifdef _OPENMP
	omp_set_num_threads(8);
#pragma omp parallel for
#endif
	{
		if(repetation==0)
		{

#ifdef _OPENMP
#pragma omp for
#endif
			for(int i=0;i<h_num;i++)
			{
				//if(PART[i].toFEM==ON)	mi=mh;
				//else
				//{
				//	mi=ms;

				//}
				//if(HYPER[i].fw==1)
				//{
				//	HYPER[i].half_p[A_X]=0;
				//	HYPER[i].half_p[A_Y]=0;
				//	HYPER[i].half_p[A_Z]=0;//
				//	//位置座標の更新
				//	PART[i].r[A_X]=PART[i].q0[A_X];
				//	PART[i].r[A_Y]=PART[i].q0[A_Y];
				//	PART[i].r[A_Z]=PART[i].q0[A_Z];

				//}
				//else
				//{
				double p_half_p[DIMENSION]={0,0,0};
				int Ni=HYPER[i].N;
				for(int j=0;j<Ni;j++)
				{		
					int k=HYPER[i].NEI[j];
					p_half_p[A_X]+=(HYPER[k].stress_n[0][0]-HYPER[k].lambda)*HYPER1[k*h_num+i].DgDq_n[0]+HYPER[k].stress_n[0][1]*HYPER1[k*h_num+i].DgDq_n[1]+HYPER[k].stress_n[0][2]*HYPER1[k*h_num+i].DgDq_n[2];
					p_half_p[A_Y]+=HYPER[k].stress_n[1][0]*HYPER1[k*h_num+i].DgDq_n[0]+(HYPER[k].stress_n[1][1]-HYPER[k].lambda)*HYPER1[k*h_num+i].DgDq_n[1]+HYPER[k].stress_n[1][2]*HYPER1[k*h_num+i].DgDq_n[2];
					p_half_p[A_Z]+=HYPER[k].stress_n[2][0]*HYPER1[k*h_num+i].DgDq_n[0]+HYPER[k].stress_n[2][1]*HYPER1[k*h_num+i].DgDq_n[1]+(HYPER[k].stress_n[2][2]-HYPER[k].lambda)*HYPER1[k*h_num+i].DgDq_n[2];
				}//jに関するfor文の終わり	
				p_half_p[A_X]+=(HYPER[i].stress_n[0][0]-HYPER[i].lambda)*HYPER1[i*h_num+i].DgDq_n[0]+HYPER[i].stress_n[0][1]*HYPER1[i*h_num+i].DgDq_n[1]+HYPER[i].stress_n[0][2]*HYPER1[i*h_num+i].DgDq_n[2];
				p_half_p[A_Y]+=HYPER[i].stress_n[1][0]*HYPER1[i*h_num+i].DgDq_n[0]+(HYPER[i].stress_n[1][1]-HYPER[i].lambda)*HYPER1[i*h_num+i].DgDq_n[1]+HYPER[i].stress_n[1][2]*HYPER1[i*h_num+i].DgDq_n[2];
				p_half_p[A_Z]+=HYPER[i].stress_n[2][0]*HYPER1[i*h_num+i].DgDq_n[0]+HYPER[i].stress_n[2][1]*HYPER1[i*h_num+i].DgDq_n[1]+(HYPER[i].stress_n[2][2]-HYPER[i].lambda)*HYPER1[i*h_num+i].DgDq_n[2];

				//重力項

				/*if(flag_G==ON)	*/p_half_p[A_Z]-=G*mi;
				//粘性項
				//if(flag_vis==ON)
				{
					p_half_p[A_X]+=HYPER[i].vis_force[A_X];
					p_half_p[A_Y]+=HYPER[i].vis_force[A_Y];
					p_half_p[A_Z]+=HYPER[i].vis_force[A_Z];
				}
				//磁力項
				if(flag_FEM==ON&&PART[i].toFEM==ON)
				{
					p_half_p[A_X]+=V*F[A_X][i];//density;
					p_half_p[A_Y]+=V*F[A_Y][i];//density;
					p_half_p[A_Z]+=V*F[A_Z][i];//density;
				}


				//half_pの更新
				HYPER[i].half_p[A_X]=HYPER[i].p_n[A_X]+Dt*0.5*p_half_p[A_X];
				HYPER[i].half_p[A_Y]=HYPER[i].p_n[A_Y]+Dt*0.5*p_half_p[A_Y];
				HYPER[i].half_p[A_Z]=HYPER[i].p_n[A_Z]+Dt*0.5*p_half_p[A_Z];//
				PART[i].r[A_X]=HYPER[i].q_n[A_X]+Dt*HYPER[i].half_p[A_X]/mi;
				PART[i].r[A_Y]=HYPER[i].q_n[A_Y]+Dt*HYPER[i].half_p[A_Y]/mi;
				PART[i].r[A_Z]=HYPER[i].q_n[A_Z]+Dt*HYPER[i].half_p[A_Z]/mi;

				//位置座標の更新
				//if(i<Nw)
				//{
				//	HYPER[i].half_p[A_Z]=0.;
				//	PART[i].r[A_Z]=PART[i].q0[A_Z];
				//}

			}
			//}

		}
		else
		{

#ifdef _OPENMP
#pragma omp for
#endif
			for(int i=0;i<h_num;i++)
			{
				//if(PART[i].toFEM==ON)	mi=mh;
				//else
				//{
				//	mi=ms;

				//}
				//if(HYPER[i].fw==1)
				//{
				//	HYPER[i].p[A_X]=0;
				//	HYPER[i].p[A_Y]=0;
				//	HYPER[i].p[A_Z]=0;//
				//}
				//else
				//{
				double p_half_p[DIMENSION]={0,0,0};
				int Ni=HYPER[i].N;
				for(int j=0;j<Ni;j++)
				{		
					int k=HYPER[i].NEI[j];
					p_half_p[A_X]+=(HYPER[k].stress[0][0]-HYPER[k].lambda)*HYPER1[k*h_num+i].DgDq[0]+HYPER[k].stress[0][1]*HYPER1[k*h_num+i].DgDq[1]+HYPER[k].stress[0][2]*HYPER1[k*h_num+i].DgDq[2];
					p_half_p[A_Y]+=HYPER[k].stress[1][0]*HYPER1[k*h_num+i].DgDq[0]+(HYPER[k].stress[1][1]-HYPER[k].lambda)*HYPER1[k*h_num+i].DgDq[1]+HYPER[k].stress[1][2]*HYPER1[k*h_num+i].DgDq[2];
					p_half_p[A_Z]+=HYPER[k].stress[2][0]*HYPER1[k*h_num+i].DgDq[0]+HYPER[k].stress[2][1]*HYPER1[k*h_num+i].DgDq[1]+(HYPER[k].stress[2][2]-HYPER[k].lambda)*HYPER1[k*h_num+i].DgDq[2];
				}//jに関するfor文の終わり	
				p_half_p[A_X]+=(HYPER[i].stress[0][0]-HYPER[i].lambda)*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[0][1]*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[0][2]*HYPER1[i*h_num+i].DgDq[2];
				p_half_p[A_Y]+=HYPER[i].stress[1][0]*HYPER1[i*h_num+i].DgDq[0]+(HYPER[i].stress[1][1]-HYPER[i].lambda)*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[1][2]*HYPER1[i*h_num+i].DgDq[2];
				p_half_p[A_Z]+=HYPER[i].stress[2][0]*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[2][1]*HYPER1[i*h_num+i].DgDq[1]+(HYPER[i].stress[2][2]-HYPER[i].lambda)*HYPER1[i*h_num+i].DgDq[2];

				//重力項
				/*if(flag_G==ON)	*/p_half_p[A_Z]-=G*mi;
				//粘性項
				//if(flag_vis==ON)
				{
					p_half_p[A_X]+=HYPER[i].vis_force[A_X];
					p_half_p[A_Y]+=HYPER[i].vis_force[A_Y];
					p_half_p[A_Z]+=HYPER[i].vis_force[A_Z];
				}
				//磁力項
				if(flag_FEM==ON&&PART[i].toFEM==ON)
				{
					p_half_p[A_X]+=V*F[A_X][i];//density;
					p_half_p[A_Y]+=V*F[A_Y][i];//density;
					p_half_p[A_Z]+=V*F[A_Z][i];//density;
				}

				//運動量の更新
				HYPER[i].p[A_X]=HYPER[i].half_p[A_X]+Dt*0.5*p_half_p[A_X];
				HYPER[i].p[A_Y]=HYPER[i].half_p[A_Y]+Dt*0.5*p_half_p[A_Y];
				HYPER[i].p[A_Z]=HYPER[i].half_p[A_Z]+Dt*0.5*p_half_p[A_Z];////
				//if(i<Nw)	HYPER[i].p[A_Z]=0.;

				//速度の更新
				PART[i].u[A_X]=HYPER[i].half_p[A_X]/mi;
				PART[i].u[A_Y]=HYPER[i].half_p[A_Y]/mi;
				PART[i].u[A_Z]=HYPER[i].half_p[A_Z]/mi;
				//角運動量の更新
				HYPER[i].ang_p[A_X]=PART[i].r[A_Y]*HYPER[i].p[A_Z]-PART[i].r[A_Z]*HYPER[i].p[A_Y];
				HYPER[i].ang_p[A_Y]=PART[i].r[A_Z]*HYPER[i].p[A_X]-PART[i].r[A_X]*HYPER[i].p[A_Z];
				HYPER[i].ang_p[A_Z]=PART[i].r[A_X]*HYPER[i].p[A_Y]-PART[i].r[A_Y]*HYPER[i].p[A_X];
			}
		}//iに関するfor文の終わり
		//}

		//	cout<<"----------OK"<<endl;
		//}
	}
}

void calc_F(mpsconfig &CON, vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1)
{
	//	cout<<"Fi計算";
	////Fiの更新
#ifdef _OPENMP
	omp_set_num_threads(8);
#pragma omp parallel for
#endif
	{
		int h_num=HYPER.size();
		double V=get_volume(&CON);

		double **p_Fi=new double *[DIMENSION];
		for(int D=0;D<DIMENSION;D++)	p_Fi[D]=new double[DIMENSION];

#ifdef _OPENMP
#pragma omp for
#endif
		for(int i=0;i<h_num;i++)
		{
			double fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
			//Fiの計算

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
			//Jの計算
			double J=calc_det3(p_Fi);
			//	for(int i=0;i<h_num;i++)	cout<<"J["<<i<<"]="<<J<<endl;
			HYPER[i].J=J;

			//t_inverse_Fiの計算
			inverse(p_Fi,DIMENSION);
			HYPER[i].t_inverse_Fi[0][0]=p_Fi[0][0];	HYPER[i].t_inverse_Fi[0][1]=p_Fi[1][0];	HYPER[i].t_inverse_Fi[0][2]=p_Fi[2][0];
			HYPER[i].t_inverse_Fi[1][0]=p_Fi[0][1];	HYPER[i].t_inverse_Fi[1][1]=p_Fi[1][1];	HYPER[i].t_inverse_Fi[1][2]=p_Fi[2][1];
			HYPER[i].t_inverse_Fi[2][0]=p_Fi[0][2];	HYPER[i].t_inverse_Fi[2][1]=p_Fi[1][2];	HYPER[i].t_inverse_Fi[2][2]=p_Fi[2][2];


			for(int j=0;j<Ni;j++)
			{			
				int k=HYPER[i].NEI[j];
				HYPER1[i*h_num+k].DgDq[A_X]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_X][0]*HYPER1[k*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_X][1]*HYPER1[k*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_X][2]*HYPER1[k*h_num+i].n0ij[2]);
				HYPER1[i*h_num+k].DgDq[A_Y]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_Y][0]*HYPER1[k*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_Y][1]*HYPER1[k*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_Y][2]*HYPER1[k*h_num+i].n0ij[2]);
				HYPER1[i*h_num+k].DgDq[A_Z]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_Z][0]*HYPER1[k*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_Z][1]*HYPER1[k*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_Z][2]*HYPER1[k*h_num+i].n0ij[2]);
				//cout<<"i"<<i<<"j"<<k<<"	"<<HYPER1[k*h_num+i].DgDq[A_X]<<","<<HYPER1[k*h_num+i].DgDq[A_Y]<<","<<HYPER1[k*h_num+i].DgDq[A_Z]<<endl;

			}
			HYPER1[i*h_num+i].DgDq[A_X]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_X][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_X][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_X][2]*HYPER1[i*h_num+i].n0ij[2]);
			HYPER1[i*h_num+i].DgDq[A_Y]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_Y][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_Y][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_Y][2]*HYPER1[i*h_num+i].n0ij[2]);
			HYPER1[i*h_num+i].DgDq[A_Z]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_Z][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_Z][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_Z][2]*HYPER1[i*h_num+i].n0ij[2]);
		}
		//	cout<<"----------OK"<<endl;

		for(int D=0;D<DIMENSION;D++)	delete[]	p_Fi[D];
		delete[]	p_Fi;
	}

	//void calc_pi(mpsconfig &CON,vector<hyperelastic> &HYPER)
	//{
	//	int h_num=HYPER.size();
	//	double Dt=CON.get_dt();
	//	double V=get_volume(&CON);
	//	double mi=V*CON.get_hyper_density();
	//	double c10=CON.get_c10();
	//	double c01=CON.get_c01();
	//	double nG[DIMENSION]={0,0,1};
	//
	//	/////////////DgDq, Stress, W, S, dSdc計算
	//	double d_Fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	//	
	//	double dC[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	//	double dC2[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	//
	//	double trace_dC=0;
	//	double trace_dC2=0;
	//	double Ic=0;
	//	double IIc=0;
	//
	//	double **in_Ci=new double *[DIMENSION];
	//	for(int D=0;D<DIMENSION;D++)	in_Ci[D]=new double[DIMENSION];
	//
	//	double first_term[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	//
	//	double S[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	//
	//
	//	for(int j=0;j<h_num;j++)
	//	{	
	//
	//		double J=HYPER[j].J;	
	//		if(J<0){
	//			d_Fi[0][0]=-1/pow(-J,1/3.)*HYPER[j].Fi[0][0];	d_Fi[0][1]=-1/pow(-J,1/3.)*HYPER[j].Fi[0][1];	d_Fi[0][2]=-1/pow(-J,1/3.)*HYPER[j].Fi[0][2];
	//			d_Fi[1][0]=-1/pow(-J,1/3.)*HYPER[j].Fi[1][0];	d_Fi[1][1]=-1/pow(-J,1/3.)*HYPER[j].Fi[1][1];	d_Fi[1][2]=-1/pow(-J,1/3.)*HYPER[j].Fi[1][2];
	//			d_Fi[2][0]=-1/pow(-J,1/3.)*HYPER[j].Fi[2][0];	d_Fi[2][1]=-1/pow(-J,1/3.)*HYPER[j].Fi[2][1];	d_Fi[2][2]=-1/pow(-J,1/3.)*HYPER[j].Fi[2][2];
	//		}
	//		else
	//		{
	//			d_Fi[0][0]=1/pow(J,1/3.)*HYPER[j].Fi[0][0];	d_Fi[0][1]=1/pow(J,1/3.)*HYPER[j].Fi[0][1];	d_Fi[0][2]=1/pow(J,1/3.)*HYPER[j].Fi[0][2];
	//			d_Fi[1][0]=1/pow(J,1/3.)*HYPER[j].Fi[1][0];	d_Fi[1][1]=1/pow(J,1/3.)*HYPER[j].Fi[1][1];	d_Fi[1][2]=1/pow(J,1/3.)*HYPER[j].Fi[1][2];
	//			d_Fi[2][0]=1/pow(J,1/3.)*HYPER[j].Fi[2][0];	d_Fi[2][1]=1/pow(J,1/3.)*HYPER[j].Fi[2][1];	d_Fi[2][2]=1/pow(J,1/3.)*HYPER[j].Fi[2][2];
	//		}
	//		
	//		//cout<<"dF"<<j<<"="<<d_Fi[0][0]<<","<<d_Fi[0][1]<<","<<d_Fi[0][2]<<endl;
	//		//cout<<d_Fi[1][0]<<","<<d_Fi[1][1]<<","<<d_Fi[1][2]<<endl;
	//		//cout<<d_Fi[2][0]<<","<<d_Fi[2][1]<<","<<d_Fi[2][2]<<endl;
	//
	//		dC[0][0]=d_Fi[0][0]*d_Fi[0][0]+d_Fi[1][0]*d_Fi[1][0]+d_Fi[2][0]*d_Fi[2][0];
	//		dC[0][1]=d_Fi[0][0]*d_Fi[0][1]+d_Fi[1][0]*d_Fi[1][1]+d_Fi[2][0]*d_Fi[2][1];
	//		dC[0][2]=d_Fi[0][0]*d_Fi[0][2]+d_Fi[1][0]*d_Fi[1][2]+d_Fi[2][0]*d_Fi[2][2];
	//		dC[1][0]=d_Fi[0][1]*d_Fi[0][0]+d_Fi[1][1]*d_Fi[1][0]+d_Fi[2][1]*d_Fi[2][0];
	//		dC[1][1]=d_Fi[0][1]*d_Fi[0][1]+d_Fi[1][1]*d_Fi[1][1]+d_Fi[2][1]*d_Fi[2][1];
	//		dC[1][2]=d_Fi[0][1]*d_Fi[0][2]+d_Fi[1][1]*d_Fi[1][2]+d_Fi[2][1]*d_Fi[2][2];
	//		dC[2][0]=d_Fi[0][2]*d_Fi[0][0]+d_Fi[1][2]*d_Fi[1][0]+d_Fi[2][2]*d_Fi[2][0];
	//		dC[2][1]=d_Fi[0][2]*d_Fi[0][1]+d_Fi[1][2]*d_Fi[1][1]+d_Fi[2][2]*d_Fi[2][1];
	//		dC[2][2]=d_Fi[0][2]*d_Fi[0][2]+d_Fi[1][2]*d_Fi[1][2]+d_Fi[2][2]*d_Fi[2][2];
	//
	//		dC2[0][0]=dC[0][0]*dC[0][0]+dC[0][1]*dC[1][0]+dC[0][2]*dC[2][0];
	//		dC2[0][1]=dC[0][0]*dC[0][1]+dC[0][1]*dC[1][1]+dC[0][2]*dC[2][1];
	//		dC2[0][2]=dC[0][0]*dC[0][2]+dC[0][1]*dC[1][2]+dC[0][2]*dC[2][2];
	//		dC2[1][0]=dC[1][0]*dC[0][0]+dC[1][1]*dC[1][0]+dC[1][2]*dC[2][0];
	//		dC2[1][1]=dC[1][0]*dC[0][1]+dC[1][1]*dC[1][1]+dC[1][2]*dC[2][1];
	//		dC2[1][2]=dC[1][0]*dC[0][2]+dC[1][1]*dC[1][2]+dC[1][2]*dC[2][2];
	//		dC2[2][0]=dC[2][0]*dC[0][0]+dC[2][1]*dC[1][0]+dC[2][2]*dC[2][0];
	//		dC2[2][1]=dC[2][0]*dC[0][1]+dC[2][1]*dC[1][1]+dC[2][2]*dC[2][1];
	//		dC2[2][2]=dC[2][0]*dC[0][2]+dC[2][1]*dC[1][2]+dC[2][2]*dC[2][2];
	//
	//		//cout<<"dC"<<j<<"="<<dC[0][0]<<","<<dC[0][1]<<","<<dC[0][2]<<endl;
	//		//cout<<dC[1][0]<<","<<dC[1][1]<<","<<dC[1][2]<<endl;
	//		//cout<<dC[2][0]<<","<<dC[2][1]<<","<<dC[2][2]<<endl;
	//		//
	//		//cout<<"dC2"<<j<<"="<<dC2[0][0]<<","<<dC2[0][1]<<","<<dC2[0][2]<<endl;
	//		//cout<<dC2[1][0]<<","<<dC2[1][1]<<","<<dC2[1][2]<<endl;
	//		//cout<<dC2[2][0]<<","<<dC2[2][1]<<","<<dC2[2][2]<<endl;
	//
	//
	//		trace_dC=dC[0][0]+dC[1][1]+dC[2][2];
	//		trace_dC2=dC2[0][0]+dC2[1][1]+dC2[2][2];
	//
	//		//cout<<"trace_dC="<<trace_dC<<", trace_dC2="<<trace_dC2;
	//
	//		Ic=trace_dC;
	//		IIc=0.50*(trace_dC*trace_dC-trace_dC2);
	//
	//		//cout<<" ,Ic="<<Ic<<", IIc="<<IIc<<endl;
	//	
	//		
	//		////////S計算
	//		in_Ci[0][0]=dC[0][0];	in_Ci[0][1]=dC[0][1];	in_Ci[0][2]=dC[0][2];
	//		in_Ci[1][0]=dC[1][0];	in_Ci[1][1]=dC[1][1];	in_Ci[1][2]=dC[1][2];
	//		in_Ci[2][0]=dC[2][0];	in_Ci[2][1]=dC[2][1];	in_Ci[2][2]=dC[2][2];
	//
	//		inverse(in_Ci,DIMENSION);
	//
	//
	//		//cout<<"in_Ci"<<j<<"="<<in_Ci[0][0]<<","<<in_Ci[0][1]<<","<<in_Ci[0][2]<<endl;
	//		//cout<<in_Ci[1][0]<<","<<in_Ci[1][1]<<","<<in_Ci[1][2]<<endl;
	//		//cout<<in_Ci[2][0]<<","<<in_Ci[2][1]<<","<<in_Ci[2][2]<<endl;
	//
	//		first_term[0][0]=1-1/3.*Ic*in_Ci[0][0];	first_term[0][1]=-1/3.*Ic*in_Ci[0][1];	first_term[0][2]=-1/3.*Ic*in_Ci[0][2];
	//		first_term[1][0]=-1/3.*Ic*in_Ci[1][0];	first_term[1][1]=1-1/3.*Ic*in_Ci[1][1];	first_term[1][2]=-1/3.*Ic*in_Ci[1][0];
	//		first_term[2][0]=-1/3.*Ic*in_Ci[2][0];	first_term[2][1]=-1/3.*Ic*in_Ci[2][1];	first_term[2][2]=1-1/3.*Ic*in_Ci[2][2];
	//		
	//		//cout<<"1/3.*Ic*in_Ci[0][0]="<<1/3.*Ic*in_Ci[0][0]<<endl;
	//		//cout<<"in_Ci[0][0]="<<in_Ci[0][0]<<endl;
	//		//cout<<"Ic="<<Ic<<endl;
	//		//cout<<"Ic*in_Ci[0][0]="<<Ic*in_Ci[0][0]<<endl;
	//		//cout<<"first_term"<<j<<"="<<first_term[0][0]<<","<<first_term[0][1]<<","<<first_term[0][2]<<endl;
	//		//cout<<first_term[1][0]<<","<<first_term[1][1]<<","<<first_term[1][2]<<endl;
	//		//cout<<first_term[2][0]<<","<<first_term[2][1]<<","<<first_term[2][2]<<endl;
	//
	//		if(J<0)
	//		{
	//			S[0][0]=-2*1/pow(-J,2/3.)*( (c10+c01*Ic)*first_term[0][0]-c01*(dC[0][0]-1/3.*trace_dC2*in_Ci[0][0]) );
	//			S[0][1]=-2*1/pow(-J,2/3.)*( (c10+c01*Ic)*first_term[0][1]-c01*(dC[0][1]-1/3.*trace_dC2*in_Ci[0][1]) );
	//			S[0][2]=-2*1/pow(-J,2/3.)*( (c10+c01*Ic)*first_term[0][2]-c01*(dC[0][2]-1/3.*trace_dC2*in_Ci[0][2]) );
	//			S[1][0]=-2*1/pow(-J,2/3.)*( (c10+c01*Ic)*first_term[1][0]-c01*(dC[1][0]-1/3.*trace_dC2*in_Ci[1][0]) );
	//			S[1][1]=-2*1/pow(-J,2/3.)*( (c10+c01*Ic)*first_term[1][1]-c01*(dC[1][1]-1/3.*trace_dC2*in_Ci[1][1]) );
	//			S[1][2]=-2*1/pow(-J,2/3.)*( (c10+c01*Ic)*first_term[1][2]-c01*(dC[1][2]-1/3.*trace_dC2*in_Ci[1][2]) );
	//			S[2][0]=-2*1/pow(-J,2/3.)*( (c10+c01*Ic)*first_term[2][0]-c01*(dC[2][0]-1/3.*trace_dC2*in_Ci[2][0]) );
	//			S[2][1]=-2*1/pow(-J,2/3.)*( (c10+c01*Ic)*first_term[2][1]-c01*(dC[2][1]-1/3.*trace_dC2*in_Ci[2][1]) );		
	//			S[2][2]=-2*1/pow(-J,2/3.)*( (c10+c01*Ic)*first_term[2][2]-c01*(dC[2][2]-1/3.*trace_dC2*in_Ci[2][2]) );
	//		}
	//		else
	//		{
	//			S[0][0]=2*1/pow(J,2/3.)*( (c10+c01*Ic)*first_term[0][0]-c01*(dC[0][0]-1/3.*trace_dC2*in_Ci[0][0]) );
	//			S[0][1]=2*1/pow(J,2/3.)*( (c10+c01*Ic)*first_term[0][1]-c01*(dC[0][1]-1/3.*trace_dC2*in_Ci[0][1]) );
	//			S[0][2]=2*1/pow(J,2/3.)*( (c10+c01*Ic)*first_term[0][2]-c01*(dC[0][2]-1/3.*trace_dC2*in_Ci[0][2]) );
	//			S[1][0]=2*1/pow(J,2/3.)*( (c10+c01*Ic)*first_term[1][0]-c01*(dC[1][0]-1/3.*trace_dC2*in_Ci[1][0]) );
	//			S[1][1]=2*1/pow(J,2/3.)*( (c10+c01*Ic)*first_term[1][1]-c01*(dC[1][1]-1/3.*trace_dC2*in_Ci[1][1]) );
	//			S[1][2]=2*1/pow(J,2/3.)*( (c10+c01*Ic)*first_term[1][2]-c01*(dC[1][2]-1/3.*trace_dC2*in_Ci[1][2]) );
	//			S[2][0]=2*1/pow(J,2/3.)*( (c10+c01*Ic)*first_term[2][0]-c01*(dC[2][0]-1/3.*trace_dC2*in_Ci[2][0]) );
	//			S[2][1]=2*1/pow(J,2/3.)*( (c10+c01*Ic)*first_term[2][1]-c01*(dC[2][1]-1/3.*trace_dC2*in_Ci[2][1]) );		
	//			S[2][2]=2*1/pow(J,2/3.)*( (c10+c01*Ic)*first_term[2][2]-c01*(dC[2][2]-1/3.*trace_dC2*in_Ci[2][2]) );
	//		}
	//		//cout<<"S"<<j<<"="<<S[0][0]<<","<<S[0][1]<<","<<S[0][2]<<endl;
	//		//cout<<S[1][0]<<","<<S[1][1]<<","<<S[1][2]<<endl;
	//		//cout<<S[2][0]<<","<<S[2][1]<<","<<S[2][2]<<endl;
	//
	//
	//
	//		HYPER[j].pi[0][0]=HYPER[j].Fi[0][0]*S[0][0]+HYPER[j].Fi[0][1]*S[1][0]+HYPER[j].Fi[0][2]*S[2][0];
	//		HYPER[j].pi[0][1]=HYPER[j].Fi[0][0]*S[0][1]+HYPER[j].Fi[0][1]*S[1][1]+HYPER[j].Fi[0][2]*S[2][1];
	//		HYPER[j].pi[0][2]=HYPER[j].Fi[0][0]*S[0][2]+HYPER[j].Fi[0][1]*S[1][2]+HYPER[j].Fi[0][2]*S[2][2];
	//		HYPER[j].pi[1][0]=HYPER[j].Fi[1][0]*S[0][0]+HYPER[j].Fi[1][1]*S[1][0]+HYPER[j].Fi[1][2]*S[2][0];
	//		HYPER[j].pi[1][1]=HYPER[j].Fi[1][0]*S[0][1]+HYPER[j].Fi[1][1]*S[1][1]+HYPER[j].Fi[1][2]*S[2][1];
	//		HYPER[j].pi[1][2]=HYPER[j].Fi[1][0]*S[0][2]+HYPER[j].Fi[1][1]*S[1][2]+HYPER[j].Fi[1][2]*S[2][2];
	//		HYPER[j].pi[2][0]=HYPER[j].Fi[2][0]*S[0][0]+HYPER[j].Fi[2][1]*S[1][0]+HYPER[j].Fi[2][2]*S[2][0];
	//		HYPER[j].pi[2][1]=HYPER[j].Fi[2][0]*S[0][1]+HYPER[j].Fi[2][1]*S[1][1]+HYPER[j].Fi[2][2]*S[2][1];
	//		HYPER[j].pi[2][2]=HYPER[j].Fi[2][0]*S[0][2]+HYPER[j].Fi[2][1]*S[1][2]+HYPER[j].Fi[2][2]*S[2][2];
	//	}
	//	for(int D=0;D<DIMENSION;D++)	delete[]	in_Ci[D];
	//	delete[]	in_Ci;
	//}
}



void calc_stress(mpsconfig &CON,vector<hyperelastic> &HYPER)
{

#ifdef _OPENMP
	omp_set_num_threads(8);
#pragma omp parallel for
#endif
	{
		//	cout<<"応力計算";
		//y = 1.61632525555788E-04x6 - 1.55947675768484E-02x5 + 6.04913409143031E-01x4 - 1.21769666879106E+01x3 + 1.37392898012643E+02x2 - 8.80879613720886E+02x + 4.42829780193931E+03
		//double c10=0.;
		//c10=CON.get_c10();//4621.012789;//4484.5883;	//加嶋さんの実験値

		int h_num=HYPER.size();

		double **d_Fi=new double *[DIMENSION];
		for(int D=0;D<DIMENSION;D++)	d_Fi[D]=new double [DIMENSION];

		double c10=CON.get_c10();
		double c01=CON.get_c01();

		double b[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
		double bb[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
		double dC[DIMENSION]={0,0,0};

#ifdef _OPENMP
#pragma omp for
#endif
		for(int j=0;j<h_num;j++)
		{	
			double J=HYPER[j].J;
			if(J<0){
				d_Fi[0][0]=-1/pow(-J,1/3.)*HYPER[j].Fi[0][0];	d_Fi[0][1]=-1/pow(-J,1/3.)*HYPER[j].Fi[0][1];	d_Fi[0][2]=-1/pow(-J,1/3.)*HYPER[j].Fi[0][2];
				d_Fi[1][0]=-1/pow(-J,1/3.)*HYPER[j].Fi[1][0];	d_Fi[1][1]=-1/pow(-J,1/3.)*HYPER[j].Fi[1][1];	d_Fi[1][2]=-1/pow(-J,1/3.)*HYPER[j].Fi[1][2];
				d_Fi[2][0]=-1/pow(-J,1/3.)*HYPER[j].Fi[2][0];	d_Fi[2][1]=-1/pow(-J,1/3.)*HYPER[j].Fi[2][1];	d_Fi[2][2]=-1/pow(-J,1/3.)*HYPER[j].Fi[2][2];
			}
			else
			{
				d_Fi[0][0]=1/pow(J,1/3.)*HYPER[j].Fi[0][0];	d_Fi[0][1]=1/pow(J,1/3.)*HYPER[j].Fi[0][1];	d_Fi[0][2]=1/pow(J,1/3.)*HYPER[j].Fi[0][2];
				d_Fi[1][0]=1/pow(J,1/3.)*HYPER[j].Fi[1][0];	d_Fi[1][1]=1/pow(J,1/3.)*HYPER[j].Fi[1][1];	d_Fi[1][2]=1/pow(J,1/3.)*HYPER[j].Fi[1][2];
				d_Fi[2][0]=1/pow(J,1/3.)*HYPER[j].Fi[2][0];	d_Fi[2][1]=1/pow(J,1/3.)*HYPER[j].Fi[2][1];	d_Fi[2][2]=1/pow(J,1/3.)*HYPER[j].Fi[2][2];
			}

			dC[0]=d_Fi[0][0]*d_Fi[0][0]+d_Fi[1][0]*d_Fi[1][0]+d_Fi[2][0]*d_Fi[2][0];
			dC[1]=d_Fi[0][1]*d_Fi[0][1]+d_Fi[1][1]*d_Fi[1][1]+d_Fi[2][1]*d_Fi[2][1];
			dC[2]=d_Fi[0][2]*d_Fi[0][2]+d_Fi[1][2]*d_Fi[1][2]+d_Fi[2][2]*d_Fi[2][2];


			b[0][0]=d_Fi[0][0]*d_Fi[0][0]+d_Fi[0][1]*d_Fi[0][1]+d_Fi[0][2]*d_Fi[0][2];
			b[0][1]=d_Fi[0][0]*d_Fi[1][0]+d_Fi[0][1]*d_Fi[1][1]+d_Fi[0][2]*d_Fi[1][2];
			b[0][2]=d_Fi[0][0]*d_Fi[2][0]+d_Fi[0][1]*d_Fi[2][1]+d_Fi[0][2]*d_Fi[2][2];
			b[1][0]=d_Fi[1][0]*d_Fi[0][0]+d_Fi[1][1]*d_Fi[0][1]+d_Fi[1][2]*d_Fi[0][2];
			b[1][1]=d_Fi[1][0]*d_Fi[1][0]+d_Fi[1][1]*d_Fi[1][1]+d_Fi[1][2]*d_Fi[1][2];
			b[1][2]=d_Fi[1][0]*d_Fi[2][0]+d_Fi[1][1]*d_Fi[2][1]+d_Fi[1][2]*d_Fi[2][2];
			b[2][0]=d_Fi[2][0]*d_Fi[0][0]+d_Fi[2][1]*d_Fi[0][1]+d_Fi[2][2]*d_Fi[0][2];
			b[2][1]=d_Fi[2][0]*d_Fi[1][0]+d_Fi[2][1]*d_Fi[1][1]+d_Fi[2][2]*d_Fi[1][2];
			b[2][2]=d_Fi[2][0]*d_Fi[2][0]+d_Fi[2][1]*d_Fi[2][1]+d_Fi[2][2]*d_Fi[2][2];

			double trace_b=b[0][0]+b[1][1]+b[2][2];
			double trace_dC=dC[0]+dC[1]+dC[2];
			//c10=(1.61632525555788e-4)*pow(trace_dC,6.)-(1.55947675768484e-2)*pow(trace_dC,5.)+0.604913409143031*pow(trace_dC,4.)-12.1769666879106*pow(trace_dC,3.)+137.392898012643*pow(trace_dC,2.)-880.879613720886*trace_dC + 4428.29780193931;

			HYPER[j].stress[0][0]=1./J*c10*(b[0][0]-1.0/3.0*trace_b);
			HYPER[j].stress[0][1]=1./J*c10*b[0][1];
			HYPER[j].stress[0][2]=1./J*c10*b[0][2];
			HYPER[j].stress[1][0]=1./J*c10*b[1][0];
			HYPER[j].stress[1][1]=1./J*c10*(b[1][1]-1.0/3.0*trace_b);
			HYPER[j].stress[1][2]=1./J*c10*b[1][2];
			HYPER[j].stress[2][0]=1./J*c10*b[2][0];
			HYPER[j].stress[2][1]=1./J*c10*b[2][1];
			HYPER[j].stress[2][2]=1./J*c10*(b[2][2]-1.0/3.0*trace_b);
			//double J=HYPER[j].J;
			//if(J<0){
			//	d_Fi[0][0]=-1/pow(-J,1/3.)*HYPER[j].Fi[0][0];	d_Fi[0][1]=-1/pow(-J,1/3.)*HYPER[j].Fi[0][1];	d_Fi[0][2]=-1/pow(-J,1/3.)*HYPER[j].Fi[0][2];
			//	d_Fi[1][0]=-1/pow(-J,1/3.)*HYPER[j].Fi[1][0];	d_Fi[1][1]=-1/pow(-J,1/3.)*HYPER[j].Fi[1][1];	d_Fi[1][2]=-1/pow(-J,1/3.)*HYPER[j].Fi[1][2];
			//	d_Fi[2][0]=-1/pow(-J,1/3.)*HYPER[j].Fi[2][0];	d_Fi[2][1]=-1/pow(-J,1/3.)*HYPER[j].Fi[2][1];	d_Fi[2][2]=-1/pow(-J,1/3.)*HYPER[j].Fi[2][2];
			//}
			//else
			//{
			//	d_Fi[0][0]=1/pow(J,1/3.)*HYPER[j].Fi[0][0];	d_Fi[0][1]=1/pow(J,1/3.)*HYPER[j].Fi[0][1];	d_Fi[0][2]=1/pow(J,1/3.)*HYPER[j].Fi[0][2];
			//	d_Fi[1][0]=1/pow(J,1/3.)*HYPER[j].Fi[1][0];	d_Fi[1][1]=1/pow(J,1/3.)*HYPER[j].Fi[1][1];	d_Fi[1][2]=1/pow(J,1/3.)*HYPER[j].Fi[1][2];
			//	d_Fi[2][0]=1/pow(J,1/3.)*HYPER[j].Fi[2][0];	d_Fi[2][1]=1/pow(J,1/3.)*HYPER[j].Fi[2][1];	d_Fi[2][2]=1/pow(J,1/3.)*HYPER[j].Fi[2][2];
			//}

			//b[0][0]=d_Fi[0][0]*d_Fi[0][0]+d_Fi[0][1]*d_Fi[0][1]+d_Fi[0][2]*d_Fi[0][2];
			//b[0][1]=d_Fi[0][0]*d_Fi[1][0]+d_Fi[0][1]*d_Fi[1][1]+d_Fi[0][2]*d_Fi[1][2];
			//b[0][2]=d_Fi[0][0]*d_Fi[2][0]+d_Fi[0][1]*d_Fi[2][1]+d_Fi[0][2]*d_Fi[2][2];
			//b[1][0]=d_Fi[1][0]*d_Fi[0][0]+d_Fi[1][1]*d_Fi[0][1]+d_Fi[1][2]*d_Fi[0][2];
			//b[1][1]=d_Fi[1][0]*d_Fi[1][0]+d_Fi[1][1]*d_Fi[1][1]+d_Fi[1][2]*d_Fi[1][2];
			//b[1][2]=d_Fi[1][0]*d_Fi[2][0]+d_Fi[1][1]*d_Fi[2][1]+d_Fi[1][2]*d_Fi[2][2];
			//b[2][0]=d_Fi[2][0]*d_Fi[0][0]+d_Fi[2][1]*d_Fi[0][1]+d_Fi[2][2]*d_Fi[0][2];
			//b[2][1]=d_Fi[2][0]*d_Fi[1][0]+d_Fi[2][1]*d_Fi[1][1]+d_Fi[2][2]*d_Fi[1][2];
			//b[2][2]=d_Fi[2][0]*d_Fi[2][0]+d_Fi[2][1]*d_Fi[2][1]+d_Fi[2][2]*d_Fi[2][2];

			//bb[0][0]=b[0][0]*b[0][0]+b[0][1]*b[1][0]+b[0][2]*b[2][0];
			//bb[0][1]=b[0][0]*b[0][1]+b[0][1]*b[1][1]+b[0][2]*b[2][1];
			//bb[0][2]=b[0][0]*b[0][2]+b[0][1]*b[1][2]+b[0][2]*b[2][2];
			//bb[1][0]=b[1][0]*b[0][0]+b[1][1]*b[1][0]+b[1][2]*b[2][0];
			//bb[1][1]=b[1][0]*b[0][1]+b[1][1]*b[1][1]+b[1][2]*b[2][1];
			//bb[1][2]=b[1][0]*b[0][2]+b[1][1]*b[1][2]+b[1][2]*b[2][2];
			//bb[2][0]=b[2][0]*b[0][0]+b[2][1]*b[1][0]+b[2][2]*b[2][0];
			//bb[2][1]=b[2][0]*b[0][1]+b[2][1]*b[1][1]+b[2][2]*b[2][1];
			//bb[2][2]=b[2][0]*b[0][2]+b[2][1]*b[1][2]+b[2][2]*b[2][2];

			//double trace_b=b[0][0]+b[1][1]+b[2][2];
			//double trace_bb=bb[0][0]+bb[1][1]+bb[2][2];

			//HYPER[j].stress[0][0]=2./J*((c10+c01*trace_b)*(b[0][0]-1.0/3.0*trace_b)-c01*(bb[0][0]-1.0/3.0*trace_bb));
			//HYPER[j].stress[0][1]=2./J*((c10+c01*trace_b)*b[0][1]-c01*bb[0][1]);
			//HYPER[j].stress[0][2]=2./J*((c10+c01*trace_b)*b[0][2]-c01*bb[0][2]);
			//HYPER[j].stress[1][0]=2./J*((c10+c01*trace_b)*b[1][0]-c01*bb[1][0]);
			//HYPER[j].stress[1][1]=2./J*((c10+c01*trace_b)*(b[1][1]-1.0/3.0*trace_b)-c01*(bb[1][1]-1.0/3.0*trace_bb));
			//HYPER[j].stress[1][2]=2./J*((c10+c01*trace_b)*b[1][2]-c01*bb[1][2]);
			//HYPER[j].stress[2][0]=2./J*((c10+c01*trace_b)*b[2][0]-c01*bb[2][0]);
			//HYPER[j].stress[2][1]=2./J*((c10+c01*trace_b)*b[2][1]-c01*bb[2][1]);
			//HYPER[j].stress[2][2]=2./J*((c10+c01*trace_b)*(b[2][2]-1.0/3.0*trace_b)-c01*(bb[2][2]-1.0/3.0*trace_bb));
		}

		for(int D=0;D>DIMENSION;D++)	delete[]	d_Fi[D];
		delete[]	d_Fi;

		//	cout<<"----------OK"<<endl;
	}
}

void calc_differential_p(mpsconfig &CON,vector<mpselastic>PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,double **F)
{
	//	cout<<"運動量微分値計算";

#ifdef _OPENMP
	omp_set_num_threads(8);
#pragma omp parallel for
#endif
	{
		int h_num=HYPER.size();
		int flag_vis=CON.get_flag_vis();
		int flag_G=CON.get_flag_G();
		int flag_FEM=CON.get_FEM_flag();
		double Dt=CON.get_dt();
		double V=get_volume(&CON);
		double mh=V*CON.get_hyper_density();
		double ms=V*CON.get_silicone_density();
		double density=CON.get_hyper_density();
		double G=9.8;

#ifdef _OPENMP
#pragma omp for
#endif
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

			//重力影響
			//if(flag_G==ON)
			//{
			/*if(PART[i].toFEM==ON)	*/HYPER[i].differential_p[A_Z]-=Dt*0.5*G*mh;
			//else
			//{
			//	HYPER[i].differential_p[A_Z]-=Dt*0.5*G*ms;
			//}
			//if(flag_vis==ON)
			{
				HYPER[i].differential_p[A_X]+=Dt*0.5*HYPER[i].vis_force[A_X];
				HYPER[i].differential_p[A_Y]+=Dt*0.5*HYPER[i].vis_force[A_Y];
				HYPER[i].differential_p[A_Z]+=Dt*0.5*HYPER[i].vis_force[A_Z];
			}
			//if(flag_FEM==ON && PART[i].toFEM==ON)
			{
				HYPER[i].differential_p[A_X]+=Dt*0.5*F[A_X][i]*V;//density;
				HYPER[i].differential_p[A_Y]+=Dt*0.5*F[A_Y][i]*V;//density;
				HYPER[i].differential_p[A_Z]+=Dt*0.5*F[A_Z][i]*V;//density;
			}
			//}
			//	cout<<"----------OK"<<endl;
		}
	}
}
void renew_lambda(mpsconfig &CON,vector<mpselastic>PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int t)
{
	clock_t t3=clock();

#ifdef _OPENMP
	omp_set_num_threads(8);
#pragma omp parallel for
#endif
	{
		//	cout<<"Lambda計算";

		int h_num=HYPER.size();
		double le=CON.get_distancebp();
		double Dt=CON.get_dt();
		double V=get_volume(&CON);
		double mh=V*CON.get_hyper_density();
		double ms=V*CON.get_silicone_density();

		double *N_Left=new double[h_num*h_num];
		double *N_Right=new double[h_num];
#ifdef _OPENMP
#pragma omp for
#endif
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
		}//jに関するfor文の終わり
		}//iに関するfor文の終わり
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

#ifdef _OPENMP
#pragma omp for
#endif
		for(int k=0;k<h_num;k++)
		{
			//	if(HYPER[k].flag_wall==OFF)
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
				}//jに関するfor文の終わり
				N_Right[k]=N_right+HYPER[k].differential_p[0]*HYPER1[k*h_num+k].DgDq[0]+HYPER[k].differential_p[1]*HYPER1[k*h_num+k].DgDq[1]+HYPER[k].differential_p[2]*HYPER1[k*h_num+k].DgDq[2];//1/mk*N_right;	
				N_Left[k*h_num+k]+=Dt*0.5*(HYPER1[k*h_num+k].DgDq[0]*HYPER1[k*h_num+k].DgDq[0]+HYPER1[k*h_num+k].DgDq[1]*HYPER1[k*h_num+k].DgDq[1]+HYPER1[k*h_num+k].DgDq[2]*HYPER1[k*h_num+k].DgDq[2]);
			}
			/*		if(HYPER[k].flag_wall==ON)
			{
			int Nk=HYPER[k].N;
			double N_right=0;
			for(int i=0;i<Nk;i++)
			{
			int in=HYPER[k].NEI[i];
			for(int j=0;j<Nk;j++)
			{
			int jn=HYPER[k].NEI[j];
			N_Left[in*h_num+jn]+=Dt*0.5*(HYPER1[in*h_num+k].DgDq[0]*HYPER1[jn*h_num+k].DgDq[0]+HYPER1[in*h_num+k].DgDq[1]*HYPER1[jn*h_num+k].DgDq[1]-HYPER1[in*h_num+k].DgDq[2]*HYPER1[jn*h_num+k].DgDq[2]);
			}
			N_Left[in*h_num+k]+=Dt*0.5*(HYPER1[in*h_num+k].DgDq[0]*HYPER1[k*h_num+k].DgDq[0]+HYPER1[in*h_num+k].DgDq[1]*HYPER1[k*h_num+k].DgDq[1]-HYPER1[in*h_num+k].DgDq[2]*HYPER1[k*h_num+k].DgDq[2]);
			N_Left[k*h_num+in]+=Dt*0.5*(HYPER1[k*h_num+k].DgDq[0]*HYPER1[in*h_num+k].DgDq[0]+HYPER1[k*h_num+k].DgDq[1]*HYPER1[in*h_num+k].DgDq[1]-HYPER1[k*h_num+k].DgDq[2]*HYPER1[in*h_num+k].DgDq[2]);
			N_right+=HYPER[in].differential_p[0]*HYPER1[k*h_num+in].DgDq[0]+HYPER[in].differential_p[1]*HYPER1[k*h_num+in].DgDq[1]-HYPER[in].differential_p[2]*HYPER1[k*h_num+in].DgDq[2];
			}//jに関するfor文の終わり
			N_Right[k]=N_right+HYPER[k].differential_p[0]*HYPER1[k*h_num+k].DgDq[0]+HYPER[k].differential_p[1]*HYPER1[k*h_num+k].DgDq[1]-HYPER[k].differential_p[2]*HYPER1[k*h_num+k].DgDq[2];//1/mk*N_right;	
			N_Left[k*h_num+k]+=Dt*0.5*(HYPER1[k*h_num+k].DgDq[0]*HYPER1[k*h_num+k].DgDq[0]+HYPER1[k*h_num+k].DgDq[1]*HYPER1[k*h_num+k].DgDq[1]-HYPER1[k*h_num+k].DgDq[2]*HYPER1[k*h_num+k].DgDq[2]);
			}//*/
		}//iに関するfor文の終わり


		/*
		for(int i=0;i<h_num;i++)	cout<<"N_Right["<<i<<"]="<<N_Right[i]<<endl;
		cout<<endl;
		for(int i=0;i<h_num;i++)
		{
		for(int j=0;j<h_num;j++)	cout<<"N_Left["<<i<<"]["<<j<<"]="<<N_Left[i*h_num+j]<<endl;
		cout<<endl;
		}*/

		//lambdaを求める

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


		/*	for(int i=0;i<h_num;i++)
		{
		if(PART[i].r[A_Z]<0)
		{
		int Ni=HYPER[i].N;
		double Np_Right_part=0;
		for(int j=0;j<Ni;j++)
		{
		int jn=HYPER[i].NEI[j];
		N_Left[i*h_num+jn]=HYPER1[jn*h_num+i].DgDq[A_Z]*Dt/2;
		Np_Right_part+=Dt/2*(HYPER[jn].stress[A_Z][A_X]*HYPER1[jn*h_num+i].DgDq[A_X]+HYPER[jn].stress[A_Z][A_Y]*HYPER1[jn*h_num+i].DgDq[A_Y]+HYPER[jn].stress[A_Z][A_Z]*HYPER1[jn*h_num+i].DgDq[A_Z]);
		}
		N_Right[i]=Np_Right_part+HYPER[i].half_p[A_Z];
		}
		}//*/
		gauss(N_Left,N_Right,h_num);
		//double ep=CON.get_FEMCGep();
		//GaussSeidelvh(N_Left,h_num,N_Right,ep);
#ifdef _OPENMP
#pragma omp for
#endif
		for(int i=0;i<h_num;i++)	HYPER[i].lambda=N_Right[i];//*/

		delete [] N_Left;
		delete [] N_Right;

		/*
		int all_ind_num=0;

		int *b_ind=new int [h_num*h_num];
		double *b_val=new double [h_num*h_num];
		int *ptr=new int [h_num+1];
		double *X=new double [h_num];

		for(int i=0;i<h_num;i++)
		{
		ptr[i+1]=0;
		X[i]=0;
		for(int j=0;j<h_num;j++)
		{
		b_ind[i*h_num+j]=0;
		b_val[i*h_num+j]=0;
		}
		}

		ptr[0]=0;
		for(int i=0;i<h_num;i++)
		{
		int ind_num=0;
		for(int j=0;j<h_num;j++)
		{
		if(N_Left[i*h_num+j]!=0)
		{
		b_ind[all_ind_num+ind_num]=j;	
		b_val[all_ind_num+ind_num]=N_Left[i*h_num+j];
		ind_num++;
		ptr[i+1]=ind_num+ptr[i];
		}			
		}
		//cout<<"ptr"<<i+1<<" "<<ptr[i+1]<<endl;
		all_ind_num+=ind_num;			
		}


		int *ind=new int [all_ind_num];
		double *val=new double [all_ind_num];
		for(int i=0;i<all_ind_num;i++)
		{
		ind[i]=b_ind[i];
		val[i]=b_val[i];
		}
		delete[] b_ind;
		delete[] b_val;

		BiCGStab2_method(&CON,val,ind,ptr,h_num,N_Right,all_ind_num,X);
		//iccg2(&CON,val,ind,ptr,h_num,N_Right,all_ind_num,X);
		//CG3D(&CON,val,ind,ptr,h_num,N_Right,all_ind_num,X);

		delete[] ind;
		delete[] val;
		delete[] ptr;

		for(int i=0;i<h_num;i++)	HYPER[i].lambda=X[i];

		delete[] X;//*/

		//	for(int i=0;i<h_num;i++)	cout<<"lambda["<<i<<"]="<<HYPER[i].lambda<<endl;

		ofstream t_loge("time_log.dat", ios::app);
		clock_t t4=clock();
		t_loge<<"step="<<t<<", time="<<1000*(t4-t3)/CLOCKS_PER_SEC<<"[e-3sec]"<<endl;
		t_loge.close();
	}

	//	cout<<"----------OK"<<endl;
}

void previous_strage(vector<mpselastic>PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int h_num)
{
#ifdef _OPENMP
	omp_set_num_threads(8);
#pragma omp parallel for
#endif
	{
#ifdef _OPENMP
#pragma omp for
#endif
		for(int i=0;i<h_num;i++)
		{
			HYPER[i].q_n[A_X]=PART[i].r[A_X];	HYPER[i].q_n[A_Y]=PART[i].r[A_Y];	HYPER[i].q_n[A_Z]=PART[i].r[A_Z];
			HYPER[i].p_n[A_X]=HYPER[i].p[A_X];	HYPER[i].p_n[A_Y]=HYPER[i].p[A_Y];	HYPER[i].p_n[A_Z]=HYPER[i].p[A_Z];
			HYPER[i].ph_n[A_X]=HYPER[i].half_p[A_X];	HYPER[i].ph_n[A_Y]=HYPER[i].half_p[A_Y];	HYPER[i].ph_n[A_Z]=HYPER[i].half_p[A_Z];
			HYPER[i].W_n=HYPER[i].W;
			HYPER[i].stress_n[A_X][A_X]=HYPER[i].stress[A_X][A_X];	HYPER[i].stress_n[A_X][A_Y]=HYPER[i].stress[A_X][A_Y];	HYPER[i].stress_n[A_X][A_Z]=HYPER[i].stress[A_X][A_Z];
			HYPER[i].stress_n[A_Y][A_X]=HYPER[i].stress[A_Y][A_X];	HYPER[i].stress_n[A_Y][A_Y]=HYPER[i].stress[A_Y][A_Y];	HYPER[i].stress_n[A_Y][A_Z]=HYPER[i].stress[A_Y][A_Z];
			HYPER[i].stress_n[A_Z][A_X]=HYPER[i].stress[A_Z][A_X];	HYPER[i].stress_n[A_Z][A_Y]=HYPER[i].stress[A_Z][A_Y];	HYPER[i].stress_n[A_Z][A_Z]=HYPER[i].stress[A_Z][A_Z];
			//HYPER[i].pi_n[A_X][A_X]=HYPER[i].pi[A_X][A_X];	HYPER[i].pi_n[A_X][A_Y]=HYPER[i].pi[A_X][A_Y];	HYPER[i].pi_n[A_X][A_Z]=HYPER[i].pi[A_X][A_Z];
			//HYPER[i].pi_n[A_Y][A_X]=HYPER[i].pi[A_Y][A_X];	HYPER[i].pi_n[A_Y][A_Y]=HYPER[i].pi[A_Y][A_Y];	HYPER[i].pi_n[A_Y][A_Z]=HYPER[i].pi[A_Y][A_Z];
			//HYPER[i].pi_n[A_Z][A_X]=HYPER[i].pi[A_Z][A_X];	HYPER[i].pi_n[A_Z][A_Y]=HYPER[i].pi[A_Z][A_Y];	HYPER[i].pi_n[A_Z][A_Z]=HYPER[i].pi[A_Z][A_Z];

			int Ni=HYPER[i].N;
			for(int j=0;j<Ni;j++)
			{
				int k=HYPER[i].NEI[j];
				HYPER1[k*h_num+i].DgDq_n[A_X]=HYPER1[k*h_num+i].DgDq[A_X];	HYPER1[k*h_num+i].DgDq_n[A_Y]=HYPER1[k*h_num+i].DgDq[A_Y];	HYPER1[k*h_num+i].DgDq_n[A_Z]=HYPER1[k*h_num+i].DgDq[A_Z];
			}
			HYPER1[i*h_num+i].DgDq_n[A_X]=HYPER1[i*h_num+i].DgDq[A_X];	HYPER1[i*h_num+i].DgDq_n[A_Y]=HYPER1[i*h_num+i].DgDq[A_Y];	HYPER1[i*h_num+i].DgDq_n[A_Z]=HYPER1[i*h_num+i].DgDq[A_Z];
		}
	}
}

void calc_W(mpsconfig &CON,vector<hyperelastic> &HYPER,int h_num)
{
	//	cout<<"弾性ポテンシャル計算";

	double d_Fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	//double c10=0.;
	//c10=CON.get_c10();	//4621.012789;//4484.5883;	//加嶋さんの実験値
	double c10=CON.get_c10();
	double c01=CON.get_c01();



	double dC[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double dC2[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

#ifdef _OPENMP
	omp_set_num_threads(8);
#pragma omp parallel for
#endif
	{

#ifdef _OPENMP
#pragma omp for
#endif
		for(int i=0;i<h_num;i++)
		{
			//double J=HYPER[i].J;
			//if(J<0){
			//	d_Fi[0][0]=-1/pow(-J,1/3)*HYPER[i].Fi[0][0];	d_Fi[0][1]=-1/pow(-J,1/3)*HYPER[i].Fi[0][1];	d_Fi[0][2]=-1/pow(-J,1/3)*HYPER[i].Fi[0][2];
			//	d_Fi[1][0]=-1/pow(-J,1/3)*HYPER[i].Fi[1][0];	d_Fi[1][1]=-1/pow(-J,1/3)*HYPER[i].Fi[1][1];	d_Fi[1][2]=-1/pow(-J,1/3)*HYPER[i].Fi[1][2];
			//	d_Fi[2][0]=-1/pow(-J,1/3)*HYPER[i].Fi[2][0];	d_Fi[2][1]=-1/pow(-J,1/3)*HYPER[i].Fi[2][1];	d_Fi[2][2]=-1/pow(-J,1/3)*HYPER[i].Fi[2][2];
			//}
			//else
			//{
			//	d_Fi[0][0]=1/pow(J,1/3)*HYPER[i].Fi[0][0];	d_Fi[0][1]=1/pow(J,1/3)*HYPER[i].Fi[0][1];	d_Fi[0][2]=1/pow(J,1/3)*HYPER[i].Fi[0][2];
			//	d_Fi[1][0]=1/pow(J,1/3)*HYPER[i].Fi[1][0];	d_Fi[1][1]=1/pow(J,1/3)*HYPER[i].Fi[1][1];	d_Fi[1][2]=1/pow(J,1/3)*HYPER[i].Fi[1][2];
			//	d_Fi[2][0]=1/pow(J,1/3)*HYPER[i].Fi[2][0];	d_Fi[2][1]=1/pow(J,1/3)*HYPER[i].Fi[2][1];	d_Fi[2][2]=1/pow(J,1/3)*HYPER[i].Fi[2][2];
			//}

			//dC[0][0]=d_Fi[0][0]*d_Fi[0][0]+d_Fi[1][0]*d_Fi[1][0]+d_Fi[2][0]*d_Fi[2][0];
			//dC[0][1]=d_Fi[0][0]*d_Fi[0][1]+d_Fi[1][0]*d_Fi[1][1]+d_Fi[2][0]*d_Fi[2][1];
			//dC[0][2]=d_Fi[0][0]*d_Fi[0][2]+d_Fi[1][0]*d_Fi[1][2]+d_Fi[2][0]*d_Fi[2][2];
			//dC[1][0]=d_Fi[0][1]*d_Fi[0][0]+d_Fi[1][1]*d_Fi[1][0]+d_Fi[2][1]*d_Fi[2][0];
			//dC[1][1]=d_Fi[0][1]*d_Fi[0][1]+d_Fi[1][1]*d_Fi[1][1]+d_Fi[2][1]*d_Fi[2][1];
			//dC[1][2]=d_Fi[0][1]*d_Fi[0][2]+d_Fi[1][1]*d_Fi[1][2]+d_Fi[2][1]*d_Fi[2][2];
			//dC[2][0]=d_Fi[0][2]*d_Fi[0][0]+d_Fi[1][2]*d_Fi[1][0]+d_Fi[2][2]*d_Fi[2][0];
			//dC[2][1]=d_Fi[0][2]*d_Fi[0][1]+d_Fi[1][2]*d_Fi[1][1]+d_Fi[2][2]*d_Fi[2][1];
			//dC[2][2]=d_Fi[0][2]*d_Fi[0][2]+d_Fi[1][2]*d_Fi[1][2]+d_Fi[2][2]*d_Fi[2][2];

			//dC2[0][0]=dC[A_X][0]*dC[0][A_X]+dC[A_X][1]*dC[1][A_X]+dC[A_X][2]*dC[2][A_X];
			//dC2[0][1]=dC[A_X][0]*dC[0][A_Y]+dC[A_X][1]*dC[1][A_Y]+dC[A_X][2]*dC[2][A_Y];
			//dC2[0][2]=dC[A_X][0]*dC[0][A_Z]+dC[A_X][1]*dC[1][A_Z]+dC[A_X][2]*dC[2][A_Z];
			//dC2[1][0]=dC[A_Y][0]*dC[0][A_X]+dC[A_Y][1]*dC[1][A_X]+dC[A_Y][2]*dC[2][A_X];
			//dC2[1][1]=dC[A_Y][0]*dC[0][A_Y]+dC[A_Y][1]*dC[1][A_Y]+dC[A_Y][2]*dC[2][A_Y];
			//dC2[1][2]=dC[A_Y][0]*dC[0][A_Z]+dC[A_Y][1]*dC[1][A_Z]+dC[A_Y][2]*dC[2][A_Z];
			//dC2[2][0]=dC[A_Z][0]*dC[0][A_X]+dC[A_Z][1]*dC[1][A_X]+dC[A_Z][2]*dC[2][A_X];
			//dC2[2][1]=dC[A_Z][0]*dC[0][A_Y]+dC[A_Z][1]*dC[1][A_Y]+dC[A_Z][2]*dC[2][A_Y];
			//dC2[2][2]=dC[A_Z][0]*dC[0][A_Z]+dC[A_Z][1]*dC[1][A_Z]+dC[A_Z][2]*dC[2][A_Z];

			//double trace_dC=dC[0][0]+dC[1][1]+dC[2][2];
			//double trace_dC2=dC2[0][0]+dC2[1][1]+dC2[2][2];

			//double Ic=trace_dC;
			//double IIc=0.50*(trace_dC*trace_dC-trace_dC2);
			//HYPER[i].W=c10*(Ic-3)+c01*(IIc-3);
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

			double trace_dC=dC[0][0]+dC[1][1]+dC[2][2];

			//c10=(1.61632525555788e-4)*pow(trace_dC,6.)-(1.55947675768484e-2)*pow(trace_dC,5.)+0.604913409143031*pow(trace_dC,4.)-12.1769666879106*pow(trace_dC,3.)+137.392898012643*pow(trace_dC,2.)-880.879613720886*trace_dC + 4428.29780193931;
			double Ic=trace_dC;
			HYPER[i].W=0.5*c10*(Ic-3);
		}
		//	cout<<"----------OK"<<endl;
	}
}

//圧力など粒子の持つ情報をコンター図で表示する関数
void momentum_movie_AVS(mpsconfig &CON,int t,vector<mpselastic> PART,vector<hyperelastic> HYPER,double **F)
{
	//参考にしている書式はmicroAVSのヘルプであなたのデータは？→「非構造格子型データ（アスキー）の書式」
	int h_num=HYPER.size();
	double TIME=CON.get_step()*CON.get_dt();
	double le=CON.get_distancebp();
	int STEP=CON.get_step()/CON.get_interval()+1;		//出力する総ステップ数
	int step;
	int n=3;	//出力次元
	double mh=get_volume(&CON)*CON.get_hyper_density();
	double ms=get_volume(&CON)*CON.get_silicone_density();

	if(CON.get_interval()==1)	step=t/CON.get_interval();
	else step=t/CON.get_interval()+1;


	if(t==1) 
	{
		ofstream fp("momentum.inp", ios::trunc);			
		fp<<STEP<<endl;//総ステップ数
		fp<<"data_geom"<<endl;
		fp.close();
	}

	//mainファイル書き込み
	ofstream fp("momentum.inp",ios :: app);
	fp<<"step"<<step<<" TIME="<<TIME<<endl;

	//fp<<step<<endl;

	//fp<<"data_geom"<<endl;
	//fp<<"step1"<<endl;
	//fp<<"step"<<t/CON->get_interval()+1<<" TIME="<<TIME<<endl;

	//if(n==3)
	//{
	fp<<h_num<<" "<<h_num<<endl;	//節点数と要素数出力

	//節点番号とその座標の出力 
	for(int i=0;i<h_num;i++) fp<<i<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;

	//要素番号と要素形状の種類、そして要素を構成する節点番号出力
	for(int i=0;i<h_num;i++)	fp<<i<<"  0 pt "<<i<<endl;

	//fp<<"2 3"<<endl;//節点の情報量が2で、要素の情報量が3ということ。
	fp<<"10 0"<<endl;//節点の情報量が8で、要素の情報量が0ということ。
	fp<<"10 1 1 1 1 1 1 1 1 1 1"<<endl;	//この行の詳細はヘルプを参照
	//fp<<"8 1 1 1 1 1 1 1 1"<<endl;	//この行の詳細はヘルプを参照
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

	//各節点の情報値入力
	for(int i=0;i<h_num;i++)
	{
		fp<<i<<" "<<HYPER[i].p[A_X]<<" "<<HYPER[i].p[A_Y]<<" "<<HYPER[i].p[A_Z]<<" "<<HYPER[i].lambda<<" "<<F[A_X][i]<<" "<<F[A_Y][i]<<" "<<F[A_Z][i]<<" "<<HYPER[i].vis_force[A_X]<<" "<<HYPER[i].vis_force[A_Y]<<" "<<HYPER[i].vis_force[A_Z]<<endl;
		//fp<<i<<" "<<NODE[i].depth<<" "<<NODE[i].L/le<<" "<<NODE[i].potential<<" "<<NODE[i].Fs<<" "<<endl;
	}
	//}
	//else if(n==2)
	//{
	//	int num=0;
	//	int *Nxz=new int [h_num];
	//	fp<<h_num<<" "<<h_num<<endl;	//節点数と要素数出力

	//	//節点番号とその座標の出力 
	//	for(int i=0;i<h_num;i++)
	//	{
	//		if(PART[i].q0[A_Y]<1/2*le&&PART[i].q0[A_Y]>-1/2*le)
	//		{
	//			fp<<i<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
	//			Nxz[num]=i;
	//			num++;
	//		}
	//	}
	//	//要素番号と要素形状の種類、そして要素を構成する節点番号出力
	//	for(int i=0;i<num;i++)	fp<<i<<"  0 pt "<<i<<endl;

	//	//fp<<"2 3"<<endl;//節点の情報量が2で、要素の情報量が3ということ。
	//	//		fp<<"10 0"<<endl;//節点の情報量が8で、要素の情報量が0ということ。
	//	//		fp<<"10 1 1 1 1 1 1 1 1 1 1"<<endl;	//この行の詳細はヘルプを参照
	//	//fp<<"8 1 1 1 1 1 1 1 1"<<endl;	//この行の詳細はヘルプを参照
	//	fp<<"7 0"<<endl;
	//	fp<<"7 1 1 1 1 1 1 1"<<endl;
	//	fp<<"v_x,"<<endl;
	//	fp<<"v_y,"<<endl;
	//	fp<<"v_z,"<<endl;
	//	fp<<"lambda,"<<endl;
	//	fp<<"F_x,"<<endl;
	//	fp<<"F_y,"<<endl;
	//	fp<<"F_z,"<<endl;
	//	/*		fp<<"v_x,"<<endl;
	//	fp<<"v_y,"<<endl;
	//	fp<<"v_z,"<<endl;*/
	//	//fp<<"P,N/m^2"<<endl;
	//	//fp<<"value1,??"<<endl;

	//	//各節点の情報値入力
	//	for(int i=0;i<num;i++)
	//	{
	//		int j=Nxz[i];
	//		if(PART[i].toFEM==ON)
	//		{
	//			fp<<j<<" "<<1/mh*HYPER[j].p[A_X]<<" "<<1/mh*HYPER[j].p[A_Y]<<" "<<1/mh*HYPER[j].p[A_Z]<<" "<<HYPER[j].lambda<<" "<<F[A_X][j]<<" "<<F[A_Y][j]<<" "<<F[A_Z][j]<<endl;
	//		}
	//		else
	//		{
	//			fp<<j<<" "<<1/ms*HYPER[j].p[A_X]<<" "<<1/ms*HYPER[j].p[A_Y]<<" "<<1/ms*HYPER[j].p[A_Z]<<" "<<HYPER[j].lambda<<" "<<F[A_X][j]<<" "<<F[A_Y][j]<<" "<<F[A_Z][j]<<endl;

	//		}
	//		//fp<<j<<" "<<HYPER[j].p[A_X]<<" "<<HYPER[j].p[A_Y]<<" "<<HYPER[j].p[A_Z]<<" "<<HYPER[j].lambda<<" "<<F[A_X][j]<<" "<<F[A_Y][j]<<" "<<F[A_Z][j]<<" "<<HYPER[j].vjs_force[A_X]<<" "<<HYPER[j].vjs_force[A_Y]<<" "<<HYPER[j].vjs_force[A_Z]<<endl;
	//		//fp<<i<<" "<<NODE[i].depth<<" "<<NODE[i].L/le<<" "<<NODE[i].potential<<" "<<NODE[i].Fs<<" "<<endl;
	//	}
	//	delete[] Nxz;
	//}
	fp.close();
}

void output_hyper_data(vector<mpselastic> PART,vector<hyperelastic> HYPER,vector<hyperelastic2> HYPER1,int t,double V)
{
	int h_num=HYPER.size();


	ofstream j("J.csv", ios::app);
	ofstream lam("lambda.csv", ios::app);
	//	ofstream p("P.csv", ios::app);
	//ofstream p_an("P_ave_norm.csv", ios::app);
	//ofstream h("model_height.csv", ios::app);
	ofstream fg("g.csv", ios::app);

	//	if(t==1)	p<<"h_p,,,d_p,,,p\n";
	//	p<<"t"<<t<<endl;
	double sum_lam=0;
	double sum_j=0;
	double sum_g=0.;
	double p_sum=0,d_p_sum=0,h_p_sum=0;
	double h_min=PART[0].r[A_Z];
	double h_max=PART[0].r[A_Z];
	for(int i=0;i<h_num;i++)
	{
		lam<<HYPER[i].lambda<<",";
		j<<HYPER[i].J<<",";
		fg<<V*(1-HYPER[i].J)<<",";
		sum_lam+=HYPER[i].lambda;
		sum_j+=HYPER[i].J;
		sum_g+=V*(1-HYPER[i].J);
		//		p<<HYPER[i].half_p[A_X]<<","<<HYPER[i].half_p[A_Y]<<","<<HYPER[i].half_p[A_Z]<<","<<HYPER[i].differential_p[A_X]<<","<<HYPER[i].differential_p[A_Y]<<","<<HYPER[i].differential_p[A_Z]<<","<<HYPER[i].p[A_X]<<","<<HYPER[i].p[A_Y]<<","<<HYPER[i].p[A_Z]<<endl;
		//p_sum+=sqrt(HYPER[i].p[A_X]*HYPER[i].p[A_X]+HYPER[i].p[A_Y]*HYPER[i].p[A_Y]+HYPER[i].p[A_Z]*HYPER[i].p[A_Z]);
		//d_p_sum+=sqrt(HYPER[i].differential_p[A_X]*HYPER[i].differential_p[A_X]+HYPER[i].differential_p[A_Y]*HYPER[i].differential_p[A_Y]+HYPER[i].differential_p[A_Z]*HYPER[i].differential_p[A_Z]);
		//h_p_sum+=sqrt(HYPER[i].half_p[A_X]*HYPER[i].half_p[A_X]+HYPER[i].half_p[A_Y]*HYPER[i].half_p[A_Y]+HYPER[i].half_p[A_Z]*HYPER[i].half_p[A_Z]);
		//if(h_min>PART[i].r[A_Z])	h_min=PART[i].r[A_Z];
		//if(h_max<PART[i].r[A_Z])	h_max=PART[i].r[A_Z];
	}
	fg<<sum_g/h_num<<endl;
	lam<<sum_lam/h_num<<endl;
	j<<sum_j/h_num<<endl;
	//	p<<endl;
	//if(t==1)
	//{
	//	//p_an<<"p_sum,d_p_sum,h_p_sum\n";
	//	h<<"h_min,h_max\n";
	//}
	////p_an<<p_sum/h_num<<","<<d_p_sum/h_num<<","<<h_p_sum/h_num<<endl;
	//h<<h_min<<","<<h_max<<endl;

	lam.close();
	j.close();
	//	p.close();
	////p_an.close();
	//h.close();
	fg.close();


	/*if(t==1)
	{
	////計算した各定数の出力
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

	/*	ofstream d_p("d_P.csv", ios::app);
	ofstream h_p("h_P.csv", ios::app);*/
	//stringstream ss_dgdq;
	//ss_dgdq<<"./DgDq/DgDq"<<t<<".csv";
	//ofstream dg(ss_dgdq.str());
	//stringstream ss_n0;
	//ss_n0<<"./n0/n0"<<t<<".csv";
	//ofstream fn0(ss_n0.str());

	//dg<<"t,"<<t<<endl;
	//fn0<<"t,"<<t<<endl;
	//for(int i=0;i<h_num;i++)
	//{
	//	int N=HYPER[i].N;
	//	dg<<",,"<<i<<endl;
	//	fn0<<",,"<<i<<endl;
	//	for(int j=0;j<N;j++)
	//	{
	//		int nei=HYPER[i].NEI[j];
	//		dg<<",,,"<<nei<<","<<HYPER1[i*h_num+nei].DgDq[A_X]<<","<<HYPER1[i*h_num+nei].DgDq[A_Y]<<","<<HYPER1[i*h_num+nei].DgDq[A_Z]<<endl;
	//		fn0<<",,,"<<nei<<","<<HYPER1[nei*h_num+i].n0ij[A_X]<<","<<HYPER1[nei*h_num+i].n0ij[A_Y]<<","<<HYPER1[nei*h_num+i].n0ij[A_Z]<<endl;
	//	}
	//	dg<<",,,"<<i<<","<<HYPER1[i*h_num+i].DgDq[A_X]<<","<<HYPER1[i*h_num+i].DgDq[A_Y]<<","<<HYPER1[i*h_num+i].DgDq[A_Z]<<endl;
	//	fn0<<",,,"<<i<<","<<HYPER1[i*h_num+i].n0ij[A_X]<<","<<HYPER1[i*h_num+i].n0ij[A_Y]<<","<<HYPER1[i*h_num+i].n0ij[A_Z]<<endl;
	//}
	//dg.close();
	//fn0.close();



	//if(t==1)
	//{
	//	ofstream fsdg("stress_DgDq.csv",ios::trunc);
	//	ofstream fpi("pi.csv",ios::trunc);
	//	ofstream fpn("pi_n0.csv",ios::trunc);
	//	ofstream fst("stress.csv",ios::trunc);
	//	fsdg.close();
	//	fpn.close();
	//	fpi.close();
	//	fst.close();
	//}

	//ofstream fsdg("stress_DgDq.csv",ios::app);
	////ofstream fpi("pi.csv",ios::app);
	//ofstream fst("stress.csv",ios::app);
	////ofstream fpn("pi_n0.csv",ios::app);


	//fsdg<<"t,"<<t<<endl;
	////fpn<<"t,"<<t<<endl;
	////fpi<<"t,"<<t<<endl;
	//fst<<"t,"<<t<<endl;
	//double sdgdq[DIMENSION]={0,0,0};
	//double pin0[DIMENSION]={0,0,0};
	//for(int i=0;i<h_num;i++)
	//{
	//	int N=HYPER[i].N;
	//	fsdg<<",,"<<i<<endl;
	//	//fpn<<",,"<<i<<endl;
	//	for(int j=0;j<N;j++)
	//	{
	//		int jn=HYPER[i].NEI[j];
	//		sdgdq[A_X]=HYPER[jn].stress[A_X][0]*HYPER1[jn*h_num+i].DgDq[0]+HYPER[jn].stress[A_X][1]*HYPER1[jn*h_num+i].DgDq[1]+HYPER[jn].stress[A_X][2]*HYPER1[jn*h_num+i].DgDq[2];
	//		sdgdq[A_Y]=HYPER[jn].stress[A_Y][0]*HYPER1[jn*h_num+i].DgDq[0]+HYPER[jn].stress[A_Y][1]*HYPER1[jn*h_num+i].DgDq[1]+HYPER[jn].stress[A_Y][2]*HYPER1[jn*h_num+i].DgDq[2];
	//		sdgdq[A_Z]=HYPER[jn].stress[A_Z][0]*HYPER1[jn*h_num+i].DgDq[0]+HYPER[jn].stress[A_Z][1]*HYPER1[jn*h_num+i].DgDq[1]+HYPER[jn].stress[A_Z][2]*HYPER1[jn*h_num+i].DgDq[2];

	//		fsdg<<",,,"<<jn<<","<<sdgdq[A_X]<<","<<sdgdq[A_Y]<<","<<sdgdq[A_Z]<<endl;

	//		//pin0[A_X]=HYPER[jn].pi[A_X][0]*HYPER1[i*h_num+jn].n0ij[0]+HYPER[jn].pi[A_X][1]*HYPER1[i*h_num+jn].n0ij[1]+HYPER[jn].pi[A_X][2]*HYPER1[i*h_num+jn].n0ij[2];
	//		//pin0[A_Y]=HYPER[jn].pi[A_Y][0]*HYPER1[i*h_num+jn].n0ij[0]+HYPER[jn].pi[A_Y][1]*HYPER1[i*h_num+jn].n0ij[1]+HYPER[jn].pi[A_Y][2]*HYPER1[i*h_num+jn].n0ij[2];
	//		//pin0[A_Z]=HYPER[jn].pi[A_Z][0]*HYPER1[i*h_num+jn].n0ij[0]+HYPER[jn].pi[A_Z][1]*HYPER1[i*h_num+jn].n0ij[1]+HYPER[jn].pi[A_Z][2]*HYPER1[i*h_num+jn].n0ij[2];
	//		//fpn<<",,,"<<jn<<","<<pin0[A_X]<<","<<pin0[A_Y]<<","<<pin0[A_Z]<<endl;
	//	}
	//	sdgdq[A_X]=HYPER[i].stress[A_X][0]*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[A_X][1]*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[A_X][2]*HYPER1[i*h_num+i].DgDq[2];
	//	sdgdq[A_Y]=HYPER[i].stress[A_Y][0]*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[A_Y][1]*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[A_Y][2]*HYPER1[i*h_num+i].DgDq[2];
	//	sdgdq[A_Z]=HYPER[i].stress[A_Z][0]*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[A_Z][1]*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[A_Z][2]*HYPER1[i*h_num+i].DgDq[2];
	//	fsdg<<",,,"<<i<<","<<sdgdq[A_X]<<","<<sdgdq[A_Y]<<","<<sdgdq[A_Z]<<endl;

	//	//pin0[A_X]=HYPER[i].pi[A_X][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].pi[A_X][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].pi[A_X][2]*HYPER1[i*h_num+i].n0ij[2];
	//	//pin0[A_Y]=HYPER[i].pi[A_Y][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].pi[A_Y][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].pi[A_Y][2]*HYPER1[i*h_num+i].n0ij[2];
	//	//pin0[A_Z]=HYPER[i].pi[A_Z][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].pi[A_Z][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].pi[A_Z][2]*HYPER1[i*h_num+i].n0ij[2];
	//	//fpn<<",,,"<<i<<","<<pin0[A_X]<<","<<pin0[A_Y]<<","<<pin0[A_Z]<<endl;

	//	//fpi<<",,"<<i<<","<<HYPER[i].pi[A_X][A_X]<<","<<HYPER[i].pi[A_X][A_Y]<<","<<HYPER[i].pi[A_X][A_Z]<<endl;
	//	//fpi<<",,,"<<HYPER[i].pi[A_Y][A_X]<<","<<HYPER[i].pi[A_Y][A_Y]<<","<<HYPER[i].pi[A_Y][A_Z]<<endl;
	//	//fpi<<",,,"<<HYPER[i].pi[A_Z][A_X]<<","<<HYPER[i].pi[A_Z][A_Y]<<","<<HYPER[i].pi[A_Z][A_Z]<<endl;
	//	fst<<",,"<<i<<","<<HYPER[i].stress[A_X][A_X]<<","<<HYPER[i].stress[A_X][A_Y]<<","<<HYPER[i].stress[A_X][A_Z]<<endl;
	//	fst<<",,,"<<HYPER[i].stress[A_Y][A_X]<<","<<HYPER[i].stress[A_Y][A_Y]<<","<<HYPER[i].stress[A_Y][A_Z]<<endl;
	//	fst<<",,,"<<HYPER[i].stress[A_Z][A_X]<<","<<HYPER[i].stress[A_Z][A_Y]<<","<<HYPER[i].stress[A_Z][A_Z]<<endl;
	//}
	//fsdg.close();
	////fpi.close();
	////fpn.close();
	//fst.close();

	////
	////	ofstream r("r_CG.csv", ios::app);
	////	double r_cg_x=0;
	////	double r_cg_y=0;
	////	double r_cg_z=0;
	//
	//	for(int i=0;i<h_num;i++)
	//	{
	//		r_cg_x+=PART[i].r[A_X];
	//		r_cg_y+=PART[i].r[A_Y];
	//		r_cg_z+=PART[i].r[A_Z];
	//	}
	//
	//	r<<r_cg_x/h_num<<","<<r_cg_y/h_num<<","<<r_cg_z/h_num<<endl;
	//	r.close();

	/*	ofstream J("J.csv", ios::app);
	ofstream ti_Fi("ti_Fi.csv", ios::app);
	ofstream Fi("Fi.csv", ios::app);*/
	/*
	ofstream lam("lambda.csv", ios::app);
	lam<<HYPER[182].lambda<<","<<HYPER[911].lambda<<","<<HYPER[1640].lambda<<","<<HYPER[2369].lambda<<endl;
	lam.close();

	ofstream j("J.csv", ios::app);
	j<<HYPER[182].J<<","<<HYPER[911].J<<","<<HYPER[1640].J<<","<<HYPER[2369].J<<endl;
	j.close();


	ofstream p("P.csv", ios::app);
	p<<sqrt(HYPER[182].p[A_X]*HYPER[182].p[A_X]+HYPER[182].p[A_Z]*HYPER[182].p[A_Z])<<",";
	p<<sqrt(HYPER[911].p[A_X]*HYPER[911].p[A_X]+HYPER[911].p[A_Z]*HYPER[911].p[A_Z])<<",";
	p<<sqrt(HYPER[1640].p[A_X]*HYPER[1640].p[A_X]+HYPER[1640].p[A_Z]*HYPER[1640].p[A_Z])<<",";	
	p<<sqrt(HYPER[2369].p[A_X]*HYPER[2369].p[A_X]+HYPER[2369].p[A_Z]*HYPER[2369].p[A_Z])<<endl;
	p.close();

	/*
	if(t==1)
	{
	p<<"t"<<",";		
	d_p<<"t"<<",";
	h_p<<"t"<<",";/
	lam<<"t"<<",";
	stress<<"t"<<",";
	ti_Fi<<"t"<<",";
	Fi<<"t"<<",";
	J<<"t"<<",";//

	for(int i=0;i<h_num;i++)
	{
	p<<i<<","<<","<<",";
	lam<<i<<",";
	d_p<<i<<","<<","<<",";
	h_p<<i<<","<<","<<",";
	stress<<i<<","<<","<<",";
	ti_Fi<<i<<","<<","<<",";
	Fi<<i<<","<<","<<",";
	J<<i<<",";
	}
	p<<endl;
	lam<<endl;
	d_p<<endl;
	h_p<<endl;
	stress<<endl;
	ti_Fi<<endl;
	Fi<<endl;
	J<<endl;
	}

	p<<t<<",";
	lam<<t<<",";
	d_p<<t<<",";
	h_p<<t<<",";
	stress<<t<<",";
	ti_Fi<<t<<",";
	Fi<<t<<",";
	J<<t<<",";
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

	}
	stress<<endl<<",";
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
	}

	p<<endl;
	lam<<endl;
	d_p<<endl;
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
	J.close();
	p.close();
	lam.close();*/
}

void output_newton_data(vector<mpselastic> PART,vector<hyperelastic> HYPER,vector<hyperelastic2> HYPER1,double *fx,double *DfDx,int hyper_number,int count, int t)
{
	int h_num=hyper_number;

	stringstream ss_r;
	ss_r<<"./Newton_raphson/position"<<t<<"_"<<count<<".csv";
	stringstream ss_df;
	ss_df<<"./Newton_raphson/dfdx "<<t<<"_"<<count<<".csv";
	stringstream ss_fx;
	ss_fx<<"./Newton_raphson/fx"<<t<<"_"<<count<<".csv";
	stringstream ss_dg;
	ss_dg<<"./Newton_raphson/DgDq"<<t<<"_"<<count<<".csv";

	ofstream r(ss_r.str());
	ofstream sfx(ss_fx.str());
	ofstream Df(ss_df.str());
	ofstream fdg(ss_dg.str());

	//double sum_fx=0;
	for(int i=0;i<h_num;i++)
	{
		r<<PART[i].r[A_X]<<","<<PART[i].r[A_Y]<<","<<PART[i].r[A_Z]<<endl;

		//sum_fx+=fx[i];
		sfx<<fx[i]<<endl;
		for(int j=0;j<h_num;j++) Df<<DfDx[i*h_num+j]<<",";
		Df<<endl;
		fdg<<i<<endl;
		int Ni=HYPER[i].N;
		for(int j=0;j<Ni;j++)
		{
			int jn=HYPER[i].NEI[j];
			fdg<<","<<jn<<","<<HYPER1[i*h_num+jn].DgDq[A_X]<<","<<HYPER1[i*h_num+jn].DgDq[A_Y]<<","<<HYPER1[i*h_num+jn].DgDq[A_Z]<<endl;
		}
	}
	//sfx<<count<<","<<sum_fx/h_num<<endl;

	r.close();
	sfx.close();
	Df.close();
	fdg.close();
}

void output_newton_data2(double E, double *XX, double *fx,int hyper_number, int count, int t)
{
	int h_num=hyper_number;
	stringstream ss_E;
	ss_E<<"./Newton_raphson/E"<<t<<".csv";

	stringstream ss_lam;
	ss_lam<<"./Newton_raphson/lambda"<<t<<"_"<<count<<".csv";

	stringstream ss_d;
	ss_d<<"./Newton_raphson/d"<<t<<"_"<<count<<".csv";

	if(count==1)
	{
		ofstream init0(ss_E.str(), ios::trunc);

		init0.close();
	}

	ofstream e(ss_E.str(), ios::app);
	ofstream lam(ss_lam.str());
	ofstream fd(ss_d.str());

	e<<count<<","<<E<<endl;

	double sum_lam=0;
	for(int i=0;i<h_num;i++)
	{
		lam<<XX[i]<<endl;
		fd<<fx[i]<<endl;
	}

	e.close();
	lam.close();
	fd.close();
}

void output_energy(mpsconfig CON, vector<mpselastic> PART, vector<hyperelastic> HYPER,int t)
{
	//	cout<<"弾性ポテンシャル計算";
	int h_num=HYPER.size();
	int p_num=PART.size();
	//double d_Fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	//double c10=0.;//CON.get_c10();
	//c10=4484.5883;	//加嶋さんの実験値
	////double c01=CON.get_c01();
	////vector<double>	W;
	////W.reserve(h_num);
	////for(int i=0;i<h_num;i++)	W.emplace_back(0);
	//double dC[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	//double dC2[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double G=9.8;
	/*
	for(int j=0;j<4;j++)
	{
	if(j==0)	i=182;
	else if(j==1) i=911;
	else if(j==2) i=1640;
	else if(j==3) i=2369;

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
	}*/

	//for(int i=0;i<h_num;i++)
	//{
	//	double J=HYPER[i].J;
	//	if(J<0){
	//		d_Fi[0][0]=-1/pow(-J,1/3)*HYPER[i].Fi[0][0];	d_Fi[0][1]=-1/pow(-J,1/3)*HYPER[i].Fi[0][1];	d_Fi[0][2]=-1/pow(-J,1/3)*HYPER[i].Fi[0][2];
	//		d_Fi[1][0]=-1/pow(-J,1/3)*HYPER[i].Fi[1][0];	d_Fi[1][1]=-1/pow(-J,1/3)*HYPER[i].Fi[1][1];	d_Fi[1][2]=-1/pow(-J,1/3)*HYPER[i].Fi[1][2];
	//		d_Fi[2][0]=-1/pow(-J,1/3)*HYPER[i].Fi[2][0];	d_Fi[2][1]=-1/pow(-J,1/3)*HYPER[i].Fi[2][1];	d_Fi[2][2]=-1/pow(-J,1/3)*HYPER[i].Fi[2][2];
	//	}
	//	else
	//	{
	//		d_Fi[0][0]=1/pow(J,1/3)*HYPER[i].Fi[0][0];	d_Fi[0][1]=1/pow(J,1/3)*HYPER[i].Fi[0][1];	d_Fi[0][2]=1/pow(J,1/3)*HYPER[i].Fi[0][2];
	//		d_Fi[1][0]=1/pow(J,1/3)*HYPER[i].Fi[1][0];	d_Fi[1][1]=1/pow(J,1/3)*HYPER[i].Fi[1][1];	d_Fi[1][2]=1/pow(J,1/3)*HYPER[i].Fi[1][2];
	//		d_Fi[2][0]=1/pow(J,1/3)*HYPER[i].Fi[2][0];	d_Fi[2][1]=1/pow(J,1/3)*HYPER[i].Fi[2][1];	d_Fi[2][2]=1/pow(J,1/3)*HYPER[i].Fi[2][2];
	//	}

	//	dC[0][0]=d_Fi[0][0]*d_Fi[0][0]+d_Fi[1][0]*d_Fi[1][0]+d_Fi[2][0]*d_Fi[2][0];
	//	dC[0][1]=d_Fi[0][0]*d_Fi[0][1]+d_Fi[1][0]*d_Fi[1][1]+d_Fi[2][0]*d_Fi[2][1];
	//	dC[0][2]=d_Fi[0][0]*d_Fi[0][2]+d_Fi[1][0]*d_Fi[1][2]+d_Fi[2][0]*d_Fi[2][2];
	//	dC[1][0]=d_Fi[0][1]*d_Fi[0][0]+d_Fi[1][1]*d_Fi[1][0]+d_Fi[2][1]*d_Fi[2][0];
	//	dC[1][1]=d_Fi[0][1]*d_Fi[0][1]+d_Fi[1][1]*d_Fi[1][1]+d_Fi[2][1]*d_Fi[2][1];
	//	dC[1][2]=d_Fi[0][1]*d_Fi[0][2]+d_Fi[1][1]*d_Fi[1][2]+d_Fi[2][1]*d_Fi[2][2];
	//	dC[2][0]=d_Fi[0][2]*d_Fi[0][0]+d_Fi[1][2]*d_Fi[1][0]+d_Fi[2][2]*d_Fi[2][0];
	//	dC[2][1]=d_Fi[0][2]*d_Fi[0][1]+d_Fi[1][2]*d_Fi[1][1]+d_Fi[2][2]*d_Fi[2][1];
	//	dC[2][2]=d_Fi[0][2]*d_Fi[0][2]+d_Fi[1][2]*d_Fi[1][2]+d_Fi[2][2]*d_Fi[2][2];

	//	dC2[0][0]=dC[A_X][0]*dC[0][A_X]+dC[A_X][1]*dC[1][A_X]+dC[A_X][2]*dC[2][A_X];
	//	dC2[0][1]=dC[A_X][0]*dC[0][A_Y]+dC[A_X][1]*dC[1][A_Y]+dC[A_X][2]*dC[2][A_Y];
	//	dC2[0][2]=dC[A_X][0]*dC[0][A_Z]+dC[A_X][1]*dC[1][A_Z]+dC[A_X][2]*dC[2][A_Z];
	//	dC2[1][0]=dC[A_Y][0]*dC[0][A_X]+dC[A_Y][1]*dC[1][A_X]+dC[A_Y][2]*dC[2][A_X];
	//	dC2[1][1]=dC[A_Y][0]*dC[0][A_Y]+dC[A_Y][1]*dC[1][A_Y]+dC[A_Y][2]*dC[2][A_Y];
	//	dC2[1][2]=dC[A_Y][0]*dC[0][A_Z]+dC[A_Y][1]*dC[1][A_Z]+dC[A_Y][2]*dC[2][A_Z];
	//	dC2[2][0]=dC[A_Z][0]*dC[0][A_X]+dC[A_Z][1]*dC[1][A_X]+dC[A_Z][2]*dC[2][A_X];
	//	dC2[2][1]=dC[A_Z][0]*dC[0][A_Y]+dC[A_Z][1]*dC[1][A_Y]+dC[A_Z][2]*dC[2][A_Y];
	//	dC2[2][2]=dC[A_Z][0]*dC[0][A_Z]+dC[A_Z][1]*dC[1][A_Z]+dC[A_Z][2]*dC[2][A_Z];

	//	double trace_dC=dC[0][0]+dC[1][1]+dC[2][2];
	//	double trace_dC2=dC2[0][0]+dC2[1][1]+dC2[2][2];

	//	double Ic=trace_dC;
	//	double IIc=0.50*(trace_dC*trace_dC-trace_dC2);
	//	W[i]=c10*(Ic-3)+c01*(IIc-3);
	//}//*/
	//	cout<<"----------OK"<<endl;

	/*	ofstream e("E.csv", ios::app);
	ofstream e_T("E_T.csv", ios::app);
	ofstream e_g("E_g.csv", ios::app);
	ofstream e_W("E_W.csv", ios::app);
	ofstream e_lam("E_lam.csv", ios::app);//*/
	ofstream e_sum("E_sum.csv", ios::app);

	stringstream ss_e;
	ss_e<<"./Energy/E"<<t<<".csv";
	ofstream fe(ss_e.str());

	double V=get_volume(&CON);
	double mh=V*CON.get_hyper_density();
	double ms=V*CON.get_silicone_density();
	double mi=0.;
	/*	double sum_e_T=0;
	double sum_e_g=0;
	double sum_e_lam=0;
	double sum_e_W=0;
	double sum_e=0;
	double vv=0;
	double energy=0;

	for(int j=0;j<4;j++)
	{
	if(j==0)	i=182;
	else if(j==1) i=911;
	else if(j==2) i=1640;
	else if(j==3) i=2369;

	vv=HYPER[i].p[0]*HYPER[i].p[0]+HYPER[i].p[1]*HYPER[i].p[1]+HYPER[i].p[2]*HYPER[i].p[2];
	energy=0.5/mi*vv+W[i]*V+HYPER[i].lambda*(1-HYPER[i].J)*V;
	//energy=0.5/mi*vv+mi*9.8*PART[i].r[A_Z]+W[i]*V+HYPER[i].lambda*(1-HYPER[i].J)*V;
	e<<energy<<",";
	e_T<<0.5/mi*vv<<",";
	e_g<<mi*9.8*PART[i].r[A_Z]<<",";
	e_W<<W[i]*V<<",";
	e_lam<<HYPER[i].lambda*(1-HYPER[i].J)*V<<",";
	}*/

	/*	if(t==1)
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
	e_lam<<t<<",";*/

	double sum_e_T=0;
	double sum_e_g=0;
	double sum_e_lam=0;
	double sum_e_W=0;
	double sum_e=0;
	double vv=0;
	double energy=0;

	for(int i=0;i<h_num;i++)
	{
		if(PART[i].toFEM==ON)
		{
			mi=mh;
		}
		else
		{
			mi=ms;
		}
		vv=HYPER[i].p[0]*HYPER[i].p[0]+HYPER[i].p[1]*HYPER[i].p[1]+HYPER[i].p[2]*HYPER[i].p[2];
		//energy=0.5/mi*vv+W[i]*V+HYPER[i].lambda*(1-HYPER[i].J)*V;
		if(CON.get_flag_G()==ON)	energy=0.5/mi*vv+mi*G*PART[i].r[A_Z]+HYPER[i].W*V+HYPER[i].lambda*(1-HYPER[i].J)*V;
		if(CON.get_flag_G()==OFF)	energy=0.5/mi*vv+HYPER[i].W*V+HYPER[i].lambda*(1-HYPER[i].J)*V;
		sum_e_T+=0.5/mi*vv;
		sum_e_g+=mi*G*PART[i].r[A_Z];
		sum_e_lam+=HYPER[i].lambda*(1-HYPER[i].J)*V;
		sum_e_W+=HYPER[i].W*V;
		sum_e+=energy;
		fe<<energy<<","<<0.5/mi*vv<<","<<mi*G*PART[i].r[A_Z]<<","<<HYPER[i].W*V<<","<<HYPER[i].lambda*(1-HYPER[i].J)*V<<endl;
		/*	e<<energy<<",";
		e_T<<0.5/mi*vv<<",";
		e_g<<mi*G*PART[i].r[A_Z]<<",";
		e_W<<W[i]*V<<",";
		e_lam<<HYPER[i].lambda*(1-HYPER[i].J)*V<<",";//*/
	}
	fe.close();
	/*	e<<endl;
	e_T<<endl;
	e_g<<endl;
	e_W<<endl;
	e_lam<<endl;*/


	if(t==1)	e_sum<<"E,E_T,E_g,E_W,E_lam\n";
	e_sum<<sum_e<<","<<sum_e_T<<","<<sum_e_g<<","<<sum_e_W<<","<<sum_e_lam<<endl;//*/

	/*	e.close();
	e_T.close();
	e_g.close();
	e_W.close();
	e_lam.close();//*/
	e_sum.close();
}

//detを求める関数　※自己流のため自信な
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

//逆行列を求める関数 function.hの圧力計算用関数に類似した名前の関数があったのですみわけ
void calc_inverse_matrix_for_NR(int N, double *a)
{
	//N=未知数
	double buf=0;

	double *inv_a=new  double[N*N];					//逆行列格納

	//単位行列を作る
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<N;j++)
		{
			if(i==j) inv_a[j*N+i]=1;
			else		inv_a[j*N+i]=0;
		}
	}

	//掃き出し法
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
			cout<<"特異行列である"<<endl;
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

void BiCGStab2_method(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X)
{
	ofstream t("time_log.dat",ios::app);
	t<<clock()*CLOCKS_PER_SEC<<"	";

	//val :非ゼロ要素の値
	//ind:非ゼロ要素の列番号格納配列
	//ptr:各行の要素がvalの何番目からはじまるのかを格納
	//X[n]:解

	int count=0;
	double pk=1;
	double E=1;//誤差
	double alp,beta,rr,w,ita;

	double *r=new double [pn];	//残差
	double *P=new double [pn];	//探索ベクトル
	double *Pj=new double [pn];
	double *rj=new double [pn];
	double *AP=new double [pn];
	double *AtPj=new double [pn];
	double *e=new double [pn];
	double *Ae=new double [pn];
	double *y=new double [pn];
	double *u=new double [pn];
	double *W=new double [pn];
	double *Z=new double [pn];
	double *e_pre=new double [pn];
	beta=0;
	w=0;
	ita=0;

	for(int n=0;n<pn;n++) //初期値
	{
		X[n]=0;
		r[n]=B[n];
	}

	for(int n=0;n<pn;n++)
	{
		P[n]=r[n];//初期化
		rj[n]=r[n];
		Pj[n]=rj[n];
		AP[n]=0;
		y[n]=0;
		u[n]=0;
		W[n]=0;
		Z[n]=0;
		e[n]=0;
		e_pre[n]=0;//1ステップ前のe[]
	}
	double rr0=0;
	for(int n=0;n<pn;n++) rr0+=r[n]*r[n];
	cout<<"rr0="<<rr0<<endl;
	cout<<"BiCGstab2法スタート  -----未知数="<<pn<<"  ---";
	unsigned int time=GetTickCount();
	double ep=CON->get_FEMCGep();//収束判定

	while(E>ep)// EP=CON->get_CGep();//収束判定(convergence test)
	{
		count++;

		for(int n=0;n<pn;n++)
		{
			//P[n]=r[n]+beta*(P[n]-w*AP[n]);
			P[n]=r[n]+beta*(P[n]-u[n]);
		} 

		////pk(rとrjの内積)を求める
		pk=0;
		for(int n=0;n<pn;n++) pk+=r[n]*rj[n];

		//////////////alpを求める
		for(int n=0;n<pn;n++)
		{      
			AP[n]=0;
			//cout<<"ptrn"<<ptr[n]<<"ptrn+1"<<ptr[n+1]<<endl;
			for(int m=ptr[n];m<ptr[n+1];m++) AP[n]+=val[m]*P[ind[m]];
			//for(int m=0;m<pn;m++) AP[n]+=A[n][m]*P[m];
		}
		double APrj=0;
		for(int n=0;n<pn;n++)  APrj+=rj[n]*AP[n];
		alp=pk/APrj;
		//cout<<"alp="<<alp<<" APrj="<<APrj<<endl;
		//////////////////////

		for(int n=0;n<pn;n++) y[n]=e[n]-r[n]-alp*W[n]+alp*AP[n];

		for(int n=0;n<pn;n++)
		{
			e_pre[n]=e[n];//変更前の値を記憶
			e[n]=r[n]-alp*AP[n];
		}

		for(int n=0;n<pn;n++)
		{
			Ae[n]=0;
			//for(int m=0;m<pn;m++) Ae[n]+=A[n][m]*e[m];
			for(int m=ptr[n];m<ptr[n+1];m++) Ae[n]+=val[m]*e[ind[m]];
		}

		if(count%2!=0)
		{
			double e_dot_ae  = 0;
			double ae_dot_ae = 0;
			for(int n=0;n<pn;n++)
			{
				e_dot_ae+=e[n]*Ae[n];
				ae_dot_ae+=Ae[n]*Ae[n];
			}
			w = e_dot_ae / ae_dot_ae;
			ita=0;
		}
		else
		{
			double e_dot_ae  = 0;
			double y_dot_y=0;
			double y_dot_e=0;
			double Ae_dot_y=0;
			double Ae_dot_Ae=0;
			for(int n=0;n<pn;n++) e_dot_ae+=e[n]*Ae[n];
			for(int n=0;n<pn;n++) y_dot_y+=y[n]*y[n];
			for(int n=0;n<pn;n++) y_dot_e+=y[n]*e[n];
			for(int n=0;n<pn;n++) Ae_dot_y+=Ae[n]*y[n];
			for(int n=0;n<pn;n++) Ae_dot_Ae+=Ae[n]*Ae[n];

			double co=Ae_dot_Ae*y_dot_y-Ae_dot_y*Ae_dot_y;

			w=y_dot_y*e_dot_ae-y_dot_e*Ae_dot_y;
			w/=co;

			ita=Ae_dot_Ae*y_dot_e-Ae_dot_y*e_dot_ae;
			ita/=co;
		}


		for(int n=0;n<pn;n++) 
		{
			u[n]=w*AP[n]+ita*(e_pre[n]-r[n]+beta*u[n]);
			Z[n]=w*r[n]+ita*Z[n]-alp*u[n];
			X[n]+=alp*P[n]+Z[n];
			r[n]=e[n]-ita*y[n]-w*Ae[n];
		}

		///beta
		beta=0;
		for(int n=0;n<pn;n++) beta+=r[n]*rj[n];
		beta/=w*APrj;

		//cout<<"beta="<<beta<<endl;
		//////////////////////

		//W[n]
		for(int n=0;n<pn;n++) W[n]=Ae[n]+beta*AP[n];

		//////////////////誤差
		rr=0;
		for(int n=0;n<pn;n++) rr+=r[n]*r[n];
		//E=rr/rr0;
		E=sqrt(rr);
		if(count==1||count%10==0)	cout<<"E="<<E<<" count="<<count<<endl;
		////////////////////////
	}

	cout<<"反復回数="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;

	delete [] r;
	delete [] Pj;
	delete [] rj;
	delete [] AP;
	delete [] AtPj;
	delete [] P;
	delete [] e;
	delete [] Ae;
	delete [] y;
	delete [] u;
	delete [] W;
	delete [] Z;
	delete [] e_pre;

	t<<clock()*CLOCKS_PER_SEC<<"	";
	t.close();
}

void iccg2(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X)
{
	//val :ゼロ要素の値
	//ind:非ゼロ要素の列番号格納配列
	//ptr:各行の要素がvalの何番目からはじまるのかを格納
	//X[n]:解

	double accel=0.87;//CON->get_CGaccl();//加速ファクタ

	int num2=0;//対角成分を含む、下三角行列だけを考慮にいれた非ゼロ要素数
	for(int k=0;k<pn;k++) for(int m=ptr[k];m<ptr[k+1];m++) if(ind[m]<=k) num2++;	
	if(num2!=(number-pn)/2+pn) cout<<"ERROR"<<endl;

	double *val2=new double [num2];
	int *ind2 = new int [num2];
	int *ptr2 = new int [pn+1];
	double *r=new double [pn];
	double *P=new double [pn];

	num2=0;
	for(int k=0;k<pn;k++)
	{	
		ptr2[k]=num2;
		r[k]=0;
		P[k]=0;
		for(int m=ptr[k];m<ptr[k+1];m++)///k行目の非０要素
		{
			if(ind[m]<=k)
			{
				val2[num2]=val[m];
				ind2[num2]=ind[m];
				if(ind[m]==k) val2[num2]*=accel;//加速ﾌｧｸﾀ
				num2++;
			}
		}
	}
	ptr2[pn]=num2;//これをしておかないと、最後に(int m=ptr2[k];m<ptr2[k+1];m++)みたいなことができない

	int *NUM = new int [pn];
	for(int k=0;k<pn;k++) NUM[k]=0;

	for(int k=0;k<pn;k++)
	{	
		for(int m=ptr2[k];m<ptr2[k+1];m++)///k行目の非０要素
		{
			int J=ind2[m];
			NUM[J]=NUM[J]+1;
		}
	}
	double **VAL=new double *[pn];//ゼロ要素の値
	for(int i=0;i<pn;i++) VAL[i]=new double [NUM[i]];
	int **IND = new int *[pn];//非ゼロ要素の行番号格納配列
	for(int i=0;i<pn;i++) IND[i]=new int [NUM[i]];

	/////////////////////////////////////iccg法
	double alp,beta;
	double rLDLt_r;
	double E=1;//誤差
	double *AP = new double [pn];
	double *y=new double [pn];
	double *LDLt_r= new double [pn];
	double *D1 = new double [pn];//D行列

	/////不完全コレスキｰ分解
	for(int k=0;k<pn;k++)
	{	
		for(int m=ptr2[k];m<ptr2[k+1];m++)///k行目の非０要素
		{
			int i=ind2[m];//列番号
			if(i==0)
			{
				val2[m]=val2[m];
				if(val2[m]<0.0001 &&val2[m]>=0) val2[m]=0.0001;
				else if(val2[m]>-0.0001 &&val2[m]<=0) val2[m]=-0.0001;

			}
			if(i>0 && i<k)
			{
				double sum=0;

				for(int j=ptr2[k];j<m;j++)
				{	
					for(int J=ptr2[i];J<ptr2[i+1];J++)//i行目のなかから列の一致するものを探している。少し手間か？
					{
						if(ind2[J]==ind2[j]) sum+=val2[j]*val2[J]*D1[ind2[j]];
					}
				}
				val2[m]=val2[m]-sum;
			}
			if(i==k)
			{
				double sum=0;
				for(int j=ptr2[k];j<m;j++) sum+=val2[j]*val2[j]*D1[ind2[j]];
				val2[m]=val2[m]-sum;
				if(val2[m]<0.0001 &&val2[m]>=0) val2[m]=0.0001;
				else if(val2[m]>-0.0001 &&val2[m]<=0) val2[m]=-0.0001;
				D1[k]=1/val2[m];
				//if(val2[m]>0) cout<<"EE"<<endl;
			}
		}
	}    
	///不完全コレスキー分解完了/////////*/

	///列を基準にした配列に値を代入
	for(int k=0;k<pn;k++) NUM[k]=0;
	for(int k=0;k<pn;k++)
	{	
		for(int m=ptr2[k];m<ptr2[k+1];m++)///k行目の非０要素
		{
			int J=ind2[m];
			VAL[J][NUM[J]]=val2[m];
			IND[J][NUM[J]]=k;
			NUM[J]=NUM[J]+1;
		}
	}////////*/

	/////////////////y[i]をもとめ、LDLt_r[i]をもとめる。
	for(int i=0;i<pn;i++)
	{
		if(i==0) y[0]=r[0]/val2[0]; //式（3.77） 
		else
		{
			double sum=0;
			/////////        
			for(int m=ptr2[i];m<ptr2[i+1]-1;m++) sum+=val2[m]*y[ind2[m]];//式（3.78）
			int m=ptr2[i+1]-1;
			y[i]=(r[i]-sum)/val2[m];
		}
	}////y[i]がもとまった。
	for(int i=pn-1;i>=0;i--)
	{
		double sum=0;
		for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
		LDLt_r[i]=y[i]-D1[i]*sum;	
	}
	/////////////////*/

	for(int n=0;n<pn;n++){ P[n]=LDLt_r[n];
	cout<<"P"<<n<<"	"<<P[n]<<endl;
	}

	cout<<"ICCG法:未知数="<<pn<<" ---";
	unsigned int time=GetTickCount();
	int count=0;
	double ep=CON->get_FEMCGep();//収束判定
	rLDLt_r=0;
	for(int n=0;n<pn;n++) rLDLt_r+=r[n]*LDLt_r[n];//最初のrLDLt_rだけここで求める
	while(E>ep)
	{
		count++;
		if(count==pn) cout<<"count=pn"<<endl;
		//////////////alpを求める
		double PAP=0;
#pragma omp parallel for reduction(+:PAP)
		for(int n=0;n<pn;n++)
		{
			AP[n]=0;
			for(int m=ptr[n];m<ptr[n+1];m++) AP[n]+=val[m]*P[ind[m]];
			PAP+=P[n]*AP[n];
		}

		//for(int n=0;n<pn;n++)  PAP+=P[n]*AP[n];
		alp=rLDLt_r/PAP;
		cout<<"alp="<<alp<<endl;
		//////////////////////
		E=0;
#pragma omp parallel for reduction(+:E)
		for(int n=0;n<pn;n++)
		{
			X[n]+=alp*P[n];// X(k+1)=X(k)+alp*P 更新後の場所
			r[n]-=alp*AP[n];// r=r-alp*AP       更新後の残差
			E+=r[n]*r[n];						//更新後の誤差
		}
		E=sqrt(E);
		cout<<"E="<<E<<endl;
		////////////////////////

		///////////////////////beta
		beta=1.0/rLDLt_r;
		rLDLt_r=0;

		/////////////////y[i]をもとめ、LDLt_r[i]をもとめる。
		for(int i=0;i<pn;i++)
		{
			if(i==0) y[0]=r[0]/val2[0]; //式（3.77） 新
			else
			{
				double sum=0;
				for(int m=ptr2[i];m<ptr2[i+1]-1;m++)//対角成分は除くからptr[i+1]-1
				{
					sum+=val2[m]*y[ind2[m]];//式（3.78）
				}
				int m=ptr2[i+1]-1;
				y[i]=(r[i]-sum)/val2[m];
			}
		}////y[i]がもとまった。

		/////////LDLt_r[i]を求める
		for(int i=pn-1;i>=0;i--)
		{
			double sum=0;
			for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];

			LDLt_r[i]=y[i]-D1[i]*sum;	
		}
		/////////////////*/

		for(int n=0;n<pn;n++) rLDLt_r+=r[n]*LDLt_r[n];

		beta=beta*rLDLt_r;
		/////////////////*/

		///////////////////// P=r+beta*P
		for(int n=0;n<pn;n++) P[n]=LDLt_r[n]+beta*P[n];//iccg
		if(count==1||count%10==0)	cout<<"E="<<E<<" count="<<count<<endl;
	}
	cout<<"反復回数="<<count<<" time="<<(GetTickCount()-time)*0.001<<"/";

	delete [] AP;

	delete [] y;
	delete [] LDLt_r;
	delete [] D1;

	delete [] val2;
	delete [] ind2;
	delete [] ptr2;
	delete [] r;
	delete [] P;
	for(int i=0;i<pn;i++)
	{
		delete [] VAL[i];
		delete [] IND[i];
	}
	delete [] VAL;
	delete [] IND;
	delete [] NUM;
	//*count2=count;//反復回数を格納して返す
}

void CG3D(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X)
{
	cout<<"CG法:未知数="<<pn<<" ---";
	unsigned int time=GetTickCount();
	int count=0;
	double ep=CON->get_FEMCGep();//収束判定

	double alp,beta;
	double E=1;//誤差
	double *AP = new double [pn];
	double *r=new double [pn];
	double *P=new double [pn];

	while(E>ep)
	{
		count++;
		if(count==pn) cout<<"count=pn"<<endl;
		//////////////alpを求める
		double PAP=0;
		double Pr=0;
#pragma omp parallel for reduction(+:PAP)
		for(int n=0;n<pn;n++)
		{
			AP[n]=0;
			for(int m=ptr[n];m<ptr[n+1];m++) AP[n]+=val[m]*P[ind[m]];
			PAP+=P[n]*AP[n];
			Pr+=P[n]*r[n];
		}

		//for(int n=0;n<pn;n++)  PAP+=P[n]*AP[n];
		alp=Pr/PAP;
		cout<<"alp="<<alp<<endl;
		//////////////////////
		E=0;
#pragma omp parallel for reduction(+:E)
		for(int n=0;n<pn;n++)
		{
			X[n]+=alp*P[n];// X(k+1)=X(k)+alp*P 更新後の場所
			r[n]-=alp*AP[n];// r=r-alp*AP       更新後の残差
			E+=r[n]*r[n];						//更新後の誤差
		}
		E=sqrt(E);
		if(count==1||count%10==0)	cout<<"E="<<E<<" count="<<count<<endl;
		////////////////////////

		///////////////////////beta
		double rAP=0;
		for(int n=0;n<pn;n++)	rAP+=r[n]*AP[n];
		beta=-rAP/PAP;
		cout<<"beta="<<beta<<endl;
		for(int n=0;n<pn;n++)	P[n]=r[n]+beta*P[n];

	}
	cout<<"反復回数="<<count<<" time="<<(GetTickCount()-time)*0.001<<"/";

	delete [] AP;
	delete [] r;
	delete [] P;

}

void GaussSeidelvh(double *A, int pn, double *b,double ep)
{

	double tmp;
	double e = 0.0;
	double *x=new double [pn];
	for(int k = 0; k < 10000; k++){
		e = 0.0;
		for(int i = 0; i < pn; i++){
			tmp = x[i];
			x[i] = b[i];
			for(int j = 0; j <pn; j++){
				x[i] -= (j != i ? A[i*pn+j] * x[j] : 0.0);
			}
			x[i] /= A[i*pn+i];
			e+=fabs(tmp - x[i]);
		}
		if(e <= ep)break;
	}
	double temp[3];
	int part = 0;
	for(int i = 0; i < pn; i++){
		b[i] = x[i];
	}
}
