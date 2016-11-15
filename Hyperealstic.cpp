#include "stdafx.h"	

void calc_half_p(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,bool repetation,double **F);
void renew_lambda(mpsconfig &CON,vector<mpselastic>PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int t);
void calc_differential_p(mpsconfig &CON,vector<mpselastic>PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,double **F);
void calc_transposed_inverse_matrix(double **M,bool transport,bool inversion);
double calc_det(double **M,int N);
double calc_det3(double **M);
void calc_stress(mpsconfig &CON,vector<hyperelastic> &HYPER);
void calc_constant(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1);
void calc_inverse_matrix_for_NR(int N, double *a);
void newton_raphson(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int t,double **F);
void newton_raphson2(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int t,double **F);
void calc_F(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1);
void calc_newton_function(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,double *lambda,double *fx,double *DfDx,int hyper_number,int count,int t,double **F);
void calc_newton_function2(mpsconfig &CON,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,double *old_hpz,double *lambda, double *fx,double *DfDx,int hyper_number,int count,int t,double **F);
void calc_newton_function3(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,double *lambda,double *fx,double *DfDx,int hyper_number,int count,int t,double **F);
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
void output_newton_data3(double *w_fx, double *w_DfDx, double *n_rx, double *n_ry, double *n_rz, double *part_DfDx, int hyper_number, int count, int t);
void output_newton2_data1(vector<hyperelastic> HYPER, double *fx, double *DfDx, int hyper_number, int count, int t);
void calc_gravity(mpsconfig CON,vector<hyperelastic> &HYPER,int hyper_number);
void calculation_vec_norm(vector<mpselastic> PART, vector<hyperelastic> &HYPER, int hyper_number,int particle_number,int t);
void output_energy(mpsconfig CON, vector<mpselastic> PART, vector<hyperelastic> HYPER,int t);
void contact_judge(mpsconfig &CON, vector<mpselastic> PART,vector<hyperelastic> &HYPER,double max_h,int t);
void BiCGStab2_method(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X);
void iccg2(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X);
void CG3D(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X);
void GaussSeidelvh(double *A, int pn, double *b,double ep);
void QP(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int t,double **F);
void calc_wg(mpsconfig &CON, vector<mpselastic> PART, vector<hyperelastic> &HYPER, vector<hyperelastic2> &HYPER1, double **dq, double *w, double *g, double *J, double **F, double **ti_F, double **S);

void calc_hyper(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t,double **F)
{	

	ofstream time("time_log.dat",ios::app);
	clock_t	start_t=clock();
	time<<start_t*CLOCKS_PER_SEC<<"	";

	int h_num=HYPER.size();
	double Dt=CON.get_dt();
	double mi=CON.get_hyper_density()*get_volume(&CON);

	int switch_w=OFF;

	cout<<"h_num="<<h_num<<endl;
	cout<<"Hypercalculation starts."<<endl;


	if(t==1)
	{
		ofstream init0("lambda_n1.csv",ios::trunc);
		ofstream init1("lambda_p1.csv",ios::trunc);
		ofstream init3("lambda_n2.csv",ios::trunc);
		init0.close();
		init1.close();
		init3.close();

		for(int i=0;i<h_num;i++)	for(int D=0;D<DIMENSION;D++)	PART[i].q0[D]=0;
		for(int i=0;i<h_num;i++)	for(int D=0;D<DIMENSION;D++)	PART[i].q0[D]=PART[i].r[D];
		calc_constant(CON,PART,HYPER,HYPER1);
		calc_stress(CON,HYPER);
	}
	
	
	if(CON.get_flag_vis()==ON)	calc_vis_f(CON,PART,HYPER,t);

	if(t==1 || t%CON.get_interval()==0)
	{
		output_hyper_data(PART,HYPER,HYPER1,t);
		momentum_movie_AVS(CON,t,PART,HYPER,F);
		output_energy(CON,PART,HYPER,t);
	}

	newton_raphson(CON,PART,HYPER,HYPER1,t,F);
	//newton_raphson2(CON,PART,HYPER,HYPER1,t,F);

	int f_wall=CON.get_flag_wall();
	//�ǌv�Z��
	if(f_wall==OFF)
	{
		calc_half_p(CON,PART,HYPER,HYPER1,0,F);
		calc_F(CON,PART,HYPER,HYPER1);
		calc_stress(CON,HYPER);
		calc_differential_p(CON,PART,HYPER,HYPER1,F);
		renew_lambda(CON,PART,HYPER,HYPER1,t);
		calc_half_p(CON,PART,HYPER,HYPER1,1,F);
	}
	//�ǌv�Z�L
	else
	{
		double *old_r_z=new double [h_num];
		int *Nw_n=new int [h_num];
		for(int i=0;i<h_num;i++)
		{
			Nw_n[i]=0;
			old_r_z[i]=PART[i].r[A_Z];
		}

		calc_half_p(CON,PART,HYPER,HYPER1,0,F);

		double *old_hpz=new double [h_num];

		int Nw=0;
		for(int i=0;i<h_num;i++)
		{
			old_hpz[i]=HYPER[i].half_p[A_Z];
			if(PART[i].r[A_Z]<0)
			{
				PART[i].r[A_Z]=0;
				HYPER[i].half_p[A_Z]=-1*old_r_z[i]/Dt*mi;
				Nw_n[Nw]=i;
				Nw++;
			}
		}
		

		calc_F(CON,PART,HYPER,HYPER1);
		calc_stress(CON,HYPER);	
		calc_differential_p(CON,PART,HYPER,HYPER1,F);

		if(Nw>0)
		{
			ofstream fs0("lambda_n1.csv",ios::app);

			stringstream ss;
			ss<<"./Position/position_before_r_changed"<<t<<".csv";
			ofstream fs(ss.str());

			stringstream ss3;
			ss3<<"./Half_P/half_p_before_r_changed"<<t<<".csv";
			ofstream fs3(ss3.str());

			stringstream ss5;
			ss5<<"./P/P_before_p_changed"<<t<<".csv";
			ofstream fs5(ss5.str());


			for(int i=0;i<h_num;i++)	
			{
				fs0<<HYPER[i].lambda<<",";
				fs<<i<<","<<PART[i].r[A_X]<<","<<PART[i].r[A_Y]<<","<<old_r_z[i]<<","<<PART[i].r[A_X]<<","<<PART[i].r[A_Y]<<","<<PART[i].r[A_Z]<<endl;
				fs3<<i<<","<<HYPER[i].half_p[A_X]<<","<<HYPER[i].half_p[A_Y]<<","<<old_hpz[i]<<","<<HYPER[i].half_p[A_X]<<","<<HYPER[i].half_p[A_Y]<<","<<HYPER[i].half_p[A_Z]<<endl;
				fs5<<i<<","<<HYPER[i].p[A_X]<<","<<HYPER[i].p[A_Y]<<","<<HYPER[i].p[A_Z]<<endl;
			}
			fs0<<endl;

			fs0.close();
			fs.close();
			fs3.close();
			fs5.close();
		}

		renew_lambda(CON,PART,HYPER,HYPER1,t);
		calc_half_p(CON,PART,HYPER,HYPER1,1,F);

		if(Nw>0)
		{
			ofstream fs1("lambda_p1.csv",ios::app);
			for(int i=0;i<h_num;i++)	fs1<<HYPER[i].lambda<<",";
			fs1<<endl;
			fs1.close();
		}


		delete[]	old_r_z;

		for(int i=0;i<Nw;i++)
		{
			int j=Nw_n[i];
			double p_norm=sqrt(HYPER[j].p[A_X]*HYPER[j].p[A_X]+HYPER[j].p[A_Y]*HYPER[j].p[A_Y]+HYPER[j].p[A_Z]*HYPER[j].p[A_Z]);
			double p_vector[DIMENSION]={HYPER[j].p[A_X]/p_norm,HYPER[j].p[A_Y]/p_norm,HYPER[j].p[A_Z]/p_norm};
			double E=0.5/mi*(HYPER[j].p[A_X]*HYPER[j].p[A_X]+HYPER[j].p[A_Y]*HYPER[j].p[A_Y]+HYPER[j].p[A_Z]*HYPER[j].p[A_Z])+(0.5/mi*old_hpz[j]*old_hpz[j]-0.5/mi*HYPER[j].half_p[A_Z]*HYPER[j].half_p[A_Z]);
			HYPER[j].p[A_X]=mi*p_vector[A_X]*sqrt(2/mi*E);
			HYPER[j].p[A_Y]=mi*p_vector[A_Y]*sqrt(2/mi*E);
			HYPER[j].p[A_Z]=mi*p_vector[A_Z]*sqrt(2/mi*E);	//*/ //c1_v1_wall_dl0.5_2
		}//*/
		delete[]	Nw_n;


			//�J��Ԃ��v�Z	
		if(Nw>0)
		{
			stringstream ss6;
			ss6<<"./P/P_after_p_changed"<<t<<".csv";
			ofstream fs6(ss6.str());
			for(int i=0;i<h_num;i++)	fs6<<i<<","<<HYPER[i].p[A_X]<<","<<HYPER[i].p[A_Y]<<","<<HYPER[i].p[A_Z]<<endl;
			fs6.close();

			double ep=1e-5;
			double dX=1;
			double *X=new double [h_num];
			double *old_X=new double [h_num];
			double *fx=new double [h_num];
			double *DfDx=new double [h_num*h_num];

			int count=0;

			for(int i=0;i<h_num;i++)
			{
				X[i]=1;
				old_X[i]=1;
				fx[i]=0;
				for(int j=0;j<h_num;j++)	DfDx[i*h_num+j]=0;
			}

			while(dX>ep)
			{
				count++;

				int way=1;
				//�P���ȌJ��Ԃ��v�Z
				/*
				for(int i=0;i<h_num;i++)
				{
					old_X[i]=HYPER[i].lambda;

					HYPER[i].half_p[A_X]+=HYPER[i].p[A_X]-pn_x[i];
					HYPER[i].half_p[A_Y]+=HYPER[i].p[A_Y]-pn_y[i];
					HYPER[i].half_p[A_Z]+=HYPER[i].p[A_Z]-pn_z[i];
	
					pn_x[i]=HYPER[i].p[A_X];
					pn_y[i]=HYPER[i].p[A_Y];
					pn_z[i]=HYPER[i].p[A_Z];
				}

				calc_differential_p(CON,PART,HYPER,HYPER1,F);
				renew_lambda(CON,PART,HYPER,HYPER1,t);
				calc_half_p(CON,PART,HYPER,HYPER1,1,F);//*/
			
				//Newton-Raphson�@

				for(int i=0;i<h_num;i++)
				{
					old_X[i]=X[i];
					fx[i]=0;
					for(int j=0;j<h_num;j++)	DfDx[i*h_num+j]=0;
				}

				calc_newton_function2(CON,HYPER,HYPER1,old_hpz,X,fx,DfDx,h_num,count,t,F);
				gauss(DfDx,fx,h_num);
				//double ep=CON.get_FEMCGep();
				///*GaussSeidelvh(DfDx,h_num,fx,ep);

				for(int i=0;i<h_num;i++)	X[i]-=fx[i];//*0.5*mi/(Dt*Dt)*V*fx[i];//

				dX=0;
				double dX_d=0;
				for(int i=0;i<h_num;i++)
				{
					dX_d+=fabs(old_X[i]);
					dX+=fabs(X[i]-old_X[i]);
				}
				dX/=dX_d;

				//�o��
				stringstream ss_E;
				ss_E<<"./Wall/E"<<t<<".csv";	
				stringstream ss_lam;
				ss_lam<<"./Wall/lambda"<<t<<".csv";		
	
				if(count==1)
				{
					ofstream init0(ss_E.str(), ios::trunc);
					ofstream init1(ss_lam.str(), ios::trunc);
					init0.close();
					init1.close();
				}
				ofstream e(ss_E.str(), ios::app);
				ofstream lam(ss_lam.str(), ios::app);

				e<<count<<","<<dX<<endl;
				lam<<count<<",";
				double sum_lam=0;	
				for(int i=0;i<h_num;i++)
				{
					lam<<X[i]<<",";	//����
					sum_lam+=X[i];
				}
				lam<<sum_lam/h_num<<","<<h_num<<endl;	//����
				e.close();
				lam.close();

				if(count%100==1)	cout<<"count"<<count<<" ,E"<<dX<<endl;
				if(count>1000)	break;
			}
		
			ofstream fs2("lambda_n2.csv",ios::app);

			for(int i=0;i<h_num;i++)
			{
				fs2<<X[i]<<",";	//������X[i]�����ׂăG���[
				HYPER[i].lambda=X[i];
			}
			fs2<<endl;
			fs2.close();

			ofstream fsw("Convergence_rate.csv", ios::app);
			fsw<<count<<","<<dX<<endl;
			fsw.close();
	
			delete[]	X;
			delete[]	old_X;
			delete[]	fx;
			delete[]	DfDx;
		}//*/
		delete[]	old_hpz;
	
	}
	
		calc_half_p(CON,PART,HYPER,HYPER1,1,F);

	//	for(int i=0;i<p_num;i++)	cout<<"renew_p_x"<<i<<"="<<HYPER[i].p[A_X]<<endl;
	cout<<"Hypercalculation ends."<<endl;

	clock_t end_t=clock();
	time<<end_t*CLOCKS_PER_SEC<<"	";
	time.close();
}

void calc_constant(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1)
{
//	cout<<"�����l�v�Z";

	double le=CON.get_distancebp();
	double r=CON.get_h_dis();
	double Dt=CON.get_dt();
	double V=get_volume(&CON);	//�l�����K�v��������Ȃ�
	double mi=V*CON.get_hyper_density();
	int h_num=HYPER.size();
	int model=CON.get_model_number();
//	cout<<"V"<<V<<endl;

	ofstream fc("constant.dat");
	fc<<"le"<<","<<le<<endl;
	fc<<"r"<<","<<r<<endl;
	fc<<"Dt"<<","<<Dt<<endl;
	fc<<"Volume"<<","<<V<<endl;
	fc<<"mi"<<","<<mi<<endl;
	fc.close();


	//��ʂ�Z���W�Ɨ��q���̒T��
	////�����^����
	//if(CON.get_flag_G()==OFF && CON.get_model_number()!=31)	for(int i=0;i<h_num;i++)	HYPER[i].p[A_Z]=-1.0*mi;
	
	//�Ȃ��˂���
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
				HYPER[i].pnd+=wiin;
				N++;
			}
			else	wiin=0;

			HYPER1[i*h_num+j].wiin=wiin;
		}
		HYPER[i].N=N;
		HYPER[i].pnd0=N;
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
//		cout<<"J["<<i<<"]="<<HYPER[i].J<<endl;

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
//			cout<<"DgDq["<<k<<","<<i<<"]="<<HYPER1[k*h_num+i].DgDq[A_X]<<","<<HYPER1[k*h_num+i].DgDq[A_Y]<<","<<HYPER1[k*h_num+i].DgDq[A_Z]<<endl;
		}
		HYPER1[i*h_num+i].DgDq[A_X]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_X][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_X][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_X][2]*HYPER1[i*h_num+i].n0ij[2]);
		HYPER1[i*h_num+i].DgDq[A_Y]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_Y][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_Y][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_Y][2]*HYPER1[i*h_num+i].n0ij[2]);
		HYPER1[i*h_num+i].DgDq[A_Z]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_Z][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_Z][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_Z][2]*HYPER1[i*h_num+i].n0ij[2]);
//		cout<<"DgDq["<<i<<","<<i<<"]="<<HYPER1[i*h_num+i].DgDq[A_X]<<","<<HYPER1[i*h_num+i].DgDq[A_Y]<<","<<HYPER1[i*h_num+i].DgDq[A_Z]<<endl;
	}	
//	cout<<"----------OK"<<endl;


 }


/////�j���[�g�����t�\���@ 
void newton_raphson(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int t,double **F)
{
	/////fx(N*N�s��̊e����)��1�����z��Ŋi�[�A(i,j)�����Ȃ�[j*N+i]�ŎQ��
	/////DfDx(N*N�s��̊e����)��1�����z��Ŋi�[�A(i,j)�����Ȃ�[j*N+i]�ŎQ��

	int calc_type=1;//�j���[�g���t�\���̔������@ 0:�Δ������̋t�s������̂܂܋��߂�@1:���`�������𗘗p

	//OPENMP���I���ł���Ώ������܂��Bthreads��pc�ň�����ő�X���b�h��
	#ifdef _OPENMP
    printf("OpenMP : On, threads = %d\n", omp_get_max_threads());
	#endif

	//�ő�X���b�h���Ōv�Z��������Ɖߕ��ׂ�CPU�����ɔM���Ȃ邽�߁A���񉻐����w�肵�Ă���
	//�Ȃ��A�ő�X���b�h��12�̂Ƃ���8�`10���炢���ڈ�
	omp_set_num_threads(8);
	//pn=2;//test,�Ƃ肠����2���łƂ��邩�ǂ����m�F 
	//////////////////�@f1(x1,x2) = x1^2 + x2^2 -5 = 0 f2(x1,x2) = x1^2/9+ x2^2 -1 = 0  http://homepage1.nifty.com/gfk/excel_newton_ren.htm

	int h_num=HYPER.size();
	double *fx=new double [h_num];//�֐��l�B
	double *DfDx=new double [h_num*h_num];//�֐��̕Δ����l�B
	double *XX=new double [h_num];//���݂̉��B	
	double *XX_old=new double [h_num];//1�X�e�b�v�O�̉��B
	double ep=1e-5;//��������
	double E=1;//���݂̌덷
	int count=0;//������
	double d;
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	double Dt=CON.get_interval();
	double sum=0;
	double E_old=0;
	int dec_flag=ON;

	#pragma omp parallel for
	for(int i=0;i<h_num;i++)
	{
		XX[i]=1;
		XX_old[i]=0;
		fx[i]=0;
		for(int j=0;j<h_num;j++)	DfDx[i*h_num+j]=0;
	}
	ofstream time("time_log.dat",ios::app);
	clock_t start_t=clock();
	time<<start_t*CLOCKS_PER_SEC<<"	";

	//	for(int i=0; i<N; i++) XX[i]=1;///�����l��^����B�Ƃ肠����1��
	cout<<"NR�@�J�n";
	while(E>ep)
	{
		count++;
		for(int i=0; i<h_num; i++)
		{
			XX_old[i]=XX[i];	//�����L��
			//HYPER[i].flag_wall=OFF;			
		}


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
			//double ep=CON.get_FEMCGep();
			///*GaussSeidelvh(DfDx,h_num,fx,ep);
			for(int i=0;i<h_num;i++)	XX[i]-=fx[i];//*0.5*mi/(Dt*Dt)*V*fx[i];//*/
		}

		//�덷�̕]��
		E_old=E;
		E=0;
		sum=0;
		for(int i=0; i<h_num; i++)
		{
			E+=fabs(XX[i]-XX_old[i]);
			//sum+=fabs(XX[i]);	//��Ό덷�ŕ]��
		}
		//E/=sum;


		if(count==1 || count%500==0)
		{		
			/*
			cout<<"XX_old["<<i<<"]-d=X["<<i<<"]	";
			cout<<XX_old[i]<<" - ";
			cout<<d<<" = ";
			cout<<XX[i]<<endl;*/

			cout<<"������	"<<count<<" E="<<E<<endl;
			output_newton_data2(E,XX,h_num,count,t);

		}
		if(count>CON.get_nr())	break;
		else if(dec_flag==ON)	if(E_old-E<0)	break;	
	}

	ofstream fs("Newton_Convergence_rate.csv", ios::app);
	fs<<count<<","<<E<<endl;
	fs.close();

//	end=clock();
//	newton_t=(end-start)/CLOCKS_PER_SEC;

	cout<<"��������";
	#pragma omp parallel for
	for(int i=0;i<h_num;i++) HYPER[i].lambda=XX[i];

//	for(int i=0;i<h_num;i++)	cout<<"lambda["<<i<<"]="<<HYPER[i].lambda<<endl;
	delete[]	fx;
	delete[]	DfDx;
	delete[]	XX;
	delete[]	XX_old;
	cout<<"---------- OK"<<endl;
}

void newton_raphson2(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int t,double **F)
{
	/////fx(N*N�s��̊e����)��1�����z��Ŋi�[�A(i,j)�����Ȃ�[j*N+i]�ŎQ��
	/////DfDx(N*N�s��̊e����)��1�����z��Ŋi�[�A(i,j)�����Ȃ�[j*N+i]�ŎQ��

	int calc_type=1;//�j���[�g���t�\���̔������@ 0:�Δ������̋t�s������̂܂܋��߂�@1:���`�������𗘗p

	//OPENMP���I���ł���Ώ������܂��Bthreads��pc�ň�����ő�X���b�h��
	#ifdef _OPENMP
    printf("OpenMP : On, threads = %d\n", omp_get_max_threads());
	#endif

	//�ő�X���b�h���Ōv�Z��������Ɖߕ��ׂ�CPU�����ɔM���Ȃ邽�߁A���񉻐����w�肵�Ă���
	//�Ȃ��A�ő�X���b�h��12�̂Ƃ���8�`10���炢���ڈ�
	omp_set_num_threads(8);
	//pn=2;//test,�Ƃ肠����2���łƂ��邩�ǂ����m�F 
	//////////////////�@f1(x1,x2) = x1^2 + x2^2 -5 = 0 f2(x1,x2) = x1^2/9+ x2^2 -1 = 0  http://homepage1.nifty.com/gfk/excel_newton_ren.htm

	int h_num=HYPER.size();
	double *fx=new double [h_num];//�֐��l�B
	double *DfDx=new double [h_num*h_num];//�֐��̕Δ����l�B
	double *XX=new double [h_num];//���݂̉��B	
	double *XX_old=new double [h_num];//1�X�e�b�v�O�̉��B
	double ep=1e-5;//��������
	double E=1;//���݂̌덷
	int count=0;//������
	double d;
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	double Dt=CON.get_interval();
	double sum=0;
	double E_old=0;
	int dec_flag=ON;

	#pragma omp parallel for
	for(int i=0;i<h_num;i++)
	{
		XX[i]=1;
		XX_old[i]=0;
		fx[i]=0;
		for(int j=0;j<h_num;j++)	DfDx[i*h_num+j]=0;
	}
	ofstream time("time_log.dat",ios::app);
	clock_t start_t=clock();
	time<<start_t*CLOCKS_PER_SEC<<"	";

	//	for(int i=0; i<N; i++) XX[i]=1;///�����l��^����B�Ƃ肠����1��
	cout<<"NR�@�J�n";
	while(E>ep)
	{
		count++;
		for(int i=0; i<h_num; i++)
		{
			XX_old[i]=XX[i];	//�����L��
			//HYPER[i].flag_wall=OFF;			
		}


//		if(count==1)	for(int i=0;i<h_num;i++)	for(int j=0;j<h_num;j++)	for(int D=0;D<DIMENSION;D++)	HYPER1[i*h_num+j].newton_DgDq[D]=HYPER1[i*N+j].DgDq[D];

		calc_newton_function3(CON,PART,HYPER,HYPER1,XX,fx,DfDx,h_num,count,t,F);


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
			//double ep=CON.get_FEMCGep();
			///*GaussSeidelvh(DfDx,h_num,fx,ep);
			for(int i=0;i<h_num;i++)	XX[i]-=fx[i];//*0.5*mi/(Dt*Dt)*V*fx[i];//*/
		}

		//�덷�̕]��
		E_old=E;
		E=0;
		sum=0;
		for(int i=0; i<h_num; i++)
		{
			E+=fabs(XX[i]-XX_old[i]);
			//sum+=fabs(XX[i]);	//��Ό덷�ŕ]��
		}
		//E/=sum;


//		if(count==1 || count%2500==0)
		{		
			/*
			cout<<"XX_old["<<i<<"]-d=X["<<i<<"]	";
			cout<<XX_old[i]<<" - ";
			cout<<d<<" = ";
			cout<<XX[i]<<endl;*/

			cout<<"������	"<<count<<" E="<<E<<endl;
			output_newton_data2(E,XX,h_num,count,t);

		}
		if(count>CON.get_nr())	break;
		else if(dec_flag==ON)	if(E_old-E<0)	break;	
	}

	ofstream fs("Newton_Convergence_rate.csv", ios::app);
	fs<<count<<","<<E<<endl;
	fs.close();

//	end=clock();
//	newton_t=(end-start)/CLOCKS_PER_SEC;

	cout<<"��������";
	#pragma omp parallel for
	for(int i=0;i<h_num;i++) HYPER[i].lambda=XX[i];

//	for(int i=0;i<h_num;i++)	cout<<"lambda["<<i<<"]="<<HYPER[i].lambda<<endl;
	delete[]	fx;
	delete[]	DfDx;
	delete[]	XX;
	delete[]	XX_old;
	cout<<"---------- OK"<<endl;
}

void calc_newton_function(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,double *lambda,double *fx,double *DfDx,int hyper_number,int count,int t,double **F)
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
	double G=9.8;
	double *n_rx=new double[h_num];
	double *n_ry=new double[h_num];
	double *n_rz=new double[h_num];

	double **n_DgDq_x=new double *[h_num];
	double **n_DgDq_y=new double *[h_num];
	double **n_DgDq_z=new double *[h_num];
	for(int i=0;i<h_num;i++)
	{
		n_rx[i]=PART[i].r[A_X];
		n_ry[i]=PART[i].r[A_Y];
		n_rz[i]=PART[i].r[A_Z];
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
/*		if(model_num==30||model_num==23)
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
		else*/
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
			if(flag_G==ON)	p_half_p[A_Z]-=G*mi;
			//�S�����̉e��
			if(flag_vis==ON)
			{
				p_half_p[A_X]+=HYPER[i].vis_force[A_X];
				p_half_p[A_Y]+=HYPER[i].vis_force[A_Y];
				p_half_p[A_Z]+=HYPER[i].vis_force[A_Z];
			}
			//����̍l��
			if(flag_FEM==ON || PART[i].toFEM==ON)
			{
				p_half_p[A_X]+=F[A_X][i]*V;//density;
				p_half_p[A_Y]+=F[A_Y][i]*V;//density;
				p_half_p[A_Z]+=F[A_Z][i]*V;//density;
			}
//		if(n_rz[i]>0)	//way2
		{
			//�ʒu���W�̌v�Z
			n_rx[i]+=Dt*(HYPER[i].p[A_X]+Dt*0.5*p_half_p[A_X])/mi;
			n_ry[i]+=Dt*(HYPER[i].p[A_Y]+Dt*0.5*p_half_p[A_Y])/mi;
			n_rz[i]+=Dt*(HYPER[i].p[A_Z]+Dt*0.5*p_half_p[A_Z])/mi;
//			cout<<"r["<<i<<"]="<<n_rx[i]<<", "<<n_ry[i]<<", "<<n_rz[i]<<endl;
			if(CON.get_flag_wall()==ON)
			{
				if(n_rz[i]<0)
				{
					n_rz[i]=0;
				}			
			}
	//		if(n_rz[i]<=0)	n_rz[i]*=-1;
		}
/*		else
		{
			n_rx[i]+=Dt*(HYPER[i].p[A_X]+Dt*0.5*p_half_p[A_X])/mi;	//way2
			n_ry[i]+=Dt*(HYPER[i].p[A_Y]+Dt*0.5*p_half_p[A_Y])/mi;	//way2
			cout<<"�ǐN�����q"<<i<<"lambda="<<lambda[i]<<endl;	//way1
		}*/
	}



	////DgDq��fx�̍X�V
	int p_num=PART.size();
	double r=CON.get_h_dis();
	for(int i=0;i<h_num;i++)
	{
		//Fi�̌v�Z
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
			fi[0][0]+=w*(n_rx[inn]-n_rx[i])*HYPER1[i*h_num+inn].aiin[A_X];
			fi[0][1]+=w*(n_rx[inn]-n_rx[i])*HYPER1[i*h_num+inn].aiin[A_Y];
			fi[0][2]+=w*(n_rx[inn]-n_rx[i])*HYPER1[i*h_num+inn].aiin[A_Z];
			fi[1][0]+=w*(n_ry[inn]-n_ry[i])*HYPER1[i*h_num+inn].aiin[A_X];
			fi[1][1]+=w*(n_ry[inn]-n_ry[i])*HYPER1[i*h_num+inn].aiin[A_Y];
			fi[1][2]+=w*(n_ry[inn]-n_ry[i])*HYPER1[i*h_num+inn].aiin[A_Z];
			fi[2][0]+=w*(n_rz[inn]-n_rz[i])*HYPER1[i*h_num+inn].aiin[A_X];
			fi[2][1]+=w*(n_rz[inn]-n_rz[i])*HYPER1[i*h_num+inn].aiin[A_Y];
			fi[2][2]+=w*(n_rz[inn]-n_rz[i])*HYPER1[i*h_num+inn].aiin[A_Z];

/*			double dis=sqrt((n_rx[inn]-n_rx[i])*(n_rx[inn]-n_rx[i])+(n_ry[inn]-n_ry[i])*(n_ry[inn]-n_ry[i])+(n_rz[inn]-n_rz[i])*(n_rz[inn]-n_rz[i]));
			pnd+=kernel4(r,dis);*/
/*			if(n_rz[inn]<0)
			{
				//double dis=sqrt((n_rx[inn]-n_rx[i])*(n_rx[inn]-n_rx[i])+(n_ry[inn]-n_ry[i])*(n_ry[inn]-n_ry[i])+(n_rz[inn]-n_rz[i])*(n_rz[inn]-n_rz[i]));
				double dis=sqrt(HYPER1[i*h_num+inn].aiin[A_X]*HYPER1[i*h_num+inn].aiin[A_X]+HYPER1[i*h_num+inn].aiin[A_Y]*HYPER1[i*h_num+inn].aiin[A_Y]+HYPER1[i*h_num+inn].aiin[A_Z]*HYPER1[i*h_num+inn].aiin[A_Z]);
				
				double d_wij=-r/dis/dis;
//				double wij=kernel4(r,dis);
				part_g+=w/HYPER[i].pnd;

				n_DgDq_x[i][inn]+=V/HYPER[i].pnd*d_wij*HYPER1[i*h_num+inn].aiin[A_X]/dis;
				n_DgDq_y[i][inn]+=V/HYPER[i].pnd*d_wij*HYPER1[i*h_num+inn].aiin[A_Y]/dis;
				n_DgDq_z[i][inn]+=V/HYPER[i].pnd*d_wij*HYPER1[i*h_num+inn].aiin[A_Z]/dis;

				n_DgDq_x[i][i]+=-V/HYPER[i].pnd*d_wij*HYPER1[i*h_num+inn].aiin[A_X]/dis;
				n_DgDq_y[i][i]+=-V/HYPER[i].pnd*d_wij*HYPER1[i*h_num+inn].aiin[A_Y]/dis;
				n_DgDq_z[i][i]+=-V/HYPER[i].pnd*d_wij*HYPER1[i*h_num+inn].aiin[A_Z]/dis;
			}*/
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

//		cout<<"F"<<i<<"="<<p_Fi[0][0]<<", "<<p_Fi[0][1]<<", "<<p_Fi[0][2]<<endl;
//		cout<<p_Fi[1][0]<<", "<<p_Fi[1][1]<<", "<<p_Fi[1][2]<<endl;
//		cout<<p_Fi[2][0]<<", "<<p_Fi[2][1]<<", "<<p_Fi[2][2]<<endl;

		//J�̌v�Z
		double J=calc_det3(p_Fi);
//		cout<<"J"<<i<<"="<<J<<endl;

/*		double part_dgdq_x=0;
		double part_dgdq_y=0;
		double part_dgdq_z=0;
		double part_g=0;
		for(int k=h_num;k<p_num;k++)
		{
			double dis=sqrt((PART[k].r[A_X]-n_rx[i])*(PART[k].r[A_X]-n_rx[i])+(PART[k].r[A_Y]-n_ry[i])*(PART[k].r[A_Y]-n_ry[i])+(PART[k].r[A_Z]-n_rz[i])*(PART[k].r[A_Z]-n_rz[i]));
			if(dis<r)
			{
				double wij=kernel4(r,dis);
				double d_wij=-r/dis/dis;
				part_g+=wij/HYPER[i].pnd;
				part_dgdq_x+=-V/HYPER[i].pnd*d_wij*(PART[k].r[A_X]-n_rx[i])/dis;
				part_dgdq_y+=-V/HYPER[i].pnd*d_wij*(PART[k].r[A_Y]-n_ry[i])/dis;
				part_dgdq_z+=-V/HYPER[i].pnd*d_wij*(PART[k].r[A_Z]-n_rz[i])/dis;
			}
		}//*/
		//fx�̌v�Z
//		fx[i]=V*(1-J+part_g);//1-J;//
		fx[i]=V*(1-J);//1-J;//
//		cout<<"fx"<<i<<"="<<V*(1-J)<<endl;

		//t_inverse_Fi�̌v�Z
		inverse(p_Fi,DIMENSION);

		//DgDq�̌v�Z
		for(int j=0;j<Ni;j++)
		{
			int k=HYPER[i].NEI[j];
			n_DgDq_x[i][k]=J*(p_Fi[0][0]*HYPER1[k*h_num+i].n0ij[0]+p_Fi[1][0]*HYPER1[k*h_num+i].n0ij[1]+p_Fi[2][0]*HYPER1[k*h_num+i].n0ij[2]);
			n_DgDq_y[i][k]=J*(p_Fi[0][1]*HYPER1[k*h_num+i].n0ij[0]+p_Fi[1][1]*HYPER1[k*h_num+i].n0ij[1]+p_Fi[2][1]*HYPER1[k*h_num+i].n0ij[2]);
			n_DgDq_z[i][k]=J*(p_Fi[0][2]*HYPER1[k*h_num+i].n0ij[0]+p_Fi[1][2]*HYPER1[k*h_num+i].n0ij[1]+p_Fi[2][2]*HYPER1[k*h_num+i].n0ij[2]);
//		cout<<"n_DgDq["<<i<<"*p_num+"<<k<<"]={"<<n_DgDq_x[i][k]<<","<<n_DgDq_y[i][k]<<","<<n_DgDq_z[i][k]<<"}"<<endl;

		}
		n_DgDq_x[i][i]=J*(p_Fi[0][0]*HYPER1[i*h_num+i].n0ij[0]+p_Fi[1][0]*HYPER1[i*h_num+i].n0ij[1]+p_Fi[2][0]*HYPER1[i*h_num+i].n0ij[2]);
		n_DgDq_y[i][i]=J*(p_Fi[0][1]*HYPER1[i*h_num+i].n0ij[0]+p_Fi[1][1]*HYPER1[i*h_num+i].n0ij[1]+p_Fi[2][1]*HYPER1[i*h_num+i].n0ij[2]);
		n_DgDq_z[i][i]=J*(p_Fi[0][2]*HYPER1[i*h_num+i].n0ij[0]+p_Fi[1][2]*HYPER1[i*h_num+i].n0ij[1]+p_Fi[2][2]*HYPER1[i*h_num+i].n0ij[2]);

//		cout<<"n_DgDq["<<i<<"*p_num+"<<i<<"]={"<<n_DgDq_x[i][i]<<","<<n_DgDq_y[i][i]<<","<<n_DgDq_z[i][i]<<"}"<<endl;
//		cout<<"fx["<<i<<"]="<<fx[i]<<endl;
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

	/*
	double *part_fx=new double [h_num];
	double *w_DfDx=new double [h_num*h_num];
	double *w_fx=new double [h_num];
	int *Nw=new int [h_num];

	for(int i=0;i<h_num;i++)
	{
		for(int j=0;j<h_num;j++)	w_DfDx[i*h_num+j]=0;
		w_fx[i]=0;
		Nw[i]=0;
	}

	int number=0;
	for(int i=0;i<h_num;i++)
	{
		part_fx[i]=0;
		if(n_rz[i]<0)
		{
			int Ni=HYPER[i].N;
			DfDx[i*h_num+i]=1;
			for(int j=0;j<Ni;j++)
			{
				Nw[number]=i;
				int jn=HYPER[i].NEI[j];
				w_DfDx[i*h_num+jn]=Dt/2*HYPER1[jn*h_num+i].DgDq[A_Z];
				part_fx[i]+=Dt/2*(HYPER[jn].stress[A_Z][A_X]*HYPER1[jn*h_num+i].DgDq[A_X]+HYPER[jn].stress[A_Z][A_Y]*HYPER1[jn*h_num+i].DgDq[A_Y]+HYPER[jn].stress[A_Z][A_Z]*HYPER1[jn*h_num+i].DgDq[A_Z]);
				DfDx[i*h_num+jn]=0;
			}
			w_fx[i]=HYPER[i].p[A_Z]+part_fx[i];
			number++;
		}
//		w_DfDx[i*h_num+i]=1;
	}
	if(number>0)
	{

		output_newton_data3(w_fx,w_DfDx,n_rx,n_ry,n_rz,part_fx,h_num,count,t);		
		gauss(w_DfDx,w_fx,h_num);

		for(int i=0;i<number;i++)
		{
			int in=Nw[i];
			fx[in]=w_fx[in];
		}
	}
	delete[]	Nw;
	delete[]	w_DfDx;
	delete[]	w_fx;
	delete[]	part_fx;//*/

	output_newton_data1(fx,DfDx,n_rx,n_ry,n_rz,h_num,count,t);
	



	////�o��
//	if(count%200==0 && count>CON.get_nr()/2)
//	if(t==1||t%CON.get_interval()==0)	if(count%200==0||count==1)	
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

void calc_newton_function2(mpsconfig &CON,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,double *old_hpz, double *lambda, double *fx,double *DfDx,int hyper_number,int count,int t,double **F)
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
	double G=9.8;



	////�ʒu���W�̍X�V	
	for(int i=0;i<h_num;i++)
	{
		//half_p�̌v�Z
		double part_p[DIMENSION]={0,0,0};

		int Ni=HYPER[i].N;
		for(int j=0;j<Ni;j++)
		{	
			int k=HYPER[i].NEI[j];
			part_p[A_X]+=(HYPER[k].stress[0][0]-lambda[k])*HYPER1[k*h_num+i].DgDq[0]+HYPER[k].stress[0][1]*HYPER1[k*h_num+i].DgDq[1]+HYPER[k].stress[0][2]*HYPER1[k*h_num+i].DgDq[2];
			part_p[A_Y]+=HYPER[k].stress[1][0]*HYPER1[k*h_num+i].DgDq[0]+(HYPER[k].stress[1][1]-lambda[k])*HYPER1[k*h_num+i].DgDq[1]+HYPER[k].stress[1][2]*HYPER1[k*h_num+i].DgDq[2];
			part_p[A_Z]+=HYPER[k].stress[2][0]*HYPER1[k*h_num+i].DgDq[0]+HYPER[k].stress[2][1]*HYPER1[k*h_num+i].DgDq[1]+(HYPER[k].stress[2][2]-lambda[k])*HYPER1[k*h_num+i].DgDq[2];
		}//j�Ɋւ���for���̏I���
		part_p[A_X]+=(HYPER[i].stress[0][0]-lambda[i])*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[0][1]*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[0][2]*HYPER1[i*h_num+i].DgDq[2];
		part_p[A_Y]+=HYPER[i].stress[1][0]*HYPER1[i*h_num+i].DgDq[0]+(HYPER[i].stress[1][1]-lambda[i])*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[1][2]*HYPER1[i*h_num+i].DgDq[2];
		part_p[A_Z]+=HYPER[i].stress[2][0]*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[2][1]*HYPER1[i*h_num+i].DgDq[1]+(HYPER[i].stress[2][2]-lambda[i])*HYPER1[i*h_num+i].DgDq[2];
		
		
		//�d�͂̉e��
		if(flag_G==ON)	part_p[A_Z]-=G*mi;
		//�S�����̉e��
		if(flag_vis==ON)
		{
			part_p[A_X]+=HYPER[i].vis_force[A_X];
			part_p[A_Y]+=HYPER[i].vis_force[A_Y];
			part_p[A_Z]+=HYPER[i].vis_force[A_Z];
		}
		//����̍l��
		if(flag_FEM==ON)
		{
			part_p[A_X]+=F[A_X][i]*V;//density;
			part_p[A_Y]+=F[A_Y][i]*V;//density;
			part_p[A_Z]+=F[A_Z][i]*V;//density;
		}
		//�ʒu���W�̌v�Z
		HYPER[i].p[A_X]=HYPER[i].half_p[A_X]+Dt*0.5*part_p[A_X];
		HYPER[i].p[A_Y]=HYPER[i].half_p[A_Y]+Dt*0.5*part_p[A_Y];
		HYPER[i].p[A_Z]=HYPER[i].half_p[A_Z]+Dt*0.5*part_p[A_Z];	
//		HYPER[i].p[A_Z]=old_hpz[i]+Dt*0.5*part_p[A_Z];	
//			cout<<"r["<<i<<"]="<<n_rx[i]<<", "<<n_ry[i]<<", "<<n_rz[i]<<endl;
	}

	////DgDq��fx�̍X�V
	double r=CON.get_h_dis();


	//calc_newton1�Ɠ������@
	/*for(int k=0;k<h_num;k++)
	{
		int Nk=HYPER[k].N;
		for(int i=0;i<Nk;i++)
		{
			int in=HYPER[k].NEI[i];
			for(int j=0;j<Nk;j++)
			{
				int jn=HYPER[k].NEI[j];
				DfDx[in*h_num+jn]-=Dt*0.5/mi*(HYPER1[in*h_num+k].DgDq[A_X]*HYPER1[jn*h_num+k].DgDq[A_X]+HYPER1[in*h_num+k].DgDq[A_Y]*HYPER1[jn*h_num+k].DgDq[A_Y]+HYPER1[in*h_num+k].DgDq[A_Z]*HYPER1[jn*h_num+k].DgDq[A_Z]);
			}
			DfDx[k*h_num+in]-=Dt*0.5/mi*(HYPER1[k*h_num+k].DgDq[A_X]*HYPER1[in*h_num+k].DgDq[A_X]+HYPER1[k*h_num+k].DgDq[A_Y]*HYPER1[in*h_num+k].DgDq[A_Y]+HYPER1[k*h_num+k].DgDq[A_Z]*HYPER1[in*h_num+k].DgDq[A_Z]);
			DfDx[in*h_num+k]-=Dt*0.5/mi*(HYPER1[in*h_num+k].DgDq[A_X]*HYPER1[k*h_num+k].DgDq[A_X]+HYPER1[in*h_num+k].DgDq[A_Y]*HYPER1[k*h_num+k].DgDq[A_Y]+HYPER1[in*h_num+k].DgDq[A_Z]*HYPER1[k*h_num+k].DgDq[A_Z]);
		}
		DfDx[k*h_num+k]-=Dt*Dt*0.5/mi*(HYPER1[k*h_num+k].DgDq[A_X]*HYPER1[k*h_num+k].DgDq[A_X]+HYPER1[k*h_num+k].DgDq[A_Y]*HYPER1[k*h_num+k].DgDq[A_Y]+HYPER1[k*h_num+k].DgDq[A_Z]*HYPER1[k*h_num+k].DgDq[A_Z]);		
		for(int l=0;l<h_num;l++)	fx[k]+=1/mi*(HYPER1[k*h_num+l].DgDq[A_X]*HYPER[l].p[A_X]+HYPER1[k*h_num+l].DgDq[A_Y]*HYPER[l].p[A_Y]+HYPER1[k*h_num+l].DgDq[A_Z]*HYPER[l].p[A_Z]);
	}//*/

	//simplest	
	for(int i=0;i<h_num;i++)
	{
		for(int j=0;j<h_num;j++)
		{
			fx[i]+=HYPER1[i*h_num+j].DgDq[A_X]*HYPER[j].p[A_X]+HYPER1[i*h_num+j].DgDq[A_Y]*HYPER[j].p[A_Y]+HYPER1[i*h_num+j].DgDq[A_Z]*HYPER[j].p[A_Z];
			for(int k=0;k<h_num;k++)	DfDx[i*h_num+j]-=Dt*0.5/mi*(HYPER1[i*h_num+k].DgDq[A_X]*HYPER1[j*h_num+k].DgDq[A_X]+HYPER1[i*h_num+k].DgDq[A_Y]*HYPER1[j*h_num+k].DgDq[A_Y]+HYPER1[i*h_num+k].DgDq[A_Z]*HYPER1[j*h_num+k].DgDq[A_Z]);
		}
		if(fx[i]!=0)	fx[i]/=mi;
	}//*/

	output_newton2_data1(HYPER,fx,DfDx,h_num,count,t);
	

	////�o��
//	if(count%200==0 && count>CON.get_nr()/2)
//	if(t==1||t%CON.get_interval()==0)	if(count%200==0||count==1)	
//	if(t==1||t%(10*CON.get_interval())==0)	if(count%200==0||count==1)	output_newton_data1(fx,DfDx,n_rx,n_ry,n_rz,h_num,count,t);

	ofstream t_loge("time_log_newton.dat", ios::app);
	clock_t t4=clock();
	t_loge<<"step="<<t<<", count="<<count<<", time="<<1000*(t4-t3)/CLOCKS_PER_SEC<<"[e-3sec]"<<endl;
	t_loge.close();


}

void calc_newton_function3(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,double *lambda,double *fx,double *DfDx,int hyper_number,int count,int t,double **F)
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
	double G=9.8;
	
	double *n_rx=new double[h_num];
	double *n_ry=new double[h_num];
	double *n_rz=new double[h_num];

	double **n_DgDq_x=new double *[h_num];
	double **n_DgDq_y=new double *[h_num];
	double **n_DgDq_z=new double *[h_num];

	int *num_w=new int [h_num];

	for(int i=0;i<h_num;i++)
	{
		HYPER[i].flag_wall=OFF;
		n_rx[i]=PART[i].r[A_X];
		n_ry[i]=PART[i].r[A_Y];
		n_rz[i]=PART[i].r[A_Z];
		n_DgDq_x[i]=new double [h_num];
		n_DgDq_y[i]=new double [h_num];
		n_DgDq_z[i]=new double [h_num];
		num_w[i]=0;
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
	
	int Nw=0;

	double **p_Fi=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	p_Fi[D]=new double [DIMENSION];

	////�ʒu���W�̍X�V	
	for(int i=0;i<h_num;i++)
	{
/*		if(model_num==30||model_num==23)
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
		else*/
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
			if(flag_G==ON)	p_half_p[A_Z]-=G*mi;
			//�S�����̉e��
			if(flag_vis==ON)
			{
				p_half_p[A_X]+=HYPER[i].vis_force[A_X];
				p_half_p[A_Y]+=HYPER[i].vis_force[A_Y];
				p_half_p[A_Z]+=HYPER[i].vis_force[A_Z];
			}
			//����̍l��
			if(flag_FEM==ON || PART[i].toFEM==ON)
			{
				p_half_p[A_X]+=F[A_X][i]*V;//density;
				p_half_p[A_Y]+=F[A_Y][i]*V;//density;
				p_half_p[A_Z]+=F[A_Z][i]*V;//density;
			}
//		if(n_rz[i]>0)	//way2
		{
			//�ʒu���W�̌v�Z
			n_rx[i]+=Dt*(HYPER[i].p[A_X]+Dt*0.5*p_half_p[A_X])/mi;
			n_ry[i]+=Dt*(HYPER[i].p[A_Y]+Dt*0.5*p_half_p[A_Y])/mi;
			n_rz[i]+=Dt*(HYPER[i].p[A_Z]+Dt*0.5*p_half_p[A_Z])/mi;
//			cout<<"r["<<i<<"]="<<n_rx[i]<<", "<<n_ry[i]<<", "<<n_rz[i]<<endl;
			if(CON.get_flag_wall()==ON)
			{
				if(n_rz[i]<0)
				{
					n_rz[i]=0;
					HYPER[i].flag_wall=ON;
				//	fx[i]=PART[i].r[A_Z]+Dt*(HYPER[i].p[A_Z]+Dt*0.5*p_half_p[A_Z])/mi;
					num_w[Nw]=i;
					Nw++;
				}			
			}
	//		if(n_rz[i]<=0)	n_rz[i]*=-1;
		}
/*		else
		{
			n_rx[i]+=Dt*(HYPER[i].p[A_X]+Dt*0.5*p_half_p[A_X])/mi;	//way2
			n_ry[i]+=Dt*(HYPER[i].p[A_Y]+Dt*0.5*p_half_p[A_Y])/mi;	//way2
			cout<<"�ǐN�����q"<<i<<"lambda="<<lambda[i]<<endl;	//way1
		}*/
	}



	////DgDq��fx�̍X�V
	int p_num=PART.size();
	double r=CON.get_h_dis();
	for(int i=0;i<h_num;i++)
	{
		//Fi�̌v�Z
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
			fi[0][0]+=w*(n_rx[inn]-n_rx[i])*HYPER1[i*h_num+inn].aiin[A_X];
			fi[0][1]+=w*(n_rx[inn]-n_rx[i])*HYPER1[i*h_num+inn].aiin[A_Y];
			fi[0][2]+=w*(n_rx[inn]-n_rx[i])*HYPER1[i*h_num+inn].aiin[A_Z];
			fi[1][0]+=w*(n_ry[inn]-n_ry[i])*HYPER1[i*h_num+inn].aiin[A_X];
			fi[1][1]+=w*(n_ry[inn]-n_ry[i])*HYPER1[i*h_num+inn].aiin[A_Y];
			fi[1][2]+=w*(n_ry[inn]-n_ry[i])*HYPER1[i*h_num+inn].aiin[A_Z];
			fi[2][0]+=w*(n_rz[inn]-n_rz[i])*HYPER1[i*h_num+inn].aiin[A_X];
			fi[2][1]+=w*(n_rz[inn]-n_rz[i])*HYPER1[i*h_num+inn].aiin[A_Y];
			fi[2][2]+=w*(n_rz[inn]-n_rz[i])*HYPER1[i*h_num+inn].aiin[A_Z];

/*			double dis=sqrt((n_rx[inn]-n_rx[i])*(n_rx[inn]-n_rx[i])+(n_ry[inn]-n_ry[i])*(n_ry[inn]-n_ry[i])+(n_rz[inn]-n_rz[i])*(n_rz[inn]-n_rz[i]));
			pnd+=kernel4(r,dis);*/
/*			if(n_rz[inn]<0)
			{
				//double dis=sqrt((n_rx[inn]-n_rx[i])*(n_rx[inn]-n_rx[i])+(n_ry[inn]-n_ry[i])*(n_ry[inn]-n_ry[i])+(n_rz[inn]-n_rz[i])*(n_rz[inn]-n_rz[i]));
				double dis=sqrt(HYPER1[i*h_num+inn].aiin[A_X]*HYPER1[i*h_num+inn].aiin[A_X]+HYPER1[i*h_num+inn].aiin[A_Y]*HYPER1[i*h_num+inn].aiin[A_Y]+HYPER1[i*h_num+inn].aiin[A_Z]*HYPER1[i*h_num+inn].aiin[A_Z]);
				
				double d_wij=-r/dis/dis;
//				double wij=kernel4(r,dis);
				part_g+=w/HYPER[i].pnd;

				n_DgDq_x[i][inn]+=V/HYPER[i].pnd*d_wij*HYPER1[i*h_num+inn].aiin[A_X]/dis;
				n_DgDq_y[i][inn]+=V/HYPER[i].pnd*d_wij*HYPER1[i*h_num+inn].aiin[A_Y]/dis;
				n_DgDq_z[i][inn]+=V/HYPER[i].pnd*d_wij*HYPER1[i*h_num+inn].aiin[A_Z]/dis;

				n_DgDq_x[i][i]+=-V/HYPER[i].pnd*d_wij*HYPER1[i*h_num+inn].aiin[A_X]/dis;
				n_DgDq_y[i][i]+=-V/HYPER[i].pnd*d_wij*HYPER1[i*h_num+inn].aiin[A_Y]/dis;
				n_DgDq_z[i][i]+=-V/HYPER[i].pnd*d_wij*HYPER1[i*h_num+inn].aiin[A_Z]/dis;
			}*/
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

//		cout<<"F"<<i<<"="<<p_Fi[0][0]<<", "<<p_Fi[0][1]<<", "<<p_Fi[0][2]<<endl;
//		cout<<p_Fi[1][0]<<", "<<p_Fi[1][1]<<", "<<p_Fi[1][2]<<endl;
//		cout<<p_Fi[2][0]<<", "<<p_Fi[2][1]<<", "<<p_Fi[2][2]<<endl;

		//J�̌v�Z
		double J=calc_det3(p_Fi);
//		cout<<"J"<<i<<"="<<J<<endl;

/*		double part_dgdq_x=0;
		double part_dgdq_y=0;
		double part_dgdq_z=0;
		double part_g=0;
		for(int k=h_num;k<p_num;k++)
		{
			double dis=sqrt((PART[k].r[A_X]-n_rx[i])*(PART[k].r[A_X]-n_rx[i])+(PART[k].r[A_Y]-n_ry[i])*(PART[k].r[A_Y]-n_ry[i])+(PART[k].r[A_Z]-n_rz[i])*(PART[k].r[A_Z]-n_rz[i]));
			if(dis<r)
			{
				double wij=kernel4(r,dis);
				double d_wij=-r/dis/dis;
				part_g+=wij/HYPER[i].pnd;
				part_dgdq_x+=-V/HYPER[i].pnd*d_wij*(PART[k].r[A_X]-n_rx[i])/dis;
				part_dgdq_y+=-V/HYPER[i].pnd*d_wij*(PART[k].r[A_Y]-n_ry[i])/dis;
				part_dgdq_z+=-V/HYPER[i].pnd*d_wij*(PART[k].r[A_Z]-n_rz[i])/dis;
			}
		}//*/
		//fx�̌v�Z
//		fx[i]=V*(1-J+part_g);//1-J;//
		fx[i]=V*(1-J);//1-J;//
//		cout<<"fx"<<i<<"="<<V*(1-J)<<endl;

		//t_inverse_Fi�̌v�Z
		inverse(p_Fi,DIMENSION);

		//DgDq�̌v�Z
		for(int j=0;j<Ni;j++)
		{
			int k=HYPER[i].NEI[j];
			n_DgDq_x[i][k]=J*(p_Fi[0][0]*HYPER1[k*h_num+i].n0ij[0]+p_Fi[1][0]*HYPER1[k*h_num+i].n0ij[1]+p_Fi[2][0]*HYPER1[k*h_num+i].n0ij[2]);
			n_DgDq_y[i][k]=J*(p_Fi[0][1]*HYPER1[k*h_num+i].n0ij[0]+p_Fi[1][1]*HYPER1[k*h_num+i].n0ij[1]+p_Fi[2][1]*HYPER1[k*h_num+i].n0ij[2]);
			n_DgDq_z[i][k]=J*(p_Fi[0][2]*HYPER1[k*h_num+i].n0ij[0]+p_Fi[1][2]*HYPER1[k*h_num+i].n0ij[1]+p_Fi[2][2]*HYPER1[k*h_num+i].n0ij[2]);
//		cout<<"n_DgDq["<<i<<"*p_num+"<<k<<"]={"<<n_DgDq_x[i][k]<<","<<n_DgDq_y[i][k]<<","<<n_DgDq_z[i][k]<<"}"<<endl;

		}
		n_DgDq_x[i][i]=J*(p_Fi[0][0]*HYPER1[i*h_num+i].n0ij[0]+p_Fi[1][0]*HYPER1[i*h_num+i].n0ij[1]+p_Fi[2][0]*HYPER1[i*h_num+i].n0ij[2]);
		n_DgDq_y[i][i]=J*(p_Fi[0][1]*HYPER1[i*h_num+i].n0ij[0]+p_Fi[1][1]*HYPER1[i*h_num+i].n0ij[1]+p_Fi[2][1]*HYPER1[i*h_num+i].n0ij[2]);
		n_DgDq_z[i][i]=J*(p_Fi[0][2]*HYPER1[i*h_num+i].n0ij[0]+p_Fi[1][2]*HYPER1[i*h_num+i].n0ij[1]+p_Fi[2][2]*HYPER1[i*h_num+i].n0ij[2]);

//		cout<<"n_DgDq["<<i<<"*p_num+"<<i<<"]={"<<n_DgDq_x[i][i]<<","<<n_DgDq_y[i][i]<<","<<n_DgDq_z[i][i]<<"}"<<endl;
//		cout<<"fx["<<i<<"]="<<fx[i]<<endl;
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

	for(int i=0;i<Nw;i++)
	{
		int in=num_w[i];
		int Ni=HYPER[in].N;
		for(int k=0;k<Ni;k++)
		{
			int	kn=HYPER[in].NEI[k]; 
			DfDx[in*h_num+kn]=-Dt*Dt*0.5/mi*HYPER1[kn*h_num+in].DgDq[A_Z];
		}
	}

	/*
	double *part_fx=new double [h_num];
	double *w_DfDx=new double [h_num*h_num];
	double *w_fx=new double [h_num];
	int *Nw=new int [h_num];

	for(int i=0;i<h_num;i++)
	{
		for(int j=0;j<h_num;j++)	w_DfDx[i*h_num+j]=0;
		w_fx[i]=0;
		Nw[i]=0;
	}

	int number=0;
	for(int i=0;i<h_num;i++)
	{
		part_fx[i]=0;
		if(n_rz[i]<0)
		{
			int Ni=HYPER[i].N;
			DfDx[i*h_num+i]=1;
			for(int j=0;j<Ni;j++)
			{
				Nw[number]=i;
				int jn=HYPER[i].NEI[j];
				w_DfDx[i*h_num+jn]=Dt/2*HYPER1[jn*h_num+i].DgDq[A_Z];
				part_fx[i]+=Dt/2*(HYPER[jn].stress[A_Z][A_X]*HYPER1[jn*h_num+i].DgDq[A_X]+HYPER[jn].stress[A_Z][A_Y]*HYPER1[jn*h_num+i].DgDq[A_Y]+HYPER[jn].stress[A_Z][A_Z]*HYPER1[jn*h_num+i].DgDq[A_Z]);
				DfDx[i*h_num+jn]=0;
			}
			w_fx[i]=HYPER[i].p[A_Z]+part_fx[i];
			number++;
		}
//		w_DfDx[i*h_num+i]=1;
	}
	if(number>0)
	{

		output_newton_data3(w_fx,w_DfDx,n_rx,n_ry,n_rz,part_fx,h_num,count,t);		
		gauss(w_DfDx,w_fx,h_num);

		for(int i=0;i<number;i++)
		{
			int in=Nw[i];
			fx[in]=w_fx[in];
		}
	}
	delete[]	Nw;
	delete[]	w_DfDx;
	delete[]	w_fx;
	delete[]	part_fx;//*/

	output_newton_data1(fx,DfDx,n_rx,n_ry,n_rz,h_num,count,t);
	



	////�o��
//	if(count%200==0 && count>CON.get_nr()/2)
//	if(t==1||t%CON.get_interval()==0)	if(count%200==0||count==1)	
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
//	if(repetation==0)	cout<<"���̉^���ʁ��ʒu���W�v�Z";
//	else	cout<<"�^���ʌv�Z";

	int h_num=HYPER.size();
	double Dt=CON.get_dt();
	double le=CON.get_distancebp();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	double G=9.8;
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
		if(flag_G==ON)	p_half_p[A_Z]-=G*mi;
		//�S����
		if(flag_vis==ON)
		{
			p_half_p[A_X]+=HYPER[i].vis_force[A_X];
			p_half_p[A_Y]+=HYPER[i].vis_force[A_Y];
			p_half_p[A_Z]+=HYPER[i].vis_force[A_Z];
		}
		//���͍�
		if(flag_FEM==ON&&PART[i].toFEM==ON)
		{
			p_half_p[A_X]+=V*F[A_X][i];//density;
			p_half_p[A_Y]+=V*F[A_Y][i];//density;
			p_half_p[A_Z]+=V*F[A_Z][i];//density;
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
			/*if(model_num==30||model_num==23)
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
			else*/
//			if(PART[i].r[A_Z]>0)	//way3
			{
				PART[i].r[A_X]+=Dt*HYPER[i].half_p[A_X]/mi;
				PART[i].r[A_Y]+=Dt*HYPER[i].half_p[A_Y]/mi;
				PART[i].r[A_Z]+=Dt*HYPER[i].half_p[A_Z]/mi;
			}
/*			else
			{
				PART[i].r[A_X]+=Dt*HYPER[i].half_p[A_X]/mi;	//way3
				PART[i].r[A_Y]+=Dt*HYPER[i].half_p[A_Y]/mi;	//way3
			}*/

		}
		else
		{
			//�^���ʂ̍X�V
			HYPER[i].p[A_X]=HYPER[i].half_p[A_X]+Dt*0.5*p_half_p[A_X];
			HYPER[i].p[A_Y]=HYPER[i].half_p[A_Y]+Dt*0.5*p_half_p[A_Y];
			HYPER[i].p[A_Z]=HYPER[i].half_p[A_Z]+Dt*0.5*p_half_p[A_Z];////
/*			if(PART[i].r[A_Z]<0)
			{
				HYPER[i].p[A_Z]*=-1;
			}//*/

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
	
//	cout<<"----------OK"<<endl;
}

void calc_F(mpsconfig &CON, vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1)
{
//	cout<<"Fi�v�Z";
	////Fi�̍X�V
	int h_num=HYPER.size();
	double V=get_volume(&CON);

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
	int p_num=PART.size();
	double r=CON.get_h_dis();
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

/*			if(PART[k].r[A_Z]<0)
			{
//				double dis=sqrt((PART[k].r[A_X]-PART[i].r[A_X])*(PART[k].r[A_X]-PART[i].r[A_X])+(PART[k].r[A_Y]-PART[i].r[A_Y])*(PART[k].r[A_Y]-PART[i].r[A_Y])+(PART[k].r[A_Z]-PART[i].r[A_Z])*(PART[k].r[A_Z]-PART[i].r[A_Z]));
				double dis=sqrt(HYPER1[i*h_num+k].aiin[A_X]*HYPER1[i*h_num+k].aiin[A_X]+HYPER1[i*h_num+k].aiin[A_Y]*HYPER1[i*h_num+k].aiin[A_Y]+HYPER1[i*h_num+k].aiin[A_Z]*HYPER1[i*h_num+k].aiin[A_Z]);
				double d_wij=-r/dis/dis;

				HYPER1[k*h_num+i].DgDq[A_X]+=V/HYPER[i].pnd*d_wij*(PART[k].r[A_X]-PART[i].r[A_X])/dis;
				HYPER1[k*h_num+i].DgDq[A_Y]+=V/HYPER[i].pnd*d_wij*(PART[k].r[A_Y]-PART[i].r[A_Y])/dis;
				HYPER1[k*h_num+i].DgDq[A_Z]+=V/HYPER[i].pnd*d_wij*(PART[k].r[A_Z]-PART[i].r[A_Z])/dis;

				HYPER1[i*h_num+i].DgDq[A_X]+=-V/HYPER[i].pnd*d_wij*(PART[k].r[A_X]-PART[i].r[A_X])/dis;
				HYPER1[i*h_num+i].DgDq[A_Y]+=-V/HYPER[i].pnd*d_wij*(PART[k].r[A_Y]-PART[i].r[A_Y])/dis;
				HYPER1[i*h_num+i].DgDq[A_Z]+=-V/HYPER[i].pnd*d_wij*(PART[k].r[A_Z]-PART[i].r[A_Z])/dis;
				HYPER1[k*h_num+i].DgDq[A_X]+=V/HYPER[i].pnd*d_wij*HYPER1[i*h_num+k].aiin[A_X]/dis;
				HYPER1[k*h_num+i].DgDq[A_Y]+=V/HYPER[i].pnd*d_wij*HYPER1[i*h_num+k].aiin[A_Y]/dis;
				HYPER1[k*h_num+i].DgDq[A_Z]+=V/HYPER[i].pnd*d_wij*HYPER1[i*h_num+k].aiin[A_Z]/dis;

				HYPER1[i*h_num+i].DgDq[A_X]+=-V/HYPER[i].pnd*d_wij*HYPER1[i*h_num+k].aiin[A_X]/dis;
				HYPER1[i*h_num+i].DgDq[A_Y]+=-V/HYPER[i].pnd*d_wij*HYPER1[i*h_num+k].aiin[A_Y]/dis;
				HYPER1[i*h_num+i].DgDq[A_Z]+=-V/HYPER[i].pnd*d_wij*HYPER1[i*h_num+k].aiin[A_Z]/dis;

			}//*/
		}
		HYPER1[i*h_num+i].DgDq[A_X]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_X][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_X][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_X][2]*HYPER1[i*h_num+i].n0ij[2]);
		HYPER1[i*h_num+i].DgDq[A_Y]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_Y][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_Y][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_Y][2]*HYPER1[i*h_num+i].n0ij[2]);
		HYPER1[i*h_num+i].DgDq[A_Z]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_Z][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_Z][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_Z][2]*HYPER1[i*h_num+i].n0ij[2]);

	/*	double part_dgdq_x=0;
		double part_dgdq_y=0;
		double part_dgdq_z=0;
		for(int k=h_num;k<p_num;k++)
		{
			double dis=sqrt((PART[k].r[A_X]-PART[i].r[A_X])*(PART[k].r[A_X]-PART[i].r[A_X])+(PART[k].r[A_Y]-PART[i].r[A_Y])*(PART[k].r[A_Y]-PART[i].r[A_Y])+(PART[k].r[A_Z]-PART[i].r[A_Z])*(PART[k].r[A_Z]-PART[i].r[A_Z]));
			if(dis<r)
			{
				double d_wij=-r/dis/dis;
				part_dgdq_x+=-V/HYPER[i].pnd*d_wij*(PART[k].r[A_X]-PART[i].r[A_X])/dis;
				part_dgdq_y+=-V/HYPER[i].pnd*d_wij*(PART[k].r[A_Y]-PART[i].r[A_Y])/dis;
				part_dgdq_z+=-V/HYPER[i].pnd*d_wij*(PART[k].r[A_Z]-PART[i].r[A_Z])/dis;
			}
		}
		HYPER1[i*h_num+i].DgDq[A_X]+=part_dgdq_x;
		HYPER1[i*h_num+i].DgDq[A_Y]+=part_dgdq_y;
		HYPER1[i*h_num+i].DgDq[A_Z]+=part_dgdq_z;//*/

		//cout<<"i"<<i<<"j"<<i<<"	"<<HYPER1[i*h_num+i].DgDq[A_X]<<","<<HYPER1[i*h_num+i].DgDq[A_Y]<<","<<HYPER1[i*h_num+i].DgDq[A_Z]<<endl;
	}
//	cout<<"----------OK"<<endl;

}

void calc_stress(mpsconfig &CON,vector<hyperelastic> &HYPER)
{
//	cout<<"���͌v�Z";
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

//	cout<<"----------OK"<<endl;
}

void calc_differential_p(mpsconfig &CON,vector<mpselastic>PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,double **F)
{
//	cout<<"�^���ʔ����l�v�Z";

	int h_num=HYPER.size();
	int flag_vis=CON.get_flag_vis();
	int flag_G=CON.get_flag_G();
	int flag_FEM=CON.get_FEM_flag();
	double Dt=CON.get_dt();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	double density=CON.get_hyper_density();
	double G=9.8;
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
		if(flag_G==ON)	HYPER[i].differential_p[A_Z]-=Dt*0.5*G*mi;
		if(flag_vis==ON)
		{
			HYPER[i].differential_p[A_X]+=Dt*0.5*HYPER[i].vis_force[A_X];
			HYPER[i].differential_p[A_Y]+=Dt*0.5*HYPER[i].vis_force[A_Y];
			HYPER[i].differential_p[A_Z]+=Dt*0.5*HYPER[i].vis_force[A_Z];
		}
		if(flag_FEM==ON && PART[i].toFEM==ON)
		{
			HYPER[i].differential_p[A_X]+=Dt*0.5*F[A_X][i]*V;//density;
			HYPER[i].differential_p[A_Y]+=Dt*0.5*F[A_Y][i]*V;//density;
			HYPER[i].differential_p[A_Z]+=Dt*0.5*F[A_Z][i]*V;//density;
		}
	}
//	cout<<"----------OK"<<endl;

}

void renew_lambda(mpsconfig &CON,vector<mpselastic>PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int t)
{
	clock_t t3=clock();

//	cout<<"Lambda�v�Z";

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
			}//j�Ɋւ���for���̏I���
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
			}//j�Ɋւ���for���̏I���
			N_Right[k]=N_right+HYPER[k].differential_p[0]*HYPER1[k*h_num+k].DgDq[0]+HYPER[k].differential_p[1]*HYPER1[k*h_num+k].DgDq[1]-HYPER[k].differential_p[2]*HYPER1[k*h_num+k].DgDq[2];//1/mk*N_right;	
			N_Left[k*h_num+k]+=Dt*0.5*(HYPER1[k*h_num+k].DgDq[0]*HYPER1[k*h_num+k].DgDq[0]+HYPER1[k*h_num+k].DgDq[1]*HYPER1[k*h_num+k].DgDq[1]-HYPER1[k*h_num+k].DgDq[2]*HYPER1[k*h_num+k].DgDq[2]);
		}//*/
	}//i�Ɋւ���for���̏I���


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


//	cout<<"----------OK"<<endl;
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
	int n=3;	//�o�͎���
	double mi=get_volume(&CON)*CON.get_hyper_density();

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
	
	if(n==3)
	{
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
	}
	else if(n==2)
	{
		int num=0;
		int *Nxz=new int [h_num];
		fp<<h_num<<" "<<h_num<<endl;	//�ߓ_���Ɨv�f���o��
	
		//�ߓ_�ԍ��Ƃ��̍��W�̏o�� 
		for(int i=0;i<h_num;i++)
		{
			if(PART[i].q0[A_Y]<1/2*le&&PART[i].q0[A_Y]>-1/2*le)
			{
				fp<<i<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
				Nxz[num]=i;
				num++;
			}
		}
		//�v�f�ԍ��Ɨv�f�`��̎�ށA�����ėv�f���\������ߓ_�ԍ��o��
		for(int i=0;i<num;i++)	fp<<i<<"  0 pt "<<i<<endl;

		//fp<<"2 3"<<endl;//�ߓ_�̏��ʂ�2�ŁA�v�f�̏��ʂ�3�Ƃ������ƁB
//		fp<<"10 0"<<endl;//�ߓ_�̏��ʂ�8�ŁA�v�f�̏��ʂ�0�Ƃ������ƁB
//		fp<<"10 1 1 1 1 1 1 1 1 1 1"<<endl;	//���̍s�̏ڍׂ̓w���v���Q��
		//fp<<"8 1 1 1 1 1 1 1 1"<<endl;	//���̍s�̏ڍׂ̓w���v���Q��
		fp<<"7 0"<<endl;
		fp<<"7 1 1 1 1 1 1 1"<<endl;
		fp<<"v_x,"<<endl;
		fp<<"v_y,"<<endl;
		fp<<"v_z,"<<endl;
		fp<<"lambda,"<<endl;
		fp<<"F_x,"<<endl;
		fp<<"F_y,"<<endl;
		fp<<"F_z,"<<endl;
/*		fp<<"v_x,"<<endl;
		fp<<"v_y,"<<endl;
		fp<<"v_z,"<<endl;*/
		//fp<<"P,N/m^2"<<endl;
		//fp<<"value1,??"<<endl;

		//�e�ߓ_�̏��l����
		for(int i=0;i<num;i++)
		{
			int j=Nxz[i];
			fp<<j<<" "<<1/mi*HYPER[j].p[A_X]<<" "<<1/mi*HYPER[j].p[A_Y]<<" "<<1/mi*HYPER[j].p[A_Z]<<" "<<HYPER[j].lambda<<" "<<F[A_X][j]<<" "<<F[A_Y][j]<<" "<<F[A_Z][j]<<endl;
			//fp<<j<<" "<<HYPER[j].p[A_X]<<" "<<HYPER[j].p[A_Y]<<" "<<HYPER[j].p[A_Z]<<" "<<HYPER[j].lambda<<" "<<F[A_X][j]<<" "<<F[A_Y][j]<<" "<<F[A_Z][j]<<" "<<HYPER[j].vjs_force[A_X]<<" "<<HYPER[j].vjs_force[A_Y]<<" "<<HYPER[j].vjs_force[A_Z]<<endl;
			//fp<<i<<" "<<NODE[i].depth<<" "<<NODE[i].L/le<<" "<<NODE[i].potential<<" "<<NODE[i].Fs<<" "<<endl;
		}
		delete[] Nxz;
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


	ofstream j("J.csv", ios::app);
	ofstream lam("lambda.csv", ios::app);
//	ofstream p("P.csv", ios::app);
	ofstream p_an("P_ave_norm.csv", ios::app);
//	ofstream h("model_height.csv", ios::app);

//	if(t==1)	p<<"h_p,,,d_p,,,p\n";
//	p<<"t"<<t<<endl;
	double sum_lam=0;
	double sum_j=0;
	double p_sum=0,d_p_sum=0,h_p_sum=0;
	double h_min=PART[0].r[A_Z];
	double h_max=PART[0].r[A_Z];
	for(int i=0;i<h_num;i++)
	{
		lam<<HYPER[i].lambda<<",";
		j<<HYPER[i].J<<",";
		sum_lam+=HYPER[i].lambda;
		sum_j+=HYPER[i].J;
//		p<<HYPER[i].half_p[A_X]<<","<<HYPER[i].half_p[A_Y]<<","<<HYPER[i].half_p[A_Z]<<","<<HYPER[i].differential_p[A_X]<<","<<HYPER[i].differential_p[A_Y]<<","<<HYPER[i].differential_p[A_Z]<<","<<HYPER[i].p[A_X]<<","<<HYPER[i].p[A_Y]<<","<<HYPER[i].p[A_Z]<<endl;
		p_sum+=sqrt(HYPER[i].p[A_X]*HYPER[i].p[A_X]+HYPER[i].p[A_Y]*HYPER[i].p[A_Y]+HYPER[i].p[A_Z]*HYPER[i].p[A_Z]);
		d_p_sum+=sqrt(HYPER[i].differential_p[A_X]*HYPER[i].differential_p[A_X]+HYPER[i].differential_p[A_Y]*HYPER[i].differential_p[A_Y]+HYPER[i].differential_p[A_Z]*HYPER[i].differential_p[A_Z]);
		h_p_sum+=sqrt(HYPER[i].half_p[A_X]*HYPER[i].half_p[A_X]+HYPER[i].half_p[A_Y]*HYPER[i].half_p[A_Y]+HYPER[i].half_p[A_Z]*HYPER[i].half_p[A_Z]);
//		if(h_min>PART[i].r[A_Z])	h_min=PART[i].r[A_Z];
//		if(h_max<PART[i].r[A_Z])	h_max=PART[i].r[A_Z];
	}
	
	lam<<sum_lam/h_num<<endl;
	j<<sum_j/h_num<<endl;
//	p<<endl;
	if(t==1)
	{
		p_an<<"p_sum,d_p_sum,h_p_sum\n";
//		h<<"h_min,h_max\n";
	}
	p_an<<p_sum/h_num<<","<<d_p_sum/h_num<<","<<h_p_sum/h_num<<endl;
//	h<<h_min<<","<<h_max<<endl;

	lam.close();
	j.close();
//	p.close();
	p_an.close();
//	h.close();

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

//	stringstream ss_dgdq;
//	ss_dgdq<<"./DgDq/DgDq"<<t<<".csv";
//	ofstream dg(ss_dgdq.str());
/*	ofstream d_p("d_P.csv", ios::app);
	ofstream h_p("h_P.csv", ios::app);*/
	stringstream ss_stress;
	ss_stress<<"./Stress/stress"<<t<<".csv";
	ofstream fs(ss_stress.str());
	for(int i=0;i<h_num;i++)
	{
		int N=HYPER[i].N;
//		dg<<i;
		for(int j=0;j<N;j++)
		{
			int nei=HYPER[i].NEI[j];
	//		dg<<","<<nei<<","<<HYPER1[i*h_num+nei].DgDq[A_X]<<","<<HYPER1[i*h_num+nei].DgDq[A_Y]<<","<<HYPER1[i*h_num+nei].DgDq[A_Z]<<","<<HYPER[nei].p[A_X]<<","<<HYPER[nei].p[A_Y]<<","<<HYPER[nei].p[A_Z]<<endl;
		}
		fs<<i<<","<<HYPER[i].stress[A_X][A_X]<<","<<HYPER[i].stress[A_X][A_Y]<<","<<HYPER[i].stress[A_X][A_Z]<<endl;
		fs<<","<<HYPER[i].stress[A_Y][A_X]<<","<<HYPER[i].stress[A_Y][A_Y]<<","<<HYPER[i].stress[A_Y][A_Z]<<endl;
		fs<<","<<HYPER[i].stress[A_Z][A_X]<<","<<HYPER[i].stress[A_Z][A_Y]<<","<<HYPER[i].stress[A_Z][A_Z]<<endl;
	}
//	dg.close();
	fs.close();

	ofstream r("r_CG.csv", ios::app);
	double r_cg_x=0;
	double r_cg_y=0;
	double r_cg_z=0;

	for(int i=0;i<h_num;i++)
	{
		r_cg_x+=PART[i].r[A_X];
		r_cg_y+=PART[i].r[A_Y];
		r_cg_z+=PART[i].r[A_Z];
	}

	r<<r_cg_x/h_num<<","<<r_cg_y/h_num<<","<<r_cg_z/h_num<<endl;
	r.close();

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

void output_newton_data1(double *fx, double *DfDx, double *n_rx, double *n_ry, double *n_rz, int hyper_number,int count, int t)
{
	int h_num=hyper_number;

	stringstream ss_r;
	ss_r<<"./Newton_raphson/position"<<t<<".csv";
/*	stringstream ss_Df;
	ss_Df<<"./Newton_raphson/DfDx "<<t<<".csv";*/
	stringstream ss_fx;
	ss_fx<<"./Newton_raphson/ave_fx"<<t<<".csv";
		
	if(count==1)
	{
		ofstream init0(ss_r.str(), ios::trunc);
		//ofstream init1(ss_Df.str(), ios::trunc);
		ofstream init3(ss_fx.str(), ios::trunc);
	
		init0.close();
		//init1.close();
		init3.close();
	}

	ofstream r(ss_r.str(), ios::app);
	ofstream sfx(ss_fx.str(), ios::app);

	//ofstream Df(ss_Df.str(), ios::app);

	r<<t<<endl;
	//Df<<t<<endl;
	double sum_fx=0;
	for(int i=0;i<h_num;i++)
	{
		r<<n_rx[i]<<","<<n_ry[i]<<","<<n_rz[i]<<endl;
		sum_fx+=fx[i];
		//for(int j=0;j<h_num;j++) Df<<DfDx[i*h_num+j]<<",";
		//Df<<endl;
	}
	sfx<<count<<","<<sum_fx/h_num<<endl;

	r.close();
	sfx.close();
	//Df.close();
}


void output_newton_data2(double E, double *XX, int hyper_number, int count, int t)
{
	int h_num=hyper_number;
	stringstream ss_E;
	ss_E<<"./Newton_raphson/E"<<t<<".csv";
	
	stringstream ss_lam;
	ss_lam<<"./Newton_raphson/lambda"<<t<<".csv";
		
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

	double sum_lam=0;
	lam<<count<<",";
	for(int i=0;i<h_num;i++)
	{
		lam<<XX[i]<<",";
		sum_lam+=XX[i];
	}
	lam<<sum_lam/h_num<<endl;

	e.close();
	lam.close();
}
void output_newton_data3(double *w_fx, double *w_DfDx, double *n_rx, double *n_ry, double *n_rz,double *part_fx, int hyper_number,int count, int t)
{
	int h_num=hyper_number;

	stringstream ss_r;
	ss_r<<"./Newton_raphson/position"<<t<<".csv";
	stringstream ss_Df;
	ss_Df<<"./Newton_raphson/w_DfDx "<<t<<".csv";
	stringstream ss_fx;
	ss_fx<<"./Newton_raphson/w_fx"<<t<<".csv";
		
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

	r<<t<<endl;
	sfx<<"fx"<<","<<"part_fx"<<endl;
	Df<<t<<endl;
	for(int i=0;i<h_num;i++)
	{
		r<<n_rx[i]<<","<<n_ry[i]<<","<<n_rz[i]<<endl;
		sfx<<w_fx[i]<<","<<part_fx[i]<<endl;
		for(int j=0;j<h_num;j++) Df<<w_DfDx[i*h_num+j]<<",";
		Df<<endl;

	}
	r<<endl;
	sfx<<endl;

	r.close();
	sfx.close();

	Df.close();
}

void output_newton2_data1(vector<hyperelastic>HYPER, double *fx, double *DfDx, int hyper_number, int count, int t)
{
	int h_num=hyper_number;

	stringstream ss_p;
	ss_p<<"./Newton_raphson2/P"<<t<<".csv";
	stringstream ss_fx;
	ss_fx<<"./Newton_raphson2/ave_fx"<<t<<".csv";
	
	stringstream ss_DfDx;
	ss_DfDx<<"./Newton_raphson2/Dfdx_t"<<t<<"_count"<<count<<".csv";
	
	if(count==1)
	{
		ofstream init0(ss_p.str(), ios::trunc);
		ofstream init2(ss_fx.str(), ios::trunc);
	
		init0.close();
		init2.close();
	}

	ofstream fp(ss_p.str(), ios::app);
	ofstream ffx(ss_fx.str(), ios::app);
	ofstream fDfDx(ss_DfDx.str(),ios::trunc);

	if(count==1)
	{
		fp<<t;
		ffx<<t;
	}
	fp<<","<<count;	
	int n=0;
	double sum_fx=0;
	ffx<<","<<count<<",";
	for(int i=0;i<h_num;i++)
	{
		if(n==0)	fp<<","<<HYPER[i].p[A_X]<<","<<HYPER[i].p[A_Y]<<","<<HYPER[i].p[A_Z]<<endl;
		else
		{
			fp<<","<<","<<HYPER[i].p[A_X]<<","<<HYPER[i].p[A_Y]<<","<<HYPER[i].p[A_Z]<<endl;
		}
		for(int j=0;j<h_num;j++)	fDfDx<<DfDx[i*h_num+j]<<",";
		fDfDx<<endl;
		ffx<<fx[i]<<",";
		sum_fx+=fx[i];
		n++;
	}
	ffx<<sum_fx/h_num<<","<<h_num<<endl;

	fp.close();
	ffx.close();
	fDfDx.close();
}

void output_energy(mpsconfig CON, vector<mpselastic> PART, vector<hyperelastic> HYPER,int t)
{
//	cout<<"�e���|�e���V�����v�Z";
	int h_num=HYPER.size();
	int p_num=PART.size();
	double d_Fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double c10=CON.get_c10();
	double c01=CON.get_c01();
	vector<double>	W;
	W.reserve(h_num);
	for(int i=0;i<h_num;i++)	W.emplace_back(0);
	double dC[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double dC2[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
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
	}//*/
//	cout<<"----------OK"<<endl;

/*	ofstream e("E.csv", ios::app);
	ofstream e_T("E_T.csv", ios::app);
	ofstream e_g("E_g.csv", ios::app);
	ofstream e_W("E_W.csv", ios::app);
	ofstream e_lam("E_lam.csv", ios::app);//*/
	ofstream e_sum("E_sum.csv", ios::app);

	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
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
		vv=HYPER[i].p[0]*HYPER[i].p[0]+HYPER[i].p[1]*HYPER[i].p[1]+HYPER[i].p[2]*HYPER[i].p[2];
		//energy=0.5/mi*vv+W[i]*V+HYPER[i].lambda*(1-HYPER[i].J)*V;
		if(CON.get_flag_G()==ON)	energy=0.5/mi*vv+mi*G*PART[i].r[A_Z]+W[i]*V+HYPER[i].lambda*(1-HYPER[i].J)*V;
		if(CON.get_flag_G()==OFF)	energy=0.5/mi*vv+W[i]*V+HYPER[i].lambda*(1-HYPER[i].J)*V;
		sum_e_T+=0.5/mi*vv;
		sum_e_g+=mi*G*PART[i].r[A_Z];
		sum_e_lam+=HYPER[i].lambda*(1-HYPER[i].J)*V;
		sum_e_W+=W[i]*V;
		sum_e+=energy;
	/*	e<<energy<<",";
		e_T<<0.5/mi*vv<<",";
		e_g<<mi*G*PART[i].r[A_Z]<<",";
		e_W<<W[i]*V<<",";
		e_lam<<HYPER[i].lambda*(1-HYPER[i].J)*V<<",";//*/
	}
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


void BiCGStab2_method(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X)
{
	ofstream t("time_log.dat",ios::app);
	t<<clock()*CLOCKS_PER_SEC<<"	";

	//val :��[���v�f�̒l
	//ind:��[���v�f�̗�ԍ��i�[�z��
	//ptr:�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
	//X[n]:��
	
	int count=0;
	double pk=1;
	double E=1;//�덷
	double alp,beta,rr,w,ita;

	double *r=new double [pn];	//�c��
	double *P=new double [pn];	//�T���x�N�g��
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

	for(int n=0;n<pn;n++) //�����l
	{
		X[n]=0;
		r[n]=B[n];
	}

	for(int n=0;n<pn;n++)
	{
		P[n]=r[n];//������
		rj[n]=r[n];
		Pj[n]=rj[n];
		AP[n]=0;
		y[n]=0;
		u[n]=0;
		W[n]=0;
		Z[n]=0;
		e[n]=0;
		e_pre[n]=0;//1�X�e�b�v�O��e[]
	}
	double rr0=0;
	for(int n=0;n<pn;n++) rr0+=r[n]*r[n];
	cout<<"rr0="<<rr0<<endl;
	 cout<<"BiCGstab2�@�X�^�[�g  -----���m��="<<pn<<"  ---";
	 unsigned int time=GetTickCount();
	double ep=CON->get_FEMCGep();//��������

	while(E>ep)// EP=CON->get_CGep();//��������(convergence test)
	{
		count++;

		for(int n=0;n<pn;n++)
		{
			//P[n]=r[n]+beta*(P[n]-w*AP[n]);
			P[n]=r[n]+beta*(P[n]-u[n]);
		} 

		////pk(r��rj�̓���)�����߂�
		pk=0;
		for(int n=0;n<pn;n++) pk+=r[n]*rj[n];

		//////////////alp�����߂�
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
			e_pre[n]=e[n];//�ύX�O�̒l���L��
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
		
		//////////////////�덷
		rr=0;
		for(int n=0;n<pn;n++) rr+=r[n]*r[n];
		//E=rr/rr0;
		E=sqrt(rr);
		if(count==1||count%10==0)	cout<<"E="<<E<<" count="<<count<<endl;
		////////////////////////
	}
	
	cout<<"������="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;

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
	//val :�[���v�f�̒l
	//ind:��[���v�f�̗�ԍ��i�[�z��
	//ptr:�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
	//X[n]:��

	double accel=0.87;//CON->get_CGaccl();//�����t�@�N�^
	
	int num2=0;//�Ίp�������܂ށA���O�p�s�񂾂����l���ɂ��ꂽ��[���v�f��
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
	    for(int m=ptr[k];m<ptr[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
			if(ind[m]<=k)
			{
				val2[num2]=val[m];
				ind2[num2]=ind[m];
				if(ind[m]==k) val2[num2]*=accel;//����̧��
				num2++;
			}
		}
	}
	ptr2[pn]=num2;//��������Ă����Ȃ��ƁA�Ō��(int m=ptr2[k];m<ptr2[k+1];m++)�݂����Ȃ��Ƃ��ł��Ȃ�

	int *NUM = new int [pn];
	for(int k=0;k<pn;k++) NUM[k]=0;
	
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
			int J=ind2[m];
			NUM[J]=NUM[J]+1;
		}
	}
	double **VAL=new double *[pn];//�[���v�f�̒l
	for(int i=0;i<pn;i++) VAL[i]=new double [NUM[i]];
	int **IND = new int *[pn];//��[���v�f�̍s�ԍ��i�[�z��
	for(int i=0;i<pn;i++) IND[i]=new int [NUM[i]];
	
	/////////////////////////////////////iccg�@
	double alp,beta;
	double rLDLt_r;
	double E=1;//�덷
	double *AP = new double [pn];
	double *y=new double [pn];
	double *LDLt_r= new double [pn];
	double *D1 = new double [pn];//D�s��
	
	/////�s���S�R���X�L�����
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
	        int i=ind2[m];//��ԍ�
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
					for(int J=ptr2[i];J<ptr2[i+1];J++)//i�s�ڂ̂Ȃ������̈�v������̂�T���Ă���B������Ԃ��H
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
	///�s���S�R���X�L�[��������/////////*/

	///�����ɂ����z��ɒl����
	for(int k=0;k<pn;k++) NUM[k]=0;
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
			int J=ind2[m];
			VAL[J][NUM[J]]=val2[m];
			IND[J][NUM[J]]=k;
			NUM[J]=NUM[J]+1;
		}
	}////////*/

	/////////////////y[i]�����Ƃ߁ALDLt_r[i]�����Ƃ߂�B
	for(int i=0;i<pn;i++)
	{
		if(i==0) y[0]=r[0]/val2[0]; //���i3.77�j 
		else
		{
		    double sum=0;
		    /////////        
		    for(int m=ptr2[i];m<ptr2[i+1]-1;m++) sum+=val2[m]*y[ind2[m]];//���i3.78�j
		    int m=ptr2[i+1]-1;
		    y[i]=(r[i]-sum)/val2[m];
		}
	}////y[i]�����Ƃ܂����B
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

	cout<<"ICCG�@:���m��="<<pn<<" ---";
	unsigned int time=GetTickCount();
	int count=0;
	double ep=CON->get_FEMCGep();//��������
	rLDLt_r=0;
	for(int n=0;n<pn;n++) rLDLt_r+=r[n]*LDLt_r[n];//�ŏ���rLDLt_r���������ŋ��߂�
	while(E>ep)
	{
		count++;
		if(count==pn) cout<<"count=pn"<<endl;
		//////////////alp�����߂�
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
			X[n]+=alp*P[n];// X(k+1)=X(k)+alp*P �X�V��̏ꏊ
			r[n]-=alp*AP[n];// r=r-alp*AP       �X�V��̎c��
			E+=r[n]*r[n];						//�X�V��̌덷
		}
		E=sqrt(E);
		cout<<"E="<<E<<endl;
		////////////////////////
		
		///////////////////////beta
		beta=1.0/rLDLt_r;
		rLDLt_r=0;
		
        /////////////////y[i]�����Ƃ߁ALDLt_r[i]�����Ƃ߂�B
		for(int i=0;i<pn;i++)
		{
			if(i==0) y[0]=r[0]/val2[0]; //���i3.77�j �V
			else
			{
			    double sum=0;
			    for(int m=ptr2[i];m<ptr2[i+1]-1;m++)//�Ίp�����͏�������ptr[i+1]-1
			    {
			        sum+=val2[m]*y[ind2[m]];//���i3.78�j
			    }
			    int m=ptr2[i+1]-1;
			    y[i]=(r[i]-sum)/val2[m];
			}
		}////y[i]�����Ƃ܂����B
	
		/////////LDLt_r[i]�����߂�
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
	cout<<"������="<<count<<" time="<<(GetTickCount()-time)*0.001<<"/";
		
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
	//*count2=count;//�����񐔂��i�[���ĕԂ�
}

void CG3D(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X)
{
	cout<<"CG�@:���m��="<<pn<<" ---";
	unsigned int time=GetTickCount();
	int count=0;
	double ep=CON->get_FEMCGep();//��������

	double alp,beta;
	double E=1;//�덷
	double *AP = new double [pn];
	double *r=new double [pn];
	double *P=new double [pn];

	while(E>ep)
	{
		count++;
		if(count==pn) cout<<"count=pn"<<endl;
		//////////////alp�����߂�
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
			X[n]+=alp*P[n];// X(k+1)=X(k)+alp*P �X�V��̏ꏊ
			r[n]-=alp*AP[n];// r=r-alp*AP       �X�V��̎c��
			E+=r[n]*r[n];						//�X�V��̌덷
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
	cout<<"������="<<count<<" time="<<(GetTickCount()-time)*0.001<<"/";
			
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

void QP(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int t)
{
	int h_num=HYPER.size();

	double **dp=new double *[h_num];
	double **old_dp=new double *[h_num];

	double **dq=new double *[h_num];
	double **old_dq=new double *[h_num];

	double Tr=0;
	double **dpTr=new double *[h_num];
	double **dqTr=new double *[h_num];

	double *J=new double [h_num];
	double **F=new double *[h_num];
	double **ti_F=new double *[h_num];
	double **S=new double *[h_num];


	double V=get_volume(&CON);
	double mi=CON.get_h_dis()*V;
	double Dt=CON.get_dt();

	double *w=new double [h_num];
	double *g=new double [h_num];
	double *h=new double [h_num];

	double *seta_g=new double [h_num];
	double *seta_h=new double [h_num];

	double a[DIMENSION]={0,0,0};
	double n[DIMENSION]={0,0,1};

	int Nx=h_num*6;
	double **d_p=new double *[h_num];	
	double **d_q=new double *[h_num];	
	double *B=new double [Nx*Nx];
	double *Nr=new double [Nx];

	for(int i=0;i<h_num;i++)
	{
		dp[i]=new double [DIMENSION];
		old_dp[i]=new double [DIMENSION];
	
		dq[i]=new double [DIMENSION];
		old_dq[i]=new double [DIMENSION];

		dpTr[i]=new double [DIMENSION];
		dqTr[i]=new double [DIMENSION];

		F[i]=new double [DIMENSION*DIMENSION];
		ti_F[i]=new double [DIMENSION*DIMENSION];

		S[i]=new double [DIMENSION*DIMENSION];
		
		d_p[i]=new double [DIMENSION];
		d_q[i]=new double [DIMENSION];
	}

	//������
	for(int i=0; i<h_num; i++)
	{
		dp[i][A_X]=0;	dp[i][A_Y]=0;	dp[i][A_Z]=0;
		dq[i][A_X]=0;	dq[i][A_Y]=0;	dq[i][A_Z]=0;
		
		old_dp[i][A_X]=0;	old_dp[i][A_Y]=0;	old_dp[i][A_Z]=0;
		old_dq[i][A_X]=0;	old_dq[i][A_Y]=0;	old_dq[i][A_Z]=0;
		
		dpTr[i][A_X]=0;	dpTr[i][A_Y]=0;	dpTr[i][A_Z]=0;
		dqTr[i][A_X]=0;	dqTr[i][A_Y]=0;	dqTr[i][A_Z]=0;

		F[i][A_X*DIMENSION+A_X]=0;	F[i][A_X*DIMENSION+A_Y]=0;	F[i][A_X*DIMENSION+A_Z]=0;	
		F[i][A_Y*DIMENSION+A_X]=0;	F[i][A_Y*DIMENSION+A_Y]=0;	F[i][A_Y*DIMENSION+A_Z]=0;	
		F[i][A_Z*DIMENSION+A_X]=0;	F[i][A_Z*DIMENSION+A_Y]=0;	F[i][A_Z*DIMENSION+A_Z]=0;	

		ti_F[i][A_X*DIMENSION+A_X]=0;	ti_F[i][A_X*DIMENSION+A_Y]=0;	ti_F[i][A_X*DIMENSION+A_Z]=0;	
		ti_F[i][A_Y*DIMENSION+A_X]=0;	ti_F[i][A_Y*DIMENSION+A_Y]=0;	ti_F[i][A_Y*DIMENSION+A_Z]=0;	
		ti_F[i][A_Z*DIMENSION+A_X]=0;	ti_F[i][A_Z*DIMENSION+A_Y]=0;	ti_F[i][A_Z*DIMENSION+A_Z]=0;	

		S[i][A_X*DIMENSION+A_X]=0;	S[i][A_X*DIMENSION+A_Y]=0;	S[i][A_X*DIMENSION+A_Z]=0;	
		S[i][A_Y*DIMENSION+A_X]=0;	S[i][A_Y*DIMENSION+A_Y]=0;	S[i][A_Y*DIMENSION+A_Z]=0;	
		S[i][A_Z*DIMENSION+A_X]=0;	S[i][A_Z*DIMENSION+A_Y]=0;	S[i][A_Z*DIMENSION+A_Z]=0;	

		d_p[i][A_X]=0;	d_p[i][A_Y]=0;	d_p[i][A_Z]=0;
		d_q[i][A_X]=0;	d_q[i][A_Y]=0;	d_q[i][A_Z]=0;
		
		w[i]=0;
		g[i]=0;
		h[i]=0;
		seta_g[i]=0;
		seta_h[i]=0;
		J[i]=0;
	}
	for(int i=0;i<h_num*6;i++)
	{
		for(int j=0;j<h_num*6;j++)
		{
			if(j==i)	B[i*h_num*6+j]=1;
			else
			{
				B[i*h_num*6+j]=0;
			}			
		}
		Nr[i]=0;
	}




	double r=0;
	double ep=1e-10;

	double E_min=1;
	int count_min=0;

	//�����v�Z�J�n
	while(E_min>ep)
	{
		count_min++;

		for(int i=0;i<h_num;i++)
		{
			old_dp[i][A_X]=dp[i][A_X];	old_dp[i][A_Y]=dp[i][A_Y];	old_dp[i][A_Z]=dp[i][A_Z];
			old_dq[i][A_X]=dq[i][A_X];	old_dq[i][A_Y]=dq[i][A_Y];	old_dq[i][A_Z]=dq[i][A_Z];

			d_p[i][A_X]=0;	d_p[i][A_Y]=0;	d_p[i][A_Z]=0;
			d_q[i][A_X]=0;	d_q[i][A_Y]=0;	d_q[i][A_Z]=0;		
		}

		for(int i=0;i<Nx;i++)
		{
			for(int j=0;j<Nx;j++)
			{
				if(j==i)	B[i*Nx+j]=1;
				else
				{
					B[i*Nx+j]=0;
				}			
			}
			Nr[i]=0;
		}

		double E=1;
		int count=0;

		calc_wg(CON, PART, HYPER, HYPER1, dq, w, g, J, F, ti_F, S);	//�ꏊ���������Ȃ�����

		Tr=0;
		for(int i=0;i<h_num;i++)
		{
			Tr+=0.5/mi*( 2*Dt*( HYPER[i].p[A_X]*dp[i][A_X]+HYPER[i].p[A_Y]*dp[i][A_Y]+HYPER[i].p[A_Z]*dp[i][A_Z] ) + Dt*Dt*( dp[i][A_X]*dp[i][A_X]+dp[i][A_Y]*dp[i][A_Y]+dp[i][A_Z]*dp[i][A_Z] ) )
				-V*(w[i]-HYPER[i].w)
				+0.5*r*(g[i]+seta_g[i])*(g[i]+seta_g[i]);			

			h[i]=-V*( (PART[i].r[A_X]+Dt*dp[i][A_X] - a[A_X])*n[A_X] + (PART[i].r[A_Y]+Dt*dp[i][A_Y] - a[A_Y])*n[A_Y] +(PART[i].r[A_Z]+Dt*dp[i][A_Z] - a[A_Z])*n[A_Z] );	
			if(h[i]+seta_h[i]>0)	Tr+=0.5*r*(h[i]+seta_h[i])*(h[i]+seta_h[i]);
		}

		for(int i=0;i<h_num;i++)
		{
			dpTr[i][A_X]=Dt/mi*(HYPER[i].p[A_X]+dp[i][A_X]);
			dpTr[i][A_Y]=Dt/mi*(HYPER[i].p[A_X]+dp[i][A_Y]);
			dpTr[i][A_Z]=Dt/mi*(HYPER[i].p[A_X]+dp[i][A_Z]);

	
			dqTr[i][A_X]=0;	dqTr[i][A_Y]=0;	dqTr[i][A_Z]=0;
			for(int j=0;j<h_num;j++)
			{
				//w�̕Δ�������
				dqTr[i][A_X]=-Dt*( ( F[j][A_X*DIMENSION+A_X]*S[j][A_X*DIMENSION+A_X] + F[j][A_X*DIMENSION+A_Y]*S[j][A_Y*DIMENSION+A_X] + F[j][A_X*DIMENSION+A_Z]*S[j][A_Z*DIMENSION+A_X] )*HYPER1[i*h_num+j].n0ij[A_X]
				+( F[j][A_X*DIMENSION+A_X]*S[j][A_X*DIMENSION+A_Y] + F[j][A_X*DIMENSION+A_Y]*S[j][A_Y*DIMENSION+A_Y] + F[j][A_X*DIMENSION+A_Z]*S[j][A_Z*DIMENSION+A_Y] )*HYPER1[i*h_num+j].n0ij[A_Y]
				+( F[j][A_X*DIMENSION+A_X]*S[j][A_X*DIMENSION+A_Z] + F[j][A_X*DIMENSION+A_Y]*S[j][A_Y*DIMENSION+A_Z] + F[j][A_X*DIMENSION+A_Z]*S[j][A_Z*DIMENSION+A_Z] )*HYPER1[i*h_num+j].n0ij[A_Z] );

				dqTr[i][A_Y]=-Dt*( ( F[j][A_Y*DIMENSION+A_X]*S[j][A_X*DIMENSION+A_X] + F[j][A_Y*DIMENSION+A_Y]*S[j][A_Y*DIMENSION+A_X] + F[j][A_Y*DIMENSION+A_Z]*S[j][A_Z*DIMENSION+A_X] )*HYPER1[i*h_num+j].n0ij[A_X]
				+( F[j][A_Y*DIMENSION+A_X]*S[j][A_X*DIMENSION+A_Y] + F[j][A_Y*DIMENSION+A_Y]*S[j][A_Y*DIMENSION+A_Y] + F[j][A_Y*DIMENSION+A_Z]*S[j][A_Z*DIMENSION+A_Y] )*HYPER1[i*h_num+j].n0ij[A_Y]
				+( F[j][A_Y*DIMENSION+A_X]*S[j][A_X*DIMENSION+A_Z] + F[j][A_Y*DIMENSION+A_Y]*S[j][A_Y*DIMENSION+A_Z] + F[j][A_Y*DIMENSION+A_Z]*S[j][A_Z*DIMENSION+A_Z] )*HYPER1[i*h_num+j].n0ij[A_Z] );

				dqTr[i][A_Z]=-Dt*( ( F[j][A_Z*DIMENSION+A_X]*S[j][A_X*DIMENSION+A_X] + F[j][A_Z*DIMENSION+A_Y]*S[j][A_Y*DIMENSION+A_X] + F[j][A_Z*DIMENSION+A_Z]*S[j][A_Z*DIMENSION+A_X] )*HYPER1[i*h_num+j].n0ij[A_X]
				+( F[j][A_Z*DIMENSION+A_X]*S[j][A_X*DIMENSION+A_Y] + F[j][A_Z*DIMENSION+A_Y]*S[j][A_Y*DIMENSION+A_Y] + F[j][A_Z*DIMENSION+A_Z]*S[j][A_Z*DIMENSION+A_Y] )*HYPER1[i*h_num+j].n0ij[A_Y]
				+( F[j][A_Z*DIMENSION+A_X]*S[j][A_X*DIMENSION+A_Z] + F[j][A_Z*DIMENSION+A_Y]*S[j][A_Y*DIMENSION+A_Z] + F[j][A_Z*DIMENSION+A_Z]*S[j][A_Z*DIMENSION+A_Z] )*HYPER1[i*h_num+j].n0ij[A_Z] );

				//���̕Δ�������
				dqTr[i][A_X]=-r*Dt*J[j]*( ti_F[j][A_X*DIMENSION+A_X]*HYPER1[i*h_num+j].n0ij[A_X] + ti_F[j][A_X*DIMENSION+A_Y]*HYPER1[i*h_num+j].n0ij[A_Y] + ti_F[j][A_X*DIMENSION+A_Z]*HYPER1[i*h_num+j].n0ij[A_Z] )*(g[j]+seta_g[j]);
				dqTr[i][A_Y]=-r*Dt*J[j]*( ti_F[j][A_Y*DIMENSION+A_X]*HYPER1[i*h_num+j].n0ij[A_X] + ti_F[j][A_Y*DIMENSION+A_Y]*HYPER1[i*h_num+j].n0ij[A_Y] + ti_F[j][A_Y*DIMENSION+A_Z]*HYPER1[i*h_num+j].n0ij[A_Z] )*(g[j]+seta_g[j]);
				dqTr[i][A_Z]=-r*Dt*J[j]*( ti_F[j][A_Z*DIMENSION+A_X]*HYPER1[i*h_num+j].n0ij[A_X] + ti_F[j][A_Z*DIMENSION+A_Y]*HYPER1[i*h_num+j].n0ij[A_Y] + ti_F[j][A_Z*DIMENSION+A_Z]*HYPER1[i*h_num+j].n0ij[A_Z] )*(g[j]+seta_g[j]);
			}
			if(h[i]+seta_h[i]>0)
			{
				dqTr[i][A_X]=-r*V*Dt*n[A_X]*(h[i]+seta_h[i]);
				dqTr[i][A_Y]=-r*V*Dt*n[A_Y]*(h[i]+seta_h[i]);
				dqTr[i][A_Z]=-r*V*Dt*n[A_Z]*(h[i]+seta_h[i]);
			}
		
			E+=sqrt(dpTr[i][A_X]*dpTr[i][A_X] + dpTr[i][A_Y]*dpTr[i][A_Y] + dpTr[i][A_Z]*dpTr[i][A_Z] + dqTr[i][A_X]*dqTr[i][A_X] + dqTr[i][A_Y]*dqTr[i][A_Y] + dqTr[i][A_Z]*dqTr[i][A_Z]);
		}
		cout<<"E0="<<E<<endl;
		if(E<ep)	break;
		else
		{
			while(E>ep)
			{
				count++;

				double **dp_k=new double *[h_num];
				double **dq_k=new double *[h_num];
				double **dpTr_k=new double *[h_num];
				double **dqTr_k=new double *[h_num];

				double **dp_a=new double *[h_num];
				double **dq_a=new double *[h_num];
	
				for(int i=0;i<h_num;i++)
				{
					dp_k[i]=new double [DIMENSION];
					dq_k[i]=new double [DIMENSION];
					dpTr_k[i]=new double [DIMENSION];
					dqTr_k[i]=new double [DIMENSION];
					dp_a[i]=new double [DIMENSION];
					dq_a[i]=new double [DIMENSION];
				}

				double *w_a=new double [h_num];
				double *g_a=new double [h_num];
				double *h_a=new double [h_num];
			
				double Tr_min=Tr;
				double a_min=1e-3;
				for(int i=0;i<1000;i++)
				{
					double alpha=(i+1)*1e-3;

					double Tr_a=0;
					for(int i=0;i<h_num;i++)
					{
						dp_a[i][A_X]=dp[i][A_X]+d_p[i][A_X]*alpha;
						dp_a[i][A_Y]=dp[i][A_Y]+d_p[i][A_Y]*alpha;
						dp_a[i][A_Z]=dp[i][A_Z]+d_p[i][A_Z]*alpha;
						dq_a[i][A_X]=dq[i][A_X]+d_q[i][A_X]*alpha;
						dq_a[i][A_Y]=dq[i][A_Y]+d_q[i][A_Y]*alpha;
						dq_a[i][A_Z]=dq[i][A_Z]+d_q[i][A_Z]*alpha;
						w_a[i]=0;
						g_a[i]=0;
					}

					calc_wg(CON, PART, HYPER, HYPER1, dq_a, w_a, g_a, J, F, ti_F, S);�@	//�ꏊ���������Ȃ�����
					
					for(int i=0;i<h_num;i++)
					{
						Tr_a+=0.5/mi*( 2*Dt*( HYPER[i].p[A_X]*dp_a[i][A_X]+HYPER[i].p[A_Y]*dp_a[i][A_Y]+HYPER[i].p[A_Z]*dp_a[i][A_Z] ) + Dt*Dt*( dp_a[i][A_X]*dp_a[i][A_X]+dp_a[i][A_Y]*dp_a[i][A_Y]+dp_a[i][A_Z]*dp_a[i][A_Z] ) )
							-V*(w_a[i]-HYPER[i].w)	+0.5*r*(g_a[i]+seta_g[i])*(g_a[i]+seta_g[i]);			

						h_a[i]=-V*( (PART[i].r[A_X]+Dt*dp_a[i][A_X] - a[A_X])*n[A_X] + (PART[i].r[A_Y]+Dt*dp_a[i][A_Y] - a[A_Y])*n[A_Y] +(PART[i].r[A_Z]+Dt*dp_a[i][A_Z] - a[A_Z])*n[A_Z] );	
						if(h_a[i]+seta_h[i]>0)	Tr+=0.5*r*(h_a[i]+seta_h[i])*(h_a[i]+seta_h[i]);
					}

					if(Tr_a<Tr_min)
					{
						Tr_min=Tr_a;
						a_min=alpha;
					}
				}
				cout<<"Tr"<<count<<"="<<Tr_min<<", alpha="<<a_min<<endl;
				
				for(int i=0;i<Nx;i++)
				{
					dp[i][A_X]+=d_q[i][A_X]*a_min;
					dp[i][A_Y]+=d_q[i][A_Y]*a_min;
					dp[i][A_Z]+=d_q[i][A_Z]*a_min;
					dq[i][A_X]+=d_q[i][A_X]*a_min;
					dq[i][A_Y]+=d_q[i][A_Y]*a_min;
					dq[i][A_Z]+=d_q[i][A_Z]*a_min;
				}


				calc_wg(CON, PART, HYPER, HYPER1, dq, w, g, J, F, ti_F, S);	//�ꏊ���������Ȃ�����

				Tr=0;
				for(int i=0;i<h_num;i++)
				{
					Tr+=0.5/mi*( 2*Dt*( HYPER[i].p[A_X]*dp[i][A_X]+HYPER[i].p[A_Y]*dp[i][A_Y]+HYPER[i].p[A_Z]*dp[i][A_Z] ) + Dt*Dt*( dp[i][A_X]*dp[i][A_X]+dp[i][A_Y]*dp[i][A_Y]+dp[i][A_Z]*dp[i][A_Z] ) )
						-V*(w[i]-HYPER[i].w)
						+0.5*r*(g[i]+seta_g[i])*(g[i]+seta_g[i]);			

					h[i]=-V*( (PART[i].r[A_X]+Dt*dp[i][A_X] - a[A_X])*n[A_X] + (PART[i].r[A_Y]+Dt*dp[i][A_Y] - a[A_Y])*n[A_Y] +(PART[i].r[A_Z]+Dt*dp[i][A_Z] - a[A_Z])*n[A_Z] );	
					if(h[i]+seta_h[i]>0)	Tr+=0.5*r*(h[i]+seta_h[i])*(h[i]+seta_h[i]);
				}

				for(int i=0;i<h_num;i++)
				{
					dpTr[i][A_X]=Dt/mi*(HYPER[i].p[A_X]+dp[i][A_X]);
					dpTr[i][A_Y]=Dt/mi*(HYPER[i].p[A_X]+dp[i][A_Y]);
					dpTr[i][A_Z]=Dt/mi*(HYPER[i].p[A_X]+dp[i][A_Z]);

	
					dqTr[i][A_X]=0;	dqTr[i][A_Y]=0;	dqTr[i][A_Z]=0;
					for(int j=0;j<h_num;j++)
					{
						//w�̕Δ�������
						dqTr[i][A_X]=-Dt*( ( F[j][A_X*DIMENSION+A_X]*S[j][A_X*DIMENSION+A_X] + F[j][A_X*DIMENSION+A_Y]*S[j][A_Y*DIMENSION+A_X] + F[j][A_X*DIMENSION+A_Z]*S[j][A_Z*DIMENSION+A_X] )*HYPER1[i*h_num+j].n0ij[A_X]
						+( F[j][A_X*DIMENSION+A_X]*S[j][A_X*DIMENSION+A_Y] + F[j][A_X*DIMENSION+A_Y]*S[j][A_Y*DIMENSION+A_Y] + F[j][A_X*DIMENSION+A_Z]*S[j][A_Z*DIMENSION+A_Y] )*HYPER1[i*h_num+j].n0ij[A_Y]
						+( F[j][A_X*DIMENSION+A_X]*S[j][A_X*DIMENSION+A_Z] + F[j][A_X*DIMENSION+A_Y]*S[j][A_Y*DIMENSION+A_Z] + F[j][A_X*DIMENSION+A_Z]*S[j][A_Z*DIMENSION+A_Z] )*HYPER1[i*h_num+j].n0ij[A_Z] );

						dqTr[i][A_Y]=-Dt*( ( F[j][A_Y*DIMENSION+A_X]*S[j][A_X*DIMENSION+A_X] + F[j][A_Y*DIMENSION+A_Y]*S[j][A_Y*DIMENSION+A_X] + F[j][A_Y*DIMENSION+A_Z]*S[j][A_Z*DIMENSION+A_X] )*HYPER1[i*h_num+j].n0ij[A_X]
						+( F[j][A_Y*DIMENSION+A_X]*S[j][A_X*DIMENSION+A_Y] + F[j][A_Y*DIMENSION+A_Y]*S[j][A_Y*DIMENSION+A_Y] + F[j][A_Y*DIMENSION+A_Z]*S[j][A_Z*DIMENSION+A_Y] )*HYPER1[i*h_num+j].n0ij[A_Y]
						+( F[j][A_Y*DIMENSION+A_X]*S[j][A_X*DIMENSION+A_Z] + F[j][A_Y*DIMENSION+A_Y]*S[j][A_Y*DIMENSION+A_Z] + F[j][A_Y*DIMENSION+A_Z]*S[j][A_Z*DIMENSION+A_Z] )*HYPER1[i*h_num+j].n0ij[A_Z] );

						dqTr[i][A_Z]=-Dt*( ( F[j][A_Z*DIMENSION+A_X]*S[j][A_X*DIMENSION+A_X] + F[j][A_Z*DIMENSION+A_Y]*S[j][A_Y*DIMENSION+A_X] + F[j][A_Z*DIMENSION+A_Z]*S[j][A_Z*DIMENSION+A_X] )*HYPER1[i*h_num+j].n0ij[A_X]
						+( F[j][A_Z*DIMENSION+A_X]*S[j][A_X*DIMENSION+A_Y] + F[j][A_Z*DIMENSION+A_Y]*S[j][A_Y*DIMENSION+A_Y] + F[j][A_Z*DIMENSION+A_Z]*S[j][A_Z*DIMENSION+A_Y] )*HYPER1[i*h_num+j].n0ij[A_Y]
						+( F[j][A_Z*DIMENSION+A_X]*S[j][A_X*DIMENSION+A_Z] + F[j][A_Z*DIMENSION+A_Y]*S[j][A_Y*DIMENSION+A_Z] + F[j][A_Z*DIMENSION+A_Z]*S[j][A_Z*DIMENSION+A_Z] )*HYPER1[i*h_num+j].n0ij[A_Z] );

						//���̕Δ�������
						dqTr[i][A_X]=-r*Dt*J[j]*( ti_F[j][A_X*DIMENSION+A_X]*HYPER1[i*h_num+j].n0ij[A_X] + ti_F[j][A_X*DIMENSION+A_Y]*HYPER1[i*h_num+j].n0ij[A_Y] + ti_F[j][A_X*DIMENSION+A_Z]*HYPER1[i*h_num+j].n0ij[A_Z] )*(g[j]+seta_g[j]);
						dqTr[i][A_Y]=-r*Dt*J[j]*( ti_F[j][A_Y*DIMENSION+A_X]*HYPER1[i*h_num+j].n0ij[A_X] + ti_F[j][A_Y*DIMENSION+A_Y]*HYPER1[i*h_num+j].n0ij[A_Y] + ti_F[j][A_Y*DIMENSION+A_Z]*HYPER1[i*h_num+j].n0ij[A_Z] )*(g[j]+seta_g[j]);
						dqTr[i][A_Z]=-r*Dt*J[j]*( ti_F[j][A_Z*DIMENSION+A_X]*HYPER1[i*h_num+j].n0ij[A_X] + ti_F[j][A_Z*DIMENSION+A_Y]*HYPER1[i*h_num+j].n0ij[A_Y] + ti_F[j][A_Z*DIMENSION+A_Z]*HYPER1[i*h_num+j].n0ij[A_Z] )*(g[j]+seta_g[j]);
					}
					if(h[i]+seta_h[i]>0)
					{
						dqTr[i][A_X]=-r*V*Dt*n[A_X]*(h[i]+seta_h[i]);
						dqTr[i][A_Y]=-r*V*Dt*n[A_Y]*(h[i]+seta_h[i]);
						dqTr[i][A_Z]=-r*V*Dt*n[A_Z]*(h[i]+seta_h[i]);
					}
		
					E+=sqrt(dpTr[i][A_X]*dpTr[i][A_X] + dpTr[i][A_Y]*dpTr[i][A_Y] + dpTr[i][A_Z]*dpTr[i][A_Z] + dqTr[i][A_X]*dqTr[i][A_X] + dqTr[i][A_Y]*dqTr[i][A_Y] + dqTr[i][A_Z]*dqTr[i][A_Z]);
				}
				if(count%100==0)	cout<<"E"<<count<<"="<<E<<endl;
				if(E<ep)	break;

				double *s=new double [Nx];
				double *y=new double [Nx];
				double *sB=new double [Nx];
				double *Bs=new double [Nx];
				double *BssB=new double [Nx*Nx];

				for(int i=0;i<h_num;i++)
				{
					s[i]=0;	s[i+h_num]=0; s[i+2*h_num]=0; s[i+3*h_num]=0; s[i+4*h_num]=0; s[i+5*h_num]=0;
					y[i]=0;	y[i+h_num]=0; y[i+2*h_num]=0;	y[i+3*h_num]=0; y[i+4*h_num]=0; y[i+5*h_num]=0;
					sB[i]=0; sB[i+h_num]=0; sB[i+2*h_num]=0; sB[i+3*h_num]=0; sB[i+4*h_num]=0; sB[i+5*h_num]=0;
					Bs[i]=0; Bs[i+h_num]=0; Bs[i+2*h_num]=0; Bs[i+3*h_num]=0; Bs[i+4*h_num]=0; Bs[i+5*h_num]=0;

					for(int j=0;j<Nx;j++)
					{
						BssB[i*Nx+j]=0;
						BssB[(i+h_num)*Nx+j+h_num]=0;
						BssB[(i+2*h_num)*Nx+j+2*h_num]=0;
						BssB[(i+3*h_num)*Nx+j+3*h_num]=0;
						BssB[(i+4*h_num)*Nx+j+4*h_num]=0;
						BssB[(i+5*h_num)*Nx+j+5*h_num]=0;
					}
				}



				double beta=0;
				double sigma=0;

				for(int i=0;i<Nx;i++)
				{
					s[i]=dp[i][A_X]-dp_k[i][A_X];			s[i+h_num]=dp[i][A_Y]-dp_k[i][A_Y];		s[i+2*h_num]=dp[i][A_Z]-dp_k[i][A_Z];
					s[i+3*h_num]=dq[i][A_X]-dq_k[i][A_X];	s[i+4*h_num]=dq[i][A_Y]-dq_k[i][A_Y];	s[i+5*h_num]=dq[i][A_Z]-dq_k[i][A_Z];
	
					y[i]=dpTr[i][A_X]-dpTr_k[i][A_X];			y[i+h_num]=dpTr[i][A_Y]-dpTr_k[i][A_Y];			y[i+2*h_num]=dpTr[i][A_Z]-dpTr_k[i][A_Z];
					y[i+3*h_num]=dqTr[i][A_X]-dqTr_k[i][A_X];	y[i+4*h_num]=dqTr[i][A_Y]-dqTr_k[i][A_Y];		y[i+5*h_num]=dqTr[i][A_Z]-dqTr_k[i][A_Z];
					
					beta+=s[i]*y[i]+s[i+h_num]*y[i+h_num]+s[i+2*h_num]*y[i+2*h_num]+s[i+3*h_num]*y[i+3*h_num]+s[i+4*h_num]*y[i+4*h_num]+s[i+5*h_num]*y[i+5*h_num];
	
					sB[i]=0; sB[i+h_num]=0; sB[i+2*h_num]=0; sB[i+3*h_num]=0; sB[i+4*h_num]=0; sB[i+5*h_num]=0;
					for(int j=0;j<h_num;j++)
					{
						sB[i]+=s[j]*B[j*Nx+i]+s[(j+h_num)]*B[(j+h_num)*Nx+i]+s[(j+2*h_num)]*B[(j+2*h_num)*Nx+i]+s[(j+3*h_num)]*B[(j+3*h_num)*Nx+i]+s[(j+4*h_num)]*B[(j+4*h_num)*Nx+i]+s[(j+5*h_num)]*B[(j+5*h_num)*Nx+i];
						sB[i+h_num]+=s[j]*B[j*Nx+i+h_num]+s[(j+h_num)]*B[(j+h_num)*Nx+i+h_num]+s[(j+2*h_num)]*B[(j+2*h_num)*Nx+i+h_num]+s[(j+3*h_num)]*B[(j+3*h_num)*Nx+i+h_num]+s[(j+4*h_num)]*B[(j+4*h_num)*Nx+i+h_num]+s[(j+5*h_num)]*B[(j+5*h_num)*Nx+i+h_num];
						sB[i+2*h_num]+=s[j]*B[j*Nx+i+2*h_num]+s[(j+h_num)]*B[(j+h_num)*Nx+i+2*h_num]+s[(j+2*h_num)]*B[(j+2*h_num)*Nx+i+2*h_num]+s[(j+3*h_num)]*B[(j+3*h_num)*Nx+i+2*h_num]+s[(j+4*h_num)]*B[(j+4*h_num)*Nx+i+2*h_num]+s[(j+5*h_num)]*B[(j+5*h_num)*Nx+i+2*h_num];
						sB[i+3*h_num]+=s[j]*B[j*Nx+i+3*h_num]+s[(j+h_num)]*B[(j+h_num)*Nx+i+3*h_num]+s[(j+2*h_num)]*B[(j+2*h_num)*Nx+i+3*h_num]+s[(j+3*h_num)]*B[(j+3*h_num)*Nx+i+3*h_num]+s[(j+4*h_num)]*B[(j+4*h_num)*Nx+i+3*h_num]+s[(j+5*h_num)]*B[(j+5*h_num)*Nx+i+3*h_num];
						sB[i+4*h_num]+=s[j]*B[j*Nx+i+4*h_num]+s[(j+h_num)]*B[(j+h_num)*Nx+i+4*h_num]+s[(j+2*h_num)]*B[(j+2*h_num)*Nx+i+4*h_num]+s[(j+3*h_num)]*B[(j+3*h_num)*Nx+i+4*h_num]+s[(j+4*h_num)]*B[(j+4*h_num)*Nx+i+4*h_num]+s[(j+5*h_num)]*B[(j+5*h_num)*Nx+i+4*h_num];
						sB[i+5*h_num]+=s[j]*B[j*Nx+i+5*h_num]+s[(j+h_num)]*B[(j+h_num)*Nx+i+5*h_num]+s[(j+2*h_num)]*B[(j+2*h_num)*Nx+i+5*h_num]+s[(j+3*h_num)]*B[(j+3*h_num)*Nx+i+5*h_num]+s[(j+4*h_num)]*B[(j+4*h_num)*Nx+i+5*h_num]+s[(j+5*h_num)]*B[(j+5*h_num)*Nx+i+5*h_num];
					}
					sigma+=sB[i]*s[i]+sB[i+h_num]*s[i+h_num]+sB[i+2*h_num]*s[i+2*h_num]+sB[i+3*h_num]*s[i+3*h_num]+sB[i+4*h_num]*s[i+4*h_num]+sB[i+5*h_num]*s[i+5*h_num];
				}
				cout<<"beta"<<count<<"="<<beta<<", sigma"<<count<<"="<<sigma<<endl;

				if(beta>0)//(beta>=0.2*sigma)
				{
					for(int i=0;i<h_num;i++)
					{
						Bs[i]=0; Bs[i+h_num]=0; Bs[i+2*h_num]=0; Bs[i+3*h_num]=0; Bs[i+4*h_num]=0; Bs[i+5*h_num]=0;
						for(int j=0;j<h_num;j++)
						{
							Bs[i]+=			B[i*Nx+j]*s[j]			+B[i*Nx+j+h_num]*s[j+h_num]				+B[i*Nx+j+2*h_num]*s[j+2*h_num]				+B[i*Nx+j+3*h_num]*s[j+3*h_num]				+B[i*Nx+j+4*h_num]*s[j+4*h_num]				+B[i*Nx+j+5*h_num]*s[j+5*h_num];
							Bs[i+h_num]+=	B[(i+h_num)*Nx+j]*s[j]	+B[(i+h_num)*Nx+j+h_num]*s[j+h_num]		+B[(i+h_num)*Nx+j+2*h_num]*s[j+2*h_num]		+B[(i+h_num)*Nx+j+3*h_num]*s[j+3*h_num]		+B[(i+h_num)*Nx+j+4*h_num]*s[j+4*h_num]		+B[(i+h_num)*Nx+j+5*h_num]*s[j+5*h_num];
							Bs[i+2*h_num]+=	B[(i+2*h_num)*Nx+j]*s[j]+B[(i+2*h_num)*Nx+j+h_num]*s[j+h_num]	+B[(i+2*h_num)*Nx+j+2*h_num]*s[j+2*h_num]	+B[(i+2*h_num)*Nx+j+3*h_num]*s[j+3*h_num]	+B[(i+2*h_num)*Nx+j+4*h_num]*s[j+4*h_num]	+B[(i+2*h_num)*Nx+j+5*h_num]*s[j+5*h_num];
							Bs[i+3*h_num]+=	B[(i+3*h_num)*Nx+j]*s[j]+B[(i+3*h_num)*Nx+j+h_num]*s[j+h_num]	+B[(i+3*h_num)*Nx+j+2*h_num]*s[j+2*h_num]	+B[(i+3*h_num)*Nx+j+3*h_num]*s[j+3*h_num]	+B[(i+3*h_num)*Nx+j+4*h_num]*s[j+4*h_num]	+B[(i+3*h_num)*Nx+j+5*h_num]*s[j+5*h_num];
							Bs[i+4*h_num]+=	B[(i+4*h_num)*Nx+j]*s[j]+B[(i+4*h_num)*Nx+j+h_num]*s[j+h_num]	+B[(i+4*h_num)*Nx+j+2*h_num]*s[j+2*h_num]	+B[(i+4*h_num)*Nx+j+3*h_num]*s[j+3*h_num]	+B[(i+4*h_num)*Nx+j+4*h_num]*s[j+4*h_num]	+B[(i+4*h_num)*Nx+j+5*h_num]*s[j+5*h_num];
							Bs[i+5*h_num]+=	B[(i+5*h_num)*Nx+j]*s[j]+B[(i+5*h_num)*Nx+j+h_num]*s[j+h_num]	+B[(i+5*h_num)*Nx+j+2*h_num]*s[j+2*h_num]	+B[(i+5*h_num)*Nx+j+3*h_num]*s[j+3*h_num]	+B[(i+5*h_num)*Nx+j+4*h_num]*s[j+4*h_num]	+B[(i+5*h_num)*Nx+j+5*h_num]*s[j+5*h_num];
						}
					}
					for(int i=0;i<Nx;i++)
					{
						for(int j=0;j<Nx;j++)
						{
							BssB[i*Nx+j]=0;
							for(int k=0;k<Nx;k++)	BssB[i*Nx+j]+=Bs[i]*s[k]*B[k*Nx+j];
							B[i*Nx+j]+=1/beta*y[i]*y[j]-1/sigma*BssB[i*Nx+j];
						}
					}
				}

				for(int i=0;i<h_num;i++)
				{
					delete[]	dp_k[i];
					delete[]	dq_k[i];
					delete[]	dpTr_k[i];
					delete[]	dqTr_k[i];
					delete[]	dp_a[i];
					delete[]	dq_a[i];
				}
				delete[]	dp_k;
				delete[]	dq_k;
				delete[]	dpTr_k;
				delete[]	dqTr_k;
				delete[]	dp_a;
				delete[]	dq_a;
				delete[]	w_a;
				delete[]	g_a;
				delete[]	h_a;

				delete[]	s;
				delete[]	y;
				delete[]	sB;
				delete[]	Bs;
				delete[]	BssB;
			}
		}

		for(int i=0;i<h_num;i++)
		{
			E_min+=sqrt((old_dp[i][A_X]-dp[i][A_X])*(old_dp[i][A_X]-dp[i][A_X]) + (old_dp[i][A_Y]-dp[i][A_Y])*(old_dp[i][A_Y]-dp[i][A_Y]) + (old_dp[i][A_Z]-dp[i][A_Z])*(old_dp[i][A_Z]-dp[i][A_Z])
				+ (old_dq[i][A_X]-dq[i][A_X])*(old_dq[i][A_X]-dq[i][A_X]) + (old_dq[i][A_Y]-dq[i][A_Y])*(old_dq[i][A_Y]-dq[i][A_Y]) + (old_dq[i][A_Z]-dq[i][A_Z])*(old_dq[i][A_Z]-dq[i][A_Z]));

			seta_g[i]+=g[i];
			seta_h[i]+=h[i];
		}
		if(E_min<ep*1000)	r*=4;
	}


	for(int i=0;i<h_num;i++)
	{
		HYPER[i].p[A_X]+=Dt*dp[i][A_X];
		HYPER[i].p[A_Y]+=Dt*dp[i][A_Y];
		HYPER[i].p[A_Z]+=Dt*dp[i][A_Z];
		PART[i].r[A_X]+=Dt*dq[i][A_X];
		PART[i].r[A_Y]+=Dt*dq[i][A_Y];
		PART[i].r[A_Z]+=Dt*dq[i][A_Z];
	}


	for(int i=0;i<h_num;i++)
	{
		delete[]	dp[i];
		delete[]	old_dp[i];

		delete[]	dq[i];
		delete[]	old_dq[i];

		delete[]	dpTr[i];
		delete[]	dqTr[i];
		delete[]	F[i];
		delete[]	ti_F[i];
		delete[]	S[i];
		delete[]	d_p[i];
		delete[]	d_q[i];

	}

	delete[]	dp;
	delete[]	old_dp;
	delete[]	dq;
	delete[]	old_dq;
	delete[]	dpTr;
	delete[]	dqTr;
	delete[]	J;
	delete[]	F;
	delete[]	ti_F;
	delete[]	S;
	delete[]	w;
	delete[]	g;
	delete[]	h;
	delete[]	seta_g;
	delete[]	seta_h;
	delete[]	d_p;
	delete[]	d_q;
	delete[]	B;
	delete[]	Nr;
}

void calc_wg(mpsconfig &CON, vector<mpselastic> PART, vector<hyperelastic> &HYPER, vector<hyperelastic2> &HYPER1, double **dq, double *w, double *g, double *J, double **F, double **ti_F, double **S)
{

//	cout<<"Fi�v�Z";
	////Fi�̍X�V
	int h_num=HYPER.size();
	double V=get_volume(&CON);
	double mi=V*CON.get_h_dis();
	double Dt=CON.get_dt();

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
			
			fi[0][0]+=w*(PART[inn].r[A_X]+Dt*dq[inn][A_X] - PART[i].r[A_X]-Dt*dq[i][A_X])*a[A_X];
			fi[0][1]+=w*(PART[inn].r[A_X]+Dt*dq[inn][A_X] - PART[i].r[A_X]-Dt*dq[i][A_X])*a[A_Y];
			fi[0][2]+=w*(PART[inn].r[A_X]+Dt*dq[inn][A_X] - PART[i].r[A_X]-Dt*dq[i][A_X])*a[A_Z];
			
			fi[1][0]+=w*(PART[inn].r[A_Y]+Dt*dq[inn][A_Y] - PART[i].r[A_Y]-Dt*dq[i][A_Y])*a[A_X];
			fi[1][1]+=w*(PART[inn].r[A_Y]+Dt*dq[inn][A_Y] - PART[i].r[A_Y]-Dt*dq[i][A_Y])*a[A_Y];
			fi[1][2]+=w*(PART[inn].r[A_Y]+Dt*dq[inn][A_Y] - PART[i].r[A_Y]-Dt*dq[i][A_Y])*a[A_Z];

			fi[2][0]+=w*(PART[inn].r[A_Z]+Dt*dq[inn][A_Z] - PART[i].r[A_Z]-Dt*dq[i][A_Z])*a[A_X];
			fi[2][1]+=w*(PART[inn].r[A_Z]+Dt*dq[inn][A_Z] - PART[i].r[A_Z]-Dt*dq[i][A_Z])*a[A_Y];
			fi[2][2]+=w*(PART[inn].r[A_Z]+Dt*dq[inn][A_Z] - PART[i].r[A_Z]-Dt*dq[i][A_Z])*a[A_Z];
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
		
		F[i][0*DIMENSION+0]=p_Fi[0][0];	F[i][0*DIMENSION+1]=p_Fi[0][1];	F[i][0*DIMENSION+2]=p_Fi[0][2];	
		F[i][1*DIMENSION+0]=p_Fi[1][0];	F[i][1*DIMENSION+1]=p_Fi[1][1];	F[i][1*DIMENSION+2]=p_Fi[1][2];	
		F[i][2*DIMENSION+0]=p_Fi[2][0];	F[i][2*DIMENSION+1]=p_Fi[2][1];	F[i][2*DIMENSION+2]=p_Fi[2][2];	
	
		//J�̌v�Z
		J[i]=calc_det3(p_Fi);
	//	for(int i=0;i<h_num;i++)	cout<<"J["<<i<<"]="<<J<<endl;
		g[i]=V*(1-J[i]);
		//t_inverse_Fi�̌v�Z
		inverse(p_Fi,DIMENSION);
		ti_F[i][0*DIMENSION+0]=p_Fi[0][0];	ti_F[i][0*DIMENSION+1]=p_Fi[1][0];	ti_F[i][0*DIMENSION+2]=p_Fi[2][0];
		ti_F[i][1*DIMENSION+0]=p_Fi[0][1];	ti_F[i][1*DIMENSION+1]=p_Fi[1][1];	ti_F[i][1*DIMENSION+2]=p_Fi[2][1];
		ti_F[i][2*DIMENSION+0]=p_Fi[0][2];	ti_F[i][2*DIMENSION+1]=p_Fi[1][2];	ti_F[i][2*DIMENSION+2]=p_Fi[2][2];
	}
	//	cout<<"----------OK"<<endl;

	for(int D=0;D>DIMENSION;D++)	delete[]	p_Fi[D];
	delete[]	p_Fi;




	double c10=CON.get_c10();
	double c01=CON.get_c01();

	double **dC=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	dC[D]=new double [DIMENSION];

	double dC2[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

	
	for(int i=0;i<h_num;i++)
	{
		double dFi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

		if(J[i]<0){
			dFi[0][0]=-1/pow(-J[i],1/3)*F[i][0*DIMENSION+0];	dFi[0][1]=-1/pow(-J[i],1/3)*F[i][0*DIMENSION+1];	dFi[0][2]=-1/pow(-J[i],1/3)*F[i][0*DIMENSION+2];
			dFi[1][0]=-1/pow(-J[i],1/3)*F[i][1*DIMENSION+0];	dFi[1][1]=-1/pow(-J[i],1/3)*F[i][1*DIMENSION+1];	dFi[1][2]=-1/pow(-J[i],1/3)*F[i][1*DIMENSION+2];
			dFi[2][0]=-1/pow(-J[i],1/3)*F[i][2*DIMENSION+0];	dFi[2][1]=-1/pow(-J[i],1/3)*F[i][2*DIMENSION+1];	dFi[2][2]=-1/pow(-J[i],1/3)*F[i][2*DIMENSION+2];
		}
		else
		{
			dFi[0][0]=1/pow(J[i],1/3)*F[i][0*DIMENSION+0];	dFi[0][1]=1/pow(J[i],1/3)*F[i][0*DIMENSION+1];	dFi[0][2]=1/pow(J[i],1/3)*F[i][0*DIMENSION+2];
			dFi[1][0]=1/pow(J[i],1/3)*F[i][1*DIMENSION+0];	dFi[1][1]=1/pow(J[i],1/3)*F[i][1*DIMENSION+1];	dFi[1][2]=1/pow(J[i],1/3)*F[i][1*DIMENSION+2];
			dFi[2][0]=1/pow(J[i],1/3)*F[i][2*DIMENSION+0];	dFi[2][1]=1/pow(J[i],1/3)*F[i][2*DIMENSION+1];	dFi[2][2]=1/pow(J[i],1/3)*F[i][2*DIMENSION+2];
		}


		dC[0][0]=dFi[0][0]*dFi[0][0]+dFi[1][0]*dFi[1][0]+dFi[2][0]*dFi[2][0];
		dC[0][1]=dFi[0][0]*dFi[0][1]+dFi[1][0]*dFi[1][1]+dFi[2][0]*dFi[2][1];
		dC[0][2]=dFi[0][0]*dFi[0][2]+dFi[1][0]*dFi[1][2]+dFi[2][0]*dFi[2][2];
		dC[1][0]=dFi[0][1]*dFi[0][0]+dFi[1][1]*dFi[1][0]+dFi[2][1]*dFi[2][0];
		dC[1][1]=dFi[0][1]*dFi[0][1]+dFi[1][1]*dFi[1][1]+dFi[2][1]*dFi[2][1];
		dC[1][2]=dFi[0][1]*dFi[0][2]+dFi[1][1]*dFi[1][2]+dFi[2][1]*dFi[2][2];
		dC[2][0]=dFi[0][2]*dFi[0][0]+dFi[1][2]*dFi[1][0]+dFi[2][2]*dFi[2][0];
		dC[2][1]=dFi[0][2]*dFi[0][1]+dFi[1][2]*dFi[1][1]+dFi[2][2]*dFi[2][1];
		dC[2][2]=dFi[0][2]*dFi[0][2]+dFi[1][2]*dFi[1][2]+dFi[2][2]*dFi[2][2];


		double trace_dC=dC[0][0]+dC[1][1]+dC[2][2];

		dC2[0][0]=dC[A_X][0]*dC[0][A_X]+dC[A_X][1]*dC[1][A_X]+dC[A_X][2]*dC[2][A_X];
		dC2[0][1]=dC[A_X][0]*dC[0][A_Y]+dC[A_X][1]*dC[1][A_Y]+dC[A_X][2]*dC[2][A_Y];
		dC2[0][2]=dC[A_X][0]*dC[0][A_Z]+dC[A_X][1]*dC[1][A_Z]+dC[A_X][2]*dC[2][A_Z];
		dC2[1][0]=dC[A_Y][0]*dC[0][A_X]+dC[A_Y][1]*dC[1][A_X]+dC[A_Y][2]*dC[2][A_X];
		dC2[1][1]=dC[A_Y][0]*dC[0][A_Y]+dC[A_Y][1]*dC[1][A_Y]+dC[A_Y][2]*dC[2][A_Y];
		dC2[1][2]=dC[A_Y][0]*dC[0][A_Z]+dC[A_Y][1]*dC[1][A_Z]+dC[A_Y][2]*dC[2][A_Z];
		dC2[2][0]=dC[A_Z][0]*dC[0][A_X]+dC[A_Z][1]*dC[1][A_X]+dC[A_Z][2]*dC[2][A_X];
		dC2[2][1]=dC[A_Z][0]*dC[0][A_Y]+dC[A_Z][1]*dC[1][A_Y]+dC[A_Z][2]*dC[2][A_Y];
		dC2[2][2]=dC[A_Z][0]*dC[0][A_Z]+dC[A_Z][1]*dC[1][A_Z]+dC[A_Z][2]*dC[2][A_Z];

		double trace_dC2=dC2[0][0]+dC2[1][1]+dC2[2][2];
		double Ic=trace_dC;
		double IIc=0.50*(trace_dC*trace_dC-trace_dC2);

		w[i]=c10*(Ic-3)+c01*(IIc-3);

		double C[DIMENSION][DIMENSION]={{dC[0][0],dC[0][1],dC[0][2]},{dC[1][0],dC[1][1],dC[1][2]},{dC[2][0],dC[2][1],dC[2][2]}};

		inverse(dC,DIMENSION);

		if(J[i]<0){
			S[i][0*DIMENSION+0]=-1/pow(-J[i],2/3)*( (c10+c01*Ic)*(1-1/3*trace_dC*dC[0][0]) -c01*(C[0][0]-1/3*trace_dC2*dC[0][0]) );
			S[i][0*DIMENSION+1]=-1/pow(-J[i],2/3)*( (c10+c01*Ic)*( -1/3*trace_dC*dC[0][1]) -c01*(C[0][1]-1/3*trace_dC2*dC[0][1]) );
			S[i][0*DIMENSION+2]=-1/pow(-J[i],2/3)*( (c10+c01*Ic)*( -1/3*trace_dC*dC[0][2]) -c01*(C[0][2]-1/3*trace_dC2*dC[0][2]) );

			S[i][1*DIMENSION+0]=-1/pow(-J[i],2/3)*( (c10+c01*Ic)*( -1/3*trace_dC*dC[1][0]) -c01*(C[1][0]-1/3*trace_dC2*dC[1][0]) );
			S[i][1*DIMENSION+1]=-1/pow(-J[i],2/3)*( (c10+c01*Ic)*(1-1/3*trace_dC*dC[1][1]) -c01*(C[1][1]-1/3*trace_dC2*dC[1][1]) );
			S[i][1*DIMENSION+2]=-1/pow(-J[i],2/3)*( (c10+c01*Ic)*( -1/3*trace_dC*dC[1][2]) -c01*(C[1][2]-1/3*trace_dC2*dC[1][2]) );

			S[i][2*DIMENSION+0]=-1/pow(-J[i],2/3)*( (c10+c01*Ic)*( -1/3*trace_dC*dC[2][0]) -c01*(C[2][0]-1/3*trace_dC2*dC[2][0]) );
			S[i][2*DIMENSION+1]=-1/pow(-J[i],2/3)*( (c10+c01*Ic)*( -1/3*trace_dC*dC[2][1]) -c01*(C[2][1]-1/3*trace_dC2*dC[2][1]) );
			S[i][2*DIMENSION+2]=-1/pow(-J[i],2/3)*( (c10+c01*Ic)*(1-1/3*trace_dC*dC[2][2]) -c01*(C[2][2]-1/3*trace_dC2*dC[2][2]) );
		}
		else
		{
			S[i][A_X*DIMENSION+0]=1/pow(J[i],2/3)*( (c10+c01*Ic)*(1-1/3*trace_dC*dC[0][0]) -c01*(C[0][0]-1/3*trace_dC2*dC[0][0]) );
			S[i][A_X*DIMENSION+1]=1/pow(J[i],2/3)*( (c10+c01*Ic)*( -1/3*trace_dC*dC[0][1]) -c01*(C[0][1]-1/3*trace_dC2*dC[0][1]) );
			S[i][A_X*DIMENSION+2]=1/pow(J[i],2/3)*( (c10+c01*Ic)*( -1/3*trace_dC*dC[0][2]) -c01*(C[0][2]-1/3*trace_dC2*dC[0][2]) );

			S[i][1*DIMENSION+0]=1/pow(J[i],2/3)*( (c10+c01*Ic)*( -1/3*trace_dC*dC[1][0]) -c01*(C[1][0]-1/3*trace_dC2*dC[1][0]) );
			S[i][1*DIMENSION+1]=1/pow(J[i],2/3)*( (c10+c01*Ic)*(1-1/3*trace_dC*dC[1][1]) -c01*(C[1][1]-1/3*trace_dC2*dC[1][1]) );
			S[i][1*DIMENSION+2]=1/pow(J[i],2/3)*( (c10+c01*Ic)*( -1/3*trace_dC*dC[1][2]) -c01*(C[1][2]-1/3*trace_dC2*dC[1][2]) );

			S[i][2*DIMENSION+0]=1/pow(J[i],2/3)*( (c10+c01*Ic)*( -1/3*trace_dC*dC[2][0]) -c01*(C[2][0]-1/3*trace_dC2*dC[2][0]) );
			S[i][2*DIMENSION+1]=1/pow(J[i],2/3)*( (c10+c01*Ic)*( -1/3*trace_dC*dC[2][1]) -c01*(C[2][1]-1/3*trace_dC2*dC[2][1]) );
			S[i][2*DIMENSION+2]=1/pow(J[i],2/3)*( (c10+c01*Ic)*(1-1/3*trace_dC*dC[2][2]) -c01*(C[2][2]-1/3*trace_dC2*dC[2][2]) );
		}
	}
	for(int D=0;D>DIMENSION;D++)	delete[]	dC[D];
	delete[]	dC;
}


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

