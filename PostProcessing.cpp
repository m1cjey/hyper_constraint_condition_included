#include "stdafx.h"
#include "Micro_AVS.h"

void post_processing(mpsconfig &CON, vector<mpselastic> &PART, elastic &ELAST, int fluid_number,int particle_number,double dt,double Umax, int t,double TIME, double **F)
{
	double le=CON.get_distancebp();
	unsigned int timeA=GetTickCount();	//�v�Z�J�n����

	cout<<"�e�����ʊJ�n: ";
	
	//AVS�ɗ��q�f�[�^�o��
	if(t==1 || t%CON.get_interval()==0) particle_movie_AVS(CON,PART,ELAST,fluid_number,particle_number,t,TIME);
		cout<<"AVS�o�͊���"<<endl;
		
	//AVS2�Ɏ������x�A���[�����c�͏o��
		if((t==1 || t%CON.get_EM_interval()==0) && CON.get_FEM_flag()==ON) particle_movie_AVS2(CON,t,PART,TIME,fluid_number,particle_number,F);
		cout<<"AVS2�o�͊���\n";

	//���x���v���b�g
	if(CON.get_current_step()==1 || CON.get_current_step()%CON.get_interval()==0) plot_speed(CON ,PART,particle_number,fluid_number);
		cout<<"���x�v���b�g����"<<endl;

	//���͂ɂ������x���v���b�g
	if(CON.get_current_step()==1 || CON.get_current_step()%CON.get_interval()==0) plot_pressure_acceleration(CON, PART);
		cout<<"�����x�v���b�g����"<<endl;

	//���͕��z���v���b�g
	if(CON.get_flag_ELAST()==ON)	if(CON.get_current_step()==1 || CON.get_current_step()%CON.get_interval()==0) plot_stress_distribution(CON, PART);
		cout<<"���͕��z�v���b�g����"<<endl;

	//����f�͂��v���b�g
	if(CON.get_flag_ELAST()==ON)	if(CON.get_current_step()==1 || CON.get_current_step()%CON.get_interval()==0) plot_shear_force(CON, PART);
		cout<<"����f�����x�v���b�g����"<<endl;

	//�Ђ��ݑ��x�ɂ������x���v���b�g
	if(CON.get_flag_ELAST()==ON)	if(CON.get_current_step()==1 || CON.get_current_step()%CON.get_interval()==0) plot_strainrate_acceleration(CON, PART);
		cout<<"�Ђ��݉����x�v���b�g����"<<endl;

	//�Ђ��ݗ����v���b�g
	if(CON.get_current_step()==1 || CON.get_current_step()%CON.get_interval()==0) plot_distortion_rate(CON, PART);
		cout<<"�Ђ��ݗ��v���b�g����"<<endl;

	//�c�������x���v���b�g
	if(CON.get_flag_ELAST()==ON)	if(CON.get_current_step()==1 || CON.get_current_step()%CON.get_interval()==0) plot_residual_acceleration(CON, PART, F);
		cout<<"�c�������x�v���b�g����"<<endl;

	//���W�v���b�g
	ofstream gnu1("0.dat");//��͏I����̑S���q���W���v���b�g
	ofstream gnu2("suf.dat");//�\�ʗ��q�������v���b�g
	
	if(CON.get_dimension()==2)
	{
		for(int i=0;i<particle_number;i++)
		{ 
			gnu1<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<PART[i].r[A_Z]<<endl;
			if(PART[i].surface==ON)
			{
				gnu2<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<PART[i].r[A_Z]<<endl;
			}
		}
	}
	else if(CON.get_dimension()==3)
	{
		for(int i=0;i<particle_number;i++)
		{ 
			if(PART[i].surface==ON && (PART[i].type==FLUID || PART[i].type==ELASTIC || PART[i].type==MAGELAST||PART[i].type==HYPERELAST))//���̂����\���������Ƃ�
			{
				gnu2<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<PART[i].r[A_Z]<<endl;
			}
		}
	}
	gnu1.close();
	gnu2.close();
	//*/
		
	if(CON.get_flag_ELAST()==ON)
	{
		//���ϗ��q���x&���͂�\��
		double ave_n0=0;
		double ave_P=0;
		int count=0;
		for(int i=0;i<fluid_number;i++) 
		{
			if(PART[i].surface==OFF)
			{
				ave_n0+=PART[i].PND; //���ϗ��q�����x
				ave_P+=fabs(PART[i].P); //���ψ���
				count++;
			}
		}
		if(count!=0){ 
			ave_n0/=count;
			ave_P/=count;
		}

		cout<<"average n0="<<ave_n0<<" average P="<<ave_P<<"/"<<CON.get_ave_P_for_FEM_flag()<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
	
		check_FEM_flag(CON, ELAST, ave_P);
	}
	CON.set_current_step(t);
}

void plot_pressure_acceleration(mpsconfig &CON, vector<mpselastic> &PART)
{
	double le=CON.get_distancebp()*0.5;
	double times=1;//CON.get_pressure_times();0.001
	int dimension=CON.get_dimension();
	int face=CON.get_speed_face();			//3D��͎���pressure.dat�̏o�͖� 0=YZ���� 1=XZ 2=XY
	double face_p=CON.get_speed_face_p();	//3D��͎���pressure.dat�̏o�͖ʂ̍��W
	int d1,d2,d3;							//3D��͎��̏o�͂ɕK�v�Ȏ���
	double xmax=-100;						//�o�͗��q�̍ő剡���W
	double ymax=-100;						//�o�͗��q�̍ő�c���W
	
	stringstream ss;
	ss<<"./Pressure"<<"/pressure_accel"<<CON.get_current_step()<<".dat";
	string filename=ss.str();
	
	ofstream vec(filename);

	if(vec.fail()){
		system("mkdir Pressure");
		ofstream vec(filename);//�Ď��s
		if(vec.fail()){
			cout<<"Pressure�f�B���N�g�����쐬�ł��܂���ł���"<<endl;
			exit(1);
		}
	}
	stringstream ss2;
	ss2<<"./Pressure"<<"/pressure3D_"<<CON.get_current_step()<<".dat";
	string filename2=ss2.str();
	
	ofstream ve(filename2);

	if(ve.fail()){
		system("mkdir Pressure");
		ofstream ve(filename2);//�Ď��s
		if(ve.fail()){
			cout<<"Pressure�f�B���N�g�����쐬�ł��܂���ł���"<<endl;
			exit(1);
		}
	}
	
	if(dimension==3)
	{
		//int d1,d2;				//�o�͂ɕK�v�Ȏ���
		if(face==0){
			d1=A_Y; d2=A_Z; d3=A_X;
		}
		else if(face==1)
		{
			d1=A_X; d2=A_Z; d3=A_Y;
		}

		for(int i=0;i<PART.size();i++)
		{
			if(PART[i].r[face]>face_p-le && PART[i].r[face]<face_p+le)
			{
				double x=PART[i].r[d1];
				double z=PART[i].r[d2];
				double u=PART[i].get_pressure_accel(d1);
				double w=PART[i].get_pressure_accel(d2);
				vec<<x<<" "<<z<<" "<<u*times<<" "<<w*times<<endl;
				ve<<x<<" "<<z<<" "<<PART[i].P<<endl;
				if(x>xmax) xmax=x;
				if(z>ymax) ymax=z;
			}
		}
	}
	xmax+=4*le;//�}����o���ʒu��ی��������ď����΂߂Ɉړ�
	ymax+=4*le;
	if(CON.get_legend_pressure()>0) vec<<xmax<<" "<<ymax<<" "<<CON.get_legend_pressure()*times<<" "<<0*times<<endl;//�Ō�ɖ}��o��
	vec.close();
	ve.close();
}

void plot_stress_distribution(mpsconfig &CON, vector<mpselastic> &PART)
{
	Micro_AVS Micro_avs;
	//�\���i�q�f�[�^�t�@�C��
	int timestep=CON.get_current_step();
	double den=CON.get_density();

	stringstream sstre;
	sstre<<"./Stress/Stress"<<timestep<<".fld";
	string stForce=sstre.str();

	ofstream fout3(stForce);
	if(fout3.fail()){
		system("mkdir Stress");
		ofstream fout3(stForce);
		if(fout3.fail()){
			cout<<"./Stress�t�H���_���J���܂���ł���"<<endl;
			exit(1);
		}
	}

	fout3 << "# AVS field file" << endl;
	fout3 << "ndim=1" << endl;
	fout3 << "dim1=" << PART.size() <<endl;
	fout3 << "nspace=3" << endl;
	fout3 << "veclen=3" << endl;
	fout3 << "data=float" << endl;
	fout3 << "field=irregular" << endl;
	fout3 << "label=e-x e-y e-z" << endl << endl;
	fout3 << "variable 1 file=./Stress"<<timestep<<" filetype=ascii skip=1 offset=0 stride=6" << endl;
	fout3 << "variable 2 file=./Stress"<<timestep<<" filetype=ascii skip=1 offset=1 stride=6" << endl;
	fout3 << "variable 3 file=./Stress"<<timestep<<" filetype=ascii skip=1 offset=2 stride=6" << endl;
	fout3 << "coord    1 file=./Stress"<<timestep<<" filetype=ascii skip=1 offset=3 stride=6" << endl;
	fout3 << "coord    2 file=./Stress"<<timestep<<" filetype=ascii skip=1 offset=4 stride=6" << endl;
	fout3 << "coord    3 file=./Stress"<<timestep<<" filetype=ascii skip=1 offset=5 stride=6" << endl;
	fout3.close();

	//�f�[�^�t�@�C��
	stringstream ss;
	ss<<"./Stress/Stress"<<timestep;
	string filename=ss.str();

	ofstream fout(filename);
	if(fout.fail()){
		cout<<"./Stress�t�H���_���J���܂���ł���"<<endl;
		exit(1);
	}

	fout<<"e-x e-y e-z x y z"<<endl;
	for(size_t i=0;i<PART.size();i++)
    {
		if(PART[i].type!=WALL){
//			Micro_avs.make_list(PART[i].r[A_X], PART[i].r[A_Y], PART[i].r[A_Z], PART[i].get_stress_accel(A_X)*den, PART[i].get_stress_accel(A_Y)*den, PART[i].get_stress_accel(A_Z)*den);
//		fout<<F[A_X][i]*times<<" "<<F[A_Y][i]*times<<" "<<F[A_Z][i]*times<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
		fout<<PART[i].get_stress_accel(A_X)*den<<" "<<PART[i].get_stress_accel(A_Y)*den<<" "<<PART[i].get_stress_accel(A_Z)*den<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
		}
	}
//	Micro_avs.Output_vector_MicroAVS("Stress",CON.get_current_step());
	fout.close();
}


void plot_shear_force(mpsconfig &CON, vector<mpselastic> &PART)
{
	Micro_AVS Micro_avs;
	double le=CON.get_distancebp()*0.5;
	double mass=CON.get_particle_mass();
	double times=CON.get_pressure_times();
	int dimension=CON.get_dimension();
	int face=CON.get_speed_face();			//3D��͎���pressure.dat�̏o�͖� 0=YZ���� 1=XZ 2=XY
	double face_p=CON.get_speed_face_p();	//3D��͎���pressure.dat�̏o�͖ʂ̍��W
	int d1,d2,d3;							//3D��͎��̏o�͂ɕK�v�Ȏ���
	double xmax=-100;						//�o�͗��q�̍ő剡���W
	double ymax=-100;						//�o�͗��q�̍ő�c���W
	
/*	stringstream ss;
	ss<<"./Shear"<<"/shear"<<CON.get_current_step()<<".dat";
	string filename=ss.str();
	
	ofstream vec(filename);

	if(vec.fail()){
		system("mkdir Shear");
		ofstream vec(filename);//�Ď��s
		if(vec.fail()){
			cout<<"Shear�f�B���N�g�����쐬�ł��܂���ł���"<<endl;
			exit(1);
		}
	}*/
	
	if(dimension==3)
	{
		//int d1,d2;				//�o�͂ɕK�v�Ȏ���
		if(face==0){
			d1=A_Y; d2=A_Z; d3=A_X;
		}
		else if(face==1)
		{
			d1=A_X; d2=A_Z; d3=A_Y;
		}
		double u_max=0;
		double w_max=0;
		for(int i=0;i<PART.size();i++)
		{
	//		if(PART[i].r[face]>face_p-le && PART[i].r[face]<face_p+le)
	//		{
				Micro_avs.make_list(PART[i].r[A_X],PART[i].r[A_Y],PART[i].r[A_Z],PART[i].get_stress_accel(A_X)*mass,PART[i].get_stress_accel(A_Y)*mass,PART[i].get_stress_accel(A_Z)*mass);
	/*			double x=PART[i].r[d1];
				double z=PART[i].r[d2];
				double u=PART[i].get_stress_accel(d1)*mass;
				double w=PART[i].get_stress_accel(d2)*mass;
				
				if(u_max<u)u_max=u;
				if(w_max<w)w_max=w;*/

	//			vec<<x<<" "<<z<<" "<<u<<" "<<w<<endl;
				
		//	}

		}
	}
	Micro_avs.Output_vector_MicroAVS("Shear",CON.get_current_step());
/*	xmax+=4*le;//�}����o���ʒu��ی��������ď����΂߂Ɉړ�
	ymax+=4*le;
	if(CON.get_legend_pressure()>0) vec<<xmax<<" "<<ymax<<" "<<CON.get_legend_pressure()*times<<" "<<0*times<<endl;//�Ō�ɖ}��o��
	vec.close();*/
}

void plot_strainrate_acceleration(mpsconfig &CON, vector<mpselastic> &PART)
{
	Micro_AVS Micro_avs;
	double le=CON.get_distancebp()*0.5;
	double times=CON.get_pressure_times();
	int dimension=CON.get_dimension();
	int face=CON.get_speed_face();			//3D��͎���pressure.dat�̏o�͖� 0=YZ���� 1=XZ 2=XY
	double face_p=CON.get_speed_face_p();	//3D��͎���pressure.dat�̏o�͖ʂ̍��W
	int d1,d2,d3;							//3D��͎��̏o�͂ɕK�v�Ȏ���
	double xmax=-100;						//�o�͗��q�̍ő剡���W
	double ymax=-100;						//�o�͗��q�̍ő�c���W

	double g[3]={0.0, 0.0, 0.0};
	if(dimension==2) g[A_Y]=CON.get_g(); 
	if(dimension==3) g[A_Z]=CON.get_g();

	double mass=CON.get_particle_mass();

//	stringstream ss;
//	ss<<"./StrainRate"<<"/strainrate"<<CON.get_current_step()<<".dat";
//	string filename=ss.str();
	
//	ofstream vec(filename);

//	if(vec.fail()){
//		system("mkdir StrainRate");
//		ofstream vec(filename);//�Ď��s
//		if(vec.fail()){
//			cout<<"StrainRate�f�B���N�g�����쐬�ł��܂���ł���"<<endl;
//			exit(1);
//		}
//	}
	
	if(dimension==3)
	{
		//int d1,d2;				//�o�͂ɕK�v�Ȏ���
		if(face==0){
			d1=A_Y; d2=A_Z; d3=A_X;
		}
		else if(face==1)
		{
			d1=A_X; d2=A_Z; d3=A_Y;
		}

		for(int i=0;i<PART.size();i++)
		{
			if(PART[i].r[face]>face_p-le && PART[i].r[face]<face_p+le)
			{
				Micro_avs.make_list(PART[i].r[d1],PART[i].r[d2],PART[i].get_stress_visco_accel(d1),PART[i].get_stress_visco_accel(d2));
/*				double x=PART[i].r[d1];
				double z=PART[i].r[d2];
				double u=PART[i].get_stress_visco_accel(d1);
				double w=PART[i].get_stress_visco_accel(d2);
				vec<<x<<" "<<z<<" "<<u*times<<" "<<w*times<<endl;
				if(x>xmax) xmax=x;
				if(z>ymax) ymax=z;*/
			}
		}
	}
	Micro_avs.Output_vector_MicroAVS("StrainRate",CON.get_current_step());
//	xmax+=4*le;//�}����o���ʒu��ی��������ď����΂߂Ɉړ�
//	ymax+=4*le;
//	if(CON.get_legend_pressure()>0) vec<<xmax<<" "<<ymax<<" "<<CON.get_legend_pressure()*times<<" "<<0*times<<endl;//�Ō�ɖ}��o��
//	vec.close();
}

void plot_distortion_rate(mpsconfig &CON, vector<mpselastic> &PART)
{
	Micro_AVS Micro_avs;
	double le=CON.get_distancebp()*0.5;
	double times=CON.get_pressure_times();
	int dimension=CON.get_dimension();
	int face=CON.get_speed_face();			//3D��͎���pressure.dat�̏o�͖� 0=YZ���� 1=XZ 2=XY
	double face_p=CON.get_speed_face_p();	//3D��͎���pressure.dat�̏o�͖ʂ̍��W
	int d1,d2,d3;							//3D��͎��̏o�͂ɕK�v�Ȏ���
	double xmax=-100;						//�o�͗��q�̍ő剡���W
	double ymax=-100;						//�o�͗��q�̍ő�c���W

	double g[3]={0.0, 0.0, 0.0};
	if(dimension==2) g[A_Y]=CON.get_g(); 
	if(dimension==3) g[A_Z]=CON.get_g();

	double mass=CON.get_particle_mass();

/*	stringstream ss;
	ss<<"./distortionRate"<<"/distortionRate"<<CON.get_current_step()<<".dat";
	string filename=ss.str();
	
	ofstream vec(filename);

	if(vec.fail()){
		system("mkdir distortionRate");
		ofstream vec(filename);//�Ď��s
		if(vec.fail()){
			cout<<"distortionRate�f�B���N�g�����쐬�ł��܂���ł���"<<endl;
			exit(1);
		}
	}*/
	
	if(dimension==3)
	{
		//int d1,d2;				//�o�͂ɕK�v�Ȏ���
		if(face==0){
			d1=A_Y; d2=A_Z; d3=A_X;
		}
		else if(face==1)
		{
			d1=A_X; d2=A_Z; d3=A_Y;
		}

		for(int i=0;i<PART.size();i++)
		{
			if(PART[i].r[face]>face_p-le && PART[i].r[face]<face_p+le)
			{
				Micro_avs.make_list(PART[i].r[d1],PART[i].r[d2],PART[i].volumetric_strain,PART[i].get_stress_visco_accel(d2));
/*				double x=PART[i].r[d1];
				double z=PART[i].r[d2];
				double u=PART[i].volumetric_strain;
				double w=PART[i].get_stress_visco_accel(d2);
	//			vec<<x<<" "<<z<<" "<<u*times<<" "<<w*times<<endl;
				vec<<x<<" "<<z<<" "<<u<<endl;
				if(x>xmax) xmax=x;
				if(z>ymax) ymax=z;*/
			}
		}
	}
	Micro_avs.Output_vector_MicroAVS("distortionRate",CON.get_current_step());
//	vec.close();
}

void plot_residual_acceleration(mpsconfig &CON, vector<mpselastic> &PART, double **F)
{
	double le=CON.get_distancebp()*0.5;
	double times=CON.get_pressure_times();
	int dimension=CON.get_dimension();
	int face=CON.get_speed_face();			//3D��͎���pressure.dat�̏o�͖� 0=YZ���� 1=XZ 2=XY
	double face_p=CON.get_speed_face_p();	//3D��͎���pressure.dat�̏o�͖ʂ̍��W
	int d1,d2,d3;							//3D��͎��̏o�͂ɕK�v�Ȏ���
	double xmax=-100;						//�o�͗��q�̍ő剡���W
	double ymax=-100;						//�o�͗��q�̍ő�c���W

	double g[3]={0.0, 0.0, 0.0};
	if(dimension==2) g[A_Y]=CON.get_g(); 
	if(dimension==3) g[A_Z]=CON.get_g();

	double mass=CON.get_particle_mass();

	stringstream ss;
	ss<<"./Residual"<<"/residual"<<CON.get_current_step()<<".dat";
	string filename=ss.str();
	
	ofstream vec(filename);

	if(vec.fail()){
		system("mkdir Residual");
		ofstream vec(filename);//�Ď��s
		if(vec.fail()){
			cout<<"Residual�f�B���N�g�����쐬�ł��܂���ł���"<<endl;
			exit(1);
		}
	}
	
	if(dimension==3)
	{
		//int d1,d2;				//�o�͂ɕK�v�Ȏ���
		if(face==0){
			d1=A_Y; d2=A_Z; d3=A_X;
		}
		else if(face==1)
		{
			d1=A_X; d2=A_Z; d3=A_Y;
		}

		for(int i=0;i<PART.size();i++)
		{
			if(PART[i].r[face]>face_p-le && PART[i].r[face]<face_p+le)
			{
				double x=PART[i].r[d1];
				double z=PART[i].r[d2];
				double u=PART[i].get_stress_accel(d1)+PART[i].get_pressure_accel(d1)+PART[i].get_stress_visco_accel(d1)+g[d1]+F[d1][i]/mass;
				double w=PART[i].get_stress_accel(d2)+PART[i].get_pressure_accel(d2)+PART[i].get_stress_visco_accel(d2)+g[d2]+F[d2][i]/mass;
				vec<<x<<" "<<z<<" "<<u*times<<" "<<w*times<<endl;
				if(x>xmax) xmax=x;
				if(z>ymax) ymax=z;
			}
		}
	}
	xmax+=4*le;//�}����o���ʒu��ی��������ď����΂߂Ɉړ�
	ymax+=4*le;
	if(CON.get_legend_pressure()>0) vec<<xmax<<" "<<ymax<<" "<<CON.get_legend_pressure()*times<<" "<<0*times<<endl;//�Ō�ɖ}��o��
	vec.close();
}

void check_FEM_flag(mpsconfig &CON, elastic &ELAST, double ave_P)
{
	if(CON.get_current_time()<CON.get_dt()){//time��double�Ȃ̂�==�͎g��Ȃ�
		ofstream fout("aveave_P_history.txt", ios::out);
		fout<<CON.get_current_time()<<"\t"<<ave_P<<endl;
		fout.close();
	}else{
		ofstream fout("aveave_P_history.txt", ios::app);
		fout<<CON.get_current_time()<<"\t"<<ave_P<<endl;
		fout.close();
	}

	//���͂��w�肵���l�ɂȂ��FEM�X�C�b�`ON
	if(CON.get_FEM_flag()==false){
		if(ave_P>CON.get_ave_P_for_FEM_flag()){
			CON.set_FEM_flag(true);//�ǂ������ɂ��Ƃ��E�E�E
			ELAST.set_FEM_switch(true);
		}
	}
}

//quickMPS�p�|�X�g�����֐��i36�s�j
void post_processing3(mpsconfig &CON,vector<mpselastic> &PART,int fluid_number,int particle_number,int t,double TIME)
{
	
	//restart�p�t�@�C���o��
	if(t==CON.get_step() || t%CON.get_autosave()==0 )
	{
		//restart�p�ɗ��q���Ɨ��q�f�[�^���L�^
		ofstream hoge1("number.dat");
		hoge1<<particle_number<<endl;
		hoge1<<TIME<<endl;
		hoge1.close();
		
		FILE *hoge6;
		hoge6=fopen("restart_input.dat","w");//mps_input.dat�Ƃ͋�ʂ��Ȃ��ƁArestart�����s�����Ƃ�����
		if(hoge6==NULL){
			cout<<"�t�@�C���I�[�v���G���["<<endl; 
			exit(1);
		}
		for(int i=0;i<particle_number;i++)//���̉�͗��q�����ɋL�q
		{
			fprintf( hoge6, "%d\t",i);
			fprintf( hoge6, "%5.10f\t",PART[i].r[A_X]);
			fprintf( hoge6, "%5.10f\t",PART[i].r[A_Y]);
			fprintf( hoge6, "%5.10f\t",PART[i].r[A_Z]);
			fprintf( hoge6, "%5.10f\t",PART[i].u[A_X]);  //���xx����
			fprintf( hoge6, "%5.10f\t",PART[i].u[A_Y]);  //���xy����
			fprintf( hoge6, "%5.10f\t",PART[i].u[A_Z]);  //���xz����
			fprintf( hoge6, "%5.10f\t",PART[i].P); //����
			fprintf( hoge6, "%5.10f\t",PART[i].h); //�G���^���s�[
			fprintf( hoge6, "%5.10f\t",PART[i].val);
			fprintf( hoge6, "%d\n",PART[i].type);
			fprintf( hoge6, "%d\n",PART[i].materialID);
			fprintf( hoge6, "%d\n",PART[i].surface);
			fprintf( hoge6, "%d\n",PART[i].toFEM);
		}
		fclose(hoge6);
		//////////////////////////////////////////////
	}
}

//���x�v���b�g�֐��i200�s�j
void plot_speed(mpsconfig &CON ,vector<mpselastic> &PART,int particle_number,int fluid_number)
{
	Micro_AVS Micro_avs;
	double le=CON.get_distancebp()*0.5;
	double times=CON.get_speedtimes();
	int d=CON.get_dimension();
	int NUM;								//AVS�ɏo�͂��闱�q��
	int startID=0;							//�ŏ��ɏo�͂��闱�q��id
	int num=0;								//���������ϐ�
	int face=CON.get_speed_face();			//3D��͎���speed.dat�̏o�͖� 0=YZ���� 1=XZ 2=XY
	double face_p=CON.get_speed_face_p();	//3D��͎���speed.dat�̏o�͖ʂ̍��W
	int d1,d2,d3;								//3D��͎��̏o�͂ɕK�v�Ȏ���
	double xmax=-100;						//�o�͗��q�̍ő剡���W
	double ymax=-100;						//�o�͗��q�̍ő�c���W
	
	//AVS�o�͗��q��NUM�v�Z
	if(CON.get_speed_plot_particle()==1) NUM=particle_number;	//�S���q�o��
	else if(CON.get_speed_plot_particle()==2) NUM=fluid_number;//���̗��q�̂ݏo��
	else if(CON.get_speed_plot_particle()==3)					//�Ǘ��q�̂ݏo��
	{
		NUM=particle_number; 
		startID=fluid_number;
	}

	{
/*		stringstream ss;
		ss<<"./Speed"<<"/speed"<<CON.get_current_step()<<".dat";
		string filename=ss.str();
	
		ofstream vec(filename);//��Α��x

		if(vec.fail()){
			system("mkdir Speed");
			ofstream vec(filename);//�Ď��s
			if(vec.fail()){
				cout<<"Speed�f�B���N�g�����쐬�ł��܂���ł���"<<endl;
				exit(1);
			}
		}*/

		if(d==2)
		{
			for(int i=startID;i<NUM;i++)
			{
				Micro_avs.make_list(PART[i].r[A_X],PART[i].r[A_Y],PART[i].u[A_X],PART[i].u[A_Y]);
	/*			vec<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].u[A_X]*times<<" "<<PART[i].u[A_Y]*times<<endl;
				if(PART[i].r[A_X]>xmax) xmax=PART[i].r[A_X];
				if(PART[i].r[A_Y]>ymax) ymax=PART[i].r[A_Y];*/
			}
		}
		else if(d==3)
		{
			//int d1,d2;				//�o�͂ɕK�v�Ȏ���
			if(face==0) {d1=A_Y; d2=A_Z; d3=A_X;}
			else if(face==1) {d1=A_X; d2=A_Z; d3=A_Y;}
			if(CON.get_ax_sym_modify()==OFF)
			{
				for(int i=startID;i<NUM;i++)
				{
					if(PART[i].r[face]>face_p-le && PART[i].r[face]<face_p+le)
					{
						Micro_avs.make_list(PART[i].r[d1],PART[i].r[d2],PART[i].u[d1],PART[i].u[d2]);
					/*	double x=PART[i].r[d1];
						double z=PART[i].r[d2];
						double u=PART[i].u[d1];
						double w=PART[i].u[d2];
						vec<<x<<" "<<z<<" "<<u*times<<" "<<w*times<<endl;
						if(x>xmax) xmax=x;
						if(z>ymax) ymax=z;*/
					}
				}
			}
			else if(CON.get_ax_sym_modify()==ON)
			{
				for(int i=startID;i<NUM;i++)
				{
					if(PART[i].type!=WALL){
						if(PART[i].r[face]>face_p-le && PART[i].r[face]<face_p+le)
						{
							double x=PART[i].r[d1];//�o�͂Ɋ֗^���鎲 �֋X��A�ϐ�����x�ƂȂ��Ă��邪�A�����Ƃ͌���Ȃ����Ƃɒ���
							double z=PART[i].r[d2];//�o�͂Ɋ֗^���鎲
							double u=PART[i].u[d1];//�o�͂Ɋ֗^���鎲
							double w=PART[i].u[d2];//�o�͂Ɋ֗^���鎲

							double y=PART[i].r[d3];//�o�͂Ɋ֗^���Ȃ���
							double v=PART[i].u[d3];//�o�͂Ɋ֗^���Ȃ���

							double r=sqrt(x*x+y*y);//���_����̋���
			
							if(r>0)
							{
								double SIN,COS;
								if(x>0)
								{
									SIN=-y/r;
									COS=x/r;
								}
								if(x<0)
								{
									SIN=y/r;
									COS=-x/r;
								}
								double x2=COS*x-SIN*y;//��]��̍��W�@�~�����̂�x2�̂݁By2�͂���Ȃ�
								double u2=COS*u-SIN*v;//��]��̑��x�@�~�����̂�u2�̂݁Bv2�͂���Ȃ�
								Micro_avs.make_list(x2,z,u2,w);
						//		vec<<x2<<" "<<z<<" "<<u2*times<<" "<<w*times<<endl;
								//if(SIN*x+COS*y!=0) cout<<SIN*x+COS*y<<endl;
								if(x2>xmax) xmax=x2;
								if(z>ymax) ymax=z;
							}
							else Micro_avs.make_list(x,z,u,w);
						}
					}
				}
			}
		}
		Micro_avs.Output_vector_MicroAVS("Speed",CON.get_current_step());
/*		xmax+=4*le;//�}����o���ʒu��ی��������ď����΂߂Ɉړ�
		ymax+=4*le;
//		if(CON.get_legend_speed()>0) vec<<xmax<<" "<<ymax<<" "<<CON.get_legend_speed()*times<<" "<<0*times<<endl;//�Ō�ɖ}��o��
		vec.close();*/
	}
	/////////////////////////////

	/////////�d�S�ɑ΂��鑊�Α��x�o��
	if(CON.get_relative_speed()==ON)
	{
		double U=0;								//���ϑ��x
		double V=0;
		ofstream vec2("relative_speed.dat");

		///���ϑ��x�v�Z
		if(d==2) for(int i=0;i<fluid_number;i++) {U+=PART[i].u[A_X]; V+=PART[i].u[A_Y];}
		else if(d==3) for(int i=0;i<fluid_number;i++) {U+=PART[i].u[d1]; V+=PART[i].u[d2];}//d1,d2�͂��łɂ��Ƃ܂��Ă���
		if(fluid_number>0) {U/=fluid_number; V/=fluid_number;}
		/////////////////////

		if(d==2)
		{
			for(int i=0;i<fluid_number;i++)
			{
				double x=PART[i].r[A_X];
				double y=PART[i].r[A_Y];
				double u=PART[i].u[A_X]-U;
				double v=PART[i].u[A_Y]-V;
				vec2<<x<<" "<<y<<" "<<u*times<<" "<<v*times<<endl;
			}
		}
		else if(d==3)
		{
			for(int i=0;i<fluid_number;i++)
			{
				if(PART[i].r[face]>face_p-le && PART[i].r[face]<face_p+le)
				{
					double x=PART[i].r[d1];
					double z=PART[i].r[d2];
					double u=PART[i].u[d1]-U;
					double w=PART[i].u[d2]-V;
					vec2<<x<<" "<<z<<" "<<u*times<<" "<<w*times<<endl;
				}
			}
		}
		vec2.close();
	}/////////////////////


	if(CON.get_flat_speed_plot()==ON && d==3) //���������̑��x�o��
	{
		ofstream flat("flat_speed.dat");		
		double face_p2=CON.get_flat_speed_p();
		//for(int i=startID;i<NUM;i++)
		for(int i=0;i<fluid_number;i++)
		{
			if(PART[i].r[A_Z]>face_p2-le && PART[i].r[A_Z]<face_p2+le) flat<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].u[A_X]*times<<" "<<PART[i].u[A_Y]*times<<endl;
		}
		flat.close();
	}

	if(CON.get_speed_AVS()==ON && d==3)
	{
		num=0;
		for(int i=startID;i<NUM;i++) if(PART[i].r[face]<face_p) num++;
		
		ofstream fout2("speed_dist.fld");
		fout2 << "# AVS field file" << endl;
		fout2 << "ndim=1" << endl;
		fout2 << "dim1=" << num <<endl;
		fout2 << "nspace=3" << endl;
		fout2 << "veclen=3" << endl;
		fout2 << "data=float" << endl;
		fout2 << "field=irregular" << endl;
		fout2 << "label=e-x e-y e-z" << endl << endl;
		fout2 << "variable 1 file=./speedvec filetype=ascii skip=1 offset=0 stride=6" << endl;
		fout2 << "variable 2 file=./speedvec filetype=ascii skip=1 offset=1 stride=6" << endl;
		fout2 << "variable 3 file=./speedvec filetype=ascii skip=1 offset=2 stride=6" << endl;
		fout2 << "coord    1 file=./speedvec filetype=ascii skip=1 offset=3 stride=6" << endl;
		fout2 << "coord    2 file=./speedvec filetype=ascii skip=1 offset=4 stride=6" << endl;
		fout2 << "coord    3 file=./speedvec filetype=ascii skip=1 offset=5 stride=6" << endl;
		fout2.close();

		ofstream fout("speedvec");
		fout<<"e-x e-y e-z x y z"<<endl;
		for(int i=startID;i<NUM;i++)
		{
			if(PART[i].r[face]<face_p)
			{
				fout<<PART[i].u[A_X]*times<<" "<<PART[i].u[A_Y]*times<<" "<<PART[i].u[A_Z]*times<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
			}
		}
		fout.close();
	}
}

//microAVS�p���q����o�͊֐��i400�s�j
void particle_movie_AVS(mpsconfig &CON,vector<mpselastic> &PART, elastic &ELAST, int fluid_number,int particle_number,int t,double T)
{
	//t:�^�C���X�e�b�v�@T:��������
	double le=CON.get_distancebp();
	double times=CON.get_P_size_AVS();					//�o�͂��闱�q�̃T�C�Y(le�̉��{��)
	bool cut=OFF;	//�f�ʕ\���t���O

	if(t==1)
	{
		ofstream fout("particle_movie.mgf");

		fout<<"# Micro AVS Geom:2.00"<<endl;
		fout<<CON.get_step()/CON.get_interval()+1<<endl;//microAVS�ɏo�͂��鑍�X�e�b�v���B�t�@�C���o�͂�CON.get_interval()���1��ƍŏ��ɍs���B
		
		fout.close();
	}

	ofstream avs("particle_movie.mgf",ios :: app);

	if(CON.get_interval()==1) 
	{
		avs<<"step"<<t<<endl;
	}

	else
		avs<<"step"<<t/CON.get_interval()+1<<endl;
		avs<<"sphere"<<endl;
		avs<<"time="<<T<<endl;
		avs<<"color"<<endl;


	double red,green,blue;	//���q�̐F��\������3���F
	
	////////////////////////////���q�̓���������\��//////////////////////////				
	if(CON.get_AVS()==0)	
	{	
		double le=CON.get_distancebp();
		double mass=CON.get_particle_mass();


	////////////////////////////�\�����闱�q�����v�Z/////////////////////////
		int num=0;//�\�����闱�q��
		int num2=0;

		if(CON.get_dimension()==2) num=particle_number;//2�����ł͑S���q��\��
		else if(CON.get_dimension()==3){
			for(int i=0;i<PART.size();i++){
		//		if(PART[i].surface==ON) num++;
				num++;
			}
		}
		avs<<num<<endl;
		/////////////////////////////////

		if(CON.get_dimension()==2)
		{
			for(int i=0;i<PART.size();i++)
			{	
				double width=CON.get_max_pressure()-CON.get_min_pressure();//���x�̕�
				double P=0.0;
				for(int D=0;D<3;D++) P+=PART[i].get_pressure_accel(D)*PART[i].get_pressure_accel(D);

				P=sqrt(P);

				double level=(P-CON.get_min_pressure())/width;
				red=level*level; //�Ȃ���悵�Ă���H
				green=-4*level*(level-1);//�^�񒆂�1�̕�����
				blue=1-red;

				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
				
				avs<<CON.get_distancebp()*times<<" ";//���q�̑傫���o��
	
				avs<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
			}
		}
		else if(CON.get_dimension()==3)
		{

			for(int i=0;i<PART.size();i++)
			{
				if(PART[i].type==MAGELAST || PART[i].type==MAGELAST2)
				{
					red=0;
					green=1;//�^�񒆂�1�̕�����
					blue=0;
				}
				else if(PART[i].type==ELASTIC)
				{
					red=0;
					green=0;
					blue=1;//�^�񒆂�1�̕�����
				}
				else if(PART[i].type==WALL)
				{
					red=1;//1
					green=0;//0
					blue=0;//�^�񒆂�1�̕����� 0
				}//*/
				else if(PART[i].type==TERMINAL1 || PART[i].type==TERMINAL2)
				{
					red=1;
					green=1;
					blue=0;//�^�񒆂�1�̕�����
				}
				else if(PART[i].type==HYPERELAST)
				{
					red=0;
					green=1;
					blue=1;
				}

				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
				
				avs<<CON.get_distancebp()*times<<" ";//���q�̑傫���o��
	
				avs<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
			}
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////���q�̈��͂��R���^�[�\��//////////////////////////////////
	else if(CON.get_AVS()==1)
	{
		int num=0;//�\�����闱�q��
		num=particle_number;
		
		avs<<num<<endl;

		///////////////���͂̍ő�l�ƍŏ��l�����Ƃ߂�//////////
		double maxP=0;//
		double minP=0;
		for(int i=0;i<particle_number;i++)
		{
			if(PART[i].P>maxP) maxP=PART[i].P;
			if(PART[i].P<minP) minP=PART[i].P;
		}
		/////maxP,minP�����܂����B
		//////////////////////////////////////////////////////

		double width=maxP-minP;
		if(width<0.0001) width=0.0001;
			
		if(CON.get_dimension()==2)
		{
			for(int i=0;i<particle_number;i++)
			{
				double level=(PART[i].P-minP)/width;
				red=level*level;
				green=-4*level*(level-1);//�^�񒆂�1�̕�����
				blue=1-red;
				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
				
				avs<<CON.get_distancebp()*times<<" ";//���q�̑傫���o��
	
				avs<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
			}
		}
		else if(CON.get_dimension()==3)
		{
			for(int i=0;i<particle_number;i++)//2013-02-26 fluid_number��particle_number
			{
				double level=(PART[i].P-minP)/width;
				red=level*level;
				green=-4*level*(level-1);//�^��(level==1/2)��1�̕�����
				blue=1-red;
				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
				
				avs<<CON.get_distancebp()*times<<" ";//���q�̑傫���o��
	
				avs<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
			}
		}
		
	}
	else if(CON.get_AVS()==2)	//���q�̉��x���R���^�[�\��
	{
		//////////////////////////////////�\�����闱�q�����v�Z
		int num=0;//�\�����闱�q��
		if(CON.get_dimension()==2) num=particle_number;//2�����ł͑S���q��\��
		else if(CON.get_dimension()==3) num=fluid_number;//3�����ł͗��̗��q�݂̂�\��
		
		avs<<num<<endl;
		/////////////////////////////////

		double le=CON.get_distancebp();
		double mass=CON.get_particle_mass();
		double T;//���qi�̉��x
		double width=CON.get_maxT()-CON.get_minT();//���x�̕�

		if(CON.get_dimension()==2)
		{
			for(int i=0;i<particle_number;i++)
			{
				double hs0=mass*CON.get_Cp()*CON.get_MP();//�Z���J�n�_�̃G���^���s�[
				double hs1=hs0+CON.get_latent_H()*mass;//�Z���I���_�̃G���^���s�[
				if(PART[i].h<hs0) T=PART[i].h/mass/CON.get_Cp();
				else if(hs0<=PART[i].h && PART[i].h<=hs1) T=CON.get_MP();
				else if(hs1<PART[i].h) T=CON.get_MP()+(PART[i].h-hs1)/mass/CON.get_Cp();
				double level=(T-CON.get_minT())/width;
				red=level*level;
				green=-4*level*(level-1);//�^�񒆂�1�̕�����
				blue=1-red;

				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
				
				avs<<CON.get_distancebp()*times<<" ";//���q�̑傫���o��
	
				avs<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
			}
		}
		else if(CON.get_dimension()==3)
		{
			for(int i=0;i<particle_number;i++)
			{
				double hs0=mass*CON.get_Cp()*CON.get_MP();//�Z���J�n�_�̃G���^���s�[
				double hs1=hs0+CON.get_latent_H()*mass;//�Z���I���_�̃G���^���s�[
				if(PART[i].h<hs0) T=PART[i].h/mass/CON.get_Cp();
				else if(hs0<=PART[i].h && PART[i].h<=hs1) T=CON.get_MP();
				else if(hs1<PART[i].h) T=CON.get_MP()+(PART[i].h-hs1)/mass/CON.get_Cp();
				double level=(T-CON.get_minT())/width;
				red=level*level;
				green=-4*level*(level-1);//�^�񒆂�1�̕�����
				blue=1-red;

				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
				
				avs<<CON.get_distancebp()*times<<" ";//���q�̑傫���o��
	
				avs<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
			}
		}
	}
	else if(CON.get_AVS()==3)//�Ǘ��q�͔�\��
	{
		//////////////////////////////////�\�����闱�q�����v�Z
		int num=fluid_number;//�\�����闱�q��
		avs<<num<<endl;
		/////////////////////////////////

		////���q�o��
		if(CON.get_dimension()==2) cout<<"error in AVS() 2D�͔�Ή�"<<endl;      
		else if(CON.get_dimension()==3)
		{
			for(int i=0;i<fluid_number;i++)
			{
				if((PART[i].type==FLUID || PART[i].type==ELASTIC ) && PART[i].surface==OFF) {red=0;green=0;blue=1;}
				else if((PART[i].type==FLUID || PART[i].type==ELASTIC) && PART[i].surface==ON) {red=0;green=0.5;blue=0.5;}
	        	else {red=0.5;green=0.5;blue=0;}//�Ǘ��q
			    
				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
			
				avs<<CON.get_distancebp()*times<<" ";//���q�̑傫���o��

				avs<<red<<" "<<green<<" "<<blue<<endl;//�F�o��        
			}
		}
	}
	else if(CON.get_AVS()==4)//���̕\�ʗ��q�̂�
	{
		//////////////////////////////////�\�����闱�q�����v�Z
		int num=0;//�\�����闱�q��
		if(CON.get_dimension()==2) num=particle_number;//2�����ł͑S���q��\��
		else if(CON.get_dimension()==3) for(int i=0;i<particle_number;i++) if((PART[i].type==FLUID || PART[i].type==ELASTIC) && PART[i].surface==ON) num++;//3�����̏ꍇ�A�������͕̂\�����Ȃ�	
		avs<<num<<endl;
		/////////////////////////////////

		////���q�o��
		if(CON.get_dimension()==2)
		{
			for(int i=0;i<particle_number;i++)
			{
				if((PART[i].type==FLUID || PART[i].type==ELASTIC) && PART[i].surface==OFF) {red=0;green=0;blue=1;}
				else if((PART[i].type==FLUID || PART[i].type==ELASTIC) && PART[i].surface==ON) {red=0;green=0.5;blue=0.5;}
	        	else {red=0.5;green=0.5;blue=0;}//�Ǘ��q
			   
				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
			
				avs<<CON.get_distancebp()/2<<" ";//���q�̑傫���o��

				avs<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
			}
			
		}
		else if(CON.get_dimension()==3)
		{
			for(int i=0;i<fluid_number;i++)
			{
				if(PART[i].surface==ON) 
				{
					if((PART[i].type==FLUID || PART[i].type==ELASTIC) && PART[i].surface==OFF) {red=0;green=0;blue=1;}
					else if((PART[i].type==FLUID || PART[i].type==ELASTIC) && PART[i].surface==ON) {red=0;green=0.5;blue=0.5;}
	        		else {red=0.5;green=0.5;blue=0;}//�Ǘ��q
			    
					avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
			
					avs<<CON.get_distancebp()/2<<" ";//���q�̑傫���o��

					avs<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
				}        
			}
		}
	}
	if(CON.get_AVS()==5)	//�ǂ�����\��.�����̊m�F�p
	{

	}
	else if(CON.get_AVS()==6){}	//����̗��q�̓���������\���B���q�̑I���̓v���O�����𒼐ڕύX���邵���Ȃ�
	


	avs.close();
		////////////////////

	if(CON.get_dimension()==3 && CON.get_AVS()==0)//�f�ʕ\���t���O
	{
		if(CON.get_cut_x()==ON)
		{	
			if(t==1)
			{
				ofstream fout("half_x_particle_movie.mgf");

				fout<<"# Micro AVS Geom:2.00"<<endl;
				fout<<CON.get_step()/CON.get_interval()+1<<endl;//microAVS�ɏo�͂��鑍�X�e�b�v���B�t�@�C���o�͂�CON.get_interval()���1��ƍŏ��ɍs���B
				fout.close();
			}

			ofstream avs2("half_x_particle_movie.mgf",ios :: app);
	
			if(CON.get_interval()==1) avs2<<"step"<<t<<endl;

			else
				avs2<<"step"<<t/CON.get_interval()+1<<endl;
				avs2<<"sphere"<<endl;
				avs2<<"time="<<T<<endl;
				avs2<<"color"<<endl;
		

			double red,green,blue;	//���q�̐F��\������3���F
	
			////////////////////////////���q�̓���������\��//////////////////////////				

			double le=CON.get_distancebp();
			double mass=CON.get_particle_mass();
			double width=CON.get_max_pressure()-CON.get_min_pressure();//���x�̕�

		////////////////////////////�\�����闱�q�����v�Z/////////////////////////
			int num=0;//�\�����闱�q��


			for(int i=0;i<PART.size();i++)
			{
				if(PART[i].r[A_Y]>=-0.00025 && PART[i].surface==ON)
				{
					num++;
				}
			}

			avs2<<num<<endl;
			/////////////////////////////////

			double P=0.0;
			for(int i=0;i<PART.size();i++)
			{
				if(PART[i].r[A_Y]>=-0.00025 && PART[i].surface==ON){
				//���͂̌v�Z

				for(int D=0;D<3;D++) P+=PART[i].get_pressure_accel(D)*PART[i].get_pressure_accel(D);

				P=sqrt(P);

				double level=(P-CON.get_min_pressure())/width;
				if(PART[i].type==MAGELAST || PART[i].type==MAGELAST2){
				red=0;
				green=1;//�^�񒆂�1�̕�����
				blue=0;
				}
				else if(PART[i].type==ELASTIC){
				red=0;
				green=0;
				blue=1;//�^�񒆂�1�̕�����
				}
				else if(PART[i].type==WALL ){
				red=1;
				green=0;
				blue=0;//�^�񒆂�1�̕�����
				}
				else if(PART[i].type==TERMINAL1 || PART[i].type==TERMINAL2){
				red=1;
				green=1;
				blue=0;//�^�񒆂�1�̕�����
				}
				else if(PART[i].type==HYPERELAST)
				{
					red=0;
					green=1;
					blue=1;
				}

				avs2<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
		
				avs2<<CON.get_distancebp()*times<<" ";//���q�̑傫���o��

				avs2<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
				}
			}
			avs2.close();
		}

		if(CON.get_cut_y()==ON)  	//�f�ʕ\���t���O
		{	
			if(t==1)
			{
				ofstream fout("half_y_particle_movie.mgf");

				fout<<"# Micro AVS Geom:2.00"<<endl;
				fout<<CON.get_step()/CON.get_interval()+1<<endl;//microAVS�ɏo�͂��鑍�X�e�b�v���B�t�@�C���o�͂�CON.get_interval()���1��ƍŏ��ɍs���B
				fout.close();
			}

			ofstream avs3("half_y_particle_movie.mgf",ios :: app);
	
			if(CON.get_interval()==1) avs3<<"step"<<t<<endl;

			else
				avs3<<"step"<<t/CON.get_interval()+1<<endl;
				avs3<<"sphere"<<endl;
				avs3<<"time="<<T<<endl;
				avs3<<"color"<<endl;
		

			double red,green,blue;	//���q�̐F��\������3���F
	
			////////////////////////////���q�̓���������\��//////////////////////////				

			double le=CON.get_distancebp();
			double mass=CON.get_particle_mass();
			double width=CON.get_max_pressure()-CON.get_min_pressure();//���x�̕�

		////////////////////////////�\�����闱�q�����v�Z/////////////////////////
			int num=0;//�\�����闱�q��


			for(int i=0;i<PART.size();i++)
			{
				if(PART[i].r[A_X]>=-0.00025 && PART[i].surface==ON)
				{
					num++;
				}
			}

			avs3<<num<<endl;
			/////////////////////////////////

			double P=0.0;
			for(int i=0;i<PART.size();i++)
			{
				if(PART[i].r[A_X]>=-0.00025 && PART[i].surface==ON){
				//���͂̌v�Z

				for(int D=0;D<3;D++) P+=PART[i].get_pressure_accel(D)*PART[i].get_pressure_accel(D);

				P=sqrt(P);

				double level=(P-CON.get_min_pressure())/width;
				if(PART[i].type==MAGELAST || PART[i].type==MAGELAST2){
				red=0;
				green=1;//�^�񒆂�1�̕�����
				blue=0;
				}
				else if(PART[i].type==ELASTIC){
				red=0;
				green=0;
				blue=1;//�^�񒆂�1�̕�����
				}
				else if(PART[i].type==WALL ){
				red=1;
				green=0;
				blue=0;//�^�񒆂�1�̕�����
				}
				else if(PART[i].type==TERMINAL1 || PART[i].type==TERMINAL2){
				red=1;
				green=1;
				blue=0;//�^�񒆂�1�̕�����
				}
				else if(PART[i].type==HYPERELAST)
				{
					red=0;
					green=1;
					blue=1;
				}

				avs3<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
		
				avs3<<CON.get_distancebp()*times<<" ";//���q�̑傫���o��

				avs3<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
				}
			}
			avs3.close();
		}
	}
}





//���͂ȂǗ��q�̎������R���^�[�}�ŕ\������֐�
void particle_movie_AVS2(mpsconfig &CON,int t,vector<mpselastic> &PART,double TIME,int fluid_number, int particle_number, double **F)
{
	//�Q�l�ɂ��Ă��鏑����microAVS�̃w���v�ł��Ȃ��̃f�[�^�́H���u��\���i�q�^�f�[�^�i�A�X�L�[�j�̏����v
	
	double le=CON.get_distancebp();
	int STEP=CON.get_step()/CON.get_EM_interval()+1;		//�o�͂��鑍�X�e�b�v��
	int *input=new int[particle_number];			//input[]=ON�Ȃ�o��

	if(t==1) 
	{
		ofstream fp("particle_movie2.inp");			
		fp<<STEP<<endl;						//microAVS�ɏo�͂��鑍�ï�ߐ��B�t�@�C���o�͂�CON->get_interval()���1��ƍŏ��ɍs���B
		fp<<"data_geom"<<endl;
		fp.close();
	}


	//main�t�@�C����������
	ofstream fp("particle_movie2.inp",ios :: app);
	fp<<"step"<<t/CON.get_EM_interval()+1<<" TIME="<<TIME<<endl;

	//�o�͗��q���Z�o
	int num=0;//�\�����闱�q��
	if(CON.get_dimension()==2)
	{
		num=particle_number;//2�����ł͑S���q��\��
		for(int i=0;i<particle_number;i++) input[i]=ON;
	}
	else if(CON.get_dimension()==3)
	{
		for(int i=0;i<particle_number;i++)
		{
			if(PART[i].type!=WALL || PART[i].type!=ELASTIC)
			{
				input[i]=ON;
				num++;//3�����̏ꍇ�A�������͕̂\�����Ȃ�
			}
			else input[i]=OFF;
		}
	}

	fp<<num<<" "<<num<<endl;	//�ߓ_���Ɨv�f���o�� ���̏ꍇ�͗����Ƃ��ߓ_��������Ă���

	//�ߓ_�ԍ��Ƃ��̍��W�̏o�� 
	int count=0;
	for(int i=0;i<particle_number;i++)
	{
//		if(input[i]==ON)
		{
			 fp<<count<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
			 count++;
		}
	}

	//�v�f�ԍ��Ɨv�f�`��̎�ށA�����ėv�f���\������ߓ_�ԍ��o��
	for(int i=0;i<num;i++)
	{
		fp<<i<<"  0 pt "<<i<<endl;
	}

	//fp<<"2 3"<<endl;//�ߓ_�̏��ʂ�2�ŁA�v�f�̏��ʂ�3�Ƃ������ƁB
	fp<<"5 0"<<endl;//�ߓ_�̏��ʂ�8�ŁA�v�f�̏��ʂ�0�Ƃ������ƁB
	fp<<"5 1 1 1 1 1"<<endl;	//���̍s�̏ڍׂ̓w���v���Q��
	fp<<"surface,"<<endl;
	fp<<"type,\n";
	fp<<"MF_x,"<<endl;
	fp<<"MF_y,"<<endl;
	fp<<"MF_z,"<<endl;


	//�e�ߓ_�̏��l����
	count=0;
	for(int i=0;i<particle_number;i++)
	{
//		if(input[i]==ON)
		{
			fp<<count<<" "<<PART[i].surface<<" "<<PART[i].type<<" "<<F[A_X][i]<<" "<<F[A_Y][i]<<" "<<F[A_Z][i]<<"\n";
			count++;
		}
	}

	
	fp.close();

	delete [] input;

	
}