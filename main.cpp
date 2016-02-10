
#include "stdafx.h"		//��v�ȃw�b�_�[�t�@�C���͂܂Ƃ߂Ă��̂Ȃ��B
//only main function

#include "Model.h"
#include "Check.h"
#include "Rigidbody.h"
#include "Micro_AVS.h"
void file_initialization();
void Make_STL();//STL�p�������킹�֐�
//////////////////////
int _tmain(int argc, _TCHAR* argv[])
{
	mpsconfig CON;				//��͏����i�[�N���X
	//CON.Set_initial_config();	//��͏�������


	clock_t t1=clock();	//�o�ߎ��Ԃ�b�ŕ\������ɂ́ACLOCKS_PER_SEC�Ŋ���
	
    int particle_number=0;	//�S���q��
	int fluid_number=0;		//���̗��q���i�e���̂̏ꍇ�͒e���̗��q���j
	int out;				//fluid_number<=i<out��INWALL�Q�ŁAout<=i��OUTWALL�Q�ɂȂ�B �Ȃ��ĂȂ��Ȃ玩���ł��Ƃŕ��ёւ���
	int count;				//�J�E���g�p�ϐ�
	int order_sw=OFF;		//���q����ёւ���K�v�����邩�ǂ����̃t���O
	int t=0;				//�X�e�b�v��
	int dimension=CON.get_dimension();
	int model_set=CON.get_model_set_way();
    double dt=CON.get_dt();
    int STEP=CON.get_step();
    double g[DIMENSION]={0,0,0};
    if(dimension==2) g[A_Y]=CON.get_g();
    if(dimension==3) g[A_Z]=CON.get_g();



    //MPS�ɂ�����萔�v�Z
	double n0=initial_pnd(CON.get_re(), dimension, model_set);		//�������q���xn0
    double n0_4=initial_pnd(CON.get_re4(), dimension, model_set);	//freeon�p�������q���x
    double TIME=0; //��͎���
	double Umax=0; //�ő呬�x
	double mindis; //CFL�̍Œᗱ�q�ԋ��� minimum distance

	cout<<"�������q�����x n0= "<<setprecision(10)<<n0<<endl;
	
	//FEM�p��class�쐬
	vector<point3D> NODE_FEM3D;
	vector<element3D> ELEM_FEM3D;
	int node_FEM3D=0;
	int nelm_FEM3D=0;	//�S�ߓ_��,�v�f��

	//�������q�z�u�������݁@restart=ON�̏ꍇ�͗��q���ǂݍ���
	initial_model_input(&CON, &particle_number, &TIME);
//	Model model;
//	particle_number=model.Model_set();

	//INDEX�֌W
    int *INDEX=new int[CON.get_number_of_mesh()];	//�e�i�q�Ɋ܂܂�闱�q�����i�[(�i�q�̐������z�񂪕K�v)
    cout<<"X_mesh="<<CON.get_X_mesh()<<" Y_mesh="<<CON.get_Y_mesh()<<" Z_mesh="<<CON.get_Z_mesh()<<endl;
    cout<<"number_of_mesh="<<CON.get_number_of_mesh()<<endl;
	Check check;
	check.Out_put_config();
	//vector�Ŋm�ۂ���ƃA���S���Y�������p�ł���
	//�������E�E�Epush_back()�ł̓R���X�g���N�^�͓����Ȃ��I�I�i�R�s�[�R���X�g���N�^�������I�I�jang�͕s��I�I
	//����ȍ~��particle_number=PART.size()�ł悢
	mpselastic PART0;
	vector<mpselastic> PART3;
	vector<mpselastic> PART2;
	vector<mpselastic> PART1;	//���q�z���tyoe���ɕ��בւ���B�ꎞ�ۊ�
	vector<mpselastic> PART;
	Rigidbody rigids0;
	vector<Rigidbody> rigids;//����
	cout<<"�x�N�g���쐬����"<<endl;

	for(int i=0;i<2;i++){
		rigids.push_back(rigids0);
	}
	cout<<"1"<<endl;
//	PART.reserve(20000);
	for(int i=0;i<particle_number;i++){
		PART1.push_back(PART0);
		PART.push_back(PART0);
	}

	////////////STL/////////
//	Make_STL();
	///////////////////////

	//���q�f�[�^���t�@�C������ǂݎ��
	//PART�̏����l�̓t�@�C������ǂݍ��ށE�E�E�R���X�g���N�^�͓����Ȃ��I�I
	//�����ŗ��q�ԍ��̐��������Ă��� MRE���V���R�[�������̑��̏���
	input_particle_data(&CON, PART1,PART, 1);	//�Ō�̈�����1��n������initial�f�[�^��ǂݎ��A����ȊO�Ȃ�O�X�e�b�v�f�[�^��ǂݎ��
	Micro_AVS avs;
	for(int i=0;i<PART.size();i++){
		avs.make_list(PART[i].r[A_X],PART[i].r[A_Y],PART[i].r[A_Z],0,0,0);
	}
	avs.Output_mgf_MicroAVS("check_particle",1);
/*	///////////////z���W�����L�����Ă���/////////�ő�ړ������p
	for(int i=0;i<particle_number;i++){
		PART1[i].r[A_Z]=PART[i].r[A_Z];
	}
	/////////////////////////////////////////////*/
	//�e���q�����J�E���g or ���ёւ��D���قƂ�ǈӖ����Ȃ��D���s���Ȃ��Ă�OK(2012/02/21)
	calc_numbers_of_particles_and_change_the_order(&CON, PART, &fluid_number,&out, &order_sw);

	
	//���̃N���X�ɗ��q���i�[
	for(int i=0;i<particle_number;i++){
		if(PART[i].type==TERMINAL1){
			PART2.push_back(PART[i]);
			}		
		if(PART[i].type==TERMINAL2){
			PART3.push_back(PART[i]);
			}
	}
	rigids[0].Get_initial_particle(PART2);
	rigids[1].Get_initial_particle(PART3);//*/

	PART2.clear();
	PART3.clear();
	//�������q�z�u�ŉ�͗̈�O�ɗ��q���Ȃ����`�F�b�N
	check_initial_position(&CON, PART);

	//FEM3D
	double **F=new double*[DIMENSION];//F[D][i]�ƂȂ��Ă��邪�����OK??
	for(int D=0;D<DIMENSION;D++){
		F[D]=new double [(unsigned)PART.size()]; //�e���q�ɓ����d����//particle_number��initial_model_input()�ŋ��܂�D
		for(int i=0;i<PART.size();i++)	F[D][i]=0.0; //������
	}

	//�z��͂̑O��reloadINDEX	
	//���qID�X�V�E���q�����x�X�V
	reload_INDEX(CON,PART, INDEX);//�i�q���̗��q���X�V
	count=0;
	int **MESH0=new int *[CON.get_number_of_mesh()];
	for(int i=0;i<CON.get_number_of_mesh();i++)
	{       
		count+=INDEX[i];
		MESH0[i]=new int [INDEX[i]];
	}
	reload_INDEX2(&CON, PART, MESH0);

	//�\�ʔ���i�e���v�Z�̏ꍇ�ŏ���������OK�I�j
	freeon(CON, PART, n0_4, INDEX, MESH0, &mindis, t);

	for(int i=0;i<CON.get_number_of_mesh();i++) delete [] MESH0[i];
	delete [] MESH0;

	//�������������ɍs���̂ł��̈ʒu(ePND[i]=PART[i].PND�Ȃ̂ł��̈ʒu�łȂ����)
	elastic ELAST(PART);
	for(int i=0;i<PART.size();i++) PART[i].initialize_particles(ELAST, t);
	ELAST.set_ground_position(CON.get_ground_position());
	//�v���v���Z�X�I��

	//file initialization
	file_initialization();
	/////////////////////
	ofstream pt("PART_model.dat");
	for(int i=0;i<PART.size();i++)
	{
		if(PART[i].r[A_X]>(0-CON.get_distancebp()/2) && PART[i].r[A_X]>(0+CON.get_distancebp()/2))
		{
			pt<<PART[i].r[A_Y]<<PART[i].r[A_Z]<<endl;
		}
	}

	pt.close();
	double L=0;
	double W=100;
	double P=0;
	int m=0;
	int wait=0;	//�҂��X�e�b�v
	double limP=0;
	bool ff=0;

	//�\�ʂ����\��
	if(bool cat=ON)
	{
		for(int i=0 ;i<PART.size();i++){
			PART[i].surface=OFF;	//�\�ʏ�����
			if(PART[i].PND<=18){
				PART[i].surface=ON;
			}
		}
	}

	//�v�Z�X�^�[�g
	for(t=1;t<=STEP;t++)
	{	
		cout<<"�z��� start:step="<<t<<", ���q��="<<particle_number<<", dt="<<dt<<endl;

		//�z��͂̑O��reloadINDEX�i���q�͈ړ�����̂Ŗ��X�e�b�v���s����K�v������j
		reload_INDEX(CON, PART, INDEX);	//�i�q���̗��q���X�V
		cout<<"reload_INDEX_OK"<<endl;
		//MESH��mpsconfig�̃����o�֐��ɂ���̂��ǂ��Bnew/delete�댯�Ȃ̂�shared_ptr��vector���g���ׂ�
		int **MESH = new int *[CON.get_number_of_mesh()];	//get_number_of_mesh() {return (int)((maxX-minX)/(distancebp*dx)*(maxY-minY)/(distancebp*dx)*(maxZ-minZ)/(distancebp*dx)+0.001);}//�i�q���FX_mesh*Y_mesh*Z_mesh
		count=0;
		for(int i=0;i<CON.get_number_of_mesh();i++)
		{       
			count+=INDEX[i];
			MESH[i]=new int [INDEX[i]];
		}
		if(count>PART.size()) cout<<"INDEX error ���q������?"<<endl;
		if(count<PART.size()) cout<<"INDEX error ���q������?"<<endl;
		reload_INDEX2(&CON, PART, MESH);

		//�z���
		unsigned timeA=GetTickCount();

		//�e�����a���̗��q���J�E���g�{mindis�̌v�Z
		//�e���v�Z�p�̊֐����g��
		cout<<"freeon_start"<<endl;
		freeon(ELAST, PART, n0_4, INDEX, MESH, &mindis, t); //�\�ʔ���

		cout<<"�z��͑O�̗��q�ˑ��֌W�v�Z�I���@�|�|time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

		//CFL�����E�g�U�������ɂ��dt�̌���
		//�V���v���e�B�b�N�X�L�[�����g���ꍇ�͎��ԉ𑜓x�ς͓K�p���Ă͂Ȃ�Ȃ�
//		if(OFF==ELAST.get_symp_flag()) courant_elastic(&CON, PART, fluid_number, t, &dt, mindis, Umax, g);
//		courant_elastic(&CON, PART, fluid_number, t, &dt, mindis, Umax, g);
		//�e���̌v�Z�����������Ă���dt��ύX����

		//model1����1311�Amodel7����731�Amodel11����1781��ICCG�@�Ō덷�̌v�Z�l�����U���A��͂��j�]����

		if(CON.get_model_number()==1)
		{
			if((CON.get_avoid_step()==0 || CON.get_current_step()!=CON.get_avoid_step())
				&& (CON.get_avoid_step2()==0 || CON.get_current_step()!=CON.get_avoid_step2())
				&& (CON.get_avoid_step3()==0 || CON.get_current_step()!=CON.get_avoid_step3())
				&& (CON.get_avoid_step4()==0 || CON.get_current_step()!=CON.get_avoid_step4())
				&& (CON.get_avoid_step5()==0 || CON.get_current_step()!=CON.get_avoid_step5())
				&& (CON.get_avoid_step6()==0 || CON.get_current_step()!=CON.get_avoid_step6())
				&& (CON.get_avoid_step7()==0 || CON.get_current_step()!=CON.get_avoid_step7()))
			{
				if(CON.get_FEM_flag()==true && t>wait)
				{

					//ELAST����CON�̊֐��Ăяo����̂œ�n���v��Ȃ��E�E�E
					if(ELAST.get_FEM_switch()==true && CON.get_mesh_input()!=2)
					{
						CON.change_step_size(); //����Ȋ֐��͍��������������Ȃ̂ŗv��Ȃ�
						ELAST.change_step_size();
						if(t==1 || (t-1)%CON.get_EM_interval()==0)
						{
							//����f���[�j����
							FEM3D(CON, PART, F, &node_FEM3D, &nelm_FEM3D, NODE_FEM3D, ELEM_FEM3D, t, TIME, fluid_number);
						}
					}
					else if(CON.get_mesh_input()==2)
					{
						if(t==1 || (t-1)%CON.get_EM_interval()==0)
						{	
							//TetGen�ɂ�郁�b�V������
							TetGenInterface(CON, PART, F, fluid_number, dt, t, particle_number, n0, TIME);
						}	
					}
				}
			}
		}

		else if(CON.get_model_number()==7)
		{
//			if(CON.get_current_step()!=2620)
			if(CON.get_FEM_flag()==true && t>wait)
			{

				//ELAST����CON�̊֐��Ăяo����̂œ�n���v��Ȃ��E�E�E
				if(ELAST.get_FEM_switch()==true && CON.get_mesh_input()!=2)
				{
					CON.change_step_size(); //����Ȋ֐��͍��������������Ȃ̂ŗv��Ȃ�
					ELAST.change_step_size();
					if(t==1 || (t-1)%CON.get_EM_interval()==0)
					{
						//����f���[�j����
						FEM3D(CON, PART, F, &node_FEM3D, &nelm_FEM3D, NODE_FEM3D, ELEM_FEM3D, t, TIME, fluid_number);
					}
				}
				else if(CON.get_mesh_input()==2)
				{
					if(t==1 || (t-1)%CON.get_EM_interval()==0)
					{	
						//TetGen�ɂ�郁�b�V������
						TetGenInterface(CON, PART, F, fluid_number, dt, t, particle_number, n0, TIME);
					}	
				}
			}
		}

		else if(CON.get_model_number()==11)
		{
			if((CON.get_avoid_step()==0 || CON.get_current_step()!=CON.get_avoid_step())
				&& (CON.get_avoid_step2()==0 || CON.get_current_step()!=CON.get_avoid_step2())
				&& (CON.get_avoid_step3()==0 || CON.get_current_step()!=CON.get_avoid_step3())
				&& (CON.get_avoid_step4()==0 || CON.get_current_step()!=CON.get_avoid_step4())
				&& (CON.get_avoid_step5()==0 || CON.get_current_step()!=CON.get_avoid_step5())
				&& (CON.get_avoid_step6()==0 || CON.get_current_step()!=CON.get_avoid_step6())
				&& (CON.get_avoid_step7()==0 || CON.get_current_step()!=CON.get_avoid_step7()))
			{
				if(CON.get_FEM_flag()==true && t>wait)
				{

					//ELAST����CON�̊֐��Ăяo����̂œ�n���v��Ȃ��E�E�E
					if(ELAST.get_FEM_switch()==true && CON.get_mesh_input()!=2)
					{
						CON.change_step_size(); //����Ȋ֐��͍��������������Ȃ̂ŗv��Ȃ�
						ELAST.change_step_size();
						if(t==1 || (t-1)%CON.get_EM_interval()==0)
						{
							//����f���[�j����
							FEM3D(CON, PART, F, &node_FEM3D, &nelm_FEM3D, NODE_FEM3D, ELEM_FEM3D, t, TIME, fluid_number);
						}
					}
					else if(CON.get_mesh_input()==2)
					{
						if(t==1 || (t-1)%CON.get_EM_interval()==0)
						{	
							//TetGen�ɂ�郁�b�V������
							TetGenInterface(CON, PART, F, fluid_number, dt, t, particle_number, n0, TIME);
						}	
					}
				}
			}
		}

		else if(CON.get_model_number()!=1 && CON.get_model_number()!=7 && CON.get_model_number()!=11)
		{
			if(CON.get_FEM_flag()==true && t>wait)
			{

				//ELAST����CON�̊֐��Ăяo����̂œ�n���v��Ȃ��E�E�E
				if(ELAST.get_FEM_switch()==true && CON.get_mesh_input()!=2)
				{
					CON.change_step_size(); //����Ȋ֐��͍��������������Ȃ̂ŗv��Ȃ�
					ELAST.change_step_size();
					if(t==1 || (t-1)%CON.get_EM_interval()==0)
					{
						//����f���[�j����
						FEM3D(CON, PART, F, &node_FEM3D, &nelm_FEM3D, NODE_FEM3D, ELEM_FEM3D, t, TIME, fluid_number);
					}
				}
				else if(CON.get_mesh_input()==2)
				{
					if(t==1 || (t-1)%CON.get_EM_interval()==0)
					{	
						//TetGen�ɂ�郁�b�V������
						TetGenInterface(CON, PART, F, fluid_number, dt, t, particle_number, n0, TIME);
					}	
				}
			}
		}


			//���q�ړ��v�Z
			calc_elastic(PART, ELAST, t, F);

			cout<<"�z��͏I�� umax="<<sqrt(Umax)<<"  limit U="<<0.2*mindis/dt<<endl;

			//���q���������̂�INDEX�X�V
			for(int i=0;i<CON.get_number_of_mesh();i++) delete [] MESH[i];
			delete [] MESH;

			double umax2=0.0;//�ő呬�x
			for(int i=0;i<PART.size();i++)
			{ 
				double speed=0.0;//���q���x
				for(int D=0;D<DIMENSION;D++) speed+=PART[i].u[D]*PART[i].u[D];
				if(speed>umax2) umax2=speed;
				
			}
			cout<<"umax="<<sqrt(umax2)<<endl;

			if(t>=wait){
			TIME+=dt;///���ԍX�V
			CON.set_current_time(TIME);
			}
			cout<<"��͕�������="<<TIME<<"/ "<<(CON.get_step()*CON.get_dt())<<endl;

			//�|�X�g�����F
			post_processing(CON, PART, ELAST, particle_number, particle_number, dt, Umax, t, TIME,F); //�e�����ʏo�́��N�[�������ɂ��dt����&microAVS�o��
			post_processing3(CON, PART, particle_number, particle_number, t, TIME); //restart�p�t�@�C������
	//		order_sw=check_position(&CON, PART, particle_number, &particle_number); //�̈�O�̗��q������ //�̈�O���q�����m����΁Aorder_sw=ON�ɂȂ�
			clock_t t3=clock();
			cout<<"CPU time="<<(t3-t1)/CLOCKS_PER_SEC<<"[sec]"<<endl;
			ofstream t_log("time_log.dat", ios::app);

			t_log<<"step="<<t<<", time="<<(t3-t1)/CLOCKS_PER_SEC<<"[sec]"<<endl;
			t_log.close();
			//�ŉ��ʂ̈��͕\��
	/*	/////////////////////////////////////////////
		double side=100;
		int side_num=0;
		for(int i=0;i<PART.size();i++)
		{
			if(PART[i].type==ELASTIC || PART[i].type==MAGELAST)
			{
				if(PART[i].r[A_Z]<side) {
					side=PART[i].r[A_Z];
					side_num=i;
				}
			}
		}
		cout<<"���ʂ̈���="<<PART[side_num].P<<" ,ID="<<side_num<<endl;
		//////////////////////////////////////////////*/
			/////////////////�ő�ړ������o��////////////
			if(CON.get_model_number()==7 && t>=wait){
				if(t==wait || t==1){
					L=PART[2753].r[A_Z];
				}
				else if(t%10==0 && t>wait){
				ofstream longz("longZ.dat", ios::app);
				double adis=0;
				 adis=(L-PART[2753].r[A_Z]);
		//		cout<<PART1[95].r[A_Z]<<" "<<PART[95].r[A_Z];
				longz<<TIME<<" "<<adis*1000<<endl;//[mm]
				longz.close();
				}
			}
			else if(CON.get_model_number()==10){
				if(t==1){
					L=PART[CON.Get_length()].r[A_Z];
				}
				else if(t%1000==0){
					ofstream longz("longZ.dat", ios::app);
					double adis=0;
					adis=(L-PART[CON.Get_length()].r[A_Z]);
					longz<<TIME<<" "<<adis*1000<<endl;//[mm]
				longz.close();
				}
			}
			/////////////////////////////////////////////*/
	/*		if(CON.get_model_number()==4){
				int point[4];
				model.Get_point(point[0],point[1],point[2],point[3]);
				if(t==1){
					L=PART[point[1]].r[A_Z]-PART[point[0]].r[A_Z];
					W=PART[point[3]].r[A_Y]-PART[point[2]].r[A_Y];
					ELAST.set_poise_flag(ON);
				}
			
				else{
				double dL=0.0;
				double e=0.0;
				double dW=0.0;
				double de=0.0;
		
				dL=(PART[point[1]].r[A_Z]-PART[point[0]].r[A_Z])-L;
				dW=(PART[point[3]].r[A_Y]-PART[point[2]].r[A_Y])-W;
				if(dL!=0)e=dL/L;
				if(dW!=0)de=fabs(dW/W);
				P=de/e;
				cout<<de<<","<<e<<endl;
				cout<<limP<<","<<P<<endl;
				cout<<floor(limP*1000000)<<","<<floor(P*1000000)<<endl;
				if(floor(limP*1000000)==floor(P*1000000) && floor(P*1000000)>100000){
					if(ff==1){
				ofstream hiz("hizumi-poa.dat", ios::app);
				hiz<<P<<" "<<e*100<<endl;
				hiz.close();
				ofstream yiz("hizumi-yang_stress.dat", ios::app);
				double stress=0;
				double pressure=0;
				double nensei=0;
				int partnum=0;
				for(int i=1452;i<=1572;i++){
					stress+=fabs(PART[i].get_stress_accel(A_Z))*ELAST.get_density();
					partnum++;
				}
			//	stress+=fabs(PART[1512].get_stress_accel(A_Z))*ELAST.get_density();
				stress/=partnum;
				yiz<<stress/(e*100)<<" "<<e*100<<endl;
				yiz.close();
				ELAST.set_poise_flag(ON);
					}
					ff=1;
			}
				else {
					limP=P;
					ELAST.set_poise_flag(OFF);
					ff=0;
				}
			}
			}//*/
		
		
		

			check.Courant_condition(PART);

			cout<<endl;
		
	}
	//���[�v�I��

	for(int D=0;D<DIMENSION;D++){
//		delete [] laplacian[D];
		delete [] F[D];
	}

	delete [] INDEX;

	clock_t t2=clock();
	cout<<"CPU last time="<<(t2-t1)/CLOCKS_PER_SEC<<"[sec]"<<endl;
//	MessageBeep(MB_ICONEXCLAMATION);//��Ƃ̏I����m�点��BEEP�� �}���O�����O�G���[���N����̂ŃR�����g�A�E�g

	return 0;
}

//�d�݊֐��@��`�����������I�I
double kernel(double r,double dis)
{
	return (dis<r) ? r/dis-1: 0;
//	return r/dis-1;
	//return r*r/(dis*dis)-1;
	//return r*r*r/(dis*dis*dis)-1;
	//return (1-dis/r)*(1-dis/r);
	//return 1;
}

//�d�݊֐��Q
double kernel2(double r,double dis,double d)
{
	return r/dis-1;
    //return r*r*r/(dis*dis*dis)-1;
	//return r*r*r*r/(dis*dis*dis*dis)-1;
	//return pow(r,d)/pow(dis,d);
}
//�d�݊֐�3�@�Ǘp
double kernel3(double r,double dis)
{
	double ndis=0;
	ndis=dis-(r*0.1);
	return (dis<r) ? 5*(r/ndis-1): 0;
}

//�������q���x�̌v�Z�i90�s�j
double initial_pnd(double r,int dimension,int calc_type)
{
	int size = (int)(r+1);//�v�Z�̈�
	double dis;//����
	double pnd=0;
	int count=0;
	if(dimension==2)
	{
		if(calc_type==0)				//�����z�u�Ƃ��Đ����i�q���Ƃ����ꍇ
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-size;j<=size;j++)
				{
					dis=sqrt((double)(i*i+j*j));
					if(dis!=0 && dis<=r )
					{
						pnd+=kernel(r,dis);
						count++;
					}			
				}
			}
		}
		if(calc_type==1)				//�����z�u�Ƃ��čז��Z�@���Ƃ����ꍇ
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-2*size;j<=2*size;j++)
				{
					double jj=j*sqrt(3.0)/2;
					double ii=(double)i;
					if(j%2!=0) ii+=0.5;//j����Ȃ�ii��0.5�i�q�������炷
					dis=sqrt(ii*ii+jj*jj);
					if(dis!=0 && dis<=r )
					{
						pnd+=kernel(r,dis);
						count++;
					}			
				}
			}
		}
	}
	else if(dimension==3)
	{
		if(calc_type==0)				//�����z�u�Ƃ��Đ����i�q���Ƃ����ꍇ
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-size;j<=size;j++)
				{
					for(int k=-size;k<=size;k++)
					{
						dis=sqrt((double)(i*i+j*j+k*k));
						if(dis!=0 && dis<=r )
						{
							pnd+=kernel(r,dis);
							count++;
						}
					}			
				}
			}
		}
		if(calc_type==1)				//�����z�u�Ƃ��čז��Z�@���Ƃ����ꍇ
		{
			for(int i=-2*size;i<=2*size;i++)
			{
				for(int j=-2*size;j<=2*size;j++)
				{
					for(int k=-2*size;k<=2*size;k++)
					{
						double ii=(double)i;
						double jj=j*sqrt(3.0)/2;
						double kk=k*sqrt(2.0/3);
						if(j%2!=0) ii+=0.5;//j����Ȃ�ii��0.5�i�q�������炷
						if(k%2!=0) {ii+=0.5; jj+=sqrt(3.0)/6;}//k����Ȃ�ii��jj�����炷
						dis=sqrt(ii*ii+jj*jj+kk*kk);
						if(dis!=0 && dis<=r )
						{
							pnd+=kernel(r,dis);
							count++;
						}
					}			
				}
			}
		}
	}
	cout<<"n0��count="<<count<<endl;
	return pnd;  
}

//���v���V�A���p�ϐ��Ɍv�Z�֐��i108�s�j
double calclambda(mpsconfig &CON)
{
	
	int dimension=CON.get_dimension();	//��͎���
	int Ini_place=CON.get_model_set_way();	//�������q�z�u���@�@0=���� 1=�ז�
	double R=CON.get_re2();			//���v���V�A���p�e�����a
	int size = (int)(R+1);//�v�Z�̈�
	int count=0;
	double dis;//����
	double w;
	double pnd=0;
	double lam=0;
	if(dimension==2)
	{
		if(Ini_place==0)				//�����z�u�Ƃ��Đ����i�q���Ƃ����ꍇ
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-size;j<=size;j++)
				{
					dis=sqrt((double)(i*i+j*j));
					if(dis!=0 && dis<=R)
					{
						double length=dis*CON.get_distancebp();
						w=kernel(R,dis);
						pnd+=w;
						lam+=length*length*w;
						count++;
					}
				}				
			}
		}
		else if(Ini_place==1)				//�����z�u�Ƃ��čז��Z�@���Ƃ����ꍇ
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-2*size;j<=2*size;j++)
				{
					double jj=j*sqrt(3.0)/2;
					double ii=(double)i;
					if(j%2!=0) ii+=0.5;//j����Ȃ�ii��0.5�i�q�������炷
					dis=sqrt(ii*ii+jj*jj);
					if(dis!=0 && dis<=R )
					{
						double length=dis*CON.get_distancebp();
						w=kernel(R,dis);
						pnd+=w;
						lam+=length*length*w;
						count++;
					}			
				}
			}
		}
	}
	else if(dimension==3)
	{
		if(Ini_place==0)				//�����z�u�Ƃ��Đ����i�q���Ƃ����ꍇ
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-size;j<=size;j++)
				{
					for(int k=-size;k<=size;k++)
					{
						dis=sqrt((double)(i*i+j*j+k*k));
						if(dis!=0 && dis<=R)
						{
							double length=dis*CON.get_distancebp();
							w=kernel(R,dis);
							pnd+=w;
							lam+=length*length*w;
							count++;
						}
					}
				}				
			}
		}
		else if(Ini_place==1)				//�����z�u�Ƃ��čז��Z�@���Ƃ����ꍇ
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-2*size;j<=2*size;j++)
				{
					for(int k=-2*size;k<=2*size;k++)
					{
						double ii=(double)i;
						double jj=j*sqrt(3.0)/2;
						double kk=k*sqrt(2.0/3);
						if(j%2!=0) ii+=0.5;//j����Ȃ�ii��0.5�i�q�������炷
						if(k%2!=0) {ii+=0.5; jj+=sqrt(3.0)/6;}//k����Ȃ�ii��jj�����炷
						dis=sqrt(ii*ii+jj*jj+kk*kk);
						if(dis!=0 && dis<=R )
						{
							double length=dis*CON.get_distancebp();
							w=kernel(R,dis);
							pnd+=w;
							lam+=length*length*w;
							count++;
						}
					}			
				}
			}
		}
	}
	lam/=pnd;
	cout<<"�ɂ�count="<<count<<endl;
	return lam;  
}

//���q�f�[�^�ǂݎ��֐��i52�s�j
void input_particle_data(mpsconfig *CON, vector<mpselastic> &PART1, vector<mpselastic> &PART, int t)
{
	int p=0;
	if(t==1)//�ŏ���initial_input.dat����ǂݍ���
	{
		ifstream fin("initial_input.dat");
		if(!fin) cout<<"\"initial_input.dat\" cannot be opened."<<endl;
		fin.unsetf(ifstream::dec);
		fin.setf(ifstream::skipws);
		for(int i=0;i<PART1.size();i++)
		{       
			fin>>PART1[i].ID;
			for(int D=0;D<DIMENSION;D++) fin>>PART1[i].r[D];
			for(int D=0;D<DIMENSION;D++) fin>>PART1[i].u[D];
			for(int D=0;D<DIMENSION;D++) fin>>PART1[i].PAcc[D];
			fin>>PART1[i].P;
			fin>>PART1[i].h;
			fin>>PART1[i].val;
			fin>>PART1[i].type;
			fin>>PART1[i].materialID;	//����͂Ȃ�̈Ӗ�������H�H
			fin>>PART1[i].surface;
			fin>>PART1[i].toFEM;
			PART1[i].dir_Pem=0;			//������
			PART1[i].dir_Pst=0;
			for(int D=0;D<3;D++) PART1[i].ang[D]=0.0;
			PART1[i].ang[3]=1.0;
			for(int D=0;D<3;D++) PART1[i].ang_u[D]=0.0;
		}
		fin.close();	
	}
	///////////////////////*/
	cout<<"���q���בւ��J�n"<<endl;
	if(t!=1) //2STEP�ȍ~��mps_input.dat����ǂݍ���
	{
		ifstream fin("restart_input.dat");
		if(!fin) cout<<"cannot open restart_input.dat"<<endl;
		fin.unsetf(ifstream::dec);
		fin.setf(ifstream::skipws);
		for(int i=0;i<PART.size();i++)
		{
			fin>>PART[i].ID;
			for(int D=0;D<DIMENSION;D++) fin>>PART[i].r[D];
			for(int D=0;D<DIMENSION;D++) fin>>PART[i].u[D];
			for(int D=0;D<DIMENSION;D++) fin>>PART1[i].PAcc[D];
			fin>>PART[i].P;
			fin>>PART[i].h;
			fin>>PART[i].val;
			fin>>PART[i].type;
			fin>>PART[i].materialID;
			fin>>PART[i].surface;
			fin>>PART[i].toFEM;
			PART[i].dir_Pem=0;
			PART[i].dir_Pst=0;
	//		for(int D=0;D<DIMENSION;D++) PART[i].u[D]=0; //�ʒu�����~�����̂ő��x�͂O 
		}
		fin.close();
	}
	bool num=false;
	////���q���בւ�/////
	if(num==true){
	for(int i=0;i<PART1.size();i++)
	{
	 if(PART1[i].type==MAGELAST) 
	 {
		 swap(PART[p],PART1[i]);
		 p++;
	 }
	}
	for(int i=0;i<PART1.size();i++)
	{
	 if(PART1[i].type==MAGELAST2) 
	 {
		 swap(PART[p],PART1[i]);
		 p++;
	 }
	}
	for(int i=0;i<PART1.size();i++)
	{
	 if(PART1[i].type==ELASTIC) 
	 {
		 swap(PART[p],PART1[i]);
		 p++;
	 }
	}
	for(int i=0;i<PART1.size();i++)
	{
	 if(PART1[i].type==TERMINAL1 || PART1[i].type==TERMINAL2 || PART1[i].type==WALL) 
	 {
		 swap(PART[p],PART1[i]);
		 p++;
		 
	 }
	}
	cout<<"���בւ��I��"<<endl;
	//////////////////////////
	}
	else if(num==false){
		for(int i=0;i<PART1.size();i++)
		{
			PART[i]=PART1[i];
		}
	}
}

//���q���J�E���g�֐� & ���ёւ��i28�s�j
//��{�I�Ƀp�[�c�̏��Ԃ�set_initial_placement�Ōv�Z���Ă���̂ŕ��ёւ���K�v�͂Ȃ�
void calc_numbers_of_particles_and_change_the_order(mpsconfig *CON,vector <mpselastic> &PART, int *fluid_number,int *out,int *order_sw)
{
	//�e���q�����J�E���g
	int elastic_num=0;						//�G���X�g�}�[���q
	int magelast_num=0;
	int wall_num=0;
	int solid_num=0;
	int out_num=0;
	int magelast_num2=0;


	for(int i=0;i<PART.size();i++) //map��OK
	{
		if(PART[i].type==ELASTIC) elastic_num++;
		else if(PART[i].type==MAGELAST) magelast_num++;
		else if(PART[i].type==WALL) wall_num++;
		else if(PART[i].type==TERMINAL1) wall_num++;
		else if(PART[i].type==TERMINAL2) wall_num++;
		else if(PART[i].type=MAGELAST2) magelast_num2++;
//		else if(PART[i].type==FLUID) fluid_num++;
//		else if(PART[i].type==INWALL) inwall_num++;
	}
	
	cout<<"elastic: "<<elastic_num<<", magelast: "<<magelast_num+magelast_num2<<", solid: "<<wall_num<<endl;
	out_num=elastic_num+magelast_num+magelast_num2+solid_num+wall_num;	//fluid_number<=i<out��INWALL�Q�ŁAout<=i��OUTWALL�Q�ɂȂ�B
	*out=out_num;
	*fluid_number=elastic_num+magelast_num+magelast_num2;
	if(out_num!=PART.size()){
		cout<<"\n�I���q�����v�G���[�I"<<endl;
		exit(EXIT_FAILURE);
	}

	//���ёւ�
/*	mpselastic PART_temp;			//���ёւ��p���q�N���X
	for(int i=0;i<fluid_num;i++)
	{
		if(PART[i].type!=ELASTIC)//if(PART[i].type!=FLUID)
		{
			cout<<"���ёւ��K�v����"<<endl;
		}
	}
*/	*order_sw=OFF;
}

//INDEX�X�V�֐��i17�s�j
//INDEX: �e�i�q�̂Ȃ��Ɋ܂܂�闱�q���B
//INDEX�𐔂���B�i�q�ԍ��͂Q�����ł͍�����0�BX�����ɂ�{�P�ŁA�E��ōő�(X_mesh*Y_mesh)�B�R�����ł͂y�����ɂ������Ă���
void reload_INDEX(mpsconfig &CON, vector<mpselastic> &PART, int *INDEX)
{       	
	
	int X, Y, Z;	//X, Y, Z�����ɉ��ڂ̊i�q�� 
	int lattice_number;		//���qi���܂ފi�q�̔ԍ�
	double width=CON.get_distancebp()*CON.get_dx();		//�i�q��
//	cout<<"error_check"<<endl;
	for(int i=0;i<CON.get_number_of_mesh();i++) INDEX[i]=0; //
	for(int i=0;i<PART.size();i++)
	{
		//�̈�O���q
		if(!(PART[i].r[A_X]>CON.get_minX() && PART[i].r[A_X]<CON.get_maxX())) cout<<"X="<<PART[i].r[A_X]<<", i="<<i<<endl;
		else if(!(PART[i].r[A_Y]>CON.get_minY() && PART[i].r[A_Y]<CON.get_maxY())) cout<<"Y="<<PART[i].r[A_Y]<<", i="<<i<<endl;
//		else if(!(PART[i].r[A_Z]>CON->get_minZ() && PART[i].r[A_Z]<CON->get_maxZ())) cout<<"Z="<<PART[i].r[A_Z]<<", i="<<i<<endl; 

		//////////////////

		X=(int)((PART[i].r[A_X]-CON.get_minX())/width);
		Y=(int)((PART[i].r[A_Y]-CON.get_minY())/width);
		Z=(int)((PART[i].r[A_Z]-CON.get_minZ())/width);

		lattice_number=(Z*CON.get_X_mesh()*CON.get_Y_mesh())+(Y*CON.get_X_mesh())+X;

        PART[i].index=lattice_number; //���qi�����Ԗڂ̊i�q�Ɋ܂܂�Ă��邩���v�Z
		INDEX[lattice_number]++; //�i�q�Ɋ܂܂�闱�q�����v�Z
	}
//	cout<<"error_check_end"<<endl;
}

//INDEX�X�V�֐����̂Q MESH�ɗ��q�ԍ����i�[����i15�s�j
//MESH: int **MESH = new int *[CON.get_number_of_mesh()];	
//get_number_of_mesh() {return (int)((maxX-minX)/(distancebp*dx)*(maxY-minY)/(distancebp*dx)*(maxZ-minZ)/(distancebp*dx)+0.001);} //�i�q���FX_mesh*Y_mesh*Z_mesh
void reload_INDEX2(mpsconfig *CON, vector<mpselastic> &PART, int **MESH)
{
	int number_of_mesh = CON->get_number_of_mesh();
	int *count = new int [number_of_mesh];
	for(int i=0;i<number_of_mesh;i++) count[i]=0;//������

	for(int i=0;i<PART.size();i++)
	{
		int lattice_number=PART[i].index; //���qi�����Ԗڂ̊i�q�Ɋ܂܂�Ă��邩���擾
        
		MESH[lattice_number][count[lattice_number]]=i;
		count[lattice_number]++;
	}
	delete [] count;
}

//�@���x�N�g���쐬�֐��i30�s�j
void direct_f(mpsconfig &CON,vector<mpselastic> &PART,int i,double *direct[DIMENSION])
{
	
	double R=CON.get_re3()*CON.get_distancebp();//�@���x�N�g���v�Z�ɗ��p����e�����a

	double px=PART[i].r[A_X]+CON.get_distancebp();//x+le
	double mx=PART[i].r[A_X]-CON.get_distancebp();//x-le
	double py=PART[i].r[A_Y]+CON.get_distancebp();//y+le
	double my=PART[i].r[A_Y]-CON.get_distancebp();//y-le
	double pz=PART[i].r[A_Z]+CON.get_distancebp();//z+le
	double mz=PART[i].r[A_Z]-CON.get_distancebp();//z-le
	
	double pnd_px=pnd_for_direct(CON,PART,px,PART[i].r[A_Y],PART[i].r[A_Z],R,i);
	double pnd_mx=pnd_for_direct(CON,PART,mx,PART[i].r[A_Y],PART[i].r[A_Z],R,i);
	double pnd_py=pnd_for_direct(CON,PART,PART[i].r[A_X],py,PART[i].r[A_Z],R,i);
	double pnd_my=pnd_for_direct(CON,PART,PART[i].r[A_X],my,PART[i].r[A_Z],R,i);
	double pnd_pz=pnd_for_direct(CON,PART,PART[i].r[A_X],PART[i].r[A_Y],pz,R,i);
	double pnd_mz=pnd_for_direct(CON,PART,PART[i].r[A_X],PART[i].r[A_Y],mz,R,i);

	direct[A_X][i]=(pnd_px-pnd_mx)/(2*CON.get_distancebp());
	direct[A_Y][i]=(pnd_py-pnd_my)/(2*CON.get_distancebp());
	direct[A_Z][i]=(pnd_pz-pnd_mz)/(2*CON.get_distancebp());
	
	double a=sqrt(direct[A_X][i]*direct[A_X][i]+direct[A_Y][i]*direct[A_Y][i]+direct[A_Z][i]*direct[A_Z][i]);
	if(a!=0)
	{ 
		direct[A_X][i]/=a;
		direct[A_Y][i]/=a;
		direct[A_Z][i]/=a;
	}
}

//direct�p���q�����x����i25�s�j
double pnd_for_direct(mpsconfig &CON,vector<mpselastic> &PART,double x,double y,double z,double R,int i)
{
	//���qi�̈ʒu�Ƃ͂����̂ŁAMESH���g�p�����ق����悢
	//���ȏ��ł͌����������Ă��邪�A�����ł͏d�݊֐���p����B
	//R=CON->get_re3()*CON.get_distancebp();
	double spnd=0;

	for(int k=0;k<PART[i].N3;k++)
	{       
		int j=PART[i].NEI3[k];
		double X=PART[j].r[A_X]-x;
		double Y=PART[j].r[A_Y]-y;
		double Z=PART[j].r[A_Z]-z;
		double dis=sqrt(X*X+Y*Y+Z*Z);
		//if(dis<R) spnd++;   //���ȏ��ǂ���
		if(dis<R)
		{
			double w=(1-dis/R)*(1-dis/R);
			spnd+=w;
		}
	}
	return spnd;
}

//�N�[������
void courant_elastic(vector<mpselastic> &PART, int fluid_number, int t, double *dt, double mindis, double Umax, double *g)
{
	mpsconfig CON;
	double CFL=CON.get_courant();	//�N�[����������
	double CFL2=0;
	double le=mindis;				//�ŒZ���q�ԋ���
	double E=CON.get_E_m();		//MAGELAST�̃����O��
	double v=CON.get_v_m();
	double lambda=v*E/((1.0+v)*(1.0-2.0*v));
	double mu=E/(2.0*(1.0+v));
//	double density=CON->Get_density();
	double factor=20.0;
	double uu;
	double cfl2;
	ofstream aaa("CFL.dat", ios::app);
	for(int i=0;i<PART.size();i++){
		if(PART[i].type==MAGELAST || PART[i].type==MAGELAST2)CFL=sqrt((lambda+2.0*mu)/CON.Get_MRE_density());
		else if(PART[i].type==ELASTIC)CFL=sqrt((lambda+2.0*mu)/CON.Get_Silicone_density());
		uu=pow(pow(PART[i].u[A_X],2)+pow(PART[i].u[A_Y],2)+pow(PART[i].u[A_Z],2),0.5);
		cfl2=(*dt*uu)/0.001;
		if(cfl2>CFL2) CFL2=cfl2; 
	}
	aaa<<"le="<<le<<", dt= "<<*dt<<", CFL= "<<CFL<<", CFL2="<<CFL2<<", new dt="<<le/CFL<<endl;
	aaa.close();
/*	///////////�N�[������//////////////////
	if(CFL>0)
	{      
		double newdt=*dt;	//�V����dt
		CFL=sqrt((lambda+2.0*mu)/density);
	//	if(Umax!=0) newdt=CFL*le/Umax;
		if(Umax!=0) newdt=(le/CFL)/factor;
		if(newdt>CON->get_dt()) *dt=CON->get_dt();
		else *dt=newdt;
		
		if(*dt!=CON->get_dt()) cout<<"CFL�����ɂ��dt�X�V dt="<<*dt<<endl;
	}
	///�g�U���̐��m�Ȓ�`�𒲂ׂď����Ȃ���
	if(CON->get_vis()!=0 && CON->get_vis_calc_type()==POSITIVE)
	{
		if(*dt>0.25*le*le/CON->get_vis()) cout<<"�g�U���ᔽ"<<endl;
	}
	/////////////////////////////*/
}

//CG�@�i90�s�j
void CG_method(mpsconfig *CON,double *r,double *P,double *AP,double *val,int *ind,int *ptr,int pn,double *X,int *countN,double EP)
{
	cout<<"CG�@�X�^�[�g------";
	//unsigned int timeCG=GetTickCount();
	int count=0;
	double rr=0;
	double E=1;//�덷
	double alp,beta;
	for(int n=0;n<pn;n++) rr+=r[n]*r[n];
	
	if(CON->get_omp_P()==OFF)//�ʏ��
	{
		while(E>EP)// EP=CON->get_CGep();//��������(convergence test)
		{
			count++;
			//////////////alp�����߂�
			for(int n=0;n<pn;n++)
			{      
				AP[n]=0;
				for(int j=ptr[n];j<ptr[n+1];j++) AP[n]+=val[j]*P[ind[j]];
			}
			double PAP=0;
			for(int n=0;n<pn;n++)  PAP+=P[n]*AP[n];
			alp=rr/PAP;
		//	cout<<"alp="<<alp<<" rr="<<rr<<" PAP="<<PAP<<endl;
			//////////////////////
		
			//////////////// ���X�V�@X(k+1)=X(k)+alp*P
			for(int n=0;n<pn;n++) X[n]+=alp*P[n];
			//////////////////////////////
			
			//////////////// r=r-alp*AP
			for(int n=0;n<pn;n++) r[n]-=alp*AP[n];
			/////////////////////////////
			
			///////////////////////beta
			beta=1.0/rr;
			rr=0;
			for(int n=0;n<pn;n++) rr+=r[n]*r[n];
			beta=beta*rr;
			///////////////////////

			//////////////////�덷
			E=sqrt(rr);
			//cout<<"E="<<E<<endl;
			////////////////////////
			
			///////////////////// P=r+beta*P
			for(int n=0;n<pn;n++) P[n]=r[n]+beta*P[n];
		}
	}
	else if(CON->get_omp_P()==ON)//openMP���g�p����ꍇ
	{
		while(E>EP)
		{
			count++;
			//////////////alp�����߂�
			double PAP=0;
			#pragma omp parallel for reduction(+:PAP)
			for(int n=0;n<pn;n++)
			{
				AP[n]=0;
				for(int j=ptr[n];j<ptr[n+1];j++) AP[n]+=val[j]*P[ind[j]];
				PAP+=P[n]*AP[n];
			}
			alp=rr/PAP;
			//////////////////////

			//////////////
			E=0;//�덷
			beta=1.0/rr;
			rr=0;
			#pragma omp parallel for reduction(+:rr)
			for(int n=0;n<pn;n++) 
			{
				X[n]+=alp*P[n];// ���X�V�@X(k+1)=X(k)+alp*P
				r[n]-=alp*AP[n];// r=r-alp*AP
				rr+=r[n]*r[n];
			}
			E=sqrt(rr);
			//cout<<"E="<<E<<endl;
		
			beta=beta*rr;///beta
			
			for(int n=0;n<pn;n++) P[n]=r[n]+beta*P[n];/// P=r+beta*P
		}
	}
	*countN=count;//�����񐔂�n��
}

//ICCG�@(230�s)
void iccg(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X,double *r,double *P,double EP,int *count2)
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

	num2=0;
	for(int k=0;k<pn;k++)
	{	
		ptr2[k]=num2;
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
	
	for(int n=0;n<pn;n++) P[n]=LDLt_r[n];

	cout<<"ICCG�@:���m��="<<pn<<" ---";
	unsigned int time=GetTickCount();
	int count=0;
	double ep=EP;//��������
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
		//cout<<"alp="<<alp<<endl;
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
		//cout<<"E="<<E<<endl;
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
	}
	//cout<<"������="<<count<<" time="<<(GetTickCount()-time)*0.001<<"/";
		
	delete [] AP;

	delete [] y;
	delete [] LDLt_r;
	delete [] D1;

	delete [] val2;
	delete [] ind2;
	delete [] ptr2;

	for(int i=0;i<pn;i++)
	{
		delete [] VAL[i];
		delete [] IND[i];
	}
	delete [] VAL;
	delete [] IND;
	delete [] NUM;
	*count2=count;//�����񐔂��i�[���ĕԂ�
}


//�K�E�X�̏����@�i24�s�j
//���͍ŏI�I��B�̂Ȃ���
void gauss(double *matrix,double *B,int N)
{
	for(int k=0;k<N;k++)
	{
		double akk=matrix[k*N+k];
		
		for(int i=0;i<N;i++)
		{
			if(i!=k)
			{
				double A=matrix[i*N+k]/akk;
				//for(int j=0;j<N;j++)
				for(int j=k;j<N;j++)
				{					
					matrix[i*N+j]-=A*matrix[k*N+j];
				}
				B[i]-=A*B[k];
			}
		}
	}
	for(int k=0;k<N;k++) B[k]/=matrix[k*N+k];

}

//���̑��x����шʒu����i46�s�j
void renewal_u_and_r_in_positive(vector<mpselastic> &PART,int fluid_number,int t,double dt,double *Umax,double **potential,double **laplacian,double *g,double **previous_Un,double **F)
{
	mpsconfig CON;
	double U=0;						//�ő呬�x
	double vis;//=CON->Get_vis();
	double mass=CON.get_particle_mass();	//���q�̎���
	int d=CON.get_dimension();
	int sw=CON.get_temporary_r();	//�n�m�Ȃ牼�̈ʒu���v�Z����

	double *old_U[DIMENSION];
	for(int D=0;D<DIMENSION;D++) old_U[D]=new double [fluid_number];//�ύX�O�̑��x���L�����Ă���

	if(t==1) for(int i=0;i<fluid_number;i++)for(int D=0;D<DIMENSION;D++) previous_Un[D][i]=0;//t=1�̂Ƃ��͏�����      
			

	//potential[D][i]���ꍇ�ɂ���Ă̓[���ɏ���������
	if(CON.get_dir_for_P()==1 || CON.get_dir_for_P()==3) //�\�ʗ��q�̕\�ʒ��͈͂��͒l�Ƃ��Čv�Z����Ă���̂�,�����ł͍l�����Ȃ��悤����������
	{
		for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++) potential[D][i]=0;
	}
	//////////////////////*/

	/////////////���x�X�V
	for(int i=0;i<fluid_number;i++)
	{        
		double speed=0;//���q���x
		for(int D=0;D<d;D++)
		{   
			if(PART[i].type==ELASTIC) vis=CON.Get_Silicone_vis();
			else if(PART[i].type==MAGELAST || PART[i].type==MAGELAST2) vis=CON.Get_MRE_vis();
			old_U[D][i]=PART[i].u[D];
			
			PART[i].u[D]+=dt*(vis *laplacian[D][i]+potential[D][i]+g[D]+F[D][i]/mass);
			
			//PART[i].u[D]=previous_Un[D][i]+dt*(vis*laplacian[D][i]+potential[D][i]+g[D]);//�^�Ƃі@
			speed+=PART[i].u[D]*PART[i].u[D];
		}
		if(speed>U) U=speed;
	}
	*Umax=U;

	//�ʒu�X�V
	if(sw==ON) for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++) PART[i].r[D]+=dt*0.5*(PART[i].u[D]+old_U[D][i]);//��`��
	
}

//���x���U�v�Z�֐��i27�s�j
double divergence(mpsconfig *CON,vector<mpselastic> &PART,int i,double n0)
{
    double W=0;										//���q�����x
    double R=CON->get_distancebp()*CON->get_re();	//�e�����a
    double div=0;									//���U�̒l

	for(int k=0;k<PART[i].N;k++)
    {    
        int j=PART[i].NEI[k]; 
        double X=PART[j].r[A_X]-PART[i].r[A_X];
		double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
		double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
		double dis=sqrt(X*X+Y*Y+Z*Z);
			       
		double w=kernel(R,dis);
		
		div+=(PART[j].u[A_X]-PART[i].u[A_X])*X*w/(dis*dis);
		div+=(PART[j].u[A_Y]-PART[i].u[A_Y])*Y*w/(dis*dis);
		div+=(PART[j].u[A_Z]-PART[i].u[A_Z])*Z*w/(dis*dis);
		W+=w;
    }
    if(W!=0)
	{
		div*=CON->get_dimension()/W;
	}
    return div;
}


///���x���U�v�Z�֐�(WLSM�@)�i198�s�j
///WLSM=Weighed Least Square Method�F�d�ݕt���ŏ����@
double divergence2(mpsconfig *CON,vector<mpselastic> &PART,int i)
{
    double div=0;//���U�̒l
	double le=CON->get_distancebp();
	
	double r=CON->get_re();
	double R=r*le;
	int d=CON->get_dimension();
	int N=0;					//�W���s��̌�
	int order=1;				//�ߎ��Ȗʂ̵��ް�B 1=���` 2=��
    
	//�W���s��̑傫���̌���
	if(d==2)
	{
		if(order==1) N=2;
		else if(order==2) N=5;
	}
	else if(d==3)
	{
		if(order==1) N=4;
		else if(order==2) N=10;
	}
	////////////////////////////////

	double *matrix=new double [N*N];	//N�~N�̌W���s��
	double *B1=new double [N];			//N�̉��s��
	double *B2=new double [N];			//N�̉��s��
	double *B3=new double [N];			//N�̉��s��

	for(int n=0;n<N*N;n++) matrix[n]=0;	//������
	for(int n=0;n<N;n++) {B1[n]=0;B2[n]=0;B3[n]=0;}

	if(d==2 && order==1)				//�񎟌�
	{
		for(int k=0;k<PART[i].N;k++)
		{
			int j=PART[i].NEI[k];
			
			double X=(PART[j].r[A_X]-PART[i].r[A_X])/le;// le�Ŋ���̂͑ł��؂�덷�h�~
			double Y=(PART[j].r[A_Y]-PART[i].r[A_Y])/le;
			double U=(PART[j].u[A_X]-PART[i].u[A_X]);
			double V=(PART[j].u[A_Y]-PART[i].u[A_Y]);
			double dis=sqrt(X*X+Y*Y);
					
			double w=1;
			if(dis>1) w=1/(dis*dis*dis*dis);
					
			matrix[0]+=X*X*w;			//��Xjwj
			matrix[1]+=X*Y*w;		//��XjYjwj
			matrix[3]+=Y*Y*w;			//��Yjwj
				
			B1[0]+=U*X*w;//��ujXjwj
			B1[1]+=U*Y*w;//��ujYjwj
			B2[0]+=V*X*w;//��vjXjwj
			B2[1]+=V*Y*w;//��vjYjwj
		}
			
		matrix[2]=matrix[1];		//��XjYjwj
		for(int n=0;n<N;n++)
		{
			B1[n]/=le;//�ł��؂�덷�h�~
			B2[n]/=le;//�ł��؂�덷�h�~
		}

		double determinant=matrix[0]*matrix[3]-matrix[1]*matrix[2];//�s��
			
		double dudx=(B1[0]*matrix[3]-matrix[1]*B1[1])/determinant;
		double dvdy=(B2[1]*matrix[0]-matrix[2]*B2[0])/determinant;

		div=dudx+dvdy;		
	}
	if(d==2 && order==2)//�񎟌�2����
	{
		for(int k=0;k<PART[i].N;k++)
		{
			int j=PART[i].NEI[k];
			
			double X=(PART[j].r[A_X]-PART[i].r[A_X]);// le�Ŋ���̂͑ł��؂�덷�h�~
			double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
			double U=(PART[j].u[A_X]-PART[i].u[A_X]);
			double V=(PART[j].u[A_Y]-PART[i].u[A_Y]);
			double dis=sqrt(X*X+Y*Y);
					
			double w=1;
			//if(dis>1) w=r*r*r*r/(dis*dis*dis*dis);
			if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);
					
			matrix[0]+=X*X*w;			//��Xjwj
			matrix[1]+=X*Y*w;		//��XjYjwj
			matrix[2]+=X*X*X*w;			//��Xj^3wj
			matrix[3]+=X*X*Y*w;			//��Xj^2Yjwj
			matrix[4]+=X*Y*Y*w;			//��XjYj^2wj

			matrix[6]+=Y*Y*w;			//��Yj^2wj
			matrix[9]+=Y*Y*Y*w;			//��Yj^3wj

			matrix[12]+=X*X*X*X*w;			//��Xj^4wj
			matrix[13]+=X*X*X*Y*w;			//��Xj^3Yjwj
			matrix[14]+=X*X*Y*Y*w;			//��Xj^2Yj^2wj
	
			matrix[19]+=X*Y*Y*Y*w;			//��XjYj^3wj

			matrix[24]+=Y*Y*Y*Y*w;			//��Yj^4wj

				
			B1[0]+=U*X*w;//��ujXjwj
			B1[1]+=U*Y*w;//��ujYjwj
			B1[2]+=U*X*X*w;//��ujXj^2wj
			B1[3]+=U*X*Y*w;//��ujXjYjwj
			B1[4]+=U*Y*Y*w;//��ujYj^2wj

			B2[0]+=V*X*w;//��vjXjwj
			B2[1]+=V*Y*w;//��vjYjwj
			B2[2]+=V*X*X*w;//��vjXj^2wj
			B2[3]+=V*X*Y*w;//��vjXjYjwj
			B2[4]+=V*Y*Y*w;//��vjYj^2wj
		}
			
		matrix[5]=matrix[1];		//��XjYjwj
		matrix[7]=matrix[3];		//��Xj^2Yjwj
		matrix[8]=matrix[4];		//��XjYj^2wj
		matrix[10]=matrix[2];		//��Xj^3Yjwj
		matrix[11]=matrix[3];		//��Xj^2Yjwj
		matrix[15]=matrix[3];		//��Xj^2Yjwj
		matrix[16]=matrix[4];		//��XjYj^2wj
		matrix[17]=matrix[13];		//��Xj^3Yjwj
		matrix[18]=matrix[14];		//��Xj^2Yj^2wj
		matrix[20]=matrix[4];		//��XjYj^2wj
		matrix[21]=matrix[9];		//��Yj^3wj
		matrix[22]=matrix[14];		//��Xj^2Yj^2wj
		matrix[23]=matrix[19];		//��XjYj^3wj

		///�ۂߌ덷�h�~
		for(int n=0;n<N;n++)
		{
			for(int m=0;m<N;m++) matrix[n*N+m]*=1e7;
			B1[n]*=1e7;
			B2[n]*=1e7;
		}//*/

		double dudx=0;
		double dvdy=0;
		return_X_for5N(matrix,N,B1,B2,&dudx,&dvdy);//5���A���������́A���P�C�Q��Ԃ��֐�
		
		div=(dudx+dvdy);
			
	}
	else if(d==3 && order==1)//3����1����
	{
		if(PART[i].N>5)
		{
			//div=calc_WLSM_divu_D3_order1(CON,PART,matrix,B1,B2,B3,i,3);//�Ō�̈����͖��m��
			div=calc_WLSM_divu_D3_order1_2(CON,PART,matrix,B1,B2,B3,i,4);//�Ō�̈����͖��m��
		}
		else 
		{
			div=0;
			//cout<<"�x�� �ߗח��q����5�ȉ�("<<PART[i].N<<")�ł�"<<endl;
		}
	}
	else if(d==3 && order==2)//3����2����
	{
		//P=Pi+a��x+b��y+c��z+d��x2+e��y2+f��z2+g��x��y+h��y��z+i��z��x�Ƃ����ƁA
		///�W���s���
		///   ����x2      ����x��y    ����x��z    ����x3       ����x��y2    ����x��z2    ����x2��y     ����x��y��z  ����x2��z     a = ����x��P  
		///   ����x��y    ����y2      ����y��z    ����x2��y    ����y3       ����y��z2    ����x��y2     ����y2��z    ����x��y��z   b = ����y��P
		///   ����x��z    ����y��z    ����z2      ����x2��z    ����y2��z    ����z3       ����x��y��z   ����y��z2    ����x��z2     c = ����z��P
		///   ����x3      ����x2��y   ����x2��z   ����x4       ����x2��y2   ����x2��z2   ����x3��y     ����x2��y��z ����x3��z     d = ����x2��P
		///   ����x��y2   ����y3      ����y2��z   ����x2��y2   ����y4       ����y2��z2   ����x��y3     ����y3��z    ����x��y2��z  e = ����y2��P
		///   ����x��z2   ����y��z2   ����z3      ����x2��z2   ����y2��z2   ����z4       ����x��y��z2  ����y��z3    ����x��z3     f = ����z2��P
		///   ����x2��y   ����x��y2   ����x��y��z ����x3��y    ����x��y3    ����x��y��z2 ����x2��y2    ����x��y2��z ����x2��y��z  g = ����x��y��P
		///   ����x��y��z ����y2��z   ����y��z2   ����x2��y��z ����y3��z    ����y��z3    ����x��y2��z  ����y2��z2   ����x��y��z2  h = ����y��z��P
		///   ����x2��z   ����x��y��z ����x��z2   ����x3��z    ����x��y2��z  ����x��z3   ����x2��y ��z ����x��y��z2 ����x2��z2    g = ����x��z��P
		
		if(PART[i].N>8)
		{
			//div=calc_WLSM_divu_D3_order2(CON,PART,matrix,B1,B2,B3,i,9);
			div=calc_WLSM_divu_D3_order2_2(CON,PART,matrix,B1,B2,B3,i,10);
		}
		else if(PART[i].N>5)
		{
			//div=calc_WLSM_divu_D3_order1(CON,PART,matrix,B1,B2,B3,i,3);//���̏ꍇ�A�Ō�̈�����N=3��n�����Ƃɒ���
			div=calc_WLSM_divu_D3_order1_2(CON,PART,matrix,B1,B2,B3,i,4);//�Ō�̈����͖��m��
		}
		else 
		{
			div=0;
			//cout<<"�x�� �ߗח��q����9�ȉ�("<<PART[i].N<<")�ł�"<<endl;
		}
	}

	delete [] matrix;
	delete [] B1;
	delete [] B2;
	delete [] B3;

    return div;
}

//5���̘A���������̉�1,2��Ԃ��֐��i24�s�j
void return_X_for5N(double *matrix,int N,double *B1,double *B2,double *dudx,double *dudy)
{
	double a11,a12,a13,a14,a15,a21,a22,a23,a24,a25,a31,a32,a33,a34,a35,a41,a42,a43,a44,a45,a51,a52,a53,a54,a55;
	double b1,b2,b3,b4,b5;
	double c1,c2,c3,c4,c5;

	a11=matrix[0];a12=matrix[1];a13=matrix[2];a14=matrix[3];a15=matrix[4];
	a21=matrix[5];a22=matrix[6];a23=matrix[7];a24=matrix[8];a25=matrix[9];
	a31=matrix[10];a32=matrix[11];a33=matrix[12];a34=matrix[13];a35=matrix[14];
	a41=matrix[15];a42=matrix[16];a43=matrix[17];a44=matrix[18];a45=matrix[19];
	a51=matrix[20];a52=matrix[21];a53=matrix[22];a54=matrix[23];a55=matrix[24];

	b1=B1[0];b2=B1[1];b3=B1[2];b4=B1[3];b5=B1[4];
	c1=B2[0];c2=B2[1];c3=B2[2];c4=B2[3];c5=B2[4];
	
	double determinant=(a11*a22*a33*a44*a55-a11*a22*a33*a45*a54-a11*a22*a34*a43*a55+a11*a22*a34*a45*a53+a11*a22*a35*a43*a54-a11*a22*a35*a44*a53-a11*a23*a32*a44*a55+a11*a23*a32*a45*a54+a11*a23*a34*a42*a55-a11*a23*a34*a45*a52-a11*a23*a35*a42*a54+a11*a23*a35*a44*a52+a11*a24*a32*a43*a55-a11*a24*a32*a45*a53-a11*a24*a33*a42*a55+a11*a24*a33*a45*a52+a11*a24*a35*a42*a53-a11*a24*a35*a43*a52-a11*a25*a32*a43*a54+a11*a25*a32*a44*a53+a11*a25*a33*a42*a54-a11*a25*a33*a44*a52-a11*a25*a34*a42*a53+a11*a25*a34*a43*a52-a12*a21*a33*a44*a55+a12*a21*a33*a45*a54+a12*a21*a34*a43*a55-a12*a21*a34*a45*a53-a12*a21*a35*a43*a54+a12*a21*a35*a44*a53+a12*a23*a31*a44*a55-a12*a23*a31*a45*a54-a12*a23*a34*a41*a55+a12*a23*a34*a45*a51+a12*a23*a35*a41*a54-a12*a23*a35*a44*a51-a12*a24*a31*a43*a55+a12*a24*a31*a45*a53+a12*a24*a33*a41*a55-a12*a24*a33*a45*a51-a12*a24*a35*a41*a53+a12*a24*a35*a43*a51+a12*a25*a31*a43*a54-a12*a25*a31*a44*a53-a12*a25*a33*a41*a54+a12*a25*a33*a44*a51+a12*a25*a34*a41*a53-a12*a25*a34*a43*a51+a13*a21*a32*a44*a55-a13*a21*a32*a45*a54-a13*a21*a34*a42*a55+a13*a21*a34*a45*a52+a13*a21*a35*a42*a54-a13*a21*a35*a44*a52-a13*a22*a31*a44*a55+a13*a22*a31*a45*a54+a13*a22*a34*a41*a55-a13*a22*a34*a45*a51-a13*a22*a35*a41*a54+a13*a22*a35*a44*a51+a13*a24*a31*a42*a55-a13*a24*a31*a45*a52-a13*a24*a32*a41*a55+a13*a24*a32*a45*a51+a13*a24*a35*a41*a52-a13*a24*a35*a42*a51-a13*a25*a31*a42*a54+a13*a25*a31*a44*a52+a13*a25*a32*a41*a54-a13*a25*a32*a44*a51-a13*a25*a34*a41*a52+a13*a25*a34*a42*a51-a14*a21*a32*a43*a55+a14*a21*a32*a45*a53+a14*a21*a33*a42*a55-a14*a21*a33*a45*a52-a14*a21*a35*a42*a53+a14*a21*a35*a43*a52+a14*a22*a31*a43*a55-a14*a22*a31*a45*a53-a14*a22*a33*a41*a55+a14*a22*a33*a45*a51+a14*a22*a35*a41*a53-a14*a22*a35*a43*a51-a14*a23*a31*a42*a55+a14*a23*a31*a45*a52+a14*a23*a32*a41*a55-a14*a23*a32*a45*a51-a14*a23*a35*a41*a52+a14*a23*a35*a42*a51+a14*a25*a31*a42*a53-a14*a25*a31*a43*a52-a14*a25*a32*a41*a53+a14*a25*a32*a43*a51+a14*a25*a33*a41*a52-a14*a25*a33*a42*a51+a15*a21*a32*a43*a54-a15*a21*a32*a44*a53-a15*a21*a33*a42*a54+a15*a21*a33*a44*a52+a15*a21*a34*a42*a53-a15*a21*a34*a43*a52-a15*a22*a31*a43*a54+a15*a22*a31*a44*a53+a15*a22*a33*a41*a54-a15*a22*a33*a44*a51-a15*a22*a34*a41*a53+a15*a22*a34*a43*a51+a15*a23*a31*a42*a54-a15*a23*a31*a44*a52-a15*a23*a32*a41*a54+a15*a23*a32*a44*a51+a15*a23*a34*a41*a52-a15*a23*a34*a42*a51-a15*a24*a31*a42*a53+a15*a24*a31*a43*a52+a15*a24*a32*a41*a53-a15*a24*a32*a43*a51-a15*a24*a33*a41*a52+a15*a24*a33*a42*a51);
	
	*dudx=(b1*a22*a33*a44*a55-b1*a22*a33*a45*a54-b1*a22*a34*a43*a55+b1*a22*a34*a45*a53+b1*a22*a35*a43*a54-b1*a22*a35*a44*a53-b1*a23*a32*a44*a55+b1*a23*a32*a45*a54+b1*a23*a34*a42*a55-b1*a23*a34*a45*a52-b1*a23*a35*a42*a54+b1*a23*a35*a44*a52+b1*a24*a32*a43*a55-b1*a24*a32*a45*a53-b1*a24*a33*a42*a55+b1*a24*a33*a45*a52+b1*a24*a35*a42*a53-b1*a24*a35*a43*a52-b1*a25*a32*a43*a54+b1*a25*a32*a44*a53+b1*a25*a33*a42*a54-b1*a25*a33*a44*a52-b1*a25*a34*a42*a53+b1*a25*a34*a43*a52-a12*b2*a33*a44*a55+a12*b2*a33*a45*a54+a12*b2*a34*a43*a55-a12*b2*a34*a45*a53-a12*b2*a35*a43*a54+a12*b2*a35*a44*a53+a12*a23*b3*a44*a55-a12*a23*b3*a45*a54-a12*a23*a34*b4*a55+a12*a23*a34*a45*b5+a12*a23*a35*b4*a54-a12*a23*a35*a44*b5-a12*a24*b3*a43*a55+a12*a24*b3*a45*a53+a12*a24*a33*b4*a55-a12*a24*a33*a45*b5-a12*a24*a35*b4*a53+a12*a24*a35*a43*b5+a12*a25*b3*a43*a54-a12*a25*b3*a44*a53-a12*a25*a33*b4*a54+a12*a25*a33*a44*b5+a12*a25*a34*b4*a53-a12*a25*a34*a43*b5+a13*b2*a32*a44*a55-a13*b2*a32*a45*a54-a13*b2*a34*a42*a55+a13*b2*a34*a45*a52+a13*b2*a35*a42*a54-a13*b2*a35*a44*a52-a13*a22*b3*a44*a55+a13*a22*b3*a45*a54+a13*a22*a34*b4*a55-a13*a22*a34*a45*b5-a13*a22*a35*b4*a54+a13*a22*a35*a44*b5+a13*a24*b3*a42*a55-a13*a24*b3*a45*a52-a13*a24*a32*b4*a55+a13*a24*a32*a45*b5+a13*a24*a35*b4*a52-a13*a24*a35*a42*b5-a13*a25*b3*a42*a54+a13*a25*b3*a44*a52+a13*a25*a32*b4*a54-a13*a25*a32*a44*b5-a13*a25*a34*b4*a52+a13*a25*a34*a42*b5-a14*b2*a32*a43*a55+a14*b2*a32*a45*a53+a14*b2*a33*a42*a55-a14*b2*a33*a45*a52-a14*b2*a35*a42*a53+a14*b2*a35*a43*a52+a14*a22*b3*a43*a55-a14*a22*b3*a45*a53-a14*a22*a33*b4*a55+a14*a22*a33*a45*b5+a14*a22*a35*b4*a53-a14*a22*a35*a43*b5-a14*a23*b3*a42*a55+a14*a23*b3*a45*a52+a14*a23*a32*b4*a55-a14*a23*a32*a45*b5-a14*a23*a35*b4*a52+a14*a23*a35*a42*b5+a14*a25*b3*a42*a53-a14*a25*b3*a43*a52-a14*a25*a32*b4*a53+a14*a25*a32*a43*b5+a14*a25*a33*b4*a52-a14*a25*a33*a42*b5+a15*b2*a32*a43*a54-a15*b2*a32*a44*a53-a15*b2*a33*a42*a54+a15*b2*a33*a44*a52+a15*b2*a34*a42*a53-a15*b2*a34*a43*a52-a15*a22*b3*a43*a54+a15*a22*b3*a44*a53+a15*a22*a33*b4*a54-a15*a22*a33*a44*b5-a15*a22*a34*b4*a53+a15*a22*a34*a43*b5+a15*a23*b3*a42*a54-a15*a23*b3*a44*a52-a15*a23*a32*b4*a54+a15*a23*a32*a44*b5+a15*a23*a34*b4*a52-a15*a23*a34*a42*b5-a15*a24*b3*a42*a53+a15*a24*b3*a43*a52+a15*a24*a32*b4*a53-a15*a24*a32*a43*b5-a15*a24*a33*b4*a52+a15*a24*a33*a42*b5)/determinant;
	*dudy=(a11*c2*a33*a44*a55-a11*c2*a33*a45*a54-a11*c2*a34*a43*a55+a11*c2*a34*a45*a53+a11*c2*a35*a43*a54-a11*c2*a35*a44*a53-a11*a23*c3*a44*a55+a11*a23*c3*a45*a54+a11*a23*a34*c4*a55-a11*a23*a34*a45*c5-a11*a23*a35*c4*a54+a11*a23*a35*a44*c5+a11*a24*c3*a43*a55-a11*a24*c3*a45*a53-a11*a24*a33*c4*a55+a11*a24*a33*a45*c5+a11*a24*a35*c4*a53-a11*a24*a35*a43*c5-a11*a25*c3*a43*a54+a11*a25*c3*a44*a53+a11*a25*a33*c4*a54-a11*a25*a33*a44*c5-a11*a25*a34*c4*a53+a11*a25*a34*a43*c5-c1*a21*a33*a44*a55+c1*a21*a33*a45*a54+c1*a21*a34*a43*a55-c1*a21*a34*a45*a53-c1*a21*a35*a43*a54+c1*a21*a35*a44*a53+c1*a23*a31*a44*a55-c1*a23*a31*a45*a54-c1*a23*a34*a41*a55+c1*a23*a34*a45*a51+c1*a23*a35*a41*a54-c1*a23*a35*a44*a51-c1*a24*a31*a43*a55+c1*a24*a31*a45*a53+c1*a24*a33*a41*a55-c1*a24*a33*a45*a51-c1*a24*a35*a41*a53+c1*a24*a35*a43*a51+c1*a25*a31*a43*a54-c1*a25*a31*a44*a53-c1*a25*a33*a41*a54+c1*a25*a33*a44*a51+c1*a25*a34*a41*a53-c1*a25*a34*a43*a51+a13*a21*c3*a44*a55-a13*a21*c3*a45*a54-a13*a21*a34*c4*a55+a13*a21*a34*a45*c5+a13*a21*a35*c4*a54-a13*a21*a35*a44*c5-a13*c2*a31*a44*a55+a13*c2*a31*a45*a54+a13*c2*a34*a41*a55-a13*c2*a34*a45*a51-a13*c2*a35*a41*a54+a13*c2*a35*a44*a51+a13*a24*a31*c4*a55-a13*a24*a31*a45*c5-a13*a24*c3*a41*a55+a13*a24*c3*a45*a51+a13*a24*a35*a41*c5-a13*a24*a35*c4*a51-a13*a25*a31*c4*a54+a13*a25*a31*a44*c5+a13*a25*c3*a41*a54-a13*a25*c3*a44*a51-a13*a25*a34*a41*c5+a13*a25*a34*c4*a51-a14*a21*c3*a43*a55+a14*a21*c3*a45*a53+a14*a21*a33*c4*a55-a14*a21*a33*a45*c5-a14*a21*a35*c4*a53+a14*a21*a35*a43*c5+a14*c2*a31*a43*a55-a14*c2*a31*a45*a53-a14*c2*a33*a41*a55+a14*c2*a33*a45*a51+a14*c2*a35*a41*a53-a14*c2*a35*a43*a51-a14*a23*a31*c4*a55+a14*a23*a31*a45*c5+a14*a23*c3*a41*a55-a14*a23*c3*a45*a51-a14*a23*a35*a41*c5+a14*a23*a35*c4*a51+a14*a25*a31*c4*a53-a14*a25*a31*a43*c5-a14*a25*c3*a41*a53+a14*a25*c3*a43*a51+a14*a25*a33*a41*c5-a14*a25*a33*c4*a51+a15*a21*c3*a43*a54-a15*a21*c3*a44*a53-a15*a21*a33*c4*a54+a15*a21*a33*a44*c5+a15*a21*a34*c4*a53-a15*a21*a34*a43*c5-a15*c2*a31*a43*a54+a15*c2*a31*a44*a53+a15*c2*a33*a41*a54-a15*c2*a33*a44*a51-a15*c2*a34*a41*a53+a15*c2*a34*a43*a51+a15*a23*a31*c4*a54-a15*a23*a31*a44*c5-a15*a23*c3*a41*a54+a15*a23*c3*a44*a51+a15*a23*a34*a41*c5-a15*a23*a34*c4*a51-a15*a24*a31*c4*a53+a15*a24*a31*a43*c5+a15*a24*c3*a41*a53-a15*a24*c3*a43*a51-a15*a24*a33*a41*c5+a15*a24*a33*c4*a51)/determinant;
	
}

//divergence2�ɂ�����A3����1���ߎ����s���֐��i169�s�j
double calc_WLSM_divu_D3_order1(mpsconfig *CON,vector<mpselastic> &PART,double *matrix, double *B1,double *B2,double *B3,int i,int N)
{
	///�W���s���
	///   ����x2    ����x��y  ����x��z  a = ����x��f  
	///  ����x��y    ����y2   ����y��z  b = ����y��f 
	///  ����x��z   ����y��z   ����z2   c = ����z��f 

	double le=CON->get_distancebp();
	double matrix_val[9];			//matrix�̒l�͈�xgauss()�Ŏg�p����ƒl���ς��̂ŁA�ύX�O��matrix��ۑ�����
	double *weight=new double [PART[i].N];	//���ӗ��q�̏d�݂��i�[����B���ƂŌ덷�]���̂Ƃ��Čv�Z���s�v�ɂȂ�B

	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
			
		double X=(PART[j].r[A_X]-PART[i].r[A_X]);// le�Ŋ���̂͑ł��؂�덷�h�~
		double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
		double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
		double U=(PART[j].u[A_X]-PART[i].u[A_X]);
		double V=(PART[j].u[A_Y]-PART[i].u[A_Y]);
		double W=(PART[j].u[A_Z]-PART[i].u[A_Z]);
		double dis=sqrt(X*X+Y*Y+Z*Z);
					
		double w=1;
		if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);
		weight[k]=w;
					
		matrix[0]+=X*X*w;			//��Xj^2wj
		matrix[1]+=X*Y*w;		//��XjYjwj
		matrix[2]+=X*Z*w;		//��XjZjwj
					
		matrix[4]+=Y*Y*w;			//��Yj^2wj
		matrix[5]+=Y*Z*w;		//��YjZjwj

		matrix[8]+=Z*Z*w;			//��Zj^2wj
			
		B1[0]+=U*X*w;//��fjXjwj
		B1[1]+=U*Y*w;//��fjYjwj
		B1[2]+=U*Z*w;//��fjZjwj

		B2[0]+=V*X*w;//��fjXjwj
		B2[1]+=V*Y*w;//��fjYjwj
		B2[2]+=V*Z*w;//��fjZjwj

		B3[0]+=W*X*w;//��fjXjwj
		B3[1]+=W*Y*w;//��fjYjwj
		B3[2]+=W*Z*w;//��fjZjwj
	}
			
	matrix[3]=matrix[1];		//��XjYjwj
	matrix[6]=matrix[2];		//��XjZjwj
	matrix[7]=matrix[5];		//��YjZjwj

	for(int L=0;L<9;L++) matrix_val[L]=matrix[L];//�s��̒l��ۑ�

	/*double dudx=0;//�������̂ق����኱�����B���ǌ덷�]���������Ȃ�K�E�X
	double dvdy=0;
	double dwdz=0;
	double determinant=(matrix[0]*matrix[4]*matrix[8]-matrix[0]*matrix[5]*matrix[7]-matrix[1]*matrix[3]*matrix[8]+matrix[1]*matrix[5]*matrix[6]+matrix[2]*matrix[3]*matrix[7]-matrix[2]*matrix[4]*matrix[6]);//�s��
			
	dudx=(B1[0]*matrix[4]*matrix[8]-B1[0]*matrix[5]*matrix[7]-matrix[1]*B1[1]*matrix[8]+matrix[1]*matrix[5]*B1[2]+matrix[2]*B1[1]*matrix[7]-matrix[2]*matrix[4]*B1[2])/determinant;
	dvdy=(matrix[0]*B2[1]*matrix[8]-matrix[0]*matrix[5]*B2[2]-B2[0]*matrix[3]*matrix[8]+B2[0]*matrix[5]*matrix[6]+matrix[2]*matrix[3]*B2[2]-matrix[2]*B2[1]*matrix[6])/determinant;
	dwdz=(matrix[0]*matrix[4]*B3[2]-matrix[0]*B3[1]*matrix[7]-matrix[1]*matrix[3]*B3[2]+matrix[1]*B3[1]*matrix[6]+B3[0]*matrix[3]*matrix[7]-B3[0]*matrix[4]*matrix[6])/determinant;
	*/	

	//�s����K�E�X�̏����@�ŉ����@����B�Ɋi�[�����
	int Xflag=OFF;//�K�E�X�̏����@�����邩���Ȃ����B���s�񂪃[���Ȃ�K�E�X�̏����@�����炾�߁B
	int Yflag=OFF;
	int Zflag=OFF;
	for(int k=0;k<3;k++)
	{
		if(B1[k]!=0) Xflag=ON;//B1[k]1�����ׂă[���Ȃ�Xflag��OFF�̂܂܁B
		if(B2[k]!=0) Yflag=ON;//B1[k]1�����ׂă[���Ȃ�Xflag��OFF�̂܂܁B
		if(B3[k]!=0) Zflag=ON;//B1[k]1�����ׂă[���Ȃ�Xflag��OFF�̂܂܁B
	}

	if(Xflag==ON)
	{
		gauss(matrix,B1,N);
		for(int L=0;L<9;L++) matrix[L]=matrix_val[L];//matrix��ύX�O�ɖ߂�
	}
	else B1[0]=0; //flag��OFF�Ȃ�ǂ݂̂�dudx�̓[���B
	
	if(Yflag==ON)
	{
		gauss(matrix,B2,N);
		for(int L=0;L<9;L++) matrix[L]=matrix_val[L];//matrix��ύX�O�ɖ߂�
	}
	else B2[1]=0;	//flag��OFF�Ȃ�ǂ݂̂�dvdy�̓[���B
	
	if(Zflag==ON) gauss(matrix,B3,N);
	else B3[2]=0;	//flag��OFF�Ȃ�ǂ݂̂�dwdz�̓[���B

	//�v�Z�I��

	//�덷�𒲍�
	double Q[3]={0,0,0};//�e�����̌덷
	double W=0;//�d�݂̑��a
	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
		double X=(PART[j].r[A_X]-PART[i].r[A_X]);
		double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
		double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
		double dU=(PART[j].u[A_X]-PART[i].u[A_X]);
		double dV=(PART[j].u[A_Y]-PART[i].u[A_Y]);
		double dW=(PART[j].u[A_Z]-PART[i].u[A_Z]);
		double w=weight[k];
		W+=w;
		double err[3];
		err[A_X]=B1[0]*X+B1[1]*Y+B1[2]*Z-dU;
		err[A_Y]=B2[0]*X+B2[1]*Y+B2[2]*Z-dV;
		err[A_Z]=B3[0]*X+B3[1]*Y+B3[2]*Z-dW;
		Q[A_X]+=err[A_X]*err[A_X]*w;
		Q[A_Y]+=err[A_Y]*err[A_Y]*w;
		Q[A_Z]+=err[A_Z]*err[A_Z]*w;
	}
	for(int D=0;D<3;D++)
	{
		if(W>0) Q[D]/=W;
		Q[D]=sqrt(Q[D]);
		double speed=sqrt(PART[i].u[D]*PART[i].u[D]);
		if(speed>1e-6) Q[D]/=speed;
		else Q[D]=0;
	}
	//if(Q[A_X]>10 ||Q[A_Y]>10||Q[A_Z]>10) cout<<"Q("<<i<<")="<<Q[A_X]<<" "<<Q[A_Y]<<" "<<Q[A_Z]<<endl;
	/////*/

	double dudx=B1[0];
	double dvdy=B2[1];
	double dwdz=B3[2];

	double div=(dudx+dvdy+dwdz);

	delete [] weight;

	return div;
}

//divergence2�ɂ�����A3����1���ߎ����s���֐�ver.2�i256�s�j
double calc_WLSM_divu_D3_order1_2(mpsconfig *CON,vector<mpselastic> &PART,double *matrix, double *B1,double *B2,double *B3,int i,int N)
{
	///�W���s���
	///   ����x2    ����x��y  ����x��z ����x a = ����xfj  
	///  ����x��y    ����y2   ����y��z ����y b = ����yfj 
	///  ����x��z   ����y��z  ����z2  ����z  c = ����zfj 
	///  ����x      ����y     ����z     ��1  d = ��fj

	double le=CON->get_distancebp();
	double matrix_val[16];			//matrix�̒l�͈�xgauss()�Ŏg�p����ƒl���ς��̂ŁA�ύX�O��matrix��ۑ�����
	double *weight=new double [PART[i].N];	//���ӗ��q�̏d�݂��i�[����B���ƂŌ덷�]���̂Ƃ��Čv�Z���s�v�ɂȂ�B

	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
			
		double X=(PART[j].r[A_X]-PART[i].r[A_X]);// le�Ŋ���̂͑ł��؂�덷�h�~
		double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
		double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
		//double U=(PART[j].u[A_X]-PART[i].u[A_X]);
		//double V=(PART[j].u[A_Y]-PART[i].u[A_Y]);
		//double W=(PART[j].u[A_Z]-PART[i].u[A_Z]);
		double U=(PART[j].u[A_X]);
		double V=(PART[j].u[A_Y]);
		double W=(PART[j].u[A_Z]);
		double dis=sqrt(X*X+Y*Y+Z*Z);
					
		double w=1;
		if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);
		weight[k]=w;
					
		matrix[0]+=X*X*w;			//��Xj^2wj
		matrix[1]+=X*Y*w;		//��XjYjwj
		matrix[2]+=X*Z*w;		//��XjZjwj
		matrix[3]+=X*w;
					
		matrix[5]+=Y*Y*w;			//��Yj^2wj
		matrix[6]+=Y*Z*w;		//��YjZjwj
		matrix[7]+=Y*w;

		matrix[10]+=Z*Z*w;			//��Zj^2wj
		matrix[11]+=Z*w;

		matrix[15]+=w;
			
		B1[0]+=U*X*w;//��fjXjwj
		B1[1]+=U*Y*w;//��fjYjwj
		B1[2]+=U*Z*w;//��fjZjwj
		B1[3]+=U*w;//��fjwj

		B2[0]+=V*X*w;//��fjXjwj
		B2[1]+=V*Y*w;//��fjYjwj
		B2[2]+=V*Z*w;//��fjZjwj
		B2[3]+=V*w;//��fjwj

		B3[0]+=W*X*w;//��fjXjwj
		B3[1]+=W*Y*w;//��fjYjwj
		B3[2]+=W*Z*w;//��fjZjwj
		B3[3]+=W*w;//��fjwj
	}
			
	matrix[4]=matrix[1];	
	matrix[8]=matrix[2];		
	matrix[9]=matrix[6];
	matrix[12]=matrix[3];
	matrix[13]=matrix[7];
	matrix[14]=matrix[11];

	matrix[15]+=1;//�������g
	B1[3]+=PART[i].u[A_X];
	B2[3]+=PART[i].u[A_Y];
	B3[3]+=PART[i].u[A_Z];

	for(int L=0;L<16;L++) matrix_val[L]=matrix[L];//�s��̒l��ۑ�

	//�s����K�E�X�̏����@�ŉ����@����B�Ɋi�[�����
	int Xflag=OFF;//�K�E�X�̏����@�����邩���Ȃ����B���s�񂪃[���Ȃ�K�E�X�̏����@�����炾�߁B
	int Yflag=OFF;
	int Zflag=OFF;
	for(int k=0;k<4;k++)
	{
		if(B1[k]!=0) Xflag=ON;//B1[k]1�����ׂă[���Ȃ�Xflag��OFF�̂܂܁B
		if(B2[k]!=0) Yflag=ON;//B2[k]1�����ׂă[���Ȃ�Xflag��OFF�̂܂܁B
		if(B3[k]!=0) Zflag=ON;//B3[k]1�����ׂă[���Ȃ�Xflag��OFF�̂܂܁B
	}

	if(Xflag==ON)
	{
		gauss(matrix,B1,N);
		for(int L=0;L<16;L++) matrix[L]=matrix_val[L];//matrix��ύX�O�ɖ߂�
	}
	else for(int k=0;k<4;k++) B1[k]=0; //flag��OFF�Ȃ�ǂ݂̂�dudx�̓[���B
	
	if(Yflag==ON)
	{
		gauss(matrix,B2,N);
		for(int L=0;L<16;L++) matrix[L]=matrix_val[L];//matrix��ύX�O�ɖ߂�
	}
	else for(int k=0;k<4;k++) B2[k]=0;	//flag��OFF�Ȃ�ǂ݂̂�dvdy�̓[���B
	
	if(Zflag==ON) gauss(matrix,B3,N);
	else for(int k=0;k<4;k++) B3[k]=0;	//flag��OFF�Ȃ�ǂ݂̂�dwdz�̓[���B

	//�v�Z�I��

	//�덷�𒲍�
	double Q[3]={0,0,0};//�e�����̌덷
	double err[3];
	double W=1;//�d�݂̑��a
	err[A_X]=B1[3]-PART[i].u[A_X];//���g�̌덷
	err[A_Y]=B2[3]-PART[i].u[A_Y];
	err[A_Z]=B3[3]-PART[i].u[A_Z];
	Q[A_X]+=err[A_X]*err[A_X];
	Q[A_Y]+=err[A_Y]*err[A_Y];
	Q[A_Z]+=err[A_Z]*err[A_Z];
	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
		double X=(PART[j].r[A_X]-PART[i].r[A_X]);
		double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
		double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
		double U=(PART[j].u[A_X]);
		double V=(PART[j].u[A_Y]);
		double W=(PART[j].u[A_Z]);
		double w=weight[k];
		W+=w;
		
		err[A_X]=B1[0]*X+B1[1]*Y+B1[2]*Z+B1[3]-U;
		err[A_Y]=B2[0]*X+B2[1]*Y+B2[2]*Z+B2[3]-V;
		err[A_Z]=B3[0]*X+B3[1]*Y+B3[2]*Z+B3[3]-W;
		Q[A_X]+=err[A_X]*err[A_X]*w;
		Q[A_Y]+=err[A_Y]*err[A_Y]*w;
		Q[A_Z]+=err[A_Z]*err[A_Z]*w;
	}
	for(int D=0;D<3;D++)
	{
		if(W>0) Q[D]/=W;
		Q[D]=sqrt(Q[D]);
		double speed=sqrt(PART[i].u[D]*PART[i].u[D]);
		if(speed>1e-6) Q[D]/=speed;
		else Q[D]=0;
	}
	//if(Q[A_X]>10 ||Q[A_Y]>10||Q[A_Z]>10) cout<<"Q("<<i<<")="<<Q[A_X]<<" "<<Q[A_Y]<<" "<<Q[A_Z]<<endl;
	/////*/

	double dudx=B1[0];
	double dvdy=B2[1];
	double dwdz=B3[2];

	double div=(dudx+dvdy+dwdz);

	delete [] weight;

	return div;
}

//divergence2�ɂ�����A3����2���ߎ����s���֐��i227�s�j
double calc_WLSM_divu_D3_order2(mpsconfig *CON,vector<mpselastic> &PART,double *matrix, double *B1,double *B2,double *B3,int i,int N)
{
	double le=CON->get_distancebp();
	double matrix_val[81];			//matrix�̒l�͈�xgauss()�Ŏg�p����ƒl���ς��̂ŁA�ύX�O��matrix��ۑ�����
	double *weight=new double [PART[i].N];	//���ӗ��q�̏d�݂��i�[����B���ƂŌ덷�]���̂Ƃ��Čv�Z���s�v�ɂȂ�B

	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
		
		double X=(PART[j].r[A_X]-PART[i].r[A_X]);
		double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
		double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
		double U=(PART[j].u[A_X]-PART[i].u[A_X]);
		double V=(PART[j].u[A_Y]-PART[i].u[A_Y]);
		double W=(PART[j].u[A_Z]-PART[i].u[A_Z]);
		double dis=sqrt(X*X+Y*Y+Z*Z);
						
		double w=1;
		if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);//���̎�����dis=R�̂Ƃ��d�݂��[���ɋ߂��Ȃ�
		weight[k]=w;

		matrix[0]+=X*X*w;		//��Xj^2wj
		matrix[1]+=X*Y*w;		//��XjYjwj
		matrix[2]+=X*Z*w;		//��XjZjwj
		matrix[3]+=X*X*X*w;		//��Xj^3wj
		matrix[4]+=X*Y*Y*w;		//��XjYj^2wj
		matrix[5]+=X*Z*Z*w;		//��XjZj^2wj
		matrix[6]+=X*X*Y*w;		//��Xj^2Yjwj
		matrix[7]+=X*Y*Z*w;		//��XjYjZjwj
		matrix[8]+=X*X*Z*w;		//��Xj^2Zjwj
	
		matrix[10]+=Y*Y*w;		
		matrix[11]+=Y*Z*w;		
		matrix[12]+=X*X*Y*w;		
		matrix[13]+=Y*Y*Y*w;		
		matrix[14]+=Y*Z*Z*w;
		matrix[15]+=X*Y*Y*w;
		matrix[16]+=Y*Y*Z*w;
						
		matrix[20]+=Z*Z*w;			
		matrix[23]+=Z*Z*Z*w;		
					
		matrix[30]+=X*X*X*X*w;
		matrix[31]+=X*X*Y*Y*w;
		matrix[32]+=X*X*Z*Z*w;	
		matrix[33]+=X*X*X*Y*w;	
		matrix[34]+=X*X*Y*Z*w;	
		matrix[35]+=X*X*X*Z*w;	
					
		matrix[40]+=Y*Y*Y*Y*w;
		matrix[41]+=Y*Y*Z*Z*w;
		matrix[42]+=X*Y*Y*Y*w;
		matrix[43]+=Y*Y*Y*Z*w;
		matrix[44]+=X*Y*Y*Z*w;

		matrix[50]+=Z*Z*Z*Z*w;	//6�s��
		matrix[51]+=X*Y*Z*Z*w;
		matrix[52]+=Y*Z*Z*Z*w;
		matrix[53]+=X*Z*Z*Z*w;

		//7�`9�s�ڂ͂��ׂĊ����̗v�f����]�p���\


		B1[0]+=U*X*w;		//a
		B1[1]+=U*Y*w;		//b
		B1[2]+=U*Z*w;		//c
		B1[3]+=U*X*X*w;		//d
		B1[4]+=U*Y*Y*w;		//e
		B1[5]+=U*Z*Z*w;		//f
		B1[6]+=U*X*Y*w;		//g
		B1[7]+=U*Y*Z*w;		//h
		B1[8]+=U*X*Z*w;		//i

		B2[0]+=V*X*w;		//a
		B2[1]+=V*Y*w;		//b
		B2[2]+=V*Z*w;		//c
		B2[3]+=V*X*X*w;		//d
		B2[4]+=V*Y*Y*w;		//e
		B2[5]+=V*Z*Z*w;		//f
		B2[6]+=V*X*Y*w;		//g
		B2[7]+=V*Y*Z*w;		//h
		B2[8]+=V*X*Z*w;		//i

		B3[0]+=W*X*w;		//a
		B3[1]+=W*Y*w;		//b
		B3[2]+=W*Z*w;		//c
		B3[3]+=W*X*X*w;		//d
		B3[4]+=W*Y*Y*w;		//e
		B3[5]+=W*Z*Z*w;		//f
		B3[6]+=W*X*Y*w;		//g
		B3[7]+=W*Y*Z*w;		//h
		B3[8]+=W*X*Z*w;		//i
		
	}
	matrix[9]=matrix[1];		//��XjYjwj
	matrix[17]=matrix[7];

	matrix[18]=matrix[2];
	matrix[19]=matrix[11];
	matrix[21]=matrix[8];
	matrix[22]=matrix[16];
	matrix[24]=matrix[7];
	matrix[25]=matrix[14];
	matrix[26]=matrix[5];

	matrix[27]=matrix[3];
	matrix[28]=matrix[12];
	matrix[29]=matrix[21];

	matrix[36]=matrix[4];
	matrix[37]=matrix[13];
	matrix[38]=matrix[22];
	matrix[39]=matrix[31];

	matrix[45]=matrix[5];
	matrix[46]=matrix[14];
	matrix[47]=matrix[23];
	matrix[48]=matrix[32];
	matrix[49]=matrix[41];

	matrix[54]=matrix[6];
	matrix[55]=matrix[15];
	matrix[56]=matrix[24];
	matrix[57]=matrix[33];
	matrix[58]=matrix[42];
	matrix[59]=matrix[51];
	matrix[60]=matrix[31];
	matrix[61]=matrix[44];
	matrix[62]=matrix[34];

	matrix[63]=matrix[7];
	matrix[64]=matrix[16];
	matrix[65]=matrix[25];
	matrix[66]=matrix[34];
	matrix[67]=matrix[43];
	matrix[68]=matrix[52];
	matrix[69]=matrix[61];
	matrix[70]=matrix[41];
	matrix[71]=matrix[51];

	matrix[72]=matrix[8];
	matrix[73]=matrix[17];
	matrix[74]=matrix[26];
	matrix[75]=matrix[35];
	matrix[76]=matrix[44];
	matrix[77]=matrix[53];
	matrix[78]=matrix[62];
	matrix[79]=matrix[71];
	matrix[80]=matrix[32];

	for(int L=0;L<81;L++) matrix_val[L]=matrix[L];//�s��̒l��ۑ�
			
	//�s����K�E�X�̏����@�ŉ����@����B�Ɋi�[�����
	int Xflag=OFF;
	int Yflag=OFF;//�K�E�X�̏����@�����邩���Ȃ����B���s�񂪃[���Ȃ�K�E�X�̏����@�����炾�߁B
	int Zflag=OFF;
	for(int k=0;k<9;k++)
	{
		if(B1[k]!=0) Xflag=ON;//B1[k]1�����ׂă[���Ȃ�Xflag��OFF�̂܂܁B
		if(B2[k]!=0) Yflag=ON;//B1[k]1�����ׂă[���Ȃ�Xflag��OFF�̂܂܁B
		if(B3[k]!=0) Zflag=ON;//B1[k]1�����ׂă[���Ȃ�Xflag��OFF�̂܂܁B
	}

	if(Xflag==ON)
	{
		gauss(matrix,B1,N);
		for(int L=0;L<81;L++) matrix[L]=matrix_val[L];//�s��̒l�߂�
	}
	else {for(int k=0;k<N;k++) B1[k]=0;} //flag��OFF�Ȃ�ǂ݂̂�dudx�̓[���B

	if(Yflag==ON) 
	{
		gauss(matrix,B2,N);
		for(int L=0;L<81;L++) matrix[L]=matrix_val[L];//�s��̒l�߂�
	}
	else {for(int k=0;k<N;k++) B2[k]=0;}	//flag��OFF�Ȃ�ǂ݂̂�dvdy�̓[���B

	if(Zflag==ON) gauss(matrix,B3,N);
	else {for(int k=0;k<N;k++) B3[k]=0;}	//flag��OFF�Ȃ�ǂ݂̂�dwdz�̓[���B

	////�v�Z�I��

	//�덷�𒲍�
	double Q[3]={0,0,0};//�e�����̕W���΍�
	double W=0;		//�d�݂̑��a
	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
		double X=(PART[j].r[A_X]-PART[i].r[A_X]);
		double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
		double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
		double dU=(PART[j].u[A_X]-PART[i].u[A_X]);
		double dV=(PART[j].u[A_Y]-PART[i].u[A_Y]);
		double dW=(PART[j].u[A_Z]-PART[i].u[A_Z]);
		double w=weight[k];
		W+=w;
		double err[3];
		err[A_X]=B1[0]*X+B1[1]*Y+B1[2]*Z+B1[3]*X*X+B1[4]*Y*Y+B1[5]*Z*Z+B1[6]*X*Y+B1[7]*Y*Z+B1[8]*X*Z-dU;
		err[A_Y]=B2[0]*X+B2[1]*Y+B2[2]*Z+B2[3]*X*X+B2[4]*Y*Y+B2[5]*Z*Z+B2[6]*X*Y+B2[7]*Y*Z+B2[8]*X*Z-dV;
		err[A_Z]=B3[0]*X+B3[1]*Y+B3[2]*Z+B3[3]*X*X+B3[4]*Y*Y+B3[5]*Z*Z+B3[6]*X*Y+B3[7]*Y*Z+B3[8]*X*Z-dW;
		Q[A_X]+=err[A_X]*err[A_X]*w;
		Q[A_Y]+=err[A_Y]*err[A_Y]*w;
		Q[A_Z]+=err[A_Z]*err[A_Z]*w;
	}
	for(int D=0;D<3;D++)
	{
		if(W>0) Q[D]/=W;
		Q[D]=sqrt(Q[D]);
		double speed=sqrt(PART[i].u[D]*PART[i].u[D]);
		if(speed>1e-6) Q[D]/=speed;
		else Q[D]=0;
	}
	//if(Q[A_X]>10 ||Q[A_Y]>10||Q[A_Z]>10) cout<<"Q("<<i<<")="<<Q[A_X]<<" "<<Q[A_Y]<<" "<<Q[A_Z]<<endl;
	//////


	double dudx=B1[0];
	double dvdy=B2[1];
	double dwdz=B3[2];
			
	double div=dudx+dvdy+dwdz;

	delete [] weight;

	return div;
}

//divergence2�ɂ�����A3����2���ߎ����s���֐�ver.2 ���m�����ЂƂ����i233�s�j
double calc_WLSM_divu_D3_order2_2(mpsconfig *CON,vector<mpselastic> &PART,double *matrix, double *B1,double *B2,double *B3,int i,int N)
{
	//P=a��x+b��y+c��z+d��x2+e��y2+f��z2+g��x��y+h��y��z+i��z��x+P�Ƃ����ƁA
	///�W���s���
	///   ����x2      ����x��y    ����x��z    ����x3       ����x��y2    ����x��z2    ����x2��y     ����x��y��z  ����x2��z    ����x    a = ����xf  
	///   ����x��y    ����y2      ����y��z    ����x2��y    ����y3       ����y��z2    ����x��y2     ����y2��z    ����x��y��z  ����y    b = ����yf
	///   ����x��z    ����y��z    ����z2      ����x2��z    ����y2��z    ����z3       ����x��y��z   ����y��z2    ����x��z2    ����z    c = ����zf
	///   ����x3      ����x2��y   ����x2��z   ����x4       ����x2��y2   ����x2��z2   ����x3��y     ����x2��y��z ����x3��z    ����x2   d = ����x2f
	///   ����x��y2   ����y3      ����y2��z   ����x2��y2   ����y4       ����y2��z2   ����x��y3     ����y3��z    ����x��y2��z ����y2   e = ����y2f
	///   ����x��z2   ����y��z2   ����z3      ����x2��z2   ����y2��z2   ����z4       ����x��y��z2  ����y��z3    ����x��z3	 ����z2   f = ����z2f
	///   ����x2��y   ����x��y2   ����x��y��z ����x3��y    ����x��y3    ����x��y��z2 ����x2��y2    ����x��y2��z ����x2��y��z ����x��y g = ����x��yf
	///   ����x��y��z ����y2��z   ����y��z2   ����x2��y��z ����y3��z    ����y��z3    ����x��y2��z  ����y2��z2   ����x��y��z2 ����y��z h = ����y��zf
	///   ����x2��z   ����x��y��z ����x��z2   ����x3��z    ����x��y2��z  ����x��z3   ����x2��y ��z ����x��y��z2 ����x2��z2   ����z��x i = ����x��zf
	///   ����x       ����y ����z ����x2      ����y2       ����z        ����x��y     ����y��z      ����x��z     ����z��x     ��1      P = ��fj

	double le=CON->get_distancebp();
	double matrix_val[100];			//matrix�̒l�͈�xgauss()�Ŏg�p����ƒl���ς��̂ŁA�ύX�O��matrix��ۑ�����
	double *weight=new double [PART[i].N];	//���ӗ��q�̏d�݂��i�[����B���ƂŌ덷�]���̂Ƃ��Čv�Z���s�v�ɂȂ�B

	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
		
		double X=(PART[j].r[A_X]-PART[i].r[A_X]);
		double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
		double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
		double U=(PART[j].u[A_X]-PART[i].u[A_X]);
		double V=(PART[j].u[A_Y]-PART[i].u[A_Y]);
		double W=(PART[j].u[A_Z]-PART[i].u[A_Z]);
		double dis=sqrt(X*X+Y*Y+Z*Z);
						
		double w=1;
		if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);//���̎�����dis=R�̂Ƃ��d�݂��[���ɋ߂��Ȃ�
		weight[k]=w;

		matrix[0]+=X*X*w;		//��Xj^2wj
		matrix[1]+=X*Y*w;		//��XjYjwj
		matrix[2]+=X*Z*w;		//��XjZjwj
		matrix[3]+=X*X*X*w;		//��Xj^3wj
		matrix[4]+=X*Y*Y*w;		//��XjYj^2wj
		matrix[5]+=X*Z*Z*w;		//��XjZj^2wj
		matrix[6]+=X*X*Y*w;		//��Xj^2Yjwj
		matrix[7]+=X*Y*Z*w;		//��XjYjZjwj
		matrix[8]+=X*X*Z*w;		//��Xj^2Zjwj
		matrix[9]+=X*w;		//��Xj^2Zjwj
	
		matrix[11]+=Y*Y*w;		
		matrix[12]+=Y*Z*w;		
		matrix[13]+=X*X*Y*w;		
		matrix[14]+=Y*Y*Y*w;		
		matrix[15]+=Y*Z*Z*w;
		matrix[16]+=X*Y*Y*w;
		matrix[17]+=Y*Y*Z*w;
		matrix[19]+=Y*w;
					
		matrix[22]+=Z*Z*w;			
		matrix[23]+=X*X*Z*w;
		matrix[24]+=Y*Y*Z*w;
		matrix[25]+=Y*Y*Y*w;
		matrix[29]+=Z*w;
					
		matrix[33]+=X*X*X*X*w;
		matrix[34]+=X*X*Y*Y*w;
		matrix[35]+=X*X*Z*Z*w;	
		matrix[36]+=X*X*X*Y*w;	
		matrix[37]+=X*X*Y*Z*w;	
		matrix[38]+=X*X*X*Z*w;	
					
		matrix[44]+=Y*Y*Y*Y*w;
		matrix[45]+=Y*Y*Z*Z*w;
		matrix[46]+=X*Y*Y*Y*w;
		matrix[47]+=Y*Y*Y*Z*w;
		matrix[48]+=X*Y*Y*Z*w;

		matrix[55]+=Z*Z*Z*Z*w;	//6�s��
		matrix[56]+=X*Y*Z*Z*w;
		matrix[57]+=Y*Z*Z*Z*w;
		matrix[58]+=X*Z*Z*Z*w;

		matrix[99]+=w;
		//7�`9�s�ڂ͂��ׂĊ����̗v�f����]�p���\

		B1[0]+=U*X*w;		//a
		B1[1]+=U*Y*w;		//b
		B1[2]+=U*Z*w;		//c
		B1[3]+=U*X*X*w;		//d
		B1[4]+=U*Y*Y*w;		//e
		B1[5]+=U*Z*Z*w;		//f
		B1[6]+=U*X*Y*w;		//g
		B1[7]+=U*Y*Z*w;		//h
		B1[8]+=U*X*Z*w;		//i
		B1[9]+=U*w;

		B2[0]+=V*X*w;		//a
		B2[1]+=V*Y*w;		//b
		B2[2]+=V*Z*w;		//c
		B2[3]+=V*X*X*w;		//d
		B2[4]+=V*Y*Y*w;		//e
		B2[5]+=V*Z*Z*w;		//f
		B2[6]+=V*X*Y*w;		//g
		B2[7]+=V*Y*Z*w;		//h
		B2[8]+=V*X*Z*w;		//i
		B2[9]+=V*w;	

		B3[0]+=W*X*w;		//a
		B3[1]+=W*Y*w;		//b
		B3[2]+=W*Z*w;		//c
		B3[3]+=W*X*X*w;		//d
		B3[4]+=W*Y*Y*w;		//e
		B3[5]+=W*Z*Z*w;		//f
		B3[6]+=W*X*Y*w;		//g
		B3[7]+=W*Y*Z*w;		//h
		B3[8]+=W*X*Z*w;		//i
		B3[9]+=W*w;
	}
	matrix[10]=matrix[1];
	matrix[18]=matrix[7];

	matrix[20]=matrix[2];
	matrix[21]=matrix[12];
	matrix[24]=matrix[16];
	matrix[26]=matrix[7];
	matrix[27]=matrix[15];
	matrix[28]=matrix[5];

	for(int k=0;k<=2;k++) matrix[30+k]=matrix[3+10*k];//30�`32�v�f
	matrix[39]=matrix[0];

	for(int k=0;k<=3;k++) matrix[40+k]=matrix[4+10*k];//40�`43�v�f
	matrix[49]=matrix[11];

	for(int k=0;k<=4;k++) matrix[50+k]=matrix[5+10*k];//50�`54�v�f
	matrix[59]=matrix[22];

	for(int k=0;k<=5;k++) matrix[60+k]=matrix[6+10*k];//60�`65�v�f
	matrix[66]=matrix[34];
	matrix[67]=matrix[48];
	matrix[68]=matrix[37];
	matrix[69]=matrix[1];

	for(int k=0;k<=6;k++) matrix[70+k]=matrix[7+10*k];//70�`76�v�f
	matrix[77]=matrix[54];
	matrix[78]=matrix[56];
	matrix[79]=matrix[12];
	
	for(int k=0;k<=7;k++) matrix[80+k]=matrix[8+10*k];//80�`87�v�f
	matrix[88]=matrix[35];
	matrix[89]=matrix[20];

	for(int k=0;k<=8;k++) matrix[90+k]=matrix[9+10*k];//90�`98�v�f

	matrix[99]+=1;//���g
	B1[9]+=PART[i].u[A_X];
	B2[9]+=PART[i].u[A_Y];
	B3[9]+=PART[i].u[A_Z];

	for(int L=0;L<100;L++) matrix_val[L]=matrix[L];//�s��̒l��ۑ�
			
	//�s����K�E�X�̏����@�ŉ����@����B�Ɋi�[�����
	int Xflag=OFF;
	int Yflag=OFF;//�K�E�X�̏����@�����邩���Ȃ����B���s�񂪃[���Ȃ�K�E�X�̏����@�����炾�߁B
	int Zflag=OFF;
	for(int k=0;k<10;k++)
	{
		if(B1[k]!=0) Xflag=ON;//B1[k]1�����ׂă[���Ȃ�Xflag��OFF�̂܂܁B
		if(B2[k]!=0) Yflag=ON;//B1[k]1�����ׂă[���Ȃ�Xflag��OFF�̂܂܁B
		if(B3[k]!=0) Zflag=ON;//B1[k]1�����ׂă[���Ȃ�Xflag��OFF�̂܂܁B
	}

	if(Xflag==ON)
	{
		gauss(matrix,B1,N);
		for(int L=0;L<100;L++) matrix[L]=matrix_val[L];//�s��̒l�߂�
	}
	else {for(int k=0;k<N;k++) B1[k]=0;} //flag��OFF�Ȃ�ǂ݂̂�dudx�̓[���B

	if(Yflag==ON) 
	{
		gauss(matrix,B2,N);
		for(int L=0;L<100;L++) matrix[L]=matrix_val[L];//�s��̒l�߂�
	}
	else {for(int k=0;k<N;k++) B2[k]=0;}	//flag��OFF�Ȃ�ǂ݂̂�dvdy�̓[���B

	if(Zflag==ON) gauss(matrix,B3,N);
	else {for(int k=0;k<N;k++) B3[k]=0;}	//flag��OFF�Ȃ�ǂ݂̂�dwdz�̓[���B

	////�v�Z�I��

	//�덷�𒲍�
	double Q[3]={0,0,0};//�e�����̕W���΍�
	double W=0;		//�d�݂̑��a
	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
		double X=(PART[j].r[A_X]-PART[i].r[A_X]);
		double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
		double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
		double dU=(PART[j].u[A_X]-PART[i].u[A_X]);
		double dV=(PART[j].u[A_Y]-PART[i].u[A_Y]);
		double dW=(PART[j].u[A_Z]-PART[i].u[A_Z]);
		double w=weight[k];
		W+=w;
		double err[3];
		err[A_X]=B1[0]*X+B1[1]*Y+B1[2]*Z+B1[3]*X*X+B1[4]*Y*Y+B1[5]*Z*Z+B1[6]*X*Y+B1[7]*Y*Z+B1[8]*X*Z-dU;
		err[A_Y]=B2[0]*X+B2[1]*Y+B2[2]*Z+B2[3]*X*X+B2[4]*Y*Y+B2[5]*Z*Z+B2[6]*X*Y+B2[7]*Y*Z+B2[8]*X*Z-dV;
		err[A_Z]=B3[0]*X+B3[1]*Y+B3[2]*Z+B3[3]*X*X+B3[4]*Y*Y+B3[5]*Z*Z+B3[6]*X*Y+B3[7]*Y*Z+B3[8]*X*Z-dW;
		Q[A_X]+=err[A_X]*err[A_X]*w;
		Q[A_Y]+=err[A_Y]*err[A_Y]*w;
		Q[A_Z]+=err[A_Z]*err[A_Z]*w;
	}
	for(int D=0;D<3;D++)
	{
		if(W>0) Q[D]/=W;
		Q[D]=sqrt(Q[D]);
		double speed=sqrt(PART[i].u[D]*PART[i].u[D]);
		if(speed>1e-6) Q[D]/=speed;
		else Q[D]=0;
	}
	//if(Q[A_X]>10 ||Q[A_Y]>10||Q[A_Z]>10) cout<<"Q("<<i<<")="<<Q[A_X]<<" "<<Q[A_Y]<<" "<<Q[A_Z]<<endl;
	//////


	double dudx=B1[0];
	double dvdy=B2[1];
	double dwdz=B3[2];
			
	double div=dudx+dvdy+dwdz;

	delete [] weight;

	return div;
}

//���q����͗̈�̊O�ɂłĂ��Ȃ����`�F�b�N�i70�s�j
int check_position(mpsconfig *CON,vector<mpselastic> &PART,int fluid_number,int *particle_number)
{
	int sw=OFF;	//�{�֐��ŕԂ��l�BOFF�Ȃ痱�q���ɕω��Ȃ��BON�Ȃ�ω������Ƃ�����
	int num=0;	//���ł��闱�q��
	double le=CON->get_distancebp();
	int *flag=new int [fluid_number];	//ON�Ȃ�̈�� OFF�Ȃ�̈�O

	double Xmax=CON->get_maxX(); double Xmin=CON->get_minX();
	double Ymax=CON->get_maxY(); double Ymin=CON->get_minY();
	double Zmax=CON->get_maxZ(); double Zmin=CON->get_minZ();

	double dx=CON->get_dx()*le;	//�i�q��

	Xmax-=dx; Ymax-=dx; Zmax-=2*dx;		//�ی���������1�i�q�������ɋ��E���Ƃ�B������O���Ȃ痱�q������
	Xmin+=dx; Ymin+=dx; Zmin+=2*dx;

	vector<mpselastic>::iterator p,p0;//�����q
	p0=PART.begin();
	
	for(int i=0;i<fluid_number;i++) 
	{
		flag[i]=ON;
		if(PART[i].r[A_X]<Xmin || PART[i].r[A_X]>Xmax) flag[i]=OFF;
		else if(PART[i].r[A_Y]<Ymin || PART[i].r[A_Y]>Ymax) flag[i]=OFF;
		else if(PART[i].r[A_Z]<Zmin || PART[i].r[A_Z]>Zmax) flag[i]=OFF;
	}//flag[i]�����܂���
	
	for(int i=0;i<fluid_number;i++) if(PART[i].N<=4) flag[i]=OFF;//���ӗ��q�������Ȃ����q���폜?

	for(int i=0;i<fluid_number;i++) if(flag[i]==OFF) num++;


	if(num>0)//�̈�O���q�����m�����Ȃ�
	{
		sw=ON;
		int *erase_id=new int[num];
		int count=0;
		for(int i=0;i<fluid_number;i++)
		{
			if(flag[i]==OFF)
			{
				erase_id[count]=i;//�����ׂ�id���L��
				count++;
			}
			
		}
	
		for(int i=0;i<num;i++)
		{
			p=PART.begin();
			p+=erase_id[i];
			
			PART.erase(p);
			for(int j=i+1;j<num;j++)
			{
				erase_id[j]=erase_id[j]-1;//1�l��������
			}
		}

		delete [] erase_id;

		//id�����ꂽ����߂�
		for(int i=0;i<(int) PART.size();i++) if(PART[i].ID!=i) PART[i].ID=i; 
	}
	
	if(sw==ON)
	{
		cout<<"�̈�O���q��T�m "<<num<<"�̗��q������--���݂̗��q��="<<PART.size()<<endl;
		*particle_number=*particle_number-num;
	}
	
	delete [] flag;

	return sw;
}

//���q���߂Â�������̂�h���֐��i63�s�j
void modify_position(mpsconfig *CON,vector<mpselastic> &PART,int fluid_number,double dt)
{
	double le=CON->get_distancebp();

	for(int i=0;i<fluid_number;i++)
	{
		//if(PART[i].surface==ON)
		{
			double mindis=le;
			//mindis=100;			//�����̒Z�����̂��������邾���łȂ��A�������̂�le�ɂ������Ƃ��͂�����
			int J=i;			//�ŋߐڗ��q
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				//if(PART[j].type==FLUID)// && j>i)
				{
					//if(PART[j].surface==ON)
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
						if(dis<mindis)
						{
							mindis=dis;
							J=j;
						}
					}
				}
			}
			if(J!=i && (PART[J].type==FLUID || PART[J].type==ELASTIC || PART[J].type==MAGELAST || PART[J].type==MAGELAST2))//le���ߐڂ��Ă��闬�̗��q���������Ȃ�
			{
				double L=le-mindis;//�J���ׂ�����
				double dL[DIMENSION];
				for(int D=0;D<DIMENSION;D++) dL[D]=PART[J].r[D]-PART[i].r[D];

				for(int D=0;D<DIMENSION;D++)
				{
					double dU=0.5*L/dt;	//�ω����ׂ����x
				//	PART[J].u[D]+=dL[D]/mindis*dU;
					PART[J].r[D]+=dL[D]/mindis*dU*dt;

				//	PART[i].u[D]-=dL[D]/mindis*dU;
					PART[i].r[D]-=dL[D]/mindis*dU*dt;
				}				
			}
			else if(J!=i && (PART[J].type!=FLUID || PART[J].type!=ELASTIC || PART[J].type!=MAGELAST || PART[J].type!=MAGELAST2))//le���ߐڂ��Ă���Ǘ��q���������Ȃ�
			{
				double L=le-mindis;//�J���ׂ�����
				double dL[DIMENSION];
				for(int D=0;D<DIMENSION;D++) dL[D]=PART[J].r[D]-PART[i].r[D];

				for(int D=0;D<DIMENSION;D++)
				{
					double dU=L/dt;	//�ω����ׂ����x
				//	PART[i].u[D]-=dL[D]/mindis*dU;
					PART[i].r[D]-=dL[D]/mindis*dU*dt;
				}				
			}
		}
	}
}


//���q�̐όv�Z�֐��i18�s�j
double get_volume(mpsconfig *CON)
{
	double V=0;//�̐�
	double le=CON->get_distancebp();
	if(CON->get_model_set_way()==0)	//�����i�q�̂Ƃ�
	{
		if(CON->get_dimension()==2){V=le*le;}
		else	{V=le*le*le;}
	}
	else if(CON->get_model_set_way()==1)	//�ז��i�q�̂Ƃ�
	{
		if(CON->get_dimension()==2){V=sqrt(3.0)/2*le*le;}
		else V=le*le*le/sqrt(2.0);
	}	
	else cout<<"���f���̐ςݕ����s��ł� �̐ς��v�Z�ł��܂���"<<endl;
	return V;
}

//���q�����x�̂��ꑪ��֐��i18�s�j
void output_particle_density(mpsconfig *CON,vector<mpselastic> &PART,int fluid_number,double n0,int particle_number,int t)
{
	///�������q�����x���z���o�͂���

	ofstream fp("initial_n0.dat");
	double le=CON->get_distancebp();

	if(CON->get_dimension()==2) {for(int i=0;i<particle_number;i++) fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].PND<<endl;}
	else if(CON->get_dimension()==3)
	{
		for(int i=0;i<particle_number;i++)
		{
			if(PART[i].r[A_Y]>-0.5*le && PART[i].r[A_Y]<0.5*le) fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<PART[i].PND<<endl;
		}
	}
	fp.close();
}


//�����z�u�Ɖ�͗̈�̊֌W�����֐��i16�s�j
void check_initial_position(mpsconfig *CON,vector<mpselastic> &PART)
{
	double region[3][2];
	region[A_X][0]=CON->get_minX(); region[A_X][1]=CON->get_maxX(); 
	region[A_Y][0]=CON->get_minY(); region[A_Y][1]=CON->get_maxY(); 
	region[A_Z][0]=CON->get_minZ(); region[A_Z][1]=CON->get_maxZ(); 

	for(int i=0;i<PART.size();i++) 
	{
		int flag=OFF;
		for(int D=0;D<3;D++) if(PART[i].r[D]>=region[D][1] || PART[i].r[D]<=region[D][0]) flag=ON;
		
		if(flag==ON)	cout<<"error:��͗̈�O�ɏ������q(id="<<i<<")�����݂��܂� x,y,z="<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
	}
}


//25�s
void FEM3D(mpsconfig &CON, vector<mpselastic> &PART, double **F, int *static_node, int *static_nelm, vector<point3D> &static_NODE, vector<element3D> &static_ELEM, int t, double TIME, int fluid_number)
{
	size_t particle_number=PART.size();
	double dt=CON.get_dt();
	if(CON.get_EM_interval()==1 || t==1 || (t-1)%CON.get_EM_interval()==0){
		FEM3D_calculation(CON, static_node, static_nelm, static_NODE, static_ELEM, t, TIME, PART, fluid_number, particle_number, dt, F);
	}else{ //���X�e�b�v���Ɉ�x�����d����v�Z���s���ꍇ
	
		if(CON.get_dir_for_P()!=2 && CON.get_dir_for_P()!=3)//�d���͂����̓f�B���N���l�̂Ƃ��͂����ł͌v�Z���Ȃ�
		{
			double F_val;
			ifstream fin5("FEM_interval.dat");
			if(!fin5) cout<<"cannot open FEM_interval.dat"<<endl;
			fin5.unsetf(ifstream::dec);
			fin5.setf(ifstream::skipws);
			////�d���͓ǂݎ��
			for(int i=0;i<PART.size();i++)
			{
				if(PART[i].type!=WALL){
				for(int D=0;D<3;D++)
				{
					fin5>>F_val;
					F[D][i]+=F_val;//F[D][i]�͐錾��������Ƀ[���Z�b�g���Ă���
				}
			}
			}
			fin5.close();
		}
	}
}

/*�Q�l
///�d���͌v�Z�֐�
void calc_electro_magnetic_force(mpsconfig &CON,vector<mpselastic> &PART,int fluid_number,double dt,int t,int particle_number,int *INDEX,int **MESH,double n0,double TIME)
{
	/////////�d������
	if(CON.get_mesher()==0)	//�g�삳��쐬���b�V���[
	{
		if(CON.get_dimention()==2)
		{
			if(CON.get_current_step()==1 || CON.get_current_step()%CON.get_FEM_interval()==0)
			{
				if((t==1 && CON.get_FEM_first_step()==ON) || CON.get_current_step()==1 || CON.get_current_step()%CON.get_FEM_interval()==0)
				{
					//�d���͏�����
					for(int i=0;i<particle_number;i++)	for(int D=0;D<CON.get_dimention();D++)	PART[i].eforce[D]=0;
					////�܂�output_node_data�֐��ŗ��̏���FEM�p�ɏo�͂���B���̂Ƃ�MPS���璊�o�����ő嗱�q����300�܂łƒ�`�B
					int N=0;///FEM�ɓ]�����闱�q��(���Ō��΂�闱�q��)
					int NUM[1000];///�\�ʗ��qID�i�[
					int trans[1000];///NUM[i]�ɑΉ�����NODE�ԍ��i�[�@���q��node�ւ̕ϊ��z��
					cout<<"FEM�p���b�V���쐬�J�n"<<endl;
					int node_number=0;//�����炪�w�肷��ߓ_��
			    
					output_node_data(CON,PART,particle_number,INDEX,MESH,&node_number,NUM,trans,&N,fluid_number,n0);
			    
					FEM(CON,PART,node_number,particle_number,N,NUM,trans,fluid_number,INDEX,MESH,dt,TIME);
					cout<<"�d���͌v�Z�I��"<<endl;
				}
			}
		}
		else if(CON.get_dimention()==3)
		{
			if(CON.get_current_step()==1 || CON.get_current_step()%CON.get_FEM_interval()==0)
			{
				if((t==1 && CON.get_FEM_first_step()==ON) || CON.get_current_step()==1 || CON.get_current_step()%CON.get_FEM_interval()==0)
				{
					//�d���͏�����
					for(int i=0;i<particle_number;i++)	for(int D=0;D<CON.get_dimention();D++)	PART[i].eforce[D]=0;
					int N=0;							//FEM�ɓ]�����闬�̗��q��
					//int *TRANS=new int[fluid_number+1]; //�ߓ_i��TRANS[i]�Ԗڂ̗��q�ɑ���
					vector<int> TRANS(fluid_number+1);	//TetGen�̓����ɂ��vector�ɑΉ�
					
					cout<<"FEM�p���b�V���쐬�J�n----";
					unsigned int timeF=GetTickCount();	//�v�Z�J�n����
					int node_number=0;					//�����炪�w�肷��ߓ_��
					//�Ód�����ɂ��Ă�MPSTOFEM3D_nanoe2�������Ă���(2011.5.2�쐬)
					//if(CON.get_model_number()==14) MPSTOFEM3D_nanoe(CON,PART,particle_number,INDEX,MESH,&node_number,TRANS,&N,fluid_number,n0);
					if(CON.get_model_number()==14) MPSTOFEM3D_nanoe2(CON,PART,particle_number,INDEX,MESH,&node_number,TRANS,&N,fluid_number,n0);
					else if(CON.get_model_number()==15) MPSTOFEM3D_ferrofluid(CON,PART,particle_number,INDEX,MESH,&node_number,TRANS,&N,fluid_number,n0);//����
					else if(CON.get_model_number()==16) MPSTOFEM3D_levitation(CON,PART,particle_number,INDEX,MESH,&node_number,TRANS,&N,fluid_number,n0);//����
//void MPSTOFEM3D_levitation(mpsconfig &CON,vector<mpselastic> &PART,int particle_number,int *INDEX,int **MESH,int *node_number,vector<int> &TRANS,int *N,int fluid_number,double n0)
					
					cout<<"ok  time="<<(GetTickCount()-timeF)*0.001<<"[sec]"<<endl;
					FEM3D(CON,PART,node_number,N,TRANS,fluid_number,particle_number,dt,TIME,t);
					//delete [] TRANS;
				}
			}
		}
	}
	else if(CON.get_mesher()==1)	//TetGen�ɂ�郁�b�V������
	{
		if(CON.get_dimention()==3)
		{
			if(CON.get_current_step()==1 || CON.get_current_step()%CON.get_FEM_interval()==0)
			{
				if((t==1 && CON.get_FEM_first_step()==ON) || CON.get_current_step()==1 || CON.get_current_step()%CON.get_FEM_interval()==0)
				{
					//delaun3D�֌W�̕ϐ��͑S�ĕs�v
					FEM3D_for_TetGen(CON,PART,fluid_number,particle_number);
				}
			}
		}
	}
	//�Ód�����̏ꍇ
	//���q�̈ړ��ɂ��A�����ɐÓd�͂����ꂽ��A�\�ʂ̐Ód�͂�0�������肷��̂�h������
	if(CON.get_FEM_calc_type()==1 && CON.get_model_number()==14)
	{
		double le=CON.get_distancebp();
		double R=CON.get_re()*le;
		for(int i=0;i<fluid_number;i++)
		{
			if(PART[i].type==BOFLUID)//�\�ʂ̏ꍇ
			{
				double F=sqrt(PART[i].eforce[A_X]*PART[i].eforce[A_X]+PART[i].eforce[A_Y]*PART[i].eforce[A_Y]+PART[i].eforce[A_Z]*PART[i].eforce[A_Z]);		
				if(fabs(F)<1.0E-16)
				{
					double newF[DIMENTION];
					for(int D=0;D<DIMENTION;D++)	newF[D]=0;
					int num=0;//���ӗ��q��
					//double W=0;//�d�݂̑��a
					
					for(int k=0;k<PART[i].N;k++)
					{
						int j=PART[i].NEI[k];
						if(PART[j].type==BOFLUID)
						{
							double X=PART[j].r[A_X]-PART[i].r[A_X];
							double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
							double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
							//double dis=sqrt(X*X+Y*Y+Z*Z);
							//double w=kernel(R,dis);
							for(int D=0;D<DIMENTION;D++)	newF[D]+=PART[j].eforce[D];
							//for(int D=0;D<DIMENTION;D++)	newF[D]+=PART[j].eforce[D]*w;
							num+=1;
							//W+=w;
						}
					}
					if(num>0)	for(int D=0;D<DIMENTION;D++)	newF[D]/=num;
					//if(W>0)	for(int D=0;D<DIMENTION;D++)	newF[D]/=W;
					for(int D=0;D<DIMENTION;D++)	PART[i].eforce[D]=newF[D];//�K�p
				}
			}
			if(PART[i].type==FRFLUID)//�����̏ꍇ
			{
				double F=sqrt(PART[i].eforce[A_X]*PART[i].eforce[A_X]+PART[i].eforce[A_Y]*PART[i].eforce[A_Y]+PART[i].eforce[A_Z]*PART[i].eforce[A_Z]);		
				if(fabs(F)>1.0E-16)//�Ód�͂������Ă���ꍇ�͎���̕\�ʗ��q�ɕ��z���āA������0�ɂ���B
				{
					int num=0;//���ӗ��q��
					for(int k=0;k<PART[i].N;k++)
					{
						int j=PART[i].NEI[k];
						if(PART[j].type==BOFLUID)	num++;
					}//�\�ʂ̗��q�������܂���
					for(int k=0;k<PART[i].N;k++)
					{
						int j=PART[i].NEI[k];
						if(PART[j].type==BOFLUID)
						{
							for(int D=0;D<DIMENTION;D++)	PART[j].eforce[D]+=PART[i].eforce[D]/num;	//���ӂ̕\�ʗ��q���Ŋ������Ód�͂�������
						}
					}
					
				}
				for(int D=0;D<DIMENTION;D++)	PART[i].eforce[D]=0;
			}
		}
	}//
	//���q�̓d����[N]�����͂̃f�B���N���l�ɂ���
	if(CON.get_FEM_calc_type()==1)
    {
		double le=CON.get_distancebp();
		//�Ód�͂�MPS�̈��̓f�B���N���l�Ƃ���ꍇ
		if(CON.get_dir_for_P()==2 || CON.get_dir_for_P()==3)
		{
			if(CON.get_eleforce()==1 || CON.get_eleforce()==2)
			{
				//eleforce==3(�\�ʗ�)��eleforce==4(divT)�͂��ꂼ��̊֐����ŕ\�ʐς�������Ƌ��߂�dirP_em���ިظڒl���i�[�ς݂Ȃ̂ŁA�����ł͂����v�Z���Ȃ��B
				
				double S;	//���q1���S������\�ʐ�
				if(CON.get_model_set_way()==0)	S=le*le;
				if(CON.get_model_set_way()==1)	S=sqrt(3.0)/2*le*le;
				double *direct[DIMENTION];	//����g���ĂȂ���ˁH
				for(int D=0;D<DIMENTION;D++) direct[D]=new double [fluid_number];
				for(int i=0;i<fluid_number;i++)
				{
					if(PART[i].type==BOFLUID)  direct_f(CON,PART,i,direct);
					else  for(int D=0;D<DIMENTION;D++) direct[D][i]=0;
				}
				for(int i=0;i<fluid_number;i++)
				{
					double fs=0;//�\�ʗ�
					if(PART[i].type==BOFLUID)//�������̂̏ꍇ��fs=0�Ƃ���
					{
						//�\�ʒ��͓͂����������ƂȂ��Ă���̂ŁA�O�����@���Ōv�Z���Ă����d���͂ɂ��ẮA�ިظڒl�ɂ���Ƃ��̓}�C�i�X��t����K�v������D
						//fs=(PART[i].eforce[A_X]*direct[A_X][i]+PART[i].eforce[A_Y]*direct[A_Y][i]+PART[i].eforce[A_Z]*direct[A_Z][i])/S;//�\�ʗ͓͂������@���x�N�g���Ƃ̓���
						fs=-sqrt((PART[i].eforce[A_X]*PART[i].eforce[A_X]+PART[i].eforce[A_Y]*PART[i].eforce[A_Y]+PART[i].eforce[A_Z]*PART[i].eforce[A_Z]))/S;//�{���͖@���x�N�g���Ƃ̓��ς��Ƃ�Ȃ��Ƃ��߂����ǁA�@���x�N�g���̐��x���悭�Ȃ��̂ŁA���͋��E�����͉����Ȃ��Ă��܂��B�Ȃ̂ł��̂悤�ɂ���B
						PART[i].dirP_em=fs;
					}
				}
				for(int D=0;D<DIMENTION;D++) delete [] direct[D];
			}
		}
    }//////////
	
	//�d���̓v���b�g
	plot_F(CON, PART, fluid_number);
	if(CON.get_F_interval()>0)
	{
		if(CON.get_current_step()==1 || CON.get_current_step()%CON.get_F_interval()==0) plot_F_log(CON,PART,fluid_number);
	}
	//�d����AVS�t�@�C���o��
	if(CON.get_avs_eforce_interval()>0)
	{
		if(CON.get_current_step()==1 || CON.get_current_step()%CON.get_avs_eforce_interval()==0) plot_avs_eforce(CON,PART,fluid_number);
	}
}
*/

void initial_model_input(mpsconfig *CON, int *particle_number, double *TIME)
{
	if(CON->get_restart()==OFF)
	{
		//���f���쐬�ƑS���q���̌v�Z�B
		set_initial_placement_using_MD(CON, particle_number);//set_initial_placement_using_MD(&CON,&particle_number);//���q���͊w�ɂ�郂�f���Z�b�g�B�v�Z���Ԃ͂����邪�A�\�ʂ����ꂢ�ɕ\��
	}
    else if(CON->get_restart()==ON)
    {
		ifstream fin("number.dat");
		if(!fin) cout<<"number.dat cannot be opened" <<endl;
		fin.unsetf(ifstream::dec);
		fin.setf(ifstream::skipws);
    
		fin>>*particle_number;
		TIME=0;//fin>>*TIME;
		fin.close();
		cout<<"�O��ɂ������͂������p��"<<endl;
    }
	
    cout<<"���q����"<<*particle_number<<endl;
}

//TetGen�ɂ�郁�b�V������
void TetGenInterface(mpsconfig &CON, vector<mpselastic> &PART, double **F, int fluid_number, double dt, int t, int particle_number, double n0, double TIME)
{
	
	/////////////////////////////////////////////////
	if(CON.get_EM_interval()==1 || t==1 || (t-1)%CON.get_EM_interval()==0){
			usingTetGen(CON,PART,F,fluid_number,particle_number,t,TIME); //delaun3D�֌W�̕ϐ��͑S�ĕs�v

			if(CON.get_Ferror()==ON){ //�G���[�Ȃ�
				double F_val;
			ifstream fin5("FEM_interval.dat");
			if(!fin5) cout<<"cannot open FEM_interval.dat"<<endl;
			fin5.unsetf(ifstream::dec);
			fin5.setf(ifstream::skipws);
			////�d���͓ǂݎ��
			for(int i=0;i<PART.size();i++)
			{
				if(PART[i].type==MAGELAST){//MAGELAST
				for(int D=0;D<3;D++)
				{
					fin5>>F_val;
					F[D][i]+=F_val;//F[D][i]�͐錾��������Ƀ[���Z�b�g���Ă���
				}
			}
			}
			fin5.close();			
			}
	}else{ //���X�e�b�v���Ɉ�x�����d����v�Z���s���ꍇ
	
		if(CON.get_dir_for_P()!=2 && CON.get_dir_for_P()!=3)//�d���͂����̓f�B���N���l�̂Ƃ��͂����ł͌v�Z���Ȃ�
		{
			double F_val;
			ifstream fin5("FEM_interval.dat");
			if(!fin5) cout<<"cannot open FEM_interval.dat"<<endl;
			fin5.unsetf(ifstream::dec);
			fin5.setf(ifstream::skipws);
			////�d���͓ǂݎ��
			for(int i=0;i<PART.size();i++)
			{
				if(PART[i].type==MAGELAST){//MAGELAST
				for(int D=0;D<3;D++)
				{
					fin5>>F_val;
					F[D][i]+=F_val;//F[D][i]�͐錾��������Ƀ[���Z�b�g���Ă���
				}
			}
			}
			fin5.close();
		}
	}
	//////////////////////////////////////////////////*/
/*	//�d���̓v���b�g
	plot_F(CON, PART, fluid_number,t);
	if(CON.get_F_interval()>0)
	{
		if(t==1 || t%CON.get_F_interval()==0) plot_F_log(CON,PART,fluid_number,t);
	}
	//�d����AVS�t�@�C���o��
	if(CON.get_avs_eforce_interval()>0)
	{
		if(t==1 || t%CON.get_avs_eforce_interval()==0) plot_avs_eforce(CON,PART,fluid_number,t);
	}//*/
}

void file_initialization()
{
	//////////////////////////records//////////////////////////////////////
	//hamiltonian file
	ofstream init1("./Elastic/hamiltonian.dat", ios::trunc);
	ofstream init2("./Elastic/kinetic_energy.dat", ios::trunc);
	ofstream init3("./Elastic/elastic_energy.dat", ios::trunc);
	ofstream init4("./Elastic/potential_energy.dat", ios::trunc);
	
	ofstream init5("aveave_P_history.txt", ios::trunc);
	ofstream init6("node.dat", ios::trunc);
	ofstream init7("time_log.dat", ios::trunc);
	ofstream init8("A-t.dat", ios::trunc);
	ofstream init9("longZ.dat", ios::trunc);
	system("rmdir /s /q Mesh");
	system("rmdir /s /q Lorentz");
	system("rmdir /s /q Speed");
	system("rmdir /s /q FluxContour");
	system("rmdir /s /q FluxAVS");
	
	system("rmdir /s /q Residual");
	system("rmdir /s /q Pressure");
	system("rmdir /s /q Current");
	/////////////////////////make file///////////////////////////////////////
	system("mkdir Mesh");
	system("mkdir Lorentz");
	system("mkdir Speed");
	system("mkdir FluxContour");
	system("mkdir FluxAVS");

	system("mkdir Residual");
	system("mkdir Pressure");
	system("mkdir Current");




	//close file
	init1.close();
	init2.close();
	init3.close();
	init4.close();
	init5.close();
	init6.close();
	init7.close();
	init8.close();
	init9.close();
}

void Make_STL(){
			double normal[368][3];
			double cooda[368][3];
			double coodb[368][3];
			double coodc[368][3];
			double normal_no[368][3];
			double cood1[368][3];
			double cood2[368][3];
			double cood3[368][3];
		/////////////////coordinate file///////////////
		ifstream fin1("cood.txt");
		if(!fin1) cout<<"cood.txt cannot be opened."<<endl;
		fin1.unsetf(ifstream::dec);
		fin1.setf(ifstream::skipws);

		for(int i=0;i<368;i++)
		{       
			for(int D=0;D<3;D++) fin1>>normal[i][D];
			for(int D=0;D<3;D++) fin1>>cooda[i][D];
			for(int D=0;D<3;D++) fin1>>coodb[i][D];
			for(int D=0;D<3;D++) fin1>>coodc[i][D];
			
		}
		fin1.close();	
		//////////////////////////////////////////////////
		/////////////////non coordinate file/////////////
		ifstream fin2("noncood.txt");
		if(!fin2) cout<<"noncood.txt cannot be opened."<<endl;
		fin2.unsetf(ifstream::dec);
		fin2.setf(ifstream::skipws);
		for(int i=0;i<368;i++)
		{       
			for(int D=0;D<3;D++) fin2>>normal_no[i][D];
			for(int D=0;D<3;D++) fin2>>cood1[i][D];
			for(int D=0;D<3;D++) fin2>>cood2[i][D];
			for(int D=0;D<3;D++) fin2>>cood3[i][D];
			
		}
		fin2.close();	
		//////////////////////////////////////////////////

		///////////////�������킹/////////////////////////
		for(int i=0;i<368;i++){
			for(int D=0;D<3;D++) cood1[i][D]+=cooda[i][D];
			for(int D=0;D<3;D++) cood2[i][D]+=coodb[i][D];
			for(int D=0;D<3;D++) cood3[i][D]+=coodc[i][D];
		}
		/////////////////////////////////////////////////

		////////////////STL�t�@�C���o��//////////////////
		ofstream fout1("COIL.STL");
		if(!fout1) cout<<"COIL.STL cannot be opened."<<endl;
		fout1.precision(6);
	
		fout1<<"COIL"<<endl;
		for(int i=0;i<368;i++)
		{       
			fout1<<"   facet normal";
			for(int D=0;D<3;D++) fout1<<" "<<normal[i][D];
			fout1<<endl<<"      outer loop"<<endl;
			fout1<<"         vertex";
			for(int D=0;D<3;D++) fout1<<" "<<cooda[i][D];
			fout1<<endl<<"         vertex";
			for(int D=0;D<3;D++) fout1<<" "<<coodb[i][D];
			fout1<<endl<<"         vertex";
			for(int D=0;D<3;D++) fout1<<" "<<coodc[i][D];
			fout1<<endl<<"      endloop"<<endl;
			fout1<<"   endfacet"<<endl;
			
		}
		fout1.close();	
		/////////////////////////////////////////////////
}