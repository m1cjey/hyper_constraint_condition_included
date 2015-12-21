#include "stdafx.h"	
//#include <fstream>

//�R���X�g���N�^
//��͏���
mpsconfig::mpsconfig()
{
/*	ifstream fin("config.txt");
	string buffer;
	stringstream ss;
	if(fin.fail()){
		cout<<"�R���t�B�O�t�@�C�����J���܂���ł���";
		exit(1);
	}
	//boost::format�Ŏ�������̂����N�ŃX�}�[�g�E�E�E
	//get_char()�Ŏ�������#(�R�����g��)��ǂݔ�΂��Ă��ǂ�
	getline(fin, buffer);
	ss<<buffer;
	ss>>step;
	fin.close();*/
	
	step=100000;				//�Sstep��	step=20000;//40000;	//30000;//10000;;	//79*20+1;
	switch_FEM=true;		//FEM�����s���邩���Ȃ��� false
	nonlinear_elastic=false;	//�e���̔���`�v�Z���邩true
	switch_vis=OFF;			//�S�����v�Z���邩���Ȃ����E�E�E����͂��Ƃŏ���
	FEMCG=2;				//FEM�ɂ�����s���@ 0:CG 1:ICCG 2:����ICCG 3:MRTR 4:ICMRTR

//	dt= (switch_FEM==OFF) ? 1.0e-5: 5.0e-6; //0.0001;�s����v���I 0.00001:����(Courant���l����) //Cf. dt_for_FEM=0.000001/2;
	dt=1.0e-5;
	dt_for_FEM=1.0e-3;
	//FEM����0.000001�Ŏ~�܂�E�E�E
	current_step=2;
	current_time=0.0;
	dimension=3;

	interval=1; //10	//particle_movie.mgf�̏o�͊Ԋu�B2�ȏ�̐����ɂ��邱��
	EM_interval=1;//1	//�d����v�Z�����X�e�b�v�Ɉ��s�����B�ʏ��1�ɐݒ�
	motion_interval=1;	//�^��������������Ɉ�������
	
	//���̈��͈ȏ�ɂȂ�����FEM�X�^�[�g cf. PostProcessing.cpp
	//�����l�������ς��ăM���M����_���Ɠ��B���Ȃ��\��������E�E�E
	ave_P_for_FEM_flag=8000000000;//80.0;//75.0;//70.0;

//���f��
	model_number=23;			//4:�������莎���� 7:MRE�A�N�`���G�[�^ 12:����
	model_set_way=1;		//model���Z�b�g������@�@0=�����i�q 1=MD

//���f���P,11��p
/*	R1=9*0.001;
	avoid_step=0;
	avoid_step2=2425;
	avoid_step3=2435;
	avoid_step4=2445;
	avoid_step5=5935;
	avoid_step6=8410;
	avoid_step7=10710;*/

//�o�͓���
	flag_cut_x_movie=OFF;
	flag_cut_y_movie=OFF;

//�d���͌v�Z
	region_shape=1;//1;		//��͗̈�`��@0=������ 1=�~��
	EM_method=3;			//�d����̉�@ 0=OFF 1=FEM 2=BEM 3=���CӰ��Ė@ 4=FEM2
	EM_calc_type=2;			//0=�f���[�j�̂� 1=�d�� 2=�Î��� 3=������ 4=����
//	EM_interval=1;//1		//�d����v�Z�����X�e�b�v�Ɉ��s�����B�ʏ��1�ɐݒ�
	//��͗̈�
	XR=0.2;
	XL=-0.2;
	YU=0.2;
	YD=-0.2;
	/*
	XR=0.1;//0.01;		
	XL=-0.1;//-0.01;
	YU=0.1;//0.01;
	YD=-0.1;//-0.01;*/
	//�~���̈�͂��ꂾ�����߂�
/*	ZU=distancebp*20;
	ZD=-distancebp*20;
	RU=distancebp*10;*/
	
	//FRMcheck�p	15/2/10
	ZU=2.5;//0.10; //0.2
	ZD=-2.5;//0.10; //0.2 				//�t�H -0.01 �R�C��:-0.15 ���:-0.0002
	RU=2.0;//0.10;//0.1;				//��͗̈悪�~���`�ƂȂ�Ƃ��̂��̔��a

//���̂̕����l
	MRE_density=1826;          //water:997.04  �G�^�m�[��:798[kg/m3]
	Silicone_density=980;
	flag_modify_density=OFF;	//���x�⏞���邩�ǂ���
	nensei=1.00; //[Pa�Es]//water:0.001 �G�^�m�[��:0.001084 nensei 8.0
	sigma=0.07196;			//water:0.07196 �G�^�m�[��:0.02361 SUS404:1.85 �\�ʒ��͌W��
	Cp=640/10;     			//water:4.2[kJ/(kg�EK)] �|:800J SUS404:645J/(kgK)�@
	k=28;       			//water:0.6[W/(m�EK)] //�M�`����
	latent_H=0;        		//water:334000J/kg  �|209.3J�@���M�i�G���^���s�[�j
	MP=273;		//�Z�_[K]
	CTE=0.002;//2.1e-4;//2.1e-4;	//���c���W�� SUS410:10.4e-6  ��:2.1e-4

//�e���̂̕����l
	/////MAGELAST/////
	E_m=20000.0;//1000; //200000		//�����O���G500kPa���ƈ��͂��ُ�Ȓl�ɂȂ��Ď~�܂�(dt�̂������Ǝv����)
	v_m=0.35;//0.49	//�|�A�\����
	///////////////////
	/////ELASTIC//////
	E_e=18500.0;
	v_e=0.28;
	///////////////////
//���q�z�u�p
	fluidwidth=20; //30;//40//15[��]	//fluidwidth=20*2;
	distancebp=2.5e-3;///0.001/2;//0.005; //distancebp=0.0125;[mm]
	wlength=2;
	height=0.0;//0.005;    

//��͗̈�
	maxX=0.2;//0.1;	//0.1/2;	//1
	minX=-0.2;//-0.1;	//-0.1/2;
	maxY=0.2;//0.1;	//0.1/2;	//0.4;
	minY=-0.2;//-0.1;	//-0.1/2;	//-0.6; //-1.0
	maxZ=0.2;//0.1;	//0.1/2;	//0.3;
	minZ=-0.2;//-0.1;	//-0.1/2;	//-0.6;  //index�̊֌W��AZ�����ɂ͗]�T�������ƁB

	//FEMcheck�p15/2/10
/*	maxX=0.2;	//0.1/2;	//
	minX=-0.2;	//-0.1/2;
	maxY=0.2;	//0.1/2;	//0.4;
	minY=-0.2;	//-0.1/2;	//-0.6; //-1.0
	maxZ=0.2;	//0.1/2;	//0.3;
	minZ=-0.2;	//-0.1/2;	//-0.6;  //index�̊֌W��AZ�����ɂ͗]�T�������ƁB*/


//���q�@�p�p�����[�^
	re=2.1;//2.1; //���z�E���U�p�̌v�Z�Ɏg��
	re2=2.1; //laplacian�Ɏg��
	re3=3.0; //�\�ʔ���
	re4=3.1; //
	re_elastic=re;//�e���̌v�Z�p
	beta=0.8;
	dx=4;					//index�p�i�q���B�ő嗱�q���a���܂ސ����ɂ��邱��
    times=1;
	
//�\�ʒ��͊֌W                                      
	surface_tension=0;      //0=OFF 1=�`��
	smooth=OFF;				//�X���[�W���O�@1=ON 0=OFF
	SFT=0;					//1=ON 0=OFF �\�ʒ���(surface tension)�̉��x(T)�ˑ����X�C�b�`
	smn=1;					//�X���[�W���O�̉�
	smth_sumP=0;			//curv�̃X���[�W���O�� 0�Ȃ�OFF
    wall_direct=1;			//for ver2~5. BDWALL�̖@���x�N�g���v�Z���@�@0=���� 1=�ʏ� 2=�������� 3=��������
	non_slip=OFF;			//for ver.1 �S������ ON�Ȃ�BDWALL��ɉ��z�I�ɉt�̂����݂���ݒ�B���̍ۂ�BDWALL�̔z�u�ɗ���
	suf_modify=ON;			//�ȗ��v�Z�œ���ꂽ�ȖʂɈ�v����悤�ɗ��q�ʒu���C�����邩�A���Ȃ���
	interpolate_curv=ON;	//�ȗ����v�Z����Ȃ��������q�̋ȗ����A���ӗ��q�̒l�����Ԃ��邩�A���Ȃ���

//���͌v�Z
	iteration_count=0;//1	//���͌v�Z���s�������� �ʏ��1�ɐݒ�	0��FEM	
	solution=1;             //0=CG method 1=ICCG 2=ICCG2�̓������ߖ�
	B_of_P=1;				//���͌v�Z�̉��s�� 0=���ȏ� 1=���x���U 2=0+1 3=���x���U2 4=3+PND 5=(ni-nk)+(nk-n0)
	w_div=100;				//���U�̏d��
	HL_sw=OFF;				//���v���V�A���ɑ΂������̗��U�����s�� 1=ON 0=OFF
	dir_for_P=0;			//���͌v�Z�̍ہA�\�ʂɔ�[����Dirichlet�l�������邩�A���Ȃ��� 0=OFF 1=�\�ʒ��� 2=�d���� 3=����
	initialP=ON;            //���͌v�Z���ɏ����l�Ƃ��Č����͂������邩�ǂ����@1=ON 0=OFF
	CGep=1e-3;				//���͌v�Z���ɂ������������
	omp_P=OFF;				//���͌v�Z����CG�@��openMP���g�p���邩 1=ON 0=OFF
	negativeP=ON;			//�����@1=���� 0=�s����
	set_P=OFF;              //0=OFF 1=�Ð��� 2=�֐�
	Pgrad=3;				//���͌��z�v�Z�@ 1=���ȏ� 2=�\�ʂ̂ݖ@�� 3=2+minP=0�̂Ƃ��\�ʔ��� 4:WLSM 5=CMPS 6=Pj+Pi 7=MPS-AS 
	minP=ON;				//�ŏ����̓X�C�b�` 1=ON 0=OFF
	ave_Pdirect=OFF;		//���͌��z�p�@���x�N�g���̕��ω� 1=ON 0=OFF
	Pgrad_order=2;			//���͌��zver.4�ɂ����鐸�x 1=���` 2=2��
	artP_sw=OFF;			//�l�H���͌v�Z�t���O 0=OFF 1=artP���� 2=Monaghan
	artP=100;				//�l�H���͒l �ʏ��0�ɐݒ�
	Pgrad_times=20;			//���͌��z���v���b�g����ۂ́A�\�ʒ��͂ɑ΂���{�� �ʏ��1�ɐݒ�
	P_AVS=10000;				//microAVS�p�̈��̓t�@�C�����o�͂���step�Ԋu�B0�Ȃ�OFF
	
//���x��֌W
	T_field=0;              //���x���́@�@1=ON 0=OFF
	insulate=1;             //�ǂƂ̒f�M��� 0=�f�M�@1=��f�M
	T_laplacian=0;          //���x�̃��v���V�A���@0=���ȏ��@1=��[i] 2=���U�E���z
	wall_density=7850;		//�ǂ̖��x[kg/m^3]
	wall_Cp=460;			//�ǂ̔�M[J/(kg�EK)]
	wall_k=24;				//�ǂ̔M�`����[W/(m�EK)]
	roomT=313;				//����(20��) int�^�ł悢
	air_cool=ON;			//��C�Ƃ̔M�`�B���l�����邩���Ȃ��� 1=ON 0=OFF
	T_expansion=OFF;		//�M�c���v�Z 1=ON 0=OFF
	buoyant=OFF;			//����(���x��Boussinesq�ߎ�)  1=ON 0=OFF
	TplotZ=0.004;			//3D��͂ɂ����āAXY���ʂ̉��x���o�͂���Ƃ��̂y���W
	T_AVS=20;				//microAVS�p�̉��x�t�@�C�����o�͂���step�Ԋu�B0�Ȃ�OFF

//���b�V��
	mesh_input=2;//2;				//mesh�̍쐬���@ 0:���� 1:Magnet 2:TetGen
	remesh_sw=OFF; //OFF;2012/02/23	//ON:remesh�̈���l�����A��������remesh OFF:FULL�̈�𖈉񕪊�
	boxalpha=2.0;				//�X�[�p�[�{�b�N�X�̑傫�� 1�`10�̎����ɂ��邱�ƁB�ʏ��2�ɐݒ�
	fine=ON;//ON2012/07/08		//�ߓ_�̎����ǉ����s�����ǂ��� 0:OFF 1:�Ӄx�[�X 2:�d�S�x�[�X(���݂͏���)
	co_fine=3.0;				//�ߓ_�����ǉ��ɂ�����W��.�Ӄx�[�X�̏ꍇ�͍ŏ��ӂƍő�ӂ̔䗦�̂������l
	add_points=50000;//50000	//�����ߓ_�ǉ���
	poly_memory=100000;//50000	//poly3D()�Ŋm�ۂ��郁�����̐�
	air_layer=2;				//���̂̎��͂ɋ�C�w�𐶐�����w��(0�Ȃ琶�����Ȃ�)
	layer_depth=0.5;//0.5(2012/03/03);		//��C�w�̕��B�������q�ԋ����̉��{��
	mesh_output_interval=1;//1;	//���b�V�������A�L���v�f�@�̃X�e�b�v�ɑ΂��āA���X�e�b�v�Ɉ�x�o�͂��邩�B
	FEMCGep=5.0e-5;			//1.3e-4PICCG ���X5.0e-6	//�{�v���O�����ł͌��ݕs�g�p15/5/24
	MRTRep=6.8e-4;		//MRTR(&ICMRTR)�@�̎�������	//�{�v���O�����ł͌��ݕs�g�p15/5/24
	FEM_calc_type=2;	//15/5/24	//3;		//0=OFF 1=�d�� 2=���� 3=����(�Q�d��) 4=���� 5=����`�Î���
	ele_type=1;				//(mesher=0�̏ꍇ) �v�f���� 0:�ߓ_�v�f 1:�ӗv�f

//����v�Z
	J_input_way=0;			//�d�����x������@ 0:OFF 1:���� 2:�\�t�g
	J0=0.0;//100;//180000000;		//�����d�����x[A/m2]
	I0=450.0;//200; mA?
	RP=2.0;//2//1.5;//1.28;	//�䓧����
	ele_conduc=1e7;			//�d�C�`����
	Hz=10;			//���g��
	div_Hz=4;				//�P�����̕�����(��͐��x) 4�̔{�������� 40?
	jheat=0;				//�Q�d���ɂ�锭�M���l�����邩�@0=OFF 1=ON
	m_force=1;//1			//�d���͌v�Z���� 1=�ߓ_�͖@ 2=kelvin 3=�ϕ��� 4=divT(�}�N�X�E�F���̉��̓e���\��) 5=VAG 6=�ϕ����ߓ_�͖@ 7=MC
	NLBHI=0;				//�̐ϗ͂ɂ����āA�v�f�a����v�f�g�����߂�ۂɔ���`�����l�����邩�A���Ȃ���(non linier B H inverter)
	NLMH=OFF;				//�l�̎Z�o�ɔ���`�����l�����邩�A���Ȃ���
	magnet_H=1.62577821*distancebp;			//�i�v���΂̍���0.005
	magnet_r=1.62577821*distancebp;//0.01	//�i�v���΂̔��a0.005 �@�@�@�@�@�@�@�@�@�@�@//J_input_way=2:���a�ł͂Ȃ����a�AJ_put_way=0:���a�@�Ǝv����B
	magnet_Z=-magnet_H-distancebp;//-8*distancebp; //-45*0.0005-0.005; //-(fluidwidth)*distancebp-0.01; //-0.0125*0.8		//�i�v���΂̒��S��Z���W 42*0.0005 //-magnet_H/2-(15*distancebp+0.002) ���f��5:-0.035
	magnet_angle=0.0;			//�i�v���΂̒������� 0�Ȃ�+Z�����ƂȂ�B��������p�x���������Ȃ�A���̊p�x[deg]����͂���
	magnet_B=0.5;//0.145;//1.20;			//�i�v���΂̋���[T] Avector3D()�Ŏw��
	magnetic_layers=1;	//�i�v���Ύ��ӂ̋�C�w�̐�1�w�͂��łɂ��� 1+
	uniform_B=0.00; //0.01;	//��l����̑傫��[T]?
	B_times=1;//0.1			//�t�@�C���o�͂���ۂ́A�������x�̔{��
	plot_B_type=1;			//�������x�o�̓^�C�v 1=�x�N�g�� 2=�X�J���[ 0=OFF

	//2012-11-02 0:38 �ӗv�f�𕪊򂳂���I�I�I���ꂪ�Ȃ��ƃf�B���N���l���S�ă[���I�I
	//�R���X�g���N�^�Œl������`�Ŏ��s����Ă��Ȃ��֐���������̂ł́H�H
	uniform_B_sw=OFF;		//��͗̈撆�Ɉ�l����𔭐������邩�ۂ� 0=OFF 1=ON�@

//BEM�֌W
	BEM_elm_type=CONSTANT;	//�v�f�^�C�v 0:���v�f 1:���`�v�f
	BEM_M_solver=2;			//���E�v�f�@�ɂ�����s���@ 0:�K�E�X�̏����@ 1:BiCGStab�@ 2:BiCGStab2�@
	gauss_N=7;				//Gauss�ϕ��̕]���_�� 3,4,7�̂ǂꂩ
	FEM_elm_type=1;			//�v�f�^�C�v 0:�ߓ_�v�f 1:�ӗv�f
	FEM_smn=1;				//�d���̓X���[�W���O�񐔁@0�Ȃ�OFF �}�C�i�X�Ȃ�\�ʂ̂� �v���X�Ȃ�������B
	max_DN=25000;           //Dirichlet�^���E�������Ƃ�ő�ߓ_(��)��
//	FEMCG=1;				//FEM�ɂ�����s���@ 0:CG 1:ICCG 2:����ICCG	//�Ȃ��R�����g�A�E�g����Ă���̂�15/5/24
	CGaccl=1.3;				//ICCG�@�ɂ���������t�@�N�^�@1�̂Ƃ��t�@�N�^OFF
	EMCGep=1.0e-12;			//�d�����ICCG�̎�������1e-5  �A�N�`���G�[�^2.5e-5�@ICCG
	FEMtimes=5;				//�d���͂��v���b�g����ۂ́A�\�ʒ��͂ɑ΂���{�� �ʏ��1�ɐݒ�
	legend_F=2e-7;			//F.dat�̖}��ɏo�͂����[N]
	tree_sw=ON;				//BEM��tree�@���g�����A�g��Ȃ���
	tree_mesh_size=4;		//tree�@�Ŏg�p����ő僌�x���̃Z���̃T�C�Y�Ble�̉��{��
	p_max=3;				//tree�@�̖��������̍ő區���B0����6�̐����ɑΉ�
	plot_F_type=1;			//F.dat�̏o�̓^�C�v  0=[N]�\�� 1=[Pa]�\��
	legend_F_Pa=300;		//F.dat�̖}��ɏo�͂���P��[Pa]

//FEMver2�֌W
	surface_depth=1;		//�\�ʂ̌��݂̔���(�����le�����������̂��{���̌���(�̔���))
	
//�d�E�v�Z
	V=40000;				//�d�E�v�Z�p�d���@�����ЂƂ͂O�Ƃ����Ă���B
	V_step=100;//1500;
	r_perm=80;				//��U�d��
	V_con=0;				//�d�������@0:�p���X�@1:���j�A�@2:���萔
	initial_V=2000;			//�d���������j�A�̂Ƃ��̏����d��
	E_times=1e-10;			//�d�E�o�͔{�� 0�Ȃ�o�͂��Ȃ�
	eleforce=4;				//�Ód�͌v�Z���@ 1:�ߓ_�͖@ 2:�ϕ��� 3:�\�ʗ� 4:divT
	charge=0;				//�d�׍l�� 0=OFF 1=�d�ז��x 2=�N�[������
	plot_E_type=1;			//�d�E�o�̓^�C�v 1=�x�N�g�� 2=�X�J���[ 0=OFF
	
//�e��X�C�b�`�@�ǂ̂悤�Ȍv�Z���l�����邩�A���Ȃ���
	g=-9.8;					//-9.8
	restart=OFF;				//1=ON 0=OFF
	autosave=100;			//�I�[�g�Z�[�u�Ԋu�B�����ɂ������Ƃ��͑傫�Ȑ����������Ă���
	courant=1.0;			//Courant�������@0�Ȃ�n�e�e
	modify_position=OFF;		//���q�ԋ�����le��菬�����ꍇ������C�����邩�A���Ȃ���
	vis_calc_type=1;		//�S�����v�Z��@�@0=POSITIVE=�z��@ 1=NEGATIVE=�A��@
	wall_adheision=2;		//0=�t���[�X���b�v 1=�m���X���b�v  2=2*�m���X���b�v //�����ł́F2
	laplacian=0;//0;        //0=���ȏ��@1=��[i] 2=���U�E���z //�����ł�2
	vis_solver=0;			//�S�������A��͂ŉ����ۂ̍s��\���o�[ 0:CG 1:ICCG
	initial_u=OFF;			//�S�������A�I�ɉ����ۂɁA�����l�Ƃ��Č��ݑ��x����͂��邩�A���Ȃ���
	temporary_r=OFF;		//�z��͌�ɉ��̈ʒu���v�Z���邩���Ȃ����B 1=ON 0=OFF �ʏ��ON�ɐݒ�
	fix_center=0;			//1=ON 0=OFF
	freeon=1;				//���q�ˑ��֌W�֐��@1:���񉻉\ 2:����s�� 4:GPU
	freeon3sw=1;			//freeon3���v�Z���邩���Ȃ��� 1=ON 0=OFF
	surface_judge2=OFF;		//surface_judge2()���g�p���邩���Ȃ���  1=ON 0=OFF
	move_prtcl=OFF;			//�ړ����q���l�����邩���Ȃ��� 1=ON 0=OFF
	move_u_dirct=-3;//-3;//2;			//�ړ����q���ړ���������@���݂́}X����=�}1,�}Y����=�}2,�}Z����=�}3
	move_speed=1.5e-3;//1.5e-3;//1e-3;//8.333*1e-3;//12.5*1e-3;	//�ړ����q�̈ړ����x[m/s]
	check_something=OFF;		//check_something()�����s���邩���Ȃ��� 1=ON 0=OFF

	
//���x�v���b�g�ϐ�
	speed_plot_particle=2;	//���x���v���b�g���闱�q�̎�� 1=���ׂ� 2=fluid 3=��
	speedtimes=5.0e-3;		//���x�v���b�g���́A���W�ɑ΂��鑬�x�̔{��
	speed_face=0;			//3D��͎���speed.dat�̏o�͖� 0=YZ���� 1=XZ
	speed_face_p=0.0;		//3D��͎���speed.dat�̏o�͖ʂ̍��W
	ax_sym_modify=OFF;		//3D����speed.dat�Ɋւ��āA���Ώ̂ɂ��o�͏C�����s�����ۂ��@1=ON 0=OF
	flat_speed_plot=OFF;	//���������̑��x(XY��)���v���b�g���邩���Ȃ���1=ON 0=OFF
	flat_speed_p=0.004;		//flat_speed.dat�̏o�͖ʂ̍��W
	relative_speed=OFF;		//�d�S�ɑ΂��鑊�Α��x���o�͂��邩���Ȃ��� 1=ON 0=OFF
	speed_AVS=OFF;			//microAVS�ɂ��3D���x���z�o�͂��邩���Ȃ��� 1=ON 0=OFF
	legend_speed=0.1;		//speed.dat�̖}��ɏo�͂��鑬�x[m/s]
	set_zero_speed=OFF;		//restart���ɑ��x���[���Z�b�g���邩���Ȃ���  1=ON 0=OFF
	
//���͂ɂ������x�̃v���b�g�ϐ�
	pressure_face=1;			//3D��͎���speed.dat�̏o�͖� 0=YZ���� 1=XZ
	pressure_face_p=0.0;		//3D��͎���speed.dat�̏o�͖ʂ̍��W
	pressure_times=1.0e-6;
	legend_pressure=0.1;

//GPU�֌W
	M_form=ELL;				//CG�@�ɂ�����W���s��̊i�[���@ CSR_scl,CSR_vec,ELL��3���T�|�[�g
	MAX_thread=512;			//�ЂƂ�SM������̍ő�X���b�h���@�ӂ���512

//�t�@�C���o�͕ϐ� 
//	interval=50; //10		//2�ȏ�̐����ɂ��邱��
	AVS=0;                  //0:���ʁ@1:���́@2:���x 3:�ǔ�\�� 4�F�\�ʂ̂� 5:�� 6:���� 
	P_size_AVS=1;			//�o�͂��闱�q�̃T�C�Y(le�̉��{��)
	maxT=1000;
	minT=293;

//AVS�֘A�o�̓t�@�C��
	F_interval=2;			//F.dat�̃��O�����X�e�b�v���ɏo�͂��邩
	avs_eforce_interval=0;			//AVS�d���̓t�@�C�������X�e�b�v���ɏo�͂��邩
	avs_mesh1_interval=200;			//AVS�d�ʃt�@�C��(�f��)�����X�e�b�v���ɏo�͂��邩
	avs_mesh2_interval=200;			//AVS���b�V���t�@�C�L��(�f��)�����X�e�b�v���ɏo�͂��邩
	avs_mesh3_interval=200;			//AVS�\���b�h���f�������X�e�b�v���ɏo�͂��邩

	times_Pa=3E-7;	//gnuplot�p�t�@�C���̃x�N�g�������̔{�� �ʐϗ�[Pa]

	max_pressure=10000;
	min_pressure=0;


//���e���v�Z 
	flag_ELAST=OFF;
	flag_HYPER=ON;
	flag_GRAVITY=ON;
	hyper_density=2970;          //water:997.04  �G�^�m�[��:798[kg/m3]
	c10=30000;//30000;
	c01=20000;//20000;
	flag_wall=OFF;
	h_dis=1.9*distancebp;
	h_vis=1;
	flag_vis=OFF;
	nr_time=1000;	//15/2/8
}


//����: �ŏ��ɕ�����̂͗��q�������̑̐ς����B���Ƃ͏����z�u�ɋ����Ė��x���␳�����̂ł����������
//�[�U��: ���� 68%�@�Ŗ� 74%
double mpsconfig::get_particle_mass()
{
	double mass;
	if(model_set_way==0)//�����i�q�̂Ƃ��[�U����PI/6=52.3%
	{
		if(dimension==2){
			mass=MRE_density*distancebp*distancebp;
		}else{
			mass=(0.523598775)*MRE_density*distancebp*distancebp*distancebp;
		}
	}
	else if(model_set_way==1)//�Ŗ��i�q�̂Ƃ�(�[�U����sqrt(2)*PI/6=0.74048049)
	{
		if(dimension==2){
			mass=MRE_density*sqrt(3.0)/4*distancebp*distancebp;
		}else{
			//mass=density*sqrt(2.0)/12*distancebp*distancebp*distancebp;
			//mass=density*distancebp*distancebp*distancebp/sqrt(2.0);
			//(4*PI/3)/2^3=4.1887902/8=0.523598776
			//(sqrt(2)*PI/6)*(4*PI/3)/2^3=0.387714678
		//	mass=(0.74048049)*distancebp*distancebp*distancebp*density;
			mass=distancebp*distancebp*distancebp*MRE_density;
		}
	}
	else{
		cout<<"���f���̐ςݕ����s��ł� ���ʂ��v�Z�ł��܂���"<<endl;
		exit(1);
	}

	return mass;
}