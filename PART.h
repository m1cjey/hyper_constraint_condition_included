#ifndef MPSPARTICLE //ifndef�ŃR���p�C�����Ȃ�
#define MPSPARTICLE

class vector3D
{
	double component[3];
public:
	vector3D(double *r0);
	void operator=(double vc[3]);
	double *get_comp(){return component;}
	double get_comp(int D){return component[D];}
};

class mpsparticle
{
public:
	unsigned ID;
	double r[DIMENSION];
	double initial_r[DIMENSION];
	double q0[DIMENSION];

	double u[DIMENSION];

	double r_temp[DIMENSION];
	double u_temp[DIMENSION];

	double n[DIMENSION];	//�P�ʖ@���x�N�g��(TetGen�p)�@�ǂ��Ō���H

	double ang[DIMENSION+1];//�p�x�E�N�H�[�^�j�I��
	double ang_u[DIMENSION];
	double ang_temp[DIMENSION+1];//�p�x�E�N�H�[�^�j�I��
	double ang_u_temp[DIMENSION];
	double P;				//����
	double volumetric_strain;//�̐ςЂ���

	double p[DIMENSION];	//������

	bool contact;

	double dirP_st;	//�\�ʒ��͂ɂ�鈳�̓f�B���N���l
	double dirP_em;	//�d���͂ɂ�鈳�̓f�B���N���l

	double h;				//�G���^���s�[
	double PND;				//re��p�������q�����x
	double PND2;
	double PND0;		
	double val;					//���̓s�x�K���Ȓl�̊i�[�ɗ��p
	double eforce[DIMENSION];	//�d����

	int fly;				//��U�^�C�v(TOUCH, GROUP, ISOLATION)
	int type;				//FLUID INWALL OUTWALL
	int materialID;			//�ގ��ԍ� 1�Ȃ�density,2�Ȃ�density2�Ȃǂ��g�p
	int surface;			//0:���� 1:�\��
	int index;				//�i�[����Ă���i�q�̔ԍ�
	int N;					//re���ɑ��݂�����ӗ��q���E�E�E���ꂪpublic�ɂȂ��Ă���̂œY���ɂ��������Ă��܂�����I�I�I
	int N2;					//re2���ɑ��݂�����ӗ��q��
	int N3;					//re3���ɑ��݂�����ӗ��q��
	int NEI[300];			//re���ɑ��݂�����ӗ��q�ԍ�
	int NEI2[300];
	int NEI3[450];
	double PAcc[3];			//[murao] Permanent acceleration

	unsigned int toFEM;		//FEM�ɓ]�����邩�ǂ����BON or OFF
	double dir_Pst;			//�\�ʒ��͂ɂ�鈳��
	double dir_Pem;			//�d���͂ɂ�鈳��

	double ave_lambda;		//���q�̕��σ�
	//PND=particle_number_density
	void set_PND(const double pnd){PND=pnd;}
	double get_PND() const{return PND;}
	double F[DIMENSION];
};

//�e���̌v�Z�p�h���N���X
class mpselastic: public mpsparticle
{
	double youngs_modulus;	//�����O��
	double poisson_ratio;	//�|�A�\����
	

	//�v�͎��ӂɂ��闱�q��ID���������Ă����Ηǂ��̂����A���̓s�x�ʒu�֌W���v�Z���Ă�����ʓ|
	vector<int> initial_neighboursID;	//�����̎��ӗ��qID
	vector<int> current_neighboursID;	//���݂̎��ӗ��qID

	//�ʒu���W
	vector<double> initial_neighbourX;	//freeon�ŏ���������/ size==initial_neighboursN
	vector<double> initial_neighbourY;
	vector<double> initial_neighbourZ;
	vector<double> current_neighbourX;	//freeon�ŏ���������
	vector<double> current_neighbourY;
	vector<double> current_neighbourZ;

	vector<double> initial_distancebps;
	vector<double> current_distancebps;

	vector<vector3D> r0_ij;
	vector<vector3D> r0_ji;

	//�͊w�p�����[�^(i�ɑ�������̂�elastic�N���X���炱����ֈړ�����)
	double elastic_density;					//�e���̖̂��x�iPND���狁�߂�)

	double stress_accel[DIMENSION];			//���͉����x
	double stress_visco_accel[DIMENSION];	//����f���͉����x
	double pressure_accel[DIMENSION];		//���͉����x
	double total_accel[DIMENSION];			//�����x���v
	
	bool stop_on_floor; //���Ƃ̐ڐG����
	bool acceleration_upward;//�����x�̔���i����������������j

public:

	mpselastic();

	//�C�e���[�^��߂�l�Ɏg�����Ƃ͉\�H

	vector<int>& get_initial_neighboursID() {return initial_neighboursID;}	//�Q�Ƃ�Ԃ��ق��������͂��i����initial_neighboursID�ɃA�N�Z�X���Ă���̂Ɠ����ɂȂ�̂Łj
	vector<int>& get_current_neighboursID() {return current_neighboursID;}	//�Q�Ƃ�Ԃ����������͂��i����initial_neighboursID�ɃA�N�Z�X���Ă���̂Ɠ����ɂȂ�̂Łj

	void set_current_neighboursID(const int j){current_neighboursID.push_back(j);}//current_neighboursID.size()==N
	void reset_current_neighboursID(){current_neighboursID.clear();}

	void get_initial_neighbours_position(int k0, double &X0, double &Y0, double &Z0){X0=initial_neighbourX[k0]; Y0=initial_neighbourY[k0]; Z0=initial_neighbourZ[k0];}
	void get_initial_neighbours_position(int k0, double *r0){r0[0]=initial_neighbourX[k0]; r0[1]=initial_neighbourY[k0]; r0[2]=initial_neighbourZ[k0];}
	void get_current_neighbours_position(int k, double &X, double &Y, double &Z){X=current_neighbourX[k]; Y=current_neighbourY[k]; Z=current_neighbourZ[k];}
	void get_current_neighbours_position(int k, double *r){r[0]=current_neighbourX[k]; r[1]=current_neighbourY[k]; r[2]=current_neighbourZ[k];}
	void set_current_neighbours_position(double X,double Y, double Z){current_neighbourX.push_back(X); current_neighbourY.push_back(Y); current_neighbourZ.push_back(Z);}
	void reset_current_neighbours_position(){current_neighbourX.clear(); current_neighbourY.clear(); current_neighbourZ.clear();}

	void set_current_distancebps(const double dis){current_distancebps.push_back(dis);}//���ꂢ��Ȃ��Binitial_distancebps�����ŏ\��
	double get_initial_distancebps(const int k){return initial_distancebps[k];}
	double get_current_distancebps(const int k){return current_distancebps[k];}
	vector<double>& get_initial_distancebps(){return initial_distancebps;}
	vector<double>& get_current_distancebps(){return current_distancebps;}

	void set_r0_ij(double *r0);
	void set_r0_ji(double *r0);
	vector<vector3D>& get_r0_ij(){return r0_ij;}
	vector<vector3D>& get_r0_ji(){return r0_ji;}

	double get_density() const{return elastic_density;}
	void set_density(const double d){elastic_density=d;}
	void reset_current_distancebps(){current_distancebps.clear();}//claer()���Ȃ��Ɣz�񂪑���������

	double get_pressure_accel(int D){return pressure_accel[D];}
	void set_pressure(const double *pressure){for(int D=0;D<3;D++) pressure_accel[D]=pressure[D];}
	
	double get_stress_accel(int D){return stress_accel[D];}
	double get_stress_visco_accel(int D){return stress_visco_accel[D];}
/*	void reset_stress(){for(int D=0;D<3;D++){ 
		stress_accel[D]=0;
		stress_visco_accel[D]=0;
	}
	}//*/

	void add_pressure(const double* press){for(int D=0;D<3;D++) pressure_accel[D]+=press[D];}

	void add_stress_accel(const double* st){for(int D=0;D<3;D++) stress_accel[D]+=st[D];}
	void add_stress_visco_accel(const double* sv){for(int D=0;D<3;D++) stress_visco_accel[D]+=sv[D];}

	void mul_pressure(const double coef){for(int D=0;D<3;D++) pressure_accel[D]*=coef;}

	void mul_stress_accel(const double coef){for(int D=0;D<3;D++) stress_accel[D]*=coef;}
	void mul_stress_visco_accel(const double coef){for(int D=0;D<3;D++) stress_visco_accel[D]*=coef;}

	void initialize_particles(elastic &ELAST, const int t);//�����elastic�̃����o�ŗǂ��̂ł́H
	void reset_particle_acceleration();

	double get_total_accel(int D){return total_accel[D];}
	void set_total_accel(const double *acceleration){for(int D=0;D<3;D++) total_accel[D]=acceleration[D];}


	//�e���萔�̃A�N�Z�T
	double get_youngs_modulus(){return youngs_modulus;}
	double get_poisson_ratio(){return poisson_ratio;}
	void set_youngs_modulus(double E){youngs_modulus=E;}
	void set_poisson_ratio(double nu){poisson_ratio=nu;}

	//������
	bool get_stop_on_floor(){return stop_on_floor;}
	void set_stop_on_floor(bool flag){stop_on_floor=flag;}

	//�����x�̔���
	bool get_acceleration_upward(){return acceleration_upward;}
	void set_acceleration_upward(bool flag){acceleration_upward=flag;}

	//�����x�`�F�b�N
	void check_acceleration();

	//�e���W���`�F�b�N
	void check_elastic_constant();

	//map�ɓ����ƃ\�[�g����Ă��܂��̂Œ���
	//�ڐG����Ŏg��
	//freeon��neighbour_distance�����߂�initialize_particles��map�ɃR�s�[����
	map<int, double> distance0;		//�����̎��ӗ��q�Ƃ̋���(size��neighbours0), ID���L�[�ɂ��ė��q�ԋ������Q�Ƃ��� //NEI�̒u��������}��
	map<int, double> distance;		//���݂̎��ӗ��q�Ƃ̋���(size��neighbours), ID���L�[�ɂ��ė��q�ԋ������Q�Ƃ���
	
};

#endif