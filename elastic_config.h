//#include "stdafx.h"

//////////////////////////////////////////////
//�w�b�_�t�@�C���ɂ̓\�[�X�̎��Ԃ͏����Ȃ��I���d��`����Ă��܂����߁I�I�I
//�w�b�_�̒��g��extern�錾����̂��肾����ʓI�łȂ�
//�w�b�_�t�@�C���Ƀ\�[�X�̎��̂�������Ă���ꍇ�́A�C���N���[�h�������ׂẴ\�[�X�t�@�C���ɓ����֐����W�J����A���̂��o���܂��B
//////////////////////////////////////////////

//�e���̂̃R���t�B�O��񂾂������B
//mpsparticle, mpselastic�̏������͕ʓr�R���X�g���N�^�ŏ���

class elastic: public mpsconfig
{

	//�e���̌v�Z�p�ϐ��B�����������낢��Ɗi�[���Ă����K�v������ �l�̓R���X�g���N�^�Ōv�Z
	double mag_youngs_modulus;
	double mag_shear_modulus;
	double mag_poisson_ratio;
	double elas_youngs_modulus;
	double elas_poisson_ratio;
	double elas_shear_modulus;
	double mag_lambda;
	double elas_lambda;
	double le;
	double r;
	double mass;
	double inertia;
	double vis;

	//�G�l���M�[�v�Z�p�̕ϐ�
	double hamiltonian;
	double potential;
	double kinetic_energy;
	double elastic_energy;
	double elastic_energy1;
	double elastic_energy2;

	double last_elastic_energy;
	double last_kinetic_energy;
	double last_potential;
	
	//�e��t���O
	bool pivot_check;	//�N�H�[�^�j�I���v�Z�̃K�E�X�����@��pivot�I�����邩�ǂ���
	bool symplectic;
	bool FEM_flag;
	bool nonlinear;

	int symplectic_order;

public:
	elastic(vector<mpselastic> &PART);

	//�|�C���^��Ԃ����set�֐����g���ق�����Ԃ�������Ȃ��idouble **normal�����K�v�����邽�߁j

	double get_mag_youngs_modulus(){return mag_youngs_modulus;}
	double get_elas_youngs_modulus(){return elas_youngs_modulus;}
	double get_mag_shear_modulus(){return mag_shear_modulus;}
	double get_elas_shear_modulus(){return elas_shear_modulus;}
	double get_mag_poisson_ratio(){return mag_poisson_ratio;}
	double get_elas_poisson_ratio(){return elas_poisson_ratio;}
	double get_mag_lambda(){return mag_lambda;}
	double get_elas_lambda(){return elas_lambda;}
	double get_le(){return le;}
	double get_r(){return r;}
	double get_mass(){return mass;}
	double get_inertia(){return inertia;}

	double get_hamiltonian(){return hamiltonian;}
	double get_elastic_energy(){return elastic_energy;}
	double get_elastic_energy1(){return elastic_energy1;}
	double get_elastic_energy2(){return elastic_energy2;}
	double get_kinetic_energy(){return kinetic_energy;}
	double get_potential_energy(){return potential;}

	void set_hamiltonian(double H){hamiltonian=H;}
	void set_elastic_energy(double EE){elastic_energy=EE;}
	void set_elastic_energy1(double EE1){elastic_energy1=EE1;}
	void set_elastic_energy2(double EE2){elastic_energy2=EE2;}
	void set_kinetic(double KE){kinetic_energy=KE;}
	void set_potential(double PE){potential=PE;}

	double get_last_elastic_energy(){return last_elastic_energy;}
	double get_last_kinetic_energy(){return last_kinetic_energy;}
	double get_last_potential(){return last_potential;}
	void set_last_elastic_energy(double EE){last_elastic_energy=EE;}
	void set_last_kinetic_energy(double KE){last_kinetic_energy=KE;}
	void set_last_potential(double PE){last_potential=PE;}

	bool get_pivot_check(){return pivot_check;}
	bool get_symp_flag(){return symplectic;}
	bool get_FEM_switch(){return FEM_flag;}
//	bool get_nonlinear_flag(){return nonlinear;} //nonlinear�̃A�N�Z�T��mpsconfig�Œ�`���Ă���

	void set_symp_flag(bool symp){symplectic=symp;}
	void set_FEM_switch(bool flag){FEM_flag=flag;}
//	void set_nonlinear_flag(bool flag){nonlinear=flag;}
	int get_symp_order(){return symplectic_order;}

	void check_elastic_config();
	void reset_energy();
	
};

void mapcheck(vector<mpselastic> &PART, int i);//map�̃L�[�ƒl���`�F�b�N����
