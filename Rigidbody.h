#pragma once

/*
Rigidbody�N���X

���̗p�N���X
*/
#include "define.h"

class Rigidbody
{
	double Rg[DIMENSION];	//���̂̏����d�S���W
	double R[DIMENSION];	//���̂̏d�S���W
	double V[DIMENSION];	//���̂̕��ϑ��x
	double theta[DIMENSION];//���̂̉�]�j
	double theta2[DIMENSION];//�������番����
	double I0[3][3];
	double I[3][3];	//���̂̃��[�����g�e���\��	
	double J[3][3];	//I�̋t�s��
//	double F[3];	//���̂ɓ�����
//	double L[3];	//�p�^����
	int rigid_num;			//�����̌ő̔ԍ�

	vector<mpselastic> PARTr;//���̂��\�����闱�q
	vector<double> ri[3];	//���S�Ɨ��qri�̑��΋���
	double Qj[4];		//�N�I�[�^�j�I��
	double Wj[3];		//�p���x
	
	mpsconfig CON;
public:
	Rigidbody();
	virtual ~Rigidbody();
	//�d�S�擾�֐�
	virtual void Get_initial_particle(vector<mpselastic> PARTg);
		virtual void _Get_initial_center();
		virtual void _Get_initial_moment();
		virtual void _Get_new_center();
		virtual void _Get_velocity();
//	virtual void Get_center_update(vector<mpselastic> &PARTup);

	virtual void Renew_part_r_v(vector<mpselastic> &PARTg);//���q�̈ʒu���ړ�
	virtual void Get_rigid_move(vector<mpselastic> &PART);
		virtual void  _Get_I(double (&E)[3][3]);
		virtual void _Rigid_force(vector<mpselastic> &PART);
		virtual void _Calc_quaternion(double q[4],double* ri);
	//�N�I�[�^�j�I����
	virtual void _Quaternion_crossproduct(double A[4],double B[4],double (&Ans)[4]);
	//�e���\����
	virtual void _Tensorial_product(double A[3][3],double B[3][3],double (&Ans)[3][3]);
	//�t�e���\��
	virtual void _Inverse_matrix(double A[3][3],double (&Ans)[3][3]);
};

