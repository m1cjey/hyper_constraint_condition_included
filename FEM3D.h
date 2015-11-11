//#include "stdafx.h"
#pragma once
class point3D//�ߓ_�N���X
{
public:
	double r[DIMENSION];
	int material; //�ގ�
	int boundary_condition; //���E���� 0=���ʁ@1,2=�Œ苫�E
	int particleID;			//�Ή����闱�q�ԍ� ���݂��Ȃ��Ƃ���-1���i�[
	int remesh;				//�����b�V���̈�ɑ����邩�A���Ȃ���
};

class element3D//�v�f�N���X
{
public:
    int node[5];//�v�f���\�z����node�̔ԍ�  �����v�܂�聨�Ă��؂�̏��Ɋi�[
	int elm[5];//�v�f�Ɨאڂ���v�f�ԍ� //�\�ʂ�0�Ɗi�[
	int sides[6+1];//�v�f���\������Ӕԍ�
	double volume;//�̐ς̂U�{
	double r[DIMENSION];//�{���m�C�_���W(�O�ڋ����S)
	double RR;//�O�ڋ����a�̓��
	int map;//�}�b�s���O
	int material;//�ގ�
};

class sides//�ӃN���X
{
public:
	int node[2+1];//�ӂ��\������ߓ_�ԍ��i�[
	int boundary_condition; //���E���� 0=���ʁ@1,2=�Œ苫�E
};

void VOLT3D(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double *V,int *jnb,double TIME,vector<mpselastic> &PART,int fluid_number,int **nei,double *RP);
void ELECTRO3D(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double *V,double **Ee);
void potential_calculation(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,int *jnb,double TIME,vector<mpselastic> &PART,int fluid_number,int **nei,double **F);
void conductive_approximation(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **Ee,int *jnb,int **nei,double *RP,vector<mpselastic> &PART,double **F,int fluid_number);
void import_J0_density(mpsconfig &CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **current);
void check_J0(mpsconfig &CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **current);
void Avector3D_node_eddy2(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,double **A,int *jnb,int *depth,int **nei2,int *branch_num,double **old_A,double II,double dt,int mps_num,double *V,int t);
void calc_transitional_EM_field(mpsconfig &CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,int *jnb,double dt,double TIME,int t,int **nei,int KTE,vector<mpselastic> &PART,int fluid_number,double **F,int KTJ);
void Bflux3D_node(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **A,double **B);
int locate3D(vector<point3D> &NODE,vector<element3D> &ELEM,int nelm,double xp,double yp,double zp);

void plot_magnetic_flux_density(mpsconfig &CON, vector<mpselastic> &PART, vector<point3D> &NODE, vector<element3D> &ELEM, int nelm, double **B, int t); //�������x���v���b�g

void poly3D2(vector<point3D> &NODE,vector<element3D> &ELEM,int *iv,int *kv,int ip,int *nelm,mpsconfig *CON);
void smoothingFs3D(mpsconfig &CON,vector<mpselastic> &PART,int fluid_number,double *Fs);
void smoothingF3D(mpsconfig &CON,vector<mpselastic> &PART,int fluid_number,double *F[3],int t);

//�ӗv�f�����֐�
int make_edge_element(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,int *jnb,int **nei,sides *SIDE,int *branch_num,int **nei2,int KTE);
void calc_static_magnetic_field(mpsconfig &CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,int *jnb,double dt,double TIME,int t,int **nei,int KTE,vector<mpselastic> &PART,int fluid_number,double **F,int KTJ);
void calc_variable_magnetic_field(mpsconfig &CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,int *jnb,double dt,double TIME,int t,int **nei,int KTE,vector<mpselastic> &PART,int fluid_number,double **F);
void set_boundary_condition3D_edge(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,int *dn,int *NN,double *PHAT,double *A);
void Avector3D(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,double *A,int *jnb,int *branch_num,double **current,double *RP);
void Avector3D2(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,double *A,int *jnb,int *branch_num,double **current,double *RP,double II,int *depth,int t);
void Bflux3D_side(mpsconfig &CON, vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,double *A,double **B,double *RP);

//��
void direct_divT3D(mpsconfig &CON,vector<mpselastic> &PART,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **Be,int fluid_number,double **F,double *RP,int *jnb,int **nei);
void kelvin_force3D(mpsconfig &CON,vector<mpselastic> &PART,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **Be,int fluid_number,double **F,double *RP,int*jnb,int **nei);
void coroid_pole(mpsconfig &CON,vector<mpselastic> &PART,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **Be,int fluid_number,double **F,double *RP,int*jnb,int **nei);
void NODE_F3D(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **Ee,int *jnb,int **nei,double *RP,vector<mpselastic> &PART,double **F,int fluid_number);

//�d�����x�v�Z�֐�(�ӗv�f�p)
void calc_current(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,int *jnb,int *branch_num,double **current,int *depth,double II);
//�ӗv�f�d�����x�v�Z�֐�
void denryu_side(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,double *T,double **current,double J1);
