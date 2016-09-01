#include "stdafx.h"	//��v�ȃw�b�_�[�t�@�C���͂܂Ƃ߂Ă��̂Ȃ��B
//#include"define.h"	//#define �i�[


/////BEM2D�p�̃N���X

class point2D//�ߓ_�N���X
{
public:
	double r[2];
	int boundary_condition; //���E���� 0=Diric 1=Neumn
	double potential;		//�|�e���V����
	double slop1;			//�@�������̌��z
	double slop2;			//�@�������̌��z
	double C;				//���p
	double L;				//���v�f�̏ꍇ�̂ݎg�p�B���v�f�̏ꍇ�ANODE[n].L�Őߓ_n�𒆓_�Ƃ�����v�f�̗v�f�̒�����������
	int particle;			//�Ή����闱�q�ԍ��i�[ �Ή�������̂��Ȃ��ꍇ�́|1���_�~�[�Ƃ��Ċi�[
};

class element2D//�v�f�N���X
{
public:
    int node[2];//�v�f���\�z����node�̔ԍ� 
	double r[2];		//���_�̍��W
	double L;//����
	double direct[2];		//�@���x�N�g��
	int boundary_condition; //���E���� 0=Diric 1=Neumn ���v�f�̎��Ɏg�p
	int map;//�}�b�s���O
	int material;//�ގ�
};

class REGION//�̈�N���X
{
public:
	int start;	//start����v�Z�_�ԍ�
	int end;
};

////////////


/////BEM3D�p�̃N���X/////////
class BEMpoint3D//�ߓ_�N���X
{
public:
	double r[3];
	int boundary_condition; //���E���� 0=Diric 1=Neumn
	double potential;		//�|�e���V����
	double slop1;			//�@�������̌��z
	double slop2;			//�@�������̌��z
	double C;				//���p
	int particle;			//�Ή����闱�q�ԍ��i�[ �Ή�������̂��Ȃ��ꍇ�́|1���_�~�[�Ƃ��Ċi�[
};

class BEMelement3D//�v�f�N���X
{
public:
    int node[3];//�v�f���\�z����node�̔ԍ� 
	double r[3];		//���_�̍��W
	double S;		//�ʐ�
	double direct[3];		//�@���x�N�g��
	int boundary_condition; //���E���� 0=Diric 1=Neumn ���v�f�̎��Ɏg�p
	int map;//�}�b�s���O
	int material;//�ގ�
};
