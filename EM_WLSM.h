#include "stdafx.h"	//��v�ȃw�b�_�[�t�@�C���͂܂Ƃ߂Ă��̂Ȃ��B
#include"define.h"	//#define �i�[

class WLSM_point3D//�ߓ_�N���X
{
public:
	double r[DIMENSION];
	int material; //�ގ�
	int boundary_condition; //���E���� 0=���ʁ@1,2=�Œ苫�E -1 ���R���E����(���z���[��)
	int particleID;			//�Ή����闱�q�ԍ� ���݂��Ȃ��Ƃ���-1���i�[
	int remesh;				//�����b�V���̈�ɑ����邩�A���Ȃ���
	int surface;			//���̂̕\�ʂł��邩�A�Ȃ����B1=ON 0=OFF, ���E�����Ƃ͕�
	int depth;				//�E�ʂ���̐[���@�E�ʂ�0�Ƃ���
	int calced;				//calced=ON�Ȃ�A���v���V�A���̂��߂̎��ӗ��q�̊�^���͌v�Z�ς݂Ȃ̂ŁA���̂܂܎g�p����΂悢�BOFF�Ȃ�Čv�Z���ċ��߂�
	double rp;				//��U�d���������͔䓧����

	int index;//�i�[����Ă���i�q�̔ԍ�
	int N2;   //re2���ɑ��݂�����ӗ��q��
	vector<int> NEI2;
	double L;				//���ϗ��q�ԋ���
	double R;				//�e�����a
	vector<double> fai;		//���v���V�A���𗣎U������̂ɕK�v�Ȏ��ӗ��q�̏d��

	double E[3];			//�d�E�܂��͎��E�x�N�g��
	double potential;		
	double F[3];			//�d����
	double Fs;				//�d������[N/m^2]
	double normal[3];		//�O�����@���x�N�g��

	int root;				//��C�v�Z�_�Ɋւ��āA���v�Z�_�̖@����萶�����ꂽ�Ȃ�root=FLUID,�Ƃ������ɁA�������̍ގ����i�[����B��͗̈�Ȃ��C�B�@��C�v�Z�_�ȊO��root�͎Q�Ƃ��Ȃ�(-1���i�[)
	
	int flag;				//�K���Ȏg���񂵗p�ϐ�
};

class INDEX//�i�q
{
public:
	vector<int> node;		//�i�[����ߓ_�ԍ�
};