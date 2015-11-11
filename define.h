#pragma once

#define DIMENSION 3
#define A_X 0
#define A_Y 1
#define A_Z 2
#define ON 1
#define OFF 0
#define PI 3.14159265358979323846264338327950288
#define POSITIVE 0
#define NEGATIVE 1
#define ERR 1.0e-12
#define EPS 1.0e-12

//////////////���q�f��///////////////
//���[���O���̐����Œu��
//���ꗱ�q 000_099
#define GHOST 000 //�v�Z���Ȃ����q
#define FACE_P 001//�\��

//�C�̗��q 100_199
#define GAS 100
#define AIR 101   //��C

//���̗��q 200_299
#define FLUID 200
#define WATER 201

//�ő̗��q 300_399
#define SOLID 300 //�ő�
#define WALL 301 //�Œ��
#define ELASTIC 302//�e����
#define MAGELAST 303//�����G���X�g�}�[
#define MAGNET 304
#define IRON 305
#define COIL 306
#define ELECTRODE 307	//�d��
#define PLATE 308
#define ELECTRODE1 309	//�d��1(nanoe�~���d��)
#define ELECTRODE2 310	//�d��2(nanoe���d��)
#define TERMINAL1 311	//�������莎���p�̒[�q
#define TERMINAL2 312	//�������莎���p�̒[�q
#define MAGELAST2 313 //�����G���X�g�}�[�Q
#define HYPERELAST 314

//GPU CG�@�ɂ�����W���s��̊i�[���@
#define CSR_scl 0
#define CSR_vec 1
#define ELL 2

#define Diric 0
#define Neumn 1
#define BOTH  2
#define CONSTANT 0
#define LINER 1

#define SIGNIFY 18		//�t�@�C���o�̗͂L������

#define UNDEFINED -1	//����`