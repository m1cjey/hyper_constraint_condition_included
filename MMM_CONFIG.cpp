#include "stdafx.h"	//��v�ȃw�b�_�[�t�@�C���͂܂Ƃ߂Ă��̂Ȃ��B
#include"define.h"	//#define �i�[
#include"MMM_CONFIG.h"	//class CON��`

MMM_config::MMM_config()
{
	///////��͏���
	Hf_type=1;		//�O������@0:��l���� 1:����
	Hf_H = 0.03 / (4 * PI * 1e-07);		//�O������H�̋���
	isWLSM = 0;
	eForce = 0;             //�d����(0:���i��,1:Kelvin��)
	force_t = 1;            //�d���͌v�Z��@(0:���͂̂�,1:�̐ϗ͂̂�)

	P_interval = 20;
}