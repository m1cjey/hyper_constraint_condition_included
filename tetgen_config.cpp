#include "stdafx.h"

tetgen_config::tetgen_config()
{
	//�Ód�������u�̐��@
	radius_column=0.0005;		//�~���d�ɔ��a
	length_column=0.02;			//�~���d�ɒ���
	height_plate=0.0075;		//���d�ɍ���
	length_plate=0.02;			//���d�Ɉ�Ӓ���
	thickness_plate=0.0005;		//���d�Ɍ���
	length_base=0.02;			//�y���Ӓ���
	thickness_base=0.005;		//�y�����
	
	//���b�V���̑e��(���E��1�v�f�̕ӂ̒������ǂꂮ�炢�ɂ��邩)
	fine_air=0.01;				//��C�̈拫�E�̑e��		0.01
	fine_plate_t=0.00025;		//���d�Ɍ��ݕ����̑e��	0.00025
	fine_plate_L=0.00050;		//���d�ɕ��ʕ����̑e��	0.0005
	fine_column_L=0.00015;		//�~���d�ɂ̒��������̑e��	0.00015
	fine_base=0.0010;			//�y��\�ʂ̑e��			0.0010

	//���H�\�ʂ̃��b�V���w�̐ݒ�
	num_layer_out=1;	//���̊O�����b�V���w��
	num_layer_in=0;		//���̓������b�V���w��
	thick_layer=0.3;	//���E���b�V��1�w�̌���(le�̉��{��)

	//�������̗v�f���폜����臒l 
	del_length=3.0;		//le�̉��{�ȏ�̕ӂ����v�f��������2.0
}

tetgen_config::tetgen_config(mpsconfig &CON)
{
	//���ΐ��@ cf. MPSTOFEM_MRE()<-MPS_TO_FEM3D.cpp
	magnet_height=CON.get_magnet_H();
	magnet_radius=CON.get_magnet_r();

	//���b�V���̑e��(���E��1�v�f�̕ӂ̒������ǂꂮ�炢�ɂ��邩)
	fine_air=0.01;				//��C�̈拫�E�̑e��		0.01
	fine_plate_t=0.00025;		//���d�Ɍ��ݕ����̑e��	0.00025
	fine_plate_L=0.00050;		//���d�ɕ��ʕ����̑e��	0.0005
	fine_column_L=0.00015;		//�~���d�ɂ̒��������̑e��	0.00015
	fine_base=0.0010;			//�y��\�ʂ̑e��			0.0010

	//���H�\�ʂ̃��b�V���w�̐ݒ�
	num_layer_out=1;	//���̊O�����b�V���w��
	num_layer_in=0;		//���̓������b�V���w��
	thick_layer=0.3;	//���E���b�V��1�w�̌���(le�̉��{��)

	//�������̗v�f���폜����臒l 
	del_length=2.0;		//le�̉��{�ȏ�̕ӂ����v�f��������

}
