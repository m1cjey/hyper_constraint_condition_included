///////////////////////����///////////////////////
//TetView�Ńt�@�C�����J���Ă���ƃv���Z�X���I������
//output�t�@�C���̓ǂݏ������o���Ȃ��̂Ńv���O�������i�܂Ȃ��I
//�K��TetView�̃v���Z�X���I�����Ă���.exe�����s���邱�ƁI�I

//TetGen�Ăяo���֐�  TetGen�̓����

#include "stdafx.h"


//TetGen����  ���f���ɂ�蕪��
void tetgen_function::call_TetGen(mpsconfig &CON, vector<mpselastic> &PART, int fluid_number, int particle_number, vector<tetgen_node> &NODEall, vector<tetgen_facet> &FACEall, vector<tetgen_element> &ELEMall, vector<int> &TRANS)
{
	
	cout<<"TetGen�ɂ�郁�b�V�������J�n"<<endl;
	clock_t t1=clock();	//�N���b�N���擾


	//////���f�����ɕ���/////////////////////////////
	TetGen_elastic(CON,PART,fluid_number,particle_number,NODEall,FACEall,ELEMall,TRANS); //�����G���X�g�}�[�p
	/////////////////////////////////////////////////


	clock_t t2=clock();	//�N���b�N���擾�@
	cout<<"���b�V����������  CPU time="<<(t2-t1)/CLOCKS_PER_SEC<<endl;
}


//////�����艺�ɁA��肽�����f���̃��b�V�������v���O�������L�q���Ă�������/////////////////////////////////////////////

//�e���̗p�@���b�V������
void tetgen_function::TetGen_elastic(mpsconfig &CON, vector<mpselastic> &PART, int fluid_number, int particle_number, vector<tetgen_node> &NODEall, vector<tetgen_facet> &FACEall, vector<tetgen_element> &ELEMall, vector<int> &TRANS)
{
	//tetgenio�N���X�錾�Ə�����
	tetgenio in, out;
	in.initialize();
	out.initialize();

	//TetGen�ݒ�N���X
	tetgen_config TET;

	//�S�̔z�񏉊���
	NODEall.clear();
	FACEall.clear();
	//�e���i�p�z��i�g���񂵁j
	vector<tetgen_node> NODE;
	vector<tetgen_facet> FACE;


	////////////�G���X�g�}�[�̈�////////////
	NODE.clear();
	FACE.clear();
	SetElastBoundary(CON, PART, TET, NODE, FACE);
	SetTRANS(NODE, TRANS);	//���q���ߓ_�ɑΉ����鎞�Ɏ��s�H�H�H
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, MAGELAST);//*


	//////////////�d���Η̈�/////////////
	NODE.clear();
	FACE.clear();
	SetMagnetBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, MAGNET);//*/
	/////////////�R�C���̈�//////////////
	/*NODE.clear();
	FACE.clear();
	SetCOILBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, COIL);*/
	////////////�S�S�̈�///////////////////*
/*NODE.clear();
	FACE.clear();
	SetIRONBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, IRON);*/

	////////////��C�̈�////////////
	NODE.clear();
	FACE.clear();
	SetAirBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);

/*	//////////////��C�̈�ǉ��ߓ_////////////
	NODE.clear();
	FACE.clear();
	SetAirFineBoundary(CON, PART, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);
	//*/
	

	//.node�t�@�C���쐬
	MakeNodeFile(CON, NODEall, "all_boundary.node");
	//.poly�t�@�C���̏o��
	MakePolyFile(CON, TET, NODEall, FACEall, "all_boundary.poly");

	cout<<"poly�t�@�C���ǂݍ���"<<endl;
	in.load_poly("all_boundary");	//.poly��ǂݍ��񂾂玩���I��.node���ǂݍ��ނ炵��
	cout<<"���������J�n"<<endl;
	tetrahedralize("pq1.1a1.0e-4AYYn", &in, &out);	//1.1�����ł͐؂�Ȃ� �f�t�H���g��rqa1.1AYYn  pq1.3a1.67e-7AYYn
	out.save_nodes("output");
	out.save_elements("output");
	//out.save_faces("output");
	//out.save_neighbors("output");
	//*/
	cout<<"�ގ��̏C���J�n\n";
	//cout<<"�ގ��C��"<<endl;
	//ModifyAttribute(CON, TET, NODEall, ELEMall, in, out);
	ModifyAttribute_tetgenio(CON, TET, NODEall, ELEMall, in, out);	//tetgenio�𒼐ڕҏW

	//�o��
	cout<<"TetView�p�t�@�C���o��"<<endl;
	//MakeNodeFile(CON, NODEall, "output_final.node");
	//MakeElemFile(CON, ELEMall, "output_final.ele");
	out.save_nodes("output_final");
	out.save_elements("output_final");

	//�ߓ_�E�v�f�f�[�^�擾
	cout<<"TetGen���ߓ_�E�v�f�f�[�^�擾"<<endl;
	GetPointList(NODEall, in, out);
	GetTetrahedronList_full(ELEMall, in, out);
}