#pragma once
/*Micro_AVS�t�@�C���o�̓N���X*/
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>

class Micro_AVS{
	std::string Filename;
	int Step;
	int Dimension; 

	std::vector<double> Xr;
	std::vector<double> Yr;
	std::vector<double> Zr;

	std::vector<double> Xdata;
	std::vector<double> Ydata;
	std::vector<double> Zdata;

	static std::vector<std::string> Filelist;
	static bool fast_flag;
public:
	Micro_AVS::Micro_AVS();

	//Micro�QAVS�ŏo�͂���l�̃��X�g�쐬
	virtual void make_list(double xr,double yr,double xdata,double ydata);
	virtual void make_list(double xr,double yr,double zr,double xdata,double ydata,double zdata);
	//Micro_AVS�o��
	virtual void Output_vector_MicroAVS(std::string fname,int step);
	virtual void Output_mgf_MicroAVS(std::string fname,int max_step);
	//���܂łɍ쐬�����t�H���_�����擾
	virtual void Get_folder_name();

protected:
	//�x�N�g���t�@�C��
	//2����
	virtual void _Wright_vector_2d();
	//3����
	virtual void _Wright_vector_3d();
	
	//mgf�t�@�C��
	//2����
//	virtual void _Wright_mgf_2d(); //���̂͂Ȃ�
	//3����
	virtual void _Wright_mgf_3d();

	//�V�����쐬�����t�H���_�����L������
	virtual void _Foldername_storing();




};