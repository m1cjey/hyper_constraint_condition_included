#pragma once
//
//���f���N���X
//

#include "Particle.h"	//�p������include
class Model : Particle
{
	vector<Particle> PART;
	vector<point3D> node;
	vector<element3D> nelm;
	Particle PART0;					//���������ꂽParticle�N���X�̃f�[�^
	ofstream fq;
	int point[4];

public:
	Model();
	~Model();
	//////////���q���͊w�֐�//////////
	virtual void MD_2D(double le,int BstartID,int beforeN,int newN);
	virtual void MD_3D(double le,int BstartID,int beforeN,int newN,double r,double region[3][2]);
		
	//////////�~���쐬�֐�//////////
	virtual int Set_cylinder(int initial_ID, double le, double R, double height);	
	//////////���쐬�֐�//////////
	virtual int Set_sphere_function(int initial_ID, double le, double R);
		//�~�쐬�֐�
		virtual int Set_circle_edge(int initial_number, double le, double R);
		//�~�����[�U�֐�
		virtual int Set_circle_in(int number, double le, double R, int initial_number);
		//�~�����[�U6�s�[�X�Ώ̊֐�
		virtual int Set_circle_in_using_6_pieces(int number,double le,double R,int  initial_ID);
		//���쐬�֐�
		virtual int Set_sphere(int number,double le,double R,int flag);
		//�~���\�ʍ��W�쐬
		virtual int Set_cylinder_face(int number,double le,double R,double height, int circle_start_id);
		//�~�������[�U�֐�
		virtual int Set_cylinder_in(int number,double le, double R, double height, int circle_start_id);
		//�����𕪊����邳���̍œK�ȕ������ƕ��������̎Z�o�֐�
		virtual void Set_calc_N_and_L(double dis,double le,int *N,double *L);
		//�~���������v�Z�֐�
		virtual int Set_calc_division_N_circle(double dis,double le);

	//////////�t�@�C���o�͊֐�//////////
	void writedata(ofstream &fp, int id, double x, double y,double z, int type,int materialID,int surface,double val,double vx,double vy,double vz,double P,double h,int toBEM);
	
	//////////���f���Z�b�g�֐�//////////
	virtual int Model_set();
		//////////�����f���֐�//////////
		virtual int Set_sphere_model();
		//////////�~�����f���֐�//////////
		virtual int Set_cylinder_model();
		//////////�����̃��f���֐�//////////
		virtual int Set_cube_model();
		//////////�������莎�����f��//////////
		virtual int Set_tensiontest_model();
			virtual void Set_point(int low_ID,int top_ID,int lefy_ID,int right_ID);
			virtual void Get_point(int &a,int &b,int &c,int &d);
		//////////�A�N�`���G�[�^���f���֐�//////////
		virtual int Set_actuator_model();
		//////////�x���`�}�[�N��胂�f���֐�//////////
		virtual int Set_benchmark_model();
	
};