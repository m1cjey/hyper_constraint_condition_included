#pragma once
//
//���q�N���X
//

class Particle 
{
	double X;	//���q���W
	double Y;
	double Z;

	double R;	//���q���a
	
public:
	

	virtual void Set(double x, double y, double z, double r);	//���W�Ɣ��a���Z�b�g
	virtual void Add(double x_add, double y_add, double z_add);	//�����Z
	virtual void Multiply(double x_mul, double y_mul, double z_mul);//�Ϗ��Z
	virtual double Get_X();
	virtual double Get_Y();
	virtual double Get_Z();
	
	
};