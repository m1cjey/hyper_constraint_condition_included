//#define MPSPARTICLE
//#ifndef MPSPARTICLE //#infdef��
#include "stdafx.h"

//freeon�ł���Ă��邱�Ƃ͂��ׂĂ̗��q�̈ʒu��񂪑���Ȃ��Ǝ��s�ł��Ȃ��̂ŁA�R���X�g���N�^��C�j�V�����C�U�ł͎��s�ł��Ȃ�
//��̓I�ɂ�"PART[j]"�̏�񂪕K�v�Ȍv�Z�̓����o�֐��ɂł��Ȃ�
//ELAST�̃R���X�g���N�^�Ŏ��s�H

vector3D::vector3D(double *r0)
{
	component[0]=r0[0];
	component[1]=r0[1];
	component[2]=r0[2];	
}

void vector3D::operator=(double vc[3])
{
	component[0]=vc[0];
	component[1]=vc[1];
	component[2]=vc[2];
}

mpselastic::mpselastic():youngs_modulus(0),poisson_ratio(0)
{
	stop_on_floor=false;

	//�p�x�������E�p���x������
	for(int D=0;D<3;D++)
	{
		ang[D]=0.0;
		ang_u[D]=0.0;
	}
	ang[3]=1.0;
	cout<<"mpselastic constructor run"<<endl;
}

//freeon�̌�Ŏ��s����
void mpselastic::initialize_particles(elastic &ELAST, const int t)
{
	try{
		if(t==0)//�ی�p��if
		{
			for(int D=0;D<3;D++) initial_r[D]=r[D];
			PND0=PND;
			set_density(ELAST.get_density());
			initial_neighboursID=current_neighboursID;
			initial_neighbourX=current_neighbourX;
			initial_neighbourY=current_neighbourY;
			initial_neighbourZ=current_neighbourZ;

			initial_distancebps=current_distancebps;
			distance0=distance;//map�̃R�s�[

			if(N!=initial_neighboursID.size()) throw initial_neighboursID.size();
			else if(N!=initial_neighbourX.size()) throw initial_neighbourX.size();
			else if(N!=initial_neighbourY.size()) throw initial_neighbourY.size();
			else if(N!=initial_neighbourZ.size()) throw initial_neighbourZ.size();
			else if(N!=initial_distancebps.size()) throw initial_distancebps.size();

//			cout<<"particle initializer worked out"<<endl;
		}else{
			throw t;
		}
	}
	catch(size_t N)
	{
		if(t!=0)
		{
			cout<<"t="<<t<<"; EXIT_FAILURE"<<endl;
			exit(1);
		}else{
			cout<<"Neighbouring particles do not balance each other.\n";
			cout<<"original N:"<<initial_neighboursID.size()<<", different N"<<N<<endl;
			exit(1);
		}
	}

}

//vector<mpselastic>::PART[i]��������
void mpselastic::reset_particle_acceleration()
{
	P=0.0;
	contact=false;

	for(int D=0;D<3;D++)
	{
		pressure_accel[D]=0.0;//���ꂪ�Ȃ��Ɛ�΂ɔ��͂����������Ȃ�I
		stress_accel[D]=0.0;
		stress_visco_accel[D]=0.0;
		total_accel[D]=0.0;
	}

	//r0�̓N���A���Ă����Ȃ���mpselastic::set_r0_ij(double *r0)��push_back����Ă��������Ȃ�
	r0_ij.clear();
	r0_ji.clear();
}

void mpselastic::set_r0_ij(double *r0)
{
	vector3D vc(r0);
	r0_ij.push_back(vc);
}

void mpselastic::set_r0_ji(double *r0)
{
	vector3D vc(r0);
	r0_ji.push_back(vc);
	//vector<vector3D> r0_ji;
}

void mpselastic::check_acceleration()
{
//	cout<<"show acceleration"<<endl;
	cout<<"stress= "<<stress_accel[2]<<", st_visco= "<<stress_visco_accel[2]<<", P= "<<pressure_accel[2]<<endl;
}

void mpselastic::check_elastic_constant()
{
	cout<<"youngs_modulus: "<<youngs_modulus<<", poisson_ratio: "<<poisson_ratio<<endl;
}


//#endif