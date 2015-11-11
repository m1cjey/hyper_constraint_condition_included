#include "stdafx.h"
//dt�Ȃǂ̃O���[�o���ȏ����̕ύX��CON�ōs�����ƁI�I

elastic::elastic(vector<mpselastic> &PART)
{
	cout<<"�e���̌v�Z�p�̏��������v�Z--";

	//�`�F�b�N�p�t���O
	symplectic=ON;
	symplectic_order=1;
	FEM_flag=get_FEM_flag();//�ŏ���CON�̐ݒ�ɓ���

	pivot_check=OFF;

	
	mag_youngs_modulus=get_E_m();
	mag_poisson_ratio=get_v_m();
	elas_youngs_modulus=get_E_e();
	elas_poisson_ratio=get_v_e();
	le=get_distancebp();
	mass=get_particle_mass();
	
	mag_shear_modulus=mag_youngs_modulus/(2.0*(1.0+mag_poisson_ratio));
	mag_lambda=(mag_poisson_ratio*mag_youngs_modulus)/((1.0+mag_poisson_ratio)*(1.0-2.0*mag_poisson_ratio));

	elas_shear_modulus=elas_youngs_modulus/(2.0*(1.0+elas_poisson_ratio));
	elas_lambda=(elas_poisson_ratio*elas_youngs_modulus)/((1.0+elas_poisson_ratio)*(1.0-2.0*elas_poisson_ratio));

	inertia=mass*le*le/6;
	r=get_re_elastic()*le;

	if(get_modify_density()==ON)
	{
		if(get_model_set_way()==0) modify_density(0.523598775);
		else if(get_model_set_way()==1) modify_density(sqrt(2.0)*PI/6);//modify_density(0.724)��0.74�I�I
		else{cout<<"model set error"<<endl; exit(1);}
	}

	//�G�l���M�[�v�Z
	hamiltonian=0.0;
	potential=0.0;
	kinetic_energy=0.0;
	elastic_energy=0.0;

	last_elastic_energy=0.0;
	last_kinetic_energy=0.0;
	last_potential=0.0;

	cout<<"Lame constant lambda="<<mag_lambda<<" shear_modulus="<<mag_shear_modulus<<endl;
	cout<<"FEM flag: "<<boolalpha<<FEM_flag<<endl;
}

void elastic::check_elastic_config()
{
	cout<<"lambda ="<<mag_lambda<<endl;
	cout<<"shear modulus ="<<mag_shear_modulus<<endl;
	cout<<"poisson's ratio ="<<mag_poisson_ratio<<endl;
	cout<<"le ="<<le<<endl;
	cout<<"r ="<<r<<endl;
}

void elastic::reset_energy()
{
	hamiltonian=0.0;
	potential=0.0;
	kinetic_energy=0.0;
	elastic_energy=0.0;
}
	
void mapcheck(vector<mpselastic> &PART, int i)
{
	cout<<"PART["<<i<<"].N: "<<PART[i].N<<", PART["<<i<<"].distance0.size(): "<<PART[i].distance0.size()<<endl;
	map<int, double>::iterator it=PART[i].distance0.begin();
	int count=0;
	while(it!=PART[i].distance0.end())
	{
		cout<<count<<"; key: "<<(*it).first<<", value: "<<(*it).second<<endl;
		it++;
		count++;
	}
	char ch;
	cin>>ch;
}