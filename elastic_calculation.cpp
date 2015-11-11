#include "stdafx.h"
#include "Rigidbody.h"
#include "Micro_AVS.h"
//#include<boost\array.hpp>
//#include <boost\shared_ptr.hpp>
using namespace std;

//�S�e���̌v�Z
//�G�l���M�[�͊֐�������Čv�Z������
void calc_elastic(vector<mpselastic> &PART, elastic &ELAST, int t, double **F)
{
	cout<<"���q�ړ��v�Z�J�n"<<endl;

	int dimension=ELAST.get_dimension();

/*	if(t==1){
		rigid_calc(PART, ELAST);
	}*/

	if(ELAST.get_modify_density()==ON) calc_modified_density(PART, ELAST); //���x�̏C��

	if(ELAST.get_nonlinear_elastic_flag()==false)
	{
		calc_pre_velocity_and_position(PART, ELAST, F);//�O��̃X�e�b�v�̓��́i�����x�j��p���ĉ��̈ʒu�E���x���v�Z

		#pragma omp parallel for
		//���̂̌v�Z
//		calc_rigid(PART,rigids);

		for(int i=0;i<PART.size();i++) PART[i].reset_particle_acceleration(); //�����x�̃��Z�b�g

		calc_quaternion_using_vectorSTL(ELAST, PART);
		calc_angular_velocity_using_vectorSTL(ELAST, PART);

		//�����ʒu�x�N�g���̉�]�i���݈ʒu�x�N�g���͉�]�����Ȃ��Ă悢�j
		calc_r0_ij(PART); 

		//
		calc_accel_for_3D(PART, ELAST); 

		//���͂ƐڐG���͂ɂ������x���v�Z
		contact_judge(PART, ELAST);

		//�e�����x�̍��v
		modify_acceleration(PART, ELAST, F);

		//���̈ʒu�E���x���狁�߂����́i�����x�j��u, r���C��
		calc_post_velocity_and_position(PART, ELAST,t);




/*		for(int i=0;i<PART.size();i++){
		cout<<"i="<<i<<", rx="<<PART[i].r[A_X]<<", ry="<<PART[i].r[A_Y]<<", rz="<<PART[i].r[A_Z]<<endl;//
	}//*/
		//�n�~���g�j�A���̃t�@�C���o��
		calc_hamiltonian(ELAST, t);
	}
	else
	{
		calc_pre_velocity_and_position(PART, ELAST, F);//�O��̃X�e�b�v�̓��́i�����x�j��p���ĉ��̈ʒu�E���x���v�Z

		#pragma omp parallel for
		for(int i=0;i<PART.size();i++) PART[i].reset_particle_acceleration(); //�����x�̃��Z�b�g

		calc_quaternion_using_vectorSTL(ELAST, PART);
		calc_angular_velocity_using_vectorSTL(ELAST, PART);

		//�����ʒu�x�N�g���̉�]
		calc_r0_ij(PART); 


		calc_nonlinear_accel_for_3D_test(PART, ELAST);

		//���͂ƐڐG���͂ɂ������x���v�Z
		contact_judge(PART, ELAST);	

		//�e�����x�̍��v
		modify_acceleration(PART, ELAST, F);

		//���̈ʒu�E���x���狁�߂����́i�����x�j��u, r���C��
		calc_post_velocity_and_position(PART, ELAST,t);

		//���̂̈ʒu�C��
//		calc_rigid(PART,rigids);

		//�n�~���g�j�A���̃t�@�C���o��
		calc_hamiltonian(ELAST, t);
	}
}

void calc_accel_for_3D(vector<mpselastic> &PART, elastic &ELAST)
{
	double dimension=static_cast<double>(ELAST.get_dimension());
	double mag_shear_modulus=ELAST.get_mag_shear_modulus();
	double mag_lambda=ELAST.get_mag_lambda();
	double elas_shear_modulus=ELAST.get_elas_shear_modulus();
	double elas_lambda=ELAST.get_elas_lambda();
	double mass=ELAST.get_mass();
	double density=ELAST.get_density();
	double re=ELAST.get_r(); //����͊���re*le
	double vis=ELAST.get_nensei(); //�S�x�i���S�x�ł͂Ȃ��j
	double g=ELAST.get_g();

	double KE=0.0;	//�^���G�l���M�[
	double EE1=0.0;	//�Ђ��݃G�l���M�[
	double EE2=0.0;	//�̐ςЂ��݂ɂ��
	double PE=0.0;	//�|�e���V����
	double ground=ELAST.get_ground_position(); //���̂����W�̎擾

	cout<<"ground: "<<ground<<endl;
	double mindis=100;
	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{   
		//WALL�͏��O���Ă͂����Ȃ�

		double EE1_temp=0.0;

		//�^���G�l���M�[�̍X�V
		if(PART[i].type==ELASTIC || PART[i].type==MAGELAST){
			for(int D=0;D<3;D++) KE+=PART[i].u_temp[D]*PART[i].u_temp[D];
		}

		//�ʒu�|�e���V�����̍X�V
		if(PART[i].type==ELASTIC || PART[i].type==MAGELAST){
			PE+=(PART[i].r[2]-ground);
		}

		int neighboursN0=static_cast<int>(PART[i].get_initial_neighboursID().size());

		//////////////////�����z�u�Ŏ��ӂɂ�����̂�����ǐՂ���E�E�E/////////////////////////
		for(int k=0;k<neighboursN0;k++)
		{
			int j=PART[i].get_initial_neighboursID()[k];

			//i����݂������z�u�ł̑��΍��W
			double r_ij_Init[3], r_ji_Init[3]; 
			for(int D=0;D<3;D++){
				r_ij_Init[D]=PART[j].initial_r[D]-PART[i].initial_r[D];
				r_ji_Init[D]=-r_ij_Init[D];
			}

			//���ݔz�u�ł̑��΍��W
			double r_ij[3], r_ji[3];
			for(int D=0;D<3;D++)
			{
				r_ij[D]=PART[j].r_temp[D]-PART[i].r_temp[D];
				r_ji[D]=-r_ij[D];
			}

			double dis=0.0; for(int D=0;D<3;D++){dis+=r_ij[D]*r_ij[D];} dis=sqrt(dis); //���ݗ��q�ԋ���
			double dis0=0.0; for(int D=0;D<3;D++){dis0+=r_ij_Init[D]*r_ij_Init[D];} dis0=sqrt(dis0); //�������q�ԋ���
			if(dis<=mindis) mindis=dis;
			double r_ij_zero[3], r_ji_zero[3];

			//calc_r0_ij(PART)��"���݈ʒu�̃N�H�[�^�j�I�����g����"��]�����������ʒu�x�N�g�����擾
			for(int D=0;D<3;D++) r_ij_zero[D]=PART[i].get_r0_ij()[k].get_comp()[D];
			for(int D=0;D<3;D++) r_ji_zero[D]=PART[i].get_r0_ji()[k].get_comp()[D];

			double w=kernel(re, dis0);

			//��]��̏����z�u���Έʒu�̒P�ʃx�N�g��
			double n_ij_zero[3], n_ji_zero[3];
			for(int D=0;D<3;D++) 
			{
				n_ij_zero[D]=r_ij_zero[D]/dis0;
				n_ji_zero[D]=r_ji_zero[D]/dis0;
			}

			//���ݔz�u���Έʒu�P�ʃx�N�g��
			double n_ij_current[3], n_ji_current[3];
			for(int D=0;D<3;D++) 
			{
				n_ij_current[D]=r_ij[D]/dis;
				n_ji_current[D]=r_ji[D]/dis;
			}

			//�ψʃx�N�g��
			double U_ij[3], U_ji[3]; 
			for(int D=0;D<3;D++)
			{
				U_ij[D]=r_ij[D]-r_ij_zero[D];
				U_ji[D]=r_ji[D]-r_ji_zero[D];
			}

			//�Ђ��݃x�N�g��
			double E_ij[3], E_ji[3];
			for(int D=0;D<3;D++)
			{
				E_ij[D]=U_ij[D]/dis0;
				E_ji[D]=U_ji[D]/dis0;
			}
				
			//���͌v�Z�̏����E�E�Ee_vol_i�̃�[]�i�ψʂ̔��U�j���v�Z
			//�����Ђ��݂̘a�����i���̓e���\���̑Ίp�����̘a�ɏd�ݕt���j
			double volumetric_strain=0.0;
			for(int D=0;D<3;D++) volumetric_strain+=E_ij[D]*n_ij_current[D];
//			for(int D=0;D<3;D++) PART[i].P+=E_ij[D]*n_ij_zero[D]*w; 
			volumetric_strain*=w;
			PART[i].P+=volumetric_strain; //���ł͂Ȃ��̂Œ��ӁI

			//�Ђ��݃G�l���M�[
			for(int D=0;D<3;D++) EE1_temp+=E_ij[D]*E_ij[D]*w;

			//���̓x�N�g���i��舵����ύX�E�E�E�d�݊֐��͗��qi�ŏ�ɉ��d���ς���B�������Ȃ��ƈ��͂̏d�݊֐��Ƃ̐����������Ȃ��j
			double sigma_ij[3], sigma_ji[3];
			for(int D=0;D<3;D++)
			{
				if(PART[i].type==MAGELAST){
				sigma_ij[D]=mag_shear_modulus*w*(E_ij[D]-E_ji[D])/dis0;
				}
				else if(PART[i].type==ELASTIC){
					sigma_ij[D]=elas_shear_modulus*w*(E_ij[D]-E_ji[D])/dis0;
				}
			}
			PART[i].add_stress_accel(sigma_ij);

//			if(fabs((dis-dis0)/dis0)>1.04) w=0;//�j�����

			//�Ђ��ݑ��x�̌v�Z
			double strain_vi[3], strain_vj[3];

			//i����݂��Ђ��ݑ��x
			strain_vi[0]=((PART[j].u_temp[0]-PART[i].u_temp[0])-(PART[i].ang_u[1]*r_ij[2]-PART[i].ang_u[2]*r_ij[1]))/dis; //X������̂Ђ��ݑ��x ��������ˉe����K�v����H�H
			strain_vi[1]=((PART[j].u_temp[1]-PART[i].u_temp[1])-(PART[i].ang_u[2]*r_ij[0]-PART[i].ang_u[0]*r_ij[2]))/dis; //Y������̂Ђ��ݑ��x
			strain_vi[2]=((PART[j].u_temp[2]-PART[i].u_temp[2])-(PART[i].ang_u[0]*r_ij[1]-PART[i].ang_u[1]*r_ij[0]))/dis; //Z������̂Ђ��ݑ��x
			//anglar_u1..3�͔z��ɏ�������

			//j����݂��Ђ��ݑ��x
			strain_vj[0]=((PART[i].u_temp[0]-PART[j].u_temp[0])-(PART[j].ang_u[1]*r_ji[2]-PART[j].ang_u[2]*r_ji[1]))/dis; //X������̂Ђ��ݑ��x ��������ˉe����K�v����H�H
			strain_vj[1]=((PART[i].u_temp[1]-PART[j].u_temp[1])-(PART[j].ang_u[2]*r_ji[0]-PART[j].ang_u[0]*r_ji[2]))/dis; //Y������̂Ђ��ݑ��x
			strain_vj[2]=((PART[i].u_temp[2]-PART[j].u_temp[2])-(PART[j].ang_u[0]*r_ji[1]-PART[j].ang_u[1]*r_ji[0]))/dis; //Z������̂Ђ��ݑ��x
			
			double sigma_v_ij[3], sigma_v_ji[3];
			for(int D=0;D<3;D++)
			{
				sigma_v_ij[D]=vis*w*(strain_vi[D]-strain_vj[D])/dis;
//				sigma_v_ji[D]=2*vis*w*strain_vj[D]/PART[j].PND/dis;
//				sigma_v_ij[D]-=sigma_v_ji[D];
			}

			PART[i].add_stress_visco_accel(sigma_v_ij);
		}//for(int k=0;k<neighboursN0;k++)���[�v�I�� �������q�z�u�ł̎��͂ɂ��闱�q
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//		PART[i].P*=dimension/PART[i].PND;//�̐ςЂ��݂����߂�ꂽ�i����̓G�l���M�[�v�Z����O�ɋ��߂Ă������Ɓj
		PART[i].P*=dimension/PART[i].PND0;//�̐ςЂ��݂����߂�ꂽ�i����̓G�l���M�[�v�Z����O�ɋ��߂Ă������Ɓj

//		double coef=dimension/PART[i].get_density()/PART[i].PND;
		double coef=dimension/PART[i].get_density()/PART[i].PND0;

		PART[i].mul_stress_accel(coef);
		PART[i].mul_stress_visco_accel(coef);

		//�Ђ��݃G�l���M�[
		if(PART[i].type==ELASTIC || PART[i].type==MAGELAST || PART[i].type==TERMINAL1 || PART[i].type==TERMINAL2){
//			EE1+=EE1_temp/PART[i].PND/PART[i].get_density();//��񍀂̃� ���x�͂��ꂼ��Ⴄ
	//		if(PART[i].type==MAGELAST || PART[i].type==TERMINAL1 || PART[i].type==TERMINAL2){
			if(PART[i].type==MAGELAST){
			EE1+=(EE1_temp/PART[i].PND0/PART[i].get_density())*mag_shear_modulus*mass*dimension;//��񍀂̃� ���x�͂��ꂼ��Ⴄ
			EE2+=((PART[i].P*PART[i].P)/PART[i].get_density())*0.5*mass*mag_lambda;//��O���̃� ���x�͂��ꂼ��Ⴄ
			}
			else if(PART[i].type==ELASTIC){
			EE1+=(EE1_temp/PART[i].PND0/PART[i].get_density())*elas_shear_modulus*mass*dimension;//��񍀂̃� ���x�͂��ꂼ��Ⴄ
			EE2+=((PART[i].P*PART[i].P)/PART[i].get_density())*0.5*mass*elas_lambda;//��O���̃� ���x�͂��ꂼ��Ⴄ
			}
		}

//		if(PART[i].type==MAGELAST || PART[i].type==TERMINAL1 || PART[i].type==TERMINAL2){
		if(PART[i].type==MAGELAST || PART[i].type==TERMINAL1 || PART[i].type==TERMINAL2){
		PART[i].P*=-mag_lambda;//�̐ςЂ��݂ɂ�鈳�͂����߂�ꂽ
		}
		else if(PART[i].type==ELASTIC){
		PART[i].P*=-elas_lambda;
		}
		else if(PART[i].type==WALL){
		PART[i].P*=-elas_lambda;
		}
//		����Œu����������WALL���Ōv�Z���ď�̃u���b�N�Ŋ֐��ŗL���𔻒肵�ĕ������t�]����������m���ł́E�E�E
		
		//���q�����x�����������ꍇ��P��u��������
		if(PART[i].PND>PART[i].PND0)
		{
			if(PART[i].type==MAGELAST || PART[i].type==TERMINAL1 || PART[i].type==TERMINAL2)
			{
				double contact_pressure=mag_lambda*(PART[i].PND-PART[i].PND0)/PART[i].PND0;
				if(contact_pressure>PART[i].P)
				{
					PART[i].P=contact_pressure;//���x������`�Ȃ̂ɂ����ŗ��q�����x�����`�Ƃ��ėǂ��H��OK�in0�ƃς͔��j
				}
				PART[i].contact=true; //�t���O�𗧂ĂĂ���
			}
			else if(PART[i].type==ELASTIC){
			double contact_pressure=elas_lambda*(PART[i].PND-PART[i].PND0)/PART[i].PND0;
			if(contact_pressure>PART[i].P)
			{
				PART[i].P=contact_pressure;//���x������`�Ȃ̂ɂ����ŗ��q�����x�����`�Ƃ��ėǂ��H��OK�in0�ƃς͔��j
			}
			PART[i].contact=true; //�t���O�𗧂ĂĂ���
			}
			else if(PART[i].type==WALL){
			double contact_pressure=elas_lambda*(PART[i].PND-PART[i].PND0)/PART[i].PND0;
			if(contact_pressure>PART[i].P)
			{
				PART[i].P=contact_pressure;//���x������`�Ȃ̂ɂ����ŗ��q�����x�����`�Ƃ��ėǂ��H��OK�in0�ƃς͔��j
			}
			PART[i].contact=true; //�t���O�𗧂ĂĂ���
			}
		}


	}//for(i=0;i<PART.size();i++)���[�v�I��
	cout<<"�ŏ����q�ԋ���"<<mindis<<endl;
	KE*=0.5*mass;
	PE*=-g*mass;
	
//	EE1*=shear_modulus*mass*dimension;
//	EE2*=0.5*mass*lambda;
	ELAST.set_kinetic(KE);
	ELAST.set_elastic_energy(EE1+EE2);//�S�n�̒e���G�l���M�[
	ELAST.set_potential(PE);
	ELAST.set_hamiltonian(KE+EE1+EE2+PE);
}

void calc_nonlinear_accel_for_3D_test(vector<mpselastic> &PART, elastic &ELAST)
{
	double dimension=static_cast<double>(ELAST.get_dimension());
	double mag_shear_modulus=ELAST.get_mag_shear_modulus();
	double elas_shear_modulus=ELAST.get_elas_shear_modulus();
	double mass=ELAST.get_mass();
	double density=ELAST.get_density();
	double re=ELAST.get_r(); //����͊���re*le
	double vis=ELAST.get_nensei(); //�S�x�i���S�x�ł͂Ȃ��j
	double g=ELAST.get_g();

//	vector<double> lambda;

	double KE=0.0;	//�^���G�l���M�[
	double EE1=0.0;	//�Ђ��݃G�l���M�[
	double EE2=0.0;	//�̐ςЂ��݂ɂ��
	double PE=0.0;	//�|�e���V����
	double ground=ELAST.get_ground_position(); //���̂����W�̎擾


//	cout<<"ground: "<<ground<<endl;
	
	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{   
		if(PART[i].type==ELASTIC || PART[i].type==MAGELAST || PART[i].type==TERMINAL1 || PART[i].type==TERMINAL2)//WALL�͒e���ό`���Ȃ��̂œ��͌v�Z�Ɋ܂߂Ȃ�
		{

			double EE1_temp=0.0;
			double EE2_temp=0.0;
			double volumetric_strain=0.0;//�̐ςЂ��݁ii���ƂɈقȂ�j

			//�^���G�l���M�[�̍X�V
			if(PART[i].type==ELASTIC || PART[i].type==MAGELAST){
				for(int D=0;D<3;D++) KE+=PART[i].u[D]*PART[i].u[D];
			}

			//�ʒu�|�e���V�����̍X�V
			if(PART[i].type==ELASTIC || PART[i].type==MAGELAST){
				PE+=(PART[i].r[2]-ground);
			}

			int neighboursN0=static_cast<int>(PART[i].get_initial_neighboursID().size());

			//���ڂ������qi�Ǝ��ӗ��qj�Ƃ́u�֌W�v���v�Z
			//�����z�u�Ŏ��ӂɂ�����̂�����ǐՂ���E�E�E
			for(int k=0;k<neighboursN0;k++) //pressure�ƍ��킹�Ȃ��Ƌ����H
			{
				int j=PART[i].get_initial_neighboursID()[k];

				//i����݂������z�u�ł̑��΍��W
				double r_ij_Init[3], r_ji_Init[3]; 
				for(int D=0;D<3;D++){
					r_ij_Init[D]=PART[j].initial_r[D]-PART[i].initial_r[D];
					r_ji_Init[D]=-r_ij_Init[D];
				}

				//���ݔz�u�ł̑��΍��W
				double r_ij[3], r_ji[3];
				for(int D=0;D<3;D++)
				{
					r_ij[D]=PART[j].r_temp[D]-PART[i].r_temp[D];
					r_ji[D]=-r_ij[D];
				}

				double dis=0.0; for(int D=0;D<3;D++){dis+=r_ij[D]*r_ij[D];} dis=sqrt(dis); //���ݗ��q�ԋ���
				double dis0=0.0; for(int D=0;D<3;D++){dis0+=r_ij_Init[D]*r_ij_Init[D];} dis0=sqrt(dis0); //�������q�ԋ���
			
				double r_ij_zero[3], r_ji_zero[3];

				//calc_r0_ij(PART)��"���݈ʒu�̃N�H�[�^�j�I�����g����"��]�����������ʒu�x�N�g�����擾
				for(int D=0;D<3;D++) r_ij_zero[D]=PART[i].get_r0_ij()[k].get_comp()[D];
				for(int D=0;D<3;D++) r_ji_zero[D]=PART[i].get_r0_ji()[k].get_comp()[D];

				double w=kernel(re, dis0);

				//��]��̏����z�u���Έʒu�̒P�ʃx�N�g��
				double n_ij_zero[3], n_ji_zero[3];
				for(int D=0;D<3;D++) 
				{
					n_ij_zero[D]=r_ij_zero[D]/dis0;
					n_ji_zero[D]=r_ji_zero[D]/dis0;
				}

				//���ݔz�u���Έʒu�P�ʃx�N�g��
				double n_ij_current[3], n_ji_current[3];
				for(int D=0;D<3;D++) 
				{
					n_ij_current[D]=r_ij[D]/dis;
					n_ji_current[D]=r_ji[D]/dis;
				}

				//�ψʃx�N�g��
				double U_ij[3], U_ji[3]; 
				for(int D=0;D<3;D++)
				{
					U_ij[D]=r_ij[D]-r_ij_zero[D];
					U_ji[D]=r_ji[D]-r_ji_zero[D];
				}

				//�Ђ��݃x�N�g��
				double E_ij[3], E_ji[3];
				for(int D=0;D<3;D++)
				{
					E_ij[D]=U_ij[D]/dis0;
					E_ji[D]=U_ji[D]/dis0;
				}

				//�����O���ƃ|�A�\����̌v�Z
				double En_ij[3];
				double volumetric_temp=0.0;
				for(int D=0;D<3;D++) volumetric_temp+=E_ij[D]*n_ij_current[D];	//E_ij��n_ij_current�̓���
				for(int D=0;D<3;D++) En_ij[D]=volumetric_temp*n_ij_current[D];	//�Ђ��݃x�N�g����r_ij�����֐��ˉe
				double ns=sqrt(En_ij[0]*En_ij[0]+En_ij[1]*En_ij[1]+En_ij[2]*En_ij[2]); //�c�Ђ���
				ns*=100;	//[%]
				//�V���R�[���̏ꍇ
				//youngs_modulus=-6.11619E-11*x^9+1.66522E-8*x^8-1.93253E-6*x^7+1.2459E-4*x^6-4.87141E-3*x^5+0.118275*x^4-1.75252*x^3+14.9237*x^2-64.3273*x+123.411;
				double elas_youngs_modulus=(((((((-0.000004092500180)*ns+0.000860882692091)*ns-0.073155040824811)*ns+3.215486271963979)*ns-77.712461766914529)*ns+1012.872595336284500)*ns-6495.847221601891800)*ns+18812.474059104068000;
				//-1.09176E-12*x^10+1.56272E-10*x^9-9.70758E-9*x^8+3.42745E-7*x^7-7.56378E-6*x^6+0.00010807*x^5+-0.00100264*x^4+5.93234E-3*x^3-0.0219534*x^2+0.0537726*x^1+0.061154
				double elas_poisson_ratio=(((((((+0.000000000918578)*ns-0.000000090882723)*ns+0.000003665373271)*ns-0.000078101785064)*ns+0.000960394839085)*ns-0.007112271340844)*ns+0.033194941391016)*ns+0.069899017991200
;

				//MRE�̏ꍇ
				//youngs_modulus=-4.91157E-8*x^7+9.42896E-6*x^6-7.31665E-4*x^5+2.95145E-2*x^4-0.663258*x^3+8.25586*x^2-53.0674*x+190.48
				double mag_youngs_modulus=(((((((-0.000002526899825)*ns+0.000526467419950)*ns-0.044510412964222)*ns+1.962528096212004)*ns-48.344643329754824)*ns+664.648197882870590)*ns-4894.930152042784400)*ns+23399.990870638816000;
				//-9.6723E-12*x^10+1.22758E-9*x^9-6.73685E-8*x^8+2.09375E-7*x^7-4.05552E-5*x^6+5.08152E-4*x^5-0.00414923*x^4+0.0218461*x^3-0.072922*x^2+0.151892*x^1+0.17713
				double mag_poisson_ratio=(((((((+0.000000006114133)*ns-0.000000557177269)*ns+0.000020762913954)*ns-0.000408515956157)*ns+0.004578221540415)*ns-0.029448326315701)*ns+0.103023742143433)*ns+0.192042424941781;

				if(i==1511 && j==1513){
					ofstream pr("poason.dat", ios::app);
					cout<<"youngs_modulus: "<<mag_youngs_modulus<<", poisson_ratio: "<<mag_poisson_ratio<<endl;
					pr<<"youngs_modulus="<<mag_youngs_modulus<<", poisson_ratio="<<mag_poisson_ratio<<", strain="<<ns<<endl;
					pr.close();
				}
				//�����萔mu
				double elas_shear_modulus=elas_youngs_modulus/(2.0*(1.0+elas_poisson_ratio));
				double mag_shear_modulus=mag_youngs_modulus/(2.0*(1.0+mag_poisson_ratio));
				//�����萔lambda
				double elas_lambda_temp=(elas_poisson_ratio*elas_youngs_modulus)/((1.0+elas_poisson_ratio)*(1.0-2.0*elas_poisson_ratio));
				double mag_lambda_temp=(mag_poisson_ratio*mag_youngs_modulus)/((1.0+mag_poisson_ratio)*(1.0-2.0*mag_poisson_ratio));
				
	//			cout<<"shear_modulus: "<<shear_modulus<<", lambda: "<<lambda<<endl;

				//���͌v�Z�̏����E�E�Ee_vol_i�̃�[]�i�ψʂ̔��U�j���v�Z
				//�����Ђ��݂̘a�����i���̓e���\���̑Ίp�����̘a�ɏd�ݕt���j
				if(PART[i].type==ELASTIC){
				volumetric_temp*=-elas_lambda_temp*w;
				PART[i].P+=volumetric_temp;
				PART[i].ave_lambda+=elas_lambda_temp;
				//�Ђ��݃G�l���M�[(�d�݂͂���Ȃ�)
				for(int D=0;D<3;D++) EE1_temp+=E_ij[D]*E_ij[D]*w;
				EE2_temp+=elas_lambda_temp*volumetric_temp*volumetric_temp;

				//���̓x�N�g���i��舵����ύX�E�E�E�d�݊֐��͗��qi�ŏ�ɉ��d���ς���Belse, ���͂̏d�݊֐��Ƃ̐����������Ȃ��j
				double sigma_ij[3], sigma_ji[3];
				for(int D=0;D<3;D++) sigma_ij[D]=elas_shear_modulus*w*(E_ij[D]-E_ji[D])/dis0;
				PART[i].add_stress_accel(sigma_ij);
				}
				///////////////////////////////////////////////////

				else if(PART[i].type==MAGELAST|| PART[i].type==TERMINAL1 || PART[i].type==TERMINAL2){
				volumetric_temp*=-mag_lambda_temp*w;
				PART[i].P+=volumetric_temp;
				PART[i].ave_lambda+=mag_lambda_temp;
				//�Ђ��݃G�l���M�[(�d�݂͂���Ȃ�)
				for(int D=0;D<3;D++) EE1_temp+=E_ij[D]*E_ij[D]*w;
		//		EE2_temp+=mag_lambda_temp*volumetric_temp*volumetric_temp;
				EE2_temp+=mag_lambda_temp*volumetric_temp*volumetric_temp/w*w; //�d�݂����Ȃ����

				//���̓x�N�g���i��舵����ύX�E�E�E�d�݊֐��͗��qi�ŏ�ɉ��d���ς���Belse, ���͂̏d�݊֐��Ƃ̐����������Ȃ��j
				double sigma_ij[3], sigma_ji[3];
				for(int D=0;D<3;D++) sigma_ij[D]=mag_shear_modulus*w*(E_ij[D]-E_ji[D])/dis0;
				PART[i].add_stress_accel(sigma_ij);
				}
				

	//			if(fabs((dis-dis0)/dis0)>1.04) w=0;//�j�����

				//�Ђ��ݑ��x�̌v�Z
				double strain_vi[3], strain_vj[3];

				//i����݂��Ђ��ݑ��x
				strain_vi[0]=((PART[j].u[0]-PART[i].u[0])-(PART[i].ang_u[1]*r_ij[2]-PART[i].ang_u[2]*r_ij[1]))/dis; //X������̂Ђ��ݑ��x ��������ˉe����K�v����H�H
				strain_vi[1]=((PART[j].u[1]-PART[i].u[1])-(PART[i].ang_u[2]*r_ij[0]-PART[i].ang_u[0]*r_ij[2]))/dis; //Y������̂Ђ��ݑ��x
				strain_vi[2]=((PART[j].u[2]-PART[i].u[2])-(PART[i].ang_u[0]*r_ij[1]-PART[i].ang_u[1]*r_ij[0]))/dis; //Z������̂Ђ��ݑ��x
				//anglar_u1..3�͔z��ɏ�������

				//j����݂��Ђ��ݑ��x
				strain_vj[0]=((PART[i].u[0]-PART[j].u[0])-(PART[j].ang_u[1]*r_ji[2]-PART[j].ang_u[2]*r_ji[1]))/dis; //X������̂Ђ��ݑ��x ��������ˉe����K�v����H�H
				strain_vj[1]=((PART[i].u[1]-PART[j].u[1])-(PART[j].ang_u[2]*r_ji[0]-PART[j].ang_u[0]*r_ji[2]))/dis; //Y������̂Ђ��ݑ��x
				strain_vj[2]=((PART[i].u[2]-PART[j].u[2])-(PART[j].ang_u[0]*r_ji[1]-PART[j].ang_u[1]*r_ji[0]))/dis; //Z������̂Ђ��ݑ��x
			
				//�S�����̓x�N�g��
				double sigma_v_ij[3], sigma_v_ji[3];
				for(int D=0;D<3;D++)
				{
					sigma_v_ij[D]=vis*w*(strain_vi[D]-strain_vj[D])/dis;
	//				sigma_v_ji[D]=vis*w*strain_vj[D]/PART[j].PND/dis;
	//				sigma_v_ij[D]-=sigma_v_ji[D];
				}

				PART[i].add_stress_visco_accel(sigma_v_ij);
			}//for(int k=0;k<neighboursN0;k++)���[�v�I��

			//�ɂ̕���
			PART[i].ave_lambda/=neighboursN0;	//�����܂ŏ������q�Ԃ̃�
			//����
			PART[i].P*=dimension/PART[i].PND0;

	//		double coef=dimension/PART[i].get_density()/PART[i].PND;
			double coef=dimension/PART[i].get_density()/PART[i].PND0;

			PART[i].mul_stress_accel(coef);
			PART[i].mul_stress_visco_accel(coef);

			//�Ђ��݃G�l���M�[
			if(PART[i].type==ELASTIC || PART[i].type==MAGELAST || PART[i].type==TERMINAL1 || PART[i].type==TERMINAL2){
				if(PART[i].type==MAGELAST){
				EE1+=(EE1_temp/PART[i].PND0/PART[i].get_density())*mag_shear_modulus*mass*dimension;//��񍀂̃� ���x�͂��ꂼ��Ⴄ
				}
				else if(PART[i].type==ELASTIC){
				EE1+=(EE1_temp/PART[i].PND0/PART[i].get_density())*elas_shear_modulus*mass*dimension;//��񍀂̃� ���x�͂��ꂼ��Ⴄ
				}
				EE2+=EE2_temp/PART[i].get_density();//��O���̃� ���x�͂��ꂼ��Ⴄ
			}

			//���q�����x�����������ꍇ��P��u��������
	

			if(PART[i].PND>PART[i].PND0)
		{
			if(PART[i].type==MAGELAST || PART[i].type==TERMINAL1||PART[i].type==TERMINAL2){
			double contact_pressure=ELAST.get_mag_lambda()*(PART[i].PND-PART[i].PND0)/PART[i].PND0;
			if(contact_pressure>PART[i].P)
			{
				PART[i].P=contact_pressure;//���x������`�Ȃ̂ɂ����ŗ��q�����x�����`�Ƃ��ėǂ��H��OK�in0�ƃς͔��j
			}
			PART[i].contact=true; //�t���O�𗧂ĂĂ���
			}
			else if(PART[i].type==ELASTIC){
			double contact_pressure=ELAST.get_elas_lambda()*(PART[i].PND-PART[i].PND0)/PART[i].PND0;
			if(contact_pressure>PART[i].P)
			{
				PART[i].P=contact_pressure;//���x������`�Ȃ̂ɂ����ŗ��q�����x�����`�Ƃ��ėǂ��H��OK�in0�ƃς͔��j
			}
			PART[i].contact=true; //�t���O�𗧂ĂĂ���
			}
			
		}
		}
	}//for(i=0;i<PART.size();i++)���[�v�I��

	//�ǂ���̔���
	#pragma omp parallel for
	for(int i=0;i<PART.size();i++){
		if(PART[i].type==WALL){
			int neighbours=static_cast<int>(PART[i].get_current_neighboursID().size());
			for(int k=0;k<neighbours;k++) //pressure�ƍ��킹�Ȃ��Ƌ����H
			{
				int j=PART[i].get_current_neighboursID()[k];
				if(PART[i].type!=WALL){ //�ǂ̉e�����a���ɓ������f�ނ̈��͂̍��v��ǂ���̔��͂Ƃ���B
				PART[i].P+=PART[j].P;
				}
			}
		}
	}
	

	KE*=0.5*mass;
	PE*=-g*mass;
//	EE1*=shear_modulus*mass*dimension;
	EE2*=0.5*mass;
	ELAST.set_kinetic(KE);
	ELAST.set_elastic_energy(EE1+EE2);//�S�n�̒e���G�l���M�[
	ELAST.set_elastic_energy1(EE1);
	ELAST.set_elastic_energy2(EE2);
	ELAST.set_potential(PE);
	ELAST.set_hamiltonian(KE+EE1+EE2+PE);
}

void contact_judge(vector<mpselastic> &PART, elastic &ELAST)
{
	//�A���S���Y��
	// 0. i���ӂ̗��q�����x�����������ꍇ�C�e�����a���ɂ��闱�q��T�����C�ȉ����s��
	// 1. �u�ڐG�̉\�������闱�q�v�i(PART[j].PND>PART[j].PND0)���^�H�j�𒲂ׂ�
	// 2. ���͂�u��
	// 3. �����z�u�̗��q�Əd�����Ȃ��悤�ɐڐG�̉\�������闱�q�Ƃ̊Ԃŗ͂��v�Z����
	double dimension=static_cast<double>(ELAST.get_dimension());
	double re=ELAST.get_r();
	double le=ELAST.get_distancebp();
	double ground=ELAST.get_ground_position();

	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{
//		if(PART[i].type==ELASTIC || PART[i].type==MAGELAST)
		{
			//�����z�u�̗��q�Ƃ͕��ʂɈ��͌��z���v�Z�i���͌v�Z�j
			if(PART[i].contact==false)
			{
				size_t neighbourN0=PART[i].get_initial_neighboursID().size();//�����z�u�ł̎��ӗ��q�����擾
				double qi[4]={PART[i].ang[0], PART[i].ang[1], PART[i].ang[2], PART[i].ang[3]};

				for(int k=0;k<neighbourN0;k++)
				{
					int j=PART[i].get_initial_neighboursID()[k];

					//���ݔz�u�ł̑��΍��W
					double r_ij[3], r_ji[3];
					for(int D=0;D<3;D++)
					{
						r_ij[D]=PART[j].r_temp[D]-PART[i].r_temp[D];
						r_ji[D]=-r_ij[D];
					}

					//���ݗ��q�ԋ���
					double dis=0.0; for(int D=0;D<3;D++){dis+=r_ij[D]*r_ij[D];} dis=sqrt(dis);

					//�e���̂̓��͂ł͏d�ݕt����dis0���g���ׂ��i�����z�u����̕ό`���d�v�Ȃ̂Łj
					double r_ij_Init[3], r_ji_Init[3]; 
					for(int D=0;D<3;D++)
					{
						r_ij_Init[D]=PART[j].initial_r[D]-PART[i].initial_r[D];
						r_ji_Init[D]=-r_ij_Init[D];
					}

					double r_ij_zero[3], r_ji_zero[3];
					//calc_r0_ij(PART)��"���݈ʒu�̃N�H�[�^�j�I�����g����"��]�����������ʒu�x�N�g�����擾
					for(int D=0;D<3;D++) r_ij_zero[D]=PART[i].get_r0_ij()[k].get_comp()[D];
					for(int D=0;D<3;D++) r_ji_zero[D]=PART[i].get_r0_ji()[k].get_comp()[D];

					double dis0=0.0; for(int D=0;D<3;D++){dis0+=r_ij_Init[D]*r_ij_Init[D];} dis0=sqrt(dis0); //�������q�ԋ���

					double w=kernel(re, dis0);

					double press_accel[3];
					for(int D=0;D<3;D++){
						//�d�݊֐��̉��d���ς̂Ƃ����ύX�Bi�̗��q�����x�Ōv�Z����
						//�d�݂͏����̋������g�������z�͌��ݔz�u�̕����g��
						//press_accel[D]=(PART[i].P*r_ij[D]-PART[j].P*r_ji[D])*w/dis0/dis0;
						//���z��dis0���g���Ă͂����Ȃ��H
						press_accel[D]=(PART[i].P*r_ij_zero[D]-PART[j].P*r_ji_zero[D])*w/dis/dis;
					}

					PART[i].add_pressure(press_accel);
				}
	//			double coef=-dimension/PART[i].get_density()/PART[i].PND;
				double coef=-dimension/PART[i].get_density()/PART[i].PND0;
				PART[i].mul_pressure(coef);
			}
			else if(PART[i].contact==true)
			{
				//���݈ʒu�ł̎��ӗ��q�����擾
				size_t neighbourN=PART[i].get_current_neighboursID().size();
				double press_accel_temp[3]={0.0, 0.0, 0.0};

				for(int k=0;k<neighbourN;k++)
				{
					int j=PART[i].get_current_neighboursID()[k];

					//�����z�u�ł͋߂��ɑ��݂��Ȃ����q�Ƃ̈��͌��z���v�Z����
					//���ݔz�u�ł̑��΍��W
					double r_ij[3], r_ji[3];
					for(int D=0;D<3;D++)
					{
						r_ij[D]=PART[j].r_temp[D]-PART[i].r_temp[D];
						r_ji[D]=-r_ij[D];
					}

					double dis=0.0; for(int D=0;D<3;D++){dis+=r_ij[D]*r_ij[D];} dis=sqrt(dis); //���ݗ��q�ԋ���

					double w=kernel(re, dis);
	//				if(PART[j].type!=WALL){
					for(int D=0;D<3;D++){
						press_accel_temp[D]+=(PART[i].P*r_ij[D]-PART[j].P*r_ji[D])*w/dis/dis;
					}
	//				}
				}//for(int k=0;k<neighbourN;k++)

				double coef=-dimension/PART[i].get_density()/PART[i].PND;//����͒e���ό`�Ƃ͌���Ȃ��̂�PND�Ŋ���
				for(int D=0;D<3;D++) press_accel_temp[D]*=coef;

				PART[i].add_pressure(press_accel_temp);
			}
//			if(PART[i].r_temp[A_Z]<ground+0.5*le) PART[i].set_stop_on_floor(true);  //�ڐG�m�F�t���O
//			else PART[i].set_stop_on_floor(false);
		}
	}
}

void modify_acceleration(vector<mpselastic> &PART, elastic &ELAST, double **F)
{
	//�����x���[���ɂ��Ă����x�����̓[���ɂȂ�Ȃ��i�����œ��������邽�߁j
	//�����x�[�����͂̒ނ荇��
	//�����ꑬ�x�̔����l���[���ɂȂ�E�E�E�H
	//�Ȃ�Ȃ���΃[���ɂȂ�悤�ȏ�����t��

	double dt=ELAST.get_dt();//OK?�m�F����E�E�E
	double mass=ELAST.get_mass();
	double g[3]={0.0, 0.0, 0.0};
	double total_accel[3];
	double dimension=static_cast<double>(ELAST.get_dimension());
	double re=ELAST.get_r();
	double le=ELAST.get_distancebp();

	g[A_Z]=ELAST.get_g();

	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{
		for(int D=0;D<3;D++) total_accel[D]=(PART[i].get_stress_accel(D)+PART[i].get_stress_visco_accel(D)+PART[i].get_pressure_accel(D)+PART[i].PAcc[D]+g[D]+F[D][i]/mass);

		//�����������x�����������ꍇ������͉ߍS���ł́H
		if((PART[i].get_stop_on_floor()==true) && (total_accel[A_Z]<0.0))
		{
			total_accel[A_Z]+=-2.0*PART[i].u_temp[A_Z]/dt;//���S�e���Փ˂�����
			PART[i].set_total_accel(total_accel);
			PART[i].set_acceleration_upward(false);
		}
		else
		{
			PART[i].set_total_accel(total_accel);
			PART[i].set_acceleration_upward(true);
		}
	}
}

void calc_pre_velocity_and_position(vector<mpselastic> &PART, elastic &ELAST, double **F)
{
	//F[D][i]����͗��q���ɑ΂���x, y, z�������ꂼ��̐���������(�ߓ_�͖@�𗘗p���ׂ�) //F�̒P�ʂ�[N]!!!
	int dimension=ELAST.get_dimension();
	double dt=ELAST.get_dt();
	double mass=ELAST.get_mass();
	double g[3]={0.0, 0.0, 0.0};
	int symplectic=ELAST.get_symp_flag();
	int symplectic_order=ELAST.get_symp_order();
	
	if(dimension==2) g[A_Y]=ELAST.get_g(); 
	if(dimension==3) g[A_Z]=ELAST.get_g();

	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{
		//���̈ʒu�E���̑��x�̍X�V
		for(int D=0;D<3;D++)
		{
			PART[i].r_temp[D]=PART[i].r[D];
			PART[i].u_temp[D]=PART[i].u[D];
		}
		if(PART[i].type==ELASTIC || PART[i].type==MAGELAST|| PART[i].type==TERMINAL1 || PART[i].type==TERMINAL2)
		{
			if(dimension==3)
			{
				//3D�ł̓V���v���N�e�B�b�N�����𖞂����Ȃ��̂ɓK�p���Ă����I�I�C���[�@�Ōv�Z����ƁE�E�E
				if(symplectic==OFF)
				{
					for(int D=0;D<dimension;D++){
//							double u_temp=PART[i].u[D];
//							PART[i].u[D]+=dt*(PART[i].get_normal(D)+PART[i].get_shear(D)+PART[i].get_pressure(D)+PART[i].get_normal_visco(D)+PART[i].get_shear_visco(D)+ELAST.get_P_visco_stress(D, i)+g[D]+F[D][i]/mass);
//							PART[i].u[D]+=dt*(PART[i].get_normal(D)+PART[i].get_shear(D)+PART[i].get_pressure(D)+PART[i].get_normal_visco(D)+PART[i].get_shear_visco(D)+g[D]+F[D][i]/mass);
					}
				}
				else
				{
//					check_velocity_and_position(PART, i, mass, g, F);

					for(int D=0;D<dimension;D++){
				//		PART[i].u_temp[D]+=dt*(PART[i].get_stress_accel(D)+PART[i].get_pressure_accel(D)+PART[i].get_stress_visco_accel(D)+g[D]+F[D][i]/mass);
				//		PART[i].r_temp[D]+=dt*PART[i].u[D];
						PART[i].u_temp[D]=PART[i].u[D];
						PART[i].r_temp[D]=PART[i].r[D];
					}
				}
			}
		}
	}
}

void calc_post_velocity_and_position(vector<mpselastic> &PART, elastic &ELAST,int t)
{
	cout<<"calc_post_velocity_and_position"<<endl;
	//F[D][i]����͗��q���ɑ΂���x, y, z�������ꂼ��̐���������(�ߓ_�͖@�𗘗p���ׂ�) //F�̒P�ʂ�[N]!!!
	int dimension=ELAST.get_dimension();
	double dt=ELAST.get_dt();
	double density=ELAST.get_density();
	double mass=ELAST.get_mass();
	double g[3]={0.0, 0.0, 0.0};
	int symplectic=ELAST.get_symp_flag();
	int symplectic_order=ELAST.get_symp_order();

	double accleration_epsilon=1.0e-12;
	double ground=ELAST.get_ground_position();
	double le=ELAST.get_distancebp();

	if(dimension==2) g[A_Y]=ELAST.get_g(); 
	if(dimension==3) g[A_Z]=ELAST.get_g();

	int tei=0;
	int tei2=120;
	int ue=2904;

if(ELAST.get_model_number()==4){
	if(ELAST.get_poise_flag()==ON){	//���߂̕��͂ł��邾���������������肽��
		for(int i=tei;i<tei2;i++){
			PART[i].u[A_Z]=-0.05;//-0.05
			PART[i+ue].u[A_Z]=0.05;//0.05 
		}
	}
	else {
		for(int i=tei;i<tei2;i++){
			PART[i].u[A_Z]=0.0;
			PART[i+ue].u[A_Z]=0.0;
		}
	}
}
	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{
		if(PART[i].type!=WALL && PART[i].type!=TERMINAL1 && PART[i].type!=TERMINAL2){	//�����Ȃ����q
			if(symplectic==OFF)
			{
				for(int D=0;D<dimension;D++){
	//						double u_temp=PART[i].u[D];
	//						PART[i].u[D]+=dt*(PART[i].get_normal(D)+PART[i].get_shear(D)+PART[i].get_pressure(D)+PART[i].get_normal_visco(D)+PART[i].get_shear_visco(D)+ELAST.get_P_visco_stress(D, i)+g[D]+F[D][i]/mass);
	//						PART[i].u[D]+=dt*(PART[i].get_normal(D)+PART[i].get_shear(D)+PART[i].get_pressure(D)+PART[i].get_normal_visco(D)+PART[i].get_shear_visco(D)+g[D]+F[D][i]/mass);
				}
			}
			else
			{
	//			check_velocity_and_position(PART, i, mass, g, F);

				//���x�̍X�V
				if(ELAST.get_model_number()==4){
				if(!((i>=tei && i<tei2) || (i>=ue && i<ue+tei2))){
				for(int D=0;D<dimension;D++) PART[i].u[D]+=dt*PART[i].get_total_accel(D);
				}
				}
				else {
					for(int D=0;D<dimension;D++) PART[i].u[D]+=dt*PART[i].get_total_accel(D);
				}
				//���x�̏C��
				//�����x���[���ɂ��Ă������œ���������̂ł��̏ꍇ�͑��x���[���ɂ���
				if(PART[i].get_stop_on_floor()==true)//���̋߂��ɂȂ��Ȃ���ɉ������Ȃ�
				{
						//�^���ʕۑ����𖞂����悤��F=��mv/��t�����߂�
						//���S�e���Փ˂����肷��
						//cout<<"PART["<<i<<"].u[A_Z]="<<PART[i].u[A_Z]<<endl;
						//PART[i].u[A_Z]*=0.0;����ł�OK
						PART[i].u[A_Z]*=-1.0;
				}
				
				//�ʒu�̍X�V	
				for(int D=0;D<dimension;D++) {
					PART[i].r[D]+=dt*PART[i].u[D];//���̃X�e�b�v�̈ʒu���i�[
				}
			}
		}
	}
}

void calc_rigid(vector<mpselastic> &PART,vector<Rigidbody> &rigids){
	
	vector<mpselastic> PARTa;
	vector<mpselastic> PARTb;
	int j=0, k=0;
	bool calcf=false;
	//���݂̗��q�ʒu���ړ�
	for(int i=0;i<PART.size();i++){
		if(PART[i].type==TERMINAL1){
			PARTa.push_back(PART[i]);
			calcf=true;
			}		
		else if(PART[i].type==TERMINAL2){
			PARTb.push_back(PART[i]);
			calcf=true;
			}
	}
	if(calcf==true){
		cout<<"���̌v�Z�J�n"<<endl;
		rigids[0].Renew_part_r_v(PARTa);
//		rigids[1].Renew_part_r_v(PARTb);
	//�S�̗��q���X�g�ɍ��̗��q��߂�
	//���q���i�[���鏇�Ԃ͌����Ă������ԂƓ����͂�
		cout<<"���̗��q��S���q���X�g�ɏ㏑��"<<endl;
	for(int i=0;i<PART.size();i++){
		if(PART[i].type==TERMINAL1){
			PART[i].r[A_X]=PARTa[j].r[A_X];
			PART[i].r[A_Y]=PARTa[j].r[A_Y];
			PART[i].r[A_Z]=PARTa[j].r[A_Z];
			j++;
			}		
		else if(PART[i].type==TERMINAL2){
			PART[i].r[A_X]=PARTb[k].r[A_X];
			PART[i].r[A_Y]=PARTb[k].r[A_Y];
			PART[i].r[A_Z]=PARTb[k].r[A_Z];
			k++;
			}
	}
	//���̂̏d�S���ړ�
	cout<<"�d�S�ړ��v�Z�J�n"<<endl;
	rigids[0].Get_rigid_move(PART);
//	rigids[1].Get_rigid_move(PART);

	Micro_AVS avs;
	for(int i=0;i<PART.size();i++){
		avs.make_list(PART[i].r[A_X],PART[i].r[A_Y],PART[i].r[A_Z],0,0,0);
	}
	avs.Output_mgf_MicroAVS("move_particle",1);
	cout<<"���̌v�Z�I��"<<endl;
	}
//	getchar();
	PARTa.clear();
	PARTb.clear();
	
}

//�N�H�[�^�j�I���ɂ��set_r0_ij()�̌v�Z
//�������ӗ��q�����ׂ�Ri�ŉ�]�����ĕێ�����E�E�E�g���܂킵�����̂�
void calc_r0_ij(vector<mpselastic> &PART)
{
	//PART[i].ang[D]�͂��łɌv�Z����Ă���Ƃ����O��Ł�calc_quaternion(), calc_angular_velocity();
	for(int i=0;i<PART.size();i++)
	{
//		if(PART[i].type==ELASTIC)
		{
			size_t neighboursN0=PART[i].get_initial_neighboursID().size();
			double qi[4]={PART[i].ang[0], PART[i].ang[1], PART[i].ang[2], PART[i].ang[3]};

			for(int k=0;k<neighboursN0;k++)
			{
				int j=PART[i].get_initial_neighboursID()[k];
				double qj[4]={PART[j].ang[0], PART[j].ang[1], PART[j].ang[2], PART[j].ang[3]};

				double r_ij_Init[3]; for(int D=0;D<3;D++) r_ij_Init[D]=PART[j].initial_r[D]-PART[i].initial_r[D];
				double r_ji_Init[3]; for(int D=0;D<3;D++) r_ji_Init[D]=-r_ij_Init[D];
				double r0_ij[3], r0_ji[3];

				rotate_r0(r_ij_Init, qi, r0_ij); //qi��p���ĉ�]
				rotate_r0(r_ji_Init, qj, r0_ji); //qj��p���ĉ�]

				PART[i].set_r0_ij(r0_ij);
				PART[i].set_r0_ji(r0_ji);
			}
		}
	}
//	cout<<"��]�����̌v�Z����"<<endl;
}

//�N�H�[�^�j�I����p������]�s��̌v�Z
void rotate_r0(double const *rInit, double const *q, double *result)
{
	result[0]=rInit[0]*(1-2*q[1]*q[1]-2*q[2]*q[2])+rInit[1]*(2*q[0]*q[1]-2*q[3]*q[2])+rInit[2]*(2*q[0]*q[2]+2*q[3]*q[1]);
	result[1]=rInit[0]*(2*q[0]*q[1]+2*q[3]*q[2])+rInit[1]*(1-2*q[0]*q[0]-2*q[2]*q[2])+rInit[2]*(2*q[1]*q[2]-2*q[3]*q[0]);
	result[2]=rInit[0]*(2*q[0]*q[2]-2*q[3]*q[1])+rInit[1]*(2*q[1]*q[2]+2*q[3]*q[0])+rInit[2]*(1-2*q[0]*q[0]-2*q[1]*q[1]);
}

//�p���x�̍X�V�E�E�ELU�����Ƀ|�C���^�łȂ�vector�̎Q�Ƃ�n��
//shared_ptr�ł��ǂ�
void calc_angular_velocity_using_vectorSTL(elastic &ELAST, vector<mpselastic> &PART)
{
	double re=ELAST.get_r();
	
	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{
		//�s����[���ɏ�����
		vector<vector<double>> A(3, vector<double>(3, 0.0));
		vector<double> b(3, 0.0);
		
		size_t neighboursN0=PART[i].get_initial_neighboursID().size(); //���̕����Ӗ����ʂ�₷��

//			if(PART[i].type==static_cast<int>(ELASTIC)) //WALL�̊p���x�͌v�Z���Ȃ�
		{
			for(int k=0;k<neighboursN0;k++)
			{
				int j=PART[i].get_initial_neighboursID()[k];//�Q�Ƃ�����Ă����ق�������

				double X=PART[j].r_temp[A_X]-PART[i].r_temp[A_X];
				double Y=PART[j].r_temp[A_Y]-PART[i].r_temp[A_Y];
				double Z=PART[j].r_temp[A_Z]-PART[i].r_temp[A_Z];
			
				double vX=PART[j].u_temp[A_X]-PART[i].u_temp[A_X];
				double vY=PART[j].u_temp[A_Y]-PART[i].u_temp[A_Y];
				double vZ=PART[j].u_temp[A_Z]-PART[i].u_temp[A_Z];

				double dis0=PART[i].get_initial_distancebps()[k];	//�������q�ԋ���
//				double dis=PART[i].get_current_distancebps()[k]; //NG�I�e�����a���ɂ��闱�q��j�Ƃ͌���Ȃ��I�E�E�E���̊֐��K�v�H
				double dis=sqrt(X*X+Y*Y+Z*Z);//�����OK�Bi��j�Ƃ̊֌W���l���Ă���̂�

				double w=kernel(re, dis0);

				//�W���s��
				A[0][0]+=(Y*Y+Z*Z)*w/dis/dis;
				A[1][1]+=(Z*Z+X*X)*w/dis/dis;
				A[2][2]+=(X*X+Y*Y)*w/dis/dis;			//A[2][2]=(X*X+Y*Y)*w/dis/dis;�{���������E�E�E2012-09-13
				A[0][1]-=(X*Y)*w/dis/dis;
				A[0][2]-=(Z*Y)*w/dis/dis;
				A[1][2]-=(Y*Z)*w/dis/dis;

				//�E�Ӄx�N�g��
				b[0]+=(Y*vZ-Z*vY)*w/dis/dis;
				b[1]+=(Z*vX-X*vZ)*w/dis/dis;
				b[2]+=(X*vY-Y*vX)*w/dis/dis;
//				vector_product(r, v, b);				//���̕�������+=�����v����Ȃ��ib�����̓s�x���Z�b�g�����j
//				for(int D=0;D<3;D++) b[D]*=w/dis/dis;

			}

			//�v�f���Ώ̂Ȃ̂Ő����グ���I����Ă�����
			A[1][0]=A[0][1];
			A[2][0]=A[0][2];
			A[2][1]=A[1][2];

			//LU����@��A��=b�������B����b�Ɋi�[�����
			lu_decomposition(A, b);

			for(int D=0;D<3;D++)
			{
				PART[i].ang_u_temp[D]=b[D];
//				cout<<"PART["<<i<<"].ang_u["<<D<<"]="<<PART[i].ang_u[D]<<" ";
			}
//			cout<<endl;
		}
	}

//	cout<<"�p���x�v�Z����"<<endl;
}

//����v�Z�p�E�E�ELU�����Ƀ|�C���^�łȂ�vector��n��
void calc_quaternion_using_vectorSTL(elastic &ELAST, vector<mpselastic> &PART)
{
	double re=ELAST.get_r();
	double r0[3];

	const int MAX_ITERATION=1000;		//�ő�J��Ԃ���
	const double EPSILON=1.0e-6;//pow(10.0, -12); //�g�������X�B�}�V���C�v�V������DBL_EPSILON�œ��� ��1.0e-16
	int particles_over_max_iteration=0; //MAX_ITERATION�𒴂������q�����J�E���g����

	//�N�H�[�^�j�I���p�̌W���s��p�֐��|�C���^�e�[�u��
	double (* const rotation[])(const double *r, const double *r0, const vector<double> &q)={rotx, roty, rotz};
	double (* const calc_jacobi_matrix[3][4])(const double *r, const double *r0, const vector<double> &q)={
		{rotx_x, rotx_y, rotx_z, rotx_s},
		{roty_x, roty_y, roty_z, roty_s},
		{rotz_x, rotz_y, rotz_z, rotz_s},
	};
	double (* const jacobi_norm[4])(const vector<double> &q)={rotn_x, rotn_y, rotn_z, rotn_s};

	vector<vector<double>> jacobi_matrix(4, vector<double>(4, 0.0));
	vector<double> d(4, 0.0);
	vector<double> qi(4);

	for(int i=0;i<PART.size();i++)
	{
		{
			int itr=0; //�����񐔂̃��Z�b�g
			size_t neighboursN0=PART[i].get_initial_neighboursID().size();//re���Ɋ܂܂�鏉�����q���E�E�EPART[i].N�͏����ł͂Ȃ��I
			
			for(int j=0;j<4;j++){
				for(int k=0;k<4;k++){
					jacobi_matrix[j][k]=0.0;
				}
				d[j]=0.0;
				qi[j]=PART[i].ang[j];
			}

			do{
				for(int k=0;k<neighboursN0;k++)//re���Řa�����
				{
					//j: ID, k: �z��̓Y���I�I ID�Ɨ��q�̓Y���̑ΏƂɂ͋C�����邱�ƁI�I
					int j=PART[i].get_initial_neighboursID()[k];

					if(j!=i)//����if�͂���Ȃ�
					{
						//i����݂����΍��W�E�E�E���݂̉e�����a���̗��q�ł͂Ȃ��A�����ɋߖT�ɂ��������q��ID��p����B�]����PART[i].get_current_position();�͖��Ӗ��I�I�I
						double X=PART[j].r_temp[A_X]-PART[i].r_temp[A_X];
						double Y=PART[j].r_temp[A_Y]-PART[i].r_temp[A_Y];
						double Z=PART[j].r_temp[A_Z]-PART[i].r_temp[A_Z];

						double dis0=PART[i].get_initial_distancebps()[k];//r_init//�����[i][j]�ł͂Ȃ�[i][k]�I�I
						PART[i].get_initial_neighbours_position(k, r0);

						double w=kernel(re, dis0);

						//�E�Ӄx�N�g���̍쐬
						double r[3]={X, Y, Z};//r_ij
						for(int D=0;D<3;D++) d[D]-=rotation[D](r, r0, qi)*w/dis0/dis0;//���ɒ��ӁI(-1)��Y��Ȃ��悤��.�a�����O�ɏ������I

						//���R�r�s��̍쐬(���P���W���̘a��n���悤�ɂ���I)
						for(int row=0;row<3;row++)
							for(int col=0;col<4;col++)
								jacobi_matrix[row][col]+=calc_jacobi_matrix[row][col](r, r0, qi)*w/dis0/dis0;
					}
				}

				//�K�i������
				d[3]=-rotn(qi);
				for(int col=0;col<4;col++) jacobi_matrix[3][col]=jacobi_norm[col](qi); //�K�i�������̕Δ���
				
				lu_decomposition(jacobi_matrix, d);
				for(int D=0;D<4;D++) qi[D]+=d[D]; //qi[]�̍X�V
		
				//�����񐔂̃C���N�������g
				itr++;
			}while(check_q_norm(qi, EPSILON) && itr<MAX_ITERATION);

			if(itr==MAX_ITERATION){
				particles_over_max_iteration++;
			}else{
				for(int D=0;D<4;D++) PART[i].ang[D]=qi[D];
//				std::cout<<"PART["<<i<<"], finished, iteration: "<<itr<<", matrix converged"<<std::endl;
//				cout<<"q convergence: "<<boolalpha<<check_q_residue(qi, EPSILON)<<"; ";
//				for(int m=0;m<4;m++) cout<<"q["<<m<<"]="<<qi[m]<<" ";
//				cout<<"iteration: "<<itr<<endl;
			}
		}
	}
	cout<<"particles_over_max_iteration: "<<particles_over_max_iteration<<endl;
}

//��]�s��ƋK�i�������̕Γ��֐�
inline double rotx(double *r, double *r0, double *q)
{
	return (r[1]*(2*(r0[0]*(q[2]*q[0]-q[1]*q[3])+r0[1]*(q[1]*q[2]+q[0]*q[3]))+r0[2]*(q[2]*q[2]+q[3]*q[3]-q[0]*q[0]-q[1]*q[1]))-r[2]*(2*(r0[0]*(q[0]*q[1]+q[2]*q[3])+r0[2]*(q[1]*q[2]-q[0]*q[3]))+r0[1]*(q[1]*q[1]+q[3]*q[3]-q[2]*q[2]-q[0]*q[0])));
}

inline double roty(double *r, double *r0, double *q)
{
	return (r[2]*(2*(r0[1]*(q[0]*q[1]-q[2]*q[3])+r0[2]*(q[2]*q[0]+q[1]*q[3]))+r0[0]*(q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2]))-r[0]*(2*(r0[0]*(q[2]*q[0]-q[1]*q[3])+r0[1]*(q[1]*q[2]+q[0]*q[3]))+r0[2]*(q[2]*q[2]+q[3]*q[3]-q[0]*q[0]-q[1]*q[1])));
}

inline double rotz(double *r, double *r0, double *q)
{
	return (r[0]*(2*(r0[0]*(q[0]*q[1]+q[2]*q[3])+r0[2]*(q[1]*q[2]-q[0]*q[3]))+r0[1]*(q[1]*q[1]+q[3]*q[3]-q[2]*q[2]-q[0]*q[0]))-r[1]*(2*(r0[1]*(q[0]*q[1]-q[2]*q[3])+r0[2]*(q[2]*q[0]+q[1]*q[3]))+r0[0]*(q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2])));
}

inline double rotn(double *q)
{
	double norm=0.0;
	for(int i=0;i<4;i++) norm+=q[i]*q[i];
	return norm-1.0;
}

//��]�s��ƋK�i�������̕Γ��֐�(�I�[�o�[���[�h)
inline double rotx(const double *r, const double *r0, const vector<double> &q)
{
	return (r[1]*(2*(r0[0]*(q[2]*q[0]-q[1]*q[3])+r0[1]*(q[1]*q[2]+q[0]*q[3]))+r0[2]*(q[2]*q[2]+q[3]*q[3]-q[0]*q[0]-q[1]*q[1]))-r[2]*(2*(r0[0]*(q[0]*q[1]+q[2]*q[3])+r0[2]*(q[1]*q[2]-q[0]*q[3]))+r0[1]*(q[1]*q[1]+q[3]*q[3]-q[2]*q[2]-q[0]*q[0])));
}

inline double roty(const double *r, const double *r0, const vector<double> &q)
{
	return (r[2]*(2*(r0[1]*(q[0]*q[1]-q[2]*q[3])+r0[2]*(q[2]*q[0]+q[1]*q[3]))+r0[0]*(q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2]))-r[0]*(2*(r0[0]*(q[2]*q[0]-q[1]*q[3])+r0[1]*(q[1]*q[2]+q[0]*q[3]))+r0[2]*(q[2]*q[2]+q[3]*q[3]-q[0]*q[0]-q[1]*q[1])));
}

inline double rotz(const double *r, const double *r0, const vector<double> &q)
{
	return (r[0]*(2*(r0[0]*(q[0]*q[1]+q[2]*q[3])+r0[2]*(q[1]*q[2]-q[0]*q[3]))+r0[1]*(q[1]*q[1]+q[3]*q[3]-q[2]*q[2]-q[0]*q[0]))-r[1]*(2*(r0[1]*(q[0]*q[1]-q[2]*q[3])+r0[2]*(q[2]*q[0]+q[1]*q[3]))+r0[0]*(q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2])));
}

inline double rotn(const vector<double> &q)
{
	double norm=0.0;
	for(int i=0;i<4;i++) norm+=q[i]*q[i];
	return norm-1.0;
}

//�N�H�[�^�j�I���̎c�����`�F�b�N����
bool check_q_norm(const double *q, const double EPSILON)
{
	double norm=0.0;
	for(unsigned i=0;i<4;i++) norm+=q[i]*q[i];
	return (fabs(1.0-norm)>EPSILON) ? true: false;
}

bool check_q_norm(const vector<double> &q, const double EPSILON)
{
	double norm=0.0;
	for(unsigned i=0;i<4;i++) norm+=q[i]*q[i];
	return (fabs(1.0-norm)>EPSILON) ? true: false;
}

//�G�l���M�[�̊m�F
void calc_hamiltonian(elastic &ELAST, int t)
{
//	double hamiltonian=ELAST.get_elastic_energy()+ELAST.get_kinetic_energy()+ELAST.get_potential_energy();
	ofstream fout1("./Elastic/hamiltonian.dat", ios::app);
	ofstream fout2("./Elastic/kinetic_energy.dat", ios::app);
	ofstream fout3("./Elastic/elastic_energy.dat", ios::app);
	ofstream fout4("./Elastic/elastic_energy1.dat", ios::app);
	ofstream fout5("./Elastic/elastic_energy2.dat", ios::app);
	ofstream fout6("./Elastic/potential_energy.dat", ios::app);
	
	if(fout1.fail() || fout2.fail() || fout3.fail() || fout4.fail() || fout5.fail() || fout6.fail()){
		system("mkdir Elastic");
		ofstream fout1("./Elastic/hamiltonian.dat", ios::app);
		ofstream fout2("./Elastic/kinetic_energy.dat", ios::app);
		ofstream fout3("./Elastic/elastic_energy.dat", ios::app);
		ofstream fout4("./Elastic/elastic_energy1.dat", ios::app);
		ofstream fout5("./Elastic/elastic_energy2.dat", ios::app);
		ofstream fout6("./Elastic/potential_energy.dat", ios::app);

		if(fout1.fail() || fout2.fail() || fout3.fail() || fout4.fail())//�Ď��s
		{
			cout<<"�G�l���M�[�m�F�p�̃f�B���N�g�����J���܂���"<<endl;
			cout<<"�v���O�������I�����܂�"<<endl;
			exit(1);
		}
	}

	fout1<<t<<"\t"<<ELAST.get_hamiltonian()<<endl;
	fout2<<t<<"\t"<<ELAST.get_kinetic_energy()<<endl;
	fout3<<t<<"\t"<<ELAST.get_elastic_energy()<<endl;
	fout4<<t<<"\t"<<ELAST.get_elastic_energy1()<<endl;
	fout5<<t<<"\t"<<ELAST.get_elastic_energy2()<<endl;
	fout6<<t<<"\t"<<ELAST.get_potential_energy()<<endl;

	fout1.close();
	fout2.close();
	fout3.close();
	fout4.close();
	fout5.close();
	fout6.close();

	//ofstream fout5("difference_of_EE.dat", ios::app);
	//fout5<<t<<"\t"<<(ELAST.get_elastic_energy()-ELAST.get_last_elastic_energy())<<endl;
	//fout5.close()

/*
	if(ELAST.get_FEM_flag()==ON)
	{
		//�{����ELASTIC�ȊO�ɂ��G�l���M�[���v�Z���ׂ�
		//<0���Ƃ�����Ƃ����덷�ŃX�C�b�`������̂�<
//		if((ELAST.get_elastic_energy()-ELAST.get_last_elastic_energy())<0)
		if((ELAST.get_elastic_energy()-ELAST.get_last_elastic_energy())<-1e-10)
		{
			ELAST.set_FEM_switch(ON);
			cout<<"step: "<<t<<"FEM switch: "<<boolalpha<<ELAST.get_FEM_switch()<<endl;
		}else{
			cout<<"FEM switch: "<<boolalpha<<ELAST.get_FEM_switch()<<endl;
		}
	}
*/
	ELAST.set_last_elastic_energy(ELAST.get_elastic_energy());
	ELAST.set_last_kinetic_energy(ELAST.get_kinetic_energy());
	ELAST.set_last_potential(ELAST.get_potential_energy());

}

//�x�N�g���ς̌v�Z
void vector_product(double *a, double *b, double *result)
{
	result[0]=a[A_Y]*b[A_Z]-a[A_Z]*b[A_Y];
	result[1]=a[A_Z]*b[A_X]-a[A_X]*b[A_Z];
	result[2]=a[A_X]*b[A_Y]-a[A_Y]*b[A_X];
}

//�����x�̊m�F
void check_velocity_and_position(vector<mpselastic> &PART, elastic &ELAST, const int i, double **F)
{
	double mass=ELAST.get_mass();
	double g[3]={0.0, 0.0, 0.0};
	if(2==ELAST.get_dimension()) g[1]=ELAST.get_g();
	else g[2]=g[1]=ELAST.get_g();

	cout<<"checking 3d velocity and position"<<endl;
	for(int D=0;D<3;D++) cout<<"PART["<<i<<"].u["<<D<<"]="<<PART[i].u[D]<<" ";
	cout<<endl;
	for(int D=0;D<3;D++) cout<<"PART["<<i<<"].r["<<D<<"]="<<PART[i].r[D]<<" ";
	cout<<endl;
	for(int D=0;D<3;D++) cout<<"PART["<<i<<"].pr["<<D<<"]="<<PART[i].get_pressure_accel(D)<<" ";
	cout<<endl;
	for(int D=0;D<3;D++) cout<<"g["<<D<<"]="<<g[D]<<" ";
	cout<<endl;
	for(int D=0;D<3;D++) cout<<"F["<<D<<"]["<<i<<"]"<<F[D][i]<<" ";
	cout<<endl;
	for(int D=0;D<3;D++) cout<<"F["<<D<<"]["<<i<<"]/mass="<<F[D][i]/mass<<" ";
	cout<<endl;
}

//���x�̏C��
void calc_modified_density(vector<mpselastic> &PART, elastic &ELAST)
{
	//�̐ςƗ��q����p�������̂�freeon()�Ŏ��s
	for(int i=0;i<PART.size();i++)
	{
//		cout<<"initial density: "<<PART[i].get_density()<<endl;
		double ratio=PART[i].PND/PART[i].PND0;
		PART[i].set_density(ELAST.get_density()*ratio);
//		cout<<"current density: "<<PART[i].get_density()<<endl;
	}
}

void rigid_calc(vector<mpselastic> &PART, elastic &ELAST)
{
	for(int i=0;i<PART.size();i++)
	{
		if(PART[i].type==TERMINAL1 || PART[i].type==TERMINAL2) //���̗��q
		{

		}
	}
}

//���l�����ɂ�郄�R�r�s��
//���S�����ߎ��B�œK�ȍ��ݕ��̑��݂ɒ��ӁI�I�I
void numerical_jacobian(double **J, double *r, double *r0, double *q)
{
	//�֐��|�C���^�e�[�u��
	double (* const rotation[])(double *r, double *r0, double *q)={rotx, roty, rotz};

	const double h=1.0e-3;//���ݕ��B�A�_�v�e�B�u�ɕύX���ׂ�
	double result0=0.0;

	//�[���ɂȂ�ꍇ�͏���
	result0=((rotation[0](r, r0, q))-(rotation[0](r, r0, q)))/(2*h);

}

void calc_nonlinear_elastic(vector<mpselastic> &PART, elastic &ELAST, int t, double **F)
{
	vector<vector<double>> residual_acceleration(3, vector<double>(PART.size(), 0.0)); //�c�������x

	//������
	for(int D=0;D<ELAST.get_dimension();D++)
	{
		for(int i=0;i<PART.size();i++)
		{
			residual_acceleration[D][i]=F[D][i];
		}
	}

	while(true)
	{	
	//	calc_quaternion(ELAST, PART); //�N�H�[�^�j�I���̌v�Z//OK
		calc_quaternion_using_vectorSTL(ELAST, PART);
	//	calc_angular_velocity(ELAST, PART); //�p���x�̌v�Z//OK
		calc_angular_velocity_using_vectorSTL(ELAST, PART);
		calc_r0_ij(PART); //�����ʒu�x�N�g���̉�]
		calc_nonlinear_accel_for_3D(PART, ELAST);
		calc_pressure_and_contact(PART, ELAST);	//���͂ƐڐG���͂ɂ������x���v�Z
		PART[2].check_acceleration(); //���q�̉����x���`�F�b�N

		calc_residual_acceleration(PART, ELAST, residual_acceleration, t, F);

		calc_nonlinear_velocity_and_position(PART, ELAST, residual_acceleration);
		calc_hamiltonian(ELAST, t);
	}
}

void calc_nonlinear_accel_for_3D(vector<mpselastic> &PART, elastic &ELAST)
{
	double dimension=static_cast<double>(ELAST.get_dimension());
	double mass=ELAST.get_mass();
	double density=ELAST.get_density();
	double re=ELAST.get_r(); //����͊���re*le
	double vis=ELAST.get_nensei(); //�S�x�i���S�x�ł͂Ȃ��j
	double g=ELAST.get_g();

	double KE=0.0;	//�^���G�l���M�[
	double EE1=0.0;	//�Ђ��݃G�l���M�[
	double EE2=0.0;	//�̐ςЂ��݂ɂ��
	double PE=0.0;	//�|�e���V����
	double ground=ELAST.get_ground_position(); //���̂����W�̎擾

	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{   
		double lambda_average=0.0; //lambda�̕��ϒl
		double particle_density=0.0; //���q�����x�E�E�E�ڐG����Ŏg�� 

		//WALL�͒e���ό`���Ȃ��̂ŏ��O
//		if(PART[i].type==ELASTIC || PART[i].type==INELASTIC || PART[i].type==BOELASTIC)
		{
			double EE1_temp=0.0;
			double EE2_temp=0.0;
			double volumetric_strain=0.0; //�̐ςЂ���

			//�^���G�l���M�[�̍X�V
			for(int D=0;D<3;D++) KE+=PART[i].u[D]*PART[i].u[D];

			//�ʒu�|�e���V�����̍X�V
			PE+=(PART[i].r[2]-ground);

			int neighboursN0=static_cast<int>(PART[i].get_initial_neighboursID().size());

			for(int k=0;k<neighboursN0;k++) //pressure�ƍ��킹�Ȃ��Ƌ����H
			{

				int j=PART[i].get_initial_neighboursID()[k];

				//�����z�u��NONELASTIC�Ȃ��̂��߂��ɂ���ꍇ�͏����𕪂���
				//���ꂪ�Ȃ��ƃG�l���M�[���x�����������Ȃ�͂��E�E�E
				if(PART[j].type==ELASTIC)
				{				
					//i����݂������z�u�ł̑��΍��W
					double r_ij_Init[3], r_ji_Init[3]; 
					for(int D=0;D<3;D++){
						r_ij_Init[D]=PART[j].initial_r[D]-PART[i].initial_r[D];
						r_ji_Init[D]=-r_ij_Init[D];
					}

					//���ݔz�u�ł̑��΍��W
					double r_ij[3], r_ji[3];
					for(int D=0;D<3;D++)
					{
						r_ij[D]=PART[j].r[D]-PART[i].r[D];
						r_ji[D]=-r_ij[D];
					}

					double dis=0.0; for(int D=0;D<3;D++){dis+=r_ij[D]*r_ij[D];} dis=sqrt(dis); //���ݗ��q�ԋ���
					double dis0=0.0; for(int D=0;D<3;D++){dis0+=r_ij_Init[D]*r_ij_Init[D];} dis0=sqrt(dis0); //�������q�ԋ���
			
					//��]�s��Ri�ɂ��r_ij_Init�̉�]�E�E�E�ȍ~��r_ij_zero���g��
					//qi�͎��X���X�ς��̂Ŏg���񂵂ł��Ȃ��E�E�EPART�̃p�����[�^�Ƃ��ė^����̂��ǂ�
					double r_ij_zero[3], r_ji_zero[3];

					for(int D=0;D<3;D++) r_ij_zero[D]=PART[i].get_r0_ij()[k].get_comp()[D];
					for(int D=0;D<3;D++) r_ji_zero[D]=PART[i].get_r0_ji()[k].get_comp()[D];

				//�d�݊֐�
					double w=kernel(re, dis0);

				//���q�����x�̍X�V
					particle_density+=w;

					//��]��̏����z�u���Έʒu�̒P�ʃx�N�g��
					double n_ij[3], n_ji[3];
					for(int D=0;D<3;D++) 
					{
						n_ij[D]=r_ij_zero[D]/dis0;
						n_ji[D]=r_ji_zero[D]/dis0;
					}

					//�ψʃx�N�g��
					double U_ij[3], U_ji[3]; 
					for(int D=0;D<3;D++)
					{
						U_ij[D]=r_ij[D]-r_ij_zero[D];
						U_ji[D]=r_ji[D]-r_ji_zero[D];
					}

				//�Ђ��݃x�N�g��
					double E_ij[3], E_ji[3];
					for(int D=0;D<3;D++)
					{
						E_ij[D]=U_ij[D]/dis0;
						E_ji[D]=U_ji[D]/dis0;
					}

				//�����O���ƃ|�A�\����̌v�Z
					double En_ij[3], En_ji[3];//�c�Ђ��݃x�N�g��
					double naiseki=0.0;
					for(int D=0;D<3;D++) naiseki+=E_ij[D]*r_ij[D];
					for(int D=0;D<3;D++) En_ij[D]=naiseki*r_ij[D]/dis/dis;	//�Ђ��݃x�N�g����r_ij�����֐��ˉe
					double ns=sqrt(En_ij[0]*En_ij[0]+En_ij[1]*En_ij[1]+En_ij[2]*En_ij[2]); //�c�Ђ���

					//�V���R�[���̏ꍇ
//					youngs_modulus=-6.11619E-11*x^9+1.66522E-8*x^8-1.93253E-6*x^7+1.2459E-4*x^6-4.87141E-3*x^5+0.118275*x^4-1.75252*x^3+14.9237*x^2-64.3273*x+123.411;
					double youngs_modulus=(((((((((-6.11619E-11*ns+1.66522E-8)*ns-1.93253E-6)*ns+1.2459E-4)*ns-4.87141E-3)*ns+0.118275)*ns-1.75252)*ns+14.9237)*ns-64.3273)*ns+123.411)*100;
//					-9.6723E-12*x^10+1.22758E-9*x^9-6.73685E-8*x^8+2.09375E-7*x^7-4.05552E-5*x^6+5.08152E-4*x^5-0.00414923*x^4+0.0218461*x^3-0.072922*x^2+0.151892*x^1+0.17713
					double poisson_ratio=(((((((((-9.6723E-12*ns+1.22758E-9)*ns-6.73685E-8)*ns+2.09375E-7)*ns-4.05552E-5)*ns+5.08152E-4)*ns-0.00414923)*ns+0.0218461)*ns-0.072922)*ns+0.151892)*ns+0.17713;

//					cout<<"youngs_modulus: "<<youngs_modulus<<", poisson_ratio: "<<poisson_ratio<<endl;

					//�����萔mu
					double shear_modulus=youngs_modulus/(2.0*(1.0+poisson_ratio));
					//�����萔lambda
					double lambda=(poisson_ratio*youngs_modulus)/((1.0+poisson_ratio)*(1.0-2.0*poisson_ratio));

//					cout<<"shear_modulus: "<<shear_modulus<<", lambda: "<<lambda<<endl;

				//���͌v�Z�E�E�Ee_vol_i�̃�[]�i�ψʂ̔��U�j���v�Z
					for(int D=0;D<3;D++) volumetric_strain+=E_ij[D]*n_ij[D]*w; //�̐ςЂ��݌v�Z�̏���
					PART[i].P+=volumetric_strain*lambda;
					EE2_temp+=PART[i].P*volumetric_strain;
					//PART[i].P*=w;�E�E�E���̂����ꂾ�Ƃ��܂������Ȃ��E�E�Ew���ݏ�ő������̂œ�����O((()*w+smt)*w+smt)*w...

					lambda_average+=lambda;

				//�Ђ��݃G�l���M�[
					for(int D=0;D<3;D++) EE1_temp+=E_ij[D]*E_ij[D]*w*shear_modulus;

				//���̓x�N�g��
					double sigma_ij[3], sigma_ji[3];
					for(int D=0;D<3;D++)
					{
						sigma_ij[D]=2*shear_modulus*w*E_ij[D]/PART[i].PND/dis0;
						sigma_ji[D]=2*shear_modulus*w*E_ji[D]/PART[j].PND/dis0;

						sigma_ij[D]-=sigma_ji[D];
					}
					PART[i].add_stress_accel(sigma_ij);

	//				if(fabs((dis-dis0)/dis0)>1.04) w=0;//�j�����

					//�Ђ��ݑ��x�̌v�Z
					double strain_vi[3], strain_vj[3];

					//i����݂��Ђ��ݑ��x
					strain_vi[0]=((PART[j].u[0]-PART[i].u[0])-(PART[i].ang_u[1]*r_ij[2]-PART[i].ang_u[2]*r_ij[1]))/dis; //X������̂Ђ��ݑ��x ��������ˉe����K�v����H�H
					strain_vi[1]=((PART[j].u[1]-PART[i].u[1])-(PART[i].ang_u[2]*r_ij[0]-PART[i].ang_u[0]*r_ij[2]))/dis; //Y������̂Ђ��ݑ��x
					strain_vi[2]=((PART[j].u[2]-PART[i].u[2])-(PART[i].ang_u[0]*r_ij[1]-PART[i].ang_u[1]*r_ij[0]))/dis; //Z������̂Ђ��ݑ��x
					//anglar_u1..3�͔z��ɏ�������

					//j����݂��Ђ��ݑ��x
					strain_vj[0]=((PART[i].u[0]-PART[j].u[0])-(PART[j].ang_u[1]*r_ji[2]-PART[j].ang_u[2]*r_ji[1]))/dis; //X������̂Ђ��ݑ��x ��������ˉe����K�v����H�H
					strain_vj[1]=((PART[i].u[1]-PART[j].u[1])-(PART[j].ang_u[2]*r_ji[0]-PART[j].ang_u[0]*r_ji[2]))/dis; //Y������̂Ђ��ݑ��x
					strain_vj[2]=((PART[i].u[2]-PART[j].u[2])-(PART[j].ang_u[0]*r_ji[1]-PART[j].ang_u[1]*r_ji[0]))/dis; //Z������̂Ђ��ݑ��x
			
					double sigma_v_ij[3], sigma_v_ji[3];
					for(int D=0;D<3;D++)
					{
						sigma_v_ij[D]=2*vis*w*strain_vi[D]/PART[i].PND/dis;
						sigma_v_ji[D]=2*vis*w*strain_vj[D]/PART[j].PND/dis;

						sigma_v_ij[D]-=sigma_v_ji[D];
					}

					PART[i].add_stress_visco_accel(sigma_v_ij);
				}
			}//for(int k=0;k<neighboursN0;k++)���[�v�I��

			PART[i].P*=-dimension/PART[i].PND;//�̐ςЂ��݂ɂ�鈳�͂����߂�ꂽ�i����̓G�l���M�[�v�Z����O�ɋ��߂Ă������Ɓj
			volumetric_strain*=-dimension/PART[i].PND;//�̐ςЂ��݂����߂�ꂽ�i����̓G�l���M�[�v�Z����O�ɋ��߂Ă������Ɓj
			PART[i].volumetric_strain=volumetric_strain;
			//�ڐG����i���ǂ��܂ށj�E�E�E�����Ō��������Ȃ��ƃG���[�H
			//�����x�v�Z��ELASTIC�����l����΂悢�̂�if(.type==ELASTIC)�̃��[�v�̒��ɓ���Ă��ς��Ȃ�
//			double n00=PART[i].PND0;
//			if(PART[i].P<(PART[i].PND-n00)/n00) PART[i].P=(PART[i].PND-n00)/n00; //���q�����x�̑�����������Βu������

			double coef=dimension/PART[i].get_density();

			PART[i].mul_stress_accel(coef);
			PART[i].mul_stress_visco_accel(coef);

			//�Ђ��݃G�l���M�[
//			EE1+=(EE1_temp)*(shear_modulus*mass*dimension/PART[i].PND/density);	//��񍀂̃�
//			EE2+=0.5*lambda*mass*(PART[i].P*PART[i].P)/density;//��O���̃�
			EE1+=EE1_temp/PART[i].PND/PART[i].get_density();//��񍀂̃� ���x�͂��ꂼ��Ⴄ
			EE2+=EE2_temp/PART[i].get_density();//��O���̃� ���x�͂��ꂼ��Ⴄ

			/****�G�l���M�[�v�Z�ɂ�lambda_average���g���ׂ��E�E�E�H****/

//			PART[i].P*=lambda;//�̐ςЂ��݂ɂ�鈳�͂����߂�ꂽ
		}//if(PART[i].type==ELASTIC || PART[i].type==INELASTIC ||PART[i].type==BOELASTIC)�I��

		//���ꎩ�͓̂����Ȃ���WALL�������x�i���́j��L���Ă���
		//���q���̑������l�����Ȃ���΂Ȃ�Ȃ��̂ł����if�̊O�ɒu��
		//�E�E�EPART[i].P�͂��ׂĂ̗��q���l������K�v�����邪�A�G�l���M�[�v�Z�ɓ����ׂ��ł͂Ȃ�
		double n00=PART[i].PND0;
//2012-11-26 ���͂Ɨ��q�����x���r���Ă����E�E�E�I
//		if(PART[i].P<(PART[i].PND-n00)/n00) PART[i].P=lambda*(PART[i].PND-n00)/n00; //���q�����x�̑�����������Βu������

		lambda_average/=static_cast<int>(PART[i].get_initial_neighboursID().size());

		//���q�����x�̑�����������Βu������
		//����P�͒e���̂��ǂ��󂯂Ă��鈳�͂ł��邱�Ƃɒ���

		if(PART[i].P<lambda_average*((PART[i].PND-n00)/n00)) PART[i].P=lambda_average*(PART[i].PND-n00)/n00; //�E�E�E���l�I�ɉ����N���Ă��邩���ׂ�I�I

	}//for(i=0;i<PART.size();i++)

	KE*=0.5*mass;
	PE*=-g*mass;
	EE1*=mass*dimension;
	EE2*=0.5*mass;
	ELAST.set_kinetic(KE);
	ELAST.set_elastic_energy(EE1+EE2);//�S�n�̒e���G�l���M�[
	ELAST.set_potential(PE);
	ELAST.set_hamiltonian(KE+EE1+EE2+PE);

}

void calc_residual_acceleration(vector<mpselastic> &PART, elastic &ELAST, vector<vector<double>> &res_accel, int t, double **F)
{
	double mass=ELAST.get_mass();
	double g[3]={0.0, 0.0, 0.0};
	g[A_Z]=ELAST.get_g();

	for(int i=0;i<PART.size();i++)
	{
		if(PART[i].type==ELASTIC || PART[i].type==MAGELAST)
		{
			for(int D=0;D<DIMENSION;D++) res_accel[D][i]=PART[i].get_stress_accel(D)+PART[i].get_pressure_accel(D)+PART[i].get_stress_visco_accel(D)+g[D]+F[D][i]/mass;
		}
	}
}

void calc_nonlinear_velocity_and_position(vector<mpselastic> &PART, elastic &ELAST, vector<vector<double>> &residual_acceleration)
{
	//F[D][i]����͗��q���ɑ΂���x, y, z�������ꂼ��̐���������(�ߓ_�͖@�𗘗p���ׂ�) //F�̒P�ʂ�[N]!!!
	double dt=ELAST.get_dt();

	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{
		if(PART[i].type==ELASTIC || PART[i].type==MAGELAST)
		{
//			check_velocity_and_position(PART, i, mass, g, F);

			for(int D=0;D<DIMENSION;D++){
				PART[i].u[D]+=dt*residual_acceleration[D][i];
				PART[i].r[D]+=dt*PART[i].u[D];
			}
		}
	}
}

/*********�e�X�g�̈�***********/
/*********�e�X�g�̈�***********/
/*********�e�X�g�̈�***********/
/*********�e�X�g�̈�***********/
/*********�e�X�g�̈�***********/
/*********�e�X�g�̈�***********/
/*********�e�X�g�̈�***********/
/*********�e�X�g�̈�***********/
/*********�e�X�g�̈�***********/


//�p���x�̍X�V�E�E�E�X�P�[�����O������i�Ƃ�������v�Z�Œ��ډ����Ă��ǂ��̂ł́H�j
void calc_angular_velocity(elastic &ELAST, vector<mpselastic> &PART)
{
	double re=ELAST.get_r();
	bool pivot_check=ELAST.get_pivot_check();

	//�Ώ̌W���s��
	double **A=new double*[3]; for(int D=0;D<3;D++) A[D]=new double[3];
	//�E�Ӄx�N�g��
	double b[3];

	#pragma omp parallel
	{
		for(int i=0;i<PART.size();i++)
		{
			//�s��̏�����
			for(int m=0;m<3;m++)
			{
				for(int n=0;n<3;n++) A[m][n]=0.0;
				b[m]=0.0;
			}

			size_t neighboursN0=PART[i].get_initial_neighboursID().size(); //���̕����Ӗ����ʂ�₷��

			if(PART[i].type==(int)ELASTIC) //WALL�̊p���x�͌v�Z���Ȃ�
			{
				for(int k=0;k<neighboursN0;k++)
				{
					int j=PART[i].get_initial_neighboursID()[k];//�Q�Ƃ�����Ă����ق�������

					double X=PART[j].r[A_X]-PART[i].r[A_X];
					double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
					double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
			
					double vX=PART[j].u[A_X]-PART[i].u[A_X];
					double vY=PART[j].u[A_Y]-PART[i].u[A_Y];
					double vZ=PART[j].u[A_Z]-PART[i].u[A_Z];

	//				double r[3]={X, Y, Z};
	//				double v[3]={vX, vY, vZ};

					double dis0=PART[i].get_initial_distancebps()[k];	//�������q�ԋ���
	//				double dis=PART[i].get_current_distancebps()[k]; //NG�I�e�����a���ɂ��闱�q��j�Ƃ͌���Ȃ��I�E�E�E���̊֐��K�v�H
					double dis=sqrt(X*X+Y*Y+Z*Z);//�����OK�Bi��j�Ƃ̊֌W���l���Ă���̂�

					double w=kernel(re, dis0);

					//�W���s��
					A[0][0]+=(Y*Y+Z*Z)*w/dis/dis;
					A[1][1]+=(Z*Z+X*X)*w/dis/dis;
					A[2][2]+=(X*X+Y*Y)*w/dis/dis;			//A[2][2]=(X*X+Y*Y)*w/dis/dis;�{���������E�E�E2012-09-13
					A[0][1]-=(X*Y)*w/dis/dis;
					A[0][2]-=(Z*Y)*w/dis/dis;
					A[1][2]-=(Y*Z)*w/dis/dis;
	//				A[0][1]=A[1][0]=-(X*Y)*w/dis/dis;		//�����=-�ł͂Ȃ��E�E�E
	//				A[0][2]=A[2][0]=-(Z*Y)*w/dis/dis;
	//				A[1][2]=A[2][1]=-(Y*Z)*w/dis/dis;

					//�E�Ӄx�N�g��
					b[0]+=(Y*vZ-Z*vY)*w/dis/dis;
					b[1]+=(Z*vX-X*vZ)*w/dis/dis;
					b[2]+=(X*vY-Y*vX)*w/dis/dis;
	//				vector_product(r, v, b);				//���̕�������+=�����v����Ȃ��ib�����̓s�x���Z�b�g�����j
	//				for(int D=0;D<3;D++) b[D]*=w/dis/dis;

				}

				//�v�f���Ώ̂Ȃ̂Ő����グ���I����Ă�����
				A[1][0]=A[0][1];
				A[2][0]=A[0][2];
				A[2][1]=A[1][2];

				//LU����@��A��=b�������B����b�Ɋi�[�����
				lu_decomposition(A, b, 3, pivot_check);
	//			cout<<"angular velocity: "<<endl;
				for(int D=0;D<3;D++)
				{
					PART[i].ang_u[D]=b[D];
	//				cout<<"PART["<<i<<"].ang_u["<<D<<"]="<<PART[i].ang_u[D]<<" ";
				}
	//			cout<<endl;
			}
		}
	}

	for(int D=0;D<3;D++) delete [] A[D];
	delete [] A;
}

//�N�H�[�^�j�I���̌v�Z�E�E�E�K�i�������̓����Ă͂����Ȃ�
//|q[i]|=<1�Ȃ̂ŁA�l�����������Ȃ�����N�H�[�^�j�I���̒�`����l���ċ����I�Ƀ��Z�b�g���邱�Ƃ��L���E�E�E�H
void calc_quaternion(elastic &ELAST, vector<mpselastic> &PART)
{
	double re=ELAST.get_r();
	double r0[3];
	bool pivot_check=ELAST.get_pivot_check();

	const int MAX_ITERATION=1000;		//�ő�J��Ԃ���
	const double EPSILON=1.0e-6;//pow(10.0, -12);	//�g�������X�B�}�V���C�v�V������DBL_EPSILON�œ��� ��1.0e-16
	int particles_over_max_iteration=0;

	//�N�H�[�^�j�I���p�̌W���s��p�֐��|�C���^�e�[�u��
	double (* const rotation[])(double *r, double *r0, double *q)={rotx, roty, rotz};
	double (* const calc_jacobi_matrix[3][4])(double *r, double *r0, double *q)={
		{rotx_x, rotx_y, rotx_z, rotx_s},
		{roty_x, roty_y, roty_z, roty_s},
		{rotz_x, rotz_y, rotz_z, rotz_s},
	};
	double (* const jacobi_norm[4])(double *q)={rotn_x, rotn_y, rotn_z, rotn_s};
	double d[4];//�E�Ӄx�N�g��
	double** jacobi_matrix=new double*[4];//���R�r�s��

	for(int D=0;D<4;D++) jacobi_matrix[D]=new double[4];

	for(int i=0;i<PART.size();i++)
	{
		if(PART[i].type==ELASTIC || PART[i].type==MAGELAST) //WALL�̊p���x�͌v�Z���Ȃ�
		{
			int itr=0; //�����񐔂̃��Z�b�g
			size_t neighboursN0=PART[i].get_initial_neighboursID().size();//re���Ɋ܂܂�鏉�����q���E�E�EPART[i].N�͏����ł͂Ȃ��I
			double qi[4]={PART[i].ang[0], PART[i].ang[1], PART[i].ang[2], PART[i].ang[3]};

//			for(int m=0;m<4;m++) cout<<"itr: "<<itr<<", q["<<m<<"]="<<qi[m]<<endl;
//			cout<<"PART["<<i<<"].N="<<PART[i].N<<", neighboursN0="<<neighboursN0<<endl;

//			cout<<"q�����l: "<<endl;
			for(int m=0;m<4;m++){
				d[m]=0.0;
				for(int n=0;n<4;n++) jacobi_matrix[m][n]=0.0;//jacobi�s��̏�����
//				cout<<"q["<<m<<"]="<<qi[m]<<" ";
			}
//			cout<<endl;

			do{
//				for(int j=0;j<neighbours0;j++)//re���Řa�����B�����ID�Əƍ����Ă��Ȃ������P���ȃG���[�I�I�I
				for(int k=0;k<neighboursN0;k++)//re���Řa�����
				{
					//j: ID, k: �z��̓Y���I�I ID�Ɨ��q�̓Y���̑ΏƂɂ͋C�����邱�ƁI�I
					int j=PART[i].get_initial_neighboursID()[k];

					if(j!=i)//����if�͂���Ȃ�
					{
						//i����݂����΍��W�E�E�E���݂̉e�����a���̗��q�ł͂Ȃ��A�����ɋߖT�ɂ��������q��ID��p����B�]����PART[i].get_current_position();�͖��Ӗ��I�I�I
						double X=PART[j].r_temp[A_X]-PART[i].r_temp[A_X];
						double Y=PART[j].r_temp[A_Y]-PART[i].r_temp[A_Y];
						double Z=PART[j].r_temp[A_Z]-PART[i].r_temp[A_Z];

						double dis0=PART[i].get_initial_distancebps()[k];//r_init//�����[i][j]�ł͂Ȃ�[i][k]�I�I
						PART[i].get_initial_neighbours_position(k, r0);

						double w=kernel(re, dis0);

//						cout<<"dis0="<<dis0<<", r0[0]="<<r0[0]<<", r0[1]="<<r0[1]<<", r0[2]="<<r0[2]<<endl;

						//�E�Ӄx�N�g���̍쐬
						double r[3]={X, Y, Z};//r_ij
						for(int D=0;D<3;D++) d[D]-=rotation[D](r, r0, qi)*w/dis0/dis0;//���ɒ��ӁI(-1)��Y��Ȃ��悤��.�a�����O�ɏ������I

						//���R�r�s��̍쐬(���P���W���̘a��n���悤�ɂ���I)
						for(int row=0;row<3;row++)
							for(int col=0;col<4;col++)
								jacobi_matrix[row][col]+=calc_jacobi_matrix[row][col](r, r0, qi)*w/dis0/dis0;
					}
				}

				//�K�i������
				d[3]=-rotn(qi);
				for(int col=0;col<4;col++) jacobi_matrix[3][col]=jacobi_norm[col](qi); //�K�i�������̕Δ���
				
				lu_decomposition(jacobi_matrix, d, 4, pivot_check);
//				pivot_gauss(jacobi_matrix, d, 4, pivot_check); //���R�r�A���͕K��J(1, 1)>0�Ȃ̂�pivot�I���Ȃ��̕�������������
				for(int D=0;D<4;D++) qi[D]+=d[D]; //qi[]�̍X�V
		
				//�����񐔂̃C���N�������g

				itr++;
			}while(check_q_norm(qi, EPSILON) && itr<MAX_ITERATION);
//			}while(vector_norm2(d, 0, 4)>EPSILON && itr<MAX_ITERATION);//�����͖{���Ȃ�c����p����ׂ�

			if(itr==MAX_ITERATION){
				particles_over_max_iteration++;
			}else{
				for(int D=0;D<4;D++) PART[i].ang[D]=qi[D];
//				std::cout<<"PART["<<i<<"], finished, iteration: "<<itr<<", matrix converged"<<std::endl;
//				cout<<"q convergence: "<<boolalpha<<check_q_residue(qi, EPSILON)<<"; ";
//				for(int m=0;m<4;m++) cout<<"q["<<m<<"]="<<qi[m]<<" ";
//				cout<<"iteration: "<<itr<<endl;
			}
		}
	}

	cout<<"particles_over_max_iteration: "<<particles_over_max_iteration<<endl;
	for(int D=0;D<4;D++) delete [] jacobi_matrix[D];
	delete [] jacobi_matrix;
}

void calc_velocity_and_position(vector<mpselastic> &PART, elastic &ELAST, double **F)
{
	//F[D][i]����͗��q���ɑ΂���x, y, z�������ꂼ��̐���������(�ߓ_�͖@�𗘗p���ׂ�) //F�̒P�ʂ�[N]!!!
	int dimension=ELAST.get_dimension();
	double dt=ELAST.get_dt();
	double density=ELAST.get_density();
	double mass=ELAST.get_mass();
	double g[3]={0.0, 0.0, 0.0};
	int symplectic=ELAST.get_symp_flag();
	int symplectic_order=ELAST.get_symp_order();
	
	if(dimension==2) g[A_Y]=ELAST.get_g(); 
	if(dimension==3) g[A_Z]=ELAST.get_g();

//	cout<<"mass= "<<mass<<endl;
//	cout<<"density= "<<density<<endl;
//	cout<<"PART[2].density= "<<PART[2].get_density()<<endl;

	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{
		if(PART[i].type==ELASTIC || PART[i].type==MAGELAST)
		{
			if(dimension==3)
			{
				//3D�ł̓V���v���N�e�B�b�N�����𖞂����Ȃ��̂ɓK�p���Ă����I�I�C���[�@�Ōv�Z����ƁE�E�E
				if(symplectic==OFF)
				{
					for(int D=0;D<dimension;D++){
//							double u_temp=PART[i].u[D];
//							PART[i].u[D]+=dt*(PART[i].get_normal(D)+PART[i].get_shear(D)+PART[i].get_pressure(D)+PART[i].get_normal_visco(D)+PART[i].get_shear_visco(D)+ELAST.get_P_visco_stress(D, i)+g[D]+F[D][i]/mass);
//							PART[i].u[D]+=dt*(PART[i].get_normal(D)+PART[i].get_shear(D)+PART[i].get_pressure(D)+PART[i].get_normal_visco(D)+PART[i].get_shear_visco(D)+g[D]+F[D][i]/mass);
					}
				}
				else
				{
//					check_velocity_and_position(PART, i, mass, g, F);

					for(int D=0;D<dimension;D++){
						PART[i].u[D]+=dt*(PART[i].get_stress_accel(D)+PART[i].get_pressure_accel(D)+PART[i].get_stress_visco_accel(D)+g[D]+F[D][i]/mass);
						PART[i].r[D]+=dt*PART[i].u[D];
					}
				}
			}
		}
	}
}

void calc_contact(vector<mpselastic> &PART, elastic &ELAST)
{
	double dimension=static_cast<double>(ELAST.get_dimension());
	double le=ELAST.get_le();
	double re=ELAST.get_r();
	double density=ELAST.get_density();

	vector<int>::iterator ip;
	map<int, double>::iterator mp; //map�T���p��iterator

	//main��reload_INDEX2()���s���Ă����Ȃ���PART[i]�̏�񂪂߂��Ⴍ����ɂȂ�

	//�c���E���k�����x�ƐڐG���͂̌v�Z
	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{   
//		countOK=countNOT=0;
		size_t neighbourN=PART[i].get_current_neighboursID().size();//���݈ʒu�ł̎��ӗ��q�����擾
		double qi[4]={PART[i].ang[0], PART[i].ang[1], PART[i].ang[2], PART[i].ang[3]};//���݈ʒu�ł�i�̃N�H�[�^�j�I�����擾

		//����re���ɑ��݂��Ȃ����w=0�Ȃ̂Ō���0�Bfind�ŒT���Ηǂ�
		for(int k=0;k<neighbourN;k++)
		{
			int j=PART[i].get_current_neighboursID()[k];//���ݎ��ӂɂ��闱�q��ID���擾

			//���ݔz�u�ł̑��΍��W
			double r_ij[3], r_ji[3];
			for(int D=0;D<3;D++)
			{
				r_ij[D]=PART[j].r[D]-PART[i].r[D];
				r_ji[D]=-r_ij[D];
			}

			double dis=0.0; for(int D=0;D<3;D++){dis+=r_ij[D]*r_ij[D];} dis=sqrt(dis); //���ݗ��q�ԋ���

			//�����z�u��re���ɂ��������q��ID��T��
			ip=find(PART[i].get_initial_neighboursID().begin(), PART[i].get_initial_neighboursID().end(), j);
			
			//���qj�������z�u��re���ɂ���ꍇ
			if(ip!=PART[i].get_initial_neighboursID().end())//���������ꍇ�i���ݎ��ӂɂ͏����z�u��ID������j
			{
				//i����݂������z�u�ł̑��΍��W
				//i��ELASTIC�������łȂ����ŏꍇ����
				if(PART[i].type==ELASTIC || PART[i].type==MAGELAST)
				{
					if(PART[j].type==ELASTIC)//�ŏ�����ގ��ŕ�����ƒe���̂ǂ����̏Փ˂ɑΉ��ł��Ȃ��Ȃ�
					{
						//�e���̂̓��͂ł͏d�ݕt����dis0���g���ׂ��i�����z�u����̕ό`���d�v�Ȃ̂Łj
						double r_ij_Init[3], r_ji_Init[3]; 
						for(int D=0;D<3;D++)
						{
							r_ij_Init[D]=PART[j].initial_r[D]-PART[i].initial_r[D];
							r_ji_Init[D]=-r_ij_Init[D];
						}

						//��]�s��Ri�ɂ��r_ij_Init�̉�]�E�E�E�ȍ~��r_ij_zero���g��
						//qi�͎��X���X�ς��̂Ŏg���񂵂ł��Ȃ��E�E�EPART�̃p�����[�^�Ƃ��ė^����̂��ǂ�
						double r_ij_zero[3], r_ji_zero[3];
						double qj[4]={PART[j].ang[0], PART[j].ang[1], PART[j].ang[2], PART[j].ang[3]};

						rotate_r0(r_ij_Init, qi, r_ij_zero);
						rotate_r0(r_ji_Init, qj, r_ji_zero);

						double dis0=0.0; for(int D=0;D<3;D++){dis0+=r_ij_Init[D]*r_ij_Init[D];} dis0=sqrt(dis0); //�������q�ԋ���

						double w=kernel(re, dis0);

						double press_accel[3];
						for(int D=0;D<3;D++)
							//�d�݊֐��̉��d���ς̂Ƃ����ύX�Bi�̗��q�����x�Ōv�Z����
							press_accel[D]=(PART[i].P*r_ij_zero[D]-PART[j].P*r_ji_zero[D])*w/dis0/dis0;
							//�d�݂͏����̋������g�������z�͌��ݔz�u�̕����g��
	//						press_accel[D]=((PART[i].P*r_ij[D]/PART[i].PND)-(PART[j].P*r_ji[D]/PART[j].PND))*w/dis/dis;
						//dis0^2�ŏ����Ă���Ƃ������Ƃ́C���ς��l���Ă���\��������D

						PART[i].add_pressure(press_accel);

					}else if(PART[j].type==WALL){
					
						//�ό`���Ȃ����̂ɑ΂��Ă͏d�ݕt����dis���g���ׂ��i�����z�u����̕ό`�͐�������R�͂̌����������d�v�Ȃ̂Łj
						//���E�����̓X�e�b�v���Ƃɕς��
						//�ڐG����ꍇ�������߂Â��̂�

						double w=kernel(re, dis);

						double press_accel[3];
						for(int D=0;D<3;D++)
							press_accel[D]=(PART[i].P*r_ij[D]-PART[j].P*r_ji[D])*w/dis/dis;

						PART[i].add_pressure(press_accel);//���͂̉����x��������
					}				
				}
				else if(PART[i].type==WALL)
				{
					if(PART[j].type==ELASTIC)//�ŏ�����ގ��ŕ�����ƒe���̂ǂ����̏Փ˂ɑΉ��ł��Ȃ��Ȃ�
					{
						double w=kernel(re, dis);

						double press_accel[3];
						for(int D=0;D<3;D++)
							//�d�݊֐��̉��d���ς̂Ƃ����ύX�Bi�̗��q�����x�Ōv�Z����
							press_accel[D]=(PART[i].P*r_ij[D]-PART[j].P*r_ji[D])*w/PART[i].PND/dis/dis;
							//�d�݂͏����̋������g�������z�͌��ݔz�u�̕����g��
	//						press_accel[D]=((PART[i].P*r_ij[D]/PART[i].PND)-(PART[j].P*r_ji[D]/PART[j].PND))*w/dis/dis;
						//dis0^2�ŏ����Ă���Ƃ������Ƃ́C���ς��l���Ă���\��������D

						PART[i].add_pressure(press_accel);

					}		
				}
			}
			else//�����z�u��re���ɗ��qj���Ȃ��ꍇ�i���ӂɗ��q���ڋ߂��Ă���or����Ă������j
			{
				//���qi��ELASTIC��WALL���ŕ�����
				if(PART[i].type==ELASTIC || PART[i].type==MAGELAST)
				{
					if(PART[j].type==ELASTIC)//�ŏ�����ގ��ŕ�����ƒe���̂ǂ����̏Փ˂ɑΉ��ł��Ȃ��Ȃ�
					{
						double r_ij_Init[3], r_ji_Init[3]; 
						for(int D=0;D<3;D++)
						{
							r_ij_Init[D]=PART[j].initial_r[D]-PART[i].initial_r[D];
							r_ji_Init[D]=-r_ij_Init[D];
						}

						double dis0=0.0; for(int D=0;D<3;D++){dis0+=r_ij_Init[D]*r_ij_Init[D];} dis0=sqrt(dis0); //�������q�ԋ���

						//��]�s��Ri�ɂ��r_ij_Init�̉�]�E�E�E�ȍ~��r_ij_zero���g��
						//qi�͎��X���X�ς��̂Ŏg���񂵂ł��Ȃ��E�E�EPART�̃p�����[�^�Ƃ��ė^����̂��ǂ�
						double r_ij_zero[3], r_ji_zero[3];
						double qj[4]={PART[j].ang[0], PART[j].ang[1], PART[j].ang[2], PART[j].ang[3]};

						rotate_r0(r_ij_Init, qi, r_ij_zero);
						rotate_r0(r_ji_Init, qj, r_ji_zero);
			
						double w=kernel(re, dis0);

						double press_accel[3];
						for(int D=0;D<3;D++)
							//���͌��z�������x�Ƃ��Čv�Z���Ă��邱�Ƃɑ�������i���q�����̒P�ʃx�N�g�������Ɍ��z�x�N�g����������j
	//						press_accel[D]=((PART[i].P*r_ij_zero[D]/PART[i].PND)-(PART[j].P*r_ji_zero[D]/PART[j].PND))*w/dis0/dis0;
							press_accel[D]=(PART[i].P*r_ij_zero[D]-PART[j].P*r_ji_zero[D])*w/dis0/dis0;

						PART[i].add_pressure(press_accel);

					}else if(PART[j].type==WALL){//��Ɠ��l

						double w=kernel(re, dis);

						double press_accel[3];
						for(int D=0;D<3;D++)
							//���͖͂c���Ƌt�����ɂ�����͂��E�E�E
							press_accel[D]=(PART[i].P*r_ij[D]-PART[j].P*r_ji[D])*w/dis/dis;

						PART[i].add_pressure(press_accel);
					}
				}
				else if(PART[i].type==WALL)
				{
					if(PART[j].type==ELASTIC)//�ŏ�����ގ��ŕ�����ƒe���̂ǂ����̏Փ˂ɑΉ��ł��Ȃ��Ȃ�
					{
						double w=kernel(re, dis);

						double press_accel[3];
						for(int D=0;D<3;D++)
							//���͌��z�������x�Ƃ��Čv�Z���Ă��邱�Ƃɑ�������i���q�����̒P�ʃx�N�g�������Ɍ��z�x�N�g����������j
	//						press_accel[D]=((PART[i].P*r_ij_zero[D]/PART[i].PND)-(PART[j].P*r_ji_zero[D]/PART[j].PND))*w/dis0/dis0;
							press_accel[D]=(PART[i].P*r_ij[D]-PART[j].P*r_ji[D])*w/dis/dis;

						PART[i].add_pressure(press_accel);

					}
				}
			}

		}//for(int k=0;k<neighbourN;k++)�E�E�Ek���[�v�I��

		double coef=-dimension/PART[i].get_density()/PART[i].PND;
		PART[i].mul_pressure(coef);
//		ELAST.set_P_visco_stress(D, i, (ELAST.get_P_visco_stress(D, i)*dimension*(-1)/ePND[i]/density));//n0�łȂ�n[i]���g���Ă���E�E�E
	}
}


void calc_accel_for_3D_ver_3(vector<mpselastic> &PART, elastic &ELAST)
{
//	cout<<"3D analysis start"<<endl;
	double dimension=static_cast<double>(ELAST.get_dimension());
	double mag_shear_modulus=ELAST.get_mag_shear_modulus();
	double elas_shear_modulus=ELAST.get_elas_shear_modulus();
	double mag_lambda=ELAST.get_mag_lambda();
	double elas_lambda=ELAST.get_elas_lambda();
	double mass=ELAST.get_mass();
	double density=ELAST.get_density();
	double re=ELAST.get_r(); //����͊���re*le
	double vis=ELAST.get_nensei(); //�S�x�i���S�x�ł͂Ȃ��j
	double g=ELAST.get_g();

	double KE=0.0;	//�^���G�l���M�[
	double EE1=0.0;	//�Ђ��݃G�l���M�[
	double EE2=0.0;	//�̐ςЂ��݂ɂ��
	double PE=0.0;	//�|�e���V����
	double ground=ELAST.get_ground_position(); //���̂����W�̎擾

	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{   
		//WALL�͒e���ό`���Ȃ��̂ŏ��O
		if(PART[i].type==ELASTIC || PART[i].type==MAGELAST)
		{
			double EE1_temp=0.0;

			//�^���G�l���M�[�̍X�V
			for(int D=0;D<3;D++) KE+=PART[i].u[D]*PART[i].u[D];

			//�ʒu�|�e���V�����̍X�V
			PE+=(PART[i].r[2]-ground);

			int neighboursN0=static_cast<int>(PART[i].get_initial_neighboursID().size());

			for(int k=0;k<neighboursN0;k++) //pressure�ƍ��킹�Ȃ��Ƌ����H
			{
				int j=PART[i].get_initial_neighboursID()[k];

				//�����z�u��NONELASTIC�Ȃ��̂��߂��ɂ���ꍇ�͏����𕪂���
				//���ӗ��q��WALL�ł��ǂ��̂ł́E�E�E����if����Ȃ��E�E�E(2012-11-29)
//				if(PART[j].type==ELASTIC)
				{				
					//i����݂������z�u�ł̑��΍��W
					double r_ij_Init[3], r_ji_Init[3]; 
					for(int D=0;D<3;D++){
						r_ij_Init[D]=PART[j].initial_r[D]-PART[i].initial_r[D];
						r_ji_Init[D]=-r_ij_Init[D];
					}

					//���ݔz�u�ł̑��΍��W
					double r_ij[3], r_ji[3];
					for(int D=0;D<3;D++)
					{
						r_ij[D]=PART[j].r[D]-PART[i].r[D];
						r_ji[D]=-r_ij[D];
					}

					double dis=0.0; for(int D=0;D<3;D++){dis+=r_ij[D]*r_ij[D];} dis=sqrt(dis); //���ݗ��q�ԋ���
					double dis0=0.0; for(int D=0;D<3;D++){dis0+=r_ij_Init[D]*r_ij_Init[D];} dis0=sqrt(dis0); //�������q�ԋ���
			
					double r_ij_zero[3], r_ji_zero[3];

					//calc_r0_ij(PART)��"���݈ʒu�̃N�H�[�^�j�I�����g����"��]�����������ʒu�x�N�g�����擾
					for(int D=0;D<3;D++) r_ij_zero[D]=PART[i].get_r0_ij()[k].get_comp()[D];
					for(int D=0;D<3;D++) r_ji_zero[D]=PART[i].get_r0_ji()[k].get_comp()[D];

					double w=kernel(re, dis0);

					//��]��̏����z�u���Έʒu�̒P�ʃx�N�g��
					double n_ij[3], n_ji[3];
					for(int D=0;D<3;D++) 
					{
						n_ij[D]=r_ij_zero[D]/dis0;
						n_ji[D]=r_ji_zero[D]/dis0;
					}

					//�ψʃx�N�g��
					double U_ij[3], U_ji[3]; 
					for(int D=0;D<3;D++)
					{
						U_ij[D]=r_ij[D]-r_ij_zero[D];
						U_ji[D]=r_ji[D]-r_ji_zero[D];
					}

					//�Ђ��݃x�N�g��
					double E_ij[3], E_ji[3];
					for(int D=0;D<3;D++)
					{
						E_ij[D]=U_ij[D]/dis0;
						E_ji[D]=U_ji[D]/dis0;
					}
				
					//���͌v�Z�̏����E�E�Ee_vol_i�̃�[]�i�ψʂ̔��U�j���v�Z
					//�����Ђ��݂̘a�����i���̓e���\���̑Ίp�����̘a�ɏd�ݕt���j
					for(int D=0;D<3;D++) PART[i].P+=E_ij[D]*n_ij[D]*w; 

					//�Ђ��݃G�l���M�[
					for(int D=0;D<3;D++) EE1_temp+=E_ij[D]*E_ij[D]*w;

					//���̓x�N�g��
					double sigma_ij[3], sigma_ji[3];
					for(int D=0;D<3;D++)
					{
						if(PART[i].type==MAGELAST){
						sigma_ij[D]=2*mag_shear_modulus*w*E_ij[D]/PART[i].PND/dis0;
						sigma_ji[D]=2*mag_shear_modulus*w*E_ji[D]/PART[j].PND/dis0;
						}
						else if(PART[i].type==ELASTIC){
						sigma_ij[D]=2*elas_shear_modulus*w*E_ij[D]/PART[i].PND/dis0;
						sigma_ji[D]=2*elas_shear_modulus*w*E_ji[D]/PART[j].PND/dis0;
						}

						sigma_ij[D]-=sigma_ji[D];
					}
					PART[i].add_stress_accel(sigma_ij);

	//				if(fabs((dis-dis0)/dis0)>1.04) w=0;//�j�����

					//�Ђ��ݑ��x�̌v�Z
					double strain_vi[3], strain_vj[3];

					//i����݂��Ђ��ݑ��x
					strain_vi[0]=((PART[j].u[0]-PART[i].u[0])-(PART[i].ang_u[1]*r_ij[2]-PART[i].ang_u[2]*r_ij[1]))/dis; //X������̂Ђ��ݑ��x ��������ˉe����K�v����H�H
					strain_vi[1]=((PART[j].u[1]-PART[i].u[1])-(PART[i].ang_u[2]*r_ij[0]-PART[i].ang_u[0]*r_ij[2]))/dis; //Y������̂Ђ��ݑ��x
					strain_vi[2]=((PART[j].u[2]-PART[i].u[2])-(PART[i].ang_u[0]*r_ij[1]-PART[i].ang_u[1]*r_ij[0]))/dis; //Z������̂Ђ��ݑ��x
					//anglar_u1..3�͔z��ɏ�������

					//j����݂��Ђ��ݑ��x
					strain_vj[0]=((PART[i].u[0]-PART[j].u[0])-(PART[j].ang_u[1]*r_ji[2]-PART[j].ang_u[2]*r_ji[1]))/dis; //X������̂Ђ��ݑ��x ��������ˉe����K�v����H�H
					strain_vj[1]=((PART[i].u[1]-PART[j].u[1])-(PART[j].ang_u[2]*r_ji[0]-PART[j].ang_u[0]*r_ji[2]))/dis; //Y������̂Ђ��ݑ��x
					strain_vj[2]=((PART[i].u[2]-PART[j].u[2])-(PART[j].ang_u[0]*r_ji[1]-PART[j].ang_u[1]*r_ji[0]))/dis; //Z������̂Ђ��ݑ��x
			
					double sigma_v_ij[3], sigma_v_ji[3];
					for(int D=0;D<3;D++)
					{
						sigma_v_ij[D]=2*vis*w*strain_vi[D]/PART[i].PND/dis;
						sigma_v_ji[D]=2*vis*w*strain_vj[D]/PART[j].PND/dis;

						sigma_v_ij[D]-=sigma_v_ji[D];
					}

					PART[i].add_stress_visco_accel(sigma_v_ij);
				}//if(PART[j].type==ELASTIC)�I���E�E�E������WALL����̂���f
			}//for(int k=0;k<neighboursN0;k++)���[�v�I��

			PART[i].P*=-dimension/PART[i].PND;//�̐ςЂ��݂����߂�ꂽ�i����̓G�l���M�[�v�Z����O�ɋ��߂Ă������Ɓj

			//�ڐG����i���ǂ��܂ށj�E�E�E�����Ō��������Ȃ��ƃG���[�H
			//�����x�v�Z��ELASTIC�����l����΂悢�̂�if(.type==ELASTIC)�̃��[�v�̒��ɓ���Ă��ς��Ȃ�
			double n00=PART[i].PND0;
//			if(PART[i].P<(PART[i].PND-n00)/n00) PART[i].P=(PART[i].PND-n00)/n00; //���q�����x�̑�����������Βu������

			double coef=dimension/PART[i].get_density();

			PART[i].mul_stress_accel(coef);
			PART[i].mul_stress_visco_accel(coef);

			//�Ђ��݃G�l���M�[
//			EE1+=(EE1_temp)*(shear_modulus*mass*dimension/PART[i].PND/density);	//��񍀂̃�
//			EE2+=0.5*lambda*mass*(PART[i].P*PART[i].P)/density;//��O���̃�
			if(PART[i].type==MAGELAST){
			EE1+=(EE1_temp/PART[i].PND/PART[i].get_density())*mag_shear_modulus*mass*dimension;//��񍀂̃� ���x�͂��ꂼ��Ⴄ
			EE2+=((PART[i].P*PART[i].P)/PART[i].get_density())*0.5*mass*mag_lambda;//��O���̃� ���x�͂��ꂼ��Ⴄ

			PART[i].P*=mag_lambda;//�̐ςЂ��݂ɂ�鈳�͂����߂�ꂽ
			}
			else if(PART[i].type==ELASTIC){
			EE1+=(EE1_temp/PART[i].PND/PART[i].get_density())*elas_shear_modulus*mass*dimension;//��񍀂̃� ���x�͂��ꂼ��Ⴄ
			EE2+=((PART[i].P*PART[i].P)/PART[i].get_density())*0.5*mass*elas_lambda;//��O���̃� ���x�͂��ꂼ��Ⴄ

			PART[i].P*=elas_lambda;//�̐ςЂ��݂ɂ�鈳�͂����߂�ꂽ
			}
		}//if(PART[i].type==ELASTIC)�I��

		//���ꎩ�͓̂����Ȃ���WALL�������x�i���́j��L���Ă���
		//�E�E�EPART[i].P�͂��ׂĂ̗��q���l������K�v�����邪�A�G�l���M�[�v�Z�ɓ����ׂ��ł͂Ȃ�
		double n00=PART[i].PND0;
//		if(PART[i].P<(PART[i].PND-n00)/n00) PART[i].P=lambda*(PART[i].PND-n00)/n00; //���q�����x�̑�����������Βu������
		//2012-11-26 ���͂Ɨ��q�����x���r���Ă����E�E�E�I
		if(PART[i].type==MAGELAST){
			if(PART[i].P<mag_lambda*((PART[i].PND-n00)/n00)) PART[i].P=mag_lambda*(PART[i].PND-n00)/n00; //���q�����x�̑�����������Βu������
		}
		else if(PART[i].type==ELASTIC){
			if(PART[i].P<elas_lambda*((PART[i].PND-n00)/n00)) PART[i].P=elas_lambda*(PART[i].PND-n00)/n00;
		}

	}//for(i=0;i<PART.size();i++)���[�v�I��

	KE*=0.5*mass;
	PE*=-g*mass;
//	EE1*=shear_modulus*mass*dimension;
//	EE2*=0.5*mass*lambda;
	ELAST.set_kinetic(KE);
	ELAST.set_elastic_energy(EE1+EE2);//�S�n�̒e���G�l���M�[
	ELAST.set_potential(PE);
	ELAST.set_hamiltonian(KE+EE1+EE2+PE);
}

void calc_accel_for_3D_ver_2(vector<mpselastic> &PART, elastic &ELAST)
{
//	cout<<"3D analysis start"<<endl;
	double dimension=static_cast<double>(ELAST.get_dimension());
	double mag_shear_modulus=ELAST.get_mag_shear_modulus();
	double elas_shear_modulus=ELAST.get_elas_shear_modulus();
	double mag_lambda=ELAST.get_mag_lambda();
	double elas_lambda=ELAST.get_elas_lambda();
	double mass=ELAST.get_mass();
	double density=ELAST.get_density();
	double re=ELAST.get_r(); //����͊���re*le
	double vis=ELAST.get_nensei(); //�S�x�i���S�x�ł͂Ȃ��j
	double g=ELAST.get_g();

	double KE=0.0;	//�^���G�l���M�[
	double EE1=0.0;	//�Ђ��݃G�l���M�[
	double EE2=0.0;	//�̐ςЂ��݂ɂ��
	double PE=0.0;	//�|�e���V����
	double ground=ELAST.get_ground_position(); //���̂����W�̎擾

	for(int i=0;i<PART.size();i++)
	{   
		//WALL�͒e���ό`���Ȃ��̂ŏ��O
		if(PART[i].type==ELASTIC || PART[i].type==MAGELAST)
		{
			double EE1_temp=0.0;

			//�^���G�l���M�[�̍X�V
			for(int D=0;D<3;D++) KE+=PART[i].u[D]*PART[i].u[D];

			//�ʒu�|�e���V�����̍X�V
			PE+=(PART[i].r[2]-ground);

			double qi[4]={PART[i].ang[0], PART[i].ang[1], PART[i].ang[2], PART[i].ang[3]}; //i�̃N�H�[�^�j�I��
			size_t neighboursN0=PART[i].get_initial_neighboursID().size();

			for(int k=0;k<neighboursN0;k++) //pressure�ƍ��킹�Ȃ��Ƌ����H
			{
				int j=PART[i].get_initial_neighboursID()[k];

				double qj[4]={PART[j].ang[0], PART[j].ang[1], PART[j].ang[2], PART[j].ang[3]}; //j�̃N�H�[�^�j�I��

//				PART[i].get_initial_neighbours_position(k, X0, Y0, Z0);�@���_�������Ƃ���
//				double dis0=PART[i].get_initial_distancebps()[k];//�������q�ԋ���

				//i����݂������z�u�ł̑��΍��W
				double r_ij_Init[3], r_ji_Init[3]; 
				for(int D=0;D<3;D++){
					r_ij_Init[D]=PART[j].initial_r[D]-PART[i].initial_r[D];
					r_ji_Init[D]=-r_ij_Init[D];
				}

				//���ݔz�u�ł̑��΍��W
				double r_ij[3], r_ji[3];
				for(int D=0;D<3;D++)
				{
					r_ij[D]=PART[j].r[D]-PART[i].r[D];
					r_ji[D]=-r_ij[D];
				}

				double dis=0.0; for(int D=0;D<3;D++){dis+=r_ij[D]*r_ij[D];} dis=sqrt(dis); //���ݗ��q�ԋ���
				double dis0=0.0; for(int D=0;D<3;D++){dis0+=r_ij_Init[D]*r_ij_Init[D];} dis0=sqrt(dis0); //�������q�ԋ���
			
				//��]�s��Ri�ɂ��r_ij_Init�̉�]�E�E�E�ȍ~��r_ij_zero���g��
				//qi�͎��X���X�ς��̂Ŏg���񂵂ł��Ȃ��E�E�EPART�̃p�����[�^�Ƃ��ė^����̂��ǂ�
				double r_ij_zero[3], r_ji_zero[3];
				rotate_r0(r_ij_Init, qi, r_ij_zero);
				rotate_r0(r_ji_Init, qj, r_ji_zero);

				double w=kernel(re, dis0);

				//��]��̏����z�u���Έʒu�̒P�ʃx�N�g��
				double n_ij[3], n_ji[3];
				for(int D=0;D<3;D++) 
				{
					n_ij[D]=r_ij_zero[D]/dis0;
					n_ji[D]=r_ji_zero[D]/dis0;
				}

				//�ψʃx�N�g��
				double U_ij[3], U_ji[3]; 
				for(int D=0;D<3;D++)
				{
					U_ij[D]=r_ij[D]-r_ij_zero[D];
					U_ji[D]=r_ji[D]-r_ji_zero[D];
				}

				//�Ђ��݃x�N�g��
				double E_ij[3], E_ji[3];
				for(int D=0;D<3;D++)
				{
					E_ij[D]=U_ij[D]/dis0;
					E_ji[D]=U_ji[D]/dis0;
				}
				
				//���͌v�Z�̏����E�E�Ee_vol_i�̃�[]�i�ψʂ̔��U�j���v�Z
				for(int D=0;D<3;D++) PART[i].P+=E_ij[D]*n_ij[D]*w; 
//				PART[i].P*=w;

				//�Ђ��݃G�l���M�[
				for(int D=0;D<3;D++) EE1_temp+=E_ij[D]*E_ij[D]*w;
//				EE1_temp*=w;

				//���̓x�N�g��
				double sigma_ij[3], sigma_ji[3];
				for(int D=0;D<3;D++)
				{
					if(PART[i].type==MAGELAST){
					sigma_ij[D]=2*mag_shear_modulus*w*E_ij[D]/PART[i].PND/dis0;
					sigma_ji[D]=2*mag_shear_modulus*w*E_ji[D]/PART[i].PND/dis0;
					}
					else if(PART[i].type==ELASTIC){
					sigma_ji[D]=2*elas_shear_modulus*w*E_ji[D]/PART[j].PND/dis0;
					sigma_ij[D]=2*elas_shear_modulus*w*E_ij[D]/PART[j].PND/dis0;
					}

					sigma_ij[D]-=sigma_ji[D];
				}
				PART[i].add_stress_accel(sigma_ij);

//				if(fabs((dis-dis0)/dis0)>1.04) w=0;//�j�����

				//�Ђ��ݑ��x�̌v�Z
				double strain_vi[3], strain_vj[3];

				//i����݂��Ђ��ݑ��x
				strain_vi[0]=((PART[j].u[0]-PART[i].u[0])-(PART[i].ang_u[1]*r_ij[2]-PART[i].ang_u[2]*r_ij[1]))/dis; //X������̂Ђ��ݑ��x ��������ˉe����K�v����H�H
				strain_vi[1]=((PART[j].u[1]-PART[i].u[1])-(PART[i].ang_u[2]*r_ij[0]-PART[i].ang_u[0]*r_ij[2]))/dis; //Y������̂Ђ��ݑ��x
				strain_vi[2]=((PART[j].u[2]-PART[i].u[2])-(PART[i].ang_u[0]*r_ij[1]-PART[i].ang_u[1]*r_ij[0]))/dis; //Z������̂Ђ��ݑ��x
				//anglar_u1..3�͔z��ɏ�������

				//j����݂��Ђ��ݑ��x
				strain_vj[0]=((PART[i].u[0]-PART[j].u[0])-(PART[j].ang_u[1]*r_ji[2]-PART[j].ang_u[2]*r_ji[1]))/dis; //X������̂Ђ��ݑ��x ��������ˉe����K�v����H�H
				strain_vj[1]=((PART[i].u[1]-PART[j].u[1])-(PART[j].ang_u[2]*r_ji[0]-PART[j].ang_u[0]*r_ji[2]))/dis; //Y������̂Ђ��ݑ��x
				strain_vj[2]=((PART[i].u[2]-PART[j].u[2])-(PART[j].ang_u[0]*r_ji[1]-PART[j].ang_u[1]*r_ji[0]))/dis; //Z������̂Ђ��ݑ��x
			
				double sigma_v_ij[3], sigma_v_ji[3];
				for(int D=0;D<3;D++)
				{
					sigma_v_ij[D]=2*vis*w*strain_vi[D]/PART[i].PND/dis;
					sigma_v_ji[D]=2*vis*w*strain_vj[D]/PART[j].PND/dis;

					sigma_v_ij[D]-=sigma_v_ji[D];
				}

				PART[i].add_stress_visco_accel(sigma_v_ij);
			}//for(int k=0;k<neighboursN0;k++)���[�v�I��

			PART[i].P*=-dimension/PART[i].PND;//�̐ςЂ��݂����߂�ꂽ�i����̓G�l���M�[�v�Z����O�ɋ��߂Ă������Ɓj

			//�ڐG����i���ǂ��܂ށj�E�E�E�����Ō��������Ȃ��ƃG���[�H
			//�����x�v�Z��ELASTIC�����l����΂悢�̂�if(.type==ELASTIC)�̃��[�v�̒��ɓ���Ă��ς��Ȃ�
			double n00=PART[i].PND0;
			if(PART[i].P<(PART[i].PND-n00)/n00) PART[i].P=(PART[i].PND-n00)/n00; //���q�����x�̑�����������Βu������

			//ELASTIC�łȂ����̂��v�Z���Ă���
			double coef=dimension/PART[i].get_density();

			PART[i].mul_stress_accel(coef);
			PART[i].mul_stress_visco_accel(coef);

			//�Ђ��݃G�l���M�[
//			EE1+=(EE1_temp)*(shear_modulus*mass*dimension/PART[i].PND/density);	//��񍀂̃�
//			EE2+=0.5*lambda*mass*(PART[i].P*PART[i].P)/density;//��O���̃�
			
			
			if(PART[i].type==MAGELAST){
				EE1+=(EE1_temp/PART[i].PND/PART[i].get_density())*mag_shear_modulus*mass*dimension;//��񍀂̃� ���x�͂��ꂼ��Ⴄ
			EE2+=((PART[i].P*PART[i].P)/PART[i].get_density())*0.5*mass*mag_lambda;//��O���̃� ���x�͂��ꂼ��Ⴄ
			PART[i].P*=mag_lambda;//�̐ςЂ��݂ɂ�鈳�͂����߂�ꂽ
			}
			else if(PART[i].type==ELASTIC){
				EE1+=(EE1_temp/PART[i].PND/PART[i].get_density())*elas_shear_modulus*mass*dimension;//��񍀂̃� ���x�͂��ꂼ��Ⴄ
			EE2+=((PART[i].P*PART[i].P)/PART[i].get_density())*0.5*mass*elas_lambda;//��O���̃� ���x�͂��ꂼ��Ⴄ
			PART[i].P*=elas_lambda;
			}
		}//if(PART[i].type==ELASTIC || PART[i].type==INELASTIC ||PART[i].type==BOELASTIC)�I��

		//���ꎩ�͓̂����Ȃ���WALL�������x�i���́j��L���Ă���
		//�E�E�EPART[i].P�͂��ׂĂ̗��q���l������K�v�����邪�A�G�l���M�[�v�Z�ɓ�����ׂ��ł͂Ȃ�
		//

	}//for(i=0;i<PART.size();i++)���[�v�I��

	KE*=0.5*mass;
	PE*=-g*mass;
//	EE1*=shear_modulus*mass*dimension;
//	EE2*=0.5*mass*lambda;
	ELAST.set_kinetic(KE);
	ELAST.set_elastic_energy(EE1+EE2);//�S�n�̒e���G�l���M�[
	ELAST.set_elastic_energy1(EE1);
	ELAST.set_elastic_energy2(EE2);
	ELAST.set_potential(PE);
	ELAST.set_hamiltonian(KE+EE1+EE2+PE);
}


void calc_pressure_and_contact_ver_3(vector<mpselastic> &PART, elastic &ELAST)
{
	double dimension=static_cast<double>(ELAST.get_dimension());
	double le=ELAST.get_le();
	double re=ELAST.get_r();
	double density=ELAST.get_density();

	//main��reload_INDEX2()���s���Ă����Ȃ���PART[i]�̏�񂪂߂��Ⴍ����ɂȂ�

	//�c���E���k�����x�ƐڐG���͂̌v�Z
	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{
		//�����z�u�̂��́{
		//(PART[i].PND>PART[i].PND0)�ƂȂ������q���v�Z����
		size_t neighbourN0=PART[i].get_initial_neighboursID().size();//���݈ʒu�ł̎��ӗ��q�����擾
		double qi[4]={PART[i].ang[0], PART[i].ang[1], PART[i].ang[2], PART[i].ang[3]};//���݈ʒu�ł�i�̃N�H�[�^�j�I�����擾

		//�����z�u�ɂ��邩�ǂ����͋C�ɂ��Ȃ�
		for(int k=0;k<neighbourN0;k++)
		{
			int j=PART[i].get_current_neighboursID()[k];//���ݎ��ӂɂ��闱�q��ID���擾

			//���ݔz�u�ł̑��΍��W
			double r_ij[3], r_ji[3];
			for(int D=0;D<3;D++)
			{
				r_ij[D]=PART[j].r[D]-PART[i].r[D];
				r_ji[D]=-r_ij[D];
			}

			double dis=0.0; for(int D=0;D<3;D++){dis+=r_ij[D]*r_ij[D];} dis=sqrt(dis); //���ݗ��q�ԋ���

			//i����݂������z�u�ł̑��΍��W
			//i��ELASTIC�������łȂ����ŏꍇ����
			if(PART[i].type==ELASTIC || PART[i].type==MAGELAST)
			{
				if(PART[j].type==ELASTIC)//�ŏ�����ގ��ŕ�����ƒe���̂ǂ����̏Փ˂ɑΉ��ł��Ȃ��Ȃ�
				{
					//�e���̂̓��͂ł͏d�ݕt����dis0���g���ׂ��i�����z�u����̕ό`���d�v�Ȃ̂Łj
					double r_ij_Init[3], r_ji_Init[3]; 
					for(int D=0;D<3;D++)
					{
						r_ij_Init[D]=PART[j].initial_r[D]-PART[i].initial_r[D];
						r_ji_Init[D]=-r_ij_Init[D];
					}

					double dis0=0.0; for(int D=0;D<3;D++){dis0+=r_ij_Init[D]*r_ij_Init[D];} dis0=sqrt(dis0); //�������q�ԋ���

					double w=kernel(re, dis0);

					double press_accel[3];
					for(int D=0;D<3;D++){
						//�d�݊֐��̉��d���ς̂Ƃ����ύX�Bi�̗��q�����x�Ōv�Z����
						//�d�݂͏����̋������g�������z�͌��ݔz�u�̕����g��
						//press_accel[D]=(PART[i].P*r_ij[D]-PART[j].P*r_ji[D])*w/dis0/dis0;�E�E�E���z��dis0���g���Ă͂����Ȃ��I�I�I
						press_accel[D]=(PART[i].P*r_ij[D]-PART[j].P*r_ji[D])*w/dis/dis;
					}

					PART[i].add_pressure(press_accel);

				}else if(PART[j].type==WALL){
					
					//�ό`���Ȃ����̂ɑ΂��Ă͏d�ݕt����dis���g���ׂ��i�����z�u����̕ό`�͐�������R�͂̌����������d�v�Ȃ̂Łj
					//�ڐG����ꍇ�������߂Â��̂ŋ��E�����̓X�e�b�v���Ƃɕς��

					double w=kernel(re, dis);

					double press_accel[3];
					for(int D=0;D<3;D++){
						press_accel[D]=(PART[i].P*r_ij[D]-PART[j].P*r_ji[D])*w/dis/dis;
					}

					PART[i].add_pressure(press_accel);//���͂̉����x��������
				}				
			}
			else if(PART[i].type==WALL)
			{
				if(PART[j].type==ELASTIC)//�ŏ�����ގ��ŕ�����ƒe���̂ǂ����̏Փ˂ɑΉ��ł��Ȃ��Ȃ�
				{
					double w=kernel(re, dis);

					double press_accel[3];
					for(int D=0;D<3;D++){
						//�d�݊֐��̉��d���ς̂Ƃ����ύX�Bi�̗��q�����x�Ōv�Z����
						//�d�݂͏����̋������g�������z�͌��ݔz�u�̕����g��
						press_accel[D]=(PART[i].P*r_ij[D]-PART[j].P*r_ji[D])*w/dis/dis;
					}
					PART[i].add_pressure(press_accel);
				}		
			}
		}//for(int k=0;k<neighbourN;k++)�E�E�Ek���[�v�I��



		double coef=-dimension/PART[i].get_density()/PART[i].PND;
		PART[i].mul_pressure(coef);
//		ELAST.set_P_visco_stress(D, i, (ELAST.get_P_visco_stress(D, i)*dimension*(-1)/ePND[i]/density));//n0�łȂ�n[i]���g���Ă���E�E�E
	}//	for(int i=0;i<PART.size();i++)
}

void calc_pressure_and_contact(vector<mpselastic> &PART, elastic &ELAST)
{
	double dimension=static_cast<double>(ELAST.get_dimension());
	double le=ELAST.get_le();
	double re=ELAST.get_r();
	double density=ELAST.get_density();

	vector<int>::iterator ip;
	map<int, double>::iterator mp; //map�T���p��iterator

	//main��reload_INDEX2()���s���Ă����Ȃ���PART[i]�̏�񂪂߂��Ⴍ����ɂȂ�

	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{   
//		countOK=countNOT=0;
		size_t neighbourN=PART[i].get_current_neighboursID().size();//���݈ʒu�ł̎��ӗ��q�����擾
		double qi[4]={PART[i].ang[0], PART[i].ang[1], PART[i].ang[2], PART[i].ang[3]};//���݈ʒu�ł�i�̃N�H�[�^�j�I�����擾

		//����re���ɑ��݂��Ȃ����w=0�Ȃ̂Ō���0�Bfind�ŒT���Ηǂ�
		for(int k=0;k<neighbourN;k++)
		{
			int j=PART[i].get_current_neighboursID()[k];//���ݎ��ӂɂ��闱�q��ID���擾

			//���ݔz�u�ł̑��΍��W
			double r_ij[3], r_ji[3];
			for(int D=0;D<3;D++)
			{
				r_ij[D]=PART[j].r[D]-PART[i].r[D];
				r_ji[D]=-r_ij[D];
			}

			double dis=0.0; for(int D=0;D<3;D++){dis+=r_ij[D]*r_ij[D];} dis=sqrt(dis); //���ݗ��q�ԋ���

			//�����z�u��re���ɂ��������q��ID��T��
			ip=find(PART[i].get_initial_neighboursID().begin(), PART[i].get_initial_neighboursID().end(), j);
			
			//�����z�u��re���ɗ��qj������ꍇ
			if(ip!=PART[i].get_initial_neighboursID().end())//���������ꍇ�i���ݎ��ӂɂ͏����z�u��ID������j
			{
				//i����݂������z�u�ł̑��΍��W
				if(PART[j].type==ELASTIC)//�ŏ�����ގ��ŕ�����ƒe���̂ǂ����̏Փ˂ɑΉ��ł��Ȃ��Ȃ�
				{
					double r_ij_Init[3], r_ji_Init[3]; 
					for(int D=0;D<3;D++)
					{
						r_ij_Init[D]=PART[j].initial_r[D]-PART[i].initial_r[D];
						r_ji_Init[D]=-r_ij_Init[D];
					}

					//��]�s��Ri�ɂ��r_ij_Init�̉�]�E�E�E�ȍ~��r_ij_zero���g��
					//qi�͎��X���X�ς��̂Ŏg���񂵂ł��Ȃ��E�E�EPART�̃p�����[�^�Ƃ��ė^����̂��ǂ�
					double r_ij_zero[3], r_ji_zero[3];
					double qj[4]={PART[j].ang[0], PART[j].ang[1], PART[j].ang[2], PART[j].ang[3]};

					rotate_r0(r_ij_Init, qi, r_ij_zero);
					rotate_r0(r_ji_Init, qj, r_ji_zero);

					double dis0=0.0; for(int D=0;D<3;D++){dis0+=r_ij_Init[D]*r_ij_Init[D];} dis0=sqrt(dis0); //�������q�ԋ���
			
					double w=kernel(re, dis0);

					double press_accel[3];
					for(int D=0;D<3;D++)
						press_accel[D]=((PART[i].P*r_ij_zero[D]/PART[i].PND)-(PART[j].P*r_ji_zero[D]/PART[j].PND))*w/dis0/dis0;
					//dis0^2�ŏ����Ă���Ƃ������Ƃ́C���ς��l���Ă���\��������D

					PART[i].add_pressure(press_accel);

				}else if(PART[j].type==WALL){

					double w=kernel(re, dis);

					double press_accel[3];
					for(int D=0;D<3;D++)
						press_accel[D]=((PART[i].P*r_ij[D]/PART[i].PND)-(PART[j].P*r_ji[D]/PART[j].PND))*w/dis/dis;

					PART[i].add_pressure(press_accel);//���͂̉����x��������
				}
			}
			else//�����z�u��re���ɗ��qj���Ȃ��ꍇ
			{
				if(PART[j].type==ELASTIC)//�ŏ�����ގ��ŕ�����ƒe���̂ǂ����̏Փ˂ɑΉ��ł��Ȃ��Ȃ�
				{
					double r_ij_Init[3], r_ji_Init[3]; 
					for(int D=0;D<3;D++)
					{
						r_ij_Init[D]=PART[j].initial_r[D]-PART[i].initial_r[D];
						r_ji_Init[D]=-r_ij_Init[D];
					}

					double dis0=0.0; for(int D=0;D<3;D++){dis0+=r_ij_Init[D]*r_ij_Init[D];} dis0=sqrt(dis0); //�������q�ԋ���

					//��]�s��Ri�ɂ��r_ij_Init�̉�]�E�E�E�ȍ~��r_ij_zero���g��
					//qi�͎��X���X�ς��̂Ŏg���񂵂ł��Ȃ��E�E�EPART�̃p�����[�^�Ƃ��ė^����̂��ǂ�
					double r_ij_zero[3], r_ji_zero[3];
					double qj[4]={PART[j].ang[0], PART[j].ang[1], PART[j].ang[2], PART[j].ang[3]};

					rotate_r0(r_ij_Init, qi, r_ij_zero);
					rotate_r0(r_ji_Init, qj, r_ji_zero);
			
					double w=kernel(re, dis0);

					double press_accel[3];
					for(int D=0;D<3;D++)
						press_accel[D]=((PART[i].P*r_ij_zero[D]/PART[i].PND)-(PART[j].P*r_ji_zero[D]/PART[j].PND))*w/dis0/dis0;

					PART[i].add_pressure(press_accel);

				}else if(PART[j].type==WALL){

					double w=kernel(re, dis);

					double press_accel[3];
					for(int D=0;D<3;D++)
						press_accel[D]=((PART[i].P*r_ij[D]/PART[i].PND)-(PART[j].P*r_ji[D]/PART[j].PND))*w/dis/dis;

					PART[i].add_pressure(press_accel);
				}
			}

		}//for(int k=0;k<neighbourN;k++)�E�E�Ek���[�v�I��

		double coef=-dimension/PART[i].get_density();
		PART[i].mul_pressure(coef);
	}
}

//�p�x�v�Z ������X�V���Ȃ��Ƃ������Ȃ��Ƃ��N����E�E�E
void calc_angular_velocity_and_angle(vector<mpselastic> &PART, elastic &ELAST)
{
	unsigned dimension=ELAST.get_dimension();
	double dt=ELAST.get_dt();

	//�p�x
	if(dimension==2)
	{
		for(int i=0;i<PART.size();i++) PART[i].ang[0]+=dt*PART[i].ang_u[0];
	}
	if(dimension==3)
	{	
		for(int i=0;i<PART.size();i++){
			double omega_norm=sqrt(PART[i].ang_u[0]*PART[i].ang_u[0]+PART[i].ang_u[1]*PART[i].ang_u[1]+PART[i].ang_u[2]*PART[i].ang_u[2]);

			double v[3]={0.0, 0.0, 0.0};//��]���x�N�g��

			if(omega_norm!=0){
				for(int D=0;D<DIMENSION; D++) v[D]=PART[i].ang_u[D]/omega_norm;
			}
			//�N�H�[�^�j�I���C����
			double theta=dt*omega_norm;//(4.33)
			
			double q[4];
			for(int D=0;D<DIMENSION;D++) q[D]=v[D]*sin(theta/2.0);
			q[3]=cos(theta/2.0);
		    
			///�p�x�C��
			double qi[4];
			for(int D=0;D<4;D++) qi[D]=PART[i].ang[D];//���qi�̃N�H�[�^�j�I��
			
			//OK??
			PART[i].ang[0]=q[3]*qi[0]+q[0]*qi[3]+q[1]*qi[2]-q[2]*qi[1];
			PART[i].ang[1]=q[3]*qi[1]+q[1]*qi[3]+q[2]*qi[0]-q[0]*qi[2];
			PART[i].ang[2]=q[3]*qi[2]+q[2]*qi[3]+q[0]*qi[1]-q[1]*qi[0];
			PART[i].ang[3]=q[3]*qi[3]-q[0]*qi[0]-q[1]*qi[1]-q[2]*qi[2];
		}
	}
}

//�p���x�̍X�V
void calc_torque_for_2D(vector<mpselastic> &PART, elastic &ELAST)
{
	double mass=ELAST.get_mass();
	double mag_shear_modulus=ELAST.get_mag_shear_modulus();
	double elas_shear_modulus=ELAST.get_elas_shear_modulus();
	double inertia=ELAST.get_inertia();
	double r=ELAST.get_r();
	double dt=ELAST.get_dt();

	double density=ELAST.get_density();
	double dimension=static_cast<double>(ELAST.get_dimension());
	double force[2];

	for(int i=0;i<PART.size();i++)
	{   
		size_t neighboursN0=PART[i].get_initial_neighboursID().size();
		if(PART[i].type==ELASTIC || PART[i].type==MAGELAST)
		{   
			//�d�v�I�d�݊֐��͉e�����a���̗��q�ɍ�p������INUM�̃��[�v�́u�e�����a���̉��Z�v�Ƃ����Ӗ�������Ii��j�̑��ݍ�p���v�Z
			for(int k=0;k<neighboursN0;k++) //���qi�Ɨ��qi�̉e�����a�Ɋ܂܂�闱�q�̑��ݍ�p���v�Z NUM[i]=PART[i].N;(�����z�u�ɂ�����e���q�̎��ӗ��q��) i!=j�͖�������Ă���I�I
			{
				int j=PART[i].get_initial_neighboursID()[k];
				double X0[2];
				double dis0=0.0;
				PART[i].get_initial_neighbours_position(k, X0); X0[A_Z]=0.0;
				dis0=PART[i].get_initial_distancebps()[k];
				
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double dis=sqrt(X*X+Y*Y);

				if(dis==0) cout<<"dis=0 "<<i<<" "<<j<<endl;

				double theta=(PART[i].ang[0]+PART[j].ang[0])/2;//���Ίp�x (4.8)
				double RX0=X0[0]*cos(theta)-X0[1]*sin(theta);//��]�s��~r0 �������q����]
				double RY0=X0[0]*sin(theta)+X0[1]*cos(theta);//��]�s��~r0

				//���Εψʂ̌v�Z: u_ij=r_ij-R*r0_ij
				double U_X=X-RX0;//�ψʂ�X���� (4.6)
				double U_Y=Y-RY0;//�ψʂ�Y����
			
				double orth_project=(U_X*X+U_Y*Y)/(dis*dis);

				double U_n_X=orth_project*X; //�ψʂ�r(ij)�ɕ��s��"X"����
				double U_n_Y=orth_project*Y; //�ψʂ�r(ij)�ɕ��s��"Y"����

				double U_s_X=U_X-U_n_X; //�ψʂ�r(ij)�ɐ�����"X"����
				double U_s_Y=U_Y-U_n_Y; //�ψʂ�r(ij)�ɐ�����"Y"����

				double EsX=U_s_X/dis0; //�Ђ��݂�r(ij)�ɐ�����"X"����
				double EsY=U_s_Y/dis0; //�Ђ��݂�r(ij)�ɐ�����"Y"����
			
				double r1=ELAST.get_r();
				double w=kernel(r1,dis0);
				double constantF=0;
				if(PART[i].type==MAGELAST){
				constantF=(4*mass*dimension*mag_shear_modulus*w)/(PART[i].get_density()*PART[i].PND*dis0);
				}
				else if(PART[i].type==ELASTIC){
				constantF=(4*mass*dimension*elas_shear_modulus*w)/(PART[i].get_density()*PART[i].PND*dis0);
				}
				force[0]=constantF*EsX;
				force[1]=constantF*EsY;
				double torque=Y*force[0]-X*force[1];

				PART[i].ang_u[0]+=dt*(-0.5)*torque/inertia;
				PART[j].ang_u[0]+=dt*(-0.5)*torque/inertia;
			}
		}
	}
}

//GNUplot�p���̓t�@�C��
void output_stress_for_GNU(vector<mpselastic> &PART, elastic &ELAST, int t)
{
	int dimension=ELAST.get_dimension();
	double le=ELAST.get_distancebp();
	double magni=ELAST.get_times()*le*le;

	if(t==1 || !(t%ELAST.get_interval())){
			char FileName1[128], FileName2[128];

		//�t�@�C�����̏���
			sprintf_s(FileName1, "./stress/normal%04d.dat", t);
			sprintf_s(FileName2, "./stress/shear%04d.dat", t);

			ofstream fout1(FileName1);	//normal
			ofstream fout2(FileName2);	//shear

//			if(dimension==2) for(int i=0;i<PART.size();i++) fout1<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<PART[i].get_normal(A_X, i)*magni<<" "<<PART[i].get_normal(A_Y, i)*magni<<endl;
//			else if(dimension==3) for(int i=0;i<PART.size();i++) if(PART[i].r[A_Y]<le && PART[i].r[A_Y]>-le) fout1<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<PART[i].get_normal(A_X, i)*magni<<"\t"<<PART[i].get_normal(A_Z, i)*magni<<endl;

//			if(dimension==2) for(int i=0;i<PART.size();i++) fout2<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<PART[i].get_shear(A_X, i)*magni<<" "<<PART[i].get_shear(A_Y, i)*magni<<endl;
//			else if(dimension==3) for(int i=0;i<PART.size();i++) if(PART[i].r[A_Y]<le && PART[i].r[A_Y]>-le) fout2<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<PART[i].get_shear(A_X, i)*magni<<"\t"<<PART[i].get_shear(A_Z, i)*magni<<endl;

			fout1.close();
			fout2.close();
	}
}

void quaternion(double P[], double axis[3], double angle){
	double Q[4]={cos(angle/2),axis[0]*sin(angle/2),axis[1]*sin(angle/2),axis[2]*sin(angle/2)};
	double R[4]={cos(angle/2),-axis[0]*sin(angle/2),-axis[1]*sin(angle/2),-axis[2]*sin(angle/2)};
	double QP[4]={-(Q[1]*R[0]+Q[2]*R[1]+Q[3]*R[2])};
}