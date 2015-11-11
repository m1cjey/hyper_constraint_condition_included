#include "stdafx.h"	

#define  A_R 0
#define  A_t 1 //��
		
using namespace std;

void make_cube_region(mpsconfig &CON,vector<point3D> &NODE,int *node, int *divN,double regionX[2],double regionY[2],double regionZ[2]);//�����̉�͗̈�쐬�֐�
void make_cylinder_region(mpsconfig &CON,vector<point3D> &NODE,int *node, int *divN,double* regionR,double* regionZ);//�~����͗̈�쐬�֐�

void MPSTOFEM3D_MRE(mpsconfig &CON,int *node_num,vector<point3D> &NODE,vector<mpselastic> &PART, int fluid_number, int particle_number);//MRE
void MPSTOFEM3D_droplet(mpsconfig &CON,int *node_num,vector<point3D> &NODE,vector<mpselastic> &PART, int fluid_number, int particle_number);//�t�H
void MPSTOFEM3D_nanoe(mpsconfig &CON,int *node_num,vector<point3D> &NODE,vector<mpselastic> &PART, int fluid_number, int particle_number);

//MPS_TO_FEM3Dmain�֐�
void MPS_TO_FEM3D(mpsconfig &CON, int *node_num, vector<point3D> &NODE, vector<mpselastic> &PART, int fluid_number, int particle_number)
{
	int model=CON.get_model_number();

	if(model==3) MPSTOFEM3D_droplet(CON, node_num, NODE,PART, fluid_number, particle_number);//�t�H

	if(model>=5 || model<=8 || model==1 || model==11) MPSTOFEM3D_MRE(CON, node_num, NODE, PART, fluid_number, particle_number); 

	if(model==14) MPSTOFEM3D_nanoe(CON,node_num,NODE,PART, fluid_number, particle_number);//�t�H

	//ofstream fp("node_c.dat");
	//for(int i=1;i<=*node_num;i++) fp<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<endl;
	//fp.close();
	cout<<"MPSTOFEM3D���� �ߓ_��="<<*node_num<<endl;
}

void MPSTOFEM3D_MRE(mpsconfig &CON, int *node_num, vector<point3D> &NODE, vector<mpselastic> &PART, int fluid_number, int particle_number)
{
	//���q�ʒu���W��ߓ_�v�f���W�Ƃ��ĕϊ�����B
	//�G���X�g�}�[����C�̈�ւ̋����ߓ_�E���E�����ݒ聨���΂̏���
	int num=0;//�ߓ_��
	double le=CON.get_distancebp();
    double err=1e-10; //1e-10;
	
	point3D NODE0;
	NODE.clear();
	NODE.push_back(NODE0);//NODE�͐ߓ_�ԍ�1����X�^�[�g���邩��A�����łЂƂm�ۂ��Ă���

    for(unsigned i=0;i<PART.size();i++)
    {
		if(PART[i].type!=WALL){
		//if(PART[i].surface==ON || i%4==0)
		if(PART[i].toFEM==ON)//��͒��AFEM�ɓn�����q�ԍ��𓝈ꂵ���� 2012/05/22
		{
			num++;
			NODE.push_back(NODE0);
			for(int D=0;D<3;D++) NODE[num].r[D]=PART[i].r[D];
			NODE[num].boundary_condition=0;	//���m��
			NODE[num].material=PART[i].type;		//FLUID;
			NODE[num].particleID=i;			//�ߓ_i�ɑΉ����闱�q�ԍ���i
			NODE[num].remesh=ON;			//�����b�V��ON
		}
/*		else{
			num++;
			NODE.push_back(NODE0);
			for(int D=0;D<3;D++) NODE[num].r[D]=PART[i].r[D];
			NODE[num].boundary_condition=0;	//���m��
			NODE[num].material=AIR;		//FLUID;
			NODE[num].particleID=i;			//�ߓ_i�ɑΉ����闱�q�ԍ���i
			NODE[num].remesh=ON;			//�����b�V��ON
		}//*/
		}
    }////////////*/

/*	////�\�ʗ��q����D�e���̂ł͕s�v�H
	for(int i=fluid_number;i<particle_number;i++)
    {
		if(PART[i].type==WALL && PART[i].surface==OFF)//���̂Ɛڂ���INWALL
		{
			num++;
			NODE.push_back(NODE0);
			for(int D=0;D<3;D++) NODE[num].r[D]=PART[i].r[D];
			NODE[num].boundary_condition=0;
			NODE[num].material=ELASTIC;
			NODE[num].particleID=-1;		//�ߓ_i�ɑΉ����闱�q�͑��݂��Ȃ�
			NODE[num].remesh=ON;			//�����b�V��ON
		}
    }*/

	//�~����C�̈�쐬
	if(CON.get_region_shape()==1)
	{
		int divN[3];
		if(CON.get_model_number()==23)	//FEM�Œǉ�����ߓ_�������炷15/2/4
		{
			divN[A_R]=20;
			divN[A_t]=20;
			divN[A_Z]=20;

		}
		else
		{
			divN[A_R]=20;
			divN[A_t]=20;
			divN[A_Z]=20;
		}
		double regionR[2]={0.0, CON.get_RU()};
		double regionZ[2]={CON.get_ZD(), CON.get_ZU()};
		
		int previous_points=num;//�����_�ł̐ߓ_��
		
		make_cylinder_region(CON,NODE,&num,divN,regionR,regionZ);//��͗̈�͉~��

		//���ɁA���܍쐬������͗̈�̏��������ɐߓ_��ǉ�����B�����FINE3D�ŋ��E�����ʂɂ悯���Ȑߓ_���ǉ�����邱�Ƃ�h�����߂ł���B
		double dR=(regionR[1]-regionR[0])/divN[A_R];				//��قǍ쐬������͗̈�̃��b�V������dR
		double dZ=(regionZ[1]-regionZ[0])/divN[A_Z];				//��قǍ쐬������͗̈�̃��b�V������dZ
		//�ȉ��̂悤�ɒ�`�������@�̔��̈���쐬����΁A��͗̈�ƐV�������Ƃ̌��Ԃɂ͗ǎ��Ȑ��l�ʑ̂ɋ߂����b�V�����쐬�����

		regionR[0]=0.0;
		regionR[1]=CON.get_RU()-dR;
		regionZ[0]=CON.get_ZD()+dZ; 
		regionZ[1]=CON.get_ZU()-dZ;

		make_cylinder_region(CON,NODE,&num,divN,regionR,regionZ);//�����̉~���쐬

		//�Œ苫�E�̐ݒ�
		for(int i=previous_points+1;i<=num;i++)
		{
			double X=NODE[i].r[A_X];
			double Y=NODE[i].r[A_Y];
			double Z=NODE[i].r[A_Z];
			double R=sqrt(X*X+Y*Y);

			if(Z>CON.get_ZU()-err) NODE[i].boundary_condition=1;
			else if(Z<CON.get_ZD()+err) NODE[i].boundary_condition=1;
			else if(R>CON.get_RU()-err) NODE[i].boundary_condition=1;
			else NODE[i].boundary_condition=0;
			NODE[i].remesh=OFF;
		}

		//remesh�̈�쐬�E�E�E����if�ŕ�����ׂ��ł́H
		regionR[0]=0.0;
		regionR[1]=CON.get_RU()*0.7;
		regionZ[0]=CON.get_ZD()*0.7;
		regionZ[1]=CON.get_ZU()*0.7;
		previous_points=num;//�����_�ł̐ߓ_��

		make_cylinder_region(CON,NODE,&num,divN,regionR,regionZ);//��͗̈�͉~��

		for(int i=previous_points+1;i<=num;i++)
		{
			NODE[i].boundary_condition=0;			//���E����
			NODE[i].remesh=ON;
		}

	}else if(CON.get_region_shape()==0){
		cout<<"�����̉�͗̈�͖�����"<<endl;
		cout<<"�v�Z�I��"<<endl;
		exit(1);
	}

	//����
	point3D NODE01;
	int divN[3];
	if(CON.get_model_number()==23)	//�ߓ_�������炷15/2/4
	{
		divN[A_R]=20;//���a����������
		divN[A_t]=20;//�p�x�����������i�����p�`�ŋߎ����邩�j
		divN[A_Z]=20;//��������������
	}
	else
	{
		divN[A_R]=20;//���a����������
		divN[A_t]=20;//�p�x�����������i�����p�`�ŋߎ����邩�j
		divN[A_Z]=20;//��������������
	}
	double Rmin=0.0;
	double Rmax=CON.get_magnet_r();

	double Zmin=CON.get_magnet_Z()-0.5*CON.get_magnet_H();
	double Zmax=CON.get_magnet_Z()+0.5*CON.get_magnet_H();
	double divL[3];
	divL[A_R]=(Rmax-Rmin)/divN[A_R];
	divL[A_t]=(2*PI)/divN[A_t];
	divL[A_Z]=(Zmax-Zmin)/divN[A_Z];
	double theta2=CON.get_magnet_angle()*PI*2/360;//���΂̉�]�p�x
	if(CON.get_EM_calc_type()==2){	//�i�v���ΐÎ���
	for(int n=0;n<=divN[A_R];n++)
	{
		if(n==0)//���S
		{
			for(int k=0;k<=divN[A_Z];k++)
			{
				num++;
				NODE.push_back(NODE01);
				NODE[num].r[A_X]=0;
				NODE[num].r[A_Y]=0;
				NODE[num].r[A_Z]=Zmin+divL[A_Z]*k;
				NODE[num].material=MAGNET;
				NODE[num].particleID=-1;					//�Ή����闱�q�����݂��Ȃ�����-1���i�[
				NODE[num].boundary_condition=0;			//���E����
				NODE[num].remesh=OFF;

				if(fabs(theta2)>0.0175)//double��0�Ɣ�r���Ă͂����Ȃ��i��v����킯���Ȃ��j
				{
					double XX=NODE[num].r[A_X];
					double ZZ=NODE[num].r[A_Z]-CON.get_magnet_Z();
					double newX=XX*cos(theta2)-ZZ*sin(theta2);
					double newZ=XX*sin(theta2)+ZZ*cos(theta2)+CON.get_magnet_Z();
					NODE[num].r[A_X]=newX;
					NODE[num].r[A_Z]=newZ;
				}
			}
		}
		else
		{
			for(int m=0;m<divN[A_t];m++)
			{
				for(int k=0;k<=divN[A_Z];k++)
				{
					num++;
					NODE.push_back(NODE01);
					double r=divL[A_R]*n;
					double theta=divL[A_t]*m;
					NODE[num].r[A_X]=r*cos(theta);
					NODE[num].r[A_Y]=r*sin(theta);
					NODE[num].r[A_Z]=Zmin+divL[A_Z]*k;
					NODE[num].material=MAGNET;
					NODE[num].particleID=-1;					//�Ή����闱�q�����݂��Ȃ�����-1���i�[
					NODE[num].boundary_condition=0;			//���E����
					NODE[num].remesh=OFF;

					if(fabs(theta2)>0.0175)//���΂��X���Ă���ꍇ(1deg���傫���ꍇ)�E�E�Edouble��0�Ɣ�r���Ă͂����Ȃ�
					//if(theta2!=0)����̓_���I�I
					{
						double XX=NODE[num].r[A_X];
						double ZZ=NODE[num].r[A_Z]-CON.get_magnet_Z();
						double newX=XX*cos(theta2)-ZZ*sin(theta2);
						double newZ=XX*sin(theta2)+ZZ*cos(theta2)+CON.get_magnet_Z();
						NODE[num].r[A_X]=newX;
						NODE[num].r[A_Z]=newZ;
					}
				}
			}
		}
	}
	}
	else if(CON.get_EM_calc_type()==3){	//�d���΁@������
	for(int n=0;n<=divN[A_R];n++)
	{
		if(n==0)//���S
		{
			for(int k=0;k<=divN[A_Z];k++)
			{
				num++;
				NODE.push_back(NODE01);
				NODE[num].r[A_X]=0;
				NODE[num].r[A_Y]=0;
				NODE[num].r[A_Z]=Zmin+divL[A_Z]*k;
				NODE[num].material=IRON;
				NODE[num].particleID=-1;					//�Ή����闱�q�����݂��Ȃ�����-1���i�[
				NODE[num].boundary_condition=0;			//���E����
				NODE[num].remesh=OFF;

				if(fabs(theta2)>0.0175)//double��0�Ɣ�r���Ă͂����Ȃ��i��v����킯���Ȃ��j
				{
					double XX=NODE[num].r[A_X];
					double ZZ=NODE[num].r[A_Z]-CON.get_magnet_Z();
					double newX=XX*cos(theta2)-ZZ*sin(theta2);
					double newZ=XX*sin(theta2)+ZZ*cos(theta2)+CON.get_magnet_Z();
					NODE[num].r[A_X]=newX;
					NODE[num].r[A_Z]=newZ;
				}
			}
		}
		else
		{
			for(int m=0;m<divN[A_t];m++)
			{
				for(int k=0;k<=divN[A_Z];k++)
				{
					num++;
					NODE.push_back(NODE01);
					double r=divL[A_R]*n;
					double theta=divL[A_t]*m;
					NODE[num].r[A_X]=r*cos(theta);
					NODE[num].r[A_Y]=r*sin(theta);
					NODE[num].r[A_Z]=Zmin+divL[A_Z]*k;
					if(n<6){
					NODE[num].material=IRON;	
					NODE[num].boundary_condition=0;			//���E����
					}
					else {
					NODE[num].material=COIL;
					 if((n==divN[A_R]) && (k%3==0))NODE[num].boundary_condition=23;//�O������
					else if((n==divN[A_R]) && (k%3==1))NODE[num].boundary_condition=22;
					else if((n==divN[A_R]) && (k%3==2))NODE[num].boundary_condition=21;
		/*			else if((n==6) && (k%3==0)) NODE[num].boundary_condition=21;	//�S�S�Ƃ̋��E��
					else if((n==6) && (k%3==1))NODE[num].boundary_condition=22;
					else if((n==6) && (k%3==2))NODE[num].boundary_condition=23;  */
					else NODE[num].boundary_condition=0;			//���E����
					}
					NODE[num].particleID=-1;					//�Ή����闱�q�����݂��Ȃ�����-1���i�[
					NODE[num].remesh=OFF;

					if(fabs(theta2)>0.0175)//���΂��X���Ă���ꍇ(1deg���傫���ꍇ)�E�E�Edouble��0�Ɣ�r���Ă͂����Ȃ�
					//if(theta2!=0)����̓_���I�I
					{
						double XX=NODE[num].r[A_X];
						double ZZ=NODE[num].r[A_Z]-CON.get_magnet_Z();
						double newX=XX*cos(theta2)-ZZ*sin(theta2);
						double newZ=XX*sin(theta2)+ZZ*cos(theta2)+CON.get_magnet_Z();
						NODE[num].r[A_X]=newX;
						NODE[num].r[A_Z]=newZ;
					}
				}
			}
		}
	}
	}
	//���΂̊O���ɋ�C�w�ݒu
	//�S�̂̋��E�Ɣ��Ȃ��悤�ɋC�����邱��
	double regionR[2];
	double regionZ[2];
	int number_of_layers=1+CON.get_magnetic_layers();
	divN[2]=22;
	regionR[0]=0.0;//�ύX�Ȃ�
	
	for(int layers=1;layers<=number_of_layers;layers++)
	{
		regionR[1]=Rmax+(divL[A_R])*layers;
		regionZ[0]=Zmin-(divL[A_Z])*layers/2;
		regionZ[1]=Zmax+(divL[A_Z])*layers/2;
		int before=num;//�����_�ł̐ߓ_��

		make_cylinder_region(CON,NODE,&num,divN,regionR,regionZ);//��͗̈�͉~��
		
		for(int i=before+1;i<=num;i++)
		{
			NODE[i].boundary_condition=0;			//���E����
			NODE[i].remesh=OFF;
		}

		if(fabs(theta2)>0.0175)//���΂��X���Ă���ꍇ(1deg���傫���ꍇ)�E�E�Edouble��0�Ɣ�r���Ă͂����Ȃ�
		{
			for(int i=before +1;i<=num;i++)
			{
				double XX=NODE[i].r[A_X];
				double ZZ=NODE[i].r[A_Z]-CON.get_magnet_Z();
				double newX=XX*cos(theta2)-ZZ*sin(theta2);
				double newZ=XX*sin(theta2)+ZZ*cos(theta2)+CON.get_magnet_Z();
				NODE[i].r[A_X]=newX;
				NODE[i].r[A_Z]=newZ;
			}
		}
	}

	/*//���΂̊O���ɋ�C�w�ݒu
	double regionR[2]={0,Rmax+divL[A_R]};//{0,Rmax+divL[A_R]/2};
	double regionZ[2]={Zmin-divL[A_Z],Zmax+divL[A_Z]};//{Zmin-divL[A_Z]/2,Zmax+divL[A_Z]/2};
	int count=num;//�����_�ł̐ߓ_��
		
	make_cylinder_region(CON,NODE,&num,divN,regionR,regionZ);//��͗̈�͉~��

	for(int i=count+1;i<=num;i++)
	{
		NODE[i].boundary_condition=0;			//���E����
		NODE[i].remesh=OFF;
	}

	if(theta2!=0)
	{
		for(int i=count+1;i<=num;i++)
		{
			double XX=NODE[i].r[A_X];
			double ZZ=NODE[i].r[A_Z]-CON.get_magnet_Z();
			double newX=XX*cos(theta2)-ZZ*sin(theta2);
			double newZ=XX*sin(theta2)+ZZ*cos(theta2)+CON.get_magnet_Z();
			NODE[i].r[A_X]=newX;
			NODE[i].r[A_Z]=newZ;
		}
	}*/

	*node_num=num;
}

//�t�H
void MPSTOFEM3D_droplet(mpsconfig &CON,int *node_num,vector<point3D> &NODE,vector<mpselastic> &PART, int fluid_number, int particle_number)
{
	int num=0;						//�ߓ_��
	double le=CON.get_distancebp();
    double err=1e-10;
	
	point3D NODE0;
	NODE.clear();
	NODE.push_back(NODE0);			//NODE�͐ߓ_�ԍ�1����X�^�[�g���邩��A�����łЂƂm�ۂ��Ă���

	////���̗��q�o��
    for(int i=0;i<fluid_number;i++)
    {
		//if(PART[i].surface==ON || i%4==0)
		//if(PART[i].toFEM==ON)//��͒��AFEM�ɓn�����q�ԍ��𓝈ꂵ����
		{
			num++;
			NODE.push_back(NODE0);
			for(int D=0;D<3;D++) NODE[num].r[D]=PART[i].r[D];
			NODE[num].boundary_condition=0;
			NODE[num].material=FLUID;
			NODE[num].particleID=i;				//�ߓ_i�ɑΉ����闱�q�ԍ���i
			NODE[num].remesh=ON;			//�����b�V��ON
		}
    }////////////*/

	/*///////�\�ʗ��q�̖@���޸��
    double *direct[DIMENTION];
    for(int D=0;D<DIMENTION;D++) direct[D]=new double [fluid_number];
    
  //  double inpoint[3];	//FEM�ɑ�������ߓ_���W(X,Y,Z)
	int count=0;			
	double Le=CON.get_distancebp()*0.5;
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].surface==ON)
		{
			count++;
			//if(count%2==0)//�Ԉ���
			{
				direct_f(CON,PART,i,direct);
				for(int D=0;D<DIMENTION;D++) direct[D][i]*=-1;//�O�����@���޸�ق��~�������畄����������
				
				for(int n=1;n<=4;n+=1)
				{   
					num++;
					NODE.push_back(NODE0);
					for(int D=0;D<3;D++) NODE[num].r[D]=PART[i].r[D]+direct[D][i]*Le*n;
					NODE[num].boundary_condition=0;
					NODE[num].material=AIR;
					NODE[num].particleID=-1;				
					NODE[num].remesh=ON;			//�����b�V��ON
				}
			}
		}
	}
	for(int D=0;D<DIMENTION;D++) delete [] direct[D];
	//////*/

	if(CON.get_region_shape()==0)
	{
		int divN[3];
		divN[A_X]=10;
		divN[A_Y]=10;
		divN[A_Z]=10;
		double regionX[2]={CON.get_XL(),CON.get_XR()};
		double regionY[2]={CON.get_YD(),CON.get_YU()};
		double regionZ[2]={CON.get_ZD(),CON.get_ZU()};
		int count=num;//�����_�ł̐ߓ_��
		
		make_cube_region(CON,NODE,&num,divN,regionX,regionY,regionZ);//��͗̈�͒�����

		//���ɁA���܍쐬������͗̈�̏��������ɐߓ_��ǉ�����B�����FINE3D�ŋ��E�����ʂɂ悯���Ȑߓ_���ǉ�����邱�Ƃ�h�����߂ł���B
		double dX=(regionX[1]-regionX[0])/divN[A_X];				//��قǍ쐬������͗̈�̃��b�V������dX
		double dY=(regionY[1]-regionY[0])/divN[A_Y];				//��قǍ쐬������͗̈�̃��b�V������dY
		double dZ=(regionZ[1]-regionZ[0])/divN[A_Z];				//��قǍ쐬������͗̈�̃��b�V������dZ
		//�ȉ��̂悤�ɒ�`�������@�̔��̈���쐬����΁A��͗̈�ƐV�������Ƃ̌��Ԃɂ͗ǎ��Ȑ��l�ʑ̂ɋ߂����b�V�����쐬�����
		regionX[0]=CON.get_XL()+dX*sqrt(3.0)*0.5; regionX[1]=CON.get_XR()-dX*sqrt(3.0)*0.5;
		regionY[0]=CON.get_YD()+dY*sqrt(3.0)*0.5; regionY[1]=CON.get_YU()-dY*sqrt(3.0)*0.5;
		regionZ[0]=CON.get_ZD()+dZ*sqrt(3.0)*0.5; regionZ[1]=CON.get_ZU()-dZ*sqrt(3.0)*0.5;

		make_cube_region(CON,NODE,&num,divN,regionX,regionY,regionZ);//�����̔��쐬

		for(int i=count+1;i<=num;i++)				//���E����
		{
			double Z=NODE[i].r[A_Z];
			if(Z>CON.get_ZU()-err) NODE[i].boundary_condition=2;
			else if(Z<CON.get_ZD()+err) NODE[i].boundary_condition=1;
			else NODE[i].boundary_condition=0;
			NODE[i].remesh=OFF;
		}

		//remesh�̈�쐬
		regionX[0]=CON.get_XL()*0.5; regionX[1]=CON.get_XR()*0.5;
		regionY[0]=CON.get_YD()*0.5; regionY[1]=CON.get_YU()*0.5;
		regionZ[0]=CON.get_ZD()*0.7; regionZ[1]=CON.get_ZU()*0.7;
		count=num;//�����_�ł̐ߓ_��
		make_cube_region(CON,NODE,&num,divN,regionX,regionY,regionZ);//��͗̈�͒�����

		for(int i=count+1;i<=num;i++)
		{
			NODE[i].boundary_condition=0;			//���E����
			NODE[i].remesh=ON;
		}
		
	}
	else if(CON.get_region_shape()==1)
	{
		int divN[3];
		divN[A_R]=10;
		divN[A_t]=30;
		divN[A_Z]=10;
		double regionR[2]={0,CON.get_RU()};
		double regionZ[2]={CON.get_ZD(),CON.get_ZU()};
		int count=num;//�����_�ł̐ߓ_��
		
		make_cylinder_region(CON,NODE,&num,divN,regionR,regionZ);//��͗̈�͉~��

		//���ɁA���܍쐬������͗̈�̏��������ɐߓ_��ǉ�����B�����FINE3D�ŋ��E�����ʂɂ悯���Ȑߓ_���ǉ�����邱�Ƃ�h�����߂ł���B
		double dR=(regionR[1]-regionR[0])/divN[A_R];				//��قǍ쐬������͗̈�̃��b�V������dR
		double dZ=(regionZ[1]-regionZ[0])/divN[A_Z];				//��قǍ쐬������͗̈�̃��b�V������dZ
		//�ȉ��̂悤�ɒ�`�������@�̔��̈���쐬����΁A��͗̈�ƐV�������Ƃ̌��Ԃɂ͗ǎ��Ȑ��l�ʑ̂ɋ߂����b�V�����쐬�����
		regionR[0]=0; regionR[1]=CON.get_RU()-dR;
		regionZ[0]=CON.get_ZD()+dZ; regionZ[1]=CON.get_ZU()-dZ;

		make_cylinder_region(CON,NODE,&num,divN,regionR,regionZ);//�����̉~���쐬

		for(int i=count+1;i<=num;i++)				//���E����
		{
			double X=NODE[i].r[A_X];
			double Y=NODE[i].r[A_Y];
			double Z=NODE[i].r[A_Z];
			double R=sqrt(X*X+Y*Y);
			if(Z>CON.get_ZU()-err) NODE[i].boundary_condition=1;
			else if(Z<CON.get_ZD()+err) NODE[i].boundary_condition=1;
			else if(R>CON.get_RU()-err) NODE[i].boundary_condition=1;
			else NODE[i].boundary_condition=0;
			NODE[i].remesh=OFF;
		}

		//remesh�̈�쐬
		regionR[0]=0; regionR[1]=CON.get_RU()*0.7;
		regionZ[0]=CON.get_ZD()*0.7; regionZ[1]=CON.get_ZU()*0.7;
		count=num;//�����_�ł̐ߓ_��
		make_cylinder_region(CON,NODE,&num,divN,regionR,regionZ);//��͗̈�͉~��

		for(int i=count+1;i<=num;i++)
		{
			NODE[i].boundary_condition=0;			//���E����
			NODE[i].remesh=ON;
		}
	}///*/

	

	

	*node_num=num;

}

///3D�p�Ód����
void MPSTOFEM3D_nanoe(mpsconfig &CON,int *node_num,vector<point3D> &NODE,vector<mpselastic> &PART, int fluid_number, int particle_number)
{
    double le=CON.get_distancebp();	//�������q�ԋ���
    double P_R=0.0005;					//�d�ɔ��a
//    double dX,dY,dZ;					//�e�������̕�����
//    int Nx,Ny,Nz;						//�e�������̕�����
    double ex_r=0.0006;					//���̂�����Ƒz�肵�Ă�̈攼�a
//    double hight;						//���̂�����Ƒz�肵�Ă鍂��
	int TOUCH=0;						//�ڐG
	int UNTOUCH=1;						//��ڐG
	double err=1e-10;
	double beta=CON.get_beta();
	int divN[3];

	int num=0;						//�ߓ_��
	
	point3D NODE0;
	NODE.clear();
	NODE.push_back(NODE0);			//NODE�͐ߓ_�ԍ�1����X�^�[�g���邩��A�����łЂƂm�ۂ��Ă���

	////���̗��q�o��
    for(int i=0;i<fluid_number;i++)
    {
       // if(PART[i].surface==ON || i%4==0)
		//if(PART[i].toFEM==ON)//��͒��AFEM�ɓn�����q�ԍ��𓝈ꂵ����
		{
			num++;
			NODE.push_back(NODE0);
			for(int D=0;D<3;D++) NODE[num].r[D]=PART[i].r[D];
			NODE[num].boundary_condition=0;
			NODE[num].material=FLUID;
			NODE[num].particleID=i;				//�ߓ_i�ɑΉ����闱�q�ԍ���i
			NODE[num].remesh=ON;			//�����b�V��ON
		}
    }////////////*/

    ///�d�ɗ��q�o��
	double minZ=0;		//�o�͂���Ǘ��q�̍ŏ�����
    for(int i=fluid_number;i<particle_number;i++)
    {   
       // if(PART[i].type==BDWALL || PART[i].type==INWALL || i%5==0)
		{
			num++;
			NODE.push_back(NODE0);
			for(int D=0;D<3;D++) NODE[num].r[D]=PART[i].r[D];
			NODE[num].boundary_condition=0;
			NODE[num].material=ELECTRODE;
			NODE[num].particleID=i;				//�ߓ_i�ɑΉ����闱�q�ԍ���i
			NODE[num].remesh=ON;			//�����b�V��ON*/
			if(minZ>PART[i].r[A_Z]) minZ=PART[i].r[A_Z];
        }
    }////////////////*/

	//���d�ɍ쐬
	divN[A_R]=10;
	divN[A_t]=20;
	divN[A_Z]=200;
	double regionR[2]={0,P_R+le};
	double regionZ[2]={CON.get_ZD(),minZ-3*le};		//�쐬����~���̍���
	int count=num;//�����_�ł̐ߓ_��
		
	make_cylinder_region(CON,NODE,&num,divN,regionR,regionZ);//��͗̈�͉~��

	for(int i=count+1;i<=num;i++)				//���E����
	{
		NODE[i].material=ELECTRODE;
		NODE[i].boundary_condition=1;
		NODE[i].remesh=OFF;
	}
    
    
   /* double WX=CON.get_XR()-CON.get_XL();///��͗̈�X������
    double WY=CON.get_YU()-CON.get_YD();///��͗̈�Y������
    double WZ=CON.get_ZU()-CON.get_ZD();///��͗̈�Z������
    
    /////��͗̈�� (���̂͂��Ȃ��Ƒz�肵�Ă���̈�)
    
	//��C�̈�1
    Nx=40;//120;//15;
    Ny=40;//120;//15;
    Nz=75*2;//60;//75;//15;
	double Xmin1=CON.get_XL()*3/20; double Xmax1=CON.get_XR()*3/20;//��C�̈�1�̍ŏ��E�ő�X���W
	double Ymin1=CON.get_YD()*3/20; double Ymax1=CON.get_YU()*3/20;//��C�̈�1�̍ŏ��E�ő�Y���W
	double Zmin1=CON.get_ZD(); double Zmax1=CON.get_ZU()*0.24;//��C�̈�1�̍ŏ��E�ő�Z���W
	dX=(Xmax1-Xmin1)/Nx;
    dY=(Ymax1-Ymin1)/Ny;
    dZ=(Zmax1-Zmin1)/Nz;

	
    for(double X=Xmin1;X<=Xmax1+err;X+=dX)
    {   
        for(double Y=Ymin1;Y<=Ymax1+err;Y+=dY)
		{   
			for(double Z=Zmin1;Z<=Zmax1+err;Z+=dZ)
			{
				if(X*X+Y*Y>ex_r*ex_r || Z>=hight)
				{
					num++;
					int BC=0;///���E����
					if(Z>=CON.get_ZU()-err) BC=2;
					if(Z<CON.get_ZD()+err) BC=1;
					fout<<X<<" "<<Y<<" "<<Z<<" "<<AIR<<" "<<BC<<endl;//�Ō�̐����͋��E����
					fprintf(fp,"%lf %lf %lf %d %d\n",X,Y,Z,AIR,BC);//�Ō�̐����͋��E����
				}
			} 
		}
    }////////*/

	/*/��C�̈�2
	Nx=120/2;//15;
    Ny=120/2;//15;
    Nz=75/2;//15;
    dX=WX/Nx;
    dY=WY/Ny;
    dZ=WZ/Nz;
	
    for(double X=CON.get_XL();X<=CON.get_XR()+err;X+=dX)
    {   
        for(double Y=CON.get_YD();Y<=CON.get_YU()+err;Y+=dY)
		{   
			for(double Z=CON.get_ZD();Z<=CON.get_ZU()+err;Z+=dZ)
			{
				int flag=1;//flag=1�Ȃ��C�̈�1�@2�Ȃ��C�̈�2
				if(X>Xmax1+le || X<Xmin1-le) flag=2;
				if(Y>Ymax1+le || Y<Ymin1-le) flag=2;
				if(Z>Zmax1+le || Z<Zmin1-le) flag=2;
				if(flag==2)
				{
					num++;
					int BC=0;///���E����
					if(Z>=CON.get_ZU()-err) BC=2;
					if(Z<CON.get_ZD()+err) BC=1;
					fout<<X<<" "<<Y<<" "<<Z<<" "<<AIR<<" "<<BC<<endl;//�Ō�̐����͋��E����
					fprintf(fp,"%lf %lf %lf %d %d\n",X,Y,Z,AIR,BC);//�Ō�̐����͋��E����
				}
			} 
		}
    }////////*/

	*node_num=num;
   
}


//�����̉�͗̈�쐬�֐�
void make_cube_region(mpsconfig &CON,vector<point3D> &NODE,int *node, int *divN,double regionX[2],double regionY[2],double regionZ[2])
{
	int node_num=*node;
	//divN[3];								//�e�ӂ̕�����
	
	double Xmin=regionX[0];				//��͗̈�
	double Xmax=regionX[1];
	double Ymin=regionY[0];
	double Ymax=regionY[1];
	double Zmin=regionZ[0];
	double Zmax=regionZ[1];

	double divL[3];								//������
	divL[A_X]=(Xmax-Xmin)/divN[A_X];
	divL[A_Y]=(Ymax-Ymin)/divN[A_Y];
	divL[A_Z]=(Zmax-Zmin)/divN[A_Z];

	point3D NODE01;

	//���
	for(int n=0;n<=divN[A_X];n++)
	{
		for(int m=0;m<=divN[A_Y];m++)
		{
			node_num++;
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmin+n*divL[A_X];
			NODE[node_num].r[A_Y]=Ymin+m*divL[A_Y];
			NODE[node_num].r[A_Z]=Zmin;					//��͗̈�̒��
			NODE[node_num].material=AIR;
			NODE[node_num].particleID=-1;					//�Ή����闱�q�����݂��Ȃ�����-1���i�[
		}
	}

	//���
	for(int n=0;n<=divN[A_X];n++)
	{
		for(int m=0;m<=divN[A_Y];m++)
		{
			node_num++;
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmin+n*divL[A_X];
			NODE[node_num].r[A_Y]=Ymin+m*divL[A_Y];
			NODE[node_num].r[A_Z]=Zmax;					//��͗̈�̏��
			NODE[node_num].material=AIR;
			NODE[node_num].particleID=-1;					//�Ή����闱�q�����݂��Ȃ�����-1���i�[
		}
	}////

	
	//����Y
	for(int n=0;n<=divN[A_X];n++)
	{
		for(int m=1;m<divN[A_Z];m++)
		{
			node_num++;
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmin+n*divL[A_X];
			NODE[node_num].r[A_Y]=Ymin;
			NODE[node_num].r[A_Z]=Zmin+m*divL[A_Z];		
			NODE[node_num].material=AIR;
			NODE[node_num].particleID=-1;					//�Ή����闱�q�����݂��Ȃ�����-1���i�[
		}
	}
	for(int n=0;n<=divN[A_X];n++)
	{
		for(int m=1;m<divN[A_Z];m++)
		{
			node_num++;
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmin+n*divL[A_X];
			NODE[node_num].r[A_Y]=Ymax;
			NODE[node_num].r[A_Z]=Zmin+m*divL[A_Z];	
			NODE[node_num].material=AIR;
			NODE[node_num].particleID=-1;					//�Ή����闱�q�����݂��Ȃ�����-1���i�[
		}
	}

	//����
	for(int n=1;n<divN[A_Y];n++)
	{
		for(int m=1;m<divN[A_Z];m++)
		{
			node_num++;
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmin;
			NODE[node_num].r[A_Y]=Ymin+n*divL[A_Y];
			NODE[node_num].r[A_Z]=Zmin+m*divL[A_Z];	
			NODE[node_num].material=AIR;
			NODE[node_num].particleID=-1;					//�Ή����闱�q�����݂��Ȃ�����-1���i�[
		}
	}
	for(int n=1;n<divN[A_Y];n++)
	{
		for(int m=1;m<divN[A_Z];m++)
		{
			node_num++;
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmax;
			NODE[node_num].r[A_Y]=Ymin+n*divL[A_Y];
			NODE[node_num].r[A_Z]=Zmin+m*divL[A_Z];	
			NODE[node_num].material=AIR;
			NODE[node_num].particleID=-1;					//�Ή����闱�q�����݂��Ȃ�����-1���i�[
		}
	}
	*node=node_num;
}

//�~����͗̈�쐬�֐�
void make_cylinder_region(mpsconfig &CON,vector<point3D> &NODE,int *node, int *divN, double* regionR, double* regionZ)
{
	int node_num=*node;
	//divN[3];					//�e�ӂ̕�����
	
	double Rmin=0;				//��͗̈�
	double Rmax=regionR[1];
	double Zmin=regionZ[0];
	double Zmax=regionZ[1];

	double divL[3];//������
	divL[A_R]=(Rmax-Rmin)/divN[A_R];
	divL[A_t]=(2*PI)/divN[A_t];
	divL[A_Z]=(Zmax-Zmin)/divN[A_Z];

	point3D NODE01;

	/////////////////////���
	node_num++;					//���S�_
	NODE.push_back(NODE01);
	NODE[node_num].r[A_X]=0;
	NODE[node_num].r[A_Y]=0;
	NODE[node_num].r[A_Z]=Zmin;					//��͗̈�̒��
	NODE[node_num].material=AIR;
	NODE[node_num].particleID=-1;

	for(int n=1;n<=divN[A_R];n++)
	{
		for(int m=0;m<divN[A_t];m++)
		{
			node_num++;
			NODE.push_back(NODE01);
			double r=divL[A_R]*n;
			double theta=divL[A_t]*m;
			NODE[node_num].r[A_X]=r*cos(theta);
			NODE[node_num].r[A_Y]=r*sin(theta);
			NODE[node_num].r[A_Z]=Zmin;					//��͗̈�̒��
			NODE[node_num].material=AIR;
			NODE[node_num].particleID=-1;					//�Ή����闱�q�����݂��Ȃ�����-1���i�[
		}
	}
	//////////////////////////////////���
	node_num++;					//���S�_
	NODE.push_back(NODE01);
	NODE[node_num].r[A_X]=0;
	NODE[node_num].r[A_Y]=0;
	NODE[node_num].r[A_Z]=Zmax;					//��͗̈�̒��
	NODE[node_num].material=AIR;
	NODE[node_num].particleID=-1;	
	for(int n=1;n<=divN[A_R];n++)
	{
		for(int m=0;m<divN[A_t];m++)
		{
			node_num++;
			NODE.push_back(NODE01);
			double r=divL[A_R]*n;
			double theta=divL[A_t]*m;
			NODE[node_num].r[A_X]=r*cos(theta);
			NODE[node_num].r[A_Y]=r*sin(theta);
			NODE[node_num].r[A_Z]=Zmax;					//��͗̈�̒��
			NODE[node_num].material=AIR;
			NODE[node_num].particleID=-1;					//�Ή����闱�q�����݂��Ȃ�����-1���i�[
		}
	}

	//����
	double RR=divL[A_R]*divN[A_R];
	for(int n=1;n<divN[A_Z];n++)
	{
		for(int m=0;m<divN[A_t];m++)
		{
			node_num++;
			NODE.push_back(NODE01);
			double theta=divL[A_t]*m;
			NODE[node_num].r[A_X]=RR*cos(theta);
			NODE[node_num].r[A_Y]=RR*sin(theta);;
			NODE[node_num].r[A_Z]=Zmin+n*divL[A_Z];		
			NODE[node_num].material=AIR;
			NODE[node_num].particleID=-1;					//�Ή����闱�q�����݂��Ȃ�����-1���i�[
		}
	}
	*node=node_num;
}