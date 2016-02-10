#include"stdafx.h"	
using namespace std;

//�\�ʔ���֐�
//***********�n�b�V���T������Α����Ȃ�H�H�H
//�e���̌v�Z�p�ɃJ�X�^�}�C�Y�i2012-11-05�j
void freeon(elastic &ELAST, vector<mpselastic> &PART, double n0_4,int *INDEX,int **MESH, double *mindis, int t)
{
    //����:freeon�֐��ł͊e���q�̗��q�����x�����߂�ƂƂ��ɁA����ꂽ���x����\�ʔ�����s���B�܂��A�Œᗱ�q�ԋ��������łɂ��Ƃ߂Ă���
	//freeon2���x�����A���̂Ԃ���񉻂��e�ՁB
	double re=ELAST.get_re()*ELAST.get_distancebp();	//get_re()�̒l�����̂܂܎g���Ă̓_���I�I
	double mass=ELAST.get_particle_mass();
//	double volume=4*3.141592653*pow(re, 3)/3;
	int SIZE=ELAST.get_X_mesh()*ELAST.get_Y_mesh();
	
	//���x�����߁A�\�ʔ�����s���B�\�ʂȂ�P=0�ɂ���
	//omp_set_num_threads(8);//�X���b�h���w��

	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{
		PART[i].PND=0;
		PART[i].PND2=0;
		PART[i].N=0;
		PART[i].N2=0;
		PART[i].reset_current_distancebps();
		PART[i].reset_current_neighboursID();
		PART[i].reset_current_neighbours_position();
		PART[i].distance.clear();

		////���q�����x����
		double pnd=0.0;	//�S���q���l���ɂ��ꂽ���q�����x
		double pnd2=0;	//�e���̌v�Z�p���q�����x
		double pnd4=0;	//�\�ʔ���p
		double sum=0;
		int N=0;		//���ӗ��q��
		int N0=0;
			
		//��ԂɌŗL��ID(INDEX, MESH)�Ɨ��q�ɌŗL��ID(PART[i].index, PART[i].NEI[N])���ƍ����ĒT������
		for(int I=PART[i].index-1;I<=PART[i].index+1;I++)
		{
			if(PART[i].index>=SIZE && PART[i].index<ELAST.get_number_of_mesh()-SIZE)
			{       
				for(int J=-1*ELAST.get_X_mesh();J<=ELAST.get_X_mesh();J+=ELAST.get_X_mesh())
				{
					for(int K=-1*SIZE;K<=SIZE;K+=SIZE)
					{
						for(int L=0;L<INDEX[I+J+K];L++)
						{       
							int j=MESH[I+J+K][L];	//(I+J+K)�Ԗڂ̊i�q�Ɋ܂܂��L�Ԗڂ̗��q��ID���擾����
							if(j!=i)				//���qID�����ݒ��ڂ��Ă������(i)�ƈقȂ�΋������v�Z����
							{ 
								double X=PART[j].r[A_X]-PART[i].r[A_X];
								double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
								double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
								double dis=sqrt(X*X+Y*Y+Z*Z);
						
								if(dis<re)	
								{
									pnd+=kernel(re, dis);	//�v�Z�����������e�����a���ł���Ώd�݊֐����v�Z����
									PART[i].NEI[N]=j;		//i�̉e�����a�ɑ��݂�����ӗ��q�Ƃ���NEI[N]��ID��o�^����BNEI[N]�Ɋi�[�����ID�͂ǂ̂悤�ȏ��ɂȂ�H
									PART[i].distance.insert(pair<int, double>(j, dis)); //���ӗ��q�Ƃ̋���
									PART[i].set_current_neighboursID(j);	//�ł��d�v
									PART[i].set_current_neighbours_position(X, Y, Z);
									PART[i].set_current_distancebps(dis);	//���ꂢ��Ȃ��E�E�E
									N++;
									
								}
							}
						}
					}
				}	
			}
		}
		PART[i].PND=pnd;	//i�̉e�����a���Ɋ܂܂�闱�q�����x���擾
		PART[i].N=N;		//i�̉e�����a���Ɋ܂܂����ӗ��q�����擾

//			PART[i].set_density(N*mass/volume);

//			cout<<"mass: "<<mass<<", N: "<<N<<", volume: "<<volume<<endl;
//			cout<<"density_i: "<<PART[i].get_density()<<", compensated density: "<<ELAST.get_density()<<endl;

	}//�S���q�̕\�ʔ��肪�I��

//�ŏ����q�ԋ��������Ƃ߂�

	if(ELAST.get_symp_flag()==ON)
	{
		double min=ELAST.get_distancebp();//�ŏ����q�ԋ���
		int type1=0;
		int type2=0;
		double X1,Y1,Z1;

		//�����������Ă��ꂢ��Ȃ��E�E�E�R�����g�A�E�g
		//�S�Ă̗��q�𑍓���Ō������Ă���E�E�Ebsearch���g���Č���������
		for(int i=0;i<PART.size();i++)
		{
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
				if(dis<min)
				{
					type1=PART[i].type;
					type2=PART[j].type;
					min=dis;
					X1=PART[i].r[A_X];
					Y1=PART[i].r[A_Y];
					Z1=PART[i].r[A_Z];
				}
			}
		}//�ŏ����q�ԋ��������Ƃ܂���
		cout<<"mindis="<<min<<" between "<<type1<<" & "<<type2<<endl;
		//cout<<"(x,y)="<<X1<<" , "<<Y1<<endl;
		*mindis=min;
	}
}

void freeon(mpsconfig &CON, vector<mpselastic> &PART, double n0_4,int *INDEX,int **MESH,double *mindis, int t)
{
    //����:freeon�֐��ł͊e���q�̗��q�����x�����߂�ƂƂ��ɁA����ꂽ���x����\�ʔ�����s���B�܂��A�Œᗱ�q�ԋ��������łɂ��Ƃ߂Ă���
	//freeon2���x�����A���̂Ԃ���񉻂��e�ՁB
	double re=CON.get_re()*CON.get_distancebp();	//get_re()�̒l�����̂܂܎g���Ă̓_���I�I
	double mass=CON.get_particle_mass();
//	double volume=4*3.141592653*pow(re, 3)/3;
	int SIZE=CON.get_X_mesh()*CON.get_Y_mesh();
	
	//���x�����߁A�\�ʔ�����s���B�\�ʂȂ�P=0�ɂ���
	//omp_set_num_threads(8);//�X���b�h���w��

	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{
		PART[i].PND=0;
		PART[i].PND2=0;
		PART[i].N=0;
		PART[i].N2=0;
		PART[i].reset_current_distancebps();
		PART[i].reset_current_neighboursID();
		PART[i].reset_current_neighbours_position();
		PART[i].distance.clear();

		////���q�����x����
		double pnd=0;	//�S���q���l���ɂ��ꂽ���q�����x
		double pnd2=0;	//�e���̌v�Z�p���q�����x
		double pnd4=0;	//�\�ʔ���p
		double sum=0;
		int N=0;		//���ӗ��q��
			
		//��ԂɌŗL��ID(INDEX, MESH)�Ɨ��q�ɌŗL��ID(PART[i].index, PART[i].NEI[N])���ƍ����ĒT������
		for(int I=PART[i].index-1;I<=PART[i].index+1;I++)
		{
			if(PART[i].index>=SIZE && PART[i].index<CON.get_number_of_mesh()-SIZE)
			{       
				for(int J=-1*CON.get_X_mesh();J<=CON.get_X_mesh();J+=CON.get_X_mesh())
				{
					for(int K=-1*SIZE;K<=SIZE;K+=SIZE)
					{
						for(int L=0;L<INDEX[I+J+K];L++)
						{       
							int j=MESH[I+J+K][L];	//(I+J+K)�Ԗڂ̊i�q�Ɋ܂܂��L�Ԗڂ̗��q��ID���擾����
							if(j!=i)				//���qID�����ݒ��ڂ��Ă������(i)�ƈقȂ�΋������v�Z����
							{ 
								double X=PART[j].r[A_X]-PART[i].r[A_X];
								double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
								double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
								double dis=sqrt(X*X+Y*Y+Z*Z);
						
								if(dis<re)	
								{
									if(PART[j].type!=WALL){	//�Ǘ��q�͏������ӗ��q�ɉ����Ȃ�
									pnd+=kernel(re,dis);	//�v�Z�����������e�����a���ł���Ώd�݊֐����v�Z����
									PART[i].NEI[N]=j;		//i�̉e�����a�ɑ��݂�����ӗ��q�Ƃ���NEI[N]��ID��o�^����BNEI[N]�Ɋi�[�����ID�͂ǂ̂悤�ȏ��ɂȂ�H
									PART[i].distance.insert(pair<int, double>(j, dis)); //���ӗ��q�Ƃ̋���
									PART[i].set_current_neighboursID(j);	//�ł��d�v
									PART[i].set_current_neighbours_position(X, Y, Z);
									PART[i].set_current_distancebps(dis);	//���ꂢ��Ȃ��E�E�E
									N++;
									}
									else {
									pnd+=kernel(re,dis);	//�v�Z�����������e�����a���ł���Ώd�݊֐����v�Z����
									PART[i].NEI[N]=j;		//i�̉e�����a�ɑ��݂�����ӗ��q�Ƃ���NEI[N]��ID��o�^����BNEI[N]�Ɋi�[�����ID�͂ǂ̂悤�ȏ��ɂȂ�H
									PART[i].distance.insert(pair<int, double>(j, dis)); //���ӗ��q�Ƃ̋���
									PART[i].set_current_neighboursID(j);	//�ł��d�v
									PART[i].set_current_neighbours_position(X, Y, Z);
									PART[i].set_current_distancebps(dis);	//���ꂢ��Ȃ��E�E�E
									N++;		//�Ǘ��q�͏������ӗ��q�ɉ����Ȃ�
									}//*/
								}
							}
						}
					}
				}	
			}
		}
		PART[i].PND=pnd;	//i�̉e�����a���Ɋ܂܂�闱�q�����x���擾
		PART[i].N=N;		//i�̉e�����a���Ɋ܂܂����ӗ��q�����擾


	}//�S���q�̕\�ʔ��肪�I��

//�ŏ����q�ԋ��������Ƃ߂�

	double min=CON.get_distancebp();//�ŏ����q�ԋ���
	int type1=0,type2=0;
	double X1,Y1,Z1;

	//�����������Ă��ꂢ��Ȃ��E�E�E�R�����g�A�E�g
	//�S�Ă̗��q�𑍓���Ō������Ă���E�E�Ebsearch���g���Č���������
	for(int i=0;i<PART.size();i++)
	{
		for(int k=0;k<PART[i].N;k++)
		{
			int j=PART[i].NEI[k];
			double X=PART[j].r[A_X]-PART[i].r[A_X];
			double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
			double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
			double dis=sqrt(X*X+Y*Y+Z*Z);
			if(dis<min)
			{
				type1=PART[i].type;
				type2=PART[j].type;
				min=dis;
				X1=PART[i].r[A_X];
				Y1=PART[i].r[A_Y];
				Z1=PART[i].r[A_Z];
			}
		}
	}//�ŏ����q�ԋ��������Ƃ܂���
	cout<<"mindis="<<min<<" between "<<type1<<" & "<<type2<<endl;
	//cout<<"(x,y)="<<X1<<" , "<<Y1<<endl;
	*mindis=min;
	
}

///�\�ʔ���֐�ver.2
void freeon2(mpsconfig &CON,vector<mpselastic> &PART,int particle_number,double n0_4,int *INDEX,int **MESH,double *mindis,int fluid_number,int out)
{
	///freeon2�֐��̐���:flag1[i]�̓����ɂ�荂�����B�������A1CPU�Ȃ瑁�����A���CPU�ɂ����񉻂͌�����
	double le=CON.get_distancebp();
	int SIZE=CON.get_X_mesh()*CON.get_Y_mesh();
	double d=2;
	if(CON.get_dimension()==3) d=3;

	if(CON.get_T_field()==ON && CON.get_insulate()==1)
	{
		out=particle_number;	//��f�M�̉��x����v�Z����ۂ́AOUTWALL�����ӗ��q���i�[���Ă����K�v������B�����out=particle_number�Ƃ��邱�Ƃł���ɑΏ�
	}

	double *PND4=new double [out];//�\�ʔ���p���q�����x
	int *flag1=new int [out];		//�����t���O�B0:�������@1:�����ς�
	///������
	#pragma omp parallel for
	for(int i=0;i<out;i++)
	{
		PART[i].PND=0;//������
	    PART[i].PND2=0;
	    PART[i].N=0;
	    PART[i].N2=0;
	    PART[i].N3=0;

		PND4[i]=0;
		flag1[i]=0;
	}
	
	//���x�����߁A�\�ʔ�����s���B�\�ʂȂ�P=0�ɂ���
	for(int i=0;i<out;i++)//OUTWALL�ȊO��CFD���q�B//OUTWALL�̗��q�����x�Ȃǂ͂���Ȃ�
	{    
		////���q�����x����
		if(PART[i].index>=SIZE && PART[i].index<CON.get_number_of_mesh()-SIZE)
		{       
			for(int I=PART[i].index-1;I<=PART[i].index+1;I++)
			{       
				for(int J=-1*CON.get_X_mesh();J<=CON.get_X_mesh();J+=CON.get_X_mesh())
				{
					for(int K=-1*SIZE;K<=SIZE;K+=SIZE)
					{
						for(int L=0;L<INDEX[I+J+K];L++)
						{       
							int j=MESH[I+J+K][L];
							if(j<out)
							{
								if(flag1[j]==0 && j!=i)//�܂��������ĂȂ��Ȃ�
								{
									double X=PART[j].r[A_X]-PART[i].r[A_X];
									double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
									double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
									double dis=sqrt(X*X+Y*Y+Z*Z);
							
									if(dis<=CON.get_re()*le)//���z�E���U
									{       
										double r=CON.get_re()*le;
										double w=kernel(r,dis);
										PART[i].PND+=w;
										PART[j].PND+=w;
										PART[i].NEI[PART[i].N]=j;
										PART[j].NEI[PART[j].N]=i;
										PART[i].N++;
										PART[j].N++;
									}
									if(dis<=CON.get_re2()*le)
									{       
										double r=CON.get_re2()*le;
										double w=kernel(r,dis);
							
										PART[i].PND2+=w;
										PART[j].PND2+=w;
										
										PART[i].NEI2[PART[i].N2]=j;
										PART[j].NEI2[PART[j].N2]=i;
										PART[i].N2++;
										PART[j].N2++;
									}
									if(dis<=CON.get_re3()*le)//�\�ʒ���re3
									{       
										PART[i].NEI3[PART[i].N3]=j;
										PART[j].NEI3[PART[j].N3]=i;
										PART[i].N3++;
										PART[j].N3++;
									}
									if(dis<=CON.get_re4()*le)
									{       
									    double r=CON.get_re4()*le;
										double w=kernel(r,dis);
										PND4[i]+=w;
										PND4[j]+=w;
									}
								}
							}
							else
							{
								double X=PART[j].r[A_X]-PART[i].r[A_X];
								double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
								double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
								double dis=sqrt(X*X+Y*Y+Z*Z);
							
								if(dis<=CON.get_re()*le)//���z�E���U
								{       
									double r=CON.get_re()*le;
									double w=kernel(r,dis);
									PART[i].PND+=w;
									PART[i].NEI[PART[i].N]=j;
									PART[i].N++;
								}
								if(dis<=CON.get_re2()*le)
								{       
									double r=CON.get_re2()*le;
									double w=kernel(r,dis);
									PART[i].PND2+=w;
									PART[i].NEI2[PART[i].N2]=j;
									PART[i].N2++;
								}
								if(dis<=CON.get_re3()*le)//�\�ʒ���re3
								{       
									PART[i].NEI3[PART[i].N3]=j;
									PART[i].N3++;
								}
								if(dis<=CON.get_re4()*le)
								{       
								    double r=CON.get_re4()*le;
									double w=kernel(r,dis);
									PND4[i]+=w;
								}
							}
						}
					}
				}
			}
		}
		flag1[i]=1;//�����I��

		if(PART[i].N3>450) cout<<"error in freeon over. Change more than "<<PART[i].N3<<" coods:"<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
		////////////////////

		if(PART[i].type==ELASTIC)
//		if(PART[i].type==FLUID)
		{
			if(PND4[i]<n0_4*CON.get_beta())//���ȉ��Ȃ�
			{
				PART[i].surface=ON;//�\�ʗ��q�Ƃ���
				PART[i].P=0;//���E�����H
			}
			else PART[i].surface=OFF;
		}
		else if(PART[i].type==WALL)
//		else if(PART[i].type==INWALL)
		{
			if(PND4[i]<n0_4*CON.get_beta())//���ȉ��Ȃ�
			{
				PART[i].surface=ON;//�Ǖ\�ʗ��q�Ƃ���
				PART[i].P=0;
			}
			else  PART[i].surface=OFF;
		}
		else if(PART[i].type==SOLID)
//		else if(PART[i].type==INWALL)
		{
			if(PND4[i]<n0_4*CON.get_beta())//���ȉ��Ȃ�
			{
				PART[i].surface=ON;//�Ǖ\�ʗ��q�Ƃ���
				PART[i].P=0;
			}
			else  PART[i].surface=OFF;
		}
	}
	for(int i=out;i<particle_number;i++)
	{
		PART[i].N=0;
		PART[i].N2=0;
		PART[i].N3=0;
	}

	////�Œᗱ�q�ԋ��������Ƃ߂�
	double min0=CON.get_distancebp();//�Œᗱ�q�ԋ���
	int type1,type2,surface1,surface2;
	double X1,Y1,Z1;
	for(int i=0;i<fluid_number;i++)
	{
		for(int k=0;k<PART[i].N;k++)
		{
			int j=PART[i].NEI[k];
			double X=PART[j].r[A_X]-PART[i].r[A_X];
			double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
			double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
			double dis=sqrt(X*X+Y*Y+Z*Z);
			if(dis<min0)
			{
				type1=PART[i].type;
				surface1=PART[i].surface;
				type2=PART[j].type;
				surface2=PART[j].surface;
				min0=dis;
				X1=PART[i].r[A_X];
				Y1=PART[i].r[A_Y];
				Z1=PART[i].r[A_Z];
			}
		}
	}/////�Œᗱ�q�ԋ��������Ƃ܂���
	cout<<"mindis="<<min0<<" between "<<type1<<" "<<surface1<<" & "<<type2<<" "<<surface2<<endl;
	/////////*/
	*mindis=min0;

	delete [] PND4;
	delete [] flag1;
}

//�\�ʔ���֐�ver.2
void surface_judge2(mpsconfig &CON,vector<mpselastic> &PART,int fluid_number)
{
	//�]���̕\�ʔ���ɉ����āA�����ł���ɂӂ邢�ɂ�����B
	//���qi�̖@���x�N�g���ƁA���qj�����޸�قƂ̊p�x������l�𒴂��Ă�����\�ʂł͂Ȃ��Ɣ��f
	double *direct[DIMENSION];
    for(int D=0;D<DIMENSION;D++) direct[D]=new double [fluid_number];
		
	
    //////�@���޸�ٌv�Z
    for(int i=0;i<fluid_number;i++)
    {
        if(PART[i].surface==ON)  direct_f(CON,PART,i,direct);
        else  for(int D=0;D<3;D++) direct[D][i]=0;
	}

	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].surface==ON)//���̕\�ʗ��q���{���ɕ\�ʗ��q���肤�邩���肷��
		{
			int flag=ON;
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				if(j<fluid_number)
				{
					double X=PART[j].r[A_X]-PART[i].r[A_X];
					double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
					double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
					double dis=sqrt(X*X+Y*Y+Z*Z);
					double nx=X/dis;//i���q����j���q�ւƌ������P���޸��
					double ny=Y/dis;
					double nz=Z/dis;
					double inp=nx*direct[A_X][i]+ny*direct[A_Y][i]+nz*direct[A_Z][i];//�@���޸�قƒP���޸�قƂ̓���
					if(inp<-0.5) flag=OFF;//���ς�-0.5�̂��̂��ЂƂł�����Γ������̗��q�B���ς�-0.5�Ƃ������Ƃ́A���̊p�x��120�x�����e�͈͂Ƃ�������
				}
			}
			PART[i].surface=flag;
		}
	}

	for(int D=0;D<DIMENSION;D++) delete [] direct[D];

}

//freeon�֐� ver.3 ���q�����x�̂ݍČv�Z�B���q�\���q�֌W�͕ω����Ȃ��Ɖ��肵�Ă���
void freeon3(mpsconfig &CON,vector<mpselastic> &PART,int particle_number,int out)
{
	double d=2;
	if(CON.get_dimension()==3) d=3;
	double R=CON.get_re()*CON.get_distancebp();
	double R2=CON.get_re2()*CON.get_distancebp();
	for(int i=0;i<out;i++)
	{
		double W=0;
		for(int k=0;k<PART[i].N;k++)
		{
			int j=PART[i].NEI[k];
			double X=PART[j].r[A_X]-PART[i].r[A_X];
			double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
			double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
			double dis=sqrt(X*X+Y*Y+Z*Z);
			double w=kernel(R,dis);
			//double w=kernel2(R,dis,d);
			W+=w;
		}
		PART[i].PND=W;
	}
	if(CON.get_re()!=CON.get_re2())
	{
		for(int i=0;i<out;i++)
		{
			double W=0;
			for(int k=0;k<PART[i].N2;k++)
			{
				int j=PART[i].NEI2[k];
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
				double w=kernel(R2,dis);
				//double w=kernel2(R2,dis,d);
				W+=w;
			}
			PART[i].PND2=W;
		}
	}
	else for(int i=0;i<out;i++) PART[i].PND2=PART[i].PND;
}