#include "stdafx.h"
#include "Model.h"

#define FULL 1
#define HALF 2
#define HALFD 3
#define HALFD_shell 4

Model :: Model()
	:fq("initial_input.dat")
{
	
}

Model :: ~Model()
{
	fq.close();
}

////////////////////////���f���Z�b�g//////////////////////
int Model :: Model_set(){
	////���������ꂽmpsconfig�N���X�̃f�[�^
	mpsconfig CON;		
	////////////////////////////////////
	int number;
	if(CON.get_model_number()==0);
	else if(CON.get_model_number()==1 || CON.get_model_number()==11) number=Model :: Set_sphere_model() ;
	else if(CON.get_model_number()==2) number=Model :: Set_cylinder_model() ;
	else if(CON.get_model_number()==3) number=Model :: Set_cube_model() ;
	else if(CON.get_model_number()==4) number=Model :: Set_tensiontest_model();
	else if(CON.get_model_number()==5) number=Model :: Set_actuator_model();
	else if(CON.get_model_number()==6) number=Model :: Set_benchmark_model();
	else {cout<<"non-existent Model number"<<endl;
	exit(0);}
	return (number);
}

int Model :: Set_sphere_model()
{
	////���������ꂽmpsconfig�N���X�̃f�[�^
	mpsconfig CON;		
	////////////////////////////////////
	double le=CON.get_distancebp();
	int initial_ID=0;
	int number=0;
	double R=CON.get_fluidwidth()*le*0.5;			//�쐬����~�̔��a
	double Zg=CON.get_height();					//���̒��S����
	

			//���쐬
			int flag=FULL;
			number=Model :: Set_sphere_function(initial_ID,  le,  R);
		
		//������񏑂�����
		//for(int i=0;i<number;i++) writedata2(fq,i,X[i],Y[i],Z[i],FRFLUID,1,0,Y[i]*40,0,0,0,0);
		double vx=0;
		double vy=0;
		double vz=0;
		
			for(int i=0;i<number;i++)
			{
				double r=sqrt(PART[i].Get_X()*PART[i].Get_X()+PART[i].Get_Y()*PART[i].Get_Y());
				double a=0.2;
		//		a=500*3;//�N�[��������̂Ƃ�
				a=200;	//�V��
				
				
				vx=-0.5*a*PART[i].Get_X();
				vy=-0.5*a*PART[i].Get_Y();
				vz=a*PART[i].Get_Z();
				 
				writedata(fq,  i,PART[i].Get_X(),PART[i].Get_Y(),PART[i].Get_Z()+Zg,ELASTIC,1, OFF,0, vx,vy, vz,0, 0, 0);
		}
			return (number);
}

int Model :: Set_cylinder_model()
{
	return (0);
}

int Model :: Set_cube_model()
{
	mpsconfig CON;					//���������ꂽmpsconfig�N���X�̃f�[�^
	double le=CON.get_distancebp();
	double z=10.0;
	double y=25.0;
	double x=25.0;
	double zw=3.0;
	double yw=31.0;
	double xw=31.0;
	int num=0;

	for(int i=0;i<x;i++){
		for(int j=0;j<y;j++){
			for(int k=0;k<z;k++){
				PART.push_back(PART0);
				PART[num].Set(i*le, j*le, k*le, le);
				num++;
			}
		}
	}
	int elasnum=num;  //1�`�̐� �Ǘ��qID��1�Ԗ�

/*	for(int i=0;i<xw;i++){
		for(int j=0;j<yw;j++){
			for(int k=0;k<zw;k++){
				PART.push_back(PART0);
				PART[num].Set(i*le, j*le, k*le, le);
				num++;
			}
		}
	}*/
	for(int i=0;i<elasnum;i++) writedata(fq,i,PART[i].Get_X()-((x-1)*le)/2.0,PART[i].Get_Y()-((y-1)*le)/2.0,PART[i].Get_Z(),MAGELAST,1,0,0,0,0,0,0,0,1);//���q��MAGELAST
//	for(int i=elasnum;i<num;i++) writedata(fq,i,PART[i].Get_X()-(xw-1)*le/2.0,PART[i].Get_Y()-(yw-1)*le/2.0,PART[i].Get_Z()-(2.0+zw)*le,WALL,1,0,0,0,0,0,0,0,0);//���q��WALL
	return num;
}
int Model :: Set_tensiontest_model()
{
	mpsconfig CON;					//���������ꂽmpsconfig�N���X�̃f�[�^
	int length=25;	//[mm]��������[�q�Ƃ͕�35
	int width=11;//15
	int depth=11;
	int tan=0;//�[�q�������q5
	int number=0;//�������q��
	double accel=0.005;
	double le=CON.get_distancebp();
	le*=0.5;//*(11.0/12.0)

	int low_ID=(((width*depth)-1)/2)+(width*depth*(8-1));
	int left_ID=(width*depth)*(length-1)/2+(width*2)+(width-1)/2;
	int right_ID=(width*depth)*(length-1)/2+(width*(width-3))+(width-1)/2;
	int top_ID=(((width*depth)-1)/2)+(width*depth*(18-1));
	Model :: Set_point(low_ID,top_ID,left_ID,right_ID);
	cout<<low_ID<<","<<top_ID<<","<<left_ID<<","<<right_ID<<endl;
	////////////���[�q////////////
	for(int i=0;i<tan;i++){
		for(int j=0;j<width;j++){
			for(int k=0;k<depth;k++){
				PART.push_back(PART0);
				PART[number].Set(k*le, j*le, i*le, le);
				number++;
			}
		}
	}
	int low=number;
	//////////////////////////////
	/////////////MRE����////////////
	for(int i=tan;i<length+tan;i++){
		for(int j=0;j<width;j++){
			for(int k=0;k<depth;k++){
				PART.push_back(PART0);
				PART[number].Set(k*le, j*le, i*le, le);
				number++;
			}
		}
	}
	int mid=number;
	/////////////////////////////////
	//////////////��[�q/////////////
	for(int i=length+tan;i<length+tan*2;i++){
		for(int j=0;j<width;j++){
			for(int k=0;k<depth;k++){
				PART.push_back(PART0);
				PART[number].Set(k*le, j*le, i*le, le);
				number++;
			}
		}
	}
	////////////////////////////////
		
	for(int i=0;i<number;i++){ 
		if(i<low)writedata(fq,i,PART[i].Get_X()-depth*le/2,PART[i].Get_Y()-width*le/2,PART[i].Get_Z()-(length+10)*le/2,WALL,0, OFF,0, 0, 0, -accel, 0, 0, 0);
		else if(i>=low && i<mid)writedata(fq,i,PART[i].Get_X()-depth*le/2,PART[i].Get_Y()-width*le/2,PART[i].Get_Z()-(length+10)*le/2,MAGELAST,0, OFF,0, 0,0, 0, 0, 0, 0);
		else if(i>=mid)writedata(fq,i,PART[i].Get_X()-depth*le/2,PART[i].Get_Y()-width*le/2,PART[i].Get_Z()-(length+10)*le/2,WALL,0, OFF,0, 0, 0, accel, 0, 0, 0);
	}
	return number;
}

///////////////�A�N�`���G�[�^���f��/////////////////
int Model :: Set_actuator_model()	
{
	//////�A�N�`���G�[�^�̐��@///////
	double elast_R=0.012;	//���a
	double elast_H=0.03;
	double MRE_R=0.015;
	double MRE_H=0.04;
	////////////////////////////////

	////���������ꂽmpsconfig�N���X�̃f�[�^
	mpsconfig CON;		
	////////////////////////////////////

	//////���q�̒��a////
	double le=CON.get_distancebp();
	////////////////////

	//////���q���i�[//////
	int total_number;	//���̊֐��ł̑S���q��
	int start_ID=0;		//�������q��
	/////////////////////

	//////�~���`�쐬/////
	total_number=Model :: Set_cylinder(start_ID, le, MRE_R, MRE_H);
	////////////////////

	for(int i=start_ID;i<total_number;i++)
		{
			if(pow(PART[i].Get_X(),2)+pow(PART[i].Get_Y(),2)<=pow(0.01,2) && PART[i].Get_Z()>0.01 && PART[i].Get_Z()<0.02+0.01){
					writedata(fq,i,PART[i].Get_X(),PART[i].Get_Y(),PART[i].Get_Z()-(MRE_H/2),WALL,1,0,0,0,0,0,0,0,0);//���q��WALL
			}
			else if(pow(PART[i].Get_X(),2)+pow(PART[i].Get_Y(),2)<=pow(elast_R,2) && PART[i].Get_Z()>0.005 && PART[i].Get_Z()<elast_H+0.0045){
				writedata(fq,i,PART[i].Get_X(),PART[i].Get_Y(),PART[i].Get_Z()-(MRE_H/2),ELASTIC,1,0,0,0,0,0,0,0,0);//���q��ELASTIC
			}
			else writedata(fq,i,PART[i].Get_X(),PART[i].Get_Y(),PART[i].Get_Z()-(MRE_H/2),MAGELAST,1,0,0,0,0,0,0,0,1);//���q��MAGELAST
			
		}

	return total_number;

}

int Model :: Set_benchmark_model()
{
	double bord_r=0.015;
	double bord_h=0.005;
	double elast_r=0.012;
	double elast_h=0.010;
	int total_number;
	int start_ID=0;
	mpsconfig CON;					//���������ꂽmpsconfig�N���X�̃f�[�^
	double le=CON.get_distancebp();

	total_number=Model :: Set_cylinder(start_ID, le, bord_r, bord_h);
	for(int i=start_ID;i<total_number;i++) writedata(fq,i,PART[i].Get_X(),PART[i].Get_Y(),PART[i].Get_Z()-(elast_h/2+bord_h),WALL,1,0,0,0,0,0,0,0,0);//���q��WALL

	start_ID=total_number;
	total_number+=Model :: Set_cylinder(start_ID, le, elast_r, elast_h);
	for(int j=start_ID;j<total_number;j++) writedata(fq,j,PART[j].Get_X(),PART[j].Get_Y(),PART[j].Get_Z()-(elast_h/2),ELASTIC,1,0,0,0,0,0,0,0,0);//���q��ELASTIC

/*	start_ID=total_number;
	total_number+=Model :: Set_cylinder(start_ID, le, bord_r, bord_h);
	for(int k=start_ID;k<total_number;k++) writedata(fq,k,PART[k].Get_X(),PART[k].Get_Y(),PART[k].Get_Z()+(elast_h/2),WALL,1,0,0,0,0,0,0,0,0);//���q��WALL
	*/
	return total_number;
}

///////////////���q�̏����O���t�@�C���ɏo��/////////////////
void Model :: writedata(ofstream &fp, int id, double x, double y,double z, int type,int materialID,int surface,double val,double vx,double vy,double vz,double P,double h,int toBEM)
{
	//�t�@�C���o��
	fp<<id<<"\t"<<x<<"\t"<<y<<"\t"<<z<<"\t"<<vx<<"\t"<<vy<<"\t"<<vz<<"\t"<<P<<"\t"<<h<<"\t"<<val<<"\t"<<type<<"\t"<<materialID<<"\t"<<surface<<"\t"<<toBEM<<"\t"<<endl;
}

int Model :: Set_sphere_function(int initial_ID, double le, double R)
{
	int number=0;
	int flag=FULL;
	cout<<"Set_circle_edge�J�n"<<endl;
	number+=Set_circle_edge(initial_ID, le, R);
	cout<<"Set_circle_in_using_6_pieces�J�n"<<endl;
	number+=Set_circle_in_using_6_pieces(number,le, R, initial_ID);
	cout<<"Set_sphere�J�n"<<endl;
	number+=Set_sphere(number,le, R, flag);
	return number;
}

int Model :: Set_cylinder(int initial_ID, double le, double R, double height)	//�~���쐬�֐�
{
	int number=0;
	cout<<"Set_circle_edge�J�n"<<endl;
	number+=Set_circle_edge(initial_ID, le, R);
	cout<<"Set_circle_in�J�n"<<endl;
	number+=Set_circle_in(number, le, R, initial_ID);
	cout<<"Set_cylinder_face�J�n"<<endl;
	number+=Set_cylinder_face(number, le, R, height, initial_ID);
	cout<<"Set_cylinder_in�J�n"<<endl;
	number+=Set_cylinder_in(number, le, R, height, initial_ID);
	cout<<"�I��"<<endl;
	return number;
}

int Model :: Set_circle_edge(int initial_number, double le, double R)	//�~�O���֐�
{
	int N=(int)((2*PI*R)/le);//�~���̕�����
	
	double L=2*PI*R/N;				//���q�ԋ���
	double theta=L/R;				//���q��z�u����p�x

	for(int n=initial_number;n<N;n++)
	{
		PART.push_back(PART0);
		PART[n].Set(R*cos(theta*n), R*sin(theta*n), 0, le);
	}
	return N;
}

int Model :: Set_circle_in(int number, double le, double R, int initial_number)	//�~�����[�U�֐�
{
	//set_circle_in()�ƈႢ�A60�x�������v�Z���A�����6�º�߰���ĉ~���\������B���ԒZ�k�Ɣz�u�̋ϓ������ړI
	//edge_startID����(edge_lastID-1)�܂ł̗��q���A�~�̊O�����\�����闱�q�ɊY������
	//vector�^�z��͎Q�Ɠn�����Ă���Bvector<double> *X�ł͂Ȃ�vector<double> &X�ł��邱�Ƃɒ��ӁB����Ŋe�z��͒ʏ�ʂ�Ɏd�l�\�B�A���[���Z�q������Ȃ�
	//�Q�Ɠn���łȂ��ʏ�̂�肩���ł��������ǁA���̏ꍇ�A�Ⴆ��a=X[5]�Ə����Ă����p�ł��Ȃ��Ba=(*X)[5]�ȂǂƂ��Ȃ���΂Ȃ�Ȃ�

	int newN=0;					//�̊֐��ŐV�����ǉ����闱�q��
	int beforeN=number;		//���̊֐��Ăяo�����ɂ����闱�q��
	

	double A=sqrt(3.0)/2;		//�悭�g���W��

	int half_WX=(int)(R/le)+1;  //�~���\���܂ގl�p�`��z�肷��B���̎l�p�`�̕�*0.5
	int half_WY=(int)(R/(le*A))+1;  //�~���\���܂ސ��l�p�`��z�肷��B���̎l�p�`�̕�*0.5
	double R2=R-le*0.5;				//���������߂̔��a��ݒ�

	//////////////////////�~��6�ɂ킯�邽�߂̒���6�𐶐�

	double temp_R_num=R/le;			//���a�����ɐݒu����w���́x���q��
	int    R_num=(int)temp_R_num;	//�^�̗��q���@�Ƃ肠�������̗��q���̐����ԂƂ���B�����ŁAtemp_R_num>R_num���������Ă���B
	double difference=temp_R_num-R_num;	//���̐��Ɛ^�̐��̍�
	
	//���̐��Ɛ^�̐��̍���0.5�܂łȂ�A�^�̐���N�Ƃ���B0.5�ȏ�Ȃ�N+1�Ƃ���
	if(difference>0.5) R_num++;
	double L=R/R_num;					//���q�ԋ���

	PART.push_back(PART0);
	PART[beforeN].Set(0, 0, 0, le);						//���S���q��ǉ�
	
	newN++;
	for(int k=0;k<6;k++)				//6�̒�����loop
	{
		double theta=PI/3*k;			//�����̊p�x
		for(int n=1;n<R_num;n++)		//���S���q�ƍŊO�����q�͂������邩��A�����ł�loop�͂�����J�E���g���Ȃ�
		{
			double r=L*n;				//���S����̋���
			PART.push_back(PART0);
			PART[beforeN+newN].Set(r*cos(theta), r*sin(theta), 0, le);
			newN++;
		}
	}	
	beforeN+=newN;	//�����܂ł̗��q��
	newN=0;			//���ꂩ�瑝���闱�q��
	/////////////////����6�𐶐�����
	//���������ʒu �������ŏ���1�s�[�X�̂�
	for(int i=0;i<=half_WX;i++)
	{
		for(int j=1;j<=half_WY;j++)
		{
			double jj=j*le*A;
			double ii=i*le;
			if(j%2!=0) ii+=0.5*le;//j����Ȃ�ii��0.5�i�q�������炷
			if(ii*ii+jj*jj<R2*R2)
			{
				if(jj<sqrt(3.0)*ii-le)		//�ŏ��̃s�[�X�̎΂ߐ����Ⴂ�̈�ɐݒu�B�������������肬��͂܂����̂ŁA�ی���-le���Ă���
				{
					PART.push_back(PART0);
					PART[beforeN+newN].Set(ii, jj, 0, le);	
					newN++;
				}
			}
		}
	}////////////////////////
	//���q���͊w�ɂ��ʒu���œK��
	//MD_2D(X,Y,Z,le,0,beforeN,beforeN,newN);
	MD_2D(le,initial_number,beforeN,newN);
	int M_N=newN;
	///�������q����������6�º�߰
	for(int angle=1;angle<6;angle++)
	{
		double theta=PI/3*angle;//��]����p�x
		for(int k=0;k<M_N;k++)
		{
			int i=beforeN+k;
			double x=cos(theta)*PART[i].Get_X()-sin(theta)*PART[i].Get_Y();//��]��̍��W
			double y=sin(theta)*PART[i].Get_X()+cos(theta)*PART[i].Get_Y();

			PART.push_back(PART0);
			PART[beforeN+newN].Set(x, y, 0, le);
			newN++;
		}
	}///////////////////*/
	return (beforeN+newN-number);	//�Œ�g���q���{�������q��
}

int Model:: Set_circle_in_using_6_pieces(int number,double le,double R,int  initial_ID)
{
	//set_circle_in()�ƈႢ�A60�x�������v�Z���A�����6�º�߰���ĉ~���\������B���ԒZ�k�Ɣz�u�̋ϓ������ړI
	//edge_startID����(edge_lastID-1)�܂ł̗��q���A�~�̊O�����\�����闱�q�ɊY������
	//vector�^�z��͎Q�Ɠn�����Ă���Bvector<double> *X�ł͂Ȃ�vector<double> &X�ł��邱�Ƃɒ��ӁB����Ŋe�z��͒ʏ�ʂ�Ɏd�l�\�B�A���[���Z�q������Ȃ�
	//�Q�Ɠn���łȂ��ʏ�̂�肩���ł��������ǁA���̏ꍇ�A�Ⴆ��a=X[5]�Ə����Ă����p�ł��Ȃ��Ba=(*X)[5]�ȂǂƂ��Ȃ���΂Ȃ�Ȃ�

	int newN=0;					//�̊֐��ŐV�����ǉ����闱�q��
	int beforeN=number;		//���̊֐��Ăяo�����ɂ����闱�q��

	double A=sqrt(3.0)/2;		//�悭�g���W��

	int half_WX=(int)(R/le)+1;  //�~���\���܂ގl�p�`��z�肷��B���̎l�p�`�̕�*0.5
	int half_WY=(int)(R/(le*A))+1;  //�~���\���܂ސ��l�p�`��z�肷��B���̎l�p�`�̕�*0.5
	double R2=R-le*0.5;				//���������߂̔��a��ݒ�


	//////////////////////�~��6�ɂ킯�邽�߂̒���6�𐶐�

	double temp_R_num=R/le;			//���a�����ɐݒu����w���́x���q��
	int    R_num=(int)temp_R_num;	//�^�̗��q���@�Ƃ肠�������̗��q���̐����ԂƂ���B�����ŁAtemp_R_num>R_num���������Ă���B
	double difference=temp_R_num-R_num;	//���̐��Ɛ^�̐��̍�
	
	//���̐��Ɛ^�̐��̍���0.5�܂łȂ�A�^�̐���N�Ƃ���B0.5�ȏ�Ȃ�N+1�Ƃ���
	if(difference>0.5) R_num++;
	double L=R/R_num;					//���q�ԋ���
	PART.push_back(PART0);
	PART[newN].Set(0,0,0,le);	//���S���q��ǉ�			
	newN++;

	for(int k=0;k<6;k++)				//6�̒�����loop
	{
		double theta=PI/3*k;			//�����̊p�x
		for(int n=1;n<R_num;n++)		//���S���q�ƍŊO�����q�͂������邩��A�����ł�loop�͂�����J�E���g���Ȃ�
		{
			double r=L*n;				//���S����̋���
			PART.push_back(PART0);
			PART[newN].Set(r*cos(theta),r*sin(theta),0,le);	//���S���q��ǉ�			
			newN++;
		}
	}
	number=number+newN;
	newN=0;
	beforeN=number;
	/////////////////����6�𐶐�����

	//���������ʒu �������ŏ���1�s�[�X�̂�
	for(int i=0;i<=half_WX;i++)
	{
		for(int j=1;j<=half_WY;j++)
		{
			double jj=j*le*A;
			double ii=i*le;
			if(j%2!=0) ii+=0.5*le;//j����Ȃ�ii��0.5�i�q�������炷
			if(ii*ii+jj*jj<R2*R2)
			{
				if(jj<sqrt(3.0)*ii-le)		//�ŏ��̃s�[�X�̎΂ߐ����Ⴂ�̈�ɐݒu�B�������������肬��͂܂����̂ŁA�ی���-le���Ă���
				{
					PART.push_back(PART0);
					PART[newN].Set(ii,jj,0,le);	//���S���q��ǉ�			
					newN++;
				}
			}
		}
	}////////////////////////

	//���q���͊w�ɂ��ʒu���œK��
	//MD_2D(X,Y,Z,le,0,beforeN,beforeN,newN);
	MD_2D(le, initial_ID,beforeN,newN);

	///�������q����������6�º�߰
	for(int angle=1;angle<6;angle++)
	{
		double theta=PI/3*angle;//��]����p�x
		for(int k=0;k<newN;k++)
		{
			int i=beforeN+k;
			double x=cos(theta)*PART[i].Get_X()-sin(theta)*PART[i].Get_Y();//��]��̍��W
			double y=sin(theta)*PART[i].Get_X()+cos(theta)*PART[i].Get_Y();

			PART.push_back(PART0);
			PART[newN].Set(x,y,0,le);	//���S���q��ǉ�			
			newN++;
			}
	}///////////////////*/

	return (number+newN*6);//newN�͂ЂƂ̃s�[�X���̗��q����\���Ă��邩�炱���ł�6�{
}


//���aR�̋��쐬�֐�
int Model :: Set_sphere(int number,double le,double R,int flag)
{
	//�܂��͔��������B���̂��߂ɂ͔����\�ʂ��쐬����K�v������B
	int newN=0;					//�̊֐��ŐV�����ǉ����闱�q��
	int beforeN=number;		//���̊֐��Ăяo�����ɂ����闱�q��
	int beforeN2=number;		//�֐��Ăяo�����̗��q���B���̊֐��̍Ō�܂ŋL�����Ă���

	double A=sqrt(3.0)/2;				//�悭�g���W��
	double B=sqrt(2.0/3);						////�悭�g���W��
	int half_WX=(int)(R/le)+1;		 //�����\���܂ގl�p�`��z�肷��B���̎l�p�`�̕�*0.5
	int half_WY=(int)(R/(le*A))+1;  //�����\���܂ސ��l�p�`��z�肷��B���̎l�p�`�̕�*0.5
	int half_WZ=(int)(R/(le*B))+1;  //�����\���܂ސ��l�p�`��z�肷��B���̎l�p�`�̕�*0.5
	double R2=R-0.5*le;				//���������߂̔��a��ݒ�

	///////////�����\��
	int Nt;						//���\�ʂ́A�ƕ����̕�����
	double Lt;					//���\�ʂ́A�ƕ����̕�������
	Model::Set_calc_N_and_L(PI/2*R,le*A,&Nt,&Lt);//���~�̕������͋����E��ǂ���ł��悢
	double d_theta=Lt/R;		//�ʂ̒�����Lt�ɂȂ�p�x

	for(int k=0;k<Nt;k++)//loop��k<Nt�ŏI��点��BNt�ɊY������Ƃ���͂��łɐݒu�ς�
	{
		double THETA=k*d_theta;	//��
		double r=R*sin(THETA);	//���̍����ɂ�����~�̔��a
		double round=2*PI*r;//���̍����ɂ�����~��

		int Nf=Model :: Set_calc_division_N_circle(round,le);//���\�ʂ́A�ƕ����̕�����
		double Lf=round/Nf;						//���\�ʂ́A�ƕ����̕�������
		double d_fai=Lf/r;						//�ʂ̒�����Lf�ɂȂ�p�x
		
		for(int i=0;i<Nf;i++)
		{
			double fai=d_fai*i;
			if(Nt%2==0)
			{
				if(k%2!=0) fai+=0.5*d_fai;//Nt�������Ȃ�A�쐬�ς݂̉~�Ɛڂ���Ƃ��͊�ԖځB����Ċ�����炷
			}
			else
			{
				if(k%2==0) fai+=0.5*d_fai;//Nt����Ȃ�A�쐬�ς݂̉~�Ɛڂ���Ƃ��͋����ԖځB����Ċ�����炷
			}
			double x=r*cos(fai);
			double y=r*sin(fai);
			double z=R*cos(THETA);
			PART.push_back(PART0);
			PART[newN].Set(x,y,z,le);			
			newN++;
		}
	}
	if(Nt%2!=0)//Nt����̂Ƃ��́A����ɗ��q���u����Ȃ���΂Ȃ�Ȃ��B���������loop�͂��ꂪ�s�\�B����Ă����Œǉ�
	{
		PART.push_back(PART0);
		PART[newN].Set(0,0,R,le);			
		newN++;
	}
	//////////////////////////////////////

	number=number+newN;

	if(flag!=HALFD_shell)
	{
	newN=0;					//�̊֐��ŐV�����ǉ����闱�q��
	beforeN=number;		//���̊֐��Ăяo�����ɂ����闱�q��
	
	//���������ʒu
	for(int i=-half_WX;i<=half_WX;i++)
	{
		for(int j=-half_WY;j<=half_WY;j++)
		{
			for(int k=1;k<half_WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//j����Ȃ�ii��0.5�i�q�������炷
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//k����Ȃ�ii��jj�����炷
				if(ii*ii+jj*jj+kk*kk<R2*R2*0.6*0.6)//���a*0.7���̗��q�͈ʒu�Œ�
				{
					PART.push_back(PART0);
					PART[newN].Set(ii,jj,kk,le);			
					newN++;
				}
			}
		}
	}
	number=number+newN;

	newN=0;					//�̊֐��ŐV�����ǉ����闱�q��
	beforeN=number;		//���̊֐��Ăяo�����ɂ����闱�q��
	///////////////////////*/

	

	//���������ʒu
	for(int i=-half_WX;i<=half_WX;i++)
	{
		for(int j=-half_WY;j<=half_WY;j++)
		{
			for(int k=1;k<half_WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//j����Ȃ�ii��0.5�i�q�������炷
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//k����Ȃ�ii��jj�����炷
				//if(ii*ii+jj*jj+kk*kk<R2*R2)
				if(ii*ii+jj*jj+kk*kk<R2*R2 && ii*ii+jj*jj+kk*kk>=R2*R2*0.6*0.6)
				{
					PART.push_back(PART0);
					PART[newN].Set(ii,jj,kk,le);			
					newN++;
				}
			}
		}
	}///////////////////////*/
	

	//���q���͊w�ɂ��ʒu���œK��
	double r=1.5;
	double rigion[3][2];	//��͗̈�
	rigion[A_X][0]=-1.2*R; rigion[A_X][1]=1.2*R;
	rigion[A_Y][0]=-1.2*R; rigion[A_Y][1]=1.2*R;
	rigion[A_Z][0]=-0.2*R; rigion[A_Z][1]=1.2*R;

	MD_3D(le,beforeN2,beforeN,newN,r,rigion);

	number=number+newN;
	}

	///�㔼�����������ֺ�߰����������
	if(flag==FULL)
	{
		newN=0;					//�V�����ǉ����闱�q��
		beforeN=number;		//���̎��ɂ����闱�q��

		for(int k=0;k<beforeN;k++)
		{
			if(PART[k].Get_Z()>0.4*le)
			{
				newN++;
				PART.push_back(PART0);
				PART[newN].Set(PART[k].Get_X(),PART[k].Get_Y(),-PART[k].Get_Z(),le);			
				newN++;
			}
		}
		number=number+newN;
	}///////////////////////////////
	if(flag==HALFD || flag==HALFD_shell)		//���������ق����Ƃ��ɁA�������㔼�����㉺���]������
	{
		for(int k=beforeN2;k<number;k++) PART[k].Multiply(1,1,-1);
	}
	return (number);
}

void Model :: Set_point(int low_ID,int top_ID,int lefy_ID,int right_ID)
{
	
	point[0]=low_ID;
	point[1]=top_ID;
	point[2]=lefy_ID;
	point[3]=right_ID;
	
}

void Model :: Get_point(int &a,int &b,int &c,int &d)
{
	a=point[0];
	b=point[1];
	c=point[2];
	d=point[3];
}

void Model :: MD_2D(double le,int BstartID,int beforeN,int newN)
{
	//���q���͊w�ɂ��newN�̗��q�̈ʒu���œK���@ID��BstartID����BendID�܂ł̂͋��E���q�Ȃ̂œ������Ȃ�

	double region[2][2];	//��͗̈�

	/////////////////////��͗̈�̌���
	region[A_X][0]=100; region[A_X][1]=-100;
	region[A_Y][0]=100; region[A_Y][1]=-100;
	for(int i=BstartID;i<beforeN;i++)
	{
		if(PART[i].Get_X()<region[A_X][0]) region[A_X][0]=PART[i].Get_X();
		else if(PART[i].Get_X()>region[A_X][1]) region[A_X][1]=PART[i].Get_X();

		if(PART[i].Get_Y()<region[A_Y][0]) region[A_Y][0]=PART[i].Get_Y();
		else if(PART[i].Get_Y()>region[A_Y][1]) region[A_Y][1]=PART[i].Get_Y();
	}
	for(int D=0;D<2;D++)
	{
		region[D][0]-=5*le;	//�����̈���L�߂ɂƂ�
		region[D][1]+=5*le;
	}//////////////////////////

	//�p�����[�^
	double k0=1;
	double r=1.5;
	double dt=0.001;
	
	//�͂�ax^3+bx^2+d�̎����̗p�B����[Bubble Mesh Automated Triangular Meshing of Non-Manifold Geometry by Sphere Packing]���Q��
	double a=(r+1)/(2*r*r-r-1)*k0/(le*le);
	double b=-0.5*k0/le-1.5*a*le;
	double d=-a*le*le*le-b*le*le;
	/////////////
	int lastN=beforeN+newN;
	vector<double> Fx(newN);	//�e���q�ɓ���X������
	vector<double> Fy(newN);	//�e���q�ɓ���Y������
	vector<double> Ax(newN,0);	//X���������x
	vector<double> Ay(newN,0);	//Y���������x
	vector<double> U(newN,0);	//X�������x
	vector<double> V(newN,0);	//Y�������x
	vector<double> visX(newN);	//X�����S���W��
	vector<double> visY(newN);	//Y�����S���W��
	//�v�Z�̍������̂��߂Ɋi�q���`�� ��͕���r*le�Ŋ���؂��Ƃ͌���Ȃ��̂ŁA�͂ݏo�����Ƃ���͐؂�̂āB�Ȃ̂Ŋe���Ƃ����̕����ɂ͗]�T��������
	double grid_width=le*((int)(r+1));								//�i�q�̕��Br���܂ސ���*le
	int grid_sizeX=(int)((region[A_X][1]-region[A_X][0])/grid_width);	//X�����̊i�q�̌�
	int grid_sizeY=(int)((region[A_Y][1]-region[A_Y][0])/grid_width);
	int plane_SIZE=grid_sizeX*grid_sizeY;
	int *index=new int[newN];									//�e�������q���܂ފi�q�ԍ�
//	cout<<"�������ۂ�"<<endl;
//	vector<int> *MESH=new vector<int>[plane_SIZE];				//�e���b�V���Ɋi�[����闱�qID�i�[
	vector<vector<int> > MESH;
	MESH.resize(plane_SIZE);
	for(int i=BstartID;i<beforeN;i++)	//�܂��͋��E���q���i�q�Ɋi�[
	{
		int xn=(int)((PART[i].Get_X()-region[A_X][0])/grid_width);	//X�����ɉ��ڂ̊i�q�� 
		int yn=(int)((PART[i].Get_Y()-region[A_Y][0])/grid_width);	//Y�����ɉ��ڂ̊i�q��
		int number=yn*grid_sizeX+xn;					//���qi���܂ފi�q�̔ԍ�
		MESH[number].push_back(i);
	}
	for(int k=0;k<newN;k++)	//���ɓ������q���i�[
	{
		int i=beforeN+k;
		int xn=(int)((PART[i].Get_X()-region[A_X][0])/grid_width);//X�����ɉ��ڂ̊i�q�� 
		int yn=(int)((PART[i].Get_Y()-region[A_Y][0])/grid_width);//Y�����ɉ��ڂ̊i�q��
		int number=yn*grid_sizeX+xn;					//���qi���܂ފi�q�̔ԍ�
		MESH[number].push_back(i);
		index[k]=number;
	}//////////////////////////////////////////

	//�v�Z�J�n
	for(int t=0;t<100;t++)
	{
		if(t%10==0 &&t>0)//MESH����蒼��
		{
			//�܂���MESH����x�j�󂷂�B
			for(int n=0;n<plane_SIZE;n++)
			{
				size_t size=MESH[n].size();
				for(int k=0;k<size;k++) MESH[n].pop_back();
			}
			
			for(int i=BstartID;i<beforeN;i++)	//�܂��͋��E���q���i�q�Ɋi�[
			{
				int xn=(int)((PART[i].Get_X()-region[A_X][0])/grid_width);	//X�����ɉ��ڂ̊i�q�� 
				int yn=(int)((PART[i].Get_Y()-region[A_Y][0])/grid_width);	//Y�����ɉ��ڂ̊i�q��
				int number=yn*grid_sizeX+xn;					//���qi���܂ފi�q�̔ԍ�
				MESH[number].push_back(i);
			}
			for(int k=0;k<newN;k++)	//���ɓ������q���i�[
			{
				int i=beforeN+k;
				int xn=(int)((PART[i].Get_X()-region[A_X][0])/grid_width);//X�����ɉ��ڂ̊i�q�� 
				int yn=(int)((PART[i].Get_Y()-region[A_Y][0])/grid_width);//Y�����ɉ��ڂ̊i�q��
				int number=yn*grid_sizeX+xn;					//���qi���܂ފi�q�̔ԍ�
				MESH[number].push_back(i);
				index[k]=number;
			}
		}////////////

		for(int k=0;k<newN;k++)
		{
			Fx[k]=0; Fy[k]=0;					//������
			int i=beforeN+k;					//�Ή����闱�q�ԍ�
			double kx=0;						//X�����o�l�W��
			double ky=0;
			int G_id=index[k];				//�i�[����i�q�ԍ�
			for(int II=G_id-1;II<=G_id+1;II++)
			{       
				for(int JJ=-1*grid_sizeX;JJ<=grid_sizeX;JJ+=grid_sizeX)
				{
					int M_id=II+JJ;
					for(int L=0;L<MESH[M_id].size();L++)
					{
						int j=MESH[M_id][L];
						double x=PART[j].Get_X()-PART[i].Get_X();
						double y=PART[j].Get_Y()-PART[i].Get_Y();
						double dis=sqrt(x*x+y*y);
						if(dis<r*le && dis!=0)			//����loop�͎������g���ʉ߂��邩��Adis!=0�͕K�v
						{
							double F=a*dis*dis*dis+b*dis*dis+d;
							Fx[k]-=F*x/dis;					//F�̒l�����̂Ƃ��͐˗͂Ȃ̂ŁA-=�ɂ���
							Fy[k]-=F*y/dis;
							double K=3*a*dis*dis+2*b*dis;//�o�l�W���@�͂̎��̔����ɑ���
							K=sqrt(K*K);					//���̒l���~�����B�����畉�̂Ƃ��ɔ����Đ��ɕϊ�
							kx+=K*x*x/(dis*dis);			//k���e�����ɕ��z�B�����ŁA��ɐ��̗ʂ����z�����悤��x*x/(dis*dis)�ƂȂ��Ă���
							ky+=K*y*y/(dis*dis);
						}
					}
				}
			}
			visX[k]=1.414*sqrt(kx);//���̂悤�Ɋe�������̔S���W�������߂�B�����u�������f���ɂ�鎩�����b�V�������vP6�Q�ƁB���������ʂ�1�Ƃ��Ă���B
			visY[k]=1.414*sqrt(ky);
			Ax[k]=(Fx[k]-visX[k]*U[k]);
			Ay[k]=(Fy[k]-visY[k]*V[k]);
		}//�e���q�̉����x�����܂����B
		
		if(t==0)	//�ŏ��̽ï�ߎ���dt������
		{
			double MaxAccel=0;
			for(int k=0;k<newN;k++)
			{
				double accel2=Ax[k]*Ax[k]+Ay[k]*Ay[k];
				if(accel2>MaxAccel) MaxAccel=accel2;
			}
			MaxAccel=sqrt(MaxAccel);//�ő�����x�����܂���
			dt=sqrt(0.02*le/MaxAccel);
		}

		for(int k=0;k<newN;k++)//���x�ƈʒu�̍X�V
		{
			int i=beforeN+k;
			double u=U[k];
			double v=V[k];
			U[k]+=dt*Ax[k];
			V[k]+=dt*Ay[k];
			PART[i].Add(dt*(U[k]+u)*0.5, dt*(V[k]+v)*0.5, 0);	
		}

		//�ċߐڋ�����le�ȉ��̏ꍇ�͂�����C��
		for(int k=0;k<newN;k++)
		{
			int i=beforeN+k;					//�Ή����闱�q�ԍ�
			int G_id=index[k];				//�i�[����i�q�ԍ�
			double mindis=le;
			int J=k;						//�ŋߐڋ����̑��藱�q
			for(int II=G_id-1;II<=G_id+1;II++)
			{       
				for(int JJ=-1*grid_sizeX;JJ<=grid_sizeX;JJ+=grid_sizeX)
				{
					int M_id=II+JJ;
					for(int L=0;L<MESH[M_id].size();L++)
					{
						int j=MESH[M_id][L];
						double x=PART[j].Get_X()-PART[i].Get_X();
						double y=PART[j].Get_Y()-PART[i].Get_Y();
						double dis=sqrt(x*x+y*y);
						if(dis<mindis && i!=j)
						{
							mindis=dis;
							J=j;
						}
					}
				}
			}
			if(J!=i && J<beforeN)//le���ߐڂ��Ă��鑊�肪���E���q�Ȃ�
			{
				double L=le-mindis;//�J���ׂ�����
				double dX=PART[J].Get_X()-PART[i].Get_X();
				double dY=PART[J].Get_Y()-PART[i].Get_Y();
				PART[i].Add(-(dX/mindis*L), -(dY/mindis*L), 0);		
			}
			else if(J!=i && J>=beforeN)//le���ߐڂ��Ă��鑊�肪�������q�Ȃ�
			{
				double L=0.5*(le-mindis);//�J���ׂ�����
				double dX=PART[J].Get_X()-PART[i].Get_X();
				double dY=PART[J].Get_Y()-PART[i].Get_Y();
				PART[i].Add(-(dX/mindis*L), -(dY/mindis*L), 0);	
				PART[J].Add(dX/mindis*L, dY/mindis*L, 0);	
			}
		}//////////*/
	}/////MD�I��

	delete [] index;
//	delete [] MESH;
}

int Model :: Set_cylinder_face(int number,double le,double R,double height,int circle_start_id)
{
//���aR,����height�̉~���̖ʂ��쐬����B���������̊֐��Ăяo�����ɂ����āA���łɉ��ʂ̉~(Z=0)�͍쐬�ς݂Ƃ���
	//top_flag=ON�Ȃ�~����ʂ��쐬����BOFF�Ȃ炵�Ȃ����A���ʂ����͍쐬����B
	int beforeN=number;
	int circle_end_id=number;
	int newN=0;

	int Nv;				//�����̕�����
	double dL_V;		//�����̕�������
	double A=sqrt(3.0)/2;		//�悭�g���W��
	/////////��������///////////////////////////
	double temp_N=height/(le*A);			//���̕������Ble�Ŋ���؂ꂽ���Ԃ������ǁA�����������Ȃ��Ƃ�������
	int Ns=(int) temp_N;				//�^�̕�����
	double difference=temp_N-Ns;		//���Ɛ^�̍�
	if(difference>0.5) Ns++;
	dL_V=height/Ns;			//���q�̋���
	Nv=Ns;
	//////////////////////////////////////////////

	int Nr=(int)((2*PI*R)/le);//�~���̕�����
	double Lr=2*PI*R/Nr;				//�~����������

	double gap=0.4*le;				//�ӂ��肬��ɓ������q��z�u���Ȃ��悤�A���Ԃ�݂���

	///////////////////////////////////����
	for(int i=0;i<Nr;i++)
	{
		for(int j=1;j<Nv;j++)//j=0,j=Nv�͉��ʁA��ʂɊY������̂ł����ł͂ʂ���
		{
			double jj=j*dL_V;
			double ii=i*Lr;
			if(j%2!=0) ii+=0.5*Lr;//j����Ȃ�ii��0.5�i�q�������炷
			if(ii<2*PI*R-gap)
			{
				if(jj<height-gap)	
				{
					double theta=2*PI*(ii/(2*PI*R));
					PART.push_back(PART0);
					PART[beforeN+newN].Set(R*cos(theta), R*sin(theta), jj, le);
					newN++;
				}
			}
		}
	}

	////////////////////////


	//��ʍ쐬�i���ʂ̃R�s�[�j������Nv����Ȃ��ʂ͉��ʂƔ��i�q����Ȃ���΂Ȃ�Ȃ�
	beforeN+=newN;
//	cout<<PART.size()<<" "<<beforeN<<endl;
	newN=0;
	if(Nv%2==0)	//�����Ȃ炻�̂܂ܺ�߰
	{
		
			for(int i=circle_start_id;i<circle_end_id;i++)
			{
				PART.push_back(PART0);
				PART[beforeN+newN].Set(PART[i].Get_X(), PART[i].Get_Y(), height, le);
				newN++;
			}
	}
	else
	{
		double d_theta=0.5*Lr/R;//���̔����p�x������]������B
			for(int i=circle_start_id;i<circle_end_id;i++)
			{
				double X2=PART[i].Get_X()*cos(d_theta)-PART[i].Get_Y()*sin(d_theta);//��]��̍��W
				double Y2=PART[i].Get_X()*sin(d_theta)+PART[i].Get_Y()*cos(d_theta);
				PART.push_back(PART0);
				PART[beforeN+newN].Set(X2, Y2, height, le);
				newN++;
			}		
	}
	return beforeN+newN-number;
}

int Model :: Set_cylinder_in(int number,double le, double R, double height, int start_id)
{
	//���aR,����height�̉~���������쐬����B���̊֐��Ăяo������0<=i<number�̗��q�ŉ~���\�ʂ��`������Ă���Ƃ���B
	int newN=0;					//�̊֐��ŐV�����ǉ����闱�q��
	int beforeN=number;		//���̊֐��Ăяo�����ɂ����闱�q��

	double A=sqrt(3.0)/2;				//�悭�g���W��
	double B=sqrt(2.0/3);						////�悭�g���W��
	int WX=(int)(R/le)+1;		 //�����\���܂ގl�p�`��z�肷��B���̎l�p�`�̕�*0.5
	int WY=(int)(R/(le*A))+1;  //�����\���܂ސ��l�p�`��z�肷��B���̎l�p�`�̕�*0.5
	int WZ=(int)(height/(le*B))+1;  //�����\���܂ސ��l�p�`��z�肷��B���̎l�p�`�̕�*0.5
	double R2=R-0.3*le;				//���������߂̔��a��ݒ�
	double height2=height-0.3*le;				//���������߂̍�����ݒ�

	//���������ʒu
	for(int i=-WX;i<=WX;i++)
	{
		for(int j=-WY;j<=WY;j++)
		{
			for(int k=1;k<=WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//j����Ȃ�ii��0.5�i�q�������炷
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//k����Ȃ�ii��jj�����炷
				if(ii*ii+jj*jj<R2*R2)
				{
					if(kk<height2)
					{
						PART.push_back(PART0);
					//	cout<<PART.size()<<"   "<<beforeN+newN<<endl;
						PART[beforeN+newN].Set(ii, jj, kk, le);
						newN++;
					}
				}
			}
		}
	}///////////////////////*/
	//���q���͊w�ɂ��ʒu���œK��
	double r=1.5;
	double rigion[3][2];	//��͗̈�
	rigion[A_X][0]=-1.2*R; rigion[A_X][1]=1.2*R;
	rigion[A_Y][0]=-1.2*R; rigion[A_Y][1]=1.2*R;
	rigion[A_Z][0]=-5*le;  rigion[A_Z][1]=height*1.5;

	MD_3D(le,start_id,beforeN,newN,r,rigion);  //0��stert_id�ɕύX
	return newN;
}

void Model :: MD_3D(double le,int BstartID,int beforeN,int newN,double r,double region[3][2])
{
//���q���͊w�ɂ��newN�̗��q�̈ʒu���œK���@ID��BstartID����BendID�܂ł̂͋��E���q�Ȃ̂œ������Ȃ�
	double k0=1;
	double dt=0.001;
	int BendID=beforeN;
	
	//�͂�ax^3+bx^2+d�̎����̗p�B����[Bubble Mesh Automated Triangular Meshing of Non-Manifold Geometry by Sphere Packing]���Q��
	double a=(r+1)/(2*r*r-r-1)*k0/(le*le);
	double b=-0.5*k0/le-1.5*a*le;
	double d=-a*le*le*le-b*le*le;
	/////////////
	int lastN=beforeN+newN;

	//cout<<"F="<<a*le*le*le+b*le*le+d<<" "<<a*1.5*le*1.5*le*1.5*le+b*1.5*le*1.5*le+d<<endl;

	vector<double> Fx(newN);	//�e���q�ɓ���X������
	vector<double> Fy(newN);	//�e���q�ɓ���Y������
	vector<double> Fz(newN);	//�e���q�ɓ���Z������
	vector<double> Ax(newN,0);	//X���������x
	vector<double> Ay(newN,0);	//Y���������x
	vector<double> Az(newN,0);	//Z���������x
	vector<double> U(newN,0);	//X�������x
	vector<double> V(newN,0);	//Y�������x
	vector<double> W(newN,0);	//Z�������x
	vector<double> visX(newN);	//X�����S���W��
	vector<double> visY(newN);	//Y�����S���W��
	vector<double> visZ(newN);	//Y�����S���W��
	vector<double> KX(newN);	//X�����o�l�W��
	vector<double> KY(newN);	//Y�����o�l�W��
	vector<double> KZ(newN);	//Y�����o�l�W��

	//�v�Z�̍������̂��߂Ɋi�q���`�� ��͕���r*le�Ŋ���؂��Ƃ͌���Ȃ��̂ŁA�͂ݏo�����Ƃ���͐؂�̂āB�Ȃ̂Ŋe���Ƃ����̕����ɂ͗]�T��������
	double grid_width=le*((int)(r+1));								//�i�q�̕��Br���܂ސ���*le
	int grid_sizeX=(int)((region[A_X][1]-region[A_X][0])/grid_width);	//X�����̊i�q�̌�
	int grid_sizeY=(int)((region[A_Y][1]-region[A_Y][0])/grid_width);
	int grid_sizeZ=(int)((region[A_Z][1]-region[A_Z][0])/grid_width);
	int grid_SIZE=grid_sizeX*grid_sizeY*grid_sizeZ;
	int plane_SIZE=grid_sizeX*grid_sizeY;
	int *index=new int[newN];									//�e�������q���܂ފi�q�ԍ�
	vector<int> *MESH=new vector<int>[grid_SIZE];				//�e���b�V���Ɋi�[����闱�qID�i�[

	for(int i=BstartID;i<BendID;i++)	//�܂��͋��E���q���i�q�Ɋi�[
	{
		int xn=(int)((PART[i].Get_X()-region[A_X][0])/grid_width);//X�����ɉ��ڂ̊i�q�� 
		int yn=(int)((PART[i].Get_Y()-region[A_Y][0])/grid_width);//Y�����ɉ��ڂ̊i�q��
		int zn=(int)((PART[i].Get_Z()-region[A_Z][0])/grid_width);//Z�����ɉ��ڂ̊i�q��
		int number=zn*grid_sizeX*grid_sizeY+yn*grid_sizeX+xn;//���qi���܂ފi�q�̔ԍ�
		MESH[number].push_back(i);
	}
	for(int k=0;k<newN;k++)	//���ɓ������q���i�[
	{
		int i=beforeN+k;
		int xn=(int)((PART[i].Get_X()-region[A_X][0])/grid_width);//X�����ɉ��ڂ̊i�q�� 
		int yn=(int)((PART[i].Get_Y()-region[A_Y][0])/grid_width);//Y�����ɉ��ڂ̊i�q��
		int zn=(int)((PART[i].Get_Z()-region[A_Z][0])/grid_width);//Z�����ɉ��ڂ̊i�q��
		int number=zn*grid_sizeX*grid_sizeY+yn*grid_sizeX+xn;//���qi���܂ފi�q�̔ԍ�
		MESH[number].push_back(i);
		index[k]=number;
	}

	//�v�Z�J�n
	for(int t=0;t<100;t++)
	{
		if(t%10==0 && t>0)
		{
			//MESH����x�j�󂷂�B
			for(int n=0;n<grid_SIZE;n++)
			{
				size_t size=MESH[n].size();
				for(int k=0;k<size;k++) MESH[n].pop_back();
			}
			
			for(int i=BstartID;i<BendID;i++)	//�܂��͋��E���q���i�q�Ɋi�[
			{
				int xn=(int)((PART[i].Get_X()-region[A_X][0])/grid_width);//X�����ɉ��ڂ̊i�q�� 
				int yn=(int)((PART[i].Get_Y()-region[A_Y][0])/grid_width);//Y�����ɉ��ڂ̊i�q��
				int zn=(int)((PART[i].Get_Z()-region[A_Z][0])/grid_width);//Z�����ɉ��ڂ̊i�q��
				int number=zn*grid_sizeX*grid_sizeY+yn*grid_sizeX+xn;//���qi���܂ފi�q�̔ԍ�
				MESH[number].push_back(i);
			}
			for(int k=0;k<newN;k++)	//���ɓ������q���i�[
			{
				int i=beforeN+k;
				int xn=(int)((PART[i].Get_X()-region[A_X][0])/grid_width);//X�����ɉ��ڂ̊i�q�� 
				int yn=(int)((PART[i].Get_Y()-region[A_Y][0])/grid_width);//Y�����ɉ��ڂ̊i�q��
				int zn=(int)((PART[i].Get_Z()-region[A_Z][0])/grid_width);//Z�����ɉ��ڂ̊i�q��
				int number=zn*grid_sizeX*grid_sizeY+yn*grid_sizeX+xn;//���qi���܂ފi�q�̔ԍ�
				MESH[number].push_back(i);
				index[k]=number;
			}
		}

		for(int k=0;k<newN;k++)
		{
			Fx[k]=0; Fy[k]=0, Fz[k]=0;			//������
			KX[k]=0;KY[k]=0; KZ[k]=0;			//�o�l�W��
		}

		for(int k=0;k<newN;k++)
		{
			int i=beforeN+k;					//�Ή����闱�q�ԍ�
			int G_id=index[k];				//�i�[����i�q�ԍ�
			for(int II=G_id-1;II<=G_id+1;II++)
			{       
				for(int JJ=-1*grid_sizeX;JJ<=grid_sizeX;JJ+=grid_sizeX)
				{
					for(int KK=-1*plane_SIZE;KK<=plane_SIZE;KK+=plane_SIZE)
					{
						int M_id=II+JJ+KK;
						for(int L=0;L<MESH[M_id].size();L++)
						{
							int j=MESH[M_id][L];
							if(j>=beforeN && j>i)	//�����̈���ł���i���傫�Ȕԍ��Ȃ�
							{
								int J=j-beforeN;	//newN���ł̔ԍ�
								double x=PART[j].Get_X()-PART[i].Get_X();
								double y=PART[j].Get_Y()-PART[i].Get_Y();
								double z=PART[j].Get_Z()-PART[i].Get_Z();
								double dis=sqrt(x*x+y*y+z*z);
								if(dis<r*le)			//����loop�͎������g���ʉ߂��邩��Adis!=0�͕K�v
								{
									double F=a*dis*dis*dis+b*dis*dis+d;
									Fx[k]-=F*x/dis;					//F�̒l�����̂Ƃ��͐˗͂Ȃ̂ŁA-=�ɂ���
									Fy[k]-=F*y/dis;
									Fz[k]-=F*z/dis;
									Fx[J]+=F*x/dis;					//���藱�q�̗͂������Ōv�Z�B�����͔��]������
									Fy[J]+=F*y/dis;
									Fz[J]+=F*z/dis;
									double K=3*a*dis*dis+2*b*dis;//�o�l�W���@�͂̎��̔����ɑ���
									K=sqrt(K*K);					//���̒l���~�����B�����畉�̂Ƃ��ɔ����Đ��ɕϊ�
									KX[k]+=K*x*x/(dis*dis);			//k���e�����ɕ��z�B�����ŁA��ɐ��̗ʂ����z�����悤��x*x/(dis*dis)�ƂȂ��Ă���
									KY[k]+=K*y*y/(dis*dis);
									KZ[k]+=K*z*z/(dis*dis);
									KX[J]+=K*x*x/(dis*dis);			//k�𑊎藱�q�ɂ����z
									KY[J]+=K*y*y/(dis*dis);
									KZ[J]+=K*z*z/(dis*dis);
								}
							}
							if(j<BendID && j>=BstartID)
							{
								double x=PART[j].Get_X()-PART[i].Get_X();
								double y=PART[j].Get_Y()-PART[i].Get_Y();
								double z=PART[j].Get_Z()-PART[i].Get_Z();
								double dis=sqrt(x*x+y*y+z*z);
								if(dis<r*le && dis>0)			//����loop�͎������g�͒ʉ߂��Ȃ��Adis!=0�͕s�v
								{
									double F=a*dis*dis*dis+b*dis*dis+d;
									Fx[k]-=F*x/dis;					//F�̒l�����̂Ƃ��͐˗͂Ȃ̂ŁA-=�ɂ���
									Fy[k]-=F*y/dis;
									Fz[k]-=F*z/dis;
									double K=3*a*dis*dis+2*b*dis;//�o�l�W���@�͂̎��̔����ɑ���
									K=sqrt(K*K);					//���̒l���~�����B�����畉�̂Ƃ��ɔ����Đ��ɕϊ�
									KX[k]+=K*x*x/(dis*dis);			//k���e�����ɕ��z�B�����ŁA��ɐ��̗ʂ����z�����悤��x*x/(dis*dis)�ƂȂ��Ă���
									KY[k]+=K*y*y/(dis*dis);
									KZ[k]+=K*z*z/(dis*dis);
								}
							}
						}
					}
				}
			}
			//visX[k]=1.414*sqrt(KX[k]);//���̂悤�Ɋe�������̔S���W�������߂�B�����u�������f���ɂ�鎩�����b�V�������vP6�Q�ƁB���������ʂ�1�Ƃ��Ă���B
			//visY[k]=1.414*sqrt(KY[k]);
			//visZ[k]=1.414*sqrt(KZ[k]);
			visX[k]=1.414*sqrt(KX[k]);//���̂悤�Ɋe�������̔S���W�������߂�B�����u�������f���ɂ�鎩�����b�V�������vP6�Q�ƁB���������ʂ�1�Ƃ��Ă���B
			visY[k]=1.414*sqrt(KY[k]);
			visZ[k]=1.414*sqrt(KZ[k]);
			Ax[k]=(Fx[k]-visX[k]*U[k]);
			Ay[k]=(Fy[k]-visY[k]*V[k]);
			Az[k]=(Fz[k]-visZ[k]*W[k]);
		}//�e���q�̉����x�����܂����B
		
		if(t==0)	//�ŏ��̽ï�ߎ���dt������
		{
			double MaxAccel=0;
			for(int k=0;k<newN;k++)
			{
				double accel2=Ax[k]*Ax[k]+Ay[k]*Ay[k]+Az[k]*Az[k];
				if(accel2>MaxAccel) MaxAccel=accel2;
			}
			MaxAccel=sqrt(MaxAccel);//�ő�����x�����܂���
			if(MaxAccel!=0)
			{
				dt=sqrt(0.02*le/MaxAccel);
			}
		}

		for(int k=0;k<newN;k++)//���x�ƈʒu�̍X�V
		{
			int i=beforeN+k;
			double u=U[k];
			double v=V[k];
			double w=W[k];
			U[k]+=dt*Ax[k];
			V[k]+=dt*Ay[k];
			W[k]+=dt*Az[k];
			PART[i].Add(dt*(U[k]+u)*0.5, dt*(V[k]+v)*0.5, dt*(W[k]+w)*0.5);
		}

		//�ċߐڋ�����le�ȉ��̏ꍇ�͂�����C��
		for(int k=0;k<newN;k++)
		{
			int i=beforeN+k;					//�Ή����闱�q�ԍ�
			int G_id=index[k];				//�i�[����i�q�ԍ�
			double mindis=le;
			int J=k;						//�ŋߐڋ����̑��藱�q
			for(int II=G_id-1;II<=G_id+1;II++)
			{       
				for(int JJ=-1*grid_sizeX;JJ<=grid_sizeX;JJ+=grid_sizeX)
				{
					for(int KK=-1*plane_SIZE;KK<=plane_SIZE;KK+=plane_SIZE)
					{
						int M_id=II+JJ+KK;
						for(int L=0;L<MESH[M_id].size();L++)
						{
							int j=MESH[M_id][L];
							double x=PART[j].Get_X()-PART[i].Get_X();
							double y=PART[j].Get_Y()-PART[i].Get_Y();
							double z=PART[j].Get_Z()-PART[i].Get_Z();
							double dis=sqrt(x*x+y*y+z*z);
							if(dis<mindis && i!=j)
							{
								mindis=dis;
								J=j;
							}
						}
					}
				}
			}
			if(J!=i && J<beforeN)//le���ߐڂ��Ă��鑊�肪���E���q�Ȃ�
			{
				double L=le-mindis;//�J���ׂ�����
				double dX=PART[J].Get_X()-PART[i].Get_X();
				double dY=PART[J].Get_Y()-PART[i].Get_Y();
				double dZ=PART[J].Get_Z()-PART[i].Get_Z();
				PART[i].Add(-dX/mindis*L, -dY/mindis*L, -dZ/mindis*L);
			}
			else if(J!=i && J>=beforeN)//le���ߐڂ��Ă��鑊�肪�������q�Ȃ�
			{
				double L=0.5*(le-mindis);//�J���ׂ�����
				double dX=PART[J].Get_X()-PART[i].Get_X();
				double dY=PART[J].Get_Y()-PART[i].Get_Y();
				double dZ=PART[J].Get_Z()-PART[i].Get_Z();
				PART[i].Add(-dX/mindis*L, -dY/mindis*L, -dZ/mindis*L);
				PART[J].Add(dX/mindis*L, dY/mindis*L, dZ/mindis*L);
			}
		}//////////*/
	}/////MD�I��

	delete [] index;
	delete [] MESH;
}

//�����𕪊����邳���̍œK�ȕ������ƕ��������̎Z�o�֐�
void Model :: Set_calc_N_and_L(double dis,double le,int *N,double *L)
{
	double temp_N=dis/le;			//���̕������Ble�Ŋ���؂ꂽ���Ԃ������ǁA�����������Ȃ��Ƃ�������
	int Ns=(int) temp_N;				//�^�̕�����
	double difference=temp_N-Ns;		//���Ɛ^�̍�
	if(difference>0.5) Ns++;
	*L=dis/Ns;			//���q�̋���
	*N=Ns;
}

//�~���������v�Z�֐�
int Model :: Set_calc_division_N_circle(double dis,double le)
{
	//�Ώ̐����l�������A�~�����L�q���闱�q�͋����łȂ���΂Ȃ�Ȃ��B�����瑼�̕ӕ������Ƃ͈�������������
	//dis:�������鋗��(�~��)
	double temp_num=dis/le;		//�~�O���ɐݒu����w���́x���q���B�������O�������܂�le�Ŋ���؂��Ƃ͌���Ȃ�

	int N1=(int)(temp_num/2);
	N1*=2;							
	int N2=N1+2;					//temp_num��N1��N2�̊Ԃɂ���B������N1,N2�͋���

	double dif1=temp_num-N1;		//�eN�Ƃ̍�
	double dif2=N2-temp_num;
	int N=N1;						//������������
	if(dif2<dif1) N=N2;				//���̏���������N�Ƃ��č̗p����B

	return N;
}