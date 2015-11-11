#include "stdafx.h"	
#include "PostProcess.h"

#define FULL 3
#define REMESH 4
#define FULL_IMPORT 5

using namespace std;

/*********�p����
�ÓI�v�f�F�����b�V�����ɍĔz�u���Ȃ��v�f�i������Ԃ̂��̂��g��������j
���I�v�f�F�����b�V�����ɍĔz�u����v�f
*****************/
//�v�f�̐[������֐�
void set_depth(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,int *depth,int KTE);
void non_linear_Avector3D(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,double *V,int *jnb,int *branch_num,double **current,double *RP,double II,int *depth,double **Be,int t);
//�d���x�N�g��
void carrent_vector(mpsconfig &CON, vector<point3D> &NODE, vector<element3D> &ELEM, int nelm, double **B,int t);
void FEM3D_calculation(mpsconfig &CON, int *static_node, int *static_nelm, vector<point3D> &static_NODE, vector<element3D> &static_ELEM, int t,double TIME, vector<mpselastic> &PART, int fluid_number, int particle_number, double dt, double **F)
{
	int delaun_flag;//�f���[�j�������s�����A�s��Ȃ���

	if(CON.get_mesh_input()==0)//MPSTOFEM�ɂ��ߓ_�A�v�f����
	{
		if(CON.get_remesh_sw()==OFF) delaun_flag=FULL;	//remesh�̈��z�肹���A��ɂ��ׂĂ��f���[�j����
		else if(CON.get_remesh_sw()==ON)
		{
			if(t==1) delaun_flag=FULL;	//�S���f�����f���[�j����
			else delaun_flag=REMESH;	//remesh�̈�̂݃f���[�j����
		}
	}

	if(delaun_flag==FULL) cout<<"FULL �f���[�j�������s"<<endl;
	else if(delaun_flag==REMESH) cout<<"remesh�̈�̂݃f���[�j�������s"<<endl;
	else if(delaun_flag==FULL_IMPORT) cout<<"Magnet�����t�@�C�����v�f��񓙓ǂݍ���"<<endl;

	vector<point3D> NODE;
	vector<element3D> ELEM;
	int node=0;					//�S�ߓ_��
	int nelm=0;					//�S�v�f��
	int KTJ;					//�ő�ߓ_��
	int KTE;					//�ő�v�f���@3���������Ƃ߂�
	const double err=1.0e-14;			//�덷����̂������l�E�E�E1e-14�}�V���C�v�V������茵�����̂ł́H
	
	if(delaun_flag==FULL)
	{
		//���q�z�u���ߓ_�z�u�����
		MPS_TO_FEM3D(CON, &node, NODE, PART, fluid_number, particle_number);

		KTJ=node;
		if(CON.get_fine()!=OFF) KTJ+=CON.get_add_points();
		KTE=12*KTJ;
	}else if(delaun_flag==REMESH){
		node=(int)static_NODE.size()-1;			//�ÓI�ߓ_��
		nelm=(int)static_ELEM.size()-1;
		KTJ=node+particle_number;				//���̂��Ɠ��I�ߓ_(����)���i�[���Ȃ��Ƃ����Ȃ�����AKTJ�𑝉�
		if(CON.get_fine()!=OFF) KTJ+=CON.get_add_points();
		KTE=12*KTJ;
	}

	if(delaun_flag==FULL)//���ׂĂ��f���[�j����
	{
		point3D NODE0;
		element3D ELEM0;
		for(int i=KTJ+8;i>node;i--) NODE.push_back(NODE0);//10^5�Ƃ�int�̍ő�l������̂ł́H�H�H	//����͋�C�ߓ_���̒ǉ��H15/2/4
		for(int i=0;i<KTE;i++) ELEM.push_back(ELEM0);//�z����m��
		int FINE_sw=CON.get_fine();//�����ĕ����X�C�b�`

		//�f���[�j����
		delaun3D_main(CON,NODE,ELEM, KTJ, KTE,&node,&nelm, FINE_sw);
		cout<<"�v�f��="<<nelm<<" �ߓ_����"<<node<<endl;

		//���b�V�������̃|�X�g����
		double *val=new double[KTJ+1];
		for(int i=1;i<=node;i++) val[i]=NODE[i].material;
//		if(t==1 || t%(CON.get_EM_interval()*CON.get_mesh_output_interval())==0) 
		{
			data_avs(node,nelm,NODE,ELEM,KTJ,val,CON);
			data_avs2(CON,node,nelm,NODE,ELEM,KTJ,t);//�f�ʐ}
			data_avs3(node,nelm,NODE,ELEM,CON);//�ގ�	
		}

		if(CON.get_remesh_sw()==ON)//remesh�̈��z�肷��Ȃ�A�ÓI�ߓ_�E�v�f�����L�����Ă���
		{
			//NODE,ELEM�̂����A�����Ȃ��v�f�A�ߓ_������static_NODE,staticELEM�Ɋi�[����
			static_ELEM.clear();
			static_NODE.clear();
			memorize_static_NODE_and_ELEM(CON,NODE,ELEM,static_NODE,static_ELEM, node, nelm);

			/*/�`�F�b�N
			int snode=(int) static_NODE.size()-1;
			int snelm=(int) static_ELEM.size()-1;
			for(int i=1;i<=snode;i++) val[i]=static_NODE[i].remesh;
			data_avs2(CON,snode,snelm,static_NODE,static_ELEM,KTJ,t);//�f�ʐ}*/
		}
		delete [] val;
	}

//���b�V�����؂�Ă���͂̌v�Z�܂�

	//�ߓ_-�v�f�֌W
    int *jnb=new int[node+1];//�e�ߓ_�ɗאڂ���v�f���i�[
    set_jnb3D(NODE,ELEM,node,nelm,jnb);

	int **nei=new int* [node+1];//�e�ߓ_�̎��ӗv�f�ԍ��i�[
    for(int i=1;i<=node;i++) nei[i]=new int [jnb[i]+1];
    set_nei3D(NODE,ELEM,node,nelm,jnb,nei);

	for(int i=1;i<=node;i++)
	{
		if(jnb[i]==0)
		{
			//cout<<"jnb=0 i="<<i<<" material="<<NODE[i].material<<" particle="<<NODE[i].particleID<<endl;
			NODE[i].boundary_condition=1;//���E�������f�B���N���^�ɂ��邱�ƂŁAICCG�ɎQ�������Ȃ�
			//if(NODE[i].material==FLUID) NODE[i].particleID=-1;//���̐ߓ_����������ꍇ�A�K�v�Ȓl���v�Z����Ȃ����߁A���q�Ƀt�B�[�h�o�b�N�ł��Ȃ�(���Ă͂����Ȃ�)
			if(NODE[i].material==FLUID || NODE[i].material==ELASTIC) if(NODE[i].particleID>=0) cout<<"suf="<<PART[NODE[i].particleID].surface<<endl; 
		}
	}

	//FEM
	if(CON.get_EM_calc_type()==1 || CON.get_EM_calc_type()==4) potential_calculation(CON,NODE,ELEM, node, nelm,jnb, TIME,PART, fluid_number,nei,F);
	if(CON.get_EM_calc_type()==2) calc_static_magnetic_field(CON, node, nelm,NODE,ELEM,jnb, dt,TIME, t,nei, KTE,PART,fluid_number,F,KTJ);
	if(CON.get_EM_calc_type()==3) calc_transitional_EM_field(CON, node, nelm,NODE,ELEM,jnb, dt, TIME,t,nei, KTE,PART,fluid_number,F,KTJ);
	if(CON.get_EM_calc_type()==5) calc_variable_magnetic_field(CON, node, nelm,NODE,ELEM,jnb, dt,TIME, t,nei, KTE,PART,fluid_number,F);

	//FEM_interval�ɂ��v�Z���ԒZ�k�̏ꍇ
    if(CON.get_EM_interval()>1)
    {
		//ofstream bb("FEM_interval.dat");///�d���͏o�́@���X�e�b�v�͂����ǂݎ��
        FILE *b=fopen("FEM_interval.dat","w");///�d���͏o�́@���X�e�b�v�͂����ǂݎ��
		for(int i=0;i<fluid_number;i++)
		{
			//bb<<F[A_X][i]<<" "<<F[A_Y][i]<<" "<<F[A_Z][i]<<endl;
			fprintf(b,"%1.15lf %1.15lf %1.15lf\n",F[A_X][i],F[A_Y][i],F[A_Z][i]);
		}
		fclose(b);
		//bb.close();
    }

	delete [] jnb;
    for(int i=1;i<=node;i++) delete [] nei[i];
    delete [] nei;
}

//�Q�d�����l���ɂ��ꂽ�ߓ_�v�f�̃x�N�g���|�e���V�����v�Z�֐�ver.2 �N�[�����Q�[�W��K�p���č��ӂ����v���V�A���Ƃ���
void Avector3D_node_eddy2(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,double **A,int *jnb,int *depth,int **nei2,int *branch_num,double **old_A,double II,double dt,int mps_num,double *V,int t)
{	cout<<"�Q�d�����l���ɂ��ꂽ�ߓ_�v�f�ɂ���޸�����ݼ�ٌv�Z�J�n"<<endl;

	double u0=4*PI*0.0000001;			//��C�̓�����
    double v0=1/u0;						//���C��R��
	double j0x,j0y,j0z;					//�d�����x
	double sigma=CON.get_ele_conduc();	//�d�C�`����
	
	double *current[3];					//�e�v�f�̓d�����x[A/m3]�i�[
	for(int D=0;D<3;D++) current[D]=new double [nelm+1];

	//////////////////////��ق��Ȃ���Γd���v�Z���s���K�v���Ȃ��B�����ź�ِߓ_������������
	int coil_node_num=0;	//��ِߓ_��
	for(int i=1;i<=node;i++) if(NODE[i].material==COIL) coil_node_num++;
	////////////////////////*/

	int *save_bound=new int [node+1];//�{�֐��͉��x���Ăяo�����̂ŁA�d���Ɋւ��鋫�E���������S�ɏ����ł��Ȃ��B�����ŕۑ�����
	
	if(coil_node_num>0)///��ِߓ_������Ȃ�d�����x�v�Z
	{
		for(int i=1;i<=node;i++) save_bound[i]=NODE[i].boundary_condition;//���E������ۑ�

		if(II!=0)//�d���l����[���Ȃ�d���v�Z
		{
			calc_current(CON,NODE,ELEM,SIDE,node,nelm,side_num,jnb,branch_num,current,depth,II);
			cout<<"�d���v�Z����"<<endl;
		}
		else	//�[���Ȃ珉����
		{
			for(int i=1;i<=nelm;i++) if(ELEM[i].material==COIL) for(int D=0;D<3;D++) current[D][i]=0;
		}
		///�d�����E�����̏������i�����d���͂��Ƃ܂��Ă��邩�炢��Ȃ�)
		for(int i=1;i<=side_num;i++) if(SIDE[i].boundary_condition>=10) SIDE[i].boundary_condition=0;
		for(int i=1;i<=node;i++) if(NODE[i].boundary_condition>=10) NODE[i].boundary_condition=0;
	}
	carrent_vector(CON, NODE, ELEM, nelm, current,t);

    int NN=0;//�f�B���N���^���E�ߓ_��
    int *dn=new int [node+1]; //�e�ߓ_���f�B���N���^���E�̉��Ԗڂ��B�f�B���N���^���E��łȂ��Ȃ�node+1���i�[
	double *PHAT[3];
	for(int D=0;D<3;D++) PHAT[D]=new double [CON.get_max_DN()];//�f�B���N���^���l
    
    ///�f�B���N���^���E��������
    for(int i=1;i<=node;i++)
    {
        if(jnb[i]==0)
		//if(NODE[i].boundary_condition==3)
		{    
			//cout<<"Jnb["<<i<<"]=0"<<endl;
	        dn[i]=NN;
	        PHAT[A_X][NN]=0;
			PHAT[A_Y][NN]=0;
			PHAT[A_Z][NN]=0;
			A[A_X][i]=0;
			A[A_Y][i]=0;
			A[A_Z][i]=0;
	        NN++;
		}
		else if(NODE[i].boundary_condition==2)
		{    
	        dn[i]=NN;
	        PHAT[A_X][NN]=0;
			PHAT[A_Y][NN]=0;
			PHAT[A_Z][NN]=0;
			A[A_X][i]=0;
			A[A_Y][i]=0;
			A[A_Z][i]=0;
	        NN++;
		}
		else if(NODE[i].boundary_condition==1)
		{    
	        dn[i]=NN;
	        PHAT[A_X][NN]=0;
			PHAT[A_Y][NN]=0;
			PHAT[A_Z][NN]=0;
			A[A_X][i]=0;
			A[A_Y][i]=0;
			A[A_Z][i]=0;
	        NN++;
		}
		else
		{
			dn[i]=node+1;
		}
    }/////////////*/
    cout<<"�ިظڐ���"<<3*NN<<endl;//���ۂɌv�Z���Ȃ��Ă����ިظڌ^�ߓ_����3*NN�ł���
	    
	///�`�Ӗ@�ɂ͓���(����)�ߓ_���̐������d�ʂɊւ��関�m�������݂���B����𐔂���

	int conducter_num=0;//����(����)�ߓ_��
	for(int i=1;i<=node;i++) if( NODE[i].material==WATER && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
	//////////////
    int pn=3*(node-NN)+conducter_num;///���m��
    int *ppn=new int [pn];		///�s���n�Ԗڂ̐ߓ_�͐ߓ_�ԍ�ppn[n]
    int *npp=new int [node+1];	///�e�ߓ_���s��̉��Ԗڂɂ����邩�B�ިظڌ^�̏ꍇ��pn+1���i�[
    int num=0; 
    for(int i=1;i<=node;i++)//�s���A1x,A1y,A1z,A1��,A2x,A2y�E�E�E�̏��Ɋi�[
    {
        if(NODE[i].boundary_condition==0 && jnb[i]!=0)
		//if(NODE[i].boundary_condition==0)
		{
			ppn[num]=i;
			ppn[num+1]=-1;
			ppn[num+2]=-2;
			npp[i]=num;
			num+=3;
			if( NODE[i].material==WATER)
			{
				//cout<<"�ߓ_�ԍ�"<<i<<endl;
				ppn[num]=-3;
				num++;
			}
		}
		else npp[i]=pn+1;
    }
    
    ////�s��̍ő啝�v�Z  �����Ƃ����̑傫���̂̓ӂ̗v�f���낤�Ɖ��肵�Ă���B�m���ł͂Ȃ����Ƃɒ���
    int mat_w=0;
	for(int i=1;i<=node;i++) if(branch_num[i]+1>mat_w && NODE[i].material==WATER) mat_w=branch_num[i]+1;
	mat_w=4*mat_w;//X,Y,Z,�Ӑ����Ƃ���̂Ł~4 //���ӂ����׼�݂ɂ��Ă��A�ӂɊւ��Ă͂S�����̂܂܂ɒ���
    //////
//	mat_w=200;//��c�̏ꍇ�Q�d��������Ȃ�
    
    ////�z��m��
    double **G=new double *[pn+1];///�S�̍s��
    for(int i=1;i<=pn;i++) G[i]=new double [mat_w+1];
    int **ROW=new int *[pn+1]; ///�e�s�́A��[���v�f�̗�ԍ��L��
    for(int i=1;i<=pn;i++) ROW[i]=new int [mat_w+1];
    int *NUM=new int [pn+1]; ///�e�s�́A��[���v�f��
    
    for(int i=1;i<=pn;i++)//������
    {
        NUM[i]=0;
        for(int j=1;j<=mat_w;j++)
		{
			G[i][j]=0;
			ROW[i][j]=0;
		}
    }
    double *B=new double [pn];//���s��
    for(int i=0;i<pn;i++) B[i]=0;//������
    ////
    
    /////////�S�̍s����쐬����
    unsigned int time=GetTickCount();
    int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int J,J2,J3,flag,flag2,flag3;
    for(int je=1;je<=nelm;je++)
    {   
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
		}
		
		///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
	
		double delta6=ELEM[je].volume;//�̐ς�6�{
	
		delta6=1/delta6/6;//�v�Z�ɕK�v�Ȃ̂�1/(36V)�Ȃ̂ł����Ŕ��]���Ă���(1/36V^2)*V=1/36V
	
		double delta=ELEM[je].volume/6;//�{���̑̐�
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
			int m=j%4+1;
			int n=m%4+1;
	    
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
			if(i%2!=0)//i����Ȃ�
			{
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
		
		if(ELEM[je].material==COIL)
		{
			j0x=current[A_X][je];
			j0y=current[A_Y][je];
			j0z=current[A_Z][je];
		}
		else
		{
			j0x=0;
			j0y=0;
			j0z=0;
		}

		///�䓧����
		double v=v0;
		if(ELEM[je].material==WATER) v/=CON.get_RP();
//		if(ELEM[je].material==IRON) v/=5000;
		double rp=u0;
		if(ELEM[je].material==WATER) rp*=CON.get_RP();
//		if(ELEM[je].material==IRON) rp*=5000;
		
		for(int n=1;n<=4;n++)
		{
			if(NODE[N[n]].boundary_condition==0)
			{
				/////X����

				int I=npp[N[n]]+1;
				//�e�a�͂����ł܂Ƃ߂Čv�Z����
				B[I-1]+=delta/4*j0x;//x,y,z�͂P�������
				B[I]+=delta/4*j0y;//
				B[I+1]+=delta/4*j0z;
				
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						J=npp[N[m]]+1;	//Aix�̍�
						flag=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aix�̍�
							    flag=1;
							}
						}
						if(flag==0)//��������i�����݂��Ȃ���������
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aix�̍�
						    ROW[I][H]=J;
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_X][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}

				/////Y����
				I=npp[N[n]]+1+1;/////Y����
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						J=npp[N[m]]+1;	//Aix�̍�
						J2=J+1;			//Aiy�̍�
						flag=0;
						flag2=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J2==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aiy�̍�
							    flag2=1;
							}
						}
						if(flag2==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=J2;	
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Y][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}
				/////*/

				/////Z����
				I=npp[N[n]]+2+1;
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						J=npp[N[m]]+1;	//Aix�̍�
						J2=J+1;			//Aiy�̍�
						J3=J2+1;		//Aiz�̍�
						flag3=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J3==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aiz�̍�
							    flag3=1;
							}
						}
						if(flag3==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=J3;
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Z][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}
			}
		}
	}
	
	//////�Q�d�����v�Z
	int J4,flag4;
	
	for(int je=1;je<=nelm;je++)
    {   
		if(ELEM[je].material==WATER)//MAGELAST�ɂ͉Q�d��������Ȃ�
		{
			
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
			double Xs=0;//�d�S���W
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
		
			///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
			//ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
	
			double delta6=ELEM[je].volume;//�̐ς�6�{
	
			delta6=1/delta6/6;//�v�Z�ɕK�v�Ȃ̂�1/(36V)�Ȃ̂ł����Ŕ��]���Ă���(1/36V^2)*V=1/36V
	
			double delta=ELEM[je].volume/6;//�{���̑̐�
	
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
				int m=j%4+1;
				int n=m%4+1;
	    
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
				if(i%2!=0)//i����Ȃ�
				{
					b[i]*=-1.0;
					c[i]*=-1.0;
					d[i]*=-1.0;
					e[i]*=-1.0;
				}
			}
			/////////
	
			double Sx=0;//�Q�d�����̂����A�P�ï�ߑO���޸�����ݼ�ق��l�����邱�Ƃœ�����A���޸��
			double Sy=0;
			double Sz=0;

			double co=delta/20.0/dt*sigma;//��V/(20dt)�@���x���v�Z���邱�ƂɂȂ�̂ŌW��������
			
			double co2=delta6*sigma;//   ��/(36V)
			double co3=sigma*dt*delta6;

			for(int n=1;n<=4;n++)
			{
				if(NODE[N[n]].boundary_condition==0)
				{
					/////X����

					int I=npp[N[n]]+1;
					Sx=0;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aix�̍�
							J4=J+3;			//�ӂ̍�
							flag=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)//Aix�̍�
							{
								if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									if(I==J) G[I][h]+=co*2;
									else G[I][h]+=co;
									flag=1;
								}
							
								//�ӂ̍�
								if(J4==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
									flag4=1;
								}///
							}
							if(flag==0)//��������i�����݂��Ȃ���������(�����ł͂���Ȃ͂��Ȃ�����)
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==J) G[I][H]+=co*2;
								else G[I][H]+=co;
								ROW[I][H]=J;
							}
							if(flag4==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
								ROW[I][H]=J4;
							}

							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							if(I==J) Sx+=2*Ax;
							else Sx+=Ax;
							
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					B[I-1]+=co*Sx;//�x�z��������X�������瓾����B�̒l

					/////Y����

					I=npp[N[n]]+1+1;/////Y����
					Sy=0;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aix�̍�
							J2=J+1;			//Aiy�̍�
							J4=J+3;			//�ӂ̍�
							flag2=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)
							{
								if(J2==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									if(I==J2) G[I][h]+=co*2;//Aiy�̍�
									else G[I][h]+=co;
									flag2=1;
								}

								//�ӂ̍�
								if(J4==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
									flag4=1;
								}
							}
							if(flag2==0)
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								if(I==J2) G[I][H]+=co*2;
								else G[I][H]+=co;
								ROW[I][H]=J2;
							}
							if(flag4==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
								ROW[I][H]=J4;
							}//*/

							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							double Ay=old_A[A_Y][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							if(I==J2) Sy+=2*Ay;
							else Sy+=Ay;
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							int NN=dn[N[m]];
							//B[I-1]-=-c[n]*d[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+e[n]*e[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}////////*/
					B[I-1]+=co*Sy;//�x�z��������X�������瓾����B�̒l

					/////Z����
					Sz=0;
					I=npp[N[n]]+2+1;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aix�̍�
							J2=J+1;			//Aiy�̍�
							J3=J2+1;		//Aiz�̍�
							J4=J+3;			//�ӂ̍�
							flag=0;
							flag2=0;
							flag3=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Aiz�̍�
								if(J3==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									if(I==J3) G[I][h]+=co*2;
									else G[I][h]+=co;
									flag3=1;
								}
								//�ӂ̍�
								if(J4==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
									flag4=1;
								}
							}
							
							if(flag3==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    if(I==J3) G[I][H]+=co*2;
								else G[I][H]+=co;
							    ROW[I][H]=J3;
							}
							
							if(flag4==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
								ROW[I][H]=J4;
							}//*/
							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							double Az=old_A[A_Z][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							if(I==J3) Sz+=2*Az;
							else Sz+=Az;
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							int NN=dn[N[m]];
							//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}////////*/
					B[I-1]+=co*Sz;//�x�z��������X�������瓾����B�̒l

					/////��
					
					I=npp[N[n]]+3+1;
					Sx=0;
					Sy=0;
					Sz=0;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aix�̍�
							J2=J+1;			//Aiy�̍�
							J3=J2+1;		//Aiz�̍�
							J4=J+3;			//�ӂ̍�
							flag=0;
							flag2=0;
							flag3=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Aix�̍�
								if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									///co2*c[n]*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])�Ƃ����ӂ���c[n]��O�Ɏ����Ă���ƁA�ł��؂�덷�̊֌W�Ŕ�Ώ̂ƔF�������̂ŁA��̎��Ɠ��l�ɍŌ�ɂ����Ă���
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
									flag=1;
								}
								//Aiy�̍�
								if(J2==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
									flag2=1;
								}
								//Aiz�̍�
								if(J3==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
									flag3=1;
								}
								//�ӂ̍�
								if(J4==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
									flag4=1;
								}
							}
							if(flag==0)
							{   
								//cout<<I<<" "<<J<<" co2="<<co2<<"  c[m]="<<c[m]<<" "<<co2*c[m]<<" B  n,m="<<n<<" "<<m<<endl;
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    //G[I][H]+=co2*c[m];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
							    ROW[I][H]=J;
							}
							if(flag2==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
							    ROW[I][H]=J2;
							}
							if(flag3==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
							    ROW[I][H]=J3;
							}
							
							if(flag4==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
								ROW[I][H]=J4;
							}
							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							double Ay=old_A[A_Y][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							double Az=old_A[A_Z][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							Sx+=Ax*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*c[n];
							Sy+=Ay*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*d[n];
							Sz+=Az*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*e[n];
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							int NN=dn[N[m]];
							//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}//////*/
					B[I-1]+=co2*(Sx+Sy+Sz);
				}
			}
		}
	}//////*/
	
    for(int D=0;D<3;D++) delete [] current[D];
    int number=0;//�s��̔�[���v�f��
    for(int i=1;i<=pn;i++) number+=NUM[i];
    
    ///���̂܂܂ł�G[i][j]��[j]�̏��������ɕ���ł��Ȃ��̂ŁAG��ROW����ёւ���
    for(int i=1;i<=pn;i++)
    {
        double tempG;
		int tempR;
		for(int j=2;j<=NUM[i];j++)
		{
		    for(int m=1;m<j;m++)
		    {
		        if(ROW[i][j]<ROW[i][m])
				{
				    tempG=G[i][m];
				    tempR=ROW[i][m];
					G[i][m]=G[i][j];
					ROW[i][m]=ROW[i][j];
					G[i][j]=tempG;
					ROW[i][j]=tempR;
				}
			}
		}
    }//////////*/

	/*for(int i=1;i<=node;i++)
	{
		if(NODE[i].boundary_condition==0)
		{
			//if(branch_num[i]==0) cout<<i<<endl;
			int I=npp[i]+1;
			//if((branch_num[i]+1)*4!=NUM[I]) cout<<branch_num[i]<<" "<<NUM[I]<<" i="<<i<<" I="<<I<<endl;
		}
	}
	//for(int k=1;k<=NUM[4];k++) if(ppn[ROW[4][k]-1]>0) cout<<ppn[ROW[4][k]-1]<<" ";
	//cout<<endl;
	*/

	///�s��̎��ۂ̍ő啝�����߂�
	int maxN=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>maxN) maxN=NUM[i];
	cout<<"�ő啝�F"<<maxN<<"/"<<mat_w;
	
	ofstream fout("matrix.dat");
	for(int i=1;i<=pn;i++)
	{
		for(int j=1;j<=NUM[i];j++)
		{
			fout<<G[i][j]<<" ";
			//cout<<ROW[i][j]<<" ";
		}
		fout<<endl;
		//cout<<endl;
	}
	fout.close();/////////*/

	///�Ώ̐��`�F�b�N
	for(int i=1;i<=pn;i++)
	{
		for(int j=1;j<=NUM[i];j++)
		{
			int J=ROW[i][j];
			for(int k=1;k<=NUM[J];k++)
			{
				if(ROW[J][k]==i)
				{
					if(G[i][j]!=G[J][k]) 
					{
						cout<<"��Ώ�"<<G[i][j]<<" "<<G[J][k]<<endl;
						G[i][j]=G[J][k];
					}
				}
			}
		}
	}///////////*/
	
    
    double *val = new double [number];
    int *ind = new int [number];//��[���v�f�̗�ԍ��i�[�z��
    int *ptr = new int [pn+1];//�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
    
    /////////////////////val,ind ,ptr�ɒl���i�[
    int index=0;
    for(int n=0;n<pn;n++)
    {
        ptr[n]=index;
		for(int m=1;m<=NUM[n+1];m++)
		{
			val[index]=G[n+1][m];
			ind[index]=ROW[n+1][m]-1;
			index++;
		}
    }	    
    ptr[pn]=number;
    ////////////////////*/
    
    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
    
    cout<<"  �s��쐬�I��  time="<<GetTickCount()-time<<endl;
    
	//CG�@���s
	double *XX=new double [pn];//�s��̓����i�[
//    if(CON.get_FEMCG()==0) CG3D_Avector(val,ind,ptr,pn,ppn,B,A);//CG�@���s
//	else 
//	{
		if(CON.get_FEMCG()==1) ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
		else if(CON.get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
		else if(CON.get_FEMCG()==3) MRTR(CON,B,pn,XX,val,ind,ptr);
		else if(CON.get_FEMCG()==4) ICMRTR(CON,B,pn,XX,val,ind,ptr);
		for(int n=0;n<pn;n++)
		{
			int i=ppn[n];
			if(i>0) A[A_X][i]=XX[n];
			else if(i==-1) A[A_Y][ppn[n-1]]=XX[n];
			else if(i==-2) A[A_Z][ppn[n-2]]=XX[n];
			else if(i==-3) V[ppn[n-3]]=XX[n];//�d�ʃ�
		}	
//	}
	delete [] XX;



	/*//old_A�̍X�V //�Q�d���v�Z�̂Ƃ����dA/dt���v�Z����������A�����ł�old_A�͂��̂܂܂ɂ��Ă���
	for(int i=1;i<=mps_num;i++)
	{
		for(int D=0;D<3;D++) old_A[D][i]=A[D][i];
	}///*/

	///old_A��̧�ُo��(���̗��q�̂ݏo�́B�ŏ���mps_num�Ԗڂ܂ł����̗��q)
	ofstream g("old_A.dat");
	for(int i=1;i<=mps_num;i++) g<<old_A[A_X][i]<<" "<<old_A[A_Y][i]<<" "<<old_A[A_Z][i]<<endl;
	g.close();
	
	if(coil_node_num>0)//���E���������Ƃɖ߂�(���̽ï�߂ōĂѓd�����v�Z���邩������Ȃ�����)
	{	
		for(int i=1;i<=node;i++) NODE[i].boundary_condition=save_bound[i];
	}

    delete [] dn;
    for(int D=0;D<3;D++) delete [] PHAT[D];
    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;

	delete [] save_bound;
	
}

//�d�ʁA���ʂȂǂ̃|�e���V�����v�Z�֐�
void potential_calculation(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,int *jnb,double TIME,vector<mpselastic> &PART,int fluid_number,int **nei,double **F)
{
	double *V=new double [node+1];	//potential
    
    double *Ee[3];					//�v�f���d�E or �v�f�����E
    for(int D=0;D<3;D++) Ee[D]=new double [nelm+1];
	double *RP=new double [nelm+1];	//�e�v�f�̓������܂��͗U�d���i�[�i���݂ł͗U�d���͎g���Ă��Ȃ��j

	double fluid_rp=CON.get_r_perm();		//�U�d��
	if(CON.get_EM_calc_type()==4) fluid_rp=CON.get_RP();	//������
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==FLUID || ELEM[i].material==MAGELAST || ELEM[i].material==MAGELAST2) RP[i]=fluid_rp;
		else if(ELEM[i].material==ELASTIC) RP[i]=1;
		else if(ELEM[i].material==AIR) RP[i]=1;
		else if(ELEM[i].material==COIL)
		{
			if(CON.get_EM_calc_type()==1) RP[i]=1000;//�{���͓��̂����灇
			else if(CON.get_EM_calc_type()==4) RP[i]=1;//������
		}
		else cout<<"error:�ގ��ɑ΂��A�U�d�����邢�͓��������s��"<<endl;
	}

	//�d�ʉ���
	VOLT3D(CON,NODE,ELEM, node, nelm,V,jnb, TIME,PART, fluid_number,nei,RP);

	//�d�E�v�Z
	ELECTRO3D(CON,NODE,ELEM, node, nelm,V,Ee);

	//�Ód�͌v�Z�֐� �\�ʂɓ������͂�P=0.5epE^2 ����͗U�d�̂𓱑̋ߎ���������
	//conductive_approximation(CON,NODE,ELEM, node, nelm,Ee, jnb, nei, RP,PART,F,fluid_number);

	

	delete [] V;
    for(int D=0;D<3;D++) delete [] Ee[D];
	delete [] RP;
}

///�d�ʌv�Z�֐�
void VOLT3D(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double *V,int *jnb,double TIME,vector<mpselastic> &PART,int fluid_number,int **nei,double *RP)
{
    double V1=0;
    double V2=CON.get_V();			//�d��
	double u0=0.000001256;			//�^��̓�����
	double ep0=8.854e-12;			//�^��̗U�d��
	double le=CON.get_distancebp();//���q�ԋ���
	unsigned timeA=GetTickCount();	//�v�Z�J�n����

	////////�d�ʌ���
	if(CON.get_EM_calc_type()==1)//�d��v�Z�Ȃ�
	{
		double tau=CON.get_dt()*CON.get_V_step();//���萔
		if(CON.get_V_con()==2)
		{
			V2=CON.get_V()*(1-exp(-TIME/tau));//�d�ʂ��w���֐��I�ɑ���������
		}
		else if(CON.get_V_con()==1)
		{
			V2=(CON.get_V()-CON.get_initial_V())*TIME/(CON.get_dt()*CON.get_V_step())+CON.get_initial_V();//�d�ʂ𒼐��I�ɑ���������
		}
		if(V2==0) V2=1;//0���ƃG���[�ɂȂ邩��1�ɂ���
		else if(V2>CON.get_V()) V2=CON.get_V();
		cout<<"�d�ʌv�Z�J�n V="<<V2<<" ";
	}
	else if(CON.get_EM_calc_type()==4)//���ʌv�Z
	{
		double R=1;//�䗦
		V1=0;
		if(CON.get_uniform_B_sw()==OFF) V2=CON.get_magnet_B()/u0*CON.get_magnet_H();
		if(CON.get_uniform_B_sw()==ON)  V2=CON.get_uniform_B()/u0*(CON.get_ZU()-CON.get_ZD());//��l����
		if(CON.get_V_con()==1) 
		{
			V2=(V2)*TIME/(CON.get_dt()*CON.get_V_step());
			R=TIME/(CON.get_dt()*CON.get_V_step());
			if(V2>CON.get_magnet_B()/u0*CON.get_magnet_H()) {V2=CON.get_magnet_B()/u0*CON.get_magnet_H();R=1;}
		}
		if(V2==0) V2=1;//0���ƃG���[�ɂȂ邩��1�ɂ���
		
		cout<<"���ʌv�Z�J�n V2="<<V2<<"("<<R*CON.get_magnet_B()<<"T)"<<endl;
	}
	////////////////////////////

	///���ʌv�Z�̏ꍇ�A��͗̈�̋��E���������R���E�����ɂ���B�i�����ŁAMPSTOFEM�̒i�K�ł��̂悤�ɂ�����FINE3D�����܂������Ȃ��j
	if(CON.get_EM_calc_type()==4)//���ʌv�Z
	{
		if(CON.get_uniform_B_sw()==OFF)//�ʏ�͂�����
		{
			for(int i=1;i<=nelm;i++)
			{
				if(ELEM[i].material==AIR)
				{
					for(int k=1;k<=4;k++)
					{
						int kelm=ELEM[i].elm[k];
						if(kelm==0)
						{
							int j=k%4+1;//ielm��jelm�̐ڂ���O�p�`���\������ߓ_�ԍ� 
							int m=4-(k-1)/2*2;
							int n=3-(k/2%2)*2;
							NODE[ELEM[i].node[j]].boundary_condition=0;//���E�����̖���
							NODE[ELEM[i].node[m]].boundary_condition=0;
							NODE[ELEM[i].node[n]].boundary_condition=0;
						}
					}
				}
			}
		}
		if(CON.get_uniform_B_sw()==ON)//��l����̏ꍇ
		{
			double ZU=CON.get_ZU();		//��͗̈��[
			double ZD=CON.get_ZD();		//��͗̈扺�[
			double err=1e-14;
			for(int i=1;i<=node;i++)
			{
				if(NODE[i].r[A_Z]>ZU-err) NODE[i].boundary_condition=2;//��[
				else if(NODE[i].r[A_Z]<ZD+err) NODE[i].boundary_condition=1;//���[
				else NODE[i].boundary_condition=0;
			}
		}
	}////////*/


    int NN=0;					//�f�B���N���^���E�ߓ_��
    int *dn=new int [node+1];	//�e�ߓ_���f�B���N���^���E�̉��Ԗڂ��B�f�B���N���^���E��łȂ��Ȃ�node+1���i�[
    double *PHAT=new double [CON.get_max_DN()];//�f�B���N���^���l
    
    ///�f�B���N���^���E��������
    for(int i=1;i<=node;i++)
    {
		if(NODE[i].boundary_condition==1)//��ŋ��E���������������Ă邩��A���̕���else if�ɂ��邱��
		{    
	        dn[i]=NN;
	        PHAT[NN]=V1;
	        V[i]=V1;
	        NN++;
		}
		else if(NODE[i].boundary_condition==2)
		{   
	        dn[i]=NN;
	        PHAT[NN]=V2;
	        V[i]=V2;
	        NN++;
		}
		else dn[i]=node+1;
    }/////////////*/
    cout<<"�f�B���N������"<<NN<<" ";
	    
    int pn=node-NN;///���m��
    int *ppn=new int [pn];		//�s���n�Ԗڂ̖��m���͐ߓ_�ԍ�ppn[n]�ɑ���
    int *npp=new int [node+1];	//i�Ԗڂ̐ߓ_�͍s���npp[i]�Ԗڂɑ����B�f�B���N���^�̏ꍇ�͍s��ɓ����Ă��Ȃ��̂�pn+1���i�[
    int num=0; 
    for(int i=1;i<=node;i++)
    {
        if(NODE[i].boundary_condition==0)	//���m���Ȃ�
		{
			ppn[num]=i;
			npp[i]=num;
			num++;
		}
		else npp[i]=pn+1;
    }
    
    ////�s��̍ő啝�v�Z
    int mat_w=0;
	mat_w=calc_matrix_width(CON,NODE,ELEM,node,nelm,jnb,nei);//������������������Ƃ߂�B�����ǎ኱�̌v�Z�R�X�g�ɂȂ�B
    //////
    
    ////�z��m��
    double **G=new double *[pn+1];						///�S�̍s��
    for(int i=1;i<=pn;i++) G[i]=new double [mat_w+1];
    int **ROW=new int *[pn+1];							///�e�s�́A��[���v�f�̗�ԍ��L��
    for(int i=1;i<=pn;i++) ROW[i]=new int [mat_w+1];
    int *NUM=new int [pn+1];							///�e�s�́A��[���v�f��
    
    for(int i=1;i<=pn;i++)//������
    {
        NUM[i]=0;
        for(int j=1;j<=mat_w;j++)
		{
			G[i][j]=0;
			ROW[i][j]=0;
		}
    }
    double *B=new double [pn];//���s��
    for(int i=0;i<pn;i++) B[i]=0;//������
    ////
	
	/*//���s��a�ɐߓ_�̓d�ׂ���
	if(CON.get_charge()==1 && CON.get_FEM_calc_type()==1)
	{
		for(int k=1;k<=WATER_N;k++)//���߂�WATER_N�̐ߓ_��TRANS[k]�̗��q�ɑ���
		{
			if(NODE[k].boundary_condition==0)//���m���Ȃ�B�ɂ��̐ߓ_�̗̈悪���邩����
			{
				int i=TRANS[k];//�ߓ_k�ɑ������闱�q�ԍ�
				for(int n=1;n<=jnb[k];n++)
				{
					int jelm=nei[k][n];//�ߓ_k�̗אڂ���v�f
					int N[5];
					for(int j=1;j<=4;j++) N[j]=ELEM[jelm].node[j];
		
					///�f���[�j�����̍ۂɋ��߂��̐ς́A�X�[�p�[�{�b�N�X���Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
					ELEM[jelm].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
					B[k]+=PART[i].angle/ep0/4*ELEM[jelm].volume;
					//if(PART[i].angle!=0)cout<<"EE"<<endl;
				}
				
			}
		}
	}///////////*/
	
	double *charge=new double [nelm+1];
	for(int n=1;n<=nelm;n++) charge[n]=0;
	
    /////////�S�̍s����쐬����
    
    int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
    for(int je=1;je<=nelm;je++)
    {   
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
		}
		
		///�f���[�j�����̍ۂɋ��߂��̐ς́A�X�[�p�[�{�b�N�X���Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
		
		double delta6=ELEM[je].volume;//�̐ς�6�{
		
		delta6=1/delta6;//�v�Z�ɕK�v�Ȃ̂͋t���Ȃ̂ł����Ŕ��]���Ă���
	
		double delta=ELEM[je].volume/6;//�{���̑̐�
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
			int m=j%4+1;
			int n=m%4+1;
	    
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]))*delta6;//i�������̂Ƃ�
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]))*delta6;//i�������̂Ƃ�
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]))*delta6;//i�������̂Ƃ�
			if(i%2!=0)//i����Ȃ�
			{
				c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
	
		/////��U�d�����`
		double ep=RP[je];
		/*if(ELEM[je].material==FLUID)
		{
			if(CON.get_EM_calc_type()==1) ep=CON.get_r_perm();
			else if(CON.get_EM_calc_type()==4) ep=CON.get_RP();//���ʌv�Z������䓧����
		}*/
	
		////�v�f�}�g���N�X�쐬�J�n
		for(int i=1;i<=4;i++)
		{
			if(NODE[N[i]].boundary_condition==0)///���m�Ȃ�
			{   
				int I=npp[N[i]]+1;///�ߓ_N[i]�͍s���I�Ԗ�
				for(int j=1;j<=4;j++)
				{					
					int J=npp[N[j]]+1;///�ߓ_N[j]�͍s���J�Ԗ�
					if(NODE[N[j]].boundary_condition==0)///���m�Ȃ�
					{   
						int flag=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=ep*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*delta;
								//if(I==85839 && h==1) cout<<h<<" "<<G[I][h]<<" "<<delta<<" "<<(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])<<endl;
								flag=1;
							}
						}
						if(flag==0)
						{   
							NUM[I]=NUM[I]+1;
							int H=NUM[I];
			    
							G[I][H]+=ep*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*delta;
							ROW[I][H]=J;
						}
					}
					else //N[j]���f�B���N���^���E�ߓ_�Ȃ�
					{
						int n=dn[N[j]];
						B[I-1]-=ep*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*delta*PHAT[n];
					}
				}
				//if(CON.get_charge()==ON) B[I-1]+=charge[je]/4*delta;
			}
		}   
    }
    ///////////////////////*/
    
	///���ۂ̍ő啝���v�Z
	double realwidth=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>realwidth) realwidth=NUM[i];
    
    int number=0;//�s��̔�[���v�f��
    for(int i=1;i<=pn;i++) number+=NUM[i];

    ///���̂܂܂ł�G[i][j]��[j]�̏��������ɕ���ł��Ȃ��̂ŁAG��ROW����ёւ���
	arrange_matrix(pn,NUM,ROW,G);

	//�Ώ̐��`�F�b�N
	check_matrix_symmetry(pn,NUM,ROW,G);

    double *val = new double [number];
    int *ind = new int [number];//��[���v�f�̗�ԍ��i�[�z��
    int *ptr = new int [pn+1];//�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
    
    /////////////////////val,ind ,ptr�ɒl���i�[
    int index=0;
    for(int n=0;n<pn;n++)
    {
        ptr[n]=index;
		for(int m=1;m<=NUM[n+1];m++)
		{
			val[index]=G[n+1][m];
			ind[index]=ROW[n+1][m]-1;
			index++;
		}
    }	    
    ptr[pn]=number;
    ////////////////////*/
    
    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
    
	cout<<"�s��쐬�I�� ��:"<<realwidth<<"/"<<mat_w<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
    
	///////////////////////�s��v�Z�J�n
	
	double *XX=new double [pn];//�s��̓����i�[
	if(CON.get_FEMCG()==1) ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON.get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON.get_FEMCG()==3) MRTR(CON,B,pn,XX,val,ind,ptr);
	else if(CON.get_FEMCG()==4) ICMRTR(CON,B,pn,XX,val,ind,ptr);
	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		V[i]=XX[n];//�d�ʃ�
	}
	delete [] XX;
	////////////////////////////
    
    ///////�d�ʂ��t�@�C���o��
	ofstream fp("V.dat");
    for(int i=1;i<=node;i++)
    {
        if(NODE[i].r[A_Y]<2*le && NODE[i].r[A_Y]>-2*le) fp<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Z]<<" "<<V[i]<<endl;
    }
	fp.close(); 
    
    ////////////////////////

    delete [] dn;
    delete [] PHAT;
    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;

	delete [] charge;
}

//�s��̕��v�Z�֐�
int calc_matrix_width(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,int *jnb,int **nei)
{
	//�W���s��̂Ȃ��ŁAi�s�ڂɊ܂܂��񂺂�v�f���́A�v�Z�_i����̂т�ӂ̐�(���Ȃ킿�ߗ׌v�Z�_��)+1(���g)�ł���B����̍ő�l�����߂�
	///�z��m��
		
	int *temp_check=new int[node+1];	//�ꎞ�I�������z��
		
	///������
	for(int i=1;i<=node;i++) temp_check[i]=0;
	/////////////////////

	int maxwidth=0;

	//if(CON.get_FEM_calc_type()==1 || CON.get_FEM_calc_type()==4 || CON.get_iteration_count()==0)
	{
		for(int i=1;i<=node;i++)
		{
			if(NODE[i].boundary_condition==0)
			{
				int width=0;
				vector<int> NEI2;//�e�ߓ_�̗אڂ���ߓ_�ԍ��i�[(nei�͐ߓ_-�v�f�ł��邪�ANEI2�͐ߓ_-�ߓ_)
				for(int k=1;k<=jnb[i];k++)
				{
					int jelm=nei[i][k];//�ߓ_i�ɗאڂ���v�f�ԍ�
					for(int j=1;j<=4;j++)
					{
						int p=ELEM[jelm].node[j];
						if(p!=i && temp_check[p]==0 && NODE[p].boundary_condition==0)
						{	
							width++;
							NEI2.push_back(p);
							temp_check[p]=1;//��������
						}
					}
				}
				for(size_t k=0;k<NEI2.size();k++) temp_check[NEI2[k]]=0;//������
				if(width>maxwidth) maxwidth=width;
			}
		}
	}

	delete [] temp_check;

	return maxwidth+1;//width�͎�������L�т�ӂ̐��ɓ������B�s��̕��͎������g�������width+1
}

//�s����ёւ��֐�
void arrange_matrix(int pn,int *NUM,int **ROW,double **G)
{
	///G[i][j]��[j]�̏��������ɕ���ł��Ȃ��̂ŁAG��ROW����ёւ���
	
    for(int i=1;i<=pn;i++)
    {
        double tempG;
		int tempR;
		for(int j=2;j<=NUM[i];j++)
		{
			for(int m=1;m<j;m++)
			{
				if(ROW[i][j]<ROW[i][m])
				{
					tempG=G[i][m];
					tempR=ROW[i][m];
					G[i][m]=G[i][j];


					ROW[i][m]=ROW[i][j];
					G[i][j]=tempG;
					ROW[i][j]=tempR;
				}
			}
		}
    }///////////
}

//ICCG�@
void ICCG3D2(mpsconfig &CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X)
{
	//val :�[���v�f�̒l
	//ind:��[���v�f�̗�ԍ��i�[�z��
	//ptr:�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
	//X[n]:��

	int count;						//�����グ�ϐ�
	double accel=CON.get_CGaccl();	//�����t�@�N�^
	
	int num2=0;//�Ίp�������܂ށA���O�p�s�񂾂����l���ɂ��ꂽ��[���v�f��
	for(int k=0;k<pn;k++) for(int m=ptr[k];m<ptr[k+1];m++) if(ind[m]<=k) num2++;	
	if(num2!=(number-pn)/2+pn) cout<<"ERROR"<<endl;
	
	double *matrix=new double[num2];//�W���s���ۑ�(�Ίp�������܂މ��O�p�s��) ��[���v�f������1���z��Ƃ��ĕۑ�
	int *ind2 = new int [num2];
	int *ptr2 = new int [pn+1];

	num2=0;
	for(int k=0;k<pn;k++)
	{	
		ptr2[k]=num2;
	    for(int m=ptr[k];m<ptr[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
			if(ind[m]<=k)
			{
				matrix[num2]=val[m];
				ind2[num2]=ind[m];
				num2++;
			}
		}
	}
	ptr2[pn]=num2;//��������Ă����Ȃ��ƁA�Ō��(int m=ptr2[k];m<ptr2[k+1];m++)�݂����Ȃ��Ƃ��ł��Ȃ�

	int *NUM = new int [pn];			//������ɂ݂��A�e��̗v�f��
	for(int k=0;k<pn;k++) NUM[k]=0;
	
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
			int J=ind2[m];
			NUM[J]=NUM[J]+1;
		}
	}
	double **VAL=new double *[pn];						//�[���v�f�̒l VAL[i][k]��i���k�Ԗڂ̔�[���v�f
	for(int i=0;i<pn;i++) VAL[i]=new double [NUM[i]];
	int **IND = new int *[pn];							//��[���v�f�̍s�ԍ��i�[�z��
	for(int i=0;i<pn;i++) IND[i]=new int [NUM[i]];
	
	/////////////////////////////////////iccg�@
	double alp,beta;
	double rLDLt_r;
	double E=1;//�덷
    double *r=new double[pn];
	double *AP = new double [pn];
	double *P = new double [pn];
	double *y=new double [pn];
	double *LDLt_r= new double [pn];
	double *L=new double[num2];//�s���S�R���X�L�[������̉��O�p�s��L�i�[
	double *D1 = new double [pn];//D�s��
	
	/////�s���S�R���X�L�[����
	Incomplete_Cholesky_Decomposition(CON,L,D1,matrix,ptr2,ind2,pn,num2);//L��D1�ɒl���������܂��
	
	delete [] matrix;

	///�����ɂ����z��ɒl����
	for(int k=0;k<pn;k++) NUM[k]=0;
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
			int J=ind2[m];
			VAL[J][NUM[J]]=L[m];
			IND[J][NUM[J]]=k;
			NUM[J]=NUM[J]+1;
		}
	}////////*/

	for(int n=0;n<pn;n++) //�����l
	{
		X[n]=0;
		r[n]=B[n];
	}

	/////////////////y[i]�����Ƃ߁ALDLt_r[i]�����Ƃ߂�B
	for(int i=0;i<pn;i++)
	{
		if(i==0) y[0]=r[0]/L[0]; //���i3.77�j 
		else
		{
		    double sum=0;
			for(int m=ptr2[i];m<ptr2[i+1]-1;m++) sum+=L[m]*y[ind2[m]];//���i3.78�j
		    int m=ptr2[i+1]-1;
			y[i]=(r[i]-sum)/L[m];
		}
	}////y[i]�����Ƃ܂����B
	for(int i=pn-1;i>=0;i--)
	{
	    double sum=0;
		for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
	    LDLt_r[i]=y[i]-D1[i]*sum;	
	}
	/////////////////*/
	
	for(int n=0;n<pn;n++) P[n]=LDLt_r[n];

    cout<<"ICCG�@�X�^�[�g-----���m��="<<pn<<"  ---";
	unsigned int time=GetTickCount();
	count=0;
	double ep=CON.get_EMCGep();//��������
	while(E>ep)
	{
		count++;
		if(count==pn) cout<<"count=pn"<<endl;
		//////////////alp�����߂�
		rLDLt_r=0;

		for(int n=0;n<pn;n++) rLDLt_r+=r[n]*LDLt_r[n];
		
		for(int n=0;n<pn;n++)
		{
			AP[n]=0;
			for(int m=ptr[n];m<ptr[n+1];m++) AP[n]+=val[m]*P[ind[m]];
		}
		double PAP=0;
		for(int n=0;n<pn;n++)  PAP+=P[n]*AP[n];
		alp=rLDLt_r/PAP;
		//cout<<"alp="<<alp<<" PAP="<<PAP<<endl;
		//////////////////////
	
		//////////////// X(k+1)=X(k)+alp*P
		for(int n=0;n<pn;n++) X[n]+=alp*P[n];
		
		//////////////// r=r-alp*AP
		for(int n=0;n<pn;n++) r[n]-=alp*AP[n];
		/////////////////////////////
		
		//////////////////�덷
		E=0;
		for(int n=0;n<pn;n++) E+=r[n]*r[n];
		
		E=sqrt(E);
		//cout<<"E="<<E<<endl;
		////////////////////////
		
		///////////////////////beta
		beta=1.0/rLDLt_r;
		rLDLt_r=0;
		
        /////////////////y[i]�����Ƃ߁ALDLt_r[i]�����Ƃ߂�B
		for(int i=0;i<pn;i++)
		{
			if(i==0) y[0]=r[0]/L[0]; //���i3.77�j �V
			else
			{
			    double sum=0;
			    for(int m=ptr2[i];m<ptr2[i+1]-1;m++)//�Ίp�����͏�������ptr[i+1]-1
			    {
					sum+=L[m]*y[ind2[m]];//���i3.78�j
			    }
			    int m=ptr2[i+1]-1;
				y[i]=(r[i]-sum)/L[m];
			}
		}////y[i]�����Ƃ܂����B
	
		/////////LDLt_r[i]�����߂�
		for(int i=pn-1;i>=0;i--)
		{
		    double sum=0;
			for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
			
		    LDLt_r[i]=y[i]-D1[i]*sum;	
		}
		/////////////////*/
	
		for(int n=0;n<pn;n++) rLDLt_r+=r[n]*LDLt_r[n];
		
		beta=beta*rLDLt_r;
		/////////////////*/
		
		///////////////////// P=r+beta*P
		for(int n=0;n<pn;n++) P[n]=LDLt_r[n]+beta*P[n];//iccg
	}
	cout<<"������="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;
	
	
    delete [] r;
	delete [] AP;
	delete [] P;

	delete [] y;
	delete [] LDLt_r;
	delete [] D1;

	delete [] ind2;
	delete [] ptr2;

	for(int i=0;i<pn;i++)
	{
		delete [] VAL[i];
		delete [] IND[i];
	}
	delete [] VAL;
	delete [] IND;
	delete [] NUM;

	delete [] L;
}



//�s���S�R���X�L�[����
void Incomplete_Cholesky_Decomposition(mpsconfig &CON,double *L,double *D1,double *matrix,int *ptr2,int *ind2,int pn,int num2)
{
	//matrix[i]:�W���s��̂����A�Ίp���܂މ��O�p���i�[����Ă���
	//ptr2[k]:matrix[i]�ɑ΂���A�N�Z�X�z��
	//L[i],D[i]:�s���S�R���X�L�[�����i�[�z��
	//pn:���m��
	//num2:matrix[i]�̗v�f��

	cout<<"�s���S�R���X�L�[�����J�n---";
	unsigned int timeC=GetTickCount();
	double accel=CON.get_CGaccl();		//�����W��
	double accel2=0;					//�����v�Z���ꂽ�����W��
	double maxvalue=0;					//matrix���̍ő�l
	int Line=1;							//�ő�l�����s�ԍ�
	int Column=1;						//�ő�l������ԍ�
	
	//1���:�܂��͕��ʂɕs���S�R���X�L�[�������s���B������D1[i]<0�Ȃ�����W�����C�����A2��ڂ��s��
	int onemoreflag=OFF;
	L[0]=matrix[0]*accel;
	D1[0]=1/L[0];			//L[0]�͌W���s��̒l�ƕς��Ȃ��̂ŁA�[���ł��Ȃ���Ε��ł��Ȃ��Ɨ\�z�����
	if(D1[0]<0) cout<<"error in �s���S�R���X�L�[���� D1[0]<0"<<endl;
	maxvalue=D1[0]*D1[0];
	for(int k=1;k<pn;k++)	//k=0�͂����ς񂾁Bk>=1����loop
	{	
		for(int m=ptr2[k];m<ptr2[k+1]-1;m++)///0����k-1�܂ł̗v�f
		{
			int i=ind2[m];//��ԍ�
		
			double sum=0;
			for(int j=ptr2[k];j<m;j++) for(int J=ptr2[i];J<ptr2[i+1];J++) if(ind2[J]==ind2[j]) sum+=L[j]*L[J]*D1[ind2[j]];//i�s�ڂ̂Ȃ������̈�v������̂�T���Ă���B������Ԃ��H
			L[m]=matrix[m]-sum;//i=0�̂Ƃ���sum=0�ƂȂ莩���I��L[m]=matrix[m]�ƂȂ�
			if(L[m]*L[m]>maxvalue)
			{
				maxvalue=L[m]*L[m];
				Line=k; Column=i;
			}
		}
		int m=ptr2[k+1]-1;				//����m�͕K���Ίp�����ɊY������
				
		double akk=matrix[m];			//�W���s��̑Ίp����.�l��ۑ�.
		double sum=0;
		for(int j=ptr2[k];j<m;j++) sum+=L[j]*L[j]*D1[ind2[j]];
		L[m]=matrix[m]*accel-sum;		//�����W����������
		if(L[m]*L[m]<1e-8)L[m]=0.0001;	//Lkk�����ɏ������Ƃ��A���̂悤�ɂ���B�����ŁA�ǂ݂̂�Lkk<0�Ȃ�������Â炢����ALkk�����̋ɏ��l�ł����̒l�ɂ��Ă��܂��Ă悢�̂ł́H
		D1[k]=1/L[m];
		
		if(L[m]<0)						//L[m]<0�̂Ƃ�D1[m]<0�ƂȂ�A�������ɒ[�ɒx���Ȃ�.�����h�����߂ɉ����W����傫������
		{
			accel2=sum/akk*1.1;
			if(accel2>accel) accel=accel2;
			onemoreflag=ON;				//���̃X�C�b�`��ON�Ȃ�A�s���S�R���X�L�[������������x�s���B
		}
		if(L[m]*L[m]>maxvalue)
		{
			maxvalue=L[m]*L[m];
			Line=k; Column=k;
		}
	}
	if(Line!=Column) cout<<"value="<<maxvalue<<" Line="<<Line<<" Column="<<Column<<endl;
	if(onemoreflag==ON)
	{
		L[0]=matrix[0]*accel;
		D1[0]=1/L[0];			//L[0]�͌W���s��̒l�ƕς��Ȃ��̂ŁA�[���ł��Ȃ���Ε��ł��Ȃ��Ɨ\�z�����
		for(int k=1;k<pn;k++)	//k=0�͂����ς񂾁Bk>=1����loop
		{	
			for(int m=ptr2[k];m<ptr2[k+1]-1;m++)///0����k-1�܂ł̗v�f
			{
				int i=ind2[m];//��ԍ�
		
				double sum=0;		//accel�̕ύX�ɂ��L��D1�̒l���ς�����̂ŁAsum���v�Z���Ȃ����B
				for(int j=ptr2[k];j<m;j++) for(int J=ptr2[i];J<ptr2[i+1];J++) if(ind2[J]==ind2[j]) sum+=L[j]*L[J]*D1[ind2[j]];//i�s�ڂ̂Ȃ������̈�v������̂�T���Ă���B������Ԃ��H
				L[m]=matrix[m]-sum;//i=0�̂Ƃ���sum=0�ƂȂ莩���I��L[m]=matrix[m]�ƂȂ�
			}
			int m=ptr2[k+1]-1;				//����m�͕K���Ίp�����ɊY������
				
			double akk=matrix[m];			//�W���s��̑Ίp����.�l��ۑ�.
			double sum=0;
			for(int j=ptr2[k];j<m;j++) sum+=L[j]*L[j]*D1[ind2[j]];
			L[m]=matrix[m]*accel-sum;		//�����W����������
			if(L[m]*L[m]<1e-8)L[m]=0.0001;	//Lkk�����ɏ������Ƃ��A���̂悤�ɂ���B�����ŁA�ǂ݂̂�Lkk<0�Ȃ�������Â炢����ALkk�����̋ɏ��l�ł����̒l�ɂ��Ă��܂��Ă悢�̂ł́H
			D1[k]=1/L[m];
		}
	}	
	cout<<"�����W��="<<accel<<" time="<<(GetTickCount()-timeC)*0.001<<"[sec]"<<endl;
	///�s���S�R���X�L�[��������/////////*/
}

//�s��̑Ώ̐������֐�
void check_matrix_symmetry(int pn,int *NUM,int **ROW,double **G)
{
	ofstream sym("symmetry_out.csv");
	for(int i=1;i<=pn;i++)
	{
		int flag_sym=0;
		for(int j=1;j<=NUM[i];j++)
		{
			int J=ROW[i][j];
			int flag=0;
			for(int k=1;k<=NUM[J];k++)
			{
				if(ROW[J][k]==i)
				{
					flag=1;
					if(G[i][j]!=G[J][k])
					{
						flag_sym=ON;
						cout<<"�Ώ̐��G���[ "<<i<<" "<<G[i][j]<<" "<<G[J][k]<<endl;							
					}
				}
			}
			if(flag==0)cout<<"DDD"<<endl;
		}
		if(flag_sym==ON)	sym<<i<<endl;
	}
	sym.close();
}

//�d�E�v�Z�֐�
void ELECTRO3D(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double *V,double **Ee)
{
    //�@Ee:�v�f�̓d�E or �v�f�̎��EB
    
	if(CON.get_EM_calc_type()==1) cout<<"�d�E�v�Z-------";
	if(CON.get_EM_calc_type()==4) cout<<"���E�v�Z(����)------";
	unsigned timeA=GetTickCount();		//�v�Z�J�n����
	double ep0=8.854e-12;				//�^��̗U�d��
	double u0=12.56e-7;					//�^��̓�����
    
    int N[4+1];							//�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];

	//double maxE=0;
	//double Xc[3];
    ///�v�f�̓d�E�E���E�����Ƃ߂�
    for(int je=1;je<=nelm;je++)
    {
        for(int j=1;j<=4;j++)
		{
			N[j]=ELEM[je].node[j];
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
		}
	
		double delta6=ELEM[je].volume;//�̐ς�6�{  �̐ς͓d�ʋ��߂�Ƃ��Ɍv�Z���Ă���
	
		delta6=1/delta6;//�v�Z�ɕK�v�Ȃ̂͋t���Ȃ̂ł����Ŕ��]���Ă���
		
		///�W��c,d,e�v�Z
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
			int m=j%4+1;
			int n=m%4+1;
	    
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]))*delta6;//i�������̂Ƃ�
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]))*delta6;//i�������̂Ƃ�
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]))*delta6;//i�������̂Ƃ�
			if(i%2!=0)//i����Ȃ�
			{
				c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}   
		}
		/////////
    
		Ee[A_X][je]=c[1]*V[N[1]]+c[2]*V[N[2]]+c[3]*V[N[3]]+c[4]*V[N[4]];
        Ee[A_Y][je]=d[1]*V[N[1]]+d[2]*V[N[2]]+d[3]*V[N[3]]+d[4]*V[N[4]];
		Ee[A_Z][je]=e[1]*V[N[1]]+e[2]*V[N[2]]+e[3]*V[N[3]]+e[4]*V[N[4]];
	
		for(int D=0;D<3;D++) Ee[D][je]*=-1;//�d�E�͓d�ʌ��z�̃}�C�i�X

		double EE=sqrt(Ee[A_X][je]*Ee[A_X][je]+Ee[A_Y][je]*Ee[A_Y][je]+Ee[A_Z][je]*Ee[A_Z][je]);
		/*if(EE>maxE) 
		{
			maxE=EE;
			Xc[A_X]=(X[1]+X[2]+X[3]+X[4])*0.25;
			Xc[A_Y]=(Y[1]+Y[2]+Y[3]+Y[4])*0.25;
			Xc[A_Z]=(Z[1]+Z[2]+Z[3]+Z[4])*0.25;
		}*/
    }///�d�E�v�Z�I��

	//cout<<"MAX="<<maxE<<" x,y,z="<<Xc[A_X]<<" "<<Xc[A_Y]<<" "<<Xc[A_Z]<<endl;
    
	/*if(CON.get_FEM_calc_type()==4)//���ʂɂ�鎥���͂̏ꍇ
	{
		///E[D][je]�Ɋi�[����Ă���̂�H�̒l�Ȃ̂ŁA����ɓ�������������
		
		for(int je=1;je<=nelm;je++)
		{
			if(ELEM[je].material==WATER) for(int D=0;D<3;D++) Ee[D][je]*=u0*CON.get_RP();
			else for(int D=0;D<3;D++) Ee[D][je]*=u0;
		}
	}///////////////////////*/

    /////�t�@�C���o��
	if(CON.get_EM_calc_type()==1 && CON.get_E_times()>0)
	{
		ofstream fp("E.dat");
		double times=CON.get_E_times();
		double le=CON.get_distancebp();
		double Xmin=CON.get_XL()+le; double Xmax=CON.get_XR()-le;//��͗̈�
		double Zmin=CON.get_ZD()+le; double Zmax=CON.get_ZU()-le;
		double dx=3*le;
		int Nx=(int)((Xmax-Xmin)/dx);//�e�����̕�����
		int Nz=(int)((Zmax-Zmin)/dx);
		int search=nelm;//locate�֐��ōŏ��ɒT������v�f�ԍ�

		if(CON.get_plot_E_type()==1)
		{
			for(int n=0;n<Nx;n++)
			{
				for(int m=0;m<Nz;m++)
				{
					double xp=dx*n+Xmin;//�o�͂���_�̍��W
					double yp=0;
					double zp=dx*m+Zmin;
					int loc=locate3D(NODE,ELEM,search,xp,yp,zp);

					fp<<xp<<" "<<zp<<" "<<-Ee[A_X][loc]*times<<" "<<-Ee[A_Z][loc]*times<<endl;///�����炢����킴�Ɓ{�|���]���ďo��
					search=loc;
				}
			}
		}
		else if(CON.get_plot_E_type()==2)//�X�J���[
		{
			for(int n=0;n<Nx;n++)
			{
				for(int m=0;m<Nz;m++)
				{
					double xp=dx*n+Xmin;//�o�͂���_�̍��W
					double yp=0;
					double zp=dx*m+Zmin;
					int loc=locate3D(NODE,ELEM,search,xp,yp,zp);

					fp<<xp<<" "<<zp<<" "<<sqrt(Ee[A_X][loc]*Ee[A_X][loc]+Ee[A_Z][loc]*Ee[A_Z][loc])<<endl;
					
					//fp<<xp<<" "<<zp<<" "<<sqrt(Ee[A_X][loc]*Ee[A_X][loc]+Ee[A_Y][loc]*Ee[A_Y][loc]+Ee[A_Z][loc]*Ee[A_Z][loc])<<endl;
					search=loc;
				}
			}
		}
		fp.close();
	}
	if(CON.get_EM_calc_type()==4)
	{
		ofstream fp("H.dat");
		double times=CON.get_B_times()*u0;//u0: ������
		double Le=CON.get_distancebp();
		int plot_type=CON.get_plot_B_type();	//�o�̓^�C�v�@1=�x�N�g���@2=�X�J���[
		for(int i=1;i<=nelm;i++)
		{
			if(ELEM[i].material!=MAGNET)
			{
				double X=0;
				double Y=0;
				double Z=0;
				for(int j=1;j<=4;j++)
				{
					X+=NODE[ELEM[i].node[j]].r[A_X]*0.25;
					Y+=NODE[ELEM[i].node[j]].r[A_Y]*0.25;
					Z+=NODE[ELEM[i].node[j]].r[A_Z]*0.25;
				}
				
				if(Y>-2*Le && Y<2*Le)
				{
					if(plot_type==1) fp<<X<<" "<<Z<<" "<<Ee[A_X][i]*times<<" "<<Ee[A_Z][i]*times<<endl;
					if(plot_type==2) fp<<X<<" "<<Z<<" "<<sqrt(Ee[A_X][i]*Ee[A_X][i]+Ee[A_Y][i]*Ee[A_Y][i]+Ee[A_Z][i]*Ee[A_Z][i])<<endl;
				}
			}
		}
		fp.close();
	}
    ////////*/
	cout<<"ok time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
}

//�Î���v�Z�֐�
void calc_static_magnetic_field(mpsconfig &CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,int *jnb,double dt,double TIME,int t,int **nei,int KTE,vector<mpselastic> &PART,int fluid_number,double **F,int KTJ)
{
	Post post;
	sides *SIDE=new sides[KTE];			//�ӃN���X
	int *branch_num=new int[node+1];	//�e�ߓ_���אڂ���ߓ_��(�e�ߓ_���牄�т�ӂ̐�)
	int max=1000;						//�ߓ_�ɗאڂ���ő�ߓ_���E�E�Evector�ɂ���ׂ�
	int **nei2=new int* [node+1];		//�e�ߓ_�̗אڂ���ߓ_�ԍ��i�[(nei�͐ߓ_-�v�f�ł��邪�Anei2�͐ߓ_-�ߓ_)
	for(int i=1;i<=node;i++) nei2[i]=new int [max];
	int side_num=0;						//�S�Ӑ��i�[
		
//�ӗv�f���� (�ߓ_�v�f���g�p����ꍇ�ł��A�d�����x�����߂�Ƃ��ɕӗv�f���ق���)
//side��edge�ɗp��𓝈ꂷ�邱��
	side_num=make_edge_element(CON,NODE,ELEM,node,nelm,jnb,nei,SIDE,branch_num,nei2,KTE);

	double *current[3];//�e�v�f�̓d�����x[A/m^2]
	for(int D=0;D<3;D++) current[D]=new double [nelm+1];

	double *Be[3];//�v�f���d�E or �v�f�����E
    for(int D=0;D<3;D++) Be[D]=new double [nelm+1];
	double *RP=new double [nelm+1];	//�e�v�f�̓�����

//����������
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==FLUID || ELEM[i].material==MAGELAST || ELEM[i].material==MAGELAST2) RP[i]=CON.get_RP();
		else if(ELEM[i].material==ELASTIC) RP[i]=1;
		else if(ELEM[i].material==AIR) RP[i]=1;
		else if(ELEM[i].material==COIL) RP[i]=1;
		else if(ELEM[i].material==MAGNET) RP[i]=1;
		else if(ELEM[i].material==IRON) RP[i]=2500;
		else cout<<"error:�ގ��ɑ΂��A���������s�� �ގ�="<<ELEM[i].material<<endl;
	}

//�d�����x����
	if(CON.get_J_input_way()==2)//�d�����x�𑼂̃\�t�g����ǂݍ���
	{
		import_J0_density(CON, node, nelm,NODE,ELEM,current);
		check_J0(CON, node, nelm,NODE,ELEM,current);
	}
	else if(CON.get_J_input_way()==1)
	{
		//������
	}

	else if(CON.get_J_input_way()==0) for(int i=0;i<=nelm;i++) for(int D=0;D<3;D++) current[D][i]=0.0;//�����d�������݂��Ȃ����珉����

	if(CON.get_FEM_elm_type()==1)//�ӗv�f
	{
		double *A=new double [side_num+1];//�x�N�g���|�e���V����
		
		//�x�N�g���|�e���V�������v�Z����
		if(CON.get_EM_calc_type()==2) Avector3D(CON,NODE,ELEM,SIDE,node,nelm,side_num,A,jnb,branch_num,current,RP);	

		//�������x���v�Z����
		Bflux3D_side(CON,NODE,ELEM,SIDE,node,nelm,side_num,A,Be,RP);

		//�������x���v���b�g
		if(post.get_plot_B()==true){
			plot_magnetic_flux_density(CON, PART, NODE, ELEM, nelm, Be, t);
			C_Fluix_data_avs(node,nelm,NODE,ELEM,KTJ,Be,CON,t);
		}
		
		//�͂̌v�Z
		if(CON.get_m_force()==1) NODE_F3D(CON,NODE,ELEM,node,nelm,Be,jnb,nei,RP,PART,F,fluid_number);
		if(CON.get_m_force()==2) kelvin_force3D(CON,PART,NODE,ELEM, node, nelm,Be, fluid_number,F,RP,jnb,nei);
		if(CON.get_m_force()==4) direct_divT3D(CON,PART,NODE,ELEM, node, nelm,Be, fluid_number,F,RP,jnb,nei);

		//�X���[�W���O�Ɨ͂̃v���b�g
		if(post.get_plot_NF()==true){
			smoothingF3D(CON,PART,fluid_number,F,t);
		}

		delete [] A;
	}

	delete [] RP;
	for(int D=0;D<3;D++)
	{
		delete [] current[D];
		delete [] Be[D];
	}
	for(int i=1;i<=node;i++) delete [] nei2[i];
	delete [] nei2;
	delete [] branch_num;
	delete [] SIDE;
}

//�Î���v�Z�֐�
void calc_variable_magnetic_field(mpsconfig &CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,int *jnb,double dt,double TIME,int t,int **nei,int KTE,vector<mpselastic> &PART,int fluid_number,double **F)
{
	Post post;
	sides *SIDE=new sides[KTE];			//�ӃN���X
	int *branch_num=new int[node+1];	//�e�ߓ_���אڂ���ߓ_��(�e�ߓ_���牄�т�ӂ̐�)
	int max=1000;						//�ߓ_�ɗאڂ���ő�ߓ_���E�E�Evector�ɂ���ׂ�
	int **nei2=new int* [node+1];		//�e�ߓ_�̗אڂ���ߓ_�ԍ��i�[(nei�͐ߓ_-�v�f�ł��邪�Anei2�͐ߓ_-�ߓ_)
	for(int i=1;i<=node;i++) nei2[i]=new int [max];
	int side_num=0;						//�S�Ӑ��i�[
		
//�ӗv�f���� (�ߓ_�v�f���g�p����ꍇ�ł��A�d�����x�����߂�Ƃ��ɕӗv�f���ق���)
//side��edge�ɗp��𓝈ꂷ�邱��
	side_num=make_edge_element(CON,NODE,ELEM,node,nelm,jnb,nei,SIDE,branch_num,nei2,KTE);

	double *current[3];//�e�v�f�̓d�����x[A/m^2]
	for(int D=0;D<3;D++) current[D]=new double [nelm+1];

	double *Be[3];//�v�f���d�E or �v�f�����E
    for(int D=0;D<3;D++) Be[D]=new double [nelm+1];
	double *RP=new double [nelm+1];	//�e�v�f�̓�����

//����������
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==FLUID || ELEM[i].material==MAGELAST || ELEM[i].material==MAGELAST2) RP[i]=CON.get_RP();
		else if(ELEM[i].material==ELASTIC) RP[i]=1;
		else if(ELEM[i].material==AIR) RP[i]=1;
		else if(ELEM[i].material==COIL) RP[i]=1;
		else if(ELEM[i].material==MAGNET) RP[i]=1;
		else if(ELEM[i].material==IRON) RP[i]=2500;
		else cout<<"error:�ގ��ɑ΂��A���������s�� �v�ǉ�"<<endl;
	}

//�d�����x����
	if(CON.get_J_input_way()==2)//�d�����x�𑼂̃\�t�g����ǂݍ���
	{
		import_J0_density(CON, node, nelm,NODE,ELEM,current);
		check_J0(CON, node, nelm,NODE,ELEM,current);
	}
	else if(CON.get_J_input_way()==1)
	{
		//������
	}

	else if(CON.get_J_input_way()==0) for(int i=0;i<=nelm;i++) for(int D=0;D<3;D++) current[D][i]=0.0;//�����d�������݂��Ȃ����珉����

	if(CON.get_FEM_elm_type()==1)//�ӗv�f
	{
		double *A=new double [side_num+1];//�x�N�g���|�e���V����
		
		//�x�N�g���|�e���V�������v�Z����
		if(CON.get_EM_calc_type()==2) Avector3D(CON,NODE,ELEM,SIDE,node,nelm,side_num,A,jnb,branch_num,current,RP);

		//�������x���v�Z����
		Bflux3D_side(CON,NODE,ELEM,SIDE,node,nelm,side_num,A,Be,RP);

		//�������x���v���b�g
		if(post.get_plot_B()==true){
			plot_magnetic_flux_density(CON, PART, NODE, ELEM, nelm, Be, t);
		}

		//�͂̌v�Z
		if(CON.get_m_force()==1) NODE_F3D(CON,NODE,ELEM,node,nelm,Be,jnb,nei,RP,PART,F,fluid_number);
		if(CON.get_m_force()==2) kelvin_force3D(CON,PART,NODE,ELEM, node, nelm,Be, fluid_number,F,RP,jnb,nei);
		if(CON.get_m_force()==4) direct_divT3D(CON,PART,NODE,ELEM, node, nelm,Be, fluid_number,F,RP,jnb,nei);

		//�X���[�W���O�Ɨ͂̃v���b�g
		if(post.get_plot_NF()==true){
			smoothingF3D(CON,PART,fluid_number,F,t);
		}

		delete [] A;
	}

	delete [] RP;
	for(int D=0;D<3;D++)
	{
		delete [] current[D];
		delete [] Be[D];
	}
	for(int i=1;i<=node;i++) delete [] nei2[i];
	delete [] nei2;
	delete [] branch_num;
	delete [] SIDE;
}


//�d���̓X���[�W���O�֐��E�E�E�e�����a���ɂ�����ӗ��q�̓d���͂̕��ϒl��^���Ă���
void smoothingF3D(mpsconfig &CON, vector<mpselastic> &PART,int fluid_number,double *F[3],int t)
{

	double le=CON.get_distancebp();
    double ep=8.854e-12;///�^��̗U�d���B
    double *newF[3];
    for(int D=0;D<3;D++) newF[D]=new double [PART.size()];

	if(CON.get_FEM_smn()>0)
	{
		for(int n=0;n<CON.get_FEM_smn();n++)//������P�ȏ�ɂ���K�v����H�H�H
		{
		    for(size_t i=0;i<PART.size();i++) 
		    {  
				if(PART[i].type!=WALL && PART[i].type==MAGELAST)
				{
					for(int D=0;D<3;D++) newF[D][i]=F[D][i];
					int num=1; //�������g���J�E���g���邩��1
					for(int k=0;k<PART[i].N;k++)
					{       
						int j=PART[i].NEI[k];
						if(PART[j].type==FLUID || PART[j].type==MAGELAST)
						{
							num++;
							for(int D=0;D<3;D++) newF[D][i]+=F[D][j];
						}
					}
					for(int D=0;D<3;D++) newF[D][i]/=num;
				}			
			} 
			for(size_t i=0;i<PART.size();i++)
			{
				if(PART[i].type!=WALL && PART[i].type==MAGELAST)	for(int D=0;D<3;D++) F[D][i]=newF[D][i];
			}
		}
	}
	else if(CON.get_FEM_smn()<0)//�\�ʂ݂̂ŃX���[�W���O
	{
		int N=-1*CON.get_FEM_smn();
		for(int n=0;n<N;n++)
		{
			for(size_t i=0;i<PART.size();i++) 
			{  
				if(PART[i].type!=WALL)
				{
					for(int D=0;D<3;D++) newF[D][i]=F[D][i];
					if(PART[i].surface==ON)
					{
						int num=1; //�������g���J�E���g���邩��1
						for(int k=0;k<PART[i].N;k++)
						{       
							int j=PART[i].NEI[k];
							if(PART[j].surface==ON && (PART[j].type==FLUID || PART[j].type==MAGELAST))
							{
								num++;
								for(int D=0;D<3;D++) newF[D][i]+=F[D][j];
							}
						}
						for(int D=0;D<3;D++) newF[D][i]/=num;
					}
				}
			} 
			for(size_t i=0;i<PART.size();i++){
				if(PART[i].type!=WALL){
				for(int D=0;D<3;D++) F[D][i]=newF[D][i];
				}
			}
		}
	}

    for(int D=0;D<3;D++) delete [] newF[D];
    /////////////////////////////////////*/
    
    //�d���͂��v���b�g, GNUplot�p

	ofstream fp("F.dat");
	double xmax=-100;//�o�͗��q�̍ő剡���W
	double ymax=-100;//�o�͗��q�̍ő�c���W
	double mass=CON.get_particle_mass();
	double times=CON.get_times()/mass*le*le*CON.get_FEMtimes();
    for(size_t i=0;i<PART.size();i++)
    {
		if(PART[i].type!=WALL){
		if(PART[i].r[A_Y]>-le*0.5 && PART[i].r[A_Y]<+le*0.5)
		{
			fp<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<F[A_X][i]*times<<"\t"<<F[A_Z][i]*times<<endl;
			if(PART[i].r[A_X]>xmax) xmax=PART[i].r[A_X];
			if(PART[i].r[A_Z]>ymax) ymax=PART[i].r[A_Z];
		}
		}
    }
	xmax+=4*le;//�}����o���ʒu��ی��������ď����΂߂Ɉړ�
	ymax+=4*le;
	if(CON.get_legend_F()>0) fp<<xmax<<" "<<ymax<<" "<<CON.get_legend_F()*times<<" "<<0*times<<endl;//�Ō�ɖ}��o��
	fp.close();

	int BOnum=0;
	for(size_t i=0;i<PART.size();i++) if(PART[i].surface==ON && PART[i].type!=WALL) BOnum++;

	//�\���i�q�f�[�^�t�@�C��
	int timestep=CON.get_current_step();
	stringstream sstr;
	sstr<<"./Lorentz/Lorentz"<<timestep<<".fld";
	string NodalForce=sstr.str();

	ofstream fout2(NodalForce);
	if(fout2.fail()){
		system("mkdir Lorentz");
		ofstream fout2(NodalForce);
		if(fout2.fail()){
			cout<<"./Lorentz�t�H���_���J���܂���ł���"<<endl;
			exit(1);
		}
	}

	fout2 << "# AVS field file" << endl;
	fout2 << "ndim=1" << endl;
	fout2 << "dim1=" << fluid_number <<endl;
	fout2 << "nspace=3" << endl;
	fout2 << "veclen=3" << endl;
	fout2 << "data=float" << endl;
	fout2 << "field=irregular" << endl;
	fout2 << "label=e-x e-y e-z" << endl << endl;
	fout2 << "variable 1 file=./Lorentz"<<timestep<<" filetype=ascii skip=1 offset=0 stride=6" << endl;
	fout2 << "variable 2 file=./Lorentz"<<timestep<<" filetype=ascii skip=1 offset=1 stride=6" << endl;
	fout2 << "variable 3 file=./Lorentz"<<timestep<<" filetype=ascii skip=1 offset=2 stride=6" << endl;
	fout2 << "coord    1 file=./Lorentz"<<timestep<<" filetype=ascii skip=1 offset=3 stride=6" << endl;
	fout2 << "coord    2 file=./Lorentz"<<timestep<<" filetype=ascii skip=1 offset=4 stride=6" << endl;
	fout2 << "coord    3 file=./Lorentz"<<timestep<<" filetype=ascii skip=1 offset=5 stride=6" << endl;
	fout2.close();

	//�f�[�^�t�@�C��
	stringstream ss;
	ss<<"./Lorentz/Lorentz"<<timestep;
	string filename=ss.str();

	ofstream fout(filename);
	if(fout.fail()){
		cout<<"./Lorentz�t�H���_���J���܂���ł���"<<endl;
		exit(1);
	}

	fout<<"e-x e-y e-z x y z"<<endl;
	for(size_t i=0;i<PART.size();i++)
    {
		if(PART[i].type!=WALL){
//		fout<<F[A_X][i]*times<<" "<<F[A_Y][i]*times<<" "<<F[A_Z][i]*times<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
		fout<<F[A_X][i]<<" "<<F[A_Y][i]<<" "<<F[A_Z][i]<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
		}
	}
	fout.close();
}

//�ӗv�f�����֐�
int make_edge_element(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,int *jnb,int **nei,sides *SIDE,int *branch_num,int **nei2,int KTE)
{
	unsigned timeA=GetTickCount();		//�v�Z�J�n����

	//�ӗv�f����
	int *check=new int [node+1];		//�e�ߓ_������������������ǂ���
	int *temp_check=new int[node+1];	//�ꎞ�I�������z��
	int max=1000;						//�ߓ_�ɗאڂ���ő�ߓ_��
	int *ele_side_num=new int[nelm+1];	//�e�v�f�̑扽�ӂ܂ł����Ƃ܂��Ă��邩
	int side_num=0;						//�S�Ӑ��i�[
		
	///������
	for(int i=1;i<=node;i++) 
	{
		check[i]=0;
		temp_check[i]=0;
		branch_num[i]=0;//������
	}
	for(int i=1;i<=nelm;i++) ele_side_num[i]=0;
	/////////////////////

	//�Ӕԍ��ƕ�-�ߓ_��񐶐�
	for(int i=1;i<=node;i++)
	{
		for(int k=1;k<=jnb[i];k++)
		{
			int jelm=nei[i][k];//�ߓ_i�ɗאڂ���v�f�ԍ�
			for(int j=1;j<=4;j++)
			{
				int p=ELEM[jelm].node[j];//p�͉��H�H
				if(p!=i && temp_check[p]==0)
				{
					branch_num[i]=branch_num[i]+1;//�ߓ_i�̗אڂ���ߓ_�����v���X
					temp_check[p]=1;//��������
					nei2[i][branch_num[i]]=p;

					if(check[p]==0)
					{	
						side_num++;//�ӗv�f���X�V
						SIDE[side_num].node[1]=i;//���̃A���S���Y�����ɂ����ẮA���i<p�ł���.�Ȃ��Ȃ�i��菬���Ȕԍ��̐ߓ_�͂��ł�chek=1������
						SIDE[side_num].node[2]=p;
						///�ߓ_�x�[�X�̋��E������Ӄx�[�X�Ɋg��
						if(NODE[i].boundary_condition==NODE[p].boundary_condition){
							SIDE[side_num].boundary_condition=NODE[i].boundary_condition;//���[���������E�����Ȃ�A���̕ӂ����̋��E�����ɏ]���B���ɗ������m���ł����Ȃ�
						}else{
							SIDE[side_num].boundary_condition=0;//����ȊO�͖��m��
						}
					}
				}
			}
		}
		check[i]=1;
		for(int k=1;k<=branch_num[i];k++) temp_check[nei2[i][k]]=0;//������
//		cout<<"for in i= "<<i<<"/"<<node<<"node, NODE["<<i<<"].attr="<<NODE[i].material<<endl; //�Ȃ���material��FLUID=0�������Ă���E�E�E
	}

	//�v�f-�ӏ�񐶐�
	cout<<"ELEM[].sides�i�["<<endl;	
	for(int i=1;i<=side_num;i++)
	{
		int ia=SIDE[i].node[1];//��i���\������ߓ_�ԍ�(ia<ib)
		int ib=SIDE[i].node[2];
		for(int k=1;k<=jnb[ia];k++)
		{
			int jelm=nei[ia][k];//�ߓ_ia�̗אڂ���v�f�ԍ�
			for(int j=1;j<=4;j++)
			{
				if(ELEM[jelm].node[j]==ib)
				{
					ele_side_num[jelm]=ele_side_num[jelm]+1;
					int a=ele_side_num[jelm];
					ELEM[jelm].sides[a]=i;
				}
			}
		}
	}////////////

	for(int i=1;i<=side_num;i++) SIDE[i].boundary_condition=0;		//�������E�E�E���̈ʒu�̓_���ł́H

	delete [] check;
	delete [] temp_check;
	delete [] ele_side_num;

	std::cout<<"�Ӑ�="<<side_num<<" �ő�Ӑ�="<<KTE<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

	return side_num;
}

//������v�Z�֐�
void calc_transitional_EM_field(mpsconfig &CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,int *jnb,double dt,double TIME,int t,int **nei,int KTE,vector<mpselastic> &PART,int fluid_number,double **F,int KTJ)
{
	sides *SIDE=new sides[KTE];			//�ӃN���X
	int *branch_num=new int[node+1];	//�e�ߓ_���אڂ���ߓ_��(�e�ߓ_���牄�т�ӂ̐�)
	int max=1000;						//�ߓ_�ɗאڂ���ő�ߓ_��
	int **nei2=new int* [node+1];		//�e�ߓ_�̗אڂ���ߓ_�ԍ��i�[(nei�͐ߓ_-�v�f�ł��邪�Anei2�͐ߓ_-�ߓ_)
		for(int i=1;i<=node;i++) nei2[i]=new int [max];
		int side_num=0;						//�S�Ӑ��i�[
	int *depth=new int [KTE+1];//�e�v�f�̐[���i�[

	Post post;
	double *V=new double[node+1];
	for(int i=1;i<=node;i++) V[i]=0;				//�������@V�͉Q�d���̓d�ʊi�[�Ɏg��
	double *current[3];								//�e�v�f�̓d�����x[A/m2]�i�[
	for(int D=0;D<3;D++) current[D]=new double [nelm+1];
	double *Be[3];					//�v�f���d�E or �v�f�����E
    for(int D=0;D<3;D++) Be[D]=new double [nelm+1];
	double *RP=new double [nelm+1];	//�e�v�f�̓�����

	side_num=make_edge_element(CON,NODE,ELEM,node,nelm,jnb,nei,SIDE,branch_num,nei2,KTE);
	//�v�f�[������
	set_depth(CON,NODE,ELEM,node,nelm,depth,KTE);

	//����������
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==FLUID || ELEM[i].material==MAGELAST || ELEM[i].material==MAGELAST2) RP[i]=CON.get_RP();	
		else if(ELEM[i].material==ELASTIC) RP[i]=1;
		else if(ELEM[i].material==AIR) RP[i]=1;
		else if(ELEM[i].material==COIL) RP[i]=1;
		else if(ELEM[i].material==IRON) RP[i]=2500;
		else cout<<"error:�ގ��ɑ΂��A���������s��	�ގ�="<<ELEM[i].material<<endl;
	}

	if(CON.get_J_input_way()==2)					//�d�����x�𑼂̃\�t�g����ǂݍ���
	{
		import_J0_density(CON, node, nelm,NODE,ELEM,current);
		check_J0(CON, node, nelm,NODE,ELEM,current);
	}
	else if(CON.get_J_input_way()==1)
	{
		//������
	}
	//////////////////��ٓd���l�����肷��
				int nm=450;			//��������
				double I0=CON.get_I0();			//�d���̐U��
				double f=CON.get_Hz();				//�𗬂̎��g��[Hz]
			//	double II=I0*cos(2*3.1415*f*TIME)+900;		//�d���g�`
				double II=I0*(-cos(2*3.1415*f*TIME))+450;		//0����
			//	double II=I0*sin(2*3.1415*f*TIME);		//SIN�g�d���g�`
			//	double II=I0+450;//����
/*				/////////////////////////////
				double II;
				if(TIME<(1/100.0)/2.0){
					II=I0*(-cos(2*3.1415*100*TIME))+450;
				}
				else II=I0+450;
				if(II==0)II=0.00000001;
			//	cout<<"Current="<<II<<endl;
				/////////////////////////////*/
		///////////////////////////////////////////
				ofstream Hz("A-t.dat", ios::app);
				Hz<<TIME<<" "<<II<<endl;
				Hz.close();
	int check=1;
	if(CON.get_FEM_calc_type()==3)//������
	{
		if(check==0){//�ړ_�v�f
			double *A[3];									//�ߓ_�ɂ�����x�N�g���|�e���V����
			for(int D=0;D<3;D++) A[D]=new double [node+1];
			double *old_A[3];		
			for(int D=0;D<3;D++) old_A[D]=new double [node+1]; //1step�O�̃x�N�g���|�e���V����
				
			/////old_A�ɒl���i�[
			if(t==1 && CON.get_restart()==OFF)				//�ŏ��̃X�e�b�v�̓[���ŏ�����
			{
				for(int i=1;i<=node;i++) for(int D=0;D<3;D++) old_A[D][i]=0;	//������
			}
			else//����ȊO�̓t�@�C������ǂݍ���
			{
				ifstream old("old_A.dat");
				if(!old) cout<<"cannot open old_A.dat"<<endl;
				old.unsetf(ifstream::dec);
				old.setf(ifstream::skipws);

				for(int i=1;i<=node;i++) for(int D=0;D<3;D++) old>>old_A[D][i];
					
				old.close();
			}
		
			Avector3D_node_eddy2(CON,NODE,ELEM,SIDE,node,nelm,side_num,A,jnb,depth,nei2,branch_num,old_A,II,dt,fluid_number,V,t);
			//�������x����
			Bflux3D_node(CON,NODE,ELEM,node,nelm,A,Be);
			/////////////////////////////////////////////////////
			//�������x���v���b�g
			if(post.get_plot_B()==true){
				plot_magnetic_flux_density(CON, PART, NODE, ELEM, nelm, Be, t);
				C_Fluix_data_avs(node,nelm,NODE,ELEM,KTJ,Be,CON,t);
			
			}

			//�͂̌v�Z
			if(CON.get_m_force()==1) NODE_F3D(CON,NODE,ELEM,node,nelm,Be,jnb,nei,RP,PART,F,fluid_number);
			if(CON.get_m_force()==2) kelvin_force3D(CON,PART,NODE,ELEM, node, nelm,Be, fluid_number,F,RP,jnb,nei);
			if(CON.get_m_force()==4) direct_divT3D(CON,PART,NODE,ELEM, node, nelm,Be, fluid_number,F,RP,jnb,nei);

			//�X���[�W���O�Ɨ͂̃v���b�g
			if(post.get_plot_NF()==true){
				smoothingF3D(CON,PART,fluid_number,F,t);
			}
			/////////////////////////////////////////////////////

			for(int D=0;D<3;D++)
			{
				delete [] A[D];
				delete [] old_A[D];
			}
		}
		else if(check==1){//�ӗv�f

			double *A=new double [side_num+1];//�x�N�g���|�e���V����
		
			//�x�N�g���|�e���V�������v�Z����
			Avector3D2(CON,NODE,ELEM,SIDE,node,nelm,side_num,A,jnb,branch_num,current,RP,II,depth,t);//�ӗv�f
//			non_linear_Avector3D(CON,NODE,ELEM,SIDE,node,nelm,side_num,A,jnb,branch_num,current,RP, II,depth,Be, t);
			//�������x���v�Z����
			Bflux3D_side(CON,NODE,ELEM,SIDE,node,nelm,side_num,A,Be,RP);

			//�������x���v���b�g
			if(post.get_plot_B()==true){
				plot_magnetic_flux_density(CON, PART, NODE, ELEM, nelm, Be, t);
				C_Fluix_data_avs(node,nelm,NODE,ELEM,KTJ,Be,CON,t);
			}
		
			//�͂̌v�Z
			if(CON.get_m_force()==1) NODE_F3D(CON,NODE,ELEM,node,nelm,Be,jnb,nei,RP,PART,F,fluid_number);
			if(CON.get_m_force()==2) kelvin_force3D(CON,PART,NODE,ELEM, node, nelm,Be, fluid_number,F,RP,jnb,nei);
			if(CON.get_m_force()==4) direct_divT3D(CON,PART,NODE,ELEM, node, nelm,Be, fluid_number,F,RP,jnb,nei);

			//�X���[�W���O�Ɨ͂̃v���b�g
			if(post.get_plot_NF()==true){
				smoothingF3D(CON,PART,fluid_number,F,t);
			}

			delete [] A;
		
		}
	}

	delete [] V;
	delete [] RP;
	delete [] depth;
	delete [] SIDE;
	delete [] branch_num;
	for(int i=1;i<=node;i++)delete [] nei2[i];
	delete [] nei2;
	for(int D=0;D<3;D++)
	{
		delete [] current[D];
		delete [] Be[D];
	}
}

//�d�����x�ǂݍ��݊֐�
void import_J0_density(mpsconfig &CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **current)
{
	cout<<"current_density.txt��苭���d�����z��ǂݍ���--";
	int id;

	ifstream fp("current_density.dat");
	if(!fp) cout<<"cannot open current_density.dat"<<endl;
	fp.unsetf(ifstream::dec);
	fp.setf(ifstream::skipws);

	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==COIL)						//current_density.dat�ɂ́A�R�C���v�f�̂ݏo�͂���Ă���d�l�ɂ��Ă�������
		{												//���̏ꍇ�A�R�C���v�f�͐ÓI�v�f�łȂ���΂Ȃ�Ȃ��B���I�Ȃ瑼�̃\�t�g����ǂݍ��߂Ȃ�
			fp>>id;
			for(int D=0;D<3;D++) fp>>current[D][i];		//
		}
		else for(int D=0;D<3;D++) current[D][i]=0;
	}
	fp.close();
	cout<<"ok"<<endl;
}

//�d�����x�o�͊֐�
void check_J0(mpsconfig &CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **current)
{
	int coil_num=0;							//�R�C���v�f��
	for(int i=1;i<=nelm;i++) if(ELEM[i].material==COIL) coil_num++;

	ofstream fout2("J0.fld");
	fout2 << "# AVS field file" << endl;
	fout2 << "ndim=1" << endl;
	//fout2 << "dim1=" << fluid_number <<endl;
	fout2 << "dim1=" << coil_num <<endl;
	fout2 << "nspace=3" << endl;
	fout2 << "veclen=3" << endl;
	fout2 << "data=float" << endl;
	fout2 << "field=irregular" << endl;
	fout2 << "label=e-x e-y e-z" << endl << endl;
	fout2 << "variable 1 file=./J0 filetype=ascii skip=1 offset=0 stride=6" << endl;
	fout2 << "variable 2 file=./J0 filetype=ascii skip=1 offset=1 stride=6" << endl;
	fout2 << "variable 3 file=./J0 filetype=ascii skip=1 offset=2 stride=6" << endl;
	fout2 << "coord    1 file=./J0 filetype=ascii skip=1 offset=3 stride=6" << endl;
	fout2 << "coord    2 file=./J0 filetype=ascii skip=1 offset=4 stride=6" << endl;
	fout2 << "coord    3 file=./J0 filetype=ascii skip=1 offset=5 stride=6" << endl;
	fout2.close();

	ofstream fout("J0");

	fout<<"e-x e-y e-z x y z"<<endl;
	double times=1e-12;
	for(int i=1;i<=nelm;i++)
    {
		if(ELEM[i].material==COIL)
		{
			double r[3]={0,0,0};
			for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
			fout<<current[A_X][i]*times<<" "<<current[A_Y][i]*times<<" "<<current[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
			
		}
	}
	fout.close();
}

//�ߓ_�v�f�������x�v�Z�֐�
void Bflux3D_node(mpsconfig &CON, vector<point3D> &NODE, vector<element3D> &ELEM, int node,int nelm, double **A, double **B)
{
	cout<<"�������x�v�Z----------";

	int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];

	double times=CON.get_B_times();
	double le=CON.get_distancebp();
	int plot_type=CON.get_plot_B_type();//1:�x�N�g���@2:�X�J���[
	//ofstream fp2("test.dat");
    for(int je=1;je<=nelm;je++)
    {   
		for(int D=0;D<3;D++) B[D][je]=0;
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
		double delta6=ELEM[je].volume;//�̐ς�6�{

		delta6=1/delta6;

		double Xs=0;
		double Ys=0;
		double Zs=0;
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
			Xs+=X[j];Ys+=Y[j];Zs+=Z[j];
		}
		Xs/=4;Ys/=4;Zs/=4;
		//�W���쐬
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
			int m=j%4+1;
			int n=m%4+1;
	    
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
			if(i%2!=0)//i����Ȃ�
			{
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////

		for(int i=1;i<=4;i++)
		{
			B[A_X][je]+=d[i]*A[A_Z][N[i]]-e[i]*A[A_Y][N[i]];
			B[A_Y][je]+=e[i]*A[A_X][N[i]]-c[i]*A[A_Z][N[i]];
			B[A_Z][je]+=c[i]*A[A_Y][N[i]]-d[i]*A[A_X][N[i]];
		}

		for(int D=0;D<3;D++) B[D][je]*=delta6;

		//if(Ys>-5*le && Ys<5*le) fp2<<Xs<<" "<<Zs<<" "<<B[A_X][je]*100<<" "<<B[A_Z][je]*100<<endl;

	}
	cout<<"ok"<<endl;
//	fp2.close();

	

	//�������x�o��
	ofstream fp("Bflux.dat");
	
	double Xmin=CON.get_XL()+le; double Xmax=CON.get_XR()-le;//��͗̈�
	double Zmin=CON.get_ZD()+le; double Zmax=CON.get_ZU()-le;
	double Rmax=CON.get_RU()-le;
	
	double dx=5*le;
	int Nx=(int)((Xmax-Xmin)/dx);//�e�����̕�����
	int Nr=(int)(2*Rmax/dx);
	int Nz=(int)((Zmax-Zmin)/dx);
	int search=nelm;//locate�֐��ōŏ��ɒT������v�f�ԍ�
	
	if(CON.get_region_shape()==1) //�~���̈�̂Ƃ��͕ϐ������������ď���
	{
		Nx=Nr;
		Xmin=-Rmax;
	}


	if(plot_type==1)		//�x�N�g���\��
	{
		for(int n=0;n<Nx;n++)
		{
			for(int m=0;m<Nz;m++)
			{
				double xp=dx*n+Xmin;//�o�͂���_�̍��W
				double yp=0;
				double zp=dx*m+Zmin;
				int loc=locate3D(NODE,ELEM,search,xp,yp,zp);
				fp<<xp<<" "<<zp<<" "<<B[A_X][loc]*times<<" "<<B[A_Z][loc]*times<<endl;
				search=loc;
				//if(loc==0) cout<<"EE"<<endl;
			}
		}
	}
	else if(plot_type==2)	//�X�J���[�\��
	{
		for(int n=0;n<Nx;n++)
		{
			for(int m=0;m<Nz;m++)
			{
				double xp=dx*n+Xmin;//�o�͂���_�̍��W
				double yp=0;
				double zp=dx*m+Zmin;
				int loc=locate3D(NODE,ELEM,search,xp,yp,zp);
				double BB=sqrt(B[A_X][loc]*B[A_X][loc]+B[A_Y][loc]*B[A_Y][loc]+B[A_Z][loc]*B[A_Z][loc]);
				fp<<xp<<" "<<zp<<" "<<BB<<endl;
				search=loc;
			}
		}
	}
	fp.close();//*/

	//���̕t�߂����s�b�N�A�b�v����eforce
	int count=0;
	for(int i=1;i<=nelm;i++)
    {
		double r[3]={0,0,0};
		for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
		if(r[A_Y]<=0 &&  r[A_Y]>-0.03 && -0.03<r[A_X] && r[A_X]<0.03 &&  0.12<r[A_Z] && r[A_Z]<0.18) 
		count++;
	}


	ofstream fout2("eforce.fld");
	fout2 << "# AVS field file" << endl;
	fout2 << "ndim=1" << endl;
	//fout2 << "dim1=" << fluid_number <<endl;
	fout2 << "dim1=" << count/*nelm*/ <<endl;
	fout2 << "nspace=3" << endl;
	fout2 << "veclen=3" << endl;
	fout2 << "data=float" << endl;
	fout2 << "field=irregular" << endl;
	fout2 << "label=e-x e-y e-z" << endl << endl;
	fout2 << "variable 1 file=./eforce filetype=ascii skip=1 offset=0 stride=6" << endl;
	fout2 << "variable 2 file=./eforce filetype=ascii skip=1 offset=1 stride=6" << endl;
	fout2 << "variable 3 file=./eforce filetype=ascii skip=1 offset=2 stride=6" << endl;
	fout2 << "coord    1 file=./eforce filetype=ascii skip=1 offset=3 stride=6" << endl;
	fout2 << "coord    2 file=./eforce filetype=ascii skip=1 offset=4 stride=6" << endl;
	fout2 << "coord    3 file=./eforce filetype=ascii skip=1 offset=5 stride=6" << endl;
	fout2.close();

	ofstream fout("eforce");
	fout<<"e-x e-y e-z x y z"<<endl;
	for(int i=1;i<=nelm;i++)
    {
		double r[3]={0,0,0};
		for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
		if(r[A_Y]<=0 &&  r[A_Y]>-0.03 && -0.03<r[A_X] && r[A_X]<0.03 &&  0.12<r[A_Z] && r[A_Z]<0.18)
		{
			fout<<B[A_X][i]*times<<" "<<B[A_Y][i]*times<<" "<<B[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
		}
	}
	fout.close();

}

//�\�ʗ̓X���[�W���O�֐�
void smoothingFs3D(mpsconfig &CON, vector<mpselastic> &PART,int fluid_number,double *Fs)
{
	double le=CON.get_distancebp();
    double *newFs=new double [fluid_number];

	int N=abs(CON.get_FEM_smn());//������
	for(int n=0;n<N;n++)
	{
		for(int i=0;i<fluid_number;i++) 
		{  
			newFs[i]=Fs[i];
			if(PART[i].surface==ON)
			{
				int num=1; //�������g���J�E���g���邩��1
				for(int k=0;k<PART[i].N;k++)
				{       
					int j=PART[i].NEI[k];
					if(PART[j].surface==ON && (PART[j].type==FLUID || PART[j].type==ELASTIC))
					{
						num++;
						newFs[i]+=Fs[j];
					}
				}
				newFs[i]/=num;
			}
		} 
		for(int i=0;i<fluid_number;i++)  Fs[i]=newFs[i];
	}

    delete [] newFs;
}

//�x�N�g���|�e���V�����v�Z�֐�(�ӗv�f�p)
void Avector3D(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,double *A,int *jnb,int *branch_num,double **current,double *RP)
{
	cout<<"�x�N�g���|�e���V�����v�Z�J�n---"<<endl;

	double u0=4*PI*0.0000001;	//�^��̓�����
    double v0=1/u0;				//���C��R��
	double j0x,j0y,j0z;			//�d�����x[A/m^3]
	unsigned timeA=GetTickCount();

	//���΂̒�������������
	double MA=CON.get_magnet_angle();
	double magnet_direction[3]={-sin(MA*2*PI/360),0,cos(MA*2*PI/360)};
	double Mx=CON.get_magnet_B()*magnet_direction[A_X];
	double My=CON.get_magnet_B()*magnet_direction[A_Y];
	double Mz=CON.get_magnet_B()*magnet_direction[A_Z];

	//���͕ǂɌŒ苫�E������ݒ�
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==AIR)//���̕\�ʂɗאڂ����C�v�f
		{
			for(int j=1;j<=4;j++)
			{
				int kelm=ELEM[i].elm[j];
				if(kelm==0)
				{
					//��j�ʂɗאڂ���O�p�`�ƌ����������Ă��钸�_�́A��j�Ԑߓ_�ł���
					int p=ELEM[i].node[j];//���E�O�p�`�ɑ����Ȃ��ߓ_
					for(int k=1;k<=6;k++)
					{
						int iside=ELEM[i].sides[k];
						int ia=SIDE[iside].node[1];
						int ib=SIDE[iside].node[2];
						if(ia!=p && ib!=p) SIDE[iside].boundary_condition=1;//�ߓ_p���܂܂Ȃ��ӂ͋��E��
						else SIDE[iside].boundary_condition=0;
					}
				}
			}
			
		}
	}
	////////////////*/

    int NN=0;//�f�B���N���^���E�Ӑ�
    int *dn=new int [side_num+1]; //�e�ӂ��f�B���N���^���E�̉��Ԗڂ��B�f�B���N���^���E��łȂ��Ȃ�side_num+1���i�[
    double *PHAT=new double [CON.get_max_DN()];//�f�B���N���^���l
    
    //�f�B���N���^���E�������́E�E�E�����łȂ���NN���X�V����Ă��Ȃ�
	//NN���v�Z����
	set_boundary_condition3D_edge(CON,NODE,ELEM,SIDE,node,nelm,side_num,dn,&NN,PHAT,A);

	cout<<"�f�B���N����="<<NN<<endl;
	//*/
    
    int pn=side_num-NN;				//���m��
    int *ppn=new int [pn];			//�s���n�Ԗڂ͕�ppn[n]
    int *npp=new int [side_num+1];	//�e�ӂ��s��̉��Ԗڂɂ����邩�B�f�B���N���^�̏ꍇ��pn+1���i�[
    int num=0;

    for(int i=1;i<=side_num;i++)
    {
        if(SIDE[i].boundary_condition==0)//���m��
		{
			ppn[num]=i;
			npp[i]=num;
			num++;
		}
		else npp[i]=pn+1;
    }
   
	cout<<"���m��: "<<pn<<endl;

    //�s��̍ő啝�v�Z
    int mat_w=0;
	int *nume=new int [side_num+1];				//SIDE[i]���܂ޗv�f��
	int *wid= new int [pn+1];

	for(int i=1;i<=side_num;i++) nume[i]=0;
	for(int i=1;i<=nelm;i++)
	{
		for(int j=1;j<=6;j++)
		{
			int edge=ELEM[i].sides[j];
			nume[edge]=nume[edge]+1;
		}
	}						
	//nume[i]�����Ƃ܂���

	for(int i=1;i<=side_num;i++)
	{	//�l����:����ӎ���̗v�f�Q���l����B�܂��A���̓��̈����Ԃɂ������i�K�ŁA�ӂ̐���6�B�Ȍ�A�v�f�������邲�Ƃ�3�ӂ�������B����Ď����ƂȂ�
		int width=6+3*(nume[i]-1);//�ꍇ�ɂ���Ă͂����菬�����l�ɂȂ邩������Ȃ��B���Ǒ������������m�ۂ���Ԃ�ɂ͖��Ȃ�
		if(width>mat_w) mat_w=width;
		if(npp[i]<pn) wid[npp[i]+1]=width;
	}
	delete [] nume;
	///////////////////////////////////////////

//	cout<<"�z��m�ۂ̑O�܂Ő���, pn: "<<pn<<endl;
	//for(int i=0;i<=pn;i++) cout<<"wid["<<i<<"]="<<wid[i]<<endl;

    //�z��m��
    double **G=new double*[pn+1];//�S�̍s��
	for(int i=1;i<=pn;i++){
		G[i]=new double[wid[i]+1];
	}
//	cout<<"�z��m�ې����@G"<<endl;

    int **ROW=new int *[pn+1]; //�e�s�́A��[���v�f�̗�ԍ��L��
	for(int i=1;i<=pn;i++) ROW[i]=new int [wid[i]+1];
//	cout<<"�z��m�ې����@ROW"<<endl;

    int *NUM=new int [pn+1]; //�e�s�́A��[���v�f��
    
//	cout<<"�z��m�ې����@NUM"<<endl;

    for(int i=1;i<=pn;i++)//������
    {
        NUM[i]=0;
        //for(int j=1;j<=mat_w;j++)
		for(int j=1;j<=wid[i];j++)
		{
			G[i][j]=0;
			ROW[i][j]=0;
		}
    }
    double *B=new double [pn];//���s��
    for(int i=0;i<pn;i++) B[i]=0;//������
    //

   
//�S�̍s����쐬����
    
    int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int table[6+1][2+1];//�v�f���\������ӂƂ��̗v�f���ߓ_�ԍ��֌W
	//cout<<"�s��쐬�J�n ";
    for(int je=1;je<=nelm;je++)
    {   
		//�Ӂ|�ߓ_�e�[�u���쐬
		for(int i=1;i<=6;i++)
		{
			int iside=ELEM[je].sides[i];
			int ia=SIDE[iside].node[1];
			int ib=SIDE[iside].node[2];
			for(int j=1;j<=4;j++)
			{
				if(ELEM[je].node[j]==ia) table[i][1]=j;
				else if(ELEM[je].node[j]==ib) table[i][2]=j;
			}
		}////////////////
	
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
		double Xs=0;//�v�f�̏d�S���W
		double Ys=0;
		double Zs=0;
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
			Xs+=X[j];
			Ys+=Y[j];
			Zs+=Z[j];
		}
		Xs/=4;Ys/=4;Zs/=4;
		////////////////////////////

		//�f���[�j�����̍ۂɋ��߂��̐ς́A�X�[�p�[�{�b�N�X���Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
	
		double delta6=ELEM[je].volume;//�̐ς�6�{
		
		delta6=1/delta6;//�v�Z�ɕK�v�Ȃ̂͋t���Ȃ̂ł����Ŕ��]���Ă���
	
		double delta=ELEM[je].volume/6;//�{���̑̐�
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;//i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
			int m=j%4+1;
			int n=m%4+1;
	    
			b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);//i�������̂Ƃ�
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
			if(i%2!=0)//i����Ȃ�
			{
				b[i]*=-1;
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////	

		double rp=RP[je];
	//	if(ELEM[je].material==FLUID) rp=CON.get_RP();//�䓧����
		//�v�f�}�g���N�X�쐬�J�n
		for(int i=1;i<=6;i++)
		{	
			int iside=ELEM[je].sides[i];//�v�fje�̕Ӕԍ�
			if(SIDE[iside].boundary_condition==0)//�f�B���N���łȂ��B�܂薢�m�Ȃ�
			{   
				int I=npp[iside]+1;///��iside�͍s���I�Ԗ�
				int I1=SIDE[iside].node[1];//iside���\������2�_
				int I2=SIDE[iside].node[2];
				int k1=table[i][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
				int k2=table[i][2];
			    for(int j=1;j<=6;j++)
				{
					int jside=ELEM[je].sides[j];
					
					int u1=table[j][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
					int u2=table[j][2];
						
					if(SIDE[jside].boundary_condition==0)//�f�B���N���łȂ��B�܂薢�m�Ȃ�
					{   
						int J=npp[jside]+1;///��jside�͍s���J�Ԗ�
						int flag=0;
						
						int J1=SIDE[jside].node[1];//jside���\������2�_
						int J2=SIDE[jside].node[2];
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6/rp;
								flag=1;
							}
						}
						if(flag==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    
						    G[I][H]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6/rp;
							ROW[I][H]=J;
						}
					}
					////
					else //jside���f�B���N���^���E�ߓ_�Ȃ�
					{
					    int n=dn[jside];
						B[I-1]-=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6*PHAT[n]/rp;
					}//////////*/
				}

				///B[I-1]���v�Z����i���s��j
				if(ELEM[je].material==COIL)
				{
					j0x=current[A_X][je];
					j0y=current[A_Y][je];
					j0z=current[A_Z][je];
					B[I-1]+=u0*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*j0x*delta6/6;
					B[I-1]+=u0*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*j0y*delta6/6;
					B[I-1]+=u0*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*j0z*delta6/6;
				}//////////
				else if(ELEM[je].material==MAGNET)	B[I-1]+=1.0/18.0/delta*((d[k1]*e[k2]-e[k1]*d[k2])*Mx+(e[k1]*c[k2]-c[k1]*e[k2])*My+(c[k1]*d[k2]-d[k1]*c[k2])*Mz);
			}
		}   	
    }
    ///////////////////////*/
	

    int number=0;//�s��̔�[���v�f��
    for(int i=1;i<=pn;i++) number+=NUM[i];
    
    //�s��̎��ۂ̍ő啝�����߂�
	int maxN=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>maxN) maxN=NUM[i];
	
	//���̂܂܂ł�G[i][j]��[j]�̏��������ɕ���ł��Ȃ��̂ŁAG��ROW����ёւ���
	arrange_matrix(pn,NUM,ROW,G);

	//�Ώ̐��`�F�b�N
	check_matrix_symmetry(pn,NUM,ROW,G);

    double *val = new double [number];
    int *ind = new int [number];//��[���v�f�̗�ԍ��i�[�z��
    int *ptr = new int [pn+1];//�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
    
    /////////////////////val,ind ,ptr�ɒl���i�[
    int index=0;
    for(int n=0;n<pn;n++)
    {
        ptr[n]=index;
		for(int m=1;m<=NUM[n+1];m++)
		{
			val[index]=G[n+1][m];
			ind[index]=ROW[n+1][m]-1;
			index++;
		}
    }	    
    ptr[pn]=number;
    ////////////////////*/
    
    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;

	delete [] wid;
    
	cout<<"�s��쐬 ���F"<<maxN<<"/"<<mat_w<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
    

//�s��v�Z�J�n
	double *XX=new double [pn];//�s��̓����i�[
	if(CON.get_FEMCG()==1) ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON.get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON.get_FEMCG()==3) MRTR(CON,B,pn,XX,val,ind,ptr);
	else if(CON.get_FEMCG()==4) ICMRTR(CON,B,pn,XX,val,ind,ptr);
	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		A[i]=XX[n];
		//cout<<A[i]<<endl;
	}
	delete [] XX;
	///////////////////////////*/
    
    delete [] dn;
    delete [] PHAT;
    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;
}
//�x�N�g���|�e���V�����v�Z�֐�(�ӗv�f�p)
void Avector3D2(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,double *A,int *jnb,int *branch_num,double **current,double *RP,double II,int *depth,int t)
{
	cout<<"�x�N�g���|�e���V�����v�Z�J�n---"<<endl;

	double u0=4*PI*0.0000001;	//�^��̓�����
    double v0=1/u0;				//���C��R��
	double j0x,j0y,j0z;			//�d�����x[A/m^3]
	unsigned timeA=GetTickCount();

	//���΂̒�������������
	double MA=CON.get_magnet_angle();
	double magnet_direction[3]={-sin(MA*2*PI/360),0,cos(MA*2*PI/360)};
	double Mx=CON.get_magnet_B()*magnet_direction[A_X];
	double My=CON.get_magnet_B()*magnet_direction[A_Y];
	double Mz=CON.get_magnet_B()*magnet_direction[A_Z];
	
	//////////////////////////////////////////////////////*/
	//////////////////////��ق��Ȃ���Γd���v�Z���s���K�v���Ȃ��B�����ź�ِߓ_������������
	int coil_node_num=0;	//��ِߓ_��
	for(int i=1;i<=node;i++) if(NODE[i].material==COIL) coil_node_num++;
	////////////////////////*/

	int *save_bound=new int [node+1];//�{�֐��͉��x���Ăяo�����̂ŁA�d���Ɋւ��鋫�E���������S�ɏ����ł��Ȃ��B�����ŕۑ�����
	
	if(coil_node_num>0)///��ِߓ_������Ȃ�d�����x�v�Z
	{
		for(int i=1;i<=node;i++) save_bound[i]=NODE[i].boundary_condition;//���E������ۑ�

		if(II!=0)//�d���l����[���Ȃ�d���v�Z
		{
			calc_current(CON,NODE,ELEM,SIDE,node,nelm,side_num,jnb,branch_num,current,depth,II);
			cout<<"�d���v�Z����"<<endl;
		}
		else	//�[���Ȃ珉����
		{
			for(int i=1;i<=nelm;i++) if(ELEM[i].material==COIL) for(int D=0;D<3;D++) current[D][i]=0;
		}
		///�d�����E�����̏������i�����d���͂��Ƃ܂��Ă��邩�炢��Ȃ�)
		for(int i=1;i<=side_num;i++) if(SIDE[i].boundary_condition>=10) SIDE[i].boundary_condition=0;
		for(int i=1;i<=node;i++) if(NODE[i].boundary_condition>=10) NODE[i].boundary_condition=0;
	}
	carrent_vector(CON, NODE, ELEM, nelm, current,t);
	//���͕ǂɌŒ苫�E������ݒ�
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==AIR)//���̕\�ʂɗאڂ����C�v�f
		{
			for(int j=1;j<=4;j++)
			{
				int kelm=ELEM[i].elm[j];
				if(kelm==0)
				{
					//��j�ʂɗאڂ���O�p�`�ƌ����������Ă��钸�_�́A��j�Ԑߓ_�ł���
					int p=ELEM[i].node[j];//���E�O�p�`�ɑ����Ȃ��ߓ_
					for(int k=1;k<=6;k++)
					{
						int iside=ELEM[i].sides[k];
						int ia=SIDE[iside].node[1];
						int ib=SIDE[iside].node[2];
						if(ia!=p && ib!=p) SIDE[iside].boundary_condition=1;//�ߓ_p���܂܂Ȃ��ӂ͋��E��
						else SIDE[iside].boundary_condition=0;
					}
				}
			}
			
		}
	}
	////////////////*/

    int NN=0;//�f�B���N���^���E�Ӑ�
    int *dn=new int [side_num+1]; //�e�ӂ��f�B���N���^���E�̉��Ԗڂ��B�f�B���N���^���E��łȂ��Ȃ�side_num+1���i�[
    double *PHAT=new double [CON.get_max_DN()];//�f�B���N���^���l
    
    //�f�B���N���^���E�������́E�E�E�����łȂ���NN���X�V����Ă��Ȃ�
	//NN���v�Z����
	set_boundary_condition3D_edge(CON,NODE,ELEM,SIDE,node,nelm,side_num,dn,&NN,PHAT,A);

	cout<<"�f�B���N����="<<NN<<endl;
	//*/
    
    int pn=side_num-NN;				//���m��
    int *ppn=new int [pn];			//�s���n�Ԗڂ͕�ppn[n]
    int *npp=new int [side_num+1];	//�e�ӂ��s��̉��Ԗڂɂ����邩�B�f�B���N���^�̏ꍇ��pn+1���i�[
    int num=0;

    for(int i=1;i<=side_num;i++)
    {
        if(SIDE[i].boundary_condition==0)//���m��
		{
			ppn[num]=i;
			npp[i]=num;
			num++;
		}
		else npp[i]=pn+1;
    }
   
	cout<<"���m��: "<<pn<<endl;

    //�s��̍ő啝�v�Z
    int mat_w=0;
	int *nume=new int [side_num+1];				//SIDE[i]���܂ޗv�f��
	int *wid= new int [pn+1];

	for(int i=1;i<=side_num;i++) nume[i]=0;
	for(int i=1;i<=nelm;i++)
	{
		for(int j=1;j<=6;j++)
		{
			int edge=ELEM[i].sides[j];
			nume[edge]=nume[edge]+1;
		}
	}						
	//nume[i]�����Ƃ܂���

	for(int i=1;i<=side_num;i++)
	{	//�l����:����ӎ���̗v�f�Q���l����B�܂��A���̓��̈����Ԃɂ������i�K�ŁA�ӂ̐���4�B�Ȍ�A�v�f�������邲�Ƃ�3�ӂ�������B����Ď����ƂȂ�
		int width=6+3*(nume[i]-1);//�ꍇ�ɂ���Ă͂����菬�����l�ɂȂ邩������Ȃ��B���Ǒ������������m�ۂ���Ԃ�ɂ͖��Ȃ�
		if(width>mat_w) mat_w=width;
		if(npp[i]<pn) wid[npp[i]+1]=width;
	}
	delete [] nume;
	///////////////////////////////////////////

//	cout<<"�z��m�ۂ̑O�܂Ő���, pn: "<<pn<<endl;
	//for(int i=0;i<=pn;i++) cout<<"wid["<<i<<"]="<<wid[i]<<endl;

    //�z��m��
    double **G=new double*[pn+1];//�S�̍s��
	for(int i=1;i<=pn;i++){
		G[i]=new double[wid[i]+1];
	}
//	cout<<"�z��m�ې����@G"<<endl;

    int **ROW=new int *[pn+1]; //�e�s�́A��[���v�f�̗�ԍ��L��
	for(int i=1;i<=pn;i++) ROW[i]=new int [wid[i]+1];
//	cout<<"�z��m�ې����@ROW"<<endl;

    int *NUM=new int [pn+1]; //�e�s�́A��[���v�f��
    
//	cout<<"�z��m�ې����@NUM"<<endl;

    for(int i=1;i<=pn;i++)//������
    {
        NUM[i]=0;
        //for(int j=1;j<=mat_w;j++)
		for(int j=1;j<=wid[i];j++)
		{
			G[i][j]=0;
			ROW[i][j]=0;
		}
    }
    double *B=new double [pn];//���s��
    for(int i=0;i<pn;i++) B[i]=0;//������
    //

   
//�S�̍s����쐬����
    
    int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int table[6+1][2+1];//�v�f���\������ӂƂ��̗v�f���ߓ_�ԍ��֌W
	//cout<<"�s��쐬�J�n ";
    for(int je=1;je<=nelm;je++)
    {   
		//�Ӂ|�ߓ_�e�[�u���쐬
		for(int i=1;i<=6;i++)
		{
			int iside=ELEM[je].sides[i];
			int ia=SIDE[iside].node[1];
			int ib=SIDE[iside].node[2];
			for(int j=1;j<=4;j++)
			{
				if(ELEM[je].node[j]==ia) table[i][1]=j;
				else if(ELEM[je].node[j]==ib) table[i][2]=j;
			}
		}////////////////
	
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
		double Xs=0;//�v�f�̏d�S���W
		double Ys=0;
		double Zs=0;
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
			Xs+=X[j];
			Ys+=Y[j];
			Zs+=Z[j];
		}
		Xs/=4;Ys/=4;Zs/=4;
		////////////////////////////

		//�f���[�j�����̍ۂɋ��߂��̐ς́A�X�[�p�[�{�b�N�X���Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
	
		double delta6=ELEM[je].volume;//�̐ς�6�{
		
		delta6=1/delta6;//�v�Z�ɕK�v�Ȃ̂͋t���Ȃ̂ł����Ŕ��]���Ă���
	
		double delta=ELEM[je].volume/6;//�{���̑̐�
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;//i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
			int m=j%4+1;
			int n=m%4+1;
	    
			b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);//i�������̂Ƃ�
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
			if(i%2!=0)//i����Ȃ�
			{
				b[i]*=-1;
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////	

		double rp=RP[je];
	//	if(ELEM[je].material==FLUID) rp=CON.get_RP();//�䓧����
		//�v�f�}�g���N�X�쐬�J�n
		for(int i=1;i<=6;i++)
		{	
			int iside=ELEM[je].sides[i];//�v�fje�̕Ӕԍ�
			if(SIDE[iside].boundary_condition==0)//�f�B���N���łȂ��B�܂薢�m�Ȃ�
			{   
				int I=npp[iside]+1;///��iside�͍s���I�Ԗ�
				int I1=SIDE[iside].node[1];//iside���\������2�_
				int I2=SIDE[iside].node[2];
				int k1=table[i][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
				int k2=table[i][2];
			    for(int j=1;j<=6;j++)
				{
					int jside=ELEM[je].sides[j];
					
					int u1=table[j][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
					int u2=table[j][2];
						
					if(SIDE[jside].boundary_condition==0)//�f�B���N���łȂ��B�܂薢�m�Ȃ�
					{   
						int J=npp[jside]+1;///��jside�͍s���J�Ԗ�
						int flag=0;
						
						int J1=SIDE[jside].node[1];//jside���\������2�_
						int J2=SIDE[jside].node[2];
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6/rp;
								flag=1;
							}
						}
						if(flag==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    
						    G[I][H]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6/rp;
							ROW[I][H]=J;
						}
					}
					////
					else //jside���f�B���N���^���E�ߓ_�Ȃ�
					{
					    int n=dn[jside];
						B[I-1]-=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6*PHAT[n]/rp;
					}//////////*/
				}

				///B[I-1]���v�Z����i���s��j
				if(ELEM[je].material==COIL)
				{
					j0x=current[A_X][je];
					j0y=current[A_Y][je];
					j0z=current[A_Z][je];
					B[I-1]+=u0*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*j0x*delta6/6;
					B[I-1]+=u0*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*j0y*delta6/6;
					B[I-1]+=u0*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*j0z*delta6/6;
				}//////////
				else if(ELEM[je].material==MAGNET)
				{
					B[I-1]+=1.0/18.0/delta*((d[k1]*e[k2]-e[k1]*d[k2])*Mx+(e[k1]*c[k2]-c[k1]*e[k2])*My+(c[k1]*d[k2]-d[k1]*c[k2])*Mz);
				}
			}
		}   	
    }
    ///////////////////////*/
	

    int number=0;//�s��̔�[���v�f��
    for(int i=1;i<=pn;i++) number+=NUM[i];
    
    //�s��̎��ۂ̍ő啝�����߂�
	int maxN=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>maxN) maxN=NUM[i];
	
	//���̂܂܂ł�G[i][j]��[j]�̏��������ɕ���ł��Ȃ��̂ŁAG��ROW����ёւ���
	arrange_matrix(pn,NUM,ROW,G);

	//�Ώ̐��`�F�b�N
	check_matrix_symmetry(pn,NUM,ROW,G);

    double *val = new double [number];
    int *ind = new int [number];//��[���v�f�̗�ԍ��i�[�z��
    int *ptr = new int [pn+1];//�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
    
    /////////////////////val,ind ,ptr�ɒl���i�[
    int index=0;
    for(int n=0;n<pn;n++)
    {
        ptr[n]=index;
		for(int m=1;m<=NUM[n+1];m++)
		{
			val[index]=G[n+1][m];
			ind[index]=ROW[n+1][m]-1;
			index++;
		}
    }	    
    ptr[pn]=number;
    ////////////////////*/
    
    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;

	delete [] wid;
    
	cout<<"�s��쐬 ���F"<<maxN<<"/"<<mat_w<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
    

//�s��v�Z�J�n
	double *XX=new double [pn];//�s��̓����i�[
	if(CON.get_FEMCG()==1) ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON.get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON.get_FEMCG()==3) MRTR(CON,B,pn,XX,val,ind,ptr);
	else if(CON.get_FEMCG()==4) ICMRTR(CON,B,pn,XX,val,ind,ptr);
	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		A[i]=XX[n];
		//cout<<A[i]<<endl;
	}
	delete [] XX;
	///////////////////////////*/
    
    delete [] dn;
    delete [] PHAT;
    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;
}
void Avector3D_OK(mpsconfig &CON,vector <point3D> &NODE,vector <element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,double *A,int *jnb,int *branch_num,int *depth)
{
	cout<<"�x�N�g���|�e���V�����v�Z�J�n---"<<endl;

	double u0=4*PI*0.0000001;	//�^��̓�����
    double v0=1/u0;				//���C��R��
	double j0x,j0y,j0z;			//�d�����x[A/m^3]
	unsigned timeA=GetTickCount();

	//���΂̒�������������
	double MA=CON.get_magnet_angle();
	double magnet_direction[3]={-sin(MA*2*PI/360),0,cos(MA*2*PI/360)};
	double Mx=CON.get_magnet_B()*magnet_direction[A_X];
	double My=CON.get_magnet_B()*magnet_direction[A_Y];
	double Mz=CON.get_magnet_B()*magnet_direction[A_Z];
	

	//////////////////////////////////////////////////////*/

	///���͕ǂɌŒ苫�E������ݒ�
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==AIR)//���̕\�ʂɗאڂ����C�v�f
		{
			for(int j=1;j<=4;j++)
			{
				int kelm=ELEM[i].elm[j];
				if(kelm==0)
				{
					///��j�ʂɗאڂ���O�p�`�ƌ����������Ă��钸�_�́A��j�Ԑߓ_�ł���
					int p=ELEM[i].node[j];//���E�O�p�`�ɑ����Ȃ��ߓ_
					for(int k=1;k<=6;k++)
					{
						int iside=ELEM[i].sides[k];
						int ia=SIDE[iside].node[1];
						int ib=SIDE[iside].node[2];
						if(ia!=p && ib!=p) SIDE[iside].boundary_condition=1;//�ߓ_p���܂܂Ȃ��ӂ͋��E��
						else SIDE[iside].boundary_condition=0;
					}
				}
			}
			
		}
	}
	////////////////*/

    int NN=0;//�f�B���N���^���E�Ӑ�
    int *dn=new int [side_num+1]; //�e�ӂ��f�B���N���^���E�̉��Ԗڂ��B�f�B���N���^���E��łȂ��Ȃ�side_num+1���i�[
    double *PHAT=new double [CON.get_max_DN()];//�f�B���N���^���l
    
    ///�f�B���N���^���E��������
	set_boundary_condition3D_edge(CON,NODE,ELEM,SIDE,node,nelm,side_num,dn,&NN,PHAT,A);

	cout<<"�f�B���N����="<<NN;
	/////////////*/
    
	    
    int pn=side_num-NN;				///���m��
    int *ppn=new int [pn];			//�s���n�Ԗڂ͕�ppn[n]
    int *npp=new int [side_num+1];	///�e�ӂ��s��̉��Ԗڂɂ����邩�B�f�B���N���^�̏ꍇ��pn+1���i�[
    int num=0; 
    for(int i=1;i<=side_num;i++)
    {
        if(SIDE[i].boundary_condition==0)//���m��
		{
			ppn[num]=i;
			npp[i]=num;
			num++;
		}
		else npp[i]=pn+1;
    }
   
    ////�s��̍ő啝�v�Z
    int mat_w=0;
	int *nume=new int [side_num+1];				//SIDE[i]���܂ޗv�f��
	int *wid= new int [pn+1];

	for(int i=1;i<=side_num;i++) nume[i]=0;
	for(int i=1;i<=nelm;i++)
	{
		for(int j=1;j<=6;j++)
		{
			int edge=ELEM[i].sides[j];
			nume[edge]=nume[edge]+1;
		}
	}											///nume[i]�����Ƃ܂���
	for(int i=1;i<=side_num;i++)
	{	//�l����:����ӎ���̗v�f�Q���l����B�܂��A���̓��̈����Ԃɂ������i�K�ŁA�ӂ̐���4�B�Ȍ�A�v�f�������邲�Ƃ�3�ӂ�������B����Ď����ƂȂ�
		int width=6+3*(nume[i]-1);//�ꍇ�ɂ���Ă͂����菬�����l�ɂȂ邩������Ȃ��B���Ǒ������������m�ۂ���Ԃ�ɂ͖��Ȃ�
		if(width>mat_w) mat_w=width;
		if(npp[i]<pn) wid[npp[i]+1]=width;
	}
	delete [] nume;
	///////////////////////////////////////////

	
    ////�z��m��
    double **G=new double *[pn+1];///�S�̍s��
	for(int i=1;i<=pn;i++) G[i]=new double [wid[i]+1];
    int **ROW=new int *[pn+1]; ///�e�s�́A��[���v�f�̗�ԍ��L��
	for(int i=1;i<=pn;i++) ROW[i]=new int [wid[i]+1];
    int *NUM=new int [pn+1]; ///�e�s�́A��[���v�f��
    
    for(int i=1;i<=pn;i++)//������
    {
        NUM[i]=0;
        //for(int j=1;j<=mat_w;j++)
		for(int j=1;j<=wid[i];j++)
		{
			G[i][j]=0;
			ROW[i][j]=0;
		}
    }
    double *B=new double [pn];//���s��
    for(int i=0;i<pn;i++) B[i]=0;//������
    ////
    
    /////////�S�̍s����쐬����
    
    int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int table[6+1][2+1];//�v�f���\������ӂƂ��̗v�f���ߓ_�ԍ��֌W
	//cout<<"�s��쐬�J�n ";
    for(int je=1;je<=nelm;je++)
    {   
		//if(ELEM[je].material!=7) cout<<ELEM[je].material<<endl;
		//�Ӂ|�ߓ_�e�[�u���쐬
		for(int i=1;i<=6;i++)
		{
			int iside=ELEM[je].sides[i];
			int ia=SIDE[iside].node[1];
			int ib=SIDE[iside].node[2];
			for(int j=1;j<=4;j++)
			{
				if(ELEM[je].node[j]==ia) table[i][1]=j;
				else if(ELEM[je].node[j]==ib) table[i][2]=j;
			}
		}////////////////
	
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
		double Xs=0;//�v�f�̏d�S���W
		double Ys=0;
		double Zs=0;
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
			Xs+=X[j];
			Ys+=Y[j];
			Zs+=Z[j];
		}
		Xs/=4;Ys/=4;Zs/=4;
		////////////////////////////

		///�f���[�j�����̍ۂɋ��߂��̐ς́A�X�[�p�[�{�b�N�X���Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
	
		double delta6=ELEM[je].volume;//�̐ς�6�{
		
		delta6=1/delta6;//�v�Z�ɕK�v�Ȃ̂͋t���Ȃ̂ł����Ŕ��]���Ă���
	
		double delta=ELEM[je].volume/6;//�{���̑̐�
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
			int m=j%4+1;
			int n=m%4+1;
	    
			b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);//i�������̂Ƃ�
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
			if(i%2!=0)//i����Ȃ�
			{
				b[i]*=-1;
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
	
		
		double rp=CON.get_RP();
	//	if(ELEM[je].material==FLUID) rp=CON.get_RP();//�䓧����
		////�v�f�}�g���N�X�쐬�J�n
		for(int i=1;i<=6;i++)
		{	
			int iside=ELEM[je].sides[i];//�v�fje�̕Ӕԍ�
			if(SIDE[iside].boundary_condition==0)///�f�B���N���łȂ��B�܂薢�m�Ȃ�
			{   
				int I=npp[iside]+1;///��iside�͍s���I�Ԗ�
				int I1=SIDE[iside].node[1];//iside���\������2�_
				int I2=SIDE[iside].node[2];
				int k1=table[i][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
				int k2=table[i][2];
			    for(int j=1;j<=6;j++)
				{
					int jside=ELEM[je].sides[j];
					
					int u1=table[j][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
					int u2=table[j][2];
						
					if(SIDE[jside].boundary_condition==0)///�f�B���N���łȂ��B�܂薢�m�Ȃ�
					{   
						int J=npp[jside]+1;///��jside�͍s���J�Ԗ�
						int flag=0;
						
						int J1=SIDE[jside].node[1];//jside���\������2�_
						int J2=SIDE[jside].node[2];
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6/rp;
								flag=1;
							}
						}
						if(flag==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    
						    G[I][H]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6/rp;
							ROW[I][H]=J;
						}
					}
					////
					else //jside���f�B���N���^���E�ߓ_�Ȃ�
					{
					    int n=dn[jside];
						B[I-1]-=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6*PHAT[n]/rp;
					}//////////*/
				}
				double* current[3]={0,0,0};
				///B[I-1]���v�Z����
				if(ELEM[je].material==COIL)
				{
					j0x=current[A_X][je];
					j0y=current[A_Y][je];
					j0z=current[A_Z][je];
					B[I-1]+=u0*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*j0x*delta6/6;
					B[I-1]+=u0*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*j0y*delta6/6;
					B[I-1]+=u0*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*j0z*delta6/6;
				}//////////
				else if(ELEM[je].material==MAGNET)
				{
					//cout<<Mz<<endl;
					B[I-1]+=1.0/18.0/delta*((d[k1]*e[k2]-e[k1]*d[k2])*Mx+(e[k1]*c[k2]-c[k1]*e[k2])*My+(c[k1]*d[k2]-d[k1]*c[k2])*Mz);
				}
			}
		}   	
    }
    ///////////////////////*/
	

    int number=0;//�s��̔�[���v�f��
    for(int i=1;i<=pn;i++) number+=NUM[i];
    
    ///�s��̎��ۂ̍ő啝�����߂�
	int maxN=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>maxN) maxN=NUM[i];
	

	///���̂܂܂ł�G[i][j]��[j]�̏��������ɕ���ł��Ȃ��̂ŁAG��ROW����ёւ���
	arrange_matrix(pn,NUM,ROW,G);

	//�Ώ̐��`�F�b�N
	check_matrix_symmetry(pn,NUM,ROW,G);

    double *val = new double [number];
    int *ind = new int [number];//��[���v�f�̗�ԍ��i�[�z��
    int *ptr = new int [pn+1];//�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
    
    /////////////////////val,ind ,ptr�ɒl���i�[
    int index=0;
    for(int n=0;n<pn;n++)
    {
        ptr[n]=index;
		for(int m=1;m<=NUM[n+1];m++)
		{
			val[index]=G[n+1][m];
			ind[index]=ROW[n+1][m]-1;
			index++;
		}
    }	    
    ptr[pn]=number;
    ////////////////////*/
    
    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;

	delete [] wid;
    
	cout<<"�s��쐬 ���F"<<maxN<<"/"<<mat_w<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
    

	///////////////////////�s��v�Z�J�n
	double *XX=new double [pn];//�s��̓����i�[
	if(CON.get_FEMCG()==1) ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON.get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON.get_FEMCG()==3) MRTR(CON,B,pn,XX,val,ind,ptr);
	else if(CON.get_FEMCG()==4) ICMRTR(CON,B,pn,XX,val,ind,ptr);
	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		A[i]=XX[n];
		//cout<<A[i]<<endl;
	}
	delete [] XX;
	///////////////////////////*/
    
    
    delete [] dn;
    delete [] PHAT;
    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;
}

//���E�����K�p�֐�(�R�c�ӗv�f�p)
void set_boundary_condition3D_edge(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,int *dn,int *NN,double *PHAT,double *A)
{
	int N=0;	//�����グ�ϐ�

	//�f�B���N���^���E��������
	if(CON.get_uniform_B_sw()==OFF)
	{
		for(int i=1;i<=side_num;i++)
		{
		    if(SIDE[i].boundary_condition==1)
			{
		        dn[i]=N;
		        PHAT[N]=0;
		        A[i]=0;
		        N++;
			}
			else if(SIDE[i].boundary_condition==2)
			{   
				/*int ia=SIDE[i].node[1];
				int ib=SIDE[i].node[2];
				A[i]=(NODE[ib].r[A_Z]-NODE[ia].r[A_Z])*0;
		        dn[i]=NN;
		        PHAT[NN]=(NODE[ib].r[A_Z]-NODE[ia].r[A_Z])*0;
		        //A[i]=0;
		        NN++;*/
				dn[i]=N;
		        PHAT[N]=0.0;
		        A[i]=0.0;
		        N++;
			}
			else
			{
				dn[i]=side_num+1;
			}
		}
	}
	if(CON.get_uniform_B_sw()==ON)	//��͗̈�S�̂Ɉ�l�����^����ꍇ
	{
		double B=CON.get_uniform_B();//��l����̑傫��[�s�n
		double R[3];					//�ӂ̒����i�[(X,Y,Z����)
		double r[3];					//�ӂ̒��_���W�i�[
		double err=1e-14;
		
		for(int i=1;i<=side_num;i++)
		{
		    if(SIDE[i].boundary_condition==1)
			{   
				int ia=SIDE[i].node[1];
				int ib=SIDE[i].node[2];
				for(int D=0;D<3;D++)
				{
					R[D]=NODE[ib].r[D]-NODE[ia].r[D];
					r[D]=(NODE[ib].r[D]+NODE[ia].r[D])*0.5;//���_
				}

			//	if(r[A_Z]<CON.get_ZU()-err && r[A_Z]>CON.get_ZD()+err)
				{

					double L=sqrt(R[A_X]*R[A_X]+R[A_Y]*R[A_Y]+R[A_Z]*R[A_Z]);//�ӂ̒���
				
					A[i]=0.5*(-B*r[A_Y]*R[A_X]+B*r[A_X]*R[A_Y]);	//�Ȃ������Ȃ�̂��̓X�g�[�N�X�̒藝�𗘗p����΂킩��B
					dn[i]=N;
					PHAT[N]=0.5*(-B*r[A_Y]*R[A_X]+B*r[A_X]*R[A_Y]);
					N++;
				}
				//else 
				//{
				////	SIDE[i].boundary_condition=0;//���R���E����
				//	dn[i]=side_num+1;
				//}
			}
			else dn[i]=side_num+1;
		}
	}

	*NN=N;
}

//�ӗv�f�������x�v�Z�֐�
void Bflux3D_side(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,double *A,double **B,double *RP)
{
	cout<<"�������x�v�Z�J�n----";
	unsigned timeA=GetTickCount();
	int plot_type=CON.get_plot_B_type();	//�t�@�C���o�͌`���@1=�x�N�g�� 2=�X�J���[
	double times=CON.get_B_times();		//�t�@�C���o�͎��̔{��
	double u0=4*PI*1e-7;					//�^��̓�����

	double *Xg=new double [nelm+1];			//�v�f�̏d�S���W
	double *Yg=new double [nelm+1];
	double *Zg=new double [nelm+1];


	//#pragma omp parallel for
	for(int je=1;je<=nelm;je++)
    {   
		//cout<<je<<" "<<omp_get_thread_num()<<endl;//�ei�̌v�Z��S�����Ă���X���b�h�ԍ��o��
		int N[4+1];								//�v�f�̊e�ߓ_�ԍ��i�[
		double X[4+1];
		double Y[4+1];
		double Z[4+1];
		double c[4+1];
		double d[4+1];
		double e[4+1];
		int table[6+1][2+1];					//�v�f���\������ӂƂ��̗v�f���ߓ_�ԍ��֌W

		//�Ӂ|�ߓ_�e�[�u���쐬
		for(int i=1;i<=6;i++)
		{
			int iside=ELEM[je].sides[i];
			int ia=SIDE[iside].node[1];
			int ib=SIDE[iside].node[2];
			for(int j=1;j<=4;j++)
			{
				if(ELEM[je].node[j]==ia) table[i][1]=j;
				else if(ELEM[je].node[j]==ib) table[i][2]=j;
			}
		}////////////////
	
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
		double Xs=0;//�v�f�̏d�S���W
		double Ys=0;
		double Zs=0;
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
			Xs+=X[j]*0.25;
			Ys+=Y[j]*0.25;
			Zs+=Z[j]*0.25;
		}
		Xg[je]=Xs; Yg[je]=Ys; Zg[je]=Zs;	//�d�S���
		////////////////////////////

		double delta6=ELEM[je].volume;//�̐ς�6�{(�������̐ς̒l�͂��łɃx�N�g���|�e���V���������߂�ۂɌv�Z���Ă���)
	
		delta6=1/delta6;//�v�Z�ɕK�v�Ȃ̂͋t���Ȃ̂ł����Ŕ��]���Ă���
	
		double delta=ELEM[je].volume/6;//�{���̑̐�
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
			int m=j%4+1;
			int n=m%4+1;
	    
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
			if(i%2!=0)//i����Ȃ�
			{
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////

		for(int D=0;D<3;D++) B[D][je]=0;//������

		for(int i=1;i<=6;i++)
		{
			int s=ELEM[je].sides[i];
			int k1=table[i][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
			int k2=table[i][2];

			B[A_X][je]+=(d[k1]*e[k2]-e[k1]*d[k2])*A[s]; //c=x, d=y, e=z�ɑΉ�, B=rotA���v�Z���Ă���
			B[A_Y][je]+=(e[k1]*c[k2]-c[k1]*e[k2])*A[s];
			B[A_Z][je]+=(c[k1]*d[k2]-d[k1]*c[k2])*A[s];
			//cout<<A[s]<<endl;
		}

		double BB=0;
		for(int D=0;D<3;D++)
		{
			B[D][je]*=delta6*delta6*2; //��(3.9)
		//	BB+=B[D][je]*B[D][je];
			//cout<<B[D][je]<<endl;
		}
		//cout<<sqrt(BB)<<"[T]"<<endl;
	}

	//�������x�o��
	string filename[4];
	filename[0]="vectorBflux";
	filename[1]="vectorField";
	filename[2]="scalarBflux";
	filename[3]="scalarField";

	for(int suffix=0;suffix<4;suffix++)
	{
		stringstream ss;

		if(suffix%2==0){
			ss<<"./Bflux/"<<filename[suffix]<<CON.get_current_step()<<".dat";
		}else{
			ss<<"./Field/"<<filename[suffix]<<CON.get_current_step()<<".dat";
		}
		filename[suffix]=ss.str();
	}

	ofstream vectorflux(filename[0]);
	ofstream vectorfield(filename[1]);
	ofstream scalarflux(filename[2]);
	ofstream scalarfield(filename[3]);

	double le=CON.get_distancebp();
	double Xmin=CON.get_XL()+le; double Xmax=CON.get_XR()-le; //��͗̈�(���E���傤�ǂ��ƁA���̓_���܂ޗv�f���ł��؂�덷�Ō�����Ȃ��Ƃ�������̂�le�����ی�������)
	double Zmin=CON.get_ZD()+le; double Zmax=CON.get_ZU()-le;

	if(CON.get_region_shape()==1)//�~���`�̉�͗̈�Ȃ�
	{
		Xmin=-1*CON.get_RU()+le; 
		Xmax=CON.get_RU()-le;
	}
	
	double dx=3*le;
	int Nx=(int)((Xmax-Xmin)/dx);								//�e�����̕�����
	int Nz=(int)((Zmax-Zmin)/dx);
	int search=nelm;												//locate�֐��ōŏ��ɒT������v�f�ԍ�
	
//	if(plot_type==1)	//�x�N�g���\��
	{
		for(int n=0;n<Nx;n++)
		{
			for(int m=0;m<Nz;m++)
			{
				double xp=dx*n+Xmin;							//�o�͂���_�̍��W
				double yp=0;
				double zp=dx*m+Zmin;
				int loc=locate3D(NODE,ELEM,search,xp,yp,zp);
				vectorflux<<xp<<"\t"<<zp<<"\t"<<B[A_X][loc]*times<<"\t"<<B[A_Z][loc]*times<<endl;
				vectorfield<<xp<<"\t"<<zp<<"\t"<<B[A_X][loc]*times/(u0*RP[loc])<<"\t"<<B[A_X][loc]*times/(u0*RP[loc])<<endl;
				search=loc;
			}
		}
	}
//	else if(plot_type==2)										//�X�J���[�\��
	{
		for(int n=0;n<Nx;n++)
		{
			for(int m=0;m<Nz;m++)
			{
				double xp=dx*n+Xmin;							//�o�͂���_�̍��W
				double yp=0;
				double zp=dx*m+Zmin;
				int loc=locate3D(NODE,ELEM,search,xp,yp,zp);
				double BB=sqrt(B[A_X][loc]*B[A_X][loc]+B[A_Y][loc]*B[A_Y][loc]+B[A_Z][loc]*B[A_Z][loc]);
				scalarflux<<xp<<"\t"<<zp<<"\t"<<BB<<endl;
				scalarfield<<xp<<"\t"<<zp<<"\t"<<BB/(u0*RP[loc])<<endl;
				search=loc;
			}
		}
	}

	vectorflux.close();
	scalarflux.close();
	vectorfield.close();
	scalarfield.close();

	delete [] Xg;
	delete [] Yg;
	delete [] Zg;
	cout<<"ok time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
		
}

//�}�N�X�E�F���e���\���̐ϗ͌v�Z�֐�
void direct_divT3D(mpsconfig &CON, vector<mpselastic> &PART,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **Be,int fluid_number,double **F,double *RP,int *jnb,int **nei)
{
    cout<<"div(T)�̒��ڕ]���ɂ��d���͌v�Z--";
	unsigned int timeA=GetTickCount();
    double le=CON.get_distancebp();
	double R=CON.get_re()*le;
    double V=CON.get_particle_volume();//���q�̑̐�
    double u0=1.257e-6;				//�^��̓�����
	int N=4;
	int order=1;
	double Fz=0;
    for(int i=0;i<fluid_number;i++) for(int D=0;D<3;D++) F[D][i]=0;//������
	double *matrix=new double [N*N];	//N�~N�̌W���s��
	double *matrix_val=new double [N*N];	//matrix�̕ۑ��p
	double *Bxx=new double [N];	//N�̉��s��
	double *Bxy=new double [N];	//N�̉��s��
	double *Bxz=new double [N];	//N�̉��s��
	double *Byx=new double [N];	//N�̉��s��
	double *Byy=new double [N];	//N�̉��s��
	double *Byz=new double [N];	//N�̉��s��
	double *Bzx=new double [N];	//N�̉��s��
	double *Bzy=new double [N];	//N�̉��s��
	double *Bzz=new double [N];	//N�̉��s��

	

	if(order==1)//3����1����
	{
		///�W���s���
		///   ����x2    ����x��y  ����x��z ����x  a = ����x��f  
		///  ����x��y   ����y2    ����y��z ����y  b = ����y��f 
		///  ����x��z   ����y��z  ����z2   ����z  c = ����z��f 
		///  ����x      ����y     ����z    ��1    d = ����f

		double Ri[3];//�ߓ_i�̍��W�i�[
		double Rg[3];	//�v�f�̏d�S�i�[
		//for(int i=1;i<=NN;i++)//���̐ߓ_
		for(int i=1;i<=node;i++)//���̐ߓ_
		{
			if(NODE[i].material==FLUID || NODE[i].material==ELASTIC || NODE[i].material==MAGELAST)
			{
				for(int n=0;n<N*N;n++) matrix[n]=0;		//������
				for(int n=0;n<N;n++) 
				{
					Bxx[n]=0;  Bxy[n]=0;  Bxz[n]=0;
					Byx[n]=0;  Byy[n]=0;  Byz[n]=0;
					Bzx[n]=0;  Bzy[n]=0;  Bzz[n]=0;
				}
				for(int D=0;D<3;D++) Ri[D]=NODE[i].r[D];
				//if(jnb[i]<3) cout<<i<<" "<<jnb[i]<<endl;
				if(jnb[i]>3)
				{
					for(int j=1;j<=jnb[i];j++)//�ߓ_i��jnb[n]�̗v�f�ɗאڂ��Ă���
					{
						int jelm=nei[i][j];
						double u=RP[jelm]*u0;
						///�}�N�X�E�F���̉��̓e���\��
						double Txx=(Be[A_X][jelm]*Be[A_X][jelm]-Be[A_Y][jelm]*Be[A_Y][jelm]-Be[A_Z][jelm]*Be[A_Z][jelm])/(2*u);
						double Txy=2*Be[A_X][jelm]*Be[A_Y][jelm]/(2*u);
						double Txz=2*Be[A_X][jelm]*Be[A_Z][jelm]/(2*u);
						double Tyx=Txy;
						double Tyy=(-Be[A_X][jelm]*Be[A_X][jelm]+Be[A_Y][jelm]*Be[A_Y][jelm]-Be[A_Z][jelm]*Be[A_Z][jelm])/(2*u);
						double Tyz=2*Be[A_Y][jelm]*Be[A_Z][jelm]/(2*u);
						double Tzx=Txz;
						double Tzy=Tyz;
						double Tzz=(-Be[A_X][jelm]*Be[A_X][jelm]-Be[A_Y][jelm]*Be[A_Y][jelm]+Be[A_Z][jelm]*Be[A_Z][jelm])/(2*u);
						//////////
						for(int D=0;D<3;D++) Rg[D]=0;
						for(int k=1;k<=4;k++)
						{
							int node=ELEM[jelm].node[k];//jelm���\������ߓ_�ԍ�
							for(int D=0;D<3;D++) Rg[D]+=NODE[node].r[D]*0.25;
						}//Rg[D]�����Ƃ܂����B
	
						double X=Rg[A_X]-Ri[A_X];
						double Y=Rg[A_Y]-Ri[A_Y];
						double Z=Rg[A_Z]-Ri[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
					
						double w=1;
						if(dis>le) w=le*le/(dis*dis);
		
						matrix[0]+=X*X*w;		
						matrix[1]+=X*Y*w;	
						matrix[2]+=X*Z*w;	
						matrix[3]+=X*w;
	
						matrix[5]+=Y*Y*w;		
						matrix[6]+=Y*Z*w;
						matrix[7]+=Y*w;
	
						matrix[10]+=Z*Z*w;
						matrix[11]+=Z*w;
	
						matrix[15]+=w;
	
						Bxx[0]+=X*w*Txx;  Bxy[0]+=X*w*Txy;  Bxz[0]+=X*w*Txz;
						Byx[0]+=X*w*Tyx;  Byy[0]+=X*w*Tyy;  Byz[0]+=X*w*Tyz;
						Bzx[0]+=X*w*Tzx;  Bzy[0]+=X*w*Tzy;  Bzz[0]+=X*w*Tzz;
	
						Bxx[1]+=Y*w*Txx;  Bxy[1]+=Y*w*Txy;  Bxz[1]+=Y*w*Txz;
						Byx[1]+=Y*w*Tyx;  Byy[1]+=Y*w*Tyy;  Byz[1]+=Y*w*Tyz;
						Bzx[1]+=Y*w*Tzx;  Bzy[1]+=Y*w*Tzy;  Bzz[1]+=Y*w*Tzz;
	
						Bxx[2]+=Z*w*Txx;  Bxy[2]+=Z*w*Txy;  Bxz[2]+=Z*w*Txz;
						Byx[2]+=Z*w*Tyx;  Byy[2]+=Z*w*Tyy;  Byz[2]+=Z*w*Tyz;
						Bzx[2]+=Z*w*Tzx;  Bzy[2]+=Z*w*Tzy;  Bzz[2]+=Z*w*Tzz;
	
						Bxx[3]+=w*Txx;  Bxy[3]+=w*Txy;  Bxz[3]+=w*Txz;
						Byx[3]+=w*Tyx;  Byy[3]+=w*Tyy;  Byz[3]+=w*Tyz;
						Bzx[3]+=w*Tzx;  Bzy[3]+=w*Tzy;  Bzz[3]+=w*Tzz;
					
					}	
					matrix[4]=matrix[1];
					matrix[8]=matrix[2];
					matrix[9]=matrix[6];
					matrix[12]=matrix[3];
					matrix[13]=matrix[7];
					matrix[14]=matrix[11];

					for(int n=0;n<N*N;n++) matrix_val[n]=matrix[n];//�l��ۑ�
	
					//�K�E�X�̏����@�ŉ���
					gauss(matrix,Bxx,N);
					for(int n=0;n<N*N;n++) matrix[n]=matrix_val[n];//�l��߂�
					gauss(matrix,Bxy,N);
					for(int n=0;n<N*N;n++) matrix[n]=matrix_val[n];//�l��߂�
					gauss(matrix,Bxz,N);
					for(int n=0;n<N*N;n++) matrix[n]=matrix_val[n];//�l��߂�
	
					gauss(matrix,Byx,N);
					for(int n=0;n<N*N;n++) matrix[n]=matrix_val[n];//�l��߂�
					gauss(matrix,Byy,N);
					for(int n=0;n<N*N;n++) matrix[n]=matrix_val[n];//�l��߂�
					gauss(matrix,Byz,N);
					for(int n=0;n<N*N;n++) matrix[n]=matrix_val[n];//�l��߂�
	
					gauss(matrix,Bzx,N);
					for(int n=0;n<N*N;n++) matrix[n]=matrix_val[n];//�l��߂�
					gauss(matrix,Bzy,N);
					for(int n=0;n<N*N;n++) matrix[n]=matrix_val[n];//�l��߂�
					gauss(matrix,Bzz,N);
	
					double fx=Bxx[0]+Bxy[1]+Bxz[2];//�̐ϗ�fx:(divT)��X����
					double fy=Byx[0]+Byy[1]+Byz[2];//�̐ϗ�fy:(divT)��X����
					double fz=Bzx[0]+Bzy[1]+Bzz[2];//�̐ϗ�fz:(divT)��X����
	
					int p=NODE[i].particleID;//�Ή����闱�q�ԍ�
					if(p>=0)
					{
						F[A_X][p]=fx*V;
						F[A_Y][p]=fy*V;//�P�ʂ�[N]
						F[A_Z][p]=fz*V;
		
						
					}
					Fz+=fz*V;
				}
			}

		}
	}
	int *check=new int[fluid_number];//check[i]=ON�Ȃ炻�̗��q�͑Ή�����ߓ_�����݂���Ƃ�������
	for(int i=0;i<fluid_number;i++) check[i]=OFF;
	for(int i=1;i<=node;i++)
	{
		if(NODE[i].material==FLUID || NODE[i].material==ELASTIC || NODE[i].material==MAGELAST)//���̐ߓ_
		{
			int p=NODE[i].particleID;//���q�ԍ�
			if(p>=0) check[p]=ON;
		}
	}
	for(int i=0;i<fluid_number;i++)
	{
		if(check[i]==OFF)//�Ή�����ߓ_���Ȃ��̂ŁA�͂��[���̂܂܂ȗ��q
		{
			double W=0;
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				if(j<fluid_number)
				{
					if(check[j]==ON)
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);//���q�Ԃ̋���
						double w=R/dis-1;//�d�݊֐�
						W+=w;
						for(int D=0;D<3;D++) F[D][i]+=F[D][j]*w;
					}
				}
			}
			if(W!=0) for(int D=0;D<3;D++) F[D][i]/=W;
		}
	}

	delete [] matrix;
	delete [] matrix_val;
	delete [] Bxx;
	delete [] Bxy;
	delete [] Bxz;
	delete [] Byx;
	delete [] Byy;
	delete [] Byz;
	delete [] Bzx;
	delete [] Bzy;
	delete [] Bzz;

	delete [] check;

	cout<<"OK Fz="<<Fz<<"[N] time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
} 

//�P���r���͌v�Z�֐�
void kelvin_force3D(mpsconfig &CON, vector<mpselastic> &PART,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **Be,int fluid_number,double **F,double *RP,int*jnb,int **nei)
{
	cout<<"kelvin�͂ɂ��d���͌v�Z--";
	//kelvin�͂́�u0(M�E��)H dv
	unsigned int timeA=GetTickCount();
    double le=CON.get_distancebp();
	double R=CON.get_re()*le;
    double V=CON.get_particle_volume();//���q�̑̐�
    double u0=1.257e-6;				//�^��̓�����
	double kai=CON.get_RP()-1;
	double Ms=14700;					//�O�a����[A/m]
	double kai0=1.172;					//���C����
	double gamma=3*kai0/Ms;

	int N=4;						//���m��
	int order=1;					//�ߎ����x
	double Fz=0;
	double Fx=0;
    for(int i=0;i<fluid_number;i++) for(int D=0;D<3;D++) F[D][i]=0;//������

	double *matrix=new double [N*N];	//N�~N�̌W���s��
	double *matrix_val=new double [N*N];	//matrix�̕ۑ��p
	double *Bx=new double [N];	//N�̉��s��
	double *By=new double [N];	//N�̉��s��		//���Ƃ���By��Hy�̔��������߂�̂Ɏg��
	double *Bz=new double [N];	//N�̉��s��
	double *Bh=new double [N];	//N�̉��s��		//H�̑傫���̌��z�����߂�̂Ɏg��

	double H[3];

	double *direct[DIMENSION];
    for(int D=0;D<DIMENSION;D++) direct[D]=new double [fluid_number];
	double *fs=new double [fluid_number];//�\�ʐߓ_�̕\�ʗ͊i�[
	int *num=new int [fluid_number];//�\�ʐߓ_�̕\�ʗ͊i�[
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].surface==ON)  direct_f(CON,PART,i,direct);
		else  for(int D=0;D<DIMENSION;D++) direct[D][i]=0;
		fs[i]=0;
		num[i]=0;
	}

	
	//H�̌��z��WLSM�ŋ��߂�
	if(order==1)//3����1����
	{
		///�W���s���
		///   ����x2    ����x��y  ����x��z ����x  a = ����x��f  
		///  ����x��y   ����y2    ����y��z ����y  b = ����y��f 
		///  ����x��z   ����y��z  ����z2   ����z  c = ����z��f 
		///  ����x      ����y     ����z    ��1    d = ����f

		double Ri[3];//�ߓ_i�̍��W�i�[
		double Rg[3];	//�v�f�̏d�S�i�[
		for(int i=1;i<=node;i++)
		{
			if(NODE[i].material==FLUID || NODE[i].material==ELASTIC)
			{
				for(int n=0;n<N*N;n++) matrix[n]=0;		//������
				for(int n=0;n<N;n++)  {Bx[n]=0; By[n]=0; Bz[n]=0; Bh[n]=0;}
				
				for(int D=0;D<3;D++) Ri[D]=NODE[i].r[D];
				if(jnb[i]>5)
				{
					double Hi[3]={0,0,0};//�ߓ_i�ɂ�����H�i�[
					double HHi=0;		////�ߓ_i�ɂ�����H�̑傫���i�[
					for(int j=1;j<=jnb[i];j++)//�ߓ_i��jnb[n]�̗v�f�ɗאڂ��Ă���
					{
						int jelm=nei[i][j];
						double u=RP[jelm]*u0;
						double HHj=0;			//jelm�ɂ����鎥��̑傫��
						for(int D=0;D<3;D++) 
						{
							H[D]=Be[D][jelm]/u;//jelm�ɂ�����g�[�^���Ȏ���@�{���͊O�����ꂶ��Ȃ��Ƃ����Ȃ��炵�����ǁH
							Hi[D]+=H[D];
							HHj+=H[D]*H[D];
						}
						HHj=sqrt(HHj);
						HHi+=HHj;
						
						for(int D=0;D<3;D++) Rg[D]=0;
						for(int k=1;k<=4;k++)
						{
							int node=ELEM[jelm].node[k];//jelm���\������ߓ_�ԍ�
							for(int D=0;D<3;D++) Rg[D]+=NODE[node].r[D]*0.25;
						}//jelm�̏d�SRg[D]�����Ƃ܂����B
	
						double X=Rg[A_X]-Ri[A_X];
						double Y=Rg[A_Y]-Ri[A_Y];
						double Z=Rg[A_Z]-Ri[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
					
						double w=1;
						if(dis>le) w=le*le/(dis*dis);
		
						matrix[0]+=X*X*w;		
						matrix[1]+=X*Y*w;	
						matrix[2]+=X*Z*w;	
						matrix[3]+=X*w;
	
						matrix[5]+=Y*Y*w;		
						matrix[6]+=Y*Z*w;
						matrix[7]+=Y*w;
	
						matrix[10]+=Z*Z*w;
						matrix[11]+=Z*w;
	
						matrix[15]+=w;
	
						Bx[0]+=X*w*H[A_X];
						By[0]+=X*w*H[A_Y];
						Bz[0]+=X*w*H[A_Z]; 
						Bh[0]+=X*w*HHj;

						Bx[1]+=Y*w*H[A_X];
						By[1]+=Y*w*H[A_Y]; 
						Bz[1]+=Y*w*H[A_Z];
						Bh[1]+=Y*w*HHj;
	
						Bx[2]+=Z*w*H[A_X]; 
						By[2]+=Z*w*H[A_Y]; 
						Bz[2]+=Z*w*H[A_Z]; 
						Bh[2]+=Z*w*HHj;

						Bx[3]+=w*H[A_X];
						By[3]+=w*H[A_Y]; 
						Bz[3]+=w*H[A_Z]; 
						Bh[3]+=w*HHj;	
					}
					matrix[4]=matrix[1];
					matrix[8]=matrix[2];
					matrix[9]=matrix[6];
					matrix[12]=matrix[3];
					matrix[13]=matrix[7];
					matrix[14]=matrix[11];

					for(int D=0;D<3;D++) Hi[D]/=jnb[i];//�ߓ_i�ɂ�����H
					HHi/=jnb[i];
				
					matrix[15]+=1;//���g
					Bx[3]+=Hi[A_X];
					By[3]+=Hi[A_Y];
					Bz[3]+=Hi[A_Z];
					Bh[3]+=HHi;

					for(int n=0;n<N*N;n++) matrix_val[n]=matrix[n];//�l��ۑ�

					//�K�E�X�̏����@�ŉ���
					gauss(matrix,Bx,N);
				
					for(int n=0;n<N*N;n++) matrix[n]=matrix_val[n];//�l��߂�
				
					gauss(matrix,By,N);
				
					for(int n=0;n<N*N;n++) matrix[n]=matrix_val[n];//�l��߂�
				
					gauss(matrix,Bz,N);

					for(int n=0;n<N*N;n++) matrix[n]=matrix_val[n];//�l��߂�
				
					gauss(matrix,Bh,N);
				
					double Hix=Bx[3];//�ߓ_i�ɂ�����Hx
					double Hiy=By[3];//�ߓ_i�ɂ�����Hy
					double Hiz=Bz[3];//�ߓ_i�ɂ�����Hz
					HHi=Bh[3];		//WLSM�ɂ�苁�܂���H���i�[���Ȃ���

					double Mx,My,Mz,MM;
					if(CON.get_NLMH()==OFF)
					{
						MM=kai*HHi;
						Mx=kai*Hix;
						My=kai*Hiy;
						Mz=kai*Hiz;
					}
					else if(CON.get_NLMH()==ON)
					{
						double HH2=sqrt(Hix*Hix+Hiy*Hiy+Hiz*Hiz);//H�̑傫��
						MM=Ms*(cosh(gamma*HH2)/sinh(gamma*HH2)-1/(gamma*HH2));//�����W���o���̎��ɂ��M�̑傫��
						Mx=MM*Hix/HH2;
						My=MM*Hiy/HH2;
						Mz=MM*Hiz/HH2;
					}

				//	double fx=u0*(Mx*Bx[0]+My*Bx[1]+Mz*Bx[2]);//�̐ϗ�fx
				//	double fy=u0*(Mx*By[0]+My*By[1]+Mz*By[2]);//�̐ϗ�fy
				//	double fz=u0*(Mx*Bz[0]+My*Bz[1]+Mz*Bz[2]);//�̐ϗ�fz
	
					double fx=u0*MM*Bh[0];//�̐ϗ�fx
					double fy=u0*MM*Bh[1];//�̐ϗ�fy
					double fz=u0*MM*Bh[2];//�̐ϗ�fz
					
					int p=NODE[i].particleID;//�Ή����闱�q�ԍ�
				
					if(p>=0)//�Ή����闱�q������Ȃ�F[D][p]�ɑ��
					{
						if(PART[p].surface==OFF)//�������̗��q�̂ݗ͂��v�Z
						{
							F[A_X][p]=fx*V;
							F[A_Y][p]=fy*V;//�P�ʂ�[N]
							F[A_Z][p]=fz*V;
						}
					}
					Fz+=fz*V;
				}
			}
		}
	}///////////*/

	//Pm=u0��MdH�̌��z��WLSM�ŋ��߂�
	if(order==1)//3����1����
	{
		///�W���s���
		///   ����x2    ����x��y  ����x��z ����x  a = ����x��f  
		///  ����x��y   ����y2    ����y��z ����y  b = ����y��f 
		///  ����x��z   ����y��z  ����z2   ����z  c = ����z��f 
		///  ����x      ����y     ����z    ��1    d = ����f

		double Ri[3];//�ߓ_i�̍��W�i�[
		double Rg[3];	//�v�f�̏d�S�i�[
		for(int i=1;i<=node;i++)
		{
			if(NODE[i].material==FLUID || NODE[i].material==ELASTIC)
			{
				for(int n=0;n<N*N;n++) matrix[n]=0;		//������
				for(int n=0;n<N;n++)  {Bx[n]=0; By[n]=0; Bz[n]=0;}
				
				for(int D=0;D<3;D++) Ri[D]=NODE[i].r[D];
				if(jnb[i]>5)
				{
					double Hi[3]={0,0,0};//�ߓ_i�ɂ�����H�i�[
					int num_nei=0;			//�v�Z�Ɋ�^����v�f���J�E���g(��C�͏���)
					for(int j=1;j<=jnb[i];j++)//�ߓ_i��jnb[n]�̗v�f�ɗאڂ��Ă���
					{
						int jelm=nei[i][j];
						//if(ELEM[jelm].material==WATER)
						{
							num_nei++;
							double u=RP[jelm]*u0;
							for(int D=0;D<3;D++) 
							{
								H[D]=Be[D][jelm]/u;//jelm�ɂ�����g�[�^���Ȏ���@�{���͊O�����ꂶ��Ȃ��Ƃ����Ȃ��炵�����ǁH
								Hi[D]+=H[D];
							}
							double HHj=sqrt(H[A_X]*H[A_X]+H[A_Y]*H[A_Y]+H[A_Z]*H[A_Z]);//jelm��H�̑傫��
							double integ_of_M=0;//u0��MdH
							//if(ELEM[jelm].material==WATER)
							{
								if(CON.get_NLMH()==ON)	integ_of_M=u0*Ms/gamma*(log(sinh(gamma*HHj))-log(gamma*HHj));//�����W���o���̎�
								else integ_of_M=u0*kai*HHj*HHj*0.5;
							}
							
							for(int D=0;D<3;D++) Rg[D]=0;
							for(int k=1;k<=4;k++)
							{
								int node=ELEM[jelm].node[k];//jelm���\������ߓ_�ԍ�
								for(int D=0;D<3;D++) Rg[D]+=NODE[node].r[D]*0.25;
							}//jelm�̏d�SRg[D]�����Ƃ܂����B
		
							double X=Rg[A_X]-Ri[A_X];
							double Y=Rg[A_Y]-Ri[A_Y];
							double Z=Rg[A_Z]-Ri[A_Z];
							double dis=sqrt(X*X+Y*Y+Z*Z);
						
							double w=1;
							if(dis>le) w=le*le/(dis*dis);
			
							matrix[0]+=X*X*w;		
							matrix[1]+=X*Y*w;	
							matrix[2]+=X*Z*w;	
							matrix[3]+=X*w;
		
							matrix[5]+=Y*Y*w;		
							matrix[6]+=Y*Z*w;
							matrix[7]+=Y*w;
	
							matrix[10]+=Z*Z*w;
							matrix[11]+=Z*w;
		
							matrix[15]+=w;
		
							Bx[0]+=X*w*integ_of_M;
		
							Bx[1]+=Y*w*integ_of_M;
		
							Bx[2]+=Z*w*integ_of_M; 
		
							Bx[3]+=w*integ_of_M;
						}
					}
					matrix[4]=matrix[1];
					matrix[8]=matrix[2];
					matrix[9]=matrix[6];
					matrix[12]=matrix[3];
					matrix[13]=matrix[7];
					matrix[14]=matrix[11];
		
					if(num_nei>0) for(int D=0;D<3;D++) Hi[D]/=num_nei;//�ߓ_i�ɂ�����H
					double HHi=sqrt(Hi[A_X]*Hi[A_X]+Hi[A_Y]*Hi[A_Y]+Hi[A_Z]*Hi[A_Z]);//i��H�̑傫��
					double integ_of_Mi;//i��u0��MdH
					if(CON.get_NLMH()==ON) integ_of_Mi=u0*Ms/gamma*(log(sinh(gamma*HHi))-log(gamma*HHi));//�����W���o���̎�
					else integ_of_Mi=0.5*u0*HHi*HHi*kai;
	
					matrix[15]+=1;//���g
					Bx[3]+=integ_of_Mi;
	
					for(int n=0;n<N*N;n++) matrix_val[n]=matrix[n];//�l��ۑ�
	
					//�K�E�X�̏����@�ŉ���
					gauss(matrix,Bx,N);
					
					for(int n=0;n<N*N;n++) matrix[n]=matrix_val[n];//�l��߂�
					
					integ_of_Mi=Bx[3];//�ߓ_i�ɂ�����u0��MdH
	
					double fx=Bx[0];//�̐ϗ�fx
					double fy=Bx[1];//�̐ϗ�fy
					double fz=Bx[2];//�̐ϗ�fz
				
					int p=NODE[i].particleID;
					
					if(p>=0)//�Ή����闱�q������Ȃ�F[D][p]�ɑ��
					{
						if(PART[p].surface==OFF)//�������̗��q�̂ݗ͂��v�Z
						{
							F[A_X][p]-=fx*V;
							F[A_Y][p]-=fy*V;//�P�ʂ�[N]
							F[A_Z][p]-=fz*V;
						}
					}
					Fz+=fz*V;
				}
			}
		}
	}///////////////////*/

	//�\�ʐϕ����s
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==FLUID || ELEM[i].material==ELASTIC)
		{
			double u=RP[i]*u0;
			for(int j=1;j<=4;j++)
			{
				int jelm=ELEM[i].elm[j];
				if(ELEM[jelm].material==AIR)
				{
					int ia=ELEM[i].node[j%4+1];//���_���̌����Ɉʒu����O�p�`�̒��_ ia��ib��ic�͂�����݂Ĕ����v�܂��
					int ib=ELEM[i].node[4-(j-1)/2*2];
					int ic=ELEM[i].node[3-(j/2%2)*2];

					double iaic[3];//ia��ic�̃x�N�g�������i�[
					double iaib[3];//ia��ib�̃x�N�g�������i�[
					for(int D=0;D<3;D++)
					{
						iaic[D]=NODE[ic].r[D]-NODE[ia].r[D];
						iaib[D]=NODE[ib].r[D]-NODE[ia].r[D];
					}
					///0.5(iaic�~iaib)�͎O�p�`ia,ib,ic�̖ʐς̑傫���������A�����͊O�����̃x�N�g���ƂȂ�
					double S[3];//��L�̃x�N�g�������i�[
					S[A_X]=0.5*(iaic[A_Y]*iaib[A_Z]-iaic[A_Z]*iaib[A_Y]);
					S[A_Y]=0.5*(iaic[A_Z]*iaib[A_X]-iaic[A_X]*iaib[A_Z]);
					S[A_Z]=0.5*(iaic[A_X]*iaib[A_Y]-iaic[A_Y]*iaib[A_X]);
					double SS=sqrt(S[A_X]*S[A_X]+S[A_Y]*S[A_Y]+S[A_Z]*S[A_Z]);
					////�ʐ�S�����Ƃ܂���

					double n[3]={S[A_X]/SS,S[A_Y]/SS,S[A_Z]/SS};//�@���x�N�g��
					double H[3]={Be[A_X][i]/u,Be[A_Y][i]/u,Be[A_Z][i]/u};//H
					double ave_B[3]={(Be[A_X][i]+Be[A_X][jelm])*0.5,(Be[A_Y][i]+Be[A_Y][jelm])*0.5,(Be[A_Z][i]+Be[A_Z][jelm])*0.5};//�E�ʂ�B
					double H2[3]={ave_B[A_X]/u0,ave_B[A_Y]/u0,ave_B[A_Z]/u0};//��C����H
					double H1[3]={ave_B[A_X]/u,ave_B[A_Y]/u,ave_B[A_Z]/u};//���̑���H
					double Hn2=H2[A_X]*n[A_X]+H2[A_Y]*n[A_Y]+H2[A_Z]*n[A_Z];//H2�̖@������
					double Hn1=H1[A_X]*n[A_X]+H1[A_Y]*n[A_Y]+H1[A_Z]*n[A_Z];
					double Mn=Hn2-Hn1;//ferrohydrodynamics�̎�5.21b�Q��
					double BB=sqrt(ave_B[A_X]*ave_B[A_X]+ave_B[A_Y]*ave_B[A_Y]+ave_B[A_Z]*ave_B[A_Z]);//�E�ʂ�B�̑傫��
					double Bn=ave_B[A_X]*n[A_X]+ave_B[A_Y]*n[A_Y]+ave_B[A_Z]*n[A_Z];
					
					double M[3];//����
					for(int D=0;D<3;D++) M[D]=0.5*(H1[D]+H2[D])*kai;
					double MM=sqrt(M[A_X]*M[A_X]+M[A_Y]*M[A_Y]+M[A_Z]*M[A_Z]);//�����̑傫��
					double HH=MM/kai;
					/*if(CON.get_NLMH()==OFF)
					{
						MM=kai*HH;
						for(int D=0;D<3;D++) M[D]=H[D]*kai;
					}
					else
					{
						MM=Ms*(cosh(gamma*HH)/sinh(gamma*HH)-1/(gamma*HH));//�����W���o���̎��ɂ��M�̑傫��
						for(int D=0;D<3;D++) M[D]=MM*H[D]/HH;
					}*/
					//double Mn=0;
					//for(int D=0;D<3;D++) Mn+=M[D]*n[D];//�@�������̎���
					//cout<<Mn<<endl;*/
					double Fs[3]={0,0,0};
					for(int D=0;D<3;D++) Fs[D]=0.5*u0*Mn*Mn*SS*n[D];//�\�ʂɓ����� //�P�ʂ�[N]
					//for(int D=0;D<3;D++) Fs[D]+=0.5*u0*MM*MM*SS*n[D];//???
					//for(int D=0;D<3;D++) Fs[D]+=0.5*u0*kai*HH*HH*SS*n[D];//???

					double integ_of_Mi;//u0��MdH
					if(CON.get_NLMH()==ON) integ_of_Mi=u0*Ms/gamma*(log(sinh(gamma*HH))-log(gamma*HH));//�����W���o���̎�
					else integ_of_Mi=0.5*u0*HH*HH*kai;
					for(int D=0;D<3;D++) Fs[D]+=integ_of_Mi*SS*n[D];//�\�ʂɓ����� //�P�ʂ�[N]

					Fz+=Fs[A_Z];

					//����ꂽFs�𒸓_�Ɋ��蓖�Ă�
					int node[3]={ia,ib,ic};//���_�̔ԍ����L��
					for(int k=0;k<3;k++)
					{
						int node_ID=node[k];//�ߓ_�ԍ�
						if(NODE[node_ID].particleID>=0)//���̐ߓ_�Ȃ�
						{
							int p=NODE[node_ID].particleID;
							for(int D=0;D<3;D++) F[D][p]+=Fs[D]/3;
							for(int D=0;D<3;D++) fs[p]+=Fs[D]/SS*n[D];
							num[p]++;
						}
					}
				}
			}
		}
	}///*/

	//���C�o�Ȏq�Ԃ̗��q�ԗ͂��v�Z
	//coroid_pole(CON,PART,NODE,ELEM, node, nelm,Be, fluid_number,F,RP,jnb,nei);

	int *check=new int[fluid_number];//check[i]=ON�Ȃ炻�̗��q�͑Ή�����ߓ_�����݂���Ƃ�������
	for(int i=0;i<fluid_number;i++) check[i]=OFF;
	for(int I=1;I<=node;I++)
	{
		if(NODE[I].material==FLUID || NODE[I].material==ELASTIC)
		{	
			int i=NODE[I].particleID;//���q�ԍ�
			if(i>=0) check[i]=ON;
		}
	}
	
	for(int i=0;i<fluid_number;i++)
	{
		if(check[i]==OFF)//�Ή�����ߓ_���Ȃ��̂ŁA�͂��[���̂܂܂ȗ��q
		{
			double W=0;
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				if(j<fluid_number)
				{
					if(check[j]==ON)
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);//���q�Ԃ̋���
						double w=R/dis-1;//�d�݊֐�
						W+=w;
						for(int D=0;D<3;D++) F[D][i]+=F[D][j]*w;
					}
				}
			}
			if(W!=0) for(int D=0;D<3;D++) F[D][i]/=W;
		}
	}

	if(CON.get_dir_for_P()==2 ||CON.get_dir_for_P()==3 )
    {
		for(int i=0;i<fluid_number;i++)
		{
			double val=0;//�\�ʗ�
			if(PART[i].surface==ON)
			{
				if(num[i]>0) val=fs[i]/num[i];
				else	//�f���[�j������jnb=0�ƂȂ������q�͂����ł�num[i]=0�ƂȂ��Ă��܂��f�B���N���l�����܂�Ȃ��B���������Ƃ��͎��ӗ��q���l���Ԃ���
				{
					double count=0;
					for(int k=0;k<PART[i].N;k++)
					{
						int j=PART[i].NEI[k];
						if(j<fluid_number && PART[j].surface==ON)
						{
							if(num[j]>0)
							{
								val+=fs[j];
								count+=1;
							}
						}
					}
					if(count>0) val/=count;		//������count=0�ƂȂ�悤�Ȃ炻�̗��q�͊��S�ɌǗ����q�Ȃ̂ŁA�f�B���N���l�̓[�����i�[�����Ă����΂悢
				}
				for(int D=0;D<3;D++) F[D][i]=0;//���̓fިظڂƂ��ēd���͂��g�p����̂ł����ł͏�����
			}
			PART[i].dir_Pem=-val;
        }
	}
	ofstream fg("Fss.dat");
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].surface==ON)
		{
			fg<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].dir_Pem*-1<<endl;
		}
	}

	delete [] matrix;
	delete [] matrix_val;
	delete [] Bx;
	delete [] By;
	delete [] Bz;
	delete [] Bh;

	delete [] check;
	for(int D=0;D<DIMENSION;D++) delete [] direct[D];
	delete [] fs;
	delete [] num;

	cout<<"OK Fz="<<Fz<<"[N] time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
	cout<<"Fx="<<Fx<<endl;
}

void coroid_pole(mpsconfig &CON, vector<mpselastic> &PART,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **Be,int fluid_number,double **F,double *RP,int*jnb,int **nei)
{
	double le=CON.get_distancebp();
    double V=CON.get_particle_volume();//���q�̑̐�
	double u0=1.257e-6;	
	double u=CON.get_RP()*u0;
	double kai=CON.get_RP()-1;

	double *M[DIMENSION];
	for(int D=0;D<DIMENSION;D++) M[D]=new double [fluid_number];//�e���q�ʒu�ł�M
	double *F_di[DIMENSION];
	for(int D=0;D<DIMENSION;D++) F_di[D]=new double [fluid_number];//���̊֐��Ōv�Z�����F
	double *B_di[DIMENSION];
	for(int D=0;D<DIMENSION;D++) B_di[D]=new double [fluid_number];//���̊֐��Ōv�Z�����B
	int *check=new int[fluid_number];//check=ON�̗��q�͑Ή�����ߓ_�����݂���B�䂦�ɑΉ�����M����ɓ���Bcheck=OFF�Ȃ���͂���M�̒l���Ԃ��Ȃ���΂Ȃ�Ȃ�

	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<DIMENSION;D++)
		{
			M[D][i]=0;//������
			F_di[D][i]=0;//������
			B_di[D][i]=0;//������
		}
		check[i]=OFF;
	}

	for(int i=1;i<=node;i++)
	{
		if(NODE[i].particleID>=0)//�Ή����闱�q������Ȃ�
		{
			int ip=NODE[i].particleID;
			int num=0;
			check[ip]=ON;				//�Ή�����ߓ_�����݂���
			for(int j=1;j<=jnb[i];j++)
			{
				int jelm=nei[i][j];//�אڂ���v�f�ԍ�
				if(ELEM[jelm].material==FLUID || ELEM[jelm].material==ELASTIC)
				{
					double H[3];
					for(int D=0;D<3;D++) H[D]=Be[D][jelm]/u;
					for(int D=0;D<3;D++) M[D][ip]+=H[D]*kai;
					num++;
				}
			}
			if(num>0) for(int D=0;D<3;D++) M[D][ip]/=num;
		}
	}

	for(int i=0;i<fluid_number;i++)
	{
		if(check[i]==OFF)
		{
			int num=0;
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				if(PART[j].type==FLUID || PART[j].type==ELASTIC) if(check[j]==ON) for(int D=0;D<3;D++) M[D][i]+=M[D][j];
			}
			if(num>0) for(int D=0;D<3;D++) M[D][i]/=num;
			check[i]=ON;
		}
	}//�S���q��M�����܂���

	for(int i=0;i<fluid_number;i++) for(int D=0;D<3;D++) M[D][i]*=V*0.32;//���C���[�����gm�ɒu������

	//�͂��v�Z
	double co=u0/(4*PI);//�悭�g���W��
	double f[3]={0,0,0};
	double L=le/10;
	for(int i=0;i<fluid_number;i++)
	{
		for(int k=0;k<PART[i].N3;k++)
		{
			int j=PART[i].NEI3[k];
			if(PART[j].type==FLUID || PART[j].type==ELASTIC)
			{	
				/*double dX[3];
				for(int D=0;D<3;D++) dX[D]=PART[i].r[D]-PART[j].r[D];//j����i�֌������x�N�g��
				double r=sqrt(dX[A_X]*dX[A_X]+dX[A_Y]*dX[A_Y]+dX[A_Z]*dX[A_Z]);//i��j�̋���
			
				double MMi=sqrt(M[A_X][i]*M[A_X][i]+M[A_Y][i]*M[A_Y][i]+M[A_Z][i]*M[A_Z][i]);//Mi�̑傫��
				double MMj=sqrt(M[A_X][j]*M[A_X][j]+M[A_Y][j]*M[A_Y][j]+M[A_Z][j]*M[A_Z][j]);//Mj�̑傫��
				double cos1=(dX[A_X]*M[A_X][i]+dX[A_Y]*M[A_Y][i]+dX[A_Z]*M[A_Z][i])/(r*MMi);		//r��Mi�̂Ȃ��p�x
				double cos2=(dX[A_X]*M[A_X][j]+dX[A_Y]*M[A_Y][j]+dX[A_Z]*M[A_Z][j])/(r*MMj);		//r��Mj�̂Ȃ��p�x
				double cos3=(M[A_X][i]*M[A_X][j]+M[A_Y][i]*M[A_Y][j]+M[A_Z][i]*M[A_Z][j])/(MMj*MMi);		//Mj��Mi�̂Ȃ��p�x
				for(int D=0;D<3;D++) F_di[D][i]+=-co*3*MMi*MMj/pow(r,4.0)*(3*cos1*cos2-cos3)*dX[D]/r;
				//for(int D=0;D<3;D++) F_di[D][i]+=co*3*MMi*MMj/pow(r,4.0)*(cos1+cos2-5*cos1*cos2+cos3)*dX[D]/r;
				for(int D=0;D<3;D++) f[D]+=-co*3*MMi*MMj/pow(r,4.0)*(3*cos1*cos2-cos3)*dX[D]/r;*/

				double Xnp[6];
				double Ynp[6];
				double Znp[6];
				double Bx[6],By[6],Bz[6];
				Xnp[0]=PART[i].r[A_X]-L; Xnp[1]=PART[i].r[A_X]+L; Xnp[2]=PART[i].r[A_X]; Xnp[3]=PART[i].r[A_X]; Xnp[4]=PART[i].r[A_X]; Xnp[5]=PART[i].r[A_X];
				Ynp[0]=PART[i].r[A_Y]; Ynp[1]=PART[i].r[A_Y]; Ynp[2]=PART[i].r[A_Y]-L; Ynp[3]=PART[i].r[A_Y]+L; Ynp[4]=PART[i].r[A_Y]; Ynp[5]=PART[i].r[A_Y]; 
				Znp[0]=PART[i].r[A_Z]; Znp[1]=PART[i].r[A_Z]; Znp[2]=PART[i].r[A_Z]; Znp[3]=PART[i].r[A_Z]; Znp[4]=PART[i].r[A_Z]-L; Znp[5]=PART[i].r[A_Z]+L; 

				for(int k=0;k<6;k++)
				{
					double dX[3];
					dX[A_X]=Xnp[k]-PART[j].r[A_X];//j����i�֌������x�N�g��
					dX[A_Y]=Ynp[k]-PART[j].r[A_Y];//j����i�֌������x�N�g��
					dX[A_Z]=Znp[k]-PART[j].r[A_Z];//j����i�֌������x�N�g��
					double r=sqrt(dX[A_X]*dX[A_X]+dX[A_Y]*dX[A_Y]+dX[A_Z]*dX[A_Z]);//i��j�̋���
				
					Bx[k]=co*3*(dX[A_X]*M[A_X][j]+dX[A_Y]*M[A_Y][j]+dX[A_Z]*M[A_Z][j])*dX[A_X]/pow(r,5.0)-co*M[A_X][j]/pow(r,3.0);
					By[k]=co*3*(dX[A_X]*M[A_X][j]+dX[A_Y]*M[A_Y][j]+dX[A_Z]*M[A_Z][j])*dX[A_Y]/pow(r,5.0)-co*M[A_Y][j]/pow(r,3.0);
					Bz[k]=co*3*(dX[A_X]*M[A_X][j]+dX[A_Y]*M[A_Y][j]+dX[A_Z]*M[A_Z][j])*dX[A_Z]/pow(r,5.0)-co*M[A_Z][j]/pow(r,3.0);
				}

				F_di[A_X][i]+=M[A_X][i]*(Bx[1]-Bx[0])/(2*L)+M[A_Y][i]*(Bx[3]-Bx[2])/(2*L)+M[A_Z][i]*(Bx[5]-Bx[4])/(2*L);
				F_di[A_Y][i]+=M[A_X][i]*(By[1]-By[0])/(2*L)+M[A_Y][i]*(By[3]-By[2])/(2*L)+M[A_Z][i]*(By[5]-By[4])/(2*L);
				F_di[A_Z][i]+=M[A_X][i]*(Bz[1]-Bz[0])/(2*L)+M[A_Y][i]*(Bz[3]-Bz[2])/(2*L)+M[A_Z][i]*(Bz[5]-Bz[4])/(2*L);

				//B�����߂�
				double dX[3];
				for(int D=0;D<3;D++) dX[D]=PART[i].r[D]-PART[j].r[D];//j����i�֌������x�N�g��
				double r=sqrt(dX[A_X]*dX[A_X]+dX[A_Y]*dX[A_Y]+dX[A_Z]*dX[A_Z]);//i��j�̋���
				
				for(int D=0;D<3;D++) B_di[D][i]+=co*3*(dX[A_X]*M[A_X][j]+dX[A_Y]*M[A_Y][j]+dX[A_Z]*M[A_Z][j])*dX[D]/pow(r,5.0)-co*M[D][j]/pow(r,3.0);

			}
		}
	}
	//F[D][i]�ɉ�����
	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<3;D++)
		{
			//F[D][i]+= F_di[D][i];
			f[D]+=F_di[D][i];
		}
	}

	cout<<f[A_X]<<" "<<f[A_Y]<<" "<<f[A_Z]<<endl;

	ofstream fp("M.dat");
	ofstream fq("F_di.dat");
	ofstream fr("B_di.dat");
	double times=1000;
	double timesF=CON.get_times()/CON.get_density()/le*CON.get_FEMtimes();
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].r[A_Y]<le && PART[i].r[A_Y]>-le)
		{
			fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<M[A_X][i]*times<<" "<<M[A_Z][i]*times<<endl;
			fq<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<F_di[A_X][i]*times<<" "<<F_di[A_Z][i]*times<<endl;
			fr<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<B_di[A_X][i]*CON.get_B_times()<<" "<<B_di[A_Z][i]*CON.get_B_times()<<endl;
		}
	}
	fp.close();
	fq.close();
	fr.close();

	for(int D=0;D<DIMENSION;D++) delete [] M[D];
	for(int D=0;D<DIMENSION;D++) delete [] F_di[D];
	for(int D=0;D<DIMENSION;D++) delete [] B_di[D];
	delete [] check;
}

//�ߓ_�͖@�v�Z�֐�
void NODE_F3D(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **Ee,int *jnb,int **nei,double *RP,vector<mpselastic> &PART,double **F,int fluid_number)
{
	cout<<"�ߓ_�͖@�ɂ��d���͌v�Z--------";
	double ep0=8.854e-12;	//�^��̗U�d���B
	double u0=12.5e-7;		//�^��̓�����
	int N[4+1];				//�v�f�̊e�ߓ_�ԍ��i�[
	double X[4+1];
	double Y[4+1];
	double Z[4+1];
	double Fz=0;			//Z�����̍���[N]
	double Fx=0;
	unsigned timeA=GetTickCount();//�v�Z�J�n����

	double *Fn[3];
	for(int D=0;D<3;D++) Fn[D]=new double [node+1];//

	for(int i=0;i<=node;i++) for(int D=0;D<3;D++) Fn[D][i]=0;//������

	if(CON.get_EM_calc_type()==1)//�Ód�͌v�Z
	{
		for(int I=1;I<=node;I++)
		{
			if(NODE[I].material==FLUID || NODE[I].material==MAGELAST  || NODE[I].material==ELASTIC || NODE[I].material==IRON)
			{			
				for(int k=1;k<=jnb[I];k++)
				{
					int jelm=nei[I][k];//�ߓ_I���אڂ���v�f�ԍ�
				
					double ep=ep0;
					if(ELEM[jelm].material==FLUID) ep*=CON.get_r_perm();

					//�}�N�X�E�F���̉��̓e���\��
					double Txx=0.5*ep*(Ee[A_X][jelm]*Ee[A_X][jelm]-Ee[A_Y][jelm]*Ee[A_Y][jelm]-Ee[A_Z][jelm]*Ee[A_Z][jelm]);
					double Txy=ep*Ee[A_X][jelm]*Ee[A_Y][jelm];
					double Txz=ep*Ee[A_X][jelm]*Ee[A_Z][jelm];
					double Tyx=Txy;
					double Tyy=0.5*ep*(-Ee[A_X][jelm]*Ee[A_X][jelm]+Ee[A_Y][jelm]*Ee[A_Y][jelm]-Ee[A_Z][jelm]*Ee[A_Z][jelm]);
					double Tyz=ep*Ee[A_Y][jelm]*Ee[A_Z][jelm];
					double Tzx=Txz;
					double Tzy=Tyz;
					double Tzz=0.5*ep*(-Ee[A_X][jelm]*Ee[A_X][jelm]-Ee[A_Y][jelm]*Ee[A_Y][jelm]+Ee[A_Z][jelm]*Ee[A_Z][jelm]);
					//////////
		   
					/////�W��c,d,e�v�Z
					for(int j=1;j<=4;j++)
					{
						N[j]=ELEM[jelm].node[j];
						X[j]=NODE[N[j]].r[A_X];
						Y[j]=NODE[N[j]].r[A_Y];
						Z[j]=NODE[N[j]].r[A_Z];
					}

					int i=0;///�ߓ_i�͗v�fjelm�̑�j�Ԗڂ̐ߓ_
					for(int j=1;j<=4;j++) if(ELEM[jelm].node[j]==I) i=j;
					if(i==0) cout<<"ERROR IN Fn"<<endl;

					int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
					int m=j%4+1;
					int n=m%4+1;
					//delta6�͑��E�����̂ł���Ȃ�
					double c=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
					double d=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
					double e=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
					if(i%2!=0)//i����Ȃ�
					{
						c*=-1;
						d*=-1;
						e*=-1;
					}
					///////////////
		
					Fn[A_X][I]+=Txx*c+Txy*d+Txz*e;
					Fn[A_Y][I]+=Tyx*c+Tyy*d+Tyz*e;
					Fn[A_Z][I]+=Tzx*c+Tzy*d+Tzz*e;
				}
				for(int D=0;D<3;D++) Fn[D][I]*=-1.00000000000/6.00000000000;
				Fz+=Fn[A_Z][I];
				
			}
		}
	}
	else if(CON.get_EM_calc_type()==4)//���ʉ�� Ee�ɂ͂g���i�[����Ă��邱�Ƃɒ���
	{
		for(int I=1;I<=node;I++)
		{
			if(NODE[I].material==FLUID || NODE[I].material==MAGELAST || NODE[I].material==ELASTIC || NODE[I].material==IRON)
			{
				for(int k=1;k<=jnb[I];k++)
				{
					int jelm=nei[I][k];//�ߓ_i���אڂ���v�f�ԍ�
					//if(ELEM[jelm].material==AIR)
					{
						double u=u0*RP[jelm];

						//�}�N�X�E�F���̉��̓e���\��
						double Txx=0.5*u*(Ee[A_X][jelm]*Ee[A_X][jelm]-Ee[A_Y][jelm]*Ee[A_Y][jelm]-Ee[A_Z][jelm]*Ee[A_Z][jelm]);
						double Txy=u*Ee[A_X][jelm]*Ee[A_Y][jelm];
						double Txz=u*Ee[A_X][jelm]*Ee[A_Z][jelm];
						double Tyx=Txy;
						double Tyy=0.5*u*(-Ee[A_X][jelm]*Ee[A_X][jelm]+Ee[A_Y][jelm]*Ee[A_Y][jelm]-Ee[A_Z][jelm]*Ee[A_Z][jelm]);
						double Tyz=u*Ee[A_Y][jelm]*Ee[A_Z][jelm];
						double Tzx=Txz;
						double Tzy=Tyz;
						double Tzz=0.5*u*(-Ee[A_X][jelm]*Ee[A_X][jelm]-Ee[A_Y][jelm]*Ee[A_Y][jelm]+Ee[A_Z][jelm]*Ee[A_Z][jelm]);
						//*/
						
					//�W��c,d,e�v�Z
						for(int j=1;j<=4;j++)
						{
							N[j]=ELEM[jelm].node[j];
							X[j]=NODE[N[j]].r[A_X];
							Y[j]=NODE[N[j]].r[A_Y];
							Z[j]=NODE[N[j]].r[A_Z];
						}
						int i=0;///�ߓ_i�͗v�fjelm�̑�J�Ԗڂ̐ߓ_
						for(int j=1;j<=4;j++) if(ELEM[jelm].node[j]==I) i=j;

						//if(i==0) cout<<"ERROR IN Fn"<<endl;
						int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
						int m=j%4+1;
						int n=m%4+1;
						//delta6�͑��E�����̂ł���Ȃ�
						double c=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
						double d=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
						double e=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
						if( i & 1 )//i����Ȃ�
						{
							c*=-1;
							d*=-1;
							e*=-1;
						}
						///////////////
			
						
						
						Fn[A_X][I]+=(Txx*c+Txy*d+Txz*e);
						Fn[A_Y][I]+=(Tyx*c+Tyy*d+Tyz*e);
						Fn[A_Z][I]+=(Tzx*c+Tzy*d+Tzz*e);
						//cout<<(Tzx*c+Tzy*d+Tzz*e)<<" "<<c<<" "<<d<<" "<<e<<endl;
					}
				}
				for(int D=0;D<3;D++) Fn[D][I]*=-1.00000000000/6.00000000000;
				//if(Fn[A_Z][I]<0) Fn[A_Z][I]=0;
				Fz+=Fn[A_Z][I];
			}
		}
	}
	else //������
	{
	for(int I=1;I<=node;I++)
	{
		if(NODE[I].material==FLUID || NODE[I].material==ELASTIC || NODE[I].material==MAGELAST || NODE[I].material==IRON)
			{
				for(int k=1;k<=jnb[I];k++)
				{
					int jelm=nei[I][k];//�ߓ_i���אڂ���v�f�ԍ�
						
					//if(ELEM[jelm].material==AIR){
					//�}�N�X�E�F���̉��̓e���\��
					double Txx=(Ee[A_X][jelm]*Ee[A_X][jelm]-Ee[A_Y][jelm]*Ee[A_Y][jelm]-Ee[A_Z][jelm]*Ee[A_Z][jelm]);
					double Txy=2*Ee[A_X][jelm]*Ee[A_Y][jelm];
					double Txz=2*Ee[A_X][jelm]*Ee[A_Z][jelm];
					double Tyx=Txy;
					double Tyy=(-Ee[A_X][jelm]*Ee[A_X][jelm]+Ee[A_Y][jelm]*Ee[A_Y][jelm]-Ee[A_Z][jelm]*Ee[A_Z][jelm]);
					double Tyz=2*Ee[A_Y][jelm]*Ee[A_Z][jelm];
					double Tzx=Txz;
					double Tzy=Tyz;
					double Tzz=(-Ee[A_X][jelm]*Ee[A_X][jelm]-Ee[A_Y][jelm]*Ee[A_Y][jelm]+Ee[A_Z][jelm]*Ee[A_Z][jelm]);
					//////////
					
					//�W��c,d,e�v�Z
					for(int j=1;j<=4;j++)
					{
						N[j]=ELEM[jelm].node[j];
						X[j]=NODE[N[j]].r[A_X];
						Y[j]=NODE[N[j]].r[A_Y];
						Z[j]=NODE[N[j]].r[A_Z];
					}

					int i=0;//�ߓ_i�͗v�fjelm�̑�j�Ԗڂ̐ߓ_
					for(int j=1;j<=4;j++) if(ELEM[jelm].node[j]==I) i=j;
					int j=i%4+1;//i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
					int m=j%4+1;
					int n=m%4+1;

				//delta6�͑��E�����̂ł���Ȃ�
					double c=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
					double d=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
					double e=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
					if(i & 1)//i����Ȃ�
					{
						c*=-1;
						d*=-1;
						e*=-1;
					}
					///////////////
		
					double u=RP[jelm];
						
					Fn[A_X][I]+=(Txx*c+Txy*d+Txz*e)/(2*u);//u0�͂��̉��ł����Ă�
					Fn[A_Y][I]+=(Tyx*c+Tyy*d+Tyz*e)/(2*u);
					Fn[A_Z][I]+=(Tzx*c+Tzy*d+Tzz*e)/(2*u);
					//}
				}
				for(int D=0;D<3;D++) Fn[D][I]*=-1.00000000000/(6.00000000000*u0);
				//if(Fn[A_Z][I]<0) Fn[A_Z][I]=0;

				Fz+=Fn[A_Z][I];
				if(NODE[I].r[A_X]<0) Fx+= Fn[A_X][I];//����͉��H�H

			}//if(NODE[I].material==FLUID || NODE[I].material==ELASTIC)
		}
	}
	//cout<<"OK Fz="<<Fz<<"[N] time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

	int *bound=new int[node+1];
	for(int n=1;n<=node;n++)
	{
		bound[n]=OFF;
		if(NODE[n].material==FLUID || NODE[n].material==ELASTIC || NODE[n].material==MAGELAST|| NODE[n].material==IRON)
		{
			for(int k=1;k<=jnb[n];k++)
			{
				int jelm=nei[n][k];//�ߓ_i���אڂ���v�f�ԍ�
				if(ELEM[jelm].material==AIR) bound[n]=ON;
			}
		}
	}
	
	

	

	//�t�@�C���o��
	ofstream fp("Fn.dat");
	double le=CON.get_distancebp();
	double times=CON.get_times()/CON.get_density()/le*CON.get_FEMtimes();

	for(int i=1;i<=node;i++)//���̐ߓ_�̂ݏo��
	{
		if(NODE[i].material==FLUID || NODE[i].material==ELASTIC || NODE[i].material==IRON) if(NODE[i].r[A_Y]>-le*0.5&& NODE[i].r[A_Y]<le*0.5) fp<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Z]<<" "<<Fn[A_X][i]*times<<" "<<Fn[A_Z][i]*times<<endl;	
	}
	fp.close();//


	for(int D=0;D<3;D++) delete [] Fn[D];
	cout<<"OK Fz="<<Fz<<"[N] time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

	if(2*Fz!=Fz+Fz){			//�G���[���ł���Ferror��ON
		CON.set_Ferror(ON);
		cout<<"Fz�ŃG���[�A�O��̓d���͂��g�p���܂�";
	}
	else {
		cout<<"not error";
		CON.set_Ferror(OFF);	//�G���[���Ȃ��Ȃ�
		///////////F�X�V////////////
		for(int n=1;n<=node;n++)
		{
			if(bound[n]==ON)
			{
				int i=NODE[n].particleID;		
				if(i>=0) for(int D=0;D<3;D++) F[D][i]=Fn[D][n];
			}
		}
		//�ߓ_�͖@�̗͂̒P�ʂ�[N]�Ȃ̂ŁA�������炳��ɗ��q���ߓ_�̗��q�Ɋւ��Ă��͂����Ƃ߂Ă��K�v�͂Ȃ�
		cout<<"!"<<endl;
	}
	delete [] bound;
	cout<<"ok1"<<endl;
}

//�������x���v���b�g(�����쐬)
void plot_magnetic_flux_density(mpsconfig &CON, vector<mpselastic> &PART, vector<point3D> &NODE, vector<element3D> &ELEM, int nelm, double **B, int t)
{

	//�������x�X���[�W���O�֐�
	double le=CON.get_distancebp();
//�Ƃ肠�����X���[�W���O�������ɂ��

	/*  double *newB[3];
	for(int D=0;D<3;D++) newB[D]=new double [nelm+1];
	if(CON.get_FEM_smn()>0)
	{
		for(int n=0;n<CON.get_FEM_smn();n++)
		{
		    for(int i=1;i<=nelm;i++) 
		    {  
		        for(int D=0;D<3;D++) newB[D][i]=B[D][i];
				int num=1; //�������g���J�E���g���邩��1
				for(int k=0;k<PART[i].N;k++)
				{       
					int j=PART[i].NEI[k];
					if(PART[j].type==FLUID || PART[j].type==ELASTIC)
					{
						num++;
						for(int D=0;D<3;D++) newB[D][i]+=B[D][j];
					}
				}
				for(int D=0;D<3;D++) newB[D][i]/=num;
		    } 
		    for(int i=1;i<nelm;i++) for(int D=0;D<3;D++) B[D][i]=newB[D][i];
		}
	}
	else if(CON.get_FEM_smn()<0)//�\�ʂ݂̂ŃX���[�W���O
	{
		int N=-1*CON.get_FEM_smn();
		for(int n=0;n<N;n++)
		{
			for(int i=1;i<=nelm;i++) 
			{  
			    for(int D=0;D<3;D++) newB[D][i]=B[D][i];
				if(PART[i].surface==ON)
				{
					int num=1; //�������g���J�E���g���邩��1
					for(int k=0;k<PART[i].N;k++)
					{       
						int j=PART[i].NEI[k];
						if(PART[j].surface==ON && (PART[j].type==FLUID || PART[j].type==ELASTIC))
						{
							num++;
							for(int D=0;D<3;D++) newB[D][i]+=B[D][j];
						}
					}
					for(int D=0;D<3;D++) newB[D][i]/=num;
				}
			} 
			for(int i=1;i<=nelm;i++) for(int D=0;D<3;D++) B[D][i]=newB[D][i];
		}
	}
    for(int D=0;D<3;D++) delete [] newB[D];
*/
	ofstream fp("./FluxAVS/Flux.dat");
	double xmax=-100;						//�o�͗��q�̍ő剡���W
	double ymax=-100;						//�o�͗��q�̍ő�c���W
	double times=1.0; //=CON.get_times()/CON.get_particle_mass();*le*le*CON.get_FEMtimes(); //?
	
	double **center_of_element=new double* [nelm+1];
	for(int i=1;i<=nelm;i++) center_of_element[i]=new double[3];

	//�d�S�̌v�Z�E�E�Enelm�S�ďo��
	for(int i=1;i<=nelm;i++)
	{
		for(int D=0;D<3;D++)
		{
			center_of_element[i][D]=0.0;
			for(int j=1;j<=4;j++) center_of_element[i][D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
		}
	}

	for(int i=1;i<=nelm;i++)
    {
		if(center_of_element[i][A_Y]>-le*0.5 && center_of_element[i][A_Y]<+le*0.5)
		{
			fp<<center_of_element[i][A_X]<<"\t"<<center_of_element[i][A_Z]<<"\t"<<B[A_X][i]*times<<"\t"<<B[A_Z][i]*times<<endl;
			if(center_of_element[i][A_X]>xmax) xmax=center_of_element[i][A_X];
			if(center_of_element[i][A_Z]>ymax) ymax=center_of_element[i][A_Z];
		}
    }
	xmax+=4*le;//�}����o���ʒu��ی��������ď����΂߂Ɉړ�
	ymax+=4*le;
	if(CON.get_B_times()>0) fp<<xmax<<" "<<ymax<<" "<<CON.get_B_times()*times<<" "<<0*times<<endl;//�Ō�ɖ}��o��
	fp.close();////

	//////////////////////////////////////
//	if(t==1 || t%(CON.get_EM_interval()*2)==0)
//	{
	int timestep=CON.get_current_step();
	int half_nelm=0;
	for(int i=1;i<=nelm;i++){				//X-Z�ʂɕ\������f�ʂ܂ł̗v�f��
		if(center_of_element[i][A_Y]<0.0) half_nelm++;
	}
	///////////////////////////////////////�������x�x�N�g���t�@�C��///////////////////////////////////////////////
		stringstream ssfd;
		ssfd<<"./FluxAVS/FluxDensity"<<t<<".fld";
		string FluxDensity=ssfd.str();	

		ofstream ffd(FluxDensity);
		if(ffd.fail()){
			system("mkdir FluxAVS");
		ofstream ffd(FluxDensity);
		if(ffd.fail()){
			cout<<"�t�@�C�����J���܂���ł����BFluxAVS�t�H���_�����邩�ǂ����m�F���Ă�������"<<endl;
			exit(1);
			}
		}

		ffd << "# AVS field file" << endl;
		ffd << "ndim=1" << endl;
		ffd << "dim1=" << half_nelm <<endl;
		ffd << "nspace=3" << endl;
		ffd << "veclen=3" << endl;
		ffd << "data=float" << endl;
		ffd << "field=irregular" << endl;
		ffd << "label=e-x e-y e-z" << endl << endl;
		ffd << "variable 1 file=./FluxDensity"<<t<<" filetype=ascii skip=1 offset=0 stride=6" << endl;
		ffd << "variable 2 file=./FluxDensity"<<t<<" filetype=ascii skip=1 offset=1 stride=6" << endl;
		ffd << "variable 3 file=./FluxDensity"<<t<<" filetype=ascii skip=1 offset=2 stride=6" << endl;
		ffd << "coord    1 file=./FluxDensity"<<t<<" filetype=ascii skip=1 offset=3 stride=6" << endl;
		ffd << "coord    2 file=./FluxDensity"<<t<<" filetype=ascii skip=1 offset=4 stride=6" << endl;
		ffd << "coord    3 file=./FluxDensity"<<t<<" filetype=ascii skip=1 offset=5 stride=6" << endl;

		ffd.close();

		stringstream ssfdd;
		ssfdd<<"./FluxAVS/FluxDensity"<<t;
		string FluxDensitydata=ssfdd.str();

		ofstream fout(FluxDensitydata);
		if(fout.fail()){
			cout<<"�f�[�^�t�@�C�����J���܂���ł����B"<<endl;
			exit(1);
		}

		fout<<"e-x e-y e-z x y z"<<endl;
		for(int i=1;i<=nelm;i++)
		{
			if(center_of_element[i][A_Y]<0){
				fout<<B[A_X][i]<<" "<<B[A_Y][i]<<" "<<B[A_Z][i]<<" "<<center_of_element[i][A_X]<<" "<<center_of_element[i][A_Y]<<" "<<center_of_element[i][A_Z]<<endl;
			}
		}
		fout.close();
		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////�������x�x�N�g���f�ʐ}///////////////////////////////////////
		//�ǂݍ��ݗp�f�[�^�t�@�C���̍쐬
		stringstream sscfd;
		sscfd<<"FluxAVS/C_FluxDensity"<<t;
		string CutFluxDensity=sscfd.str();

		ofstream flux_numbered(CutFluxDensity);
		if(flux_numbered.fail()){
			cout<<"�t�@�C�����J���܂���ł����BFluxAVS�t�H���_�����邩�ǂ����m�F���Ă�������2"<<endl;
			exit(1);
		}

		flux_numbered<<"e-x e-y e-z x y z"<<endl;

		int counter=0;//�f�[�^�t�@�C���̃f�[�^�_��
		for(int i=1;i<=nelm;i++)
		{	
			if(center_of_element[i][A_Y]>-le*0.5 && center_of_element[i][A_Y]<le*0.5){
				flux_numbered<<B[A_X][i]<<" "<<B[A_Y][i]<<" "<<B[A_Z][i]<<" "<<center_of_element[i][A_X]<<" "<<center_of_element[i][A_Y]<<" "<<center_of_element[i][A_Z]<<endl;
				counter++;
			}
		}
		flux_numbered.close();

		//�\���p�\���i�q�^�f�[�^�t�@�C���̍쐬
		sscfd.clear();
		string filename=CutFluxDensity+".fld";
		ofstream flux(filename);

		flux << "# AVS field file" << endl;
		flux << "ndim=1" << endl;
		flux << "dim1=" << counter <<endl;
		flux << "nspace=3" << endl;
		flux << "veclen=3" << endl;
		flux << "data=float" << endl;
		flux << "field=irregular" << endl;
		flux << "label=e-x e-y e-z" << endl << endl;
		flux << "variable 1 file=./C_FluxDensity"<<t<<" filetype=ascii skip=1 offset=0 stride=6" << endl;
		flux << "variable 2 file=./C_FluxDensity"<<t<<" filetype=ascii skip=1 offset=1 stride=6" << endl;
		flux << "variable 3 file=./C_FluxDensity"<<t<<" filetype=ascii skip=1 offset=2 stride=6" << endl;
		flux << "coord    1 file=./C_FluxDensity"<<t<<" filetype=ascii skip=1 offset=3 stride=6" << endl;
		flux << "coord	  2 file=./C_FluxDensity"<<t<<" filetype=ascii skip=1 offset=4 stride=6" << endl;
		flux << "coord    3 file=./C_FluxDensity"<<t<<" filetype=ascii skip=1 offset=5 stride=6" << endl;

		flux.close();
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
/*	//�\���i�q�f�[�^�t�@�C��
	///////////////////////////////////�������x�R���^�[�}///////////////////////////////////////////////////////
	stringstream sstr;
	sstr<<"./FluxContour/FluxContour"<<t<<".fld";
	string FluxContour=sstr.str();
	
	ofstream ffc(FluxContour);
	if(ffc.fail()){
		system("mkdir FluxContour");
		ofstream ffc(FluxContour);
		if(ffc.fail()){
		cout<<"./FluxContour�t�H���_���J���܂���ł���"<<endl;
		exit(1);
		}
	}
	
	ffc << "# AVS field file" << endl;
	ffc << "ndim=1" << endl;
	ffc << "dim1=" << half_nelm <<endl;
	ffc << "nspace=3" << endl;
	ffc << "veclen=1" << endl;
	ffc << "data=float" << endl;
	ffc << "field=irregular" << endl;
	ffc << "label=FluxContour" << endl << endl; //e-x e-y e-z
	ffc << "variable 1 file=./FluxContour"<<t<<" filetype=ascii skip=1 offset=0 stride=4" << endl;
	ffc << "coord    1 file=./FluxContour"<<t<<" filetype=ascii skip=1 offset=1 stride=4" << endl;
	ffc << "coord    2 file=./FluxContour"<<t<<" filetype=ascii skip=1 offset=2 stride=4" << endl;
	ffc << "coord    3 file=./FluxContour"<<t<<" filetype=ascii skip=1 offset=3 stride=4" << endl;
	ffc.close();
	//�f�[�^�t�@�C��
	stringstream ssfc;
	ssfc<<"./FluxContour/FluxContour"<<t;
	string datafluxcontour=ssfc.str();
	ofstream fout3(datafluxcontour);
	if(fout3.fail()){
		cout<<"./FluxContour�t�H���_���J���܂���ł���"<<endl;
		exit(1);
	}
	fout3<<"e-x e-y e-z x y z"<<endl;
	for(int i=1;i<=nelm;i++)
    {
		if(center_of_element[i][A_X]<0){
//		fout<<F[A_X][i]*times<<" "<<F[A_Y][i]*times<<" "<<F[A_Z][i]*times<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
		fout3<<sqrt(B[A_X][i]*B[A_X][i]+B[A_Y][i]*B[A_Y][i]+B[A_Z][i]*B[A_Z][i])<<" "<<center_of_element[i][A_X]<<" "<<center_of_element[i][A_Y]<<" "<<center_of_element[i][A_Z]<<endl;
		}
	}
	fout3.close();*/
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	}
	for(int i=1;i<=nelm;i++) delete [] center_of_element[i];
	delete [] center_of_element;
}

//TetGen�p �O�����L���v�f�@
void usingTetGen(mpsconfig &CON,vector<mpselastic> &PART, double **F, int fluid_number,int particle_number,int t,double TIME)
{
	double dt=CON.get_dt();
//	double TIME=t*dt;
	//////////////////////////////////////////////////////////////////////////////////////////
	//TetGen�p�z��쐬
	vector<tetgen_node> NODEall;
	vector<tetgen_facet> FACEall;
	vector<tetgen_element> ELEMall;

	//TRANS[i]�ɂ́A�ߓ_�ԍ�i�ɑΉ����闱�q�ԍ����i�[����B
	//FEM3D.cpp �ł́A�ߓ_�ԍ���1����n�܂�̂ŁATRANS[0]�ɂ͐錾���-1������B
	vector<int> TRANS;
	TRANS.push_back(-1);	//TRANS[0]�ɂ�-1

    ///TetGen�ɂ�郁�b�V������////////////////////////////////////////////////////

	tetgen_function TETFUNC;
    TETFUNC.call_TetGen(CON,PART,fluid_number,particle_number,NODEall,FACEall,ELEMall,TRANS);
	
	///////////////////////////////////////////////////////////////////////////////

	int N=static_cast<int>(TRANS.size())-1;	//FEM�ߓ_�Ɋ܂܂�闱�q��(0�Ԗڂ�����TRANS[]�̒���)
	int node=static_cast<int>(NODEall.size()); //�ߓ_��
	int nelm=static_cast<int>(ELEMall.size()); //�v�f��
	int KTJ=node; //�ő�ߓ_��
	int KTE=nelm*6; //�ő�v�f��
	
	int *depth=new int [KTE+1];			//�e�v�f�̐[���i�[
//	vector<int> depth(KTE+1);

	//�z��m��
	vector<point3D> NODE(NODEall.size()+1);
//	NODE.reserve(NODEall.size()+1);
	vector<element3D> ELEM(ELEMall.size()+1);
//	ELEM.reserve(ELEMall.size()+1);
	
	cout<<"NODE.size(): "<<NODE.size()<<", ELEM.size(): "<<ELEM.size()<<endl;
	//TetGen�̐ߓ_�E�v�f�f�[�^���擾	�ߓ_�ԍ���1���炷
	//�ߓ_�f�[�^
//	int num_magnet_attr=0;
    for(int i=0;i<node;i++)
    {
		NODE[i+1].r[A_X]=NODEall[i].r[A_X];
		NODE[i+1].r[A_Y]=NODEall[i].r[A_Y];
		NODE[i+1].r[A_Z]=NODEall[i].r[A_Z];
		NODE[i+1].material=static_cast<int>(NODEall[i].attribute);
		NODE[i+1].boundary_condition=static_cast<int>(NODEall[i].boundary);
		NODE[i].particleID=-1; //-1�ŏ�����
//		if(static_cast<int>(NODEall[i].attribute==MAGNET)) num_magnet_attr++;
    }

	////////////���q�ƑΉ�����ړ_�͂����ŗ��q�ԍ�������///////////////
	for(int i=0;i<TRANS.size();i++){
		NODE[i].particleID=TRANS[i];
	}
	//////////////////////////////////////////////////////////////////////

//	cout<<"num_magnet_attr: "<<num_magnet_attr<<endl;
	
	//�v�f�f�[�^
    for(int i=0;i<nelm;i++)
	{
		//�ގ�
		ELEM[i+1].material=ELEMall[i].attribute;

		//�\���ߓ_
		for(int n=0;n<4;n++) ELEM[i+1].node[n+1]=ELEMall[i].node[n]+1;
		
		//�v�f-�v�f�֌W
		for(int n=0;n<4;n++) 
		{
			if(ELEMall[i].nei_elem[n]==-1){
				ELEM[i+1].elm[n+1]=0;	//�ߗחv�f�Ȃ��̏ꍇTetGen��-1��Ԃ��Ă��邽��0�ɏC������
			}else{
				ELEM[i+1].elm[n+1]=ELEMall[i].nei_elem[n]+1;
			}
		}
	}


	//���b�V�������̃|�X�g����
		double *val=new double[KTJ+1];
		for(int i=1;i<=node;i++) val[i]=NODE[i].material;
//		if(t==1 || t%(CON.get_EM_interval()*CON.get_mesh_output_interval())==0) 
		{
//			data_avs(node,nelm,NODE,ELEM,KTJ,val,CON);
			data_avs2(CON,node,nelm,NODE,ELEM,KTJ,t);//�f�ʐ}
			data_avs3(node,nelm,NODE,ELEM,CON);//�ގ�
		}

		delete [] val;
	//double min_volume=10;
/*	//�̐ς���ъO�ڋ��p�����[�^�v�Z
    for(int i=1;i<=nelm;i++)
	{
		//4�̐ߓ_
		int ia=ELEM[i].node[1];
		int ib=ELEM[i].node[2];
		int ic=ELEM[i].node[3];
		int ip=ELEM[i].node[4];
		//�v�f�̐όv�Z
		ELEM[i].volume=volume3D(NODE,ia,ib,ic,ip);
		
		//�O�ڋ����S���W����є��a�v�Z
		sphere3D(NODE,ELEM,ia,ib,ic,ip,i);
	}
	//cout<<"�ŏ��̐�="<<min_volume<<endl;*/
	////////////////////////////////////////////////////////////////
	//���b�V�����؂�Ă���͂̌v�Z�܂�

	//�ߓ_-�v�f�֌W
    int *jnb=new int[node+1];//�e�ߓ_�ɗאڂ���v�f���i�[
    set_jnb3D(NODE,ELEM,node,nelm,jnb);

	int **nei=new int* [node+1];//�e�ߓ_�̎��ӗv�f�ԍ��i�[
    for(int i=1;i<=node;i++) nei[i]=new int [jnb[i]+1];
    set_nei3D(NODE,ELEM,node,nelm,jnb,nei);

	for(int i=1;i<=node;i++)
	{
		if(jnb[i]==0)
		{
			//cout<<"jnb=0 i="<<i<<" material="<<NODE[i].material<<" particle="<<NODE[i].particleID<<endl;
			NODE[i].boundary_condition=1;//���E�������f�B���N���^�ɂ��邱�ƂŁAICCG�ɎQ�������Ȃ�
			//if(NODE[i].material==FLUID) NODE[i].particleID=-1;//���̐ߓ_����������ꍇ�A�K�v�Ȓl���v�Z����Ȃ����߁A���q�Ƀt�B�[�h�o�b�N�ł��Ȃ�(���Ă͂����Ȃ�)
			if(NODE[i].material==FLUID || NODE[i].material==ELASTIC) if(NODE[i].particleID>=0) cout<<"suf="<<PART[NODE[i].particleID].surface<<endl;
		}
	}

	//FEM
	if(CON.get_EM_calc_type()==1 || CON.get_EM_calc_type()==4) potential_calculation(CON,NODE,ELEM, node, nelm,jnb, TIME,PART, fluid_number,nei,F);
	if(CON.get_EM_calc_type()==2) calc_static_magnetic_field(CON, node, nelm,NODE,ELEM,jnb, dt,TIME, t,nei, KTE,PART,fluid_number,F,KTJ);
	if(CON.get_EM_calc_type()==3) calc_transitional_EM_field(CON, node, nelm,NODE,ELEM,jnb, dt, TIME,t,nei, KTE,PART,fluid_number,F,KTJ);
	if(CON.get_EM_calc_type()==5) calc_variable_magnetic_field(CON, node, nelm,NODE,ELEM,jnb, dt,TIME, t,nei, KTE,PART,fluid_number,F);
	//FEM_interval�ɂ��v�Z���ԒZ�k�̏ꍇ
    if(CON.get_EM_interval()>1 && CON.get_Ferror()==OFF) //����[���Ȃ��Ȃ�d���͂��t�@�C���ɏ�������
    {
		//ofstream bb("FEM_interval.dat");///�d���͏o�́@���X�e�b�v�͂����ǂݎ��
        FILE *b=fopen("FEM_interval.dat","w");///�d���͏o�́@���X�e�b�v�͂����ǂݎ��
		for(int i=0;i<fluid_number;i++)
		{
			//bb<<F[A_X][i]<<" "<<F[A_Y][i]<<" "<<F[A_Z][i]<<endl;
			fprintf(b,"%1.15lf %1.15lf %1.15lf\n",F[A_X][i],F[A_Y][i],F[A_Z][i]);
		}
		fclose(b);
		//bb.close();
    }

	delete [] jnb;
    for(int i=1;i<=node;i++) delete [] nei[i];
    delete [] nei;

/*	///////////////////////////////////////////////////////////////
    //�ߓ_-�v�f�֌W
	int *jnb=new int[node+1]; //�e�ߓ_�ɗאڂ���v�f���i�[
	set_jnb3D(NODE, ELEM, node, nelm, jnb); //delaun3D��set_jnb3D�̃o�O�΍􂪏����Ă���
	int **nei=new int*[node+1]; //�e�ߓ_�̎��ӗv�f�ԍ��i�[
	for(int i=1;i<=node;i++) nei[i]=new int[jnb[i]+1]; 
    set_nei3D(NODE, ELEM, node, nelm, jnb, nei);
	cout<<"�v�f��="<<nelm<<" �ߓ_��="<<node<<endl;
//���b�V����������
	//�\�ʃ��b�V���o��
	if(t==1 || t%CON.get_avs_mesh3_interval()==0)
	{
		t_data_avs3(node,nelm,NODE,ELEM,CON,t);//�ގ�
	}
    //�����܂ł�OK, 2012-10-30, 14:14
	//���f�����̐ݒ�(���E�����Ȃ�)
	//�Ód����
	if(CON.get_model_number()==14)
	{
	    for(int i=0;i<node;i++)
		{
			//i��1����Ă��邱�Ƃɒ���
			if(NODEall[i].boundary==ELECTRODE1)	NODE[i+1].boundary_condition=1;	//�~���d�ɂ���ѓy��
			else if(NODEall[i].boundary==ELECTRODE2) NODE[i+1].boundary_condition=2;	//���d��
			else NODE[i+1].boundary_condition=0;	//���̑��̐ߓ_�͖��m��
		}
		//�v�f-�v�f�֌W�ɂ���U����
	//	if(CON.get_fly_judge()==2)	fly_judge_FEMelement(CON,PART,NODE,ELEM,node,nelm,TRANS);
	}
	//�����G���X�g�}�[�̋��E�����w��E�E�E���ۂ̋��E�����́H
	//AIR�́H
	if(CON.get_model_number()==5)
	{
	    for(int i=0;i<node;i++)
		{
			//i��1����Ă��邱�Ƃɒ���			
			if(NODEall[i].boundary==MAGNET){
				NODE[i+1].boundary_condition=1;//1;	//����
				NODE[i+1].particleID=-1;
			}else{
				NODE[i+1].boundary_condition=0;	//���̑��̐ߓ_�͖��m���E�E�EAIR��
			}
		}
	}
	//set_material()��delaun3D.cpp�@���ꂪ�Ȃ��Ɨv�f�ގ���`�ł��Ȃ��H�H�H
//�L���v�f�@�v�Z�J�n
	//�z��m��
	//����v�Z�̏���
	double *V=new double[node+1];
	//vector<double> V(node+1);	//�d��
	//vector<vector<double> > Eb(3,vector<double>(nelm+1));//�v�f���d�E or �v�f�����E	
    double *RP=new double[nelm+1];//�e�v�f�̓������܂��͗U�d���i�[�i���݂ł͗U�d���͎g���Ă��Ȃ��j
	
	//FEM�����s�����Ƃ��ł����b�V�����݂�邩��A�ŏ����������ŏo��
	if(t==1){
		//data_avs(node,nelm,NODE,ELEM,KTJ,V,CON);
		//data_avs2(CON,node,nelm,NODE,ELEM,KTJ,t);
	}
	//������
	if(CON.get_FEM_calc_type()==2 || CON.get_FEM_calc_type()==3 || CON.get_FEM_calc_type()==5)
	{
		sides* SIDE=new sides[KTE]; //�ӃN���X[KTE+1]�ł͂Ȃ�
		int* branch_num=new int[node+1]; //�e�ߓ_���אڂ���ߓ_��(�e�ߓ_���牄�т�ӂ̐�)		
		int max=1000; //�ߓ_�ɗאڂ���ő�ߓ_��		
		int** nei2=new int*[node+1];
		for(int i=1;i<=node;i++) nei2[i]=new int[max];//�C�R�[���ɒ���
		int side_num=0;	//�S�Ӑ��i�[
		//�ӗv�f���� (�ߓ_�v�f���g�p����ꍇ�ł��A�d�����x�����߂�Ƃ��ɕӗv�f���ق���)
		side_num=make_edge_element(CON,NODE,ELEM,node,nelm,jnb,nei,SIDE,branch_num,nei2,KTE);
		cout<<"�ӗv�f�����I��"<<endl;
		//����v�Z�̏���
		double *B[3]; //�v�f�����E
		for(int D=0;D<3;D++) B[D]=new double [nelm+1];
		double** current=new double*[DIMENSION];
		for(int D=0;D<3;D++) current[D]=new double[nelm+1];
		if(CON.get_J_input_way()==0){
			for(int D=0;D<3;D++){
				for(int i=0;i<=nelm;i++){
					current[D][i]=0.0;//�����d�������݂��Ȃ����珉����
				}
			}
		}
		for(int n=1;n<=nelm;n++)
		{
			RP[n]=1;
			if(ELEM[n].material==ELASTIC) RP[n]=CON.get_RP();
		}
		if(CON.get_ele_type()==0){
			cout<<"�ߓ_�v�f�ɂ�鎥��v�Z�͖���`�I"<<endl;
			exit(EXIT_FAILURE);
		}
		if(CON.get_ele_type()==1)//�ӗv�f
		{
			double *A=new double[side_num+1];//�x�N�g���|�e���V����
		
			//�x�N�g���|�e���V��������
			if(CON.get_FEM_calc_type()==2) Avector3D(CON,NODE,ELEM,SIDE,node,nelm,side_num,A,jnb,branch_num,current,RP);
			if(CON.get_FEM_calc_type()==2) Avector3D_OK(CON, NODE, ELEM, SIDE, node, nelm, side_num, A, jnb, branch_num, depth);
			//else if(CON.get_FEM_calc_type()==5) non_linear_Avector3D(CON,NODE,ELEM,SIDE,node,nelm,side_num,A,jnb,branch_num,RP,Eb);//����`
			
			//�������x����
			Bflux3D_side(CON,NODE,ELEM,SIDE,node,nelm,side_num,A,B,RP);
	
			//�d���͉���
			if(CON.get_m_force()==0) NODE_F3D(CON,NODE,ELEM,node,nelm,B,jnb,nei,RP,PART,F,fluid_number);//�ߓ_�͖@
//			else if(CON.get_m_force()==1) VolumeForce3D(&CON,PART,NODE,ELEM,nelm,B,fluid_number,RP,jnb,nei,N,TRANS,particle_number);
//			else if(CON.get_m_force()==2) kelvin_force3D(&CON,PART,NODE,ELEM,nelm,B,fluid_number,RP,N,TRANS,jnb,nei);
//			else if(CON.get_m_force()==3) integral_surface_F3D(&CON,NODE,ELEM,node,nelm,B,jnb,nei,RP,PART,N,TRANS,fluid_number,depth);
//			else if(CON.get_m_force()==4) direct_divT3D(&CON,PART,NODE,ELEM,nelm,B,fluid_number,RP,jnb,nei,N,TRANS);
//			else if(CON.get_m_force()==5) virtual_air_gap_F3D(&CON,NODE,ELEM,node,nelm,B,jnb,nei,RP,PART,N,TRANS,fluid_number,depth);
//			else if(CON.get_m_force()==6) NODE_F3D_with_integral(&CON,NODE,ELEM,node,nelm,B,jnb,nei,RP,PART,N,TRANS,fluid_number);//�ߓ_�͖@
//			else if(CON.get_m_force()==7) Magnetic_Charge_Method_F3D(&CON,NODE,ELEM,node,nelm,B,jnb,nei,RP,PART,N,TRANS,fluid_number);
			delete [] A;
		}
		delete [] SIDE;
		delete [] branch_num;
		for(int i=0;i<=node;i++) delete [] nei2[i];
		delete [] nei2;
		for(int D=0;D<3;D++) delete [] current[D];
		delete [] current;
		for(int D=0;D<3;D++) delete [] B[D];
	}
	//�d���̓X���[�W���O
	smoothingF3D(CON,PART,fluid_number,F,t);
	//�f�ʃ��b�V���o��
	if(t==1 || t%CON.get_avs_mesh2_interval()==0)
	{
		if(CON.get_FEM_calc_type()==1) data_avs(node,nelm,NODE,ELEM,KTJ,V,CON);
		//data_avs2(CON,node,nelm,NODE,ELEM,KTJ,V,t);
		data_avs2(CON,node,nelm,NODE,ELEM,KTJ,t); //������V�Ȃǂ��w�肷��Βl���m�F�ł���
		//data_avs2_movie(CON,node,nelm,NODE,ELEM,KTJ,V,t);
	}
	//�\�ʃ��b�V���o��
	if(t==1 || t%CON.get_avs_mesh3_interval()==0)
	{
		if(CON.get_FEM_calc_type()==1) data_avs(node,nelm,NODE,ELEM,KTJ,V,CON);
//		data_avs3(node,nelm,NODE,ELEM,CON,t);//�ގ�
	}
 
	delete [] V;
	delete [] RP;
    delete [] jnb;
	for(int i=1;i<=node;i++) delete [] nei[i];
    delete [] nei;*/
	delete [] depth;
	cout<<"ok2"<<endl;
}

//�d���̓v���b�g�֐�
void plot_F(mpsconfig &CON, vector<mpselastic> &PART, int fluid_number, int t)
{
	double le=CON.get_distancebp();
	double xmax=-100;						//�o�͗��q�̍ő剡���W
	double ymax=-100;						//�o�͗��q�̍ő�c���W
	double times=CON.get_times()/CON.get_density()/le*CON.get_FEMtimes();
	double times_Pa=CON.get_times_Pa()*CON.get_FEMtimes();

	double Fz_sum=0;						//�d���͂�Z���������̍��v
	double Fall_sum=0;						//�d���̓x�N�g���̒����̍��v
	
	ofstream fp("F.dat");

	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].r[A_Y]>-le*0.5 && PART[i].r[A_Y]<le*0.5)
		{
			if(CON.get_plot_F_type()==0) fp<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<PART[i].eforce[A_X]*times<<"\t"<<PART[i].eforce[A_Z]*times<<endl;//[N]�\��
//			if(CON.get_plot_F_type()==1) fp<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<PART[i].dirP_em*PART[i].n[A_X]*times_Pa<<"\t"<<PART[i].dirP_em*PART[i].n[A_Z]*times_Pa<<endl;//[Pa]�\��
			if(PART[i].r[A_X]>xmax) xmax=PART[i].r[A_X];
			if(PART[i].r[A_Z]>ymax) ymax=PART[i].r[A_Z];
		}
//		if(PART[i].fly==TOUCH)	Fz_sum+=PART[i].eforce[A_Z];	//Z���������𑫂��Ă���
//		if(PART[i].type==BOFLUID)	Fall_sum+=sqrt(PART[i].eforce[A_X]*PART[i].eforce[A_X]+PART[i].eforce[A_Y]*PART[i].eforce[A_Y]+PART[i].eforce[A_Z]*PART[i].eforce[A_Z]);
	}
	
	//�}��   �o�͈ʒu��ی��������ď����΂߂Ɉړ�
	xmax+=4*le;
	ymax+=4*le;
	if(CON.get_plot_F_type()==0)	if(CON.get_legend_F()>0) fp<<xmax<<" "<<ymax<<" "<<CON.get_legend_F()*times<<" "<<0*times<<endl;//�Ō�ɖ}��o�� [N]�\��
	if(CON.get_plot_F_type()==1)	if(CON.get_legend_F_Pa()>0) fp<<xmax<<" "<<ymax<<" "<<CON.get_legend_F_Pa()*times_Pa<<" "<<0*times_Pa<<endl;//�Ō�ɖ}��o�� [Pa]�\��
	fp.close();/////


	//Z�����̍��͂��v���b�g
	if(t==1)//1�X�e�b�v�ڂ̓t�@�C�������Z�b�g
	{
		ofstream freset("Fz.dat");
		freset.close();
	}
	ofstream fout("Fz.dat",ios::app);
	fout<<t*CON.get_dt()<<setprecision(SIGNIFY)<<scientific<<"\t"<<Fz_sum<<endl;
	fout.close();

	//Z�����̍��͂��v���b�g
	if(t==1)//1�X�e�b�v�ڂ̓t�@�C�������Z�b�g
	{
		ofstream allreset("Fall.dat");
		allreset.close();
	}
	ofstream all("Fall.dat",ios::app);
	all<<t*CON.get_dt()<<setprecision(SIGNIFY)<<scientific<<"\t"<<Fall_sum<<endl;
	all.close();
	//�m�F
	cout<<"�Ód�͂̑傫���̑��a Fall="<<Fall_sum<<endl;

}

//�d���̓v���b�g�֐�(���O)
void plot_F_log(mpsconfig &CON, vector<mpselastic> &PART, int fluid_number,int t)
{
	double xmax=-100;						//�o�͗��q�̍ő剡���W
	double ymax=-100;						//�o�͗��q�̍ő�c���W
	double le=CON.get_distancebp();
	double times=CON.get_times()/CON.get_density()/le*CON.get_FEMtimes();
	double times_Pa=CON.get_times_Pa();
	
	char filename[20];
	sprintf_s(filename,"F%d.dat", t);
	ofstream fp(filename);

	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].r[A_Y]>-le*0.5 && PART[i].r[A_Y]<le*0.5)
		{
			if(CON.get_plot_F_type()==0)	fp<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<PART[i].eforce[A_X]*times<<"\t"<<PART[i].eforce[A_Z]*times<<endl;//[N]�\��
			if(CON.get_plot_F_type()==1)	fp<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<PART[i].dirP_em*PART[i].n[A_X]*times_Pa<<"\t"<<PART[i].dirP_em*PART[i].n[A_Z]*times_Pa<<endl;//[Pa]�\��
			if(PART[i].r[A_X]>xmax) xmax=PART[i].r[A_X];
			if(PART[i].r[A_Z]>ymax) ymax=PART[i].r[A_Z];
		}
    }
	
	//�}��   �o�͈ʒu��ی��������ď����΂߂Ɉړ�
	xmax+=4*le;
	ymax+=4*le;
	if(CON.get_plot_F_type()==0)	if(CON.get_legend_F()>0) fp<<xmax<<" "<<ymax<<" "<<CON.get_legend_F()*times<<" "<<0*times<<endl;//�Ō�ɖ}��o�� [N]�\��
	if(CON.get_plot_F_type()==1)	if(CON.get_legend_F_Pa()>0) fp<<xmax<<" "<<ymax<<" "<<CON.get_legend_F_Pa()*times_Pa<<" "<<0*times_Pa<<endl;//�Ō�ɖ}��o�� [Pa]�\��
	fp.close();/////
}


//OpenMP�ɂ��ICCG�@
void parallel_ICCG3D2(mpsconfig &CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X)
{
	//val :�[���v�f�̒l
	//ind:��[���v�f�̗�ԍ��i�[�z��
	//ptr:�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
	//X[n]:��

	double accel=CON.get_CGaccl();//�����t�@�N�^
	
	int num2=0;//�Ίp�������܂ށA���O�p�s�񂾂����l���ɂ��ꂽ��[���v�f��
	for(int k=0;k<pn;k++) for(int m=ptr[k];m<ptr[k+1];m++) if(ind[m]<=k) num2++;	
	if(num2!=(number-pn)/2+pn) cout<<"ERROR"<<endl;
	
	double *val2=new double [num2];
	int *ind2 = new int [num2];
	int *ptr2 = new int [pn+1];

	num2=0;
	for(int k=0;k<pn;k++)
	{	
		ptr2[k]=num2;
	    for(int m=ptr[k];m<ptr[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
			if(ind[m]<=k)
			{
				val2[num2]=val[m];
				ind2[num2]=ind[m];
				if(ind[m]==k) val2[num2]*=accel;//����̧��
				num2++;
			}
		}
	}
	ptr2[pn]=num2;//��������Ă����Ȃ��ƁA�Ō��(int m=ptr2[k];m<ptr2[k+1];m++)�݂����Ȃ��Ƃ��ł��Ȃ�

	int *NUM = new int [pn];
	for(int k=0;k<pn;k++) NUM[k]=0;
	
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
			int J=ind2[m];
			NUM[J]=NUM[J]+1;
		}
	}
	double **VAL=new double *[pn];//�[���v�f�̒l
	for(int i=0;i<pn;i++) VAL[i]=new double [NUM[i]];
	int **IND = new int *[pn];//��[���v�f�̍s�ԍ��i�[�z��
	for(int i=0;i<pn;i++) IND[i]=new int [NUM[i]];
	
	/////////////////////////////////////iccg�@
	double alp,beta;
	double rLDLt_r;
	double E=1;//�덷
    double *r=new double[pn];
	for(int n=0;n<pn;n++) X[n]=0;
	
	double *AP = new double [pn];
	double *P = new double [pn];
	double *y=new double [pn];
	double *LDLt_r= new double [pn];
	double *D1 = new double [pn];//D�s��
	
	/////�s���S�R���X�L�����
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
	        int i=ind2[m];//��ԍ�
	        if(i==0)
			{
				val2[m]=val2[m];
				if(val2[m]<0.0001 &&val2[m]>=0) val2[m]=0.0001;
				else if(val2[m]>-0.0001 &&val2[m]<=0) val2[m]=-0.0001;
		    
			}
			if(i>0 && i<k)
			{
				double sum=0;
				
				for(int j=ptr2[k];j<m;j++)
				{	
					for(int J=ptr2[i];J<ptr2[i+1];J++)//i�s�ڂ̂Ȃ������̈�v������̂�T���Ă���B������Ԃ��H
					{
						if(ind2[J]==ind2[j]) sum+=val2[j]*val2[J]*D1[ind2[j]];
					}
				}
				val2[m]=val2[m]-sum;
			}
			if(i==k)
			{
				double sum=0;
				for(int j=ptr2[k];j<m;j++) sum+=val2[j]*val2[j]*D1[ind2[j]];
				val2[m]=val2[m]-sum;
				if(val2[m]<0.0001 &&val2[m]>=0) val2[m]=0.0001;
				else if(val2[m]>-0.0001 &&val2[m]<=0) val2[m]=-0.0001;
				D1[k]=1/val2[m];
            }
	    }
	}    
	///�s���S�R���X�L�[��������/////////*/

	///�����ɂ����z��ɒl����
	for(int k=0;k<pn;k++) NUM[k]=0;
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
			int J=ind2[m];
			VAL[J][NUM[J]]=val2[m];
			IND[J][NUM[J]]=k;
			NUM[J]=NUM[J]+1;
		}
	}////////*/

	for(int n=0;n<pn;n++) //�����l
	{
		X[n]=0;
		r[n]=B[n];
	}

	/////////////////y[i]�����Ƃ߁ALDLt_r[i]�����Ƃ߂�B
	for(int i=0;i<pn;i++)
	{
		if(i==0) y[0]=r[0]/val2[0]; //���i3.77�j 
		else
		{
		    double sum=0;
		    /////////        
		    for(int m=ptr2[i];m<ptr2[i+1]-1;m++) sum+=val2[m]*y[ind2[m]];//���i3.78�j
		    int m=ptr2[i+1]-1;
		    y[i]=(r[i]-sum)/val2[m];
		}
	}////y[i]�����Ƃ܂����B
	for(int i=pn-1;i>=0;i--)
	{
	    double sum=0;
		for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
	    LDLt_r[i]=y[i]-D1[i]*sum;	
	}
	/////////////////*/
	
	for(int n=0;n<pn;n++) P[n]=LDLt_r[n];

    cout<<"parallel_ICCG�@�X�^�[�g  -----���m��="<<pn<<"  ---"<<endl;
	ofstream Eee("PICCG.dat", ios::trunc);
	Eee.close();
	unsigned int time=GetTickCount();
	int count=0;
	double lowE=100;
	double ep=CON.get_FEMCGep();//��������
	while(E>ep)
	{
//		if(E<=lowE) lowE=E;
		count++;
		if(count==pn){
			cout<<"count>pn E="<<E<<"lowE="<<lowE<<endl;
//			ep=lowE;
		}
		//////////////alp�����߂�
		rLDLt_r=0;
		double PAP=0;
		//for(int n=0;n<pn;n++) rLDLt_r+=r[n]*LDLt_r[n];//sirial
		#pragma omp parallel for reduction(+:rLDLt_r) reduction(+:PAP)
		for(int n=0;n<pn;n++)
		{	//printf("%d\n",omp_get_thread_num());
			rLDLt_r+=r[n]*LDLt_r[n];
			AP[n]=0;
			for(int m=ptr[n];m<ptr[n+1];m++) AP[n]+=val[m]*P[ind[m]];
			PAP+=P[n]*AP[n];
		}
		
		//for(int n=0;n<pn;n++)  PAP+=P[n]*AP[n];//sirial
		alp=rLDLt_r/PAP;
		
//		cout<<"alp="<<alp<<endl;
		//////////////////////
	
		//////////////// X(k+1)=X(k)+alp*P
		E=0;
		#pragma omp parallel for reduction(+:E)
		for(int n=0;n<pn;n++) 
		{
			X[n]+=alp*P[n];	// X(k+1)=X(k)+alp*P
			r[n]-=alp*AP[n];// r=r-alp*AP
			E+=r[n]*r[n];	//�덷
		}
		E=sqrt(E);
		//////////////////////////////
		
		//////////////// r=r-alp*AP
		//for(int n=0;n<pn;n++) r[n]-=alp*AP[n];
		/////////////////////////////
		
		/*/////////////////�덷
		E=0;
		for(int n=0;n<pn;n++) E+=r[n]*r[n];		
		E=sqrt(E);
		//cout<<"E="<<E<<endl;
		printf("�c�� E=%lf          \r",E);
		///////////////////////*/
		ofstream Ee("PICCG.dat", ios::app);
		Ee<<count<<" "<<E<<endl;
		Ee.close();
		///////////////////////beta
		beta=1.0/rLDLt_r;
		rLDLt_r=0;
		
        /////////////////y[i]�����Ƃ߁ALDLt_r[i]�����Ƃ߂�B
		for(int i=0;i<pn;i++)
		{
			if(i==0) y[0]=r[0]/val2[0]; //���i3.77�j �V
			else
			{
			    double sum=0;
			    for(int m=ptr2[i];m<ptr2[i+1]-1;m++)//�Ίp�����͏�������ptr[i+1]-1
			    {
			        sum+=val2[m]*y[ind2[m]];//���i3.78�j
			    }
			    int m=ptr2[i+1]-1;
			    y[i]=(r[i]-sum)/val2[m];
			}
		}////y[i]�����Ƃ܂����B
	
		/////////LDLt_r[i]�����߂�
		for(int i=pn-1;i>=0;i--)
		{
		    double sum=0;
			for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
			
		    LDLt_r[i]=y[i]-D1[i]*sum;	
		}
		/////////////////*/
	
		for(int n=0;n<pn;n++) rLDLt_r+=r[n]*LDLt_r[n];
		
		beta=beta*rLDLt_r;
		/////////////////*/
		
		///////////////////// P=r+beta*P
		for(int n=0;n<pn;n++) P[n]=LDLt_r[n]+beta*P[n];//iccg
	}
	cout<<"������="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;
	
	
    delete [] r;
	delete [] AP;
	delete [] P;

	delete [] y;
	delete [] LDLt_r;
	delete [] D1;

	delete [] val2;
	delete [] ind2;
	delete [] ptr2;

	for(int i=0;i<pn;i++)
	{
		delete [] VAL[i];
		delete [] IND[i];
	}
	delete [] VAL;
	delete [] IND;
	delete [] NUM;
}

///�d�����x�v�Z�֐�(�ӗv�f�p)
void calc_current(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,int *jnb,int *branch_num,double **current,int *depth,double II)
{
	//II:�d��[A]

	cout<<"�d�����x�v�Z�J�n"<<endl;
	
	double p=1.68e-8;//���̓d�C��R��[��m]

	int side_num2=0;//�R�C�����\������Ӑ�
	for(int i=1;i<=side_num;i++)
	{
		int ia=SIDE[i].node[1];
		int ib=SIDE[i].node[2];
		if(NODE[ia].material==COIL && NODE[ib].material==COIL) side_num2++;
	}///side_num2�����Ƃ܂���

	int *side_id=new int [side_num2+1];//�R�C�����\������Ӕԍ��i�[
	side_num2=0;
	for(int i=1;i<=side_num;i++)
	{
		int ia=SIDE[i].node[1];
		int ib=SIDE[i].node[2];
		if(NODE[ia].material==COIL && NODE[ib].material==COIL)
		{
			side_num2++;
			side_id[side_num2]=i;//�Ӕԍ�i���i�[
		}
	}///�R�C�����\������Ӕԍ����Ǘ�
	//////////////////////////��ٕӂ̓d���Ɋւ��鋫�E������ݒ�
	for(int i=1;i<=nelm;i++)
	{
		if(depth[i]==1)//���̕\�ʂɗאڂ����C�v�f
		{
			for(int j=1;j<=4;j++)
			{
				int kelm=ELEM[i].elm[j];
				if(kelm!=0)
				{
					if(ELEM[kelm].material==COIL) 
					{
						///��j�ʂɗאڂ���O�p�`�ƌ����������Ă��钸�_�́A��j�Ԑߓ_�ł���
						int p=ELEM[i].node[j];//���E�O�p�ɑ����Ȃ��ߓ_
						for(int k=1;k<=6;k++)
						{
							
							int iside=ELEM[i].sides[k];
							//cout<<"iside="<<iside<<endl;
							int ia=SIDE[iside].node[1];
							int ib=SIDE[iside].node[2];
							//cout<<"ia="<<ia<<"ib="<<ib<<endl;
							if(ib<ia)
							{
								int temp=ia;
								ia=ib;
								ib=temp;
							}///����ŕK��ia<ib�ƂȂ���
							
							if(ia!=p && ib!=p)//��iside�̓R�C�����E��Ƃ�������
							{
								if(NODE[ia].boundary_condition==11 || NODE[ib].boundary_condition==11)
								{
									SIDE[iside].boundary_condition=11;//11���ЂƂł��܂񂾕ӂ͎��R���E����
									//cout<<(NODE[ia].r[A_X]+NODE[ib].r[A_X])/2<<" "<<(NODE[ia].r[A_Y]+NODE[ib].r[A_Y])/2<<" "<<(NODE[ia].r[A_Z]+NODE[ib].r[A_Z])/2<<endl;
								}
								else
								{
									SIDE[iside].boundary_condition=10;//T=0�ƂȂ��
								}
								//T=0�̖ʂ���ю��R���E�����ʂɌŒ苫�E�ӂ�ݒu����
								if((NODE[ia].boundary_condition==21 && NODE[ib].boundary_condition==22) ||(NODE[ia].boundary_condition==22 && NODE[ib].boundary_condition==23)||(NODE[ia].boundary_condition==23 && NODE[ib].boundary_condition==21)) SIDE[iside].boundary_condition=21;//�ӂ�21��22�̕���
								else if((NODE[ia].boundary_condition==22 && NODE[ib].boundary_condition==21) ||(NODE[ia].boundary_condition==21 && NODE[ib].boundary_condition==23)||(NODE[ia].boundary_condition==23 && NODE[ib].boundary_condition==22)) SIDE[iside].boundary_condition=22;//�ӂ�22��21�̕���
								
							}
						}
					}
				}
			}
		}
		if(depth[i]==5)//���̕\�ʂɗאڂ����C�v�f
		{
			for(int j=1;j<=4;j++)
			{
				int kelm=ELEM[i].elm[j];
				if(kelm!=0)
				{
					if(ELEM[kelm].material==COIL) 
					{
						///��j�ʂɗאڂ���O�p�`�ƌ����������Ă��钸�_�́A��j�Ԑߓ_�ł���
						int p=ELEM[i].node[j];//���E�O�p�ɑ����Ȃ��ߓ_
						for(int k=1;k<=6;k++)
						{
							
							int iside=ELEM[i].sides[k];
							//cout<<"iside="<<iside<<endl;
							int ia=SIDE[iside].node[1];
							int ib=SIDE[iside].node[2];
							//cout<<"ia="<<ia<<"ib="<<ib<<endl;
							if(ib<ia)
							{
								int temp=ia;
								ia=ib;
								ib=temp;
							}///����ŕK��ia<ib�ƂȂ���
							
							if(ia!=p && ib!=p)//��iside�̓R�C�����E��Ƃ�������
							{
								if(NODE[ia].boundary_condition==11 || NODE[ib].boundary_condition==11)
								{
									SIDE[iside].boundary_condition=11;//11���ЂƂł��܂񂾕ӂ͎��R���E����
									//cout<<(NODE[ia].r[A_X]+NODE[ib].r[A_X])/2<<" "<<(NODE[ia].r[A_Y]+NODE[ib].r[A_Y])/2<<" "<<(NODE[ia].r[A_Z]+NODE[ib].r[A_Z])/2<<endl;
								}
								else
								{
									SIDE[iside].boundary_condition=10;//T=0�ƂȂ��
								}
								//T=0�̖ʂ���ю��R���E�����ʂɌŒ苫�E�ӂ�ݒu����
								if((NODE[ia].boundary_condition==21 && NODE[ib].boundary_condition==22) ||(NODE[ia].boundary_condition==22 && NODE[ib].boundary_condition==23)||(NODE[ia].boundary_condition==23 && NODE[ib].boundary_condition==21)) SIDE[iside].boundary_condition=23;//�ӂ�21��22�̕���
								else if((NODE[ia].boundary_condition==22 && NODE[ib].boundary_condition==21) ||(NODE[ia].boundary_condition==21 && NODE[ib].boundary_condition==23)||(NODE[ia].boundary_condition==23 && NODE[ib].boundary_condition==22)) SIDE[iside].boundary_condition=24;//�ӂ�22��21�̕���
								
							}
						}
					}
				}
			}
		}


	}//�R�C���[�ʂɎ��R���E�A���ʂ�T=0�̌Œ苫�E��~�����B
	////////////////*/

	///���E�����o�́@���܂������Ȃ��Ƃ��ɂ݂�
	//data_avs_J_boundary(node,nelm,NODE,ELEM);

	for(int k=1;k<=side_num2;k++)
    {
		int i=side_id[k];
		//���R���E�����𖢒m������������
        if(SIDE[i].boundary_condition==11) SIDE[i].boundary_condition=0;
	}
	
	double *T=new double [side_num+1];//�d���x�N�g�����ݼ��
    int NN=0;//�f�B���N���^���E�Ӑ�
    int *dn=new int [side_num+1]; //�e�ӂ��f�B���N���^���E�̉��Ԗڂ��B�f�B���N���^���E��łȂ��Ȃ�side_num2+1���i�[
    double *PHAT=new double [CON.get_max_DN()];//�f�B���N���^���l
    double J1=II/0.0000896;//89.6
	for(int k=1;k<=side_num;k++)T[k]=0;
    ///�f�B���N���^���E��������
    for(int k=1;k<=side_num2;k++)
    {
		int i=side_id[k];
        if(SIDE[i].boundary_condition==10)
		{
			dn[i]=NN;//i�Ԗڂ̕ӂ�NN�Ԗڂ̃f�B���N�����E��
	        PHAT[NN]=0;
	        T[i]=0;
	        NN++;
		}
		else if(SIDE[i].boundary_condition==21)
		{    
			int ia=SIDE[i].node[1];
			int ib=SIDE[i].node[2];
			double X=NODE[ia].r[A_X]-NODE[ib].r[A_X];
			double Y=NODE[ia].r[A_Y]-NODE[ib].r[A_Y];
			double Z=NODE[ia].r[A_Z]-NODE[ib].r[A_Z];
			double L=sqrt(X*X+Y*Y+Z*Z);//�ӂ̒���
	        dn[i]=NN;
	        PHAT[NN]=II;
	        T[i]=II;
	        NN++;
		}
		else if(SIDE[i].boundary_condition==22)
		{   
			int ia=SIDE[i].node[1];
			int ib=SIDE[i].node[2];
			double X=NODE[ia].r[A_X]-NODE[ib].r[A_X];
			double Y=NODE[ia].r[A_Y]-NODE[ib].r[A_Y];
			double Z=NODE[ia].r[A_Z]-NODE[ib].r[A_Z];
			double L=sqrt(X*X+Y*Y+Z*Z);//�ӂ̒���
	        dn[i]=NN;
	        PHAT[NN]=-II;
	        T[i]=-II;
	        NN++;
		}
		else if(SIDE[i].boundary_condition==23)
		{    
			int ia=SIDE[i].node[1];
			int ib=SIDE[i].node[2];
			double X=NODE[ia].r[A_X]-NODE[ib].r[A_X];
			double Y=NODE[ia].r[A_Y]-NODE[ib].r[A_Y];
			double Z=NODE[ia].r[A_Z]-NODE[ib].r[A_Z];
			double L=sqrt(X*X+Y*Y+Z*Z);//�ӂ̒���
	        dn[i]=NN;
	        PHAT[NN]=II;
	        T[i]=II;
	        NN++;
		}
		else if(SIDE[i].boundary_condition==24)
		{   
			int ia=SIDE[i].node[1];
			int ib=SIDE[i].node[2];
			double X=NODE[ia].r[A_X]-NODE[ib].r[A_X];
			double Y=NODE[ia].r[A_Y]-NODE[ib].r[A_Y];
			double Z=NODE[ia].r[A_Z]-NODE[ib].r[A_Z];
			double L=sqrt(X*X+Y*Y+Z*Z);//�ӂ̒���
	        dn[i]=NN;
	        PHAT[NN]=-II;
	        T[i]=-II;
	        NN++;
		}
		else dn[i]=side_num2+1;
    }
	cout<<"�ިظڐ���"<<NN<<endl;
	/////////////*/

	
    //////int pn=side_num-NN;				///���m��
	int pn=side_num2-NN;				///���m��
    int *ppn=new int [pn];			//�s���n�Ԗڂ͕�ppn[n]
    int *npp=new int [side_num+1];	///�e�ӂ��s��̉��Ԗڂɂ����邩�B�ިظڌ^�̏ꍇ��pn+1���i�[
    int num=0; 
    for(int k=1;k<=side_num2;k++)
    {
		int i=side_id[k];
        if(SIDE[i].boundary_condition==0)//���m�� 
		{
			ppn[num]=i;
			npp[i]=num;
			num++;
		}
		else npp[i]=pn+1;
    }
    
    ////�s��̍ő啝�v�Z
    int mat_w=0;
	for(int k=1;k<=side_num2;k++)
	{
		int i=side_id[k];
		int ia=SIDE[i].node[1];
		int ib=SIDE[i].node[2];
		int width=branch_num[ia]+branch_num[ib];//�s��̕�
		if(width>mat_w) mat_w=width;
	}
	//mat_w*=5;
	////////////
	
    ////�z��m��
    double **G=new double *[pn+1];///�S�̍s��
    for(int i=1;i<=pn;i++) G[i]=new double [mat_w+1];
    int **ROW=new int *[pn+1]; ///�e�s�́A��[���v�f�̗�ԍ��L��
    for(int i=1;i<=pn;i++) ROW[i]=new int [mat_w+1];
    int *NUM=new int [pn+1]; ///�e�s�́A��[���v�f��
    
    for(int i=1;i<=pn;i++)//������
    {
        NUM[i]=0;
        for(int j=1;j<=mat_w;j++)
		{
			G[i][j]=0;
			ROW[i][j]=0;
		}
    }
    double *B=new double [pn];//���s��
    for(int i=0;i<pn;i++) B[i]=0;//������
    ////
    
    /////////�S�̍s����쐬����
    
    int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int table[6+1][2+1];//�v�f���\������ӂƂ��̗v�f���ߓ_�ԍ��֌W
	
    for(int je=1;je<=nelm;je++)
    {   
		if(ELEM[je].material==COIL)
		{
			//�Ӂ|�ߓ_ð��ٍ쐬
			for(int i=1;i<=6;i++)
			{
				int iside=ELEM[je].sides[i];
				int ia=SIDE[iside].node[1];
				int ib=SIDE[iside].node[2];
				for(int j=1;j<=4;j++)
				{
					if(ELEM[je].node[j]==ia) table[i][1]=j;
					else if(ELEM[je].node[j]==ib) table[i][2]=j;
				}
			}////////////////
		
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
			
			double Xs=0;//�v�f�̏d�S���W
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j];
				Ys+=Y[j];
				Zs+=Z[j];
			}
			Xs/=4;Ys/=4;Zs/=4;
			////////////////////////////
	
			///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
	        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
		
			double delta6=ELEM[je].volume;//�̐ς�6�{
		
			delta6=1/delta6;//�v�Z�ɕK�v�Ȃ̂͋t���Ȃ̂ł����Ŕ��]���Ă���
		
			double delta=ELEM[je].volume/6;//�{���̑̐�
		
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
				int m=j%4+1;
				int n=m%4+1;
		    
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);//i�������̂Ƃ�
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
				if(i%2!=0)//i����Ȃ�
				{
					b[i]*=-1;
				    c[i]*=-1;
					d[i]*=-1;
					e[i]*=-1;
				}
			}
			/////////
		
			////�v�f��ظ��쐬�J�n
			for(int i=1;i<=6;i++)
			{	
				int iside=ELEM[je].sides[i];//�v�fje�̕Ӕԍ�
				if(SIDE[iside].boundary_condition==0)///�ިظڂłȂ��B�܂薢�m�Ȃ�
				{   
					int I=npp[iside]+1;///��iside�͍s���I�Ԗ�
					int I1=SIDE[iside].node[1];//iside���\������2�_
					int I2=SIDE[iside].node[2];
					int k1=table[i][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
					int k2=table[i][2];
				    for(int j=1;j<=6;j++)
					{
						int jside=ELEM[je].sides[j];
						
						int u1=table[j][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
						int u2=table[j][2];
							
						if(SIDE[jside].boundary_condition==0)///�ިظڂłȂ��B�܂薢�m�Ȃ�
						{   
							int J=npp[jside]+1;///��jside�͍s���J�Ԗ�
							int flag=0;
							//if(J<=I){
							int J1=SIDE[jside].node[1];//jside���\������2�_
							int J2=SIDE[jside].node[2];
							for(int h=1;h<=NUM[I];h++)
							{
								if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2/3*p*delta6*delta6*delta6;
								    flag=1;
								}
							}
							if(flag==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    
								G[I][H]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2/3*p*delta6*delta6*delta6;
							    ROW[I][H]=J;
							}
							//}
						}
						////
						else //jside���ިظڌ^���E�ߓ_�Ȃ�
						{
						    int n=dn[jside];
						    B[I-1]-=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2/3*p*delta6*delta6*delta6*PHAT[n];
						}//////////*/
					}
				}
			}
		}
    }
    ///////////////////////*/
	

    int number=0;//�s��̔�[���v�f��
    for(int i=1;i<=pn;i++) number+=NUM[i];
   
    ///���̂܂܂ł�G[i][j]��[j]�̏��������ɕ���ł��Ȃ��̂ŁAG��ROW����ёւ���
    for(int i=1;i<=pn;i++)
    {
        double tempG;
		int tempR;
		for(int j=2;j<=NUM[i];j++)
		{
		    for(int m=1;m<j;m++)
		    {
		        if(ROW[i][j]<ROW[i][m])
				{
				    tempG=G[i][m];
				    tempR=ROW[i][m];
					G[i][m]=G[i][j];
					ROW[i][m]=ROW[i][j];
					G[i][j]=tempG;
					ROW[i][j]=tempR;
				}
			}
		}
    }///////////

	///�Ώ̐��`�F�b�N
	for(int i=1;i<=pn;i++)
	{
		for(int j=1;j<=NUM[i];j++)
		{
			int J=ROW[i][j];
			for(int k=1;k<=NUM[J];k++) if(ROW[J][k]==i) if(G[i][j]!=G[J][k])
			{
				cout<<"matrix isn't symmetric   "<<G[i][j]<<" "<<G[J][k]<<endl;
			}
		}
	}///////////*/

    double *val = new double [number];
    int *ind = new int [number];//��[���v�f�̗�ԍ��i�[�z��
    int *ptr = new int [pn+1];//�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
    
    /////////////////////val,ind ,ptr�ɒl���i�[
    int index=0;
    for(int n=0;n<pn;n++)
    {
        ptr[n]=index;
		for(int m=1;m<=NUM[n+1];m++)
		{
			val[index]=G[n+1][m];
			ind[index]=ROW[n+1][m]-1;
			index++;
		}
    }	    
    ptr[pn]=number;
    ////////////////////*/
    
    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
    
    cout<<"�s��쐬�I��  "<<endl;
    
    double *XX=new double [pn];//�s��̓����i�[
	//CG3D(val,ind,ptr,pn,ppn,B,T);//CG�@���s
	//ICCG3D(val,ind,ptr,pn,ppn,B,T,number);//ICCG�@���s
	ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	//MRTR(CON,B,pn,XX,val,ind,ptr);
    ///////////
	for(int n=0;n<pn;n++)
			{
				int i=ppn[n];
				T[i]=XX[n];
			}
			delete [] XX;
	denryu_side(CON,NODE,ELEM,SIDE,node,nelm,side_num,T,current,J1);
    
	delete [] side_id;
	delete [] T;
    ///////////////////////*/
    delete [] dn;
    delete [] PHAT;
    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;
}
////�ӗv�f�d�����x�v�Z�֐�
void denryu_side(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,double *T,double **current,double J1)
{
	int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int table[6+1][2+1];//�v�f���\������ӂƂ��̗v�f���ߓ_�ԍ��֌W

	ofstream fp("j.dat");
	for(int je=1;je<=nelm;je++)
    {   
		if(ELEM[je].material==COIL)
		{
			//�Ӂ|�ߓ_ð��ٍ쐬
			for(int i=1;i<=6;i++)
			{
				int iside=ELEM[je].sides[i];
				int ia=SIDE[iside].node[1];
				int ib=SIDE[iside].node[2];
				for(int j=1;j<=4;j++)
				{
					if(ELEM[je].node[j]==ia) table[i][1]=j;
					else if(ELEM[je].node[j]==ib) table[i][2]=j;
				}
			}////////////////
	
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
			double Xs=0;//�v�f�̏d�S���W
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j];
				Ys+=Y[j];
				Zs+=Z[j];
			}
			Xs/=4;
			Ys/=4;
			Zs/=4;
			////////////////////////////
	
			double delta6=ELEM[je].volume;//�̐ς�6�{(�������̐ς̒l�͂��łɃx�N�g���|�e���V���������߂�ۂɌv�Z���Ă���)
		
			delta6=1/delta6;//�v�Z�ɕK�v�Ȃ̂͋t���Ȃ̂ł����Ŕ��]���Ă���
		
			double delta=ELEM[je].volume/6;//�{���̑̐�
		
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
				int m=j%4+1;
				int n=m%4+1;
		    
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
				if(i%2!=0)//i����Ȃ�
				{
				    c[i]*=-1;
					d[i]*=-1;
					e[i]*=-1;
				}
			}
			/////////
	
			for(int D=0;D<3;D++) current[D][je]=0;//������
	
			for(int i=1;i<=6;i++)
			{
				int s=ELEM[je].sides[i];
				int k1=table[i][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
				int k2=table[i][2];
	
				current[A_X][je]+=(d[k1]*e[k2]-e[k1]*d[k2])*T[s];
				current[A_Y][je]+=(e[k1]*c[k2]-c[k1]*e[k2])*T[s];
				current[A_Z][je]+=(c[k1]*d[k2]-d[k1]*c[k2])*T[s];
				
			}
	
			for(int D=0;D<3;D++) current[D][je]*=delta6*delta6*2;
			///////////////////////////////////////////////////////////
			double Along=sqrt(pow(current[A_X][je],2)+pow(current[A_Y][je],2)+pow(current[A_Z][je],2));
			current[A_X][je]=current[A_X][je]*J1/Along;
			current[A_Y][je]=current[A_Y][je]*J1/Along;
			current[A_Z][je]=current[A_Z][je]*J1/Along;

			//////////////////////////////////////////////////////////*/
			//if(Zs>0 && Zs<0.0001) fp<<Xs<<" "<<Ys<<" "<<current[A_X][je]*1e-12<<" "<<current[A_Y][je]*1e-12<<endl;
			if(Zs>0.0005 && Zs<0.001)
			fp<<Xs<<" "<<Ys<<" "<<current[A_X][je]*1e-12<<" "<<current[A_Y][je]*1e-12<<endl;		
		}
		
	}	
	fp.close();
	
}

//�v�f�̐[������֐�
void set_depth(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,int *depth,int KTE)
{
	//for(int i=1;i<=KTE;i++) depth[i]=0;
	for(int i=1;i<=nelm;i++) depth[i]=0;

	for(int i=1;i<=nelm;i++)//�����ߓ_���ЂƂł��܂ދ�C�v�f�͋��E�v�f(�[��1)�ƒ�`����
	{
		if(ELEM[i].material==AIR)//�v�fi����C
		{
			for(int j=1;j<=4;j++)
			{
				int ia=ELEM[i].node[j];
				if(NODE[ia].material!=AIR)
				{
					depth[i]=1;//�����v�f�ɐڂ����C�v�f��depth���P�ƒ�`
				}
			}
		}
		else if(ELEM[i].material==IRON)	//�R�C���p
		{
			{
			for(int j=1;j<=4;j++)
			{
				int ia=ELEM[i].node[j];
				if(NODE[ia].material==COIL)
				{
					depth[i]=5;//�����v�f�ɐڂ����C�v�f��depth���P�ƒ�`
				}
			}
		}
		}
	}///�[���P�̗v�f�����Ƃ܂���(ϸ���ِϕ��ʂ͐[��1�̗v�f���[��2�Ɛڂ���ʂƂ���΂悢)
	
	int count=10;
	int A=1;//���ݒ��ڂ���[��
	while(count>0)
	{
		count=0;//�V�����[�������Ƃ܂����v�f��
		for(int i=1;i<=nelm;i++)
		{
			if(ELEM[i].material==AIR)
			{
				if(depth[i]==A)
				{
					for(int j=1;j<=4;j++)
					{
						int jelm=ELEM[i].elm[j];//i�̗אڂ���v�f
						if(jelm!=0)
						{
							if(ELEM[jelm].material==AIR && depth[jelm]==0)
							{
								depth[jelm]=A+1;
								count++;
							}
						}
					}
				}
			}
		}
		A++;
	}////�Q�ȍ~�̐[���Ɋւ��ẮA�����ł̒�`�Ƃ̂��Ɏ������������Ƃ��̒�`�ł͏����قȂ�̂Œ��ӁB
	///check
	for(int i=1;i<=nelm;i++) if(ELEM[i].material==AIR) if(depth[i]==0) cout<<"�v�f�[��error"<<endl;
}

void carrent_vector(mpsconfig &CON, vector<point3D> &NODE, vector<element3D> &ELEM, int nelm, double **B,int t)
{
	double **center_of_element=new double* [nelm+1];
	for(int i=1;i<=nelm;i++) center_of_element[i]=new double[3];
	//�d�S�̌v�Z�E�E�Enelm�S�ďo��
	for(int i=1;i<=nelm;i++)
	{
		for(int D=0;D<3;D++)
		{
			center_of_element[i][D]=0.0;
			for(int j=1;j<=4;j++) center_of_element[i][D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
		}
	}

	int timestep=CON.get_current_step();
	int half_nelm=0;
	for(int i=1;i<=nelm;i++){				//Y-Z�ʂɕ\������f�ʂ܂ł̗v�f��
		if(center_of_element[i][A_X]<0.0) half_nelm++;
	}
//	if(t==1 || t%10==0){
		stringstream sscv;
		sscv<<"./Current/Cvector"<<t<<".fld";
		string Cvector=sscv.str();	
		ofstream fcv(Cvector);
		if(fcv.fail()){
			system("mkdir Current");
		ofstream fcv(Cvector);
		if(fcv.fail()){
			cout<<"�t�@�C�����J���܂���ł����BCurrent�t�H���_�����邩�ǂ����m�F���Ă�������"<<endl;
			exit(1);
			}
		}

		fcv << "# AVS field file" << endl;
		fcv << "ndim=1" << endl;
		fcv << "dim1=" << half_nelm <<endl;
		fcv << "nspace=3" << endl;
		fcv << "veclen=3" << endl;
		fcv << "data=float" << endl;
		fcv << "field=irregular" << endl;
		fcv << "label=e-x e-y e-z" << endl << endl;
		fcv << "variable 1 file=./Cvector"<<t<<" filetype=ascii skip=1 offset=0 stride=6" << endl;
		fcv << "variable 2 file=./Cvector"<<t<<" filetype=ascii skip=1 offset=1 stride=6" << endl;
		fcv << "variable 3 file=./Cvector"<<t<<" filetype=ascii skip=1 offset=2 stride=6" << endl;
		fcv << "coord    1 file=./Cvector"<<t<<" filetype=ascii skip=1 offset=3 stride=6" << endl;
		fcv << "coord    2 file=./Cvector"<<t<<" filetype=ascii skip=1 offset=4 stride=6" << endl;
		fcv << "coord    3 file=./Cvector"<<t<<" filetype=ascii skip=1 offset=5 stride=6" << endl;

		fcv.close();

		stringstream sscvd;
		sscvd<<"./Current/Cvector"<<t;
		string Cvectordata=sscvd.str();

		ofstream fout(Cvectordata);
		if(fout.fail()){
			cout<<"�f�[�^�t�@�C�����J���܂���ł����B"<<endl;
			exit(1);
		}

		fout<<"e-x e-y e-z x y z"<<endl;
		for(int i=1;i<=nelm;i++)
		{
			if(center_of_element[i][A_X]<0){
				fout<<B[A_X][i]/1000000<<" "<<B[A_Y][i]/1000000<<" "<<B[A_Z][i]/1000000<<" "<<center_of_element[i][A_X]<<" "<<center_of_element[i][A_Y]<<" "<<center_of_element[i][A_Z]<<endl;
			}
		}
		fout.close();
//	}//*/
		for(int i=1;i<=nelm;i++) delete [] center_of_element[i];
		delete [] center_of_element;
}
///����`�޸�����ݼ�ٌv�Z�֐�(�ӗv�f�p)(����)
void non_linear_Avector3D(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,double *V,int *jnb,int *branch_num,double **current,double *RP,double II,int *depth,double **Be,int t)
{
	//////////////////////////////////�d�����x�x�N�g���v�Z///////////////////////////////////////////////
	//////////////////////��ق��Ȃ���Γd���v�Z���s���K�v���Ȃ��B�����ź�ِߓ_������������
	int coil_node_num=0;	//��ِߓ_��
	for(int i=1;i<=node;i++) if(NODE[i].material==COIL) coil_node_num++;
	////////////////////////*/

	int *save_bound=new int [node+1];//�{�֐��͉��x���Ăяo�����̂ŁA�d���Ɋւ��鋫�E���������S�ɏ����ł��Ȃ��B�����ŕۑ�����
	
	if(coil_node_num>0)///��ِߓ_������Ȃ�d�����x�v�Z
	{
		for(int i=1;i<=node;i++) save_bound[i]=NODE[i].boundary_condition;//���E������ۑ�

		if(II!=0)//�d���l����[���Ȃ�d���v�Z
		{
			calc_current(CON,NODE,ELEM,SIDE,node,nelm,side_num,jnb,branch_num,current,depth,II);
			cout<<"�d���v�Z����"<<endl;
		}
		else	//�[���Ȃ珉����
		{
			for(int i=1;i<=nelm;i++) if(ELEM[i].material==COIL) for(int D=0;D<3;D++) current[D][i]=0;
		}
		///�d�����E�����̏������i�����d���͂��Ƃ܂��Ă��邩�炢��Ȃ�)
		for(int i=1;i<=side_num;i++) if(SIDE[i].boundary_condition>=10) SIDE[i].boundary_condition=0;
		for(int i=1;i<=node;i++) if(NODE[i].boundary_condition>=10) NODE[i].boundary_condition=0;
	}
	carrent_vector(CON, NODE, ELEM, nelm, current,t);
	///////////////////////////////////////////////////////////////////////////////////////////////////


	//�����ł́AV[i]�͊e�ӂ��޸�����ݼ�ق̔����ω���A��\�����Ƃɒ��ӁB
	cout<<"����`�޸�����ݼ�ٌv�Z�J�n"<<endl;
	double u0=4*PI*0.0000001;	//�^��̓�����
    double v0=1/u0;///���C��R��
	double j0x, j0y, j0z;//�d�����x[A/m^3]
	double MA=CON.get_magnet_angle();
	double magnet_direction[3]={-sin(MA*2*PI/360),0,cos(MA*2*PI/360)};
	double Mx=CON.get_magnet_B()*magnet_direction[A_X];
	double My=CON.get_magnet_B()*magnet_direction[A_Y];
	double Mz=CON.get_magnet_B()*magnet_direction[A_Z];
	int graphtype=2;

	//���͕ǂɌŒ苫�E������ݒ�
	for(int i=1;i<=nelm; i++){
		if(ELEM[i].material==AIR){//���̕\�ʂɗאڂ����C�v�f
			for(int j=1; j<=4; j++){
				int kelm=ELEM[i].elm[j];
				if(kelm==0){
					//��j�ʂɗאڂ���O�p�`�ƌ����������Ă��钸�_�́C��j�Ԑړ_�ł���
					int p=ELEM[i].node[j];//���E�O�p�`�ɑ����Ȃ��ړ_
					for(int k=1;k<=6;k++){
						int iside=ELEM[i].sides[k];
						int ia=SIDE[iside].node[1];
						int ib=SIDE[iside].node[2];
						if(ia!=p && ib!=p) SIDE[iside].boundary_condition=1;//�ړ_p���܂܂Ȃ��ӂ͋��E��
						else SIDE[iside].boundary_condition=0;
					}
				}
			}
		}
	}

	//////////////////////////////////////�A�L�}�⊮�ɂ��B-H�Ȑ��쐬/////////////////////////////////
	////////////////////����`������
	double H[14];//Nn+4
	double M[14];
	double Bflux[14];
	int Nn=10;		//�f�[�^��
	if(graphtype==0)//�T�C�g
	{
		H[2]=0.0;
		H[3]=17391;
		H[4]=43478;
		H[5]=78260;
		H[6]=173913;
		H[7]=391304;
		M[2]=0.0;
		M[3]=0.01;
		M[4]=0.015;
		M[5]=0.02;
		M[6]=0.025;
		M[7]=0.03;
		for(int n=2;n<Nn+2;n++) Bflux[n]=u0*H[n]+M[n];
	}
	else if(graphtype==1)//�����W���o��
	{
		double Ms=14700;			//�O�a����[A/m]
		double kai0=1.172;			//���C����
		double gam=3*kai0/Ms;
		H[2]=0.0;
		H[3]=2000;
		H[4]=4000;
		H[5]=6000;
		H[6]=8000;
		H[7]=12000;
		M[2]=0.0;
		for(int k=3;k<=7;k++) M[k]=Ms*(1/tanh(gam*H[k])-1/(gam*H[k]));
		for(int n=2;n<Nn+2;n++) Bflux[n]=u0*H[n]+u0*M[n];
	}
	else if(graphtype==2)
	{
		H[2]=29.083;
		H[3]=26119.17;
		H[4]=88160.77;
		H[5]=212216.5;
		H[6]=428597.7;
		H[7]=800722.5;
		H[8]=1020144;
		H[9]=1212705;
		H[10]=1483964;
		H[11]=1556321;
		Bflux[2]=0.0;
		Bflux[3]=0.102173;
		Bflux[4]=0.389118;
		Bflux[5]=0.859961;
		Bflux[6]=1.40583;
		Bflux[7]=1.930318;
		Bflux[8]=2.268766;
		Bflux[9]=2.523716;
		Bflux[10]=2.881576;
		Bflux[11]=2.976717;
	}

	///���[��2�_�ǉ�
	Bflux[0]=Bflux[2]-(Bflux[4]-Bflux[2]);
	Bflux[1]=Bflux[2]-(Bflux[3]-Bflux[2]);
	H[0]=H[2]*(Bflux[0]-Bflux[3])*(Bflux[0]-Bflux[4])/((Bflux[2]-Bflux[3])*(Bflux[2]-Bflux[4]))+H[3]*(Bflux[0]-Bflux[2])*(Bflux[0]-Bflux[4])/((Bflux[3]-Bflux[2])*(Bflux[3]-Bflux[4]))+H[4]*(Bflux[0]-Bflux[2])*(Bflux[0]-Bflux[3])/((Bflux[4]-Bflux[2])*(Bflux[4]-Bflux[3]));
	H[1]=H[2]*(Bflux[1]-Bflux[3])*(Bflux[1]-Bflux[4])/((Bflux[2]-Bflux[3])*(Bflux[2]-Bflux[4]))+H[3]*(Bflux[1]-Bflux[2])*(Bflux[1]-Bflux[4])/((Bflux[3]-Bflux[2])*(Bflux[3]-Bflux[4]))+H[4]*(Bflux[1]-Bflux[2])*(Bflux[1]-Bflux[3])/((Bflux[4]-Bflux[2])*(Bflux[4]-Bflux[3]));

	Bflux[Nn+2]=Bflux[Nn+1]+(Bflux[Nn+1]-Bflux[Nn]);
	Bflux[Nn+3]=Bflux[Nn+1]+(Bflux[Nn+1]-Bflux[Nn-1]);
	H[Nn+2]=H[Nn-1]*(Bflux[Nn+2]-Bflux[Nn])*(Bflux[Nn+2]-Bflux[Nn+1])/((Bflux[Nn-1]-Bflux[Nn])*(Bflux[Nn-1]-Bflux[Nn+1]))+H[Nn]*(Bflux[Nn+2]-Bflux[Nn-1])*(Bflux[Nn+2]-Bflux[Nn+1])/((Bflux[Nn]-Bflux[Nn-1])*(Bflux[Nn]-Bflux[Nn+1]))+H[Nn+1]*(Bflux[Nn+2]-Bflux[Nn-1])*(Bflux[Nn+2]-Bflux[Nn])/((Bflux[Nn+1]-Bflux[Nn-1])*(Bflux[Nn+1]-Bflux[Nn]));
	H[Nn+3]=H[Nn-1]*(Bflux[Nn+3]-Bflux[Nn])*(Bflux[Nn+3]-Bflux[Nn+1])/((Bflux[Nn-1]-Bflux[Nn])*(Bflux[Nn-1]-Bflux[Nn+1]))+H[Nn]*(Bflux[Nn+3]-Bflux[Nn-1])*(Bflux[Nn+3]-Bflux[Nn+1])/((Bflux[Nn]-Bflux[Nn-1])*(Bflux[Nn]-Bflux[Nn+1]))+H[Nn+1]*(Bflux[Nn+3]-Bflux[Nn-1])*(Bflux[Nn+3]-Bflux[Nn])/((Bflux[Nn+1]-Bflux[Nn-1])*(Bflux[Nn+1]-Bflux[Nn]));
	////////////

	/*/////////////////�Ώ̐����l���ĕ␳
	H[0]=-H[4];
	H[1]=-H[3];
	Bflux[0]=-Bflux[4];
	Bflux[1]=-Bflux[3];
	//////////////////////*/

	/////���[�U�[�̓��͂����_(B-H�Ȑ����v���b�g
	ofstream fp3("B-H.dat");
	for(int i=0;i<Nn+4;i++) fp3<<H[i]<<" "<<Bflux[i]<<endl;
	fp3.close();
	//////////////////*/

	///�e��Ԃ�H-B�Ȑ��̌X�����v�Z
	double *mm=new double [Nn+3];//�e��Ԃ̌X���B
	for(int i=0;i<Nn+3;i++) mm[i]=(H[i+1]-H[i])/(Bflux[i+1]-Bflux[i]);

	///���[�U�[�̓��͓_�ɂ�����H-B�Ȑ��̌X��,�܂莥�C��R���v�Z
	double *vm=new double[Nn+4];//�X��
	for(int i=2;i<Nn+2;i++)
	{
		double a1=sqrt((mm[i+1]-mm[i])*(mm[i+1]-mm[i]));
		double a2=sqrt((mm[i-1]-mm[i-2])*(mm[i-1]-mm[i-2]));
		double A0=a1+a2;

		if(A0<0.0001) vm[i]=0.5*(mm[i-1]+mm[i]);
		else vm[i]=(a1*mm[i-1]+a2*mm[i])/A0;
		
	}///////////

	ofstream fout("H-B_akima.dat");//�A�L�}��Ԃɂ��H-B�Ȑ���ۯ�
	for(double i=Bflux[2];i<=Bflux[Nn+1];i+=0.001)
	{
		for(int n=2;n<=Nn+1;n++)
		{
			if(i>=Bflux[n] && i<Bflux[n+1])
			{
				double h=Bflux[n+1]-Bflux[n];
				double aj0=H[n];
				double aj1=vm[n];
				double aj2=(3*mm[n]-2*vm[n]-vm[n+1])/h;
				double aj3=(vm[n]+vm[n+1]-2*mm[n])/(h*h);
				double val=aj0+aj1*(i-Bflux[n])+aj2*(i-Bflux[n])*(i-Bflux[n])+aj3*(i-Bflux[n])*(i-Bflux[n])*(i-Bflux[n]);
				double val2=aj1+2*aj2*(i-Bflux[n])+3*aj3*(i-Bflux[n])*(i-Bflux[n]);
				fout<<i<<" "<<val<<endl;
			}
		}
	}
	fout.close();
	///////////////////////////*/

	//���[2�_��vm�����߂�
	vm[0]=vm[2]*(Bflux[0]-Bflux[3])*(Bflux[0]-Bflux[4])/((Bflux[2]-Bflux[3])*(Bflux[2]-Bflux[4]))+vm[3]*(Bflux[0]-Bflux[2])*(Bflux[0]-Bflux[4])/((Bflux[3]-Bflux[2])*(Bflux[3]-Bflux[4]))+vm[4]*(Bflux[0]-Bflux[2])*(Bflux[0]-Bflux[3])/((Bflux[4]-Bflux[2])*(Bflux[4]-Bflux[3]));
	vm[1]=vm[2]*(Bflux[1]-Bflux[3])*(Bflux[1]-Bflux[4])/((Bflux[2]-Bflux[3])*(Bflux[2]-Bflux[4]))+vm[3]*(Bflux[1]-Bflux[2])*(Bflux[1]-Bflux[4])/((Bflux[3]-Bflux[2])*(Bflux[3]-Bflux[4]))+vm[4]*(Bflux[1]-Bflux[2])*(Bflux[1]-Bflux[3])/((Bflux[4]-Bflux[2])*(Bflux[4]-Bflux[3]));
	vm[Nn+2]=vm[Nn-1]*(Bflux[Nn+2]-Bflux[Nn])*(Bflux[Nn+2]-Bflux[Nn+1])/((Bflux[Nn-1]-Bflux[Nn])*(Bflux[Nn-1]-Bflux[Nn+1]))+vm[Nn]*(Bflux[Nn+2]-Bflux[Nn-1])*(Bflux[Nn+2]-Bflux[Nn+1])/((Bflux[Nn]-Bflux[Nn-1])*(Bflux[Nn]-Bflux[Nn+1]))+vm[Nn+1]*(Bflux[Nn+2]-Bflux[Nn-1])*(Bflux[Nn+2]-Bflux[Nn])/((Bflux[Nn+1]-Bflux[Nn-1])*(Bflux[Nn+1]-Bflux[Nn]));
	vm[Nn+3]=vm[Nn-1]*(Bflux[Nn+3]-Bflux[Nn])*(Bflux[Nn+3]-Bflux[Nn+1])/((Bflux[Nn-1]-Bflux[Nn])*(Bflux[Nn-1]-Bflux[Nn+1]))+vm[Nn]*(Bflux[Nn+3]-Bflux[Nn-1])*(Bflux[Nn+3]-Bflux[Nn+1])/((Bflux[Nn]-Bflux[Nn-1])*(Bflux[Nn]-Bflux[Nn+1]))+vm[Nn+1]*(Bflux[Nn+3]-Bflux[Nn-1])*(Bflux[Nn+3]-Bflux[Nn])/((Bflux[Nn+1]-Bflux[Nn-1])*(Bflux[Nn+1]-Bflux[Nn]));
	
	//���[�U�[�̓��͓_�ɂ�����v-B�Ȑ���ۯ�
	ofstream fp("v-B.dat");
	for(int i=0;i<Nn+4;i++) fp<<Bflux[i]<<" "<<vm[i]<<endl;
	fp.close();
	/////////////*/

	///�e��Ԃ�v-B�Ȑ��̌X�����v�Z
	double *mm2=new double [Nn+3];//�e��Ԃ̌X���B
	for(int i=0;i<Nn+3;i++) mm2[i]=(vm[i+1]-vm[i])/(Bflux[i+1]-Bflux[i]);

	///���[�U�[�̓��͓_�ɂ�����v-B�Ȑ��̌X���v�Z
	double *t2=new double[Nn+4];//�X��
	for(int i=2;i<Nn+2;i++)
	{
		double a1=sqrt((mm2[i+1]-mm2[i])*(mm2[i+1]-mm2[i]));
		double a2=sqrt((mm2[i-1]-mm2[i-2])*(mm2[i-1]-mm2[i-2]));
		double A=a1+a2;

		if(A<0.0001) t2[i]=0.5*(mm2[i-1]+mm2[i]);
		else t2[i]=(a1*mm2[i-1]+a2*mm2[i])/A;
	}///////////

	/////�A�L�}��Ԃɂ��v-B�Ȑ���ۯ�
	ofstream fout2("v-B_akima.dat");
	for(double i=Bflux[2];i<=Bflux[Nn+1];i+=0.001)
	{
		for(int n=2;n<=Nn+1;n++)
		{
			if(i>=Bflux[n] && i<Bflux[n+1])
			{
				double h=Bflux[n+1]-Bflux[n];
				double aj0=vm[n];
				double aj1=t2[n];
				double aj2=(3*mm2[n]-2*t2[n]-t2[n+1])/h;
				double aj3=(t2[n]+t2[n+1]-2*mm2[n])/(h*h);
				double val=aj0+aj1*(i-Bflux[n])+aj2*(i-Bflux[n])*(i-Bflux[n])+aj3*(i-Bflux[n])*(i-Bflux[n])*(i-Bflux[n]);
				fout2<<i<<" "<<val<<endl;
			}
		}
	}
	fout2.close();
	/////////*/

	double *Bflux2=new double [Nn+4];//�a^2�i�[
	for(int i=0;i<Nn+4;i++) Bflux2[i]=Bflux[i]*Bflux[i];

	///�e��Ԃ�v-B^2�Ȑ��̌X�����v�Z
	double *mm3=new double [Nn+3];//�e��Ԃ̌X���B
	for(int i=0;i<Nn+3;i++) mm3[i]=(vm[i+1]-vm[i])/(Bflux2[i+1]-Bflux2[i]);

	///���[�U�[�̓��͓_�ɂ�����v-B^2�Ȑ��̌X����v/��B^2�v�Z
	double *dvdB2=new double[Nn+4];//�X��
	for(int i=2;i<Nn+2;i++)
	{
		double a1=sqrt((mm3[i+1]-mm3[i])*(mm3[i+1]-mm3[i]));
		double a2=sqrt((mm3[i-1]-mm3[i-2])*(mm3[i-1]-mm3[i-2]));
		double A=a1+a2;

		if(A<0.0001) dvdB2[i]=0.5*(mm3[i-1]+mm3[i]);
		else dvdB2[i]=(a1*mm3[i-1]+a2*mm3[i])/A;
	}///////////

	/////�A�L�}��Ԃɂ��v-B^2�Ȑ���ۯ�
	ofstream fout3("v-B2_akima.dat");
	for(double i=Bflux2[2];i<=Bflux2[Nn+1];i+=0.0001)
	{
		for(int n=2;n<=Nn+1;n++)
		{
			if(i>=Bflux2[n] && i<Bflux2[n+1])
			{
				double h=Bflux2[n+1]-Bflux2[n];
				double aj0=vm[n];
				double aj1=dvdB2[n];
				double aj2=(3*mm3[n]-2*dvdB2[n]-dvdB2[n+1])/h;
				double aj3=(dvdB2[n]+dvdB2[n+1]-2*mm3[n])/(h*h);
				double val=aj0+aj1*(i-Bflux2[n])+aj2*(i-Bflux2[n])*(i-Bflux2[n])+aj3*(i-Bflux2[n])*(i-Bflux2[n])*(i-Bflux2[n]);
				fout3<<i<<" "<<val<<endl;
			}
		}
	}
	fout3.close();
	/////////*/

	//���[2�_��dvdB2�����߂�
	dvdB2[0]=dvdB2[2]*(Bflux2[0]-Bflux2[3])*(Bflux2[0]-Bflux2[4])/((Bflux2[2]-Bflux2[3])*(Bflux2[2]-Bflux2[4]))+dvdB2[3]*(Bflux2[0]-Bflux2[2])*(Bflux2[0]-Bflux2[4])/((Bflux2[3]-Bflux2[2])*(Bflux2[3]-Bflux2[4]))+dvdB2[4]*(Bflux2[0]-Bflux2[2])*(Bflux2[0]-Bflux2[3])/((Bflux2[4]-Bflux2[2])*(Bflux2[4]-Bflux2[3]));
	dvdB2[1]=dvdB2[2]*(Bflux2[1]-Bflux2[3])*(Bflux2[1]-Bflux2[4])/((Bflux2[2]-Bflux2[3])*(Bflux2[2]-Bflux2[4]))+dvdB2[3]*(Bflux2[1]-Bflux2[2])*(Bflux2[1]-Bflux2[4])/((Bflux2[3]-Bflux2[2])*(Bflux2[3]-Bflux2[4]))+dvdB2[4]*(Bflux2[1]-Bflux2[2])*(Bflux2[1]-Bflux2[3])/((Bflux2[4]-Bflux2[2])*(Bflux2[4]-Bflux2[3]));
	dvdB2[Nn+2]=dvdB2[Nn-1]*(Bflux2[Nn+2]-Bflux2[Nn])*(Bflux2[Nn+2]-Bflux2[Nn+1])/((Bflux2[Nn-1]-Bflux2[Nn])*(Bflux2[Nn-1]-Bflux2[Nn+1]))+dvdB2[Nn]*(Bflux2[Nn+2]-Bflux2[Nn-1])*(Bflux2[Nn+2]-Bflux2[Nn+1])/((Bflux2[Nn]-Bflux2[Nn-1])*(Bflux2[Nn]-Bflux2[Nn+1]))+dvdB2[Nn+1]*(Bflux2[Nn+2]-Bflux2[Nn-1])*(Bflux2[Nn+2]-Bflux2[Nn])/((Bflux2[Nn+1]-Bflux2[Nn-1])*(Bflux2[Nn+1]-Bflux2[Nn]));
	dvdB2[Nn+3]=dvdB2[Nn-1]*(Bflux2[Nn+3]-Bflux2[Nn])*(Bflux2[Nn+3]-Bflux2[Nn+1])/((Bflux2[Nn-1]-Bflux2[Nn])*(Bflux2[Nn-1]-Bflux2[Nn+1]))+dvdB2[Nn]*(Bflux2[Nn+3]-Bflux2[Nn-1])*(Bflux2[Nn+3]-Bflux2[Nn+1])/((Bflux2[Nn]-Bflux2[Nn-1])*(Bflux2[Nn]-Bflux2[Nn+1]))+dvdB2[Nn+1]*(Bflux2[Nn+3]-Bflux2[Nn-1])*(Bflux2[Nn+3]-Bflux2[Nn])/((Bflux2[Nn+1]-Bflux2[Nn-1])*(Bflux2[Nn+1]-Bflux2[Nn]));
	
	///�e��Ԃ́�v/��B^2-B^2�Ȑ��̌X�����v�Z
	double *mm4=new double [Nn+3];//�e��Ԃ̌X���B
	for(int i=0;i<Nn+3;i++) mm4[i]=(dvdB2[i+1]-dvdB2[i])/(Bflux2[i+1]-Bflux2[i]);

	///���[�U�[�̓��͓_�ɂ������v/��B^2-B^2�Ȑ��̌X���v�Z
	double *t4=new double[Nn+4];//�X��
	for(int i=2;i<Nn+2;i++)
	{
		double a1=sqrt((mm4[i+1]-mm4[i])*(mm4[i+1]-mm4[i]));
		double a2=sqrt((mm4[i-1]-mm4[i-2])*(mm4[i-1]-mm4[i-2]));
		double A=a1+a2;

		if(A<0.0001) t4[i]=0.5*(mm4[i-1]+mm4[i]);
		else t4[i]=(a1*mm4[i-1]+a2*mm4[i])/A;
	}///////////

	/////�A�L�}��Ԃɂ���v/��B^2-B^2�Ȑ���ۯ�
	ofstream fout4("dvdB2-B2_akima.dat");
	for(double i=Bflux2[2];i<=Bflux2[Nn+1];i+=0.001)
	{
		for(int n=2;n<=Nn+1;n++)
		{
			if(i>=Bflux2[n] && i<Bflux2[n+1])
			{
				double h=Bflux2[n+1]-Bflux2[n];
				double aj0=dvdB2[n];
				double aj1=t4[n];
				double aj2=(3*mm4[n]-2*t4[n]-t4[n+1])/h;
				double aj3=(t4[n]+t4[n+1]-2*mm4[n])/(h*h);
				double val=aj0+aj1*(i-Bflux2[n])+aj2*(i-Bflux2[n])*(i-Bflux2[n])+aj3*(i-Bflux2[n])*(i-Bflux2[n])*(i-Bflux2[n]);
				fout4<<i<<" "<<val<<endl;
			}
		}
	}
	fout4.close();
	/////////*/
	double *A=new double [side_num+1];
    int NN=0;//�f�B���N���^���E�Ӑ�
    int *dn=new int [side_num+1]; //�e�ӂ��f�B���N���^���E�̉��Ԗڂ��B�f�B���N���^���E��łȂ��Ȃ�side_num+1���i�[
    double *PHAT=new double [CON.get_max_DN()];//�f�B���N���^���l
	//NN���v�Z����
	set_boundary_condition3D_edge(CON,NODE,ELEM,SIDE,node,nelm,side_num,dn,&NN,PHAT,A);

	cout<<"�ިظڐ���"<<NN<<endl;
	/////////////*/
    
	    
    int pn=side_num-NN;				///���m��
    int *ppn=new int [pn];			//�s���n�Ԗڂ͕�ppn[n]
    int *npp=new int [side_num+1];	///�e�ӂ��s��̉��Ԗڂɂ����邩�B�ިظڌ^�̏ꍇ��pn+1���i�[
    int num=0; 
    for(int i=1;i<=side_num;i++)
    {
        if(SIDE[i].boundary_condition==0)//���m��
		{
			ppn[num]=i;
			npp[i]=num;
			num++;
			A[i]=0;//������
		}
		else npp[i]=pn+1;
    }
	cout<<"���m��: "<<pn<<endl;
    cout<<"�s�񕝌���@";
    ////�s��̍ő啝�v�Z
    int mat_w=0;
	int *nume=new int[side_num+1];
	int *wid=new int [pn+1];

	for(int i=1;i<=side_num;i++) nume[i]=0;
	for(int i=1;i<=nelm;i++)
	{
		for(int j=1;j<=6;j++)
		{
			int edge=ELEM[i].sides[j];
			nume[edge]=nume[edge]+1;
		}
	}
	for(int i=1;i<=side_num;i++)
	{
		int width=6+3*(nume[i]-1);
		if(width>mat_w) mat_w=width;
		if(npp[i]<pn) wid[npp[i]+1]=width;
	}

	delete [] nume;	
	//mat_w*=5;
	////////////
	
    ////�z��m��
    double **G=new double *[pn+1];///�S�̍s��
    for(int i=1;i<=pn;i++) G[i]=new double [wid[i]+1];
    int **ROW=new int *[pn+1]; ///�e�s�́A��[���v�f�̗�ԍ��L��
    for(int i=1;i<=pn;i++) ROW[i]=new int [wid[i]+1];
    int *NUM=new int [pn+1]; ///�e�s�́A��[���v�f��
    double *B=new double [pn];//���s��
    ////

	double *vB2=new double [nelm+1];//�e�v�f��dv/dB^2�̒l�i�[
	for(int je=1;je<=nelm;je++) 
	{
		vB2[je]=0;	//������
		RP[je]=1;	//������
		if(ELEM[je].material==IRON) RP[je]=2500;
		else if(ELEM[je].material==MAGELAST) RP[je]=CON.get_RP();
	}
    
	
    /////////�S�̍s����쐬����
    
    int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int table[6+1][2+1];		//�v�f���\������ӂƂ��̗v�f���ߓ_�ԍ��֌W
	cout<<"�s��쐬�J�n    ";
	int ENDFLAG=OFF;
	int loop_count=0;			//�����񐔁B������ȏ�ɂȂ�����v�Z��؂�グ��

	
	double err2;				//����`���[�v��������덷
	double old_err;
	double CG=0.01; //2e-7;				//��������(convergence)
	double CGtype=2;			//���������̃^�C�v�@1:V/A 2:A
	while(ENDFLAG==OFF)
	{
		for(int i=1;i<=pn;i++){//������
			NUM[i]=0;
			for(int j=1;j<=wid[i];j++){
				G[i][j]=0;
				ROW[i][j]=0;
			}
		}
		for(int i=0;i<pn;i++) B[i]=0;//������
		loop_count++;

		for(int je=1;je<=nelm;je++)
		{   
			//�Ӂ|�ߓ_ð��ٍ쐬
			for(int i=1;i<=6;i++)
			{
				int iside=ELEM[je].sides[i];
				int ia=SIDE[iside].node[1];
				int ib=SIDE[iside].node[2];
				for(int j=1;j<=4;j++)
				{
					if(ELEM[je].node[j]==ia) table[i][1]=j;
					else if(ELEM[je].node[j]==ib) table[i][2]=j;
				}
			}////////////////
	
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
			double Xs=0;//�v�f�̏d�S���W
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
			////////////////////////////

			///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B�������P�x���߂���Q��ڂ���͌v�Z���Ȃ��Ă���
			if(loop_count==1) ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
	
			double delta6=ELEM[je].volume;//�̐ς�6�{
		
			delta6=1/delta6;//�v�Z�ɕK�v�Ȃ̂͋t���Ȃ̂ł����Ŕ��]���Ă���
	
			double delta=ELEM[je].volume/6;//�{���̑̐�
	
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
				int m=j%4+1;
				int n=m%4+1;
	    
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);//i�������̂Ƃ�
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
				if(i%2!=0)//i����Ȃ�
				{
					b[i]*=-1;
					c[i]*=-1;
					d[i]*=-1;
					e[i]*=-1;
				}
			}
			/////////
	
			double v=v0;
			v/=RP[je]; //v�͎��C������Ȃ̂ɒ���
			double rp=RP[je];
/*			////�v�f��ظ��쐬�J�n
			for(int i=1;i<=6;i++)
			{	
				int iside=ELEM[je].sides[i];//�v�fje�̕Ӕԍ�
				if(SIDE[iside].boundary_condition==0)///�ިظڂłȂ��B�܂薢�m�Ȃ�
				{   
					int I=npp[iside]+1;///��iside�͍s���I�Ԗ�
					int I1=SIDE[iside].node[1];//iside���\������2�_
					int I2=SIDE[iside].node[2];
					int k1=table[i][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
					int k2=table[i][2];
					////���s��a���쐬
					for(int j=1;j<=6;j++)
					{
						int jside=ELEM[je].sides[j];
						int u1=table[j][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
						int u2=table[j][2];
					//	B[I-1]-=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2/3*v*delta6*delta6*delta6*A[jside];
						B[I-1]-=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2/3*delta6*delta6*delta6*A[jside]/rp;
					}
					if(ELEM[je].material==MAGNET)
					{
						B[I-1]+=1.0/18/delta*((d[k1]*e[k2]-e[k1]*d[k2])*Mx+(e[k1]*c[k2]-c[k1]*e[k2])*My+(c[k1]*d[k2]-d[k1]*c[k2])*Mz);
						//B[I-1]+=1/18/delta*((d[k1]*e[k2]-e[k1]*d[k2])*Mx+(e[k1]*c[k2]-c[k1]*e[k2])*My+(c[k1]*d[k2]-d[k1]*c[k2])*Mz);
					}
					else if(ELEM[je].material==COIL)
				{
					j0x=current[A_X][je];
					j0y=current[A_Y][je];
					j0z=current[A_Z][je];
					B[I-1]+=((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*j0x*delta6/6;
					B[I-1]+=((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*j0y*delta6/6;
					B[I-1]+=((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*j0z*delta6/6;
				}//////////
					/////////////////////////////////////
					
					for(int j=1;j<=6;j++)
					{
						int jside=ELEM[je].sides[j];
					
						int m1=table[j][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
						int m2=table[j][2];//
						
						if(SIDE[jside].boundary_condition==0)///�ިظڂłȂ��B�܂薢�m�Ȃ�
						{   
							int J=npp[jside]+1;///��jside�͍s���J�Ԗ�
							int flag=0;
						
							int J1=SIDE[jside].node[1];//jside���\������2�_
							int J2=SIDE[jside].node[2];
							for(int h=1;h<=NUM[I];h++)
							{
								if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									if(ELEM[je].material==MAGELAST)
									{
										double Uie=0;//p.106(4.88)�́Au�Ɋւ��釔�̍��v�l
										double Uje=0;//p.106(4.88)�́Ak�Ɋւ��釔�̍��v�l
										for(int u=1;u<=6;u++)
										{
											int uside=ELEM[je].sides[u];
											int u1=table[u][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
											int u2=table[u][2];
											Uie+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*A[uside];
											Uje+=(d[k1]*e[k2]-e[k1]*d[k2])*(d[m1]*e[m2]-e[m1]*d[m2])*A[uside];
										}
										double add=2/(18*18*delta*delta*delta*delta)*vB2[je]*Uje;
										G[I][h]+=add*Uie;
										//G[I][h]+=add*Uie*RP[je]*u0;
									}
					//				G[I][h]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[m1]*e[m2]-e[m1]*d[m2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[m1]*c[m2]-c[m1]*e[m2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[m1]*d[m2]-d[m1]*c[m2]))*2/3*v*delta6*delta6*delta6;//(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.106(4.88)�Ɠ���.)
									G[I][h]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[m1]*e[m2]-e[m1]*d[m2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[m1]*c[m2]-c[m1]*e[m2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[m1]*d[m2]-d[m1]*c[m2]))*2/3*delta6*delta6*delta6/rp;//(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.106(4.88)�Ɠ���.)
									flag=1;
								}
							}
							if(flag==0)
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								if(ELEM[je].material==MAGELAST)
								{
									double Uie=0;//p.106(4.88)�́Au�Ɋւ��釔�̍��v�l
									double Uje=0;//p.106(4.88)�́Ak�Ɋւ��釔�̍��v�l
									for(int u=1;u<=6;u++)
									{
										int uside=ELEM[je].sides[u];
										int u1=table[u][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
										int u2=table[u][2];
										Uie+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*A[uside];
										Uje+=(d[k1]*e[k2]-e[k1]*d[k2])*(d[m1]*e[m2]-e[m1]*d[m2])*A[uside];
									}
									double add=2/(18*18*delta*delta*delta*delta)*vB2[je]*Uje;
									G[I][H]+=add*Uie;
									//G[I][H]+=add*Uie*RP[je]*u0;
								}
						//		G[I][H]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[m1]*e[m2]-e[m1]*d[m2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[m1]*c[m2]-c[m1]*e[m2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[m1]*d[m2]-d[m1]*c[m2]))*2/3*v*delta6*delta6*delta6;//(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.106(4.88)�Ɠ���.)					    
								G[I][H]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[m1]*e[m2]-e[m1]*d[m2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[m1]*c[m2]-c[m1]*e[m2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[m1]*d[m2]-d[m1]*c[m2]))*2/3*delta6*delta6*delta6/rp;//(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.106(4.88)�Ɠ���.)
								ROW[I][H]=J;
							}
						}
						////
						else //jside���ިظڌ^���E�ߓ_�Ȃ� ���܂͌Œ�l���O�Ȃ̂Ōv�Z���Ȃ����A�����I�ɂ͂ǂ��ɂ����邱��
						{
							//int n=dn[jside];
							//B[I-1]-=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2/3*v*delta6*delta6*delta6*PHAT[n];
						}///////////
					}
				}
			}   	
		}
		///////////////////////*///�m�F
			for(int i=1;i<=6;i++)
		{	
			int iside=ELEM[je].sides[i];//�v�fje�̕Ӕԍ�
			if(SIDE[iside].boundary_condition==0)//�f�B���N���łȂ��B�܂薢�m�Ȃ�
			{   
				int I=npp[iside]+1;///��iside�͍s���I�Ԗ�
				int I1=SIDE[iside].node[1];//iside���\������2�_
				int I2=SIDE[iside].node[2];
				int k1=table[i][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
				int k2=table[i][2];

	//			if(loop_count=1){
			    for(int j=1;j<=6;j++)
				{
					int jside=ELEM[je].sides[j];
					
					int u1=table[j][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
					int u2=table[j][2];
						
					if(SIDE[jside].boundary_condition==0)//�f�B���N���łȂ��B�܂薢�m�Ȃ�
					{   
						int J=npp[jside]+1;///��jside�͍s���J�Ԗ�
						int flag=0;
						
						int J1=SIDE[jside].node[1];//jside���\������2�_
						int J2=SIDE[jside].node[2];
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
							//	G[I][h]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6*v;
								G[I][h]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6/rp;
								flag=1;
							}
						}
						if(flag==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    
						   // G[I][H]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6*v;
							G[I][H]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6/rp;
							ROW[I][H]=J;
						}
					}
					////
					else //jside���f�B���N���^���E�ߓ_�Ȃ�
					{
					    int n=dn[jside];
					//	B[I-1]-=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6*PHAT[n]*v;
						B[I-1]-=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6*PHAT[n]/rp;
					}//////////
				}

				///B[I-1]���v�Z����i���s��j
				if(ELEM[je].material==COIL)
				{
					j0x=current[A_X][je];
					j0y=current[A_Y][je];
					j0z=current[A_Z][je];
					B[I-1]+=u0*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*j0x*delta6/6;
					B[I-1]+=u0*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*j0y*delta6/6;
					B[I-1]+=u0*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*j0z*delta6/6;
				}//////////
				else if(ELEM[je].material==MAGNET)
				{
					B[I-1]+=1.0/18.0/delta*((d[k1]*e[k2]-e[k1]*d[k2])*Mx+(e[k1]*c[k2]-c[k1]*e[k2])*My+(c[k1]*d[k2]-d[k1]*c[k2])*Mz);
				}
//				}
//				else if(loop_count>1){
				////////////

	/*			////���s��a���쐬
					for(int j=1;j<=6;j++)
					{
						int jside=ELEM[je].sides[j];
						int u1=table[j][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
						int u2=table[j][2];
					//	B[I-1]-=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2/3*v*delta6*delta6*delta6*A[jside];
						B[I-1]-=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2/3*delta6*delta6*delta6*A[jside]/rp;
					}
					if(ELEM[je].material==MAGNET)
					{
						B[I-1]+=1.0/18/delta*((d[k1]*e[k2]-e[k1]*d[k2])*Mx+(e[k1]*c[k2]-c[k1]*e[k2])*My+(c[k1]*d[k2]-d[k1]*c[k2])*Mz);
						//B[I-1]+=1/18/delta*((d[k1]*e[k2]-e[k1]*d[k2])*Mx+(e[k1]*c[k2]-c[k1]*e[k2])*My+(c[k1]*d[k2]-d[k1]*c[k2])*Mz);
					}
					else if(ELEM[je].material==COIL)
				{
					j0x=current[A_X][je];
					j0y=current[A_Y][je];
					j0z=current[A_Z][je];
					B[I-1]+=((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*j0x*delta6/6;
					B[I-1]+=((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*j0y*delta6/6;
					B[I-1]+=((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*j0z*delta6/6;
				}//////////*/

				for(int j=1;j<=6;j++)
					{
						int jside=ELEM[je].sides[j];
					
						int m1=table[j][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
						int m2=table[j][2];//
						
						if(SIDE[jside].boundary_condition==0)///�ިظڂłȂ��B�܂薢�m�Ȃ�
						{   
							int J=npp[jside]+1;///��jside�͍s���J�Ԗ�
							int flag=0;
						
							int J1=SIDE[jside].node[1];//jside���\������2�_
							int J2=SIDE[jside].node[2];
							for(int h=1;h<=NUM[I];h++)
							{
								if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									if(ELEM[je].material==MAGELAST)
									{
										double Uie=0;//p.106(4.88)�́Au�Ɋւ��釔�̍��v�l
										double Uje=0;//p.106(4.88)�́Ak�Ɋւ��釔�̍��v�l
										for(int u=1;u<=6;u++)
										{
											int uside=ELEM[je].sides[u];
											int u1=table[u][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
											int u2=table[u][2];
											Uie+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*A[uside];
											Uje+=(d[k1]*e[k2]-e[k1]*d[k2])*(d[m1]*e[m2]-e[m1]*d[m2])*A[uside];
										}
										double add=2/(18*18*delta*delta*delta*delta)*vB2[je]*Uje;
										G[I][h]+=add*Uie;
										//G[I][h]+=add*Uie*RP[je]*u0;
									}
							//		G[I][h]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[m1]*e[m2]-e[m1]*d[m2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[m1]*c[m2]-c[m1]*e[m2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[m1]*d[m2]-d[m1]*c[m2]))*2/3*v*delta6*delta6*delta6;//(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.106(4.88)�Ɠ���.)
									G[I][h]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[m1]*e[m2]-e[m1]*d[m2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[m1]*c[m2]-c[m1]*e[m2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[m1]*d[m2]-d[m1]*c[m2]))*2/3*delta6*delta6*delta6/rp;//(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.106(4.88)�Ɠ���.)
									flag=1;
								}
							}
							if(flag==0)
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								if(ELEM[je].material==MAGELAST)
								{
									double Uie=0;//p.106(4.88)�́Au�Ɋւ��釔�̍��v�l
									double Uje=0;//p.106(4.88)�́Ak�Ɋւ��釔�̍��v�l
									for(int u=1;u<=6;u++)
									{
										int uside=ELEM[je].sides[u];
										int u1=table[u][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
										int u2=table[u][2];
										Uie+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*A[uside];
										Uje+=(d[k1]*e[k2]-e[k1]*d[k2])*(d[m1]*e[m2]-e[m1]*d[m2])*A[uside];
									}
									double add=2/(18*18*delta*delta*delta*delta)*vB2[je]*Uje;
									G[I][H]+=add*Uie;
									//G[I][H]+=add*Uie*RP[je]*u0;
								}
						//		G[I][H]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[m1]*e[m2]-e[m1]*d[m2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[m1]*c[m2]-c[m1]*e[m2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[m1]*d[m2]-d[m1]*c[m2]))*2/3*v*delta6*delta6*delta6;//(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.106(4.88)�Ɠ���.)					    
								G[I][H]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[m1]*e[m2]-e[m1]*d[m2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[m1]*c[m2]-c[m1]*e[m2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[m1]*d[m2]-d[m1]*c[m2]))*2/3*delta6*delta6*delta6/rp;//(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.106(4.88)�Ɠ���.)
								ROW[I][H]=J;
							}
						}
						////
						else //jside���ިظڌ^���E�ߓ_�Ȃ� ���܂͌Œ�l���O�Ȃ̂Ōv�Z���Ȃ����A�����I�ɂ͂ǂ��ɂ����邱��
						{
							//int n=dn[jside];
							//B[I-1]-=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2/3*v*delta6*delta6*delta6*PHAT[n];
						}///////////
					}//
		//		}
			}
		}   	
    }
    ///////////////////////*/

		int number=0;//�s��̔�[���v�f��
		for(int i=1;i<=pn;i++) number+=NUM[i];
    
		///�s��̎��ۂ̍ő啝�����߂�
		if(loop_count==1)
		{
			int maxN=0;
			for(int i=1;i<=pn;i++) if(NUM[i]>maxN) maxN=NUM[i];
			cout<<"�ő啝�F"<<maxN<<"/"<<mat_w<<endl;
			mat_w=maxN; //�ő啝�̏�������
		}

		///���̂܂܂ł�G[i][j]��[j]�̏��������ɕ���ł��Ȃ��̂ŁAG��ROW����ёւ���
		arrange_matrix(pn,NUM,ROW,G);

		//�Ώ̐��`�F�b�N
		check_matrix_symmetry(pn,NUM,ROW,G);
    
		ofstream fout("matrix.dat");
		for(int i=1;i<=100;i++)
		{
			for(int j=1;j<=NUM[i];j++)
			{
				fout<<G[i][j]<<" ";
				//cout<<ROW[i][j]<<" ";
			}
			fout<<endl;
			//cout<<endl;
		}
		fout.close();/////////*/

		double *val = new double [number];
		int *ind = new int [number];//��[���v�f�̗�ԍ��i�[�z��
		int *ptr = new int [pn+1];//�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
    
		/////////////////////val,ind ,ptr�ɒl���i�[
		int index=0;
		for(int n=0;n<pn;n++)
		{
			ptr[n]=index;
			for(int m=1;m<=NUM[n+1];m++)
			{
				val[index]=G[n+1][m];
				ind[index]=ROW[n+1][m]-1;
				index++;
			}
		}	    
		ptr[pn]=number;
		////////////////////*/

		///////////////////////�s��v�Z�J�n
//		if(CON.get_FEMCG()==0) ICCG3D(val,ind,ptr,pn,ppn,B,A);
//		else
//		{
			double *XX=new double [pn];//�s��̓����i�[
			if(CON.get_FEMCG()==1) ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
			else if(CON.get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
			else if(CON.get_FEMCG()==3) MRTR(CON,B,pn,XX,val,ind,ptr);
			else if(CON.get_FEMCG()==4) ICMRTR(CON,B,pn,XX,val,ind,ptr);
	
			for(int n=0;n<pn;n++)
			{
				int i=ppn[n];
				V[i]=XX[n];
			}
			delete [] XX;
//		}///////////////////////////*/
    
		//�޸�����ݼ�ٍX�V
		err2=0;//�덷��������
		if(loop_count==1)
		{
			for(int i=1;i<=side_num;i++)
			{
				if(SIDE[i].boundary_condition==0)///�ިظڂłȂ��B�܂薢�m�Ȃ�
				{
					double a=A[i];
					A[i]+=V[i];
					if(CGtype==1)
					{
						if(A[i]!=0) err2+=sqrt((V[i]/A[i])*(V[i]/A[i]));
					}
					else if(CGtype==2) err2+=sqrt((A[i]-a)*(A[i]-a));
				}
			}
		}
		else if(loop_count>1 && CGtype==1)
		{
			double *A2=new double [side_num+1];
			
			double alpha=2;//�����W��
			double err3=0;
			int breakflag=OFF;
			double tempA=alpha;
			double temperr=old_err;
			for(int n=1;n<=10;n++)
			{
				for(int i=1;i<=side_num;i++) A2[i]=A[i];
				err3=0;
				alpha*=0.5;
				for(int i=1;i<=side_num;i++)
				{
					if(SIDE[i].boundary_condition==0)///�ިظڂłȂ��B�܂薢�m�Ȃ�
					{
						A2[i]+=alpha*V[i];
						if(A2[i]!=0) err3+=sqrt((V[i]/A2[i])*(V[i]/A2[i]));
					}
				}
				if(err3<temperr)
				{
					tempA=alpha;
					temperr=err3;
				}
			}///�������Ƃ܂���
			alpha=tempA;
			cout<<loop_count<<"��ڂ�alpha="<<alpha;
			//���̃����g���Ď��ۂɂ`[i]���X�V
			for(int i=1;i<=side_num;i++)
			{
				if(SIDE[i].boundary_condition==0)///�ިظڂłȂ��B�܂薢�m�Ȃ�
				{
					A[i]+=alpha*V[i];
					if(A[i]!=0) err2+=sqrt((V[i]/A[i])*(V[i]/A[i]));
				}
			}
			delete [] A2;
		}//////////*/
		else if(loop_count>1 && CGtype==2)
		{
			for(int i=1;i<=side_num;i++)
			{
				if(SIDE[i].boundary_condition==0)///�ިظڂłȂ��B�܂薢�m�Ȃ�
				{
					double a=A[i];
					A[i]+=V[i];
				//	err2+=sqrt((A[i]-a)*(A[i]-a));
					err2+=sqrt(V[i]*V[i]);
				}
			}
		}//////////*/
		
		
		if(err2<CG || loop_count>5) ENDFLAG=ON;
		cout<<"  �덷="<<err2<<endl;
		old_err=err2;//�O���err��ۑ�
		//if(CGtype==2 && loop_count==1) CG=err2*1e-6;//�덷����Ƃ���A[i]�݂̂��g���Ƃ��́A���قɂ��œK�Ȍ덷���肪�ς��̂ŁA���̂悤�ɒ�`����

		////////////////////////////�V���������������Ƃ߂邽�߂́A�����̗v�f�̎��������Ƃ߂�
		if(ENDFLAG==OFF)
		{
			for(int je=1;je<=nelm;je++)
			{   
				if(ELEM[je].material==MAGELAST)
				{
					//�Ӂ|�ߓ_ð��ٍ쐬
					for(int i=1;i<=6;i++)
					{
						int iside=ELEM[je].sides[i];
						int ia=SIDE[iside].node[1];
						int ib=SIDE[iside].node[2];
						for(int j=1;j<=4;j++)
						{
							if(ELEM[je].node[j]==ia) table[i][1]=j;
							else if(ELEM[je].node[j]==ib) table[i][2]=j;
						}
					}////////////////
			
					for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
					
					double Xs=0;//�v�f�̏d�S���W
					double Ys=0;
					double Zs=0;
					for(int j=1;j<=4;j++)
					{
						X[j]=NODE[N[j]].r[A_X];
						Y[j]=NODE[N[j]].r[A_Y];
						Z[j]=NODE[N[j]].r[A_Z];
						Xs+=X[j];
						Ys+=Y[j];
						Zs+=Z[j];
					}
					Xs/=4;Ys/=4;Zs/=4;
					////////////////////////////
		
					double delta6=ELEM[je].volume;//�̐ς�6�{(�������̐ς̒l�͂��łɃx�N�g���|�e���V���������߂�ۂɌv�Z���Ă���)
			
					delta6=1/delta6;//�v�Z�ɕK�v�Ȃ̂͋t���Ȃ̂ł����Ŕ��]���Ă���
			
					double delta=ELEM[je].volume/6;//�{���̑̐�
			
					for(int i=1;i<=4;i++)
					{
						int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
						int m=j%4+1;
						int n=m%4+1;
			    
						c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
						d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
						e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
						if(i%2!=0)//i����Ȃ�
						{
						    c[i]*=-1;
							d[i]*=-1;
							e[i]*=-1;
						}
					}
					/////////
		
					for(int D=0;D<3;D++) Be[D][je]=0;//������
		
					for(int i=1;i<=6;i++)
					{
						int s=ELEM[je].sides[i];
						int k1=table[i][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
						int k2=table[i][2];
			
						Be[A_X][je]+=(d[k1]*e[k2]-e[k1]*d[k2])*A[s];
						Be[A_Y][je]+=(e[k1]*c[k2]-c[k1]*e[k2])*A[s];
						Be[A_Z][je]+=(c[k1]*d[k2]-d[k1]*c[k2])*A[s];
					}
			
					for(int D=0;D<3;D++) Be[D][je]*=delta6*delta6*2;
			
					double BB=0;//�������x�̂Q��
					for(int D=0;D<3;D++) BB+=Be[D][je]*Be[D][je];
//					cout<<"je="<<je<<", B="<<sqrt(BB)<<endl;
					/////�A�L�}��Ԃɂ��v-B^2�Ȑ�����v�f��RP[i]�����Ƃ߂�
					int flagB=0;
					for(int n=2;n<=Nn+1;n++)
					{
						if(BB>=Bflux2[n] && BB<Bflux2[n+1])
						{
							double h=Bflux2[n+1]-Bflux2[n];
							double aj0=vm[n];
							double aj1=dvdB2[n];
							double aj2=(3*mm3[n]-2*dvdB2[n]-dvdB2[n+1])/h;
							double aj3=(dvdB2[n]+dvdB2[n+1]-2*mm3[n])/(h*h);
							double val=aj0+aj1*(BB-Bflux2[n])+aj2*(BB-Bflux2[n])*(BB-Bflux2[n])+aj3*(BB-Bflux2[n])*(BB-Bflux2[n])*(BB-Bflux2[n]);
							RP[je]=1/val/u0;
//							cout<<"je="<<je<<",RP="<<RP[je]<<"aj0="<<aj0<<"aj1="<<aj1<<"aj2="<<aj2<<"aj3="<<aj3<<"val="<<val<<endl;
							double bj0=dvdB2[n];
							double bj1=t4[n];
							double bj2=(3*mm4[n]-2*t4[n]-t4[n+1])/h;
							double bj3=(t4[n]+t4[n+1]-2*mm4[n])/(h*h);
							double val2=bj0+bj1*(BB-Bflux2[n])+bj2*(BB-Bflux2[n])*(BB-Bflux2[n])+bj3*(BB-Bflux2[n])*(BB-Bflux2[n])*(BB-Bflux2[n]);
							vB2[je]=val2;
							flagB=1;
						}
						
					}
					if(flagB==0) cout<<"�������x���z��̈�O�ł� B="<<sqrt(BB)<<endl;
					
				}
			}
		}
		delete [] val;
		delete [] ind;
		delete [] ptr;
	}
	///�{���v���O������V[i]��ߓ_i���޸�����ݼ�قƂ��Čv�Z����̂ŁA������V[i]=A[i]�Ƃ��Ēl��n���Ă���
	for(int i=1;i<=side_num;i++) V[i]=A[i];
	

    ///////////
	ofstream fp2("V.dat");
    for(int i=1;i<=node;i++)
    {
		if(NODE[i].boundary_condition==1) fp2<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<endl;
    }
    fp2.close();
	ofstream fp4("RP.dat");
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==MAGELAST)
		{
			double R[3]={0,0,0};
			for(int j=1;j<=4;j++) for(int D=0;D<3;D++) R[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
		
			if(R[1]<0.0015 && R[1]>0) 
			{
				fp4<<R[0]<<" "<<R[2]<<" "<<RP[i]<<endl;
			}
		}
	}
	fp4.close();

    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
	delete [] wid;
    ///////////////////////*/
    delete [] dn;
    delete [] PHAT;
    delete [] npp;
    delete [] ppn;
    delete [] B;
	delete [] A;
    
    

	delete [] vB2;

	delete [] mm;
	delete [] vm;
	delete [] mm2;
	delete [] t2;
	delete [] Bflux2;
	delete [] mm3;
	delete [] dvdB2;
	delete [] mm4;
	delete [] t4;
}
//ICMRTR�@//
void ICMRTR(mpsconfig &CON,double *B,int pn,double *X,double *val,int *ind,int *ptr)
{
	////////////////////////////////////////////////////////////////
	//val :�[���v�f�̒l
	//ind:��[���v�f�̗�ԍ��i�[�z��
	//ptr:�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
	//X[n]:��
	int count=0;
	double E=1;//�덷
	double BB=0;
	double rr=0;
		
	double alp;
	double beta;
	double zeta;
	double zeta_old;
	double eta;
	double nu;

	double Ar_r;//cAr_r: (cAr,r)(����)���w���B(u,v)=��((cu)*v)
	double Ar_Ar;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
	double y_Ar;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
	double Ar_y;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)

	double w_r;
	double v_w;
	double y_w;
	double w_y;

	double *r= new double[pn];
	double *Ar = new double [pn];  //   _ 
	double *P = new double [pn];
	double *y = new double [pn];
	double *u = new double [pn];
	double *v = new double [pn];
	double *w = new double [pn];
	double *z = new double [pn];

	zeta=0;
	zeta_old=0;
	eta=0;
	nu=1;

	ofstream c("convergence.dat");

	for(int n=0;n<pn;n++)//�����l
	{
		X[n]=0;
		r[n]=B[n];
		P[n]=r[n];
		y[n]=0;
	}

	//////////�O����//////////
	double accel_re=CON.get_CGaccl();
	double accel;	//�����t�@�N�^
	accel= accel_re;

	int num2=0;//�Ίp�������܂ށA���O�p�s�񂾂����l���ɂ��ꂽ��[���v�f��
	for(int k=0;k<pn;k++) for(int m=ptr[k];m<ptr[k+1];m++) if(ind[m]<=k) num2++;

	double *val2=new double [num2];//�W���s���ۑ�(�Ίp�������܂މ��O�p�s��) ��[���v�f������1���z��Ƃ��ĕۑ�
	int *ind2 = new int [num2];
	int *ptr2 = new int [pn+1];

	num2=0;
	double one=1;
	double sum=0;
	for(int k=0;k<pn;k++)
	{	
		ptr2[k]=num2;
	    for(int m=ptr[k];m<ptr[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
			if(ind[m]<=k)
			{
				val2[num2]=val[m];
				ind2[num2]=ind[m];
				if(ind[m]==k) val2[num2]*=accel;
				num2++;
			}
		}
	}
	ptr2[pn]=num2;//��������Ă����Ȃ��ƁA�Ō��(int m=ptr2[k];m<ptr2[k+1];m++)�݂����Ȃ��Ƃ��ł��Ȃ�

	int *NUM = new int [pn];			//������ɂ݂��A�e��̗v�f��
	for(int k=0;k<pn;k++) NUM[k]=0;
	
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
			int J=ind2[m];
			NUM[J]=NUM[J]+1;
		}
	}
	double **VAL=new double *[pn];						//�[���v�f�̒l VAL[i][k]��i���k�Ԗڂ̔�[���v�f
	for(int i=0;i<pn;i++) VAL[i]=new double [NUM[i]];
	int **IND = new int *[pn];							//��[���v�f�̍s�ԍ��i�[�z��
	for(int i=0;i<pn;i++) IND[i]=new int [NUM[i]];
	
	/////////////////////////////////////iccg�@
	double rLDLt_r;
	double *y2= new double [pn];
	double *LDLt_r= new double [pn];
//	double *L=new double[num2];//�s���S�R���X�L�[������̉��O�p�s��L�i�[
	double *D1 = new double [pn];//D�s��
	/////�s���S�R���X�L�[����
//	Incomplete_Cholesky_Decomposition(CON,L,D1,val2,ptr2,ind2,pn,num2);//L��D1�ɒl���������܂��
	for(int k=0;k<pn;k++)
	{
		for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�0�v�f
		{
			int i=ind2[m];//��ԍ�
			if(i==0)
			{
				val2[m]=val2[m];
			}
			if(i>0 && i<k)
			{
				sum=0;
				for(int j=ptr2[k];j<m;j++)
				{
					for(int J=ptr2[i];J<ptr2[i+1];J++)
					{
						if(ind2[J]==ind2[j]) sum+=val2[j]*val2[J]*D1[ind2[j]];
					}
				}
				val2[m]=val2[m]-sum;
			}
			if(i==k)
			{
				sum=0;
				for(int j=ptr2[k];j<m;j++) sum+=val2[j]*val2[j]*D1[ind2[j]];
				val2[m]=val2[m]-sum;
				D1[k]=one/val2[m];
			}
		}
	}
	//�s���S�R���X�L�[��������///////
//	delete [] val2;
	///�����ɂ����z��ɒl����
	for(int k=0;k<pn;k++) NUM[k]=0;
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
			int J=ind2[m];
			VAL[J][NUM[J]]=val2[m];
			IND[J][NUM[J]]=k;
			NUM[J]=NUM[J]+1;
		}
	}////////*/
	
	/////////////////y[i]�����Ƃ߁ALDLt_r[i]�����Ƃ߂�B
	for(int i=0;i<pn;i++)
	{
		if(i==0) y2[0]=r[0]/val2[0]; //���i3.77�j 
		else
		{
		    double sum=0;
			for(int m=ptr2[i];m<ptr2[i+1]-1;m++) sum+=val2[m]*y2[ind2[m]];//���i3.78�j
		    int m=ptr2[i+1]-1;
			y2[i]=(r[i]-sum)/val2[m];
		}
	}////y[i]�����Ƃ܂����B
	for(int i=pn-1;i>=0;i--)
	{
	    double sum=0;
		for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
	    LDLt_r[i]=y2[i]-D1[i]*sum;	
	}
	/////////////////*/

	for(int n=0;n<pn;n++) u[n]=LDLt_r[n];
	////////////////////////////////////////////////////////////////
	for(int n=0;n<pn;n++) BB+=abs(B[n])*abs(B[n]);
	ofstream ICEEe("ICMRTR.dat", ios::trunc);
	ICEEe.close();
	cout<<"ICMRTR�@�X�^�[�g-----���m��="<<pn<<"  ---";
	unsigned int time=GetTickCount();
	if(CON.get_omp_P()==OFF)//�ʏ��
	{
		//while(E>CON->get_FEMCGep())// EP=CON->get_CGep();//��������(convergence test)
		while(E>CON.get_MRTRep())// CON.get_MRTRep()//��������(convergence test)
		{
			if(count==pn) cout<<"count=pn"<<endl;
			ofstream ICEE("ICMRTR.dat", ios::app);
			ICEE<<count<<" "<<E<<endl;
			ICEE.close();
			count++;
			/////v�̌v�Z
			for(int n=0;n<pn;n++)
			{
				v[n]=0;
				for(int j=ptr[n];j<ptr[n+1];j++)
				{
					v[n]+=val[j]*u[ind[j]];
				}
			}
		
			/////Ar,Ar_c�v�Z
			for(int n=0;n<pn;n++)
			{    
				Ar[n]=0.0;
				for(int j=ptr[n];j<ptr[n+1];j++)
				{
					Ar[n]+=val[j]*u[ind[j]];
				}
			}
		
			//////w�̌v�Z
			//////////////////////////y[i]�����߁ALDLt_r[i]�����߂�B
			for(int i=0;i<pn;i++)
			{
				if(i==0) y2[0]=v[0]/val2[0]; //��(3.77)
				else
				{
					sum=0;
					for(int m=ptr2[i];m<ptr2[i+1]-1;m++) sum+=val2[m]*y2[ind2[m]];//���i3.78�j
				    int m=ptr2[i+1]-1;
						y2[i]=(v[i]-sum)/val2[m];
				}
			}////y[i]�����Ƃ܂����B
		
			for(int i=pn-1;i>=0;i--)
			{
				double sum=0;
				for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
				LDLt_r[i]=y2[i]-D1[i]*sum;	
			}
			/////////////////*/
			for(int n=0;n<pn;n++) w[n]=LDLt_r[n];

			w_r=0;
			v_w=0;
			y_w=0;
			w_y=0;
		
			//���ς̌v�Z
			for(int n=0;n<pn;n++)
			{
				w_r+=w[n]*r[n];
				v_w+=v[n]*w[n];
				y_w+=y[n]*w[n];
				w_y+=w[n]*y[n];
			}	

			zeta=nu*w_r/(nu*v_w-y_w*w_y);
			eta=-y_w*w_r/(nu*v_w-y_w*w_y);
			//////nu�̌v�Z//////////////
			nu=0.0;
			for(int n=0;n<pn;n++) nu+=zeta*r[n]*w[n];
			
			//////p�̌v�Z//////////////
			for(int n=0;n<pn;n++) P[n]=r[n] + (zeta_old*eta/zeta)*P[n];
			zeta_old=zeta;

			///////X�̌v�Z////////////
			for(int n=0;n<pn;n++) X[n]+=zeta*P[n];
	
			//////y�̌v�Z////////////
			for(int n=0;n<pn;n++) y[n]=eta*y[n]+zeta*Ar[n];
	
			//////r�̌v�Z////////////
			for(int n=0;n<pn;n++) r[n]-=y[n];

			//////z�̌v�Z///////////
			for(int n=0;n<pn;n++) z[n]=eta*z[n]+zeta*w[n];

			//////u�̌v�Z///////////
			for(int n=0;n<pn;n++) u[n]-=z[n];

			/////////////////*/
			//�덷�]��
			rr=0;
			for(int n=0;n<pn;n++) rr+=abs(r[n])*abs(r[n]);
			
			E=sqrt(rr/BB);
			c<<count<<" "<<E<<endl;
			if(count==30000) break;
	
			if(count>pn)
			{
				count=1;		//�����񐔂��P�ɂ��Ċ֐�����E�o�B�����񐔂��P������A�G���[���Ƃ킩��
				E=0;
			}
		}
	}
	cout<<"�I�� ������="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;
	
	delete [] LDLt_r;
	delete [] D1;

	delete [] val2;
	delete [] ind2;
	delete [] ptr2;

	for(int i=0;i<pn;i++)
	{
		delete [] VAL[i];
		delete [] IND[i];
	}
	delete [] VAL;
	delete [] IND;
	delete [] NUM;
	delete [] y2;

	delete [] r;
	delete [] Ar;
	delete [] P;
	delete [] y;
	delete [] u;
	delete [] v;
	delete [] w;
	delete [] z;
	
	/////////////////
	c.close();
}
//MRTR�@//
void MRTR(mpsconfig &CON,double *B,int pn,double *X,double *val,int *ind,int *ptr)
{
	
	//unsigned int timeCG=GetTickCount();
	int count=0;
	double EP=CON.get_MRTRep();
	//complex<double> rr=(0,0);
	double E=1;//�덷
	double BB=0;
	double rr=0;
		
	double zeta;
	double zeta_old;
	double eta;
	double v;

	double Ar_r;//cAr_r: (cAr,r)(����)���w���B(u,v)=��((cu)*v)
	double Ar_Ar;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
	double y_Ar;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
	double Ar_y;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)

    double *r= new double[pn];
	double *Ar = new double [pn];  //   _ _
	
	double *P = new double [pn];
	double *y = new double [pn];

	zeta=0;
	zeta_old=0;
	eta=0.0;
	//v=complex<double>(1.0,0.0);//1��ڂ̔����ɂ����āA0�łȂ����Ƃɒ���
	v=1.0;

	ofstream c("convergence.dat");

	for(int n=0;n<pn;n++) //�����l
	{
		X[n]=0.0;			//�ߎ���
		r[n]=B[n];			//�E�Ӎs��@Ax=b
		P[n]=r[n];			//r�͎c��
		y[n]=0.0;
	}

	for(int n=0;n<pn;n++) BB+=r[n]*r[n];
	//BB=sqrt(BB);
	//cout<<"||r0||2="<<sqrt(BB)<<endl;
	//cout<<"r[0]="<<r[0]<<" r[1]="<<r[1]<<endl;
	ofstream EEe("MRTR.dat", ios::trunc);
	EEe.close();
	cout<<"MRTR�@�X�^�[�g-----���m��="<<pn<<"  ---";
	unsigned int time=GetTickCount();
	if(CON.get_omp_P()==OFF)//�ʏ��
	{
		//while(E>CON->get_FEMCGep())// EP=CON->get_CGep();//��������(convergence test)
		while(E>EP)// EP=CON->get_CGep();//��������(convergence test)
		{
			//if(count==pn) cout<<"count=pn"<<endl;
			count++;
			ofstream EE("MRTR.dat", ios::app);
			EE<<count<<" "<<E<<endl;
			EE.close();
			for(int n=0;n<pn;n++)//Ar,Ar_c�v�Z
			{    
				Ar[n]=0.0;
				for(int j=ptr[n];j<ptr[n+1];j++)
				{
					Ar[n]+=val[j]*r[ind[j]];
				}
			}

			Ar_r=0.0;//cAr_r: (cAr,r)(����)���w���B(u,v)=��((cu)*v)
			Ar_Ar=0.0;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
			y_Ar=0.0;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
			Ar_y=0.0;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
	
			//#pragma omp parallel for
			for(int n=0;n<pn;n++)  
			{
				Ar_r+=Ar[n]*r[n];
				Ar_Ar+=Ar[n]*Ar[n];
				y_Ar+=y[n]*Ar[n];
				Ar_y+=Ar[n]*y[n];
			}

			zeta=v*Ar_r/(v*Ar_Ar-y_Ar*Ar_y);
			eta=-y_Ar*Ar_r/(v*Ar_Ar-y_Ar*Ar_y);
		//	cout<<zeta<<" "<<eta<<endl;

			//////v�̌v�Z//////////////
			v=0.0;
			for(int n=0;n<pn;n++) v+=zeta*Ar[n]*r[n];

			//////p�̌v�Z//////////////
			for(int n=0;n<pn;n++) P[n]=r[n] + (zeta_old*eta/zeta)*P[n];
			zeta_old=zeta;

			///////X�̌v�Z////////////
			for(int n=0;n<pn;n++) X[n]+=zeta*P[n];
	
			//////y�̌v�Z////////////
			for(int n=0;n<pn;n++)
			{
				y[n]=eta*y[n]+zeta*Ar[n];
				//cout<<y[n]<<endl;
			}

			//////r�̌v�Z////////////
			for(int n=0;n<pn;n++) r[n]-=y[n];

			//�덷�]��
			rr=0;
			for(int n=0;n<pn;n++) rr+=r[n]*r[n];
			//for(int n=0;n<pn;n++) rr+=conj(r[n])*conj(r[n]);
			
			E=sqrt(rr/BB);
			//cout<<E<<endl;
			c<<count<<" "<<E<<endl;
			//if(count%1000==0) cout<<"������="<<count<<" E="<<E<<endl;
	
			if(count>pn)
			{
				count=1;		//�����񐔂��P�ɂ��Ċ֐�����E�o�B�����񐔂��P������A�G���[���Ƃ킩��
				E=0;
			}
			if(count>=10000) break;
		}
	}
	cout<<"�I�� ������="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;
	
	delete [] r;
	delete [] Ar;
	
	delete [] P;
	delete [] y;
	c.close();
}
////ic�t��MRTR�@
//void cs_ICMRTR(mpsconfig *CON,complex<double> *val,int *ind,int *ptr,int pn,complex<double> *B,int number,complex<double> *X)
//{
//	
//	//unsigned int timeCG=GetTickCount();
//	int count=0;
//	//complex<double> rr=(0,0);
//	double E=1;//�덷
//	double BB=0;
//	double rr=0;
//	
//	complex<double> alp;
//	complex<double> beta;
//	complex<double> zeta;
//	complex<double> zeta_old;
//	complex<double> eta;
//	complex<double> nu;
//
//	complex<double> cAr_r;//cAr_r: (cAr,r)(����)���w���B(u,v)=��((cu)*v)
//	complex<double> cAr_Ar;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
//	complex<double> cy_Ar;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
//	complex<double> cAr_y;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
//
//	complex<double> cw_r;//cAr_r: (cAr,r)(����)���w���B(u,v)=��((cu)*v)
//	complex<double> cv_w;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
//	complex<double> cy_w;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
//	complex<double> cw_y;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
//
//    complex<double> *r= new complex<double>[pn];
//	complex<double> *Ar = new complex<double> [pn];  //   _ _
//	complex<double> *cAr = new complex<double> [pn];//cAr:A r(���f�������o�[�ŕ\�L)
//	complex<double> *P = new complex<double> [pn];
//	complex<double> *y = new complex<double> [pn];
//	complex<double> *u = new complex<double> [pn];
//	complex<double> *v = new complex<double> [pn];
//	complex<double> *w = new complex<double> [pn];
//	complex<double> *z = new complex<double> [pn];
//
//	zeta=complex<double>(0.0,0.0);
//	zeta_old=complex<double>(0.0,0.0);
//	eta=complex<double>(0.0,0.0);
//	nu=complex<double>(1.0,0.0);//1��ڂ̔����ɂ����āA0�łȂ����Ƃɒ���
//
//	ofstream c("convergence.dat");
//
//	for(int n=0;n<pn;n++) //�����l
//	{
//		X[n]=complex<double> (0.0,0.0);
//		r[n]=B[n];
//		P[n]=r[n];
//		y[n]=complex<double> (0.0,0.0);
//		u[n]=complex<double> (0.0,0.0);
//		v[n]=complex<double> (0.0,0.0);
//		w[n]=complex<double> (0.0,0.0);
//		z[n]=complex<double> (0.0,0.0);
//	}
//
//	//////�O����/////
//	double accel_re=CON->get_CGaccl();
//	complex<double> accel;//CON->get_CGaccl();//�����t�@�N�^
//	accel=complex<double> (accel_re,0);
//	double accel2=0.001;//���f�V�t�g�p�����t�@�N�^
//	
//	int num2=0;//�Ίp�������܂ށA���O�p�s�񂾂����l���ɂ��ꂽ��[���v�f��
//	for(int k=0;k<pn;k++) for(int m=ptr[k];m<ptr[k+1];m++) if(ind[m]<=k) num2++;	
//	if(num2!=(number-pn)/2+pn) cout<<"ERROR"<<endl;
//	
//	complex<double> *val2=new complex<double> [num2];
//	int *ind2 = new int [num2];
//	int *ptr2 = new int [pn+1];
//
//	num2=0;
//
//	complex<double> one;
//	one=complex<double> (1.0,0);
//	complex<double> Im;
//	Im=complex<double> (0.0,1.0);
//	complex<double> sum;
//	sum=complex<double> (0.0,0.0);
//
//	for(int k=0;k<pn;k++)
//	{	
//		ptr2[k]=num2;
//	    for(int m=ptr[k];m<ptr[k+1];m++)///k�s�ڂ̔�O�v�f
//	    {
//			if(ind[m]<=k)
//			{
//				val2[num2]=val[m];
//				ind2[num2]=ind[m];
//				if(ind[m]==k) val2[num2]*=accel;//����̧��
//				//if(ind[m]==k) val2[num2]+=accel2*Im;//���f�V�t�g
//				num2++;
//			}
//		}
//	}
//	ptr2[pn]=num2;//��������Ă����Ȃ��ƁA�Ō��(int m=ptr2[k];m<ptr2[k+1];m++)�݂����Ȃ��Ƃ��ł��Ȃ�
//
//	int *NUM = new int [pn];
//	for(int k=0;k<pn;k++) NUM[k]=0;
//	
//	for(int k=0;k<pn;k++)
//	{	
//	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�O�v�f
//	    {
//			int J=ind2[m];
//			NUM[J]=NUM[J]+1;
//		}
//	}
//	complex<double> **VAL=new complex<double> *[pn];//�[���v�f�̒l
//	for(int i=0;i<pn;i++) VAL[i]=new complex<double> [NUM[i]];
//	int **IND = new int *[pn];//��[���v�f�̍s�ԍ��i�[�z��
//	for(int i=0;i<pn;i++) IND[i]=new int [NUM[i]];
//	
//	/////////////////////////////////////iccg�@
//	complex<double> rLDLt_r;
//	complex<double> *y2=new complex<double> [pn];
//	complex<double> *LDLt_r= new complex<double> [pn];
//	complex<double> *D1 = new complex<double> [pn];//D�s��
//	
//	/////�s���S�R���X�L�����
//	for(int k=0;k<pn;k++)
//	{	
//	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�O�v�f
//	    {
//	        int i=ind2[m];//��ԍ�
//	        if(i==0)
//			{
//				val2[m]=val2[m];
//				//if(val2[m]<0.0001 &&val2[m]>=0) val2[m]=0.0001;
//				//else if(val2[m]>-0.0001 &&val2[m]<=0) val2[m]=-0.0001;
//		    
//			}
//			if(i>0 && i<k)
//			{
//				sum=complex<double> (0,0);
//				
//				for(int j=ptr2[k];j<m;j++)
//				{	
//					for(int J=ptr2[i];J<ptr2[i+1];J++)//i�s�ڂ̂Ȃ������̈�v������̂�T���Ă���B������Ԃ��H
//					{
//						if(ind2[J]==ind2[j]) sum+=val2[j]*val2[J]*D1[ind2[j]];
//					}
//				}
//				val2[m]=val2[m]-sum;
//			}
//			if(i==k)
//			{
//				sum=complex<double> (0,0);
//				for(int j=ptr2[k];j<m;j++) sum+=val2[j]*val2[j]*D1[ind2[j]];
//				val2[m]=val2[m]-sum;
//				//if(val2[m]<0.0001 &&val2[m]>=0) val2[m]=0.0001;
//				//else if(val2[m]>-0.0001 &&val2[m]<=0) val2[m]=-0.0001;
//				//if(val2[m].real()<0.0001 &&val2[m].real()>=0) val2[m].real(0.0001);
//				//else if(val2[m].real()>-0.0001 &&val2[m].real()<=0) val2[m].real(-0.0001);
//				//if(val2[m].imag()<0.0001 &&val2[m].imag()>=0) val2[m].imag(0.0001);
//				//else if(val2[m].imag()>-0.0001 &&val2[m].imag()<=0) val2[m].imag(-0.0001);
//				
//				D1[k]=one/val2[m];
//				//if(val2[m]>0) cout<<"EE"<<endl;
//            }
//	    }
//	}    
//	///�s���S�ڽ����������/////////*/
//
//	///�����ɂ����z��ɒl����
//	for(int k=0;k<pn;k++) NUM[k]=0;
//	for(int k=0;k<pn;k++)
//	{	
//	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�O�v�f
//	    {
//			int J=ind2[m];
//			VAL[J][NUM[J]]=val2[m];
//			IND[J][NUM[J]]=k;
//			NUM[J]=NUM[J]+1;
//		}
//	}////////*/
//
//	/////////////////y[i]�����Ƃ߁ALDLt_r[i]�����Ƃ߂�B
//	for(int i=0;i<pn;i++)
//	{
//		if(i==0) y2[0]=r[0]/val2[0]; //���i3.77�j 
//		else
//		{
//		    sum=complex<double> (0,0);
//		    /////////        
//		    for(int m=ptr2[i];m<ptr2[i+1]-1;m++) sum+=val2[m]*y2[ind2[m]];//���i3.78�j
//		    int m=ptr2[i+1]-1;
//		    y2[i]=(r[i]-sum)/val2[m];
//		}
//	}////y[i]�����Ƃ܂����B
//	for(int i=pn-1;i>=0;i--)
//	{
//	    sum=complex<double> (0,0);
//		for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
//	    LDLt_r[i]=y2[i]-D1[i]*sum;	
//	}
//	/////////////////*///LDLt��M-1�ɑ����H�܂�ALDLt_r��u0�Ɠ�����
//	
//	for(int n=0;n<pn;n++) u[n]=LDLt_r[n];
//		/////////////////*/
//	//////////////////////
//
//	//for(int n=0;n<pn;n++) BB+=norm(B[n]);
//	for(int n=0;n<pn;n++) BB+=abs(B[n])*abs(B[n]);
//	//BB=sqrt(BB);
//	cout<<"||r0||2="<<sqrt(BB)<<endl;
//
//
//	 cout<<"cs_ICMRTR�@�X�^�[�g-----���m��="<<pn<<"  ---";
//	unsigned int time=GetTickCount();
//	if(CON->get_omp_P()==OFF)//�ʏ��
//	{
//		//while(E>CON->get_FEMCGep())// EP=CON->get_CGep();//��������(convergence test)
//		while(E>CON->get_FEMCGep())// EP=CON->get_CGep();//��������(convergence test)
//		{
//			if(count==pn) cout<<"count=pn"<<endl;
//			count++;
//
//			/////v�̌v�Z
//			for(int n=0;n<pn;n++)
//			{    
//				v[n]=complex<double>(0.0,0.0);//v�ɑ���
//				for(int j=ptr[n];j<ptr[n+1];j++)
//				{
//					v[n]+=val[j]*u[ind[j]];
//				}
//			}
//
//			/*///
//			for(int n=0;n<pn;n++)//Ar,Ar_c�v�Z
//			{    
//				Ar[n]=complex<double>(0.0,0.0);//v�ɑ���
//				cAr[n]=complex<double>(0.0,0.0);//conj(v)�ɑ���
//				for(int j=ptr[n];j<ptr[n+1];j++)
//				{
//					Ar[n]+=val[j]*u[ind[j]];
//					cAr[n]+=conj(val[j]) * conj(u[ind[j]]);
//				}
//			}
//			///*/
//
//			////w�̌v�Z
//
//			/////////////////y[i]�����Ƃ߁ALDLt_r[i]�����Ƃ߂�B
//			for(int i=0;i<pn;i++)
//			{
//				if(i==0) y2[0]=v[0]/val2[0]; //���i3.77�j 
//				else
//				{
//					sum=complex<double> (0,0);
//					/////////        
//					for(int m=ptr2[i];m<ptr2[i+1]-1;m++) sum+=val2[m]*y2[ind2[m]];//���i3.78�j
//					int m=ptr2[i+1]-1;
//					y2[i]=(v[i]-sum)/val2[m];
//				}
//			}////y[i]�����Ƃ܂����B
//			for(int i=pn-1;i>=0;i--)
//			{
//				sum=complex<double> (0,0);
//				for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
//				LDLt_r[i]=y2[i]-D1[i]*sum;	
//			}
//			/////////////////*///LDLt��M-1�ɑ����H�܂�ALDLt_r��u0�Ɠ����� ���l�ɁALDL_v��w�H
//	
//			for(int n=0;n<pn;n++) w[n]=LDLt_r[n];
//
//			cw_r=complex<double>(0.0,0.0);//cAr_r: (cAr,r)(����)���w���B(u,v)=��((cu)*v)
//			cv_w=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
//			cy_w=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
//			cw_y=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
//	
//			//���ς̌v�Z
//			for(int n=0;n<pn;n++)  
//			{
//				cw_r+=w[n]*r[n];
//				cv_w+=v[n]*w[n];
//				cy_w+=y[n]*w[n];
//				cw_y+=w[n]*y[n];
//			}
//
//			zeta=nu*cw_r/(nu*cv_w-cy_w*cw_y);
//			eta=-cy_w*cw_r/(nu*cv_w-cy_w*cw_y);
//
//			//////nu�̌v�Z//////////////
//			nu=complex<double>(0.0,0.0);
//			for(int n=0;n<pn;n++) nu+=zeta*r[n]*w[n];
//
//			//////p�̌v�Z//////////////
//			for(int n=0;n<pn;n++) P[n]=u[n] + (zeta_old*eta/zeta)*P[n];
//			zeta_old=zeta;
//
//			///////X�̌v�Z////////////
//			for(int n=0;n<pn;n++) X[n]+=zeta*P[n];
//	
//			//////y�̌v�Z////////////
//			for(int n=0;n<pn;n++) y[n]=eta*y[n]+zeta*v[n];
//
//			//////r�̌v�Z////////////
//			for(int n=0;n<pn;n++) r[n]-=y[n];
//
//			//////z�̌v�Z////////////
//			for(int n=0;n<pn;n++) z[n]=eta*z[n]+zeta*w[n];
//
//			//////u�̌v�Z////////////
//			for(int n=0;n<pn;n++) u[n]-=z[n];
//
//			//�덷�]��
//			rr=0;
//			//for(int n=0;n<pn;n++) rr+=norm(r[n]);
//			for(int n=0;n<pn;n++) rr+=abs(r[n])*abs(r[n]);
//			//for(int n=0;n<pn;n++) rr+=conj(r[n])*conj(r[n]);
//			
//			E=sqrt(rr/BB);
//			c<<count<<" "<<E<<endl;
//			if(count%1000==0) cout<<"������="<<count<<" E="<<E<<endl;
//	
//			
//		}
//	}
//
//	cout<<"�I�� ������="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;
//	
//	delete [] LDLt_r;
//	delete [] D1;
//
//	delete [] val2;
//	delete [] ind2;
//	delete [] ptr2;
//
//	for(int i=0;i<pn;i++)
//	{
//		delete [] VAL[i];
//		delete [] IND[i];
//	}
//	delete [] VAL;
//	delete [] IND;
//	delete [] NUM;
//	delete [] y2;
//
//	delete [] r;
//	delete [] Ar;
//	delete [] cAr;	
//	delete [] P;
//	delete [] y;
//	delete [] u;
//	delete [] v;
//	delete [] w;
//	delete [] z;
//
//	c.close();
//}

