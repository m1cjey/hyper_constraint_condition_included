#include "stdafx.h"	
//���b�V���������~�܂�ꍇ���炩�̗�O�����@�\���K�v�E�E�E
//�����炭poly3D_for_FINE3D()�Ŏ~�܂��Ă���

void box3D(vector<point3D> &NODE,vector<element3D> &ELEM,double alpha,int *nelm,double rax,double ray,double raz,int KTJ);
double volume3D(vector<point3D> &NODE,int ia,int ib,int ic,int ip);
void sphere3D(vector<point3D> &NODE,vector<element3D> &ELEM,int ia,int ib,int ic,int ip,int i);
int locate3D(vector<point3D> &NODE,vector<element3D> &ELEM,int nelm,double xp,double yp,double zp);
int locate3D(vector<point3D> &NODE,vector<element3D> &ELEM,int nelm,double xp,double yp,double zp,int i);
int iface3D(vector<element3D> &ELEM,int ielm,int jelm);
void qsorti3D(vector<point3D> &NODE,vector<element3D> &ELEM,int n,int *list);
void poly3D(vector<point3D> &NODE,vector<element3D> &ELEM,int *iv,int *kv,int ip,int *nelm,mpsconfig &CON);
void remove3D(vector<point3D> &NODE,vector<element3D> &ELEM,int *nelm,int iv,int *kv);
void fill3D(vector<point3D> &NODE,vector<element3D> &ELEM,int nelm);
void set_material(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm);
int poly3D_for_FINE3D(vector<point3D> &NODE,vector<element3D> &ELEM,int *iv,int *kv,int ip,int *nelm,mpsconfig &CON);
void FINE3D(vector<point3D> &NODE,vector<element3D> &ELEM,int KTJ,int KTE,int *node,int *nelm,mpsconfig &CON,double rrm,int startID);
int make_air_layer(vector<point3D> &NODE,vector<element3D> &ELEM,int *nelm,mpsconfig &CON,int *node_num,int KTE,double rrm,int KTJ);
//���b�V���`����t�@�C���ɏo��
void check_meshshape(int node,int *nelm,vector <point3D> &NODE,vector <element3D> &ELEM,mpsconfig &CON);

//�f���[�j����main�֐�
void delaun3D_main(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int KTJ,int KTE,int *node_num,int *nelm,int FINE_sw)
{
	/////////////�ߓ_���W�̐��K��
    double xmin=NODE[1].r[A_X];
    double ymin=NODE[1].r[A_Y];
    double zmin=NODE[1].r[A_Z];
    double xmax=xmin;
    double ymax=ymin;
    double zmax=zmin;

    //���W�̍ő�A�ŏ��l�����߂�
    for(int i=2;i<=*node_num;i++)
    {
        if(NODE[i].r[A_X]<xmin) xmin=NODE[i].r[A_X];
		else if(NODE[i].r[A_X]>xmax) xmax=NODE[i].r[A_X];
	
		if(NODE[i].r[A_Y]<ymin) ymin=NODE[i].r[A_Y];
		else if(NODE[i].r[A_Y]>ymax) ymax=NODE[i].r[A_Y];
	
		if(NODE[i].r[A_Z]<zmin) zmin=NODE[i].r[A_Z];
		else if(NODE[i].r[A_Z]>zmax) zmax=NODE[i].r[A_Z];
    }
    ////

    double rax=sqrt((xmax-xmin)*(xmax-xmin));	///X�������̐��@
    double ray=sqrt((ymax-ymin)*(ymax-ymin));	///Y�������̐��@
    double raz=sqrt((zmax-zmin)*(zmax-zmin));	///Z�������̐��@
    double rmax=rax;		///�ő吡�@
    if(ray>rmax) rmax=ray;
    if(raz>rmax) rmax=raz;      //������else�ɂ�����_�� �@

    //���W�ϊ�
    double rrm=1.000000/rmax;///�������������������邱�ƂŁA���l�덷�����点��E�E�H
    for(int i=1;i<=*node_num;i++)
    {   //   A/B�Ƃ����v�Z�������Ƃ��A�`�̒l�ɂ���Ĕ�����1/B�Ƃ����{�����������Ă���̂ł͂Ȃ����ƍl���āA���̂悤�ȏ������ɂ��Ă���
        NODE[i].r[A_X]=(NODE[i].r[A_X]-xmin)*rrm;
		NODE[i].r[A_Y]=(NODE[i].r[A_Y]-ymin)*rrm;
		NODE[i].r[A_Z]=(NODE[i].r[A_Z]-zmin)*rrm;
    }
    rax*=rrm;
    ray*=rrm;
    raz*=rrm;
    /////
	
	delaun3D(CON,NODE,ELEM, KTJ, KTE, rax, ray, raz,node_num,nelm, FINE_sw,rrm);
	//���b�V����������

	//���W�����ɖ߂�
    for(int i=1;i<=*node_num;i++)
    {
        NODE[i].r[A_X]=rmax*NODE[i].r[A_X]+xmin;
		NODE[i].r[A_Y]=rmax*NODE[i].r[A_Y]+ymin;
		NODE[i].r[A_Z]=rmax*NODE[i].r[A_Z]+zmin;
    }
}

//�f���[�j����
void delaun3D(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int KTJ,int KTE,double rax,double ray,double raz,int *node_num,int *nelm,int FINE_sw,double rrm)
{   
	cout<<"�f���[�j�����J�n----";

	unsigned int timeA=GetTickCount();	//�v�Z�J�n����
	double alpha=CON.get_boxalpha();	//�X�[�p�[�{�b�N�X�̑傫�������߂�W����
    int node=*node_num;


    for(int i=1;i<KTE;i++) ELEM[i].map=0;//������

    //�X�[�p�[�{�b�N�X�ւ̈ړ��ʂ��v�Z(�X�[�p�[�{�b�N�X�̑傫����rax*alpha,ray*alpha,raz*alpha)
    double xbar=0.500000*(alpha-1.000000)*rax;
    double ybar=0.500000*(alpha-1.000000)*ray;
    double zbar=0.500000*(alpha-1.000000)*raz;

    //�X�[�p�[�{�b�N�X�ֈړ� ����Ȃ镽�s�ړ��ł����āA�{���͕ω����Ă��Ȃ�
    for(int i=1;i<=node;i++)
    {
        NODE[i].r[A_X]+=xbar;
		NODE[i].r[A_Y]+=ybar;
		NODE[i].r[A_Z]+=zbar;
    }
    ///////////

    //�X�[�p�[�{�b�N�X��6�̎l�ʑ̂ɕ�������
    box3D(NODE,ELEM,alpha,nelm,rax,ray,raz,KTJ);//�v�f��6��������A�e�v�f�̑̐ρE�O�ڋ����S���W�����Ƃ܂���

    //�����ߓ_�𓱓����Ă���
    int *kv=new int[KTE];//�V�ߓ_���O�ڋ��ɂӂ��ޗv�f�Q
    int *istack=new int[KTE];//�ꎞ�z��

	int *flag=new int [node+1];
	for(int i=1;i<=node;i++) flag[i]=ON;
	int *jnb=new int[node+1];///�e�ߓ_�ɗאڂ���v�f���i�[ jnb=0�̋����ߓ_�͍ē���
	int num_ON=node;//�ǉ��ł��Ă��Ȃ��ߓ_�̐�
	int *unbreak=new int[KTE];//�j��ł��Ȃ��v�f
	for(int i=0;i<=KTE;i++) unbreak[i]=OFF;
//�����܂ł͂ł��Ă���B

/*	ofstream nod("node.dat",ios::app);
	for(int o=1;o<=node;o++)
	{
		nod<<NODE[o].r[A_X]<<", "<<NODE[o].r[A_Y]<<", "<<NODE[o].r[A_Z]<<endl;
	}
	nod.close();*/
	int asd=0;
	while(num_ON!=0)
	{
		for(int i=1;i<=node;i++)
		{   
			if(flag[i]==ON && (NODE[i].material==IRON || NODE[i].material==COIL || NODE[i].material==MAGNET))
			{
				int ip=i;
				double xp=NODE[ip].r[A_X];//��������ߓ_�̍��W
				double yp=NODE[ip].r[A_Y];
				double zp=NODE[ip].r[A_Z];
				
				//�V�ߓ_���܂ޗv�f�̒T��
				int loc=locate3D(NODE,ELEM,*nelm,xp,yp,zp,i);
			
				//////////�O�ڋ����ɐV�ߓ_���܂ޗv�f�̒��o
				int iv=0;
				int msk=0;
		
				iv++;//�O�ڋ����ɐV�ߓ_���܂ޗv�f��
				kv[iv]=loc;
				ELEM[loc].map=1;//map��1�̗v�f�́A�O�ڋ��ɐߓ_i���܂ނ��ǂ����������ς݂Ƃ�������
				msk++;
				istack[msk]=loc;
				
				while(msk!=0)
				{   
					int isk=istack[msk];//���ܒ��ڂ��Ă���v�f�̔ԍ�
					msk--;
					for(int j=1;j<=4;j++)
					{
						int jelm=ELEM[isk].elm[j];//isk�Ɛڂ���v�f
						if(unbreak[jelm]==OFF){//�j�󂵂��Ⴞ�߂ȖʂłȂ��Ȃ�
							if(jelm!=0)//���ꂪ�\�ʂłȂ��Ȃ�
							{
								if(ELEM[jelm].map!=1) //�܂��������ĂȂ��Ȃ�
								{   
				         
									double rad=ELEM[jelm].RR*(1.000000+ERR);//�O�ڋ����a�̂Q��
									double dst=ELEM[jelm].r[A_X]*ELEM[jelm].r[A_X]-2.000000*ELEM[jelm].r[A_X]*xp+xp*xp;///�O�ڋ����S�ƐV�ߓ_�̋���
								
									if(dst<rad)///������dst>rad�Ȃ��΂ɊO�ڋ��Ɋ܂܂Ȃ�
									{
										dst+=ELEM[jelm].r[A_Y]*ELEM[jelm].r[A_Y]-2.000000*ELEM[jelm].r[A_Y]*yp+yp*yp;
										if(dst<rad)
										{
											dst+=ELEM[jelm].r[A_Z]*ELEM[jelm].r[A_Z]-2.000000*ELEM[jelm].r[A_Z]*zp+zp*zp;
											if(dst<rad)//�O�ڋ����Ɋ܂�
											{
													iv++;//�O�ڋ����ɐV�ߓ_���܂ޗv�f�����{�P
													kv[iv]=jelm;//���X�g�ɂ����
													ELEM[jelm].map=1;//�O�ڋ��ɐV�_�܂ނƂ������邵
													msk++;
													istack[msk]=jelm;
											}
										}
									}
								}
							}
						}
					}
				}//�O�ڋ����ɐV�ߓ_���܂ޗv�f��iv�ƁA���̗v�f�ԍ�kv[iv]�����Ƃ܂���
				/////////////////
		
				int n0=*nelm;
				////����ꂽ���ʑ̂��l�ʑ̂ɕ�������
				poly3D(NODE,ELEM,&iv,kv,ip,nelm,CON); 
			}
		}
		/////////////////////////////////////////////
		for(int i=1;i<=node;i++)
		{   
			if(flag[i]==ON && (NODE[i].material==AIR ))
			{
				int ip=i;
				double xp=NODE[ip].r[A_X];//��������ߓ_�̍��W
				double yp=NODE[ip].r[A_Y];
				double zp=NODE[ip].r[A_Z];
				
				//�V�ߓ_���܂ޗv�f�̒T��
				int loc=locate3D(NODE,ELEM,*nelm,xp,yp,zp,i);
			
				//////////�O�ڋ����ɐV�ߓ_���܂ޗv�f�̒��o
				int iv=0;
				int msk=0;
		
				iv++;//�O�ڋ����ɐV�ߓ_���܂ޗv�f��
				kv[iv]=loc;
				ELEM[loc].map=1;//map��1�̗v�f�́A�O�ڋ��ɐߓ_i���܂ނ��ǂ����������ς݂Ƃ�������
				msk++;
				istack[msk]=loc;
				
				while(msk!=0)
				{   
					int isk=istack[msk];//���ܒ��ڂ��Ă���v�f�̔ԍ�
					msk--;
					for(int j=1;j<=4;j++)
					{
						int jelm=ELEM[isk].elm[j];//isk�Ɛڂ���v�f
						if(unbreak[jelm]==OFF){//�j�󂵂��Ⴞ�߂ȕӂłȂ��Ȃ�
							if(jelm!=0)//���ꂪ�\�ʂłȂ��Ȃ�
							{
								if(ELEM[jelm].map!=1) //�܂��������ĂȂ��Ȃ�
								{   
				         
									double rad=ELEM[jelm].RR*(1.000000+ERR);//�O�ڋ����a�̂Q��
									double dst=ELEM[jelm].r[A_X]*ELEM[jelm].r[A_X]-2.000000*ELEM[jelm].r[A_X]*xp+xp*xp;///�O�ڋ����S�ƐV�ߓ_�̋���
								
									if(dst<rad)///������dst>rad�Ȃ��΂ɊO�ڋ��Ɋ܂܂Ȃ�
									{
										dst+=ELEM[jelm].r[A_Y]*ELEM[jelm].r[A_Y]-2.000000*ELEM[jelm].r[A_Y]*yp+yp*yp;
										if(dst<rad)
										{
											dst+=ELEM[jelm].r[A_Z]*ELEM[jelm].r[A_Z]-2.000000*ELEM[jelm].r[A_Z]*zp+zp*zp;
											if(dst<rad)//�O�ڋ����Ɋ܂�
											{
												int o=ELEM[jelm].node[1];//�v�f���`������ړ_
												int p=ELEM[jelm].node[2];
												int q=ELEM[jelm].node[3];
												int r=ELEM[jelm].node[4];
												if(!(NODE[o].material==NODE[p].material && NODE[p].material==NODE[q].material && NODE[q].material==NODE[r].material &&(NODE[r].material==COIL || NODE[r].material==IRON|| NODE[r].material==FLUID || NODE[r].material==MAGNET ))){
												iv++;//�O�ڋ����ɐV�ߓ_���܂ޗv�f�����{�P
												kv[iv]=jelm;//���X�g�ɂ����
												ELEM[jelm].map=1;//�O�ڋ��ɐV�_�܂ނƂ������邵
												msk++;
												istack[msk]=jelm;
												}
											}
										}
									}
								}
							}
						}
					}
				}//�O�ڋ����ɐV�ߓ_���܂ޗv�f��iv�ƁA���̗v�f�ԍ�kv[iv]�����Ƃ܂���
				/////////////////
		
				int n0=*nelm;
				////����ꂽ���ʑ̂��l�ʑ̂ɕ�������
				poly3D(NODE,ELEM,&iv,kv,ip,nelm,CON); 
			}
		}
		/////////////////////////////////////////////
		for(int i=1;i<=node;i++)
		{   
//			if(flag[i]==ON && (NODE[i].material==MAGELAST || NODE[i].material==FLUID || NODE[i].material==HYPERELAST))	//HYPERELAST�ǉ�15/2/4
			if(flag[i]==ON && (NODE[i].material==MAGELAST || NODE[i].material==FLUID))	//HYPERELAST�ǉ�15/2/4
			{
				int ip=i;
				double xp=NODE[ip].r[A_X];//��������ߓ_�̍��W
				double yp=NODE[ip].r[A_Y];
				double zp=NODE[ip].r[A_Z];
				
				//�V�ߓ_���܂ޗv�f�̒T��
				int loc=locate3D(NODE,ELEM,*nelm,xp,yp,zp,i);
			
				//////////�O�ڋ����ɐV�ߓ_���܂ޗv�f�̒��o
				int iv=0;
				int msk=0;
		
				iv++;//�O�ڋ����ɐV�ߓ_���܂ޗv�f��
				kv[iv]=loc;
				ELEM[loc].map=1;//map��1�̗v�f�́A�O�ڋ��ɐߓ_i���܂ނ��ǂ����������ς݂Ƃ�������
				msk++;
				istack[msk]=loc;
				
				while(msk!=0)
				{   
					int isk=istack[msk];//���ܒ��ڂ��Ă���v�f�̔ԍ�
					msk--;
					for(int j=1;j<=4;j++)
					{
						int jelm=ELEM[isk].elm[j];//isk�Ɛڂ���v�f
						if(unbreak[jelm]==OFF){//�j�󂵂��Ⴞ�߂ȕӂłȂ��Ȃ�
							if(jelm!=0)//���ꂪ�\�ʂłȂ��Ȃ�
							{
								if(ELEM[jelm].map!=1) //�܂��������ĂȂ��Ȃ�
								{   
				         
									double rad=ELEM[jelm].RR*(1.000000+ERR);//�O�ڋ����a�̂Q��
									double dst=ELEM[jelm].r[A_X]*ELEM[jelm].r[A_X]-2.000000*ELEM[jelm].r[A_X]*xp+xp*xp;///�O�ڋ����S�ƐV�ߓ_�̋���
								
									if(dst<rad)///������dst>rad�Ȃ��΂ɊO�ڋ��Ɋ܂܂Ȃ�
									{
										dst+=ELEM[jelm].r[A_Y]*ELEM[jelm].r[A_Y]-2.000000*ELEM[jelm].r[A_Y]*yp+yp*yp;
										if(dst<rad)
										{
											dst+=ELEM[jelm].r[A_Z]*ELEM[jelm].r[A_Z]-2.000000*ELEM[jelm].r[A_Z]*zp+zp*zp;
											if(dst<rad)//�O�ڋ����Ɋ܂�
											{
									//				int o=ELEM[jelm].node[1];//�v�f���`������ړ_
									//				int p=ELEM[jelm].node[2];
									//				int q=ELEM[jelm].node[3];
									//				int r=ELEM[jelm].node[4];
									//				if(!(NODE[o].material==AIR || NODE[p].material==AIR || NODE[q].material==AIR || NODE[r].material==AIR)){//��蒼�����Ɏ������番������ɍs���̂�h��
													iv++;//�O�ڋ����ɐV�ߓ_���܂ޗv�f�����{�P
													kv[iv]=jelm;//���X�g�ɂ����
													ELEM[jelm].map=1;//�O�ڋ��ɐV�_�܂ނƂ������邵
													msk++;
													istack[msk]=jelm;
									//			}
											}
										}
									}
								}
							}
						}
					}
				}//�O�ڋ����ɐV�ߓ_���܂ޗv�f��iv�ƁA���̗v�f�ԍ�kv[iv]�����Ƃ܂���
				/////////////////
		
				int n0=*nelm;
				////����ꂽ���ʑ̂��l�ʑ̂ɕ�������
				poly3D(NODE,ELEM,&iv,kv,ip,nelm,CON); 
			}
		}
		////////////////////////////////////////////
		
		/////////////////////////////////////////////////////////////////
		//set_jnb3D(NODE,ELEM,node,*nelm,jnb);		//set_jnb3D()���g���ƁA�v�f���\������ߓ_���X�[�p�[�{�b�N�̒��_�������Ƃ��Ƀo�O��B�����牺�̂悤�Ɍv�Z
		for(int i=1;i<=node;i++) jnb[i]=0;//������
		for(int je=1;je<=*nelm;je++)
		{
			for(int j=1;j<=4;j++)
			{
				int n=ELEM[je].node[j];
				if(n<=node) jnb[n]++;
			}
		}
		
		for(int i=1;i<=node;i++)
		{
			if(flag[i]==ON) num_ON--;
			flag[i]=OFF;
					
			if(jnb[i]==0)
			{
				flag[i]=ON;
				num_ON++;
			}
		}
	//	cout<<"num_ON="<<num_ON<<endl;
		asd++;
		if(asd>=1000) num_ON=0;
	}
	
	delete [] flag;
	delete [] jnb;

    //�X�[�p�[�{�b�N�X�𒸓_�ɂ��l�ʑ̂���������
    //������iv��kv�̖�ڂ�ς���Biv�͏��������l�ʑ̐��ł���Akv�͂��̃��X�g�B
    int iv=0;
    for(int i=1;i<=*nelm;i++)  kv[i]=0;//������
    
    for(int i=1;i<=*nelm;i++)
    {   //�ǂꂩ�ЂƂł��X�[�p�[�{�b�N�X�̒��_���܂܂�Ă�����
        if(ELEM[i].node[1]>KTJ || ELEM[i].node[2]>KTJ||ELEM[i].node[3]>KTJ||ELEM[i].node[4]>KTJ)
		{
			iv++;
			kv[iv]=i;
		}
    }
	int before_nelm=*nelm;//�����O�v�f���L��

  //8

	cout<<"ok �����O�v�f��="<<before_nelm<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

    remove3D(NODE,ELEM,nelm,iv,kv);
    ///////////////////
    
    //�v�f�����܂���������Ă��邩�`�F�b�N����
    fill3D(NODE,ELEM,*nelm);

	//�v�f�ގ�����
    set_material(CON,NODE,ELEM,node,*nelm);
	check_meshshape(node,nelm,NODE,ELEM,CON);//�ގ�
	//set_depth(CON,NODE,ELEM,node,*nelm,depth,KTE);//�v�f�[������

    ///���W�����ɖ߂�  ���������S�ɖ߂��̂�main��
    for(int i=1;i<=node;i++)
    {
        NODE[i].r[A_X]-=xbar;
		NODE[i].r[A_Y]-=ybar;
		NODE[i].r[A_Z]-=zbar;
    }

    delete [] kv;
    delete [] istack;
	delete [] unbreak;
	////////////////////////////////////////////���b�V���̎����ו���
	if(FINE_sw==ON) FINE3D(NODE,ELEM,KTJ,KTE,node_num,nelm,CON,rrm,1);

}

//�X�[�p�[�{�b�N�X�쐬�֐�
void box3D(vector<point3D> &NODE,vector<element3D> &ELEM,double alpha,int *nelm,double rax,double ray,double raz,int KTJ)
{
	//�X�[�p�[�{�b�N�X�̑傫����rax*alpha,ray*alpha,raz*alpha

    *nelm=6;
    double xone=alpha*rax;//�ő���W
    double yone=alpha*ray;
    double zone=alpha*raz;
//	cout<<"xone="<<xone<<", yone="<<yone<<", zone="<<zone<<endl;
    /////�l�ʑ̏��쐬
    //�ߓ_���W
    NODE[KTJ+1].r[A_X]=0.000000;
    NODE[KTJ+1].r[A_Y]=0.000000;
    NODE[KTJ+1].r[A_Z]=0.000000;

    NODE[KTJ+2].r[A_X]=xone;
    NODE[KTJ+2].r[A_Y]=0.000000;
    NODE[KTJ+2].r[A_Z]=0.000000;
    
    NODE[KTJ+3].r[A_X]=xone;
    NODE[KTJ+3].r[A_Y]=yone;
    NODE[KTJ+3].r[A_Z]=0.000000;

    NODE[KTJ+4].r[A_X]=0.000000;
    NODE[KTJ+4].r[A_Y]=yone;
    NODE[KTJ+4].r[A_Z]=0.000000;

    NODE[KTJ+5].r[A_X]=0.000000;
    NODE[KTJ+5].r[A_Y]=0.000000;
    NODE[KTJ+5].r[A_Z]=zone;

    NODE[KTJ+6].r[A_X]=xone;
    NODE[KTJ+6].r[A_Y]=0.000000;
    NODE[KTJ+6].r[A_Z]=zone;

    NODE[KTJ+7].r[A_X]=xone;
    NODE[KTJ+7].r[A_Y]=yone;
    NODE[KTJ+7].r[A_Z]=zone;

    NODE[KTJ+8].r[A_X]=0.000000;
    NODE[KTJ+8].r[A_Y]=yone;
	NODE[KTJ+8].r[A_Z]=zone;
    ///////

    ////�v�f�\�ߓ_���W
    ELEM[1].node[1]=KTJ+2;
    ELEM[1].node[2]=KTJ+7;
    ELEM[1].node[3]=KTJ+5;
    ELEM[1].node[4]=KTJ+6;
    
    ELEM[2].node[1]=KTJ+1;
    ELEM[2].node[2]=KTJ+2;
    ELEM[2].node[3]=KTJ+3;
    ELEM[2].node[4]=KTJ+5;
    
    ELEM[3].node[1]=KTJ+2;
    ELEM[3].node[2]=KTJ+3;
    ELEM[3].node[3]=KTJ+5;
    ELEM[3].node[4]=KTJ+7;
    
    ELEM[4].node[1]=KTJ+5;
    ELEM[4].node[2]=KTJ+4;
    ELEM[4].node[3]=KTJ+8;
    ELEM[4].node[4]=KTJ+7;
    
    ELEM[5].node[1]=KTJ+1;
    ELEM[5].node[2]=KTJ+3;
    ELEM[5].node[3]=KTJ+4;
    ELEM[5].node[4]=KTJ+5;
    
    ELEM[6].node[1]=KTJ+4;
    ELEM[6].node[2]=KTJ+3;
    ELEM[6].node[3]=KTJ+7;
    ELEM[6].node[4]=KTJ+5;
    /////////

    ////�v�f�\�v�f���
    ELEM[1].elm[1]=0;
    ELEM[1].elm[2]=0;
    ELEM[1].elm[3]=0;
    ELEM[1].elm[4]=3;
    
    ELEM[2].elm[1]=3;
    ELEM[2].elm[2]=5;
    ELEM[2].elm[3]=0;
    ELEM[2].elm[4]=0;
    
    ELEM[3].elm[1]=6;
    ELEM[3].elm[2]=1;
    ELEM[3].elm[3]=0;
    ELEM[3].elm[4]=2;
    
    ELEM[4].elm[1]=0;
    ELEM[4].elm[2]=0;
    ELEM[4].elm[3]=6;
    ELEM[4].elm[4]=0;
    
    ELEM[5].elm[1]=6;
    ELEM[5].elm[2]=0;
    ELEM[5].elm[3]=2;
    ELEM[5].elm[4]=0;
    
    ELEM[6].elm[1]=3;
    ELEM[6].elm[2]=4;
    ELEM[6].elm[3]=5;
    ELEM[6].elm[4]=0;
    /////

 

    for(int i=1;i<=6;i++)
    {
        int ia=ELEM[i].node[1];
		int ib=ELEM[i].node[2];
		int ic=ELEM[i].node[3];
		int ip=ELEM[i].node[4];
		ELEM[i].volume=volume3D(NODE,ia,ib,ic,ip);//�̐ς�6�{�ł��邱�Ƃɒ���
		sphere3D(NODE,ELEM,ia,ib,ic,ip,i);
    }
}

///�̐όv�Z�֐�
double volume3D(vector<point3D> &NODE,int ia,int ib,int ic,int ip)
{
    ////4���_ia,ib,ic,ip����̐ς����Ƃ߂�B�������̐ς�6�{�̒l��Ԃ�(�ǂ����v�Z�ɕK�v�Ȃ̂͂u���U������)
    
    double xa=NODE[ia].r[A_X];
    double ya=NODE[ia].r[A_Y];
    double za=NODE[ia].r[A_Z];
    
    double xb=NODE[ib].r[A_X];
    double yb=NODE[ib].r[A_Y];
    double zb=NODE[ib].r[A_Z];
    
    double xc=NODE[ic].r[A_X];
    double yc=NODE[ic].r[A_Y];
    double zc=NODE[ic].r[A_Z];
    
    double xp=NODE[ip].r[A_X];
    double yp=NODE[ip].r[A_Y];
    double zp=NODE[ip].r[A_Z];
    
    double va=xb*yc*zp+xa*ya*zp+xb*ya*za+xa*yc*za-(xb*yc*za+xa*ya*za+xb*ya*zp+xa*yc*zp);
    double vb=yb*zc*xp+ya*za*xp+yb*za*xa+ya*zc*xa-(yb*zc*xa+ya*za*xa+yb*za*xp+ya*zc*xp);
    double vc=zb*xc*yp+za*xa*yp+zb*xa*ya+za*xc*ya-(zb*xc*ya+za*xa*ya+zb*xa*yp+za*xc*yp);
    
    double wa=xb*zc*ya+xa*za*ya+xb*za*yp+xa*zc*yp-(xb*zc*yp+xa*za*yp+xb*za*ya+xa*zc*ya);
    double wb=yb*xc*za+ya*xa*za+yb*xa*zp+ya*xc*zp-(yb*xc*zp+ya*xa*zp+yb*xa*za+ya*xc*za);
    double wc=zb*yc*xa+za*ya*xa+zb*ya*xp+za*yc*xp-(zb*yc*xp+za*ya*xp+zb*ya*xa+za*yc*xa);
    
    double volume=va+vb+vc+wa+wb+wc;
    
    return volume;
}

///�O�ڋ����Ұ��v�Z�֐�
void sphere3D(vector<point3D> &NODE,vector<element3D> &ELEM,int ia,int ib,int ic,int ip,int i)
{
    double xa=NODE[ia].r[A_X];
    double ya=NODE[ia].r[A_Y];
    double za=NODE[ia].r[A_Z];
    
    double xb=NODE[ib].r[A_X];
    double yb=NODE[ib].r[A_Y];
    double zb=NODE[ib].r[A_Z];
    
    double xc=NODE[ic].r[A_X];
    double yc=NODE[ic].r[A_Y];
    double zc=NODE[ic].r[A_Z];
    
    double xp=NODE[ip].r[A_X];
    double yp=NODE[ip].r[A_Y];
    double zp=NODE[ip].r[A_Z];
    
    double p11=yc*zp+ya*za+yp*za+ya*zc-(yc*za+ya*zp+yp*zc+ya*za);
    double p12=xp*zc+xa*za+xc*za+xa*zp-(xp*za+xa*zc+xc*zp+xa*za);
    double p13=xc*yp+xa*ya+xp*ya+xa*yc-(xc*ya+xa*yp+xp*yc+xa*ya);
    double p21=yp*zb+ya*za+yb*za+ya*zp-(yp*za+ya*zb+yb*zp+ya*za);
    double p22=xb*zp+xa*za+xp*za+xa*zb-(xb*za+xa*zp+xp*zb+xa*za);
    double p23=xp*yb+xa*ya+xb*ya+xa*yp-(xp*ya+xa*yb+xb*yp+xa*ya);
    double p31=yb*zc+ya*za+yc*za+ya*zb-(yb*za+ya*zc+yc*zb+ya*za);
    double p32=xc*zb+xa*za+xb*za+xa*zc-(xc*za+xa*zb+xb*zc+xa*za);
    double p33=xb*yc+xa*ya+xc*ya+xa*yb-(xb*ya+xa*yc+xc*yb+xa*ya);
    
    double xyza=xa*xa+ya*ya+za*za;
    double aa=0.5000000*(xb*xb+yb*yb+zb*zb-xyza);
    double bb=0.5000000*(xc*xc+yc*yc+zc*zc-xyza);
    double cc=0.5000000*(xp*xp+yp*yp+zp*zp-xyza);
    
    double xx=p11*aa+p21*bb+p31*cc;
    double yy=p12*aa+p22*bb+p32*cc;
    double zz=p13*aa+p23*bb+p33*cc;
    
    double determ=ELEM[i].volume;//�̐ς̂U�{
    double xv=xx/determ; 
    double yv=yy/determ;
    double zv=zz/determ;
    
    ELEM[i].r[A_X]=xv; //�O�ڋ����S���W
    ELEM[i].r[A_Y]=yv;
    ELEM[i].r[A_Z]=zv;
    
    ELEM[i].RR=xa*xa+xv*xv+ya*ya+yv*yv+za*za+zv*zv-2.000000*(xa*xx+ya*yy+za*zz)/determ;
}

///�ߓ_����v�f�T���֐�
int locate3D(vector<point3D> &NODE,vector<element3D> &ELEM,int nelm,double xp,double yp,double zp)
{
    int itet=nelm;//��ԍŌ�ɐ������ꂽ�v�f�ԍ�
    int flag=1;
	int count=0;
	int maxsize=ELEM.size();
	double err=ERR;
    while(flag==1)
    {
        flag=0;
        for(int n=1;n<=4;n++)//��P�`��S�ʂ𒲂ׂ�
        {
            if(flag==0)
			{
                int i=ELEM[itet].node[n%4+1];
				int j=ELEM[itet].node[4-(n-1)/2*2];
				int k=ELEM[itet].node[3-n/2%2*2];
	
				double xi=NODE[i].r[A_X];//�ߓ_i�̍��W
				double yi=NODE[i].r[A_Y];
				double zi=NODE[i].r[A_Z];
	
				double xj=NODE[j].r[A_X];//�ߓ_j�̍��W
				double yj=NODE[j].r[A_Y];
				double zj=NODE[j].r[A_Z];
	
				double xk=NODE[k].r[A_X];//�ߓ_k�̍��W
				double yk=NODE[k].r[A_Y];
				double zk=NODE[k].r[A_Z];
	
				double a=yi*zj+yj*zk+yk*zi-(yi*zk+yj*zi+yk*zj);
				double b=zi*xj+zj*xk+zk*xi-(zi*xk+zj*xi+zk*xj);
				double c=xi*yj+xj*yk+xk*yi-(xi*yk+xj*yi+xk*yj);
				double d=-a*xi-b*yi-c*zi;
				if(a*xp+b*yp+c*zp+d<-err)
				{   //�O��
					itet=ELEM[itet].elm[n];//itet�̑�n�ʂɗאڂ��Ă���v�f�����̌����ΏۂƂ���
					flag=1;
					count++;
					if(count>maxsize){ 
						cout<<"locate�Ŗ������[�vcount="<<count<<" �T���v�f="<<itet<<"�ڋߋ�����"<<a*xp+b*yp+c*zp+d<<endl;
						cout<<"xp="<<xp<<", yp="<<yp<<", zp="<<zp<<", err"<<err<<endl;
					}
					if(itet==0) return 0;	//�̈�̊O�ɂł��ꍇ��0��Ԃ��B�ʏ�̃f���[�j�ł͂��肦�Ȃ��󋵂����A���Ƃ���remesh�̈�ɂ������C�w�������ɂ́A�@���x�N�g����̐ߓ_�ʒu���ÓI�v�f�̈�ɂ͂���Ƃ��̂悤�ȃG���[�ƂȂ�
				}
			}
        }
		
    }
    return itet;
}

///�ߓ_����v�f�T���֐��ǉ��ړ_��ID�t
int locate3D(vector<point3D> &NODE,vector<element3D> &ELEM,int nelm,double xp,double yp,double zp,int n)
{
    int itet=nelm;//��ԍŌ�ɐ������ꂽ�v�f�ԍ�
    int flag=1;
	int count=0;
	int maxsize=ELEM.size();
	double err=ERR;
    while(flag==1)
    {
        flag=0;
        for(int n=1;n<=4;n++)//��P�`��S�ʂ𒲂ׂ�
        {
            if(flag==0)
			{
                int i=ELEM[itet].node[n%4+1];
				int j=ELEM[itet].node[4-(n-1)/2*2];
				int k=ELEM[itet].node[3-n/2%2*2];
	
				double xi=NODE[i].r[A_X];//�ߓ_i�̍��W
				double yi=NODE[i].r[A_Y];
				double zi=NODE[i].r[A_Z];
	
				double xj=NODE[j].r[A_X];//�ߓ_j�̍��W
				double yj=NODE[j].r[A_Y];
				double zj=NODE[j].r[A_Z];
	
				double xk=NODE[k].r[A_X];//�ߓ_k�̍��W
				double yk=NODE[k].r[A_Y];
				double zk=NODE[k].r[A_Z];
	
				double a=yi*zj+yj*zk+yk*zi-(yi*zk+yj*zi+yk*zj);
				double b=zi*xj+zj*xk+zk*xi-(zi*xk+zj*xi+zk*xj);
				double c=xi*yj+xj*yk+xk*yi-(xi*yk+xj*yi+xk*yj);
				double d=-a*xi-b*yi-c*zi;
				if(a*xp+b*yp+c*zp+d<-err)
				{   //�O��
					itet=ELEM[itet].elm[n];//itet�̑�n�ʂɗאڂ��Ă���v�f�����̌����ΏۂƂ���
					flag=1;
					count++;
					if(count>maxsize){ 
						cout<<"locate�Ŗ������[�vcount="<<count<<" �T���v�f="<<itet<<"�ڋߋ�����"<<a*xp+b*yp+c*zp+d<<endl;
						cout<<"xp="<<xp<<", yp="<<yp<<", zp="<<zp<<", err��"<<err<<"�ǉ��v�f�ގ���"<<NODE[n].material<<endl;
					}
					if(itet==0) return 0;	//�̈�̊O�ɂł��ꍇ��0��Ԃ��B�ʏ�̃f���[�j�ł͂��肦�Ȃ��󋵂����A���Ƃ���remesh�̈�ɂ������C�w�������ɂ́A�@���x�N�g����̐ߓ_�ʒu���ÓI�v�f�̈�ɂ͂���Ƃ��̂悤�ȃG���[�ƂȂ�
				}
			}
        }
		
    }
    return itet;
}

//���ʑ̕����֐�
void poly3D(vector<point3D> &NODE,vector<element3D> &ELEM,int *iv,int *kv,int ip,int *nelm,mpsconfig &CON)
{   
	////���̊֐��̑O�̒i�K�ŁA�V�ߓ_���O�ڋ��Ɋ܂ގl�ʑ̐�iv�Ƃ��̗v�f�ԍ�kv[iv]�����Ƃ܂��Ă���B
	////�܂��A�V�ߓ_���O�ڋ��Ɋ܂ގl�ʑ̂�ELEM[i].map=1�ƂȂ��Ă���
    
	int memory=CON.get_poly_memory();		//���I�Ɋm�ۂ��郁������
	
    int ix=0;				//�\�ʂ̐� ��ʂ�iv=ix�Ƃ͂Ȃ�Ȃ�.�ЂƂ̗v�f�������̕\�ʂ�S���̂ŁB(iv<=ix�ł���)
	int *imen[3+1];
	for(int i=1;i<=3;i++) imen[i]=new int [memory]; //���ʑ̕\�ʎO�p�`�̐ߓ_�ԍ��i�[
    int *jmen=new int [memory];		//���ʑ̕\�ʎO�p�`�ɗאڂ���l�ʑ̔ԍ��i�[
    int *kmen=new int [memory];		//���ʑ̕\�ʎO�p�`�ɗאڂ���l�ʑ̗̂אږʔԍ� (����͑扽�ʂŎ����Ɛڂ��Ă��邩)
    double *vol=new double [memory];		//���ʑ̂̑̐ς̂U�{
    
	//if(ip==7884) data_avs4(CON,ip,NODE,ELEM, *iv,kv);

    ///////���ʑ̕\�ʐ�ix�Ƃ�����\�����钸�_�����Ƃ߂�
    int flag=0;
    while(flag==0)
    {   
        ix=0;
        flag=1;
        for(int i=1;i<=*iv;i++)//*iv�͐V�_���O�ڋ����Ɋ܂ޗv�f�̐�
        {
			if(flag==1)
			{
        		int ielm=kv[i];//�V�_���O�ڋ����Ɋ܂ޗv�f�ԍ�
				for(int j=1;j<=4;j++)
				{
					if(flag==1)
					{
						int jelm=ELEM[ielm].elm[j];//ielm�v�f�ɗאڂ���v�f�ԍ�
						int ia=ELEM[ielm].node[j%4+1];//ielm��jelm�̐ڂ���O�p�`���\������ߓ_�ԍ�  ia,ib,ic,ip�𒸓_�Ƃ��ĐV�����v�f�������
						int ib=ELEM[ielm].node[4-(j-1)/2*2];
						int ic=ELEM[ielm].node[3-(j/2%2)*2];
	                
						if(jelm==0)//�\�ʂȂ�ielm�͕\�ʗv�f�ł���Ƃ킩��
						{
							ix++;
							imen[1][ix]=ia;	//���ʑ̕\�ʎO�p�`�̒��_�ԍ�
							imen[2][ix]=ib;
							imen[3][ix]=ic;
							jmen[ix]=0;		//���ʑ̕\�ʎO�p�`�ɗאڂ���v�f�ԍ�
							kmen[ix]=0;		//���ʑ̕\�ʎO�p�`�ɗאڂ���v�f�̗אږʔԍ�  �����0������A�����ł�0����
		
							vol[ix]=volume3D(NODE,ia,ib,ic,ip);
						}
						else if(ELEM[jelm].map==0)//ielm�ɐڂ���jelm�͐V�_���O�ڋ��Ɋ܂܂Ȃ��B�܂�ielm��jelm�̋��E�͑��ʑ̕\�ʂƂ�������
						{
							ix++;
							imen[1][ix]=ia;
							imen[2][ix]=ib;
							imen[3][ix]=ic;
							jmen[ix]=jelm;
							kmen[ix]=iface3D(ELEM,ielm,jelm);//jelm��ielm�ɐڂ���ʔԍ�
							
							if(ix>=memory)cout<<" ix>memory"<<endl;
							vol[ix]=volume3D(NODE,ia,ib,ic,ip);
							if(vol[ix]<=ERR)//���ȏ��}3.10�̂悤�ɁA�s�K�؂Ȏl�ʑ̂�z�肵�Ă��܂��A���ʑ̐ς����ɂȂ�B���̏ꍇ�A���ȏ��ɏ����Ă���Ƃ���A�ŏ�������Ȃ���(goto 10)
							{   							
								*iv=(*iv)+1;//���X�g�𑝂₷ *iv++�͂��߁H
								
								kv[*iv]=jelm;//�{���͊O�ڋ����ɐV�ߓ_���܂܂Ȃ����A���ʑ̕������\�ɂ��邽�߂Ƀ��X�g�ɂ����
								
								ELEM[jelm].map=1;
								flag=0;
								//if(ip==7884) cout<<jelm<<endl;
							}
						}
					}
				}
			}
        }
    }

    //////���ʑ̕\�ʐ�ix�Ƃ�����\�����钸�_�����Ƃ܂���
  // if(ip==7884) data_avs4(CON,ip,NODE,ELEM, *iv,kv);
    /////////////�\�ʂ̒��_�ƐV�ߓ_���Ȃ���.����Ƒ��ʑ̂���ix�̗v�f�����������

    int ibound=ix;//ix�̑���B�܂�\�ʂ̐���\��
    for(int i=*(iv)+1;i<=ibound;i++) ///iv�̗v�f���j�󂳂�ĐV����ibound�̗v�f����������邩��A������v�f����ibound-iv
    {   
        *nelm=*nelm+1;		//*nelm++�Ƃ����������ł͂��߁H
        kv[i]=*nelm;		//�V�ߓ_���O�ڋ��Ɋ܂ޗv�f���X�g��������B
		ELEM[*nelm].map=1;	//���������v�f�͕K���V�ߓ_���O�ڋ��Ɋ܂�(�Ƃ��Ă���).(���̍s�Ӗ��Ȃ��Ȃ��H)
    }
    
    for(int i=1;i<=ibound;i++) ELEM[kv[i]].map=0;//�}�b�s���O�̏�����
    
    for(int i=1;i<=ibound;i++)//�v�f��񐶐�
    {   
        int ielm=kv[i];
		double determ=vol[i];
		int ia=imen[1][i];
		int ib=imen[2][i];
		int ic=imen[3][i];
		ELEM[ielm].node[1]=ia;
		ELEM[ielm].node[2]=ib;
		ELEM[ielm].node[3]=ic;
		ELEM[ielm].node[4]=ip;		//�V�_�͂S�Ԗڂƒ�`	
		ELEM[ielm].elm[4]=jmen[i];	
		if(jmen[i]!=0) ELEM[jmen[i]].elm[kmen[i]]=ielm;
		///����ł����H
		ELEM[ielm].volume=vol[i];
	
		sphere3D(NODE,ELEM,ia,ib,ic,ip,ielm);//�O�ڋ��̒��S�i�{���m�C�_)�Ɣ��a�̓����v�Z
    }
    ///////////////////
    
    //�v�f�\�v�f�֌W�C��/////��̏����ő�4�ʂŐڂ���v�f�ԍ��͂킩���Ă���̂ŁA�c������߂�
	//						�����ŁA1�`3�ʂ͑��ʑ̂��\������v�f�Ƃ̋��E�ʂł��邱�Ƃɒ���
    ix=0;
    for(int i=1;i<=ibound;i++)
    {
        int ielm=kv[i];
		//if(ip==49893) if(ielm==11564) cout<<"RR"<<endl;
		for(int j=1;j<=3;j++)//ELEM[ielm].elm[4]�͂��łɂ��Ƃ܂�������A����ȊO�����Ƃ߂�
		{
			///ELEM[ielm].node[4]=ip�ł���
			ELEM[ielm].elm[j]=-1;				//������
			int ia=ELEM[ielm].node[j%3+1];		//j=1,2,3�̂Ƃ��A2,3,1�̏�
			int ib=ELEM[ielm].node[(j%3+1)%3+1];//j=1,2,3�̂Ƃ��A3,1,2�̏�
			int flag=0;
			for(int k=1;k<=ix;k++)
			{
				if(flag==0)
				{
					int ja=imen[1][k];
					int jb=imen[2][k];
					if(ia==ja && ib==jb)//�ߓ_����v������
					{
					    ELEM[ielm].elm[j]=jmen[k];//���炩����ؽĂ��Ă����������i�[
					    ELEM[jmen[k]].elm[kmen[k]]=ielm;
					    imen[1][k]=imen[1][ix];//k�Ԗڂ̏��͂����s�v�B�Ȃ̂Ŕz��̈�ԍŌ�̏���k�Ԗڂɂ����Ă��āA����܂ł̏��͔j������
					    imen[2][k]=imen[2][ix];
					    jmen[k]=jmen[ix];
					    kmen[k]=kmen[ix];
					    ix--;	//�҂��Ӑ�����
						flag=1;	//ELEM[ielm].elm[j]�͂��Ƃ܂����̂ŁA���̃l�X�g�ɓ���K�v�͂Ȃ��̂�flag=1
					}
				}
			}
			if(flag==0)
			{
			    ix++;			//�����ł�ix�́A[�אڊ֌W�𖞂����v�f]���܂��Ă���[��]�̐���\���B
				imen[1][ix]=ib;	//�����̐ߓ_�̕��т��L�������A�ʂ̗v�f�����̕��т𖞂����̂�҂Bib��ia�̕��т��t�ɂ��Ă��邱�Ƃɒ���
				imen[2][ix]=ia;
				jmen[ix]=ielm;
				kmen[ix]=j;
			}
		}
    }///�v�f-�v�f�֌W�C������

    /////////�V���ɍ��ꂽ�l�ʑ̂̌����A�Â����̂�菭�Ȃ��Ȃ����ꍇ(���ʑ̂��\�������Ƃ��A�ǂ̖ʂ����ʑ̂̋��E�ł͂Ȃ������v�f�����݂���Ƃ�)
    if(*iv>ibound)
    {   
		//cout<<ip<<endl;
        int ir=*(iv)-ibound;	//ir�̗v�f���폜����Ȃ��Ă͂Ȃ�Ȃ�.���������ʂɍ폜�����̂ł́A�v�f�z���[��]�������邱�ƂɂȂ�

		for(int i=1;i<=ir;i++)
		{
			kv[i]=kv[ibound+i];
			ELEM[kv[i]].map=kv[i];//map�Ƃ��ĂƂ肠�������̗v�f�ԍ����i�[
		}
		///��ő������map�̒l(�v�f�ԍ�)�����������ɂȂ�т�����B���������ۂɕ��ёւ��̎w�W�ɂ��Ă���̂�kv[i]�ł���Bkv[i]�Ɨv�f�ԍ��͂����ł͓�����
		qsorti3D(NODE,ELEM,ir,kv);
	
		for(int i=1;i<=ir;i++)
		{   
			int ielm=kv[ir-i+1];//i��������ɂ�Air-i+1��ir����1���������Ȃ��Ă���
			ELEM[ielm].map=0;	//��̂ق���map�̏��������s���Ă��邪�A����͕��ёւ��ɕK�v����������ŁA���яI��������炱���ł͏�����
			if(ielm!=*nelm)		//ielm�Ԗڂ̗v�f�ɁAnelm�Ԗڂ̗v�f�����㏑��
			{   
			    ELEM[ielm].r[A_X]=ELEM[*nelm].r[A_X];
				ELEM[ielm].r[A_Y]=ELEM[*nelm].r[A_Y];
				ELEM[ielm].r[A_Z]=ELEM[*nelm].r[A_Z];
				ELEM[ielm].RR=ELEM[*nelm].RR;
				for(int j=1;j<=4;j++)
				{   
				    ELEM[ielm].node[j]=ELEM[*nelm].node[j];
				    int jelm=ELEM[*nelm].elm[j];
				    ELEM[ielm].elm[j]=jelm;
				    if(jelm!=0)
				    {   
				        int N=iface3D(ELEM,*nelm,jelm);
						ELEM[jelm].elm[N]=ielm;
				    }
				}
			}
			*nelm=*nelm-1;
		}
    }///////////////////*/

	delete [] jmen;
	delete [] kmen;
	delete [] vol;

	for(int i=1;i<=3;i++) delete [] imen[i];
}		

//�אږʔԍ��v�Z�֐�
int iface3D(vector<element3D> &ELEM,int ielm,int jelm)
{
    int iface=-1;//return����ϐ�
    for(int n=1;n<=4;n++)
    {
        if(ELEM[jelm].elm[n]==ielm)
		{
			iface=n;
			// break;
		}
    }
    if(iface<0) cout<<"error in inface �v�f�ԍ�="<<jelm<<endl;
    return iface;
}

//���בւ��֐�
void qsorti3D(vector<point3D> &NODE,vector<element3D> &ELEM,int n,int *list)
{
    //list=kv[i]
    int maxstk=100;//�p�����[�^
    
    int ll=1;
    int lr=n;//���ёւ����
    int istk=0;
    
    int *ilst=new int[maxstk];
    int *irst=new int[maxstk];
    
    int flag2=0;
    while(flag2==0)  //10 continue
    {
        flag2=1;
        while(ll<lr)
        {
            int nl=ll;
			int nr=lr;
			int lm=(ll+lr)/2;
			int iguess=ELEM[list[lm]].map;
	
			int flag=0;
			while(flag==0)
			{   
				flag=1;
				while(ELEM[list[nl]].map<iguess)
				{
					nl++;
				}
	
				while(iguess<ELEM[list[nr]].map)
				{
					nr--;
				}
	
				if(nl<nr-1)
				{
					int ltemp=list[nl];
					list[nl]=list[nr];
					list[nr]=ltemp;
					nl++;
					nr--;
				flag=0;
				}
			}
	
			if(nl<=nr)
            {
                if(nl<nr)
				{
					int ltemp=list[nl];
					list[nl]=list[nr];
					list[nr]=ltemp;
				}
				nl++;
				nr--;
            }
    
            istk++;
			if(istk>maxstk) cout<<"ERROR(qsorti) istk>maxstk"<<endl;
	
			if(nr<lm)
			{
				ilst[istk]=nl;
				irst[istk]=lr;
				lr=nr;
			}
			else
			{
				ilst[istk]=ll;
				irst[istk]=nr;
				ll=nl;
			}
	
        }
        if(istk!=0)
        {
            ll=ilst[istk];
			lr=irst[istk];
			istk--;
			flag2=0;
        }
    }
    delete [] ilst;
    delete [] irst;
}

//�s�v�v�f�����֐�
void remove3D(vector<point3D> &NODE,vector<element3D> &ELEM,int *nelm,int iv,int *kv)
{
	//iv:�����v�f���@kv[i]:�����v�f�ԍ�

    int m=0;///��������Ȃ��v�f�� �v�Z�̓s����2��ޕK�v
    int n=0;///��������Ȃ��v�f��
    for(int i=1;i<=*nelm;i++) ELEM[i].map=1;	//������ �����ł�map�͏�����̗v�f�ԍ���\��
    
    ///////�e�v�f��map(�V�����v�f�ԍ�)�����Ƃ߂�

    for(int i=1;i<=iv;i++) ELEM[kv[i]].map=0;//���������̂�����0�ԖڂƂȂ�
    
    for(int i=1;i<=*nelm;i++)
    {
        if(ELEM[i].map!=0)//��������Ȃ��Ȃ�
		{
			m++;
			ELEM[i].map=m;
		}
    }/////////
    
    ////�e���̂Ђ����@
    for(int i=1;i<=*nelm;i++)
    {
        if(ELEM[i].map!=0)//��������Ȃ��Ȃ�
		{
			n++;
			ELEM[n].r[A_X]=ELEM[i].r[A_X];
			ELEM[n].r[A_Y]=ELEM[i].r[A_Y];
			ELEM[n].r[A_Z]=ELEM[i].r[A_Z];
			ELEM[n].RR=ELEM[i].RR;
			for(int ia=1;ia<=4;ia++)
			{
				//�v�f-�ߓ_���̃R�s�[
				ELEM[n].node[ia]=ELEM[i].node[ia];

				//�v�f-�v�f���̃R�s�[
				if(ELEM[i].elm[ia]==0) ELEM[n].elm[ia]=0;	//�\�ʂƐڂ��Ă���Ȃ炻�̏���P���ɃR�s�[
				else //�����łȂ��Ȃ�A�v�f�ԍ��̕ω����l���ɂ��ꂽ�R�s�[���s��Ȃ��Ƃ����Ȃ�
				{
					ELEM[n].elm[ia]=ELEM[ELEM[i].elm[ia]].map;//�����ŁA�X�[�p�[�{�b�N�X���_���܂ޗv�f(�s�v�v�f)��map=0������A�s�v�v�f�ɗאڂ���v�f�͗אڗv�f�Ƃ��ă[�����i�[�����
				}
			} 
		}
    }///////////
    
	//�s�v�v�f�̏�񏉊���
    for(int i=n+1;i<=*nelm;i++)
    {
        ELEM[i].r[A_X]=0.00000;
		ELEM[i].r[A_Y]=0.00000;
		ELEM[i].r[A_Z]=0.00000;
		ELEM[i].RR=0.00000;
		for(int ia=1;ia<=4;ia++)
		{
			ELEM[i].node[ia]=0;
			ELEM[i].elm[ia]=0;
		}
    }
    
    *nelm=*nelm-iv; 
}

//�v�f�����m�F�֐�
void fill3D(vector<point3D> &NODE,vector<element3D> &ELEM,int nelm)
{
    for(int i=1;i<=nelm;i++)
    {
        int ielm=i;
		for(int j=1;j<=4;j++)
		{
			int ia=ELEM[ielm].node[j%4+1];
			int ib=ELEM[ielm].node[4-(j-1)/2*2];
			int ic=ELEM[ielm].node[3-j/2%2*2];
			int jelm=ELEM[ielm].elm[j];
			
			if(jelm!=0 && ielm<jelm)
			{
				int k=iface3D(ELEM,ielm,jelm);
				int ja=ELEM[jelm].node[k%4+1];
				int jb=ELEM[jelm].node[4-(k-1)/2*2];
				int jc=ELEM[jelm].node[3-k/2%2*2];
				int flag=0;
				if(ia==jc && ib==jb && ic==ja) flag=1;
				if(ib==jc && ic==jb && ia==ja) flag=1;
				if(ic==jc && ia==jb && ib==ja) flag=1;
				if(flag==0) cout<<"ERROR IN FILL"<<endl;
			}
		}
    }
}

//�v�f�ގ�����֐�
void set_material(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm)
{
    for(int i=1;i<=nelm;i++)
    {
        int M1=NODE[ELEM[i].node[1]].material;
		int M2=NODE[ELEM[i].node[2]].material;
		int M3=NODE[ELEM[i].node[3]].material;
		int M4=NODE[ELEM[i].node[4]].material;

		if(M1==HYPERELAST)cout<<"M1"<<" "<<"HYPERELAST"<<endl;
		if(M2==HYPERELAST)cout<<"M2"<<" "<<"HYPERELAST"<<endl;
		if(M3==HYPERELAST)cout<<"M3"<<" "<<"HYPERELAST"<<endl;
		if(M4==HYPERELAST)cout<<"M4"<<" "<<"HYPERELAST"<<endl;
	
		///4���_���ׂĂ������ގ��Ȃ�v�f������ɂȂ炤�B
		///�ЂƂł��قȂ��Ă������C�ƒ�`
		if(M1==M2 && M2==M3 && M3==M4) ELEM[i].material=M1;
		else if(M1!=AIR && M2!=AIR && M3!=AIR && M4!=AIR) 
		{
			if(M1!=ELASTIC && M2!=ELASTIC && M3!=ELASTIC && M4!=ELASTIC)
			{
				if(M1!=HYPERELAST && M2!=HYPERELAST && M3!=HYPERELAST && M4!=HYPERELAST)	ELEM[i].material=MAGELAST;	//15/2/10

/*				if(M1!=MAGELAST && M2!=MAGELAST && M3!=MAGELAST && M4!=MAGELAST)
				{
					if(M1!=HYPERELAST && M2!=HYPERELAST && M3!=HYPERELAST && M4!=HYPERELAST)	ELEM[i].material=MAGELAST;
					else if (true)
					{
						ELEM[i].material=HYPERELAST;
					}
				}*/

			}
			//ELEM[i].material=HYPERELAST��ǉ�15/2/4
			//if(CON.get_model_number()==15) ELEM[i].material=AIR;
			else ELEM[i].material=ELASTIC; //else ELEM[i].material=FLUID;
			if(M1==IRON || M2==IRON || M3==IRON || M4==IRON) ELEM[i].material=IRON;//�R�C���̗v�f�̓R�C���ړ_�̓����ɂȂ�悤�ɂ���
		}
		else ELEM[i].material=AIR;
    }
}

//���b�V���f�[�^�o�͊֐�
void data_avs(int node_number,int nelm,vector <point3D> &NODE,vector <element3D> &ELEM,int ktj,double *val,mpsconfig &CON)
{
	stringstream ss;
	ss<<"./Mesh"<<"/mesh1_"<<CON.get_current_step()<<".inp";
	string filename=ss.str();
	
	ofstream fq(filename);
	if(fq.fail()){
		system("mkdir Mesh");
		ofstream fq(filename);
		if(fq.fail()){
		cout<<"���b�V���t�@�C�����J���܂���"<<endl;
		exit(1);
		}
	}

    cout<<"Now writing-----";	
	
	//mesh1
	fq<<"# Micro AVS"<<endl;
	fq<<"1"<<endl;
	fq<<"data"<<endl;
	fq<<"step1"<<endl;
	//fq<<node_number+8<<" "<<nelm<<endl;//�������̓X�[�p�[�{�b�N�X���\������Ƃ�
	fq<<node_number<<" "<<nelm<<endl;

	/*for(i=ktj+1;i<=ktj+8;i++) //�������̓X�[�p�[�{�b�N�X���\������Ƃ�
	{
		fq<<i<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<endl;
		//fprintf(fp, "\t%d\t%f\t%f\t%f\n", i, NODE[i].r[A_X], NODE[i].r[A_Y], NODE[i].r[A_Z]);
		
	} */

	//�ߓ_�ԍ��Ƃ��̍��W�̏o��
	for(int i=1;i<=node_number;i++) fq<<i<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<endl;
	
	//�v�f�ԍ��Ɨv�f�`��̎�ށA�����ėv�f���\������ߓ_�ԍ��o��
	for(int i=1;i<=nelm;i++)
	{
		fq<<i<<" t0 tet ";
		for(int j=1;j<=4;j++) fq<<ELEM[i].node[j]<<" ";
		fq<<endl;
	}
	
	fq<<"1 0"<<endl;
	fq<<"1 1"<<endl;
	fq<<"element, e"<<endl;

    /*for(int i=ktj+1;i<=ktj+8;i++) //�������̓X�[�p�[�{�b�N�X���\������Ƃ�
	{
		fq<<i<<" 0"<<endl;
	}*/

	//�e�ߓ_�̒l���
	for(int i=1;i<=node_number;i++) fq<<i<<" "<<val[i]<<endl;		//�o�͂������X�J���[�l��val�ɓ��͂��Ă�������
	
	cout<<"OK"<<endl;
	fq.close();
}
//���b�V���f�[�^�o�͊֐�
void C_Fluix_data_avs(int node,int nelm,vector <point3D> &NODE,vector <element3D> &ELEM,int ktj,double **B,mpsconfig &CON,int t)
{
	stringstream ss;
	ss<<"./FluxContour"<<"/FluxC_"<<t<<".inp";
	string filename=ss.str();
	ofstream ffc(filename);
	if(ffc.fail()){
		system("mkdir FluxContour");
		ofstream ffc(filename);
		if(ffc.fail()){
		cout<<"./FluxContour�t�H���_���J���܂���ł���"<<endl;
		exit(1);
		}
	}
	//�Q�l�ɂ��Ă��鏑����microAVS�̃w���v�ł��Ȃ��̃f�[�^�́H���u��\���i�q�^�f�[�^�i�A�X�L�[�j�̏����v
	int step=CON.get_step()/(CON.get_EM_interval()*CON.get_mesh_output_interval())+1;				//�o�͂��鑍�X�e�b�v��
	int nelm2=0;			//�o�͂���v�f��
	int *node_flag=new int[node+1];			//�e�ߓ_���o�͂��邩�A���Ȃ���
	int node_on=0;							//�o�͂���ߓ_��

	for(int i=1;i<=nelm;i++) ELEM[i].map=0;//ϯ��ݸޏ�����(map�͗v�f�쐬�̍ۂɗ��p���Ă������A�����ł�������؂��g��)
	for(int i=0;i<=node;i++) node_flag[i]=OFF;		//������

	for(int i=1;i<=nelm;i++)
	{
		double Y=0;//�v�f�̏d�S��y���W
		for(int j=1;j<=4;j++) Y+=NODE[ELEM[i].node[j]].r[A_Y]*0.25;
		if(Y<0)
		{
			nelm2++;
			ELEM[i].map=1;//�o�͂���Ƃ�����
			for(int j=1;j<=4;j++) node_flag[ELEM[i].node[j]]=ON;//�֌W����ߓ_�̃t���O��ON
		}
	}
	///�o�͂���ߓ_�Ɨv�f�A����т��̐������Ƃ܂���

	for(int i=1;i<=node;i++) if(node_flag[i]=ON) node_on++;
	
		ofstream fout(filename);

		fout<<1<<endl;
		fout<<"data_geom"<<endl;
		fout<<"step1"<<endl;
		fout.close();

	ofstream fp(filename,ios :: app);
//	if(CON.get_EM_interval()==1) fp<<"step"<<t/(CON.get_EM_interval()*CON.get_mesh_output_interval())<<endl;
//	else if(CON.get_EM_interval()>1) fp<<"step"<<t/(CON.get_EM_interval()*CON.get_mesh_output_interval())+1<<endl;
	
	fp<<node_on<<" "<<nelm2<<endl;	//�ߓ_�͖��֌W�̂��̂��o�͂��Ă������ǁA���ꂾ�ƃt�@�C�����d���Ȃ邩��A�K�v�Ȑߓ_�����o�� 
	
	
	
	//�ߓ_�ԍ��Ƃ��̍��W�̏o��
	for(int i=1;i<=node;i++)//�ߓ_�͖��֌W�̂��̂��o�͂��Ă������ǁA���ꂾ�ƃt�@�C�����d���Ȃ邩��A�K�v�Ȑߓ_�����o�� 
	{
		//fp<<i<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<endl;
		if(node_flag[i]==ON) fp<<i<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<endl;
	}
	
	
	//�v�f�ԍ��Ɨv�f�`��̎�ށA�����ėv�f���\������ߓ_�ԍ��o��
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].map==1)
		{
			fp<<i<<"  0 tet ";
			for(int j=1;j<=4;j++)	fp<<ELEM[i].node[j]<<" ";
			fp<<endl;
		}
	}

	fp<<"0 1"<<endl;//�ߓ_�̏��ʂ��[���ŁA�v�f�̏��ʂ�num_info�Ƃ������ƁB
	fp<<"1 1"<<endl;
	fp<<"material, material"<<endl;

	//�e�ߓ_�̒l���
	//for(int i=1;i<=node_number;i++) fp<<i<<" "<<NODE[i].material<<endl;
	//for(int i=1;i<=node;i++) fp<<i<<" "<<val[i]<<endl;

	for(int i=1;i<=nelm;i++) if(ELEM[i].map==1) fp<<i<<"  "<<sqrt(B[A_X][i]*B[A_X][i]+B[A_Y][i]*B[A_Y][i]+B[A_Z][i]*B[A_Z][i])<<endl;
	
	cout<<"OK"<<endl;
	fp.close();
	delete [] node_flag;
}
//���b�V���f�[�^�o�͊֐�ver.2 ����f�ʂ̃��b�V�����݂����Ƃ��Ɏg��
void data_avs2(mpsconfig &CON,int node,int nelm,vector <point3D> &NODE,vector <element3D> &ELEM,int ktj,int t)
{
	//�Q�l�ɂ��Ă��鏑����microAVS�̃w���v�ł��Ȃ��̃f�[�^�́H���u��\���i�q�^�f�[�^�i�A�X�L�[�j�̏����v
	cout<<"Now ver.2 writing-----";
	int step=CON.get_step()/(CON.get_EM_interval()*CON.get_mesh_output_interval())+1;				//�o�͂��鑍�X�e�b�v��
	int nelm2=0;			//�o�͂���v�f��
	int *node_flag=new int[node+1];			//�e�ߓ_���o�͂��邩�A���Ȃ���
	int node_on=0;							//�o�͂���ߓ_��
	
	for(int i=1;i<=nelm;i++) ELEM[i].map=0;//ϯ��ݸޏ�����(map�͗v�f�쐬�̍ۂɗ��p���Ă������A�����ł�������؂��g��)
	for(int i=0;i<=node;i++) node_flag[i]=OFF;		//������

	for(int i=1;i<=nelm;i++)
	{
		double Y=0;//�v�f�̏d�S��y���W
		for(int j=1;j<=4;j++) Y+=NODE[ELEM[i].node[j]].r[A_Y]*0.25;
		if(Y<0)
		{
			nelm2++;
			ELEM[i].map=1;//�o�͂���Ƃ�����
			for(int j=1;j<=4;j++) node_flag[ELEM[i].node[j]]=ON;//�֌W����ߓ_�̃t���O��ON
		}
	}
	///�o�͂���ߓ_�Ɨv�f�A����т��̐������Ƃ܂���

	for(int i=1;i<=node;i++) if(node_flag[i]=ON) node_on++;

	
		stringstream ss;
		ss<<"./Mesh"<<"/mesh2_"<<CON.get_current_step()<<".inp";
		string filename=ss.str();
	
		ofstream fp(filename);
		if(fp.fail()){
		system("mkdir Mesh");
		ofstream fp(filename);
		if(fp.fail()){
		cout<<"���b�V���t�@�C�����J���܂���"<<endl;
		exit(1);
		}
	}

		fp<<step<<endl;
		fp<<"data_geom"<<endl;
		
	
	

	
//	if(CON.get_EM_interval()==1) fp<<"step"<<t/(CON.get_EM_interval()*CON.get_mesh_output_interval())<<endl;
//	else if(CON.get_EM_interval()>1) fp<<"step"<<t/(CON.get_EM_interval()*CON.get_mesh_output_interval())+1<<endl;
	fp<<"step"<<1<<endl;
	fp<<node_on<<" "<<nelm2<<endl;	//�ߓ_�͖��֌W�̂��̂��o�͂��Ă������ǁA���ꂾ�ƃt�@�C�����d���Ȃ邩��A�K�v�Ȑߓ_�����o�� 
	
	
	
	//�ߓ_�ԍ��Ƃ��̍��W�̏o��
	for(int i=1;i<=node;i++)//�ߓ_�͖��֌W�̂��̂��o�͂��Ă������ǁA���ꂾ�ƃt�@�C�����d���Ȃ邩��A�K�v�Ȑߓ_�����o�� 
	{
		//fp<<i<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<endl;
		if(node_flag[i]==ON) fp<<i<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<endl;
	}
	
	
	//�v�f�ԍ��Ɨv�f�`��̎�ށA�����ėv�f���\������ߓ_�ԍ��o��
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].map==1)
		{
			fp<<i<<"  0 tet ";
			for(int j=1;j<=4;j++)	fp<<ELEM[i].node[j]<<" ";
			fp<<endl;
		}
	}

	fp<<"0 1"<<endl;//�ߓ_�̏��ʂ��[���ŁA�v�f�̏��ʂ�num_info�Ƃ������ƁB
	fp<<"1 1"<<endl;
	fp<<"material, material"<<endl;

	//�e�ߓ_�̒l���
	//for(int i=1;i<=node_number;i++) fp<<i<<" "<<NODE[i].material<<endl;
	//for(int i=1;i<=node;i++) fp<<i<<" "<<val[i]<<endl;
	int mate=0;
	for(int i=1;i<=nelm;i++){
		if(ELEM[i].map==1){ 
			if(ELEM[i].material==AIR) mate=0;	//��
			else if(ELEM[i].material==MAGELAST) mate=3;	//��
			else if(ELEM[i].material==COIL) mate=4;	//��
			else if(ELEM[i].material==IRON) mate=5;	//��
			else if(ELEM[i].material==MAGELAST2) mate=6;
//			else if(ELEM[i].material==HYPERELAST) mate=3;//HYPERELAST�ǉ�15/2/4
			fp<<i<<"  "<<mate<<endl;
		}
	}
	
	cout<<"OK"<<endl;
	fp.close();
	delete [] node_flag;
}

//���b�V���f�[�^�o�͊֐�(�ގ�ver�j
void data_avs3(int node,int nelm,vector <point3D> &NODE,vector <element3D> &ELEM,mpsconfig &CON)
{
	stringstream ss;
	ss<<"./Mesh"<<"/mesh3_"<<CON.get_current_step()<<".mgf";
	string filename=ss.str();
	
	ofstream fout(filename);
	if(fout.fail()){
		system("mkdir Mesh");
		ofstream fout(filename);
		if(fout.fail()){
			cout<<filename<<"���J���܂���B�t�H���_���m�F���Ă�������"<<endl;
			exit(1);
		}
	}

	cout<<"Now ver.3 writing-----";
	int step=1;
	fout <<"# Micro AVS Geom:2.00" << endl;
	//fout << "step" << step<< endl;
	fout << "disjoint polygon" << endl;
	fout <<"element" << endl;
	fout <<"facet" << endl;
	fout <<"color" << endl;

	int count=0;//�����\�ʐ�
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material!=AIR)
		{
			for(int j=1;j<=4;j++)
			{
				int jelm=ELEM[i].elm[j];
				if(ELEM[jelm].material==AIR) count++;
			}
		}
	}///�\�ʐ������Ƃ܂���
	
	fout <<count<< endl;//�v�f���o��
	double red, green, blue;
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material!=AIR)
		{
			for(int j=1;j<=4;j++)
			{
				int jelm=ELEM[i].elm[j];
				if(ELEM[jelm].material==AIR)
				{
					if(ELEM[i].material==FLUID ||  ELEM[i].material==MAGELAST || ELEM[i].material==MAGELAST2 ||ELEM[i].material==HYPERELAST)
					{
						red   = 0.0;
						green = 1.0;
						blue  = 0.0;
					}
					else if(ELEM[i].material==ELECTRODE || ELEM[i].material==MAGNET || ELEM[i].material==COIL)
					{
						red   = 1.0;
						green = 0.0;
						blue  = 0.0;
					}
					else if(ELEM[i].material==ELASTIC)
					{
						red   = 0.0;
						green = 0.0;
						blue  = 1.0;
					}
					else if(ELEM[i].material==WALL)
					{
						red=0.0;
						green=1.0;
						blue=1.0;
					}
					
					fout <<"3"<< endl;
					int ia=ELEM[i].node[j%4+1];//���_���̌����Ɉʒu����O�p�`�̒��_
					int ib=ELEM[i].node[4-(j-1)/2*2];
					int ic=ELEM[i].node[3-(j/2%2)*2];
					for(int D=0;D<3;D++) fout << NODE[ia].r[D]<<" ";
					fout << red << " " << green << " " << blue;
					fout << endl;

					for(int D=0;D<3;D++) fout << NODE[ib].r[D]<<" ";
					fout << red << " " << green << " " << blue;
					fout << endl;

					for(int D=0;D<3;D++) fout << NODE[ic].r[D]<<" ";
					fout << red << " " << green << " " << blue;
					fout << endl;
				}
			}
		}
	}////*/

	fout.close();
	cout<<"OK"<<endl;
}

//���b�V���f�[�^�o�͊֐�ver.4 �������̃��b�V�����݂����Ƃ��Ɏg��
void data_avs4(mpsconfig &CON,int node_number,vector <point3D> &NODE,vector <element3D> &ELEM,int nelm,int *kv)
{
	//nelm:�o�͂���v�f��
	//kv[]:�o�͂���v�f�ԍ�
	ofstream fp("mesh4.inp");
	
	fp<<"# Micro AVS"<<endl;
	fp<<"1"<<endl;
	fp<<"data"<<endl;
	fp<<"step1"<<endl;
	fp<<node_number<<" "<<nelm<<endl;
	
	//�ߓ_�ԍ��Ƃ��̍��W�̏o��
	for(int i=1;i<=node_number;i++)//�ߓ_�͖��֌W�̂��̂��o�͂���΂��� 
	{
		fp<<i<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<endl;
	}
	
	
	//�v�f�ԍ��Ɨv�f�`��̎�ށA�����ėv�f���\������ߓ_�ԍ��o��
	for(int k=1;k<=nelm;k++)
	{
		int i=kv[k];
		
		fp<<i<<"  0 tet ";
		for(int j=1;j<=4;j++)	fp<<ELEM[i].node[j]<<" ";
		fp<<endl;
	}
	
	fp<<"1 0"<<endl;
	fp<<"1 1"<<endl;
	fp<<"element, e"<<endl;

	//�e�ߓ_�̒l���
	for(int i=1;i<=node_number;i++) fp<<i<<" "<<NODE[i].material<<endl;
	//for(int i=1;i<=node_number;i++) fp<<i<<" "<<val[i]<<endl;
	
	fp.close();


	ofstream fout("test.mgf");

	fout<<"# Micro AVS Geom:2.00"<<endl;
	fout<<"1"<<endl;//microAVS�ɏo�͂��鑍�ï�ߐ��B�t�@�C���o�͂�CON.get_interval()���1��ƍŏ��ɍs���B
		
	fout<<"step"<<1<<endl;
	fout<<"sphere"<<endl;
	fout<<"time="<<1<<endl;
	fout<<"color"<<endl;

	double red,green,blue;	//���q�̐F��\������3���F
	
	int num=2;		//�o�͐ߓ_��
		
	fout<<num<<endl;
	red=0.5;green=0.5;blue=0;
		
	fout<<NODE[node_number].r[A_X]<<" "<<NODE[node_number].r[A_Y]<<" "<<NODE[node_number].r[A_Z]<<" ";//���W�o��
			
	fout<<CON.get_distancebp()<<" "<<red<<" "<<green<<" "<<blue<<endl;//�F�o��

	red=1;green=0;blue=0;
		
	fout<<NODE[7883].r[A_X]<<" "<<NODE[7883].r[A_Y]<<" "<<NODE[7883].r[A_Z]<<" ";//���W�o��

	fout<<CON.get_distancebp()<<" "<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
	
	fout.close();
}

///�������b�V�������֐�
void FINE3D(vector<point3D> &NODE,vector<element3D> &ELEM,int KTJ,int KTE,int *node,int *nelm,mpsconfig &CON,double rrm,int startID)
{
	
	cout<<"�v�f�̍ĕ������s �ߓ_��:"<<*node<<" �v�f��:"<<*nelm<<endl;

	unsigned int timeA=GetTickCount();	//�v�Z�J�n����

	/////�v�f�̑̐ςƃ{���m�C�_�����Ƃ߂Ă���(�X�[�p�[�{�b�N�X�̂��ƂŊ֐����Ăяo����Ă��Ή��ł���悤��)
	for(int i=1;i<=*nelm;i++)
    {   
		int ia=ELEM[i].node[1];
		int ib=ELEM[i].node[2];
		int ic=ELEM[i].node[3];
		int ip=ELEM[i].node[4];
		
		ELEM[i].volume=volume3D(NODE,ia,ib,ic,ip);//�̐ς�6�{�ł��邱�Ƃɒ���
	
		sphere3D(NODE,ELEM,ia,ib,ic,ip,i);//�O�ڋ��̒��S�i�{���m�C�_)�Ɣ��a�̓����v�Z
    }
	//////////////*/

	//��C�w����
	int layer_node_num=0;
	layer_node_num=make_air_layer(NODE,ELEM,nelm,CON,node, KTE,rrm,KTJ);	//layer_node_num�͂��̊֐��ɂ�葝�������ߓ_��

	int Mnum=10;								//�c���Ɣ��肳��A�ĕ������ꂽ�v�f��
	int newnode=0;								//�V�������������ߓ_��
	int limit_num=CON.get_add_points()-layer_node_num;//�����ł���ő�ߓ_��
	//if(limit_num<0) limit_num=0;
	int *kv=new int[KTE];						//�V�ߓ_���O�ڋ��ɂӂ��ޗv�f�Q
    int *istack=new int[KTE];					//�ꎞ�z��

	int MESH;
	double limit=CON.get_co_fine();	//���ް��̏ꍇ�͍ŏ��ӂƍő�ӂ̔䗦�̂������l
	int remesh;
	
	if(CON.get_fine()==1 && limit_num>0)//���ް�
	{
		while(Mnum>0)
		{
			Mnum=0;			
		
			///������
			for(int i=1;i<KTE;i++) ELEM[i].map=0;
	
			///�����ߓ_�𓱓����Ă���
	    
			MESH=*nelm;
		
			for(int je=startID;je<=MESH;je++)
			{
				//cout<<je<<endl;
				if(newnode<limit_num)
				{
					//if(ELEM[je].material==AIR && depth[je]>=1)
					if(ELEM[je].material==AIR)
					{
						int flag3=0;//1�Ȃ�V�ߓ_�𓱓�����ӏ������܂����Ƃ�������
						double r[3];//�V�ߓ_�̍��W�i�[
						for(int j=1;j<=4;j++)
						{
							if(flag3==0 && ELEM[je].elm[j]!=0)
							{
								int ia=ELEM[je].node[j%4+1];//ielm��jelm�̐ڂ���O�p�`���\������ߓ_�ԍ�  ia,ib,ic,ip�𒸓_�Ƃ��ĐV�����v�f�������
								int ib=ELEM[je].node[4-(j-1)/2*2];
								int ic=ELEM[je].node[3-(j/2%2)*2];
		
								double iaib=0;//��ia,ib�̋���
								double iaic=0;//��ia,ic�̋���
								double ibic=0;//��ib,ic�̋���
								for(int D=0;D<3;D++)
								{
									iaib+=(NODE[ib].r[D]-NODE[ia].r[D])*(NODE[ib].r[D]-NODE[ia].r[D]);
									iaic+=(NODE[ia].r[D]-NODE[ic].r[D])*(NODE[ia].r[D]-NODE[ic].r[D]);
									ibic+=(NODE[ib].r[D]-NODE[ic].r[D])*(NODE[ib].r[D]-NODE[ic].r[D]);
								}
								iaib=sqrt(iaib);
								iaic=sqrt(iaic);
								ibic=sqrt(ibic);
								double minL=iaib;//�ŏ��Ӓ���
								double maxL=iaib;//�ő�Ӓ���
								int N[2]={ia,ib};//�ő�ӂ��\������ߓ_�ԍ�
								if(iaic<minL) minL=iaic;
								else if(iaic>maxL) 
								{
									maxL=iaic;
									N[0]=ia;
									N[1]=ic;
								}
								if(ibic<minL) minL=ibic;
								else if(ibic>maxL)
								{
									maxL=ibic;
									N[0]=ib;
									N[1]=ic;
								}
								if(maxL>limit*minL) //�ő咷�����ŏ���limit�{�ȏ゠���
								{
									int n1=N[0];//�ő�ӂ��\������ߓ_�ԍ�
									int n2=N[1];
									if(NODE[n1].boundary_condition==0 && NODE[n2].boundary_condition==0)//���m���Ȃ�
									{
										flag3=1;
										for(int D=0;D<3;D++) r[D]=0.5*(NODE[n1].r[D]+NODE[n2].r[D]);	//�V�_�͍ő�ӂ̒��_�ɐݒu
										
										if(NODE[n1].remesh==NODE[n2].remesh) remesh=NODE[n1].remesh;	//remesh�������Ȃ炻��������p���B�قȂ�Ȃ炻��͔�remesh�̈���̐ߓ_�Ȃ̂�remesh=OFF
										else remesh=OFF;
									}
								}
							}
						}
						if(flag3==1)
						{
							Mnum++;
							newnode++;		//�V�_������
							*node=*node+1;
							
							int ip=*node;
							for(int D=0;D<3;D++) NODE[ip].r[D]=r[D];
							NODE[ip].boundary_condition=0;
							NODE[ip].material=AIR;
							NODE[ip].remesh=remesh;
							NODE[ip].particleID=-1;		//�Ή����闱�q�͑��݂��Ȃ�
							double xp=NODE[ip].r[A_X];//��������ߓ_�̍��W
							double yp=NODE[ip].r[A_Y];
							double zp=NODE[ip].r[A_Z];
			
							///�V�ߓ_���܂ޗv�f�̒T��
							
							int loc=je;//���炩��je�͐V�_���O�ڋ��Ɋ܂�
							
							//////////�O�ڋ����ɐV�ߓ_���܂ޗv�f�̒��o
							int iv=0;
							int msk=0;
			
							iv++;//�O�ڋ����ɐV�ߓ_���܂ޗv�f��
							kv[iv]=loc;
							ELEM[loc].map=1;//map��1�̗v�f�́A�O�ڋ��ɐߓ_i���܂ނƂ�������
							msk++;
							istack[msk]=loc;
			
							while(msk!=0)
							{   
								int isk=istack[msk];//���ܒ��ڂ��Ă���v�f�̔ԍ�
								msk--;
								for(int j=1;j<=4;j++)
								{
									int jelm=ELEM[isk].elm[j];//isk�Ɛڂ���v�f
									
									if(jelm!=0)//���ꂪ�\�ʂłȂ��Ȃ�
									{
										if(ELEM[jelm].map==0) //�܂��������ĂȂ��Ȃ�
										{   
											double rad=ELEM[jelm].RR*(1.000000+ERR);//�O�ڋ����a�̂Q��
											double dst=ELEM[jelm].r[A_X]*ELEM[jelm].r[A_X]-2.000000*ELEM[jelm].r[A_X]*xp+xp*xp;///�O�ڋ����S�ƐV�ߓ_�̋���
											if(dst<rad)///������dst>rad�Ȃ��΂ɊO�ڋ��Ɋ܂܂Ȃ�
											{
												dst+=ELEM[jelm].r[A_Y]*ELEM[jelm].r[A_Y]-2.000000*ELEM[jelm].r[A_Y]*yp+yp*yp;
												if(dst<rad)
												{
													dst+=ELEM[jelm].r[A_Z]*ELEM[jelm].r[A_Z]-2.000000*ELEM[jelm].r[A_Z]*zp+zp*zp;
													if(dst<rad)//�O�ڋ����Ɋ܂�
													{
														if(ELEM[jelm].material==AIR )	//���͉̂󂵂Ăق����Ȃ����炱���Ɋ܂߂Ȃ�
														{
															iv++;//�O�ڋ����ɐV�ߓ_���܂ޗv�f�����{�P
															kv[iv]=jelm;//���X�g�ɂ����
															msk++;
															istack[msk]=jelm;
															ELEM[jelm].map=1;//jelm�͊O�ڋ����ɐV�ߓ_���܂�
														}
													}
												}
											}
										}
									}
								}
							}//�V�_���O�ڋ����Ɋ܂ޗv�f��iv�Ƃ��̗v�f�ԍ�kv[iv]�����Ƃ܂���

							/////////////////
							int NN=*nelm;//�����O�̗v�f��
		
							////����ꂽ���ʑ̂��l�ʑ̂ɕ������A�V�����v�f�𐶐�����B�܂��A�֐����ŐV�v�f�̍ގ������肷��
							int FLAG=poly3D_for_FINE3D(NODE,ELEM,&iv,kv,ip,nelm,CON); 
							if(FLAG==OFF)		//poly�����s�����ꍇ
							{
								*node=*node-1;	//�V�ߓ_������������߂�
								Mnum--;
								newnode--;
							}
						}
					}
				}
				else Mnum=0;//newnode��limit_num�ɓ������Ȃ�����Mnum=0�ɂ���ٰ�߂���E�o
			}
		}
	}
	cout<<"�ߓ_+="<<newnode<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
	

	delete [] kv;
	delete [] istack;
}

//FINE�p���ʑ̕����֐�
int poly3D_for_FINE3D(vector<point3D> &NODE,vector<element3D> &ELEM,int *iv,int *kv,int ip,int *nelm,mpsconfig &CON)
{   

    ////���̊֐��̑O�̒i�K�ŁA�V�ߓ_���O�ڋ��Ɋ܂ގl�ʑ̐�iv�Ƃ��̗v�f�ԍ�kv[iv]�����Ƃ܂��Ă���B
    ////�܂��A�V�ߓ_���O�ڋ��Ɋ܂ގl�ʑ̂�ELEM[i].map=1�ƂȂ��Ă���

    int ix=0;				//�\�ʂ̐�. ��ʂ�iv=ix�Ƃ͂Ȃ�Ȃ�.�ЂƂ̗v�f�������̕\�ʂ�S���̂ŁB(iv<=ix�ł���)
    int imen[10000][3+1];	//���ʑ̕\�ʎO�p�`�̐ߓ_�ԍ��i�[
    int jmen[10000];		//���ʑ̕\�ʎO�p�`�ɗאڂ���l�ʑ̔ԍ��i�[
    int kmen[10000];		//���ʑ̕\�ʎO�p�`�ɗאڂ���l�ʑ̗̂אږʔԍ� (����͑扽�ʂŎ����Ɛڂ��Ă��邩)
    double vol[10000];		//���ʑ̂̑̐ς̂U�{
    
    ///////���ʑ̕\�ʐ�ix�Ƃ�����\�����钸�_�����Ƃ߂�
    int flag=0;
    while(flag==0)
    {   
        ix=0;
        flag=1;
        for(int i=1;i<=*iv;i++)//*iv�͐V�_���O�ڋ����Ɋ܂ޗv�f�̐�
        {
			if(flag==1)
			{
        		int ielm=kv[i];//�V�_���O�ڋ����Ɋ܂ޗv�f�ԍ�
				for(int j=1;j<=4;j++)
				{
					if(flag==1)
					{
						int jelm=ELEM[ielm].elm[j];		//ielm�v�f�ɗאڂ���v�f�ԍ�
						int ia=ELEM[ielm].node[j%4+1];	//ielm��jelm�̐ڂ���O�p�`���\������ߓ_�ԍ�  ia,ib,ic,ip�𒸓_�Ƃ��ĐV�����v�f�������
						int ib=ELEM[ielm].node[4-(j-1)/2*2];
						int ic=ELEM[ielm].node[3-(j/2%2)*2];
							
						if(jelm==0)//�\�ʂȂ�ielm�͕\�ʗv�f�ł���Ƃ킩��
						{
							ix++;
							imen[ix][1]=ia;	//���ʑ̕\�ʎO�p�`�̒��_�ԍ�
							imen[ix][2]=ib;
							imen[ix][3]=ic;
							jmen[ix]=0;		//���ʑ̕\�ʎO�p�`�ɗאڂ���v�f�ԍ�
							kmen[ix]=0;		//���ʑ̕\�ʎO�p�`�ɗאڂ���v�f�̗אږʔԍ�  �����0������A�����ł�0����
		
							vol[ix]=volume3D(NODE,ia,ib,ic,ip);

							if(vol[ix]<ERR)//���̍s���������ق������܂������H�C�̂����H
							{
								
								for(int ii=1;ii<=*iv;ii++) ELEM[kv[ii]].map=0;//�����ŏ��������Ă����Ȃ��Ǝ��̐ߓ_�̂Ƃ��o�O��
								return OFF;//�V�ߓ_������������߂�
								
							}
						}
						else if(ELEM[jelm].map==0)//ielm�ɐڂ���jelm�͐V�_���O�ڋ��Ɋ܂܂Ȃ��B�܂�ielm��jelm�̋��E�͑��ʑ̕\�ʂƂ�������
						{
							ix++;
							imen[ix][1]=ia;
							imen[ix][2]=ib;
							imen[ix][3]=ic;
							jmen[ix]=jelm;
							kmen[ix]=iface3D(ELEM,ielm,jelm);//jelm��ielm�ɐڂ���ʔԍ�
							if(ix>=10000)cout<<"ix>10000"<<endl;

//2012-09-28�ύX			if(kmen[ix]==-1) cout<<"error ielm/jelm="<<ielm<<"/"<<jelm<<"  ELEM[jelm].elm="<<ELEM[jelm].elm[1]<<","<<ELEM[jelm].elm[2]<<","<<ELEM[jelm].elm[3]<<","<<ELEM[jelm].elm[4]<<endl;

							if(kmen[ix]==-1)
							{	//�Ȃ��iface�G���[�ɂȂ邩�悭�킩��Ȃ�
								for(int ii=1;ii<=*iv;ii++) ELEM[kv[ii]].map=0;//�����ŏ��������Ă����Ȃ��Ǝ��̐ߓ_�̂Ƃ��o�O��
								return OFF;//�V�ߓ_������������߂�
							}
							vol[ix]=volume3D(NODE,ia,ib,ic,ip);
							//if(ELEM[jelm].material!=AIR || depth[jelm]==1)//���͔̂j�󂵂Ăق����Ȃ�+�[��1������B���������ꂷ��Ƃ��܂�iface���G���[�ɂȂ�B�Ȃ�ŁH
							
							//�[��1�̋�C�������Ă悢�Ƃ��͉�������(nanoe�͂n�m ferrofluid�͂n�e�e�j
							//if(CON.get_model_number()!=15)
							{
								if(ELEM[jelm].material==COIL || ELEM[jelm].material==IRON)//���͔̂j�󂵂Ăق����Ȃ�����
								{	
									for(int ii=1;ii<=*iv;ii++) ELEM[kv[ii]].map=0;//�����ŏ��������Ă����Ȃ��Ǝ��̐ߓ_�̂Ƃ��o�O��
									return OFF;//�V�ߓ_������������߂�
								}
							}//////*/

							/*if(ELEM[jelm].material==FLUID)//���͔̂j�󂵂Ăق����Ȃ�����
							{	//�Ȃ��iface�G���[�ɂȂ邩�悭�킩��Ȃ�
								for(int ii=1;ii<=*iv;ii++) ELEM[kv[ii]].map=0;//�����ŏ��������Ă����Ȃ��Ǝ��̐ߓ_�̂Ƃ��o�O��
								return OFF;//�V�ߓ_������������߂�
							}///*/

							if(ELEM[jelm].material==MAGNET)//���͔̂j�󂵂Ăق����Ȃ�����
							{	//�Ȃ��iface�G���[�ɂȂ邩�悭�킩��Ȃ�
								for(int ii=1;ii<=*iv;ii++) ELEM[kv[ii]].map=0;//�����ŏ��������Ă����Ȃ��Ǝ��̐ߓ_�̂Ƃ��o�O��
								return OFF;//�V�ߓ_������������߂�
							}//////*/
							if(vol[ix]<=ERR)//���ȏ��}3.10�̂悤�ɁA�s�K�؂Ȏl�ʑ̂�z�肵�Ă��܂��A���ʑ̐ς����ɂȂ�B���̏ꍇ�A���ȏ��ɏ����Ă���Ƃ���A�ŏ�������Ȃ���(goto 10)
							{   
								if(ELEM[jelm].material==AIR)
								{
									*iv=*iv+1;//���X�g�𑝂₷ *iv++�͂��߁H
								
									kv[*iv]=jelm;//�{���͊O�ڋ����ɐV�ߓ_���܂܂Ȃ����A���ʑ̕������\�ɂ��邽�߂Ƀ��X�g�ɂ����
									ELEM[jelm].map=1;
									flag=0;
								}
								else	//���͉̂󂵂Ăق����Ȃ����烊�X�g�ɂ���Ȃ�
								{
									for(int ii=1;ii<=*iv;ii++) ELEM[kv[ii]].map=0;//�����ŏ��������Ă����Ȃ��Ǝ��̐ߓ_�̂Ƃ��o�O��
									return OFF;//�V�ߓ_������������߂�
								}
							}
						}
					}
				}
			}
        }
    }
    //////���ʑ̕\�ʐ�ix�Ƃ�����\�����钸�_�����Ƃ܂���
   
    /////////////�\�ʂ̒��_�ƐV�ߓ_���Ȃ���.����Ƒ��ʑ̂���ix�̗v�f�����������

    int ibound=ix;//ix�̑���B�܂�\�ʂ̐���\��
    
	///�̐ς�0�ɂȂ�Ȃ��߂�
	for(int i=1;i<=ibound;i++)
	{
		if(vol[i]<=0)
		{
			if(vol[i]==0) cout<<"vol=0"<<" ip="<<ip<<" i="<<i<<endl;
			else if(vol[i]<0) cout<<"vol<=0"<<endl;
			//�Ȃ��0�ɂȂ邩�悭�킩��Ȃ�
			for(int ii=1;ii<=*iv;ii++) ELEM[kv[ii]].map=0;//�����ŏ��������Ă����Ȃ��Ǝ��̐ߓ_�̂Ƃ��o�O��
			return OFF;//�V�ߓ_������������߂�
		}
	}////

	int nelm0=*nelm;//�ύX�O�̗v�f�����L��

    for(int i=*(iv)+1;i<=ibound;i++) ///iv�̗v�f���j�󂳂�ĐV����ibound�̗v�f����������邩��A������v�f����ibound-iv
    {   
        *nelm=*nelm+1;		//*nelm++�Ƃ����������ł͂��߁H
        kv[i]=*nelm;		//�V�ߓ_���܂ޗv�f���X�g��������B
		ELEM[*nelm].map=1;	//���������v�f�͕K���V�ߓ_���O�ڋ��Ɋ܂�(�Ƃ��Ă���).(���̍s�Ӗ��Ȃ��Ȃ��H)
		
    }
    
    for(int i=1;i<=ibound;i++) ELEM[kv[i]].map=0;//�}�b�s���O�̏�����

	/*///�V����ibound�̗v�f�̐[�������Ƃ߂�
	int aveDEP=0;//���ϐ[��
	for(int i=1;i<=ibound;i++)
	{
		aveDEP+=depth[kv[i]];
	}
	aveDEP/=ibound;		*/							////���ϐ[�������Ƃ܂���

    for(int i=1;i<=ibound;i++)//�v�f��񐶐�
    {   
        int ielm=kv[i];
	//	double determ=vol[i];
		int ia=imen[i][1];
		int ib=imen[i][2];
		int ic=imen[i][3];
		ELEM[ielm].node[1]=ia;
		ELEM[ielm].node[2]=ib;
		ELEM[ielm].node[3]=ic;
		ELEM[ielm].node[4]=ip;//�V�_�͂S�Ԗڂƒ�`	
		ELEM[ielm].elm[4]=jmen[i];
		if(jmen[i]!=0) ELEM[jmen[i]].elm[kmen[i]]=ielm;
		ELEM[ielm].volume=vol[i];
		
	
		sphere3D(NODE,ELEM,ia,ib,ic,ip,ielm);//�O�ڋ��̒��S�i�{���m�C�_)�Ɣ��a�̓����v�Z
		
		ELEM[ielm].material=AIR;//�ʏ��poly�֐��Ƃ̑���_�B�����ōގ�����C�ƌ��肷��

		//depth[ielm]=aveDEP;//���̐[���Ƃ��āA���ϐ[������
    }
    ///////////////////
    
    //�v�f-�v�f�֌W�C��/////////��̏����ő�4�ʂŐڂ���v�f�ԍ��͂킩���Ă���̂ŁA�c������߂�
	//						�����ŁA1�`3�ʂ͑��ʑ̂��\������v�f�Ƃ̋��E�ʂł��邱�Ƃɒ���
    ix=0;
    for(int i=1;i<=ibound;i++)
    {
        int ielm=kv[i];
		for(int j=1;j<=3;j++)//ELEM[ielm].elm[4]�͂��łɂ��Ƃ܂�������A����ȊO�����Ƃ߂�
		{
			///ELEM[ielm].node[4]=ip�ł���
			int ia=ELEM[ielm].node[j%3+1];		//j=1,2,3�̂Ƃ��A2,3,1�̏�
			int ib=ELEM[ielm].node[(j%3+1)%3+1];//j=1,2,3�̂Ƃ��A3,1,2�̏�
			int flag=0;
			for(int k=1;k<=ix;k++)
			{
				if(flag==0)
				{
					int ja=imen[k][1];
					int jb=imen[k][2];
					if(ia==ja && ib==jb)//�ߓ_����v������
					{
						ELEM[ielm].elm[j]=jmen[k];//���炩����ؽĂ��Ă����������i�[
						ELEM[jmen[k]].elm[kmen[k]]=ielm;
						imen[k][1]=imen[ix][1];		//k�Ԗڂ̏��͂����s�v�B�Ȃ̂Ŕz��̈�ԍŌ�̏���k�Ԗڂɂ����Ă��āA����܂ł̏��͔j������
						imen[k][2]=imen[ix][2];
						jmen[k]=jmen[ix];
						kmen[k]=kmen[ix];
						ix--;						//�҂��Ӑ�����
						flag=1;						//ELEM[ielm].elm[j]�͂��Ƃ܂����̂ŁA���̃l�X�g�ɓ���K�v�͂Ȃ��̂�flag=1
					}
				}
			}
			if(flag==0)
			{
				ix++;			//�����ł�ix�́A[�אڊ֌W�𖞂����v�f]���܂��Ă���[��]�̐���\���B
				imen[ix][1]=ib;	//�����̐ߓ_�̕��т��L�������A�ʂ̗v�f�����̕��т𖞂����̂�҂Bib��ia�̕��т��t�ɂ��Ă��邱�Ƃɒ���
				imen[ix][2]=ia;
				jmen[ix]=ielm;
				kmen[ix]=j;
			}
		}
    }///�v�f-�v�f�֌W�C������

	/*/////////////////////////////////////////////////�V�v�f�̐[������(�[���ω����A���ɂȂ�悤�ɁA��������܂Ōv�Z���J��Ԃ�)
	int count2=10;	//�[�����C�����ꂽ�v�f��
	while(count2>0)
	{
		count2=0;	//������
		for(int i=1;i<=ibound;i++)
		{
		    int ielm=kv[i];
			int flag=0;		//���ꂪ1�Ȃ畨�̐ߓ_�������Ă���Ƃ������ƁB���̏ꍇ�͐[����1�ɐݒ肷��
			for(int j=1;j<=4;j++)
			{
				int ia=ELEM[ielm].node[j];
				if(NODE[ia].material!=AIR)
				{
					flag=1;			//���E�v�f�ł���Ƃ킩��
					if(depth[ielm]!=1) count2++;
					depth[ielm]=1;	//���E�v�f�͐[��1�ɐݒ�
				}
			}
			if(flag==0)						//���E�v�f�łȂ��̂Ȃ�
			{
				int smallDEP=depth[ielm];	//�אڂ���v�f�̂Ȃ��ł̍ŏ��[�������Ƃ߂�
				for(int j=1;j<=4;j++)		//ELEM[ielm].elm[4]�͑��ʑ̂̊O�̗v�f�ł���
				{
					int jelm=ELEM[ielm].elm[j];
					if(jelm!=0) if(depth[jelm]<smallDEP) smallDEP=depth[jelm];
				}
				///�v�f�̐[�����A���͂̍ŏ��[��+1�ƒ�`�B���̂ɗאڂ���ꍇ�͉��̎����[����1�ƂȂ�
				if(depth[ielm]!=smallDEP+1) count2++;
				depth[ielm]=smallDEP+1;//���ʑ̓��v�f��depth�������؂�ɂ��Ƃ߂đ��v���H(���̂ɗאڂ���ꍇ�͎����[����1�ƂȂ邩��)
			}
		}
	}

	//���ʑ̂̂ЂƂO���̗v�f�̐[�����C��
	for(int i=1;i<=ibound;i++)
	{
		int ielm=kv[i];
		int jelm=ELEM[ielm].elm[4];//ELEM[ielm].elm[4]�͑��ʑ̂̊O�̗v�f�ł���
		///jelm��ielm���[���ꍇ�A�[����ielm���1�����傫������B���Ƃ��Ƃ����������Ȃ牺�̎����v�Z���Ă��ω����Ȃ�.
		//�܂��Ajelm�̂ق����󂩂�����(���̏ꍇ��-1�̂͂�)�A�[����ielm�Ɠ����ł���ꍇ�͂��̂܂܂ł悢
		if(jelm!=0) if(depth[jelm]>depth[ielm]) depth[jelm]=depth[ielm]+1;
	}
	////////////////////�[���C������*/

	//���׼�ݖ@���s
	//LAPLAS1(NODE,ELEM,ip,ip,kv,ibound);


    /////////�V���ɍ��ꂽ�l�ʑ̂̌����A�Â����̂�菭�Ȃ��Ȃ����ꍇ(���ʑ̂��\�������Ƃ��A�ǂ̖ʂ����ʑ̂̋��E�ł͂Ȃ������v�f�����݂���Ƃ�)
    if(*iv>ibound)
    {   
        int ir=*(iv)-ibound;		//ir�̗v�f���폜����Ȃ��Ă͂Ȃ�Ȃ�.���������ʂɍ폜�����̂ł́A�v�f�z���[��]�������邱�ƂɂȂ�
		
		for(int i=1;i<=ir;i++)
		{
			kv[i]=kv[ibound+i];
			ELEM[kv[i]].map=kv[i];	//map�Ƃ��ĂƂ肠�������̗v�f�ԍ����i�[
		}
		///��ő������map�̒l�����������ɂȂ�т�����B���̓s����A��ł�kv[1],kv[2]�E�E�E�ƒl���߂Ă���
		qsorti3D(NODE,ELEM,ir,kv);
	
		for(int i=1;i<=ir;i++)
		{   
			int ielm=kv[ir-i+1];//i��������ɂ�Air-i+1��ir����1���������Ȃ��Ă���
			ELEM[ielm].map=0;	//��̂ق���map�̏��������s���Ă��邪�Aiv>ibound�̏ꍇ�Air�̗v�f�͏���������Ă��Ȃ��̂ŁA�����ŏ�����
	    
			if(ielm!=*nelm)		//ielm�Ԗڂ̗v�f�ɁAnelm�Ԗڂ̗v�f�����㏑��
			{   
				ELEM[ielm].r[A_X]=ELEM[*nelm].r[A_X];
				ELEM[ielm].r[A_Y]=ELEM[*nelm].r[A_Y];
				ELEM[ielm].r[A_Z]=ELEM[*nelm].r[A_Z];
				ELEM[ielm].RR=ELEM[*nelm].RR;
				ELEM[ielm].material=ELEM[*nelm].material;	//���ʂ�poly�Ƃ̑���_�B�ގ��������p��
				//depth[ielm]=depth[*nelm];					//���ʂ�poly�Ƃ̑���_�B�[���������p��
				for(int j=1;j<=4;j++)
				{   
					ELEM[ielm].node[j]=ELEM[*nelm].node[j];
					int jelm=ELEM[*nelm].elm[j];
					ELEM[ielm].elm[j]=jelm;
					if(jelm!=0)
					{   
						int N=iface3D(ELEM,*nelm,jelm);
						ELEM[jelm].elm[N]=ielm;
					}
				}
			}
			*nelm=*nelm-1;
		}
    } ///////////////////*/

	return ON;//poly�����������Ƃ������邵��Ԃ�
}		


//��C�w�����֐�
int make_air_layer(vector<point3D> &NODE,vector<element3D> &ELEM,int *nelm,mpsconfig &CON,int *node_num,int KTE,double rrm,int KTJ)
{
	//��C�w���쐬����֐� return �Ƃ��Đߓ_��������Ԃ�
	cout<<"��C�w����--";
	int node=*node_num;		//���݂̐ߓ_��
	int nelm0=*nelm;		//���݂̗v�f��
	int limit_num=CON.get_add_points();		//�ǉ��ł���ő�ߓ_��
	ELEM[0].material=AIR;	//ELEM[0]�̍ގ��͒�܂��Ă��Ȃ�����A�A�N�Z�X������G���[�ɂȂ邩���B�Ȃ̂ł����œK���ɂ���Ă����B

	int *nei_num=new int[node+1];			//�ߗו\�ʗv�f��
	
	double *direct[DIMENSION];				//�@���x�N�g���쐬
    for(int D=0;D<DIMENSION;D++) direct[D]=new double [node+1];//�O�����@���x�N�g��

	for(int i=0;i<=node;i++)
	{
		nei_num[i]=0;
		for(int D=0;D<DIMENSION;D++) direct[D][i]=0;
	}

	for(int je=1;je<=nelm0;je++)
	{
//		if(ELEM[je].material==FLUID || ELEM[je].material==ELASTIC || ELEM[je].material==MAGELAST||ELEM[je].material==HYPERELAST)	//HYPERELAST�ǉ�15/2/4
		if(ELEM[je].material==FLUID || ELEM[je].material==ELASTIC || ELEM[je].material==MAGELAST)	//HYPERELAST�ǉ�15/2/4
		{
			for(int j=1;j<=4;j++)
			{
				int jelm=ELEM[je].elm[j];
				if(jelm!=0 && ELEM[jelm].material==AIR)//��C�Ɛڂ����ʂȂ�
				{
					int ia=ELEM[je].node[j%4+1];//je��jelm�̐ڂ���O�p�`���\������ߓ_�ԍ�  ia,ib,ic,ip�𒸓_�Ƃ��ĐV�����v�f�������
					int ib=ELEM[je].node[4-(j-1)/2*2];
					int ic=ELEM[je].node[3-(j/2%2)*2];

					double iaic[3];//ia��ic���޸�ِ����i�[
					double iaib[3];//ia��ib���޸�ِ����i�[
					for(int D=0;D<3;D++)
					{
						iaic[D]=NODE[ic].r[D]-NODE[ia].r[D];
						iaib[D]=NODE[ib].r[D]-NODE[ia].r[D];
					}
					///0.5(iaic�~iaib)�͎O�p�`ia,ib,ic�̖ʐς̑傫���������A�����͊O�������޸�قƂȂ�
					double S[3];//��L���޸�ِ����i�[
					S[A_X]=0.5*(iaic[A_Y]*iaib[A_Z]-iaic[A_Z]*iaib[A_Y]);
					S[A_Y]=0.5*(iaic[A_Z]*iaib[A_X]-iaic[A_X]*iaib[A_Z]);
					S[A_Z]=0.5*(iaic[A_X]*iaib[A_Y]-iaic[A_Y]*iaib[A_X]);

					double SS=sqrt(S[A_X]*S[A_X]+S[A_Y]*S[A_Y]+S[A_Z]*S[A_Z]);
					
					////�ʐ�S�����Ƃ܂���

					double n[3]={S[A_X]/SS,S[A_Y]/SS,S[A_Z]/SS};//�O�����P�ʖ@���޸��
					
					for(int D=0;D<3;D++)
					{
						direct[D][ia]+=n[D];
						direct[D][ib]+=n[D];
						direct[D][ic]+=n[D];
					}
					nei_num[ia]++;
					nei_num[ib]++;
					nei_num[ic]++;
				}
			}
		}
	}

	//�@���x�N�g���𐳋K��
	for(int i=1;i<=node;i++)
	{
		if(nei_num[i]>0)//�\�ʐߓ_�Ȃ�
		{
			double val=0;				//�@���x�N�g���̑傫��
			for(int D=0;D<DIMENSION;D++) val+=direct[D][i]*direct[D][i];
			val=sqrt(val);
			for(int D=0;D<DIMENSION;D++) direct[D][i]/=val;		//���K��
		}
	}
	//�P�ʖ@���x�N�g�������܂���
	

	//�P�ʖ@���x�N�g����̐V�ߓ_��p���čĕ������s
	
	int newnode=0;							//�V�������������ߓ_��
	int *kv=new int[KTE];					//�V�ߓ_���O�ڋ��ɂӂ��ޗv�f�Q
    int *istack=new int[KTE];				//�ꎞ�z��

	int MESH;								//locate3D�ōŏ��ɒT������v�f�ԍ�
	int add_num=CON.get_air_layer();							//��C�w�����w�������邩
	double dL=CON.get_distancebp()*rrm*CON.get_layer_depth();	//��C�w�̕�
		
	///������
	for(int i=1;i<KTE;i++) ELEM[i].map=0;
	
	///�����ߓ_�𓱓����Ă���
	    
	MESH=*nelm;

	for(int i=1;i<=node;i++)
	{
		if(nei_num[i]>0) //�\�ʐߓ_�Ȃ�
		{
			for(int k=1;k<=add_num;k++)
			{
				if(newnode<=limit_num)
				{
					
					double r[3];//�V�ߓ_�̍��W�i�[
					for(int D=0;D<3;D++) r[D]=NODE[i].r[D]+direct[D][i]*dL*k;
					double xp=r[A_X];
					double yp=r[A_Y];
					double zp=r[A_Z];
					int flag3=OFF;//ON�Ȃ�V�ߓ_�𓱓����邩�ǂ��������܂����Ƃ�������
					int loc=locate3D(NODE,ELEM,MESH,r[A_X],r[A_Y],r[A_Z]);//�V�_���܂ޗv�f��T�� �v�f�ԍ�[MESH]����T������
					if(loc>0)	//remesh�̈�̂݃f���[�j�������ɂ́A�V�_���ÓI�v�f�̈�ɐN������ꍇ���l������B���̂Ƃ���loc=0�ƂȂ��Ă��܂��B
					{
						if(ELEM[loc].material==AIR) flag3=ON;	//�V�_���܂ޗv�f����C�v�f�Ȃ�FINE���s
						MESH=loc;								//loc=0�̂Ƃ���MESH=loc����ƁA������locate3D�ł��Ȃ��̂ŁA���̕������̃J�b�R�̊O�ɏo���Ȃ��悤����
					}
					if(flag3==ON )
					{
						
						newnode++;		//�V�_������
						*node_num=*node_num+1;
						int ip=*node_num;
						for(int D=0;D<3;D++) NODE[ip].r[D]=r[D];
						NODE[ip].boundary_condition=0;
						NODE[ip].material=AIR;
						NODE[ip].remesh=ON;			//�\������Aremesh�̈�ɑ��݂��邱�ƂɂȂ�B
						NODE[ip].particleID=-1;		//�Ή����闱�q�͑��݂���
								
						//////////�O�ڋ����ɐV�ߓ_���܂ޗv�f�̒��o
						int iv=0;
						int msk=0;
				
						iv++;//�O�ڋ����ɐV�ߓ_���܂ޗv�f��
						kv[iv]=loc;
						ELEM[loc].map=1;//map��1�̗v�f�́A�O�ڋ��ɐߓ_i���܂ނƂ�������
						msk++;
						istack[msk]=loc;
						
						while(msk!=0)
						{   
							int isk=istack[msk];//���ܒ��ڂ��Ă���v�f�̔ԍ�
							msk--;
							for(int j=1;j<=4;j++)
							{
								int jelm=ELEM[isk].elm[j];//isk�Ɛڂ���v�f
								if(jelm!=0)//���ꂪ�\�ʂłȂ��Ȃ�
								{
									if(ELEM[jelm].map==0) //�܂��������ĂȂ��Ȃ�
									{   
										double rad=ELEM[jelm].RR*(1.000000+ERR);//�O�ڋ����a�̂Q��
										double dst=ELEM[jelm].r[A_X]*ELEM[jelm].r[A_X]-2.000000*ELEM[jelm].r[A_X]*xp+xp*xp;///�O�ڋ����S�ƐV�ߓ_�̋���
										if(dst<rad)///������dst>rad�Ȃ��΂ɊO�ڋ��Ɋ܂܂Ȃ�
										{
											dst+=ELEM[jelm].r[A_Y]*ELEM[jelm].r[A_Y]-2.000000*ELEM[jelm].r[A_Y]*yp+yp*yp;
											if(dst<rad)
											{
												dst+=ELEM[jelm].r[A_Z]*ELEM[jelm].r[A_Z]-2.000000*ELEM[jelm].r[A_Z]*zp+zp*zp;
												if(dst<rad)//�O�ڋ����Ɋ܂�
												{
													
												
													iv++;//�O�ڋ����ɐV�ߓ_���܂ޗv�f�����{�P
													kv[iv]=jelm;//���X�g�ɂ����
													msk++;
													istack[msk]=jelm;
													ELEM[jelm].map=1;//jelm�͊O�ڋ����ɐV�ߓ_���܂�
				
												}
											}
										}
									}
								}
							}
						}//�V�_���O�ڋ����Ɋ܂ޗv�f��iv�Ƃ��̗v�f�ԍ�kv[iv]�����Ƃ܂���

						/////////////////
						int NN=*nelm;//�����O�̗v�f��
						
						////����ꂽ���ʑ̂��l�ʑ̂ɕ������A�V�����v�f�𐶐�����B�܂��A�֐����ŐV�v�f�̍ގ������肷��
						int FLAG=poly3D_for_FINE3D(NODE,ELEM,&iv,kv,ip,nelm,CON); 
	
						if(FLAG==OFF)		//poly�����s�����ꍇ
						{
							*node_num=*node_num-1;	//�V�ߓ_������������߂�
							newnode--;
						}
						
					}
				}
			}
		}
		
	}
	
	delete [] nei_num;
	for(int D=0;D<DIMENSION;D++) delete [] direct[D];
	delete [] kv;			
    delete [] istack;
	cout<<"����  �ߓ_��+="<<newnode<<endl;
	return newnode;
}

//�ÓI�v�f�ۑ��֐�
void memorize_static_NODE_and_ELEM(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,vector<point3D> &static_NODE,vector<element3D> &static_ELEM,int node,int nelm)
{
	//NODE,ELEM�̂����A�����Ȃ��v�f�A�ߓ_������static_NODE,staticELEM�Ɋi�[����
	cout<<"�ÓI�ߓ_�E�v�f���쐬--";
	int STATIC=1;
	int DYNAMIC=2;
	int snode=0;		//�ÓI�ߓ_��
	int snelm=0;		//�ÓI�v�f��

	int *node_map=new int[node+1];					//�e�ߓ_���ÓI���ǂ����𔻒f����
	int *newID=new int[node+1];						//�e�ߓ_��static�Ƃ��ĉ��Ԗڂ̐ߓ_�ɕύX�ɂȂ�����.DYNAMIC�Ȃ�_�~�[�Ƃ���-1���i�[
	int *new_elemID=new int[nelm+1];				//�e�v�f��static�Ƃ��ĉ��Ԗڂ̗v�f�ɕύX�ɂȂ�����.DYNAMIC�Ȃ�_�~�[�Ƃ���-1���i�[

	for(int i=1;i<=nelm;i++)
	{
		new_elemID[i]=-1;
		ELEM[i].map=DYNAMIC;	//������ 
	}
	for(int i=0;i<=node;i++)
	{
		node_map[i]=DYNAMIC;	//������
		newID[i]=-1;
	}

	//�e�v�f���ÓI�����I���𔻒f����B���ʂ�.map�Ɋi�[
	for(int i=1;i<=nelm;i++)
	{
		for(int j=1;j<=4;j++)
		{
			int ip=ELEM[i].node[j];
			if(NODE[ip].remesh==OFF) ELEM[i].map=STATIC;//non-remesh�ߓ_���ЂƂł��܂�ł���΁A����͐ÓI�v�f
		}
	}////////////

	//�e�ߓ_�̐ÓI�E���I�𔻒f
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].map==STATIC)
		{
			for(int j=1;j<=4;j++) node_map[ELEM[i].node[j]]=STATIC;//�ÓI�v�f���\������ߓ_���ÓI�Ɣ��f
		}
	}

	//static_NODE�ɏo��
	point3D NODE0;
	static_NODE.push_back(NODE0);		//�ߓ_�ԍ���1����Ȃ̂ŁA�����ň��A[0]�̔z����m�ۂ������Ă���
	for(int i=1;i<=node;i++)
	{
		if(node_map[i]==STATIC)			//�ÓI�Ȑߓ_���̂�static_NODE�Ɋi�[
		{
			snode++;
			static_NODE.push_back(NODE0);
			for(int D=0;D<3;D++) static_NODE[snode].r[D]=NODE[i].r[D];
			static_NODE[snode].boundary_condition=NODE[i].boundary_condition;
			static_NODE[snode].material=NODE[i].material;
			static_NODE[snode].particleID=NODE[i].particleID;
			static_NODE[snode].remesh=NODE[i].remesh;

			newID[i]=snode;							//�ߓ_i�͐ߓ_snode�ɂȂ���
		}
	}/////////

	//////static_ELEM�쐬
	element3D ELEM0;
	static_ELEM.push_back(ELEM0);
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].map==STATIC)
		{
			snelm++;
			static_ELEM.push_back(ELEM0);
			for(int j=1;j<=4;j++)
			{
				static_ELEM[snelm].node[j]=newID[ELEM[i].node[j]];//�ߓ_�ԍ����V�����Ȃ��Ă��邱�Ƃɒ���
				static_ELEM[snelm].elm[j]=ELEM[i].elm[j];			//�v�f�ԍ��͌Â��܂܁B���ƂŏC������
			}
			static_ELEM[snelm].volume=ELEM[i].volume;
			static_ELEM[snelm].material=ELEM[i].material;
			static_ELEM[snelm].RR=ELEM[i].RR;
			static_ELEM[snelm].map=ELEM[i].map;
			for(int D=0;D<3;D++) static_ELEM[snelm].r[D]=ELEM[i].r[D];
			//�ӏ��͂ǂ�����H �ӂ�newID���Ȃ��Ƃ����Ȃ�
			new_elemID[i]=snelm;								//i�Ԗڂ̗v�f��snelm�ԖڂɂȂ���
		}
	}

	//�v�f�[�v�f���쐬
	for(int i=1;i<=snelm;i++)
	{
		for(int j=1;j<=4;j++)
		{
			int jelm=static_ELEM[i].elm[j];//�i�[����Ă���v�f�ԍ��B�Â��̂ŐV�����̂ɂ�����
			if(jelm!=0)						//�[���̂Ƃ��͏���������K�v�Ȃ�
			{
				if(new_elemID[jelm]!=-1) static_ELEM[i].elm[j]=new_elemID[jelm];
				else
				{
					static_ELEM[i].elm[j]=0;//���I�v�f�Ɛڂ��Ă���ꍇ�A�Ƃ肠�����[�����i�[
					//cout<<i<<endl;
				}
			}
		}
	}

	delete [] node_map;
	delete [] newID;
	delete [] new_elemID;
	cout<<"ok"<<endl;
}

//�ߓ_�אڗv�f���v�Z�֐�
void set_jnb3D(vector<point3D> &NODE, vector<element3D> &ELEM,int node,int nelm, int *jnb)
{
    for(int i=1;i<=node;i++) jnb[i]=0;//������
    for(int je=1;je<=nelm;je++)
    {
        for(int j=1;j<=4;j++) jnb[ELEM[je].node[j]]=jnb[ELEM[je].node[j]]+1;
    }
}

////�ߓ_�אڗv�f�ԍ��i�[�֐�
void set_nei3D(vector<point3D> &NODE, vector<element3D> &ELEM,int node,int nelm,int *jnb,int **nei)
{
    int *num=new int [node+1];///���������ϐ�
    for(int i=1;i<=node;i++)
    {
        num[i]=0;//������
        for(int j=1;j<=jnb[i];j++) nei[i][j]=0;//������
    }
    
    for(int je=1;je<=nelm;je++)
    {
        for(int j=1;j<=4;j++)
		{
			num[ELEM[je].node[j]]=num[ELEM[je].node[j]]+1;///set_jnb3D�ł���Ă邵�A��x��Ԃ��E�E�E
			nei[ELEM[je].node[j]][num[ELEM[je].node[j]]]=je;
        }
    }
    /*//check
    for(int i=1;i<=node;i++)
    {
        for(int k=1;k<=jnb[i];k++)
		{
			int jelm=nei[i][k];
			int flag=0;
			for(int j=1;j<=4;j++) if(ELEM[jelm].node[j]==i) flag=1;
			if(flag==0) cout<<"EE"<<endl;
		}
    }/////////*/
    
    delete [] num;
}  

void poly3D2(vector<point3D> &NODE,vector<element3D> &ELEM,int *iv,int *kv,int ip,int *nelm,mpsconfig &CON)
{   
    ////���̊֐��̑O�̒i�K�ŁA�V�ߓ_���O�ڋ��Ɋ܂ގl�ʑ̐�iv�Ƃ��̗v�f�ԍ�kv[iv]�����Ƃ܂��Ă���B
	////�܂��A�V�ߓ_���O�ڋ��Ɋ܂ގl�ʑ̂�ELEM[i].map=1�ƂȂ��Ă���
    
	int memory=CON.get_poly_memory();		//���I�Ɋm�ۂ��郁������
	
    int ix=0;				//�\�ʂ̐� ��ʂ�iv=ix�Ƃ͂Ȃ�Ȃ�.�ЂƂ̗v�f�������̕\�ʂ�S���̂ŁB(iv<=ix�ł���)
	int *imen[3+1];
	for(int i=1;i<=3;i++) imen[i]=new int [memory]; //���ʑ̕\�ʎO�p�`�̐ߓ_�ԍ��i�[
    int *jmen=new int [memory];		//���ʑ̕\�ʎO�p�`�ɗאڂ���l�ʑ̔ԍ��i�[
    int *kmen=new int [memory];		//���ʑ̕\�ʎO�p�`�ɗאڂ���l�ʑ̗̂אږʔԍ� (����͑扽�ʂŎ����Ɛڂ��Ă��邩)
    double *vol=new double [memory];		//���ʑ̂̑̐ς̂U�{
    
    ///////���ʑ̕\�ʐ�ix�Ƃ�����\�����钸�_�����Ƃ߂�
    int flag=0;
    while(flag==0)
    {   
        ix=0;
        flag=1;
        for(int i=1;i<=*iv;i++)//*iv�͐V�_���O�ڋ����Ɋ܂ޗv�f�̐�
        {
			if(flag==1)
			{
        		int ielm=kv[i];//�V�_���O�ڋ����Ɋ܂ޗv�f�ԍ�
				for(int j=1;j<=4;j++)
				{
					if(flag==1)
					{
						int jelm=ELEM[ielm].elm[j];//ielm�v�f�ɗאڂ���v�f�ԍ�
						int ia=ELEM[ielm].node[j%4+1];//ielm��jelm�̐ڂ���O�p�`���\������ߓ_�ԍ�  ia,ib,ic,ip�𒸓_�Ƃ��ĐV�����v�f�������
						int ib=ELEM[ielm].node[4-(j-1)/2*2];
						int ic=ELEM[ielm].node[3-(j/2%2)*2];
	                
						if(ix>49998) cout<<"pinch"<<endl;
						if(jelm==0)//�\�ʂȂ�ielm�͕\�ʗv�f�ł���Ƃ킩��
						{
							ix++;
							imen[1][ix]=ia;	//���ʑ̕\�ʎO�p�`�̒��_�ԍ�
							imen[2][ix]=ib;
							imen[3][ix]=ic;
							jmen[ix]=0;		//���ʑ̕\�ʎO�p�`�ɗאڂ���v�f�ԍ�
							kmen[ix]=0;		//���ʑ̕\�ʎO�p�`�ɗאڂ���v�f�̗אږʔԍ�  �����0������A�����ł�0����
		
							vol[ix]=volume3D(NODE,ia,ib,ic,ip);
						}
						else if(ELEM[jelm].map==0)//ielm�ɐڂ���jelm�͐V�_���O�ڋ��Ɋ܂܂Ȃ��B�܂�ielm��jelm�̋��E�͑��ʑ̕\�ʂƂ�������
						{
							ix++;
						
							imen[1][ix]=ia;
							imen[2][ix]=ib;
							imen[3][ix]=ic;
							jmen[ix]=jelm;
							kmen[ix]=iface3D(ELEM,ielm,jelm);//jelm��ielm�ɐڂ���ʔԍ�
							
							if(ix>=memory)cout<<" ix>memory"<<endl;
							vol[ix]=volume3D(NODE,ia,ib,ic,ip);
							if(vol[ix]<=ERR)//���ȏ��}3.10�̂悤�ɁA�s�K�؂Ȏl�ʑ̂�z�肵�Ă��܂��A���ʑ̐ς����ɂȂ�B���̏ꍇ�A���ȏ��ɏ����Ă���Ƃ���A�ŏ�������Ȃ���(goto 10)
							{   
								*iv=(*iv)+1;//���X�g�𑝂₷ *iv++�͂��߁H
								
								kv[*iv]=jelm;//�{���͊O�ڋ����ɐV�ߓ_���܂܂Ȃ����A���ʑ̕������\�ɂ��邽�߂Ƀ��X�g�ɂ����
								
								ELEM[jelm].map=1;
								flag=0;
							}
						}
					}
				}
			}
        }
    }
    //////���ʑ̕\�ʐ�ix�Ƃ�����\�����钸�_�����Ƃ܂���
   
    /////////////�\�ʂ̒��_�ƐV�ߓ_���Ȃ���.����Ƒ��ʑ̂���ix�̗v�f�����������

    int ibound=ix;//ix�̑���B�܂�\�ʂ̐���\��
    
    for(int i=*(iv)+1;i<=ibound;i++) ///iv�̗v�f���j�󂳂�ĐV����ibound�̗v�f����������邩��A������v�f����ibound-iv
    {   
        *nelm=*nelm+1;		//*nelm++�Ƃ����������ł͂��߁H
        kv[i]=*nelm;		//�V�ߓ_���O�ڋ��Ɋ܂ޗv�f���X�g��������B
		ELEM[*nelm].map=1;	//���������v�f�͕K���V�ߓ_���O�ڋ��Ɋ܂�(�Ƃ��Ă���).(���̍s�Ӗ��Ȃ��Ȃ��H)
    }
    
    for(int i=1;i<=ibound;i++) ELEM[kv[i]].map=0;//�}�b�s���O�̏�����
    
    for(int i=1;i<=ibound;i++)//�v�f��񐶐�
    {   
        int ielm=kv[i];
		double determ=vol[i];
		
		int ia=imen[1][i];
		int ib=imen[2][i];
		int ic=imen[3][i];
		ELEM[ielm].node[1]=ia;
		ELEM[ielm].node[2]=ib;
		ELEM[ielm].node[3]=ic;
		ELEM[ielm].node[4]=ip;		//�V�_�͂S�Ԗڂƒ�`	
		ELEM[ielm].elm[4]=jmen[i];	
		if(jmen[i]!=0) ELEM[jmen[i]].elm[kmen[i]]=ielm;
		///����ł����H
		ELEM[ielm].volume=vol[i];
		ELEM[ielm].material=AIR;
	
		sphere3D(NODE,ELEM,ia,ib,ic,ip,ielm);//�O�ڋ��̒��S�i�{���m�C�_)�Ɣ��a�̓����v�Z
    }
    ///////////////////
    
    //�v�f�\�v�f�֌W�C��/////��̏����ő�4�ʂŐڂ���v�f�ԍ��͂킩���Ă���̂ŁA�c������߂�
	//						�����ŁA1�`3�ʂ͑��ʑ̂��\������v�f�Ƃ̋��E�ʂł��邱�Ƃɒ���
    ix=0;
    for(int i=1;i<=ibound;i++)
    {
        int ielm=kv[i];
		for(int j=1;j<=3;j++)//ELEM[ielm].elm[4]�͂��łɂ��Ƃ܂�������A����ȊO�����Ƃ߂�
		{
			///ELEM[ielm].node[4]=ip�ł���
			ELEM[ielm].elm[j]=-1;				//������
			int ia=ELEM[ielm].node[j%3+1];		//j=1,2,3�̂Ƃ��A2,3,1�̏�
			int ib=ELEM[ielm].node[(j%3+1)%3+1];//j=1,2,3�̂Ƃ��A3,1,2�̏�
			int flag=0;
			for(int k=1;k<=ix;k++)
			{
				if(flag==0)
				{
					int ja=imen[1][k];
					int jb=imen[2][k];
					if(ia==ja && ib==jb)//�ߓ_����v������
					{
					    ELEM[ielm].elm[j]=jmen[k];//���炩����ؽĂ��Ă����������i�[
					    ELEM[jmen[k]].elm[kmen[k]]=ielm;
					    imen[1][k]=imen[1][ix];//k�Ԗڂ̏��͂����s�v�B�Ȃ̂Ŕz��̈�ԍŌ�̏���k�Ԗڂɂ����Ă��āA����܂ł̏��͔j������
					    imen[2][k]=imen[2][ix];
					    jmen[k]=jmen[ix];
					    kmen[k]=kmen[ix];
					    ix--;	//�҂��Ӑ�����
						flag=1;	//ELEM[ielm].elm[j]�͂��Ƃ܂����̂ŁA���̃l�X�g�ɓ���K�v�͂Ȃ��̂�flag=1
					}
				}
			}
			if(flag==0)
			{
			    ix++;			//�����ł�ix�́A[�אڊ֌W�𖞂����v�f]���܂��Ă���[��]�̐���\���B
				imen[1][ix]=ib;	//�����̐ߓ_�̕��т��L�������A�ʂ̗v�f�����̕��т𖞂����̂�҂Bib��ia�̕��т��t�ɂ��Ă��邱�Ƃɒ���
				imen[2][ix]=ia;
				jmen[ix]=ielm;
				kmen[ix]=j;
			}
		}
    }///�v�f-�v�f�֌W�C������

    /////////�V���ɍ��ꂽ�l�ʑ̂̌����A�Â����̂�菭�Ȃ��Ȃ����ꍇ(���ʑ̂��\�������Ƃ��A�ǂ̖ʂ����ʑ̂̋��E�ł͂Ȃ������v�f�����݂���Ƃ�)
    if(*iv>ibound)
    {   
	
        int ir=*(iv)-ibound;	//ir�̗v�f���폜����Ȃ��Ă͂Ȃ�Ȃ�.���������ʂɍ폜�����̂ł́A�v�f�z���[��]�������邱�ƂɂȂ�

		for(int i=1;i<=ir;i++)
		{
			kv[i]=kv[ibound+i];
			ELEM[kv[i]].map=kv[i];//map�Ƃ��ĂƂ肠�������̗v�f�ԍ����i�[
		}
		///��ő������map�̒l(�v�f�ԍ�)�����������ɂȂ�т�����B���̓s����A��ł�kv[1],kv[2]�E�E�E�ƒl���߂Ă���
		qsorti3D(NODE,ELEM,ir,kv);
	
		for(int i=1;i<=ir;i++)
		{   
			int ielm=kv[ir-i+1];//i��������ɂ�Air-i+1��ir����1���������Ȃ��Ă���
			ELEM[ielm].map=0;	//��̂ق���map�̏��������s���Ă��邪�Aiv>ibound�̏ꍇ�Air�̗v�f�͏���������Ă��Ȃ��̂ŁA�����ŏ�����
			
			if(ielm!=*nelm)		//ielm�Ԗڂ̗v�f�ɁAnelm�Ԗڂ̗v�f�����㏑��
			{   
			    ELEM[ielm].r[A_X]=ELEM[*nelm].r[A_X];
				ELEM[ielm].r[A_Y]=ELEM[*nelm].r[A_Y];
				ELEM[ielm].r[A_Z]=ELEM[*nelm].r[A_Z];
				ELEM[ielm].RR=ELEM[*nelm].RR;
				for(int j=1;j<=4;j++)
				{   
				    ELEM[ielm].node[j]=ELEM[*nelm].node[j];
				    int jelm=ELEM[*nelm].elm[j];
				    ELEM[ielm].elm[j]=jelm;
				    if(jelm!=0)
				    {   
				        int N=iface3D(ELEM,*nelm,jelm);
						ELEM[jelm].elm[N]=ielm;
						//�ގ��́H�H�H
				    }
				}
			}
			*nelm=*nelm-1;
		}
    }///////////////////*/

	delete [] jmen;
	delete [] kmen;
	delete [] vol;

	for(int i=1;i<=3;i++) delete [] imen[i];
	
}

//���b�V���f�[�^�o�͊֐�(�ގ�ver�j
void t_data_avs3(int node,int nelm,vector <point3D> &NODE,vector <element3D> &ELEM,mpsconfig &CON,int t)
{
	ofstream fout("mesh3.mgf");
	cout<<"Now ver.3 writing-----";
	int step=1;
	fout <<"# Micro AVS Geom:2.00" << endl;
	//fout << "step" << step<< endl;
	fout << "disjoint polygon" << endl;
	fout <<"element" << endl;
	fout <<"facet" << endl;
	fout <<"color" << endl;

	int count=0;//�����\�ʐ�
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material!=AIR)
		{
			for(int j=1;j<=4;j++)
			{
				int jelm=ELEM[i].elm[j];
				if(ELEM[jelm].material==AIR) count++;
			}
		}
	}///�\�ʐ������Ƃ܂���
	
	fout <<count<< endl;//�v�f���o��
	double red, green, blue;
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material!=AIR)
		{
			for(int j=1;j<=4;j++)
			{
				int jelm=ELEM[i].elm[j];
				if(ELEM[jelm].material==AIR)//���ӗ��q����C�Ȃ�
				{
					if(ELEM[i].material==FLUID || ELEM[i].material==MAGELAST)
					{
						red   = 0.0;
						green = 0.0;
						blue  = 1.0;
					}
					else if(ELEM[i].material==ELECTRODE || ELEM[i].material==MAGNET || ELEM[i].material==COIL)
					{
						red   = 1.0;
						green = 0.0;
						blue  = 0.0;
					}
					else if(ELEM[i].material==IRON)
					{
						red   = 0.0;
						green = 1.0;
						blue  = 1.0;
					}
					else
					{
						red   = 0.0;
						green = 1.0;
						blue  = 0.0;
					}
					
					fout <<"3"<< endl;
					int ia=ELEM[i].node[j%4+1];//���_���̌����Ɉʒu����O�p�`�̒��_
					int ib=ELEM[i].node[4-(j-1)/2*2];
					int ic=ELEM[i].node[3-(j/2%2)*2];
					for(int D=0;D<3;D++) fout << NODE[ia].r[D]<<" ";
					fout << red << " " << green << " " << blue;
					fout << endl;

					for(int D=0;D<3;D++) fout << NODE[ib].r[D]<<" ";
					fout << red << " " << green << " " << blue;
					fout << endl;

					for(int D=0;D<3;D++) fout << NODE[ic].r[D]<<" ";
					fout << red << " " << green << " " << blue;
					fout << endl;
				}
			}
		}
	}////*/

	fout.close();
	cout<<"OK"<<endl;
}
//check�p���b�V���f�[�^�o�͊֐�(�ގ�ver�j
void check_meshshape(int node,int *nel,vector <point3D> &NODE,vector <element3D> &ELEM,mpsconfig &CON)
{
	stringstream ss;
	ss<<"./Check"<<"/meshshape_"<<CON.get_current_step()<<".mgf";
	string filename=ss.str();
	
	ofstream fout(filename);
	if(fout.fail()){
		system("mkdir Check");
		ofstream fout(filename);
		if(fout.fail()){
			cout<<filename<<"���J���܂���B�t�H���_���m�F���Ă�������"<<endl;
			exit(1);
		}
	}
	int nelm=*nel;
	cout<<"Now ver.3 writing-----";
	int step=1;
	fout <<"# Micro AVS Geom:2.00" << endl;
	//fout << "step" << step<< endl;
	fout << "disjoint polygon" << endl;
	fout <<"element" << endl;
	fout <<"facet" << endl;
	fout <<"color" << endl;

	int count=0;//�����\�ʐ�
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==MAGNET || ELEM[i].material==COIL)//�\���������\�ʂ̍ގ�
		{
			for(int j=1;j<=4;j++)
			{
				int jelm=ELEM[i].elm[j];
				if(ELEM[jelm].material==AIR) count++;
			}
		}
	}///�\�ʐ������Ƃ܂���
	
	fout <<count<< endl;//�v�f���o��
	double red, green, blue;
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==MAGNET || ELEM[i].material==COIL)//�\���������ގ�
		{
			for(int j=1;j<=4;j++)
			{
				int jelm=ELEM[i].elm[j];
				if(ELEM[jelm].material==AIR)
				{
//					if(ELEM[i].material==FLUID || ELEM[i].material==ELASTIC || ELEM[i].material==MAGELAST || ELEM[i].material==HYPERELAST)//�I������HYPERELAST�ǉ�15/2/4
					if(ELEM[i].material==FLUID || ELEM[i].material==ELASTIC || ELEM[i].material==MAGELAST)//�I������HYPERELAST�ǉ�15/2/4
					{
						red   = 0.0;
						green = 0.0;
						blue  = 1.0;
					}
					else if(ELEM[i].material==ELECTRODE || ELEM[i].material==MAGNET)
					{
						red   = 1.0;
						green = 0.0;
						blue  = 0.0;
					}
					else
					{
						red   = 0.0;
						green = 1.0;
						blue  = 0.0;
					}
					
					fout <<"3"<< endl;
					int ia=ELEM[i].node[j%4+1];//���_���̌����Ɉʒu����O�p�`�̒��_
					int ib=ELEM[i].node[4-(j-1)/2*2];
					int ic=ELEM[i].node[3-(j/2%2)*2];
					for(int D=0;D<3;D++) fout << NODE[ia].r[D]<<" ";
					fout << red << " " << green << " " << blue;
					fout << endl;

					for(int D=0;D<3;D++) fout << NODE[ib].r[D]<<" ";
					fout << red << " " << green << " " << blue;
					fout << endl;

					for(int D=0;D<3;D++) fout << NODE[ic].r[D]<<" ";
					fout << red << " " << green << " " << blue;
					fout << endl;
				}
			}
		}
	}////*/

	fout.close();
	cout<<"OK"<<endl;
}