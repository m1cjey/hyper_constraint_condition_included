//#faset��#facet�̊ԈႢ�H�H�H
//�t�@�C���ǂݍ��݂���łł��Ă���H�H�H

#include "stdafx.h"

#define  A_R 0
#define  A_t 1 //��


using namespace std;
//.node�f�[�^�擾�֐�
void tetgen_function::GetPointList(vector<tetgen_node> &NODE, tetgenio &in, tetgenio &out)
{
	NODE.clear();

	tetgen_node temp;
	for(int i=0;i<out.numberofpoints;i++)
	{
		temp.id=i+out.firstnumber;
		for(int n=0;n<3;n++)	temp.r[n]=out.pointlist[i*3+n];
		temp.attribute=(int)out.pointattributelist[i];
		temp.boundary=out.pointmarkerlist[i];

		NODE.push_back(temp);
	}
}


//.ele�f�[�^�擾�֐�(�ȈՔ�)
void tetgen_function::GetTetrahedronList(vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	ELEM.clear();

	tetgen_element temp;
	for(int i=0;i<out.numberoftetrahedra;i++)
	{
		temp.id=i+out.firstnumber;
		for(int n=0;n<4;n++)	temp.node[n]=out.tetrahedronlist[i*4+n];
		temp.attribute=0; //���炭�A����`�ő������ƕςȐ��l�������ăo�O��̂ł����ł�0�ɂ��Ƃ�

		ELEM.push_back(temp);
	}
}
//.ele�f�[�^�擾�֐�(MAGELAST�p)
void tetgen_function::GetMTetrahedronList(vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	ELEM.clear();

	tetgen_element temp;
	for(int i=0;i<out.numberoftetrahedra;i++)
	{
		temp.id=i+out.firstnumber;
		for(int n=0;n<4;n++)	temp.node[n]=out.tetrahedronlist[i*4+n];
		temp.attribute=MAGELAST; //���炭�A����`�ő������ƕςȐ��l�������ăo�O��̂ł����ł�0�ɂ��Ƃ�

		ELEM.push_back(temp);
	}
}

void tetgen_function::Geteleattribute(vector<tetgen_node> &NODE,vector<tetgen_element> &ELEM,tetgenio &out)
{
 for(int i=0;i<out.numberoftetrahedra;i++)
    {
        int M1=(int)out.pointattributelist[out.tetrahedronlist[i*4+0]];
		int M2=(int)out.pointattributelist[out.tetrahedronlist[i*4+1]];
		int M3=(int)out.pointattributelist[out.tetrahedronlist[i*4+2]];
		int M4=(int)out.pointattributelist[out.tetrahedronlist[i*4+3]];
	
		///4���_���ׂĂ������ގ��Ȃ�v�f������ɂȂ炤�B
		///�ЂƂł��قȂ��Ă������C�ƒ�`
		if(M1==M2 && M2==M3 && M3==M4) out.tetrahedronattributelist[i]=M1;
		else if(M1!=AIR && M2!=AIR && M3!=AIR && M4!=AIR) 
		{
			if(M1!=ELASTIC && M2!=ELASTIC && M3!=ELASTIC && M4!=ELASTIC) out.tetrahedronattributelist[i]=MAGELAST;
			//if(CON.get_model_number()==15) ELEM[i].material=AIR;
			else out.tetrahedronattributelist[i]=ELASTIC; //else ELEM[i].material=FLUID;
			if(M1==IRON || M2==IRON || M3==IRON || M4==IRON) out.tetrahedronattributelist[i]=IRON;//�R�C���̗v�f�̓R�C���ړ_�̓����ɂȂ�悤�ɂ���
		}
		else out.tetrahedronattributelist[i]=AIR;
    }
}


//.ele�f�[�^�擾�֐�(�ގ��E�v�f�v�f�֌W�܂�)
void tetgen_function::GetTetrahedronList_full(vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	ELEM.clear();

	tetgen_element temp;
	for(int i=0;i<out.numberoftetrahedra;i++)
	{
		temp.id=i+out.firstnumber;										//�ߓ_�ԍ�
		for(int n=0;n<4;n++) temp.node[n]=out.tetrahedronlist[i*4+n];	//�\���ߓ_
		for(int n=0;n<4;n++) temp.nei_elem[n]=out.neighborlist[i*4+n];	//�v�f-�v�f�֌W
		temp.attribute=(int)out.tetrahedronattributelist[i];				//�ގ�
		//temp.volume=out.tetrahedronvolumelist[i];							//�̐�(�擾�ł��Ȃ�)

		ELEM.push_back(temp);
	}
}


//.face�f�[�^�擾�֐�
void tetgen_function::GetFacetList(vector<tetgen_facet> &FACE, tetgenio &in, tetgenio &out, int boundary)
{
	//��boundarymarker�͈����Ƃ��ė^���Ă���BPLC���[�h�ō�������b�V���ł͂Ȃ����ߋ��E���o�͂���Ȃ��B

	FACE.clear();

	tetgen_facet temp;
	for(int i=0;i<out.numberoftrifaces;i++)
	{
		temp.id=i+out.firstnumber;
		for(int n=0;n<3;n++)	temp.node[n]=out.trifacelist[i*3+n];
	//	temp.boundary=out.trifacemarkerlist[i];
		temp.boundary=boundary;

		FACE.push_back(temp);
	}//*/

	//for(int i=0;i<3;i++)
	//for(int n=0;n<3;n++)	cout<<out.trifacelist[i*3+n]<<endl;
}


//.node�t�@�C���쐬�֐�
void tetgen_function::MakeNodeFile(mpsconfig &CON, vector<tetgen_node> &NODE, char *filename)
{
	//cout<<filename<<" �o��"<<endl;

	ofstream fout(filename);
	
	fout<<"#node"<<endl;
	fout<<(int)NODE.size()<<" "<<"3"<<" "<<"1"<<" "<<"1"<<endl;
	for(int i=0;i<(int)NODE.size();i++)
	{
		fout<<NODE[i].id<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<" "<<NODE[i].attribute<<" "<<NODE[i].boundary<<endl;
	}

	fout.close();
}


//.node�t�@�C���쐬�֐�
void tetgen_function::MakeNodeFile_NonAttributeAndBoundary(mpsconfig &CON, vector<tetgen_node> &NODE, char *filename)
{
	cout<<filename<<" �o��"<<endl;

	ofstream fout(filename);
	
	fout<<"#node"<<endl;
	fout<<(int)NODE.size()<<" "<<"3"<<" "<<"0"<<" "<<"0"<<endl;
	for(int i=0;i<(int)NODE.size();i++)
	{
		fout<<NODE[i].id<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<endl;
	}

	fout.close();
}


//.ele�t�@�C���쐬�֐�
void tetgen_function::MakeElemFile(mpsconfig &CON, vector<tetgen_element> &ELEM, char *filename)
{
	ofstream fout(filename);

	fout<<(int)ELEM.size()<<"\t"<<"4"<<"\t"<<"1"<<endl;
	for(int i=0;i<(int)ELEM.size();i++)
	{
		fout<<ELEM[i].id<<"\t";
		for(int n=0;n<4;n++)	fout<<ELEM[i].node[n]<<"\t";
		fout<<ELEM[i].attribute<<endl;
	}

	fout.close();
}


//.face�t�@�C���쐬�֐�
void tetgen_function::MakeFaceFile(mpsconfig &CON, vector<tetgen_facet> &FACE, char *filename)
{
	ofstream fout(filename);

	fout<<(int)FACE.size()<<" "<<"1"<<endl;
	for(int i=0;i<(int)FACE.size();i++)
	{
		fout<<FACE[i].id;
		for(int n=0;n<3;n++)	fout<<" "<<FACE[i].node[n];
		fout<<" "<<FACE[i].boundary;
		fout<<endl;
	}

	fout.clear();
}


//.poly�t�@�C���쐬�֐�
void tetgen_function::MakePolyFile(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODE, vector<tetgen_facet> &FACE, char *filename)
{
	//cout<<filename<<" �o��"<<endl;

	ofstream fout(filename);

	//node list (�����ł͏o�͂��Ȃ�)
	fout<<"#node"<<endl;
	fout<<"0"<<" "<<"3"<<" "<<"1"<<" "<<"1"<<endl;
	//fout<<(int)NODE.size()<<" "<<"3"<<" "<<"0"<<" "<<"0"<<endl;
	//for(int i=0;i<(int)NODE.size();i++)
	//{
		//fout<<NODE[i].id<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<" "<<NODE[i].attribute<<" "<<NODE[i].boundary<<endl;
	//}

	//faset list
	fout<<"#faset"<<endl;
//	fout<<"#facet"<<endl;
	fout<<(int)FACE.size()<<" "<<"1"<<endl;
	for(int i=0;i<(int)FACE.size();i++)
	{
		fout<<"1"<<" "<<"0"<<" "<<FACE[i].boundary<<endl;
		//fout<<"1"<<" "<<FACE[i].boundary<<endl;
		//fout<<"3"<<" "<<FACE[i].node[A_X]<<" "<<FACE[i].node[A_Y]<<" "<<FACE[i].node[A_Z]<<" "<<FACE[i].boundary<<endl;
		fout<<"3"<<" "<<FACE[i].node[A_X]<<" "<<FACE[i].node[A_Y]<<" "<<FACE[i].node[A_Z]<<endl;
	}

	//hole list
	fout<<"#hole"<<endl;
	fout<<"0"<<endl;
	fout<<endl;

	////////////////////////////region attribute�̌��� (�z��Ɋi�[���Ă���o�͂���)/////////////////////////////
	//�ގ��̎w��́C���E���ɂ����_�̍��W�����߁C�����̍ގ����w�肷�邱�ƂŁC�������E���ɂ���v�f���S�Ă��̍ގ��ɂȂ�D
	//���͕���𔺂��C�ގ��̎w�肪����ł��邽�߁C�����ł͍s��Ȃ��D
	//��̍ގ��̏C���ɂ����āC����`�ƂȂ��Ă���v�f�𐅗v�f�Ƃ���D

	vector<region_attribute_list> REGION;
	region_attribute_list temp;
	temp.id=0;
	temp.region_number=0;
	temp.region_attribute=0;
	
	if(CON.get_model_number()==2)
	{
		//��C
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0;
		temp.r[A_Z]=CON.get_ZU()*0.9;	//��͗̈�����9���̂Ƃ���
		temp.region_number=AIR;
		temp.region_attribute=AIR;
		REGION.push_back(temp);

		//MAGNET
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0;
		temp.r[A_Z]=CON.get_magnet_Z();	//���΂̒��S�_
		temp.region_number=MAGNET;
		temp.region_attribute=MAGNET;
		REGION.push_back(temp);
	}

	//�����G���X�g�}�[
	if(CON.get_model_number()==6)
	{
		//��C
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0;
		temp.r[A_Z]=CON.get_ZU()*0.9;	//��͗̈�����9���̂Ƃ���
		temp.region_number=AIR;
		temp.region_attribute=AIR;
		REGION.push_back(temp);

		//��C(����)
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=CON.get_magnet_r();
		temp.r[A_Z]=CON.get_magnet_Z()+CON.get_magnet_H()/2+0.0001;	//���΂̏㕔
		temp.region_number=AIR;
		temp.region_attribute=AIR;
		REGION.push_back(temp);

		//�R�C��
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0;
		temp.r[A_Z]=CON.get_magnet_Z();	//���΂̒��S�_
		temp.region_number=COIL;
		temp.region_attribute=COIL;
		REGION.push_back(temp);
		int i_no;
		for(int i=0;i<NODE.size();i++){
			if(NODE[i].part_no==1)  i_no=i;//95
		}
		//MAGELAST
		temp.id+=1;
		temp.r[A_X]=NODE[i_no].r[A_X]-0.0001; //MAGELAST�̃m�[�h���ŏ��ɒǉ������
		temp.r[A_Y]=NODE[i_no].r[A_Y];
		temp.r[A_Z]=NODE[i_no].r[A_Z]+0.0001;
		temp.region_number=MAGELAST;
		temp.region_attribute=MAGELAST;
		REGION.push_back(temp);
	}
	//�����G���X�g�}�[
	if(CON.get_model_number()==7 && CON.get_model_number()==1 && CON.get_model_number()==11) 
	{
		//��C
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0;
		temp.r[A_Z]=CON.get_ZU()*0.9;	//��͗̈�����9���̂Ƃ���
		temp.region_number=AIR;
		temp.region_attribute=AIR;
		REGION.push_back(temp);

		//��C2
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0;
		temp.r[A_Z]=CON.get_magnet_Z()-CON.get_magnet_H()/2-0.001;	//���΂̒��S�_
		temp.region_number=AIR;
		temp.region_attribute=AIR;
		REGION.push_back(temp);

		//�R�C��
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0;
		temp.r[A_Z]=CON.get_magnet_Z();	//���΂̒��S�_
		temp.region_number=COIL;
		temp.region_attribute=COIL;
		REGION.push_back(temp);

		int i_no;
		for(int i=0;i<NODE.size();i++)
		{
			if(NODE[i].part_no==1) i_no=i;
		}
		//MAGELAST

		temp.id+=1;
		temp.r[A_X]=NODE[i_no].r[A_X]-0.0001; //MAGELAST�̃m�[�h���ŏ��ɒǉ������
		temp.r[A_Y]=NODE[i_no].r[A_Y];
		temp.r[A_Z]=NODE[i_no].r[A_Z]+0.0001;
		temp.region_number=MAGELAST;
		temp.region_attribute=MAGELAST;
		REGION.push_back(temp);
	}
	if(CON.get_model_number()==8)
	{
		//��C
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0;
		temp.r[A_Z]=CON.get_ZU()*0.9;	//��͗̈�����9���̂Ƃ���
		temp.region_number=AIR;
		temp.region_attribute=AIR;
		REGION.push_back(temp);

		//�R�C��
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0;
		temp.r[A_Z]=CON.get_magnet_Z();	//���΂̒��S�_
		temp.region_number=COIL;
		temp.region_attribute=COIL;
		REGION.push_back(temp);

		int i_no;
		for(int i=0;i<NODE.size();i++){
			if(NODE[i].part_no==2)  i_no=i;//
		}
/*		//MAGELAST
		temp.id+=1;
		temp.r[A_X]=NODE[i_no].r[A_X]-0.00001; //MAGELAST�̃m�[�h���ŏ��ɒǉ������
		temp.r[A_Y]=NODE[i_no].r[A_Y];
		temp.r[A_Z]=NODE[i_no].r[A_Z]+0.00001;
		temp.region_number=MAGELAST;
		temp.region_attribute=MAGELAST;
		REGION.push_back(temp);//*/
	}

	if(CON.get_model_number()==23 || CON.get_model_number()==24)
	{
		double le=CON.get_distancebp();
		//��C
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0;
		temp.r[A_Z]=CON.get_ZU()*0.9;	//��͗̈�����9���̂Ƃ���
		temp.region_number=AIR;
		temp.region_attribute=AIR;
		REGION.push_back(temp);

		//MAGNET
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0;
		temp.r[A_Z]=CON.get_magnet_Z();	//���΂̒��S�_
		temp.region_number=MAGNET;
		temp.region_attribute=MAGNET;
		REGION.push_back(temp);

		int i_no;
		for(int i=0;i<NODE.size();i++)	if(NODE[i].part_no==0)  i_no=i;

		//MAGELAST
		temp.id+=1;
		temp.r[A_X]=NODE[i_no].r[A_X]+le; //MAGELAST�̃m�[�h���ŏ��ɒǉ������
		temp.r[A_Y]=NODE[i_no].r[A_Y]+le;
		temp.r[A_Z]=NODE[i_no].r[A_Z]+le;
		temp.region_number=MAGELAST;
		temp.region_attribute=MAGELAST;
		REGION.push_back(temp);//
	}
	////////////////////////////////////////////////////////////////////*/

	//region attribute list
	fout<<"#region attribute"<<endl;
	fout<<(int)REGION.size()<<endl;
	for(int i=0;i<(int)REGION.size();i++)
	{
		fout<<REGION[i].id<<" "<<REGION[i].r[A_X]<<" "<<REGION[i].r[A_Y]<<" "<<REGION[i].r[A_Z]<<" "<<REGION[i].region_number<<" "<<REGION[i].region_attribute<<endl;
	}

	fout.close();
}

//.smesh�t�@�C���쐬�֐�
void tetgen_function::MakeSmeshFile(mpsconfig &CON, vector<tetgen_facet> &FACE, char *filename)
{
	double le=CON.get_distancebp();

	ofstream fsmesh(filename);

	//node list (�����ł͏o�͂��Ȃ�)
	fsmesh<<"#node"<<endl;
	fsmesh<<"0"<<" "<<"3"<<" "<<"1"<<" "<<"1"<<endl;
	//fsmesh<<(int)NODE.size()<<" "<<"3"<<" "<<"0"<<" "<<"0"<<endl;
	//for(int i=0;i<(int)NODE.size();i++)
	//{
		//fsmesh<<NODE[i].id<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<" "<<NODE[i].attribute<<" "<<NODE[i].boundary<<endl;
	//}

	//faset list
	fsmesh<<"#faset"<<endl;
//	fsmesh<<"#facet"<<endl;
	fsmesh<<(int)FACE.size()<<" "<<"1"<<endl;
	for(int i=0;i<(int)FACE.size();i++)
	{
		fsmesh<<"3"<<" "<<FACE[i].node[A_X]<<" "<<FACE[i].node[A_Y]<<" "<<FACE[i].node[A_Z]<<" "<<FACE[i].boundary<<endl;
	}

	//hole list
	fsmesh<<"#hole"<<endl;
	fsmesh<<"0"<<endl;
	fsmesh<<endl;

	//region attribute list
	fsmesh<<"#region attribute"<<endl;
	fsmesh<<"2"<<endl;
	fsmesh<<"1"<<" "<<"0"<<" "<<"0"<<" "<<le*2<<" "<<"1"<<" "<<"1"<<endl;
	fsmesh<<"2"<<" "<<"0"<<" "<<"0"<<" "<<-le*2<<" "<<"2"<<" "<<"2"<<endl;
	
	fsmesh.close();
}


//�e���̋��E�ʍ쐬
void tetgen_function::SetElastBoundary(mpsconfig &CON, vector<mpselastic> &PART, tetgen_config &TET, vector<tetgen_node> &NODEe, vector<tetgen_facet> &FACEe)
{
	cout<<"�����G���X�g�}�[���E�쐬"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	vector<tetgen_element> ELEMe;
	vector<int> trans;

	tetgen_node temp;
	temp.id=0;

	double le=CON.get_distancebp();	//���q�ԋ���
	int type;
	int part_no;
	
	//�����G���X�g�}�[
/*		for(int i=0;i<(int)PART.size();i++)
		{
			if(PART[i].type==MAGELAST)
			{
				part_no=i;
				for(int d=0;d<3;d++) temp.r[d]=PART[i].r[d];
				type=PART[i].type;
				temp.boundary=0;	
				if(PART[i].surface==0) temp.attribute=MAGELAST;
				else if(PART[i].surface==1) temp.attribute=FACE_P;
				trans.push_back(part_no);
				NODEe.push_back(temp);
				temp.id+=1;
			}
			
		}*/
	//�����G���X�g�}�[
	int count=0;
		for(int i=0;i<(int)PART.size();i++)
		{
			if(PART[i].type==MAGELAST)
			{			
				part_no=i;
				for(int d=0;d<3;d++) temp.r[d]=PART[i].r[d];
				type=PART[i].type;
				temp.boundary=0;	
				if(PART[i].surface==0)
				{
					temp.attribute=MAGELAST;
					count++;
				}
				else if(PART[i].surface==1) 
				{
					temp.attribute=FACE_P;
					//temp.boundary=1;
				}
				trans.push_back(part_no);
				NODEe.push_back(temp);
				temp.id+=1;
			}
		}
		cout<<"count"<<count<<endl;
	//node�t�@�C���쐬
	cout<<"MREnode�쐬-----";
	MakeNodeFile(CON, NODEe, "MAGELAST.node");
	cout<<"OK"<<endl;

	//.node�t�@�C���ǂݎ��
	in.load_node("MAGELAST");

	//�܂��͗��̐ߓ_�݂̂ŕ���
	cout<<"MRE���b�V������-----";
	tetrahedralize("", &in, &out);
	cout<<"OK"<<endl;
		// i �v�f���֐ߓ_��ǉ�(���̏������Ă��@�\���Ă���)
		// f .face�t�@�C���ɋ��E�ł͂Ȃ��ʂ��܂߂�
		// e .edge�t�@�C���̏o��(ON�ɂ���ƂȂ����~�܂��Ă��܂�)
		// n .neigh�t�@�C���̏o��

	//�o��
	out.save_nodes("MAGELAST_whole");
	out.save_elements("MAGELAST_whole");



	//////////////////�����܂łŗ��̐ߓ_�݂̂��g���āA���ׂĂ̗v�f���q�������ʂȃ��b�V�����ł���(fluid_whole�Ŋm�F�\)*/


	///////////////�s�v�ȗv�f�̍폜

	//.node�̎擾
	GetPointList(NODEe, in, out);
	//.ele�̎擾
	GetMTetrahedronList(ELEMe, in, out);

	//�����v�f�̏���
	DelThinTetrahedron(CON, TET, NODEe, ELEMe, in, out);

	//�ߓ_-�v�f�֌W
	SetRelation_NodeElem(CON, NODEe, ELEMe);
	//�v�f-�v�f�֌W
	SetRelation_ElemElem(CON, NODEe, ELEMe);
	//�v�f-�v�f�֌W���e���̕\�ʎ擾
	GetFacetList_from_neigh(CON, ELEMe, FACEe);

	//���q�ԍ���NODEe���ɑ��
	for(int i=0;i<NODEe.size();i++)	NODEe[i].part_no=trans[i];

	//�\�ʂ��\������ߓ_��I�����C�z��ԍ����l�߂� ���̓��������q�̐ߓ_���g���ꍇ�̓R�����g�A�E�g
	//SelectFaceNode(CON, NODEe, FACEe);

	/////////////////�v�f�m�F�p�t�@�C��///////////////////////////////////
	out.save_nodes("boundary_MAGELAST");	//fluid.2.node�Ɠ����t�@�C��
	MakeElemFile(CON, ELEMe, "boundary_MAGELAST.ele");
	MakeFaceFile(CON, FACEe, "boundary_MAGELAST.face");
	////////////////�����܂łŃG���X�g�}�[�̃��b�V�����؂ꂽ//////////////////////*/
//	out.save_elements("boundary_MAGELAST");
//	out.save_faces("boundary_MAGELAST");
}


//�����v�f�̏���
void tetgen_function::DelThinTetrahedron(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	cout<<"�s�v�ȗv�f�̍폜  ";

	double le=CON.get_distancebp();
	double delL;
	int flag=0;

	//del_length.dat������Γǂݍ���
	ifstream del("del_length.dat");
	if(!del)
	{
		cout<<"tetgen_config���폜����Ӓ���������  ";
		delL=TET.del_length;
	}
	else
	{
		cout<<"del_length.dat���폜����Ӓ���������  ";
		del>>delL;
		flag=1;		//�t�@�C������ǂݍ��񂾂�t���OON
	}
	del.close();

	if(flag==1)//�t�@�C���̐�����߂��Ă���
	{
		ofstream del2("del_length.dat");
		del2<<TET.del_length<<endl;
		del2.close();
	}

	cout<<"del_length="<<delL<<endl;


	cout<<"�����O�v�f��: "<<(int)ELEM.size()<<endl;

	int del_count=0;
	int i=0;
	while(i<(int)ELEM.size())
	{
		int flag1=UNDEFINED;
		int flag2=UNDEFINED;
		int flag3=UNDEFINED;	//1�ō폜
		int del=OFF;
		int count=0;

		/*//4�_���\�ʐߓ_�ō\������Ă���΃t���O1ON
		count=0;
		for(int n=0;n<4;n++)
		{
			if(NODE[ELEM[i].node[n]].attribute==FACE_P)	count+=1;
		}
		if(count==4)	flag1=ON;
		else			flag1=OFF;//*/

		/*//4�_�������ߓ_�ō\������Ă���΃t���O2ON
		count=0;
		for(int n=0;n<4;n++)
		{
			if(NODE[ELEM[i].node[n]].boundary==FRFLUID)	count+=1;
		}
		if(count==4)	flag2=ON;
		else			flag2=OFF;//*/

		//��ł������ӂ�����΃t���OON
		for(int n1=0;n1<4;n1++)
		{
			for(int n2=n1+1;n2<4;n2++)
			{
				double dis=Distance(NODE[ELEM[i].node[n1]+out.firstnumber], NODE[ELEM[i].node[n2]+out.firstnumber]);
				if(dis>le*delL)	flag3=ON;
			}
		}//*/

		if(flag3==ON)	del=ON;
		if(flag1==ON)	del=ON;
		//if(flag3==2)	del=ON;

		//�폜
		if(del==ON)
		{
			vector<tetgen_element>::iterator it=ELEM.begin();	//�C�e���[�^������
			it+=i;				//i���w��
			it=ELEM.erase(it);	//�폜���ăC�e���[�^��Ԃ�
			del_count++;
		}
		else i++;
	}

	//�v�f�ԍ��̐U�蒼��
	for(int i=0;i<(int)ELEM.size();i++)	ELEM[i].id=i+out.firstnumber;

	cout<<"������v�f��: "<<(int)ELEM.size()<<endl;
	cout<<del_count<<"�̗v�f���폜 ----------OK"<<endl;
}


//�s�v�ȗv�f�̏���(�O���_�~�[�ߓ_�@�p)
void tetgen_function::DelTetrahedron_OutsideDummy(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	cout<<"�s�v�ȗv�f�̍폜"<<endl;
	cout<<"�����O�v�f��: "<<(int)ELEM.size()<<endl;

	double le=CON.get_distancebp();

	int del_count=0;
	int i=0;
	while(i<(int)ELEM.size())
	{
		int flag=0;	//1�ō폜
		
		//1�ł��_�~�[(��C)�ߓ_������΃t���OON
		int count=0;
		for(int n=0;n<4;n++)
		{
			if(NODE[ELEM[i].node[n]].boundary==AIR)
			{
				flag=1;
				break;
			}
		}//*/

		/*//4�_���\�ʐߓ_�ō\������Ă���΃t���OON
		int count=0;
		for(int n=0;n<4;n++)
		{
			if(NODE[ELEM[i].node[n]].boundary==BOFLUID)	count+=1;
		}
		if(count==4)	flag=1;//*/


		if(flag==1)
		{
			vector<tetgen_element>::iterator it=ELEM.begin();	//�C�e���[�^������
			it+=i;				//i���w��
			it=ELEM.erase(it);	//�폜���ăC�e���[�^��Ԃ�
			del_count++;
		}
		else i++;
	}

	//�v�f�ԍ��̐U�蒼��
	for(int i=0;i<(int)ELEM.size();i++)	ELEM[i].id=i+out.firstnumber;

	cout<<"������v�f��: "<<(int)ELEM.size()<<endl;
	cout<<del_count<<"�̗v�f���폜 ----------OK"<<endl;
}


//�ߓ_-�v�f�֌W(tetgenio��edge���X�g����擾) �����O��tetgenio����f�[�^���擾���Ă�������
void tetgen_function::SetRelation_NodeNode(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	cout<<"�ߓ_-�ߓ_�֌W";

	//�ꉞ������
	for(int i=0;i<(int)NODE.size();i++)
	{
		NODE[i].nei_node.clear();
	}

	for(int i=0;i<out.numberofedges;i++)
	{
		int node1=out.edgelist[i*2+0];
		int node2=out.edgelist[i*2+1];
		
		NODE[node1].nei_node.push_back(node2);
		NODE[node2].nei_node.push_back(node1);
	}

	cout<<"----------OK"<<endl;
}


//�ߓ_-�v�f�֌W
void tetgen_function::SetRelation_NodeElem(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM)
{
	cout<<"�ߓ_-�v�f�֌W";

	//�ꉞ������
	for(int i=0;i<(int)NODE.size();i++)
	{
		NODE[i].nei_node.clear();
		NODE[i].nei_elem.clear();
	}

	for(int i=0;i<(int)ELEM.size();i++)
	{
		for(int n=0;n<4;n++)
		{
			NODE[ELEM[i].node[n]].nei_elem.push_back(i);	//�v�f��ǉ�
		}
	}

	/*//�o�� ����эő吔�E�ŏ����̏o��
	//int max=0;
	//int min=(int)NODE[0].nei_elem.size();

	ofstream fout("neigh_node-elem.dat");
	for(int i=0;i<(int)NODE.size();i++)
	{
		fout<<NODE[i].id;
		for(int n=0;n<(int)NODE[i].nei_elem.size();n++)
		{
			fout<<" "<<NODE[i].nei_elem[n];
		}
		fout<<endl;

		//�ő�ŏ��X�V
		//if(max<(int)NODE[i].nei_elem.size())	max=(int)NODE[i].nei_elem.size();
		//if(min>(int)NODE[i].nei_elem.size())	min=(int)NODE[i].nei_elem.size();
	}
	fout.clear();

	//cout<<"�ő吔: "<<max<<endl;
	//cout<<"�ŏ���: "<<min<<endl;
	//*/

	cout<<"----------OK"<<endl;
}


//�v�f-�v�f�֌W
void tetgen_function::SetRelation_ElemElem(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM)
{
	cout<<"�v�f-�v�f�֌W";
	
	vector<int> nei_all;
	
	//������
	for(int i=0;i<(int)ELEM.size();i++)
	{
		for(int n=0;n<4;n++)	ELEM[i].nei_elem[n]=-2;	//����`��-2�Ƃ��Ă���
	}

	for(int i=0;i<(int)ELEM.size();i++)
	{
		//4�̐ߓ_�̐ߓ_-�v�f�֌W�ɂ���v�f���i�[����
		nei_all.clear();
		
		for(int n=0;n<4;n++)
		{
			for(int j=0;j<(int)NODE[ELEM[i].node[n]].nei_elem.size();j++)
			{
				int elem=NODE[ELEM[i].node[n]].nei_elem[j];
				if(elem!=i)	nei_all.push_back(elem);	//�v�fi�͏���
			}
		}//nei_all�Ɋi�[����(�����v�f�ԍ����܂܂�Ă���\������)

		//�m�F�p
		//if(i==10000)	for(int j=0;j<(int)nei_all.size();j++)	cout<<nei_all[j]<<endl;

		//�ʂ�T��
		for(int ni=0;ni<4;ni++)
		{
			//�܂��͗v�fi�̖ʂ��w��
			int face[3];	//3�ԍ��Ŗʂ��w��
			int c=0;//�����グ�ϐ�

			for(int f=0;f<4;f++)
			{
				if(ni!=f)
				{
					face[c]=ELEM[i].node[f];
					c++;
				}
			}//face[3]��n�Ԗڂ̖ʂ��i�[

			//�ʂ�T��
			int correct_nei=-1;
			for(int j=0;j<(int)nei_all.size();j++)
			{
				int count=0;//���̃J�E���g��3�ɂȂ�Ίm��

				for(int nj=0;nj<4;nj++)
				{
					int node_j=ELEM[nei_all[j]].node[nj];
					for(int f=0;f<3;f++)
					{
						if(node_j==face[f])	count++;
					}	
				}
				if(count==3)
				{
					correct_nei=nei_all[j];
					break;
				}
			}//����������Ȃ�������correcr_nei�ɂ�-1�������Ă���

			ELEM[i].nei_elem[ni]=correct_nei;
		}
	}//*/

	/*//�o��
	ofstream fout("neigh.dat");
	for(int i=0;i<(int)ELEM.size();i++)
	{
		fout<<ELEM[i].id;
		for(int n=0;n<4;n++)
		{
			fout<<" "<<ELEM[i].nei_elem[n];		
		}
		fout<<endl;
	}
	fout.clear();
	//*/

	cout<<"----------OK"<<endl;
}


//�v�f������̗��̕\�ʒ�`
void tetgen_function::GetFacetList_from_neigh(mpsconfig &CON, vector<tetgen_element> &ELEM, vector<tetgen_facet> &FACE)
{
	cout<<"�v�f������̗��̕\�ʒ�`";

	//������
	FACE.clear();
	
	int id=0;	//id�p(�\�ʂ̐�)
	
	for(int i=0;i<(int)ELEM.size();i++)
	{
		for(int n=0;n<4;n++)
		{
			if(ELEM[i].nei_elem[n]==-1)//�Ζʂ����݂��Ȃ����\��
			{
				tetgen_facet temp;	
				int c=0;

				for(int f=0;f<4;f++)
				{
					if(f!=n)
					{
						temp.node[c]=ELEM[i].node[f];
						c++;
					}
				}

				temp.id=id;
				temp.boundary=MAGELAST;
				FACE.push_back(temp);
				id++;
			}
		}
	}

	cout<<"----------OK"<<endl;
}


//�\�ʐߓ_�ȊO���폜�C�\�ʐߓ_�ɔԍ���U��Ȃ���
void tetgen_function::SelectFaceNode(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_facet> &FACE)
{
	//boundary���t���O�Ɏg�킹�Ă��炤�B
	//boundary==-1�͕\�ʂ��\�����Ă��Ȃ��ߓ_���폜
	//boundary==-2�͕\�ʂ��\�����Ă���ߓ_���V�����ߓ_�ԍ���^����
	
	//�Ƃ肠����������
	for(int i=0;i<(int)NODE.size();i++)	NODE[i].boundary=-1;
	
	//FACE�ɂ���ߓ_��flag��-2��
	for(int i=0;i<(int)FACE.size();i++)
	{
		for(int n=0;n<3;n++)
		{
			NODE[FACE[i].node[n]].boundary=-2;
		}
	}

	//�V�����ߓ_�ԍ��̌���
	int id=0;
	for(int i=0;i<(int)NODE.size();i++)
	{
		if(NODE[i].boundary==-2)
		{
			NODE[i].boundary=id;
			id++;
		}
	}
	//�\�ʐߓ_�ɂ�boundary�ɐV�����ߓ_�ԍ�������B�����ߓ_�ɂ�-1������

	//FACE�̍\���ߓ_�ԍ��̕ϊ�
	for(int i=0;i<(int)FACE.size();i++)
	{
		for(int n=0;n<3;n++)
		{
			FACE[i].node[n]=NODE[FACE[i].node[n]].boundary;	//������-1��-2�ƂȂ���̂͂Ȃ��͂��B����Ώ�̏������Ԉ���Ă���
		}
	}

	//�����ߓ_�̍폜 NODE�̐ߓ_�ԍ��̕ϊ� boundary�����ɖ߂�
	int k=0;
	while(k<(int)NODE.size())
	{
		if(NODE[k].boundary==-1)
		{
			vector<tetgen_node>::iterator it=NODE.begin();	//�C�e���[�^������
			it+=k;				//k�Ԗڂ��w��
			it=NODE.erase(it);	//�폜���ăC�e���[�^��Ԃ�
		}
		else
		{
			NODE[k].id=NODE[k].boundary;
			NODE[k].boundary=WATER;
			k++;
		}
	}

}


//�_�~�[�ߓ_���폜
void tetgen_function::DelDummyNode(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_facet> &FACE, int num_dummy)
{
	//NODE�̒��ɂ͑O���ɗ��̕\�ʐߓ_�C�㔼�Ƀ_�~�[�ߓ_���ł܂��Ċi�[����Ă���̂ŁC�㔼�̃_�~�[�ߓ_�̕����݂̂������΂悢
	//�_�~�[�v�f�������Ă���\�ʃf�[�^���擾���Ă���̂ŁC�\�ʂ��\������ߓ_�ԍ��̕ύX�͕s�v

	//�_�~�[�ϐ��̐�����popback�Ŗ����̗v�f�������
	for(int i=0;i<num_dummy;i++)
	{
		NODE.pop_back();
	}
}


//��C���E�ʍ쐬
void tetgen_function::SetAirBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEa, vector<tetgen_facet> &FACEa)
{
	cout<<"��C���E�쐬"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

//	vector<tetgen_element> ELEMa;

	tetgen_node temp_n;
	temp_n.id=0;
	temp_n.attribute=AIR;
	temp_n.boundary=0;

	//����������  1�������I�t�Z�b�g���Ă���
	double dL=TET.fine_air;
	int lx = int((CON.get_XL()-0.1*dL)/dL);
	int ux = int((CON.get_XR()+0.1*dL)/dL);
	int ly = int((CON.get_YD()-0.1*dL)/dL);
	int uy = int((CON.get_YU()+0.1*dL)/dL);
	int lz = int((CON.get_ZD()-0.1*dL)/dL);
	int uz = int((CON.get_ZU()+0.1*dL)/dL);

	//�ړ_�f�[�^�쐬
	//�Ód����
	if(CON.get_model_number()==14)
	{
		for(int z=lz;z<=uz;z++)
		{
			for(int y=ly;y<=uy;y++)
			{
				for(int x=lx;x<=ux;x++)
				{
					if(x==lx || x==ux || y==ly || y==uy || z==lz || z==uz)	//��͗̈�̒[�ɗ����Ƃ��ߓ_��u��
					{
						temp_n.r[A_X]=x*dL;
						temp_n.r[A_Y]=y*dL;
						temp_n.r[A_Z]=z*dL;
						NODEa.push_back(temp_n);
						temp_n.id+=1;
					}
				}
			}
		}
	}//*/

	//���C���[����
	if(CON.get_model_number()==20)
	{
		temp_n.boundary=ELECTRODE1;

		for(int z=lz;z<=uz;z++)
		{
			for(int y=ly;y<=uy;y++)
			{
				for(int x=lx;x<=ux;x++)
				{
					if(x==lx || x==ux || y==ly || y==uy || z==lz || z==uz)	//��͗̈�̒[�ɗ����Ƃ��ߓ_��u��
					{
						temp_n.r[A_X]=x*dL;
						temp_n.r[A_Y]=y*dL;
						temp_n.r[A_Z]=z*dL;
						NODEa.push_back(temp_n);
						temp_n.id+=1;
					}
				}
			}
		}
	}

	//�����G���X�g�}�[
	//��͗̈�̒[�ɂ������ɐߓ_��u���Ă��邪for���[�v�g��Ȃ��Ă��ǂ��̂ł́H�H�H
	//���Ɖ~���`�ɒu�����ق����ǂ�
	
	{
//		temp_n.boundary=ELECTRODE1;

		//for(int z=lz;z<=uz;z++)
		//{
		//	for(int y=ly;y<=uy;y++)
		//	{
		//		for(int x=lx;x<=ux;x++)
		//		{
		//			if(x==lx || x==ux || y==ly || y==uy || z==lz || z==uz)	//��͗̈�̒[�ɗ����Ƃ��ߓ_��u��
		//			{
		//				temp_n.r[A_X]=x*dL;
		//				temp_n.r[A_Y]=y*dL;
		//				temp_n.r[A_Z]=z*dL;
		//				NODEa.push_back(temp_n);
		//				temp_n.id+=1;
		//			}
		//		}
		//	}
		//}
		/////////////////////////////////////////////////////////////////////////////
		int divN[3];
		divN[A_R]=20;
		divN[A_t]=20;
		divN[A_Z]=20;
		double regionR[2]={0.0, CON.get_RU()};
		double regionZ[2]={CON.get_ZD(), CON.get_ZU()};
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
					//���S�_
	temp_n.r[A_X]=0;
	temp_n.r[A_Y]=0;
	temp_n.r[A_Z]=Zmin;					//��͗̈�̒��
	NODEa.push_back(temp_n);
	temp_n.id+=1;

	for(int n=1;n<=divN[A_R];n++)
	{
		for(int m=0;m<divN[A_t];m++)
		{
			double r=divL[A_R]*n;
			double theta=divL[A_t]*m;
			temp_n.r[A_X]=r*cos(theta);
			temp_n.r[A_Y]=r*sin(theta);
			temp_n.r[A_Z]=Zmin;					//��͗̈�̒��
			NODEa.push_back(temp_n);
			temp_n.id+=1;				//�Ή����闱�q�����݂��Ȃ�����-1���i�[
		}
	}
	//////////////////////////////////���
					//���S�_
	temp_n.r[A_X]=0;
	temp_n.r[A_Y]=0;
	temp_n.r[A_Z]=Zmax;					//��͗̈�̒��
	NODEa.push_back(temp_n);
			temp_n.id+=1;	
	for(int n=1;n<=divN[A_R];n++)
	{
		for(int m=0;m<divN[A_t];m++)
		{
			double r=divL[A_R]*n;
			double theta=divL[A_t]*m;
			temp_n.r[A_X]=r*cos(theta);
			temp_n.r[A_Y]=r*sin(theta);
			temp_n.r[A_Z]=Zmax;					//��͗̈�̒��
			NODEa.push_back(temp_n);
			temp_n.id+=1;					//�Ή����闱�q�����݂��Ȃ�����-1���i�[
		}
	}

	//����
	double RR=divL[A_R]*divN[A_R];
	for(int n=1;n<divN[A_Z];n++)
	{
		for(int m=0;m<divN[A_t];m++)
		{
			double theta=divL[A_t]*m;
			temp_n.r[A_X]=RR*cos(theta);
			temp_n.r[A_Y]=RR*sin(theta);;
			temp_n.r[A_Z]=Zmin+n*divL[A_Z];		
			NODEa.push_back(temp_n);
			temp_n.id+=1;						//�Ή����闱�q�����݂��Ȃ�����-1���i�[
		}
	}

		/////////////////////////////////////////////////////////////////////////////
	}

	//node�t�@�C���쐬
	MakeNodeFile(CON, NODEa, "NODEa1.node");
	//.ele�̎擾
//	GetTetrahedronList(ELEMa, in, out);

	in.load_node("NODEa1");
	tetrahedralize("", &in, &out);
	out.save_nodes("boundary_air1");
	out.save_elements("boundary_air1");
	out.save_faces("boundary_air1");
	//*/

	//���E�ʃf�[�^�擾
	GetFacetList(FACEa, in, out, AIR);
	//.face�t�@�C���쐬
	//MakeFaceFile(CON, FACEa1, "FACEa1.face");
	//�����v�f�̏���
//	DelThinTetrahedron(CON, TET, NODEa, ELEMa, in, out);

	//�ߓ_-�v�f�֌W
//	SetRelation_NodeElem(CON, NODEa, ELEMa);
	//�v�f-�v�f�֌W
//	SetRelation_ElemElem(CON, NODEa, ELEMa);
	//�v�f-�v�f�֌W���e���̕\�ʎ擾
//	GetFacetList_from_neigh(CON, ELEMa, FACEa);
	

}


//���H�E�G���X�g�}�[�\�ʕt�ߒǉ��ߓ_
void tetgen_function::SetAirFineBoundary(mpsconfig &CON, vector<mpselastic> &PART, tetgen_config &TET, vector<tetgen_node> &NODEa, vector<tetgen_facet> &FACEa)
{
	//���H�̕\�ʂ̖@�������ɉ��w���̃��b�V���w���쐬����

	cout<<"�d���Ε\�ʕt�ߒǉ��ߓ_"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp;
	temp.id=0;
	temp.attribute=AIR;
	temp.boundary=0;

	//����������
	double le=CON.get_distancebp();
	double R=CON.get_magnet_r()+5*le;		//�~�����a
	double uz=CON.get_magnet_Z()+CON.get_magnet_H()/2+5*le;		//�~���ő卂��
	double lz=CON.get_magnet_Z()-CON.get_magnet_H()/2-5*le;		//�~���ŏ����� -5*le
	double dx=le;					//�~�������������b�V���e��

		//�㉺���ʂŃe�[�p�[��t���������͕�����ׂ������܂Ƃ߂Ă����������E�E�E

	//���S
	temp.r[A_X]=0;
	temp.r[A_Y]=0;	
	temp.r[A_Z]=uz;
	temp.boundary=0;
	NODEa.push_back(temp);
	temp.id+=1;

	//���
	for(double r=dx;r<R+0.1*dx;r+=dx)
	{
		int nr=static_cast<int>(2.0*PI*r/dx);//�����p�`�ŋߎ����邩�i�p�x�����̕������j
		double d_theta=360.0/static_cast<double>(nr);

		for(double theta=0.0;theta<360-0.1*d_theta;theta+=d_theta)
		{
			temp.r[A_X]=r*sin(theta*PI/180.0);
			temp.r[A_Y]=r*cos(theta*PI/180.0);	
			temp.r[A_Z]=uz;
			temp.boundary=0;
			NODEa.push_back(temp);
			temp.id+=1;
		}
	}

	//����	//R*=1.005;	//�킸���ɑ�������(���b�V�����q����̂�h������)�E�E�E����͂����
	for(double z=uz-dx;z>lz;z-=dx)
	{
			int nr=static_cast<int>(2.0*PI*R/dx);//r=dR�̎���nr=7
			double d_theta=360.0/static_cast<double>(nr);

			for(double theta=0.0;theta<360.0-0.1*d_theta;theta+=d_theta)
			{
					temp.r[A_X]=R*sin(theta*PI/180.0);
					temp.r[A_Y]=R*cos(theta*PI/180.0);
					temp.r[A_Z]=z;
					temp.boundary=0;
					NODEa.push_back(temp);
					temp.id+=1;
			}	
	}

	//����
	

	for(double r=dx;r<R+0.1*dx;r+=dx)//	for(double r=R;r>le-0.1*le;r-=le)
	{
		int nr=static_cast<int>(2.0*PI*r/dx);
		double d_theta=360.0/static_cast<double>(nr);

		for(double theta=0.0;theta<360.0-0.1*d_theta;theta+=d_theta)
		{
				temp.r[A_X]=r*sin(theta*PI/180.0);
				temp.r[A_Y]=r*cos(theta*PI/180.0);
				temp.r[A_Z]=lz;
				temp.boundary=0;
				NODEa.push_back(temp);
				temp.id+=1;
		}
	}//*/

	//.stl�t�@�C���ǂݎ��
//	in.load_stl("COIL");
	//node�t�@�C���쐬
	MakeNodeFile(CON, NODEa, "air2.node");
	//.node�t�@�C���ǂݎ��
	in.load_node("air2");

	//�܂��͗��̐ߓ_�݂̂ŕ���
	tetrahedralize("", &in, &out);
		// i �v�f���֐ߓ_��ǉ�(���̏������Ă��@�\���Ă���)
		// f .face�t�@�C���ɋ��E�ł͂Ȃ��ʂ��܂߂�
		// e .edge�t�@�C���̏o��(ON�ɂ���ƂȂ����~�܂��Ă��܂�)
		// n .neigh�t�@�C���̏o��
	//�o��
	out.save_nodes("air2_whole");
	out.save_elements("air2_whole");
	out.save_faces("air2_whole");

	GetFacetList(FACEa, in, out, COIL);
}


//���d�ɋ��E�ʍ쐬
void tetgen_function::SetPlateElectrodeBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEp, vector<tetgen_facet> &FACEp)
{
	cout<<"���d�ɋ��E�쐬"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp;
	temp.id=0;
	temp.attribute=ELECTRODE2;
	temp.boundary=ELECTRODE2;

	//����������
	double	h = TET.height_plate;						//������
	double	dh = TET.fine_plate_t;						//���ݕ������b�V���e��
	int		nh = int((TET.thickness_plate+0.1*dh)/dh);	//���ݕ���������
	double	dL = TET.fine_plate_L;						//xy�������b�V���e��
	int		nL = int((TET.length_plate/2+0.1*dL)/dL);	//xy����������(�Б�)

	//�ړ_�f�[�^�쐬
	for(int z=0;z<=nh;z++)
	{
		for(int y=-nL;y<=nL;y++)
		{
			for(int x=-nL;x<=nL;x++)
			{
				if(x==-nL || x==nL || y==-nL || y==nL || z==0 || z==nh)	//���d�ɗ̈�̒[�ɗ����Ƃ��ߓ_��u��
				{
					temp.r[A_X]=x*dL;
					temp.r[A_Y]=y*dL;
					temp.r[A_Z]=z*dh+h;
					//temp.attribute=ELECTRODE;
					//temp.boundary=ELECTRODE;
					NODEp.push_back(temp);
					temp.id+=1;
				}
			}
		}
	}//*/

	//node�t�@�C���쐬
	MakeNodeFile(CON, NODEp, "NODEp.node");

	in.load_node("NODEp");
	tetrahedralize("", &in, &out);
	out.save_nodes("boundary_plate");
	out.save_elements("boundary_plate");
	out.save_faces("boundary_plate");
	//*/

	//���E�ʃf�[�^�擾
	GetFacetList(FACEp, in, out, ELECTRODE2);
	//.face�t�@�C���쐬
	//MakeFaceFile(CON, FACEp, "FACEp.face");
}


//�~���d�ɋ��E�ʍ쐬
void tetgen_function::SetColumnElectrodeBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEc, vector<tetgen_facet> &FACEc)
{
	cout<<"�~���d�ɋ��E�쐬"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp;
	temp.id=0;
	temp.attribute=ELECTRODE1;
	temp.boundary=ELECTRODE1;

	//����������
	double le=CON.get_distancebp();
	double rc=TET.radius_column;	//�~�����a
	double L=TET.length_column;		//�~������
	double dL=TET.fine_column_L;	//�~�������������b�V���e��
	//double z0=0;					//��ʂ̈ʒu

	
	//���
	//���S
	temp.r[A_X]=0;
	temp.r[A_Y]=0;	
	temp.r[A_Z]=0;
	NODEc.push_back(temp);
	temp.id+=1;

	for(double r=le;r<rc+0.1*le;r+=le)
	{
		int nr=int(2.0*PI*r/le);
		double d_theta=360.0/(double)nr;

		for(double theta=0;theta<360-0.1*d_theta;theta+=d_theta)
		{
			//double theta=nt*d_theta;
			temp.r[A_X]=r*sin(theta*PI/180.0);
			temp.r[A_Y]=r*cos(theta*PI/180.0);	
			temp.r[A_Z]=0;
			NODEc.push_back(temp);
			temp.id+=1;
		}
	}


	//����
	int flag=0;
	double dz=le;
	double z=0;
	rc*=1.005;	//�킸���ɑ�������(���b�V�����q����̂�h������)
	int nr=int(2.0*PI*rc/(2*le));//le��2�{����Ă��邱�Ƃɒ���
	double d_theta=360.0/(double)nr;
	
	while(1)
	{
		//z�����ւ̈ړ�
		if(flag==0)			dz*=1.05;
		else if(flag==1)	dz=dL;
		z-=dz;

		if(dz>dL)	flag=1;
		if(z<-L+le)
		{
			z=-L+le;
			//flag=2;
			break;
		}

		for(double theta=0;theta<360-0.1*d_theta;theta+=d_theta)
		{
			temp.r[A_X]=rc*sin(theta*PI/180.0);
			temp.r[A_Y]=rc*cos(theta*PI/180.0);
			temp.r[A_Z]=z;
			NODEc.push_back(temp);
			temp.id+=1;
		}

		//if(flag==2)	break;
	}

	//����
	//���S
	temp.r[A_X]=0;
	temp.r[A_Y]=0;	
	temp.r[A_Z]=-L+le;
	NODEc.push_back(temp);
	temp.id+=1;

	//for(double r=le;r<rc+0.1*le;r+=le)
	for(double r=rc;r>le-0.1*le;r-=2*le)//le��2�{����Ă��邱�Ƃɒ���
	{
		int nr=int(2.0*PI*r/(2*le));//le��2�{����Ă��邱�Ƃɒ���
		double d_theta=360.0/(double)nr;

		for(double theta=0;theta<360-0.1*d_theta;theta+=d_theta)
		{
			temp.r[A_X]=r*sin(theta*PI/180.0);
			temp.r[A_Y]=r*cos(theta*PI/180.0);	
			temp.r[A_Z]=-L+le;
			NODEc.push_back(temp);
			temp.id+=1;
		}
	}//*/


	//node�t�@�C���쐬
	MakeNodeFile(CON, NODEc, "NODEc.node");

	in.load_node("NODEc");
	tetrahedralize("", &in, &out);
	out.save_nodes("boundary_column");
	out.save_elements("boundary_column");
	out.save_faces("boundary_column");
	//*/

	//���E�ʃf�[�^�擾
	GetFacetList(FACEc, in, out, ELECTRODE1);
	//.face�t�@�C���쐬
	//MakeFaceFile(CON, FACEc, "FACEc.face");
}


//�y�䋫�E�ʍ쐬
void tetgen_function::SetBaseBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEb, vector<tetgen_facet> &FACEb)
{
	cout<<"�y�䋫�E�쐬"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp;
	temp.id=0;
	temp.attribute=ELECTRODE1;
	temp.boundary=ELECTRODE1;

	//����������
	double	le = CON.get_distancebp();
	double	rc = TET.radius_column;
	double	h = TET.length_column;						//������
	double	dh = TET.fine_base;						//���ݕ������b�V���e��
	int		nh = int((TET.thickness_base+0.1*dh)/dh);	//���ݕ���������
	double	L = TET.length_base;
	double	dL = TET.fine_base;						//xy�������b�V���e��
	int		nL = int((TET.length_base/2+0.1*dL)/dL);	//xy����������(�Б�)


	//��ʐڑ�����
	temp.r[A_X]=0;
	temp.r[A_Y]=0;	
	temp.r[A_Z]=-L;
	NODEb.push_back(temp);
	temp.id+=1;

	rc*=1.005;//�~���ɍ��킹�Ă킸���ɑ�������
	for(double r=rc;r>le-0.1*le;r-=2*le)//le��2�{����Ă��邱�Ƃɒ���
	{
		int nr=int(2.0*PI*r/(2*le));//le��2�{����Ă��邱�Ƃɒ���
		double d_theta=360.0/(double)nr;

		for(double theta=0;theta<360-0.1*d_theta;theta+=d_theta)
		{
			temp.r[A_X]=r*sin(theta*PI/180.0);
			temp.r[A_Y]=r*cos(theta*PI/180.0);	
			temp.r[A_Z]=-L;
			NODEb.push_back(temp);
			temp.id+=1;
		}
	}//*/

	//��ʕ\��
	double r=rc+2*le;
	double s=2*le;
	while(r<sqrt(2.0)*L/2+dL)
	{
		int nr=int(2.0*PI*r/s);
		double d_theta=360.0/(double)nr;

		for(double theta=0;theta<360-0.1*d_theta;theta+=d_theta)
		{
			temp.r[A_X]=r*sin(theta*PI/180.0);
			temp.r[A_Y]=r*cos(theta*PI/180.0);	
			temp.r[A_Z]=-L;
			if(fabs(temp.r[A_X])<L/2 && fabs(temp.r[A_Y])<L/2)
			{
				NODEb.push_back(temp);
				temp.id+=1;
			}
		}
		r*=1.15;
		s*=1.15;
	}//*/

	//���̑��̖�
	for(int z=0;z<=nh;z++)
	{
		for(int y=-nL;y<=nL;y++)
		{
			for(int x=-nL;x<=nL;x++)
			{
				if(x==-nL || x==nL || y==-nL || y==nL || z==nh)	//���d�ɗ̈�̒[�ɗ����Ƃ��ߓ_��u��
				{
					//if(sqrt(temp.r[A_X]*temp.r[A_X]+temp.r[A_Y]*temp.r[A_Y])>rc+le || z>0)
					{
						temp.r[A_X]=x*dL;
						temp.r[A_Y]=y*dL;
						temp.r[A_Z]=-h-z*dL;
						NODEb.push_back(temp);
						temp.id+=1;
					}
				}
			}
		}
	}//*/

	//node�t�@�C���쐬
	MakeNodeFile(CON, NODEb, "NODEb.node");

	in.load_node("NODEb");
	tetrahedralize("", &in, &out);
	out.save_nodes("boundary_base");
	out.save_elements("boundary_base");
	out.save_faces("boundary_base");
	//*/

	//���E�ʃf�[�^�擾
	GetFacetList(FACEb, in, out, ELECTRODE1);
	//.face�t�@�C���쐬
	//MakeFaceFile(CON, FACEp, "FACEp.face");
}

void tetgen_function::SetMagnetBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEc, vector<tetgen_facet> &FACEc)
{
	cout<<"���΋��E�쐬"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp;
	temp.id=0;
	temp.attribute=MAGNET;
	temp.boundary=MAGNET;

	//����������
	double le=CON.get_distancebp();
	double R=CON.get_magnet_r();	//���Δ��a
	double L=CON.get_magnet_H();	//���΍���
	double Zmin=CON.get_magnet_Z()-0.5*CON.get_magnet_H();
	double Zmax=CON.get_magnet_Z()+0.5*CON.get_magnet_H();
	double dL=(Zmax-Zmin)/10;	//���Β����������b�V���e��
	double dR=le;//���a�����̕����P��

	//�㉺���ʂŃe�[�p�[��t���������͕�����ׂ������܂Ƃ߂Ă����������E�E�E

	//���S
	temp.r[A_X]=0;
	temp.r[A_Y]=0;	
	temp.r[A_Z]=Zmax;
	NODEc.push_back(temp);
	temp.id+=1;

	//���
	for(double r=dR;r<R+0.1*dR;r+=dR)
	{
		int nr=static_cast<int>(2.0*PI*r/dR);//�����p�`�ŋߎ����邩�i�p�x�����̕������j
		double d_theta=360.0/static_cast<double>(nr);

		for(double theta=0.0;theta<360-0.1*d_theta;theta+=d_theta)
		{
			temp.r[A_X]=r*sin(theta*PI/180.0);
			temp.r[A_Y]=r*cos(theta*PI/180.0);	
			temp.r[A_Z]=Zmax;
			NODEc.push_back(temp);
			temp.id+=1;
		}
	}

	//����	//R*=1.005;	//�킸���ɑ�������(���b�V�����q����̂�h������)�E�E�E����͂����
	for(double z=Zmax-dL;z>Zmin;z-=dL)
	{
		//���S�i�Y��Ȃ��j
		temp.r[A_X]=0;
		temp.r[A_Y]=0;	
		temp.r[A_Z]=z;
		NODEc.push_back(temp);
		temp.id+=1;

		for(double r=dR;r<R+0.1*dR;r+=dR)
		{
			int nr=static_cast<int>(2.0*PI*r/dR);//r=dR�̎���nr=7
			double d_theta=360.0/static_cast<double>(nr);

			for(double theta=0.0;theta<360.0-0.1*d_theta;theta+=d_theta)
			{
				temp.r[A_X]=r*sin(theta*PI/180.0);
				temp.r[A_Y]=r*cos(theta*PI/180.0);
				temp.r[A_Z]=z;
				NODEc.push_back(temp);
				temp.id+=1;
			}	
		}
	}

	//����
	//���S
	temp.r[A_X]=0;
	temp.r[A_Y]=0;	
	temp.r[A_Z]=Zmin;
	NODEc.push_back(temp);
	temp.id+=1;

	for(double r=dR;r<R+0.1*dR;r+=dR)//	for(double r=R;r>le-0.1*le;r-=le)
	{
		int nr=static_cast<int>(2.0*PI*r/dR);
		double d_theta=360.0/static_cast<double>(nr);

		for(double theta=0.0;theta<360.0-0.1*d_theta;theta+=d_theta)
		{
			temp.r[A_X]=r*sin(theta*PI/180.0);
			temp.r[A_Y]=r*cos(theta*PI/180.0);	
			temp.r[A_Z]=Zmin;
			NODEc.push_back(temp);
			temp.id+=1;
		}
	}//*/

	//node�t�@�C���쐬
	MakeNodeFile(CON, NODEc, "NODEc.node");

	in.load_node("NODEc");
	tetrahedralize("", &in, &out);
	out.save_nodes("boundary_column");
	out.save_elements("boundary_column");
	out.save_faces("boundary_column");
	//*/

	//���E�ʃf�[�^�擾
	GetFacetList(FACEc, in, out, MAGNET);
	//.face�t�@�C���쐬
	//MakeFaceFile(CON, FACEc, "FACEc.face");
}

void tetgen_function::SetCOILBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEc, vector<tetgen_facet> &FACEc)
{
	cout<<"�R�C�����E�쐬"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

//	vector<tetgen_element> ELEMc;
	tetgen_node temp;
	temp.id=0;

	//����������
	double le=CON.get_distancebp();
	double R=CON.get_magnet_r();	//���Δ��a
	double L=CON.get_magnet_H();	//���΍���
	double Zmin=CON.get_magnet_Z()-0.5*CON.get_magnet_H();
	double Zmax=CON.get_magnet_Z()+0.5*CON.get_magnet_H();
	double dL=(Zmax-Zmin)/20;	//���Β����������b�V���e��
	double dR=0.001;//���a�����̕����P��

	//�㉺���ʂŃe�[�p�[��t���������͕�����ׂ������܂Ƃ߂Ă����������E�E�E

	//���S
	temp.r[A_X]=0;
	temp.r[A_Y]=0;	
	temp.r[A_Z]=Zmax;
	temp.attribute=IRON;
	temp.boundary=0;
	NODEc.push_back(temp);
	temp.id+=1;

	//���
	for(double r=dR;r<R+0.1*dR;r+=dR)
	{
		int nr=static_cast<int>(2.0*PI*r/dR);//�����p�`�ŋߎ����邩�i�p�x�����̕������j
		double d_theta=360.0/static_cast<double>(nr);

		for(double theta=0.0;theta<360-0.1*d_theta;theta+=d_theta)
		{
			if(r<=dR*4){
			temp.r[A_X]=r*sin(theta*PI/180.0);
			temp.r[A_Y]=r*cos(theta*PI/180.0);	
			temp.r[A_Z]=Zmax;
			temp.attribute=IRON; //IRON
			temp.boundary=0;
			NODEc.push_back(temp);
			temp.id+=1;
			}
			else {
			temp.r[A_X]=r*sin(theta*PI/180.0);
			temp.r[A_Y]=r*cos(theta*PI/180.0);	
			temp.r[A_Z]=Zmax;
			temp.attribute=COIL; //COIL
			if(r>=R-0.1*dR){	//(r==R)���ƌ덷�Œʂ�Ȃ�
				int a=int((Zmax-Zmin)/dL);		//�R�C���̕\�ʂɃf�B���N���`���E������}��
				if(a%3==2) temp.boundary=23;
				else if(a%3==1) temp.boundary=22;
				else if(a%3==0) temp.boundary=21;
			}
			else temp.boundary=0;
			NODEc.push_back(temp);
			temp.id+=1;
			}
		}
	}

	//����	//R*=1.005;	//�킸���ɑ�������(���b�V�����q����̂�h������)�E�E�E����͂����
	for(double z=Zmax-dL;z>Zmin;z-=dL)
	{
		//���S�i�Y��Ȃ��j
		temp.r[A_X]=0;
		temp.r[A_Y]=0;	
		temp.r[A_Z]=z;
		temp.attribute=IRON; //IRON
		temp.boundary=0;
		NODEc.push_back(temp);
		temp.id+=1;

		for(double r=dR;r<R+0.1*dR;r+=dR)
		{
			int nr=static_cast<int>(2.0*PI*r/dR);//r=dR�̎���nr=7
			double d_theta=360.0/static_cast<double>(nr);

			for(double theta=0.0;theta<360.0-0.1*d_theta;theta+=d_theta)
			{
				if(r<=dR*4){
					temp.r[A_X]=r*sin(theta*PI/180.0);
					temp.r[A_Y]=r*cos(theta*PI/180.0);
					temp.r[A_Z]=z;
					temp.attribute=IRON; //IRON
					temp.boundary=0;
					NODEc.push_back(temp);
					temp.id+=1;
				}
				else {
					temp.r[A_X]=r*sin(theta*PI/180.0);
					temp.r[A_Y]=r*cos(theta*PI/180.0);	
					temp.r[A_Z]=z;
					temp.attribute=COIL; //COIL
					if(r>=R-0.1*dR){
						int a=int((z-Zmin)/dL);
						if(a%3==2) {temp.boundary=23;}
						else if(a%3==1) {temp.boundary=22;}
						else if(a%3==0) {temp.boundary=21;}
					}
					else {temp.boundary=0;}
					NODEc.push_back(temp);
					temp.id+=1;
				}
			}	
		}
	}

	//����
	//���S
	temp.r[A_X]=0;
	temp.r[A_Y]=0;	
	temp.r[A_Z]=Zmin;
	temp.attribute=IRON; //IRON
	temp.boundary=0;
	NODEc.push_back(temp);
	temp.id+=1;

	for(double r=dR;r<R+0.1*dR;r+=dR)//	for(double r=R;r>le-0.1*le;r-=le)
	{
		int nr=static_cast<int>(2.0*PI*r/dR);
		double d_theta=360.0/static_cast<double>(nr);

		for(double theta=0.0;theta<360.0-0.1*d_theta;theta+=d_theta)
		{
			if(r<=dR*4){
				temp.r[A_X]=r*sin(theta*PI/180.0);
				temp.r[A_Y]=r*cos(theta*PI/180.0);
				temp.r[A_Z]=Zmin;
				temp.attribute=IRON; //IRON
				temp.boundary=0;
				NODEc.push_back(temp);
				temp.id+=1;
				}
			else {
				temp.r[A_X]=r*sin(theta*PI/180.0);
				temp.r[A_Y]=r*cos(theta*PI/180.0);	
				temp.r[A_Z]=Zmin;
				temp.attribute=COIL; //COIL
				if(r>=R-0.1*dR) {temp.boundary=21;}
				else {temp.boundary=0;}
				NODEc.push_back(temp);
				temp.id+=1;
			}
		}
	}//*/

	//.stl�t�@�C���ǂݎ��
//	in.load_stl("COIL");
	//node�t�@�C���쐬
	MakeNodeFile(CON, NODEc, "COIL.node");
	//.node�t�@�C���ǂݎ��
	in.load_node("COIL");

	//�܂��͗��̐ߓ_�݂̂ŕ���
	tetrahedralize("", &in, &out);
		// i �v�f���֐ߓ_��ǉ�(���̏������Ă��@�\���Ă���)
		// f .face�t�@�C���ɋ��E�ł͂Ȃ��ʂ��܂߂�
		// e .edge�t�@�C���̏o��(ON�ɂ���ƂȂ����~�܂��Ă��܂�)
		// n .neigh�t�@�C���̏o��
	//�o��
	out.save_nodes("COIL_whole");
	out.save_elements("COIL_whole");
	out.save_faces("COIL_whole");
	//////////////////�����܂łŗ��̐ߓ_�݂̂��g���āA���ׂĂ̗v�f���q�������ʂȃ��b�V�����ł���(fluid_whole�Ŋm�F�\)*/
	//�v�f�ގ��̌���
//	Geteleattribute(NODEc,ELEMc, out);
	///////////////�s�v�ȗv�f�̍폜

	//.node�̎擾
//	GetPointList(NODEc, in, out);
	//.ele�̎擾
//	GetTetrahedronList(ELEMc, in, out);
	//�����v�f�̏���
//	DelThinTetrahedron(CON, TET, NODEc, ELEMc, in, out);

	//�ߓ_-�v�f�֌W
//	SetRelation_NodeElem(CON, NODEc, ELEMc);
	//�v�f-�v�f�֌W
//	SetRelation_ElemElem(CON, NODEc, ELEMc);
	//���E�ʃf�[�^�擾
	GetFacetList(FACEc, in, out, COIL);

/*	/////////////////�v�f�m�F�p�t�@�C��///////////////////////////////////
	out.save_nodes("boundary_COIL");	//fluid.2.node�Ɠ����t�@�C��
	MakeElemFile(CON, ELEMc, "boundary_COIL.ele");
	MakeFaceFile(CON, FACEc, "boundary_COIL.face");
	////////////////�����܂łŃG���X�g�}�[�̃��b�V�����؂ꂽ//////////////////////*/

	//���E�ʃf�[�^�擾
//	GetFacetList(FACEc, in, out, MAGNET);
}

void tetgen_function::SetIRONBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEi, vector<tetgen_facet> &FACEi)
{
	cout<<"�S�S���E�쐬"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	vector<tetgen_element> ELEMi;
	tetgen_node temp;
	temp.id=0;

	//����������
	int divN[3];
	divN[A_R]=20;//���a����������
	divN[A_t]=20;//�p�x�����������i�����p�`�ŋߎ����邩�j
	divN[A_Z]=20;//��������������
	double Rmin=0.0;
	double Rmax=CON.get_magnet_r();

	double Zmin=CON.get_magnet_Z()-0.5*CON.get_magnet_H();
	double Zmax=CON.get_magnet_Z()+0.5*CON.get_magnet_H();
	double divL[3];
	divL[A_R]=(Rmax-Rmin)/divN[A_R];
	divL[A_t]=(2*PI)/divN[A_t];
	divL[A_Z]=(Zmax-Zmin)/divN[A_Z];
	double theta2=CON.get_magnet_angle()*PI*2/360;//���΂̉�]�p�x

	//�㉺���ʂŃe�[�p�[��t���������͕�����ׂ������܂Ƃ߂Ă����������E�E�E

	for(int n=0;n<5;n++)
	{
		if(n==0)//���S
		{
			for(int k=0;k<=divN[A_Z];k++)
			{
				temp.r[A_X]=0;
				temp.r[A_Y]=0;
				temp.r[A_Z]=Zmin+divL[A_Z]*k;
				temp.attribute=IRON;				
				temp.boundary=0;			//���E����
			
				if(fabs(theta2)>0.0175)//double��0�Ɣ�r���Ă͂����Ȃ��i��v����킯���Ȃ��j
				{
					double XX=temp.r[A_X];
					double ZZ=temp.r[A_Z]-CON.get_magnet_Z();
					double newX=XX*cos(theta2)-ZZ*sin(theta2);
					double newZ=XX*sin(theta2)+ZZ*cos(theta2)+CON.get_magnet_Z();
					temp.r[A_X]=newX;
					temp.r[A_Z]=newZ;
				}
				NODEi.push_back(temp);
				temp.id+=1;	
			}
		}
		else 
		{
			for(int m=0;m<divN[A_t];m++)
			{
				for(int k=0;k<=divN[A_Z];k++)
				{
					double r=divL[A_R]*n;
					double theta=divL[A_t]*m;
					temp.r[A_X]=r*cos(theta);
					temp.r[A_Y]=r*sin(theta);
					temp.r[A_Z]=Zmin+divL[A_Z]*k;
					
					temp.attribute=IRON;	
					temp.boundary=0;			//���E����

		//			temp.attribute=IRON;
		//			 if((n==divN[A_R]) && (k%3==0))temp.boundary=23;//�O������
		//			else if((n==divN[A_R]) && (k%3==1))temp.boundary=22;
		//			else if((n==divN[A_R]) && (k%3==2))temp.boundary=21;
		/*			else if((n==6) && (k%3==0)) NODE[num].boundary_condition=21;	//�S�S�Ƃ̋��E��
					else if((n==6) && (k%3==1))NODE[num].boundary_condition=22;
					else if((n==6) && (k%3==2))NODE[num].boundary_condition=23;  */
					/*temp.boundary=0;*/			//���E����
	
					if(fabs(theta2)>0.0175)//���΂��X���Ă���ꍇ(1deg���傫���ꍇ)�E�E�Edouble��0�Ɣ�r���Ă͂����Ȃ�
					//if(theta2!=0)����̓_���I�I
					{
						double XX=temp.r[A_X];
						double ZZ=temp.r[A_Z]-CON.get_magnet_Z();
						double newX=XX*cos(theta2)-ZZ*sin(theta2);
						double newZ=XX*sin(theta2)+ZZ*cos(theta2)+CON.get_magnet_Z();
						temp.r[A_X]=newX;
						temp.r[A_Z]=newZ;
					}
					NODEi.push_back(temp);
					temp.id+=1;					//ID�i�[
				}
			}
		}
	}
	//node�t�@�C���쐬
	MakeNodeFile(CON, NODEi, "IRON.node");

	//.node�t�@�C���ǂݎ��
	in.load_node("IRON");

	//�܂��͗��̐ߓ_�݂̂ŕ���
	tetrahedralize("", &in, &out);
		// i �v�f���֐ߓ_��ǉ�(���̏������Ă��@�\���Ă���)
		// f .face�t�@�C���ɋ��E�ł͂Ȃ��ʂ��܂߂�
		// e .edge�t�@�C���̏o��(ON�ɂ���ƂȂ����~�܂��Ă��܂�)
		// n .neigh�t�@�C���̏o��

	//�o��
	out.save_nodes("IRON_whole");
	out.save_elements("IRON_whole");

	//////////////////�����܂łŗ��̐ߓ_�݂̂��g���āA���ׂĂ̗v�f���q�������ʂȃ��b�V�����ł���(fluid_whole�Ŋm�F�\)*/
	//�v�f�ގ��̌���
//	Geteleattribute(NODEc,ELEMc, out);
	///////////////�s�v�ȗv�f�̍폜

	//.node�̎擾
	GetPointList(NODEi, in, out);
	//.ele�̎擾
	GetTetrahedronList(ELEMi, in, out);

	//�����v�f�̏���
//	DelThinTetrahedron(CON, TET, NODEc, ELEMc, in, out);

	//�ߓ_-�v�f�֌W
	SetRelation_NodeElem(CON, NODEi, ELEMi);
	//�v�f-�v�f�֌W
	SetRelation_ElemElem(CON, NODEi, ELEMi);
	//���E�ʃf�[�^�擾
	GetFacetList(FACEi, in, out, IRON);

	/////////////////�v�f�m�F�p�t�@�C��///////////////////////////////////
	out.save_nodes("boundary_IRON");	//fluid.2.node�Ɠ����t�@�C��
	MakeElemFile(CON, ELEMi, "boundary_IRON.ele");
	MakeFaceFile(CON, FACEi, "boundary_IRON.face");
	////////////////�����܂łŃG���X�g�}�[�̃��b�V�����؂ꂽ//////////////////////*/

	//���E�ʃf�[�^�擾
//	GetFacetList(FACEc, in, out, MAGNET);
}



//���E�ߓ_�E���E�ʃf�[�^�̌���
void tetgen_function::UniteBoundaryData(mpsconfig &CON, 
					   vector<tetgen_node> &NODE, vector<tetgen_node> &NODEa1, vector<tetgen_node> &NODEa2, vector<tetgen_node> &NODEp, vector<tetgen_node> &NODEc, vector<tetgen_node> &NODEb, vector<tetgen_node> &NODEw, 
					   vector<tetgen_facet> &FACE, vector<tetgen_facet> &FACEa1, vector<tetgen_facet> &FACEa2, vector<tetgen_facet> &FACEp, vector<tetgen_facet> &FACEc, vector<tetgen_facet> &FACEb, vector<tetgen_facet> &FACEw, 
					   vector<int> &TRANS)
{
	cout<<"���E�ߓ_�E���E�ʃf�[�^�̌���"<<endl;

	NODE.clear();
	FACE.clear();
	tetgen_node temp_n;
	tetgen_facet temp_f;

	temp_n.boundary=0;
	temp_f.boundary=0;

	int offset_n=0;	//�ߓ_�ԍ��̃I�t�Z�b�g��
	int offset_f=0;	//�\�ʔԍ��̃I�t�Z�b�g��


	//���H���E
	for(int i=0;i<(int)NODEw.size();i++)//�ߓ_
	{
		TRANS.push_back(NODEw[i].part_no);	//TRANS�ɗ��q�ԍ����i�[(boundary�ɗ��q�ԍ����i�[���Ă���)

		temp_n=NODEw[i];
		temp_n.id+=offset_n;
		//temp_n.boundary=WATER;
		NODE.push_back(temp_n);
	}
	for(int i=0;i<(int)FACEw.size();i++)//�\��
	{
		temp_f=FACEw[i];
		for(int n=0;n<3;n++)	temp_f.node[n]+=offset_n;
		temp_f.id+=offset_f;
		temp_f.boundary=WATER;
		FACE.push_back(temp_f);
	}

	//�I�t�Z�b�g�ʍX�V
	offset_n+=(int)NODEw.size();
	offset_f+=(int)FACEw.size();

	//��C���E
	for(int i=0;i<(int)NODEa1.size();i++)//�ߓ_
	{
		temp_n=NODEa1[i];
		temp_n.id+=offset_n;
		//temp_n.boundary=AIR;
		NODE.push_back(temp_n);
	}
	for(int i=0;i<(int)FACEa1.size();i++)//�\��
	{
		temp_f=FACEa1[i];
		for(int n=0;n<3;n++)	temp_f.node[n]+=offset_n;
		temp_f.id+=offset_f;
		temp_f.boundary=AIR;
		FACE.push_back(temp_f);
	}

	//�I�t�Z�b�g�ʍX�V
	offset_n+=(int)NODEa1.size();
	offset_f+=(int)FACEa1.size();

	//��C���E ���𑜓x��
	for(int i=0;i<(int)NODEa2.size();i++)//�ߓ_
	{
		temp_n=NODEa2[i];
		temp_n.id+=offset_n;
		//temp_n.boundary=AIR;
		NODE.push_back(temp_n);
	}
	for(int i=0;i<(int)FACEa2.size();i++)//�\��
	{
		temp_f=FACEa2[i];
		for(int n=0;n<3;n++)	temp_f.node[n]+=offset_n;
		temp_f.id+=offset_f;
		temp_f.boundary=AIR;
		FACE.push_back(temp_f);
	}

	//�I�t�Z�b�g�ʍX�V
	offset_n+=(int)NODEa2.size();
	offset_f+=(int)FACEa2.size();

	//���d�ɋ��E
	for(int i=0;i<(int)NODEp.size();i++)//�ߓ_
	{
		temp_n=NODEp[i];
		temp_n.id+=offset_n;
		//temp_n.boundary=ELECTRODE2;
		NODE.push_back(temp_n);
	}
	for(int i=0;i<(int)FACEp.size();i++)//�\��
	{
		temp_f=FACEp[i];
		for(int n=0;n<3;n++)	temp_f.node[n]+=offset_n;
		temp_f.id+=offset_f;
		temp_f.boundary=ELECTRODE2;
		FACE.push_back(temp_f);
	}

	//�I�t�Z�b�g�ʍX�V
	offset_n+=(int)NODEp.size();
	offset_f+=(int)FACEp.size();

	//�~���d�ɋ��E
	for(int i=0;i<(int)NODEc.size();i++)//�ߓ_
	{
		temp_n=NODEc[i];
		temp_n.id+=offset_n;
		//temp_n.boundary=ELECTRODE1;
		NODE.push_back(temp_n);
	}
	for(int i=0;i<(int)FACEc.size();i++)//�\��
	{
		temp_f=FACEc[i];
		for(int n=0;n<3;n++)	temp_f.node[n]+=offset_n;
		temp_f.id+=offset_f;
		temp_f.boundary=ELECTRODE1;
		FACE.push_back(temp_f);
	}

	//�I�t�Z�b�g�ʍX�V
	offset_n+=(int)NODEc.size();
	offset_f+=(int)FACEc.size();

	//�y�䋫�E
	for(int i=0;i<(int)NODEb.size();i++)//�ߓ_
	{
		temp_n=NODEb[i];
		temp_n.id+=offset_n;
		//temp_n.boundary=ELECTRODE1;
		NODE.push_back(temp_n);
	}
	for(int i=0;i<(int)FACEb.size();i++)//�\��
	{
		temp_f=FACEb[i];
		for(int n=0;n<3;n++)	temp_f.node[n]+=offset_n;
		temp_f.id+=offset_f;
		temp_f.boundary=ELECTRODE1;
		FACE.push_back(temp_f);
	}

	//�I�t�Z�b�g�ʍX�V
	//offset_n+=(int)NODEb.size();
	//offset_f+=(int)FACEb.size();

	//cout<<"----------OK"<<endl;
}

//���E�ߓ_�E���E�ʃf�[�^�̒ǉ�
void tetgen_function::AddBoundaryData(mpsconfig &CON, vector<tetgen_node> &NODEall, vector<tetgen_facet> &FACEall, vector<tetgen_node> &NODE, vector<tetgen_facet> &FACE, int attribute)
{
	//NODE,FACE�Ɋi�[����Ă���e���i�̃f�[�^���CNODEall,FACEall�Ɋi�[���Ă����D

	tetgen_node temp_n;
	tetgen_facet temp_f;

	temp_n.boundary=0;
	temp_f.boundary=0;

	int offset_n=(int)NODEall.size();	//�ߓ_�ԍ��̃I�t�Z�b�g��
	int offset_f=(int)FACEall.size();	//�\�ʔԍ��̃I�t�Z�b�g��


	//�ߓ_�̒ǉ�
	for(int i=0;i<(int)NODE.size();i++)//�ߓ_
	{
		temp_n=NODE[i];
		temp_n.id+=offset_n;
		temp_n.boundary=NODE[i].boundary;//�Ȃ����R�����g�A�E�g����Ă����E�E�E
		NODEall.push_back(temp_n);
	}

	//�ʂ̒ǉ�
	for(int i=0;i<(int)FACE.size();i++)//�\��
	{
		temp_f=FACE[i];
		for(int n=0;n<3;n++) temp_f.node[n]+=offset_n;
		temp_f.id+=offset_f;
		temp_f.boundary=FACE[i].boundary;
		FACEall.push_back(temp_f);
	}
}

//TRANS[]�̊i�[
void tetgen_function::SetTRANS(vector<tetgen_node> &NODE, vector<int> &TRANS)
{
	//TRANS[i]�ɂ́A�ߓ_�ԍ�i�ɑΉ����闱�q�ԍ����i�[����B
	//FEM3D.cpp �ł́A�ߓ_�ԍ���1����n�܂�̂ŁATRANS[0]�ɂ͐錾���-1������B�i�����ł͊��ɓ����Ă���j

	for(int i=0;i<(int)NODE.size();i++)
	{
		TRANS.push_back(NODE[i].part_no);	//TRANS�ɗ��q�ԍ����i�[
	}
}


//�ގ��̌���i�܂��������j
void tetgen_function::DecisionAttribute(vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	double M[4];	//�ߓ_�̍ގ����i�[
	
	for(int i=0;i<(int)ELEM.size();i++)
	{	
		if(ELEM[i].attribute==AIR)
		{
			for(int n=0;n<4;n++)	M[n]=NODE[ELEM[i].node[n]].attribute;

			if(M[0]==AIR || M[1]==AIR || M[2]==AIR || M[3]==AIR)	//1�ł���C�ߓ_������΋�C�v�f
			{
				ELEM[i].attribute=AIR;
				out.tetrahedronattributelist[i]=AIR;
			}
			else
			{
				ELEM[i].attribute=WATER;	//�c��͐�
				out.tetrahedronattributelist[i]=WATER;
			}
		}
	}
}


//�ގ��̏C��  ���̊֐����g���Ƃ��͒��O��tetgenio����f�[�^���擾���Ă�������
void tetgen_function::ModifyAttribute(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	//tetgenio����f�[�^���擾�ς݂Ɖ��肷��
	//�ߓ_-�ߓ_�֌W�������Ă���Ɖ��肷��

	/*double le=CON.get_distancebp();
	double rc=TET.radius_column;
	double L=TET.length_column;*/

	//�v�f�̍ގ��̏C��  ���̎��_�ł͐��v�f�ƂȂ�̈�͍ގ��ԍ���0(����`)�ƂȂ��Ă���
	for(int i=0;i<(int)ELEM.size();i++)
	{
		//����`(0)�̍ގ��𐅂ɂ���
		if(ELEM[i].attribute==0)
		{
			ELEM[i].attribute=WATER;
		}
	}

	//�ߓ_�̍ގ���v�f�̍ގ��ƍ��킹��

	//�܂��S����C�ɂ���
	for(int i=0;i<(int)ELEM.size();i++)
	{
		NODE[i].attribute=AIR;
	}
	//���ߓ_�̌���
	for(int i=0;i<(int)ELEM.size();i++)
	{
		if(ELEM[i].attribute==WATER)
		{
			for(int n=0;n<4;n++)	NODE[ELEM[i].node[n]].attribute=ELEM[i].attribute;
		}
	}
	//�d�ɐߓ_�̌���
	for(int i=0;i<(int)ELEM.size();i++)
	{
		if(ELEM[i].attribute==WATER)
		{
			for(int n=0;n<4;n++)	NODE[ELEM[i].node[n]].attribute=ELEM[i].attribute;
		}
	}

	//�ߓ_-�ߓ_�֌W���ߓ_�̏C��
	for(int i=0;i<(int)NODE.size();i++)
	{
		if(NODE[i].attribute==AIR)
		{
			int num_air=0;
			int num_ele=0;
		}
	}
}


//�ގ��̏C��  tetgenio�𒼐ڕҏW
void tetgen_function::ModifyAttribute_tetgenio(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	double le=CON.get_distancebp();
	double rc=TET.radius_column;
	double L=TET.length_column;

	//�ߓ_�����}�����ꂽ�ߓ_�̍ގ����Ƃ肠������C�Ƃ���
	for(int i=0;i<out.numberofpoints;i++)
	{
		if(out.pointattributelist[i]==0)
		{
			out.pointattributelist[i]=AIR;
		}
	}//*/

	//�v�f�̍ގ��̏C��  ���̎��_�ł͐��v�f�ƂȂ�̈�͍ގ��ԍ����f�t�H���g��0(����`)�ƂȂ��Ă���
	if(CON.get_model_number()==14)
	{
		for(int i=0;i<out.numberoftetrahedra;i++)
		{
			//����`�̍ގ��𐅂ɂ���
			if(out.tetrahedronattributelist[i]==0)
			{
				out.tetrahedronattributelist[i]=WATER;//0���������̂�WATER�ɂȂ�
			}
		}
	}

	if(CON.get_model_number()==2)
	{
			for(int i=0;i<out.numberofpoints;i++){
	if(out.pointattributelist[i]==FACE_P) out.pointattributelist[i]=MAGELAST;//�v�f��߂�
	}
		for(int i=0;i<out.numberoftetrahedra;i++)
		{
			/////////////////�v�f�����_�̗v�f/////////////////////////
			int ai=(int)out.pointattributelist[out.tetrahedronlist[4*i]];
			int bi=(int)out.pointattributelist[out.tetrahedronlist[4*i+1]];
			int ci=(int)out.pointattributelist[out.tetrahedronlist[4*i+2]];
			int di=(int)out.pointattributelist[out.tetrahedronlist[4*i+3]];
			////////////////////////////////////////////////////////////
			
			if(ai==bi && bi==ci && ci==di) out.tetrahedronattributelist[i]=di;//���ׂē����f�ނȂ�v�f�����̑f��
			
	//		if(ai==AIR || bi==AIR || ci==AIR || di==AIR) out.tetrahedronattributelist[i]=AIR;//�ЂƂł���C�ړ_���܂�ł���Ȃ��C	//����������MRE���ɋ�C�v�f���ł���
			
		}
		for(int i=0;i<out.numberoftetrahedra;i++)
		{
			if(out.tetrahedronattributelist[i]==0) out.tetrahedronattributelist[i]=MAGELAST;
		}
	}
	else
	{
		for(int i=0;i<out.numberofpoints;i++){
	if(out.pointattributelist[i]==FACE_P) out.pointattributelist[i]=MAGELAST;//�v�f��߂�
	}
		for(int i=0;i<out.numberoftetrahedra;i++)
		{
			/////////////////�v�f�����_�̗v�f/////////////////////////
			int ai=(int)out.pointattributelist[out.tetrahedronlist[4*i]];
			int bi=(int)out.pointattributelist[out.tetrahedronlist[4*i+1]];
			int ci=(int)out.pointattributelist[out.tetrahedronlist[4*i+2]];
			int di=(int)out.pointattributelist[out.tetrahedronlist[4*i+3]];
			////////////////////////////////////////////////////////////
			
			if(ai==bi && bi==ci && ci==di) out.tetrahedronattributelist[i]=di;//���ׂē����f�ނȂ�v�f�����̑f��
			
			else if((ai==IRON || bi==IRON || ci==IRON || di==IRON) && !(ai==AIR || bi==AIR || ci==AIR || di==AIR)) {
				out.tetrahedronattributelist[i]=IRON;//�ЂƂł��S�ړ_���܂�ł���Ȃ�,����C���ЂƂ��܂�łȂ��Ȃ�S
			}
			//if(ai==AIR || bi==AIR || ci==AIR || di==AIR) out.tetrahedronattributelist[i]=AIR;//�ЂƂł���C�ړ_���܂�ł���Ȃ��C	//����������MRE���ɋ�C�v�f���ł���
			
		}
		for(int i=0;i<out.numberoftetrahedra;i++)
		{
			if(out.tetrahedronattributelist[i]==0) out.tetrahedronattributelist[i]=AIR;
		}
/*		//�F�̒����@���������ɐԁE�΁E�E���F ���ꂷ��Ɖ�͂����Ȃ��@���b�V���ޗ��ɂ��̒l���g�p����邽��
		for(int i=0;i<out.numberoftetrahedra;i++)
		{
			if(out.tetrahedronattributelist[i]==COIL) out.tetrahedronattributelist[i]=1;
			else if(out.tetrahedronattributelist[i]==MAGELAST) out.tetrahedronattributelist[i]=2;
			else if(out.tetrahedronattributelist[i]==AIR) out.tetrahedronattributelist[i]=3;
			else if(out.tetrahedronattributelist[i]==IRON) out.tetrahedronattributelist[i]=4;	
		}*/
	}

	//�����ǉ����ꂽ�ߓ_�̂�attribute��boundary_marker�̏C��
	//���̎��_�ł͗��̓�����d�ɓ�����d�ɓ����̐ߓ_�͋�C�ߓ_�ɂȂ��Ă���̂ł��ꂼ��̍ގ��ɏC������
	//���d�ɂ̈�ԏ�̐ߓ_�͓d�ɂ̐ߓ_�Ƃ��邽�߁A��ɐ��̐ߓ_���߂Ă���d�ɐߓ_�����߂�
	
	if(CON.get_model_number()==2)
	{
		for(int i=0;i<out.numberoftetrahedra;i++)
		{
			if(out.tetrahedronattributelist[i]==MAGNET)
			{
				for(int n=0;n<4;n++)
				{
					out.pointattributelist[out.tetrahedronlist[i*4+n]]=MAGNET;
					out.pointmarkerlist[out.tetrahedronlist[i*4+n]]=MAGNET;
				}
			//	//�ގ���ELECTRODE�ɖ߂��Ă���
			//	out.tetrahedronattributelist[i]=ELECTRODE;
			}
		}
		
	}
}


//�v�f�̍ו���
void tetgen_function::FineElement(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	tetgenio add;
	add.initialize();

	vector<tetgen_node> NODEadd;
	tetgen_node temp;
	temp.id=0;
	temp.attribute=0;
	temp.boundary=0;

	for(int i=0;i<out.numberoftetrahedra;i++)
	{
		if(out.tetrahedronattributelist[i]==AIR)
		{
			double r[3]={0,0,0};
			//for(int d=0;d<3;d++)	r[d]=out.vpointlist[i*3+d];
			for(int d=0;d<3;d++)
			{
				for(int n=0;n<4;n++)	r[d]+=out.pointlist[3*(out.tetrahedronlist[4*i+n])+d];
				r[d]/=4;
			}

			if(fabs(r[A_X])<0.0005 && fabs(r[A_Y])<0.0005 && fabs(r[A_Z])<0.0005)
			{
				temp.r[A_X]=r[A_X];
				temp.r[A_Y]=r[A_Y];	
				temp.r[A_Z]=r[A_Z];
				NODEadd.push_back(temp);
				temp.id+=1;
			}
		}
	}

	MakeNodeFile(CON, NODEadd, "output-a.node");
}//*/

//�ߓ_�ԋ����v�Z�֐�
double tetgen_function::Distance(tetgen_node &point1, tetgen_node &point2)
{
	double dis=0;

	for(int d=0;d<3;d++)	dis+=(point2.r[d]-point1.r[d])*(point2.r[d]-point1.r[d]);

	return sqrt(dis);
}


//�v�f�d�S���W�v�Z�֐�
void tetgen_function::CalcBarycentricElement(mpsconfig&, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM)
{
	double r[3]={0,0,0};

	for(int i=0;i<(int)ELEM.size();i++)
	{
		for(int d=0;d<3;d++)	for(int n=0;n<4;n++)	r[d]+=NODE[ELEM[i].node[n]].r[d];
		for(int d=0;d<3;d++)	ELEM[i].g[d]=r[d]/4;
	}
}