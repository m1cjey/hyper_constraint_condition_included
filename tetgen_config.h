#ifndef CLASSH
#define CLASSH


//TetGen�̐ݒ�
class tetgen_config
{
public:
	tetgen_config();
	tetgen_config(mpsconfig &CON);

	//���ΐ��@
	double magnet_height;
	double magnet_radius;

	//�Ód�������u�̐��@
	double radius_column;		//�d�ɔ��a
	double length_column;		//�d�ɒ���
	double height_plate;		//������
	double length_plate;		//���̈�ӂ̒���
	double thickness_plate;		//���̌���
	double length_base;			//�y��̈�ӂ̒���
	double thickness_base;		//�y��̌���

	//���b�V���̑e��(���E��1�v�f�̕ӂ̒������ǂꂮ�炢�ɂ��邩)
	double fine_air;
	double fine_plate_t;
	double fine_plate_L;
	double fine_column_L;
	double fine_base;

	//���H�\�ʂ̃��b�V���w�̐ݒ�
	int num_layer_out;			//���̊O�����b�V���w��
	int num_layer_in;			//���̓������b�V���w��
	double thick_layer;			//���E���b�V��1�w�̌���
	
	//�����v�f���폜����臒l
	double del_length;			//le�̉��{�ȏ�̕ӂ����v�f��������
};

//�m�[�h
class tetgen_node
{
public:
	int id;
	double r[3];
	int attribute;//����
	int boundary;//���E
	vector<int> nei_node;
	vector<int> nei_elem;
	int part_no;
};

//�t�@�Z�b�g
class tetgen_facet
{
public:
	int id;
	int node[3];
	int boundary;
};

//�v�f
class tetgen_element
{
public:
	int id;
	int node[4];
	double g[3];	//�d�S���W
	int attribute;
	double volume;
	int nei_elem[4];
};

//
class region_attribute_list
{
public:
	int id;
	double r[3];
	int region_number;
	int region_attribute;
};


#endif