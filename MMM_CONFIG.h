////�V�������q��ǉ������Ƃ�������Ȃ��Ƃ����Ȃ��͇̂@u_laplacian_f���ݽد�߇AP_gradient�BP_gradient2�Ccalc_Temperature��k[i]
////�D freeon
#ifndef mmmconfig
#define mmmconfig
class MMM_config
{       
	///////��͏���

	
	
	int Hf_type;	//�O������@0:��l���� 1:����
	double Hf_H;			//�O������H�̋���
	int isWLSM;             //�d�݂��ŏ����@=0, �Z�p����=1
	int eForce;             //�d���͌v�Z�@(���i��=0�Ckelvin��=1)
	int force_t;			//�d���̓^�C�v(����=0, �̐ϗ�=1)

	int P_interval;         //���̓R���^�[�o�͊Ԋu

	
	public:
	MMM_config();
	int		get_Hf_type()	{return Hf_type;}
	int		get_Hf_H()		{return Hf_H;}
	int		get_isWLSM()	{return isWLSM;}
	int		get_eForce()	{return eForce;}
	int		get_force_t()	{return force_t;}

	int		get_P_interval(){return P_interval;}

};


#endif
