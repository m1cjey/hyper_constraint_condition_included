////�V�������q��ǉ������Ƃ�������Ȃ��Ƃ����Ȃ��͇̂@u_laplacian_f�̃m���X���b�v�AP_gradient�BP_gradient2�Ccalc_Temperature��k[i]
////�D freeon
#ifndef config
#define config

using namespace std;

class mpsconfig
{    
//��͏���

	double dt;//���ԍ��ݕ�
	double dt_for_FEM; //FEM�p���ԍ���
	int step;
	int current_step;
	double current_time;

	int avoid_step;
	int avoid_step2;
	int avoid_step3;
	int avoid_step4;
	int avoid_step5;
	int avoid_step6;
	int avoid_step7;

	bool flag_cut_x_movie;
	bool flag_cut_y_movie;

	int dimension;//����
	double maxX;//��͗̈�
	double minX;
	double maxY;
	double minY;
	double maxZ;
	double minZ;
	
//�e���́i���́j�̕����l
	
	double MRE_density;		//MRE���x
	double Silicone_density;	//Silicone���x
	double nensei;		//�S���W��
	double sigma;		//�\�ʒ��͌W��
	double Cp;			//�舳��M
	double k;			//�M�`����
	double latent_H;	//���M(latent_heat)
	double MP;			//�Z�_(melting point)
	double CTE;			//�M�c���W��[1/K] 
	double E_m;			//�����O��
	double v_m;			//�|�A�\����
	double E_e;
	double v_e;
	bool nonlinear_elastic;
	
//���q�z�u�p
	
	int fluidwidth;		//���̂̑�\����
	double distancebp;	//�������q�ԋ���
	double wlength;		//���E�̕ǂ̋����B���̂̉��{��
	double height;		//���̂̏d�S����
	double ground;		//����ʂ̈ʒu
	
//���q�@�p�p�����[�^
	
	double re;		//��ʓI�ȗ��q���are
	double re2;		//���v���V�A���p��re
	double re3;		//�\�ʒ��͗p��re
	double re4;		//freeon
	double re_elastic;
	double beta;	//��
	int dx;			//index�p�i�q���B�ő嗱�q���a���܂ސ����ɂ��邱��
	double times;	//�x�N�g���v���b�g�p�̔{��

	int motion_interval;
	
	//�\�ʒ��͊֌W
	
	int surface_tension; //�\�ʒ��̓X�C�b�`�@0=OFF 1=���ȏ� 2=���q�Ԉ���
	int smooth;			//�X���[�W���O�X�C�b�`�@1=ON 0=OFF
	int SFT;				//1=ON 0=OFF �\�ʒ���(surface tension)�̉��x(T)�ˑ����X�C�b�`
	int smn;				//�X���[�W���O�̉�
	int smth_sumP;    //for ver3. curv�̃X���[�W���O 0=OFF 1=ON
	int wall_direct;		//for ver2. BDWALL�̖@���x�N�g���v�Z���@�@0=���� 1=�ʏ� 2=�������� 3=��������
	int non_slip;		//for ver.5 �S������ ON�Ȃ�BDWALL��ɉ��z�I�ɉt�̂����݂���ݒ�B���̍ۂ�BDWALL�̔z�u�ɗ���
	int suf_modify;		//�ȗ��v�Z�œ���ꂽ�ȖʂɈ�v����悤�ɗ��q�ʒu���C�����邩�A���Ȃ���
	int interpolate_curv;//�ȗ����v�Z����Ȃ��������q�̋ȗ����A���ӗ��q�̒l������
	
	//���͊֌W

	int iteration_count;//���͌v�Z���s���� �ʏ��1�ɐݒ�
	int solution;		//�A���������̉���@0=CG method 1=ICCG method
	int B_of_P;			//���͌v�Z�̉��s�� 0=���ȏ� 1=���x���U 2=0+1
	double w_div;		//���U�̏d�� 
	int HL_sw;			//���v���V�A���ɑ΂������̗��U�����s�� 1=ON 0=OFF
	int	dir_for_P;		//���͌v�Z�̍ہA�\�ʂɔ�[���̃f�B���N���l�������邩�A���Ȃ���
	int initialP;		//���͌v�Z���ɏ����l�Ƃ��Č����͂������邩�ǂ����@1=ON 0=OFF
	double CGep;		//���͌v�Z���ɂ������������
	int	omp_P;			//���͌v�Z����CG�@��openMP���g�p���邩 1=ON 0=OFF
	int negativeP;		//�����@1=���� 0=�s����
	int set_P;			//���͋�������X�C�b�` 1=ON 0=OFF
	int Pgrad;			//���͌��z�v�Z�@ 0=���ȏ� 1=�\�ʂ̂ݖ@��
	int minP;			//�ŏ����̓X�C�b�` 1=ON 0=OFF
	int ave_Pdirect;	//���͌��z�p�@���x�N�g���̕��ω� 1=ON 0=OFF
	int Pgrad_order;	//���͌��zver.4�ɂ����鐸�x 1=���` 2=2��
	int artP_sw;		//�l�H���͌v�Z�t���O 0=OFF 1=artP���� 2=Monaghan
	double artP;		//�l�H���͒l[Pa]
	double Pgrad_times;	//���͌��z���v���b�g����ۂ́A�\�ʒ��͂ɑ΂���{�� �ʏ��1�ɐݒ�
	int P_AVS;		//microAVS�p�̈��̓t�@�C�����o�͂���step�Ԋu�B0�Ȃ�OFF

	//���x��֌W
	
	int T_field;			//���x���́@�@1=ON 0=OFF
	int insulate;			//�ǂƂ̒f�M��� 0=�f�M�@1=��f�M
	int T_laplacian;		//���x�̃��v���V�A���B�@0=���ȏ��@1=��[i] 2=���U�E���z
	double wall_density;	//�ǂ̖��x[kg/m^3]
	double wall_Cp;			//�ǂ̔�M[J/(kg�EK)]
	double wall_k;			//�ǂ̔M�`����[W/(m�EK)]
	double roomT;			//���� int�^�ł悢
	int air_cool;			//��C�Ƃ̔M�`�B���l�����邩���Ȃ��� 1=ON 0=OFF
	int T_expansion;		//�M�c���v�Z�@1=ON 0=OFF
	int buoyant;			//����(���x�̃u�V�l�X�N�ߎ�)  1=ON 0=OFF
	double TplotZ;			//3D��͂ɂ����āAXY���ʂ̉��x���o�͂���Ƃ��̂y���W
	int T_AVS;				//microAVS�p�̉��x�t�@�C�����o�͂���step�Ԋu�B0�Ȃ�OFF

	//�d����̉�@

	int	EM_method;		//�d����̉�@ 0=OFF 1=FEM 2=BEM 3=���C���[�����g�@
	int	EM_calc_type;	//0=OFF 1=�d�� 2=���� 3=����(�Q�d��)
	int	EM_interval; //�d����v�Z�����X�e�b�v�Ɉ��s�����B�ʏ��1�ɐݒ�
	int region_shape;	//��͗̈�`��@0=������ 1=�~��
	double XL;			//��͗̈捶�[
	double XR;    //��͗̈�E�[
	double YU;    //��͗̈��[
	double YD;    //��͗̈扺�[
	double ZU;			//��͗̈�Z������[
	double ZD;			//��͗̈�Z�������[
	double RU;			//��͗̈悪�~���`�ƂȂ�Ƃ��̂��̔��a

	//���b�V��
	int mesh_input;		//mesh�̍쐬���@ 0:���� 1:Magnet 2:TetGen
	int remesh_sw;		//ON:remesh�̈���l�����A��������remesh OFF:FULL�̈�𖈉񕪊�
	double boxalpha;		//�X�[�p�[�{�b�N�X�̑傫�� 1�`10�̎����ɂ��邱�ƁB�ʏ��2�ɐݒ�
	int fine;			//�ߓ_�̎����ǉ����s�����ǂ��� 0:OFF 1:ON
	double co_fine;			//�ߓ_�����ǉ��ɂ�����W��.�傫���قǕ��̂��痣�ꂽ��e���Ȃ�
	int add_points;		//�����ߓ_�ǉ���
	int poly_memory;		//poly3D()�Ŋm�ۂ��郁�����̐�
	int air_layer;		//���̂̎��͂ɋ�C�w�𐶐�����w��(0�Ȃ琶�����Ȃ�)
	double layer_depth;		//��C�w�̕��B�������q�ԋ����̉��{��
	int mesh_output_interval;//���b�V�������A�L���v�f�@�̃X�e�b�v�ɑ΂��āA���X�e�b�v�Ɉ�x�o�͂��邩�B
	double FEMCGep;
	double MRTRep;

	//BEM�֌W
	int BEM_elm_type;	//�v�f�^�C�v 0:�ߓ_�v�f 1:�ӗv�f
	int BEM_M_solver;	//���E�v�f�@�ɂ�����s���@ 0:�K�E�X�̏����@ 1:BiCGStab�@
	int gauss_N;			//Gauss�ϕ��̕]���_�� 3,4,7�̂ǂꂩ
	int tree_sw;			//BEM��tree�@���g�����A�g��Ȃ���
	double tree_mesh_size;	//tree�@�Ŏg�p����ő僌�x���̃Z���̃T�C�Y�Ble�̉��{��
	int p_max;			//tree�@�̖��������̍ő區���B0����6�̐����ɑΉ�

	//FEM�֌W
	int FEM_elm_type;	//�v�f�^�C�v 0:�ߓ_�v�f 1:�ӗv�f
	int FEM_smn;			//�d���̓X���[�W���O�񐔁@0�Ȃ�OFF
	int max_DN;			//�f�B���N���^���E�������Ƃ�ő�ߓ_��
	int FEMCG;			//FEM�ɂ�����s���@ 0:CG 1:ICCG
	double CGaccl;			//CG,ICCG�@�ɂ���������t�@�N�^�@CGaccelerator
	double EMCGep;			//�d�����ICCG�̎�������
	double FEMtimes;		//�d���͂��v���b�g����ۂ́A�\�ʒ��͂ɑ΂���{�� �ʏ��1�ɐݒ�
	double legend_F;		//F.dat�̖}��ɏo�͂����[N]
	int plot_F_type;	//F.dat�̏o�̓^�C�v  0=[N]�\�� 1=[Pa]�\��
	double legend_F_Pa;	//F.dat�̖}��ɏo�͂����[Pa]

	///////FEM ver2�֌W
	double surface_depth;			//�\�ʂ̌��݂̔���(�����le�����������̂��{���̌���(�̔���))
	
	///////�d�E�v�Z
	double V;			//�d��
	int V_step;			//�d������J�n�X�e�b�v
	double r_perm;		//��U�d�� (relative permittivity)
	int V_con;			//�d�������@0:�p���X�@1:���j�A�@2:���萔
	double initial_V;	//�d���������j�A�̂Ƃ��̏����d��
	double E_times;		//�t�@�C���o�͂���ۂ́A�d�E�̔{��
	int eleforce;		//�Ód�͌v�Z���@ 1:�ߓ_�͖@ 2:�ϕ���
	int charge;			//�d�ׂ��l�����邩���Ȃ����B0=OFF 1=ON
	int plot_E_type;	//�d�E�o�̓^�C�v 1=�x�N�g�� 2=�X�J���[

	//����v�Z
	
	int J_input_way;		//�d�����x������@ 0:���� 1:�\�t�g
	double J0;				//�����d�����x[A/m2]
	double I0;				//�����d���l[A]
	double RP;				//�䓧����(Relative Permeability)
	double ele_conduc;		//�d�C�`����
	int Hz;					//�𗬂̎��g��
	int div_Hz;				//�P�����̕�����(��͐��x)
	int jheat;				//�Q�d���ɂ�锭�M���l�����邩�@0=OFF 1=ON
	int	m_force;			//�d���͌v�Z���� 0=�}�N�X�E�F���̉��� 1=�̐ϗ�
	int NLBHI;				//�̐ϗ͂ɂ����āA�v�f�a����v�f�g�����߂�ۂɔ���`�����l�����邩�A���Ȃ���(non linier B H inverter)
	int NLMH;				//�l�̎Z�o�ɔ���`�����l�����邩�A���Ȃ���
	double magnet_H;		//�i�v���΂̑傫��
	double magnet_B;		//�i�v���΂̋���[T]
	double magnet_angle;	//�i�v���΂̒������� 0�Ȃ�+Z�����ƂȂ�B��������p�x���������Ȃ�A���̊p�x[deg]����͂���
	double magnet_r;		//�i�v���΂̔��a
	double magnet_Z;		//�i�v���΂̒��S��Z���W
	int uniform_B_sw;		//��͗̈撆�Ɉ�l����𔭐������邩�ۂ� 0=OFF 1=ON
	int magnetic_layers;	//�i�v���Ύ��ӂ̋�C�w�̐�
	double uniform_B;		//��l����̑傫��[T]
	double B_times;			//�t�@�C���o�͂���ۂ́A�������x�̔{��
	int plot_B_type;		//�������x�o�̓^�C�v 1=�x�N�g�� 2=�X�J���[

	int FEM_calc_type;		//0=OFF 1=�d�� 2=���� 3=����(�Q�d��)
	int ele_type;
	//�e��X�C�b�`�@�ǂ̂悤�Ȍv�Z���l�����邩�A���Ȃ���
	
	double g;					//�d�͉����x�@�ʏ��-9.8�Ƃ��邱��
	int restart;				//restart�p�X�C�b�`
	int autosave;			//�I�[�g�Z�[�u�Ԋu�B�����ɂ������Ƃ��͑傫�Ȑ����������Ă���
	double courant;			//�N�[���������� 0�Ȃ�OFF
	int modify_position;	//���q�ԋ�����le��菬�����ꍇ������C�����邩�A���Ȃ���
	int vis_calc_type;		//�S�����v�Z��@�@0=�z��@ 1=�A��@
	int wall_adheision;		//�ǂ̔S����ԁ@1=�m���X���b�v 0=�ذ�د��
	int laplacian;			//���v���V�A���B�@0=���ȏ��@1=��[i] 2=���U�E���z
	int	vis_solver;			//�S�������A��͂ŉ����ۂ̍s����ް 0:CG 1:ICCG
	int initial_u;			//�S�������A�I�ɉ����ۂɁA�����l�Ƃ��Č��ݑ��x����͂��邩�A���Ȃ���
	int temporary_r;		//�z��͌�ɉ��̈ʒu���v�Z���邩���Ȃ����B 1=ON 0=OFF �ʏ��ON�ɐݒ�
	
	int fix_center;			//1=ON 0=OFF
	int freeon;				//���q�ˑ��֌W�֐��@1:���񉻉\ 2:����s��
	int freeon3sw;			//freeon3���v�Z���邩���Ȃ��� 1=ON 0=OFF
	int surface_judge2;		//surface_judge2()���g�p���邩���Ȃ���  1=ON 0=OFF
	int move_prtcl;			//�ړ����q���l�����邩���Ȃ��� 1=ON 0=OFF
	int move_u_dirct;		//�ړ����q���ړ���������@���݂́}X����=�}1,�}Y����=�}2,�}Z����=�}3
	double move_speed;		//�ړ����q�̈ړ����x[m/s]
	int check_something;	//check_something()�����s���邩���Ȃ��� 1=ON 0=OFF
	
	int model_number;
	int model_set_way;		//model���Z�b�g������@�@0=�����i�q 1=MD

	double R1;

	bool switch_FEM;			//FEM�����s���邩���Ȃ��� OFF: 0, ON: 1
	bool switch_vis;			//�S�e���v�Z���邩���Ȃ��� OFF: 0, ON: 1

	///////���x�v���b�g�ϐ�
	int speed_plot_particle;	//���x���v���b�g���闱�q�̎�� 1=���ׂ� 2=fluid
	double speedtimes;			//���x�v���b�g���́A���W�ɑ΂��鑬�x�̔{��
	int speed_face;				//speed.dat�̏o�͖� 0=YZ���� 1=XZ 2=XY
	double speed_face_p;		//3D��͎���speed.dat�̏o�͖ʂ̍��W
	int ax_sym_modify;			//3D����speed.dat�Ɋւ��āA���Ώ̂ɂ��o�͏C�����s�����ۂ��@1=ON 0=OF
	int flat_speed_plot;		//���������̑��x���v���b�g���邩���Ȃ���
	double flat_speed_p;		//flat_speed.dat�̏o�͖ʂ̍��W
	int relative_speed;			//�d�S�ɑ΂��鑊�Α��x���o�͂��邩���Ȃ��� 1=ON 0=OFF
	int speed_AVS;				//microAVS�ɂ��3D���x���z�o�͂��邩���Ȃ��� 1=ON 0=OFF
	double legend_speed;		//speed.dat�̖}��ɏo�͂��鑬�x[m/s]
	int set_zero_speed;			//restart���ɑ��x���[���Z�b�g���邩���Ȃ���  1=ON 0=OFF

	//���̓`�F�b�N�p
	int pressure_face;
	double pressure_face_p;
	double pressure_times;
	double legend_pressure;		//speed.dat�̖}��ɏo�͂��鑬�x[m/s]

	///////GPU�֌W
	int M_form;				//CG�@�ɂ�����W���s��̊i�[���@ CSR_scl,CSR_vec,ELL��3���߰�
	int MAX_thread;			//�ЂƂ�SM������̍ő�X���b�h���@�ӂ���512
	
	//�t�@�C���o�͕ϐ�
	
	int interval;			//AVS�p�X�e�b�v�Ԋu
	int AVS;				//0:���ʁ@1:���́@2:���x

	//////AVS�t�@�C���o��
	int avs_eforce_interval;//AVS�d���̓t�@�C���������1��o�͂��邩 0=�o�͂��Ȃ�
	int F_interval;//�d����F.dat�����X�e�b�v���ɏo�͂��邩
	int avs_mesh1_interval;	//AVS�d�ʃt�@�C��(�f��)�������1��o�͂��邩 0=�o�͂��Ȃ�
	int avs_mesh2_interval;	//AVS���b�V���t�@�C��(�f��)�������1��o�͂��邩 0=�o�͂��Ȃ�
	int avs_mesh3_interval;	//AVS���b�V���t�@�C��(�ގ�)�������1��o�͂��邩 0=�o�͂��Ȃ�

	double P_size_AVS;		//�o�͂��闱�q�̃T�C�Y(le�̉��{��)
	double maxT;			//AVS(2)�ɂ�����ő剷�x
	double minT;			//AVS(2)�ɂ�����ŏ����x

	double max_pressure;
	double min_pressure;

	bool flag_modify_density;

	double times_Pa;

	double ave_P_for_FEM_flag;
	bool Ferror;
	bool poise_flag;
	double length;


	//���e���v�Z
	int flag_ELAST;
	int flag_HYPER;
	double hyper_density;
	double c01;
	double c10;
	int flag_wall;
	double r_z_wall;
	double h_dis;
	double h_vis;
	int flag_vis;
	int tension_test;	//�������莎����͗p15/2/8
	int nr_time;

public:
	mpsconfig();
	void Set_length(double len){length=len;}
	double Get_length(){return length;}
	double get_dt(){return dt;}
	int get_step(){return step;}
	int get_current_step(){return current_step;}
	void set_current_step(int cs){current_step=cs;}
	double get_current_time(){return current_time;}
	void set_current_time(const double ct){current_time=ct;}
	int get_avoid_step(){return avoid_step;}
	int get_avoid_step2(){return avoid_step2;}
	int get_avoid_step3(){return avoid_step3;}
	int get_avoid_step4(){return avoid_step4;}
	int get_avoid_step5(){return avoid_step5;}
	int get_avoid_step6(){return avoid_step6;}
	int get_avoid_step7(){return avoid_step7;}

	bool get_cut_x(){return flag_cut_x_movie;}
	bool get_cut_y(){return flag_cut_y_movie;}

	int get_dimension(){return dimension;}
	double get_maxX(){return maxX;}
	double get_minX(){return minX;}
	double get_maxY(){return maxY;}
	double get_minY(){return minY;}
	double get_maxZ(){return maxZ;}
	double get_minZ(){return minZ;}
	
	double get_density(){return MRE_density;}
	double Get_MRE_density(){return MRE_density;}
	double Get_Silicone_density(){return Silicone_density;}
	bool get_modify_density(){return flag_modify_density;}
	double get_nensei(){return nensei;}
	double get_sigma(){return sigma;}
	double get_vis(){return nensei/MRE_density;}	//���S���W��
	double Get_MRE_vis(){return nensei/MRE_density;}
	double Get_Silicone_vis(){return nensei/Silicone_density;}
	double get_Cp(){return Cp;}
	double get_k(){return k;}
	double get_latent_H(){return latent_H;}
	double get_MP(){return MP;}
	double get_CTE(){return CTE;}
	double get_E_m(){return E_m;}
	double get_v_m(){return v_m;}
	double get_E_e(){return E_e;}
	double get_v_e(){return v_e;}

	bool get_Ferror(){return Ferror;}//Fe�G���[
	void set_Ferror(bool ero){Ferror=ero;}
	int get_fluidwidth(){return fluidwidth;}
	double get_distancebp(){return distancebp;}
	double get_wlength(){return wlength;}
	double get_height(){return height;}
	
	double get_re(){return re;}
	double get_re2(){return re2;}
	double get_re3(){return re3;}
	double get_re4(){return re4;}
	double get_re_elastic(){return re_elastic;}
	double get_beta(){return beta;}
	int get_dx(){return dx;}
	double get_times(){return times;}

	int get_motion_interval(){return motion_interval; }

	int get_X_mesh(){return (int)((maxX-minX)/(distancebp*dx)+0.00000000000001);}		//X�������̊i�q�� �ۂߌ덷��h�����߂�0.001�𑫂��Ă���i�H�j
	int get_Y_mesh(){return (int)((maxY-minY)/(distancebp*dx)+0.00000000000001);}		//Y�������̊i�q�� �ۂߌ덷��h�����߂�0.001�𑫂��Ă���
	int get_Z_mesh(){return (int)((maxZ-minZ)/(distancebp*dx)+0.00000000000001);}		//Z�������̊i�q�� �ۂߌ덷��h�����߂�0.001�𑫂��Ă���
	int get_number_of_mesh(){return (int)((maxX-minX)/(distancebp*dx)*(maxY-minY)/(distancebp*dx)*(maxZ-minZ)/(distancebp*dx)+0.001);}//�i�q���FX_mesh*Y_mesh*Z_mesh
	
	int get_surface_tension(){return surface_tension;}
	int get_smooth(){return smooth;}
	int get_SFT(){return SFT;}
	int get_smn(){return smn;}
	int get_smth_sumP(){return smth_sumP;}
	int get_wall_direct(){return wall_direct;}
	int get_non_slip(){return non_slip;}
	int get_suf_modify(){return suf_modify;}
	int get_interpolate_curv(){return interpolate_curv;}

	int get_iteration_count(){return iteration_count;}
	int get_
		(){return Pgrad;}
	int get_ave_Pdirect(){return ave_Pdirect;}
	int get_solution(){return solution;}
	int get_B_of_P(){return B_of_P;}
	double get_w_div(){return w_div;}
	int get_HL_sw(){return HL_sw;}
	int get_dir_for_P(){return dir_for_P;}
	int get_initialP(){return initialP;}
	double get_CGep(){return CGep;}
	int get_omp_P(){return omp_P;}
	int get_negativeP(){return negativeP;}
	int get_minP(){return minP;}
	int get_set_P(){return set_P;}
	int get_artP_sw(){return artP_sw;}
	double get_artP(){return artP;}
	double get_Pgrad_times(){return Pgrad_times;}
	int get_Pgrad_order(){return Pgrad_order;}
	int get_P_AVS(){return P_AVS;}
	
	int get_T_field(){return T_field;}
	int get_insulate(){return insulate;}
	int get_T_laplacian(){return T_laplacian;}
	double get_wall_density(){return wall_density;}
	double get_wall_Cp(){return wall_Cp;}
	double get_wall_k(){return wall_k;}
	double get_roomT(){return roomT;}
	int get_air_cool(){return air_cool;}
	int get_T_expansion(){return T_expansion;}
	int get_buoyant(){return buoyant;}
	double get_TplotZ(){return TplotZ;}
	int get_T_AVS(){return T_AVS;}

	int get_EM_method(){return EM_method;}
	int get_EM_calc_type(){return EM_calc_type;}
	int get_EM_interval(){return EM_interval;}
	int get_region_shape(){return region_shape;}
	double get_XL(){return XL;}
	double get_XR(){return XR;}
	double get_YU(){return YU;}
	double get_YD(){return YD;}
	double get_ZU(){return ZU;}
	double get_ZD(){return ZD;}
	double get_RU(){return RU;}

	int get_mesh_input(){return mesh_input;}
	int get_remesh_sw(){return remesh_sw;}
	double get_boxalpha(){return boxalpha;}
	int get_fine(){return fine;}
	double get_co_fine(){return co_fine;}
	int get_add_points(){return add_points;}
	int get_poly_memory(){return poly_memory;}
	int get_air_layer(){return air_layer;}
	double get_layer_depth(){return layer_depth;}
	int get_mesh_output_interval(){return mesh_output_interval;}
	
	int get_BEM_elm_type(){return BEM_elm_type;}
	int get_BEM_M_solver(){return BEM_M_solver;}
	int get_gauss_N(){return gauss_N;}
	int get_FEM_elm_type(){return FEM_elm_type;}
	int get_FEM_smn(){return FEM_smn;}
	int get_max_DN(){return max_DN;}
	int get_FEMCG(){return FEMCG;}
	double get_CGaccl(){return CGaccl;}
	double get_EMCGep(){return EMCGep;}
	double get_FEMtimes(){return FEMtimes;}
	double get_legend_F(){return legend_F;}
	int get_tree_sw(){return tree_sw;}
	double get_tree_mesh_size(){return tree_mesh_size;}
	int get_p_max(){return p_max;}
	int    get_plot_F_type(){return plot_F_type;}
	double get_legend_F_Pa(){return legend_F_Pa;}

	double get_surface_depth(){return surface_depth;}
	
	double get_V(){return V;}
	int get_V_step(){return V_step;}
	double get_r_perm(){return r_perm;}
	int get_V_con(){return V_con;}
	double get_initial_V(){return initial_V;}
	double get_E_times(){return E_times;}
	int	get_eleforce(){return eleforce;}
	int get_charge(){return charge;}
	int get_plot_E_type()	{return plot_E_type;}
	
	int get_J_input_way(){return J_input_way;}
	double get_J0(){return J0;}
	double get_I0(){return I0;}
	double get_RP(){return RP;}
	double get_ele_conduc(){return ele_conduc;}
	int get_Hz(){return Hz;}
	int get_div_Hz(){return div_Hz;}
	int get_jheat(){return jheat;}
	int get_m_force(){return m_force;}
	int get_NLBHI(){return NLBHI;}
	int get_NLMH(){return NLMH;}
	double get_magnet_H(){return magnet_H;}
	double get_magnet_r(){return magnet_r;}
	double get_magnet_Z(){return magnet_Z;}
	double get_magnet_angle(){return magnet_angle;}
	double get_magnet_B(){return magnet_B;}
	int get_uniform_B_sw(){return uniform_B_sw;}
	int get_magnetic_layers(){return magnetic_layers;}
	double get_uniform_B(){return uniform_B;}
	double get_B_times(){return B_times;}
	int get_plot_B_type(){return plot_B_type;}

    int get_FEM_calc_type()	{return FEM_calc_type;}
	int get_ele_type(){return ele_type;}
	double get_FEMCGep() {return FEMCGep;}
	double get_MRTRep() {return MRTRep;}

	double get_g(){return g;}
	int get_restart(){return restart;}
	int get_autosave(){return autosave;}
	double get_courant(){return courant;}
	int get_modify_position(){return modify_position;}
	int get_vis_calc_type(){return vis_calc_type;}
	int get_wall_adheision(){return wall_adheision;}
	int get_laplacian(){return laplacian;}
	int get_vis_solver(){return vis_solver;}
	int get_initial_u(){return initial_u;}
	int get_temporary_r(){return temporary_r;}
	
	int get_fix_center(){return fix_center;}
	int get_freeon(){return freeon;}
	int get_freeon3sw(){return freeon3sw;}
	int get_surface_judge2(){return surface_judge2;}
	int get_move_prtcl(){return move_prtcl;}
	int get_move_u_dirct(){return move_u_dirct;}
	double get_move_speed(){return move_speed;}
	int get_check_something(){return check_something;}
	int get_set_zero_speed(){return set_zero_speed;}

	int get_model_number(){return model_number;}
	int get_model_set_way(){return model_set_way;}

	double get_R1(){return R1;}

	bool get_FEM_flag(){return switch_FEM;}
	void set_FEM_flag(bool sw){switch_FEM=sw;}
	bool get_vis_flag(){return switch_vis;}

	int get_speed_plot_particle(){return speed_plot_particle;}
	double get_speedtimes(){return speedtimes;}
	int get_speed_face(){return speed_face;}
	double get_speed_face_p(){return speed_face_p;}
	int get_ax_sym_modify(){return ax_sym_modify;}
	int get_flat_speed_plot(){return flat_speed_plot;}
	double get_flat_speed_p(){return flat_speed_p;}
	int get_relative_speed(){return relative_speed;}
	int get_speed_AVS(){return speed_AVS;}
	double get_legend_speed(){return legend_speed;}

	int get_pressure_face(){return pressure_face;}			//3D��͎���speed.dat�̏o�͖� 0=YZ���� 1=XZ
	double get_pressure_face_p(){return pressure_face_p;}		//3D��͎���speed.dat�̏o�͖ʂ̍��W
	double get_pressure_times(){return pressure_times;}
	double get_legend_pressure(){return legend_pressure;}


	int get_M_form(){return M_form;}
	int get_MAX_thread(){return MAX_thread;}
	
	int get_interval(){return interval;}
	int get_AVS(){return AVS;}
	double get_P_size_AVS()	{return P_size_AVS;}

	double get_maxT(){return maxT;}
	double get_minT(){return minT;}
	double get_max_pressure(){return max_pressure;}
	double get_min_pressure(){return min_pressure;}

	double get_particle_mass();
	void modify_density(double coefficient){MRE_density*=coefficient;}

	double get_particle_volume(){ return (get_particle_mass()/get_density());}

	int get_F_interval(){return F_interval;}

	int get_avs_eforce_interval(){return avs_eforce_interval;}
	int get_avs_mesh1_interval(){return avs_mesh1_interval;}
	int	get_avs_mesh2_interval(){return avs_mesh2_interval;}
	int	get_avs_mesh3_interval(){return avs_mesh3_interval;}

	double get_times_Pa(){return times_Pa;}

	double get_ave_P_for_FEM_flag(){return ave_P_for_FEM_flag;}

	bool get_nonlinear_elastic_flag(){return nonlinear_elastic;}
	void set_nonlinear_elastic_flag(bool flag){nonlinear_elastic=flag;}

	void change_step_size(){dt=dt_for_FEM;} //���ԍ��ݕ���FEM�v�Z�p�ɕύX����

	double get_ground_position(){return ground;}
	void set_ground_position(double g){ground=g;}

	bool get_poise_flag(){return poise_flag;}
	void set_poise_flag(bool flg){poise_flag=flg;}



	//���e���v�Z
	int get_flag_ELAST(){return flag_ELAST;}
	int get_flag_HYPER(){return flag_HYPER;}
	double get_hyper_density(){return hyper_density;}
	double get_c10(){return c10;}
	double get_c01(){return c01;}
	int get_flag_wall(){return flag_wall;}
	double get_h_dis(){return h_dis;}
	double get_h_viscousity(){return h_vis;}
	int get_flag_vis(){return flag_vis;}
	int get_nr(){return nr_time;}
};

class elastic; //�O���錾 ���ꂪ�Ȃ��ƃG���[���o��?

#endif