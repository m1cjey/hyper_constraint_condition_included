////initial for MPS
using namespace std;

double kernel(double r,double dis); //�d�݊֐�
double kernel2(double r,double dis,double d); //�d�݊֐�
double kernel3(double r,double dis);
double kernel4(double r,double dis);
double initial_pnd(double r,int dimension,int calc_type); //�����z�u�ɂ�����e���q�̗��q�����x
double calclambda(mpsconfig &CON); //���v���V�A���p�Ɍv�Z�֐�
void input_particle_data(mpsconfig *CON, vector<mpselastic> &PART1, vector<mpselastic> &PART, int t); //���q�f�[�^�ǂݎ��֐�
void courant_elastic(mpsconfig *CON,vector<mpselastic> &PART,int fluid_number,int t,double *dt,double mindis,double Umax,double *g); //
double get_volume(mpsconfig *CON);
void output_particle_density(mpsconfig *CON,vector<mpselastic> &PART,int fluid_number,double n0,int particle_number,int t);
void check_initial_position(mpsconfig *CON,vector<mpselastic> &PART);
///
void culan(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int t,double *dt,double mindis,double Umax);
///

//���q���J�E���g�֐� & ���ёւ�
void calc_numbers_of_particles_and_change_the_order(mpsconfig *CON, vector <mpselastic> &PART, int *fluid_number,int *hyper_number,int *magnetic_number,int *out,int *order_sw);	//main_another�d�l	15/2/10
//void calc_numbers_of_particles_and_change_the_order(mpsconfig *CON, vector <mpselastic> &PART, int *fluid_number,int *out,int *order_sw);	//main�d�l	15/2/10
int check_position(mpsconfig *CON,vector<mpselastic> &PART, int fluid_number, int *particle_number);

//���q�@�ɂ����郂�f���쐬
void writedata(ofstream fp, int id, double x, double y,double z, int type,int surface,double val, ofstream gnu, double vx, double vy, double vz, double P, double h, int toBEM);
void set_initial_placement_using_MD(mpsconfig *CON,int *particle_number);

//�i�q(index)�����֌W
void reload_INDEX(mpsconfig &CON, vector<mpselastic> &PART, int *INDEX);
void reload_INDEX2(mpsconfig *CON, vector<mpselastic> &PART, int **MESH);

//�\�ʔ���
void freeon(mpsconfig &CON, vector<mpselastic> &PART, double n0_4,int *INDEX,int **MESH,double *mindis, int t);
void freeon(elastic &ELAST, vector<mpselastic> &PART, double n0_4,int *INDEX,int **MESH, double *mindis, int t);
void freeon2(mpsconfig &CON,vector<mpselastic> &PART,int particle_number,double n0_4,int *INDEX,int **MESH,double *mindis,int fluid_number,int out);

void calc_neighbor_relation(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,double n0_4,int fluid_number,int out);
void surface_judge2(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number);
void surface_judge2_old(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number);
void surface_judge2_new(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number);
void surface_judge3(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number);

//�@���x�N�g��
double pnd_for_direct(mpsconfig &CON,vector<mpselastic> &PART,double x,double y,double z,double R,int i);
void direct_f(mpsconfig &CON,vector<mpselastic> &PART,int i,double *direct[DIMENSION]);

//�t�@�C���o��
void plot_speed(mpsconfig &CON ,vector<mpselastic> &PART,int particle_number,int fluid_number);
void post_processing(mpsconfig &CON,vector<mpselastic> &PART, elastic &ELAST, int fluid_number,int particle_number,double dt,double Umax,int t,double TIME,double **F);
void particle_movie_AVS(mpsconfig &CON,vector<mpselastic> &PART,elastic &ELAST, int fluid_number,int particle_number,int t,double T);
void particle_movie_AVS2(mpsconfig &CON,int t,vector<mpselastic> &PART,double TIME,int fluid_number, int particle_number, double **F);
void post_processing3(mpsconfig &CON,vector<mpselastic> &PART,int fluid_number,int particle_number,int t,double TIME);
void plot_pressure_acceleration(mpsconfig &CON, vector<mpselastic> &PART);
void plot_stress_distribution(mpsconfig &CON, vector<mpselastic> &PART);
void plot_shear_force(mpsconfig &CON, vector<mpselastic> &PART);
void plot_strainrate_acceleration(mpsconfig &CON, vector<mpselastic> &PART);
void plot_distortion_rate(mpsconfig &CON, vector<mpselastic> &PART);
void plot_residual_acceleration(mpsconfig &CON, vector<mpselastic> &PART, double **F);
///
void output_alldata_AVS(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,double dt,double Umax,double mindis,int t,double TIME, unsigned int time0,int count_avs);
///


//�s���@
void CG_method(mpsconfig *CON,double *r,double *P,double *AP,double *val,int *ind,int *ptr,int pn,double *X,int *countN,double EP);
void iccg(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X,double *r,double *P,double EP,int *count2);
void gauss(double *matrix,double *B,int N);
void jacobi(double **matrix,double *B,int N);

//�S��
void visterm_negative(mpsconfig *CON,vector<mpselastic> &PART,double *laplacian[DIMENSION],double n0,double lambda,int fluid_number,int particle_number,double dt,int t); //�A�I�Ȍv�Z
void calc_viscous_term(mpsconfig *CON,vector<mpselastic> &PART,int fluid_number,int particle_number,double dt,double N0,double *laplacian[DIMENSION],double lambda,int t); //�S�����̌v�Z�֐�
///
void u_laplacian_f(mpsconfig *CON,vector<mpsparticle> &PART,double *laplacian[DIMENSION],double n0,double lamda,int fluid_number,double dt);
///

//�\�ʒ���
void smoothing(mpsconfig *CON,vector<mpselastic> &PART,int particle_number,int fluid_number,double *potential[DIMENSION],double n0);
void surface_tension1(mpsconfig *CON,vector<mpselastic> &PART,int fluid_number,double *potential[DIMENSION],int particle_number);
///
void calc_surface_tension(mpsconfig *CON,vector<mpselastic> &PART,int fluid_number,double dt,int particle_number,double n0,double **potential);
void plot_ST(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *potential[DIMENSION],int t);
void plot_ST_each(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *potential[DIMENSION],int t);
void plot_PND(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *PND4,int t);
void plot_PND_each(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *PND4,int t);
void plot_F(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *F[DIMENSION],int t);
void plot_speed_each(mpsconfig *CON ,vector<mpsparticle> &PART,int particle_number,int fluid_number,int t);
///

/////���̑��x����шʒu����
void renewal_u_and_r_in_positive(mpsconfig *CON,vector<mpselastic> &PART,int fluid_number,int t,double dt,double *Umax,double **potential,double **laplacian,double *g,double **previous_Un,double **F);

//���͌v�Z
void negative1(mpsconfig *CON,vector<mpselastic> &PART,int fluid_number,int particle_number,int out,int t,double dt,double lambda,double N0,double *PND2,double n0,double **Un);
void pressure(mpsconfig *CON,vector<mpselastic> &PART,int fluid_number,int particle_number,int out,double dt,int t,double lambda,double N0,double *PND2,double n0);
void calc_P_main(mpsconfig *CON,vector<mpselastic> &PART,int particle_number,double lambda,double N0,double dt,double *PND2,double n0,int fluid_number,int out,int pn);
void set_N0_and_PND2(mpsconfig *CON,vector<mpselastic> &PART,int particle_number,int fluid_number,double *N0,int out);
void set_Dirichlet_P(mpsconfig *CON,vector<mpselastic> &PART,int fluid_number,int flag,double *Dirichlet_P);
void output_dirichlet_vector_files(mpsconfig *CON,vector<mpselastic> &PART,int fluid_number,int flag,double *Dirichlet_P);
double Dndt(mpsconfig *CON,vector<mpselastic> &PART,int i,double n0);
void plot_P(mpsconfig *CON ,vector<mpselastic> &PART,int particle_number,int t,int fluid_number);
void output_pressuer_avs(mpsconfig *CON,vector<mpselastic> &PART,int t,int particle_number,int fluid_number);
void calc_Pgradient(mpsconfig *CON,vector<mpselastic> &PART,int particle_number,int fluid_number,double **reU,double n0,double dt);
void set_minP(mpsconfig *CON,vector<mpselastic> &PART,int fluid_number,double *minP);
void plot_Pgradient(mpsconfig *CON,vector<mpselastic> &PART,int fluid_number,double **Pgrad);
void P_gradient3(mpsconfig *CON,vector<mpselastic> &PART,int fluid_number,double *direct[DIMENSION],double dt,double **reU);
void calc_WLSM_P_D3_order1(mpsconfig *CON,vector<mpselastic> &PART,double *matrix, double *B,double **P_grad,int i,int N,double *minP);
void calc_WLSM_P_D3_order2(mpsconfig *CON,vector<mpselastic> &PART,double *matrix, double *B,double **P_grad,int i,int N,double *minP);
void P_gradient_MPS(mpsconfig *CON,vector<mpselastic> &PART,int i,double n0,double *P_grad);
void P_gradient4(mpsconfig *CON,vector<mpselastic> &PART,double dt,int fluid_number,double *reU[3]);
///
void negative1_twice(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,int out,int t,double dt,double lamda,double N0,double *PND2,double n0,double **Un,double n0_4);
void plot_P_each(mpsconfig *CON ,vector<mpsparticle> &PART,int particle_number,int t,int fluid_number);
void modify_reU(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,int t,double dt,double n0,double *udiv, double **reU);
void calc_inverse_matrix(mpsconfig *CON,vector<mpsparticle> &PART,int N, double *matrix);

///

//�A��͌�̑��x�E���W�X�V
void modify_u_and_x_after_Pcalc(mpsconfig *CON,vector<mpselastic> &PART,int fluid_number,double **reU,double dt,double **Un);


//���x���U
double divergence(mpsconfig *CON,vector<mpselastic> &PART,int i,double n0);
double divergence2(mpsconfig *CON,vector<mpselastic> &PART,int i);
void return_X_for5N(double *matrix,int N,double *B1,double *B2,double *dudx,double *dudy);
double calc_WLSM_divu_D3_order1(mpsconfig *CON,vector<mpselastic> &PART,double *matrix, double *B1,double *B2,double *B3,int i,int N);
double calc_WLSM_divu_D3_order1_2(mpsconfig *CON,vector<mpselastic> &PART,double *matrix, double *B1,double *B2,double *B3,int i,int N);
double calc_WLSM_divu_D3_order2(mpsconfig *CON,vector<mpselastic> &PART,double *matrix, double *B1,double *B2,double *B3,int i,int N);
double calc_WLSM_divu_D3_order2_2(mpsconfig *CON,vector<mpselastic> &PART,double *matrix, double *B1,double *B2,double *B3,int i,int N);

//���q�@���̑�
void modify_position(mpsconfig *CON,vector<mpselastic> &PART,int fluid_number,double dt);

//���x��
void plot_T(mpsconfig *CON ,vector<mpselastic> &PART,int particle_number,double *T,double height);
void calc_Temperature(mpsconfig *CON,vector<mpselastic> &PART,int fluid_number,int particle_number,double n0,double lambda,double dt,int t);
void output_temperature_avs(mpsconfig *CON,vector<mpselastic> &PART,int t,int particle_number,int fluid_number,double *T,double height);

void move_particle(mpsconfig *CON,vector<mpselastic> &PART,int particle_number,int fluid_number,double dt);

//��Ώ̍s���@
void BiCGStab_method(mpsconfig *CON,double *r,int pn,double *X,int *countN,double EP,double **A);
//void BiCGStab2_method(mpsconfig *CON,double *r,int pn,double *X,int *countN,double EP,double **A);

//�����x���֐�
double get_Legendre(int l,int m,double x);
double Legendre(int l,int m,double x);
void get_derivative_of_spherical_harmonics_with_theta(int l, int m, double COStheta, double fai,double *real, double *imaginary,double SIN2);
double Legendre_small_L(int l,double x);
double Legendre_big_L(int l,double x);
double Legendre_derivative_small(int l,int m,double x);

//�t�@�C���������ݗp
void WriteFile();

//�f���[�j����
void delaun3D(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int KTJ,int KTE,double rax,double ray,double raz,int *node_num,int *nelm,int FINE_sw,double rrm);
void data_avs(int node_number,int nelm,vector <point3D> &NODE,vector <element3D> &ELEM,int ktj,double *val,mpsconfig &CON);
void C_Fluix_data_avs(int node_number,int nelm,vector <point3D> &NODE,vector <element3D> &ELEM,int ktj,double **B,mpsconfig &CON,int t);
void data_avs2(mpsconfig &CON,int node,int nelm,vector <point3D> &NODE,vector <element3D> &ELEM,int ktj,int t);
void data_avs3(int node,int nelm,vector <point3D> &NODE,vector <element3D> &ELEM,mpsconfig &CON);
void t_data_avs3(int node,int nelm,vector <point3D> &NODE,vector <element3D> &ELEM,mpsconfig &CON,int t);
void data_avs4(mpsconfig &CON,int node_number,vector <point3D> &NODE,vector <element3D> &ELEM,int nelm,int *kv);
void delaun3D_main(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int KTJ,int KTE,int *node_num,int *nelm,int FINE_sw);
void memorize_static_NODE_and_ELEM(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,vector<point3D> &static_NODE,vector<element3D> &static_ELEM,int node,int nelm);
void set_jnb3D(vector<point3D> &NODE, vector<element3D> &ELEM,int node,int nelm,int *jnb);
void set_nei3D(vector<point3D> &NODE, vector<element3D> &ELEM,int node,int nelm,int *jnb,int **nei);
double volume3D(vector<point3D> &NODE,int ia,int ib,int ic,int ip);
void sphere3D(vector<point3D> &NODE,vector<element3D> &ELEM,int ia,int ib,int ic,int ip,int i);
void poly3D(vector<point3D> &NODE,vector<element3D> &ELEM,int *iv,int *kv,int ip,int *nelm,mpsconfig *CON);
void FINE3D(vector<point3D> &NODE,vector<element3D> &ELEM,int KTJ,int KTE,int *node,int *nelm,mpsconfig *CON,double rrm,int startID);
void fill3D(vector<point3D> &NODE,vector<element3D> &ELEM,int nelm);

//MPSTOFEM3D
void MPS_TO_FEM3D(mpsconfig &CON, int *node_num, vector<point3D> &NODE, vector<mpselastic> &PART, int fluid_number, int particle_number);
void make_dense_cube_region(mpsconfig *CON,vector<point3D> &NODE,int *node, int *divN,double regionX[2],double regionY[2],double regionZ[2]);

//FEM3D
void FEM3D(mpsconfig &CON, vector<mpselastic> &PART, double **F, int *static_node, int *static_nelm,vector<point3D> &static_NODE, vector<element3D> &static_ELEM, int t, double TIME, int fluid_number);
void FEM3D_calculation(mpsconfig &CON,int *static_node,int *static_nelm,vector<point3D> &static_NODE, vector<element3D> &static_ELEM,int t,double TIME,vector<mpselastic> &PART,int fluid_number,int particle_number,double dt,double **F);
int calc_matrix_width(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,int *jnb,int **nei);
void arrange_matrix(int pn,int *NUM,int **ROW,double **G);
void Incomplete_Cholesky_Decomposition(mpsconfig &CON,double *L,double *D1,double *matrix,int *ptr2,int *ind2,int pn,int num2);
void ICCG3D2(mpsconfig &CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X);
void parallel_ICCG3D2(mpsconfig &CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X);
void check_matrix_symmetry(int pn,int *NUM,int **ROW,double **G);
void MRTR(mpsconfig &CON,double *B,int pn,double *X,double *val,int *ind,int *ptr);
void ICMRTR(mpsconfig &CON,double *B,int pn,double *X,double *val,int *ind,int *ptr);

void initial_model_input(mpsconfig *CON, int *particle_number, double *TIME);

//TetGen�p
//�߂�l���w�肵�Ă�����Ɠ��������ǂ����m�F���Ă��ǂ�
void TetGenInterface(mpsconfig &CON, vector<mpselastic> &PART, double **F, int fluid_number,double dt,int t,int particle_number,double n0,double TIME);
void usingTetGen(mpsconfig &CON,vector<mpselastic> &PART, double **F, int fluid_number,int particle_number,int t,double TIME);
void plot_F(mpsconfig &CON, vector<mpselastic> &PART, int fluid_number, int t);
void plot_F_log(mpsconfig &CON, vector<mpselastic> &PART, int fluid_number,int t);

//deraun3D�p
void data_avs3(int node,int nelm,vector <point3D> &NODE,vector <element3D> &ELEM,mpsconfig &CON,int t);
void data_avs4(mpsconfig *CON,int node_number,vector <point3D> &NODE,vector <element3D> &ELEM,int nelm,int *kv);

//FEM�p�̃t���O
void check_FEM_flag(mpsconfig &CON, elastic &ELAST, double ave_P);


//���e��
void calc_hyper(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t,double **F);


//�S�����v�Z
void calc_vis_f(mpsconfig &CON,vector<mpselastic>PART,vector<hyperelastic>&HYPER,vector<hyperelastic2>HYPER1,int rigid_number,int t);
void calc_spl_f(mpsconfig &CON,vector<mpselastic>PART,vector<hyperelastic2>&HYPER1,int hyper_number);

//���͏o�͊֐�
void force_movie_AVS(mpsconfig *CON,int t,vector<mpselastic> &PART,int particle_number,double m);

//���C���[�����g�@
void Magnetic_Moment_Method(mpsconfig &CON,vector<mpselastic>&PART,double **F,double n0,double lamda,int fluid_number,int particle_number, double current_time, int t);
void Magnetic_Moment_Methodv2(mpsconfig &CON,vector<mpselastic> &PART,double **F,double n0,double lamda,int fluid_number,int particle_number, double current_time, int t);
