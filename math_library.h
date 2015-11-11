//�A���ꎟ�������̉�@
void pivot_gauss(double **A, double *b, const int size, bool pivot_check);//�K�E�X�̏����@
void lu_decomposition(double **A, double *b, size_t size, bool pivot_check);//LU������@
void lu_decomposition(vector<vector<double>> &A, vector<double> &b);

//Choresky����
void modified_cholesky_decomposition(double **A, double *b, int size);

//�������z�@
void conjugate_gradient(double **A, double *b, double *x, int size);
void matrix_vector_product(double **A, double *b, double *x, int row, int col);

//�x�N�g���m����
double vector_norm1(double*v, const int start, const int end);
double vector_norm2(double*v, const int start, const int end);
