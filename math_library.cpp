#include"stdafx.h"

/********************�x�N�g���m�����̌v�Z********************/

//1�m����
double vector_norm1(double *d, int start, int end)
{
	double norm=0.0;
	
	for(int i=start;i<end;i++) norm+=fabs(d[i]);
	//for(int i=m;i<n;i++) norm+=fabs(d[i]);
	//cout<<norm<<endl;
//	cout<<"finished, n:"<<n<<endl;

	return norm;
}

//2�m����
double vector_norm2(double *d, int start, int end)
{
	double sum=0.0;

	for(int i=start;i<end;i++) sum+=d[i]*d[i];
//	cout<<sqrt(sum)<<endl;

	return sqrt(sum);
}

/********************�A���ꎟ�������̌v�Z********************/

//�K�E�X�̏����@
void pivot_gauss(double **A, double *b, int size, bool pivot_check)
{
	int i, j, k, ip;
	double alpha, temp;
	double amax, eps=pow(2.0, -12.0);

	for(k=0;k<size-1;k++)
	{
		amax=fabs(A[k][k]);
		ip=k;

		//�s�{�b�g�I���ƍs����
		if(pivot_check==ON)
		{
			for(i=k+1;i<size;i++)//�s�{�b�g�T��
			{
				if(fabs(A[i][k])>amax)
				{
					amax=fabs(A[i][k]);
					ip=i;
				}
			}

			if(amax<eps){cout<<"�s��͐����ł͂Ȃ��B�v�Z���f"<<endl; return;}	//�������̔���

			//�s����
			if(ip!=k)//�s�{�b�g�̍s�����ڂ��Ă���s�Ɠ����s�Ȃ�����̕K�v�Ȃ�
			{
				for(j=k;j<size;j++)//�e�񂲂ƂɌ�����Ɓij�͗�ԍ��j
				{
					temp=A[k][j];
					A[k][j]=A[ip][j];
					A[ip][j]=temp;
				}
				temp=b[k];//temp=x[k];
				b[k]=b[ip];//x[k]=x[ip];
				b[ip]=temp;//x[ip]=temp;
			}
		}

		//�O�i�����i��k���k+1�s�ȉ��͑S��0�ɂȂ�...�͂��j
		for(i=k+1;i<size;i++)
		{
			alpha=-1*A[i][k]/A[k][k];
			for(j=k+1;j<size;j++)
			{
				A[i][j]+=alpha*A[k][j];
			}
			b[i]+=alpha*b[k];
		}
	}

	//��ޑ��
	b[size-1]/=A[size-1][size-1];
	for(k=size-2;k>=0;k--)//k=size-1�͊��ɋ��߂Ă���
	{
		temp=b[k];
		for(j=k+1;j<size;j++)
		{
			temp=temp-A[k][j]*b[j];
		}
		b[k]=temp/A[k][k];
	}
}

//LU������p���� size:=ROW(�����s�������)
void lu_decomposition(double **A, double *b, size_t size, bool pivot_check)
{
	int i, j, k, ip;
	double alpha, temp;
	double amax;
	int *p=new int[size];

	for(k=0;k<size-1;k++)
	{
		//�s�{�b�g�I��
		amax=fabs(A[k][k]);//�Ίp����
		ip=k;

		for(i=k+1;i<size;i++)
		{
			if(fabs(A[i][k])>amax)
			{
				amax=fabs(A[i][k]);
				ip=i;
			}
		}

		//�������̔���
		if(amax<EPS){
			cout<<"�s��͐����ł͂Ȃ��B�v�Z�I��"<<endl; 
			return;
		}

		p[k]=ip; //ip��z��p[i]�ɕۑ�

		//�s����
		if(ip!=k)
		{
			for(j=k;j<size;j++)
			{
				temp=A[k][j];
				A[k][j]=A[ip][j];
				A[ip][j]=temp;
			}
		}

		//�O�i����
		for(i=k+1;i<size;i++)
		{
			alpha=-A[i][k]/A[k][k];
			A[i][k]=alpha;
			for(j=k+1;j<size;j++)
			{
				A[i][j]+=alpha*A[k][j];
			}
		}
	}

	for(k=0;k<size-1;k++)
	{
		temp=b[k];
		b[k]=b[p[k]];//???
		b[p[k]]=temp;
		//�O�i���
		for(i=k+1;i<size;i++)
		{
			b[i]+=A[i][k]*b[k];
		}
	}

	//��ޑ���E�E�EA[size-1][size-1]!=0�͕ۏ؂���Ă���H�H�H
	b[size-1]/=A[size-1][size-1];
	for(k=size-2;k>=0;k--)
	{
		temp=b[k];
		for(j=k+1;j<size;j++)
		{
			temp-=A[k][j]*b[j];
		}
		b[k]=temp/A[k][k];
	}

	delete [] p;
}

//LU������p����(vector���g��)
void lu_decomposition(vector<vector<double>> &A, vector<double> &b)
{
	int i, j, k, ip;
	size_t size=b.size();
	double alpha, temp;
	double amax;
	vector<int> p(size); //�錾�Ɠ����ɃT�C�Y���m�ۂ��Ă���

	for(k=0;k<size-1;k++)
	{
		//�s�{�b�g�I��
		amax=fabs(A[k][k]);//�Ίp����
		ip=k;

		for(i=k+1;i<size;i++)
		{
			if(fabs(A[i][k])>amax)
			{
				amax=fabs(A[i][k]);
				ip=i;
			}
		}

		//�������̔���
		if(amax<EPS){
			cout<<"�s��͐����ł͂Ȃ��B�v�Z�I��"<<endl; 
			exit(1);
		}

		p[k]=ip; //ip��z��p[i]�ɕۑ�

		//�s����
		if(ip!=k)
		{
			for(j=k;j<size;j++)
			{
				temp=A[k][j];
				A[k][j]=A[ip][j];
				A[ip][j]=temp;
			}
		}

		//�O�i����
		for(i=k+1;i<size;i++)
		{
			alpha=-A[i][k]/A[k][k];
			A[i][k]=alpha;
			for(j=k+1;j<size;j++)
			{
				A[i][j]+=alpha*A[k][j];
			}
		}
	}

	for(k=0;k<size-1;k++)
	{
		temp=b[k];
		b[k]=b[p[k]];//???
		b[p[k]]=temp;
		//�O�i���
		for(i=k+1;i<size;i++)
		{
			b[i]+=A[i][k]*b[k];
		}
	}

	//��ޑ���E�E�EA[size-1][size-1]!=0�͕ۏ؂���Ă���H�H�H
	b[size-1]/=A[size-1][size-1];
	for(k=size-2;k>=0;k--)
	{
		temp=b[k];
		for(j=k+1;j<size;j++)
		{
			temp-=A[k][j]*b[j];
		}
		b[k]=temp/A[k][k];
	}
}


//Cholesky����@��p����
//�Ώ̐���l�s��ɂ̂ݓK�p�\
void modified_cholesky_decomposition(double **A, double *b, const int size)
{
	//�C���R���X�L�[����
	int i, j, k;
	double tmp;

	for(i=2;i<=size;i++)
	{
		for(j=1;j<=i-1;j++)
		{
			tmp=0.0;
			for(k=1;k<=j-1;k++)
			{
				tmp+=A[i][k]*A[k][k]*A[j][k];
			}
			A[i][j]=(A[i][j]-tmp)/A[i][j];
		}
		tmp=0.0;
		for(k=1;k<=j-1;k++)
		{
			tmp+=A[i][k]*A[i][k]*A[k][k];
		}
		A[i][i]=A[i][i]-tmp;
	}

	//�R���X�L�[�����ɂ��s���@
	//LDy=b
	b[1]=b[1]/A[1][1];
	for(i=2;i<=size;i++)
	{
		tmp=0.0;
		for(j=1;j<=i-1;j++)
		{
			tmp+=A[j][j]*A[i][j]*b[j];
		}
		b[i]=(b[i]-tmp)/A[i][i];
	}

	//L^t x=y
	for(i=size-1;i>=1;i--)
	{
		tmp=0.0;
		for(j=i+1;j<=size;j++)
		{
			tmp+=A[j][i]*b[j];
		}
		b[i]=b[i]-tmp;
	}
	//*b�͉��x�N�g��
}

//CG�@�i�������z�@�j
void conjugate_gradient(double **A, double *b, double *x, const int size)
{
	const double EPSILON=1e-8;
	const int MAX_ITERATION=100;
	double *r=new double[size];
	double *p=new double[size];
	double *tmp=new double[size];
	double alpha, beta, work;

	int iteration=0;//������

	//tmp=A*b���v�Z
	matrix_vector_product(A, b, tmp, size, size);

	for(int i=0;i<size;i++)
	{
		p[i]=b[i]-tmp[i];
		r[i]=p[i];
	}

	do{
		//alpha�̌v�Z
		matrix_vector_product(A, p, tmp, size, size);
		
		work=alpha=0.0;
		for(int i=0;i<size;i++) work+=p[i]*tmp[i];//����
		for(int i=0;i<size;i++) alpha+=p[i]*r[i];
		alpha/=work;

		//x_(k+1)��r_(k+1)�̌v�Z
		for(int i=0;i<size;i++) x[i]+=alpha*p[i];
		for(int i=0;i<size;i++) r[i]-=alpha*tmp[i];

		//��������
		if(vector_norm1(r, 0, size)<EPSILON){
			break;//do-while�𔲂���
		}else{
			iteration++;
			beta=0.0;
			for(int i=0;i<size;i++) beta+=r[i]*tmp[i];
			beta*=-1;
			beta/=work;
			for(int i=0;i<size;i++) p[i]=r[i]+beta*p[i];
		}

	}while(iteration<MAX_ITERATION);

	if(iteration==MAX_ITERATION){
		std::cout<<"Iteration reached MAX_ITERATION: "<<MAX_ITERATION<<endl;
		//try-catch
	}else{
		//std::cout<<"Iteration: "<<iteration<<endl;
	}

	delete [] r;
	delete [] p;
	delete [] tmp;
}

//�s��ƃx�N�g���̐�
void matrix_vector_product(double **A, double *b, double *x, int row, int col)
{
	for(int i=0;i<row;i++)
	{
		double temp=0.0;
		for(int j=0;j<col;j++)
		{
			temp+=A[i][j]*b[j];
		}
		x[i]=temp;
	}
}