class Post
{
	bool plot_B;	//�������x�o��
	bool plot_NF;	//�ߓ_�͏o��
public:
	Post();
	bool get_plot_B(){return plot_B;}
	bool get_plot_NF(){return plot_NF;}
};

//�\���p�t���O�̐ݒ�
Post::Post()
{
	plot_B=true;
	plot_NF=true;
}