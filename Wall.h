#pragma once
//
//���q�N���X
//

class Wall
{
	double box_node[8][3];
	double box_eleID[6][3];
	
public:
	
	////////////�{�b�N�X���b�V����/////////////
	virtual void Set_boxmeshuWall(double size[3], double core[3]);
	///////////////////////////////////////////
	//////////////////���ʕ�///////////////////
	virtual void Set_flatWall();
	///////////////////////////////////////////
	
};