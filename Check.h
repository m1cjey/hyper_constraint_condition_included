#pragma once
/*
Check�N���X

�ォ�猩�������Ƃ��ɂ킩��悤�ɉ�͏����E�v�Z�����L�^����N���X

*/
class Check
{
	string checkf;

public:
	Check();
	virtual ~Check();
	////////////////�N�����擾//////////////////
	virtual string Set_y_m_d();

	///////////��͏����t�@�C���o��/////////////
	virtual void Out_put_config();

	/////////////�N�[�������̊m�F///////////////
	virtual void Courant_condition(vector<mpselastic> &PART);
};