#include "fdct_wrapping_inc.hpp"
#include "nummat.hpp"
#include <iostream>
#include "fdct_wrapping.hpp"
#include "fdct_wrapping_inline.hpp"
#include "math.h"
using std::ostream;
using std::istream;
using namespace std;
using namespace fdct_wrapping_ns;

void segy2curvelet(CpxNumMat& x, float **mysegy,int m,int n );//����������ת���������任�����ܹ����������
void curvelet2segy(CpxNumMat& x, float **mysegy,int m,int n );//segy2curvelet�ķ��任
/*
������ֵ����
���룺�����任ϡ��c
	  ��ֵ�ٷֱ�th_p
�������ֵ����
*/
double *SearchThreshold(vector<vector <CpxNumMat> > c,double th_p);
int * TotalElements(vector<vector <CpxNumMat> > c);//���������任ϵ��c�е���Ԫ�ظ���
