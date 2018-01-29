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

void segy2curvelet(CpxNumMat& x, float **mysegy,int m,int n );//将地震数据转换成曲波变换函数能够处理的类型
void curvelet2segy(CpxNumMat& x, float **mysegy,int m,int n );//segy2curvelet的反变换
/*
查找阈值函数
输入：曲波变换稀疏c
	  阈值百分比th_p
输出：阈值数组
*/
double *SearchThreshold(vector<vector <CpxNumMat> > c,double th_p);
int * TotalElements(vector<vector <CpxNumMat> > c);//计算曲波变换系数c中的总元素个数
