#include "fdct_wrapping_inc.hpp"
#include "nummat.hpp"
#include <iostream>
#include "fdct_wrapping.hpp"
#include "fdct_wrapping_inline.hpp"
#include "math.h"
#include "myfunction.h"
using std::ostream;
using std::istream;
using namespace std;
using namespace fdct_wrapping_ns;
/*
ȫ�����С������
���룺
u_t���ɴ�ϵ����ʾ������
curcoef_t,�����任ϵ���еĴ�ϵ��
m��n�������ݳߴ�
dt������
time����������
�����
�����result
*/
float ** TV_function(float ** u_t, vector< vector<CpxNumMat> > curcoef_t,int m=16,int n=16,double dt=0.002,int time=5);
