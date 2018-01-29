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
全变分最小化函数
输入：
u_t，由大系数表示的数据
curcoef_t,曲波变换系数中的大系数
m，n地震数据尺寸
dt：步长
time：迭代次数
输出：
最后结果result
*/
float ** TV_function(float ** u_t, vector< vector<CpxNumMat> > curcoef_t,int m=16,int n=16,double dt=0.002,int time=5);
