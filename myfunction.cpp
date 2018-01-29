//求最小值
#include "myfunction.h"


// 地震数据转换成曲波能处理的复数数据
void segy2curvelet(CpxNumMat& x, float **mysegy,int m,int n )
{
	
    for(int i=0; i<m; i++)
    for(int j=0; j<n; j++)
	x(i,j)=cpx(mysegy[i][j],0);
}

void curvelet2segy(CpxNumMat& x,float **mysegy,int m,int n )
{
	cpx* cpx_tmp;
    cpx c_tmp;
	for(int  i=0; i<m; i++)
    for(int j=0; j<n; j++)
	{   cpx_tmp=x.data();
		c_tmp=cpx_tmp[j*x.m()+i];
		mysegy[i][j]=float(c_tmp.real());

	}

}



int * TotalElements(vector<vector <CpxNumMat> > c)
{
	int *temp;
	temp=new int[c.end()-c.begin()];
	for (int s=0;s<c.end()-c.begin();s++)
	{temp[s]=0;
		for(int w=0;w<c[s].end()-c[s].begin();w++)
			temp[s]+=c[s][w].m()*c[s][w].n();
	}
	return temp;
}

double *SearchThreshold(vector<vector <CpxNumMat> > c,double th_p)
{
	int num=0;int num2=0;
	double *th;
	th=new double[c.end()-c.begin()];
	double *temp;
	int *te;
	te=TotalElements(c);
	for (int s=0;s<c.end()-c.begin();s++)
	{	
		num=0;
		temp=new double[te[s]];
		for(int w=0;w<c[s].end()-c[s].begin();w++)
			for(int i=0;i<c[s][w].m();i++)
				for(int j=0;j<c[s][w].n();j++)
				{
					temp[num]=sqrt(norm(c[s][w].data()[i*c[s][w].n()+j]));num++;num2++;
				}
	 
			
	 sort(temp,temp+te[s]);
	 th[s]=temp[int(ceil((te[s]-1)*(1-th_p)))];
	 //cout<<th[s]<<endl;
	 delete temp;
	
	
	
	}

return th;

}
