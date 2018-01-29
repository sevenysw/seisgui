#include "readwrite.h"
float ** Read2dSegy(char filename[],int m,int n)
{
	float **data;
	data=new float *[m];
	for (int i=0;i<m;i++)
		data[i]=new float [n];
	READ(filename,data,m,n,0,1);
	float **data2;
	data2=new float *[m];
	for (int i=0;i<m;i++)
		data2[i]=new float [n];
	for(int i=0;i<m;i++)
		for(int j=0;j<n;j++)
			data2[i][j]=float(data[i][j]);
	return data2;
}
void Write2dSegy(char filename[],float **data2,int m,int n)
{
	float **data;
	data=new float *[m];
	for (int i=0;i<m;i++)
		data[i]=new float [n];
	for(int i=0;i<m;i++)
		for(int j=0;j<n;j++)
			data[i][j]=float(data2[i][j]);
	WRITE(filename,data,m,n,0,1);

}