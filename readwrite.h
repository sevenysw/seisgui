#ifndef R_W
#define R_W
#include <iostream>
#include <fstream>
#include "susgy.h"
using namespace std;

int READ(char filename[],float **data,int nx,int n0,int volumehead,int tracehead)
////filename(�ļ���),data�����ݣ�,int nx�������꣩,int n0�������꣩,int volumehead����ͷ��һ��Ϊ0����ʾ�޾�ͷ��Ϣ��,int tracehead����ͷ��һ��Ϊ1����ʾ�е�ͷ��Ϣ��)

{
	float *volhead=new float [900];
	segy sgy;
	sgy.ntr=nx;
    sgy.ns=n0;
	sgy.dt=1000;
	int i;
	ifstream infile(filename,ios_base::binary);
	if(infile)
	{
		i=0;
		if(volumehead==1)infile.read((char*)volhead,3600);
		while(!infile.eof()&&i<nx)
		{
			if(tracehead==1)infile.read((char*)&sgy,sizeof(segy));
			infile.read((char*)data[i],sizeof(float)*n0);	
			i++;
		}
		infile.close();
		return 1;
	}
	else
	{
		cout<<"cannot open input file!!!!!!!"<<endl;
		return 0;
	}
}
int WRITE(char filename[],float **data,int nx,int n0,int volumehead,int tracehead)
{
	float *volhead=new float [900];
	segy sgy;
	sgy.ntr=nx;
    sgy.ns=n0;
	sgy.dt=1000;
	int i;
	ofstream outfile(filename,ios_base::binary);
	if(outfile)
	{
		i=0;
		if(volumehead==1)outfile.write((char*)volhead,3600);
		while(i<nx)
		{
			if(tracehead==1)outfile.write((char*)&sgy,sizeof(segy));
			outfile.write((char*)data[i],sizeof(float)*n0);	
			i++;
		}
		outfile.close();
		return 1;
	}
	else
	{
		cout<<"cannot open output file!!!!!!!"<<endl;
		return 0;
	}
	delete [] volhead;
}
#endif