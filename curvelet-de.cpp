/*Curvelet+TV             */

/* Written by HIT MathGeo */

/* 2017.09                */



//include system headers
#include <fstream>
#include <iostream>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <stdio.h>

//include user headers
//#include "opencv2/imgproc/imgproc.hpp"
//#include "opencv2/highgui/highgui.hpp"
#include "susgy.h"
#include "as_wave.h"
#include "fdct_wrapping.hpp"
#include "fdct_wrapping_inline.hpp"
#include "myfunction.h"
#include "TV_function.h"


//using namespace
using namespace fdct_wrapping_ns;
//using namespace cv;
using namespace std;

//initialize functions
void curvelet_denoise(float *datac,float *datas, int ns, int ntr,int nbscales,int nbangles_coarse,float ac);
double *SearchThreshold(vector<vector <CpxNumMat> > c,double th_p);
void WriteSegy(float *array, char outputname[],int ns, int ntr, float dt, float dx);//ntr:trace number. ns:sample number.
/*++++++++++++++++++++++++++++++++++++++++++++++++++*/


//int maxshottrc = 200;/*maximum trace number in a shot gather;
//IMPORTANT NOTICE: must be greater than shot traces number!
//*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*IMPORTANT:
Input SEG-Y trace data should be in single-precision IBM data format;
*/
const int samplebyte = 4; /* float IBM: 4
float IEEE: 4
32-bit fixed: 4
16-bit fixed: 2
*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++*/

void stack_mov_4(char * buf) {
/*stack_mov_4(char *buf):
4 bytes of char buf[4], change the stack memory sequence;
used for reading binary file.*/

	char buf_temp;

	buf_temp = buf[0];
	buf[0] = buf[3];
	buf[3] = buf_temp;

	buf_temp = buf[1];
	buf[1] = buf[2];
	buf[2] = buf_temp;
}

void stack_mov_2(char * buf) {
/*stack_mov_2(char *buf):
2 bytes of char buf[2], change the stack memory sequence;
used for reading binary file.*/

	char buf_temp;

	buf_temp = buf[0];
	buf[0] = buf[1];
	buf[1] = buf_temp;
}

float singleIBM2float(char *buf) {
/*Transform 32-bits single precision IBM to float;

For 32-bit IBM,
V=pow(-1.0,S) * M * pow(16.0,E-64);
M= C+F = 0+F = F;
Important:
# include <math.h>
*/

	unsigned char b1 = buf[0];
	unsigned char b2 = buf[1];
	unsigned char b3 = buf[2];
	unsigned char b4 = buf[3];

	double s = (b1 & 0x80) >> 7;
	double e = b1 & 0x7f;
	double f = (b2 << 16) + (b3 << 8) + b4;
double m = f / 16777216.0; /* 2^24 */
//      cout<<" s "<<s<<" e "<<e<<" f "<<f<<" m "<<m<<" ";

	int int_s = s;
	int int_e = e;

	if (s == 0 && e == 0 && f == 0) {
		return 0.0;
	} else {
		return (float) pow(-1.0, int_s) * m * pow(16.0, int_e - 64);
	}
}

int  float2singleIBM(char * buf, float value) {

//	char * buf = new char[4];//MEMORY LEAKAGE
/*NOTE: Here the pointer buf pointed memory won't be freed after this function are executed,
unless we use a free synstax. AND it's won't be freed in the main function unless we do
free or the entire program has been executed done.*/
	long sign;
	long exp;
	double mantissa;
	long fraction;

	int int_sign;
	int int_exp;

	sign = (value < 0 ? 1 : 0);
	int_sign = sign;
	value = 1.0 * value * pow(-1.0, int_sign);

	if (value > 0) {

		exp = (long) (log(value) / (4.0 * log(2.0)) + 65);
		int_exp = exp;
		mantissa = value * 1.0 / pow(16.0, int_exp - 64);
		fraction = (long) (mantissa * 16777216.0);

	} else {
		exp = 0;
		fraction = 0;
	}

buf[0] = (sign << 7) & 0x80 | exp & 0x7f; //NOTE: | stands for "OR";
buf[1] = (fraction >> 16) & 0xff;
buf[2] = (fraction >> 8) & 0xff;
buf[3] = fraction & 0xff;

return 0;
}


char string_Line[100];
char * input_segy_file = new char [100];
char output_segy_file[100];
int  index_key;
int maxshottrc;
int sample_num;
int iter;
int ns,ntr;
void readpara(void){


	std::ifstream fpin("inputpara.txt", std::ios::in);
//	if(fn_in.tellg() == -1){
//		cout << " 参数文件     inputpara.txt 不存在 " << endl ;
//       return false ;
//	}
//	else
	{
		cout << "Paramter file: inputpara.txt "<< endl ;
	}
	int rheader=0;
	while( !fpin.eof())
	{
		fpin.getline (string_Line,100);
		rheader++;

		if(rheader==5)
		{
// input_segy_file=string_Line;//这个地方不知道怎么改快捷???这样子只能读一个首字符
			strcpy(input_segy_file, string_Line);
			cout << "The name of  SEGYIN " << input_segy_file << endl ;
		}
		if(rheader==8)
		{
//output_segy_file= string_Line;
			strcpy(output_segy_file, string_Line);

			cout << "The name of出SEGYOUT " << output_segy_file << endl ;
		}
		if(rheader==11)
		{
			index_key =atoi(string_Line);
			cout << "输The Key TraceHeader in position  " << index_key << endl ;
		}
		if(rheader==14)
		{
maxshottrc =atoi(string_Line);//将字符转换成数字
cout <<"The Max Trace Num for all shots " << maxshottrc << endl ;
}

if(rheader==17)
{
	sample_num =atoi(string_Line)  ;
	cout << "The number of sampling points " << sample_num << endl ;
}
if(rheader==20)
{
	iter =atoi(string_Line)  ;
	cout << "Maxium  Number of Iterations of  TV  "<< iter << endl ;
}

}
fpin.close();
}

double tmp[5];
double tmp_double;
float dt;
int tmp_int, tmp_int2;
int inlineno, xlineno;
int findflag, existflag;
int shotnumcount, shottrccount, shottrccount_judge, alltrccount;
int identnum_p, identnum_n;
/*--------------------------------------------------------------------*/
/*	Shot gather trace headers are stored in bufshotheader[M][K] in binary format;
M=shottrccount;		K=240;
Shot gather trace data are stored in shotdata_trans[M][L] in double format;
M=shottrccount;		L=sample_num;
*/
//Input the data which want to be processed
//double **a;
//a=new double*[shottrccount];
char buf3600[3600];
char buf240[240];
char buf180[180];
char buf4[4];
char buf52[52];
char buf2[2];
char **bufshotheader;
char **bufshotdata;
double **shotdata_trans;		 //根据数据*/读取道数和采样点、采样率
//open file
ifstream filein;
ofstream fileout;

int readheader(void){


//       filein.open(input_segy_file, ios::in | ios::binary); 
	filein.open("133_shot.sgy", ios::in | ios::binary);
//	fileout.open(output_segy_file, ios::out | ios::binary); 
	fileout.open("133_shot_out.sgy", ios::out | ios::binary);
	if(filein.tellg() == -1){
		cout<<" Error: Failed to open file"<<endl;
		return 1;
	}

	cout << "The name of  SEGYIN " << input_segy_file << endl ;
//Reading sgyfile's fileheader: 3200byte EBDIC & 400byte binary header;
//Writing to output;
	filein.read(buf3600, 3600);
	fileout.write(buf3600, 3600);

	cout << "Now reading sgy file one shot by oneshot. Please wait..." << endl
	<< endl;
	cout<<"in.peek   "<<filein.peek()<<"  in.tellg  "<<filein.tellg()<<"  out.tellp  "<<fileout.tellp()<<endl;
//	getchar();
	shottrccount_judge = 0;
	alltrccount = 0;
	shotnumcount = 0;
	identnum_n = 0;
}

void readdata(void){
	int i, j, k, l, ns,ntr;


	cout<<"testtesttest"<<endl;

//define dynamic array to store one shot gather;
	bufshotheader = new char *[maxshottrc];
	for (i = 0; i < maxshottrc; i++) {
		bufshotheader[i] = new char[240];
	}
	bufshotdata = new char *[maxshottrc];
	for (i = 0; i < maxshottrc; i++) {
		bufshotdata[i] = new char[sample_num * samplebyte];
	}

//Processing on shot gather;  分选单炮，并
	shottrccount = 0;
	for (k = 0; k < maxshottrc; k++) {

//in case of insuffcient maxshottrc, and not completed. NEED to be improved here.
//			if(shotnumcount != 0 && shottrccount_judge != maxshottrc && identnum_n == identnum_p){
//				cout<<" Error: Defined Maximum Shot Traces Number (maxshottrc) is too small! "<<endl;
//				exit(0);
//			}

//in case of file end;
		if (filein.peek() == EOF) {
			break;
		}

		identnum_p = identnum_n;

//Read trace headers;
		filein.read(bufshotheader[shottrccount], 240);

//Extract Shot Number;
for (i = 0; i < 4; i++) {/*Shot Number: 9-12byte  道头：FFID在9-12*/
		buf4[i] = bufshotheader[shottrccount][index_key-1 + i];
	}
	stack_mov_4(buf4);
	identnum_n = *((int *) (buf4));


for (i = 0; i < 2; i++) {/*采样率: 9-12byte  道头：FFID在9-12*/
	buf2[i] = bufshotheader[shottrccount][116 + i];
}

stack_mov_2(buf2);
dt = *((float *) (buf2))/1000000;

if (k != 0 && identnum_n != identnum_p) {
	filein.seekg(-240, ios::cur);
	break;
}
//Read trace values;
filein.read(bufshotdata[shottrccount], sample_num * samplebyte);
shottrccount++;
alltrccount++;
//			cout<<" Traces_all "<<alltrccount<<" identnum_p "<<identnum_p<<" identnum_n "<<identnum_n<<endl;
//			getchar();
}
shottrccount_judge = shottrccount;

// define dynamic shot data array, to store shot data in double;
shotdata_trans = new double *[shottrccount];
for (i = 0; i < shottrccount; i++) {
	shotdata_trans[i] = new double[sample_num];
}
// transform bufshotdata[][](single IBM) to shotdata_trans[][](double);
for (i = 0; i < shottrccount; i++) {
	for (j = 0; j < sample_num; j++) {

		for (k = 0; k < 4; k++) {
			buf4[k] = bufshotdata[i][j * samplebyte + k];
		}
		shotdata_trans[i][j] = (double) singleIBM2float(buf4);
//				cout<<shotdata_trans[i][j]<<endl;
//				getchar();
	}
}

//FREE MEMORY
for (i = 0; i < maxshottrc; i++) {
	delete[] bufshotdata[i];
}
delete[] bufshotdata;

cout << "Reading shot " << identnum_p << " completed! " << endl;
cout << "Now start denoising process of shot " << identnum_p << " ..."
<< endl;
}
float *datas,*dataf;
void denoise(void){

	int i, j, k, l;
// float **data2d;
// Mat src,des;
//     ifstream infile;
//      infile.open("2dgatherlb.segy",ios::binary);
	ntr =shottrccount;   
	ns = sample_num; 
	cout<<"trace="<<ntr<<"samplenumber="<<ns<<endl; 
//	cout<<ntr<<"   "<<ns<<endl;getchar(); 
	datas = (float *) alloc1(ntr*ns,sizeof(float));
	dataf = (float *) alloc1(ntr*ns,sizeof(float));
	for(i=0;i<ntr*ns;i++) {
		datas[i]=0.0;
		dataf[i]=0.0;
	}
	srand48( (long)time(NULL) );
// int nbscales=(int)(log10((float) ntr)/log10(2)-3);
	int nbscales=6;
	cout<<"scales="<<nbscales<<endl;
//      int nbangles_coarse=((int)(pow(2.0,(floor(nbscales/2.0)))));
	int nbangles_coarse=8;
//  printf("nbangles_coarse=%d\n",nbangles_coarse);
//read segy data and add noise
//     infile.seekg(0,ios::beg);    
	for(i=0;i<ntr;i++){
		for(j=0;j<ns;j++){
			datas[i*ns+j]= shotdata_trans[i][j];
		}
	}

//show noise data with opencv
/*  src = Mat(ntr, ns, CV_32F, datas);
normalize(src,des,1,0,NORM_MINMAX,-1);
namedWindow( "noise", WINDOW_NORMAL );
imshow( "noise", des );*/



//curvelet denoise
//ERROR: Segmentation fault!
	float ac=0.0;
	curvelet_denoise(dataf,datas,ns,ntr, nbscales, nbangles_coarse,ac);

//show denoised data with opencv
/* src = Mat(ntr, ns, CV_32F, dataf);
normalize(src,des,1,0,NORM_MINMAX,-1);
namedWindow( "test", WINDOW_NORMAL );
imshow( "test", des );
waitKey(0);*/

	for(i=0;i<ntr;i++){
		for(j=0;j<ns;j++){
			shotdata_trans[i][j]=dataf[i*ns+j];
		}
	}
}

void writedata(void){
	int i, j, k, l;
	cout << "Processing finished, and now write shot " << identnum_p
	<< " to disk..." << endl;
// transform shotdata_trans[][](double) to bufshotdata[][](single IBM);
// redefine of bufshotdata;

	bufshotdata = new char *[shottrccount];
	for (i = 0; i < shottrccount; i++) {
		bufshotdata[i] = new char[sample_num * samplebyte];
	}
//transform process;
	for (i = 0; i < shottrccount; i++) {
		for (j = 0; j < sample_num; j++) {

			float2singleIBM( buf4, (float) shotdata_trans[i][j]);
			for (k = 0; k < 4; k++) {
				bufshotdata[i][j * samplebyte + k] = buf4[k];
				buf4[k]=0x00;
			}

		}
	}
//write to outputfile;
	for (k = 0; k < shottrccount; k++) {
//Write trace headers;
		fileout.write(bufshotheader[k], 240);
//Write trace values; 
		fileout.write(bufshotdata[k], sample_num * samplebyte);
	}

//FREE MEMORY
	for (i = 0; i < maxshottrc; i++) {
		delete[] bufshotheader[i];
	}
	delete[] bufshotheader;

	for (i = 0; i < shottrccount; i++) {
		delete[] bufshotdata[i];
	}
	delete[] bufshotdata;

	for (i = 0; i < shottrccount; i++) {
		delete[] shotdata_trans[i];
	}
	delete[] shotdata_trans;

//		if ((alltrccount)%10000==1){
	cout << " Traces_all " << alltrccount << "   Shot " << identnum_p
	<< " Shot count " << shotnumcount + 1 << endl;
//		getchar();
//cout<<setiosflags(ios::fixed)<<setw(5)<<setprecision(1);
//cout<<100.0*countsgy/counttxt<<"% Processed..."<<endl;
//		}

	shotnumcount++;
}

//main function
int main(int argc, char* argv[])
{
	cout<<"Hello"<<endl;
	clock_t ck0,ck1;
	ck0=clock();
//stage
/*++++++++++++++++++++++读入参数文本++++++++++++++++++*/
	readpara();
	readheader();


	while (filein.peek() != EOF) {
		readdata();
		denoise();
		writedata();
	}

//shotnumcount++; 

	cout << "Processing finished. Now closing files..." << endl << endl;

	filein.close();
	fileout.close();

//		cout << "***************Completed!*****************" << endl << endl;
	return 0;


}

void curvelet_denoise(float *dataf,float *datas, int ns, int ntr,int nbscales,int nbangles_coarse,float ac)
{
//preparate data as complex data 
	int m=ntr;
	int n=ns;
//int nbscales=6;
// int nbangles_coarse=8;
	CpxNumMat e(m,n);
	CpxNumMat x(m,n);
	for(int i=0; i<m; i++)
		for(int j=0; j<n; j++)
			x(i,j) = cpx(datas[i*n+j] , 0);

//forward curvelet transform
vector< vector<CpxNumMat> > c;  //vector<int> extra;
fdct_wrapping(m, n, nbscales, nbangles_coarse, ac, x, c);
vector < vector<CpxNumMat> > c2(c);

//compute threshold, reserve big value
double temp1, temp2, temp;
double *myth;

myth=new double[c.end()-c.begin()];
myth=SearchThreshold(c,0.2);


for (int s=0;s<(c.size());s++)
	for(int w=0;w<(c[s].size());w++)
		for(int i=0;i<c[s][w].m();i++)
			for(int j=0;j<c[s][w].n();j++)
			{
				temp1=c[s][w](i,j).real();
				temp2=c[s][w](i,j).imag();
				temp=sqrt(temp1*temp1+temp2*temp2);
				c2[s][w](i,j)=cpx(temp1*(temp>myth[s]),temp2*(temp>myth[s]));
			}



//inverse transform 
			CpxNumMat y(x); clear(y);
			ifdct_wrapping(m, n, nbscales, nbangles_coarse, ac, c2, y);

//1D to 2D
			int i,j;
			float **data2d,**data2dd;
			data2d=new float *[ntr];
			data2dd=new float *[ntr];
			for (i=0;i<ntr;i++)

			{ 
				data2d[i]=new float [ns];
				for(j=0;j<ns;j++)
				{
					data2d[i][j]=datas[i*ns+j];
				}
			}

//return
			curvelet2segy(y,data2d,m,n);

//Tv
			double dt=0.002;
			int iter=1;
			data2dd=TV_function(data2d,c2,m,n,dt,iter);



//2D to 1D
			for (i=0;i<ntr;i++)
				for(j=0;j<ns;j++)
				{
					dataf[i*ns+j]=data2dd[i][j];
				}


			}

void WriteSegy(float *array, char outputname[],int ns, int ntr, float dt, float dx)//ntr:trace num. ns:sample num.
{
	int ii,jj;
	segy TraceHead;
	TraceHead.dt=(unsigned short)(dt*1e6);
	TraceHead.ns=(unsigned short)ns;
	FILE *pSegy=fopen(outputname,"wb");
	fseek(pSegy,3600L,0);
	for(ii=0;ii<ntr;ii++)
	{
		TraceHead.tracf=ii+1;
		TraceHead.offset=ii*(int)dx;
		fwrite(&TraceHead,sizeof(segy),1,pSegy);
		for(jj=0;jj<ns;jj++)
		{
			fwrite(&array[ii*ns+jj],sizeof(float),1,pSegy);
		}
		cout<<endl;
	}
	fclose(pSegy);
}

