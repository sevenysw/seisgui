#include <gtk/gtk.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "susgy.h"
#include "as_wave.h"
//include "myrw.h"
#include "ReadSegy.h"
#include "mymath.h"
#include "fdct_wrapping.hpp"
#include "fdct_wrapping_inline.hpp"
#include "myfunction.h"
#include "TV_function.h"

using namespace fdct_wrapping_ns;
//using namespace cv;
using namespace std;


GtkWidget *window_rd2d = NULL;
GtkWidget *image = NULL;

GtkWidget *image2 = NULL;

const int samplebyte = 4; /* float IBM: 4

float IEEE: 4
32-bit fixed: 4
16-bit fixed: 2
*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++*/

float tmp[5];
float tmp_float;
float dt;
int tmp_int, tmp_int2;
int inlineno, xlineno;
int findflag, existflag;
int shotnumcount, shottrccount, shottrccount_judge, alltrccount;
int identnum_p, identnum_n;
/*--------------------------------------------------------------------*/
/*  Shot gather trace headers are stored in bufshotheader[M][K] in binary format;
M=shottrccount;   K=240;
Shot gather trace data are stored in shotdata_trans[M][L] in float format;
M=shottrccount;   L=sample_num;
*/
//Input the data which want to be processed
//float **a;
//a=new float*[shottrccount];
char buf3600[3600];
char buf240[240];
char buf180[180];
char buf4[4];
char buf52[52];
char buf2[2];
char **bufshotheader;
char **bufshotdata;
float **shotdata_trans;    //根据数据*/读取道数和采样点、采样率
//open file
ifstream filein;
ofstream fileout;


void stack_mov_4(char * buf) ;
void stack_mov_2(char * buf) ;
float singleIBM2float(char *buf) ;
int  float2singleIBM(char * buf, float value) ;


char string_Line[100];
char * input_segy_file = new char [100];
char output_segy_file[100];
int  index_key;
int maxshottrc;
int sample_num;
int iter;

float ** data_u;

void readpara(void);
int readheader(void);
float ** readdata(void);
void curvelet_denoise(float *dataf,float *datas, int ns, int ntr,int nbscales,int nbangles_coarse,float ac);
void writedata(void);

static void put_pixel (GdkPixbuf *pixbuf, int x, int y, guchar red, guchar green, guchar blue, guchar alpha) 
{
    int width, height, rowstride, n_channels; 
    guchar *pixels, *p;   
    n_channels = gdk_pixbuf_get_n_channels (pixbuf);  
    g_assert (gdk_pixbuf_get_colorspace (pixbuf) == GDK_COLORSPACE_RGB); 
    g_assert (gdk_pixbuf_get_bits_per_sample (pixbuf) == 8); 
    g_assert (gdk_pixbuf_get_has_alpha (pixbuf));   
    g_assert (n_channels == 4);   
    width = gdk_pixbuf_get_width (pixbuf);  
    height = gdk_pixbuf_get_height (pixbuf); 
   g_assert (x >= 0 && x < width);   
   g_assert (y >= 0 && y < height);  
   rowstride = gdk_pixbuf_get_rowstride (pixbuf); 
   pixels = gdk_pixbuf_get_pixels (pixbuf);    
   p = pixels + y * rowstride + x * n_channels;  
   p[0] = red;  
   p[1] = green; 
   p[2] = blue;  
   p[3] = alpha;
}

void show_data(int flag){

    //load data    
    GdkPixbuf *image_pixbuf;
    int width, height, nchannels, rowstride;
    float max_u, min_u;

    image_pixbuf = gtk_image_get_pixbuf ((GtkImage *)image);
    width = gdk_pixbuf_get_width (image_pixbuf);
    height = gdk_pixbuf_get_height (image_pixbuf);
    nchannels = gdk_pixbuf_get_n_channels(image_pixbuf);
    rowstride = gdk_pixbuf_get_rowstride (image_pixbuf);
    printf("width:%d, hight:%d, nchannels:%d, rowstride:%d\n\n", width,height,nchannels,rowstride);

    max_u = mymax(data_u,shottrccount,sample_num)*0.5;
    min_u = mymin(data_u,shottrccount,sample_num)*0.5;

    printf("\nns:%d, ntr:%d\n", sample_num,shottrccount);
    int i = 0;
    int j = 0;
    int ii,jj;
    guchar *pixels, *p;
    pixels = gdk_pixbuf_get_pixels (image_pixbuf);
    for(i=0; i<height; i++){
        for(j=0; j<width; j++){
            p = pixels + i * rowstride + j * nchannels;
            ii = floor(1.0*i*shottrccount/height);
            jj = floor(1.0*j*sample_num/width);
            
            p[0] = int((data_u[ii][jj]-min_u)/(max_u-min_u)*255);
            p[1] = p[0];//255-p[1];
            p[2] = p[0];//255-p[2];
        }
    }
    if (flag == 1 )
    	gtk_image_set_from_pixbuf (GTK_IMAGE(image),image_pixbuf);
    if (flag == 2)
    	gtk_image_set_from_pixbuf (GTK_IMAGE(image2),image_pixbuf);
}


void denoise(void){

  int i, j;
  int ns, ntr;
  float *datas,*dataf;
  ntr =shottrccount;   
  ns = sample_num; 

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
      datas[i*ns+j]= data_u[i][j];
    }
  }


  float ac=0.0;
  curvelet_denoise(dataf,datas,ns,ntr, nbscales, nbangles_coarse,ac);

  for(i=0;i<ntr;i++){
    for(j=0;j<ns;j++){
      data_u[i][j]=dataf[i*ns+j];
    }
  }
}

void save_clicked(GtkWidget *widget, gpointer data) { 

	GtkWidget *dialog;
	gchar *filename;
	gint res;

	dialog=gtk_file_chooser_dialog_new("SelectFile",NULL,GTK_FILE_CHOOSER_ACTION_SAVE,GTK_STOCK_CANCEL,GTK_RESPONSE_CANCEL,GTK_STOCK_OK,GTK_RESPONSE_ACCEPT,NULL);

	res = gtk_dialog_run (GTK_DIALOG (dialog));
	if (res == GTK_RESPONSE_ACCEPT)
	  {
	    char *filename;
	    GtkFileChooser *chooser = GTK_FILE_CHOOSER (dialog);
	    filename = gtk_file_chooser_get_filename (chooser);
	    printf ("The segy file is:\n" );
	    printf (filename);
	    printf ("\n");
	    g_free (filename);
	  }

	gtk_widget_destroy (dialog);
    writedata();

}

void begin_clicked(GtkWidget *widget, gpointer data) { 

    denoise();
    show_data(2);

}

void load_clicked(GtkWidget *widget, gpointer data) { 

	GtkWidget *dialog;
	gchar *filename;
	gint res;

	dialog=gtk_file_chooser_dialog_new("SelectFile",NULL,GTK_FILE_CHOOSER_ACTION_OPEN,GTK_STOCK_CANCEL,GTK_RESPONSE_CANCEL,GTK_STOCK_OK,GTK_RESPONSE_ACCEPT,NULL);

	res = gtk_dialog_run (GTK_DIALOG (dialog));
	if (res == GTK_RESPONSE_ACCEPT)
	  {
	    char *filename;
	    GtkFileChooser *chooser = GTK_FILE_CHOOSER (dialog);
	    filename = gtk_file_chooser_get_filename (chooser);
	    printf ("The segy file is:\n" );
	    printf (filename);
	    printf ("\n");
	    g_free (filename);
	  }

	gtk_widget_destroy (dialog);
    
    readpara();
    readheader();    
    data_u = readdata();
    show_data(1);

}
void rd2d_clicked(GtkWidget *widget, gpointer data) {  
  
  GtkWidget *button_load;
  GtkWidget *button_begin;
  GtkWidget *button_save;
  GtkWidget *button_stop;


  GtkWidget *hbox;
  GtkWidget *vbox_left;
  GtkWidget *vbox_right;

//window set up

  window_rd2d = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  gtk_window_set_position(GTK_WINDOW(window_rd2d), GTK_WIN_POS_CENTER);
  gtk_window_set_default_size(GTK_WINDOW(window_rd2d), 500, 250);
  gtk_window_set_title(GTK_WINDOW(window_rd2d), "Random 2D");
  gtk_container_set_border_width(GTK_CONTAINER(window_rd2d), 5);

// data setup

  image = gtk_image_new_from_file("redrock.png");
  image2 = gtk_image_new_from_file("redrock.png");

// layout setup

  hbox = gtk_hbox_new(TRUE, 1);
  vbox_left = gtk_vbox_new(TRUE, 1);
  vbox_right = gtk_vbox_new(TRUE, 1);
  gtk_container_add(GTK_CONTAINER(window_rd2d), hbox);

// button setup

  button_load = gtk_button_new_with_label("Load");
  button_begin = gtk_button_new_with_label("Begin");
  button_stop = gtk_button_new_with_label("Stop");
  button_save = gtk_button_new_with_label("Save");

  gtk_box_pack_start(GTK_BOX(hbox), vbox_left, TRUE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(hbox), vbox_right, TRUE, TRUE, 0);

  gtk_box_pack_start(GTK_BOX(vbox_left), image, TRUE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(vbox_left), image2, TRUE, TRUE, 0);

  gtk_box_pack_start(GTK_BOX(vbox_right), button_load, TRUE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(vbox_right), button_begin, TRUE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(vbox_right), button_stop, TRUE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(vbox_right), button_save, TRUE, TRUE, 0);

  g_signal_connect(G_OBJECT(button_load), "clicked", 
      G_CALLBACK(load_clicked), NULL);
  g_signal_connect(G_OBJECT(button_begin), "clicked", 
      G_CALLBACK(begin_clicked), NULL);
  g_signal_connect(G_OBJECT(button_save), "clicked", 
      G_CALLBACK(save_clicked), NULL);

  g_print("clicked\n");

  gtk_widget_show_all(window_rd2d);


}

int main(int argc, char *argv[]) {

  GtkWidget *window;
  GtkWidget *vbox;

  GtkWidget *rd2d;
  GtkWidget *gr2d;
  GtkWidget *rd3d;
  GtkWidget *gr3d;

  gtk_init(&argc, &argv);

  window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  gtk_window_set_position(GTK_WINDOW(window), GTK_WIN_POS_CENTER);
  gtk_window_set_default_size(GTK_WINDOW(window), 230, 250);
  gtk_window_set_title(GTK_WINDOW(window), "CS Denoising");
  gtk_container_set_border_width(GTK_CONTAINER(window), 5);

  vbox = gtk_vbox_new(TRUE, 1);
  gtk_container_add(GTK_CONTAINER(window), vbox);

  rd2d = gtk_button_new_with_label("Random 2D");
  gr2d = gtk_button_new_with_label("Groundroll 2D");
  rd3d = gtk_button_new_with_label("Random 3D");
  gr3d = gtk_button_new_with_label("Groundroll 3D");

  gtk_box_pack_start(GTK_BOX(vbox), rd2d, TRUE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(vbox), gr2d, TRUE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(vbox), rd3d, TRUE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(vbox), gr3d, TRUE, TRUE, 0);


  g_signal_connect(G_OBJECT(rd2d), "clicked", 
      G_CALLBACK(rd2d_clicked), NULL);

  g_signal_connect(G_OBJECT(window), "destroy",
        G_CALLBACK(gtk_main_quit), G_OBJECT(window));

  gtk_widget_show_all(window);

  gtk_main();

  return 0;
}


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

  float s = (b1 & 0x80) >> 7;
  float e = b1 & 0x7f;
  float f = (b2 << 16) + (b3 << 8) + b4;
  float m = f / 16777216.0; /* 2^24 */
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

//  char * buf = new char[4];//MEMORY LEAKAGE
/*NOTE: Here the pointer buf pointed memory won't be freed after this function are executed,
unless we use a free synstax. AND it's won't be freed in the main function unless we do
free or the entire program has been executed done.*/
  long sign;
  long exp;
  float mantissa;
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


void readpara(void){


  std::ifstream fpin("inputpara.txt", std::ios::in);
//  if(fn_in.tellg() == -1){
//    cout << " 参数文件     inputpara.txt 不存在 " << endl ;
//       return false ;
//  }
//  else
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

      cout << "The name of SEGYOUT " << output_segy_file << endl ;
    }
    if(rheader==11)
    {
      index_key =atoi(string_Line);
      cout << "The Key TraceHeader in position  " << index_key << endl ;
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



int readheader(void){

  
//       filein.open(input_segy_file, ios::in | ios::binary); 
  filein.open("133_shot.sgy", ios::in | ios::binary);
//  fileout.open(output_segy_file, ios::out | ios::binary); 
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
//  getchar();
  shottrccount_judge = 0;
  alltrccount = 0;
  shotnumcount = 0;
  identnum_n = 0;
}



float ** readdata(void){
  
  int i, j, k, l;
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
//      cout<<" Traces_all "<<alltrccount<<" identnum_p "<<identnum_p<<" identnum_n "<<identnum_n<<endl;
//      getchar();
}
shottrccount_judge = shottrccount;

// define dynamic shot data array, to store shot data in float;
shotdata_trans = new float *[shottrccount];
for (i = 0; i < shottrccount; i++) {
  shotdata_trans[i] = new float[sample_num];
}
// transform bufshotdata[][](single IBM) to shotdata_trans[][](float);

for (i = 0; i < shottrccount; i++) {
  for (j = 0; j < sample_num; j++) {

    for (k = 0; k < 4; k++) {
      buf4[k] = bufshotdata[i][j * samplebyte + k];
    }
    shotdata_trans[i][j] = (float) singleIBM2float(buf4);
//        cout<<shotdata_trans[i][j]<<endl;
//        getchar();
  }
}

//FREE MEMORY
for (i = 0; i < maxshottrc; i++) {
  delete[] bufshotdata[i];
}
delete[] bufshotdata;

cout << "Reading shot " << identnum_p << " completed! " << endl;

return shotdata_trans;

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
myth=SearchThreshold(c,0.1);


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
          data2d[i][j]=0;//datas[i*ns+j];
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

			float2singleIBM( buf4, (float) data_u[i][j]);
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
		delete[] data_u[i];
	}
	delete[] data_u;

//		if ((alltrccount)%10000==1){
	cout << " Traces_all " << alltrccount << "   Shot " << identnum_p
	<< " Shot count " << shotnumcount + 1 << endl;
//		getchar();
//cout<<setiosflags(ios::fixed)<<setw(5)<<setprecision(1);
//cout<<100.0*countsgy/counttxt<<"% Processed..."<<endl;
//		}

	shotnumcount++;
}