
#include "TV_function.h"
const int shot_num=1;  

//const double th=0.0005;
const double th=0.005;
float ** TV_function(float ** u_t, vector< vector<CpxNumMat> > curcoef_t,int m,int n,double dt,int time)
{
	int nbscales=6;int nbangles_coarse=8;
	//nbscales=int(ceil(log10(double(min(m,n)))/log10(double(2)))-3);
	int ac=0;
	double b=0.00000000001;
	float **u;
	u=new float *[m+2];
	for(int i=0;i<(m+2);i++)
	{
		u[i]=new float [n+2];
		for(int j=0;j<(n+2);j++)
			u[i][j]=0;
		
	}

		float **g;
	g=new float *[m];
	for(int i=0;i<m;i++)
	{
		g[i]=new float [n];
		for(int j=0;j<n;j++)
		g[i][j]=0;
	}

	double temp=0;double temp1=0;double temp2=0;
	vector< vector<CpxNumMat> > c;  //vector<int> extra;
	CpxNumMat x(m,n);
    
	CpxNumMat y(x); clear(y);
	for (int i=1;i<(m+1);i++)
		for(int j=1;j<(n+1);j++)
		{
			u[i][j]=u_t[i-1][j-1];
			g[i-1][j-1]=0;
		}
	for(int kk=1;kk<=time;kk++)
		//compute g(k)
		{
			printf("kk=%d\n",kk);
			for(int i=1;i<(m+1);i++)
				for(int j=1;j<(n+1);j++)
					g[i-1][j-1]=(2*u[i][j]-u[i+1][j]-u[i][j+1])*pow((pow(u[i][j+1]-u[i][j],2.0)+pow(u[i+1][j]-u[i][j],2.0)+b),-0.5)+
								(u[i][j]-u[i-1][j])*pow((pow(u[i-1][j+1]-u[i-1][j],2.0)+pow(u[i-1][j]-u[i][j],2.0)+b),-0.5)+
								(u[i][j]-u[i][j-1])*pow((pow(u[i+1][j-1]-u[i][j-1],2.0)+pow(u[i][j]-u[i][j-1],2.0)+b),-0.5);
						
			
	segy2curvelet(x,g,m,n);
	fdct_wrapping(m, n, nbscales, nbangles_coarse, ac, x, c);
	cpx cpx_tmp;
	//有错误用的是curcoef而不是c
	for(int s=0;s<(c.end()-c.begin()-1);s++)
		for(int w=0;w<(c[s].end()-c[s].begin());w++)
			for(int i=0;i<c[s][w].m();i++)
			  for(int j=0;j<c[s][w].n();j++)
			  { cpx_tmp=curcoef_t[s][w](i,j);
				temp1=cpx_tmp.real();
				temp2=cpx_tmp.imag();
				//temp=sqrt(temp1*temp1+temp2+temp2);
				temp=sqrt(norm(cpx_tmp));
			    c[s][w].data()[i*c[s][w].n()+j] =cpx(temp1*(th>temp),temp2*(th>temp));
			    
			  }
	
	ifdct_wrapping(m, n, nbscales, nbangles_coarse, ac, c, y);
	curvelet2segy(y,g,m,n);

	for (int i=1;i<(m+1);i++)
		for(int j=1;j<(n+1);j++)
		{
			u[i][j]=u[i][j]-dt*g[i-1][j-1];
			
		}


	}
	
	for (int i=1;i<(m+1);i++)
	for(int j=1;j<(n+1);j++)
	{
		u_t[i-1][j-1]=u[i][j];
		
	}



	return u_t;

}
