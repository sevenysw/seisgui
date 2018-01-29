#define PI 3.1415926 

/* allocate a 3-d array */
//void ***alloc3 (size_t n1, size_t n2, size_t n3, size_t size)
void ***alloc3 (int n1, int n2, int n3, int size)
{
        int i3,i2;
        void ***p;

        if ((p=(void***)malloc(n3*sizeof(void**)))==NULL)
                return NULL;
        if ((p[0]=(void**)malloc(n3*n2*sizeof(void*)))==NULL) {
                free(p);
                return NULL;
        }
        if ((p[0][0]=(void*)malloc(n3*n2*n1*size))==NULL) {
                free(p[0]);
                free(p);
                return NULL;
        }

        for (i3=0; i3<n3; i3++) {
                p[i3] = p[0]+n2*i3;
                for (i2=0; i2<n2; i2++)
                        p[i3][i2] = (char*)p[0][0]+size*n1*(i2+n2*i3);
        }
        return p;
}

/* free a 3-d array */
void free3 (void ***p)
{
        free(p[0][0]);
        free(p[0]);
        free(p);
}

/* allocate a 2-d array */
void **alloc2 (int n1, int n2, int size)
{
        int i2;
        void **p;

        if ((p=(void**)malloc(n2*sizeof(void*)))==NULL) 
                return NULL;
        if ((p[0]=(void*)malloc(n2*n1*size))==NULL) {
                free(p);
                return NULL;
        }
        for (i2=0; i2<n2; i2++)
                p[i2] = (char*)p[0]+size*n1*i2;
        return p;
}

/* free a 2-d array */
void free2 (float **p)
{
        free(p[0]);
        free(p);
}

/* allocate a 1-d array */
void *alloc1 (int n1, int size)
{
	void *p;

	if ((p=(void *)malloc(n1*size))==NULL)
		return NULL;
	return p;
}

/* free a 1-d array */
void free1 (void *p)
{
	free(p);
}


/************************************************************************
* tforce_ricker -- Compute the	time response of a source function as
*	a Ricker wavelet with peak frequency "fpeak" Hz.	
*************************************************************************
*************************************************************************
* Input: 
*	int nt		number of time step
*	float dt 	time step
*	float fpeak	peak frequency of the Ricker wavelet
*	
*************************************************************************
* Output: float *tforce		source array
*
*************************************************************************/

void tforce_ricker(int nt, float *tforce, float dt, float fpeak)
{
  int	it;
  float	t1, t0; 

  t0=1.0/fpeak;
 // t0 = (float)nt*dt/2.0;

  for (it=0; it<nt; it++) 
  {
    t1=it*dt;
    tforce[it] = exp(-PI*PI*fpeak*fpeak*(t1-t0)*(t1-t0))*
      (1.0-2.*PI*PI*fpeak*fpeak*(t1-t0)*(t1-t0));
  }
}



typedef struct
{
       float real;
       float imag;
} COMPLEX;
/**********************************************************************

DFT.C - SOURCE CODE FOR DISCRETE FOURIER TRANSFORM FUNCTIONS

dft     Discrete Fourier Transform
idft    Inverse Discrete Fourier Transform
fft     In-place radix 2 decimation in time FFT
ifft    In-place radix 2 decimation in time inverse FFT
log2    Base 2 logarithm
***********************************************************************/

/***********************************************************************

dft - Discrete Fourier Transform
  
This function performs a straight DFT of N points on an array of
complex numbers whose first member is pointed to by Datain.  The
output is placed in an array pointed to by Dataout.

void dft(COMPLEX *Datain, COMPLEX *Dataout, int N)

*************************************************************************/

void dft(COMPLEX *Datain, COMPLEX *Dataout, int N)
{
    int i,k,n,p;
    static int nstore = 0;      /* store N for future use */
    static COMPLEX *cf;         /* coefficient storage */
    COMPLEX *cfptr,*Dinptr;
    float arg;

/* Create the coefficients if N has changed */

    if(N != nstore) {
        if(nstore != 0) free((char *) cf);    /* free previous */

        cf = (COMPLEX  *) calloc(N, sizeof(COMPLEX));
        if (!cf) {
            printf("\nUnable to allocate memory for coefficients.\n");
            exit(1);
        }

        arg = 8.0*atan(1.0)/N;
        for (i=0 ; i<N ; i++) {
            cf[i].real = (float)cos(arg*i);
            cf[i].imag = -(float)sin(arg*i);
        }
    }

/* Perform the DFT calculation */

    /*printf("\n");*/
    for (k=0 ; k<N ; k++) {

        Dinptr = Datain;
        Dataout->real = Dinptr->real;
        Dataout->imag = Dinptr->imag;
        Dinptr++;
        for (n=1; n<N; n++) {

        p = (int)((long)n*k % N);
            cfptr = cf + p;         /* pointer to cf modulo N */

            Dataout->real += Dinptr->real * cfptr->real
                             - Dinptr->imag * cfptr->imag;

            Dataout->imag += Dinptr->real * cfptr->imag
                             + Dinptr->imag * cfptr->real;
            Dinptr++;
        }
        /*if (k % 32 == 31) printf("*");*/
        if (k % 32 == 31) {        };
        Dataout++;          /* next output */
    }
    /*printf("\n");*/
}

/***********************************************************************

idft - Inverse Discrete Fourier Transform
  
This function performs an inverse DFT of N points on an array of
COMPLEX numbers whose first member is pointed to by Datain.  The
output is placed in an array pointed to by Dataout.
It returns nothing.

void idft(COMPLEX *Datain, COMPLEX *Dataout, int N)

*************************************************************************/

void idft(COMPLEX *Datain, COMPLEX *Dataout, int N)
{
    int i,k,n,p;
    static int nstore = 0;      /* store N for future use */
    static COMPLEX *cf;         /* coefficient storage */
    COMPLEX *cfptr,*Dinptr;
    float arg;

/* Create the coefficients if N has changed */

    if(N != nstore) {
        if(nstore != 0) free((char *) cf);    /* free previous */

        cf = (COMPLEX  *) calloc(N, sizeof(COMPLEX));
        if (cf == 0) {
            printf("\nUnable to allocate memory for coefficients.\n");
            exit(1);
        }

/* scale stored values by 1/N */
        arg = 8.0*atan(1.0)/N;
        for (i=0 ; i<N ; i++) {
            cf[i].real = (float)(cos(arg*i)/(float)N);
            cf[i].imag = (float)(sin(arg*i)/(float)N);
        }
    }

/* Perform the DFT calculation */

    /*printf("\n");*/
    for (k=0 ; k<N ; k++) {

        Dinptr = Datain;
        Dataout->real = Dinptr->real * cf[0].real;
        Dataout->imag = Dinptr->imag * cf[0].real;
        Dinptr++;
        for (n=1; n<N; n++) {

        p = (int)((long)n*k % N);
            cfptr = cf + p;         /* pointer to cf modulo N */

            Dataout->real += Dinptr->real * cfptr->real
                             - Dinptr->imag * cfptr->imag;

            Dataout->imag += Dinptr->real * cfptr->imag
                             + Dinptr->imag * cfptr->real;
            Dinptr++;
        }
        /*if (k % 32 == 31) printf("*");*/
        Dataout++;          /* next output */
    }
    /*printf("\n");*/
}

/**************************************************************************

fft - In-place radix 2 decimation in time FFT

Requires pointer to COMPLEX array, x and power of 2 size of FFT, m
(size of FFT = 2**m).  Places FFT output on top of input COMPLEX array.

void fft(COMPLEX *x, int m)

*************************************************************************/

void fft(COMPLEX *x, int m)
{
    static COMPLEX *w;           /* used to store the w COMPLEX array */
    static int mstore = 0;       /* stores m for future reference */
    static int n = 1;            /* length of fft stored for future */

    COMPLEX u,temp,tm;
    COMPLEX *xi,*xip,*xj,*wptr;

    int i,j,k,l,le,windex;

    float arg,w_real,w_imag,wrecur_real,wrecur_imag,wtemp_real;

    if(m != mstore) {

/* free previously allocated storage and set new m */

        if(mstore != 0) free(w);
        mstore = m;
        if(m == 0) return;       /* if m=0 then done */

/* n = 2**m = fft length */

        n = 1 << m;
        le = n/2;

/* allocate the storage for w */

        w = (COMPLEX *) calloc(le-1,sizeof(COMPLEX));
        if(!w) {
            printf("\nUnable to allocate COMPLEX W array\n");
            exit(1);
        }

/* calculate the w values recursively */

        arg = 4.0*atan(1.0)/le;         /* PI/le calculation */
        wrecur_real = w_real = cos(arg);
        wrecur_imag = w_imag = -sin(arg);
        xj = w;
        for (j = 1 ; j < le ; j++) {
            xj->real = (float)wrecur_real;
            xj->imag = (float)wrecur_imag;
            xj++;
            wtemp_real = wrecur_real*w_real - wrecur_imag*w_imag;
            wrecur_imag = wrecur_real*w_imag + wrecur_imag*w_real;
            wrecur_real = wtemp_real;
        }
    }

/* start fft */

    le = n;
    windex = 1;
    for (l = 0 ; l < m ; l++) {
        le = le/2;

/* first iteration with no multiplies */

        for(i = 0 ; i < n ; i = i + 2*le) {
            xi = x + i;
            xip = xi + le;
            temp.real = xi->real + xip->real;
            temp.imag = xi->imag + xip->imag;
            xip->real = xi->real - xip->real;
            xip->imag = xi->imag - xip->imag;
            *xi = temp;
        }

/* remaining iterations use stored w */

        wptr = w + windex - 1;
        for (j = 1 ; j < le ; j++) {
            u = *wptr;
            for (i = j ; i < n ; i = i + 2*le) {
                xi = x + i;
                xip = xi + le;
                temp.real = xi->real + xip->real;
                temp.imag = xi->imag + xip->imag;
                tm.real = xi->real - xip->real;
                tm.imag = xi->imag - xip->imag;             
                xip->real = tm.real*u.real - tm.imag*u.imag;
                xip->imag = tm.real*u.imag + tm.imag*u.real;
                *xi = temp;
            }
            wptr = wptr + windex;
        }
        windex = 2*windex;
    }            

/* rearrange data by bit reversing */

    j = 0;
    for (i = 1 ; i < (n-1) ; i++) {
        k = n/2;
        while(k <= j) {
            j = j - k;
            k = k/2;
        }
        j = j + k;
        if (i < j) {
            xi = x + i;
            xj = x + j;
            temp = *xj;
            *xj = *xi;
            *xi = temp;
        }
    }
}

/**************************************************************************

ifft - In-place radix 2 decimation in time inverse FFT

Requires pointer to complex array, x and power of 2 size of FFT, m
(size of FFT = 2**m).  Places inverse FFT output on top of input
frequency domain COMPLEX array.

void ifft(COMPLEX *x, int m)

*************************************************************************/

void ifft(COMPLEX *x, int m)
{
    static COMPLEX *w;           /* used to store the w complex array */
    static int mstore = 0;       /* stores m for future reference */
    static int n = 1;            /* length of ifft stored for future */

    COMPLEX u,temp,tm;
    COMPLEX *xi,*xip,*xj,*wptr;

    int i,j,k,l,le,windex;

    float arg,w_real,w_imag,wrecur_real,wrecur_imag,wtemp_real;
    float scale;

    if(m != mstore) {

/* free previously allocated storage and set new m */

        if(mstore != 0) free(w);
        mstore = m;
        if(m == 0) return;       /* if m=0 then done */

/* n = 2**m = inverse fft length */

        n = 1 << m;
        le = n/2;

/* allocate the storage for w */

        w = (COMPLEX *) calloc(le-1,sizeof(COMPLEX));
        if(!w) {
            printf("\nUnable to allocate complex W array\n");
            exit(1);
        }

/* calculate the w values recursively */

        arg = 4.0*atan(1.0)/le;         /* PI/le calculation */
        wrecur_real = w_real = cos(arg);
        wrecur_imag = w_imag = sin(arg);  /* opposite sign from fft */
        xj = w;
        for (j = 1 ; j < le ; j++) {
            xj->real = (float)wrecur_real;
            xj->imag = (float)wrecur_imag;
            xj++;
            wtemp_real = wrecur_real*w_real - wrecur_imag*w_imag;
            wrecur_imag = wrecur_real*w_imag + wrecur_imag*w_real;
            wrecur_real = wtemp_real;
        }
    }

/* start inverse fft */

    le = n;
    windex = 1;
    for (l = 0 ; l < m ; l++) {
        le = le/2;

/* first iteration with no multiplies */

        for(i = 0 ; i < n ; i = i + 2*le) {
            xi = x + i;
            xip = xi + le;
            temp.real = xi->real + xip->real;
            temp.imag = xi->imag + xip->imag;
            xip->real = xi->real - xip->real;
            xip->imag = xi->imag - xip->imag;
            *xi = temp;
        }

/* remaining iterations use stored w */

        wptr = w + windex - 1;
        for (j = 1 ; j < le ; j++) {
            u = *wptr;
            for (i = j ; i < n ; i = i + 2*le) {
                xi = x + i;
                xip = xi + le;
                temp.real = xi->real + xip->real;
                temp.imag = xi->imag + xip->imag;
                tm.real = xi->real - xip->real;
                tm.imag = xi->imag - xip->imag;             
                xip->real = tm.real*u.real - tm.imag*u.imag;
                xip->imag = tm.real*u.imag + tm.imag*u.real;
                *xi = temp;
            }
            wptr = wptr + windex;
        }
        windex = 2*windex;
    }            

/* rearrange data by bit reversing */

    j = 0;
    for (i = 1 ; i < (n-1) ; i++) {
        k = n/2;
        while(k <= j) {
            j = j - k;
            k = k/2;
        }
        j = j + k;
        if (i < j) {
            xi = x + i;
            xj = x + j;
            temp = *xj;
            *xj = *xi;
            *xi = temp;
        }
    }

/* scale all results by 1/n */
    scale = (float)(1.0/n);
    for(i = 0 ; i < n ; i++) {
        x->real = scale*x->real;
        x->imag = scale*x->imag;
        x++;
    }
}






/***********************************************************************
log2 - base 2 logarithm

Returns base 2 log such that i = 2**ans where ans = log2(i).
if log2(i) is between two values, the larger is returned.

int log2(unsigned int x)

*************************************************************************/

int log2(unsigned int x)
{
    unsigned int mask,i;

    if(x == 0) return(-1);      /* zero is an error, return -1 */

    x--;                        /* get the max index, x-1 */

    for(mask = 1 , i = 0 ; ; mask *= 2 , i++) {
        if(x == 0) return(i);   /* return log2 if all zero */
        x = x & (~mask);        /* AND off a bit */
    }
}

/***********************************************************************
FFT_2D
************************************************************************/
void FFT_2D(COMPLEX **x, int mi, int mj, COMPLEX *w)
{
	 //static COMPLEX *w;
	 int i, j, ni, nj;

	 ni = 1<<mi;
	 nj = 1<<mj;

	 for(i=0;i<ni;i++) {
		 fft(x[i],mj);
	 }

/* for(i=0;i<2;i++) {
		 for(j=0;j<nj;j++) {
		    printf("sss i=%d,j=%d,f2darry=%f\n",i,j,x[i][j].real);
		 }
		 getchar();
	}
	 printf("2\n");

	 for(j=0; j<nj; j++) {
	    printf("i=%d,f2darry=%f\n",j, x[0][j].real);
	}
	 getchar();*/

	 //w = (COMPLEX *) calloc(ni,sizeof(COMPLEX));
	 for(j=0;j<nj;j++) {
		 for(i=0;i<ni;i++) {
			 w[i].real = (float)x[i][j].real;
			 w[i].imag = (float)x[i][j].imag;

		 }
		 fft(w,mi);
		 for(i=0;i<ni;i++) {
			 x[i][j].real = (float)w[i].real;
			 x[i][j].imag = (float)w[i].imag;

		 }
	 }

    /* for(i=0;i<ni;i++) {
		 for(j=0;j<nj;j++) {
		    printf("i=%d,j=%d,f2darry=%f\n",i,j,x[i][j].real);
		 }
		 getchar();
	}*/
}

/********************************************************************
IFFT_2D
********************************************************************/
void IFFT_2D(COMPLEX **x, int mi, int mj, COMPLEX *w)

{
	 int i, j, ni, nj;

	 ni = 1<<mi;
	 nj = 1<<mj;

	 for(j=0;j<nj;j++) {
		 for(i=0;i<ni;i++) {
			 w[i].real = (float)x[i][j].real;
			 w[i].imag = (float)x[i][j].imag;

		 }
		 ifft(w,mi);
		 for(i=0;i<ni;i++) {
			 x[i][j].real = (float)w[i].real;
			 x[i][j].imag = (float)w[i].imag;

		 }  
	 }

	 for(i=0;i<ni;i++) {
		 ifft(x[i],mj);		 
	 }
 
}



	 




