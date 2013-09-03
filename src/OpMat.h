/*
 *  OpMat.h
 *
 *  Created by Alejandro Quiroz on 3/20/11.
 *
 */

#include "libMatInv.h"  

static double maxarg1,maxarg2; 
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
(maxarg1) : (maxarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

/*======================================================
 =======================================================
 Funcion para extraer una submatriz
 ========================================================
 ======================================================*/

void SubMat(double **data,int n, int *indices, int lenInd, double **Ndata)
{
	int i,j,k;
	for(i=0;i<n;i++)
	{
		for(j=0;j<lenInd;j++)
		{
			k=indices[j];
			Ndata[i][j]=data[i][k];
		}
		
	}
}

/*======================================================
========================================================
 Funcion para calcular la suma de una matriz
========================================================
======================================================*/

void MatSum(double **Xmat, double **Ymat, int rX, int cX, double **Result)
{
	int i, j;	
	for (i=0; i<rX; i++) 
	{
		for (j=0; j<cX; j++) 
		{
			Result[i][j]=Xmat[i][j]+Ymat[i][j];
		}
	}
}

/*======================================================
 =======================================================
 Funcion para calcular la resta de una matriz
========================================================
======================================================*/

void MatRes(double **Xmat, double **Ymat, int rX, int cX,double **Result)
{
	int i, j;
	for (i=0; i<rX; i++) 
	{
		for (j=0; j<cX; j++) 
		{
			Result[i][j]=Xmat[i][j]-Ymat[i][j];
		}
	}
}

/*======================================================
========================================================
 Funcion para calcular la transpuesta de una matriz
=========================================================
=======================================================*/
void MatTrans(double **XMat, int rX, int cX,double **Resu)
{
	int i, j;
	
	for(i=0;i<cX;i++)
	{
		for(j=0;j<rX;j++)
		{
			Resu[i][j]=XMat[j][i];
		}
	}
}


/*======================================================
========================================================
 Funcion para calcular la multiplicacion de matrices
=========================================================
=======================================================*/

void MatMult(double **Xmat, double **Ymat, int rX, int cX, int cY, double **Result)
{
	int i, j, k;
	double suma;
	
	for (i=0; i<rX; i++) 
	{
		for (j=0; j<cY; j++) 
		{
			suma=0;
			for (k=0; k<cX; k++) 
			{
				suma=suma+Xmat[i][k]*Ymat[k][j];
			}
			Result[i][j]=suma;	 
		}
	}
}
/*======================================================
========================================================
 Funcion para calcular la inversa de una matriz
========================================================
======================================================*/

void MatInv(double **A, int n, double **Result)
{
	int i,j;
	int *indx;
	double *col;
	
	col=array1M(n);
	indx=arrayInt1M(n);
	
	LUdcmp(A,n,indx);
	for(j=0;j<n;j++)
	{
		for(i=0;i<n;i++)
		{
			col[i]=0;
		}
		col[j]=1;
		LUbksb(A, n, indx, col);
		for(i=0;i<n;i++)
		{
			Result[i][j]=col[i];
		}
	}
	free(col);
	free(indx);	
}

/*======================================================
========================================================
 Funcion para calcular Cholesky
========================================================
======================================================*/

void MatChol(double **A,int n)
{
	int i,j,k;
	double sum;
	double *p;
	
	p=array1(n);
	
	for(i=0;i<n;i++)
	{
		for(j=i;j<n;j++)
		{
			sum=A[i][j];
			for(k=i-1;k>=0;k--)
			{
				sum=sum-A[i][k]*A[j][k];
			}
			if(i==j)
			{
				p[i]=sqrt(sum);
			}else{
				A[j][i]=sum/p[i];
			}
		}
	}
	
	for(i=0;i<n;i++)
	{
		A[i][i]=p[i];
		for(j=i+1;j<n;j++)
		{
			A[i][j]=0;
		}
	}
	//free(p);	
}


void MatCholDc(double **A,double *p,int n)
{
	int i,j,k;
	double sum;
	
	for(i=0;i<n;i++)
	{
		for(j=i;j<n;j++)
		{
            sum=A[i][j];
            for(k=i-1;k>=0;k--)
            {
                sum=sum-A[i][k]*A[j][k];
            }
            if(i==j)
            {
                p[i]=sqrt(sum);
            }else{
                A[j][i]=sum/p[i];
            }
        }
	}
}

void MatCholSl(double **A, double *p,double **b, double **x, int n)
{
    int i,k,m;
    double sum;
    
    for(i=0;i<n;i++)
    {
        sum=b[i][0];
        for(k=i-1;k>=0;k--)
        {
            sum=sum-A[i][k]*x[k][0];
        }
        x[i][0]=sum/p[i];
    }
    m=n-1;
    for(i=m;i>=0;i--)
    {
        sum=x[i][0];
        for(k=i+1;k<n;k++)
        {
            sum=sum-A[k][i]*x[k][0];
        }
        x[i][0]=sum/p[i];
    }
}



/*======================================================
========================================================
 Funcion para calcular QR
========================================================
======================================================*/
void MatQRdcmp(double **A, int n, double *c, double *d, int *sing)
{
	int i,j,k;
	double scale,sigma,sum,tau;
	
	*sing=0;
	for(k=0;k<n-1;k++)
	{
		scale=0;
		for(i=k;i<n;i++)
		{
			scale=FMAX(scale,fabs(A[i][k]));
		}
		if(scale==0)
		{
			*sing=1;
			c[k]=d[k]=0;
		}else{
			for(i=k;i<n;i++)
			{
				A[i][k]/=scale;
			}
			sum=0;
			for(i=k;i<n;i++)
			{
				sum+=SQR(A[i][k]);
			}
			sigma=SIGN(sqrt(sum),A[k][k]);
			A[k][k]+=sigma;
			c[k]=sigma*A[k][k];
			d[k]=-scale*sigma;
			for(j=k+1;j<n;j++)
			{
				sum=0;
				for(i=k;i<n;i++)
				{
					sum+=A[i][k]*A[i][j];
				}
				tau=sum/c[k];
				for(i=k;i<n;i++)
				{
					A[i][j]-=tau*A[i][k];
				}
			}
		}
	}
	d[n-1]=A[n-1][n-1];
	if(d[n-1]==0)
		*sing=1;
	
	for(i=0;i<n;i++)
		A[i][i]=d[i];
	
}
