/*
 *  libMatInv.h
 *
 *  Created by Alejandro Quiroz on 3/20/11.
 *
 */

#include "libmat.h"  //Las librerias que hacen las matrices y la impresion de los arrelgos.
#define TINY 1.0e-20  // Se definio para LU y la inversa de las matrices.

void LUdcmp(double **A, int n, int *indx)
{
	int i, j, k, imax;
	double big,temp,sum,dum;
	double *vv;
	
	vv=array1(n);
	
    for(i=0;i<n;i++)
	{
		big=0;
		for(j=0;j<n;j++)
		{
			temp=fabs(A[i][j]);
			if(temp>big) 
			{
				big=temp;
			}
		}
		vv[i]=1/big;	   
	}
	
	
	for(j=0;j<n;j++)
	{
		for(i=0;i<j;i++)
		{
			sum=A[i][j];
			for(k=0;k<i;k++)
			{
				sum=sum-A[i][k]*A[k][j];	
			}
			A[i][j]=sum;
		}
		
		big=0;
	    for(i=j;i<n;i++)
	    {
			sum=A[i][j];
			for(k=0;k<j;k++)
			{
				sum=sum-A[i][k]*A[k][j];	
			}
			A[i][j]=sum;
			dum=vv[i]*fabs(sum);
			if(dum>=big)
			{
				big=dum;
				imax=i;
			}
		}
		if(j!=imax)
		{
			for(k=0;k<n;k++)
			{
				dum=A[imax][k];
				A[imax][k]=A[j][k];
				A[j][k]=dum;
			}
			vv[imax]=vv[j];
		}
		
		indx[j]=imax;
		if(A[j][j]==0)
		{
			A[j][j]=TINY;	
		}
		
		if(j!=n)
		{
			dum=1/(A[j][j]);
			for(i=j+1;i<n;i++)
			{
				A[i][j]*=dum;
			}
		}
		
	}
}


void LUbksb(double **A, int n, int *indx, double *b)
{
	int i,ii,ip,j,contador; 
	double sum;
	
    ii=0;
	contador=0;
	
	for(i=0;i<n;i++)
	{
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if(contador!=0){
			for(j=ii;j<=i-1;j++) sum-=A[i][j]*b[j];
		}
		else if (sum!=0) {
			contador=1;
			ii=i;
		}
		b[i]=sum;
	}
	for(i=n-1;i>=0;i--)
	{
		sum=b[i];
		for(j=i+1;j<n;j++)
		{
			sum-=A[i][j]*b[j];
		}
		b[i]=sum/A[i][i];
	}
}
