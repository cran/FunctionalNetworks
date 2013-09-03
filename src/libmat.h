//
//  libmat.h
//  GibbsCS
//
//  Created by Alejandro Quiroz on 1/9/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/PrtUtil.h>


int *arrayInt1(int row);
double *array1(int row);
int **arrayInt2(int row, int col);
double **array2(int row, int col);
int *array1Intsrce(int *srce, int row);
int **array2Intsrce(int *srce, int row, int col);
double *array1srce(double *srce, int row);
double **array2srce(double *srce, int row, int col);
double prodArray1(double *v1, double *v2, int m);
void array2toArray1(double **a, double *v, int nrow, int cols);
double *array2toVec(double **a, double *v ,int nrow, int cols);


/***************************************************************/

int *arrayInt1(int row)
{
	int *mat;
	
	mat=(int*)R_alloc(row,sizeof(int));
	return mat;
}


int *arrayInt1M(int row)
{
	int *mat;
	
	mat=(int*)malloc(row*sizeof(int));
    //if (mat==NULL) {
	//	printf("Falla al asignar memoria a los renglones del arreglo\n");
	//	exit(0);
	//}
	return mat;
}


double *array1(int row)
{
	double *mat;
	
	mat=(double*)R_alloc(row,sizeof(double));
	return mat;
}


double *array1M(int row)
{
	double *mat;
	
	mat=(double*)malloc(row*sizeof(double));
	//if (mat==NULL) {
	//	printf("Falla al asignar memoria a los renglones del arreglo\n");
	//	exit(0);
	//}
	return mat;
}

int **arrayInt2(int row, int col)
{
	int *ptr1, **mat;
	int i;
	
	ptr1=(int*)R_alloc(row*col,sizeof(int));
	mat=(int**)R_alloc(row,sizeof(int*));
	for (i=0; i<row; i++)
		mat[i]=ptr1+(i*col);
	return mat;
}

double **array2(int row, int col)
{
	double *ptr1, **mat;
	int i;
	
	ptr1=(double*)R_alloc(row*col,sizeof(double));
	mat=(double**)R_alloc(row,sizeof(double*));
	for (i=0; i<row; i++)
		mat[i]=ptr1+(i*col);
	return mat;
}

double **array2M(int row, int col)
{
	double *ptr1, **mat;
	int i;
	
	ptr1=(double*)malloc(row*col*sizeof(double));
	//if (ptr1==NULL) {
	//	printf("Falla al asignar memoria a la matriz\n");
	//	exit(0);
	//}
	mat=(double**)malloc(row*sizeof(double*));
	//if (mat==NULL) {
	//	printf("Falla al asignar memoria a los renglones de la matriz\n");
	//	exit(0);
	//}
	for (i=0; i<row; i++)
		mat[i]=ptr1+(i*col);
	return mat;
}

void freeArray(double **array)
{
	free(array[0]);
	free(array);
}


int *array1Intsrce(int *srce, int row)
{
	int *ptr1;
	int i,k;
	ptr1=(int*)R_alloc(row,sizeof(int));
	
	k=0;
    for(i=0; i<row; i++)
    {
        ptr1[i]=srce[k];
        k++;
    }
	return ptr1;
}

int **array2Intsrce(int *srce, int row, int col)
{
	int *ptr1, **mat;
	int i, j, k;
	
	ptr1=(int*)R_alloc(row*col,sizeof(int));
	mat=(int**)R_alloc(row,sizeof(int*));
	for (i=0; i<row; i++)
		mat[i]=ptr1+(i*col);
	
	k=0;
	for(j=0; j<col; j++)
		for(i=0; i<row; i++){
			mat[i][j]=srce[k];
			k++;
		}
	return mat;
}

double *array1srce(double *srce, int row)
{
	double *ptr1;
	int i,k;
	ptr1=(double*)R_alloc(row,sizeof(double));
	
	k=0;
    for(i=0; i<row; i++)
    {
        ptr1[i]=srce[k];
        k++;
    }
	return ptr1;
}

double **array2srce(double *srce, int row, int col)
{
	double *ptr1, **mat;
	int i, j, k;
	
	ptr1=(double*)R_alloc(row*col,sizeof(double));
	mat=(double**)R_alloc(row,sizeof(double*));
	for (i=0; i<row; i++)
		mat[i]=ptr1+(i*col);
	
	k=0;
	for(j=0; j<col; j++)
		for(i=0; i<row; i++){
			mat[i][j]=srce[k];
			k++;
		}
	return mat;
}


double prodArray1(double *v1, double *v2, int m)
{
	double resul;
	int i;
	
	resul=0;
	for (i=0; i<m; i++)
		resul+=(v2[i]*v1[i]);
	return resul;
}

void array2toArray1(double **a, double *v, int nrow, int cols)
{
	int j;
	
	for(j=0; j<cols; j++)
		v[j]=a[nrow][j];
}

double *array2toVec(double **a, double *v ,int nrow, int cols)
{
	int i,j;
	int contador;
	
	contador=0;
	for(i=0;i<nrow;i++)
	{
    	for(j=0; j<cols; j++)
		{
			v[contador]=a[i][j];
			contador=contador+1;
		}
	}
	return v;
}
