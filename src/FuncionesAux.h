//
//  FuncionesAux.h
//  Red2ndAttmpt
//
//  Created by Alejandro Quiroz on 2/15/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include "OpMat.h"



//----------------------------------------------------------------
//------------ Tabla: Indices e indices de interes 0 -------------
//------------------------------------------
/*int *TablaIndices(int *vector,int tam_vector,int indice)
{
    int i,cuantos;
    double max;
    
    max=-30000;
    
    for(i=0;i<tam_vector;i++)
    {
        if(vector[i]>)
        
    }
}*/

//-------------------------------------
//----- min identifier INT -------------
//-------------------------------------
int minimo_ind_int(int num_elements, int *a)
{
    int i,indice;
    int minimo;
    minimo=a[0];
    indice=0;
    for(i=1; i<num_elements; i++)
    {
        if (a[i]<minimo)
        {
            minimo=a[i];
            indice=i;
        }
    }
    return(indice);
}

//-------------------------------------
//----- min identifier -------------
//-------------------------------------
int minimo_ind(int num_elements, double *a)
{
    int i,indice;
    double minimo;
    minimo=a[0];
    indice=0;
    for(i=1; i<num_elements; i++)
    {
        if (a[i]<minimo)
        {
            minimo=a[i];
            indice=i;
        }
    }
    return(indice);
}

//-------------------------------------
//----- Maximo identifier -------------
//-------------------------------------
int maximo_ind(int num_elements, double *a)
{
    int i,indice;
    double max;
    max=-320000;
    indice=0;
    for (i=0; i<num_elements; i++)
    {
        if (a[i]>max)
        {
            max=a[i];
            indice=i;
        }
    }
    return(indice);
}

//-------------------------------------
//---------------- Maximo -------------
//-------------------------------------
double maximo(int num_elements, double *a)
{
    int i; 
    double max;
    max=-320000;
    for (i=0; i<num_elements; i++)
    {
        if (a[i]>max)
        {
            max=a[i];
        }
    }
    return(max);
}

//-------------------------------------
//---- Sample from a vector -----------
//-------------------------------------

int *SampleFun(int cuantos, int tamNoCeros, int tamVec, int *vector, int *resultado)
{
    int i,j,auxTam;
    int *xx,*yy;
    double alea;
    
    xx=arrayInt1M(tamNoCeros);
    yy=arrayInt1M(cuantos);
    auxTam=tamNoCeros;
    
    for(i=0;i<tamNoCeros;i++)
    {
        xx[i]=i;
    }
    
    for(i=0;i<cuantos;i++)
    {
        GetRNGstate();
        alea=runif(0,1);
        PutRNGstate();
        
        j=(int)(auxTam*alea);
        yy[i]=xx[j];
        xx[j]=xx[--auxTam];
    }
    
    for(i=0;i<cuantos;i++)
    {
        resultado[i]=vector[yy[i]];
    }
    
    free(xx);
    free(yy);
    return resultado;
}

int *SamplePesoFun(int cuantos, int tamNoCeros, int tamVec, int *vector,double *pesos, int *resultado)
{
    int i,j,auxTam;
    int *xx,*yy;
    double alea,suma;
    
    xx=arrayInt1M(tamNoCeros);
    yy=arrayInt1M(cuantos);
    auxTam=tamNoCeros;
    
    suma=0;
    for(i=0;i<tamNoCeros;i++)
    {
        xx[i]=i;
        suma=suma+pesos[i];
    }
    
    for(i=0;i<cuantos;i++)
    {
        GetRNGstate();
        alea=runif(0,suma);
        PutRNGstate();
       
        for(j=0;j<auxTam;j++)
        {
            if(alea<pesos[j])
            {
                yy[i]=xx[j];
                xx[j]=xx[--auxTam];
                suma=suma-pesos[j];
                pesos[j]=pesos[auxTam];
                break;
            }else{
                alea=alea-pesos[j];
            }
        }
    }
    
    for(i=0;i<cuantos;i++)
    {
        resultado[i]=vector[yy[i]];
    }
    
    free(xx);
    free(yy);
    return resultado;
}



//-----------------------------------------------------------
//---- Da el Id del padre potencial con mas votos -----------
//-----------------------------------------------------------

int *GrpPtnclPdrsIdsMAXRnd(int obj,int nrgrsg,int ncgrsg,int ncgsgr,int **grsg,int **GenRedSrc,int **gsgr,int *resultado)
{
    int i,j,k,l;
    int ncgrsgaux,ncgsgraux,cntsGensObj,cntsGensSrc,aux,aux2,cntsGrpsSrc,max_voto_ind,cuenta;
    double max_voto;
    
    ncgrsgaux=ncgrsg-1;
    ncgsgraux=ncgsgr-1;
    cntsGensObj=grsg[obj][ncgrsgaux];     // Cuantos genes tiene asociadio el GrpObj
    
    cntsGrpsSrc=0;
    for(i=0;i<cntsGensObj;i++)
    {
        aux=grsg[obj][i]-1;               // GenObj asociado al GrpObj 
        cntsGensSrc=GenRedSrc[aux][3];    // Cuantos genes tiene como GenSrc el GenObj asociadio al GrpObj
        for(j=0;j<cntsGensSrc;j++)
        {
            aux2=GenRedSrc[aux][j]-1;     // Obtenemos el GenSrc para cada GenObj asociadio al GrpObj
            cntsGrpsSrc=cntsGrpsSrc+gsgr[aux2][ncgsgraux];   //Guardamos el numero de Grps que tiene asociado ese GenSrc.
        }
    }
    
    for(i=0;i<7;i++)
    {
        resultado[i]=0;
    }
    if(cntsGrpsSrc!=0)
    {
        double votoaux;
        int *dndGrpSrc;
        double *votos;
        
        dndGrpSrc=arrayInt1M(cntsGrpsSrc);    //Aqui Guardaremos todos los Grps Src
        l=0;
        for(i=0;i<cntsGensObj;i++)
        {
            aux=grsg[obj][i]-1;               // GenObj asociado al GrpObj    
            cntsGensSrc=GenRedSrc[aux][3];    // Cuantos genes tiene como GenSrc el GenObj asociadio al GrpObj
            for(j=0;j<cntsGensSrc;j++)
            {
                aux2=GenRedSrc[aux][j]-1;     // Obtenemos el GenSrc para cada GenObj asociadio al GrpObj      
                for(k=0;k<gsgr[aux2][ncgsgraux];k++)
                {
                    dndGrpSrc[l]=gsgr[aux2][k];   // Aqui guardamos el GrpSrc asociado al GenSrc
                    l=l+1;
                }
            }
        }
        
        votos=array1M(nrgrsg);
        for(i=0;i<nrgrsg;i++)
        {
            votos[i]=0;
        }
        for(i=0;i<cntsGrpsSrc;i++)
        {
            aux=dndGrpSrc[i]-1;
            votoaux=1;
            votos[aux]=votos[aux]+votoaux;
        }
        
        cuenta=0;
        for(i=0;i<3;i++)
        {
            max_voto_ind=maximo_ind(nrgrsg,votos);
            max_voto=maximo(nrgrsg,votos);
            if(max_voto!=0)
            {
                resultado[i]=max_voto_ind+1;
                j=i+4;
                resultado[j]=max_voto;
                votos[max_voto_ind]=0;
                cuenta=cuenta+1;
            }else{
                break;
            }
        }
        resultado[3]=cuenta;
        
        free(dndGrpSrc);
        free(votos);
    }
    
    return resultado;
}

int *GrpPtnclPdrsIds(int obj,int nrgrsg,int ncgrsg,int ncgsgr,int **grsg,int **GenRedSrc,int **gsgr,int *resultado)
{
    int i,j,k,l;
    int ncgrsgaux,ncgsgraux,cntsGensObj,cntsGensSrc,aux,aux2,cntsGrpsSrc,max_voto_ind,cuenta;
    double max_voto;
    
    ncgrsgaux=ncgrsg-1;
    ncgsgraux=ncgsgr-1;
    cntsGensObj=grsg[obj][ncgrsgaux];     // Cuantos genes tiene asociadio el GrpObj
    
    cntsGrpsSrc=0;
    for(i=0;i<cntsGensObj;i++)
    {
        aux=grsg[obj][i]-1;               // GenObj asociado al GrpObj 
        cntsGensSrc=GenRedSrc[aux][3];    // Cuantos genes tiene como GenSrc el GenObj asociadio al GrpObj
        for(j=0;j<cntsGensSrc;j++)
        {
            aux2=GenRedSrc[aux][j]-1;     // Obtenemos el GenSrc para cada GenObj asociadio al GrpObj
            cntsGrpsSrc=cntsGrpsSrc+gsgr[aux2][ncgsgraux];   //Guardamos el numero de Grps que tiene asociado ese GenSrc.
        }
    }
    
    for(i=0;i<4;i++)
    {
        resultado[i]=0;
    }
    if(cntsGrpsSrc!=0)
    {
        double votoaux;
        int *dndGrpSrc;
        double *votos;
        
        dndGrpSrc=arrayInt1M(cntsGrpsSrc);    //Aqui Guardaremos todos los Grps Src
        l=0;
        for(i=0;i<cntsGensObj;i++)
        {
            aux=grsg[obj][i]-1;               // GenObj asociado al GrpObj    
            cntsGensSrc=GenRedSrc[aux][3];    // Cuantos genes tiene como GenSrc el GenObj asociadio al GrpObj
            for(j=0;j<cntsGensSrc;j++)
            {
                aux2=GenRedSrc[aux][j]-1;     // Obtenemos el GenSrc para cada GenObj asociadio al GrpObj      
                for(k=0;k<gsgr[aux2][ncgsgraux];k++)
                {
                    dndGrpSrc[l]=gsgr[aux2][k];   // Aqui guardamos el GrpSrc asociado al GenSrc
                    l=l+1;
                }
            }
        }
        
        votos=array1M(nrgrsg);
        for(i=0;i<nrgrsg;i++)
        {
            votos[i]=0;
        }
        for(i=0;i<cntsGrpsSrc;i++)
        {
            aux=dndGrpSrc[i]-1;
            votoaux=1;
            votos[aux]=votos[aux]+votoaux;
        }
        
        cuenta=0;
        for(i=0;i<3;i++)
        {
            max_voto_ind=maximo_ind(nrgrsg,votos);
            max_voto=maximo(nrgrsg,votos);
            if(max_voto!=0)
            {
                resultado[i]=max_voto_ind+1;
                votos[max_voto_ind]=0;
                cuenta=cuenta+1;
            }else{
                break;
            }
        }
        resultado[3]=cuenta;
        
        free(dndGrpSrc);
        free(votos);
    }
     
    return resultado;
}


int *GrpPtnclPdrsIdsALEA(int obj,int nrgrsg,int ncgrsg,int ncgsgr,int **grsg,int **GenRedSrc,int **gsgr,int *resultado)
{
    int i,j,k;
    int ncgrsgaux,ncgsgraux,cntsGensObj,cntsGens,cntsGensSrc,aux,aux2,cntsUncs,aux3;
    int *votos;
    
    int minCnsts,aleatorio,ss;
    double uu;
    
    votos=arrayInt1M(nrgrsg);
    for(i=0;i<nrgrsg;i++)
    {
        votos[i]=0;
    }
    
    ncgrsgaux=ncgrsg-1;
    ncgsgraux=ncgsgr-1;
    cntsGensObj=grsg[obj][ncgrsgaux];     // Cuantos genes tiene asociadio el GrpObj
    cntsGens=0;
    for(i=0;i<cntsGensObj;i++)
    {
        aux=grsg[obj][i]-1;               // GenObj asociado al GrpObj 
        cntsGensSrc=GenRedSrc[aux][3];    // Cuantos genes tiene como GenSrc el GenObj asociadio al GrpObj
        for(j=0;j<cntsGensSrc;j++)
        {
            aux2=GenRedSrc[aux][j]-1;     // Obtenemos el GenSrc para cada GenObj asociadio al GrpObj
            for(k=0;k<gsgr[aux2][ncgsgraux];k++)
            {
                aux3=gsgr[aux2][k]-1;
                votos[aux3]=votos[aux3]+1;
            }
        }
        cntsGens=cntsGens+cntsGensSrc;
    }
    
    cntsUncs=0;
    for(i=0;i<nrgrsg;i++)
    {
        if(votos[i]!=0)
        {
            cntsUncs=cntsUncs+1;
        }
    }
    for(i=0;i<4;i++)
    {
        resultado[i]=0;
    }
    
    if(cntsGens!=0)
    {
        int *GensUncs;
        
        GensUncs=arrayInt1M(cntsUncs);
        j=0;
        for(i=0;i<nrgrsg;i++)
        {
            if(votos[i]!=0)
            {
                GensUncs[j]=i+1;
                j=j+1;
            }
        }
        
        minCnsts=3;
        if(cntsUncs<3)
        {
            minCnsts=cntsUncs;
        }
        
        for(i=0;i<minCnsts;i++)
        {
            ss=0;
            while(ss==0)
            {
                GetRNGstate();
                uu=runif(0,cntsUncs);
                PutRNGstate();
                if(uu==cntsUncs)
                {
                    aleatorio=uu-1;
                }else{
                    aleatorio=uu;
                }
                
                if(GensUncs[aleatorio]!=0)
                {
                    resultado[i]=GensUncs[aleatorio];
                    ss=1;
                    GensUncs[aleatorio]=0;
                }
            }
        }
        resultado[3]=minCnsts;
        free(GensUncs);
    }
    
    free(votos);
    return resultado;
}


int *GrpPtnclPdrsIdsALEAw(int obj,int nrgrsg,int ncgrsg,int ncgsgr,int **grsg,int **GenRedSrc,int **gsgr,int *resultado)
{
    int i,j,k;
    int ncgrsgaux,ncgsgraux,cntsGensObj,cntsGens,cntsGensSrc,aux,aux2,cntsUncs,aux3;
    int *votos;
    
    int minCnsts;
    double uu;
    
    int sumaVotos;
    
    
    votos=arrayInt1M(nrgrsg);
    for(i=0;i<nrgrsg;i++)
    {
        votos[i]=0;
    }
    
    ncgrsgaux=ncgrsg-1;
    ncgsgraux=ncgsgr-1;
    cntsGensObj=grsg[obj][ncgrsgaux];     // Cuantos genes tiene asociadio el GrpObj
    cntsGens=0;
    for(i=0;i<cntsGensObj;i++)
    {
        aux=grsg[obj][i]-1;               // GenObj asociado al GrpObj 
        cntsGensSrc=GenRedSrc[aux][3];    // Cuantos genes tiene como GenSrc el GenObj asociadio al GrpObj
        for(j=0;j<cntsGensSrc;j++)
        {
            aux2=GenRedSrc[aux][j]-1;     // Obtenemos el GenSrc para cada GenObj asociadio al GrpObj
            for(k=0;k<gsgr[aux2][ncgsgraux];k++)
            {
                aux3=gsgr[aux2][k]-1;
                votos[aux3]=votos[aux3]+1;
            }
        }
        cntsGens=cntsGens+cntsGensSrc;
    }
    
    cntsUncs=0;
    sumaVotos=0;
    for(i=0;i<nrgrsg;i++)
    {
        if(votos[i]!=0)
        {
            cntsUncs=cntsUncs+1;
            sumaVotos=sumaVotos+votos[i];
        }
    }
    for(i=0;i<4;i++)
    {
        resultado[i]=0;
    }
    
    if(cntsGens!=0)
    {
        int *GensUncs,*PesosUncs;
        
        GensUncs=arrayInt1M(cntsUncs);
        PesosUncs=arrayInt1M(cntsUncs);
        j=0;
        for(i=0;i<nrgrsg;i++)
        {
            if(votos[i]!=0)
            {
                GensUncs[j]=i+1;
                PesosUncs[j]=votos[i];
                j=j+1;
            }
        }
        
        minCnsts=3;
        if(cntsUncs<3)
        {
            minCnsts=cntsUncs;
        }
        
        for(i=0;i<minCnsts;i++)
        {
            GetRNGstate();
            uu=runif(0,sumaVotos);
            PutRNGstate();
            for(j=0;j<cntsUncs;j++)
            {
                if(uu<PesosUncs[j])
                {
                    resultado[i]=GensUncs[j];
                    sumaVotos=sumaVotos-PesosUncs[j];
                    GensUncs[j]=0;
                    PesosUncs[j]=0;
                    break;
                }else{
                    uu=uu-PesosUncs[j];
                }
            }
        }
        resultado[3]=minCnsts;
        free(GensUncs);
        free(PesosUncs);
    }
    
    free(votos);
    return resultado;
}

int *GrpPtnclPdrsIdsALEAwLinks(int obj,int nrgrsg,int nrgsgr,int ncgrsg,int ncgsgr,int **grsg,int **GenRedSrc,int **gsgr,int *resultado)
{
    int i,j,k;
    int ncgrsgaux,ncgsgraux,cntsGensObj,cntsGens,cntsGensSrc,aux,aux2,cntsUncs,aux3,aux4;
    int *votosGrps,*votosGens;
    
    int minCnsts;
    double uu;
    
    int sumaVotos;
    
    
    votosGrps=arrayInt1M(nrgrsg);
    for(i=0;i<nrgrsg;i++)
    {
        votosGrps[i]=0;
    }
    votosGens=arrayInt1M(nrgsgr);    //<----- aqui falta nrgsgr
    for(i=0;i<nrgsgr;i++)
    {
        votosGens[i]=0;
    }
    
    ncgrsgaux=ncgrsg-1;
    ncgsgraux=ncgsgr-1;
    cntsGensObj=grsg[obj][ncgrsgaux];     // Cuantos genes tiene asociadio el GrpObj
    for(i=0;i<cntsGensObj;i++)
    {
        aux=grsg[obj][i]-1;               // GenObj asociado al GrpObj 
        cntsGensSrc=GenRedSrc[aux][3];    // Cuantos genes tiene como GenSrc el GenObj asociadio al GrpObj
        for(j=0;j<cntsGensSrc;j++)
        {
            aux2=GenRedSrc[aux][j]-1;     // Obtenemos el GenSrc para cada GenObj asociadio al GrpObj
            votosGens[aux2]=votosGens[aux2]+1; // Guardamos,el numero de links que tiene cada GenSrc
        }
    }
    
    cntsGens=0;
    for(i=0;i<cntsGensObj;i++)
    {
        aux=grsg[obj][i]-1;               // GenObj asociado al GrpObj 
        cntsGensSrc=GenRedSrc[aux][3];    // Cuantos genes tiene como GenSrc el GenObj asociadio al GrpObj
        for(j=0;j<cntsGensSrc;j++)
        {
            aux2=GenRedSrc[aux][j]-1;     // Obtenemos el GenSrc para cada GenObj asociadio al GrpObj
            for(k=0;k<gsgr[aux2][ncgsgraux];k++)
            {
                aux3=gsgr[aux2][k]-1;
                aux4=votosGens[aux2];
                votosGrps[aux3]=votosGrps[aux3]+aux4; // Los votos seran proporcionales al numero de outgoing links que tiene en total cada grupo debido a los genes.
            }
        }
        cntsGens=cntsGens+cntsGensSrc;
    }
    
    
    cntsUncs=0;
    sumaVotos=0;
    for(i=0;i<nrgrsg;i++)
    {
        if(votosGrps[i]!=0)
        {
            cntsUncs=cntsUncs+1;
            sumaVotos=sumaVotos+votosGrps[i];
        }
    }
    for(i=0;i<4;i++)
    {
        resultado[i]=0;
    }
    
    if(cntsGens!=0)
    {
        int *GensUncs,*PesosUncs;
        
        GensUncs=arrayInt1M(cntsUncs);
        PesosUncs=arrayInt1M(cntsUncs);
        j=0;
        for(i=0;i<nrgrsg;i++)
        {
            if(votosGrps[i]!=0)
            {
                GensUncs[j]=i+1;
                PesosUncs[j]=votosGrps[i];
                j=j+1;
            }
        }
        
        minCnsts=3;
        if(cntsUncs<3)
        {
            minCnsts=cntsUncs;
        }
        
        for(i=0;i<minCnsts;i++)
        {
            GetRNGstate();
            uu=runif(0,sumaVotos);
            PutRNGstate();
            for(j=0;j<cntsUncs;j++)
            {
                if(uu<PesosUncs[j])
                {
                    resultado[i]=GensUncs[j];
                    sumaVotos=sumaVotos-PesosUncs[j];
                    GensUncs[j]=0;
                    PesosUncs[j]=0;
                    break;
                }else{
                    uu=uu-PesosUncs[j];
                }
            }
        }
        resultado[3]=minCnsts;
        free(GensUncs);
        free(PesosUncs);
    }
    
    free(votosGrps);
    free(votosGens);
    return resultado;
}


int *GenPtnclPdrsIds(int obj,int nrgsgr,int ncgsgr,int ncgrsg,int **gsgr,int **GrpRedSrc,int **grsg,int *resultado)
{
    int i,j,k;
    int cntsGrpsObj,cntsGrps,cntsGrpsSrc,ncgsgraux,ncgrsgaux,aux,aux2,cntsUncs,aux3;
    int *votos;
    
    votos=arrayInt1M(nrgsgr);
    for(i=0;i<nrgsgr;i++)
    {
        votos[i]=0;
    }
    
    ncgsgraux=ncgsgr-1;
    ncgrsgaux=ncgrsg-1;
    
    cntsGrpsObj=gsgr[obj][ncgsgraux];     // Cuantos grupos tiene asociadio el GenObj
    cntsGrps=0;
    for(i=0;i<cntsGrpsObj;i++)
    {
        aux=gsgr[obj][i]-1;
        cntsGrpsSrc=GrpRedSrc[aux][3];
        for(j=0;j<cntsGrpsSrc;j++)
        {
            aux2=GrpRedSrc[aux][j]-1;
            for(k=0;k<grsg[aux2][ncgrsgaux];k++)
            {
                aux3=grsg[aux2][k]-1;
                votos[aux3]=votos[aux3]+1;
             }
        }
        cntsGrps=cntsGrps+cntsGrpsSrc;
    }
    
    cntsUncs=0;
    for(i=0;i<nrgsgr;i++)
    {
        if(votos[i]!=0)
        {
            cntsUncs=cntsUncs+1;
        }
    }
    
    for(i=0;i<4;i++)
    {
        resultado[i]=0;
    }
    
    if(cntsGrps!=0)
    {
        int *GensUncs,*mstrGen,*mstrGen2;
        GensUncs=arrayInt1M(cntsUncs);
        j=0;
        for(i=0;i<nrgsgr;i++)
        {
            if(votos[i]!=0)
            {
                GensUncs[j]=i;
                j=j+1;
            }
        }
        
        mstrGen=arrayInt1M(cntsUncs);
        mstrGen2=arrayInt1M(cntsUncs);
        mstrGen=SampleFun(3,cntsUncs,cntsUncs,GensUncs,mstrGen2);
        for(i=0;i<3;i++)
        {
            resultado[i]=mstrGen[i]+1;
        }
        resultado[3]=3;
        
        free(GensUncs);
        free(mstrGen);
    }
    
    free(votos);
    
    return resultado;
}
    
int *GenPtnclPdrsIdsMAX(int obj,int nrgsgr,int ncgsgr,int ncgrsg,int **gsgr,int **GrpRedSrc,int **grsg,int *resultado)
{
    int i,j,k,l;
    int cntsGrpsObj,cntsGrpsSrc,cntsGenSrc,ncgsgraux,ncgrsgaux,aux,aux2,cuenta,max_voto_ind;
    double max_voto;
    
    ncgsgraux=ncgsgr-1;
    ncgrsgaux=ncgrsg-1;
    
    cntsGrpsObj=gsgr[obj][ncgsgraux];     // Cuantos grupos tiene asociadio el GenObj
    cntsGenSrc=0;
    for(i=0;i<cntsGrpsObj;i++)
    {
        aux=gsgr[obj][i]-1;
        cntsGrpsSrc=GrpRedSrc[aux][3];
        for(j=0;j<cntsGrpsSrc;j++)
        {
            aux2=GrpRedSrc[aux][j]-1;
            cntsGenSrc=cntsGenSrc+grsg[aux2][ncgrsgaux];
        }
    }
    
    for(i=0;i<4;i++)
    {
        resultado[i]=0;
    }
    
    if(cntsGenSrc!=0)
    {
        int *dndGenSrc;
        double *votos;
        
        dndGenSrc=arrayInt1M(cntsGenSrc);
        l=0;
        for(i=0;i<cntsGrpsObj;i++)
        {
            aux=gsgr[obj][i]-1;
            cntsGrpsSrc=GrpRedSrc[aux][3];
            for(j=0;j<cntsGrpsSrc;j++)
            {
                aux2=GrpRedSrc[aux][j]-1;
                for(k=0;k<grsg[aux2][ncgrsgaux];k++)
                {
                    dndGenSrc[l]=grsg[aux2][k];
                    l=l+1;
                }
            }
        }
        
        votos=array1M(nrgsgr);
        for(i=0;i<nrgsgr;i++)
        {
            votos[i]=0;
        }
        for(i=0;i<cntsGenSrc;i++)
        {
            aux=dndGenSrc[i]-1;
            votos[aux]=votos[aux]+1;
        }
        
        cuenta=0;
        for(i=0;i<3;i++)
        {
            max_voto_ind=maximo_ind(nrgsgr,votos);
            max_voto=maximo(nrgsgr,votos);
            if(max_voto!=0)
            {
                resultado[i]=max_voto_ind+1;
                votos[max_voto_ind]=0;
                cuenta=cuenta+1;
            }else{
                break;
            }
        }
        resultado[3]=cuenta;
        
        free(dndGenSrc);
        free(votos);
    }
    
    return resultado;
}


int *GenPtnclPdrsIdsMAXRnd(int obj,int nrgsgr,int ncgsgr,int ncgrsg,int **gsgr,int **GrpRedSrc,int **grsg,int *resultado)
{
    int i,j,k,l;
    int cntsGrpsObj,cntsGrpsSrc,cntsGenSrc,ncgsgraux,ncgrsgaux,aux,aux2,cuenta,max_voto_ind;
    double max_voto;
    
    ncgsgraux=ncgsgr-1;
    ncgrsgaux=ncgrsg-1;
    
    cntsGrpsObj=gsgr[obj][ncgsgraux];     // Cuantos grupos tiene asociadio el GenObj
    cntsGenSrc=0;
    for(i=0;i<cntsGrpsObj;i++)
    {
        aux=gsgr[obj][i]-1;
        cntsGrpsSrc=GrpRedSrc[aux][3];
        for(j=0;j<cntsGrpsSrc;j++)
        {
            aux2=GrpRedSrc[aux][j]-1;
            cntsGenSrc=cntsGenSrc+grsg[aux2][ncgrsgaux];
        }
    }
    
    for(i=0;i<7;i++)
    {
        resultado[i]=0;
    }
    
    if(cntsGenSrc!=0)
    {
        int *dndGenSrc;
        double *votos;
        
        dndGenSrc=arrayInt1M(cntsGenSrc);
        l=0;
        for(i=0;i<cntsGrpsObj;i++)
        {
            aux=gsgr[obj][i]-1;
            cntsGrpsSrc=GrpRedSrc[aux][3];
            for(j=0;j<cntsGrpsSrc;j++)
            {
                aux2=GrpRedSrc[aux][j]-1;
                for(k=0;k<grsg[aux2][ncgrsgaux];k++)
                {
                    dndGenSrc[l]=grsg[aux2][k];
                    l=l+1;
                }
            }
        }
        
        votos=array1M(nrgsgr);
        for(i=0;i<nrgsgr;i++)
        {
            votos[i]=0;
        }
        for(i=0;i<cntsGenSrc;i++)
        {
            aux=dndGenSrc[i]-1;
            votos[aux]=votos[aux]+1;
        }
        
        cuenta=0;
        for(i=0;i<3;i++)
        {
            max_voto_ind=maximo_ind(nrgsgr,votos);
            max_voto=maximo(nrgsgr,votos);
            if(max_voto!=0)
            {
                resultado[i]=max_voto_ind+1;
                j=i+4;
                resultado[j]=max_voto;
                votos[max_voto_ind]=0;
                cuenta=cuenta+1;
            }else{
                break;
            }
        }
        resultado[3]=cuenta;
        
        free(dndGenSrc);
        free(votos);
    }
    
    return resultado;
}
//--------------------------------------------------------------------------
//---- Da un vector de 4, con los ids de los potenciales padres y tam ------
//--------------------------------------------------------------------------
int *GrpSrcRedPropFunMAXRnd(int *GrpsPtncls,int *Pesos,int *resultado)
{
    int i,aleatorios;
    int *x,*y,*z;
    int *x_p,*y_p,*z_p;
    double *x_pesos;
    
    x=arrayInt1M(4);
    y=arrayInt1M(4);
    z=arrayInt1M(4);
    
    x_p=arrayInt1M(3);
    x_pesos=array1M(3);
    y_p=arrayInt1M(3);
    z_p=arrayInt1M(3);
    
    for(i=0;i<4;i++)
    {
        x[i]=i;
    }
    
    for(i=0;i<3;i++)
    {
        x_p[i]=GrpsPtncls[i];
        x_pesos[i]=Pesos[i];
    }
    
    y=SampleFun(1,4,4,x,z);
    aleatorios=y[0];
    
    for(i=0;i<4;i++)
    {
        resultado[i]=0;
    }
    
    if((aleatorios!=0) & (GrpsPtncls[0]!=0))
    {
        if(aleatorios>GrpsPtncls[3])
        {
            aleatorios=GrpsPtncls[3];
        }
        y_p=SamplePesoFun(aleatorios,GrpsPtncls[3],3,x_p,x_pesos,z_p);
        for(i=0;i<aleatorios;i++)
        {
            resultado[i]=y_p[i];
        }
        resultado[3]=aleatorios;
    }
    
    free(x);
    free(y);
    
    free(y_p);
    free(x_p);
    free(x_pesos);
    
    return resultado;
}


int *GrpSrcRedPropFun(int *GrpsPtncls,int *resultado)
{
    int i,aleatorios;
    double uu;
    
    GetRNGstate();
    uu=runif(0,4);
    PutRNGstate();
    
    if(uu==4)
    {
        uu=3;
    }else{
        aleatorios=uu;
    }
    
    for(i=0;i<4;i++)
    {
        resultado[i]=0;
    }
    
    if((aleatorios!=0) & (GrpsPtncls[3]!=0))
    {
        if(aleatorios>GrpsPtncls[3])
        {
            aleatorios=GrpsPtncls[3];
        }
        for(i=0;i<aleatorios;i++)
        {
            resultado[i]=GrpsPtncls[i];
        }
        resultado[3]=aleatorios;
    }
    return resultado;
}

int *GenSrcRedPropFun(int *GrpsPtncls,int *resultado)
{
    int i,aleatorios;
    int *x,*y,*z;
    
    x=arrayInt1M(4);
    y=arrayInt1M(4);
    z=arrayInt1M(4);
    for(i=0;i<4;i++)
    {
        x[i]=i;
    }
    y=SampleFun(1,4,4,x,z);
    aleatorios=y[0];
    
    for(i=0;i<4;i++)
    {
        resultado[i]=0;
    }

    if((aleatorios!=0) & (GrpsPtncls[0]!=0))
    {
        for(i=0;i<aleatorios;i++)
        {
            resultado[i]=GrpsPtncls[i];
        }
        resultado[3]=aleatorios;
    }
    
    free(x);
    free(y);
    
    return resultado;
}

int *GenSrcRedPropFunMAX(int *GrpsPtncls,int *resultado)
{
    int i,aleatorios;
    int *x,*y,*z;
    
    x=arrayInt1M(4);
    y=arrayInt1M(4);
    z=arrayInt1M(4);
    for(i=0;i<4;i++)
    {
        x[i]=i;
    }
    y=SampleFun(1,4,4,x,z);
    aleatorios=y[0];
    
    for(i=0;i<4;i++)
    {
        resultado[i]=0;
    }
    
    if((aleatorios!=0) & (GrpsPtncls[3]!=0))
    {
        if(aleatorios>GrpsPtncls[3])
        {
            aleatorios=GrpsPtncls[3];
        }
        for(i=0;i<aleatorios;i++)
        {
            resultado[i]=GrpsPtncls[i];
        }
        resultado[3]=aleatorios;
    }
    
    free(x);
    free(y);
    
    return resultado;
}

int *GenSrcRedPropFunMAXRnd(int *GrpsPtncls,int *Pesos,int *resultado)
{
    int i,aleatorios;
    int *x,*y,*z;
    int *x_p,*y_p,*z_p;
    double *x_pesos;
    
    x=arrayInt1M(4);
    y=arrayInt1M(4);
    z=arrayInt1M(4);
    
    x_p=arrayInt1M(3);
    x_pesos=array1M(3);
    y_p=arrayInt1M(3);
    z_p=arrayInt1M(3);
    
    for(i=0;i<4;i++)
    {
        x[i]=i;
    }
    
    for(i=0;i<3;i++)
    {
        x_p[i]=GrpsPtncls[i];
        x_pesos[i]=Pesos[i];
    }
    
    y=SampleFun(1,4,4,x,z);
    aleatorios=y[0];
    
    for(i=0;i<4;i++)
    {
        resultado[i]=0;
    }
    
    if((aleatorios!=0) & (GrpsPtncls[0]!=0))
    {
        if(aleatorios>GrpsPtncls[3])
        {
            aleatorios=GrpsPtncls[3];
        }
        y_p=SamplePesoFun(aleatorios,GrpsPtncls[3],3,x_p,x_pesos,z_p);
        for(i=0;i<aleatorios;i++)
        {
            resultado[i]=y_p[i];
        }
        resultado[3]=aleatorios;
    }
    
    free(x);
    free(y);
    
    free(y_p);
    free(x_p);
    free(x_pesos);
    
    return resultado;
}

//--------------------------------------------------------------------------
//---- Da un vector de 4, con los ids de los potenciales padres y tam ------
//--------------------------------------------------------------------------
double lmBIC_Chol(double **datos,int *IndT0,int *IndT1,int *auxSrce,int auxObj,int TamT, int TamSrc)
{
    int i,j,k,auxTamX,auxX,auxT0,auxT1;
    double SSR,aux,lenY,lenBeta,loglenY,BIC;
    double **X,**Y,**XT,**XTX,**XTY,**Beta,**YPred,**Residuales; 
    double *pp;
    
    auxTamX=TamSrc+1;
    X=array2M(TamT,auxTamX);
    Y=array2M(TamT,1);
    XT=array2M(auxTamX,TamT);
    XTX=array2M(auxTamX,auxTamX);
    XTY=array2M(auxTamX,1);
    pp=array1M(auxTamX);
    Beta=array2M(auxTamX,1);
    YPred=array2M(TamT,1);
    Residuales=array2M(TamT,1);
    
    for(i=1;i<auxTamX;i++)
    {
        j=i-1;
        auxX=auxSrce[j]-1;
        for(k=0;k<TamT;k++)
        {
            auxT0=IndT0[k]-1;
            X[k][i]=datos[auxX][auxT0];
        }
    }
    
    for(i=0;i<TamT;i++)
    {
        X[i][0]=1;
        auxT1=IndT1[i]-1;
        Y[i][0]=datos[auxObj][auxT1];        
    }
    
    //---------------------------------------------------------
    //---------------------------------------------------------
    MatTrans(X, TamT, auxTamX, XT);
    MatMult(XT, X, auxTamX, TamT, auxTamX, XTX);
    MatMult(XT, Y, auxTamX, TamT, 1, XTY);
    
    MatCholDc(XTX,pp,auxTamX);
    MatCholSl(XTX,pp,XTY,Beta,auxTamX);
    
    MatMult(X, Beta, TamT, auxTamX, 1, YPred);
    MatRes(Y, YPred, TamT, 1, Residuales);
    //---------------------------------------------------------
    //---------------------------------------------------------
    
    SSR=0;
    for(i=0;i<TamT;i++)
    {
        aux=Residuales[i][0]*Residuales[i][0];
        SSR=SSR+aux;
    }
    SSR=SSR/TamT;
    SSR=log(SSR);
    
    lenY=TamT;
    lenBeta=auxTamX;
    loglenY=log(lenY);
    
    BIC=lenY*SSR+lenBeta*loglenY;
    
    freeArray(X);
    freeArray(Y);
    freeArray(XT);
    freeArray(XTX);
    freeArray(XTY);
    free(pp);
    freeArray(Beta);
    freeArray(YPred);
    freeArray(Residuales);
    
    return BIC;
}



double *lmBIC_Chol_Todos(double **datos,int *IndT0,int *IndT1,int *auxSrce,int auxObj,int TamT, double *resultado)
{
    int i,j,k,ii,renglones,columnas,auxX,auxT0,auxT1;
    int *ind_unos,*ind_no_unos,*ind_data,*aux_mult,*len_b;
    double lenY,lenBeta,loglenY,aux;
    double *pp,*SSR;
    double **X,**Y,**XT,**XTX,**XTY,**Beta,**YPred,**Residuales;
    
    renglones=TamT*8;
    columnas=20;
    
    X=array2M(renglones,columnas);
    Y=array2M(renglones,1);
    XT=array2M(columnas,renglones);
    XTX=array2M(columnas,columnas);
    XTY=array2M(columnas,1);
    pp=array1M(columnas);
    Beta=array2M(columnas,1);
    YPred=array2M(renglones,1);
    Residuales=array2M(renglones,1);
    //==============================================================================
    SSR=array1M(8);
    //==============================================================================
    ind_unos=arrayInt1M(8);
    ind_no_unos=arrayInt1M(12);
    ind_data=arrayInt1M(12);
    aux_mult=arrayInt1M(12);
    len_b=arrayInt1M(8);
    
    ind_unos[0]=0;    len_b[0]=1;
    ind_unos[1]=1;    len_b[1]=2;
    ind_unos[2]=3;    len_b[2]=2;
    ind_unos[3]=5;    len_b[3]=2;
    ind_unos[4]=7;    len_b[4]=3;
    ind_unos[5]=10;   len_b[5]=3;
    ind_unos[6]=13;   len_b[6]=3;
    ind_unos[7]=16;   len_b[7]=4;
    
    ind_no_unos[0]=2;    ind_data[0]=0;    aux_mult[0]=1;
    ind_no_unos[1]=4;    ind_data[1]=1;    aux_mult[1]=2;
    ind_no_unos[2]=6;    ind_data[2]=2;    aux_mult[2]=3;
    ind_no_unos[3]=8;    ind_data[3]=0;    aux_mult[3]=4;
    ind_no_unos[4]=9;    ind_data[4]=1;    aux_mult[4]=4;
    ind_no_unos[5]=11;   ind_data[5]=0;    aux_mult[5]=5;
    ind_no_unos[6]=12;   ind_data[6]=2;    aux_mult[6]=5;
    ind_no_unos[7]=14;   ind_data[7]=1;    aux_mult[7]=6;
    ind_no_unos[8]=15;   ind_data[8]=2;    aux_mult[8]=6;
    ind_no_unos[9]=17;   ind_data[9]=0;    aux_mult[9]=7;
    ind_no_unos[10]=18;  ind_data[10]=1;   aux_mult[10]=7;
    ind_no_unos[11]=19;  ind_data[11]=2;   aux_mult[11]=7;
    
    for(i=0;i<renglones;i++)
    {
        for(j=0;j<columnas;j++)
            X[i][j]=0;
    }
    

    k=0;
    for(i=0;i<8;i++)
    {
        for(j=0;j<TamT;j++)
        {
            X[k][ind_unos[i]]=1;
            auxT1=IndT1[j]-1;
            Y[k][0]=datos[auxObj][auxT1];
            k=k+1;
        }
    }
    for(i=0;i<12;i++)
    {
        auxX=auxSrce[ind_data[i]]-1;
        for(j=0;j<TamT;j++)
        {
            k=j+(TamT*aux_mult[i]);
            auxT0=IndT0[j]-1;
            X[k][ind_no_unos[i]]=datos[auxX][auxT0];
        }
    }
    
    MatTrans(X, renglones, columnas, XT);
    MatMult(XT, X, columnas, renglones, columnas, XTX);
    MatMult(XT, Y, columnas, renglones, 1, XTY);
    MatCholDc(XTX,pp,columnas);
    MatCholSl(XTX,pp,XTY,Beta,columnas);
    MatMult(X, Beta, renglones, columnas, 1, YPred);
    MatRes(Y, YPred, renglones, 1, Residuales);
    
    //---------------------------------------------------------
    //---------------------------------------------------------
    k=0;
    for(i=0;i<8;i++)
    {
        SSR[i]=0;
        for(j=0;j<TamT;j++)
        {
            aux=Residuales[k][0]*Residuales[k][0];
            SSR[i]=SSR[i]+aux;
            k=k+1;
        }
        SSR[i]=SSR[i]/TamT;
        SSR[i]=log(SSR[i]);
        
        lenY=TamT;
        lenBeta=len_b[i];
        loglenY=log(lenY);
        
        resultado[i]=lenY*SSR[i]+lenBeta*loglenY;
    }
    
    freeArray(X);
    freeArray(Y);
    freeArray(XT);
    freeArray(XTX);
    freeArray(XTY);
    free(pp);
    freeArray(Beta);
    freeArray(YPred);
    freeArray(Residuales);
    //==============================================================================
    free(SSR);
    //==============================================================================
    free(ind_unos);
    free(ind_no_unos);
    free(ind_data);
    free(aux_mult);
    free(len_b);
    
    return resultado;
}




int **RedProm(int tamRed, int *NoCeros, int **RedIter, int **MatResult)
{
    int i,j,k;
    
    for(i=0;i<tamRed;i++)
    {
        for(j=0;j<NoCeros[i];j++)
        {
            k=RedIter[i][j]-1;
            MatResult[i][k]=MatResult[i][k]+1;
        }
    }    
    return MatResult;
}



int **RedPromPreBIC(int tamRed, int *NoCeros, int *RedIter, int **MatResult)
{
    int i,j,k;
    
    for(i=0;i<tamRed;i++)
    {
        if(NoCeros[i]==1)
        {
            k=RedIter[i]-1;
            MatResult[i][k]=MatResult[i][k]+1;
        }
    }
    return MatResult;
}
