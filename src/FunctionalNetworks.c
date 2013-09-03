//
//  main.c
//  RedPaper2
//
//  Created by Alejandro Quiroz on 6/25/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include "FuncionesAux.h"
#include <time.h>


void Red(int *fcnets,int *nrfcnet,int *gnets,int *nrgnet,int *procagens,int *nrpag,int *ncpag,int *gensaproc,int *nrgap,int *ncgap,double *data,int *nrdata,int *ncdata,double *datafc,int *nrdatafc,int *ncdatafc,int *dondeaffy,int *tamaffy,int *day,int *tamday,int *sim,int *burnin,int *matResG,int *matResFC,double *bicresgen,double *bicres,double *ssrg,double *ssrfc,double *bicPreGen0, double *bicPreGen1, double *bicPreFc0,double *bicPreFc1)
{
    
    //---------------------------
    //--- Variables del Input ---
    //---------------------------
    int *FCNetS,*GNetS,**ProcaGens,**GensaProc;
    int *DondeAffy,*Dia;
    int nrFCNet,nrGNet,nrPaG,ncPaG,nrGaP,ncGaP,nrDatos,ncDatos,nrDatosFC,ncDatosFC,TamAffy,TamDia,nsim,BurnIn;
    double **Datos,**DatosFC;
    
    double *BIC_Pre_Gen_0,**BIC_Pre_Gen_1,*BIC_Pre_Fc_0,**BIC_Pre_Fc_1;
    
    nrFCNet=*nrfcnet;
    FCNetS=array1Intsrce(fcnets,nrFCNet);
    nrGNet=*nrgnet;
    GNetS=array1Intsrce(gnets,nrGNet);
    nrPaG=*nrpag;
    ncPaG=*ncpag;
    ProcaGens=array2Intsrce(procagens,nrPaG,ncPaG);
    nrGaP=*nrgap;
    ncGaP=*ncgap;
    GensaProc=array2Intsrce(gensaproc,nrGaP,ncGaP);
    nrDatos=*nrdata;
    ncDatos=*ncdata;
    Datos=array2srce(data,nrDatos,ncDatos);
    nrDatosFC=*nrdatafc;
    ncDatosFC=*ncdatafc;
    DatosFC=array2srce(datafc,nrDatosFC,ncDatosFC);
    TamAffy=*tamaffy;
    DondeAffy=array1Intsrce(dondeaffy,TamAffy);
    TamDia=*tamday;
    Dia=array1Intsrce(day,TamDia);
    nsim=*sim;
    BurnIn=*burnin;
    
    BIC_Pre_Gen_0=array1srce(bicPreGen0,nrGNet);
    BIC_Pre_Gen_1=array2srce(bicPreGen1,nrGNet,nrGNet);
    
    BIC_Pre_Fc_0=array1srce(bicPreFc0,nrFCNet);
    BIC_Pre_Fc_1=array2srce(bicPreFc1,nrFCNet,nrFCNet);
    
    //=============================================================================
    //=====   Variables del proceso
    //=============================================================================
    int kk,j,i,k,ii,ww,vv,ncGaPaux,ncPaGaux,auxproc,auxproc_i,aux_fc_pred,indice,tamauxvec,PtnclsPredGns,AffyProp,Affy0,AffyObj,comparacion_cual,NoCerosSrc_aux_k,bandera_gen,bb,cc,contador_nsim,contador_nsim_print;
    int *NoCerosSrcProp,*NoCerosSrcFC,*NoCerosSrc0,*tabla_genes,*NumPredGeneAux,*GNetSProp,*tabla_genes_fc,*FCNetSProp,*vec_nsim;
    double suma_bic_gene,suma_bic0_gene,suma_pesos,uu,Bic0,comparacion,BicGNProp,suma_g,suma_3g,suma_4g,TamDiaD,aa,g_abs,nrGNet_d,nrFCNet_d;
    double *tabla_pesos,*BicPropVec,*comparacion_vec,*BicSSRGProp,*BicSSRG0;
    
    int tamaux,aux,aux_k,aux_gene,PtnclsPredFC,FCProp,FC0,FCObj,NoCerosSrcFC_aux,cntsIgls,renglon,columna,contador;
    int *NoCerosSrcFCProp;
    int **MatRedFC,**MatRedGen;
    double BicFCProp,g_affi,g_affi_s;
    double *tabla_pesos_fc,*BicSSRFCProp,*BicSSRFC0;
    
    
    NoCerosSrcProp=arrayInt1(nrGNet);
    tabla_pesos=array1(nrGNet);
    ncPaGaux=ncPaG-1;
    ncGaPaux=ncGaP-1;
    
    NoCerosSrc0=arrayInt1(nrGNet);
    tabla_genes=arrayInt1(nrGNet);
    for(i=0;i<nrGNet;i++)
    {
        tabla_genes[i]=i+1;
        NoCerosSrc0[i]=1;
    }
    
    NoCerosSrcFC=arrayInt1(nrFCNet);
    tabla_genes_fc=arrayInt1(nrFCNet);
    for(i=0;i<nrFCNet;i++)
    {
        tabla_genes_fc[i]=i+1;
        NoCerosSrcFC[i]=1;
    }
    
    NumPredGeneAux=arrayInt1(2);
    BicPropVec=array1(2);
    comparacion_vec=array1(2);
    GNetSProp=arrayInt1(nrGNet);
    BicSSRGProp=array1(nrGNet);
    BicSSRG0=array1(nrGNet);
    TamDiaD=TamDia;
    nrGNet_d=nrGNet;
    nrFCNet_d=nrFCNet;
    
    FCNetSProp=arrayInt1(nrFCNet);
    
    vec_nsim=arrayInt1(100);
    bb=nsim/100;
    cc=(bb+0.5);
    for(i=0;i<99;i++)
        vec_nsim[i]=(cc*(i+1));
    vec_nsim[99]=nsim-1;
    
    contador_nsim=0;
    contador_nsim_print=1;
    
    NoCerosSrcFCProp=arrayInt1(nrFCNet);
    tabla_pesos_fc=array1(nrFCNet);
    
    BicSSRFCProp=array1(nrFCNet);
    BicSSRFC0=array1(nrFCNet);
    MatRedFC=arrayInt2(nrFCNet,nrFCNet);
    MatRedGen=arrayInt2(nrGNet,nrGNet);
    
    for(i=0;i<nrFCNet;i++)
    {
        for(j=0;j<nrFCNet;j++)
            MatRedFC[i][j]=0;
    }
    for(i=0;i<nrGNet;i++)
    {
        for(j=0;j<nrGNet;j++)
            MatRedGen[i][j]=0;
    }
    
    for(kk=0;kk<nsim;kk++)
    {
        suma_bic_gene=0;
        suma_bic0_gene=0;
        suma_3g=0;
        suma_4g=0;
        //======================================
        //--- Parte update de la red de GEN ----
        //======================================
        for(j=0;j<nrGNet;j++)
        {
            //--- Reseteo de la tabla de pesos
            NoCerosSrcProp[j]=0;
            for(i=0;i<nrGNet;i++)
                tabla_pesos[i]=0;      // Resetamos los pesos de la tabla para la seleccion aleatoria
            tamauxvec=GensaProc[j][ncGaPaux];
            suma_pesos=0;
            for(k=0;k<tamauxvec;k++)
            {
                auxproc=GensaProc[j][k];   //   Obtenemos el proceso asociados al gen de interes
                auxproc_i=auxproc-1;           // escalamos para indices de c
                if(NoCerosSrcFC[auxproc_i]==1)
                {
                    aux_fc_pred=FCNetS[auxproc_i]-1;
                    for(ii=0;ii<ProcaGens[aux_fc_pred][ncPaGaux];ii++)
                    {
                        indice=ProcaGens[aux_fc_pred][ii]-1;      // Obtenemos el gen asociado a ese proceso
                        tabla_pesos[indice]=tabla_pesos[indice]+1;  // Registramos ese gen en los pesos
                        suma_pesos=suma_pesos+1;
                    }
                }
            }
            suma_pesos=suma_pesos-tabla_pesos[j];
            tabla_pesos[j]=0;
            //=============================================
            
            GetRNGstate();
            uu=runif(0,suma_pesos);
            PutRNGstate();
            for(ii=0;ii<nrGNet;ii++)
            {
                if(uu<tabla_pesos[ii])
                {
                    PtnclsPredGns=tabla_genes[ii];
                    suma_pesos=suma_pesos-tabla_pesos[ii];
                    tabla_pesos[ii]=0;
                    break;
                }else{
                    uu=uu-tabla_pesos[ii];
                }
            }
            AffyProp=PtnclsPredGns-1;
            if(NoCerosSrc0[j]==1)
                Affy0=GNetS[j]-1;
            AffyObj=j;
            
            //=====================================================
            if(NoCerosSrc0[j]==0)
            {
                Bic0=BIC_Pre_Gen_0[AffyObj];
            }else{
                Bic0=BIC_Pre_Gen_1[AffyObj][Affy0];
            }
            
            NumPredGeneAux[0]=0;
            BicPropVec[0]=BIC_Pre_Gen_0[AffyObj];
            comparacion_vec[0]=BicPropVec[0]-Bic0;
            
            NumPredGeneAux[1]=1;
            BicPropVec[1]=BIC_Pre_Gen_1[AffyObj][AffyProp];
            comparacion_vec[1]=BicPropVec[1]-Bic0;
            
            comparacion_cual=minimo_ind(2,comparacion_vec);
            //=============================================================
            comparacion=comparacion_vec[comparacion_cual];
            
            GNetSProp[j]=0;
            NoCerosSrc_aux_k=NoCerosSrc0[j];
            
            NoCerosSrcProp[j]=NumPredGeneAux[comparacion_cual];
            if(NumPredGeneAux[comparacion_cual]==1)
                GNetSProp[j]=(AffyProp+1);
            BicGNProp=BicPropVec[comparacion_cual];
            
            if(comparacion>0)
            {
                aa=1;
                if(comparacion<30)
                    aa=exp(comparacion)/(exp(comparacion)+1);
                GetRNGstate();
                uu=runif(0,1);
                PutRNGstate;
                if(uu<aa)
                {
                    NoCerosSrcProp[j]=NoCerosSrc_aux_k;
                    GNetSProp[j]=0;
                    if(NoCerosSrcProp[j]==1)
                        GNetSProp[j]=(Affy0+1);
                    BicGNProp=Bic0;
                }
            }
            
            BicSSRGProp[j]=BicGNProp;
            BicSSRG0[j]=Bic0;
            
            suma_bic_gene=suma_bic_gene+BicGNProp;
            suma_bic0_gene=suma_bic0_gene+Bic0;
            suma_3g=suma_3g+exp((BicSSRG0[j]-(NoCerosSrc0[j]+1)*log(TamDiaD))/TamDiaD);
            suma_4g=suma_4g+exp((BicSSRGProp[j]-(NoCerosSrcProp[j]+1)*log(TamDiaD))/TamDiaD);
        }
        comparacion=suma_bic_gene-suma_bic0_gene;
        //=========== Random sampling =================
        bandera_gen=0;
        if(comparacion>0)
        {
            aa=1;
            if(comparacion<30)
                aa=exp(comparacion)/(exp(comparacion)+1);
            GetRNGstate();
            uu=runif(0,1);
            PutRNGstate();
            if(uu<aa)
                bandera_gen=1;
        }
        if(bandera_gen==1)
        {
            bicresgen[kk]=suma_bic0_gene;
            ssrg[kk]=suma_3g;
        }else{
            bicresgen[kk]=suma_bic_gene;
            ssrg[kk]=suma_4g;
            for(i=0;i<nrGNet;i++)
            {
                GNetS[i]=GNetSProp[i];
                NoCerosSrc0[i]=NoCerosSrcProp[i];
            }
        }
        //================================================================
        //================================================================
        //================================================================
        //================================================================
        //================================================================
        //================================================================
        
        suma_bic_gene=0;
        suma_bic0_gene=0;
        for(j=0;j<nrFCNet;j++)
        {
            NoCerosSrcFCProp[j]=0;
            tamaux=ProcaGens[j][ncPaGaux];          // numero de genes del proc
            for(i=0;i<nrFCNet;i++)
                tabla_pesos_fc[i]=0;
            
            suma_pesos=0;
            for(k=0;k<tamaux;k++)
            {
                aux=ProcaGens[j][k];                // gene del proceso objetivo
                aux_k=aux-1;
                if(NoCerosSrc0[aux_k]==1)
                {
                    aux_gene=GNetS[aux_k]-1;      // genes predictores
                    for(ii=0;ii<GensaProc[aux_gene][ncGaPaux];ii++)
                    {
                        auxproc=GensaProc[aux_gene][ii]-1;  //proceso del gene predictor
                        tabla_pesos_fc[auxproc]=tabla_pesos_fc[auxproc]+1;
                        suma_pesos=suma_pesos+1;
                    }
                }
            }
            suma_pesos=suma_pesos-tabla_pesos_fc[j];
            tabla_pesos_fc[j]=0;
            //======================================
            
            GetRNGstate();
            uu=runif(0,suma_pesos);
            PutRNGstate();
            for(ii=0;ii<nrFCNet;ii++)
            {
                if(uu<tabla_pesos_fc[ii])
                {
                    PtnclsPredFC=tabla_genes_fc[ii];
                    suma_pesos=suma_pesos-tabla_pesos_fc[ii];
                    tabla_pesos_fc[ii]=0;
                    break;
                }else{
                    uu=uu-tabla_pesos_fc[ii];
                }
            }
            FCProp=PtnclsPredFC-1;
            if(NoCerosSrcFC[j]==1)
                FC0=FCNetS[j]-1;
            FCObj=j;
            //=====================================================
            if(NoCerosSrcFC[j]==0)
            {
                Bic0=BIC_Pre_Fc_0[FCObj];
            }else{
                Bic0=BIC_Pre_Fc_1[FCObj][FC0];
            }
            
            NumPredGeneAux[0]=0;
            BicPropVec[0]=BIC_Pre_Fc_0[FCObj];
            comparacion_vec[0]=BicPropVec[0]-Bic0;
            
            NumPredGeneAux[1]=1;
            BicPropVec[1]=BIC_Pre_Fc_1[FCObj][FCProp];
            comparacion_vec[1]=BicPropVec[1]-Bic0;
            
            comparacion_cual=minimo_ind(2,comparacion_vec);
            //===========================================================
            comparacion=comparacion_vec[comparacion_cual];
            
            FCNetSProp[j]=0;
            NoCerosSrcFC_aux=NoCerosSrcFC[j];
            
            NoCerosSrcFCProp[j]=NumPredGeneAux[comparacion_cual];
            if(NumPredGeneAux[comparacion_cual]==1)
                FCNetSProp[j]=(FCProp+1);
            BicFCProp=BicPropVec[comparacion_cual];
            
            if(comparacion>0)
            {
                aa=1;
                if(comparacion<30)
                    aa=exp(comparacion)/(exp(comparacion)+1);
                GetRNGstate();
                uu=runif(0,1);
                PutRNGstate();
                if(uu<aa)
                {
                    NoCerosSrcFCProp[j]=NoCerosSrcFC_aux;
                    FCNetSProp[j]=0;
                    if(NoCerosSrcFCProp[j]==1)
                        FCNetSProp[j]=(FC0+1);
                    BicFCProp=Bic0;
                }
            }
            
            BicSSRFCProp[j]=BicFCProp;
            BicSSRFC0[j]=Bic0;
            
            suma_bic_gene=suma_bic_gene+BicFCProp;
            suma_bic0_gene=suma_bic0_gene+Bic0;
            suma_3g=suma_3g+exp((BicSSRFC0[j]-(NoCerosSrcFC[j]+1)*log(TamDiaD))/TamDiaD);
            suma_4g=suma_4g+exp((BicSSRFCProp[j]-(NoCerosSrcFCProp[j]+1)*log(TamDiaD))/TamDiaD);
        }
        comparacion=suma_bic_gene-suma_bic0_gene;
        //=========== Random sampling =================
        bandera_gen=0;
        if(comparacion>0)
        {
            aa=1;
            if(comparacion<30)
                aa=exp(comparacion)/(exp(comparacion)+1);
            GetRNGstate();
            uu=runif(0,1);
            PutRNGstate();
            if(uu<aa)
                bandera_gen=1;
        }
        
        if(bandera_gen==1)
        {
            bicres[kk]=suma_bic0_gene;
            ssrfc[kk]=suma_3g;
        }else{
            bicres[kk]=suma_bic_gene;
            ssrfc[kk]=suma_4g;
            for(i=0;i<nrFCNet;i++)
            {
                FCNetS[i]=FCNetSProp[i];
                NoCerosSrcFC[i]=NoCerosSrcFCProp[i];
            }
        }
        
        if(kk>=BurnIn)
        {
            for(ww=0;ww<nrFCNet;ww++)
            {
                if(NoCerosSrcFC[ww]==1)
                {
                    k=FCNetS[ww]-1;
                    MatRedFC[ww][k]=MatRedFC[ww][k]+1;
                }
            }
            for(ww=0;ww<nrGNet;ww++)
            {
                if(NoCerosSrc0[ww]==1)
                {
                    k=GNetS[ww]-1;
                    MatRedGen[ww][k]=MatRedGen[ww][k]+1;
                }
            }
        }
        
        if(kk==vec_nsim[contador_nsim])
        {
            Rprintf("Percentage of iterations completed: %d \n",contador_nsim_print,contador_nsim_print);
            contador_nsim_print=contador_nsim_print+1;
            contador_nsim=contador_nsim+1;
        }
    }
    //==========================================
    contador=0;
    for(j=0;j<nrFCNet;j++)
    {
        for(i=0;i<nrFCNet;i++)
        {
            matResFC[contador]=MatRedFC[i][j];
            contador=contador+1;
        }
    }
    contador=0;
    for(j=0;j<nrGNet;j++)
    {
        for(i=0;i<nrGNet;i++)
        {
            matResG[contador]=MatRedGen[i][j];
            contador=contador+1;
        }
    }
}





