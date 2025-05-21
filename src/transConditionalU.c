/*
 *  transConditionalU
 *
 *  Created by Yasaman Tahernezhad on June 29, 2024.
 *
 */
#include <R_ext/Lapack.h>
#include <R_ext/Applic.h>
#include <R_ext/PrtUtil.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <R.h>
#include <Rmath.h>
#include "utility.h"
#include "transConditionalU.h"



/*********************************************************************
 *
 * 
 *
 *********************************************************************/
 
void transConditionalU(double* RsummaryBetas,int* ldJ, int* dims, int* RIndexMatrix, int* IndJ, double* Rldvec, double* RSDvec,  double* RDensityU, double* RJointBmatrix,
   double* RAllAn_i, double* Nvec, double* sqrt_det_Zcov,double *sqrt_det_Ai, double *sqrt_inv_det_An_i, double* RA0, double* RA0_inv, double* Res_alphaMatrix, int* RSharedPattern, int* numI, double* weight, double* RZcov_Inv, double* Rb0Matrix, double* Resjointbetas, double* RpostEb1, double* RpostEb2, double*numvec)
{
    int nrow_ldvec = dims[0], P = dims[1];
    int nrow_JointBmatrix = dims[2], ncol_JointBmatrix = dims[3];
    int nrow_IndexMatrix = dims[4], ncol_IndexMatrix = dims[5];
    int **IndexMatrix;
    double **JointBmatrix, **summaryBetas, **DensityU;
    double  **AllAn_i, **A0, **ldvec, **A0_inv;
    int Q = dims[10], numAlpha = dims[11], num_element = dims[12];
    int nrow_es_alphaMatrix = dims[13], ncol_es_alphaMatrix = dims[14];
    int nrow_SharedPattern = dims[17],ncol_SharedPattern = dims[18];
    int nrow_DensityU = dims[19], ncol_DensityU = dims[20];
    int ** SharedPattern;
    double **SDvec, **es_alphaMatrix, **Zcov_Inv, **b0Matrix, **esjointbetas, **postEb1, **postEb2;
    //*Rprintf("check=%e\n",0.001);
    //Rprintf("check0=%d\n",0);
    
    /* reorganize vector into matrix */
    reorg(RsummaryBetas, &summaryBetas, P, Q);
    reorg_int(RIndexMatrix, &IndexMatrix, nrow_IndexMatrix, ncol_IndexMatrix);
    reorg(RSDvec, &SDvec, P, Q);
    reorg(RDensityU, &DensityU, nrow_DensityU, ncol_DensityU);
    reorg(RpostEb1, &postEb1, nrow_DensityU, ncol_DensityU);
    reorg(RpostEb2, &postEb2, nrow_DensityU, ncol_DensityU);
    
    reorg(RJointBmatrix, &JointBmatrix, nrow_JointBmatrix, ncol_JointBmatrix);
    reorg(RAllAn_i, &AllAn_i, numAlpha, num_element);
    reorg(RA0, &A0, Q, Q);
    reorg(Rldvec, &ldvec, nrow_ldvec, Q);
    reorg(RZcov_Inv, &Zcov_Inv, Q, Q);
    reorg(Rb0Matrix, &b0Matrix, P, Q);
    reorg(Resjointbetas, &esjointbetas, P, Q);
    reorg(RA0_inv, &A0_inv, Q, Q);



    reorg(Res_alphaMatrix, &es_alphaMatrix, nrow_es_alphaMatrix, ncol_es_alphaMatrix);
    reorg_int(RSharedPattern, &SharedPattern, nrow_SharedPattern, ncol_SharedPattern);
    
   
    //Rprintf("check0=%d\n",2);
    transConditionalUC(summaryBetas, ldJ, dims, IndexMatrix, IndJ, ldvec, SDvec,   DensityU, JointBmatrix, AllAn_i, Nvec, sqrt_det_Zcov, sqrt_det_Ai,  sqrt_inv_det_An_i, A0, A0_inv, es_alphaMatrix, SharedPattern, numI, weight, Zcov_Inv, b0Matrix, esjointbetas, postEb1, postEb2, numvec);
} 


void transConditionalUC(double** summaryBetas,int* ldJ, int* dims, int** IndexMatrix, int* IndJ, double** ldvec, double** SDvec,  double** DensityU, double** JointBmatrix,
   double** AllAn_i, double* Nvec, double* sqrt_det_Zcov,double *sqrt_det_Ai, double *sqrt_inv_det_An_i, double** A0, double** A0_inv, double** es_alphaMatrix, int** SharedPattern, int* numI, double* weight, double** Zcov_Inv, double** b0Matrix, double** esjointbetas, double** postEb1, double** postEb2, double* numvec)
{
  
    int P = dims[1];
    int t1, i, j, k, j1, j2, k1;
    int NumSNP = dims[4], NumInd = dims[6];
    int Zscale = dims[7];

    int Q = dims[10], pos, aIndex, numAlpha = dims[11], num_element = dims[12];
    int nrow_JointBmatrix = dims[2], ncol_JointBmatrix = dims[3];
    int nrow_IndexMatrix = dims[4], ncol_IndexMatrix = dims[5];
    
    int nrow_SharedPattern = dims[17],ncol_SharedPattern = dims[18];

    
    double tmp0, *btmp2, bj_bar, alphai;
    double Ipart;
    btmp2 = (double *)calloc(Q, sizeof(double));

    
    int i3,j3, ncol_Matrix1 = Q, nrow_Matrix1 = Q;
    
    double *b0, *usedb0, valbab, *phiInvbetam, valbab0, *phiInvbetam0, pival = pow((2*M_PI),(-1.5));
    b0 = (double *)calloc(Q, sizeof(double));
    usedb0 = (double *)calloc(Q, sizeof(double));
    phiInvbetam = (double *)calloc(Q, sizeof(double));
    phiInvbetam0 = (double *)calloc(Q, sizeof(double));
    int k2;

    
    
    
    GetRNGstate();
   

    for(t1=0; t1<nrow_JointBmatrix; t1++){
      //Rprintf("t1=%d\n",t1);
      if(NumInd!=0){
        for(j1=0; j1<NumInd; j1++){
          j = IndJ[j1];
          for(k1=0; k1<Q; k1++){
            btmp2[k1] = summaryBetas[j][k1];
          }
          for(k1=0; k1<Q; k1++){
            tmp0 = 0.0;
            for(k2=0; k2<Q; k2++){
              tmp0 = tmp0 + A0[k1][k2]*btmp2[k2];
            }
            b0[k1] = tmp0;
            if(t1==0){
              b0Matrix[j][k1] = tmp0;
              esjointbetas[j][k1] =  btmp2[k1];
            }
          }
          for(k1=0; k1<Q; k1++){
            tmp0 = 0.0;
            for(k2=0; k2<Q; k2++){
              tmp0 = tmp0 + b0[k2]*A0_inv[k2][k1];
            }
            phiInvbetam0[k1] = tmp0;
          }
          valbab0 = 0.0;
          for(k1=0; k1<Q; k1++){
            valbab0 = valbab0 + b0[k1]*phiInvbetam0[k1];
          }
          if(valbab<0){
            error("valbab0 negative");
          }
          if(t1==0 && j1==0){
            //Rprintf("j=%d\n",j);
            //Rprintf("valbab0=%e\n",valbab0);
          }
          for(aIndex=0; aIndex<numAlpha; aIndex++){
            alphai = es_alphaMatrix[t1][aIndex];
            if(numI[aIndex]!=0){
              nrow_Matrix1 = numI[aIndex];
              ncol_Matrix1 = numI[aIndex];
              
              k = 0;
              for(i3=0 ; i3<ncol_SharedPattern; i3++){
                if(SharedPattern[aIndex][i3]==1){
                  usedb0[k] = b0[i3];
                  k = k + 1;
                }
              }
              for(i3=0; i3<nrow_Matrix1; i3++){
                tmp0 = 0.0;
                for(j3=0; j3<nrow_Matrix1; j3++){
                  tmp0 = tmp0 + AllAn_i[aIndex][(i3*nrow_Matrix1)+j3]*usedb0[j3];
                }
                phiInvbetam[i3] = tmp0;
              }
              valbab = 0.0;
              for(i3=0; i3<nrow_Matrix1; i3++){
                valbab = valbab + usedb0[i3]*phiInvbetam[i3];
              }
              if(valbab<0){
                error("valbab negative");
              }
              Ipart = alphai*pival*numvec[0]*sqrt_det_Ai[aIndex]*sqrt_inv_det_An_i[aIndex]*exp((-valbab0/2))*exp((valbab/2));
   
            }else{
              Ipart = alphai*pival*numvec[0]*exp((-valbab0/2));
            }
            if(t1==0 && j1==0){
              //Rprintf("aIndex=%d\n",aIndex);
              //Rprintf("valbab=%e\n",valbab);
              //Rprintf("Ipart=%e\n",Ipart);

            }
            DensityU[t1][((j*numAlpha) + aIndex)] = Ipart;
            if(SharedPattern[aIndex][0]==1){
              postEb1[t1][((j*numAlpha) + aIndex)] = phiInvbetam[0];
            }
            if(SharedPattern[aIndex][1]==1){
              if(numI[aIndex]==1){
                postEb2[t1][((j*numAlpha) + aIndex)] = phiInvbetam[0];
              }
              if(numI[aIndex]==2){
                postEb2[t1][((j*numAlpha) + aIndex)] = phiInvbetam[1];
              }
            }
          }
        }
      }
      for(j1=0; j1<NumSNP; j1++){
        //Rprintf("j1=%d\n",j1);
        j = IndexMatrix[j1][0];
        for(k1=0; k1<Q; k1++){
          tmp0 = 0.0;
          for(j2=IndexMatrix[j1][1]; j2<(IndexMatrix[j1][2]+1); j2++){
            tmp0 = tmp0 + ldvec[j2][k1]*JointBmatrix[t1][(ldJ[j2] + k1*P)];
          }
          btmp2[k1] = (summaryBetas[j][k1] - tmp0);
        }
        for(k1=0; k1<Q; k1++){
          tmp0 = 0.0;
          for(k2=0; k2<Q; k2++){
            tmp0 = tmp0 + A0[k1][k2]*btmp2[k2];
            //Rprintf("A0[k1][k2]=%e\n",A0[k1][k2]);
          }
          b0[k1] = tmp0;
          if(t1==0){
            b0Matrix[j][k1] = tmp0;
            esjointbetas[j][k1] =  btmp2[k1];
          }
          //Rprintf("b0[k1]=%e\n",b0[k1]);
          //Rprintf("btmp2[k1]=%e\n",btmp2[k1]);
        }
        for(k1=0; k1<Q; k1++){
          tmp0 = 0.0;
          for(k2=0; k2<Q; k2++){
            tmp0 = tmp0 + b0[k2]*A0_inv[k2][k1];
          }
          phiInvbetam0[k1] = tmp0;
        }
        valbab0 = 0.0;
        for(k1=0; k1<Q; k1++){
          valbab0 = valbab0 + b0[k1]*phiInvbetam0[k1];
        }
        if(valbab<0){
          error("valbab0 negative");
        }
        if(t1==0 && j1==0){
          //Rprintf("j=%d\n",j);
          //Rprintf("valbab0=%e\n",valbab0);
        }
        for(aIndex=0; aIndex<numAlpha; aIndex++){
          //Rprintf("aIndex=%d\n",aIndex);
          alphai = es_alphaMatrix[t1][aIndex];
      
          
          if(numI[aIndex]!=0){
            /* compute A */
            nrow_Matrix1 = numI[aIndex];
            ncol_Matrix1 = numI[aIndex];
              
            k = 0;
            //Rprintf("ncol_SharedPattern=%d\n",ncol_SharedPattern);
            for(i3=0 ; i3<ncol_SharedPattern; i3++){
              //Rprintf("SharedPattern[aIndex][i3]=%d\n",SharedPattern[aIndex][i3]);
              if(SharedPattern[aIndex][i3]==1){
                usedb0[k] = b0[i3];
                k = k + 1;
              }
            }
        
            for(i3=0; i3<nrow_Matrix1; i3++){
              tmp0 = 0.0;
              for(j3=0; j3<nrow_Matrix1; j3++){
                tmp0 = tmp0 + AllAn_i[aIndex][(i3*nrow_Matrix1)+j3]*usedb0[j3];
              }
              phiInvbetam[i3] = tmp0;
            }
            valbab = 0.0;
            for(i3=0; i3<nrow_Matrix1; i3++){
                valbab = valbab + usedb0[i3]*phiInvbetam[i3];
            }
            if(valbab<0){
              error("valbab negative");
            }
            
            //Rprintf("check=%d\n",3);
            Ipart = alphai*pival*numvec[0]*sqrt_det_Ai[aIndex]*sqrt_inv_det_An_i[aIndex]*exp((-valbab0/2))*exp((valbab/2));
   
            /*if(t1==89 && j==17){
            Rprintf("j=%d\n",j);
            Rprintf("aIndex=%d\n",aIndex);
            Rprintf("Ipart=%e\n",Ipart);
            Rprintf("valbab0=%e\n",valbab0);
            Rprintf("valbab=%e\n",valbab);
            Rprintf("alphai=%e\n",alphai);
            }*/
            
          }else{
              Ipart = alphai*pival*numvec[0]*exp((-valbab0/2));
              
          }
          if(t1==0 && j1==0){
            //Rprintf("aIndex=%d\n",aIndex);
            //Rprintf("valbab=%e\n",valbab);
            //Rprintf("Ipart=%e\n",Ipart);

          }
          
          DensityU[t1][((j*numAlpha) + aIndex)] = Ipart;
          if(SharedPattern[aIndex][0]==1){
            postEb1[t1][((j*numAlpha) + aIndex)] = phiInvbetam[0];
          }
          if(SharedPattern[aIndex][1]==1){
            if(numI[aIndex]==1){
              postEb2[t1][((j*numAlpha) + aIndex)] = phiInvbetam[0];
            }
            if(numI[aIndex]==2){
              postEb2[t1][((j*numAlpha) + aIndex)] = phiInvbetam[1];
            }
          }
        }
      }
    }
    //Rprintf("pival=%e\n",pival);
    //Rprintf("numvec[0]=%e\n",numvec[0]);

    PutRNGstate();
 }

 
