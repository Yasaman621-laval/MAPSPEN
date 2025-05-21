/*
 *  IterativeSNPC.c
 *
 *  Created by Yasaman Tahernezhad on Feb 06, 2024.
 *
 *  Last updated by Yasaman Tahernezhad on Feb 25, 2025.
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
#include "transPEN.h"



/*********************************************************************
 *
 * 
 *
 *********************************************************************/
 
void transPEN(double* RsummaryBetas,int* All_ldJ, int* dims, int* Numitervec,
     int* RAll_IndexMatrix, int* All_IndJ,
     double* All_ldvec, double* inv_summaryBetas, int*ChrIndexBeta, double*doublevec,
     double* RInit_summaryBetas, double* RSDvec,
     double* RtuningMatrix, double* RBetaMatrix, char **penalty_, int* q_nrowIndexMatrix, int* q_ldJ, int* q_IndJ, int* q_ldvec, int* All_length_wind, int* Dfq1, double* RSSvec)
{
    int NumTuning = dims[0], P = dims[1];
    int nrow_All_IndexMatrix = dims[4], ncol_All_IndexMatrix = dims[5];
    int **All_IndexMatrix;
    int ncol_tuningMatrix = dims[9], Q = dims[10], ncolBetaMatrix = dims[11];
    int nrow_SDvec = dims[16], ncol_SDvec = dims[17];
    int nrow_SSvec = dims[18], ncol_SSvec = dims[19];
    double **BetaMatrix, **tuningMatrix, **summaryBetas, **SDvec, **Init_summaryBetas;
    double **SSvec;
    
    /* reorganize vector into matrix */
    reorg(RsummaryBetas, &summaryBetas, P, Q);
    reorg(RInit_summaryBetas, &Init_summaryBetas, P, Q);
    reorg(RSDvec, &SDvec, nrow_SDvec, ncol_SDvec);
    reorg(RSSvec, &SSvec, nrow_SSvec, ncol_SSvec);
    reorg(RBetaMatrix, &BetaMatrix, NumTuning, ncolBetaMatrix);
    
    reorg(RtuningMatrix, &tuningMatrix, NumTuning, ncol_tuningMatrix);
    reorg_int(RAll_IndexMatrix, &All_IndexMatrix, nrow_All_IndexMatrix, ncol_All_IndexMatrix);
    //Rprintf("check=%d\n",1);
    transPENC(summaryBetas, All_ldJ, dims, Numitervec,
     All_IndexMatrix, All_IndJ,
     All_ldvec, inv_summaryBetas, ChrIndexBeta, doublevec,
     Init_summaryBetas, SDvec, 
     tuningMatrix, BetaMatrix, penalty_, q_nrowIndexMatrix, q_ldJ, q_IndJ, q_ldvec, All_length_wind, Dfq1, SSvec);
} 


void transPENC(double** summaryBetas, int* All_ldJ, int* dims, int* Numitervec,
     int** All_IndexMatrix, int* All_IndJ,
     double* All_ldvec, double* inv_summaryBetas, int*ChrIndexBeta, double*doublevec,
     double** Init_summaryBetas, double** SDvec, 
     double**tuningMatrix, double** BetaMatrix,
     char **penalty_,int* q_nrowIndexMatrix, int* q_ldJ, int* q_IndJ, int* q_ldvec, int* All_length_wind, int* Dfq1, double** SSvec)
{
  
    int NumTuning = dims[0], P = dims[1], NumIter = dims[2], breaking = dims[3];
    int niter, t1, i, j,j1, j2, found, **skipb;
    int NumSNP, NumInd = dims[6];
    int Zscale = dims[7];

    int ncol_tuningMatrix = dims[9], Q = dims[10], k, k1, k2, cindex, warmStart = dims[14], IniBeta = dims[15];
    int OutputForPath = dims[12], dfMax = dims[13], df_q1 = 0.0;
    char *penalty=penalty_[0];
    
    int start_IndJ, start_ldJ, start_ldvec, start_nrowIndexMatrix, ktemp;

    
    double tmp0, tmp1, tmp2, *tmpbetas, threshold = 0.0, bj_bar = 0.0, epsilon = 0.0001;
    double *sumDiffBetas, *sumBetas, *ratio, upperVal = doublevec[0], basew = doublevec[1];
    double lambda0 = 0.0, tau = doublevec[2], alpha = doublevec[3], lambda1 = 0.0;
    double *sumBetas_after, weightN = doublevec[4];
    

    
    /* allocate memory */
    tmpbetas = (double *)calloc(P, sizeof(double));
    sumDiffBetas = (double *)calloc(P, sizeof(double));
    sumBetas = (double *)calloc(P, sizeof(double));
    sumBetas_after = (double *)calloc(P, sizeof(double));
    ratio = (double *)calloc(P, sizeof(double));
    
    GetRNGstate();
    
    double **jointBmatrix, **tempBmatrix;
    
    jointBmatrix = (double **) malloc(P * sizeof(double*));
    tempBmatrix = (double **) malloc(P * sizeof(double*));
    skipb = (int **) malloc(P * sizeof(int*));
  
    jointBmatrix[0] = (double *) calloc(P*Q, sizeof(double));
    if(jointBmatrix[0] == NULL){ error("fail to allocate memory of jointBmatrix"); }
    for(j=0; j<P; j++){  
      jointBmatrix[j] = jointBmatrix[0] + j*Q; 
    }
    
    tempBmatrix[0] = (double *) calloc(P*Q, sizeof(double));
    if(tempBmatrix[0] == NULL){ error("fail to allocate memory of tempBmatrix"); }
    for(j=0; j<P; j++){  
      tempBmatrix[j] = tempBmatrix[0] + j*Q; 
    }
    
    skipb[0] = (int *) calloc(P*Q, sizeof(int));
    if(skipb[0] == NULL){ error("fail to allocate memory of skipb"); }
    for(j=0; j<P; j++){  
      skipb[j] = skipb[0] + j*Q; 
    }

    for(i=0; i<P; i++){  
      sumDiffBetas[i] = 0.0;
      for(j=0; j<Q; j++){  
        jointBmatrix[i][j] = 0.0;
        tempBmatrix[i][j] = 0.0;
        skipb[i][j] = 0;
      } 
    }
    
    //printf("(ncol_tuningMatrix-1)=%d\n",(ncol_tuningMatrix-1));
    //Rprintf("tuningMatrix[0][0]=%f\n",tuningMatrix[0][0]);
    //Rprintf("tuningMatrix[0][1]=%f\n",tuningMatrix[0][1]);
    //Rprintf("(pow(3.0,2.0)=%f\n",pow(3.0,2.0));

    
    for(t1=0; t1<NumTuning; t1++){
      //printf("t1=%d\n",t1);
      Numitervec[t1] = 0;
      
      lambda0 = tuningMatrix[t1][0];
      lambda1 = tuningMatrix[t1][1];
      
      
      if (strcmp(penalty,"LassowII")==0){
        alpha = tuningMatrix[t1][1];
      }
      
      if(warmStart==0){
        for(i=0; i<P; i++){
          for(j=0; j<Q; j++){
            jointBmatrix[i][j] = 0.0;
            tempBmatrix[i][j] = 0.0;
            skipb[i][j] = 0;
          }
          sumDiffBetas[i] = 0.0;
          sumBetas_after[i] = 0.0;
          sumBetas[i] = 0.0;
          ratio[i] = 1.0;
        }
      }
      //Rprintf("check=%d\n",2);
      if(t1 > 0){
        if(warmStart==1 && tuningMatrix[t1][(ncol_tuningMatrix-1)]==0 && tuningMatrix[t1-1][(ncol_tuningMatrix-1)]==0){
          if(Numitervec[t1-1] <= 0){
            continue;
          }
        }
      }
      //Rprintf("t1=%d\n",t1);
      //Rprintf("warmStart=%d\n",warmStart);
      //Rprintf("tuningMatrix[t1][(ncol_tuningMatrix-1)]=%f\n",tuningMatrix[t1][(ncol_tuningMatrix-1)]);
      if(warmStart==1 && tuningMatrix[t1][(ncol_tuningMatrix-1)]==1){
        //Rprintf("warm t1=%d\n",t1);
        for(i=0; i<P; i++){
          for(j=0; j<Q; j++){
            if(IniBeta==1){
              jointBmatrix[i][j] = Init_summaryBetas[i][j];
            }else{
              jointBmatrix[i][j] = 0.0;
            }
            tempBmatrix[i][j] = 0.0;
            skipb[i][j] = 0;
          }
          sumDiffBetas[i] = 0.0;
          sumBetas[i] = 0.0;
          sumBetas_after[i] = 0.0;
          ratio[i] = 1.0;
        }
      }
      //Rprintf("check=%d\n",3);
      cindex = 0;
      for(niter=0; niter<NumIter; niter++){
        //printf("niter=%d\n",niter);
        for(k1=0; k1<Q; k1++){
          //printf("k1=%d\n",k1);
          NumSNP = q_nrowIndexMatrix[k1];
          start_IndJ = 0;
          start_ldJ = 0;
          start_ldvec = 0;
          start_nrowIndexMatrix = 0;
          if(k1 > 0){
            for(ktemp=0; ktemp<k1; ktemp++){
              start_IndJ = start_IndJ + q_IndJ[ktemp];
              start_ldJ = start_ldJ + q_ldJ[ktemp];
              start_ldvec = start_ldvec + q_ldvec[ktemp];
              start_nrowIndexMatrix = start_nrowIndexMatrix + q_nrowIndexMatrix[ktemp];
            }
          }
          NumInd = All_length_wind[k1];
          //Rprintf("check=%d\n",4);
          if(NumInd!=0){
            for(j1=0; j1<NumInd; j1++){
              j = All_IndJ[j1+start_IndJ];
              if (strcmp(penalty,"Lasso")==0){
                threshold = lambda0 + pow(sumDiffBetas[j], alpha);
              }
              if (strcmp(penalty,"LassowII")==0){
                if(sumBetas[j] == 0){
                  threshold = lambda0 + (basew + sumDiffBetas[j])*alpha;
                }else{
                  threshold = lambda0 + sumDiffBetas[j]*alpha;
                }
              }
              if (strcmp(penalty,"LassowIII")==0){
                threshold = lambda0 + lambda1*(sumDiffBetas[j]+tau)/(sumBetas[j]+tau);
              }
              if (strcmp(penalty,"transw")==0){
                threshold = lambda0 + lambda1*(sumDiffBetas[j]+basew)/(sumBetas[j]+tau);
              }
              if (strcmp(penalty,"transwII")==0){
                threshold = lambda0 + lambda1*(ratio[j]+basew);
              }
              if (strcmp(penalty,"transwIII")==0){
                threshold = lambda0 + lambda1*(pow(ratio[j],basew));
              }
              if (strcmp(penalty,"ratiow")==0){
                threshold = lambda0 + (lambda1/(pow((1-ratio[j]),alpha) + tau));
              }
              if (strcmp(penalty,"ratiowS")==0){
                threshold = lambda0 + lambda1*(sumDiffBetas[j] + tau)/(sumBetas_after[j] + tau);
              }
              
              if (strcmp(penalty,"mixLOG")==0){
                  threshold = lambda0 + lambda1*(1/(sumBetas[j] + tau));
              }
              if (strcmp(penalty,"mixLOG_after")==0){
                  threshold = lambda0 + lambda1*(1/(sumBetas_after[j] + tau));
              }
              if (strcmp(penalty,"mixLOGsq")==0){
                  threshold = lambda0 + lambda1*(1/(sumBetas_after[j] + tau));
              }
              
              
              bj_bar = summaryBetas[j][k1];
              if(bj_bar!=0.0){
                if (strcmp(penalty,"ratiowII")==0){
                  threshold = lambda0 + (lambda1/(pow((1-ratio[j]),alpha) + tempBmatrix[j][k1] + tau));
                }
                if (strcmp(penalty,"ratiowIII")==0){
                  threshold = (lambda1/(pow((1-ratio[j]),alpha) + tempBmatrix[j][k1] + tau));
                }
                
                if(Zscale==1){
                  threshold = threshold*SDvec[j][k1];
                }else{
                  threshold = threshold/SSvec[j][k1];
                }
                if(bj_bar > threshold){
                  jointBmatrix[j][k1] = bj_bar - threshold;
                }else if(bj_bar < -threshold){
                  jointBmatrix[j][k1] = bj_bar + threshold;
                }else{
                  jointBmatrix[j][k1] = 0.0;
                }
              }
              if(summaryBetas[j][k1]*jointBmatrix[j][k1]<0){
                Rprintf("summaryBetas[j]=%d\n",j);
                Rprintf("summaryBetas[k1]=%d\n",k1);
                Rprintf("summaryBetas[j][k1]=%e\n",summaryBetas[j][k1]);
                Rprintf("jointBmatrix[j][k1]=%e\n",jointBmatrix[j][k1]);
                error("sign inverse");
              }
            }
          }
          //Rprintf("check=%d\n",5);
          for(j1=0; j1<NumSNP; j1++){
            //Rprintf("j1=%d\n",j1);
            j = All_IndexMatrix[j1+start_nrowIndexMatrix][0];
            if (strcmp(penalty,"Lasso")==0){
              threshold = lambda0 + pow(sumDiffBetas[j], alpha);
            }
            if (strcmp(penalty,"LassowII")==0){
              if(sumBetas[j] == 0){
                threshold = lambda0 + (basew + sumDiffBetas[j])*alpha;
              }else{
                threshold = lambda0 + sumDiffBetas[j]*alpha;
              }
            }
            if (strcmp(penalty,"LassowIII")==0){
              threshold = lambda0 + lambda1*(sumDiffBetas[j]+tau)/(sumBetas[j]+tau);
            }
            if (strcmp(penalty,"transw")==0){
              threshold = lambda0 + lambda1*(sumDiffBetas[j]+basew)/(sumBetas[j]+tau);
            }
            //printf("lambda0=%f\n",lambda0);
            //printf("lambda1=%f\n",lambda1);
            //printf("sumDiffBetas[j]=%f\n",sumDiffBetas[j]);
            //printf("basew=%f\n",basew);
            //printf("tau=%f\n",tau);
            //printf("sumBetas[j]=%f\n",sumBetas[j]);
            //printf("threshold=%f\n",threshold);
            if (strcmp(penalty,"transwII")==0){
              threshold = lambda0 + lambda1*(ratio[j]+basew);
            }
            if (strcmp(penalty,"transwIII")==0){
              threshold = lambda0 + lambda1*(pow(ratio[j],basew));
            }
            if (strcmp(penalty,"ratiow")==0){
              threshold = lambda0 + (lambda1/(pow((1-ratio[j]),alpha) + tau));
            }
            if (strcmp(penalty,"ratiowS")==0){
                threshold = lambda0 + lambda1*(sumDiffBetas[j] + tau)/(sumBetas_after[j] + tau);
            }
            if (strcmp(penalty,"mixLOG")==0){
              threshold = lambda0 + lambda1*(1/(sumBetas[j] + tau));
            }
            if (strcmp(penalty,"mixLOG_after")==0){
              threshold = lambda0 + lambda1*(1/(sumBetas_after[j] + tau));
            }
              
            if(skipb[j][k1]==0){
              if(summaryBetas[j][k1]!=0.0){
                tmp0 = 0.0;
                if(Zscale==1){
                  for(j2=All_IndexMatrix[j1+start_nrowIndexMatrix][1]; j2<(All_IndexMatrix[j1+start_nrowIndexMatrix][2]+1); j2++){
                  tmp0 = tmp0 + All_ldvec[j2+start_ldvec]*jointBmatrix[All_ldJ[j2+start_ldJ]][k1];
                  }
                }else{
                  for(j2=All_IndexMatrix[j1+start_nrowIndexMatrix][1]; j2<(All_IndexMatrix[j1+start_nrowIndexMatrix][2]+1); j2++){
                  tmp0 = tmp0 + sqrt(SSvec[j][k1])*sqrt(SSvec[All_ldJ[j2+start_ldJ]][k1])*All_ldvec[j2+start_ldvec]*jointBmatrix[All_ldJ[j2+start_ldJ]][k1];
                  }
                }
          
                if(Zscale==1){
                  bj_bar = (summaryBetas[j][k1] - tmp0);
                }else{
                  bj_bar = (SSvec[j][k1]*summaryBetas[j][k1] - tmp0)/SSvec[j][k1];
                }
                //printf("threshold=%f\n",threshold);
                //Rprintf("pow(sumDiffBetas[j], alpha)=%f\n",pow(sumDiffBetas[j], alpha));
                //Rprintf("check=%d\n",6);
                if(fabs(bj_bar)>upperVal){
                  if(breaking==1){
                    Numitervec[t1] = -1;
                    break;
                  }else{
                    bj_bar = 0.0;
                    skipb[j][k1] = 1;
                  }
                }
                if (strcmp(penalty,"ratiowII")==0){
                  threshold = lambda0 + (lambda1/(pow((1-ratio[j]),alpha) + tempBmatrix[j][k1] + tau));
                }
                if (strcmp(penalty,"ratiowIII")==0){
                  threshold = (lambda1/(pow((1-ratio[j]),alpha) + tempBmatrix[j][k1] + tau));
                }
                //Rprintf("check=%d\n",7);
                if(Zscale==1){
                  threshold = threshold*SDvec[j][k1];
                }else{
                  threshold = threshold/SSvec[j][k1];
                }
                if(bj_bar > threshold){
                  jointBmatrix[j][k1] = bj_bar - threshold;
                }else if(bj_bar < -threshold){
                  jointBmatrix[j][k1] = bj_bar + threshold;
                }else{
                  jointBmatrix[j][k1] = 0.0;
                }
              }else{
                jointBmatrix[j][k1] = 0.0;
              }  
            }else{
              jointBmatrix[j][k1] = 0.0;
            }
            //Rprintf("check=%d\n",8);
            //printf("jointBmatrix=%f\n",jointBmatrix[j][k1]);
          }
        }
        //Rprintf("check=%d\n",9);
       /**
        * check convergence
        */
              
        found = 0;
        df_q1 = 0;
        for(j=0; j<P; j++){
          for(k1=0; k1<Q; k1++){
            if(fabs(jointBmatrix[j][k1])>upperVal){
              skipb[j][k1] = 1;
              //printf("skipb%d\n",skipb[j][k1]);
            }
            if(k1==0){
              if(fabs(jointBmatrix[j][k1])!=0){
                df_q1 = df_q1 + 1;
              }
              if(df_q1 > dfMax){
                cindex = 1;
                //printf("df_q1=%d\n",df_q1);
                break;
              }
            }
            if(fabs(tempBmatrix[j][k1] - jointBmatrix[j][k1]) > epsilon){
              found = 1;
              break;
            }
          } 
        }
        //printf("df_q1=%d\n",df_q1);
        //printf("dfMax=%d\n",dfMax);
        if(cindex == 1){
          Numitervec[t1] = -2;
          //printf("Numitervec=%d\n",Numitervec[t1]);
          break;
        }
        //printf("found=%d\n",found);
        if(found==0){
          k1=0;
          df_q1 = 0;
          for(i=0; i<Q; i++){
            for(j=0; j<P; j++){
              BetaMatrix[t1][k1] = jointBmatrix[j][i];
              if(i==0){
                if(fabs(jointBmatrix[j][i])!=0){
                  df_q1 = df_q1 + 1;
                }
              }
              k1 = k1 + 1;
            } 
          } 
          Numitervec[t1] = (niter+1);
          Dfq1[t1] = df_q1;
          break;
        }
        for(j=0; j<P; j++){
          tmp1 = 0.0;
          tmp2 = 0.0;
          for(k=0; k<Q; k++){
            tempBmatrix[j][k] = jointBmatrix[j][k];
            if(weightN==0){
              tmp1 = tmp1 + fabs(jointBmatrix[j][k]);
            }else{
              if(k==0){tmp1 = tmp1 + weightN*fabs(jointBmatrix[j][k]);}
              if(k==1){tmp1 = tmp1 + (1-weightN)*fabs(jointBmatrix[j][k]);}
            }
            tmp2 = tmp2 + jointBmatrix[j][k];
          }
          sumBetas[j] = tmp1;
          sumBetas_after[j] = fabs(tmp2);
          
          tmp0 = 0.0;
          for(k1=0; k1<(Q-1); k1++){
            for(k2=(k1+1); k2<Q; k2++){
              tmp0 = tmp0 + fabs(jointBmatrix[j][k1]-jointBmatrix[j][k2]);
            }
          }
          sumDiffBetas[j] = tmp0;
          
          if (strcmp(penalty,"ratiow")==0){
            if(sumBetas[j]==0){
              ratio[j] = 1;
            }else{
              ratio[j] = sumDiffBetas[j]/sumBetas[j];
            }
          }
    
          if (strcmp(penalty,"ratiowS")==0){
            if(sumBetas_after[j]==0){
              ratio[j] = 1;
            }else{
              ratio[j] = sumDiffBetas[j]/sumBetas_after[j];
            }
          }
        }
        //printf("Numitervec1=%d\n",Numitervec[0]);
        //printf("(NumIter-1)=%d\n",(NumIter-1));
        //printf("OutputForPath=%d\n",OutputForPath);

        if(OutputForPath==1 && niter==(NumIter-1)){
          //printf("OutputForPathENTER=%d\n",OutputForPath);
          k1=0;
          for(i=0; i<Q; i++){
            for(j=0; j<P; j++){
              BetaMatrix[t1][k1] = jointBmatrix[j][i];
              k1 = k1 + 1;
            } 
          } 
        }
      }
    } 
    PutRNGstate();
    
    free(jointBmatrix[0]);
    free(tempBmatrix[0]);
    free(skipb[0]);
    free(jointBmatrix);
    free(tempBmatrix);
    free(skipb);
    free(tmpbetas);
    free(sumDiffBetas);

    
 }

 
