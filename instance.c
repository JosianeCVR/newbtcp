#include "lp.h"
#include "instance.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <glpk.h>
#include <math.h>

#define RMAIORI '1'
#define RMENORI '2'
#define RIGUAL '4'


//void loadInstance( Instance *inst,const char* filename, char* filesol){
void loadInstance( Instance *inst,const char* filename) {

    int i,j=0,cont=0;
    int *idx1;
    double *coef1;
    int *idx2;
    double *coef2;
    int *equal;
    /*FILE *arq;
    int pos =0;
    char*  linha;
    char dados[500];
    char * words;
    char name[10];
    int valor;*/
    int l=0;
    int col=0;


    LinearProgram *lp = lp_create();
    lp_read(lp,filename);

    /*number cols and rows*/
    inst->nRows = lp_rows(lp);
    inst->nCols = lp_cols(lp);
    inst->varsol = (int*) malloc((inst->nCols+1) * sizeof(int));
    inst->varwrite = (int*) malloc((inst->nCols+1) * sizeof(int));
    equal= (int*) malloc((inst->nRows+1) * sizeof(int));

    // printf("\nNCOLS: %d", inst->nCols);
    // printf("\nNROWS: %d", inst->nRows);

    //arquivo
    /*arq = fopen(filesol, "r");
    if (arq == NULL){
        printf("Erro ao abrir o arquivo \n");
        exit (EXIT_FAILURE);
    }
    linha = fgets(dados, 500, arq);

    while(linha!= NULL){
        if(l>1){
            words  = strtok(linha, " \n");
            pos=0;
            while( words!= NULL ){
                if(pos == 1){
                    strcpy(name,words);
                    col=  lp_col_index(lp,name);
                }
                else if(pos == 2){
                    valor= atoi(words);
                    inst->varwrite[col]= valor;
                }
                    words  = strtok( NULL," \n");
                    pos++;
            }
            i++;
        }
        linha = fgets(dados, 500, arq);
        l++;
    }
    fclose(arq);*/

    /*Type of row*/
    inst->rowType = (char*) malloc((inst->nRows+1) * sizeof(char));
    for( i=0; i<inst->nRows; i++) {
        if(lp_sense(lp,i)=='G') {
            inst->rowType[i+1]='1';
        }
        else if(lp_sense(lp,i)=='L') {
            inst->rowType[i+1]='2';
        }
        else if(lp_sense(lp,i)=='E') {
            inst->rowType[i+1]='4';
            equal[cont]=i+1;
            cont++;
        }
    }
    // printf("\n CONT: %d",cont);
    int num= inst->nRows+cont+1;
    inst->rowType = (char *)realloc (inst->rowType, num*sizeof(char));
    for( i=inst->nRows+1; i<=inst->nRows+cont; i++) {
        inst->rowType[i]='4';
    }
    int nRows= inst->nRows;
    inst->nRows+=cont;
    printf("\nNCOLS: %d", inst->nCols);
    printf("\nNROWS: %d", inst->nRows);
    /*for( i=1; i<=inst->nRows;i++){
        printf("\nTYPE %d: %c",i,inst->rowType[i]);
    }*/

    int nz = lp_nz(lp);
    nz+=cont;
    printf("\n NZ: %d",nz);

    //allocation of memory
    inst->obj = (double*) malloc((inst->nCols+1) * sizeof(double));
    inst->nElRow = (int*) malloc((inst->nRows+1) * sizeof(int));
    idx1 = (int*) malloc((nz+1) * sizeof(int));
    coef1 = (double*) malloc((nz+1) * sizeof(double));
    inst->idxRow = (int**) malloc((inst->nRows+1) * sizeof(int*));
    inst->coefRow =(double**) malloc((inst->nRows+1) * sizeof(double*));
    inst->rhs = (double*) malloc((inst->nRows+1) * sizeof(double));
    inst->colname = (char**) malloc((inst->nCols+1) * sizeof(char*));
    inst->Rowname = (char**) malloc((inst->nRows+1) * sizeof(char*));
    idx2=(int*) malloc((nz+1) * sizeof(int*));
    coef2=(double*) malloc((nz+1) * sizeof(double*));
    inst->nzCol=(int*) malloc((inst->nCols+1) * sizeof(int*));
    inst->idxCol=(int**) malloc((inst->nCols+1) * sizeof(int*));
    inst->coefCol=(double**) malloc((inst->nCols+1) * sizeof(double*));
    inst->newrhs = (double*) malloc((inst->nRows+1) * sizeof(double));
    inst->ncoefrow = (double**) malloc((inst->nRows+1) * sizeof(double*));
    for ( i = 1; i <= inst->nCols; i++ ) {
        inst->colname[i] = (char*) malloc((128) * sizeof(char));
    }
    for ( i = 1; i <= inst->nRows; i++ ) {
        inst->Rowname[i] = (char*) malloc((128) * sizeof(char));
    }
    //column name
    for(i=0; i<inst->nCols; i++) {
        lp_col_name(lp,i,inst->colname[i+1]);
    }
    //row name
    for(i=0; i<nRows; i++) {
        lp_row_name(lp,i,inst->Rowname[i+1]);
    }
    int a=0;
    for(i=nRows+1; i<=inst->nRows; i++) {
        char name[20];
        strcpy(name,"line");
        strcat(name,inst->Rowname[equal[a]]);
        strcpy(inst->Rowname[i],name);
        a++;
    }

    /*for(i=1;i<=inst->nRows;i++){
        printf("\n%d- ROWNAME: %s",i,inst->Rowname[i]);
    }
    for(i=1;i<=inst->nCols;i++){
        printf("\n%d- Colname: %s",i,inst->colname[i]);
    }*/

    const double *res =  lp_obj_coef(lp);
    //coefficients FO
    for( i=0; i<inst->nCols; i++) {
        inst->obj[i+1] = res[i];
    }
    /*for(i=1;i<=inst->nCols;i++){
        printf("\nOBJ: %s %.2f",inst->colname[i],inst->obj[i]);
    }*/

    //numero de elementos de cada linha
    for( i=0; i<nRows; i++) {
        inst->nElRow[i+1]= lp_row(lp,i,idx1,coef1);
    }
    for( i=1; i<=nRows; i++) {
        inst->idxRow[i] = (int*) malloc((inst->nElRow[i]+1) * sizeof(int));
        inst->coefRow[i] = (double*) malloc((inst->nElRow[i]+1) * sizeof(double));
    }
    //coeficientes e indice de cada linha
    for( i=0; i<nRows; i++) {
        inst->nElRow[i+1]= lp_row(lp,i,inst->idxRow[i+1]+1,inst->coefRow[i+1]+1);
    }
    for( i=1; i<=nRows; i++) {
        for(j=1; j<=inst->nElRow[i]; j++) {
            inst->idxRow[i][j]++;
        }
    }
    //novas linhas de restrições que são sinal =
    a=0;
    for( i=nRows+1; i<=inst->nRows; i++) {
        inst->nElRow[i] = inst->nElRow[equal[a]];
        a++;
    }
    for( i=nRows+1; i<=inst->nRows; i++) {
        inst->idxRow[i] = (int*) malloc((inst->nElRow[i]+1) * sizeof(int));
        inst->coefRow[i] = (double*) malloc((inst->nElRow[i]+1) * sizeof(double));
    }
    a=0;
    for( i=nRows+1; i<=inst->nRows; i++) {
        for(j=1; j<=inst->nElRow[i]; j++) {
            inst->idxRow[i][j] =inst->idxRow[equal[a]][j] ;
            inst->coefRow[i][j]=inst->coefRow[equal[a]][j] ;
        }
        a++;
    }
    /* printf("\n");
     for( i=1; i<=inst->nRows;i++){
         printf("%d-(%d) ",i, inst->nElRow[i]);
         for(j=1;j<=inst->nElRow[i];j++){
             printf(" %.2f  %s  ",inst->coefRow[i][j],inst->colname[inst->idxRow[i][j]] );
         }
         printf("\n");
     }*/

    //bound constraint
    for( i=0; i<nRows; i++) {
        inst->rhs[i+1]=lp_rhs(lp,i);
    }
    a=0;
    for( i=nRows+1; i<=inst->nRows; i++) {
        inst->rhs[i] =inst->rhs[equal[a]];
        a++;
    }
    /*for( i=1; i<=inst->nRows;i++){
        printf("\n%d- %.2f",i,inst->rhs[i]);
    }*/
    int *nzcol = (int*) malloc((inst->nCols+1) * sizeof(int));
    //numero de elementos que a variavel i aparece
    for( i=0; i<inst->nCols; i++) {
        inst->nzCol[i+1]= lp_col(lp,i,idx2,coef2);
        nzcol[i+1]=inst->nzCol[i+1];
    }
    for( i=nRows+1; i<=inst->nRows; i++) {
        for(j=1; j<=inst->nElRow[i]; j++) {
            inst->nzCol[inst->idxRow[i][j]]++;
            //  printf("\n%s  - %d",inst->colname[inst->idxRow[i][j]],inst->nzCol[inst->idxRow[i][j]]);
        }
    }
    // printf("\n");

    for( i=1; i<=inst->nCols; i++) {
        inst->idxCol[i] = (int*) malloc((inst->nzCol[i]+1) * sizeof(int));
        inst->coefCol[i] = (double*) malloc((inst->nzCol[i]+1) * sizeof(double));
    }
    //constraint and coefficients that this variable i
    for( i=0; i<inst->nCols; i++) {
        inst->nzCol[i+1]= lp_col(lp,i,inst->idxCol[i+1]+1,inst->coefCol[i+1]+1);
    }
    for( i=nRows+1; i<=inst->nRows; i++) {
        for(j=1; j<=inst->nElRow[i]; j++) {
            inst->nzCol[inst->idxRow[i][j]]++;
        }
    }
    for( i=1; i<=inst->nCols; i++) {
        for(j=1; j<=inst->nzCol[i]; j++) {
            inst->idxCol[i][j]++;
        }
    }

    for( i=nRows+1; i<=inst->nRows; i++) {
        for(j=1; j<=inst->nElRow[i]; j++) {
            inst->idxCol[inst->idxRow[i][j]][nzcol[inst->idxRow[i][j]]+1]=i;
            inst->coefCol[inst->idxRow[i][j]][nzcol[inst->idxRow[i][j]]+1]=inst->coefRow[i][j];
            nzcol[inst->idxRow[i][j]]++;
        }
    }
    /*printf("\n");
    for( i=1; i<=inst->nCols;i++){
        printf("%d- (%d) ",i, inst->nzCol[i]);
        printf(" %s  ",inst->colname[i] );
        for(j=1;j<=inst->nzCol[i];j++){
            printf(" %s  ",inst->Rowname[inst->idxCol[i][j]]);
        }
        printf("\n");
    }*/

    //novos coeficientes de cada linha seguindo a equação Aij Xj<= Bi
    for( i=1; i<=inst->nRows; i++) {
        inst->ncoefrow[i]= (double*)malloc((inst->nElRow[i]+1)* sizeof(double));
    }
    for( i=1; i<=nRows; i++) {
        if(inst->rowType[i]=='1') {
            for( j=1; j<=inst->nElRow[i]; j++) {
                inst->ncoefrow[i][j]= inst->coefRow[i][j]*-1;
            }
            inst->newrhs[i]= inst->rhs[i]*-1;
        }
        else {
            for( j=1; j<=inst->nElRow[i]; j++) {
                inst->ncoefrow[i][j]= inst->coefRow[i][j];
            }
            inst->newrhs[i]= inst->rhs[i];
        }
    }
    for( i=nRows+1; i<=inst->nRows; i++) {
        for(j=1; j<=inst->nElRow[i]; j++) {
            inst->ncoefrow[i][j]= inst->coefRow[i][j]*-1;
        }
        inst->newrhs[i]= inst->rhs[i]*-1;
    }
    /* printf("\n\n");
     for( i=1;i<=inst->nRows;i++){
         printf("%d  |",i);
         for( j=1;j<=inst->nElRow[i];j++){
             printf("%.2f %s   ",inst->ncoefrow[i][j],inst->colname[inst->idxRow[i][j]]);
         }
         printf("  %c   %.2f\n",inst->rowType[i],inst->newrhs[i]);
     }
     printf("\n\n");*/

    // int result = lp_optimize_as_continuous(lp);
    //printf("\nRESULT: %d",result);
    inst->relax= (double*)malloc((inst->nCols+1)* sizeof(double));
    //double *rel = lp_x(lp);

    //coefficients FO
    // for( i=0; i<inst->nCols;i++){
    //     inst->relax[i+1] = rel[i];
//   }

    free (nzcol);
    free (equal);
    free(idx1);
    free (coef1);
    free(idx2);
    free (coef2);
    free(lp);
}

void isntance_free(Instance *inst) {

    int i=0;
    free(inst->obj);
    free(inst->nElRow);
    free(inst->rhs);
    free(inst->rowType);
    free(inst->nzCol);
    free(inst->newrhs);
    free(inst->varsol);
    free(inst->varwrite);
    free(inst->relax);
    //free(inst->consType);
    for(i=1; i<=inst->nRows; i++) {
        free(inst->idxRow[i]);
    }
    free(inst->idxRow);
    for(i=1; i<=inst->nRows; i++) {
        free(inst->coefRow[i]);
    }
    free(inst->coefRow);
    for(i=1; i<=inst->nCols; i++) {
        free(inst->colname[i]);
    }
    free(inst->colname);
    for(i=1; i<=inst->nRows; i++) {
        free(inst->Rowname[i]);
    }
    free(inst->Rowname);
    for(i=1; i<=inst->nCols; i++) {
        free(inst->idxCol[i]);
    }
    free(inst->idxCol);
    for(i=1; i<=inst->nCols; i++) {
        free(inst->coefCol[i]);
    }
    free(inst->coefCol);
    for(i=1; i<=inst->nRows; i++) {
        free(inst->ncoefrow[i]);
    }
    free(inst->ncoefrow);

}

