#include "instance.h"
#include "solution.h"
#include <stdio.h>
#include <stdlib.h>
#include <glpk.h>
#include <math.h>

#define RMAIORI '1'
#define RMENORI '2'
#define RIGUAL '4'
/*Solution Starts*/
void SolutionStarts(Instance *inst, Solution *sol) {
    int i=1;

    /* alocation*/
    sol->vect = (char*) malloc((inst->nCols+1) * sizeof(char));
    sol->lhs = (double*) malloc((inst->nRows+1) * sizeof(double));
    sol->inf = (double*) malloc((inst->nRows+1) * sizeof(double));
    sol->newinf = (double*) malloc((inst->nRows+1) * sizeof(double));

    /*start variables*/

    for( i=1; i<=inst->nCols; i++) {
        sol->vect[i]= 1;
    }
    for(i=1; i<=inst->nRows; i++) {
        sol->inf[i]= 0.0;
    }
    for(i=1; i<=inst->nRows; i++) {
        sol->lhs[i]= 0.0;
    }


    for(i=1; i<=inst->nRows; i++) {
        sol->newinf[i]= 0.0;
    }
    sol->infeas=99999.0;
    sol->cost= 999999.0;
}

/*generate random vector solution*/
void GenerateRandom(const Instance *inst, Solution *sol) {
    int i;

    /*generate start solution*/
    for ( i=1; i<=inst->nCols; i++) {
        sol->vect[i]= rand() % 2;
    }
}



/*calculation the constraints according to the solution vector*/
void CalcConstraints(const Instance *inst, Solution *sol) {
    int i=1,j=1;
    double z=0.0;

    /* calcula rhs*/
    for( i=1; i<=inst->nRows; i++) {
        for( j=1; j<=inst->nElRow[i]; j++) {
            z=z + inst->coefRow[i][j]*((double)sol->vect[inst->idxRow[i][j]]);
        }
        sol->lhs[i]=z;
        z=0.0;
    }
}

/*calculate the distance of the restrictions and limits*/
void Infeas(const Instance *inst, Solution *sol) {
    double d=0.0;
    int i=1;
    sol->infeas=0.0;

    for(i=1; i<=inst->nRows; i++) {
        d=0.0;
        switch (inst->rowType[i]) {
        case RMAIORI:
            if(sol->lhs[i]>=inst->rhs[i]) {
                d=0.0;
            }
            else {
                d= fabs(inst->rhs[i]-sol->lhs[i]);
            }
            break;
        case RMENORI:
            if(sol->lhs[i]<=inst->rhs[i]) {
                d=0.0;
            }
            else {
                d= fabs(sol->lhs[i]-inst->rhs[i]);
            }
            break ;
        case RIGUAL:
            if(sol->lhs[i]==inst->rhs[i]) {
                d=0.0;
            }
            else {
                d= fabs(inst->rhs[i]-sol->lhs[i]);
            }
            break ;
        }
        sol->inf[i]=d;
        sol->infeas= sol->infeas + d;
    }
}

/*calculation  Cost*/
void CalculationCost(const Instance *inst, Solution *sol) {
    int i=1;
    sol->cost=0.0;

    for (i = 1; i <= inst->nCols; i++) {
        sol->cost = sol->cost + (inst->obj[i]* sol->vect[i]);

    }
}


void newInfeas(const Instance *inst,Solution *sol, int pos) {
    int i=0;
    double d=0.0;
    for(i=1; i<=inst->nRows; i++) {
        sol->newinf[i]= sol->inf[i];
    }
    for(i=1; i<=inst->nzCol[pos]; i++) {
        d=0.0;
        switch (inst->rowType[inst->idxCol[pos][i]]) {

        case RMAIORI:
            if(sol->lhs[inst->idxCol[pos][i]]>=inst->rhs[inst->idxCol[pos][i]]) {
                d=0.0;
            }
            else {
                d= fabs(inst->rhs[inst->idxCol[pos][i]]-sol->lhs[inst->idxCol[pos][i]]);
            }
            break;
        case RMENORI:
            if(sol->lhs[inst->idxCol[pos][i]]<=inst->rhs[inst->idxCol[pos][i]]) {
                d=0.0;
            }
            else {
                d= fabs(sol->lhs[inst->idxCol[pos][i]]-inst->rhs[inst->idxCol[pos][i]]);
            }
            break ;
        case RIGUAL:
            if(sol->lhs[inst->idxCol[pos][i]]==inst->rhs[inst->idxCol[pos][i]]) {
                d=0.0;
            }
            else {
                d= fabs(inst->rhs[inst->idxCol[pos][i]]-sol->lhs[inst->idxCol[pos][i]]);
            }
            break ;
        }
        sol->newinf[inst->idxCol[pos][i]]=d;
        sol->infeas = sol->infeas + d - sol->inf[inst->idxCol[pos][i]];
    }
    for(i=1; i<=inst->nRows; i++) {
        sol->inf[i]= sol->newinf[i];
    }
}

void CalculationNewCost(const Instance *inst, Solution *sol, int pos) {
    sol->cost = sol->cost + (inst->obj[pos]* sol->vect[pos])- inst->obj[pos] *(1- sol->vect[pos]);
}

void CopySolution(const Instance *inst, Solution *target, const Solution *source ) {
    int h=1;


    for (h = 1; h <= inst->nCols; h++) {
        target->vect[h] = source->vect[h];
    }
    for (h = 1; h <= inst->nRows; h++) {
        target->lhs[h] = source->lhs[h];
    }
    for (h = 1; h <= inst->nRows; h++) {
        target->inf[h] = source->inf[h];
    }
    target->infeas=source->infeas;
    target->cost = source->cost;

}



int SolIsBetter( const Solution *solLS, const Solution *sol) {
    int update;
    if(solLS->infeas<sol->infeas) {
        update=1;
    }
    else if(solLS->infeas==sol->infeas) {
        if(solLS->cost<sol->cost) { //muda sinal
            update=1;
        }
        else {
            update=0;
        }
    }
    else {
        update=0;
    }
    return update;
}




void liberasol(Solution *sol,const Instance *inst) {
    free(sol->vect);
    free(sol->lhs);
    free(sol->inf);
    free(sol->newinf);
}




