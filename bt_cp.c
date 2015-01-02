#include "instance.h"
#include "solution.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include "bt_cp.h"
#include "cp.h"

#define EPS 1e-5
#define MAX_DEPTH 200
#define SIZE_MAX_CONF 500
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

#define RMAIORI '1'
#define RMENORI '2'
#define RIGUAL '4'

#define FIX_AT_ZERO 0
#define FIX_AT_ONE 1
#define NO_IMPLICATION 2
#define CONFLICT 3
#define true 1
#define false 0
#define EPS   1e-5


static int n_process=0;
double tRemainder;
int nvarfix=0;
clock_t  startBTt ;

static int lev=0;


void time_start(double t) {
    tRemainder=t;
}

double time_end() {
    tRemainder = ((double) (clock() - startBTt) / CLOCKS_PER_SEC);
    return tRemainder;
}


int startBT(const Instance *inst, Bounds* bd,Bounds * bdc,CPResultStack * cprs,Solution *sol,Solution *solbest,Conflicts * c, Constraints *constr,VUnfixed * vu) {
    int i=0, n_proc=0;
    for(i=1; i<=inst->nCols; i++) {
        inst->varsol[i]=1; // preencho a solução para guiar o BT com o valor 1
    }

    printf("\nNCOLS: %d", inst->nCols);
    printf("\nNROWS: %d", inst->nRows);

    int varfix=1; //primeira variavel a ser passada quando de modo ordenado
    //int varfix= 1 + rand() % (inst.nCols- 1); //primeira variavel a ser passada quando de modo aleatorio
    int valu=1; // valor passado a cada variavel - 1 se  modo ordenado
    // int valu=inst.varsol[varfix]; // valor passado a cada variavel - deste modo quando aleatorio
    int level =0; // nivel do BT
    int control=0; // se control=0, não é inicial, se control=1, é solução inicial
    int rando=1;// se zero ativo ordenado, se 1 ativo aleatorio
    double time= 600.0; // passo o tempo a ser utilizado pelo BT

    time_start(time); // inicio o tempo que o loop gastará
    startBTt = clock(); //pego o tempo corrente

    clock_t  start = clock(); //inicia tempo somente dentro desta função
    b_t_r(inst,bd,bdc,cprs,varfix,valu,sol,solbest,c,level,constr,control,vu,rando); //chama BT com valor primeiro valor
    n_proc= n_proc + n_process; //numero de nos
    bt_clear(inst,vu); //limpa as estruturas referentes ao CP antes de passar
    cp_clear(inst,cprs,bd,bdc,constr); //limpa a estrutura referente ao BT antes de passar
    b_t_r(inst,bd,bdc,cprs,varfix,1-valu,sol,solbest,c,level,constr,control,vu,rando); //chama BT com valor segundo valor
    n_proc= n_proc + n_process; //numero de nos

    clock_t end = clock(); //termina o tempo somente dentro desta função
    double cpuTime = ((double) (end - start)) / CLOCKS_PER_SEC; //calcula o tempo gasto
    printf("\n\ndone in %.3f.\n\n", cpuTime );
    fflush(stdout); // imprime tempo gasto nesta função

    return n_proc;
}

int startBT_Direcionado_Samuel_Combinado(const Instance *inst, Bounds* bd,Bounds * bdc,CPResultStack * cprs,Solution *sol,Solution *solbest,Solution *solstart,Conflicts * c, Constraints *constr,VUnfixed * vu) {
    int i=0,level, varfix, valu, control,rando, n_proc=0;
    double time;

    clock_t  start = clock(); //inicia tempo somente dentro desta função
    clock_t  startRemainder = clock(); //inicia tempo que será o restante para o loop do-while

    //variavel ordenada com valor iniciado com 1
    for(i=1; i<=inst->nCols; i++) {
        inst->varsol[i]=1;//// preencho a solução para guiar o BT com o valor 1
    }
    varfix=1; //primeira variavel a ser passada quando de modo ordenado
    valu=1; // valor passado a cada variavel - 1 se  modo ordenado
    level =0; // nivel do BT
    control=1; // se control=0, não é inicial, se control=1, é solução inicial
    rando=0;// se zero ativo ordenado, se 1 ativo aleatorio
    time= 3.0; // passo o tempo a ser utilizado pelo BT
    time_start(time); // inicio o tempo que o loop gastará
    startBTt = clock(); //pego o tempo corrente
    b_t_r(inst,bd,bdc,cprs,varfix,valu,sol,solbest,c,level,constr,control,vu,rando); //chama BT com valor primeiro valor
    n_proc= n_proc + n_process; //numero de nos
    bt_clear(inst,vu); //limpa as estruturas referentes ao CP antes de passar
    cp_clear(inst,cprs,bd,bdc,constr); //limpa a estrutura referente ao BT antes de passar
    b_t_r(inst,bd,bdc,cprs,varfix,1-valu,sol,solbest,c,level,constr,control,vu,rando); //chama BT com valor segundo valor
    n_proc= n_proc + n_process; //numero de nos
    bt_clear(inst,vu); //limpa as estruturas referentes ao CP antes de passar
    cp_clear(inst,cprs,bd,bdc,constr); //limpa a estrutura referente ao BT antes de passar


    //variavel ordenada com valor iniciado com 0
    for(i=1; i<=inst->nCols; i++) {
        inst->varsol[i]=0;//// preencho a solução para guiar o BT com o valor 1
    }
    varfix=1; //primeira variavel a ser passada quando de modo ordenado
    valu=1; // valor passado a cada variavel - 1 se  modo ordenado
    level =0; // nivel do BT
    control=1; // se control=0, não é inicial, se control=1, é solução inicial
    rando=0;// se zero ativo ordenado, se 1 ativo aleatorio
    time= 3.0; // passo o tempo a ser utilizado pelo BT
    time_start(time); // inicio o tempo que o loop gastará
    startBTt = clock(); //pego o tempo corrente
    b_t_r(inst,bd,bdc,cprs,varfix,valu,sol,solbest,c,level,constr,control,vu,rando); //chama BT com valor primeiro valor
    n_proc= n_proc + n_process; //numero de nos
    bt_clear(inst,vu); //limpa as estruturas referentes ao CP antes de passar
    cp_clear(inst,cprs,bd,bdc,constr); //limpa a estrutura referente ao BT antes de passar
    b_t_r(inst,bd,bdc,cprs,varfix,1-valu,sol,solbest,c,level,constr,control,vu,rando); //chama BT com valor segundo valor
    n_proc= n_proc + n_process; //numero de nos
    bt_clear(inst,vu); //limpa as estruturas referentes ao CP antes de passar
    cp_clear(inst,cprs,bd,bdc,constr); //limpa a estrutura referente ao BT antes de passar


    //variavel aleatoria com valor iniciado com 1
    for(i=1; i<=inst->nCols; i++) {
        inst->varsol[i]=1;//// preencho a solução para guiar o BT com o valor 1
    }
    varfix= 1 + rand() % (inst->nCols- 1); //primeira variavel a ser passada quando de modo aleatorio
    valu=inst->varsol[varfix]; // valor passado a cada variavel - deste modo quando aleatorio
    level =0; // nivel do BT
    control=1; // se control=0, não é inicial, se control=1, é solução inicial
    rando=1;// se zero ativo ordenado, se 1 ativo aleatorio
    time= 3.0; // passo o tempo a ser utilizado pelo BT
    time_start(time); // inicio o tempo que o loop gastará
    startBTt = clock(); //pego o tempo corrente
    b_t_r(inst,bd,bdc,cprs,varfix,valu,sol,solbest,c,level,constr,control,vu,rando); //chama BT com valor primeiro valor
    n_proc= n_proc + n_process; //numero de nos
    bt_clear(inst,vu); //limpa as estruturas referentes ao CP antes de passar
    cp_clear(inst,cprs,bd,bdc,constr); //limpa a estrutura referente ao BT antes de passar
    b_t_r(inst,bd,bdc,cprs,varfix,1-valu,sol,solbest,c,level,constr,control,vu,rando); //chama BT com valor segundo valor
    n_proc= n_proc + n_process; //numero de nos
    bt_clear(inst,vu); //limpa as estruturas referentes ao CP antes de passar
    cp_clear(inst,cprs,bd,bdc,constr); //limpa a estrutura referente ao BT antes de passar

    //variavel aleatoria com valor iniciado com 0
    for(i=1; i<=inst->nCols; i++) {
        inst->varsol[i]=0;//// preencho a solução para guiar o BT com o valor 1
    }
    varfix= 1 + rand() % (inst->nCols- 1); //primeira variavel a ser passada quando de modo aleatorio
    valu=inst->varsol[varfix]; // valor passado a cada variavel - deste modo quando aleatorio
    level =0; // nivel do BT
    control=1; // se control=0, não é inicial, se control=1, é solução inicial
    rando=1;// se zero ativo ordenado, se 1 ativo aleatorio
    time= 3.0; // passo o tempo a ser utilizado pelo BT
    time_start(time); // inicio o tempo que o loop gastará
    startBTt = clock(); //pego o tempo corrente
    b_t_r(inst,bd,bdc,cprs,varfix,valu,sol,solbest,c,level,constr,control,vu,rando); //chama BT com valor primeiro valor
    n_proc= n_proc + n_process; //numero de nos
    bt_clear(inst,vu); //limpa as estruturas referentes ao CP antes de passar
    cp_clear(inst,cprs,bd,bdc,constr); //limpa a estrutura referente ao BT antes de passar
    b_t_r(inst,bd,bdc,cprs,varfix,1-valu,sol,solbest,c,level,constr,control,vu,rando); //chama BT com valor segundo valor
    n_proc= n_proc + n_process; //numero de nos
    bt_clear(inst,vu); //limpa as estruturas referentes ao CP antes de passar
    cp_clear(inst,cprs,bd,bdc,constr); //limpa a estrutura referente ao BT antes de passar


    for(i=1; i<=inst->nCols; i++) {
        solstart->vect[i]=inst->varwrite[i];//pego valor da solução lida .txt e calculo a FO e a Infactibilidade
    }

    CalcConstraints(inst,solbest);// calcula a infactibilidade e custo da melhor solução encontrada pelo BT
    Infeas(inst,solbest); // calcula a infactibilidade e custo da melhor solução encontrada pelo BT
    CalculationCost(inst,solbest); // calcula a infactibilidade e custo da melhor solução encontrada pelo BT
    CalcConstraints(inst,solstart); // calcula a infactibilidade e custo da solução lida pelo .txt
    Infeas(inst,solstart); // calcula a infactibilidade e custo da solução lida pelo .txt
    CalculationCost(inst,solstart);// calcula a infactibilidade e custo da solução lida pelo .txt
    printf("\n START: INF: %.3f Cost: %.3f", solstart->infeas,solstart->cost);
    printf("\n BEST: INF: %.3f Cost: %.3f", solbest->infeas,solbest->cost);
    if(solstart->infeas<solbest->infeas) { // se a inf da solução .txt for melhor que a inf da solbest
        for(i=1; i<=inst->nCols; i++) {
            inst->varsol[i]=solstart->vect[i]; //utlizo a solução .txt para guiar BT
            solbest->vect[i]=solstart->vect[i]; //pego o valor do txt par sol best
        }
        printf("\n START");
    }
    else if(solstart->infeas==solbest->infeas) { // se a inf da solução .txt for igual que a inf da solbest
        if(solstart->cost<solbest->cost) { //verifico se o custo da .txt é melhor que o custo da solbest
            for(i=1; i<=inst->nCols; i++) {
                inst->varsol[i]=solstart->vect[i]; //utlizo a solução .txt para guiar BT
                solbest->vect[i]=solstart->vect[i]; //pego o valor do txt par sol best
            }
            printf("\n START");
        }
        else {
            for(i=1; i<=inst->nCols; i++) {
                inst->varsol[i]=solbest->vect[i]; //utlizo a solução a solbest para guiar BT
            }
            printf("\n best");
        }
    }
    else {
        for(i=1; i<=inst->nCols; i++) {
            inst->varsol[i]=solbest->vect[i]; //utlizo a solução .txt para guiar BT
        }
        printf("\n best");
    }

    // faço o BT guiado
    time = 300.0-((double) (clock() - startRemainder) / CLOCKS_PER_SEC); // passo o tempo a ser utilizado pelo BT
    time_start(time); // inicio o tempo que o loop gastará
    startBTt = clock(); //pego o tempo corrente
    control=0; // se control=0, não é inicial, se control=1, é solução inicial
    rando=0;// se zero ativo ordenado, se 1 ativo aleatorio
    do {
        varfix=1;
        //varfix= 1 + rand() % (inst->nCols- 1); //primeira variavel a ser passada quando de modo aleatorio
        valu=inst->varsol[varfix]; // valor passado a cada variavel - deste modo quando aleatorio
        level =0; // nivel do BT
        b_t_r(inst,bd,bdc,cprs,varfix,valu,sol,solbest,c,level,constr,control,vu,rando); //chama BT com valor primeiro valor
        n_proc= n_proc + n_process; //numero de nos
        bt_clear(inst,vu); //limpa as estruturas referentes ao CP antes de passar
        cp_clear(inst,cprs,bd,bdc,constr); //limpa a estrutura referente ao BT antes de passar
        b_t_r(inst,bd,bdc,cprs,varfix,1-valu,sol,solbest,c,level,constr,control,vu,rando); //chama BT com valor segundo valor
        n_proc= n_proc + n_process; //numero de nos
        bt_clear(inst,vu); //limpa as estruturas referentes ao CP antes de passar
        cp_clear(inst,cprs,bd,bdc,constr); //limpa a estrutura referente ao BT antes de passar
        for(i=1; i<=inst->nCols; i++) {
            inst->varsol[i]=solbest->vect[i];
        }
    } while((((double) (clock() - startBTt) / CLOCKS_PER_SEC)  <  tRemainder )); // enquanto tiver tempo


    clock_t end = clock(); //termina o tempo somente dentro desta função
    double cpuTime = ((double) (end - start)) / CLOCKS_PER_SEC; //calcula o tempo gasto
    printf("\n\ndone in %.3f.\n\n", cpuTime );
    fflush(stdout); // imprime tempo gasto nesta função


    return n_process;
}

//constroi, aloca e inicia a estrutura VUnfixed que é utilizada quando as variaveis do BT são passadas de forma aleatoria
void vUfixed_create(const Instance *inst,VUnfixed* vu) {

    int i=0;
    vu->var= (int*) malloc((inst->nCols+1) * sizeof(int));
    vu->indexv= (int*) malloc((inst->nCols+1) * sizeof(int));

    for( i=1; i<=inst->nCols; i++) {
        vu->var[i]= i;
        vu->indexv[i]=i;
    }
    vu->nunfix=inst->nCols;
}

//função para limpar as estruturas antes de chamar novamente o BT
void bt_clear(const Instance *inst,VUnfixed* vu) {

    int i;

    for( i=1; i<=inst->nCols; i++) {
        vu->var[i]= i;
        vu->indexv[i]=i;
    }
    vu->nunfix=inst->nCols;

    n_process=0;
    nvarfix=0;

}
// Verifica se a nova solução corrente é melhor que a melhor encontrada
//se for melhor copia e salva
void solBetter(const Instance *inst, Bounds * bd,Solution *sol, Solution *solbest) {
    int j=0;
    for ( j=1; j<=inst->nCols; j++) {
        sol->vect[j]=(char) bd->l[j];
    }

    CalcConstraints(inst,sol);// calculo as restrições
    Infeas(inst, sol); // calcula infactibilidade
    CalculationCost(inst, sol); //calcula custo
    int update= SolIsBetter(sol,solbest);/*verificando se melhorou solução*/
    if(update==1) {
      //  printf("\nsolBest: %.4f cost: %.4f\n",sol->infeas,sol->cost);
        CopySolution(inst,solbest,sol); // se melhor copia
    }
}

//copia solução parao BT guiado por solução
void copy_sol(Instance *inst, Solution *solbest) {

    int i=1;
    for(i=1; i<=inst->nCols; i++) {
        inst->varsol[i]=solbest->vect[i];
    }

}

//função utilizada para remover variavel quando as variaveis do BT são passadas de modo aleatorio
void remove_var(const Instance *inst,VUnfixed * vu, int var) {
    //  printf("\n REMOVE VAR: %d",var);
    int position = vu->indexv[var];
    int var_realoc= vu->var[vu->nunfix];
    vu->var[position]= vu->var[vu->nunfix];
    vu->var[vu->nunfix]=0;
    vu->indexv[var_realoc]= position;
    vu->indexv[var]=0;
    vu->nunfix--;
}

//função utilizada para remover variavel quando as variaveis do BT são passadas de modo aleatorio
void add_var(const Instance *inst,VUnfixed * vu, int var) {
    // printf("\n ADD VAR: %d",var);
    vu->nunfix++;
    vu->var[vu->nunfix]= var;
    vu->indexv[var]= vu->nunfix;
}
// função back tracking
void b_t_r(const Instance *inst, Bounds * bd,Bounds * bdc,CPResultStack * cprs,int varfix, int value,Solution *sol,Solution *solbest,Conflicts * c, int level, Constraints *constr, int control,VUnfixed * vu, int rando) {


    int i,j,s,r=0,val,position,pMax=0,pos;
    int nfixC=0;
    int *upBouds; //vetor das variaveis fixadas pelo CP , utilizado para retornar
    n_process++; // incremento numero de nos processados
    nvarfix++; // incremento o numero de variaveis fixadas



    level ++; // incremento nivel do BT
    // printf("\nBT___________________Nivel: %d  - VF: %d  Name: %s  Valor: %d", level,varfix, inst->colname[varfix], value);
    int vari= varfix; // pego a variavel fixada

    if(lev<level) {
        lev=level;
    }

    //fixa valor da variavel vetores L e U
    bd->l[varfix]=value;
    bd->u[varfix]=value;
    remove_var(inst,vu,vari);
    if ((((double) (clock() - startBTt) / CLOCKS_PER_SEC)  >  tRemainder )) {
        //  if ((((double) (clock() - startBTt) / CLOCKS_PER_SEC)  >  tRemainder )||(n_process>inst->nCols*15)){
        goto TERMINATE;
    }

    //testa se já possuo solução e se possui verifica solução
    if(nvarfix==inst->nCols) { //se possui soloução chama função preenche sol e verifica
        solBetter(inst,bd,sol,solbest);
        //printf("\nCOMPLETE SOLUTION 1 INF:%.4f     COST:%.4f  NIVEL: %d",sol->infeas, sol->cost,nivel);
        if(control==1) {
            tRemainder=0.0;
        }
        goto TERMINATE;
    }//termina loop solução

    constraints_clear(inst,constr);//limpa o CP
    //cprs_clear(inst,cprs); // limpa o cprs
    //preenche as retrições a serem passadas
    for(j=1; j<=inst->nzCol[varfix]; j++) { //loop no numero de restrições que a variavel fixada aparece
        r=inst->idxCol[varfix][j]; // pega a restrição que a variavel fixada aparece
        constraint_add(inst,constr,r); // se a restrição possui variavel não fixada, adiciona restrição na estrutura contraint
    }

    bounds_copy(inst,bd,bdc);//copia o cp para o cpbt

    s= cPropagation (inst,bdc,cprs,level,constr,tRemainder,startBTt);// chama constraint propagation
    // printf("\nRetorno: %d",s);
    if(s == CONFLICT) {
        conflicts_add(inst,c,cprs);
        cprs_clear(inst,cprs); // limpa o cprs
        goto TERMINATE; // se o resultado do CP for conflito, retorna
    }
    // se possui fixações atualizo vetor l e u
    else if ((s==FIX_AT_ZERO) || (s==FIX_AT_ONE)) {
        int t=0;
        if(cprs->tfix>0) {
            upBouds = (int*) malloc((cprs->tfix+1) * sizeof(int));//aloca vetor para armazenar as variaveis fixadas
            pMax=cprs->depth; // armazena quantos nivies do CP tiveram fixações na chamada
            nfixC=cprs->tfix; // armazena numero de variaveis fixadas na chamada
            for(i=0; i<=cprs->depth; i++) { // percorro os niveis de fixações do CP
                for(j=0; j<cprs->cpr[i].nfix; j++) { // percorre as fixações em cada nivel do CP
                    position = cprs->cpr[i].idxfix[j]; // pego a variavel fixada
                    val = cprs->cpr[i].valfix[j]; // pego o valor da variavel fixada
                    bd->u[position]=val;// atualizo o upper bound da posição da variavel fixada com o valor fixado
                    upBouds[t]=position;
                    nvarfix++; // incremento o numero de variaveis fixadas
                    t++;
                    remove_var(inst,vu,position);
                }
            }
        }
        cprs_clear(inst,cprs); // limpa o cprs
    }
    // se retornou NO_IMPLICATION mas possui fixações, atualizo vetor l e u
    else if (s==NO_IMPLICATION) {
        int t=0;
        // printf(" tfix: %d",tfix);
        if(cprs->tfix>0) {
            upBouds = (int*) malloc((cprs->tfix+1) * sizeof(int));//aloca vetor para armazenar as variaveis fixadas
            pMax=cprs->depth; // armazena quantos nivies do CP tiveram fixações na chamada
            nfixC=cprs->tfix; // armazena numero de variaveis fixadas na chamada
            for(i=0; i<=cprs->depth; i++) { // percorro os niveis de fixações do CP
                for(j=0; j<cprs->cpr[i].nfix; j++) { // percorre as fixações em cada nivel do CP
                    position = cprs->cpr[i].idxfix[j]; // pego a variavel fixada
                    val = cprs->cpr[i].valfix[j]; // pego o valor da variavel fixada
                    bd->l[position]=val; // atualizo o lower bound da posição da variavel fixada com o valor fixado
                    bd->u[position]=val;// atualizo o upper bound da posição da variavel fixada com o valor fixado
                    upBouds[t]=position;
                    nvarfix++; // incremento o numero de variaveis fixadas
                    t++;
                    remove_var(inst,vu,position);
                }
            }
        }
        cprs_clear(inst,cprs); // limpa o cprs
    }
    //testa se já possuo solução e se possui verifica solução
    if(nvarfix==inst->nCols) { //se possui soloução chama função preenche sol e verifica
        solBetter(inst,bd,sol,solbest);
        if(control==1) {
            tRemainder=0.0;
        }
        goto TERMINATE;
    }//termina loop solução

    //pega proxima variavel ordenada sem fixar
    if(rando==0) {
        for ( i = 1; i <= inst->nCols; i++ ) {
            if(bd->l[i]!=bd->u[i]) { // se não possui a variavel na posição i não possui fixação
                varfix= i; //pego a variavel
                break;
            }
        }
    }
    //pega proxima variável aleatoria
    if(rando ==1) {
        if( vu->nunfix>1) {
            int tam =  vu->nunfix-1;
            pos = rand() % tam;
            pos ++;
        }
        else {
            pos =1;
        }
        varfix = vu->var[pos];
    }

    //if ((((double) (clock() - startBTt) / CLOCKS_PER_SEC)  >  tRemainder )||(n_process>inst->nCols*15)){
    if ((((double) (clock() - startBTt) / CLOCKS_PER_SEC)  >  tRemainder )) {
        goto TERMINATE;
    }

    int valu= inst->varsol[varfix];
    // printf("\nVarfix: %d",varfix);
    b_t_r(inst,bd,bdc,cprs,varfix,valu,sol,solbest,c,level,constr, control,vu,rando);// chama BT com valor 1
    b_t_r(inst,bd,bdc,cprs,varfix,1-valu,sol,solbest,c,level,constr,control,vu,rando); //chama BT com valor 0

TERMINATE:
    bd->l[vari]=0; // desfaz a variavel fixada pelo BT
    bd->u[vari]=1;// desfaz a variavel fixada pelo BT
    add_var(inst,vu,vari);//adiociona variavel na estrutura para variaveis aleatorias
    nvarfix--;

    if(pMax>0) { // se forem feitas fixações pelo CP
        for(i=0; i<nfixC; i++) { // percorro a estrutura que armazenou as fixações

            position= upBouds[i]; // pego a variavel fixada
            bd->l[position] = 0; // libero a variavel
            bd->u[position] = 1; // libera a variavel
            add_var(inst,vu,position); //adiociona variavel na estrutura para variaveis aleatorias
            nvarfix--; // decremento o numero de variaveis fixads
        }
        free(upBouds); // libera a estrutura de fixações realizadas pelo CP
    }
    return;
}

// nivel de profundidade maximo atigindo pelo BT
int level() {
    return lev;
}

int SolIsInfac( const Instance *inst, Bounds* bd,int varfix) {

    int j,r,h;
    double lhs;
    int ret=0;
    int teste=0;

    for(j=1; j<=inst->nzCol[varfix]; j++) { //loop no numero de restrições que a variavel fixada aparece

        r=inst->idxCol[varfix][j]; // pega a restrição que a variavel fixada aparece
        //  printf("\nR: %s, var: %s",inst->Rowname[r], inst->colname[varfix]);
        teste=0;
        for(h=1; h<= inst->nElRow[r]; h++) { // loop nas variaveis da restrição
            int var= inst->idxRow[r][h]; //pega uma variavel
            if(bd->l[var]!=bd->u[var]) { //verifica se não foi fixada
                teste=1;
                break;
            }
        }
        lhs=0;
        if(teste==0) {
            for(h=1; h<= inst->nElRow[r]; h++) { // loop nas variaveis da restrição
                ret=9999999;
                lhs=lhs + inst->ncoefrow[r][h]*((double)bd->l[inst->idxRow[r][h]]);
            }
            //printf("\n LHS: %.2f  RHS:  %.2f ",lhs,inst->newrhs[r]);
            if(lhs>inst->newrhs[r]) {
                //  printf("\nREstr: %s  type: %c LHS: %.2f  RHS:  %.2f ,   :",inst->Rowname[r],inst->rowType[r],lhs,inst->newrhs[r]);
                for(h=1; h<= inst->nElRow[r]; h++) { // loop nas variaveis da restrição
                    //   printf("\n%s: %d %d var %d ",inst->colname[inst->idxRow[r][h]],bd->l[inst->idxRow[r][h]],bd->u[inst->idxRow[r][h]],inst->idxRow[r][h]);
                    //   getchar();
                }
                return r;
            }
            else {
                //   printf("  Calc OK!");
            }
        }
    }
    return ret;
}

//desaloca a esturutura vUfixed
void vUfixed_free(VUnfixed* vu) {

    free (vu->var);
    free (vu->indexv);

}
