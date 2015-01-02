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

//constroi, aloca e inicia com 0 a estrutura cpResultStack - que contem o resultado das fixações
void cprs_create(const Instance *inst,CPResultStack * cprs) {
    // aloca a estrutura
    cprs->depth_max= MAX_DEPTH;
    cprs->cpr = (CPResult*) malloc((MAX_DEPTH) * sizeof(CPResult));

    int i=0,j=0;
    for( i=0; i<MAX_DEPTH; i++) {
        cprs->cpr[i].idxfix = (int*) malloc((inst->nCols+1) * sizeof(int));
        cprs->cpr[i].valfix = (int*) malloc((inst->nCols+1) * sizeof(int));
    }
    //inicializa a estrutura com 0
    for( i=0; i< MAX_DEPTH ; i++) {
        cprs->cpr[i].result= NO_IMPLICATION;
        cprs->cpr[i].nfix=0;
        for( j=1; j<inst->nCols+1; j++) {
            cprs->cpr[i].idxfix[j] = 0;
            cprs->cpr[i].valfix[j] = 0;
        }
    }
    cprs->depth=0;

}

//constroi, aloca e inicia com 0 a estrutura Constraints que contem as restrições a serem passadas para o CP
void constraints_create(const Instance *inst,Constraints * constr) {

    int i=0;
    //aloca
    constr->c_indx = (int*) malloc((inst->nRows+1) * sizeof(int));
    constr->c = (int*) malloc((inst->nRows+1) * sizeof(int));

    //inicia com 0
    for( i=1; i<inst->nRows; i++) {
        constr->c_indx[i]=0;
        constr->c[i]=0;
    }
    constr->nconst=0;
}

//constroi, aloca e inicia a estrutura CPBT - que contem os vetores de limites L e U
void bounds_create(const Instance *inst, Bounds * bd) {

    int i=0;
    // aloca
    bd->l = (int*) malloc((inst->nCols+1) * sizeof(int));
    bd->u = (int*) malloc((inst->nCols+1) * sizeof(int));
    //lower bounds=0, upper bounds=1
    for( i=1; i<=inst->nCols; i++) {
        bd->l[i]=0;
        bd->u[i]=1;
    }
}

//controi e aloca a estrutura Conflicts que contem os pares de conflitos econtrados
void conflicts_create(const Instance *inst,Conflicts * c) {
    int i,j;
    c->nconf = (int*) malloc((inst->nCols+1) * sizeof(int));
    c->ncmax = (int*) malloc((inst->nCols+1) * sizeof(int));
    c->confl = (Conflict**) malloc((inst->nCols+1) * sizeof(Conflict*));
    for(i=1; i<=inst->nCols; i++) {
        c->confl[i] = (Conflict*) malloc((SIZE_MAX_CONF+1) * sizeof(Conflict));
    }
    for(i=1; i<=inst->nCols; i++) {
        c->nconf[i]=0;
        c->ncmax[i]= SIZE_MAX_CONF;
        for(j=0; j<=SIZE_MAX_CONF; j++) {
            c->confl[i][j].idx1=99999;
            c->confl[i][j].idx2=99999;
            c->confl[i][j].val1=99999;
            c->confl[i][j].val2=99999;
        }
    }
}

//função para zerar as estruturas antes de chamar novamente  o BT com 0 no main
void cp_clear(const Instance *inst,CPResultStack * cprs,Bounds * bd,Bounds * bdc,Constraints *constr) {

    int i,j;
    for( i=0; i< MAX_DEPTH ; i++) {
        cprs->cpr[i].result= NO_IMPLICATION;
        cprs->cpr[i].nfix=0;
        for( j=1; j<inst->nCols+1; j++) {
            cprs->cpr[i].idxfix[j] = 0;
            cprs->cpr[i].valfix[j] = 0;
        }
    }
    for( i=1; i<=inst->nCols; i++) {
        bd->l[i]=0;
        bd->u[i]=1;
        bdc->l[i]=0;
        bdc->u[i]=1;
    }
    for( i=1; i<=inst->nRows; i++) {
        constr->c_indx[i]=0;
        constr->c[i]=0;
    }
    constr->nconst=0;
}
//limpa a estrutura CPBT, antes de passar para um novo CP
void constraints_clear(const Instance *inst, Constraints *constr) {
    int i=0;

    for( i=1; i<=inst->nRows; i++) {
        constr->c_indx[i]=0;
        constr->c[i]=0;
    }
    constr->nconst=0;
}

//limpa a estrutura cpResultStack, antes de passar para a estrutura CP
void cprs_clear(const Instance *inst,CPResultStack * cprs) {
    int i=0,j=0;

    for( i=0; i<= cprs->depth ; i++) {
        cprs->cpr[i].result= NO_IMPLICATION;
        cprs->cpr[i].nfix=0;
        for( j=1; j<inst->nCols+1; j++) {
            cprs->cpr[i].idxfix[j] = 0;
            cprs->cpr[i].valfix[j] = 0;
        }
    }
}

//adiciona uma nova restrição a estrutura Constraints
void constraint_add(const Instance *inst, Constraints *constr, int r) {
    constr->nconst++; //acrescento o numero de restrições
    constr->c[constr->nconst]=r; //adiciono a restrição no final do vetor
    constr->c_indx[r]= constr->nconst; //coloco o indice que a restrição foi fixada na posição da restrição
}

//remove uma restrição da estrutura Constraints
void constraint_remove(const Instance *inst, Constraints *constr, int r) {

    constr->c_indx[constr->nconst]= constr->c_indx[r];//pega o indice da restrição que será retirada
    int pos= constr->c_indx[r]; //pego o indice que a restrição trocada será colocada
    constr->c_indx[r]=0; // na posição da restrição retirada coloca 0
    constr->c[pos]= constr->c[constr->nconst]; //coloca a ultima restrição no lugar da restrição retirada
    constr->c[constr->nconst]=0; // coloca zero na ultima posição
    constr->nconst--; //incremento o numero de restrições
}

//copia uma estrutura CPBT para outra CPBT
void bounds_copy(const Instance *inst, Bounds * bd,Bounds *bdc) {
    int i=0;

    for( i=1; i<=inst->nCols; i++) {
        bdc->l[i]=bd->l[i];
        bdc->u[i]=bd->u[i];
    }
}

//Adiciona os pares de conflitos encontrados na estrutura Conflicts
void conflicts_add(const Instance *inst,Conflicts * c, CPResultStack * cprs) {
    int i,p1,p2,col1,col2,v1,v2;
    int test1=0,test2=0,j,k;

    for(i=0; i<= cprs->depth; i++) { // loop em todos os niveis de conflitos

        if(cprs->cpr[i].nfix==2) { // verifico se no nivel foram encontrados 2 conflitos

            test1=1;
            test2=1;
            p1 = cprs->cpr[i].idxfix[0]; //indice da var 1
            v1 = cprs->cpr[i].valfix[0]; // valor da var 1
            p2 = cprs->cpr[i].idxfix[1]; // undice da var 2
            v2 = cprs->cpr[i].valfix[1]; //valor da var 2
            col1= c->nconf[p1]; // numero de conflitos da var 1
            col2= c->nconf[p2]; //munero de conflitos da var 2
            if(c->nconf[p1]>c->ncmax[p1]) { // Verifico se já ultrapassou o num max de conflitos encontrados
                c->ncmax[p1]= c->ncmax[p1]*2;
                c->confl[p1] = (Conflict *)realloc (c->confl[p1], (c->ncmax[p1]+1)*sizeof(Conflict)); // se sim realoca
            }
            if(c->nconf[p2]>c->ncmax[p2]) { // Verifico se já ultrapassou o num max de conflitos encontrados
                c->ncmax[p2]= c->ncmax[p2]*2;
                c->confl[p2] = (Conflict *)realloc (c->confl[p2], (c->ncmax[p2]+1)*sizeof(Conflict));// se sim realoca
            }
            for(j=0; j<c->nconf[p1]; j++) { // verifico se já foi armazenado par de conflitos iguais
                if((c->confl[p1][j].idx2)==p2) {
                    if((c->confl[p1][j].val2)==v2) {
                        if((c->confl[p1][j].val1)==v1) {
                            test1=0;
                        }
                    }
                }
            }
            for(k=0; k<c->nconf[p2]; k++) { // verifico se já foi armazenado par de conflitos iguais
                if((c->confl[p2][k].idx2)==p1) {
                    if((c->confl[p2][k].val2)==v1) {
                        if((c->confl[p2][k].val1)==v2) {
                            test2=0;
                        }
                    }
                }
            }
            if(test1>=1) { // se não foi armazenado ainda armazeno P1
                c->confl[p1][col1].idx1=p1;
                c->confl[p1][col1].val1=v1;
                c->confl[p1][col1].idx2=p2;
                c->confl[p1][col1].val2=v2;
                c->nconf[p1]++;
            }
            if(test2>=1) { // se não foi armazenado ainda armazeno P2
                c->confl[p2][col2].idx1=p2;
                c->confl[p2][col2].val1=v2;
                c->confl[p2][col2].idx2=p1;
                c->confl[p2][col2].val2=v1;
                c->nconf[p2]++;
            }
        }
    }
}

// Adiciona as fixações realizadas no CP- recebe pos- nivel do CP, var= variavel, val = valor da variavel
void cprs_add_result(const Instance *inst,int pos, int var, int val,CPResultStack * cprs) {
    int i=0,j=0;
    if(pos>=cprs->depth_max) { // verifico se ja ultrapassou o nivel maximo
        cprs->cpr = (CPResult *)realloc (cprs->cpr, (cprs->depth_max *2)*sizeof(CPResult)); // sem sim realoco

        for( i=cprs->depth_max; i<(cprs->depth_max*2); i++) { // aloca os vetores internos
            cprs->cpr[i].idxfix = (int*) malloc((inst->nCols+1) * sizeof(int));
            cprs->cpr[i].valfix = (int*) malloc((inst->nCols+1) * sizeof(int));
        }
        //inicializo as novas alocações
        for( i=cprs->depth_max; i< (cprs->depth_max*2) ; i++) { //inicio
            cprs->cpr[i].result= NO_IMPLICATION;
            cprs->cpr[i].nfix=0;
            for( j=1; j<inst->nCols+1; j++) {
                cprs->cpr[i].idxfix[j] = 0;
                cprs->cpr[i].valfix[j] = 0;
            }
        }
        cprs->depth_max= cprs->depth_max*2;
    }
    else { // se não ultrapassou preencho
        cprs->cpr[cprs->depth].idxfix[cprs->cpr[cprs->depth].nfix]=var;//insere a variavel fixada
        cprs->cpr[cprs->depth].valfix[cprs->cpr[cprs->depth].nfix]=val; //insee o valor da variavel fixada
        cprs->cpr[cprs->depth].result= val; //insere o resultado da variavel fixada
        cprs->cpr[cprs->depth].nfix ++; //incremento o numero de fixações
    }
}

// verifica se no nivel atual esta restrição ja foi fixada, e se foi se ela foi fixada com o mesmo valor da nova chamada do evalueteBound
int cprs_check(const Instance *inst,CPResultStack * cprs, int var) {
    int h=0;
    int variable=-1;
    for(h=0; h<cprs->cpr[cprs->depth].nfix; h++) { // verifica se no nivel atual esta restrição ja foi fixada
        if(var==cprs->cpr[cprs->depth].idxfix[h]) {
            variable=h; // se estiver fixa, pego posição que foi fixada
            break;
        }
    }
    return variable;
}

//calucula o valor do limite U
double calculatedU(const Instance *inst, Bounds * bd, int r, int var, int nivel) {
    int i=0;
    int index;
    double pos=0.0;
    double lim=0.0;
    int bound;

    for(i=1; i<= inst->nElRow[r]; i++) { // loop de 1 ate o numero de elementos da restrição r
        if(inst->idxRow[r][i]!=var) { // verifica se não estou pegando a variavel var na retrição r
            index=inst->idxRow[r][i]; // pega o indice da variavel variavel da restrição r
            if(bd->l[index]==bd->u[index]) { //se a variavel ja estiver fixa, pego o limite fixado
                bound= bd->l[index]; //pego o limite
            }
            else if(inst->ncoefrow[r][i]>=0.0) { // se não foi fixada, e se for positiva o limite é 0
                bound= 0;
            }
            else { //se não foi fixada, e for negativa o limite é 1
                bound = 1;
            }
            pos = pos + (inst->ncoefrow[r][i]* ((double)bound)); //multiplica o coeficiente pelo limite e acrescenta na soma
        }
    }
    //  printf("   pos: %.3f",pos);
    lim = inst->newrhs[r]- (pos); //diminuo o limite da direita dos calculados
    return lim;// retorna o LImite Uij
}

// função que verifica se fixo em 1 ou 0 ou conflito ou não sei em que fixar
int evalueteBound (const Instance *inst, Bounds * bd, int r, int var, int nivel) {

    int bound=0;
    double lim=0.0;
    int h;

    lim=calculatedU(inst,bd,r,var,nivel); //pega o limite Uij
//   printf("   U: %f",lim);
    for(h=1; h<= inst->nElRow[r]; h++) { // percorre a restrição r

        if(inst->idxRow[r][h]==var) { //percorre a restrição para achar a variavel

            //conjunto positivos
            if(inst->ncoefrow[r][h]>=0.0) { //verifica se coeficinte da var é positivo
                //      printf("   inst->ncoefrow[r][h] %f, lim %f", inst->ncoefrow[r][h],lim);
                //   printf("\npositivos");
                if(lim<0.0) { // se limite menor que 0
                    //  printf("   CONFLITO %d",var);
                    bound=3; //conflito
                }
                else if (inst->ncoefrow[r][h] >lim+EPS) { //se não se coeficiente de var >lim
                    bound=0; //fixo var em 0
                    //    printf("   fixou %s: %d",inst->colname[var], 0);
                }
                else {
                    bound=2; // não sei em que fixar
                    //  printf("   Não coneguiu fixar");
                }

            }
            //conjunto negativos
            else { //verifica se coeficinte da var é negativo
                //  printf(" \nnegativos");
                if(inst->ncoefrow[r][h]> lim+EPS) { //verifica se coeficiente é maior que lim
                    // printf("   COMFLITO  %d:",var);
                    bound=3;
                }
                else if (lim<0) { // se não se lim<0
                    bound=1; // fixo var em 1
                    // printf("   fixou %d: %d", var, 1);
                }
                else {
                    bound=2;// se não não sei em que fixar
                    //  printf("   Não coneguiu fixar");
                }
            }
            h=inst->nElRow[r]+1; //sai do loop quando acha a var
        }
    }
    return bound; //retorna o bound calculado
}


int calc_lhs(const Instance *inst, Bounds * bd, int r) {

    int h;
    double lhs=0;
    for(h=1; h<= inst->nElRow[r]; h++) { // loop nas variaveis da restrição
        lhs=lhs + inst->ncoefrow[r][h]*((double)bd->l[inst->idxRow[r][h]]);
    }
    if(lhs>inst->newrhs[r]) {
        return CONFLICT;
    }
    else {
        return NO_IMPLICATION;
    }
}

// Função Constraint Proapagation
int cPropagation (const Instance *inst, Bounds * bd, CPResultStack * cprs, int nivel, Constraints *constr, double tRemainder ,clock_t  startBTt ) {
    //printf("\nCP____________________________________%d", constr->nconst);
    int i,j,var, variable, restr;
    int c_unfix=0,r=0;
    int bound =0;
    int new_implication=false;

    cprs->depth=0;
    cprs->tfix=0; //inicializo total de fixações com 0
    do {
        // printf("\nCP____________________________________%d", constr->nconst);
        new_implication=false; // incializo a cada loop sem novas imlicações
        cprs->depth++; // incremento o nivel
        cprs->cpr[cprs->depth].result= NO_IMPLICATION; //inicializo o nivel p
        cprs->cpr[cprs->depth].nfix= 0; //inicializo o nivel p
        for(i=1; i<=constr->nconst; i++) { // loop em todas as retrições para ver qual devo analizar
            restr= constr->c[i]; //pego a restrição que será verificada
            c_unfix=0;
            for(j=1; j<= inst->nElRow[restr]; j++) { // loop em todos os elementos que a variavel aparece
                var= inst->idxRow[restr][j]; // variavel corrente
                if(bd->l[var]!=bd->u[var]) { // verifico se não foi fixad
                    c_unfix=1; // se não foi fixada c_unfix= 1, isto é tenho queverificar LHS
                }
            }
            if(c_unfix<=0) { // verifica lhs
                int check=calc_lhs(inst,bd,restr);  // calcula se infactivel
                if(check == CONFLICT) {
                    return CONFLICT;
                }
            }
            else {
                for(j=1; j<=inst->nElRow[restr]; j++) { // loop em todos elementos de restr
                    var= inst->idxRow[restr][j]; // pego uma variavel de i
                    //  printf("\nVAR:  %s- %d",inst->colname[var],var);
                    if(bd->l[var]!=bd->u[var]) { // verifica se ja esta fixada
                        variable=-1;
                        bound = evalueteBound (inst,bd,restr,var,nivel); // chama função evoluate bound
                        variable=cprs_check(inst,cprs,var);
                        /*for(h=0;h<cprs->cpr[cprs->depth].nfix;h++){ // verifica se no nivel atual esta restrição ja foi fixada
                                if(var==cprs->cpr[cprs->depth].idxfix[h]){
                                    variable=h; // se estiver fixa, pego posição que foi fixada
                                    break;
                                }
                        }*/
                        if(variable>=0) { // se ja foi fixada verifica se o valor retornado pelo evaluate bound é o mesmo do fixado
                            if ((bound == FIX_AT_ZERO)||(bound ==FIX_AT_ONE)) {
                                if(cprs->cpr[cprs->depth].valfix[variable]!=bound) { // se for diferente retorno
                                    return CONFLICT;
                                }
                            }
                        }
                        else { // se não foi fixada no nivel atual
                            if ((bound == FIX_AT_ZERO)||(bound ==FIX_AT_ONE)) { // se retornado 0 ou 1 adiciono os valores e cprs
                                cprs_add_result(inst,cprs->depth,var,bound,cprs); //add cprs
                                cprs->tfix++; // incremento o numero de fixações
                                new_implication=true; // coloco que teve novas fixações
                            }

                            else if (bound==CONFLICT) { // se retornado 3, retorno conflito
                                cprs->cpr[cprs->depth].result= bound;
                                return CONFLICT;
                            }
                        }
                    }
                    if ((((double) (clock() - startBTt) / CLOCKS_PER_SEC)  >  tRemainder )) {
                        return bound;
                    }
                }
            }

            if ((((double) (clock() - startBTt) / CLOCKS_PER_SEC)  >  tRemainder )) {
                return bound;
            }
        }
        if(cprs->cpr[cprs->depth].nfix<=0) { // se não encontrou nenhuma implicação do nivel atual
            return NO_IMPLICATION; // retorno NO_IMPLICATION
        }
        if ((((double) (clock() - startBTt) / CLOCKS_PER_SEC)  >  tRemainder )) {
            return bound;
        }
        constraints_clear(inst,constr);//limpa as restrições
        for(i=0; i<cprs->cpr[cprs->depth].nfix; i++) { // se possui implicações atualizo vetores l e u
            int position = cprs->cpr[cprs->depth].idxfix[i]; // pego a variavel fixada
            int value = cprs->cpr[cprs->depth].valfix[i]; //pego o valor da variavel fixada
            bd->l[position]= value; // na posição da variavel fixda atualizo l com valor
            bd->u[position]= value;// na posição da variavel fixda atualizo u com valor
            for(j=1; j<=inst->nzCol[position]; j++) { //loop em todas as retrições que a variavel aparece
                r=inst->idxCol[position][j]; //retrição i que a variavel aparece
                if(constr->c_indx[r]<1) {
                    constraint_add(inst,constr,r); // adiciono a restrição
                }
            }
        }
    } while((constr->nconst>0)&&(new_implication!=false)); // enquanto houver retrições e novas implicações acontecerem

    return bound; //retorno o limite
}


//quando falta apenas a variavel fixada pelo BT em todas as restrições que a variavel aparece verifico se posso fixar ou se a fixação gera conflito
int percorre(const Instance *inst, Bounds * bd, int var, int value, int nivel) {
    int j,r,h;
    int bound;

    for(j=1; j<=inst->nzCol[var]; j++) { //percorro as retrições em que a variavel aparece
        r=inst->idxCol[var][j]; // pego a retrição
        bound= evalueteBound (inst,bd,r,var,nivel); // verifico em qual valor fixar ou se ja existe um conflito
        if((bound==FIX_AT_ZERO)||(bound==FIX_AT_ONE)) {
            if(bound!=value) { //se o valor fixado for diferente do valor retornado pelo evalueteBound
                //  printf("\nposition:  %s- value: %d Restr: %s",inst->colname[var],value,inst->Rowname[r]);
                for(h=1; h<= inst->nElRow[r]; h++) {
                    //      printf("\n%s: %d  ",inst->colname[inst->idxRow[r][h]],bd->l[h]);
                }
                //  getchar();
                return CONFLICT; //retorno conflito
            }
        }
    }
    return NO_IMPLICATION; // se OK, retorno NO_IMPLICATION
}

//Libera estrutura cpResultStack
void cprs_free(const Instance *inst,CPResultStack * cprs) {

    int i=0;
    for( i=0; i<cprs->depth_max; i++) {
        free(cprs->cpr[i].idxfix);
        free(cprs->cpr[i].valfix);
    }
    free(cprs->cpr);
}
//Libera estrutura Bounds
void bounds_free(const Instance *inst, Bounds * bd) {

    free(bd->l);
    free(bd->u);
}
//desaloca a esturutura Constraints
void contraints_free(Constraints  * constr) {

    free(constr->c);
    free(constr->c_indx);
}
//desaloca a esturutura Conflicts
void conflicts_free(const Instance *inst,Conflicts * c) {

    int i=0;
    for( i=1; i<=inst->nCols; i++) {
        free(c->confl[i]);
    }
    free(c->confl);
    free(c->ncmax);
    free(c->nconf);
}


