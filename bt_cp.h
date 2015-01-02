#ifndef BT_CP_H_INCLUDED
#define BT_CP_H_INCLUDED

//estrutura que armazena cada nivel de fixa��es do CP
typedef struct {
    int result; // resultado do EvaluateBound
    int nfix; // numero de fixa��es por nivel
    int *idxfix; // indice das variaveis fixadas
    int *valfix; //valor das variaveis fixadas
} CPResult;

//estrutura para todos os niveis de arrmazenamento de fixa��es do CP
typedef struct {
    CPResult *cpr; // nivel so cpResult
    int depth_max; // numero maximo de niveis utilizados
    int depth; // quantidade de niveis que houveram fixa��es no CP
    int tfix; // numero total de fixa��es
} CPResultStack;

//estrutura que armazena as restri��es a serem passadas para o CP
typedef struct {
    int *c_indx; // indice de onde est�o as restri��es fixdas
    int *c; // restri��es fixadas
    int nconst; // numero de restri��es fixadas
} Constraints;

//estrutura que armazena os limites L e U, as restri��es a serem passadas e os novos coeficientes
typedef struct {
    int *l; //limite inferior L
    int *u; //Limite Superior U
} Bounds;

// estrutura que armazena as variaveis que n�o foram fixadas para utilizar quando o BT pega var de forma aleatoria
typedef struct {
    int *var; //variaveis n�o fixadas
    int *indexv; // indice das posi��es das variaveis n�o fixadas
    int nunfix; // numero de variaveis n�o fixadas
} VUnfixed;

//estrutura para retornar o BT
typedef struct {
    int *j; // variavel
    int nfixC; // numero de fixa��es da restri��o
} UpdatedBounds;

// aramazena os pares de conflitos
typedef struct {
    int idx1;
    int idx2;
    int val1;
    int val2;
} Conflict;

//aramazena todos os pares de conflitos
typedef struct {
    Conflict **confl;
    int *nconf;
    int *ncmax;
} Conflicts;

void time_start(double t); //inicia o tempo
double time_end(); //retorna o tempo de processamento
int startBT(const Instance *inst, Bounds* bd,Bounds * bdc,CPResultStack * cprs,Solution *sol,Solution *solbest,Conflicts * c,Constraints *constr,VUnfixed * vu);
int startBT_Direcionado_Samuel_Combinado(const Instance *inst, Bounds* bd,Bounds * bdc,CPResultStack * cprs,Solution *sol,Solution *solbest,Solution *solstart,Conflicts * c, Constraints *constr,VUnfixed * vu);
void vUfixed_create(const Instance *inst,VUnfixed * vu);//estrutura que armazena vairiavies n�o fixadas
void bt_clear(const Instance *inst,VUnfixed* vu); //limpa as estruturas antes de passar para o BT com 0
void solBetter(const Instance *inst, Bounds * bd,Solution *sol, Solution *solbest); //verica se a solu��o � melhor que a corrente
void copy_sol(Instance *inst, Solution *solbest);//copia solu��o parao BT guiado por solu��o
void remove_var(const Instance *inst,VUnfixed * vu, int var); //remove var da estrutura VUnfixed
void add_var(const Instance *inst,VUnfixed * vu, int var);// adiciona resultdo na estrutura VUnfixed
void b_t_r(const Instance *inst, Bounds* bd,Bounds * bdc,CPResultStack * cprs,int varfix, int value,Solution *sol,Solution *solbest,Conflicts * c, int nivel, Constraints *constr, int control,VUnfixed * vu, int rando); //backtrack
int level(); // nivel de profundidade maximo atigindo pelo BT
int SolIsInfac( const Instance *inst, Bounds* bd,int var);//testa de fixa��o de uma determinda variavel deixa a instancia infactivel
void vUfixed_free(VUnfixed * vu); //libera a estrutura que armazena as variaveis que n�o foram fixadas para utilizar quando o BT pega var de forma aleatoria


#endif // BT_CP_H_INCLUDED
