#ifndef CP_H_INCLUDED
#define CP_H_INCLUDED


void cprs_create(const Instance *inst,CPResultStack * cprs);//constroi, aloca e inicia com 0 a estrutura cpResultStack - que contem o resultado das fixações
void constraints_create(const Instance *inst,Constraints * constr);// estrutura aramazena restrições a serem passadas para o CP
void bounds_create(const Instance *inst, Bounds * bd); //cria estrutura que armazena limites e restrições
void conflicts_create(const Instance *inst,Conflicts * c); //cria estrutura que adiciona pares de conflitos
void cp_clear(const Instance *inst,CPResultStack * cprs,Bounds * bd,Bounds * bdc,Constraints *constr); // zerar as estruturas do CP antes de passar para o BT
void constraints_clear(const Instance *inst, Constraints *constr); //limpa a estrutura Constraints
void cprs_clear(const Instance *inst,CPResultStack * cprs);//limpa a estrutura cpResultStack, antes de passar para a estrutura CP
void constraint_add(const Instance *inst, Constraints *constr, int r);//adiciona uma nova restrição a estrutura Constraints
void constraint_remove(const Instance *inst, Constraints *constr, int r);//remove uma restrição da estrutura Constraints
void bounds_copy(const Instance *inst, Bounds * bd,Bounds *bdc);//copia uma estrutura Bounds para outra Bounds
void conflicts_add(const Instance *inst,Conflicts * c, CPResultStack * cprs);//Adiciona os pares de conflitos encontrados na estrutura Conflicts
void cprs_add_result(const Instance *inst,int pos, int var, int val,CPResultStack * cprs);// Adiciona as fixações realizadas no CP- recebe pos- nivel do CP, var= variavel, val = valor da variavel
int cprs_check(const Instance *inst,CPResultStack * cprs, int var);// verifica se no nivel atual esta restrição ja foi fixada, e se foi se ela foi fixada com o mesmo valor da nova chamada do evalueteBound
double calculatedU(const Instance *inst, Bounds * bd, int r, int var, int nivel);//calucula o valor do limite U
int evalueteBound (const Instance *inst, Bounds * bd, int r, int var, int nivel);// função que verifica se fixo em 1 ou 0 ou conflito ou não sei em que fixar
int calc_lhs(const Instance *inst, Bounds * bd, int r);// verfica se uma restrição gera infactibilidade
int cPropagation (const Instance *inst, Bounds * bd, CPResultStack * cprs, int nivel, Constraints *constr, double tRemainder ,clock_t  startBTt );// Função Constraint Proapagation
int percorre(const Instance *inst, Bounds * bd, int var, int value, int nivel);//quando falta apenas a variavel fixada pelo BT em todas as restrições que a variavel aparece verifico se posso fixar ou se a fixação gera conflito
void cprs_free(const Instance *inst,CPResultStack * cprs);//Libera estrutura cpResultStack
void bounds_free(const Instance *inst, Bounds * bd) ;//Libera estrutura Bounds
void contraints_free(Constraints  * constr);//desaloca a esturutura Constraints
void conflicts_free(const Instance *inst,Conflicts * c);//desaloca a esturutura Connflicts
#endif // CP_H_INCLUDED
