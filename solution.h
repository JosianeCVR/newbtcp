#ifndef SOLUTION_H_INCLUDED
#define SOLUTION_H_INCLUDED

typedef struct {
    char *vect;        /* vector 0 and 1 - solution*/
    double *lhs;       /* summation constraints left*/
    double infeas;     /* distance fact*/
    double cost;       /*FO*/
    double *inf;
    double *newinf;
} Solution;

void SolutionStarts(Instance *inst, Solution *sol);
void GenerateRandom(const Instance *inst, Solution *sol);
void CalcConstraints(const Instance *inst, Solution *sol);
void Infeas(const Instance *inst, Solution *sol);
void CalculationCost(const Instance *inst, Solution *sol);
void ChangesConstraints(const Instance *inst, Solution *sol, int pos);
void newInfeas(const Instance *inst,Solution *sol, int pos);
void CalculationNewCost(const Instance *inst, Solution *sol, int pos);
void CopySolution(const Instance *inst, Solution *target, const Solution *source );
void flip(const Instance *inst, Solution *sol,int pos);
int SolIsBetter( const Solution *solLS, const Solution *sol);
void liberasol(Solution *sol,const Instance *inst);



#endif // SOLUTION_H_INCLUDED
