#ifndef INSTANCE_H_INCLUDED
#define INSTANCE_H_INCLUDED
typedef enum ContraintType
{
    SetPacking,
    SetPartitioning,
    SetCovering,
    Cardinality,
    OtherConstraint
} ContraintType;

typedef struct {
    int nRows;         /* number of rows (constraints)*/
    int nCols;         /* number of columns (variables)*/
    double *obj;       /* coefficients in the objective function*/
    int *nElRow;       /* number of elements in a row*/
    int **idxRow;       /* index of variables in a row*/
    double **coefRow;  /* coefficients of variables row*/
    char *rowType;     /* type of constraint*/
    ContraintType *consType;
    double *rhs;     /* row lower bounds*/
    int direction; /*maximar ou minimizar*/
    char **colname; /* name of column*/
    char **Rowname; /* name of line*/
    int *nzCol;   /*dizendo em quantas restrições a variável j aparece)*/
    int **idxCol;  /*(dizendo em quais restrições a variável j aparece)*/
    double **coefCol; /* (para cada restrição que j aparece o coeficiente com que ela aparece)*/
    double **ncoefrow; // novos coeficientes de cada linha seguindo a equação Aij Xj<= Bi
    double *newrhs; // novo limite da direita
    int *varsol; //valor da variavel na solução inicial
    int *varwrite; // valor da variavel lida no .txt
    double *relax;

} Instance;


//void loadInstance( Instance *inst,const char* filename, char* filesol);
void instance (Instance *inst,const char* filename);
void isntance_free(Instance *inst);
#endif // INSTANCE_H_INCLUDED
