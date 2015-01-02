#include <stdio.h>
#include <stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include "instance.h"
#include "lp.h"
#include "solution.h"
#include "bt_cp.h"
#include "cp.h"

#define MAX_DEPTH 20

int nc=0,nr=0;


int main( int argc,char **argv) {


    srand(27);
    char filename[30]; // armazena o nome da instacia
    char insol[30]; // aramzena o nome da solução da instancia
    // char filesol[30]; // armazena o nome do arquivo da solução de saida
    // char inconf[30]; //armazena o nome do arquivo dos conflitos encontrados


    /*variables to main*/
    Instance inst;
    int i,j;
    int n_process=0; // numero de nos processados
    int conflic=0; //numero de conflitos encontrados

    if(argc<2) {
        printf("\nerro");
        exit( EXIT_FAILURE);
    }
    for(i=1; i<argc; i++) {
        if(strcmp(argv[i],"-insol")==0) {
            strcpy(insol,argv[i+1]); // pego o nome que será dado ao arquivo de saida com a solução
        }
        /*if(strcmp(argv[i],"-inconf")==0){
            strcpy(inconf,argv[i+1]);
        }
        if(strcmp(argv[i],"-filesol")==0){
            strcpy(filesol,argv[i+1]); //pego o nome do arquivo .txt lido
        }*/

    }

    //arquivos de saída
    FILE * arqsol;
    FILE * arqdados;
    //FILE * arqconf;
    FILE * fo;
    FILE * infeas;
    FILE * time;
    FILE * nos;
    FILE * confli;
    FILE * nive;


    arqsol=fopen (insol,"a");
    arqdados=fopen ("xdados.txt","a");
    //arqconf=fopen (inconf,"a");

    fo=fopen ("xdfo.txt","a");
    infeas=fopen ("xdinfeas.txt","a");
    time=fopen ("xdtime.txt","a");
    nos=fopen ("xdnos.txt","a");
    confli=fopen ("xdconflic.txt","a");
    nive = fopen ("xnivel.txt","a");

    // loadInstance( &inst,argv[1],filesol);
    loadInstance( &inst,argv[1]);
    strcpy(filename,argv[1]);
    printf("\nfile: %s",filename);
    printf("\ninsol: %s",insol);

    /*allocation of variables to functions*/
    Solution sol, solbest,solstart;
    Bounds bd;
    Bounds bdc;
    CPResultStack  cprs;
    Constraints constr;
    VUnfixed vu;
    Conflicts c;

    SolutionStarts(&inst,&sol); //Inicia estrutura Solution
    SolutionStarts(&inst,&solbest); //Inicia estrutura Solution
    SolutionStarts(&inst,&solstart); //Inicia estrutura Solution
    bounds_create(&inst,&bd); //Inicia estrutura Bounds
    bounds_create(&inst,&bdc);//Inicia estrutura Bounds
    cprs_create(&inst,&cprs);//Inicia estrutura cpResultStack
    constraints_create(&inst,&constr);//Inicia estrutura Constraints
    vUfixed_create(&inst,&vu); //Inicia estrutura de retorno quando as variaveis são passadas aleatorias
    conflicts_create(&inst,&c); //cria a estrutura que armazena pares de conflitos

    fprintf(arqdados, "\nFILE NAME: %s",filename);

    clock_t start = clock();

    printf("\nbababab");

    n_process= startBT(&inst,&bd,&bdc,&cprs,&sol,&solbest,&c,&constr,&vu); //chama BT com valor 1

    // n_process= startBT_Direcionado_Samuel(&inst,&bd,&bdc,&cprs,&sol,&solbest,&c,&constr,&vu);

    //n_process= startBT_Direcionado_Samuel_Combinado(&inst,&bd,&bdc,&cprs,&sol,&solbest,&solstart,&c,&constr,&vu);



    clock_t end = clock();

    double cpuTime = ((double) (end - start)) / CLOCKS_PER_SEC; //calcula o tempo gasto
    printf("\n\ndone in %.3f.\n\n", cpuTime );
    fflush(stdout);

    CalcConstraints(&inst,&solbest);
    Infeas(&inst, &solbest);
    CalculationCost(&inst, &solbest);

    int lev= level();

    printf("\ndistance: %.1f",solbest.infeas);
    printf("\nFO: %.1f",solbest.cost);


    /*fprintf(arqconf,"\nFILE NAME: %s\n\n",filename);
    for(i=1; i<= inst.nCols;i++){
        if( c.nconf[i]>0){
            fprintf(arqconf,"\n%d_________________________\n",i);
            for(j=0; j< c.nconf[i];j++){
               fprintf(arqconf,"\nidx1:%d   val1:%d   idx2:%d   val2:%d   ",c.confl[i][j].idx1,c.confl[i][j].val1,c.confl[i][j].idx2,c.confl[i][j].val2);
            }
        }
    }*/

    for(i=1; i<= inst.nCols; i++) {
        if( c.nconf[i]>0) {
            for(j=0; j< c.nconf[i]; j++) {
                conflic ++;
            }
        }
    }
    printf("\nCONFLICT: %d",conflic);
    fprintf(confli,"\n%d",conflic);
    fprintf(arqsol, "\n__________________________________________________________\n");
    fprintf(arqsol, "\nFILE NAME: %s\n\n",filename);

    for( i=1; i<= inst.nCols; i++) {
        fprintf(arqsol, "rest%d: %s =%i\n",i,inst.colname[i],solbest.vect[i]);
    }


    fprintf(arqdados, "\nFO: %.3f",solbest.cost);
    fprintf(arqdados, "\nInfeas: %.3f",solbest.infeas);
    fprintf(arqdados, "\nNos PRocessados: %d",n_process);
    fprintf(arqdados, "\nTempo: %.3f",cpuTime);
    fprintf(arqdados,"\nConflitos: %d",conflic);
    fprintf(arqdados,"\nNivel: %d",lev);
    fprintf(arqdados, "\n________________________________________________\n");

    fprintf(fo, "%.3f\n",solbest.cost);
    fprintf(infeas, "%.3f\n",solbest.infeas);
    fprintf(time, "%.f\n",cpuTime);
    fprintf(nos, "%d\n",n_process);
    fprintf(nive, "%d\n",lev);


    liberasol(&sol,&inst);
    liberasol(&solbest,&inst);
    liberasol(&solstart,&inst);
    cprs_free(&inst,&cprs);
    bounds_free(&inst,&bd);
    bounds_free(&inst,&bdc);
    contraints_free(&constr);
    conflicts_free(&inst,&c);
    vUfixed_free(&vu);
    isntance_free(&inst);


    return EXIT_SUCCESS;
}
