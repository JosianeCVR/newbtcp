/* at least one solver need to be selected to compile this file */

#ifndef CBC
#ifndef GLPK
#ifndef CPX
#error Compile selecting one solver, i.e. include in GCC parameters: -DCBC or -DGLPK or -DCPX
#endif
#endif
#endif

#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <assert.h>
#include <cfloat>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <string>
#include <vector>
#include <map>
extern "C" {
#include "lp.h"
}
#include <omp.h>

#define ERROR( msg ) \
    fprintf( stderr, msg ); \
    abort(); \
    exit(1);

#define LP_CHECK_COL_INDEX( lp, col )\
    if ((col<0)||(col>=lp_cols(lp))) {\
        fprintf( stderr, "ERROR: Invalid column index: %d\n", col );\
        fprintf( stderr, "\tat %s:%d\n", __FILE__, __LINE__ );\
        abort();\
        exit(EXIT_FAILURE);}

#define LP_CHECK_ROW_INDEX( lp, col )\
    if ((row<0)||(row>=lp_rows(lp))) {\
        fprintf( stderr, "ERROR: Invalid row index: %d\n", row );\
        fprintf( stderr, "\tat %s:%d\n", __FILE__, __LINE__ );\
        abort();\
        exit(EXIT_FAILURE);}


#include <sstream>

#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

using namespace std;

// solver dependent includes and definitions
#ifdef GLPK
#include <glpk.h>
#endif
#ifdef CBC
#include <coin/OsiSolverInterface.hpp>
#include <coin/OsiClpSolverInterface.hpp>
#include <coin/OsiCbcSolverInterface.hpp>
#include <coin/CoinBuild.hpp>
#include <coin/CglPreProcess.hpp>
#endif
#ifdef CPX
#include <cplex.h>
static CPXENVptr LPcpxDefaultEnv = NULL;
#endif

#define EPS 1e-5

#define INT_NOT_SET INT_MIN

struct _LinearProgram {
    int mipEmphasis;

    char optAsContinuous;
    char mipPreProcess;

    double obj;
    int status;

    int nOptimizations;

    vector< double > *_x;
    vector< double > *_pi;
    vector< double > *_obj;

    vector< vector< double > > *_savedSol;
    vector< double > *_savedObj;

    map< string, int > *colNameIdx;

    std::vector< std::pair< std::string, double > > *_mipStart;
    std::vector< int > msIdx;
    std::vector< double > msVal;

    char silent;

    double solutionTime;

    /* parameters */
    int heurFPPasses;
    int heurProx;
    int maxSeconds;
    int maxSolutions; // exit after found n solutions
    int maxSavedSols;
    int maxNodes;
    int cuts;
    int printMessages;
    double absMIPGap;
    double relMIPGap;
    int parallel;

    char solOutFN[256];
    char solInFN[256];

    // defining solver dependent
    // LP storing type
#ifdef GLPK
    glp_prob *_lp;
#endif
#ifdef CBC
    OsiClpSolverInterface *_lp;
    OsiSolverInterface *osiLP;

    // if model is pre processed then this should be stored
    CglPreProcess *cglPP;

#endif
#ifdef CPX
    CPXLPptr cpxLP;
#endif // CPX
};

bool lp_str_is_num( const char *str );

#ifdef GLPK
void lp_config_glpk_params(LinearProgram *lp, glp_iocp *iocp);
#endif
#ifdef CBC
void lp_config_cbc_params(LinearProgram *lp, vector<string> &cbcOpt);
#endif
#ifdef CPX
void lp_check_for_cpx_error(CPXENVptr env, int errorCode, const char *sourceFile, int sourceLine);
void lp_config_cpx_params(LinearProgram *lp);
#endif

void lp_printRootRelaxInfo(LinearProgram *lp);

/* from names */
void lp_gen_mipstart( LinearProgram *lp )
{
    const int n = lp->_mipStart->size();
    lp->msIdx.clear();
    lp->msVal.clear();
    lp->msIdx.reserve( n );
    lp->msVal.reserve( n );
    map< string, int > nameIndex;
    char name[512];

    for (int i=0 ; (i<lp_cols(lp)) ; i++ )
        nameIndex[lp_col_name(lp,i,name)] = i;

    vector< pair< string, double > >::const_iterator mIt;
    for ( mIt=lp->_mipStart->begin() ; (mIt!=lp->_mipStart->end()) ; ++mIt )
    {
        map<string, int>::iterator it = nameIndex.find( mIt->first );
        assert( it != nameIndex.end() );
        lp->msIdx.push_back( it->second );
        lp->msVal.push_back( mIt->second );
    }

}

void lp_initialize(LinearProgram *lp);

LinearProgramPtr lp_create()
{
    LinearProgram *result = (LinearProgramPtr) malloc(sizeof(LinearProgram));
    assert(result);

    lp_initialize(result);

#ifdef GLPK
    result->_lp = glp_create_prob();
    assert(result->_lp);
#endif
#ifdef CBC
    result->_lp   = new OsiClpSolverInterface();
    result->_lp->messageHandler()->setLogLevel(0);
    result->osiLP = dynamic_cast<OsiSolverInterface *>(result->_lp);
    result->cglPP = NULL;

    result->osiLP->setIntParam(OsiNameDiscipline, 1);
#endif
#ifdef CPX
    int cpxError = 1;

    if (!LPcpxDefaultEnv) {
        LPcpxDefaultEnv = CPXopenCPLEX(&cpxError);
        if (!LPcpxDefaultEnv) {
            fprintf(stderr, "Error opening CPLEX environment. Quiting.\n");
            abort();
        }
    }

    result->cpxLP = CPXcreateprob(LPcpxDefaultEnv, &cpxError, "mip");
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
#endif // CPX

    return result;
}

char getFileType(const char *fileName)
{
    if ( strstr(fileName, ".mps") || strstr(fileName, ".MPS") || strstr(fileName, ".mps.gz") || strstr(fileName, ".MPS.GZ") )
        return 'M';

    return 'L';
}

void lp_read(LinearProgram *lp, const char *fileName)
{
    assert(lp != NULL);

#ifdef CPX
    int cpxError;
    cpxError = CPXreadcopyprob(LPcpxDefaultEnv, lp->cpxLP, fileName, NULL);
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
#else
    switch (getFileType(fileName)) {
    case 'L':
#ifdef GLPK
        glp_read_lp(lp->_lp, NULL, fileName);
#endif
#ifdef CBC
        lp->osiLP->readLp(fileName);
#endif
        break;
    case 'M':
#ifdef GLPK
        glp_read_mps(lp->_lp, GLP_MPS_FILE, NULL, fileName);
#endif
#ifdef CBC
        lp->osiLP->readMps(fileName);
#endif
        break;
    }
#endif

    lp->colNameIdx->clear();
    char colName[512];

    for (int i = 0 ; (i < lp_cols(lp)) ; ++i)
        (*lp->colNameIdx)[lp_col_name(lp, i, colName)] = i;
}

void lp_write_lp(LinearProgram *lp, const char *fileName)
{
    assert(lp != NULL);

#ifdef CPX
    int cpxError;
    cpxError = CPXwriteprob(LPcpxDefaultEnv, lp->cpxLP, fileName, "LP");
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);

    return;
#endif

    char fileType = getFileType(fileName);

    switch (fileType) {
    case 'L':
#ifdef GLPK
        glp_write_lp(lp->_lp, NULL, fileName);
#endif
#ifdef CBC
        {
            char outFile[256];
            strcpy(outFile, fileName);
            char *s = NULL;
            if ((s = strstr(outFile, ".lp"))) {
                if (s != outFile) // not at the start
                    *s = '\0';
            }
            lp->osiLP->writeLp(outFile);
        }
#endif

        break;
    case 'M':
#ifdef GLPK
        glp_write_mps(lp->_lp, GLP_MPS_FILE, NULL, fileName);
#endif
#ifdef CBC
        lp->osiLP->writeMps(fileName);
#endif

        break;
    }
}

void lp_set_direction(LinearProgram *lp, const char direction)
{
    assert(lp != NULL);

    switch (direction) {
    case LP_MIN:
#ifdef GLPK
        glp_set_obj_dir(lp->_lp, GLP_MIN);
#endif
#ifdef CBC
        lp->osiLP->setObjSense(1.0);
#endif
#ifdef CPX
        CPXchgobjsen(LPcpxDefaultEnv, lp->cpxLP, CPX_MIN);
#endif
        break;
    case LP_MAX:
#ifdef GLPK
        glp_set_obj_dir(lp->_lp, GLP_MAX);
#endif
#ifdef CBC
        lp->osiLP->setObjSense(-1.0);
#endif
#ifdef CPX
        CPXchgobjsen(LPcpxDefaultEnv, lp->cpxLP, CPX_MAX);
#endif
        break;
    default:
        fprintf(stderr, "Unknow optimization direction: %c\n", direction);
        break;
    }
}

int lp_get_direction(LinearProgram *lp)
{
    assert(lp != NULL);

#ifdef GLPK
    if (glp_get_obj_dir(lp->_lp) ==  GLP_MAX)
        return LP_MAX;
    return LP_MIN;
#endif
#ifdef CBC
    if ((fabs(lp->osiLP->getObjSense() - 1.0)) <= EPS)
        return LP_MIN;
    return LP_MAX;
#endif
#ifdef CPX
    switch (CPXgetobjsen(LPcpxDefaultEnv, lp->cpxLP)) {
    case CPX_MIN:
        return LP_MIN;
        break;
    case CPX_MAX:
        return LP_MAX;
        break;
    default:
        fprintf(stderr, "Invalid value for CPXgetobjsen.\n");
        abort();
    }
#endif
}

void lp_add_row(LinearProgram *lp, const int nz,  int *indexes, double *coefs, const char *name, char sense, const double rhs)
{
#ifdef DEBUG
    assert( lp != NULL );
    assert( nz>=0 );
    assert( indexes != NULL );
    assert( coefs != NULL );
    assert( nz <= lp_cols(lp) );

    /* checking indexes */
    for (int i = 0 ; (i < nz) ; ++i) {
        LP_CHECK_COL_INDEX(lp, indexes[i]);
    }
#endif

    sense = toupper(sense);

    char strsense[2];
    sprintf(strsense, "%c", sense);
    if (!strstr("GLE", strsense)) {
        fprintf(stderr, "LP: cannot handle sense %c.\n", sense);
        abort();
    }

#ifdef CPX
    int cpxError;

    int matBeg[] = { 0, nz };

    cpxError = CPXaddrows(LPcpxDefaultEnv, lp->cpxLP, 0, 1, nz, &rhs, &sense, matBeg, indexes, coefs, NULL, (char **) &name);
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);

    return;
#endif // CPX
#ifdef GLPK
    int currRow = lp_rows(lp) + 1;

    glp_add_rows(lp->_lp, 1);
    glp_set_row_name(lp->_lp, currRow, name);
    switch (sense) {
    case 'L' :
        glp_set_row_bnds(lp->_lp, currRow, GLP_UP, 0.0, rhs);
        break;
    case 'G' :
        glp_set_row_bnds(lp->_lp, currRow, GLP_LO, rhs, 0.0);
        break;
    case 'E' :
        glp_set_row_bnds(lp->_lp, currRow, GLP_FX, rhs, 0.0);
        break;
    }

    register int *endIdx = indexes + nz;
    register int *idxPtr;
    for (idxPtr = indexes ; (idxPtr < endIdx) ; idxPtr++)
        (*idxPtr)++;

    glp_set_mat_row(lp->_lp, currRow, nz, indexes - 1, coefs - 1);

    // restoring indexes
    for (idxPtr = indexes ; (idxPtr < endIdx) ; idxPtr++)
        (*idxPtr)--;
#endif
#ifdef CBC
    int currRow = lp_rows(lp);

    double rLB, rUB;
    switch (sense) {
    case 'E':
        rLB = rhs;
        rUB = rhs;
        break;
    case 'L':
        rLB = -COIN_DBL_MAX;
        rUB = rhs;
        break;
    case 'G':
        rLB = rhs;
        rUB = COIN_DBL_MAX;
        break;
    }
    lp->osiLP->addRow(nz, indexes, coefs, rLB, rUB);
    lp->osiLP->setRowName(currRow, name);
#endif
}

void lp_add_cols(LinearProgram *lp, const int count, double *obj, double *lb, double *ub, char *integer, char **name)
{
#ifdef DEBUG
    assert( lp != NULL );
    assert( count >1  );
    assert( name != NULL );

    if ( (lb) && (ub) )
        for ( int i=0 ; (i<count) ; i++ )
            assert(lb[i] <= ub[i]);
#endif

#ifdef CPX
    int cpxError;

    vector< char > type(count, CPX_CONTINUOUS);

    double *_lb = lb;
    double *_ub = ub;

    vector< double > rLB;
    vector< double > rUB;
    if (!_lb) {
        rLB.resize(count, 0.0);
        _lb = &rLB[0];
    }
    if (!_ub) {
        rUB.resize(count, CPX_INFBOUND);
        _ub = &rUB[0];
    }

    int nIntegers = 0;

    for (int i = 0 ; (i < count) ; i++) {
        if (_ub[i] > CPX_INFBOUND)
            _ub[i] = CPX_INFBOUND;
    }


    if (integer) {
        for (int i = 0 ; (i < count) ; ++i)
            if (integer[i]) {
                nIntegers++;
                _lb[i] = floor(_lb[i] + 0.5);

                if (_ub[i] != CPX_INFBOUND)
                    _ub[i] = floor(_ub[i] + 0.5);

                if ((fabs(_lb[i]) < EPS) && (fabs(_ub[i] - 1.0) < EPS)) {
                    //printf( "var %s binary lb %g  ub %g\n", name[i], _lb[i], _ub[i] );
                    type[i] = CPX_BINARY;
                }
                else {
                    //printf( "var %s integer lb %g  ub %g %g\n", name[i], _lb[i], _ub[i], CPX_INFBOUND );
                    type[i] = CPX_INTEGER;
                }
            }
    }

    /*
    if ((nIntegers) &&
            (CPXgetprobtype(LPcpxDefaultEnv, lp->cpxLP) == CPXPROB_LP)) {
        int cpxError = CPXchgprobtype(LPcpxDefaultEnv, lp->cpxLP, CPXPROB_MILP);
        lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
    } */

    cpxError = CPXnewcols(LPcpxDefaultEnv, lp->cpxLP, count, obj, _lb, _ub, &type[0] , name);
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
#endif // CPX
#ifdef GLPK
    register int j, cols, currCol;
    cols = lp_cols(lp);

    glp_add_cols(lp->_lp, count);

    if (name)
        for (currCol = cols + 1, j = 0 ; (j < count) ; j++, currCol++)
            glp_set_col_name(lp->_lp, currCol, name[j]);

    if (obj)
        for (currCol = cols + 1, j = 0 ; (j < count) ; j++, currCol++)
            glp_set_obj_coef(lp->_lp, currCol, obj[j]);

    if (integer)
        for (currCol = cols + 1, j = 0 ; (j < count) ; j++, currCol++)
            if (integer[j]) {
                if ((fabs(ub[j] - 1.0) <= EPS))
                    glp_set_col_kind(lp->_lp, currCol, GLP_BV);
                else
                    glp_set_col_kind(lp->_lp, currCol, GLP_IV);
            }

    for (currCol = cols + 1, j = 0 ; (j < count) ; j++, currCol++) {
        if (((ub) && (lb))) {
            if ((lb[j] != -DBL_MAX) && (ub[j] != DBL_MAX))
                glp_set_col_bnds(lp->_lp, currCol, GLP_DB, lb[j] , ub[j]);
            else if ((ub[j] == DBL_MAX) && (lb[j] != -DBL_MAX))
                glp_set_col_bnds(lp->_lp, currCol, GLP_LO, lb[j] , ub[j]);
            else if ((ub[j] != DBL_MAX) && (lb[j] == -DBL_MAX))
                glp_set_col_bnds(lp->_lp, currCol, GLP_UP, lb[j] , ub[j]);
            else if ((ub[j] == DBL_MAX) && (lb[j] == -DBL_MAX))
                glp_set_col_bnds(lp->_lp, currCol, GLP_FR, -DBL_MAX , DBL_MAX);
        }   // both LB and UB informed
        else {
            if ((!lb) && (!ub)) {
                glp_set_col_bnds(lp->_lp, currCol, GLP_LO, 0.0 , DBL_MAX);
            }  // no LB and UB informed
            else {
                if (ub) {
                    if (ub[j] == DBL_MAX)
                        glp_set_col_bnds(lp->_lp, currCol, GLP_LO, 0.0 , DBL_MAX);
                    else
                        glp_set_col_bnds(lp->_lp, currCol, GLP_DB, 0.0 , ub[j]);
                }
                else if (lb) {
                    if (lb[j] == -DBL_MAX)
                        glp_set_col_bnds(lp->_lp, currCol, GLP_FR, -DBL_MAX , DBL_MAX);
                    else
                        glp_set_col_bnds(lp->_lp, currCol, GLP_DB, lb[j] , DBL_MAX);
                }
            }  // only LB or UB informed
        }
    } // all columns

    // configuring positive, continuous variables
    // (for binary vars it is not necessary to inform
    // bounds)
    if ((integer) && (!ub))
        for (currCol = cols + 1, j = 0 ; (j < count) ; j++, currCol++)
            if (!integer[j])
                glp_set_col_bnds(lp->_lp, currCol, GLP_LO, 0.0 , DBL_MAX);
#endif
#ifdef CBC
    int cols = lp->osiLP->getNumCols();

    {
        vector< int > starts(count + 1, 0);
        lp->osiLP->addCols(count, &starts[0], NULL, NULL, lb, ub, obj);
    }


    if (integer) {
        for (int j = 0 ; (j < count) ; j++)
            if (integer[j])
                lp->osiLP->setInteger(cols + j);
    }

    for (int j = 0 ; (j < count) ; j++)
        lp->osiLP->setColName(cols + j, name[j]);
#endif
    //  updating column name index
    {
        int idxCol = lp_cols(lp) - count;
        for (int j = 0 ; (j < count) ; j++)
            (*lp->colNameIdx)[name[j]] = idxCol++;
    }
}

void lp_add_cols_same_bound(LinearProgram *lp, const int count, double *obj, double lb, double ub, char *integer, char **name)
{
    assert(lp != NULL);

    vector< double > vlb(count, lb);
    vector< double > vub(count, ub);
    lp_add_cols(lp, count, obj, &vlb[0], &vub[0], integer, name);
}

void lp_add_bin_cols( LinearProgram *lp, const int count, double *obj, char **name )
{
    vector< double > vlb(count, 0.0 );
    vector< double > vub(count, 1.0 );
    vector< char > integer( count, 1 );
    lp_add_cols( lp, count, obj, &vlb[0], &vub[0], &(integer[0]), name );
}

const double *lp_obj_coef( LinearProgram *lp )
{
    assert(lp);

#ifdef CBC
    OsiClpSolverInterface *osilp = lp->_lp;

    return osilp->getObjCoefficients();
#endif
#ifdef GLPK
    {
        int i;
        glp_prob *glp = lp->_lp;
        lp->_obj->resize( lp_cols(lp) );
        for ( i=0 ; (i<lp_cols(lp)) ; i++ )
            (*lp->_obj)[i] = glp_get_obj_coef( glp, i+1 );
    }
    return &((*lp->_obj)[0]);
#endif
#ifdef CPX
    {
        int cpxError =  CPXgetobj( LPcpxDefaultEnv, lp->_lp, &(lp->_obj[0]), 0, lp_cols(lp)-1 );
        lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
    }
    return &((*lp->_obj)[0]);
#endif

    return NULL;
}

int lp_cols(LinearProgram *lp)
{
    assert(lp != NULL);

#ifdef GLPK
    return glp_get_num_cols(lp->_lp);
#endif
#ifdef CBC
    return lp->osiLP->getNumCols();
#endif
#ifdef CPX
    return CPXgetnumcols(LPcpxDefaultEnv, lp->cpxLP);
#endif
}

int lp_rows(LinearProgram *lp)
{
    assert( lp != NULL );

#ifdef GLPK
    return glp_get_num_rows(lp->_lp);
#endif
#ifdef CBC
    return lp->osiLP->getNumRows();
#endif
#ifdef CPX
    return CPXgetnumrows(LPcpxDefaultEnv, lp->cpxLP);
#endif
}


char lp_is_integer(LinearProgram *lp, const int j)
{
    assert( lp != NULL );

    LP_CHECK_COL_INDEX(lp, j);

#ifdef GLPK
    switch (glp_get_col_kind(lp->_lp, j + 1)) {
    case GLP_IV:
    case GLP_BV:
        return 1;
        break;
    }

    return 0;
#endif
#ifdef CBC
    return lp->osiLP->isInteger(j);
#endif
#ifdef CPX
    if (CPXgetprobtype(LPcpxDefaultEnv, lp->cpxLP) == CPXPROB_LP)
        return 0;

    char colType[2];
    int cpxError = CPXgetctype(LPcpxDefaultEnv, lp->cpxLP, colType, j, j);
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);

    return ((colType[0] == CPX_BINARY) || (colType[0] == CPX_INTEGER));
#endif
}

char lp_isMIP(LinearProgram *lp)
{
    int nCols = lp_cols(lp);
    int j;

    for (j = 0 ; (j < nCols) ; j++)
        if (lp_is_integer(lp, j))
            return 1;

    return 0;
}

int lp_optimize_as_continuous(LinearProgram *lp)
{
    assert(lp != NULL);

    lp->optAsContinuous = 1;
    int status = lp_optimize(lp);
    lp->optAsContinuous = 0;

    return status;
}

int lp_optimize(LinearProgram *lp)
{
    assert(lp != NULL);

    lp->solutionTime = 0.0;

    time_t startT;
    time(&startT);

    int isMIP = 0;

    if (!lp->optAsContinuous)
        isMIP = lp_isMIP(lp);

    lp->_x->resize(lp_cols(lp));
    lp->_obj->resize(lp_cols(lp));
    lp->_pi->resize(lp_rows(lp));
    lp->_savedSol->clear();
    lp->_savedObj->clear();
    lp->obj = DBL_MAX;
    lp->status = LP_ERROR;

    fill(lp->_x->begin(), lp->_x->end(), DBL_MAX);
    fill(lp->_obj->begin(), lp->_obj->end(), DBL_MAX);
    fill(lp->_pi->begin(), lp->_pi->end(), DBL_MAX);

    lp->nOptimizations++;

    /* error handling */
    char errorMsg[256] = "";
    int errorLine = -1;

#ifdef GLPK
    {
        // solving linear relaxation
        int status = 1;

        glp_smcp parm;
        glp_init_smcp(&parm);

        if (lp->silent)
            parm.msg_lev = GLP_MSG_OFF;

        status = glp_simplex(lp->_lp, &parm);

        // checking error type
        if (status) {
            switch (status) {
            case GLP_EBADB :
                sprintf(errorMsg,
                        "Unable to start the search, because the initial basis specified in the problem object is\
    invalid—the number of basic (auxiliary and structural) variables is not the same as the\
   number of rows in the problem object.");
                errorLine = __LINE__;
                break;
            case GLP_ESING :
                sprintf(errorMsg,
                        "Unable to start the search, because the basis matrix corresponding\
    to the initial basis is singular within the working precision.");
                errorLine = __LINE__;
                break;
            case GLP_ECOND :
                sprintf(errorMsg,
                        "Unable to start the search, because the basis matrix corresponding to the initial basis is\
    ill-conditioned, i.e. its condition number is too large.");
                errorLine = __LINE__;
                break;
            case GLP_EBOUND :
                sprintf(errorMsg,
                        "Unable to start the search, because some double-bounded\
    (auxiliary or structural) variables have incorrect bounds.");
                errorLine = __LINE__;
                break;
            case GLP_EFAIL :
                sprintf(errorMsg,
                        "The search was prematurely terminated due to the solver failure.");
                errorLine = __LINE__;
                break;
            case GLP_EOBJLL :
                sprintf(errorMsg,
                        "The search was prematurely terminated, because the objective function being maximized has reached its\
    lower limit and continues decreasing (the dual simplex only).");
                errorLine = __LINE__;
                break;
            case GLP_EOBJUL :
                sprintf(errorMsg,
                        "The search was prematurely terminated, because the objective function being minimized has reached its\
    upper limit and continues increasing (the dual simplex only).");
                errorLine = __LINE__;
                break;
            case GLP_EITLIM :
                sprintf(errorMsg,
                        "The search was prematurely terminated, because the simplex iteration limit has been exceeded.");
                errorLine = __LINE__;
                break;
            case GLP_ETMLIM :
                sprintf(errorMsg,
                        "The search was prematurely terminated, because the time limit has been exceeded.");
                errorLine = __LINE__;
                break;
            case GLP_ENOPFS :
                sprintf(errorMsg,
                        "The LP problem instance has no primal feasible solution (only if the LP presolver is used).");
                errorLine = __LINE__;
                break;
            case GLP_ENODFS :
                sprintf(errorMsg,
                        "The LP problem instance has no dual feasible solution (only if the LP presolver is used).");
                errorLine = __LINE__;
                break;
            }

            goto OPTIMIZATION_ERROR;
        }
        else {
            status = glp_get_status(lp->_lp);

            switch (status) {
            case GLP_OPT:
                if (!isMIP) {
                    for (int i = 0 ; (i < lp_cols(lp)) ; ++i)
                        lp->_x->at(i) = glp_get_col_prim(lp->_lp, i + 1);
                    for (int i = 0 ; (i < lp_rows(lp)) ; ++i)
                        lp->_pi->at(i) = glp_get_row_dual(lp->_lp, i + 1);

                    lp->obj = glp_get_obj_val(lp->_lp);
                    lp->status = LP_OPTIMAL;
                    goto OPTIMIZATION_CONCLUDED;
                }
                break;
            case GLP_INFEAS:
            case GLP_NOFEAS:
                lp->status = LP_INFEASIBLE;
                goto OPTIMIZATION_CONCLUDED;
            case GLP_UNBND:
                goto OPTIMIZATION_CONCLUDED;
                lp->status = LP_UNBOUNDED;
            default:
                sprintf(errorMsg, "\n\nGLPK ERROR CALLING glp_simplex at file %s\n", __FILE__);
                errorLine = __LINE__;
                goto OPTIMIZATION_ERROR;
            }
        }

        if (isMIP) {
            {
                glp_iocp ioParams;
                lp_config_glpk_params(lp, &ioParams);
                status = glp_intopt(lp->_lp, &ioParams);
            }

            switch (status) {
            case GLP_EFAIL :
            case GLP_EROOT :
            case GLP_EBOUND :
                sprintf(errorMsg, "\n\nGLPK ERROR CALLING glp_intopt at file %s\n", __FILE__);
                errorLine = __LINE__;
                goto OPTIMIZATION_ERROR;
            }

            switch (glp_mip_status(lp->_lp)) {
            case GLP_OPT:
                lp->status = LP_OPTIMAL;
                goto GLPK_GET_MIP_SOLUTION;
            case GLP_FEAS:
                lp->status = LP_FEASIBLE;
                goto GLPK_GET_MIP_SOLUTION;
            case GLP_NOFEAS:
                lp->status = LP_INFEASIBLE;
                goto OPTIMIZATION_CONCLUDED;
            default:
                sprintf(errorMsg, "\n\nGLPK ERROR CALLING glp_mip_status at file %s\n", __FILE__);
                errorLine = __LINE__;
                goto OPTIMIZATION_ERROR;
            }

GLPK_GET_MIP_SOLUTION:
            for (int i = 0 ; (i < lp_cols(lp)) ; ++i)
                lp->_x->at(i) = glp_mip_col_val(lp->_lp, i + 1);

            lp->obj = glp_mip_obj_val(lp->_lp);
            goto OPTIMIZATION_CONCLUDED;
        }
    }
#endif
#ifdef CBC
    bool deleteLP = false;
    OsiSolverInterface *linearProgram = NULL;

    {
        lp->status = LP_ERROR;

        if (lp_isMIP(lp)) {
            linearProgram = lp->osiLP->clone();
            deleteLP = true;
        }
        else
            linearProgram = lp->osiLP;

        if (lp_isMIP(lp))
            linearProgram->initialSolve();
        else {
            if (lp->nOptimizations >= 2)
                linearProgram->resolve();
            else
                linearProgram->initialSolve();
        }

        if (linearProgram->isAbandoned()) {
            sprintf(errorMsg, "Linear program isAbandoned()\n");
            errorLine = __LINE__;
            goto CBC_OPTIMIZATION_ERROR;
        }
        if ((linearProgram->isProvenPrimalInfeasible()) || (linearProgram->isProvenDualInfeasible())) {
            lp->status = LP_INFEASIBLE;
            goto CBC_OPTIMIZATION_CONCLUDED;
        }

        if ((!isMIP) && (linearProgram->isProvenOptimal())) {
            memcpy(&((*(lp->_x))[0]) , linearProgram->getColSolution(), sizeof(double)*lp_cols(lp));
            memcpy(&((*(lp->_pi))[0]), linearProgram->getRowPrice(), sizeof(double)*lp_rows(lp));

            lp->obj = linearProgram->getObjValue();

            lp->status = LP_OPTIMAL;
            goto CBC_OPTIMIZATION_CONCLUDED;
        }

        if (isMIP) {
            // Pass to Cbc initialize defaults
            CbcModel modelA(*linearProgram);
            if (lp->_mipStart->size())
                modelA.setMIPStart( *lp->_mipStart );

            CbcModel *model = &modelA;
            {
                // filling options and calling solver
#define STR_OPT_SIZE 256
                vector<string> cbcP;
                char **cbcOptStr = NULL;
                CbcMain0(modelA);
                lp_config_cbc_params(lp, cbcP);
                cbcOptStr = (char **) malloc(sizeof(char *)*cbcP.size());
                assert(cbcOptStr);
                cbcOptStr[0] = (char *)malloc(sizeof(char) * cbcP.size() * STR_OPT_SIZE);
                assert(cbcOptStr[0]);
                for (int i = 1 ; (i < (int)cbcP.size()) ; i++) cbcOptStr[i] = cbcOptStr[i - 1] + STR_OPT_SIZE;
                for (int i = 0 ; (i < (int)cbcP.size()) ; i++) strncpy(cbcOptStr[i], cbcP[i].c_str(), STR_OPT_SIZE);

                CbcMain1(cbcP.size(), (const char **)cbcOptStr, modelA);
#undef STR_OPT_SIZE
                free(cbcOptStr[0]);
                free(cbcOptStr);
                cbcOptStr = NULL;
            }

            lp->status = LP_NO_SOL_FOUND;

            if (model->isAbandoned()) {
                sprintf(errorMsg, "Model isAbandoned()\n");
                errorLine = __LINE__;
                goto CBC_OPTIMIZATION_ERROR;
            }
            if ((model->isProvenInfeasible()) || (model->isProvenDualInfeasible())) {
                lp->status = LP_INFEASIBLE;
                goto CBC_OPTIMIZATION_CONCLUDED;
            }

            if (model->bestSolution())
                lp->status = LP_FEASIBLE;
            if (model->isProvenOptimal())
                lp->status = LP_OPTIMAL;

            if (model->bestSolution()) {
                memcpy(&((*(lp->_x))[0]), model->solver()->getColSolution(), sizeof(double)*lp_cols(lp));

                for (int i = 0; (i < model->numberSavedSolutions()) ; i++) {
                    const double *si = model->savedSolution(i);
                    vector< double > ti;
                    ti.insert(ti.end(), si, si + lp_cols(lp));

                    lp->_savedSol->push_back(ti);
                    lp->_savedObj->push_back(model->savedSolutionObjective(i));
                } // saved solution
            } // best solution

            lp->obj = model->getObjValue();
        }
    }
CBC_OPTIMIZATION_CONCLUDED:
    if (deleteLP) {
        assert(linearProgram);
        delete linearProgram;
    }

    goto OPTIMIZATION_CONCLUDED;
CBC_OPTIMIZATION_ERROR:
    if (deleteLP) {
        assert(linearProgram);
        delete linearProgram;
    }

    goto OPTIMIZATION_ERROR;
#endif
#ifdef CPX
    lp_config_cpx_params(lp);

    if ((isMIP) && (!lp->optAsContinuous)) {
        int cpxError = CPXmipopt(LPcpxDefaultEnv, lp->cpxLP);
        lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
        int solStat = CPXgetstat(LPcpxDefaultEnv, lp->cpxLP);
        bool hasSolution = false;

        lp->status = LP_NO_SOL_FOUND;

        switch (solStat) {
        case CPXMIP_OPTIMAL :
        case CPXMIP_OPTIMAL_TOL :
            hasSolution = true;
            {
                int status  = CPXgetobjval(LPcpxDefaultEnv, lp->cpxLP, &lp->obj);
                if (status) {
                    fprintf(stderr, "Could not get objval. At %s:%d.\n", __FILE__, __LINE__);
                    abort();
                }
            }

            lp->status = LP_OPTIMAL;
            break;
        case CPX_STAT_INFEASIBLE :
        case CPX_STAT_INForUNBD :
        case CPXMIP_INFEASIBLE :
        case CPXMIP_INForUNBD :
            lp->status = LP_INFEASIBLE;
            break;
        case CPX_STAT_UNBOUNDED :
        case CPXMIP_UNBOUNDED :
            lp->status = LP_UNBOUNDED;
            break;
        default: {
            int status  = CPXgetobjval(LPcpxDefaultEnv, lp->cpxLP, &lp->obj);
            if (status)
                lp->status = LP_NO_SOL_FOUND;
            else {
                hasSolution = true;
                lp->status = LP_FEASIBLE;
            }
        }
        break;
        }

        if (hasSolution) {
            //printf("get x for %d cols\n", lp_cols(lp) ); getchar();
            assert( lp->_x->size()>=lp_cols(lp) );
            int cpxError = CPXgetx(LPcpxDefaultEnv, lp->cpxLP, &((*(lp->_x))[0]), 0, lp_cols(lp) - 1);
            lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
        }

        goto OPTIMIZATION_CONCLUDED;
    }
    else {
        int status = 0;
        if (CPXgetprobtype(LPcpxDefaultEnv, lp->cpxLP) == CPXPROB_MILP) {
            int cpxError = CPXchgprobtype(LPcpxDefaultEnv, lp->cpxLP, CPXPROB_LP);
            lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
        }
        int cpxError = CPXlpopt(LPcpxDefaultEnv, lp->cpxLP);
        lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
        int solStat = CPXgetstat(LPcpxDefaultEnv, lp->cpxLP);

        switch (solStat) {
        case CPX_STAT_OPTIMAL : {
            status = CPXgetobjval(LPcpxDefaultEnv, lp->cpxLP, &lp->obj);
            if (status) {
                sprintf(errorMsg, "Could not get objval.");
                errorLine = __LINE__;
                goto OPTIMIZATION_ERROR;
            }

            int cpxError = CPXgetx(LPcpxDefaultEnv, lp->cpxLP, &((*(lp->_x))[0]), 0, lp_cols(lp) - 1);
            lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
            cpxError = CPXgetpi(LPcpxDefaultEnv, lp->cpxLP, &((*(lp->_pi))[0]), 0, lp_rows(lp) - 1);
            lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
        }

        lp->status = LP_OPTIMAL;
        goto OPTIMIZATION_CONCLUDED;
        break;
        case CPX_STAT_INFEASIBLE :
            lp->status = LP_INFEASIBLE;
            goto OPTIMIZATION_CONCLUDED;
            break;
        case CPX_STAT_INForUNBD :
            lp->status = LP_INFEASIBLE;
            goto OPTIMIZATION_CONCLUDED;
            break;
        case CPX_STAT_UNBOUNDED :
            lp->status = LP_UNBOUNDED;
            goto OPTIMIZATION_CONCLUDED;
            break;
        default :
            char statStr[256];
            CPXgetstatstring(LPcpxDefaultEnv, solStat, statStr);
            sprintf(errorMsg, "CPLEX CPXlpopt returned unhandled optimization status %s.\n", statStr);
            errorLine = __LINE__;
            goto OPTIMIZATION_ERROR;
            break;
        }
    }
#endif
OPTIMIZATION_CONCLUDED:


    time_t endT;
    time(&endT);

    lp->solutionTime = difftime( endT, startT );
    if ((strlen(lp->solOutFN)) && (lp_obj_value(lp) != DBL_MAX))
        lp_write_sol(lp, lp->solOutFN);
    return lp->status;

OPTIMIZATION_ERROR:
    fprintf(stderr, "\n\n===--->>> ERROR <<<---===\n");
    fprintf(stderr, "\tAt lp.cpp, line: %d\n", errorLine);
    fprintf(stderr, "\tMessage: %s\n", errorMsg);
    fprintf(stderr, "\tSaving LP in error.lp\n");
    lp_write_lp(lp, "error.lp");
    abort();
    exit(EXIT_FAILURE);   // if the first one does not quits...
}

double lp_obj_value(LinearProgram *lp)
{
    assert(lp != NULL);
    if (lp->nOptimizations == 0) {
        fprintf(stderr, "No optimizations have been made with this model.\n");
        abort();
    }

    return lp->obj;
}

double *lp_row_price(LinearProgram *lp)
{
    assert(lp != NULL);

    if (lp->nOptimizations == 0) {
        fprintf(stderr, "No optimizations have been made with this model.\n");
        abort();
        exit(EXIT_FAILURE);
    }

    if (lp->status != LP_OPTIMAL) {
        fprintf(stderr, "\n\nERROR: no dual solution available.\n At: %s:%d\n\n", __FILE__, __LINE__);
        abort();
        exit(EXIT_FAILURE);
    }

    return &((*(lp->_pi))[0]);
}

double *lp_x(LinearProgram *lp)
{
    assert(lp != NULL);

    if ((lp->status != LP_OPTIMAL) && (lp->status != LP_FEASIBLE)) {
        fprintf(stderr, "\n\nERROR: no solution available.\n At: %s:%d\n\n", __FILE__, __LINE__);
        abort();
        exit(EXIT_FAILURE);
    }

    if (lp->nOptimizations == 0) {
        fprintf(stderr, "No optimizations have been made with this model.\n");
        abort();
    }

    return &((*(lp->_x))[0]);
}

int lp_get_mip_emphasis(LinearProgram *lp)
{
    assert(lp != NULL);

    return lp->mipEmphasis;
}

void lp_set_mip_emphasis(LinearProgram *lp, const int mipEmphasis)
{
    assert(lp != NULL);

    lp->mipEmphasis = mipEmphasis;
}

void lp_printRootRelaxInfo(LinearProgram *lp)
{
    assert(lp != NULL);

    double obj = lp_obj_value(lp);
    printf("Root node linear relaxation info:\n");
    printf("\tObjective value: %g\n", obj);
    fflush(stdout);
}

void lp_free(LinearProgramPtr *lp)
{
    assert(lp != NULL);

#ifdef GLPK
    glp_delete_prob((*lp)->_lp);
#endif
#ifdef CBC
    if ((*lp)->cglPP == NULL) {
        delete(*lp)->_lp;
        (*lp)->osiLP = NULL;
    }
    else {
        delete(*lp)->cglPP;
        (*lp)->_lp = NULL;
        (*lp)->osiLP = NULL;
    }
#endif
#ifdef CPX
    CPXfreeprob(LPcpxDefaultEnv, &((*lp)->cpxLP));
#endif // CPX

    delete(*lp)->_x;
    delete(*lp)->_obj;
    delete(*lp)->_pi;
    delete(*lp)->_mipStart;
    delete(*lp)->_savedSol;
    delete(*lp)->_savedObj;
    delete(*lp)->colNameIdx;

    free(*lp);
    *lp = NULL;
}

int lp_col( LinearProgram *lp, int col, int *idx, double *coef )
{
    LP_CHECK_COL_INDEX( lp, col );

    int result = -INT_MAX;

#ifdef CPX
    int surplus= -INT_MAX;
    int rmatbeg[2] = { -INT_MAX, -INT_MAX };
    int cpxError = CPXgetcols(LPcpxDefaultEnv, lp->cpxLP, &result, &rmatbeg[0], idx, coef, lp_rows(lp) + 1, &surplus, row, row );
    assert( surplus>=0 );
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
#endif
#ifdef GLPK
    result = glp_get_mat_col(lp->_lp, col + 1, idx - 1, coef - 1);
    for ( int i = 0 ; (i<result) ; ++i )
        idx[i]--;
#endif
#ifdef CBC
    const CoinPackedMatrix *cpmCol =  lp->osiLP->getMatrixByCol();
    const int nzCol = cpmCol->getVectorLengths()[col];
    const CoinBigIndex *starts = cpmCol->getVectorStarts();
    const int *ridx = cpmCol->getIndices() + starts[col];
    const double *rcoef = cpmCol->getElements() + starts[col];

    for (int j = 0 ; (j < nzCol) ; ++j) {
        idx[j] = ridx[j];
        coef[j] = rcoef[j];
    }

    result = nzCol;
#endif

    return result;
}

int lp_row(LinearProgram *lp, int row, int *idx, double *coef)
{
    LP_CHECK_ROW_INDEX( lp, row );

    int result = -INT_MAX;

#ifdef CPX
    int surplus= -INT_MAX;
    int rmatbeg[2] = { -INT_MAX, -INT_MAX };
    int cpxError = CPXgetrows(LPcpxDefaultEnv, lp->cpxLP, &result, &rmatbeg[0], idx, coef, lp_cols(lp) + 1, &surplus, row, row);
    assert( surplus>=0 );
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
#endif

#ifdef CBC
    const CoinPackedMatrix *cpmRow =  lp->osiLP->getMatrixByRow();
    const int nzRow = cpmRow->getVectorLengths()[row];
    const CoinBigIndex *starts = cpmRow->getVectorStarts();
    const int *ridx = cpmRow->getIndices() + starts[row];
    const double *rcoef = cpmRow->getElements() + starts[row];
    for (int j = 0 ; (j < nzRow) ; ++j) {
        idx[j] = ridx[j];
        coef[j] = rcoef[j];
    }

    result = nzRow;
#endif
#ifdef GLPK
    result = glp_get_mat_row(lp->_lp, row + 1, idx - 1, coef - 1);
    for (int i = 0 ; (i<result) ; ++i)
        idx[i]--;
#endif

    return result;
}

double lp_rhs(LinearProgram *lp, int row)
{
    assert(lp != NULL);

#ifdef DEBUG
    assert(row >= 0);
    assert(row < lp_rows(lp));
#endif
#ifdef CPX
    double rhs;
    int cpxError = CPXgetrhs(LPcpxDefaultEnv, lp->cpxLP, &rhs, row, row);
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
    return rhs;
#endif
#ifdef CBC
    return lp->osiLP->getRightHandSide()[row];
#endif
#ifdef GLPK
    switch (glp_get_row_type(lp->_lp, row + 1)) {
    case GLP_LO:
        return glp_get_row_lb(lp->_lp, row + 1);
    case GLP_UP:
        return glp_get_row_ub(lp->_lp, row + 1);
    case GLP_FX:
        return glp_get_row_ub(lp->_lp, row + 1);
    default:
        abort();
    }
#endif
}

char lp_sense(LinearProgram *lp, int row)
{
    assert(lp != NULL);

#ifdef DEBUG
    assert(row >= 0);
    assert(row < lp_rows(lp));
#endif
#ifdef CPX
    int cpxError;
    char result;
    cpxError = CPXgetsense(LPcpxDefaultEnv, lp->cpxLP, &result, row, row);
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
    return result;
#endif
#ifdef CBC
    return lp->osiLP->getRowSense()[row];
#endif
#ifdef GLPK
    switch (glp_get_row_type(lp->_lp, row + 1)) {
    case GLP_LO:
        return 'G';
    case GLP_UP:
        return 'L';
    case GLP_FX:
        return 'E';
    default:
        abort();
    }
#endif
}

char *lp_row_name(LinearProgram *lp, int row, char *dest)
{

#ifdef DEBUG
    LP_CHECK_ROW_INDEX( lp, row );
    assert(dest);
#endif
#ifdef CPX
    int surplus = 0;
    int cpxError = CPXgetrowname(LPcpxDefaultEnv, lp->cpxLP, &dest, dest, 256, &surplus, row, row);
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
#endif
#ifdef CBC
    strcpy(dest, lp->osiLP->getRowName(row).c_str());
#endif
#ifdef GLPK
    strcpy(dest, glp_get_row_name(lp->_lp, row + 1));
#endif
    return dest;
}

char *lp_col_name(LinearProgram *lp, int col, char *dest)
{
    assert(lp != NULL);

    LP_CHECK_COL_INDEX(lp, col);
    assert(dest);

#ifdef CPX
    int surplus = 0;
    int cpxError = CPXgetcolname(LPcpxDefaultEnv, lp->cpxLP, &dest, dest, 256, &surplus, col, col);
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
#endif
#ifdef CBC
    strcpy(dest, lp->osiLP->getColName(col).c_str());
#endif
#ifdef GLPK
    strcpy(dest, glp_get_col_name(lp->_lp, col + 1));
#endif

    return dest;
}

double lp_col_lb(LinearProgram *lp, int col)
{
    assert(lp != NULL);
    LP_CHECK_COL_INDEX(lp, col);

#ifdef DEBUG
    assert(col >= 0);
    assert(col < lp_cols(lp));
#endif
#ifdef CPX
    double lb;
    int begin = col;
    int end = col;
    int cpxError = CPXgetlb(LPcpxDefaultEnv, lp->cpxLP, &lb, begin, end);
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);

    return lb;
#endif
#ifdef CBC
    return lp->osiLP->getColLower()[col];
#endif
#ifdef GLPK
    return glp_get_col_lb(lp->_lp, col + 1);
#endif

}

double lp_col_ub(LinearProgram *lp, int col)
{
    assert(lp != NULL);
    LP_CHECK_COL_INDEX(lp, col);

#ifdef DEBUG
    assert(col >= 0);
    assert(col < lp_cols(lp));
#endif
#ifdef CPX
    double ub;
    int begin = col;
    int end = col;
    int cpxError = CPXgetub(LPcpxDefaultEnv, lp->cpxLP, &ub, begin, end);
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
    return ub;
#endif
#ifdef CBC
    return lp->osiLP->getColUpper()[col];
#endif
#ifdef GLPK
    return glp_get_col_ub(lp->_lp, col + 1);
#endif
}

void lp_set_obj(LinearProgram *lp, double *obj)
{
    assert(lp != NULL);

#ifdef CPX
    vector< int > idx(lp_cols(lp), 0);
    for (int i = 0 ; (i < lp_cols(lp)) ; i++)
        idx[i] = i;

    int cpxError = CPXchgobj(LPcpxDefaultEnv, lp->cpxLP, lp_cols(lp), &idx[0], obj);
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
#endif
#ifdef CBC
    return lp->osiLP->setObjective(obj);
#endif
#ifdef GLPK
    for (int i = 0 ; (i < lp_cols(lp)) ; ++i)
        glp_set_obj_coef(lp->_lp, i + 1, obj[i]);
#endif
}

void lp_add_col(LinearProgram *lp, double obj, double lb, double ub, char integer, char *name, int nz, int *rowIdx, double *rowCoef)
{
    assert(lp != NULL);

#ifdef CPX
    char type = CPX_CONTINUOUS;
    if (integer) {
        if ((fabs(lb) < EPS) && (fabs(ub - 1.0) < EPS))
            type = CPX_BINARY;
        else
            type = CPX_INTEGER;
    }

    int matBeg[] = { 0, nz };
    int cpxError = CPXaddcols( LPcpxDefaultEnv, lp->cpxLP, 1, nz, &obj, matBeg, rowIdx, rowCoef, &lb, &ub, &name );
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);

    return;
#endif // CPX
#ifdef CBC
    int starts[] = { 0, nz };
    lp->osiLP->addCols(1, starts, rowIdx, rowCoef, &lb, &ub, &obj);
    if (integer)
        lp->osiLP->setInteger(lp_cols(lp) - 1);
    lp->osiLP->setColName(lp_cols(lp) - 1, name);
#endif
#ifdef GLPK
    int cols;

    glp_add_cols(lp->_lp, 1);
    cols = lp_cols(lp);

    if (name)
        glp_set_col_name(lp->_lp, cols, name);

    glp_set_obj_coef(lp->_lp, cols, obj);

    if (integer) {
        if ((fabs(ub - 1.0) <= EPS) && (fabs(lb) <= EPS))
            glp_set_col_kind(lp->_lp, cols, GLP_BV);
        else
            glp_set_col_kind(lp->_lp, cols, GLP_IV);
    }

    if ((lb != -DBL_MAX) && (ub != DBL_MAX))
        glp_set_col_bnds(lp->_lp, cols, GLP_DB, lb , ub);
    else if ((ub == DBL_MAX) && (lb != -DBL_MAX))
        glp_set_col_bnds(lp->_lp, cols, GLP_LO, lb , ub);
    else if ((ub != DBL_MAX) && (lb == -DBL_MAX))
        glp_set_col_bnds(lp->_lp, cols, GLP_UP, lb , ub);
    else if ((ub == DBL_MAX) && (lb == -DBL_MAX))
        glp_set_col_bnds(lp->_lp, cols, GLP_FR, -DBL_MAX , DBL_MAX);

    for (int i = 0 ; (i < nz) ; ++i)
        rowIdx[i]++;

    glp_set_mat_col(lp->_lp, cols, nz, rowIdx - 1, rowCoef - 1);

    for (int i = 0 ; (i < nz) ; ++i)
        rowIdx[i]--;
#endif
}

void lp_set_rhs(LinearProgram *lp, int row, double rhs)
{
    assert(lp != NULL);

#ifdef DEBUG
    assert(row >= 0);
    assert(row < lp_rows(lp));
#endif

#ifdef CPX
    int cpxError = CPXchgrhs(LPcpxDefaultEnv, lp->cpxLP, 1, &row, &rhs);
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
#else

    char sense = lp_sense(lp, row);

#ifdef CBC
    switch (sense) {
    case 'E':
        lp->osiLP->setRowBounds(row, rhs, rhs);
        break;
    case 'G':
        lp->osiLP->setRowBounds(row, rhs, COIN_DBL_MAX);
        break;
    case 'L':
        lp->osiLP->setRowBounds(row, -COIN_DBL_MAX, rhs);
        break;
    default:
        fprintf(stderr, "Unknow sense: %c!\n", sense);
        abort();
        exit(1);
    }
#endif
#ifdef GLPK
    switch (sense) {
    case 'E':
        glp_set_row_bnds(lp->_lp, row + 1, GLP_FX, rhs, rhs);
        break;
    case 'G':
        glp_set_row_bnds(lp->_lp, row + 1, GLP_LO, rhs, DBL_MAX);
        break;
    case 'L':
        glp_set_row_bnds(lp->_lp, row + 1, GLP_UP, -DBL_MAX, rhs);
        break;
    default :
        fprintf(stderr, "Unknow sense: %c!\n", sense);
        abort();
        exit(1);
    }
#endif
#endif
}

double lp_solution_time(LinearProgram *lp)
{
    assert(lp != NULL);

    return lp->solutionTime;
}

#ifdef CPX
void lp_check_for_cpx_error(CPXENVptr env, int errorCode, const char *sourceFile, int sourceLine)
{
    if (errorCode) {
        char errorStr[256];
        CPXgeterrorstring(LPcpxDefaultEnv, errorCode, errorStr);
        fprintf(stderr, "CPLEX Error: %s\n", errorStr);
        fprintf(stderr, "Inside LP Library - %s:%d\n\n", sourceFile, sourceLine);
        abort();
        exit(EXIT_FAILURE);
    }
}
#endif

void lp_set_cuts(LinearProgram *lp, char onOff)
{
    assert(lp != NULL);
    lp->cuts = onOff;
}

void lp_set_max_seconds(LinearProgram *lp, int _max)
{
    assert(lp != NULL);
    lp->maxSeconds = _max;
}

void lp_set_max_solutions(LinearProgram *lp, int _max)
{
    assert(lp != NULL);
    lp->maxSolutions = _max;
}


void lp_set_max_saved_sols(LinearProgram *lp, int _max)
{
    assert(lp != NULL);
    lp->maxSavedSols = _max;
}

void lp_set_max_nodes(LinearProgram *lp, int _max)
{
    assert(lp != NULL);
    lp->maxNodes = _max;
}

void lp_set_print_messages(LinearProgram *lp, char onOff)
{
    assert(lp != NULL);
    lp->printMessages = onOff;
}

void lp_set_parallel(LinearProgram *lp, char onOff)
{
    assert(lp != NULL);
    lp->parallel = onOff;
}

void lp_set_heur_proximity(LinearProgram *lp, char onOff)
{
    assert(lp != NULL);
    lp->heurProx = onOff;
}

void lp_set_heur_fp_passes(LinearProgram *lp, int passes)
{
    assert(lp != NULL);

    lp->heurFPPasses = passes;
}

void lp_config_cbc_params(LinearProgram *lp, vector<string> &cbcP)
{
    assert(lp != NULL);

    cbcP.push_back("someMIP"); // problem name
    if (lp->cuts != INT_NOT_SET) {
        if (lp->cuts)
            cbcP.push_back("-cuts=on");
        else
            cbcP.push_back("-cuts=off");
    }
    if (lp->printMessages != INT_NOT_SET) {
        if (lp->printMessages)
            cbcP.push_back("-log=1");
        else
            cbcP.push_back("-log=0");
    }

    if ( lp->parallel != INT_NOT_SET )
    {
        if ( lp->parallel )
        {
            int nthreads =1;
            //int nthreads = omp_get_num_procs();
            printf("CBC will use %d threads.\n", nthreads );
            char stropt[256];
            sprintf( stropt, "-threads=%d", nthreads );
            cbcP.push_back( stropt );
        }
        else
            cbcP.push_back( "-threads=0" );
    }

    if (lp->maxSeconds != INT_NOT_SET)
    {
        cbcP.push_back("-timeM=elapsed");
        cbcP.push_back("-seconds=" + SSTR(lp->maxSeconds));
    }
    if (lp->maxSolutions != INT_NOT_SET)
        cbcP.push_back("-maxSol=" + SSTR(lp->maxSolutions));
    if (lp->maxNodes != INT_NOT_SET)
        cbcP.push_back("-maxNodes=" + SSTR(lp->maxNodes));
    if (lp->heurFPPasses  != INT_NOT_SET)
        cbcP.push_back("-passF=" + SSTR(lp->heurFPPasses));
    if (lp->heurProx != INT_NOT_SET) {
        if (lp->heurProx)
            cbcP.push_back("-proximity=on");
        else
            cbcP.push_back("-proximity=off");
    }
    if (lp->maxSavedSols != INT_NOT_SET)
        cbcP.push_back("-maxSaved=" + SSTR(lp->maxSavedSols));

    if (strlen(lp->solInFN))
        cbcP.push_back("-mips=" + string(lp->solInFN));

    cbcP.push_back("-solve");
    cbcP.push_back("-quit");
}

void lp_set_col_bounds(LinearProgram *lp, int col, const double lb, const double ub)
{
    LP_CHECK_COL_INDEX(lp, col);

    double l = lb;
    double u = ub;
    if ( lp_is_integer(lp,col) )
    {
        if (( l!=-DBL_MAX ) && (l!=DBL_MIN))
            l = floor(lb + 0.5);
        if ( u!=DBL_MAX )
            u = floor(ub + 0.5);
    }

#ifdef CBC
    OsiSolverInterface *linearProgram = lp->osiLP;
    linearProgram->setColBounds(col, l, u);
#endif
#ifdef GLPK
    if ( fabs(l-u)<= EPS )
        glp_set_col_bnds(lp->_lp, col+1, GLP_FX, l, l);
    else if ( (l==-DBL_MAX) || (l==DBL_MIN) )
        glp_set_col_bnds( lp->_lp, col+1, GLP_UP, l, u );
    else if ( u==DBL_MAX )
        glp_set_col_bnds( lp->_lp, col+1, GLP_LO, l, u );
    else
        glp_set_col_bnds( lp->_lp, col+1, GLP_DB, l, u );
#endif
#ifdef CPX
    const int idx[] = { col };
    if ( fabs(l-u)<= EPS )
    {
        int cpxError = CPXchgbds( LPcpxDefaultEnv, lp->_lp, 1, idx, 'B', lb ) ;
        lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
    }
    else
    {
        if ( (l!=-DBL_MAX) && (l!=DBL_MIN) )
        {
            int cpxError;
            cpxError = CPXchgbds( LPcpxDefaultEnv, lp->_lp, 1, idx, 'L', lb ) ;
            lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
        }
        if ( u!=DBL_MAX )
        {
            int cpxError = CPXchgbds( LPcpxDefaultEnv, lp->_lp, 1, idx, 'U', ub ) ;
            lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
        }
    }
#endif


}

double lp_saved_sol_obj(LinearProgram *lp, int isol)
{
    assert(isol >= 0);
    assert(isol < (int)lp->_savedObj->size());
    return lp->_savedObj->at(isol);
}

double *lp_saved_sol_x(LinearProgram *lp, int isol)
{
    assert(lp != NULL);

    assert(isol >= 0);
    assert(isol < (int)lp->_savedSol->size());
    return &lp->_savedSol->at(isol)[0];
}

LinearProgram *lp_clone(LinearProgram *lp)
{
    assert(lp != NULL);

    LinearProgram *result = lp_create();

    memcpy(result, lp, sizeof(LinearProgram));

    lp_initialize(result);
    result->_x         = lp->_x;
    result->_obj       = lp->_obj;
    result->_pi        = lp->_pi;
    result->_mipStart  = lp->_mipStart;
    result->_savedSol  = lp->_savedSol;
    result->_savedObj  = lp->_savedObj;
    result->colNameIdx = lp->colNameIdx;



#ifdef CBC
    result->osiLP = lp->osiLP->clone();
    result->_lp = dynamic_cast<OsiClpSolverInterface *>(result->osiLP);
#endif
#ifdef GLPK
    result->_lp = glp_create_prob();
    glp_copy_prob(result->_lp, lp->_lp, GLP_ON);
#endif
#ifdef CPX
    int cpxError;
    result->cpxLP = CPXcreateprob(LPcpxDefaultEnv, &cpxError, "mip");
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
    //cpxError = CPXcopylp( LPcpxDefaultEnv, lp->cpxLP, lp_cols(lp), lp_rows(lp), CPXgetobjsen( lp->cpxLP), , const double * rhs, const char * sense, const int * matbeg, const int * matcnt, const int * matind, const double * matval, const double * lb, const double * ub, const double * rngval)
#endif

    return result;
}

void lp_fix_col(LinearProgram *lp, int col, double val)
{
    LP_CHECK_COL_INDEX(lp, col);

    if (lp_is_integer(lp, col))
        val = floor(val + 0.5);
#ifdef CBC
    lp->osiLP->setColBounds(col, val, val);
#endif
#ifdef CPX
    char lu = 'B';
    int cpxError = CPXchgbds(LPcpxDefaultEnv , lp->cpxLP, 1, &col, &lu, &val);
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
#endif
#ifdef GLPK
    glp_set_col_bnds(lp->_lp, col+1, GLP_FX, val, val);
#endif
}

LinearProgram *lp_pre_process(LinearProgram *lp)
{
    assert(lp != NULL);

#ifdef CBC
    LinearProgram *result = (LinearProgramPtr) malloc(sizeof(LinearProgram));

    lp_initialize(result);

    result->cglPP = new CglPreProcess();
    result->_lp = dynamic_cast<OsiClpSolverInterface *>(result->cglPP->preProcess(*(lp->_lp), false, 4));
    result->osiLP = dynamic_cast<OsiSolverInterface *>(result->_lp);
    result->_lp->setIntParam(OsiNameDiscipline, 1);
    result->_lp->messageHandler()->setLogLevel(0);

    return result;
#else
    fprintf(stderr, "\nERROR: Pre-processing not implemented for other solvers.\n");
    abort();
    exit(EXIT_FAILURE);
#endif
}

void lp_initialize(LinearProgram *lp)
{
    assert(lp != NULL);

    lp->_x = new vector< double >();
    lp->_obj = new vector< double >();
    lp->_pi = new vector< double >();
    lp->_mipStart = new std::vector< std::pair< std::string, double > >();
    lp->_savedSol = new vector< vector<double> >();
    lp->_savedObj = new vector<double>();
    lp->colNameIdx = new map< string, int >();

    lp->optAsContinuous = 0;
    lp->mipPreProcess = 0;
    lp->nOptimizations = 0;
    lp->silent = 0;
    lp->solutionTime = 0.0;
    lp->obj = DBL_MAX;
    lp->status = LP_ERROR;
    lp->absMIPGap = DBL_MAX;
    lp->relMIPGap = DBL_MAX;

    strcpy(lp->solOutFN, "");
    strcpy(lp->solInFN, "");

    /* parameters */
    lp->maxSeconds    = INT_NOT_SET;
    lp->maxSavedSols  = INT_NOT_SET;
    lp->heurFPPasses  = INT_NOT_SET;
    lp->heurProx      = INT_NOT_SET;
    lp->maxNodes      = INT_NOT_SET;
    lp->cuts          = INT_NOT_SET;
    lp->printMessages = INT_NOT_SET;
    lp->maxSolutions  = INT_NOT_SET;
    lp->mipEmphasis   = LP_ME_DEFAULT;
    lp->parallel      = INT_NOT_SET;
}

int lp_col_index(LinearProgram *lp, const char *name)
{
    assert(lp != NULL);

    map< string, int >::const_iterator mIt;
    mIt = lp->colNameIdx->find(string(name));
    if (mIt == lp->colNameIdx->end())
        return -1;

    return mIt->second;
}

#ifdef GLPK
void lp_config_glpk_params(LinearProgram *lp, glp_iocp *iocp)
{
    glp_init_iocp(iocp);
    if (lp->silent)
        iocp->msg_lev = GLP_MSG_OFF;

    if ( lp->cuts == INT_NOT_SET )
        lp->cuts = 1;

    if ( lp->cuts )
    {
        iocp->gmi_cuts  = GLP_ON;
        iocp->mir_cuts  = GLP_ON;
        iocp->cov_cuts  = GLP_ON;
        iocp->clq_cuts  = GLP_ON;
    }

    if ( lp->heurProx != INT_NOT_SET )
        iocp->ps_heur = GLP_ON;
    if ( lp->heurFPPasses != INT_NOT_SET )
        iocp->fp_heur = GLP_ON;

    /*iocp->ps_heur = GLP_ON;*/
    iocp->presolve = GLP_ON;
    iocp->br_tech = GLP_BR_PCH;

    switch (lp->mipEmphasis) {
    case LP_ME_OPTIMALITY:
        iocp->presolve  = GLP_ON;
        iocp->pp_tech   = GLP_PP_ALL;
        iocp->gmi_cuts  = GLP_ON;
        iocp->mir_cuts  = GLP_ON;
        iocp->cov_cuts  = GLP_ON;
        iocp->clq_cuts  = GLP_ON;
        iocp->fp_heur   = GLP_ON;

        break;
    }

    if (lp->maxSeconds != INT_NOT_SET)
        iocp->tm_lim = lp->maxSeconds * 1000.0;

}
#endif

#ifdef CPX
void lp_config_cpx_params(LinearProgram *lp)
{
    if (lp->maxSeconds != INT_NOT_SET)
        CPXsetdblparam(LPcpxDefaultEnv, CPX_PARAM_TILIM, lp->maxSeconds);
    if (lp->maxSolutions == 1)
        CPXsetdblparam(LPcpxDefaultEnv, CPX_PARAM_EPGAP, 1.0);
    if (lp->maxNodes != INT_NOT_SET)
        CPXsetintparam(LPcpxDefaultEnv, CPX_PARAM_NODELIM, lp->maxNodes);
    if ((lp->silent) || (!lp->printMessages))
        CPXsetintparam(LPcpxDefaultEnv, CPX_PARAM_SCRIND, CPX_OFF);
    else
        CPXsetintparam(LPcpxDefaultEnv, CPX_PARAM_SCRIND, CPX_ON);

    if ( lp->absMIPGap != DBL_MAX )
    {
        double agap;
        CPXgetdblparam( LPcpxDefaultEnv, CPXPARAM_MIP_Tolerances_AbsMIPGap, &agap );
        printf("changing absolute MIP GAP from %g to %g\n", agap, lp->absMIPGap );
        CPXsetdblparam( LPcpxDefaultEnv, CPXPARAM_MIP_Tolerances_AbsMIPGap, lp->absMIPGap );
    }

    if ( lp->relMIPGap != DBL_MAX )
    {
        double rgap;
        CPXgetdblparam( LPcpxDefaultEnv, CPXPARAM_MIP_Tolerances_MIPGap, &rgap );
        printf("changing relative MIP GAP from %g to %g\n", rgap, lp->relMIPGap );
        CPXsetdblparam( LPcpxDefaultEnv, CPXPARAM_MIP_Tolerances_MIPGap, lp->relMIPGap );
    }



    if ( (lp->_mipStart) && (lp->_mipStart->size()) )
    {
        printf("setting cpx mips start\n");
        lp_gen_mipstart(lp);
        int cpxError = CPXcopymipstart( LPcpxDefaultEnv  , lp->cpxLP, lp->msIdx.size(), &lp->msIdx[0], &lp->msVal[0]);
        lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
    }
}
#endif // CPX


void lp_write_sol(LinearProgram *lp, const char *fileName)
{
    assert(lp != NULL);

    if (lp_obj_value(lp) == DBL_MAX) {
        fprintf(stderr, "No solution to write.\n");
        abort();
        exit(EXIT_FAILURE);
    }

    FILE *fsol = fopen(fileName, "w");
    if (fsol == NULL) {
        fprintf(stderr, "Could not open file %s.\n", fileName);
        abort();
        exit(EXIT_FAILURE);
    }

    if (lp->status == LP_OPTIMAL)
        fprintf(fsol, "Optimal (within gap tolerance) - objective value %g - solution time %g\n", lp_obj_value(lp), lp->solutionTime );
    else
        fprintf(fsol, "Stopped on iterations - objective value %g - solution time %g\n", lp_obj_value(lp), lp->solutionTime );

    const int nCols = lp_cols(lp);
    const double *x = lp_x(lp);
    char cName[256];
    for (int i = 0 ; (i < nCols) ; ++i) {
        double xv = x[i];

        if (fabs(xv)<EPS)
            continue;

        if (lp_is_integer(lp, i))
            xv = floor(xv + 0.5);

        fprintf(fsol, "%d %s %g\n", i, lp_col_name(lp, i, cName), xv);
    }
    fclose(fsol);
}

void lp_parse_options(LinearProgram *lp, int argc, const char **argv)
{
    for (int i = 0 ; (i < argc) ; ++i) {
        if (argv[i][0] != '-')
            continue;

        char optLower[256];
        strncpy(optLower, argv[i], 256);
        int len = strlen(optLower);
        for (int j = 0 ; (j < len) ; ++j)
            optLower[j] = tolower(optLower[j]);

        if (strstr(optLower, "maxsec")) {
            if (i + 1 == argc) {
                fprintf(stderr, "enter number of seconds.\n");
                exit(EXIT_FAILURE);
            }

            int sec = atoi(argv[i + 1]);

            printf("> setting max seconds to %d\n", sec);

            lp_set_max_seconds(lp, sec);
            continue;
        }
        if (strstr(optLower, "maxsol")) {
            if (i + 1 == argc) {
                fprintf(stderr, "enter number of solutions.\n");
                exit(EXIT_FAILURE);
            }

            int sol = atoi(argv[i + 1]);

            printf("> setting max solutions to %d\n", sol);

            lp_set_max_solutions(lp, sol);
            continue;
        }
        if (strstr(optLower, "absgap")) {
            if (i + 1 == argc) {
                fprintf(stderr, "enter the allowed absolute MIP gap.\n");
                exit(EXIT_FAILURE);
            }

            double agap = atof(argv[i + 1]);

            lp_set_abs_mip_gap( lp, agap );
            continue;
        }
        if (strstr(optLower, "relgap")) {
            if (i + 1 == argc) {
                fprintf(stderr, "enter the relative MIP gap.\n");
                exit(EXIT_FAILURE);
            }

            double rgap = atof(argv[i + 1]);

            lp_set_rel_mip_gap( lp, rgap );
            continue;
        }
        if (strstr(optLower, "soloutfn")) {
            if (i + 1 == argc) {
                fprintf(stderr, "enter solution file name.\n");
                exit(EXIT_FAILURE);
            }

            printf("> setting solution output file name to %s\n", argv[i + 1]);

            lp_set_sol_out_file_name(lp, argv[i + 1]);
            continue;
        }
        if (strstr(optLower, "solinfn")) {
            if (i + 1 == argc) {
                fprintf(stderr, "enter solution file name.\n");
                exit(EXIT_FAILURE);
            }

            printf("> setting solution input file name to %s\n", argv[i + 1]);

            lp_set_sol_in_file_name(lp, argv[i + 1]);
            continue;
        }
    }
}

void lp_set_sol_out_file_name(LinearProgram *lp, const char *sfn)
{
    assert(lp);

    strcpy(lp->solOutFN, sfn);
}

void lp_set_sol_in_file_name(LinearProgram *lp, const char *sfn)
{
    assert(lp);

    strcpy(lp->solInFN, sfn);
}

void lp_load_mip_start(LinearProgram *lp, int count, const char **colNames, const double *colValues)
{
    lp->_mipStart->clear();
    for (int i = 0 ; (i < count) ; i++)
        lp->_mipStart->push_back(pair< string, double >(colNames[i], colValues[i]));
}

void lp_chg_obj(LinearProgram *lp, int count, int idx[], double obj[])
{
#ifdef CPX
    int cpxError = CPXchgobj(LPcpxDefaultEnv, lp->cpxLP, count, idx, obj);
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
#endif
#ifdef CBC
    for (int i=0 ; (i<count) ; i++ )
        lp->_lp->setObjCoeff( idx[i], obj[i]);
#endif
}

int lp_nz(LinearProgram *lp)
{
    assert(lp);
#ifdef CPX
    return CPXgetnumnz(LPcpxDefaultEnv, lp->cpxLP);
#endif
#ifdef CBC
    lp->_lp->getNumElements();
#endif
#ifdef GLPK
    return glp_get_num_nz(lp->_lp);
#endif

}

void lp_help_options()
{
    printf("options:\n");
    printf("\t-maxSec sec        :  specifies timelimit of 'sec' seconds\n");
    printf("\t-maxSol sol        :  search will be ended when 'sol' solutions are found.\n");
    printf("\t-absgap gap        :  set allowed the absolute MIP gap to 'gap'.\n");
    printf("\t-relgap gap        :  set allowed the relative MIP gap to 'gap'.\n");
    printf("\t-solOutFN solfname :  file name to save solution\n");
    printf("\t-solInFN  solfname :  file name to read solution\n");
}

void lp_read_mip_start( const char *fileName, char **colNames, double colValues[] )
{
#define STR_SIZE 256
    FILE *f = fopen( fileName, "r" );
    if (!f)
    {
        fprintf( stderr, "Could not open mipstart from file %s at %s:%d.\n", fileName, __FILE__, __LINE__ );
        abort();
        exit( EXIT_FAILURE );
    }
    char line[STR_SIZE];

    int nLine = 0;
    int nRCols = 0;
    while (fgets( line, STR_SIZE, f ))
    {
        ++nLine;
        char col[4][STR_SIZE];
        int nread = sscanf( line, "%s %s %s %s", col[0], col[1], col[2], col[3] );
        if (!nread)
            continue;
        /* line with variable value */
        if (strlen(col[0])&&isdigit(col[0][0])&&(nread>=3))
        {
            if (!lp_str_is_num(col[0]))
            {
                fprintf( stderr, "LP warning: reading: %s, line %d - first column in mipstart file should be numeric, ignoring.", fileName, nLine );
                continue;
            }
            if (!lp_str_is_num(col[2]))
            {
                fprintf( stderr, "LP warning: reading: %s, line %d - third column in mipstart file should be numeric, ignoring.", fileName, nLine  );
                continue;
            }

            //int idx = atoi( col[0] );
            char *name = col[1];
            double value = atof( col[2] );
            //double obj = 0.0;
//         if (nread >= 4)
//            obj = atof( col[3] );

            strncpy( colNames[nRCols], name, STR_SIZE );
            colValues[nRCols] = value;

            ++nRCols;
        }
    }

    if (nRCols)
        printf("LP: mipstart values read for %d variables.\n", nRCols );
    else
        printf("LP: no mipstart solution read from %s.\n", fileName );

    fclose(f);

    return;
#undef STR_SIZE
}

void lp_set_abs_mip_gap( LinearProgram *lp, const double _value )
{
    lp->absMIPGap = _value;
}

void lp_set_rel_mip_gap( LinearProgram *lp, const double _value )
{
    lp->relMIPGap = _value;
}

bool lp_str_is_num( const char *str )
{
    const size_t l = strlen(str);

    for ( size_t i=0 ; i<l ; ++i )
        if (!(isdigit(str[i])||(str[i]=='.')))
            return false;

    return true;
}

