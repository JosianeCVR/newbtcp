#ifndef LP_HEADER
#define LP_HEADER

#define LP_ME_DEFAULT     0
#define LP_ME_OPTIMALITY  1
#define LP_ME_FEASIBILITY 2

/* Optimization direction */
#define LP_MIN 0
#define LP_MAX 1

/* Optimization result: */
#define LP_OPTIMAL       0
#define LP_INFEASIBLE    1
#define LP_UNBOUNDED     2
#define LP_FEASIBLE      3
#define LP_INTINFEASIBLE 4
#define LP_NO_SOL_FOUND  5
#define LP_ERROR         6

typedef struct _LinearProgram LinearProgram;
typedef LinearProgram * LinearProgramPtr;

/* Model input & output */
void lp_read( LinearProgram *lp, const char *fileName );
void lp_write_lp( LinearProgram *lp, const char *fileName );
void lp_write_sol( LinearProgram *lp, const char *fileName );
void lp_load_mip_start( LinearProgram *lp, int count, const char **colNames, const double *colValues );
void lp_read_mip_start( const char *fileName, char **colNames, double colValues[] );

/* Model creation, modification and destruction */
LinearProgram *lp_create();
LinearProgram *lp_clone( LinearProgram *lp );
void lp_add_row( LinearProgram *lp, const int nz, int *indexes, double *coefs, const char *name, char sense, const double rhs );

/** @brief adds new columns
 *
 *  adds new columns to lp, specifying objective function, bounds, integrality and names
 *
 *  @param count number of columns
 *  @param obj objective function coefficients
 *  @param lb lower bounds - if NULL is specified then it is assumed that all variables have lb=0.0
 *  @param ub upper bounds - if NULL is specified then it is assumed that all variables have ub=infinity
 *  @param integer - vector of boolean values indicating if each variable is integer, if NULL all variables
 *     are assumed to be integral
 *  @param names variable names
 */
void lp_add_col( LinearProgram *lp, double obj, double lb, double ub, char integer, char *name, int nz, int *rowIdx, double *rowCoef );
void lp_add_cols( LinearProgram *lp, const int count, double *obj, double *lb, double *ub, char *integer, char **name );
void lp_add_cols_same_bound( LinearProgram *lp, const int count, double *obj, double lb, double ub, char *integer, char **name );
void lp_add_bin_cols( LinearProgram *lp, const int count, double *obj, char **name );
void lp_free( LinearProgramPtr *lp );
void lp_set_direction( LinearProgram *lp, const char direction );
int lp_get_direction( LinearProgram *lp );
void lp_set_obj( LinearProgram *lp, double *obj );
void lp_chg_obj(LinearProgram *lp, int count, int idx[], double obj[] );
void lp_set_rhs( LinearProgram *lp, int row, double rhs );
void lp_set_col_bounds( LinearProgram *lp, int col, const double lb, const double ub );
void lp_fix_col( LinearProgram *lp, int col, double val );
LinearProgram *lp_pre_process( LinearProgram *lp );

/* Model optimization, results query
   and solution methods parameters */
void lp_set_mip_emphasis( LinearProgram *lp, const int mipEmphasis );
int lp_get_mip_emphasis( LinearProgram *lp );
int lp_optimize( LinearProgram *lp );
double lp_solution_time( LinearProgram *lp );
/* if it is a mip, optimizes as a continuous problem */
int lp_optimize_as_continuous( LinearProgram *lp );
/* primal and dual solution */
double lp_obj_value( LinearProgram *lp );
double *lp_x( LinearProgram *lp );
double *lp_row_price( LinearProgram *lp );
/* multiple solutions (if available) */
int lp_num_saved_sols( LinearProgram *lp );
double lp_saved_sol_obj( LinearProgram *lp, int isol );
double *lp_saved_sol_x( LinearProgram *lp, int isol );

/* command line options */
void lp_parse_options( LinearProgram *lp, int argc, const char **argv );
void lp_help_options( );

/* parameters - input/output */
void lp_set_sol_out_file_name( LinearProgram *lp, const char *sfn );
void lp_set_sol_in_file_name( LinearProgram *lp, const char *sfn );
/* parameters - heuristics */
void lp_set_heur_proximity( LinearProgram *lp, char onOff );
void lp_set_heur_fp_passes( LinearProgram *lp, int passes );
/* parameters - cuts */
void lp_set_cuts( LinearProgram *lp, char onOff );
/* parameters - input/output */
void lp_set_print_messages( LinearProgram *lp, char onOff );
/* parameters - limits */
void lp_set_max_seconds( LinearProgram *lp, int _max );
void lp_set_max_solutions( LinearProgram *lp, int _max );
void lp_set_max_nodes( LinearProgram *lp, int _max );
void lp_set_max_saved_sols( LinearProgram *lp, int _max );
void lp_set_abs_mip_gap( LinearProgram *lp, const double _value );
void lp_set_rel_mip_gap( LinearProgram *lp, const double _value );
/* parameters - parallel */
void lp_set_parallel( LinearProgram *lp, char onOff );


/* Model query */
char lp_is_mip( LinearProgram *lp );
char lp_is_integer( LinearProgram *lp, const int j );
int lp_cols( LinearProgram *lp );
int lp_rows( LinearProgram *lp );
int lp_nz( LinearProgram *lp );
int lp_row( LinearProgram *lp, int row, int *idx, double *coef );
int lp_col( LinearProgram *lp, int col, int *idx, double *coef );
double lp_rhs( LinearProgram *lp, int row );
char lp_sense( LinearProgram *lp, int row );
char *lp_row_name( LinearProgram *lp, int row, char *dest );
char *lp_col_name( LinearProgram *lp, int col, char *dest );
double lp_col_lb( LinearProgram *lp, int col );
double lp_col_ub( LinearProgram *lp, int col );
int lp_col_index( LinearProgram *lp, const char *name );
const double *lp_obj_coef( LinearProgram *lp );

#endif

