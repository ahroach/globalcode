#define PI 3.1415926535897932384626433832795029
#define EPSILON 1e-10
#define MAXSTEPS 1000000

typedef struct {
  double eta; //Magnetic Diffusivity
  double nu; //Viscosity
  double rho; //Density
  double B0; //Applied field
  double va; //Alfven Speed
  double r1; //Inner cylinder radius
  double r2; //Outer cylinder radius
  double omega1; //Rotation frequency of inner cylinder
  double omega2; //Rotation frequency of outer cylinder
  double height; //Effective height of experiment
  double nmode; //Number of half-wavelengths in vertical direction
  double k;
  double kva;
  double m; //Azimuthal mode number
  double Pm;
  double Re;
  double Rm;
  double Ha;
  double zetabar;
  int numcells;
  int magnetic_bc; //0 for perfectly conducting B.C., 1 for insulating
} PARAMS_STRUCT;

typedef struct {
  double *r; //Pointer to array of r-values
  double *x;
  double *r2inv; //Array of 1/r^2 values.
  double *diffuse;
  double *diffuse2;
  double dx;
  int numcells;
  int is;
  int ie;
} GRID_STRUCT;

typedef struct {
  double *omega;
  double a;
  double b;
} ROTATION_STRUCT;

typedef struct {
  int filenum;
  char basefilename[256];
} OUTPUT_CONTROL;


typedef struct {
  double complex *lambda;
  double complex *z;
  double *residual;
  int nconv;
  int itersused;
} RESULTS_STRUCT;


typedef struct compressed_matrix {
  int n;
  int m;
  int kl;
  int ku;
  int lda;
  int ldab;
  double complex *A;
  double complex *B;
  double complex *Bb;
} COMPRESSED_MATRIX;

typedef struct {
  int nummodes;
  double complex sigma;
  double tol;
  int iterate;
  int maxiters;
  char which[5];
} ARPACK_CONTROL;


int get_iparam(char paramname[], char filename[]);
double get_dparam(char paramname[], char filename[]);
void get_sparam(char paramname[], char filename[], char destination[]);
void probgen(char input_file_name[], PARAMS_STRUCT *params);

void arpack_handler(char *input_file_name);
void fullmode_handler(char *input_file_name);
void batchmode_handler(char *input_file_name);

ARPACK_CONTROL *setup_arpack(char *input_file_name);
double complex find_sigma(COMPRESSED_MATRIX *matrix,
                          PARAMS_STRUCT *params, GRID_STRUCT *grid,
                          ROTATION_STRUCT *rotation,
                          ARPACK_CONTROL *arpack_params);
COMPRESSED_MATRIX *create_matrix(int numelems);
GRID_STRUCT *gridgen(PARAMS_STRUCT *params);
ROTATION_STRUCT *couette(PARAMS_STRUCT *params, GRID_STRUCT *grid);
ROTATION_STRUCT *dataprofile(PARAMS_STRUCT *params, GRID_STRUCT *grid);
ROTATION_STRUCT *shearlayer(PARAMS_STRUCT *params, GRID_STRUCT *grid,
			    double shear_width, double shear_radius);
void output(PARAMS_STRUCT *params, GRID_STRUCT *grid,
	    ROTATION_STRUCT *rotation, OUTPUT_CONTROL *output_control,
	    double complex *state, double complex gr);
void wAelem(int i, int j, COMPRESSED_MATRIX *matrix, double complex value);
void pAelem(int i, int j, COMPRESSED_MATRIX *matrix);
void wBelem(int i, int j, COMPRESSED_MATRIX *matrix, double complex value);
RESULTS_STRUCT *eigensolve(COMPRESSED_MATRIX *matrix,
			   PARAMS_STRUCT *params, GRID_STRUCT *grid,
			   ROTATION_STRUCT *rotation,
			   ARPACK_CONTROL *arpack_params);
RESULTS_STRUCT *eigensolve_full(COMPRESSED_MATRIX *matrix,
				PARAMS_STRUCT *params, GRID_STRUCT *grid,
				ROTATION_STRUCT *rotation,
				ARPACK_CONTROL *arpack_params);
int wnetcdf(PARAMS_STRUCT *params, GRID_STRUCT *grid,
            ROTATION_STRUCT *rotation, OUTPUT_CONTROL *output_control,
	    ARPACK_CONTROL *arpack_params, RESULTS_STRUCT *results);
