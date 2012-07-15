#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <netcdf.h>
#include "global.h"

typedef struct {
  int pos;
  double complex eigenvalue;
} EIGENVALUE_STRUCT;

int cmpeigenvalueg(EIGENVALUE_STRUCT *eigen1, EIGENVALUE_STRUCT *eigen2);

int wnetcdf(PARAMS_STRUCT *params, GRID_STRUCT *grid,
	    ROTATION_STRUCT *rotation, OUTPUT_CONTROL *output_control,
	    ARPACK_CONTROL *arpack_params, RESULTS_STRUCT *results) {

  int ncid; //NetCDF ID

  //Dimension IDs
  int dim_r_id; //r-dimension ID
  int dim_index_id; //index-dimension ID
  int dim_polar_id; //polar-dimension ID
  int dim_complex_id; //complex-dimension ID
  int dim_ids[3]; 

  //Variable IDs for each of the components
  int br_id;
  int bt_id;
  int bz_id;
  int vr_id;
  int vt_id;
  int vz_id;
  int pi_id;
  int r_id;
  int omega_id;
  int lambda_id;
  int residual_id;

  int ncerror;
  char output_file_path[256];
  size_t put_index[3];
  double sigma_r, sigma_i;
  double complex bz;
  double complex vz;
  int pos;

  double writeable_lambda[results->nconv][2];
  double writeable_residual[results->nconv];
  double writeable_br[results->nconv][grid->numcells][2];
  double writeable_bt[results->nconv][grid->numcells][2];
  double writeable_bz[results->nconv][grid->numcells][2];
  double writeable_vr[results->nconv][grid->numcells][2];
  double writeable_vt[results->nconv][grid->numcells][2];
  double writeable_vz[results->nconv][grid->numcells][2];
  double writeable_pi[results->nconv][grid->numcells][2];
  double writeable_r[grid->numcells];
  double writeable_omega[grid->numcells];

  size_t start[3] = {0, 0, 0};
  size_t count[3];


  //Create local pointers to these guys
  double complex *lambda = results->lambda;
  double complex *eigvecs = results->z;
  double *residual = results->residual;

  //First, lets sort these things so that we only save the eigenvectors
  //that actually have data.

  EIGENVALUE_STRUCT eigenvalues[5*grid->numcells];

  //Sort the eigenvalues
  for (int i = 0; i < results->nconv; i++) {
    eigenvalues[i].pos = i;
    eigenvalues[i].eigenvalue = lambda[i];
  }
  qsort(eigenvalues, results->nconv, sizeof(EIGENVALUE_STRUCT),
        cmpeigenvalueg);
  
  strcpy(output_file_path, output_control->basefilename);
  strcat(output_file_path, ".nc");

  

  ncerror = nc_create(output_file_path, NC_CLOBBER, &ncid);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to create netCDF file.\n");
    return 1;
  }

  ncerror = nc_def_dim(ncid, "index", results->nconv, &dim_index_id);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to create netCDF dimension 'index'.\n");
    return 1;
  }

  ncerror = nc_def_dim(ncid, "r", grid->numcells, &dim_r_id);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to create netCDF dimension 'r'.\n");
    return 1;
  }

  ncerror = nc_def_dim(ncid, "polar", 2, &dim_polar_id);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to create netCDF dimension 'polar'.\n");
    return 1;
  }

  ncerror = nc_def_dim(ncid, "complex", 2, &dim_complex_id);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to create netCDF dimension 'complex'.\n");
    return 1;
  }

  dim_ids[0] = dim_index_id;
  dim_ids[1] = dim_complex_id;

  ncerror = nc_def_var(ncid, "lambda", NC_FLOAT, 2, dim_ids, &lambda_id);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to create netCDF variable 'lambda'.\n");
    return 1;
  }

  ncerror = nc_def_var(ncid, "residual", NC_FLOAT, 1, dim_ids, &residual_id);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to create netCDF variable 'residual'.\n");
    return 1;
  }

  dim_ids[0] = dim_index_id;
  dim_ids[1] = dim_r_id;
  dim_ids[2] = dim_polar_id;

  ncerror = nc_def_var(ncid, "br", NC_FLOAT, 3, dim_ids, &br_id);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to create netCDF variable 'br'.\n");
    return 1;
  }

  ncerror = nc_def_var(ncid, "bt", NC_FLOAT, 3, dim_ids, &bt_id);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to create netCDF variable 'bt'.\n");
    return 1;
  }
  
  ncerror = nc_def_var(ncid, "bz", NC_FLOAT, 3, dim_ids, &bz_id);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to create netCDF variable 'bz'.\n");
    return 1;
  }
  
  ncerror = nc_def_var(ncid, "vr", NC_FLOAT, 3, dim_ids, &vr_id);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to create netCDF variable 'vr'.\n");
    return 1;
  }

  ncerror = nc_def_var(ncid, "vt", NC_FLOAT, 3, dim_ids, &vt_id);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to create netCDF variable 'vt'.\n");
    return 1;
  }

  ncerror = nc_def_var(ncid, "vz", NC_FLOAT, 3, dim_ids, &vz_id);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to create netCDF variable 'vz'.\n");
    return 1;
  }

  ncerror = nc_def_var(ncid, "pi", NC_FLOAT, 3, dim_ids, &pi_id);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to create netCDF variable 'pi'.\n");
    return 1;
  }

  dim_ids[0] = dim_r_id;

  ncerror = nc_def_var(ncid, "r", NC_FLOAT, 1, dim_ids, &r_id);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to create netCDF variable 'r'.\n");
    return 1;
  }

  ncerror = nc_def_var(ncid, "omega", NC_FLOAT, 1, dim_ids, &omega_id);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to create netCDF variable 'omega'.\n");
    return 1;
  }


  //Define these global attributes.

  ncerror = nc_put_att_double(ncid, NC_GLOBAL, "eta",
			      NC_FLOAT, 1, &params->eta);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to write attribute 'eta'.\n");
    return 1;
  }  

  ncerror = nc_put_att_double(ncid, NC_GLOBAL, "nu",
			      NC_FLOAT, 1, &params->nu);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to write attribute 'nu'.\n");
    return 1;
  }

  ncerror = nc_put_att_double(ncid, NC_GLOBAL, "rho",
			      NC_FLOAT, 1, &params->rho);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to write attribute 'rho'.\n");
    return 1;
  }

  ncerror = nc_put_att_double(ncid, NC_GLOBAL, "va",
			      NC_FLOAT, 1, &params->va);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to write attribute 'va'.\n");
    return 1;
  }
  
  ncerror = nc_put_att_double(ncid, NC_GLOBAL, "B0",
			      NC_FLOAT, 1, &params->B0);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to write attribute 'B0'.\n");
    return 1;
  }  

  ncerror = nc_put_att_double(ncid, NC_GLOBAL, "r1",
			      NC_FLOAT, 1, &params->r1);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to write attribute 'r1'.\n");
    return 1;
  }  

  ncerror = nc_put_att_double(ncid, NC_GLOBAL, "r2",
			      NC_FLOAT, 1, &params->r2);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to write attribute 'r2'.\n");
    return 1;
  }  

  ncerror = nc_put_att_double(ncid, NC_GLOBAL, "height",
			      NC_FLOAT, 1, &params->height);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to write attribute 'height'.\n");
    return 1;
  }  

  ncerror = nc_put_att_double(ncid, NC_GLOBAL, "nmode",
			      NC_FLOAT, 1, &params->nmode);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to write attribute 'nmode'.\n");
    return 1;
  }

  ncerror = nc_put_att_double(ncid, NC_GLOBAL, "kz",
			      NC_FLOAT, 1, &params->k);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to write attribute 'kz'.\n");
    return 1;
  }  

  ncerror = nc_put_att_double(ncid, NC_GLOBAL, "m",
			      NC_FLOAT,1, &params->m);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to write attribute 'm'.\n");
    return 1;
  }

  ncerror = nc_put_att_int(ncid, NC_GLOBAL, "numcells",
			   NC_INT, 1, &params->numcells);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to write attribute 'numcells'.\n");
    return 1;
  }

  ncerror = nc_put_att_double(ncid, NC_GLOBAL, "Pm",
			      NC_FLOAT, 1, &params->Pm);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to write attribute 'Pm'.\n");
    return 1;
  }

  ncerror = nc_put_att_double(ncid, NC_GLOBAL, "Re",
			      NC_FLOAT, 1, &params->Re);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to write attribute 'Re'.\n");
    return 1;
  }

  ncerror = nc_put_att_double(ncid, NC_GLOBAL, "Rm",
			      NC_FLOAT, 1, &params->Rm);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to write attribute 'Rm'.\n");
    return 1;
  }

  ncerror = nc_put_att_double(ncid, NC_GLOBAL, "Ha",
			      NC_FLOAT, 1, &params->Ha);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to write attribute 'Ha'.\n");
    return 1;
  }

  ncerror = nc_put_att_int(ncid, NC_GLOBAL, "arpack_modes_requested",
			   NC_INT, 1, &arpack_params->nummodes);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to write attribute 'arpack_modes_requested'.\n");
    return 1;
  }

  ncerror = nc_put_att_int(ncid, NC_GLOBAL, "arpack_modes_converged",
			   NC_INT, 1, &results->nconv);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to write attribute 'arpack_modes_converged'.\n");
    return 1;
  }

  ncerror = nc_put_att_int(ncid, NC_GLOBAL, "arpack_maxiters",
			   NC_INT, 1, &arpack_params->maxiters);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to write attribute 'arpack_maxiters'.\n");
    return 1;
  }

  ncerror = nc_put_att_int(ncid, NC_GLOBAL, "arpack_itersused",
			   NC_INT, 1, &results->itersused);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to write attribute 'arpack_iterused'.\n");
    return 1;
  }

  ncerror = nc_put_att_double(ncid, NC_GLOBAL, "arpack_tol",
			      NC_FLOAT, 1, &arpack_params->tol);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to write attribute 'arpack_tol'.\n");
    return 1;
  }

  sigma_r = creal(arpack_params->sigma);
  ncerror = nc_put_att_double(ncid, NC_GLOBAL, "arpack_sigma_r",
			      NC_FLOAT, 1, &sigma_r);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to write attribute 'arpack_sigma_r'.\n");
    return 1;
  }

  sigma_i = cimag(arpack_params->sigma);
  ncerror = nc_put_att_double(ncid, NC_GLOBAL, "arpack_sigma_i",
			      NC_FLOAT, 1, &sigma_i);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to write attribute 'arpack_sigma_i'.\n");
    return 1;
  }


  //Write units variable attributes
  ncerror = nc_put_att_text(ncid, r_id, "units", 2, "cm");
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to write attribute for 'r' units.\n");
    return 1;
  }

  ncerror = nc_put_att_text(ncid, lambda_id, "units", 3, "1/s");
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to write attribute for 'lambda' units.\n");
    return 1;
  }

  ncerror = nc_put_att_text(ncid, omega_id, "units", 3, "1/s");
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to write attribute for 'omega' units.\n");
    return 1;
  }


  ncerror = nc_put_att_text(ncid, NC_GLOBAL, "title", 46,
			    "From nonaxisymmetric MHD global stability code");
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to write attribute for 'omega' units.\n");
    return 1;
  }

  // A little extra effort to define the magnetic BC attribute

  if(params->magnetic_bc == 1) {
    ncerror = nc_put_att_text(ncid, NC_GLOBAL, "magnetic_bc", 10,
			      "Insulating");
  } else {
    ncerror = nc_put_att_text(ncid, NC_GLOBAL, "magnetic_bc", 10,
			      "Conducting");
  }
  
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to write attribute 'magnetic_bc'.\n");
    return 1;
  }  


  //Leave define mode

  ncerror = nc_enddef(ncid);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to leave define mode.\n");
    return 1;
  }

  
  //Now we just need to fill out all of these fields.

  //Get the eigenvalues (lambda) and residuals in the correct order
  //and load them in the arrays to write.
  
  for(int i=0; i < results->nconv; i++) {    
    writeable_lambda[i][0] = creal(eigenvalues[i].eigenvalue);
    writeable_lambda[i][1] = cimag(eigenvalues[i].eigenvalue);
    writeable_residual[i] = residual[eigenvalues[i].pos];
  }
  
  //Write these values to the file

  count[0] = results->nconv;
  count[1] = 2;
  nc_put_vara_double(ncid, lambda_id, start, count,
		     &writeable_lambda[0][0]);
  nc_put_vara_double(ncid, residual_id, start, count,
		     &writeable_residual[0]);



  //Load up the data arrays with results in the correct order.
   

  for(int i=0; i < results->nconv; i++) {
    pos = eigenvalues[i].pos;
    for (int j = 0; j < grid->numcells; j++) {	
      //br
      writeable_br[i][j][0] = cabs(eigvecs[5*grid->numcells*pos + 5*j + 0]);
      writeable_br[i][j][1] = carg(eigvecs[5*grid->numcells*pos + 5*j + 0]);
      writeable_bt[i][j][0] = cabs(eigvecs[5*grid->numcells*pos + 5*j + 1]);
      writeable_bt[i][j][1] = carg(eigvecs[5*grid->numcells*pos + 5*j + 1]);
      writeable_vr[i][j][0] = cabs(eigvecs[5*grid->numcells*pos + 5*j + 2]);
      writeable_vr[i][j][1] = carg(eigvecs[5*grid->numcells*pos + 5*j + 2]);
      writeable_vt[i][j][0] = cabs(eigvecs[5*grid->numcells*pos + 5*j + 3]);
      writeable_vt[i][j][1] = carg(eigvecs[5*grid->numcells*pos + 5*j + 3]);
      writeable_pi[i][j][0] = cabs(eigvecs[5*grid->numcells*pos + 5*j + 4]);
      writeable_pi[i][j][1] = carg(eigvecs[5*grid->numcells*pos + 5*j + 4]);

      //Now calculate vz and bz
      if (j == 0 || j == grid->numcells-1) {
	vz = 0.0;
	bz = 0.0;
      } else {
	vz = (I*params->m*eigvecs[5*grid->numcells*pos + 5*j + 3] 
	      + eigvecs[5*grid->numcells*pos + 5*j + 2]
	      + (eigvecs[5*grid->numcells*pos + 5*(j+1) + 2] 
		 - eigvecs[5*grid->numcells*pos + 5*(j-1) + 2])/
	      (2.0*grid->dx))/(params->k*grid->r[j]);
	bz = -(I*params->m*eigvecs[5*grid->numcells*pos + 5*j + 1] 
	       + eigvecs[5*grid->numcells*pos + 5*j + 0]
	       + (eigvecs[5*grid->numcells*pos + 5*(j+1) + 0]
		  - eigvecs[5*grid->numcells*pos + 5*(j-1) + 0])/
	       (2.0*grid->dx))/(params->k*grid->r[j]);
      }

      writeable_bz[i][j][0] = cabs(bz);
      writeable_bz[i][j][1] = carg(bz);
      writeable_vz[i][j][0] = cabs(vz);
      writeable_vz[i][j][1] = carg(vz);
    }      
  }


  //Write all of the data arrays to the file
  count[0] = results->nconv;
  count[1] = grid->numcells;
  count[2] = 2;

  ncerror = nc_put_vara_double(ncid, br_id, start, count,
			       (double *)writeable_br);
  if (ncerror != NC_NOERR) {
    if(ncerror == NC_ERANGE) {
      fprintf(stderr, "Value out of range in 'br'\n");
    } else {
      fprintf(stderr, "Failed to write value for 'br'. Error %i.\n",
	      ncerror);
      return 1;
    }
  }  


  ncerror = nc_put_vara_double(ncid, bt_id, start, count,
			       (double *)writeable_bt);
  if (ncerror != NC_NOERR) {
    if(ncerror == NC_ERANGE) {
      fprintf(stderr, "Value out of range in 'bt'\n");
    } else {
      fprintf(stderr, "Failed to write value for 'bt'. Error %i.\n",
	      ncerror);
      return 1;
    }
  }  

  ncerror = nc_put_vara_double(ncid, bz_id, start, count,
			       (double *)writeable_bz);
  if (ncerror != NC_NOERR) {
    if(ncerror == NC_ERANGE) {
      fprintf(stderr, "Value out of range in 'bz'\n");
    } else {
      fprintf(stderr, "Failed to write value for 'bz'. Error %i.\n",
	      ncerror);
      return 1;
    }
  }  


  ncerror = nc_put_vara_double(ncid, vr_id, start, count,
			       (double *)writeable_vr);
  if (ncerror != NC_NOERR) {
    if(ncerror == NC_ERANGE) {
      fprintf(stderr, "Value out of range in 'vr'\n");
    } else {
      fprintf(stderr, "Failed to write value for 'vr'. Error %i.\n",
	      ncerror);
      return 1;
    }
  }  

  ncerror = nc_put_vara_double(ncid, vt_id, start, count,
			       (double *)writeable_vt);
  if (ncerror != NC_NOERR) {
    if(ncerror == NC_ERANGE) {
      fprintf(stderr, "Value out of range in 'vt'\n");
    } else {
      fprintf(stderr, "Failed to write value for 'vt'. Error %i.\n",
	      ncerror);
      return 1;
    }
  }  

  
  ncerror = nc_put_vara_double(ncid, vz_id, (size_t[3]){0,0,0},
			       (size_t[3]){results->nconv, grid->numcells, 2},
			       (double *)writeable_vz);
  if (ncerror != NC_NOERR) {
    if(ncerror == NC_ERANGE) {
      fprintf(stderr, "Value out of range in 'vz'\n");
    } else {
      fprintf(stderr, "Failed to write value for 'vz'. Error %i.\n",
	      ncerror);
      return 1;
    }
  }  


  ncerror = nc_put_vara_double(ncid, pi_id, start, count,
			       (double *)writeable_pi);
  if (ncerror != NC_NOERR) {
    if(ncerror == NC_ERANGE) {
      fprintf(stderr, "Value out of range in 'pi'\n");
    } else {
      fprintf(stderr, "Failed to write value for 'pi'. Error %i.\n",
	      ncerror);
      return 1;
    }
  }  


  //Now write the r and omega values
  //(Finally, something that's already in order and easy to write!)

  count[0] = grid->numcells;
  nc_put_vara_double(ncid, r_id, start, count, grid->r);
  nc_put_vara_double(ncid, omega_id, start, count, rotation->omega);
   

  //Close the file
  ncerror = nc_close(ncid);
  if (ncerror != NC_NOERR) {
    fprintf(stderr, "Failed to close netCDF dataset properly.\n");
    return 1;
  }

return 0;
}


int cmpeigenvalueg(EIGENVALUE_STRUCT *eigen1, EIGENVALUE_STRUCT *eigen2)
{
  if ((creal(eigen1->eigenvalue)) > (creal(eigen2->eigenvalue))) {
    return -1;
  } else if ((creal(eigen1->eigenvalue))
             == (creal(eigen2->eigenvalue))) {
    return 0;
  } else {
    return 1;
  }
}

    
