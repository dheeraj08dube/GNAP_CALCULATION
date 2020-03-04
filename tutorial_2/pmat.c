#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include <stddef.h>
#include <time.h>
#include "myhead.h"
#define PI 3.14159265359
#define MAX_LINE_LENGTH 1000
#define D 3			/////////////please also take care of the dimension when switching from 2 to 3 dimension


double periodic (double, double, double);

extern void eigen_calculation (double *a, int N, double *evalues,
			       double *evectors);

extern void dgetrf_ (int *m, int *n, double *A, int *LDA, int *IPIV,
		     int *INFO);
extern void dgetri_ (int *n, double *A, int *LDA, int *IPIV,
		     double *WORK, int *LWORK, int *INFO);
		     
/* DSYEV prototype */
extern void dsyev( char* jobz, char* uplo, int* n, double* a, int* lda,
                double* w, double* work, int* lwork, int* info );
/* Auxiliary routines prototypes */
void print_matrix( char* desc, int m, int n, double* a, int lda , double *evalues);


double
periodic (double a, double b, double box_len)
{
  double deltap;
  deltap = (a - b) - box_len * rint ((a - b) / box_len);
  return (deltap);
}




int
main ()
{
  int zero = 0;
  static int i, j, k, l, count, pi[NP], residue_ID[NP], indicies[NP], p, p1, f, b, count_k, m, no_of_neighood[NP], rows, cols, dummy;
  rows = zero;
  cols = 0;

  static int *pi_of_neighood[NP];

  static double lx, ly, lz, NAP[NP], nap0_radius[NP];
  static double ab1, ab2, ab3;
  int abc;			//for use and throw
  static double Rx[NP], Ry[NP], Rz[NP], rx[NP], ry[NP], rz[NP], r0x, r0y, r0z,
    R0x, R0y, R0z, averageevalues[8], temp[6], romega, dist;
    static double mod_delta[NP];
  for (j = 0; j < 8; j++)
    averageevalues[j] = (double) 0.0;

  char buf[MAX_LINE_LENGTH], a[10][20] ;

  for (i = 0; i < NP; i++)
    {
      Rx[i] = 0.0;
      Ry[i] = 0.0;
      Rz[i] = 0.0;
      rx[i] = 0.0;
      ry[i] = 0.0;
      rz[i] = 0.0;
      mod_delta[i] = 0.0;
    }

  for (i = 0; i < NP; i++)
    {
      no_of_neighood[i] = 0;
    }

  FILE *f1, *f2, *f21, *f22, *f3, *f4, *f15;
  
  f1 = fopen ("xyz.dat", "r");	
  f2 = fopen ("neighbourhood.dat","w");
  f21 = fopen ("for_deri.dat","w");
  f22 = fopen ("population_neigh.dat","w");
  f3 = fopen ("PMATRIX.dat","w");
    


  for (i = 0; i < NP; i++)
    {
      pi[i] = i + 1 ;
      fscanf (f1, " %lf %lf %lf  ",  &Rx[i], &Ry[i], &Rz[i]);	
      //printf("%d %lf %lf %lf \n", pi[i] , Rx[i], Ry[i], Rz[i]);
    }
    
  fscanf (f1,"%lf %lf %lf", &lx, &ly, &lz);
  lx = lx/1.0 ;
  ly = ly/1.0 ;
  lz = lz/1.0 ;
//printf("%lf %lf %lf\n",lx,ly,lz);
  fclose (f1);




  f1 = fopen ("Romega.dat", "r");	////////please make sure that u are attaching the correct .dat file
  fscanf (f1, "%lf", &romega);
  //printf ("%lf\n", romega);
  fclose (f1);


  //printf("\n%d\n", NP);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Capturing the no of neighbour hood particles and their particle ID
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Capturing the no of neighbour hood particles and their particle ID 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Capturing the no of neighbour hood particles and their particle ID 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Capturing the no of neighbour hood particles and their particle ID 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Capturing the no of neighbour hood particles and their particle ID 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Capturing the no of neighbour hood particles and their particle ID 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Capturing the no of neighbour hood particles and their particle ID 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Capturing the no of neighbour hood particles and their particle ID  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Capturing the no of neighbour hood particles and their particle ID 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Capturing the no of neighbour hood particles and their particle ID 
  dist = 0.0;
  double *C_Matrix[NP];
  for (p = 0; p < NP; p++)
    {
    //if ( pi[p] == 1 )
    //if (pi[p] == 1 || pi[p] == 2 ||  pi[p] == 3 ||  pi[p] == 4 ||  pi[p] == 5 ||  pi[p] == 6 ||  pi[p] == 7 ||  pi[p] == 8 )
      //if (pi[p] == 2042 || pi[p] == 1854 || pi[p] == 2043 || pi[p] == 2707 || pi[p] == 1671 || pi[p] == 310 || pi[p] == 1761 || pi[p] == 87 || pi[p] == 309 || pi[p] == 493 || pi[p] == 3981 || pi[p] == 2523 || pi[p] == 3894 || pi[p] == 1760 || pi[p] == 2041 || pi[p] == 3886/* ||residue_ID[p]==25 ||residue_ID[p]==1 ||residue_ID[p]==217 ||residue_ID[p]==178 ||residue_ID[p]==291 */ )	/////////////////please change this while changing the protein
	{
	  R0x = Rx[p];
	  R0y = Ry[p];
	  R0z = Rz[p];

	  for (k = 0; k < NP; k++)
	    {

	      dist = sqrt (pow (periodic (Rx[k], R0x, lx), 2)
			   + pow (periodic (Ry[k], R0y, ly), 2)
			   + pow (periodic (Rz[k], R0z, lz), 2));


	      if (dist <= romega && (dist != 0.000))
		{
		  (no_of_neighood[p])++;
		}
	      dist = 0.0;
	    }
fprintf(f22,"%d\n",(no_of_neighood[p]));
	  pi_of_neighood[p] = calloc ((no_of_neighood[p]), sizeof (int));

	  i = 0;
		  fprintf(f21," %d ", pi[p]);
	  for (k = 0; k < NP; k++)
	    {

	      dist = sqrt (pow (periodic (Rx[k], R0x, lx), 2)
			   + pow (periodic (Ry[k], R0y, ly), 2)
			   + pow (periodic (Rz[k], R0z, lz), 2));


	      if (dist <= romega && (dist != 0.000))
		{
		  pi_of_neighood[p][i] = pi[k];
		  i++;
		  fprintf(f2," %d ", pi_of_neighood[p][i-1]);
		  fprintf(f21," %d ", pi_of_neighood[p][i-1]);
		}
	      dist = 0.0;
	    }
	  rows = (no_of_neighood[p] * D);
	  cols = (D * D);
	  C_Matrix[p] = calloc ((rows * rows), (sizeof (double)));
	  fprintf(f2,"\n\n");
	  fprintf(f21,"\n\n");
	}
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Capturing the no of neighbour hood particles and their particle ID
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Capturing the no of neighbour hood particles and their particle ID
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Capturing the no of neighbour hood particles and their particle ID
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Capturing the no of neighbour hood particles and their particle ID
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Capturing the no of neighbour hood particles and their particle ID
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Capturing the no of neighbour hood particles and their particle ID
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Capturing the no of neighbour hood particles and their particle ID
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Capturing the no of neighbour hood particles and their particle ID
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Capturing the no of neighbour hood particles and their particle ID
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Capturing the no of neighbour hood particles and their particle ID




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////calculate P matrix
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////calculate P matrix
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////calculate P matrix
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////calculate P matrix
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////calculate P matrix
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////calculate P matrix
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////calculate P matrix
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////calculate P matrix
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////calculate P matrix
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////calculate P matrix
  double *P_matrix[NP];

  for (p = 0; p < NP; p++)	//the coordinate p sweeps through various particles in a particular frame  //LOOP for PARTICLES
    {
    //if ( pi[p] == 1 )
    //if (pi[p] == 1 || pi[p] == 2 ||  pi[p] == 3 ||  pi[p] == 4 ||  pi[p] == 5 ||  pi[p] == 6 ||  pi[p] == 7 ||  pi[p] == 8 )
      //if (pi[p] == 2042 || pi[p] == 1854 || pi[p] == 2043 || pi[p] == 2707 || pi[p] == 1671 || pi[p] == 310 || pi[p] == 1761 || pi[p] == 87 || pi[p] == 309 || pi[p] == 493 || pi[p] == 3981 || pi[p] == 2523 || pi[p] == 3894 || pi[p] == 1760 || pi[p] == 2041 || pi[p] == 3886/* ||residue_ID[p]==25 ||residue_ID[p]==1 ||residue_ID[p]==217 ||residue_ID[p]==178 ||residue_ID[p]==291 */ )	/////////////////please change this while changing the protein
	{
	  R0x = Rx[p];
	  R0y = Ry[p];
	  R0z = Rz[p];
	  rows = (no_of_neighood[p] * D);
	  cols = (D * D);
	  double *Rmat = calloc ((rows * cols), (sizeof (double)));


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////calculating the R matrix
	  for (i = 0, m = 0; i < rows, m < no_of_neighood[p]; i++, m++)
	    {
	      b = pi_of_neighood[p][m];
	      *(Rmat + i * cols + 0) = periodic (Rx[b - 1], R0x, lx);
	      *(Rmat + i * cols + 1) = periodic (Ry[b - 1], R0y, ly);
	      *(Rmat + i * cols + 2) = periodic (Rz[b - 1], R0z, lz);
	      i++;
	      *(Rmat + i * cols + 3) = periodic (Rx[b - 1], R0x, lx);
	      *(Rmat + i * cols + 4) = periodic (Ry[b - 1], R0y, ly);
	      *(Rmat + i * cols + 5) = periodic (Rz[b - 1], R0z, lz);
	      i++;
	      *(Rmat + i * cols + 6) = periodic (Rx[b - 1], R0x, lx);
	      *(Rmat + i * cols + 7) = periodic (Ry[b - 1], R0y, ly);
	      *(Rmat + i * cols + 8) = periodic (Rz[b - 1], R0z, lz);
	    }


	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////calculating Rt matrix
	  double *Rmat_transpose = calloc ((cols * rows), (sizeof (double)));
	  for (i = 0; i < cols; i++)
	    {
	      for (j = 0; j < rows; j++)
		{
		  *(Rmat_transpose + i * rows + j) = *(Rmat + j * cols + i);
		}
	    }


	  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////calculating RtR matrix
	  double *Rt_R_mat = calloc ((cols * cols), (sizeof (double)));

	 // double *Rt_R_mat1 = calloc ((cols * cols), (sizeof (double)));

	  for (i = 0; i < cols; i++)
	    for (j = 0; j < cols; j++)
	      for (k = 0; k < rows; k++)
		*(Rt_R_mat + i * cols + j) =
		  *(Rt_R_mat + i * cols + j) +
		  (*(Rmat_transpose + i * rows + k)) *
		  (*(Rmat + k * cols + j));

	/*  for (i = 0; i < cols; i++)
	    for (j = 0; j < cols; j++)
	      for (k = 0; k < rows; k++)
		*(Rt_R_mat1 + i * cols + j) =
		  *(Rt_R_mat1 + i * cols + j) +
		  (*(Rmat_transpose + i * rows + k)) *
		  (*(Rmat + k * cols + j));


	    
	    printf("\n\n\n\n\n\n\n THE MATRIX\n");
	    
	    
	    
	  for (i = 0; i < cols; i++)
	  {
	    for (j = 0; j < cols; j++)
	    {
	    printf("%lf  ", *(Rt_R_mat1 + i * cols + j) );
	    }
	    printf("\n");
	    }
	    
	    */


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////inverting the RtR matrix

	  int N = cols;
	  int NN = cols * cols;
	  int *pivotArray = (int *) malloc (N * sizeof (int));	//since our matrix has three rows
	  int errorHandler;
	  double *lapackWorkspace = (double *) malloc (NN * sizeof (double));
	  dgetrf_ (&N, &N, Rt_R_mat, &N, pivotArray, &errorHandler);
	  //printf ("dgetrf eh, %d, should be zero\n", errorHandler);
	  dgetri_ (&N, Rt_R_mat, &N, pivotArray, lapackWorkspace, &NN,&errorHandler);
	  //printf ("dgetri eh, %d, should be zero\n\n", errorHandler);
	  int index;


/*
	    printf("\n\n\n\n\n\n\n IT'S inverse\n");

	  for (i = 0; i < cols; i++)
	  {
	    for (j = 0; j < cols; j++)
	    {
	    printf("%lf  ", *(Rt_R_mat + i * cols + j) );
	    }
	    printf("\n");
	    }
	    
	    

	    printf("\n\n\n\n\n\n\n \n\n");
	    */

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////calculating the R * RtR-1 matrix
	  double *R_Rt_R_mat = calloc ((rows * cols), (sizeof (double)));

	  for (i = 0; i < rows; i++)
	    for (j = 0; j < cols; j++)
	      for (k = 0; k < cols; k++)
		*(R_Rt_R_mat + i * cols + j) =
		  *(R_Rt_R_mat + i * cols + j) +
		  (*(Rmat + i * cols + k)) * (*(Rt_R_mat + k * cols + j));


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////calculating I-P matrix
	  double *I_P_mat = calloc ((rows * rows), (sizeof (double)));

	  for (i = 0; i < rows; i++)
	    for (j = 0; j < rows; j++)
	      for (k = 0; k < cols; k++)
		*(I_P_mat + i * rows + j) =
		  *(I_P_mat + i * rows + j) +
		  (*(R_Rt_R_mat + i * cols + k)) *
		  (*(Rmat_transpose + k * rows + j));


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Finally calculating the P Matrix
	  if(no_of_neighood[p]<=3)
	  {
	  P_matrix[p] = calloc ((rows * rows), (sizeof (double)));
	  }
	  else
	  {
	  P_matrix[p] = calloc ((rows * rows), (sizeof (double)));
	  for (i = 0; i < rows; i++)
	    for (j = 0; j < rows; j++)
	      *(P_matrix[p] + i * rows + j) =	(-1.0) * *(I_P_mat + i * rows + j);

	  for (i = 0; i < rows; i++)
	  {
	    *(P_matrix[p] + i * rows + i) =  *(P_matrix[p] + i * rows + i) + 1.0;
	  }
	  }



	  for (i = 0; i < rows; i++)
	    for (j = 0; j < rows; j++)
{
		  fprintf(f3," %lf  ", *(P_matrix[p] + i * rows + j) );

}

fprintf(f3,"\n\n");

	  free (Rmat);
	  free (Rmat_transpose);
	  free (Rt_R_mat);
	  free (R_Rt_R_mat);
	  free (I_P_mat);




	}
    }



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////calculate P matrix
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////calculate P matrix
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////calculate P matrix
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////calculate P matrix
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////calculate P matrix
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////calculate P matrix
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////calculate P matrix
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////calculate P matrix
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////calculate P matrix
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////calculate P matrix



}





