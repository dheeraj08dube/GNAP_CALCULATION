/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	Copyright (c) 2011-2014 The plumed team
	(see the PEOPLE file at the root of the distribution for a list of names)

	See http://www.plumed-code.org for more information.

	This file is part of plumed, version 2.

	plumed is free software: you can redistribute it and/or modify
	it under the terms of the GNU Lesser General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	plumed is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public License
	along with plumed.  If not, see <http://www.gnu.org/licenses/>.
	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Colvar.h"
#include "core/PlumedMain.h"
#include "ActionRegister.h"
#include "tools/PDB.h"
#include "core/Atoms.h"
#include <iostream>
#include <new>
#include <memory>
#include <cstring>
#include <string>
#include <array>
#include <iomanip>
#include <cmath>
#include <string>
#include <cstdio>
#include <string>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include <stddef.h>
#include <time.h>
//#include "myhead.h"
#include <fstream>
#include "Colvar.h"
#include "ActionRegister.h"
#include<stdio.h>
#include <string>
#include <cmath>
#include <cassert>
const int NP = 4795;
#define MAX_LINE_LENGTH 1000


using namespace std;

namespace PLMD
{
  namespace colvar
  {

    class GNAP:public Colvar
    {

//                      PLMD::GNAPBase* gnap;
      bool squared;
      int i, j, k, k1, kex, l, m, n, NP1, indexall;	///////////these are necessary indicies that we shall exploit time to time
      int abc, abc1, abc2;
      int rows, cols;
      //double gnap;                                                                        ////////////////////////this is the global Non-affinity parameter that we need to calculate
      double dGNAPx[NP], dGNAPy[NP], dGNAPz[NP];	/////////////////////////these are the (3 x NP) number of "analytical" derivatives that we need to extract
      double dGNAPx_val[NP], dGNAPy_val[NP], dGNAPz_val[NP];	/////////////////////////these are the (3 x NP) number of  derivatives that we need to extract  
      char buf[MAX_LINE_LENGTH];

/////////////////////necessary stuffs to be read from the input files 
      int pi[NP];		///////////////////////////particle index for the NP particles

      double Rix[NP], Riy[NP], Riz[NP];	///////////////////////////X, Y, Z coordinates of each of the NP particles in the reference frame

      double Rfx[NP], Rfy[NP], Rfz[NP];	///////////////////////////X, Y, Z coordinates of each of the NP particles in the current(md) frame

      int *pi_of_neighboour[NP];	///////////////////////////particle index of the neighbourhood of each of the ( NP ) particle

      double *P_matrix[NP];	/////////////////////////// P  matrix for each of the ( NP ) particle 

      int no_of_neigh[NP];	///////////////////////////no of particles in the neighbourhood of each of the NP particles

    public:
        GNAP (const ActionOptions &);
       ~GNAP ();
      virtual void calculate ();
      static void registerKeywords (Keywords & keys);
    };


    using namespace std;

    //+PLUMEDOC DCOLVAR GNAP
    /*

     */
    //+ENDPLUMEDOC

      PLUMED_REGISTER_ACTION (GNAP, "GNAP")
      void GNAP::registerKeywords (Keywords & keys)
    {
      Colvar::registerKeywords (keys);
      keys.add ("compulsory", "PM", "a file containing all coefficients for the P matrix.");	//File Path for coefficients
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
      //* keys.add ("compulsory", "FRAME", "a file containing all coefficients for the P matrix.");      //File Path for coefficients
      keys.add ("atoms", "ATOMS",
		"the group of atoms that you are calculating the Gyration Tensor for");
      keys.add ("compulsory", "PI_NEIGHBOUR", "a file containing neighbours for all particles.");	// File path for neighbourhoods
      keys.add ("compulsory", "NUM_NEIGHBOUR",
		"number of neighbors per particle");
      keys.add ("compulsory", "REFERENCE",
		"file containing reference structure");
    }

    GNAP::GNAP (const ActionOptions & ao):PLUMED_COLVAR_INIT (ao),
      squared (false)
    {
      std::vector < AtomNumber > atoms;
      parseAtomList ("ATOMS", atoms);
      string reference;
        parse ("REFERENCE", reference);
      string pm;
        parse ("PM", pm);
      //////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
      //*string frame;
      //*  parse ("FRAME", frame);
      string pi_neigh;
        parse ("PI_NEIGHBOUR", pi_neigh);
      string num_neigh;
        parse ("NUM_NEIGHBOUR", num_neigh);
      //TODO:
      //1. parse all arguments into strings
      //2. use those strins to open files and load values into your matrices
      //




      // Parse all arguments into files
      // read all files to do this:
      ////////////////////    reading population_neigh.dat                                                         1                     reading            no_of_neigh[NP]
      FILE *population_neigh_dat;
        population_neigh_dat = std::fopen (num_neigh.c_str (), "r");
      for (int i = 0; i < NP; i++)
	{
	  fscanf (population_neigh_dat, "%d ", &no_of_neigh[i]);
	  //population_neigh_dat >> no_of_neigh[i];

	  m = no_of_neigh[i];
	  pi_of_neighboour[i] = new int[m];

	    m = (3 * no_of_neigh[i]) * (3 * no_of_neigh[i]);
	    P_matrix[i] = new double[m];
	}

////////////////////    reading PMATRIX.dat                                                                  2                     reading               P_matrix[NP]
      FILE *PMATRIX_dat;
      PMATRIX_dat = std::fopen (pm.c_str (), "r");
      for (int l = 0; l < NP; l++)
	{
	  rows = 3 * no_of_neigh[l];
	  cols = 3 * no_of_neigh[l];
	  for (i = 0; i < rows; i++)
	    {
	      for (j = 0; j < rows; j++)
		{
		  fscanf (PMATRIX_dat, "%lf ", &P_matrix[l][i * rows + j]);
		  //PMATRIX_dat >> P_matrix[l][ i * rows + j] ;
		}
	    }
	}

////////////////////    reading neighbourhood.dat                                                             3                       reading               pi_of_neighboour[NP]
      FILE *neighbourhood_dat;
      neighbourhood_dat = std::fopen (pi_neigh.c_str (), "r");
      for (int l = 0; l < NP; l++)
	{
	  rows = no_of_neigh[l];
	  for (i = 0; i < rows; i++)
	    {
	      fscanf (neighbourhood_dat, "%d ", &pi_of_neighboour[l][i]);
	      //neighbourhood_dat >> pi_of_neighboour[l][i] ;
	    }
	}

////////////////////    reading ref.gro                                                                       4                        reading               Rix[NP], Riy[NP], Riz[NP]
      FILE *ref_gro;
      ref_gro = std::fopen (reference.c_str (), "r");
      for (l = 0; l < NP; l++)
	{
	  fscanf (ref_gro, " %*s   %*s  %d %lf %lf %lf", &pi[l], &Rix[l],
		  &Riy[l], &Riz[l]);
	}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
///////////////please remove this later on
//*
/*
FILE *frame_dat;
      frame_dat = std::fopen (frame.c_str (), "r");
for (int l = 0; l < NP; l++)
	{
	  fscanf (frame_dat, " %*s   %*s  %*s  %*s  %*s   %lf %lf %lf %*s   %*s  ",  &Rfx[l] , &Rfy[l] , &Rfz[l] );
	
	}
*/






      // 1. Get PMatrix file into P_matrix
      // 2. get Number of neighbors into no_of_neighs[NP]
      // 3. Get the list of neighbors into  PI_Neighbours
      // 4. Get reference structure loaded into matrix
      checkRead ();


      //FILE* fp=fopen(file.c_str(),"r");
      addValueWithDerivatives ();
      setNotPeriodic ();
      
      //log.printf("  method for alignment : %s \n",type.c_str() );
      if (squared)
	log.printf ("  chosen to use SQUARED option for MSD instead of GNAP\n");

      requestAtoms (atoms);
    }

    GNAP::~GNAP ()
    {
      //delete gnap;
    }




    // calculator
    void GNAP::calculate ()
    {
      std::vector < Vector > derivatives (getNumberOfAtoms ());
      std::vector < Vector > posi_final (getNumberOfAtoms ());
      double gnapcal;

      //      plumed.getAtoms()
      for (unsigned ui = 0; ui < getNumberOfAtoms (); ui++)
	{
	  //pos+=(getMass(i)/mass)*getPosition(i);
	  posi_final[ui] = getPosition (ui);
	  Rfx[ui] = posi_final[ui][0];
	  Rfy[ui] = posi_final[ui][1];
	  Rfz[ui] = posi_final[ui][2];
	  ///Rfx[ui] = Rfx[ui]/10.0;
	  //Rfy[ui] = Rfy[ui]/10.0;
	  //Rfz[ui] = Rfz[ui]/10.0;
	  //log.printf("%d %lf %lf %lf \n",(ui+1), Rfx[ui], Rfy[ui], Rfz[ui]);
	}





/////////////////////  calculating  dGNAPx
      for (l = 0; l < getNumberOfAtoms (); l++)
	{
	  //------------part->1
	  dGNAPx[l] = 0.e0;
	  rows = 3 * no_of_neigh[l];
	  cols = 3 * no_of_neigh[l];
	  double *X_all;
	  X_all = new double[rows];
	  double *YnZ_all;
	  YnZ_all = new double[rows];
	  double *expression_all;
	  expression_all = new double[rows * cols];

	  for (i = 0; i < rows; i++)
	    for (j = 0; j < cols; j++)
	      expression_all[i * rows + j] = (0.e0);

	  k = 0;
	  kex = 0;
	  for (i = 0; i < no_of_neigh[l]; i++)
	    {
	      abc = pi_of_neighboour[l][i];
	      indexall = abc - 1 ;
	      X_all[k] =	-((Rfx[indexall] - Rix[indexall]) - (Rfx[l] - Rix[l]));
	      k++;
	      X_all[k] = 0.e0;
	      k++;
	      X_all[k] = 0.e0;
	      k++;	    
	     
	      YnZ_all[kex] = 0.e0;
	      kex++;
	      YnZ_all[kex] = -((Rfy[indexall] - Riy[indexall]) - (Rfy[l] - Riy[l]));
	      kex++;
	      YnZ_all[kex] = -((Rfz[indexall] - Riz[indexall]) - (Rfz[l] - Riz[l]));
	      kex++;
	    }


	  for (i = 0; i < rows; i = i + 3)
	    {
	      for (j = 0; j < rows; j++)
		{
		  expression_all[i * rows + j] += X_all[j];
		}
	    }

	  for (j = 0; j < rows; j = j + 3)
	    {
	      for (i = 0; i < rows; i++)
		{
		  expression_all[i * rows + j] += X_all[i];
		}
	    }

	  for (i = 0; i < rows; i = i + 3)
	    {
	      for (j = 0; j < rows; j++)
		{
		  expression_all[i * rows + j] += YnZ_all[j];
		}
	    }

	  for (j = 0; j < rows; j = j + 3)
	    {
	      for (i = 0; i < rows; i++)
		{
		  expression_all[i * rows + j] += YnZ_all[i];
		}
	    }




	  for (i = 0; i < rows; i++)
	    {
	      for (j = 0; j < cols; j++)
		{
		  dGNAPx[l] += expression_all[i * rows +
					       j] * P_matrix[l][i * rows + j];
		}
	    }

	  delete[]X_all;
	  delete[]YnZ_all;
    delete[]expression_all;


	  //------------part->2
	  for (i = 0; i < no_of_neigh[l]; i++)
	    {
	      abc = pi_of_neighboour[l][i];


	      rows = 3 * no_of_neigh[abc - 1];
	      cols = 3 * no_of_neigh[abc - 1];
	      double *XnYnZ_all;
	      XnYnZ_all = new double[rows];
	      double *expression_all;
	      expression_all = new double[rows * cols];


	      for (m = 0; m < rows; m++)
		for (n = 0; n < rows; n++)
		  expression_all[m * rows + n] = 0.e0;


	      k1 = 0;
	      for (j = 0; j < no_of_neigh[abc - 1]; j++)
		{
		  abc1 = pi_of_neighboour[abc - 1][j];
		  if (abc1 == pi[l])
		    {
		      abc2 = j;
		    }

		  XnYnZ_all[k1] =
		    ((Rfx[abc1 - 1] - Rix[abc1 - 1]) -
		     (Rfx[abc - 1] - Rix[abc - 1]));
		  k1++;
		  XnYnZ_all[k1] =
		    ((Rfy[abc1 - 1] - Riy[abc1 - 1]) -
		     (Rfy[abc - 1] - Riy[abc - 1]));
		  k1++;
		  XnYnZ_all[k1] =
		    ((Rfz[abc1 - 1] - Riz[abc1 - 1]) -
		     (Rfz[abc - 1] - Riz[abc - 1]));
		  k1++;
		}


	      for (j = 0; j < rows; j++)
		{
		  expression_all[((abc2 * 3) * rows + j)] += XnYnZ_all[j];
		}

	      for (j = 0; j < rows; j++)
		{
		  expression_all[(j * rows + (abc2 * 3))] += XnYnZ_all[j];
		}




	      for (m = 0; m < rows; m++)
		{
		  for (n = 0; n < cols; n++)
		    {
		      dGNAPx[l] += expression_all[m * rows +
						   n] * P_matrix[abc -
								 1][m * rows +
								    n];
		    }
		}



	      delete[]XnYnZ_all;
	  delete[]expression_all;
	    }




//log.printf (" dgnapx = %lf \n", dGNAPx[l]);
	}








/////////////////////  calculating  dGNAPy
      for (l = 0; l < getNumberOfAtoms (); l++)
	{
	  //------------part->1
	  dGNAPy[l] = 0.e0;
	  rows = 3 * no_of_neigh[l];
	  cols = 3 * no_of_neigh[l];
	  double *Y_all;
	  Y_all = new double[rows];
	  double *XnZ_all;
	  XnZ_all = new double[rows];
	  double *expression_all;
	  expression_all = new double[rows * cols];

	  for (i = 0; i < rows; i++)
	    for (j = 0; j < cols; j++)
	      expression_all[i * rows + j] = 0.e0;

	  k = 0;
	  kex = 0;
	  for (i = 0; i < no_of_neigh[l]; i++)
	    {
	      abc = pi_of_neighboour[l][i];
	      indexall = abc - 1 ;
	      Y_all[k] = 0.e0;
	      k++;
	      Y_all[k] =	-((Rfy[indexall] - Riy[indexall]) - (Rfy[l] - Riy[l]));
	      k++;
	      Y_all[k] = 0.e0;
	      k++;
	      
	      XnZ_all[kex] =	-((Rfx[indexall] - Rix[indexall]) - (Rfx[l] - Rix[l]));
	      kex++;
	      XnZ_all[kex] = 0.e0;
	      kex++;
	      XnZ_all[kex] =	-((Rfz[indexall] - Riz[indexall]) - (Rfz[l] - Riz[l]));
	      kex++;
	    }


	  for (i = 1; i < rows; i = i + 3)
	    {
	      for (j = 0; j < rows; j++)
		{
		  expression_all[i * rows + j] += Y_all[j];
		}
	    }

	  for (j = 1; j < rows; j = j + 3)
	    {
	      for (i = 0; i < rows; i++)
		{
		  expression_all[i * rows + j] += Y_all[i];
		}
	    }

	  for (i = 1; i < rows; i = i + 3)
	    {
	      for (j = 0; j < rows; j++)
		{
		  expression_all[i * rows + j] += XnZ_all[j];
		}
	    }

	  for (j = 1; j < rows; j = j + 3)
	    {
	      for (i = 0; i < rows; i++)
		{
		  expression_all[i * rows + j] += XnZ_all[i];
		}
	    }




	  for (i = 0; i < rows; i++)
	    {
	      for (j = 0; j < cols; j++)
		{
		  dGNAPy[l] += expression_all[i * rows +
					       j] * P_matrix[l][i * rows + j];
		}
	    }

	  delete[]Y_all;
	  delete[]XnZ_all;
    delete[]expression_all;


	  //------------part->2
	  for (i = 0; i < no_of_neigh[l]; i++)
	    {
	      abc = pi_of_neighboour[l][i];


	      rows = 3 * no_of_neigh[abc - 1];
	      cols = 3 * no_of_neigh[abc - 1];
	      double *XnYnZ_all;
	      XnYnZ_all = new double[rows];
	      double *expression_all;
	      expression_all = new double[rows * cols];


	      for (m = 0; m < rows; m++)
		for (n = 0; n < rows; n++)
		  expression_all[m * rows + n] = 0.e0;


	      k1 = 0;
	      for (j = 0; j < no_of_neigh[abc - 1]; j++)
		{
		  abc1 = pi_of_neighboour[abc - 1][j];
		  if (abc1 == pi[l])
		    {
		      abc2 = j;
		    }

		  XnYnZ_all[k1] =
		    ((Rfx[abc1 - 1] - Rix[abc1 - 1]) -
		     (Rfx[abc - 1] - Rix[abc - 1]));
		  k1++;
		  XnYnZ_all[k1] =
		    ((Rfy[abc1 - 1] - Riy[abc1 - 1]) -
		     (Rfy[abc - 1] - Riy[abc - 1]));
		  k1++;
		  XnYnZ_all[k1] =
		    ((Rfz[abc1 - 1] - Riz[abc1 - 1]) -
		     (Rfz[abc - 1] - Riz[abc - 1]));
		  k1++;
		}


	      for (j = 0; j < rows; j++)
		{
		  expression_all[((abc2 * 3 + 1) * rows + j)] +=
		    XnYnZ_all[j];
		}

	      for (j = 0; j < rows; j++)
		{
		  expression_all[(j * rows + (abc2 * 3 + 1))] +=
		    XnYnZ_all[j];
		}




	      for (m = 0; m < rows; m++)
		{
		  for (n = 0; n < cols; n++)
		    {
		      dGNAPy[l] += expression_all[m * rows +
						   n] * P_matrix[abc -
								 1][m * rows +
								    n];
		    }
		}



	      delete[]XnYnZ_all;
	  delete[]expression_all;
	    }



	}





/////////////////////  calculating  dGNAPz
      for (l = 0; l < getNumberOfAtoms (); l++)
	{
	  //------------part->1
	  dGNAPz[l] = 0.e0;
	  rows = 3 * no_of_neigh[l];
	  cols = 3 * no_of_neigh[l];
	  double *Z_all;
	  Z_all = new double[rows];
	  double *XnY_all;
	  XnY_all = new double[rows];
	  double *expression_all;
	  expression_all = new double[rows * cols];

	  for (i = 0; i < rows; i++)
	    for (j = 0; j < cols; j++)
	      expression_all[i * rows + j] = 0.e0;

	  k = 0;
	  kex = 0;
	  for (i = 0; i < no_of_neigh[l]; i++)
	    {
	      abc = pi_of_neighboour[l][i];
	      indexall = abc - 1 ;
	      Z_all[k] = 0.e0;
	      k++;
	      Z_all[k] = 0.e0;
	      k++;
	      Z_all[k] =	-((Rfz[indexall] - Riz[indexall]) - (Rfz[l] - Riz[l]));
	      k++;  
	      
	      XnY_all[kex] =	-((Rfx[indexall] - Rix[indexall]) - (Rfx[l] - Rix[l]));
	      kex++;
	      XnY_all[kex] =	-((Rfy[indexall] - Riy[indexall]) - (Rfy[l] - Riy[l]));
	      kex++;
	      XnY_all[kex] = 0.e0;
	      kex++;
	    }

	  for (i = 2; i < rows; i = i + 3)
	    {
	      for (j = 0; j < rows; j++)
		{
		  expression_all[i * rows + j] += Z_all[j];
		}
	    }

	  for (j = 2; j < rows; j = j + 3)
	    {
	      for (i = 0; i < rows; i++)
		{
		  expression_all[i * rows + j] += Z_all[i];
		}
	    }

	  for (i = 2; i < rows; i = i + 3)
	    {
	      for (j = 0; j < rows; j++)
		{
		  expression_all[i * rows + j] += XnY_all[j];
		}
	    }

	  for (j = 2; j < rows; j = j + 3)
	    {
	      for (i = 0; i < rows; i++)
		{
		  expression_all[i * rows + j] += XnY_all[i];
		}
	    }



	  for (i = 0; i < rows; i++)
	    {
	      for (j = 0; j < cols; j++)
		{
		  dGNAPz[l] += expression_all[i * rows +
					       j] * P_matrix[l][i * rows + j];
		}
	    }

	  delete[]Z_all;
	  delete[]XnY_all;
    delete[]expression_all;


	  //------------part->2
	  for (i = 0; i < no_of_neigh[l]; i++)
	    {
	      abc = pi_of_neighboour[l][i];


	      rows = 3 * no_of_neigh[abc - 1];
	      cols = 3 * no_of_neigh[abc - 1];
	      double *XnYnZ_all;
	      XnYnZ_all = new double[rows];
	      double *expression_all;
	      expression_all = new double[rows * cols];


	      for (m = 0; m < rows; m++)
		for (n = 0; n < rows; n++)
		  expression_all[m * rows + n] = 0.e0;



	      k1 = 0;
	      for (j = 0; j < no_of_neigh[abc - 1]; j++)
		{
		  abc1 = pi_of_neighboour[abc - 1][j];
		  if (abc1 == pi[l])
		    {
		      abc2 = j;
		    }

		  XnYnZ_all[k1] =
		    ((Rfx[abc1 - 1] - Rix[abc1 - 1]) -
		     (Rfx[abc - 1] - Rix[abc - 1]));
		  k1++;
		  XnYnZ_all[k1] =
		    ((Rfy[abc1 - 1] - Riy[abc1 - 1]) -
		     (Rfy[abc - 1] - Riy[abc - 1]));
		  k1++;
		  XnYnZ_all[k1] =
		    ((Rfz[abc1 - 1] - Riz[abc1 - 1]) -
		     (Rfz[abc - 1] - Riz[abc - 1]));
		  k1++;
		}


	      for (j = 0; j < rows; j++)
		{
		  expression_all[((abc2 * 3 + 2) * rows + j)] +=
		    XnYnZ_all[j];
		}

	      for (j = 0; j < rows; j++)
		{
		  expression_all[(j * rows + (abc2 * 3 + 2))] +=
		    XnYnZ_all[j];
		}



	      for (m = 0; m < rows; m++)
		{
		  for (n = 0; n < cols; n++)
		    {
		      dGNAPz[l] += expression_all[m * rows +
						   n] * P_matrix[abc -
								 1][m * rows +
								    n];
		    }
		}


	      delete[]XnYnZ_all;
	  delete[]expression_all;
	    }





	}






      //Vector deriv;
      for (m = 0; m < getNumberOfAtoms (); m++)
	{
//                      deriv.clear();
	  //dGNAPx_val[m] = dGNAPx[m];
	  //dGNAPy_val[m] = dGNAPy[m];
	  //dGNAPz_val[m] = dGNAPz[m];
	  derivatives[m][0] = (dGNAPx[m]);
	  derivatives[m][1] = (dGNAPy[m]);
	  derivatives[m][2] = (dGNAPz[m]);

	  //log.printf("particle : %d \n",(m+1));
	  //log.printf("dGNAPx = %lf  \n",derivatives[m][0]);
	  //log.printf("dGNAPy = %lf  \n",derivatives[m][1]);
	  //log.printf("dGNAPz = %lf  \n\n\n",derivatives[m][2]);
	  
	  
	  
	  setAtomsDerivatives (m, derivatives[m]);
	}









///////////////////calculating the gnapcal

      gnapcal = 0.0;
      for (l = 0; l < getNumberOfAtoms (); l++)
	{
	  rows = 3 * no_of_neigh[l];
	  cols = 3 * no_of_neigh[l];
	  double *XnYnZ_all;
	  XnYnZ_all = new double[rows];



	  k1 = 0;
	  for (j = 0; j < no_of_neigh[l]; j++)
	    {
	      abc1 = pi_of_neighboour[l][j];



	      XnYnZ_all[k1] =
		((Rfx[abc1 - 1] - Rix[abc1 - 1]) - (Rfx[l] - Rix[l]));
	      k1++;
	      XnYnZ_all[k1] =
		((Rfy[abc1 - 1] - Riy[abc1 - 1]) - (Rfy[l] - Riy[l]));
	      k1++;
	      XnYnZ_all[k1] =
		((Rfz[abc1 - 1] - Riz[abc1 - 1]) - (Rfz[l] - Riz[l]));
	      k1++;
	    }




	  for (m = 0; m < rows; m++)
	    {
	      for (n = 0; n < cols; n++)
		{
		  gnapcal =
		    gnapcal + P_matrix[l][m * rows +
					  n] * XnYnZ_all[m] * XnYnZ_all[n];
		}
	    }

      

	      delete[]XnYnZ_all;

	}

      gnapcal = (double) (gnapcal / (double) getNumberOfAtoms ());

      //log.printf (" gnapc = %lf ", gnapcal);






      setValue (gnapcal);


      setBoxDerivativesNoPbc ();
    }
    




  }
}






















