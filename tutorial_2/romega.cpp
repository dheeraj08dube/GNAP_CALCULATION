#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "myhead.h"
#define MAX_LINE_LENGTH 1000
//#define N 4722  //please mention the number of atom in your protein!!!!!! HERE!!!!!!


double periodic (double, double, double);

int main()
{ 
  
    static double Rn[NP][3],rn[NP][3], lx, ly, lz;
    char amino_name[NP][6], atom_name[NP][5];
    int atom_index[NP]; 
    int N1;
    char buf[MAX_LINE_LENGTH];
    int count,price,i,j,k,p,m;
    count=0;
    FILE *f1, *testfile;



    /* {reading the coordinates of the native protein into the variable} 
      the following part is meant for intakining all the coordinates of the native structure of the protein 
      into a varible. here we basically open the file (e.g. file by the name "c_n_coordinates.dat"). then we
      make sure whether or not the file is successfully opened which is accomplished by invoking the if-else
      statements following the file-opening step. then in order to demonstrate that the values are really 
      read into the variable we get the values printed on the screen and this part where we print the coordinates
      on the screen will be later commented out.*/
    f1 = fopen ("xyz.dat", "r");

for(i=0 ; i<NP ; i++)
    {
      atom_index[i] = i + 1 ;
      fscanf (f1, " %lf %lf %lf  ",  &Rn[i][0],&Rn[i][1],&Rn[i][2]);	
      //fscanf(f1, "%*s  %s  %d %lf %lf %lf",  atom_name[i], &atom_index[i],&Rn[i][0],&Rn[i][1],&Rn[i][2]);
     
  
    }
    fscanf(f1, "%lf %lf %lf ",  &lx, &ly, &lz);
    
fclose(f1);

double x[NP],y[NP],z[NP];
int newindex[NP];
N1=0;
for(i=0 ; i<NP ; i++)
    {
     //if(  ((atom_name[i][0])=='C') || ((atom_name[i][0])=='N') )   //carbon and nitrogen
    // if(  ((atom_name[i][0])=='C')  )                         //only carbon
    // if(  ((atom_name[i][0])=='C') || ((atom_name[i][0])=='N') || ((atom_name[i][0])=='O') ) //carbon nitrogen and oxygen
        {
        x[N1] = 1.0 * Rn[i][0];
        y[N1] = 1.0 * Rn[i][1];
        z[N1] = 1.0 * Rn[i][2];
        newindex[N1] = atom_index[i];
        N1++;
        }
    }



    /*The following 3 parts are actually meant for calculating R-omega, in successive steps*/
    

    /*Part-1: Here we calculate all the distances between all the possible particles and thereby
              create a square matrix having NxN enteries. Here again the printing part will be 
              later on commented out */
              static double distance[NP][NP];
              for(j=0 ; j<N1 ; j++)
                 {
                  for(i=j ; i<N1 ; i++)
                     {
                      if(i==j)
                      distance[i][j] = 0.0;
                      else
                      {
                      distance[i][j] = sqrt( periodic(x[i],x[j],lx) * periodic(x[i],x[j],lx)
                                           + periodic(y[i],y[j],ly) * periodic(y[i],y[j],ly) 
                                           + periodic(z[i],z[j],lz) * periodic(z[i],z[j],lz) );
                                
                      distance[j][i] = sqrt( periodic(x[i],x[j],lx) * periodic(x[i],x[j],lx)
                                           + periodic(y[i],y[j],ly) * periodic(y[i],y[j],ly) 
                                           + periodic(z[i],z[j],lz) * periodic(z[i],z[j],lz) );
                                           
                     
                      } 
                     }
                 }
                


    /*Part-2: Here we visit each column in the 2d "distance" array generated and sorting each of them 
              in ascending in the downward direction.*/
                  double temp;
              for(k=0 ; k<N1 ; k++)
{
                  for (i = 0; i < N1; ++i)
    {
        for (j = i + 1; j < N1; ++j)
        {
            if (distance[i][k] > distance[j][k])
            {
                temp =  distance[i][k];
                distance[i][k] = distance[j][k];
                distance[j][k] = temp;
            }
        }
    }
}             

FILE *fx;
//fx = fopen("ROMEGA_for_all_particles.dat","w");
//for(i=0;i<N1;i++)
//fprintf(fx,"%lf\n",distance[4][i]);
   


    /*Part-4: Now we sort the Array Romega in ascending order in the downward direction and capture the largest
              value in the file Romega.*/
                  for (i = 0; i < N1; ++i)
    {
        for (j = i + 1; j < N1; ++j)
        {
            if (distance[4][i] > distance[4][j])
            {
                temp =  distance[4][i];
                distance[4][i] = distance[4][j];
                distance[4][j] = temp;
            }
        }
    }

distance[4][NP-1]=1.010e0*distance[4][NP-1];
FILE *fp;
fp = fopen("Romega.dat","w");
fprintf(fp,"%lf",distance[4][N1-1]);
fclose(fp);


}

double
periodic (double a, double b, double box_len)
{
  double delta;
  delta = (a - b) - (box_len) * rint ((a - b) / (box_len));
  return (delta);
}






