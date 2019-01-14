/* 
Ian O'Brien
CS451 = Numerical Methods
matrix solver, (w/ multiple solution columns)
this program gets user input for a coefficnet matrix and then a solution matrix with the same number of rows
but as many colums as you want so you can solve systems with same coefficents but different solutions concurrently
*/ 

/*
functionality includes handling...
-checking for over and under determined systems
-checking for diagonal dominance
-divison by zero
-ability to partial pivot 
-ability to specify many solution columns to get as many solutions to a system as you want
*/

/*
possibile remaining bugs...
-when specifying many solutions to a system and at least one of them has no solutions and some do have solutions you dont get to see solutions to systems that are consistent, program bails before that happens
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// main function
int main(){

int row, col, srow, scol, i, j, p = 0, store = 0, flag = 1;
double a[10][10], b[10][10], temp, max = 0, pivot[10], diag, sum =0, pivot2[10];

// get size of coefficnet matrix (row X column)
printf("Enter number of rows COEF MATRIX (between 1 and 10): ");
scanf("%d", &row);
printf("Enter number of columns COEF MATRIX(between 1 and 10): ");
scanf("%d", &col);


// if more unknowns than equations, bail
if (row > col)
{
  printf("\n\n\tUnderdetermined system, please enter appropriate matrix...terminating");
  return 0;
}

// if more equations than unknowns, bail
if (col > row)
{
  printf("\n\n\tOverdetermined system, please enter appropriate matrix...terminating");
  return 0;
}

// loop to get elements of matrix
printf("\nEnter elements of coefficent matrix\n\n");

// build matrix A (coeffireplcent matrix)
for(i=0; i<row; ++i)
for(j=0; j<col; ++j)
{
  printf("Enter element a%d%d: ",i+1,j+1);
  scanf("%lf",&a[i][j]);
}

// check for diagonal dominance in coefficent matrix
for(i=0;i<row;++i)
for(j=0;j<col;++j)
{
  if(i == j)
  {
      diag = abs(a[i][j]);
      ++j;
  }
  sum += abs(a[i][j]);
  if(j==col-1)
  {
    if(diag < sum)
    {
      flag = 0;
    }
    sum = 0;
  }
}
if(flag == 1)
  printf("\n\t YES, DIAGONALLY DOMINANT");
else
  printf("\n\t NO, NOT DIAGONALLY DOMINANT");

// get size of solution matrix matrix (row X column)
// must have same number of rows, but can have many columns
srow = row;
printf("\n\nEnter number of columns SOLUTION MATRIX(between 1 and 10): ");
scanf("%d", &scol);

// loop to get elements of matrix B (solution matrix)
printf("\nEnter elements of solution matrix\n\n");

for(i=0; i<srow; ++i)
for(j=0; j<scol; ++j)
{
  printf("Enter element a%d%d: ",i+1,j+1);
  scanf("%lf",&b[i][j]);
}

// print matrix A (coefficent Matrix)
printf("\n coefficnet Matrix is: \n\n");
for(i=0;i<row;++i)
for(j=0;j<col;++j)
{
  printf("%lf ",a[i][j]);
  if(j==col-1)
  {
    printf("\n\n");
  }
}

// print matrix B (solution matrix)
printf("\n Solutin Matrix is: \n\n");
for(i=0;i<srow;++i)
for(j=0;j<scol;++j)
{
  printf("%lf ",b[i][j]);
  if(j==scol-1)
  {
    printf("\n\n");
  }
}

// get both matricies in RREF while checking for inconsistentcies and division by zero
 int lead = 0; 
    while (lead < row) 
    {

        // check for inconsistent system by checking for rows of 0
        // in coefficnet matrix and non-zero in column entries in solution 
        // matrix
        for(i=0;i<row;++i)
        {
          temp = 0;
          for(j=0;j<col;++j)
          {
            temp +=  abs(a[i][j]);
            if(j==col-1)
            {
              if (temp == 0)
              {
                if (b[i][j] != 0)
                {
                  printf("\n\n\t INCONSISTENT SYTEM WITH NO SOLUTION!! OR ONE SYSTEM IS SCALAR OF OTHER");
                  return 0;
                }
              }
            }
          }
        }

        // check for need to partial pivot
        for(int p = 0; p < row; p++)
        {
          for (int c = 0; c < (scol+col); c++)
          {
            // OKAY WE NEED TO DO PARTIAL PIVOT
            if((p ==c) && (a[p][c] == 0))
            {
              //STORE SWITCH ROW IN TEMP ARRAY
              for(int z = 0; z < (scol+col); ++z)
              {
                pivot[z] = a[p][z];
                pivot2[z] = b[p][z];
              }
              //FIND ROW TO SWITCH WITH BY FINDING ROW WITH MAX ABS VALUE IN COL OF ZERO
              // itration by row on column c to check for max abs and get row
              for(int z = 0; z < row; ++z)
              {
                if(abs(a[z][c] > max))
                {
                  max = abs(a[z][c]);
                  store = z;
                }
                // now we have row z that needs to be switched with 
              }
              // store row z in original row
              for(int z = 0; z < (scol+col); ++z)
              {
                a[p][z] = a[store][z];
                b[p][z] = b[store][z];
              }
              //move origianl row to row z
              for(int z = 0; z < (scol+col); ++z)
              {
                a[store][z] = pivot[z];
                b[store][z] = pivot2[z];
              }
            }
          }        
        }

        // for each row calculate divisor and multiplier 
        for (int r = 0; r < 10; r++) 
        {
            double d, m, y, z;
            d = a[lead][lead];
            m = a[r][lead] / a[lead][lead];

            for (int c = 0; c < (srow+col); c++)
            { // for each column ...
                if (r == lead)
                {
                    a[r][c] /= d;               // make pivot = 1
                    b[r][c] /= d;               // make pivot = 1
                }
                else
                {
                    a[r][c] -= a[lead][c] * m;  // make other = 0
                    b[r][c] -= b[lead][c] * m;  // make other = 0
                }
            }
        }
        // increment the lead for the Row reduction
        lead++;
    }

// print row reduced matrix
printf("\n RREF coefficent matrix 'A' is: \n\n");
for(i=0;i<row;++i)
for(j=0;j<col;++j)
{
  printf("%lf ",a[i][j]);
  if(j==col-1)
  {
    printf("\n\n");
  }
}

// print matrix B (solution matrix)
printf("\n Solutin Matrix is: \n\n");
for(i=0;i<srow;++i)
for(j=0;j<scol;++j)
{
  printf("%lf ",b[i][j]);
  if(j==scol-1)
  {
    printf("\n\n");
  }
}

//end of program
return 0;
}
