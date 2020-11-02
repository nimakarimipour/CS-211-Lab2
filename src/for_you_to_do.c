#include "../include/for_you_to_do.h"

int get_block_size()
{
    //return the block size you'd like to use
    /*add your code here */
    return 128;
}

/**
 * 
 * this function computes LU factorization
 * for a square matrix
 * 
 * syntax 
 *  
 *  input : 
 *      A     n by n , square matrix
 *      ipiv  1 by n , vector
 *      n            , length of vector / size of matrix
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/

int mydgetrf(double *A, int *ipiv, int n)
{
    /* add your code here */
    int i = 0;
    int max_index = 0; 
    int t = 0;
    int j = 0;
    int tmp = 0;
    int k = 0;
    double max = 0;
    for (i = 0; i < n - 1; i++)
    {
        /*Pivoting*/
        max_index = i;
        max = A[i * n + i] > 0 ? A[i * n + i] : (-1 * A[i * n + i]);
        for (t = i + 1; t < n; t++)
        {
            double pos_value = A[t * n + i] > 0 ? A[t * n + i] : (-1 * A[t * n + i]);
            if (pos_value > max)
            {
                max = pos_value;
                max_index = t;
            }
        }
        if (max == 0)
        {
            return -1;
        }
        else if (max_index != i)
        {
            /*Save Pivot Information*/
            tmp = ipiv[i];
            ipiv[i] = ipiv[max_index];
            ipiv[max_index] = tmp;

            /*Swap Rows*/
            for (j = 0; j < n; j++)
            {
                double tempv;
                tempv = A[i * n + j];
                A[i * n + j] = A[max_index * n + j];
                A[max_index * n + j] = tempv;
            }
        }
        /*factorization*/
        for (j = i + 1; j < n; j++)
        {
            A[j * n + i] = (double)A[j * n + i] / A[i * n + i];
            for (k = i + 1; k < n; k++)
                A[j * n + k] = A[j * n + k] - A[j * n + i] * A[i * n + k];
        }
    }
    return 0;
}

/**
 * 
 * this function computes triangular matrix - vector solver
 * for a square matrix . according to lecture slides, this
 * function computes forward AND backward subtitution in the
 * same function.
 * 
 * syntax 
 *  
 *  input :
 *      UPLO  'L' or 'U' , denotes whether input matrix is upper
 *                         lower triangular . ( forward / backward
 *                         substitution )
 * 
 *      A     n by n     , square matrix
 * 
 *      B     1 by n     , vector
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *  
 *  output :
 *      none
 * 
 **/
void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv)
{
    double y[n];
    int i;
    int j;
    switch (UPLO)
    {
    case 'L':
        y[0] = B[ipiv[0]];
        for (i = 1; i < n; i++)
        {
            double sum = 0;
            for (j = 0; j < i; j++)
            {
                sum += y[j] * A[i * n + j];
            }
            y[i] = B[ipiv[i]] - sum;
        }
        break;

    case 'U':
        y[n - 1] = (double)B[n - 1] / A[(n - 1) * n + (n - 1)];
        for (i = n - 2; i >= 0; i--)
        {
            double sum = 0;
            for (j = i + 1; j < n; j++)
            {
                sum += y[j] * A[i * n + j];
            }
            y[i] = (B[i] - sum) / A[i * n + i];
        }
        break;
    }
    for (i = 0; i < n; i++)
    {
        B[i] = y[i];
    }
    return;
}

/**
 * 
 * Same function as what you used in lab1, cache_part4.c : optimal( ... ).
 * 
 **/
void mydgemm(double *A, double *B, double *C, int n, int matx, int maty, int b)
{
    int i = 0;
    int j = 0;
    int k = 0;
    int i_b = 0;
    int j_b = 0;
    int k_b = 0;
    for (k = 0; k < maty; k += b)
    {
        for (i = 0; i < matx; i += b)
        {
            for (j = 0; j < matx; j += b)
            { 
                for (k_b = k; k_b < k + b && k_b < maty; k_b += 2)
                {
                    for (i_b = i; i_b < i + b && i_b < matx; i_b += 2)
                    {
                        register int A_0_0 = i_b * n + k_b;
                        register int A_1_0 = A_0_0 + n;
                        register double a_0_0 = A[A_0_0];
                        register double a_0_1 = A[A_0_0 + 1];
                        register double a_1_0 = A[A_1_0];
                        register double a_1_1 = A[A_1_0 + 1];
                        for (j_b = j; j_b < j + b && j_b < matx; j_b += 2)
                        {
                            register int B_0_0 = k_b * n + j_b;
                            register int C_0_0 = i_b * n + j_b;
                            register int B_1_0 = B_0_0 + n;
                            register int C_1_0 = C_0_0 + n;
                            register double c_0_0 = C[C_0_0]; 
                            register double c_0_1 = C[C_0_0 + 1]; 
                            register double c_1_0 = C[C_1_0];
                            register double c_1_1 = C[C_1_0 + 1];
                            register double b_0_0 = B[B_0_0];
                            register double b_0_1 = B[B_0_0 + 1]; 
                            register double b_1_0 = B[B_1_0];
                            register double b_1_1 = B[B_1_0 + 1];
                            C[C_0_0] -= a_0_0 * b_0_0 + a_0_1 * b_1_0;
                            C[C_0_0 + 1] -= a_0_0 * b_0_1 + a_0_1 * b_1_1;
                            C[C_1_0] -= a_1_0 * b_0_0 + a_1_1 * b_1_0;
                            C[C_1_0 + 1] -= a_1_0 * b_0_1 + a_1_1 * b_1_1;
                        }
                    }
                }
            }
        }
    }
}

/**
 * 
 * this function computes LU decomposition
 * for a square matrix using block gepp introduced in course
 * lecture .
 * 
 * just implement the block algorithm you learned in class.
 * 
 * syntax 
 *  
 *  input :
 *      
 * 
 *      A     n by n     , square matrix
 * 
 *    
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *     
 *      b                , block size   
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf_block(double *A, int *ipiv, int n, int b)
{
    int i = 0;
    int t = 0;
    int j = 0;
    int i_b = 0;
    int last = 0;
    int max_index = 0;
    int k = 0;
    int tmp = 0;
    double max;
    for (i_b = 0; i_b < n - 1; i_b += b)
    {
        /*Partial Pivoting*/
        last = ((n - 1) > (i_b + -1 + b)) ? (i_b - 1 + b) : n - 1;
        for (i = i_b; i <= last; i++)
        {
            max_index = i;
            max = A[i * n + i] > 0 ? A[i * n + i] : (-1 * A[i * n + i]);
            for (k = i + 1; k < n; k++)
            {
                register double pos_elem = A[k * n + i] > 0 ? A[k * n + i] : (-1 * A[k * n + i]);
                if (pos_elem > max)
                {
                    max_index = k;
                    max = A[k * n + i] > 0 ? A[k * n + i] : (-1 * A[k * n + i]);
                }
            }
            if (max == 0)
            {
                return -1;
            }
            else if (max_index != i)
            {
                /*Save Pivot Infortmation*/
                tmp = ipiv[i];
                ipiv[i] = ipiv[max_index];
                ipiv[max_index] = tmp;
                /*Swap Rows*/
                for (j = 0; j < n; j++)
                {
                    double tempv;
                    tempv = A[i * n + j];
                    A[i * n + j] = A[max_index * n + j];
                    A[max_index * n + j] = tempv;
                }
            }
            /*Update columns i+1 until last*/
            for (j = i + 1; j < n; j++)
            {
                /* Since we are going to use the following values very often in indices, we can store the result in a registr */
                register int j_cross_n = j * n;
                register int i_cross_n = i * n;
                A[j_cross_n + i] = (double)A[j_cross_n + i] / A[i_cross_n + i];
                for (t = i + 1; t <= last; t++)
                {
                    A[j_cross_n + t] = A[j_cross_n + t] - A[j_cross_n + i] * A[i_cross_n + t];
                }
            }
        }
        /*inv(LL)*/
        for (i = i_b; i <= last; i++)
        {
            for (k = last + 1; k < n; k++)
            {
                double sum = 0;
                for (j = i_b; j < i; j++)
                {
                    sum += A[i * n + j] * A[j * n + k];
                }
                A[i * n + k] -= sum;
            }
        }
        /* copied and adapted from previous */
        mydgemm(&A[(last + 1) * n + i_b],
                &A[i_b * n + last + 1],
                &A[(last + 1) * n + (last + 1)],
                n,
                (n - last - 1), (last - i_b + 1),
                32);
    }
    return 0;
}
