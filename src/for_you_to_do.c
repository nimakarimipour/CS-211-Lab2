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
    int i, max_index, t, j, tmp, k;
    double max;
    for (i = 0; i < n - 1; i++)
    {
        /*Pivoting*/
        max_index = i;
        max = fabs(A[i * n + i]);
        for (t = i + 1; t < n; t++)
        {
            if (fabs(A[t * n + i]) > max)
            {
                max_index = t;
                max = fabs(A[t * n + i]);
            }
        }
        if (max == 0)
        {
            return -1;
        }
        else if (max_index != i)
        {
            /*save pivot information*/
            tmp = ipiv[i];
            ipiv[i] = ipiv[max_index];
            ipiv[max_index] = tmp;
            /*swap rows*/
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
    if (UPLO == 'L')
    {
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
    }
    else if (UPLO == 'U')
    {
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
    int i, j, k, iB, jB, kB;
    for (k = 0; k < maty; k += b)
    {
        for (i = 0; i < matx; i += b)
        {
            for (j = 0; j < matx; j += b)
            {
                for (kB = k; kB < k + b && kB < maty; kB += 2)
                {
                    for (iB = i; iB < i + b && iB < matx; iB += 2)
                    {
                        register int regA00 = iB * n + kB;
                        register int regA10 = regA00 + n;
                        register double a00 = A[regA00], a01 = A[regA00 + 1], a10 = A[regA10], a11 = A[regA10 + 1];
                        for (jB = j; jB < j + b && jB < matx; jB += 2)
                        {
                            register int regB00 = kB * n + jB, regC00 = iB * n + jB;
                            register int regB10 = regB00 + n, regC10 = regC00 + n;
                            register double b00 = B[regB00], b01 = B[regB00 + 1], b10 = B[regB10], b11 = B[regB10 + 1];
                            register double c00 = C[regC00], c01 = C[regC00 + 1], c10 = C[regC10], c11 = C[regC10 + 1];
                            C[regC00] -= a00 * b00 + a01 * b10;
                            C[regC00 + 1] -= a00 * b01 + a01 * b11;
                            C[regC10] -= a10 * b00 + a11 * b10;
                            C[regC10 + 1] -= a10 * b01 + a11 * b11;
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
