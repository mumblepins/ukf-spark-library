#include "ukf-spark-library.h"
#include <math.h>

/*
 * The MIT License (MIT)
 * 
 * Copyright (c) 2014 Daniel Sullivan (Mumblepins)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * */

UkfLib::UkfLib()
{
  return;
}

void UkfLib::ukf(
        FuncPtr f, float *x, float *P, FuncPtr hmeas, float *z, float *Q, float *R, //inputs
        int L, int m) //dimensions
{
  //x and P get modified in place

  float lambda = ALPHA * ALPHA * (L + KI) - L; //scaling factor
  float c = L + lambda; //scaling factor

  // <editor-fold defaultstate="collapsed" desc=" Memory Reservations ">
  float Wc[L * 2 + 1];
  float x1[L];
  float X1[L * (2 * L + 1)];
  float P1[L * L];
  float X2[L * (2 * L + 1)];
  float z1[m];
  float Z1[m * (2 * L + 1)];
  float P2[m * m];
  float Z2[m * (2 * L + 1)];
  float diagWc[(L * 2 + 1)*(L * 2 + 1)];
  float t1[L * (L * 2 + 1)];
  float P12[L * m];
  float K[L * m];
  float t3[m];
  float t4[L];
  float transposeP12[m * L];
  float t5[L * L];
  float transposeZ2[m * (2 * L + 1)];
  float X[L * (2 * L + 1)];
  memset(Wc, 0, sizeof (Wc));
  memset(x1, 0, sizeof (x1));
  memset(X1, 0, sizeof (X1));
  memset(P1, 0, sizeof (P1));
  memset(X2, 0, sizeof (X2));
  memset(z1, 0, sizeof (z1));
  memset(Z1, 0, sizeof (Z1));
  memset(P2, 0, sizeof (P2));
  memset(Z2, 0, sizeof (Z2));
  memset(diagWc, 0, sizeof (diagWc));
  memset(t1, 0, sizeof (t1));
  memset(P12, 0, sizeof (P12));
  memset(K, 0, sizeof (K));
  memset(t3, 0, sizeof (t3));
  memset(t4, 0, sizeof (t4));
  memset(transposeP12, 0, sizeof (transposeP12));
  memset(t5, 0, sizeof (t5));
  memset(transposeZ2, 0, sizeof (transposeZ2));
  memset(X, 0, sizeof (X));

  // </editor-fold>

  // Weights for covariance
  Wc[0] = lambda / c + (1 - ALPHA * ALPHA + BETA);
  for (int i = 1; i < L * 2 + 1; i++) {
    Wc[i] = 0.5 / c;
  }
  c = sqrt(c);

  // find sigma points
  sigmas(x, P, c, L, X);

  // unscented transform on process
  unscentedTransform(
          f, X, Wc, Q,
          L, (2 * L + 1),
          x1, X1, P1, X2);

  // unscented transform on measurements
  unscentedTransform(
          hmeas, X1, Wc, R,
          m, (2 * L + 1),
          z1, Z1, P2, Z2);

  // find transformed cross-covariance
  for (int i = 0; i < (L * 2 + 1); i++) {
    for (int j = 0; j < (L * 2 + 1); j++) {
      if (i == j)
        diagWc[i * (L * 2 + 1) + j] = Wc[i];
      else
        diagWc[i * (L * 2 + 1) + j] = 0;
    }
  }

  transpose(Z2, m, (2 * L + 1), transposeZ2);
  multiply(X2, diagWc, L, (L * 2 + 1), (L * 2 + 1), t1);
  multiply(t1, transposeZ2, L, (L * 2 + 1), m, P12);

  // find K
  invert(P2, m);
  multiply(P12, P2, L, m, m, K);

  //state update
  subtract(z, z1, m, 1, t3);
  multiply(K, t3, L, m, 1, t4);
  memset(x, 0, sizeof (x));
  add(x1, t4, L, 1, x);

  // covariance update
  transpose(P12, L, m, transposeP12);
  multiply(K, transposeP12, L, m, L, t5);
  memset(P, 0, sizeof (P));
  subtract(P1, t5, L, L, P);
}

void UkfLib::printArray(float* A, int m, int n, String label)
{
  // A = input matrix (m x n)
  int i, j;
  Serial.println();
  Serial.println(label);
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      Serial.print(A[n * i + j]);
      Serial.print("\t");
    }
    Serial.println();
  }
}

void UkfLib::unscentedTransform(
        FuncPtr f, float *X, float *Wc, float *R, // Inputs
        int n, int L, // dimensions
        float *y, float *Y, float *P, float *Y1) // Outputs
{
  /*
   * Unscented Transformation
   * Input:
   *         f: nonlinear map (maps Sn inputs to n outputs)
   *         X: sigma points (Sn x (L)) (Sn: number of state variables)
   *         Wm: weights for mean (L)
   *         Wc: weights for covraiance (L)
   *         R: additive covariance [n]
   * 
   *         n: number of outputs of f
   *         L: (1+2* Sn) (Sn: number of state variables) (note, different L from above function)
   * Output:
   *         y: transformed mean (n x 1)
   *         Y: transformed smapling points (n x L)
   *         P: transformed covariance (nxn)
   *        Y1: transformed deviations (n x L)
   * 
   *      */

  // <editor-fold defaultstate="collapsed" desc=" Memory Reservations ">

  float funcIn[(L - 1) / 2];
  float funcOut[n];
  float diagWc[L * L];
  float t1[n * L];
  float t2[n * n];
  float transposeY1[L * n];

  memset(y, 0, sizeof (y));
  memset(Y, 0, sizeof (Y));
  memset(funcIn, 0, sizeof (funcIn));
  memset(funcOut, 0, sizeof (funcOut));
  memset(diagWc, 0, sizeof (diagWc));
  memset(t1, 0, sizeof (t1));
  memset(t2, 0, sizeof (t2));
  memset(transposeY1, 0, sizeof (transposeY1));

  // </editor-fold>

  //determine transformed sampling points
  for (int k = 0; k < L; k++) {
    for (int i = 0; i < (L - 1) / 2; i++) {
      funcIn[i] = X[i * L + k];
    }
    f(funcIn, funcOut);
    for (int i = 0; i < n; i++) {
      Y[i * L + k] = funcOut[i];
    }
  }

  //determine transformed mean
  for (int i = 0; i < n; i++) {
    // row
    y[i] = Y[i * L] / 2;

    for (int j = 1; j < L; j++) {
      //column
      y[i] += Y[i * L + j] / ((L - 1)*2);
    }
  }

  // determine transformed deviations
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < L; j++) {
      Y1[i * L + j] = Y[i * L + j] - y[i];
    }
  }

  // determine transformed covariance
  // make a diagonal matrix
  for (int i = 0; i < L; i++) {
    for (int j = 0; j < L; j++) {
      if (i == j)
        diagWc[i * L + j] = Wc[i];
      else
        diagWc[i * L + j] = 0;

    }
  }
  transpose(Y1, n, L, transposeY1);
  multiply(Y1, diagWc, n, L, L, t1);
  multiply(t1, transposeY1, n, L, n, t2);
  add(t2, R, n, n, P);
}

void UkfLib::sigmas(float *x, float *P, float c, int n, float *X)
{
  /*
 % Sigma points around reference point
         % Inputs :
         % x : reference point
         % P : covariance
         % c : coefficient
         % Output :
         % X : Sigma points
   * */

  // <editor-fold defaultstate="collapsed" desc=" Memory Reservations ">

  float t1[n * n];
  float t2[n * n];
  float A[n * n];
  float Y[n * n];
  memset(t1, 0, sizeof (t1));
  memset(t2, 0, sizeof (t2));
  memset(A, 0, sizeof (A));
  memset(Y, 0, sizeof (Y));
  memset(X, 0, sizeof (X));

  // </editor-fold>

  cholesky(P, t1, n);
  transpose(t1, n, n, A);
  scale(A, n, n, c);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      Y[i * n + j] = x[i];
    }
  }
  add(Y, A, n, n, t1);
  subtract(Y, A, n, n, t2);
  for (int i = 0; i < n; i++) {
    X[i * (2 * n + 1)] = x[i];
    for (int j = 0; j < n; j++) {
      X[i * (2 * n + 1) + j + 1] = t1[i * n + j];
      X[i * (2 * n + 1) + j + n + 1] = t2[i * n + j];
    }
  }
}

void UkfLib::cholesky(float *A, float *L, int n)
{
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < (i + 1); j++) {
      float s = 0;
      for (int k = 0; k < j; k++)
        s += L[i * n + k] * L[j * n + k];
      L[i * n + j] = (i == j) ?
              sqrt(A[i * n + i] - s) :
              (1.0 / L[j * n + j] * (A[i * n + j] - s));
    }
  }
}

void UkfLib::copy(float* A, int n, int m, float* B)
{
  int i, j, k;
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      B[n * i + j] = A[n * i + j];
    }
}

void UkfLib::multiply(float* A, float* B, int m, int p, int n, float* C)
{
  //Matrix Multiplication Routine
  // C = A*B
  // A = input matrix (m x p)
  // B = input matrix (p x n)
  // m = number of rows in A
  // p = number of columns in A = number of rows in B
  // n = number of columns in B
  // C = output matrix = A*B (m x n)
  int i, j, k;
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      C[n * i + j] = 0;
      for (k = 0; k < p; k++)
        C[n * i + j] = C[n * i + j] + A[p * i + k] * B[n * k + j];
    }
}

void UkfLib::add(float* A, float* B, int m, int n, float* C)
{
  // Matrix Addition Routine
  // A = input matrix (m x n)
  // B = input matrix (m x n)
  // m = number of rows in A = number of rows in B
  // n = number of columns in A = number of columns in B
  // C = output matrix = A+B (m x n)
  int i, j;
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      C[n * i + j] = A[n * i + j] + B[n * i + j];
}

void UkfLib::subtract(float* A, float* B, int m, int n, float* C)
{
  // Matrix Subtraction Routine
  // A = input matrix (m x n)
  // B = input matrix (m x n)
  // m = number of rows in A = number of rows in B
  // n = number of columns in A = number of columns in B
  // C = output matrix = A-B (m x n)
  int i, j;
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      C[n * i + j] = A[n * i + j] - B[n * i + j];
}

void UkfLib::transpose(float* A, int m, int n, float* C)
{
  // Matrix Transpose Routine
  // A = input matrix (m x n)
  // m = number of rows in A
  // n = number of columns in A
  // C = output matrix = the transpose of A (n x m)
  int i, j;
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      C[m * j + i] = A[n * i + j];
}

void UkfLib::scale(float* A, int m, int n, float k)
{
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      A[n * i + j] = A[n * i + j] * k;
}

int UkfLib::invert(float* A, int n)
{
  //Matrix Inversion Routine
  // * This function inverts a matrix based on the Gauss Jordan method.
  // * Specifically, it uses partial pivoting to improve numeric stability.
  // * The algorithm is drawn from those presented in 
  //	 NUMERICAL RECIPES: The Art of Scientific Computing.
  // * The function returns 1 on success, 0 on failure.
  // * NOTE: The argument is ALSO the result matrix, meaning the input matrix is REPLACED
  // A = input matrix AND result matrix
  // n = number of rows = number of columns in A (n x n)
  int pivrow; // keeps track of current pivot row
  int k, i, j; // k: overall index along diagonal; i: row index; j: col index
  int pivrows[n]; // keeps track of rows swaps to undo at end
  float tmp; // used for finding max value and making column swaps

  for (k = 0; k < n; k++) {
    // find pivot row, the row with biggest entry in current column
    tmp = 0;
    for (i = k; i < n; i++) {
      if (abs(A[i * n + k]) >= tmp) // 'Avoid using other functions inside abs()?'
      {
        tmp = abs(A[i * n + k]);
        pivrow = i;
      }
    }

    // check for singular matrix
    if (A[pivrow * n + k] == 0.0f) {
      Serial.println("Inversion failed due to singular matrix");
      return 0;
    }

    // Execute pivot (row swap) if needed
    if (pivrow != k) {
      // swap row k with pivrow
      for (j = 0; j < n; j++) {
        tmp = A[k * n + j];
        A[k * n + j] = A[pivrow * n + j];
        A[pivrow * n + j] = tmp;
      }
    }
    pivrows[k] = pivrow; // record row swap (even if no swap happened)

    tmp = 1.0f / A[k * n + k]; // invert pivot element
    A[k * n + k] = 1.0f; // This element of input matrix becomes result matrix

    // Perform row reduction (divide every element by pivot)
    for (j = 0; j < n; j++) {
      A[k * n + j] = A[k * n + j] * tmp;
    }

    // Now eliminate all other entries in this column
    for (i = 0; i < n; i++) {
      if (i != k) {
        tmp = A[i * n + k];
        A[i * n + k] = 0.0f; // The other place where in matrix becomes result mat
        for (j = 0; j < n; j++) {
          A[i * n + j] = A[i * n + j] - A[k * n + j] * tmp;
        }
      }
    }
  }

  // Done, now need to undo pivot row swaps by doing column swaps in reverse order
  for (k = n - 1; k >= 0; k--) {
    if (pivrows[k] != k) {
      for (i = 0; i < n; i++) {
        tmp = A[i * n + k];
        A[i * n + k] = A[i * n + pivrows[k]];
        A[i * n + pivrows[k]] = tmp;
      }
    }
  }
  return 1;
}
