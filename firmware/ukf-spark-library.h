#ifndef _UKF_LIBRARY
#define _UKF_LIBRARY

#include "application.h"

/*
 * Unscented Karman Filter (UKF)
 * somewhat based on  Matlab code by Yi Cao
 * http://www.mathworks.com/matlabcentral/fileexchange/18217-learning-the-unscented-kalman-filter
 * 
 * Matrix math codes from http://playground.arduino.cc/Code/MatrixMath 
 * /////
 * MatrixMath library for Arduino
 * by Charlie Matlack
 * contact: eecharlie in Arduino forums
 * This library was modified from code posted by RobH45345,
 * notably including replacement of the inversion algorithm. See the book
 * NUMERICAL RECIPES: The Art of Scientific Computing.
 * /////
 * 
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

class UkfLib {
  typedef void (* FuncPtr) (float *, float *);
public:
  UkfLib();
  void ukf(
          FuncPtr f, float *x, float *P, FuncPtr hmeas, float *z, float *Q, float *R,
          int L, int m);
  void printArray(float* A, int m, int n, String label);

private:
  void unscentedTransform(
          FuncPtr f, float *X, float *Wc, float *R,
          int n, int L,
          float *y, float *Y, float *P, float *Y1);
  void sigmas(float *x, float *P, float c, int n, float *X);
  void cholesky(float* A, float* L, int n);
  void copy(float* A, int n, int m, float* B);
  void multiply(float* A, float* B, int m, int p, int n, float* C);
  void add(float* A, float* B, int m, int n, float* C);
  void subtract(float* A, float* B, int m, int n, float* C);
  void transpose(float* A, int m, int n, float* C);
  void scale(float* A, int m, int n, float k);
  int invert(float* A, int n);

  const float ALPHA = 0.1; //default, tunable
  const float KI = 0; //default, tunable
  const float BETA = 2; //default, tunable

};



#endif
