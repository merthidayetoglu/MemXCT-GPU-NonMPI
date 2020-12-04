#ifndef __VARS_H__
#define __VARS_H__ 

#include <stdio.h>
#include <cmath>
#include <complex>
#include <limits>
#include <omp.h>
#include <dragon.h>

using namespace std;

struct matrix{
  int ind;
  float len;
};

typedef enum {
    GPUMEM = 0,
    UVM_READONLY = 1,
    UVM_DIRECT = 2,
    UVM_READONLY_NVLINK = 3,
    UVM_DIRECT_NVLINK = 4,
    DRAGON_MAP = 5,
} mem_type;

void findnumpix(float, float, float*, int*);
void findpixind(float, float, float*, int*, int, int*);
void findlength(float, float, float*, float*);

void projection(float*, float*);
void backprojection(float*, float*);

int encode(unsigned short, unsigned short);
int xy2d (int n, int x, int y);
void d2xy(int n, int d, int *x, int *y);

void setup_gpu(float**,float**,float**,float**,float**,float**);
float norm_kernel(float*, int);
float dot_kernel(float*, float*, int);
void copy_kernel(float*, float*, int);
void subtract_kernel(float*, float*, float*, int);
void saxpy_kernel(float*, float*, float, float*, int);


#endif //__VARS_H__
