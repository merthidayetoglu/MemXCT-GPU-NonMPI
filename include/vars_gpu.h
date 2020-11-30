#ifndef __VARS_GPU_H__
#define __VARS_GPU_H__


__global__ void kernel_SpMV_buffered(float*, float*, short*, float*, int, int*, int*, int*, int);
__global__ void kernel_SpReduce(float*, float*, int*, int*, int);
__global__ void kernel_SpGather(float*, float*, int*, int);


#endif //__VARS_GPU_H__
