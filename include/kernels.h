#ifndef __KERNELS_H__
#define __KERNELS_H__


#include "vars.h"
//#include "vars_gpu.h"
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cuComplex.h>
#include <fstream>
#include <iostream>


extern int *raysendstart;
extern int *rayrecvstart;
extern int *raysendcount;
extern int *rayrecvcount;

extern int *rayraystart;
extern int *rayrayind;
extern int *rayrecvlist;

extern double ftime;
extern double btime;
extern double fktime;
extern double frtime;
extern double bktime;
extern double brtime;
extern double aftime;
extern double abtime;

extern int numproj;
extern int numback;

extern int proj_rownztot;
extern int *proj_rowdispl;
extern int *proj_rowindex;
extern float *proj_rowvalue;
extern int proj_blocksize;
extern int proj_numblocks;
extern int proj_blocknztot;
extern int *proj_blockdispl;
extern int *proj_blockindex;
extern float *proj_blockvalue;
extern int proj_buffsize;
extern int *proj_buffdispl;
extern int proj_buffnztot;
extern int *proj_buffmap;
extern short *proj_buffindex;
extern float *proj_buffvalue;

extern int back_rownztot;
extern int *back_rowdispl;
extern int *back_rowindex;
extern float *back_rowvalue;
extern int back_blocksize;
extern int back_numblocks;
extern int back_blocknztot;
extern int *back_blockdispl;
extern int *back_blockindex;
extern float *back_blockvalue;
extern int back_buffsize;
extern int *back_buffdispl;
extern int back_buffnztot;
extern int *back_buffmap;
extern short *back_buffindex;
extern float *back_buffvalue;

int *proj_blockdispl_d;
int *proj_buffdispl_d;
int *proj_buffmap_d;
short *proj_buffindex_d;
float *proj_buffvalue_d;
int *back_blockdispl_d;
int *back_buffdispl_d;
int *back_buffmap_d;
short *back_buffindex_d;
float *back_buffvalue_d;

int *rayraystart_d;
int *rayrayind_d;
int *rayindray_d;

float *tomogram_d;
float *sinogram_d;
float *raypart_d;
float *raybuff_d;


extern mem_type mem;

extern const char *pidxfile; 
extern const char *pvalfile; 
extern const char *bidxfile; 
extern const char *bvalfile; 


extern float *raypart;
extern float *raybuff;
extern int numpix;
extern int numray;

__global__ void kernel_SpMV_buffered(float *y, float *x, short *index, float *value, int numrow, int *blockdispl, int *buffdispl, int *buffmap, int buffsize){
  extern __shared__ float shared[];
  float reduce = 0;
  int ind;
  

  for(int buff = blockdispl[blockIdx.x]; buff < blockdispl[blockIdx.x+1]; buff++){
    for(int i = threadIdx.x; i < buffsize; i += blockDim.x)
      shared[i] = x[buffmap[buff*buffsize+i]];
    __syncthreads();
    for(int n = buffdispl[buff]; n < buffdispl[buff+1]; n++){
      ind = n*blockDim.x+threadIdx.x;
      reduce = reduce + shared[index[ind]]*value[ind];
    }
    __syncthreads();
  }
  ind = blockIdx.x*blockDim.x+threadIdx.x;
  if(ind < numrow)
    y[ind] = reduce;
}
 

void setup_gpu(float **obj,float **gra, float **dir,float **mes,float **res,float **ray){

  int device = 0;
  printf("device: %d\n",device);
  cudaSetDevice(device);
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  printf("\n");
  printf("Device Count: %d\n",deviceCount);
  //for (int dev = 0; dev < deviceCount; dev++) {
    int dev = device;
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, dev);
    printf("Device %d name: %s\n",dev,deviceProp.name);
    printf("Computational Capabilities: %d, %d\n",deviceProp.major,deviceProp.minor);
    printf("Maximum global memory size: %lu\n",deviceProp.totalGlobalMem);
    printf("Maximum constant memory size: %zu\n",deviceProp.totalConstMem);
    printf("Maximum shared memory size per block: %zu\n",deviceProp.sharedMemPerBlock);
    printf("Maximum block dimensions: %dx%dx%d\n",deviceProp.maxThreadsDim[0],deviceProp.maxThreadsDim[1],deviceProp.maxThreadsDim[2]);
    printf("Maximum grid dimensions: %dx%dx%d\n",deviceProp.maxGridSize[0],deviceProp.maxGridSize[1],deviceProp.maxGridSize[2]);
    printf("Maximum threads per block: %d\n",deviceProp.maxThreadsPerBlock);
    printf("Warp size: %d\n",deviceProp.warpSize);
    printf("\n");
  //}

  cudaMallocHost((void**)obj,sizeof(float)*numpix);
  cudaMallocHost((void**)gra,sizeof(float)*numpix);
  cudaMallocHost((void**)dir,sizeof(float)*numpix);
  cudaMallocHost((void**)mes,sizeof(float)*numray);
  cudaMallocHost((void**)res,sizeof(float)*numray);
  cudaMallocHost((void**)ray,sizeof(float)*numray);

  float projmem = 0;
  projmem = projmem + sizeof(int)/1e9*(proj_numblocks+1);
  projmem = projmem + sizeof(int)/1e9*(proj_blocknztot+1);
  projmem = projmem + sizeof(int)/1e9*(proj_blocknztot*proj_buffsize);
  projmem = projmem + sizeof(int)/1e9*(proj_buffnztot*proj_blocksize);
  projmem = projmem + sizeof(float)/1e9*(proj_buffnztot*proj_blocksize);
  //printf("PROC %d FORWARD PROJECTION MEMORY: %f GB\n",myid,projmem);

  cudaMalloc((void**)&proj_blockdispl_d,sizeof(int)*(proj_numblocks+1));
  cudaMalloc((void**)&proj_buffdispl_d,sizeof(int)*(proj_blocknztot+1));
  cudaMalloc((void**)&proj_buffmap_d,sizeof(int)*proj_blocknztot*proj_buffsize);
  cudaMemcpy(proj_blockdispl_d,proj_blockdispl,sizeof(int)*(proj_numblocks+1),cudaMemcpyHostToDevice);
  cudaMemcpy(proj_buffdispl_d,proj_buffdispl,sizeof(int)*(proj_blocknztot+1),cudaMemcpyHostToDevice);
  cudaMemcpy(proj_buffmap_d,proj_buffmap,sizeof(int)*proj_blocknztot*proj_buffsize,cudaMemcpyHostToDevice);

  cudaMalloc((void**)&back_blockdispl_d,sizeof(int)*(back_numblocks+1));
  cudaMalloc((void**)&back_buffdispl_d,sizeof(int)*(back_blocknztot+1));
  cudaMalloc((void**)&back_buffmap_d,sizeof(int)*back_blocknztot*back_buffsize);
  cudaMemcpy(back_blockdispl_d,back_blockdispl,sizeof(int)*(back_numblocks+1),cudaMemcpyHostToDevice);
  cudaMemcpy(back_buffdispl_d,back_buffdispl,sizeof(int)*(back_blocknztot+1),cudaMemcpyHostToDevice);
  cudaMemcpy(back_buffmap_d,back_buffmap,sizeof(int)*back_blocknztot*back_buffsize,cudaMemcpyHostToDevice);

  std::ifstream fppidx,fppval,fpbidx,fpbval;
  uint64_t pidx_size = sizeof(short)*proj_buffnztot*proj_blocksize; 
  uint64_t pval_size = sizeof(float)*proj_buffnztot*proj_blocksize; 
  uint64_t bidx_size = sizeof(short)*back_buffnztot*back_blocksize; 
  uint64_t bval_size = sizeof(float)*back_buffnztot*back_blocksize; 

  switch(mem){
    case GPUMEM: 
        cudaMalloc((void**)&proj_buffindex_d,pidx_size);
        cudaMalloc((void**)&proj_buffvalue_d,pval_size);
        cudaMalloc((void**)&back_buffindex_d,bidx_size);
        cudaMalloc((void**)&back_buffvalue_d,bval_size);
        cudaMemcpy(proj_buffindex_d,proj_buffindex,pidx_size,cudaMemcpyHostToDevice);
        cudaMemcpy(proj_buffvalue_d,proj_buffvalue,pval_size,cudaMemcpyHostToDevice);
        cudaMemcpy(back_buffindex_d,back_buffindex,bidx_size,cudaMemcpyHostToDevice);
        cudaMemcpy(back_buffvalue_d,back_buffvalue,bval_size,cudaMemcpyHostToDevice);
        break; 
    case UVM_READONLY:
        cudaMallocManaged((void**)&proj_buffindex_d, pidx_size);
        cudaMallocManaged((void**)&proj_buffvalue_d, pval_size);
        cudaMemAdvise(proj_buffindex_d,pidx_size , cudaMemAdviseSetReadMostly, 0);
        cudaMemAdvise(proj_buffvalue_d,pval_size , cudaMemAdviseSetReadMostly, 0);

        cudaMallocManaged((void**)&back_buffindex_d,bidx_size);
        cudaMallocManaged((void**)&back_buffvalue_d,bval_size);
        cudaMemAdvise((void**)&back_buffindex_d,bidx_size, cudaMemAdviseSetReadMostly, 0);
        cudaMemAdvise((void**)&back_buffvalue_d,bval_size, cudaMemAdviseSetReadMostly, 0);

        fppidx.open(pidxfile, std::ios::in | std::ios::binary);
        fppval.open(pvalfile, std::ios::in | std::ios::binary);
        fpbidx.open(bidxfile, std::ios::in | std::ios::binary);
        fpbval.open(bvalfile, std::ios::in | std::ios::binary);
        if(!fppidx.is_open() || !fppval.is_open() || !fpbidx.is_open() || !fpbval.is_open()){
            fprintf(stderr, "File opening failed\n");
            exit(1);
        }
        // memcpy(proj_buffindex_d,proj_buffindex,sizeof(short)*proj_buffnztot*proj_blocksize);
        // memcpy(proj_buffvalue_d,proj_buffvalue,sizeof(float)*proj_buffnztot*proj_blocksize);
        // memcpy(back_buffindex_d,back_buffindex,sizeof(short)*back_buffnztot*back_blocksize);
        // memcpy(back_buffvalue_d,back_buffvalue,sizeof(float)*back_buffnztot*back_blocksize);
        fppidx.read((char*)proj_buffindex_d, pidx_size);
        fppval.read((char*)proj_buffvalue_d, pval_size);
        fpbidx.read((char*)back_buffindex_d, bidx_size);
        fpbval.read((char*)back_buffvalue_d, bval_size);
        break; 
    case UVM_DIRECT:  
        // printf("\n\n\n\n\nENTERING UVM DIRECT\n\n\n\n\n\n");
        cudaMallocManaged((void**)&proj_buffindex_d, pidx_size);
        cudaMallocManaged((void**)&proj_buffvalue_d, pval_size);
        cudaMemAdvise(proj_buffindex_d,pidx_size , cudaMemAdviseSetAccessedBy, 0);
        cudaMemAdvise(proj_buffvalue_d,pval_size , cudaMemAdviseSetAccessedBy, 0);

        cudaMallocManaged((void**)&back_buffindex_d,bidx_size);
        cudaMallocManaged((void**)&back_buffvalue_d,bval_size);
        cudaMemAdvise(back_buffindex_d, bidx_size, cudaMemAdviseSetAccessedBy, 0);
        cudaMemAdvise(back_buffvalue_d, bval_size, cudaMemAdviseSetAccessedBy, 0);

        fppidx.open(pidxfile, std::ios::in | std::ios::binary);
        fppval.open(pvalfile, std::ios::in | std::ios::binary);
        fpbidx.open(bidxfile, std::ios::in | std::ios::binary);
        fpbval.open(bvalfile, std::ios::in | std::ios::binary);
        if(!fppidx.is_open() || !fppval.is_open() || !fpbidx.is_open() || !fpbval.is_open()){
            fprintf(stderr, "File opening failed\n");
            exit(1);
        }
        // memcpy(proj_buffindex_d,proj_buffindex,sizeof(short)*proj_buffnztot*proj_blocksize);
        // memcpy(proj_buffvalue_d,proj_buffvalue,sizeof(float)*proj_buffnztot*proj_blocksize);
        // memcpy(back_buffindex_d,back_buffindex,sizeof(short)*back_buffnztot*back_blocksize);
        // memcpy(back_buffvalue_d,back_buffvalue,sizeof(float)*back_buffnztot*back_blocksize);
        fppidx.read((char*)proj_buffindex_d, pidx_size);
        fppval.read((char*)proj_buffvalue_d, pval_size);
        fpbidx.read((char*)back_buffindex_d, bidx_size);
        fpbval.read((char*)back_buffvalue_d, bval_size);
        break;
    case DRAGON_MAP: 
        // printf("\n\n\n\n\nENTERING DRAGON\n\n\n\n\n\n");
        if((dragon_map(pidxfile, pidx_size, D_F_READ, (void**) &proj_buffindex_d)) != D_OK){
                  printf("Dragon Map Failed for pidxfile\n");
                  exit(1);
        }
        if((dragon_map(pvalfile, pval_size, D_F_READ, (void**) &proj_buffvalue_d)) != D_OK){
                  printf("Dragon Map Failed for pvalfile\n");
                  exit(1);
        }
        if((dragon_map(bidxfile, bidx_size, D_F_READ, (void**) &back_buffindex_d)) != D_OK){
                  printf("Dragon Map Failed for bidxfilen\n");
                  exit(1);
        }
        if((dragon_map(bvalfile, bval_size, D_F_READ, (void**) &back_buffvalue_d)) != D_OK){
                  printf("Dragon Map Failed for bvalfile\n");
                  exit(1);
        }
        break; 


  }

  if(mem == UVM_DIRECT || mem == UVM_READONLY){
      fppidx.close();
      fppval.close();
      fpbidx.close();
      fpbval.close();
   }
  
  /*  
  FILE *fppidx = fopen((const char*)pidxfile, "wb");
  FILE *fppval = fopen((const char*)pvalfile, "wb");

  if(fppidx == NULL) {printf("Error open file projbufffile\n");}
  if(fppval == NULL) {printf("Error open file projbufffile\n");}

  printf("Size of pindex: %llu elements need %llu\n", proj_buffnztot*proj_blocksize, sizeof(short)*proj_buffnztot*proj_blocksize);
  printf("Size of pval: %llu elements need %llu\n", proj_buffnztot*proj_blocksize, sizeof(float)*proj_buffnztot*proj_blocksize);

  fwrite((void*)proj_buffindex_d,sizeof(short),proj_buffnztot*proj_blocksize,fppidx);
  fwrite((void*)proj_buffvalue_d,sizeof(float),proj_buffnztot*proj_blocksize,fppval);
*/
  
 
/*
  FILE *fpbidx = fopen((const char*)bidxfile, "wb");
  FILE *fpbval = fopen((const char*)bvalfile, "wb");

  if(fpbidx == NULL) {printf("Error open file projbufffile\n");}
  if(fpbval == NULL) {printf("Error open file projbufffile\n");}

  printf("Size of bindex: %llu elements need %f MB\n", back_buffnztot*back_blocksize, sizeof(short)*back_buffnztot*back_blocksize/(1024*1024.0));
  printf("Size of bval: %llu elements need %llu\n", back_buffnztot*back_blocksize, sizeof(float)*back_buffnztot*back_blocksize/(1024.0*1024));

  fwrite((void*)back_buffindex_d,sizeof(short),back_buffnztot*back_blocksize,fpbidx);
  fwrite((void*)back_buffvalue_d,sizeof(float),back_buffnztot*back_blocksize,fpbval);

*/


  float backmem = 0;
  backmem = backmem + sizeof(int)/1e9*(back_numblocks+1);
  backmem = backmem + sizeof(int)/1e9*(back_blocknztot+1);
  backmem = backmem + sizeof(int)/1e9*(back_blocknztot*back_buffsize);
  backmem = backmem + sizeof(int)/1e9*(back_buffnztot*back_blocksize);
  backmem = backmem + sizeof(float)/1e9*(back_buffnztot*back_blocksize);
  //printf("PROC %d BACKPROJECTION MEMORY: %f GB\n",myid,backmem);

  printf("TOTAL GPU MEMORY: %f GB\n",projmem+backmem);

  cudaMalloc((void**)&tomogram_d,sizeof(float)*numpix);
  cudaMalloc((void**)&sinogram_d,sizeof(float)*numray);

 }

void projection(float *mes, float *obj){
  double timef = omp_get_wtime();
  {
    int blocksize = proj_blocksize;
    int numblocks = proj_numblocks;
    int buffsize = proj_buffsize;
    int *blockdispl = proj_blockdispl_d;
    int *buffdispl = proj_buffdispl_d;
    int *buffmap = proj_buffmap_d;
    short *buffindex = proj_buffindex_d;
    float *buffvalue = proj_buffvalue_d;
    cudaMemcpy(tomogram_d,obj,sizeof(float)*numpix,cudaMemcpyHostToDevice);
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start);
    kernel_SpMV_buffered<<<numblocks,blocksize,sizeof(float)*buffsize>>>(sinogram_d,tomogram_d,buffindex,buffvalue,numray,blockdispl,buffdispl,buffmap,buffsize);
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  float milliseconds = 0;
  cudaEventElapsedTime(&milliseconds, start, stop);
  fktime = fktime + milliseconds/1000;
    cudaMemcpy(mes,sinogram_d,sizeof(float)*numray,cudaMemcpyDeviceToHost);
  }
  ftime = ftime + omp_get_wtime()-timef;
  numproj++;
}

void backprojection(float *gra, float *res){
  double timeb = omp_get_wtime();
  {
    int blocksize = back_blocksize;
    int numblocks = back_numblocks;
    int buffsize = back_buffsize;
    int *blockdispl = back_blockdispl_d;
    int *buffdispl = back_buffdispl_d;
    int *buffmap = back_buffmap_d;
    short *buffindex = back_buffindex_d;
    float *buffvalue = back_buffvalue_d;
    cudaMemcpy(sinogram_d,res,sizeof(float)*numray,cudaMemcpyHostToDevice);
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start);
    kernel_SpMV_buffered<<<numblocks,blocksize,sizeof(float)*buffsize>>>(tomogram_d,sinogram_d,buffindex,buffvalue,numpix,blockdispl,buffdispl,buffmap,buffsize);
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  float milliseconds = 0;
  cudaEventElapsedTime(&milliseconds, start, stop);
  bktime = bktime + milliseconds/1000;
    cudaMemcpy(gra,tomogram_d,sizeof(float)*numpix,cudaMemcpyDeviceToHost);
  }
  btime = btime + omp_get_wtime()-timeb;
  numback++;
}

#endif //__KERNELS_H__
