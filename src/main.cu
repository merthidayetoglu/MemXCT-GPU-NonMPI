#include <vars.h>
#include <raytrace.h>
#include <kernels.h>
#include <omp.h> 

double ftime = 0;
double btime = 0;
double fktime = 0;
double frtime = 0;
double bktime = 0;
double brtime = 0;
double aftime = 0;
double abtime = 0;

int numx; //X SIZE OF DOMAIN
int numy; //Y SIZE OF DOMAIN
int numthe; //NUMBER OF THETAS
int numrho; //NUMBER OF RHOS

float xstart; //X START OF DOMAIN
float ystart; //Y START OF DOMAIN
float pixsize; //PIXEL SIZE
float rhostart; //RHO START
float raylength; //RAYLENGTH

char *sinfile; //SINOGRAM FILE
char *thefile; //THETA FILE
char *outfile; //OUTPUT FILE

int numiter;

int spatindexing;
int specindexing;

int spatsize; //SPATIAL TILE SIZE
int specsize; //SPECTRAL TILE SIZE
int numxtile; //NUMBER OF X TILES
int numytile; //NUMBER OF Y TILES
int numttile; //NUMBER OF THETA TILES
int numrtile; //NUMBER OF RHO TILES
int numpix; //NUMBER OF PIXELS (EXTENDED)
int numray; //NUMBER OF RAYS (EXTENDED)
int numspattile; //NUMBER OF SPATIAL TILES
int numspectile; //NUMBER OF SPECTRAL TILES

float *raypart;
float *raybuff;

int *raysendstart;
int *rayrecvstart;
int *raysendcount;
int *rayrecvcount;

int *rayraystart;
int *rayrayind;

int *rayrecvlist;

int proj_rownztot;
int *proj_rowdispl;
int *proj_rowindex;
float *proj_rowvalue;
int proj_blocksize;
int proj_numblocks;
int proj_blocknztot;
int *proj_blockdispl;
int *proj_blockindex;
float *proj_blockvalue;
int proj_buffsize;
int *proj_buffdispl;
int proj_buffnztot;
int *proj_buffmap;
short *proj_buffindex;
float *proj_buffvalue;

int *back_rowdispl;
int *back_rowindex;
float *back_rowvalue;
int back_blocksize;
int back_numblocks;
int back_blocknztot;
int *back_blockdispl;
int *back_blockindex;
float *back_blockvalue;
int back_buffsize;
int *back_buffdispl;
int back_buffnztot;
int *back_buffmap;
short *back_buffindex;
float *back_buffvalue;

int numproj;
int numback;

int main(int argc, char** argv){

  double timetot = omp_get_wtime();

  int numthreads;
  #pragma omp parallel
  if(omp_get_thread_num()==0)
          numthreads=omp_get_num_threads();

  //SCANNING GEOMETRY DATA
  char *chartemp;
  chartemp = getenv("NUMTHE");
  numthe = atoi(chartemp);
  chartemp = getenv("NUMRHO");
  numrho = atoi(chartemp);
  //chartemp = getenv("NUMX");
  numx = numrho;//atoi(chartemp);
  //chartemp = getenv("NUMY");
  numy = numrho;//atoi(chartemp);

  //chartemp = getenv("XSTART");
  float xstart = -numrho/2.0;//atof(chartemp);
  //chartemp = getenv("YSTART");
  float ystart = xstart;//atof(chartemp);
  chartemp = getenv("PIXSIZE");
  pixsize = atof(chartemp);
  //chartemp = getenv("RHOSTART");
  rhostart = xstart;//atof(chartemp);
  //chartemp = getenv("RAYLENGTH");
  raylength = 2*numx;// atof(chartemp);

  //chartemp = getenv("SPATINDEXING");
  spatindexing = 5;//atoi(chartemp);
  //chartemp = getenv("SPECINDEXING");
  specindexing = 5;//atoi(chartemp);

  chartemp = getenv("NUMITER");
  numiter = atoi(chartemp);

  chartemp = getenv("SPATSIZE");
  spatsize = atof(chartemp);
  chartemp = getenv("SPECSIZE");
  specsize = atof(chartemp);

  chartemp = getenv("PROJBLOCK");
  proj_blocksize = atoi(chartemp);
  chartemp = getenv("BACKBLOCK");
  back_blocksize = atoi(chartemp);
  chartemp = getenv("PROJBUFF");
  proj_buffsize = atoi(chartemp);
  chartemp = getenv("BACKBUFF");
  back_buffsize = atoi(chartemp);

  sinfile = getenv("SINFILE");
  thefile = getenv("THEFILE");
  outfile = getenv("OUTFILE");

  //FIND NUMBER OF TILES
  numxtile = numx/spatsize;
  if(numx%spatsize)numxtile++;
  numytile = numy/spatsize;
  if(numy%spatsize)numytile++;
  numrtile = numrho/specsize;
  if(numrho%specsize)numrtile++;
  numttile = numthe/specsize;
  if(numthe%specsize)numttile++;
  numspattile = numxtile*numytile;
  numspectile = numrtile*numttile;
  numpix = numspattile*pow(spatsize,2);
  numray = numspectile*pow(specsize,2);
  //PRINT DATA
  printf("NUM. THETA             : %d\n",numthe);
  printf("NUM. RHO               : %d\n",numrho);
  printf("NUM. X PIXELS          : %d\n",numx);
  printf("NUM. Y PIXELS          : %d\n",numy);
  printf("\n");
  printf("NUM. OF PIXELS         : %d\n",numx*numy);
  printf("NUM. OF RAYS           : %d\n",numrho*numthe);
  printf("\n");
  printf("NUM. ITERATIONS        : %d\n",numiter);
  printf("\n");
  printf("SPATIAL TILE SIZE      : %d\n",spatsize);
  printf("SPECTRAL TILE SIZE     : %d\n",specsize);
  printf("\n");
  printf("NUMBER OF X TILES      : %d\n",numxtile);
  printf("NUMBER OF Y TILES      : %d\n",numytile);
  printf("NUM. OF THETA TILES    : %d\n",numttile);
  printf("NUM. OF RHO TILES      : %d\n",numrtile);
  printf("\n");
  printf("NUM. SPATIAL TILES     : %d\n",numspattile);
  printf("NUM. SPECTRAL TILES    : %d\n",numspectile);
  printf("\n");
  printf("NUM. OF X PIXELS (EXT) : %d\n",numxtile*spatsize);
  printf("NUM. OF Y PIXELS (EXT) : %d\n",numytile*spatsize);
  printf("NUM. OF ANGLES (EXT)   : %d\n",numttile*specsize);
  printf("NUM. OF RHOS (EXT)     : %d\n",numrtile*specsize);
  printf("\n");
  printf("NUM. OF PIXELS (EXT)   : %d\n",numpix);
  printf("NUM. OF RAYS (EXT)     : %d\n",numray);
  printf("\n");
  printf("NUMBER OF THREADS      : %d\n",numthreads);
  printf("\n");
  printf("INTEGER: %d, FLOAT: %d, LONG: %d, SHORT: %d, POINTER: %d\n",(int)sizeof(int),(int)sizeof(float),(int)sizeof(long),(int)sizeof(short),(int)sizeof(complex<float>*));
  printf("APPROX. MEMORY     TOTAL: %e GB\n",(numthe*numrho)/1024.0/1024.0/1024.0*numrho*16.0*1.5);
  printf("\n");
  printf("X START       : %e\n",xstart);
  printf("Y START       : %e\n",ystart);
  printf("PIXEL SIZE    : %e\n",pixsize);
  printf("RHO START     : %e\n",rhostart);
  printf("RAY LENGTH    : %e\n",raylength);
  printf("\n");
  printf("SPATIAL INDEXING       : %d\n",spatindexing);
  printf("SPECTRAL INDEXING      : %d\n",specindexing);
  printf(" 1: CARTESIAN, NATURAL\n 2: CARTESIAN, TRANSPOSED\n 3: MORTON, NATURAL\n 4: MORTON, TRANSPOSED\n");
  printf("PROJECTION BLOCK SIZE      : %d\n",proj_blocksize);
  printf("BACKPROJECTION BLOCK SIZE  : %d\n",back_blocksize);
  printf("PROJECTION BUFFER SIZE      : %d KB\n",proj_buffsize);
  printf("BACKPROJECTION BUFFER SIZE  : %d KB\n",back_buffsize);
  printf("\n");
  printf("SINOGRAM FILE : %s\n",sinfile);
  printf("   THETA FILE : %s\n",thefile);
  printf("  OUTPUT FILE : %s\n",outfile);
  printf("\n");
  printf("\n");
  proj_buffsize = proj_buffsize*1024/4;
  back_buffsize = back_buffsize*1024/4;

  printf("PLACE TILES\n");
  //PLACE SPATIAL TILES
  int lspat = numxtile;
  if(numytile > lspat)lspat = numytile;
  int spatlevel = 0;
  while(true){
    if(lspat<=pow(2,spatlevel))break;
    spatlevel++;
  }
  int lspatdim = pow(2,spatlevel);
  printf("lspat %d lspatdim %d\n",lspat,lspatdim);
  complex<float> *spatlltemp = new complex<float>[lspatdim*lspatdim];
  #pragma omp parallel for
  for(int lspat = 0; lspat < lspatdim*lspatdim; lspat++)
    spatlltemp[lspat].real(xstart-1);
  #pragma omp parallel for
  for(int spat = 0; spat < numspattile; spat++){
    int ytile = spat/numxtile;
    int xtile = spat%numxtile;
    int ind = 0;
      //ind = ytile*numxtile+xtile;
      //ind = encode(xtile,ytile);
      ind = xy2d (numxtile,xtile,ytile);
    float x = xstart+xtile*spatsize*pixsize;
    float y = ystart+ytile*spatsize*pixsize;
    spatlltemp[ind]=complex<float>(x,y);
  }
  complex<float> *spatll = new complex<float>[numspattile];
  int spatcount = 0;
  for(int lspat = 0; lspat < lspatdim*lspatdim; lspat++)
    if(spatlltemp[lspat].real()>xstart-0.5){
      spatll[spatcount] = spatlltemp[lspat];
      spatcount++;
    }
  delete[] spatlltemp;
  //PLACE SPECTRAL TILES
  int lspec = numrtile;
  if(numttile > lspec)lspec = numttile;
  int speclevel = 0;
  while(true){
    if(lspec<=pow(2,speclevel))break;
    speclevel++;
  }
  int lspecdim = pow(2,speclevel);
  printf("lspec %d lspecdim %d\n",lspec,lspecdim);
  complex<float> *speclltemp = new complex<float>[lspecdim*lspecdim];
  #pragma omp parallel for
  for(int lspec = 0; lspec < lspecdim*lspecdim; lspec++)
    speclltemp[lspec].real(rhostart-1);
  #pragma omp parallel for
  for(int spec = 0; spec < numspectile; spec++){
    int thetile = spec/numrtile;
    int rhotile = spec%numrtile;
    int ind = 0;
      //ind = thetile*numrtile+rhotile;
      //ind = encode(rhotile,thetile);
      ind = xy2d(lspecdim,rhotile,thetile);
    float rho = rhostart+rhotile*specsize;
    float the = thetile*specsize*M_PI/numthe;
    speclltemp[ind]=complex<float>(rho,the);
  }
  complex<float> *specll = new complex<float>[numspectile];
  int speccount = 0;
  for(int lspec = 0; lspec < lspecdim*lspecdim; lspec++)
    if(speclltemp[lspec].real()>rhostart-0.5){
      specll[speccount] = speclltemp[lspec];
      speccount++;
    }
  delete[] speclltemp;
  printf("FILL PIXELS AND RAYS\n");
  //PLACE PIXELS
  complex<float> *pixcoor = new complex<float>[numpix];
  int *pixglobalind = new int[numpix];
  #pragma omp parallel for
  for(int pix = 0; pix < numpix; pix++){
    int tile = pix/(spatsize*spatsize);
    int pixloc = pix%(spatsize*spatsize);
    int pixlocy = pixloc/spatsize;
    int pixlocx = pixloc%spatsize;
    int ind = 0;
    if(spatindexing==1)
      ind = pixlocy*spatsize+pixlocx;
    if(spatindexing==2)
      ind = pixlocx*spatsize+pixlocy;
    if(spatindexing==3)
      ind = encode(pixlocx,pixlocy);
    if(spatindexing==4)
      ind = encode(pixlocy,pixlocx);
    if(spatindexing==5)
      //ind = encode(pixlocx,pixlocy);
      ind = xy2d(spatsize,pixlocx,pixlocy);
    ind = ind + tile*spatsize*spatsize;
    float x = spatll[tile].real()+pixsize/2+pixlocx*pixsize;
    float y = spatll[tile].imag()+pixsize/2+pixlocy*pixsize;
    pixcoor[ind] = complex<float>(x,y);
    //GLOBAL SPATIAL INDEX (EXTENDED)
    int xglobalind = (int)((x-xstart)/pixsize);
    int yglobalind = (int)((y-ystart)/pixsize);
    pixglobalind[ind] = yglobalind*numxtile*spatsize+xglobalind;
  }
 
  float *mestheta = new float[numthe];
  printf("INPUT THETA DATA: %s\n",thefile);
  FILE *thetaf = fopen(thefile,"rb");
  fread(mestheta,sizeof(float),numthe,thetaf);
  fclose(thetaf);

  //PLACE RAYS
  complex<float> *raycoor = new complex<float>[numray];
  int *rayglobalind = new int[numray];
  int *raymesind = new int[numray];
  #pragma omp parallel for
  for(int ray = 0; ray < numray; ray++){
    int tile = ray/(specsize*specsize);
    int rayloc = ray%(specsize*specsize);
    int raylocthe = rayloc/specsize;
    int raylocrho = rayloc%specsize;
    int ind = 0;
    if(specindexing==1)
      ind = raylocthe*specsize+raylocrho;
    if(specindexing==2)
      ind = raylocrho*specsize+raylocthe;
    if(specindexing==3)
      ind = encode(raylocrho,raylocthe);
    if(specindexing==4)
      ind = encode(raylocthe,raylocrho);
    if(specindexing==5)
      //ind = encode(raylocrho,raylocthe);
      ind = xy2d(specsize,raylocrho,raylocthe);
    ind = ind + tile*specsize*specsize;
    float rho = specll[tile].real()+0.5+raylocrho;
    float the = specll[tile].imag()+raylocthe*M_PI/numthe;
    //GLOBAL SPECTRAL INDEX (EXTENDED)
    int rhoglobalind = (int)(rho-rhostart);
    int theglobalind = (int)((the+(M_PI/numthe)/2)/(M_PI/numthe));
    rayglobalind[ind] = theglobalind*numrtile*specsize+rhoglobalind;
    if(theglobalind<numthe && rhoglobalind<numrho){
      raymesind[ind] = theglobalind*numrho+rhoglobalind;
      raycoor[ind] = complex<float>(rho,mestheta[theglobalind]);
      //raycoor[ind] = complex<float>(rho,the);
    }
    else{
      raycoor[ind].real(5*raylength);
      raymesind[ind] = -1;
    }
  }
  //delete[] mestheta;
  delete[] specll;
  printf("SKIP DOMAIN PARTITIONING\n"); 
  double timep = omp_get_wtime();
  double project_time = omp_get_wtime();
  printf("\nCONSTRUCT PROJECTION MATRIX\n");
  {
    int *rownz = new int[numray];
    #pragma omp parallel for schedule(dynamic,proj_blocksize)
    for(int k = 0; k < numray; k++){
      float rho = raycoor[k].real();
      float theta = raycoor[k].imag();
      rownz[k] = 0;
      for(int tile = 0; tile < numspattile; tile++){
        float domain[4];
        domain[0]=spatll[tile].real();
        domain[1]=domain[0]+spatsize*pixsize;
        domain[2]=spatll[tile].imag();
        domain[3]=domain[2]+spatsize*pixsize;
        //REMOVE SPATIAL EDGE CONDITION
        if(domain[1] > xstart+numx*pixsize)domain[1]=xstart+numx*pixsize;
        if(domain[3] > ystart+numy*pixsize)domain[3]=ystart+numy*pixsize;
        findnumpix(theta,rho,&domain[0],&rownz[k]);
      }
    }
    int rownzmax = 0;
    for(int k = 0; k < numray; k++)
      if(rownz[k]>rownzmax)rownzmax=rownz[k];
    int *rowdispl = new int[numray+1];
    rowdispl[0] = 0;
    for(int k = 1; k < numray+1; k++)
      rowdispl[k] = rowdispl[k-1]+rownz[k-1];
    delete[] rownz;
    int rownztot = rowdispl[numray];
    printf("CSR STORAGE: %d (%f GB) rownzmax %d\n",rownztot,rownztot*sizeof(float)*2/1024.0/1024.0/1024.0,rownzmax);
    int *rowindex = new int[rownztot];
    #pragma omp parallel for schedule(dynamic,proj_blocksize)
    for(int k = 0; k < numray; k++){
      float rho = raycoor[k].real();
      float theta = raycoor[k].imag();
      int start = rowdispl[k];
      for(int tile = 0; tile < numspattile; tile++){
        float domain[4];
        domain[0]=spatll[tile].real();
        domain[1]=domain[0]+spatsize*pixsize;
        domain[2]=spatll[tile].imag();
        domain[3]=domain[2]+spatsize*pixsize;
        //REMOVE SPATIAL EDGE CONDITION
        if(domain[1] > xstart+numx*pixsize)domain[1]=xstart+numx*pixsize;
        if(domain[3] > ystart+numy*pixsize)domain[3]=ystart+numy*pixsize;
        int offset = tile*spatsize*spatsize;
        int pixtemp = 0;
        findpixind(theta,rho,&domain[0],&pixtemp,offset,&rowindex[start]);
        start=start+pixtemp;
      }
    }
    proj_rownztot = rownztot;
    proj_rowdispl = rowdispl;
    proj_rowindex = rowindex;
  }
  printf("RAY-TRACING TIME: %e\n",omp_get_wtime()-project_time);
  project_time = omp_get_wtime();
  printf("CONSTRUCT BACKPROJECTION MATRIX\n");
  {
    int *csrRowInd = new int[proj_rownztot];
    int *inter = new int[(numthreads+1)*numpix];;
    int *intra = new int[proj_rownztot];
    #pragma omp parallel for
    for(int k = 0; k < numray; k++)
      for(int n = proj_rowdispl[k]; n < proj_rowdispl[k+1]; n++)
        csrRowInd[n] = k;
    #pragma omp parallel for
    for(int n = 0; n < (numthreads+1)*numpix; n++)
      inter[n] = 0;
    #pragma omp parallel for
    for(int n = 0; n < proj_rownztot; n++){
      intra[n] = inter[(omp_get_thread_num()+1)*numpix+proj_rowindex[n]];
      inter[(omp_get_thread_num()+1)*numpix+proj_rowindex[n]]++;
    }
    #pragma omp parallel for
    for(int m = 0; m < numpix; m++)
      for(int t = 1; t < numthreads+1; t++)
        inter[t*numpix+m] = inter[t*numpix+m]+inter[(t-1)*numpix+m];
    int *rowdispl = new int[numpix+1];
    rowdispl[0] = 0;
    for(int m = 1; m < numpix+1; m++)
      rowdispl[m] = rowdispl[m-1] + inter[numthreads*numpix+m-1];
    int rownzmax = 0;
    for(int k = 0; k < numpix; k++){
      int rownz = rowdispl[k+1]-rowdispl[k];
      if(rownz>rownzmax)rownzmax=rownz;
    }
    int rownztot = rowdispl[numpix];
    printf("CSR STORAGE: %d (%f GB) rownzmax %d\n",rownztot,rownztot*sizeof(float)*2/1024.0/1024.0/1024.0,rownzmax);
    int *rowindex = new int[rownztot];
    #pragma omp parallel for
    for(int n = 0; n < rownztot; n++){
      rowindex[rowdispl[proj_rowindex[n]]+
      inter[omp_get_thread_num()*numpix+
      proj_rowindex[n]]+intra[n]] = csrRowInd[n];
    }
    delete[] inter;
    delete[] intra;
    delete[] csrRowInd;
    back_rowdispl = rowdispl;
    back_rowindex = rowindex;
  }
  printf("TRANSPOSITION TIME: %e\n",omp_get_wtime()-project_time);
  project_time = omp_get_wtime();
  printf("\nBLOCKING PROJECTION MATRIX\n");
  {
    int *rowindex = proj_rowindex;
    int *rowdispl = proj_rowdispl;
    int blocksize = proj_blocksize;
    int buffsize = proj_buffsize;
    int numblocks = numray/blocksize;
    if(numray%blocksize)numblocks++;
    printf("NUMBER OF BLOCKS: %d BUFFSIZE: %d\n",numblocks,buffsize);
    int *blocknz = new int[numblocks];
    #pragma omp parallel
    {
      int *numint = new int[numpix];
      #pragma omp for schedule(dynamic)
      for(int block = 0; block < numblocks; block++){
        for(int n = 0; n < numpix; n++)
          numint[n] = 0;
        for(int m = block*blocksize; m < (block+1)*blocksize && m < numray; m++)
          for(int n = rowdispl[m]; n < rowdispl[m+1]; n++)
            numint[rowindex[n]]++;
        int count = 0;
        for(int n = 0; n < numpix; n++)
          if(numint[n])count++;
        blocknz[block] = count/buffsize;
        if(count%buffsize)blocknz[block]++;
      }
      delete[] numint;
    }
    int *blockdispl = new int[numblocks+1];
    blockdispl[0] = 0;
    for(int block = 1; block < numblocks+1; block++)
      blockdispl[block] = blockdispl[block-1] + blocknz[block-1];
    int blocknzmax = 0;
    for(int block = 0; block < numblocks; block++)
      if(blocknz[block]>blocknzmax)blocknzmax = blocknz[block];
    delete[] blocknz;
    int blocknztot = blockdispl[numblocks];
    printf("NUMBER OF BUFFERS: %d AVERAGE BUFF/BLOCK: %f MAX BUFF/BLOCK: %d\n",blocknztot,blocknztot/(float)numblocks,blocknzmax);
    printf("BUFF MAP: %d (%f GB)\n",blocknztot*buffsize,blocknztot*buffsize/1024.0/1024.0/1024.0*sizeof(float));
    int *buffmap = new int[blocknztot*buffsize];
    int *buffnz = new int[blocknztot];
    long footprint = 0;
    #pragma omp parallel
    {
      int *numint = new int[numpix];
      int buffnztemp[blocknzmax];
      #pragma omp for
      for(int n = 0; n < blocknztot*buffsize; n++)
        buffmap[n] = 0;
      #pragma omp for
      for(int n = 0; n < blocknztot; n++)
        buffnz[n] = 0;
      #pragma omp for schedule(dynamic)
      for(int block = 0; block < numblocks; block++){
        for(int n = 0; n < numpix; n++)
          numint[n] = 0;
        for(int m = block*blocksize; m < (block+1)*blocksize && m < numray; m++)
          for(int n = rowdispl[m]; n < rowdispl[m+1]; n++)
            numint[rowindex[n]]++;
        int count = 0;
        for(int n = 0; n < numpix; n++)
          if(numint[n]){
            int buffloc = count/buffsize;
            int mapind = (blockdispl[block]+buffloc)*buffsize+count%buffsize;
            buffmap[mapind] = n;
            numint[n] = buffloc;
            count++;
          }
        #pragma omp atomic
        footprint += count;
        for(int m = block*blocksize; m < (block+1)*blocksize && m < numray; m++){
          for(int n = 0; n < blocknzmax; n++)
            buffnztemp[n] = 0;
          for(int n = rowdispl[m]; n < rowdispl[m+1]; n++)
            buffnztemp[numint[rowindex[n]]]++;
          for(int buff = blockdispl[block]; buff < blockdispl[block+1]; buff++){
            int buffloc = buff-blockdispl[block];
            if(buffnztemp[buffloc]>buffnz[buff])buffnz[buff]=buffnztemp[buffloc];
          }
        }
      }
      delete[] numint;
    }
    int *buffdispl = new int[blocknztot+1];
    buffdispl[0] = 0;
    for(int buff = 1; buff < blocknztot+1; buff++)
      buffdispl[buff] = buffdispl[buff-1]+buffnz[buff-1];
    int buffnzmax = 0;
    for(int buff = 0; buff < blocknztot; buff++)
      if(buffnzmax < buffnz[buff])buffnzmax = buffnz[buff];
    delete[] buffnz;
    int buffnztot = buffdispl[blocknztot];
    printf("ELLPACK STORAGE: %ld (%f GB) buffnzmax: %d STORAGE EFFICIENCY: %f DATA REUSE: %f\n",buffnztot*(long)blocksize,buffnztot*(float)blocksize*sizeof(float)*1.5/1024.0/1024.0/1024.0,buffnzmax,proj_rownztot/(float)buffnztot/blocksize*1.5,proj_rownztot/(float)footprint);
    short *buffindex = new short[buffnztot*blocksize];
    #pragma omp parallel
    {
      int *numint = new int[numpix];
      int *numind = new int[numpix];
      int buffnztemp[blocknzmax];
      #pragma omp for
      for(int n = 0; n < buffnztot*blocksize; n++)
        buffindex[n] = -1;
      #pragma omp for schedule(dynamic)
      for(int block = 0; block < numblocks; block++){
        for(int n = 0; n < numpix; n++)
          numint[n] = 0;
        for(int m = block*blocksize; m < (block+1)*blocksize && m < numray; m++)
          for(int n = rowdispl[m]; n < rowdispl[m+1]; n++)
            numint[rowindex[n]]++;
        int count = 0;
        for(int n = 0; n < numpix; n++)
          if(numint[n]){
            int buffloc = count/buffsize;
            numint[n] = buffloc;
            numind[n] = count%buffsize;
            count++;
          }
        for(int m = block*blocksize; m < (block+1)*blocksize && m < numray; m++){
          for(int n = 0; n < blocknzmax; n++)
            buffnztemp[n] = 0;
          for(int n = rowdispl[m]; n < rowdispl[m+1]; n++){
            int buffloc = numint[rowindex[n]];
            int row = buffdispl[blockdispl[block]+buffloc]+buffnztemp[buffloc];
            buffindex[row*blocksize+m%blocksize] = numind[rowindex[n]];
            buffnztemp[buffloc]++;
          }
        }
      }
      delete[] numint;
      delete[] numind;
    }
    delete[] rowindex;
    proj_numblocks = numblocks;
    proj_blocknztot = blocknztot;
    proj_blockdispl = blockdispl;
    proj_buffnztot = buffnztot;
    proj_buffdispl = buffdispl;
    proj_buffmap = buffmap;
    proj_buffindex = buffindex;
  }
  printf("BLOCKING TIME: %e\n",omp_get_wtime()-project_time);
  project_time = omp_get_wtime();
  printf("BLOCKING BACKPROJECTION MATRIX\n");
  {
    int *rowindex = back_rowindex;
    int *rowdispl = back_rowdispl;
    int blocksize = back_blocksize;
    int buffsize = back_buffsize;
    int numblocks = numpix/blocksize;
    if(numpix%blocksize)numblocks++;
    printf("NUMBER OF BLOCKS: %d BUFFSIZE: %d\n",numblocks,buffsize);
    int *blocknz = new int[numblocks];
    #pragma omp parallel
    {
      int *numint = new int[numray];
      #pragma omp for schedule(dynamic)
      for(int block = 0; block < numblocks; block++){
        for(int n = 0; n < numray; n++)
          numint[n] = 0;
        for(int m = block*blocksize; m < (block+1)*blocksize && m < numpix; m++)
          for(int n = rowdispl[m]; n < rowdispl[m+1]; n++)
            numint[rowindex[n]]++;
        int count = 0;
        for(int n = 0; n < numray; n++)
          if(numint[n])count++;
        blocknz[block] = count/buffsize;
        if(count%buffsize)blocknz[block]++;
      }
      delete[] numint;
    }
    int *blockdispl = new int[numblocks+1];
    blockdispl[0] = 0;
    for(int block = 1; block < numblocks+1; block++)
      blockdispl[block] = blockdispl[block-1] + blocknz[block-1];
    int blocknzmax = 0;
    for(int block = 0; block < numblocks; block++)
      if(blocknz[block]>blocknzmax)blocknzmax = blocknz[block];
    delete[] blocknz;
    int blocknztot = blockdispl[numblocks];
    printf("NUMBER OF BUFFERS: %d AVERAGE BUFF/BLOCK: %f MAX BUFF/BLOCK: %d\n",blocknztot,blocknztot/(float)numblocks,blocknzmax);
    printf("BUFF MAP: %d (%f GB)\n",blocknztot*buffsize,blocknztot*buffsize/1024.0/1024.0/1024.0*sizeof(float));
    int *buffmap = new int[blocknztot*buffsize];
    int *buffnz = new int[blocknztot];
    long footprint = 0;
    #pragma omp parallel
    {
      int *numint = new int[numray];
      int buffnztemp[blocknzmax];
      #pragma omp for
      for(int n = 0; n < blocknztot*buffsize; n++)
        buffmap[n] = 0;
      #pragma omp for
      for(int n = 0; n < blocknztot; n++)
        buffnz[n] = 0;
      #pragma omp for schedule(dynamic)
      for(int block = 0; block < numblocks; block++){
        for(int n = 0; n < numray; n++)
          numint[n] = 0;
        for(int m = block*blocksize; m < (block+1)*blocksize && m < numpix; m++)
          for(int n = rowdispl[m]; n < rowdispl[m+1]; n++)
            numint[rowindex[n]]++;
        int count = 0;
        for(int n = 0; n < numray; n++)
          if(numint[n]){
            int buffloc = count/buffsize;
            int mapind = (blockdispl[block]+buffloc)*buffsize+count%buffsize;
            buffmap[mapind] = n;
            numint[n] = buffloc;
            count++;
          }
        #pragma omp atomic
        footprint += count;
        for(int m = block*blocksize; m < (block+1)*blocksize && m < numpix; m++){
          for(int n = 0; n < blocknzmax; n++)
            buffnztemp[n] = 0;
          for(int n = rowdispl[m]; n < rowdispl[m+1]; n++)
            buffnztemp[numint[rowindex[n]]]++;
          for(int buff = blockdispl[block]; buff < blockdispl[block+1]; buff++){
            int buffloc = buff-blockdispl[block];
            if(buffnztemp[buffloc]>buffnz[buff])buffnz[buff]=buffnztemp[buffloc];
          }
        }
      }
      delete[] numint;
    }
    int *buffdispl = new int[blocknztot+1];
    buffdispl[0] = 0;
    for(int buff = 1; buff < blocknztot+1; buff++)
      buffdispl[buff] = buffdispl[buff-1]+buffnz[buff-1];
    int buffnzmax = 0;
    for(int buff = 0; buff < blocknztot; buff++)
      if(buffnzmax < buffnz[buff])buffnzmax = buffnz[buff];
    delete[] buffnz;
    int buffnztot = buffdispl[blocknztot];
    printf("ELLPACK STORAGE: %d (%f GB) buffnzmax: %d STORAGE EFFICIENCY: %f DATA REUSE: %f\n",buffnztot*blocksize,buffnztot*(float)blocksize*sizeof(float)*1.5/1024.0/1024.0/1024.0,buffnzmax,proj_rownztot/(float)buffnztot/blocksize*1.5,proj_rownztot/(float)footprint);
    short *buffindex = new short[buffnztot*blocksize];
    #pragma omp parallel
    {
      int *numint = new int[numray];
      int *numind = new int[numray];
      int buffnztemp[blocknzmax];
      #pragma omp for
      for(int n = 0; n < buffnztot*blocksize; n++)
        buffindex[n] = -1;
      #pragma omp for schedule(dynamic)
      for(int block = 0; block < numblocks; block++){
        for(int n = 0; n < numray; n++)
          numint[n] = 0;
        for(int m = block*blocksize; m < (block+1)*blocksize && m < numpix; m++)
          for(int n = rowdispl[m]; n < rowdispl[m+1]; n++)
            numint[rowindex[n]]++;
        int count = 0;
        for(int n = 0; n < numray; n++)
          if(numint[n]){
            int buffloc = count/buffsize;
            numint[n] = buffloc;
            numind[n] = count%buffsize;
            count++;
          }
        for(int m = block*blocksize; m < (block+1)*blocksize && m < numpix; m++){
          for(int n = 0; n < blocknzmax; n++)
            buffnztemp[n] = 0;
          for(int n = rowdispl[m]; n < rowdispl[m+1]; n++){
            int buffloc = numint[rowindex[n]];
            int row = buffdispl[blockdispl[block]+buffloc]+buffnztemp[buffloc];
            buffindex[row*blocksize+m%blocksize] = numind[rowindex[n]];
            buffnztemp[buffloc]++;
          }
        }
      }
      delete[] numint;
      delete[] numind;
    }
    delete[] rowindex;
    back_numblocks = numblocks;
    back_blocknztot = blocknztot;
    back_blockdispl = blockdispl;
    back_buffnztot = buffnztot;
    back_buffdispl = buffdispl;
    back_buffmap = buffmap;
    back_buffindex = buffindex;
  }
  printf("BLOCKING TIME: %e\n",omp_get_wtime()-project_time);
  project_time = omp_get_wtime();
  printf("\nFILL PROJECTION MATRIX\n");
  {
    proj_buffvalue = new float[proj_buffnztot*proj_blocksize];
    #pragma omp parallel for
    for(int n = 0; n < proj_buffnztot*proj_blocksize; n++)
      proj_buffvalue[n] = 0;
    #pragma omp parallel for schedule(dynamic)
    for(int block = 0; block < proj_numblocks; block++)
      for(int k = block*proj_blocksize; k < (block+1)*proj_blocksize && k < numray; k++){
        float rho = raycoor[k].real();
        float theta = raycoor[k].imag();
        for(int buff = proj_blockdispl[block]; buff < proj_blockdispl[block+1]; buff++)
          for(int row = proj_buffdispl[buff]; row < proj_buffdispl[buff+1]; row++){
            int ind = row*proj_blocksize+k%proj_blocksize;
            if(proj_buffindex[ind]==-1)
              proj_buffindex[ind] = 0;
            else{
              int mapind = buff*proj_buffsize+proj_buffindex[ind];
              int pixind = proj_buffmap[mapind];
              float domain[4];
              domain[0]=pixcoor[pixind].real()-pixsize/2;
              domain[1]=domain[0]+pixsize;
              domain[2]=pixcoor[pixind].imag()-pixsize/2;
              domain[3]=domain[2]+pixsize;
              findlength(theta,rho,&domain[0],&proj_buffvalue[ind]);
            }
          }
      }
    #pragma omp parallel for
    for(int n = 0; n < proj_buffnztot*proj_blocksize; n++)
      if(proj_buffindex[n]==-1){
        proj_buffindex[n] = 0;
        proj_buffvalue[n] = 0;
      }
  }
  printf("TIME: %e\n",omp_get_wtime()-project_time);
  printf("FILL BACKPROJECTION MATRIX\n");
  {
    back_buffvalue = new float[back_buffnztot*back_blocksize];
    #pragma omp parallel for
    for(int n = 0; n < back_buffnztot*back_blocksize; n++)
      back_buffvalue[n] = 0;
    #pragma omp parallel for schedule(dynamic)
    for(int block = 0; block < back_numblocks; block++)
      for(int n = block*back_blocksize; n < (block+1)*back_blocksize && n < numpix; n++){
        float domain[4];
        domain[0]=pixcoor[n].real()-pixsize/2;
        domain[1]=domain[0]+pixsize;
        domain[2]=pixcoor[n].imag()-pixsize/2;
        domain[3]=domain[2]+pixsize;
        for(int buff = back_blockdispl[block]; buff < back_blockdispl[block+1]; buff++)
          for(int row = back_buffdispl[buff]; row < back_buffdispl[buff+1]; row++){
            int ind = row*back_blocksize+n%back_blocksize;
            if(back_buffindex[ind]==-1)
              back_buffindex[ind] = 0;
            else{
              int mapind = buff*back_buffsize+back_buffindex[ind];
              int rayind = back_buffmap[mapind];
              float rho = raycoor[rayind].real();
              float theta = raycoor[rayind].imag();
              findlength(theta,rho,&domain[0],&back_buffvalue[ind]);
            }
          }
      }
    #pragma omp parallel for
    for(int n = 0; n < back_buffnztot*back_blocksize; n++)
      if(back_buffindex[n]==-1){
        back_buffindex[n] = 0;
        back_buffvalue[n] = 0;
      }
  }
  printf("TIME: %e\n",omp_get_wtime()-project_time);
  printf("PREPROCESSING TIME: %e\n",omp_get_wtime()-timep);

  float *obj;
  float *mes;
  float *ray;
  float *res;
  float *gra;
  float *dir;

  setup_gpu(&obj,&gra,&dir,&mes,&res,&ray);

  numproj = 0;
  numback = 0;

  printf("INPUT MEASUREMENT DATA\n");
  float *mesdata = new float[numrho*numthe];
  FILE *dataf = fopen(sinfile,"rb");
  fread(mesdata,sizeof(float),numrho*numthe,dataf);
  fclose(dataf);
  #pragma omp parallel for
  for(int k = 0; k < numray; k++)
    if(raymesind[k]>-1)mes[k] = mesdata[raymesind[k]];
    else mes[k] = 0;
  delete[] mesdata;
  printf("INPUT ENDS\n");
  delete[] rayglobalind;
  delete[] raymesind;

  printf("GRADIENT-DESCENT OPTIMIZATION\n");
  //FILE *resf = fopen("residual.txt","w");
  //INITIALIZE OBJECT
  #pragma omp parallel for
  for(int n = 0; n < numpix; n++)
    obj[n] = 0;
  {
    double time;
    double ctime = 0;
    double wtime = 0;
    double rtime = omp_get_wtime();
    //FORWARD PROJECTION
    projection(ray,obj);
    //FIND RESIDUAL ERROR
    time = omp_get_wtime();
    subtract_kernel(res,ray,mes,numray);
    ctime = ctime + omp_get_wtime() - time;
    //FIND GRADIENT
    backprojection(gra,res);
    time = omp_get_wtime();
    float error = norm_kernel(res,numray);
    float gradnorm = norm_kernel(gra,numpix);
    printf("iter: %d error: %e gradnorm: %e\n",0,error,gradnorm);
    printf("Before save direction.\n");
    fflush(stdout);
    //SAVE DIRECTION
    copy_kernel(dir,gra,numpix);
    float oldgradnorm = gradnorm;
    ctime = ctime + omp_get_wtime() - time;
    //START ITERATIONS
    for(int iter = 1; iter <= numiter; iter++){
      //PROJECT DIRECTION
      projection(ray,dir);
      //FIND STEP SIZE
      time = omp_get_wtime();
      float temp1 = dot_kernel(res,ray,numray);
      float temp2 = norm_kernel(ray,numray);
      //STEP SIZE
      float alpha = temp1/temp2;
      saxpy_kernel(obj,obj,-1.0*alpha,dir,numpix);
      ctime = ctime + omp_get_wtime() - time;
      //FORWARD PROJECTION
      projection(ray,obj);
      //FIND RESIDUAL ERROR
      time = omp_get_wtime();
      subtract_kernel(res,ray,mes,numray);
      ctime = ctime + omp_get_wtime() - time;
      //FIND GRADIENT
      backprojection(gra,res);
      time = omp_get_wtime();
      float error = norm_kernel(res,numray);
      float gradnorm = norm_kernel(gra,numpix);
      //UPDATE DIRECTION
      printf("iter: %d error: %e gradnorm: %e\n",iter,error,gradnorm);
      printf("Before sxay kernel\n");
      fflush(stdout);
      //fprintf(resf,"%e %e\n",error,gradnorm);
      float beta = gradnorm/oldgradnorm;
      //float beta = 0;
      oldgradnorm = gradnorm;
      saxpy_kernel(dir,gra,beta,dir,numpix);
      ctime = ctime + omp_get_wtime() - time;
    }
    rtime = omp_get_wtime()-rtime;
    //fclose(resf);

    printf("recon: %e proj: %e (%e %e %e) backproj: %e (%e %e %e) ctime: %e wtime: %e\n",rtime,ftime,fktime,aftime,frtime,btime,brtime,abtime,bktime,ctime,wtime);

    printf("numproj: %d numback: %d\n",numproj,numback);
    double projflops = proj_rownztot/1e9*2*numproj/fktime;
    double backflops = proj_rownztot/1e9*2*numback/bktime;
    double totflops = proj_rownztot/1e9*2*(numproj+numback)/(fktime+bktime);
    double projbw = (proj_buffnztot*proj_blocksize/1e9*6+proj_blocknztot*proj_buffsize/1e9*4)*numproj/fktime;
    double backbw = (back_buffnztot*back_blocksize/1e9*6+back_blocknztot*back_buffsize/1e9*4)*numback/bktime;
    double totbw = ((proj_buffnztot*proj_blocksize/1e9*6+proj_blocknztot*proj_buffsize/1e9*4)*numproj+(back_buffnztot*back_blocksize/1e9*6+back_blocknztot*back_buffsize/1e9*4)*numback)/(fktime+bktime);
    printf("proj: %e s (%f GFLOPS) backproj: %e s (%f GFLOPS) GFLOPS: %f\n",fktime,projflops,bktime,backflops,totflops);
    printf("proj: %f GB/s back: %f GB/s %f GB/s\n",projbw,backbw,totbw);

  }
  float *objtemp = new float[numpix];
  #pragma omp parallel for
  for(int n = 0; n < numpix; n++)
    objtemp[pixglobalind[n]] = obj[n];
  FILE *objf = fopen(outfile,"wb");
  fwrite(objtemp,sizeof(float),numpix,objf);
  fclose(objf);
  delete[] objtemp;

  printf("Total Time: %e\n",omp_get_wtime()-timetot);

  return 0;
}

float norm_kernel(float *a, int dim){
  float reduce = 0;
  #pragma omp parallel for reduction(+:reduce)
  for(int n = 0; n < dim; n++)
    reduce += norm(a[n]);
  return reduce;
};
float dot_kernel(float *a, float *b, int dim){
  float reduce = 0;
  #pragma omp parallel for reduction(+:reduce)
  for(int n = 0; n < dim; n++)
    reduce += a[n]*b[n];
  return reduce;
};
void saxpy_kernel(float *a, float *b, float coef, float *c, int dim){
  #pragma omp parallel for
  for(int n = 0; n < dim; n++)
    a[n] = b[n] + coef*c[n];
};
void copy_kernel(float *a, float *b, int dim){
  #pragma omp parallel for
  for(int n = 0; n < dim; n++)
    a[n] = b[n];
};
void subtract_kernel(float *a, float *b, float *c, int dim){
  #pragma omp parallel for
  for(int n = 0; n < dim; n++)
    a[n] = b[n] - c[n];
};
