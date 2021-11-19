#include <cstdlib>
#include <cstdio>
#include <string.h>
#include <helper_cuda.h>

#include <math.h>

// cube length
__device__ __constant__ float L;

// whether minimum image convention periodic boundary conditions apply
__device__ __constant__ float periodic;

// for the radial distribution function
__device__ __constant__ float dr2;

// pressure via the virial theorem
__device__ float pressure;

// tau_xy of stress-eneryg tensor via virial theorem
__device__ float Pxy;

//__device__ float4
__device__ void
bodyBodyInteraction(float4 &ai, float4 bi, float4 &bj, int *RDF)
{
  float3 r;

  r.x = bj.x - bi.x;
  r.y = bj.y - bi.y;
  r.z = bj.z - bi.z;

  // Minimum image convention
  r.x += -periodic * L * lround(r.x / L);
  r.y += -periodic * L * lround(r.y / L);
  r.z += -periodic * L * lround(r.z / L);

  float distSqr = r.x * r.x + r.y * r.y + r.z * r.z;

  if (distSqr > 1e-6f)
  {
    int indr2 = min(__float2int_rd(distSqr / dr2), 255);

    if (RDF != NULL)
      atomicAdd(&RDF[indr2], 1);

    float invDist2 = 1.0f / distSqr;
    float invDist6 = invDist2 * invDist2 * invDist2;

    float s = invDist2 * (12.0f * invDist6 * invDist6 - 6.0f * invDist6);

    ai.x += r.x * s;
    ai.y += r.y * s;
    ai.z += r.z * s;
    ai.w += invDist6 * invDist6 - invDist6;
    bj.w += s * distSqr;

    atomicAdd(&Pxy, -r.x * r.y * s);
  }
}


// This is the "tile_calculation" function from the GPUG3 article.
__device__ float4 tile_force(float4 &myPos, float4 accel, int maxnumb, int *RDF)
{
  extern __shared__ float4 sharedPos[];
  int i;

#pragma unroll 8
  for (i = 0; i < maxnumb; ++i)
  {
    bodyBodyInteraction(accel, sharedPos[i], myPos, RDF);
  }

  return accel;
}

__device__ float4
computeBodyAccel(float4 &bodyPos, float4* positions, int numBodies, int maxnumblast, int *RDF)
{
  extern __shared__ float4 sharedPos[];

  float4 frc = { 0.0f, 0.0f, 0.0f, 0.0f };

  int start = 0;
  int tile = 0;
  int finish = start + numBodies - blockDim.x;

  bodyPos.w = 0.f;

  for (int i = start; i < finish; i += blockDim.x, tile++)
  {
    sharedPos[threadIdx.x] = positions[tile * blockDim.x + threadIdx.x];

    __syncthreads();
    frc = tile_force(bodyPos, frc, blockDim.x, RDF);
    __syncthreads();
  }

  if (threadIdx.x < maxnumblast)
    sharedPos[threadIdx.x] = positions[tile * blockDim.x + threadIdx.x];

  __syncthreads();

  frc = tile_force(bodyPos, frc, maxnumblast, RDF);
  __syncthreads();

  // Lennard-Jones factor
  frc.x *= 4.f;
  frc.y *= 4.f;
  frc.z *= 4.f;
  bodyPos.w *= 4.0f / 2.0f / 3.0f;

  return frc;
}

__global__ void
calculateForces(float4* Pos, float4* Force,
  int numBodies, int calcmaxnumblast, int *RDF)
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;

  int g = blockIdx.x;
  int* gmem = RDF + g * 256;
  for (int i = threadIdx.x; i < 256; i += blockDim.x)
    gmem[i] = 0;
  __syncthreads();

  int index2 = index;
  if (index2 >= numBodies) {
    index2 = numBodies - 1;
    gmem = NULL;
  }

  float4 pos = Pos[index2];
  float4 accel = computeBodyAccel(pos, Pos, numBodies, calcmaxnumblast, gmem);

  if (index < numBodies) {
    Force[index2] = accel;
    atomicAdd(&pressure, pos.w);
  }
}

__global__ void
RDFmerger(int* RDF, int blocks, int *RDFout)
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index >= 256)
    return;

  int total = 0;
  for (int j = 0; j < blocks; ++j) {
    total += RDF[index + blockDim.x * j];
  }

  RDFout[index] = total;
}

extern "C"
{

  void allocateNBodyArrays(float* vel[2], int numBodies)
  {
    // 4 floats each for alignment reasons
    unsigned int memSize = sizeof(float) * 4 * numBodies;

    checkCudaErrors(cudaMalloc((void**)&vel[0], memSize));
    checkCudaErrors(cudaMalloc((void**)&vel[1], memSize));
  }

  void allocateArray(float** dest, int number)
  {
    // 4 floats each for alignment reasons
    unsigned int memSize = sizeof(float) * 4 * number;

    checkCudaErrors(cudaMalloc((void**)dest, memSize));
  }

  void deleteNBodyArrays(float* vel[2])
  {
    checkCudaErrors(cudaFree((void**)vel[0]));
    checkCudaErrors(cudaFree((void**)vel[1]));
  }

  void deleteArray(float* arr)
  {
    checkCudaErrors(cudaFree((void**)arr));
  }

  void allocateIntArray(int** dest, int number)
  {
    // 4 floats each for alignment reasons
    unsigned int memSize = sizeof(int) * number;

    checkCudaErrors(cudaMalloc((void**)dest, memSize));
  }

  void deleteIntArray(int* arr)
  {
    checkCudaErrors(cudaFree((void**)arr));
  }

  void copyArrayFromDevice(float* host,
    const float* device,
    unsigned int pbo,
    int numBodies)
  {
    //if (pbo)
    //    checkCudaErrors(cudaGLMapBufferObject((void**)&device, pbo));
    checkCudaErrors(cudaMemcpy(host, device, numBodies * 4 * sizeof(float),
      cudaMemcpyDeviceToHost));
    //if (pbo)
    //   checkCudaErrors(cudaGLUnmapBufferObject(pbo));
  }

  void copyArrayToDevice(float* device, const float* host, int numBodies)
  {
    checkCudaErrors(cudaMemcpy(device, host, numBodies * 4 * sizeof(float),
      cudaMemcpyHostToDevice));
  }

  void registerGLBufferObject(unsigned int pbo)
  {
    //checkCudaErrors(cudaGLRegisterBufferObject(pbo));
  }

  void unregisterGLBufferObject(unsigned int pbo)
  {
    //checkCudaErrors(cudaGLUnregisterBufferObject(pbo));
  }

  void threadSync() { cudaThreadSynchronize(); }

  void
    calculateNForces(float* Pos, float* Force, float* host_pressure, float* host_Pxy,
      int numBodies, float host_L, int Lperiodic, int *host_RDF, float host_dr2, int blockSize, int q)
  {
    int sharedMemSize = blockSize * sizeof(float4);
    int gridSize = (int)((numBodies + blockSize - 1) / blockSize);
    dim3 dimGrid(gridSize);
    dim3 dimBlock(blockSize);
    int calcmaxnumblast = numBodies % blockSize;
    calcmaxnumblast = (calcmaxnumblast == 0) ? blockSize : calcmaxnumblast;

    checkCudaErrors(cudaMemcpyToSymbol(L,
      &host_L,
      sizeof(float), 0,
      cudaMemcpyHostToDevice));

    float host_periodic = static_cast<float>(Lperiodic);

    checkCudaErrors(cudaMemcpyToSymbol(periodic,
      &host_periodic,
      sizeof(float), 0,
      cudaMemcpyHostToDevice));

    *host_pressure = 0.f;
    *host_Pxy = 0.f;

    checkCudaErrors(cudaMemcpyToSymbol(pressure,
      host_pressure,
      sizeof(float), 0,
      cudaMemcpyHostToDevice));

    checkCudaErrors(cudaMemcpyToSymbol(Pxy,
      host_Pxy,
      sizeof(float), 0,
      cudaMemcpyHostToDevice));

    checkCudaErrors(cudaMemcpyToSymbol(dr2,
      &host_dr2,
      sizeof(float), 0,
      cudaMemcpyHostToDevice));

    int* RDFall;
    allocateIntArray(&RDFall, 256 * gridSize);

    calculateForces << < dimGrid, dimBlock, sharedMemSize >> >
      ((float4*)Pos, (float4*)Force,
        numBodies, calcmaxnumblast, RDFall);

    checkCudaErrors(cudaMemcpyFromSymbol(host_pressure,
      pressure,
      sizeof(float), 0,
      cudaMemcpyDeviceToHost));

    checkCudaErrors(cudaMemcpyFromSymbol(host_Pxy,
      Pxy,
      sizeof(float), 0,
      cudaMemcpyDeviceToHost));

    int* RDFtot;
    allocateIntArray(&RDFtot, 256);

    RDFmerger << <   dimGrid, dimBlock, sharedMemSize >> >
      (RDFall, gridSize, RDFtot);

    checkCudaErrors(cudaMemcpy(host_RDF, RDFtot, 256 * sizeof(int),
      cudaMemcpyDeviceToHost));

    deleteIntArray(RDFall);
    deleteIntArray(RDFtot);

    // check if kernel invocation generated an error
    getLastCudaError("Kernel execution failed");
  }


  void threadExit()
  {
    cudaThreadExit();
  }

}
