extern "C"
__global__ void recipSummation(float* data, float* recip, int len)
{
    const int y = blockIdx.y * gridDim.x * blockDim.x;
    const int x = blockIdx.x * blockDim.x;
    const int i = threadIdx.x + x + y;
    if (i < len) {
       const int j = 2 * i;
       data[j]     *= recip[i];
       data[j + 1] *= recip[i];
    }
}
