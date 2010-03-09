extern "C"
__global__ void recipSummation(float* data, float* recip, int len)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < len) {
       const int j = 2 * i;
       data[j]     *= recip[i];
       data[j + 1] *= recip[i];
    }
}
