    // OpenCL Kernel Function for element by element vector multiply
    kernel void VectorMultiply(global double* data, global const double* recip, int numElements) {

        // get index into global data array
        int iGID = get_global_id(0);

        // bound check, equivalent to the limit on a 'for' loop
        if (iGID >= numElements)  {
            return;
        }

        // Element wise multiply.
        int i2 = iGID * 2;
        data[i2] = data[i2] * recip[iGID];
        data[i2 + 1] = data[i2 + 1] * recip[iGID];
    }
