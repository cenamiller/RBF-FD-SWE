// Written by Richelle Streater, Summer 2018.
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel __attribute__((vec_type_hint(double4)))
__kernel void eval_RHS_kernel(  __global double4* K, __global const int* idx, __global const double* Dx, __global const double* Dy,
                                __global const double* Dz, __global const double* L, __global const double* x, __global const double* y,
                                __global const double* z, __global const double* f, __global const double* ghm, __global const double* p_u,
                                __global const double* p_v, __global const double* p_w, __global const double* gradghm, __global double4* F,
                                const double gh0, const int padded_Nnodes, const int padded_Nnbr, const int Nnbr, const int compute_pid_s) {
    int i = get_global_id(0);
    if (i < padded_Nnodes) {
    
        double4 H;
        double4 rhs;
    
        // Four-vectors represent u,v,w,h components of the state variable matrices
        double4 dH_dx = (double4)(0.0);
        double4 dH_dy = (double4)(0.0);
        double4 dH_dz = (double4)(0.0);
        double4 LH = (double4)(0.0);

        #pragma unroll
        for (int j = 0; j < Nnbr; j++){
            int dmid = i*(padded_Nnbr) + j; // DM id = padded_Nnbr*i + j
            H = K[idx[dmid]];
        
            // Update sums
            dH_dx += Dx[dmid] * H; // dH_m/d_n = (D^n)_(ij) (H_m)^(ij)
            dH_dy += Dy[dmid] * H;
            dH_dz += Dz[dmid] * H;
            LH += L[dmid] * H;
        }

        // Reset to get H at current node (rather than final neighbor)
        H = K[i+compute_pid_s];
        
        // Evaluate RHS of unconstrained momentum equations
        rhs = (H.s0 * dH_dx) + (H.s1 * dH_dy) + (H.s2 * dH_dz);
        rhs.s0 += f[i] * ((y[i] * H.s2) - (z[i] * H.s1)) + dH_dx.s3;
        rhs.s1 += f[i] * ((z[i] * H.s0) - (x[i] * H.s2)) + dH_dy.s3;
        rhs.s2 += f[i] * ((x[i] * H.s1) - (y[i] * H.s0)) + dH_dz.s3;
        rhs *= -1;
    
        // Project momentum equations to confine motion to the surface of the sphere and add hyperviscosity
        // Get indices of x, y, and z components of projection matrices
        int x_id = 3*i;
        int y_id = 3*i + 1;
        int z_id = 3*i + 2;
    
        // Evaluate projections and apply hyperviscosity
        F[i].s0 = (p_u[x_id] * rhs.s0) + (p_u[y_id] * rhs.s1) + (p_u[z_id] * rhs.s2);
        F[i].s1 = (p_v[x_id] * rhs.s0) + (p_v[y_id] * rhs.s1) + (p_v[z_id] * rhs.s2);
        F[i].s2 = (p_w[x_id] * rhs.s0) + (p_w[y_id] * rhs.s1) + (p_w[z_id] * rhs.s2);
    
        // Calculate the RHS of the geopotential equation
        F[i].s3 = rhs.s3 + (H.s0 * gradghm[x_id]) + (H.s1 * gradghm[y_id]) + (H.s2 * gradghm[z_id]) - (H.s3 + gh0 - ghm[i]) * (dH_dx.s0 + dH_dy.s1 + dH_dz.s2);
        F[i] += LH;
    }
    return;
}

__kernel void copy_arr_kernel(__global double* D, __global const double* F, const int arr_size) {
    int i = get_global_id(0);
    if (i >= arr_size) return;
        
    D[i] = F[i]; // D = F
    
    return;
}
    
__kernel void update_D_kernel(__global const double* F, __global double* D, const double coefficient, const int Nnodes) {
    int i = get_global_id(0);
    if (i >= Nnodes*4) return;

    D[i] += coefficient * F[i];  // D = c*F

    return;
}

__kernel void eval_K_kernel(__global double* H, __global const double* F, __global double* K,
                            const double dt, const int Nnodes, const int compute_pid_s) {
    int i = get_global_id(0);
    if (i >= Nnodes*4) return;

    K[i+compute_pid_s*4] = H[i+compute_pid_s*4] + dt * F[i];        // K = H + (coeff)*dt*F (RK coefficient handled in invocation (RK_ocl.c)

    return;
}

__kernel void update_H_kernel(	__global double* H, const double dt, __global double* D,
                                __global double* K, const int Nnodes, const int compute_pid_s) {
    int i = get_global_id(0);
    if (i >= Nnodes*4) return;
	
    H[i+compute_pid_s*4] += (dt/6.0) * D[i];        // H += dt/6*D (Add all of the RK pieces and multiply by 1/6 per RK 4th order algorithm)
    K[i+compute_pid_s*4] = H[i+compute_pid_s*4];    // Update to use H in the first round of the next time step w/o updating eval_RHS args
    return;
}
