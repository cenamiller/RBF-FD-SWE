// Written by Richelle Streater, Summer 2018.

#pragma OPENCL EXTENSION cl_khr_fp64 : enable

#define ROUNDDOWN_SIMD(N) (((N) / SIMD_LENGTH) * SIMD_LENGTH)

#if ((SIMD_LENGTH == 4) || (SIMD_LENGTH == 8)) // Set up to group with double4 or double8. See OpenCL vector types for full list of possible vectors.
#define TILE_LENGTH SIMD_LENGTH
#else
#define TILE_LENGTH 1
#endif

#ifdef USE_CFDL // If uvwh is the fast index (USE_CFDL), access normally. Otherwise, use formula to access 3d array with two variables.
#define SVM_LIN_ID(PID, SV_ID) ((4 * (PID)) + (SV_ID))
#else
#define SVM_LIN_ID(PID, SV_ID) ((ROUNDDOWN_SIMD(PID) * 4) + ((SV_ID) * SIMD_LENGTH) + ((PID) % SIMD_LENGTH))
#endif

#if !defined(USE_CFDL) && (TILE_LENGTH > 1) // If fast index on SVM matrices is nodes direction, group incoming buffers by 4 or 8
#if (TILE_LENGTH == 4)
typedef double4 vec_f_svm;
#else
typedef double8 vec_f_svm;
#endif
#else
typedef double vec_f_svm;
#endif

#if (TILE_LENGTH == 1) // Group private variables with double4 or double8, and create ACCESS function to access vectors in a loop.
typedef double vec_f;
typedef int vec_i;
#define ACCESS(VEC, E) (VEC)
#elif (TILE_LENGTH == 4)
typedef double4 vec_f;
typedef int4 vec_i;
#define ACCESS(VEC, E) ((E==0) ? ((VEC).s0): (E==1) ? ((VEC).s1): (E==2) ? ((VEC).s2): (E==3) ? ((VEC).s3): 0)
#else
typedef double8 vec_f;
typedef int8 vec_i;
#define ACCESS(VEC, E) ((E==0) ? ((VEC).s0): (E==1) ? ((VEC).s1): (E==2) ? ((VEC).s2): (E==3) ? ((VEC).s3): (E==4) ? ((VEC).s4): (E==5) ? ((VEC).s5): (E==6) ? ((VEC).s6): (E==7) ? ((VEC).s7): 0)
#endif

#if defined(USE_SFDL) && (TILE_LENGTH>1) // If fast index on diff. matrices is nodes direction, group incoming buffers by 4 or 8
#if (TILE_LENGTH == 4)
typedef double4 vec_f_dm;
typedef int4 vec_i_dm;
#else
typedef double8 vec_f_dm;
typedef int8 vec_i_dm;
#endif
#else
typedef double vec_f_dm;
typedef int vec_i_dm;
#endif

// Set accessing functions- if incoming arrays are grouped, 3d arrays act like 2d arrays- but if they are tiled and cannot be
// grouped (SIMD_LENGTH is not 4 or 8), it is necessary to use formula for accessing 3d array with 2 variables.
#if !(defined(USE_SFDL) && (TILE_LENGTH == 1))
#define DM_LIN_ID(NODE_ID, NBR_ID, NNBR) (((NODE_ID) * (NNBR)) + (NBR_ID))
#define VECT_LIN_ID(NODE_ID, COMP_ID) (((NODE_ID) * 3) + (COMP_ID))
#else
#define DM_LIN_ID(NODE_ID, NBR_ID, NNBR) ((ROUNDDOWN_SIMD(NODE_ID) * (NNBR)) + ((NBR_ID) * SIMD_LENGTH) + ((NODE_ID) % SIMD_LENGTH))
#define VECT_LIN_ID(NODE_ID, COMP_ID) ((ROUNDDOWN_SIMD(NODE_ID) * 3) + ((COMP_ID) * SIMD_LENGTH) + ((NODE_ID) % SIMD_LENGTH))
#endif

#if !defined(USE_CFDL) && (TILE_LENGTH == 1)
#define F_ID(PID, SV_ID) ((ROUNDDOWN_SIMD(PID) * 4) + ((SV_ID) * SIMD_LENGTH) + ((PID) % SIMD_LENGTH))
#else
#define F_ID(PID, SV_ID) ((4 * (PID)) + (SV_ID))
#endif

__kernel __attribute__((vec_type_hint(vec_f)))
__kernel void eval_RHS_kernel(  __global const double* K, __global const vec_i_dm* idx, __global const vec_f_dm* Dx, __global const vec_f_dm* Dy,
                                __global const vec_f_dm* Dz, __global const vec_f_dm* L, __global const vec_f_dm* x, __global const vec_f_dm* y,
                                __global const vec_f_dm* z, __global const vec_f_dm* f, __global const vec_f_dm* ghm, __global const vec_f_dm* p_u,
                                __global const vec_f_dm* p_v, __global const vec_f_dm* p_w, __global const vec_f_dm* gradghm, __global vec_f_svm* F,
                                const double gh0, const int padded_Nnodes, const int padded_Nnbr, const int Nnbr, const int compute_pid_s) {
    int i = get_global_id(0);
    if (i >= padded_Nnodes/TILE_LENGTH) return;
    
    union {int arr[TILE_LENGTH]; vec_i vec;} nbr_pid; // union allows accessing vector elements in "for" loop
    union {double arr[TILE_LENGTH]; vec_f vec;} u;
    union {double arr[TILE_LENGTH]; vec_f vec;} v;
    union {double arr[TILE_LENGTH]; vec_f vec;} w;
    union {double arr[TILE_LENGTH]; vec_f vec;} h;
    
    union {double arr[TILE_LENGTH]; vec_f vec;} Dx_temp;
    union {double arr[TILE_LENGTH]; vec_f vec;} Dy_temp;
    union {double arr[TILE_LENGTH]; vec_f vec;} Dz_temp;
    union {double arr[TILE_LENGTH]; vec_f vec;} L_temp;
    
    union {double arr[TILE_LENGTH]; vec_f vec;} x_temp;
    union {double arr[TILE_LENGTH]; vec_f vec;} y_temp;
    union {double arr[TILE_LENGTH]; vec_f vec;} z_temp;
    union {double arr[TILE_LENGTH]; vec_f vec;} f_temp;
    union {double arr[TILE_LENGTH]; vec_f vec;} ghm_temp;
    
    union {double arr[3][TILE_LENGTH]; vec_f vec[3];} p_u_temp;
    union {double arr[3][TILE_LENGTH]; vec_f vec[3];} p_v_temp;
    union {double arr[3][TILE_LENGTH]; vec_f vec[3];} p_w_temp;
    union {double arr[3][TILE_LENGTH]; vec_f vec[3];} gradghm_temp;
    
    vec_f du_dx = (vec_f)(0.0);
    vec_f dv_dx = (vec_f)(0.0);
    vec_f dw_dx = (vec_f)(0.0);
    vec_f dh_dx = (vec_f)(0.0);

    vec_f du_dy = (vec_f)(0.0);
    vec_f dv_dy = (vec_f)(0.0);
    vec_f dw_dy = (vec_f)(0.0);
    vec_f dh_dy = (vec_f)(0.0);
    
    vec_f du_dz = (vec_f)(0.0);
    vec_f dv_dz = (vec_f)(0.0);
    vec_f dw_dz = (vec_f)(0.0);
    vec_f dh_dz = (vec_f)(0.0);
    
    vec_f Lu = (vec_f)(0.0);
    vec_f Lv = (vec_f)(0.0);
    vec_f Lw = (vec_f)(0.0);
    vec_f Lh = (vec_f)(0.0);
    
    #pragma unroll
    for (int j = 0; j < Nnbr; j++){
        
        #if defined(USE_SFDL)
        nbr_pid.vec = idx[DM_LIN_ID(i, j, padded_Nnbr)];
        #else
        for (int k = 0; k < TILE_LENGTH; k++) nbr_pid.arr[k] = idx[DM_LIN_ID(TILE_LENGTH * i + k, j, padded_Nnbr)];
        #endif

        // Use ID from idx matrix to get state variables from K
        for (int k = 0; k < TILE_LENGTH; k++) u.arr[k] = K[SVM_LIN_ID(nbr_pid.arr[k], 0)];
        for (int k = 0; k < TILE_LENGTH; k++) v.arr[k] = K[SVM_LIN_ID(nbr_pid.arr[k], 1)];
        for (int k = 0; k < TILE_LENGTH; k++) w.arr[k] = K[SVM_LIN_ID(nbr_pid.arr[k], 2)];
        for (int k = 0; k < TILE_LENGTH; k++) h.arr[k] = K[SVM_LIN_ID(nbr_pid.arr[k], 3)];
        
        // Update sums- save vectors to private memory so that partial tiling (tiling either DM or SVM matrices) is possible
        #if defined(USE_SFDL)
        Dx_temp.vec = Dx[DM_LIN_ID(i, j, padded_Nnbr)];
        Dy_temp.vec = Dy[DM_LIN_ID(i, j, padded_Nnbr)];
        Dz_temp.vec = Dz[DM_LIN_ID(i, j, padded_Nnbr)];
        L_temp.vec = L[DM_LIN_ID(i, j, padded_Nnbr)];
        #else
        for (int k = 0; k < TILE_LENGTH; k++) Dx_temp.arr[k] = Dx[DM_LIN_ID(TILE_LENGTH * i + k, j, padded_Nnbr)];
        for (int k = 0; k < TILE_LENGTH; k++) Dy_temp.arr[k] = Dy[DM_LIN_ID(TILE_LENGTH * i + k, j, padded_Nnbr)];
        for (int k = 0; k < TILE_LENGTH; k++) Dz_temp.arr[k] = Dz[DM_LIN_ID(TILE_LENGTH * i + k, j, padded_Nnbr)];
        for (int k = 0; k < TILE_LENGTH; k++) L_temp.arr[k] = L[DM_LIN_ID(TILE_LENGTH * i + k, j, padded_Nnbr)];
        #endif
        
        du_dx += Dx_temp.vec * u.vec; // Summation notation: d_m/d_n = (D^n)_(ij) * (H_m)^(ij), where m = u,v,w,h;
        dv_dx += Dx_temp.vec * v.vec; // n = x,y,z; i = 0,1,2,...,nnbr, j = 0,1,2,...,nnodes
        dw_dx += Dx_temp.vec * w.vec;
        dh_dx += Dx_temp.vec * h.vec;
        
        du_dy += Dy_temp.vec * u.vec;
        dv_dy += Dy_temp.vec * v.vec;
        dw_dy += Dy_temp.vec * w.vec;
        dh_dy += Dy_temp.vec * h.vec;
        
        du_dz += Dz_temp.vec * u.vec;
        dv_dz += Dz_temp.vec * v.vec;
        dw_dz += Dz_temp.vec * w.vec;
        dh_dz += Dz_temp.vec * h.vec;
        
        Lu += L_temp.vec * u.vec;
        Lv += L_temp.vec * v.vec;
        Lw += L_temp.vec * w.vec;
        Lh += L_temp.vec * h.vec;
    }
    
    // Reset to get state variables at current node (rather than final neighbor).
    for (int k = 0; k < TILE_LENGTH; k++) u.arr[k] = K[SVM_LIN_ID(TILE_LENGTH * i + compute_pid_s + k, 0)];
    for (int k = 0; k < TILE_LENGTH; k++) v.arr[k] = K[SVM_LIN_ID(TILE_LENGTH * i + compute_pid_s + k, 1)];
    for (int k = 0; k < TILE_LENGTH; k++) w.arr[k] = K[SVM_LIN_ID(TILE_LENGTH * i + compute_pid_s + k, 2)];
    for (int k = 0; k < TILE_LENGTH; k++) h.arr[k] = K[SVM_LIN_ID(TILE_LENGTH * i + compute_pid_s + k, 3)];
    
    #if defined(USE_SFDL)
    x_temp.vec = x[i];
    y_temp.vec = y[i];
    z_temp.vec = z[i];
    f_temp.vec = f[i];
    ghm_temp.vec = ghm[i];
    #else
    for (int k = 0; k < TILE_LENGTH; k++) x_temp.arr[k] = x[TILE_LENGTH * i + k];
    for (int k = 0; k < TILE_LENGTH; k++) y_temp.arr[k] = y[TILE_LENGTH * i + k];
    for (int k = 0; k < TILE_LENGTH; k++) z_temp.arr[k] = z[TILE_LENGTH * i + k];
    for (int k = 0; k < TILE_LENGTH; k++) f_temp.arr[k] = f[TILE_LENGTH * i + k];
    for (int k = 0; k < TILE_LENGTH; k++) ghm_temp.arr[k] = ghm[TILE_LENGTH * i + k];
    #endif
    
    // Evaluate RHS of unconstrained momentum equations
    vec_f rhs_u = - ((u.vec * du_dx) + (v.vec * du_dy) + (w.vec * du_dz) + (f_temp.vec * ((y_temp.vec * w.vec) - (z_temp.vec * v.vec))) + dh_dx);
    vec_f rhs_v = - ((u.vec * dv_dx) + (v.vec * dv_dy) + (w.vec * dv_dz) + (f_temp.vec * ((z_temp.vec * u.vec) - (x_temp.vec * w.vec))) + dh_dy);
    vec_f rhs_w = - ((u.vec * dw_dx) + (v.vec * dw_dy) + (w.vec * dw_dz) + (f_temp.vec * ((x_temp.vec * v.vec) - (y_temp.vec * u.vec))) + dh_dz);
    
    // Get x, y, and z components of projection matrices
    for (int k = 0; k < 3; k++) {
        #if defined(USE_SFDL)
        p_u_temp.vec[k] = p_u[VECT_LIN_ID(i, k)];
        p_v_temp.vec[k] = p_v[VECT_LIN_ID(i, k)];
        p_w_temp.vec[k] = p_w[VECT_LIN_ID(i, k)];
        gradghm_temp.vec[k] = gradghm[VECT_LIN_ID(i, k)];
        #else
        for (int m = 0; m < TILE_LENGTH; m++) p_u_temp.arr[k][m] = p_u[VECT_LIN_ID(TILE_LENGTH * i + m, k)];
        for (int m = 0; m < TILE_LENGTH; m++) p_v_temp.arr[k][m] = p_v[VECT_LIN_ID(TILE_LENGTH * i + m, k)];
        for (int m = 0; m < TILE_LENGTH; m++) p_w_temp.arr[k][m] = p_w[VECT_LIN_ID(TILE_LENGTH * i + m, k)];
        for (int m = 0; m < TILE_LENGTH; m++) gradghm_temp.arr[k][m] = gradghm[VECT_LIN_ID(TILE_LENGTH * i + m, k)];
        #endif
    }
    
    // evaluate projections and apply hyperviscosity
    #if !defined(USE_CFDL)
    F[F_ID(i, 0)] = (p_u_temp.vec[0] * rhs_u) + (p_u_temp.vec[1] * rhs_v) + (p_u_temp.vec[2] * rhs_w) + Lu; // F_u = ...
    F[F_ID(i, 1)] = (p_v_temp.vec[0] * rhs_u) + (p_v_temp.vec[1] * rhs_v) + (p_v_temp.vec[2] * rhs_w) + Lv; // F_v = ...
    F[F_ID(i, 2)] = (p_w_temp.vec[0] * rhs_u) + (p_w_temp.vec[1] * rhs_v) + (p_w_temp.vec[2] * rhs_w) + Lw; // F_w = ...
    F[F_ID(i, 3)] = - ((u.vec * (dh_dx - gradghm_temp.vec[0])) + (v.vec * (dh_dy - gradghm_temp.vec[1])) + (w.vec * (dh_dz - gradghm_temp.vec[2])) + ((h.vec + gh0 - ghm_temp.vec) * (du_dx + dv_dy + dw_dz))) + Lh; // F_h = ...
    #else
    for (int k = 0; k < TILE_LENGTH; k++) F[F_ID(TILE_LENGTH * i + k, 0)] = (p_u_temp.arr[0][k] * ACCESS(rhs_u, k)) + (p_u_temp.arr[1][k] * ACCESS(rhs_v, k)) + (p_u_temp.arr[2][k] * ACCESS(rhs_w, k)) + ACCESS(Lu, k);
    for (int k = 0; k < TILE_LENGTH; k++) F[F_ID(TILE_LENGTH * i + k, 1)] = (p_v_temp.arr[0][k] * ACCESS(rhs_u, k)) + (p_v_temp.arr[1][k] * ACCESS(rhs_v, k)) + (p_v_temp.arr[2][k] * ACCESS(rhs_w, k)) + ACCESS(Lv, k);
    for (int k = 0; k < TILE_LENGTH; k++) F[F_ID(TILE_LENGTH * i + k, 2)] = (p_w_temp.arr[0][k] * ACCESS(rhs_u, k)) + (p_w_temp.arr[1][k] * ACCESS(rhs_v, k)) + (p_w_temp.arr[2][k] * ACCESS(rhs_w, k)) + ACCESS(Lw, k);
    for (int k = 0; k < TILE_LENGTH; k++) F[F_ID(TILE_LENGTH * i + k, 3)] = - ((u.arr[k] * (ACCESS(dh_dx, k) - gradghm_temp.arr[0][k])) + (v.arr[k] * (ACCESS(dh_dy, k) - gradghm_temp.arr[1][k])) + (w.arr[k] * (ACCESS(dh_dz, k) - gradghm_temp.arr[2][k])) + ((h.arr[k] + gh0 - ghm_temp.arr[k]) * (ACCESS(du_dx, k) + ACCESS(dv_dy, k) + ACCESS(dw_dz, k)))) + ACCESS(Lh, k);
    #endif
    return;
}

__kernel void copy_arr_kernel(__global double* D, __global const double* F, const int arr_size) {
    int i = get_global_id(0);
    if (i >= arr_size) return;
    
    D[i] = F[i]; // D = F
    
    return;
}

__kernel void update_D_kernel(__global const double* F, __global double* D, const double coefficient, const int Nnodes){
    int i = get_global_id(0);
    if (i >= Nnodes*4) return;
    
    D[i] += coefficient * F[i];  // D = c*F
    
    return;
}

__kernel void eval_K_kernel(__global double* H, __global const double* F, __global double* K,
                            const double dt, const int Nnodes, const int compute_pid_s){
    int i = get_global_id(0);
    if (i >= Nnodes*4) return;
    
    K[compute_pid_s*4 + i] = H[compute_pid_s*4 + i] + dt * F[i];        // K = H + (coeff)*dt*F (RK coefficient handled in invocation (RK_ocl.c)
    
    return;
}

__kernel void update_H_kernel(  __global double* H, const double dt, __global double* D,
                                __global double* K, const int Nnodes, const int compute_pid_s){
    int i = get_global_id(0);
    if (i >= Nnodes*4) return;
    
    H[compute_pid_s*4 + i] += (dt/6.0) * D[i]; // H += dt/6*D (Add all of the RK pieces and multiply by 1/6 per RK-45 algorithm)
    K[compute_pid_s*4 + i] = H[compute_pid_s*4 + i]; // Update to use H in the first round of the next time step w/o updating eval_RHS args
    return;
}
