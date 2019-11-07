// Written by Tuan Ta, 2016- from Gen1 RBF_SWE code. Last updated- Richelle Streater, Summer 2018.

#ifdef USE_OCL
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#include <device_setup.h>

extern sim_params_struct sim_params;

extern int mpi_rank;
extern int mpi_size;
 
// Choose from a list of possible platforms (platformID is the index in the list to choose)
void createPlatform(cl_platform_id* platform, int platformID){
    cl_uint platformCount = 0;

    // Figure out how many platforms are available
    cl_int err = clGetPlatformIDs(1, NULL, &platformCount);

    // Create array of all possible platforms
    cl_platform_id platforms[platformCount];

    // Add found platform ID's to the platforms array
    err = clGetPlatformIDs(1, platforms, NULL);

    if (err != CL_SUCCESS){
        printf("Error in clGetPlatformIDs: %s\n", checkError(err));
        exit(-1);
    }
    
    if (mpi_rank == 0) {

        // Print device info
        printf("\n\n============================================== OpenCL Device Information ===================================================\n\n");
        printf("Available platforms:\n\n");
        for (int i = 0; i < platformCount; i++){
            printf("Platform # %d\t-> \t", i);
            char buffer[10240];
            err = clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, sizeof(buffer), buffer, NULL);
            printf("%s\n", buffer);
	    }

        printf(">>> Choose platform %d\n", platformID);
    }
	
    // Select from the list of possible platforms
    *platform = platforms[platformID];
}

// Create a context, picking a device to use out of a list (deviceID is the index of the device to choose) 
void createContext(cl_context* context, cl_platform_id platform, cl_device_id* device, int deviceID){
    cl_int err = 0;
    cl_uint deviceCount = 0;

    // Find available devices
    err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, 0, NULL, &deviceCount);

    // Create array of possible devices
    cl_device_id devices[deviceCount];

    // Fill array of devices with found devices
    err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, deviceCount, devices, NULL);

    if (err != CL_SUCCESS){
        printf("Error in clGetDeviceIDs: %s\n", checkError(err));
        exit(-1);
    }
	
if (mpi_rank == 0) {

    // Print device info
    printf("\nAvailable devices:\n\n");
	
    for (int i = 0; i < deviceCount; i++){
        printf("Device # %d\t-> \t", i);
		
        // Buffer stores messages about the device
        char buffer[10240];
        err = clGetDeviceInfo(devices[i], CL_DEVICE_NAME, sizeof(buffer), buffer, NULL);
        printf("%s\n", buffer);
    }
	
    if (deviceID >= deviceCount){
        printf("Invalid device id\n");
        exit(-1);
    }
	
    printf(">>> Choose device %d\n", deviceID);
        printf("\nNumber of subdevices: %d\n", sim_params.num_compute_units);
    }
    
    cl_device_id subDevice;
    cl_device_partition_property* props = (cl_device_partition_property*) malloc(sizeof(cl_device_partition_property) *    (sim_params.num_compute_units + 3));
    
    #ifdef SPLITDEV
    props[0] = CL_DEVICE_PARTITION_BY_NAMES_EXT;
    
    // Partition by names- if device splitting fails, run in serial
    for (int i = mpi_rank; i < sim_params.num_compute_units/mpi_size + mpi_rank; i++) props[1 + i] = i;
    props[sim_params.num_compute_units + 1] = CL_PARTITION_BY_NAMES_LIST_END_EXT;
    props[sim_params.num_compute_units + 2] = 0;

    err = clCreateSubDevices(devices[deviceID], props, 1, &subDevice, NULL);
    if (err != CL_SUCCESS) {
        if (mpi_rank == 0) printf("Error in clCreateSubDevices: %s\n--> Executing in serial\n", checkError(err));
        sim_params.num_compute_units = 1;
        *device = devices[deviceID];
    }
    else *device = subDevice;

    #else
    if (mpi_rank == 0) printf("Executable is not set up to support device splitting (SPLIT_DEV=0)\n--> Executing in serial\n");
    sim_params.num_compute_units = 1;
    *device = devices[deviceID];
    #endif
	
    *context = clCreateContext(0, 1, device, NULL, NULL, &err);
    if (err != CL_SUCCESS){
        printf("Error in clCreateContext: %s\n", checkError(err));
        exit(-1);
    }
    if (mpi_rank == 0) printf("\n============================================================================================================================\n\n");
}

// Create command queue: Will give instructions about what to do and when
void createCommandQueue(cl_context context, cl_device_id device, cl_command_queue* commandQueue){
	cl_int err = 0;
	
	// Use the context and device (already defined) to form a command queue
	*commandQueue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &err);
	
	if (err != CL_SUCCESS){
		printf("Error in clCreateCommandQueue: %s\n", checkError(err));
		exit(-1);
	}	
}

// Compile OpenCL kernel from a file/kernel name
void compileKernel(const char* kernelFileName, const char* kernelName, cl_context context, cl_program* program, cl_device_id* device, cl_kernel* kernel){
    cl_int err = 0;

    // Read the OpenCL kernel from source file
    char* source = readKernelSource(kernelFileName);

    // Create a program with kernel source file
    *program = clCreateProgramWithSource(context, 1, (const char **)&source, NULL, &err);

    if (err != CL_SUCCESS){
        printf("Error in clCreateProgramWithSource: %s\n", checkError(err));
        exit(-1);
    }

    // Build the program
    err = clBuildProgram(*program, 1, device, NULL, NULL, NULL);

    if (err != CL_SUCCESS){
        printf("Error in clBuildProgram: %s\n", checkError(err));
        size_t build_log_size;
        err = clGetProgramBuildInfo(*program, *device, CL_PROGRAM_BUILD_LOG, 0, NULL, &build_log_size);

        char log[build_log_size];
        err = clGetProgramBuildInfo(*program, *device, CL_PROGRAM_BUILD_LOG, build_log_size, log, NULL);

        printf("%s\n", log);
        exit(-1);
    }

    // Create kernel
    *kernel = clCreateKernel(*program, kernelName, &err);

    if (err != CL_SUCCESS){
        printf("Error in clCreateKernel: %s\n", checkError(err));
        exit(-1);
    }
}

// Create a buffer to read into or write from a kernel- usually an array
cl_mem createBuffer(cl_context context, cl_mem_flags flags, size_t buffSize){
    cl_int err = 0;
    cl_mem buffer = clCreateBuffer(context, flags, buffSize, NULL, &err);

    if (err != CL_SUCCESS){
        printf("Error in clCreateBuffer: %s\n", checkError(err));
        exit(-1);
    }

    return buffer;
}

// Add a map of a buffer to the command queue- allows reading/writing of the buffer (select with flags)
void* mapBuffer(cl_command_queue commandQueue, cl_mem buffer, cl_map_flags flags, size_t buffSize){
	cl_int err = 0;
	
    // Add map to the command queue; make it map the ENTIRE buffer. Include blocking
    void* ptr = clEnqueueMapBuffer(commandQueue, buffer, CL_TRUE, flags, 0, buffSize, 0u, NULL, NULL, &err);
	
    if (err != CL_SUCCESS){
        printf("Error in clEnqueueMapBuffer: %s\n", checkError(err));
        exit(-1);
    }
	
    return ptr;
}

// Cleans up map of buffer when finished
void unmapBuffer(cl_command_queue commandQueue, cl_mem buffer, void* ptr){
    cl_int err = clEnqueueUnmapMemObject(commandQueue, buffer, ptr, 0u, NULL, NULL);

    if (err != CL_SUCCESS){
        printf("Error in clEnqueueUnmapMemObject: %s\n", checkError(err));
        exit(-1);
    }
}

// Waits for command queue to finish tasks
void commandQueueSync(cl_command_queue commandQueue){
    cl_int err = clFinish(commandQueue);

    if (err != CL_SUCCESS){
        printf("Error in clFinish: %s\n", checkError(err));
        exit(-1);
    }
}

// Set all inputs for the kernel that evaluates the SWE righthand side
void set_eval_RHS_args(cl_kernel kernel, LPSMD_buffers* LPSMD_buffs, cl_mem K_buff, cl_mem F_buff, fType gh0, int Nnodes, int padded_Nnbr, int Nnbr, int compute_pid_s){
    cl_int err = 0;
	
    err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void*)&(K_buff));
    err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void*)&(LPSMD_buffs->idx));
    err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void*)&(LPSMD_buffs->Dx));
    err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void*)&(LPSMD_buffs->Dy));
    err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void*)&(LPSMD_buffs->Dz));
    err = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void*)&(LPSMD_buffs->L));
    err = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void*)&(LPSMD_buffs->x));
    err = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void*)&(LPSMD_buffs->y));
    err = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void*)&(LPSMD_buffs->z));
    err = clSetKernelArg(kernel, 9, sizeof(cl_mem), (void*)&(LPSMD_buffs->f));
    err = clSetKernelArg(kernel, 10, sizeof(cl_mem), (void*)&(LPSMD_buffs->ghm));
    err = clSetKernelArg(kernel, 11, sizeof(cl_mem), (void*)&(LPSMD_buffs->p_u));
    err = clSetKernelArg(kernel, 12, sizeof(cl_mem), (void*)&(LPSMD_buffs->p_v));
    err = clSetKernelArg(kernel, 13, sizeof(cl_mem), (void*)&(LPSMD_buffs->p_w));
    err = clSetKernelArg(kernel, 14, sizeof(cl_mem), (void*)&(LPSMD_buffs->gradghm));
    err = clSetKernelArg(kernel, 15, sizeof(cl_mem), (void*)&(F_buff));
    err = clSetKernelArg(kernel, 16, sizeof(fType), (void*)&gh0);
    err = clSetKernelArg(kernel, 17, sizeof(int), (void*)&Nnodes);
    err = clSetKernelArg(kernel, 18, sizeof(int), (void*)&padded_Nnbr);
    err = clSetKernelArg(kernel, 19, sizeof(int), (void*)&Nnbr);
    err = clSetKernelArg(kernel, 20, sizeof(int), (void*)&compute_pid_s);
}
    
// Set all inputs for the kernel that copies F's contents to D's
void set_copy_arr_args(cl_kernel kernel, cl_mem D_buff, cl_mem F_buff, const int arr_size){
    cl_int err = 0;
                             
    err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void*)&D_buff);
    err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void*)&F_buff);
    err = clSetKernelArg(kernel, 2, sizeof(int), (void*)&arr_size);
}
                         
// Set all inputs for the kernel that updates D (represents differentiation matrix x velocity)
void set_update_D_args(cl_kernel kernel, cl_mem F_buff, cl_mem D_buff, const int Nnodes){
    cl_int err = 0;
                             
    err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void*)&F_buff);
    err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void*)&D_buff);
    err = clSetKernelArg(kernel, 3, sizeof(int), (void*)&Nnodes);
}
                         
                         
// Set all inputs for the kernel that calculates K in the SWE's
void set_eval_K_args(cl_kernel kernel, cl_mem H_buff, cl_mem F_buff, const int Nnodes, cl_mem K_buff, int compute_pid_s){
    cl_int err = 0;

    err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void*)&H_buff);
    err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void*)&F_buff);
    err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void*)&K_buff);
    err = clSetKernelArg(kernel, 4, sizeof(int), (void*)&Nnodes);
    err = clSetKernelArg(kernel, 5, sizeof(int), (void*)&compute_pid_s);
}

// Set arguments for the kernel that updates H in the SWE's
void set_update_H_args (cl_kernel kernel, cl_mem H_buff, cl_mem D_buff, cl_mem K_buff, const int Nnodes, const fType dt, int compute_pid_s){
    cl_int err = 0;

    err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void*)&H_buff);
    err = clSetKernelArg(kernel, 1, sizeof(fType), (void*)&dt);
    err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void*)&D_buff);
    err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void*)&K_buff);
    err = clSetKernelArg(kernel, 4, sizeof(cl_int), (void*)&Nnodes);
    err = clSetKernelArg(kernel, 5, sizeof(cl_int), (void*)&compute_pid_s);
}

// Read the kernel text into our program from the other file
char* readKernelSource(const char *fileName){

    // Load the kernel source code into *char source
    FILE *fp;
    char *source;
    char *SIMD_len_str;
    char *source_full;
    size_t size;

    fp = fopen(fileName, "r");
    if (!fp) {
        fprintf(stderr, "Failed to load kernel.\n");
        exit(1);
    }

    // Measure the size of source file to tell where to stop
    fseek(fp, 0, SEEK_END);
    size = ftell(fp);

    fseek(fp, 0, SEEK_SET);

    source = (char*)malloc(size+1);
    size = fread(source, 1, size, fp);
    source[size]='\0';
    fclose( fp );
    
    SIMD_len_str = (char*)malloc(4);
    source_full = (char*)malloc(80+strlen(source));

    // Read in SIMD length as a string
    sprintf(SIMD_len_str, "%d\n", SIMD_LENGTH);
    
    // Define the SIMD length in the kernel
    strcpy(source_full, "#define SIMD_LENGTH ");
    strcat(source_full, SIMD_len_str);
    
    // Tell the kernel which ordering the main program uses
    #ifdef USE_SFDL
    strcat(source_full, "#define USE_SFDL\n");
    #endif
    
    // CFDL- the default
    #ifdef USE_CFDL
    strcat(source_full, "#define USE_CFDL\n");
    #endif
    
    strcat(source_full, source);

    return source_full;
}

// Functions to print specific errors to the user
const char* checkError(cl_int err){
    switch (err) {
        case CL_SUCCESS:                            return "Success!";
        case CL_DEVICE_NOT_FOUND:                   return "Device not found.";
        case CL_DEVICE_NOT_AVAILABLE:               return "Device not available";
        case CL_COMPILER_NOT_AVAILABLE:             return "Compiler not available";
        case CL_MEM_OBJECT_ALLOCATION_FAILURE:      return "Memory object allocation failure";
        case CL_OUT_OF_RESOURCES:                   return "Out of resources";
        case CL_OUT_OF_HOST_MEMORY:                 return "Out of host memory";
        case CL_PROFILING_INFO_NOT_AVAILABLE:       return "Profiling information not available";
        case CL_MEM_COPY_OVERLAP:                   return "Memory copy overlap";
        case CL_IMAGE_FORMAT_MISMATCH:              return "Image format mismatch";
        case CL_IMAGE_FORMAT_NOT_SUPPORTED:         return "Image format not supported";
        case CL_BUILD_PROGRAM_FAILURE:              return "Program build failure";
        case CL_MAP_FAILURE:                        return "Map failure";
        case CL_INVALID_VALUE:                      return "Invalid value";
        case CL_INVALID_DEVICE_TYPE:                return "Invalid device type";
        case CL_INVALID_PLATFORM:                   return "Invalid platform";
        case CL_INVALID_DEVICE:                     return "Invalid device";
        case CL_INVALID_CONTEXT:                    return "Invalid context";
        case CL_INVALID_QUEUE_PROPERTIES:           return "Invalid queue properties";
        case CL_INVALID_COMMAND_QUEUE:              return "Invalid command queue";
        case CL_INVALID_HOST_PTR:                   return "Invalid host pointer";
        case CL_INVALID_MEM_OBJECT:                 return "Invalid memory object";
        case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:    return "Invalid image format descriptor";
        case CL_INVALID_IMAGE_SIZE:                 return "Invalid image size";
        case CL_INVALID_SAMPLER:                    return "Invalid sampler";
        case CL_INVALID_BINARY:                     return "Invalid binary";
        case CL_INVALID_BUILD_OPTIONS:              return "Invalid build options";
        case CL_INVALID_PROGRAM:                    return "Invalid program";
        case CL_INVALID_PROGRAM_EXECUTABLE:         return "Invalid program executable";
        case CL_INVALID_KERNEL_NAME:                return "Invalid kernel name";
        case CL_INVALID_KERNEL_DEFINITION:          return "Invalid kernel definition";
        case CL_INVALID_KERNEL:                     return "Invalid kernel";
        case CL_INVALID_ARG_INDEX:                  return "Invalid argument index";
        case CL_INVALID_ARG_VALUE:                  return "Invalid argument value";
        case CL_INVALID_ARG_SIZE:                   return "Invalid argument size";
        case CL_INVALID_KERNEL_ARGS:                return "Invalid kernel arguments";
        case CL_INVALID_WORK_DIMENSION:             return "Invalid work dimension";
        case CL_INVALID_WORK_GROUP_SIZE:            return "Invalid work group size";
        case CL_INVALID_WORK_ITEM_SIZE:             return "Invalid work item size";
        case CL_INVALID_GLOBAL_OFFSET:              return "Invalid global offset";
        case CL_INVALID_EVENT_WAIT_LIST:            return "Invalid event wait list";
        case CL_INVALID_EVENT:                      return "Invalid event";
        case CL_INVALID_OPERATION:                  return "Invalid operation";
        case CL_INVALID_GL_OBJECT:                  return "Invalid OpenGL object";
        case CL_INVALID_BUFFER_SIZE:                return "Invalid buffer size";
        case CL_INVALID_MIP_LEVEL:                  return "Invalid mip-map level";
        case CL_DEVICE_PARTITION_FAILED:            return "Device partition failed";
        default: return "Unknown";
	}
}
#endif
