
// c_  signifies it is a function that interacts with CUDA, that can be called from a .cpp file.

void c_cudaMalloc(void** devPtr, int size);		// Allocates device memory.
void c_cudaFree(void* devPtr);					// Releases device memory.

void c_registerVBO(struct cudaGraphicsResource** cuda_vbo_resource, unsigned int vbo);	// Registers VBO and returns a handle to it.
void c_unregisterVBO(struct cudaGraphicsResource* cuda_vbo_resource);					// Unregisters VBO.

void* c_mapVBO(struct cudaGraphicsResource** cuda_vbo_resource);	// Maps VBO for access by CUDA and returns a pointer to it.
void c_unmapVBO(struct cudaGraphicsResource* cuda_vbo_resource);	// Unmaps the VBO.

// Copy arrays between host and device.
void c_cudaMemcpyDeviceToHost(void* host, const void* device, struct cudaGraphicsResource** cuda_vbo_resource, int size);
void c_cudaMemcpyHostToDevice(void* device, const void* host, int size);

// Uniform Grid functions that call corresponding device kernels.
void c_calculateHashValues(int threads, int blocks, int number_of_particles, int* particle_hashes, float* positions, int grid_width, int cell_size, float scale);
//void c_sortParticles(int* particle_hashes, int* particle_IDs, int number_of_particles, float* positions, float* velocities, float* forces, float* densities, float* pressures);
void c_findCellStartEndIndices(int threads, int blocks, int number_of_particles, int* cell_start_indices, int* cell_end_indices, int* particle_hashes);

// Simulation functions that call corresponding device kernels.
void c_computeDensityPressure(int threads, int blocks, float* positions, float* densitys, float* pressures, int number_of_particles, float particle_diameter_squared, float particle_mass, float POLY6, float BOLTZMANN_CONSTANT, float REST_DENSITY, int cell_size, int grid_width, int grid_height, int* cell_start_indices, int* cell_end_indices, float scale);
void c_computeForces(int threads, int blocks, float* positions, float* velocities, float* forces, float* densities, float* pressures, int number_of_particles, float particle_diameter, float particle_mass, float particle_viscocity, float gravity_x, float gravity_y, float SPIKY_GRAD, float VISC_LAP, int cell_size, int grid_width, int grid_height, int* cell_start_indices, int* cell_end_indices, int* particle_IDs, float scale);
void c_integrate(int threads, int blocks, float scale, float* positions, float* velocities, float* forces, float* densities, int number_of_particles, float integration_timestep, float damping);

