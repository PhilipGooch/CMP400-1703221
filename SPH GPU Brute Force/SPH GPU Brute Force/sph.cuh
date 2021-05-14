
// c_  signifies it is a function that interacts with CUDA, that can be called from a .cpp file.

void c_cudaMalloc(void** devPtr, int size);		// Allocates device memory.
void c_cudaFree(void* devPtr);					// Releases device memory.

void c_registerVBO(struct cudaGraphicsResource** cuda_vbo_resource, unsigned int vbo);	// Registers VBO and returns a handle to it.
void c_unregisterVBO(struct cudaGraphicsResource* cuda_vbo_resource);					// Unregisters VBO.

void* c_mapVBO(struct cudaGraphicsResource** cuda_vbo_resource);	// Maps VBO for access by CUDA and returns a pointer to it.
void c_unmapVBO(struct cudaGraphicsResource* cuda_vbo_resource);	// Unmaps the VBO.

// Simulation functions that call corresponding device kernels.
void c_computeDensityPressure(int threads, int blocks, float* positions, float* densitys, float* pressures, int number_of_particles, float particle_diameter_squared, float particle_mass, float POLY6, float BOLTZMANN_CONSTANT, float REST_DENSITY);
void c_computeForces(int threads, int blocks, float* positions, float* velocities, float* forces, float* densities, float* pressures, int number_of_particles, float particle_diameter, float particle_mass, float particle_viscocity, float gravity_x, float gravity_y, float SPIKY_GRAD, float VISC_LAP);
void c_integrate(int threads, int blocks, float scale, float* positions, float* velocities, float* forces, float* densities, int number_of_particles, float integration_timestep, float damping);

