#include <freeglut.h>        
#include <cuda_gl_interop.h>    

#include <cstdio>					// printf()
#include <math.h>					// powf() sqrt()

__global__
void g_findCellStartEndIndices(int number_of_particles, int* cell_start_indices, int* cell_end_indices, int* particle_hashes)
{
	int i = (blockIdx.x * blockDim.x + threadIdx.x);

	if (i < number_of_particles)
	{
		// Set the cell start hash table at the first particle's cell index (hash value) to be the first particle's index (0).
		cell_start_indices[particle_hashes[0]] = 0;
		if (i > 0)
		{
			// If this particle is in a difference cell to the previous particle.
			if (particle_hashes[i] != particle_hashes[i - 1])
			{
				// Set the cell end hash table at the previous particle's cell index (hash value) to be the current particle's index.
				cell_end_indices[particle_hashes[i - 1]] = i;

				// Set the cell start hash table at this particle's cell index (hash value) to be the current particle's index.
				cell_start_indices[particle_hashes[i]] = i;
			}
			// Set the cell end hash table at the final particle's cell index (hash value) to be the last particles index + 1 (number of particles).
			cell_end_indices[particle_hashes[number_of_particles - 1]] = number_of_particles;
		}
	}
}

__device__
void d_calculateGridHash(int grid_x, int grid_y, int& hash, int grid_width)
{
	// Calculate the linear cell index from the 2D cell array.
	hash = (grid_y * grid_width) + grid_x;
}

__device__
void d_calculateGridPosition(float position_x, float position_y, int& grid_x, int& grid_y, int cell_size, float scale)
{
	// Calculate x and y position of the cell the particle is in.
	grid_x = (int)((position_x + (512 * scale)) / cell_size);
	grid_y = (int)((position_y + (256 * scale)) / cell_size);
}

__global__
void g_calculateHashValues(int number_of_particles, int* particle_hashes, float* positions, int grid_width, int cell_size, float scale)
{
	int i = (blockIdx.x * blockDim.x + threadIdx.x);

	if (i < number_of_particles)
	{
		int grid_position_x;
		int grid_position_y;
		d_calculateGridPosition(positions[i * 2], positions[i * 2 + 1], grid_position_x, grid_position_y, cell_size, scale);
		d_calculateGridHash(grid_position_x, grid_position_y, particle_hashes[i], grid_width);
	}
}

__global__
void g_computeDensityPressure(float* positions, float* densities, float* pressures, int number_of_particles, float particle_diameter_squared, float particle_mass, float POLY6, float BOLTZMANN_CONSTANT, float REST_DENSITY, int cell_size, int grid_width, int grid_height, int* cell_start_indices, int* cell_end_indices, float scale)
{
	int i = (blockIdx.x * blockDim.x + threadIdx.x);

	if (i < number_of_particles)
	{
		// Storing particle variables for readability.
		float position_x = positions[i * 2];
		float position_y = positions[i * 2 + 1];

		// Pointers as density and pressure are being manipulated.
		float* density = &densities[i];
		float* pressure = &pressures[i];

		// Calculate the grid position of the particle.
		int grid_position_x;
		int grid_position_y;
		d_calculateGridPosition(position_x, position_y, grid_position_x, grid_position_y, cell_size, scale);

		// Reset density so it can be recalculated each iteration.
		*density = 0.0f;

		// Loop over surrounding 8 cells and itself.
		for (int y = -1; y <= 1; y++)
		{
			for (int x = -1; x <= 1; x++)
			{
				// Calculate the grid position of the neighboring cell.
				int neighboring_grid_position_x = grid_position_x + x;
				int neighboring_grid_position_y = grid_position_y + y;

				// If the cell is on an edge, ignore the cell positions outside the grid.
				if (neighboring_grid_position_x < 0 || neighboring_grid_position_x >= grid_width ||	
					neighboring_grid_position_y < 0 || neighboring_grid_position_y >= grid_height)
				{
					continue;
				}

				// Calculate hash value of neighboring cell.
				int neighboring_grid_hash;
				d_calculateGridHash(neighboring_grid_position_x, neighboring_grid_position_y, neighboring_grid_hash, grid_width);

				// Get the index of the first particle in this cell.
				int cell_start_index = cell_start_indices[neighboring_grid_hash];

				// Get the index of the first particle in the next cell.
				int cell_end_index = cell_end_indices[neighboring_grid_hash];

				// Loop over particles in cell.
				for (int j = cell_start_index; j < cell_end_index; j++)	
				{
					// Storing variables for other particle for readability.
					float other_position_x = positions[j * 2];
					float other_position_y = positions[j * 2 + 1];

					// Calculate vector between the pair of particles.
					float difference_x = other_position_x - position_x;
					float difference_y = other_position_y - position_y;

					// Calculate the squared distance between the pair of particles. 
					float distance_squared = powf(difference_x, 2) + powf(difference_y, 2);

					// If the particles are overlapping.
					if (distance_squared < particle_diameter_squared)
					{
						// Add the other particle's mass, scaled by a muller kernel to this particle's density.
						*density += particle_mass * POLY6 * powf(particle_diameter_squared - distance_squared, 3.0f);
					}
				}
			}
		}
		// Calculate pressure of particle using an "equation of state", relating it's density to a given rest density.
		*pressure = BOLTZMANN_CONSTANT * (*density - REST_DENSITY);
	}
}

__global__
void g_computeForces(float* positions, float* velocities, float* forces, float* densities, float* pressures, int number_of_particles, float particle_diameter, float particle_mass, float particle_viscocity, float gravity_x, float gravity_y, float SPIKY_GRAD, float VISC_LAP, int cell_size, int grid_width, int grid_height, int* cell_start_indices, int* cell_end_indices, int* particle_IDs, float scale)
{
	int i = (blockIdx.x * blockDim.x + threadIdx.x);

	if (i < number_of_particles)
	{
		// Storing particle variables for readability.
		float position_x = positions[i * 2];
		float position_y = positions[i * 2 + 1];
		float vx = velocities[i * 2];
		float vy = velocities[i * 2 + 1];
		float density = densities[i];
		float pressure = pressures[i];

		// Pointers as force is getting manipulated.
		float* fx = &forces[i * 2];
		float* fy = &forces[i * 2 + 1];

		// Calculate the grid position of the particle.
		int grid_position_x;
		int grid_position_y;
		d_calculateGridPosition(position_x, position_y, grid_position_x, grid_position_y, cell_size, scale);

		// Reset pressure and viscocity contributions so they can be recalculated each iteration.
		float pressure_contribution_x = 0.0f;
		float pressure_contribution_y = 0.0f;
		float viscocity_contribution_x = 0.0f;
		float viscocity_contribution_y = 0.0f;

		// Loop over surrounding 8 cells and itself.
		for (int y = -1; y <= 1; y++)
		{
			for (int x = -1; x <= 1; x++)
			{
				// Calculate the grid position of the neighboring cell.
				int neighboring_grid_position_x = grid_position_x + x;
				int neighboring_grid_position_y = grid_position_y + y;

				// If the cell is on an edge, ignore the cell positions outside the grid.
				if (neighboring_grid_position_x < 0 || neighboring_grid_position_x >= grid_width ||	// added to clamp to edges instead of weird wrapping nvidia seems to be doing
					neighboring_grid_position_y < 0 || neighboring_grid_position_y >= grid_height)
				{
					continue;
				}

				// Calculate hash value of neighboring cell.
				int neighboring_grid_hash;
				d_calculateGridHash(neighboring_grid_position_x, neighboring_grid_position_y, neighboring_grid_hash, grid_width);

				// Get the index of the first particle in this cell.
				int cell_start_index = cell_start_indices[neighboring_grid_hash];

				// Get the index of the first particle in the next cell.
				int cell_end_index = cell_end_indices[neighboring_grid_hash];

				// loop over particles in cell
				for (int j = cell_start_index; j < cell_end_index; j++)
				{
					// If it is not this particle.
					if (particle_IDs[j] != particle_IDs[i])
					{
						// Storing variables for other particle for readability.
						float other_position_x = positions[j * 2];
						float other_position_y = positions[j * 2 + 1];
						float other_vx = velocities[j * 2];
						float other_vy = velocities[j * 2 + 1];
						float other_density = densities[j];
						float other_pressure = pressures[j];

						// Calculate vector between the pair of particles.
						float difference_x = other_position_x - position_x;
						float difference_y = other_position_y - position_y;

						// Calculate the distance between the pair of particles. Distance needed later so no benefit from comparing squared distances.
						float distance = sqrt(powf(difference_x, 2.f) + powf(difference_y, 2.f));

						// If particles are overlapping.
						if (distance < particle_diameter)		// <--- distance != 0 is handling case where particles are at the same position.
						{
							// Calculate the desired direction vector of this particle.
							float desired_x = -difference_x / distance;
							float desired_y = -difference_y / distance;

							// Add other particle's pressure and viscocity contributions using Navier-Stokes equations, scaled by Muller kernels.
							pressure_contribution_x += desired_x * particle_mass * (pressure + other_pressure) / (2.f * other_density) * SPIKY_GRAD * powf(particle_diameter - distance, 2.f);
							pressure_contribution_y += desired_y * particle_mass * (pressure + other_pressure) / (2.f * other_density) * SPIKY_GRAD * powf(particle_diameter - distance, 2.f);

							viscocity_contribution_x += particle_viscocity * particle_mass * (other_vx - vx) / other_density * VISC_LAP * (particle_diameter - distance);
							viscocity_contribution_y += particle_viscocity * particle_mass * (other_vy - vy) / other_density * VISC_LAP * (particle_diameter - distance);

						}
					}
				}
				
			}
		}
		// Calculate grivity contributions by multiplying gravity by this particle's calculated density.
		float gravity_contribution_x = gravity_x * density;
		float gravity_contribution_y = gravity_y * density;

		// Add all force contributions together to calculate the particle's force.
		*fx = pressure_contribution_x + viscocity_contribution_x + gravity_contribution_x;
		*fy = pressure_contribution_y + viscocity_contribution_y + gravity_contribution_y;
	}
}

__global__
void g_integrate(float scale, float* positions, float* velocities, float* forces, float* densities, int number_of_particles, float integration_timestep, float damping, float jitter)
{
	int i = (blockIdx.x * blockDim.x + threadIdx.x);

	if (i < number_of_particles)
	{
		// Storing particle variables for readability.
		float* x = &positions[i * 2];
		float* y = &positions[i * 2 + 1];
		float* vx = &velocities[i * 2];
		float* vy = &velocities[i * 2 + 1];
		float* fx = &forces[i * 2];
		float* fy = &forces[i * 2 + 1];
		float density = densities[i];

		// Update velocity of particle.
		*vx += *fx / density * integration_timestep;
		*vy += *fy / density * integration_timestep;

		// Update position of particle.
		*x += *vx * integration_timestep;
		*y += *vy * integration_timestep;

		//// Handle collision detection at edges of environment rectangle.
		// LEFT
		if (*x < -512.0f * scale)
		{
			*vx *= damping;
			*x = (-512.0f * scale) + jitter;
		}
		// RIGHT
		if (*x > 512.0f * scale)
		{
			*vx *= damping;
			*x = (512.0f * scale) - jitter;
		}
		// TOP
		if (*y > 256.0f * scale)
		{
			*vy *= damping;
			*y = (256.0f * scale) - jitter;
		}
		// BOTTOM
		if (*y < -256.0f * scale)
		{
			*vy *= damping;
			*y = (-256.0f * scale) + jitter;
		}
	}
}

void errorCheck(const int line, const char* const file)
{
    // Checks last cuda kernel for errors.
    cudaError_t error = cudaGetLastError();
    // If there was an error.
    if (error != cudaSuccess)
    {
        // Print error message with line number and file.
        printf("CUDA ERROR: %s at line %d in\n %s\n", cudaGetErrorString(error), line, file);
        // Exit the application.
        exit(EXIT_FAILURE);
    }
}

void c_cudaMalloc(void** devPtr, int size) 
{
    // Allocates device memory.
    cudaMalloc(devPtr, size);
    errorCheck(__LINE__, __FILE__);
}

void c_cudaFree(void* devPtr)
{
    // Releases device memory.
    cudaFree(devPtr);
    errorCheck(__LINE__, __FILE__);
}

void c_registerVBO(struct cudaGraphicsResource** cuda_vbo_resource, unsigned int vbo)
{
    // Registers OpenGL VBO and returns a handle to it.
    cudaGraphicsGLRegisterBuffer(cuda_vbo_resource, vbo, cudaGraphicsMapFlagsNone);
    errorCheck(__LINE__, __FILE__);
}

void c_unregisterVBO(struct cudaGraphicsResource* cuda_vbo_resource)
{
    // Unregisters VBO.
    cudaGraphicsUnregisterResource(cuda_vbo_resource);
    errorCheck(__LINE__, __FILE__);
}

void* c_mapVBO(struct cudaGraphicsResource** cuda_vbo_resource)
{
    // Maps VBO for access by CUDA.
    cudaGraphicsMapResources(1, cuda_vbo_resource);
    errorCheck(__LINE__, __FILE__);

    void* vbo_pointer;
    size_t bytes;
    // Returns a pointer to the mapped VBO.
    cudaGraphicsResourceGetMappedPointer((void**)&vbo_pointer, &bytes, *cuda_vbo_resource);
    errorCheck(__LINE__, __FILE__);

    return vbo_pointer;
}

void c_unmapVBO(struct cudaGraphicsResource* cuda_vbo_resource)
{
    // Unmaps VBO.
    cudaGraphicsUnmapResources(1, &cuda_vbo_resource, 0);
    errorCheck(__LINE__, __FILE__);
}

void c_cudaMemcpyHostToDevice(void* device, const void* host, int size)
{
	// Copies array from host to device.
	cudaMemcpy((char*)device, host, size, cudaMemcpyHostToDevice);
	errorCheck(__LINE__, __FILE__);
}

void c_cudaMemcpyDeviceToHost(void* host, const void* device, struct cudaGraphicsResource** cuda_vbo_resource, int size)
{
	// Maps Graphics Resource for access by CUDA and stores the pointer to it.
	device = c_mapVBO(cuda_vbo_resource);

	// Copy data from device to host.
	cudaMemcpy(host, device, size, cudaMemcpyDeviceToHost);
	errorCheck(__LINE__, __FILE__);

	// Unmaps Graphics Resource.
	c_unmapVBO(*cuda_vbo_resource);
}

void c_calculateHashValues(int threads, int blocks, int number_of_particles, int* particle_hashes, float* positions, int grid_width, int cell_size, float scale)
{
    g_calculateHashValues<<<threads, blocks>>>(number_of_particles, particle_hashes, positions, grid_width, cell_size, scale);
	errorCheck(__LINE__, __FILE__);

	// Sync threads before continuing.
    cudaDeviceSynchronize();
    errorCheck(__LINE__, __FILE__);
}

void c_findCellStartEndIndices(int threads, int blocks, int number_of_particles, int* cell_start_indices, int* cell_end_indices, int* particle_hashes)
{
    g_findCellStartEndIndices<<<threads, blocks>>>(number_of_particles, cell_start_indices, cell_end_indices, particle_hashes);
	errorCheck(__LINE__, __FILE__);

	// Sync threads before continuing.
    cudaDeviceSynchronize();
    errorCheck(__LINE__, __FILE__);
}

void c_computeDensityPressure(int threads, int blocks, float* positions, float* densities, float* pressures, int number_of_particles, float particle_diameter_squared, float particle_mass, float POLY6, float BOLTZMANN_CONSTANT, float REST_DENSITY, int cell_size, int grid_width, int grid_height, int* cell_start_indices, int* cell_end_indices, float scale)
{
    g_computeDensityPressure<<<blocks, threads>>>(positions, densities, pressures, number_of_particles, particle_diameter_squared, particle_mass, POLY6, BOLTZMANN_CONSTANT, REST_DENSITY, cell_size, grid_width, grid_height, cell_start_indices, cell_end_indices, scale);
    errorCheck(__LINE__, __FILE__);

	// Sync threads before continuing.
    cudaDeviceSynchronize();
    errorCheck(__LINE__, __FILE__);
}

void c_computeForces(int threads, int blocks, float* positions, float* velocities, float* forces, float* densities, float* pressures, int number_of_particles, float particle_diameter, float particle_mass, float particle_viscocity, float gravity_x, float gravity_y, float SPIKY_GRAD, float VISC_LAP, int cell_size, int grid_width, int grid_height, int* cell_start_indices, int* cell_end_indices, int* particle_IDs, float scale)
{
    g_computeForces<<<blocks, threads>>>(positions, velocities, forces, densities, pressures, number_of_particles, particle_diameter, particle_mass, particle_viscocity, gravity_x, gravity_y, SPIKY_GRAD, VISC_LAP, cell_size, grid_width, grid_height, cell_start_indices, cell_end_indices, particle_IDs, scale);
    errorCheck(__LINE__, __FILE__);

	// Sync threads before continuing.
    cudaDeviceSynchronize();
    errorCheck(__LINE__, __FILE__);
}

void c_integrate(int threads, int blocks, float scale, float* positions, float* velocities, float* forces, float* densities, int number_of_particles, float integration_timestep, float damping)
{
    // A small offest applied to particles on boundary collision to avoid stacking.
    float jitter = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    g_integrate<<<blocks, threads>>>(scale, positions, velocities, forces, densities, number_of_particles, integration_timestep, damping, jitter);
    errorCheck(__LINE__, __FILE__);

}


