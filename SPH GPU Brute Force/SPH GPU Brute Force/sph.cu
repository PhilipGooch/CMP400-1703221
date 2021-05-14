#include <freeglut.h>        
#include <cuda_gl_interop.h>    

#include <cstdio> // printf()

#include <math.h>

__global__
void g_computeDensityPressure(float* positions, float* densities, float* pressures, int number_of_particles, float particle_diameter_squared, float particle_mass, float POLY6, float BOLTZMANN_CONSTANT, float REST_DENSITY)
{
    int i = (blockIdx.x * blockDim.x + threadIdx.x);

    if (i < number_of_particles)
    {
        // Storing particle variables for readability.
        float x = positions[i * 2];
        float y = positions[i * 2 + 1];

        // Pointers as density and pressure are being manipulated.
        float* density = &densities[i];
        float* pressure = &pressures[i];

        // Reset density so it can be recalculated each iteration.
        *density = 0;

        for (int j = 0; j < number_of_particles; j++)
        {
            // Storing variables for other particle for readability.
            float other_x = positions[j * 2];
            float other_y = positions[j * 2 + 1];

            // Calculate vector between the pair of particles.
            float difference_x = other_x - x;
            float difference_y = other_y - y;

            // Calculate the squared distance between the pair of particles. 
            float distance_squared = powf(difference_x, 2) + powf(difference_y, 2);

            // If the particles are overlapping.
            if (distance_squared < particle_diameter_squared)
            {
                // Add the other particle's mass, scaled by a muller kernel to this particle's density.
                *density += particle_mass * POLY6 * powf(particle_diameter_squared - distance_squared, 3.0f);
            }
        }
        // Calculate pressure of particle using an "equation of state", relating it's density to a given rest density.
        *pressure = BOLTZMANN_CONSTANT * (*density - REST_DENSITY);
    }
}

__global__
void g_computeForces(float* positions, float* velocities, float* forces, float* densities, float* pressures, int number_of_particles, float particle_diameter, float particle_mass, float particle_viscocity, float gravity_x, float gravity_y, float SPIKY_GRAD, float VISC_LAP)
{
    int i = (blockIdx.x * blockDim.x + threadIdx.x);

    if (i < number_of_particles)
    {
        // Storing particle variables for readability.
        float x = positions[i * 2];
        float y = positions[i * 2 + 1];
        float vx = velocities[i * 2];
        float vy = velocities[i * 2 + 1];
        float density = densities[i];
        float pressure = pressures[i];

        // Pointers as force is getting manipulated.
        float* fx = &forces[i * 2];
        float* fy = &forces[i * 2 + 1];

        // Reset pressure and viscocity contributions so they can be recalculated each iteration.
        float pressure_contribution_x = 0.0f;
        float pressure_contribution_y = 0.0f;
        float viscocity_contribution_x = 0.0f;
        float viscocity_contribution_y = 0.0f;

        for (int j = 0; j < number_of_particles; j++)
        {
            // Do not compare particles if they are the same.
            if (i != j)
            {
                // Storing variables for other particle for readability.
                float other_x = positions[j * 2];
                float other_y = positions[j * 2 + 1];
                float other_vx = velocities[j * 2];
                float other_vy = velocities[j * 2 + 1];
                float other_density = densities[j];
                float other_pressure = pressures[j];

                // Calculate vector between the pair of particles.
                float difference_x = other_x - x;
                float difference_y = other_y - y;

                // Calculate the distance between the pair of particles. Distance needed later so no benefit from comparing squared distances.
                float distance = sqrt(powf(difference_x, 2.f) + powf(difference_y, 2.f));

                // If particles are overlapping.
                if (distance != 0 && distance < particle_diameter)	// <--- distance != 0 is handling case where particles are at the same position.
                {
                    // Calculate the direction vector of this particle.
                    float direction_x = -difference_x / distance;
                    float direction_y = -difference_y / distance;

                    // Add other particle's pressure and viscocity contributions using Navier-Stokes equations, scaled by Muller kernels.
                    pressure_contribution_x += direction_x * particle_mass * (pressure + other_pressure) / (2.f * other_density) * SPIKY_GRAD * powf(particle_diameter - distance, 2.f);
                    pressure_contribution_y += direction_y * particle_mass * (pressure + other_pressure) / (2.f * other_density) * SPIKY_GRAD * powf(particle_diameter - distance, 2.f);

                    viscocity_contribution_x += particle_viscocity * particle_mass * (other_vx - vx) / other_density * VISC_LAP * (particle_diameter - distance);
                    viscocity_contribution_y += particle_viscocity * particle_mass * (other_vy - vy) / other_density * VISC_LAP * (particle_diameter - distance);

                }
            }
        }
        // Calculate gravity contribution.
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

        // Handle collision detection at edges of environment rectangle.
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
    cudaGraphicsMapResources(1, cuda_vbo_resource, 0);
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

void c_computeDensityPressure(int threads, int blocks, float* positions, float* densities, float* pressures, int number_of_particles, float particle_diameter_squared, float particle_mass, float POLY6, float BOLTZMANN_CONSTANT, float REST_DENSITY)
{
    g_computeDensityPressure<<<blocks, threads>>>(positions, densities, pressures, number_of_particles, particle_diameter_squared, particle_mass, POLY6, BOLTZMANN_CONSTANT, REST_DENSITY);
    errorCheck(__LINE__, __FILE__);

    // Sync threads before continuing.
    cudaDeviceSynchronize();
    errorCheck(__LINE__, __FILE__);
}

void c_computeForces(int threads, int blocks, float* positions, float* velocities, float* forces, float* densities, float* pressures, int number_of_particles, float particle_diameter, float particle_mass, float particle_viscocity, float gravity_x, float gravity_y, float SPIKY_GRAD, float VISC_LAP)
{
    g_computeForces<<<blocks, threads>>>(positions, velocities, forces, densities, pressures, number_of_particles, particle_diameter, particle_mass, particle_viscocity, gravity_x, gravity_y, SPIKY_GRAD, VISC_LAP);
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

