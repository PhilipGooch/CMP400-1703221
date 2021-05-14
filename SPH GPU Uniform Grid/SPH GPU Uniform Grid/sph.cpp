#include "sph.h"
#include "sph.cuh"
#include "thrust.cuh"

#include <glew.h>       // glGenBuffers() glBindBuffer() glBufferData() glBufferSubData()

#include <memory.h>     // memset()
#include <random>       // rand()
#include <time.h>       // time()

#if BENCHMARKING
#include <iostream>
#include <chrono>  
using namespace std::chrono;
#endif

#define PI 3.141592654589793239f

SPH::SPH(int number_of_particles, float dam_width, float scale) :
    number_of_particles_(number_of_particles),
    dam_width_(dam_width),
    scale_(scale),
    gravity_(4200.0f * 9.8f),
    gravity_x_(0.f),
    gravity_y_(-gravity_),
    damping_(-0.5f),
    integration_timestep_(0.0008f),
    particle_diameter_(8.0f),
    particle_diameter_squared_(pow(particle_diameter_, 2)),
    particle_mass_(50.f),
    particle_viscocity_(250.f),
    REST_DENSITY_(1000.f),
    BOLTZMANN_CONSTANT_(2000.f),
    POLY6_(315.f / (65.f * PI * pow(particle_diameter_, 9.f))),
    SPIKY_GRAD_(-45.f / (PI * pow(particle_diameter_, 6.f))),
    VISC_LAP_(45.f / (PI * pow(particle_diameter_, 6.f))),
    cell_size_(particle_diameter_),
    grid_width_((1024 * scale_) / cell_size_),
    grid_height_((512 * scale_) / cell_size_),
    number_of_cells_(grid_width_* grid_height_)
#if BENCHMARKING
    ,simulations_(10),
    simulation_(0),
    iterations_(1000),
    start_iteration_(25),
    iteration_(0),
    frames_(0),
    start_time_(0),
    time_sum_(0)
#endif
{
    srand(time(NULL));

    // Execution configuration for device kernels.
    threads_ = std::min(number_of_particles_, 256);
    blocks_ = number_of_particles_ % threads_ == 0 ? number_of_particles_ / threads_ : number_of_particles_ / threads_ + 1;

    // Initialize host variables.
    host_positions_ = new float[number_of_particles * 2];
    host_velocities_ = new float[number_of_particles * 2];
    host_forces_ = new float[number_of_particles * 2];
    host_densities_ = new float[number_of_particles];
    host_pressures_ = new float[number_of_particles];
    host_hashes_ = new int[number_of_particles];
    host_IDs_ = new int[number_of_particles];
    host_cell_start_indices_ = new int[grid_width_ * grid_height_];
    host_cell_end_indices_ = new int[grid_width_ * grid_height_];
    memset(host_positions_, 0, number_of_particles * sizeof(float) * 2);
    memset(host_velocities_, 0, number_of_particles * sizeof(float) * 2);
    memset(host_forces_, 0, number_of_particles * sizeof(float) * 2);
    memset(host_densities_, 0, number_of_particles * sizeof(float));
    memset(host_pressures_, 0, number_of_particles * sizeof(float));
    memset(host_hashes_, 0, number_of_particles * sizeof(int));
    memset(host_IDs_, 0, number_of_particles * sizeof(int));
    memset(host_cell_start_indices_, 0, grid_width_ * grid_height_ * sizeof(int));
    memset(host_cell_end_indices_, 0, grid_width_ * grid_height_ * sizeof(int));

    // Create OpenGL Vertex Buffer Object for device positions.
    glGenBuffers(1, &vbo_ID_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_ID_);
    glBufferData(GL_ARRAY_BUFFER, number_of_particles * sizeof(float) * 2, 0, GL_DYNAMIC_DRAW);

    // Register VBO with CUDA.
    c_registerVBO(&cuda_vbo_resource_, vbo_ID_);

    // Initialize Device variables.
    c_cudaMalloc((void**)&device_velocities_, number_of_particles * sizeof(float) * 2);
    c_cudaMalloc((void**)&device_forces_, number_of_particles * sizeof(float) * 2);
    c_cudaMalloc((void**)&device_densities_, number_of_particles * sizeof(float));
    c_cudaMalloc((void**)&device_pressures_, number_of_particles * sizeof(float));
    c_cudaMalloc((void**)&device_IDs_, number_of_particles * sizeof(int));
    c_cudaMalloc((void**)&device_hashes_, number_of_particles * sizeof(int));
    c_cudaMalloc((void**)&device_cell_start_indices_, grid_width_ * grid_height_ * sizeof(int));
    c_cudaMalloc((void**)&device_cell_end_indices_, grid_width_ * grid_height_ * sizeof(int));

    resetParticles();

#if BENCHMARKING
    std::ofstream temp(std::to_string(number_of_particles_) + ".csv", std::ofstream::out);
    output_file_stream_.swap(temp);
    start_time_ = system_clock::now().time_since_epoch().count();
#endif
}

SPH::~SPH()
{
    // Delete dynamically allocated arrays.
    delete[] host_positions_;
    delete[] host_velocities_;
    delete[] host_forces_;
    delete[] host_densities_;
    delete[] host_pressures_;
    delete[] host_IDs_;
    delete[] host_hashes_;
    delete[] host_cell_start_indices_;
    delete[] host_cell_end_indices_;

    // Release device memory.
    c_cudaFree(device_velocities_);
    c_cudaFree(device_forces_);
    c_cudaFree(device_densities_);
    c_cudaFree(device_pressures_);
    c_cudaFree(device_IDs_);
    c_cudaFree(device_hashes_);
    c_cudaFree(device_cell_start_indices_);
    c_cudaFree(device_cell_end_indices_);

    // Release vertex buffer object.
    c_unregisterVBO(cuda_vbo_resource_);
    glDeleteBuffers(1, &vbo_ID_);
}

void SPH::setGravity(DIRECTION direction)
{
#if !BENCHMARKING
    // Manipulating gravity depending on user input.
    switch (direction)
    {
    case DIRECTION::LEFT:
        gravity_x_ = -gravity_;
        gravity_y_ = 0.0f;
        break;
    case DIRECTION::RIGHT:
        gravity_x_ = gravity_;
        gravity_y_ = 0.0f;
        break;
    case DIRECTION::UP:
        gravity_x_ = 0.0f;
        gravity_y_ = gravity_;
        break;
    case DIRECTION::DOWN:
        gravity_x_ = 0.0f;
        gravity_y_ = -gravity_;
        break;
    }
#endif
}

void SPH::update()
{
    // Map VBO to access device positions.
    float *device_positions = (float *) c_mapVBO(&cuda_vbo_resource_);

#if !BENCHMARKING

    // Calculate particle hash values. (the linear cell index of the cell it is in)
    c_calculateHashValues(threads_, blocks_, number_of_particles_, device_hashes_, device_positions, grid_width_, cell_size_, scale_);

    // Sort particles based on hash values.
    c_sortParticles(device_hashes_, device_IDs_, number_of_particles_, device_positions, device_velocities_, device_forces_, device_densities_, device_pressures_);
    
    // Store indices of first and last particle in cell in hash table. 
    c_findCellStartEndIndices(threads_, blocks_, number_of_particles_, device_cell_start_indices_, device_cell_end_indices_, device_hashes_);

    // Simulation functions.
    c_computeDensityPressure(threads_, blocks_, device_positions, device_densities_, device_pressures_, number_of_particles_, particle_diameter_squared_, particle_mass_, POLY6_, BOLTZMANN_CONSTANT_, REST_DENSITY_, cell_size_, grid_width_, grid_height_, device_cell_start_indices_, device_cell_end_indices_, scale_);
    c_computeForces(threads_, blocks_, device_positions, device_velocities_, device_forces_, device_densities_, device_pressures_, number_of_particles_, particle_diameter_, particle_mass_, particle_viscocity_, gravity_x_, gravity_y_, SPIKY_GRAD_, VISC_LAP_, cell_size_, grid_width_, grid_height_, device_cell_start_indices_, device_cell_end_indices_, device_IDs_, scale_);
    c_integrate(threads_, blocks_, scale_, device_positions, device_velocities_, device_forces_, device_densities_, number_of_particles_, integration_timestep_, damping_);

#else
    auto start = steady_clock::now();
    c_calculateHashValues(threads_, blocks_, number_of_particles_, device_hashes_, device_positions, grid_width_, cell_size_, scale_);
    c_sortParticles(device_hashes_, device_IDs_, number_of_particles_, device_positions, device_velocities_, device_forces_, device_densities_, device_pressures_);
    c_findCellStartEndIndices(threads_, blocks_, number_of_particles_, device_cell_start_indices_, device_cell_end_indices_, device_hashes_);
    c_computeDensityPressure(threads_, blocks_, device_positions, device_densities_, device_pressures_, number_of_particles_, particle_diameter_squared_, particle_mass_, POLY6_, BOLTZMANN_CONSTANT_, REST_DENSITY_, cell_size_, grid_width_, grid_height_, device_cell_start_indices_, device_cell_end_indices_, scale_);
    c_computeForces(threads_, blocks_, device_positions, device_velocities_, device_forces_, device_densities_, device_pressures_, number_of_particles_, particle_diameter_, particle_mass_, particle_viscocity_, gravity_x_, gravity_y_, SPIKY_GRAD_, VISC_LAP_, cell_size_, grid_width_, grid_height_, device_cell_start_indices_, device_cell_end_indices_, device_IDs_, scale_);
    c_integrate(threads_, blocks_, scale_, device_positions, device_velocities_, device_forces_, device_densities_, number_of_particles_, integration_timestep_, damping_);
    auto end = steady_clock::now();
    auto elapsed = duration_cast<nanoseconds>(end - start).count();
    if (iteration_ >= start_iteration_)
    {
        time_sum_ += elapsed;
    }
    if (simulation_ == simulations_)
    {
        output_file_stream_.close();
        exit(EXIT_SUCCESS);
    }
    if (++iteration_ == iterations_ + start_iteration_)
    {
        std::cout << simulation_ << std::endl;
        time_sum_ /= iterations_;
        output_file_stream_ << time_sum_ << ",";
        long long end_time = system_clock::now().time_since_epoch().count();
        float seconds = (float)(end_time - start_time_) / 10000000;
        long long fps = (iterations_ / seconds);
        output_file_stream_ << seconds  << "," << fps << ",\n";
        iteration_ = 0;
        simulation_++;
        time_sum_ = 0;
        frames_ = 0;
        resetParticles();
        start_time_ = system_clock::now().time_since_epoch().count();
    }
    frames_++;
#endif

    // Unmap ready to be mapped again next iteration.
    c_unmapVBO(cuda_vbo_resource_);
}

void SPH::resetParticles()
{
    // Width of dam in particles.
    for (int i = 0; i < number_of_particles_; i++)
    {
        // A small offest applied to x axis on initialization of particles.
        float jitter = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        host_positions_[i * 2] = (-512.0f * scale_) + (i % dam_width_) * particle_diameter_ + jitter; // x
        host_positions_[i * 2 + 1] = (-256.0f * scale_) + (i / dam_width_) * particle_diameter_;      // y
        host_velocities_[i * 2] = 0.0f;         // x
        host_velocities_[i * 2 + 1] = 0.0f;     // y
        host_forces_[i * 2] = 0.0f;             // x
        host_forces_[i * 2 + 1] = 0.0f;         // y
        host_densities_[i] = 0.0f;
        host_pressures_[i] = 0.0f;
        host_hashes_[i] = 0;
        host_IDs_[i] = i;
    }
    
    // Set data in VBO to host positions.
    glBufferSubData(GL_ARRAY_BUFFER, 0, number_of_particles_ * sizeof(float) * 2, host_positions_);
    
    // Initialize device memory.
    c_cudaMemcpyHostToDevice(device_velocities_, host_IDs_, number_of_particles_ * sizeof(int) * 2);
    c_cudaMemcpyHostToDevice(device_forces_, host_forces_, number_of_particles_ * sizeof(int) * 2);
    c_cudaMemcpyHostToDevice(device_densities_, host_densities_, number_of_particles_ * sizeof(int));
    c_cudaMemcpyHostToDevice(device_pressures_, host_pressures_, number_of_particles_ * sizeof(int));
    c_cudaMemcpyHostToDevice(device_IDs_, host_IDs_, number_of_particles_ * sizeof(int));
    c_cudaMemcpyHostToDevice(device_hashes_, host_hashes_, number_of_particles_ * sizeof(int));
    c_cudaMemcpyHostToDevice(device_cell_start_indices_, host_cell_start_indices_, grid_width_ * grid_height_ * sizeof(int));
    c_cudaMemcpyHostToDevice(device_cell_end_indices_, host_cell_end_indices_, grid_width_ * grid_height_ * sizeof(int));
}

void SPH::debug(float* device_positions)
{
    // Copy arrays from device to host.
    c_cudaMemcpyDeviceToHost(host_positions_, device_positions, 0, number_of_particles_ * sizeof(float) * 2);
    c_cudaMemcpyDeviceToHost(host_velocities_, device_velocities_, 0, number_of_particles_ * sizeof(float) * 2);
    c_cudaMemcpyDeviceToHost(host_forces_, device_forces_, 0, number_of_particles_ * sizeof(float) * 2);
    c_cudaMemcpyDeviceToHost(host_densities_, device_densities_, 0, number_of_particles_ * sizeof(float));
    c_cudaMemcpyDeviceToHost(host_pressures_, device_pressures_, 0, number_of_particles_ * sizeof(float));
    c_cudaMemcpyDeviceToHost(host_IDs_, device_IDs_, 0, number_of_particles_ * sizeof(float));
    c_cudaMemcpyDeviceToHost(host_hashes_, device_hashes_, 0, number_of_particles_ * sizeof(float));
    c_cudaMemcpyDeviceToHost(host_cell_start_indices_, device_cell_start_indices_, 0, grid_height_ * grid_width_ * sizeof(float));
    c_cudaMemcpyDeviceToHost(host_cell_end_indices_, device_cell_end_indices_, 0, grid_height_ * grid_width_ * sizeof(float));

    // Copy arrays into defined c-style arrays for inspection when debugging.
    memcpy(&position[0], host_positions_, number_of_particles_ * sizeof(float) * 2);
    memcpy(&velocity[0], host_velocities_, number_of_particles_ * sizeof(float) * 2);
    memcpy(&force[0], host_forces_, number_of_particles_ * sizeof(float) * 2);
    memcpy(&density[0], host_densities_, number_of_particles_ * sizeof(float));
    memcpy(&pressure[0], host_pressures_, number_of_particles_ * sizeof(float));
    memcpy(&ID[0], host_IDs_, number_of_particles_ * sizeof(float));
    memcpy(&hash[0], host_hashes_, number_of_particles_ * sizeof(float));
    memcpy(&start[0], host_cell_start_indices_, grid_height_ * grid_width_ * sizeof(float));
    memcpy(&end[0], host_cell_end_indices_, grid_height_ * grid_width_ * sizeof(float));
}

