#include "sph.h"
#include "sph.cuh"

#include <glew.h>

#include <memory.h>     // memset()
#include <random>       // rand()
#include <time.h>       // time()

#if BENCHMARKING
#include <iostream>
#include <chrono>  
using namespace std::chrono;
#endif

#define PI 3.141592654589793239f

SPH::SPH(int number_of_particles, float scale) :
    number_of_particles_(number_of_particles),
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
    scale_(scale)
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
    memset(host_positions_, 0, number_of_particles * sizeof(float) * 2);
    memset(host_velocities_, 0, number_of_particles * sizeof(float) * 2);
    memset(host_forces_, 0, number_of_particles * sizeof(float) * 2);
    memset(host_densities_, 0, number_of_particles * sizeof(float));
    memset(host_pressures_, 0, number_of_particles * sizeof(float));

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

    // Release device memory.
    c_cudaFree(device_velocities_);
    c_cudaFree(device_forces_);
    c_cudaFree(device_densities_);
    c_cudaFree(device_pressures_);

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

    // Simulation functions.
    c_computeDensityPressure(threads_, blocks_, device_positions, device_densities_, device_pressures_, number_of_particles_, particle_diameter_squared_, particle_mass_, POLY6_, BOLTZMANN_CONSTANT_, REST_DENSITY_);
    c_computeForces(threads_, blocks_, device_positions, device_velocities_, device_forces_, device_densities_, device_pressures_, number_of_particles_, particle_diameter_, particle_mass_, particle_viscocity_, gravity_x_, gravity_y_, SPIKY_GRAD_, VISC_LAP_);
    c_integrate(threads_, blocks_, scale_, device_positions, device_velocities_, device_forces_, device_densities_, number_of_particles_, integration_timestep_, damping_);

#else
    auto start = steady_clock::now();
    c_computeDensityPressure(threads_, blocks_, device_positions, device_densities_, device_pressures_, number_of_particles_, particle_diameter_squared_, particle_mass_, POLY6_, BOLTZMANN_CONSTANT_, REST_DENSITY_);
    c_computeForces(threads_, blocks_, device_positions, device_velocities_, device_forces_, device_densities_, device_pressures_, number_of_particles_, particle_diameter_, particle_mass_, particle_viscocity_, gravity_x_, gravity_y_, SPIKY_GRAD_, VISC_LAP_);
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
        time_sum = 0;
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
    int width = 64;
    for (int i = 0; i < number_of_particles_; i++)
    {
        // A small offest applied to x axis on initialization of particles.
        float jitter = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        host_positions_[i * 2] = (-512.0f * scale_) + (i % width) * particle_diameter_ + jitter; // x
        host_positions_[i * 2 + 1] = (-256.0f * scale_) + (i / width) * particle_diameter_;      // y
        host_velocities_[i * 2] = 0.0f;         // x
        host_velocities_[i * 2 + 1] = 0.0f;     // y
        host_forces_[i * 2] = 0.0f;             // x
        host_forces_[i * 2 + 1] = 0.0f;         // y
        host_densities_[i] = 0.0f;
        host_pressures_[i] = 0.0f;
    }
    // Set data in VBO to host positions.
    glBufferSubData(GL_ARRAY_BUFFER, 0, number_of_particles_ * sizeof(float) * 2, host_positions_);
}

