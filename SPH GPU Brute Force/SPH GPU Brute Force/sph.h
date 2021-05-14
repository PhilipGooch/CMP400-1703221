#pragma once

#define BENCHMARKING 0

#if BENCHMARKING
#include <fstream>
#endif

struct Particle
{
    float x, y, vx, vy, fx, fy, density, pressure;
};

class SPH
{
public:
    SPH(int number_of_particles, float scale);
    ~SPH();

    // Directions for gravity.
    enum DIRECTION
    {
        LEFT,
        RIGHT,
        UP,
        DOWN
    };

    void update();

    void setGravity(DIRECTION direction);

private: 
    void resetParticles();

private:
    int number_of_particles_;		    // Number of particles in the SPH simulation.
    float scale_;                       // Manipulates environment size and zoom of camera.

    float gravity_;                     // Magnitude of gravity vector.
    float gravity_x_;					// Force acting on every particle.
    float gravity_y_;
    float damping_;					    // Scaler applied to reflected velocity vectors when particles go out of bounds.
    float integration_timestep_;		// Determines how far to progress SPH simulation each frame.

    // Shared particle attributes. All particles are considered to be identical.	
    float particle_diameter_;
    float particle_diameter_squared_;
    float particle_mass_;
    float particle_viscocity_;

    // Equation of state (Ideal gas law) constants.
    float REST_DENSITY_;				// Density of a particle at rest.
    float BOLTZMANN_CONSTANT_;		    // Relates kinetic energy of particles to thermodynamic temperature of fluid.

    // Smoothing kernel constants defined by Hans-Georg Müller.
    float POLY6_;
    float SPIKY_GRAD_;
    float VISC_LAP_;

    // CPU data.
    float* host_positions_;             
    float* host_velocities_;            
    float* host_forces_;
    float* host_densities_;
    float* host_pressures_;

    // GPU data.
    float* device_velocities_;
    float* device_forces_; 
    float* device_densities_;
    float* device_pressures_;

    // CUDA vertex buffer object for particle positions.
    float *vbo_;         

    // ID of CUDA vertex buffer.
    unsigned int vbo_ID_;                               

    // CUDA VBO resource. (handles CUDA / OpenGL exchange)
    struct cudaGraphicsResource* cuda_vbo_resource_;    

    // Execution configuration for device kernels.
    int threads_;
    int blocks_;

#if BENCHMARKING
    std::ofstream output_file_stream_;
    const int simulations_;
    int simulation_;
    const int iterations_;
    const int start_iteration_;
    int iteration_;
    long long frames_;
    long long start_time_;
    long long time_sum_;
#endif
};

