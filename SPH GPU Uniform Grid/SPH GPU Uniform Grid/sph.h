#pragma once

#define BENCHMARKING 0

#if BENCHMARKING
#include <fstream>
#endif

class SPH
{
public:
    SPH(int number_of_particles, float dam_width, float scale);
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

    // Used for debugging to view values of host particle variable arrays.
    void debug(float* device_positions);

private:
    int number_of_particles_;		    // Number of particles in the SPH simulation.
    int dam_width_;                     // Width of dam in particles.
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

    // CPU data
    float* host_positions_;
    float* host_velocities_;
    float* host_forces_;
    float* host_densities_;
    float* host_pressures_;
    int* host_IDs_;
    int* host_hashes_;
    int* host_cell_start_indices_;
    int* host_cell_end_indices_;

    // GPU data
    float* device_velocities_;
    float* device_forces_;
    float* device_densities_;
    float* device_pressures_;
    int* device_IDs_;
    int* device_hashes_;
    int* device_cell_start_indices_;    // Hash table that stores the index of the first particle in each cell. 
    int* device_cell_end_indices_;      // Hash table that stores the index of the first particle in the next cell. 

    // CUDA vertex buffer object for particle positions.
    float* vbo_;         

    // ID of CUDA vertex buffer.
    unsigned int vbo_ID_;                               

    // CUDA VBO resource. (handles CUDA / OpenGL exchange)
    struct cudaGraphicsResource* cuda_vbo_resource_;    

    // Execution configuration for device kernels.
    int threads_;
    int blocks_;

    // Uniform Grid variables.
    const int cell_size_;           // square cells with a size equal to the particles' diameter.
    const int grid_width_;          // The number of cells on the x axis. The environment width divided by cell size.
    const int grid_height_;         // The number of cells on the y axis. The environment height divided by cell size.
    const int number_of_cells_;     // Number of cells in grid.
                                    
    // Debug.                       
    static const int number_of_particles = 8192;    // <--- dangerous! (must allways match actual number of particles)
    static const int grid_width = 128;
    static const int grid_height = 64;
    float position[number_of_particles * 2];
    float velocity[number_of_particles * 2];
    float force[number_of_particles * 2];
    float density[number_of_particles];
    float pressure[number_of_particles];
    int hash[number_of_particles];
    int ID[number_of_particles];
    int start[grid_width * grid_height];
    int end[grid_width * grid_height];

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

