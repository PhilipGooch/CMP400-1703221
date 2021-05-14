#pragma once
#include "Input.h"

using namespace sf;
using namespace std;

const float PI = 3.141592653589793239f;	// Max precision.

#define BENCHMARKING 0

#if BENCHMARKING
#include <fstream>
#endif

class Application
{
public:
	Application(RenderWindow* window);
	~Application();

public:
	// Application loop.
	void run();

private:

	struct Particle {					
		Particle(float _x, float _y) : position(_x, _y), velocity(0.f, 0.f), force(0.f, 0.f), density(0), pressure(0.f) {}
		Vector2f position, velocity, force;
		float density, pressure;
		int ID;						// Giving particles an ID as they get shuffled and need to identify them.
		int hash;					// Hash is the linear cell index from the 2D cell array. 
	};

	void handleInput(); 
	void update();
	void render();

	// Uniform Grid functions.
	/////////////////////////////////////////////////////////////////////////////////////

	// Calculate x and y position the of cell the particle is in.
	Vector2i calculateGridPosition(Vector2f position);

	// Calculate the linear cell index from the 2D cell array.
	int calculateGridHash(Vector2i grid_position);

	// Calculate the hash values of the particles.
	void calculateHashValues();

	// Sort particles based on low to high hash values.
	void sortParticles();

	// Store the indices of the first particle in each cell, in cell start hash table.
	// Store the indices of the first particle in the next cell, in cell end hash table.
	void findCellStartEndIndices();

	// SPH Functions.
	/////////////////////////////////////////////////////////////////////////////////////

	// Initializes particles in the form of a dam.
	void initializeParticles();			
	// Logic for particles when they go out of bounds.
	void edges(Particle& particle);	
	// Compute the density and pressure of all particles.
	void computeDensityPressure();		
	// Compute the force each particle is applying.
	void computeForces();	
	// Update positions of particles.
	void integrate();
	// Reset particles.
	void reset();						

	RenderWindow* window;		// SFML object for handling the creation of a window to work in.
	Input input;				// Input handler class.
	bool running;				// If false, the application will close.
	int environment_width;		// Environment refers to the rectangle the particles exist in.
	int environment_height;		// Environment refers to the rectangle the particles exist in.

	vector<Particle> particles;		// All particles.
		
	const int number_of_particles;		// Number of particles in the SPH simulation.
	const float gravity_magnitude;		// Used to set gravity vector in different direcctions.
	Vector2f gravity;					// Force acting on every particle.
	const float damping;				// Scaler applied to reflected velocity vectors when particles go out of bounds.
	const float integration_timestep;	// Determines how far to progress SPH simulation each frame.

	// Shared particle attributes - all particles are considered to be identical.				
	/////////////////////////////////////////////////////////////////////////////////////
	
	const float particle_diameter;			
	const float particle_diameter_squared;
	const float particle_mass;				
	const float particle_viscocity; 

	// Equation of state (Ideal gas law) constants
	/////////////////////////////////////////////////////////////////////////////////////

	const float REST_DENSITY;			// Density of a particle at rest.
	const float BOLTZMANN_CONSTANT;		// Relates kinetic energy of particles to thermodynamic temperature of fluid.

	// Smoothing kernel constants defined by Hans-Georg Müller 
	/////////////////////////////////////////////////////////////////////////////////////

	const float POLY6;
	const float SPIKY_GRAD;
	const float VISC_LAP;

	// Uniform Grid 
	/////////////////////////////////////////////////////////////////////////////////////

	const int cell_size;				// square cells with a size equal to the particles' diameter.
	const int grid_width;				// The number of cells on the x axis. The environment width divided by cell size.
	const int grid_height;				// The number of cells on the y axis. The environment height divided by cell size.
	const int number_of_cells;			// Number of cells in grid.
	vector<int> cell_start_indices;		// Hash table that stores the index of the first particle in each cell. 
	vector<int> cell_end_indices;		// Hash table that stores the index of the first particle in the next cell. 

#if BENCHMARKING
	std::ofstream output_file_stream_;
	const int simulations_;
	int simulation_;
	const int iterations_;
	const int start_iteration_;
	int iteration_;
	long long frames_;
	long long start_time_;
	long long computeDensityPressure_sum_;
	long long computeForces_sum_;
	long long integrate_sum_;
	long long time_sum_;
#endif
};

