#pragma once
#include "Input.h"

#include "circle.h"

using namespace sf;
using namespace std;

#define BENCHMARKING 0

#if BENCHMARKING
#include <fstream>
#endif

const float PI = 3.141592653589793239f;	// Max precision.

class QuadTree;

class Application
{
public:
	Application(RenderWindow* window);
	~Application();

public:
	// Application loop.
	void run();

private:
	void handleInput(); 
	void update();
	void render();
	
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

	RenderWindow* window_;		// SFML object for handling the creation of a window to work in.
	Input input_;				// Input handler class.
	bool running_;				// If false, the application will close.
	int environment_width_;		// Environment refers to the rectangle the particles exist in.
	int environment_height_;	// Environment refers to the rectangle the particles exist in.

	vector<Particle> particles_;	// All particles.
		
	const int number_of_particles_;		// Number of particles in the SPH simulation.
	float gravity_magnitude_;
	Vector2f gravity_;			// Force acting on every particle.
	const float damping_;				// Scaler applied to reflected velocity vectors when particles go out of bounds.
	const float integration_timestep_;	// Determines how far to progress SPH simulation each frame.

	// Shared particle attributes - all particles are considered to be identical.				
	/////////////////////////////////////////////////////////////////////////////////////		<---- consider making each variable of Particle object a seperate array to make more data oriented.
	
	const float particle_diameter_;			
	const float particle_radius_squared_;
	const float particle_mass_;				
	const float particle_viscocity_; 

	// Equation of state (Ideal gas law) constants.
	/////////////////////////////////////////////////////////////////////////////////////

	const float REST_DENSITY;			// Density of a particle at rest.
	const float BOLTZMANN_CONSTANT;		// Relates kinetic energy of particles to thermodynamic temperature of fluid.

	// Smoothing kernel constants defined by Hans-Georg Müller .
	/////////////////////////////////////////////////////////////////////////////////////

	const float POLY6;
	const float SPIKY_GRAD;
	const float VISC_LAP;

	QuadTree* quad_tree_;		// Quad Tree data structure.
	Circle circle;				// Shape used to queery Quad Tree. (used for returning particles in the circle)

#if BENCHMARKING
	std::ofstream output_file_stream_;
	const int simulations_;
	int simulation_;
	const int iterations_;
	const int start_iteration_;
	int iteration_;
	long long start_time_;
	long long time_sum_;
#endif
};

