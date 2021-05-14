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
	};

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
	// Reset simulation.
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
	/////////////////////////////////////////////////////////////////////////////////////		<---- consider making each variable of Particle object a seperate array to make more data oriented.
	
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

