//https://lucasschuermann.com/writing/implementing-sph-in-2d

#include "Application.h"

#if BENCHMARKING
#include <iostream>
#include <chrono>  
using namespace std::chrono;
#endif

Application::Application(RenderWindow* window) :
	window(window),
	running(true),
	environment_width(1024),
	environment_height(512),
	number_of_particles(256),
	gravity_magnitude(4200 * 9.8f),
	gravity(0.f, gravity_magnitude),
	damping(-0.5f),
	integration_timestep(0.0008f),
	particle_diameter(8.0f),
	particle_diameter_squared(pow(particle_diameter, 2)),
	particle_mass(50.f),
	particle_viscocity(250.f),
	REST_DENSITY(1000.f),
	BOLTZMANN_CONSTANT(2000.f),
	POLY6(315.f / (65.f * PI * pow(particle_diameter, 9.f))),
	SPIKY_GRAD(-45.f / (PI * pow(particle_diameter, 6.f))),
	VISC_LAP(45.f / (PI * pow(particle_diameter, 6.f)))
#if BENCHMARKING
    ,simulations_(10),
    simulation_(0),
    iterations_(1025),
    start_iteration_(25),
    iteration_(0),
    frames_(0),
    start_time_(0),
	time_sum_(0)
#endif
{
	srand(time(NULL));

#if BENCHMARKING
	std::ofstream temp(to_string(number_of_particles) + ".csv", std::ofstream::out);
	output_file_stream_.swap(temp);
#endif

	initializeParticles();

#if BENCHMARKING
	start_time_ = system_clock::now().time_since_epoch().count();
#endif
}

Application::~Application()
{
}

void Application::run()
{
	// Variables for calculating frames per second.
	Clock clock;
	Time fps_timer = Time::Zero;	
	int frames = 0;				

	// Application loop.
	while (running)
	{
		// Queery relevant Windows input messages.
		Event event;
		while (window->pollEvent(event))
		{
			switch (event.type)
			{
			// If the X symbol in top right is clicked, exit application.
			case Event::Closed:
				running = false;
				break;
			// If a key is pressed, tell input manager what key has been pressed.
			case Event::KeyPressed:
				input.setKeyDown(event.key.code);
				break;
			// If a key is released, tell input manager what key has been released.
			case Event::KeyReleased:
				input.setKeyUp(event.key.code);
				break;
			// If a mouse button is pressed, tell input manager what mouse button has been pressed.
			case Event::MouseButtonPressed:
				if (event.mouseButton.button == sf::Mouse::Left)
				{
					input.setMouseLeftDown(true);
				}
				break;
			// If a mouse button is released, tell input manager what mouse button has been released.
			case Event::MouseButtonReleased:
				if (event.mouseButton.button == sf::Mouse::Left)
				{
					input.setMouseLeftDown(false);
				}
				break;
			// If the mouse is moved, tell input manager the new coordinates of the mouse.
			case Event::MouseMoved:
				input.setMousePosition(event.mouseMove.x, event.mouseMove.y);
				break;
			default:
				break;
			}
		}
#if !BENCHMARKING
		handleInput();
#endif
		update();
		render();
	}
	window->close();
}

void Application::handleInput()
{
	// Manipulating gravity.
	if (input.getKeyDown('a' - 97))
	{
		gravity.x = -gravity_magnitude;
		gravity.y = 0;
	}
	if (input.getKeyDown('d' - 97))
	{
		gravity.x = gravity_magnitude;
		gravity.y = 0;
	}
	if (input.getKeyDown('w' - 97))
	{
		gravity.x = 0;
		gravity.y = -gravity_magnitude;
	}
	if (input.getKeyDown('s' - 97))
	{
		gravity.x = 0;
		gravity.y = gravity_magnitude;
	}
	
}

void Application::update()
{
#if !BENCHMARKING

	// SPH simulation.
	computeDensityPressure();
	computeForces();
	integrate();

#else
	auto start = steady_clock::now();
	computeDensityPressure();
	computeForces();
	integrate();
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
		time_sum_ /= frames_;
		output_file_stream_ << time_sum_ << ",";
		long long end_time = system_clock::now().time_since_epoch().count();
		float seconds = (float)(end_time - start_time_) / 10000000;
		long long fps = (iterations_ / seconds);
		output_file_stream_ << "," << fps << ",\n";
		iteration_ = 0;
		simulation_++;
		time_sum_ = 0;
		frames_ = 0;
		reset();
		start_time_ = system_clock::now().time_since_epoch().count();
	}
	frames_++;
#endif
}

void Application::render()
{
	window->clear(Color::Black);

	// Render Particles.
	for (Particle& particle : particles)
	{
		RectangleShape square;
		square.setFillColor(Color(0, 127, 255));
		square.setSize(Vector2f(4, 4));
		square.setPosition(particle.position);
		square.setOrigin(Vector2f(2, 2));
		window->draw(square);
	}

	window->display();
}

void Application::initializeParticles()
{
	// Loop up from the bottom of the screen until all particles have been initialized.
 	for (float y = environment_height; ; y -= particle_diameter)
	{
		// Loop the width of the dam.
		for (float x = 0; x < 8 * particle_diameter; x += particle_diameter)
		{
			if (particles.size() < number_of_particles)
			{
				// A small offest applied to X axis on initialization of particles.
				float jitter = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);		
				particles.push_back(Particle(x + jitter, y));
			}
			else
			{
				return;
			}
		}
	}
}

void Application::edges(Particle& particle)
{
	// A small offest applied to particles on boundary collision to avoid stacking.
	float jitter = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

	// Handle collision detection at edges of environment rectangle.
	// LEFT
	if (particle.position.x < 0.0f)
	{
		particle.velocity.x *= damping;
		particle.position.x = jitter;
	}
	// RIGHT
	if (particle.position.x > environment_width - 1)
	{
		particle.velocity.x *= damping;
		particle.position.x = environment_width - 1 - jitter;
	}
	// TOP
	if (particle.position.y < 0.0f)
	{
		particle.velocity.y *= damping;
		particle.position.y = jitter;
	}
	// BOTTOM
	if (particle.position.y > environment_height - 1)
	{
		particle.velocity.y *= damping;
		particle.position.y = environment_height - 1 - jitter;
	}
}

void Application::computeDensityPressure()
{
	for (Particle& particle : particles)
	{
		// Reset particles density so it can be recalculated.
		particle.density = 0.f;
		for (Particle& other_particle : particles)
		{
			// Calculate vector between the pair of particles.
			Vector2f difference = other_particle.position - particle.position;
			// Calculate the squared distance between the pair of particles. 
			float distance_squared = pow(difference.x, 2) + pow(difference.y, 2);
			// If the particles are overlapping.
			if (distance_squared < particle_diameter_squared)
			{
				// Add the other particle's mass, scaled by a muller kernel to this particle's density.
				particle.density += particle_mass * POLY6 * pow(particle_diameter_squared - distance_squared, 3.f);
			}
		}
		// Calculate pressure of particle using an "equation of state", relating it's density to a given rest density.
		particle.pressure = BOLTZMANN_CONSTANT * (particle.density - REST_DENSITY);
	}
}

void Application::computeForces()
{
	for (Particle& particle : particles)
	{
		// Reset pressure and viscocity calcualtions so they can be recalculated.
		Vector2f pressure_contribution(0.f, 0.f);
		Vector2f viscocity_contribution(0.f, 0.f);
		for (Particle& other_particle : particles)
		{
			// Do not compare particles if they are the same.
			if (&particle != &other_particle)
			{
				// Calculate vector between the pair of particles.
				Vector2f difference = other_particle.position - particle.position;

				// Calculate the distance between the pair of particles. Distance needed later so no benefit from comparing squared distances.
				float distance = sqrt(pow(difference.x, 2) + pow(difference.y, 2));
				// If particles are overlapping.
				if (distance != 0 && distance < particle_diameter)	// <--- distance != 0 is handling case where particles are at the same position.
				{
					// Calculate the desired direction vector of this particle. (Negative distance, normalized.)
					Vector2f desired = -Vector2f(difference.x / distance, difference.y / distance);
					// Add other particle's pressure and viscocity contributions using Navier-Stokes equations, scaled by Muller kernels.
					pressure_contribution += desired * particle_mass * (particle.pressure + other_particle.pressure) / (2 * other_particle.density) * SPIKY_GRAD * pow(particle_diameter - distance, 2.f);
					viscocity_contribution += particle_viscocity * particle_mass * (other_particle.velocity - particle.velocity) / other_particle.density * VISC_LAP * (particle_diameter - distance);
				}
			}
		}
		// Calculate grivity contributions by multiplying gravity by this particle's calculated density.
		Vector2f gravity_contribution = gravity * particle.density;
		// Add all force contributions together to calculate the particle's force.
		particle.force = pressure_contribution + viscocity_contribution + gravity_contribution;
	}
}

void Application::integrate()
{
	// Update velocity and positions of particles.
	for (Particle& particle : particles)
	{
		particle.velocity += particle.force / particle.density * integration_timestep;
		particle.position += particle.velocity * integration_timestep;
		// Handle collision detection at edges of environment rectangle.
		edges(particle);
	}
}

void Application::reset()
{
	// Reset simulation.
	particles.clear();
	initializeParticles();
}


