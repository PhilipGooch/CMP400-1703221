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
	number_of_particles(1024),
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
	VISC_LAP(45.f / (PI * pow(particle_diameter, 6.f))),
	cell_size(particle_diameter),
	grid_width((environment_width / cell_size)),
	grid_height((environment_height / cell_size)),
	number_of_cells(grid_width* grid_height)
#if BENCHMARKING
	,simulations_(10),
	simulation_(0),
	iterations_(1000),
	start_iteration_(25),
	iteration_(0),
	frames_(0),
	start_time_(0),
	computeDensityPressure_sum_(0),
	computeForces_sum_(0),
	integrate_sum_(0),
	time_sum_(0)
#endif
{
	cell_start_indices = vector<int>(number_of_cells);
	cell_end_indices = vector<int>(number_of_cells);

	initializeParticles();

#if BENCHMARKING
	std::ofstream temp(to_string(number_of_particles) + ".csv", std::ofstream::out);
	output_file_stream_.swap(temp);
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

	// Calculate particle hash values. (the linear cell index of the cell it is in)
	calculateHashValues();

	// Sort particles based on hash values.
	sortParticles();

	// Store indices of first and last particle in cell in hash table. 
	findCellStartEndIndices();

	// SPH simulation.
	computeDensityPressure();
	computeForces();
	integrate();

#else
	auto start = steady_clock::now();
	calculateHashValues();
	sortParticles();
	findCellStartEndIndices();
	computeDensityPressure();
	computeForces();
	integrate();
	auto end = steady_clock::now();
	auto elapsed = duration_cast<nanoseconds>(end - start).count();
	if (iteration_ >= start_iteration_) time_sum_ += elapsed;
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

	// Render particles.
	for (Particle& particle : particles)
	{
		RectangleShape square;
		square.setFillColor(Color(0, 127, 255));
		square.setSize(Vector2f(4, 4));
		square.setPosition(particle.position);
		square.setOrigin(Vector2f(2, 2));
		window->draw(square);
	}

	for (int i = 0; i <= grid_width; i++)
	{
		sf::Vertex line[] =
		{
			Vertex(Vector2f(i * cell_size, 0), Color(50, 50, 50)),
			Vertex(Vector2f(i * cell_size, environment_height), Color(50, 50, 50))
		};
		window->draw(line, 2, sf::Lines);
	}

	for (int i = 0; i <= grid_height; i++)
	{
		sf::Vertex line[] =
		{
			Vertex(Vector2f(0, i * cell_size), Color(50, 50, 50)),
			Vertex(Vector2f(environment_width, i * cell_size), Color(50, 50, 50))
		};
		window->draw(line, 2, sf::Lines);
	}

	window->display();
}

Vector2i Application::calculateGridPosition(Vector2f position)
{
	// Calculate x and y position of the cell the particle is in.
	return Vector2i((int)(position.x / cell_size), (int)(((environment_height - 1) - position.y) / cell_size));
}

int Application::calculateGridHash(Vector2i grid_position)
{
	// Calculate the linear cell index from the 2D cell array.
	return (grid_position.y * grid_width) + grid_position.x;
}

void Application::calculateHashValues()
{
	// Calculate hash values for each particle.
	for (Particle& particle : particles)
	{
		Vector2i grid_position = calculateGridPosition(particle.position);
		particle.hash = calculateGridHash(grid_position);
	}
}

void Application::sortParticles()
{
	// Sort particles based on low to high hash values.
	std::sort(particles.begin(), particles.end(), [](const Particle& a, const Particle& b) { return a.hash < b.hash; });
}

void Application::findCellStartEndIndices()	
{
	// Set the cell start hash table at the first particle's cell index (hash value) to be the first particle's index (0).
	cell_start_indices[particles[0].hash] = 0;

	// loop through the rest of the particles.
	for (int i = 1; i < particles.size(); i++)
	{
		// If this particle is in a difference cell to the previous particle.
		if (particles[i].hash != particles[i - 1].hash)
		{
			// Set the cell end hash table at the previous particle's cell index (hash value) to be the current particle's index.
			cell_end_indices[particles[i - 1].hash] = i;
			// Set the cell start hash table at this particle's cell index (hash value) to be the current particle's index.
			cell_start_indices[particles[i].hash] = i;
		}
	}
	// Set the cell end hash table at the final particle's cell index (hash value) to be the last particles index + 1 (number of particles).
	cell_end_indices[particles[number_of_particles - 1].hash] = number_of_particles;
}

void Application::initializeParticles()
{
	int ID = 0;
	// Loop up from the bottom of the screen until all particles have been initialized.
 	for (float y = environment_height; ; y -= particle_diameter)
	{
		// Loop the width of the dam.
		for (float x = 0; x < 32 * particle_diameter; x += particle_diameter)
		{
			if (particles.size() < number_of_particles)
			{
				// A small offest applied to X axis on initialization of particles.
				float jitter = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);		
				particles.push_back(Particle(x + jitter, y));
				particles.back().ID = ID++;							
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
		// Calculate the grid position of the particle.
		Vector2i grid_position = calculateGridPosition(particle.position);

		// Reset particle's density so it can be recalculated.
		particle.density = 0.f;

		// Loop over surrounding 8 cells and itself.
		for (int y = -1; y <= 1; y++)
		{
			for (int x = -1; x <= 1; x++)
			{
				// Calculate the grid position of the neighboring cell.
				Vector2i neighboring_grid_position = grid_position + Vector2i(x, y);

				// If the cell is on an edge, ignore the cell positions outside the grid.
				if (neighboring_grid_position.x < 0 || neighboring_grid_position.x >= grid_width ||	
					neighboring_grid_position.y < 0 || neighboring_grid_position.y >= grid_height)
				{
					continue;
				}

				// Calculate hash value of neighboring cell.
				int neighboring_grid_hash = calculateGridHash(neighboring_grid_position);

				// Get the index of the first particle in this cell.
				int start_index = cell_start_indices[neighboring_grid_hash];

				// Get the index of the first particle in the next cell.
				int end_index = cell_end_indices[neighboring_grid_hash];

				// Loop over particles in cell.
				for (int i = start_index; i < end_index; i++)	
				{
					Particle other_particle = particles[i];

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
		// Calculate the grid position of the particle.
		Vector2i grid_position = calculateGridPosition(particle.position);

		// Reset pressure and viscocity calcualtions so they can be recalculated.
		Vector2f pressure_contribution(0.f, 0.f);
		Vector2f viscocity_contribution(0.f, 0.f);

		// Loop over surrounding 8 cells and itself.
		for (int y = -1; y <= 1; y++)
		{
			for (int x = -1; x <= 1; x++)
			{
				// Calculate the grid position of the neighboring cell.
				Vector2i neighboring_grid_position = grid_position + Vector2i(x, y);

				// If the cell is on an edge, ignore the cell positions outside the grid.
				if (neighboring_grid_position.x < 0 || neighboring_grid_position.x >= grid_width ||	
					neighboring_grid_position.y < 0 || neighboring_grid_position.y >= grid_height)
				{
					continue;
				}

				// Calculate hash value of neighboring cell.
				int neighboring_grid_hash = calculateGridHash(neighboring_grid_position);

				// Get the index of the first particle in this cell.
				int start_index = cell_start_indices[neighboring_grid_hash];

				// Get the index of the first particle in the next cell.
				int end_index = cell_end_indices[neighboring_grid_hash];

				// Loop over particles in cell.
				for (int i = start_index; i < end_index; i++)	
				{
					// If it is not this particle.
					if (particles[i].ID != particle.ID)
					{
						Particle& other_particle = particles[i];

						// Calculate vector between the pair of particles.
						Vector2f difference = other_particle.position - particle.position;

						// Calculate the distance between the pair of particles. Distance needed later so no benefit from comparing squared distances.
						float distance = sqrt(pow(difference.x, 2) + pow(difference.y, 2));

						// If particles are overlapping.
						if (distance != 0 && distance < particle_diameter)	// <--- distance != 0 is handling case where particles are at the same position.
						{
							// Calculate the direction vector of this particle. (Negative distance, normalized.)
							Vector2f direction = -Vector2f(difference.x / distance, difference.y / distance);
							// Add other particle's pressure and viscocity contributions using Navier-Stokes equations, scaled by Muller kernels.
							pressure_contribution += direction * (particle_mass) * (particle.pressure + other_particle.pressure) / (2 * other_particle.density) * SPIKY_GRAD * pow(particle_diameter - distance, 2.f);
							viscocity_contribution += particle_viscocity * particle_mass * (other_particle.velocity - particle.velocity) / other_particle.density * VISC_LAP * (particle_diameter - distance);
						}
					}
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
	// Update velocity and position of particles.
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


