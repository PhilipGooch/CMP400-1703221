#include "Application.h"

#include "quad_tree.h"
#include "Particle.h"

#if BENCHMARKING
#include <chrono>  

using std::chrono::duration_cast;
using std::chrono::nanoseconds;
using std::chrono::steady_clock;
using namespace std::chrono;
#endif

Application::Application(RenderWindow* window) :
	window_(window),
	running_(true),
	environment_width_(1024),
	environment_height_(512),
	number_of_particles_(512),
	gravity_magnitude_(4200 * 9.8f),
	gravity_(0.f, gravity_magnitude_),
	damping_(-0.5f),
	integration_timestep_(0.0008f),
	particle_diameter_(8.0f),
	particle_radius_squared_(pow(particle_diameter_, 2)),
	particle_mass_(50.f),
	particle_viscocity_(250.f),
	REST_DENSITY(1000.f),
	BOLTZMANN_CONSTANT(2000.f),
	POLY6(315.f / (65.f * PI * pow(particle_diameter_, 9.f))),
	SPIKY_GRAD(-45.f / (PI * pow(particle_diameter_, 6.f))),
	VISC_LAP(45.f / (PI * pow(particle_diameter_, 6.f)))
#if BENCHMARKING
	,simulations_(10),
	simulation_(0),
	iterations_(1000),
	start_iteration_(25),
	iteration_(0),
	start_time_(0),
	time_sum_(0)
#endif
{
	srand(time(NULL));

	// Placeholder Quad Tree to save conditional check in update loop.
	quad_tree_ = new QuadTree(Rectangle(0, 0, 0, 0), 0);

	initializeParticles();

#if BENCHMARKING
	std::ofstream temp(to_string(number_of_particles_) + ".csv", std::ofstream::out);
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
	while (running_)
	{
		// Queery relevant Windows input messages.
		Event event;
		while (window_->pollEvent(event))
		{
			switch (event.type)
			{
			// If the X symbol in top right is clicked, exit application.
			case Event::Closed:
				running_ = false;
				break;
			// If a key is pressed, tell input manager what key has been pressed.
			case Event::KeyPressed:
				input_.setKeyDown(event.key.code);
				break;
			// If a key is released, tell input manager what key has been released.
			case Event::KeyReleased:
				input_.setKeyUp(event.key.code);
				break;
			// If a mouse button is pressed, tell input manager what mouse button has been pressed.
			case Event::MouseButtonPressed:
				if (event.mouseButton.button == sf::Mouse::Left)
				{
					input_.setMouseLeftDown(true);
				}
				break;
			// If a mouse button is released, tell input manager what mouse button has been released.
			case Event::MouseButtonReleased:
				if (event.mouseButton.button == sf::Mouse::Left)
				{
					input_.setMouseLeftDown(false);
				}
				break;
			// If the mouse is moved, tell input manager the new coordinates of the mouse.
			case Event::MouseMoved:
				input_.setMousePosition(event.mouseMove.x, event.mouseMove.y);
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
	window_->close();
}

void Application::handleInput()
{
	// Manipulating Gravity.
	if (input_.getKeyDown('a' - 97))
	{
		gravity_.x = -gravity_magnitude_;
		gravity_.y = 0;
	}
	if (input_.getKeyDown('d' - 97))
	{
		gravity_.x = gravity_magnitude_;
		gravity_.y = 0;
	}
	if (input_.getKeyDown('w' - 97))
	{
		gravity_.x = 0;
		gravity_.y = -gravity_magnitude_;
	}
	if (input_.getKeyDown('s' - 97))
	{
		gravity_.x = 0;
		gravity_.y = gravity_magnitude_;
	}
}

void Application::update()
{
	// Releases all memory from dynamically alocated child quad trees.
	quad_tree_->release();
	delete quad_tree_;

	// Constructing new Quad Tree every frame.
	quad_tree_ = new QuadTree(Rectangle(0, 0, max(environment_width_, environment_height_), max(environment_width_, environment_height_)), 4);
	for (Particle& particle : particles_)
	{
		quad_tree_->insert(&particle);
	}

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
		time_sum_ /= iterations_;
		output_file_stream_ << time_sum_ << ", ";
		long long end_time = system_clock::now().time_since_epoch().count();
		float seconds = (float)(end_time - start_time_) / 10000000;
		long long fps = (iterations_ / seconds);
		output_file_stream_ << "," << fps << ",\n";
		iteration_ = 0;
		simulation_++;
		time_sum_ = 0;
		reset();
		start_time_ = system_clock::now().time_since_epoch().count();
	}
#endif
}

void Application::render()
{
	window_->clear(Color::Black);

	// Renders Quad Tree visualisation.
	quad_tree_->render(window_);

	// Rendering particles.
	for (Particle& particle : particles_)
	{
		RectangleShape square;
		square.setFillColor(Color(0, 127, 255));
		square.setSize(Vector2f(4, 4));
		square.setPosition(particle.position_);
		square.setOrigin(Vector2f(2, 2));
		window_->draw(square);
	}

	window_->display();
}

void Application::initializeParticles()
{
	// Loop up from the bottom of the screen until all particles have been initialized.
 	for (float y = environment_height_; ; y -= particle_diameter_)
	{
		// Loop the width of the dam.
		for (float x = 0; x < 16 * particle_diameter_; x += particle_diameter_)
		{
			if (particles_.size() < number_of_particles_)
			{
				// A small offest applied to X axis on initialization of particles.
				float jitter = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);		
				particles_.push_back(Particle(x + jitter, y));
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
	if (particle.position_.x < 0.0f)
	{					 
		particle.velocity_.x *= damping_;
		particle.position_.x = jitter;
	}					 
	// RIGHT			 
	if (particle.position_.x > environment_width_ - 1)
	{					 
		particle.velocity_.x *= damping_;
		particle.position_.x = environment_width_ - 1 - jitter;
	}	
	// BOTTOM			 
	if (particle.position_.y < 0.0f)
	{
		particle.velocity_.y *= damping_;
		particle.position_.y = jitter;
	}
	// BOTTOM			 
	if (particle.position_.y > environment_height_ - 1)
	{					 
		particle.velocity_.y *= damping_;
		particle.position_.y = environment_height_ - 1 - jitter;
	}
}

void Application::computeDensityPressure()
{
	
	for (Particle& particle : particles_)
	{
		// Querying the quad tree to get particles in a radius around the current particle.
		circle.position_.x = particle.position_.x;
		circle.position_.y = particle.position_.y;
		circle.radius_ = particle_diameter_;
		std::vector<Particle*> other_particles = quad_tree_->query(circle);

		// Reset particles density so it can be recalculated.
		particle.density_ = 0.f;
		for (Particle* other_particle : other_particles)
		{
			// Calculate vector between the pair of particles.
			Vector2f difference = other_particle->position_ - particle.position_;
			// Calculate the squared distance between the pair of particles.
			float distance_squared = pow(difference.x, 2) + pow(difference.y, 2);
			// If the particles are overlapping.
			if (distance_squared < particle_radius_squared_)
			{
				// Add the other particle's mass, scaled by a muller kernel to this particle's density.
				particle.density_ += particle_mass_ * POLY6 * pow(particle_radius_squared_ - distance_squared, 3.f);
			}
		}
		// Calculate pressure of particle using an "equation of state", relating it's density to a given rest density.
		particle.pressure_ = BOLTZMANN_CONSTANT * (particle.density_ - REST_DENSITY);
	}
}

void Application::computeForces()
{
	for (Particle& particle : particles_)
	{
		// Querying the quad tree to get particles in a radius around the current particle.
		circle.position_.x = particle.position_.x;
		circle.position_.y = particle.position_.y;
		circle.radius_ = particle_diameter_;
		std::vector<Particle*> other_particles = quad_tree_->query(circle);

		// Reset pressure and viscocity calcualtions so they can be recalculated.
		Vector2f pressure_contribution(0.f, 0.f);
		Vector2f viscocity_contribution(0.f, 0.f);
		for (Particle* other_particle : other_particles)
		{
			// Do not compare particles if they are the same.
			if (&particle != other_particle)
			{
				// Calculate vector between the pair of particles.
				Vector2f difference = other_particle->position_ - particle.position_;
				// Calculate the distance between the pair of particles. Distance needed later so no benefit from comparing squared distances.
				float distance = sqrt(pow(difference.x, 2) + pow(difference.y, 2));
				// If particles are overlapping.
				if (distance != 0 && distance < particle_diameter_)	// <--- distance != 0 is handling case where particles are at the same position.
				{
					// Calculate the desired direction vector of this particle.
					Vector2f desired = -Vector2f(difference.x / distance, difference.y / distance);
					// Add other particle's pressure and viscocity contributions using Navier-Stokes equations, scaled by Muller kernels.
					pressure_contribution += desired * particle_mass_ * (particle.pressure_ + other_particle->pressure_) / (2.f * other_particle->density_) * SPIKY_GRAD * pow(particle_diameter_ - distance, 2.f);
					viscocity_contribution += particle_viscocity_ * particle_mass_ * (other_particle->velocity_ - particle.velocity_) / other_particle->density_ * VISC_LAP * (particle_diameter_ - distance);
				}
			}
		}
		// Calculate grivity contributions by multiplying gravity by this particle's calculated density.
		Vector2f gravity_contribution = gravity_ * particle.density_;
		// Add all force contributions together to calculate the particle's force.
		particle.force_ = pressure_contribution + viscocity_contribution + gravity_contribution;
	}
}

void Application::integrate()
{
	// Update velocity and positions of particles.
	for (Particle& particle : particles_)
	{
		particle.velocity_ += particle.force_ / particle.density_ * integration_timestep_;
		particle.position_ += particle.velocity_ * integration_timestep_;
		// Handle collision detection at edges of environment rectangle.
		edges(particle);
	}
}

void Application::reset()
{
	particles_.clear();
	initializeParticles();
}
