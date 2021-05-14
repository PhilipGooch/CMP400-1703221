#include "quad_tree.h"
#include "circle.h"
#include "particle.h"

QuadTree::QuadTree(Rectangle rectangle, int capacity) :
	rectangle_(rectangle),
	capacity_(capacity),
	divided_(false),
	top_left_(nullptr),
	top_right_(nullptr),
	bottom_left_(nullptr),
	bottom_right_(nullptr)
{
}

QuadTree::~QuadTree()
{
}

void QuadTree::insert(Particle* particle)
{
	// If particle is in the boundary of the Quad Tree node.
	if (rectangle_.contains(particle))
	{
		// If there is still space availaible in the node's array.
		if (particles_.size() < capacity_)
		{
			// Store the particle in this array.
			particles_.push_back(particle);
		}
		// If there is no space availaible in the node's array.
		else
		{
			// If node is not yet subdived.
			if (!divided_)
			{
				// Subdivide the node into 4 smaller Quad Trees.
				subdivide();
			}
			// Recursively call the insert function on each of the child nodes.
			top_left_->insert(particle);
			top_right_->insert(particle);
			bottom_left_->insert(particle);
			bottom_right_->insert(particle);
		}
	}
}

void QuadTree::subdivide()
{
	// Top left position, width and height of the Quad Tree node's boundary. 
	int x = rectangle_.position_.x;
	int y = rectangle_.position_.y;
	int w = rectangle_.dimensions_.x;
	int h = rectangle_.dimensions_.y;

	// 4 rectangles to describe the boundaries of the child Quad Tree nodes.
	Rectangle top_left    (x,         y,         w / 2, h / 2);
	Rectangle top_right   (x + w / 2, y,         w / 2, h / 2);
	Rectangle bottom_left (x,         y + h / 2, w / 2, h / 2);
	Rectangle bottom_right(x + w / 2, y + h / 2, w / 2, h / 2);

	// Creating the Quad Tree nodes.
	top_left_ = new QuadTree(top_left, capacity_);
	top_right_ = new QuadTree(top_right, capacity_);
	bottom_left_ = new QuadTree(bottom_left, capacity_);
	bottom_right_ = new QuadTree(bottom_right, capacity_);

	// Flagging current Quad Tree node as divided.
	divided_ = true;
}

std::vector<Particle*> QuadTree::query(Circle circle)
{
	// Particles inside circle.
	std::vector<Particle*> particles;

	// If the circle overlaps the Quad Tree node's boundary.
	if (rectangle_.intersects(circle))
	{
		// Loop over particles contained in this node.
		for (Particle* particle : particles_)
		{
			// If the particle is inside the circle.
			if (circle.contains(particle))
			{
				// Add the particle to the particles array to be returned.
				particles.push_back(particle);
			}
		}
		// If the node is subdivided.
		if (divided_)
		{
			// Recursively call query function on child nodes and store returned particles.
			std::vector<Particle*> top_left_particles = top_left_->query(circle);
			std::vector<Particle*> top_right_particles = top_right_->query(circle);
			std::vector<Particle*> bottom_left_particles = bottom_left_->query(circle);
			std::vector<Particle*> bottom_right_particles = bottom_right_->query(circle);

			// Add particles from child nodes to particles array to be returned.
			particles.insert(particles.end(), top_left_particles.begin(), top_left_particles.end());
			particles.insert(particles.end(), top_right_particles.begin(), top_right_particles.end());
			particles.insert(particles.end(), bottom_left_particles.begin(), bottom_left_particles.end());
			particles.insert(particles.end(), bottom_right_particles.begin(), bottom_right_particles.end());
		}
	}
	return particles;
}

std::vector<Particle*> QuadTree::query(Rectangle rectangle)
{
	// Particles inside rectangle.
	std::vector<Particle*> particles;

	// If the rectangle overlaps the Quad Tree node's boundary.
	if(rectangle_.intersects(rectangle))
	{
		// Loop over particles contained in this node.
		for (Particle* particle : particles_)
		{
			// If the particle is inside the rectangle.
			if (rectangle.contains(particle))
			{
				// Add the particle to the particles array to be returned.
				particles.push_back(particle);
			}
		}
		if (divided_)
		{
			// Recursively call query function on child nodes and store returned particles.
			std::vector<Particle*> top_left_particles = top_left_->query(rectangle);
			std::vector<Particle*> top_right_particles = top_right_->query(rectangle);
			std::vector<Particle*> bottom_left_particles = bottom_left_->query(rectangle);
			std::vector<Particle*> bottom_right_particles = bottom_right_->query(rectangle);

			// Add particles from child nodes to particles array to be returned.
			particles.insert(particles.end(), top_left_particles.begin(), top_left_particles.end());
			particles.insert(particles.end(), top_right_particles.begin(), top_right_particles.end());
			particles.insert(particles.end(), bottom_left_particles.begin(), bottom_left_particles.end());
			particles.insert(particles.end(), bottom_right_particles.begin(), bottom_right_particles.end());
		}
	}
	return particles;
}

void QuadTree::render(sf::RenderWindow* window)
{
	// Draw rectangle for boundary of the Quad Tree node.
	sf::RectangleShape rectangle;
	rectangle.setOutlineColor(sf::Color(50, 50, 50));
	rectangle.setOutlineThickness(1);
	rectangle.setFillColor(sf::Color::Transparent);
	rectangle.setSize(rectangle_.dimensions_);
	rectangle.setPosition(rectangle_.position_);
	window->draw(rectangle);

	// Recursively call render function on child nodes.
	if (divided_)
	{
		top_left_->render(window);
		top_right_->render(window);
		bottom_left_->render(window);
		bottom_right_->render(window);
	}
}

void QuadTree::release()
{
	// Recursively call release function on child nodes.
	if (divided_)
	{
		top_left_->release();
		top_right_->release();
		bottom_left_->release();
		bottom_right_->release();
	}

	// Delete child Quad Trees.
	delete top_left_;
	delete top_right_;
	delete bottom_left_;
	delete bottom_right_;
}

