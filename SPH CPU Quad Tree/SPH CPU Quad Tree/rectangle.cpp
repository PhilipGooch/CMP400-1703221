#include "rectangle.h"
#include "particle.h"
#include <math.h>

Rectangle::Rectangle()
{
}

Rectangle::Rectangle(float x, float y, float width, float height)
{
	position_.x = x;
	position_.y = y;
	dimensions_.x = width;
	dimensions_.y = height;
}

bool Rectangle::contains(Particle* particle)
{
	// Point and rectangle collision.
	return (particle->position_.x >= position_.x &&
		    particle->position_.x < position_.x + dimensions_.x &&
		    particle->position_.y >= position_.y &&
		    particle->position_.y < position_.y + dimensions_.y);
}

bool Rectangle::intersects(Rectangle rectangle)
{
	// AABB collision.
	return !(rectangle.position_.x >= position_.x + dimensions_.x ||
		     rectangle.position_.x + rectangle.dimensions_.x < position_.x ||
		     rectangle.position_.y >= position_.y + dimensions_.y ||
		     rectangle.position_.y + rectangle.dimensions_.y < position_.y);
}

// https://www.youtube.com/watch?v=62-pRVZuS5c
bool Rectangle::intersects(Circle circle)
{
	// Half width and height of boundary rectangle.
	sf::Vector2f extent = sf::Vector2f(dimensions_.x / 2, dimensions_.y / 2);

	// The point is the circle's position, relative to the centre of the rectangle.
	sf::Vector2f point = circle.position_ - (position_ + extent);

	// Signed Distance Field (SDF) function for a point and a rectangle. 
	// Distance squared from the centre of the circle to the edge of the rectangle.
	float distance_squared = pow(fmax(fabsf(point.x) - extent.x, 0), 2) + pow(fmax(fabsf(point.y) - extent.y, 0), 2);

	return distance_squared <= circle.radius_ * circle.radius_;
}