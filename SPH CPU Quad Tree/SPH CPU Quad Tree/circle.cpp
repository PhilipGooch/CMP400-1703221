#include "circle.h"
#include "particle.h"
#include <math.h>

Circle::Circle()
{
}

Circle::Circle(float x, float y, float radius)
{
	position_.x = x;
	position_.y = y;
	radius_ = radius;
}

bool Circle::contains(Particle* particle)
{
	// Point and circle collision.
	return (pow((particle->position_.x - position_.x), 2) + pow((particle->position_.y - position_.y), 2) <= pow(radius_, 2));
}

bool Circle::intersects(Circle circle)
{
	// Circle and circle collision.
	return (pow((circle.position_.x - position_.x), 2) + pow((circle.position_.y - position_.y), 2) <= pow(radius_ + circle.radius_, 2));
}
