#pragma once
#include <SFML/System/Vector2.hpp>

class Particle;

class Circle
{
public:
	Circle();
	Circle(float x, float y, float radius);

	// Checks if a particle is inside this circle.
	bool contains(Particle* particle);

	// Checks if a circle intersects with this circle
	bool intersects(Circle circle);

	sf::Vector2f position_;			// Centre of the circle.
	float radius_;					// Radius of the circle.
};

