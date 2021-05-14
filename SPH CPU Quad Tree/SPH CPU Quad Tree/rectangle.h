#pragma once
#include "circle.h"
#include <SFML/System/Vector2.hpp>

class Particle;

class Rectangle
{
public:
	Rectangle();
	Rectangle(float x, float y, float width, float height);

	// Checks if a particle is inside this rectangle.
	bool contains(Particle* particle);

	// Checks if a rectangle intersects this rectangle.
	bool intersects(Rectangle rectangle);

	// Checks if a circle intersects with this rectangle
	bool intersects(Circle circle);

	sf::Vector2f position_;			// Top left corner of rectangle.
	sf::Vector2f dimensions_;		// Dimensions of rectangle.
};

