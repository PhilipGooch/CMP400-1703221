#pragma once
#include <SFML/Graphics.hpp>

class Particle
{
public:
	Particle(float _x, float _y) : position_(_x, _y), velocity_(0.f, 0.f), force_(0.f, 0.f), density_(0), pressure_(0.f) {}
	sf::Vector2f position_, velocity_, force_;
	float density_, pressure_;
};

