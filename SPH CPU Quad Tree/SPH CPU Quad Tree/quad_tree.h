#pragma once
#include "rectangle.h"

#include <vector>
#include <SFML/Graphics.hpp>

class Particle;

class QuadTree
{
public:
	QuadTree(Rectangle rectangle, int capacity);
	~QuadTree();

	// Inserts a particle into the Quad Tree
	void insert(Particle* particle);

	// Subdivides the Quad Tree by initializing 4 smaller child quad trees. Splitting it into quarters.
	void subdivide();

	// Queries the Quad Tree. Returns the particles that fall within the given shape.
	std::vector<Particle*> query(Circle circle);
	std::vector<Particle*> query(Rectangle rectangle);

	// Renders The boundary of the Quad tree for visualization.
	void render(sf::RenderWindow* window);

	// Recursive release function. Deletes child Quad Trees.
	void release();

	Rectangle rectangle_;					// Rectangle describing boundary of the Quad Tree node.
	int capacity_;							// Number of particles that can be stored in each node.
	bool divided_;							// Flag for if the node is subdivided.
	std::vector<Particle*> particles_;		// Particles stored in this node.
	QuadTree* top_left_;					// Child Quad Tree nodes.
	QuadTree* top_right_;
	QuadTree* bottom_left_;
	QuadTree* bottom_right_;
};

