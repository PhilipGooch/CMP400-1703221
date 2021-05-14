#pragma once
#include <SFML/Graphics.hpp>

class Input
{
public:

	Input();

	// Getters and Setters
	///////////////////////////////////////////////////////////
	bool getKeyDown(int);
	void setKeyDown(int);
	void setKeyUp(int);
	inline sf::Vector2i getMousePosition() { return mouse_position; }
	inline bool getMouseLeftDown() { return mouse_left_down; }
	inline void setMousePosition(int x, int y) { mouse_position.x = x; mouse_position.y = y; }
	inline void setMouseLeftDown(bool mouseLeftDown) { mouse_left_down = mouseLeftDown; }

private:
	
	sf::Vector2i mouse_position;
	bool mouse_left_down;
	bool keys[256]{ false };		

};


