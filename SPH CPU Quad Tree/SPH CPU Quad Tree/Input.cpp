#include "Input.h"

Input::Input() :
	mouse_position(sf::Vector2i(0, 0)),
	mouse_left_down(false)
{
}

bool Input::getKeyDown(int key)
{
	return keys[key];
}

void Input::setKeyDown(int key)
{
	if (key >= 0)
	{
		keys[key] = true;
	}
}

void Input::setKeyUp(int key)
{
	if (key >= 0)
	{
		keys[key] = false;
	}
}

