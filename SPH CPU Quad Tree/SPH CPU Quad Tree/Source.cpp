#include "Application.h"

using namespace sf;

int main()
{
    RenderWindow window(VideoMode(1024, 512), "CPU Quad Tree. 512 particles.");

    Application application(&window);

    application.run();

    return 0;
}