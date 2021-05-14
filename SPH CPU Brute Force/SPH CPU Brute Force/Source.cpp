#include "Application.h"

using namespace sf;

int main()
{
    RenderWindow window(VideoMode(1024, 512), "CPU Brute Force. 256 particles.");

    Application application(&window);

    application.run();

    return 0;
}