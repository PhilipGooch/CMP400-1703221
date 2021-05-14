
// This function uses CUDAs "thrust/sort.h" and "thrust/device_ptr.h" headers.
// It takes a long time to compile so it is kept seperate.

// Sorts device particle variables based on particle hash values.
void c_sortParticles(int* particle_hashes, int* particle_IDs, int number_of_particles, float* positions, float* velocities, float* forces, float* densities, float* pressures);