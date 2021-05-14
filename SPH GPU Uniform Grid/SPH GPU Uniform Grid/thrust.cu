#include "thrust/sort.h"            // sort_by_key()
#include "thrust/device_ptr.h"      // device_ptr<>

void c_sortParticles(int* particle_hashes, int* particle_IDs, int number_of_particles, float* positions, float* velocities, float* forces, float* densities, float* pressures)
{
	// Store two copies of particle hashes. 
	int* save_particle_hashes_A;				// <--- Gets sorted each time.
	int* save_particle_hashes_B;				// <--- Stays constant.
	cudaMalloc((void**)&save_particle_hashes_A, number_of_particles * sizeof(int));
	cudaMalloc((void**)&save_particle_hashes_B, number_of_particles * sizeof(int));
	cudaMemcpy(save_particle_hashes_A, particle_hashes, number_of_particles * sizeof(int), cudaMemcpyDeviceToDevice);
	cudaMemcpy(save_particle_hashes_B, particle_hashes, number_of_particles * sizeof(int), cudaMemcpyDeviceToDevice);

	// Original particlal hashes and particle IDs get sorted in order of low to high particle hashes.
	thrust::sort_by_key(thrust::device_ptr<int>(particle_hashes),
						thrust::device_ptr<int>(particle_hashes + number_of_particles),
						thrust::device_ptr<int>(particle_IDs));

	// positions get sorted with save_particle_hashes_A.
	thrust::sort_by_key(thrust::device_ptr<int>(save_particle_hashes_A),
						thrust::device_ptr<int>(save_particle_hashes_A + number_of_particles),
						thrust::device_ptr<float2>((float2*)positions));

	// save_particle_hashes_A gets reset with save_particle_hashes_B. (this repeats for each particle variable)
	cudaMemcpy(save_particle_hashes_A, save_particle_hashes_B, number_of_particles * sizeof(int), cudaMemcpyDeviceToDevice);

	// velocities.
	thrust::sort_by_key(thrust::device_ptr<int>(save_particle_hashes_A),
						thrust::device_ptr<int>(save_particle_hashes_A + number_of_particles),
						thrust::device_ptr<float2>((float2*)velocities));
	cudaMemcpy(save_particle_hashes_A, save_particle_hashes_B, number_of_particles * sizeof(int), cudaMemcpyDeviceToDevice);

	// forces.
	thrust::sort_by_key(thrust::device_ptr<int>(save_particle_hashes_A),
						thrust::device_ptr<int>(save_particle_hashes_A + number_of_particles),
						thrust::device_ptr<float2>((float2*)forces));
	cudaMemcpy(save_particle_hashes_A, save_particle_hashes_B, number_of_particles * sizeof(int), cudaMemcpyDeviceToDevice);

	// densities.
	thrust::sort_by_key(thrust::device_ptr<int>(save_particle_hashes_A),
						thrust::device_ptr<int>(save_particle_hashes_A + number_of_particles),
						thrust::device_ptr<float>(densities));
	cudaMemcpy(save_particle_hashes_A, save_particle_hashes_B, number_of_particles * sizeof(int), cudaMemcpyDeviceToDevice);

	// pressures.
	thrust::sort_by_key(thrust::device_ptr<int>(save_particle_hashes_A),
						thrust::device_ptr<int>(save_particle_hashes_A + number_of_particles),
						thrust::device_ptr<float>(pressures));

	// Release memory.
	cudaFree(save_particle_hashes_A);
	cudaFree(save_particle_hashes_B);
}