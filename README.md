# PACO-VMP
This is the algorithm used in the paper "PACO-VMP: Parallel Ant Colony Optimization for Virtual Machine Placement". This code may freely be used for comparison with your own algorithm or other experimental purposes.

While the code should take few if any adjustments to get running, certain alterations may need to be made on a machine-by-machine basis.

For machines that are AVX2 compatible (which includes most CPUs developed in the last decade), \_\_AVX2__ should be defined in platform.h (this is done by default). Otherwise, ensure it is not defined.

## Code Compilation
It should be noted that the Boost library is required in order for the code to compile successfully, which should be included using the -I compiler flag. The -mavx2 flag should also be included if your machine is AVX2 compatible, as well as -fopenmp for parallelization and -O2 for optimization.

## Input Parameters
The algorithm has 6 input parameters, which should be included in the following order:

1. Problem file: Assuming that the specified file is contained in the "Instances" directory, this is simply the file name, i.e. VMP_C177.vmp
2. Number of Ants: This can be any integer value, although we recommend a number that is either the number of available threads on your hardware or a multiple of that number in order to get the most from the OpenMP parallelization.
3. Number of Iterations: This can be any integer value. In the experiments for the paper, we used 50.
4. Seed: A seed for the ranluxgen random number generator. We generally use 5 digit values.
5. Alpha: The ACO parameter that determines the influence of pheromone. This should be an integer value between 0 and 6, as higher numbers can lead to floating point underflows. Parameter tuning indicates that 1 is the best value for Alpha, at least for our instances.
6. Beta: The ACO parameter that determines the influence of heuristic. This should be an integer value between 0 and 6, as higher numbers can lead to floating point underflows. Parameter tuning indicates that 6 is the best value for Beta, at least for our instances.
7. Rho: The ACO parameter that determines the pheromone decay rate. This can be any floating point value between 0 and 1, with 1 leading to no pheromone decay and 0 fully decaying the pheromone after each iteration. Parameter tuning indicates that 0.8 is a good value for Rho, although the tuning was less conclusive for Rho than Alpha or Beta.
