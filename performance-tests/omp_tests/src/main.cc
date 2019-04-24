#include <iostream>
#include <omp.h>


int main(int argc, char **argv) {
    int sum = 50;

#pragma omp parallel default(shared)
    {
        #pragma omp single
        {
            std::cout << "Spawned threads: " << omp_get_num_threads() << std::endl;
        }
        #pragma omp for reduction(+:sum)
        for (int i = 0; i < 1000; i++) {
                sum++;
        }
    }

    std::cout << "Reduction: " << sum << std::endl;

    return 0;
}
