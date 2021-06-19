#include <math.h>
#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <unistd.h>

//#include "OpenMP sources/omp.h"

#ifdef _OPENMP

#include <omp.h>

#else

double omp_get_wtime() {
    return 0.0;
}

int omp_get_num_procs() {
    return 1;
}


#endif// ifdef _OPENMP


void selectionSortOfPart(double*, const int, const int);
void mergeParts(double*, const int, const int, const int);


void generateTwoArrays(double* first, double* second, int length, unsigned* seed) {
    for (int j = 0; j < length; j++) {
        first[j] = rand_r(seed) % 900 + 1;
    }

    for (int j = 0; j < length / 2; j++) {
        second[j] = 900 + rand_r(seed) % (9 * 900 + 1);
    }
}


void map(double* first, double* second, int length) {

#pragma omp parallel for
    for (int j = 0; j < length; j++) {
        first[j] = pow(sinh(first[j]), 2);
    }

    for (int j = length / 2 - 1; j > 0; j--) {
        second[j] = sqrt((second[j] + second[j - 1]) * exp(1));
    }
}


void merge(double* first, double* second, int length) {

#pragma omp parallel for shared(first, second)
    for (int j = 0; j < length / 2; j++) {
        second[j] = second[j]/first[j];
    }
}

void selectionSort(double* array, const int length) {
    int const numberOfSortingThreads = omp_get_num_procs() - 1;
    int const partSize = length / numberOfSortingThreads;
    int const lastPartSize = length - (numberOfSortingThreads) * partSize;

#pragma omp parallel sections
    {
        for (int i = 0; i < numberOfSortingThreads; i++) {
#pragma omp section
            selectionSortOfPart(array, i * partSize, partSize);
        }

#pragma omp section
        selectionSortOfPart(array, numberOfSortingThreads * partSize, lastPartSize);
    }

    mergeParts(array, partSize, length, numberOfSortingThreads + 1);
}

void reduce(double* array, double* result, const int length) {

    double minNonZero = __DBL_MAX__;
    double res = 0;

#pragma omp parallel
    {
#pragma omp for reduction(min:minNonZero)
        for (int j = 0; j < length / 2; j++) {
            double value = array[j];
            if (minNonZero > value && value != 0) {
                minNonZero = value;
            }
        }

#pragma omp for reduction(+:res)
        for (int j = 0; j < length / 2; j++) {
            if ((int)floor(array[j] / minNonZero) % 2) {
                res += sin(array[j]);
            }
        }
    }
    *result = res;
}

void selectionSortOfPart(double* array, const int offset, const int length) {

    for (int i = offset; i < offset + length - 1; i++) {
        int minIndex = i;
        int localIndex = minIndex;

#pragma omp parallel shared(minIndex) private(localIndex)
        {
#pragma omp for
            for (int j = i + 1; j < offset + length; j++) {
                if (array[j] < array[localIndex]) {
                    localIndex = j;
                }
            }

#pragma omp critical
            if (array[localIndex] < array[minIndex]) {
                minIndex = localIndex;
            }
        }

        if (minIndex != i) {
            double const temp = array[i];
            array[i] = array[minIndex];
            array[minIndex] = temp;
        }
    }
}

void mergeParts(double* array, const int partSize, const int length, const int numberOfParts) {

    selectionSortOfPart(array, 0, length);
}

void timer(double* completion) {
    int second = 0;
    while (*completion < 100) {
        sleep(1);
        second++;
        printf("%d: %f\n", second, *completion);
    }
}
int main(int argc, const char * argv[]) {
    double completion = 0;

#pragma omp parallel sections shared(completion)
    {

#pragma omp section
        {
            timer(&completion);
        }

#pragma omp section
        {
            int const length = atoi(argv[1]);
            unsigned seed;
            double results[50] = { 0 };

            double const startTime = omp_get_wtime();

            for (unsigned i = 0; i < 50; i++) {
                srand(i);
                double firstArray[length];
                double secondArray[length / 2];
                seed = i;
                generateTwoArrays(firstArray, secondArray, length, &seed);
                map(firstArray, secondArray, length);
                merge(firstArray, secondArray, length);
                selectionSort(secondArray, length / 2);
                reduce(secondArray, results + i, length / 2);
                completion += 100.0 / 50;
            }

            printf("\nN=%d. Milliseconds passed: %f\n", length, 1000 * (omp_get_wtime() - startTime));
        }
    }
}
