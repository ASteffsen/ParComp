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



int const aTask = 900;
#define NUMBER_OF_ITERATIONS 50


void generateTwoArrays(double*, double*, int, unsigned int*);
void map(double*, double*, int);
void merge(double*, double*, int);
void selectionSort(double*, const int);
void reduce(double*, double*, const int);
void selectionSortOfPart(double*, const int, const int);
//void mergeParts(double*, const int, const int);
void mergeParts(double*, const int, const int, const int);
void timer(double*);

//MARK:- main()
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

                //MARK:- Generate
                seed = i;
                generateTwoArrays(firstArray, secondArray, length, &seed);


                //MARK:- Map
                map(firstArray, secondArray, length);


                //MARK:- Merge
                merge(firstArray, secondArray, length);


                //MARK:- Sort
                selectionSort(secondArray, length / 2);


                //MARK:- Reduce
                reduce(secondArray, results + i, length / 2);

                completion += 100.0 / NUMBER_OF_ITERATIONS;
            }

            printf("\nN=%d. Milliseconds passed: %f\n", length, 1000 * (omp_get_wtime() - startTime));
        }
    }
}

/// This function fills two arrays using generation task

void generateTwoArrays(double* first, double* second, int length, unsigned* seed) {
    //force single-threaded for correct result
    for (int j = 0; j < length; j++) {
        first[j] = rand_r(seed) % aTask + 1;
    }

    for (int j = 0; j < length / 2; j++) {
        second[j] = aTask + rand_r(seed) % (9 * aTask + 1);
    }
}

/// This function maps some functions over the array

void map(double* first, double* second, int length) {

#pragma omp parallel for
    for (int j = 0; j < length; j++) {
        first[j] = pow(sinh(first[j]), 2);
    }

    for (int j = length / 2 - 1; j > 0; j--) {
        second[j] = sqrt((second[j] + second[j - 1]) * exp(1));
    }
}

/// This function merges two arrays

void merge(double* first, double* second, int length) {

#pragma omp parallel for shared(first, second)
    for (int j = 0; j < length / 2; j++) {
        second[j] = second[j]/first[j];
    }
}


/// This function sorts given array using selection sorting algorithm
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

/// This function reduces array to result value
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


/// This function sorts given array using selection sorting algorithm
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


/// This function merges parts of given array in ascending order
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
