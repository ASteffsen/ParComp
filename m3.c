#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <omp.h>
void map(long double* M1, long double* M2, int N, int TN){
    /*to each M1 use hyperbolic sinus, then ^2*/
    int N1 = N/2;
    /*each element in M2 add to previous and then *e then sqrt*/
    long double cM2[N1];
#pragma omp parallel num_threads(TN)
    {
#pragma omp for schedule(guided,5)
        for (int j = 0; j < N; j++) {
            M1[j] = pow(sinh(M1[j]), 2);
        }
        for (int j = 0; j < N1; j++) {
            cM2[j] = M2[j];
        }
        M2[0] = sqrt(M2[0] * exp(1));
        for (int j = 1; j < N1; j++) {
            M2[j] = sqrt((M2[j] + cM2[j - 1]) * exp(1));
        }
    }
}
void merge(long double* M1, long double* M2, int N, int TN){

    /*to each element with similar index M2[i]=M1[i]/M2[i]*/
#pragma omp parallel for schedule(guided,5) num_threads(TN) shared(M1, M2)
    for (int j=0; j<N; j++){
        M2[j] = M1[j]/M2[j];
    }
}
//void sort(long double* M2, int N){

    /*selection sort*/
   /* for (int i = 0; i < N - 1; i++) {
        int minIndex= i;
        int localIndex = minIndex;
        

#pragma omp parallel shared(minIndex) private(localIndex)
        {
#pragma omp for
	
            for (int j = i + 1; j < N; j++)
                if (M2[j] < M2[localIndex]) localIndex = j;

#pragma omp critical
            if (localIndex < minIndex) minIndex = localIndex;
        }  

        if (minIndex != i) {
            long double const temp = M2[i];
            M2[i] = M2[minIndex];
            M2[minIndex] = temp;
        }
    }

}*/

struct Compare { long double val; int index; };
#pragma omp declare reduction(maximum : struct Compare : omp_out = omp_in.val > omp_out.val ? omp_in : omp_out)

void sort(long double* arr, int size, int TN)
{
    for (int i = size - 1; i > 0; --i)
    {
        struct Compare max;
        max.val = arr[i];
        max.index = i;
        #pragma omp parallel for schedule(guided,5) num_threads(TN) reduction(maximum:max)
        for (int j = i - 1; j >= 0; --j)
        {
            if (arr[j] > max.val)
            {
                max.val = arr[j];
                max.index = j;
            }
        }
        int tmp = arr[i];
        arr[i] = max.val;
        arr[max.index] = tmp;
    }
}

void reduce(long double* M2, int N, int TN){

    double minNonZero = __DBL_MAX__;
    double res = 0;

#pragma omp parallel num_threads(TN)
    {
#pragma omp for schedule(guided,5) reduction(min:minNonZero)
        for (int j = 0; j < N / 2; j++) {
            double value = M2[j];
            if (minNonZero > value && value != 0) {
                minNonZero = value;
            }
        }

#pragma omp for schedule(guided,5) reduction(+:res)
        for (int j = 0; j < N / 2; j++) {
            if ((int)floor(M2[j] / minNonZero) % 2) {
                res += sin(M2[j]);
            }
        }
    }
    printf(" X= %f\n", res);
}

int main(int argc, char* argv[]) {
    printf ("0");
    int N;
    struct timeval T1,T2;
    unsigned seed;
    long delta_ms;
    N = atoi(argv[1]);
    int TN = atoi(argv[2]);
    gettimeofday(&T1, NULL);
    for (unsigned i=0; i<50; i++)
    {
        srand(i);
        seed = i;

        /*stage generate*/
        /*Create array M1 of N elements*/
        long double M1[N];
        int j;
        int A = 900;
        for (j = 0; j<N; j++){
            M1[j] = rand_r(&seed) % A;

        }
        //printf(" M1 - %f", M1[1]);
        /*Fill M1 with rand_r from 1 to A*/
        /*Create array M2 of N/2*/
        int N1 = N/2;
        long double M2[N1];

        /*Fill M2 with rand_r from A to 10*A */
        for (j=0; j<N1; j++){
            M2[j] = A + rand_r(&seed) % 10*A;
        }
        //printf(" M2 - %f", M2[1]);
        /*stage map*/
        map(M1, M2, N, TN);
        /*stage merge*/
        merge(M1,M2, N1, TN);
        /*stage sort*/
        sort(M2, N1, TN);
        /*stage reduce*/
        reduce(M2, N1, TN);


    }
    gettimeofday(&T2,NULL);
    delta_ms = 1000*(T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec)/1000;
    printf("\nN=%d. Milliseconds passed: %ld\n",N, delta_ms);

    return 0;
} 
