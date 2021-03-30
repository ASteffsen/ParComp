#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

int main(int argc, char* argv[]) {
    printf ("0");
    int N;
    struct timeval T1,T2;
    unsigned seed;
    long delta_ms;
    N = atoi(argv[1]);
    gettimeofday(&T1, NULL);
    for (unsigned i=0; i<50; i++)
    {
        srand(i);
        seed = i;
        /*stage generate*/
        /*Create array M1 of N elements*/
        double M1[N];
        int j,k;
        int A = 900;
        for (j = 0; j<N; j++){
            M1[j] = rand_r(&seed) % A;

        }
        printf(" M1 - %f", M1[1]);
        /*Fill M1 with rand_r from 1 to A*/
        /*Create array M2 of N/2*/
        int N1 = N/2;
        double M2[N1];

        /*Fill M2 with rand_r from A to 10*A */
        for (j=0; j<N1; j++){
            M2[j] = A + rand_r(&seed) % 10*A;
        }

        printf(" M2 - %f", M2[1]);
        /*stage map*/
        /*to each M1 use hyperbolic sinus, then ^2*/
        for (j=0; j<N; j++){
            M1[j] = pow(sinh(M1[j]), 2);
        }
        /*each element in M2 add to previous and then *e then sqrt*/
        double cM2[N1];
        for (j=0; j<N1; j++){
            cM2[j] = M2[j];
        }
        M2[0] = sqrt(M2[0] * exp(1));
        for (j=1; j<N1; j++){
            M2[j] = sqrt((M2[j]+cM2[j-1])*exp(1));
        }
        /*stage merge*/
        /*to each element with similar index M2[i]=M1[i]/M2[i]*/
        for (j=0; j<N1; j++){
            M2[j] = M1[j]/M2[j];
        }
        /*stage sort*/
        /*selection sort*/
        int min;
        for (j=0; j<N1; j++){
            min = i;

            for (k=j+1; k<N1; k++){
                if (M2[min]>M2[k]){
                    min = k;
                }
            }
            double t = M2[j];
            M2[j] = M2[min];
            M2[min] = t;
        }
        double minM2;
        /*cycle to find min not null*/
        k = 0;
        while (M2[k]==0){
            k++;
        }
        minM2 = M2[k];
        /*stage reduce*/
        /*count sum sin elements M2 that give odd number after delenie at min elements of M2*/
        double X=0;
        int b;
        double a;
        for (j=0; j<N1; j++){
            a = M2[j]/minM2;
            b = (int) a;
            if (b % 2 == 0){
                X = X+M2[j];
            }
        }
        printf(" X= %f\n", X);

    }
    gettimeofday(&T2,NULL);
    delta_ms = 1000*(T2.tv_sec - T1.tv_sec) + (T2.tv_sec - T1.tv_sec)/1000;
    printf("\nN=%d. Milliseconds passed: %ld\n",N, delta_ms);

    return 0;
}


