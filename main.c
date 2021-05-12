#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <ipp.h>
#include <ippvm.h>


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
        ippsSinh_64f_A50(M1, M1, N);
        ippsPowx_64f_A50(M1, 2, M1, N);
        /*each element in M2 add to previous and then *e then sqrt*/
        double M2_CP[N1];
        M2_CP[0] = 0;
        ippsCopy_64f(M2, &M2_CP[1], N1);
        ippsAdd_64f(M2, M2_CP, M2, N1);
        double val = exp(1);
        ippsMulC_64f(M2, val, M2, N1);
        /*stage merge*/
        /*to each element with similar index M2[i]=M1[i]/M2[i]*/
        ippsDiv_64f(M1, M2, M2, N1);
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
        long double X=0;
        int b;
        long double a;
        for (j=0; j<N1; j++){
            a = M2[j]/minM2;
            b = (int) a;
            if (b % 2 == 0){
                X = X+M2[j];
            }
        }
        printf(" X= %Lf\n", X);

    }
    gettimeofday(&T2,NULL);
    delta_ms = 1000*(T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec)/1000;
    printf("\nN=%d. Milliseconds passed: %ld\n",N, delta_ms);

    return 0;
}
