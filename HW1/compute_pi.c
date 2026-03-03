//
// Value of pi up to 100 decimal places (wikipedia.org):
// 3.14159 26535 89793 23846 26433 83279 50288 41971 69399 37510
//   58209 74944 59230 78164 06286 20899 86280 34825 34211 70679
//
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define MAX_THREADS 8192   // Change this for HW assignment

pthread_t p_threads[MAX_THREADS];
pthread_attr_t attr;

long int sample_points_per_thread;
unsigned int seeds[MAX_THREADS];   // separate seeds from hits (important!)

void *compute_pi(void *arg) {
    unsigned int seed = *(unsigned int *)arg;
    long int local_hits = 0;

    for (long int i = 0; i < sample_points_per_thread; i++) {
        double x = (double)rand_r(&seed) / (double)RAND_MAX;
        double y = (double)rand_r(&seed) / (double)RAND_MAX;
        double d = (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5);
        if (d < 0.25) local_hits++;

        // avoid seed *= (i+1) overflow weirdness; keep seed evolving safely
        seed = seed * 1103515245u + 12345u;
    }

    long int *ret = (long int *)malloc(sizeof(long int));
    if (ret == NULL) pthread_exit(NULL);
    *ret = local_hits;
    pthread_exit((void *)ret);
}

int main(int argc, char *argv[]) {
    struct timeval start, stop;
    double computed_pi, error_pi, total_time;

    long int total_hits = 0;
    long int sample_points;
    int num_threads;
    int status;

    if (argc != 3) {
        printf("Need two integers as input \n");
        printf("Use: <executable_name> <sample_points> <num_threads>\n");
        exit(0);
    }

    sample_points = atol(argv[argc - 2]);
    if ((num_threads = atoi(argv[argc - 1])) > MAX_THREADS) {
        printf("Maximum number of threads allowed: %d.\n", MAX_THREADS);
        exit(0);
    }

    sample_points_per_thread = sample_points / num_threads;

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    gettimeofday(&start, NULL);

    // create threads with independent seeds
    for (int i = 0; i < num_threads; i++) {
        seeds[i] = 123456u + (unsigned int)i * 17u;  // simple distinct seeds
        status = pthread_create(&p_threads[i], &attr, compute_pi, (void *)&seeds[i]);
        if (status != 0) printf("Non-zero status when creating thread # %d\n", i);
    }

    // join and sum hits
    for (int i = 0; i < num_threads; i++) {
        void *retptr = NULL;
        pthread_join(p_threads[i], &retptr);
        if (retptr != NULL) {
            total_hits += *(long int *)retptr;
            free(retptr);
        }
    }

    gettimeofday(&stop, NULL);
    total_time = (stop.tv_sec - start.tv_sec) + 0.000001 * (stop.tv_usec - start.tv_usec);

    computed_pi = (4.0 * (double)total_hits) / (double)sample_points;
    error_pi = fabs(3.14159265358979323846 - computed_pi) / 3.14159265358979323846;

    printf("Trials = %ld, Threads = %4d, pi = %14.10f, error = %8.2e, time (sec) = %8.4f\n",
           sample_points, num_threads, computed_pi, error_pi, total_time);

    pthread_attr_destroy(&attr);
    return 0;
}