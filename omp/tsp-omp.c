#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include "node.h"
#include "queue.h"

double calculate_first_LB(const double *min1, const double *min2, int cities) {
    double lower_bound;
    double total_min1 = 0, total_min2 = 0;

    for (int i = 0; i < cities; i++) {
        total_min1 = total_min1 + min1[i];
        total_min2 = total_min2 + min2[i];
    }

    lower_bound = 0.5 * (total_min1 + total_min2);

    return lower_bound;
}

double
calculate_new_LB(const double *min1, const double *min2, double lower_bound, int from, int to, double **weight_matrix) {
    double new_lower_bound;
    double cf, ct;

    if (weight_matrix[from][to] >= min2[from]) {
        cf = min2[from];
    } else {
        cf = min1[from];
    }

    if (weight_matrix[from][to] >= min2[to]) {
        ct = min2[to];
    } else {
        ct = min1[to];
    }

    new_lower_bound = lower_bound + weight_matrix[from][to] - ((cf + ct) / 2);

    return new_lower_bound;
}

void pop_first_node(int best_tour_cost, int *final_tour, int number_of_cities,
                    node *N, priority_queue_t *queue, double inf, double new_bound, int *tour2,
                    int cost, node *N2, int length, double **matrix_weights, double *min1, double *min2) {
    queue_push(queue, N);

    N = queue_pop(queue);

    // todo: rmv commented code, this will never happen for the first node
    /*if (N->lower_bound > best_tour_cost) {
        return final_tour;
    }*/

    for (int v = 1; v < number_of_cities; v++) {
        if (N->viewed[v] == 0 && v != N->currCity && matrix_weights[N->currCity][v] != inf) {
            new_bound = calculate_new_LB(min1, min2, N->lower_bound, N->currCity, v, matrix_weights);

            for (int i = 0; i < N->tourLength; i++) {
                tour2[i] = N->tour[i];
            }

            tour2[N->tourLength] = v;
            cost = N->cost + matrix_weights[N->currCity][v];

            length = N->tourLength + 1;

            N2 = create_node(tour2, length, new_bound, cost, v, number_of_cities);

            queue_push(queue, N2);
        }
    }
    free_node(N);
}

int *procedure(double **matrix_weights, int number_of_cities, double best_tour_cost, double *min1, double *min2,
               priority_queue_t *queue, int *final_tour, int limit_threads) {
    node *N, *N2, **VN;
    priority_queue_t *queue_master;
    int length, v, bestId, flag=0, stop, cnt;
    int tour2[number_of_cities], tour[number_of_cities + 1];
    double lower_bound, cost, new_bound=0, bestlb, inf = INFINITY;
    int check, finished=0, iter;

    bestlb=inf;

    queue_master = queue_create(compare_node);

    lower_bound = calculate_first_LB(min1, min2, number_of_cities);
    length = 0;
    cost = 0;
    tour[length] = 0;     // Initialize tour
    length = length + 1; 

    N = create_node(tour, length, lower_bound, cost, 0, number_of_cities);
    
    VN = (node**) malloc (sizeof(node*)*limit_threads);
    for(int i=0;i<limit_threads;i++){
        VN[i] = create_node(tour, length, lower_bound, cost, 0, number_of_cities);
    }

    pop_first_node(best_tour_cost,final_tour,number_of_cities,N,queue_master,inf,new_bound,tour2,cost,N2,length,matrix_weights,min1,min2);

    int qsize = queue_master->size/limit_threads;


#pragma omp parallel num_threads(limit_threads) default(none) private(cnt, queue, N, v, new_bound, tour2, length, cost, N2, check, iter) shared(finished, VN, limit_threads, flag, bestlb, bestId, qsize, queue_master, best_tour_cost, final_tour, min1, min2, matrix_weights, number_of_cities, inf)
    {

        queue = queue_create(compare_node);
        
        check = 0;
        cost = 0;
        length = 0;
        new_bound = 0;
        v=0;
        cnt=0;
        iter=1000;

        int tid = omp_get_thread_num();
#pragma omp critical
        {
            for(int i=0; i<qsize; i++){
                N = queue_pop(queue_master);
                queue_push(queue, N);
            }
        }
#pragma omp barrier

#pragma omp critical
        {
            if(queue_master->size!=0){
                N = queue_pop(queue_master);
                queue_push(queue, N);
            }
        }

        while (queue->size != 0) {
            
            N = queue_pop(queue);
            
            //SINCRONIZAÇÃO
            if(flag==limit_threads && finished == 0){
                #pragma omp critical
                {
                    if(N->lower_bound<bestlb){
                       
                        bestlb = N->lower_bound;
                        bestId = tid;

                    }
                }
                queue_push(queue,N);

                #pragma omp barrier

                if(bestId==tid){

                        for(int i=0; i<limit_threads;i++){
                            VN[i]=queue_pop(queue);
                        }
                }                
                #pragma omp barrier
                
                queue_push(queue, VN[tid]);
                N = queue_pop(queue);

                #pragma omp barrier
                flag=0;
                bestlb=inf;
                check=0;
            }

            if (N->lower_bound >= best_tour_cost) {
                #pragma omp atomic write
                    finished = 1;

                break;
            }
            if (N->tourLength == number_of_cities) {

#pragma omp critical
                {
                    if (N->cost + matrix_weights[N->currCity][0] < best_tour_cost && matrix_weights[N->currCity][0] != inf) {

                        best_tour_cost = N->cost + matrix_weights[N->currCity][0];
                        
                        for (int i = 0; i < (N->tourLength); i++) {
                            final_tour[i] = N->tour[i];
                        }
                        final_tour[N->tourLength] = 0;
                    }
                }
            }
            else {
                for (v = 1; v < number_of_cities; v++) {
                    cnt = cnt+1;

                    if (N->viewed[v] == 0 && v != N->currCity && matrix_weights[N->currCity][v] != inf) {
                        new_bound = calculate_new_LB(min1, min2, N->lower_bound, N->currCity, v, matrix_weights);

                        if (new_bound > best_tour_cost) {
                            continue;
                        }

                        for (int i = 0; i < N->tourLength; i++) {
                            tour2[i] = N->tour[i];
                        }

                        tour2[N->tourLength] = v;
                        cost = N->cost + matrix_weights[N->currCity][v];

                        length = N->tourLength + 1;

                        N2 = create_node(tour2, length, new_bound, cost, v, number_of_cities);

                        queue_push(queue, N2);
                    }
                }
            }
            if(cnt>iter && check == 0){
                
                #pragma omp atomic write
                flag=flag+1;
                check=1;
            }
            free_node(N);
        }
    }

    return final_tour;
}

// Calculate min1 and min2
void calculate_mins(double **matrix_weights, int number_of_cities, double *min1, double *min2) {
    int used_min_index;

    for (int i = 0; i < number_of_cities; i++) {
        for (int j = 0; j < number_of_cities; j++) {
            if (j == 0) {  // to allocate first found val as our min reference
                min1[i] = matrix_weights[i][j];
                used_min_index = j;
            } else {
                if (min1[i] > matrix_weights[i][j]) {
                    min1[i] = matrix_weights[i][j];
                    used_min_index = j;
                }
            }
        }
        for (int k = 0; k < number_of_cities; k++) {
            if (k == 0 && used_min_index != k) {
                min2[i] = matrix_weights[i][k];
            } else if (k == 0) {
                min2[i] = matrix_weights[i][k + 1];
                k++;
            } else {  // not first element (min reference)
                if (min2[i] > matrix_weights[i][k] && used_min_index != k && matrix_weights[i][k] != INFINITY) {
                    min2[i] = matrix_weights[i][k];
                }
            }
        }
    }
}

int main(int argc, char *argv[]) {
    FILE *fptr1;
    int i, j, number_of_cities, total_roads = 0;  // k, used_min_index;
    double max_value_accepted;
    double weight, inf = INFINITY;
    priority_queue_t *queue;
    double *min1, *min2;
    int *final_tour_cost;
    double final_cost = 0;
    int limit_threads = 4;  // tenho o numero de threads = ao numero de processadores que o meu computador possui!

    number_of_cities = 0;
    max_value_accepted = 0;

    if (argc != 3) {
        exit(1);
    }

    max_value_accepted = atoi(argv[2]);

    fptr1 = fopen(argv[1], "r");

    if (fptr1 == NULL) {
        printf("Error! opening 1st file");
        exit(0);
    }

    if (fscanf(fptr1, "%d %d", &number_of_cities, &total_roads) != 2)
        return 0;

    // allocate size to matrix and its elements
    double **matrix_weights = (double **)malloc(number_of_cities * sizeof(double *));
    for (int i = 0; i < number_of_cities; i++) {
        matrix_weights[i] = (double *)malloc(number_of_cities * sizeof(double));
    }

    // initialize matrix with infinities
    for (int i = 0; i < number_of_cities; i++) {
        for (int j = 0; j < number_of_cities; j++) {
            matrix_weights[i][j] = inf;
        }
    }

    // populate nodes' edge weights
    while (fscanf(fptr1, "%d %d %lf", &i, &j, &weight) != EOF) {
        matrix_weights[i][j] = weight;
        matrix_weights[j][i] = weight;
    }
    fclose(fptr1);

    // allocate mem
    min1 = malloc(number_of_cities * sizeof(double));
    min2 = malloc(number_of_cities * sizeof(double));
    final_tour_cost = calloc(number_of_cities + 1, sizeof(int)); 

//////////////////////////////////////////////////////////////////////////
    double exec_time;
    exec_time = -omp_get_wtime();

    calculate_mins(matrix_weights, number_of_cities, min1, min2);

    final_tour_cost = procedure(matrix_weights, number_of_cities, max_value_accepted, min1, min2, queue,
                                final_tour_cost, limit_threads);

    if (final_tour_cost[1] != 0) {
        for (i = 0; i < number_of_cities; i++) {
            final_cost = final_cost + matrix_weights[final_tour_cost[i]][final_tour_cost[i + 1]];
        }
    }
    if (max_value_accepted > final_cost && final_cost != 0) {
        fprintf(stdout, "%.1f\n", final_cost);
        for (i = 0; i < number_of_cities + 1; i++) {
            fprintf(stdout, "%d ", final_tour_cost[i]);
        }
        fprintf(stdout, "\n");
    } else {
        fprintf(stdout, "NO SOLUTION\n");
    }

    /////////////////////////////////////////////////////////////////////////////////
    exec_time += omp_get_wtime();
    fprintf(stderr, "%.1fs\n", exec_time);

    free(min1);
    free(min2);

    for (i = 0; i < number_of_cities; i++) {
        free(matrix_weights[i]);
    }

    free(matrix_weights);
    free(final_tour_cost);

    return 0;
}