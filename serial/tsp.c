#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "queue.h"
#include "node.h"
#include <math.h>
#include <omp.h>

double calculate_first_LB(double *min1, double *min2, int cities)
{

    double lower_bound;
    double total_min1 =0, total_min2 = 0;

    for(int i=0; i<cities; i++){
        total_min1 = total_min1 + min1[i];
        total_min2 = total_min2 + min2[i];
    }

    lower_bound = 0.5 * (total_min1 + total_min2);

    return lower_bound;
}

double calculate_new_LB(double *min1, double *min2, double lower_bound, int from, int to, double **weight_matrix)
{
    double new_lower_bound;
    double cf=0, ct=0;

    if (weight_matrix[from][to] >= min2[from]){
        cf = min2[from];
    }
    else{
        cf = min1[from];
    }

    if (weight_matrix[from][to] >= min2[to]){
        ct = min2[to];
    }
    else{
        ct = min1[to];
    }

    new_lower_bound = lower_bound + weight_matrix[from][to] - ((cf + ct) / 2);

    return new_lower_bound;
}

int* procedure(double **matrix_weights, int number_of_cities, double best_tour_cost, double *min1, double *min2, priority_queue_t *queue, int *final_tour){

    node *N, *N2;
    int length, v;
    int tour2[number_of_cities], tour[number_of_cities+1];
    double lower_bound, cost, new_bound, inf = INFINITY;

    lower_bound = calculate_first_LB(min1, min2, number_of_cities);
    length = 0;
    cost = 0;
    tour[length] = 0; //Initialize tour
    length = length +1;

    N = create_node(tour, length, lower_bound, cost, 0, number_of_cities);

    queue_push(queue, N);

    int check = 0;
    while(queue->size != 0){

        N = queue_pop(queue); 

        if(N->lower_bound >= best_tour_cost){
            return final_tour; 
        }
        if(N->tourLength == number_of_cities){
            if(N->cost + matrix_weights[N->currCity][0] < best_tour_cost && matrix_weights[N->currCity][0]!=inf){
        
                best_tour_cost = N->cost + matrix_weights[N->currCity][0];
                for(int i=0; i<(N->tourLength); i++){
                    final_tour[i] = N->tour[i];
                }
                final_tour[N->tourLength] = 0;
            }
        }
        else{   
            for(v=1; v<number_of_cities; v++){
                
                if(N->viewed[v] == 0 && v!=N->currCity && matrix_weights[N->currCity][v]!=inf){

                    if(check==1){

                        check=0;
                    }
                    new_bound = calculate_new_LB(min1, min2, N->lower_bound, N->currCity, v, matrix_weights);  

                    if(new_bound > best_tour_cost){
                        continue;
                    }

                    for(int i=0; i<N->tourLength; i++){
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
        free_node(N);
    }

    return final_tour;

}

int main(int argc, char *argv[])
{

    FILE *fptr1;
    int i, j, number_of_cities, total_roads=0, k, aux;
    double max_value_accepted;
    double weight, inf = INFINITY;
    priority_queue_t *queue;
    node *N4;
    double *min1, *min2;
    int *final_tour; 
    double final_cost=0;

    number_of_cities = 0;
    max_value_accepted = 0;

    if (argc != 3)
    {
        exit(1);
    }

    max_value_accepted = atoi(argv[2]);

    fptr1 = fopen(argv[1], "r");
    
    if (fptr1 == NULL)
    {
        printf("Error! opening 1st file");
        exit(0);
    }

    if(fscanf(fptr1, "%d %d", &number_of_cities, &total_roads) != 2)
        return 0;

    
    double **matrix_weights;

    matrix_weights = (double**)malloc(number_of_cities* sizeof(double*));

    for(int i=0; i<number_of_cities; i++){
        matrix_weights[i] = (double*)malloc(number_of_cities*sizeof(double));
    }

    //initialize matrix to inf
    for (int i = 0; i<number_of_cities; i++){
        for (int j=0; j<number_of_cities; j++){
            matrix_weights[i][j]=inf;
        }
    }

    while (fscanf(fptr1, "%d %d %lf", &i, &j, &weight) != EOF)
    {
        matrix_weights[i][j] = weight;
        matrix_weights[j][i] = weight;
    }
    
    fclose(fptr1);

    min1 = malloc(number_of_cities * sizeof (double));
    min2 = malloc(number_of_cities * sizeof (double));
    final_tour = calloc(number_of_cities+1,sizeof(int));

//////////////////////////////////////////////////////////////////////////
    double exec_time;
    exec_time = -omp_get_wtime();

    //Calculate min1 and min2
    for(i=0; i<number_of_cities; i++){
        for(j=0; j < number_of_cities; j++){
            if(j==0 ){
                min1[i] = matrix_weights[i][j];
                aux = j;
            }
            else{
                if(min1[i] > matrix_weights[i][j] ){
                    min1[i] = matrix_weights[i][j];
                    aux = j;
                }
            }
        }
        for(k=0; k<number_of_cities; k++){
            if(k==0 && aux!=k){
                min2[i] = matrix_weights[i][k];
            }
            else if(k==0){
                min2[i] = matrix_weights[i][k+1];
                k++;
            }
            else{
                if(min2[i] > matrix_weights[i][k] && aux !=k && matrix_weights[i][k]!=inf){
                    min2[i] = matrix_weights[i][k];
                }
            }
        }
    }

    queue = queue_create(compare_node);

    final_tour = procedure(matrix_weights, number_of_cities, max_value_accepted, min1, min2, queue, final_tour);

    if(final_tour[1] != 0){
        for(i=0; i<number_of_cities; i++){
            final_cost = final_cost + matrix_weights[final_tour[i]][final_tour[i+1]];
        }
    }
    if(max_value_accepted > final_cost && final_cost != 0){
        fprintf(stdout, "%.1f\n", final_cost);
        for(i=0; i<number_of_cities+1; i++){
            fprintf(stdout, "%d ", final_tour[i]);
        }
        fprintf(stdout, "\n");
    }    else{
        fprintf(stdout, "NO SOLUTION\n");
    }

/////////////////////////////////////////////////////////////////////////////////
    exec_time += omp_get_wtime();
    fprintf(stderr, "%.1fs\n", exec_time);

    while(queue->size != 0){
        N4 = queue_pop(queue);
        free_node(N4);
    }

    queue_delete(queue);
    free(min1);
    free(min2);

    for(i=0; i<number_of_cities; i++){
        free(matrix_weights[i]);
    }

    free(matrix_weights);
    free(final_tour);

    return 0;
}