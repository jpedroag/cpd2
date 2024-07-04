#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
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

node** pop_first_node(int best_tour_cost, int *final_tour, int number_of_cities,
                    node *N, priority_queue_t *queue, double inf, double new_bound, int *tour2,
                    int cost, node *N2, int length, double **matrix_weights, double *min1, double *min2, int neighbours_of_zero) {

    node **vector;
    int k=0;
    int tour[2];

    tour[0] = 0;

    vector = (node **) malloc(sizeof(node *) * neighbours_of_zero);
    for (int i = 0; i < neighbours_of_zero; i++) {
        vector[i] = create_node(tour, length, 0, cost, 0, number_of_cities);
    }

    queue_push(queue, N);

    N = queue_pop(queue);

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

            vector[k] = N2;
            k++;
        }
    }
    free_node(N);

    return vector;
}

char *encode(char *buffer, node *N){

    int pos=0;

    memcpy(buffer + pos, &(N->tourLength), sizeof(int));
    pos += sizeof(int);
    memcpy(buffer + pos, &(N->tour), sizeof(int*) );
    pos += sizeof(int*);
    memcpy(buffer + pos, N->tour, sizeof(int)*N->tourLength);
    pos += sizeof(int)*N->tourLength;
    memcpy(buffer + pos, &(N->lower_bound), sizeof(double) );
    pos += sizeof(double);
    memcpy(buffer + pos, &(N->cost), sizeof(double) );
    pos += sizeof(double);
    memcpy(buffer + pos, &(N->currCity), sizeof(int) );
    pos += sizeof(int);
    memcpy(buffer + pos, &(N->n_cities), sizeof(int));
    pos += sizeof(int);
    memcpy(buffer + pos, &(N->viewed), sizeof(int*) );
    pos += sizeof(int*);
    memcpy(buffer + pos, N->viewed, sizeof(int)*N->n_cities);

    return buffer;
}

node *decode(char *buffer){

    int pos=0;
    node *N = (node*)malloc(sizeof(node));

    N->tourLength = *(int*)(buffer + pos);
    pos += sizeof(int);
    memcpy(&(N->tour), buffer + pos, sizeof(int*));
    pos += sizeof(int*);
    N->tour = (int*)malloc(sizeof(int) * N->tourLength);
    memcpy(N->tour, buffer + pos, sizeof(int) * N->tourLength);
    pos += sizeof(int) * N->tourLength;
    N->lower_bound = *(double*)(buffer + pos);
    pos += sizeof(double);
    N->cost = *(double*)(buffer + pos);
    pos += sizeof(double);
    N->currCity = *(int*)(buffer + pos);
    pos += sizeof(int);
    N->n_cities = *(int*)(buffer + pos);
    pos += sizeof(int);
    memcpy(&(N->viewed), buffer + pos, sizeof(int*));
    pos += sizeof(int*);
    N->viewed = (int*)malloc(sizeof(int) * N->n_cities);
    memcpy(N->viewed, buffer + pos, sizeof(int) * N->n_cities);

    return N;
}

int *procedure(double **matrix_weights, int number_of_cities, double best_tour_cost, double *min1, double *min2,
               priority_queue_t *queue, int *final_tour, int size, int rank, int neighbours_of_zero) {

    node *N, *N2, **VN, **vector_neighbours_zero;
    priority_queue_t *queue_master;
    int length, v, bestId, cnt;
    int tour2[number_of_cities], tour[number_of_cities + 1];
    double lower_bound, cost, new_bound = 0, bestlb, inf = INFINITY;
    double finished = 0.0;
    int neighbours_left, iter;
    char *buffer;
    double bestvec[size-1], bestcostvec[size-1], best_cost[size-1];
    int buffer_size;
    int *final_vec;
    int lowint = -1;
    double lowfloat = -1.0;

    bestlb = inf;

    lower_bound = calculate_first_LB(min1, min2, number_of_cities);
    length = 0;
    cost = 0;
    tour[length] = 0;     // Initialize tour
    length = length + 1;

    int qsize = qsize = neighbours_of_zero / size;

    N2 = create_node(tour, length, lower_bound, cost, 0, number_of_cities);

    queue_master = queue_create(compare_node);
    N = create_node(tour, length, lower_bound, cost, 0, number_of_cities);
    vector_neighbours_zero = pop_first_node(best_tour_cost, final_tour, number_of_cities, N, queue_master, inf, new_bound, tour2, cost, N2, length, matrix_weights, min1, min2, neighbours_of_zero);

    queue = queue_create(compare_node);

    VN = (node **) malloc(sizeof(node *) * qsize);
    for (int i = 0; i < qsize; i++) {
        VN[i] = create_node(tour, length, lower_bound, cost, 0, number_of_cities);
    }

    int root = 0;
    cost = 0;
    length = 0;
    new_bound = 0;
    v = 0;
    cnt = 0;
    iter = 1000;

    for(int i=rank*qsize; i<rank*qsize + qsize;i++){
        queue_push(queue, vector_neighbours_zero[i]);
    }

    neighbours_left = neighbours_of_zero - (qsize*size);

    MPI_Barrier(MPI_COMM_WORLD);

    int i =1;
    if(rank<neighbours_left && neighbours_left != 0){

        queue_push(queue, vector_neighbours_zero[(neighbours_of_zero - neighbours_left)+ i*rank]);
    }

    int maxIter = 0;
    if(size < 5) {
        maxIter = 1000000;
    } else {
        maxIter = 10000000;
    }

    while (queue->size != 0) {

        N = queue_pop(queue);

        //SINCRONIZAÇÃO
        if (cnt > iter && finished == 0.0) {

            printf("Rank: %d SYNC!\n", rank);

            cnt=0;
            if(iter<maxIter)
                iter=iter*10;

            bestId=0;
            bestlb=N->lower_bound;

            queue_push(queue, N);

            if(rank!=0){
                if(finished==1)
                    bestlb=inf;  ///////////////////////////////////////////////////
                
                MPI_Send(&bestlb, 1, MPI_DOUBLE, root, 0, MPI_COMM_WORLD);
                MPI_Send(&best_tour_cost, 1, MPI_DOUBLE, root, 1, MPI_COMM_WORLD);
            }
            ////0 calcula o melhor LB
            else{
                for(int i=1; i<size; i++){
                    MPI_Recv(&bestvec[i-1], 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(&best_cost[i-1], 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                for(int i=1; i<size; i++){
                    if(bestvec[i-1]<bestlb){
                        bestId=i;
                        bestlb=bestvec[i-1];
                    }
                    if(bestvec[i-1]==-1.0){
                        finished=-1.0;
                        bestId=-1;
                    }
                    if(best_cost[i-1] < best_tour_cost){
                        best_tour_cost = best_cost[i-1];
                    }
                }
                for(int i=1; i<size; i++){
                    if(finished==0.0){
                        MPI_Send(&bestId, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                        MPI_Send(&best_tour_cost, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
                    }
                    else{
                        MPI_Send(&lowint, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                        MPI_Send(&best_tour_cost, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
                    }
                }
            }
            ////recebem informacao do melhor LB
            if(rank!=0){
                MPI_Recv(&bestId, 1, MPI_INT, root, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&best_tour_cost, 1, MPI_DOUBLE, root, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            ////melhor LB distribui
            if(rank==bestId){
                for(int i=0; i<size; i++){
                    if(i!=rank){
                        N2 = queue_pop(queue);
                        buffer_size = sizeof(int)*3 + sizeof(int*)*2 + sizeof(int)*N2->tourLength + sizeof(int)*(number_of_cities+2) + sizeof(double)*2;
                        buffer = (char*)malloc(buffer_size);
                        buffer = encode(buffer, N2);
                        MPI_Send(buffer, buffer_size, MPI_BYTE, i, 2, MPI_COMM_WORLD);

                    }
                }
            }
            else{
                if(bestId==-1){
                    finished=-1;
                }
                else{

                    MPI_Status status;
                    MPI_Probe(bestId, 2, MPI_COMM_WORLD, &status);
                    int count;
                    MPI_Get_count(&status, MPI_BYTE, &count);

                    char* buffer_receive = (char*)malloc(count);

                    MPI_Recv(buffer_receive, count, MPI_BYTE, bestId, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                    N2 = decode(buffer_receive);

                    queue_push(queue, N2);
                }
            }
            N = queue_pop(queue);

        }

        if (N->lower_bound >= best_tour_cost) {
            printf("Rank: %d finished!\n", rank);
            finished = -1.0;
            break;
        }

        if (N->tourLength == number_of_cities) {

            if (N->cost + matrix_weights[N->currCity][0] < best_tour_cost &&
                matrix_weights[N->currCity][0] != inf) {

                best_tour_cost = N->cost + matrix_weights[N->currCity][0];

                for (int i = 0; i < (N->tourLength); i++) {
                    final_tour[i] = N->tour[i];
                }
                final_tour[N->tourLength] = 0;

                printf("Rank: %d foundbest! btc: %.1lf\n", rank, best_tour_cost);
            }

        } else {
            for (v = 1; v < number_of_cities; v++) {

                cnt = cnt + 1;

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

        free_node(N);
    }

    if(rank!=0){
        MPI_Send(&lowfloat, 1, MPI_DOUBLE, root, 0, MPI_COMM_WORLD);
        MPI_Send(&best_tour_cost, 1, MPI_DOUBLE, root, 1, MPI_COMM_WORLD);
    }
    else{
        for(int i=1; i<size; i++){
            MPI_Send(&lowint, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&best_tour_cost, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
        }
    }

    printf("Rank: %d left while!\n", rank);

    MPI_Barrier(MPI_COMM_WORLD);

    bestId=0;
    double final_cost = 0;


    for (i = 0; i < number_of_cities; i++) {
        final_cost = final_cost + matrix_weights[final_tour[i]][final_tour[i + 1]];
    }

    if(rank!=0){
        MPI_Send(&final_cost, 1, MPI_DOUBLE, root, 7, MPI_COMM_WORLD);
    }
    else{
        best_tour_cost = final_cost;
        for(int i=1; i<size; i++){
            MPI_Recv(&bestcostvec[i-1], 1, MPI_DOUBLE, i, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        for(int i=1; i<size; i++){
            if(bestcostvec[i-1]<=best_tour_cost){
                bestId=i;
                best_tour_cost = bestcostvec[i-1];
            }
        }
        for(int i=1; i<size; i++){
            MPI_Send(&bestId, 1, MPI_INT, i, 8, MPI_COMM_WORLD);
        }
    }

    if(rank!=0)
        MPI_Recv(&bestId, 1, MPI_INT, root, 8, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    final_vec = calloc(number_of_cities + 2, sizeof(int));

    for(int i=0; i<number_of_cities+1; i++){
        final_vec[i+1] = final_tour[i];
    }

    for(int i=0; i<number_of_cities+1; i++){
        final_vec[i+1] = final_tour[i];
    }

    final_vec[0]=bestId;

    return final_vec;
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
    priority_queue_t *queue = NULL;
    double *min1, *min2;
    int *final_tour_cost;
    double final_cost = 0;
    int neighbours_of_zero = 0;
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    number_of_cities = 0;
    max_value_accepted = 0;

    if (argc != 3) {
        exit(1);
    }

    max_value_accepted = atoi(argv[2]);

    fptr1 = fopen(argv[1], "r");

    if (fptr1 == NULL) {
        exit(0);
    }

    if (fscanf(fptr1, "%d %d", &number_of_cities, &total_roads) != 2)
        return 0;

    // allocate size to matrix and its elements
    double **matrix_weights = (double **) malloc(number_of_cities * sizeof(double *));
    for (int i = 0; i < number_of_cities; i++) {
        matrix_weights[i] = (double *) malloc(number_of_cities * sizeof(double));
    }

    // initialize matrix with infinities
    for (int i = 0; i < number_of_cities; i++) {
        for (int j = 0; j < number_of_cities; j++) {
            matrix_weights[i][j] = inf;
        }
    }

    // populate nodes' edge weights
    while (fscanf(fptr1, "%d %d %lf", &i, &j, &weight) != EOF) {
        if(i==0){
            neighbours_of_zero++;
        }
        matrix_weights[i][j] = weight;
        matrix_weights[j][i] = weight;
    }
    fclose(fptr1);

    // allocate mem
    min1 = malloc(number_of_cities * sizeof(double));
    min2 = malloc(number_of_cities * sizeof(double));
    final_tour_cost = calloc(number_of_cities + 2, sizeof(int));

    calculate_mins(matrix_weights, number_of_cities, min1, min2);

    //////////////////////////////////////////////////////////////////////////
    double starttime, endtime;

    //It will only start counting the time when it passes through here
    MPI_Barrier(MPI_COMM_WORLD);
    starttime = MPI_Wtime();

    final_tour_cost = procedure(matrix_weights, number_of_cities, max_value_accepted, min1, min2, queue,
                                final_tour_cost, size, rank, neighbours_of_zero);

    if(rank==final_tour_cost[0]){
        if (final_tour_cost[2] != 0) {
            for (i = 1; i < number_of_cities+1; i++) {
                final_cost = final_cost + matrix_weights[final_tour_cost[i]][final_tour_cost[i + 1]];
            }
        }
        if (max_value_accepted > final_cost && final_cost != 0) {
            fprintf(stdout, "%.1f\n", final_cost);
            for (i = 1; i < number_of_cities + 2; i++) {
                fprintf(stdout, "%d ", final_tour_cost[i]);
            }
            fprintf(stdout, "\n");
        } else {
            fprintf(stdout, "NO SOLUTION\n");
        }
    }

    /////////////////////////////////////////////////////////////////////////////////
    MPI_Barrier(MPI_COMM_WORLD);
    endtime = MPI_Wtime();

    if(rank==0)
    fprintf(stderr, "%.1fs\n", endtime-starttime);

    free(min1);
    free(min2);

    for (i = 0; i < number_of_cities; i++) {
        free(matrix_weights[i]);
    }

    free(matrix_weights);
    free(final_tour_cost);


    MPI_Finalize();
    return 0;
}