#include "queue.h"
#include "node.h"

node* create_node(
    int *tour2,
    int tourLength2,
    double lower_bound2,
    double cost2,
    int currCity2, int n_cities)
{
    node *Node;

    Node = malloc(sizeof(node));

    Node->tour = (int*)malloc((tourLength2)*sizeof(int));
    Node->viewed = (int*)calloc(n_cities,sizeof(int));
    
    for(int i=0; i<tourLength2; i++){
        Node->tour[i] = tour2[i];
        Node->viewed[tour2[i]] = 1;
    }
    
    Node->tourLength = tourLength2;
    Node->lower_bound = lower_bound2;
    Node->cost = cost2;
    Node->currCity = currCity2;

    return Node;
}

void free_node(node *a){

    free(a->tour);
    free(a->viewed);
    free(a);
}

char compare_node(void *a, void *b){

    node *c = a;
    node *d = b;

    if(c->lower_bound > d->lower_bound){
        return 1;
    }
    else if(c->lower_bound == d->lower_bound){
        if(c->currCity > d->currCity){
            return 1;
        }
    }

    return 0;
}

void print_tour(int *tour, int length){

    printf("[");
    for(int i=0; i<length; i++){
        printf("%d ", tour[i]);
    }
    printf("]\n");
}