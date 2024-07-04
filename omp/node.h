#ifndef _TSP_NODE_H

typedef struct
{
    int *tour; // curr path
    int tourLength;
    double lower_bound;
    double cost; // current tour cost (or not)
    int currCity; // node's city
    int *viewed;
    int n_cities;
} node;

node* create_node(
    int *tour2,
    int tourLength2,
    double lower_bound2,
    double cost2,
    int currCity2, int n_cities);

void free_node(node *a);

char compare_node(void *a, void *b);

int is_on_tour(int *tour, int vertex, int length);

void print_tour(int *tour, int length);

#endif //_TSP_NODE_H