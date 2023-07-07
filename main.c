#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Definir la estructura de un punto en el espacio
typedef struct {
    double x;
    double y;
} Point;

// Función para generar un número aleatorio entre min y max
double random_double(double min, double max) {
    return min + ((double)rand() / RAND_MAX) * (max - min);
}

// Función para calcular la distancia euclidiana entre dos puntos
double euclidean_distance(Point p1, Point p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    return sqrt(dx * dx + dy * dy);
}

// Función para asignar cada punto al centroide más cercano
void assign_points(Point* points, int num_points, Point* centroids, int num_centroids, int* assignments) {
    for (int i = 0; i < num_points; i++) {
        double min_distance = INFINITY;
        int centroid_index = 0;
        
        for (int j = 0; j < num_centroids; j++) {
            double distance = euclidean_distance(points[i], centroids[j]);
            
            if (distance < min_distance) {
                min_distance = distance;
                centroid_index = j;
            }
        }
        
        assignments[i] = centroid_index;
    }
}

// Función para recalcular las posiciones de los centroides
void update_centroids(Point* points, int num_points, Point* centroids, int num_centroids, int* assignments) {
    int* counts = (int*)calloc(num_centroids, sizeof(int));
    double* sum_x = (double*)calloc(num_centroids, sizeof(double));
    double* sum_y = (double*)calloc(num_centroids, sizeof(double));
    
    for (int i = 0; i < num_points; i++) {
        int centroid_index = assignments[i];
        counts[centroid_index]++;
        sum_x[centroid_index] += points[i].x;
        sum_y[centroid_index] += points[i].y;
    }
    
    for (int i = 0; i < num_centroids; i++) {
        if (counts[i] > 0) {
            centroids[i].x = sum_x[i] / counts[i];
            centroids[i].y = sum_y[i] / counts[i];
        }
    }
    
    free(counts);
    free(sum_x);
    free(sum_y);
}

// Función para imprimir los resultados del algoritmo
void print_results(Point* points, int num_points, Point* centroids, int num_centroids, int* assignments) {
    for (int i = 0; i < num_points; i++) {
        printf("Punto %d: (%.2f, %.2f) asignado al centroide %d\n", i + 1, points[i].x, points[i].y, assignments[i]);
    }
    
    printf("\nCentroides finales:\n");
    for (int i = 0; i < num_centroids; i++) {
        printf("Centroide %d: (%.2f, %.2f)\n", i + 1, centroids[i].x, centroids[i].y);
    }
}

int main() {
    // Configuración del algoritmo
    int num_points = 10;          // Número de puntos
    int num_centroids = 2;         // Número de centroides
    int max_iterations = 3;      // Número máximo de iteraciones
    double min_value = 0.0;        // Valor mínimo para generar puntos aleatorios
    double max_value = 10.0;       // Valor máximo para generar puntos aleatorios
    
    // Creación de los puntos aleatorios
    Point* points = (Point*)malloc(num_points * sizeof(Point));
    for (int i = 0; i < num_points; i++) {
        points[i].x = random_double(min_value, max_value);
        points[i].y = random_double(min_value, max_value);
    }
    
    // Creación de los centroides iniciales
    Point* centroids = (Point*)malloc(num_centroids * sizeof(Point));
    for (int i = 0; i < num_centroids; i++) {
        centroids[i].x = random_double(min_value, max_value);
        centroids[i].y = random_double(min_value, max_value);
    }
    
    // Asignación inicial de puntos a centroides
    int* assignments = (int*)malloc(num_points * sizeof(int));
    
    // Bucle principal del algoritmo
    int iteration = 0;
    while (iteration < max_iterations) {
        assign_points(points, num_points, centroids, num_centroids, assignments);
        update_centroids(points, num_points, centroids, num_centroids, assignments);
        iteration++;
    }
    
    // Imprimir resultados
    print_results(points, num_points, centroids, num_centroids, assignments);
    
    // Liberar memoria
    free(points);
    free(centroids);
    free(assignments);
    
    return 0;
}