#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <random>

using namespace std;

// Definir la estructura de un punto en el espacio
typedef struct {
    double x;
    double y;
} Point;

// Función para generar un número aleatorio entre min y max
double random_double(double min, double max, int seed) {
    mt19937_64 drng;
    drng.seed(seed);
    uniform_real_distribution<double> dist(min, max);
    return min + (dist(drng) / max) * (max - min);
}

// Función para calcular la distancia euclidiana entre dos puntos
double euclidean_distance(Point p1, Point p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    return sqrt(dx * dx + dy * dy);
}

// Función para asignar cada punto al centroide más cercano
void assign_points(Point* points, int num_points, Point* centroids, int num_centroids, int* assignments, int iteration) {
    if(iteration < 30){
        printf("asignando puntos...\n");
    }
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
void update_centroids(Point* points, int num_points, Point* centroids, int num_centroids, int* assignments, int iteration) {
    int* counts = (int*)calloc(num_centroids, sizeof(int));
    double* sum_x = (double*)calloc(num_centroids, sizeof(double));
    double* sum_y = (double*)calloc(num_centroids, sizeof(double));
    if(iteration < 30){
        printf("Actualizando centroides...\n\n");
    }
    
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
    remove("resultados.txt");
    FILE* file = fopen("resultados.txt", "a");
    if (file == NULL) {
        printf("No se pudo crear el archivo.\n");
        return;
    }
    fprintf(file, "Punto;CoordX;CoordY;Centroide\n");
    for (int i = 0; i < num_points; i++) {
       fprintf(file, "P%d;%.2f;%.2f;%d\n", i + 1, points[i].x, points[i].y, assignments[i]);
    }
    
    printf("Guardando resultados finales.......\n");
    fprintf(file,"Guardando resultados finales.......\n");
    for (int i = 0; i < num_centroids; i++) {
        fprintf(file, "Centroide %d: (%.2f, %.2f)\n", i + 1, centroids[i].x, centroids[i].y);
    }
    fclose(file);
}

int main() {

    // Configuración del algoritmo
    int seed=1234;
    int num_points = 1000;          // Número de puntos
    int num_centroids = 2;         // Número de centroides
    int max_iterations = 30;      // Número máximo de iteraciones
    double min_value = 0.0;        // Valor mínimo para generar puntos aleatorios
    double max_value = 10.0;       // Valor máximo para generar puntos aleatorios
    
    // Creación de los puntos aleatorios
    Point* points = (Point*)malloc(num_points * sizeof(Point));
    printf("Generando datos aleatorios... \n\n");
    for (int i = 0; i < num_points; i++) {
        points[i].x = random_double(min_value, max_value, seed+i);
        points[i].y = random_double(min_value, max_value, seed+i);
    }
    
    // Creación de los centroides iniciales
    Point* centroids = (Point*)malloc(num_centroids * sizeof(Point));
    printf("Inicializando centroides... \n\n");
    for (int i = 0; i < num_centroids; i++) {
        centroids[i].x = random_double(min_value, max_value, seed-i);
        centroids[i].y = random_double(min_value, max_value, seed-i);
    }
    
    // Asignación inicial de puntos a centroides
    int* assignments = (int*)malloc(num_points * sizeof(int));

    // Medir el tiempo de ejecución
    clock_t start_time = clock();

    // Bucle principal del algoritmo
    int iteration = 0;
    while (iteration < max_iterations) {
        if(max_iterations < 30){
            printf("Iteracion %i\n", iteration+1);
        }
       
        assign_points(points, num_points, centroids, num_centroids, assignments, max_iterations);
        update_centroids(points, num_points, centroids, num_centroids, assignments, max_iterations);
        iteration++;
    }
    
    // Calcular el tiempo transcurrido en segundos
    clock_t end_time = clock();
    double execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    // Imprimir resultados
    printf("Tiempo de ejecución: %.2f segundos\n", execution_time);
    print_results(points, num_points, centroids, num_centroids, assignments);
    
    // Liberar memoria
    free(points);
    free(centroids);
    free(assignments);
    
    return 0;
}