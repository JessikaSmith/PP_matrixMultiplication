#include "mpi.h"
#include <iostream>
#include <random>
#include <algorithm>

#define SIZE = 1000
#define NUM_OF_SAMPLES = 51

using namespace std;

mt19937 mt(1729);
normal_distribution<float> dist(-1, 1);


// random matrix type - True, else zero matrix
float** GenerateMatrix(bool type) {
	float** matrix = 0;
	matrix = new float*[SIZE];
	for (int i = 0; i < SIZE; i++) {
		matrix[i] = new float[SIZE];
		for (int j = 0; j < N; j++) {
			if (type) matrix[i][j] = dist(mt);
			else matrix[i][j] = 0;
		}
	}
	return matrix;
}

void SerialMatrixMultiplication() {
	float** matrix1 = GenerateMatrix(1);
	float** matrix2 = GenerateMatrix(1);
	float** matrix3 = GenerateMatrix(0);
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			for (int k = 0; k < SIZE; k++) {
				matrix3[i][j] += matrix1[i][j] * matrix2[i][j];
			}
		}
	}
}

// Column-wize Block-Striped Decomposition
void ParallelCalculation() {
	float** matrix1 = GenerateMatrix(1);
	float** matrix2 = GenerateMatrix(1);
	float** matrix3 = GenerateMatrix(0);
	int procNum, procRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);
	MPI_Status status;
	// get the number of processes p
	// float **storage_mat = malloc(sizeof(float) * SIZE / p);
	blockSize = SIZE / procNum;
	if(proc)
}

double CountTime(void(*f)(unsigned)) {
	int size = 1000;
	// cout << "Enter matrix size n: ";
	// cin >> size;
	double t1, t2;
	double* time = new double[NUM_OF_SAMPLES];
	for (int timeSample = 0; timeSample < NUM_OF_SAMPLES; timeSample++) {
		t1 = MPI_Wtime();
		f();
		t2 = MPI_Wtime();
		time[timeSample] = t2 - t1;
	}
	sort(time, time + NUM_OF_SAMPLES);
	return time[(NUM_OF_SAMPLES - 1) / 2];
}

int main(int argc, char** argv) {

	MPI_Init(&argc, &argv);

	// creating derived column datatype
	MPI_Datatype dt_temp, dt_column;
	MPI_Type_vector(SIZE, 1, SIZE, MPI_DOUBLE, &dt_temp);
	MPI_Type_create_resized(dt_temp, 0, sizeof(double), &dt_column);
	MPI_Type_commit(&dt_column);

	// Serial computations
	// void(*func1)() = SerialMatrixMultiplication;
	// Serial time: 18.2206
	// cout << "Time for sequential calculation: " << CountTime(func1, NUM_OF_SAMPLES) << endl;

	MPI_Type_free(&dt_column)
	MPI_Finalize();
}
