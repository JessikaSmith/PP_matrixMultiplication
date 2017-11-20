#include "mpi.h"
#include <iostream>
#include <random>
#include <algorithm>

#define SIZE 1000
#define NUM_OF_SAMPLES 51
#define NUM_OF_PROC_PREDEFINED 2

MPI_Status status;

using namespace std;

mt19937 mt(1729);
normal_distribution<float> dist(-1, 1);

// random matrix type - True, else zero matrix
float** GenerateMatrix(bool type) {
	float** matrix = 0;
	matrix = new float*[SIZE];
	for (int i = 0; i < SIZE; i++) {
		matrix[i] = new float[SIZE];
		for (int j = 0; j < SIZE; j++) {
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


//double CountTime(void(*f)(unsigned)) {
//	int size = 1000;
//	// cout << "Enter matrix size n: ";
//	// cin >> size;
//	double t1, t2;
//	double* time = new double[NUM_OF_SAMPLES];
//	for (int timeSample = 0; timeSample < NUM_OF_SAMPLES; timeSample++) {
//		t1 = MPI_Wtime();
//		f();
//		t2 = MPI_Wtime();
//		time[timeSample] = t2 - t1;
//	}
//	sort(time, time + NUM_OF_SAMPLES);
//	return time[(NUM_OF_SAMPLES - 1) / 2];
//}

int main(int argc, char** argv) {


	// Serial computations
	// void(*func1)() = SerialMatrixMultiplication;
	// Serial time: 18.2206
	// cout << "Time for sequential calculation: " << CountTime(func1, NUM_OF_SAMPLES) << endl;
	
	int procNum, procRank;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);

	int numOfWorkers = procNum - 1;
	int rows = SIZE / numOfWorkers;

	MPI_Datatype dt_temp, dt_column;
	MPI_Type_vector(SIZE, rows, SIZE, MPI_FLOAT, &dt_temp);
	MPI_Type_create_resized(dt_temp, 0, sizeof(float), &dt_column);
	MPI_Type_commit(&dt_column);

	MPI_Barrier(MPI_COMM_WORLD);

	double start = MPI_Wtime();
	int offset = 0, cyclicShift, i, j, k = 0;
	float** matrix1;
	float** matrix2;
	float matrix3[SIZE][SIZE];

	if (procRank == 0) {
		matrix1 = GenerateMatrix(1);
		matrix2 = GenerateMatrix(1);

		// pass colums to the processes
		for (int r = 0; r < rows; r++) {
			offset += rows;
			MPI_Send(&matrix2[0][offset], 1, dt_column, r + 1, 0, MPI_COMM_WORLD);
		}

		// broadcast rows to the processes with cyclic shift
		// set cyclicShift
		for (int cyclicShift = 0; cyclicShift < rows; cyclicShift++) {
			for (int worker = 1; worker < numOfWorkers; worker++) {
				//MPI_Send(&cyclicShift, 1, MPI_INT, worker, 1, MPI_COMM_WORLD);
				MPI_Send(&offset, 1, MPI_INT, worker, 1, MPI_COMM_WORLD);
				MPI_Send(&matrix1[offset][0], rows*SIZE, MPI_FLOAT, worker, 1, MPI_COMM_WORLD);
				offset += rows;
			}
			for (int worker = 1; worker < numOfWorkers; worker++) {
				MPI_Recv(&offset, 1, MPI_INT, worker, 2, MPI_COMM_WORLD, &status);
				MPI_Recv(&matrix3[offset][0], rows*SIZE, MPI_FLOAT, worker, 2, MPI_COMM_WORLD, &status);

			}
		}

		// receive results in submatrix 

		// print result 
		printf("Here is the result matrix:\n");
		for (i = 0; i < SIZE; i++) {
			for (j = 0; j < SIZE; j++)
				printf("%6.2f   ", matrix3[i][j]);
			printf("\n");
		}

	}

	else {
		// get the columns and rows
		MPI_Recv(&offset, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
		MPI_Recv(&matrix2, 1, dt_column, 0, 1, MPI_COMM_WORLD, &status);
		MPI_Recv(&matrix1, rows*SIZE, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &status);

		// Calculations
		for (k = 0; k < rows; k++) {
			for (i = 0; i < rows; i++) {
				matrix3[i][k] = 0.0;
				for (j = 0; j < rows; j++) {
					matrix3[i][k] = matrix3[i][k] + matrix1[i][j] * matrix2[j][k];
				}
			}
		}

		MPI_Send(&offset, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
		MPI_Send(&matrix3, rows*rows, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
	}
	double end = MPI_Wtime();
	cout << "Computational time: " << end - start << endl;
	
	MPI_Type_free(&dt_column);
	MPI_Finalize();
}
