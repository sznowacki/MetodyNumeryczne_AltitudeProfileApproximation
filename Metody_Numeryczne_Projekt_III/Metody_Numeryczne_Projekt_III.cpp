// Metody_Numeryczne_Projekt_III.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Point.h"
#include <cstdlib>
#include <string>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

double* lagrangeFiFunction(int numberOfPoints, double x, Point** tableOfPoints)
{
	double* fiTable = new double[numberOfPoints];
	for (int i = 0; i < numberOfPoints; i++)
	{
//		//a/b, a - licznik, b - mianownik
		double a = 1, b = 1;
		for (int j = 0; j < numberOfPoints; j++)
		{
			if (i == j) {
				continue;
			}
			a = a * (x - tableOfPoints[j]->getXValue());
			b = b * (tableOfPoints[i]->getXValue() - tableOfPoints[j]->getXValue());
		}
		fiTable[i] = a / b;
		}
	return fiTable;
}

double lagrangeInterpolation(int numberOfPoints, double x, Point** tableOfPoints)
{
	double* fiTable = lagrangeFiFunction(numberOfPoints, x, tableOfPoints);
	
	double Fx = 0;
	for(int i = 0; i < numberOfPoints; i++)
	{
		Fx += tableOfPoints[i]->getYValue() * fiTable[i];
	}
	delete fiTable;
	return Fx;
}

double splineThirdDegreePolynominalFunction(Point** tableOfPoints, double* tableOfFactors, double x, int n)
{
	double a = tableOfFactors[4 * n];
	double b = tableOfFactors[4 * n + 1];
	double c = tableOfFactors[4 * n + 2];
	double d = tableOfFactors[4 * n + 3];
	double xDifference = x - tableOfPoints[n]->getXValue();
	double Sx = a + b*xDifference + c*(xDifference*xDifference) + d*(xDifference*xDifference*xDifference);

	return Sx;
}

void printMatrix(double** thirdDegreeMatrix, int numberOfPoints)
{
	int matrixSize = 4 * (numberOfPoints - 1);
	for (int j = 0; j < matrixSize; j++)
	{
		for (int i = 0; i < matrixSize; i++)
		{
			printf("%f ", thirdDegreeMatrix[j][i]);
		}
		printf("\n");
	}
}

double** createThirdDegreeMatrix(int numberOfPoints)
{
	int matrixSize = 4 * (numberOfPoints - 1);
	double** thirdDegreeMatrix = new double*[matrixSize];
	for (int i = 0; i < matrixSize; i++)
	{
		thirdDegreeMatrix[i] = new double[matrixSize];
	}

	for(int i = 0; i < matrixSize; i++)
	{
		for(int j = 0; j < matrixSize; j++)
		{
			thirdDegreeMatrix[j][i] = 0;
		}
	}
	return thirdDegreeMatrix;
}

void fillThirdDegreeMatrix(double** thirdDegreeMatrix, double h, int numberOfPoints)
{	
	//Filling upper part of matrix
	int matrixSize = 4 * (numberOfPoints - 1);

	for(int i = 0; i < 2*(numberOfPoints-1); i++)
	{
		if(i%2 == 0)
		{
			thirdDegreeMatrix[2 * i][i] = 1.0;
		}
		else if(i % 2 == 1)
		{
			thirdDegreeMatrix[2*i-2][i] = 1.0;
			thirdDegreeMatrix[2*i-1][i] = h;
			thirdDegreeMatrix[2*i][i] = h*h;
			thirdDegreeMatrix[2*i+1][i] = h*h*h;
		}
	}

	//Filling bottom part of matrix without last 2 rows
	for (int i = 0; i < 2 * (numberOfPoints - 1) - 2; i++)
	{
		if (i % 2 == 0)
		{
			thirdDegreeMatrix[2*i+1][i+2*(numberOfPoints-1)] = 1;
			thirdDegreeMatrix[2*i+2][i+ 2 * (numberOfPoints - 1)] = 2*h;
			thirdDegreeMatrix[2*i+3][i+ 2 * (numberOfPoints - 1)] = 3*h*h;
			thirdDegreeMatrix[2*i+5][i+ 2 * (numberOfPoints - 1)] = -1;
		}
		else if (i % 2 == 1)
		{
			thirdDegreeMatrix[2*i][i+ 2 * (numberOfPoints - 1)] = 2;
			thirdDegreeMatrix[2*i+1][i+ 2 * (numberOfPoints - 1)] = 6*h;
			thirdDegreeMatrix[2*i+4][i+ 2 * (numberOfPoints - 1)] = -2;
		}
	}
	//Filllig bottom part of matrix
	thirdDegreeMatrix[2][matrixSize - 2] = 2.0;
	thirdDegreeMatrix[matrixSize-2][matrixSize - 1] = 2.0;
	thirdDegreeMatrix[matrixSize-1][matrixSize - 1] = 6.0*h;
}

void transposeFillThirdDegreeMatrix(double** thirdDegreeMatrix, double h, int numberOfPoints)
{
	//Filling upper part of matrix
	int matrixSize = 4 * (numberOfPoints - 1);

	for (int i = 0; i < 2 * (numberOfPoints - 1); i++)
	{
		if (i % 2 == 0)
		{
			thirdDegreeMatrix[i][2 * i] = 1.0;
		}
		else if (i % 2 == 1)
		{
			thirdDegreeMatrix[i][2 * i - 2] = 1.0;
			thirdDegreeMatrix[i][2 * i - 1] = h;
			thirdDegreeMatrix[i][2 * i] = h*h;
			thirdDegreeMatrix[i][2 * i + 1] = h*h*h;
		}
	}

	//Filling bottom part of matrix without last 2 rows
	for (int i = 0; i < 2 * (numberOfPoints - 1) - 2; i++)
	{
		if (i % 2 == 0)
		{
			thirdDegreeMatrix[i + 2 * (numberOfPoints - 1)][2 * i + 1] = 1;
			thirdDegreeMatrix[i + 2 * (numberOfPoints - 1)][2 * i + 2] = 2 * h;
			thirdDegreeMatrix[i + 2 * (numberOfPoints - 1)][2 * i + 3] = 3 * h*h;
			thirdDegreeMatrix[i + 2 * (numberOfPoints - 1)][2 * i + 5] = -1;
		}
		else if (i % 2 == 1)
		{
			thirdDegreeMatrix[i + 2 * (numberOfPoints - 1)][2 * i] = 2;
			thirdDegreeMatrix[i + 2 * (numberOfPoints - 1)][2 * i + 1] = 6 * h;
			thirdDegreeMatrix[i + 2 * (numberOfPoints - 1)][2 * i + 4] = -2;
		}
	}
	//Filllig bottom part of matrix
	thirdDegreeMatrix[matrixSize - 2][2] = 2.0;
	thirdDegreeMatrix[matrixSize - 1][matrixSize - 2] = 2.0;
	thirdDegreeMatrix[matrixSize - 1][matrixSize - 1] = 6.0*h;
}

void fillTableOfValues(Point** tableOfPoints, double* tableOfValues, int numberOfPoints)
{
	tableOfValues[0] = tableOfPoints[0]->getYValue();
	for(int i = 1; i < numberOfPoints-1; i++)
	{
		tableOfValues[2 * i - 1] = tableOfPoints[i]->getYValue();
		tableOfValues[2 * i] = tableOfPoints[i]->getYValue();
	}
	tableOfValues[2*(numberOfPoints-1) - 1] = tableOfPoints[numberOfPoints-1]->getYValue();
}

double* createTable(int numberOfPoints)
{
	int tableSize = 4 * (numberOfPoints - 1);
	double* table = new double[tableSize];
	
	for(int i = 0; i < tableSize; i++)
	{
		table[i] = 0;
	}
	return table;
}

//LUPDecompose from wikipedia
void LUPDecompose(double **A, int N, int *P) {

	int i, j, k, imax;
	double maxA, absA;
	double* ptr = new double[N];

	for (i = 0; i <= N; i++)
		P[i] = i; //Unit permutation matrix, P[N] initialized with N

	for (i = 0; i < N; i++) {
		maxA = 0.0;
		imax = i;

		for (k = i; k < N; k++)
			if ((absA = fabs(A[k][i])) > maxA) {
				maxA = absA;
				imax = k;
			}
		
		if (imax != i) {
			//pivoting P
			j = P[i];
			P[i] = P[imax];
			P[imax] = j;

			for (int q = 0; q < N; q++) {
				ptr[q] = A[i][q];
			}

			for (int q = 0; q < N; q++) {
				A[i][q] = A[imax][q];
				A[imax][q] = ptr[q];
			}
			//counting pivots starting from N (for determinant)
			P[N]++;
		}

		for (j = i + 1; j < N; j++) {
			A[j][i] /= A[i][i];

			for (k = i + 1; k < N; k++)
				A[j][k] -= A[j][i] * A[i][k];
		}
	}
}
//LUPSolve from wikipedia
void LUPSolve(double **A, int *P, double *b, int N, double *x) {

	for (int i = 0; i < N; i++) {
		x[i] = b[P[i]];

		for (int k = 0; k < i; k++)
			x[i] -= A[i][k] * x[k];
	}

	for (int i = N - 1; i >= 0; i--) {
		for (int k = i + 1; k < N; k++)
			x[i] -= A[i][k] * x[k];

		x[i] = x[i] / A[i][i];
	}
}

int checkIntervalFunctionNumber(Point** tableOfPoints,double x, int numberOfPoints)
{
	if(x <= 1)
	{
		return 0;
	}

	int n = 0;
	for(int i = 0; i < numberOfPoints; i++)
	{
		if (x > tableOfPoints[i]->getXValue())
		{
			continue;
		}
		n = i - 1;
		break;
	}
	return n;
}

void splineFunctionsInterpolation(Point** tableOfPoints, int numberOfPoints, string fileName)
{
	double h = tableOfPoints[1]->getXValue() - tableOfPoints[0]->getXValue();
	int matrixSize = 4 * (numberOfPoints - 1);
	double** thirdDegreeMatrix = createThirdDegreeMatrix(numberOfPoints);
	double* tableOfFactors = createTable(numberOfPoints);
	double* tableOfYValues = createTable(numberOfPoints);
	int* tableOfPermutation = new int[matrixSize + 1];
	transposeFillThirdDegreeMatrix(thirdDegreeMatrix, h, numberOfPoints);
	fillTableOfValues(tableOfPoints, tableOfYValues, numberOfPoints);
	LUPDecompose(thirdDegreeMatrix, matrixSize, tableOfPermutation);
	LUPSolve(thirdDegreeMatrix, tableOfPermutation, tableOfYValues, matrixSize, tableOfFactors);
	printMatrix(thirdDegreeMatrix, numberOfPoints);
	fstream file;
	file.open(fileName, ios::out | ios::app);
	int range = int(tableOfPoints[numberOfPoints - 1]->getXValue());
	double step = tableOfPoints[1]->getXValue() - tableOfPoints[0]->getXValue();
	for (double i = 0; i <= tableOfPoints[numberOfPoints - 1]->getXValue(); i=i+5)
	{
		int n = checkIntervalFunctionNumber(tableOfPoints, i, numberOfPoints);
		double y = splineThirdDegreePolynominalFunction(tableOfPoints, tableOfFactors, i, n);
//		printf("F(%f) = %f\n", i, y);
		file << i << "," << y << "\n";
	}
	file.close();
}

int countLines(string fileName)
{
	fstream file;
	//Zliczenie linii w pliku
	file.open(fileName, ios::in);

	if (file.good() == false)
	{
		printf("File doesn't exists");
		exit(0);
	}
	string line;
	int numberOfLines = 1;

	while (getline(file, line))
	{
		numberOfLines++;
	}
	file.close();

	return numberOfLines - 1;
}

void readFromFile(string fileName, Point** tableOfPoints)
{
	fstream file;
	//Odczytanie wartosci punktow z pliku
	file.open(fileName, ios::in);
	if (file.good() == false)
	{
		printf("File doesn't exists");
		exit(0);
	}
	string distance, height;
	int j = 1;
	while (getline(file, distance, ','))
	{
		double x = atof(distance.c_str());
		getline(file, height);
		double y = atof(height.c_str());
		//printf("%d:(%f,%f)\n", j, x, y);
		tableOfPoints[j - 1] = new Point(x, y);
		j++;
	}
	file.close();
}

void saveToFileLagrange(int numberOfLines, Point** tableOfPoints, string fileName)
{
	fstream file;
	file.open(fileName, ios::out | ios::app);
	int range = int(tableOfPoints[numberOfLines - 1]->getXValue());
	double step = tableOfPoints[1]->getXValue() - tableOfPoints[0]->getXValue();
	for (double i = 0; i <= range; i=i+5)
	{
		double y = lagrangeInterpolation(numberOfLines - 1, i, tableOfPoints);
//		printf("F(%f) = %f\n", i, y);
		file << i << "," << y << "\n";
	}
	file.close();
}

int main()
{
	string filesNames[5] = { "test.txt", "GlebiaChallengera.txt", "MountEverest.txt", "SpacerniakGdansk.txt", "WielkiKanionKolorado.txt" };
	string filesToSaveNames[5] = { "test.csv", "GlebiaChallengeraLagrange.csv", "MountEverestLagrange.csv", "SpacerniakGdanskLagrange.csv", "WielkiKanionKoloradoLagrange.csv" };
	string filesToSaveNames2[5] = { "test.csv", "GlebiaChallengeraSpline.csv", "MountEverestSpline.csv", "SpacerniakGdanskSpline.csv", "WielkiKanionKoloradoSpline.csv" };
//	for (int i = 0; i < 5; i++) {
		int numberOfLines = countLines(filesNames[0]);
		Point** tableOfPoints = new Point*[numberOfLines];
		readFromFile(filesNames[0], tableOfPoints);
		saveToFileLagrange(numberOfLines, tableOfPoints, filesToSaveNames[0]);
		splineFunctionsInterpolation(tableOfPoints, numberOfLines, filesToSaveNames2[0]);
//	}
    return 0;
}

