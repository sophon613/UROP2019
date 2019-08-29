#include <iostream>
#include <cmath>
#include <Eigen/Dense>

using namespace Eigen;

MatrixXd matrixGen(int rows, int cols);
MatrixXd quaternionConv2R(double ANGLE, int x, int y, int z);
MatrixXd predefinedTransformation(MatrixXd points, MatrixXd R, int translation[3]);
MatrixXd centroidFinding(MatrixXd points, MatrixXd pointsNew);
MatrixXd findH(MatrixXd points, MatrixXd pointsNew, MatrixXd centroids);
MatrixXd singleValueDecomposition(MatrixXd H);
MatrixXd findTranslation(MatrixXd centroids, MatrixXd X);

int main() {
	int rows = 41472;
	int cols = 3;

	//Generate a matrix with random numbers(later we will use PCL)
	MatrixXd points = matrixGen(rows, cols); //Generates the points matrix

	//Quaternion conversion to R
	double ANGLE = 3.14 / 6;
	int x = 1;
	int y = 0;
	int z = 0;

	MatrixXd R = quaternionConv2R(ANGLE, x, y, z);

	//Apply arbitrary transformation to the points

	int translation[3] = { 3,2,1 };

	MatrixXd pointsNew = predefinedTransformation(points,  R, translation);

	// Find the centroids

	MatrixXd centroids = centroidFinding(points, pointsNew);

	std::cout << "\n" << centroids << std::endl;

	// Find H matrix

	MatrixXd H = findH(points, pointsNew, centroids);

	std::cout << "\n" << H << std::endl;

	MatrixXd X = singleValueDecomposition(H);

	std::cout << "\n" << X << std::endl;

	//Find the t vector

	MatrixXd t = findTranslation(centroids, X);

	std::cout << "\n" << t << std::endl;

	std::cin.get();
}