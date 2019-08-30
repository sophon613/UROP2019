#include <Eigen/Dense>
#include <iostream>
#include <cmath>

MatrixXd matrixGen(int rows, int cols);
MatrixXd quaternionConv2R(double ANGLE, int x, int y, int z);
MatrixXd predefinedTransformation(MatrixXd points, MatrixXd R, int translation[3]);
MatrixXd centroidFinding(MatrixXd points, MatrixXd pointsNew);
MatrixXd findH(MatrixXd points, MatrixXd pointsNew, MatrixXd centroids);
MatrixXd singleValueDecomposition(MatrixXd H);
MatrixXd findTranslation(MatrixXd centroids, MatrixXd X);