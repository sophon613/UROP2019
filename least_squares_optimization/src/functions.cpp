#include <iostream>
#include <cmath>
#include <Eigen/Dense>

using namespace Eigen;

//Matrix generator function

MatrixXd matrixGen(int rows, int cols) {


	MatrixXd points = MatrixXd::Random(rows, cols);

	return points;
}

//Quaternion converter

MatrixXd quaternionConv2R(double ANGLE, int x, int y, int z) {

	double q[4] = { cos(ANGLE), sin(ANGLE) * x, sin(ANGLE) * y, sin(ANGLE) * z };

	double xx = q[1] * q[1];
	double xy = q[1] * q[2];
	double xz = q[1] * q[3];
	double xw = q[1] * q[0];

	double yy = q[2] * q[2];
	double yz = q[2] * q[3];
	double yw = q[2] * q[0];

	double zz = q[3] * q[3];
	double zw = q[3] * q[0];

	double m00 = 1 - 2 * (yy + zz);
	double m01 = 2 * (xy - zw);
	double m02 = 2 * (xz + yw);

	double m10 = 2 * (xy + zw);
	double m11 = 1 - 2 * (xx + zz);
	double m12 = 2 * (yz - xw);

	double m20 = 2 * (xz - yw);
	double m21 = 2 * (yz + xw);
	double m22 = 1 - 2 * (xx + yy);

	MatrixXd R(3, 3);

	R(0, 0) = m00;
	R(0, 1) = m01;
	R(0, 2) = m02;

	R(1, 0) = m10;
	R(1, 1) = m11;
	R(1, 2) = m12;

	R(2, 0) = m20;
	R(2, 1) = m21;
	R(2, 2) = m22;

	return R;
}

//Predefined transformation
MatrixXd predefinedTransformation(MatrixXd points, MatrixXd R, int translation[3]) {
	MatrixXd pointsNew = points * R;

	for (int i = 1; i < pointsNew.rows(); i++) {
		pointsNew(i, 0) += translation[0];
		pointsNew(i, 1) -= translation[1];
		pointsNew(i, 2) += translation[2];
	}
	return pointsNew;

}
//Finding centroids of two data clouds
MatrixXd centroidFinding(MatrixXd points, MatrixXd pointsNew) {
	double sum_1x = 0;
	double sum_1y = 0;
	double sum_1z = 0;

	double sum_2x = 0;
	double sum_2y = 0;
	double sum_2z = 0;

	for (int i = 1; i < points.rows(); i++) {
		sum_1x += points(i, 0);
		sum_1y += points(i, 1);
		sum_1z += points(i, 2);

		sum_2x += pointsNew(i, 0);
		sum_2y += pointsNew(i, 1);
		sum_2z += pointsNew(i, 2);
	}
	double avg_cloud_1x = sum_1x / points.rows();
	double avg_cloud_1y = sum_1y / points.rows();
	double avg_cloud_1z = sum_1z / points.rows();

	double avg_cloud_2x = sum_2x / points.rows();
	double avg_cloud_2y = sum_2y / points.rows();
	double avg_cloud_2z = sum_2z / points.rows();

	MatrixXd centroids(2, 3);

	centroids(0, 0) = avg_cloud_1x;
	centroids(0, 1) = avg_cloud_1y;
	centroids(0, 2) = avg_cloud_1z;

	centroids(1, 0) = avg_cloud_2x;
	centroids(1, 1) = avg_cloud_2y;
	centroids(1, 2) = avg_cloud_2z;


	return centroids;

}

MatrixXd findH(MatrixXd points, MatrixXd pointsNew, MatrixXd centroids) {
	MatrixXd diff1 = MatrixXd::Constant(points.rows(), 3, 0);
	MatrixXd diff2 = MatrixXd::Constant(points.rows(), 3, 0);

	MatrixXd H = MatrixXd::Constant(3, 3, 0);
	MatrixXd a = MatrixXd::Constant(3, 1, 0);
	MatrixXd b = MatrixXd::Constant(1, 3, 0);

	for (int i = 1; i < points.rows(); i++) {
		diff1(i, 0) = points(i, 0) - centroids(0, 0); //x
		diff1(i, 1) = points(i, 1) - centroids(0, 1); // y
		diff1(i, 2) = points(i, 2) - centroids(0, 2); // z

		diff2(i, 0) = pointsNew(i, 0) - centroids(1, 0);// x
		diff2(i, 1) = pointsNew(i, 1) - centroids(1, 1);// y
		diff2(i, 2) = pointsNew(i, 2) - centroids(1, 2);// z

		a << diff1(i, 0),
			diff1(i, 1),
			diff1(i, 2);

		b << diff2(i, 0), diff2(i, 1), diff2(i, 2);
		H = H + a * b;

	}
	return H;

}

//Jacobian SVD decomposition of matrix to recover R

MatrixXd singleValueDecomposition(MatrixXd H) {
	JacobiSVD<MatrixXd> svd(H, ComputeThinU | ComputeThinV);

	MatrixXd V = svd.matrixV();
	MatrixXd U = svd.matrixU();
	MatrixXd S = svd.singularValues();

	MatrixXd X = U * V.transpose();

	MatrixXd correction(3, 3);

	correction << 1, 0, 0,
		0, 1, 0,
		0, 0, X.determinant();


	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {

			if (abs(X(i, j)) < 1e-5) {
				X(i, j) = 0.00;
			}

		}
	}

	if (X.determinant() == 1) {
		std::cout << "Congratz" << std::endl;


	}
	else if (X.determinant() == -1) {
		std::cout << "Humm, check with Admin" << std::endl;
		X = U * correction * V.transpose();
	}
	return X;
}

MatrixXd findTranslation(MatrixXd centroids, MatrixXd X) {
	MatrixXd t = MatrixXd::Constant(3, 1, 0);
	MatrixXd centroid_2_Vector = MatrixXd::Constant(3, 1, 0);
	MatrixXd centroid_1_Vector = MatrixXd::Constant(3, 1, 0);

	centroid_1_Vector << round(centroids(0, 0)),
		round(centroids(0, 1)),
		round(centroids(0, 2));

	centroid_2_Vector << round(centroids(1, 0)),
		round(centroids(1, 1)),
		round(centroids(1, 2));

	t = centroid_2_Vector - X * centroid_1_Vector;
	return t;
}




