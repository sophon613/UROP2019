#include <iostream>
#include <cmath>
#include <Eigen/Dense>
using namespace Eigen;



int main() {
	int rows = 41472;

	//Extract point cloud points as matrix, for now use a random number generator
	MatrixXd points = MatrixXd::Random(rows, 3);
	
	
	//Quaternion conversion to R
	double ANGLE = 3.14 / 6;
	int x = 1;
	int y = 0;
	int z = 0;

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

	MatrixXd R(3,3);

	R(0, 0) = m00;
	R(0, 1) = m01;
	R(0, 2) = m02;

	R(1, 0) = m10;
	R(1, 1) = m11;
	R(1, 2) = m12;

	R(2, 0) = m20;
	R(2, 1) = m21;
	R(2, 2) = m22;

	std::cout << "\n" << R << std::endl;

	//Now we apply an arbitrary transformation to this matrix

	MatrixXd pointsNew = points * R;

	for (int i = 1; i < rows; i++) {
		pointsNew(i, 0) += 3;
		pointsNew(i, 1) -= 2;
		pointsNew(i, 2) += 1;
	}

	//Add some code to visualize what is happening, view the point cloud
	std::cout << "\n" << pointsNew.rows() << std::endl;
	//Now lets find the centroid

	
	double sum_1x = 0;
	double sum_1y = 0;
	double sum_1z = 0;
	
	double sum_2x = 0;
	double sum_2y = 0;
	double sum_2z = 0;

	for (int i = 1; i < rows; i++) {
		sum_1x += points(i, 0);
		sum_1y += points(i, 1);
		sum_1z += points(i, 2);

		sum_2x += pointsNew(i, 0);
		sum_2y += pointsNew(i, 1);
		sum_2z += pointsNew(i, 2);
	}

	double avg_cloud_1x = sum_1x / rows;
	double avg_cloud_1y = sum_1y / rows;
	double avg_cloud_1z = sum_1z / rows;

	double avg_cloud_2x = sum_2x / rows;
	double avg_could_2y = sum_2y / rows;
	double avg_cloud_2z = sum_2z / rows;

	double centroid_cloud_1[3] = { avg_cloud_1x, avg_cloud_1y, avg_cloud_1z };
	double centroid_cloud_2[3] = { avg_cloud_2x,avg_could_2y, avg_cloud_2z };

	std::cout << "\n" << centroid_cloud_1 << std::endl;
	std::cout << "\n" << centroid_cloud_2 << std::endl;

	//We must now investigate the difference matrix between each point and the PCL centroid

	//Initialize

	MatrixXd diff1 = MatrixXd::Constant(rows, 3, 0);
	MatrixXd diff2 = MatrixXd::Constant(rows, 3, 0);

	MatrixXd H = MatrixXd::Constant(3, 3, 0);
	MatrixXd a = MatrixXd::Constant(3, 1, 0);
	MatrixXd b = MatrixXd::Constant(1, 3, 0);

	//Loop through the matrix of the two point clouds to form two matrix

	for (int i = 1; i < rows; i++) {
		diff1(i, 0) = points(i, 0) - centroid_cloud_1[0]; //x
		diff1(i, 1) = points(i, 1) - centroid_cloud_1[1]; // y
		diff1(i, 2) = points(i, 2) - centroid_cloud_1[2]; // z

		diff2(i, 0) = pointsNew(i, 0) - centroid_cloud_2[0];// x
		diff2(i, 1) = pointsNew(i, 1) - centroid_cloud_2[1];// y
		diff2(i, 2) = pointsNew(i, 2) - centroid_cloud_2[2];// z

		a << diff1(i, 0),
			diff1(i, 1),
			diff1(i, 2);

		b << diff2(i, 0), diff2(i, 1), diff2(i, 2);
		H = H + a * b;

	}

	std::cout << "\n" << H << std::endl;

	//SVD decomposition

	JacobiSVD<MatrixXd> svd(H, ComputeThinU | ComputeThinV);

	MatrixXd V = svd.matrixV();
	MatrixXd U = svd.matrixU();
	MatrixXd S = svd.singularValues();

	MatrixXd X = U * V.transpose();

	MatrixXd correction = MatrixXd::Constant(3, 3, 0);

	correction << 1, 0, 0,
		0, 1, 0,
		0, 0, X.determinant();

	
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {

			if (abs(X(i, j) < 1e-5)) {
				X(i, j) = 0.00;
			}
		
			
		}
	}



	if (X.determinant() == 1) {
		std::cout << "Congratz" << std::endl;

	}
	else if (X.determinant() == -1) {
		std::cout << "Humm, check with Admin"<< std::endl;
		X = U * correction * V.transpose();
	}

	MatrixXd t = MatrixXd::Constant(3, 1, 0);
	MatrixXd centroid_2_Vector = MatrixXd::Constant(3, 1, 0);
	MatrixXd centroid_1_Vector = MatrixXd::Constant(3, 1, 0);

	centroid_1_Vector << round(centroid_cloud_1[0]),
		round(centroid_cloud_1[1]),
		round(centroid_cloud_1[2]);

	centroid_2_Vector << round(centroid_cloud_2[0]),
		round(centroid_cloud_2[1]),
		round(centroid_cloud_2[2]);

	t =centroid_2_Vector - X * centroid_1_Vector;

	std::cout << "\n" << X << std::endl;
	
	std::cout << "\n" << t << std::endl;

	std::cout << "\n" << X.determinant() << std::endl;


	

	std::cin.get();


}