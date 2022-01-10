#include <iostream>
#include <Eigen/core>
#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <vector>

#include <opencv2/core.hpp>
using namespace cv;
using namespace std;
void PSVDSolveRT(const std::vector<Point3d>& pts1,
	const std::vector<Point3d>& pts2,
	Mat&R, Mat& T);


int main()
{
	cout << "hello world__" << endl;
	return 0;
}

void PSVDSolveRT(const std::vector<Point3d>& pts1,
	const std::vector<Point3d>& pts2,
	Mat& R, Mat& T)
{
	//点集中心
	Point3d p1;
	Point3d p2;

	int N = pts1.size();

	for (int i = 0; i < N; i++)
	{
		p1 += pts1[i];
		p2 += pts2[i];
	}

	//求解质心
	p1 = Point3d(Vec3d(p1) / N);
	p2 = Point3d(Vec3d(p2) / N);
	vector<Point3d> q1(N), q2(N);
	
	//去掉质心
	for (int i = 0; i < N; i++)
	{
		q1[i] = pts1[i] - p1;
		q2[i] = pts2[i] - p2;
	}

	//计算 q1*q2^T
	Eigen::Matrix3d W = Eigen::Matrix3d::Zero();
	for (int i = 0; i < N; i++)
	{
		W += Eigen::Vector3d( q1[i].x, q1[i].y, q1[i].z) * Eigen::Vector3d(q2[i].x, q2[i].y, q2[i].z).transpose();
	}

	cout << "W=" << W << endl;

	//SVD on W
	Eigen::JacobiSVD<Eigen::Matrix3d> svd(W, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Matrix3d U = svd.matrixU();
	Eigen::Matrix3d V = svd.matrixV();

	if (U.determinant() * V.determinant() < 0)
	{
		for (int x = 0; x < 3; x++)
		{
			U(x, 2) *= -1;
		}
	}

	cout << "U=" << U << endl;
	cout << "V=" << V << endl;
	
	Eigen::Matrix3d R_ = U * (V.transpose());
	Eigen::Vector3d t_ = Eigen::Vector3d(p1.x, p1.y, p1.z) - R_ * Eigen::Vector3d(p2.x, p2.y, p2.z);

	R = (Mat_<double>(3, 3) <<
		R_(0, 0), R_(0, 1), R_(0, 2),
		R_(1, 0), R_(1, 1), R_(1, 2),
		R_(2, 0), R_(2, 1), R_(2, 2)
		);
	T = (Mat_<double>(3, 1) << t_(0, 0), t_(1, 0), t_(2, 0));
}