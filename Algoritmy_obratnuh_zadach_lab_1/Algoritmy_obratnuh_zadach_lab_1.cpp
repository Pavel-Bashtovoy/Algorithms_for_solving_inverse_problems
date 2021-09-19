#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
using namespace std;
const double PI = 3.141593;
struct point3d
{
	double x, y, z;
};
double r(point3d first, point3d second)
{
	return sqrt
	(
		(first.x - second.x) * (first.x - second.x) +
		(first.y - second.y) * (first.y - second.y) +
		(first.z - second.z) * (first.z - second.z)
	);
}
void readSettings(int& max_iteration, double& eps, double& sigma, double
	& sigma_iteration, string fileName)
{
	ifstream file(fileName);
	file >> max_iteration;
	file >> eps;
	file >> sigma;
	file >> sigma_iteration;
	file.close();
}
void readModel(pair <point3d, point3d> &source, vector<pair<point3d, point3d>> &receivers, double  &I, string fileName)
{
	int count;
	point3d A, B, M, N;
	ifstream file(fileName);
	file >> count;
	for (int i = 0; i < count; i++)
	{
		file >> M.x >> M.y >> M.z >> N.x >> N.y >> N.z;
		receivers.push_back(pair<point3d, point3d>(M, N));
	}
	file >> A.x >> A.y >> A.z >> B.x >> B.y >> B.z;
	source = pair <point3d, point3d>(A, B);
	file >> I;
	file.close();
}
double potential(double I, double sigma, pair<point3d, point3d> source, pair<point3d,
	point3d> receiver)
{
	point3d A = source.first;
	point3d B = source.second;
	point3d M = receiver.first;
	point3d N = receiver.second;
	return (I / (2 * PI * sigma)) * ((1.0 / r(B, M) - 1.0 / r(A, M)) - (1.0 / r(B, N) - 1.0
		/ r(A, N)));
}
double functional(vector<double> V, vector<double> V_iteration, int size)
{
	double w, result = 0.0;
	for (int i = 0; i < size; i++)
	{
		w = 1.0 / V[i];
		result += pow((w * (V_iteration[i] - V[i])), 2);
	}
	return result;
}
double diff_sigma(double I, double sigma, pair<point3d, point3d> source, pair<point3d,
	point3d> receiver)
{
	point3d A = source.first;
	point3d B = source.second;
	point3d M = receiver.first;
	point3d N = receiver.second;
	return -(I / (2 * PI * sigma * sigma)) * ((1.0 / r(B, M) - 1.0 / r(A, M)) - (1.0 / r(B,
		N) - 1.0 / r(A, N)));
}
double GetRandomDouble(double RealMax)
{
	double g = (1.0 * rand() / RAND_MAX);
	if (g > 0.999)
	{
		g = 0.999;
	}
	return (g > 0.5 ? RealMax * g : -RealMax * g);
}
void main()
{
	srand((unsigned)time(NULL));
	int iteration, max_iteration;
	double F = 0, eps, sigma, sigma_iteration, I, w, A, B, dSigma, a, b ;

	vector<double> V, V_iteration_sigma;
	pair<point3d, point3d>source;
	vector<pair<point3d, point3d>> receivers;
	readSettings(max_iteration, eps, sigma, sigma_iteration, "settings");
	readModel(source, receivers, I, "model");
	for (int i = 0; i < receivers.size(); i++)
	{
		V.push_back(potential(I, sigma, source, receivers[i]));
		V_iteration_sigma.push_back(potential(I, sigma_iteration, source, receivers[i]) *
			(1.0 + GetRandomDouble(0.1)));
	}
	for (iteration = 1; iteration <= max_iteration; iteration++)
	{
		A = B = 0;
		for (int i = 0; i < receivers.size(); i++)
		{
			w = 1.0 / V[i];
			dSigma = diff_sigma(I, sigma_iteration, source, receivers[i]);
			A += pow((w * dSigma), 2);
			B -= pow(w, 2) * dSigma * (V_iteration_sigma[i] - V[i]);
		}
		sigma_iteration += B / A;
		for (int i = 0; i < receivers.size(); i++)
		{
			V_iteration_sigma[i] = potential(I, sigma_iteration, source, receivers[i]);
		}
		F = functional(V, V_iteration_sigma, receivers.size());
		if (F < eps)
		{
			break;
		}
	}
	if (iteration > max_iteration)
	{
		iteration = max_iteration;
	}
	cout << "Functional: "<<F<< endl << "Itearations: " << iteration << endl << endl << "Sigma: " << sigma_iteration << endl;
	system("pause");
}