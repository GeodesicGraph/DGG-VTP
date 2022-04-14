#ifndef GEODESIC_RESULT_PROCESS_H
#define GEODESIC_RESULT_PROCESS_H

#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include "geodesic_mesh_elements.h"

class GeodesicResultProcess
{
public:
	GeodesicResultProcess(){}
	~GeodesicResultProcess(){}

	bool GeodesicDistanceResultExact(std::string file1, std::string file2);
	double GeodesicDistanceApproximateError(std::string exactFile, std::string apprFile, std::string errorFile);
	void GeodesicResultTexture(std::string scrData, std::string distanceData);
	double GeodesicResultRelativeErrorValueTexture(std::string scrData, std::string errorValueFile);
	unsigned FindNearestVertexIndex(geodesic::Point3D scrVert, std::string inFile);
};

bool GeodesicResultProcess::GeodesicDistanceResultExact(std::string file1, std::string file2)
{
	double eps = 1e-6;
	std::vector<double> filesrc;
	std::ifstream fin(file1);
	if (!fin)
	{
		printf("error: unable to open input file: %s\n", file1.c_str());
		return false;
	}
	filesrc.clear();
	double res;
	while (fin >> res)
	{

		filesrc.push_back(res);
	}
	fin.close();
	fin.open(file2);
	if (!fin)
	{
		printf("error: unable to open the input file: %s\n", file2.c_str());
		return false;
	}
	int i = 0;
	//int cn = 0;
	bool IsExact = true;
	while (fin >> res)
	{
		//printf("%.8lf\n", res);
		//cn++;
		//if (cn == 10) return;
		if (fabs(res - filesrc[i])/filesrc[i] > eps)
		{
			IsExact = false;
			printf("i: %d, src: %.20lf, res: %.20lf\n", i, filesrc[i], res);
		}
		i++;
	}

	return IsExact;
}

double GeodesicResultProcess::GeodesicDistanceApproximateError(std::string exactFile, std::string apprFile, std::string errorFile)
{

	double eps = 1e-12;
	std::vector<double> exactRes;
	std::ifstream fin(exactFile);
	if (!fin)
	{
		printf("Error: unable to open input file: %s\n", exactFile.c_str());
		return -1;
	}
	exactRes.clear();
	double res;
	while (fin >> res)
	{
		exactRes.push_back(res);
	}
	fin.close();
	fin.open(apprFile);
	if (!fin)
	{
		printf("Error: unable to open input file: %s\n", apprFile.c_str());
		return -1;
	}
	int i = 0;
	double avg_relative_error = 0;
	double relative_error = 0;

	
	std::ofstream fout(errorFile.c_str());

	while (fin >> res)
	{
		if (exactRes[i] < eps)
		{
			i++;
			relative_error = 0;
			fout << std::setprecision(20) << relative_error << std::endl;
			continue;
		}
		relative_error = (res - exactRes[i]) / exactRes[i];
		if (relative_error < 0)
		{
			printf("error: appr_result is smaller than exact_result!\n");
			system("pause");
			return -1;
		}
		fout << std::setprecision(20) << relative_error << std::endl;
		avg_relative_error += relative_error;
		i++;
	}
	fin.close();
	fout.close();
	avg_relative_error /= double(i);

	printf("Relative Error: %.6lf%%\n", avg_relative_error * 100);
	return avg_relative_error * 100;
}

void GeodesicResultProcess::GeodesicResultTexture(std::string scrData, std::string distanceData)// example: scrData: bunny3k.obj, distanceData: bunny3k_result.txt
{
	std::string outFile = distanceData.substr(0, distanceData.length() - 4);
	outFile.append("_texture.obj");
	double mxDis = -1.0;
	std::vector<double> distance;
	FILE *fin = fopen(distanceData.c_str(), "r");
	if (!fin)
	{
		printf("error: unable to open the input file: %s\n", distanceData.c_str());
		return;
	}
	distance.clear();
	double dis;
	while (fscanf(fin, "%lf", &dis) != EOF)
	{
		mxDis = max(dis, mxDis);
		distance.push_back(dis);
	}
	fclose(fin);
	for (int i = 0; i < distance.size(); i++)
	{
		distance[i] = distance[i] / mxDis;
	}
	fin = fopen(scrData.c_str(), "r");
	if (!fin)
	{
		printf("error: unable to open the input file: %s\n", scrData.c_str());
		return;
	}
	FILE *fout = fopen(outFile.c_str(), "w");
	if (!fout)
	{
		printf("error: unable to open the output file: %s\n", outFile.c_str());
		return;
	}
	//int n_v, n_f;
	//fscanf(fin, "%d %d", &n_v, &n_f);
	//printf("n_v: %d, n_f: %d\n", n_v, n_f);
	double x1, y1, z1;
	int f1, f2, f3;
	char ch = ' ';
	int i = 0;
	char str[255] = { '\0' };
	while (fscanf(fin, "%c", &ch) != EOF)
	{
		if (ch == 'v')
		{
			fscanf(fin, "%lf %lf %lf", &x1, &y1, &z1);
			fprintf(fout, "v %.18lf %.18lf %.18lf\n", x1, y1, z1);
			fprintf(fout, "vt %.18lf %.18lf\n", distance[i], distance[i]);
			i++;
		}
		else if (ch == 'f')
		{
			fscanf(fin, "%d %d %d\n", &f1, &f2, &f3);
			fprintf(fout, "f %d/%d %d/%d %d/%d\n", f1, f1, f2, f2, f3, f3);
		}
		else if (ch == '\n')
		{
			continue;
		}
		else
		{
			fgets(str, 255, fin);
		}
	}
	fclose(fin);
	fclose(fout);
}

//10%, 1%, 0.1%, 0.01%, 0.001%
double GeodesicResultProcess::GeodesicResultRelativeErrorValueTexture(std::string scrData, std::string errorValueFile)
{
	double eps = 1e-12;
	std::string outFile = errorValueFile.substr(0, errorValueFile.length() - 4);
	outFile.append("_texture.obj");
	std::vector<double> texture_distance;
	FILE *fin = fopen(errorValueFile.c_str(), "r");
	if (!fin)
	{
		printf("error: unable to open the input file: %s\n", errorValueFile.c_str());
		return -1;
	}
	texture_distance.clear();
	double relativeErrorValue;
	double texture_value = 1;
	double up_boundary = log10(10); //10%
	double down_boundary = log10(0.001); //0.001%
	double mx_relative_error_percent = 0;
	while (fscanf(fin, "%lf", &relativeErrorValue) != EOF)
	{
		if (relativeErrorValue * 100 > mx_relative_error_percent)
			mx_relative_error_percent = relativeErrorValue * 100;

		if (relativeErrorValue < eps)
			relativeErrorValue = eps;

		relativeErrorValue = log10(relativeErrorValue * 100);

		if (relativeErrorValue > up_boundary)
			relativeErrorValue = up_boundary;
		
		texture_value = (relativeErrorValue - down_boundary) / (up_boundary - down_boundary);

		if (texture_value < 0.01)
			texture_value = 0.01;
		if (texture_value > 0.99)
			texture_value = 0.99;

		texture_distance.push_back(texture_value);
	}

	fclose(fin);
	fin = fopen(scrData.c_str(), "r");
	if (!fin)
	{
		printf("error: unable to open the input file: %s\n", scrData.c_str());
		return -1;
	}
	FILE *fout = fopen(outFile.c_str(), "w");
	if (!fout)
	{
		printf("error: unable to open the output file: %s\n", outFile.c_str());
		return -1;
	}
	//int n_v, n_f;
	//fscanf(fin, "%d %d", &n_v, &n_f);
	//printf("n_v: %d, n_f: %d\n", n_v, n_f);
	double x1, y1, z1;
	int f1, f2, f3;
	char ch = ' ';
	int i = 0;
	char str[255] = { '\0' };
	while (fscanf(fin, "%c", &ch) != EOF)
	{
		if (ch == 'v')
		{
			fscanf(fin, "%lf %lf %lf", &x1, &y1, &z1);
			fprintf(fout, "v %.18lf %.18lf %.18lf\n", x1, y1, z1);
			fprintf(fout, "vt %.18lf %.18lf\n", texture_distance[i], texture_distance[i]);
			i++;
		}
		else if (ch == 'f')
		{
			fscanf(fin, "%d %d %d\n", &f1, &f2, &f3);
			fprintf(fout, "f %d/%d %d/%d %d/%d\n", f1, f1, f2, f2, f3, f3);
		}
		else if (ch == '\n')
		{
			continue;
		}
		else
		{
			fgets(str, 255, fin);
		}
	}
	fclose(fin);
	fclose(fout);

	return mx_relative_error_percent;
}

unsigned GeodesicResultProcess::FindNearestVertexIndex(geodesic::Point3D scrVert, std::string inFile)
{
	FILE *fin = fopen(inFile.c_str(), "r");
	if (!fin)
	{
		printf("error: unable to open the file: %s\n", inFile.c_str());
		return 0;
	}
	double x1, y1, z1;
	geodesic::Point3D vert_it;
	double minDistance = 1e100;
	double dis = 0;
	char ch = ' ';
	unsigned min_id = 0;
	char str[255] = { '\0' };
	unsigned i = 0;
	while (fscanf(fin, "%c", &ch) != EOF)
	{
		if (ch == 'v')
		{
			fscanf(fin, "%lf %lf %lf", &vert_it.x(), &vert_it.y(), &vert_it.z());
			dis = scrVert.distance(&vert_it);
			if (dis < minDistance)
			{
				minDistance = dis;
				min_id = i;
			}
			i++;
		}
		else if (ch == '\n')
		{
			continue;
		}
		else
		{
			fgets(str, 255, fin);
		}
	}
	fclose(fin);
	return min_id;
}


#endif