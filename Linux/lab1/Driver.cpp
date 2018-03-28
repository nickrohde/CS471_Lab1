#include <iostream>
#include <chrono>
#include <random>
#include <vector>
#include <math.h>
#include <cmath>

#define M_PI       3.14159265358979323846   // pi -- stolen from corecrt_math_defines.h

using namespace std;


template <typename T> 
T getRandomNumberInRange(const T* p_MIN, const T* p_MAX)
{
	static std::random_device rd{};
	static std::mt19937 engine{ rd() };
	static std::uniform_real_distribution<T> dist{ *p_MIN, *p_MAX };

	return dist(engine);
} // end method getRandomNumberInRange


template <typename T>
vector<T>* getRandomVector(const size_t ui_SIZE, const T* p_MIN, const T* p_MAX)
{
	vector<T>* vec = new vector<T>(ui_SIZE);

	for (size_t i = 0; i < ui_SIZE; i++)
	{
		vec->at(i) = (getRandomNumberInRange(p_MIN, p_MAX));
	} // end for

	return vec;
} // end method getRandomVector

double schwefelsFunction(vector<double>* vect)
{
	double total = 0.0;
	
	for(auto& d: *vect)
	{
		double temp = -1 * d;
		double root = sqrt(abs(d));
		temp *= sin(root);
		total += temp;
	} // end for
	
	return total;
} // end method schwefelsFunction

double firstDeJongsFunction(vector<double>* vect)
{
	double total = 0.0;
	
	for(auto& d: *vect)
	{
		double temp = d * d;
		total += temp;
	} // end for
	
	return total;
} // end method firstDeJongsFunction

double rosenbrockFunction(vector<double>* vect)
{
	double total = 0.0;
	
	for(size_t i = 0; i < vect->size() - 1; i++)
	{
		double temp = 1 - vect->at(i),
			   temp2 = vect->at(i) * vect->at(i) - vect->at(i+1);
			   temp2 *= temp2;
		temp *= temp;
		temp2 *= 100;
		
		total += temp2 + temp;
	} // end for
	
	return total;
} // end method rosenbrockFunction

double rastriginFunction(vector<double>* vect)
{
	double total = 0.0;
	
	for(size_t i = 0; i < vect->size(); i++)
	{
		double temp = vect->at(i) * vect->at(i);
		double temp2 = cos((M_PI * 2) * vect->at(i));
		
		temp2 *= 10;
		temp -= temp2;
	} // end for
	
	total *= (2 * vect->size());
	
	return total;
} // end method rastriginFunction


int main(void)
{
	const size_t ui_SIZE = 10;
	const double d_MIN = 0.0;
	const double d_MAX = 10000.0;

	for (int i = 0; i < 30; i++)
	{
		vector<double>* vec = getRandomVector(ui_SIZE, &d_MIN, &d_MAX);


		delete vec;
	} // end for

	system("pause");

	return EXIT_SUCCESS;
} // end method Main                                                              