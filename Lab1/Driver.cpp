#include <iostream>
#include <chrono>
#include <random>
#include <vector>


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