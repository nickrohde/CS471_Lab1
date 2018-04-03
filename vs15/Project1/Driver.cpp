#define _USE_MATH_DEFINES 1

#include <iostream>
#include <chrono>
#include <random>
#include <vector>
#include <math.h>
#include <fstream>

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1
#endif // !EXIT_SUCCESS

//#define M_PI 3.141592653589793238462643383279
//#define M_E  2.718281828459045235360287471352


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
vector<T>* getRandomVector(const std::size_t ui_SIZE, const T* p_MIN, const T* p_MAX)
{
	vector<T>* vec = new vector<T>(ui_SIZE);

	for (std::size_t i = 0; i < ui_SIZE; i++)
	{
		vec->at(i) = (getRandomNumberInRange(p_MIN, p_MAX));
	} // end for

	return vec;
} // end method getRandomVector

double schwefelsFunction(const vector<double>* vect)
{
	double total = 0.0;
	
	// SUM[1->n]
	for(auto& d: *vect)
	{
		double	temp = -1 * d,       // -x_i
				root = sqrt(abs(d)); // sqrt( | x_i | )
		temp *= sin(root);           // sin( sqrt( | x_i | ) )
		total += temp;               // sum
	} // end for
	
	return total;
} // end method schwefelsFunction

double firstDeJongsFunction(const vector<double>* vect)
{
	double total = 0.0;
	
	// SUM[1->n]
	for(auto& d: *vect)
	{
		total += d * d; // (x_i)^2
	} // end for
	
	return total;
} // end method firstDeJongsFunction

double rosenbrockFunction(const vector<double>* vect)
{
	double total = 0.0;
	
	// SUM[1->n-1]
	for(std::size_t i = 0; i < vect->size() - 1; i++)
	{
		double temp = 1 - vect->at(i),                            // (1 - x_i)
			   temp2 = vect->at(i) * vect->at(i) - vect->at(i+1); // (x_i^2 - x_i+1)

		temp *= temp;   // (1-x_i)^2

		temp2 *= temp2; // (x_i^2 - x_i+1)^2
		temp2 *= 100;   // 100(x_i^2 - x_i+1)^2
		
		total += temp2 + temp; // 100(x_i^2 - x_i+1)^2 + (1-x_i)^2
	} // end for
	
	return total;
} // end method rosenbrockFunction

double rastriginFunction(vector<double>* vect)
{
	double total = 0.0;
	
	// SUM[1->n]
	for(std::size_t i = 0; i < vect->size(); i++)
	{
		double temp = vect->at(i) * vect->at(i);      // x_i^2
		double temp2 = cos((M_PI * 2) * vect->at(i)); // cos(2pi * x_i)
		
		temp2 *= 10;             // 10cos(2pi * x_i)
		
		total += temp - temp2;   // x_i^2 - 10cos(2pi * x_i)
	} // end for
	
	total *= (10 * vect->size()); // 2n * SUM
	
	return total;
} // end method rastriginFunction


double griewangkFunction(vector<double>* vect)
{
	double total   = 0.0, 
		   sum     = 0.0, // SUM[1->n]
		   product = 1.0; // PROD[1->n]
	
	for(std::size_t i = 0; i < vect->size(); i++)
	{
		// SUM[1->n]
		{
			double tempSum = vect->at(i) * vect->at(i); // x_i^2

			tempSum /= 4000; // (x_i)^2 / 4000

			sum += tempSum;
		} // end SUM[1->n]

		// PROD[1->n]
		{
			double tempProd = vect->at(i) / sqrt(static_cast<double>(i)); // x_i / sqrt(i)

			tempProd = cos(tempProd); // cos( x_i / sqrt(i) )

			product *= tempProd;
		} // end PROD[1->n]
	} // end for
	
	total += 1 + sum - product; // 1 + SUM - PROD
	
	return total;
} // end method griewangkFunction


double sineEnvelopeSineWaveFunction(vector<double>* vect)
{
	double total   = 0.5 * (vect->size() - 1), // 0.5(n-1) + SUM
		   sum     = 0.0;
	
	// SUM[1->n-1]
	for(std::size_t i = 0; i < vect->size() - 1; i++)
	{
		double temp2 = 0,
			   quotient = 0,
			   sumOfSquares = 0;
		
		sumOfSquares = vect->at(i) * vect->at(i);      // x_i^2
		sumOfSquares += vect->at(i+1) * vect->at(i+1); // x_i^2 + (x_i+1)^2

		sum = sumOfSquares - 0.5;     // x_i^2 + (x_i+1)^2 - 0.5
		sum = sum * sum;              // ( x_i^2 + (x_i+1)^2 - 0.5 )^2
		sum = sin(sum);               // sin( x_i^2 + (x_i+1)^2 - 0.5 )^2
		
		temp2 = sumOfSquares * 0.001; //  0.001(x_i^2 + (x_i+1)^2)
		temp2 += 1;                   //  0.001(x_i^2 + (x_i+1)^2) + 1
		temp2 = temp2 * temp2;		  // (0.001(x_i^2 + (x_i+1)^2) + 1)^2
		
		quotient = sum/temp2;         // sin( x_i^2 + (x_i+1)^2 - 0.5 )^2 / (0.001(x_i^2 + (x_i+1)^2) + 1)^2
		
		total += quotient;
	} // end for
	
	total *= -1; // - SUM
	
	return total;
} // end method sineEnvelopeSineWaveFunction


double stretchedVSineWaveFunction(vector<double>* vect)
{
	double total   = 0.0;
	
	// SUM[1->n-1]
	for(std::size_t i = 0; i < vect->size() - 1; i++)
	{
		double sumOfSquares = 0.0,
			   product = 0.0,
			   temp = 0.0,
			   temp2 = 0.0;
		
		sumOfSquares = vect->at(i) * vect->at(i);          // x_i^2 
		sumOfSquares += vect->at(i + 1) * vect->at(i + 1); // x_i^2 + (x_i+1)^2

		temp = pow(sumOfSquares, 0.25); // root_4(x_i^2 + (x_i+1)^2)
		temp2 = pow(sumOfSquares, 0.1); // root_10(x_i^2 + (x_i+1)^2)
		
		temp2 *= 50.0;					// 50( root_10( x_i^2 + (x_i+1)^2 ))

		temp2 = temp2 * temp2;          // ( 50( root_10( x_i^2 + (x_i+1)^2 )))^2
		temp2 = sin(temp2);             // sin( 50( root_10( x_i^2 + (x_i+1)^2 )))^2

		product = temp2 * temp;         // sin( 50( root_10( x_i^2 + (x_i+1)^2 )))^2 * root_4(x_i^2 + (x_i+1)^2)

		total += product + 1;			// sin( 50( root_10( x_i^2 + (x_i+1)^2 )))^2 * root_4(x_i^2 + (x_i+1)^2) + 1 
	} // end for
	
	return total;
} // end method stretchedVSineWaveFunction


double ackleysOneFunction(vector<double>* vect)
{
	const double oneOverE = 1 / pow(M_E, 0.2); // 1 / e ^ 0.2

	double total = 0.0;

	// SUM[1->n-1]
	for (std::size_t i = 0; i < vect->size() - 1; i++)
	{
		double sumOfSquares = 0,
			   temp = 0,
			   temp2 = 0;

		sumOfSquares = vect->at(i) * vect->at(i);          // x_i^2 
		sumOfSquares += vect->at(i + 1) * vect->at(i + 1); // x_i^2 + (x_i+1)^2

		sumOfSquares = sqrt(sumOfSquares); // sqrt(x_i ^ 2 + (x_i + 1) ^ 2)
		sumOfSquares *= oneOverE;          // 1/(e^0.2) * sqrt(x_i ^ 2 + (x_i + 1) ^ 2)

		temp  = cos(2 * vect->at(i));      // cos(2 * x_i)
		temp2 = sin(2 * vect->at(i + 1));  // sin(2 * x_i+1)

		temp += temp2;                     // cos(2 * x_i) + sin(2 * x_i+1)
		temp *= 3;						   // 3(cos(2 * x_i) + sin(2 * x_i+1))

		total += temp + sumOfSquares;      // 1/(e^0.2) * sqrt(x_i ^ 2 + (x_i + 1) ^ 2) + 3(cos(2 * x_i) + sin(2 * x_i+1))
	} // end for

	return total;
} // end method ackleysOneFunction


double ackleysTwoFunction(vector<double>* vect)
{
	double total = 20.0 * (vect->size() - 1.0), // 20(n-1) + SUM
		   product = 0.0;

	// SUM[1->n-1]
	for (std::size_t i = 0; i < vect->size() - 1; i++)
	{
		double sumOfSquares = 0.0,
			   temp  = 0.0,
			   temp2 = 0.0;

		temp  = cos(2 * M_PI * vect->at(i));     // cos(2pi * x_i)
		temp2 = cos(2 * M_PI * vect->at(i + 1)); // cos(2pi * x_i+1)

		temp += temp2;                           // cos(2pi * x_i) + cos(2pi * x_i+1)
		temp *= 0.5;                             // (cos(2pi * x_i) + cos(2pi * x_i+1))/2

		temp2 = exp(temp);                       // e ^ ((cos(2pi * x_i) + cos(2pi * x_i+1))/2)

		sumOfSquares = vect->at(i) * vect->at(i);          // x_i^2 
		sumOfSquares += vect->at(i + 1) * vect->at(i + 1); // x_i^2 + (x_i+1)^2
		sumOfSquares /= 2;								   // (x_i^2 + (x_i+1)^2)/2

		sumOfSquares = sqrt(sumOfSquares);                 // sqrt((x_i^2 + (x_i+1)^2)/2)
		sumOfSquares *= 0.2;                               // (sqrt((x_i^2 + (x_i+1)^2)/2))/5

		temp = exp(sumOfSquares);                          // e ^ (sqrt((x_i^2 + (x_i+1)^2)/2))/5

		total += M_E - (20 / temp) - temp2;                // e - e ^ (sqrt((x_i^2 + (x_i+1)^2)/2))/5 - e ^ ((cos(2pi * x_i) + cos(2pi * x_i+1))/2)
	} // end for

	return total;
} // end method ackleysTwoFunction


double eggHolderFunction(vector<double>* vect)
{
	double total = 0,
		   product = 0.0;

	// SUM[1->n-1]
	for (std::size_t i = 0; i < vect->size() - 1; i++)
	{
		double temp = 0,
			   temp2 = 0;

		temp = vect->at(i) - vect->at(i + 1) - 47; // (x_i - x_(i+1) - 47)
		temp = abs(temp);                          // |(x_i - x_(i+1) - 47)|

		temp = sqrt(temp);                         // sqrt(|(x_i - x_(i+1) - 47)|)
		temp = sin(temp);                          // sin( sqrt(|(x_i - x_(i+1) - 47)|))
		temp *= vect->at(i);                       // sin( sqrt(|(x_i - x_(i+1) - 47)|)) * x_i

		temp2 = (vect->at(i) / 2) + vect->at(i + 1) + 47; // (x_i/2) + x_(i+1) + 47
		temp2 = abs(temp2);                               // |(x_i/2) + x_(i+1) + 47|
		
		temp2 = sqrt(temp2);                              // sqrt( |(x_i/2) + x_(i+1) + 47| )
		temp2 = sin(temp2);                               // sin( sqrt( |(x_i/2) + x_(i+1) + 47| ))
		temp2 *= (vect->at(i+1) + 47);                    // sin( sqrt( |(x_i/2) + x_(i+1) + 47| )) * (x_(i+1) + 47)

		total += temp + temp2;							  // (sin( sqrt(|(x_i - x_(i+1) - 47)|)) * x_i) + (sin( sqrt( |(x_i/2) + x_(i+1) + 47| )) * (x_(i+1) + 47))
	} // end for

	total *= -1; // - SUM

	return total;
} // end method eggHolderFunction


double ranaFunction(vector<double>* vect)
{
	double total = 20 * (vect->size() - 1),
		   product = 0.0;

	// SUM[1->n-1]
	for (std::size_t i = 0; i < vect->size() - 1; i++)
	{
		double	sum = 0,
				sumPlusXi = 0,
				sumMinXi = 0,
				temp = 0,
				temp2 = 0;

		sum = vect->at(i + 1) + 1;     // x_(i+1) + 1

		sumPlusXi = sum + vect->at(i); //  (x_(i+1) + 1) + x_i
		sumPlusXi = abs(sumPlusXi);    // |(x_(i+1) + 1) + x_i|

		sumMinXi = sum - vect->at(i);  //  (x_(i+1) + 1) - x_i
		sumMinXi = abs(sumMinXi);      // |(x_(i+1) + 1) - x_i|

		temp = vect->at(i) * sin(sumMinXi) * cos(sumPlusXi);          // x_i * sin(|(x_(i+1) + 1) - x_i|) * cos(|(x_(i+1) + 1) + x_i|)
		temp2 = (vect->at(i+1) + 1) * cos(sumMinXi) * sin(sumPlusXi); // (x_(i+1) + 1) * cos(|(x_(i+1) + 1) - x_i|) * sin(|(x_(i+1) + 1) + x_i|)

		total += temp + temp2; // [ x_i * sin(|(x_(i+1) + 1) - x_i|) * cos(|(x_(i+1) + 1) + x_i|)] + [ (x_(i+1) + 1) * cos(|(x_(i+1) + 1) - x_i|) * sin(|(x_(i+1) + 1) + x_i|)]
	} // end for

	return total;
} // end method ranaFunction


double pathologicalFunction(vector<double>* vect)
{
	double total = 0.5 * (vect->size() - 1);

	// SUM[1->n-1]
	for (std::size_t i = 0; i < vect->size() - 1; i++)
	{
		double	temp  = 0,
				temp2 = 0;

		// 100(x_i)^2 + x_(i+1)^2
		temp = 100 * vect->at(i) * vect->at(i) + vect->at(i + 1) * vect->at(i + 1);
		temp = sqrt(temp); // sqrt( 100(x_i)^2 + x_(i+1)^2)

		temp *= temp;      // (sqrt( 100(x_i)^2 + x_(i+1)^2))^2
		temp = sin(temp);  // sin(sqrt( 100(x_i)^2 + x_(i+1)^2))^2

		temp -= 0.5;       // sin(sqrt( 100(x_i)^2 + x_(i+1)^2))^2 - 0.5

		// (x_i)^2 - 2(x_i) * x_(i+1) + (x_(i+1))^2
		temp2 = vect->at(i) * vect->at(i) - 2 * vect->at(i) * vect->at(i+1) + vect->at(i + 1) * vect->at(i + 1);
		temp2 = temp2 * temp2; // ((x_i)^2 - 2(x_i) * x_(i+1) + (x_(i+1))^2)^2
		temp2 *= 0.001;        // 0.001((x_i)^2 - 2(x_i) * x_(i+1) + (x_(i+1))^2)^2
		temp2 += 1;            // 0.001((x_i)^2 - 2(x_i) * x_(i+1) + (x_(i+1))^2)^2 + 1

		total += temp / temp2; // [ sin^2(sqrt( 100(x_i)^2 + x_(i+1)^2))^2 - 0.5 ] / [ 0.001((x_i)^2 - 2(x_i) * x_(i+1) + (x_(i+1))^2)^2 + 1 ]
	} // end for

	return total;
} // end method pathologicalFunction


double michalewiczFunction(vector<double>* vect)
{
	double total = 0;

	// SUM[1->n]
	for (std::size_t i = 0; i < vect->size(); i++)
	{
		double	temp = 0,
				temp2 = 0;

		temp = sin(vect->at(i));               // sin(x_i)
		temp2 = i * vect->at(i) * vect->at(i); // i * (x_i)^2
		temp2 /= M_PI;                         // (i * (x_i)^2) / pi

		temp2 = pow(temp2, 20);                // ((i * (x_i)^2) / pi)^20
		temp2 = sin(temp2);                    // sin((i * (x_i)^2) / pi)^20
		
		total += temp * temp2;                 // sin(x_i) * sin((i * (x_i)^2) / pi)^20
	} // end for

	total *= -1; // - SUM

	return total;
} // end method michalewiczFunction


double masterCosineWaveFunction(vector<double>* vect)
{
	double total = 0;

	// SUM[1->n-1]
	for (std::size_t i = 0; i < vect->size() - 1; i++)
	{
		double	temp = 0,
				temp2 = 0;

		// x_i^2 + 0.5(x_i)(x_i+1) + (x_i+1)^2
		temp = (vect->at(i) * vect->at(i)) + (0.5 * vect->at(i) * vect->at(i + 1)) + (vect->at(i + 1) * vect->at(i + 1));

		temp2 = temp;

		temp *= -0.125;                  // (-1/8)(x_i^2 + 0.5(x_i)(x_i+1) + (x_i+1)^2)

		temp2 = 4 * sqrt(temp2);         // 4( x_i^2 + 0.5(x_i)(x_i+1) + (x_i+1)^2 )

		total += exp(temp) * cos(temp2); // e^((-1/8)(x_i^2 + 0.5(x_i)(x_i+1) + (x_i+1)^2)) * cos( 4( x_i^2 + 0.5(x_i)(x_i+1) + (x_i+1)^2 ))
	} // end for

	total *= -1; // - SUM

	return total;
} // end method masterCosineWaveFunction


double shekelsFoxholesFunction(const vector<double>* vect, double** da_A, const std::size_t ui_M = 30)
{
	const double da_C[] = { 0.806,0.517,0.1,0.908,0.965,0.669,0.524,0.902,0.351,0.876,0.462,
							0.491,0.463,0.741,0.352,0.869,0.813,0.811,0.0828,0.964,0.789,0.360,0.369,
							0.992,0.332,0.817,0.632,0.883,0.608,0.326 };

	double total = 0;

	// SUM[1->m]
	for (std::size_t i = 0; i < ui_M; i++)
	{
		double temp = 0;

		// SUM[1->n]
		for (std::size_t j = 0; j < vect->size(); j++)
		{
			double	temp2 = 0;

			temp2 = vect->at(j) - da_A[i][j]; // x_j - a_i,j
			temp += temp2 * temp2;			  // (x_j - a_i,j)^2
		} // end for

		temp += da_C[i];   // SUM[1->n](x_j - a_i,j)^2 + C_i

		total += 1 / temp; // 1 / ( SUM[1->n](x_j - a_i,j)^2 + C_i )
	} // end for

	total *= -1;
	
	return total;
} // end method shekelsFoxholesFunction


void makeMatrix(double**& da_A)
{
	da_A[0] = new double[10] {9.681, 0.667, 4.783, 9.095, 3.517, 9.325, 6.544, 0.211, 5.122, 2.02};
	da_A[1] = new double[10] { 9.4, 2.041, 3.788, 7.931, 2.882, 2.672, 3.568, 1.284, 7.033, 7.374 };
	da_A[2] = new double[10] { 8.025, 9.152, 5.114, 7.621, 4.564, 4.711, 2.996, 6.126, 0.734, 4.982 };
	da_A[3] = new double[10] { 2.196, 0.415, 5.649, 6.979, 9.510, 9.166, 6.304, 6.054, 9.377, 1.426 };
	da_A[4] = new double[10] { 8.074, 8.777, 3.467, 1.863, 6.708, 6.349, 4.534, 0.276, 7.633, 1.567 };
	da_A[5] = new double[10] { 7.650, 5.658, 0.720, 2.764, 3.278, 5.283, 7.474, 6.274, 1.409, 8.208 };
	da_A[6] = new double[10] { 1.256, 3.605, 8.623, 6.905, 4.584, 8.133, 6.071, 6.888, 4.187, 5.448 };
	da_A[7] = new double[10] { 8.314, 2.261, 4.24, 1.781, 4.124, 0.932, 8.129, 8.658, 1.208, 5.762 };
	da_A[8] = new double[10] { 0.226, 8.858, 1.42, 0.954, 1.622, 4.698, 6.228, 9.096, 0.972, 7.637 };
	da_A[9] = new double[10] { 7.305, 2.228, 1.242, 5.928, 9.133, 1.826, 4.06, 5.204, 8.713, 8.247 };
	da_A[10] = new double[10] { 0.652, 7.027, 0.508, 4.876, 8.807, 4.632, 5.808, 6.937, 3.291, 7.016 };
	da_A[11] = new double[10] { 2.699, 3.516, 5.847, 4.119, 4.461, 7.496, 8.817, 0.69, 6.593, 9.789 };
	da_A[12] = new double[10] { 8.327, 3.897, 2.017, 9.57, 9.825, 1.15, 1.395, 3.885, 6.354, 0.109 };
	da_A[13] = new double[10] { 2.132, 7.006, 7.136, 2.641, 1.882, 5.943, 7.273, 7.691, 2.88, 0.564 };
	da_A[14] = new double[10] { 4.707, 5.579, 4.08, 0.581, 9.698, 8.542, 8.077, 8.515, 9.231, 4.67 };
	da_A[15] = new double[10] { 8.304, 7.559, 8.567, 0.322, 7.128, 8.392, 1.472, 8.524, 2.277, 7.826 };
	da_A[16] = new double[10] { 8.632, 4.409, 4.832, 5.768, 7.05, 6.715, 1.711, 4.323, 4.405, 4.591 };
	da_A[17] = new double[10] { 4.887, 9.112, 0.17, 8.967, 9.693, 9.867, 7.508, 7.77, 8.382, 6.74 };
	da_A[18] = new double[10] { 2.44, 6.686, 4.299, 1.007, 7.008, 1.427, 9.398, 8.48, 9.95, 1.675 };
	da_A[19] = new double[10] { 6.306, 8.583, 6.084, 1.138, 4.350, 3.134, 7.853, 6.061, 7.457, 2.258 };
	da_A[20] = new double[10] { 0.652, 2.343, 1.37, 0.821, 1.31, 1.063, 0.689, 8.819, 8.833, 9.07 };
	da_A[21] = new double[10] { 5.558, 1.272, 5.756, 9.857, 2.279, 2.764, 1.284, 1.677, 1.244, 1.234 };
	da_A[22] = new double[10] { 3.352, 7.549, 9.817, 9.437, 8.687, 4.167, 2.57, 6.54, 0.228, 0.027 };
	da_A[23] = new double[10] { 8.798, 0.88, 2.37, 0.168, 1.701, 3.68, 1.231, 2.39, 2.499, 0.064 };
	da_A[24] = new double[10] { 1.46, 8.057, 1.337, 7.217, 7.914, 3.615, 9.981, 9.198, 5.292, 1.224 };
	da_A[25] = new double[10] { 0.432, 8.645, 8.774, 0.249, 8.081, 7.461, 4.416, 0.652, 4.002, 4.644 };
	da_A[26] = new double[10] { 0.679, 2.8, 5.523, 3.049, 2.968, 7.225, 6.73, 4.199, 9.614, 9.229 };
	da_A[27] = new double[10] { 4.263, 1.074, 7.286, 5.599, 8.291, 5.2, 9.214, 8.272, 4.398, 4.506 };
	da_A[28] = new double[10] { 9.496, 4.83, 3.15, 8.27, 5.079, 1.231, 5.731, 9.494, 1.883, 9.732 };
	da_A[29] = new double[10] { 4.138, 2.562, 2.532, 9.661, 5.611, 5.5, 6.886, 2.341, 9.699, 6.5 };
} // end method makeMatrix


int main(void)
{
	const std::size_t ui_SIZE = 10;
	const double d_MIN = 0.0;
	const double d_MAX = M_PI;
	
	double** da_A = new double*[30];

	makeMatrix(da_A);

	ofstream results("results.csv", ios::app | ios::out);

	for (int i = 0; i < 30; i++)
	{
		vector<double>* vec = getRandomVector(ui_SIZE, &d_MIN, &d_MAX);
		double temp = 0;

		temp = schwefelsFunction(vec);

		results << "F1: " << temp << ",";

		temp = firstDeJongsFunction(vec);
		results << "F2: " << temp << ",";
		temp = rosenbrockFunction(vec);
		results << "F3: " << temp << ",";
		temp = rastriginFunction(vec);
		results << "F4: " << temp << ",";
		temp = griewangkFunction(vec);
		results << "F5: " << temp << ",";
		temp = sineEnvelopeSineWaveFunction(vec);
		results << "F6: " << temp << ",";
		temp = stretchedVSineWaveFunction(vec);
		results << "F7: " << temp << ",";
		temp = ackleysOneFunction(vec);
		results << "F8: " << temp << ",";
		temp = ackleysTwoFunction(vec);
		results << "F9: " << temp << ",";
		temp = eggHolderFunction(vec);
		results << "F10: " << temp << ",";

		temp = ranaFunction(vec);
		results << "F11: " << temp << ",";
		temp = pathologicalFunction(vec);
		results << "F12: " << temp << ",";
		temp = michalewiczFunction(vec);
		results << "F13: " << temp << ",";
		temp = masterCosineWaveFunction(vec);
		results << "F14: " << temp << ",";
		temp = shekelsFoxholesFunction(vec, da_A);
		results << "F15: " << temp << ",\n";

		delete vec;
	} // end for


	for (std::size_t i = 0; i < 30; i++)
	{
		delete[] da_A[i];
	} // end for

	delete[] da_A;

	system("pause");

	return EXIT_SUCCESS;
} // end method Main                                                              