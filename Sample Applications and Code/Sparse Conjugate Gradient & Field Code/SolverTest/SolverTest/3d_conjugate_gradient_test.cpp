#include <iostream>

#include "field_dense.h"
#include "field_sparse.h"
#include "dMatrix_Laplacian7.h"
#include "linear_solvers.h"
#include "drogon_status.h"

using namespace std;
using namespace drogon;


int main(int argc, char* argv[])
{
	//cout << "Conjugate Gradient Test Starting ... " << endl << endl;

	//size_t i;
	//size_t j;
	//size_t k;
	//size_t max_its;

	//cout << "Please enter field dimension width (use only mutiples of two, e.g. 32)" << endl;
	//cin >> i;

	//while(i%2 != 0)
	//{
	//	cout << "input is not divisible by two, please enter another value..." << endl;
	//	cin >> i;
	//}

	//cout << endl;
	//cout << "Please enter field dimension height (use only mutiples of two, e.g. 32)" << endl;
	//cin >> j;

	//while(j%2 != 0)
	//{
	//	cout << "input is not divisible by two, please enter another value..." << endl;
	//	cin >> j;
	//}

	//cout << endl;
	//cout << "Please enter field dimension depth (use only mutiples of two, e.g. 32)" << endl;
	//cin >> k;

	//while(k%2 != 0)
	//{
	//	cout << "input is not divisible by two, please enter another value..." << endl;
	//	cin >> k;
	//}

	//cout << endl;
	//cout << "Please enter maximum number of cg iterations" << endl;
	//cin >> max_its;

	//cout << endl;
	//cout << "Tile-dimension set automatically ..." << endl;
	//cout << endl;
	//cout << "Allocating memory ... " << endl;

	size_t i = 32;
	size_t j = 32;
	size_t k = 32;
	size_t max_its = 10000;

	size_t tile_dim = 4;

	FIELD_EXTENT data_window;
	data_window.minX() = 0; data_window.minY() = 0; data_window.minZ() = 0; data_window.maxX() = i; data_window.maxY() = j; data_window.maxZ() = k;

	FIELD_EXTENT simulation_window;
	simulation_window.minX() = 0; simulation_window.minY() = 0; simulation_window.minZ() = 0; simulation_window.maxX() = i; simulation_window.maxY() = j; simulation_window.maxZ() = k;

	FIELD_EXTENT display_window;
	display_window.minX() = 0; display_window.minY() = 0; display_window.minZ() = 0; display_window.maxX() = i; display_window.maxY() = j; display_window.maxZ() = k;

	field_dense<float> dense_1("dense_1", data_window, display_window, simulation_window);
	field_sparse<float> sparse_1("sparse_1", data_window, display_window, simulation_window, tile_dim);

	field_dense<double> result_2d("result_1", data_window, display_window, simulation_window);

	dMatrix_Laplacian7<float> matrix_1("matrix_1", data_window, i * j * k, 0);

	dense_1.set_field(2.58);
	sparse_1.set_field(1.0);
	result_2d.set_field(0.0);

	matrix_1.set_Diag(5.0);
	matrix_1.set_I(1.1);
	matrix_1.set_J(1.1);
	matrix_1.set_K(1.1);
	matrix_1.check_coefficients();

	cout << "Testing linear solvers ..." << endl;
	cout << endl << endl;

	linear_solvers<double, float> solver;

	drogon_status<double> stat_1 = solver.conjugate_gradient(matrix_1, result_2d, sparse_1, 1e-6, max_its);
	stat_1.print_cg_results();

	result_2d.set_field(0.0);
	drogon_status<double> stat_2 = solver.conjugate_gradient_SEQ(matrix_1, result_2d, sparse_1, 1e-6, max_its);
	stat_2.print_cg_results();

	result_2d.set_field(0.0);
	drogon_status<double> stat_3 = solver.conjugate_gradient_DENSE(matrix_1, result_2d, dense_1, 1e-6, max_its);
	stat_3.print_cg_results();

	cout << endl << endl;
	cout << "Test completed..." << endl;
	int _end = 0;
	cin >> _end;

	return 0;
}