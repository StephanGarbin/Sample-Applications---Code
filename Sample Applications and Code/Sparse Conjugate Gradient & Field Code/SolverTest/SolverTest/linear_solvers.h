#pragma once

#ifndef DROGON_LINEAR_SOLVERS_H
#define DROGON_LINEAR_SOLVERS_H

//STL

//TBB
#include "tbb\tick_count.h"

//DROGON
#include "drogon_status.h"
#include "dBLAS_1.h"
#include "dBLAS_CG.h"
#include "field_sparse.h"
#include "field_dense.h"
#include "dMatrix_Laplacian7.h"
#include "dVector.h"

namespace drogon {

	template<class SOLVER, class SYSTEM>
	class linear_solvers
	{
	public:
		linear_solvers()
		{
		}
		~linear_solvers()
		{
		}

		static drogon_status<SOLVER> conjugate_gradient(dMatrix_Laplacian7<SYSTEM>& A, field_dense<SOLVER>& x, field_sparse<SYSTEM>& b,
			SOLVER tolerance, int max_iterations)
		{
			//1. PRELIMINARIES
			//start the timer
			tick_count start = tick_count::now();

			//number of iterations taken
			size_t it = 0;

			//Initialisation of temporary vectoes
			size_t size = x.get_size();
			size_t width = x.get_dataWindow().x();
			size_t height = x.get_dataWindow().y();
			size_t depth = x.get_dataWindow().z();

			//residual
			field_dense<SOLVER> residual("residual", x.get_dataWindow_EXT(), x.get_displayWindow_EXT(), x.get_simulationWindow_EXT());

			//search direction
			field_dense<SOLVER> d("d", x.get_dataWindow_EXT(), x.get_displayWindow_EXT(), x.get_simulationWindow_EXT());

			//array to hold A * d, which is used multiple times
			field_dense<SOLVER> Ad("Ad", x.get_dataWindow_EXT(), x.get_displayWindow_EXT(), x.get_simulationWindow_EXT());

			SOLVER alpha = 0;
			SOLVER beta = 0;
			SOLVER rho_old = 0;
			SOLVER rho_new = 0;

			dVector<SOLVER, 2> temp;

			//2. FIRST STEP
			// calculate residual & direction: r = b - Ax
			dBLAS_CG<SOLVER, SYSTEM>::calculate_residual_sparse_PAR(residual, b, x, A);

			//set first search direction
			d.copy_from(residual, true);

			//set rho_old = rTr
			rho_new = dBLAS_1<SOLVER>::dot_product_PAR(residual, residual);
			rho_old = rho_new;

			//3. ITERATE
			for (it = 1; it < max_iterations; ++it)
			{
				//Calculate Ad & Ad *d
				temp[0] = dBLAS_CG<SOLVER, SYSTEM>::matrix_vector__dot_product_PAR(Ad, d, A);

				//Calculate alpha: rho_new / d SOLVER A*d
				alpha = rho_new / temp[0];

				temp = dBLAS_CG<SOLVER, SYSTEM>::add_subtr_scale_infnorm_PAR(residual, x, d, Ad, alpha);

				//check for convergence
				if (temp[1] <= tolerance)
				{
					//end the timer
					tick_count end = tick_count::now();
					return drogon_status<SOLVER>(true, start, end, it, temp[1], tolerance);
				}

				//Calculate / assign rho: rho_old = rho_new | rho_new = rTr
				rho_old = rho_new;
				rho_new = temp[0];

				//Calculate beta: beta = rho_new / rho_old
				beta = rho_new / rho_old;

				//Calculate new search direction: d = r + beta * d
				dBLAS_CG<SOLVER, SYSTEM>::add_scale_PAR(d, residual, beta);
			}
			//no solution found using the given number of iterations :(

			//We might as well recalculate the residual for debug purposes here
			temp[1] = dBLAS_1<SOLVER>::inf_norm_PAR(residual);

			//end the timer
			tick_count end = tick_count::now();
			return drogon_status<SOLVER>(false, start, end, it, temp[1], tolerance);
		}

		static drogon_status<SOLVER> conjugate_gradient_SEQ(dMatrix_Laplacian7<SYSTEM>& A, field_dense<SOLVER>& x, field_sparse<SYSTEM>& b,
			SOLVER tolerance, int max_iterations)
		{
			//1. PRELIMINARIES
			//start the timer
			tick_count start = tick_count::now();

			//number of iterations taken
			size_t it = 0;

			//Initialisation of temporary vectoes
			size_t size = x.get_size();
			size_t width = x.get_dataWindow().x();
			size_t height = x.get_dataWindow().y();
			size_t depth = x.get_dataWindow().z();

			//residual
			field_dense<SOLVER> residual("residual", x.get_dataWindow_EXT(), x.get_displayWindow_EXT(), x.get_simulationWindow_EXT());

			//search direction
			field_dense<SOLVER> d("d", x.get_dataWindow_EXT(), x.get_displayWindow_EXT(), x.get_simulationWindow_EXT());

			//array to hold A * d, which is used multiple times
			field_dense<SOLVER> Ad("Ad", x.get_dataWindow_EXT(), x.get_displayWindow_EXT(), x.get_simulationWindow_EXT());

			SOLVER alpha = 0;
			SOLVER beta = 0;
			SOLVER rho_old = 0;
			SOLVER rho_new = 0;

			dVector<SOLVER, 2> temp;

			//2. FIRST STEP
			// calculate residual & direction: r = b - Ax
			dBLAS_CG<SOLVER, SYSTEM>::calculate_residual_sparse_SEQ(residual, b, x, A);

			//set first search direction
			d.copy_from(residual, true);

			//set rho_old = rTr
			rho_new = dBLAS_1<SOLVER>::dot_product_SEQ(residual, residual);
			rho_old = rho_new;

			//3. ITERATE
			for (it = 1; it < max_iterations; ++it)
			{
				//Calculate Ad & Ad *d
				temp[0] = dBLAS_CG<SOLVER, SYSTEM>::matrix_vector__dot_product_SEQ(Ad, d, A);

				//Calculate alpha: rho_new / d SOLVER A*d
				alpha = rho_new / temp[0];

				temp = dBLAS_CG<SOLVER, SYSTEM>::add_subtr_scale_infnorm_SEQ(residual, x, d, Ad, alpha);

				//check for convergence
				if (temp[1] <= tolerance)
				{
					//end the timer
					tick_count end = tick_count::now();
					return drogon_status<SOLVER>(true, start, end, it, temp[1], tolerance);
				}

				//Calculate / assign rho: rho_old = rho_new | rho_new = rTr
				rho_old = rho_new;
				rho_new = temp[0];

				//Calculate beta: beta = rho_new / rho_old
				beta = rho_new / rho_old;

				//Calculate new search direction: d = r + beta * d
				dBLAS_CG<SOLVER, SYSTEM>::add_scale_SEQ(d, residual, beta);
			}
			//no solution found using the given number of iterations :(

			//We might as well recalculate the residual for debug purposes here
			temp[1] = dBLAS_1<SOLVER>::inf_norm_SEQ(residual);

			//free memory
			//delete residual;
			//delete d;
			//delete Ad;

			//end the timer
			tick_count end = tick_count::now();
			return drogon_status<SOLVER>(false, start, end, it, temp[1], tolerance);
		}

		static drogon_status<SOLVER> conjugate_gradient_DENSE(dMatrix_Laplacian7<SYSTEM>& A, field_dense<SOLVER>& x, field_dense<SYSTEM>& b,
			SOLVER tolerance, int max_iterations)
		{
			//1. PRELIMINARIES
			//start the timer
			tick_count start = tick_count::now();

			//number of iterations taken
			size_t it = 0;

			//Initialisation of temporary vectoes
			size_t size = x.get_size();
			size_t width = x.get_dataWindow().x();
			size_t height = x.get_dataWindow().y();
			size_t depth = x.get_dataWindow().z();

			//residual
			field_dense<SOLVER> residual("residual", x.get_dataWindow_EXT(), x.get_displayWindow_EXT(), x.get_simulationWindow_EXT());

			//search direction
			field_dense<SOLVER> d("d", x.get_dataWindow_EXT(), x.get_displayWindow_EXT(), x.get_simulationWindow_EXT());

			//array to hold A * d, which is used multiple times
			field_dense<SOLVER> Ad("Ad", x.get_dataWindow_EXT(), x.get_displayWindow_EXT(), x.get_simulationWindow_EXT());

			SOLVER alpha = 0;
			SOLVER beta = 0;
			SOLVER rho_old = 0;
			SOLVER rho_new = 0;

			dVector<SOLVER, 2> temp;

			//2. FIRST STEP
			// calculate residual & direction: r = b - Ax
			dBLAS_CG<SOLVER, SYSTEM>::calculate_residual_PAR(residual, b, x, A);

			//set first search direction
			d.copy_from(residual, true);

			//set rho_old = rTr
			rho_new = dBLAS_1<SOLVER>::dot_product_PAR(residual, residual);
			rho_old = rho_new;

			//3. ITERATE
			for (it = 1; it < max_iterations; ++it)
			{
				//Calculate Ad & Ad *d
				temp[0] = dBLAS_CG<SOLVER, SYSTEM>::matrix_vector__dot_product_PAR(Ad, d, A);

				//Calculate alpha: rho_new / d SOLVER A*d
				alpha = rho_new / temp[0];

				temp = dBLAS_CG<SOLVER, SYSTEM>::add_subtr_scale_infnorm_PAR(residual, x, d, Ad, alpha);

				//check for convergence
				if (temp[1] <= tolerance)
				{
					//end the timer
					tick_count end = tick_count::now();
					return drogon_status<SOLVER>(true, start, end, it, temp[1], tolerance);
				}

				//Calculate / assign rho: rho_old = rho_new | rho_new = rTr
				rho_old = rho_new;
				rho_new = temp[0];

				//Calculate beta: beta = rho_new / rho_old
				beta = rho_new / rho_old;

				//Calculate new search direction: d = r + beta * d
				dBLAS_CG<SOLVER, SYSTEM>::add_scale_PAR(d, residual, beta);
			}
			//no solution found using the given number of iterations :(

			//We might as well recalculate the residual for debug purposes here
			temp[1] = dBLAS_1<SOLVER>::inf_norm_PAR(residual);

			//end the timer
			tick_count end = tick_count::now();
			return drogon_status<SOLVER>(false, start, end, it, temp[1], tolerance);
		}

		static drogon_status<SOLVER> conjugate_gradient_DENSE_2D(dMatrix_Laplacian5<SYSTEM>& A, field_dense_2D<SOLVER>& x, field_dense_2D<SYSTEM>& b,
			SOLVER tolerance, int max_iterations)
		{
			//1. PRELIMINARIES
			//start the timer
			tick_count start = tick_count::now();

			//number of iterations taken
			size_t it = 0;

			//Initialisation of temporary vectoes
			size_t size = x.get_size();
			size_t width = x.get_dataWindow().x();
			size_t height = x.get_dataWindow().y();

			//residual
			field_dense_2D<SOLVER> residual("residual", x.get_dataWindow_EXT(), x.get_displayWindow_EXT(), x.get_simulationWindow_EXT());

			//search direction
			field_dense_2D<SOLVER> d("d", x.get_dataWindow_EXT(), x.get_displayWindow_EXT(), x.get_simulationWindow_EXT());

			//array to hold A * d, which is used multiple times
			field_dense_2D<SOLVER> Ad("Ad", x.get_dataWindow_EXT(), x.get_displayWindow_EXT(), x.get_simulationWindow_EXT());

			SOLVER alpha = 0;
			SOLVER beta = 0;
			SOLVER rho_old = 0;
			SOLVER rho_new = 0;

			dVector<SOLVER, 2> temp;

			//2. FIRST STEP
			// calculate residual & direction: r = b - Ax
			dBLAS_CG<SOLVER, SYSTEM>::calculate_residual_PAR_2D(residual, b, x, A);

			//set first search direction
			d.copy_from(residual, true);

			//set rho_old = rTr
			rho_new = dBLAS_1<SOLVER>::dot_product_PAR_2D(residual, residual);
			rho_old = rho_new;

			//3. ITERATE
			for (it = 1; it < max_iterations; ++it)
			{
				//Calculate Ad & Ad *d
				temp[0] = dBLAS_CG<SOLVER, SYSTEM>::matrix_vector__dot_product_PAR_2D(Ad, d, A);

				//Calculate alpha: rho_new / d SOLVER A*d
				alpha = rho_new / temp[0];

				temp = dBLAS_CG<SOLVER, SYSTEM>::add_subtr_scale_infnorm_PAR_2D(residual, x, d, Ad, alpha);

				//check for convergence
				if (temp[1] <= tolerance)
				{
					//end the timer
					tick_count end = tick_count::now();
					return drogon_status<SOLVER>(true, start, end, it, temp[1], tolerance);
				}

				//Calculate / assign rho: rho_old = rho_new | rho_new = rTr
				rho_old = rho_new;
				rho_new = temp[0];

				//Calculate beta: beta = rho_new / rho_old
				beta = rho_new / rho_old;

				//Calculate new search direction: d = r + beta * d
				dBLAS_CG<SOLVER, SYSTEM>::add_scale_PAR_2D(d, residual, beta);
			}
			//no solution found using the given number of iterations :(

			//We might as well recalculate the residual for debug purposes here
			temp[1] = dBLAS_1<SOLVER>::inf_norm_PAR_2D(residual);

			//end the timer
			tick_count end = tick_count::now();
			return drogon_status<SOLVER>(false, start, end, it, temp[1], tolerance);
		}

	private:

	};

} /*namespace drogon*/

#endif /*DROGON_LINEAR_SOLVERS_H*/