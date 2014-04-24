#pragma once

#ifndef DROGON_DBLAS_CG_H
#define DROGON_DBLAS_CG_H

//STL

//TBB
#include "tbb\parallel_for.h"
#include "tbb\parallel_reduce.h"
#include "tbb\blocked_range.h"
#include "tbb\blocked_range2d.h"
#include "tbb\blocked_range3d.h"

//DROGON
#include "dVector.h"
#include "field_dense.h"
#include "field_dense_2D.h"
#include "field_sparse.h"
#include "dMatrix_Laplacian7.h"
#include "dMatrix_Laplacian5.h"

using namespace tbb;

namespace drogon {

	template<class SOLVER, class SYSTEM>
	class dBLAS_CG
	{
	public:
		dBLAS_CG() {}
		~dBLAS_CG() {}

		//! Calculates residual for dense field systems (SEQUENTIAL)
		static void calculate_residual_SEQ(field_dense<SOLVER>& residual, field_dense<SOLVER>& b, field_dense<SOLVER>& x, dMatrix_Laplacian7<SYSTEM>& A)
		{
			size_t slice = A.get_system_width() * A.get_system_height();
			size_t width = A.get_system_width();

			//put iterations this way to avoid cache misses
			for (size_t k = 0; k < A.get_system_depth(); ++k)
			{
				for (size_t j = 0; j < A.get_system_height(); ++j)
				{
					for (size_t i = 0; i < width; ++i)
					{
						size_t index = I_2(i, j, k, width, slice);
						residual.set(index, (SOLVER)b.get(index) - matrix_vector_product_row_inner(A, x, i, j, k, index));
					}
				}
			}
		}

		//! Calculates residual for dense field systems (PARALLEL)
		static void calculate_residual_PAR(field_dense<SOLVER>& residual, field_dense<SYSTEM>& b, field_dense<SOLVER>& x, dMatrix_Laplacian7<SYSTEM>& A)
		{
			size_t depth = A.get_system_depth();
			size_t height = A.get_system_height();
			size_t width = A.get_system_width();
			size_t slice = width * height;

			parallel_for(blocked_range3d<size_t, size_t, size_t>(0, depth, 0, height, 0, width),
				[&](blocked_range3d<size_t, size_t, size_t> &r) {
				//make local copies to tell compiler to hold these in registers
				size_t l_slice = slice;
				size_t l_width = width;

				for (size_t k = r.pages().begin(); k != r.pages().end(); ++k)
				{
					for (size_t j = r.rows().begin(); j != r.rows().end(); ++j)
					{
						for (size_t i = r.cols().begin(); i != r.cols().end(); ++i)
						{
							size_t index = I_2(i, j, k, l_width, l_slice);
							residual.set(index, (SOLVER)b.get(index) - matrix_vector_product_row_inner(A, x, i, j, k, index));
						}
					}
				}
			});
		}

		//! Calculates residual for dense field systems (PARALLEL) (2D Problems)
		static void calculate_residual_PAR_2D(field_dense_2D<SOLVER>& residual, field_dense_2D<SYSTEM>& b, field_dense_2D<SOLVER>& x, dMatrix_Laplacian5<SYSTEM>& A)
		{
			size_t height = A.get_system_height();
			size_t width = A.get_system_width();

			parallel_for(blocked_range2d<size_t, size_t>(0, height, 0, width),
				[&](blocked_range2d<size_t, size_t> &r) {
				//make local copies to tell compiler to hold these in registers
				size_t l_width = width;

				for (size_t j = r.rows().begin(); j != r.rows().end(); ++j)
				{
					for (size_t i = r.cols().begin(); i != r.cols().end(); ++i)
					{
						size_t index = I_2_2D(i, j, l_width);

						residual.set(index, (SOLVER)b.get(index) - matrix_vector_product_row_inner_2D(A, x, i, j, index));
					}
				}
			});
		}

		//! Calculates residual for sparse field systems (SEQUENTIAL)
		static void calculate_residual_sparse_SEQ(field_dense<SOLVER>& residual, field_sparse<SYSTEM>& b, field_dense<SOLVER>& x, dMatrix_Laplacian7<SYSTEM>& A)
		{
			size_t tile_dim = b.get_tile_edgeLength();
			size_t slice = A.get_system_width() * A.get_system_height();
			size_t width = A.get_system_width();

			size_t system_tiles_slice = b.get_tiles_window().x() * b.get_tiles_window().y();
			size_t system_tiles_width = b.get_tiles_window().x();

			//put iterations this way to avoid cache misses
			for (size_t t_k = 0; t_k < b.get_tiles_window().z(); ++t_k)
			{
				for (size_t t_j = 0; t_j < b.get_tiles_window().y(); ++t_j)
				{
					for (size_t t_i = 0; t_i < b.get_tiles_window().x(); ++t_i)
					{
						size_t tile_index = I_2(t_i, t_j, t_k, system_tiles_width, system_tiles_slice);
						if (b.is_active(tile_index))
						{
							//inner loop
							for (size_t k = 0; k < tile_dim; ++k)
							{
								for (size_t j = 0; j < tile_dim; ++j)
								{
									for (size_t i = 0; i < tile_dim; ++i)
									{
										size_t global_i = t_i * tile_dim + i;
										size_t global_j = t_j * tile_dim + j;
										size_t global_k = t_k * tile_dim + k;

										size_t index = I_2(global_i, global_j, global_k, width, slice);
										/*residual.set(index, (SOLVER)b.get_fast(tile_index, i, j, k) - matrix_vector_product_row_inner(A, x, global_i, global_j, global_k, index));*/
										residual.set(index, (SOLVER)b.get_fast(tile_index, I_2(i, j, k, tile_dim, tile_dim)) - matrix_vector_product_row_inner(A, x, global_i, global_j, global_k, index));

										//cout << index << "\n";
									}
								}
							}
						}
						//else // VECTORISATION
						//{
						//	size_t global_i = 0;
						//	size_t global_j = 0;
						//	size_t global_k = 0;

						//	size_t index = 0;
						//	//residual.set(index, (SOLVER)b.get_fast(tile_index, I_2(i, j, k, tile_dim, tile_dim)) - matrix_vector_product_row_inner(A, x, global_i, global_j, global_k, index));
						//}
					}
				}
			}
		}

		//improved cache efficiency at the expense of more operations
		//do not use this other than for experimentative purposes!
		static void calculate_residual_sparse_SEQ_2(field_dense<SOLVER>& residual, field_sparse<SYSTEM>& b, field_dense<SOLVER>& x, dMatrix_Laplacian7<SYSTEM>& A)
		{
			size_t tile_dim = b.get_tile_edgeLength();
			size_t slice = A.get_system_width() * A.get_system_height();
			size_t width = A.get_system_width();

			//put iterations this way to avoid cache misses
			for (size_t k = 0; k < A.get_system_depth(); ++k)
			{
				for (size_t j = 0; j < A.get_system_height(); ++j)
				{
					for (size_t i = 0; i < width; ++i)
					{
						size_t t_i = i / tile_dim;
						size_t t_j = j / tile_dim;
						size_t t_k = k / tile_dim;
						size_t tile_index = I_2(t_i, t_j, t_k, width, slice);
						if (b.is_active(tile_index))
						{
							//WORK IN PROGRESS
							size_t index = 0;
							//size_t index = I_2(i, j, k, width, tile_area);
							residual.set(index, (SOLVER)b.get_fast(tile_index, I_2(t_i%tile_dim, t_j%tile_dim, t_k%tile_dim, tile_dim, tile_dim)) - matrix_vector_product_row_inner(A, x, i, j, k, index));
						}
					}
				}
			}
		}

		//! Calculates residual for sparse field systems (PARALLEL)
		static void calculate_residual_sparse_PAR(field_dense<SOLVER>& residual, field_sparse<SYSTEM>& b, field_dense<SOLVER>& x, dMatrix_Laplacian7<SYSTEM>& A)
		{
			size_t width = b.get_dataWindow().x(); // correct accessors!
			size_t tile_dim = b.get_tile_edgeLength();

			size_t system_tiles_slice = b.get_tiles_window().x() * b.get_tiles_window().y();
			size_t system_tiles_width = b.get_tiles_window().x();


			size_t slice = width * b.get_dataWindow().y();

			parallel_for(blocked_range3d<size_t, size_t, size_t>(0, b.get_tiles_window().z(), 0, b.get_tiles_window().y(), 0, b.get_tiles_window().x()),
				[&](blocked_range3d<size_t, size_t, size_t> &r) {
				//make local copies to tell compiler to hold these in registers
				size_t l_width = width;
				size_t l_slice = slice;

				size_t l_system_tiles_width = system_tiles_width;
				size_t l_system_tiles_slice = system_tiles_slice;

				for (size_t t_k = r.pages().begin(); t_k != r.pages().end(); ++t_k)
				{
					for (size_t t_j = r.rows().begin(); t_j != r.rows().end(); ++t_j)
					{
						for (size_t t_i = r.cols().begin(); t_i != r.cols().end(); ++t_i)
						{
							size_t tile_index = I_2(t_i, t_j, t_k, l_system_tiles_width, l_system_tiles_slice);
							if (b.is_active(tile_index))
							{
								//inner loop
								for (size_t k = 0; k < tile_dim; ++k)
								{
									for (size_t j = 0; j < tile_dim; ++j)
									{
										for (size_t i = 0; i < tile_dim; ++i)
										{
											size_t global_i = t_i * tile_dim + i;
											size_t global_j = t_j * tile_dim + j;
											size_t global_k = t_k * tile_dim + k;

											size_t index = I_2(global_i, global_j, global_k, l_width, l_slice);
											residual.set(index, (SOLVER)b.get_fast(tile_index, I_2(i, j, k, tile_dim, tile_dim)) - matrix_vector_product_row_inner(A, x, global_i, global_j, global_k, index));
											//cout << index << "\n";
										}
									}
								}
							}
						}
					}
				}
			});
		}

		//! Combines Matrix-vector & dot products for greater efficiency (SEQUENTIAL)
		static SOLVER matrix_vector__dot_product_SEQ(field_dense<SOLVER>& result, field_dense<SOLVER>& x, dMatrix_Laplacian7<SYSTEM>& A)
		{
			SOLVER sum = 0;
			size_t slice = A.get_system_width() * A.get_system_height();
			size_t width = A.get_system_width();

			//put iterations this way to avoid cache misses
			for (size_t k = 0; k < A.get_system_depth(); ++k)
			{
				for (size_t j = 0; j < A.get_system_height(); ++j)
				{
					for (size_t i = 0; i < width; ++i)
					{
						size_t index = I_2(i, j, k, width, slice);
						//matrix-vector op
						result.set(index, matrix_vector_product_row_inner(A, x, i, j, k, index));
						//dot product op
						sum += result.get(index) * x.get(index);
					}
				}
			}
			return sum;
		}

		//! Combines Matrix-vector & dot products for greater efficiency (PARALLEL)
		static SOLVER matrix_vector__dot_product_PAR(field_dense<SOLVER>& result, field_dense<SOLVER>& x, dMatrix_Laplacian7<SYSTEM>& A)
		{
			size_t depth = x.get_dataWindow().z();
			size_t height = x.get_dataWindow().y();
			size_t width = x.get_dataWindow().x();
			size_t slice = width * height;

			return parallel_reduce(blocked_range3d<size_t, size_t, size_t>(0, depth, 0, height, 0, width), SOLVER(0),
				[&](const blocked_range3d<size_t, size_t, size_t> &r, SOLVER init)->SOLVER {
				//make local copies to tell compiler to hold these in registers
				SOLVER local_init = init;
				size_t l_slice = slice;
				size_t l_width = width;
				for (size_t k = r.pages().begin(); k != r.pages().end(); ++k)
				{
					for (size_t j = r.rows().begin(); j != r.rows().end(); ++j)
					{
						for (size_t i = r.cols().begin(); i != r.cols().end(); ++i)
						{
							size_t index = I_2(i, j, k, l_width, l_slice);
							result.set(index, matrix_vector_product_row_inner(A, x, i, j, k, index));
							local_init += result.get(index) * x.get(index);
						}
					}
				}
				return local_init;
			}, [](SOLVER x, SOLVER y)->SOLVER {
				return x + y;
			});
		}

		//! Combines Matrix-vector & dot products for greater efficiency (PARALLEL) (2D)
		static SOLVER matrix_vector__dot_product_PAR_2D(field_dense_2D<SOLVER>& result, field_dense_2D<SOLVER>& x, dMatrix_Laplacian5<SYSTEM>& A)
		{
			size_t height = x.get_dataWindow().y();
			size_t width = x.get_dataWindow().x();

			return parallel_reduce(blocked_range2d<size_t, size_t>(0, height, 0, width), SOLVER(0),
				[&](const blocked_range2d<size_t, size_t> &r, SOLVER init)->SOLVER {
				//make local copies to tell compiler to hold these in registers
				SOLVER local_init = init;
				size_t l_width = width;
				for (size_t j = r.rows().begin(); j != r.rows().end(); ++j)
				{
					for (size_t i = r.cols().begin(); i != r.cols().end(); ++i)
					{
						size_t index = I_2_2D(i, j, l_width);
						result.set(index, matrix_vector_product_row_inner_2D(A, x, i, j, index));
						local_init += result.get(index) * x.get(index);
					}
				}
				return local_init;
			}, [](SOLVER x, SOLVER y)->SOLVER {
				return x + y;
			});
		}

		//! Combines addition, subtraction, scaling and infinity norm calculation (SEQUENTIAL)
		static dVector<SOLVER, 2> add_subtr_scale_infnorm_SEQ(field_dense<SOLVER>& residual, field_dense<SOLVER>& x, field_dense<SOLVER>& d, field_dense<SOLVER>& ad, SOLVER alpha)
		{
			dVector<SOLVER, 2> sum_inf;
			sum_inf[0] = (SOLVER)0;
			sum_inf[1] = (SOLVER)0;

			for (size_t i = 0; i < d.get_size(); ++i)
			{
				//add & scale
				x.set(i, x.get(i) + d.get(i) * alpha);

				//subtract & scale
				residual.set(i, residual.get(i) - ad.get(i) * alpha);

				//dot product
				sum_inf[0] += residual.get(i) * residual.get(i);

				//infinity norm
				sum_inf[1] = (residual.get(i) > sum_inf[1]) ? residual.get(i) : sum_inf[1];

			}
			return sum_inf;
		}

		//! Combines addition, subtraction, scaling and infinity norm calculation (PAR)
		static dVector<SOLVER, 2> add_subtr_scale_infnorm_PAR(field_dense<SOLVER>& residual, field_dense<SOLVER>& x, field_dense<SOLVER>& d, field_dense<SOLVER>& ad, SOLVER alpha)
		{
			//zero vector to initialise calculation
			dVector<SOLVER, 2> initialise;
			initialise[0] = (SOLVER)0;
			initialise[1] = (SOLVER)0;

			return parallel_reduce(blocked_range<size_t>(0, d.get_size()), initialise,
				[&](const blocked_range<size_t> &r, dVector<SOLVER, 2> init)->dVector<SOLVER, 2> {
				dVector<SOLVER, 2> local_init = init;
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					//add & scale
					x.set(i, x.get(i) + d.get(i) * alpha);

					//subtract & scale
					residual.set(i, residual.get(i) - ad.get(i) * alpha);

					//dot product
					local_init[0] += residual.get(i) * residual.get(i);

					//infinity norm
					local_init[1] = (residual.get(i) > local_init[1]) ? residual.get(i) : local_init[1];
				}
				return local_init;
			}, [](dVector<SOLVER, 2> x, dVector<SOLVER, 2> y)->dVector<SOLVER, 2> {
				dVector<SOLVER, 2> z;
				z[0] = (SOLVER)0;
				z[1] = (SOLVER)0;

				z[0] += x[0] + y[0];
				z[1] = (x[1] > y[1]) ? x[1] : y[1];
				return z;
			});
		}

		//! Combines addition, subtraction, scaling and infinity norm calculation (PAR) (2D)
		static dVector<SOLVER, 2> add_subtr_scale_infnorm_PAR_2D(field_dense_2D<SOLVER>& residual, field_dense_2D<SOLVER>& x, field_dense_2D<SOLVER>& d, field_dense_2D<SOLVER>& ad, SOLVER alpha)
		{
			//zero vector to initialise calculation
			dVector<SOLVER, 2> initialise;
			initialise[0] = (SOLVER)0;
			initialise[1] = (SOLVER)0;

			return parallel_reduce(blocked_range<size_t>(0, d.get_size()), initialise,
				[&](const blocked_range<size_t> &r, dVector<SOLVER, 2> init)->dVector<SOLVER, 2> {
				dVector<SOLVER, 2> local_init = init;
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					//add & scale
					x.set(i, x.get(i) + d.get(i) * alpha);

					//subtract & scale
					residual.set(i, residual.get(i) - ad.get(i) * alpha);

					//dot product
					local_init[0] += residual.get(i) * residual.get(i);

					//infinity norm
					local_init[1] = (residual.get(i) > local_init[1]) ? residual.get(i) : local_init[1];
				}
				return local_init;
			}, [](dVector<SOLVER, 2> x, dVector<SOLVER, 2> y)->dVector<SOLVER, 2> {
				dVector<SOLVER, 2> z;
				z[0] = (SOLVER)0;
				z[1] = (SOLVER)0;

				z[0] += x[0] + y[0];
				z[1] = (x[1] > y[1]) ? x[1] : y[1];
				return z;
			});
		}

		//! Combines addition and subtraction of two dense fields (SEQUENTIAL)
		static void add_scale_SEQ(field_dense<SOLVER>& d, field_dense<SOLVER>& residual, SOLVER beta)
		{
			for (int i = 0; i < d.get_size(); ++i)
			{
				d.set(i, residual.get(i) + d.get(i) * beta);
			}
		}

		//! Combines addition and subtraction of two dense fields (PARALLEL)
		static void add_scale_PAR(field_dense<SOLVER>& d, field_dense<SOLVER>& residual, SOLVER beta)
		{
			parallel_for(blocked_range<size_t>(0, d.get_size()),
				[&](blocked_range<size_t> &r) {
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					d.set(i, residual.get(i) + d.get(i) * beta);
				}
			});
		}

		//! Combines addition and subtraction of two dense fields (PARALLEL)
		static void add_scale_PAR_2D(field_dense_2D<SOLVER>& d, field_dense_2D<SOLVER>& residual, SOLVER beta)
		{
			parallel_for(blocked_range<size_t>(0, d.get_size()),
				[&](blocked_range<size_t> &r) {
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					d.set(i, residual.get(i) + d.get(i) * beta);
				}
			});
		}


	private:
		static inline SOLVER matrix_vector_product_row_inner(dMatrix_Laplacian7<SYSTEM>& A, field_dense<SOLVER>& x, size_t i, size_t j, size_t k, size_t index)
		{
			SOLVER sum = 0.0f;

			//Always exists
			sum += (SOLVER)A.get_Diag(index) * x.get(index);

			//cout << index << ": ";

			//Rearrange in this way for possible cache-exploitation
			if (i < A.get_system_width_MV())
			{
				sum += (SOLVER)A.get_PlusI(index) * x.get(index + 1);
				//cout << A->get_PlusI(index) * x[index + 1] << "(+I) ";

			}
			if (i > 0)
			{
				sum += (SOLVER)A.get_MinusI(index) * x.get(index - 1);
				//cout << " + " << A->get_MinusI(index) * x[index - 1] << "(-I) ";
			}

			if (j < A.get_system_height_MV())
			{
				sum += (SOLVER)A.get_PlusJ(index) * x.get(index + A.get_system_width());
				//cout << " + " << A->get_PlusJ(index) * x[index + A->get_system_width()] << "(+J) ";
			}
			if (j > 0)
			{
				sum += (SOLVER)A.get_MinusJ(index) * x.get(index - A.get_system_width());
				//cout << " + " << A->get_MinusJ(index) * x[index - A->get_system_width()] << "(-J) ";
			}

			if (k < (SOLVER)A.get_system_depth_MV())
			{
				sum += A.get_PlusK(index) * x.get(index + A.get_system_slice());
				//cout << " + " << A->get_PlusK(index) * x[index + A->get_system_width_slice()] << "(+K) ";
			}
			if (k > 0)
			{
				sum += (SOLVER)A.get_MinusK(index) * x.get(index - A.get_system_slice());
				//cout << " + " << A->get_MinusK(index) * x[index - A->get_system_slice()] << "(-K) ";
			}

			//cout << " = " << sum << "\n";
			return sum;
		}

		static inline size_t I(const size_t i, const size_t j, const size_t k, const size_t width, const size_t height)
		{
			return width*height*k + width*j + i;
		}

		static inline size_t I_2(const size_t i, const size_t j, const size_t k, const size_t width, const size_t slice)
		{
			return slice*k + width*j + i;
		}

		static inline size_t I_tile(const size_t i, const size_t j, const size_t k, const size_t tile_dim)
		{
			return tile_dim*tile_dim*k + tile_dim*j + i;
		}

		static inline size_t I_tile_2(const size_t i, const size_t j, const size_t k, const size_t tile_dim, const size_t tile_dim_sqr)
		{
			return tile_dim_sqr*k + tile_dim*j + i;
		}

		//2D VARIANTS------------------------------------------
		static inline SOLVER matrix_vector_product_row_inner_2D(dMatrix_Laplacian5<SYSTEM>& A, field_dense_2D<SOLVER>& x, size_t i, size_t j, size_t index)
		{
			SOLVER sum = 0.0f;

			//Always exists
			sum += (SOLVER)A.get_Diag(index) * x.get(index);

			//cout << index << ": ";

			//Rearrange in this way for possible cache-exploitation
			if (i < A.get_system_width_MV())
			{
				sum += (SOLVER)A.get_PlusI(index) * x.get(index + 1);
				//cout << A->get_PlusI(index) * x[index + 1] << "(+I) ";

			}
			if (i > 0)
			{
				sum += (SOLVER)A.get_MinusI(index) * x.get(index - 1);
				//cout << " + " << A->get_MinusI(index) * x[index - 1] << "(-I) ";
			}

			if (j < A.get_system_height_MV())
			{
				sum += (SOLVER)A.get_PlusJ(index) * x.get(index + A.get_system_width());
				//cout << " + " << A->get_PlusJ(index) * x[index + A->get_system_width()] << "(+J) ";
			}
			if (j > 0)
			{
				sum += (SOLVER)A.get_MinusJ(index) * x.get(index - A.get_system_width());
				//cout << " + " << A->get_MinusJ(index) * x[index - A->get_system_width()] << "(-J) ";
			}

			//cout << " = " << sum << "\n";
			return sum;
		}

		static inline size_t I_2_2D(const size_t i, const size_t j, const size_t width)
		{
			return width * j + i;
		}


	}; /*class dBLAS_CG*/

}/*namespace drogon*/

#endif /*DROGON_DBLAS_CG_H*/