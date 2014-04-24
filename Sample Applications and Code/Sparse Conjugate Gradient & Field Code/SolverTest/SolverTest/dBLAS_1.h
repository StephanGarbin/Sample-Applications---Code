#pragma once

#ifndef DROGON_DBLAS_1_H
#define DROGON_DBLAS_1_H

//STL

//TBB
#include "tbb\parallel_for.h"
#include "tbb\parallel_reduce.h"
#include "tbb\blocked_range.h"
#include "tbb\blocked_range2d.h"
#include "tbb\blocked_range3d.h"
#include "tbb\cache_aligned_allocator.h"

//DROGON
#include "field_dense.h"
#include "field_dense_2D.h"
#include "field_sparse.h"
#include "dMatrix_Laplacian7.h"

using namespace tbb;
using namespace std;

namespace drogon {

	template <class T>
	class dBLAS_1
	{
	public:
		dBLAS_1() {}
		~dBLAS_1() {}

		//! Dot product of two dense fields (SEQUENTIAL)
		static T dot_product_SEQ(field_dense<T>& a, field_dense<T>& b)
		{
			T sum = (T)0;
			for (size_t i = 0; i < a.get_size(); ++i)
			{
				sum += a.get(i) * b.get(i);
			}
			return sum;
		}

		//! Dot product of two sparse fields (SEQUENTIAL)
		static T dot_product_SEQ_sparse(field_sparse<T>& a, field_sparse<T>& b)
		{
			T sum = (T)0;
			for (int i = 0; i < a.get_num_tiles(); ++i)
			{
				//only perform calculation (& function call) if the tile is activated
				if (a.is_active(i) && b.is_active(i))
				{
					sum += dot_product_sparse_inner(a, b, i);
				}
			}

			return sum;
		}

		//! Dot product of two dense fields (PARALLEL)
		static T dot_product_PAR(field_dense<T>& a, field_dense<T>& b)
		{
			size_t n = a.get_size();

			//Make sure initial input is '(T)0' NOT '0.0f'
			//Otherwise floating point round-off error is incurred
			return parallel_reduce(blocked_range<size_t>(0, n), (T)0,
				[&](const blocked_range<size_t> &r, T init)->T {
				//make local copy to improve performance
				T local_init = init;
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					local_init += a.get(i) * b.get(i);
				}
				return local_init;
			}, [](T x, T y)->T {
				return x + y;
			});
		}

		//! Dot product of two dense fields (PARALLEL) (2D)
		static T dot_product_PAR_2D(field_dense_2D<T>& a, field_dense_2D<T>& b)
		{
			size_t n = a.get_size();

			//Make sure initial input is '(T)0' NOT '0.0f'
			//Otherwise floating point round-off error is incurred
			return parallel_reduce(blocked_range<size_t>(0, n), (T)0,
				[&](const blocked_range<size_t> &r, T init)->T {
				//make local copy to improve performance
				T local_init = init;
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					local_init += a.get(i) * b.get(i);
				}
				return local_init;
			}, [](T x, T y)->T {
				return x + y;
			});
		}

		//! Dot product of two sparse fields (PARALLEL)
		static T dot_product_PAR_sparse(field_sparse<T>& a, field_sparse<T>& b)
		{
			size_t n = a.get_num_tiles();

			//Make sure initial input is '(T)0' NOT '0.0f'
			//Otherwise floating point round-off error is incurred
			return parallel_reduce(blocked_range<size_t>(0, n), T(0),
				[&](const blocked_range<size_t> &r, T init)->T {
				T local_init = init;
				for (int i = r.begin(); i != r.end(); ++i)
				{
					//only perform calculation (& function call) if the tile is activated
					if (a.is_active(i) && b.is_active(i))
					{
						local_init += dot_product_sparse_inner(a, b, i);
					}
				}
				return local_init;
			}, [](T x, T y)->T {
				return x + y;
			});
		}

		//! Infinity norm of dense field (SEQUENTIAL)
		static T inf_norm_SEQ(field_dense<T>& a)
		{
			T infnorm = (T)0;
			for (size_t i = 0; i < a.get_size(); ++i)
			{
				infnorm = (a.get(i) > infnorm) ? a.get(i) : infnorm;
			}
			return infnorm;
		}

		//! Infinity norm of dense field (PARALLEL)
		static T inf_norm_PAR(field_dense<T>& a)
		{
			return parallel_reduce(blocked_range<size_t>(0, a.get_size()), T(0),
				[&](blocked_range<size_t> &r, T init)->T {
				T local_init = init;
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					local_init = (a.get(i) > init) ? a.get(i) : init;
				}
				return local_init;
			}, [](T x, T y) {
				return (x > y) ? x : y;
			});
		}

		//! Infinity norm of dense field (PARALLEL) (2D)
		static T inf_norm_PAR_2D(field_dense_2D<T>& a)
		{
			return parallel_reduce(blocked_range<size_t>(0, a.get_size()), T(0),
				[&](blocked_range<size_t> &r, T init)->T {
				T local_init = init;
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					local_init = (a.get(i) > init) ? a.get(i) : init;
				}
				return local_init;
			}, [](T x, T y) {
				return (x > y) ? x : y;
			});
		}


		//! Matrix-vector product for dense fields (SEQUENTIAL)
		static void matrix_vector_product_SEQ(dMatrix_Laplacian7<T>& A, field_dense<T>& x, field_dense<T>& result)
		{
			size_t slice = A.get_system_width() * A.get_system_height();
			size_t width = A.get_system_width();

			//put iterations this way to avoid cache misses
			for (size_t k = 0; k < A.get_system_depth(); ++k)
			{
				for (size_t j = 0; j < A.get_system_height(); ++j)
				{
					for (size_t i = 0; i < A.get_system_width(); ++i)
					{
						size_t index = I(i, j, k, width, slice);
						result.set(index, matrix_vector_product_row_inner(A, x, i, j, k, index));
					}
				}
			}
		}

		//! Matrix-vector product for sparse fields (SEQUENTIAL)
		static void matrix_vector_product_SEQ_sparse(dMatrix_Laplacian7<T>& A, field_sparse<T>& x, field_dense<T>& result)
		{
			size_t tiles_width = x.get_tiles_window().x();
			size_t tiles_height = x.get_tiles_window().y();
			size_t tiles_depth = x.get_tiles_window().z();
			size_t slice = A.get_system_width() * A.get_system_height();

			//loop over tiles
			for (size_t t_k = 0; t_k < tiles_depth; ++t_k)
			{
				for (size_t t_j = 0; t_j < tiles_height; ++t_j)
				{
					for (size_t t_i = 0; t_i < tiles_width; ++t_i)
					{
						size_t tile_index = I(t_i, t_j, t_k, tiles_width, tiles_height * tiles_width);
						size_t tile_dim = x.get_tile_edgeLength();
						if (x.is_active(tile_index))
						{
							//only the borders are problematic in terms of inter-tile overlap!
							for (size_t k = 0; k < tile_dim; ++k)
							{
								for (size_t j = 0; j < tile_dim; ++j)
								{
									for (size_t i = 0; i < tile_dim; ++i)
									{
										size_t global_i = t_i * tile_dim + i;
										size_t global_j = t_j * tile_dim + j;
										size_t global_k = t_k * tile_dim + k;
										size_t global_index = I(global_i, global_j, global_k, A.get_system_width(), slice);
										size_t local_index = I(i, j, k, tile_dim, tile_dim);
										//ifs are generally very quick to evaluate
										if (i == 0 | j == 0 | k == 0 | i == tile_dim | j == tile_dim | k == tile_dim)
										{
											result.set(global_index, matrix_vector_product_tile_inner_sparse(A, x, global_i, global_j, global_k, global_index, tile_index, local_index));
										}
										else
										{
											result.set(global_index, matrix_vector_product_tile_inner_sparse_amort(A, x, global_i, global_j, global_k, i, j, k, global_index, tile_index, local_index));
										}
									}
								}
							}
						}
						else
						{
							//set value in the result to zero (in case the vector has not been reset to 0)
							x.set_tile(tile_index, 0.0f);
						}
					}
				}
			}
		}

		static void matrix_vector_product_PAR(dMatrix_Laplacian7<T>& A, field_dense<T>& x, field_dense<T>& result)
		{
			size_t n = A.get_matrix_dim();
			size_t width = A.get_system_width();
			size_t height = A.get_system_height();

			parallel_for(blocked_range2d<size_t, size_t>(0, A.get_system_depth(), 0, A.get_system_height()),
				[&](blocked_range2d<size_t, size_t> &r) {
				for (size_t k = r.cols().begin(); k != r.cols().end(); ++k)
				{
					for (size_t j = r.rows().begin(); j != r.rows().end(); ++j)
					{
						//Do not parallelise innermost loop to increase cache efficiency
						for (size_t i = 0; i < width; ++i)
						{
							size_t index = I(i, j, k, width, height * width);
							result.set(index, matrix_vector_product_row_inner(A, x, i, j, k, index));
						}
					}
				}
			});
		}

	private:
		//! Calculates one pair of calculations in the dot-product
		static inline const T dot_product_sparse_inner(field_sparse<T>& a, field_sparse<T>& b, size_t tile)
		{
			T sum = (T)0;
			for (size_t i = 0; i < a.get_tile_size(); ++i)
			{
				sum += a.get_fast(tile, i) * b.get_fast(tile, i);
			}
			return sum;
		}

		//! Performs calculation of input vector with row index of the input matrix
		static inline T matrix_vector_product_row_inner(dMatrix_Laplacian7<T>& A, field_dense<T>& x, size_t i, size_t j, size_t k, size_t index)
		{
			T sum = (T)0;

			//Always exists
			sum += A.get_Diag(index) * x.get(index);

			//cout << index << ": ";

			//Rearrange in this way for possible cache-exploitation
			if (i < A.get_system_width_MV())
			{
				sum += A.get_PlusI(index) * x.get(index + 1);
				//cout << A->get_PlusI(index) * x[index + 1] << "(+I) ";

			}
			if (i > 0)
			{
				sum += A.get_MinusI(index) * x.get(index - 1);
				//cout << " + " << A->get_MinusI(index) * x[index - 1] << "(-I) ";
			}

			if (j < A.get_system_height_MV())
			{
				sum += A.get_PlusJ(index) * x.get(index + A.get_system_width());
				//cout << " + " << A->get_PlusJ(index) * x[index + A->get_system_width()] << "(+J) ";
			}
			if (j > 0)
			{
				sum += A.get_MinusJ(index) * x.get(index - A.get_system_width());
				//cout << " + " << A->get_MinusJ(index) * x[index - A->get_system_width()] << "(-J) ";
			}

			if (k < A.get_system_depth_MV())
			{
				sum += A.get_PlusK(index) * x.get(index + A.get_system_slice());
				//cout << " + " << A->get_PlusK(index) * x[index + A->get_system_width_slice()] << "(+K) ";
			}
			if (k > 0)
			{
				sum += A.get_MinusK(index) * x.get(index - A.get_system_slice());
				//cout << " + " << A->get_MinusK(index) * x[index - A->get_system_slice()] << "(-K) ";
			}

			//cout << " = " << sum << "\n";
			return sum;
		}

		//! Performs calculation of input vector with row index of the input matrix for sparse fields
		static inline T matrix_vector_product_tile_inner_sparse(dMatrix_Laplacian7<T>& A, field_sparse<T>& x, size_t global_i, size_t global_j, size_t global_k, size_t global_index, size_t tile_index, size_t local_index)
		{
			T sum = (T)0.0;

			//Always exists
			sum += A.get_Diag(global_index) * x.get_fast(tile_index, local_index);

			//cout << index << ": ";

			//Rearrange in this way for possible cache-exploitation
			if (global_i < A.get_system_width_MV())
			{
				sum += A.get_PlusI(global_index) * x.get(global_index + 1);
				//cout << A->get_PlusI(index) * x[index + 1] << "(+I) ";

			}
			if (global_i > 0)
			{
				sum += A.get_MinusI(global_index) * x.get(global_index - 1);
				//cout << " + " << A->get_MinusI(index) * x[index - 1] << "(-I) ";
			}

			if (global_j < A.get_system_height_MV())
			{
				sum += A.get_PlusJ(global_index) * x.get(global_index + A.get_system_width());
				//cout << " + " << A->get_PlusJ(index) * x[index + A->get_system_width()] << "(+J) ";
			}
			if (global_j > 0)
			{
				sum += A.get_MinusJ(global_index) * x.get(global_index - A.get_system_width());
				//cout << " + " << A->get_MinusJ(index) * x[index - A->get_system_width()] << "(-J) ";
			}

			if (global_k < A.get_system_depth_MV())
			{
				sum += A.get_PlusK(global_index) * x.get(global_index + A.get_system_slice());
				//cout << " + " << A->get_PlusK(index) * x[index + A->get_system_width_slice()] << "(+K) ";
			}
			if (global_k > 0)
			{
				sum += A.get_MinusK(global_index) * x.get(global_index - A.get_system_slice());
				//cout << " + " << A->get_MinusK(index) * x[index - A->get_system_slice()] << "(-K) ";
			}

			//cout << " = " << sum << "\n";
			return sum;
		}

		//! Performs calculation of input vector with row index of the input matrix for sparse fields. This is a faster version
		static inline T matrix_vector_product_tile_inner_sparse_amort(dMatrix_Laplacian7<T>& A, field_sparse<T>& x, size_t global_i, size_t global_j, size_t global_k, size_t local_i, size_t local_j, size_t local_k, size_t global_index, size_t tile_index, size_t local_index)
		{
			T sum = (T)0.0;
			size_t dim = x.get_tile_edgeLength();

			//Always exists
			sum += A.get_Diag(global_index) * x.get_fast(tile_index, local_index);

			//cout << index << ": ";

			//Rearrange in this way for possible cache-exploitation
			if (global_i < A.get_system_width_MV())
			{
				sum += A.get_PlusI(global_index) * x.get_fast(tile_index, I(local_i + 1, local_j, local_k, dim, dim));
				//cout << A->get_PlusI(index) * x[index + 1] << "(+I) ";

			}
			if (global_i > 0)
			{
				sum += A.get_MinusI(global_index) * x.get_fast(tile_index, I(local_i - 1, local_j, local_k, dim, dim));
				//cout << " + " << A->get_MinusI(index) * x[index - 1] << "(-I) ";
			}

			if (global_j < A.get_system_height_MV())
			{
				sum += A.get_PlusJ(global_index) * x.get_fast(tile_index, I(local_i, local_j + 1, local_k, dim, dim));
				//cout << " + " << A->get_PlusJ(index) * x[index + A->get_system_width()] << "(+J) ";
			}
			if (global_j > 0)
			{
				sum += A.get_MinusJ(global_index) * x.get_fast(tile_index, I(local_i, local_j - 1, local_k, dim, dim));
				//cout << " + " << A->get_MinusJ(index) * x[index - A->get_system_width()] << "(-J) ";
			}

			if (global_k < A.get_system_depth_MV())
			{
				sum += A.get_PlusK(global_index) * x.get_fast(tile_index, I(local_i, local_j, local_k + 1, dim, dim));
				//cout << " + " << A->get_PlusK(index) * x[index + A->get_system_width_slice()] << "(+K) ";
			}
			if (global_k > 0)
			{
				sum += A.get_MinusK(global_index) * x.get_fast(tile_index, I(local_i, local_j, local_k - 1, dim, dim));
				//cout << " + " << A->get_MinusK(index) * x[index - A->get_system_slice()] << "(-K) ";
			}

			//cout << " = " << sum << "\n";
			return sum;
		}


		//! Local index conversion
		static inline size_t I(const size_t i, const size_t j, const size_t k, const size_t width, const size_t slice)
		{
			return slice * k + width * j + i;
		}
	}; /*class dBLAS_1*/

} /*namespace drogon*/

#endif /*DROGON_DBLAS_1_H*/
