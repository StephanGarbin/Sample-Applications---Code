#pragma once

#ifndef DROGON_DMATRIX_LAPLACIAN7_H
#define DROGON_DMATRIX_LAPLACIAN7_H

//STL
#include <string>

//TBB
#include <tbb\parallel_for.h>
#include <tbb\cache_aligned_allocator.h>

//DROGON
#include "extents.h"
#include "dVector.h"

namespace drogon {

	template <class T>
	class dMatrix_Laplacian7
	{
	public:
		dMatrix_Laplacian7(string name, FIELD_EXTENT data_window, const size_t matrix_dim, const T default_value)
		{
			m_data_window = data_window;

			this->matrix_dim = matrix_dim;
			this->system_width = data_window.totalX();
			this->system_height = data_window.totalY();
			this->system_depth = data_window.totalZ();
			this->system_slice = system_width * system_height;
			this->default_value = default_value;

			this->system_width_MV = system_width - 1;
			this->system_height_MV = system_height - 1;
			this->system_depth_MV = system_depth - 1;

			//if(matrix_dim < 1e+5) //FOR TEST PURPOSES
			//{
			//	set_all(default_value);
			//}
			//else
			//{
			//	set_all_PAR(default_value);
			//}

			Diag.resize(matrix_dim, default_value);
			PlusI.resize(matrix_dim, default_value);
			PlusJ.resize(matrix_dim, default_value);
			PlusK.resize(matrix_dim, default_value);
			//PlusJ.resize(matrix_dim, default_value);
			//PlusK.resize(matrix_dim, default_value);
		}


		//! dtor
		~dMatrix_Laplacian7()
		{
			Diag.clear();
			PlusI.clear();
			PlusJ.clear();
			PlusK.clear();

			m_name.clear();
		}

		inline T get_Diag(const size_t index)
		{
			return Diag[index];
		}

		inline T get_PlusI(const size_t index)
		{
			return PlusI[index];
		}
		inline T get_PlusJ(const size_t index)
		{
			//NB vector is semi-sparse
			return PlusJ[index];
		}
		inline T get_PlusK(const size_t index)
		{
			//NB vector is semi-sparse
			return PlusK[index];
		}

		inline T get_MinusI(const size_t index)
		{
			return PlusI[index - 1];
		}
		inline T get_MinusJ(const size_t index)
		{
			//NB vector is semi-sparse
			return PlusJ[index - system_width];
		}
		inline T get_MinusK(const size_t index)
		{
			//NB vector is semi-sparse
			return PlusK[index - system_slice];
		}

		inline void set_Diag(const size_t index, const T value)
		{
			Diag[index] = value;
		}

		inline void set_PlusI(const size_t index, const T value)
		{
			PlusI[index] = value;
		}
		inline void set_PlusJ(const size_t index, const T value)
		{
			//NB vector is semi-sparse
			PlusJ[index] = value;
		}
		inline void set_PlusK(const size_t index, const T value)
		{
			//NB vector is semi-sparse
			PlusK[index] = value;
		}

		inline void set_MinusI(const size_t index, const T value)
		{
			PlusI[index - 1] = value;
		}
		inline void set_MinusJ(const size_t index, const T value)
		{
			//NB vector is semi-sparse
			PlusJ[index - system_width] = value;
		}
		inline void set_MinusK(const size_t index, const T value)
		{
			//NB vector is semi-sparse
			PlusK[index - system_slice] = value;
		}

		void set_Diag(const T value)
		{
			for (size_t i = 0; i < matrix_dim; ++i)
			{
				Diag[i] = value;
			}
		}
		void set_I(const T value)
		{
			for (size_t i = 0; i < matrix_dim; ++i)
			{
				PlusI[i] = value;
			}
		}
		void set_J(const T value)
		{
			for (size_t i = 0; i < matrix_dim - system_width; ++i)
			{
				PlusJ[i] = value;
			}
		}
		void set_K(const T value)
		{
			for (size_t i = 0; i < matrix_dim - system_slice; ++i)
			{
				PlusK[i] = value;
			}
		}

		inline size_t get_matrix_dim()
		{
			return matrix_dim;
		}
		inline size_t get_system_width()
		{
			return this->system_width;
		}
		inline size_t get_system_height()
		{
			return this->system_height;
		}
		inline size_t get_system_depth()
		{
			return this->system_depth;
		}
		inline size_t get_system_slice()
		{
			return system_slice;
		}
		inline size_t get_system_width_MV()
		{
			return system_width_MV;
		}
		inline size_t get_system_height_MV()
		{
			return system_height_MV;
		}
		inline size_t get_system_depth_MV()
		{
			return system_depth_MV;
		}

		//Diagnostics
		bool check_coefficients()
		{
			bool is_okay = true;

			for (size_t k = 0; k < system_depth; ++k)
			{
				for (size_t j = 0; j < system_height; ++j)
				{
					for (size_t i = 0; i < system_width; ++i)
					{
						T sum = 0;
						size_t index = I(i, j, k);

						if (i > 0)
						{
							sum += get_MinusI(index);
						}
						else if (i < system_width_MV)
						{
							sum += get_PlusI(index);
						}

						if (j > 0)
						{
							sum += get_MinusJ(index);
						}
						else if (j < system_height_MV)
						{
							sum += get_PlusJ(index);
						}

						if (k > 0)
						{
							sum += get_MinusK(index);
						}
						else if (k < system_depth_MV)
						{
							sum += get_PlusK(index);
						}

						//now check diagonal and correct if necessary

						if (get_Diag(index) < sum)
						{
							if (is_okay == true)
							{
								is_okay = false;
							}
							set_Diag(index, sum);
						}

						if (get_Diag(index) > sum)
						{
							set_Diag(index, sum);
						}
					}
				}
			}

			return is_okay;
		}

		std::vector<T, tbb::cache_aligned_allocator<T>> Diag;
		std::vector<T, tbb::cache_aligned_allocator<T>> PlusI;
		std::vector<T, tbb::cache_aligned_allocator<T>> PlusJ;
		std::vector<T, tbb::cache_aligned_allocator<T>> PlusK;

	private:
		void set_all(T value)
		{
			for (size_t i = 0; i < matrix_dim; ++i)
			{
				Diag[i] = value;
				PlusI[i] = value;
				if (i < matrix_dim - system_width)
				{
					PlusJ[i] = value;
				}
				if (i < matrix_dim - system_slice)
				{
					PlusK[i] = value;
				}
			}
		}
		void set_all_PAR(T value)
		{
			size_t n = matrix_dim;
			parallel_for(blocked_range<size_t>(0, n),
				[=](const blocked_range<size_t> &r) {
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					Diag[i] = value;
					PlusI[i] = value;
					if (i < matrix_dim - system_width)
					{
						PlusJ[i] = value;
					}
					if (i < matrix_dim - system_slice)
					{
						PlusK[i] = value;
					}
				}
			});
		}

		//Indexing functions
		inline size_t I(const size_t i, const size_t j, const size_t k)
		{
			return system_slice * k + system_width * j + i;
		}

		//! Data Window: The extent of data stored for this field in FIELD SPACE.
		FIELD_EXTENT m_data_window;

		//! The field's name
		string m_name;

		//width * length of the input system
		size_t system_slice;
		//absolute dimension of the matrix system
		size_t matrix_dim;

		//system dimensions
		size_t system_width;
		size_t system_height;
		size_t system_depth;


		//for BLAS_1
		size_t system_width_MV;
		size_t system_height_MV;
		size_t system_depth_MV;

		T default_value;
	};

} /*namespace drogon*/

#endif /*DROGON_DMATRIX_LAPLACIAN7_H*/