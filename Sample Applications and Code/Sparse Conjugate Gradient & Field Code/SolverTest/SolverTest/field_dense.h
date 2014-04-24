#pragma once //to stop intellisense warning...

#ifndef DROGON_FIELD_DENSE_H
#define DROGON_FIELD_DENSE_H

//STL
#include <iostream> //for debug only

//TBB

//DROGON
#include "field_base_dense.h"

namespace drogon {

	//! 3d Field dense field
	//! Internal indexing is size_t, Local and World access are float.
	template<class T_DATA>
	class field_dense : public field_base_dense<T_DATA, float>
	{
	public:
		//! ctor
		field_dense(string name, FIELD_EXTENT data_window, FIELD_EXTENT display_window, FIELD_EXTENT simulation_window) : field_base_dense((data_window.maxX() - data_window.minX()) * (data_window.maxY() - data_window.minY()) * (data_window.maxZ() - data_window.minZ()), name, data_window, display_window, simulation_window)
		{
			//use delegating constructors
			/*field_base_dense((data_window.max_x() - data_window.min_x()) * (data_window.max_y() - data_window.min_y()) * (data_window.max_z() - data_window.min_z()))*/
			//field_base(name, data_window, display_window, simulation_window);

			//allocate storage in vector
			m_DATA.resize(m_size);

			//cout << "ctor field_dense: " << m_name<< endl;
		}

		//! cpy ctor
		//field_dense(const field_dense& other)
		//{
		//	
		//}

		////! move ctor
		//field_dense(field_dense&& other)
		//{

		//}

		~field_dense()
		{
			m_DATA.clear();
			//cout << "dtor field_dense: " << m_name<< endl;
		}

		//-----------------------------------SPECIFIC TO ALL FIELDS

		//! Returns value at the given index.
		inline void get(const size_t index, T_DATA &value)
		{
			value = m_DATA[index];
		}

		//! Returns value at the given index.
		inline const T_DATA& get(const size_t index) const
		{
			return m_DATA[index];
		}

		//! Sets value at the given index. 
		inline void set(const size_t index, const T_DATA value)
		{
			m_DATA[index] = value;
		}

		//-----------------------------------SPECIFIC TO 3D FIELDS

		//! Returns value at the given index.
		inline void get(const size_t i, const size_t j, const size_t k, T_DATA &value)
		{
			value = m_DATA[m_mapping.I_dense(i, j, k)];
		}

		//! Returns value at the given index (for 3D).
		inline T_DATA get(const size_t i, const size_t j, const size_t k) const
		{
			return m_DATA[m_mapping.I_dense(i, j, k)];
		}

		//-----------------------------------SPECIFIC TO 2D FIELDS
		//! Returns value at the given index (for 2D).
		inline T_DATA get(const size_t i, const size_t j) const
		{
			return m_DATA[m_mapping.I_dense(i, j)];
		}

		//! Sets value at the given index (for 3D).
		inline void set(const size_t i, const size_t j, const size_t k, const T_DATA value)
		{
			m_DATA[m_mapping.I_dense(i, j, k)] = value;
		}

		//! Sets value at the given index (for 2D).
		inline void set(const size_t i, const size_t j, const T_DATA value)
		{
			m_DATA[m_mapping.I_dense(i, j)] = value;
		}



		//-----------------------------------FIELD DATA FUNCTIONS
		//! Sets the entire field to the input value
		void set_field(T_DATA const value)
		{
			for (size_t i = 0; i < m_size; ++i)
			{
				m_DATA[i] = value;
			}
		}

		void print_vector(const size_t start, size_t end)
		{
			end = (end > m_DATA.size()) ? m_DATA.size() : end;

			for (size_t i = start; i < end; ++i)
			{
				cout << "[" << i << "]: " << m_DATA[i] << endl;
			}
		}

		void copy_from(const field_dense& other, bool run_parallel)
		{
			m_DATA = other.m_DATA;
		}

		T_DATA max_value() const
		{
			T_DATA temp = 0;
			for (auto i : m_DATA)
			{
				temp = (i > temp) ? i : temp;
			}
			return temp;
		}

		size_t get_1D_size()
		{
			return m_DATA.size();
		}

	private:
		//-----------------------------------DATA STORAGE SPECIFIC TO DENSE FIELDS

		//! Array containing the field's data.
		vector<T_DATA> m_DATA;



	}; /*class field_dense*/

} /*namespace drogon*/

#endif /*DROGON_FIELD_DENSE_H*/