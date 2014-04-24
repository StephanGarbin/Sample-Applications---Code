#pragma once

#ifndef DROGON_FIELD_BASE_DENSE_H
#define DROGON_FIELD_BASE_DENSE_H

//STL

//TBB


//Drogon
#include "field_base.h"

namespace drogon {

	template<class T_DATA, class T_COORD>
	class field_base_dense : public field_base<T_DATA, T_COORD>
	{
	public:
		//-----------------------------------CONSTRUCTOR VERSIONS
		field_base_dense(size_t size, string name, FIELD_EXTENT data_window, FIELD_EXTENT display_window, FIELD_EXTENT simulation_window) : field_base(name, data_window, display_window, simulation_window), m_size(size)
		{
			//cout << "ctor field_base_dense:" << m_name << endl;
		}

		//! Destructor
		virtual ~field_base_dense()
		{
			//cout << "dtor field_base_dense: " << m_name << endl;
		};


		//-----------------------------------ACCESSORS SPECIFIC TO DENSE FIELDS

		//! Returns the total length of the underlying 1d vector.
		size_t get_size() const { return m_size; };

	protected:
		//-----------------------------------MISCELLANEOUS SPECIFIC TO DENSE FIELDS

		//! Total length of the underlying 1d vector.
		size_t m_size;


	}; /*class field_base_dense*/

} /*namespace drogon*/

#endif /*DROGON_FIELD_BASE_DENSE_H*/