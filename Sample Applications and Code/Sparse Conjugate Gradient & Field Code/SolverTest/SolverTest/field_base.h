#pragma once

#ifndef DROGON_FIELD_BASE_H
#define DROGON_FIELD_BASE_H

//STL
#include <string>
#include <vector>

//TBB


//Drogon
#include "math.h"
#include "dVector.h"
#include "extents.h"
#include "field_mapping.h"
//#include "transform.h"

#include <iostream>

namespace drogon {

	//-----------------------------------TYPEDEFS

	template<class T_DATA, class T_COORD>
	class field_base
	{
	public:
		//-----------------------------------CONSTRUCTORS

		//!
		field_base(string name, FIELD_EXTENT data_window, FIELD_EXTENT display_window, FIELD_EXTENT simulation_window)
		{
			m_name = name;
			m_data_window = data_window;
			m_display_window = display_window;
			m_simulation_window = simulation_window;

			//cout << "ctor field_base: " << m_name << endl;
		}

		//! Destructor
		virtual ~field_base()
		{
			//cout << "dtor field_base: " << m_name << endl;
		}

		//-----------------------------------SPECIFIC TO ALL FIELDS

		//! Returns value at the given index.
		//! The value is returned by reference. This has consistently proven to be faster than using a return statement in tests.
		virtual inline void get(const size_t index, T_DATA &value) = 0;
		//virtual inline const T_DATA& get(const size_t index) const = 0;  

		//! Sets value at the given index.
		virtual inline void set(const size_t index, const T_DATA value) = 0;

		//-----------------------------------SPECIFIC TO 3D FIELDS

		virtual inline void get(const size_t i, const size_t j, const size_t k, T_DATA &value) = 0;
		//virtual inline const T_DATA& get(const size_t i, const size_t j, const size_t k) const = 0;

		virtual inline void set(const size_t i, const size_t j, const size_t k, const T_DATA value) = 0;

		//-----------------------------------FIELD RESIZING

		//! If the field is not surrounded on all sides by a number of empty cells equal to the padding value,
		//! it will be resized in that direction.
		//virtual inline void resizeField_padding() = 0;

		//-----------------------------------FIELD MEMBER ACCESSORS

		//! Returns the field's name
		string get_name() const { return m_name; }

		//! Returns the extent of data stored in the field.
		size_t get_fieldSpace() const { return m_data_window; }

		//! Returns the extent of tiles stored in the field.
		size_t get_fieldTileSpace() const { return m_data_window_tiles; }

		//! Returns the field's default value.
		T_DATA get_default_value() const { return m_defaul_value; }

		//! Returns a handle to the field's mapping class.
		const field_mapping<T_COORD, T_COORD> get_mapping() const { return m_mapping; }

		//! Returns the data window
		//! WORK on this... highly inefficient
		dVector3_size_t get_dataWindow()
		{
			dVector3_size_t temp;
			temp[0] = m_data_window.totalX();
			temp[1] = m_data_window.totalY();
			temp[2] = m_data_window.totalZ();

			return temp;
		}

		FIELD_EXTENT get_dataWindow_EXT() const
		{
			return m_data_window;
		}

		FIELD_EXTENT get_displayWindow_EXT() const
		{
			return m_display_window;
		}

		FIELD_EXTENT get_simulationWindow_EXT() const
		{
			return m_simulation_window;
		}

		//! Returns the field's padding value for any resize operations
		bool get_resizePadding() const { return m_resizePadding; };

		//! See whether or not the field is resizeable
		bool is_frozen() const { return m_frozen; }

		//-----------------------------------FIELD RESIZING MEMBER ACCESSORS

		//void freezeResize() { m_frozen = true; }

		//void unfreezeResize() { m_frozen = false; }

		//void set_resizePadding(const size_t padding) { m_resizePadding = padding; }

		//-----------------------------------FIELD DATA FUNCTIONS
		//! Sets the entire field to the input value
		virtual void set_field(T_DATA const value) = 0;

	protected:
		//-----------------------------------DATA STORAGE

		//Handled by the derived subclasses as this differs between dense and sparse fields

		//-----------------------------------MISCELLANEOUS

		//! The field's name
		string m_name;

		//! The field's instance of the mapping class
		field_mapping<T_COORD, T_COORD> m_mapping;

		//! The field's default value. This is returned in case no data is stored for an element and is used in default
		//! initialisation.
		T_DATA m_default_value;

		//! Determines whether the field will be resized if data is written to a location not yet allocated.
		bool m_frozen;

		//! Determines by how much the field will be resized if data is written to a location not yet allocated.
		size_t m_resize_padding;

		//-----------------------------------FIELD EXTENTS

		//! Data Window: The extent of data stored for this field in FIELD SPACE.
		FIELD_EXTENT m_data_window;

		//! Display Window: The extent of data displayed for this field.
		FIELD_EXTENT m_display_window;

		//! Simulation Window: The extent of data on which simulation is run for this field.
		FIELD_EXTENT m_simulation_window;

		//-----------------------------------WORLD TRANSFORM
		//transform<T_COORD> m_world_transform;



		//-----------------------------------FRIENDS

	private:
		//prevent slicing
		field_base(const field_base& other);
		field_base& operator=(const field_base& other);


	}; /*class field_base*/

} /*namespace drogon*/

#endif /*DROGON_FIELD_BASE_H*/