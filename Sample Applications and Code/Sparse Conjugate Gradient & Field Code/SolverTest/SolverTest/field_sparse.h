#pragma once


#ifndef DROGON_FIELD_SPARSE_H
#define DROGON_FIELD_SPARSE_H

//STL

//TBB

//DROGON
#include "field_base_sparse.h"
#include "field_tile.h"

namespace drogon {

	template<class T_DATA>
	class field_sparse : public field_base_sparse<T_DATA, float>
	{
	public:
		field_sparse(string name, FIELD_EXTENT data_window, FIELD_EXTENT display_window, FIELD_EXTENT simulation_window, size_t tile_dim) : field_base_sparse(name, data_window, display_window, simulation_window, tile_dim)
		{
			m_DATA.resize(m_num_tiles);
		}

		~field_sparse()
		{
			m_DATA.clear();
		}

		//-----------------------------------SPECIFIC TO ALL SPARSE FIELDS

		//! Returns value at the given index & tile. Access is faster if FIELD TILE & TILE INDEX are known.
		inline void get_fast(const size_t field_tile, const size_t tile_index, T_DATA &value)
		{
			value = m_DATA[field_tile][tile_index];
		}

		//! Returns value at the given index & tile. Access is faster if FIELD TILE & TILE INDEX are known.
		inline const T_DATA& get_fast(const size_t field_tile, const size_t tile_index) const
		{
			return m_DATA[field_tile][tile_index];
		}

		//! Sets value at the given index & tile. Access is faster if FIELD TILE & TILE INDEX are known.
		inline void set_fast(const size_t field_tile, const size_t tile_index, const T_DATA value)
		{
			m_DATA[field_tile][tile_index] = value;
		}


		//-----------------------------------SPECIFIC TO SPARSE 3D FIELDS

		//! Returns value at the given index & tile. Access is faster if FIELD TILE & TILE INDICES are known.
		inline void get_fast(const size_t field_tile, const size_t tile_i, const size_t tile_j, const size_t tile_k, T_DATA &value)
		{
			value = m_DATA[field_tile][m_mapping.I_sparse_tileSpace(tile_i, tile_j, tile_k)];
		}

		//! Sets value at the given index & tile. Access is faster if FIELD TILE & TILE INDICES are known.
		inline void set_fast(const size_t field_tile, const size_t tile_i, const size_t tile_j, const size_t tile_k, const T_DATA value)
		{
			m_DATA[field_tile][m_mapping.I_sparse_tileSpace(tile_i, tile_j, tile_k)] = value;
		}

		//-----------------------------------SPECIFIC TO ALL FIELDS

		//! Returns value at the given index.
		inline void get(const size_t index, T_DATA &value)
		{
			//see which tile we're in
			size_t tile_index = index / m_tile_edgeLength;

			//check if tile contains any data
			if (is_active(tile_index))
			{
				//calculate local index
				size_t local_index = index%m_tile_edgeLength;
				value = m_DATA[tile_index][local_index];
			}
			else
			{
				value = m_default_value;
			}
		}

		//! Returns value at the given index.
		inline T_DATA get(const size_t index)
		{
			//see which tile we're in
			size_t tile_index = index / m_tile_size;

			//check if tile contains any data
			if (is_active(tile_index))
			{
				//calculate local index
				size_t local_index = index%m_tile_edgeLength;
				return m_DATA[tile_index][local_index];
			}
			else
			{
				return m_default_value;
			}
		}

		//! Sets value at the given index. 
		inline void set(const size_t index, const T_DATA value)
		{
			//Not sure what to do with this function yet
			size_t field_tile_space = index - ((index / m_tile_edgeLength) * m_tile_edgeLength);
			size_t tile_space = index - field_tile_space;
			m_DATA[field_tile_space][tile_space] = value;
		}

		//-----------------------------------SPECIFIC TO 3D FIELDS

		//! Returns value at the given index.
		inline void get(const size_t i, const size_t j, const size_t k, T_DATA &value)
		{
			size_t field_tile_space;
			size_t tile_space;
			m_mapping.I_sparse(i, j, k, field_tile_space, tile_space);

			value = m_DATA[field_tile_space][tile_space];
		}

		//! Sets value at the given index.
		inline void set(const size_t i, const size_t j, const size_t k, const T_DATA value)
		{
			size_t field_tile_space;
			size_t tile_space;
			m_mapping.I_sparse(i, j, k, field_tile_space, tile_space);

			m_DATA[field_tile_space][tile_space] = value;
		}

		//-----------------------------------FIELD DATA FUNCTIONS
		//! Sets the entire field to the input value
		void set_field(T_DATA const value)
		{
			for (size_t i = 0; i < m_num_tiles; ++i)
			{
				if (m_DATA[i].is_active())
				{
					m_DATA[i] = value;
				}
				else
				{
					m_DATA[i].activate(m_tile_size);
					m_DATA[i].set_tile(value);
				}
			}
		}

		//! Sets a specific tile of the field to the input value
		bool set_tile(const size_t tile_index, const T_DATA value)
		{
			if (is_active(tile_index))
			{
				for (size_t i = 0; i < m_tile_size; ++i)
				{
					m_DATA[tile_index][i] = value;
				}
				return true;
			}
			else
			{
				return false;
			}
		}

		//!determine whether a tile is activated
		inline bool is_active(const size_t index)
		{
			return m_DATA[index].is_active();
		}

	private:
		//-----------------------------------DATA STORAGE

		//! Array containing the tiles storing the field's data.
		vector<field_tile<T_DATA>> m_DATA;


	}; /*class field_sparse*/

} /*namespace drogon*/


#endif /*DROGON_FIELD_SPARSE_H*/