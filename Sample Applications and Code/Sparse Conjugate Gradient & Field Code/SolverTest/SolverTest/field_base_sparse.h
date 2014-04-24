#ifndef DROGON_FIELD_BASE_SPARSE_H
#define DROGON_FIELD_BASE_SPARSE_H

//STL

//TBB


//Drogon
#include "field_base.h"

namespace drogon {

	template<class T_DATA, class T_COORD>
	class field_base_sparse : public field_base<T_DATA, T_COORD>
	{
	public:
		//-----------------------------------CONSTRUCTOR VERSIONS

		//! Constructor for sparse 3d fields.
		field_base_sparse(string name, FIELD_EXTENT data_window, FIELD_EXTENT display_window, FIELD_EXTENT simulation_window, size_t tile_dim) : field_base(name, data_window, display_window, simulation_window)
		{
			m_tile_edgeLength = tile_dim;
			m_tile_size = tile_dim * tile_dim * tile_dim;
			m_num_tiles = (data_window.maxX() - data_window.minX()) * (data_window.maxY() - data_window.minY()) * (data_window.maxZ() - data_window.minZ()) / m_tile_size;
			m_data_window_tiles[0] = (data_window.maxX() - data_window.minX()) / tile_dim;
			m_data_window_tiles[1] = (data_window.maxY() - data_window.minY()) / tile_dim;
			m_data_window_tiles[2] = (data_window.maxZ() - data_window.minZ()) / tile_dim;
		}

		//! Destructor
		~field_base_sparse()
		{

		}


		//-----------------------------------SPECIFIC TO ALL SPARSE FIELDS

		//! Returns value at the given index & tile. Access is faster if FIELD TILE & TILE INDEX are known.
		virtual inline void get_fast(const size_t field_tile, const size_t tile_index, T_DATA &value) = 0;

		//! Sets value at the given index & tile. Access is faster if FIELD TILE & TILE INDEX are known.
		virtual inline void set_fast(const size_t field_tile, const size_t tile_index, const T_DATA value) = 0;


		//-----------------------------------SPECIFIC TO SPARSE 3D FIELDS

		//! Returns value at the given index & tile. Access is faster if FIELD TILE & TILE INDICES are known.
		virtual inline void get_fast(const size_t field_tile, const size_t tile_i, const size_t tile_j, const size_t tile_k, T_DATA &value) = 0;

		//! Sets value at the given index & tile. Access is faster if FIELD TILE & TILE INDICES are known.
		virtual inline void set_fast(const size_t field_tile, const size_t tile_i, const size_t tile_j, const size_t tile_k, const T_DATA value) = 0;


		//-----------------------------------ACCESSORS SPECIFIC TO SPARSE FIELDS

		//! Returns the total number of tiles for this field.
		inline size_t get_num_tiles() const { return m_num_tiles; }

		//! Returns the tile edge length (=dimension) for this field.
		inline size_t get_tile_edgeLength() const { return m_tile_edgeLength; }

		//! Returns the total length of the 1d vector underlying each tile.
		inline size_t get_tile_size() const { return m_tile_size; }

		//! Returns a const ref to the data_window containing tile information
		inline dVector3_size_t get_tiles_window()
		{
			return m_data_window_tiles;
		}

	protected:

		//-----------------------------------MISCELLANEOUS SPECIFIC TO SPARSE FIELDS

		//! Total number of tiles for this field.
		size_t m_num_tiles;

		//! Tile edge length (=dimension) for this field.
		size_t m_tile_edgeLength;

		//! Total length of the 1d vector underlying the each tile.
		size_t m_tile_size;

		//! Data Window: The extent of data stored for this field in FIELD TILE SPACE.
		//! This is of type dVector as it is always kept in the range from 0 - x
		dVector3_size_t m_data_window_tiles;


	}; /*class field_base_sparse*/

} /*namespace drogon*/

#endif /*DROGON_FIELD_BASE_SPARSE_H*/