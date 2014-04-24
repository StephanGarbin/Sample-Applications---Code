#pragma once

#ifndef DROGON_FIELD_MAPPING_H
#define DROGON_FIELD_MAPPING_H

//STL

//TBB

//Drogon
#include <math.h>
#include "dVector.h"
#include "extents.h"

namespace drogon {

	//! This class implements coordinate transforms to, from and between world, local and tile space.
	//! It also handles all field / tile-internal index conversion.
	template<class WORLD_T, class LOCAL_T>
	class field_mapping
	{
	public:
		field_mapping()
		{
		}

		//! Constructor for dense fields
		field_mapping(dVector3_size_t data_window)
		{
			update_dataWindow(data_window);
		}


		//! Constructor for sparse fields
		field_mapping(dVector3_size_t data_window, dVector3_size_t data_window_tiles, size_t tile_edge)
		{
			update_dataWindow(data_window);
			update_dataWindowTiles(data_window_tiles);
			update_tile_edge(tile_edge);
		}

		//! cpy ctor
		field_mapping(const field_mapping& other)
		{
			m_data_window = other.m_data_window;
			m_data_window_slice = other.m_data_window_slice;
			m_data_window_tiles = other.m_data_window_tiles;
			m_data_window_tiles_slice = other.m_data_window_tiles_slice;
		}

		//! Destructor
		~field_mapping()
		{

		}


		//-----------------------------------MAPPING

		//-----------------------------------MAPPING FROM LOCAL TO INTERNAL INDEXING - DENSE FIELDS
		//! Maps an index in local space to internal indexing for dense fields.
		//! Coordinates are automatically rounded. If this is not the desired behaviour, use the interpolator instead.
		inline void I_local(const LOCAL_T i, const LOCAL_T j, const LOCAL_T k, size_t &index) const
		{
			//By doing this, rounding occurs automatically
			size_t field_i = i * m_data_window[0];
			size_t field_j = j * m_data_window[1];
			size_t field_k = k * m_data_window[2];

			I_dense(field_i, field_j, field_k, index);
		}

		//! return by value overload
		inline size_t I_local(const LOCAL_T i, const LOCAL_T j, const LOCAL_T k)
		{
			//By doing this, rounding occurs automatically
			size_t field_i = i * m_data_window[0];
			size_t field_j = j * m_data_window[1];
			size_t field_k = k * m_data_window[2];

			return I_dense(field_i, field_j, field_k);
		}


		//-----------------------------------MAPPING FROM LOCAL TO INTERNAL INDEXING - SPARSE FIELDS

		//! Maps an index in local space to internal indexing for sparse fields.
		//! Coordinates are automatically rounded. If this is not the desired behaviour, use the interpolator instead.
		inline void I_local_sparse(const LOCAL_T i, const LOCAL_T j, const LOCAL_T k, size_t &tile, size_t &index)
		{
			//By doing this, rounding occurs automatically
			size_t field_i = i * m_data_window[0];
			size_t field_j = j * m_data_window[1];
			size_t field_k = k * m_data_window[2];

			I_sparse(field_i, field_j, field_k, tile, index);
		}

		//-----------------------------------INTERNAL INDEXING - DENSE FIELDS

		//! This function converts a 3d to a 1d FIELD index
		inline void I_dense(const size_t i, const size_t j, const size_t k, size_t &index) const
		{
			//index = m_data_window.x * m_data_window.y * k + m_data_window.x * j + i;
			//More efficient:
			index = m_data_window_slice * k + m_data_window[0] * j + i;
		}

		//! return by value overload
		//3D---------------------------
		inline size_t I_dense(const size_t i, const size_t j, const size_t k) const
		{
			return m_data_window_slice * k + m_data_window[0] * j + i;
		}

		inline size_t I_dense(const size_t i, const size_t j, const size_t k)
		{
			return m_data_window_slice * k + m_data_window[0] * j + i;
		}
		//2D---------------------------
		inline size_t I_dense(const size_t i, const size_t j) const
		{
			return m_data_window[0] * j + i;
		}

		inline size_t I_dense(const size_t i, const size_t j)
		{
			return m_data_window[0] * j + i;
		}

		//-----------------------------------INTERNAL INDEXING - SPARSE FIELDS

		//! This function converts a 3d FIELD SPACE to FIELD-TILE SPACE and TILE SPACE index.
		inline void I_sparse(float i, float j, float k, size_t &tile, size_t &index)
		{
			//Tile Space
			float tile_i;
			float tile_j;
			float tile_k;

			//Field Tile Space & Tile Space assignment
			size_t field_tile_i = modf(i, &tile_i);
			size_t field_tile_j = modf(j, &tile_j);
			size_t field_tile_k = modf(k, &tile_k);

			//Calculate 1d Field Tile Space index
			I_sparse_fieldTileSpace(field_tile_i, field_tile_j, field_tile_k, tile);

			//Calculate 1d Tile Space index
			I_sparse_tileSpace(tile_i, tile_j, tile_k, index);
		}


		//! This function converts a 3d FIELD TILE SPACE  to a 1d FIELD TILE SPACE index.
		inline void I_sparse_fieldTileSpace(size_t field_tile_i, const size_t field_tile_j, const size_t field_tile_k, size_t &tile) const
		{
			//tile = m_data_window_tiles.x * m_data_window_tiles.y * field_tile_k + m_data_window_tiles.x * field_tile_j + field_tile_i;

			//Following should be quicker as avoiding function calls
			//tile = m_data_window_tiles.DATA[0] * m_data_window_tiles.DATA[1] * field_tile_k + m_data_window_tiles.DATA[0] * field_tile_j + field_tile_i;

			//Following should be quicker as also using the precomputed field slice
			tile = m_data_window_tiles_slice * field_tile_k + m_data_window_tiles[0] * field_tile_j + field_tile_i;
		}

		//! Return by value overload
		inline size_t I_sparse_fieldTileSpace(size_t field_tile_i, const size_t field_tile_j, const size_t field_tile_k)
		{
			//Following should be quicker as also using the precomputed field slice
			return m_data_window_tiles_slice * field_tile_k + m_data_window_tiles[0] * field_tile_j + field_tile_i;
		}

		//! This function converts a 3d TILE SPACE  to a 1d TILE SPACE index.
		inline void I_sparse_tileSpace(const size_t tile_i, const size_t tile_j, const size_t tile_k, size_t &index) const
		{
			//Use the precalculated m_tile_edge_square to avoid one multiplication. This is possible as each tile is a square cube.
			index = m_tile_edge_squared * tile_k + m_tile_edge * tile_j + tile_i;
		}

		//! Return by value overload
		inline size_t I_sparse_tileSpace(const size_t tile_i, const size_t tile_j, const size_t tile_k)
		{
			//Use the precalculated m_tile_edge_square to avoid one multiplication. This is possible as each tile is a square cube.
			return m_tile_edge_squared * tile_k + m_tile_edge * tile_j + tile_i;
		}

		//-----------------------------------ACCESSORS & MUTATORS

		//! This needs to be called if the origin of the field changes.
		void update_origin(const dVector<WORLD_T, 3> origin) const
		{
			m_origin = origin;
		}

		//! This needs to be called if the data window of the field changes.
		void update_dataWindow(const extent3<size_t> data_window) const
		{
			m_data_window = data_window;
			//Update slice for indexing calculations
			m_data_window_slice = m_data_window.x() * m_data_window.y();
		}

		//! This needs to be called if the data tile window of the field changes.
		void update_dataWindowTiles(const extent3<size_t> data_window_tiles) const
		{
			m_data_window_tiles = data_window_tiles;
			//Update slice for indexing calculations
			m_data_window_tiles_slice = m_data_window_tiles.x() * m_data_window_tiles.y();
		}

		//! This needs to be called if the tile dimension is changed.
		void update_tile_edge(const size_t tile_edge) const
		{
			m_tile_edge = tile_edge;
			//Update slice for indexing calculations
			m_tile_edge_squared = m_tile_edge * m_tile_edge;
		}

	private:
		//-----------------------------------TILE INFORMATION

		//! Tile edge length
		size_t m_tile_edge;

		//! Tile edge length square. This represents a 'slice' of tiles
		size_t m_tile_edge_squared;

		//-----------------------------------FIELD EXTENTS

		//! Data Window: The extent of data stored for this field in FIELD SPACE.
		//! For the purposes of the mapping class, this is ALWAYS the entire extent of the field starting from 0,0,0. This is why it is stored as dVector.
		dVector3_size_t m_data_window;

		//! Calculation of a 3-dimensional index always involves multiplying the depth coordinate (z) times and entire slice of the system, which can be precalculated.
		//! Represents on slice of the data window (x*y). The update function above MUST therefore be called if the data window changes.
		size_t m_data_window_slice;

		//! Data Window: The extent of data stored for this field in FIELD TILE SPACE.
		//! For the purposes of the mapping class, this is ALWAYS the entire extent of the field starting from 0,0,0. This is why it is stored as dVector.
		dVector3_size_t m_data_window_tiles;

		//! Represents on slice of the data tile window (x*y). The update function above MUST therefore be called if the data tile window changes.
		size_t m_data_window_tiles_slice;

		//! The field's world space origin
		dVector<WORLD_T, 3> m_origin;

	}; /*class field_mapping*/



} /*namespace drogon*/

#endif /*DROGON_FIELD_MAPPING_H*/