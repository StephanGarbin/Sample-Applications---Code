#pragma once

#ifndef DROGON_FIELD_TILE_H
#define DROGON_FIELD_TILE_H

#include "math.h"

namespace drogon {

	template<class T>
	class field_tile
	{
	public:
		//! ctor
		field_tile() : m_size(0), m_active(false)
		{
			//do not initalise any memory for values at this point
		}

		//! ctor intial value
		field_tile(size_t size) : m_size(size), m_active(false)
		{
			//do not initalise any memory for values at this point
		}

		//! cpy ctor
		field_tile(const field_tile& other) : m_size(other.m_size), m_active(other.m_active)
		{
			// NB: everything except m_DATA has been initialised
			if (m_active)
			{
				activate();
				for (size_t i = 0; i < m_size; ++i)
				{
					m_DATA[i] = other.m_DATA[i];
				}
			}
		}

		//! move ctor
		field_tile(field_tile&& other) : m_size(other.m_size), m_active(other.m_active)
		{
			// NB: everything except m_DATA has been initialised
			if (m_active)
			{
				m_DATA = other.m_DATA;

				other.m_active = false;
				other.m_DATA = nullptr;
				other.m_size = 0;
			}
		}

		//! dtor
		~field_tile()
		{
			if (m_active)
			{
				m_active = false;
				delete[] m_DATA;
			}
			/*m_size = 0;*/
		};

		//! Activate the tiles.
		void activate(size_t size)
		{
			m_size = size;

			m_DATA = new T[m_size];
			m_active = true;
		}

		//! Deactivate the tiles.
		inline void deactivate()
		{
			m_active = false;
			delete[] m_DATA;
		}

		//! Query whether current tile is activated.
		inline bool is_active() const
		{
			return m_active;
		}

		//! deactivates tile if all values are below tolerance; returns true if tile has been deactivated
		bool clear(const T tolerance)
		{
			bool below_tolerance = true;

			for (size_t i = 0; i < m_size; ++i)
			{
				if (abs(m_DATA[i]) > tolerance)
				{
					below_tolerance = false;
					break;
				}
			}

			if (below_tolerance)
			{
				deactivate();
				return true;
			}
			else
			{
				return false;
			}
		}

		//! Sets all entries in the tile to the input value.
		void set_tile(const T value)
		{
			for (size_t i = 0; i < m_size; ++i)
			{
				m_DATA[i] = value;
			}
		}

		//OPERATORS----------------------------------------------------

		//! = op
		field_tile& operator=(const field_tile& rhs)
		{
			//we assumes the tiles are of equal dimension
			if (rhs.m_active)
			{
				for (size_t i = 0; i < m_size; ++i)
				{
					m_DATA[i] = rhs.m_DATA[i];
				}
			}
			else
			{
				if (m_active)
				{
					deactivate();
				}
			}

			return *this;
		}

		//! = move op
		field_tile& operator=(field_tile&& rhs)
		{
			if (rhs.is_active())
			{
				m_DATA = rhs.m_DATA;

				rhs.m_active = false;
				rhs.m_DATA = nullptr;
				rhs.m_size = 0;
			}
			else
			{
				if (m_active)
				{
					deactivate();
				}
			}

			return *this;
		}

		//! Member access 1d
		inline T& operator[] (const size_t index)
		{
			return m_DATA[index];
		}

		//! Member access 1d
		inline const T& operator[] (const size_t index) const
		{
			return m_DATA[index];
		}


	private:
		bool m_active;
		size_t m_size;
		T *m_DATA;

	}; /*class field_tile*/

} /*namespace drogon*/


#endif /*DROGON_FIELD_TILE_H*/