#pragma once

#ifndef DROGON_EXTENTS_H
#define DROGON_EXTENTS_H

#include "dVector.h"

namespace drogon {

	//! This struct encapsulates two 3d vectors to store the extents of a cube.
	template<class T>
	struct extent3 {

		//! Min vector.
		dVector<T, 3> m_MinVector;
		//! Max vector.
		dVector<T, 3> m_MaxVector;

		inline T& minX()
		{
			return m_MinVector.x();
		}

		inline T& minY()
		{
			return m_MinVector.y();
		}

		inline T& minZ()
		{
			return m_MinVector.z();
		}

		inline T& maxX()
		{
			return m_MaxVector.x();
		}

		inline T& maxY()
		{
			return m_MaxVector.y();
		}

		inline T& maxZ()
		{
			return m_MaxVector.z();
		}

		inline T totalX()
		{
			return m_MaxVector.x() - m_MinVector.x();
		}

		inline T totalY()
		{
			return m_MaxVector.y() - m_MinVector.y();
		}

		inline T totalZ()
		{
			return m_MaxVector.z() - m_MinVector.z();
		}

	}; //!< class extent3


	//TYPEDEFS----------------------------------
	typedef extent3<size_t> FIELD_EXTENT;


} //!< namespace drogon


#endif//!< DROGON_EXTENTS_H