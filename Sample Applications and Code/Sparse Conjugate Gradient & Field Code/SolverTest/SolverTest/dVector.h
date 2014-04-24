#pragma once

#ifndef DROGON_DVECTOR_H
#define DROGON_DVECTOR_H

using namespace std;

namespace drogon {

	template<class T, size_t SIZE>
	class dVector
	{
	public:
		//! ctor
		dVector() : m_DATA(new T[SIZE]), m_size(SIZE)
		{
			//initialise vector to 0
			for (size_t i = 0; i < m_size; ++i)
			{
				m_DATA[i] = 0;
			}
		}

		//! cpy ctor
		dVector(const dVector<T, SIZE>& other) : m_DATA(new T[other.m_size]), m_size(other.m_size)
		{
			//copy values
			for (size_t i = 0; i < m_size; ++i)
			{
				m_DATA[i] = other.m_DATA[i];
			}
		}

		//! move ctor
		dVector(dVector<T, SIZE>&& src)
		{
			//m_DATA = src.m_DATA;
			//m_size = src.m_size;

			//src.m_DATA = 0;
			//src.m_DATA = nullptr;
			//src.m_size = 0;

			//This is redundant code; just call the move assignment operator
			*this = std::move(src);
		}

		//! ctor : initial value
		explicit dVector(T initial_value) : m_DATA(new T[SIZE]), m_size(SIZE)
		{
			cout << "value ctor called" << endl;
			//initialise all entries to the provided initial value
			for (size_t i = 0; i < m_size; ++i)
			{
				m_DATA[i] = initial_value;
			}
		}

		//!dtor
		~dVector()
		{
			if (m_DATA)
			{
				//delete the ptr
				delete[] m_DATA;
			}
		}

		//OPERATOR OVERLOADS-------------------------------------------------------

		//! assingment cpy op
		dVector& operator=(const dVector<T, SIZE>& rhs)
		{
			//if(rhs.m_size >= m_size)
			//{
			for (size_t i = 0; i < m_size; ++i)
			{
				m_DATA[i] = rhs.m_DATA[i];
			}
			//}

			m_size = rhs.m_size;

			return *this;
		}

		//! assingment move op
		dVector& operator=(dVector<T, SIZE>&& rhs)
		{
			//Otherwise we created a mem leak!
			delete[] m_DATA;

			m_DATA = rhs.m_DATA;
			m_size = rhs.m_size;

			rhs.m_DATA = nullptr;

			rhs.m_size = 0;

			return *this;
		}

		//! [] op
		T& operator[](const size_t index)
		{
			return m_DATA[index];
		}

		T& operator[](const size_t index) const
		{
			return m_DATA[index];
		}

		//! += op
		template<class X>
		dVector& operator+=(const dVector<X, SIZE>& rhs)
		{
			for (size_t i = 0; i < m_size; ++i)
			{
				m_DATA[i] += rhs[i];
			}
			return *this;
		}

		//! -= op
		template<class X>
		dVector& operator-=(const dVector<X, SIZE>& rhs)
		{
			for (size_t i = 0; i < m_size; ++i)
			{
				m_DATA[i] -= rhs[i];
			}
			return *this;
		}

		//! *= op
		template<class X>
		dVector& operator*=(const dVector<X, SIZE>& rhs)
		{
			for (size_t i = 0; i < m_size; ++i)
			{
				//Do not try to access data directly here as instance with different template is inaccessible
				m_DATA[i] *= rhs[i];
			}
			return *this;
		}

		//! /= op
		template<class X>
		dVector& operator/=(const dVector<X, SIZE>& rhs)
		{
			for (size_t i = 0; i < m_size; ++i)
			{
				m_DATA[i] /= rhs[i];
			}
			return *this;
		}

		//! += op for scalar values
		dVector& operator+=(const T rhs)
		{
			for (size_t i = 0; i < m_size; ++i)
			{
				m_DATA[i] += rhs
			}
			return *this;
		}

		//! -= op for scalar values
		dVector& operator-=(const T rhs)
		{
			for (size_t i = 0; i < m_size; ++i)
			{
				m_DATA[i] -= rhs;
			}
			return *this;
		}

		//! *= op for scalar values
		dVector& operator*=(const T rhs)
		{
			for (size_t i = 0; i < m_size; ++i)
			{
				m_DATA[i] *= rhs
			}
			return *this;
		}

		//! /= op for scalar values
		dVector& operator/=(const T rhs)
		{
			for (size_t i = 0; i < m_size; ++i)
			{
				m_DATA[i] /= rhs
			}
			return *this;
		}

		//! + op
		template<class X>
		dVector operator+(const dVector<X, SIZE>& rhs)
		{
			dVector<T, SIZE> temp;

			for (size_t i = 0; i < m_size; ++i)
			{
				temp.m_DATA[i] = m_DATA[i] + rhs[i];
			}

			return temp;
		}

		//! - op
		template<class X>
		dVector operator-(const dVector<X, SIZE>& rhs)
		{
			dVector<T, SIZE> temp;

			for (size_t i = 0; i < m_size; ++i)
			{
				temp.m_DATA[i] = m_DATA[i] - rhs[i];
			}

			return temp;
		}

		//! * op
		template<class X>
		dVector operator*(const dVector<X, SIZE>& rhs)
		{
			dVector<T, SIZE> temp;

			for (size_t i = 0; i < m_size; ++i)
			{
				temp.m_DATA[i] = m_DATA[i] * rhs[i];
			}

			return temp;
		}

		//! / op
		template<class X>
		dVector operator/(const dVector<X, SIZE>& rhs)
		{
			dVector<T, SIZE> temp;

			for (size_t i = 0; i < m_size; ++i)
			{
				temp.m_DATA[i] = m_DATA[i] / rhs[i];
			}

			return temp;
		}

		//! + op for scalar values
		dVector operator+(const T rhs)
		{
			dVector<T, SIZE> temp;

			for (size_t i = 0; i < m_size; ++i)
			{
				temp.m_DATA[i] = m_DATA[i] + rhs;
			}

			return temp;
		}

		//! - op for scalar values
		dVector operator-(const T rhs)
		{
			dVector<T, SIZE> temp;

			for (size_t i = 0; i < m_size; ++i)
			{
				temp.m_DATA[i] = m_DATA[i] - rhs;
			}

			return temp;
		}

		//! * op for scalar values
		dVector operator*(const T rhs)
		{
			dVector<T, SIZE> temp;

			for (size_t i = 0; i < m_size; ++i)
			{
				temp.m_DATA[i] = m_DATA[i] * rhs;
			}

			return temp;
		}

		//! / op for scalar values
		dVector operator/(const T rhs)
		{
			dVector<T, SIZE> temp;

			for (size_t i = 0; i < m_size; ++i)
			{
				temp.m_DATA[i] = m_DATA[i] / rhs;
			}

			return temp;
		}

		//! < op
		template<class X>
		bool operator<(const dVector<X, SIZE>& rhs)
		{
			bool temp = false;

			for (size_t i = 0; i < m_size; ++i)
			{
				temp = (m_DATA[i] < rhs[i]) ? true : false;
			}

			return temp;
		}

		//! > op
		template<class X>
		bool operator>(const dVector<X, SIZE>& rhs)
		{
			bool temp = false;

			for (size_t i = 0; i < m_size; ++i)
			{
				temp = (m_DATA[i] > rhs[i]) ? true : false;
			}

			return temp;
		}

		//! <= op
		template<class X>
		bool operator<=(const dVector<X, SIZE>& rhs)
		{
			bool temp = false;

			for (size_t i = 0; i < m_size; ++i)
			{
				temp = (m_DATA[i] <= rhs[i]) ? true : false;
			}

			return temp;
		}

		//! > op
		template<class X>
		bool operator>=(const dVector<X, SIZE>& rhs)
		{
			bool temp = false;

			for (size_t i = 0; i < m_size; ++i)
			{
				temp = (m_DATA[i] >= rhs[i]) ? true : false;
			}

			return temp;
		}

		//! < op for scalar values
		bool operator<(const T rhs)
		{
			bool temp = false;

			for (size_t i = 0; i < m_size; ++i)
			{
				temp = (m_DATA[i] < rhs) ? true : false;
			}

			return temp;
		}

		//! > op for scalar values
		bool operator>(const T rhs)
		{
			bool temp = false;

			for (size_t i = 0; i < m_size; ++i)
			{
				temp = (m_DATA[i] > rhs) ? true : false;
			}

			return temp;
		}

		//! <= op for scalar values
		bool operator<=(const T rhs)
		{
			bool temp = false;

			for (size_t i = 0; i < m_size; ++i)
			{
				temp = (m_DATA[i] <= rhs) ? true : false;
			}

			return temp;
		}

		//! > op for scalar values
		bool operator>=(const T rhs)
		{
			bool temp = false;

			for (size_t i = 0; i < m_size; ++i)
			{
				temp = (m_DATA[i] >= rhs) ? true : false;
			}

			return temp;
		}

		//! % op
		template<class X>
		dVector operator%(const dVector<X, SIZE>& rhs)
		{
			dVector<T, SIZE> temp;

			for (size_t i = 0; i < m_size; ++i)
			{
				temp.m_DATA[i] = m_DATA[i] % rhs[i];
			}

			return temp;
		}

		//! % op for scalar values
		dVector operator%(const T rhs)
		{
			dVector<T, SIZE> temp;

			for (size_t i = 0; i < m_size; ++i)
			{
				temp.m_DATA[i] = m_DATA[i] % rhs;
			}

			return temp;
		}

		//BLAS ROUTINES-------------------------------------------------------

		//! dot product
		template<class X>
		T dot_product(const dVector<X, SIZE>& other)
		{
			T temp = (T)0;

			for (size_t i = 0; i < m_size; ++i)
			{
				temp += m_DATA[i] * other[i];
			}
		}

		//! cross product
		template<class X>
		dVector cross_product(const dVector<X, SIZE>& other)
		{

		}

		//SPECIAL ACCESS OPERATORS-------------------------------------------------------

		T& x()
		{
			return m_DATA[0];
		}

		T& y()
		{
			return m_DATA[1];
		}

		T& z()
		{
			return m_DATA[2];
		}

	private:
		size_t m_size;
		T* m_DATA;

	}; //class dVector

	//-----------------------------------TYPEDEFS
	typedef dVector<size_t, 2> dVector2_size_t;
	typedef dVector<float, 2> dVector2_float;
	typedef dVector<double, 2> dVector2_double;

	typedef dVector<size_t, 3> dVector3_size_t;
	typedef dVector<float, 3> dVector3_float;
	typedef dVector<double, 3> dVector3_double;

} //namespace drogon

#endif //DROGON_DVECTOR_H