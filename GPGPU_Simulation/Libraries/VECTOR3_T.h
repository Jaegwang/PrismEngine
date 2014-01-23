

#pragma once 

#include "MATH_DEFINITION.h"

class VECTOR3_T;

static VECTOR3_T operator + ( const VECTOR3_T& v0, const VECTOR3_T& v1 ) restrict(amp,cpu);
static VECTOR3_T operator - ( const VECTOR3_T& v0, const VECTOR3_T& v1 ) restrict(amp,cpu);
static VECTOR3_T operator * ( const VECTOR3_T& v0, const VECTOR3_T& v1 ) restrict(amp,cpu);
static VECTOR3_T operator / ( const VECTOR3_T& v0, const VECTOR3_T& v1 ) restrict(amp,cpu);

// unary negation
static VECTOR3_T operator - ( const VECTOR3_T& v ) restrict(amp,cpu);

// multiply and divide by scalar
static VECTOR3_T operator * ( T f, const VECTOR3_T& v ) restrict(amp,cpu);
static VECTOR3_T operator * ( const VECTOR3_T& v, T f ) restrict(amp,cpu);
static VECTOR3_T operator / ( const VECTOR3_T& v, T f ) restrict(amp,cpu);
static bool operator == ( const VECTOR3_T& v0, const VECTOR3_T& v1 ) restrict(amp,cpu);
static bool operator != ( const VECTOR3_T& v0, const VECTOR3_T& v1 ) restrict(amp,cpu);

class VECTOR3_T
{
public:

	static const VECTOR3_T ZERO;
	static const VECTOR3_T UP;
	static const VECTOR3_T RIGHT;
	static const VECTOR3_T FORWARD;

	union
	{
		struct{T x, y, z;};
		T m_elements[ 3 ];
	};

public:

	VECTOR3_T(T f=(T)0) restrict(amp,cpu)
	{
		m_elements[0] = f;
		m_elements[1] = f;
		m_elements[2] = f;
	}
	VECTOR3_T(T* v) restrict(amp,cpu)
	{
		m_elements[0] = v[0];
		m_elements[1] = v[1];
		m_elements[2] = v[2];
	}
	VECTOR3_T(T x, T y, T z) restrict(amp,cpu)
	{
		m_elements[0] = x;
		m_elements[1] = y;
		m_elements[2] = z;
	}
	VECTOR3_T(const VECTOR3_T& v) restrict(amp,cpu)
	{
		m_elements[0] = v.x;
		m_elements[1] = v.y;
		m_elements[2] = v.z;
	}

	// assignment operators
	VECTOR3_T& operator = ( const VECTOR3_T& rv ) restrict(amp,cpu)
	{
		if( this != &rv )
		{
			m_elements[0] = rv[0];
			m_elements[1] = rv[1];
			m_elements[2] = rv[2];
		}
		return *this;
	}

	// returns the ith element
	const T& operator [] ( int i ) const restrict(amp,cpu)
	{
		return m_elements[i];
	}

    T& operator [] ( int i ) restrict(amp,cpu)
	{
		return m_elements[i];
	}

	T Magnitude() const restrict(amp,cpu)  
	{
		return concurrency::fast_math::sqrt( m_elements[0] * m_elements[0] + m_elements[1] * m_elements[1] + m_elements[2] * m_elements[2] );
	}

    T MagnitudeSquared() const restrict(amp,cpu)
	{
		return ( m_elements[0] * m_elements[0] + m_elements[1] * m_elements[1] + m_elements[2] * m_elements[2] );
	}

	void Normalize() restrict(amp,cpu)
	{
		T norm = Magnitude();
//		assert(norm != (T)0);
		m_elements[0] /= norm;
		m_elements[1] /= norm;
		m_elements[2] /= norm;
	}

	VECTOR3_T Normalized() const restrict(amp,cpu)
	{
		T norm = Magnitude();
//		assert(norm != (T)0);
		return VECTOR3_T(m_elements[0] / norm, m_elements[1] / norm, m_elements[2] / norm);
	}

	void Negate() restrict(amp,cpu)
	{
		m_elements[0] = -m_elements[0];
		m_elements[1] = -m_elements[1];
		m_elements[2] = -m_elements[2];
	}

	// Utility
	operator const T* () const restrict(amp,cpu) // automatic type conversion for OpenGL
	{
		return m_elements;
	}

	operator T* () restrict(amp,cpu) // automatic type conversion for OpenGL 
	{
		return m_elements;
	}

	VECTOR3_T& operator += ( const VECTOR3_T& v ) restrict(amp,cpu)
	{
		m_elements[ 0 ] += v.m_elements[ 0 ];
		m_elements[ 1 ] += v.m_elements[ 1 ];
		m_elements[ 2 ] += v.m_elements[ 2 ];
		return *this;
	}

	VECTOR3_T& operator -= ( const VECTOR3_T& v ) restrict(amp,cpu)
	{
		m_elements[ 0 ] -= v.m_elements[ 0 ];
		m_elements[ 1 ] -= v.m_elements[ 1 ];
		m_elements[ 2 ] -= v.m_elements[ 2 ];
		return *this;
	}

	VECTOR3_T& operator *= ( T f ) restrict(amp,cpu)
	{
		m_elements[ 0 ] *= f;
		m_elements[ 1 ] *= f;
		m_elements[ 2 ] *= f;
		return *this;
	}

    static T DotProduct( const VECTOR3_T& v0, const VECTOR3_T& v1 ) restrict(amp,cpu)
	{
		return v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2];
	}

	static VECTOR3_T CrossProduct( const VECTOR3_T& v0, const VECTOR3_T& v1 ) restrict(amp,cpu)
	{
		return VECTOR3_T
			(
				v0.y * v1.z - v0.z * v1.y,
				v0.z * v1.x - v0.x * v1.z,
				v0.x * v1.y - v0.y * v1.x
			);
	}
    
    // computes the linear interpolation between v0 and v1 by alpha \in [0,1]
	// returns v0 * ( 1 - alpha ) * v1 * alpha
	static VECTOR3_T Lerp( const VECTOR3_T& v0, const VECTOR3_T& v1, T alpha ) restrict(amp,cpu)
	{
		return alpha * ( v1 - v0 ) + v0;
	}

	// computes the cubic catmull-rom interpolation between p0, p1, p2, p3
    // by t \in [0,1].  Guarantees that at t = 0, the result is p0 and
    // at p1, the result is p2.
	static VECTOR3_T CubicInterpolate( const VECTOR3_T& p0, const VECTOR3_T& p1, const VECTOR3_T& p2, const VECTOR3_T& p3, T t ) restrict(amp,cpu)
	{
		// geometric construction:
		//            t
		//   (t+1)/2     t/2
		// t+1        t	        t-1

		// bottom level
		VECTOR3_T p0p1 = VECTOR3_T::Lerp( p0, p1, t + 1 );
		VECTOR3_T p1p2 = VECTOR3_T::Lerp( p1, p2, t );
		VECTOR3_T p2p3 = VECTOR3_T::Lerp( p2, p3, t - 1 );

		// middle level
		VECTOR3_T p0p1_p1p2 = VECTOR3_T::Lerp( p0p1, p1p2, 0.5f * ( t + 1 ) );
		VECTOR3_T p1p2_p2p3 = VECTOR3_T::Lerp( p1p2, p2p3, 0.5f * t );

		// top level
		return VECTOR3_T::Lerp( p0p1_p1p2, p1p2_p2p3, t );
	}
};



VECTOR3_T operator + ( const VECTOR3_T& v0, const VECTOR3_T& v1 ) restrict(amp,cpu) 
{
    return VECTOR3_T( v0[0] + v1[0], v0[1] + v1[1], v0[2] + v1[2] );
}

VECTOR3_T operator - ( const VECTOR3_T& v0, const VECTOR3_T& v1 ) restrict(amp,cpu)  
{
    return VECTOR3_T( v0[0] - v1[0], v0[1] - v1[1], v0[2] - v1[2] );
}

VECTOR3_T operator * ( const VECTOR3_T& v0, const VECTOR3_T& v1 ) restrict(amp,cpu)  
{
    return VECTOR3_T( v0[0] * v1[0], v0[1] * v1[1], v0[2] * v1[2] );
}

VECTOR3_T operator / ( const VECTOR3_T& v0, const VECTOR3_T& v1 ) restrict(amp,cpu)  
{
    return VECTOR3_T( v0[0] / v1[0], v0[1] / v1[1], v0[2] / v1[2] );
}

VECTOR3_T operator - ( const VECTOR3_T& v ) restrict(amp,cpu)  
{
    return VECTOR3_T( -v[0], -v[1], -v[2] );
}

VECTOR3_T operator * ( T f, const VECTOR3_T& v ) restrict(amp,cpu)  
{
    return VECTOR3_T( v[0] * f, v[1] * f, v[2] * f );
}

VECTOR3_T operator * ( const VECTOR3_T& v, T f ) restrict(amp,cpu)  
{
    return VECTOR3_T( v[0] * f, v[1] * f, v[2] * f );
}

VECTOR3_T operator / ( const VECTOR3_T& v, T f ) restrict(amp,cpu)  
{
    return VECTOR3_T( v[0] / f, v[1] / f, v[2] / f );
}

bool operator == ( const VECTOR3_T& v0, const VECTOR3_T& v1 ) restrict(amp,cpu)  
{
    return( v0.x == v1.x && v0.y == v1.y && v0.z == v1.z );
}

bool operator != ( const VECTOR3_T& v0, const VECTOR3_T& v1 ) restrict(amp,cpu)  
{
    return !( v0 == v1 );
}

typedef VECTOR3_T Vector3T;
typedef VECTOR3_T Vec3T;