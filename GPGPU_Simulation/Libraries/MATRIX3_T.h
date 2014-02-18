

#pragma once
#include "MATH_DEFINITION.h"
#include "VECTOR3_T.h"

class MATRIX3_T;

static VECTOR3_T operator * ( const MATRIX3_T& m, const VECTOR3_T& v ) restrict(cpu,amp);
static MATRIX3_T operator * ( const MATRIX3_T& x, const MATRIX3_T& y ) restrict(cpu,amp);

class MATRIX3_T
{
public:

    // Fill a 3x3 matrix with "fill", default to 0.
	MATRIX3_T( T fill = 0.f ) restrict(cpu,amp)
	{
		for( int i = 0; i < 9; ++i )
		{
			m_elements[ i ] = fill;
		}
	}

	MATRIX3_T( T m00, T m01, T m02,	   T m10, T m11, T m12,     T m20, T m21, T m22 ) restrict(cpu,amp)
	{
		m_elements[ 0 ] = m00;
		m_elements[ 1 ] = m10;
		m_elements[ 2 ] = m20;

		m_elements[ 3 ] = m01;
		m_elements[ 4 ] = m11;
		m_elements[ 5 ] = m21;

		m_elements[ 6 ] = m02;
		m_elements[ 7 ] = m12;
		m_elements[ 8 ] = m22;
	}

	// setColumns = true ==> sets the columns of the matrix to be [v0 v1 v2]
	// otherwise, sets the rows
	MATRIX3_T( const VECTOR3_T& v0, const VECTOR3_T& v1, const VECTOR3_T& v2, bool setColumns = true ) restrict(cpu,amp)
	{
		if( setColumns )
		{
			SetCol( 0, v0 );
			SetCol( 1, v1 );
			SetCol( 2, v2 );
		}
		else
		{
			SetRow( 0, v0 );
			SetRow( 1, v1 );
			SetRow( 2, v2 );
		}
	}

	MATRIX3_T( const MATRIX3_T& rm ) restrict(cpu,amp) // copy constructor
	{
		for( int i = 0; i < 9; ++i )
		{
			m_elements[ i ] = rm.m_elements[i];
		}
	}
	MATRIX3_T& operator = ( const MATRIX3_T& rm ) restrict(cpu,amp) // assignment operator
	{
		if( this != &rm )
		{
			for( int i = 0; i < 9; ++i )
			{
				m_elements[ i ] = rm.m_elements[i];
			}
		}
		return *this;
	}

	const T& operator () ( int i, int j ) const restrict(cpu,amp)
	{
		return m_elements[ j * 3 + i ];
	}

	T& operator () ( int i, int j ) restrict(cpu,amp)
	{
		return m_elements[ j * 3 + i ];
	}
	
	VECTOR3_T GetRow( int i ) const restrict(cpu,amp)
	{
		return VECTOR3_T
		(
			m_elements[ i ],
			m_elements[ i + 3 ],
			m_elements[ i + 6 ]
		);
	}

	void SetRow( int i, const VECTOR3_T& v ) restrict(cpu,amp)
	{
		m_elements[ i ] = v.x;
		m_elements[ i + 3 ] = v.y;
		m_elements[ i + 6 ] = v.z;
	}

	VECTOR3_T GetCol( int j ) const restrict(cpu,amp)
	{
		int colStart = 3 * j;

		return VECTOR3_T
		(
			m_elements[ colStart ],
			m_elements[ colStart + 1 ],
			m_elements[ colStart + 2 ]			
		);
	}

	void SetCol( int j, const VECTOR3_T& v ) restrict(cpu,amp)
	{
		int colStart = 3 * j;

		m_elements[ colStart ] = v.x;
		m_elements[ colStart + 1 ] = v.y;
		m_elements[ colStart + 2 ] = v.z;
	}

	// Computes the 3x3 matrix determinant
	T Determinant() const restrict(cpu,amp)
	{
		return MATRIX3_T::Determinant3x3
		(
			m_elements[ 0 ], m_elements[ 3 ], m_elements[ 6 ],
			m_elements[ 1 ], m_elements[ 4 ], m_elements[ 7 ],
			m_elements[ 2 ], m_elements[ 5 ], m_elements[ 8 ]
		);
	}

    // Returns the inverse of this matrix.
    // if pbIsSingular is not NULL, then
    // it is set to whether or not this matrix was singular.
    // A matrix is considered singular if its determinant has
    // absolute value less than epsilon.
	MATRIX3_T Inverse( bool* pbIsSingular = 0, T epsilon = (T)0 ) const restrict(cpu,amp)
	{
		T m00 = m_elements[ 0 ];
		T m10 = m_elements[ 1 ];
		T m20 = m_elements[ 2 ];

		T m01 = m_elements[ 3 ];
		T m11 = m_elements[ 4 ];
		T m21 = m_elements[ 5 ];

		T m02 = m_elements[ 6 ];
		T m12 = m_elements[ 7 ];
		T m22 = m_elements[ 8 ];

		T cofactor00 =  determinant2x2( m11, m12, m21, m22 );
		T cofactor01 = -determinant2x2( m10, m12, m20, m22 );
		T cofactor02 =  determinant2x2( m10, m11, m20, m21 );

		T cofactor10 = -determinant2x2( m01, m02, m21, m22 );
		T cofactor11 =  determinant2x2( m00, m02, m20, m22 );
		T cofactor12 = -determinant2x2( m00, m01, m20, m21 );

		T cofactor20 =  determinant2x2( m01, m02, m11, m12 );
		T cofactor21 = -determinant2x2( m00, m02, m10, m12 );
		T cofactor22 =  determinant2x2( m00, m01, m10, m11 );

		T determinant = m00 * cofactor00 + m01 * cofactor01 + m02 * cofactor02;
	
		bool isSingular = ( concurrency::fast_math::fabs( determinant ) < epsilon );
		if( isSingular )
		{
			if( pbIsSingular != NULL )
			{
				*pbIsSingular = true;
			}
			return MATRIX3_T();
		}
		else
		{
			if( pbIsSingular != NULL )
			{
				*pbIsSingular = false;
			}

			T reciprocalDeterminant = 1.0f / determinant;

			return MATRIX3_T
			(
				cofactor00 * reciprocalDeterminant, cofactor10 * reciprocalDeterminant, cofactor20 * reciprocalDeterminant,
				cofactor01 * reciprocalDeterminant, cofactor11 * reciprocalDeterminant, cofactor21 * reciprocalDeterminant,
				cofactor02 * reciprocalDeterminant, cofactor12 * reciprocalDeterminant, cofactor22 * reciprocalDeterminant
			);
		}
	}

    // in-place transpose
	void Transpose() restrict(cpu,amp)
	{
		T temp;
		for( int i = 0; i < 2; ++i )
		{
			for( int j = i + 1; j < 3; ++j )
			{
				temp = ( *this )( i, j );
				( *this )( i, j ) = ( *this )( j, i );
				( *this )( j, i ) = temp;
			}
		}
	}

    // Returns the transpose of this matrix
	MATRIX3_T Transposed() const restrict(cpu,amp)
	{
		MATRIX3_T out;
		for( int i = 0; i < 3; ++i )
		{
			for( int j = 0; j < 3; ++j )
			{
				out( j, i ) = ( *this )( i, j );
			}
		}
		return out;
	}

	// ---- Utility ----
	operator const T* () const restrict(cpu,amp) // automatic type conversion for GL
	{
		return m_elements;
	}
	operator T* () restrict(cpu,amp) // automatic type conversion for GL
	{
		return m_elements;
	}

	static T determinant2x2( T m00, T m01,    T m10, T m11 ) restrict(cpu,amp)
	{
		return( m00 * m11 - m01 * m10 );
	}

	static T Determinant3x3( T m00, T m01, T m02,   T m10, T m11, T m12,   T m20, T m21, T m22 ) restrict(cpu,amp)
	{
		return
			(
				  m00 * ( m11 * m22 - m12 * m21 )
				- m01 * ( m10 * m22 - m12 * m20 )
				+ m02 * ( m10 * m21 - m11 * m20 )
			);
	}

    // a 3x3 matrix of all ones
	static MATRIX3_T Ones() restrict(cpu,amp)
	{
		MATRIX3_T m;
		for( int i = 0; i < 9; ++i )
		{
			m.m_elements[ i ] = 1;
		}
		return m;
	}

    // a 3x3 identity matrix
	static MATRIX3_T Identity() restrict(cpu,amp)
	{
		MATRIX3_T m;

		m( 0, 0 ) = 1;
		m( 1, 1 ) = 1;
		m( 2, 2 ) = 1;

		return m;
	}

    // Returns the rotation matrix for a rotation of
    // "radians" radians, counterclockwise about the x-axis.
	static MATRIX3_T RotateX( T radians ) restrict(cpu,amp)
	{
		T c = concurrency::fast_math::cos( radians );
		T s = concurrency::fast_math::sin( radians );

		return MATRIX3_T
		(
			1, 0, 0,
			0, c, -s,
			0, s, c
		);
	}
	
    // Returns the rotation matrix for a rotation of
    // "radians" radians, counterclockwise about the y-axis.
	static MATRIX3_T RotateY( T radians ) restrict(cpu,amp)
	{
		T c = concurrency::fast_math::cos( radians );
		T s = concurrency::fast_math::sin( radians );

		return MATRIX3_T
		(
			c, 0, s,
			0, 1, 0,
			-s, 0, c
		);
	}

	// Returns the rotation matrix for a rotation of
    // "radians" radians, counterclockwise about the z-axis.
    static MATRIX3_T RotateZ( T radians ) restrict(cpu,amp)
	{
		T c = concurrency::fast_math::cos( radians );
		T s = concurrency::fast_math::sin( radians );

		return MATRIX3_T
		(
			c, -s, 0,
			s, c, 0,
			0, 0, 1
		);
	}

    // Returns a matrix that non-uniformsly scales
    // the three axes by sx, sy, and sz.
	static MATRIX3_T Scaling( T sx, T sy, T sz ) restrict(cpu,amp)
	{
		return MATRIX3_T
		(
			sx, 0, 0,
			0, sy, 0,
			0, 0, sz
		);
	}

    // Returns a matrix that uniformly scales
    // each axis by "s".
	static MATRIX3_T UniformScaling( T s ) restrict(cpu,amp)
	{
		return MATRIX3_T
		(
			s, 0, 0,
			0, s, 0,
			0, 0, s
		);
	}

    // Returns the rotation matrix for a rotation of
    // "radians" radians, counterclockwise
    // about the axis "rDirection"
	static MATRIX3_T Rotation( const VECTOR3_T& rDirection, T radians ) restrict(cpu,amp)
	{
		VECTOR3_T normalizedDirection = rDirection.Normalized();
	
		T cosTheta = concurrency::fast_math::cos( radians );
		T sinTheta = concurrency::fast_math::sin( radians );

		T x = normalizedDirection.x;
		T y = normalizedDirection.y;
		T z = normalizedDirection.z;

		return MATRIX3_T
			(
				x * x * ( 1.0f - cosTheta ) + cosTheta,			y * x * ( 1.0f - cosTheta ) - z * sinTheta,		z * x * ( 1.0f - cosTheta ) + y * sinTheta,
				x * y * ( 1.0f - cosTheta ) + z * sinTheta,		y * y * ( 1.0f - cosTheta ) + cosTheta,			z * y * ( 1.0f - cosTheta ) - x * sinTheta,
				x * z * ( 1.0f - cosTheta ) - y * sinTheta,		y * z * ( 1.0f - cosTheta ) + x * sinTheta,		z * z * ( 1.0f - cosTheta ) + cosTheta
			);
	}

	MATRIX3_T& operator *= ( T f ) restrict(amp,cpu)
	{
		for( int i = 0; i < 9; ++i )
		{
			m_elements[i] *= f;
		}
		return *this;
	}

public:

	T m_elements[ 9 ];

};


// Matrix-Vector multiplication
// 3x3 * 3x1 ==> 3x1
static VECTOR3_T operator * (const MATRIX3_T& m, const VECTOR3_T& v) restrict(cpu, amp)
{
	VECTOR3_T output( 0, 0, 0 );

	for( int i = 0; i < 3; ++i )
	{
		for( int j = 0; j < 3; ++j )
		{
			output[ i ] += m( i, j ) * v[ j ];
		}
	}

	return output;
}

// Matrix-Matrix multiplication
static MATRIX3_T operator * (const MATRIX3_T& x, const MATRIX3_T& y) restrict(cpu, amp)
{
	MATRIX3_T product; // zeroes

	for( int i = 0; i < 3; ++i )
	{
		for( int j = 0; j < 3; ++j )
		{
			for( int k = 0; k < 3; ++k )
			{
				product( i, k ) += x( i, j ) * y( j, k );
			}
		}
	}

	return product;
}

static MATRIX3_T operator * (const MATRIX3_T& x, const T& s) restrict(cpu, amp)
{
	MATRIX3_T product; // zeroes
	for( int i = 0; i < 9; ++i )
	{
		product.m_elements[i] = x.m_elements[i]*s;
	}

	return product;
}

typedef MATRIX3_T Matrix3T;
typedef MATRIX3_T Mat3T;