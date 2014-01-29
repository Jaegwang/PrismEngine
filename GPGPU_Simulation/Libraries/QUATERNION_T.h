
#pragma once
#include "MATH_DEFINITION.h"
#include "VECTOR3_T.h"

class QUATERNION_T;

static QUATERNION_T operator + ( const QUATERNION_T& q0, const QUATERNION_T& q1 ) restrict(cpu,amp);
static QUATERNION_T operator - ( const QUATERNION_T& q0, const QUATERNION_T& q1 ) restrict(cpu,amp);
static QUATERNION_T operator * ( const QUATERNION_T& q0, const QUATERNION_T& q1 ) restrict(cpu,amp);
static QUATERNION_T operator * ( T f, const QUATERNION_T& q ) restrict(cpu,amp);
static QUATERNION_T operator * ( const QUATERNION_T& q, T f ) restrict(cpu,amp);


class QUATERNION_T
{
public:

	static const QUATERNION_T ZERO;
	static const QUATERNION_T IDENTITY;

public:

	QUATERNION_T() restrict(cpu,amp)
	{
		m_elements[ 0 ] = 0;
		m_elements[ 1 ] = 0;
		m_elements[ 2 ] = 0;
		m_elements[ 3 ] = 0;
	}

	// q = w + x * i + y * j + z * k
	QUATERNION_T( T w, T x, T y, T z ) restrict(cpu,amp)
	{
		m_elements[ 0 ] = w;
		m_elements[ 1 ] = x;
		m_elements[ 2 ] = y;
		m_elements[ 3 ] = z;
	}
		
	QUATERNION_T( const QUATERNION_T& rq ) restrict(cpu,amp) // copy constructor
	{
		m_elements[ 0 ] = rq.m_elements[ 0 ];
		m_elements[ 1 ] = rq.m_elements[ 1 ];
		m_elements[ 2 ] = rq.m_elements[ 2 ];
		m_elements[ 3 ] = rq.m_elements[ 3 ];
	}
	QUATERNION_T& operator = ( const QUATERNION_T& rq ) restrict(cpu,amp) // assignment operator	
	{
		if( this != ( &rq ) )
		{
			m_elements[ 0 ] = rq.m_elements[ 0 ];
			m_elements[ 1 ] = rq.m_elements[ 1 ];
			m_elements[ 2 ] = rq.m_elements[ 2 ];
			m_elements[ 3 ] = rq.m_elements[ 3 ];
		}
		return( *this );
	}

	// returns a quaternion with 0 real part
	QUATERNION_T( const VECTOR3_T& v ) restrict(cpu,amp)
	{
		m_elements[ 0 ] = 0;
		m_elements[ 1 ] = v[ 0 ];
		m_elements[ 2 ] = v[ 1 ];
		m_elements[ 3 ] = v[ 2 ];
	}

	// returns the ith element
	const T & operator [] ( int i ) const restrict(cpu,amp)
	{
		return m_elements[ i ];
	}

	T & operator [] ( int i ) restrict(cpu,amp)
	{
		return m_elements[ i ];
	}

	VECTOR3_T XYZ() const restrict(cpu,amp)
	{
		return VECTOR3_T(m_elements[ 1 ], m_elements[ 2 ], m_elements[ 3 ]);
	}

	T Magnitude() const restrict(cpu,amp)
	{
		return concurrency::fast_math::sqrt(MagnitudeSquared());
	}

	T MagnitudeSquared() const restrict(cpu,amp)
	{
		return
		(
			m_elements[ 0 ] * m_elements[ 0 ] +
			m_elements[ 1 ] * m_elements[ 1 ] +
			m_elements[ 2 ] * m_elements[ 2 ] +
			m_elements[ 3 ] * m_elements[ 3 ]
		);
	}
	
	void Normalize() restrict(cpu,amp)
	{
		T reciprocalAbs = 1.f / Magnitude();

		m_elements[ 0 ] *= reciprocalAbs;
		m_elements[ 1 ] *= reciprocalAbs;
		m_elements[ 2 ] *= reciprocalAbs;
		m_elements[ 3 ] *= reciprocalAbs;
	}

	QUATERNION_T Normalized() const restrict(cpu,amp)
	{
		QUATERNION_T q( *this );
		q.Normalize();
		return q;
	}

	void Conjugate() restrict(cpu,amp)
	{
		m_elements[ 1 ] = -m_elements[ 1 ];
		m_elements[ 2 ] = -m_elements[ 2 ];
		m_elements[ 3 ] = -m_elements[ 3 ];
	}

	QUATERNION_T Conjugated() const restrict(cpu,amp)
	{
		return QUATERNION_T
		(
			 m_elements[ 0 ],
			-m_elements[ 1 ],
			-m_elements[ 2 ],
			-m_elements[ 3 ]
		);
	}

	void Invert() restrict(cpu,amp)
	{
		QUATERNION_T inverse = Conjugated() * ( (T)1 / MagnitudeSquared() );

		m_elements[ 0 ] = inverse.m_elements[ 0 ];
		m_elements[ 1 ] = inverse.m_elements[ 1 ];
		m_elements[ 2 ] = inverse.m_elements[ 2 ];
		m_elements[ 3 ] = inverse.m_elements[ 3 ];
	}

	QUATERNION_T Inverse() const restrict(cpu,amp)
	{
		return Conjugated() * ( (T)1 / MagnitudeSquared() );
	}

	// log and exponential maps
	QUATERNION_T Log() const restrict(cpu,amp)
	{
		T len =
			concurrency::fast_math::sqrt
			(
				m_elements[ 1 ] * m_elements[ 1 ] +
				m_elements[ 2 ] * m_elements[ 2 ] +
				m_elements[ 3 ] * m_elements[ 3 ]
			);

		if( len < 1e-6 )
		{
			return QUATERNION_T( 0, m_elements[ 1 ], m_elements[ 2 ], m_elements[ 3 ] );
		}
		else
		{
			T coeff = concurrency::fast_math::acos( m_elements[ 0 ] ) / len;
			return QUATERNION_T( 0, m_elements[ 1 ] * coeff, m_elements[ 2 ] * coeff, m_elements[ 3 ] * coeff );
		}
	}

	QUATERNION_T Exp() const restrict(cpu,amp)
	{
		T theta =
			concurrency::fast_math::sqrt
			(
				m_elements[ 1 ] * m_elements[ 1 ] +
				m_elements[ 2 ] * m_elements[ 2 ] +
				m_elements[ 3 ] * m_elements[ 3 ]
			);

		if( theta < 1e-6 )
		{
			return QUATERNION_T( concurrency::fast_math::cos( theta ), m_elements[ 1 ], m_elements[ 2 ], m_elements[ 3 ] );
		}
		else
		{
			T coeff = concurrency::fast_math::sin( theta ) / theta;
			return QUATERNION_T( concurrency::fast_math::cos( theta ), m_elements[ 1 ] * coeff, m_elements[ 2 ] * coeff, m_elements[ 3 ] * coeff );		
		}
	}
	
	// returns unit vector for rotation and radians about the unit vector
	VECTOR3_T GetAxisAngle( T* radiansOut ) restrict(cpu,amp)
	{
		T theta = concurrency::fast_math::acos( w ) * 2;
		T vectorNorm = concurrency::fast_math::sqrt( x * x + y * y + z * z );
		T reciprocalVectorNorm = 1.f / vectorNorm;

		*radiansOut = theta;
		return VECTOR3_T
		(
			x * reciprocalVectorNorm,
			y * reciprocalVectorNorm,
			z * reciprocalVectorNorm
		);
	}

	// sets this quaternion to be a rotation of fRadians about v = < fx, fy, fz >, v need not necessarily be unit length
	void SetAxisAngle( T radians, const VECTOR3_T& axis ) restrict(cpu,amp)
	{
		m_elements[ 0 ] = concurrency::fast_math::cos( radians / 2 );

		T sinHalfTheta = concurrency::fast_math::sin( radians / 2 );
		T vectorNorm = axis.Magnitude();
		T reciprocalVectorNorm = 1.f / vectorNorm;

		m_elements[ 1 ] = axis.x * sinHalfTheta * reciprocalVectorNorm;
		m_elements[ 2 ] = axis.y * sinHalfTheta * reciprocalVectorNorm;
		m_elements[ 3 ] = axis.z * sinHalfTheta * reciprocalVectorNorm;
	}

	void SetAngularVector( VECTOR3_T& v ) restrict(cpu,amp)
	{
		const T magnitude = v.Magnitude();
		if(magnitude <= (T)1e-8)
		{
			(*this) = QUATERNION_T( 1, 0, 0, 0 );
			return;
		}

		m_elements[0] = concurrency::fast_math::cos((T)0.5*magnitude);

		VECTOR3_T axis = ((T)concurrency::fast_math::sin((T)0.5*magnitude)/magnitude)*v;

		m_elements[1] = axis[0]; m_elements[2] = axis[1]; m_elements[3] = axis[2];
	}

	// ---- Utility ----
	// quaternion dot product (a la vector)
	static T DotProduct( const QUATERNION_T& q0, const QUATERNION_T& q1 ) restrict(cpu,amp)
	{
		return
		(
			q0.w * q1.w +
			q0.x * q1.x +
			q0.y * q1.y +
			q0.z * q1.z
		);
	}
	
	// linear (stupid) interpolation
	static QUATERNION_T Lerp( const QUATERNION_T& q0, const QUATERNION_T& q1, T alpha ) restrict(cpu,amp)
	{
		return( ( q0 + alpha * ( q1 - q0 ) ).Normalized() );
	}

	// spherical linear interpolation
	static QUATERNION_T Slerp( const QUATERNION_T& a, const QUATERNION_T& b, T t, bool allowFlip = true ) restrict(cpu,amp)
	{
		T cosAngle = QUATERNION_T::DotProduct( a, b );

		T c1;
		T c2;

		// Linear interpolation for close orientations
		if( ( 1.0f - concurrency::fast_math::fabs( cosAngle ) ) < 0.01f )
		{
			c1 = 1.0f - t;
			c2 = t;
		}
		else
		{
			// Spherical interpolation
			T angle = concurrency::fast_math::acos( concurrency::fast_math::fabs( cosAngle ) );
			T sinAngle = concurrency::fast_math::sin( angle );
			c1 = concurrency::fast_math::sin( angle * ( 1.0f - t ) ) / sinAngle;
			c2 = concurrency::fast_math::sin( angle * t ) / sinAngle;
		}

		// Use the shortest path
		if( allowFlip && ( cosAngle < 0.0f ) )
		{
			c1 = -c1;
		}

		return QUATERNION_T( c1 * a[ 0 ] + c2 * b[ 0 ], c1 * a[ 1 ] + c2 * b[ 1 ], c1 * a[ 2 ] + c2 * b[ 2 ], c1 * a[ 3 ] + c2 * b[ 3 ] );
	}

	// spherical quadratic interoplation between a and b at point t
	// given quaternion tangents tanA and tanB (can be computed using squadTangent)	
	static QUATERNION_T Squad( const QUATERNION_T& a, const QUATERNION_T& tanA, const QUATERNION_T& tanB, const QUATERNION_T& b, T t ) restrict(cpu,amp)
	{
		QUATERNION_T ab = QUATERNION_T::Slerp( a, b, t );
		QUATERNION_T tangent = QUATERNION_T::Slerp( tanA, tanB, t, false );
		return QUATERNION_T::Slerp( ab, tangent, 2.0f * t * ( 1.0f - t ), false );
	}

	static QUATERNION_T CubicInterpolate( const QUATERNION_T& q0, const QUATERNION_T& q1, const QUATERNION_T& q2, const QUATERNION_T& q3, T t ) restrict(cpu,amp)
	{
		// geometric construction:
		//            t
		//   (t+1)/2     t/2
		// t+1        t	        t-1

		// bottom level
		QUATERNION_T q0q1 = QUATERNION_T::Slerp( q0, q1, t + 1 );
		QUATERNION_T q1q2 = QUATERNION_T::Slerp( q1, q2, t );
		QUATERNION_T q2q3 = QUATERNION_T::Slerp( q2, q3, t - 1 );

		// middle level
		QUATERNION_T q0q1_q1q2 = QUATERNION_T::Slerp( q0q1, q1q2, 0.5f * ( t + 1 ) );
		QUATERNION_T q1q2_q2q3 = QUATERNION_T::Slerp( q1q2, q2q3, 0.5f * t );

		// top level
		return QUATERNION_T::Slerp( q0q1_q1q2, q1q2_q2q3, t );
	}

	// Log-difference between a and b, used for squadTangent
	// returns log( a^-1 b )	
	static QUATERNION_T LogDifference( const QUATERNION_T& a, const QUATERNION_T& b ) restrict(cpu,amp)
	{
		QUATERNION_T diff = a.Inverse() * b;
		diff.Normalize();
		return diff.Log();
	}

	// Computes a tangent at center, defined by the before and after quaternions
	// Useful for squad()
	static QUATERNION_T SquadTangent( const QUATERNION_T& before, const QUATERNION_T& center, const QUATERNION_T& after ) restrict(cpu,amp)
	{
		QUATERNION_T l1 = QUATERNION_T::LogDifference( center, before );
		QUATERNION_T l2 = QUATERNION_T::LogDifference( center, after );
	
		QUATERNION_T e;
		for( int i = 0; i < 4; ++i )
		{
			e[ i ] = -0.25f * ( l1[ i ] + l2[ i ] );
		}
		e = center * ( e.Exp() );

		return e;
	}

	static QUATERNION_T FromAngularVector( const VECTOR3_T& v ) restrict(cpu,amp)
	{
		const T magnitude = v.Magnitude();
		if(magnitude <= (T)1e-8) return QUATERNION_T( 1, 0, 0, 0 );

		QUATERNION_T q;
		q.m_elements[0] = concurrency::fast_math::cos((T)0.5*magnitude);

		VECTOR3_T axis = ((T)concurrency::fast_math::sin((T)0.5*magnitude)/magnitude)*v;

		q.m_elements[1] = axis[0]; q.m_elements[2] = axis[1]; q.m_elements[3] = axis[2];

		return q;
	}

	// returns a unit quaternion that's a uniformly distributed rotation
	// given u[i] is a uniformly distributed random number in [0,1]
	// taken from Graphics Gems II
	static QUATERNION_T RandomRotation( T u0, T u1, T u2 ) restrict(cpu,amp)
	{
		T z = u0;
		T theta = static_cast< T >( 2.f * M_PI * u1 );
		T r = concurrency::fast_math::sqrt( 1.f - z * z );
		T w = static_cast< T >( M_PI * u2 );

		return QUATERNION_T
		(
			concurrency::fast_math::cos( w ),
			concurrency::fast_math::sin( w ) * concurrency::fast_math::cos( theta ) * r,
			concurrency::fast_math::sin( w ) * concurrency::fast_math::sin( theta ) * r,
			concurrency::fast_math::sin( w ) * z
		);
	}

public:

	union
	{
		struct{T w, x, y, z;};
		T m_elements[ 4 ];
	};
};


// static
const QUATERNION_T QUATERNION_T::ZERO = QUATERNION_T( 0, 0, 0, 0 );
const QUATERNION_T QUATERNION_T::IDENTITY = QUATERNION_T( 1, 0, 0, 0 );


static QUATERNION_T operator + (const QUATERNION_T& q0, const QUATERNION_T& q1) restrict(cpu, amp)
{
	return QUATERNION_T
	(
		q0.w + q1.w,
		q0.x + q1.x,
		q0.y + q1.y,
		q0.z + q1.z
	);
}

static QUATERNION_T operator - (const QUATERNION_T& q0, const QUATERNION_T& q1) restrict(cpu, amp)
{
	return QUATERNION_T
	(
		q0.w - q1.w,
		q0.x - q1.x,
		q0.y - q1.y,
		q0.z - q1.z
	);
}

static QUATERNION_T operator * (const QUATERNION_T& q0, const QUATERNION_T& q1) restrict(cpu, amp)
{
	return QUATERNION_T
	(
		q0.w * q1.w - q0.x * q1.x - q0.y * q1.y - q0.z * q1.z,
		q0.w * q1.x + q0.x * q1.w + q0.y * q1.z - q0.z * q1.y,
		q0.w * q1.y - q0.x * q1.z + q0.y * q1.w + q0.z * q1.x,
		q0.w * q1.z + q0.x * q1.y - q0.y * q1.x + q0.z * q1.w
	);
}

static QUATERNION_T operator * (T f, const QUATERNION_T& q) restrict(cpu, amp)
{
	return QUATERNION_T
	(
		f * q.w,
		f * q.x,
		f * q.y,
		f * q.z
	);
}

static QUATERNION_T operator * (const QUATERNION_T& q, T f) restrict(cpu, amp)
{
	return QUATERNION_T
	(
		f * q.w,
		f * q.x,
		f * q.y,
		f * q.z
	);
}

typedef QUATERNION_T QuaternionT;
typedef QUATERNION_T QuaterT;