
#pragma once

#include "VECTOR3_T.h"
#include "MATRIX3_T.h"
#include "QUATERNION_T.h"

// Returns the rotation matrix represented by a unit quaternion
// if q is not normalized, it is normalized first.
static MATRIX3_T MatrixFromQuaternion( const QUATERNION_T& rq ) restrict(cpu,amp)
{
	QUATERNION_T q = rq.Normalized();

	T xx = q.x * q.x;
	T yy = q.y * q.y;
	T zz = q.z * q.z;

	T xy = q.x * q.y;
	T zw = q.z * q.w;

	T xz = q.x * q.z;
	T yw = q.y * q.w;

	T yz = q.y * q.z;
	T xw = q.x * q.w;

	return MATRIX3_T
		(
			1.0f - 2.0f * ( yy + zz ),		2.0f * ( xy - zw ),				2.0f * ( xz + yw ),
			2.0f * ( xy + zw ),				1.0f - 2.0f * ( xx + zz ),		2.0f * ( yz - xw ),
			2.0f * ( xz - yw ),				2.0f * ( yz + xw ),				1.0f - 2.0f * ( xx + yy )
		);
}

static QUATERNION_T QuaternionFromRotationMatrix( const MATRIX3_T& m ) restrict(cpu,amp)
{
	T x;
	T y;
	T z;
	T w;

	// Compute one plus the trace of the matrix
	T onePlusTrace = 1.0f + m( 0, 0 ) + m( 1, 1 ) + m( 2, 2 );

	if( onePlusTrace > 1e-5 )
	{
		// Direct computation
		T s = concurrency::fast_math::sqrt( onePlusTrace ) * 2.0f;
		x = ( m( 2, 1 ) - m( 1, 2 ) ) / s;
		y = ( m( 0, 2 ) - m( 2, 0 ) ) / s;
		z = ( m( 1, 0 ) - m( 0, 1 ) ) / s;
		w = 0.25f * s;
	}
	else
	{
		// Computation depends on major diagonal term
		if( ( m( 0, 0 ) > m( 1, 1 ) ) & ( m( 0, 0 ) > m( 2, 2 ) ) )
		{
			T s = concurrency::fast_math::sqrt( 1.0f + m( 0, 0 ) - m( 1, 1 ) - m( 2, 2 ) ) * 2.0f;
			x = 0.25f * s;
			y = ( m( 0, 1 ) + m( 1, 0 ) ) / s;
			z = ( m( 0, 2 ) + m( 2, 0 ) ) / s;
			w = ( m( 1, 2 ) - m( 2, 1 ) ) / s;
		}
		else if( m( 1, 1 ) > m( 2, 2 ) )
		{
			T s = concurrency::fast_math::sqrt( 1.0f + m( 1, 1 ) - m( 0, 0 ) - m( 2, 2 ) ) * 2.0f;
			x = ( m( 0, 1 ) + m( 1, 0 ) ) / s;
			y = 0.25f * s;
			z = ( m( 1, 2 ) + m( 2, 1 ) ) / s;
			w = ( m( 0, 2 ) - m( 2, 0 ) ) / s;
		}
		else
		{
			T s = concurrency::fast_math::sqrt( 1.0f + m( 2, 2 ) - m( 0, 0 ) - m( 1, 1 ) ) * 2.0f;
			x = ( m( 0, 2 ) + m( 2, 0 ) ) / s;
			y = ( m( 1, 2 ) + m( 2, 1 ) ) / s;
			z = 0.25f * s;
			w = ( m( 0, 1 ) - m( 1, 0 ) ) / s;
		}
	}

	QUATERNION_T q( w, x, y, z );
	return q.Normalized();
}

static QUATERNION_T QuaternionFromRotatedBasis( const VECTOR3_T& x, const VECTOR3_T& y, const VECTOR3_T& z ) restrict(cpu,amp)
{
	return QuaternionFromRotationMatrix( MATRIX3_T( x, y, z ) );
}
