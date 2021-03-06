

#pragma once

#include <stdlib.h>
#include <math.h>
#include <amp_math.h>

#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

#define M_PI 3.14159265358979323846

#define ABS(a) ((a) > 0 ? (a) : -(a))

#define MIN(a, b)							((a) > (b) ? (b) : (a))
#define MIN3(a, b, c)						(MIN(MIN(a, b), (c)))
#define MIN4(a, b, c, d)					(MIN(MIN3(a, b, c), (d)))
#define MIN5(a, b, c, d, e)					(MIN(MIN4(a, b, c, d), (e)))
#define MIN7(a, b, c, d, e, f, g)			(MIN(MIN3(a, b, c), MIN4(d, e, f, g)))
#define MIN8(a, b, c, d, e, f, g, h)		(MIN(MIN7(a, b, c, d, e, f, g), h))

#define MAX(a, b)							((a) > (b) ? (a) : (b))
#define MAX3(a, b, c)						(MAX(MAX(a, b), (c)))
#define MAX4(a, b, c, d)					(MAX(MAX3(a, b, c), (d)))
#define MAX5(a, b, c, d, e)					(MAX(MAX4(a, b, c, d), (e)))
#define MAX7(a, b, c, d, e, f, g)			(MAX(MAX3(a, b, c), MAX4(d, e, f, g)))
#define MAX8(a, b, c, d, e, f, g, h)		(MAX(MAX7(a, b, c, d, e, f, g), h))

#define MIN_ABS2(a, b)						(ABS(a) > ABS(b) ? (b) : (a))
#define MAX_ABS2(a, b)						(ABS(a) > ABS(b) ? (a) : (b))

#define CLAMP(v, min, max)		((v) > (max) ? (max) : ((v) < (min) ? (min) : (v)))

#define SQUARE(a)				((a)*(a))
#define POW2(a)					((a)*(a))
#define POW3(a)					((a)*(a)*(a))

#define SIGN(a) ((a > 0) ? 1 : -1)


#define SWAP(a, b, temp) temp=a; a=b; b=temp;

inline float POW4(const float& a){const float a2 = a*a; return a2*a2;}
inline float POW5(const float& a){const float a2 = a*a; return a2*a2*a;}
inline float POW6(const float& a){const float a3 = a*a*a; return a3*a3;}
inline float POW7(const float& a){const float a3 = a*a*a; return a3*a3*a;}
inline float POW8(const float& a){const float a2 = a*a; const float a4 = a2*a2; return a4*a4;}

#define TS  float // or double
#define INT int

typedef glm::uvec2 TU2;
typedef glm::uvec3 TU3;
typedef glm::uvec4 TU4;

typedef glm::vec2 TV2;
typedef glm::vec3 TV3;
typedef glm::vec4 TV4;

typedef glm::mat2 TM2;
typedef glm::mat3 TM3;
typedef glm::mat4 TM4;

typedef glm::quat TQ;