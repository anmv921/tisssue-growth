#ifndef V_DEFS

#define V_DEFS

typedef struct {
	real x, y;
} VecR;

#define VZero(v) \
	(v).x = 0., \
	(v).y = 0.
	
#define VProd(v, c) \
	v.x *= c, \
	v.y *= c
	
#define Dist(x, y) \
	sqrt(x*x + y*y)
	
#define DotProd(v1, v2) \
	(v1).x * (v2).x + (v1).y * (v2).y

#define GetDr(DR, R1, R2) \
	DR.x = R1.x - R2.x, \
	DR.y = R1.y - R2.y

#define AddForce(f0, ffx, ffy) \
	(f0).x += ffx, \
	(f0).y += ffy

#endif

