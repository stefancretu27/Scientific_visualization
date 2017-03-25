/* 
 * Defines
 */
#ifndef GLOBALS
#define GLOBALS

#define DIM 50
#define MAX_COLOR_BANDS 256
#define max(x, y) x > y ? x : y
#define min(x, y) x < y ? x : y
#define vec_magnitude(x, y) sqrt(x*x + y*y)
#define V_SCALE_FACTOR 50
#define F_SCALE_FACTOR 100
//useful for drawing legend or colors, especially in the case of glyphs
#define MATTER_TYPE 0												
#define GLYPH_TYPE 1
#define LEGEND_TYPE 2 
//used for color bands 
#define SCALAR 0
#define VECTORIAL 1
#define GRADIENT_SCALE_FACTOR 5

#endif

/*
 * Enums 
 */
typedef enum
{
	FALSE,
	TRUE
}bool; 

//types of colormaps
typedef enum 
{
	BLACKWHITE,
	RAINBOW,
	SATURATION
}colormap;

//matter attributes
typedef enum
{
	DENSITY,													
	VELOCITY_MAGNITUDE,
	FORCE_MAGNITUDE
}matter_attribute;

//glyphs attributes
typedef enum
{
	VELOCITY,
	FORCE,
	DENSITY_GRADIENT,
	VEL_MAGN_GRADIENT
}glyphs_attribute;

//glyphs types
typedef enum
{
	HEDGEHOGS,
	CONES,
	ELLIPSES
}glyphs_type;
