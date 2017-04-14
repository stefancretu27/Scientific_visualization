#include <rfftw.h>              //the numerical simulation FFTW 

#ifndef VIS_H
#define VIS_H
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

#define MAX_SEEDS 10

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

//glyphs types
typedef enum
{
	MATTER,
	GLYPHS,
	STREAMLINES
}visualization_technique;

/*
 * Structs
 */
typedef struct int_coord int_coord;
typedef struct point point;  
typedef struct vector vector; 
typedef struct grid_cell grid_cell;

struct int_coord
{
	int a, b;
};

struct point
{
	float r, s;
};

struct vector
{
	//vector's values on Ox and Oy axes
	fftw_real x, y;
};

struct grid_cell
{
	//cel dimensions
	fftw_real width, height;
	//bottom left corner
	fftw_real x, y;
};

/*
 * Structs related functions
 */ 
void initialize_cell(grid_cell *cell, fftw_real w, fftw_real h, vector iterator, vector step);
void draw_cell(grid_cell *cell); 

/*
 * Helper functions 
 */
void clamp_value_to_01(float *value);
void draw_triangle_vertex( double x, double y, float color_value);
float direction2angle(float x, float y);
fftw_real vector_magnitude(fftw_real vx, fftw_real vy);

/* 
 * Functions implemented for step 2. 
 */
//color mapping functions
void grayscale(float value, float *R, float *G, float *B);
void rainbow(float value, float* R, float* G, float* B);
void green_saturation(float value, float *R, float *G, float *B);
void red_saturation(float value, float *R, float *G, float *B);

//functions used for parameterize the colormap
float set_color_bands(float value, unsigned short int type);
void rgb2hsv(float r, float g , float b, float *h, float *s, float *v);
void hsv2rgb(float h, float s, float v, float *r, float *g, float *b);

//scale the color map
float set_scale(float value, float x0, float x1,  float y0, float y1);
float clamp_value( float value, float min, float max);

//color legends
void draw_matter_color_legend();
void draw_glyph_color_legend(unsigned int winWidth, unsigned int winHeight);

//drawing functions
void draw_matter_fluid_smoke(fftw_real wn, fftw_real hn, fftw_real *rho, fftw_real *vx, fftw_real *vy, fftw_real *fx, fftw_real *fy);

/* 
 * Functions implemented for step 3. 
 */
void draw_hedgehog(grid_cell cell, vector data, float vec_scale);

//Draw 3D glyphs
void draw_cones(grid_cell cell, vector data, float vec_scale);	
void draw_ellipsoids(grid_cell cell, vector data, float vec_scale);

//draw the glyphs base don user's input
void draw_glyphs(grid_cell cell, vector data, float vec_scale);

//set the glyphs' color based on selected scalar field
void glyphs_scalar_color(float density, vector velocity, vector force);

//interpolation
void cell_corners_indexing(int i, int j, int *idx1, int *idx2, int *idx3, int *idx4);
fftw_real bilinear_interpolation(int idx1, int idx2, int idx3, int idx4, point data_point, fftw_real *attribute);

/* 
 * Functions implemented for step 4 
 */
void compute_density_gradient(int i, int j, vector step, point data_point, fftw_real *density, vector *gradient);
void compute_velocity_magnitude_gradient(int i, int j, vector step, point data_point, fftw_real *vx, fftw_real *vy, vector *gradient);

/* 
 * Functions implemented for step 5 
 */
void draw_seed_points(grid_cell cell);
#endif 
