/*
 * Includes 
 */

#include <stdio.h>              //for printing the help text

#include <math.h>               //for various math functions #include
#include <GL/glut.h>            //the GLUT graphics library
#include </usr/include/GL/freeglut.h>

#include <visualization.h>

/*
 * External data: They are declared here but not instantiated. The memory allocation and initialization happens in fluids.c
 */
extern unsigned int winWidth, winHeight;
extern unsigned short int no_glyphs_x, no_glyphs_y; 
 
//MATTER_GLOBALS          												//draw the matter (1) or not (0)
extern matter_attribute show_matter_attribute;							//default shown attribute
extern colormap matter_color_type;        								//default method for smoke attributes coloring
extern unsigned int matter_color_bands;									//default number of color bands
extern bool matter_scale;
extern bool matter_clamp;
extern float matter_clamp_lmin, matter_clamp_lmax, matter_scale_lmax, matter_scale_lmin;

//GLYPH_GLOBALS
extern colormap glyph_color_type;         								//method for force&velocity coloring
extern unsigned int glyph_color_bands;									//default number of color bands
extern glyphs_attribute show_glyph_attribute;
extern glyphs_type show_glyph_type;

//functions used for setting colors. They're defined and implemented in fluids.c
extern void set_colormap(float scalar_value);
extern void direction_to_color(float x, float y, float value, unsigned short int type);

//Parameters used of rocmputing streamlines
float Xmin = -5, Xmax = 5;							//Range of variable x that we sample
float Ymin = -5, Ymax = 5;							//Range of variable y that we sample
int Nx = 50, Ny	= 50;								//Number of samples taken along the x and y axes

/*
 * Struct related functions
 * Inputs: cell: the instance to grid cell structure
 * 		   w, h: width and height of a cell grid, when it is not sampled/resized
 * 		   step_x, step_y: represent how many times should the glyphs be increased on x and y, in order to cover the same area on the screen
 * 		   i, j: cell indeces on the grid. They change everytime the grid is resized and are incremented with step_x/y
 */ 
void initialize_cell(grid_cell *cell, fftw_real w, fftw_real h, vector iterator, vector step)
{
	cell->width = step.x*w;
	cell->height = step.y*h;
	cell->x = iterator.x*w;
	cell->y = iterator.y*h;
}

void draw_cell(grid_cell *cell)
{
    glLineWidth(1);																																					
	//draw rectangle
	glBegin(GL_LINE_LOOP);
	glVertex2f(cell->x, cell->y);
	glVertex2f(cell->x, cell->y + cell->height);
	glVertex2f(cell->x + cell->width, cell->y + cell->height);
	glVertex2f(cell->x + cell->width, cell->y);
	glEnd();

}

/*
 * Helper functions section
 */

fftw_real vector_normalize(vector v)
{
	return (fftw_real) sqrt(v.x*v.x + v.y*v.y);
}

fftw_real vector_magnitude(fftw_real x, fftw_real y)
{
	return (fftw_real) sqrt(x*x + y*y);
}

//Helper function to clamp inputed value to [0;1]
void clamp_value_to_01(float *value)
{
	if (*value < 0)
	{
		*value = 0;
	}
	if (*value > 1)
	{
		*value = 1;
	}
}

//helper function to set color and 2D position for a vertex 
void draw_triangle_vertex( double x, double y, float color_value)
{
	set_colormap(color_value);    
	glVertex2f(x, y);
}

//compute the angle between the vector and the Ox axis
float direction2angle(float x, float y)			
{
	float n = vector_magnitude(x, y); 
	
	//normalize the vector
	if (n<1.0e-6) 
	{
		n=1;
	} 
		
	x/=n; 
	y/=n; 

	return atan2(x,-y) * (180 / M_PI);
}

/* 
 * Colormapping functions used for step 2.
 * rainbow and grayscale are used for both matter and glyphs.
 * green_saturation is used only for matter.
 * red_saturation is used for glyphs direction (direction to color).
 * All 4 functions convert scalar float value to a color accordingly to their interpolation logic. Firstly, the value is clamped to interval [0;1]
 * Inputs: scalar value, pointers to r,g,b values that are the results of interpolation
 * Output: none
 * */

/* Implements a color palette, mapping the scalar 'value' to a rainbow color RGB. Taken form the book, section 5.2
 * It shows a blue to red transition, with other colors in between: low values are encoded to blue, high values are encoded to red. 
 */
void rainbow(float value, float* R, float* G, float* B)
{
	const float dx=0.8;
	
	clamp_value_to_01(&value);
	
	value = (6-2*dx)*value+dx;
	*R = max(0.0,(3-fabs(value-4)-fabs(value-5))/2);
	*G = max(0.0,(4-fabs(value-2)-fabs(value-4))/2);
	*B = max(0.0,(3-fabs(value-1)-fabs(value-2))/2);
}

void grayscale(float value, float *R, float *G, float *B)
{
	*R = *G = *B = value;
}

/* Implements a green-based color palette. Taken form the book sample 5.
 * It shows a black to white transition, with several green hues in between: low values are encoded to black, high values are encoded to white. 
 */
void green_saturation(float value, float *R, float *G, float *B)
{
	float r = 0, g = 1, b = 0;
	
	clamp_value_to_01(&value);
	
	if (value < 0.5)										//value in [0,0.5]: modulate the luminance from black to the base-color.
	{   
		*R = 2*value*r;
		*G = 2*value*g;
		*B = 2*value*b;
	}
	else													//value in [0.5,1]: modulate the saturation from base-color to white.
	{	
		value = 2*(value-0.5);
		*R = (1-value)*r + 1*value;
		*G = (1-value)*g + 1*value;
		*B = (1-value)*b + 1*value;
	}
}

/* Implements a red-based color palette. It is a variation of the colormap function implemented above.
 * It shows a black to white transition, with several red hues in between: low values are encoded to black, high values are encoded to white. 
 */
void red_saturation(float value, float *R, float *G, float *B)
{
	float r = 1, g = 0, b = 0;
	
	clamp_value_to_01(&value);
	
	if (value < 0.5)										//value in [0,0.5]: modulate the luminance from black to the base-color.
	{   
		*R = 2*value*r;
		*G = 2*value*g;
		*B = 2*value*b;
	}
	else													//value in [0.5,1]: modulate the saturation from base-color to white.
	{	
		value = 2*(value-0.5);
		*R = (1-value)*r + 1*value;
		*G = (1-value)*g + 1*value;
		*B = (1-value)*b + 1*value;
	}
}

/* This is like dividing the interval [0;1] (which is considered in every coloring function),
 * in a number of subintervals equal to "bands" and assign a color for each head of every interval
 * Eg: bands = 2 => 3 colors (only red, blue and green)); bands = 4 => 5 colors, one between R and G (yellow) the other one between B and G (bright blue)
 * */
float set_color_bands(float value, unsigned short int type)
{
	unsigned short int bands;
	
	if(type == SCALAR)
	{
		bands = matter_color_bands;
	}
	
	if(type == VECTORIAL)
	{
		bands = glyph_color_bands;
	}
	
	value *= bands;
	value = (int)(value);
	return value/bands;
}

/* Functions used at step 2 for converting from RGB to HSV and viceversa.
 * Their implementation is an adaptaion to C from the C++ code examples in the book section 3.6.3 
 */

//convert from rgb to hsv. Taken from the book, section 3.6.3
void rgb2hsv(float r, float g , float b, float *h, float *s, float *v)
{
	//scale input values in [0;1]
	clamp_value_to_01(&r);
	clamp_value_to_01(&g);
	clamp_value_to_01(&b);
	
	float M1 = max(g, b);
	float M = max(r, M1);
	float m1 = min(g, b);
	float m = min(r, m1);
	float d = M-m;
	*v = M; //value = max( r, g, b)
	*s = (M > 0.00001)?d/M:0; //saturation
	if (*s == 0) 
	{
		*h = 0;
	}//achromatic case , hue=0 by convention
	else
	{
		if(r == M)
		{
			*h = (g - b)/d;
		}
		else if(g == M) 
		{
			*h = 2 + (b - r)/d;
		}
		else
		{
			*h = 4 + (r - g)/d ;
		}
		
		*h /= 6 ;
		if(h < 0) 
		{
			*h += 1;
		}
	}//chromatic case

}

//convert from hsv to rgb. Taken from the book, section 3.6.3
void hsv2rgb(float h, float s, float v, float *r, float *g, float *b)
{
	//scale input values in [0;1]
	clamp_value_to_01(&h);
	clamp_value_to_01(&s);
	clamp_value_to_01(&v);
		
	int hueCase = (int) (h*6);
	float frac = 6*h - hueCase;
	float lx = v*(1 - s);
	float ly = v*(1 - s*frac);
	float lz = v*(1 - s*(1 - frac));
	
	switch(hueCase)
	{
		case 0:
		case 6: *r = v; *g = lz; *b = lx; break; // 0<hue<1/6
		case 1: *r = ly; *g = v ; *b = lx; break; // 1/6<hue<2/6
		case 2: *r = lx; *g = v ; *b = lz; break; // 2/6<hue<3/6
		case 3: *r = lx; *g = ly; *b = v; break; // 3/6<hue/4/6
		case 4: *r = lz; *g = lx; *b = v; break; // 4/6<hue<5/6
		case 5: *r = v; *g = lx; *b = ly; break; // 5/6<hue<1
	}
	
	//scale output values in [0;1]
	clamp_value_to_01(r);
	clamp_value_to_01(g);
	clamp_value_to_01(b);
}

/*
 * Function used to scale a range x0..x1 to a new range y0..y1: y = y0 + (y1 - y0) * (x - x0) / (x1 - x0)
 */
void compute_attributes_min_max(fftw_real *rho, fftw_real *vx, fftw_real *vy, fftw_real *fx, fftw_real *fy, float *min, float *max)
{
	int i;
	fftw_real magn, attr_min, attr_max;
	
	if(show_matter_attribute == DENSITY)
	{
		attr_min = attr_max = rho[0];
		for (i = 1; i < DIM * DIM; i++)
		{
			if(rho[i] < attr_min)
			{
				attr_min = rho[i];
			}
			if(rho[i] > attr_max)
			{
				attr_max = rho[i];
			}
		}
			
		*min = attr_min;
		*max = attr_max; 
	}
	else if(show_matter_attribute == VELOCITY_MAGNITUDE)
	{
		attr_min = attr_max = vector_magnitude(vx[0], vy[0]);
		for (i = 1; i < DIM * DIM; i++)
		{
			magn = vector_magnitude(vx[i], vy[i]);
			if(magn < attr_min)
			{
				attr_min = magn;
			}
			if(magn > attr_max)
			{
				attr_max = magn;
			}
		}
			
		*min = attr_min;
		*max = attr_max;
	}
	else if(show_matter_attribute == FORCE_MAGNITUDE)
	{
		attr_min = attr_max = vector_magnitude(fx[0], fy[0]);
		for (i = 1; i < DIM * DIM; i++)
		{
			magn = vector_magnitude(fx[i], fy[i]);
			if(magn < attr_min)
			{
				attr_min = magn;
			}
			if(magn > attr_max)
			{
				attr_max = magn;
			}
		}
			
		*min = attr_min;
		*max = attr_max;
	}
} 
 
float set_scale(float value, float x0, float x1, float y0, float y1)
{
	if(x1 >= 1)
		return ((y1 - y0) * (value - y0) / (x1 - x0)) + y0;
	else
		return ((x1 - x0) * (value - y0) / (y1 - y0)) + y0;
}
/*
 * Function used to clamp a vlue to a user predefined interval.
 * It can be also used as helper function to clamp glyph's positions to a cell's dimensions
 */
float clamp_value( float value, float min, float max)
{
	if( value < min)
	{
		value = min;
	}
	if( value > max)
	{
		value = max;
	}
	
	return value;
}

/* 
 * Draw fluid/smoke/matter as a triangle mesh 
 */
void draw_matter_fluid_smoke(fftw_real wn, fftw_real hn, fftw_real *rho, fftw_real *vx, fftw_real *vy, fftw_real *fx, fftw_real *fy)
{
	int i, j, idx0, idx1, idx2, idx3;
	double px0, py0, px1, py1, px2, py2, px3, py3;
		
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glBegin(GL_TRIANGLES);
		for (j = 0; j < DIM; j++)            //draw smoke
		{
			for (i = 0; i < DIM; i++)
			{
				px0 = wn + (fftw_real)i * wn;
				py0 = hn + (fftw_real)j * hn;
				idx0 = (j * DIM) + i;


				px1 = wn + (fftw_real)i * wn;
				py1 = hn + (fftw_real)(j + 1) * hn;
				idx1 = ((j + 1) * DIM) + i;


				px2 = wn + (fftw_real)(i + 1) * wn;
				py2 = hn + (fftw_real)(j + 1) * hn;
				idx2 = ((j + 1) * DIM) + (i + 1);


				px3 = wn + (fftw_real)(i + 1) * wn;
				py3 = hn + (fftw_real)j * hn;
				idx3 = (j * DIM) + (i + 1);

				if(show_matter_attribute == DENSITY)
				{
					draw_triangle_vertex( px0, py0, rho[idx0]);
					draw_triangle_vertex( px1, py1, rho[idx1]);
					draw_triangle_vertex( px2, py2, rho[idx2]);

					draw_triangle_vertex( px0, py0, rho[idx0]);
					draw_triangle_vertex( px2, py2, rho[idx2]);
					draw_triangle_vertex( px3, py3, rho[idx3]);
				}
				
				if(show_matter_attribute == VELOCITY_MAGNITUDE)
				{
					draw_triangle_vertex( px0, py0, vector_magnitude(vx[idx0], vy[idx0])*V_SCALE_FACTOR);
					draw_triangle_vertex( px1, py1, vector_magnitude(vx[idx1], vy[idx1])*V_SCALE_FACTOR);
					draw_triangle_vertex( px2, py2, vector_magnitude(vx[idx2], vy[idx2])*V_SCALE_FACTOR);

					draw_triangle_vertex( px0, py0, vector_magnitude(vx[idx0], vy[idx0])*V_SCALE_FACTOR);
					draw_triangle_vertex( px2, py2, vector_magnitude(vx[idx2], vy[idx2])*V_SCALE_FACTOR);
					draw_triangle_vertex( px3, py3, vector_magnitude(vx[idx3], vy[idx3])*V_SCALE_FACTOR);
				}
				
				if(show_matter_attribute == FORCE_MAGNITUDE)
				{
					draw_triangle_vertex( px0, py0, vector_magnitude(fx[idx0], fy[idx0])*F_SCALE_FACTOR);
					draw_triangle_vertex( px1, py1, vector_magnitude(fx[idx1], fy[idx1])*F_SCALE_FACTOR);
					draw_triangle_vertex( px2, py2, vector_magnitude(fx[idx2], fy[idx2])*F_SCALE_FACTOR);

					draw_triangle_vertex( px0, py0, vector_magnitude(fx[idx0], fy[idx0])*F_SCALE_FACTOR);
					draw_triangle_vertex( px2, py2, vector_magnitude(fx[idx2], fy[idx2])*F_SCALE_FACTOR);
					draw_triangle_vertex( px3, py3, vector_magnitude(fx[idx3], fy[idx3])*F_SCALE_FACTOR);
				}
			}
		}
	glEnd();
}

/*
 * Functions used for drawing the color legends.
 */
 
//Draw color legend for smoke/fluid/matter/scalar
void draw_matter_color_legend()
{
	unsigned short int i, val;
	const unsigned short int max_displayable_values = 32;
	float color_value, x_size, y_size;
	char buffer[10]={'\0'}, buffer2[30]={'\0'};
	
	//rectangle's length on X and Y axes
	x_size = (float) MAX_COLOR_BANDS*5;
	x_size /= matter_color_bands + 1;
	y_size = 20;
	
	//show attribute's name		
	if(show_matter_attribute == DENSITY)
	{
		sprintf(buffer2, "%s", "Scalar: density");
		glRasterPos2f(0, 40);
		glColor3f(1.0f, 1.0f, 1.0f);
		glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10,(const unsigned char*) buffer2);
	}
	else if(show_matter_attribute == VELOCITY_MAGNITUDE)
	{
		sprintf(buffer2, "%s", "Scalar: velocity magnitude");
		glRasterPos2f(0, 40);
		glColor3f(1.0f, 1.0f, 1.0f);
		glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10, (const unsigned char*)buffer2);
	}
	else if(show_matter_attribute == FORCE_MAGNITUDE)
	{
		sprintf(buffer2, "%s", "Scalar: force magnitude");
		glRasterPos2f(0, 40);
		glColor3f(1.0f, 1.0f, 1.0f);
		glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10, (const unsigned char*)buffer2);
	}
	
	for(i = 0; i <= matter_color_bands; i++)
	{
		//obtain the color value depending on the selected number of color bands
		color_value = (float) i;
		if(!matter_scale)
		{
			color_value /= matter_color_bands;
		}
		else
		{
			if(show_matter_attribute == VELOCITY_MAGNITUDE)
			{
				color_value *= 3.5;
			}
		}
	
		//draw rectangle containing the color set for the given value
		glBegin(GL_QUADS);
		set_colormap(color_value); 			//obtain the color for each value, corresponding to a color band
		glVertex2f(i*x_size, y_size);		//upper left corner
		glVertex2f((i+1)*x_size, y_size); 	//upper right corner
		glVertex2f((i+1)*x_size, 0); 		//down right corner
		glVertex2f(i*x_size, 0); 			//down left corner
		glEnd();
		
		//create a buffer to store the float value to be ouputed on the legend
		if(matter_scale)
		{
			if(show_matter_attribute == VELOCITY_MAGNITUDE)
			{
				color_value /= 3.5;
			}
			color_value = set_scale(color_value, matter_scale_lmin, matter_scale_lmax, 0, 1);
		}
	
		//apply clamping. clamp_lmin and clamp_lmax are set by the user
		if(matter_clamp)
		{
			color_value = clamp_value(color_value, matter_clamp_lmin, matter_clamp_lmax);
		}
		sprintf(buffer, "%0.3f", color_value);
		
		//display text
		//show color values
		//If there are more than 32 color bands and displayed values, the latter overlap and form a smaller color legend
		if(matter_color_bands <= max_displayable_values)
		{
			glRasterPos2f((i+0.5)*x_size, 30);
			//glClear(GL_COLOR_BUFFER_BIT);
			glColor3f(1.0f, 1.0f, 1.0f);
			glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10, (const unsigned char*) buffer);
		}
		//Thus, for color bands higher than 32 there are displayed only 32 values
		else
		{
			val = i%(matter_color_bands/max_displayable_values);
			
			if(val == 0)
			{
				glRasterPos2f((i+0.5)*x_size, 30);
				//glClear(GL_COLOR_BUFFER_BIT);
				glColor3f(1.0f, 1.0f, 1.0f);
				glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10, (const unsigned char*) buffer);
			}
		}
		//end of text displaying
	}
}

/* 
 * Draw color legend for glyphs, when direction coloring is considered
 */
void draw_glyph_color_legend(unsigned int winWidth, unsigned int winHeight)
{
	unsigned short int i;
	const unsigned short int max_displayable_values = 32;
	float color_value, x_size, y_size;
	char buffer[10]={'\0'}, buffer2[20]={'\0'};
	
	//rectangle's length on X and Y axes
	x_size = (float) MAX_COLOR_BANDS*5;
	x_size /= glyph_color_bands + 1;
	y_size = 20;
	
	if(show_glyph_attribute == VELOCITY)
	{
		sprintf(buffer2, "%s", "Vector field: velocity");
		glRasterPos2f(0, winHeight-40);
		glColor3f(1.0f, 1.0f, 1.0f);
		glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10, (const unsigned char*) buffer2);
	}
	else if(show_glyph_attribute == FORCE)
	{
		sprintf(buffer2, "%s", "Vector field: force");
		glRasterPos2f(0, winHeight-40);
		glColor3f(1.0f, 1.0f, 1.0f);
		glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10, (const unsigned char*) buffer2);
	}
	
	for(i = 0; i <= glyph_color_bands; i++)
	{
		//obtain the color value depending on the selected number of color bands
		color_value = (float) i;
		
		//create a buffer to store the float value to be ouputed on the legend
		sprintf(buffer, "%0.3f", color_value);
				
		//draw rectangle containing the color set for the given value
		glBegin(GL_QUADS);		
		direction_to_color(0, 0, color_value, LEGEND_TYPE);		//obtain the color for each value, corresponding to a color band
		glVertex2f(i*x_size, winHeight-y_size);					//down left corner
		glVertex2f((i+1)*x_size, winHeight-y_size); 			//down right corner
		glVertex2f((i+1)*x_size, winHeight); 					//up right corner
		glVertex2f(i*x_size, winHeight); 						//up left corner
		glEnd();
		
		//If there are more than 32 color bands and displayed values, the latter overlap and form a smaller color legend
		if(glyph_color_bands <= max_displayable_values)
		{
			glRasterPos2f(i*x_size, winHeight-30);
			glColor3f(0.8f, 0.8f, 0.8f);
			glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10, (const unsigned char*) buffer);
		}
		else
		{
			unsigned short int val = i%(glyph_color_bands/max_displayable_values);
			
			if(val == 0)
			{
				glRasterPos2f(i*x_size, winHeight-30);
				glColor3f(1.0f, 1.0f, 1.0f);
				glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10, (const unsigned char*) buffer);
			}
		}
	}
} 

/*
 * Functions used at step 3 for drawing glyphs
 * attribute can be either velocity or force; 
 * Inputs: x, y: hedgehogs indeces on the grid
 * 		   attribute_x, attribute_y: attribute's values on x and y axes
 * 		   cell_width, cell_height" glyph's cell's dimensions (since the glyphs can be resized, they have their own cells, appart from the fluid)
 */
//draw hedgehogs: vertical lines
void draw_hedgehog(grid_cell cell, vector data, float vec_scale)
{
	//top point's coordinates
	float top_x = cell.x + cell.width*1/2; 
	float top_y = cell.y;
	
	//base point coordinates
	float base_x = top_x;
	float base_y;
	
	//rotation angle 
	float angle = direction2angle(data.x, data.y);
	
	//compute vector magnitude and scale it
	float magnitude = vector_magnitude(data.x, data.y);
	magnitude *= vec_scale/20;
	
	//scale factors for the top point
	float scale_x;
	float scale_y;
	
	//maximum scale size is how many times  both x and y dimensions fit in the min(cell_height, cell_width)
	//minimum size is half of the min(cell_height, cell_width)			 
	if(cell.width >= cell.height)
	{
		base_y = cell.y - cell.height*3/4;
		
		//set hedgehog's height accordingly to the glyphs cell's dimensions
		scale_x = clamp_value(magnitude, cell.width/(2*(top_y-base_y)),  cell.height/(top_y-base_y));
		scale_y = clamp_value(magnitude, cell.height/(2*(top_y-base_y)), cell.height/(top_y-base_y));
	}
	else
	{
		base_y = cell.y - cell.width*3/4;
		//set hedgehog's height accordingly to the glyphs cell's dimensions
		scale_x = clamp_value(magnitude, cell.width/(2*(top_y-base_y)), cell.width/(top_y-base_y));
		scale_y = clamp_value(magnitude, cell.height/(2*(top_y-base_y)), cell.width/(top_y-base_y));
	}
	
    glPushMatrix();
    
    glTranslatef(top_x, top_y, 0);
	glRotatef(angle, 0, 0, 1.0);
	glScalef(scale_x, 0.9*scale_y, 1.0);
	glTranslatef(-1*top_x, -1*top_y, 1.0);					
	glBegin(GL_LINES);
	glVertex2f(base_x, base_y);
	glVertex2f(top_x, top_y);
	glEnd();
	
	glPopMatrix();
}

//Draw 3D glyphs
void draw_cones(grid_cell cell, vector data, float vec_scale)					
{
	//cone's height and base's radius
	float base_radius = cell.width*1/3;
	float cone_height = cell.height;
	
	//compute the base's center
	float cx = cell.x+ cell.width*1/2;
	float cy =  cell.y;
	
	//rotation angle
	float angle = direction2angle(data.x, data.y);
	
	//compute vector magnitude and scale it
	float magnitude = vector_magnitude(data.x, data.y);
	magnitude *= vec_scale/100;
	
	//scaling factors
	float base_scaling_factor;
	float height_scaling_factor;
	
	//clamp values so the cones won't overlap (too much). They do not exceed the maximum size = cell_height
	if(cell.width >= cell.height)
	{
		base_scaling_factor = clamp_value(magnitude, cell.width/(8*base_radius), cell.height/(2*base_radius));
		height_scaling_factor = clamp_value(magnitude, cell.height/(10*cone_height), cell.height/(2*cone_height));
	}
	//cones do not exceed the maximum size = cell_width
	else
	{
		base_scaling_factor = clamp_value(magnitude, cell.width/(8*base_radius), cell.width/(2*base_radius) );
		height_scaling_factor = clamp_value(magnitude, cell.height/(10*cone_height), cell.width/(2*cone_height) );
	} 				
	
    //push the current matrix on stack and perform transformations
    glPushMatrix();
    
    glTranslatef(cx, cy, 0.0);
    glRotatef(angle, 0, 0, 1);
    glRotatef(180, 0, 1, 0);
    glRotatef(90, 1, 0, 0);										
    glScalef(base_scaling_factor, height_scaling_factor, 1.0);					
	glutSolidCone(base_radius, cone_height, 30, 30);
	
	//delete matrix from stack			
	glPopMatrix();
}

void draw_ellipsoids(grid_cell cell, vector data, float vec_scale)
{
	//compute the ellipsoid's center
	float cx = cell.x + cell.width*1/2;
	float cy = cell.y + cell.height*1/2;
	
	//compute vector magnitude and scale it
	float magnitude = vector_magnitude(data.x, data.y);
	magnitude *= vec_scale/100;
	
	//rotation angle
	float angle = direction2angle(data.x, data.y);
	
	//compute scaling factors
	float scale_x;
	float scale_y;
	
	//radius of the globe
	float radius;
	
	if(cell.width >= cell.height)
	{
		radius = cell.height*1/4;
		
		scale_x = clamp_value(magnitude, cell.width/(4*radius), cell.height/(2*radius) );
		scale_y = clamp_value(magnitude, cell.height/(6*radius), cell.height/(1.75*radius) );
	}
	else
	{
		radius = cell.width*1/4;
		
		scale_x = clamp_value(magnitude, cell.width/(4*radius), cell.width/(1.75*radius) );
		scale_y = clamp_value(magnitude, cell.height/(6*radius), cell.width/(2*radius) );
	}

	//push the current matrix on stack and perform transformations
	glPushMatrix();
	
	glTranslatef(cx, cy, 0);
	glRotatef(angle, 0, 0, 1.0);
	glScalef(scale_x, scale_y, 1.0);
	glutSolidSphere(radius, 30, 30);
	
	//delete matrix from stack
	glPopMatrix();
}

//Draw 2D glyphs: Not updated, used initially for testing purpose
void draw_triangle(int x, int y, float attribute_x, float attribute_y, float cell_width, float cell_height, float vec_scale)
{
	//compute the vertices' coordinates
	float right_x = (fftw_real)x * cell_width + cell_width*1/4;
	float right_y = (fftw_real)y * cell_height + cell_height*1/4;
	float left_x = (fftw_real)x * cell_width + cell_width*3/4;
	float left_y = right_y;
	float top_x = (fftw_real)x * cell_width + cell_width*1/2; 
	float top_y;  
	
	//triangle's center coordinates
	float cx = top_x;
	float cy = (fftw_real)y * cell_height + cell_height*5/12;
	
	//compute scaling factors
	float scale_x = vec_scale * attribute_x;
	float scale_y = vec_scale * attribute_y;
	
	if(cell_width >= cell_height)
	{
		top_y = (fftw_real)y * cell_height + cell_height*3/4; 
	}
	else
	{
		top_y = (fftw_real)y * cell_height + cell_width*3/4; 
	}
	
	scale_x = clamp_value(scale_x, 0.75, 2);
	scale_y = clamp_value(scale_y, 0.75, 2);
	
	float angle = direction2angle(attribute_x, attribute_y);
	
	glClear(GL_DEPTH_BUFFER_BIT); 
    glMatrixMode(GL_MODELVIEW);																	//it is applied to object coordinates
    glPushMatrix();       																		//clones the previous matrix and puts it on the stack for applying transformations
    glTranslatef(cx, cy, 0.0);																	//translate to world coordinates
	glRotatef(angle, 0.0, 0.0, 1.0);															//rotate around Oz axis, in world coordinates (0,0,0)
	glScalef(scale_x, scale_y, 1.0);
	glTranslatef(-1*cx, -1*cy, 0.0);   					 										//translate to back to object coordinates
	//draw glyph
	glBegin(GL_TRIANGLES);
	glVertex2f( right_x, right_y);
	glVertex2f( left_x, left_y);
	glVertex2f( top_x, top_y);
	glEnd();
	
	glPopMatrix();																				//deletes matrix from stack, after drawing is done
		
}

void draw_filled_ellipse(float x, float y, float attribute_x, float attribute_y, float cell_width, float cell_height, float vec_scale)
{
	unsigned short int i;
	int num_segments = 20;
	float angle, xx, yy;
	
	//compute the ellipse's center
	float cx = (fftw_real)x * cell_width;
	float cy =  (fftw_real)y * cell_height;
	//compute ellipse's radius
	float rx = cell_width*1/4;
	float ry = cell_height*1/4;
	//compute scaling factors
	float scale_x = vec_scale * attribute_x;
	float scale_y = vec_scale * attribute_y;
	//clamp values so the ellipses won't overlap. They do not exceed the maximum size = cell_height
	scale_x = clamp_value(scale_x, 0.5, 1);
	scale_y = clamp_value(scale_y, 0.5, 2); 
	
	float theta = direction2angle(attribute_x, attribute_y);
	
	glClear(GL_DEPTH_BUFFER_BIT); 
    glMatrixMode(GL_MODELVIEW);																	//it is applied to object coordinates
    glPushMatrix();       																		//clones the previous matrix and puts it on the stack for applying transformations
    glTranslatef(cx, cy, 0.0);																	//translate to world coordinates
	glRotatef(theta, 0.0, 0.0, 1.0);																//rotate around Oz axis, in world coordinates (0,0,0)
	glScalef(scale_x, scale_y, 1.0);																
	glTranslatef(-1*cx, -1*cy, 0.0);   					 										//translate to back to object coordinates
	//draw the circle
    glBegin(GL_TRIANGLE_FAN);
    //set the center of circle 
    glVertex2f(cx, cy); 														
    for (i = 0; i <= num_segments; i++)   
    {
		angle = 2.0f * M_PI * (float)i/(float)num_segments;
		xx = rx * cosf(angle);																	//calculate the x component 
        yy = ry * sinf(angle);																	//calculate the y component  
        glVertex2f(cx + xx, cy + yy);
    }
    glEnd();
    
    glPopMatrix();																				//deletes matrix from stack, after drawing is done
}

void draw_glyphs(grid_cell cell, vector data, float vec_scale)
{	
	if(show_glyph_type == HEDGEHOGS)
	{
		glMatrixMode(GL_MODELVIEW);
		draw_hedgehog(cell, data, vec_scale);
						
	}
	else if(show_glyph_type == CONES)
	{
		glMatrixMode(GL_PROJECTION);						
		glLoadIdentity ();
		glOrtho(0, winWidth, 0, winHeight, -100, 100);
		glDisable(GL_DEPTH_TEST);
										
		glClear(GL_DEPTH_BUFFER_BIT);
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_COLOR_MATERIAL);
		
		glMatrixMode(GL_MODELVIEW);
		draw_cones(cell, data, vec_scale);
						
		//disable lightning
		glDisable(GL_LIGHTING);
		}
		else if(show_glyph_type == ELLIPSES)
		{
			glClear(GL_DEPTH_BUFFER_BIT);
			glEnable(GL_LIGHTING);
			glEnable(GL_LIGHT0);
			glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
			glEnable(GL_DEPTH_TEST);
			glEnable(GL_COLOR_MATERIAL);
		
			glMatrixMode(GL_MODELVIEW);
			draw_ellipsoids(cell, data, vec_scale);
						
			//disable lightning
			glDisable(GL_LIGHTING);
		}
}

//set the glyphs' color based on selected scalar field
void glyphs_scalar_color(float density, vector velocity, vector force)
{
	float magnitude = 0.0;
	
	if(show_matter_attribute == DENSITY)
	{
		set_colormap(density);
	}
	else if(show_matter_attribute == VELOCITY_MAGNITUDE)
	{
		magnitude = vector_magnitude( velocity.x, velocity.y);
		magnitude *= V_SCALE_FACTOR; 
						
		set_colormap(magnitude);
	}
	else if(show_matter_attribute == FORCE_MAGNITUDE)
	{
		magnitude = vector_magnitude( force.x, force.y);
		magnitude *= F_SCALE_FACTOR;
						
		set_colormap(magnitude);
	}				
}

/*
 * Interpolation
 */
 
//Linear interpolation used for the initial version of the bilinear interpolation 
fftw_real linear_interpolation(int x1, fftw_real f_x1, int x2, fftw_real f_x2, float x)
{
	fftw_real result = (fftw_real) ((x - x1)/(x2-x1)*f_x1  + (x2-x)/(x2-x1)*f_x2);
	return result;
}

/*
 * Inputs: i, j = coordinates of the grid cell in the original/initial/default grid
 * 		   x_offset, y_offset = how much a glyph moved from the original/initial position
 * 		   *attribute = array of data that needs to be interpolated
 */ 

//initial version
fftw_real bilinear_interpolation_initial(int i, int j, float x_offset, float y_offset, fftw_real *attribute)
{
	//compute cells' coordinates
	int idx1 = (j*DIM) + i;
	int idx2 = (j*DIM) + (i+1);
	int idx3 = ((j + 1)*DIM) + i;
	int idx4 = ((j + 1)*DIM) + (i+1);
	
	//compute the values of the attribute in the selected 4 points
	fftw_real f11 = attribute[idx1], f12 = attribute[idx2], f21 = attribute[idx3], f22 = attribute[idx4];
	
	//compute the point's coordinates in the given cell
	float s = i + x_offset, t  = j + y_offset;
	
	//interpolate on Ox axis
	fftw_real r1 = linear_interpolation(i, f11, i+1, f12, s);
	fftw_real r2 = linear_interpolation(i, f21, i+1, f22, s);
								
	//interpolate the results for the previous interpolations on Oy
	return linear_interpolation(j, r1, j+1, r2, t);
}

//functions used for indexing the attributes at the four corners of a cell
void cell_corners_indexing(int i, int j, cell_sample_points *cell_corners_indeces)
{
	(*cell_corners_indeces).idx1 = (j*DIM) + i;
	(*cell_corners_indeces).idx2 = (j*DIM) + (i+1);
	(*cell_corners_indeces).idx3 = ((j + 1)*DIM) + i;
	(*cell_corners_indeces).idx4 = ((j + 1)*DIM) + (i+1);
}

//currently used version based on the model provided in book, chapter 3.2
fftw_real bilinear_interpolation(cell_sample_points cell_corners_indeces, point data_point, fftw_real *attribute)
{
	//compute the values of the attribute in the selected 4 points
	fftw_real f11 = attribute[cell_corners_indeces.idx1], f12 = attribute[cell_corners_indeces.idx2], f21 = attribute[cell_corners_indeces.idx3], f22 = attribute[cell_corners_indeces.idx4];
	
	//compute the basis function values (linear basis functions)
	fftw_real q11 = (1 - data_point.r) * (1 - data_point.s); 
	fftw_real q12 = data_point.r * (1 - data_point.s); 
	fftw_real q21 = data_point.r * data_point.s;
	fftw_real q22 = (1 - data_point.r) * data_point.s; 
	
	//return the interpolation result: "weighted sum of the basis functions, where the weights are the sample values"
	return f11*q11 + f12*q12 + f21*q21 + f22*q22;
}								

/*
 * Functions used at step 4 for computing gradients 
 */

void compute_density_gradient(int i, int j, vector step, point data_point, fftw_real *density, vector *gradient)
{
	//Part1: compute gradient
	//compute cells' coordinates
	int idx1 = (j*DIM) + i;
	int idx2 = (j*DIM) + (i+1);
	int idx3 = ((j + 1)*DIM) + i;
	int idx4 = ((j + 1)*DIM) + (i+1);
	
	//compute the values in the sample points
	fftw_real f1 = density[idx1], f2 = density[idx2], f3 = density[idx3], f4 = density[idx4];
	
	(*gradient).x = (1 - data_point.r)*(f2 - f1)/step.x + data_point.r*(f4 - f3)/step.x;
	(*gradient).y = (1 - data_point.s)*(f3 - f1)/step.y + data_point.s*(f4 - f2)/step.y;
}

void compute_velocity_magnitude_gradient(int i, int j, vector step, point data_point, fftw_real *vx, fftw_real *vy, vector *gradient)
{
	//compute cells' coordinates
	int idx1 = (j*DIM) + i;
	int idx2 = (j*DIM) + (i+1);
	int idx3 = ((j + 1)*DIM) + i;
	int idx4 = ((j + 1)*DIM) + (i+1);
					
	//compute the values in the sample points
	fftw_real f1 = vector_magnitude(vx[idx1], vy[idx1]), f2 = vector_magnitude(vx[idx2], vy[idx2]);
	fftw_real f3 = vector_magnitude(vx[idx3], vy[idx3]), f4 = vector_magnitude(vx[idx4], vy[idx4]);
					
	//gradient definition from the book
	(*gradient).x = (1 - data_point.r)*(f2 - f1)/step.x + data_point.r*(f4 - f3)/step.x;
	(*gradient).y = (1 - data_point.s)*(f3 - f1)/step.y + data_point.s*(f4 - f2)/step.y;
}

/*
 * Functions used at step 5
 */
void draw_points(point current_point, vector velocity, float step)
{
	unsigned short int i;
	int num_segments = 20;
	float angle, xx, yy, radius = 2;
	fftw_real vel_magnitude;
	
	//compute the point's center
	float cx = current_point.r;
	//printf("a: %f\n", current_point.s);
	float cy = winHeight - current_point.s;
	//printf("b: %f\n", cy);
	
	//get the vector magnitude
	vel_magnitude = vector_magnitude(velocity.x, velocity.y);
	vel_magnitude *= V_SCALE_FACTOR;
	
	//set color based on velocity magnitude
	set_colormap(vel_magnitude);
	
	glClear(GL_DEPTH_BUFFER_BIT); 
    glMatrixMode(GL_MODELVIEW);																	//it is applied to object coordinates
    glPushMatrix();       																		//clones the previous matrix and puts it on the stack for applying transformations
    glBegin(GL_TRIANGLE_FAN);
    //set the center of circle 
    glVertex2f(cx, cy); 														
    for (i = 0; i <= num_segments; i++)   
    {
		angle = 2.0f * M_PI * (float)i/(float)num_segments;
		xx = radius * cosf(angle);																	//calculate the x component 
        yy = radius * sinf(angle);																	//calculate the y component  
        glVertex2f(cx + xx, cy + yy);
    }
    glEnd();
}

void draw_lines(point prev_point, point current_point, vector velocity)
{
	fftw_real vel_magnitude;
	
	//printf("%f %f %f %f\n", prev_point.r, winHeight - prev_point.s, current_point.s, winHeight - current_point.s);
	
	//get the vector magnitude
	vel_magnitude = vector_magnitude(velocity.x, velocity.y);
	vel_magnitude *= V_SCALE_FACTOR;
	
	//set color based on velocity magnitude
	set_colormap(vel_magnitude);
	
	glMatrixMode(GL_MODELVIEW);																	
    glPushMatrix();
	glBegin(GL_LINE_STRIP);
	//glBegin(GL_LINES);
	glVertex2f(prev_point.r, winHeight - prev_point.s);
	glVertex2f(current_point.r, winHeight - current_point.s);
	glEnd();
}
