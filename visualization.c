/*
 * Includes 
 */
#include <rfftw.h>              //the numerical simulation FFTW 
#include <stdio.h>              //for printing the help text

#include <math.h>               //for various math functions #include
#include <GL/glut.h>            //the GLUT graphics library

#include <fluids.h>
#include <visualization.h>

/*
 * External data: They are declared here but not instantiated. The memory allocation and initialization happens in fluids.c
 */
extern int winWidth, winHeight; 
extern float vec_scale;
 
//MATTER_GLOBALS          										//draw the matter (1) or not (0)extern 
extern matter_attribute show_matter_attribute;							//default shown attribute
extern colormap matter_color_type;        								//default method for smoke attributes coloring
extern unsigned int matter_color_bands;									//default number of color bands


//GLYPH_GLOBALS
extern colormap glyph_color_type;         								//method for force&velocity coloring
extern unsigned int glyph_color_bands;									//default number of color bands
extern glyphs_attribute show_glyph_attribute;


/*
 * Helper functions section
 */

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
	float n = vec_magnitude(x, y); 
	
	//normalize the vector
	if (n<1.0e-6) 
	{
		n=1;
	} 
		
	x/=n; 
	y/=n; 

	return atan2(x,-y) * (180 / M_PI);

	/*float cosa = x;
	float sina = y;
	float a;
	
	if (sina >= 0)
	{
		a = acos(cosa);
	}
	else
	{
		a = 2*M_PI - acos(cosa);
	}*/
		
	//return 180*a/M_PI;
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
	
	//scale output values in [0;1]
	//clamp_value_to_01(h);
	//clamp_value_to_01(s);
	//clamp_value_to_01(v);
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
float set_scale(float value, float x0, float x1,  float y0, float y1)
{
     return ((y1 - y0) * (value - x0) / (x1 - x0)) + y0;
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
 * Functions used for drawing the color legends.
 */
 
/* 
 * Draw color legend for smoke/fluid/matter/scalar
 */

void draw_matter_color_legend()
{
	unsigned short int i;
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
		glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10, buffer2);
	}
	else if(show_matter_attribute == VELOCITY_MAGNITUDE)
	{
		sprintf(buffer2, "%s", "Scalar: velocity magnitude");
		glRasterPos2f(0, 40);
		glColor3f(1.0f, 1.0f, 1.0f);
		glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10, buffer2);
	}
	else if(show_matter_attribute == FORCE_MAGNITUDE)
	{
		sprintf(buffer2, "%s", "Scalar: force magnitude");
		glRasterPos2f(0, 40);
		glColor3f(1.0f, 1.0f, 1.0f);
		glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10, buffer2);
	}
	
	for(i = 0; i <= matter_color_bands; i++)
	{
		//obtain the color value depending on the selected number of color bands
		color_value = (float) i;
		color_value /= matter_color_bands;
		
		//create a buffer to store the float value to be ouputed on the legend
		sprintf(buffer, "%0.3f", color_value);
	
		//draw rectangle containing the color set for the given value
		glBegin(GL_QUADS);
		set_colormap(color_value); 			//obtain the color for each value, corresponding to a color band
		glVertex2f(i*x_size, y_size);		//upper left corner
		glVertex2f((i+1)*x_size, y_size); 	//upper right corner
		glVertex2f((i+1)*x_size, 0); 		//down right corner
		glVertex2f(i*x_size, 0); 			//down left corner
		glEnd();
		
		//display text
		//show color values
		//If there are more than 32 color bands and displayed values, the latter overlap and form a smaller color legend
		if(matter_color_bands <= max_displayable_values)
		{
			glRasterPos2f(i*x_size, 30);
			glColor3f(1.0f, 1.0f, 1.0f);
			glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10, buffer);
		}
		//Thus, for color bands higher than 32 there are displayed only 32 values
		else
		{
			unsigned short int val = i%(matter_color_bands/max_displayable_values);
			
			if(val == 0)
			{
				glRasterPos2f(i*x_size, 30);
				glColor3f(1.0f, 1.0f, 1.0f);
				glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10, buffer);
			}
		}
		//end of text displaying
	}
}

/* 
 * Draw color legend for glyphs, when direction coloring is considered
 */
void draw_glyph_color_legend()
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
		glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10, buffer2);
	}
	else if(show_glyph_attribute == FORCE)
	{
		sprintf(buffer2, "%s", "Vector field: force");
		glRasterPos2f(0, winHeight-40);
		glColor3f(1.0f, 1.0f, 1.0f);
		glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10, buffer2);
	}
	
	for(i = 0; i <= glyph_color_bands; i++)
	{
		//obtain the color value depending on the selected number of color bands
		color_value = (float) i;
		color_value /= glyph_color_bands;
		
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
			glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10, buffer);
		}
		else
		{
			unsigned short int val = i%(glyph_color_bands/max_displayable_values);
			
			if(val == 0)
			{
				glRasterPos2f(i*x_size, winHeight-30);
				glColor3f(1.0f, 1.0f, 1.0f);
				glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10, buffer);
			}
		}
	}
} 

/*
 * Functions used for draeing glyphs
 * attribute can be either velocity or force; 
 * Inputs: x, y: hedgehogs indeces on the grid
 * 		   attribute_x, attribute_y: attribute's values on x and y axes
 * 		   cell_width, cell_height" glyph's cell's dimensions (since the glyphs can be resized, they have their own cells, appart from the fluid)
 */
//draw hedgehogs: vertical lines
void draw_hedgehog(int x, int y, float attribute_x, float attribute_y, float cell_width, float cell_height)
{
	//base point coordinates
	float base_x = (fftw_real)x * cell_width + cell_width*1/2;
	float base_y = (fftw_real)y * cell_height;
	//top point's coordinates
	float top_x = base_x; 
	float top_y;
	
	//scale factors for the top point
	float scale_x = vec_scale * attribute_x;
	float scale_y = vec_scale * attribute_y;
	
	
	//maximum scale size is how many times  both x and y dimensions fit in the min(cell_height, cell_width)
	//minimum size is half of the min(cell_height, cell_width)			 
	if(cell_width >= cell_height)
	{
		top_y = base_y + cell_height*1/2;

		//set hedgehog's height accordingly to the glyphs cell's dimensions
		scale_x = clamp_value(scale_x, cell_height/(4*cell_width*1/2), cell_height/(cell_width*1/2) );
		scale_y = clamp_value(scale_y, cell_height/(4*(top_y-base_y)), cell_height/(top_y-base_y) );
	}
	else
	{
		top_y = base_y + cell_width*1/2; 
		
		//set hedgehog's height accordingly to the glyphs cell's dimensions
		scale_x = clamp_value(scale_x, cell_width/(4*cell_width*1/2), cell_width/(cell_width*1/2) );
		scale_y = clamp_value(scale_y, cell_width/(4*(top_y-base_y)), cell_width/(top_y-base_y) );
	}

	float angle = direction2angle(attribute_x, attribute_y);
	
    glPushMatrix();
    
    glTranslatef(base_x, base_y, 0);
	glRotatef(angle, 0, 0, 1.0);
	glScalef(scale_x, scale_y, 1.0);
	glTranslatef(-1*base_x, -1*base_y, 1.0);					
	glBegin(GL_LINES);
	glVertex2f( base_x, base_y);
	glVertex2f(top_x, top_y);
	glEnd();
	
	glPopMatrix();
}
 
//Draw 2D glyphs
void draw_triangle(int x, int y, float attribute_x, float attribute_y, float cell_width, float cell_height)
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
	glRotatef(angle, 0.0, 0.0, 1.0);																//rotate around Oz axis, in world coordinates (0,0,0)
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

void draw_filled_ellipse(float x, float y, float attribute_x, float attribute_y, float cell_width, float cell_height)
{
	unsigned short int i;
	int num_segments = 20;
	float angle, xx, yy;
	
	//compute the ellipse's center
	float cx = (fftw_real)x * cell_width;
	float cy =  (fftw_real)y * cell_height;
	//compute ellipse's radius, including scaling
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

//Draw 3D glyphs
void draw_cones(int x, int y, float attribute_x, float attribute_y, float cell_width, float cell_height)					
{
	//cone's height and base's radius
	float base_radius = cell_width*1/2;
	float height = cell_height;
	
	//scaling factors
	float base_scaling_factor = vec_scale * attribute_x;
	float height_scaling_factor = vec_scale * attribute_y;
	
	//compute the base's center
	float cx = (fftw_real)x * cell_width + cell_width*1/2;
	float cy =  (fftw_real)y * cell_height;
	
	//clamp values so the ellipses won't overlap. They do not exceed the maximum size = cell_height
	if(cell_width >= cell_height)
	{
		base_scaling_factor = clamp_value(base_scaling_factor, cell_height/(4*base_radius), cell_height/(2*base_radius) );
		height_scaling_factor = clamp_value(height_scaling_factor, cell_height/(4*height), cell_height/(2*height) );
	}
	else
	{
		base_scaling_factor = clamp_value(base_scaling_factor, cell_width/(4*base_radius), cell_width/(2*base_radius) );
		height_scaling_factor = clamp_value(height_scaling_factor, cell_width/(4*height), cell_width/(2*height) );
	} 
				
	float angle = direction2angle(attribute_x, attribute_y);				
	
	//glMatrixMode(GL_MODELVIEW);
    //push the current matrix on stack and perform transformations
    glPushMatrix();
    glTranslatef(cx, cy, 0.0);
    glRotatef(angle, 0, 0, 1);
    glRotatef(180, 0, 1, 0);
    glRotatef(90, 1, 0, 0);										
    glScalef(base_scaling_factor, height_scaling_factor, 1.0);					
	glutSolidCone(base_radius, height, 30, 30);
	//delete matrix from stack			
	glPopMatrix();
}

void draw_ellipsoids(int x, int y, float attribute_x, float attribute_y, float cell_width, float cell_height)
{
	//compute the ellipsoid's center
	float cx = (fftw_real)x * cell_width + cell_width*1/2;
	float cy =  (fftw_real)y * cell_height + cell_height*1/2;
	//compute scaling factors
	float scale_x = vec_scale * attribute_x;
	float scale_y = vec_scale * attribute_y;
	//radius of the globe
	float radius;
	
	if(cell_width >= cell_height)
	{
		radius = cell_height*1/4;
		
		scale_x = clamp_value(scale_x, cell_width/(4*radius), cell_height/(2*radius) );
		scale_y = clamp_value(scale_y, cell_height/(4*radius), cell_height/(2*radius) );
	}
	else
	{
		radius = cell_width*1/4;
		
		scale_x = clamp_value(scale_x, cell_width/(4*radius), cell_width/(2*radius) );
		scale_y = clamp_value(scale_y, cell_height/(4*radius), cell_width/(2*radius) );
	}
	
	float theta = direction2angle(attribute_x, attribute_y);

	//glMatrixMode(GL_MODELVIEW);

	//push the current matrix on stack and perform transformations
	glPushMatrix();
	glTranslatef(cx, cy, 0);
	glRotatef(theta, 0, 0, 1.0);
	glScalef(scale_x, scale_y, 1.0);
	glutSolidSphere(radius, 30, 30);
	//delete matrix from stack
	glPopMatrix();
}
