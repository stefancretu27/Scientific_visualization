// Usage: Drag with the mouse to add smoke to the fluid. This will also move a "rotor" that disturbs
//        the velocity field at the mouse location. Press the indicated keys to change options
//--------------------------------------------------------------------------------------------------

#include <rfftw.h>              //the numerical simulation FFTW library
#include <stdio.h>              //for printing the help text
#include <math.h>               //for various math functions
#include <GL/glut.h>            //the GLUT graphics library

#define PI 3.1415927
#define MAX_COLOR_BANDS 256
#define max(x, y) x > y ? x : y
#define min(x, y) x < y ? x : y
#define vec_magnitude(x, y) sqrt(x*x + y*y)
#define bool short int
#define FALSE 0
#define TRUE 1

//--- SIMULATION PARAMETERS ------------------------------------------------------------------------
#define DIM 50				//size of simulation grid
double dt = 0.4;				//simulation time step
float visc = 0.001;				//fluid viscosity
fftw_real *vx, *vy;             //(vx,vy)   = velocity field at the current moment
fftw_real *vx0, *vy0;           //(vx0,vy0) = velocity field at the previous moment
fftw_real *fx, *fy;	            //(fx,fy)   = user-controlled simulation forces, steered with the mouse
fftw_real *rho, *rho0;			//smoke density at the current (rho) and previous (rho0) moment
rfftwnd_plan plan_rc, plan_cr;  //simulation domain discretization


//--- VISUALIZATION PARAMETERS ---------------------------------------------------------------------
int   winWidth, winHeight;      //size of the graphics window, in pixels
float vec_scale = 1000;			//scaling of hedgehogs
#define COLOR_BLACKWHITE 0   //different types of color mapping: black-and-white, rainbow, banded
#define COLOR_RAINBOW 1
#define COLOR_SAT 2
#define SHOW_DENSITY 0
#define SHOW_VELOCITY_MAGNITUDE 1
#define SHOW_FORCE_MAGNITUDE 2
#define V_SCALE_FACTOR 50
#define F_SCALE_FACTOR 30
#define SMOKE_TYPE 0
#define GLYPH_TYPE 1
#define LEGEND_TYPE 2  
bool draw_smoke = 0;           									//draw the smoke or not
bool draw_vecs = 1;            									//draw the vector field or not
unsigned short int scalar_color_type = COLOR_BLACKWHITE;        //method for smoke coloring
unsigned short int glyph_color_type = COLOR_BLACKWHITE;         //method for force&velocity coloring
bool frozen = FALSE;               								//toggles on/off the animation
unsigned int color_bands = 2;
unsigned short int show_attribute = SHOW_DENSITY;
bool enable_hsv = FALSE;
float scale_h = 0, scale_s = 0, scale_v = 0;
bool scalar_scale = FALSE;
float scale_lmin = 0, scale_lmax = 0;
bool scalar_clamp = FALSE;
float clamp_lmin = 0, clamp_lmax = 0;



//------ SIMULATION CODE STARTS HERE -----------------------------------------------------------------

//init_simulation: Initialize simulation data structures as a function of the grid size 'n'.
//                 Although the simulation takes place on a 2D grid, we allocate all data structures as 1D arrays,
//                 for compatibility with the FFTW numerical library.
void init_simulation(int n)
{
	int i; size_t dim;

	dim     = n * 2*(n/2+1)*sizeof(fftw_real);        //Allocate data structures
	vx       = (fftw_real*) malloc(dim);
	vy       = (fftw_real*) malloc(dim);
	vx0      = (fftw_real*) malloc(dim);
	vy0      = (fftw_real*) malloc(dim);
	dim     = n * n * sizeof(fftw_real);
	fx      = (fftw_real*) malloc(dim);
	fy      = (fftw_real*) malloc(dim);
	rho     = (fftw_real*) malloc(dim);
	rho0    = (fftw_real*) malloc(dim);
	plan_rc = rfftw2d_create_plan(n, n, FFTW_REAL_TO_COMPLEX, FFTW_IN_PLACE);
	plan_cr = rfftw2d_create_plan(n, n, FFTW_COMPLEX_TO_REAL, FFTW_IN_PLACE);

	for (i = 0; i < n * n; i++)                      //Initialize data structures to 0
	{ vx[i] = vy[i] = vx0[i] = vy0[i] = fx[i] = fy[i] = rho[i] = rho0[i] = 0.0f; }
}


//FFT: Execute the Fast Fourier Transform on the dataset 'vx'.
//     'dirfection' indicates if we do the direct (1) or inverse (-1) Fourier Transform
void FFT(int direction,void* vx)
{
	if(direction==1) rfftwnd_one_real_to_complex(plan_rc,(fftw_real*)vx,(fftw_complex*)vx);
	else             rfftwnd_one_complex_to_real(plan_cr,(fftw_complex*)vx,(fftw_real*)vx);
}

int clamp(float x)
{ return ((x)>=0.0?((int)(x)):(-((int)(1-(x))))); }

//solve: Solve (compute) one step of the fluid flow simulation
void solve(int n, fftw_real* vx, fftw_real* vy, fftw_real* vx0, fftw_real* vy0, fftw_real visc, fftw_real dt)
{
	fftw_real x, y, x0, y0, f, r, U[2], V[2], s, t;
	int i, j, i0, j0, i1, j1;

	for (i=0;i<n*n;i++)
	{
		vx[i] += dt*vx0[i];
		vx0[i] = vx[i];
		vy[i] += dt*vy0[i];
		vy0[i] = vy[i];
	}

	for ( x=0.5f/n,i=0 ; i<n ; i++,x+=1.0f/n )
		for ( y=0.5f/n,j=0 ; j<n ; j++,y+=1.0f/n )
		{
			x0 = n*(x-dt*vx0[i+n*j])-0.5f;
			y0 = n*(y-dt*vy0[i+n*j])-0.5f;
			i0 = clamp(x0); s = x0-i0;
			i0 = (n+(i0%n))%n;
			i1 = (i0+1)%n;
			j0 = clamp(y0); t = y0-j0;
			j0 = (n+(j0%n))%n;
			j1 = (j0+1)%n;
			vx[i+n*j] = (1-s)*((1-t)*vx0[i0+n*j0]+t*vx0[i0+n*j1])+s*((1-t)*vx0[i1+n*j0]+t*vx0[i1+n*j1]);
			vy[i+n*j] = (1-s)*((1-t)*vy0[i0+n*j0]+t*vy0[i0+n*j1])+s*((1-t)*vy0[i1+n*j0]+t*vy0[i1+n*j1]);
		}

	for(i=0; i<n; i++)
		for(j=0; j<n; j++)
		{
			vx0[i+(n+2)*j] = vx[i+n*j];
			vy0[i+(n+2)*j] = vy[i+n*j];
		}

	FFT(1,vx0);
	FFT(1,vy0);

	for (i=0;i<=n;i+=2)
	{
		x = 0.5f*i;
		for (j=0;j<n;j++)
		{
			y = j<=n/2 ? (fftw_real)j : (fftw_real)j-n;
			r = x*x+y*y;
			if ( r==0.0f ) continue;
			f = (fftw_real)exp(-r*dt*visc);
			U[0] = vx0[i  +(n+2)*j];
			V[0] = vy0[i  +(n+2)*j];
			U[1] = vx0[i+1+(n+2)*j];
			V[1] = vy0[i+1+(n+2)*j];

			vx0[i  +(n+2)*j] = f*((1-x*x/r)*U[0]     -x*y/r *V[0]);
			vx0[i+1+(n+2)*j] = f*((1-x*x/r)*U[1]     -x*y/r *V[1]);
			vy0[i+  (n+2)*j] = f*(  -y*x/r *U[0] + (1-y*y/r)*V[0]);
			vy0[i+1+(n+2)*j] = f*(  -y*x/r *U[1] + (1-y*y/r)*V[1]);
		}
	}

	FFT(-1,vx0);
	FFT(-1,vy0);

	f = 1.0/(n*n);
	for (i=0;i<n;i++)
		for (j=0;j<n;j++)
		{
			vx[i+n*j] = f*vx0[i+(n+2)*j];
			vy[i+n*j] = f*vy0[i+(n+2)*j];
		}
}


// diffuse_matter: This function diffuses matter that has been placed in the velocity field. It's almost identical to the
// velocity diffusion step in the function above. The input matter densities are in rho0 and the result is written into rho.
void diffuse_matter(int n, fftw_real *vx, fftw_real *vy, fftw_real *rho, fftw_real *rho0, fftw_real dt)
{
	fftw_real x, y, x0, y0, s, t;
	int i, j, i0, j0, i1, j1;

	for ( x=0.5f/n,i=0 ; i<n ; i++,x+=1.0f/n )
		for ( y=0.5f/n,j=0 ; j<n ; j++,y+=1.0f/n )
		{
			x0 = n*(x-dt*vx[i+n*j])-0.5f;
			y0 = n*(y-dt*vy[i+n*j])-0.5f;
			i0 = clamp(x0);
			s = x0-i0;
			i0 = (n+(i0%n))%n;
			i1 = (i0+1)%n;
			j0 = clamp(y0);
			t = y0-j0;
			j0 = (n+(j0%n))%n;
			j1 = (j0+1)%n;
			rho[i+n*j] = (1-s)*((1-t)*rho0[i0+n*j0]+t*rho0[i0+n*j1])+s*((1-t)*rho0[i1+n*j0]+t*rho0[i1+n*j1]);
		}
}

//set_forces: copy user-controlled forces to the force vectors that are sent to the solver.
//            Also dampen forces and matter density to get a stable simulation.
void set_forces(void)
{
	int i;
	for (i = 0; i < DIM * DIM; i++)
	{
		rho0[i]  = 0.995 * rho[i];
		fx[i] *= 0.85;
		fy[i] *= 0.85;
		vx0[i]    = fx[i];
		vy0[i]    = fy[i];
	}
}


//do_one_simulation_step: Do one complete cycle of the simulation:
//      - set_forces:
//      - solve:            read forces from the user
//      - diffuse_matter:   compute a new set of velocities
//      - gluPostRedisplay: draw a new visualization frame
void do_one_simulation_step(void)
{
	if (!frozen)
	{
		set_forces();
		solve(DIM, vx, vy, vx0, vy0, visc, dt);
		diffuse_matter(DIM, vx, vy, rho, rho0, dt);
		glutPostRedisplay();
	}
}


//------ VISUALIZATION CODE STARTS HERE -----------------------------------------------------------------
//Helper function to scale any value to [0;1]
void scale_value(float *value)
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

//rainbow: Implements a color palette, mapping the scalar 'value' to a rainbow color RGB. Taken form the book, section 5.2
void rainbow(float value, float* R, float* G, float* B)
{
	const float dx=0.8;
	scale_value(&value);
	
	value = (6-2*dx)*value+dx;
	*R = max(0.0,(3-fabs(value-4)-fabs(value-5))/2);
	*G = max(0.0,(4-fabs(value-2)-fabs(value-4))/2);
	*B = max(0.0,(3-fabs(value-1)-fabs(value-2))/2);
}

void grayscale(float value, float *R, float *G, float *B)
{
	*R = *G = *B = value;
}

//green_saturation: implements a green-based color palette. Taken form the book sample 5
void green_saturation(float value, float *R, float *G, float *B)
{
	float r = 0, g = 1, b = 0;
	
	scale_value(&value);
	
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

//red_saturation: implements a red-based color palette. Variation of the above color map
void red_saturation(float value, float *R, float *G, float *B)
{
	float r = 1, g = 0, b = 0;
	
	scale_value(&value);
	
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

//this is like dividing the interval [0;1] (which is considered in every coloring function)
//in a number of subintervals equal to "color_bands" and assign a color for each head of every interval
//color_bands = 2 => 3 colors, only RGB; color_bands = 4 => 5 colors, one between R and G (yellow) the other one between B and G (bright blue)
float set_color_bands(float value)
{
	value *= color_bands;
	value = (int)(value);
	return value/color_bands;
}

//convert from rgb to hsv. Taken from the book, section 3.6.3
void rgb2hsv(float r, float g , float b, float *h, float *s, float *v)
{
	//scale input values in [0;1]
	scale_value(&r);
	scale_value(&g);
	scale_value(&b);
	
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
	//scale_value(h);
	//scale_value(s);
	//scale_value(v);
}

//convert from hsv to rgb. Taken from the book, section 3.6.3
void hsv2rgb(float h, float s, float v, float *r, float *g, float *b)
{
	//scale input values in [0;1]
	scale_value(&h);
	scale_value(&s);
	scale_value(&v);
		
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
	scale_value(r);
	scale_value(g);
	scale_value(b);
}

//To scale a range x0..x1 to a new range y0..y1: y = y0 + (y1 - y0) * (x - x0) / (x1 - x0)
float set_scale(float value, float x0, float x1,  float y0, float y1)
{
     return ((y1 - y0) * (value - x0) / (x1 - x0)) + y0;
}

//set_colormap: Sets three different types of colormaps for smoke
void set_colormap(float scalar_value)
{
	float r, g, b, h, s, v;

	//apply color scalling. scale_lmin and scale_lmax are set by the user. The interval [0;1] is scaled to the user-inserted interval
	if(scalar_scale)
	{
		scalar_value = set_scale(scalar_value, scale_lmin, scale_lmax, 0, 1);
	}
	
	//apply clamping. clamp_lmin and clamp_lmax are set by the user
	if(scalar_clamp)
	{
		if(scalar_value < clamp_lmin)
		{
			scalar_value = clamp_lmin;
		}
		
		if(scalar_value > clamp_lmax)
		{
			scalar_value = clamp_lmax;
		}
	}

	//apply color bands
	scalar_value = set_color_bands(scalar_value);

	if (scalar_color_type == COLOR_BLACKWHITE)
	{
		grayscale(scalar_value, &r, &g, &b);
	}
	else if (scalar_color_type == COLOR_RAINBOW)
	{
		rainbow(scalar_value, &r, &g, &b);
	}
	else if (scalar_color_type == COLOR_SAT)
	{
		green_saturation(scalar_value, &r, &g, &b);
	}

	//if the user chose to insert special values for hue, saturation and brightness
	if(enable_hsv)
	{
		rgb2hsv(r, g, b, &h, &s, &v);
		// add scale values to current values
		h += scale_h;
		s += scale_s;
		v += scale_v;
		//convert to rgb
		hsv2rgb(h, s, v, &r, &g, &b);
		//set color
		glColor3f(r,g,b);
	}
	else
	{
		//set color
		glColor3f(r,g,b);
	}
}

//direction_to_color: Set the current color by mapping a direction vector (x,y) for the glyphs
void direction_to_color(float x, float y, float value, unsigned short int type)
{
	float r, g, b, h, s, v, color_value;
	
	if(type == GLYPH_TYPE)
	{
		//compute the value whose color is to be set
		color_value = atan2(y,x)/PI + 1;
	}
	if(type == LEGEND_TYPE)
	{
		color_value = value;
	}
		//apply color scalling. scale_lmin and scale_lmax are set by the user. The interval [0;1] is scaled to the user-inserted interval
	if(scalar_scale)
	{
		color_value = set_scale(color_value, scale_lmin, scale_lmax, 0, 1);
	}
	
	//apply clamping. clamp_lmin and clamp_lmax are set by the user
	if(scalar_clamp)
	{
		if(color_value < clamp_lmin)
		{
			color_value = clamp_lmin;
		}
		
		if(color_value > clamp_lmax)
		{
			color_value = clamp_lmax;
		}
	}

	//apply color bands
	color_value = set_color_bands(color_value);
	
	if (glyph_color_type == COLOR_BLACKWHITE)
	{
		grayscale(color_value, &r, &g, &b);
	}
	else if (glyph_color_type == COLOR_RAINBOW)
	{
		rainbow(color_value, &r, &g, &b);
	}
	else if (glyph_color_type == COLOR_SAT)
	{
		red_saturation(color_value, &r, &g, &b);
	}
	
	//if the user chose to insert special values for hue, saturation and brightness
	if(enable_hsv)
	{
		rgb2hsv(r, g, b, &h, &s, &v);
		// add scale values to current values
		h += scale_h;
		s += scale_s;
		v += scale_v;
		//convert to rgb
		hsv2rgb(h, s, v, &r, &g, &b);
	}

	//set the color
	glColor3f(r,g,b);
}

//draw color_legend
void draw_color_legend(unsigned short int type)
{
	unsigned short int i;
	const unsigned short int max_displayable_values = 32;
	float color_value, x_size, y_size;
	
	//rectangle's length on X and Y axes
	x_size = (float) MAX_COLOR_BANDS*5;
	x_size /= color_bands + 1;
	y_size = 20;
	
	
	for(i = 0; i <= color_bands; i++)
	{
		//obtain the color value depending on the selected number of color bands
		color_value = (float) i;
		color_value /= color_bands;
		
		//create a buffer to store the float value to be ouputed on the legend
		char buffer[10]={'\0'};
		sprintf(buffer, "%0.3f", color_value);
		
		if(type == SMOKE_TYPE)
		{
			//draw rectangle containing the color set for the given value
			glBegin(GL_QUADS);
			set_colormap(color_value); 			//obtain the color for each value, corresponding to a color band
			glVertex2f(i*x_size, y_size);		//upper left corner
			glVertex2f((i+1)*x_size, y_size); 	//upper right corner
			glVertex2f((i+1)*x_size, 0); 		//down right corner
			glVertex2f(i*x_size, 0); 			//down left corner
			glEnd();
		}
		if(type == GLYPH_TYPE)
		{
			//draw rectangle containing the color set for the given value
			glBegin(GL_QUADS);		
			direction_to_color(0, 0, color_value, LEGEND_TYPE);		//obtain the color for each value, corresponding to a color band
			glVertex2f(i*x_size, winHeight-y_size);					//down left corner
			glVertex2f((i+1)*x_size, winHeight-y_size); 			//down right corner
			glVertex2f((i+1)*x_size, winHeight); 					//up right corner
			glVertex2f(i*x_size, winHeight); 						//up left corner
			glEnd();
		}

		
		//If there are more than 32 color bands and displayed values, the latter overlap and form a smaller color legend
		if(color_bands <= max_displayable_values)
		{
			glRasterPos2f(i*x_size, 30);
			if(type == GLYPH_TYPE)
			{
				glRasterPos2f(i*x_size, winHeight-30);
			}
			glColor3f(0.8f, 0.8f, 0.8f);
			glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10, buffer);
		}
		//Thus, for color bands higher than 32 there are displayed only 32 values
		else
		{
			unsigned short int val = i%(color_bands/max_displayable_values);
			
			if(val == 0)
			{
				glRasterPos2f(i*x_size, 30);
				if(type == GLYPH_TYPE)
				{
					glRasterPos2f(i*x_size, winHeight-30);
				}
				glColor3f(0.8f, 0.8f, 0.8f);
				glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10, buffer);
			}
		}
	}
}

//helper function to set color and 2D position for a vertex 
//it removes a bit the code redundancy :D
void draw_triangle_vertex( double x, double y, float color_value)
{
	set_colormap(color_value);    
	glVertex2f(x, y);
}

//visualize: This is the main visualization function
void visualize(void)
{
	int        i, j, idx;
	fftw_real  wn = (fftw_real)winWidth / (fftw_real)(DIM + 1);   // Grid cell width
	fftw_real  hn = (fftw_real)winHeight / (fftw_real)(DIM + 1);  // Grid cell heigh

	if (draw_smoke)
	{	
		int idx0, idx1, idx2, idx3;
		double px0, py0, px1, py1, px2, py2, px3, py3;
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glBegin(GL_TRIANGLES);
		for (j = 0; j < DIM - 1; j++)            //draw smoke
		{
			for (i = 0; i < DIM - 1; i++)
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

				if(show_attribute == SHOW_DENSITY)
				{
					draw_triangle_vertex( px0, py0, rho[idx0]);
					draw_triangle_vertex( px1, py1, rho[idx1]);
					draw_triangle_vertex( px2, py2, rho[idx2]);

					draw_triangle_vertex( px0, py0, rho[idx0]);
					draw_triangle_vertex( px2, py2, rho[idx2]);
					draw_triangle_vertex( px3, py3, rho[idx3]);
				}
				
				if(show_attribute == SHOW_VELOCITY_MAGNITUDE)
				{
					draw_triangle_vertex( px0, py0, vec_magnitude(vx[idx0], vy[idx0])*V_SCALE_FACTOR);
					draw_triangle_vertex( px1, py1, vec_magnitude(vx[idx1], vy[idx1])*V_SCALE_FACTOR);
					draw_triangle_vertex( px2, py2, vec_magnitude(vx[idx2], vy[idx2])*V_SCALE_FACTOR);

					draw_triangle_vertex( px0, py0, vec_magnitude(vx[idx0], vy[idx0])*V_SCALE_FACTOR);
					draw_triangle_vertex( px2, py2, vec_magnitude(vx[idx2], vy[idx2])*V_SCALE_FACTOR);
					draw_triangle_vertex( px3, py3, vec_magnitude(vx[idx3], vy[idx3])*V_SCALE_FACTOR);
				}
				
				if(show_attribute == SHOW_FORCE_MAGNITUDE)
				{
					draw_triangle_vertex( px0, py0, vec_magnitude(fx[idx0], fy[idx0])*F_SCALE_FACTOR);
					draw_triangle_vertex( px1, py1, vec_magnitude(fx[idx1], fy[idx1])*F_SCALE_FACTOR);
					draw_triangle_vertex( px2, py2, vec_magnitude(fx[idx2], fy[idx2])*F_SCALE_FACTOR);

					draw_triangle_vertex( px0, py0, vec_magnitude(fx[idx0], fy[idx0])*F_SCALE_FACTOR);
					draw_triangle_vertex( px2, py2, vec_magnitude(fx[idx2], fy[idx2])*F_SCALE_FACTOR);
					draw_triangle_vertex( px3, py3, vec_magnitude(fx[idx3], fy[idx3])*F_SCALE_FACTOR);
				}
			}
		}
		glEnd();
		
		draw_color_legend(SMOKE_TYPE);
	}

	if (draw_vecs)		//draw velocities
	{
		glBegin(GL_LINES);				
		for (i = 0; i < DIM; i++)
			for (j = 0; j < DIM; j++)
			{
				idx = (j * DIM) + i;
				
				direction_to_color(vx[idx], vy[idx], 0, GLYPH_TYPE);
				
				glVertex2f(wn + (fftw_real)i * wn, hn + (fftw_real)j * hn);
				glVertex2f((wn + (fftw_real)i * wn) + vec_scale * vx[idx], (hn + (fftw_real)j * hn) + vec_scale * vy[idx]);
			}
		glEnd();
		
		draw_color_legend(GLYPH_TYPE);
	}
	/*
	 * else	//draw forces
	{
		glBegin(GL_LINES);				
		for (i = 0; i < DIM; i++)
			for (j = 0; j < DIM; j++)
			{
				idx = (j * DIM) + i;
				
				direction_to_color(fx[idx], fy[idx]);
				
				glVertex2f(wn + (fftw_real)i * wn, hn + (fftw_real)j * hn);
				glVertex2f((wn + (fftw_real)i * wn) + vec_scale * fx[idx], (hn + (fftw_real)j * hn) + vec_scale * fy[idx]);
			}
		glEnd();
	}
	* */
	
}


//------ INTERACTION CODE STARTS HERE -----------------------------------------------------------------

//display: Handle window redrawing events. Simply delegates to visualize().
void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	visualize();
	glFlush();
	glutSwapBuffers();
}

//reshape: Handle window resizing (reshaping) events
void reshape(int w, int h)
{
	glViewport(0.0f, 0.0f, (GLfloat)w, (GLfloat)h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, (GLdouble)w, 0.0, (GLdouble)h);
	winWidth = w; winHeight = h;
}

//keyboard: Handle key presses
void keyboard(unsigned char key, int x, int y)
{
	switch (key)
	{
		case 't': dt -= 0.001; break;
		case 'T': dt += 0.001; break; 
		case 'c': glyph_color_type++;
			if (glyph_color_type > COLOR_SAT)
			{
				glyph_color_type = COLOR_BLACKWHITE;
			}
			break;
		case 'S': vec_scale *= 1.2; break;
		case 's': vec_scale *= 0.8; break;
		case 'V': visc *= 5; break;
		case 'v': visc *= 0.2; break;
		case 'x':
			draw_smoke = 1 - draw_smoke;
			if (draw_smoke==0) draw_vecs = 1;
			break;
		case 'y':
			draw_vecs = 1 - draw_vecs;
			if (draw_vecs==0) draw_smoke = 1;
			break;
		case 'm':
			scalar_color_type++;
			if (scalar_color_type > COLOR_SAT)
			{
				scalar_color_type = COLOR_BLACKWHITE;
			}
			break;
		case 'b': color_bands *= 2;
			if(color_bands > MAX_COLOR_BANDS)
			{
				color_bands = 2;
			}
			break;
		case 'n':
			show_attribute++;
			if (show_attribute > SHOW_FORCE_MAGNITUDE)
			{
				show_attribute = SHOW_DENSITY;
			}
			break;
		case 'k': scalar_scale = 1-scalar_scale; 
			if(scalar_scale)
			{
				scalar_clamp = 0;
				printf("Set fluid color scale interval range\n");
				printf("Insert lower value and upper value:\n");
				scanf("%f %f", &scale_lmin, &scale_lmax);
			}
			break;
		case 'j': scalar_clamp = 1-scalar_clamp; 
			if(scalar_clamp)
			{
				scalar_scale = 0;
				printf("Set fluid color clamp interval range\n");
				printf("Insert lower value and upper value:\n");
				scanf("%f %f", &clamp_lmin, &clamp_lmax);
			}
			break;
		case 'l': enable_hsv = 1 - enable_hsv; 
			if(enable_hsv)
			{
				printf("Increase h,s,v values by:\n");
				scanf("%f %f %f", &scale_h, &scale_s, &scale_v);
			}
			break; 
		case 'a': frozen = 1-frozen; break;
		case 'q': exit(0);
	}
}

// drag: When the user drags with the mouse, add a force that corresponds to the direction of the mouse
//       cursor movement. Also inject some new matter into the field at the mouse location.
void drag(int mx, int my)
{
	int xi,yi,X,Y;
	double  dx, dy, len;
	static int lmx=0,lmy=0;				//remembers last mouse location

	// Compute the array index that corresponds to the cursor location
	xi = (int)clamp((double)(DIM + 1) * ((double)mx / (double)winWidth));
	yi = (int)clamp((double)(DIM + 1) * ((double)(winHeight - my) / (double)winHeight));

	X = xi;
	Y = yi;

	if (X > (DIM - 1))
		X = DIM - 1;
	if (Y > (DIM - 1))
		Y = DIM - 1;
	if (X < 0)
		X = 0;
	if (Y < 0)
		Y = 0;

	// Add force at the cursor location
	my = winHeight - my;
	dx = mx - lmx;
	dy = my - lmy;
	len = sqrt(dx * dx + dy * dy);
	if (len != 0.0)
	{
		dx *= 0.1 / len;
		dy *= 0.1 / len;
	}
	fx[Y * DIM + X] += dx;
	fy[Y * DIM + X] += dy;
	rho[Y * DIM + X] = 10.0f;
	lmx = mx;
	lmy = my;
}


//main: The main program
int main(int argc, char **argv)
{
	printf("Fluid Flow Simulation and Visualization\n");
	printf("=======================================\n");
	printf("Click and drag the mouse to steer the flow!\n");
	printf("T/t:   increase/decrease simulation timestep\n");
	printf("S/s:   increase/decrease hedgehog scaling\n");
	printf("c:     toggle direction coloring on/off\n");
	printf("V/v:   increase decrease fluid viscosity\n");
	printf("x:     toggle drawing matter on/off\n");
	printf("y:     toggle drawing hedgehogs on/off\n");
	printf("a:     toggle the animation on/off\n");
	printf("m:     toggle thru scalar coloring\n");
	printf("n:     toggle thru shown attributes (density, velocity, force) \n");
	printf("b:     change color bands for color mapping (2...256)\n");
	printf("j:     toggle the color map clamping on/off\n");
	printf("k:     toggle the color map scaling on/off\n");
	printf("l:     toggle the hsv coloring on/off\n");
	printf("q:     quit\n\n");

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize(1300,800);
	glutCreateWindow("Real-time smoke simulation and visualization");
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutIdleFunc(do_one_simulation_step);
	glutKeyboardFunc(keyboard);
	glutMotionFunc(drag);

	init_simulation(DIM);	//initialize the simulation data structures
	glutMainLoop();			//calls do_one_simulation_step, keyboard, display, drag, reshape
	return 0;
}
