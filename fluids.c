// Usage: Drag with the mouse to add smoke to the fluid. This will also move a "rotor" that disturbs
//        the velocity field at the mouse location. Press the indicated keys to change options
//--------------------------------------------------------------------------------------------------

#include <rfftw.h>              //the numerical simulation FFTW 
#include <stdio.h>              //for printing the help text

#include <math.h>               //for various math functions #include
#include <GL/glut.h>            //the GLUT graphics library

#include <fluids.h>
#include <visualization.h>

//--- SIMULATION PARAMETERS ------------------------------------------------------------------------
#define DIM 50
unsigned short int no_glyphs_x = 50;
unsigned short int no_glyphs_y = 50;
double dt = 0.4;													//simulation time step
float visc = 0.001;													//fluid viscosity
fftw_real *vx, *vy;             									//(vx,vy)   = velocity field at the current moment
fftw_real *vx0, *vy0;          										//(vx0,vy0) = velocity field at the previous moment
fftw_real *fx, *fy;	            									//(fx,fy) = user-controlled simulation forces, steered with the mouse
fftw_real *rho, *rho0;												//smoke density at the current (rho) and previous (rho0) moment
rfftwnd_plan plan_rc, plan_cr;  									//simulation domain discretization


//--- VISUALIZATION PARAMETERS ---------------------------------------------------------------------
int   winWidth, winHeight;      										//size of the graphics window, in pixels
float vec_scale = 1000;													//scaling of hedgehogs
bool frozen = FALSE;   													//toggles on/off the animation

//scalar-related global data: matter
bool draw_matter = FALSE;           									//draw the matter (1) or not (0)
matter_attribute show_matter_attribute = DENSITY;						//default shown attribute
colormap matter_color_type = BLACKWHITE;        						//default method for smoke attributes coloring
unsigned int matter_color_bands = 2;									//default number of color bands
bool matter_enable_hsv = FALSE;											//default hsv colormap
float matter_scale_h = 0, matter_scale_s = 0, matter_scale_v = 0;		//amount by which h,s,v are increased
bool matter_scale = FALSE;												//default colormap scaling
float matter_scale_lmin = 0, matter_scale_lmax = 0;						//the limits of the interval to be scaled
bool matter_clamp = FALSE;												//default colormap clamping
float matter_clamp_lmin = 0, matter_clamp_lmax = 0;						//the limits to which a value is clamped

//vector-related global data: glyphs
bool draw_vecs = TRUE;													//draw glyphs or not (if not, matter is drawn)
bool color_dir = TRUE;            										//use direction color-coding or not (if not, glyphs are colored using scalar values => 2 types of data are shown)
colormap glyph_color_type = BLACKWHITE;         						//method for force&velocity coloring
unsigned int glyph_color_bands = 2;										//default number of color bands
glyphs_attribute show_glyph_attribute = VELOCITY;  						//draw the velocity or force
glyphs_type show_glyph_type = HEDGEHOGS;
bool glyph_enable_hsv = FALSE;											//default hsv colormap
float glyph_scale_h = 0, glyph_scale_s = 0, glyph_scale_v = 0;			//amount by which h,s,v are increased
bool glyph_scale = FALSE;												//default colormap scaling
float glyph_scale_lmin = 0, glyph_scale_lmax = 0;						//the limits of the interval to be scaled
bool glyph_clamp = FALSE;												//default colormap clamping
float glyph_clamp_lmin = 0, glyph_clamp_lmax = 0;						//the limits to which a value is clamped


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
//set_colormap: Sets three different types of colormaps for matter
void set_colormap(float scalar_value)
{
	float r, g, b, h, s, v;

	//apply color scalling. scale_lmin and scale_lmax are set by the user. The interval [0;1] is scaled to the user-inserted interval
	if(matter_scale)
	{
		scalar_value = set_scale(scalar_value, matter_scale_lmin, matter_scale_lmax, 0, 1);
	}
	
	//apply clamping. clamp_lmin and clamp_lmax are set by the user
	if(matter_clamp)
	{
		scalar_value = clamp_value(scalar_value, matter_clamp_lmin, matter_clamp_lmax);
	}

	//apply color bands
	scalar_value = set_color_bands(scalar_value, SCALAR);
	
	//apply color mapping
	if (matter_color_type == BLACKWHITE)
	{
		glClearColor(0.6, 0.4, 0.2, 0);							//set background color
		grayscale(scalar_value, &r, &g, &b);
	}
	else if (matter_color_type == RAINBOW)
	{
		glClearColor(0.0, 0.0, 0.0, 0);							//set background color
		rainbow(scalar_value, &r, &g, &b);
	}
	else if (matter_color_type == SATURATION)
	{
		glClearColor(0.2, 0.0, 0.2, 0);							//set background color
		green_saturation(scalar_value, &r, &g, &b);
	}

	//if the user chose to insert special values for hue, saturation and brightness
	if(matter_enable_hsv)
	{
		rgb2hsv(r, g, b, &h, &s, &v);
		// add scale values to current values
		h += matter_scale_h;
		s += matter_scale_s;
		v += matter_scale_v;
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
	
	//if function is called for drawing glyphs
	if(type == GLYPH_TYPE)
	{
		//compute the value whose color is to be set
		color_value = atan2(y,x)/M_PI + 1;
	}
	//if function is called for drawing color legend for glyphs
	if(type == LEGEND_TYPE)
	{
		color_value = value;
	}
	
	//apply color scalling. scale_lmin and scale_lmax are set by the user. The interval [0;1] is scaled to the user-inserted interval
	if(glyph_scale)
	{
		color_value = set_scale(color_value, glyph_scale_lmin, glyph_scale_lmax, 0, 1);
	}
	
	//apply clamping. clamp_lmin and clamp_lmax are set by the user
	if(glyph_clamp)
	{
		color_value = clamp_value(color_value, glyph_clamp_lmin, glyph_clamp_lmax);
	}

	//apply color bands
	color_value = set_color_bands(color_value, VECTORIAL);
	
	if (glyph_color_type == BLACKWHITE)
	{
		glClearColor(0.6, 0.4, 0.2, 0);
		grayscale(color_value, &r, &g, &b);
	}
	else if (glyph_color_type == RAINBOW)
	{
		glClearColor(0.0, 0.0, 0.0, 0);
		rainbow(color_value, &r, &g, &b);
	}
	else if (glyph_color_type == SATURATION)
	{
		glClearColor(0.0, 0.2, 0.2, 0);
		red_saturation(color_value, &r, &g, &b);
	}
	
	//if the user chose to insert special values for hue, saturation and brightness
	if(glyph_enable_hsv)
	{
		rgb2hsv(r, g, b, &h, &s, &v);
		// add scale values to current values
		h += glyph_scale_h;
		s += glyph_scale_s;
		v += glyph_scale_v;
		//convert to rgb
		hsv2rgb(h, s, v, &r, &g, &b);
	}

	//set the color
	glColor3f(r,g,b);
}

//visualize: This is the main visualization function
void visualize(void)
{
	int        i, j, idx;
	fftw_real  wn = (fftw_real)winWidth / (fftw_real)(DIM + 1);   // Grid cell width
	fftw_real  hn = (fftw_real)winHeight / (fftw_real)(DIM + 1);  // Grid cell heigh

	if (draw_matter)
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
					draw_triangle_vertex( px0, py0, vec_magnitude(vx[idx0], vy[idx0])*V_SCALE_FACTOR);
					draw_triangle_vertex( px1, py1, vec_magnitude(vx[idx1], vy[idx1])*V_SCALE_FACTOR);
					draw_triangle_vertex( px2, py2, vec_magnitude(vx[idx2], vy[idx2])*V_SCALE_FACTOR);

					draw_triangle_vertex( px0, py0, vec_magnitude(vx[idx0], vy[idx0])*V_SCALE_FACTOR);
					draw_triangle_vertex( px2, py2, vec_magnitude(vx[idx2], vy[idx2])*V_SCALE_FACTOR);
					draw_triangle_vertex( px3, py3, vec_magnitude(vx[idx3], vy[idx3])*V_SCALE_FACTOR);
				}
				
				if(show_matter_attribute == FORCE_MAGNITUDE)
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
		
		//draw color legend for matter
		draw_matter_color_legend();
	}

	if (draw_vecs)		//draw vectors
	{
		fftw_real glyph_x_dim = (fftw_real)winWidth / (fftw_real)(no_glyphs_x);   // Grid cell width
		fftw_real glyph_y_dim = (fftw_real)winHeight / (fftw_real)(no_glyphs_y);  // Grid cell heigh
		float magnitude = 0.0;
						
		for (i = 0; i < no_glyphs_x; i++)
			for (j = 2; j < no_glyphs_y-1; j++)
			{
				idx = (j * DIM) + i;
				
				//choose vectorial attribute
				if(show_glyph_attribute == VELOCITY)
				{	
					//choose glyphs' color
					//glyphs are colored based on velocity direction
					if(color_dir == TRUE)
					{
						direction_to_color(vx[idx], vy[idx], 0, GLYPH_TYPE);
					}
					else
					{
						//set the glyphs' color based on selected scalar field
						if(show_matter_attribute == DENSITY)
						{
							set_colormap(rho[idx]);
						}
						else if(show_matter_attribute == VELOCITY_MAGNITUDE)
						{
							magnitude = vec_magnitude( vx[idx], vy[idx]);
							magnitude *= V_SCALE_FACTOR;
						
							set_colormap(magnitude);
						}
						else if(show_matter_attribute == FORCE_MAGNITUDE)
						{
							magnitude = vec_magnitude( fx[idx], fy[idx]);
							magnitude *= F_SCALE_FACTOR;
						
							set_colormap(magnitude);
						}
					}
				
					//choose glyphs' type
					if(show_glyph_type == HEDGEHOGS)
					{
						glMatrixMode(GL_MODELVIEW);
						draw_hedgehog(i, j, vx[idx], vy[idx], glyph_x_dim, glyph_y_dim);
					}
					else if(show_glyph_type == CONES)
					{
						glMatrixMode(GL_PROJECTION);						
						glLoadIdentity ();
						glOrtho(0, winWidth, 0, winHeight, -100, 100);
						glDisable(GL_DEPTH_TEST);
						glViewport(0,0,winWidth,winHeight);
						
						glClear(GL_DEPTH_BUFFER_BIT);
						glEnable(GL_LIGHTING);
						glEnable(GL_LIGHT0);
						glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
						glEnable(GL_DEPTH_TEST);
						glEnable(GL_COLOR_MATERIAL);
		
						glMatrixMode(GL_MODELVIEW);
						draw_cones(i, j, vx[idx], vy[idx], glyph_x_dim, glyph_y_dim);
						
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
						draw_ellipsoids(i, j, vx[idx], vy[idx], glyph_x_dim, glyph_y_dim);
						
						//disable lightning
						glDisable(GL_LIGHTING);
					}
				}
				
				if(show_glyph_attribute == FORCE)
				{					
					//glyphs are colored based on velocity direction
					if(color_dir == TRUE)
					{
						direction_to_color(fx[idx], fy[idx], 0, GLYPH_TYPE);
					}
					else
					{
						//set the glyphs' color based on selected scalar field
						if(show_matter_attribute == DENSITY)
						{
							set_colormap(rho[idx]);
						}
						else if(show_matter_attribute == VELOCITY_MAGNITUDE)
						{
							magnitude = vec_magnitude( vx[idx], vy[idx]);
							magnitude *= V_SCALE_FACTOR;
						
							set_colormap(magnitude);
						}
						else if(show_matter_attribute == FORCE_MAGNITUDE)
						{
							magnitude = vec_magnitude( fx[idx], fy[idx]);
							magnitude *= F_SCALE_FACTOR;
						
							set_colormap(magnitude);
						}
					}
					
					if(show_glyph_type == HEDGEHOGS)
					{
						glMatrixMode(GL_MODELVIEW);
						draw_hedgehog(i, j, fx[idx], fy[idx], glyph_x_dim, glyph_y_dim);
					}
					else if(show_glyph_type == CONES)
					{
						glMatrixMode(GL_PROJECTION);						
						glLoadIdentity ();
						glOrtho(0, winWidth, 0, winHeight, -100, 100);
						glDisable(GL_DEPTH_TEST);
						glViewport(0,0,winWidth,winHeight);
						
						glClear(GL_DEPTH_BUFFER_BIT);
						glEnable(GL_LIGHTING);
						glEnable(GL_LIGHT0);
						glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
						glEnable(GL_DEPTH_TEST);
						glEnable(GL_COLOR_MATERIAL);
		
						glMatrixMode(GL_MODELVIEW);
						draw_cones(i, j, fx[idx], fy[idx], glyph_x_dim, glyph_y_dim);
						
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
						draw_ellipsoids(i, j, fx[idx], fy[idx], glyph_x_dim, glyph_y_dim);
						
						//disable lightning
						glDisable(GL_LIGHTING);
					}
				}
			}
		
		//draw color legend for glyphs
		if(color_dir == TRUE)
		{
			draw_glyph_color_legend();
		}
		else
		{
			draw_matter_color_legend();
		}
	}
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
		case 'a': frozen = 1-frozen; break;
		case 't': dt -= 0.001; break;
		case 'T': dt += 0.001; break;
		case 'c': color_dir = 1 - color_dir; break; 
		case 'S': vec_scale *= 1.2; break;
		case 's': vec_scale *= 0.8; break;
		case 'V': visc *= 5; break;
		case 'v': visc *= 0.2; break;
		case 'x':
			draw_matter = 1 - draw_matter;
			if (draw_matter == 0) draw_vecs = 1;
			break;
		case 'y':
			draw_vecs = 1 - draw_vecs;
			if (draw_vecs==0) draw_matter = 1;
			break;
		case 'm': matter_color_type++;
			if (matter_color_type > SATURATION)
			{
				matter_color_type = BLACKWHITE;
			}
			break;
		case 'M': glyph_color_type++;
			if (glyph_color_type > SATURATION)
			{
				glyph_color_type = BLACKWHITE;
			}
			break;
		case 'b': matter_color_bands *= 2;
			if(matter_color_bands > MAX_COLOR_BANDS)
			{
				matter_color_bands = 2;
			}
			break;
		case 'B': glyph_color_bands *= 2;
			if(glyph_color_bands > MAX_COLOR_BANDS)
			{
				glyph_color_bands = 2;
			}
			break;
		case 'n':
			show_matter_attribute++;
			if (show_matter_attribute > FORCE_MAGNITUDE)
			{
				show_matter_attribute = DENSITY;
			}
			break;
		case 'N':
			show_glyph_attribute++;
			if (show_glyph_attribute > FORCE)
			{
				show_glyph_attribute = VELOCITY;
			}
			break;
		case 'k': matter_scale = 1-matter_scale; 
			if(matter_scale)
			{
				matter_clamp = 0;
				printf("Set smoke color scale interval range\n");
				printf("Insert lower value and upper value:\n");
				scanf("%f %f", &matter_scale_lmin, &matter_scale_lmax);
			}
			break;
		case 'K': glyph_scale = 1-glyph_scale; 
			if(glyph_scale)
			{
				glyph_clamp = 0;
				printf("Set glyph color scale interval range\n");
				printf("Insert lower value and upper value:\n");
				scanf("%f %f", &glyph_scale_lmin, &glyph_scale_lmax);
			}
			break;
		case 'j': matter_clamp = 1-matter_clamp; 
			if(matter_clamp)
			{
				matter_scale = 0;
				printf("Set smoke color clamp interval range\n");
				printf("Insert lower value and upper value:\n");
				scanf("%f %f", &matter_clamp_lmin, &matter_clamp_lmax);
			}
			break;
		case 'J': glyph_clamp = 1-glyph_clamp; 
			if(glyph_clamp)
			{
				glyph_scale = 0;
				printf("Set glyph color clamp interval range\n");
				printf("Insert lower value and upper value:\n");
				scanf("%f %f", &glyph_clamp_lmin, &glyph_clamp_lmax);
			}
			break;	
		case 'l': matter_enable_hsv = 1-matter_enable_hsv; 
			if(matter_enable_hsv)
			{
				printf("Increase smoke h,s,v values by:\n");
				scanf("%f %f %f", &matter_scale_h, &matter_scale_s, &matter_scale_v);
			}
			break;
		case 'L': glyph_enable_hsv = 1-glyph_enable_hsv; 
			if(glyph_enable_hsv)
			{
				printf("Increase glyph h,s,v values by:\n");
				scanf("%f %f %f", &glyph_scale_h, &glyph_scale_s, &glyph_scale_v);
			}
			break;
		case 'z':
			show_glyph_type++;
			if (show_glyph_type > ELLIPSES)
			{
				show_glyph_type = HEDGEHOGS;
			}
			break;
		case 'd': no_glyphs_x -= 5;
			if(no_glyphs_x == 30)
			{
				no_glyphs_x = 50;
			}
		case 'D': no_glyphs_y -= 5;
			if(no_glyphs_y == 30)
			{
				no_glyphs_y = 50;
			}
			break;    
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
	printf("c:     toggle direction coloring on/off\n");
	printf("S/s:   increase/decrease glyph scaling\n");
	printf("V/v:   increase decrease fluid viscosity\n");
	printf("x:     toggle drawing smoke on/off\n");
	printf("y:     toggle drawing glyphs on/off\n");
	printf("a:     toggle the animation on/off\n");
	printf("m:     toggle thru scalar coloring\n");
	printf("M:     toggle thru glyphs coloring (only when direction coloring is on)\n");
	printf("n:     toggle thru scalar shown attributes (density, velocity magnitude, force magnitude) \n");
	printf("N:     toggle thru vector shown attributes (velocity, force)\n");
	printf("b:     change color bands for smoke: 2->256\n");
	printf("B:     change color bands for glyphs: 2->256 (only when direction coloring is on)\n");
	printf("j:     toggle the smoke color map clamping on/off\n");
	printf("J:     toggle the glyph color map clamping on/off (only when direction coloring is on)\n");
	printf("k:     toggle the smoke color map scaling on/off\n");
	printf("K:     toggle the glyph color map scaling on/off (only when direction coloring is on)\n");
	printf("l:     toggle the smoke hsv coloring on/off\n");
	printf("L:     toggle the glyph hsv coloring on/off (only when direction coloring is on)\n");
	printf("z:     toggle thru glyphs types\n");
	printf("d/D:   toggle thru glyphs number on x/y axis\n");
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
