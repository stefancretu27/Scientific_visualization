// Usage: Drag with the mouse to add smoke to the fluid. This will also move a "rotor" that disturbs
//        the velocity field at the mouse location. Press the indicated keys to change options
//--------------------------------------------------------------------------------------------------

//#include <rfftw.h>              //the numerical simulation FFTW 
#include <stdio.h>              //for printing the help text

#include <math.h>               //for various math functions #include
#include <GL/glut.h>            //the GLUT graphics library

#include <visualization.h>

//--- SIMULATION PARAMETERS ------------------------------------------------------------------------

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
unsigned int winWidth, winHeight;      									//size of the graphics window, in pixels
float vec_scale = 1000;													//scaling of hedgehogs
bool frozen = FALSE;   													//toggles on/off the animation
bool draw_glyph_grid = FALSE;
visualization_technique vis_tech = MATTER;

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

//streamlines-related data
int_coord *seed;
int seeds_count = -1;
unsigned short int streamline_length_scaling_factor = 5;
streamlines_type streamline_type = POINTS;


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
	
	seed = (int_coord*)calloc(MAX_SEEDS, sizeof(int_coord));
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
		compute_attributes_min_max(rho, vx, vy, fx, fy, &matter_scale_lmin, &matter_scale_lmax);

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
		//glClearColor(0.6, 0.4, 0.2, 0);							//set background color
		grayscale(scalar_value, &r, &g, &b);
	}
	else if (matter_color_type == RAINBOW)
	{
		//glClearColor(0.0, 0.0, 0.0, 0);							//set background color
		rainbow(scalar_value, &r, &g, &b);
	}
	else if (matter_color_type == SATURATION)
	{
		//glClearColor(0.2, 0.0, 0.2, 0);							//set background color
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
void direction_to_color(vector vec, float value, unsigned short int type)
{
	float r, g, b, h, s, v, color_value;
	
	//if function is called for drawing glyphs
	if(type == GLYPH_TYPE)
	{
		//compute the value whose color is to be set
		color_value = atan2(vec.y, vec.x)/M_PI + 1;
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
		//glClearColor(0.0, 0.0, 0.0, 0);
		rainbow(color_value, &r, &g, &b);
	}
	else if (glyph_color_type == SATURATION)
	{
		//lClearColor(0.0, 0.2, 0.2, 0);
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
	//compute grid's cells' width and height, counted as pixels
	fftw_real wn = (fftw_real)winWidth / (fftw_real)(DIM + 1);
	fftw_real hn = (fftw_real)winHeight / (fftw_real)(DIM + 1);
	
	//draw smoke/fluid/matter
	if (vis_tech == MATTER)
	{	
		//draw the matter
		draw_matter_fluid_smoke(wn, hn, rho, vx, vy, fx, fy);
		
		//draw color legend for matter
		draw_matter_color_legend();
	}
	
	//draw vectors
	if (vis_tech == GLYPHS)		
	{
		//variables used to create grid cells and to iterate through them
		grid_cell cell;
		vector cell_iterator;
		cell_sample_points cell_corners_indeces;
		
		//variables used to iterate through the data and to store the data for a given cell
		int i, j, idx;
		
		//needed for storing the interpolation results for each cell
		fftw_real density;
		vector velocity, force;
		
		//needed to store the gradient computation result
		vector gradient;
		
		//compute grid's' cells' width and height, with respect to the initial data samples grid. 
		//It is a vector that shows how much the bottom left corner has moved from the original position on both axes.
		vector step;
		step.x = (fftw_real) DIM/no_glyphs_x;
		step.y = (fftw_real) DIM/no_glyphs_y;
		
		//compute the ratio with which the glyphs increase. Shows how much the lower left corner has moved from the initial position.
		point data_point;
		data_point.r = (step.x - 1);
		data_point.s = (step.y - 1);
						
		for (cell_iterator.x = 0; cell_iterator.x < DIM - 0.1; cell_iterator.x += step.x)
			for (cell_iterator.y = 0; cell_iterator.y < DIM - 0.1; cell_iterator.y += step.y)									
			{
				i = (int)cell_iterator.x;
				j = (int)cell_iterator.y;
				
				/*
				 * Step 1: initialize cell grid. It's used for positioning glyphs on the visualization window.
				 */
				initialize_cell(&cell, wn, hn, cell_iterator, step);
				
				//Eventually draw the cell grid
				if(draw_glyph_grid)
				{
					draw_cell(&cell);
				}
				
				/*
				 * Step 2: Compute the attributes values in the given cell grid. 
				 * Each grid cell displays a glyph which shows the value for the bottom left corner point, where data is measured.
				 * If #glyphs does not change => no interpolation. Otherwise the data point changes it's position and its value is computed using interpolation
				 */
				if( step.x == 1 && step.y == 1)
				{
					idx = (j * DIM) + i;
					
					density = rho[idx];
					velocity.x = vx[idx];
					velocity.y = vy[idx];
					force.x = fx[idx];
					force.y = fy[idx];
				}
				else
				{
					//It is assumed that when the glyph increase in size, they don't exceed more than 4 grid cells from the initial/original/default grid
					cell_corners_indexing(i, j, &cell_corners_indeces);
					
					//interpolate for rho
					density = bilinear_interpolation(cell_corners_indeces, data_point, rho);
					//interpolate for vx
					velocity.x = bilinear_interpolation(cell_corners_indeces, data_point, vx);
					//interpolate for vy
					velocity.y = bilinear_interpolation(cell_corners_indeces, data_point, vy);
					//interpolate for vx
					force.x = bilinear_interpolation(cell_corners_indeces, data_point, fx);
					//interpolate for vy
					force.y = bilinear_interpolation(cell_corners_indeces, data_point, fy);	
				}

				/*
				 * Step 3: choose vectorial attribute: velocity, force, gradients
				 */
				if(show_glyph_attribute == VELOCITY)
				{	
					//Part 1: choose glyphs' color
					if(color_dir == TRUE)
					{
						//glyphs are colored based on velocity direction
						direction_to_color(velocity, 0, GLYPH_TYPE);
					}
					else
					{
						//set the glyphs' color based on selected scalar field
						glyphs_scalar_color(density, velocity, force);
					}
				
					//Part 2: choose glyps' type. Don't overlap with the color legend
					if(cell_iterator.y >= 2)
					{ 
						draw_glyphs(cell, velocity, vec_scale);
					}
				}
				else if(show_glyph_attribute == FORCE)
				{	
					//Part1: choose glyps'color				
					if(color_dir == TRUE)
					{
						//glyphs are colored based on force direction
						direction_to_color(force, 0, GLYPH_TYPE);
					}
					else
					{
						//set the glyphs' color based on selected scalar field
						glyphs_scalar_color(density, velocity, force);
					}
					
					//Part2: choose glyps' type
					if(cell_iterator.y >= 2)
					{
						draw_glyphs(cell, force, vec_scale);
					}
				}
				else if(show_glyph_attribute == DENSITY_GRADIENT && color_dir == FALSE)
				{
					//Part 1: compute gradient
					compute_density_gradient(i, j, step, data_point, rho, &gradient);
					gradient.x *= 5; gradient.y *= 5;
					
					//Part 2: choose glyphs color: set the glyphs' color based on selected scalar field
					glyphs_scalar_color(density, velocity, force);
					
					//Part 3: choose glyphs' type
					if(cell_iterator.y >= 2)
					{
						draw_glyphs(cell, gradient, vec_scale);
					}
				}
				else if(show_glyph_attribute == VEL_MAGN_GRADIENT && color_dir == FALSE)
				{
					//Part 1: compute gradient
					compute_velocity_magnitude_gradient(i, j, step, data_point, vx, vy, &gradient);
					gradient.x *= 5; gradient.y *= 5;
					
					//Part 2: choose glyphs color: set the glyphs' color based on selected scalar field
					glyphs_scalar_color(density, velocity, force);
					
					//Part 3: choose glyphs' type
					if(cell_iterator.y >= 2)
					{
						draw_glyphs(cell, gradient, vec_scale);
					}
				}
			}
		
		/*
		 * Step 4: Draw color legend for glyphs
		 */
		if(color_dir == TRUE)
		{
			draw_glyph_color_legend(winWidth, winHeight);
		}
		else
		{
			draw_matter_color_legend();
		}
	}
	
	//draw streamlines
	if(vis_tech == STREAMLINES)
	{
		//variables used to create grid cells and to iterate through them
		grid_cell cell;
		//the next 2 are used to fit into cell creating functions pattern
		vector cell_iterator;
		vector step;
		step.x = 1;
		step.y = 1;
		//vector coordinates and magnitude in the given cell
		vector current_point_velocity, next_point_velocity;
		//clicked point's coordinates within a cell in [0;1]
		point point_coord_in_cell, prev_point, current_point, next_point;
		int_coord cell_coord_in_grid;
		cell_sample_points cell_corners_indeces;

		int k;
		float norm, streamline_length;  
		float step_size = 0.5 * hn, max_length = streamline_length_scaling_factor * hn; 

		for (cell_iterator.x = 0; cell_iterator.x < DIM - 0.1; cell_iterator.x += step.x)
			for (cell_iterator.y = 0; cell_iterator.y < DIM - 0.1; cell_iterator.y += step.y)
			{
				/*
				 * Step 1: initialize cell grid. It's used for positioning glyphs on the visualization window.
				 */
				initialize_cell(&cell, wn, hn, cell_iterator, step);
				
				//Eventually draw the cell grid
				if(draw_glyph_grid && cell_iterator.y >= 2)
				{
					draw_cell(&cell);
				}
				
				/*
				 * Step 2: for each selected point in the screen, identify the cell they belong to, draw the point and set color for it.
				 * Screen's origin is in the top left corner, whereas the cell grid starts from bottom left corner.
				 */ 
				for(k = 1; k <= seeds_count; k++)
				{
					//Find the cell where a clicked point is.
					if(cell.x <= seed[k].a && seed[k].a < cell.x + cell.width)
						if(cell.y - cell.height <= seed[k].b && seed[k].b < cell.y)			//as cell indexing starts from top left corner
						{
							//Compute the coordinates of the clicked point in the given cell.
							//This helps to obtain the offset of the point related to the bottom left corner of the cell.
							point_coord_in_cell.r = seed[k].a/cell.width;
							point_coord_in_cell.s = seed[k].b/cell.height;
							//Get the cell's indeces in the grid
							cell_coord_in_grid.a = (int) point_coord_in_cell.r;
							cell_coord_in_grid.b = (int) point_coord_in_cell.s;
							cell_coord_in_grid.b = DIM - cell_coord_in_grid.b - 1; 
							//Get fractional part which shows where in the cell was the point clicked. 
							point_coord_in_cell.r -= cell_coord_in_grid.a;
							point_coord_in_cell.s -= cell_coord_in_grid.b;
							
							//Compute the indeces of the sample points placed at the corners of the cell where the point was clicked
							cell_corners_indexing(cell_coord_in_grid.a, cell_coord_in_grid.b, &cell_corners_indeces);
						
							//printf("%d %d %d %d\n", cell_corners_indeces.idx1, cell_corners_indeces.idx2, cell_corners_indeces.idx3, cell_corners_indeces.idx4);
						
							//Get the velocity's value by doing piecewise linear interpolation for the given cell
							current_point_velocity.x = bilinear_interpolation(cell_corners_indeces, point_coord_in_cell, vx);
							current_point_velocity.y = bilinear_interpolation(cell_corners_indeces, point_coord_in_cell, vy);
							
							//Create the new points in the streamline starting from the current seed point
							current_point.r = (float) seed[k].a;
							current_point.s = (float) seed[k].b;
							
							//Draw seed points 
							if(streamline_type == POINTS && cell_iterator.y >= 2)
							{
								draw_points(current_point, current_point_velocity, 0);
							}
							
							//printf("%f %f %f %f %f\n", current_point.r, current_point.s, current_point_velocity.x, current_point_velocity.y, vector_normalize(current_point_velocity));
							
							if(current_point_velocity.x != 0 || current_point_velocity.y != 0)
							{
								for(streamline_length = 0; streamline_length < max_length; streamline_length += step_size)
								{
									//Firstly, normalize the vector values at the current point
									norm = vector_normalize(current_point_velocity);
									current_point_velocity.x /= norm;
									current_point_velocity.y /= norm;
								
									//Then, compute next point's coordinates using Euclidean integration
									next_point.r = current_point.r + current_point_velocity.x*step_size;
									//because on Oy axis the direction should be inverse due to the fact that the coordinates start from top-left corner and sample points from bottom left
									next_point.s = current_point.s - current_point_velocity.y*step_size;
									
									//if the next point is in the cell grid, proceed further
									if(next_point.r >= 0 && next_point.r <= 50*wn && next_point.s >= hn && next_point.s <= 49*hn)
									{
										//the next point becomes the current point so a new point can be generated. Beforehand, the current point becomes the prev point.
										//It's original coordinates are kept as the current point's coordinates and used for drawing lines
										prev_point.r = current_point.r;
										prev_point.s = current_point.s;
										current_point.r = next_point.r;
										current_point.s = next_point.s;
								
										//Compute the coordinates of the clicked point in the given cell.
										//This helps to obtain the offset of the point related to the bottom left corner of the cell.
										next_point.r /= cell.width;
										next_point.s /= cell.height;
										//get the integer part that shows the cell's indeces in the grid starting from top left corner
										cell_coord_in_grid.a = (int) next_point.r;
										cell_coord_in_grid.b = (int) next_point.s;
										cell_coord_in_grid.b = DIM - cell_coord_in_grid.b - 1; 
										//get fractional part which shows where in the cell was the point clicked. 
										next_point.r -= cell_coord_in_grid.a;
										next_point.s -= cell_coord_in_grid.b;
								
										//Compute the indeces of the sample points placed at the corners of the cell where the point was clicked
										cell_corners_indexing(cell_coord_in_grid.a, cell_coord_in_grid.b, &cell_corners_indeces);
									
										//printf("%d %d %d %d\n", cell_corners_indeces.idx1, cell_corners_indeces.idx2, cell_corners_indeces.idx3, cell_corners_indeces.idx4);
									
										//Compute the velocity's value by doing piecewise linear interpolation for the given cell
										next_point_velocity.x = bilinear_interpolation(cell_corners_indeces, next_point, vx);
										next_point_velocity.y = bilinear_interpolation(cell_corners_indeces, next_point, vy);
										
										current_point_velocity.x = next_point_velocity.x;
										current_point_velocity.y = next_point_velocity.y;

										if(streamline_type == POINTS)
										{
											draw_points(current_point, current_point_velocity, streamline_length);
										}
										if(streamline_type == LINES)
										{
											draw_lines(prev_point, current_point, current_point_velocity);
										}
									}
								}
							}
						}
				}
			}
		draw_matter_color_legend();
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
		case 'x': vis_tech++;
			if (vis_tech > STREAMLINES) 
			{
				vis_tech = MATTER;
			}
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
			if (show_glyph_attribute > VEL_MAGN_GRADIENT)
			{
				show_glyph_attribute = VELOCITY;
			}
			break;
		case 'k': matter_scale = 1-matter_scale; 
			if(matter_scale)
			{
				matter_clamp = 0;
				printf("Set smoke color scale interval range\n");
				//printf("Insert lower value and upper value:\n");
				//scanf("%f %f", &matter_scale_lmin, &matter_scale_lmax);
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
			if(no_glyphs_x == 20)
			{
				no_glyphs_x = 50;
			}
			break;
		case 'D': no_glyphs_y -= 5;
			if(no_glyphs_y == 20)
			{
				no_glyphs_y = 50;
			}
			break;
		case 'f': streamline_length_scaling_factor += 5;
			if(streamline_length_scaling_factor == 30)
			{
				streamline_length_scaling_factor = 5;
			}
			break;
		case 'g':
			streamline_type++;
			if (streamline_type > LINES)
			{
				streamline_type = POINTS;
			}
			break;
		case '`': draw_glyph_grid = 1 - draw_glyph_grid; 
				break;    
		case 'q': exit(0);
	}
}

// drag: When the user drags with the mouse, add a force that corresponds to the direction of the mouse
// cursor movement. Also inject some new matter into the field at the mouse location.
void drag(int mx, int my)
{
	//if(vis_tech != STREAMLINES)
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
}

void on_mouse_click(int button, int state, int x, int y)
{
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN && vis_tech == STREAMLINES) 
	{ 
		if(seeds_count < MAX_SEEDS)
		{
			//store the x,y value where the click happened
			seeds_count++;
			seed[seeds_count].a = x;
			seed[seeds_count].b = y;
		
			//printf("on_mouse_click: %d %d %d %d %d \n", seeds_count, seed[seeds_count].a, seed[seeds_count].b, x, y);
		}
		else
		{
			seeds_count = 0;
		}
	}
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
	printf("x:     toggle thru visualization techniques: matter/glyphs/streamlines\n");
	printf("a:     toggle the animation on/off\n");
	printf("m/M:   toggle thru coloring of matter/direction (only when direction coloring is on)\n");
	printf("b/B:   change color bands for smoke/ glyphs (only when direction coloring is on):2->256\n");
	printf("n:     toggle thru scalar attributes (density, velocity magnitude, force magnitude) \n");
	printf("N:     toggle thru vector attributes (velocity, force, rho gradient, vel magnitude gradient)\n");
	printf("j/J:   toggle the colormap clamping on/off for matter/direction (only when direction coloring is on)\n");
	printf("k/K:   toggle the colormap scaling on/off for matter/direction (only when direction coloring is on)\n");
	printf("l/L:   toggle the hsv coloring on/off for matter/direction (only when direction coloring is on)\n");
	printf("z:     toggle thru glyphs types\n");
	printf("d/D:   toggle thru glyphs number on x/y axis\n");
	printf("f:     change streamline's length\n");
	printf("g:     toggle thru streamlines' rendering type\n");
	printf("q:     quit\n\n");

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize(1300,800);
	glutCreateWindow("Real-time smoke simulation and visualization");
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutIdleFunc(do_one_simulation_step);
	glutKeyboardFunc(keyboard);
	glutMouseFunc(on_mouse_click);
	glutMotionFunc(drag);
	
	init_simulation(DIM);	//initialize the simulation data structures
	glutMainLoop();			//calls do_one_simulation_step, keyboard, display, drag, reshape
	return 0;
}
