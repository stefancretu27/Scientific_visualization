/*
 * Helper functions 
 */
void clamp_value_to_01(float *value);
void draw_triangle_vertex( double x, double y, float color_value);
float direction2angle(float x, float y);

/* 
 * Functions implemented for step 2. 
 */
//color mapping functions
void grayscale(float value, float *R, float *G, float *B);
void rainbow(float value, float* R, float* G, float* B);
void green_saturation(float value, float *R, float *G, float *B);
void red_saturation(float value, float *R, float *G, float *B);

//functions used for setting colors
extern void set_colormap(float scalar_value);
extern void direction_to_color(float x, float y, float value, unsigned short int type);

//functions used for parameterize the colormap
float set_color_bands(float value, unsigned short int type);
void rgb2hsv(float r, float g , float b, float *h, float *s, float *v);
void hsv2rgb(float h, float s, float v, float *r, float *g, float *b);

//scale the color map
float set_scale(float value, float x0, float x1,  float y0, float y1);
float clamp_value( float value, float min, float max);

//color legends
void draw_matter_color_legend();
void draw_glyph_color_legend();

/* 
 * Functions implemented for step 3. 
 */
void draw_hedgehog(int x, int y, float attribute_x, float attribute_y, float cell_width, float cell_height);
//Draw 2D glyphs
void draw_triangle(int x, int y, float attribute_x, float attribute_y, float cell_width, float cell_height);
void draw_filled_ellipse(float x, float y, float attribute_x, float attribute_y, float cell_width, float cell_height);
//Draw 3D glyphs
void draw_cones(int x, int y, float attribute_x, float attribute_y, float cell_width, float cell_height);	
void draw_ellipsoids(int x, int y, float attribute_x, float attribute_y, float cell_width, float cell_height);
