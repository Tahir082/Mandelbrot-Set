
#include "bitmap.h"

#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <pthread.h>


struct params
{
	double xmin, xmax, ymin, ymax;
	int smax, swidth, sheight;
	int width_start, width_end; //for calculating in different threads
	struct bitmap *bm;
};


int iteration_to_color( int i, int max );
int iterations_at_point( double x, double y, int max );
void compute_image( struct bitmap *bm, double xmin, double xmax, double ymin, double ymax, int max );
void *mandelthread(void *arg);

void show_help()
{
	printf("Use: mandel [options]\n");
	printf("Where options are:\n");
	printf("-m <max>    The maximum number of iterations per point. (default=1000)\n");
	printf("-x <coord>  X coordinate of image center point. (default=0)\n");
	printf("-y <coord>  Y coordinate of image center point. (default=0)\n");
	printf("-s <scale>  Scale of the image in Mandlebrot coordinates. (default=4)\n");
	printf("-W <pixels> Width of the image in pixels. (default=500)\n");
	printf("-H <pixels> Height of the image in pixels. (default=500)\n");
	printf("-o <file>   Set output file. (default=mandel.bmp)\n");
	printf("-h          Show this help text.\n");
	printf("\nSome examples are:\n");
	printf("mandel -x -0.5 -y -0.5 -s 0.2\n");
	printf("mandel -x -.38 -y -.665 -s .05 -m 100\n");
	printf("mandel -x 0.286932 -y 0.014287 -s .0005 -m 1000\n\n");
}

int main( int argc, char *argv[] )
{
	char c;

	// These are the default configuration values used
	// if no command line arguments are given.

	const char *outfile = "mandel.bmp";
	double xcenter = 0.286932;
	double ycenter = 0.014287;
	double scale = 0.0005;
	int    image_width = 4500;
	int    image_height = 4500;
	int    max = 25000;

	// For each command line argument given,
	// override the appropriate configuration value.

	while((c = getopt(argc,argv,"x:y:s:W:H:m:o:h"))!=-1) {
		switch(c) {
			case 'x':
				xcenter = atof(optarg);
				break;
			case 'y':
				ycenter = atof(optarg);
				break;
			case 's':
				scale = atof(optarg);
				break;
			case 'W':
				image_width = atoi(optarg);
				break;
			case 'H':
				image_height = atoi(optarg);
				break;
			case 'm':
				max = atoi(optarg);
				break;
			case 'o':
				outfile = optarg;
				break;
			case 'h':
				show_help();
				exit(1);
				break;
		}
	}

	// Display the configuration of the image.
	printf("mandel: x=%lf y=%lf scale=%lf max=%d outfile=%s\n",xcenter,ycenter,scale,max,outfile);

	// Create a bitmap of the appropriate size.
	struct bitmap *bm = bitmap_create(image_width,image_height);

	// Fill it with a dark blue, for debugging
	bitmap_reset(bm,MAKE_RGBA(0,0,255,0));

	// Compute the Mandelbrot image
	compute_image(bm,xcenter-scale,xcenter+scale,ycenter-scale,ycenter+scale,max);

	// Save the image in the stated file.
	if(!bitmap_save(bm,outfile)) {
		fprintf(stderr,"mandel: couldn't write to %s: %s\n",outfile,strerror(errno));
		return 1;
	}

	return 0;
}

/*
Compute an entire Mandelbrot image, writing each point to the given bitmap.
Scale the image to the range (xmin-xmax,ymin-ymax), limiting iterations to "max"
*/

void compute_image( struct bitmap *bm, double xmin, double xmax, double ymin, double ymax, int max )
{
	pthread_t th1, th2, th3, th4;
	struct params p1, p2, p3, p4;
	int width = bitmap_width(bm);
	int height = bitmap_height(bm);
	p1.swidth = width;
	p2.swidth = width;
	p3.swidth = width;
	p4.swidth = width;
	p1.sheight = height;
	p2.sheight = height;
	p3.sheight = height;
	p4.sheight = height;	
	p1.width_start = 0;
	p2.width_start = width/4;
	p3.width_start = width/2;
	p4.width_start = 3*(width)/4;
	p1.xmin = xmin;
	p1.xmax = xmax;
	p1.ymin = ymin;
	p1.ymax = ymax;
	p1.smax = max;
	p1.bm = bm;
	p1.width_end = width/4;
	p2.xmin = xmin;
	p2.xmax = xmax;
	p2.ymin = ymin;
	p2.ymax = ymax;
	p2.smax = max;
	p2.bm = bm;
	p2.width_end= width/2;
	p3.xmin = xmin;
	p3.xmax = xmax;
	p3.ymin = ymin;
	p3.ymax = ymax;
	p3.smax = max;
	p3.bm = bm;
	p3.width_end = 3*(width)/4;
	p4.xmin = xmin;
	p4.xmax = xmax;
	p4.ymin = ymin;
	p4.ymax = ymax;
	p4.smax = max;
	p4.bm = bm;
	p4.width_end = width;

	
	printf("creating thread   \n");
	pthread_create(&th1, NULL, mandelthread, &p1);
	pthread_create(&th2, NULL, mandelthread, &p2);
	pthread_create(&th3, NULL, mandelthread, &p3);
	pthread_create(&th4, NULL, mandelthread, &p4);
	pthread_join(th1, NULL);
	pthread_join(th2, NULL);
	pthread_join(th3, NULL);
	pthread_join(th4, NULL);
}


void *mandelthread(void *arg)
{
	int i,j;
	struct params *c = (struct params*)arg;



	// For every pixel in the image...

	for(j=c->width_start; j<c->width_end;j++) 
	{

		for(i=0;i<c->swidth;i++) 
		{

			// Determine the point in x,y space for that pixel.
			double x = c->xmin + i*(c->xmax - c->xmin)/c->swidth;
			double y = c->ymin + j*(c->ymax - c->ymin)/c->sheight;

			// Compute the iterations at that point.
			int iters = iterations_at_point(x,y,c->smax);

			// Set the pixel in the bitmap.
			bitmap_set(c->bm,i,j,iters);
		}
	}
}

/*
Return the number of iterations at point x, y
in the Mandelbrot space, up to a maximum of max.
*/

int iterations_at_point( double x, double y, int max )
{
	double x0 = x;
	double y0 = y;

	int iter = 0;

	while( (x*x + y*y <= 4) && iter < max ) {

		double xt = x*x - y*y + x0;
		double yt = 2*x*y + y0;

		x = xt;
		y = yt;

		iter++;
	}

	return iteration_to_color(iter,max);
}

/*
Convert a iteration number to an RGBA color.
Here, we just scale to gray with a maximum of imax.
Modify this function to make more interesting colors.
*/

int iteration_to_color( int i, int max )
{
	int color1, color2 , color3;
	
        	if (i < 255) {
            		color1 = 16*i;
            		color2 = 8*(255 - i)/max;
           		color3 = 2;
        	} else if (i < 128 ) {
            		color1 = 16*(128 - i)/max;
            		color2 = 8*i;
            		color3 = 4;
        	} else if (i < 64) {
            		color1 = 16*(i - 64);
            		color2 = 8*(64- i)/max;
            		color3 = 6;
        	} else {
            		color1 = i*16;
            		color2 = i*8;
            		color3 = 3;
        	}
	return MAKE_RGBA(color1, color2, color3, 0);
}




