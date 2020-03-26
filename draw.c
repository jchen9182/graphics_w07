#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ml6.h"
#include "display.h"
#include "draw.h"
#include "matrix.h"

/*======== void add_box() ==========
  Inputs:   struct matrix * edges
            double x
            double y
            double z
            double width
            double height
            double depth
  add the points for a rectagular prism whose
  upper-left-front corner is (x, y, z) with width,
  height and depth dimensions.
  ====================*/
void add_box( struct matrix * edges,
              double x, double y, double z,
              double width, double height, double depth ) {
    // Left Face
    add_edge(edges, x, y, z, x, y - height, z);
    add_edge(edges, x, y - height, z, x, y - height, z - depth);
    add_edge(edges, x, y - height, z - depth, x, y, z - depth);
    add_edge(edges, x, y, z - depth, x, y, z);

    // Right Face
    add_edge(edges, x + width, y, z, x + width, y - height, z);
    add_edge(edges, x + width, y - height, z, x + width, y - height, z - depth);
    add_edge(edges, x + width, y - height, z - depth, x + width, y, z - depth);
    add_edge(edges, x + width, y, z - depth, x + width, y, z);

    // Connecting Edges
    add_edge(edges, x, y, z, x + width, y, z);
    add_edge(edges, x, y - height, z, x + width, y - height, z);
    add_edge(edges, x, y - height, z - depth, x + width, y - height, z - depth);
    add_edge(edges, x, y, z - depth, x + width, y, z - depth);
}

/*======== void add_sphere() ==========
  Inputs:   struct matrix * points
            double cx
            double cy
            double cz
            double r
            int step
  adds all the points for a sphere with center (cx, cy, cz)
  and radius r using step points per circle/semicircle.
  Since edges are drawn using 2 points, add each point twice,
  or add each point and then another point 1 pixel away.
  should call generate_sphere to create the necessary points
  ====================*/
void add_sphere( struct matrix * edges,
                 double cx, double cy, double cz,
                 double r, int step ) {
    struct matrix * sphere = generate_sphere(cx, cy, cz, r, step);
    double ** matrix = sphere -> m;
    int cols = sphere -> cols;

    for (int col = 0; col < cols - 1; col++) {
        double x = matrix[0][col];
        double y = matrix[1][col];
        double z = matrix[2][col];
        
        add_edge(edges, x, y, z, x + 1, y + 1, z + 1);
    }
}

/*======== void generate_sphere() ==========
  Inputs:   struct matrix * points
            double cx
            double cy
            double cz
            double r
            int step
  Returns: Generates all the points along the surface
           of a sphere with center (cx, cy, cz) and
           radius r using step points per circle/semicircle.
           Returns a matrix of those points
  ====================*/
struct matrix * generate_sphere(double cx, double cy, double cz,
                                double r, int step ) {
    struct matrix * points = new_matrix(4, 100);
    double phi = 0;
    double theta = 0;

    for (int p = 0; p < 1000; p += step) {
        phi = (2 * M_PI) * p / 1000;

        for (int t = 0; t < 1000; t += step) {
            theta = M_PI * t / 1000;

            double x = r * cos(theta) + cx;
            double y = r * sin(theta) * cos(phi) + cy;
            double z = r * sin(theta) * sin(phi) + cz;

            add_point(points, x, y, z);
        }
    }

    return points;
}

/*======== void add_torus() ==========
  Inputs:   struct matrix * points
            double cx
            double cy
            double cz
            double r1
            double r2
            double step
  Returns:
  adds all the points required for a torus with center (cx, cy, cz),
  circle radius r1 and torus radius r2 using step points per circle.
  should call generate_torus to create the necessary points
  ====================*/
void add_torus( struct matrix * edges,
                double cx, double cy, double cz,
                double r1, double r2, int step ) {
    struct matrix * torus = generate_torus(cx, cy, cz, r1, r2, step);
    double ** matrix = torus -> m;
    int cols = torus -> cols;

    for (int col = 0; col < cols - 1; col++) {
        double x = matrix[0][col];
        double y = matrix[1][col];
        double z = matrix[2][col];
        
        add_edge(edges, x, y, z, x, y, z);
    }
}

/*======== void generate_torus() ==========
  Inputs:   struct matrix * points
            double cx
            double cy
            double cz
            double r
            int step
  Returns: Generates all the points along the surface
           of a torus with center (cx, cy, cz),
           circle radius r1 and torus radius r2 using
           step points per circle.
           Returns a matrix of those points
  ====================*/
struct matrix * generate_torus( double cx, double cy, double cz,
                                double r1, double r2, int step ) {
    struct matrix * points = new_matrix(4, 100);
    double phi = 0;
    double theta = 0;

    for (int p = 0; p < 1000; p += step) {
        phi = (2 * M_PI) * p / 1000;

        for (int t = 0; t < 1000; t += step) {
            theta = (2 * M_PI) * t / 1000;

            double x = cos(phi) * (r1 * cos(theta) + r2) + cx;
            double y = r1 * sin(theta) + cy;
            double z = -1 * sin(phi) * (r1 * cos(theta) + r2) + cz;

            add_point(points, x, y, z);
        }
    }

    return points;
}

/*======== void add_circle() ==========
  Inputs:   struct matrix * edges
            double cx
            double cy
            double r
            double step
  Adds the circle at (cx, cy) with radius r to edges
  ====================*/
void add_circle( struct matrix *edges,
                 double cx, double cy, double cz,
                 double r, double step ) {
    double angle = 0;
    double x0 = cx + r;
    double y0 = cy;

    for (double t = 0; t < 1; t += step) {
        double x1 = r * cos(angle) + cx;
        double y1 = r * sin(angle) + cy;

        add_edge(edges, x0, y0, 0, x1, y1, 0);
        x0 = x1;
        y0 = y1;
        angle += (2 * M_PI) * step;
    }
}

/*======== void add_curve() ==========
Inputs:   struct matrix *edges
         double x0
         double y0
         double x1
         double y1
         double x2
         double y2
         double x3
         double y3
         double step
         int type
Adds the curve bounded by the 4 points passsed as parameters
of type specified in type (see matrix.h for curve type constants)
to the matrix edges
====================*/
void add_curve( struct matrix *edges,
                double x0, double y0,
                double x1, double y1,
                double x2, double y2,
                double x3, double y3,
                double step, int type ) {
    struct matrix * xco = generate_curve_coefs(x0, x1, x2, x3, type);
    struct matrix * yco = generate_curve_coefs(y0, y1, y2, y3, type);
    double ** xm = xco -> m;
    double ** ym = yco -> m;
    double xold = x0;
    double yold = y0;

    for (double t = 0; t < 1; t += step) {
        double xnew, ynew;
        xnew = t * (t * (xm[0][0] * t + xm[1][0]) + xm[2][0]) + xm[3][0];
        ynew = t * (t * (ym[0][0] * t + ym[1][0]) + ym[2][0]) + ym[3][0];

        add_edge(edges, xold, yold, 0, xnew, ynew, 0);
        xold = xnew;
        yold = ynew;
    }
}

/*======== void add_point() ==========
Inputs:   struct matrix * points
int x
int y
int z
Returns:
adds point (x, y, z) to points and increment points.lastcol
if points is full, should call grow on points
====================*/
void add_point( struct matrix * points, double x, double y, double z) {
    double ** matrix = points -> m;
    int cols = points -> cols;
    int lastcol = points -> lastcol;

    if (lastcol == cols) grow_matrix(points, cols + 100);

    matrix[0][lastcol] = x;
    matrix[1][lastcol] = y;
    matrix[2][lastcol] = z;
    matrix[3][lastcol] = 1;

    points -> lastcol++;
} //end add_point

/*======== void add_edge() ==========
Inputs:   struct matrix * points
int x0, int y0, int z0, int x1, int y1, int z1
Returns:
add the line connecting (x0, y0, z0) to (x1, y1, z1) to points
should use add_point
====================*/
void add_edge( struct matrix * points,
    double x0, double y0, double z0,
    double x1, double y1, double z1) {
        add_point(points, x0, y0, z0);
        add_point(points, x1, y1, z1);
}

/*======== void draw_lines() ==========
Inputs:   struct matrix * points
screen s
color c
Returns:
Go through points 2 at a time and call draw_line to add that line
to the screen
====================*/
void draw_lines( struct matrix * points, screen s, color c) {
    int lastcol = points -> lastcol;

    if (lastcol < 2) {
        printf("Need at least 2 points to draw a line!\n");
        return;
    }

    for (int point = 0; point < lastcol - 1; point += 2) {
        int x0 = points -> m[0][point];
        int y0 = points -> m[1][point];
        int x1 = points -> m[0][point+1];
        int y1 =points -> m[1][point+1];

        draw_line(x0, y0, x1, y1, s, c);
    }
}// end draw_lines

void draw_line(int x0, int y0, int x1, int y1, screen s, color c) {
    int x, y, d, A, B;
    //swap points if going right -> left
    int xt, yt;
    if (x0 > x1) {
        xt = x0;
        yt = y0;
        x0 = x1;
        y0 = y1;
        x1 = xt;
        y1 = yt;
    }

    x = x0;
    y = y0;
    A = 2 * (y1 - y0);
    B = -2 * (x1 - x0);

    //octants 1 and 8
    if ( abs(x1 - x0) >= abs(y1 - y0) ) {

        //octant 1
        if ( A > 0 ) {

            d = A + B/2;
            while ( x < x1 ) {
                plot( s, c, x, y );
                if ( d > 0 ) {
                    y+= 1;
                    d+= B;
                }
                x++;
                d+= A;
            } //end octant 1 while
            plot( s, c, x1, y1 );
        } //end octant 1

        //octant 8
        else {
            d = A - B/2;

            while ( x < x1 ) {
                //printf("(%d, %d)\n", x, y);
                plot( s, c, x, y );
                if ( d < 0 ) {
                    y-= 1;
                    d-= B;
                }
                x++;
                d+= A;
            } //end octant 8 while
            plot( s, c, x1, y1 );
        } //end octant 8
    }//end octants 1 and 8

    //octants 2 and 7
    else {

        //octant 2
        if ( A > 0 ) {
            d = A/2 + B;

            while ( y < y1 ) {
                plot( s, c, x, y );
                if ( d < 0 ) {
                    x+= 1;
                    d+= A;
                }
                y++;
                d+= B;
            } //end octant 2 while
            plot( s, c, x1, y1 );
        } //end octant 2

        //octant 7
        else {
            d = A/2 - B;

            while ( y > y1 ) {
                plot( s, c, x, y );
                if ( d > 0 ) {
                    x+= 1;
                    d+= A;
                }
                y--;
                d-= B;
            } //end octant 7 while
            plot( s, c, x1, y1 );
        } //end octant 7
    }//end octants 2 and 7
} //end draw_line

//======== void change_color() ==========
void change_color(color * c, int r, int g, int b) {
    c->red = r;
    c->green = g;
    c->blue = b;
}