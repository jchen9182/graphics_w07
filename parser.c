#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ml6.h"
#include "display.h"
#include "draw.h"
#include "matrix.h"
#include "parser.h"


/*======== void parse_file () ==========
Inputs:   char * filename
          struct matrix * transform,
          struct matrix * pm,
          screen s
Returns:
Goes through the file named filename and performs all of the actions listed in that file.
The file follows the following format:
     Every command is a single character that takes up a line
     Any command that requires arguments must have those arguments in the second line.
     The commands are as follows:
         sphere: add a sphere to the edge matrix -
                 takes 4 arguemnts (cx, cy, cz, r)
         torus: add a torus to the edge matrix -
                takes 5 arguemnts (cx, cy, cz, r1, r2)
         box: add a rectangular prism to the edge matrix -
              takes 6 arguemnts (x, y, z, width, height, depth)
         clear: clears the edge matrix
         circle: add a circle to the edge matrix -
                 takes 4 arguments (cx, cy, cz, r)
         hermite: add a hermite curve to the edge matrix -
                  takes 8 arguments (x0, y0, x1, y1, rx0, ry0, rx1, ry1)
         bezier: add a bezier curve to the edge matrix -
                 takes 8 arguments (x0, y0, x1, y1, x2, y2, x3, y3)
         line: add a line to the edge matrix -
               takes 6 arguemnts (x0, y0, z0, x1, y1, z1)
         ident: set the transform matrix to the identity matrix -
         scale: create a scale matrix,
                then multiply the transform matrix by the scale matrix -
                takes 3 arguments (sx, sy, sz)
         move: create a translation matrix,
               then multiply the transform matrix by the translation matrix -
               takes 3 arguments (tx, ty, tz)
         rotate: create a rotation matrix,
                 then multiply the transform matrix by the rotation matrix -
                 takes 2 arguments (axis, theta) axis should be x y or z
         apply: apply the current transformation matrix to the edge matrix
         display: clear the screen, then
                  draw the lines of the edge matrix to the screen
                  display the screen
         save: clear the screen, then
               draw the lines of the edge matrix to the screen
               save the screen to a file -
               takes 1 argument (file name)
         quit: end parsing
See the file script for an example of the file format
IMPORTANT MATH NOTE:
the trig functions int math.h use radian mesure, but us normal
humans use degrees, so the file will contain degrees for rotations,
be sure to conver those degrees to radians (M_PI is the constant
for PI)
====================*/
void parse_file(char * filename,
                struct matrix * transform,
                struct matrix * edges,
                screen s) {
    FILE *f;
    char line[256];
    clear_screen(s);
    int SIZE = 500;
    color c;
    change_color(&c, 249, 179, 164);
    double line_step = .001;
    double sphere_step = 10;

    if (strcmp(filename, "stdin") == 0) f = stdin;
    else f = fopen(filename, "r");

    char lines[SIZE][256];
    int counter = 0;
    while (fgets(line, 255, f) != NULL) {
        line[strlen(line)-1]='\0';
        strcpy(lines[counter++], line);
    }
    while (counter < SIZE) {
        strcpy(lines[counter++], "\0");
    }

    for (int i = 0; i < SIZE; i++) {
        if (!strcmp(lines[i], "quit")) break;
        if (!strcmp(lines[i], "\0") && !strcmp(lines[i + 1], "\0")) break;

        else if (!strcmp(lines[i], "ident")) ident(transform);

        else if (!strcmp(lines[i], "apply")) matrix_mult(transform, edges);

        else if (!strcmp(lines[i], "line")) {
            char * args = lines[++i];
            double x0, y0, z0, x1, y1, z1;

            sscanf(args, "%le %le %le %le %le %le",
                   &x0, &y0, &z0, &x1, &y1, &z1);

            add_edge(edges, x0, y0, z0, x1, y1, z1);
        }

        else if (!strcmp(lines[i], "scale")) {
            char * args = lines[++i];
            double x, y, z;

            sscanf(args, "%le %le %le",
                   &x, &y, &z);

            struct matrix * scale = make_scale(x, y, z);
            matrix_mult(scale, transform);
        }

        else if (!strcmp(lines[i], "move")) {
            char * args = lines[++i];
            double x, y, z;

            sscanf(args, "%le %le %le",
                   &x, &y, &z);

            struct matrix * translate = make_translate(x, y, z);
            matrix_mult(translate, transform);
        }

        else if (!strcmp(lines[i], "rotate")) {
            char * args = lines[++i];

            double angle;
            char axis;
            sscanf(args, "%c %lf", &axis, &angle);
            double rad = angle * M_PI / 180;

            struct matrix * rotate;
            if ('x' == axis) rotate = make_rotX(rad);
            else if ('y' == axis) rotate = make_rotY(rad);
            else if ('z' == axis) rotate = make_rotZ(rad);

            matrix_mult(rotate, transform);
        }

        else if (!strcmp(lines[i], "circle")) {
            char * args = lines[++i];
            double cx, cy, cz, r;

            sscanf(args, "%le %le %le %le", &cx, &cy, &cz, &r);

            add_circle(edges, cx, cy, cz, r, line_step);
        }

        else if (!strcmp(lines[i], "hermite") || !strcmp(lines[i], "bezier")) {
            char * args = lines[++i];
            double x0, y0, x1, y1, x2, y2, x3, y3;
            int type;

            if (!strcmp(lines[i - 1], "hermite")) type = 0;
            else type = 1;

            sscanf(args, "%le %le %le %le %le %le %le %le",
                          &x0, &y0, &x1, &y1, &x2, &y2, &x3, &y3);

            add_curve(edges, x0, y0, x1, y1, x2, y2, x3, y3, line_step, type);
        }

        else if (!strcmp(lines[i], "sphere")) {
            char * args = lines[++i];
            double cx, cy, cz, r;

            sscanf(args, "%le %le %le %le", &cx, &cy, &cz, &r);

            add_sphere(edges, cx, cy, cz, r, sphere_step);
        }

        else if (!strcmp(lines[i], "torus")) {
            char * args = lines[++i];
            double cx, cy, cz, r1, r2;

            sscanf(args, "%le %le %le %le %le", &cx, &cy, &cz, &r1, &r2);

            add_torus(edges, cx, cy, cz, r1, r2, sphere_step);
        }

        else if (!strcmp(lines[i], "box")) {
            char * args = lines[++i];
            double x, y, z, width, height, depth;

            sscanf(args, "%le %le %le %le %le %le",
                          &x, &y, &z, &width, &height, &depth);

            add_box(edges, x, y, z, width, height, depth);
        }

        else if (!strcmp(lines[i], "clear")) {
            edges -> lastcol = 0;
        }

        else if (!strcmp(lines[i], "display")) {
            clear_screen(s);
            draw_lines(edges, s, c);
            display(s);
        }

        else if (!strcmp(lines[i], "save")) {
            clear_screen(s);
            draw_lines(edges, s, c);
            char * arg = lines[++i];
            save_extension(s, arg);
            printf("Saved as %s\n", arg);
        }

        else if (!strcmp(lines[i], "color")) {
            char * args = lines[++i];
            int r, g, b;

            sscanf(args, "%d %d %d",
                   &r, &g, &b);

            change_color(&c, r, g, b);
        }
    }
}