/*
 * Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
 *
 * NVIDIA Corporation and its licensors retain all intellectual property and
 * proprietary rights in and to this software and related documentation.
 * Any use, reproduction, disclosure, or distribution of this software
 * and related documentation without an express license agreement from
 * NVIDIA Corporation is strictly prohibited.
 *
 * Please refer to the applicable NVIDIA end user license agreement (EULA)
 * associated with this source code for terms and conditions that govern
 * your use of this NVIDIA software.
 *
 */


#ifndef __CPU_BITMAP_H__
#define __CPU_BITMAP_H__

#include "gl_helper.h"
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glext.h>
#include <GL/glx.h>

struct CPUBitmap {
    unsigned char    *pixels;
    int     x, y;
    void    *dataBlock;
    void (*bitmapExit)(void*);
}

void CPUBitmapCon(CPUBitmap* bm, int width, int height, void *d = NULL ) {
    bm->pixels = new unsigned char[width * height * 4];
    bm->x = width;
    bm->y = height;
    bm->dataBlock = d;
}

long image_size( CPUBitmap* bm ) const { return bm->x * bm->y * 4; }

void display_and_exit( CPUBitmap *bm, void(*e)(void*) = NULL ) {
    bm->bitmapExit = e;
    // a bug in the Windows GLUT implementation prevents us from
    // passing zero arguments to glutInit()
    int c=1;
    char dummy[] = "";
    char* d = dummy;
    glutInit( &c, &d);
    glutInitDisplayMode( GLUT_SINGLE | GLUT_RGBA );
    glutInitWindowSize( bm->x, bm->y );
    glutCreateWindow( "bitmap" );
    glutKeyboardFunc(Key);
    glutDisplayFunc(Draw);
    glutMainLoop();
}

// static method used for glut callbacks
static void Key(CPUBitmap *bitmap, unsigned char key, int x, int y) {
    printf("key=%d\n",key);
    switch (key) {
        case 27:
            if (bitmap->dataBlock != NULL && bitmap->bitmapExit != NULL)
                bitmap->bitmapExit( bitmap->dataBlock );
            exit(0);
    }
}

// static method used for glut callbacks
static void Draw( void ) {
    CPUBitmap*   bitmap = *(get_bitmap_ptr());
    glClearColor( 0.0, 0.0, 0.0, 1.0 );
    glClear( GL_COLOR_BUFFER_BIT );
    glDrawPixels( bitmap->x, bitmap->y, GL_RGBA, GL_UNSIGNED_BYTE, bitmap->pixels );
    glFlush();
}

#endif  // __CPU_BITMAP_H__

