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
#include "thread.h"

struct CPUBitmap {
    unsigned char *pixels;
    int  x, y;
    void *dataBlock;
    void (*bitmapExit)(void*);
    int windows[2];

    CPUBitmap( int width, int height, void *d = NULL ) {
        pixels = new unsigned char[width * height * 4];
        x = width;
        y = height;
        dataBlock = d;
    }

    ~CPUBitmap() {
        delete [] pixels;
    }

    unsigned char* get_ptr( void ) const   { return pixels; }
    long image_size( void ) const { return x * y * 4; }

    static void init(){
        int c=1;
        char dummy[] = "";
        char* d = dummy;
        glutInit( &c, &d);
        glutInitDisplayMode( GLUT_SINGLE | GLUT_RGBA );
    }

    void display_and_exit( void(*e)(void*) = NULL ) {
        CPUBitmap**   bitmap = get_bitmap_ptr();
        *bitmap = this;
        bitmapExit = e;
        // a bug in the Windows GLUT implementation prevents us from
        // passing zero arguments to glutInit()
        glutInitWindowSize( x, y );
        windows[0]=glutCreateWindow( "bitmap" );
        glutSetWindow(windows[0]);
        glutKeyboardFunc(Key);
        glutDisplayFunc(Draw);
        glutMainLoop();
        printf("after main loop\n");
    }

    // static method used for glut callbacks
    static CPUBitmap** get_bitmap_ptr( void ) {
        static CPUBitmap   *gBitmap;
        return &gBitmap;
    }

   // static method used for glut callbacks
    static void Key(unsigned char key, int x, int y) {
#if DEBUG
        printf("key %c pressed\n",key);
#endif
        switch (key) {
            case 'q':
                CPUBitmap*   bitmap = *(get_bitmap_ptr());
                if (bitmap->dataBlock != NULL && bitmap->bitmapExit != NULL){
                    bitmap->bitmapExit( bitmap->dataBlock );
                }
                glutDestroyWindow(bitmap->windows[0]);
                exit_thread();
                //exit(0);
        }
         //   pthread_cancel
    }

    // static method used for glut callbacks
    static void Draw( void ) {
        CPUBitmap*   bitmap = *(get_bitmap_ptr());
#if DEBUG
        printf("Draw: bitmap address: %p\n",bitmap);
        printf("Draw: data address: %p\n",bitmap->pixels);
#endif
        glutSetWindow(bitmap->windows[0]);
        glClearColor( 0.0, 0.0, 0.0, 1.0 );
        glClear( GL_COLOR_BUFFER_BIT );
        glDrawPixels( bitmap->x, bitmap->y, GL_RGBA, GL_UNSIGNED_BYTE, bitmap->pixels );
#if DEBUG
        FILE* f = fopen("draw.data","wb");
        fwrite(bitmap->pixels, sizeof(unsigned char), bitmap->x*bitmap->y*4, f);
        fclose(f);
#endif
        //glutSolidSphere(100,100,100);
        glFlush();
#if DEBUG
        printf("Draw: flushed\n");
#endif
    }
};

CUT_THREADPROC show_img_core(void* bitmap){
    CPUBitmap* img = ((CPUBitmap*)bitmap);
    *(img->get_bitmap_ptr())=img;
    glutInitWindowSize( img->x, img->y );
    img->windows[0]=glutCreateWindow( "bitmap" );
    glutSetWindow(img->windows[0]);
    glutKeyboardFunc(img->Key);
    glutDisplayFunc(img->Draw);
    glutMainLoop();
    printf("after main loop\n");
    CUT_THREADEND;
}

void show_img(ft* img, int w, int h, ft min, ft max){
    CPUBitmap image( w, h );
    unsigned char *ptr = image.get_ptr(), value;
#if DEBUG
    printf("show_img: width=%d, height=%d, min=%f, max=%f\n",w,h,min,max);
    printf("show_img: image address: %p\n", &image);
    printf("show_img: data address: %p\n", ptr);
#endif
    int offset;
    for(int i=0; i < w; i++) for(int j=0; j < h; j++){
        offset = j*w+i;
        value = img[offset]<min ? 0 : 
                img[offset]>max ? 255 : 
                (unsigned char)(255 * (img[offset]-min)/(max-min));
        offset = (h-1-j)*w+i;
        ptr[(offset<<2)+0] = value;
        ptr[(offset<<2)+1] = value;
        ptr[(offset<<2)+2] = value;
        ptr[(offset<<2)+3] = 0xff;
    }
    CUTThread thread;

    thread=start_thread(show_img_core, &image);
    end_thread(thread);
}

void show_img(ft* img, int w, int h, ft min){
    ft max=img[0];
    for(int i=0; i < w*h; i++){
        max=MAX(img[i],max);
    }
    show_img(img,w,h,min,max);
}

void show_img(ft* img, int w, int h){
    ft min=img[0], max=img[0];
    for(int i=0; i < w*h; i++){
        min=MIN(img[i],min);
        max=MAX(img[i],max);
    }
    show_img(img,w,h,min,max);
}


#endif  // __CPU_BITMAP_H__

