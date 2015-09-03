#ifndef __THREAD_H_
#define __THREAD_H_

#if _WIN32
    //Windows threads.
    #include <windows.h>
    #include <tchar.h>
    #include <strsafe.h>


    typedef HANDLE CUTThread;
    typedef unsigned (WINAPI *CUT_THREADROUTINE)(void *);

    #define CUT_THREADPROC unsigned WINAPI
    #define CUT_THREADEND return 0

#else
    //POSIX threads.
    #include <pthread.h>
    #include <stdlib.h>
    #include <stdio.h>

    typedef pthread_t CUTThread;
    typedef void* (*CUT_THREADROUTINE)(void *);

    #define CUT_THREADPROC void*
    #define CUT_THREADEND return 0
#endif

//Create thread.
CUTThread start_thread( CUT_THREADROUTINE, void *data );

//Wait for thread to finish.
void end_thread( CUTThread thread );

//Destroy thread.
void destroy_thread( CUTThread thread );

//Wait for multiple threads.
void wait_for_threads( const CUTThread *threads, int num );

// kill the calling thread
void exit_thread();

#if _WIN32
    //Create thread
    CUTThread start_thread(CUT_THREADROUTINE func, void *data){
        HANDLE thread=CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)func, data, 0, NULL);
        if(thread==NULL){
            printf(TEXT("CreateThread"));
            ExitProcess(3);
        }
        return thread;
    }

    //Wait for thread to finish
    void end_thread(CUTThread thread){
        WaitForSingleObject(thread, INFINITE);
        CloseHandle(thread);
    }

    //Destroy thread
    void destroy_thread( CUTThread thread ){
        TerminateThread(thread, 0);
        CloseHandle(thread);
    }

    //Wait for multiple threads
    void wait_for_threads(const CUTThread * threads, int num){
        WaitForMultipleObjects(num, threads, TRUE, INFINITE);

        for(int i = 0; i < num; i++)
            CloseHandle(threads[i]);
    }

    void exit_thread(){
        ExitThread(0);
    }

#else
    //Create thread
    CUTThread start_thread(CUT_THREADROUTINE func, void * data){
        pthread_t thread;
        int res=pthread_create(&thread, NULL, func, data);
        if (res != 0) {
            perror("Thread creation failed");
            exit(EXIT_FAILURE);
        }
        return thread;
    }

    //Wait for thread to finish
    void end_thread(CUTThread thread){
#if DEBUG
        void *thread_result;
        int res=pthread_join(thread, &thread_result);
#else
        int res=pthread_join(thread, NULL);
#endif
        if (res != 0) {
            perror("Thread join failed.");
            exit(EXIT_FAILURE);
        }
#if DEBUG
        printf("Thread joined, it returned %s\n", (char *)thread_result);
#endif
    }

    //Destroy thread
    void destroy_thread( CUTThread thread ){
        pthread_cancel(thread);
    }

    //Wait for multiple threads
    void wait_for_threads(const CUTThread * threads, int num){
        for(int i = 0; i < num; i++)
            end_thread( threads[i] );
    }

    void exit_thread(){
        int a;
        pthread_exit(&a);
    }

#endif

#endif

