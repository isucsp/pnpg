#ifndef __THREAD_H_
#define __THREAD_H_

#if _WIN32
    //Windows threads.
    #include <windows.h>

    typedef HANDLE CUTThread;
    typedef unsigned (WINAPI *CUT_THREADROUTINE)(void *);

    #define CUT_THREADPROC unsigned WINAPI
    #define  CUT_THREADEND return 0

#else
    //POSIX threads.
    #include <pthread.h>

    typedef pthread_t CUTThread;
    typedef void *(*CUT_THREADROUTINE)(void *);

    #define CUT_THREADPROC void
    #define  CUT_THREADEND
#endif

//Create thread.
CUTThread start_thread( CUT_THREADROUTINE, void *data );

//Wait for thread to finish.
void end_thread( CUTThread thread );

//Destroy thread.
void destroy_thread( CUTThread thread );

//Wait for multiple threads.
void wait_for_threads( const CUTThread *threads, int num );

#if _WIN32
    //Create thread
    CUTThread start_thread(CUT_THREADROUTINE func, void *data){
        return CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)func, data, 0, NULL);
    }

    void start_thread(CUTThread* thread, CUT_THREADROUTINE func, void *data){
        &thread=CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)func, data, 0, NULL);
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
        WaitForMultipleObjects(num, threads, true, INFINITE);

        for(int i = 0; i < num; i++)
            CloseHandle(threads[i]);
    }

#else
    //Create thread
    CUTThread start_thread(CUT_THREADROUTINE func, void * data){
        pthread_t thread;
        pthread_create(&thread, NULL, func, data);
        if (res != 0) {
            perror("Thread creation failed");
            exit(EXIT_FAILURE);
        }
        return thread;
    }

    void start_thread(CUTThread* thread, CUT_THREADROUTINE func, void * data){
        pthread_create(thread, NULL, func, data);
        if (res != 0) {
            perror("Thread creation failed");
            exit(EXIT_FAILURE);
        }
    }

    //Wait for thread to finish
    void end_thread(CUTThread thread){
        pthread_join(thread, NULL);
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

#endif

#endif

