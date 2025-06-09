/*
* Autor: Wilker Sebastian Afonso Pereira
* GitHub: https://github.com/WilkerSebastian/p1-tp2
*
* Prova 1 de Tópicos em Computação 2
*/
#ifndef __DefaultSoAAllocator_h
#define __DefaultSoAAllocator_h

#include <cstddef>

struct DefaultSoAAllocator 
{

    template <typename T>
    static T* allocate(size_t count) 
    {

        if (count == 0) 
            return nullptr;
            
        return new T[count];
    }

    template <typename T>
    static void free(T* ptr) 
    {
        delete[] ptr;
    }

};

#endif