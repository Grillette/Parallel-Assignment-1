#include "sse_test.hpp"
#include <stdlib.h>
#include <time.h>

float randFloat() {
    return (float) (rand()) / (float) (RAND_MAX);
}

// We change the our rand variables datatype so we change aswell the datatype of the
// randFloat4 function return variable
sse_float4 randFloat4() {
    sse_float4 res;
    res.elements[0] = randFloat();
    res.elements[1] = randFloat();
    res.elements[2] = randFloat();
    res.elements[3] = randFloat();
    return res;
}

void sse_test(Mesh &mesh) {
    //Not allowed to change:
    unsigned int const loadFactor = 1000;

    // We can divide the iteration by 4 because we are gonna use SSE and do 4 operations
    // at the same time
    unsigned int sizeSSE = mesh.vertexCount/4;
    // Now our vertices and rand variables use sse_float4 type
    sse_float4* vertices = new sse_float4[sizeSSE];
    sse_float4* rand1 = new sse_float4[sizeSSE];
    sse_float4* rand2 = new sse_float4[sizeSSE];
    sse_float4* rand3 = new sse_float4[sizeSSE];
    sse_float4* rand4 = new sse_float4[sizeSSE];

    srand(time(NULL));

    std::cout << "SSE_TEST: Initializing vectors... " << std::flush;
    for (unsigned int i=0; i < sizeSSE; i++) {
        vertices[i].elements[0] = mesh.vertices[i].x;
        vertices[i].elements[1] = mesh.vertices[i].y;
        vertices[i].elements[2] = mesh.vertices[i].z;
        vertices[i].elements[3] = mesh.vertices[i].w;
        rand1[i] = randFloat4();
        rand2[i] = randFloat4();
        rand3[i] = randFloat4();
        rand4[i] = randFloat4();
    }
    std::cout << "finished!"  << std::endl;
    for (unsigned int j = 0; j < loadFactor; j++) {
        std::cout << "SSE_TEST: " << (j+1) << "/" << loadFactor << " Crunching numbers on " << sizeSSE << " vertices... " << "\r" << std::flush;
        for (unsigned int i=0; i < sizeSSE; i++) {
            vertices[i].vector = vertices[i].vector + rand1[i].vector;
            vertices[i].vector = vertices[i].vector - rand2[i].vector;
            vertices[i].vector = vertices[i].vector * rand3[i].vector;
            if (rand4[i].elements != 0) {
                vertices[i].vector = vertices[i].vector / rand4[i].vector;
            }
        }
    }

    std::cout << std::endl;
}
