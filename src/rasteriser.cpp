#include "rasteriser.hpp"
#include "utilities/lodepng.h"
#include <vector>
#include <chrono>
#include <math.h>

// --- Overview ---

// I'm going to assume most of you who are reading through this file have never worked with computer graphics before.
// Let me therefore start by explaining the general idea of this file a little bit, because you'll need to at least _slightly_
// understand what is going on in this program if you are to optimise it effectively.

// We start off with two lists; one with vertices, and one with indices.
// Vertices are three-dimensional coordinates (even though in this file they have the type float4, where the last coordinate is always 1. The reason why is a long story).
// Indices represent the indices of vertices in the vertex buffer that should be connected together,
// much like you do in "connect the dots". Except here every three indices define a triangle.
// So for instance, index 0, 1, and 2 reference the three vertices of the first triangle. Indices 3, 4, and 5 the next. And so on.

// Rendering an image is done in two steps.
// First, all vertices are multiplied with a so-called "transformation matrix" in a section / program called the "vertex shader".
// How it works is not really relevant now. All you should know is that its purpose is to move objects around on the scene, and apply a "perspective" effect
// that causes objects to appear smaller the further they move away from the camera.
// The entire process only requires that all vertices are multiplied by a 4x4 matrix, using regular matrix multiplication (4x4 matrix with a 4x1 vertex)
// The second main stage involves rasterising the transformed triangles
// A "depth buffer" is used to make sure that triangles that should appear in front of others
// A "fragment shader" is run on each rendered pixel, which determines the colour the pixel should receive.
// The rendered image is stored in the so-called "framebuffer".


/**
 * Executes the vertex shader, transforms vertices and normals of the mesh object
 * @param mesh                    Mesh object with all vertices and normals
 * @param transformedVertexBuffer returned transformed vertices
 * @param transformedNormalBuffer returned transformed normals
 */
void runVertexShader( Mesh &mesh,
					  std::vector<float4> &transformedVertexBuffer,
					  std::vector<float4> &transformedNormalBuffer )
{
	// The & in front of the variable names cause the function to modify variables from the function
	// calling this one, rather than making a copy of them.

	// The matrices defined below are the ones used to transform the vertices and normals.

	// This projection matrix assumes a 16:9 aspect ratio, and an field of view (FOV) of 90 degrees.
	mat4x4 projectionMatrix(
		0.347270,   0, 			0, 		0,
		0,	  		0.617370, 	0,		0,
		0,	  		0,			-1, 	-0.2f,
		0,	  		0,			-1,		0);

	mat4x4 viewMatrix(
		0.5, 	0, 		0, 		5,
		0, 		-0.5, 	0, 		30,
		0, 		0, 		0.5, 	-55,
		0, 		0, 		0, 		1);

	mat4x4 normalMatrix(
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1);

	mat4x4 MVP = projectionMatrix * viewMatrix;

	for(unsigned int i = 0; i < transformedVertexBuffer.size(); i++) {
		float4 transformed = MVP * mesh.vertices[i];
		transformed.x /= transformed.w;
		transformed.y /= transformed.w;
		transformed.z /= transformed.w;
		transformedVertexBuffer.at(i) = transformed;
	}

	for(unsigned int j = 0; j < transformedNormalBuffer.size(); j++) {
		float4 normal;
		normal.x = mesh.normals[j].x;
		normal.y = mesh.normals[j].y;
		normal.z = mesh.normals[j].z;
		normal.w = 1;
		transformedNormalBuffer.at(j) = normalMatrix * normal;
	}
}

/**
 * Executes the fragment shader to calculate the colour of the pixel from its
 * normal
 * @param normal triangle pixel normal
 * @return		 colour of the pixel in RGBA
 */
inline void runFragmentShader(std::vector<unsigned char> &pixelColour, float3 const normal )
{
	//std::vector<unsigned char> pixelColour(4);
	const float3 lightDirection(0.0f, 0.0f, 1.0f);

	// Computing the dot product between the surface normal and a light
	// direction gives a disffuse-like reflection. It looks more than
	// good enough for a few static images.
	float colour = normal.x * lightDirection.x +
		normal.y * lightDirection.y +
		normal.z * lightDirection.z;

	// We first scale the colour value from a range between 0 and 1,
	// to between 0 and 255.
	// Since single bytes are only able to go between 0 and 255,
	// we subsequently clamp the colour to lie within that range.
	unsigned char colourByte = (unsigned char) std::min(255.0f,	std::max(colour * 255.0f, 0.0f));

	// And this writes the pixel to the pixelColor vector. The first three
	// channels are red, green, and blue. The fourth represents transparency.
	pixelColour[0] = colourByte;
	pixelColour[1] = colourByte;
	pixelColour[2] = colourByte;
	//pixelColour[3] = 255;

	// This colour vector is supposed to go into the frame buffer
	//return pixelColour;
}

/**
 * interpolates given normals using the barycentric weights
 *
 * Using float3 to return all three normals at once
 *
 * @param  n0 normal of triangle vertex
 * @param  n1 normal of triangle vertex
 * @param  n2 normal of triangle vertex
 * @param  w0 barycentric weight
 * @param  w1 barycentric weight
 * @param  w2 barycentric weight
 * @return    interpolated normal
 */
inline void interpolateNormals(float3 * res, float4 const n0,
						   float4 const n1,
					 	   float4 const n2,
					 	   float const w0,
						   float const w1,
						   float const w2 )
{

	res->x = w0 * n0.x + w1 * n1.x + w2 * n2.x;
	res->y = w0 * n0.y + w1 * n1.y + w2 * n2.y;
	res->z = w0 * n0.z + w1 * n1.z + w2 * n2.z;

}

/**
 * converts a vertex from clipping space to screen pixel coordinates
 * @param  vertex vertex
 * @param  width  screen width
 * @param  height screen height
 * @return        vertex in screen pixel coordinates
 */
inline float4 convertClippingSpace( float4 const vertex,
							 unsigned int const width,
							 unsigned int const height )
{
	float4 res;
	res.x = (vertex.x + 0.5f) * float(width);
	res.y = (vertex.y + 0.5f) * float(height);
	res.z = vertex.z;
	res.w = vertex.w;
	return res;
}

/**
 * Computing barycentric weights. This is both for determining whether the point
 * lies within the triangle, as well as interpolating coordinates needed for
 * rendering. If you need an explanation, you can find one here:
 * https://codeplea.com/triangular-interpolation
 * Using float3 to return all weights at once (x -> w0, y -> w1, z -> w2)
 *
 * @param  v0 triangle vertex
 * @param  v1 triangle vertex
 * @param  v2 triangle vertex
 * @param  x  screen pixel x-coordinate
 * @param  y  screen pixel y-coordinate
 * @return    barycentric weights of the pixel in relation to the triangle vertices
 */
inline void getTriangleBarycentricWeights(float3 &res, float4 const v0,
									  float4 const v1,
									  float4 const v2,
									  unsigned int const x,
									  unsigned int const y )
{
	float v2x = v2.x;
	float v2y = v2.y;
	float div = 1/(((v1.y - v2y) * (v0.x - v2x)) + ((v2x - v1.x) * (v0.y - v2y)));

	res.x = (((v1.y - v2y) * (x    - v2x)) + ((v2x - v1.x) * (y    - v2y))) * div;
	res.y = (((v2y - v0.y) * (x    - v2x)) + ((v0.x - v2x) * (y    - v2y))) * div;
	res.z = 1 - res.x - res.y;

}

/**
 * Calculates the triangle's pixel depth from the barycentric weights
 * @param  v0 triangle vertex
 * @param  v1 triangle vertex
 * @param  v2 triangle vertex
 * @param  w0 barycentric weight
 * @param  w1 barycentric weight
 * @param  w2 barycentric weight
 * @return    pixel depth (z)
 */
inline void getTrianglePixelDepth(float &pixelDepth, float4 const v0,
							 float4 const v1,
							 float4 const v2,
							 float const w0,
							 float const w1,
							 float const w2 )
{
	pixelDepth = w0 * v0.z + w1 * v1.z + w2 * v2.z;
}

/**
 * The main procedure which rasterises all triangles on the framebuffer
 * @param mesh                    Mesh object
 * @param transformedVertexBuffer transformed vertices from the mesh obj
 * @param transformedNormalBuffer transformed normals from the mesh obj
 * @param frameBuffer             frame buffer for the rendered image
 * @param depthBuffer             depth buffer for every pixel on the image
 * @param width                   width of the image
 * @param height                  height of the image
 */
void rasteriseTriangles( Mesh &mesh,
                         std::vector<float4> &transformedVertexBuffer,
                         std::vector<float4> &transformedNormalBuffer,
                         std::vector<unsigned char> &frameBuffer,
                         std::vector<float> &depthBuffer,
                         unsigned int width,
                         unsigned int height )
{
	// We rasterise one triangle at a time
	unsigned int triangleCount = mesh.indexCount / 3;



	// We look up those triangles here
	float4 *vertex0 = new float4();
	float4 *vertex1 = new float4();
	float4 *vertex2 = new float4();

	float4 *normal0 = new float4();
	float4 *normal1 = new float4();
	float4 *normal2 = new float4();

	float3 *interpolatedNormal = new float3();

	float weight0;
	float weight1;
	float weight2;
	float pixelDepth;
	float normalLength;
	float invNormalLength;

	float3 barycentricWeights;
	std::vector<unsigned char> pixelColour(4);

	for(unsigned int triangleIndex = 0; triangleIndex < triangleCount; triangleIndex++) {
		// '\r' returns to the beginning of the current line
		std::cout << "Rasterising triangle " << (triangleIndex + 1) << "/" << triangleCount << "\r" << std::flush;

		// As vertices are commonly reused within a model, rendering libraries use an
		// index buffer which specifies the indices of the vertices in the vertex buffer
		// which together make up the specific triangle.
		unsigned int index0 = mesh.indices[3 * triangleIndex + 0];
		unsigned int index1 = mesh.indices[3 * triangleIndex + 1];
		unsigned int index2 = mesh.indices[3 * triangleIndex + 2];

		// We look up those triangles here
		*vertex0 = transformedVertexBuffer[index0];
		*vertex1 = transformedVertexBuffer[index1];
		*vertex2 = transformedVertexBuffer[index2];

		// These triangles are still in so-called "clipping space". We first convert them
		// to screen pixel coordinates
		*vertex0 = convertClippingSpace(*vertex0, width, height);
		*vertex1 = convertClippingSpace(*vertex1, width, height);
		*vertex2 = convertClippingSpace(*vertex2, width, height);

		// We have to find the max and min for x and y for the triangle so we iterate on a smaller number of pixel
		float xMax = vertex0->x;
		if (vertex1->x > xMax) {
			xMax = vertex1->x;
		}
		if (vertex2->x > xMax) {
			xMax = vertex2->x;
		}
		float xMin = vertex0->x;
		if (vertex1->x < xMin) {
			xMin = vertex1->x;
		}
		if (vertex2->x < xMin) {
			xMin = vertex2->x;
		}
		unsigned int xMaxInt = floor(xMax)+1;
		unsigned int xMinInt = floor(xMin);

		float yMax = vertex0->y;
		if (vertex1->y > yMax) {
			yMax = vertex1->y;
		}
		if (vertex2->y > yMax) {
			yMax = vertex2->y;
		}
		float yMin = vertex0->y;
		if (vertex1->y < yMin) {
			yMin = vertex1->y;
		}
		if (vertex2->y < yMin) {
			yMin = vertex2->y;
		}
		unsigned int yMaxInt = floor(yMax)+1;
		unsigned int yMinInt = floor(yMin);

		// Read the normals belonging to each vertex
		*normal0 = transformedNormalBuffer[index0];
		*normal1 = transformedNormalBuffer[index1];
		*normal2 = transformedNormalBuffer[index2];
        // We now iterate only on a small square delimited by max and min variables
        // calculated before. The triangle we are computing is inside this square
		for(unsigned int x = xMinInt; x < xMaxInt; x++) {
			for(unsigned int y = yMinInt; y < yMaxInt; y++) {
				//Coordinate of the current pixel in the framebuffer, remember RGBA color code
				unsigned int pixelBaseCoordinate = 4 * (x + y * width);

				// Calculating the barycentric weights of the pixel in relation to the triangle
				getTriangleBarycentricWeights(barycentricWeights, *vertex0, *vertex1, *vertex2, x, y);
				weight0 = barycentricWeights.x;
				weight1 = barycentricWeights.y;
				weight2 = barycentricWeights.z;

				// Now we can determine the depth of our pixel
				getTrianglePixelDepth(pixelDepth, *vertex0, *vertex1, *vertex2, weight0, weight1, weight2);

				// Z-clipping discards pixels too close or too far from the camera
				if(pixelDepth >= -1 && pixelDepth <= 1) {
					// The weights have the nice property that if only one is negative, the pixel lies outside the triangle
					if(weight0 >= 0 && weight1 >= 0 && weight2 >= 0) {
						//Have we drawn a pixel above the current?
						if(pixelDepth < depthBuffer[y * width + x]) {
							// This pixel is going into the frame buffer,
							// save its depth to skip all next pixels underneath it
							depthBuffer[y * width + x] = pixelDepth;
							// Copy the calculated pixel colour into the frame buffer - RGBA
							/*for (unsigned int i = 0; i < pixelColour.size(); i++) {
									frameBuffer.at(pixelBaseCoordinate + i) = pixelColour.at(i);
							}*/
							interpolateNormals(interpolatedNormal,*normal0, *normal1, *normal2, weight0, weight1, weight2);

							normalLength = std::sqrt( interpolatedNormal->x * interpolatedNormal->x +
								interpolatedNormal->y * interpolatedNormal->y +
								interpolatedNormal->z * interpolatedNormal->z );
							invNormalLength = 1/normalLength;
							interpolatedNormal->x *= invNormalLength;
							interpolatedNormal->y *= invNormalLength;
							interpolatedNormal->z *= invNormalLength;
							runFragmentShader(pixelColour,*interpolatedNormal);
							frameBuffer[pixelBaseCoordinate + 0] = pixelColour[0];
							frameBuffer[pixelBaseCoordinate + 1] = pixelColour[1];
							frameBuffer[pixelBaseCoordinate + 2] = pixelColour[2];
							//frameBuffer[pixelBaseCoordinate + 3] = pixelColour[3];
						}
					}
				}
			}
		}
		// Cleanup
	}
	delete normal0;
	delete normal1;
	delete normal2;
	delete vertex0;
	delete vertex1;
	delete vertex2;
	delete interpolatedNormal;
	// finish the progress output with a new line
	std::cout << std::endl;
}

/**
 * Procedure to kick of the rasterisation process
 * @param mesh            Mesh object
 * @param outputImageFile path of the output image
 * @param width           width of the output image
 * @param height          height of the output image
 */
void rasterise(Mesh mesh, std::string outputImageFile, unsigned int width, unsigned int height) {
	// We first need to allocate some buffers.

	// The framebuffer contains the image being rendered.
	std::vector<unsigned char> frameBuffer;
	frameBuffer.resize(width * height * 4, 0);

	// The depth buffer is used to make sure that objects closer to the camera occlude/obscure objects that are behind it
	std::vector<float> depthBuffer;
	depthBuffer.resize(width * height, 1);

	// And these two buffers store vertices and normals processed by the vertex shader.
	std::vector<float4> transformedVertexBuffer;
	transformedVertexBuffer.resize(mesh.vertexCount);

	std::vector<float4> transformedNormalBuffer;
	transformedNormalBuffer.resize(mesh.vertexCount);

	for(unsigned int x = 0; x < width; x++) {
		for (unsigned int y = 0; y < height; y++) {
			frameBuffer.at(4 * ( x + y * width ) + 0) = 0;
			frameBuffer.at(4 * ( x + y * width ) + 1) = 0;
			frameBuffer.at(4 * ( x + y * width ) + 2) = 0;
			frameBuffer.at(4 * ( x + y * width ) + 3) = 255;
		}
	}

	std::cout << "Running the vertex shader... ";

	// Internal profiling
	auto start = std::chrono::high_resolution_clock::now();
	runVertexShader(mesh, transformedVertexBuffer, transformedNormalBuffer);
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "Execution of runVertexShader took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;

	start = std::chrono::high_resolution_clock::now();
	rasteriseTriangles(mesh, transformedVertexBuffer, transformedNormalBuffer, frameBuffer, depthBuffer, width, height);
	end = std::chrono::high_resolution_clock::now();
	std::cout << "Execution of rasteriseTriangles took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;

	std::cout << "complete!" << std::endl;


	std::cout << "Finished rendering!" << std::endl;

	std::cout << "Writing image to '" << outputImageFile << "'..." << std::endl;

	unsigned error = lodepng::encode(outputImageFile, frameBuffer, width, height);

	if(error)
	{
		std::cout << "An error occurred while writing the image file: " << error << ": " << lodepng_error_text(error) << std::endl;
	}
}
