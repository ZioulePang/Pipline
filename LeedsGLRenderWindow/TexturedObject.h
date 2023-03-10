///////////////////////////////////////////////////
//
//  Hamish Carr
//  September, 2020
//
//  ------------------------
//  TexturedObject.h
//  ------------------------
//  
//  Base code for rendering assignments.
//
//  Minimalist (non-optimised) code for reading and 
//  rendering an object file
//  
//  We will make some hard assumptions about input file
//  quality. We will not check for manifoldness or 
//  normal direction, &c.  And if it doesn't work on 
//  all object files, that's fine.
//
//  While I could set it up to use QImage for textures,
//  I want this code to be reusable without Qt, so I 
//  shall make a hard assumption that textures are in 
//  ASCII PPM and use my own code to read them
//  
///////////////////////////////////////////////////

// include guard for TexturedObject
#ifndef TEXTURED_OBJECT_H
#define TEXTURED_OBJECT_H

// include the C++ standard libraries we need for the header
#include <vector>
#include <iostream>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

// and also include our fake GL library
#include "LeedsGL.h"

// include the unit with Cartesian 3-vectors
#include "Cartesian3.h"
// the render parameters
#include "RenderParameters.h"
// the image class for a texture
#include "RGBAImage.h" 

class TexturedObject
    { // class TexturedObject
    public:
    // vector of vertices
    std::vector<Cartesian3> vertices;

    // vector of normals
    std::vector<Cartesian3> normals;
    
    // vector of texture coordinates (stored as triple to simplify code)
    std::vector<Cartesian3> textureCoords;

    // vector of faces
    std::vector<std::vector<unsigned long long> > faceVertices;

    // corresponding vector of normals
    std::vector<std::vector<unsigned long long> > faceNormals;
    
    // corresponding vector of texture coordinates
    std::vector<std::vector<unsigned long long> > faceTexCoords;


    //Buffers for LeedsGL
    std::vector<Homogeneous4> vertexBuffer;
    std::vector<Homogeneous4> normalsBuffer;
    std::vector<Cartesian3> texCoordsBuffer;
    std::vector<RGBAValueF> colorsBuffer;

    // RGBA Image for storing a texture
    RGBAImage texture;

    // a variable to store the texture's ID on the GPU
    GLuint textureID;

    // centre of gravity - computed after reading
    Cartesian3 centreOfGravity;

    // size of object - i.e. radius of circumscribing sphere centred at centre of gravity
    float objectSize;

    // constructor will initialise to safe values
    TexturedObject();
    
    // read routine returns true on success, failure otherwise
    bool ReadObjectStream(std::istream &geometryStream, std::istream &textureStream);

    // write routine
    void WriteObjectStream(std::ostream &geometryStream, std::ostream &textureStream);


    //routine for students to use when rendering
    void LeedsGLRender(LeedsGL& leedsGL);

    // routine to transfer assets to Leeds GL
    void TransferAssetsToLeedsGL(LeedsGL& leedsGL);

    void InitializeLeedsGLArrays();

    }; // class TexturedObject

// end of include guard for TexturedObject
#endif
