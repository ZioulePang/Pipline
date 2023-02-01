///////////////////////////////////////////////////
//
//  Hamish Carr
//  September, 2020
//
//  ------------------------
//  TexturedObject.cpp
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

// include the header file
#include "TexturedObject.h"

// include the C++ standard libraries we want
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <memory>

// include the Cartesian 3- vector class
#include "Cartesian3.h"

#define MAXIMUM_LINE_LENGTH 1024

// constructor will initialise to safe values
TexturedObject::TexturedObject()
    : centreOfGravity(0.0,0.0,0.0)
    { // TexturedObject()
    // force arrays to size 0
    vertices.resize(0);
    normals.resize(0);
    textureCoords.resize(0);
    } // TexturedObject()

// read routine returns true on success, failure otherwise
bool TexturedObject::ReadObjectStream(std::istream &geometryStream, std::istream &textureStream)
    { // ReadObjectStream()
    
    // create a read buffer
    char readBuffer[MAXIMUM_LINE_LENGTH];
    
    // the rest of this is a loop reading lines & adding them in appropriate places
    while (true)
        { // not eof
        // character to read
        char firstChar = geometryStream.get();
        // check for eof() in case we've run out
        if (geometryStream.eof())
            break;

        // otherwise, switch on the character we read
        switch (firstChar)
            { // switch on first character
            case '#':       // comment line
                // read and discard the line
                geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
                break;
                
            case 'v':       // vertex data of some type
                { // some sort of vertex data
                // retrieve another character
                char secondChar = geometryStream.get();
                
                // bail if we ran out of file
                if (geometryStream.eof())
                    break;

                // now use the second character to choose branch
                switch (secondChar)
                    { // switch on second character
                    case ' ':       // space - indicates a vertex
                        { // vertex read
                        Cartesian3 vertex;
                        geometryStream >> vertex;
                        vertices.push_back(vertex);
                        break;
                        } // vertex read
                    case 'n':       // n indicates normal vector
                        { // normal read
                        Cartesian3 normal;
                        geometryStream >> normal;
                        normals.push_back(normal);
                        break;
                        } // normal read
                    case 't':       // t indicates texture coords
                        { // tex coord
                        Cartesian3 texCoord;
                        geometryStream >> texCoord;
                        textureCoords.push_back(texCoord);
                        break;                  
                        } // tex coord
                    default:
                        break;
                    } // switch on second character 
                break;
                } // some sort of vertex data
                
            case 'f':       // face data
                { // face
                // a face can have an arbitrary number of vertices, which is a pain
                // so we will create a separate buffer to read from
                geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);

                // turn into a C++ string
                std::string lineString = std::string(readBuffer);

                // create a string stream
                std::stringstream lineParse(lineString); 

                // create vectors for the IDs (with different names from the master arrays)
                std::vector<unsigned long long> faceVertexSet;
                std::vector<unsigned long long> faceNormalSet;
                std::vector<unsigned long long> faceTexCoordSet;
                
                // now loop through the line
                while (!lineParse.eof())
                    { // lineParse isn't done
                    // the triple of vertex, normal, tex coord IDs
                    unsigned int vertexID;
                    unsigned int normalID;
                    unsigned int texCoordID;

                    // try reading them in, breaking if we hit eof
                    lineParse >> vertexID;
                    // retrieve & discard a slash
                    lineParse.get();
                    // check eof
                    if (lineParse.eof())
                        break;
                    
                    // and the tex coord
                    lineParse >> texCoordID;
                    lineParse.get();
                    if (lineParse.eof())
                        break;
                        
                    // read normal likewise
                    lineParse >> normalID;
                        
                    // if we got this far, we presumably have three valid numbers, so add them
                    // but notice that .obj uses 1-based numbering, where our arrays use 0-based
                    faceVertexSet.push_back(vertexID-1);
                    faceNormalSet.push_back(normalID-1);
                    faceTexCoordSet.push_back(texCoordID-1);
                    } // lineParse isn't done

                // as long as the face has at least three vertices, add to the master list
                if (faceVertexSet.size() > 2)
                    { // at least 3
                    faceVertices.push_back(faceVertexSet);
                    faceNormals.push_back(faceNormalSet);
                    faceTexCoords.push_back(faceTexCoordSet);
                    } // at least 3
                
                break;
                } // face
                
            // default processing: do nothing
            default:
                break;

            } // switch on first character

        } // not eof

    // compute centre of gravity
    // note that very large files may have numerical problems with this
    centreOfGravity = Cartesian3(0.0, 0.0, 0.0);

    // if there are any vertices at all
    if (vertices.size() != 0)
        { // non-empty vertex set
        // sum up all of the vertex positions
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
            centreOfGravity = centreOfGravity + vertices[vertex];
        
        // and divide through by the number to get the average position
        // also known as the barycentre
        centreOfGravity = centreOfGravity / vertices.size();

        // start with 0 radius
        objectSize = 0.0;

        // now compute the largest distance from the origin to a vertex
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
            { // per vertex
            // compute the distance from the barycentre
            float distance = (vertices[vertex] - centreOfGravity).length();         
            
            // now test for maximality
            if (distance > objectSize)
                objectSize = distance;
            } // per vertex
        } // non-empty vertex set

    // now read in the texture file

    texture.ReadPPM(textureStream);

    // return a success code
    return true;
    } // ReadObjectStream()

// write routine
void TexturedObject::WriteObjectStream(std::ostream &geometryStream, std::ostream &textureStream)
    { // WriteObjectStream()
    // output the vertex coordinates
    for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
        geometryStream << "v  " << std::fixed << vertices[vertex] << std::endl;
    geometryStream << "# " << vertices.size() << " vertices" << std::endl;
    geometryStream << std::endl;

    // and the normal vectors
    for (unsigned int normal = 0; normal < normals.size(); normal++)
        geometryStream << "vn " << std::fixed << normals[normal] << std::endl;
    geometryStream << "# " << normals.size() << " vertex normals" << std::endl;
    geometryStream << std::endl;

    // and the texture coordinates
    for (unsigned int texCoord = 0; texCoord < textureCoords.size(); texCoord++)
        geometryStream << "vt " << std::fixed << textureCoords[texCoord] << std::endl;
    geometryStream << "# " << textureCoords.size() << " texture coords" << std::endl;
    geometryStream << std::endl;

    // and the faces
    for (unsigned int face = 0; face < faceVertices.size(); face++)
        { // per face
        geometryStream << "f ";
        
        // loop through # of vertices
        for (unsigned int vertex = 0; vertex < faceVertices[face].size(); vertex++)
            geometryStream << faceVertices[face][vertex]+1 << "/" << faceTexCoords[face][vertex]+1 << "/" << faceNormals[face][vertex]+1 << " " ;
        
        geometryStream << std::endl;
        } // per face
    geometryStream << "# " << faceVertices.size() << " polygons" << std::endl;
    geometryStream << std::endl;
    
    // now output the texture
    texture.WritePPM(textureStream);
    } // WriteObjectStream()


void TexturedObject::TransferAssetsToLeedsGL(LeedsGL& leedsGL)
{
    leedsGL.texImage2D(&texture);
}

void TexturedObject::InitializeLeedsGLArrays()
{

    // loop through the faces: note that they may not be triangles, which complicates life
    for (unsigned long long face = 0; face < faceVertices.size(); face++)
        { // per face
        // on each face, treat it as a triangle fan starting with the first vertex on the face
        for (unsigned long long triangle = 0; triangle < faceVertices[face].size() - 2; triangle++)
            { // per triangle
            // now do a loop over three vertices
            for (unsigned long long vertex = 0; vertex < 3; vertex++)
                { // per vertex
                // we always use the face's vertex 0
                unsigned long long faceVertex = 0;
                // so if it isn't 0, we want to add the triangle base ID
                if (vertex != 0)
                    faceVertex = triangle + vertex;

                // now we use that ID to lookup
                Cartesian3 n
                    (
                    normals         [faceNormals    [face][faceVertex]  ].x,
                    normals         [faceNormals    [face][faceVertex]  ].y,
                    normals         [faceNormals    [face][faceVertex]  ].z
                    );

                // set the texture coordinate
                Cartesian3 uv
                    (
                    textureCoords   [faceTexCoords  [face][faceVertex]  ].x,
                    textureCoords   [faceTexCoords  [face][faceVertex]  ].y,
                    0
                    );

                //we set the color of the model to the UV coordinates, for debug purposes.
                RGBAValueF c = RGBAValueF(uv.x,uv.y,0);

                // and set the vertex position
                Homogeneous4 p
                    (
                    vertices        [faceVertices   [face][faceVertex]].x,
                    vertices        [faceVertices   [face][faceVertex]].y,
                    vertices        [faceVertices   [face][faceVertex]].z
                    );

                //add them to the buffer we want to render
                vertexBuffer.push_back(p);
                normalsBuffer.push_back(n);
                texCoordsBuffer.push_back(uv);
                colorsBuffer.push_back(c);

                } // per vertex
            } // per triangle
        } // per face

}

void TexturedObject::LeedsGLRender(LeedsGL& leedsGL)
{
    if(vertexBuffer.size() == 0)
        InitializeLeedsGLArrays();
    leedsGL.drawArrays(vertexBuffer,normalsBuffer,texCoordsBuffer,colorsBuffer,leedsGL.TRIANGLES);
}
