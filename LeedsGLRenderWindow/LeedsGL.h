#ifndef LEEDSGL_H
#define LEEDSGL_H

#include <string>
#include <map>
#include <vector>
#include <memory>
#include <cstddef>
#include "Matrix4.h"
#include "RGBAValue.h"
#include "RGBAValueF.h"
#include "RGBAImage.h"

namespace LeedsGLUtils
{
    Matrix4 calculateViewportMatrix(float cx, float cy, float width, float height);
    Matrix4 calculateProjectionFrustum(float left, float right,float bottom, float top, float near, float far);
    Matrix4 calculateProjectionOrtho(float left, float right,float bottom, float top, float near, float far);
    float distancePointLine(Cartesian3 r, Cartesian3 n, Cartesian3 p);
}

// class with vertex attributes
struct InputVertex
{
    // I wanna position normal and texture
    //so...
    Homogeneous4 position;
    Homogeneous4 normal;
    Cartesian3 texCoords;

    //I also wanna color
    RGBAValueF color;
};

struct TransformedVertex
{
    //as we need to caculate different things like light and texture
    //so we need different point's in different coordinates system
    //so I preseve position twice
    Cartesian3 position;
    Homogeneous4 normals;
    Cartesian3 texCoords;
    RGBAValueF colors;
    Cartesian3 temp_p;
    Homogeneous4 p;
};

struct Primitive
{
     std::vector<TransformedVertex> transformedVertices;
};


struct Fragment
{
    //as we need to caculate each fragment color
    //so we need row and column
    int row;
    int col;
    //and color
    RGBAValueF color;

    //finally, we need do deep test
    float depth;
};


class LeedsGL
{
public:
    LeedsGL();
    ~LeedsGL();

    //RENDERING PARAMETERS:
    void setUniform(const std::string& name,const bool value);
    void setUniform(const std::string& name, const Matrix4& mat);
    void setUniform(const std::string& name, const RGBAValueF& col);
    void setUniform(const std::string& name, const Homogeneous4& pos);
    void setUniform(const std::string& name, const float val);

    //PIPELINE CONTROL:
    void clearColor(const RGBAValueF &col);
    void clear(std::byte mask);
    void enable(const std::byte function);
    void disable(const std::byte function);
    void texImage2D(RGBAImage const* textureImage);
    void resizeBuffers(unsigned const int width,unsigned const int height);
    void lineWidth(const float width);
    void pointSize(const float size);

    //MAIN PIPELINE IMPLEMENTATION
    void drawArrays(const std::vector<Homogeneous4>& vertices,
                    const std::vector<Homogeneous4>& normals,
                    const std::vector<Cartesian3>& textureCoordinates,
                    const std::vector<RGBAValueF>& colors,std::byte mode);
    void inputAssembly(const std::vector<Homogeneous4>& vertices,
                       const std::vector<Homogeneous4>& normals,
                       const std::vector<Cartesian3>& textureCoordinates,
                       const std::vector<RGBAValueF>& colors,
                       std::vector<InputVertex>& result);
    void transformVertices(std::vector<InputVertex>& vertices,
                           std::vector<TransformedVertex>& result);
    void primitiveAssembly(std::vector<TransformedVertex>& vertices,
                           std::byte mode,
                           std::vector<Primitive>& result);
    void clipAndCull(std::vector<Primitive>& primitives,
                     std::byte mode,
                     std::vector<Primitive>& result);
    void rasterisePrimitives(std::vector<Primitive>& primitives,
                             std::byte mode,
                             std::vector<Fragment>& result);

    void rasterisePoint(const Primitive& point,
                        std::vector<Fragment>& output);
    void rasteriseLine(const Primitive& line,
                       std::vector<Fragment>& output);
    void rasteriseTriangle(const Primitive& triangle,
                           std::vector<Fragment>& output);
    void processFragments(std::vector<Fragment>& fragments);

    RGBAValue BilinearLookup(float s, float t);

    //SHADING.
    RGBAValueF CalculateLighting(const Homogeneous4& n_vcs,
                                 const Homogeneous4& v_vcs,
                                 const RGBAValueF& em,
                                 const RGBAValueF& am,
                                 const RGBAValueF& diff,
                                 const RGBAValueF& spec,
                                 float shin);
    //BUFFERS
    RGBAImage frameBuffer;
    RGBAImage depthBuffer;
    //Masks
    static const std::byte COLORMASK{1};
    static const std::byte DEPTHMASK{2};
    //Function constants
    static const std::byte DEPTHTEST{3};
    static const std::byte PERSPECTIVE{4};
    //Drawing modes

    static const std::byte POINTS{5};
    static const std::byte LINES{6};
    static const std::byte TRIANGLES{7};

private:

//uniform variables
    RGBAImage const* enabledTexture = nullptr;
    bool texturingEnabled;
    bool textureModulationEnabled;
    bool UVColourDebug;
    bool lightingEnabled;

    Matrix4 lightMatrix;
    Matrix4 viewPortMatrix;
    Matrix4 projectionMatrix;
    Matrix4 modelviewMatrix;

    Homogeneous4 lightPosition;
    RGBAValueF lightColour;

    RGBAValueF emissiveMaterial;
    RGBAValueF ambientMaterial;
    RGBAValueF diffuseMaterial;
    RGBAValueF specularMaterial;
    float shininessMaterial;

//global states
    float rasterizedLineWidth;
    float rasterizedPointSize;
    RGBAValueF bufferClearColor;
    bool depthTestEnabled;
    bool perspective;

//Queues
    std::vector<InputVertex> inputQueue;
    std::vector<TransformedVertex> transformedQueue;
    std::vector<Primitive> primitivesQueue;
    std::vector<Primitive> clippedPrimitivesQueue;
    std::vector<Fragment> fragmentQueue;

};

#endif // LEEDSGL_H
