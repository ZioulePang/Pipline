#include "LeedsGL.h"
#include <math.h>
#include "omp.h"

#define NAME(variable) (#variable)
using namespace std;
LeedsGL::LeedsGL()
{
    //TODO: Initialize all variables that should be initialized!
}

LeedsGL::~LeedsGL()
{
    //TODO: Free any resources you allocate.
    enabledTexture = nullptr;
}

void LeedsGL::clear(byte mask)
{
    //clear color and depth value
    if(mask == COLORMASK)
    {
        //if color on , clear all of the color into the color we want
        for(int i = 0; i < frameBuffer.width*frameBuffer.height;i++)
        {
            frameBuffer.block[i] =  bufferClearColor;
        }
    }
    if(mask == DEPTHTEST)
    {
        //if depth on ,vic visa
        clear(COLORMASK);

    }
    inputQueue.clear();
    transformedQueue.clear();
    primitivesQueue.clear();
    clippedPrimitivesQueue.clear();
    fragmentQueue.clear();
}

void LeedsGL::setUniform(const string& name,const bool value)
{
    //uniform variables : boolean
    if(NAME(lightingEnabled) == name) lightingEnabled = value;
    if(NAME(textureModulationEnabled) == name) textureModulationEnabled = value;
    if(NAME(texturingEnabled) == name) texturingEnabled = value;
    if(NAME(UVColorDebug) == name) UVColourDebug = value;
}

void LeedsGL::setUniform(const string &name, const Matrix4 &mat)
{
    //uniform variables : matrix
    if(NAME(lightMatrix) == name) lightMatrix = mat;
    if(NAME(viewPortMatrix) == name) viewPortMatrix = mat;
    if(NAME(projectionMatrix) == name) projectionMatrix = mat;
    if(NAME(modelviewMatrix) == name) modelviewMatrix = mat;
}

void LeedsGL::setUniform(const string &name, const RGBAValueF &col)
{
    //uniform variables : RGBAcolour
    if(NAME(lightColour) == name) lightColour = col;
    if(NAME(emissiveMaterial) == name) emissiveMaterial = col;
    if(NAME(ambientMaterial) == name) ambientMaterial = col;
    if(NAME(diffuseMaterial) == name) diffuseMaterial = col;
    if(NAME(specularMaterial) == name) specularMaterial = col;
}

void LeedsGL::setUniform(const std::string &name, const float val)
{
    //uniform variables : float
    if(NAME(shininessMaterial) == name) shininessMaterial = val;
}

void LeedsGL::setUniform(const std::string& name, const Homogeneous4& pos)
{
    //uniform variables : positions
    if(NAME(lightPosition) == name) lightPosition = pos;
}

void LeedsGL::clearColor(const RGBAValueF &col)
{
    //let all of the color equal to the color we want
    bufferClearColor = col;
}

void LeedsGL::resizeBuffers(unsigned const int width, unsigned const int height)
{
    //set buffer's column and row
    frameBuffer.Resize(width, height);
}

Matrix4 LeedsGLUtils::calculateViewportMatrix(float cx, float cy, float width, float height)
{
   Matrix4 m = Matrix4();
   m.SetTranslation(Cartesian3(width/2,height/2,0));

   //caculate viewport matrix
   Matrix4 sc;
   sc.SetScale(width,height,1);
   sc = m * sc;
   return sc;
}

Matrix4 LeedsGLUtils::calculateProjectionOrtho(float left, float right, float bottom, float top, float near, float far)
{
    //right or left handedness may have effects on other parts of your code,
    //such as shading, clipping, culling, etc.
    near = -near;
    far = -far;
    Matrix4 m = Matrix4();
    Matrix4 scale;
    //caculate ottho matarix
    scale.SetIdentity();
    m[0][0] = 2.0/(right - left);
    m[1][1] = 2.0/(top - bottom);
    m[2][2] = 1.0/(far - near);
    m[3][3] = 1.0;

    //caculate scale matrix
    scale[0][3] = -1.0 * (right + left)/2.0;
    scale[1][3] = -1.0 * (top + bottom)/2.0;
    scale[3][2] = 1.0f;
    scale[2][3] = -1.0 * (near);
    scale[3][3] = 0.0f;

    return m * scale;
}

Matrix4 LeedsGLUtils::calculateProjectionFrustum(float left, float right, float bottom, float top, float near, float far)
{
    //right or left handedness may have effects on other parts of your code,
    //such as shading, clipping, culling, etc.
    near = -near;
    far = -far;

    //caculate projection matrix;
    Matrix4 t = Matrix4();
    t[0][0] = near/right;
    t[1][1] = near/top;
    t[2][2] = far / (far - near);
    t[3][2] = 1.0f;
    t[2][3] = -1.0f * far * near / (far - near);
    return t;
}

void LeedsGL::texImage2D(RGBAImage const* textureImage)
{
    if(textureImage!= NULL)
    {
        //set texture
         enabledTexture = textureImage;
    }

}

void LeedsGL::enable(const std::byte function)
{
    switch (function)
    {
    case DEPTHTEST:
        //if depth on , disabled
        depthTestEnabled = true;
        break;
    case PERSPECTIVE:
        //if perspective on , disabled
        perspective = true;
        break;
    }
}

void LeedsGL::disable(const std::byte function)
{
    switch (function)
    {
    case DEPTHTEST:
        //if depth on , disabled
        depthTestEnabled = false;
        break;
    case PERSPECTIVE:
        //if perspective on , disabled
        perspective = false;
        break;
    }
}

void LeedsGL::lineWidth(const float width)
{
    //set line width
    rasterizedLineWidth = width;
}

void LeedsGL::pointSize(const float size)
{
    // set point size
    rasterizedPointSize = size;
}

void LeedsGL::drawArrays(const std::vector<Homogeneous4>& vertices, const std::vector<Homogeneous4>& normals, const std::vector<Cartesian3>& textureCoordinates, const std::vector<RGBAValueF>& colors,std::byte mode)
{
    //Calls the whole pipeline, step by step.
    inputAssembly(vertices,normals,textureCoordinates,colors,inputQueue);
    transformVertices(inputQueue,transformedQueue);
    inputQueue.clear();
    primitiveAssembly(transformedQueue,mode,primitivesQueue);
    transformedQueue.clear();
    clipAndCull(primitivesQueue,mode,clippedPrimitivesQueue);
    primitivesQueue.clear();
    rasterisePrimitives(clippedPrimitivesQueue,mode,fragmentQueue);
    processFragments(fragmentQueue);
    //fragmentQueue.clear();
}

void LeedsGL::inputAssembly(const std::vector<Homogeneous4> &vertices, const std::vector<Homogeneous4> &normals, const std::vector<Cartesian3> &textureCoordinates, const std::vector<RGBAValueF> &colors, std::vector<InputVertex> &result)
{
    //This function should combine this disjoint information into a series of InputVertex to be processed
    //by the next step.

    //I wanna input vertex has positions, normals, textures,etc.
    //so...
    int i = 0;
    result.resize(vertices.size());

   //set Parallelization
#pragma omp parallel for schedule(dynamic)
    //put verteices into input vertex
    for (i = 0;i<vertices.size();i++)
    {
       InputVertex input;
       input.position = vertices[i];
       result[i] = input;
    }
    //set Parallelization
#pragma omp parallel for schedule(dynamic)
    //put normals into input vertex
    for (i = 0;i<normals.size();i++)
    {
        Homogeneous4 temp_normal = Homogeneous4(normals[i].x,normals[i].y,normals[i].z,0);
        result[i].normal = temp_normal;
    }
    //set Parallelization
#pragma omp parallel for schedule(dynamic)
    //put texture into input vertex
    for (i = 0;i<textureCoordinates.size();i++)
    {
        result[i].texCoords = textureCoordinates[i];
    }
    //set Parallelization
#pragma omp parallel for schedule(dynamic)
    //put color into input vertex
    for(i = 0;i<colors.size();i++)
    {
        result[i].color = colors[i];
    }
}

void LeedsGL::transformVertices(std::vector<InputVertex> &vertices, std::vector<TransformedVertex>& result)
{
   //Also pass all the necessary information to the next steps of the pipeline.
   //You should check the slides to decide which is the appropriate coordinate system to transform them to.

    //change to left hand cor
    Matrix4 left;
    left.SetIdentity();
    left[2][2] = -1;

    //as I wanna transformedQueue has transformed vertices information, so I need to do matrix things
    //as follow...
    result.resize(vertices.size());
    //set Parallelization
#pragma omp parallel for schedule(dynamic)
    for (int i = 0;i<vertices.size();i++)
    {
        TransformedVertex t;;
        //mvp transformation
        t.p = viewPortMatrix * projectionMatrix * modelviewMatrix * vertices[i].position;
        t.position = t.p.Point();
        //as I need to caculate baricentirc coordinates in vcs, so...
        t.temp_p = ( left * modelviewMatrix * vertices[i].position).Point();
        t.normals = (left * modelviewMatrix* vertices[i].normal);
        //put texture into transformedQueue
        t.texCoords = vertices[i].texCoords;

        //put color into transformedQueue
        t.colors = vertices[i].color;
        result[i] = t;
    }
}

void LeedsGL::primitiveAssembly(std::vector<TransformedVertex> &vertices,std::byte mode, std::vector<Primitive>& result)
{
    //As there have three modes in this project
    //so....
    switch (mode)
    {
    case POINTS:
        result.resize(vertices.size());
        //set Parallelization
        #pragma omp parallel for schedule(dynamic)
        //if it is a point ,just put it in to queue
        for (int i = 0;i<vertices.size();i++)
        {
            Primitive p;
            p.transformedVertices.push_back(vertices[i]);
            result[i] = p;
        }
        break;
    case LINES:
        result.resize(vertices.size()/2);
        //set Parallelization
        #pragma omp parallel for schedule(dynamic)
        // if it is a line, it has two points
        //so we need to push twice
        for (int i = 0;i<vertices.size();i+=2)
        {
            Primitive p;
            p.transformedVertices.push_back(vertices[i]);
            p.transformedVertices.push_back(vertices[i + 1]);
            result[i/2] = p;
        }
        break;
    case TRIANGLES:
        result.resize(vertices.size()/3);
        //set Parallelization
        #pragma omp parallel for schedule(dynamic)
        // if it is a triangle, it has three points
        //so we need to push trible times
        for (int i = 0;i<vertices.size();i+=3)
        {
            Primitive p;
            p.transformedVertices.push_back(vertices[i]);
            p.transformedVertices.push_back(vertices[i + 1]);
            p.transformedVertices.push_back(vertices[i + 2]);
            result[i/3] = p;
        }
        break;
    }

}

void LeedsGL::clipAndCull(std::vector<Primitive>& primitives,std::byte mode, std::vector<Primitive>& result)
{
    //TODO: Implement clipping and culling. Should have a different behavior for each type of primitive.
    //Pay attention to what type of projection you are using, as your clipping planes will be different.
    //If you choose to skip this step as it is one of your last tasks, just return all the same primitives.

    //let's see eye position is 0,0,1
    //so...
    Cartesian3 eye = Cartesian3(0,0,1);

    //it is as same as primitiveAssmbly, it has three types of mode
    switch (mode)
    {
    case POINTS:
        //Point cliping just caculate one piont positon and edge
        for (int i = 0;i < primitives.size();i++)
        {
            TransformedVertex vertex0 = primitives[i].transformedVertices[0];

            //x larger then right edge
            if(vertex0.position.x > 1) vertex0.position = Cartesian3(vertex0.p.w,vertex0.position.y,vertex0.position.z);
            //x smaller then left edge
            else if(vertex0.position.x < -1 ) vertex0.position = Cartesian3(vertex0.p.w * -1,vertex0.position.y,vertex0.position.z);

            //y larger then top edge
            if(vertex0.position.y > 1) vertex0.position = Cartesian3(vertex0.position.x,vertex0.p.w,vertex0.position.z);
            //y smaller then bottom edge
            else if(vertex0.position.y < -1) vertex0.position = Cartesian3(vertex0.position.x,vertex0.p.w * -1,vertex0.position.z);

            //z larger then back edge
            if(vertex0.position.z > 1) vertex0.position = Cartesian3(vertex0.position.x,vertex0.position.y,vertex0.p.w);
            //y smaller then front edge
            else if(vertex0.position.z < 0) vertex0.position = Cartesian3(vertex0.position.x,vertex0.position.y * -1,0);
        }
        break;
    case LINES:
        //Line cliping
        //we need to caculate two points' position and edge line.
        //if points' position over the egdes, we need two clip then and modify their position
        for (int i = 0;i < primitives.size();i++)
        {
            TransformedVertex vertex0 = primitives[i].transformedVertices[0];
            TransformedVertex vertex1 = primitives[i].transformedVertices[1];

            float y_right = ((-1 - vertex1.position.x) / (vertex0.position.x - vertex1.position.x)) * (vertex0.position.y - vertex1.position.y) + vertex1.position.y;
            float y_left = ((1 - vertex1.position.x) / (vertex0.position.x - vertex1.position.x)) * (vertex0.position.y - vertex1.position.y) + vertex1.position.y;

            float x_top = ((1 - vertex1.position.y) / (vertex0.position.y - vertex1.position.y)) * (vertex0.position.x - vertex1.position.x) + vertex1.position.x;
            float x_bottom = ((-1 - vertex1.position.y) / (vertex0.position.y - vertex1.position.y)) * (vertex0.position.x - vertex1.position.x) + vertex1.position.x;

            float z_front = ((0 - vertex1.position.x) / (vertex0.position.x - vertex1.position.x)) * (vertex0.position.z - vertex1.position.z) + vertex1.position.z;
            float z_back = ((1 - vertex1.position.x) / (vertex0.position.x - vertex1.position.x)) * (vertex0.position.z - vertex1.position.z) + vertex1.position.z;


            //x larger then right edge
            if(vertex0.position.x > 1) vertex0.position = Cartesian3(vertex0.p.w,y_left,vertex0.position.z);
            //x smaller then left edge
            else if(vertex0.position.x < -1 ) vertex0.position = Cartesian3(vertex0.p.w * -1,y_right,vertex0.position.z);

            //y larger then top edge
            if(vertex0.position.y > 1) vertex0.position = Cartesian3(x_top,vertex0.p.w,vertex0.position.z);
            //y larger then top edge
            else if(vertex0.position.y < -1 ) vertex0.position = Cartesian3(x_bottom,vertex0.p.w * -1,vertex0.position.z);

            //z larger then back edge
            if(vertex0.position.z > vertex0.p.w) vertex0.position = Cartesian3(vertex0.position.x,vertex0.position.y,vertex0.p.w);
            //y smaller then front edge
            else if(vertex0.position.z < 0) vertex0.position = Cartesian3(vertex0.position.x,vertex0.position.y * -1,0);

            //x larger then right edge
            if(vertex1.position.x >1) vertex1.position = Cartesian3(vertex1.p.w,y_left,vertex1.position.z);
            //x smaller then left edge
            else if(vertex1.position.x < -1 ) vertex1.position = Cartesian3(vertex1.p.w * -1,y_right,vertex1.position.z);

            //y larger then top edge
            if(vertex1.position.y >1) vertex1.position = Cartesian3(x_top,vertex1.p.w,vertex1.position.z);
            //y larger then top edge
            else if(vertex1.position.y < -1 ) vertex1.position = Cartesian3(x_bottom,vertex1.p.w * -1,vertex1.position.z);

            //z larger then back edge
            if(vertex1.position.z > 1) vertex1.position = Cartesian3(vertex1.position.x,vertex1.position.y,vertex1.p.w);
            //z smaller then front edge
            else if(vertex1.position.z < 0) vertex1.position = Cartesian3(vertex1.position.x,vertex1.position.y * -1,0);
        }



        break;
    case TRIANGLES:
        //Triangle cliping



        //Triangle culling
        //if this triangle is back to our eye, then cull them
        //so.....
        for (int i = 0;i < primitives.size();i++)
        {
            TransformedVertex vertex0 = primitives[i].transformedVertices[0];
            TransformedVertex vertex1 = primitives[i].transformedVertices[1];
            TransformedVertex vertex2 = primitives[i].transformedVertices[2];

            Cartesian3 v0v1 = vertex1.position-vertex0.position;
            Cartesian3 v1v2 = vertex2.position-vertex1.position;
            Cartesian3 v2v0 = vertex0.position-vertex2.position;

            //caculate direction
            Cartesian3 normal1 = v0v1.cross(v1v2);
            Cartesian3 normal2 = v1v2.cross(v2v0);
            Primitive temp;
            //if back to us
            // then...
            if(normal1.unit().dot(eye) != -1 && normal2.unit().dot(Cartesian3(0,0,1)) != -1)
            {
                temp = primitives[i];
                result.push_back(temp);
            }

        }
        break;
    }

    result = primitives;
}

void LeedsGL::rasterisePrimitives(std::vector<Primitive> &primitives, std::byte mode, std::vector <Fragment>& results)
{
    //as this project has three modes
    //so....
    switch (mode)
    {
    case POINTS:

        //if it is a point, just rasterise it
        for (int i = 0;i<primitives.size();i++)
        {
            //raterise points
            rasterisePoint(primitives[i],results);
        }
        break;
    case LINES:
        //if it is a line, just rasterise it
        for (int i = 0;i<primitives.size();i++)
        {
            //raterise lines
            rasteriseLine(primitives[i],results);
        }
        break;
    case TRIANGLES:
        //if it is a triangle, just rasterise it
        for (int i = 0;i<primitives.size();i++)
        {
            //raterise triangles
            rasteriseTriangle(primitives[i],results);
        }
        break;
    }

}

void LeedsGL::rasterisePoint(const Primitive &point,std::vector<Fragment>& output)
{
    //gain the point we push in
    TransformedVertex v = point.transformedVertices[0];

    //we try to use a cube cover it
    //so we use four points to search it position
    float minX = v.position.x-rasterizedPointSize/2, maxX = v.position.x+rasterizedPointSize/2;
    float minY = v.position.y-rasterizedPointSize/2, maxY = v.position.y+rasterizedPointSize/2;
    if(minX<0) minX = 0;
    if(maxX>frameBuffer.width) maxX = frameBuffer.width;
    if(minY<0) minY = 0;
    if(maxY>frameBuffer.height) maxY = frameBuffer.height;

    //we have got the position
    //then we need to caculate it's information
    for (int i = int(minX); i <= maxX; i++)
    {
        for (int j = int(minY); j <= maxY; j++)
        {
             Fragment rasterFragment;
             rasterFragment.col = i;
             rasterFragment.row = j;
             rasterFragment.color = v.colors;
             rasterFragment.depth= v.temp_p.z;
             output.push_back(rasterFragment);
        }
    }
}

void LeedsGL::rasteriseLine(const Primitive &line, std::vector<Fragment> &output)
{
    //gain the points we needed
     TransformedVertex vertex1 = line.transformedVertices[0];
     TransformedVertex vertex2 = line.transformedVertices[1];

     Cartesian3 v1 = vertex1.position;
     Cartesian3 v2 = vertex2.position;

     //as we need to use depth so we need a vcs position
     Cartesian3 temp_v1 = vertex1.temp_p;
     Cartesian3 temp_v2 = vertex2.temp_p;


     //caculate length
     Cartesian3 v = v2 - v1;
     Cartesian3 n = Cartesian3(v.y,-1 * v.x,0);
     n = n/n.length();

     //caculate color
     Cartesian3 color1 = Cartesian3(vertex1.colors.red,vertex1.colors.green,vertex1.colors.blue);
     Cartesian3 color2 = Cartesian3(vertex2.colors.red,vertex2.colors.green,vertex2.colors.blue);
     Cartesian3 total_c = color2 - color1;
     RGBAValueF temp_c = RGBAValueF(total_c.x,total_c.y,total_c.z,1);

     //we have got the position
     //then we need to caculate it's information
     for (int i = 0;i<rasterizedLineWidth/2;i++)
     {
         //increase each point's width
         //so we use normal to enhance their length
         Cartesian3 up_v1 = v1 + n * i;
         Cartesian3 low_v1 = v1 - n * i;

         //increase each point's width
         //so we use normal to enhance their length
         Cartesian3 up_v2 = v2 + n * i;
         Cartesian3 low_v2 = v2 - n * i;

         Cartesian3 temp_up_v1 = temp_v1 + n * i;
         Cartesian3 temp_low_v1 = temp_v1 - n * i;

         Cartesian3 temp_up_v2 = temp_v2 + n * i;
         Cartesian3 temp_low_v2 = temp_v2 - n * i;

         //as we need to caculate every points' position between these two points
         //so....
         for (float t = 0;t<1.0;t+=0.001)
         {
             //caculate these two points information
             Fragment f,f1,f2;
             f.col = (v1.x + (v2.x-v1.x) * t);
             f.row = (v1.y + (v2.y-v1.y) * t);
             f.color = vertex1.colors + operator*(t,temp_c);
             f.depth = (vertex1.temp_p.z + (vertex2.temp_p.z-vertex1.temp_p.z) * t);

             //caculate lower lines information
             f1.col = (up_v1.x + (up_v2.x-up_v1.x) * t);
             f1.row = (up_v1.y + (up_v2.y-up_v1.y) * t);
             f1.color = vertex1.colors + operator*(t,temp_c);
             f.depth = (vertex1.temp_p.z + (vertex2.temp_p.z-vertex1.temp_p.z) * t);

             //caculate upper lines information
             f2.col = (low_v1.x + (low_v2.x-low_v1.x) * t);
             f2.row = (low_v1.y + (low_v2.y-low_v1.y) * t);
             f2.color = vertex1.colors + operator*(t,temp_c);
             f.depth = (vertex1.temp_p.z + (vertex2.temp_p.z-vertex1.temp_p.z) * t);

             //push them back
             output.push_back(f);
             output.push_back(f1);
             output.push_back(f2);
         }
     }

}

float LeedsGLUtils::distancePointLine(Cartesian3 r, Cartesian3 n, Cartesian3 p)
{
    //assumes n is normalized
       return n.dot(r) - n.dot(p);
}

RGBAValueF LeedsGL::CalculateLighting(const Homogeneous4& n_vcs, const Homogeneous4& v_vcs,const RGBAValueF& em, const RGBAValueF& am, const RGBAValueF& diff, const RGBAValueF& spec, float shin)
{
    Matrix4 left;
    left.SetIdentity();
    left[2][2] = -1;

    if((n_vcs.x == 0.0f && n_vcs.y == 0.0f && n_vcs.z == 0.0f)) // we shouldn't try shading if there are no normals
    return RGBAValueF();

    Cartesian3 lightVector;
    Cartesian3 unitNormal = n_vcs.Vector().unit();


    //Directional Light
    Homogeneous4 lp = left * lightMatrix * lightPosition;
    if(abs(lp.w - 0) < std::numeric_limits<float>::epsilon())
        lightVector = Cartesian3(0,0,0) - lp.Vector().unit();
    else //point light
        lightVector = (lp - v_vcs).Vector().unit();
    Cartesian3 eyeVector = perspective? -1*v_vcs.Point(): Cartesian3(0,0,1);
    Cartesian3 bisector = (lightVector + eyeVector).unit();

    RGBAValueF emissive = em;

    RGBAValueF ambient =  am.modulate(ambientMaterial);
    float dDot = unitNormal.dot(lightVector);
    dDot = dDot <0? 0: dDot;

    RGBAValueF diffuse = dDot * diff.modulate(lightColour);

    float sDot = unitNormal.dot(bisector);


    sDot = sDot <0? 0:sDot;
    sDot = pow(sDot,shin);
    //sDot = ((f.shininess+2)/8.0f)*sDot*dDot;
    sDot = dDot>0? sDot : 0;
    sDot = sDot* dDot *(shin+2)/2*float(M_PI);

    Cartesian3 fs = Cartesian3(spec.red,spec.green,spec.blue);

    Cartesian3 air = Cartesian3(1,1,1);
    Cartesian3 a = (air - fs);
    Cartesian3 b = (air+fs);
    Cartesian3 r0 = Cartesian3(a.x/b.x,a.y/b.y,a.z/b.z);

    r0 = Cartesian3(r0.x*r0.x,r0.y*r0.y,r0.z*r0.z);
    Cartesian3 rschlick = r0 + (air-r0) * powf((1-bisector.dot(lightVector)),5);

    RGBAValueF updatedSpecular = RGBAValueF(rschlick.x,rschlick.y,rschlick.z,1);
    RGBAValueF specular = sDot * updatedSpecular.modulate(lightColour);
    return emissive + ambient + diffuse + specular;
}

void LeedsGL::rasteriseTriangle(const Primitive &triangle, std::vector<Fragment> &output)
{
    TransformedVertex vertex0 = triangle.transformedVertices[0];
    TransformedVertex vertex1 = triangle.transformedVertices[1];
    TransformedVertex vertex2 = triangle.transformedVertices[2];

    // compute a bounding box that starts inverted to frame size
    // clipping will happen in the raster loop proper
    float minX = frameBuffer.width, maxX = 0.0;
    float minY = frameBuffer.height, maxY = 0.0;

    // test against all vertices
    if (vertex0.position.x < minX) minX = vertex0.position.x;
    if (vertex0.position.x > maxX) maxX = vertex0.position.x;
    if (vertex0.position.y < minY) minY = vertex0.position.y;
    if (vertex0.position.y > maxY) maxY = vertex0.position.y;

    if (vertex1.position.x < minX) minX = vertex1.position.x;
    if (vertex1.position.x > maxX) maxX = vertex1.position.x;
    if (vertex1.position.y < minY) minY = vertex1.position.y;
    if (vertex1.position.y > maxY) maxY = vertex1.position.y;

    if (vertex2.position.x < minX) minX = vertex2.position.x;
    if (vertex2.position.x > maxX) maxX = vertex2.position.x;
    if (vertex2.position.y < minY) minY = vertex2.position.y;
    if (vertex2.position.y > maxY) maxY = vertex2.position.y;

    Cartesian3 v0 = Cartesian3(vertex0.position.x,vertex0.position.y,0);
    Cartesian3 v1 = Cartesian3(vertex1.position.x,vertex1.position.y,0);
    Cartesian3 v2 = Cartesian3(vertex2.position.x,vertex2.position.y,0);

    Cartesian3 v0v1 = v1-v0;
   Cartesian3 n_v0v1 = Cartesian3(-v0v1.y,v0v1.x,0);
   Cartesian3 v1v2 = v2-v1;
   Cartesian3 n_v1v2 = Cartesian3(-v1v2.y,v1v2.x,0);
   Cartesian3 v2v0 = v0-v2;
   Cartesian3 n_v2v0 = Cartesian3(-v2v0.y,v2v0.x,0);
   float dAlpha =  LeedsGLUtils::distancePointLine(v0, n_v1v2,v1);
   float dBeta = LeedsGLUtils::distancePointLine(v1, n_v2v0,v2);
   float dGamma = LeedsGLUtils::distancePointLine(v2,n_v0v1,v0);

   if (abs(dAlpha-0)<std::numeric_limits<float>::epsilon() ||
       abs(dBeta-0)<std::numeric_limits<float>::epsilon() ||
       abs(dGamma-0)<std::numeric_limits<float>::epsilon())
        return;

    // create a fragment for reuse
    Fragment rasterFragment;

    //we use some parameters to preserve light color and texture color
    RGBAValueF blinn_phong;
    RGBAValueF texture;

    // loop through the pixels in the bounding box
    for (rasterFragment.row = int(minY); rasterFragment.row <= maxY; rasterFragment.row++)
    { // per row
        // this is here so that clipping works correctly
        if (rasterFragment.row < 0) continue;
        if (rasterFragment.row >= int(frameBuffer.height)) continue;
        for (rasterFragment.col = int(minX); rasterFragment.col <= maxX; rasterFragment.col++)
        { // per pixel
            // this is also for correct clipping
            if (rasterFragment.col < 0) continue;
            if (rasterFragment.col >= int(frameBuffer.width)) continue;

            // the pixel in cartesian format
            Cartesian3 pixel(rasterFragment.col+0.5f, rasterFragment.row+0.5f, 0.0f);

            // right - we have a pixel inside the frame buffer AND the bounding box
            // note we *COULD* compute gamma = 1.0 - alpha - beta instead
            float alpha = LeedsGLUtils::distancePointLine(pixel,n_v1v2,v1) / dAlpha;
            float beta = LeedsGLUtils::distancePointLine(pixel,n_v2v0,v2)/dBeta;
            float gamma = LeedsGLUtils::distancePointLine(pixel,n_v0v1,v0)/dGamma;

            // now perform the half-plane test
            if ((alpha < 0.0f) || (beta < 0.0f) || (gamma < 0.0f))
                continue;

            //caculate each point's depth vaule
            float w = 1.0f / ((alpha/vertex0.temp_p.z) + (beta/vertex1.temp_p.z)+ (gamma/vertex2.temp_p.z));
            rasterFragment.depth= w;

            //caculate baricentric normals , positions and uv coordinates
            Homogeneous4 bc_n = vertex0.normals * alpha + vertex1.normals * beta + vertex2.normals * gamma;
            //caculate baricentric normals , positions and uv coordinates
            Homogeneous4 bc_p = vertex0.temp_p * alpha + vertex1.temp_p * beta + vertex2.temp_p * gamma;
            //caculate baricentric normals , positions and uv coordinates
            Homogeneous4 bc = vertex0.position * alpha + vertex1.position * beta + vertex2.position * gamma;
            //caculate baricentric normals , positions and uv coordinates
            Cartesian3 uv = ((vertex0.texCoords/vertex0.temp_p.z) * alpha + (vertex1.texCoords/vertex1.temp_p.z) * beta + (vertex2.texCoords/vertex2.temp_p.z) * gamma) * w;

            // if we click texture ,light and modulation
            if (texturingEnabled && lightingEnabled && textureModulationEnabled)
            {
                //do texture part
                // ensure pointer is not null
                if(enabledTexture != NULL)
                {
                    //then we caculate light color and preserve it
                    blinn_phong = LeedsGL::CalculateLighting(bc_n,bc_p,emissiveMaterial,ambientMaterial,diffuseMaterial,specularMaterial,shininessMaterial);

                    //use linear interpolation to caculte each part texture color
                    int i = uv.x * enabledTexture->width;                    // truncates s to get i
                    int j = uv.y * enabledTexture->height;                    // truncates t to get j
                    float sParm = uv.x * enabledTexture->width - i;          // compute s for interpolation
                    float tParm = uv.y * enabledTexture->height - j;          // compute t for interpolation

                    // grab four nearest texel colours
                    RGBAValueF colour00 = RGBAValueF(float((*enabledTexture)[i][j].red)/255.0f,float((*enabledTexture)[i][j].green)/255.0f,float((*enabledTexture)[i][j].blue)/255.0f);
                    RGBAValueF colour01 = RGBAValueF(float((*enabledTexture)[i][j+1].red)/255.0f,float((*enabledTexture)[i][j+1].green)/255.0f,float((*enabledTexture)[i][j+1].blue)/255.0f);
                    RGBAValueF colour10 = RGBAValueF(float((*enabledTexture)[i + 1][j].red)/255.0f,float((*enabledTexture)[i + 1][j].green)/255.0f,float((*enabledTexture)[i + 1][j].blue)/255.0f);
                    RGBAValueF colour11 = RGBAValueF(float((*enabledTexture)[i + 1][j+1].red)/255.0f,float((*enabledTexture)[i + 1][j+1].green)/255.0f,float((*enabledTexture)[i + 1][j+1].blue)/255.0f);

                    Cartesian3 color_00 = Cartesian3(colour00.red,colour00.green,colour00.blue);
                    Cartesian3 color_01 = Cartesian3(colour01.red,colour01.green,colour01.blue);
                    Cartesian3 color_10 = Cartesian3(colour10.red,colour10.green,colour10.blue);
                    Cartesian3 color_11 = Cartesian3(colour11.red,colour11.green,colour11.blue);

                    // compute colours on edges
                    Cartesian3 colour0 = color_00 + tParm * (color_01 - color_00);
                    Cartesian3 colour1 = color_10 + tParm * (color_11 - color_10);

                    Cartesian3 result = colour1 + sParm * (colour1 - colour0);

                    // compute colour for interpolated texel
                    texture = RGBAValueF(result.x,result.y,result.z);

                    //finally,we modulate light and texture
                    rasterFragment.color = blinn_phong.modulate(texture);
                }

            }
            //if we only click ligth on
            else if(lightingEnabled)
            {
                //caculate light color
                blinn_phong = LeedsGL::CalculateLighting(bc_n,bc_p,emissiveMaterial,ambientMaterial,diffuseMaterial,specularMaterial,shininessMaterial);
                rasterFragment.color = blinn_phong;

            }
            //if we only click texture
            else if(texturingEnabled)
            {
                // ensure pointer is not null
                if(enabledTexture != NULL)
                {
                    //use linear interpolation to caculte each part texture color
                    int i = uv.x * enabledTexture->width;                    // truncates s to get i
                    int j = uv.y * enabledTexture->height;                    // truncates t to get j
                    float sParm = uv.x * enabledTexture->width - i;          // compute s for interpolation
                    float tParm = uv.y * enabledTexture->height - j;          // compute t for interpolation

                    // grab four nearest texel colours
                    RGBAValueF colour00 = RGBAValueF(float((*enabledTexture)[i][j].red)/255.0f,float((*enabledTexture)[i][j].green)/255.0f,float((*enabledTexture)[i][j].blue)/255.0f);
                    RGBAValueF colour01 = RGBAValueF(float((*enabledTexture)[i][j+1].red)/255.0f,float((*enabledTexture)[i][j+1].green)/255.0f,float((*enabledTexture)[i][j+1].blue)/255.0f);
                    RGBAValueF colour10 = RGBAValueF(float((*enabledTexture)[i + 1][j].red)/255.0f,float((*enabledTexture)[i + 1][j].green)/255.0f,float((*enabledTexture)[i + 1][j].blue)/255.0f);
                    RGBAValueF colour11 = RGBAValueF(float((*enabledTexture)[i + 1][j+1].red)/255.0f,float((*enabledTexture)[i + 1][j+1].green)/255.0f,float((*enabledTexture)[i + 1][j+1].blue)/255.0f);

                    Cartesian3 color_00 = Cartesian3(colour00.red,colour00.green,colour00.blue);
                    Cartesian3 color_01 = Cartesian3(colour01.red,colour01.green,colour01.blue);
                    Cartesian3 color_10 = Cartesian3(colour10.red,colour10.green,colour10.blue);
                    Cartesian3 color_11 = Cartesian3(colour11.red,colour11.green,colour11.blue);

                    // compute colours on edges
                    Cartesian3 colour0 = color_00 + tParm * (color_01 - color_00);
                    Cartesian3 colour1 = color_10 + tParm * (color_11 - color_10);

                    Cartesian3 result = colour1 + sParm * (colour1 - colour0);

                    // compute colour for interpolated texel
                    texture = RGBAValueF(result.x,result.y,result.z);
                    rasterFragment.color = texture;
                }

            }
            //if we only click UVcolor
            else if(UVColourDebug)
            {
                //caculate baricentric color
                rasterFragment.color =RGBAValueF(uv.x,uv.y,w);
            }
            else
            {
                RGBAValueF result_color = operator*(alpha,vertex0.colors) +  operator*(beta,vertex1.colors) + operator*(gamma,vertex2.colors);
                rasterFragment.color = result_color;
            }

            output.push_back(rasterFragment);
        } // per pixel
    } // per row

}

void LeedsGL::processFragments(std::vector<Fragment> &fragments)
{
   //Depth test should go here. We don't explicitly have a pre or post fragment stage.
   //Consider the "shading" as the fragment stage. Decide if the depth test should go before or after, and justify.

    //set depth buffer equals to 1000
    float depthBuffer[frameBuffer.height][frameBuffer.width];

    //set each depth euqals to a very large number
    for (int i = 0;i<frameBuffer.height;i++)
    {
        for (int j = 0;j<frameBuffer.width;j++)
        {
            depthBuffer[i][j] = 1000.0f;
        }
    }

    //do loop for each fragments
    for (int i = 0;i<fragments.size();i++)
    {
        //if our model our of boundary , we do not render them
        if(fragments[i].row <0 || fragments[i].row>frameBuffer.height || fragments[i].col<0 || fragments[i].col >frameBuffer.width) continue;

        //use depth test
        if(depthTestEnabled)
        {
            //when our depth smaller than depth buffer value
            if(depthBuffer[fragments[i].row][fragments[i].col] > fragments[i].depth)
            {
                //replace color
                frameBuffer[fragments[i].row][fragments[i].col] = fragments[i].color;
                //replace depth value;
                depthBuffer[fragments[i].row][fragments[i].col] = fragments[i].depth;
            }
        }
        else
        {
             frameBuffer[fragments[i].row][fragments[i].col] = fragments[i].color;
        }
    }


}

