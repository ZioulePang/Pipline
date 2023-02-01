#include "LeedsGLRenderWidget.h"

LeedsGLRenderWidget::LeedsGLRenderWidget(TexturedObject      *newTexturedObject,
                                         // the render parameters to use
                                         RenderParameters    *newRenderParameters,
                                         // parent widget in visual hierarchy
                                         QWidget             *parent)

// start by calling inherited constructor with parent widget's pointer
    :QOpenGLWidget (parent),
// then store the pointers that were passed in
    texturedObject(newTexturedObject),
    renderParameters(newRenderParameters)
{
    // leaves nothing to put into the constructor body
}

// destructor
LeedsGLRenderWidget::~LeedsGLRenderWidget()
{
// all of our pointers are to data owned by another class
// so we have no responsibility for destruction
// and OpenGL cleanup is taken care of by Qt
}

void LeedsGLRenderWidget::initializeGL()
{
    //Set clear color
    leedsGL.clearColor(RGBAValueF(0.8f, 0.8f, 0.6f, 1.0f));
    texturedObject->TransferAssetsToLeedsGL(leedsGL);

    axisColors.push_back(RGBAValueF(0.4f,0.0f,0.0f));
    axisCoords.push_back(Cartesian3(0.0f,0.0f,0.0f));
    axisColors.push_back(RGBAValueF(1.0f,0.0f,0.0f));
    axisCoords.push_back(Cartesian3(1.0f,0.0f,0.0f));

    // Y axis is green
    axisColors.push_back(RGBAValueF(0.0f,0.4f,0.0f));
    axisCoords.push_back(Cartesian3(0.0f,0.0f,0.0f));
    axisColors.push_back(RGBAValueF(0.0f,1.0f,0.0f));
    axisCoords.push_back(Cartesian3(0.0f,1.0f,0.0f));

    // Z axis is blue
    axisColors.push_back(RGBAValueF(0.0f,0.0f,0.4f));
    axisCoords.push_back(Cartesian3(0.0f,0.0f,0.0f));
    axisColors.push_back(RGBAValueF(0.0f,0.0f,1.0f));
    axisCoords.push_back(Cartesian3(0.0f,0.0f,1.0f));
}

// called every time the widget is resized
void LeedsGLRenderWidget::resizeGL(int w, int h)
{

    w = std::max(w,0);
    h = std::max(h,0);

    // resize the render images
    leedsGL.resizeBuffers(static_cast<unsigned int>(w), static_cast<unsigned int>(h));
    // reset the Matrices
    Matrix4 VP =  LeedsGLUtils::calculateViewportMatrix(0, 0, w, h);
    Matrix4 P;
    // compute the aspect ratio of the widget
    float aspectRatio =  float(w) / float(h);

    if(renderParameters->perspective){
        //We need this for lighting!
        leedsGL.enable(LeedsGL::PERSPECTIVE);
        if (aspectRatio > 1.0f)
            P = LeedsGLUtils::calculateProjectionFrustum(-0.01f*aspectRatio, 0.01f*aspectRatio, -0.01f, 0.01f, 0.1f, 10.0f);
        // otherwise, make left & right -1.0 & 1.0
        else
            P =  LeedsGLUtils::calculateProjectionFrustum(-0.01f, 0.01f, -0.01f/aspectRatio, 0.01f/aspectRatio, 0.1f, 10.0f);
    }else{
        //We need this for lighting!
        leedsGL.disable(LeedsGL::PERSPECTIVE);
        if (aspectRatio > 1.0f)
            P =  LeedsGLUtils::calculateProjectionOrtho(-aspectRatio, aspectRatio, -1.0f, 1.0f, -1.0f, 1.0f);
        // otherwise, make left & right -1.0 & 1.0
        else
            P =  LeedsGLUtils::calculateProjectionOrtho(-1.0f, 1.0f, -1.0f/aspectRatio, 1.0f/aspectRatio, -1.0f, 1.0f);
    }

    leedsGL.setUniform("viewPortMatrix",VP);
    leedsGL.setUniform("projectionMatrix",P);

}
void LeedsGLRenderWidget::paintGL()
{
    // set background colour to white
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);
    // call the paintFakeGL() routine to prepare the image
    paintLeedsGL();

    // and display the image
    glDrawPixels(int(leedsGL.frameBuffer.width), int(leedsGL.frameBuffer.height), GL_RGBA, GL_UNSIGNED_BYTE, leedsGL.frameBuffer.block);


}

void LeedsGLRenderWidget::paintLeedsGL()
    { // FakeGLRenderWidget::paintFakeGL()

    Matrix4 P;

    float aspectRatio =  float(leedsGL.frameBuffer.width) / float(leedsGL.frameBuffer.height);
    if(renderParameters->perspective){
        leedsGL.enable(LeedsGL::PERSPECTIVE);
        if (aspectRatio > 1.0f)
            P = LeedsGLUtils::calculateProjectionFrustum(-0.01f*aspectRatio, 0.01f*aspectRatio, -0.01f, 0.01f, 0.01f, 10.0f);
        // otherwise, make left & right -1.0 & 1.0
        else
            P =  LeedsGLUtils::calculateProjectionFrustum(-0.01f, 0.01f, -0.01f/aspectRatio, 0.01f/aspectRatio, 0.01f, 10.0f);
    }else{
        leedsGL.disable(LeedsGL::PERSPECTIVE);
        if (aspectRatio > 1.0f)
            P =  LeedsGLUtils::calculateProjectionOrtho(-aspectRatio, aspectRatio, -1.0f, 1.0f, 0.0f, 1.0f);
        // otherwise, make left & right -1.0 & 1.0
        else
            P =  LeedsGLUtils::calculateProjectionOrtho(-1.0f, 1.0f, -1.0f/aspectRatio, 1.0f/aspectRatio, 0.0f, 1.0f);
    }
    leedsGL.setUniform("projectionMatrix",P);

    if (renderParameters->depthTestOn)
        leedsGL.enable(LeedsGL::DEPTHTEST);
    else
        leedsGL.disable(LeedsGL::DEPTHTEST);

    // clear the buffer
    leedsGL.clear(LeedsGL::COLORMASK | (renderParameters->depthTestOn ? LeedsGL::DEPTHMASK: std::byte{0}));

    // if lighting is turned on
    if (renderParameters->useLighting)
        { // use lighting
        // make sure lighting is on
        leedsGL.setUniform("lightingEnabled",true);

        // set light position first, pushing/popping matrix so that it the transformation does
        // not affect the position of the geometric object
        Matrix4 mt;
        mt.SetIdentity();
        //Im expecting a left hand coordinate system so
        mt = mt * renderParameters->lightMatrix;
        leedsGL.setUniform("lightMatrix",mt);
        // and set them in OpenGL
        leedsGL.setUniform("lightColour", RGBAValueF(1.0f,1.0f,1.0f));
        Homogeneous4 LP = Homogeneous4(renderParameters->lightPosition[0],renderParameters->lightPosition[1],renderParameters->lightPosition[2],renderParameters->lightPosition[3]);
        leedsGL.setUniform("lightPosition",LP);
     } // use lighting
    else {
        leedsGL.setUniform("lightingEnabled",false);

    }
    // translate by the visual translation
    Matrix4 M;
    M.SetTranslation(Cartesian3(renderParameters->xTranslate,renderParameters->yTranslate,renderParameters->zTranslate-4.0f));
    //I'm expecting a left hand coordinate system, so:
    Matrix4 R = renderParameters->rotationMatrix;
    M = M * R;


    leedsGL.setUniform("modelviewMatrix",M);

    leedsGL.pointSize(8.0f);
    leedsGL.drawArrays(std::vector<Homogeneous4>{Homogeneous4(0,0,0)},std::vector<Homogeneous4>(),std::vector<Cartesian3>(),std::vector<RGBAValueF>{RGBAValueF(0.3f,0.5f,0.5f)},leedsGL.POINTS);

    // now we start using the render parameters
    if (renderParameters->showAxes)
        { // show axes
        // start with lighting turned off
        leedsGL.setUniform("ligthingEnabled",false);

        // set the lines to be obvious in width
        leedsGL.lineWidth(4.0f);

        leedsGL.drawArrays(axisCoords,std::vector<Homogeneous4>(),std::vector<Cartesian3>(),axisColors,LeedsGL::LINES);
        // reset lighting on if needed
        if (renderParameters->useLighting)
            leedsGL.setUniform("ligthingEnabled",true);
        } // show axes

    // tell the object to draw itself. We set things up before, and the object will call its arrays.
    if (renderParameters->showObject)
    {
        // if we have texturing enabled . . .
        if (renderParameters->texturedRendering)
            { // textures enabled
            // enable textures
            leedsGL.setUniform("texturingEnabled",true);
            // use our other flag to specify replace or modulate
            if (renderParameters->textureModulation)
                leedsGL.setUniform("textureModulationEnabled",true);
            else
                leedsGL.setUniform("textureModulationEnabled",false);
            } // textures enabled
        else
            { // textures disabled
            // make sure that they are disabled
            leedsGL.setUniform("texturingEnabled",false);
            } // textures disabled

        // Scale defaults to the zoom setting
        float scale = 1.0f;
        // if object scaling is requested, apply it as well
        if (renderParameters->scaleObject)
           scale /= texturedObject->objectSize;

        //  now scale everything
        Matrix4 scaleMatrix;
        scaleMatrix.SetScale(scale,scale,scale);

        // apply the translation to the centre of the object if requested
        Matrix4 centerMatrix;
        centerMatrix.SetIdentity();
        if (renderParameters->centreObject){
            Cartesian3 cog = texturedObject->centreOfGravity;
            centerMatrix.SetTranslation(-scale * cog);
        }

        //Now we set the final transformation matrix!
        Matrix4 modelview;
        modelview.SetIdentity();
        modelview = M * scaleMatrix * centerMatrix * modelview;
        leedsGL.setUniform("modelviewMatrix",modelview);

        // we assume a single material for the entire object
        leedsGL.setUniform("emissiveMaterial",RGBAValueF(renderParameters->emissiveMaterial,renderParameters->emissiveMaterial,renderParameters->emissiveMaterial));
        leedsGL.setUniform("ambientMaterial",RGBAValueF(renderParameters->ambientMaterial,renderParameters->ambientMaterial,renderParameters->ambientMaterial));
        leedsGL.setUniform("diffuseMaterial",RGBAValueF(renderParameters->diffuseMaterial,renderParameters->diffuseMaterial,renderParameters->diffuseMaterial));
        leedsGL.setUniform("specularMaterial",RGBAValueF(renderParameters->specularMaterial,renderParameters->specularMaterial,renderParameters->specularMaterial));
        leedsGL.setUniform("shininessMaterial",renderParameters->specularExponent);

        leedsGL.setUniform("UVColorDebug",renderParameters->mapUVWToRGB);

        texturedObject->LeedsGLRender(leedsGL);

        // if we have texturing enabled, turn texturing back off
        leedsGL.setUniform("texturingEnabled",false);
    }

    } // FakeGLRenderWidget::paintFakeGL()


// mouse-handling
void LeedsGLRenderWidget::mousePressEvent(QMouseEvent *event)
    { // FakeGLRenderWidget::mousePressEvent()
    // store the button for future reference
    auto whichButton = event->button();
    // scale the event to the nominal unit sphere in the widget:
    // find the minimum of height & width
    float size = (width() > height()) ? height() : width();
    // scale both coordinates from that
    float x = (2.0f * event->x() - size) / size;
    float y = (size - 2.0f * event->y() ) / size;


    // and we want to force mouse buttons to allow shift-click to be the same as right-click
    auto modifiers = event->modifiers();

    // shift-click (any) counts as right click
    if (modifiers & Qt::ShiftModifier)
        whichButton = Qt::RightButton;

    // send signal to the controller for detailed processing
    emit BeginScaledDrag(int(whichButton), x,y);
    } // FakeGLRenderWidget::mousePressEvent()

void LeedsGLRenderWidget::mouseMoveEvent(QMouseEvent *event)
    { // FakeGLRenderWidget::mouseMoveEvent()
    // scale the event to the nominal unit sphere in the widget:
    // find the minimum of height & width
    float size = (width() > height()) ? height() : width();
    // scale both coordinates from that
    float x = (2.0f * event->x() - size) / size;
    float y = (size - 2.0f * event->y() ) / size;

    // send signal to the controller for detailed processing
    emit ContinueScaledDrag(x,y);
    } // FakeGLRenderWidget::mouseMoveEvent()

void LeedsGLRenderWidget::mouseReleaseEvent(QMouseEvent *event)
    { // FakeGLRenderWidget::mouseReleaseEvent()
    // scale the event to the nominal unit sphere in the widget:
    // find the minimum of height & width
    float size = (width() > height()) ? height() : width();
    // scale both coordinates from that
    float x = (2.0f * event->x() - size) / size;
    float y = (size - 2.0f * event->y() ) / size;

    //std::cout << fakeGL << std::endl;

    // send signal to the controller for detailed processing
    emit EndScaledDrag(x,y);
    } // FakeGLRenderWidget::mouseReleaseEvent()
