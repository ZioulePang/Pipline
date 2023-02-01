/////////////////////////////////////////////////////////////////
//
//  University of Leeds
//  COMP 5812M Foundations of Modelling & Rendering
//  User Interface for Coursework
//
//  September, 2020
//
//  -----------------------------
//  Render Window
//  -----------------------------
//  
//  The render window class is really just a container
//  for tracking the visual hierarchy.  While it is
//  entirely possible to use Qt Creator, I try to avoid
//  over-commitment to it because I need to write code in
//  multiple environments, some of which are not well-suited
//  to IDEs in general, let alone Qt Creator
//
//  Also, by providing sample code, the didactic purpose of 
//  showing how things fit together is better served.
//
/////////////////////////////////////////////////////////////////

// include guard
#ifndef _RENDER_WINDOW_H
#define _RENDER_WINDOW_H

// include standard Qt widgets
#include <QtWidgets>

// include a custom arcball widget 
#include "ArcBallWidget.h"
// include the custom render widget
// and a second widget which the student will edit
#include "LeedsGLRenderWidget.h"

// a window that displays an geometric model with controls
class RenderWindow : public QWidget
    { // class RenderWindow
    private:
    // the geometric object being shown
    TexturedObject              *texturedObject;
    
    // the values set in the interface
    RenderParameters            *renderParameters;

    // window layout    
    QGridLayout                 *windowLayout;

    // custom widgets
    ArcBallWidget               *modelRotator;
    ArcBallWidget               *lightRotator;
    LeedsGLRenderWidget          *leedsGLRenderWidget;

    // standard widgets
    // check boxes to control render options
    QCheckBox                   *depthTestBox;
    QCheckBox                   *lightingBox;
    QCheckBox                   *texturedRenderingBox;
    QCheckBox                   *textureModulationBox;

    // check boxes for modelling options
    QCheckBox                   *showAxesBox;
    QCheckBox                   *showObjectBox;
    QCheckBox                   *centreObjectBox;
    QCheckBox                   *scaleObjectBox;
    QCheckBox                   *mapUVWToRGBBox;
    QCheckBox                   *perspectiveBox;


    // sliders for spatial manipulation
    QSlider                     *xTranslateSlider;
    // we want one slider under each widget
    QSlider                     *yTranslateSlider;
    QSlider                     *zTranslateSlider;

    // sliders for setting lighting parameters
    QSlider                     *emissiveLightSlider;
    QSlider                     *ambientLightSlider;
    QSlider                     *diffuseLightSlider;
    QSlider                     *specularLightSlider;
    QSlider                     *specularExponentSlider;

    // labels for sliders & arcballs
    QLabel                      *modelRotatorLabel;
    QLabel                      *lightRotatorLabel;
    QLabel                      *yTranslateLabel;
    QLabel                      *zTranslateLabel;
    QLabel                      *emissiveLightLabel;
    QLabel                      *diffuseLightLabel;
    QLabel                      *specularExponentLabel;
    QLabel                      *ambientLightLabel;
    QLabel                      *specularLightLabel;

    public:
    // constructor
    RenderWindow
        (
        // the object to be rendered
        TexturedObject          *newTexturedObject, 
        // the model object storing render parameters
        RenderParameters        *newRenderParameters,
        // the title for the window (with default value)
        const char              *windowName = "Object Renderer"
        );  
    
    // routine to reset the interface
    // sets every visual control to match the model
    // gets called by the controller after each change in the model
    void ResetInterface();

    // declare the render controller class a friend so it can access the UI elements
    friend class RenderController;

    }; // class RenderWindow

// end of include guard
#endif
