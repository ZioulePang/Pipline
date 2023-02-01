/////////////////////////////////////////////////////////////////
//
//  University of Leeds
//  COMP 5812M Foundations of Modelling & Rendering
//  User Interface for Coursework
//
//  September, 2020
//
//  -----------------------------
//  Render Parameters
//  -----------------------------
//  
//  This is part of the "model" in the MVC paradigm
//  We separate out the render parameters from the object rendered
//
/////////////////////////////////////////////////////////////////

// include guard
#ifndef _RENDER_PARS_H
#define _RENDER_PARS_H

#include "Matrix4.h"

// class for the render parameters
class RenderParameters
    { // class RenderParameters
    public:
    
    // we have a widget with an arcball which stores the rotation
    // we'll be lazy and leave it there instead of moving it here
    
    // we store x & y translations
    float xTranslate, yTranslate;

    // and a zoom scale
    float zTranslate;
    
    // we have the position of the light
    float lightPosition[4];
    
    // we will want two homogeneous matrices holding the rotations for the
    // model and the light.
    Matrix4 rotationMatrix;
    Matrix4 lightMatrix;
    
    // and the various lighting parameters
    float emissiveMaterial;
    float ambientMaterial;
    float diffuseMaterial;
    float specularMaterial;
    float specularExponent;
    
    // and the booleans
    bool useLighting;
    bool texturedRendering;
    bool textureModulation;
    bool depthTestOn;

    bool showAxes;
    bool showObject;
    bool centreObject;
    bool perspective;
    bool scaleObject;
    bool mapUVWToRGB;

    // constructor
    RenderParameters()
        :
        xTranslate(0.0f),
        yTranslate(0.0f),
        zTranslate(0.0f),
        emissiveMaterial(0.0f),
        ambientMaterial(0.2f),
        diffuseMaterial(0.6f),
        specularMaterial(0.3f),
        specularExponent(4.0),
        useLighting(false),
        texturedRendering(false),
        textureModulation(false),
        depthTestOn(false),
        showAxes(false),
        showObject(false),
        centreObject(false),
        perspective(true),
        scaleObject(false),
        mapUVWToRGB(false)
        { // constructor
        
        // start the lighting at the viewer's direction
        lightPosition[0] = 0.0f;
        lightPosition[1] = 0.0f;
        lightPosition[2] = 1.0f;
        lightPosition[3] = 0.0f;

        // because we are paranoid, we will initialise the matrices to the identity
        rotationMatrix.SetIdentity();
        lightMatrix.SetIdentity();
        } // constructor

    // accessor for scaledXTranslate

    }; // class RenderParameters

// now define some macros for bounds on parameters
#define TRANSLATE_MIN -1.0f
#define TRANSLATE_MAX 1.0f


#define LIGHTING_MIN 0.0f
#define LIGHTING_MAX 1.0f

#define SPECULAR_EXPONENT_LOG_MIN -2.0f
#define SPECULAR_EXPONENT_LOG_MAX 2.0f
#define SPECULAR_EXPONENT_MIN 0.01f
#define SPECULAR_EXPONENT_MAX 100.0f

// this is to scale to/from integer values
#define PARAMETER_SCALING 100


// end of include guard
#endif
