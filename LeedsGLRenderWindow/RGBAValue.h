//////////////////////////////////////////////////////////////////////
//
//  University of Leeds
//  COMP 5812M Foundations of Modelling & Rendering
//  User Interface for Coursework
//
//  September, 2020
//
//  ------------------------
//  RGBAValue.h
//  ------------------------
//  
//  A minimal class for a colour in RGB space
//  
///////////////////////////////////////////////////

#ifndef RGBAVALUE_H
#define RGBAVALUE_H

#include <iostream>
#include <cstddef>
#include "RGBAValueF.h"

class RGBAValue
    { // RGBAValue
    public:
    // just a container for the components
    std::byte red, green, blue, alpha;

    // default constructor
    RGBAValue();

    // value constructor with default alpha
    RGBAValue(std::byte Red, std::byte Green, std::byte Blue, std::byte Alpha = std::byte{0xff});


    // value constructor with floats
    // values outside 0.0-255.0 are clamped
    RGBAValue(float Red, float Green, float Blue, float Alpha = 1.0f);

    // copy constructor
    RGBAValue(const RGBAValue &other);

    // copy constructor from float
    RGBAValue(const RGBAValueF &other);

    
    // colour modulation routine:
    // NB: this routine scales each component by 1/255.0, multiplies then inverts
    RGBAValue modulate(const RGBAValue &right) const;
    
    }; // RGBAValue

// convenience routines for scalar multiplication and addition
RGBAValue operator *(float scalar, const RGBAValue &colour);
RGBAValue operator +(const RGBAValue &left, const RGBAValue &right);
bool operator ==(const RGBAValue &left, const RGBAValue &right);

// stream input
// WARNING - always sets alpha to 255
std::istream & operator >> (std::istream &inStream, RGBAValue &value);

// stream output
// WARNING - does not output alpha
std::ostream & operator << (std::ostream &outStream, const RGBAValue &value);

#endif
