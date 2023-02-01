//////////////////////////////////////////////////////////////////////
//
//  University of Leeds
//  COMP 5812M Foundations of Modelling & Rendering
//  User Interface for Coursework
//
//  September, 2022
//
//  ------------------------
//  RGBAValue.h
//  ------------------------
//
//  A minimal class for a colour in RGB space
//
///////////////////////////////////////////////////

#ifndef RGBAVALUEF_H
#define RGBAVALUEF_H

#include <iostream>
#include <cstddef>

class RGBAValueF
    {
    public:
    // just a container for the components
    float red, green, blue, alpha;

    // default constructor
    RGBAValueF();

    // value constructor with floats
    // values outside 0.0-255.0 are clamped
    RGBAValueF(float Red, float Green, float Blue, float Alpha = 1.0f);

    // copy constructor
    RGBAValueF(const RGBAValueF &other);

    // colour modulation routine:
    RGBAValueF modulate(const RGBAValueF &right) const;

    }; // RGBAValue


// convenience routines for scalar multiplication and addition
RGBAValueF operator *(float scalar, const RGBAValueF &colour);
RGBAValueF operator +(const RGBAValueF &left, const RGBAValueF &right);
bool operator ==(const RGBAValueF &left, const RGBAValueF &right);

// stream input
std::istream & operator >> (std::istream &inStream, RGBAValueF &value);

// stream output
std::ostream & operator << (std::ostream &outStream, const RGBAValueF &value);

#endif
