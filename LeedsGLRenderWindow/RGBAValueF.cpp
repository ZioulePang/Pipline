//////////////////////////////////////////////////////////////////////
//
//  University of Leeds
//  COMP 5812M Foundations of Modelling & Rendering
//  User Interface for Coursework
//
//  September, 2022
//
//  ------------------------
//  RGBAValueF.h
//  ------------------------
//
//  A minimal class for a colour in RGB space
//
///////////////////////////////////////////////////

#include "RGBAValueF.h"

using namespace std;

// default constructor
RGBAValueF::RGBAValueF()
    :
    red{0},
    green{0},
    blue{0},
    alpha{0}
    { // constructor
    } // constructor



// value constructor with floats
// values outside 0.0-1.0 are clamped
RGBAValueF::RGBAValueF(float Red, float Green, float Blue, float Alpha)
    { // constructor
    if (Red > 1.0f) Red = 1.0f;
    if (Red < 0.0f) Red = 0.0f;

    red = Red;

    if (Green > 1.0f) Green = 1.0f;
    if (Green < 0.0f) Green = 0.0f;
    green = Green;

    if (Blue > 1.0f) Blue = 1.0f;
    if (Blue < 0.0f) Blue = 0.0f;
    blue = Blue;

    if (Alpha > 1.0f) Alpha = 1.0f;
    if (Alpha < 0.0f) Alpha = 0.0f;
    alpha = Alpha;
    } // constructor

// copy constructor
RGBAValueF::RGBAValueF(const RGBAValueF &other)
    :
    red(other.red),
    green(other.green),
    blue(other.blue),
    alpha(other.alpha)
    { // copy constructor
    } // copy constructor

// convenience routines for scalar multiplication and addition
RGBAValueF operator *(float scalar, const RGBAValueF &colour)
    { // operator * ()
    // feed the new values into a constructor & return
    return RGBAValueF(  scalar * colour.red,
                        scalar * colour.green,
                        scalar * colour.blue,
                        scalar * colour.alpha);
    } // operator * ()

RGBAValueF operator +(const RGBAValueF &left, const RGBAValueF &right)
    { // operator + ()
    // compute new values, then use constructor to scale & clamp
    return RGBAValueF(  (left.red   + right.red),
                        (left.green + right.green),
                        (left.blue  + right.blue),
                        (left.alpha + right.alpha));
    } // operator + ()

bool operator ==(const RGBAValueF &left, const RGBAValueF &right)
    { // operator  ()
    return    abs(left.red - right.red) < std::numeric_limits<float>::epsilon() &&
              abs(left.green - right.green) < std::numeric_limits<float>::epsilon() &&
              abs(left.blue - right.blue) < std::numeric_limits<float>::epsilon()  &&
              abs(left.alpha - right.alpha) < std::numeric_limits<float>::epsilon();
    } // operator  ()

// colour modulation routine:
// NB: this routine scales each component by 1/255.0, multiplies then inverts
RGBAValueF RGBAValueF::modulate(const RGBAValueF &right) const
    {
    return RGBAValueF(red * right.red, green * right.green, blue * right.blue, alpha * right.alpha);
    } // modulate()

// stream input
std::istream & operator >> (std::istream &inStream, RGBAValueF &value)
    { // stream input
    float component;
    inStream >> value.red;
    inStream >> value.green;
    inStream >> value.blue;
    inStream >> value.alpha;
    return inStream;
    } // stream input

// stream output
// WARNING - does not output alpha
std::ostream & operator << (std::ostream &outStream, const RGBAValueF &value)
    { // stream output
    outStream << value.red << " ";
    outStream << value.green<< " ";
    outStream << value.blue << " ";
    outStream << value.alpha;
    return outStream;

} // stream output
