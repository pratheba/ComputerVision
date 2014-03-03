#ifndef BILINEARINTERPOLATIONCLASS_H
#define BILINEARINTERPOLATIONCLASS_H

#include <iostream>
#include "math.h"
#include <cmath>
#include <QDebug>

struct Position {
    int XPixelPos;
    int YPixelPos;

    Position():XPixelPos(-1), YPixelPos(-1){}
    Position(int X, int Y):XPixelPos(X), YPixelPos(Y){}
    Position(const Position& position_):XPixelPos(position_.XPixelPos), YPixelPos(position_.YPixelPos){}
};

class BilinearInterpolationClass
{
public:
    BilinearInterpolationClass (const BilinearInterpolationClass&) = delete;
    BilinearInterpolationClass& operator=(const BilinearInterpolationClass&) = delete;

    static BilinearInterpolationClass* getInstance(double* image_, int width_);

    void release();

    double GetBilinearInterpolatedPixelValue(double colPixel, double rowPixel);

private:
    BilinearInterpolationClass();
    BilinearInterpolationClass(double *image_, int width_);
    static BilinearInterpolationClass* bilinearInterpolation;
    ~BilinearInterpolationClass();

    double*     image;
    int         width;
    double      colPixel, rowPixel;
    double*     pixelWeight;
    Position*   position;

    void    GetNeighbouringPixelPositionAndWeight();
    double  GetPixelValueByInterpolation();
    double  GetIntermediatePixelValueByInterpolation(double influenceOfPixel1, double influenceOfPixel2,
                                                  double pixel1Weight, double pixel2Weight);
    void    GetPixelWeight(const Position& pixelPosition, double* pixelWeight);
    int     GetPixelPositionForEdges(int currentPixel, int edgesize);

};

#endif // BILINEARINTERPOLATIONCLASS_H
