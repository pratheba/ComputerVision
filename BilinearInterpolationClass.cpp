#include "bilinearinterpolationClass.h"

BilinearInterpolationClass* BilinearInterpolationClass::bilinearInterpolation = NULL;


BilinearInterpolationClass::BilinearInterpolationClass(): image(NULL), width(0)
{
}

BilinearInterpolationClass::BilinearInterpolationClass(double* image_, int width_):image(image_), width(width_) {
     pixelWeight = new double[4];
     position = new Position[4];
}

BilinearInterpolationClass::~BilinearInterpolationClass() {
    if(image != NULL) { delete image; image = NULL; }
    if(position != NULL) { delete position; position = NULL; }
    if(pixelWeight != NULL) { delete pixelWeight; pixelWeight = NULL; }
}

BilinearInterpolationClass* BilinearInterpolationClass::getInstance(double* image_, int width_) {
    if(bilinearInterpolation == NULL)
        bilinearInterpolation = new BilinearInterpolationClass(image_, width_);
    return bilinearInterpolation;
}

void BilinearInterpolationClass::release() {
    if(bilinearInterpolation != NULL) {
    delete bilinearInterpolation;
    bilinearInterpolation = NULL;
    }
}

int BilinearInterpolationClass::GetPixelPositionForEdges(int currentPixel, int edgesize) {
    if(currentPixel < 0)
        return std::abs(currentPixel);
    else if(currentPixel > (edgesize))
        return std::abs(edgesize-(currentPixel - edgesize));
    else
        return edgesize;
}

void BilinearInterpolationClass::GetPixelWeight(const Position& pixelPosition, double* pixelWeight) {

    int edgePixelCol = pixelPosition.YPixelPos;
    int edgePixelRow = pixelPosition.XPixelPos;

    if((edgePixelCol) < 0 || (edgePixelCol)> (width-1))
        edgePixelCol = GetPixelPositionForEdges(edgePixelCol,width-1);
   // if((edgePixelRow) < 0 || (edgePixelRow)> (image->height()-1))
     //   edgePixelRow = GetPixelPositionForEdges(edgePixelRow,image->height()-1);

    *pixelWeight =  image[edgePixelRow* width + edgePixelCol];
}


double BilinearInterpolationClass::GetIntermediatePixelValueByInterpolation(double influenceOfPixel1, double influenceOfPixel2,
                                              double pixel1Weight, double pixel2Weight) {

    double pixelIntensity = influenceOfPixel1 * pixel2Weight + influenceOfPixel2 * pixel1Weight;
    return pixelIntensity;
}


double BilinearInterpolationClass::GetPixelValueByInterpolation() {
    int colPixel1 = (int) (floor(colPixel));
    int rowPixel1 = (int) (floor(rowPixel));

    double influenceOfPixelX1 =  (double)(rowPixel) - (double)(rowPixel1);
    double influenceOfPixelY1 =  (double)(colPixel) - (double)(colPixel1);
    double influenceOfPixelX2 =  (double)1 - influenceOfPixelX1;
    double influenceOfPixelY2 =  (double)1 - influenceOfPixelY1;

    double XY1 = GetIntermediatePixelValueByInterpolation(influenceOfPixelX1, influenceOfPixelX2, pixelWeight[0], pixelWeight[2]);
    double XY2 = GetIntermediatePixelValueByInterpolation(influenceOfPixelX1, influenceOfPixelX2, pixelWeight[1], pixelWeight[3]);

    double XY =  GetIntermediatePixelValueByInterpolation(influenceOfPixelY1, influenceOfPixelY2, XY1, XY2);

    return XY;

}

void BilinearInterpolationClass::GetNeighbouringPixelPositionAndWeight() {
    int colPixel1 = (int) ((colPixel));
    int rowPixel1 = (int) ((rowPixel));

    position[0] = Position(rowPixel1,colPixel1);
    GetPixelWeight(Position(rowPixel1,colPixel1), &pixelWeight[0]);

    position[1] = Position(rowPixel1,colPixel1+1);
    GetPixelWeight(Position(rowPixel1,colPixel1+1), &pixelWeight[1]);

    position[2] = Position(rowPixel1+1,colPixel1);
    GetPixelWeight(Position(rowPixel1+1,colPixel1), &pixelWeight[2]);

    position[3] = Position(rowPixel1+1,colPixel1+1);
    GetPixelWeight(Position(rowPixel1+1,colPixel1+1), &pixelWeight[3]);

}


double BilinearInterpolationClass::GetBilinearInterpolatedPixelValue(double colPixel_, double rowPixel_) { // X - column , Y - row

    colPixel    =   colPixel_;
    rowPixel    =   rowPixel_;

    GetNeighbouringPixelPositionAndWeight();
    return(GetPixelValueByInterpolation());
}
