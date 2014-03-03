#include "FeatureDetection.h"
#include <iostream>

FeatureDetection::FeatureDetection():averageFace(NULL), averageBackground(NULL), averagePatch(NULL),
    numberOfFaces(0), numberOfBackground(0), numberOfPixels(0), patchSize(0)
{
}



void FeatureDetection::InitializeAveragePatches(int patchsize)
{
    patchSize           =   patchsize;
    numberOfPixels      =   patchsize * patchsize;
    averageFace         =   new double[numberOfPixels];
    averageBackground   =   new double[numberOfPixels];

    for(int pixel = 0; pixel < numberOfPixels; pixel++) {
        averageFace[pixel]          =   0.0;
        averageBackground[pixel]    =   0.0;
    }

}

void FeatureDetection::AddPatches(double *trainingData,int dataIndex) {

    for(int rowPixel = 0; rowPixel < patchSize; rowPixel++)
        for(int colPixel = 0; colPixel < patchSize; colPixel++)
            averagePatch[rowPixel * patchSize + colPixel] +=
                    trainingData[dataIndex * numberOfPixels + rowPixel*patchSize + colPixel];
}

void FeatureDetection::AddPatches(double *trainingData, int *trainingLabel, int dataIndex) {
    if(trainingLabel[dataIndex]==1) {
        averagePatch    =   averageFace;
        numberOfFaces++;
        facesIndex.push_back(dataIndex);
    }
    else {
        averagePatch    =   averageBackground;
        numberOfBackground++;
        backgroundIndex.push_back(dataIndex);
    }

    AddPatches(trainingData, dataIndex);
}

void FeatureDetection::NormalizeData() {
    for(int pixel = 0; pixel < numberOfPixels; pixel++) {
        averageFace[pixel]          /=  (double)numberOfFaces;
        averageBackground[pixel]    /=  (double)numberOfBackground;
    }
}

void FeatureDetection::GetAveragePatch(double *trainingData, int *trainingLabel, int numTrainingExamples, int patchSize) {

    InitializeAveragePatches(patchSize);

    for(int dataIndex = 0 ; dataIndex < numTrainingExamples; dataIndex ++)
            AddPatches(trainingData, trainingLabel, dataIndex);

    NormalizeData();
}

void FeatureDetection::DisplayAverageImage(QImage* displayImage) {

    int numberOfPatchColumn =   displayImage->width() / patchSize;
    int numberOfPatchRow    =   displayImage->height() / patchSize;



    std::vector<int>::iterator iterBegin    =   facesIndex.begin();
    std::vector<int>::iterator iterEnd      =   facesIndex.end();

    qDebug() << displayImage->offset() << numberOfPatchColumn << numberOfPatchRow <<"::"<< facesIndex.size();

    //for( ; iterBegin != iterEnd; iterBegin++) {
       // int rowindex    =   (iterBegin / numberOfPatchColumn) - 1;
        //int colindex    =   (iterBegin % numberOfPatchColumn)

        for(int rowPixel = 0; rowPixel < patchSize; rowPixel++)
            for(int colPixel = 0; colPixel < patchSize; colPixel++) {
                double color = averageFace[rowPixel * patchSize + colPixel];
                displayImage->setPixel(colPixel, rowPixel, qRgb(color, color, color));
            }

    //*displayImage = tempimage->copy(0,0,patchSize, patchSize);
}
