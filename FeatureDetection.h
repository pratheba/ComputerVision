#ifndef FEATUREDETECTION_H
#define FEATUREDETECTION_H

#include <QMainWindow>
#include<QDebug>

class FeatureDetection
{
public:
    FeatureDetection();
    ~FeatureDetection();
    void GetAveragePatch(double *trainingData, int *trainingLabel, int numTrainingExamples, int patchSize);
    double* GetAverageFace() { return averageFace; }
    double* GetAverageBackground() { return averageBackground; }
    void    DisplayAverageImage(QImage* image);

private:
    double  *averageFace;
    double  *averageBackground;
    double  *averagePatch;
    int     numberOfFaces;
    int     numberOfBackground;
    int     numberOfPixels;
    int     patchSize;
    std::vector<int>facesIndex;
    std::vector<int>backgroundIndex;

    void InitializeAveragePatches(int patchSize);
    void AddPatches(double *trainingData,int dataIndex);
    void AddPatches(double *trainingData, int *trainingLabel, int dataIndex);
    void NormalizeData();
};


#endif // FEATUREDETECTION_H
