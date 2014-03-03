#include "bestclassifier.h"

BestClassifier* BestClassifier::bestClassifier = NULL;
int* BestClassifier::trainingLabel =    NULL;
int BestClassifier::numTrainingExamples =   0;

BestClassifier::BestClassifier(): dataWeights(NULL),featureSortIdx(NULL),features(NULL) {
}

 BestClassifier::~BestClassifier() {
     if(dataWeights != NULL) { delete dataWeights; dataWeights = NULL; }
     if(featureSortIdx != NULL) { delete featureSortIdx; featureSortIdx = NULL; }
     if(features != NULL) { delete features; features = NULL; }
 }

 void BestClassifier::release() {
     delete trainingLabel;
     delete bestClassifier;
     bestClassifier = NULL;
 }

BestClassifier* BestClassifier::getInstance(int* trainingLabel_, int numTrainingExamples_) {
    if(bestClassifier == NULL) {
        bestClassifier      =   new BestClassifier();
        trainingLabel       =   trainingLabel_;
        numTrainingExamples =   numTrainingExamples_;
    }
    return bestClassifier;
}

void BestClassifier::Initialize(int *featureSortIdx_, double *features_, double* dataWeights_) {
    featureSortIdx  =   featureSortIdx_;
    features        =   features_;
    dataWeights     =   dataWeights_;

    totalPositiveWeight = 0;
    totalNegativeWeight = 0;
    positiveWeightBeforeSample = 0;
    negativeWeightBeforeSample = 0;

    classifierParameter.error = 9999999.0;
    classifierParameter.polarity = -1;
    classifierParameter.threshold = -9999999.0;
    classifierParameter.weight = 1000000.0;

    for(int dataIndex = 0; dataIndex < numTrainingExamples; dataIndex++) {
        if(trainingLabel[dataIndex] == 0)
            totalNegativeWeight += dataWeights[dataIndex];
        else
            totalPositiveWeight += dataWeights[dataIndex];
    }

}

void BestClassifier::UpdateWeightsbeforeSample() {

    if(trainingLabel[currentIndex] == 0)
        negativeWeightBeforeSample += dataWeights[currentIndex];
    else
        positiveWeightBeforeSample += dataWeights[currentIndex];
}


void BestClassifier::FindErrorForTheta() {
    UpdateWeightsbeforeSample();

    double  FaceAboveThresholdError     =   positiveWeightBeforeSample + ( totalNegativeWeight - negativeWeightBeforeSample);
    double  FaceBelowThresholdError     =   negativeWeightBeforeSample + ( totalPositiveWeight - positiveWeightBeforeSample);
    double  currentError                =   std::min(FaceAboveThresholdError, FaceBelowThresholdError);

    if( currentError <= classifierParameter.error ) {
            classifierParameter.error   =   currentError;
        if(FaceAboveThresholdError < FaceBelowThresholdError)
            classifierParameter.polarity = 1;
        else
            classifierParameter.polarity = 0;

        classifierParameter.weight     =    log((1- currentError)/(currentError));
        classifierParameter.threshold  =   (features[currentIndex]);// + features[featureSortIdx[currentIndex+1]])/(double)2
    }

}

ClassifierParameter BestClassifier::FindOptimalClassifierParameter() {

    for(int dataIndex = 0; dataIndex < numTrainingExamples; dataIndex++) {
        currentIndex    =       featureSortIdx[dataIndex];
        FindErrorForTheta();
    }

    return classifierParameter;
}
