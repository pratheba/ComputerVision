#ifndef BESTCLASSIFIER_H
#define BESTCLASSIFIER_H

#include <iostream>
#include <map>
#include <math.h>
#include <QDebug>
#include <cmath>

struct ClassifierParameter {
    int polarity;
    double weight;
    double threshold;
    double error;

    ClassifierParameter():polarity(-1), weight(100000.0), threshold(-999999.0), error(999999.0){}
    ClassifierParameter(const ClassifierParameter& classifierParameter_) : polarity(classifierParameter_.polarity),
        weight(classifierParameter_.weight), threshold(classifierParameter_.threshold), error(classifierParameter_.error){}
    ClassifierParameter(int polarity_, double weight_, double threshold_, double error_): polarity(polarity_),
        weight(weight_), threshold(threshold_), error(error_){}

};

class BestClassifier
{
public:
    static BestClassifier* getInstance(int* trainingLabel_, int numTrainingExamples_);
    void        Initialize(int *featureSortIdx, double *features,double* dataWeights);
    void        release();
    ClassifierParameter FindOptimalClassifierParameter();
    ClassifierParameter classifierParameter;



private:
    void        UpdateWeightsbeforeSample();
    void        FindErrorForTheta();
    void        calculatePostiveNegative();


    BestClassifier();
     BestClassifier(double* dataWeights_);
    static BestClassifier* bestClassifier;
    ~BestClassifier();



#pragma region membervariables
    double              totalPositiveWeight;
    double              totalNegativeWeight;
    double              positiveWeightBeforeSample;
    double              negativeWeightBeforeSample;
    int                 currentIndex;
    static int          *trainingLabel;
    static int          numTrainingExamples;
    double              *dataWeights;
    int                 *featureSortIdx;
    double              *features;

#pragma endregion membervariables
};

#endif // BESTCLASSIFIER_H
