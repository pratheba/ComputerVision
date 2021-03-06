
#include "mainwindow.h"
#include "math.h"
#include "ui_mainwindow.h"
#include <QtGui>

#include <iostream>
#include "FeatureDetection.h"
#include "bilinearinterpolationClass.h"
#include "bestclassifier.h"



BilinearInterpolationClass* bilinearInterpolation = NULL;
BestClassifier*             bestClassifierClass = NULL;

FILE *featureValue = fopen("featureValue.txt", "w");

/*******************************************************************************
    The following are helper routines with code already written.
    The routines you'll need to write for the assignment are below.
*******************************************************************************/

/*******************************************************************************
    Open the training dataset
        posdirectory - Directory containing face images
        negdirectory - Directory containing non-face images
        trainingData - Array used to store the data
        trainingLabel - Label assigned to training data (1 = face, 0 = non-face)
        numTrainingExamples - Number of training examples
        patchSize - Size of training patches
*******************************************************************************/
void MainWindow::OpenDataSet(QDir posdirectory, QDir negdirectory, double *trainingData, int *trainingLabel, int numTrainingExamples, int patchSize)
{
    int i, c, r;
    QStringList imgNames;
    QImage inImage;
    QRgb pixel;

    imgNames = posdirectory.entryList();

    int idx = 0;

    for(i=0;i<imgNames.length();i++)
        if(idx < numTrainingExamples/2)
    {
        // use "\\" for windows machine
        inImage.load(posdirectory.absolutePath() + "/" + imgNames.at(i));

        if(!(inImage.isNull()))
        {
            for(r=0;r<patchSize;r++)
                for(c=0;c<patchSize;c++)
                {
                    pixel = inImage.pixel(c, r);
                    trainingData[idx*patchSize*patchSize + r*patchSize + c] = (double) qGreen(pixel);
                }

            trainingLabel[idx] = 1;

            idx++;
        }
    }

    imgNames = negdirectory.entryList();

    for(i=0;i<imgNames.length();i++)
        if(idx < numTrainingExamples)
    {
        // use "\\" for windows machine
        inImage.load(negdirectory.absolutePath() + "/" + imgNames.at(i));

        if(!(inImage.isNull()))
        {
            for(r=0;r<patchSize;r++)
                for(c=0;c<patchSize;c++)
                {
                    pixel = inImage.pixel(c, r);
                    trainingData[idx*patchSize*patchSize + r*patchSize + c] = (double) qGreen(pixel);
                }

            trainingLabel[idx] = 0;

            idx++;
        }
    }
}

/*******************************************************************************
    DisplayTrainingDataset - Display example patches from training dataset
        displayImage - Display image
        trainingData - Array used to store the data
        trainingLabel - Label assigned to training data (1 = face, 0 = non-face)
        numTrainingExamples - Number of training examples
        patchSize - Size of training patches
*******************************************************************************/
void MainWindow::DisplayTrainingDataset(QImage *displayImage, double *trainingData, int *trainingLabel, int numTrainingExamples, int patchSize)
{
    int w = displayImage->width();
    int h = displayImage->height();
    int r, c;
    int rOffset = 0;
    int cOffset = 0;
    bool inBounds = true;
    int ct = 0;

    while(inBounds)
    {
        int idx = rand()%numTrainingExamples;

        for(r=0;r<patchSize;r++)
            for(c=0;c<patchSize;c++)
            {
                if(trainingLabel[idx] == 1)
                {
                    int val = (int) trainingData[idx*patchSize*patchSize + r*patchSize + c];
                    displayImage->setPixel(c + cOffset, r + rOffset, qRgb(val, val, val));

                }
                else
                {
                    int val = (int) trainingData[idx*patchSize*patchSize + r*patchSize + c];
                    displayImage->setPixel(c + cOffset, r + rOffset, qRgb(val, val, val));
                }
            }

        cOffset += patchSize;

        if(cOffset + patchSize >= w)
        {
            cOffset = 0;
            rOffset += patchSize;

            if(rOffset + patchSize >= h)
                inBounds = false;
        }

        ct++;
    }
}

/*******************************************************************************
    SaveClassifier - Save the computed AdaBoost classifier
        fileName - Name of file
*******************************************************************************/
void MainWindow::SaveClassifier(QString fileName)
{
   int i, j;
   FILE *out;

   out = fopen(fileName.toLatin1(), "w");

   fprintf(out, "%d\n", m_NumWeakClassifiers);

   for(i=0;i<m_NumWeakClassifiers;i++)
   {
       fprintf(out, "%d\n", m_WeakClassifiers[i].m_NumBoxes);

       for(j=0;j<m_WeakClassifiers[i].m_NumBoxes;j++)
           fprintf(out, "%lf\t%lf\t%lf\t%lf\t%lf\n", m_WeakClassifiers[i].m_BoxSign[j], m_WeakClassifiers[i].m_Box[j][0][0], m_WeakClassifiers[i].m_Box[j][0][1],
                   m_WeakClassifiers[i].m_Box[j][1][0], m_WeakClassifiers[i].m_Box[j][1][1]);

       fprintf(out, "%lf\n", m_WeakClassifiers[i].m_Area);
       fprintf(out, "%lf\n", m_WeakClassifiers[i].m_Polarity);
       fprintf(out, "%lf\n", m_WeakClassifiers[i].m_Threshold);
       fprintf(out, "%lf\n", m_WeakClassifiers[i].m_Weight);
   }

   fclose(out);
}

/*******************************************************************************
    OpenClassifier - Open the computed AdaBoost classifier
        fileName - Name of file
*******************************************************************************/
void MainWindow::OpenClassifier(QString fileName)
{
    int i, j;
    FILE *in;

    in = fopen(fileName.toLatin1(), "r");

    fscanf(in, "%d\n", &m_NumWeakClassifiers);
    m_WeakClassifiers = new CWeakClassifiers [m_NumWeakClassifiers];

    for(i=0;i<m_NumWeakClassifiers;i++)
    {
        fscanf(in, "%d\n", &(m_WeakClassifiers[i].m_NumBoxes));
        m_WeakClassifiers[i].m_Box = new double [m_WeakClassifiers[i].m_NumBoxes][2][2];
        m_WeakClassifiers[i].m_BoxSign = new double [m_WeakClassifiers[i].m_NumBoxes];

        for(j=0;j<m_WeakClassifiers[i].m_NumBoxes;j++)
            fscanf(in, "%lf\t%lf\t%lf\t%lf\t%lf\n", &(m_WeakClassifiers[i].m_BoxSign[j]), &(m_WeakClassifiers[i].m_Box[j][0][0]), &(m_WeakClassifiers[i].m_Box[j][0][1]),
                    &(m_WeakClassifiers[i].m_Box[j][1][0]), &(m_WeakClassifiers[i].m_Box[j][1][1]));

        fscanf(in, "%lf\n", &(m_WeakClassifiers[i].m_Area));
        fscanf(in, "%lf\n", &(m_WeakClassifiers[i].m_Polarity));
        fscanf(in, "%lf\n", &(m_WeakClassifiers[i].m_Threshold));
        fscanf(in, "%lf\n", &(m_WeakClassifiers[i].m_Weight));
    }

    fclose(in);

}

/*******************************************************************************
    DisplayClassifiers - Display the Haar wavelets for the classifier
        displayImage - Display image
        weakClassifiers - The weak classifiers used in AdaBoost
        numWeakClassifiers - Number of weak classifiers
*******************************************************************************/
void MainWindow::DisplayClassifiers(QImage *displayImage, CWeakClassifiers *weakClassifiers, int numWeakClassifiers)
{
    int w = displayImage->width();
    int h = displayImage->height();
    int i, j, r, c;
    int rOffset = 0;
    int cOffset = 0;
    int size = 50;
    bool inBounds = true;

    displayImage->fill(qRgb(0,0,0));

    for(i=0;i<numWeakClassifiers & inBounds;i++)
    {
        for(r=0;r<size;r++)
            for(c=0;c<size;c++)
            {
                 displayImage->setPixel(c + cOffset, r + rOffset, qRgb(128, 128, 128));
            }

        for(j=0;j<weakClassifiers[i].m_NumBoxes;j++)
            for(r=(int) ((double) size*weakClassifiers[i].m_Box[j][0][1]);r<(int) ((double) size*weakClassifiers[i].m_Box[j][1][1]);r++)
                for(c=(int) ((double) size*weakClassifiers[i].m_Box[j][0][0]);c<(int) ((double) size*weakClassifiers[i].m_Box[j][1][0]);c++)
                {
                    if(weakClassifiers[i].m_BoxSign[j] > 0.0)
                        displayImage->setPixel(c + cOffset, r + rOffset, qRgb(255, 255, 255));
                    else
                        displayImage->setPixel(c + cOffset, r + rOffset, qRgb(0, 0, 0));
                }

        cOffset += size+1;

        if(cOffset + size >= w)
        {
            cOffset = 0;
            rOffset += size + 1;

            if(rOffset + size >= h)
                inBounds = false;
        }
    }
}

/*******************************************************************************
    DisplayIntegralImage - Display the integral image
        displayImage - Display image
        integralImage - Output integral image
        w, h - Width and height of image
*******************************************************************************/
void MainWindow::DisplayIntegralImage(QImage *displayImage, double *integralImage, int w, int h)
{
    int r, c;
    double maxVal = integralImage[(h-1)*w + w-1];

    for(r=0;r<h;r++)
        for(c=0;c<w;c++)
        {
            int val = (int) (255.0*integralImage[r*w + c]/maxVal);

            displayImage->setPixel(c, r, qRgb(val, val, val));
        }
}

/*******************************************************************************
    InitializeFeatures - Randomly initialize the candidate weak classifiers
        weakClassifiers - Candidate weak classifiers
        numWeakClassifiers - Number of candidate weak classifiers
*******************************************************************************/
void MainWindow::InitializeFeatures(CWeakClassifiers *weakClassifiers, int numWeakClassifiers)
{
    int i;

    for(i=0;i<numWeakClassifiers;i++)
    {
        double x, y, w, h;

        // We don't know these values yet, so just initialize to 0
        weakClassifiers[i].m_Polarity = 0.0;
        weakClassifiers[i].m_Threshold = 0.0;
        weakClassifiers[i].m_Weight = 0.0;

        // The Haar wavelet's corners can range in the area of 0.02 to 0.98, with a minimum size of 0.25
        // We limit the range to [0.2, 0.98], instead of [0, 1] so we don't need to worry about checking
        // out of bounds errors later on, i.e. in the BilinearInterpolation function.

        // x position of box and width
        w = 0.25 + 0.71*(double) rand()/(double) RAND_MAX;
        x = 0.02 + (0.96 - w)*(double) rand()/(double) RAND_MAX;

        // y position of box and height
        h = 0.25 + 0.71*(double) rand()/(double) RAND_MAX;
        y = 0.02 + (0.96 - h)*(double) rand()/(double) RAND_MAX;

        int boxType = rand()%4;

        if(boxType == 0)
        {
            // Vertical boxes
            weakClassifiers[i].m_NumBoxes = 2;
            weakClassifiers[i].m_Box = new double [weakClassifiers[i].m_NumBoxes][2][2];
            weakClassifiers[i].m_BoxSign = new double [weakClassifiers[i].m_NumBoxes];

            weakClassifiers[i].m_BoxSign[0] = 1.0;
            weakClassifiers[i].m_Box[0][0][0] = x;
            weakClassifiers[i].m_Box[0][0][1] = y;
            weakClassifiers[i].m_Box[0][1][0] = x + w/2;
            weakClassifiers[i].m_Box[0][1][1] = y + h;

            weakClassifiers[i].m_BoxSign[1] = -1.0;
            weakClassifiers[i].m_Box[1][0][0] = x + w/2;
            weakClassifiers[i].m_Box[1][0][1] = y;
            weakClassifiers[i].m_Box[1][1][0] = x + w;
            weakClassifiers[i].m_Box[1][1][1] = y + h;
        }

        if(boxType == 1)
        {
            // 2 Horizontal boxes
            weakClassifiers[i].m_NumBoxes = 2;
            weakClassifiers[i].m_Box = new double [weakClassifiers[i].m_NumBoxes][2][2];
            weakClassifiers[i].m_BoxSign = new double [weakClassifiers[i].m_NumBoxes];

            weakClassifiers[i].m_BoxSign[0] = 1.0;
            weakClassifiers[i].m_Box[0][0][0] = x;
            weakClassifiers[i].m_Box[0][0][1] = y;
            weakClassifiers[i].m_Box[0][1][0] = x + w;
            weakClassifiers[i].m_Box[0][1][1] = y + h/2;

            weakClassifiers[i].m_BoxSign[1] = -1.0;
            weakClassifiers[i].m_Box[1][0][0] = x;
            weakClassifiers[i].m_Box[1][0][1] = y + h/2;
            weakClassifiers[i].m_Box[1][1][0] = x + w;
            weakClassifiers[i].m_Box[1][1][1] = y + h;
        }

        if(boxType == 2)
        {
            // 3 Vertical boxes
            weakClassifiers[i].m_NumBoxes = 3;
            weakClassifiers[i].m_Box = new double [weakClassifiers[i].m_NumBoxes][2][2];
            weakClassifiers[i].m_BoxSign = new double [weakClassifiers[i].m_NumBoxes];

            weakClassifiers[i].m_BoxSign[0] = 1.0;
            weakClassifiers[i].m_Box[0][0][0] = x;
            weakClassifiers[i].m_Box[0][0][1] = y;
            weakClassifiers[i].m_Box[0][1][0] = x + w/3;
            weakClassifiers[i].m_Box[0][1][1] = y + h;

            weakClassifiers[i].m_BoxSign[1] = -2.0;
            weakClassifiers[i].m_Box[1][0][0] = x + w/3;
            weakClassifiers[i].m_Box[1][0][1] = y;
            weakClassifiers[i].m_Box[1][1][0] = x + 2*w/3;
            weakClassifiers[i].m_Box[1][1][1] = y + h;

            weakClassifiers[i].m_BoxSign[2] = 1.0;
            weakClassifiers[i].m_Box[2][0][0] = x + 2*w/3;
            weakClassifiers[i].m_Box[2][0][1] = y;
            weakClassifiers[i].m_Box[2][1][0] = x + w;
            weakClassifiers[i].m_Box[2][1][1] = y + h;
        }

        if(boxType == 4)
        {
            // 4 diagonal boxes
            weakClassifiers[i].m_NumBoxes = 4;
            weakClassifiers[i].m_Box = new double [weakClassifiers[i].m_NumBoxes][2][2];
            weakClassifiers[i].m_BoxSign = new double [weakClassifiers[i].m_NumBoxes];

            weakClassifiers[i].m_BoxSign[0] = 1.0;
            weakClassifiers[i].m_Box[0][0][0] = x;
            weakClassifiers[i].m_Box[0][0][1] = y;
            weakClassifiers[i].m_Box[0][1][0] = x + w/2;
            weakClassifiers[i].m_Box[0][1][1] = y + h/2;

            weakClassifiers[i].m_BoxSign[1] = -1.0;
            weakClassifiers[i].m_Box[1][0][0] = x + w/2;
            weakClassifiers[i].m_Box[1][0][1] = y;
            weakClassifiers[i].m_Box[1][1][0] = x + w;
            weakClassifiers[i].m_Box[1][1][1] = y + h/2;

            weakClassifiers[i].m_BoxSign[2] = 1.0;
            weakClassifiers[i].m_Box[2][0][0] = x ;
            weakClassifiers[i].m_Box[2][0][1] = y + h/2;
            weakClassifiers[i].m_Box[2][1][0] = x + w/2;
            weakClassifiers[i].m_Box[2][1][1] = y + h;

            weakClassifiers[i].m_BoxSign[3] = -1.0;
            weakClassifiers[i].m_Box[3][0][0] = x + w/2;
            weakClassifiers[i].m_Box[3][0][1] = y + h/2;
            weakClassifiers[i].m_Box[3][1][0] = x + w;
            weakClassifiers[i].m_Box[3][1][1] = y + h;
        }

        weakClassifiers[i].m_Area = w*h;
    }
}

/*******************************************************************************
    ConvertColorToDouble - Simple helper function to convert from RGB to double
        image - Input image
        dImage - Output double image
        w, h - Image width and height
*******************************************************************************/
void MainWindow::ConvertColorToDouble(QImage image, double *dImage, int w, int h)
{
    QRgb pixel;
    int r, c;

    for(r=0;r<h;r++)
        for(c=0;c<w;c++)
        {
            pixel = image.pixel(c, r);
            dImage[r*w + c] = qGreen(pixel);
        }
}

/*******************************************************************************
    ComputeTrainingSetFeatures - Compute all of the features for the training dataset
        trainingData - Array used to store the data
        features - Array holding feature values
        numTrainingExamples - Number of training examples
        patchSize - Size of training patches
        weakClassifiers - Candidate weak classifiers
        numWeakClassifiers - Number of candidate weak classifiers
*******************************************************************************/
void MainWindow::ComputeTrainingSetFeatures(double *trainingData, double *features,
                                int numTrainingExamples, int patchSize, CWeakClassifiers *weakClassifiers, int numWeakClassifiers)
{
    int i;
    double *integralImage = new double [patchSize*patchSize];

    for(i=0;i<numTrainingExamples;i++)
    {
        // Compute features for training examples

        // First compute the integral image for each patch
        IntegralImage(&(trainingData[i*patchSize*patchSize]), integralImage, patchSize, patchSize);

        // Compute the Haar wavelets
        ComputeFeatures(integralImage, 0, 0, patchSize, &(features[i*numWeakClassifiers]), weakClassifiers, numWeakClassifiers, patchSize);
    }


    // We shouldn't need the training data anymore so let's delete it.
    delete [] trainingData;

    delete [] integralImage;

    fclose(featureValue);


}

/*******************************************************************************
    DisplayFeatures - Display the computed features (green = faces, red = background)
        displayImage - Display image
        features - Array holding feature values
        trainingLabel - Label assigned to training data (1 = face, 0 = non-face)
        numFeatures - Number of features
        numTrainingExamples - Number of training examples
*******************************************************************************/
void MainWindow::DisplayFeatures(QImage *displayImage, double *features, int *trainingLabel, int numFeatures, int numTrainingExamples)
{
    int r, c;
    int w = displayImage->width();
    int h = displayImage->height();
    int posCt = 0;
    int negCt = 0;

    double mean = 0.0;
    double meanCt = 0.0;

    for(r=0;r<numTrainingExamples;r+=10)
    {
        for(c=0;c<numFeatures;c++)
        {
            mean += fabs(features[r*numFeatures + c]);
            meanCt++;
        }
    }

    mean /= meanCt;

    for(r=0;r<numTrainingExamples;r++)
    {
        if(trainingLabel[r] == 1 && posCt < h/2)
        {
            for(c=0;c<numFeatures;c++)
                if(c < w)
            {
                int val = 255.0*(features[r*numFeatures + c]/(4.0*mean)) + 128.0;
                val = min(255, max(0, val));

                displayImage->setPixel(c, posCt, qRgb(0, val, 0));
            }

            posCt++;
        }

        if(trainingLabel[r] == 0 && negCt < h/2)
        {
            for(c=0;c<numFeatures;c++)
                if(c < w)
            {
                int val = 255.0*(features[r*numFeatures + c]/(4.0*mean)) + 128.0;
                val = min(255, max(0, val));

                displayImage->setPixel(c, negCt + h/2, qRgb(val, 0, 0));
            }

            negCt++;
        }
    }

}

/*******************************************************************************
    AdaBoost - Computes and AdaBoost classifier using a set of candidate weak classifiers
        features - Array of feature values pre-computed for the training dataset
        trainingLabel - Ground truth labels for the training examples (1 = face, 0 = background)
        numTrainingExamples - Number of training examples
        candidateWeakClassifiers - Set of candidate weak classifiers
        numCandidateWeakClassifiers - Number of candidate weak classifiers
        weakClassifiers - Set of weak classifiers selected by AdaBoost
        numWeakClassifiers - Number of selected weak classifiers
*******************************************************************************/
void MainWindow::AdaBoost(double *features, int *trainingLabel, int numTrainingExamples,
              CWeakClassifiers *candidateWeakClassifiers, int numCandidateWeakClassifiers, CWeakClassifiers *weakClassifiers, int numWeakClassifiers)
{
    FILE *out;
    out = fopen("AdaBoost.txt", "w");

    double *scores = new double [numTrainingExamples];
    double weightSum = 0.0;
    int *featureSortIdx = new int [numTrainingExamples*numCandidateWeakClassifiers];
    double *featureTranspose = new double [numTrainingExamples*numCandidateWeakClassifiers];

    // Record the classification socres for each training example
    memset(scores, 0, numTrainingExamples*sizeof(double));

    int i, j;
    // The weighting for each training example
    double *dataWeights = new double [numTrainingExamples];

    // Begin with uniform weighting
    for(i=0;i<numTrainingExamples;i++)
        dataWeights[i] = 1.0/(double) (numTrainingExamples);


    // Let's sort the feature values for each candidate weak classifier
    for(i=0;i<numCandidateWeakClassifiers;i++)
    {
        QMap<double, int> featureSort;
        QMap<double, int>::const_iterator iterator;


        for(j=0;j<numTrainingExamples;j++)
        {
            featureSort.insertMulti(features[j*numCandidateWeakClassifiers + i], j);

            // For ease later on we'll store a transposed version of the feature array
            featureTranspose[i*numTrainingExamples + j] = features[j*numCandidateWeakClassifiers + i];
        }

        j = 0;
        iterator = featureSort.constBegin();
        // Let's remember the indices of the sorted features for later.
        while (iterator != featureSort.constEnd())
        {
            featureSortIdx[i*numTrainingExamples + j] = iterator.value();
            iterator++;
            j++;
        }
    }

    // We shouldn't need the features anymore so let's delete it.
    delete [] features;


    // Find a set of weak classifiers using AdaBoost
    for(i=0;i<numWeakClassifiers;i++)
    {
        double bestError = 99999.0;
        int bestIdx = 0;

        // For each potential weak classifier
        for(j=0;j<numCandidateWeakClassifiers;j++)
        {
            CWeakClassifiers bestClassifier;

            // Find the best threshold, polarity and weight for the candidate weak classifier
            double error = FindBestClassifier(&(featureSortIdx[j*numTrainingExamples]),
                                              &(featureTranspose[j*numTrainingExamples]),
                                              trainingLabel, dataWeights, numTrainingExamples,
                                              candidateWeakClassifiers[j], &bestClassifier);

            // Is this the best classifier found so far?
            if(error < bestError)
            {
                bestError = error;
                bestIdx = j;

                // Remember the best classifier
                bestClassifier.copy(&(weakClassifiers[i]));
            }
        }

        // Given the best weak classifier found, update the weighting of the training data.
        UpdateDataWeights(&(featureTranspose[bestIdx*numTrainingExamples]), trainingLabel, weakClassifiers[i], dataWeights, numTrainingExamples);

        // Let's compute the current error for the training dataset
        weightSum += weakClassifiers[i].m_Weight;
        double error = 0.0;
        for(j=0;j<numTrainingExamples;j++)
        {
            if(featureTranspose[bestIdx*numTrainingExamples + j] > weakClassifiers[i].m_Threshold)
            {
                scores[j] += weakClassifiers[i].m_Weight*weakClassifiers[i].m_Polarity;
            }
            else
            {
                scores[j] += weakClassifiers[i].m_Weight*(1.0 - weakClassifiers[i].m_Polarity);
            }

            if((scores[j] > 0.5*weightSum && trainingLabel[j] == 0) ||
                    (scores[j] < 0.5*weightSum && trainingLabel[j] == 1))
                error++;
        }

        // Output information that you might find useful for debugging
        fprintf(out, "Count: %d\tIdx: %d\tWeight: %lf\tError: %lf\n", i, bestIdx,
                weakClassifiers[i].m_Weight, error/(double) numTrainingExamples);
        fflush(out);
    }

    delete [] dataWeights;
    delete [] scores;
    delete [] featureSortIdx;
    delete [] featureTranspose;

    fclose(out);

   // bestClassifierClass->release();

}

/*******************************************************************************
    FindFaces - Find faces in an image
        weakClassifiers - Set of weak classifiers
        numWeakClassifiers - Number of weak classifiers
        threshold - Classifier must be above Threshold to return detected face.
        minScale, maxScale - Minimum and maximum scale to search for faces.
        faceDetections - Set of face detections
        displayImage - Display image showing detected faces.
*******************************************************************************/
void MainWindow::FindFaces(QImage inImage, CWeakClassifiers *weakClassifiers, int numWeakClassifiers, double threshold, double minScale, double maxScale,
                           QMap<double, CDetection> *faceDetections, QImage *displayImage)
{
    int w = inImage.width();
    int h = inImage.height();
    double *integralImage = new double [w*h];
    double *dImage = new double [w*h];
    double scaleMulti = 1.26;
    double scale;
    int r, c;

    ConvertColorToDouble(inImage, dImage, w, h);
    // Compute the integral image
    IntegralImage(dImage, integralImage, w, h);

    // Serach in scale space
    for(scale=minScale;scale<maxScale;scale*=scaleMulti)
    {
        // Find size of bounding box, and the step size between neighboring bounding boxes.
        int faceSize = (int) scale;
        int stepSize = max(2, faceSize/8);

        // For every possible position
        for(r=0;r<h-faceSize;r+=stepSize)
            for(c=0;c<w-faceSize;c+=stepSize)
            {
                // Compute the score of the classifier
                double score = ClassifyBox(integralImage, c, r, faceSize, weakClassifiers, numWeakClassifiers, w);

                // Is the score above threshold?
                if(score > threshold)
                {
                    CDetection detection;
                    detection.m_Score = score;
                    detection.m_Scale = scale;
                    detection.m_X = (double) c;
                    detection.m_Y = (double) r;

                    // Remember the detection
                    faceDetections->insertMulti(score, detection);
                }

            }
    }

    // Draw face bounding boxes
    DrawFace(displayImage, faceDetections);

    delete [] dImage;
    delete [] integralImage;
}

/*******************************************************************************
    DrawFace - Draw the detected faces.
        displayImage - Display image
        faceDetections - Set of face detections
*******************************************************************************/
void MainWindow::DrawFace(QImage *displayImage, QMap<double, CDetection> *faceDetections)
{
    int r, c;
    QMap<double, CDetection>::const_iterator iterator = faceDetections->constBegin();

    while(iterator != faceDetections->constEnd())
    {
        CDetection detection = iterator.value();
        int c0 = (int) detection.m_X;
        int r0 = (int) detection.m_Y;
        int size = (int) detection.m_Scale;

        for(r=r0;r<r0+size;r++)
            displayImage->setPixel(c0, r, qRgb(255, 0, 0));

        for(r=r0;r<r0+size;r++)
            displayImage->setPixel(c0 + size, r, qRgb(255, 0, 0));

        for(c=c0;c<c0+size;c++)
            displayImage->setPixel(c, r0, qRgb(255, 0, 0));

        for(c=c0;c<c0+size;c++)
            displayImage->setPixel(c, r0 + size, qRgb(255, 0, 0));

        iterator++;
    }

}


/*******************************************************************************
*******************************************************************************
*******************************************************************************

    The routines you need to implement are below

*******************************************************************************
*******************************************************************************
*******************************************************************************/

/*******************************************************************************int *featureSortIdx, double *features,
    DisplayAverageFace - Display the average face and non-face image
        displayImage - Display image, draw the average images on this image
        trainingData - Array used to store the data
        trainingLabel - Label assigned to training data (1 = face, 0 = non-face)
        numTrainingExamples - Number of training examples
        patchSize - Size of training patches in one dimension (patches have patchSize*patchSize pixels)
*******************************************************************************/

void MainWindow::DisplayAverageFace(QImage *displayImage, double *trainingData, int *trainingLabel, int numTrainingExamples, int patchSize)
{
    // Add your code here.
    FeatureDetection* featureDetection   =   new FeatureDetection();
    featureDetection->GetAveragePatch(trainingData, trainingLabel, numTrainingExamples, patchSize);
    featureDetection->DisplayAverageImage(displayImage);
}

/*******************************************************************************
    IntegralImage - Compute the integral image
        image - Input double image
        integralImage - Output integral image
        w, h - Width and height of image
*******************************************************************************/
void MainWindow::IntegralImage(double *image, double *integralImage, int width, int height)
{
    // Add your code here.
    // every pixel should be sum of pixel to it upper left.

    double** SumOfColumns = new double*[height];
    for(int row = 0; row < height; row++)
        SumOfColumns[row] = new double[width+1];

    for(int row = 0; row < height; row++)
        SumOfColumns[row][0] = 0.0;

    for(int colPixel = 1; colPixel <= width; colPixel++) {
        SumOfColumns[0][colPixel] = image[0*width + colPixel] + SumOfColumns[0][colPixel-1];
        integralImage[0*width + colPixel-1] = SumOfColumns[0][colPixel];
    }

    for(int rowPixel = 1; rowPixel < height; rowPixel++ )
        for(int colPixel = 0; colPixel < width; colPixel++) {
            SumOfColumns[rowPixel][colPixel+1] = image[rowPixel*width + colPixel] + SumOfColumns[rowPixel][colPixel];
            integralImage[rowPixel*width + colPixel] = integralImage[(rowPixel-1)* width + colPixel] + SumOfColumns[rowPixel][colPixel+1];
        }

    for(int row = 0; row < height; row++)
        delete[] SumOfColumns[row];
}

/*******************************************************************************
    SumBox - Helper function for SumBox - standard bilinear interpolation
        image - image
        x, y - Position to interpolate
        w - Width of image (integralImage)
*******************************************************************************/
double MainWindow::BilinearInterpolation(double *image, double colPixel, double rowPixel, int width)
{
    // Add your code here (or cut and paste from a previous assignment.)
    return(bilinearInterpolation->GetBilinearInterpolatedPixelValue(colPixel, rowPixel));
}

/*******************************************************************************
    SumBox - Helper function for ComputeFeatures - compute the sum of the pixels within a box.
        integralImage - integral image
        x0, y0 - Upper lefthand corner of box
        x1, y1 - Lower righthand corner of box
        w - Width of image (integralImage)
*******************************************************************************/
double MainWindow::SumBox(double *integralImage, double upperLeftCol, double upperLeftRow, double bottomRightCol, double bottomRightRow, int width)
{
    // Add your code here, use BilinearInterpolation as a helper function.

    double sum = 0.0;

    double A    =   BilinearInterpolation(integralImage, upperLeftCol, upperLeftRow, width);
    double B    =   BilinearInterpolation(integralImage, bottomRightCol, upperLeftRow, width);
    double C    =   BilinearInterpolation(integralImage, upperLeftCol, bottomRightRow, width);
    double D    =   BilinearInterpolation(integralImage, bottomRightCol, bottomRightRow, width);

    sum =   A + D - (B + C);

    return sum;
}

/*******************************************************************************
    ComputeFeatures - Compute all of the features for a specific bounding box
        integralImage - integral image
        c0, r0 - position of upper lefthand corner of bounding box
        size - Size of bounding box
        features - Array for storing computed feature values, access using features[i] for all i less than numWeakClassifiers.
        weakClassifiers - Weak classifiers
        numWeakClassifiers - Number of weak classifiers
        w - Width of image (integralImage)
*******************************************************************************/
void MainWindow::ComputeFeatures(double *integralImage, int c0, int r0, int size, double *features, CWeakClassifiers *weakClassifiers, int numWeakClassifiers, int width)
{
    int i, j;
    double upperLeftCol, upperLeftRow, bottomRightCol, bottomRightRow;
    bilinearInterpolation = BilinearInterpolationClass::getInstance(integralImage,width);

    for(i=0;i<numWeakClassifiers;i++)
    {
        features[i] = 0.0;

        for(j=0;j<weakClassifiers[i].m_NumBoxes;j++)
        {
            // Add your code to compute the sum of the pixels within each box weakClassifiers[i].m_Box[j]
            upperLeftCol     =   weakClassifiers[i].m_Box[j][0][0] * (double)size + (double)c0;
            upperLeftRow     =   weakClassifiers[i].m_Box[j][0][1] * (double)size + (double)r0;
            bottomRightCol   =   weakClassifiers[i].m_Box[j][1][0] * (double)size + (double)c0;
            bottomRightRow   =   weakClassifiers[i].m_Box[j][1][1] * (double)size + (double)r0;

            double sum = SumBox(integralImage, upperLeftCol, upperLeftRow, bottomRightCol, bottomRightRow, width);

            // Store the final feature value
            features[i] += weakClassifiers[i].m_BoxSign[j]*sum/((double) (size*size));
        }
    }

}

/*******************************************************************************
    FindBestClassifier - AdaBoost helper function.  Find the best threshold for the candidate classifier
        featureSortIdx - Indexes of the training examples sorted based on the feature responses (lowest to highest)
                Use these indices to index into the other arrays, ifeatures.e. features, trainingLabel, dataWeights.
        features - Array of feature values for the candidate classifier
        trainingLabel - Ground truth labels for the training examples (1 = face, 0 = background)
        dataWeights - Weights used to weight each training example
        numTrainingExamples - Number of training examples
        candidateWeakClassifier - Candidate classifier
        bestClassifier - Returned best classifier (updated threshold, weight and polarity)
*******************************************************************************/

double MainWindow::FindBestClassifier(int *featureSortIdx, double *features, int *trainingLabel, double *dataWeights,
                                      int numTrainingExamples, CWeakClassifiers candidateWeakClassifier, CWeakClassifiers *bestClassifier)
{
    double bestError = 99999999.0;
    candidateWeakClassifier.copy(bestClassifier); // Copy the weak classifiers params



    //BestClassifier* bestClassifierClass  =   BestClassifier::getInstance(trainingLabel, numTrainingExamples);
    bestClassifierClass  =   BestClassifier::getInstance(trainingLabel, numTrainingExamples);
    bestClassifierClass->Initialize(featureSortIdx, features, dataWeights);
    bestClassifierClass->FindOptimalClassifierParameter();


    bestError                   =   bestClassifierClass->classifierParameter.error;
    bestClassifier->m_Polarity  =   bestClassifierClass->classifierParameter.polarity;
    bestClassifier->m_Threshold =   bestClassifierClass->classifierParameter.threshold;
    bestClassifier->m_Weight    =   bestClassifierClass->classifierParameter.weight;

    // Once you find the best weak classifier, you'll need to update the following member variables:
    //      bestClassifier->m_Polarity
    //      bestClassifier->m_Threshold
    //      bestClassifier->m_Weight - this is the alpha value in the course notes

    return bestError;

}

/*******************************************************************************
    UpdateDataWeights - AdaBoost helper function.  Updates the weighting of the training examples
        features - Array of feature values for the candidate classifier
        trainingLabel - Ground truth labels for the training examples (1 = face, 0 = background)
        weakClassifier - A weak classifier
        dataWeights - Weights used to weight each training example.  These are teh weights updated.
        numTrainingExamples - Number of training examples
*******************************************************************************/
void MainWindow::UpdateDataWeights(double *features, int *trainingLabel, CWeakClassifiers weakClassifier, double *dataWeights, int numTrainingExamples)
{
    // Add you code here.
    double updatedWeightSum = 0.0;

    int e = -1;
     double beta     =   std::exp(-weakClassifier.m_Weight);
    for(int dataIndex = 0; dataIndex < numTrainingExamples; dataIndex++) {
        if(weakClassifier.m_Polarity == 0) {
            if(features[dataIndex] <= weakClassifier.m_Threshold)
                e = trainingLabel[dataIndex] == 1 ? 0 : 1;
            else
                e = trainingLabel[dataIndex] == 0 ? 0 : 1;
        } else {
            if(features[dataIndex] > weakClassifier.m_Threshold)
                e = trainingLabel[dataIndex] == 1 ? 0 : 1;
            else
                e = trainingLabel[dataIndex] == 0 ? 0 : 1;
        }

        dataWeights[dataIndex] = dataWeights[dataIndex]* std::pow(beta,(1-e));
        updatedWeightSum += dataWeights[dataIndex];
    }

    for (int dataIndex = 0; dataIndex < numTrainingExamples; ++dataIndex) {
           dataWeights[dataIndex] /= updatedWeightSum;
     }
}

/*******************************************************************************
    ClassifyBox - FindFaces helper function.  Return classification score for bounding box
        integralImage - integral image
        c0, r0 - position of upper lefthand corner of bounding box
        size - Size of bounding box
        weakClassifiers - Weak classifiers
        numWeakClassifiers - Number of weak classifiers
        w - Width of image (integralImage)
*******************************************************************************/
double MainWindow::ClassifyBox(double *integralImage, int c0, int r0, int size, CWeakClassifiers *weakClassifiers, int numWeakClassifiers, int width)
{
    // Add your code here.
    double score        =   0.0;
    double* features    =   new double[numWeakClassifiers];
    ComputeFeatures(integralImage, c0, r0, size, features, weakClassifiers, numWeakClassifiers, width);

    for(int index = 0; index < numWeakClassifiers; index++) {
        if(features[index] > weakClassifiers[index].m_Threshold)
            score +=    weakClassifiers[index].m_Weight * weakClassifiers[index].m_Polarity;
        else
            score +=    weakClassifiers[index].m_Weight * ( 1.0 - weakClassifiers[index].m_Polarity);
    }

    delete[] features;
    return score;
}

/*******************************************************************************
    NMS - Non-maximal suppression of face detections (neighboring face detections must be beyond
                xyThreshold AND scaleThreshold in position and scale respectivitely.)
        faceDetections - Set of face detections
        xyThreshold - Minimum distance in position between neighboring detections
        scaleThreshold - Minimum distance in scale between neighboring detections
        displayImage - Display image
*******************************************************************************/


void MainWindow::NMS(QMap<double, CDetection> *faceDetections, double xyThreshold, double scaleThreshold, QImage *displayImage)
{
    QMap<double, CDetection>::const_iterator iterator = faceDetections->constBegin();
    // Store the final set of face detections in finalFaceDetections
    QMap<double, CDetection> finalFaceDetections;



    // This is how you iterate through all the faces detections (lowest face detection score first.)
    while(iterator != faceDetections->constEnd())
    {

        // Add your code here.
        if(iterator == faceDetections->constBegin()) {
            finalFaceDetections.insertMulti(iterator.key(), iterator.value());
        } else {
            // loop through all the finalfaceDetections
            // if scale or XY - current facedetection < threshold
            // find highest score
            // remove from finalFacedetection and add new

             iter = finalFaceDetections.begin();
             bool shouldIncludeValue    =   true;
             for( iter ;iter != finalFaceDetections.end();) {
                 if(((abs((*iterator).m_Scale - (*iter).m_Scale)) < scaleThreshold)
                     &&(((abs((*iterator).m_X - (*iter).m_X)) < xyThreshold) &&
                         ((abs((*iterator).m_Y - (*iter).m_Y)) < xyThreshold))) {
                     if((*iter).m_Score < (*iterator).m_Score) {
                         finalFaceDetections.erase(iter++);
                         //finalFaceDetections.insertMulti(iter.key(), iter.value());
                         shouldIncludeValue = true;
                     } else {
                         iter++;
                         shouldIncludeValue = false;
                     }
                 } else
                     iter++;
             }
             if(shouldIncludeValue)
                 finalFaceDetections.insertMulti(iterator.key(), iterator.value());
        }
        // Add a face detection to finalFaceDetections using:
        // finalFaceDetections.insertMulti(iterator.key(), iterator.value());

        iterator++;
    }

    DrawFace(displayImage, &finalFaceDetections);

}
