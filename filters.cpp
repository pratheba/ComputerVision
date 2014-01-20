#include "mainwindow.h"
#include "math.h"
#include "ui_mainwindow.h"
#include "time.h"
#include <QtGui>
#include <iostream>

// Include utility class
//#include "utility.h"

/***********************************************************************
  This is the only file you need to change for your assignment.  The
  other files control the UI (in case you want to make changes.)
************************************************************************/


// The first four functions provide example code to help get you started

// Convert an image to grey-scale
void MainWindow::BlackWhiteImage(QImage *image)
{
    int r, c;
    QRgb pixel;

    for(r=0;r<image->height();r++)
    {
        for(c=0;c<image->width();c++)
        {
            pixel = image->pixel(c, r);
            double red = (double) qRed(pixel);
            double green = (double) qGreen(pixel);
            double blue = (double) qBlue(pixel);

            // Compute intensity from colors - these are common weights
            double intensity = 0.3*red + 0.6*green + 0.1*blue;

            image->setPixel(c, r, qRgb( (int) intensity, (int) intensity, (int) intensity));
        }
    }
}

// Add random noise to the image
void MainWindow::AddNoise(QImage *image, double mag, bool colorNoise)
{
    int r, c;
    QRgb pixel;
    int noiseMag = mag;
    noiseMag *= 2;

    for(r=0;r<image->height();r++)
    {
        for(c=0;c<image->width();c++)
        {
            pixel = image->pixel(c, r);
            int red = qRed(pixel);
            int green = qGreen(pixel);
            int blue = qBlue(pixel);

            // If colorNoise, add color independently to each channel
            if(colorNoise)
            {
                red += rand()%noiseMag - noiseMag/2;
                green += rand()%noiseMag - noiseMag/2;
                blue += rand()%noiseMag - noiseMag/2;
            }
            // otherwise add the same amount of noise to each channel
            else
            {
                int noise = rand()%noiseMag - noiseMag/2;

                red += noise;
                green += noise;
                blue += noise;
            }

            // Make sure we don't over or under saturate
            red = min(255, max(0, red));
            green = min(255, max(0, green));
            blue = min(255, max(0, blue));

            image->setPixel(c, r, qRgb( red, green, blue));
        }
    }
}

// Here is an example of blurring an image using a mean or box filter with the specified radius.
// This could be implemented using separable filters to make it much more efficient, but it is not.
void MainWindow::MeanBlurImage(QImage *image, int radius)
{
    if(radius == 0)
        return;

    int r, c, rd, cd, i;
    QRgb pixel;

    // This is the size of the kernel
    int size = 2*radius + 1;

    // Create a buffer image so we're not reading and writing to the same image during filtering.
    QImage buffer;
    int w = image->width();
    int h = image->height();

    // This creates an image of size (w + 2*radius, h + 2*radius) with black borders.
    // This could be improved by filling the pixels using a different padding technique (reflected, fixed, etc.)
    buffer = image->copy(-radius, -radius, w + 2*radius, h + 2*radius);

    // Compute kernel to convolve with the image.
    double *kernel = new double [size*size];

    for(i=0;i<size*size;i++)
    {
        kernel[i] = 1.0;
    }

    // Make sure kernel sums to 1
    double denom = 0.000001;
    for(i=0;i<size*size;i++)
        denom += kernel[i];
    for(i=0;i<size*size;i++)
        kernel[i] /= denom;

    // For each pixel in the image...
    for(r=0;r<h;r++)
    {
        for(c=0;c<w;c++)
        {
            double rgb[3];

            rgb[0] = 0.0;
            rgb[1] = 0.0;
            rgb[2] = 0.0;

            // Convolve the kernel at each pixel
            for(rd=-radius;rd<=radius;rd++)
                for(cd=-radius;cd<=radius;cd++)
                {
                     // Get the pixel value
                     pixel = buffer.pixel(c + cd + radius, r + rd + radius);

                     // Get the value of the kernel
                     double weight = kernel[(rd + radius)*size + cd + radius];

                     rgb[0] += weight*(double) qRed(pixel);
                     rgb[1] += weight*(double) qGreen(pixel);
                     rgb[2] += weight*(double) qBlue(pixel);
                }

            // Store mean pixel in the image to be returned.
            image->setPixel(c, r, qRgb((int) floor(rgb[0] + 0.5), (int) floor(rgb[1] + 0.5), (int) floor(rgb[2] + 0.5)));
        }
    }

    // Clean up.
    delete [] kernel;
}

// Downsample the image by 1/2
void MainWindow::HalfImage(QImage &image)
{
    QImage buffer;
    int w = image.width();
    int h = image.height();
    int r, c;

    buffer = image.copy();

    // Reduce the image size.
    image = QImage(w/2, h/2, QImage::Format_RGB32);

    // Copy every other pixel
    for(r=0;r<h/2;r++)
        for(c=0;c<w/2;c++)
        {
             image.setPixel(c, r, buffer.pixel(c*2, r*2));
        }
}

#pragma region Utility
struct ImageMetaData {
    int Width;
    int Height;
    int radius;
    ImageMetaData():Width(0),Height(0),radius(0){}
    ImageMetaData(int width_, int height_, int radius_):Width(width_),Height(height_),radius(radius_){}
    ImageMetaData(const ImageMetaData& metaData):Width(metaData.Width),Height(metaData.Height),radius(metaData.radius){}
};

void ResetValues(double* imageRGB, double** tempStorageHorizConvolution, int size ) {
     
    if(imageRGB != NULL) {
     imageRGB[0] = 0.0;
     imageRGB[1] = 0.0;
     imageRGB[2] = 0.0;
    }

    if((tempStorageHorizConvolution != NULL) && size > 0) {
     for(int i=0; i< size; i++) 
         for(int j=0; j<3; j++)
             tempStorageHorizConvolution[i][j] = 0.0;
    }
}


QImage InitializeOutputImageBuffer(const QImage* image, double sigma, ImageMetaData& metaData)
{
    metaData.radius            =   (int)(3*sigma);
    metaData.Width             =   image->width();
    metaData.Height            =   image->height();
    double* kernel              =  NULL;

    QImage GaussImagebuffer     =   image->copy(-metaData.radius,-metaData.radius, 
        metaData.Width + 2*metaData.radius, metaData.Height + 2*metaData.radius);
    return (GaussImagebuffer);
}
#pragma endregion Utility

/////////////////////////////////////////////// GAUSSIAN BLUR ///////////////////////////////////////////////////
#pragma region GaussianBlurImage

//struct GaussBlurImage {
//	QImage* GaussImageBuffer;
//    double* GaussianKernel;
//
//    GaussBlurImage():GaussianKernel(NULL),GaussImageBuffer(NULL){}
//	GaussBlurImage(const GaussImage& GaussImage_): GaussImageBuffer(GaussImage_.GaussImageBuffer){}
//    GaussBlurImage(QImage* GaussImageBuffer_):GaussImageBuffer(GaussImageBuffer_){}
//};



struct Position {
    int XPixelPos;
    int YPixelPos;

    Position():XPixelPos(-1), YPixelPos(-1){}
    Position(int X, int Y):XPixelPos(X), YPixelPos(Y){}
    Position(const Position& position_):XPixelPos(position_.XPixelPos), YPixelPos(position_.YPixelPos){}
};



#pragma region GaussianKernel

struct GaussianKernel {
    int KernelSize;
    double variance;
    double* Kernel;

    GaussianKernel():KernelSize(0), variance(0), Kernel(NULL){}
    GaussianKernel(int kernelSize_, double variance_, double* Kernel): 
        KernelSize(kernelSize_), variance(variance_), Kernel(Kernel){}
    GaussianKernel(const GaussianKernel& gaussianKernel_): 
        KernelSize(gaussianKernel_.KernelSize), variance(gaussianKernel_.variance), Kernel(gaussianKernel_.Kernel){}
};

void InitializeGaussianKernel(const int radius, GaussianKernel& gaussianKernel) {
    gaussianKernel.KernelSize   =       2 * radius + 1;
    gaussianKernel.Kernel       =       new double[gaussianKernel.KernelSize * gaussianKernel.KernelSize];
}

void LoadvaluesOfkernelComponents(const int radius,GaussianKernel& gaussianKernel)
{
    double kernelConstant = (double)1 / (double)(2 * M_PI * powf((float)gaussianKernel.variance,(float)2));

    for(int row = -radius; row <= radius; ++row) {
        for(int col = -radius; col <= radius; ++col) {
            int Xvalue      =   (int)powf(float(row ), (float)2);
            int Yvalue      =   (int)powf(float(col ), (float)2);
            double expValue    =   (double)(-1)*((Xvalue + Yvalue)/(2*(powf((float)gaussianKernel.variance,(float)2))));

            gaussianKernel.Kernel[(row + radius)*gaussianKernel.KernelSize + col + radius] = kernelConstant * exp(expValue);
        }   }
}

void NormalizeKernel(const GaussianKernel& kernel) {
    double epsilon = 0.000001;
    double magnitude = epsilon;
    for(int pos = 0; pos < kernel.KernelSize * kernel.KernelSize ; pos++) {
         magnitude +=   kernel.Kernel[pos];
    }

    for(int pos = 0; pos < kernel.KernelSize * kernel.KernelSize ; pos++)
        kernel.Kernel[pos] /= magnitude;
}
    
void GetGaussianKernel(const int radius,GaussianKernel& gaussianKernel) {
    InitializeGaussianKernel(radius, gaussianKernel);
    LoadvaluesOfkernelComponents(radius, gaussianKernel);
    NormalizeKernel(gaussianKernel);
}
#pragma endregion GaussianKernel

double* ApplyFilterToPixelAndReturnRGB(const ImageMetaData& metaData, Position position, const QImage& outputImageBuffer, const GaussianKernel& kernel)
{
    int radius                = metaData.radius;
    double* imageRGB          = new double[3]();
    QRgb pixel;

    for(int row = -radius; row <= radius; ++row) {
        for(int col = -radius; col <= radius; ++col) {
            pixel           =   outputImageBuffer.pixel(position.YPixelPos + col + radius, position.XPixelPos + row + radius);
            double  weight  =   kernel.Kernel[(row + radius)* kernel.KernelSize + (col + radius)];

            imageRGB[0]     +=  weight*(double) qRed(pixel);
            imageRGB[1]     +=  weight*(double) qGreen(pixel);
            imageRGB[2]     +=  weight*(double) qBlue(pixel);
        }
    }

        return imageRGB;
}

void ConvolveImagewithGaussianKernel(const ImageMetaData& metaData, QImage* inputImage, QImage& outputImageBuffer,GaussianKernel& kernel) {
   
    double* imageRGB = new double[3];
    int radius = metaData.radius;
    QRgb pixel;

    for(int rowPixel = 0; rowPixel < metaData.Height; rowPixel++) {
        for(int colPixel = 0; colPixel < metaData.Width; colPixel++) {
            //imageRGB = ApplyFilterToPixelAndReturnRGB(metaData, Position(rowPixel, colPixel), outputImageBuffer, kernel);

            ResetValues(imageRGB,NULL,0);
             for(int row = -radius; row <= radius; ++row) {
        for(int col = -radius; col <= radius; ++col) {

            pixel           =   outputImageBuffer.pixel(colPixel + col + radius, rowPixel + row + radius);
            double  weight  =   kernel.Kernel[(row + radius)* kernel.KernelSize + (col + radius)];

            imageRGB[0]     +=  weight*(double) qRed(pixel);
            imageRGB[1]     +=  weight*(double) qGreen(pixel);
            imageRGB[2]     +=  weight*(double) qBlue(pixel);
        }
    }
            // Store mean pixel in the image to be returned.
            inputImage->setPixel(colPixel, rowPixel, qRgb((int) floor(imageRGB[0] + 0.5), 
                (int) floor(imageRGB[1] + 0.5), (int) floor(imageRGB[2] + 0.5)));
            }
}

        delete[] imageRGB;
        imageRGB = NULL;
}

void MainWindow::GaussianBlurImage(QImage *image, double sigma)
{
    // Add your code here.  Look at MeanBlurImage to get yourself started.

	// Assume radius to be one.
     clock_t start, end;
    ImageMetaData metaData(0,0,0);
    GaussianKernel gaussianKernel(0,sigma,NULL);

	QImage OutPutImageBuffer        = InitializeOutputImageBuffer(image,sigma, metaData);
    GetGaussianKernel(metaData.radius, gaussianKernel);

    start = clock();
    ConvolveImagewithGaussianKernel(metaData, image,OutPutImageBuffer,gaussianKernel);
    end = clock();
    qDebug() << "Time required for execution: "   << (double)(end-start)/CLOCKS_PER_SEC ;
}

#pragma endregion GaussianBlurImage

/////////////////////////////////////////////// SEPERABLE GAUSSIAN BLUR ///////////////////////////////////////////
#pragma region SeperableGaussianBlurImage

void InitializeSeperableGaussianKernel(const int radius, GaussianKernel& gaussianKernel) {
    gaussianKernel.KernelSize   =       2 * radius + 1;
    gaussianKernel.Kernel       =       new double[gaussianKernel.KernelSize];
}

void LoadvaluesOfSeperablekernelComponents(const int radius,GaussianKernel& gaussianKernel)
{
    double kernelConstant = (double)1 / (double)(sqrtf(2 * M_PI) * (float)gaussianKernel.variance);

    for(int pos = -radius; pos <= radius; ++pos) {
            int value                               =   (int)powf(float(pos), (float)2);
            double expValue                         =   (double)(-1)*((value)/(2*(powf((float)gaussianKernel.variance,(float)2))));
            gaussianKernel.Kernel[(pos + radius)]   =   kernelConstant * exp(expValue);
    }
}

void NormalizeSeperableKernel(const GaussianKernel& kernel) {
    double epsilon = 0.000001;
    double magnitude = epsilon;

    for(int pos = 0; pos < kernel.KernelSize ; pos++) 
         magnitude +=   kernel.Kernel[pos];

    for(int pos = 0; pos < kernel.KernelSize ; pos++)
        kernel.Kernel[pos] /= magnitude;
}

void GetSeperableGaussianKernel(const int radius,GaussianKernel& gaussianKernel) {
    InitializeSeperableGaussianKernel(radius, gaussianKernel);
    LoadvaluesOfSeperablekernelComponents(radius, gaussianKernel);
    NormalizeSeperableKernel(gaussianKernel);
}

//typedef double* (*callback_function)(const QImage&,const GaussianKernel&, Position,int ); // type for conciseness

//void ApplyVerticalFilterToPixelAndReturnRGB(double* imageRGB, double** tempStorageHorizConvolution, const GaussianKernel& kernel, int radius ) {
//    
//    for(int pos = -radius; pos <= radius; pos++) {
//        imageRGB[0] += (tempStorageHorizConvolution[pos+radius][0] )*kernel.Kernel[pos + radius];
//        imageRGB[1] += (tempStorageHorizConvolution[pos+radius][1] )*kernel.Kernel[pos + radius];
//        imageRGB[2] += (tempStorageHorizConvolution[pos+radius][2] )*kernel.Kernel[pos + radius];
//    }
//}

//void ApplyHorizontalFilterToPixelAndReturnRGB(const QImage& outputImageBuffer,double** tempStorageHorizConvolution,const GaussianKernel& kernel, Position position, int radius) {
//    QRgb pixel;
//
//    for(int row = -radius; row <= radius; ++row) 
//        for(int col = -radius; col <= radius; ++col) {
//            pixel           =   outputImageBuffer.pixel(position.YPixelPos + col + radius,position.XPixelPos +row + radius); 
//            double  weight  =   kernel.Kernel[col + radius];
//            
//            tempStorageHorizConvolution[row+radius][0]     +=  weight*(double) qRed(pixel);
//            tempStorageHorizConvolution[row+radius][1]     +=  weight*(double) qGreen(pixel);
//            tempStorageHorizConvolution[row+radius][2]     +=  weight*(double) qBlue(pixel);
//        }    
//}

//double* ApplyFilterToPixelSeperableAndReturnRGB(const ImageMetaData& metaData, callback_function sample, Position position, const QImage& outputImageBuffer, const GaussianKernel& kernel)
//{
//    int radius                = metaData.radius;
//    double* imageRGB          = new double[3]();
//    QRgb pixel;
//
//    for(int row = -radius; row <= radius; ++row) {
//       for(int col = -radius; col <= radius; ++col) {
//            //pixel           =   sample()
//            //pixel           =   outputImageBuffer.pixel(position.YPixelPos + col + radius, position.XPixelPos + row + radius);
//            double  weight  =   kernel.Kernel[col + radius];
//
//            imageRGB[0]     +=  weight*(double) qRed(pixel);
//            imageRGB[1]     +=  weight*(double) qGreen(pixel);
//            imageRGB[2]     +=  weight*(double) qBlue(pixel);
//        }
//
//    }
//      
//    //for(int pos = -radius; pos <= radius; ++pos) {
//    //        pixel           =   sample(outputImageBuffer, position, radius, pos );//outputImageBuffer.pixel(position.YPixelPos , position.XPixelPos + pos + radius);
//    //        double  weight  =   kernel.Kernel[pos + radius];
//
//    //        imageRGB[0]     +=  weight*(double) qRed(pixel);
//    //        imageRGB[1]     +=  weight*(double) qGreen(pixel);
//    //        imageRGB[2]     +=  weight*(double) qBlue(pixel);
//    //    }
//        return imageRGB;
//}



// Reflected pixel
// edgeSize is either width or height of image (minus) 1 
int GetPixelPositionForEdges(int currentPixel, int edgesize) {
    if(currentPixel < 0)
        return abs(currentPixel);
    else if(currentPixel > (edgesize))
        return abs(edgesize-(currentPixel - edgesize));
}

void  ApplyVerticalFilterToPixel(QImage* inputImage, QImage* tempImage,const ImageMetaData& metaData,const GaussianKernel& kernel ) {
    QRgb pixel;
    int edgePixel;
    double* imageRGB = new double[3];
    int radius = metaData.radius;

    for(int rowPixel = 0; rowPixel < metaData.Height; rowPixel++) { // X
        for(int colPixel = 0; colPixel < metaData.Width; colPixel++) { // Y
            
            ResetValues(imageRGB, NULL, 0);
             for(int row = -radius; row <= radius; ++row) {

                 if((rowPixel+row) < 0 || (rowPixel+row)> (metaData.Height-1)) {
                     edgePixel = GetPixelPositionForEdges(rowPixel+row,metaData.Height-1);
                     pixel = tempImage->pixel(colPixel, edgePixel);
                 }
                 else
                     pixel = tempImage->pixel(colPixel,rowPixel+row);

                 double weight = kernel.Kernel[row+radius];

                 imageRGB[0] += (double)qRed(pixel) * weight;
                 imageRGB[1] += (double)qGreen(pixel) * weight;
                 imageRGB[2] += (double)qBlue(pixel) * weight;
             }
                  inputImage->setPixel(colPixel, rowPixel, qRgb((int) floor(imageRGB[0] + 0.5), 
                (int) floor(imageRGB[1] + 0.5), (int) floor(imageRGB[2] + 0.5))); 
        }
    }
}

void ApplyHorizontalFilterToPixel(QImage* inputImage, QImage* tempImage,const ImageMetaData& metaData,const GaussianKernel& kernel ) {
    QRgb pixel;
    int edgePixel;
    double* imageRGB = new double[3];
    int radius = metaData.radius;

    for(int rowPixel = 0; rowPixel < metaData.Height; rowPixel++) { // X
        for(int colPixel = 0; colPixel < metaData.Width; colPixel++) { // Y
            
            ResetValues(imageRGB, NULL, 0);
             for(int col = -radius; col <= radius; ++col) {

                 if((colPixel+col) < 0 || (colPixel+col)> (metaData.Width-1)) {
                     edgePixel = GetPixelPositionForEdges(colPixel+col,metaData.Width-1);
                     pixel = inputImage->pixel(edgePixel, rowPixel);
                 }
                 else
                     pixel = inputImage->pixel(colPixel+col,rowPixel);

                 double weight = kernel.Kernel[col+radius];

                 imageRGB[0] += (double)qRed(pixel) * weight;
                 imageRGB[1] += (double)qGreen(pixel) * weight;
                 imageRGB[2] += (double)qBlue(pixel) * weight;
             }
                  tempImage->setPixel(colPixel, rowPixel, qRgb((int) floor(imageRGB[0] ), 
                (int) floor(imageRGB[1] ), (int) floor(imageRGB[2] ))); 
        }
    }
}

void ConvolveImagewithSeperableGaussianKernel(const ImageMetaData& metaData, QImage* inputImage, QImage& outputImageBuffer,GaussianKernel& kernel) {
     
   QImage tempImage = inputImage->copy();
   
    ApplyHorizontalFilterToPixel(inputImage, &tempImage, metaData,kernel);
    ApplyVerticalFilterToPixel(inputImage,&tempImage,metaData, kernel);
}

void MainWindow::SeparableGaussianBlurImage(QImage *image, double sigma)
{
    // Add your code here.  Done right, you should be able to copy most of the code from GaussianBlurImage.
    clock_t start, end;
    ImageMetaData metaData(0,0,0);
    GaussianKernel gaussianKernel(0,sigma,NULL);
   
    QImage OutPutImageBuffer        =       InitializeOutputImageBuffer(image,sigma,metaData);
    GetSeperableGaussianKernel(metaData.radius, gaussianKernel);

    start = clock();
    ConvolveImagewithSeperableGaussianKernel(metaData, image,OutPutImageBuffer,gaussianKernel);
    end = clock();
    qDebug() << "Time required for execution: "   << (double)(end-start)/CLOCKS_PER_SEC ;
}

#pragma endregion SeperableGaussianBlurImage

///////////////////////////////////////////// DERIVATIVE IMAGE ////////////////////////////////////////////////////
#pragma region DerivativeImage
#pragma region FirstDerivativeImage

double** GetSobelFilter(int radius) {
    double** SobelFilter = new double* [2*radius +1];
    for(int i = 0; i < (2*radius+1); i++) {
        SobelFilter[i] = new double[3];
        for(int j=0; j<3;j++)
            SobelFilter[i][j] = 0.0;
    }

    // change from default to adapt to filter of any size.

    SobelFilter[0][0] = -1; SobelFilter[1][0] = -2; SobelFilter[2][0] = -1; 
    SobelFilter[0][1] = -1; SobelFilter[1][1] = -2; SobelFilter[2][1] = -1;
    SobelFilter[0][2] = -1; SobelFilter[1][2] = -2; SobelFilter[2][2] = -1;

    return SobelFilter;
}

void GetIntensityDifferenceBetweenneighbourPixelsinXDirection(const QImage& OutPutImageBuffer,double** temporaryRowArray, int radius, int rowPixel) {
    QRgb pixel;

    for(int colPixel = radius; colPixel < OutPutImageBuffer.width() - radius; colPixel++) {

            // If we want to include the differentiation filter size > 3
            // double intensityDifference = 0.0;
            // int XDerivativeFilter[2*radius];// = {-1, 0, 1};
            // for( int pos = -radius; pos <= radius; pos++) 
            //      intensityDifference += OutPutImageBuffer.pixel(colPixel+pos) * XDerivativeFilter (pos+radius); 
           
            pixel = OutPutImageBuffer.pixel(colPixel+1,rowPixel) - OutPutImageBuffer.pixel(colPixel -1, rowPixel);
            temporaryRowArray[colPixel-radius][0] = qRed(pixel);
            temporaryRowArray[colPixel-radius][1] = qGreen(pixel);
            temporaryRowArray[colPixel-radius][2] = qBlue(pixel);
        }

}

void SetFirstDerivativePixelIinXToImage(QImage& image,double** temporaryRowArray, int rowPixel, int width) {
      for(int colPixel = 0; colPixel < width ; colPixel++) {
            image.setPixel(colPixel,rowPixel,qRgb((int) floor(temporaryRowArray[colPixel][0] + 128), 
                (int) floor(temporaryRowArray[colPixel][1] + 128), (int) floor(temporaryRowArray[colPixel][2] + 128)));

            temporaryRowArray[colPixel][0] = 0.0;
            temporaryRowArray[colPixel][1] = 0.0;
            temporaryRowArray[colPixel][2] = 0.0;
      }
}

void ApplyFirstDerivativeinXDirectionToImage(QImage& OutPutImageBuffer,QImage& image,const int radius) {

    QRgb pixel;

    int width                       =   image.width();
    int height                      =   image.height();
    double** temporaryRowArray      =   new double*[width]();

    for(int i = 0; i < width ; i++) 
        temporaryRowArray[i] = new double[3];
    
    ResetValues(NULL, temporaryRowArray, width);

    for(int rowPixel = radius; rowPixel < height+radius; rowPixel++) {
        GetIntensityDifferenceBetweenneighbourPixelsinXDirection(OutPutImageBuffer,temporaryRowArray, radius, rowPixel);
        SetFirstDerivativePixelIinXToImage(image, temporaryRowArray, rowPixel - radius, width);
    }

     for(int i = 0; i < width ; i++) 
        delete[] temporaryRowArray[i];
}

void MainWindow::FirstDerivImage(QImage *image, double sigma)
{
    // Add your code here.
    ImageMetaData metaData(0,0,0);
    GaussianKernel gaussianKernel(0,sigma,NULL);
   
    QImage OutPutImageBuffer        =       InitializeOutputImageBuffer(image,sigma,metaData);
    ApplyFirstDerivativeinXDirectionToImage(OutPutImageBuffer, *image,metaData.radius);
    SeparableGaussianBlurImage(image,sigma);
}

#pragma endregion FirstDerivativeImage

#pragma region SecondDerivativeImage

void SetSecondDerivativePixelIinYToImage(QImage& image,double** temporaryRowArray, int colPixel, int height) {
      for(int rowPixel = 0; rowPixel < height ; rowPixel++) {
            image.setPixel(colPixel,rowPixel,qRgb((int) floor(temporaryRowArray[rowPixel][0] + 128), 
                (int) floor(temporaryRowArray[rowPixel][1] + 128), (int) floor(temporaryRowArray[rowPixel][2] + 128)));

            temporaryRowArray[rowPixel][0] = 0.0;
            temporaryRowArray[rowPixel][1] = 0.0;
            temporaryRowArray[rowPixel][2] = 0.0;
      }
}

void GetSecondOrderDiffOfPixelsinYDirection(const QImage& OutPutImageBuffer,double** temporaryRowArray, int radius, int colPixel) {
    QRgb pixel;

    for(int rowPixel = radius; rowPixel < OutPutImageBuffer.height() - radius; rowPixel++) {

            // If we want to include the differentiation filter size > 3
            // double intensityDifference = 0.0;
            // int XDerivativeFilter[2*radius];// = {-1, 0, 1};
            // for( int pos = -radius; pos <= radius; pos++) 
            //      intensityDifference += OutPutImageBuffer.pixel(colPixel+pos) * XDerivativeFilter (pos+radius); 
           
        pixel = OutPutImageBuffer.pixel(colPixel,rowPixel + 1) + OutPutImageBuffer.pixel(colPixel, rowPixel-1) - 2*OutPutImageBuffer.pixel(colPixel,rowPixel);
            temporaryRowArray[rowPixel-radius][0] = qRed(pixel);
            temporaryRowArray[rowPixel-radius][1] = qGreen(pixel);
            temporaryRowArray[rowPixel-radius][2] = qBlue(pixel);
        }

}

void ApplySecondDerivativeinYDirectionToImage(QImage& OutPutImageBuffer,QImage& image,const int radius) {

    QRgb pixel;

    int width                       =   image.width();
    int height                      =   image.height();
    double** temporaryRowArray      =   new double*[height]();

    for(int i = 0; i < height ; i++) 
        temporaryRowArray[i] = new double[3];
    
    ResetValues(NULL, temporaryRowArray, height);

    for(int colPixel = radius; colPixel < width+radius; colPixel++) {
        GetSecondOrderDiffOfPixelsinYDirection(OutPutImageBuffer,temporaryRowArray, radius, colPixel);
        SetSecondDerivativePixelIinYToImage(image, temporaryRowArray, colPixel - radius, height);
    }

     for(int i = 0; i < height ; i++) 
        delete[] temporaryRowArray[i];
}

void GetSecondOrderDiffOfPixelsinXDirection(const QImage& OutPutImageBuffer,double** temporaryRowArray, int radius, int rowPixel) {
    QRgb pixel;

    for(int colPixel = radius; colPixel < OutPutImageBuffer.width() - radius; colPixel++) {

            // If we want to include the differentiation filter size > 3
            // double intensityDifference = 0.0;
            // int XDerivativeFilter[2*radius];// = {-1, 0, 1};
            // for( int pos = -radius; pos <= radius; pos++) 
            //      intensityDifference += OutPutImageBuffer.pixel(colPixel+pos) * XDerivativeFilter (pos+radius); 
           
            pixel = OutPutImageBuffer.pixel(colPixel+1,rowPixel) + OutPutImageBuffer.pixel(colPixel -1, rowPixel) - 2*OutPutImageBuffer.pixel(colPixel,rowPixel);
            temporaryRowArray[colPixel-radius][0] = qRed(pixel);
            temporaryRowArray[colPixel-radius][1] = qGreen(pixel);
            temporaryRowArray[colPixel-radius][2] = qBlue(pixel);
        }

}

void SetSecondDerivativePixelIinXToImage(QImage& image,double** temporaryRowArray, int rowPixel, int width) {
      for(int colPixel = 0; colPixel < width ; colPixel++) {
            image.setPixel(colPixel,rowPixel,qRgb((int) floor(temporaryRowArray[colPixel][0] + 128), 
                (int) floor(temporaryRowArray[colPixel][1] + 128), (int) floor(temporaryRowArray[colPixel][2] + 128)));

            temporaryRowArray[colPixel][0] = 0.0;
            temporaryRowArray[colPixel][1] = 0.0;
            temporaryRowArray[colPixel][2] = 0.0;
      }
}

void ApplySecondDerivativeinXDirectionToImage(QImage& OutPutImageBuffer,QImage& image,const int radius) {

    QRgb pixel;

    int width                       =   image.width();
    int height                      =   image.height();
    double** temporaryRowArray      =   new double*[width]();

    for(int i = 0; i < width ; i++) 
        temporaryRowArray[i] = new double[3];
    
    ResetValues(NULL, temporaryRowArray, width);

    for(int rowPixel = radius; rowPixel < height+radius; rowPixel++) {
        GetSecondOrderDiffOfPixelsinXDirection(OutPutImageBuffer,temporaryRowArray, radius, rowPixel);
        SetSecondDerivativePixelIinXToImage(image, temporaryRowArray, rowPixel - radius, width);
    }

     for(int i = 0; i < width ; i++) 
        delete[] temporaryRowArray[i];
}

QImage ReplicateImage(const QImage* originalImage) {
    QImage replicatedImage     =   originalImage->copy(0, 0, originalImage->width() , originalImage->height());
    return replicatedImage;
}

void LaplacianOfGaussian(QImage& image,  const QImage& FirstDerivativeImage, const QImage& SecondDerivativeImage, const ImageMetaData& metaData) {
    
    int width = metaData.Width;
    int height = metaData.Height;

    QRgb pixelFirstDerivImage;
    QRgb pixelSecondDerivImage;
    double* temporaryColorBuffer = new double[3];

    ResetValues(temporaryColorBuffer, NULL, 0);

    for(int rowPixel = 0; rowPixel < height; rowPixel++) {
        for(int colPixel =0; colPixel < width; colPixel++) {
            pixelFirstDerivImage = FirstDerivativeImage.pixel(colPixel,rowPixel);
            pixelSecondDerivImage = SecondDerivativeImage.pixel(colPixel, rowPixel);

            temporaryColorBuffer[0] = qRed(pixelFirstDerivImage) + qRed(pixelSecondDerivImage);
            temporaryColorBuffer[1] = qGreen(pixelFirstDerivImage) + qGreen(pixelSecondDerivImage);
            temporaryColorBuffer[2] = qBlue(pixelFirstDerivImage) + qBlue(pixelSecondDerivImage);

             image.setPixel(colPixel,rowPixel,qRgb((int) floor(temporaryColorBuffer[0] + 128), 
                (int) floor(temporaryColorBuffer[1] + 128), (int) floor(temporaryColorBuffer[2] + 128)));

             ResetValues(temporaryColorBuffer, NULL, 0);
        }
    }

    delete[] temporaryColorBuffer;
}

// This is implemented with 3x3 Laplacian operator.
// Change to dynamic filter
void MainWindow::SecondDerivImage(QImage *image, double sigma)
{
    QImage FirstDerivativeImage, SecondDerivativeImage ;
    ImageMetaData metaData(0,0,0);
    GaussianKernel gaussianKernel(0,sigma,NULL);
   
    FirstDerivativeImage = ReplicateImage(image);
    SecondDerivativeImage = ReplicateImage(image);

    QImage OutputImageBuffer        =       InitializeOutputImageBuffer(image,sigma,metaData);
    ApplySecondDerivativeinXDirectionToImage(OutputImageBuffer, FirstDerivativeImage,metaData.radius); 
    SeparableGaussianBlurImage(&FirstDerivativeImage,sigma);
    
    OutputImageBuffer               =       InitializeOutputImageBuffer(image,sigma,metaData);
    ApplySecondDerivativeinYDirectionToImage(OutputImageBuffer, SecondDerivativeImage,metaData.radius); 
    SeparableGaussianBlurImage(&SecondDerivativeImage,sigma);
 
    LaplacianOfGaussian(*image, FirstDerivativeImage, SecondDerivativeImage, metaData);

}

#pragma endregion SecondDerivativeImage

#pragma endregion DerivativeImage

/////////////////////////////////////////// SHARPEN IMAGE /////////////////////////////////////////////////////////
#pragma region SharpernImage
void ApplySubtractionFilterToImage(QImage *image, QImage* derivativeImage, double alpha) {
    
    int width       =   image->width();
    int height      =   image->height();
    QRgb imagepixel;
    QRgb derivativePixel;
    double* imageRGB          = new double[3]();

    for(int rowPixel = 0; rowPixel < height; rowPixel++)
        for(int colPixel = 0; colPixel < width; colPixel++) {
            imagepixel =  image->pixel(colPixel,rowPixel);
            derivativePixel = derivativeImage->pixel(colPixel,rowPixel);
            imageRGB[0] = (double)(qRed(imagepixel)      - (alpha * (double)(qRed(derivativePixel))));
            imageRGB[1] = (double)(qGreen(imagepixel)    - (alpha * (double)(qGreen(derivativePixel))));
            imageRGB[2] = (double)(qBlue(imagepixel)     - (alpha * (double)(qBlue(derivativePixel))));

            /*imageRGB[0] = min(255.0, max(0.0, imageRGB[0]));
            imageRGB[1] = min(255.0, max(0.0, imageRGB[1]));
            imageRGB[2] = min(255.0, max(0.0, imageRGB[2]));*/
               /* imageRGB[0] = min(255.0, max(0.0, imageRGB[0]+128));
                imageRGB[1] = min(255.0, max(0.0, imageRGB[1]+128));
                imageRGB[2] = min(255.0, max(0.0, imageRGB[2]+128));
               */

            image->setPixel(colPixel, rowPixel, qRgb((int) floor(imageRGB[0] +128), 
                (int) floor(imageRGB[1] +128), (int) floor(imageRGB[2] +128)));

            ResetValues(imageRGB, NULL, 0);
        }

}

void MainWindow::SharpenImage(QImage *image, double sigma, double alpha)
{
    // Add your code here.  It's probably easiest to call SecondDerivImage as a helper function.
    int radius = (int)(3 * sigma);
    QImage DerivativeImage     =   image->copy(0, 0, image->width() , image->height());
    SecondDerivImage(&DerivativeImage,sigma);
    ApplySubtractionFilterToImage(image, &DerivativeImage, alpha);
}

#pragma endregion SharpernImage


void MainWindow::BilateralImage(QImage *image, double sigmaS, double sigmaI)
{
    // Add your code here.  Should be similar to GaussianBlurImage.
}


struct SobelKernel {
    double* XSobelKernel;
    double* YSobelKernel;
    int KernelSize;

    SobelKernel():KernelSize(0), XSobelKernel(NULL), YSobelKernel(NULL){}
    SobelKernel(int kernelSize_,  double* XKernel_, double* YKernel_): 
        KernelSize(kernelSize_), XSobelKernel(XKernel_), YSobelKernel(YKernel_) {}
    SobelKernel(const SobelKernel& SobelKernel_): 
        KernelSize(SobelKernel_.KernelSize), XSobelKernel(SobelKernel_.XSobelKernel),YSobelKernel(SobelKernel_.YSobelKernel){}
};

void InitializeSobelKernel(SobelKernel& sobelKernel) {
    sobelKernel.KernelSize   =       3;
    sobelKernel.XSobelKernel       =       new double[3];
    sobelKernel.YSobelKernel = new double[3];
}

void LoadvaluesOfkernelComponents(SobelKernel& sobelKernel)
{
    double XKernel[3] = {1,0,-1};
    double YKernel[3] = {1,2,1};

    sobelKernel.XSobelKernel    =   XKernel;
    sobelKernel.YSobelKernel    =   YKernel;
}

struct SobelImagePixel {
    double magnitude;
    double orientation;
    double Gx;
    double Gy;

    SobelImagePixel(){}
    SobelImagePixel(double magnitude_, double orientation_):magnitude(magnitude_), orientation(orientation_){}
};


void GetSobelKernel(SobelKernel& sobelKernel) {

    InitializeSobelKernel(sobelKernel);
    LoadvaluesOfkernelComponents(sobelKernel);
}

void ApplyVerticalFilterToPixelAndReturnRGB(double* imageRGB, double** tempStorageHorizConvolution, const double* kernel, int radius ) {
    
    for(int pos = -radius; pos <= radius; pos++) {
        imageRGB[0] += (tempStorageHorizConvolution[pos+radius][0] )* kernel[pos + radius];
        imageRGB[1] += (tempStorageHorizConvolution[pos+radius][1] )* kernel[pos + radius];
        imageRGB[2] += (tempStorageHorizConvolution[pos+radius][2] )* kernel[pos + radius];
    }
}

void ApplyHorizontalFilterToPixelAndReturnRGB(const QImage& outputImageBuffer,double** tempStorageHorizConvolution,const double* kernel, Position position, int radius) {
    QRgb pixel;

    for(int row = -radius; row <= radius; ++row) 
        for(int col = -radius; col <= radius; ++col) {
            pixel           =   outputImageBuffer.pixel(position.YPixelPos + col + radius,position.XPixelPos +row + radius); 
            double  weight  =   kernel[col + radius];
            
            tempStorageHorizConvolution[row+radius][0]     +=  weight*(double) qRed(pixel);
            tempStorageHorizConvolution[row+radius][1]     +=  weight*(double) qGreen(pixel);
            tempStorageHorizConvolution[row+radius][2]     +=  weight*(double) qBlue(pixel);
        }    
}

void ApplySobelFilterToImage(const ImageMetaData& metaData, QImage* inputImage, QImage& outputImageBuffer,SobelKernel& kernel) {
   
    double* imageRGB                        =   new double[3]();
    int radius                              =   metaData.radius;
    double** tempStorageHorizConvolution    =   new double*[kernel.KernelSize]();

    double XKernel[3] = {1, 0, -1};
    double Ykernel[3] = {1, 2, 1};

   

    double Gx, Gy, mag, orient;

    for(int i=0; i<kernel.KernelSize; i++) 
        tempStorageHorizConvolution[i] = new double[3];
    

    for(int rowPixel = 0; rowPixel < metaData.Height; rowPixel++) { // X
        for(int colPixel = 0; colPixel < metaData.Width; colPixel++) { // Y

            ResetValues(imageRGB, tempStorageHorizConvolution, kernel.KernelSize);
            ApplyHorizontalFilterToPixelAndReturnRGB(outputImageBuffer, tempStorageHorizConvolution, XKernel, Position(rowPixel,colPixel),radius);
            ApplyVerticalFilterToPixelAndReturnRGB(imageRGB, tempStorageHorizConvolution, Ykernel, radius);

            Gx = (double) floor(imageRGB[0] + imageRGB[1] + imageRGB[2]);
                Gx = Gx/(double)3;

            ResetValues(imageRGB, tempStorageHorizConvolution, kernel.KernelSize);
            ApplyHorizontalFilterToPixelAndReturnRGB(outputImageBuffer, tempStorageHorizConvolution, Ykernel, Position(rowPixel,colPixel),radius);
            ApplyVerticalFilterToPixelAndReturnRGB(imageRGB, tempStorageHorizConvolution, XKernel, radius);

                Gy = (double) floor(imageRGB[0] + imageRGB[1] + imageRGB[2]);
                Gy = Gy / (double)3;

             
                mag = sqrt(Gx*Gx + Gy*Gy);
                orient = atan2( Gy , Gx );
               

            double red = (sin(orient) + 1.0)/2.0;
             double green = (cos(orient) + 1.0)/2.0;
             double blue = 1.0 - red - green;
             
             red *= mag*4.0;
             green *= mag*4.0;
             blue *= mag*4.0;

    // Make sure the pixel values range from 0 to 255
    red = min(255.0, max(0.0, red));
    green = min(255.0, max(0.0, green));
    blue = min(255.0, max(0.0, blue));

    inputImage->setPixel(colPixel, rowPixel, qRgb( (int) (red), (int) (green), (int) (blue)));

           
        }    } 
    
    delete imageRGB;
    imageRGB = NULL;
    
    for(int i=0; i< kernel.KernelSize ; i++)
        delete tempStorageHorizConvolution[i];
}
void MainWindow::SobelImage(QImage *image)
{

    // Add your code here.
    double XSobelKernel[3][3] = {{1,0,-1},{2,0,-2},{1,0,-1}};
    double YSobelKernel[3][3] = {{1,2,1},{0,0,0},{-1,-2,-1}};

    double sigma = double(1) / double(3);

    ImageMetaData metaData(0,0,0); 
    SobelKernel sobelKernel(0,NULL,NULL);
     GaussianKernel gaussianKernel(0,(double)1/(double)3,NULL);
   
    QImage OutPutImageBuffer        =       InitializeOutputImageBuffer(image,sigma,metaData);
    GetSobelKernel(sobelKernel);
    ApplySobelFilterToImage(metaData, image, OutPutImageBuffer, sobelKernel);

    /***********************************************************************
      When displaying the orientation image I
      recommend the following:

    double mag; // magnitude of the gradient
    double orien; // orientation of the gradient

    double red = (sin(orien) + 1.0)/2.0;
    double green = (cos(orien) + 1.0)/2.0;
    double blue = 1.0 - red - green;

    red *= mag*4.0;
    green *= mag*4.0;
    blue *= mag*4.0;

    // Make sure the pixel values range from 0 to 255
    red = min(255.0, max(0.0, red));
    green = min(255.0, max(0.0, green));
    blue = min(255.0, max(0.0, blue));

    image->setPixel(c, r, qRgb( (int) (red), (int) (green), (int) (blue)));

    ************************************************************************/
}

void GetPixelWeight(QImage* image, const Position& pixelPosition, double* pixelWeight) {

    int edgePixelCol = pixelPosition.YPixelPos;
    int edgePixelRow = pixelPosition.XPixelPos;

    if((edgePixelCol) < 0 || (edgePixelCol)> (image->width()-1))
        edgePixelCol = GetPixelPositionForEdges(edgePixelCol,image->width()-1);
    if((edgePixelRow) < 0 || (edgePixelRow)> (image->height()-1)) 
        edgePixelRow = GetPixelPositionForEdges(edgePixelRow,image->height()-1);
                 

    QRgb pixel =  image->pixel(edgePixelCol, edgePixelRow);
         
    pixelWeight[0] = qRed(pixel);
    pixelWeight[1] = qGreen(pixel);
    pixelWeight[2] = qBlue(pixel);
}

void GetNeighbouringPixelPositionAndWeight(QImage *image, double colPixel, double rowPixel, Position* position, double** pixelWeight) {
    int colPixel1 = (int) (floor(colPixel));
    int rowPixel1 = (int) (floor(rowPixel));

    position[0] = Position(rowPixel1,colPixel1);
    GetPixelWeight(image, Position(rowPixel1,colPixel1), pixelWeight[0]);

    position[1] = Position(rowPixel1,colPixel1+1);
    GetPixelWeight(image, Position(rowPixel1,colPixel1+1), pixelWeight[1]);

    position[2] = Position(rowPixel1+1,colPixel1);
    GetPixelWeight(image, Position(rowPixel1+1,colPixel1), pixelWeight[2]);

    position[3] = Position(rowPixel1+1,colPixel1+1);
    GetPixelWeight(image, Position(rowPixel1+1,colPixel1+1), pixelWeight[3]);
}

double* GetIntermediatePixelValueByInterpolation(double influenceOfPixel1, double influenceOfPixel2, 
                                              double* pixel1Weight, double* pixel2Weight) {
    
    double* imageRGB = new double[3];
    imageRGB[0] = influenceOfPixel1 * pixel2Weight[0] + influenceOfPixel2 * pixel1Weight[0];
    imageRGB[1] = influenceOfPixel1 * pixel2Weight[1] + influenceOfPixel2 * pixel1Weight[1];
    imageRGB[2] = influenceOfPixel1 * pixel2Weight[2] + influenceOfPixel2 * pixel1Weight[2];

    return imageRGB;
}

double* GetPixelValueByInterpolation(double colPixel, double rowPixel, double** pixelWeight) {
    int colPixel1 = (int) (floor(colPixel));
    int rowPixel1 = (int) (floor(rowPixel));

   
    double* imageRGB = new double[3];

    double influenceOfPixelX1 =  (double)(rowPixel) - (double)(rowPixel1);
    double influenceOfPixelY1 =  (double)(colPixel) - (double)(colPixel1);
    double influenceOfPixelX2 =  (double)1 - influenceOfPixelX1;
    double influenceOfPixelY2 =  (double)1 - influenceOfPixelY1;

    double* XY1 = GetIntermediatePixelValueByInterpolation(influenceOfPixelX1, influenceOfPixelX2, pixelWeight[0], pixelWeight[2]);
    double* XY2 = GetIntermediatePixelValueByInterpolation(influenceOfPixelX1, influenceOfPixelX2, pixelWeight[1], pixelWeight[3]);

    double* XY =  GetIntermediatePixelValueByInterpolation(influenceOfPixelY1, influenceOfPixelY2, XY1, XY2);

    return XY;
 
}

void MainWindow::BilinearInterpolation(QImage *image, double colPixel, double rowPixel, double rgb[3])
{
    // Add your code here.  Return the RGB values for the pixel at location (x,y) in double rgb[3].
    Position* position = new Position[4];
    double** pixelWeight = new double*[4](); 
    for(int i=0; i< 4; i++)
        pixelWeight[i] = new double[3];
       
    if( colPixel < 0 || rowPixel < 0 || colPixel > (image->width() -1) || rowPixel > (image->height() -1)) {
        rgb[0] = 0; rgb[1] =0; rgb[2]=0; 
        return;
    }

    GetNeighbouringPixelPositionAndWeight(image, colPixel, rowPixel, position, pixelWeight);
    double* imageRGB = GetPixelValueByInterpolation(colPixel, rowPixel, pixelWeight);

    for(int i=0; i< 3; i++)
        rgb[i] = imageRGB[i];
}

// Here is some sample code for rotating an image.  I assume orien is in degrees.

void MainWindow::RotateImage(QImage *image, double orien)
{
    int r, c;
    QRgb pixel;
    QImage buffer;
    int w = image->width();
    int h = image->height();
    double radians = -2.0*3.141*orien/360.0;

    //double radians = 3.141*orien/180.0;

    buffer = image->copy();

    pixel = qRgb(0, 0, 0);
    image->fill(pixel);

    for(r=0;r<h;r++)
    {
        std::cout << "r = " << r << std::endl;
        for(c=0;c<w;c++)
        {
            double rgb[3];
            double x0, y0;
            double x1, y1;

            // Rotate around the center of the image.
            x0 = (double) (c - w/2); // x = column
            y0 = (double) (r - h/2); // y = row

            // Rotate using rotation matrix
            x1 = x0*cos(radians) - y0*sin(radians);
            y1 = x0*sin(radians) + y0*cos(radians);

            x1 += (double) (w/2); // x1 = colPixel
            y1 += (double) (h/2); // y1 = rowPixel

            BilinearInterpolation(&buffer, x1, y1, rgb);

            image->setPixel(c, r, qRgb((int) floor(rgb[0] + 0.5), (int) floor(rgb[1] + 0.5), (int) floor(rgb[2] + 0.5)));
        }
    }

}



void GetNonMaximumSupression(QImage* image, double row, double col, int rowPixel, int colPixel, double* magnitude, double thres ) {

       // Add your code here.  Return the RGB values for the pixel at location (x,y) in double rgb[3].
    Position* position = new Position[4];
    double rgb[3] = {0.0};
    double** pixelWeight = new double*[4](); 
    for(int i=0; i< 4; i++)
        pixelWeight[i] = new double[3];

     double* E0imageRGB;
      double* E1imageRGB;
       
    if( colPixel < 0 || rowPixel < 0 || colPixel > (image->width() -1) || rowPixel > (image->height() -1)) {
        rgb[0] = 0; rgb[1] =0; rgb[2]=0; 
        return;
    }

    GetNeighbouringPixelPositionAndWeight(image, colPixel + col, rowPixel + row, position, pixelWeight);
    E0imageRGB = GetPixelValueByInterpolation(colPixel + col, rowPixel + row, pixelWeight);

     GetNeighbouringPixelPositionAndWeight(image, colPixel -col, rowPixel -row, position, pixelWeight);
     E1imageRGB = GetPixelValueByInterpolation(colPixel-col, rowPixel-row, pixelWeight);

     double E0magnitude = (E0imageRGB[0]+E0imageRGB[1]+E0imageRGB[2])/3;
    double E1magnitude = (E1imageRGB[0]+E1imageRGB[1]+E1imageRGB[2])/3;
   
    if((*magnitude > E0magnitude) && (*magnitude > E1magnitude) && (*magnitude > thres))
        *magnitude = 255;
    else 
        *magnitude = 0;
       
}

void findPeaksSobel(const ImageMetaData& metaData, QImage* inputImage, QImage& outputImageBuffer,SobelKernel& kernel,double thres) {
     double* imageRGB                        =   new double[3]();
    int radius                              =   1;
    double** tempStorageHorizConvolution    =   new double*[kernel.KernelSize](); 

    double XKernel[3] = {1, 0, -1};
    double Ykernel[3] = {1, 2, 1};

    double Gx, Gy; 
    double** magnitude = new double*[metaData.Height];
    int** orientation = new int*[metaData.Height];

    double mag, orient;

    for(int i=0; i< metaData.Height; i++){
        magnitude[i] = new double[metaData.Width];
        orientation[i] = new int[metaData.Width];
    }

    for(int i=0; i<kernel.KernelSize; i++) 
        tempStorageHorizConvolution[i] = new double[3];
    
    for(int rowPixel = 0; rowPixel < metaData.Height; rowPixel++) { // X
        for(int colPixel = 0; colPixel < metaData.Width; colPixel++) { // Y

            double red = 0.0;
             double green = 0.0;
             double blue = 0.0;

            ResetValues(imageRGB, tempStorageHorizConvolution, kernel.KernelSize);
            ApplyHorizontalFilterToPixelAndReturnRGB(outputImageBuffer, tempStorageHorizConvolution, XKernel, Position(rowPixel,colPixel),radius);
            ApplyVerticalFilterToPixelAndReturnRGB(imageRGB, tempStorageHorizConvolution, Ykernel, radius);

            Gx = (int) floor(imageRGB[0] + imageRGB[1] + imageRGB[2]);
           /* red = abs(imageRGB[0]);
            green = abs(imageRGB[1]);
            blue = abs(imageRGB[2]);*/

            Gx = Gx/3;

            ResetValues(imageRGB, tempStorageHorizConvolution, kernel.KernelSize);
            ApplyHorizontalFilterToPixelAndReturnRGB(outputImageBuffer, tempStorageHorizConvolution, Ykernel, Position(rowPixel,colPixel),radius);
            ApplyVerticalFilterToPixelAndReturnRGB(imageRGB, tempStorageHorizConvolution, XKernel, radius);

           /* red += abs(imageRGB[0]);
            green += abs(imageRGB[1]);
            blue += abs(imageRGB[2]);*/

            Gy = (int) floor(imageRGB[0] + imageRGB[1] + imageRGB[2]);
            Gy = Gy / 3;

             
                
            mag = sqrt(Gx*Gx + Gy*Gy);
            magnitude[rowPixel][colPixel] = mag;

            orient = atan2( Gy , Gx );

            double value = Gx + Gy;
           /* double red = (sin(orient) + 1.0)/2.0;
            double green = (cos(orient) + 1.0)/2.0;
            double blue = 1.0 - red - green;*/
             
            /*red *= mag*4.0;
            green *= mag*4.0;
            blue *= mag*4.0;*/

            // Make sure the pixel values range from 0 to 255
            /*red = min(255.0, max(0.0, red));
            green = min(255.0, max(0.0, green));
            blue = min(255.0, max(0.0, blue));*/

            orient = (int)(180 + (orient / M_PI * 180.0));   

            if(orient >= 0 && orient <=45)
                orientation[rowPixel][colPixel] = 0;
            else if(orient > 45 && orient <= 90)
                orientation[rowPixel][colPixel] = 90;
            else if(orient > 90 && orient <= 135)
                orientation[rowPixel][colPixel] = 135;
            else if(orient > 135 && orient <= 180)
                orientation[rowPixel][colPixel] = 180;

             //double red = mag;//*4.0;
             //double green = mag;//*4.0;
             //double blue = mag;//*4.0;

    // Make sure the pixel values range from 0 to 255
            red = min(255.0, max(0.0, red));
            green = min(255.0, max(0.0, green));
            blue = min(255.0, max(0.0, blue));

    //inputImage->setPixel(colPixel, rowPixel, qRgb( (int) (red), (int) (green), (int) (blue)));
     inputImage->setPixel(colPixel, rowPixel, qRgb( (int) (value), (int) (value), (int) (value)));


        }    
}
    delete imageRGB;
    imageRGB = NULL;
    
    for(int i=0; i< kernel.KernelSize ; i++)
        delete tempStorageHorizConvolution[i];

    bool isEdge = false;

     for(int rowPixel = 0; rowPixel < metaData.Height; rowPixel++) { // X
        for(int colPixel = 0; colPixel < metaData.Width; colPixel++) { // Y

            switch(orientation[rowPixel][colPixel]) {
            case 0 :
                GetNonMaximumSupression(inputImage,0,1, rowPixel, colPixel, &magnitude[rowPixel][colPixel], thres);
                break;
            case 45 :
                GetNonMaximumSupression(inputImage,(double)-1/(double)sqrtf(2),(double)1/(double)sqrtf(2), rowPixel, colPixel,&magnitude[rowPixel][colPixel], thres);
                break;
            case 90:
                GetNonMaximumSupression(inputImage,1,0, rowPixel, colPixel, &magnitude[rowPixel][colPixel], thres);
                break;
            case 135:
                GetNonMaximumSupression(inputImage,(double)-1/(double)sqrtf(2),(double)-1/(double)sqrtf(2), rowPixel, colPixel, &magnitude[rowPixel][colPixel], thres);
                break;
            }
            
        }
     }

      for(int rowPixel = 0; rowPixel < metaData.Height; rowPixel++) { // X
        for(int colPixel = 0; colPixel < metaData.Width; colPixel++) { // Y
            int value = magnitude[rowPixel][colPixel];
            inputImage->setPixel(colPixel, rowPixel, qRgb(value, value, value));
        }
      }
}

void MainWindow::FindPeaksImage(QImage *image, double thres)
{
     double XSobelKernel[3][3] = {{1,0,-1},{2,0,-2},{1,0,-1}};
    double YSobelKernel[3][3] = {{1,2,1},{0,0,0},{-1,-2,-1}};

    double sigma = double(1) / double(3);

    ImageMetaData metaData(0,0,0); 
    SobelKernel sobelKernel(0,NULL,NULL);
     GaussianKernel gaussianKernel(0,(double)1/(double)3,NULL);
   
    QImage OutPutImageBuffer        =       InitializeOutputImageBuffer(image,sigma,metaData);
    GetSobelKernel(sobelKernel);
    //GaussianBlurImage(image,2);
    ApplySobelFilterToImage(metaData, image, OutPutImageBuffer, sobelKernel);
    OutPutImageBuffer        =       InitializeOutputImageBuffer(image,sigma,metaData);
    findPeaksSobel(metaData,image,OutPutImageBuffer, sobelKernel,thres);
    // Add your code here.
}

void MainWindow::MedianImage(QImage *image, int radius)
{
    // Add your code here
}

void MainWindow::HoughImage(QImage *image)
{
    // Add your code here
}

void MainWindow::CrazyImage(QImage *image)
{
    // Add your code here
}
