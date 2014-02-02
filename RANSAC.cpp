#include "mainwindow.h"
#include "math.h"
#include "ui_mainwindow.h"
#include <QtGui>
#include "Matrix.h"

/*******************************************************************************
    The following are helper routines with code already written.
    The routines you'll need to write for the assignment are below.
*******************************************************************************/

/*******************************************************************************
Draw detected Harris corners
    interestPts - interest points
    numInterestsPts - number of interest points
    imageDisplay - image used for drawing

    Draws a red cross on top of detected corners
*******************************************************************************/
void MainWindow::DrawInterestPoints(CIntPt *interestPts, int numInterestsPts, QImage &imageDisplay)
{
   int i;
   int r, c, rd, cd;
   int w = imageDisplay.width();
   int h = imageDisplay.height();

   for(i=0;i<numInterestsPts;i++)
   {
       c = (int) interestPts[i].m_X;
       r = (int) interestPts[i].m_Y;

       for(rd=-2;rd<=2;rd++)
           if(r+rd >= 0 && r+rd < h && c >= 0 && c < w)
               imageDisplay.setPixel(c, r + rd, qRgb(255, 0, 0));

       for(cd=-2;cd<=2;cd++)
           if(r >= 0 && r < h && c + cd >= 0 && c + cd < w)
               imageDisplay.setPixel(c + cd, r, qRgb(255, 0, 0));
   }
}

/*******************************************************************************
Compute interest point descriptors
    image - input image
    interestPts - array of interest points
    numInterestsPts - number of interest points

    If the descriptor cannot be computed, i.e. it's too close to the boundary of
    the image, its descriptor length will be set to 0.

    I've implemented a very simple 8 dimensional descriptor.  Feel free to
    improve upon this.
*******************************************************************************/
void MainWindow::ComputeDescriptors(QImage image, CIntPt *interestPts, int numInterestsPts)
{
    int r, c, cd, rd, i, j;
    int w = image.width();
    int h = image.height();
    double *buffer = new double [w*h];
    QRgb pixel;

    // Descriptor parameters
    double sigma = 2.0;
    int rad = 4;

    // Computer descriptors from green channel
    for(r=0;r<h;r++)
       for(c=0;c<w;c++)
        {
            pixel = image.pixel(c, r);
            buffer[r*w + c] = (double) qGreen(pixel);
        }

    // Blur
    SeparableGaussianBlurImage(buffer, w, h, sigma);

    // Compute the desciptor from the difference between the point sampled at its center
    // and eight points sampled around it.
    for(i=0;i<numInterestsPts;i++)
    {
        int c = (int) interestPts[i].m_X;
        int r = (int) interestPts[i].m_Y;

        if(c >= rad && c < w - rad && r >= rad && r < h - rad)
        {
            double centerValue = buffer[(r)*w + c];
            int j = 0;

            for(rd=-1;rd<=1;rd++)
                for(cd=-1;cd<=1;cd++)
                    if(rd != 0 || cd != 0)
                {
                    interestPts[i].m_Desc[j] = buffer[(r + rd*rad)*w + c + cd*rad] - centerValue;
                    j++;
                }

            interestPts[i].m_DescSize = DESC_SIZE;
        }
        else
        {
            interestPts[i].m_DescSize = 0;
        }
    }

    delete [] buffer;
}

/*******************************************************************************
Draw matches between images
    matches - matching points
    numMatches - number of matching points
    image1Display - image to draw matches
    image2Display - image to draw matches

    Draws a green line between matches
*******************************************************************************/
void MainWindow::DrawMatches(CMatches *matches, int numMatches, QImage &image1Display, QImage &image2Display)
{
    int i;
    // Show matches on image
    QPainter painter;
    painter.begin(&image1Display);
    QColor green(0, 250, 0);
    QColor red(250, 0, 0);

    for(i=0;i<numMatches;i++)
    {
        painter.setPen(green);
        painter.drawLine((int) matches[i].m_X1, (int) matches[i].m_Y1, (int) matches[i].m_X2, (int) matches[i].m_Y2);
        painter.setPen(red);
        painter.drawEllipse((int) matches[i].m_X1-1, (int) matches[i].m_Y1-1, 3, 3);
    }

    QPainter painter2;
    painter2.begin(&image2Display);
    painter2.setPen(green);

    for(i=0;i<numMatches;i++)
    {
        painter2.setPen(green);
        painter2.drawLine((int) matches[i].m_X1, (int) matches[i].m_Y1, (int) matches[i].m_X2, (int) matches[i].m_Y2);
        painter2.setPen(red);
        painter2.drawEllipse((int) matches[i].m_X2-1, (int) matches[i].m_Y2-1, 3, 3);
    }

}


/*******************************************************************************
Given a set of matches computes the "best fitting" homography
    matches - matching points
    numMatches - number of matching points
    h - returned homography
    isForward - direction of the projection (true = image1 -> image2, false = image2 -> image1)
*******************************************************************************/
bool MainWindow::ComputeHomography(CMatches *matches, int numMatches, double h[3][3], bool isForward)
{
    int error;
    int nEq=numMatches*2;

    dmat M=newdmat(0,nEq,0,7,&error);
    dmat a=newdmat(0,7,0,0,&error);
    dmat b=newdmat(0,nEq,0,0,&error);

    double x0, y0, x1, y1;

    for (int i=0;i<nEq/2;i++)
    {
        if(isForward == false)
        {
            x0 = matches[i].m_X1;
            y0 = matches[i].m_Y1;
            x1 = matches[i].m_X2;
            y1 = matches[i].m_Y2;
        }
        else
        {
            x0 = matches[i].m_X2;
            y0 = matches[i].m_Y2;
            x1 = matches[i].m_X1;
            y1 = matches[i].m_Y1;
        }


        //Eq 1 for corrpoint
        M.el[i*2][0]=x1;
        M.el[i*2][1]=y1;
        M.el[i*2][2]=1;
        M.el[i*2][3]=0;
        M.el[i*2][4]=0;
        M.el[i*2][5]=0;
        M.el[i*2][6]=(x1*x0*-1);
        M.el[i*2][7]=(y1*x0*-1);

        b.el[i*2][0]=x0;
        //Eq 2 for corrpoint
        M.el[i*2+1][0]=0;
        M.el[i*2+1][1]=0;
        M.el[i*2+1][2]=0;
        M.el[i*2+1][3]=x1;
        M.el[i*2+1][4]=y1;
        M.el[i*2+1][5]=1;
        M.el[i*2+1][6]=(x1*y0*-1);
        M.el[i*2+1][7]=(y1*y0*-1);

        b.el[i*2+1][0]=y0;

    }
    int ret=solve_system (M,a,b);
    if (ret!=0)
    {
        freemat(M);
        freemat(a);
        freemat(b);

        return false;
    }
    else
    {
        h[0][0]= a.el[0][0];
        h[0][1]= a.el[1][0];
        h[0][2]= a.el[2][0];

        h[1][0]= a.el[3][0];
        h[1][1]= a.el[4][0];
        h[1][2]= a.el[5][0];

        h[2][0]= a.el[6][0];
        h[2][1]= a.el[7][0];
        h[2][2]= 1;
    }

    freemat(M);
    freemat(a);
    freemat(b);

    return true;
}


/*******************************************************************************
*******************************************************************************
*******************************************************************************

    The routines you need to implement are below

*******************************************************************************
*******************************************************************************
*******************************************************************************/


/*******************************************************************************
Blur a single channel floating point image with a Gaussian.
    image - input and output image
    w - image width
    h - image height
    sigma - standard deviation of Gaussian

    This code should be very similar to the code you wrote for assignment 1.
*******************************************************************************/

#pragma region UTILITY

struct ImageMetaData {
    int Width;
    int Height;
    int radius;
    ImageMetaData():Width(0),Height(0),radius(0){}
    ImageMetaData(int width_, int height_, int radius_):Width(width_),Height(height_),radius(radius_){}
    ImageMetaData(const ImageMetaData& metaData):Width(metaData.Width),Height(metaData.Height),radius(metaData.radius){}
};

struct Position {
    int XPixelPos;
    int YPixelPos;

    Position():XPixelPos(-1), YPixelPos(-1){}
    Position(int X, int Y):XPixelPos(X), YPixelPos(Y){}
    Position(const Position& position_):XPixelPos(position_.XPixelPos), YPixelPos(position_.YPixelPos){}
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

int GetPixelPositionForEdges(int currentPixel, int edgesize) {
    if(currentPixel < 0)
        return abs(currentPixel);
    else if(currentPixel > (edgesize))
        return abs(edgesize-(currentPixel - edgesize));
    else
        return edgesize;
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
#pragma endregion UTILITY

#pragma region SEPERABLEGAUSSIANBLURIMAGE 

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

void  ApplyVerticalFilterToPixel(double* inputImage, double* tempImage,const ImageMetaData& metaData,const GaussianKernel& kernel ) {
    double pixel;
    int edgePixel;
    double imageRGB =0.0;

    int radius = metaData.radius;

    for(int rowPixel = 0; rowPixel < metaData.Height; rowPixel++) { // X
        for(int colPixel = 0; colPixel < metaData.Width; colPixel++) { // Y
            
            imageRGB = 0.0;
             for(int row = -radius; row <= radius; ++row) {

                 if((rowPixel+row) < 0 || (rowPixel+row)> (metaData.Height-1)) {
                     edgePixel = GetPixelPositionForEdges(rowPixel+row,metaData.Height-1);
                     pixel = tempImage[edgePixel * metaData.Width + colPixel];
                 }
                 else
                     pixel = tempImage[(rowPixel+row) * metaData.Width + colPixel];
                 double weight = kernel.Kernel[row+radius];
                 imageRGB += pixel * weight;
             }
             //qDebug() << imageRGB;
              //imageRGB = min(255.0, max(0.0, imageRGB));
             inputImage[rowPixel * metaData.Width + colPixel] = imageRGB;
        }
    }
}

void ApplyHorizontalFilterToPixel(double* inputImage,double* tempImage, const ImageMetaData& metaData,const GaussianKernel& kernel ) {
    double pixel;
    int edgePixel;
    double imageRGB = 0.0;
    int radius = metaData.radius;

    for(int rowPixel = 0; rowPixel < metaData.Height; rowPixel++) { // X
        for(int colPixel = 0; colPixel < metaData.Width; colPixel++) { // Y
            imageRGB = 0.0;
            for(int col = -radius; col <= radius; ++col) {
                if((colPixel+col) < 0 || (colPixel+col)> (metaData.Width-1)) {
                    edgePixel = GetPixelPositionForEdges(colPixel+col,metaData.Width-1);
                    pixel = inputImage[rowPixel * metaData.Width + edgePixel];
                }
                else 
                    pixel = inputImage[rowPixel * metaData.Width + colPixel +col];

                 double weight = kernel.Kernel[col+radius];
                 imageRGB  += (double)(pixel) * weight;
             }

            //qDebug() << imageRGB;
           //imageRGB = min(255.0, max(0.0, imageRGB));
            tempImage[rowPixel * metaData.Width + colPixel] = imageRGB;
        }
    }
}

void ConvolveImagewithSeperableGaussianKernel(const ImageMetaData& metaData, double* inputImage,GaussianKernel& kernel) {
     
    int sizeOfTempImage = (metaData.Width )* (metaData.Height );
    double* tempImage = new double[sizeOfTempImage];

    ApplyHorizontalFilterToPixel(inputImage, tempImage, metaData,kernel);
    ApplyVerticalFilterToPixel(inputImage,tempImage,metaData, kernel);
}

void MainWindow::SeparableGaussianBlurImage(double *image, int w, int h, double sigma)
{
    // Add your code here

    // To access the pixel (c,r), use image[r*width + c].

    ImageMetaData metaData(w,h,3*sigma);
    GaussianKernel gaussianKernel(0, sigma, NULL);
    GetSeperableGaussianKernel(metaData.radius, gaussianKernel);
    ConvolveImagewithSeperableGaussianKernel(metaData, image, gaussianKernel);
}

#pragma endregion SEPERABLEGAUSSIANBLURIMAGEE
/*******************************************************************************
Detect Harris corners.
    image - input image
    sigma - standard deviation of Gaussian used to blur corner detector
    thres - Threshold for detecting corners
    interestPts - returned interest points
    numInterestsPts - number of interest points returned
    imageDisplay - image returned to display (for debugging)
*******************************************************************************/
/*******COVARIANCE MATRIX *********************************/
/*   [XX   XY]   */
/*   [XY   YY]   */
/****************************************/


int ComputeHarrisResponse(const ImageMetaData& metaData,double* XDerivativeBuffer, double* YDerivativeBuffer, double* XYDerivativeBuffer,
                           double thres, CIntPt **interestPts) {

    double XX, YY, XY;
    double* determinantOfCovarianceMatrix = new double[metaData.Width * metaData.Height];
    double* traceCovarianceMatrix         = new double[metaData.Width * metaData.Height];
    double* harrisResponse                = new double[metaData.Width * metaData.Height];

    double epsilon = 0.00001;
    int numInterestsPts = 0;
    double sumOfValues = 0.0;
    double minValue = (double)MAXINT;
    double maxValue = (double)MININT;

    double** ValueOfHarrisResponse = new double*[metaData.Height];
    for(int i = 0; i < metaData.Height; i++)
        ValueOfHarrisResponse[i] = new double[metaData.Width];

    std::vector<int> interestPtsPosition;
  

    for(int rowPixel = 0; rowPixel < metaData.Height; rowPixel++) 
        for(int colPixel = 0; colPixel < metaData.Width; colPixel++) {
            int pixelPosition = rowPixel * metaData.Width + colPixel;

            XX = XDerivativeBuffer[pixelPosition];
            YY = YDerivativeBuffer[pixelPosition];
            XY = XYDerivativeBuffer[pixelPosition];

            determinantOfCovarianceMatrix[pixelPosition] = (XX * YY) - (XY * XY);
            traceCovarianceMatrix[pixelPosition] = XX + YY;
            harrisResponse[pixelPosition] =  determinantOfCovarianceMatrix[pixelPosition] / (traceCovarianceMatrix[pixelPosition] + epsilon);

            if(minValue > harrisResponse[pixelPosition])
                minValue = harrisResponse[pixelPosition];
            if(maxValue < harrisResponse[pixelPosition])
                maxValue = harrisResponse[pixelPosition];

            //sumOfValues += harrisResponse[pixelPosition];
        }

        // Normalize values

        //double ValueOfHarrisResponse = 0.0;
        double differenceOfValue = maxValue - minValue;

        std::vector<Position>XYPosition;

        for(int rowPixel = 0; rowPixel < metaData.Height; rowPixel++) {
            for(int colPixel = 0; colPixel < metaData.Width; colPixel++) {
                int pixelPosition = rowPixel * metaData.Width + colPixel;
            //ValueOfHarrisResponse = 0.0;
            //ValueOfHarrisResponse = (((255) * ( harrisResponse[pixelPosition] - minValue)) / differenceOfValue) + minValue;

                ValueOfHarrisResponse[rowPixel][colPixel] = (((255) * ( harrisResponse[pixelPosition] - minValue)) / differenceOfValue) + minValue;
            }
        }

   
        // Compute Non maximum supression
/*
        for(int rowPixel = 0; rowPixel < metaData.Height; rowPixel++) {
            for(int colPixel = 0; colPixel < metaData.Width; colPixel++) {
*/
        // Leaving the first/last row and first/last column
         for(int rowPixel = 1; rowPixel < metaData.Height-1; rowPixel++) {
            for(int colPixel = 1; colPixel < metaData.Width-1; colPixel++) {
               
                if(ValueOfHarrisResponse[rowPixel][colPixel] > thres) {

                    double val = ValueOfHarrisResponse[rowPixel][colPixel];

                    if( (val > ValueOfHarrisResponse[rowPixel-1][colPixel-1]) &&
                        (val > ValueOfHarrisResponse[rowPixel-1][colPixel]) &&
                        (val > ValueOfHarrisResponse[rowPixel-1][colPixel+1]) &&
                        (val > ValueOfHarrisResponse[rowPixel][colPixel-1]) &&
                        (val > ValueOfHarrisResponse[rowPixel][colPixel+1]) &&
                        (val > ValueOfHarrisResponse[rowPixel+1][colPixel-1]) &&
                        (val > ValueOfHarrisResponse[rowPixel+1][colPixel]) &&
                        (val > ValueOfHarrisResponse[rowPixel+1][colPixel+1])) {
                            
                            XYPosition.push_back(Position(rowPixel, colPixel));
                            numInterestsPts++;
                    }}
                    // check around the 8 pixels surrounding it.
                }
            }
        




            //harrisResponse[pixelPosition] /= (double)sumOfValues;

            //float f2 = max(0.0, min(255.0, harrisResponse[pixelPosition]*255.0));
            //float b = floor(f2 == 1.0 ? 255 : f2 * 256.0);

   //         if(ValueOfHarrisResponse > thres) {
   //         //if(f2 > thres) {
   //             interestPtsPosition.push_back(pixelPosition);
   //             XYPosition.push_back(Position(rowPixel, colPixel));
   //             numInterestsPts++;
   //         }
   //     }
   //     
   *interestPts = new CIntPt [numInterestsPts];
   int count = 0;
   //// Access the values using: (*interestPts)[i].m_X = 5.0;

   /*for(int rowPixel = 0; rowPixel < metaData.Height; rowPixel++) 
       for(int colPixel = 0; colPixel < metaData.Width; colPixel++) {
           int pixelPosition = rowPixel * metaData.Width + colPixel;
           if(harrisResponse[pixelPosition] > thres) {
               (*interestPts)[count].m_Y = rowPixel;
               (*interestPts)[count].m_X = colPixel;

               count++;
           }*/
   ////           //Find peaks in the response that are above the threshold "thres", and store the interest point locations in "interestPts. "
   ////     }

   for(int num = 0; num < numInterestsPts; num++) {
       (*interestPts)[num].m_Y = XYPosition.at(num).XPixelPos;
       (*interestPts)[num].m_X = XYPosition.at(num).YPixelPos;
   }

       return numInterestsPts;
}

void ComputeSumOfDerivativesWithGaussianWindow(const ImageMetaData& metaData,double* XDerivativeBuffer, 
                                  double* YDerivativeBuffer, double* XYDerivativeBuffer, double sigma) {

    double* sumOfXDerivative = new double[metaData.Width * metaData.Height];
    double* sumOfYDerivative = new double[metaData.Width * metaData.Height];
    double* sumOfXYDerivative = new double[metaData.Width * metaData.Height];

    //SeparableGaussianBlurImage(XDerivativeBuffer, metaData.Width, metaData.Height, sigma);



}

void ApplyDerivativeFilterToImage(const ImageMetaData& metaData, double* inputImage, double* XDerivativeBuffer, 
                                  double* YDerivativeBuffer, double* XYDerivativeBuffer) {
   
    double* imageRGB                        =   new double[3]();
    int radius                              =   metaData.radius;
    double Gx, Gy, mag, orient;

    double XKernel[3][3] = {{-1, 0, 1}, {-1, 0, 1},{-1, 0, 1}};
    double YKernel[3][3] = {{-1, -1, -1}, {0, 0, 0},{1, 1, 1}};

//    double *XDerivativeBuffer = new double[metaData.Width * metaData.Height];
 //   double *YDerivativeBuffer = new double[metaData.Width * metaData.Height];
 //   double *XYDerivativeBuffer = new double[metaData.Width * metaData.Height];

//    double* tempImage = new double[metaData.Width * metaData.Height];
    double pixel;

    int rowEdgePixel, colEdgePixel;
    
    for(int rowPixel = 0; rowPixel < metaData.Height; rowPixel++) { // X
        for(int colPixel = 0; colPixel < metaData.Width; colPixel++) { // Y
            Gx = 0.0; Gy = 0.0;

            for(int row = -1; row < 2; row++) {
                rowEdgePixel = rowPixel + row;
                if((rowEdgePixel) < 0 || (rowEdgePixel) > (metaData.Height-1))
                    rowEdgePixel = GetPixelPositionForEdges(rowEdgePixel,metaData.Height-1);
                
                for(int col = -1; col < 2; col++) {
                    colEdgePixel = colPixel + col;              
                    if((colEdgePixel) < 0 || (colEdgePixel) > (metaData.Width-1)) 
                        colEdgePixel = GetPixelPositionForEdges(colEdgePixel,metaData.Width-1);

                    pixel = inputImage[rowEdgePixel* metaData.Width + colEdgePixel] ;
                    
                    double weight = XKernel[row+1][col+1];
                    Gx += (double)(pixel) * weight;
                    
                    weight = YKernel[row+1][col+1];
                    Gy += (double)(pixel) * weight;
                }
            }

            //qDebug() << Gx << Gy ;
            XDerivativeBuffer[rowPixel * metaData.Width + colPixel] = Gx * Gx;
            YDerivativeBuffer[rowPixel * metaData.Width + colPixel] = Gy * Gy;
            XYDerivativeBuffer[rowPixel * metaData.Width + colPixel] = Gx * Gy;
        }
      }           
}

void MainWindow::HarrisCornerDetector(QImage image, double sigma, double thres, CIntPt **interestPts, int &numInterestsPts, QImage &imageDisplay)
{
    int rowPixel, colPixel;
    int width = image.width();
    int height = image.height();
    double *buffer = new double [width*height];
    QRgb pixel;

    numInterestsPts = 0;

    // Compute the corner response using just the green channel
    for(rowPixel = 0; rowPixel < height; rowPixel++)
       for(colPixel =0; colPixel <width; colPixel++) {
           pixel = image.pixel(colPixel, rowPixel);
           buffer[rowPixel * width + colPixel] = (double) qGreen(pixel);
       }

    // Write your Harris corner detection code here.

    // Apply the X and Y derivative of Image.
       ImageMetaData metaData(width,height, 3*sigma);

       double *XDerivativeBuffer = new double[width * height];
       double *YDerivativeBuffer = new double[width * height];
       double *XYDerivativeBuffer = new double[width * height];

       ApplyDerivativeFilterToImage(metaData, buffer, XDerivativeBuffer, YDerivativeBuffer, XYDerivativeBuffer);

       // ApplyGaussian filter
        SeparableGaussianBlurImage(XDerivativeBuffer, metaData.Width, metaData.Height, sigma);
        SeparableGaussianBlurImage(YDerivativeBuffer, metaData.Width, metaData.Height, sigma);
        SeparableGaussianBlurImage(XYDerivativeBuffer, metaData.Width, metaData.Height, sigma);

        int numpts = ComputeHarrisResponse(metaData, XDerivativeBuffer, YDerivativeBuffer, XYDerivativeBuffer,thres, interestPts);
       // change it to accomodate first and last row / first and last column

      /* for(rowPixel = 1; rowPixel < height-1; rowPixel++) {
           int currentRow = rowPixel * width ;
           for(colPixel = 1; colPixel < width-1; colPixel++) {

               XDerivativeBuffer[ currentRow + colPixel] =  buffer[currentRow + colPixel + 1] - buffer[currentRow + colPixel - 1] ;
               YDerivativeBuffer[ currentRow + colPixel] =  buffer[currentRow + width + colPixel] - buffer[currentRow -width + colPixel] ;
               XYDerivativeBuffer[currentRow + colPixel] =  XDerivativeBuffer [ currentRow + colPixel] * YDerivativeBuffer [ currentRow + colPixel];
           }
       }*/

    // Once you uknow the number of interest points allocate an array as follows:
    // *interestPts = new CIntPt [numInterestsPts];
    // Access the values using: (*interestPts)[i].m_X = 5.0;
    //
    // The position of the interest point is (m_X, m_Y)
    // The descriptor of the interest point is stored in m_Desc
    // The length of the descriptor is m_DescSize, if m_DescSize = 0, then it is not valid.

    // Once you are done finding the interest points, display them on the image

        numInterestsPts = numpts;
    DrawInterestPoints(*interestPts, numInterestsPts, imageDisplay);

    //delete [] buffer;
}



void FindClosestDescriptor(CIntPt *interestPts1, int numInterestsPts1,CIntPt *interestPts2, int numInterestsPts2, CMatches **matches) {

    // Using L1 Norm

    int numMatches = MAXINT;
    double l1norm = 0.0;
   // int* matchingPoints = new int[numInterestsPts1];
    int matchingPoint = -1;

    int numberOfMatching = 0;

     // Once you uknow the number of matches allocate an array as follows:
     *matches = new CMatches [numInterestsPts1];
    //
    // The position of the interest point in iamge 1 is (m_X1, m_Y1)
    // The position of the interest point in image 2 is (m_X2, m_Y2)


   for(int num1 = 0; num1 < numInterestsPts1 ; num1++) {
        numMatches = MAXINT;
       

       for(int num2 = 0; num2 < numInterestsPts2; num2++) {
            l1norm = 0.0;
           for(int descPts = 0; descPts < 8; descPts++)
               l1norm += abs(interestPts1[num1].m_Desc[descPts] - interestPts2[num2].m_Desc[descPts]);

           if(l1norm < numMatches) {
               numMatches = l1norm;
               matchingPoint = num2;
              // matchingPoints[num1] = num2;
           }
       }

       (*matches)[num1].m_X1 = interestPts1[num1].m_X;
       (*matches)[num1].m_Y1 = interestPts1[num1].m_Y;
       (*matches)[num1].m_X2 = interestPts2[matchingPoint].m_X;//[matchingPoints[num1]].m_X;
       (*matches)[num1].m_Y2 = interestPts2[matchingPoint].m_Y;//[matchingPoints[num1]].m_Y;
   }

   
}

/*******************************************************************************
Find matching interest points between images.
    image1 - first input image
    interestPts1 - interest points corresponding to image 1
    numInterestsPts1 - number of interest points in image 1
    image2 - second input image
    interestPts2 - interest points corresponding to image 2
    numInterestsPts2 - number of interest points in image 2
    matches - set of matching points to be returned
    numMatches - number of matching points returned
    image1Display - image used to display matches
    image2Display - image used to display matches
*******************************************************************************/
void MainWindow::MatchInterestPoints(QImage image1, CIntPt *interestPts1, int numInterestsPts1,
                             QImage image2, CIntPt *interestPts2, int numInterestsPts2,
                             CMatches **matches, int &numMatches, QImage &image1Display, QImage &image2Display)
{
    numMatches = 0;

    // Compute the descriptors for each interest point.
    // You can access the descriptor for each interest point using interestPts1[i].m_Desc[j].
    // If interestPts1[i].m_DescSize = 0, it was not able to compute a descriptor for that point
    ComputeDescriptors(image1, interestPts1, numInterestsPts1);
    ComputeDescriptors(image2, interestPts2, numInterestsPts2);

    // Add your code here for finding the best matches for each point.

    // Find Closest Descriptor

    FindClosestDescriptor(interestPts1, numInterestsPts1,  interestPts2, numInterestsPts2, matches);
    // Once you uknow the number of matches allocate an array as follows:
    // *matches = new CMatches [numMatches];
    //
    // The position of the interest point in iamge 1 is (m_X1, m_Y1)
    // The position of the interest point in image 2 is (m_X2, m_Y2)

    numMatches = numInterestsPts1;
    // Draw the matches
    DrawMatches(*matches, numMatches, image1Display, image2Display);
}

/*******************************************************************************
Project a point (x1, y1) using the homography transformation h
    (x1, y1) - input point
    (x2, y2) - returned point
    h - input homography used to project point
*******************************************************************************/
void MainWindow::Project(double x1, double y1, double &x2, double &y2, double h[3][3])
{
    // Add your code here.

    x2 = (x1*h[0][0] + y1*h[0][1] + h[0][2])/(x1*h[2][0] + y1*h[2][1] + h[2][2]);
    y2 = (x1*h[1][0] + y1*h[1][1] + h[1][2])/(x1*h[2][0] + y1*h[2][1] + h[2][2]);
}

/*******************************************************************************
Count the number of inliers given a homography.  This is a helper function for RANSAC.
    h - input homography used to project points (image1 -> image2
    matches - array of matching points
    numMatches - number of matchs in the array
    inlierThreshold - maximum distance between points that are considered to be inliers

    Returns the total number of inliers.
*******************************************************************************/
int MainWindow::ComputeInlierCount(double h[3][3], CMatches *matches, int numMatches, double inlierThreshold)
{
    // Add your code here.
    double projectedX1, projectedY1;
    int totalInliers = 0;
    double distance = 0.0;

    for(int num = 0; num < numMatches; num++) {
        Project(matches[num].m_X1, matches[num].m_Y1, projectedX1, projectedY1, h);
        double x1diff = projectedX1 - matches[num].m_X2;
        double y1diff = projectedY1 - matches[num].m_Y2;
        distance = sqrt((x1diff * x1diff) + (y1diff * y1diff));

        if(distance < inlierThreshold)
            totalInliers++;
    }

    return totalInliers;
}

template <class T>
bool contains(const std::vector<T> &vec, const T &value)
{
    return std::find(vec.begin(), vec.end(), value) != vec.end();
}

/*******************************************************************************
Compute homography transformation between images using RANSAC.
    matches - set of matching points between images
    numMatches - number of matching points
    numIterations - number of iterations to run RANSAC
    inlierThreshold - maximum distance between points that are considered to be inliers
    hom - returned homography transformation (image1 -> image2)
    homInv - returned inverse homography transformation (image2 -> image1)
    image1Display - image used to display matches
    image2Display - image used to display matches
*******************************************************************************/
void MainWindow::RANSAC(CMatches *matches, int numMatches, int numIterations, double inlierThreshold,
                        double hom[3][3], double homInv[3][3], QImage &image1Display, QImage &image2Display)
{
    // Add your code here.

    std::vector<int> randomNumber;
    CMatches *subsetmatches = new CMatches[4];

    double bestHomography[3][3];
    int numOfInliers = 0.0;
    // Picking 4 random numbers

    for(int iter = 0; iter < numIterations; iter++) {
        randomNumber.clear();
        
        while(1) {
            int random =  (rand() % (numMatches -1));
            if(randomNumber.size() < 4) { 
                if(contains(randomNumber, random))
                    continue;
                else 
                    randomNumber.push_back(random);
            } else
                break;
        }
  

    // Computing Homography
    
    for(int random = 0; random < randomNumber.size(); random++) {
        //subsetmatches[random] = 0.0;
        subsetmatches[random] = matches[randomNumber[random]];
    }

    ComputeHomography(subsetmatches, 4, hom, true);
    int currInlier = ComputeInlierCount(hom, matches, numMatches, inlierThreshold);
    if( currInlier > numOfInliers) {
        for(int i=0; i< 3; i++)
            for(int j=0; j< 3; j++)
                bestHomography[i][j] = hom[i][j];
        numOfInliers = currInlier;
    }
    }

     numOfInliers = ComputeInlierCount(bestHomography, matches, numMatches, inlierThreshold);
	
    // Add your code here.
    double projectedX1, projectedY1;
    int totalInliers = 0;
    double distance = 0.0;

    CMatches *inliers = new CMatches[numOfInliers];

    int inliercount = 0;

    for(int num = 0; num < numMatches; num++) {
        Project(matches[num].m_X1, matches[num].m_Y1, projectedX1, projectedY1, bestHomography);
        double x1diff = projectedX1 - matches[num].m_X2;
        double y1diff = projectedY1 - matches[num].m_Y2;
        distance = sqrt((x1diff * x1diff) + (y1diff * y1diff));

        if(distance < inlierThreshold) {
            inliers[inliercount] = matches[num];
            inliercount++;
        }
            totalInliers++;
    }

    ComputeHomography(inliers, inliercount, hom, true);
    ComputeHomography(inliers, inliercount, homInv, false);

    // After you're done computing the inliers, display the corresponding matches.
    DrawMatches(inliers, inliercount, image1Display, image2Display);

}

/*******************************************************************************
Bilinearly interpolate image (helper function for Stitch)
    image - input image
    (x, y) - location to interpolate
    rgb - returned color values

    You can just copy code from previous assignment.
*******************************************************************************/

double* GetIntermediatePixelValueByInterpolation(double influenceOfPixel1, double influenceOfPixel2, 
                                              double* pixel1Weight, double* pixel2Weight) {
    
    double* imageRGB = new double[3];
    imageRGB[0] = influenceOfPixel1 * pixel2Weight[0] + influenceOfPixel2 * pixel1Weight[0];
    imageRGB[1] = influenceOfPixel1 * pixel2Weight[1] + influenceOfPixel2 * pixel1Weight[1];
    imageRGB[2] = influenceOfPixel1 * pixel2Weight[2] + influenceOfPixel2 * pixel1Weight[2];

    return imageRGB;
}

void GetPixelValueByInterpolation(double colPixel, double rowPixel, double** pixelWeight, double rgb[3]) {
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

   rgb[0] = XY[0];
   rgb[1] = XY[1];
   rgb[2] = XY[2];
 
}

void GetNeighbouringPixelPositionAndWeight(QImage *image, double colPixel, double rowPixel, Position* position, double** pixelWeight) {
    int colPixel1 = (int) ((colPixel));
    int rowPixel1 = (int) ((rowPixel));

    position[0] = Position(rowPixel1,colPixel1);
    GetPixelWeight(image, Position(rowPixel1,colPixel1), pixelWeight[0]);

    position[1] = Position(rowPixel1,colPixel1+1);
    GetPixelWeight(image, Position(rowPixel1,colPixel1+1), pixelWeight[1]);

    position[2] = Position(rowPixel1+1,colPixel1);
    GetPixelWeight(image, Position(rowPixel1+1,colPixel1), pixelWeight[2]);

    position[3] = Position(rowPixel1+1,colPixel1+1);
    GetPixelWeight(image, Position(rowPixel1+1,colPixel1+1), pixelWeight[3]);
}

bool MainWindow::BilinearInterpolation(QImage *image, double colPixel, double rowPixel, double rgb[3])
{
    // Add your code here.
     Position* position = new Position[4];
    double** pixelWeight = new double*[4](); 
    for(int i=0; i< 4; i++)
        pixelWeight[i] = new double[3];
       
    if( colPixel < 0 || rowPixel < 0 || colPixel > (image->width() -1) || rowPixel > (image->height() -1)) {
        rgb[0] = 0; rgb[1] =0; rgb[2]=0; 
       return true;
    }

    GetNeighbouringPixelPositionAndWeight(image, colPixel, rowPixel, position, pixelWeight);
    GetPixelValueByInterpolation(colPixel, rowPixel, pixelWeight, rgb);


    return true;
}


/*******************************************************************************
Stitch together two images using the homography transformation
    image1 - first input image
    image2 - second input image
    hom - homography transformation (image1 -> image2)
    homInv - inverse homography transformation (image2 -> image1)
    stitchedImage - returned stitched image
*******************************************************************************/



void MainWindow::Stitch(QImage image1, QImage image2, double hom[3][3], double homInv[3][3], QImage &stitchedImage)
{
    // Width and height of stitchedImage
    int ws = 0;
    int hs = 0;

    // Add your code to compute ws and hs here.

    int image2Width = image2.width();
    int image2Height = image2.height();

    int image1Width = image1.width();
    int image1Height = image1.height();

    // project four corners of image2 onto image1;
    double projectedCorner[4][2];
    Project(0,0,projectedCorner[0][0],projectedCorner[0][1], homInv); // x - column , y - row
    Project(image2Width-1,0,projectedCorner[1][0],projectedCorner[1][1], homInv);
    Project(0,image2Height-1,projectedCorner[2][0],projectedCorner[2][1], homInv);
    Project(image2Width-1,image2Height-1,projectedCorner[3][0],projectedCorner[3][1], homInv);

   // compare each corner of image1, with all corners of image 2.

    double min_height = (double)min(0,min((int)floor(projectedCorner[0][1]), min((int)floor(projectedCorner[1][1]), min((int)floor(projectedCorner[2][1]), (int)floor(projectedCorner[3][1])))));
    double max_height = (double)max(image1Height, max((int)ceil(projectedCorner[0][1]), max((int)ceil(projectedCorner[1][1]), max((int)ceil(projectedCorner[2][1]), (int)ceil(projectedCorner[3][1])))));
    double min_width  = (double)min(0,min((int)floor(projectedCorner[0][0]), min((int)floor(projectedCorner[1][0]), min((int)floor(projectedCorner[2][0]), (int)floor(projectedCorner[3][0])))));
    double max_width  = (double)max(image1Width,max((int)ceil(projectedCorner[0][0]), max((int)ceil(projectedCorner[1][0]), max((int)ceil(projectedCorner[2][0]), (int)ceil(projectedCorner[3][0])))));

   /* double min_height = (double)min((int)floor(projectedCorner[0][1]), min((int)floor(projectedCorner[1][1]), min((int)floor(projectedCorner[2][1]), (int)floor(projectedCorner[3][1]))));
    double max_height = (double)max((int)ceil(projectedCorner[0][1]), max((int)ceil(projectedCorner[1][1]), max((int)ceil(projectedCorner[2][1]), (int)ceil(projectedCorner[3][1]))));
    double min_width  = (double)min((int)floor(projectedCorner[0][0]), min((int)floor(projectedCorner[1][0]), min((int)floor(projectedCorner[2][0]), (int)floor(projectedCorner[3][0]))));
    double max_width  = (double)max((int)ceil(projectedCorner[0][0]), max((int)ceil(projectedCorner[1][0]), max((int)ceil(projectedCorner[2][0]), (int)ceil(projectedCorner[3][0]))));
*/


    ws = abs(max_width - min_width);
    hs = abs(max_height - min_height);


    stitchedImage = QImage(ws, hs, QImage::Format_RGB32);
    stitchedImage.fill(qRgb(0,0,0));

    // Add you code to warp image1 and image2 to stitchedImage here.
    double Xprojectedvalue;
    double Yprojectedvalue;
  

    // copying image1 onto stiched imge

    /*for(int rowPixel = 0; rowPixel < image1Height; rowPixel++) {
        for(int colPixel = 0; colPixel < image1Width; colPixel++) {
            QRgb pixel = image1.pixel(colPixel, rowPixel);
            stitchedImage.setPixel(colPixel+abs(min_width), rowPixel+abs(min_height), qRgb((int)floor(qRed(pixel)), (int)floor(qGreen(pixel)), (int)floor(qBlue(pixel))));
        }
    }*/

    //// projecting stitched image onto image2
    double rgb[3];
    for(int rowPixel = 0; rowPixel < hs; rowPixel++) {
        for(int colPixel = 0; colPixel < ws; colPixel++) {
            Project(colPixel, rowPixel, Xprojectedvalue, Yprojectedvalue, hom);
         
            if(Xprojectedvalue >= 0 && Xprojectedvalue < image2Width && Yprojectedvalue >=0 && Yprojectedvalue < image2Height) {
                BilinearInterpolation(&image2, Xprojectedvalue,Yprojectedvalue, rgb);
                if(rowPixel > min_height)
                    stitchedImage.setPixel(colPixel, rowPixel, qRgb((int)floor(rgb[0]), (int)floor(rgb[1]), (int)floor(rgb[2])));
                else
                    stitchedImage.setPixel(colPixel  , rowPixel  + abs(min_width), qRgb((int)floor(rgb[0]), (int)floor(rgb[1]), (int)floor(rgb[2])));  
            }
        }
      }



}

