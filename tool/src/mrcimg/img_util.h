#ifndef IMG_UTIL_H__
#define IMG_UTIL_H__

#include "dataf/keypoint.h"
#include "util/exception.h"
#include <vector>
#include <fstream>
// #include "highgui.h"
// #include "cxcore.h"
// #include "cv.h"
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/core_c.h>
#include <opencv2/opencv.hpp>

namespace util {

// static void DrawFeature(IplImage* img, const util::feature& feat);
// 
// inline void DrawX(IplImage* img, int x, int y)
// {
//     CvScalar color = CV_RGB(255, 255, 255);
//     int r = 3;
//     cvLine(img, cvPoint(x - r, y - r), cvPoint(x + r, y + r), color, 3, 8, 0);
//     cvLine(img, cvPoint(x + r, y - r), cvPoint(x - r, y + r), color, 3, 8, 0);
// }
// 
inline void DrawLine(cv::Mat& img, const util::_point& pt1, const util::_point& pt2)
{
    cv::Scalar color = CV_RGB(255, 255, 255);
    cv::line(img, cv::Point(pt1.x, pt1.y), cv::Point(pt2.x, pt2.y), color);
}
// 
// inline void DrawFeatures(IplImage* img, const std::vector<util::feature>& feats)
// {
// //     CvScalar color = CV_RGB(255, 255, 255);
// 
//     for(int i = 0; i < feats.size(); i++) {
//         DrawFeature(img, feats[i]);
//     }
// }
// 
// static void DrawFeature(IplImage* img, const util::feature& feat)
// {
//     CvScalar color = CV_RGB(255, 255, 255);
// 
//     int len, hlen, blen, start_x, start_y, end_x, end_y, h1_x, h1_y, h2_x, h2_y;
//     double scl, ori;
//     double scale = 5.0;
//     double hscale = 0.75;
//     CvPoint start, end, h1, h2;
// 
//     /* compute points for an arrow scaled and rotated by feat's scl and ori */
//     start_x = cvRound(feat.kpt.x);
//     start_y = cvRound(feat.kpt.y);
//     scl = feat.kpt.sigma;
//     ori = feat.kpt.orient;
//     len = cvRound(scl * scale);
//     hlen = cvRound(scl * hscale);
//     blen = len - hlen;
//     end_x = cvRound(len *  cos(ori)) + start_x;
//     end_y = cvRound(len * -sin(ori)) + start_y;
//     h1_x = cvRound(blen *  cos(ori + CV_PI / 18.0)) + start_x;
//     h1_y = cvRound(blen * -sin(ori + CV_PI / 18.0)) + start_y;
//     h2_x = cvRound(blen *  cos(ori - CV_PI / 18.0)) + start_x;
//     h2_y = cvRound(blen * -sin(ori - CV_PI / 18.0)) + start_y;
//     start = cvPoint(start_x, start_y);
//     end = cvPoint(end_x, end_y);
//     h1 = cvPoint(h1_x, h1_y);
//     h2 = cvPoint(h2_x, h2_y);
// 
//     cvLine(img, start, end, color, 1, 8, 0);
//     cvLine(img, end, h1, color, 1, 8, 0);
//     cvLine(img, end, h2, color, 1, 8, 0);
// }

inline void SaveImage(const cv::Mat& img, const char* filename) {
    EX_TRACE("Save Image %s\n", filename)
    cv::Mat copy = cv::Mat(img.size(), CV_8UC1);
	
	for(int y = 0; y < img.size().height; y++){
		float* src = (float*)(img.data+y*img.step);
		char* start = (char*)(copy.data+y*copy.step);
		for(int x = 0; x < img.size().width; x++){
			*start++ = (*src++)*255;
		}
	}
	
	cv::imwrite(filename, copy);
   // cvSaveImage(filename, copy);
}

// typedef struct _IplImage
// {
//     int  nSize;             /* sizeof(IplImage) */
//     int  ID;                /* version (=0)*/
//     int  nChannels;         /* Most of OpenCV functions support 1,2,3 or 4 channels */
//     int  alphaChannel;      /* Ignored by OpenCV */
//     int  depth;             /* Pixel depth in bits: IPL_DEPTH_8U, IPL_DEPTH_8S, IPL_DEPTH_16S,
//                                IPL_DEPTH_32S, IPL_DEPTH_32F and IPL_DEPTH_64F are supported.  */
//     char colorModel[4];     /* Ignored by OpenCV */
//     char channelSeq[4];     /* ditto */
//     int  dataOrder;         /* 0 - interleaved color channels, 1 - separate color channels.
//                                cvCreateImage can only create interleaved images */
//     int  origin;            /* 0 - top-left origin,
//                                1 - bottom-left origin (Windows bitmaps style).  */
//     int  align;             /* Alignment of image rows (4 or 8).
//                                OpenCV ignores it and uses widthStep instead.    */
//     int  width;             /* Image width in pixels.                           */
//     int  height;            /* Image height in pixels.                          */
//     struct _IplROI *roi;    /* Image ROI. If NULL, the whole image is selected. */
//     struct _IplImage *maskROI;      /* Must be NULL. */
//     void  *imageId;                 /* "           " */
//     struct _IplTileInfo *tileInfo;  /* "           " */
//     int  imageSize;         /* Image data size in bytes
//                                (==image->height*image->widthStep
//                                in case of interleaved data)*/
//     char *imageData;        /* Pointer to aligned image data.         */
//     int  widthStep;         /* Size of aligned image row in bytes.    */
//     int  BorderMode[4];     /* Ignored by OpenCV.                     */
//     int  BorderConst[4];    /* Ditto.                                 */
//     char *imageDataOrigin;  /* Pointer to very origin of image data
//                                (not necessarily aligned) -
//                                needed for correct deallocation */    std::cout<<"test:"<<img.channels<<std::endl;
// }
// IplImage;

inline void SeriesSaveToFile(const cv::Mat& img, const char* filename){
	std::ofstream ofs;
	ofs.open(filename, std::ofstream::binary);
    int width = img.size().width;
    int height = img.size().height;
    int depth = img.depth();
    int channels = img.channels();
    short em_size = img.elemSize();
    short type = img.type();

    ofs.write(reinterpret_cast<const char*>(&width), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&height), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&depth), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&channels), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&em_size), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&type), sizeof(int));
	ofs.write(reinterpret_cast<const char*>(img.data), sizeof(char) * em_size * width * height);
	ofs.close();
}

inline void SeriesReadFromFile(cv::Mat* img, const char* filename){
	std::ifstream ifs;
	ifs.open(filename, std::ofstream::binary);
	int channl(0);
    int width(0);
    int height(0);
    int depth(0);
    short em_size;
    short type;
    
    ifs.read(reinterpret_cast<char*>(&width), sizeof(int));
    ifs.read(reinterpret_cast<char*>(&height), sizeof(int));
    ifs.read(reinterpret_cast<char*>(&depth), sizeof(int));
    ifs.read(reinterpret_cast<char*>(&channl), sizeof(int));
    ifs.read(reinterpret_cast<char*>(&em_size), sizeof(int));
    ifs.read(reinterpret_cast<char*>(&type), sizeof(int));
    *img = cv::Mat(width, height, type);
    ifs.read((char *)img->data, em_size * width * height);
}

inline cv::Mat GetSubImage(cv::Mat& image, cv::Rect roi)
{
    cv::Mat result;
	
	if(roi.x < 0 || roi.x+roi.width >= image.size().width || roi.y < 0 || roi.y+roi.height >= image.size().height){
		return cv::Mat();
       // return NULL;
	}
	
	cv::Mat roilimage=image(roi);
	if(roilimage.size().width != roi.width || roilimage.size().height != roi.height){
		return cv::Mat();
        //return NULL;
	}
	result = cv::Mat(cv::Size(roi.width, roi.height), image.depth(), image.channels());
    result=roilimage.clone();
    
    return result;
}

inline void ScaleImage(cv::Mat& image, float scale)
{
	if(scale == 1){
		return;
	}
	
	//IplImage* cpy = cvCreateImage(cvSize(image->width*scale, image->height*scale), image->depth, image->nChannels);
	cv::Size size(image.size().width*scale, image.size().height*scale);
	cv::Mat cpy(size, image.depth(), image.channels());
	//cvResize(image, cpy, CV_INTER_CUBIC);
	cv::resize(image, cpy, size, 0, 0, cv::INTER_CUBIC);
    //cvReleaseImage(&image);
	image = cpy;
}

inline void Reversal(cv::Mat& img)
{
    double minimum, maximum;
    cv::minMaxLoc(img, &minimum, &maximum);
    double p = maximum + minimum;
    for(int y = 0; y < img.size().height; y++){
		float* ptr = (float*)(img.data+y*img.step);
		for(int x = 0; x < img.size().width; x++){
			*ptr =p-*ptr;
			ptr++;
		}
	}
}

// 	cvErode( img, img, NULL, 5);
//      cvDilate(img, img, NULL, 5);

inline void MedianSmooth(cv::Mat& img)
{
    //cvSmooth(img, img, CV_MEDIAN, 5);
    cv::medianBlur(img, img, 5);
}

#define HISTOGRAM_BIN			1024//	256
// inline void HistogramStretch(IplImage *img, int binsize = 1024)
// {
//     cv::Mat mat(img, 0);
//     double minimum, maximum;
//     cvMinMaxLoc(img, &minimum, &maximum);
//     minimum = minimum<0?minimum:0;
//     maximum = maximum>1?maximum:1;
//     float range[] = {minimum, maximum};
//     const float* ranges = {range};
//     cv::Mat hist;
//     const int channel = 0;
// 	
//     cv::calcHist(&mat, 1, &channel, cv::Mat(), hist, 1, &binsize, &ranges);	/*warning: if has some problem, look here */
//     cv::normalize(hist, hist, 1, 0, CV_L1);
//     float* vlu = (float*)hist.data;
// 
//     for(int i = 1; i < binsize; i++) {
//         vlu[i] += vlu[i-1];
//     }
// 
//     float delta = (maximum-minimum)/(binsize);
// 
//     float res = delta/2/(maximum-minimum);
//     for(int x = 0; x < img->height; x++) {
//         for(int y = 0; y < img->width; y++) {
//             int loc = (int)((CV_IMAGE_ELEM(img, float, x, y)-minimum)/delta);
//             CV_IMAGE_ELEM(img, float, x, y) = (maximum-minimum)*(vlu[loc]+delta)+minimum;
//         }
//     }
// }
#undef HISTOGRAM_BIN


inline void ConvertTo1(cv::Mat& img, bool cutoff = false)
{
    double minimum, maximum;
    cv::Scalar mean, sdv;
    cv::meanStdDev(img, mean, sdv);
    cv::minMaxLoc(img, &minimum, &maximum);
    
    double m,sd=0;

#define CUTOFF (3)

    maximum = (maximum > mean.val[0]+CUTOFF*sdv.val[0])? mean.val[0]+CUTOFF*sdv.val[0] : maximum;
    minimum = (minimum > mean.val[0]-CUTOFF*sdv.val[0])? minimum : mean.val[0]-CUTOFF*sdv.val[0];
    
	
	float span;
	(maximum-minimum)==0?span = 1:(span = maximum-minimum);

	for(int y = 0; y < img.size().height; y++){
		float* ptr = (float*)(img.ptr()+y*img.step);
		for(int x = 0; x < img.size().width; x++){
			if(cutoff){
				if(fabs(*ptr-mean.val[0]) <= CUTOFF*sdv.val[0]){
					*ptr = (*ptr-minimum)/span;
				}
				else {
					if(*ptr > mean.val[0]) {
						*ptr = 1;
					}
					else {
						*ptr = 0;
					}
				}
			}
			else{
				*ptr = (*ptr-minimum)/span;
			}
			ptr++;
		}
	}
	
#undef CUTOFF
}

// inline void CopyToUChar256(const IplImage* img, IplImage** cpy, int w, int b, bool cutoff = false)
// {
//     *cpy = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_8U, 1);
//     double minimum, maximum;
//     CvScalar mean, sdv;
//     cvAvgSdv(img, &mean, &sdv);
//     cvMinMaxLoc(img, &minimum, &maximum);
// 
// #define CUTOFF (3)//(4)//(5)
// 
//     maximum = (maximum > mean.val[0]+CUTOFF*sdv.val[0])? mean.val[0]+CUTOFF*sdv.val[0] : maximum;
//     minimum = (minimum > mean.val[0]-CUTOFF*sdv.val[0])? minimum : mean.val[0]-CUTOFF*sdv.val[0];
//     
//     for(int y = 0; y < img->height; y++){
// 		float* ptr = (float*)(img->imageData+y*img->widthStep);
// 		uchar* dst = (uchar*)((*cpy)->imageData+y*(*cpy)->widthStep);
// 		for(int x = 0; x < img->width; x++){
// 			if(cutoff){
// 				if(fabs(*ptr-mean.val[0]) <= CUTOFF*sdv.val[0]) {
// 					*dst = (uchar)((*ptr-minimum)*(w-b)/((maximum-minimum)==0?1:(maximum-minimum))+b);
// 					*dst = *dst > w ? w : *dst;
// 				}
// 				else {
// 					if(*ptr > mean.val[0]) {
// 						*dst = w;
// 					}
// 					else {
// 						*dst = b;
// 					}
// 				}
// 			}
// 			else{
// 				*dst = (uchar)((*ptr-minimum)*(w-b)/((maximum-minimum)==0?1:(maximum-minimum))+b);
// 				*dst = *dst > w ? w : *dst;
// 			}
// 			ptr++;
// 			dst++;
// 		}
// 	}
// 
// #undef CUTOFF
// }
// 
// inline void LaplaceSharpen(IplImage* img)
// {
//     CvMat* kernel;
//     kernel = cvCreateMat(3,3,CV_32F);
//     cvmSet(kernel, 0, 0, 0);
//     cvmSet(kernel, 0, 1, -1);
//     cvmSet(kernel, 0, 2, 0);
//     cvmSet(kernel, 1, 0, -1);
//     cvmSet(kernel, 1, 1, 5);
//     cvmSet(kernel, 1, 2, -1);
//     cvmSet(kernel, 2, 0, 0);
//     cvmSet(kernel, 2, 1, -1);
//     cvmSet(kernel, 2, 2, 0);
// 
//     cvFilter2D(img, img, kernel);
// }

}

#endif
