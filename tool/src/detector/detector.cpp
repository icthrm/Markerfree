#include "detector.h"
#include "mrcimg/img_util.h"
#include "ransac/xform.h"
#include <util/exception.h>
#include <iostream>
#include <stack>
#include <algorithm>
#include <cstdio>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>
#include <malloc.h>

Detector Detector::detector;

//CV_IMAGE_ELEM
static inline float pixval32f( cv::Mat& img, int r, int c )
{
    return ( (float*)(img.data + img.step*r))[c];
}

static void DrawFiducialMarkerPositions(cv::Mat& img, float diameter, const std::vector<util::point2d>& fids)
{
	cv::Scalar color = CV_RGB(255, 0, 0);
	int r = int(diameter+.5), t = int(.1*r+.5);
	if(t <= 0){
		t = 1;
	}
	for(int k = 0; k < fids.size(); k++){
		int x = int(fids[k].x+.5), y = int(fids[k].y+.5);
		cv::circle(img, cv::Point(x, y), r, color, t);
	}	
}


Detector::Detector():fid_diameter(DEFAULT_FIDD)
{

}

Detector::~Detector()
{
	if(!fid_avgtmplt.empty()){
        fid_avgtmplt.release();
// 		cvReleaseImage(&fid_avgtmplt);
	}
	if(!fid_tmplt.empty()){
        fid_tmplt.release();
// 		cvReleaseImage(&fid_tmplt);
	}
}

void Detector::GenerateDiameterSpace(std::vector<float>& diavec, float dia_begin, float dia_end, int bin_size)
{
	diavec.clear();
	float step = (dia_end-dia_begin)/(bin_size-1);
	for(int i = 0; i < bin_size; i++){
		diavec.push_back(dia_begin+i*step);
	}
}


void Detector::CalculateCorralation(const cv::Mat& src, const cv::Mat& tmplt, cv::Mat* dst)
{
    *dst=cv::Mat::zeros(cv::Size(src.size().width-tmplt.size().width+1, src.size().height-tmplt.size().height+1), CV_32FC1);
    //*dst = cvCreateImage(cvSize(src.size().width-tmplt.size().width+1, src.size().height-tmplt.size().height+1), IPL_DEPTH_32F, 1);
    
    cv::matchTemplate(src, tmplt, *dst, cv::TM_CCORR_NORMED);//CV_TM_CCOEFF_NORMED);//
}

void Detector::FindMaxPeak(cv::Mat& corr, int seed_x, int seed_y, int idiameter, int* peak_x, int* peak_y)
{
	*peak_x = seed_x, *peak_y = seed_y;
    float max_value = corr.at<float>(seed_y, seed_x);
	//float max_value = CV_IMAGE_ELEM(corr, float, seed_y, seed_x);
	for(int x = seed_x; x < seed_x+idiameter && x < corr.size().width; x++){
		for(int y = seed_y; y < seed_y+idiameter && y < corr.size().height; y++){
			if(corr.at<float>(y, x) > max_value){
				//max_value = CV_IMAGE_ELEM(corr, float, y, x);
                max_value = corr.at<float>(y, x);
				*peak_x = x;
				*peak_y = y;
			}
		}
	}
}


float Detector::GetAverageOfMarkerPatch(cv::Mat& img, int seed_x, int seed_y, int diameter)
{
	int iradius = int(diameter*.5+.5);
	int iradius2 = iradius*iradius*.81;
	if(iradius2 == 0){
		iradius2 = iradius;
	}
	
	int centre_x = seed_x+iradius, centre_y = seed_y+iradius;
	float value = 0;
	int count = 0;
	for(int x = seed_x; x < seed_x+diameter && x < img.size().width; x++){
		for(int y = seed_y; y < seed_y+diameter && y < img.size().height; y++){
			if((x-centre_x)*(x-centre_x)+(y-centre_y)*(y-centre_y) < iradius2){
// 				value += CV_IMAGE_ELEM(img, float, y, x);
                value += img.at<float>(y, x);
				count++;
			}
		}
	}
	
	return value/count;
}


void Detector::FindFiducialMarkerPositions(cv::Mat& img, const cv::Mat& tmplt, float diameter, std::vector< util::point2d >& fids, float* positive_ratio, bool limit_contrast, bool forsubtomo)
{
#define DIA_MIN_THRE		0.36
#define DIA_MAX_THRE		1.414
	
	fid_diameter = diameter;
	
	cv::Mat corr;
	
	CalculateCorralation(img, tmplt, &corr);

// 	if(forsubtomo)
// 	{
// 		std::cout<<corr<<std::endl;
// 	}

// 	std::cout<<corr<<std::endl;
    
	util::ConvertTo1(corr, false);//true);
// 	cvNormalize(corr, corr, 1, 0, CV_MINMAX); 
// 	util::ConvertTo1(corr, true);
// 	util::SaveImage(corr, "corr.pgm");
	
	cv::Scalar corr_avg, corr_std;
	//cvAvgSdv(corr, &corr_avg, &corr_std);
    cv::meanStdDev(corr, corr_avg, corr_std);
	
	cv::Scalar pixel_avg, pixel_std;
	//cvAvgSdv(img, &pixel_avg, &pixel_std);
    cv::meanStdDev(img, pixel_avg, pixel_std);

	if(forsubtomo)
	{
		std::cout<<corr_avg<<std::endl;
	}
    

#define CUTOFF_THRE       0//2.5 
	
	if(sampling < 2500){
		corr_threshold = 0;//corr_avg.val[0]+CUTOFF_THRE*corr_std.val[0];
		pixel_threshold = 0;//pixel_avg.val[0]+CUTOFF_THRE*pixel_std.val[0];
	}
	else{
		corr_threshold = corr_avg.val[0]+CUTOFF_THRE*corr_std.val[0];
		pixel_threshold = pixel_avg.val[0]+CUTOFF_THRE*pixel_std.val[0];
	}

	if(forsubtomo)
	{
		std::cout<<"corr_threshold: "<<corr_threshold<<"    "<<pixel_threshold<<std::endl;
	}

#undef 	CUTOFF_THRE
// 	corr_real_threshold = 0;
// 	pixel_real_threshold = 0;
	
	std::vector<ImgRef> refvec;
	int idia = int(diameter), iradius = int(diameter*.5); 
	int sqrt2idia = int(idia*1.414);
	
	for(int x = 0; x < corr.size().width; x += sqrt2idia){
		for(int y = 0; y < corr.size().height; y += sqrt2idia){
			int peak_x, peak_y;
			FindMaxPeak(corr, x, y, sqrt2idia, &peak_x, &peak_y);
//             std::cout<<"x:"<<x<<"y:"<<y<<std::endl;
//             std::cout<<"peak_x:"<<peak_x<<"peak_y:"<<peak_y<<std::endl;
			//if(corr.at<float>(peak_y, peak_x) > corr_threshold && CV_IMAGE_ELEM(img, float, peak_y+iradius, peak_x+iradius) > pixel_threshold){
            if(corr.at<float>(peak_y, peak_x) > corr_threshold && img.at<float>(peak_y+iradius, peak_x+iradius) > pixel_threshold){
				refvec.push_back(ImgRef(corr, peak_x, peak_y));
			}
		}
	}
// 	std::cout<<refvec.size()<<std::endl;
	std::vector<util::point2d> raw_fids;
	std::vector<std::pair<float, float> > scores;
	GetRawPositions(img, refvec, diameter, raw_fids, scores);

// 	if(test){
// 		for(int i = 0; i < scores.size(); i++){
// 			std::cout<<scores[i].first<<" "<<scores[i].second<<std::endl;
// 		}
// 		
// 		for(int i = 0; i < scores.size(); i++){
// 			std::cout<<scores[i].first*scores[i].second<<std::endl;
// 		}
// 		
// 		std::cout<<raw_fids.size()<<std::endl;
// 	}
	
	std::vector<util::point2d> new_fids;
	
 	RefineFiducialMarkersByGaussianDistribution(raw_fids, scores, new_fids);
	
	std::cout<<"Stage 1:"<<new_fids.size()<<std::endl;
	
// 	DrawFiducialMarkerPositions(img, new_fids);
// 	util::SaveImage(img, "img.pgm");
// 	
 	raw_fids.clear();
 	refvec.clear();
// 	
    cv::Mat img_cpy = img.clone();
 	//IplImage* img_cpy = cvCloneImage(img);
// 	
 	std::vector<util::point2d> candidates;
	for(int i = 0; i < new_fids.size(); i++){
		util::RECT region;
		RegionGrow(img_cpy, diameter, new_fids[i], region);
		float w = region.right-region.left;
		float h = region.bottom-region.top;
		if((w < DIA_MAX_THRE*diameter && h < DIA_MAX_THRE*diameter) && (!limit_contrast || (w >= DIA_MIN_THRE*diameter && h >= DIA_MIN_THRE*diameter))){  // && w >= DIA_MIN_THRE*diameter && h >= DIA_MIN_THRE*diameter){
			candidates.push_back(new_fids[i]);
		}
	}
	
	for(int i = 0; i < candidates.size(); i++){
		util::point2d fid;
		if(GetCenter(img, candidates[i], diameter, fid)){
			fids.push_back(fid);
		}
	}
// 	DrawFiducialMarkerPositions(img, fids);
// 	util::SaveImage(img, "img.pgm");
	if(forsubtomo){
		int magnif = 1;
		
		if(diameter <= 16){
			magnif = 16;
		}
		else if(diameter <= 32){
			magnif = 8;
		}
		else if(diameter <= 128){
			magnif = 4;
		}
		else if(diameter <= 150){
			magnif = 4;
		}
		else{
			magnif = 2;
		}
// 		std::cout<<tmplt<<std::endl;
		RefinePositionsForSubTomo(img, tmplt, fids, diameter, magnif);		//WARNING
	}
	std::cout<<"Stage 2:"<<fids.size()<<std::endl;
	*positive_ratio = (float)fids.size()/new_fids.size();
	
	//cvReleaseImage(&img_cpy);

#undef DIA_MIN_THRE	
#undef 	DIA_THRE
}


void Detector::RefineFiducialMarkersByGaussianDistribution(const std::vector< util::point2d >& raw_fids, 
														   const std::vector< std::pair< float, float > >& scores, std::vector< util::point2d >& new_fids)
{
	new_fids.clear();
	std::vector<float> sscores;
	for(int i = 0; i < scores.size(); i++){
		sscores.push_back(scores[i].first*scores[i].second);
	}
	float avg = 0, stdev = 0;
	for(int i = 0; i < sscores.size(); i++){
		avg += sscores[i];
		stdev += sscores[i]*sscores[i];
	}
	avg /= sscores.size();
	stdev = sqrt(stdev/sscores.size()-avg*avg);

// #define CUTOFF_THRE  2.5//2
// 	float thre = avg+CUTOFF_THRE*stdev;
// #undef 	CUTOFF_THRE
	float thre = avg+(sampling>2500?2.5:2)*stdev;
	
	for(int i = 0; i < raw_fids.size(); i++){
		if(sscores[i] > thre){
			new_fids.push_back(raw_fids[i]);
		}
	}
}


void Detector::GetRawPositions(cv::Mat& img, std::vector<ImgRef>& refvec, float diameter, 
							   std::vector<util::point2d>& raw_fids, std::vector<std::pair<float, float> >& scores)
{
	std::sort(refvec.begin(), refvec.end(), std::greater<ImgRef>());
	
	int dim = 2;
	int nk = 5;
	float radius = diameter/2;
	int nPts = refvec.size();
    ANNpoint qryPt = annAllocPt(dim);
    ANNpointArray dataPts = annAllocPts(nPts, dim);
    ANNidxArray nNIdx = new ANNidx[nk]; // near neighbor indices
    ANNdistArray dists = new ANNdist[nk]; // near neighbor distances
    ANNkd_tree* kdTree; // search structure

    for(int i = 0; i < refvec.size(); i++){
		dataPts[i][0] = refvec[i].x;
		dataPts[i][1] = refvec[i].y;
	}
	
	kdTree = new ANNkd_tree(dataPts, nPts, dim);
	float dists_thre = diameter;//*1.414;
	
	for(int i = 0; i < refvec.size(); i++){
		if(refvec[i].Value() < 0){
			continue;
		}
		
		qryPt[0] = refvec[i].x;
		qryPt[1] = refvec[i].y;
		kdTree->annkSearch(qryPt, nk, nNIdx, dists, 0.5);
		for(int j = 1; j < nk; j++){
			if(dists[j] < dists_thre){
				refvec[nNIdx[j]].Value() = -1;
 			}
		}
		
		raw_fids.push_back(util::point2d(refvec[i].x+radius, refvec[i].y+radius));
		std::pair<float, float> score;
		score.first = refvec[i].Value();
		
// 		int peak_x, peak_y; float half_radius = radius*.5;
// 		FindMaxPeak(img, refvec[i].x+half_radius, refvec[i].y+half_radius, half_radius, &peak_x, &peak_y);
// 		score.second = CV_IMAGE_ELEM(img, float, peak_y, peak_x);
		score.second = GetAverageOfMarkerPatch(img, refvec[i].x, refvec[i].y, int(diameter+.5));
		scores.push_back(score);
	}
	
	delete [] nNIdx;
    delete [] dists;
    delete kdTree;
    annDeallocPts(dataPts);
    annDeallocPt(qryPt);
}

bool Detector::RegionGrow(cv::Mat& img, float diameter, const util::point2d& seed, util::RECT& region)
{
	int idia = int(diameter*1.2533+.5);		//1.2533 = sqrt(pi/2)
	
	cv::Rect roi;
	roi.x = int(seed.x - idia*.5+.5);
	roi.y = int(seed.y - idia*.5+.5);
	roi.width = idia;
	roi.height = idia;
	cv::Mat sub = util::GetSubImage(img, roi);
    
	if(sub.empty()){
		return false;
	}
	
	cv::Scalar pixel_avg, pixel_std;
	cv::meanStdDev(sub, pixel_avg, pixel_std);
	
// 	std::cout<<pixel_avg.val[0]<<" "<<pixel_std.val[0]<<" "<<CV_IMAGE_ELEM(img, float, int(seed.y), int(seed.x))<<std::endl;
	
	
	int seed_y = seed.y;
	int seed_x = seed.x;
// 	float max_value = ( (float*)(img->imageData + img->widthStep*seed_y))[seed_x];//seed_y, seed_x);

    double d = pixel_avg.val[0]+.5*pixel_std.val[0];//max_value*RG_THRE;
    //std::cout<<"d:"<<d<<std::endl;
    std::stack<util::point2d> seedd;
    seedd.push(seed);

    region.left = seed.x;
    region.right = seed.x;
    region.top = seed.y;
    region.bottom = seed.y;

    int width = img.size().width;
    int height = img.size().height;
    while(!seedd.empty()){
		util::point2d point = seedd.top();
        seedd.pop();
		int intx = int(point.x);
		int inty = int(point.y);
		
		if(!(intx > 0 && intx < width-1 && inty > 0 && inty < height-1)){
			break;
		}
		
		((float*)(img.data + inty*img.step))[intx] = 0;
		float value = ((float*)(img.data + inty*img.step))[intx-1];		//(x-1, y)
		
		if(value > d){
			if(point.x-1 < region.left){
				region.left = point.x-1;
			}
			util::point2d tmp(point.x-1, point.y);
			seedd.push(tmp);
		}

		value = ((float*)(img.data + inty*img.step))[intx+1];			//(x+1, y)
		if(value > d){
			if(point.x+1 > region.right){
				region.right = point.x+1;
			}
			util::point2d tmp(point.x+1, point.y);
			seedd.push(tmp);
		}

		value = ((float*)(img.data + (inty-1)*img.step))[intx];		//(x, y-1)
		if(value > d){
			if(point.y-1 < region.top){
				region.top = point.y-1;
			}
			util::point2d tmp(point.x, point.y-1);
			seedd.push(tmp);
		}

		value = ((float*)(img.data + (inty+1)*img.step))[intx];		//(x, y+1)
		if(value > d){
			if(point.y+1 > region.bottom){
				region.bottom = point.y+1;
			}
			util::point2d tmp(point.x, point.y+1);
			seedd.push(tmp);
		}
    }
    
    if(region.left < 0 || region.right >= width || region.top < 0 || region.bottom >= height){
        return false;
	}

    return true;

#undef BG_FID_THRE
#undef RG_THRE
}

bool Detector::GetCenter(cv::Mat& img, const util::point2d& seed, float diameter, util::point2d& centre)
{
	int idia = int(diameter+.5);		//1.2533 = sqrt(pi/2)
	
	cv::Rect roi;
	roi.x = int(seed.x - idia*.5+.5);
	roi.y = int(seed.y - idia*.5+.5);
	roi.width = idia;
	roi.height = idia;
	cv::Mat sub = util::GetSubImage(img, roi);
	
	if(sub.empty()){
		return false;
	}
	
	std::vector<float> values;
	for(int i = 0; i < sub.size().width; i++){
		for(int j = 0; j < sub.size().height; j++){
			float value = ((float*)(sub.data + j*sub.step))[i];
			values.push_back(value);
		}
	}
	std::sort(values.begin(), values.end(), std::greater<float>());
    float d = values[values.size()*.785*.81];		//0.785 = pi/4
	
	float sx = 0, sy = 0, sum = 0;
	
	for(int i = 0; i < sub.size().width; i++){
		for(int j = 0; j < sub.size().height; j++){
			float value = ((float*)(sub.data + j*sub.step))[i];
            if(value > d){
                sx += i*value;
                sy +=j*value;
                sum += value;
            }
		}
	}
	
	
	if(sum != 0){
		sx /= sum;
		sy /= sum;
		centre.x = roi.x+sx;
		centre.y = roi.y+sy;
		return true;
	}
	return false;	
}


void Detector::CreateTemplate(cv::Mat* tmplt, float diameter)
{
#define TMP_ROUND(x)		int(x+0.5)
#define TEMPLATE_RATIO		16
	float radius = diameter*.5;
	int dia = TMP_ROUND(radius*TEMPLATE_RATIO);
    *tmplt = cv::Mat::zeros(cv::Size(TMP_ROUND(diameter+1), TMP_ROUND(diameter+1)), CV_32FC1);
    cv::Mat tmp(cv::Size(2*dia+1, 2*dia+1), CV_32FC1, cv::Scalar::all(-.5f));
    cv::circle(tmp, cv::Point(dia, dia), dia, CV_RGB(.5f, .5f, .5f), -1, 8);//CV_RGB(1, 1, 1), -1, 8);//
	cv::resize(tmp, *tmplt, cv::Size(TMP_ROUND(diameter+1), TMP_ROUND(diameter+1)), 0, 0, cv::INTER_CUBIC);
    
}



void Detector::GaussianSmoothBasedOnSampling(cv::Mat& img, int sampling)
{
	if(sampling > 2500){
		
	}
	else if(sampling > 1250){
		//cvSmooth(img, img, CV_GAUSSIAN, 3);
        cv::GaussianBlur(img, img, cv::Size(3, 3), 0, 0);
	}
	else{
		cv::GaussianBlur(img, img, cv::Size(5, 5), 0, 0);
        //cvSmooth(img, img, CV_GAUSSIAN, 5);
	}
}

void Detector::CreateAverageTemplate(cv::Mat& img, const std::vector< util::point2d >& fids, float diameter, cv::Mat* tmplt, bool zoom_out, int magnif)
{
	int TEMPLATE_BOUND = 4;
	float radius = diameter*.5;
	int dia = TMP_ROUND(radius*magnif);
// 	*tmplt = cvCreateImage(cvSize(TMP_ROUND(diameter+11), TMP_ROUND(diameter+11)), IPL_DEPTH_32F, 1);
// 	std::cout << img << std::endl;
    cv::Mat tmp(cv::Size(2*dia+1, 2*dia+1), CV_32FC1, cv::Scalar::all(-0.0f));
	//IplImage* tmp = cvCreateImage(cvSize(2*dia+1, 2*dia+1), IPL_DEPTH_32F, 1);
	//cvFillImage(tmp, 0.0f);

    for(int i = 0; i < fids.size(); i++){
		cv::Rect rect;
		rect.x = TMP_ROUND(fids[i].x-radius-TEMPLATE_BOUND);
		rect.y = TMP_ROUND(fids[i].y-radius-TEMPLATE_BOUND);
		rect.width = TMP_ROUND(diameter+TEMPLATE_BOUND*2);
		rect.height = TMP_ROUND(diameter+TEMPLATE_BOUND*2);
		
		*tmplt = util::GetSubImage(img, rect);
		if(tmplt->empty()){
			continue;
		}
        
		cv::Mat tmp_patch=cv::Mat::zeros(cv::Size(rect.width*magnif, rect.height*magnif), CV_32FC1);

// 		IplImage* tmp_patch = cvCreateImage(cvSize(rect.width*magnif, rect.height*magnif), IPL_DEPTH_32F, 1);
        //cvResize(*tmplt, tmp_patch, CV_INTER_CUBIC);
 		cv::resize(*tmplt, tmp_patch, tmp_patch.size(), 0, 0, cv::INTER_CUBIC);

		int fid_xr = int(magnif*fids[i].x)-rect.x*magnif;
		int fid_yr = int(magnif*fids[i].y)-rect.y*magnif;
		rect.x = fid_xr-dia;
		rect.y = fid_yr-dia;
		rect.width = 2*dia+1;
		rect.height = 2*dia+1;
		*tmplt = util::GetSubImage(tmp_patch, rect);
		cv::add(*tmplt, tmp, tmp);
    }

    
//     util::Reversal(tmp);
// 	util::ConvertTo1(tmp, false);
// 	util::SaveImage(tmp, "tmp.pgm");
    tmp.convertTo(tmp, -1, 1.f/(int)fids.size());
	if(!zoom_out){
		*tmplt = cv::Mat::zeros(cv::Size(TMP_ROUND(diameter+1), TMP_ROUND(diameter+1)), CV_32FC1);
		cv::resize(tmp, *tmplt, cv::Size(TMP_ROUND(diameter+1), TMP_ROUND(diameter+1)), 0, 0, cv::INTER_CUBIC);
	}
	else{
		*tmplt = tmp;
	}
}


void Detector::RefinePositionsForSubTomo(cv::Mat& img, const cv::Mat& tmplt, std::vector< util::point2d >& fids, float diameter, int magnif)
{
	float radius = diameter*.5;
	int dia = TMP_ROUND(radius*magnif);
	
	int TEMPLATE_BOUND = diameter*.1;
	if(TEMPLATE_BOUND<=3){
		TEMPLATE_BOUND = 3;
	}
// 	std::cout << img << std::endl;
	cv::Mat cpytmplt = cv::Mat(cv::Size(2*dia+1, 2*dia+1), tmplt.depth(), tmplt.channels());
	cv::resize(tmplt, cpytmplt, cv::Size(2*dia+1, 2*dia+1), 0, 0, cv::INTER_CUBIC);
	while(true){
		cv::Mat avgtmplt;
		CreateAverageTemplate(img, fids, diameter, &avgtmplt, true, magnif);
// 		std::cout << cpytmplt << std::endl;
		cv::add(cpytmplt, avgtmplt, avgtmplt);
		std::vector< util::point2d > finfids;
		float res = 0;
// 		std::cout << avgtmplt.size() << std::endl;
		for(int i = 0; i < fids.size(); i++){
// 			std::cout << fids[i].x<<"  "<< fids[i].y<< std::endl;
			cv::Rect rect;
			rect.x = TMP_ROUND(fids[i].x-radius-TEMPLATE_BOUND);
			rect.y = TMP_ROUND(fids[i].y-radius-TEMPLATE_BOUND);
			rect.width = TMP_ROUND(diameter+TEMPLATE_BOUND*2);
			rect.height = TMP_ROUND(diameter+TEMPLATE_BOUND*2);
// 			std::cout << rect.x<<"  "<< rect.y<<"  "<< rect.width<<"  "<< rect.height<<"  " << std::endl;
			cv::Mat fidtmplt = util::GetSubImage(img, rect);
			if(fidtmplt.empty()){
				continue;
			}
// 			std::cout << fidtmplt << std::endl;
			cv::Mat tmp_patch = cv::Mat::zeros(cv::Size(rect.width*magnif, rect.height*magnif), CV_32FC1);
			cv::resize(fidtmplt, tmp_patch, cv::Size(rect.width*magnif, rect.height*magnif), 0, 0, cv::INTER_CUBIC);
// 			std::cout << tmp_patch << std::endl;
			cv::Mat corr;
			CalculateCorralation(tmp_patch, avgtmplt, &corr);
// 			std::cout << tmp_patch << std::endl;
			int peak_x, peak_y;
			FindMaxPeak(corr, 0, 0, corr.size().width, &peak_x, &peak_y);
			peak_x += dia; peak_y += dia;

// 			std::cout <<"peak_x:"<< peak_x<<"peak_y:"<<peak_y << std::endl;
// 			std::cout<<(float)peak_x/magnif<<std::endl;
			util::point2d pt;
			pt.x = rect.x+(float)peak_x/magnif;
			pt.y = rect.y+(float)peak_y/magnif;
			
			res += sqrt((fids[i].x-pt.x)*(fids[i].x-pt.x)+ (fids[i].y-pt.y)*(fids[i].y-pt.y));
// 			std::cout << "***********************"<<res<<std::endl;
			finfids.push_back(pt);

// 			cvReleaseImage(&corr);
// 			cvReleaseImage(&tmp_patch);
		}
		
// 		cvReleaseImage(&avgtmplt);
		fids = finfids;
// 		std::cout << "***********************"<<res<<std::endl;
		res /= fids.size();
// 		std::cout<<res<<std::endl;
		
		if(res < 0.173){
			break;
		}
	}
	
// 	cvReleaseImage(&cpytmplt);
	
#undef TMP_ROUND
#undef TEMPLATE_RATIO	
}


float Detector::InitEstimateDiameter(const cv::Mat& img)
{
    
	cv::Size size(img.size().width, img.size().height);
    cv::Mat cpy = cv::Mat::zeros(size, CV_32FC1);
    //cv::Mat cpy = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_32F, 1);
    cv::resize(img, cpy, size, 0, 0, cv::INTER_CUBIC);
	//cvResize(img, cpy, CV_INTER_CUBIC);
    cv::GaussianBlur(cpy, cpy, cv::Size(9, 9), 0, 0);
	//cvSmooth(cpy, cpy, CV_GAUSSIAN, 9);
	
	std::vector<float> whv;
	
	int count = 0;
    //return 5;
	
	while(count < 6){
		int max_i, max_j;
		float max_value = -9999;
		for(int i = 0; i < cpy.size().width; i++){
			for(int j = 0; j < cpy.size().height; j++){
				float cpy_ij = ((float*)(cpy.ptr() + j*cpy.step))[i];
				if(cpy_ij > max_value){
					max_value = cpy_ij;
					max_i = i; max_j = j;
				}
			}
		}
		
		util::RECT region;
		BlindRegionGrow(cpy, max_i, max_j, region);
		
		float width = region.right-region.left;
		float height = region.bottom-region.top;
		
		if(width < 0.7*height || 0.7*width > height){
			continue;
		}
		
		if(region.left < img.size().width*.025 || region.right > img.size().width*.975 || region.top < img.size().width*.025 || region.bottom > img.size().width*.975){
			continue;
		}
		
		
		whv.push_back(width*.5+height*.5);
		
		std::cout<<region.right-region.left<<"\t"<<region.bottom-region.top<<std::endl;
		
		count++;
	}
	
	//cvReleaseImage(&cpy);
	
	std::sort(whv.begin(), whv.end());
	
 	return (whv[2]+whv[3])*.5;
}

void Detector::BlindRegionGrow(cv::Mat& img, int seed_x, int seed_y, util::RECT& region)
{
	float max_value = ( (float*)(img.ptr() + img.step*seed_y))[seed_x];//seed_y, seed_x);

    double d = (max_value*.65+.5*.35);//max_value*RG_THRE;
    std::stack<util::point2d> seedd;
	util::point2d seed;
	seed.x = seed_x; seed.y = seed_y;
	
    seedd.push(seed);


    region.left = seed.x;
    region.right = seed.x;
    region.top = seed.y;
    region.bottom = seed.y;

    int width = img.size().width;
    int height = img.size().height;
    while(!seedd.empty()){
		util::point2d point = seedd.top();
        seedd.pop();
		int intx = int(point.x);
		int inty = int(point.y);
		
		if(!(intx >= 0 && intx <= width-1 && inty >= 0 && inty <= height-1)){
			continue;
		}
		
		((float*)(img.ptr() + inty*img.step))[intx] = 0;
		
		float value = 0;
		if(intx-1 >= 0){
			value= ((float*)(img.ptr() + inty*img.step))[intx-1];		//(x-1, y)
		}
		
		if(value > d){
			if(point.x-1 < region.left){
				region.left = point.x-1;
			}
			util::point2d tmp(point.x-1, point.y);
			seedd.push(tmp);
		}

		if(intx+1 <= width-1){
			value = ((float*)(img.ptr() + inty*img.step))[intx+1];			//(x+1, y)
		}
		
		if(value > d){
			if(point.x+1 > region.right){
				region.right = point.x+1;
			}
			util::point2d tmp(point.x+1, point.y);
			seedd.push(tmp);
		}

		if(inty-1 >= 0){
			value = ((float*)(img.ptr() + (inty-1)*img.step))[intx];		//(x, y-1)
		}
		
		if(value > d){
			if(point.y-1 < region.top){
				region.top = point.y-1;
			}
			util::point2d tmp(point.x, point.y-1);
			seedd.push(tmp);
		}

		if(inty+1 <= height-1){
			value = ((float*)(img.ptr() + (inty+1)*img.step))[intx];		//(x, y+1)
		}
		
		if(value > d){
			if(point.y+1 > region.bottom){
				region.bottom = point.y+1;
			}
			util::point2d tmp(point.x, point.y+1);
			seedd.push(tmp);
		}
    }
}

bool Detector::DetermineDiameter(const cv::Mat& img, float diameter, float* refined_diameter, cv::Mat* avgtmplt)
{
#define OUTPUT_THRESHOLD	0.5
	
	EX_TRACE("Input diameter value = %.2f\n", diameter)
	
	float dbdia_begin;
	float dbdia_end;
	int bin_size;
	
	if(diameter < 0){
		EX_TRACE("Program will estimate the initial value of diameter for the user!\n")
		diameter = InitEstimateDiameter(img);
		diameter *= .9;
		EX_TRACE("Initial estimation of the diameter = %.2f\n", diameter)
		
		if(diameter < 30){
			dbdia_begin = diameter*1;
			dbdia_end = diameter*2;
			bin_size = 11;
		}
		else{
			dbdia_begin = diameter*1;
			dbdia_end = diameter*1.8;
			bin_size = int(diameter*.4)+1;
		}
	}
	else if(diameter < 20){
		dbdia_begin = diameter*.5;
		dbdia_end = diameter*1.5;
		bin_size = 11;
	}
	else if(diameter >= 20 && diameter < 30){					//diameter >= 30
		dbdia_begin = diameter-10;
		dbdia_end = diameter+10;
		bin_size = 11;
	}
	else if(diameter >= 30 && diameter < 40){					//diameter >= 30
		dbdia_begin = diameter-15;
		dbdia_end = diameter+15;
		bin_size = 16;
	}
	else if(diameter >= 40 && diameter < 60){
		dbdia_begin = diameter-18;
		dbdia_end = diameter+18;
		bin_size = 19;
	}
	else if(diameter >= 60 && diameter < 80){
		dbdia_begin = diameter-26;
		dbdia_end = diameter+26;
		bin_size = 27;
	}
	else if(diameter >= 80){
		dbdia_begin = diameter-32;
		dbdia_end = diameter+32;
		bin_size = 33;
	}
	
	EX_TIME_BEGIN("Refining the fiducial marker diameter...")

	while(true){
		std::vector<float> diavec;
		std::vector<float> ccvec;
		std::vector<cv::Mat> tmpltvec;
		GenerateDiameterSpace(diavec, dbdia_begin, dbdia_end, bin_size);
//         for (auto it = diavec.begin(); it != diavec.end(); it++)
//             std::cout << *it << " ";
        
 		bin_size = 11;
// 		diavec[0] = 8;
		for(int i = 0; i < diavec.size(); i++){
			std::vector<util::point2d> fids;
			cv::Mat manutmplt;
			CreateTemplate(&manutmplt, diavec[i]);
			
            cv::Size size(img.size().width, img.size().height);
            cv::Mat cpy=cv::Mat::zeros(size, CV_32FC1);
			//IplImage* cpy = cvCreateImage(cvSize(img.size().width, img.size().height), IPL_DEPTH_32F, 1);
		// 	cvSmooth(img, cpy, CV_GAUSSIAN, 5);
			cv::resize(img, cpy, size, 0, 0, cv::INTER_CUBIC);
			
			int csampling = int(cpy.size().width*cpy.size().height/(diavec[i]*diavec[i])*.5);		//in case of electorn noise in large scale
            

			GaussianSmoothBasedOnSampling(cpy, csampling);
			sampling = csampling;
			float positive_value;
			FindFiducialMarkerPositions(cpy, manutmplt, diavec[i], fids, &positive_value, true, false);
            
			
			cv::Mat tmplt, ccv;
			
			int magnif;
			if(diavec[i] <= 16){
				magnif = 16;
			}
			else if(diavec[i] <= 32){
				magnif = 8;
			}
			else if(diavec[i] <= 128){
				magnif = 4;
			}
			else if(diavec[i] <= 150){
				magnif = 4;
			}
			else{
				magnif = 2;
			}
			
			CreateAverageTemplate(cpy, fids, diavec[i], &tmplt, false);
    
  			CalculateCorralation(manutmplt, tmplt, &ccv);
			float ncc = pixval32f(ccv, 0, 0);
// 			float ncc = GetCorrelationScore(cpy, fids, tmplt);
			if(tmplt.size().width > 16){
				ccvec.push_back(ncc);//*positive_value);//(positive_value > 0.7 ? 1 : sqrt(positive_value)));
			}
			else{
				ccvec.push_back(ncc*positive_value);
			}
			
// 			util::SaveImage(tmplt, "tmplt.pgm");
			std::cout<<"diameter(value "<<diavec[i]<<"): cc score = "<<ncc<<", acceptive ratio = "<<positive_value<<std::endl;
			
			tmpltvec.push_back(tmplt);
		}
		
		int max_idx = 0;
		float max_cc = ccvec[0];
		for(int i = 1; i < ccvec.size(); i++){
			if(ccvec[i] > max_cc){
				max_idx = i;
				max_cc = ccvec[i];
			}
		}
		
		if(max_idx > 0 && max_idx < diavec.size()-1){
			dbdia_begin = diavec[max_idx-1];
			dbdia_end = diavec[max_idx+1];
		}
		else if(max_idx == 0){
			dbdia_begin = diavec[0];
			dbdia_end = diavec[1];
		}
		else{
			dbdia_begin = diavec[diavec.size()-2];
			dbdia_end = diavec[diavec.size()-1];
		}
		
		if(dbdia_end-dbdia_begin < OUTPUT_THRESHOLD){
			*avgtmplt = tmpltvec[max_idx];
// 			for(int i = 0; i < max_idx-1; i++){
// 				//cvReleaseImage(&(tmpltvec[i]));
// 			}
// 			for(int i = max_idx+1; i < tmpltvec.size(); i++){
// 				//cvReleaseImage(&(tmpltvec[i]));
// 			}
			break;
		}
		
// 		for(int i = 0; i < tmpltvec.size(); i++){
// 			//cvReleaseImage(&(tmpltvec[i]));
// 		}
 	}
	
	*refined_diameter = (dbdia_begin+dbdia_end)*.5;
	
	std::cout<<"Totally sampling number will be "<<int(img.size().width*img.size().height/ (*refined_diameter* *refined_diameter)*.5)<<std::endl;
 	EX_TIME_END("Refined fiducial marker diameter: value = %.2f", *refined_diameter)
	
#undef OUTPUT_THRESHOLD
}

void Detector::SetFidDiameter(float __fid_diameter)
{
	fid_diameter = __fid_diameter;
}

void Detector::SetSampling(cv::Mat& img, float __fid_diameter)
{
	sampling = int(img.size().width*img.size().height/(__fid_diameter*__fid_diameter)*.5);
}

void Detector::SetFidAverageTemplate(const cv::Mat& __fid_tmplt)
{
    //fid_tmplt.copyTo(__fid_tmplt);
    fid_tmplt=__fid_tmplt.clone();
	//fid_tmplt = cvCloneImage(__fid_tmplt);
}

void Detector::Process(cv::Mat& img, std::vector< util::point2d >& fids, bool use_avgtmplt)
{
	EX_TIME_BEGIN("\nCompute Positions of Fiducial Markers (X Y)")

	cv::Mat tmplt;
	if(!use_avgtmplt){
		CreateTemplate(&tmplt, fid_diameter);
	}
	else{
        tmplt=fid_avgtmplt.clone();
	}
	
	GaussianSmoothBasedOnSampling(img, sampling);
	float p_v;
	FindFiducialMarkerPositions(img, tmplt, fid_diameter, fids, &p_v, false, true);
	
    EX_TIME_END("Totally found %ld fiducial markers", fids.size())
}

float Detector::FindLocalPeak(cv::Mat& img, const cv::Mat& tmplt, const util::point2d& seed, float roi_radio, util::point2d* loc)
{
	cv::Rect rect;
	rect.x = seed.x-roi_radio;
	rect.y = seed.y-roi_radio;
	rect.width = roi_radio*2;
	rect.height = roi_radio*2;
	
	cv::Mat reg = util::GetSubImage(img, rect);
	if(reg.empty()){
		return -1;
	}
	
	cv::Mat corr;
	CalculateCorralation(reg, tmplt, &corr);
// 	util::ConvertTo1(corr, false);//true);
	
	int peak_x, peak_y;
	FindMaxPeak(corr, 0, 0, corr.size().width, &peak_x, &peak_y);
// 	float max_value = CV_IMAGE_ELEM(corr, float, peak_y, peak_x);
    float max_value = corr.at<float>(peak_y, peak_y);
	
	float radius = tmplt.size().width*.5;
	loc->x = rect.x+peak_x+radius;
	loc->y = rect.y+peak_y+radius;
	
	util::point2d rloc;
	
	bool success = GetCenter(img, *loc, (float)tmplt.size().width, rloc);
	*loc = rloc;
	
	if(success){
		return max_value;
	}
	else{
		return -1;
	}
}

void Detector::InitLocalDetector(float diameter, bool use_avgtmplt)
{
	if(use_avgtmplt){
		util::SeriesReadFromFile(&detector.fid_avgtmplt, "avgtmplt");
		
        detector.fid_tmplt = detector.fid_avgtmplt.clone();
		//detector.fid_tmplt = cvCloneImage(detector.fid_avgtmplt);
		
		detector.SetFidDiameter(detector.fid_avgtmplt.size().width);				//replace the fid_diameter in detector by average template
	}
	else{
		detector.SetFidDiameter(diameter);
		detector.CreateTemplate(&(detector.fid_avgtmplt), diameter);
	}
}

float Detector::LocalMarkerDetect(cv::Mat& img, const util::point2d& seed, util::point2d* loc)
{
#define ROI_RADIUS_RATIO		.6f
// 	std::cout<<detector.fid_tmplt<<std::endl;
	return detector.FindLocalPeak(img, detector.fid_tmplt, seed, detector.fid_diameter*ROI_RADIUS_RATIO, loc);
}


void Detector::PredictMissFiducial(const util::point2d& seed, HSetVector& hsets, int z_idx1, int z_idx2, std::vector< util::point2d >& slocs)
{
	slocs.clear();
	h_set* txset = hsets.GetHSetWithIdx(z_idx1, z_idx2);
			
	if(txset){
		for(int m = 0; m < txset->size(); m++){
			cv::Mat txm = txset->h[m];
			cv::Point2f tmp, p;
//             CvPoint2D32f tmp, p;
			p.x = seed.x; p.y = seed.y;
			util::_point tmpp;
			tmp = persp_xform_pt(p, txm);
			tmpp.x = tmp.x; tmpp.y = tmp.y;
			slocs.push_back(tmpp);
		}
	}
	else{
		txset = hsets.GetHSetWithIdx(z_idx2, z_idx1);
		
		if(!txset){
			return;
		}
		
		for(int m = 0; m < txset->size(); m++){
			cv::Mat txm = txset->h[m];
			cv::Mat txinv = cv::Mat::zeros(3, 3, CV_64FC1);
			cv::invert(txm, txinv);
			cv::Point2f tmp, p;
			p.x = seed.x; p.y = seed.y;
			util::_point tmpp;
			tmp = persp_xform_pt(p, txinv);
			tmpp.x = tmp.x; tmpp.y = tmp.y;
			slocs.push_back(tmpp);
		}
	}
}

void Detector::LocalDetectorMain(util::MrcStack& mrcs, util::TrackSpace& trackspace, HSetVector& hsets, float diameter, util::FiducialStack* fidsk, util::ImgMatchVector* imvector)
{
	fidsk->ReSize(mrcs.Size());
    fidsk->SetWxH(mrcs.Width(), mrcs.Height());
	
	int zidx = trackspace.Size()/2;
	Detector::InitLocalDetector(diameter);
	
	std::vector<LocalSeed>* lsearchinfos = new std::vector<LocalSeed>[trackspace.Size()];
 	
	for(util::TrackSpace::Iterator itr = trackspace.Z_Iterator(zidx); !itr.IsNULL(); itr++){
		util::TrackSpace::TrackNode node_itr = util::TrackSpace::TrackNode(itr);
		
		for(; !node_itr.IsNULL(); node_itr++){
			util::TrackSpace::TrackNode tmp = node_itr;
			tmp++;

			if(tmp.IsNULL()){
				if(node_itr.Z()+1 == trackspace.Size()){
					continue;
				}
				LocalSeed lsed;
				lsed.z_bp.x = node_itr.X(); lsed.z_bp.y = node_itr.Y();
				lsed.z_begin = node_itr.Z(); lsed.z_next = -1;
				lsed.is_final = false;
				lsed.z_img_idx = node_itr.Z()+1;
				lsearchinfos[lsed.z_img_idx].push_back(lsed);
			}
			else if(tmp.Z() != node_itr.Z()+1){
				LocalSeed lsed;
				lsed.z_bp.x = node_itr.X(); lsed.z_bp.y = node_itr.Y();
				lsed.z_np.x = tmp.X(); lsed.z_np.y = tmp.Y(); 
				lsed.z_begin = node_itr.Z(); lsed.z_next = tmp.Z();
				if(lsed.z_begin+2 != lsed.z_next){
					lsed.is_final = false;
				}
				else{
					lsed.is_final = true;
				}
				lsed.z_img_idx = node_itr.Z()+1;
				
				lsearchinfos[lsed.z_img_idx].push_back(lsed);
			}
		}
		
		node_itr = util::TrackSpace::TrackNode(itr);
		
		for(; !node_itr.IsNULL(); node_itr--){
			util::TrackSpace::TrackNode tmp = node_itr;
			tmp--;

			if(tmp.IsNULL()){
				if(node_itr.Z() == 0){
					continue;
				}
				LocalSeed lsed;
				lsed.z_bp.x = node_itr.X(); lsed.z_bp.y = node_itr.Y();
				lsed.z_begin = node_itr.Z(); lsed.z_next = -1;
				lsed.is_final = false;
				lsed.z_img_idx = node_itr.Z()-1;
				lsearchinfos[lsed.z_img_idx].push_back(lsed);
			}
			
			else if(tmp.Z() != node_itr.Z()-1){
				LocalSeed lsed;
				lsed.z_bp.x = node_itr.X(); lsed.z_bp.y = node_itr.Y();
				lsed.z_np.x = tmp.X(); lsed.z_np.y = tmp.Y(); 
				lsed.z_begin = node_itr.Z(); lsed.z_next = tmp.Z();
				if(lsed.z_begin-2 != lsed.z_next){
					lsed.is_final = false;
				}
				else{
					lsed.is_final = true;
				}
				lsed.z_img_idx = node_itr.Z()-1;
				
				lsearchinfos[lsed.z_img_idx].push_back(lsed);
			}
		}
	}
 	
	for(int i = 0/*0 zidx+1*/; i < trackspace.Size(); i++){
		cv::Mat slice = mrcs.GetStackImage(i);
		util::Reversal(slice);
		util::ConvertTo1(slice);
		util::MedianSmooth(slice);
		
		for(int j = 0; j < lsearchinfos[i].size(); j++){
			LocalSeed& lsed = lsearchinfos[i][j];
			
			std::vector<util::point2d> slocs;
			
			PredictMissFiducial(lsed.z_bp, hsets, lsed.z_begin, lsed.z_img_idx, slocs);
			
			h_set* txset = hsets.GetHSetWithIdx(lsed.z_begin, lsed.z_img_idx);
			
			util::point2d nxloc;
			float max_score = -999;
			for(int m = 0; m < slocs.size(); m++){
				util::point2d p;
				float score = LocalMarkerDetect(slice, slocs[m], &p);
				if(score > max_score){
					max_score = score;
					nxloc = p;
				}
			}
			
			if(max_score > 0.85){								//WARNING threshold
				fidsk->V(i).push_back(nxloc);
			}
			else{
				continue;
			}
            
            util::img_match* imatch;
			bool noex;
			
			if(lsed.is_final){
				noex = true;
				imatch = imvector->GetMatchSetWithIdx(lsed.z_begin, lsed.z_img_idx, noex);
				if(!imatch){
					imvector->MallocNewMatch();
					imatch = &((*imvector)[imvector->Size()-1]);
					imatch->idx1 = lsed.z_begin; imatch->idx2 = lsed.z_img_idx; 
				}
				if(noex){
					imatch->pairs.push_back(std::make_pair(lsed.z_bp, nxloc));
				}
				else{
					imatch->pairs.push_back(std::make_pair(nxloc, lsed.z_bp));
				}
				
				noex = true;
				imatch = imvector->GetMatchSetWithIdx(lsed.z_img_idx, lsed.z_next, noex);
				if(!imatch){
					imvector->MallocNewMatch();
					imatch = &((*imvector)[imvector->Size()-1]);
					imatch->idx1 = lsed.z_img_idx; imatch->idx2 = lsed.z_next; 
				}
				if(noex){
					imatch->pairs.push_back(std::make_pair(nxloc, lsed.z_np));
				}
				else{
					imatch->pairs.push_back(std::make_pair(lsed.z_np, nxloc));
				}
			}
			else if(lsed.z_img_idx == trackspace.Size()-1 || lsed.z_img_idx == 0){
				noex = true;
				imatch = imvector->GetMatchSetWithIdx(lsed.z_begin, lsed.z_img_idx, noex);
				if(!imatch){
					imvector->MallocNewMatch();
					imatch = &((*imvector)[imvector->Size()-1]);
					imatch->idx1 = lsed.z_begin; imatch->idx2 = lsed.z_img_idx; 
				}
				if(noex){
					imatch->pairs.push_back(std::make_pair(lsed.z_bp, nxloc));
				}
				else{
					imatch->pairs.push_back(std::make_pair(nxloc, lsed.z_bp));
				}
				
				imatch->pairs.push_back(std::make_pair(lsed.z_bp, nxloc));
			}
			else{
				float sign = 1;
				
				if(lsed.z_begin < zidx){
					sign = -1;
				}
				
				for(int step = 1; step <= 2 && lsed.z_img_idx+step < trackspace.Size() && lsed.z_img_idx-step >= 0; step++){				
					PredictMissFiducial(nxloc, hsets, lsed.z_img_idx, lsed.z_img_idx+step*sign, slocs);
					std::vector<util::point2d> candidates;
					for(int m = 0; m < slocs.size(); m++){
						for(util::TrackSpace::Iterator itr = trackspace.Z_Iterator(lsed.z_img_idx+step*sign); !itr.IsNULL(); itr++){
							float xdelt = slocs[m].x-itr.X(), ydelt = slocs[m].y-itr.Y();
							if(xdelt*xdelt+ydelt*ydelt < 1.5*1.5){
								candidates.push_back(util::point2d(itr.X(), itr.Y()));
								break;
							}
						}
					}

					if(candidates.size() == 1){
						noex = true;
						imatch = imvector->GetMatchSetWithIdx(lsed.z_begin, lsed.z_img_idx, noex);
						if(!imatch){
							imvector->MallocNewMatch();
							imatch = &((*imvector)[imvector->Size()-1]);
							imatch->idx1 = lsed.z_begin; imatch->idx2 = lsed.z_img_idx; 
						}
						if(noex){
							imatch->pairs.push_back(std::make_pair(lsed.z_bp, nxloc));
						}
						else{
							imatch->pairs.push_back(std::make_pair(nxloc, lsed.z_bp));
						}
						
						noex = true;
						imatch = imvector->GetMatchSetWithIdx(lsed.z_img_idx, lsed.z_img_idx+step*sign, noex);
						if(!imatch){
							imvector->MallocNewMatch();
							imatch = &((*imvector)[imvector->Size()-1]);
							imatch->idx1 = lsed.z_img_idx; imatch->idx2 = lsed.z_img_idx+step*sign; 
						}
						if(noex){
							imatch->pairs.push_back(std::make_pair(nxloc, candidates[0]));
						}
						else{
							imatch->pairs.push_back(std::make_pair(candidates[0], nxloc));
						}
					}
				}
			}
		}
		
// 		cvReleaseImage(&slice);
	}
        delete []lsearchinfos;
}

void Detector::DetectorMain(util::MrcStack& mrcr, util::FiducialStack* fidsk, float& diameter, float ratio)
{
    fidsk->ReSize(mrcr.Size());
 	fidsk->SetRatio(ratio);
    fidsk->SetWxH(mrcr.Width()*ratio, mrcr.Height()*ratio);
 	diameter *= ratio;

 	EX_TRACE("MRC Image rescale ratio: %.2f\n", ratio)
    	std::cout<<"ok"<<std::endl;
 	cv::Mat mid_slice = mrcr.GetStackImage(mrcr.Size()/2)-1;
// 	cv::Mat mid_slice = mrcr.GetStackImage(0);
 	util::ScaleImage(mid_slice, ratio);
 	util::Reversal(mid_slice);
 	util::ConvertTo1(mid_slice, false);
 	util::MedianSmooth(mid_slice);
// 	std::cout<<mid_slice<<std::endl;
// 	
 	float refined_diameter;
 	cv::Mat avgtmplt;
 	detector.DetermineDiameter(mid_slice, diameter, &refined_diameter, &avgtmplt);
 	detector.SetFidDiameter(refined_diameter);
 	detector.SetSampling(mid_slice, refined_diameter);
 	detector.SetFidAverageTemplate(avgtmplt);
 	diameter = refined_diameter;
 	detector.fid_avgtmplt = avgtmplt.clone();

 	util::ConvertTo1(detector.fid_avgtmplt, false);
 	
 	util::SeriesSaveToFile(avgtmplt, "avgtmplt");
 	util::Reversal(avgtmplt);
 	util::ConvertTo1(avgtmplt, false);
 	util::SaveImage(avgtmplt, "avgtmplt(reversal).pgm");
// 	std::cout<<avgtmplt<<std::endl;
    
// 	for(int i = 12; i < mrcr.Size()-12; i++){
    for(int i = 0; i < mrcr.Size(); i++){
        EX_TIME_BEGIN("%sProcessing MRC[%d]", _DASH, i)
		
        cv::Mat slice = mrcr.GetStackImage(i);
		util::ScaleImage(slice, ratio);
		util::Reversal(slice);
        util::ConvertTo1(slice, false);
        util::MedianSmooth(slice);
//         util::HistogramStretch(slice);

// 	std::stringstream ss;
// 	ss<<"histo"<<i<<".pgm";
// 	util::SaveImage(slice, ss.str().c_str());
// std::cout << slice << std::endl;
        std::vector<util::point2d>& fids = fidsk->V(i);
        detector.Process(slice, fids, true);

        EX_TIME_END("Processing MRC[%d]", i)
    }

}

void Detector::Test(util::MrcStack& mrcr, const util::FiducialStack& fidsk, float diameter, const char* folder)
{
    EX_TIME_BEGIN("%s\nDetector Testing", _DASH)
	
	if(diameter < 0){
		diameter = detector.fid_diameter;
	}

    for(int i = 0; i < mrcr.Size(); i++){
        cv::Mat slice = mrcr.GetStackImage(i);

        cv::Mat init = cv::Mat::zeros(cv::Size(slice.size().width*fidsk.Ratio(), slice.size().height*fidsk.Ratio()), CV_32FC1);
        cv::resize(slice, init, cv::Size(slice.size().width*fidsk.Ratio(), slice.size().height*fidsk.Ratio()), 0, 0, cv::INTER_CUBIC);
        slice = init;
        util::ConvertTo1(slice, true);

        const std::vector<util::point2d>& fids = fidsk.V(i);

		diameter = diameter;//*fidsk.Ratio();			!WARNING
        DrawFiducialMarkerPositions(slice, diameter, fids);

        if(access(folder,0) == -1){		//create file folder
            mkdir(folder,0777);
        }
        std::ostringstream oss;
		 oss <<"/home/xzh/文档/markerauto/mk_gai/build_tt/bin/fids"<<"("<<i<<").pgm";
//         oss <<folder<<"/"<<mrcr.Name()<<"("<<i<<").pgm";
        try {
            util::SaveImage(slice, oss.str().c_str());
        } catch(ex::Exception& e){
            EX_TRACE("%s\n", e.Msg())
        }
    }
    EX_TIME_END("Detector Testing")
}



















