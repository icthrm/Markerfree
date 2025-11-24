#include "match_core.h"
#include <opencv2/opencv.hpp>
#include <iomanip>
#include <sys/stat.h>
#include <unistd.h>
#include <levmar.h>
#include "matrix/matrix.h"

void ModelMatch::PairMatch(std::vector< util::point2d >& __fids1, std::vector< util::point2d >& __fids2, double beta1, double beta2, double dist_err_tol, 
						   std::vector<std::pair<util::point2d, util::point2d> >& matchvec, std::vector<cv::Mat>& hset)
{
	Ran4PEstimator ran4PEstimator(__fids1, __fids2);
	ran4PEstimator.SetSupplementData(beta1, beta2);
	ran4PEstimator.AffineTransformEstimation(dist_err_tol, matchvec, hset, false, true);	//ran4p use half diameter
	
	if(!ran4PEstimator.UmFids1().size()|!ran4PEstimator.UmFids2().size()){		//there are no obvious warp
	}
	else{
		CPDEstimator cpdEstimator(ran4PEstimator.UmFids2(), ran4PEstimator.UmFids1());		//the first is fixed and the second is moving; therefore, it is inverse with 4p method
		cpdEstimator.PointDriftEstimation(dist_err_tol, matchvec, hset, false, false);		//CPD use full diameter
	}
}


void ModelMatch::MatchMain(util::FiducialStack& fstack,const std::vector<float>& angles, util::ImgMatchVector* imvector, HSetVector* hsetvec, float dist_err_tol, bool eigenlimit, bool do_test)
{
	EX_TRACE("Used distance threshold = %.2f\n", dist_err_tol)
	
	#define STEP_ARRAY_SIZE			2
    assert(STEP_ARRAY_SIZE < fstack.Size());

    imvector->Clear();
	hsetvec->Clear();
	
    int step_length[STEP_ARRAY_SIZE];
    for(int i = 0; i < STEP_ARRAY_SIZE; i++){
        step_length[i] = i+1;
    }
	
    int turn = 0;
	
	int sampling = fstack.Width()*fstack.Height()/(dist_err_tol*dist_err_tol*4)*.5;
	
    while(turn < STEP_ARRAY_SIZE){
        for(int i = 0; i+step_length[turn] < fstack.Size(); i++){
            int idx1 = i;//14;//6;//
            int idx2 = i+step_length[turn];//16;//7
            
//             if(!do_test){
// 				if(fstack.V(idx1).size() < fstack.V(idx2).size()){		//this strategy is for Ran4PEstimator
// 					int tmp = idx1;
// 					idx1 = idx2;
// 					idx2 = tmp;
// 				}
// 			}
// 			else{
// 				if(fstack.V(idx1).size() > fstack.V(idx2).size()){		//this strategy is for CPDEstimator
// 					int tmp = idx1;
// 					idx1 = idx2;
// 					idx2 = tmp;
// 				}
// 			}

            if(fstack.V(idx1).size() > fstack.V(idx2).size()){		//this strategy is for Ran4PEstimator
                    int tmp = idx1;
                    idx1 = idx2;
                    idx2 = tmp;
            }
			
            EX_TIME_BEGIN("#\nMatching Point Set (MRC[%d] & MRC[%d])", idx1, idx2)
            util::img_match& imatch = imvector->MallocNewMatch();
            imatch.idx1 = idx1;
            imatch.idx2 = idx2;
			
			h_set& hset = hsetvec->MallocNewHSet();
			hset.idx1 = idx1;
			hset.idx2 = idx2;
            
            PairMatch(fstack.V(idx1), fstack.V(idx2), angles[idx1], angles[idx2], dist_err_tol, imatch.pairs, hset.h);
			
// 			if(do_test){
//                 PairMatch(fstack.V(idx1), fstack.V(idx2), angles[idx1], angles[idx2], dist_err_tol, imatch.pairs, hset.h);
// 				CPDEstimator cpdEstimator(fstack.V(idx1), fstack.V(idx2));		//the first is fixed and th second is moving
// 				cpdEstimator.GlobalMatch(1.7*dist_err_tol, imatch.pairs, hset.h, false, eigenlimit);		//CPD use full diameter
//                 
// 				if(imatch.size() == 0){
// 					Ran4PEstimator ran4PEstimator(fstack.V(idx1), fstack.V(idx2));
// 					ran4PEstimator.GlobalMatch(dist_err_tol, imatch.pairs, hset.h, false, eigenlimit);	//ran4p use half diameter
// 				}
// 				if(sampling < 2500){
// 					EX_TRACE("Low sampling --> Distribution correction is turned on.\n");
// 					cpdEstimator.DistributionCorrection(4*dist_err_tol, imatch.pairs);
// 				}
// 			}
// 			else{
// 				Ran4PEstimator ran4PEstimator(fstack.V(idx1), fstack.V(idx2));
// 				ran4PEstimator.GlobalMatch(dist_err_tol, imatch.pairs, hset.h, false, eigenlimit);	//ran4p use half diameter
// 				if(sampling < 2500){
// 					EX_TRACE("Low sampling --> Distribution correction is turned on.\n");
// 					ran4PEstimator.DistributionCorrection(4*dist_err_tol, imatch.pairs);
// 				}
// 			}
			//if(imatch.size() < 5 || (imatch.size() < fstack.V(idx1).size()*CPDEstimator::COVERAGE_REFRESH_RATIO*.5 &&  
				//imatch.size() < fstack.V(idx2).size()*CPDEstimator::COVERAGE_REFRESH_RATIO*.5) ){
			if(imatch.size() < 5 || (imatch.size() < 9 && imatch.size() < fstack.V(idx1).size()*CPDEstimator::COVERAGE_REFRESH_RATIO*.25 &&  
				imatch.size() < fstack.V(idx2).size()*CPDEstimator::COVERAGE_REFRESH_RATIO*.25) ){
				imatch.pairs.clear();
				hset.h.clear();
			}

            EX_TIME_END("MRC[%d] & MRC[%d]: %ld/(%ld,%ld) Pairs found", idx1, idx2, imatch.size(), fstack.V(idx1).size(), fstack.V(idx2).size())
        }
        turn++;
    }
/*	
	Ran4PEstimator ran4PEstimator(fstack.V(49), fstack.V(50));
	std::vector<std::pair<util::point2d, util::point2d> > matchvec;
	CvMat* H;
	ran4PEstimator.RansacMatch(DIS_ERR_TOL, matchvec, H, true);*/
}

void ModelMatch::DrawMatch(cv::Mat& canvas, const cv::Mat& img1, const cv::Mat& img2, const std::vector<std::pair<util::point2d, util::point2d> >& vpair)
{
//     cv::Zero(canvas);
    canvas.setTo(0);
    cv::Mat roilimage = canvas(cv::Rect(0, 0, img1.size().width, img1.size().height));
//     cvSetImageROI(canvas, cvRect(0, 0, img1->width, img1->height));
    cv::add(img1, roilimage, roilimage, NULL);
    cv::Mat roilimage2 = canvas(cv::Rect(img1.size().width, 0, img1.size().width+img2.size().width, img2.size().height));
//     cvSetImageROI(canvas, cvRect(img1->width, 0, img1->width+img2->width, img2->height));
    cv::add(img2, roilimage2, roilimage2, NULL);
    for(int i = 0; i < vpair.size(); i++) {
        cv::Scalar color = CV_RGB(255, 255, 255);
//         util::DrawX(canvas, vpair[i].first.x, vpair[i].first.y);
//         util::DrawX(canvas, img1->width+vpair[i].second.x, vpair[i].second.y);
		util::DrawLine(canvas, vpair[i].first, util::_point(img1.size().width+vpair[i].second.x, vpair[i].second.y));
    }
}

void ModelMatch::Test(util::MrcStack& mrcr, const util::ImgMatchVector& imvector, float ratio, const char* folder)
{
    EX_TIME_BEGIN("%sMatch Testing", _DASH)

    cv::Mat p[2], canvas, tmp;
//     IplImage* p[2], * canvas, *tmp;
    canvas = cv::Mat::zeros(cv::Size(mrcr.Width()*2*ratio, mrcr.Height()*ratio), CV_32FC1);
//     canvas = cvCreateImage(cvSize(mrcr.Width()*2*ratio, mrcr.Height()*ratio), IPL_DEPTH_32F, 1);

    for(int i = 0; i < imvector.Size(); i++) {
        tmp = mrcr.GetStackImage(imvector[i].idx1);
        p[0] = cv::Mat::zeros(cv::Size(tmp.size().width*ratio, tmp.size().height*ratio), CV_32FC1);
        cv::resize(tmp, p[0], cv::Size(tmp.size().width*ratio, tmp.size().height*ratio), 0 , 0, cv::INTER_CUBIC);
        util::ConvertTo1(p[0], true);

        tmp = mrcr.GetStackImage(imvector[i].idx2);
        p[1] = cv::Mat::zeros(cv::Size(tmp.size().width*ratio, tmp.size().height*ratio), CV_32FC1);
        cv::resize(tmp, p[1], cv::Size(tmp.size().width*ratio, tmp.size().height*ratio), 0 , 0, cv::INTER_CUBIC);
        util::ConvertTo1(p[1], true);

		std::vector<std::pair<util::point2d, util::point2d> > tmppairs;
		for(int j = 0; j < imvector[i].pairs.size(); j++){
			util::point2d pair1, pair2;
			pair1.x = imvector[i].pairs[j].first.x*ratio;
			pair1.y = imvector[i].pairs[j].first.y*ratio;
			pair2.x = imvector[i].pairs[j].second.x*ratio;
			pair2.y = imvector[i].pairs[j].second.y*ratio;
			tmppairs.push_back(std::make_pair(pair1, pair2));
		}
		
        DrawMatch(canvas, p[0], p[1], tmppairs);
        if(access(folder,0) == -1) {		//create file folder
            mkdir(folder,0777);
        }
        std::ostringstream oss;
        oss <<folder<<"/"<<"("<<imvector[i].idx1<<")&("<<imvector[i].idx2<<").pgm";
        try {
            util::SaveImage(canvas, oss.str().c_str());
        } catch(ex::Exception& e) {
            EX_TRACE("%s\n", e.Msg())
        }
    }
    EX_TIME_END("Match Testing")
}

void HSetVector::WriteTransforms(const h_set& hset, std::ostream& out)
{
	out<<hset.idx1<<"\t"<<hset.idx2<<"\t"<<hset.size()<<std::endl;
	for(int i = 0; i < hset.size(); i++){
		for(int m = 0; m < 3; m++){
			for(int n = 0; n < 3; n++){
                out<<hset.h[i].at<double>(m, n)<<"\t";
// 				out<<cvGetReal2D(hset.h[i], m, n)<<"\t";
			}
			out<<std::endl;
		}
		out<<std::endl<<std::endl;
	}
}

void HSetVector::ReadTransforms(h_set* hset, std::istream& in)
{
	int hsize;
	in>>hset->idx1>>hset->idx2>>hsize;
	for(int i = 0; i < hsize; i++){
		cv::Mat h = cv::Mat::zeros(3, 3, CV_64FC1);
		for(int m = 0; m < 3; m++){
			for(int n = 0; n < 3; n++){
				double val;
				in>>val;
                h.at<double>(m, n) = val;
// 				cvmSet(h, m, n, val);
			}
		}
		hset->h.push_back(h);
	}
}

void HSetVector::WriteVectorByFolder(const char* folderpath) const
{
    if(access(folderpath,0) == -1) {		//create file folder
        mkdir(folderpath,0777);
    }
    std::ostringstream ooo;
    ooo <<folderpath<<"/attributes";
    std::ofstream out(ooo.str().c_str());
    out<<"Z:"<<hset_vector->size()<<"\n";
    out.close();
    for(int i = 0; i < hset_vector->size(); i++) {
        std::ostringstream oss;
        oss <<folderpath<<"/"<<i<<".txt";
        std::ofstream o(oss.str().c_str());
        try {
            WriteTransforms((*hset_vector)[i], o);
        } catch(ex::Exception& e) {
            EX_TRACE("%s\n", e.Msg())
        }
    }
}

void HSetVector::ReadVectorByFolder(const char* folderpath)
{
    Clear();
    std::cout <<std::setprecision(8)<<std::endl;
    std::ostringstream iii;
    iii <<folderpath<<"/attributes";
    std::ifstream in(iii.str().c_str());
    char ch;
    int _size;
    in>>ch>>ch>>_size;
    in.close();

    for(int i = 0; i < _size; i++) {
        std::ostringstream oss;
        oss <<folderpath<<"/"<<i<<".txt";
        std::ifstream in(oss.str().c_str());
        if(!in.good()) {
            ex::EX_THROW("Can't Open File");
        }
        h_set& hset = MallocNewHSet();
        ReadTransforms(&hset, in);
        in.close();
    }
}
