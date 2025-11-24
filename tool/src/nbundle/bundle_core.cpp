#include "bundle_core.h"
#include "params.h"
#include "compiler.h"
#include <math.h>
#include "dataf/calibration.h"
#include <vector>
#include "matrix/matrix.h"
#include "micros.h"
#include "util/exception.h"
#include "triangulate.h"
#include <unistd.h>

#ifndef ABS
#define ABS(x) ( ( (x) < 0 )? -(x) : (x) )
#endif

using namespace bundle;

void PPBundleApp::InitGeometryData(util::TrackSpace& trackspace)
{
    std::vector<ImageData>& image_data = data.m_image_data;
    std::vector<TrackData>& track_data = data.m_track_data;
    MatchTable& matches = data.m_matches;

    image_data.clear();
    for(int i = 0; i < trackspace.Size(); i++){
        ImageData data;
        data.m_angle_in_mrc = trackspace.Z_Angle(i);
        image_data.push_back(data);
    }
    // Create the match table
    matches = MatchTable(trackspace.Size());

    for(int i = 0; i < trackspace.Size(); i++){
        util::TrackSpace::Iterator itr = trackspace.Z_Iterator(i);
        int num = 0;
        while(!itr.IsNULL()){
            image_data[i].AddKey(itr.X(), itr.Y());
            itr.SetExtraAsZVectorIndex(num);
            num++;
            itr++;
        }
    }

    for(int i = 0; i < trackspace.Size(); i++){
        util::TrackSpace::Iterator itr = trackspace.Z_Iterator(i);

        while(!itr.IsNULL()){
            TrackData new_track;

            util::TrackSpace::TrackNode node_itr = util::TrackSpace::TrackNode(itr);

            if(node_itr.IsBegin()){
                std::vector<ImageKey>& key_v = new_track.m_views;

                for(; !node_itr.IsNULL(); node_itr++){
                    key_v.push_back(ImageKey(node_itr.Z(), (int)node_itr.VecIndex()));
                }

                if(key_v.size() >=2){
                    track_data.push_back(new_track);
                }
// 		EX_TRACE("Point with %d projections found\n", (int)key_v.size());

                key_v.clear();
            }

            itr++;
        }
    }

    for(int i = 0; i < track_data.size(); i++){
        int num_features = (int) track_data[i].m_views.size();

        for (int j = 0; j < num_features; j++){
            int img1 = track_data[i].m_views[j].first;
            int key1 = track_data[i].m_views[j].second;

            image_data[img1].m_visible_points.push_back(i);
            image_data[img1].m_visible_keys.push_back(key1);
            image_data[img1].m_keys[key1].m_track = i;
        }
    }

    /* Set match flags and full fill the matchtable*/
    for(int i = 0; i < track_data.size(); i++){
        TrackData &track = track_data[i];
        int num_views = (int) track.m_views.size();

        for (int j = 0; j < num_views; j++){
            int img1 = track.m_views[j].first;
            int k1 = track.m_views[j].second;

            assert(img1 >= 0 && img1 < image_data.size());

            for (int k = j+1; k < num_views; k++){
                int img2 = track.m_views[k].first;
                int k2 = track.m_views[k].second;

                assert(img2 >= 0 && img2 < image_data.size());

                matches.SetMatch(MatchIndex(img1, img2));
                matches.GetMatchList(MatchIndex(img1, img2)).push_back(KeypointMatch(k1, k2));

                matches.SetMatch(MatchIndex(img2, img1));
                matches.GetMatchList(MatchIndex(img2, img1)).push_back(KeypointMatch(k2, k1));
            }
        }
    }

    EX_NOTE(PrintTrackBin();)

    trackspace.Release();
}

void PPBundleApp::InitBoundary(int _width, int _height)
{
    width = _width;
    height = _height;
}

PPBundleApp::PPBundleApp(){}

PPBundleApp::~PPBundleApp()
{
    delete [] result;
}

void PPBundleApp::InitMotParam(float angle, double* params) const
{
	PInitating::InitSR123T12(params, 1.0/PInitating::f, 0, angle, 0, 0, 0);
}

void PPBundleApp::InitCameraParam(float angle, mx::pproj_params* camera) const
{
    double mot[6];
    InitMotParam(angle, mot);
	mx::PProjPCopyFormMotParams(camera, mot);
}

void PPBundleApp::InitCameraParam(double alpha, double beta, double gamma, double t0, double t1, mx::pproj_params* camera)
{
	double params[6];
	PInitating::InitSR123T12(params, 1.0/PInitating::f, alpha, beta, gamma, t0, t1);
	mx::PProjPCopyFormMotParams(camera, params);
}


v3_t PPBundleApp::Triangulate(int num_points, const v2_t* pv, const double* mots, double* error_out)
{
	int num_eqs = 2*num_points;
    int num_vars = 3;

    double *As = new double[num_eqs*num_vars];
    double *bs = new double[num_eqs];
    double *x = new double[num_vars];

    v3_t r;

    for(int i = 0; i < num_points; i++){
		const double* mot = mots+6*i;
		double s, alpha, beta, gamma, t0, t1;
		s = mot[0];
		alpha = mot[1];
		beta = mot[2];
		gamma = mot[3];
		t0 = mot[4];
		t1 = mot[5];

		double cos_alpha = cos(alpha);
		double sin_alpha = sin(alpha);
		double cos_beta = cos(beta);
		double sin_beta = sin(beta);
		double cos_gamma = cos(gamma);
		double sin_gamma = sin(gamma);

        double* A = As+6*i;
        double* b = bs+2*i;
        A[0] = cos_beta;
        A[1] = sin_alpha*sin_beta;
        A[2] = -cos_alpha*sin_beta;
        A[3] = 0;
        A[4] = cos_alpha;
        A[5] = sin_alpha;

        b[0] = s*(cos_gamma*Vx(pv[i])+sin_gamma*Vy(pv[i])+t0);
        b[1] = s*(-sin_gamma*Vx(pv[i])+cos_gamma*Vy(pv[i])+t1);
    }

    /* Find the least squares result */
    dgelsy_driver(As, bs, x, num_eqs, num_vars, 1);

    if(error_out != NULL){
        double error = 0.0;
        for(int i = 0; i < num_points; i++){
            double dx, dy;
            double pp[2];
			PSFMDriver::Project(mots+6*i, x, pp);

            dx = pp[0] - Vx(pv[i]);
            dy = pp[1] - Vy(pv[i]);
            error += dx*dx+dy*dy;
        }
        error = sqrt(error / num_points);

        *error_out = error;
    }

    r = v3_new(x[0], x[1], x[2]);

    delete [] As;
    delete [] bs;
    delete [] x;

    return r;
}

v3_t PPBundleApp::Triangulate(v2_t p, v2_t q, const mx::pproj_params& c1, const mx::pproj_params& c2, double* proj_error)
{
	v2_t pv[2] = {p, q};
	double mots[12];

	mx::MotCopyFormPProjParams(mots, c1);
	mx::MotCopyFormPProjParams(mots+6, c2);

    return Triangulate(2, pv, mots, proj_error);
}

/** @brief Triangulate a subtrack */
v3_t PPBundleApp::TriangulateNViews(const ImageKeyVector& views, int* added_order, mx::pproj_params* cameras, double& error)
{
    int num_views = (int) views.size();

    v2_t* pv = new v2_t[num_views];
    double* mots = new double[6*num_views];

    int num_valid = 0;

    for(int i = 0; i < num_views; i++){
        mx::pproj_params* cam = NULL;

        int camera_idx = views[i].first;
        int image_idx = added_order[camera_idx];

        int key_idx = views[i].second;
        Keypoint& key = data.GetKey(image_idx, key_idx);

        pv[i] = v2_new(key.m_x, key.m_y);

        cam = cameras + camera_idx;
		mx::MotCopyFormPProjParams(mots+6*i, *cam);
        num_valid++;
    }

    v3_t pt = Triangulate(num_valid, pv, mots, &error);

    delete [] pv;
    delete [] mots;

    return pt;
}

void PPBundleApp::Process(std::vector<mx::pproj_params>& are_cam, float percent)
{
	double sgamma[data.NumImages()];
	double ty[data.NumImages()];
	double tx[data.NumImages()];
	memset(sgamma, 0, sizeof(double)*data.NumImages());
	memset(ty, 0, sizeof(double)*data.NumImages());
	memset(tx, 0, sizeof(double)*data.NumImages());

	EX_PRINT("Gross Estimation of Global Rotation and Y Shift\n")

	GRTxyRefiner grtxyrefiner;

	*sgamma = 0;
	grtxyrefiner.Run(data, percent, sgamma, tx, ty);
// 	EX_PRINT("Finished (Global Rotation: %f)\n",gamma/M_PI*180)

// 	for(int i = 0; i < data.NumImages(); i++){
// 		std::cout<<cos(gamma)<<"\t"<<sin(gamma)<<"\t"<<-sin(gamma)<<"\t"<<cos(gamma)<<"\t"<<tx[i]<<"\t"<<ty[i]<<std::endl;
// 	}
//     for (const auto& alignment : are_cam) {
//         std::cout << "s: " << alignment.s << ", alpha: " << alignment.alpha
//                   << ", beta: " << alignment.beta << ", gamma: " << alignment.gamma<<", TX: "<<alignment.t0<< ", TY: "<<alignment.t1 << std::endl;
//     }

	mx::pproj_params* camparams = new mx::pproj_params[data.NumImages()];		//as the order of added_order
	mx::pproj_params* camparams_are = new mx::pproj_params[data.NumImages()];		//as the order of added_order
	int* added_order = new int[data.NumImages()];
    int* added_order2 = new int[data.NumImages()];
	int max_pts =(int)data.m_track_data.size();
	v3_t* points = new v3_t[max_pts];
    v3_t* points_are = new v3_t[max_pts];
	int pt_count = 0;


    for(int i = 0; i < data.NumImages(); i++){
// 		InitCameraParam(0, data.OriAngles(i), sgamma[i]/M_PI*180, tx[i], ty[i], camparams+i);
//         InitCameraParam(0, are_cam[i].beta, are_cam[i].gamma, are_cam[i].t0, are_cam[i].t1, camparams+i);
        InitCameraParam(0, are_cam[i].beta, 0, 0, 0, camparams+i);
        // std::cout<<"["<<(camparams+i)->s<<", "<<(camparams+i)->alpha<<", "<<(camparams+i)->beta
            // <<", "<<(camparams+i)->gamma<<", "<<(camparams+i)->t0<<", "<<(camparams+i)->t1<<"]"<<std::endl;
		added_order[i] = i;
	}

	for(int i = 0; i < data.NumImages(); i++){
// 		InitCameraParam(0, are_cam[i].beta, are_cam[i].gamma, are_cam[i].t0, are_cam[i].t1, camparams_are+i);
        InitCameraParam(0, are_cam[i].beta, 0, 0, 0, camparams_are+i);
        // std::cout<<"camparams_are:  "<<"["<<(camparams_are+i)->s<<", "<<(camparams_are+i)->alpha<<", "<<(camparams_are+i)->beta
            // <<", "<<(camparams_are+i)->gamma<<", "<<(camparams_are+i)->t0<<", "<<(camparams_are+i)->t1<<"]"<<std::endl;
		added_order2[i] = i;
	}




	std::vector<TrackData>& track_data = data.m_track_data;
	int num_tracks_total =(int)data.m_track_data.size();
	std::vector<ImageKeyVector> pt_views;

	int length_thre = percent*data.NumImages();
	for(int i = 0; i < track_data.size(); i++){
		TrackData& track = track_data[i];
		if(track.m_views.size() < length_thre){
			continue;
		}

		double error;
		v3_t pt = TriangulateNViews(track.m_views, added_order, camparams, error);
        v3_t pt2 = TriangulateNViews(track.m_views, added_order2, camparams_are, error);

// 		std::cout<<error<<"\t"<<track.m_views.size()<<std::endl;

		points[pt_count] = pt;
        points_are[pt_count] = pt2;
        // std::cout<<"***************************************"<<std::endl;
        // std::cout<<pt.p[0]<<"  "<<pt.p[1]<<"  "<<pt.p[2]<<"  "<<std::endl;
        // std::cout<<pt2.p[0]<<"  "<<pt2.p[1]<<"  "<<pt2.p[2]<<"  "<<std::endl;
        // std::cout<<"***************************************"<<std::endl;
		pt_views.push_back(track.m_views);

		for(int j = 0; j < track.m_views.size(); j++){
			data.GetKey(track.m_views[j].first, track.m_views[j].second).m_extra = pt_count;
		}

		track.m_extra = pt_count;

		pt_count++;
	}

	EX_PRINT("%d Structure Points used\n", pt_count)

//     std::string res_file1 = "/home/xzh/桌面/W.txt";
//     std::ofstream outputfile_res1(res_file1);
//
//     for(int i=0;i<pt_count;i++)
//     {
//         outputfile_res1<<points[i].p[0]<<std::endl;
//         outputfile_res1<<points[i].p[1]<<std::endl;
//         outputfile_res1<<points[i].p[2]<<std::endl;
// //         outputfile_res<<"\n";
//     }
//
//     outputfile_res1.close();

// 	for(int i = 0; i < pt_count; i++){
// 		v3_print(points[i]);
// 	}

	driver.SetGlobMask(0, 1, 1, 0, 0, 0, 0);
// 	driver.Run(added_order, camparams, 0, data.NumImages(), 0, points, pt_count, 0, pt_views, data, false);
    driver.Run_are(added_order2, camparams_are, 0, data.NumImages(), 0, points_are, pt_count, 0, pt_views, data, false);

	//re-set centre of the reconstructed point

    double avg=0;

// 	driver.SetGlobMask(0, 0, 0, 0, 0, 0, 0);
// 	avg=driver.Run(added_order, camparams, 0, data.NumImages(), 0, points, pt_count, 0, pt_views, data);
//
//     std::cout<<avg<<std::endl;

    result = new mx::pproj_params[data.NumImages()];

    for(int i = 0; i < data.NumImages(); i++){
        result[i].valid = false;
    }

    for(int i = 0; i < data.NumImages(); i++){
        int img = added_order[i];
        if(!data.m_image_data[img].m_ignore_in_bundle){
            result[img] = camparams[i];
            result[img].valid = true;
        }
        else {
            result[img].valid = false;
        }
    }

    EX_NOTE(PrintRepErrorBin(camparams, data.NumImages(), added_order, points);)

    /* Points */
    for(int i = 0; i < pt_count; i++){
        /* Check if the point is visible in any view */
        if((int)pt_views[i].size()== 0)
            continue;/* Invisible */

        PointData pdata;
        pdata.m_pos[0] = Vx(points[i]);
        pdata.m_pos[1] = Vy(points[i]);
        pdata.m_pos[2] = Vz(points[i]);

#ifdef FULL_INFO
        for(int j = 0; j <(int)pt_views[i].size(); j++){
            int v = pt_views[i][j].first;
            int vnew = added_order[v];
            pdata.m_views.push_back(ImageKey(vnew, pt_views[i][j].second));
        }
#endif
        m_point_data.push_back(pdata);
    }

    delete [] added_order;
    delete [] camparams;
    delete [] points;
//     std::ofstream out("camera_C.txt");
//     for(int i = 0; i < data.NumImages(); i++){
// 	PrintMxCofCamera(camparams[i], out);
//     }
//     out.close();
}

void PPBundleApp::PrintCameras(const char* filename) const
{
    std::ofstream out(filename);
    for(int i = 0; i < data.NumImages(); i++){
        if(result[i].valid){
            out<<"Camera ["<<i<<"]\n";
            out<<result[i].s<<"\t"<<result[i].alpha<<"\t"<<result[i].beta<<"\t"<<result[i].gamma<<"\t"<<result[i].t0<<"\t"<<result[i].t1<<std::endl;
        }
        else {
            out<<"Camera ["<<i<<"] Invalid"<<std::endl;
        }
    }

    out.close();
}

void PPBundleApp::PrintPoints(const char* filename) const
{
    std::ofstream out(filename);
    for(int i = 0; i < m_point_data.size(); i++){
        out<<m_point_data[i].m_pos[0]<<" "<<m_point_data[i].m_pos[1]<<" "<<m_point_data[i].m_pos[2]<<std::endl;
    }

    out.close();
}

void PPBundleApp::PrintTrackBin(const char* filename) const
{
    int bin[data.NumImages()+1];
    memset(bin, 0, sizeof(int)*(data.NumImages()+1));
    for(int i = 0; i < data.m_track_data.size(); i++){
        bin[data.m_track_data[i].m_views.size()]++;
    }
    std::ofstream out(filename);
    for(int i = data.NumImages(); i >= 2; i--){
        out<<i<<"\t"<<bin[i]<<std::endl;
    }
    out.close();
}

void PPBundleApp::PrintRepErrorBin(mx::pproj_params* camparams, int ncams, const int added_order[],
                                 const v3_t* points, const char* foldername) const
{
    if(access(foldername,0) == -1){		//create file folder
        mkdir(foldername,0777);
    }
    std::ostringstream ooo;
    ooo <<foldername<<"/global.hist";
    std::ofstream out(ooo.str().c_str());

    int hist[81];
    float hist_gap = 0.05;
    int total = 0;
    memset(hist, 0, sizeof(int)*81);

    for(int i = 0; i < ncams; i++){
        double mot[6];
		mx::MotCopyFormPProjParams(mot, camparams[i]);

        int num_keys = data.GetNumKeys(added_order[i]);

        int num_pts_proj = 0;
        for(int j = 0; j < num_keys; j++){
            if(data.GetKey(added_order[i], j).m_extra >= 0){
                num_pts_proj++;
            }
        }

        double *dists = new double[num_pts_proj];
        int pt_count = 0;

        std::vector<Keypoint>::const_iterator iter;

        for(iter = data.m_image_data[added_order[i]].m_keys.begin(); iter != data.m_image_data[added_order[i]].m_keys.end(); iter++){
            const Keypoint &key = *iter;

            if(key.m_extra >= 0){
                int pt_idx = key.m_extra;
                double pr[2];

				PSFMDriver::Project(mot, points[pt_idx].p, pr);

                double dx = pr[0]-key.m_x;
                double dy = pr[1]-key.m_y;

                double dist = sqrt(dx * dx + dy * dy);
                int glo_num = round(dist/hist_gap) <= 80 ? round(dist/hist_gap) : 80;
                hist[glo_num]++;
                total++;
// 		global_dist.push_back(dist);

                dists[pt_count] = dist;

                pt_count++;
            }
        }

#define NUM_ERROR_BINS 10
        std::sort(dists, dists+num_pts_proj);

        double pr_min = dists[0];
        double pr_max = dists[num_pts_proj-1];
        double pr_step =(pr_max - pr_min)/NUM_ERROR_BINS;

        /* Break histogram into 10 bins */
        std::ostringstream oss;
        oss <<foldername<<"/cam"<<added_order[i]<<".hist";
        std::ofstream o(oss.str().c_str());

        int idx_count = 0;
        for(int i = 0; i < NUM_ERROR_BINS; i++){
            double max = pr_min +(i+1)* pr_step;
            int start = idx_count;

            while(idx_count < num_pts_proj && dists[idx_count] <= max){
                idx_count++;
            }

            int bin_size = idx_count - start;
            o<<"   ["<<max-pr_step<<"~"<<max<<"]: \t"<<bin_size<<"\t"<<bin_size/(double)num_pts_proj<<std::endl;
        }
        o.close();

        delete [] dists;
    }

    for(int i = 0; i < 81; i++){
        out<<"   ["<<i*hist_gap<<"]: \t"<<hist[i]<<"\t"<<hist[i]/(double)total<<std::endl;
    }
    out.close();
}

void PPBundleApp::PrintCamerasAsIMOD(const std::vector<mx::pproj_params>& cameras, const float rotation, const float scale, const char* xffilename,
									 const char* xanglefilename, const char* yanglefilename, const char* exfilename)
{
    std::ofstream out(xffilename);
	std::ofstream xagout(xanglefilename);
    std::ofstream yagout(yanglefilename);
    std::ofstream exceptf(exfilename);
    // std::cout<<" cameras.size():"<<cameras.size()<<std::endl;
    for(int i = 0; i < cameras.size(); i++){
        // std::cout << "[" << cameras[i].s << ", " << cameras[i].alpha << ", " << cameras[i].beta << ", " << cameras[i].gamma << ", " << cameras[i].t0 << ", " << cameras[i].t1 << "]" << std::endl;

//         cameras[i].beta = cameras[i].beta * 180/M_PI*180;

        double cos_x_rot, sin_x_rot;
        cos_x_rot = cos(cameras[i].gamma);
        sin_x_rot = sin(cameras[i].gamma);
        double cos_y_rot, sin_y_rot;
        cos_y_rot = cos(cameras[i].gamma);
        sin_y_rot = sin(cameras[i].gamma);
        double A[4] = {cos_x_rot, sin_y_rot, -sin_x_rot, cos_y_rot};

		double rscale = cameras[i].s*PInitating::f;


//         std::cout<<"rotation: "<<rotation<<std::endl;

//         if(cameras[i].valid) {
        if(1) {
            out<<(A[0]*cos(rotation)+A[1]*sin(rotation))*rscale
               <<"\t"<<(-A[0]*sin(rotation)+A[1]*cos(rotation))*rscale
               <<"\t"<<(A[2]*cos(rotation)+A[3]*sin(rotation))*rscale
               <<"\t"<<(-A[2]*sin(rotation)+A[3]*cos(rotation))*rscale
               <<"\t"<<cameras[i].t0*scale*rscale<<"\t"<<cameras[i].t1*scale*rscale<<std::endl;
        }
        else {
            out<<1<<"\t"<<0<<"\t"<<0<<"\t"<<1<<"\t"<<0<<"\t"<<0<<std::endl;
            exceptf<<i<<" "<<std::endl;
        }
        xagout<<cameras[i].alpha/M_PI*180<<std::endl;
        yagout<<cameras[i].beta/M_PI*180<<std::endl;
    }

    out.close();
    yagout.close();
	exceptf.close();
}

void PPBundleApp::DumpOutCameras(std::vector< mx::pproj_params >* cameras) const
{
    for(int i = 0; i < data.NumImages(); i++){
        cameras->push_back(result[i]);
    }
}

void PPBundleApp::DumpOutPoints(std::vector< v3_t >* points) const
{
    for(int i = 0; i < m_point_data.size(); i++){
        v3_t tmp;
        memcpy(tmp.p, m_point_data[i].m_pos, sizeof(double)*3);

        points->push_back(tmp);
    }
}

void PPBundleApp::BundleMain(util::TrackSpace& tspace_, int width_, int height_, std::vector<mx::pproj_params>& are_cam, std::vector<mx::pproj_params>* cameras, std::vector<v3_t>* points)
{
    PPBundleApp app;
    app.InitGeometryData(tspace_);
    app.InitBoundary(width_, height_);
    app.Process(are_cam, 0.7);
    app.DumpOutCameras(cameras);
    app.DumpOutPoints(points);
    app.PrintCameras();
    app.PrintPoints();
}

void PPBundleApp::Test()
{

}
