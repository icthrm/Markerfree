#include "opts.h"
#include <iostream>
#include <fstream>
#include "dataf/dataf.h"
#include "modelmatch/match_core.h"
#include "detector/detector.h"
#include "nbundle/bundle_core.h"

// #define TESTDEVELOP

using namespace std;
// using namespace ann_1_1_char;

struct Alignment {
    int SEC;
    double ROT;
    double TX;
    double TY;
    double TILT;
};

struct Alignment_are {
    int SEC;
    double ROT;
    double TX;
    double TY;
    double TILT;
    double GMAG;
    double SMEAN;
    double SFIT;
    double SCALE;
    double BASE;
};

int main(int argc, char **argv)
{

//     std::ifstream file("/home/xzh/下载/data/aretomo对齐数据/BBb_align.txt");
    int a=1;
    struct options opts;

    opts.diameter = -1;
    opts.verbose = 0;
    opts.rotation_angle = 0;
    opts.testmode = false;

    if (GetOpts(argc, argv, &opts) <= 0)
    {
        EX_TRACE("***WRONG INPUT.\n");
        return -1;
    }

    std::vector<mx::pproj_params> are_cam;
    if(1){ //liu
        // std::ifstream file("/home/xzh/下载/╩²╛▌/nmar20024_align_de.txt");
        // std::ifstream file("/home/xzh/下载/data/aretomo对齐数据/BBb_align.txt");
        std::ifstream file(opts.txtinput);
        if (!file.is_open())
        {
            std::cerr << "Error: Cannot open file " << opts.txtinput << std::endl;
            return -1;
        }
        //         std::ifstream file("/home/xzh/桌面/去噪/对齐/nmar_001.txt");
        // std::ifstream file("/home/xzh/下载/待测试数据/参数文件/多/10016_align_multi.txt");
        // std::ifstream file("/home/xzh/下载/待测试数据/参数文件/多/nmar.txt");
        // std::ifstream file("/home/xzh/下载/待测试数据/参数文件/多/10007_align_multi.txt");
        //  std::ifstream file("/home/xzh/下载/待测试数据/new/10643/b3tilt51_align.txt");
        //          std::ifstream file("/home/xzh/下载/ETdata/BBb_joint/BBb_200_ali/patch_ali199.txt");
        // std::ifstream file("/home/xzh/桌面/科研项目/无标记点/ver_filter/build/bin/110001/110001_fin.txt");
        // std::ifstream file("/home/xzh/下载/data/aretomo对齐数据/110001_align_tmp.txt");
        // std::ifstream file("/home/xzh/下载/待测试数据/110001_liu.txt");
        //         std::ifstream file("/home/xzh/下载/待测试数据/110001_ali_001.txt");
        // std::ifstream file("/home/xzh/桌面/科研项目/无标记点/sim_result/对齐结果/sim_001.txt");
        // std::ifstream file("/home/xzh/下载/待测试数据/new/10643/b3tilt51_align.txt");
        std::string line;
        std::vector<Alignment> alignments;

//         std::vector<mx::pproj_params> are_cam;

        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') {
                continue;
            }

//             if (line.find("angleoffset") != std::string::npos ||
//             line.find("Output Bin") != std::string::npos ||
//             line.find("Elapsed time") != std::string::npos) {
//                 break;  // 如果遇到这些行，停止读取
//             }

            std::istringstream iss(line);
            Alignment alignment;
            mx::pproj_params are_alignment;
            iss >> alignment.SEC >> alignment.ROT >> alignment.TX >> alignment.TY >> alignment.TILT;
            are_alignment.s=1;
            are_alignment.alpha=0;
            are_alignment.beta=alignment.TILT;
            are_alignment.gamma=alignment.ROT;
            are_alignment.t0=alignment.TX;
            are_alignment.t1=alignment.TY;
            are_cam.push_back(are_alignment);
            alignments.push_back(alignment);
        }
    }


    cout<<CV_VERSION<<endl;
//     EX_TRACE("\nMARKERAUTO --version 1.6.3\n")



    
    vector<float> angles;
    if(!util::ReadAnglesByName(opts.inputangle, &angles)) {
        std::cout<<"Can't open tilt angle file."<<endl;
        return -1;
    }

    int tmp_index=0;
    // for (auto& alignment : are_cam) {
    //     alignment.beta = angles[tmp_index];
    //     tmp_index++;
    //     std::cout << "s: " << alignment.s << ", alpha: " << alignment.alpha
    //               << ", beta: " << alignment.beta << ", gamma: " << alignment.gamma<<", TX: "<<alignment.t0<< ", TY: "<<alignment.t1 << std::endl;
    // }
    
    util::MrcStack mrcs;
    mrcs.Open(opts.input);

    util::FiducialStack fidstack;
//  fidstack.ReadFidsByFile("fids.txt");
#ifndef TESTDEVELOP
    EX_TIME_BEGIN("\n%sDo DetectorMain", _WAVE)
    Detector::DetectorMain(mrcs, &fidstack, opts.diameter, 1);
    fidstack.WriteFidsByFile("fids.txt");

    if(opts.verbose >= 1) {
        Detector::Test(mrcs, fidstack, opts.diameter);
    }

    EX_TIME_END("Do DetectorMain")
#else
    fidstack.ReadFidsByFile("fids.txt");
#endif
    int sum=0;
    for(int i=0;i<fidstack.Size();i++)
    {
        std::vector<util::point2d>& fids = fidstack.V(i);
        sum=sum+fids.size();
    }
    std::cout<<sum<<std::endl;

    util::ImgMatchVector imvector;
    HSetVector hset;
// #undef TESTDEVELOP
#ifndef TESTDEVELOP
    cv::Mat tmplt;
    util::SeriesReadFromFile(&tmplt, "avgtmplt");
//     std::cout<<tmplt.size().width<<std::endl;
    EX_TIME_BEGIN("\n%sDo MatchMain", _WAVE)
    ModelMatch::MatchMain(fidstack, angles, &imvector, &hset, 0.5*tmplt.size().width, true, opts.testmode);	//0.85
    imvector.WriteVectorByFolder("matches");
    hset.WriteVectorByFolder("transmx");


    if(opts.verbose >= 1) {
        ModelMatch::Test(mrcs, imvector, 0.5f, "matches_ill");
    }

    EX_TIME_END("Do MatchMain")
#else
    imvector.ReadVectorByFolder("matches");
    hset.ReadVectorByFolder("transmx");
#endif
    
    util::TrackSpace trackspace;
    trackspace.Create(imvector, angles);
    int num = 0;
        for(int i = 0; i < trackspace.Size(); i++){
        util::TrackSpace::Iterator itr = trackspace.Z_Iterator(i);
                while(!itr.IsNULL()){
            //image_data[i].AddKey(itr.X(), itr.Y());
            //itr.SetExtraAsZVectorIndex(num);
            //num++;
            num++;
            itr++;
        }
        }
    std::cout<<num<<std::endl;

    util::FiducialStack addedfsk;
    util::ImgMatchVector addedimv;
    
#ifndef TESTDEVELOP
    Detector::LocalDetectorMain(mrcs, trackspace, hset, -1, &addedfsk, &addedimv);
    addedfsk.WriteFidsByFile("addfids.txt");
    addedimv.WriteVectorByFolder("addimvec");
#else
    addedfsk.ReadFidsByFile("addfids.txt");
    addedimv.ReadVectorByFolder("addimvec");
#endif
// 	Detector::Test(mrcs, addedfsk, opts.diameter, "addedfsk");
// 	ModelMatch::Test(mrcs, addedimv, 0.5f, "addedimv_ill");

    trackspace.InsertMatchVector(addedimv);

    trackspace.CoordinateTransform(fidstack.Width(), fidstack.Height());

    if(opts.rotation_angle < -0.01 || opts.rotation_angle > 0.01) {
        EX_TRACE("Do pre-rotation of series...\n")
        trackspace.PreRotate(-DEG2RAD(opts.rotation_angle));
    }
    
    std::vector<mx::pproj_params> cameras;
    std::vector<v3_t> points;
    EX_TIME_BEGIN("\n%sDo BundleMain", _WAVE)
    PPBundleApp::BundleMain(trackspace, 0/*1024*//*featsk.Width()*/, 0/*1024*//*featsk.Height()*/, are_cam, &cameras, &points);

    std::vector<mx::pproj_params> cameras_tmp;

    for(const auto& camera : are_cam){
         mx::pproj_params camera_tmp;
         camera_tmp.s = camera.s;
         camera_tmp.alpha = camera.alpha*M_PI/180;
         camera_tmp.beta = camera.beta*M_PI/180;
         camera_tmp.gamma = camera.gamma*M_PI/180;
         camera_tmp.t0 = camera.t0;
         camera_tmp.t1 = camera.t1;
         cameras_tmp.push_back(camera_tmp);
    }


    // PPBundleApp::PrintCamerasAsIMOD(cameras, -DEG2RAD(opts.rotation_angle), 1, opts.outputxf, "xtiltangle.txt", opts.outputangle, "invalid.txt");
    EX_TIME_END("Do BundleMain")
    
    mrcs.Close();
}
