#include "sfm_driver.h"
#include "img_projs.h"
#include <cmath>
#include <iostream>
#include "util/exception.h"
#include "util/matrix.h"
#include "matrix/matrix.h"
#include "util/qsort.h"
#include <cfloat>
#include <matrix/vector.h>
#include "geometry_data.h"
#include <ext/algorithm>
#include <numeric>
#include <vector>
#include <lm.h>
#include "ceres/ceres.h"
#include <fstream>
// #include <ext/hash_map>

#define MAXITER 		800

using namespace mx;

bundle::PSFMDriver::PSFMDriver()
{
    Initialize();
}
bundle::PSFMDriver::~PSFMDriver() {}

void bundle::PSFMDriver::Initialize()
{
    /* call sparse LM routine */
    opts[0]=SBA_INIT_MU;
    opts[1]=SBA_STOP_THRESH;
    opts[2]=SBA_STOP_THRESH;
    opts[3]=SBA_STOP_THRESH;
    //opts[3]=0.05*numprojs; // uncomment to force termination if the average reprojection error drops below 0.05
    opts[4]=0.0;
    //opts[4]=1E-12; // uncomment to force termination if the relative reduction in the RMS reprojection error drops below 1E-05
}

#ifndef MIN_POINTS
#define MIN_POINTS		6
#endif

void bundle::PSFMDriver::Project(const double param[6], const double M[3], double n[2])
{
    calcImgParallelProj(param, M, n);
}

void bundle::PSFMDriver::SetGlobMask(bool s, bool alpha, bool beta, bool gamma, bool t0, bool t1, bool p)
{
	glob.cpmask[0] = s;
	glob.cpmask[1] = alpha;
	glob.cpmask[2] = beta;
	glob.cpmask[3] = gamma;
	glob.cpmask[4] = t0;
	glob.cpmask[5] = t1;
	glob.pmask = p;
}

static bool pairCompare(const std::pair<double, std::pair<int, int> >& firstElem, const std::pair<double, std::pair<int, int> >& secondElem){
  return firstElem.first > secondElem.first;
}

// 定义残差函数
struct Residual {
    Residual(double observed_x, double observed_y) : observed_x(observed_x), observed_y(observed_y) {}

    template<typename T>
    bool operator()(const T* const point, T* residuals) const {
        residuals[0] = point[0] - T(observed_x);
        residuals[1] = point[1] - T(observed_y);
        return true;
    }

private:
    double observed_x;
    double observed_y;
};

// 定义相机模型

struct CameraModel {

    template<typename T>

    bool operator()(const T* const camera, const T* const point, T* residuals) const {
//     bool operator()(const T* const camera, const T* const point, T* residuals) const {
//         // 从相机参数中提取焦距和主点
//         const T& fx = camera[0];
//         const T& fy = camera[1];
//         const T& cx = camera[2];
//         const T& cy = camera[3];
//         // 计算投影坐标
//         T x = point[0];
//         T y = point[1];
//         T z = point[2];
//         T u = fx * x / z + cx;
//         T v = fy * y / z + cy;
//         // 计算残差
//         residuals[0] = u - T(observed_x);
//         residuals[1] = v - T(observed_y);

        T X, Y, Z;
        X=point[0];
        Y=point[1];
        Z=point[2];

//         mx::pproj_params cam;
//         cam.s=camera[0];
//         cam.alpha=camera[1];
//         cam.beta=camera[2];
//         cam.gamma=camera[3];
//         cam.t0=camera[4];
//         cam.t1=camera[5];
        const T& s=camera[0];
        const T& alpha=camera[1];
        const T& beta=camera[2];
        const T& gamma=camera[3];
        const T& t0=camera[4];
        const T& t1=camera[5];

        T point3d[3], pr[2];
        point3d[0]=X;point3d[1]=Y;point3d[2]=Z;
//         bundle::PSFMDriver::Project(cam.mot, point3d, pr);

        T cos_alpha = cos(alpha);
        T sin_alpha = sin(alpha);
        T cos_beta = cos(beta);
        T sin_beta = sin(beta);
        T cos_gamma = cos(gamma);
        T sin_gamma = sin(gamma);

        T tm1, tm2;
        tm1 = (cos_beta*X+sin_alpha*sin_beta*Y-cos_alpha*sin_beta*Z)/s-t0;
        tm2 = (cos_alpha*Y+sin_alpha*Z)/s-t1;
        pr[0] = cos_gamma*tm1-sin_gamma*tm2;
        pr[1] = sin_gamma*tm1+cos_gamma*tm2;

                // 计算残差
        residuals[0] = pr[0] - T(observed_x);
        residuals[1] = pr[1] - T(observed_y);

        return true;

    }

    double observed_x;
    double observed_y;

};

/** the first index of @c pt_views is the index of @c added_order*/
double bundle::PSFMDriver::Run_are(int* added_order, mx::pproj_params* pparams, const int start_camera, const int ncams, const int nconcam,
                              v3_t* init_pts, const int n3Dpts, const int ncon3Dpts, std::vector<ImageKeyVector>& pt_views, GeometryData& data, bool remove_outliers)
{
    int total_removed_points = 0;
    int num_outliers = 0;

    double dist_total = 0.0;
    int num_dists = 0;
	double global_error;

    int *remap = new int [n3Dpts];		//store the index of 3D points in nz_pts, if no, store -1
    double* motstruct = new double[ncams*cnp + n3Dpts*pnp];
    double* nz_pts = motstruct+ncams*cnp;

    char* vmask = new char[n3Dpts*ncams];

    int num_projections = 0;
    for(int i = 0; i < n3Dpts; i++) {
        num_projections += (int)pt_views[i].size();
    }

    double* projections = new double[mnp*num_projections];
	int totalkeys[ncams];
	int outacc[ncams];

    int num_3d_pts;

    int arr_idx = 0;
    int nz_count = 0;
    memset(totalkeys, 0, sizeof(int)*ncams);

    /* Set up the vmask and projections */
    memset(vmask, 0, sizeof(char)*n3Dpts*ncams);

    int fixed_pt_num = ncon3Dpts;

    PosePointParametersBlock states;
    PosePointParametersBlock states2;
//    PosePointParametersBlock init_states;
    states.create(ncams, n3Dpts);
    states2.create(ncams, n3Dpts);


//     delete [] vmask;

//     for(int i = 0; i < data.NumImages(); i++){
//         std::cout<<"["<<(pparams+i)->s<<", "<<(pparams+i)->alpha<<", "<<(pparams+i)->beta
//             <<", "<<(pparams+i)->gamma<<", "<<(pparams+i)->t0<<", "<<(pparams+i)->t1<<"]"<<std::endl;
// 	}


    for(int i = 0; i < n3Dpts; i++) {
        int num_views =(int)pt_views[i].size();

        if(num_views > 0) {
            for(int j = 0; j < num_views; j++) {
                int c = pt_views[i][j].first;
                int v = added_order[c];
                int k = pt_views[i][j].second;
                vmask[nz_count * ncams + c] = 1;
                totalkeys[c]++;										//patch
                projections[2 * arr_idx + 0] = data.GetKey(v,k).m_x;
                projections[2 * arr_idx + 1] = data.GetKey(v,k).m_y;

                arr_idx++;
            }

            remap[i] = nz_count;
            memcpy(&nz_pts[nz_count*3], init_pts[i].p, sizeof(double)*3);
            nz_count++;
        }
        else {
            if(i < ncon3Dpts) {
                fixed_pt_num--;
            }
                remap[i] = -1;
        }
    }

//     double* motparams = motstruct;     //相机参数与三维点存储
//     for(int i = 0; i < ncams; i++) {
//         mx::MotCopyFormPProjParams(motparams, pparams[i+start_camera]);
//         motparams += cnp;
//     }
//

//     for(int i = 0; i < data.NumImages(); i++){
//         std::cout<<"["<<(pparams+i)->s<<", "<<(pparams+i)->alpha<<", "<<(pparams+i)->beta
//             <<", "<<(pparams+i)->gamma<<", "<<(pparams+i)->t0<<", "<<(pparams+i)->t1<<"]"<<std::endl;
// 	}


    for(int i=0;i<ncams;i++)
    {
        double *cameras_true (states.pose(i));
//        double *cameras=A;
        cameras_true[0]=(pparams+i)->s;
        cameras_true[1]=(pparams+i)->alpha;
        cameras_true[2]=(pparams+i)->beta;
        cameras_true[3]=(pparams+i)->gamma;
        cameras_true[4]=(pparams+i)->t0;
        cameras_true[5]=(pparams+i)->t1;
//         std::cout<<cameras_true[0]<<"   "<<cameras_true[1]<<"   "<<cameras_true[2]<<"   "<<cameras_true[3]<<std::endl;

    }


    for(int i=0;i<ncams;i++)
    {
        double *cameras_true (states2.pose(i));
//        double *cameras=A;
        cameras_true[0]=(pparams+i)->s;
        cameras_true[1]=(pparams+i)->alpha;
        cameras_true[2]=(pparams+i)->beta;
        cameras_true[3]=(pparams+i)->gamma;
        cameras_true[4]=(pparams+i)->t0;
        cameras_true[5]=(pparams+i)->t1;
//         std::cout<<cameras_true[0]<<"   "<<cameras_true[1]<<"   "<<cameras_true[2]<<"   "<<cameras_true[3]<<std::endl;

    }

//     for(int i=0;i<ncams;i++){
//         std::cout<<states.pose(i)[0]<<"   "<<states.point(i)[1]<<"   "<<states.point(i)[2]<<"   "<<states.point(i)[3]<<"   "<<states.point(i)[4]<<"   "<<states.point(i)[5]<<std::endl;
//     }

    for(int i=0;i<n3Dpts;i++)
    {
        Eigen::Map<Eigen::Vector3d> true_pt(states.point(i));  //类似于引用地址
        true_pt = Eigen::Vector3d(init_pts[i].p[0],
                                  init_pts[i].p[1],
                                  init_pts[i].p[2]);
    }
// std::cout<<init_pts[1].p[0]<<"  "<<init_pts[1].p[1]<<"  "<<init_pts[1].p[2]<<"  "<<std::endl;
//     std::cout<<states.point(1)[0]<<"   "<<states.point(1)[1]<<"   "<<states.point(1)[2]<<std::endl;
    // for(int i=0;i<n3Dpts;i++){
    //     std::cout<<states.point(i)[0]<<"   "<<states.point(i)[1]<<"   "<<states.point(i)[2]<<std::endl;
    // }

    // 创建Ceres Solver优化问题
    ceres::Problem problem;


    std::vector<std::array<double, 3>> Points(n3Dpts);
    for(int i = 0; i < n3Dpts; i++) {
        Points[i] = {init_pts[i].p[0], init_pts[i].p[1], init_pts[i].p[2]};
    }


    std::vector<std::array<double, 6>> Cameras(ncams);
    for(int i = 0; i < ncams; i++) {
//         Cameras[i] = {pparams[i].s, pparams[i].alpha, pparams[i].beta, pparams[i].gamma, pparams[i].t0, pparams[i].t1};
        Cameras[i] = {1,0,pparams[i].beta,0,0,0};
    }

    for (auto& point : Points) {

        problem.AddParameterBlock(point.data(), 3);

    }

   
    int gap=0;
    double x_o,y_o;
//     double point_tmp[3] = {0.0, 0.0, 0.0};  // 初始三维点坐标
    for(int i=0; i<n3Dpts; i++){
        double point[3] = {init_pts[i].p[0], init_pts[i].p[1], init_pts[i].p[2]};  // 初始三维点坐标
        for(int j=0;j<ncams; j++){
            if(vmask[i*ncams+j]){

                x_o=projections[i*(ncams*2)+2*j+0-gap*2];
                y_o=projections[i*(ncams*2)+2*j+1-gap*2];
                problem.AddParameterBlock(point, 3);

                ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<CameraModel, 2, 6, 3>(
                    new CameraModel {x_o, y_o}
                );
//                 problem.AddResidualBlock(cost_function, nullptr, states.pose(j), states.point(i));
                problem.AddResidualBlock(cost_function, nullptr, Cameras[j].data(), Points[i].data());
//                 problem.SetParameterBlockConstant(states.pose(j));

            }
            else{
                gap++;
            }
        }
    }

    // 设置相机参数为常量（不优化）

    for (auto& camera : Cameras) {

        problem.SetParameterBlockConstant(camera.data());

    }

    // 设置优化选项
    ceres::Solver::Options options;
    options.max_num_iterations = 100;
    options.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY;

    // 运行优化
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

//     输出优化结果
//     std::cout << summary.FullReport() << "\n";
    std::cout << summary.BriefReport() << "\n";
//     std::cout << "Optimized 3D point: " << point[0] << ", " << point[1] << ", " << point[2] << "\n";
    // std::cout << "********************优化后三维*************************" << "\n";
//     for(int i=0;i<n3Dpts;i++){
//         std::cout<<states.point(i)[0]<<"   "<<states.point(i)[1]<<"   "<<states.point(i)[2]<<std::endl;
//     }
//
//         for(int i=0;i<ncams;i++){
//         std::cout<<states.pose(i)[0]<<"   "<<states.pose(i)[1]<<"   "<<states.pose(i)[2]<<"   "<<states.pose(i)[3]<<"   "<<states.pose(i)[4]<<"   "<<states.pose(i)[5]<<std::endl;
//     }

    // for (int i = 0; i < n3Dpts; ++i) {

    //     std::cout << "Optimized 3D point " << i+1 << ": "

    //               << Points[i][0] << ", " << Points[i][1] << ", " << Points[i][2] << "\n";

    // }



    gap = 0;
    std::vector<double> f, u;
    std::vector<util::point2d> img_2d;
    std::vector<double> residual_e;
    int arr_num=0;
    for(int i=0;i<n3Dpts;i++)
    {
        double res = 0;
        double X[3], pr[2];
        util::point2d p;
//         Eigen::Vector3d point_tmp=Eigen::Vector3d(init_pts[i].p[0],init_pts[i].p[1],init_pts[i].p[2]);
//         X[0]=init_pts[i].p[0]; X[1]=init_pts[i].p[1]; X[2]=init_pts[i].p[2];
        Eigen::Vector3d point_tmp=Eigen::Vector3d(Points[i][0],Points[i][1],Points[i][2]);
        X[0]=point_tmp[0]; X[1]=point_tmp[1]; X[2]=point_tmp[2];
//        std::cout<<"Point"<<i<<std::endl;
        for(int j=0;j<ncams;j++)
        {
            double x_o,y_o;
            mx::pproj_params cam_tmp;
//             cam_tmp.s= states.pose(j)[0];
//             cam_tmp.alpha= states.pose(j)[1];
//             cam_tmp.beta= states.pose(j)[2];
//             cam_tmp.gamma= states.pose(j)[3];
//             cam_tmp.t0= states.pose(j)[4];
//             cam_tmp.t1= states.pose(j)[5];

            cam_tmp.s= Cameras[j][0];
            cam_tmp.alpha= Cameras[j][1];
            cam_tmp.beta= Cameras[j][2];
            cam_tmp.gamma= Cameras[j][3];
            cam_tmp.t0= Cameras[j][4];
            cam_tmp.t1= Cameras[j][5];
            Project(cam_tmp.mot, X, pr);
//             p = cam_tmp.Project(point_tmp);
//            std::cout<<"cam_tmp.gamma:"<<cam_tmp.gamma<<" cam_tmp.t0:"<<cam_tmp.t0<<std::endl;
            if(vmask[i*ncams+j]){
                x_o=projections[i*(ncams*2)+2*j+0-gap*2];
                y_o=projections[i*(ncams*2)+2*j+1-gap*2];

//                 util::point2d point_tmp;
//                 point_tmp.x=x_o;
//                 point_tmp.y=y_o;
//                 img_2d.push_back(point_tmp);

//                std::cout<<"x_o:"<<x_o<<" y_o:"<<y_o<<std::endl;
//                std::cout<<"pr[0]:"<<pr[0]<<" pr[1]:"<<pr[1]<<std::endl;
//                if(y_o == pr[1])
//                {
//                    std::cout<<"x_o:"<<x_o<<" y_o:"<<y_o<<std::endl;
//                std::cout<<"pr[0]:"<<pr[0]<<" pr[1]:"<<pr[1]<<std::endl;
//             }
                u.push_back(x_o);
                u.push_back(y_o);

                f.push_back(pr[0]);
                f.push_back(pr[1]);

                double dx, dy;
                dx=pr[0]-x_o;
                dy=pr[1]-y_o;
                res = sqrt(dx * dx + dy * dy);
                // std::cout<<"dx:"<<dx<<" dy:"<<dy<<std::endl;
                residual_e.push_back(sqrt(dx * dx + dy * dy));
                arr_num++;
            }
            else{
                f.push_back(0.0);
                f.push_back(0.0);
                u.push_back(0.0);
                u.push_back(0.0);
                gap++;
            }
        }
    }


    double residual = 0;
    //     std::cout<<f.size()<<std::endl;
    for (int i = 0; i < f.size(); i=i+1) {
        double delt;
        delt=sqrt(f[i]*f[i]);
        residual+=delt;
    }
    double avg_res=0, std_sum=0, std=0;

//     sort(residual_e.begin(), residual_e.end());

    double residual_tmp=0;
    int arr_num_tmp=1500;
    for (int i = 0; i < residual_e.size(); i=i+1) {
//     for (int i = 0; i < arr_num_tmp; i=i+1) {
        double delt;
        delt=residual_e[i];
//         std::cout<<"delt:"<<delt<<std::endl;
        residual_tmp=residual_tmp+delt;
    }

    avg_res=residual_tmp/arr_num;
    for (int i = 0; i < residual_e.size(); i=i+1) {
//     for (int i = 0; i < arr_num_tmp; i=i+1) {
        double delt=0;
        delt=residual_e[i]-avg_res;
        std_sum=delt*delt+std_sum;
    }

    std = sqrt(std_sum/arr_num);
// std = sqrt(std_sum/32000);

    std::cout<<std::endl;
//     std::cout<<"ResErr:"<<residual/(arr_num)<<std::endl;
    std::cout<<"ResErr:"<<residual_tmp/(arr_num)<<"   std:"<<std<<std::endl;
//     std::cout<<"ResErr_new:"<<residual_tmp/(32000)<<"   std:"<<std<<std::endl;
    // std::cout<<"residual:"<<residual<<"        arr_num:"<<arr_num<<std::endl;



//     std::cout<<"ResErr_new:"<<residual_tmp/(arr_num_tmp)<<"   std:"<<std<<std::endl;
// //     std::cout<<"ResErr_new:"<<residual_tmp/(32000)<<"   std:"<<std<<std::endl;
//     std::cout<<"residual:"<<residual<<"        arr_num:"<<arr_num_tmp<<std::endl;
















//     double init_residual=0;
//     double max=0;
//     std::vector<double> delt_2000;
//     std::cout<<"u.size():"<<u.size()<<std::endl;
//     for (int i = 0; i < u.size(); i=i+1) {
//         double delt;
// //         if(i == 2000) break;
//         delt=sqrt((u[i]-f[i])*(u[i]-f[i]));
// //         if(delt!=0){
//         if(1){
//             delt_2000.push_back(delt);
//         }
// //         else{
// //             std::cout<<"      0:"<<u[i]<<"     "<<f[i]<<std::endl;
// //         }
// //         if(delt_2000.size()==1000) break;
// //         delt_2000.push_back(delt);
// //         std::cout<<"(u[i]-f[i]):"<<u[i]<<"     "<<f[i]<<std::endl;
// //         std::cout<<"(u[i]-f[i]):"<<delt<<std::endl;
// //         if(delt>max) max=delt;
// //        delt=abs(u[i]-f[i]);
// //         init_residual=init_residual+delt;
//     }
//
//     sort(delt_2000.begin(), delt_2000.end());
//
//     std::cout<<"u.size():"<<delt_2000.size()<<std::endl;
//     int arr_vmask=0;
//
//     for(int i=0;i<n3Dpts*ncams;i++)
//     {
//         if(vmask[i]) arr_vmask++;
//     }
//
// //     int arr_num=arr_vmask;
//     int arr_num=delt_2000.size();
//     std::cout<<"arr_num:"<<arr_num<<std::endl;
//
//     for (int i = 0; i < delt_2000.size(); i=i+1) {
// //         if(i==arr_num) break;
// //         std::cout<<"delt_2000[i]:"<<delt_2000[i]<<std::endl;
//         init_residual=init_residual+delt_2000[i];
//     }
//
// //     init_residual = accumulate(delt_2000.begin(), delt_2000.end(), 0);
// //     std::cout<<"max:"<<max<<std::endl;
//
// //     int arr_num=0, arr_vmask=0;
//     std::cout<<"arr_vmask:"<<arr_vmask<<std::endl;
//
// //     std::cout<<"init_ResErr:"<<init_residual/(arr_vmask)<<std::endl;
//     std::cout<<"init_ResErr:"<<init_residual/arr_num<<std::endl;

}

/** the first index of @c pt_views is the index of @c added_order*/
double bundle::PSFMDriver::Run(int* added_order, mx::pproj_params* pparams, const int start_camera, const int ncams, const int nconcam,
                              v3_t* init_pts, const int n3Dpts, const int ncon3Dpts, std::vector<ImageKeyVector>& pt_views, GeometryData& data, bool remove_outliers)
{
    int total_removed_points = 0;
    int num_outliers = 0;

    double dist_total = 0.0;
    int num_dists = 0;
	double global_error;

    int *remap = new int [n3Dpts];		//store the index of 3D points in nz_pts, if no, store -1
    double* motstruct = new double[ncams*cnp + n3Dpts*pnp];
    double* nz_pts = motstruct+ncams*cnp;

    char* vmask = new char[n3Dpts*ncams];

    int num_projections = 0;
    for(int i = 0; i < n3Dpts; i++) {
        num_projections += (int)pt_views[i].size();
    }

    double* projections = new double[mnp*num_projections];
	int totalkeys[ncams];
	int outacc[ncams];

    int num_3d_pts;

//             std::string res_file1 = "/home/xzh/桌面/W.txt";
//     std::ofstream outputfile_res1(res_file1);
//
//     for(int i=0;i<ncams*cnp + n3Dpts*pnp;i++)
//     {
//         outputfile_res1<<motparams[i]<<std::endl;
// //         outputfile_res<<"\n";
//     }
//
//     outputfile_res1.close();


    do {
        if((num_3d_pts = n3Dpts - total_removed_points) < MIN_POINTS) {
            EX_PRINT("# Too few points remaining, exiting!\n")

            dist_total = DBL_MAX;
            break;
        }

        int arr_idx = 0;
        int nz_count = 0;
		memset(totalkeys, 0, sizeof(int)*ncams);

        /* Set up the vmask and projections */
        memset(vmask, 0, sizeof(char)*num_3d_pts*ncams);

        int fixed_pt_num = ncon3Dpts;

        for(int i = 0; i < n3Dpts; i++) {
            int num_views =(int)pt_views[i].size();

            if(num_views > 0) {
                for(int j = 0; j < num_views; j++) {
                    int c = pt_views[i][j].first;
                    int v = added_order[c];
                    int k = pt_views[i][j].second;
                    vmask[nz_count * ncams + c] = 1;
					totalkeys[c]++;										//patch
                    projections[2 * arr_idx + 0] = data.GetKey(v,k).m_x;
                    projections[2 * arr_idx + 1] = data.GetKey(v,k).m_y;

                    arr_idx++;
                }

                remap[i] = nz_count;
                memcpy(&nz_pts[nz_count*3], init_pts[i].p, sizeof(double)*3);
                nz_count++;
            }
            else {
                if(i < ncon3Dpts) {
                    fixed_pt_num--;
                }
                remap[i] = -1;
            }
        }

        double* motparams = motstruct;
        for(int i = 0; i < ncams; i++) {
            mx::MotCopyFormPProjParams(motparams, pparams[i+start_camera]);
            motparams += cnp;
        }

        dist_total = 0.0;
        num_dists = 0;

        // ###


        double *Point3D=new double[n3Dpts*3];

        for(int i=0;i<n3Dpts;i++)
        {
            Point3D[3*i]=motstruct[ncams*6+3*i];
            Point3D[3*i+1]=motstruct[ncams*6+3*i+1];
            Point3D[3*i+2]=motstruct[ncams*6+3*i+2];
        }

        mx::pproj_params* camparams = new mx::pproj_params[ncams];


        for(int i=0;i<ncams;i++)
        {
            mx::pproj_params cam_tmp;
            cam_tmp.s=motstruct[i*6]; cam_tmp.alpha=motstruct[i*6+1]; cam_tmp.beta=motstruct[i*6+2];
            cam_tmp.gamma=motstruct[i*6+3]; cam_tmp.t0=motstruct[i*6+4]; cam_tmp.t1=motstruct[i*6+5];
            camparams[i]=cam_tmp;
        }

        int  gap = 0;
        std::vector<double> f, u;
        std::vector<util::point2d> img_2d;
        for(int i=0;i<n3Dpts;i++)
        {
            double X[3], pr[2];
            util::point2d p;
    //         Eigen::Vector3d point_tmp=Eigen::Vector3d(init_pts[i].p[0],init_pts[i].p[1],init_pts[i].p[2]);
    //         X[0]=init_pts[i].p[0]; X[1]=init_pts[i].p[1]; X[2]=init_pts[i].p[2];
            Eigen::Vector3d point_tmp=Eigen::Vector3d(Point3D[3*i],Point3D[3*i+1],Point3D[3*i+2]);
            X[0]=Point3D[3*i]; X[1]=Point3D[3*i+1]; X[2]=Point3D[3*i+2];
    //        std::cout<<"Point"<<i<<std::endl;
            for(int j=0;j<ncams;j++)
            {
                double x_o,y_o;
//                 mx::pproj_params cam_tmp;
                mx::pproj_params cam_tmp = camparams[j];
                Project(cam_tmp.mot, X, pr);
    //             p = cam_tmp.Project(point_tmp);
    //            std::cout<<"x:"<<p.x<<" y:"<<p.y<<std::endl;
                if(vmask[i*ncams+j]){
                    x_o=projections[i*(ncams*2)+2*j+0-gap*2];
                    y_o=projections[i*(ncams*2)+2*j+1-gap*2];

    //                std::cout<<"x_o:"<<x_o<<" y_o:"<<y_o<<std::endl;
                    u.push_back(x_o);
                    u.push_back(y_o);

                    f.push_back(pr[0]);
                    f.push_back(pr[1]);
                }
                else{
                    f.push_back(0.0);
                    f.push_back(0.0);
                    u.push_back(0.0);
                    u.push_back(0.0);
                    gap++;
                }
            }
        }

        double init_residual=0;
        double max=0;
        for (int i = 0; i < u.size(); i=i+1) {
            double delt;
            delt=sqrt((u[i]-f[i])*(u[i]-f[i]));
            init_residual=init_residual+delt;
        }

        int arr_num=0, arr_vmask=0;
        for(int i=0;i<n3Dpts*ncams;i++)
        {
            if(vmask[i]) arr_vmask++;
        }
        std::cout<<"arr_vmask:"<<arr_vmask<<std::endl;

        std::cout<<"init_ResErr:"<<init_residual/(arr_vmask)<<std::endl;





        //###










        EX_BEGIN_CLOCK()
        Run(motstruct, projections, vmask, num_3d_pts, ncams, nconcam, fixed_pt_num);
        EX_END_CLOCK()
        EX_TRACE("# SFM using %d 3D pts(%d fixed), %d frames(%d fixed) and %d image projections(%g p/p), error %g [initial %g](elapse: %ld)\n",
                 num_3d_pts, ncon3Dpts, ncams, nconcam, arr_idx, ((double)arr_idx)/num_3d_pts, sqrt(info[1]/arr_idx)>1?sqrt(info[1]/arr_idx):info[1]/arr_idx,
				 sqrt(info[0]/arr_idx)>1?sqrt(info[0]/arr_idx):info[0]/arr_idx, EX_ELAPSE());

        double avg_res=0, std_sum=0, std=0;
        std::vector<double> res;

		global_error = sqrt(info[1]/arr_idx);

        motparams = motstruct;
        for(int i = 0; i < ncams; i++){
            mx::PProjPCopyFormMotParams(&pparams[i+start_camera], motparams);
            motparams += cnp;
        }

        std::vector<std::pair<int, int> > outliers;
		std::vector<std::pair<double, std::pair<int, int> > > outidx;

        if(!remove_outliers){
            goto end;
        }

        for(int i = 0; i < ncams; i++){
            double params[6];
			mx::MotCopyFormPProjParams(params, pparams[i+start_camera]);

            int num_keys = data.GetNumKeys(added_order[i]);

            int num_pts_proj = 0;
            for(int j = 0; j < num_keys; j++) {
                if(data.GetKey(added_order[i], j).m_extra >= 0) {
                    num_pts_proj++;
                }
            }

            double *dists = new double[num_pts_proj];
            int pt_count = 0;

            std::vector<Keypoint>::iterator iter;

            for(iter = data.m_image_data[added_order[i]].m_keys.begin(); iter != data.m_image_data[added_order[i]].m_keys.end(); iter++) {
                const Keypoint &key = *iter;

                if(key.m_extra >= 0) {
                    int pt_idx = key.m_extra;
                    double X[3], pr[2];
                    memcpy(X, &nz_pts[remap[pt_idx]*3], sizeof(double)*3);

                    Project(params, X, pr);

                    double dx = pr[0]-key.m_x;
                    double dy = pr[1]-key.m_y;

                    double dist = sqrt(dx * dx + dy * dy);
                    dist_total += dist;
                    num_dists++;

                    dists[pt_count] = dist;

                    pt_count++;

                    res.push_back(dist);
                }
            }

            /* Estimate the median of the distances */
            double med = kth_element_copy(num_pts_proj, int(0.5/*0.7/*0.8/* 0.9 */* num_pts_proj), dists);

#define NUM_STDDEV 2//2.0//3.0//6.0
            double thresh = 1.2 * NUM_STDDEV * med;/* k * stddev */
            thresh = CLAMP(thresh, min_proj_error_threshold, max_proj_error_threshold);
			if(global_error > max_proj_error_threshold){		//considering that the global error; large noise
				thresh = global_error;
			}

            /* Compute the average reprojection error for this camera */
            double sum = 0.0;
            for(int j = 0; j < num_pts_proj; j++){
                sum += dists[j];
            }

            double avg = sum/num_pts_proj;
            EX_PRINT("# Camera %d[%d] (%d pts), mean error: %0.3f [median %0.3f(0.7 quantile %0.3f), error thresh %0.3f]\n",
					 i, added_order[i], num_pts_proj, avg, med, kth_element_copy(num_pts_proj, int(0.7 * num_pts_proj), dists), thresh);
                    // i, added_order[i], num_pts_proj, avg, kth_element_copy(num_pts_proj, int(0.5 * num_pts_proj), dists), med, thresh);

            pt_count = 0;
			outidx.clear();
            for(int j = 0; j < num_keys; j++) {
                int pt_idx = data.GetKey(added_order[i],j).m_extra;

                if(pt_idx < 0)
                    continue;

                if(dists[pt_count] > thresh){ //|| dists[pt_count] > max_proj_error_threshold) {
//                     if(dists[pt_count] > 10000){ //|| dists[pt_count] > max_proj_error_threshold) {
                    /* Remove this point from consideration */
//                     outliers.push_back(std::pair<int, int>(pt_idx, i));
					outidx.push_back(std::make_pair(dists[pt_count], std::pair<int, int>(pt_idx, i)));		//? sort of outliers??? WARNING
                }
                pt_count++;
            }

			std::sort(outidx.begin(), outidx.end(), pairCompare);

			for(int i = 0; i < outidx.size(); i++){
				outliers.push_back(outidx[i].second);
			}

            delete [] dists;
        }
        avg_res=dist_total/num_dists;
            std::cout<<avg_res<<"      "<<num_dists<<std::endl;
            for(int i=0;i<res.size();i++)
            {
                double delt=0;
                delt=res[i]-avg_res;
                std_sum=delt*delt+std_sum;
            }
            std = sqrt(std_sum/num_dists);
            std::cout<<"ResErr_new:"<<avg_res<<"   std:"<<std<<std::endl;

        /* Remove outlying points */
        num_outliers = 0;
		memset(outacc, 0, sizeof(int)*ncams);
        for(int i = 0; i <(int)outliers.size(); i++) {
            int idx = outliers[i].first;

            if(idx < ncon3Dpts) {
                continue;
            }

            if(!pt_views[idx].size()) {
                continue;
            }

            for(ImageKeyVector::iterator itr = pt_views[idx].begin(); itr != pt_views[idx].end(); ) {
                int v = (*itr).first;
                int k = (*itr).second;
                if(v == outliers[i].second && totalkeys[added_order[v]]-outacc[added_order[v]] > 6){
                    if(data.GetKey(added_order[v], k).m_extra != idx) {
                        EX_ERROR("Error!  Entry for(%d,%d) should be %d, but is %d\n",
                                 added_order[v], k, idx, data.GetKey(added_order[v], k).m_extra);
                    }
                    data.GetKey(added_order[v], k).m_extra = -2;
                    pt_views[idx].erase(itr);
                    num_outliers++;
					outacc[added_order[v]]++;
                    break;
                }
                else {
                    itr++;
                }
            }

            for(int i = 0; i < ncams; i++){
				totalkeys[i] -= outacc[i];
            }

            if(pt_views[idx].size() < 2) {
                for(ImageKeyVector::iterator itr = pt_views[idx].begin(); itr != pt_views[idx].end(); itr++) {
                    int v = (*itr).first;
                    int k = (*itr).second;

                    data.GetKey(added_order[v], k).m_extra = -2;
                }

                pt_views[idx].clear();
                total_removed_points++;
            }
        }

        outliers.clear();

        EX_PRINT("# Removing %d outliers\n", num_outliers);

end:
        for(int i = 0; i < n3Dpts; i++) {
            if(remap[i] != -1) {
                memcpy(init_pts[i].p, &nz_pts[remap[i]*3], sizeof(double)*3);
            }
        }

    } while(num_outliers > 0);

    delete [] vmask;
    delete [] projections;

    delete [] remap;
    delete [] motstruct;

    return dist_total/num_dists;
}

v2_t bundle::PSFMDriver::Project(const mx::pproj_params& pparam, v3_t pt)
{
	v2_t p2;
    Project(pparam.mot, pt.p, p2.p);

    return p2;
}

void bundle::PSFMDriver::RefineCameraParameters(v3_t* points, v2_t* projs, int num_points, mx::pproj_params* pparams)			//
{
    double* motstruct = new double[cnp + num_points*pnp];
    double* nz_pts = motstruct+cnp;

    char* vmask = new char[num_points];
    double* projections = new double[mnp*num_points];

    memset(vmask, 0, sizeof(char)*num_points);

    double* proj = projections;
    double* z_3dpt = nz_pts;
    for(int i = 0; i < num_points; i++) {
        memcpy(proj, projs[i].p, sizeof(double)*2);
        proj += 2;
        memcpy(z_3dpt, points[i].p, sizeof(double)*3);
        z_3dpt +=3;
        vmask[i] = 1;
    }

    mx::MotCopyFormPProjParams(motstruct, *pparams);

    EX_BEGIN_CLOCK()
    Run(motstruct, projections, vmask, num_points, 1, 0, num_points);
    EX_END_CLOCK()
    EX_TRACE("# CameraRefine using %d pts, %d frames(%d fixed) and %d image projections, error %g [initial %g](elapse: %ld)\n",
             num_points, 1, 0, num_points, info[1]/num_points, info[0]/num_points, EX_ELAPSE());

    delete [] motstruct;
    delete [] vmask;
    delete [] projections;

	mx::PProjPCopyFormMotParams(pparams, motstruct);

    return;
}

int bundle::PSFMDriver::Run(double* motstruct, double* imgpts, char* vmask,
                           const int n3Dpts, const int ncams, const int nconcam, const int ncon3Dpts)
{
    int nvars = ncams*cnp+n3Dpts*pnp;
    int n = sba_motstr_levmar_x(n3Dpts, ncon3Dpts, ncams, nconcam, vmask, motstruct,
                                cnp, pnp, imgpts, NULL, mnp, img_ParallelProj_x,
                                img_ParallelProj_jac_x, &glob, MAXITER, verbose, opts, info);


    return n;
}

