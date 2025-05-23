#pragma once
#include <stdio.h>
#include <math.h>
#include <vector>
#include <cufft.h>
#include "Cufft1&2.h"
#include "../mrc/mrcstack.h"

enum MAXIN {
  MIN, //代表最小
  MAX //代表最大
};

struct AlignParam  //如果要加入每张图像的倾斜角数据angles，想办法把它从向量变成数组
{
    float* shiftX;
    float* shiftY;
    float* rotate;
    float angleOffset;
};

struct SysInfo {
  int id;
  int procs;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int namelen;
};

int GetFrameIdxFromTilt(int nz, std::vector<float> p_angles, float fTilt);
int FindRefIndex(int nz, std::vector<float> p_angles, int z); //找到相邻两张图像中倾角绝对值比本身小的那个,也就是更靠近零倾斜的那个
const double D2R = 4.0 * atan(1.0) / 180.0;

class PreProcess
{
  public:
    PreProcess();
    ~PreProcess();
    void SetPositive();  //负责像素点变正值
    void MassNormalization();  //负责图像归一化
  public:
    MrcStackM* rawstack;
    MrcStackM* stack;
    std::vector<float> angles;
  private:
    float FindMinAndMax(float *proj, int width, int height, MAXIN maxin);
    float mCalcMean(int i);
    void mScale(int i, float proj);
  private:
    int mstart[2];
    int msize[2];
    int nz;
};

class CalcCC
{
  public:
    CalcCC();
    ~CalcCC();
    void Setup(int* BinSize, int Factor);
    void SetSize(int* BinSize);
    void DoIt(float* pfRefImg, float* pfImg, float fRefTilt, float fTilt, float fTiltAxis);
    float CalcMoment(float* proj, int width, int height, int padwidth, int Exponent);
    void Norm2D(float* proj, int Padwidth, int Padheight, float mean, float std);
    void mNormalize(float* proj);
    void Stretch(float* inproj, int* piSize, bool bPadded, float dStretch, float TiltAxis, float *outproj, bool Randfill);
    void RoundEdge(float* Proj, int* piSize, bool bPadded, float fPower, float* MaskCent, float* MaskSize);
    void preforCC(cufftComplex* RefProj, cufftComplex* fProj, float Factor);
    void getshift(float &shiftX, float &shiftY, float binX, float binY);
    // void ProjCC(cufftComplex* RefProj, cufftComplex* fProj, float Factor);
  private:
    void PadProj(float* proj, float* padproj);
    void getCC(cufftComplex* RefProj, cufftComplex* fProj, float Factor, float* m_pfXcfImg);
    float FindPeak(float* CC);
  private:
    int BinSize[2];
    int PadSize[2];
    int CmpSize[2]; 
    int m_Factor;
    cufft2D m_fft;
    cufft2D in_fft;
    cufftComplex* PadRefProj;
    cufftComplex* PadfProj;
    cufftComplex* stretchProj;
    float fshiftX;
    float fshiftY;
};

class SaveParam
{
  public:
    SaveParam();
    ~SaveParam();
    void GetParam(int bin);
  public:
    MrcStackM* output;
    AlignParam* param;
    std::vector<float> angles;
    int iBin;
    float time;
    const char* outfile;
  private:
    void SaveHeader();
    void SaveAllParam();
    void SaveTime();
    void CloseFile();
  private:
    void* m_pvFile;
    int mbin;
};