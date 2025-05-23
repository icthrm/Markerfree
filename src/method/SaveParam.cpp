#include "Util.h"
#include "../mrc/mrcstack.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

SaveParam::SaveParam()
{
}

SaveParam::~SaveParam()
{
}

void SaveParam::GetParam(int bin)
{
	mbin = bin;
    char filename[256] = {'\0'};
    strcpy(filename, outfile);
    char* pcOutSlash = strrchr(filename, '/');
	if(pcOutSlash == 0L) printf("Error output path");
    char* pcMrc = strstr(filename, ".mrc");
    if(pcMrc == 0L) 
        strcat(filename, ".txt");
    else 
        strcpy(pcMrc, ".txt");
    FILE* pFile = fopen(filename, "wt");
    if(pFile == 0L)
	{	printf("Unable to open %s.\n", filename);
		printf("Alignment data will not be saved\n\n");
		return;
	}
    m_pvFile = pFile;
    SaveHeader();
    SaveAllParam();
    SaveTime();
    CloseFile();
}

void SaveParam::SaveHeader()
{
	FILE* pFile = (FILE*)m_pvFile;
	fprintf(pFile, "# AreTomo Alignment\n");
	fprintf(pFile, "# RawSize = %d %d %d\n", output->X(),output->Y(),output->Z());
	//---------------------------------------------------
}

void SaveParam::SaveAllParam()
{
    if(param == 0L) return;
	FILE* pFile = (FILE*)m_pvFile;
	fprintf( pFile, "# SEC       ROT          "
	   "TX        TY         TILT\n");
	//--------------------------------------------------------------------
	float afShift[] = {0.0f, 0.0f};
	for(int i=0; i<output->Z(); i++)
	{
		float fTilt = angles[i];
		float fTiltAxis = param->rotate[i];
		afShift[0] = param->shiftX[i] * mbin;
        afShift[1] = param->shiftY[i] * mbin;
		fprintf( pFile, "%5d   %9.4f  %9.3f  %9.3f   "
		   "%8.2f\n", i, fTiltAxis, 
		   afShift[0], afShift[1], fTilt);
	}
    fprintf(pFile, "angleoffset: %f\n", param->angleOffset);
	fprintf(pFile, "Output Bin: %d\n", iBin);
}

void SaveParam::SaveTime()
{
    FILE* pFile = (FILE*)m_pvFile;
    fprintf(pFile, "Elapsed time: %f seconds\n", time);
}

void SaveParam::CloseFile()
{
	if(m_pvFile == 0L) return;
	fclose((FILE*)m_pvFile);
	m_pvFile = 0L;
}
