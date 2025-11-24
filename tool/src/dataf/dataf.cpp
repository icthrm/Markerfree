#include "dataf.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <iostream>
#include "keypoint.h"
#include <sys/stat.h>
#include <unistd.h>
#include <iomanip>


void util::FiducialStack::WriteFidsByFile(const char* filename) const {
    std::ofstream out(filename);
    if(!out.good()) {
        ex::EX_THROW("Can't Create File");
    }
    std::cout <<std::setprecision(8)<<std::endl;
    out<<width<<"\t"<<height<<"\t"<<size<<"\t"<<ratio<<std::endl;
    for(int i = 0 ; i < size ; i++) {
        out<<"frame "<<i<<std::endl;
        out<<"num = "<<vfidp[i].size() <<std::endl;
        for(int j = 0 ; j < vfidp[i].size() ; j++) {
            float x = vfidp[i][j].x;
            float y = vfidp[i][j].y;
            out<<x<<"\t"<<y<<std::endl;
        }
    }
    return;
}

bool util::FiducialStack::ReadFidsByFile(const char* filename) {
    std::string s;
    std::stringstream ss;
    std::string first_str;
    char ch;

    std::ifstream fin(filename);

    if(!fin.good()) {
        std::cout<<"Unable to open "<<filename<<std::endl;
        return false;
    }

    getline(fin , s);

    ss.str("");
    ss<<s;
    ss>>width>>height>>size>>ratio;

    if(vfidp) {
        delete [] vfidp;
    }
    vfidp = new std::vector<util::point2d>[size];

    while(getline(fin ,s)) {
        ss.clear();
        ss<<s;
        ss>>first_str;
        int feats_num;
        if(first_str == "frame") {
            int frame;
            ss>>frame;
            getline(fin ,s);
            ss.clear();
            ss<<s;
            ss>>first_str;

            if(first_str == "num") {
                ss>>ch>>feats_num;
                for(int i=0; i< feats_num ; i++) {
                    util::point2d pt;
                    fin>>pt.x>>pt.y;
                    vfidp[frame].push_back(pt);
                }
            }
        }
    }
    return true;
}

void util::ImgMatchVector::PrintPairs(int index, std::ostream& o) const
{
    o<<(*match_vector)[index].idx1<<" "<<(*match_vector)[index].idx2<<std::endl;
    for(int i = 0; i < (*match_vector)[index].pairs.size(); i++) {
        o<<"("<<((*match_vector)[index].pairs)[i].first.x<<","<<((*match_vector)[index].pairs)[i].first.y<<")&"
         <<"("<<((*match_vector)[index].pairs)[i].second.x<<","<<((*match_vector)[index].pairs)[i].second.y<<")\n";
    }
    o<<"#";
}

void util::ImgMatchVector::CoordinateTransform(int width, int height)
{
    for(int index = 0; index < (*match_vector).size(); index++) {
        for(int i = 0; i < (*match_vector)[index].pairs.size(); i++) {
            ((*match_vector)[index].pairs)[i].first.x -= width/2;
            ((*match_vector)[index].pairs)[i].first.y -= height/2;
            ((*match_vector)[index].pairs)[i].second.x -= width/2;
            ((*match_vector)[index].pairs)[i].second.y -= height/2;
        }
    }
}

void util::ImgMatchVector::PreRotate(float angle)
{
    for(int index = 0; index < (*match_vector).size(); index++) {
        for(int i = 0; i < (*match_vector)[index].pairs.size(); i++) {
            ((*match_vector)[index].pairs)[i].first.x = cos(angle)*((*match_vector)[index].pairs)[i].first.x-sin(angle)*((*match_vector)[index].pairs)[i].first.y;
            ((*match_vector)[index].pairs)[i].first.y = sin(angle)*((*match_vector)[index].pairs)[i].first.x+cos(angle)*((*match_vector)[index].pairs)[i].first.y;
            ((*match_vector)[index].pairs)[i].second.x = cos(angle)*((*match_vector)[index].pairs)[i].second.x-sin(angle)*((*match_vector)[index].pairs)[i].second.y;
            ((*match_vector)[index].pairs)[i].second.y = sin(angle)*((*match_vector)[index].pairs)[i].second.x+cos(angle)*((*match_vector)[index].pairs)[i].second.y;
        }
    }
}

void util::ImgMatchVector::WriteVectorByFolder(const char* folderpath) const
{
    if(access(folderpath,0) == -1) {		//create file folder
        mkdir(folderpath,0777);
    }
    std::ostringstream ooo;
    ooo <<folderpath<<"/attributes";
    std::ofstream out(ooo.str().c_str());
    out<<"Z:"<<match_vector->size()<<"\n";
    out.close();
    for(int i = 0; i < match_vector->size(); i++) {
        std::ostringstream oss;
        oss <<folderpath<<"/"<<i<<".txt";
        std::ofstream o(oss.str().c_str());
        try {
            PrintPairs(i, o);
        } catch(ex::Exception& e) {
            EX_TRACE("%s\n", e.Msg())
        }
    }
}

util::NodePlane::_node* util::NodePlane::AddNode(const util::_point& point)
{
	_node newnode(point);
	newnode.ref = information;
	bool found = false;
	
	std::list<_node>::iterator itr;
	
	if(pre_visited != node_list.begin() && *pre_visited > newnode){
		itr = node_list.begin();
	}
	else{
		itr = pre_visited;
	}
	
	for(; itr != node_list.end(); itr++){
		if(*itr == newnode){
			pre_visited = itr;
			found = true;
			break;
		}
		if(*itr > newnode){
			break;
		}
	}
	
	if(!found){
		pre_visited = node_list.insert(itr, newnode);
	}
	return &(*pre_visited);
}

void util::TrackSpace::CoordinateTransform(int width, int height)
{
	float w_2 = width*.5, h_2 = height*.5;
    for(int i = 0; i < size; i++) {
        NodePlane::NodeReader reader = xyplane_array[i].GetReader();
        NodePlane::_node* pnode;
        while((pnode = reader.Next()) != NULL) {
			pnode->p.x -= w_2;
			pnode->p.y -= h_2;
        }
    }
}

void util::TrackSpace::PreRotate(float angle)
{
	float cosa = cos(angle), sina = sin(angle);
	for(int i = 0; i < size; i++) {
        NodePlane::NodeReader reader = xyplane_array[i].GetReader();
        NodePlane::_node* pnode;
        while((pnode = reader.Next()) != NULL) {
			_point p = pnode->p;
			pnode->p.x = cosa*p.x-sina*p.y;
			pnode->p.y = sina*p.x+cosa*p.y;
        }
    }
}
