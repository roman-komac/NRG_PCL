#include <string>
#include <iostream>
#include "ObjectSegmentation.hpp"
#include "CSCSegmentation.cpp"

using namespace std;
using namespace pcl;

int main(int argv, char** argc){

	PointCloud<PointXYZRGBNormal>::Ptr cloud (new PointCloud<PointXYZRGBNormal>);
	PointCloud<PointXYZRGBNormal>::Ptr subsampledCloud (new PointCloud<PointXYZRGBNormal>);
	PointCloud<PointXYZRGBNormal>::Ptr segmentedCloud (new PointCloud<PointXYZRGBNormal>);
	string base("");

	if(argv < 2){
		cout << "No point cloud file specified." << endl;
		cout << "Please input as first argument." << endl;
		cout << "Example:" << endl;
		cout << string(argc[0]) << " pointcloud.ply" << endl;
		return -1;
	}

	cout << "Reading point cloud file: " << string(argc[1]) << endl;
	if(io::loadPLYFile<pcl::PointXYZRGBNormal> ( string(argc[1]), *cloud) == -1){
		cerr << "Error reading point cloud" << endl;
		return -1;
	}

	cout << "Point cloud size before subsampling: " << cloud->size() << " points" << endl;

	UniformSampling<PointXYZRGBNormal>::Ptr sampler(new UniformSampling<PointXYZRGBNormal>);
	sampler->setInputCloud(cloud);
	sampler->setRadiusSearch(3);
	sampler->filter(*subsampledCloud);

	cout << "Point cloud size after subsampling: " << subsampledCloud->size() << " points" << endl;

	double middlex = 0;
	double middley = 0;
	double middlez = 0;
	PointCloud<PointXYZRGBNormal>::iterator cit;
	for(cit = subsampledCloud->begin(); cit != subsampledCloud->end(); cit++){
		middlex += cit->x;
		middley += cit->y;
		middlez += cit->z;
	}

	middlex /= subsampledCloud->size();
	middley /= subsampledCloud->size();
	middlez /= subsampledCloud->size();

	for(cit = subsampledCloud->begin(); cit != subsampledCloud->end(); cit++){
		cit->x -= middlex;
		cit->y -= middley;
		cit->z -= middlez;
	}

	cout << "Centered point cloud" << endl;

	PointCloud<PointXYZRGB>::Ptr subCloudRGB (new PointCloud<PointXYZRGB>);
	PointCloud<PointXYZRGBNormal>::Ptr segCloudRGB (new PointCloud<PointXYZRGBNormal>);

	PointCloud<PointXYZRGB> cloud_xyzrgb;
	copyPointCloud(*subsampledCloud, *subCloudRGB);

	CSCSegmentator segmentator;
	segmentator.setInputCloud(subsampledCloud);
	//Manj k 100 so posamezne krpe
	segmentator.setMinimumInliers(100);
	segmentator.segment(segCloudRGB);
	segmentator.saveSegments("pointcloud_save.ply");
	segmentator.viewSegments();

}