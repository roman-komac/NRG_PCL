#pragma once

/* std library */
#include <string>
#include <iostream>
#include <type_traits>

/* pcl library */
#include <pcl/point_cloud.h>
#include <pcl/PCLPointCloud2.h>
#include <pcl/filters/uniform_sampling.h>
#include <pcl/cloud_iterator.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/io/ply_io.h>

class Segmentator{
	protected:
	pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr pcld;
	pcl::PointCloud<pcl::PointXYZL>::Ptr final;
	std::vector<unsigned int> labels;
	public:
	void setInputCloud(pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr& pcldinput){
		pcld = pcldinput;
		labels = std::vector<unsigned int>(pcld->size());
	}
	void segment(pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr& pcldinput){

	}
	void saveSegments(std::string fname){
		pcl::io::savePLYFileBinary(fname, *final);
	}
	void viewSegments(){

	}
};