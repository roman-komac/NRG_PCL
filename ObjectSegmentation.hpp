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

/* templated trio (same as std::pair<>) but holds three values */
template <class c1, class c2, class c3>
struct trio{
	c1 first;
	c2 second;
	c3 third;
	trio(c1 fCl, c2 sCl, c3 tCl) : first(fCl), second(sCl), third(tCl) {}
};

class Segmentator{
	protected:
	pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr pcld;
	std::vector<unsigned int> labels;
	public:
	void setInputCloud(pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr& pcldinput){
		pcld = pcldinput;
		labels = std::vector<unsigned int>(pcld->size());
	}
	void segment(pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr& pcldinput){

	}
};