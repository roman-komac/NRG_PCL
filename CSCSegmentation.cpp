/* std library */
#include <map>
#include <list>
#include <limits>
#include <algorithm>
#include <math.h>
#include <unistd.h>
#include <chrono>
#include <thread>
#include <sstream>

/* pcl library */
#include <pcl/ModelCoefficients.h>
#include <pcl/point_types.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/console/parse.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/visualization/pcl_plotter.h>
#include <pcl/segmentation/supervoxel_clustering.h>

/* Eigen library */
#include <Eigen/SparseCore>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>

/* vtk library */
#include <vtkPolyLine.h>

/* Project files*/
#include "ObjectSegmentation.hpp"

/* CONSTRAINED SPECTRAL CLUSTERING */
class CSCSegmentator : public Segmentator{
	private:
	unsigned int minimum_ui = 1;
	double DELTA_ERROR = 0.000001;
	double rbf_param = 1.0;
	double cutThreshold = 100;
	double wtreshanglesqr = (0.29289 * 0.29289) + (0.70711 * 0.70711); /* 45 degrees difference between vectors squared */
	double wtresh = lin_rbf(wtreshanglesqr, rbf_param);
	static bool sortfunct(std::pair< double,int32_t > p1, std::pair< double,int32_t > p2){p1.first < p2.first;}
	pcl::PointCloud<pcl::PointXYZL> finalCLD;
	public:
	CSCSegmentator() : Segmentator(){
		
	}

	/* in paper: Algorithm 1 */
	std::list < std::multimap<uint32_t,uint32_t> > SSC(std::multimap<uint32_t,uint32_t> graph, std::map<uint32_t, pcl::Supervoxel<pcl::PointXYZRGB>::Ptr > clusters){
		
		// C <- 0     C: final set of clusters
		std::list < std::multimap<uint32_t,uint32_t> >C;

		// Q <- {G}
		std::list < std::multimap<uint32_t,uint32_t> >Q;
		Q.push_back(graph);

		// while sizeof(Q) > 0 do
		while(Q.size() > 0){

			// S <- pop(Q)
			std::multimap<uint32_t,uint32_t> S = Q.front();
			std::map< uint32_t,uint32_t > relabelled, relabelledT;
			/* relabel because the labels may not be sequential */
			relabel(S, relabelled, relabelledT);
			Eigen::SparseMatrix<double> edgeWeight = calc_weights(S, clusters, relabelled, relabelled.size());
			
			// L <- BUILDLAPLACIAN(S)
			Eigen::SparseMatrix<double> L(relabelled.size(), relabelled.size());
			L.reserve(4*clusters.size());
			signedLaplacian(edgeWeight, L);
			Eigen::MatrixXd denseL = Eigen::MatrixXd(L);

			// eigenvalues,eigenvectors <- EIGENDECOMPOSITION(L)
			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(denseL);
			/* transform to vector for further use */
			std::vector< std::vector<double> > eigenvecs;
			EigenToStd(es.eigenvectors(), eigenvecs);

			// A <- {S}
			std::list < std::multimap<uint32_t,uint32_t> >A;
			A.push_back(S);

			// i <- 0      i: cut counter
			int i = 0;

			int cnt = 0;
			Eigen::MatrixXd rvt = es.eigenvalues();
			std::vector< std::vector<double> > evls;
			

			for(int l = 0; l < rvt.rows(); l++){
				if(rvt.coeffRef(l) > DELTA_ERROR){
					cout << rvt.coeffRef(l) << endl;
					cnt = l;
					break;
				}
			}


			// while sizeof(A) > 0 do
			while(A.size() > 0){

				// SG <- pop(A)
				std::multimap<uint32_t,uint32_t> SG = A.front();
				std::multimap<uint32_t,uint32_t>::iterator mmitr;
				mmitr = SG.begin();

				// cut <- FINDBESTCUT(SG, eigenvectors)
				/* returns eigenvector index and position in it */
				for(int z = 0; z < rvt.rows(); z++){
					if(rvt.coeffRef(z) > 0){
						cout << rvt.coeffRef(z) << endl;
						cnt = z;
						break;
					}
				}
				Eigen::MatrixXd xd = es.eigenvectors().block(0, cnt, es.eigenvectors().rows(), std::min(cnt+5, (int)es.eigenvectors().cols() - cnt));
				std::pair<int,int> cut = bestCut(SG, xd, relabelled, relabelledT, edgeWeight);		

				std::vector< std::pair<double, double> > data;

				std::vector< std::pair<int,double> > rsrt = argsort(xd, cut.first);
				for(int f = 0; f < rsrt.size(); f++){
					data.push_back(std::pair<double, double>(f, rsrt[f].second));
				}
				cout << "cut at: " << xd.coeffRef(rsrt[cut.second].first, cut.first) << endl;

				pcl::visualization::PCLPlotter * plotter = new pcl::visualization::PCLPlotter ();
				plotter->addPlotData(data);
				plotter->plot();

				//cout << es.eigenvectors().coeffRef(cut.second, cut.first) << " split" << endl;

				//cout << "cututility: " << CutUtility(SG, cut, relabelled, es.eigenvectors(), edgeWeight) << endl;
				
				if(CutUtility(SG, cut, relabelled, es.eigenvectors(), edgeWeight) > 0.0){
					std::pair< std::multimap<uint32_t,uint32_t>, std::multimap<uint32_t,uint32_t> > prCC = cutGraph(SG, cut, es.eigenvectors(), relabelled, relabelledT);
					if(prCC.first.size() == 0)
						C.push_back(prCC.second);
					else if(prCC.second.size() == 0)
						C.push_back(prCC.first);
					else{
						Q.push_back(prCC.first);
						Q.push_back(prCC.second);
					}

					//cout << "pushed A" << endl;
					i++;
				} else {
					if(i > 0){
						Q.push_back(SG);
						//cout << "pushed Q" << endl;
					}
				}

				A.pop_front();
			}
			if(i == 0){
				C.push_back(S);
				//cout << "pushed C" << endl;
			}

			Q.pop_front();
		}

		return C;
	}

	std::pair< std::multimap<uint32_t,uint32_t>, std::multimap<uint32_t,uint32_t> > cutGraph(std::multimap<uint32_t,uint32_t>& SG, std::pair<int,int> cut, const Eigen::MatrixXd& eigenvecs, std::map< uint32_t,uint32_t >& relabelled, std::map< uint32_t,uint32_t >& relabelledT){
		std::multimap<uint32_t,uint32_t> C1;
		std::multimap<uint32_t,uint32_t> C2;
		std::vector<uint32_t> idx1, idx2;
		for(uint32_t i = 0; i < (uint32_t)eigenvecs.rows(); i++){
				if(eigenvecs.coeffRef(i, cut.first) < eigenvecs.coeffRef(cut.second, cut.first)+DELTA_ERROR){
					idx1.push_back(i);	
				} else {
					idx2.push_back(i);
				}
		}

		std::multimap<uint32_t,uint32_t>::iterator mmit;
		for(uint32_t i = 0; i < idx1.size(); i++){
			std::pair< std::multimap<uint32_t,uint32_t>::iterator, std::multimap<uint32_t,uint32_t>::iterator > prsg = SG.equal_range(relabelledT[idx1[i]]);
			for(mmit = prsg.first; mmit != prsg.second; mmit++){
				if(invec(relabelled[mmit->second], idx1)){
					C1.insert(std::pair<uint32_t, uint32_t>(mmit->first, mmit->second));
				}
			}
		}
		
		for(uint32_t i = 0; i < idx2.size(); i++){
			std::pair< std::multimap<uint32_t,uint32_t>::iterator, std::multimap<uint32_t,uint32_t>::iterator > prsg = SG.equal_range(relabelledT[idx2[i]]);
			for(mmit = prsg.first; mmit != prsg.second; mmit++){
				if(invec(relabelled[mmit->second], idx2)){
					C2.insert(std::pair<uint32_t, uint32_t>(mmit->first, mmit->second));
				}
			}
		}
		return std::pair< std::multimap<uint32_t,uint32_t>, std::multimap<uint32_t,uint32_t> >(C1, C2);
	}

	bool invec(uint32_t val, std::vector<uint32_t>& vec){
		for(int i = 0; i < vec.size(); i++){
			if(vec[i] == val)
				return true;
		}
		return false;
	}

	double CutUtility(std::multimap<uint32_t,uint32_t>& C, std::pair<int,int>& cut, std::map< uint32_t,uint32_t >& relabelled, const Eigen::MatrixXd& eigenvecs, Eigen::SparseMatrix<double>& edgeWeight){
		std::multimap<uint32_t,uint32_t>::iterator mmit;
		double totalCL = 0.0;
		double totalCut = DELTA_ERROR;
		for(mmit = C.begin(); mmit != C.end(); mmit++){
			if(relabelled[mmit->first] < cut.second && relabelled[mmit->second] >= cut.second){
				double ew = edgeWeight.coeffRef(relabelled[mmit->first], relabelled[mmit->second]);
				totalCL += cannot_link(ew, 0.99)? 1 : 0;
				totalCut += ew;
			}
		}

		return totalCL/totalCut;
	}

	std::pair<int,int> bestCut(std::multimap<uint32_t,uint32_t>& C, const Eigen::MatrixXd& eigenvecs, std::map< uint32_t,uint32_t >& relabelled, std::map< uint32_t,uint32_t >& relabelledT, Eigen::SparseMatrix<double>& edgeWeight){
		
		// cuti, cutv, cuts <- 0, 0, 0
		int cuti = 0;
		int cutv = 0;
		double cuts = 0.0;

		// for e in [1 .. k] do
		for(int i = 0; i < eigenvecs.cols(); i++){

			// w <- 0
			double totalw = 0.0;

			// cl <- 0
			int totalcl = 0;

			// for v in argsort(ve) do
			std::vector< std::pair<int,double> > rsrt = argsort(eigenvecs, i);
			for(int v = 0; v < rsrt.size(); v++){

				double inC = 0.0;
				double outC = 0.0;
				for(int a = 0; a <= v; a++){
					for(int b = 0; b <= v; b++){
						if(a != b && isInMultiMap(C, rsrt[a].first, rsrt[b].first)){
							inC += edgeWeight.coeffRef(rsrt[a].first, rsrt[b].first);
						}
					}	
				}
				for(int a = v+1; a < rsrt.size(); a++){
					for(int b = v+1; b < rsrt.size(); b++){
						if(a != b && isInMultiMap(C, rsrt[a].first, rsrt[b].first)){
							outC += edgeWeight.coeffRef(rsrt[a].first, rsrt[b].first);
						}
					}	
				}
				
				/* 2 instead of 1 because we have counted them twice */
				double objectC = 2/inC + 2/outC;

				double scr = signedCut(v, rsrt, edgeWeight, relabelledT, C)*objectC;

				// if (cl/w) > cuts then
				if( scr > cuts ){

					// cuti, cutv, cuts <- e, v, (cl/w)
					cuti = i;
					cutv = v;
					cuts = scr;

				}
			}
		}

		int idxx = 0;
		std::vector< std::pair<int,double> > rsrt = argsort(eigenvecs, cuti);
		double minimprob = std::numeric_limits<double>::max();
		std::multimap<uint32_t, uint32_t>::iterator MMit;
		std::pair< std::multimap<uint32_t, uint32_t>::iterator, std::multimap<uint32_t, uint32_t>::iterator > pos;


		return std::pair<int,int>(cuti, rsrt[cutv].first);
	}

	double signedCut(int idx, std::vector< std::pair<int,double> >& srtd, Eigen::SparseMatrix<double>& edgeWeight, std::map< uint32_t,uint32_t >& relabelledT, std::multimap<uint32_t,uint32_t>& C){
		double scc = 0;
		for(int h = 0; h <= idx; h++){
			for(int i = idx+1; i < srtd.size(); i++){
				if(isInMultiMap(C, relabelledT[srtd[h].first], relabelledT[srtd[i].first])){
					scc += 2 * edgeWeight.coeffRef(srtd[h].first, srtd[i].first);
				}
			}
			 
		}
		return scc;
	}

	bool isInMultiMap(std::multimap<uint32_t,uint32_t>& mm, uint32_t key, uint32_t val){
		std::multimap<uint32_t, uint32_t>::iterator MMit;
		std::pair< std::multimap<uint32_t, uint32_t>::iterator, std::multimap<uint32_t, uint32_t>::iterator > pos = mm.equal_range(key);

		for(MMit = pos.first; MMit != pos.second; MMit++){
			if(MMit->second == val)
				return true;
		}
		return false;
	}

	std::vector< std::pair<int, double> > argsort(const Eigen::MatrixXd& eigenvecs, int idx){
		std::vector< std::pair<int,double> > sorted(eigenvecs.rows());
		for(int i = 0; i < eigenvecs.rows(); i++){
			sorted[i] = std::pair<int, double>(i, eigenvecs.coeffRef(i, idx));
		}

		std::sort(sorted.begin(), sorted.end(),
			[](const std::pair<int, double> & a, const std::pair<int, double> & b) -> bool
		{ 
		    return a.second < b.second; 
		});

		return sorted;
	}

	void EigenToStd(const Eigen::MatrixXd& eigenvecs, std::vector< std::vector<double> >& stdvecs){
		for(int i = 0; i < eigenvecs.cols(); i++){
			stdvecs.push_back(std::vector<double>(eigenvecs.rows()));
			for(int j = 0; j < eigenvecs.rows(); j++){
				/* CoeffRef is row,column */
				stdvecs[i][j] = eigenvecs.coeffRef(j, i);
			}
		}
	}

	bool cannot_link(double ew, double wtr){
		return ew < wtr;
	}

	void relabel(std::multimap<uint32_t,uint32_t> graph, std::map< uint32_t,uint32_t >& lmp, std::map< uint32_t,uint32_t >& lmpT){
		std::multimap<uint32_t,uint32_t>::const_iterator cit;
		std::map< uint32_t,uint32_t >::iterator mit;
		std::pair< std::map<uint32_t,uint32_t>::iterator, bool> ret;
		uint32_t lbl = 0;
		for(cit = graph.begin(); cit != graph.end(); cit++){
			ret = lmp.insert( std::pair<uint32_t,uint32_t>(cit->first, lbl));
			if(ret.second)
				lbl++;
		}

		for(mit = lmp.begin(); mit != lmp.end(); mit++){
			lmpT.insert(std::pair<uint32_t,uint32_t>(mit->second, mit->first));
		}

	}

	Eigen::SparseMatrix<double> calc_weights(std::multimap<uint32_t,uint32_t> supervoxel_adjacency, std::map<uint32_t, pcl::Supervoxel<pcl::PointXYZRGB>::Ptr > supervoxel_clusters, std::map< uint32_t,uint32_t >& lmp, int sz){
		Eigen::SparseMatrix<double> edgeWeight(sz, sz);
		//Povprecna rezervacija prostora
		edgeWeight.reserve(4*supervoxel_clusters.size());
		//Polnjenje matrike povezav grafa

		std::multimap<uint32_t,uint32_t>::iterator label_itr;
		label_itr = supervoxel_adjacency.begin();
		while (label_itr != supervoxel_adjacency.end())
		{
			uint32_t supervoxel_label = label_itr->first;
			pcl::Supervoxel<pcl::PointXYZRGB>::Ptr supervoxel = supervoxel_clusters.at(supervoxel_label);
			pcl::PointNormal pn1;
			supervoxel->getCentroidPointNormal(pn1);

			std::multimap<uint32_t,uint32_t>::iterator adjacent_itr = supervoxel_adjacency.equal_range (supervoxel_label).first;
			std::map<uint32_t, bool> convexity;
			for (; adjacent_itr != supervoxel_adjacency.equal_range(supervoxel_label).second; ++adjacent_itr)
			{
				pcl::Supervoxel<pcl::PointXYZRGB>::Ptr neighbor_supervoxel = supervoxel_clusters.at(adjacent_itr->second);
				pcl::PointNormal pn2;
				neighbor_supervoxel->getCentroidPointNormal(pn2);


					if(convex_edge(pn1, pn2)){
						edgeWeight.coeffRef(lmp[supervoxel_label], lmp[adjacent_itr->second]) = 1;
					} else {
						
						edgeWeight.coeffRef(lmp[supervoxel_label], lmp[adjacent_itr->second]) = lin_rbf(normal_distance_sqr(pn1, pn2), rbf_param);
					}
					if(cannot_link(edgeWeight.coeffRef(lmp[supervoxel_label], lmp[adjacent_itr->second]), wtresh)){
						edgeWeight.coeffRef(lmp[supervoxel_label], lmp[adjacent_itr->second]) = -1;
					}

			}
			label_itr = supervoxel_adjacency.upper_bound (supervoxel_label);
		}

		return edgeWeight;
	}




	/* sign, vkljucuje relativno napako
	*/
	int scalar_sign(pcl::PointNormal& pn1, pcl::PointNormal& pn2){
		float scalar = (pn2.x - pn1.x) * pn2.normal[0] + (pn2.y - pn1.y) * pn2.normal[1] + (pn2.z - pn1.z) * pn2.normal[2];
		return (scalar < -DELTA_ERROR)? -1 : ((scalar > DELTA_ERROR)? 1 : 0);
	}

	/* metrika za konveksnost roba
	*/
	bool convex_edge(pcl::PointNormal& pn1, pcl::PointNormal& pn2){
		return (pn2.x - pn1.x) * (pn2.normal[0] - pn1.normal[0]) + 
		       (pn2.y - pn1.y) * (pn2.normal[1] - pn1.normal[1]) +
		       (pn2.z - pn1.z) * (pn2.normal[2] - pn1.normal[2]) > 0;
	}

	/* evklidska razdalja
	*/
	double euclidean_distance(pcl::PointNormal& pn1, pcl::PointNormal& pn2){
		return sqrt((pn2.x-pn1.x)*(pn2.y-pn1.y)*(pn2.z-pn1.z));
	}

	/* kvadratna evklidska razdalja
	*/
	double euclidean_distance_sqr(pcl::PointNormal& pn1, pcl::PointNormal& pn2){
		return (pn2.x-pn1.x)*(pn2.y-pn1.y)*(pn2.z-pn1.z);
	}

	/* kvadratna evklidska razdalja
	*/
	double normal_distance_sqr(pcl::PointNormal& pn1, pcl::PointNormal& pn2){
		double f = (pn2.normal[0] - pn1.normal[0]);
		double s = (pn2.normal[1] - pn1.normal[1]);
		double t = (pn2.normal[2] - pn1.normal[2]);
		return f*f + s*s + t*t;
	}

	/* krizni produkt normal
	*/
	void normal_cross_product(pcl::PointNormal& pn1, pcl::PointNormal& pn2, double* vec){
		vec[0] = pn1.normal[1] * pn2.normal[2] - pn1.normal[2] * pn2.normal[1];
		vec[1] = pn1.normal[2] * pn2.normal[0] - pn1.normal[0] * pn2.normal[2];
		vec[2] = pn1.normal[0] * pn2.normal[1] - pn1.normal[1] * pn2.normal[0];
	}

	/* produkt vektorjev normal
	*/
	double normal_dot(pcl::PointNormal& pn1, pcl::PointNormal& pn2){
		return pn2.normal[0] * pn1.normal[0] + pn2.normal[1] * pn1.normal[1] + pn2.normal[2] * pn1.normal[2];
	}

	/* produkt vektorjev normal
	*/
	double normal_dot(pcl::PointXYZ& pn1, pcl::PointXYZ& pn2){
		return pn2.x * pn1.x + pn2.y * pn1.y + pn2.z * pn1.z;
	}

	/* krizni produkt normal
	*/
	double normal_angle(pcl::PointXYZ& pn1, pcl::PointXYZ& pn2){
		return acos(normal_dot(pn1, pn2));
	}

	/* rotiranje tocke za 180°
	*/
	pcl::PointXYZ& nflip(pcl::PointXYZ& pt){
		pt.x = -pt.x;
		pt.y = -pt.y;
		pt.z = -pt.z;
		return pt;
	}

	/* kvadratna inverzna radial basis funkcija
	*/
	static double rbf(double c, double par = 1.0){
		return 1.0 / (1.0 + par*c*c);
	}

	/* linearna inverzna radial basis funkcija
	*/
	static double lin_rbf(double c, double par = 1.0){
		return 1.0 / (1.0 + par*c);
	}

	/* Za izris povezav med supervoksli, izdelava poligonske mreze
	*/
	void polygonize (pcl::PointXYZRGBA &supervoxel_center,
                                  pcl::PointCloud<pcl::PointXYZRGBA> &adjacent_supervoxel_centers,
                                  std::string supervoxel_name,
                                  boost::shared_ptr<pcl::visualization::PCLVisualizer> & viewer)
	{
		vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New ();
		vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New ();
		vtkSmartPointer<vtkPolyLine> polyLine = vtkSmartPointer<vtkPolyLine>::New ();

		pcl::PointCloud<pcl::PointXYZRGBA>::iterator adjacent_itr = adjacent_supervoxel_centers.begin ();
		for ( ; adjacent_itr != adjacent_supervoxel_centers.end (); ++adjacent_itr)
		{
			points->InsertNextPoint (supervoxel_center.data);
			points->InsertNextPoint (adjacent_itr->data);
		}
		vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New ();
		polyData->SetPoints (points);
		polyLine->GetPointIds  ()->SetNumberOfIds(points->GetNumberOfPoints ());
		for(unsigned int i = 0; i < points->GetNumberOfPoints (); i++)
			polyLine->GetPointIds ()->SetId (i,i);
		cells->InsertNextCell (polyLine);
		polyData->SetLines (cells);
		viewer->addModelFromPolyData (polyData,supervoxel_name);
	}

	/* Izracun predznacene Laplaceove matrike
	*/
	void signedLaplacian(Eigen::SparseMatrix<double> A, Eigen::SparseMatrix<double>& L){
		//construct D matrix
		//ni vazna cols/rows, ker je kvadratna
		Eigen::SparseMatrix<double> D(A.rows(), A.cols());
		for(int i = 0; i < A.cols(); i++){
			D.coeffRef(i, i) = 0.0;
		}
		for (int k = 0; k < A.outerSize(); ++k){
			for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it){
				D.coeffRef(it.col(), it.col()) += abs(it.value());
			}
		}

		L = D - A;
	}

	void Laplacian(Eigen::SparseMatrix<double> A, Eigen::SparseMatrix<double>& L){
		//construct D matrix
		//ni vazna cols/rows, ker je kvadratna
		Eigen::SparseMatrix<double> D(A.rows(), A.cols());
		for(int i = 0; i < A.cols(); i++){
			D.coeffRef(i, i) = 0.0;
		}
		for (int k = 0; k < A.outerSize(); ++k){
			for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it){
				D.coeffRef(it.col(), it.col()) += it.value();
			}
		}

		L = D - A;
	}

	void segment(pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr& pcldinput){
		unsigned int cntr = 0;

		//supervokselizacija
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr subCloudRGB (new pcl::PointCloud<pcl::PointXYZRGB>);
		pcl::PointCloud<pcl::Normal>::Ptr nrms (new pcl::PointCloud<pcl::Normal>);
		copyPointCloud(*pcld, *subCloudRGB);
		copyPointCloud(*pcld, *nrms);

		//supervokeslizacijski parametri
		float voxel_resolution = 2.0f;
		float seed_resolution = 64.0f;
		float color_importance = 0.2f;
		float spatial_importance = 0.4f;
		float normal_importance = 1.0f;

		pcl::SupervoxelClustering<pcl::PointXYZRGB> super (voxel_resolution, seed_resolution);
		
		//False ce point cloud ne vsebuje podatka o globini(npr. zajem s Kinectom)
		super.setUseSingleCameraTransform(false);
		super.setInputCloud(subCloudRGB);
		super.setNormalCloud(nrms);
		super.setColorImportance(color_importance);
		super.setSpatialImportance(spatial_importance);
		super.setNormalImportance(normal_importance);

		std::map<uint32_t, pcl::Supervoxel<pcl::PointXYZRGB>::Ptr > supervoxel_clusters;

		super.extract(supervoxel_clusters);
		cout << "Found " << supervoxel_clusters.size() << " supervoxels" << endl;
		//super.refineSupervoxels(2, supervoxel_clusters);

		boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer ("3D Viewer"));
		viewer->setBackgroundColor(0, 0, 0);
		


		pcl::PointCloud<pcl::PointXYZRGB>::Ptr voxel_centroid_cloud = super.getVoxelCentroidCloud();
		//viewer->addPointCloud(voxel_centroid_cloud, "voxel centroids");
		//viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE,2.0, "voxel centroids");
		//viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_OPACITY,0.95, "voxel centroids");

		
		pcl::PointCloud<pcl::PointXYZL>::Ptr labeled_cloud = super.getLabeledCloud();
		viewer->addPointCloud(labeled_cloud, "labeled cld");
		while (!viewer->wasStopped ())
		{
		viewer->spinOnce (100);
		}

		pcl::PointCloud<pcl::PointXYZLNormal>::Ptr nrmlabeled(new pcl::PointCloud<pcl::PointXYZLNormal>);
		pcl::concatenateFields(*labeled_cloud, *nrms, *nrmlabeled);
		//viewer->addPointCloud(labeled_cloud, "labeled voxels");
		//viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_OPACITY,0.8, "labeled voxels");
		

		pcl::PointCloud<pcl::PointNormal>::Ptr sv_normal_cloud = super.makeSupervoxelNormalCloud(supervoxel_clusters);
		//viewer->addPointCloudNormals<pcl::PointNormal> (sv_normal_cloud, 1, 5.0f, "supervoxel_normals");

		//cout << "Getting supervoxel adjacency" << endl;
		std::multimap<uint32_t, uint32_t> supervoxel_adjacency;
		super.getSupervoxelAdjacency(supervoxel_adjacency);
		
		//Graf konveksnosti
		std::map<uint32_t, bool> convexity;
		


		//Polnjenje grafa konveksnosti
		std::multimap<uint32_t,uint32_t>::iterator label_itr = supervoxel_adjacency.begin();
		
		while (label_itr != supervoxel_adjacency.end())
		{
			uint32_t supervoxel_label = label_itr->first;
			pcl::Supervoxel<pcl::PointXYZRGB>::Ptr supervoxel = supervoxel_clusters.at(supervoxel_label);
			pcl::PointNormal pn1;
			supervoxel->getCentroidPointNormal(pn1);

			std::multimap<uint32_t,uint32_t>::iterator adjacent_itr = supervoxel_adjacency.equal_range (supervoxel_label).first;
			std::map<uint32_t, bool> convexity;
			int connects = 0;
			for (; adjacent_itr != supervoxel_adjacency.equal_range(supervoxel_label).second; ++adjacent_itr)
			{
				pcl::Supervoxel<pcl::PointXYZRGB>::Ptr neighbor_supervoxel = supervoxel_clusters.at (adjacent_itr->second);
				pcl::PointNormal pn2;
				neighbor_supervoxel->getCentroidPointNormal(pn2);
				connects += scalar_sign(pn1, pn2);
			}
			convexity.insert(std::pair<uint32_t, bool>(supervoxel_label, connects > 0));
			label_itr = supervoxel_adjacency.upper_bound (supervoxel_label);
		}
		//velikost = last index+1
		uint32_t maxlabel = super.getMaxLabel()+1;



		std::list < std::multimap<uint32_t,uint32_t> > scsc = SSC(supervoxel_adjacency, supervoxel_clusters);
		//cout << "before remapping" << endl;

		pcl::PointCloud<pcl::PointXYZL>::iterator pcit;

		std::list < std::multimap<uint32_t,uint32_t> >::iterator lmmit;
		uint32_t labeler = 0;
		
		std::multimap<uint32_t,uint32_t>::iterator mmiter;
		for(lmmit = scsc.begin(); lmmit != scsc.end(); lmmit++){
			mmiter = (*lmmit).begin();
			while( mmiter != (*lmmit).end()){
				for(pcit = labeled_cloud->begin(); pcit != labeled_cloud->end(); pcit++){
					if(pcit->label == mmiter->first){
						pcl::PointXYZL clpt;
						clpt.x = pcit->x;
						clpt.y = pcit->y;
						clpt.z = pcit->z;
						clpt.label = labeler;
						finalCLD.push_back(clpt);
					}
				}
				mmiter = (*lmmit).upper_bound(mmiter->first);
				
			}
			labeler++;
		}


	}
	void setMinimumInliers(unsigned int ui){
		minimum_ui = ui;
	}
	void saveSegments(std::string fname){
		pcl::PointCloud<pcl::PointXYZL>::iterator litr;
		std::map<uint32_t, pcl::PointCloud<pcl::PointXYZL> > mp;
		for(litr = finalCLD.begin(); litr != finalCLD.end(); litr++){
			if(mp.find(litr->label) == mp.end()){
				mp[litr->label] = pcl::PointCloud<pcl::PointXYZL>();
			}
			mp[litr->label].push_back(*litr);
		}
		std::map<uint32_t, pcl::PointCloud<pcl::PointXYZL> >::iterator mitr;
		for(mitr = mp.begin(); mitr != mp.end(); mitr++){
			std::ostringstream ss;
			ss << mitr->first << fname;
			pcl::io::savePLYFileBinary(ss.str(), mitr->second);
		}
		
	}
	void viewSegments(){
		boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer ("3D Viewer"));
		viewer->setBackgroundColor(0, 0, 0);
		pcl::PointCloud<pcl::PointXYZL>::Ptr pr(&finalCLD);
		viewer->addPointCloud(pr, "labeled voxels");

		while (!viewer->wasStopped ())
		{
		viewer->spinOnce (100);
		}
		return;
	}
};



