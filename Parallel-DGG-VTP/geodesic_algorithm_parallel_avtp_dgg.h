// Parallel-VTP 

#ifndef GEODESIC_ALGORITHM_PARALLEL_AVTP_DGG_H
#define GEODESIC_ALGORITHM_PARALLEL_AVTP_DGG_H

#include "stdafx.h"
#include "geodesic_memory.h"
#include "geodesic_algorithm_base.h"
#include "geodesic_algorithm_exact_elements.h"
#include "FWP_Queue.h"

namespace geodesic{

	class GeodesicAlgorithmParallelAVTPDGG :public GeodesicAlgorithmBase
	{
	public:
		GeodesicAlgorithmParallelAVTPDGG(geodesic::Mesh* mesh, unsigned num_procs, unsigned k_concurrent, double eps, double side_fan_const) :
			GeodesicAlgorithmBase(mesh),
			m_vertex_FWP_queue(mesh->vertices().size(), 1),
			num_procs(num_procs),
			mx_concurrent_list_num(k_concurrent),
			epsilon(eps),
			side_fan_angle(side_fan_const*sqrt(eps))
		{};
		~GeodesicAlgorithmParallelAVTPDGG(){};


		void clear(); 
		void propagate(unsigned source);
		void print_statistics();


	private:


		void initialize_propagation_data();
		void create_pseudo_source_windows(vertex_pointer &v, bool UpdataQueue);
		void erase_from_FWP_queue(vertex_pointer v);

		//Rule 2: pairwindow checking
		void check_with_vertices(list_pointer &list);
		windows_state check_between_two_windows(interval_pointer &w1, interval_pointer &w2, Triangle &tri);
		void pairwise_windows_checking(list_pointer &list, Triangle &tri);

		//Rule 1:
		void find_separating_point(list_pointer &list, Triangle &tri, Vertex &top_t);
		void propagate_windows_to_two_edges(list_pointer &list, Triangle &tri, IntervalList &wl_left, IntervalList &wl_right);
		
		void updateDistanceBound(interval_pointer &w); // update the distance bound for a window based on FastDGG Rules..
		void propagate_from_pseudo_source(vertex_pointer pseudo_source, bool inside_traversed_area, face_pointer to_face);
		void create_flat_v_fan_shape(vertex_pointer &flat_v, bool UpdateFIFOQueue); // propagate fan from spehrical / euclidena vertex (i.e. "flat" when the faces around the vertex are unfolded)
		void updateDistance(vertex_pointer&pseudo_source);

		unsigned int m_source;

		unsigned prevPhase;
		FWP_Queue<vertex_pointer> m_vertex_FWP_queue;
		
		unsigned from_list, to_list;
		std::vector<std::vector<list_pointer>> list_queue;
		//list_queue[0]: from_queue, list_queue[1]: next_queue.

		//for propagate lists concurrently, update vertices and merge lists sequentially
		std::vector<Triangle> tris;
		std::vector<IntervalList> wl_lefts;
		std::vector<IntervalList> wl_rights;
		std::vector<Vertex> top_ts;
		std::vector<int> isNotBoundary;
		std::vector<vertex_pointer> vertextoprop;
		unsigned estimate_wavefront;
		
		//K window lists selection
		unsigned mx_concurrent_list_num;
		unsigned concurrent_num;
		//double tau_value;

		//----------for statistics-----
		unsigned num_procs;
		unsigned num_propagated_vertices;
		unsigned m_iterations;
		double list_concurrent_time;
		double epsilon;
		double side_fan_angle;
		double clock_neededoperations = 0;
		double clockfanshape = 0;
	};

	inline void GeodesicAlgorithmParallelAVTPDGG::clear()
	{
		GeodesicAlgorithmBase::clear();

		
		//limit_list_num = 5000;
		from_list = 0;
		to_list = 0;
		list_queue.resize(2);//list_queue[0]: from_queue, list_queue[1]: next_queue.
		list_queue[0].clear();
		list_queue[1].clear();
		estimate_wavefront = 0;
		tris.resize(estimate_wavefront);
		wl_lefts.resize(estimate_wavefront);
		wl_rights.resize(estimate_wavefront);
		top_ts.resize(estimate_wavefront);
		isNotBoundary.resize(estimate_wavefront, 0);

	}

	inline void GeodesicAlgorithmParallelAVTPDGG::propagate_from_pseudo_source(vertex_pointer pseudo_source, bool inside_traversed_area, face_pointer to_face){

		//m_windows_pseudo++;
		//clock_t start_clock = clock();
		face_pointer face_it = to_face;
		edge_pointer edge_it = face_it->opposite_edge(pseudo_source);

		list_pointer list = (edge_it->adjacent_faces()[0] == to_face) ? list = interval_list_0(edge_it) : list = interval_list_1(edge_it);

		// create a window
		interval_pointer candidate = new Interval;
		//if (start < 0 || end > edge_it->length())return;
		/*if (start_vertex->id() == list->start_vertex()->id()){
		candidate->start() = start;
		candidate->stop() = end;
		}
		else{
		candidate->start() = edge_it->length() - end;
		candidate->stop() = edge_it->length() - start;
		}*/
		candidate->start() = 0;
		candidate->stop() = edge_it->length();

		candidate->d() = pseudo_source->geodesic_distance();
		double angle = face_it->vertex_angle(list->start_vertex());
		double length = face_it->next_edge(edge_it, list->start_vertex())->length();
		candidate->pseudo_x() = cos(angle) * length;
		candidate->pseudo_y() = -sin(angle) * length;

		edge_pointer to_start_vertex = face_it->next_edge(edge_it, list->start_vertex());
		candidate->left_CH_distance() = to_start_vertex->length();
		candidate->right_CH_distance() = face_it->next_edge(to_start_vertex, pseudo_source)->length();
		//if (candidate->start() == 0)
		//		candidate->left_cosine_CH_angle() = 1;
		/*else{
		double vector_start_x, vector_start_y, vector_left_CH_x,vector_left_CH_y;
		vector_start_x = candidate->start() - candidate->pseudo_x();
		vector_left_CH_x = 0 - candidate->pseudo_x();
		vector_left_CH_y = vector_start_y = 0 - candidate->pseudo_y();
		candidate->left_cosine_CH_angle() = compute_cosine_vector(vector_start_x,vector_start_y,vector_left_CH_x,vector_left_CH_y);
		}*/
		//if (candidate->stop() == edge_it->length())
		//	candidate->right_cosine_CH_angle() = 1;
		/*else{
		double vector_stop_x, vector_stop_y, vector_right_CH_x, vector_right_CH_y;
		vector_stop_x = candidate->stop() - candidate->pseudo_x();
		vector_right_CH_x = edge_it->length() - candidate->pseudo_x();
		vector_right_CH_y = vector_stop_y = 0 - candidate->pseudo_y();
		candidate->right_cosine_CH_angle() = compute_cosine_vector(vector_stop_x, vector_stop_y, vector_right_CH_x, vector_right_CH_y);
		}*/



		//		candidate->origin_vertex() = pseudo_source;
		candidate->distance_bound() = 1e10;
		//updateDistanceBound(candidate);
		// insert into list
		list->push_back(candidate);

		// Statistics
		++m_windows_wavefront;
		if (m_windows_peak < m_windows_wavefront)
			m_windows_peak = m_windows_wavefront;
		//clock_t stop_clock = clock();
		//clockpseudoprop += static_cast<double>(stop_clock - start_clock);
	}

	inline void GeodesicAlgorithmParallelAVTPDGG::updateDistance(vertex_pointer&pseudo_source){
		for (unsigned i = 0; i < pseudo_source->adjacent_edges().size(); ++i)
		{
			edge_pointer edge_it = pseudo_source->adjacent_edges()[i];
			vertex_pointer vert_it = edge_it->opposite_vertex(pseudo_source);
			double distance = pseudo_source->geodesic_distance() + edge_it->length();
			if (distance < vert_it->geodesic_distance())
			{
				erase_from_FWP_queue(vert_it);
				vert_it->geodesic_distance() = distance;

				if (vert_it->state() == Vertex::OUTSIDE)
					vert_it->state() = Vertex::FRONT;

				vert_it->incident_face() = edge_it->adjacent_faces()[0];
				edge_pointer next_edge = vert_it->incident_face()->next_edge(edge_it, pseudo_source);
				//vert_it->incident_point() = (next_edge->v0() == pseudo_source) ? 0 : next_edge->length();
				vert_it->edge_start_vertex() = edge_it;
				vert_it->p_start_angle() = 0;
				vert_it->ptrInQueue = m_vertex_FWP_queue.push(vert_it, vert_it->geodesic_distance());

			}
		}
	}

	inline void GeodesicAlgorithmParallelAVTPDGG::create_flat_v_fan_shape(vertex_pointer &pseudo_source, bool inside_traversed_area){

		clock_t start = clock();

		double totalanglearound = pseudo_source->sum_angle_around_v();
		if (totalanglearound <= (M_PI))return;
		face_pointer incident_face = pseudo_source->incident_face();
		edge_pointer edge_on_start_vertex = pseudo_source->edge_start_vertex();
		/*short int index = incident_face->face_index_around_vertex().at(pseudo_source->id());
		face_pointer nexttoincidentface = pseudo_source->edge_start_vertex()->opposite_face(incident_face);
		short int indexminusone, indexplusone;

		if (index == pseudo_source->adjacent_faces().size() - 1){
		indexplusone = 0;
		}
		else indexplusone = index + 1;
		if (nexttoincidentface->id() == pseudo_source->adjacent_faces()[indexplusone]->id()){
		pseudo_source->p_start_angle() = incident_face->vertex_angle(pseudo_source) - pseudo_source->p_start_angle();
		//printf("Wrong rotation\n");
		}
		//else printf("completelu mistaken\n");

		double angleuntilgeod = angle_sum_until(pseudo_source->id(), index) + pseudo_source->p_start_angle();
		double anglefanone, anglefantwo, angledirectone, angledirecttwo;
		anglefanone = M_PI - side_fan_angle; anglefantwo = -anglefanone + totalanglearound;
		angledirectone = M_PI+side_fan_angle; angledirecttwo = -angledirectone + totalanglearound;
		double anglestarttemp = min({ anglefanone, anglefantwo, angledirectone, angledirecttwo });
		double anglestart = anglestarttemp + angleuntilgeod;
		if (anglestart > pseudo_source->sum_angle_around_v())anglestart -= pseudo_source->sum_angle_around_v();
		vector<double>& angle_around_pseudo = one_vertex_angle_around(pseudo_source->id());
		index = lower_bound(angle_around_pseudo.begin(), angle_around_pseudo.end(), anglestart) - angle_around_pseudo.begin();
		double angleend = max({anglefanone, angledirectone, anglefantwo, angledirecttwo}) - anglestarttemp + anglestart;

		while (anglestart < angleend){
		//printf("%lf %lf %d\n", anglestart, angleend, index);
		if (index >= pseudo_source->adjacent_faces().size()){
		index = 0;
		}
		else if (index < 0)printf("ERror indexing\n");
		face_pointer to_face = pseudo_source->adjacent_faces()[index];
		propagate_from_pseudo_source(pseudo_source, inside_traversed_area, to_face);
		anglestart += to_face->vertex_angle(pseudo_source);
		index++;


		}
		for (unsigned i = 0; i < pseudo_source->adjacent_edges().size(); ++i)
		{
		edge_pointer   edge_it = pseudo_source->adjacent_edges()[i];
		vertex_pointer vert_it = edge_it->opposite_vertex(pseudo_source);

		double distance = pseudo_source->geodesic_distance() + edge_it->length();
		if (distance < vert_it->geodesic_distance())
		{
		m_vertex_queue.erase(vert_it);

		vert_it->geodesic_distance() = distance;
		if (vert_it->state() == Vertex::OUTSIDE)
		vert_it->state() = Vertex::FRONT;

		vert_it->incident_face() = edge_it->adjacent_faces()[0];
		edge_pointer next_edge = vert_it->incident_face()->next_edge(edge_it, pseudo_source);
		vert_it->incident_point() = (next_edge->v0() == pseudo_source) ? 0 : next_edge->length();
		vert_it->p_start_angle() = 0;
		vert_it->edge_start_vertex() = edge_it;
		m_vertex_queue.insert(vert_it);
		}
		}*/


		//cin.get();
		edge_pointer edge_on_end_vertex = incident_face->next_edge(edge_on_start_vertex, pseudo_source);
		double p_start_len = edge_on_start_vertex->length();
		double sum_angle_start_dir = pseudo_source->p_start_angle();
		double sum_angle_end_dir = incident_face->vertex_angle(pseudo_source) - sum_angle_start_dir;

		//assert(0 <= sum_angle_start_dir <= (1 + 1e-10)*incident_face->vertex_angle(pseudo_source));

		double sum_angle = sum_angle_start_dir;
		face_pointer incident_it = incident_face;
		double total_angle_around = pseudo_source->sum_angle_around_v();
		if (total_angle_around <= (M_PI - side_fan_angle)){
			//create_pseudo_source_windows(pseudo_source, inside_traversed_area);
			return;
		}
		double start_fan, stop_fan;
		start_fan = M_PI - side_fan_angle;

		stop_fan = total_angle_around - (M_PI - side_fan_angle);

		if (start_fan > stop_fan){
			std::swap(start_fan, stop_fan);
			if (stop_fan < M_PI + side_fan_angle)
				stop_fan = M_PI + side_fan_angle;
			if (start_fan > total_angle_around - (M_PI + side_fan_angle))
				start_fan = total_angle_around - (M_PI + side_fan_angle);
		}
		double stop_fan_opposite_directions = total_angle_around - stop_fan;

		while (sum_angle <= start_fan){
			incident_it = edge_on_start_vertex->opposite_face(incident_it);
			sum_angle += incident_it->vertex_angle(pseudo_source);
			edge_on_start_vertex = incident_it->next_edge(edge_on_start_vertex, pseudo_source);

		}
		edge_on_start_vertex = incident_it->next_edge(edge_on_start_vertex, pseudo_source);
		face_pointer f_from_start_vertex = incident_it;

		double angle_f_start = f_from_start_vertex->vertex_angle(pseudo_source) - (sum_angle - start_fan);



		incident_it = incident_face;
		sum_angle = sum_angle_end_dir;
		while (sum_angle <= stop_fan_opposite_directions){
			incident_it = edge_on_end_vertex->opposite_face(incident_it);
			sum_angle += incident_it->vertex_angle(pseudo_source);
			edge_on_end_vertex = incident_it->next_edge(edge_on_end_vertex, pseudo_source);

		}
		face_pointer f_from_end_vertex = incident_it;
		edge_on_end_vertex = incident_it->next_edge(edge_on_end_vertex, pseudo_source);
		double angle_f_end = f_from_end_vertex->vertex_angle(pseudo_source) - (sum_angle - stop_fan_opposite_directions);

		// update vertices around pseudo_source

		


		// propagate from here....
		if (f_from_end_vertex->id() == f_from_start_vertex->id()){
			vertex_pointer start_vertex = edge_on_start_vertex->opposite_vertex(pseudo_source);
			vertex_pointer end_vertex = edge_on_end_vertex->opposite_vertex(pseudo_source);

			//double mid_start_angle = M_PI - angle_f_start - f_from_start_vertex->vertex_angle(start_vertex);
			//double start = edge_on_start_vertex->length()*sin(angle_f_start) / sin(mid_start_angle);
			//double mid_end_angle = M_PI - angle_f_end - f_from_start_vertex->vertex_angle(end_vertex);
			//double end = edge_on_end_vertex->length()*sin(angle_f_end) / sin(mid_end_angle);
			//edge_pointer edge_opp = f_from_start_vertex->opposite_edge(pseudo_source);
			//end = edge_opp->length() - end;
			//assert((1 + 1e-10)*start <= end);

			propagate_from_pseudo_source(pseudo_source, inside_traversed_area, f_from_start_vertex);
		}
		else{
			vertex_pointer start_vertex = edge_on_start_vertex->opposite_vertex(pseudo_source);
			//double mid_start_angle = M_PI - angle_f_start - f_from_start_vertex->vertex_angle(start_vertex);
			//double start = edge_on_start_vertex->length()*sin(angle_f_start) / sin(mid_start_angle);
			edge_pointer edge_opp = f_from_start_vertex->opposite_edge(pseudo_source);
			//double end = edge_opp->length();
			//assert((1 + 1e-10)*start <= end);
			//if (start*(1 + 1e-6) < end)
			propagate_from_pseudo_source(pseudo_source, inside_traversed_area, f_from_start_vertex);

			incident_it = f_from_start_vertex;
			edge_on_start_vertex = incident_it->next_edge(edge_on_start_vertex, pseudo_source);
			incident_it = edge_on_start_vertex->opposite_face(incident_it);
			while (incident_it->id() != f_from_end_vertex->id()){
				// propagate from start to end face

				edge_pointer edge_it = incident_it->opposite_edge(pseudo_source);
				double end = edge_it->length();
				propagate_from_pseudo_source(pseudo_source, inside_traversed_area, incident_it);
				edge_on_start_vertex = incident_it->next_edge(edge_on_start_vertex, pseudo_source);
				incident_it = edge_on_start_vertex->opposite_face(incident_it);

			}
			start_vertex = edge_on_end_vertex->opposite_vertex(pseudo_source);
			//double mid_end_angle = M_PI - angle_f_end - f_from_end_vertex->vertex_angle(start_vertex);
			//start = edge_on_end_vertex->length()*sin(angle_f_end) / sin(mid_end_angle);
			//end = f_from_end_vertex->opposite_edge(pseudo_source)->length();
			propagate_from_pseudo_source(pseudo_source, inside_traversed_area, f_from_end_vertex);

		}
				clock_t stop = clock();
				clockfanshape += static_cast<double>(stop - start);

	}

	inline void GeodesicAlgorithmParallelAVTPDGG::updateDistanceBound(interval_pointer& w){

		//clock_t start = clock();
		double sourceXCoord = w->pseudo_x();
		double sourceYCoord = w->pseudo_y();

		//double edgelength = model.Edge(w.indexOfCurEdge).length;
		double windowLeftXCoord = w->start();
		double windowRightXCoord = w->stop();
		double RightXDifference = windowRightXCoord - sourceXCoord; // difference between point b1 to source in term of x coordinate
		double LeftXDifference = windowLeftXCoord - sourceXCoord; // difference between point b0 to source in term of x coordinate
		/*	double tanRightDiff = abs(RightXDifference) / abs(sourceYCoord);
		double tanLeftDIff = abs(LeftXDifference) / abs(sourceYCoord);
		double tanLeftCH = sqrtf(1 - w->left_cosine_CH_angle()*w->left_cosine_CH_angle())/w->left_cosine_CH_angle();
		*/
		double LeftDistanceSquared = (LeftXDifference)*(LeftXDifference)+sourceYCoord*sourceYCoord; // (distance of source to b0)^2.
		double RightDistanceSquared = (RightXDifference)*(RightXDifference)+sourceYCoord*sourceYCoord; // (distance of source to b1)^2.
		double isoscelesBottomSquared;

		double RightDistance;
		//if (w->right_CH_distance() == w.rightLen)RightDistance = w.rightLen;
		RightDistance = sqrtf(RightDistanceSquared);
		double LeftDistance;
		//if (w.distanceLeftCHVertex == w.leftLen) LeftDistance = w.leftLen;
		LeftDistance = sqrtf(LeftDistanceSquared);

		double triangledistance;
		if (RightDistanceSquared < LeftDistanceSquared){ // extend right hand side of the window since it's the one that's shorter
			float ratio = (LeftDistance / RightDistance);
			double isoscelesRightProjectionX = sourceXCoord + ratio*RightXDifference;
			double isoscelesRightProjectionY = sourceYCoord - ratio*(sourceYCoord);
			isoscelesBottomSquared = (isoscelesRightProjectionX - windowLeftXCoord)*(isoscelesRightProjectionX - windowLeftXCoord)
				+ isoscelesRightProjectionY*isoscelesRightProjectionY; // (isosceles triangle base length)^2
			//isoscelesAltitude = sqrtf(LeftDistanceSquared - 0.25*isoscelesBottomSquared); // calculate triangle altitude
			triangledistance = LeftDistance;
		}
		else{ // extend left hand side of the window since it's the one that's shorter
			float ratio = RightDistance / LeftDistance;
			double isoscelesLeftProjectionX = sourceXCoord + ratio*LeftXDifference;
			double isoscelesLeftProjectionY = sourceYCoord - ratio*(sourceYCoord);
			isoscelesBottomSquared = (isoscelesLeftProjectionX - windowRightXCoord)*(isoscelesLeftProjectionX - windowRightXCoord)
				+ isoscelesLeftProjectionY*isoscelesLeftProjectionY; // (isosceles triangle base length)^2
			//isoscelesAltitude = sqrtf(RightDistanceSquared - 0.25*isoscelesBottomSquared); // calculate triangle altitude
			triangledistance = RightDistance;
		}
		//double bound = isoscelesBottomSquared / (8 * epsilon*isoscelesAltitude) + isoscelesAltitude;
		//if (bound < w.maxDistance) w.maxDistance = bound;
		if (isoscelesBottomSquared >= 2 * triangledistance*triangledistance)return; // check if angle > 90 degrees
		double CHdistance_ratio, largerdistance;
		if (w->left_CH_distance() < w->right_CH_distance()){
			CHdistance_ratio = w->left_CH_distance() / w->right_CH_distance();
			largerdistance = w->right_CH_distance();
		}
		else{
			CHdistance_ratio = w->right_CH_distance() / w->left_CH_distance();
			largerdistance = w->left_CH_distance();
		}
		double beta, bottom, oneminratio;
		double ScaleFactor = 2.25;
		oneminratio = 1 / (1 - CHdistance_ratio);
		bottom = (largerdistance / triangledistance)*sqrtf(isoscelesBottomSquared);
		if (0.99 < CHdistance_ratio && CHdistance_ratio < 1.01){ beta = 0.5; }
		else{
			beta = oneminratio - sqrtf(oneminratio*(oneminratio - 1) - largerdistance*(largerdistance - bottom) * 2 * ScaleFactor* epsilon / (bottom*bottom));
		}
		if (beta <= 0 || beta >= 1 || isnan(beta))beta = 1;
		double beta_far, beta_near;
		double triangle_altitude, cosine_half_winopening, sine_half_winopening;
		double bound;
		double intersect_left, intersect_right;
		

		bound = (beta)*(beta)*CHdistance_ratio*bottom*bottom / (2 * ScaleFactor* epsilon*(largerdistance - bottom)) + CHdistance_ratio*largerdistance;
		if (bound < w->distance_bound() && largerdistance > bottom)w->distance_bound() = bound;
		//clock_t stop = clock();
		//clockcountbound += static_cast<double>(stop)-static_cast<double>(start);
	}

	inline void GeodesicAlgorithmParallelAVTPDGG::initialize_propagation_data()
	{
		clear();
		clock_neededoperations = 0;
		//initialize source's parameters
		vertex_pointer source = &(this->mesh()->vertices()[m_source]);
		source->geodesic_distance() = 0;
		source->state() = Vertex::INSIDE;


		//for fwp's 'bucket' structure
		unsigned int initialBinSize = 0;
		initialBinSize = source->adjacent_edges().size();
		m_vertex_FWP_queue.setParameters1(binWidth, initialBinSize);

		//initialize windows around pai
		create_pseudo_source_windows(source, false);

		m_vertex_FWP_queue.setParameters2(Kmin, Kmax, step);
		m_vertex_FWP_queue.setBinRange();
		m_vertex_FWP_queue.setPhase(1);
	}

	inline void GeodesicAlgorithmParallelAVTPDGG::create_pseudo_source_windows(vertex_pointer &pseudo_source, bool inside_traversed_area)
	{
		for (unsigned i = 0; i < pseudo_source->adjacent_edges().size(); i++)
		{
			edge_pointer edge_it = pseudo_source->adjacent_edges()[i];
			vertex_pointer vert_it = edge_it->opposite_vertex(pseudo_source);
			double distance = pseudo_source->geodesic_distance() + edge_it->length();
			if (distance < vert_it->geodesic_distance())
			{
				erase_from_FWP_queue(vert_it);
				vert_it->geodesic_distance() = distance;

				if (vert_it->state() == Vertex::OUTSIDE)
					vert_it->state() = Vertex::FRONT;

				vert_it->incident_face() = edge_it->adjacent_faces()[0];
				//edge_pointer next_edge = vert_it->incident_face()->next_edge(edge_it, pseudo_source);
				edge_pointer next_edge = vert_it->incident_face()->next_edge(edge_it, pseudo_source);
//				vert_it->incident_point() = (next_edge->v0() == pseudo_source) ? 0 : next_edge->length();
				vert_it->edge_start_vertex() = edge_it;
				vert_it->p_start_angle() = 0;
				vert_it->ptrInQueue = m_vertex_FWP_queue.push(vert_it, vert_it->geodesic_distance());
				
			}
		}
		for (unsigned i = 0; i < pseudo_source->adjacent_faces().size(); i++)
		{
			face_pointer face_it = pseudo_source->adjacent_faces()[i];
			edge_pointer edge_it = face_it->opposite_edge(pseudo_source);
			list_pointer list = (edge_it->adjacent_faces()[0] == face_it) ? list = interval_list_0(edge_it) : list = interval_list_1(edge_it);

			//create a window
			interval_pointer candidate = new Interval;
			candidate->start() = 0;
			candidate->stop() = edge_it->length();
			candidate->d() = pseudo_source->geodesic_distance();
			double angle = face_it->vertex_angle(list->start_vertex());
			double length = face_it->next_edge(edge_it, list->start_vertex())->length();
			candidate->pseudo_x() = cos(angle)*length;
			candidate->pseudo_y() = -sin(angle)*length;
			edge_pointer left_edge = face_it->next_edge(edge_it, list->start_vertex());
			candidate->left_CH_distance() = left_edge->length();
			//			candidate->left_cosine_CH_angle() = 1;
			candidate->right_CH_distance() = face_it->next_edge(left_edge, pseudo_source)->length();
			candidate->distance_bound() = 1e10;
			//insert into list
			list->push_back(candidate);

			// Statistics
			++m_windows_wavefront;
			if (m_windows_peak < m_windows_wavefront)
				m_windows_peak = m_windows_wavefront;
		}
	}

	inline void GeodesicAlgorithmParallelAVTPDGG::erase_from_FWP_queue(vertex_pointer v)
	{
		if (v->ptrInQueue)
		{
			m_vertex_FWP_queue.remove(v->ptrInQueue);
			v->ptrInQueue = NULL;
		}
	}

	//Rule 2: pairwise windows checking : case 12 ICH strategy
	inline void GeodesicAlgorithmParallelAVTPDGG::check_with_vertices(list_pointer &list)
	{
		if (list->empty()) return;
		interval_pointer iter = list->begin();
		interval_pointer iter_t;
		while ((!list->empty()) && (iter != NULL))
		{
			interval_pointer &w = iter;
			
			bool w_survive = true;

			if (w->stop() - w->start() < Global_eps)
				w_survive = false;

			edge_pointer e = list->edge();
			vertex_pointer v1 = list->start_vertex();
			vertex_pointer v2 = e->opposite_vertex(v1);
			double d1 = GEODESIC_INF;
			
			d1 = w->d() + sqrt((w->stop() - w->pseudo_x()) * (w->stop() - w->pseudo_x()) + w->pseudo_y() * w->pseudo_y());
			if (v1->geodesic_distance() + w->stop() + Global_eps < d1)
				w_survive = false;

			d1 = w->d() + sqrt((w->start() - w->pseudo_x()) * (w->start() - w->pseudo_x()) + w->pseudo_y() * w->pseudo_y());
			if (v2->geodesic_distance() + e->length() - w->start() + Global_eps < d1)
				w_survive = false;

			if (w->compute_min_distance() > w->distance_bound())w_survive = false;

			iter_t = iter;
			iter = iter->next();
			if (!w_survive)
			{
				list->erase(iter_t);
				delete iter_t;
				--m_windows_wavefront;
			}
		}
	}

	//implement the discussed 6 cases as follow for simplicity. incompletely 6 cases
	inline windows_state GeodesicAlgorithmParallelAVTPDGG::check_between_two_windows(interval_pointer &w1, interval_pointer &w2, Triangle &tri)
	{
		double NUMERCIAL_EPSILON = 1 - 1e-12;
		// we implement the discussed 6 cases as follows for simplicity

		if ((w1->start() >= w2->start()) && (w1->start() <= w2->stop())) // w1->start
		{
			double Intersect_X, Intersect_Y;

			// judge the order of the two windows
			CalculateIntersectionPoint(w2->pseudo_x(), w2->pseudo_y(), w1->start(), 0, w1->pseudo_x(), w1->pseudo_y(), w1->stop(), 0, Intersect_X, Intersect_Y);

			if ((Intersect_Y < 0) && (Intersect_Y >= w1->pseudo_y()) && (Intersect_Y >= w2->pseudo_y()))
			{
				double d1, d2;
				d1 = w1->d() + sqrt((w1->start() - w1->pseudo_x()) * (w1->start() - w1->pseudo_x()) + (w1->pseudo_y()) * (w1->pseudo_y()));
				d2 = w2->d() + sqrt((w1->start() - w2->pseudo_x()) * (w1->start() - w2->pseudo_x()) + (w2->pseudo_y()) * (w2->pseudo_y()));

				if (d2 < d1 * NUMERCIAL_EPSILON)
					return w1_invalid;
				if (d1 < d2 * NUMERCIAL_EPSILON)
					w2->start() = w1->start();
			}
		}

		if ((w1->stop() >= w2->start()) && (w1->stop() <= w2->stop())) // w1->stop
		{
			double Intersect_X, Intersect_Y;

			// judge the order of the two windows
			CalculateIntersectionPoint(w2->pseudo_x(), w2->pseudo_y(), w1->stop(), 0, w1->pseudo_x(), w1->pseudo_y(), w1->start(), 0, Intersect_X, Intersect_Y);

			if ((Intersect_Y < 0) && (Intersect_Y >= w1->pseudo_y()) && (Intersect_Y >= w2->pseudo_y()))
			{
				double d1, d2;
				d1 = w1->d() + sqrt((w1->stop() - w1->pseudo_x()) * (w1->stop() - w1->pseudo_x()) + (w1->pseudo_y()) * (w1->pseudo_y()));
				d2 = w2->d() + sqrt((w1->stop() - w2->pseudo_x()) * (w1->stop() - w2->pseudo_x()) + (w2->pseudo_y()) * (w2->pseudo_y()));

				if (d2 < d1 * NUMERCIAL_EPSILON)
					return w1_invalid;
				if (d1 < d2 * NUMERCIAL_EPSILON)
					w2->stop() = w1->stop();
			}
		}

		if ((w2->start() >= w1->start()) && (w2->start() <= w1->stop())) // w2->start
		{
			double Intersect_X, Intersect_Y;

			// judge the previous order of the two windows
			CalculateIntersectionPoint(w1->pseudo_x(), w1->pseudo_y(), w2->start(), 0, w2->pseudo_x(), w2->pseudo_y(), w2->stop(), 0, Intersect_X, Intersect_Y);

			if ((Intersect_Y < 0) && (Intersect_Y >= w1->pseudo_y()) && (Intersect_Y >= w2->pseudo_y()))
			{
				double d1, d2;
				d1 = w1->d() + sqrt((w2->start() - w1->pseudo_x()) * (w2->start() - w1->pseudo_x()) + (w1->pseudo_y()) * (w1->pseudo_y()));
				d2 = w2->d() + sqrt((w2->start() - w2->pseudo_x()) * (w2->start() - w2->pseudo_x()) + (w2->pseudo_y()) * (w2->pseudo_y()));

				if (d1 < d2 * NUMERCIAL_EPSILON)
					return w2_invalid;
				if (d2 < d1 * NUMERCIAL_EPSILON)
					w1->start() = w2->start();
			}
		}

		if ((w2->stop() >= w1->start()) && (w2->stop() <= w1->stop())) // w2->stop
		{
			double Intersect_X, Intersect_Y;

			// judge the previous order of the two windows
			CalculateIntersectionPoint(w1->pseudo_x(), w1->pseudo_y(), w2->stop(), 0, w2->pseudo_x(), w2->pseudo_y(), w2->start(), 0, Intersect_X, Intersect_Y);

			if ((Intersect_Y < 0) && (Intersect_Y >= w1->pseudo_y()) && (Intersect_Y >= w2->pseudo_y()))
			{
				double d1, d2;
				d1 = w1->d() + sqrt((w2->stop() - w1->pseudo_x()) * (w2->stop() - w1->pseudo_x()) + (w1->pseudo_y()) * (w1->pseudo_y()));
				d2 = w2->d() + sqrt((w2->stop() - w2->pseudo_x()) * (w2->stop() - w2->pseudo_x()) + (w2->pseudo_y()) * (w2->pseudo_y()));

				if (d1 < d2 * NUMERCIAL_EPSILON)
					return w2_invalid;
				if (d2 < d1 * NUMERCIAL_EPSILON)
					w1->stop() = w2->stop();
			}
		}

		if (w1->start() >= w2->stop())
		{
			double Intersect_X, Intersect_Y;

			// judge the previous order of the two windows
			CalculateIntersectionPoint(w1->pseudo_x(), w1->pseudo_y(), w1->start(), 0, w2->pseudo_x(), w2->pseudo_y(), w2->stop(), 0, Intersect_X, Intersect_Y);

			face_pointer f = tri.bottom_edge->opposite_face(tri.face);
			edge_pointer e = f->next_edge(tri.bottom_edge, tri.left_vertex);
			double angle = f->vertex_angle(tri.left_vertex);
			double Cx = e->length() * cos(angle);
			double Cy = e->length() * -sin(angle);

			if ((PointInTriangle(Intersect_X, Intersect_Y, tri.bottom_edge->length(), Cx, Cy))
				&& (Intersect_Y < 0) && (Intersect_Y >= w1->pseudo_y()) && (Intersect_Y >= w2->pseudo_y()))
			{
				double d1, d2;
				d1 = w1->d() + sqrt((Intersect_X - w1->pseudo_x()) * (Intersect_X - w1->pseudo_x()) + (Intersect_Y - w1->pseudo_y()) * (Intersect_Y - w1->pseudo_y()));
				d2 = w2->d() + sqrt((Intersect_X - w2->pseudo_x()) * (Intersect_X - w2->pseudo_x()) + (Intersect_Y - w2->pseudo_y()) * (Intersect_Y - w2->pseudo_y()));

				if (d1 < d2 * NUMERCIAL_EPSILON)
					return w2_invalid;
				if (d2 < d1 * NUMERCIAL_EPSILON)
					return w1_invalid;
			}
		}

		if (w2->start() >= w1->stop())
		{
			double Intersect_X, Intersect_Y;

			// judge the previous order of the two windows
			CalculateIntersectionPoint(w2->pseudo_x(), w2->pseudo_y(), w2->start(), 0, w1->pseudo_x(), w1->pseudo_y(), w1->stop(), 0, Intersect_X, Intersect_Y);

			face_pointer f = tri.bottom_edge->opposite_face(tri.face);
			edge_pointer e = f->next_edge(tri.bottom_edge, tri.left_vertex);
			double angle = f->vertex_angle(tri.left_vertex);
			double Cx = e->length() * cos(angle);
			double Cy = e->length() * -sin(angle);

			if ((PointInTriangle(Intersect_X, Intersect_Y, tri.bottom_edge->length(), Cx, Cy))
				&& (Intersect_Y < 0) && (Intersect_Y >= w1->pseudo_y()) && (Intersect_Y >= w2->pseudo_y()))
			{
				double d1, d2;
				d1 = w1->d() + sqrt((Intersect_X - w1->pseudo_x()) * (Intersect_X - w1->pseudo_x()) + (Intersect_Y - w1->pseudo_y()) * (Intersect_Y - w1->pseudo_y()));
				d2 = w2->d() + sqrt((Intersect_X - w2->pseudo_x()) * (Intersect_X - w2->pseudo_x()) + (Intersect_Y - w2->pseudo_y()) * (Intersect_Y - w2->pseudo_y()));

				if (d1 < d2 - NUMERCIAL_EPSILON)
					return w2_invalid;
				if (d2 < d1 - NUMERCIAL_EPSILON)
					return w1_invalid;
			}
		}

		return both_valid;
	}

	inline void GeodesicAlgorithmParallelAVTPDGG::pairwise_windows_checking(list_pointer &list, Triangle &tri)
	{
		if (list->empty()) return;

		interval_pointer iter = list->begin();
		interval_pointer next, iter_t;
		next = iter->next();

		while ((!list->empty()) && (next != NULL))
		{
			windows_state ws = check_between_two_windows(iter, next, tri);
			switch (ws)
			{
			case geodesic::w1_invalid:
				iter_t = iter;
				if (iter == list->begin())
				{
					iter = iter->next();
				}
				else
				{
					iter = iter->previous();
				}
				list->erase(iter_t);
				delete iter_t;
				--m_windows_wavefront;
				break;

			case geodesic::w2_invalid:
				list->erase(next);
				delete next;
				--m_windows_wavefront;
				break;
			case geodesic::both_valid:
				iter = iter->next();
				break;
			default:
				break;
			}
			next = iter->next();
		}
	}

	//Rule 1
	inline void GeodesicAlgorithmParallelAVTPDGG::find_separating_point(list_pointer &list, Triangle &tri, Vertex &top_t)
	{
		const double LOCAL_EPSILON = 1e-20 * list->edge()->length(); // numerical issue

		double L = tri.left_edge->length();
		double top_x = L * cos(tri.left_alpha);
		double top_y = L * sin(tri.left_alpha);

		memcpy(&top_t, tri.top_vertex, sizeof(Vertex));
		
		top_t.geodesic_distance() = GEODESIC_INF;

		interval_pointer iter = list->begin();

		double wlist_sp = 0;
		double wlist_pseudo_x = 0;
		double wlist_pseudo_y = 0;
		bool updated = false;
		while (iter != NULL)
		{
			interval_pointer &w = iter;


			double w_sp = w->pseudo_x() - w->pseudo_y() * ((top_x - w->pseudo_x()) / (top_y - w->pseudo_y()));
			double distance = GEODESIC_INF;

			// shortest path from the window
			if ((w_sp - w->start() > LOCAL_EPSILON) && (w_sp - w->stop() < -LOCAL_EPSILON))
			{
				distance = w->d() + sqrt((top_x - w->pseudo_x()) * (top_x - w->pseudo_x()) + (top_y - w->pseudo_y()) * (top_y - w->pseudo_y()));
				//w->shortest_distance() = distance;
			}
			else if (w_sp - w->start() <= LOCAL_EPSILON)
			{
				distance = w->d() + sqrt((top_x - w->start()) * (top_x - w->start()) + top_y * top_y) + sqrt((w->start() - w->pseudo_x()) * (w->start() - w->pseudo_x()) + w->pseudo_y() * w->pseudo_y());
				//w->shortest_distance() = distance;
				w_sp = w->start();
			}
			else if (w_sp - w->stop() >= -LOCAL_EPSILON)
			{
				distance = w->d() + sqrt((top_x - w->stop()) * (top_x - w->stop()) + top_y * top_y) + sqrt((w->stop() - w->pseudo_x()) * (w->stop() - w->pseudo_x()) + w->pseudo_y() * w->pseudo_y());
				//w->shortest_distance() = distance;
				w_sp = w->stop();
			}


			// update information at top_t
			if (distance < top_t.geodesic_distance())
			{
				top_t.geodesic_distance() = distance;
				top_t.incident_face() = tri.face;
				//top_t.incident_point() = (list->start_vertex() == list->edge()->v0()) ? w_sp : list->edge()->length() - w_sp;
				wlist_sp = w_sp;
				wlist_pseudo_x = w->pseudo_x();
				wlist_pseudo_y = w->pseudo_y();
				updated = true;
			}
			w->sp() = w_sp;

			iter = iter->next();
		}
		if (updated){
			double middle_edge_len = sqrtf((top_x - wlist_sp)*(top_x - wlist_sp) + top_y*top_y);
			double p_start_angle = acos((L*L + middle_edge_len*middle_edge_len - wlist_sp*wlist_sp) / (2 * L*middle_edge_len));
			top_t.p_start_angle() = p_start_angle;
			top_t.edge_start_vertex() = tri.left_edge;
		}

		list->sp() = wlist_sp;
		list->pseudo_x() = wlist_pseudo_x;
		list->pseudo_y() = wlist_pseudo_y;
	}

	inline void GeodesicAlgorithmParallelAVTPDGG::propagate_windows_to_two_edges(list_pointer &list, Triangle &tri, IntervalList &wl_left, IntervalList &wl_right)
	{
		const double LOCAL_EPSILON = 1e-8 * list->edge()->length(); // numerical issue

		interval_pointer iter = list->begin();
		interval_pointer iter_t;

		enum PropagationDirection
		{
			LEFT,
			RIGHT,
			BOTH
		};

		PropagationDirection direction;

		while (!list->empty() && (iter != NULL))
		{
			interval_pointer &w = iter;
			assert(w->start() <= w->stop());


			if (w->sp() < list->sp() - LOCAL_EPSILON)
			{
				// only propagate to left edge
				double Intersect_X, Intersect_Y;

				// judge the positions of the two windows
				CalculateIntersectionPoint(list->pseudo_x(), list->pseudo_y(), list->sp(), 0, w->pseudo_x(), w->pseudo_y(), w->stop(), 0, Intersect_X, Intersect_Y);
				if ((w->stop() < list->sp()) || ((Intersect_Y <= 0) && (Intersect_Y >= list->pseudo_y()) && (Intersect_Y >= w->pseudo_y())))
				{
					direction = LEFT;
				}
				
				else direction = BOTH;
			}
			else if (w->sp() > list->sp() + LOCAL_EPSILON)
			{
				// only propagate to right edge
				double Intersect_X, Intersect_Y;

				// judge the positions of the two windows
				CalculateIntersectionPoint(list->pseudo_x(), list->pseudo_y(), list->sp(), 0, w->pseudo_x(), w->pseudo_y(), w->start(), 0, Intersect_X, Intersect_Y);
				if ((w->start() > list->sp()) || ((Intersect_Y <= 0) && (Intersect_Y >= list->pseudo_y()) && (Intersect_Y >= w->pseudo_y())))
				{
					direction = RIGHT;
				}
				
				else direction = BOTH;
			}
			else
			{
				// propagate to both edges
				direction = BOTH;
			}
			
			
			bool ValidPropagation;
			interval_pointer right_w;

			switch (direction) {
			case LEFT:
				ValidPropagation = compute_propagated_parameters(w->pseudo_x(),
					w->pseudo_y(),
					w->start(),
					w->stop(),
					tri.left_alpha,
					tri.left_edge->length(),
					w,
					w->d());

				iter_t = iter->next();
				if (ValidPropagation)
				{
					if (w->stop() == tri.left_edge->length()){
						w->right_CH_distance() = sqrtf((w->stop() - w->pseudo_x())*(w->stop() - w->pseudo_x()) + w->pseudo_y()*w->pseudo_y());
						//w->right_cosine_CH_angle() = 1;
						// update distance bound here..
						updateDistanceBound(w);
					}
					list->erase(w);
					wl_left.push_back(w);
				}
				else
				{
					list->erase(w);
					delete w;
					--m_windows_wavefront;
				}
				iter = iter_t;
				break;

			case RIGHT:
				ValidPropagation = compute_propagated_parameters(tri.bottom_edge->length() - w->pseudo_x(),
					w->pseudo_y(),
					tri.bottom_edge->length() - w->stop(),
					tri.bottom_edge->length() - w->start(),
					tri.right_alpha,
					tri.right_edge->length(),
					w,
					w->d());

				iter_t = iter->next();
				if (ValidPropagation)
				{
					double length = tri.right_edge->length(); // invert window
					double start = length - w->stop();
					w->stop() = length - w->start();
					w->start() = start;
					w->pseudo_x() = length - w->pseudo_x();

					if (w->start() == 0){
						w->left_CH_distance() = sqrtf((0 - w->pseudo_x())*(0 - w->pseudo_x()) + w->pseudo_y()*w->pseudo_y());
						//w->left_cosine_CH_angle() = 1;
						// update distance bound here..
						updateDistanceBound(w);
					}

					list->erase(w);
					wl_right.push_back(w);
				}
				else
				{
					list->erase(w);
					delete w;
					--m_windows_wavefront;
				}
				iter = iter_t;
				break;

			case BOTH:
				right_w = new Interval;
				memcpy(right_w, w, sizeof(Interval));

				ValidPropagation = compute_propagated_parameters(w->pseudo_x(),
					w->pseudo_y(),
					w->start(),
					w->stop(),
					tri.face->vertex_angle(tri.left_vertex),
					tri.left_edge->length(),
					w,
					w->d());

				iter_t = iter->next();
				if (ValidPropagation)
				{
					if (w->stop() == tri.left_edge->length()){
						w->right_CH_distance() = sqrtf((w->stop() - w->pseudo_x())*(w->stop() - w->pseudo_x()) + w->pseudo_y()*w->pseudo_y());
						//	w->right_cosine_CH_angle() = 1;
						// update distance bound here..
						updateDistanceBound(w);
					}
					list->erase(w);
					wl_left.push_back(w);
				}
				else
				{
					list->erase(w);
					delete w;
					--m_windows_wavefront;
				}
				iter = iter_t;

				ValidPropagation = compute_propagated_parameters(tri.bottom_edge->length() - right_w->pseudo_x(),
					right_w->pseudo_y(),
					tri.bottom_edge->length() - right_w->stop(),
					tri.bottom_edge->length() - right_w->start(),
					tri.face->vertex_angle(tri.right_vertex),
					tri.right_edge->length(),
					right_w,
					right_w->d());

				if (ValidPropagation)
				{
					// invert window
					double length = tri.right_edge->length();
					double start = length - right_w->stop();
					right_w->stop() = length - right_w->start();
					right_w->start() = start;
					right_w->pseudo_x() = length - right_w->pseudo_x();

					if (right_w->start() == 0){
						right_w->left_CH_distance() = sqrtf((0 - right_w->pseudo_x())*(0 - right_w->pseudo_x()) + right_w->pseudo_y()*right_w->pseudo_y());
						//	right_w->left_cosine_CH_angle() = 1;
						// update distance bound here..
						updateDistanceBound(right_w);
					}

					wl_right.push_back(right_w);

					++m_windows_wavefront;
					if (m_windows_peak < m_windows_wavefront)
						m_windows_peak = m_windows_wavefront;
				}
				else
				{
					delete right_w;
				}
				break;

			default:
				break;
			}
		}
	}

	inline void GeodesicAlgorithmParallelAVTPDGG::propagate(unsigned source)
	{
		for (int i = 0; i < mesh()->vertices().size(); i++){
			Vertex& v = mesh()->vertices()[i];
			vertex_pointer p_v = &v;
			double sum_angle = 0;
			//std::vector<double> temp;
			for (int j = 0; j<v.adjacent_faces().size(); j++){
				sum_angle += v.adjacent_faces()[j]->vertex_angle(p_v);

			}
			//sum_vertex_angles.push_back(temp);
			v.sum_angle_around_v() = sum_angle;

		}

		
		printf("num_procs: %d\n", num_procs);
		m_source = source;
		//printf("m_windows_wavefront : %d\n", m_windows_wavefront);
		//printf("m_windows_peak : %d\n", m_windows_peak);
		initialize_propagation_data();
		//printf("after init, m_windows_wavefront : %d\n", m_windows_wavefront);
		//printf("m_windows_peak : %d\n", m_windows_peak);
		num_propagated_vertices = 0;

		prevPhase = 1;
		m_iterations = 0;
		double proportion = 0;
		from_list = 0;
		double push_list_time = 0;
		list_concurrent_time = 0;
		//double st = omp_get_wtime();
		tbb::tick_count st = tbb::tick_count::now();
		double push_vertex_time = 0;
		int count_push_list = 0;
		while (!m_vertex_FWP_queue.empty())
		{
			
			if (prevPhase > 1)
			{
				prevPhase++;
				m_vertex_FWP_queue.updateParamenters();
				
			}
			concurrent_num = 0;
			vertex_pointer vert = NULL;
			clock_t start = clock();
			// pop suitable num of window lists
			int size = m_vertex_FWP_queue.size();
			int count = 0;
			vertextoprop.clear();
			do{
				//1) if too many window lists, break
				if (list_queue[from_list].size() > mx_concurrent_list_num)
				{
					break;
				}

				m_vertex_FWP_queue.updatePhase();
			
				//2) if pop the vertex at next wavefront, break
				if (prevPhase < m_vertex_FWP_queue.getPhase())
				{
					prevPhase++;
					//printf("Proportion : %lf\n", count/(double)size);
					proportion += count / (double)size;
					break;
				}
				
				vert = m_vertex_FWP_queue.pop_for_parallel();
				count++;
				vert->ptrInQueue = NULL;
				
				concurrent_num++;
				
				// for statistics
				num_propagated_vertices++;

				// update wavefront: update adjacent vertices' state into 'FRONT'
				vert->state() = Vertex::INSIDE;
				/*for (unsigned i = 0; i < vert->adjacent_edges().size(); i++)
				{
					vertex_pointer vert_it = vert->adjacent_edges()[i]->opposite_vertex(vert);
					if (vert_it->state() == Vertex::OUTSIDE) vert_it->state() = Vertex::FRONT;
				}*/

				//handle saddle vertex
				//if (vert->saddle_or_boundary()) create_pseudo_source_windows(vert, false, from_list);
				updateDistance(vert);
				vertextoprop.push_back(vert);
				//printf("after flat_v_fan, m_windows_wavefront : %d\n", m_windows_wavefront);
				//printf("m_windows_peak : %d\n", m_windows_peak);
				//push window lists on the wavefront incident to 'vert' into 'list_queue[from_list]'
				clock_t start_l = clock();
				for (unsigned i = 0; i < vert->adjacent_edges().size(); i++)
				{
					edge_pointer edge_it = vert->adjacent_edges()[i];
					if (!interval_list_0(edge_it)->empty() && !interval_list_0(edge_it)->readyPropagate())
					{
						list_queue[from_list].push_back(interval_list_0(edge_it));
						interval_list_0(edge_it)->readyPropagate() = true;
						count_push_list++;
					}
					if (!interval_list_1(edge_it)->empty() && !interval_list_1(edge_it)->readyPropagate())
					{
						list_queue[from_list].push_back(interval_list_1(edge_it));
						interval_list_1(edge_it)->readyPropagate() = true;
						count_push_list++;
					}
				}
				for (unsigned i = 0; i < vert->adjacent_faces().size(); i++)
				{
					edge_pointer edge_it = vert->adjacent_faces()[i]->opposite_edge(vert);
					vertex_pointer vert_0 = edge_it->v0();
					vertex_pointer vert_1 = edge_it->v1();
					if ((edge_it->adjacent_faces().size() < 2) || (vert_0->state() == Vertex::INSIDE) || (vert_1->state() == Vertex::INSIDE))
					{
						if (!interval_list_0(edge_it)->empty() && !interval_list_0(edge_it)->readyPropagate())
						{
							list_queue[from_list].push_back(interval_list_0(edge_it));
							interval_list_0(edge_it)->readyPropagate() = true;
							count_push_list++;
						}
						if (!interval_list_1(edge_it)->empty() && !interval_list_1(edge_it)->readyPropagate())
						{
							list_queue[from_list].push_back(interval_list_1(edge_it));
							interval_list_1(edge_it)->readyPropagate() = true;
							count_push_list++;
						}
					}
				}
				clock_t stop_l = clock();
				push_list_time += static_cast<double>(stop_l - start_l);
			} while (!m_vertex_FWP_queue.empty());

			clock_t stop = clock();
			clock_neededoperations += static_cast<double>(stop - start);
			if (concurrent_num < 1) continue;
			tbb::tick_count st_concurrent_list = tbb::tick_count::now();
			tbb::parallel_for(tbb::blocked_range<size_t>(0, vertextoprop.size()), [&](const tbb::blocked_range<size_t>& r)
			{
				for (size_t i = r.begin(); i != r.end(); i++)
				{
					//create_flat_v_fan_shape(vertextoprop[i], false);
					create_pseudo_source_windows(vertextoprop[i], false);
				}
			});
			tbb::tick_count nd_concurrent_list = tbb::tick_count::now();
			list_concurrent_time += (nd_concurrent_list - st_concurrent_list).seconds();

			//for propagate lists concurrently, update vertices and merge lists sequentially
			
			int iter = 0;
			while (!list_queue[from_list].empty())
			{
				iter++;
				if (estimate_wavefront < list_queue[from_list].size())
				{
					estimate_wavefront = list_queue[from_list].size();
					tris.resize(estimate_wavefront);
					wl_lefts.resize(estimate_wavefront);
					wl_rights.resize(estimate_wavefront);
					top_ts.resize(estimate_wavefront);
					isNotBoundary.resize(estimate_wavefront, 0);
				}
				tbb::tick_count st_concurrent_list = tbb::tick_count::now();
				tbb::parallel_for(tbb::blocked_range<size_t>(0, list_queue[from_list].size(), 5), [&](const tbb::blocked_range<size_t>& r)
				{
					for (size_t i = r.begin(); i != r.end(); i++)
					{
						list_pointer list = list_queue[from_list][i];
						bool is_boundary = calculate_triangle_parameters(list, tris[i]);
						if (!is_boundary)
						{
							isNotBoundary[i] = 1;
							check_with_vertices(list);
							//pairwise_windows_checking(list, tris[i]);

							find_separating_point(list, tris[i], top_ts[i]);
							wl_lefts[i].clear(); wl_rights[i].clear();
							propagate_windows_to_two_edges(list, tris[i], wl_lefts[i], wl_rights[i]);
							list->readyPropagate() = false;
						}
					}
				});
				tbb::tick_count nd_concurrent_list = tbb::tick_count::now();
				list_concurrent_time += (nd_concurrent_list - st_concurrent_list).seconds();

				to_list = from_list ^ 1;
				list_queue[to_list].clear();
				//printf("after prop windows to 2 edges, m_windows_wavefront : %d\n", m_windows_wavefront);
				//printf("m_windows_peak : %d\n", m_windows_peak);
				//sequentially: update vertex and merge lists
				for (unsigned i = 0; i < list_queue[from_list].size(); i++)
				{
					if (isNotBoundary[i])
					{
						isNotBoundary[i] = 0;
						//update vertex
						clock_t start_q = clock();
						if (top_ts[i].geodesic_distance() < tris[i].top_vertex->geodesic_distance())
						{

							//if (tris[i].top_vertex->state() == Vertex::FRONT) //
								erase_from_FWP_queue(tris[i].top_vertex);
								
							memcpy(tris[i].top_vertex, &top_ts[i], sizeof(Vertex));
							//clock_t stop_push_v = clock();
							//push_vertex_time += static_cast<double>(stop_push_v - start_push_v);
							//if (tris[i].top_vertex->state() == Vertex::FRONT)//
								tris[i].top_vertex->ptrInQueue = m_vertex_FWP_queue.push(tris[i].top_vertex, tris[i].top_vertex->geodesic_distance());

							
						}
						clock_t stop_q = clock();
						clock_neededoperations += static_cast<double>(stop_q - start_q);
						push_vertex_time += static_cast<double>(stop_q - start_q);
						//merge windows lists
						clock_t start_l = clock();
						if (!wl_lefts[i].empty())
						{
							if (!tris[i].left_list->empty())
							{
								tris[i].left_list->begin()->previous() = wl_lefts[i].end();
								wl_lefts[i].end()->next() = tris[i].left_list->begin();
								tris[i].left_list->begin() = wl_lefts[i].begin();
							}
							else
							{
								tris[i].left_list->begin() = wl_lefts[i].begin();
								tris[i].left_list->end() = wl_lefts[i].end();
							}

							if ((!tris[i].left_list->readyPropagate())
								&& ((tris[i].left_edge->v0()->state() == Vertex::INSIDE) || (tris[i].left_edge->v1()->state() == Vertex::INSIDE)) && (!tris[i].left_list->empty()))
							{
								list_queue[to_list].push_back(tris[i].left_list);
								tris[i].left_list->readyPropagate() = true;
								count_push_list++;
							}
							else if ((!tris[i].left_list->empty()) && (!tris[i].left_list->readyPrinted()))
							{
								tris[i].left_list->readyPrinted() = true;
								count_push_list++;
							}
						}

						if (!wl_rights[i].empty())
						{
							if (!tris[i].right_list->empty())
							{
								tris[i].right_list->end()->next() = wl_rights[i].begin();
								wl_rights[i].begin()->previous() = tris[i].right_list->end();
								tris[i].right_list->end() = wl_rights[i].end();
							}
							else
							{
								tris[i].right_list->begin() = wl_rights[i].begin();
								tris[i].right_list->end() = wl_rights[i].end();
							}

							if ((!tris[i].right_list->readyPropagate())
								&& ((tris[i].right_edge->v0()->state() == Vertex::INSIDE) || (tris[i].right_edge->v1()->state() == Vertex::INSIDE)) && (!tris[i].right_list->empty()))
							{
								list_queue[to_list].push_back(tris[i].right_list);
								tris[i].right_list->readyPropagate() = true;
								count_push_list++;
							}
							else if ((!tris[i].right_list->empty()) && (!tris[i].right_list->readyPrinted()))
							{
								tris[i].right_list->readyPrinted() = true;
								count_push_list++;
							}
						}
						clock_t stop_l = clock();
						push_list_time += static_cast<double>(stop_l - start_l);
					}

				}
				
				from_list = to_list;
			}
			//printf("No Iterations : %d\n", iter);
			m_iterations++;
		}
		
		//double nd = omp_get_wtime();
		tbb::tick_count nd = tbb::tick_count::now();
		m_time_consumed = (nd - st).seconds();
		proportion /= m_iterations;
		printf("Push list time : %lf\n",push_list_time/CLOCKS_PER_SEC);
		printf("List pushed : %d\n",count_push_list);
		printf("Proportions average vertex propagated : %lf\n", proportion);
		printf("push vertex time; %lf\n", push_vertex_time / CLOCKS_PER_SEC);
	}

	inline void GeodesicAlgorithmParallelAVTPDGG::print_statistics()
	{
		printf("\n-------------------------- for statistics --------------------------\n");
		
		printf("Time in total: %.6lf\n", m_time_consumed);
		printf("Processing list concurrently cost %.6lf\n", list_concurrent_time);
		printf("Time for needed sequentials : %lf\n", clock_neededoperations/CLOCKS_PER_SEC);
		double memory = sizeof(Interval);
		double used_memory = memory * m_windows_peak / 1e6;
		printf("Time for fan-shape : %lf\n",clockfanshape/CLOCKS_PER_SEC);
		printf("Num of propagated vertices: %d\n", num_propagated_vertices);
		printf("Num of propagated windows: %lld\n", m_windows_propagation);
		printf("Peak number of windows on wavefront: %d\n", m_windows_peak);
		printf("Used memory: %.6lf\n", used_memory);
		printf("Max num of window list: %d\n", estimate_wavefront);
		printf("Iterations: %d\n\n", m_iterations);

		//output data to .csv(excel file)
		FILE* file = fopen("Parallel_FWP_VTP.csv", "a");
		fprintf(file, "%.6lf, %.6lf, %d, %lld, %d,%.6lf,%d, %d, %d\n", m_time_consumed, list_concurrent_time, num_propagated_vertices, m_windows_propagation,
			m_windows_peak, used_memory, m_iterations, mx_concurrent_list_num, estimate_wavefront);
		fclose(file);
	}

}

#endif