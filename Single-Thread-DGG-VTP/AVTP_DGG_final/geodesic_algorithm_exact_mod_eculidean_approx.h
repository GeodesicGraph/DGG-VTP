#ifndef GEODESIC_ALGORITHM_EXACT
#define GEODESIC_ALGORITHM_EXACT

#include "stdafx.h"
#include "geodesic_memory.h"
#include "geodesic_algorithm_base.h"
#include "geodesic_algorithm_exact_elements.h"
//#include <unordered_set>
#pragma warning(disable : 4996)
namespace geodesic {

	class GeodesicAlgorithmExact : public GeodesicAlgorithmBase
	{
	public:

		// basic functions related to class
		GeodesicAlgorithmExact(geodesic::Mesh* mesh, double eps, double tol, int sidefan_const) :
			GeodesicAlgorithmBase(mesh)
		//	m_memory_allocator(mesh->edges().size(), mesh->edges().size())
		{
			epsilon = eps;
			tolerance = tol;
			side_fan_angle = sidefan_const*(sqrt(epsilon));
		};
		~GeodesicAlgorithmExact() {};
		void clear() {
			GeodesicAlgorithmBase::clear();
		}
			// m_memory_allocator.clear(); };

		// main entry
		void propagate(unsigned source);

		// print the resulting statistics
		void print_statistics();

	private:

		// simple functions
		void updateDistanceBound(interval_pointer &w); // update the distance bound for a window based on FastDGG Rules..
		void initialize_propagation_data();
		double compute_cosine_vector(double v1_x, double v1_y, double v2_x, double v2_y);
		void create_pseudo_source_windows(vertex_pointer &v, bool UpdateFIFOQueue);
		void propagate_from_pseudo_source(vertex_pointer pseudo_source, bool inside_traversed_area, face_pointer to_face);
		void create_pseudo_source_fan_shape(vertex_pointer &v, bool UpdateFIFOQueue); // propagate from saddle vertex - exact algo

		void create_saddle_fan_shape(vertex_pointer &saddle, bool UpdateFIFOQueue); // propagate fan from saddle vertex with fan - approx algo
		void create_flat_v_fan_shape(vertex_pointer &flat_v, bool UpdateFIFOQueue); // propagate fan from spehrical / euclidena vertex (i.e. "flat" when the faces around the vertex are unfolded)

		void erase_from_queue(vertex_pointer v);

		// propagate a windows list (Rule 1)
		void find_separating_point(list_pointer &list); // find the separating point of the windows and the list
		void propagate_windows_to_two_edges(list_pointer &list); // propagates windows to two edges accross a triangle face

		// pairwise windows checking (Rule 2)
		void check_with_vertices(list_pointer &list);
		windows_state check_between_two_windows(interval_pointer &w1, interval_pointer &w2); // Check two neighbouring crossing windows on same edge
		void pairwise_windows_checking(list_pointer &list); // Check crossing windows on same edge
		double clockcountbound;
		double clockpseudoprop;
		// main operation
		void propagate_one_windows_list(list_pointer &list);

		// member variables
		//std::multiset<vertex_pointer, Vertex> m_vertex_queue; // priority queue for vertices
		std::set<vertex_pointer, Vertex> m_vertex_queue;
		std::queue<list_pointer> m_list_queue;                // FIFO queue for lists
		MemoryAllocator<Interval> m_memory_allocator;		  // quickly allocate and deallocate intervals 
		std::vector<std::vector<double>> sum_vertex_angles;
		std::vector<std::vector<face_pointer>>face_pointers_around_vertex;
		double side_fan_angle; // the fan circular sector angle size located on the outside of the left & right of original fan region.
		double epsilon;
		double tolerance;
		unsigned m_source;
	};

	//----------------- simple functions ---------------------
	inline void GeodesicAlgorithmExact::initialize_propagation_data()
	{
		clear();

		// initialize source's parameters
		vertex_pointer source = &(this->mesh()->vertices()[m_source]);
		source->geodesic_distance() = 0;
		source->state() = Vertex::INSIDE;

		// initialize windows around source
		create_pseudo_source_windows(source, false);
	}

	inline void GeodesicAlgorithmExact::erase_from_queue(vertex_pointer v)
	{
		assert(m_vertex_queue.count(v) <= 1);

		std::multiset<vertex_pointer, Vertex>::iterator it = m_vertex_queue.find(v);
		if (it != m_vertex_queue.end())
			m_vertex_queue.erase(it);
	}

	inline double GeodesicAlgorithmExact::compute_cosine_vector(double v1_x, double v1_y, double v2_x, double v2_y){

		double dotproduct = v1_x*v2_x + v1_y*v2_y;
		double squaredlengthsproduct = (v1_x*v1_x + v1_y*v1_y)*(v2_x*v2_x + v2_y*v2_y);
		return dotproduct / sqrtf(squaredlengthsproduct);
	}

	inline void GeodesicAlgorithmExact::propagate_from_pseudo_source(vertex_pointer pseudo_source, bool inside_traversed_area, face_pointer to_face){

		m_windows_pseudo++;
		clock_t start_clock = clock();
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
		updateDistanceBound(candidate);
		// insert into list
		list->push_back(candidate);

		// push into M_LIST_QUEUE if inside traversed area
		if ((inside_traversed_area) && ((edge_it->v0()->state() != Vertex::FRONT) || (edge_it->v1()->state() != Vertex::FRONT)))
			m_list_queue.push(list);

		// Statistics
		++m_windows_wavefront;
		if (m_windows_peak < m_windows_wavefront)
			m_windows_peak = m_windows_wavefront;
		clock_t stop_clock = clock();
		//clockpseudoprop += static_cast<double>(stop_clock - start_clock);
	}

	inline void GeodesicAlgorithmExact::create_flat_v_fan_shape(vertex_pointer &pseudo_source, bool inside_traversed_area){

		clock_t start_fan_clock = clock();

		double totalanglearound = pseudo_source->sum_angle_around_v();
		if (totalanglearound <= (0.9*M_PI)){
			for (unsigned i = 0; i < pseudo_source->adjacent_edges().size(); ++i)
			{
				edge_pointer   edge_it = pseudo_source->adjacent_edges()[i];
				vertex_pointer vert_it = edge_it->opposite_vertex(pseudo_source);

				double distance = pseudo_source->geodesic_distance() + edge_it->length();
				if (distance < vert_it->geodesic_distance())
				{
					//clock_t start_q = clock();
					m_vertex_queue.erase(vert_it);
					//clock_t stop_q = clock();
					//clockpseudoprop += static_cast<double>(stop_q - start_q);
					vert_it->geodesic_distance() = distance;
					if (vert_it->state() == Vertex::OUTSIDE)
						vert_it->state() = Vertex::FRONT;

					vert_it->incident_face() = edge_it->adjacent_faces()[0];
					edge_pointer next_edge = vert_it->incident_face()->next_edge(edge_it, pseudo_source);
					vert_it->incident_point() = (next_edge->v0() == pseudo_source) ? 0 : next_edge->length();
					vert_it->p_start_angle() = 0;
					vert_it->edge_start_vertex() = edge_it;
					//clock_t start = clock();
					m_vertex_queue.insert(vert_it);
					//clock_t stop = clock();
					//clockpseudoprop += static_cast<double>(stop - start);
				}
			}
			
			return; }
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

		assert(0 <= sum_angle_start_dir <= (1 + 1e-10)*incident_face->vertex_angle(pseudo_source));

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
		
		for (unsigned i = 0; i < pseudo_source->adjacent_edges().size(); ++i)
		{
			edge_pointer   edge_it = pseudo_source->adjacent_edges()[i];
			vertex_pointer vert_it = edge_it->opposite_vertex(pseudo_source);

			double distance = pseudo_source->geodesic_distance() + edge_it->length();
			if (distance < vert_it->geodesic_distance())
			{
				//clock_t start_q = clock();
				m_vertex_queue.erase(vert_it);
				//clock_t stop_q = clock();
				//clockpseudoprop += static_cast<double>(stop_q - start_q);
				vert_it->geodesic_distance() = distance;
				if (vert_it->state() == Vertex::OUTSIDE)
					vert_it->state() = Vertex::FRONT;

				vert_it->incident_face() = edge_it->adjacent_faces()[0];
				edge_pointer next_edge = vert_it->incident_face()->next_edge(edge_it, pseudo_source);
				vert_it->incident_point() = (next_edge->v0() == pseudo_source) ? 0 : next_edge->length();
				vert_it->p_start_angle() = 0;
				vert_it->edge_start_vertex() = edge_it;
				//clock_t start = clock();
				m_vertex_queue.insert(vert_it);
				//clock_t stop = clock();
				//clockpseudoprop += static_cast<double>(stop - start);
			}
		}
		
		
		// propagate from here....
		if (f_from_end_vertex->id() == f_from_start_vertex->id()){
			vertex_pointer start_vertex = edge_on_start_vertex->opposite_vertex(pseudo_source);
			vertex_pointer end_vertex = edge_on_end_vertex->opposite_vertex(pseudo_source);
			//double mid_start_angle = M_PI - angle_f_start - f_from_start_vertex->vertex_angle(start_vertex);
			//double start = edge_on_start_vertex->length()*sin(angle_f_start) / sin(mid_start_angle);
			//double mid_end_angle = M_PI - angle_f_end - f_from_start_vertex->vertex_angle(end_vertex);
			//double end = edge_on_end_vertex->length()*sin(angle_f_end) / sin(mid_end_angle);
			edge_pointer edge_opp = f_from_start_vertex->opposite_edge(pseudo_source);
			//end = edge_opp->length() - end;
			//assert((1 + 1e-10)*start <= end);

			propagate_from_pseudo_source(pseudo_source, inside_traversed_area, f_from_start_vertex);
		}
		else{
			vertex_pointer start_vertex = edge_on_start_vertex->opposite_vertex(pseudo_source);
			/*double mid_start_angle = M_PI - angle_f_start - f_from_start_vertex->vertex_angle(start_vertex);
			double start = edge_on_start_vertex->length()*sin(angle_f_start) / sin(mid_start_angle);*/
			edge_pointer edge_opp = f_from_start_vertex->opposite_edge(pseudo_source);
			double end = edge_opp->length();
//			assert((1 + 1e-10)*start <= end);
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
			end = f_from_end_vertex->opposite_edge(pseudo_source)->length();
			propagate_from_pseudo_source(pseudo_source, inside_traversed_area, f_from_end_vertex);

		}
		clock_t stop_fan_clock = clock();
		clockpseudoprop += static_cast<double>(stop_fan_clock - start_fan_clock);
	}

	inline void GeodesicAlgorithmExact::create_pseudo_source_fan_shape(vertex_pointer &pseudo_source, bool inside_traversed_area){

		/*double sum_vertex_angle = 0;
		for (int i = 0; i < pseudo_source->adjacent_faces().size(); i++){
		sum_vertex_angle += pseudo_source->adjacent_faces()[i]->vertex_angle(pseudo_source);
		}
		if (sum_vertex_angle <= (1+1e-9)*2 * M_PI)return;

		edge_pointer edge_on_start_vertex = pseudo_source->edge_start_vertex();
		face_pointer incident_face = pseudo_source->incident_face();
		edge_pointer edge_on_end_vertex = incident_face->next_edge(edge_on_start_vertex, pseudo_source);
		double p_start_len = edge_on_start_vertex->length();
		double sum_angle_start_dir = pseudo_source->p_start_angle();
		double sum_angle_end_dir = incident_face->vertex_angle(pseudo_source) - sum_angle_start_dir;

		assert(0 <= sum_angle_start_dir <= (1+1e-10)*incident_face->vertex_angle(pseudo_source));

		double sum_angle = sum_angle_start_dir;
		face_pointer incident_it = incident_face;

		while (sum_angle <= M_PI){
		incident_it = edge_on_start_vertex->opposite_face(incident_it);
		sum_angle += incident_it->vertex_angle(pseudo_source);
		edge_on_start_vertex = incident_it->next_edge(edge_on_start_vertex, pseudo_source);

		}
		edge_on_start_vertex = incident_it->next_edge(edge_on_start_vertex, pseudo_source);
		face_pointer f_from_start_vertex = incident_it;

		double angle_f_start = f_from_start_vertex->vertex_angle(pseudo_source) - (sum_angle - M_PI);



		incident_it = incident_face;
		sum_angle = sum_angle_end_dir;
		while (sum_angle <= M_PI){
		incident_it = edge_on_end_vertex->opposite_face(incident_it);
		sum_angle += incident_it->vertex_angle(pseudo_source);
		edge_on_end_vertex = incident_it->next_edge(edge_on_end_vertex, pseudo_source);

		}
		face_pointer f_from_end_vertex = incident_it;
		edge_on_end_vertex = incident_it->next_edge(edge_on_end_vertex, pseudo_source);
		double angle_f_end = f_from_end_vertex->vertex_angle(pseudo_source) - (sum_angle - M_PI);*/


		// update vertices around pseudo_source

		/*for (unsigned i = 0; i < pseudo_source->adjacent_edges().size(); ++i)
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
		}

		// propagate from here....
		if (f_from_end_vertex->id() == f_from_start_vertex->id()){
		vertex_pointer start_vertex = edge_on_start_vertex->opposite_vertex(pseudo_source);
		vertex_pointer end_vertex = edge_on_end_vertex->opposite_vertex(pseudo_source);
		double mid_start_angle = M_PI - angle_f_start - f_from_start_vertex->vertex_angle(start_vertex);
		double start = edge_on_start_vertex->length()*sin(angle_f_start) / sin(mid_start_angle);
		double mid_end_angle = M_PI - angle_f_end - f_from_start_vertex->vertex_angle(end_vertex);
		double end = edge_on_end_vertex->length()*sin(angle_f_end) / sin(mid_end_angle);
		edge_pointer edge_opp = f_from_start_vertex->opposite_edge(pseudo_source);
		end = edge_opp->length() - end;
		assert((1+1e-10)*start <= end);
		if (start*(1 + 1e-6) < end)
		propagate_from_pseudo_source(pseudo_source, inside_traversed_area, f_from_start_vertex, start_vertex, start, end);
		}
		else{
		vertex_pointer start_vertex = edge_on_start_vertex->opposite_vertex(pseudo_source);
		double mid_start_angle = M_PI - angle_f_start - f_from_start_vertex->vertex_angle(start_vertex);
		double start = edge_on_start_vertex->length()*sin(angle_f_start) / sin(mid_start_angle);
		edge_pointer edge_opp = f_from_start_vertex->opposite_edge(pseudo_source);
		double end = edge_opp->length();
		assert((1 + 1e-10)*start <= end);
		if (start*(1 + 1e-6) < end)
		propagate_from_pseudo_source(pseudo_source, inside_traversed_area, f_from_start_vertex, start_vertex, start, end);

		incident_it = f_from_start_vertex;
		edge_on_start_vertex = incident_it->next_edge(edge_on_start_vertex, pseudo_source);
		incident_it = edge_on_start_vertex->opposite_face(incident_it);
		while (incident_it->id() != f_from_end_vertex->id()){
		// propagate from start to end face

		edge_pointer edge_it = incident_it->opposite_edge(pseudo_source);
		double end = edge_it->length();
		propagate_from_pseudo_source(pseudo_source, inside_traversed_area, incident_it, edge_it->v0(), 0, end);
		edge_on_start_vertex = incident_it->next_edge(edge_on_start_vertex, pseudo_source);
		incident_it = edge_on_start_vertex->opposite_face(incident_it);

		}
		start_vertex = edge_on_end_vertex->opposite_vertex(pseudo_source);
		double mid_end_angle = M_PI - angle_f_end - f_from_end_vertex->vertex_angle(start_vertex);
		start = edge_on_end_vertex->length()*sin(angle_f_end) / sin(mid_end_angle);
		end = f_from_end_vertex->opposite_edge(pseudo_source)->length();
		propagate_from_pseudo_source(pseudo_source, inside_traversed_area, f_from_end_vertex, start_vertex, start, end);

		}*/

	}


	inline void GeodesicAlgorithmExact::create_pseudo_source_windows(vertex_pointer &pseudo_source, bool inside_traversed_area)
	{
		// update vertices around pseudo_source

		for (unsigned i = 0; i < pseudo_source->adjacent_edges().size(); ++i)
		{
			edge_pointer   edge_it = pseudo_source->adjacent_edges()[i];
			vertex_pointer vert_it = edge_it->opposite_vertex(pseudo_source);

			double distance = pseudo_source->geodesic_distance() + edge_it->length();
			if (distance < vert_it->geodesic_distance())
			{
				clock_t start = clock();
				m_vertex_queue.erase(vert_it);
				clock_t stop = clock();
				clockpseudoprop += static_cast<double>(stop - start);
				vert_it->geodesic_distance() = distance;
				if (vert_it->state() == Vertex::OUTSIDE)
					vert_it->state() = Vertex::FRONT;

				vert_it->incident_face() = edge_it->adjacent_faces()[0];
				edge_pointer next_edge = vert_it->incident_face()->next_edge(edge_it, pseudo_source);
				vert_it->incident_point() = (next_edge->v0() == pseudo_source) ? 0 : next_edge->length();
				vert_it->edge_start_vertex() = edge_it;
				vert_it->p_start_angle() = 0;
				clock_t start_q = clock();
				m_vertex_queue.insert(vert_it);
				clock_t stop_q = clock();
				clockpseudoprop += static_cast<double>(stop_q - start_q);
			}
		}

		face_pointer f_incident = pseudo_source->incident_face();



		/*unsigned incident_index = -5;
		for (unsigned i = 0; i < pseudo_source->adjacent_faces().size(); i++){
		if (f_incident == pseudo_source->adjacent_faces()[i]){
		incident_index = i;
		break;
		}
		}
		unsigned neighbour_plus, neighbour_minus, neighbour_plus_two,neighbour_minus_two;
		if (incident_index == pseudo_source->adjacent_faces().size() - 1){
		neighbour_plus = 0;
		}
		else neighbour_plus = incident_index + 1;

		if (incident_index == 0)
		neighbour_minus = pseudo_source->adjacent_faces().size() - 1;
		else neighbour_minus = incident_index - 1;


		if (neighbour_plus == pseudo_source->adjacent_faces().size() - 1){
		neighbour_plus_two = 0;
		}
		else neighbour_plus_two = neighbour_plus + 1;

		if (neighbour_minus == 0)
		neighbour_minus_two = pseudo_source->adjacent_faces().size() - 1;
		else neighbour_minus_two = neighbour_minus - 1;
		*/
		// update pseudo_source windows around pseudo_source
		double sum_angle = 0;
		//if (pseudo_source->adjacent_faces().size() < 6)
		//printf("less than 6 faces\n");
		for (unsigned i = 0; i < pseudo_source->adjacent_faces().size(); ++i)
		{
			//sum_angle += pseudo_source->adjacent_faces()[i]->vertex_angle(pseudo_source);
			/*if (m_source != pseudo_source->id()){
			if (i == incident_index || i == neighbour_minus || i == neighbour_plus )continue;
			if (pseudo_source->adjacent_faces().size() >= 6 && (i == neighbour_minus_two || i == neighbour_plus_two))continue;
			}*/
			//else{
			//printf("First source\n");
			//}
			face_pointer face_it = pseudo_source->adjacent_faces()[i];
			edge_pointer edge_it = face_it->opposite_edge(pseudo_source);

			list_pointer list = (edge_it->adjacent_faces()[0] == face_it) ? list = interval_list_0(edge_it) : list = interval_list_1(edge_it);

			// create a window
			interval_pointer candidate = new Interval;
			candidate->distance_bound() = 1e10;
			candidate->start() = 0;
			candidate->stop() = edge_it->length();
			candidate->d() = pseudo_source->geodesic_distance();
			double angle = face_it->vertex_angle(list->start_vertex());
			double length = face_it->next_edge(edge_it, list->start_vertex())->length();
			candidate->pseudo_x() = cos(angle) * length;
			candidate->pseudo_y() = -sin(angle) * length;
		//	candidate->origin_vertex() = pseudo_source;
			edge_pointer left_edge = face_it->next_edge(edge_it, list->start_vertex());
			candidate->left_CH_distance() = left_edge->length();
//			candidate->left_cosine_CH_angle() = 1;
			candidate->right_CH_distance() = face_it->next_edge(left_edge, pseudo_source)->length();
	//		candidate->right_cosine_CH_angle() = 1;

			updateDistanceBound(candidate);
			// insert into list
			list->push_back(candidate);

			// push into M_LIST_QUEUE if inside traversed area
			if ((inside_traversed_area) && ((edge_it->v0()->state() != Vertex::FRONT) || (edge_it->v1()->state() != Vertex::FRONT)))
				m_list_queue.push(list);

			// Statistics
			++m_windows_wavefront;
			if (m_windows_peak < m_windows_wavefront)
				m_windows_peak = m_windows_wavefront;

		}
		//printf("print angle %lf\n",sum_angle*180/3.14);
	}

	//----------------- propagate a windows list (Rule 1) ---------------------
	inline void GeodesicAlgorithmExact::find_separating_point(list_pointer &list)
	{
		const double LOCAL_EPSILON = 1e-20 * list->edge()->length(); // numerical issue

		double L = Tri.left_edge->length();
		double top_x = L * cos(Tri.left_alpha);
		double top_y = L * sin(Tri.left_alpha);

		Vertex top_t; // temporary top_vertex
		memcpy(&top_t, Tri.top_vertex, sizeof(Vertex));
		top_t.geodesic_distance() = GEODESIC_INF;

		interval_pointer iter = list->begin();

		double wlist_sp = 0;
		double wlist_pseudo_x = 0;
		double wlist_pseudo_y = 0;
		bool directtovertex = false;
		while (iter != NULL)
		{
			interval_pointer &w = iter;

			double w_sp = w->pseudo_x() - w->pseudo_y() * ((top_x - w->pseudo_x()) / (top_y - w->pseudo_y()));
			double distance = GEODESIC_INF;
			bool tovertex;
			// shortest path from the window
			if ((w_sp - w->start() > LOCAL_EPSILON) && (w_sp - w->stop() < -LOCAL_EPSILON))
			{
				distance = w->d() + sqrt((top_x - w->pseudo_x()) * (top_x - w->pseudo_x()) + (top_y - w->pseudo_y()) * (top_y - w->pseudo_y()));
				w->shortest_distance() = distance;
			}
			else if (w_sp - w->start() <= LOCAL_EPSILON)
			{
				distance = w->d() + sqrt((top_x - w->start()) * (top_x - w->start()) + top_y * top_y) + sqrt((w->start() - w->pseudo_x()) * (w->start() - w->pseudo_x()) + w->pseudo_y() * w->pseudo_y());
				w->shortest_distance() = distance;
				w_sp = w->start();
			}
			else if (w_sp - w->stop() >= -LOCAL_EPSILON)
			{
				distance = w->d() + sqrt((top_x - w->stop()) * (top_x - w->stop()) + top_y * top_y) + sqrt((w->stop() - w->pseudo_x()) * (w->stop() - w->pseudo_x()) + w->pseudo_y() * w->pseudo_y());
				w->shortest_distance() = distance;
				w_sp = w->stop();
			}

			// update information at top_t
			if (distance < top_t.geodesic_distance())
			{

				//double differs = top_t.geodesic_distance() / distance;
				//if (differs > 1.01 && top_t.geodesic_distance() <10)
				//printf("%lf %lf \n",distance, differs);
				top_t.geodesic_distance() = distance;
				top_t.incident_face() = Tri.face;
				list->start_vertex();
				top_t.incident_point() = (list->start_vertex() == list->edge()->v0()) ? w_sp : list->edge()->length() - w_sp;
				wlist_sp = w_sp;
				wlist_pseudo_x = w->pseudo_x();
				wlist_pseudo_y = w->pseudo_y();
			}
			w->sp() = w_sp;

			iter = iter->next();
		}
		//const double EPS = 
		// update top_vertex and M_VERTEX_QUEUE
		if ((1 + 1e-10)*top_t.geodesic_distance() < Tri.top_vertex->geodesic_distance())
		{
			double middle_edge_len = sqrtf((top_x - wlist_sp)*(top_x - wlist_sp) + top_y*top_y);
			double p_start_angle = acos((L*L + middle_edge_len*middle_edge_len - wlist_sp*wlist_sp) / (2 * L*middle_edge_len));
			top_t.p_start_angle() = p_start_angle;
			top_t.edge_start_vertex() = Tri.left_edge;
			clock_t start_q = clock();
			if (Tri.top_vertex->state() == Vertex::FRONT) erase_from_queue(Tri.top_vertex);
			memcpy(Tri.top_vertex, &top_t, sizeof(Vertex));
			if (Tri.top_vertex->state() == Vertex::FRONT) m_vertex_queue.insert(Tri.top_vertex);
			clock_t stop_q = clock();
			clockpseudoprop += static_cast<double>(stop_q - start_q);

			if ((Tri.top_vertex->state() == Vertex::INSIDE)){
				//Tri.left_edge->opposite_face(top_t.incident_face());
				create_flat_v_fan_shape(Tri.top_vertex, true); // handle saddle vertex
			}


		}

		list->sp() = wlist_sp;
		list->pseudo_x() = wlist_pseudo_x;
		list->pseudo_y() = wlist_pseudo_y;
	}

	inline void GeodesicAlgorithmExact::propagate_windows_to_two_edges(list_pointer &list)
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

			/*if (w->compute_min_distance() > w->distance_bound()){
			iter = iter->next();
			list->erase(w);
			//delete w;
			//--m_windows_wavefront;

			continue;
			}*/

			assert(w->start() <= w->stop());

			if (w->sp() < list->sp() - LOCAL_EPSILON)
			{
				// only propagate to left edge
				double Intersect_X, Intersect_Y;

				// judge the positions of the two windows
				CalculateIntersectionPoint(list->pseudo_x(), list->pseudo_y(), list->sp(), 0, w->pseudo_x(), w->pseudo_y(), w->stop(), 0, Intersect_X, Intersect_Y);
				if ((w->stop() < list->sp()) || ((Intersect_Y <= 0) && (Intersect_Y >= list->pseudo_y()) && (Intersect_Y >= w->pseudo_y())))
				{
					direction = PropagationDirection::LEFT;
				}
				else
				{
					direction = PropagationDirection::BOTH;
				}
			}
			else if (w->sp() > list->sp() + LOCAL_EPSILON)
			{
				// only propagate to right edge
				double Intersect_X, Intersect_Y;

				// judge the positions of the two windows
				CalculateIntersectionPoint(list->pseudo_x(), list->pseudo_y(), list->sp(), 0, w->pseudo_x(), w->pseudo_y(), w->start(), 0, Intersect_X, Intersect_Y);
				if ((w->start() > list->sp()) || ((Intersect_Y <= 0) && (Intersect_Y >= list->pseudo_y()) && (Intersect_Y >= w->pseudo_y())))
				{
					direction = PropagationDirection::RIGHT;
				}
				else
				{
					direction = PropagationDirection::BOTH;
				}
			}
			else
			{
				// propagate to both edges
				direction = PropagationDirection::BOTH;
			}

			bool ValidPropagation;
			interval_pointer right_w;

			switch (direction) {
			case PropagationDirection::LEFT:
				ValidPropagation = compute_propagated_parameters(w->pseudo_x(),
					w->pseudo_y(),
					w->start(),
					w->stop(),
					Tri.left_alpha,
					Tri.left_edge->length(),
					w,
					w->d());

				iter_t = iter->next();
				if (ValidPropagation)
				{
					if (w->stop() == Tri.left_edge->length()){
						w->right_CH_distance() = sqrtf((w->stop() - w->pseudo_x())*(w->stop() - w->pseudo_x()) + w->pseudo_y()*w->pseudo_y());
						//w->right_cosine_CH_angle() = 1;
						// update distance bound here..
						updateDistanceBound(w);
					}
					/*else{
					double rightVertexCosine = compute_cosine_vector(w->stop() - w->pseudo_x(), -w->pseudo_y(), Tri.left_edge->length() - w->pseudo_x(), -w->pseudo_y());
					if (rightVertexCosine > w->right_cosine_CH_angle()){
					w->right_cosine_CH_angle() = rightVertexCosine;
					w->right_CH_distance() = sqrtf((Tri.left_edge->length() - w->pseudo_x())*(Tri.left_edge->length() - w->pseudo_x()) + w->pseudo_y()*w->pseudo_y());
					// update distance bound
					updateDistanceBound(w);
					}
					}*/
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

			case PropagationDirection::RIGHT:
				ValidPropagation = compute_propagated_parameters(Tri.bottom_edge->length() - w->pseudo_x(),
					w->pseudo_y(),
					Tri.bottom_edge->length() - w->stop(),
					Tri.bottom_edge->length() - w->start(),
					Tri.right_alpha,
					Tri.right_edge->length(),
					w,
					w->d());

				iter_t = iter->next();
				if (ValidPropagation)
				{
					double length = Tri.right_edge->length(); // invert window
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
					/*else{
					double leftVertexCosine = compute_cosine_vector(w->start() - w->pseudo_x(), -w->pseudo_y(), 0 - w->pseudo_x(), -w->pseudo_y());
					if (leftVertexCosine > w->left_cosine_CH_angle()){
					w->left_cosine_CH_angle() = leftVertexCosine;
					w->left_CH_distance() = sqrtf((0 - w->pseudo_x())*(0 - w->pseudo_x()) + w->pseudo_y()*w->pseudo_y());
					// update distance bound
					updateDistanceBound(w);
					}
					}*/



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

			case PropagationDirection::BOTH:
				right_w = new Interval;
				memcpy(right_w, w, sizeof(Interval));

				ValidPropagation = compute_propagated_parameters(w->pseudo_x(),
					w->pseudo_y(),
					w->start(),
					w->stop(),
					Tri.face->vertex_angle(Tri.left_vertex),
					Tri.left_edge->length(),
					w,
					w->d());

				iter_t = iter->next();
				if (ValidPropagation)
				{
					if (w->stop() == Tri.left_edge->length()){
						w->right_CH_distance() = sqrtf((w->stop() - w->pseudo_x())*(w->stop() - w->pseudo_x()) + w->pseudo_y()*w->pseudo_y());
					//	w->right_cosine_CH_angle() = 1;
						// update distance bound here..
						updateDistanceBound(w);
					}
					/*else{
					double rightVertexCosine = compute_cosine_vector(w->stop() - w->pseudo_x(), -w->pseudo_y(), Tri.left_edge->length() - w->pseudo_x(), -w->pseudo_y());
					if (rightVertexCosine > w->right_cosine_CH_angle()){
					w->right_cosine_CH_angle() = rightVertexCosine;
					w->right_CH_distance() = sqrtf((Tri.left_edge->length() - w->pseudo_x())*(Tri.left_edge->length() - w->pseudo_x()) + w->pseudo_y()*w->pseudo_y());

					// update distance bound
					updateDistanceBound(w);
					}
					}*/

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

				ValidPropagation = compute_propagated_parameters(Tri.bottom_edge->length() - right_w->pseudo_x(),
					right_w->pseudo_y(),
					Tri.bottom_edge->length() - right_w->stop(),
					Tri.bottom_edge->length() - right_w->start(),
					Tri.face->vertex_angle(Tri.right_vertex),
					Tri.right_edge->length(),
					right_w,
					right_w->d());

				if (ValidPropagation)
				{
					// invert window
					double length = Tri.right_edge->length();
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
					/*else{
					double leftVertexCosine = compute_cosine_vector(right_w->start() - right_w->pseudo_x(), -right_w->pseudo_y(), 0 - right_w->pseudo_x(), -right_w->pseudo_y());
					if (leftVertexCosine > right_w->left_cosine_CH_angle()){
					right_w->left_cosine_CH_angle() = leftVertexCosine;
					right_w->left_CH_distance() = sqrtf((0 - right_w->pseudo_x())*(0 - right_w->pseudo_x()) + right_w->pseudo_y()*right_w->pseudo_y());
					// update distance bound
					updateDistanceBound(right_w);
					}
					}*/

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

	//----------------- pairwise windows checking (Rule 2) ----------------------
	inline void GeodesicAlgorithmExact::check_with_vertices(list_pointer &list)
	{
		if (list->empty()) return;
		//return; // try to remove, what's happen
		interval_pointer iter = list->begin();
		interval_pointer iter_t;

		while ((!list->empty()) && (iter != NULL))
		{
			interval_pointer &w = iter;
			bool w_survive = true;
			if (w->compute_min_distance() > w->distance_bound())w_survive = false;
			edge_pointer   e = list->edge();
			vertex_pointer v1 = list->start_vertex();
			vertex_pointer v2 = e->opposite_vertex(v1);
			double d1 = GEODESIC_INF;

			d1 = w->d() + sqrt((w->stop() - w->pseudo_x()) * (w->stop() - w->pseudo_x()) + w->pseudo_y() * w->pseudo_y());
			if (v1->geodesic_distance() + w->stop() < d1)
				w_survive = false;

			d1 = w->d() + sqrt((w->start() - w->pseudo_x()) * (w->start() - w->pseudo_x()) + w->pseudo_y() * w->pseudo_y());
			if (v2->geodesic_distance() + e->length() - w->start() < d1)
				w_survive = false;
			/*if (w->d() > w->origin_vertex()->geodesic_distance())
			w_survive = false;*/
			//v1->
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

	inline windows_state GeodesicAlgorithmExact::check_between_two_windows(interval_pointer &w1, interval_pointer &w2)
	{
		double NUMERCIAL_EPSILON = 1 - 1e-12;
		// we implement the discussed 6 cases as follows for simplicity

		if ((w1->start() >= w2->start()) && (w1->start() <= w2->stop())) // w1->start
		{
			double Intersect_X, Intersect_Y;

			// judge the order of the two windows
			CalculateIntersectionPoint(w2->pseudo_x(), w2->pseudo_y(), w1->start(), 0, w1->pseudo_x(), w1->pseudo_y(), w1->stop(), 0, Intersect_X, Intersect_Y);

			if ((Intersect_Y <= 0) && (Intersect_Y >= w1->pseudo_y()) && (Intersect_Y >= w2->pseudo_y()))
			{
				double d1, d2;
				d1 = w1->d() + sqrt((w1->start() - w1->pseudo_x()) * (w1->start() - w1->pseudo_x()) + (w1->pseudo_y()) * (w1->pseudo_y()));
				d2 = w2->d() + sqrt((w1->start() - w2->pseudo_x()) * (w1->start() - w2->pseudo_x()) + (w2->pseudo_y()) * (w2->pseudo_y()));

				if (d2 < d1 * NUMERCIAL_EPSILON)
					return w1_invalid;
				/*if (d1 < d2 * NUMERCIAL_EPSILON)
				w2->start() = w1->start();*/
			}
		}

		if ((w1->stop() >= w2->start()) && (w1->stop() <= w2->stop())) // w1->stop
		{
			double Intersect_X, Intersect_Y;

			// judge the order of the two windows
			CalculateIntersectionPoint(w2->pseudo_x(), w2->pseudo_y(), w1->stop(), 0, w1->pseudo_x(), w1->pseudo_y(), w1->start(), 0, Intersect_X, Intersect_Y);

			if ((Intersect_Y <= 0) && (Intersect_Y >= w1->pseudo_y()) && (Intersect_Y >= w2->pseudo_y()))
			{
				double d1, d2;
				d1 = w1->d() + sqrt((w1->stop() - w1->pseudo_x()) * (w1->stop() - w1->pseudo_x()) + (w1->pseudo_y()) * (w1->pseudo_y()));
				d2 = w2->d() + sqrt((w1->stop() - w2->pseudo_x()) * (w1->stop() - w2->pseudo_x()) + (w2->pseudo_y()) * (w2->pseudo_y()));

				if (d2 < d1 * NUMERCIAL_EPSILON)
					return w1_invalid;
				/*if (d1 < d2 * NUMERCIAL_EPSILON)
				w2->stop() = w1->stop();*/
			}
		}

		if ((w2->start() >= w1->start()) && (w2->start() <= w1->stop())) // w2->start
		{
			double Intersect_X, Intersect_Y;

			// judge the previous order of the two windows
			CalculateIntersectionPoint(w1->pseudo_x(), w1->pseudo_y(), w2->start(), 0, w2->pseudo_x(), w2->pseudo_y(), w2->stop(), 0, Intersect_X, Intersect_Y);

			if ((Intersect_Y <= 0) && (Intersect_Y >= w1->pseudo_y()) && (Intersect_Y >= w2->pseudo_y()))
			{
				double d1, d2;
				d1 = w1->d() + sqrt((w2->start() - w1->pseudo_x()) * (w2->start() - w1->pseudo_x()) + (w1->pseudo_y()) * (w1->pseudo_y()));
				d2 = w2->d() + sqrt((w2->start() - w2->pseudo_x()) * (w2->start() - w2->pseudo_x()) + (w2->pseudo_y()) * (w2->pseudo_y()));

				if (d1 < d2 * NUMERCIAL_EPSILON)
					return w2_invalid;
				/*if (d2 < d1 * NUMERCIAL_EPSILON)
				w1->start() = w2->start();*/
			}
		}

		if ((w2->stop() >= w1->start()) && (w2->stop() <= w1->stop())) // w2->stop
		{
			double Intersect_X, Intersect_Y;

			// judge the previous order of the two windows
			CalculateIntersectionPoint(w1->pseudo_x(), w1->pseudo_y(), w2->stop(), 0, w2->pseudo_x(), w2->pseudo_y(), w2->start(), 0, Intersect_X, Intersect_Y);

			if ((Intersect_Y <= 0) && (Intersect_Y >= w1->pseudo_y()) && (Intersect_Y >= w2->pseudo_y()))
			{
				double d1, d2;
				d1 = w1->d() + sqrt((w2->stop() - w1->pseudo_x()) * (w2->stop() - w1->pseudo_x()) + (w1->pseudo_y()) * (w1->pseudo_y()));
				d2 = w2->d() + sqrt((w2->stop() - w2->pseudo_x()) * (w2->stop() - w2->pseudo_x()) + (w2->pseudo_y()) * (w2->pseudo_y()));

				if (d1 < d2 * NUMERCIAL_EPSILON)
					return w2_invalid;
				/*if (d2 < d1 * NUMERCIAL_EPSILON)
				w1->stop() = w2->stop();*/
			}
		}

		if (w1->start() >= w2->stop())
		{
			double Intersect_X, Intersect_Y;

			// judge the previous order of the two windows
			CalculateIntersectionPoint(w1->pseudo_x(), w1->pseudo_y(), w1->start(), 0, w2->pseudo_x(), w2->pseudo_y(), w2->stop(), 0, Intersect_X, Intersect_Y);

			face_pointer f = Tri.bottom_edge->opposite_face(Tri.face);
			edge_pointer e = f->next_edge(Tri.bottom_edge, Tri.left_vertex);
			double angle = f->vertex_angle(Tri.left_vertex);
			double Cx = e->length() * cos(angle);
			double Cy = e->length() * -sin(angle);

			if ((PointInTriangle(Intersect_X, Intersect_Y, Tri.bottom_edge->length(), Cx, Cy))
				&& (Intersect_Y <= 0) && (Intersect_Y >= w1->pseudo_y()) && (Intersect_Y >= w2->pseudo_y()))
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

			face_pointer f = Tri.bottom_edge->opposite_face(Tri.face);
			edge_pointer e = f->next_edge(Tri.bottom_edge, Tri.left_vertex);
			double angle = f->vertex_angle(Tri.left_vertex);
			double Cx = e->length() * cos(angle);
			double Cy = e->length() * -sin(angle);

			if ((PointInTriangle(Intersect_X, Intersect_Y, Tri.bottom_edge->length(), Cx, Cy))
				&& (Intersect_Y <= 0) && (Intersect_Y >= w1->pseudo_y()) && (Intersect_Y >= w2->pseudo_y()))
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

	inline void GeodesicAlgorithmExact::pairwise_windows_checking(list_pointer &list)
	{
		if (list->empty()) return;

		interval_pointer iter = list->begin();
		interval_pointer next, iter_t;

		next = iter->next();

		// traverse successive pairs of windows
		while ((!list->empty()) && (next != NULL))
		{
			windows_state ws = check_between_two_windows(iter, next);

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

	//-- compute distance bound based on Fast DGG -------------------
	inline void GeodesicAlgorithmExact::updateDistanceBound(interval_pointer& w){

		clock_t start = clock();
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
		//bool bound_updated = false;
		// calculate beta_near and beta_far
		/*if (w->left_cosine_CH_angle() != 1 || w->right_cosine_CH_angle() != 1){
		triangle_altitude = sqrtf(largerdistance*largerdistance - 0.25*bottom*bottom);
		double tan_half_winopening = (0.5*bottom)/triangle_altitude;
		// calculate beta_left
		if (w->left_cosine_CH_angle() != 1){
		// calculate tan(half_winopening - left_ch_angle)
		double left_sine_CH_angle = sqrtf(1 - w->left_cosine_CH_angle()*w->left_cosine_CH_angle());
		double tan_CH_angle = left_sine_CH_angle / w->left_cosine_CH_angle();
		double tan_difference = (tan_half_winopening - tan_CH_angle) / (1 + tan_half_winopening*tan_CH_angle);
		// calculate position of window intersection with bottom
		intersect_left = 0.5*bottom - tan_difference*triangle_altitude;

		}
		//calculate beta right
		if (w->right_cosine_CH_angle() != 1){
		// calculate tan(half_winopening - left_ch_angle)
		double right_sine_CH_angle = sqrtf(1 - w->right_cosine_CH_angle()*w->right_cosine_CH_angle());
		double tan_CH_angle = right_sine_CH_angle / w->right_cosine_CH_angle();
		double tan_difference = (tan_half_winopening - tan_CH_angle) / (1 + tan_half_winopening*tan_CH_angle);
		// calculate position of window intersection with bottom
		intersect_right = 0.5*bottom - tan_difference*triangle_altitude;
		}
		if (w->left_CH_distance() > w->right_CH_distance()){
		beta_far = 1 - intersect_left / bottom;
		beta_near = intersect_right / bottom;
		}
		else{
		beta_far = 1 - intersect_right / bottom;
		beta_near = intersect_left / bottom;
		}
		assert(beta_near <= (1+1e-9)*beta_far);
		if (beta_near > beta){
		CHdistance_ratio = 1;
		beta = 1-beta_near;
		}
		else if (beta_far < beta){
		//CHdistance_ratio = 1;
		beta = beta_far;
		}
		}*/

		bound = (beta)*(beta)*CHdistance_ratio*bottom*bottom / (2 * ScaleFactor* epsilon*(largerdistance - bottom)) + CHdistance_ratio*largerdistance;
		if (bound < w->distance_bound() && largerdistance > bottom)w->distance_bound() = bound;
		clock_t stop = clock();
		clockcountbound += static_cast<double>(stop)-static_cast<double>(start);
	}

	//------------------------- main operation ----------------------------
	inline void GeodesicAlgorithmExact::propagate_one_windows_list(list_pointer &list)
	{
		if (list->empty()) return;

		if (list->edge()->adjacent_faces().size() > 1)
		{
			// Rule 2: pairwise windows checking
			check_with_vertices(list);
		//	pairwise_windows_checking(list);

			// Rule 1: "One Angle Two Sides"
			find_separating_point(list);
			propagate_windows_to_two_edges(list);

		}
	}

	//-------------------------- main entry --------------------------
	inline void GeodesicAlgorithmExact::propagate(unsigned source)
	{
		int saddle = 0;
		clockcountbound = 0;
		clockpseudoprop = 0;

		//printf("Loading angles...\n");

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
			if (sum_angle < 1.9*M_PI){

				saddle++;
			}

		}
		printf("Sperical vertex with strange angle : %d\n", saddle);
		// initialization
		printf("Epsilon : %.12lf  Tolerance : %.12lf   side fan angle : %lf\n", epsilon, tolerance, side_fan_angle);
		m_windows_pseudo = 0;
		m_source = source;
		initialize_propagation_data();

		clock_t start = clock();
		Vertex& src = mesh()->vertices()[m_source];
		int count = 0;
		double average_err = 0;
		while (!m_vertex_queue.empty())
		{
			// (1) pop a vertex from M_VERTEX_QUEUE
			clock_t start_q = clock();
			vertex_pointer vert = *m_vertex_queue.begin();
			count++;
			if (count < 51) {
				//Point3D& src_xyz = src.xyz;
				//Point3D& vert_xyz = vert->xyz;
				double dist = sqrt((src.x() - vert->x())*(src.x() - vert->x()) + (src.y() - vert->y())*(src.y() - vert->y()) + (src.z() - vert->z())*(src.z() - vert->z()));
				double err = abs((dist - vert->geodesic_distance()) / vert->geodesic_distance());
				average_err += err;
			}
			else if (count == 51) {
				average_err /= 50;
				cout << "Error using straight line :" << average_err << endl;
			}

			//  vert->xyz;
			
			m_vertex_queue.erase(m_vertex_queue.begin());
			

			// (2) update wavefront
			vert->state() = Vertex::INSIDE;
			for (unsigned i = 0; i < vert->adjacent_edges().size(); ++i)
			{
				vertex_pointer vert_it = vert->adjacent_edges()[i]->opposite_vertex(vert);
				if (vert_it->state() == Vertex::OUTSIDE) vert_it->state() = Vertex::FRONT;
			}
			clock_t stop_q = clock();
			clockpseudoprop += static_cast<double>(stop_q - start_q);
			//printf("vertex id : %d\n",vert->id());
			// (3) handle saddle vertex
			if (vert->id() != m_source) create_flat_v_fan_shape(vert, false);
			 start_q = clock();
			// (4) push window lists on the wavefront incident to v into M_LIST_QUEUE
			for (unsigned i = 0; i < vert->adjacent_edges().size(); ++i)
			{
				edge_pointer edge_it = vert->adjacent_edges()[i];
				if (!interval_list_0(edge_it)->empty()) m_list_queue.push(interval_list_0(edge_it));
				if (!interval_list_1(edge_it)->empty()) m_list_queue.push(interval_list_1(edge_it));
			}

			for (unsigned i = 0; i < vert->adjacent_faces().size(); ++i)
			{
				edge_pointer   edge_it = vert->adjacent_faces()[i]->opposite_edge(vert);
				vertex_pointer vert_it = (edge_it->adjacent_faces().size() < 2) ? NULL : edge_it->opposite_face(vert->adjacent_faces()[i])->opposite_vertex(edge_it);
				if (edge_it->adjacent_faces().size() < 2 || vert_it->state() != Vertex::OUTSIDE)
				{
					if (!interval_list_0(edge_it)->empty()) m_list_queue.push(interval_list_0(edge_it));
					if (!interval_list_1(edge_it)->empty()) m_list_queue.push(interval_list_1(edge_it));
				}
			}
			 stop_q = clock();
			clockpseudoprop += static_cast<double>(stop_q - start_q);
			// (5) propagate window lists in a FIFO order
			while (!m_list_queue.empty())
			{
				// pop an list from M_LIST_QUEUE
				list_pointer list = m_list_queue.front();
				m_list_queue.pop();

				bool is_boundary = calculate_triangle_parameters(list, Tri);

				if (!is_boundary)
				{
					// propagate the window list using Rule 1 and 2
					wl_left.clear(); wl_right.clear();
					propagate_one_windows_list(list);
					 start_q = clock();
					// merge windows lists
					if (!wl_left.empty())
					{
						// in VTP, both "PrimeMerge" and "SecondMerge" connect window lists in an order-free way
						if (!Tri.left_list->empty())
						{
							Tri.left_list->begin()->previous() = wl_left.end();
							wl_left.end()->next() = Tri.left_list->begin();
							Tri.left_list->begin() = wl_left.begin();
						}
						else
						{
							Tri.left_list->begin() = wl_left.begin();
							Tri.left_list->end() = wl_left.end();
						}

						// push updated list into M_LIST_QUEUE
						if (((Tri.left_edge->v0()->state() == Vertex::INSIDE) || (Tri.left_edge->v1()->state() == Vertex::INSIDE)) && (!Tri.left_list->empty()))
							m_list_queue.push(Tri.left_list);
					}

					if (!wl_right.empty())
					{
						// in VTP, both "PrimeMerge" and "SecondMerge" connect window lists in an order-free way
						if (!Tri.right_list->empty())
						{
							Tri.right_list->end()->next() = wl_right.begin();
							wl_right.begin()->previous() = Tri.right_list->end();
							Tri.right_list->end() = wl_right.end();
						}
						else
						{
							Tri.right_list->begin() = wl_right.begin();
							Tri.right_list->end() = wl_right.end();
						}

						// push updated list into M_LIST_QUEUE
						if (((Tri.right_edge->v0()->state() == Vertex::INSIDE) || (Tri.right_edge->v1()->state() == Vertex::INSIDE)) && (!Tri.right_list->empty()))
							m_list_queue.push(Tri.right_list);
						
					}
					stop_q = clock();
					clockpseudoprop += static_cast<double>(stop_q - start_q);
				}

				list->clear();
			}

			// statistics
			if (m_vertex_queue.size() > m_queue_max_size)
				m_queue_max_size = m_vertex_queue.size();
		}

		clock_t stop = clock();
		m_time_consumed = (static_cast<double>(stop)-static_cast<double>(start)) / CLOCKS_PER_SEC;
		printf("time for updateDistance : %lf\n", clockcountbound / CLOCKS_PER_SEC);
		printf("time for vertex queue, fan-shape propagation, merge list, push list : %lf\n", clockpseudoprop / CLOCKS_PER_SEC);
	}

	//---------------------- print statistics --------------------------
	inline void GeodesicAlgorithmExact::print_statistics()
	{
		GeodesicAlgorithmBase::print_statistics();

		double memory = sizeof(Interval);

		//std::cout << std::endl;
		std::cout << "Peak number of intervals on wave-front " << m_windows_peak << std::endl;
		std::cout << "uses about " << memory * m_windows_peak / 1e6 << "MB of memory" << std::endl;
		std::cout << "total interval propagation number " << m_windows_propagation << std::endl;
		std::cout << "maximum interval queue size is " << m_queue_max_size << std::endl;
		std::cout << "number interval pseudo " << m_windows_pseudo << std::endl;

		FILE* file = fopen("AVTP_DGG_Stats.csv", "a");
		double peakmem = memory * m_windows_peak / 1e6;
		fprintf(file, "%.6lf, %lld, %d,%.6lf, %.6lf, %d, \n", m_time_consumed, m_windows_propagation,
			m_windows_peak, peakmem, clockpseudoprop / CLOCKS_PER_SEC, m_windows_pseudo);
		fclose(file);
	}

}		//geodesic

#endif
