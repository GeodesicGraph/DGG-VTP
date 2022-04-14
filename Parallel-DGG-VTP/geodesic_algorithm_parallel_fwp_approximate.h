//Parallel Approximate VTP

#ifndef GEODESIC_ALGORITHM_PARALLEL_FWP_APPROXIMATE_H
#define GEODESIC_ALGORITHM_PARALLEL_FWP_APPROXIMATE_H

#include "stdafx.h"
#include "geodesic_memory.h"
#include "geodesic_algorithm_base.h"
#include "geodesic_algorithm_exact_elements.h"
#include "FWP_Queue.h"
#include "omp.h"

namespace geodesic{

	class GeodesicAlgorithmParallelFWPApproximate : public GeodesicAlgorithmBase
	{
	public:
		GeodesicAlgorithmParallelFWPApproximate(geodesic::Mesh* mesh, double lambda, unsigned num_procs, unsigned k_concurrent) :
			GeodesicAlgorithmBase(mesh),
			lambda(lambda),
			m_vertex_FWP_queue(mesh->vertices().size(), 1),
			num_procs(num_procs),
			mx_concurrent_list_num(k_concurrent)
		{};
		~GeodesicAlgorithmParallelFWPApproximate(){};

		//clear
		void clear(){

			//for fwp 'bucket' structure
			m_iterations = 0;
			m_fwp_queue_max_size = 0;

			//for concurrently propagating window lists
			list_queue.resize(2);
			list_queue[0].clear();
			list_queue[1].clear();
			estimate_wavefront = 0;
			tris.resize(estimate_wavefront);
			wl_lefts.resize(estimate_wavefront);
			wl_rights.resize(estimate_wavefront);
			top_ts.resize(estimate_wavefront);
			isNotBoundary.resize(estimate_wavefront, 0);

			// for Approximate_VTP
			radius = lambda * this->m_mesh->avg_edge();
			guard = 8 * this->m_mesh->avg_edge();
			list_queue_appr.clear();

			//for statistics
			num_propagated_vertices = 0;
			num_propagated_local_exact_vertex = 0;
			num_propagated_pseudo_source_windows = 0;
			num_rings = 0;


			list_concurrent_time = 0;

			GeodesicAlgorithmBase::clear();
		};

		//main entry
		void propagate(unsigned source);


		//print the resulting statistics
		void print_statistics();

	private:

		// simple functions
		void initialize_propagation_data();
		void create_pseudo_source_windows(vertex_pointer &v, bool UpdateFIFOQueue, unsigned list_index);
		void erase_from_queue(vertex_pointer v);

		//for fwp's structure
		void erase_from_FWP_queue(vertex_pointer v);

		// propagate a windows list (Rule 1)
		void find_separating_point(list_pointer &list, Triangle &tri, Vertex &top_t); // find the separating point of the windows and the list
		void propagate_windows_to_two_edges(list_pointer &list, Triangle &tri, IntervalList &wl_left, IntervalList &wl_right); // propagates windows to two edges accross a triangle face

		// pairwise windows checking (Rule 2)
		void check_with_vertices(list_pointer &list);
		windows_state check_between_two_windows(interval_pointer &w1, interval_pointer &w2, Triangle &tri); // Check two neighbouring crossing windows on same edge
		void pairwise_windows_checking(list_pointer &list, Triangle &tri); // Check crossing windows on same edge

		// approximate operation (Du Jie)
		void propagate_local_exact_pseudo_source_at_wavefront(unsigned list_index);

		std::set<vertex_pointer, Vertex> m_vertex_queue;

		//for fwp 'bucket' structure
		FWP_Queue<vertex_pointer> m_vertex_FWP_queue;
		unsigned m_iterations;
		unsigned m_fwp_queue_max_size;

		std::queue<vertex_pointer> m_local_exact_vertex_queue;
		std::queue<list_pointer> m_list_queue;
		unsigned m_source;

		//parameter for approximate: lambda times 
		double lambda;
		double radius;
		double guard;
		double radius_guard;
		std::vector<list_pointer> list_queue_appr;

		//parameter for parallel; parameter for selecting suitable num of window lists
		unsigned num_procs;
		unsigned prevPhase;
		unsigned concurrent_num;
		unsigned mx_concurrent_list_num;
		unsigned from_list, to_list;
		std::vector<std::vector<list_pointer>> list_queue;
		
		//list_queue[0]: from_queue, list_queue[1]: next_queue;
		
		//for propagate lists concurrently, update vertices and merge lists sequentially
		unsigned estimate_wavefront;
		std::vector<int> isNotBoundary;
		std::vector<Triangle> tris;
		std::vector<IntervalList> wl_lefts;
		std::vector<IntervalList> wl_rights;
		std::vector<Vertex> top_ts;

		//for statistics
		unsigned num_propagated_vertices;
		unsigned num_propagated_local_exact_vertex;
		unsigned num_propagated_pseudo_source_windows;
		unsigned num_rings;
		double list_concurrent_time;
	};

	//----------------------simple functions-------------------
	inline void GeodesicAlgorithmParallelFWPApproximate::initialize_propagation_data()
	{
		clear();

		// initialize source's parameters
		vertex_pointer source = &(this->mesh()->vertices()[m_source]);
		source->geodesic_distance() = 0;
		source->state() = Vertex::INSIDE;

		//for fwp's 'bucket' structure
		unsigned int intialBinSize = 0;
		intialBinSize = source->adjacent_edges().size();
		m_vertex_FWP_queue.setParameters1(binWidth, intialBinSize);

		// initialize windows around source
		create_pseudo_source_windows(source, false, 0);

		m_vertex_FWP_queue.setParameters2(Kmin, Kmax, step);
		m_vertex_FWP_queue.setBinRange();
		m_vertex_FWP_queue.setPhase(1);
	}

	inline void GeodesicAlgorithmParallelFWPApproximate::erase_from_queue(vertex_pointer v)
	{
		assert(m_vertex_queue.count(v) <= 1);

		std::set<vertex_pointer, Vertex>::iterator it = m_vertex_queue.find(v);

		if (it != m_vertex_queue.end())
		{
			m_vertex_queue.erase(it);

		}
	}

	inline void GeodesicAlgorithmParallelFWPApproximate::erase_from_FWP_queue(vertex_pointer v)
	{
		if (v->ptrInQueue)
		{
			m_vertex_FWP_queue.remove(v->ptrInQueue);
			v->ptrInQueue = NULL;
		}
	}

	inline void GeodesicAlgorithmParallelFWPApproximate::create_pseudo_source_windows(vertex_pointer &pseudo_source, bool inside_traversed_area, unsigned list_index)
	{
		// update vertices around pseudo_source
		for (unsigned i = 0; i < pseudo_source->adjacent_edges().size(); ++i)
		{
			edge_pointer   edge_it = pseudo_source->adjacent_edges()[i];
			vertex_pointer vert_it = edge_it->opposite_vertex(pseudo_source);

			double distance = pseudo_source->geodesic_distance() + edge_it->length();
			if (distance < vert_it->geodesic_distance())
			{

				//if (vert_it->state() == Vertex::FRONT)//
					erase_from_FWP_queue(vert_it);
				//m_vertex_queue.erase(vert_it);

				vert_it->geodesic_distance() = distance;

				if (vert_it->state() == Vertex::OUTSIDE)
					vert_it->state() = Vertex::FRONT;

				//if (vert_it->state() == Vertex::FRONT)//
					vert_it->ptrInQueue = m_vertex_FWP_queue.push(vert_it, vert_it->geodesic_distance());
				//m_vertex_queue.insert(vert_it);

				vert_it->incident_face() = edge_it->adjacent_faces()[0];
				edge_pointer next_edge = vert_it->incident_face()->next_edge(edge_it, pseudo_source);
				vert_it->incident_point() = (next_edge->v0() == pseudo_source) ? 0 : next_edge->length();

			}
		}

		// update pseudo_source windows around pseudo_source
		for (unsigned i = 0; i < pseudo_source->adjacent_faces().size(); ++i)
		{
			face_pointer face_it = pseudo_source->adjacent_faces()[i];
			edge_pointer edge_it = face_it->opposite_edge(pseudo_source);
			list_pointer list = (edge_it->adjacent_faces()[0] == face_it) ? list = interval_list_0(edge_it) : list = interval_list_1(edge_it);

			num_propagated_pseudo_source_windows++;
			// create a window
			interval_pointer candidate = new Interval;

			//for debug
//			candidate->pseudo_id() = pseudo_source->id();
	//		candidate->mother_list_id() = -1;

			candidate->start() = 0;
			candidate->stop() = edge_it->length();
			candidate->d() = pseudo_source->geodesic_distance();
			double angle = face_it->vertex_angle(list->start_vertex());
			double length = face_it->next_edge(edge_it, list->start_vertex())->length();
			candidate->pseudo_x() = cos(angle) * length;
			candidate->pseudo_y() = -sin(angle) * length;

			// insert into list
			list->push_back(candidate);

			// push into M_LIST_QUEUE if inside traversed area
			/*
			if ((inside_traversed_area) && ((edge_it->v0()->state() != Vertex::FRONT) || (edge_it->v1()->state() != Vertex::FRONT)))
			{
				if (!list->readyPropagate())
				{
					list_queue[list_index].push_back(list);
					list->readyPropagate() = true;
				}
			}*/

			// Statistics
			++m_windows_wavefront;
			if (m_windows_peak < m_windows_wavefront)
				m_windows_peak = m_windows_wavefront;

		}
	}

	//----------------- propagate a windows list (Rule 1) modified by (Du Jie) ---------------------
	inline void GeodesicAlgorithmParallelFWPApproximate::find_separating_point(list_pointer &list, Triangle &tri, Vertex &top_t)
	{
		const double LOCAL_EPSILON = 1e-20 * list->edge()->length(); // numerical issue

		double L = tri.left_edge->length();
		double top_x = L * cos(tri.left_alpha);
		double top_y = L * sin(tri.left_alpha);

		//Vertex top_t; // temporary top_vertex
		memcpy(&top_t, tri.top_vertex, sizeof(Vertex));
		top_t.geodesic_distance() = GEODESIC_INF;

		interval_pointer iter = list->begin();

		double wlist_sp = 0;
		double wlist_pseudo_x = 0;
		double wlist_pseudo_y = 0;

		while (iter != NULL)
		{
			interval_pointer &w = iter;

			double w_sp = w->pseudo_x() - w->pseudo_y() * ((top_x - w->pseudo_x()) / (top_y - w->pseudo_y()));
			double distance = GEODESIC_INF;

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
				top_t.geodesic_distance() = distance;
				top_t.incident_face() = tri.face;
				top_t.incident_point() = (list->start_vertex() == list->edge()->v0()) ? w_sp : list->edge()->length() - w_sp;
				wlist_sp = w_sp;
				wlist_pseudo_x = w->pseudo_x();
				wlist_pseudo_y = w->pseudo_y();
			}
			w->sp() = w_sp;

			iter = iter->next();
		}

		list->sp() = wlist_sp;
		list->pseudo_x() = wlist_pseudo_x;
		list->pseudo_y() = wlist_pseudo_y;
	}

	inline void GeodesicAlgorithmParallelFWPApproximate::propagate_windows_to_two_edges(list_pointer &list, Triangle &tri, IntervalList &wl_left, IntervalList &wl_right)
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
				else
				{
					direction = BOTH;
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
					direction = RIGHT;
				}
				else
				{
					direction = BOTH;
				}
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
//					w->mother_list_id() = list->edge()->id();//
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

//					w->mother_list_id() = list->edge()->id();//
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
//					w->mother_list_id() = list->edge()->id();//
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

//					right_w->mother_list_id() = list->edge()->id();//
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
	inline void GeodesicAlgorithmParallelFWPApproximate::check_with_vertices(list_pointer &list)
	{
		if (list->empty()) return;

		interval_pointer iter = list->begin();
		interval_pointer iter_t;
		while ((!list->empty()) && (iter != NULL))
		{
			interval_pointer &w = iter;
			bool w_survive = true;

			edge_pointer   e = list->edge();
			vertex_pointer v1 = list->start_vertex();
			vertex_pointer v2 = e->opposite_vertex(v1);

			if (w->stop() - w->start() < 1e-8)
				w_survive = false;
			else
			{
				double d1 = GEODESIC_INF;

				d1 = w->d() + sqrt((w->stop() - w->pseudo_x()) * (w->stop() - w->pseudo_x()) + w->pseudo_y() * w->pseudo_y());
				if (v1->geodesic_distance() + w->stop() < d1)
					w_survive = false;

				d1 = w->d() + sqrt((w->start() - w->pseudo_x()) * (w->start() - w->pseudo_x()) + w->pseudo_y() * w->pseudo_y());
				if (v2->geodesic_distance() + e->length() - w->start() < d1)
					w_survive = false;
			}


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

	inline windows_state GeodesicAlgorithmParallelFWPApproximate::check_between_two_windows(interval_pointer &w1, interval_pointer &w2, Triangle &tri)
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
				if (d1 < d2 * NUMERCIAL_EPSILON)
					w2->start() = w1->start();
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
				if (d1 < d2 * NUMERCIAL_EPSILON)
					w2->stop() = w1->stop();
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
				if (d2 < d1 * NUMERCIAL_EPSILON)
					w1->start() = w2->start();
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

			face_pointer f = tri.bottom_edge->opposite_face(tri.face);
			edge_pointer e = f->next_edge(tri.bottom_edge, tri.left_vertex);
			double angle = f->vertex_angle(tri.left_vertex);
			double Cx = e->length() * cos(angle);
			double Cy = e->length() * -sin(angle);

			if ((PointInTriangle(Intersect_X, Intersect_Y, tri.bottom_edge->length(), Cx, Cy))
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

	inline void GeodesicAlgorithmParallelFWPApproximate::pairwise_windows_checking(list_pointer &list, Triangle &tri)
	{
		if (list->empty()) return;

		interval_pointer iter = list->begin();
		interval_pointer next, iter_t;

		next = iter->next();

		// traverse successive pairs of windows
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


	//----------------- approximate operation (Du Jie)----------------------------
	inline void GeodesicAlgorithmParallelFWPApproximate::propagate_local_exact_pseudo_source_at_wavefront(unsigned list_index)
	{
		//  ------------concurrent part---------------
		while (!m_vertex_FWP_queue.empty())
		{
			vertex_pointer vert = m_vertex_FWP_queue.top();
			if (vert->geodesic_distance() > radius_guard)
				break;
			m_vertex_FWP_queue.remove_without_record(vert->ptrInQueue);
			vert->ptrInQueue = NULL;
			vert->saddle_or_boundary() = true;
			//vert->state() = Vertex::EXACT;
			m_local_exact_vertex_queue.push(vert);
			num_propagated_local_exact_vertex++;
			for (unsigned i = 0; i < vert->adjacent_edges().size(); i++)
			{
				edge_pointer edge_it = vert->adjacent_edges()[i];
				vertex_pointer vert_it = edge_it->opposite_vertex(vert);
				
				if (vert_it->geodesic_distance() > radius_guard)//
					continue;//

				// prepare for deleting window lists
				if (!interval_list_0(edge_it)->empty() && !interval_list_0(edge_it)->readyDelete())
				{
					list_queue_appr.push_back(interval_list_0(edge_it));
					interval_list_0(edge_it)->readyDelete() = true;
				}
				if (!interval_list_1(edge_it)->empty() && !interval_list_1(edge_it)->readyDelete())
				{
					list_queue_appr.push_back(interval_list_1(edge_it));
					interval_list_1(edge_it)->readyDelete() = true;
				}
			}
		}

		
		tbb::tick_count st_concurrent_list = tbb::tick_count::now();
		tbb::parallel_for(tbb::blocked_range<size_t>(0, list_queue_appr.size()), [&](const tbb::blocked_range<size_t>& r)
		{
			for (size_t i = r.begin(); i != r.end(); i++)
			{
				list_pointer list = list_queue_appr[i];
				interval_pointer iter = list->begin();
				interval_pointer iter_t;
				while (iter)
				{
					iter_t = iter;
					iter = iter->next();
					list->erase(iter_t);
					delete iter_t;
					--m_windows_wavefront;
				}
				list->readyDelete() = false;
			}
		});
		tbb::tick_count nd_concurrent_list = tbb::tick_count::now();
		list_concurrent_time += (nd_concurrent_list - st_concurrent_list).seconds();

		list_queue_appr.clear();

		while (!m_local_exact_vertex_queue.empty())
		{
			vertex_pointer vert = m_local_exact_vertex_queue.front();
			m_local_exact_vertex_queue.pop();
			vert->ptrInQueue = m_vertex_FWP_queue.push(vert, vert->geodesic_distance());
		}

	}

	//----------------------main entry (Du Jie)---------------------
	inline void GeodesicAlgorithmParallelFWPApproximate::propagate(unsigned source)
	{
		printf("num_procs: %d\n", num_procs);
		//initialization
		m_source = source;
		initialize_propagation_data();
		prevPhase = 1;
		from_list = 0;

		//clock_t start = clock();
		//double st = omp_get_wtime();
		tbb::tick_count st = tbb::tick_count::now();
		//main propagation
		while (!m_vertex_FWP_queue.empty())
		{

			if (prevPhase > 1)
			{
				prevPhase++;
				m_vertex_FWP_queue.updateParamenters();
			}

			concurrent_num = 0;
			vertex_pointer vert = NULL;


			// pop suitable num of window lists
			do{
				//1) if too many window lists, break;
				if (list_queue[from_list].size() > mx_concurrent_list_num)
				{
					break;
				}
				m_vertex_FWP_queue.updatePhase();

				//2) if pop the vertex at next wavefront, break
				if (prevPhase < m_vertex_FWP_queue.getPhase())
				{
					prevPhase++;
					break;
				}
				// (1) pop a vertex from M_VERTEX_FWP_QUEUE
				vert = m_vertex_FWP_queue.top_for_parallel();
				
				//propagate local exact vertices
				if (vert->geodesic_distance() > radius)
				{
					m_vertex_FWP_queue.setParametersToInitial();
					vert = m_vertex_FWP_queue.top_for_parallel();

					//3) if the prepared local exact vertex is not the shortest vertex, set Parameters to initial and break,
					//   until next wavefront arrive to propagate the shorter vertices.
					if (vert->geodesic_distance() < radius)
					{
						//prevPhase = m_vertex_FWP_queue.getPhase();
						if (m_vertex_FWP_queue.getPhase() == 1)
							m_vertex_FWP_queue.setPhase(2);
						break;
					}

					//4) if arrive the prepared exact vertex for appr_VTP, break;
					if (concurrent_num > 0)
					{
						break;
					}
					
					num_rings++;
					radius_guard = radius + guard;

					//for approximate_VTP algorithm
					propagate_local_exact_pseudo_source_at_wavefront(from_list); // also can be concurrent
					
					if (m_vertex_FWP_queue.empty())
						break;


					radius += lambda * this->m_mesh->avg_edge();
					
					break;
				}

				//prepare for concurrently propagating window list
				// (1) pop a vertex from M_VERTEX_FWP_QUEUE
				m_vertex_FWP_queue.remove_without_record(vert->ptrInQueue);
				vert->ptrInQueue = NULL;

				num_propagated_vertices++;

				//for recording how many vertices regarding to concurrently propagating window lists
				concurrent_num++;

				// (2) update wavefront
				vert->state() = Vertex::INSIDE;
				for (unsigned i = 0; i < vert->adjacent_edges().size(); ++i)
				{
					vertex_pointer vert_it = vert->adjacent_edges()[i]->opposite_vertex(vert);
					if (vert_it->state() == Vertex::OUTSIDE) vert_it->state() = Vertex::FRONT;
				}

				// (3) handle saddle vertex
				if (vert->saddle_or_boundary()) create_pseudo_source_windows(vert, false, from_list);

				// (4) push window lists on the wavefront incident to 'vert' into 'list_queue[list_layer & 1]'
				for (unsigned i = 0; i < vert->adjacent_edges().size(); ++i)
				{
					edge_pointer edge_it = vert->adjacent_edges()[i];
					if (!interval_list_0(edge_it)->empty() && !interval_list_0(edge_it)->readyPropagate())
					{
						list_queue[from_list].push_back(interval_list_0(edge_it));
						interval_list_0(edge_it)->readyPropagate() = true;
					}
					if (!interval_list_1(edge_it)->empty() && !interval_list_1(edge_it)->readyPropagate())
					{
						list_queue[from_list].push_back(interval_list_1(edge_it));
						interval_list_1(edge_it)->readyPropagate() = true;
					}
				}

				for (unsigned i = 0; i < vert->adjacent_faces().size(); ++i)
				{
					edge_pointer   edge_it = vert->adjacent_faces()[i]->opposite_edge(vert);
					vertex_pointer vert_0 = edge_it->v0();
					vertex_pointer vert_1 = edge_it->v1();
					if ((edge_it->adjacent_faces().size() < 2) || (vert_0->state() == Vertex::INSIDE) || (vert_1->state() == Vertex::INSIDE))
					{
						if (!interval_list_0(edge_it)->empty() && !interval_list_0(edge_it)->readyPropagate())
						{
							list_queue[from_list].push_back(interval_list_0(edge_it));
							interval_list_0(edge_it)->readyPropagate() = true;
						}
						if (!interval_list_1(edge_it)->empty() && !interval_list_1(edge_it)->readyPropagate())
						{
							list_queue[from_list].push_back(interval_list_1(edge_it));
							interval_list_1(edge_it)->readyPropagate() = true;
						}
					}
				}

			} while (!m_vertex_FWP_queue.empty());

			if (concurrent_num < 1) continue;


			// (5) propagate window lists in 'list_queue[from_list]'
			while (!list_queue[from_list].empty())
			{
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
				tbb::parallel_for(tbb::blocked_range<size_t>(0, list_queue[from_list].size()), [&](const tbb::blocked_range<size_t>& r)
				{
					for (size_t i = r.begin(); i != r.end(); i++)
					{
						list_pointer list = list_queue[from_list][i];
						bool is_boundary = calculate_triangle_parameters(list, tris[i]);
						if (!is_boundary)
						{
							isNotBoundary[i] = 1;
							check_with_vertices(list);
							pairwise_windows_checking(list, tris[i]);

							find_separating_point(list, tris[i], top_ts[i]);
							wl_lefts[i].clear(); wl_rights[i].clear();
							propagate_windows_to_two_edges(list, tris[i], wl_lefts[i], wl_rights[i]);
							list->readyPropagate() = false;
						}
						//list_queue[from_list][i]->clear();
					}
				});
				tbb::tick_count nd_concurrent_list = tbb::tick_count::now();
				list_concurrent_time += (nd_concurrent_list - st_concurrent_list).seconds();
				
				to_list = from_list ^ 1;
				list_queue[to_list].clear();

				//sequentially: update vertex and merge lists
				for (unsigned i = 0; i < list_queue[from_list].size(); i++)
				{
					if (isNotBoundary[i])
					{
						isNotBoundary[i] = 0;
						//update vertex
						if (top_ts[i].geodesic_distance() < tris[i].top_vertex->geodesic_distance())
						{

							//if (tris[i].top_vertex->state() == Vertex::FRONT) //
								erase_from_FWP_queue(tris[i].top_vertex);

							memcpy(tris[i].top_vertex, &top_ts[i], sizeof(Vertex));

							//if (tris[i].top_vertex->state() == Vertex::FRONT)//
								tris[i].top_vertex->ptrInQueue = m_vertex_FWP_queue.push(tris[i].top_vertex, tris[i].top_vertex->geodesic_distance());

						}

						//merge windows lists
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
							}
						}
					}

				}
				from_list = to_list;
			}
			m_iterations++;
		}
		tbb::tick_count nd = tbb::tick_count::now();
		m_time_consumed = (nd - st).seconds();
		//double nd = omp_get_wtime();
		//m_time_consumed = nd - st;
	}

	
	//---------------------- print statistics --------------------------
	inline void GeodesicAlgorithmParallelFWPApproximate::print_statistics()
	{
		//GeodesicAlgorithmBase::print_statistics();

		double memory = sizeof(Interval);
		double used_memory = m_windows_peak*memory / 1e6;

		printf("\n");
		printf("--------------------Parallel-Appr-VTP statistics-------------------\n");
		printf("Propagation total time: %.6lf, concurrent part time: %.6lf\n", m_time_consumed, list_concurrent_time);

		printf("Peak number of windows on wavefront: %d\n", m_windows_peak);
		printf("Used memory: %.6lf\n", used_memory);

		printf("The number of propagated windows: %d\n", m_windows_propagation);
		printf("The number of propagated vertices: %d\n", num_propagated_vertices);
		printf("The number of propagated exact vertices: %d\n", num_propagated_local_exact_vertex);
		printf("The number of propagated pseudo windows: %d\n", num_propagated_pseudo_source_windows);
		printf("The number of rings: %d\n", num_rings);
		printf("The number of max_list_num: %d\n", estimate_wavefront);

		num_propagated_vertices += num_propagated_local_exact_vertex;

		//output result
		FILE* file = fopen("Parallel_FWP_Appr_VTP.csv", "a");
		fprintf(file, "%.6lf, %.6lf, %d, %d, %d, %d, %.6lf, %d, %d, %d\n", m_time_consumed, list_concurrent_time, num_rings, num_propagated_vertices, 
			m_windows_propagation, m_windows_peak, used_memory, m_iterations, mx_concurrent_list_num, estimate_wavefront);

		fclose(file);

	}


}//geodesic


#endif