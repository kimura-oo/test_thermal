
#include "shapefunc.h"
#include "../libBB/std.h"
#include "../libBB/calc.h"


int BBFE_std_surface_hex1st_get_surface_node(
		bool*    node_is_on_surface,
		int      total_num_nodes,
		double** x,
		int      total_num_elems,
		int**    conn)
{
	double** norm;
	norm = BB_std_calloc_2d_double(norm, total_num_nodes, 3);

	for(int e=0; e<total_num_elems; e++) {
		for(int i=0; i<6; i++) {
			int surf_conn[4];
			BBFE_std_shapefunc_hex1st_get_surface(
					surf_conn, i);

			double ans[3];  double vec_1[3];  double vec_2[3];
			int nid_0 = conn[e][ surf_conn[0] ];
			int nid_1 = conn[e][ surf_conn[1] ];
			int nid_2 = conn[e][ surf_conn[2] ];
			int nid_3 = conn[e][ surf_conn[3] ];
			for(int d=0; d<3; d++) {
				vec_1[d] = x[ nid_1 ][d] - x[ nid_0 ][d];
				vec_2[d] = x[ nid_2 ][d] - x[ nid_0 ][d];
			}
			BB_calc_vec3d_cross(ans, vec_1, vec_2);
			BB_calc_vec3d_normal_vec(ans);

			for(int d=0; d<3; d++) {
				norm[ nid_0 ][d] += ans[d];
				norm[ nid_1 ][d] += ans[d];
				norm[ nid_2 ][d] += ans[d];
				norm[ nid_3 ][d] += ans[d];
			}

		}
	}
	
	int num_surface_nodes = 0;
	for(int i=0; i<(total_num_nodes); i++) {
		double len = BB_calc_vec3d_length(norm[i]);
		if(len > 1.0) {
			num_surface_nodes++;
			node_is_on_surface[i] = true;
		}
	}
	
	BB_std_free_2d_double(norm, total_num_nodes, 3);

	return num_surface_nodes;
}


int BBFE_std_surface_tet1st_get_surface_node(
		bool*    node_is_on_surface,
		int      total_num_nodes,
		double** x,
		int      total_num_elems,
		int**    conn)
{
	double** norm;
	norm = BB_std_calloc_2d_double(norm, total_num_nodes, 3);

	for(int e=0; e<total_num_elems; e++) {
		for(int i=0; i<4; i++) {
			int surf_conn[3];
			BBFE_std_shapefunc_tet1st_get_surface(
					surf_conn, i);

			double ans[3];  double vec_1[3];  double vec_2[3];
			int nid_0 = conn[e][ surf_conn[0] ];
			int nid_1 = conn[e][ surf_conn[1] ];
			int nid_2 = conn[e][ surf_conn[2] ];
			for(int d=0; d<3; d++) {
				vec_1[d] = x[ nid_1 ][d] - x[ nid_0 ][d];
				vec_2[d] = x[ nid_2 ][d] - x[ nid_0 ][d];
			}
			BB_calc_vec3d_cross(ans, vec_1, vec_2);
			BB_calc_vec3d_normal_vec(ans);

			for(int d=0; d<3; d++) {
				norm[ nid_0 ][d] += ans[d];
				norm[ nid_1 ][d] += ans[d];
				norm[ nid_2 ][d] += ans[d];
			}

		}
	}
	
	int num_surface_nodes = 0;
	for(int i=0; i<(total_num_nodes); i++) {
		double len = BB_calc_vec3d_length(norm[i]);
		if(len > 1.0) {
			num_surface_nodes++;
			node_is_on_surface[i] = true;
		}
	}
	
	BB_std_free_2d_double(norm, total_num_nodes, 3);

	return num_surface_nodes;
}
