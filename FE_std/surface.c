
#include "surface.h"
#include "shapefunc.h"
#include "BB/std.h"
#include "BB/calc.h"
#include "mapping.h"

#include <math.h>
#include <stdbool.h>

void BBFE_set_surface_get_outward_normal_vector(
		const int local_num_nodes,
		double**  local_x,
		double*   dN_dxi,
		double*   dN_det,
		double    n_vec[3])
{
	// calc tangent vectors on the surface
	double dx_dxi[3];  double dx_det[3];
	BBFE_std_mapping_vector3d(
			dx_dxi, local_num_nodes, local_x, dN_dxi);
	BBFE_std_mapping_vector3d(
			dx_det, local_num_nodes, local_x, dN_det);

	// calc normal vector on the surface
	BB_calc_vec3d_cross(n_vec, dx_dxi, dx_det);
}


int BBFE_std_surface_get_num_surfs_in_elem(
		int  local_num_nodes)
{
	int num_surfs_in_elem = 0;
	switch(local_num_nodes) {
		case 4 : // tet 1st
		case 10: // tet 2nd
			num_surfs_in_elem = 4;
			break;

		case 8 : // hex 1st
		case 27: // hex 2nd
			num_surfs_in_elem = 6;
			break;
	}

	return num_surfs_in_elem;
}


int BBFE_std_surface_get_num_nodes_on_surf(
		int  local_num_nodes)
{
	int num_nodes_in_surf = 0;
	switch(local_num_nodes) {
		case 4 : // tet 1st
			num_nodes_in_surf = 3;
			break;
		case 10: // tet 2nd
			num_nodes_in_surf = 6;
			break;

		case 8 : // hex 1st
			num_nodes_in_surf = 4;
			break;
		case 27: // hex 2nd
			num_nodes_in_surf = 9;
			break;
	}

	return num_nodes_in_surf;
}


int BBFE_std_surface_get_surface_node_3d(
		bool*    node_is_on_surface,
		int      total_num_nodes,
		double** x,
		int      total_num_elems,
		int      local_num_nodes,
		int**    conn)
{
	int num_surface_nodes = 0;
	int num_nodes_on_surf = 0;
	
	switch(local_num_nodes) {
		case 4 : // tet 1st
			num_nodes_on_surf = 3;
			num_surface_nodes =
				BBFE_std_surface_tet1st_get_surface_node(
						node_is_on_surface,
						total_num_nodes,
						x,
						total_num_elems,
						conn);
			break;
		
		case 10: // tet 2nd
			num_nodes_on_surf = 6;
			num_surface_nodes =
				BBFE_std_surface_tet2nd_get_surface_node(
						node_is_on_surface,
						total_num_nodes,
						x,
						total_num_elems,
						conn);
			break;


		case 8 : // hex 1st
		case 27: // hex 2nd
			num_nodes_on_surf = 4;
			num_surface_nodes = 
				BBFE_std_surface_hex1st_get_surface_node(
						node_is_on_surface,
						total_num_nodes,
						x,
						total_num_elems,
						conn);
			break;

	}

	return num_surface_nodes;
}


int BBFE_std_surface_get_surface(
		bool**   surf_is_on_surface, //[total_num_elems][num_surfs]
		bool*    node_is_on_surface,
		int      total_num_elems,
		int      local_num_nodes,
		int**    conn)
{
	int num_surfs = 0;

	switch(local_num_nodes) {
		case 4:
		case 10:
			num_surfs = BBFE_std_surface_tet1st_get_surface(
					surf_is_on_surface,
					node_is_on_surface,
					total_num_elems,
					conn);

			break;

		case 8:
		case 27:
			num_surfs = BBFE_std_surface_hex1st_get_surface(
					surf_is_on_surface,
					node_is_on_surface,
					total_num_elems,
					conn);

			break;
	}

	return num_surfs;
}


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
		if(len > 0.1) {
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
		if(len > 0.1) {
			num_surface_nodes++;
			node_is_on_surface[i] = true;
		}
	}
	
	BB_std_free_2d_double(norm, total_num_nodes, 3);

	return num_surface_nodes;
}


int BBFE_std_surface_tet2nd_get_surface_node(
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
			int surf_conn[6];
			BBFE_std_shapefunc_tet2nd_get_surface(
					surf_conn, i);

			double ans[3];  double vec_1[3];  double vec_2[3];
			int nid_0 = conn[e][ surf_conn[0] ];
			int nid_1 = conn[e][ surf_conn[1] ];
			int nid_2 = conn[e][ surf_conn[2] ];
			int nid_3 = conn[e][ surf_conn[3] ];
			int nid_4 = conn[e][ surf_conn[4] ];
			int nid_5 = conn[e][ surf_conn[5] ];
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
				norm[ nid_4 ][d] += ans[d];
				norm[ nid_5 ][d] += ans[d];
			}

		}
	}
	
	int num_surface_nodes = 0;
	for(int i=0; i<(total_num_nodes); i++) {
		double len = BB_calc_vec3d_length(norm[i]);
		if(len > 0.1) {
			num_surface_nodes++;
			node_is_on_surface[i] = true;
		}
	}
	
	BB_std_free_2d_double(norm, total_num_nodes, 3);

	return num_surface_nodes;
}


int BBFE_std_surface_tet1st_get_surface(
		bool**   surf_is_on_surface, //[total_num_elems][num_surfs]
		bool*    node_is_on_surface,
		int      total_num_elems,
		int**    conn)
{
	int total_num_surfs = 0;

	const int num_surf_nodes = 3;
	const int num_surfs = 4;

	for(int e=0; e<total_num_elems; e++) {

		for(int s=0; s<num_surfs; s++) {
			int loc_conn[ num_surf_nodes ];
			int nid[ num_surf_nodes ];
			
			BBFE_std_shapefunc_tet1st_get_surface(loc_conn, s);
			for(int i=0; i<num_surf_nodes; i++) {
				nid[i] = conn[e][ loc_conn[i] ];
			}

			int count = 0;
			for(int i=0; i<num_surf_nodes; i++) {
				if( node_is_on_surface[ nid[i] ] ) {
					count++;
				}
			}

			if(count == num_surf_nodes) {
				surf_is_on_surface[e][s] = true;
			}
			else {
				surf_is_on_surface[e][s] = false;
			}
		}

	}

	for(int e=0; e<total_num_elems; e++) {

		for(int s=0; s<num_surfs; s++) {

			if( surf_is_on_surface[e][s] ) {
				int loc_conn[ num_surf_nodes ];
				int nid[ num_surf_nodes ];
				
				BBFE_std_shapefunc_tet1st_get_surface(loc_conn, s);
				for(int i=0; i<num_surf_nodes; i++) {
					nid[i] = conn[e][ loc_conn[i] ];
				}

				if( BBFE_std_surface_tet1st_search_identical_surface(
						nid, total_num_elems, conn, e) ) {
					surf_is_on_surface[e][s] = false;
				}
				else {
					surf_is_on_surface[e][s] = true;
					total_num_surfs++;
				}
			}

		}
	}

	return total_num_surfs;
}


int BBFE_std_surface_hex1st_get_surface(
		bool**   surf_is_on_surface, //[total_num_elems][num_surfs]
		bool*    node_is_on_surface,
		int      total_num_elems,
		int**    conn)
{
	int total_num_surfs = 0;

	const int num_surf_nodes = 4;
	const int num_surfs = 6;

	for(int e=0; e<total_num_elems; e++) {

		for(int s=0; s<num_surfs; s++) {
			int loc_conn[ num_surf_nodes ];
			int nid[ num_surf_nodes ];
			
			BBFE_std_shapefunc_hex1st_get_surface(loc_conn, s);
			for(int i=0; i<num_surf_nodes; i++) {
				nid[i] = conn[e][ loc_conn[i] ];
			}

			int count = 0;
			for(int i=0; i<num_surf_nodes; i++) {
				if( node_is_on_surface[ nid[i] ] ) {
					count++;
				}
			}

			if(count == num_surf_nodes) {
				surf_is_on_surface[e][s] = true;
			}
			else {
				surf_is_on_surface[e][s] = false;
			}
		}

	}

	for(int e=0; e<total_num_elems; e++) {

		for(int s=0; s<num_surfs; s++) {

			if( surf_is_on_surface[e][s] ) {
				int loc_conn[ num_surf_nodes ];
				int nid[ num_surf_nodes ];
				
				BBFE_std_shapefunc_hex1st_get_surface(loc_conn, s);
				for(int i=0; i<num_surf_nodes; i++) {
					nid[i] = conn[e][ loc_conn[i] ];
				}

				if( BBFE_std_surface_hex1st_search_identical_surface(
						nid, total_num_elems, conn, e) ) {
					surf_is_on_surface[e][s] = false;
				}
				else {
					surf_is_on_surface[e][s] = true;
					total_num_surfs++;
				}
			}

		}
	}

	return total_num_surfs;
}


static bool identical_elements(
		int* n1,
		int* n2,
		int  num)
{
	bool identical[num];

	for(int i=0; i<num; i++) {
		identical[i] = false;
	}

	for(int i=0; i<num; i++) {
		for(int j=0; j<num; j++) {
			if( n1[i] == n2[j] ) {
				identical[i] = true;
			}
		}
	}

	for(int i=0; i<num; i++) {
		if(!identical[i]) { return false;}
	}

	return true;
}


bool BBFE_std_surface_tet1st_search_identical_surface(
		int*  surf_nodes,
		int   total_num_elems,
		int** conn,
		int   elem_num)
{
	const int num_surf_nodes = 3;
	const int num_surfs = 4;
	
	for(int e=0; e<total_num_elems; e++) {
		if( e == elem_num) { continue; }

		for(int s=0; s<num_surfs; s++) {
			int loc_conn[ num_surf_nodes ];
			int nid[ num_surf_nodes ];
			
			BBFE_std_shapefunc_tet1st_get_surface(loc_conn, s);
			for(int i=0; i<num_surf_nodes; i++) {
				nid[i] = conn[e][ loc_conn[i] ];
			}
			
			if( identical_elements(surf_nodes, nid, num_surf_nodes) ) {
				return true;
			}
		}
	}

	return false;
}


bool BBFE_std_surface_hex1st_search_identical_surface(
		int*  surf_nodes,
		int   total_num_elems,
		int** conn,
		int   elem_num)
{
	const int num_surf_nodes = 4;
	const int num_surfs = 6;
	
	for(int e=0; e<total_num_elems; e++) {
		if( e == elem_num) { continue; }

		for(int s=0; s<num_surfs; s++) {
			int loc_conn[ num_surf_nodes ];
			int nid[ num_surf_nodes ];
			
			BBFE_std_shapefunc_hex1st_get_surface(loc_conn, s);
			for(int i=0; i<num_surf_nodes; i++) {
				nid[i] = conn[e][ loc_conn[i] ];
			}
			
			if( identical_elements(surf_nodes, nid, num_surf_nodes) ) {
				return true;
			}
		}
	}

	return false;
}
