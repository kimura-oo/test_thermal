#pragma once

#include <stdbool.h>

void BBFE_set_surface_get_outward_normal_vector(
		const int local_num_nodes,
		double**  local_x,
		double*   dN_dxi,
		double*   dN_det,
		double    n_vec[3]);

int BBFE_std_surface_get_num_surfs_in_elem(
		int  local_num_nodes);

int BBFE_std_surface_get_num_nodes_on_surf(
		int  local_num_nodes);

int BBFE_std_surface_get_surface_node_3d(
		bool*    node_is_on_surface,
		int      total_num_nodes,
		double** x,
		int      total_num_elems,
		int      local_num_nodes,
		int**    conn);

int BBFE_std_surface_get_surface(
		bool**   surf_is_on_surface, //[total_num_elems][num_surfs]
		bool*    node_is_on_surface,
		int      total_num_elems,
		int      local_num_nodes,
		int**    conn);

int BBFE_std_surface_hex1st_get_surface_node(
		bool*    node_is_on_surface,
		int      total_num_nodes,
		double** x,
		int      total_num_elems,
		int**    conn);

int BBFE_std_surface_tet1st_get_surface_node(
		bool*    node_is_on_surface,
		int      total_num_nodes,
		double** x,
		int      total_num_elems,
		int**    conn);

int BBFE_std_surface_tet2nd_get_surface_node(
		bool*    node_is_on_surface,
		int      total_num_nodes,
		double** x,
		int      total_num_elems,
		int**    conn);

int BBFE_std_surface_tet1st_get_surface(
		bool**   surf_is_on_surface, //[total_num_elems][num_surfs]
		bool*    node_is_on_surface,
		int      total_num_elems,
		int**    conn);

int BBFE_std_surface_hex1st_get_surface(
		bool**   surf_is_on_surface, //[total_num_elems][num_surfs]
		bool*    node_is_on_surface,
		int      total_num_elems,
		int**    conn);

bool BBFE_std_surface_tet1st_search_identical_surface(
		int*  surf_nodes,
		int   total_num_elems,
		int** conn,
		int   elem_num);

bool BBFE_std_surface_hex1st_search_identical_surface(
		int*  surf_nodes,
		int   total_num_elems,
		int** conn,
		int   elem_num);
