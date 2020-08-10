#pragma once


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
