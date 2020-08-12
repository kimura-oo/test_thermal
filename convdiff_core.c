
#include "convdiff_core.h"


const char* BBFE_convdiff_get_directory_name(
		int         argc,
		char*       argv[],
		const char* codename)
{
	const char* dir_name;

	if(argc < 2) { dir_name = "."; }
	else         { dir_name = argv[1]; }

	printf("%s Main directory: %s\n", codename, dir_name);

	return dir_name;
}


void BBFE_convdiff_pre(
		FE_SYSTEM* sys,
		int        argc,
		char*      argv[],
		int        num_integ_points_each_axis,
		bool       manufactured_solution)
{
	sys->cond.directory = 
		BBFE_convdiff_get_directory_name(argc, argv, CODENAME);	

	int n_axis = num_integ_points_each_axis;

	BBFE_sys_read_node(
			&(sys->fe),
			INPUT_FILENAME_NODE,
			sys->cond.directory);
	BBFE_sys_read_elem(
			&(sys->fe),
			INPUT_FILENAME_ELEM,
			sys->cond.directory,
			n_axis*n_axis*n_axis);

	BBFE_sys_memory_allocation_integ(
			&(sys->basis),
			n_axis*n_axis*n_axis,
			3);
	BBFE_sys_memory_allocation_shapefunc(
			&(sys->basis),
			sys->fe.local_num_nodes,
			1,
			n_axis*n_axis*n_axis);

	BBFE_convdiff_memory_allocation_nodal_values(
			&(sys->vals),
			sys->fe.total_num_nodes);

	if(manufactured_solution) {
		switch( sys->fe.local_num_nodes ) {
			case 4:
				BBFE_manusol_overwrite_bc_file_tet(
						&(sys->fe),
						BLOCK_SIZE, 
						INPUT_FILENAME_D_BC,
						sys->cond.directory);
				break;

			case 8:
				BBFE_manusol_overwrite_bc_file_hex(
						&(sys->fe),
						BLOCK_SIZE, 
						INPUT_FILENAME_D_BC,
						sys->cond.directory);
				break;
		}
	}

	BBFE_sys_read_Dirichlet_bc(
			&(sys->bc),
			INPUT_FILENAME_D_BC,
			sys->cond.directory,
			sys->fe.total_num_nodes);

	BBFE_convdiff_set_basis(
			&(sys->basis),
			sys->fe.local_num_nodes,
			n_axis);

	monolis_initialize(&(sys->monolis));
	monolis_get_nonzero_pattern(
			&(sys->monolis),
			sys->fe.total_num_nodes,
			sys->fe.local_num_nodes,
			1,
			sys->fe.total_num_elems,
			sys->fe.conn);
}


void BBFE_convdiff_memory_allocation_nodal_values(
		NODAL_VALUES*   vals,
		const int       total_num_nodes)
{
	vals->T        = BB_std_calloc_1d_double(vals->T,     total_num_nodes);
	vals->error    = BB_std_calloc_1d_double(vals->error, total_num_nodes);
	vals->theo_sol = BB_std_calloc_1d_double(vals->error, total_num_nodes);
}


void BBFE_convdiff_set_basis(
		FE_3D_BASIS*  basis,
		int           local_num_nodes,
		int           num_integ_points_each_axis)
{
	switch( local_num_nodes ) {
		case 4:
			basis->num_integ_points = 
				BBFE_std_integ_tet_set_arbitrary_points(
						num_integ_points_each_axis,
						basis->integ_point,
						basis->integ_weight);

			for(int i=0; i<(basis->num_integ_points); i++) {
				BBFE_std_shapefunc_tet1st_get_val(
						basis->integ_point[i],
						basis->N[i]);

				BBFE_std_shapefunc_tet1st_get_derivative(
						basis->integ_point[i],
						basis->dN_dxi[i],
						basis->dN_det[i],
						basis->dN_dze[i]);
			}
			printf("%s Element type: 1st-order tetrahedron.\n", CODENAME);
			break;

		case 8:
			basis->num_integ_points = 
				BBFE_std_integ_hex_set_arbitrary_points(
						num_integ_points_each_axis,
						basis->integ_point,
						basis->integ_weight);

			for(int i=0; i<(basis->num_integ_points); i++) {
				BBFE_std_shapefunc_hex1st_get_val(
						basis->integ_point[i],
						basis->N[i]);

				BBFE_std_shapefunc_hex1st_get_derivative(
						basis->integ_point[i],
						basis->dN_dxi[i],
						basis->dN_det[i],
						basis->dN_dze[i]);
			}
			printf("%s Element type: 1st-order hexahedron.\n", CODENAME);
			break;
	}
	printf("%s The number of integration points: %d\n", CODENAME, basis->num_integ_points);
}


void BBFE_convdiff_finalize(
		FE_SYSTEM* sys)
{
	BBFE_sys_memory_free_integ(&(sys->basis), 3);
	BBFE_sys_memory_free_shapefunc(&(sys->basis));

	BBFE_sys_memory_free_node(&(sys->fe), 3);
	BBFE_sys_memory_free_elem(&(sys->fe), sys->basis.num_integ_points, 3);
	
	BBFE_sys_memory_free_Dirichlet_bc(&(sys->bc), sys->fe.total_num_nodes, BLOCK_SIZE);
}
