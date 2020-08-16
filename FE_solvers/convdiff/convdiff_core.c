
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
		BBFE_DATA*    fe,
		BBFE_BASIS*   basis,
		BBFE_BC*      bc,
		MONOLIS*      monolis,
		int           argc,
		char*         argv[],
		const char*   directory,
		int           num_integ_points_each_axis,
		bool          manufactured_solution)
{
	int n_axis = num_integ_points_each_axis;

	BBFE_sys_read_node(
			fe,
			INPUT_FILENAME_NODE,
			directory);
	BBFE_sys_read_elem(
			fe,
			INPUT_FILENAME_ELEM,
			directory,
			n_axis*n_axis*n_axis);

	BBFE_sys_memory_allocation_integ(
			basis,
			n_axis*n_axis*n_axis,
			3);
	BBFE_sys_memory_allocation_shapefunc(
			basis,
			fe->local_num_nodes,
			1,
			n_axis*n_axis*n_axis);

	BBFE_sys_read_Dirichlet_bc(
			bc,
			INPUT_FILENAME_D_BC,
			directory,
			fe->total_num_nodes,
			BLOCK_SIZE);

	BBFE_convdiff_set_basis(
			basis,
			fe->local_num_nodes,
			n_axis);

	monolis_initialize(monolis);
	monolis_get_nonzero_pattern(
			monolis,
			fe->total_num_nodes,
			fe->local_num_nodes,
			1,
			fe->total_num_elems,
			fe->conn);
}


void BBFE_convdiff_set_basis(
		BBFE_BASIS*   basis,
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
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		BBFE_BC*     bc)
{
	BBFE_sys_memory_free_integ(basis, 3);
	BBFE_sys_memory_free_shapefunc(basis);

	BBFE_sys_memory_free_node(fe, 3);
	BBFE_sys_memory_free_elem(fe, basis->num_integ_points, 3);
	
	BBFE_sys_memory_free_Dirichlet_bc(bc, fe->total_num_nodes, BLOCK_SIZE);
}
