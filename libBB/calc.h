#pragma once

void BB_calc_void();

/**********************************************************
 * 3D vector
 **********************************************************/
// rule: BB_calc_vec3d_
double BB_calc_vec3d_dot(
		const double  vec_1[3],
		const double  vec_2[3]);

void BB_calc_vec3d_cross(
		double        ans[3],
		const double  vec_1[3],
		const double  vec_2[3]);

double BB_calc_vec3d_length(
		const double vec[3]);

double BB_calc_vec3d_distance(
		const double x1[3],
		const double x2[3]);

void BB_calc_vec3d_normal_vec(
		double vec[3]);

void BB_calc_vec3d_copy(
		double       vec_copied[3],
		const double vec_orig[3]);
/**********************************************************
 * 3D matrix
 **********************************************************/
// rule: BB_calc_mat3d_
double BB_calc_mat3d_determinant(
		double mat[3][3]);

void BB_calc_mat3d_inverse(
		double mat[3][3], 
		double det_mat,   
		double invMat[3][3]);

void BB_calc_mat3d_copy(
		double mat_copied[3][3],
		double mat_orig[3][3]);
