#include <math.h>


void BB_calc_void(){}

/**********************************************************
 * 3D vector
 **********************************************************/
double BB_calc_vec3d_dot(
		const double  vec_1[3],
		const double  vec_2[3])
{
	return( vec_1[0]*vec_2[0] + vec_1[1]*vec_2[1] + vec_1[2]*vec_2[2] );
}


void BB_calc_vec3d_cross(
		double        ans[3],
		const double  vec_1[3],
		const double  vec_2[3])
{
	ans[0] = vec_1[1]*vec_2[2] - vec_1[2]*vec_2[1];
	ans[1] = vec_1[2]*vec_2[0] - vec_1[0]*vec_2[2];
	ans[2] = vec_1[0]*vec_2[1] - vec_1[1]*vec_2[0];
}


double BB_calc_vec3d_length(
		const double vec[3])
{
	double len = vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];
	return( sqrt(len) );
}


double BB_calc_vec3d_distance(
		const double x1[3],
		const double x2[3])
{
	return ( 
			sqrt( 
				(x1[0]-x2[0])*(x1[0]-x2[0]) + 
				(x1[1]-x2[1])*(x1[1]-x2[1]) +
				(x1[2]-x2[2])*(x1[2]-x2[2]) 
				) 
			);
}


void BB_calc_vec3d_normal_vec(
		double vec[3])
{
	double len = BB_calc_vec3d_length(vec);
	if(len == 0.0) {
		return;
	}
	else {
		for(int d=0; d<3; d++) {
			vec[d] = vec[d]/len;
		}
	}
}


void BB_calc_vec3d_copy(
		double       vec_copied[3],
		const double vec_orig[3])
{
	vec_copied[0] = vec_orig[0];
	vec_copied[1] = vec_orig[1];
	vec_copied[2] = vec_orig[2];
}


/**********************************************************
 * 3D matrix
 **********************************************************/
double BB_calc_mat3d_determinant(
		double mat[3][3])
{
	return ( mat[0][0]*mat[1][1]*mat[2][2] +
			mat[0][1]*mat[1][2]*mat[2][0] +
			mat[0][2]*mat[2][1]*mat[1][0] -
			mat[0][2]*mat[1][1]*mat[2][0] -
			mat[0][1]*mat[1][0]*mat[2][2] -
			mat[0][0]*mat[2][1]*mat[1][2]   );
}


void BB_calc_mat3d_inverse(
		double mat[3][3], 
		double det_mat,   
		double invMat[3][3])
{
	double invDet = 1.0/det_mat;

	invMat[0][0] = invDet * (mat[1][1]*mat[2][2] - mat[1][2]*mat[2][1]);
	invMat[0][1] = invDet * (mat[0][2]*mat[2][1] - mat[0][1]*mat[2][2]);
	invMat[0][2] = invDet * (mat[0][1]*mat[1][2] - mat[0][2]*mat[1][1]);

	invMat[1][0] = invDet * (mat[1][2]*mat[2][0] - mat[1][0]*mat[2][2]);
	invMat[1][1] = invDet * (mat[0][0]*mat[2][2] - mat[0][2]*mat[2][0]);
	invMat[1][2] = invDet * (mat[0][2]*mat[1][0] - mat[0][0]*mat[1][2]);

	invMat[2][0] = invDet * (mat[1][0]*mat[2][1] - mat[1][1]*mat[2][0]);
	invMat[2][1] = invDet * (mat[0][1]*mat[2][0] - mat[0][0]*mat[2][1]);
	invMat[2][2] = invDet * (mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0]);
}


void BB_calc_mat3d_copy(
		double mat_copied[3][3],
		double mat_orig[3][3])
{
	mat_copied[0][0] = mat_orig[0][0];
	mat_copied[0][1] = mat_orig[0][1];
	mat_copied[0][2] = mat_orig[0][2];

	mat_copied[1][0] = mat_orig[1][0];
	mat_copied[1][1] = mat_orig[1][1];
	mat_copied[1][2] = mat_orig[1][2];

	mat_copied[2][0] = mat_orig[2][0];
	mat_copied[2][1] = mat_orig[2][1];
	mat_copied[2][2] = mat_orig[2][2];
}
