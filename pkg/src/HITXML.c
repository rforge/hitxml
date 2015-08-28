#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <Rcpp.h>

// TODO: Use OpenMP to speed up
// TODO: Use SEXP rather than truely C objects
void XML_add_spot(const int* n,
					const double pixel_x_mm[],
					const double pixel_y_mm[],
					const double* spot_x_mm,
					const double* spot_y_mm,
					const double* focus_x_mm,
					const double* focus_y_mm,
					const double* particles,
					double pixel_value[])

	 {
		double sigmaX_mm = *focus_x_mm / (2 * sqrt(2 * log(2)));
		double sigmaY_mm = *focus_y_mm / (2 * sqrt(2 * log(2)));
		double sigX = sigmaX_mm * sigmaX_mm;
		double sigY = sigmaY_mm * sigmaY_mm;
		double prefactor = *particles / (2 * M_PI * sigmaX_mm * sigmaY_mm);
		int i;
		double dist_x_mm2, dist_y_mm2;
		for (i = 0; i < *n; i++){
			dist_x_mm2     = (pixel_x_mm[i] - *spot_x_mm) * (pixel_x_mm[i] - *spot_x_mm);
			dist_y_mm2     = (pixel_y_mm[i] - *spot_y_mm) * (pixel_y_mm [i] - *spot_y_mm);
			//TODO: insert distance criterion: ii            <- ((m[,1] - df$x.mm[i])/df$focusX.mm[i])^2 + ((m[,2] - df$y.mm[i])/df$focusY.mm[i])^2 <= 1
			pixel_value[i] +=  prefactor * exp((-1) * 0.5 * ((dist_x_mm2 / sigX) + (dist_y_mm2 / sigY)));
		}
	 }

void XML_PBR (const int *n_df,
		const double df_x[],
		const double df_y[],
		const double df_focusX[],
		const double df_focusY[],
		const double df_particles[],
		const int* n_m,
		const double m_x[],
		const double m_y[],
		double m_value[])
	{
		int i;
		for (i = 0; i < *n_m; i++){
			m_value[i] = 0;
		}

		for (i = 0; i < *n_df; i++)
		{
			XML_add_spot(n_m,
						 m_x,
						 m_y,
						 &(df_x[i]),
						 &(df_y[i]),
						 &(df_focusX[i]),
						 &(df_focusY[i]),
						 &(df_particles[i]),
						 m_value);
		}
}
