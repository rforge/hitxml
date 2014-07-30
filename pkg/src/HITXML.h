#ifndef HIT_XML_H_
#define HIT_XML_H_

#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <math.h>

void XML_add_spot(const int* n,
					const double pixel_x_mm[],
					const double pixel_y_mm[],
					const double* spot_x_mm,
					const double* spot_y_mm,
					const double* focus_x_mm,
					const double* focus_y_mm,
					const double* particles,
					double pixel_value[]);

void XML_PBR (const int *n_df,
		const double df_x[],
		const double df_y[],
		const double df_focusX[],
		const double df_focusY[],
		const double df_particles[],
		const int* n_m,
		const double m_x[],
		const double m_y[],
		double m_value[]);

#endif /* HIT_XML_H_ */
