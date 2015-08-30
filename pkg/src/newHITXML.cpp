#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector XML_PBR_new(DataFrame beam_spot_grid,
                          NumericMatrix fluenceMatrix){

    NumericVector fluence(fluenceMatrix.nrow());
    
    NumericVector sigmaX_mm = beam_spot_grid["focus.X.FWHM.mm"];
    sigmaX_mm     = sigmaX_mm / (2 * sqrt(2 * log(2)));
    NumericVector sigX     = sigmaX_mm * sigmaX_mm;
 
    NumericVector sigmaY_mm = beam_spot_grid["focus.Y.FWHM.mm"];
    sigmaY_mm     = sigmaY_mm / (2 * sqrt(2 * log(2)));
    NumericVector sigY     = sigmaY_mm * sigmaY_mm;

    NumericVector prefactor = beam_spot_grid["N.particles"];
    prefactor     = prefactor / (2 * M_PI * sigmaX_mm * sigmaY_mm);
    
    NumericVector spot_x_mm = beam_spot_grid["x.mm"];
    NumericVector spot_y_mm = beam_spot_grid["y.mm"];
    
    
    for (int i = 0; i < beam_spot_grid.nrows(); i++){
     NumericVector dist_x_mm = fluenceMatrix(_,0) - spot_x_mm[i];
     NumericVector dist_y_mm = fluenceMatrix(_,1) - spot_y_mm[i];
     
     //TODO: insert distance criterion: ii            <- ((m[,1] - df$x.mm[i])/df$focusX.mm[i])^2 + ((m[,2] - df$y.mm[i])/df$focusY.mm[i])^2 <= 1
     fluence = fluence + prefactor[i] * exp(-0.5 * (dist_x_mm * dist_x_mm/ sigX[i] + dist_y_mm * dist_y_mm / sigY[i]));
    }
  
  return(fluence);
  }

RcppExport void XML_add_spot(const int* n,
                  const double pixel_x_mm[],
                  const double pixel_y_mm[],
                  const double* spot_x_mm,
                  const double* spot_y_mm,
                  const double* focus_x_mm,
                  const double* focus_y_mm,
                  const double* particles,
                  double pixel_value[])

{
  double sigmaX_mm = * focus_x_mm / (2 * sqrt(2 * log(2)));
  double sigmaY_mm = * focus_y_mm / (2 * sqrt(2 * log(2)));
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

RcppExport void XML_PBR (const int *n_df,
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
// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
