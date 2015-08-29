#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector XML_PBR_new(DataFrame beam_spot_grid,
                          NumericMatrix fluenceMatrix){
    NumericVector fluence(fluenceMatrix.nrow());
    
    NumericVector sigmaX_mm = beam_spot_grid["focus.X.FWHM.mm"];
    sigmaX_mm     = sigmaX_mm * 1.0 / (2 * sqrt(2 * log(2)));
    NumericVector sigX     = sigmaX_mm * sigmaX_mm;
 
    NumericVector sigmaY_mm = beam_spot_grid["focus.Y.FWHM.mm"];
    sigmaY_mm     = sigmaY_mm * 1.0 / (2 * sqrt(2 * log(2)));
    NumericVector sigY     = sigmaY_mm * sigmaY_mm;

    NumericVector prefactor = beam_spot_grid["particles"] / (2 * M_PI * sigmaX_mm * sigmaY_mm);
    
    NumericVector spot_x_mm = beam_spot_grid["x.mm"];
    NumericVector spot_y_mm = beam_spot_grid["y.mm"];
    
    
    for (int i = 0; i < beam_spot_grid.nrows(); i++){
     NumericVector dist_x_mm = fluenceMatrix(1,_) - spot_x_mm[i];
     NumericVector dist_y_mm = fluenceMatrix(2,_) - spot_y_mm[i];
     //     //TODO: insert distance criterion: ii            <- ((m[,1] - df$x.mm[i])/df$focusX.mm[i])^2 + ((m[,2] - df$y.mm[i])/df$focusY.mm[i])^2 <= 1
 
     fluence =  fluence + prefactor[i] * exp((-0.5 * ((dist_x_mm * dist_x_mm/ sigX[i]) + (dist_y_mm * dist_y_mm / sigY[i]))));
   }
  
  return(fluence);
  }

  
// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
XML_add_spot_new(1:5,
                 6:10,
                 1,
                 2,
                 3,
                 4,
                 5)
*/
