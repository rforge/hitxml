#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 10;
}

// [[Rcpp::export]]
double XML_add_spot_new(NumericVector pixel_x_mm,
  NumericVector pixel_y_mm,
  double spot_x_mm,
  double spot_y_mm,
  double focus_x_mm,
  double focus_y_mm,
  double particles) {

  double sigmaX_mm = focus_x_mm / (2 * sqrt(2 * log(2)));
  double sigmaY_mm = focus_y_mm / (2 * sqrt(2 * log(2)));
  double sigX = sigmaX_mm * sigmaX_mm;
  double sigY = sigmaY_mm * sigmaY_mm;
  double prefactor = particles / (2 * M_PI * sigmaX_mm * sigmaY_mm);

  NumericVector dist_x_mm = (pixel_x_mm - spot_x_mm);
  NumericVector dist_y_mm = (pixel_y_mm - spot_y_mm);
//     //TODO: insert distance criterion: ii            <- ((m[,1] - df$x.mm[i])/df$focusX.mm[i])^2 + ((m[,2] - df$y.mm[i])/df$focusY.mm[i])^2 <= 1
  NumericVector pixel_value =  prefactor * exp((-0.5 * ((dist_x_mm * dist_x_mm/ sigX) + (dist_y_mm * dist_y_mm / sigY))));

  return sum(pixel_value);
  }
  
  
  
  
// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
