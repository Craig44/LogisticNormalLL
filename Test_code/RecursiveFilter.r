## testing the recurisive Filter function written in C++ 
## Author C.Marsh
## date 11/10/09
Acf = c(0.4,0.2)
Acf = c(0.5,0.6)
Acf = c(0.5,0.4)
nBin = 16
xx = rep(0.0, nBin)

round(filter(xx, Acf, "recursive", init = rev(Acf)),4)

x = vector();
Rec_filter = function(Acf, xx, init_vals) {
  ## My C++ version
  x[1] = rev(init_vals)[1];
  x[2] = rev(init_vals)[2];
  for (i in 2:(length(xx) + 1)) {
    if (i == 1) {
      x[i] = x[i] * Acf[1] + x[i] * Acf[2];
    } else if ( i == 2) {
      x[i] =   x[i - 1] *Acf[1]  + x[i] *  Acf[2] 
    } else {  
      x[i] = x[i - 1] *  Acf[1] + x[i - 2] * Acf[2] ;
    } 
  }
  x[-1]
}


Rec_filter(Acf = c(0.3508772, 0.5001754), xx, init_vals = c(0.7020007, 0.7464914))