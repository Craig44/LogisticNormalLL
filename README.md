# LogisticNormalLL
A repo that looks at implementing the logistic normal negative log likelihood for composition data. I have implemented it across two platforms (other than R where it exists) so that it 
would be useful for others. The two platforms are C++ that depends on Eigen thirdparty library and TMB. To look at comparisons see the R scripts which test the functionality 
across a range of age based compositions. This work follows Francis (2014) Replacing the multinomial in stock assessment models: A first step. Fisheries Research, 151, 70-84. 
I would like to thank Chris Francis for sharing his R-code with me which is the back bone for this repo.


##Using this code in TMB

Copy this [file](https://github.com/Craig44/LogisticNormalLL/blob/master/TMB/Logistic_normal_likelihood.hpp) into your TMB .cpp and include the following line in your model.cpp
```#include "Logistic_normal_likelihood.hpp"```. This will allow you to call all the functions in this script. The main two that I call are ```NLLlogistnorm(..)``` to evaluate the negative
likelihood and ```logis_norm_stand_resids(...)``` to calculate standardized residuals from the distribution. To see these functions in usage you can look at the ```model.cpp``` found in 
this [directory](https://github.com/Craig44/LogisticNormalLL/tree/master/TMB)


