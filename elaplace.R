rm(list=ls())

#expectiles gaussiens

#install.packages("devtools")

#require(devtools)

#install_version("expectreg",version="0.39")

#install.packages("VGAM")

#install.packages("expectreg")

#install.packages("Rcpp")

#install.packages("RcppArmadillo")

#install.packages("microbenchmark")

library(Rcpp)

library(VGAM)

library(expectreg)

library(RcppArmadillo)

library(microbenchmark)

cppFunction(depends="RcppArmadillo",plugins="cpp11",code="
            arma::vec explaplace(arma::vec& alphas, const double mean=0.0, const double sigma=1.0, const double lambda=1.0) {
              arma::vec alphabis=(abs((2*alphas)-1.0)+1)/2.0;
              arma::vec e(arma::size(alphas));
              arma::vec e1(arma::size(alphas));
              arma::vec ecart(arma::size(alphas));
              ecart.ones();
              e.zeros();
              alphabis.elem( find(alphabis<0) ).zeros();
              while(ecart.max()>0.0000001){
                e1=((1.0-2.0*alphabis)%exp(-sqrt(2.0/lambda)*e)%(sqrt(lambda/2.0)+e))/((2.0*alphabis-1.0)%(2-exp(-sqrt(2.0/lambda)*e))-2*alphabis);
                ecart=abs(e1-e);
                e=e1;
              }
              e.elem( find(alphabis<0.0000000000000001) ).fill(arma::datum::inf);
              e.elem( find(alphabis<0.0000000000000001) )=e.elem( find(alphabis<0.0000000000000001) )*(-1.0);
              e.elem( find(alphas<0) ).fill(arma::datum::nan);
              e.elem( find(alphas>1) ).fill(arma::datum::nan);
              if(alphas.max()>1 || alphas.min()<0){
                std::cout << \"Warning : Values of alpha must be between 0 and 1 \" << std::endl;
              }
              return(((e*=sigma)%=sign(alphas-0.5))+=mean);
            }
            ") #multiniveaux avec critère d'arrêt

lambda=1

microbenchmark(times=1000,qnorm(0.95),qenorm(0.95),enorm(0.95),explaplace(0.95))

microbenchmark(times=1000,qnorm(seq(0.05,0.95,0.05)),qenorm(seq(0.05,0.95,0.05)),enorm(seq(0.05,0.95,0.05)),explaplace(seq(0.05,0.95,0.05)))

plot(sign(1/2-seq(0.01,0.99,0.01))*sqrt(lambda/2)*log(1-abs(2*(1-seq(0.01,0.99,0.01))-1))~seq(0.01,0.99,0.01),type="l")

points(explaplace(seq(0.01,0.99,0.01))~seq(0.01,0.99,0.01),type="l",col="red")


