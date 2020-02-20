library(Rcpp)

library(RcppArmadillo)

cppFunction(depends="RcppArmadillo",plugins="cpp11",code="
            arma::vec explog(arma::vec& alphas, const double mean=0.0, const double sigma=1.0) {
              arma::vec alphabis=(abs((2*alphas)-1.0)+1)/2.0;
              arma::vec e(arma::size(alphas));
              arma::vec e1(arma::size(alphas));
              arma::vec ecart(arma::size(alphas));
              ecart.ones();
              e.zeros();
              alphabis.elem( find(alphabis<0) ).zeros();
              while(ecart.max()>0.0000001){
                e1=(1.0-2.0*alphabis)%((exp(e)+1)%log(exp(e)+1)-e%exp(e))/((alphabis-1)%exp(e)-alphabis);
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
            ") 