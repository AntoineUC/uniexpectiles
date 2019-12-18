library(Rcpp)

library(RcppArmadillo)

cppFunction(depends="RcppArmadillo",plugins="cpp11",code="
            arma::vec epareto(arma::vec& alphas, double location=1.0, double shape=2) {
              arma::vec alphabis=alphas;
              arma::vec e(arma::size(alphas));
              arma::vec e1(arma::size(alphas));
              arma::vec ecart(arma::size(alphas));
              ecart.fill(location);
              e.ones();
              alphas.elem( find(alphas<0) ).fill(0.5);
              alphas.elem( find(alphas>0.999999999) ).fill(0.5);
              while(ecart.max()>0.0000001){
                e1=shape*location/(shape-1)+((pow(location,shape)/(shape-1))*((2*alphas-1)/(1-alphas)))%pow(e,1-shape);
                ecart=abs(e1-e);
                e=e1;
              }
              e.elem( find(alphas<0.0000000000000001) ).fill(location);
              e.elem( find(alphabis>0.999999999) ).fill(arma::datum::inf);
              e.elem( find(alphabis<0) ).fill(arma::datum::nan);
              e.elem( find(alphabis>1) ).fill(arma::datum::nan);
              if(alphabis.max()>1 || alphabis.min()<0){
                std::cout << \"Warning : Values of alpha must be between 0 and 1 \" << std::endl;
              }
              return(e);
            }
            ") 
