#ifndef SVM_H
#define SVM_H

#include <linalg.h>

template <typename T>
void sdca(const Vector<T>& y, const Matrix<T>& X, Matrix<T>& W, const Vector<T>& tablambda, const T eps, const int max_iter, const int minibatch,const bool random_cycle) {
   const int n = y.n();
   const int p = X.m();
   const int nlambda=tablambda.n();
   W.resize(p,nlambda);
   W.setZeros();
   Vector<T> alpha(n);
   alpha.setZeros();
   Vector<T> z(n);
   z.setZeros();
   Vector<T> normX;
   X.norm_2sq_cols(normX);
   Vector<T> w;
   W.refCol(0,w);
   Vector<T> tmp;
   Vector<T> xi;
   Vector<T> ys;
   Matrix<T> Xs;
   Vector<T> grad;
   Vector<T> alphas;
   Vector<T> normXs;

   cout << "Problem size: p x n: " << p << " " << n << endl;
   for (int jj = 0; jj<nlambda; ++jj) {
      const T lambda=tablambda[jj];
      cout << "*********************" << endl;
      cout << "Processes Lambda " << lambda << endl;
      if (jj > 0) {
         W.refCol(jj,w);
         W.refCol(jj-1,tmp);
         w.copy(tmp);
         w.scal(tablambda[jj-1]/lambda);
      }
      for (int ii = 0; ii<max_iter; ++ii) {
         if (ii > 0 && (ii % (10*n/minibatch)) == 0) {
            T primalA=0;
            X.multTrans(w,z);
            for (int kk=0; kk<n; ++kk) primalA += MAX(0,1-y[kk]*z[kk]);
            primalA /= n;
            T primalB = w.nrm2sq();
            const T primal=primalA + 0.5*lambda*primalB;
            T dualA=-0.5*(lambda)*primalB;
            T dualB= y.dot(alpha)/n;
            const T dual = dualA + dualB;
            cout << "Iteration: " << ii << ", primal: " << primal << ", dual: " << dual << ", gap: " << (primal-dual) << endl;
            if ((primal - dual) < eps) break;
         }

         if (minibatch == 1) {
            const int ind = random_cycle ? random() % n : ii % n;
            const T yi=y[ind];
            X.refCol(ind,xi);
            const T A = normX[ind]/(n*lambda);
            const T deltaT = (yi-xi.dot(w))/A;
            const T delta=yi*MAX(0,MIN(1,yi*(alpha[ind]+deltaT)))-alpha[ind];
            alpha[ind]+=delta;
            w.add(xi,delta/(lambda*n));
         } else {
            const int ind = random() % n;
            const int sizebatch= MIN(minibatch,n-ind);
            y.refSubVec(ind,sizebatch,ys);
            normX.refSubVec(ind,sizebatch,normXs);
            X.refSubMat(ind,sizebatch,Xs);
            alpha.refSubVec(ind,sizebatch,alphas);
            Xs.mult(alpha,tmp);
            Xs.multTrans(tmp,grad);
            grad.scal(-T(1.0)/((lambda*n)*n));
            grad.add(ys,T(1.0)/n);
            const T L = sqrt(normXs.sum())/((lambda*n)*n);
            alphas.add(grad,-T(1.0)/L);
            for (int kk=0; kk<alphas.n(); ++kk) 
               if (ys[kk] > 0) {
                  alphas[kk]=MIN(MAX(alphas[kk],-1),0);
               } else {
                  alphas[kk]=MIN(MAX(alphas[kk],0),1);
               }
         }
      }
   }
}

#endif
