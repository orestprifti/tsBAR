// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

// -----------------------------------------------------------------------------
// UPDATE FUNCTION: T = 2*sqrt(lambda)
// -----------------------------------------------------------------------------
inline double updateCoefficient(double b, double la, double tol) {
    double T = 2.0 * std::sqrt(la);
    if (std::fabs(b) >= T) {
        double disc = b * b - T * T;
        if (disc < tol) disc = 0.0;
        return 0.5 * (b + std::copysign(1.0, b) * std::sqrt(disc));
    } else {
        return 0.0;
    }
}

// -----------------------------------------------------------------------------
// GLMNET (beta initialization)
// -----------------------------------------------------------------------------
NumericMatrix ic_glmnet(const NumericMatrix& y,
                    const NumericMatrix& x,
                    int n,
                    int p,
                    int np,
					std::string crit = "ebic",
                    double alpha = 0,
                    bool intercept = true,
                    bool standardize = true) {
	NumericMatrix coefficients(np,p);
	Function glmnet("glmnet");
	Function as_matrix("as.matrix");
	Function as_vector("as.vector");   	
	for (int i = 0; i < p; i++) {
		NumericVector yi = y(_,i);
		List model = glmnet(
			x, yi,
			Named("alpha") = alpha,
			Named("intercept") = intercept,
			Named("standardize") = standardize
		);
		S4 beta_S4 = model["beta"];
		NumericMatrix coefs = as_matrix(beta_S4);
		NumericVector a0 = as_vector(model["a0"]);
		int ncol = coefs.ncol();
		int nrow = coefs.nrow();
		NumericVector mse(ncol);
		NumericVector nvar(ncol);
		NumericVector criteria(ncol);
		for (int k = 0; k < ncol; ++k) {
			double sse = 0.0;
			int n_nonzero = 0;
			for (int r = 0; r < n; ++r) {
				double yhat = a0[k];
				for (int c = 0; c < nrow; ++c) yhat += x(r,c) * coefs(c,k);
				double e = yi[r] - yhat;
				sse += e * e;
			}
			mse[k] = sse / n;
			for (int c = 0; c < nrow; ++c) if (coefs(c,k) != 0.0) n_nonzero++;
			nvar[k] = n_nonzero;
			if (crit == "aic") {
				criteria[k] = n * std::log(mse[k]) + 2 * nvar[k];
			}
			else if (crit == "bic") {
				criteria[k] = n * std::log(mse[k]) + nvar[k] * std::log(n);
			}
			else if (crit == "ebic") {
				criteria[k] = n * log(mse[k]) + nvar[k] * log(n) + (2 * nvar[k] * 0.5 * log(np));
			}
			else {
				stop("Not a valid information criterion. Use: aic, bic, ebic.");
			}
		}
		int best_lambda = 0;
		double best_crit = criteria[0];
		for (int k = 1; k < ncol; ++k) {
			if (criteria[k] < best_crit) {
				best_crit = criteria[k];
				best_lambda = k;
			}
		}
		coefficients(_,i) = coefs(_,best_lambda);
	}
   
	return coefficients;
}

// -----------------------------------------------------------------------------
// BROKEN ADAPTIVE RIDGE
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List bar_ccd_var(
    const NumericMatrix& x,        
    const NumericMatrix& y,
    double lambda_max_input = 0.0,
    double lambda_min_input = 0.0,
    double c_lambda_min = 0.001,
    bool flg_lambda_min = true,
	bool intercept = true,
	bool standardize = true,
	double alpha_init = 0,
	Nullable<NumericMatrix> penalty_factor = R_NilValue,
	int nlambda = 100,
	std::string crit = "ebic",
    double epsilon = 1e-4,
    double thresh = 1e-9,
    int max_itt = 100000,
    const double tol = 1e-12, 
	const double sqrt_tol = 1e-12
) {
	NumericMatrix yy = clone(y);
	NumericMatrix xx = clone(x);
    int n = xx.nrow();
    int np = xx.ncol();
    int p = yy.ncol();
	if (n <= np) epsilon = 1e-2;
	//if (n <= np) c_lambda_min = 0.1;
	NumericVector ym(p,0.0);
	NumericVector xm(np,0.0);
	if (intercept) {
		for (int j = 0; j < p; ++j) {
			ym[j] = mean(yy(_,j));
			for (int i = 0; i < n; ++i) yy(i,j) -= ym[j];
		}
		for (int j = 0; j < np; ++j) {
			xm[j] = mean(xx(_,j));
			for (int i = 0; i < n; ++i) xx(i,j) -= xm[j];
		}
	}
	NumericVector ys(p,1.0);
	for (int j = 0; j < p; ++j) {
		if (sd(yy(_,j)) < 0.1) {
			ys[j] = sd(yy(_,j));
			for (int i = 0; i < n; ++i) yy[i] /= ys[j];
		}
	}
	NumericVector xs(np);
	if (standardize) {
		for (int j = 0; j < np; ++j) {
			xs[j] = sd(xx(_,j));
			if (xs[j] == 0.0) xs[j] = 1.0;
			for (int i = 0; i < n; ++i) xx(i,j) /= xs[j];
		}
	} else {
		std::fill(xs.begin(), xs.end(), 1.0);
	}
	NumericMatrix w_mat;
    if (penalty_factor.isNotNull()) {
		w_mat = NumericMatrix(penalty_factor);
		if (w_mat.nrow() != np || w_mat.ncol() != p) stop("penalty.factor must have dimensions ncol(x) x ncol(y)");
	} else {
		w_mat = NumericMatrix(np, p);
		std::fill(w_mat.begin(), w_mat.end(), 1.0);
	}
    NumericMatrix b_init = ic_glmnet(y, x, n, p, np, crit, alpha_init, intercept, standardize);
    NumericMatrix beta(np,p);
    NumericVector lambda(nlambda);
    NumericVector criteria(nlambda);
    NumericVector vec_lambda_max(p);
    NumericVector vec_lambda_min(p);
    NumericVector lambda_opt(p);
    IntegerVector vec_best_lambda(p);
    List List_beta(p);

	NumericMatrix vec_max_itt(p, nlambda);
    // Loop over response variables
    for (int col = 0; col < p; ++col) {
		
        std::vector<double> ycol(n);
        for (int i = 0; i < n; ++i) ycol[i] = yy(i,col);
		NumericVector w_vec = w_mat(_,col);
		
		// Penalty factor
		std::vector<double> w(np, 1.0);
		if (penalty_factor.isNotNull()) {
			double sum_w = 0.0;
			for (int j = 0; j < np; ++j) {
				double abs_w = std::fabs(w_vec[j]);
				double wj = std::max(abs_w, 0.0);
				w[j] = wj;
				sum_w += wj;
			}
		    // Normalize weights w.r.t. to the number of variables np
			for (int j = 0; j < np; ++j) w[j] = w[j] * ((double)np / sum_w);
		}
		
        std::vector<double> b(np), bm(np), bm_old(np), res(n);
        for (int j = 0; j < np; ++j) {
            b[j] = b_init(j, col);
            bm[j] = b[j];
        }

        // Initial gradient
        std::vector<double> g(np, 0.0);
        for (int j = 0; j < np; ++j) {
            double dot = 0.0;
            for (int i = 0; i < n; ++i) dot += xx(i,j) * ycol[i];
            g[j] = dot / n;
        }

        // Lambda sequence
        double lambda_max = lambda_max_input;
        if (lambda_max == 0.0) {
			double maxval = 0.0;
            for (int j = 0; j < np; ++j) {
                double xy = 0.0;
				double yjs = sd(yy(_,col));
                for (int i = 0; i < n; ++i) xy += xx(i, j) * (ycol[i] / yjs);
                double sq = (xy * xy) / (n * 4.0);
                if (sq > maxval) maxval = sq;
            }
            lambda_max = maxval;
        }
        vec_lambda_max[col] = lambda_max;

        double minval;
        if (flg_lambda_min && lambda_min_input == 0.0) {
            minval = c_lambda_min * (std::log(np) / n);
        } else if (flg_lambda_min && lambda_min_input != 0.0) {
            stop("Conflicting inputs: flg_lambda_min = TRUE and lambda_min > 0 are not allowed.");
        } else if (lambda_min_input != 0.0) {
            minval = lambda_min_input;
        } else {
            minval = epsilon * lambda_max;
        }
        vec_lambda_min[col] = minval;

        double logmax = std::log(lambda_max);
        double logmin = std::log(minval);
        for (int k = 0; k < nlambda; ++k) {
			if (nlambda==1) {
				lambda[0] = lambda_max_input;
			} else {
				lambda[k] = std::exp(logmax - ((logmax - logmin) * k / (nlambda - 1)));
			}
		}

        NumericMatrix b_tmp(np, nlambda);
		
        // ---------------------------------------------------------------------
        // Lambda loop
        // ---------------------------------------------------------------------
		bool lambda_max_found = false;
        for (int k = 0; k < nlambda; ++k) {
			
            double la = lambda[k];
			
            // Warm start
            for (int j = 0; j < np; ++j) bm[j] = b[j];
	
			// Strong rule
			std::vector<char> strong_set(np, 0);
            for (int j = 0; j < np; ++j) {
                double Tj = 2.0 * std::sqrt(la * w[j]);
                if (std::fabs(g[j]) >= Tj - tol) strong_set[j] = 1;
            }

            // Active set
            std::vector<int> active_idx;
            active_idx.reserve(np);
            for (int j = 0; j < np; ++j)
                if (std::fabs(b[j]) > 0.0 || strong_set[j]) active_idx.push_back(j);

            if (active_idx.empty()) {
                int best = 0;
                double gm = std::fabs(g[0]);
                for (int j = 1; j < np; ++j)
                    if (std::fabs(g[j]) > gm) { gm = std::fabs(g[j]); best = j; }
                active_idx.push_back(best);
            }

            // Order active set by b_init descending
            std::sort(active_idx.begin(), active_idx.end(), [&](int a, int bb) {
                return std::fabs(bm[a]) > std::fabs(bm[bb]);
            });

            // -----------------------------------------------------------------
            // COORDINATE DESCENT
            // -----------------------------------------------------------------
            for (int itt = 0; itt < max_itt; ++itt) {
				for (int j = 0; j < np; ++j) bm_old[j] = bm[j];
				
                // Compute residual
                for (int i = 0; i < n; ++i) {
                    double sum = 0.0;
                    for (int idx = 0; idx < (int)active_idx.size(); ++idx) {
                        int j = active_idx[idx];
                        sum += xx(i, j) * bm_old[j];
                    }
                    res[i] = ycol[i] - sum;
                }
                // Coordinate updates
                for (int idx = 0; idx < (int)active_idx.size(); ++idx) {
                    int j = active_idx[idx];
                    double dotprod = 0.0;
                    for (int i = 0; i < n; ++i) {
                        double rji = res[i] + xx(i,j) * bm_old[j];
                        dotprod += xx(i,j) * rji;
                    }
                    double b_j = dotprod / n;
                    bm[j] = updateCoefficient(b_j, la * w[j], sqrt_tol);
					double delta = bm[j] - bm_old[j];
					// Residual update in-place			   
					for (int i = 0; i < n; ++i) res[i] -= xx(i,j) * delta;
                }
                // Convergence check
                double diff = 0.0;
                for (int j = 0; j < np; ++j) {
                    double d = bm[j] - bm_old[j];
                    diff += d * d;
                }
                if (std::sqrt(diff) < thresh) break;
                // KKT every iterations
                if (true) {//itt % 2 == 1  
					for (int i = 0; i < n; ++i) {
						double sum = 0.0;
						for (int j = 0; j < np; ++j) sum += xx(i,j) * bm[j];
						res[i] = ycol[i] - sum;
					}
					bool new_added = false;
					// Boolean map for 'in_active' coordinates
					for (int j = 0; j < np; ++j) {
						bool in_active = false;
						for (int idx = 0; idx < (int)active_idx.size(); ++idx)
							if (active_idx[idx] == j) { in_active = true; break; }
						if (in_active) continue;
						double dot = 0.0;
						for (int i = 0; i < n; ++i) dot += xx(i,j) * res[i];
						double gfull = std::fabs(dot / n);
						// KKT
						double Tj = 2.0 * std::sqrt(la * w[j]);
						if (gfull >= Tj - tol) {
							active_idx.push_back(j);															 
							new_added = true;
						}
					}
					if (new_added) continue;  
                }
				vec_max_itt(col, k) = itt;
            }
            // Solution for specific lambda
			int nvar = 0;
			double a0_tmp = ym[col];		
			for (int j = 0; j < np; ++j) {
				b_tmp(j,k) = bm[j] * (ys[col] / xs[j]);
				a0_tmp -= xm[j] * b_tmp(j,k);
			}
			double mse = 0.0;
			for (int i = 0; i < n; ++i) {
				double yhat = a0_tmp;
				for (int j = 0; j < np; ++j) yhat += x(i,j) * b_tmp(j,k);
				double e = y(i,col) - yhat;
				mse += e * e;
			}
			mse /= n;
			for (int j = 0; j < np; ++j) if (std::fabs(b_tmp(j,k)) > 0.0) ++nvar;
			if (crit == "aic") {
				criteria[k] = n * std::log(mse) + 2 * nvar;
			}
			else if (crit == "bic") {
				criteria[k] = n * std::log(mse) + nvar * std::log(n);
			}
			else if (crit == "ebic") {
				criteria[k] = n * log(mse) + nvar * log(n) + (2 * nvar * 0.5 * log(np));
			}
			else {
				stop("Not a valid information criterion. Use: aic, bic, ebic.");
			}
			if (!lambda_max_found && nvar > 0) {
				lambda_max_found = true;
				double new_lambda_max = lambda[k-1];
				double new_minval = 0.0;
				if (flg_lambda_min && lambda_min_input == 0.0) {
					new_minval = c_lambda_min * (std::log(np) / n);			   
				} else if (flg_lambda_min && lambda_min_input != 0.0) {
					stop("Conflicting inputs: flg_lambda_min = TRUE and lambda_min > 0 are not allowed.");
				} else if (lambda_min_input != 0.0) {
					new_minval = lambda_min_input;
				} else {
					new_minval = epsilon * new_lambda_max;
				}
				vec_lambda_min[col] = new_minval;
				double logmax = std::log(new_lambda_max);
				double logmin = std::log(new_minval);
				for (int h = 0; h < nlambda; ++h) lambda[h] = std::exp(logmax - ((logmax - logmin) * h / (nlambda - 1)));
				b_tmp = NumericMatrix(np,nlambda);
				criteria = NumericVector(nlambda);
				k = -1;
				continue;
			}
        } // End lambda loop

        // Best lambda by criterion
        int best_lambda = 0;
        double best_crit = criteria[0];
        for (int k = 1; k < nlambda; ++k) {
            if (R_finite(criteria[k]) && criteria[k] < best_crit) {
                best_crit = criteria[k];
                best_lambda = k;
            }
        }
        vec_best_lambda[col] = best_lambda;
        lambda_opt[col] = lambda[best_lambda];
        for (int j = 0; j < np; ++j) beta(j,col) = b_tmp(j,best_lambda);
        List_beta[col] = b_tmp;
    } // End response variables loop
	
	// a0 + beta
	NumericVector a0(p);
	NumericMatrix coef(np+1, p);
	for (int col = 0; col < p; ++col) {
		a0[col] = ym[col];		
		for (int j = 0; j < np; ++j) {
			a0[col] -= xm[j] * beta(j,col);
			coef(j+1,col) = beta(j,col);
		}
		coef(0,col) = a0[col];
	}

    return List::create(
		_["a0"] = a0,
        _["beta"] = beta,
		_["coefficients"] = coef,
        _["lambda.opt"] = lambda_opt,
        _["best.lambda.index"] = vec_best_lambda,
        _["lambda.sequence"] = lambda,
		_["lambda.max"] = vec_lambda_max,
		_["lambda.min"] = vec_lambda_min,
        _["beta.path"] = List_beta,
		_["n.iterations"]=vec_max_itt
    );
}