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
NumericVector ic_glmnet(const NumericVector& y,
                    const NumericMatrix& x,
                    const int& n,
                    const int& p,
					const std::string& crit = "ebic",				   
                    const double& alpha = 0,
                    const bool& intercept = true,
                    const bool& standardize = true) {
	Function glmnet("glmnet");
	Function as_matrix("as.matrix");
	Function as_vector("as.vector");
	List model = glmnet(
	  x, y,
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
			double e = y[r] - yhat;
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
			criteria[k] = n * log(mse[k]) + nvar[k] * log(n) + (2 * nvar[k] * 0.5 * log(p));
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
	NumericVector coefficients = coefs(_,best_lambda);
   
  return coefficients;
}

// -----------------------------------------------------------------------------
// BROKEN ADAPTIVE RIDGE
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
List bar_ccd_reg(
    const NumericMatrix& x,        
    const NumericVector& y,
    double lambda_max_input = 0.0,
    double lambda_min_input = 0.0,
    double c_lambda_min = 0.001,
    bool flg_lambda_min = true,
    bool intercept = true,
	bool standardize = true,
	double alpha_init = 0,
	Nullable<NumericVector> penalty_factor = R_NilValue,
	int nlambda = 100,
	std::string crit = "ebic",
    double epsilon = 1e-4,
    double thresh = 1e-9,
    int max_itt = 100000,
    const double tol = 1e-12, 
	const double sqrt_tol = 1e-12
) {
	NumericVector yy = clone(y);
	NumericMatrix xx = clone(x);
    int n = xx.nrow();
    int p = xx.ncol();
	if (n <= p) epsilon = 1e-2;
	//if (n <= p) c_lambda_min = 0.1;
	NumericVector xm(p);
	double ym = 0.0;
	if (intercept) {
		ym = mean(yy);
		for (int i = 0; i < n; ++i) yy[i] -= ym;
		for (int j = 0; j < p; ++j) {
			xm[j] = mean(xx(_,j));
			for (int i = 0; i < n; ++i) xx(i,j) -= xm[j];
		}
	}
	double ys = 1.0;
	if (sd(yy) < 0.1) {
		ys = sd(yy);
		for (int i = 0; i < n; ++i) yy[i] /= ys;
	}
	NumericVector xs(p);
	if (standardize) {
		for (int j = 0; j < p; ++j) {
			xs[j] = sd(xx(_,j));
			if (xs[j] == 0.0) xs[j] = 1.0;
			for (int i = 0; i < n; ++i) xx(i,j) /= xs[j];
		}
	} else {
		std::fill(xs.begin(), xs.end(), 1.0);
	}
	NumericVector w_vec;
    if (penalty_factor.isNotNull()) {
		w_vec = NumericVector(penalty_factor);
		if (w_vec.size() != p) stop("penalty.factor must have length ncol(x)");
	} else {
		w_vec = NumericVector(p,1.0);
	}
	
	NumericVector b_init = ic_glmnet(y, x, n, p, crit, alpha_init, intercept, standardize);
    NumericVector lambda(nlambda);
    NumericVector criteria(nlambda);
	NumericVector vec_max_itt(nlambda);
        
	std::vector<double> beta(p), b(p), bm(p), bm_old(p), res(n);
	for (int j = 0; j < p; ++j) {
		b[j] = b_init[j];
		bm[j] = b[j];
	}
	
	// Penalty factor
	std::vector<double> w(p,1.0);
	if (penalty_factor.isNotNull()) { 
		double sum_w = 0.0;
		for (int j = 0; j < p; ++j) {
			double abs_w = std::fabs(w_vec[j]);
			double wj = std::max(abs_w,0.0);
			w[j] = wj;
			sum_w += wj;
		}
		// Normalize weights w.r.t. to the number of variables p
		for (int j = 0; j < p; ++j) w[j] = w[j] * ((double)p / sum_w);		
	}

	// Initial gradient
	std::vector<double> g(p,0.0);
	for (int j = 0; j < p; ++j) {
		double dot = 0.0;
		for (int i = 0; i < n; ++i) dot += xx(i,j) * yy[i];
		g[j] = dot / n;
	}

	// Lambda sequence
	double lambda_max = lambda_max_input;
	if (lambda_max == 0.0) {
		double maxval = 0.0;
		for (int j = 0; j < p; ++j) {
			double xy = 0.0;
			double yjs = sd(yy);
			for (int i = 0; i < n; ++i) xy += xx(i, j) * (yy[i] / yjs);
			double sq = (xy * xy) / (n * 4.0);
			if (sq > maxval) maxval = sq;
		}
		lambda_max = maxval;
	}

	double minval;
	if (flg_lambda_min && lambda_min_input == 0.0) {
		minval = c_lambda_min * (std::log(p) / n);
	} else if (flg_lambda_min && lambda_min_input != 0.0) {
		stop("Conflicting inputs: flg_lambda_min = TRUE and lambda_min > 0 are not allowed.");
	} else if (lambda_min_input != 0.0) {
		minval = lambda_min_input;
	} else {
		minval = epsilon * lambda_max;
	}

	double logmax = std::log(lambda_max);
	double logmin = std::log(minval);
	for (int k = 0; k < nlambda; ++k) {
		if (nlambda==1) {
			lambda[0] = lambda_max_input;
		} else {
			lambda[k] = std::exp(logmax - ((logmax - logmin) * k / (nlambda - 1)));
		}
	}

	NumericMatrix b_tmp(p, nlambda);
	
	// ---------------------------------------------------------------------
	// Lambda loop
	// ---------------------------------------------------------------------
	bool lambda_max_found = false;
	for (int k = 0; k < nlambda; ++k) {
		double la = lambda[k];
		
		// Warm start
		for (int j = 0; j < p; ++j) bm[j] = b[j];

		// Strong rule
		std::vector<char> strong_set(p, 0);
		for (int j = 0; j < p; ++j) {
			double Tj = 2.0 * std::sqrt(la * w[j]);
			if (std::fabs(g[j]) >= Tj - tol) strong_set[j] = 1;
		}

		// Active set
		std::vector<int> active_idx;
		active_idx.reserve(p);
		for (int j = 0; j < p; ++j)
			if (std::fabs(b[j]) > 0.0 || strong_set[j]) active_idx.push_back(j);

		if (active_idx.empty()) {
			int best = 0;
			double gm = std::fabs(g[0]);
			for (int j = 1; j < p; ++j)
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
			for (int j = 0; j < p; ++j) bm_old[j] = bm[j];
			
			// Compute residual
			for (int i = 0; i < n; ++i) {
				double sum = 0.0;
				for (int idx = 0; idx < (int)active_idx.size(); ++idx) {
					int j = active_idx[idx];
					sum += xx(i,j) * bm_old[j];
				}
				res[i] = yy[i] - sum;
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
			for (int j = 0; j < p; ++j) {
				double d = bm[j] - bm_old[j];
				diff += d * d;
			}
			if (std::sqrt(diff) < thresh) break;

			// KKT every iterations
			if (true) { //itt % 2 == 1					  
				for (int i = 0; i < n; ++i) {
					double sum = 0.0;
					for (int j = 0; j < p; ++j) sum += xx(i,j) * bm[j];
					res[i] = yy[i] - sum;
				}
				bool new_added = false;
				// Boolean map for 'in_active' coordinates
				for (int j = 0; j < p; ++j) {
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
			vec_max_itt[k] = itt;
		}
        // Solution for specific lambda
		int nvar = 0;
		double a0_tmp = ym;		
		for (int j = 0; j < p; ++j) {
			b_tmp(j,k) = bm[j] * (ys / xs[j]);
			a0_tmp -= xm[j] * b_tmp(j,k);
		}
		double mse = 0.0;
		for (int i = 0; i < n; ++i) {
			double yhat = a0_tmp;
			for (int j = 0; j < p; ++j) yhat += x(i,j) * b_tmp(j,k);
			double e = y[i] - yhat;
			mse += e * e;
		}
		mse /= n;
		for (int j = 0; j < p; ++j)
			if (std::fabs(b_tmp(j,k)) > 0.0) ++nvar;
			
		if (crit == "aic") {
			criteria[k] = n * std::log(mse) + 2 * nvar;
		}
		else if (crit == "bic") {
			criteria[k] = n * std::log(mse) + nvar * std::log(n);
		}
		else if (crit == "ebic") {
			criteria[k] = n * log(mse) + nvar * log(n) + (2 * nvar * 0.5 * log(p));
		}
		else {
			stop("Not a valid information criterion. Use: aic, bic, ebic.");																						   
		}
		if (!lambda_max_found && nvar > 0) {
			lambda_max_found = true;
			double new_lambda_max = lambda[k-1];
   
			double new_minval = 0.0;
			if (flg_lambda_min && lambda_min_input == 0.0) {
				new_minval = c_lambda_min * (std::log(p) / n);			   
			} else if (flg_lambda_min && lambda_min_input != 0.0) {
				stop("Conflicting inputs: flg_lambda_min = TRUE and lambda_min > 0 are not allowed.");
			} else if (lambda_min_input != 0.0) {
				new_minval = lambda_min_input;
			} else {
				new_minval = epsilon * new_lambda_max;
			}
			double logmax = std::log(new_lambda_max);
			double logmin = std::log(new_minval);
			for (int h = 0; h < nlambda; ++h) 
				lambda[h] = std::exp(logmax - ((logmax - logmin) * h / (nlambda - 1)));

			b_tmp = NumericMatrix(p, nlambda);
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
	
	for (int j = 0; j < p; ++j) beta[j] = b_tmp(j, best_lambda);

	// a0 + beta
	double a0 = ym;
	NumericVector coef(p+1);					  
	for (int j = 0; j < p; ++j) {
		a0 -= xm[j] * beta[j];
		coef[j+1] = beta[j];
	}
	coef[0] = a0;
	
    return List::create(
		_["a0"] = a0,
        _["beta"] = beta,
		_["coefficients"] = coef,
        _["lambda.opt"] = lambda[best_lambda],
        _["best.lambda.index"] = best_lambda,
        _["lambda.sequence"] = lambda,
		_["lambda.max"] = lambda_max,
		_["lambda.min"] = minval,
        _["beta.path"] = b_tmp,
		_["n.iterations"]=vec_max_itt,
		_["beta.init"] = b_init
    );
}
