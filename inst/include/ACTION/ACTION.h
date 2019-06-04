#ifndef ACTION_H
#define ACTION_H

#include <omp.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <cmath>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <errno.h>
#include <sys/stat.h>
#include <getopt.h>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <unordered_map>
#include <ctime>


#include <RcppArmadillo.h>
using namespace arma;
using namespace std;

#define DEFAULT_PCA_DIM 50

// Expression reduction methods
#define PCA_only 0
#define ACTIONplusPCA 1

// Generic structure returned by all reduction methods
struct Projection {
	mat S_r;
	mat V;
	vec lambda;
	vec exp_var;
};


struct SPA_results {
	uvec selected_columns;
	vec column_norms;
};


struct ACTION_results {
	vector<uvec> selected_cols;
	vector<mat> H;
	vector<mat> C;
};

namespace ACTION {
	SPA_results SPA(mat M, int k); // Solve convex-NMF using the Successive Projection Algorithm (SPA)
	void simplexRegression(mat &A, mat &B, double *X_ptr); // min_{X} (|| AX - B ||) s.t. simplex constraint
	field<mat> AA (mat X, mat Z0); // Robust archetypal analysis method
	ACTION_results runACTION(mat S_r, int k_min, int k_max, int numThreads); // Main ACTION function	
	Projection reduceGeneExpression(sp_mat &expression, int reduced_dim, int method, int iter);
}

#endif
