#include <actionetcore.h>
#include <vptree.h>

namespace ACTIONetcore {
	mat computeFullDist(mat &H_stacked, int thread_no=-1, int verbose = 1) {	
		if(verbose > 0)
			printf("Building full distance matrix (returns full distance matrix)\n");
			
		int sample_no = H_stacked.n_cols;
		int archetype_no = H_stacked.n_rows;;

		mat D(sample_no, sample_no);
		
		H_stacked.transform( [](double val) { return (val < 0?0:val); } );
		mat H_stacked_norm = normalise(H_stacked, 1);		
		
		mat logProb = log(H_stacked_norm);
		logProb.transform( [](double val) { return (isinf(val)?0:val); } );
		
		vec entropies = -trans(sum(H_stacked_norm % logProb, 0));
				
		double scaler = 1.0 / (2.0*log(2.0));

		int perc = 1;
		int total_counts = 0;
		#pragma omp parallel for num_threads(thread_no) 			
		for(int v = 0; v < sample_no; v++) {
			total_counts ++;
			if(round(100*(double)total_counts / sample_no) > perc) {
				if(verbose > 0)
					printf("%d %%\n", perc);
				perc++;
			}			
			
			vec p = H_stacked_norm.col(v);					
			mat logM(H_stacked_norm.n_rows, sample_no);
			for(int c = 0; c < sample_no; c++) {
				logM.col(c) = log(0.5*(p + H_stacked_norm.col(c)));
			}		
			logM.transform( [](double val) { return ((isinf(val) || isnan(val))?0:val); } );
					
			vec DpM = trans(-p.t()*logM - entropies(v));
			vec DQM = trans(-sum(H_stacked_norm % logM, 0)) - entropies;
			
			vec JS_div = scaler*(DpM + DQM);
			JS_div(v) = 0;
		
			D.col(v) = sqrt(JS_div); // Sqrt of JS Div, not JS Div, is a metric
		}		
		D.transform( [](double val) { return (isnan(val)?0:val); } );
		D.transform( [](double val) { val = val < 0? 0:val; val = 1 < val? 1:val; return (val); } ); // Make sure that scores are normalized properly
		
		if(verbose > 0)
			printf("done\n");
			
		return D;
	}	
	
	
	sp_mat computeNearestDist(mat &H_stacked, int kNN, int thread_no=-1) {	
		printf("Building distance matrix of the %d nearest neighbors of each node (returns sparse distance matrix)\n", kNN);
		
		double epsilon = 1e-10;
		int sample_no = H_stacked.n_cols;
		
		if(kNN >= sample_no || kNN <= 1)
			kNN = min(30, sample_no-1);
			
		umat subs(2, kNN*sample_no);
		vec vv(kNN*sample_no);

		H_stacked.transform( [](double val) { return (val < 0?0:val); } );

		mat H_stacked_norm = normalise(H_stacked, 1, 0); // make the norm (sum) of each column 1			
		int archetype_no = H_stacked_norm.n_rows;
		
		// build ball tree on set
		VpTree<DataPoint, JSDiv_sqrt_distance>* tree = new VpTree<DataPoint, JSDiv_sqrt_distance>();
		std::vector<DataPoint> samples(sample_no); //, DataPoint(archetype_no, -1, data));
		for (int i = 0; i < sample_no; i++) {
			samples[i] = DataPoint(archetype_no, i, H_stacked_norm.colptr(i));
			//(H_stacked_norm.col(i)).print("col");
		}
		tree -> create(samples, 0);
		
		
		int perc = 1;
		int total_counts = 1;
		#pragma omp parallel num_threads(thread_no) 
		{
			vector< vector<int> > ind_arr(sample_no, std::vector<int>(kNN+1));
			vector< vector<double> > dist_arr(sample_no, std::vector<double>(kNN+1));
			#pragma omp for
			for (int v = 0; v < sample_no; v++) {
				total_counts ++;
				if(round(100*(double)total_counts / sample_no) > perc) {
					printf("%d %%\n", perc);
					perc++;
				}

				
				tree -> search(samples[v], kNN+1, &ind_arr[v], &dist_arr[v]);
				
				int base = v*kNN;
				for(int i = 1; i <= kNN; i++) {
					double d = dist_arr[v][i]; // To match with old scores
					d = d < epsilon?epsilon:d;
						
					subs(0, base + i-1) = ind_arr[v][i]-1;
					subs(1, base + i-1) = v;
					vv(base + i-1) = d;			
				}				
			}
		}
		samples.clear();
		delete tree;
			
		sp_mat D(subs, vv, sample_no, sample_no);	

		return(D);
	}
	

	field<mat> computeNearestDist_edgeList(mat &H_stacked, int kNN, int thread_no=-1) {	
		printf("Building distance matrix of the %d nearest neighbors of each node (returns edge list)\n", kNN);
		
		double epsilon = 1e-10;
		int sample_no = H_stacked.n_cols;
		
		if(kNN >= sample_no || kNN <= 1)
			kNN = min(30, sample_no-1);
			
		umat subs(2, kNN*sample_no);
		vec vv(kNN*sample_no);

		H_stacked.transform( [](double val) { return (val < 0?0:val); } );
		mat H_stacked_norm = normalise(H_stacked, 1, 0); // make the norm (sum) of each column 1			
		
		int archetype_no = H_stacked_norm.n_rows;
		
		// build ball tree on set
		VpTree<DataPoint, JSDiv_sqrt_distance>* tree = new VpTree<DataPoint, JSDiv_sqrt_distance>();
		std::vector<DataPoint> samples(sample_no); //, DataPoint(archetype_no, -1, data));
		for (int i = 0; i < sample_no; i++) {
			samples[i] = DataPoint(archetype_no, i, H_stacked_norm.colptr(i));
			//(H_stacked_norm.col(i)).print("col");
		}
		tree -> create(samples, 0);
		
		
		mat idx = zeros(sample_no, kNN+1);
		mat dist = zeros(sample_no, kNN+1);

		int perc = 1;
		int total_counts = 0;
		#pragma omp parallel num_threads(thread_no) 			
		{
			vector< vector<int> > ind_arr(sample_no, std::vector<int>(kNN+1));
			vector< vector<double> > dist_arr(sample_no, std::vector<double>(kNN+1));
			#pragma omp for
			for (int v = 0; v < sample_no; v++) {
				total_counts ++;
				if(round(100*(double)total_counts / sample_no) > perc) {
					printf("%d %%\n", perc);
					perc++;
				}
				
				tree -> search(samples[v], kNN+1, &ind_arr[v], &dist_arr[v]);
				idx(v, 0) = v+1;
				dist(v, 0) = 0;
							
				int base = v*kNN;
				for(int i = 1; i <= kNN; i++) {
					double d = dist_arr[v][i]; 
					d = d < epsilon?epsilon:d;
					
					#pragma omp atomic write 
					idx(v, i) = ind_arr[v][i];

					#pragma omp atomic write 
					dist(v, i) = d;
				}				
			}
		}			
		
		samples.clear();
		delete tree;
			
		field<mat> output(2);
		output(0) = idx;
		output(1) = dist;
		
		return(output);
	}


	
	field<sp_mat> buildACTIONet(mat &H_stacked, int kNN, int thread_no=-1) {	
		printf("Building ACTIONet\n");
		int sample_no = H_stacked.n_cols;		
		if(kNN >= sample_no || kNN < 1)
			kNN = min(30, sample_no-1);

		sp_mat D = computeNearestDist(H_stacked, kNN, thread_no);

		sp_mat G = D;
		for(sp_mat::iterator it = G.begin(); it != G.end(); ++it) {
		  (*it) = 1.0 - (*it);
		}					
		
		
		field<sp_mat> output(2);
		output(0) = sqrt(G % trans(G));
		output(1) = G;
		
		return(output);	
	}


	// k^{*}-Nearest Neighbors: From Global to Local (NIPS 2016)
	field<sp_mat> buildAdaptiveACTIONet(mat &H_stacked, double LC = 1.0, int thread_no=-1)	{

		printf("Building adaptive ACTIONet\n");
		
		
		double kappa = 3.0;
		int sample_no = H_stacked.n_cols;		
		int kNN = min(sample_no-1, (int)(kappa*round(sqrt(sample_no)))); // start with uniform k=sqrt(N) ["Pattern Classification" book by Duda et al.]

		field<mat> res = computeNearestDist_edgeList(H_stacked, kNN, thread_no);
		printf("Done precomputing distances\n");fflush(stdout);
		
		mat idx = res(0);
		mat dist = res(1);
		
		mat beta = LC*dist;
		vec beta_sum = zeros(sample_no);
		vec beta_sq_sum = zeros(sample_no);
		
		printf("Computing Lambda\n");
		mat lambda(size(beta));
		register int k;
		for(k = 1; k <= kNN; k++) {
			beta_sum += beta.col(k);
			beta_sq_sum += square(beta.col(k));
			
			lambda.col(k) = (1/(double)k) * ( beta_sum + sqrt(k + square(beta_sum) - k*beta_sq_sum) );
		}
		printf("done\n");fflush(stdout);

		sp_mat D(sample_no, sample_no);
		//# pragma omp parallel for num_threads(thread_no) 
		printf("Computing D\n");
		lambda = trans(lambda);
		vec node_lambda = zeros(sample_no);
		beta = trans(beta);
		for(int v = 0; v < sample_no; v++) {
			vec delta = lambda.col(v) - beta.col(v);			

			//printf("sum = %.2f, len = %d, idx = %d\n", sum(delta < 0), shift_idx.n);
			uvec shift_idx = find(delta <= 0, 1, "first");
			int neighbor_no = shift_idx.n_elem == 0?kNN:(shift_idx(0)-1);
			
			node_lambda(v) = lambda(neighbor_no, v);
			
			rowvec v_idx  = idx.row(v);
			rowvec v_dist = dist.row(v); // Shall we use lambda instead?			
			//vec v_dist = lambda(neighbor_no) - beta(span(0, neighbor_no), v);
			//v_dist /= max(v_dist);
			for (int i = 1; i <= neighbor_no; i++) {				
				int src = v_idx(i)-1;
				int dst = v;
				if(src > D.n_rows-1 || src < 0) {
					printf("v = %d, i = %d, src = %d, nrows = %d\n", v, i, src, D.n_rows);
					continue;
				}
				if(dst > D.n_cols-1 || dst < 0) {
					printf("v = %d, i = %d, dst = %d, ncols = %d\n", v, i, dst, D.n_rows);
					continue;
				}
				D(src, dst) = v_dist(i);
			}
		}
		printf("done\n");fflush(stdout);
		
		
		/*		
		sp_mat G = D;
		double epsilon = 1e-16;
		for(sp_mat::iterator it = G.begin(); it != G.end(); ++it) {
			(*it) = std::max(epsilon, node_lambda(it.col()) - LC*(*it));
		}			
		*/
		
		printf("Computing back to similarities\n");
		// Old method of computing similarities
		sp_mat G = D;
		for(sp_mat::iterator it = G.begin(); it != G.end(); ++it) {
		  (*it) = 1.0 - (*it);
		}					
		printf("done\n"); fflush(stdout);
		
		
		field<sp_mat> output(2);
		output(0) = sqrt(G % trans(G));
		output(1) = normalise(G, 1, 0);
		
		return(output);		
	}


	sp_mat smoothKNN(sp_mat D, int thread_no = -1) {		
		double epsilon = 1e-6;

		int nV = D.n_rows;
		sp_mat G = D;			
		
		#pragma omp parallel for num_threads(thread_no) 			
		for(int i = 0; i < nV; i++) {
			sp_mat v = D.col(i);
			
			vec vals = nonzeros(v);	
			if(vals.n_elem == 0)
				continue;
				
			double rho = min(vals);
			vec negated_shifted_vals = -(vals - rho); 
			double target = log2(vals.n_elem);
			
			// Binary search to find optimal sigma
			double sigma = 1.0;
			double lo = 0.0;
			double hi = DBL_MAX;
			
			int j;
			for(j = 0; j < 64; j ++) {
				double obj = sum(exp(negated_shifted_vals / sigma));

				if (abs(obj - target) < epsilon) {
					break;
				}

				if (target < obj) {
					hi = sigma;
					sigma = 0.5 * (lo + hi);
				}
				else {
					lo = sigma;
					if (hi == DBL_MAX) {
						sigma *= 2;
					}
					else {
						sigma = 0.5 * (lo + hi);
					}
				}				
			}
			
			double obj = sum(exp(negated_shifted_vals / sigma));			
			//printf("%d- rho = %.3f, degree = %d, log2(k) = %.2e, sigma = %.2e, residual = %.2e, iters = %d\n", i, rho, vals.n_elem, target, sigma, abs(obj - target), j);
			
			for(sp_mat::col_iterator it = G.begin_col(i); it != G.end_col(i); ++it) {
				*it = max(1e-16, exp( -max(0.0, (*it) - rho ) / sigma ));
			}			
		}
		return(G);
	}
}






