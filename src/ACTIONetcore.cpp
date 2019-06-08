#include <actionetcore.h>
#include <RcppArmadillo.h>

#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

#define ARMA_USE_CXX11_RNG


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat computeFullDist(mat &H_stacked, int thread_no = 4, int verbose = 1) {  	
    mat D = ACTIONetcore::computeFullDist(H_stacked, thread_no, verbose);
    
    return D;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
sp_mat computeNearestDist(mat &H_stacked, double kNN, int thread_no = 4) {  
    sp_mat D = ACTIONetcore::computeNearestDist(H_stacked, kNN, thread_no);
    
    return D;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
sp_mat smoothKNN(sp_mat D, int thread_no = 4) {  
    sp_mat G = ACTIONetcore::smoothKNN(D, thread_no);
    
    return G;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List computeNearestDist_edgeList(mat &H_stacked, double kNN, int thread_no = 4) {  	
    field<mat> NN = ACTIONetcore::computeNearestDist_edgeList(H_stacked, kNN, thread_no);
    
	List out_list;		
	out_list["idx"] = NN(0);
	out_list["dist"] = NN(1);

    return out_list;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List buildACTIONet(mat &H_stacked, int kNN = 30, int thread_no = 8) {  
    field<sp_mat> res = ACTIONetcore::buildACTIONet(H_stacked, kNN, thread_no);

	List out_list;		
	out_list["ACTIONet"] = res(0);
	out_list["ACTIONet_asym"] = res(1);

    return out_list;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List buildAdaptiveACTIONet(mat &H_stacked, double LC = 1.0, double epsilon = 0.0, int thread_no = 8) {  
		
    field<sp_mat> res = ACTIONetcore::buildAdaptiveACTIONet(H_stacked, LC, epsilon, thread_no);

	List out_list;		
	out_list["ACTIONet"] = res(0);
	out_list["ACTIONet_asym"] = res(1);

    return out_list;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List layoutACTIONet(sp_mat G,
	mat S_r,
	int compactness_level=1.0,
	unsigned int n_epochs = 500,
	int thread_no = 8) {

	field<mat> res = ACTIONetcore::layoutACTIONet(G, S_r, compactness_level, n_epochs, thread_no);
    
	List out_list;		
	out_list["coordinates"] = res(0);
	out_list["coordinates_3D"] = res(1);
	out_list["colors"] = res(2);

    return out_list;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List reconstructArchetypes(sp_mat &S, const List& C_trace, const List& H_trace, double z_threshold = 1.0) {
	
	int n_list = H_trace.size();
	vector<mat> C_trace_vec(n_list+1);
	vector<mat> H_trace_vec(n_list+1);
	for (int i = 0; i < n_list; i++) {
		if(Rf_isNull(H_trace[i])) {
			continue;
		}
		C_trace_vec[i+1] = (as<mat>(C_trace[i]));
		H_trace_vec[i+1] = (as<mat>(H_trace[i]));
	}

	ACTIONetcore::multilevel_archetypal_decomposition results = ACTIONetcore::reconstructArchetypes(S, C_trace_vec, H_trace_vec, z_threshold = -1.0);

	List out_list;		
	
	
	for(int i = 0; i < results.selected_archs.n_elem; i++) results.selected_archs[i]++;
	out_list["selected_archs"] = results.selected_archs; 

	for(int i = 0; i < results.landmark_cells.n_elem; i++) results.landmark_cells[i]++;
	out_list["landmark_cells"] = results.landmark_cells; 

	out_list["C_stacked"] = results.C_stacked;
	out_list["H_stacked"] = results.H_stacked;
	
	out_list["archetype_profile"] = results.archetype_profile;

	out_list["backbone"] = results.backbone;

    return out_list;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat extractArchetypeAssociatedSamples(sp_mat &G, mat &H_stacked, double alpha = 0.85) {

	mat scores = ACTIONetcore::extractArchetypeAssociatedSamples(G, H_stacked, alpha);
	
	return scores;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat assessFeatureSets(sp_mat &S, List index_sets, int rand_perm = 100) {
	
	int n_list = index_sets.size();
	field<uvec> feature_sets(n_list);
	
	for (int i = 0; i < n_list; i++) {
		vec v = as<vec>(index_sets[i]);
		uvec u(v.n_elem);
		for(int j = 0; j < v.n_elem; j++) {
			u(j) = v(j) - 1;
		}
		feature_sets(i) = u;
	}
	
	mat scores = ACTIONetcore::assessFeatureSets(S, feature_sets, rand_perm);		

    return scores;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat assessFeatureSets_archs(mat &archetype_profile, List index_sets, int rand_perm = 100) {
	
	int n_list = index_sets.size();
	field<uvec> feature_sets(n_list);
	
	for (int i = 0; i < n_list; i++) {
		vec v = as<vec>(index_sets[i]);
		uvec u(v.n_elem);
		for(int j = 0; j < v.n_elem; j++) {
			u(j) = v(j) - 1;
		}
		feature_sets(i) = u;
	}
	
	mat scores = ACTIONetcore::assessFeatureSets_archs(archetype_profile, feature_sets, rand_perm);		

    return scores;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat assessFeatureSets_decoupled(mat archetype_profile, mat H_stacked, List index_sets, int rand_perm = 100) {
	
	int n_list = index_sets.size();
	field<uvec> feature_sets(n_list);
	
	for (int i = 0; i < n_list; i++) {
		vec v = as<vec>(index_sets[i]);
		uvec u(v.n_elem);
		for(int j = 0; j < v.n_elem; j++) {
			u(j) = v(j) - 1;
		}
		feature_sets(i) = u;
	}
	
	mat scores = ACTIONetcore::assessFeatureSets_decoupled(archetype_profile, H_stacked, feature_sets, rand_perm);	

    return scores;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List computeAutocorrelation (sp_mat &G, mat &scores, int rand_perm = 100, int num_shuffles = 10000) {
	
	field<vec> results = ACTIONetcore::computeAutocorrelation (G, scores, rand_perm, num_shuffles);
	

	List out_list;		
	out_list["C"] = results(0);
	out_list["mu"] = results(1);
	out_list["sigma"] = results(2);
	out_list["Z"] = results(3);
	
    return out_list;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat phenotypeEnrichment(mat &H_stacked, mat &phenotype_associations, int rand_perm_no) {	
	mat Z = ACTIONetcore::phenotypeEnrichment (H_stacked, phenotype_associations, rand_perm_no);

    return Z;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat MWM(mat &G) {	
	mat G_matched = ACTIONetcore::MWM(G);

    return G_matched;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat batchPR(sp_mat &G, mat U, double alpha = 0.85, int thread_no = 8, double tol = 1e-6) {
	
	mat U_smoothed = ACTIONetcore::batchPR(G, U, alpha, thread_no, tol);
	
    return U_smoothed;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec sweepcut(sp_mat A, vec s) {
    vec conductance = ACTIONetcore::sweepcut(A, s);

    return conductance;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
sp_mat mergeArchetypes(mat C_stacked, mat H_stacked) {
	sp_mat results = ACTIONetcore::mergeArchetypes(C_stacked, H_stacked);

    return results;	
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec signed_cluster(sp_mat A, double resolution_parameter = 1.0) {
	vec clusters = ACTIONetcore::signed_cluster(A, resolution_parameter);

    return clusters;	
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec unsigned_cluster(sp_mat A, double resolution_parameter = 1.0) {
	vec clusters = ACTIONetcore::unsigned_cluster(A, resolution_parameter);

    return clusters;	
}
