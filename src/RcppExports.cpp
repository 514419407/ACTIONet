// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/ACTIONet.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// reduceGeneExpression
List reduceGeneExpression(sp_mat expression, int reduced_dim, int method, int iters);
RcppExport SEXP _ACTIONet_reduceGeneExpression(SEXP expressionSEXP, SEXP reduced_dimSEXP, SEXP methodSEXP, SEXP itersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< sp_mat >::type expression(expressionSEXP);
    Rcpp::traits::input_parameter< int >::type reduced_dim(reduced_dimSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< int >::type iters(itersSEXP);
    rcpp_result_gen = Rcpp::wrap(reduceGeneExpression(expression, reduced_dim, method, iters));
    return rcpp_result_gen;
END_RCPP
}
// runsimplexRegression
mat runsimplexRegression(mat A, mat B);
RcppExport SEXP _ACTIONet_runsimplexRegression(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< mat >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(runsimplexRegression(A, B));
    return rcpp_result_gen;
END_RCPP
}
// runSPA
mat runSPA(mat A, int k);
RcppExport SEXP _ACTIONet_runSPA(SEXP ASEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(runSPA(A, k));
    return rcpp_result_gen;
END_RCPP
}
// runAA
List runAA(mat A, mat W0);
RcppExport SEXP _ACTIONet_runAA(SEXP ASEXP, SEXP W0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< mat >::type W0(W0SEXP);
    rcpp_result_gen = Rcpp::wrap(runAA(A, W0));
    return rcpp_result_gen;
END_RCPP
}
// runACTION
List runACTION(mat S_r, int k_min, int k_max, int thread_no);
RcppExport SEXP _ACTIONet_runACTION(SEXP S_rSEXP, SEXP k_minSEXP, SEXP k_maxSEXP, SEXP thread_noSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type S_r(S_rSEXP);
    Rcpp::traits::input_parameter< int >::type k_min(k_minSEXP);
    Rcpp::traits::input_parameter< int >::type k_max(k_maxSEXP);
    Rcpp::traits::input_parameter< int >::type thread_no(thread_noSEXP);
    rcpp_result_gen = Rcpp::wrap(runACTION(S_r, k_min, k_max, thread_no));
    return rcpp_result_gen;
END_RCPP
}
// computeFullDist
mat computeFullDist(mat& H_stacked, int thread_no, int verbose);
RcppExport SEXP _ACTIONet_computeFullDist(SEXP H_stackedSEXP, SEXP thread_noSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type H_stacked(H_stackedSEXP);
    Rcpp::traits::input_parameter< int >::type thread_no(thread_noSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(computeFullDist(H_stacked, thread_no, verbose));
    return rcpp_result_gen;
END_RCPP
}
// computeNearestDist
sp_mat computeNearestDist(mat& H_stacked, double kNN, int thread_no);
RcppExport SEXP _ACTIONet_computeNearestDist(SEXP H_stackedSEXP, SEXP kNNSEXP, SEXP thread_noSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type H_stacked(H_stackedSEXP);
    Rcpp::traits::input_parameter< double >::type kNN(kNNSEXP);
    Rcpp::traits::input_parameter< int >::type thread_no(thread_noSEXP);
    rcpp_result_gen = Rcpp::wrap(computeNearestDist(H_stacked, kNN, thread_no));
    return rcpp_result_gen;
END_RCPP
}
// smoothKNN
sp_mat smoothKNN(sp_mat D, int thread_no);
RcppExport SEXP _ACTIONet_smoothKNN(SEXP DSEXP, SEXP thread_noSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< sp_mat >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type thread_no(thread_noSEXP);
    rcpp_result_gen = Rcpp::wrap(smoothKNN(D, thread_no));
    return rcpp_result_gen;
END_RCPP
}
// computeNearestDist_edgeList
List computeNearestDist_edgeList(mat& H_stacked, double kNN, int thread_no);
RcppExport SEXP _ACTIONet_computeNearestDist_edgeList(SEXP H_stackedSEXP, SEXP kNNSEXP, SEXP thread_noSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type H_stacked(H_stackedSEXP);
    Rcpp::traits::input_parameter< double >::type kNN(kNNSEXP);
    Rcpp::traits::input_parameter< int >::type thread_no(thread_noSEXP);
    rcpp_result_gen = Rcpp::wrap(computeNearestDist_edgeList(H_stacked, kNN, thread_no));
    return rcpp_result_gen;
END_RCPP
}
// buildACTIONet
List buildACTIONet(mat& H_stacked, int kNN, int thread_no);
RcppExport SEXP _ACTIONet_buildACTIONet(SEXP H_stackedSEXP, SEXP kNNSEXP, SEXP thread_noSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type H_stacked(H_stackedSEXP);
    Rcpp::traits::input_parameter< int >::type kNN(kNNSEXP);
    Rcpp::traits::input_parameter< int >::type thread_no(thread_noSEXP);
    rcpp_result_gen = Rcpp::wrap(buildACTIONet(H_stacked, kNN, thread_no));
    return rcpp_result_gen;
END_RCPP
}
// buildAdaptiveACTIONet
List buildAdaptiveACTIONet(mat& H_stacked, double LC, double epsilon, int thread_no);
RcppExport SEXP _ACTIONet_buildAdaptiveACTIONet(SEXP H_stackedSEXP, SEXP LCSEXP, SEXP epsilonSEXP, SEXP thread_noSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type H_stacked(H_stackedSEXP);
    Rcpp::traits::input_parameter< double >::type LC(LCSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< int >::type thread_no(thread_noSEXP);
    rcpp_result_gen = Rcpp::wrap(buildAdaptiveACTIONet(H_stacked, LC, epsilon, thread_no));
    return rcpp_result_gen;
END_RCPP
}
// layoutACTIONet
List layoutACTIONet(sp_mat G, mat S_r, int compactness_level, unsigned int n_epochs, int thread_no);
RcppExport SEXP _ACTIONet_layoutACTIONet(SEXP GSEXP, SEXP S_rSEXP, SEXP compactness_levelSEXP, SEXP n_epochsSEXP, SEXP thread_noSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< sp_mat >::type G(GSEXP);
    Rcpp::traits::input_parameter< mat >::type S_r(S_rSEXP);
    Rcpp::traits::input_parameter< int >::type compactness_level(compactness_levelSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type n_epochs(n_epochsSEXP);
    Rcpp::traits::input_parameter< int >::type thread_no(thread_noSEXP);
    rcpp_result_gen = Rcpp::wrap(layoutACTIONet(G, S_r, compactness_level, n_epochs, thread_no));
    return rcpp_result_gen;
END_RCPP
}
// reconstructArchetypes
List reconstructArchetypes(sp_mat& S, const List& C_trace, const List& H_trace, double z_threshold);
RcppExport SEXP _ACTIONet_reconstructArchetypes(SEXP SSEXP, SEXP C_traceSEXP, SEXP H_traceSEXP, SEXP z_thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< sp_mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const List& >::type C_trace(C_traceSEXP);
    Rcpp::traits::input_parameter< const List& >::type H_trace(H_traceSEXP);
    Rcpp::traits::input_parameter< double >::type z_threshold(z_thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(reconstructArchetypes(S, C_trace, H_trace, z_threshold));
    return rcpp_result_gen;
END_RCPP
}
// extractArchetypeAssociatedSamples
mat extractArchetypeAssociatedSamples(sp_mat& G, mat& H_stacked, double alpha);
RcppExport SEXP _ACTIONet_extractArchetypeAssociatedSamples(SEXP GSEXP, SEXP H_stackedSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< sp_mat& >::type G(GSEXP);
    Rcpp::traits::input_parameter< mat& >::type H_stacked(H_stackedSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(extractArchetypeAssociatedSamples(G, H_stacked, alpha));
    return rcpp_result_gen;
END_RCPP
}
// assessFeatureSets
mat assessFeatureSets(sp_mat& S, List index_sets, int rand_perm);
RcppExport SEXP _ACTIONet_assessFeatureSets(SEXP SSEXP, SEXP index_setsSEXP, SEXP rand_permSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< sp_mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< List >::type index_sets(index_setsSEXP);
    Rcpp::traits::input_parameter< int >::type rand_perm(rand_permSEXP);
    rcpp_result_gen = Rcpp::wrap(assessFeatureSets(S, index_sets, rand_perm));
    return rcpp_result_gen;
END_RCPP
}
// assessFeatureSets_archs
mat assessFeatureSets_archs(mat& archetype_profile, List index_sets, int rand_perm);
RcppExport SEXP _ACTIONet_assessFeatureSets_archs(SEXP archetype_profileSEXP, SEXP index_setsSEXP, SEXP rand_permSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type archetype_profile(archetype_profileSEXP);
    Rcpp::traits::input_parameter< List >::type index_sets(index_setsSEXP);
    Rcpp::traits::input_parameter< int >::type rand_perm(rand_permSEXP);
    rcpp_result_gen = Rcpp::wrap(assessFeatureSets_archs(archetype_profile, index_sets, rand_perm));
    return rcpp_result_gen;
END_RCPP
}
// assessFeatureSets_decoupled
mat assessFeatureSets_decoupled(mat archetype_profile, mat H_stacked, List index_sets, int rand_perm);
RcppExport SEXP _ACTIONet_assessFeatureSets_decoupled(SEXP archetype_profileSEXP, SEXP H_stackedSEXP, SEXP index_setsSEXP, SEXP rand_permSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type archetype_profile(archetype_profileSEXP);
    Rcpp::traits::input_parameter< mat >::type H_stacked(H_stackedSEXP);
    Rcpp::traits::input_parameter< List >::type index_sets(index_setsSEXP);
    Rcpp::traits::input_parameter< int >::type rand_perm(rand_permSEXP);
    rcpp_result_gen = Rcpp::wrap(assessFeatureSets_decoupled(archetype_profile, H_stacked, index_sets, rand_perm));
    return rcpp_result_gen;
END_RCPP
}
// computeAutocorrelation
List computeAutocorrelation(sp_mat& G, mat& scores, int rand_perm, int num_shuffles);
RcppExport SEXP _ACTIONet_computeAutocorrelation(SEXP GSEXP, SEXP scoresSEXP, SEXP rand_permSEXP, SEXP num_shufflesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< sp_mat& >::type G(GSEXP);
    Rcpp::traits::input_parameter< mat& >::type scores(scoresSEXP);
    Rcpp::traits::input_parameter< int >::type rand_perm(rand_permSEXP);
    Rcpp::traits::input_parameter< int >::type num_shuffles(num_shufflesSEXP);
    rcpp_result_gen = Rcpp::wrap(computeAutocorrelation(G, scores, rand_perm, num_shuffles));
    return rcpp_result_gen;
END_RCPP
}
// phenotypeEnrichment
mat phenotypeEnrichment(mat& H_stacked, mat& phenotype_associations, int rand_perm_no);
RcppExport SEXP _ACTIONet_phenotypeEnrichment(SEXP H_stackedSEXP, SEXP phenotype_associationsSEXP, SEXP rand_perm_noSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type H_stacked(H_stackedSEXP);
    Rcpp::traits::input_parameter< mat& >::type phenotype_associations(phenotype_associationsSEXP);
    Rcpp::traits::input_parameter< int >::type rand_perm_no(rand_perm_noSEXP);
    rcpp_result_gen = Rcpp::wrap(phenotypeEnrichment(H_stacked, phenotype_associations, rand_perm_no));
    return rcpp_result_gen;
END_RCPP
}
// MWM
mat MWM(mat& G);
RcppExport SEXP _ACTIONet_MWM(SEXP GSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type G(GSEXP);
    rcpp_result_gen = Rcpp::wrap(MWM(G));
    return rcpp_result_gen;
END_RCPP
}
// batchPR
mat batchPR(sp_mat& G, mat U, double alpha, int thread_no, double tol);
RcppExport SEXP _ACTIONet_batchPR(SEXP GSEXP, SEXP USEXP, SEXP alphaSEXP, SEXP thread_noSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< sp_mat& >::type G(GSEXP);
    Rcpp::traits::input_parameter< mat >::type U(USEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type thread_no(thread_noSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(batchPR(G, U, alpha, thread_no, tol));
    return rcpp_result_gen;
END_RCPP
}
// sweepcut
vec sweepcut(sp_mat A, vec s);
RcppExport SEXP _ACTIONet_sweepcut(SEXP ASEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< sp_mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< vec >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(sweepcut(A, s));
    return rcpp_result_gen;
END_RCPP
}
// mergeArchetypes
sp_mat mergeArchetypes(mat C_stacked, mat H_stacked);
RcppExport SEXP _ACTIONet_mergeArchetypes(SEXP C_stackedSEXP, SEXP H_stackedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type C_stacked(C_stackedSEXP);
    Rcpp::traits::input_parameter< mat >::type H_stacked(H_stackedSEXP);
    rcpp_result_gen = Rcpp::wrap(mergeArchetypes(C_stacked, H_stacked));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ACTIONet_reduceGeneExpression", (DL_FUNC) &_ACTIONet_reduceGeneExpression, 4},
    {"_ACTIONet_runsimplexRegression", (DL_FUNC) &_ACTIONet_runsimplexRegression, 2},
    {"_ACTIONet_runSPA", (DL_FUNC) &_ACTIONet_runSPA, 2},
    {"_ACTIONet_runAA", (DL_FUNC) &_ACTIONet_runAA, 2},
    {"_ACTIONet_runACTION", (DL_FUNC) &_ACTIONet_runACTION, 4},
    {"_ACTIONet_computeFullDist", (DL_FUNC) &_ACTIONet_computeFullDist, 3},
    {"_ACTIONet_computeNearestDist", (DL_FUNC) &_ACTIONet_computeNearestDist, 3},
    {"_ACTIONet_smoothKNN", (DL_FUNC) &_ACTIONet_smoothKNN, 2},
    {"_ACTIONet_computeNearestDist_edgeList", (DL_FUNC) &_ACTIONet_computeNearestDist_edgeList, 3},
    {"_ACTIONet_buildACTIONet", (DL_FUNC) &_ACTIONet_buildACTIONet, 3},
    {"_ACTIONet_buildAdaptiveACTIONet", (DL_FUNC) &_ACTIONet_buildAdaptiveACTIONet, 4},
    {"_ACTIONet_layoutACTIONet", (DL_FUNC) &_ACTIONet_layoutACTIONet, 5},
    {"_ACTIONet_reconstructArchetypes", (DL_FUNC) &_ACTIONet_reconstructArchetypes, 4},
    {"_ACTIONet_extractArchetypeAssociatedSamples", (DL_FUNC) &_ACTIONet_extractArchetypeAssociatedSamples, 3},
    {"_ACTIONet_assessFeatureSets", (DL_FUNC) &_ACTIONet_assessFeatureSets, 3},
    {"_ACTIONet_assessFeatureSets_archs", (DL_FUNC) &_ACTIONet_assessFeatureSets_archs, 3},
    {"_ACTIONet_assessFeatureSets_decoupled", (DL_FUNC) &_ACTIONet_assessFeatureSets_decoupled, 4},
    {"_ACTIONet_computeAutocorrelation", (DL_FUNC) &_ACTIONet_computeAutocorrelation, 4},
    {"_ACTIONet_phenotypeEnrichment", (DL_FUNC) &_ACTIONet_phenotypeEnrichment, 3},
    {"_ACTIONet_MWM", (DL_FUNC) &_ACTIONet_MWM, 1},
    {"_ACTIONet_batchPR", (DL_FUNC) &_ACTIONet_batchPR, 5},
    {"_ACTIONet_sweepcut", (DL_FUNC) &_ACTIONet_sweepcut, 2},
    {"_ACTIONet_mergeArchetypes", (DL_FUNC) &_ACTIONet_mergeArchetypes, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_ACTIONet(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
