#include <actionetcore.h>

#include <GraphHelper.h>
#include <RBConfigurationVertexPartition.h>
#include <CPMVertexPartition.h>
#include <Optimiser.h>

namespace ACTIONetcore {
	vec signed_cluster(sp_mat A, double resolution_parameter) {
		int nV = A.n_rows;
		int nE = A.n_nonzero;
		
		vec clusters = zeros(nV);
		
		printf("[signed] Clustering graph with %d vertices and %d edges\n", nV, nE);
				
		igraph_t g;
		igraph_vector_t v;
		igraph_vector_init(&v, 2*nE);
		vector<double> edge_weights(nE);


		int idx = 0;

		sp_mat::iterator it     = A.begin();
		sp_mat::iterator it_end = A.end();

		for(; it != it_end; ++it) {
			edge_weights[idx] = (*it);
			VECTOR(v)[2*idx] = it.row();
			VECTOR(v)[2*idx+1] = it.col();
			
			idx++;
		}
		igraph_create(&g, &v, nV, 0);

		Graph *G = new Graph(&g, edge_weights);
		//printf("idx = %d, V = %d, E = %d, total_weight = %f, weighted = %d, directed = %d\n", idx, G->vcount(), G->ecount(), G->total_weight(), G->is_weighted(), G->is_directed());
		

		CPMVertexPartition* partition = new CPMVertexPartition(G, resolution_parameter);

		Optimiser *opt = new Optimiser();
		opt->optimise_partition(partition);
	
		for(int i = 0; i < nV; i++) {
			clusters(i) = partition->membership(i)+1;
		}
		
	    partition->destructor_delete_graph = true;
		delete(G);
		igraph_vector_destroy(&v);
		igraph_destroy(&g);
		delete(opt);
		
		return(clusters);
	}
	
	vec unsigned_cluster(sp_mat A, double resolution_parameter) {
		int nV = A.n_rows;
		int nE = A.n_nonzero;
		
		vec clusters = zeros(nV);
		
		printf("[unsigned] Clustering graph with %d vertices and %d edges\n", nV, nE);
				
		igraph_t g;
		igraph_vector_t v;
		igraph_vector_init(&v, 2*nE);
		vector<double> edge_weights(nE);


		int idx = 0;

		sp_mat::iterator it     = A.begin();
		sp_mat::iterator it_end = A.end();

		for(; it != it_end; ++it) {
			edge_weights[idx] = (*it);
			VECTOR(v)[2*idx] = it.row();
			VECTOR(v)[2*idx+1] = it.col();
			
			idx++;
		}
		igraph_create(&g, &v, nV, 0);

		Graph *G = new Graph(&g, edge_weights);
		//printf("idx = %d, V = %d, E = %d, total_weight = %f, weighted = %d, directed = %d\n", idx, G->vcount(), G->ecount(), G->total_weight(), G->is_weighted(), G->is_directed());
		

		RBConfigurationVertexPartition*partition = new RBConfigurationVertexPartition(G, resolution_parameter);

		Optimiser *opt = new Optimiser();
		opt->optimise_partition(partition);
	
		for(int i = 0; i < nV; i++) {
			clusters(i) = partition->membership(i)+1;
		}
		
	    partition->destructor_delete_graph = true;
		delete(G);
		igraph_vector_destroy(&v);
		igraph_destroy(&g);
		delete(opt);
		
		return(clusters);
	}
	
}
