#include <actionetcore.h>
mat zscore(mat A) {	
	rowvec mu = mean(A, 0);
	rowvec sigma = stddev(A, 0);

	
	for(int j = 0; j < A.n_cols; j++) {
		A.col(j) = (A.col(j) - mu(j)) / sigma(j);
	}
	
	return A;
}

namespace ACTIONetcore {	
	field<mat> layoutACTIONet(sp_mat &G,
		mat &S_r,
		int compactness_level = 50,
		unsigned int n_epochs = 500,
		int thread_no = 8) { 

		unsigned int nV = G.n_rows;
		
		mat initial_coordinates = S_r.rows(regspace<uvec>(0, 1));		
		
		// Convert back from similarity to distance, and then smooth them using the UMAP framework
		/*
		for(sp_mat::iterator it = G.begin(); it != G.end(); ++it) {
		  (*it) = 1.0 - (*it);
		}
		*/			
		printf("Running layout with: compactness=%d, # epochs = %d\n", compactness_level, n_epochs);
		
		G.for_each( [](sp_mat::elem_type& val) { val = 1.0 - val; } );
		
		G = smoothKNN(G, thread_no);


		if(compactness_level < 0 || compactness_level > 100)
			compactness_level = 50;
			
		double a_param = UMAP_A[compactness_level];
		double b_param = UMAP_B[compactness_level];

		// linearized list of edges (1-simplices)
		unsigned int nE = G.n_nonzero;
		vector<unsigned int> positive_head(nE);
		vector<unsigned int> positive_tail(nE);
		vector<double> epochs_per_sample(nE);		
		
		int i = 0;
		double w_max = max(max(G));
		for(sp_mat::iterator it = G.begin(); it != G.end(); ++ it) {
			epochs_per_sample[i] = w_max / (*it); // Higher the weight of the edge, the more likely it is to be sampled (inversely proportional to epochs_per_sample)
			positive_head[i] = it.row();
			positive_tail[i++] = it.col();
		}		
		
		// Initial coordinates of vertices (0-simplices)
		mat Ct = initial_coordinates;		
		vector<double> head_vec(Ct.memptr(), Ct.memptr()+Ct.n_elem);		
		vector<double> tail_vec(head_vec);		
		
		
		printf("Computing 2D layout ... "); fflush(stdout);		
		// Stores linearized coordinates [x1, y1, x2, y2, ...]
		vector<double> result;
		const apumap_gradient gradient(a_param, b_param, GAMMA);
		result = optimize_layout(
			gradient, head_vec, tail_vec, positive_head, positive_tail,
			n_epochs, nV, epochs_per_sample, LEARNING_RATE,
			NEGATIVE_SAMPLE_RATE, UMAP_SEED);	

		mat coordinates(result.data(), 2, nV);
		coordinates = zscore(trans(coordinates));		
		printf("Done\n"); fflush(stdout);

		
		/****************************
		 *  Now compute node colors *
		 ***************************/	
		Ct = trans(join_horiz(coordinates, zscore(trans(S_r.row(2)))));

		head_vec.clear(); 
		head_vec.resize(Ct.n_elem);
		std::copy(Ct.memptr(), Ct.memptr() + Ct.n_elem, head_vec.begin());
		
		tail_vec = head_vec;		
		
		printf("Compute 3D layout ... "); fflush(stdout);
		result.clear();
		result = optimize_layout(
			gradient, head_vec, tail_vec, positive_head, positive_tail,
			n_epochs, nV, epochs_per_sample, LEARNING_RATE,
			NEGATIVE_SAMPLE_RATE, UMAP_SEED);	

		mat coordinates_3D(result.data(), 3, nV);	
		coordinates_3D = zscore(trans(coordinates_3D));		
		printf("Done\n"); fflush(stdout);
	  
	  
		printf("Estimating node colors ... "); fflush(stdout);
		vec a = 128*coordinates_3D.col(0) / max(abs(coordinates_3D.col(0)));
		vec b = 128*coordinates_3D.col(1) / max(abs(coordinates_3D.col(1)));
		vec L = coordinates_3D.col(2);
		L = 25.0 + 50.0*(L - min(L)) / (max(L) - min(L));

		double r_channel, g_channel, b_channel;
		mat RGB_colors = zeros(nV, 3);
		for(int i = 0; i < nV; i++) {
			Lab2Rgb(&r_channel, &g_channel, &b_channel, L(i), a(i), b(i));			

			RGB_colors(i, 0) = min(1.0, max(0.0, r_channel));
			RGB_colors(i, 1) = min(1.0, max(0.0, g_channel));
			RGB_colors(i, 2) = min(1.0, max(0.0, b_channel));			
		}
		
		printf("done\n");
		
		field<mat> res(3);
		res(0) = coordinates;
		res(1) = coordinates_3D;
		res(2) = RGB_colors;
		
		return res;	  		
	}	
}





