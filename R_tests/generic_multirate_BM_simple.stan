functions {
  
  // BEAST of a function that calculates the likelihood of the trait data, marginalized across trees, given the current parameters
  real prune(matrix intra_cor, vector intra_var, matrix[] inter_cor, vector[] inter_var, // parameters
             int N, int K, int S, int T, int E, int C,
             int[] n_obs, int[] ind_obs, matrix Y,
             int[] n_inc, int[] pos_inc, vector r, int[] s,
             int[] n_des, int[] ind_des, int[] pos_des, int[,] ind_tip, int[] ind_root, int[] ind_prune,
             int[] n_dim, int[] ind_dim, int[] pos_dim, int[] code_obs, int[] code_tip, int[] code_nod){
               
               matrix[K, K] intra_cov; // intra-tip variance-covariance matrix
               matrix[K, K] code_intra_P[C]; // intra-tip precision matrix for each observation code
               matrix[K, K] inter_cov[S]; // inter-tip variance covariance matrices for each regime (AKA rate matrices)
               matrix[K, E] X; // expected trait values for each node given descendants (i.e., "expectations")
               matrix[K, K] P[E]; // precision matrix of expectations for each node  (i.e., "precisions")
               vector[E] R; // log scaling constants for partial likelihoods of expectations for each node (i.e., "scalars")
               vector[T] L;
               real max_L;
               int obs_counter; // counter to keep track of position along raw observation indices vector (ind_obs)
               
               // INITIALIZE TIPS
               // TODO: add ability to specify a priori intra-tip variance-covariance matrices for certain tips
               // TODO: you could even allow a priori specification for specific traits, as long as you assume no correlation with estimated variance-covariance matrices, I think...
               // TODO: otherwise impossible correlation structures could be sampled...though maybe it wouldn't be too bad an issue?
               // form intra-tip variance-covariance matrix
               intra_cov = quad_form_diag(intra_cor, sqrt(intra_var));
               // invert subsets of intra-tip variance-covariance matrix according to observation codes (avoids redundant inversions in loop below)
               for(i in 1:C){
                 int tmp_pos;
                 int tmp_k;
                 tmp_pos = pos_dim[i];
                 tmp_k = n_dim[i];
                 code_intra_P[i] = rep_matrix(0, K, K);
                 if(tmp_k){
                   int tmp_k_ind[tmp_k];
                   tmp_k_ind = segment(ind_dim, tmp_pos, tmp_k);
                   code_intra_P[i, tmp_k_ind, tmp_k_ind] = inverse_spd(intra_cov[tmp_k_ind, tmp_k_ind]);
                 }
               }
               // get expectation, precision, and scalar for each tip
               obs_counter = 1;
               for(i in 1:N){
                 // get indices and initialize output quantities
                 int tmp_n;
                 int ind_nod[T];
                 matrix[K, K] cur_P; // current precision for ENTIRE tip
                 vector[K] cur_X; // current expectation for ENTIRE tip
                 real cur_R; // current scalar for ENTIRE tip
                 vector[K] sum_PY; // sum of raw observations, weighted by their precisions
                 real sum_k; // involved in scalar calculations
                 tmp_n = n_obs[i];
                 ind_nod = ind_tip[, i];
                 cur_P = rep_matrix(0, K, K);
                 cur_X = rep_vector(0, K);
                 cur_R = 0;
                 sum_PY = rep_vector(0, K);
                 sum_k = 0;
                 if(tmp_n){
                   int tmp_n_ind[tmp_n];
                   tmp_n_ind = segment(ind_obs, obs_counter, tmp_n);
                   // do necessary calculations for each raw observation
                   for(j in tmp_n_ind){
                     int tmp_code;
                     int tmp_k_pos;
                     int tmp_k;
                     tmp_code = code_obs[j];
                     tmp_k_pos = pos_dim[tmp_code];
                     tmp_k = n_dim[tmp_code];
                     { // raw observations always have at least 1 informative dimension
                       int tmp_k_ind[tmp_k];
                       matrix[tmp_k, tmp_k] tmp_P; // temporary precision for raw observation
                       vector[tmp_k] tmp_Y; // temporary observed trait values
                       vector[tmp_k] tmp_PY; // temporary precision right multiplied by temporary observed trait values
                       tmp_k_ind = segment(ind_dim, tmp_k_pos, tmp_k);
                       tmp_P = code_intra_P[tmp_code, tmp_k_ind, tmp_k_ind];
                       cur_P[tmp_k_ind, tmp_k_ind] += tmp_P;
                       tmp_Y = Y[tmp_k_ind, j];
                       tmp_PY = tmp_P * tmp_Y;
                       sum_PY[tmp_k_ind] += tmp_PY;
                       sum_k -= tmp_k;
                       cur_R += log_determinant(tmp_P) - tmp_Y' * tmp_PY;
                     }
                   }
                   // do necessary calculations for ENTIRE tip
                   {
                     int tmp_code;
                     int tmp_k_pos;
                     int tmp_k;
                     tmp_code = code_tip[i];
                     tmp_k_pos = pos_dim[tmp_code];
                     tmp_k = n_dim[tmp_code];
                     {
                       int tmp_k_ind[tmp_k];
                       matrix[tmp_k, tmp_k] tmp_P; // precision for ENTIRE tip, subsetted to informative dimensions
                       vector[tmp_k] tmp_X; // expectation for ENTIRE tip, subsetted to informative dimensions
                       tmp_k_ind = segment(ind_dim, tmp_k_pos, tmp_k);
                       tmp_P = cur_P[tmp_k_ind, tmp_k_ind];
                       tmp_X = mdivide_left_spd(tmp_P, sum_PY[tmp_k_ind]);
                       cur_X[tmp_k_ind] = tmp_X;
                       sum_k += tmp_k;
                       cur_R += tmp_X' * tmp_P * tmp_X - log_determinant(tmp_P);
                     }
                   }
                 }
                 // insert output quantities into full expectation, precision, and scalar arrays
                 X[, ind_nod] = rep_matrix(cur_X, T);
                 P[ind_nod] = rep_array(cur_P, T);
                 R[ind_nod] = rep_vector(0.5 * (cur_R + sum_k * log(2 * pi())), T);
                 // increment observation counter
                 obs_counter += tmp_n;
               }
               
               // PRUNING ALGORITHM
               // form rate matrices
               for(i in 1:S){
                 inter_cov[i] = quad_form_diag(inter_cor[i], sqrt(inter_var[i]));
               }
               for(i in ind_prune){
                 // get indices and initialize output quantities
                 int tmp_n;
                 int tmp_n_pos;
                 matrix[K, K] cur_P; // current precision for focal node
                 vector[K] cur_X; // current expectation for focal node
                 real cur_R; // current scalar for focal node
                 vector[K] sum_PX; // sum of descendant expectations, weighted by their "propagated" precisions (see below)
                 real sum_k; // involved in scalar calculations
                 tmp_n = n_des[i];
                 tmp_n_pos = pos_des[i];
                 cur_P = rep_matrix(0, K, K);
                 cur_X = rep_vector(0, K);
                 sum_PX = rep_vector(0, K);
                 sum_k = 0;
                 { // every node in pruning sequence has at least 1 descendant
                   int tmp_n_ind[tmp_n];
                   tmp_n_ind = segment(ind_des, tmp_n_pos, tmp_n);
                   // initialize scalar with sum of descendant scalars
                   cur_R = 2 * sum(R[tmp_n_ind]);
                   // do necessary calculations for each descendant
                   // involves "propagating" precisions down each descendant edge to account for trait change along edges
                   for(j in tmp_n_ind){
                     int tmp_code;
                     int tmp_k_pos;
                     int tmp_k;
                     int tmp_m_pos;
                     int tmp_m;
                     tmp_code = code_nod[j];
                     tmp_k_pos = pos_dim[tmp_code];
                     tmp_k = n_dim[tmp_code];
                     tmp_m_pos = pos_inc[j];
                     tmp_m = n_inc[j];
                     if(tmp_k){
                       int tmp_k_ind[tmp_k];
                       vector[tmp_m] tmp_r;
                       int tmp_s[tmp_m];
                       vector[S] sum_r; // summed rate scalars multiplied by time durations in each state for a given edge
                       matrix[tmp_k, tmp_k] tmp_V; // temporary variance-covariance matrix of trait change along desendant edge
                       matrix[tmp_k, tmp_k] tmp_P; // temporary post-propagation precision for expectation
                       vector[tmp_k] tmp_X; // temporary expectation
                       vector[tmp_k] tmp_PX; // temporary precision right multiplied by temporary expectation
                       tmp_r = segment(r, tmp_m_pos, tmp_m);
                       tmp_s = segment(s, tmp_m_pos, tmp_m);
                       sum_r = rep_vector(0, S);
                       for(l in 1:tmp_m){ // skipping conventional "k" index to avoid confusion with tmp_k
                         sum_r[tmp_s[l]] += tmp_r[l];
                       }
                       tmp_k_ind = segment(ind_dim, tmp_k_pos, tmp_k);
                       tmp_V = rep_matrix(0, tmp_k, tmp_k);
                       for(l in 1:S){ // skipping conventional "k" index to avoid confusion with tmp_k
                         tmp_V += sum_r[l] * inter_cov[l, tmp_k_ind, tmp_k_ind];
                       }
                       tmp_P = inverse_spd(inverse_spd(P[j, tmp_k_ind, tmp_k_ind]) + tmp_V);
                       cur_P[tmp_k_ind, tmp_k_ind] += tmp_P;
                       tmp_X = X[tmp_k_ind, j];
                       tmp_PX = tmp_P * tmp_X;
                       sum_PX[tmp_k_ind] += tmp_PX;
                       sum_k -= tmp_k;
                       cur_R += log_determinant(tmp_P) - tmp_X' * tmp_PX;
                     }
                   }
                   // do necessary calculations for focal node
                   {
                     int tmp_code;
                     int tmp_k_pos;
                     int tmp_k;
                     tmp_code = code_nod[i];
                     tmp_k_pos = pos_dim[tmp_code];
                     tmp_k = n_dim[tmp_code];
                     if(tmp_k){
                       int tmp_k_ind[tmp_k];
                       matrix[tmp_k, tmp_k] tmp_P; // precision for focal node, subsetted to informative dimensions
                       vector[tmp_k] tmp_X; // expectation for focal node, subsetted to informative dimensions
                       tmp_k_ind = segment(ind_dim, tmp_k_pos, tmp_k);
                       tmp_P = cur_P[tmp_k_ind, tmp_k_ind];
                       tmp_X = mdivide_left_spd(tmp_P, sum_PX[tmp_k_ind]);
                       cur_X[tmp_k_ind] = tmp_X;
                       sum_k += tmp_k;
                       cur_R += tmp_X' * tmp_P * tmp_X - log_determinant(tmp_P);
                     }
                   }
                 }
                 // insert output quantities into full expectation, precision, and scalar arrays
                 X[, i] = cur_X;
                 P[i] = cur_P;
                 R[i] = 0.5 * (cur_R + sum_k * log(2 * pi()));
               }
               
               // MARGINALIZE OVER TREES
               // if not combining with prior on root, I belieive the integrated likelihood for each tree is simply given by scalars...
               // get scalar at each root node
               L = R[ind_root];
               // do log-sum-exp trick to sum together the exponentiated likelihoods
               max_L = max(L);
               L = exp(L - max_L);
               // Fitzjohn prior --> weight each tree by the probability it yielded the observed data, given the parameters
               // given by the sum of the likelihoods, weighted by the proportion they contribute to the "overall" likelihood
               // i.e., sum(L * L / sum(L)) = 1/sum(L) * sum(L^2) = log(sum(L^2)) - log(sum(L))
               // but need to multiply by max_L since L was divided by it to prevent underflow
               return log(dot_self(L)) - log(sum(L)) + max_L;
               // // after some testing, a more "naive" marginalization approach MAY yield better correlation estimates...
               // return log_sum_exp(R[ind_root]) - log(T);
             }
             
}

data {
  
  // basic quantities
  int N; // number of tips (MUST be same for each tree)
  int K; // number of traits
  int S; //numer of regimes
  int T; // number of trees
  int E; // number of edges across all trees
  int C; // number of "observation codes"
  
  // trait info
  int n_obs[N]; // vector of number of raw observations for each tip
  int ind_obs[sum(n_obs)]; // vector of raw observation indices for each tip
  matrix[K, sum(n_obs)] Y; // matrix of observed trait data (0 if unobserved)
  
  // topological info
  int n_inc[E]; // vector of number of increments for each edge
  int pos_inc[E]; // vector of starting position along increment vectors for each edge
  vector[sum(n_inc)] t; // vector of time durations of increments
  vector[sum(n_inc)] z; // vector of "explanatory trait" values which are transformed to rate scalars
  int s[sum(n_inc)]; // vector of states of increments
  int n_des[E]; // vector of number of descendants for each edge for each tree
  int ind_des[sum(n_des)]; //vector of descendant edge indices for each edge
  int pos_des[E]; // vector of starting positions along ind_des for each edge
  int ind_tip[T, N]; // matrix of tip edge indicators for each tree
  int ind_root[T]; // vector of root edge indicators for each tree
  int ind_prune[E - N * T]; // vector of edge indices indicating order to prune in
  
  // missing data info
  int n_dim[C]; // vector of number of informative dimensions for each observation code
  int ind_dim[sum(n_dim)]; // vector of informative dimension indices for each observation code
  int pos_dim[C]; // vector of starting positions along ind_dim for each observation code
  int code_obs[sum(n_obs)]; // observation codes for raw observations
  int code_tip[N]; // observation codes for tips
  int code_nod[E]; // observation codes for all nodes, indexed by ancestral edge
  
}

transformed data {
  
  vector[K] fixed_inter_var[S]; // "base" inter-tip variance vectors (see below) which fixes to the first variance to 1
  // this prevents unidentifiability with rate scalars
  
  for(i in 1:S){
    fixed_inter_var[i] = rep_vector(1, K);
  }
  
}

parameters {
  
  // intra-tip stuff (sometimes termed "phenotypic")
  corr_matrix[K] intra_cor; // intra-tip correlation matrix
  vector<lower = 0>[K] intra_var; // intra-tip variances
  
  // inter-tip stuff (sometimes termed "evolutionary")
  corr_matrix[K] inter_cor[S]; // inter-tip correlation matrices
  vector<lower = 0>[K - 1] unfixed_intervar_1; // unfixed inter-tip variance vectors (sometimes variance in this context is termed "rate") for 1st state
  vector<lower = 0>[K] unfixed_intervar_2[S - 1]; // unifixed rates for other states
  // rate for first trait in first state will always be set to 1 to prevent unidenfitiability
  
  // rate scalar stuff
  real<lower = 0> a; // intercept of rate scalar-explanatory trait relationship
  real b; // slope of rate scalar-explanatory trait relationship
  real<lower = 0> c; // offset of rate scalar-explanatory trait relationship
  
}

transformed parameters {
  
  vector[K] inter_var[S]; // compete inter-tip variance vectors
  vector[sum(n_inc)] r; // rate scalars
  
  inter_var[1, 1] = 1;
  if(K > 1){
    inter_var[1, 2:K] = unfixed_intervar_1;
  }
  if(S > 1){
    inter_var[2:S] = unfixed_intervar_2;
  }
  
  r = t .* (a * exp(b * z) + c);

}

model {
  
  // just call pruning function with current parameters and all input data
  target += prune(intra_cor, intra_var, inter_cor, inter_var,
                  N, K, S, T, E, C,
                  n_obs, ind_obs, Y,
                  n_inc, pos_inc, r, s, 
                  n_des, ind_des, pos_des, ind_tip, ind_root, ind_prune,
                  n_dim, ind_dim, pos_dim, code_obs, code_tip, code_nod);
                  
}

