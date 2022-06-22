functions {
  
  // BEAST of a function that calculates the likelihood of the trait data, marginalized across trees, given the current parameters
  real prune(real intra_var, vector alpha, vector beta, vector gamma, // parameters
             int N, int K, int T, int E, int C,
             int[] n_obs, int[] ind_obs, matrix Y,
             int[] n_inc, int[] pos_inc, vector t, matrix U, matrix W, matrix Z,
             int[] n_des, int[] ind_des, int[] pos_des, int[,] ind_tip, int[] ind_root, int[] ind_prune,
             int[] n_dim, int[] ind_dim, int[] pos_dim, int[] code_obs, int[] code_tip, int[] code_nod){
               
               vector[sum(n_inc)] r; // rate scalars
               vector[E] X; // expected trait values for each node given descendants (i.e., "expectations")
               vector[E] P; // precision matrix of expectations for each node  (i.e., "precisions")
               vector[E] R; // log scaling constants for partial likelihoods of expectations for each node (i.e., "scalars")
               
               // vector[T] L;
               // real max_L;
               
               int obs_counter; // counter to keep track of position along raw observation indices vector (ind_obs)
               
               // GET RATE SCALARS
               r = t .* (exp(U * alpha) .* inv_logit(W * beta) + exp(Z * gamma));
               
               // INITIALIZE TIPS
               // get expectation, precision, and scalar for each tip
               // could be made more efficient probably
               obs_counter = 1;
               for(i in 1:N){
                 // get indices and initialize output quantities
                 int tmp_n;
                 int ind_nod[T];
                 real cur_P; // current precision for ENTIRE tip
                 real cur_X; // current expectation for ENTIRE tip
                 real cur_R; // current scalar for ENTIRE tip
                 real sum_k; // involved in scalar calculations
                 tmp_n = n_obs[i];
                 ind_nod = ind_tip[, i];
                 cur_P = 0;
                 cur_X = 0;
                 cur_R = 0;
                 sum_k = 0;
                 if(tmp_n){
                    // this stuff could benefit from some precalculation
                   int tmp_n_ind[tmp_n];
                   vector[tmp_n] tmp_Y;
                   real sum_Y;
                   tmp_n_ind = segment(ind_obs, obs_counter, tmp_n);
                   tmp_Y = Y[1, tmp_n_ind]';
                   sum_Y = sum(tmp_Y);
                   cur_P = tmp_n / intra_var;
                   cur_X = sum_Y / tmp_n;
                   sum_k = 1 - tmp_n;
                   cur_R = sum_k * log(intra_var) - log(tmp_n) + (sum_Y * cur_X  - dot_self(tmp_Y)) / intra_var;
                 }
                 // insert output quantities into full expectation, precision, and scalar arrays
                 X[ind_nod] = rep_vector(cur_X, T);
                 P[ind_nod] = rep_vector(cur_P, T);
                 R[ind_nod] = rep_vector(0.5 * (cur_R + sum_k * log(2 * pi())), T);
                 // increment observation counter
                 obs_counter += tmp_n;
               }
               
               // PRUNING ALGORITHM
               for(i in ind_prune){
                 // get indices and initialize output quantities
                 int tmp_n;
                 int tmp_n_pos;
                 real cur_P; // current precision for focal node
                 real cur_X; // current expectation for focal node
                 real cur_R; // current scalar for focal node
                 real sum_PX; // sum of descendant expectations, weighted by their "propagated" precisions (see below)
                 real sum_k; // involved in scalar calculations
                 tmp_n = n_des[i];
                 tmp_n_pos = pos_des[i];
                 cur_P = 0;
                 cur_X = 0;
                 sum_PX = 0;
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
                     int tmp_k;
                     tmp_code = code_nod[j];
                     tmp_k = n_dim[tmp_code];
                     if(tmp_k){
                       //  below is not necessary, but will necessary again in some form for multivariate case!
                       // int tmp_k_pos;
                       // vector[tmp_m] tmp_r;
                       // int tmp_s[tmp_m];
                       // vector[S] sum_r; // summed rate scalars multiplied by time durations in each state for a given edge
                       // tmp_k_pos = pos_dim[tmp_code];
                       int tmp_m_pos;
                       int tmp_m;
                       real tmp_V; // temporary variance-covariance matrix of trait change along desendant edge
                       real tmp_P; // temporary post-propagation precision for expectation
                       real tmp_X; // temporary expectation
                       real tmp_PX; // temporary precision right multiplied by temporary expectation
                       tmp_m_pos = pos_inc[j];
                       tmp_m = n_inc[j];
                       tmp_V = sum(segment(r, tmp_m_pos, tmp_m));
                       tmp_P = 1 / (1 / P[j] + tmp_V);
                       cur_P += tmp_P;
                       tmp_X = X[j];
                       tmp_PX = tmp_P * tmp_X;
                       sum_PX += tmp_PX;
                       sum_k -= tmp_k;
                       cur_R += log(tmp_P) - tmp_X * tmp_PX;
                     }
                   }
                   // do necessary calculations for focal node
                   {
                     int tmp_code;
                     int tmp_k;
                     tmp_code = code_nod[i];
                     tmp_k = n_dim[tmp_code];
                     if(tmp_k){
                       // int tmp_k_pos;
                       // tmp_k_pos = pos_dim[tmp_code];
                       cur_X = sum_PX / cur_P;
                       sum_k += tmp_k;
                       cur_R += cur_X^2 * cur_P - log(cur_P);
                     }
                   }
                 }
                 // insert output quantities into full expectation, precision, and scalar arrays
                 X[i] = cur_X;
                 P[i] = cur_P;
                 R[i] = 0.5 * (cur_R + sum_k * log(2 * pi()));
               }
               
               // MARGINALIZE OVER TREES
               // if not combining with prior on root, I belieive the integrated likelihood for each tree is simply given by scalars...
               
               // // get scalar at each root node
               // L = R[ind_root];
               // // do log-sum-exp trick to sum together the exponentiated likelihoods
               // max_L = max(L);
               // L = exp(L - max_L);
               // // Fitzjohn prior --> weight each tree by the probability it yielded the observed data, given the parameters
               // // given by the sum of the likelihoods, weighted by the proportion they contribute to the "overall" likelihood
               // // i.e., sum(L * L / sum(L)) = 1/sum(L) * sum(L^2) = log(sum(L^2)) - log(sum(L))
               // // but need to multiply by max_L since L was divided by it to prevent underflow
               // return log(dot_self(L)) - log(sum(L)) + max_L;
               
               // after some testing, a more "naive" marginalization approach MAY yield better correlation estimates...
               return log_sum_exp(R[ind_root]) - log(T);
             }
             
}

data {
  
  // basic quantities
  int N; // number of tips (MUST be same for each tree)
  int K; // number of traits
  int T; // number of trees
  int E; // number of edges across all trees
  int C; // number of "observation codes"
  int A; // numnber of columns in design matrix for rate scalar maximum (alpha; threshold model only)
  int B; // number of columns in design matrix for rate scalar variation (beta)
  int G; // number of columns in design matrix for min rate scalar minimum (gamma)
  
  // trait info
  int n_obs[N]; // vector of number of raw observations for each tip
  int ind_obs[sum(n_obs)]; // vector of raw observation indices for each tip
  matrix[K, sum(n_obs)] Y; // matrix of observed trait data (0 if unobserved)
  
  // topological info
  int n_inc[E]; // vector of number of increments for each edge
  int pos_inc[E]; // vector of starting position along increment vectors for each edge
  vector[sum(n_inc)] t; // vector of time durations of increments
  matrix[sum(n_inc), A] U; // design matrix for rate scalar maximum
  matrix[sum(n_inc), B] W; // design matrix for rate scalar variation
  matrix[sum(n_inc), G] Z; // design matrix for rate scalar minimum
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

parameters {
  
  real<lower = 0> intra_var;
  
  // note that you will need multiple design matrices/betas for multiple traits
  vector[A] alpha; // rate scalar maximum coefficients ("scalars")
  vector[B] beta; // rate scalar variation coefficients ("slopes")
  vector[G] gamma; // rate scalar minimum coefficients ("offsets")
  
}

model {
  
  // just call pruning function with current parameters and all input data
  target += prune(intra_var, alpha, beta, gamma,
                  N, K, T, E, C,
                  n_obs, ind_obs, Y,
                  n_inc, pos_inc, t, U, W, Z,
                  n_des, ind_des, pos_des, ind_tip, ind_root, ind_prune,
                  n_dim, ind_dim, pos_dim, code_obs, code_tip, code_nod);
                  
}

