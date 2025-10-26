data{
  int K;   // the number of occupancy parameters
  int L;   // the number of detection parameters
  int N;   // the number of sites
  int NJ;  // the total number of observations
  int S;   // the number of unique combinations of 1s and 0s

  int obs[N];         // the number of replicate surveys at each site
  int start[N];       // index of the first replicate survey at each site
  matrix[S, K] x[N];  // design matrix

  // detection model covariates
  // vector[NJ] trl;
  vector[NJ] prot;

  // indicators of whether each species was detected at least once at a site
  int I1[N];
  int I2[N];
  int I3[N];
  int I4[N];

  // detection / non-detection data over all sites and replicate surveys
  vector[NJ] Y1;
  vector[NJ] Y2;
  vector[NJ] Y3;
  vector[NJ] Y4;
}

parameters{
  // detection model covariates
  vector[1] a1;
  vector[1] a2;
  vector[L] a3;
  vector[L] a4;

  real a1_c1;
  real a1_c2;
  real a2_c1;
  real a2_c2;
  real a2_c3;

  vector[K] beta;  // occupancy model covariates
}

model{

  vector[S] psi[N];   // probability of each combination of 1s and 0s
  vector[S] prob[N];  // psi * probability of detection history
  vector[S] z[N];     // log contribution of each site to the likelihood

  // probability of observing the detection history at each site
  vector[N] cd3;
  vector[N] cd4;

  vector[N] cd1_c11;
  vector[N] cd1_c10;
  vector[N] cd1_c01;
  vector[N] cd1_c00;
  vector[N] cd2_c111;
  vector[N] cd2_c110;
  vector[N] cd2_c100;
  vector[N] cd2_c010;
  vector[N] cd2_c011;
  vector[N] cd2_c001;
  vector[N] cd2_c101;
  vector[N] cd2_c000;

  beta ~ logistic(0, 1);
  a1 ~ logistic(0, 1);
  a2 ~ logistic(0, 1);
  a3 ~ logistic(0, 1);
  a4 ~ logistic(0, 1);

  a1_c1 ~ logistic(0, 1);
  a1_c2 ~ logistic(0, 1);
  a2_c1 ~ logistic(0, 1);
  a2_c2 ~ logistic(0, 1);
  a2_c3 ~ logistic(0, 1);
  
  for(i in 1:N){

    // Stan does not support ragged indexing.  This section accommodates 
    // different numbers of replicate surveys at each site

    // elements of the detection model design matrix
    vector[obs[i]] one;
    matrix[obs[i],1] w;
    matrix[obs[i],L] w2;
 
    // detection probability at each replicate survey
    vector[obs[i]] lp3;
    vector[obs[i]] lp4;

    vector[obs[i]] lp1_c11;
    vector[obs[i]] lp1_c10;
    vector[obs[i]] lp1_c01;
    vector[obs[i]] lp1_c00;
    vector[obs[i]] lp2_c111;
    vector[obs[i]] lp2_c110;
    vector[obs[i]] lp2_c100;
    vector[obs[i]] lp2_c010;
    vector[obs[i]] lp2_c011;
    vector[obs[i]] lp2_c001;
    vector[obs[i]] lp2_c101;
    vector[obs[i]] lp2_c000;

    // detection history of each species at each replicate survey
    vector[obs[i]] y1;
    vector[obs[i]] y2;
    vector[obs[i]] y3;
    vector[obs[i]] y4;

    y1 <- segment(Y1, start[i], obs[i]);
    y2 <- segment(Y2, start[i], obs[i]);
    y3 <- segment(Y3, start[i], obs[i]);
    y4 <- segment(Y4, start[i], obs[i]);

    one <- rep_vector(1,obs[i]);
    w <- rep_matrix(one,1);
    w2<-append_col(one,segment(prot, start[i], obs[i]));
  
    lp3 <- exp(w2 * a3) ./ (1 + exp(w2 * a3));
    lp4 <- exp(w2 * a4) ./ (1 + exp(w2 * a4));
    
    lp1_c11 <- exp(w * a1 + a1_c1 + a1_c2) ./ (1 + exp(w * a1 + a1_c1 + a1_c2));
    lp1_c10 <- exp(w * a1 + a1_c1) ./ (1 + exp(w * a1 + a1_c1));
    lp1_c01 <- exp(w * a1 + a1_c2) ./ (1 + exp(w * a1 + a1_c2));
    lp1_c00 <- exp(w * a1) ./ (1 + exp(w * a1));
    lp2_c111 <- exp(w * a2 + a2_c1 + a2_c2 + a2_c3) ./ (1 + exp(w * a2 + a2_c1 + a2_c2 + a2_c3));
    lp2_c110 <- exp(w * a2 + a2_c1 + a2_c2) ./ (1 + exp(w * a2 + a2_c1 + a2_c2));
    lp2_c100 <- exp(w * a2 + a2_c1) ./ (1 + exp(w * a2 + a2_c1));
    lp2_c010 <- exp(w * a2 + a2_c2) ./ (1 + exp(w * a2 + a2_c2));
    lp2_c001 <- exp(w * a2 + a2_c3) ./ (1 + exp(w * a2 + a2_c3));
    lp2_c101 <- exp(w * a2 + a2_c1 + a2_c3) ./ (1 + exp(w * a2 + a2_c1 + a2_c3));
    lp2_c011 <- exp(w * a2 + a2_c2 + a2_c3) ./ (1 + exp(w * a2 + a2_c2 + a2_c3));
    lp2_c000 <- exp(w * a2) ./ (1 + exp(w * a2));
    
    cd3[i] <- exp(sum(y3 .* log(lp3) + (1 - y3) .* log(1 - lp3)));
    cd4[i] <- exp(sum(y4 .* log(lp4) + (1 - y4) .* log(1 - lp4)));

    cd1_c11[i] <- exp(sum(y1 .* log(lp1_c11) + (1 - y1) .* log(1 - lp1_c11)));
    cd1_c10[i] <- exp(sum(y1 .* log(lp1_c10) + (1 - y1) .* log(1 - lp1_c10)));
    cd1_c01[i] <- exp(sum(y1 .* log(lp1_c01) + (1 - y1) .* log(1 - lp1_c01)));
    cd1_c00[i] <- exp(sum(y1 .* log(lp1_c00) + (1 - y1) .* log(1 - lp1_c00)));
    cd2_c111[i] <- exp(sum(y1 .* log(lp2_c111) + (1 - y1) .* log(1 - lp2_c111)));
    cd2_c110[i] <- exp(sum(y1 .* log(lp2_c110) + (1 - y1) .* log(1 - lp2_c110)));
    cd2_c100[i] <- exp(sum(y1 .* log(lp2_c100) + (1 - y1) .* log(1 - lp2_c100)));
    cd2_c010[i] <- exp(sum(y1 .* log(lp2_c010) + (1 - y1) .* log(1 - lp2_c010)));
    cd2_c001[i] <- exp(sum(y1 .* log(lp2_c001) + (1 - y1) .* log(1 - lp2_c001)));
    cd2_c011[i] <- exp(sum(y1 .* log(lp2_c011) + (1 - y1) .* log(1 - lp2_c011)));
    cd2_c101[i] <- exp(sum(y1 .* log(lp2_c101) + (1 - y1) .* log(1 - lp2_c101)));
    cd2_c000[i] <- exp(sum(y1 .* log(lp2_c000) + (1 - y1) .* log(1 - lp2_c000)));

    psi[i] <- softmax(x[i] * beta);

    prob[i][1] <- psi[i][1] * cd1_c11[i] * cd2_c111[i] * cd3[i] * cd4[i];
    prob[i][2] <- psi[i][2] * cd1_c10[i] * cd2_c110[i] * cd3[i];
    prob[i][3] <- psi[i][3] * cd1_c01[i] * cd2_c101[i] * cd4[i];
    prob[i][4] <- psi[i][4] * cd1_c00[i] * cd2_c100[i];
    prob[i][5] <- psi[i][5] * cd1_c11[i] * cd3[i] * cd4[i];
    prob[i][6] <- psi[i][6] * cd1_c10[i] * cd3[i];
    prob[i][7] <- psi[i][7] * cd1_c01[i] * cd4[i];
    prob[i][8] <- psi[i][8] * cd1_c00[i];
    prob[i][9] <- psi[i][9] * cd2_c011[i] * cd3[i] * cd4[i];
    prob[i][10] <- psi[i][10] * cd2_c010[i] * cd3[i];
    prob[i][11] <- psi[i][11] * cd2_c001[i] * cd4[i];
    prob[i][12] <- psi[i][12] * cd2_c000[i];
    prob[i][13] <- psi[i][13] * cd3[i] * cd4[i];
    prob[i][14] <- psi[i][14] * cd3[i];
    prob[i][15] <- psi[i][15] * cd4[i];
    prob[i][16] <- psi[i][16];

    z[i][1] <- I1[i] * I2[i] * I3[i] * I4[i] * log(prob[i][1]);
    z[i][2] <- I1[i] * I2[i] * I3[i] * (1 - I4[i]) *
                log(prob[i][1] + prob[i][2]);
    z[i][3] <- I1[i] * I2[i] * (1 - I3[i]) * I4[i] *
                log(prob[i][1] + prob[i][3]);
    z[i][4] <- I1[i] * I2[i] * (1 - I3[i]) * (1 - I4[i]) *
                log(prob[i][1] + prob[i][2] + prob[i][3] + prob[i][4]);
    z[i][5] <- I1[i] * (1 - I2[i]) * I3[i] * I4[i] *
                log(prob[i][1] + prob[i][5]);
    z[i][6] <- I1[i] * (1 - I2[i]) * I3[i] * (1 - I4[i]) *
                log(prob[i][1] + prob[i][2] + prob[i][5] + prob[i][6]);
    z[i][7] <- I1[i] * (1 - I2[i]) * (1 - I3[i]) * I4[i] *
                log(prob[i][1] + prob[i][3] + prob[i][5] + prob[i][7]);
    z[i][8] <- I1[i] * (1 - I2[i]) * (1 - I3[i]) * (1 - I4[i]) *
                log(prob[i][1] + prob[i][2] + prob[i][3] + prob[i][4] +
                    prob[i][5] + prob[i][6] + prob[i][7] + prob[i][8]);
    z[i][9] <- (1 - I1[i]) * I2[i] * I3[i] * I4[i] *
                log(prob[i][1] + prob[i][9]);
    z[i][10] <- (1 - I1[i]) * I2[i] * I3[i] * (1 - I4[i]) *
                  log(prob[i][1] + prob[i][2] + prob[i][9] + prob[i][10]);
    z[i][11] <- (1 - I1[i]) * I2[i] * (1 - I3[i]) * I4[i] *
                  log(prob[i][1] + prob[i][3] + prob[i][9] + prob[i][11]);
    z[i][12] <- (1 - I1[i]) * I2[i] * (1 - I3[i]) * (1 - I4[i]) *
                  log(prob[i][1] + prob[i][2] + prob[i][3] + prob[i][4] +
                      prob[i][9] + prob[i][10] + prob[i][11] + prob[i][12]);
    z[i][13] <- (1 - I1[i]) * (1 - I2[i]) * I3[i] * I4[i] *
                  log(prob[i][1] + prob[i][5] + prob[i][9] + prob[i][13]);
    z[i][14] <- (1 - I1[i]) * (1 - I2[i]) * I3[i] * (1 - I4[i]) *
                  log(prob[i][1] + prob[i][2] + prob[i][5] + prob[i][6] +
                      prob[i][9] + prob[i][10] + prob[i][13] + prob[i][14]);
    z[i][15] <- (1 - I1[i]) * (1 - I2[i]) * (1 - I3[i]) * I4[i] *
                  log(prob[i][1] + prob[i][3] + prob[i][5] + prob[i][7] +
                      prob[i][9] + prob[i][11] + prob[i][13] + prob[i][15]);
    z[i][16] <- (1 - I1[i]) * (1 - I2[i]) * (1 - I3[i]) * (1 - I4[i]) *
                  log(sum(prob[i]));

    // likelihood not included in Stan; this allows specification of
    // non-standard likelihoods
    increment_log_prob(sum(z[i]));
  }
}
generated quantities{
  vector[S] psi[N];   // probability of each combination of 1s and 0s
  vector[S] prob[N];  // psi * probability of detection history
  vector[S] z[N];     // log contribution of each site to the likelihood

  // probability of observing the detection history at each site
  vector[N] cd3;
  vector[N] cd4;

  vector[N] cd1_c11;
  vector[N] cd1_c10;
  vector[N] cd1_c01;
  vector[N] cd1_c00;
  vector[N] cd2_c111;
  vector[N] cd2_c110;
  vector[N] cd2_c100;
  vector[N] cd2_c010;
  vector[N] cd2_c011;
  vector[N] cd2_c001;
  vector[N] cd2_c101;
  vector[N] cd2_c000;
  
  vector[N] ll;  // log likelihood
  
  for(i in 1:N){

    // Stan does not support ragged indexing.  This section accommodates 
    // different numbers of replicate surveys at each site

    // elements of the detection model design matrix
    vector[obs[i]] one;
    matrix[obs[i],1] w;
    matrix[obs[i],L] w2;
 
    // detection probability at each replicate survey
    vector[obs[i]] lp3;
    vector[obs[i]] lp4;

    vector[obs[i]] lp1_c11;
    vector[obs[i]] lp1_c10;
    vector[obs[i]] lp1_c01;
    vector[obs[i]] lp1_c00;
    vector[obs[i]] lp2_c111;
    vector[obs[i]] lp2_c110;
    vector[obs[i]] lp2_c100;
    vector[obs[i]] lp2_c010;
    vector[obs[i]] lp2_c011;
    vector[obs[i]] lp2_c001;
    vector[obs[i]] lp2_c101;
    vector[obs[i]] lp2_c000;

    // detection history of each species at each replicate survey
    vector[obs[i]] y1;
    vector[obs[i]] y2;
    vector[obs[i]] y3;
    vector[obs[i]] y4;

    y1 <- segment(Y1, start[i], obs[i]);
    y2 <- segment(Y2, start[i], obs[i]);
    y3 <- segment(Y3, start[i], obs[i]);
    y4 <- segment(Y4, start[i], obs[i]);

    one <- rep_vector(1,obs[i]);
    w <- rep_matrix(one,1);
    w2<-append_col(one,segment(prot, start[i], obs[i]));
  
    lp3 <- exp(w2 * a3) ./ (1 + exp(w2 * a3));
    lp4 <- exp(w2 * a4) ./ (1 + exp(w2 * a4));
    
    lp1_c11 <- exp(w * a1 + a1_c1 + a1_c2) ./ (1 + exp(w * a1 + a1_c1 + a1_c2));
    lp1_c10 <- exp(w * a1 + a1_c1) ./ (1 + exp(w * a1 + a1_c1));
    lp1_c01 <- exp(w * a1 + a1_c2) ./ (1 + exp(w * a1 + a1_c2));
    lp1_c00 <- exp(w * a1) ./ (1 + exp(w * a1));
    lp2_c111 <- exp(w * a2 + a2_c1 + a2_c2 + a2_c3) ./ (1 + exp(w * a2 + a2_c1 + a2_c2 + a2_c3));
    lp2_c110 <- exp(w * a2 + a2_c1 + a2_c2) ./ (1 + exp(w * a2 + a2_c1 + a2_c2));
    lp2_c100 <- exp(w * a2 + a2_c1) ./ (1 + exp(w * a2 + a2_c1));
    lp2_c010 <- exp(w * a2 + a2_c2) ./ (1 + exp(w * a2 + a2_c2));
    lp2_c001 <- exp(w * a2 + a2_c3) ./ (1 + exp(w * a2 + a2_c3));
    lp2_c101 <- exp(w * a2 + a2_c1 + a2_c3) ./ (1 + exp(w * a2 + a2_c1 + a2_c3));
    lp2_c011 <- exp(w * a2 + a2_c2 + a2_c3) ./ (1 + exp(w * a2 + a2_c2 + a2_c3));
    lp2_c000 <- exp(w * a2) ./ (1 + exp(w * a2));
    
    cd3[i] <- exp(sum(y3 .* log(lp3) + (1 - y3) .* log(1 - lp3)));
    cd4[i] <- exp(sum(y4 .* log(lp4) + (1 - y4) .* log(1 - lp4)));

    cd1_c11[i] <- exp(sum(y1 .* log(lp1_c11) + (1 - y1) .* log(1 - lp1_c11)));
    cd1_c10[i] <- exp(sum(y1 .* log(lp1_c10) + (1 - y1) .* log(1 - lp1_c10)));
    cd1_c01[i] <- exp(sum(y1 .* log(lp1_c01) + (1 - y1) .* log(1 - lp1_c01)));
    cd1_c00[i] <- exp(sum(y1 .* log(lp1_c00) + (1 - y1) .* log(1 - lp1_c00)));
    cd2_c111[i] <- exp(sum(y1 .* log(lp2_c111) + (1 - y1) .* log(1 - lp2_c111)));
    cd2_c110[i] <- exp(sum(y1 .* log(lp2_c110) + (1 - y1) .* log(1 - lp2_c110)));
    cd2_c100[i] <- exp(sum(y1 .* log(lp2_c100) + (1 - y1) .* log(1 - lp2_c100)));
    cd2_c010[i] <- exp(sum(y1 .* log(lp2_c010) + (1 - y1) .* log(1 - lp2_c010)));
    cd2_c001[i] <- exp(sum(y1 .* log(lp2_c001) + (1 - y1) .* log(1 - lp2_c001)));
    cd2_c011[i] <- exp(sum(y1 .* log(lp2_c011) + (1 - y1) .* log(1 - lp2_c011)));
    cd2_c101[i] <- exp(sum(y1 .* log(lp2_c101) + (1 - y1) .* log(1 - lp2_c101)));
    cd2_c000[i] <- exp(sum(y1 .* log(lp2_c000) + (1 - y1) .* log(1 - lp2_c000)));

    psi[i] <- softmax(x[i] * beta);

    prob[i][1] <- psi[i][1] * cd1_c11[i] * cd2_c111[i] * cd3[i] * cd4[i];
    prob[i][2] <- psi[i][2] * cd1_c10[i] * cd2_c110[i] * cd3[i];
    prob[i][3] <- psi[i][3] * cd1_c01[i] * cd2_c101[i] * cd4[i];
    prob[i][4] <- psi[i][4] * cd1_c00[i] * cd2_c100[i];
    prob[i][5] <- psi[i][5] * cd1_c11[i] * cd3[i] * cd4[i];
    prob[i][6] <- psi[i][6] * cd1_c10[i] * cd3[i];
    prob[i][7] <- psi[i][7] * cd1_c01[i] * cd4[i];
    prob[i][8] <- psi[i][8] * cd1_c00[i];
    prob[i][9] <- psi[i][9] * cd2_c011[i] * cd3[i] * cd4[i];
    prob[i][10] <- psi[i][10] * cd2_c010[i] * cd3[i];
    prob[i][11] <- psi[i][11] * cd2_c001[i] * cd4[i];
    prob[i][12] <- psi[i][12] * cd2_c000[i];
    prob[i][13] <- psi[i][13] * cd3[i] * cd4[i];
    prob[i][14] <- psi[i][14] * cd3[i];
    prob[i][15] <- psi[i][15] * cd4[i];
    prob[i][16] <- psi[i][16];

    z[i][1] <- I1[i] * I2[i] * I3[i] * I4[i] * log(prob[i][1]);
    z[i][2] <- I1[i] * I2[i] * I3[i] * (1 - I4[i]) *
                log(prob[i][1] + prob[i][2]);
    z[i][3] <- I1[i] * I2[i] * (1 - I3[i]) * I4[i] *
                log(prob[i][1] + prob[i][3]);
    z[i][4] <- I1[i] * I2[i] * (1 - I3[i]) * (1 - I4[i]) *
                log(prob[i][1] + prob[i][2] + prob[i][3] + prob[i][4]);
    z[i][5] <- I1[i] * (1 - I2[i]) * I3[i] * I4[i] *
                log(prob[i][1] + prob[i][5]);
    z[i][6] <- I1[i] * (1 - I2[i]) * I3[i] * (1 - I4[i]) *
                log(prob[i][1] + prob[i][2] + prob[i][5] + prob[i][6]);
    z[i][7] <- I1[i] * (1 - I2[i]) * (1 - I3[i]) * I4[i] *
                log(prob[i][1] + prob[i][3] + prob[i][5] + prob[i][7]);
    z[i][8] <- I1[i] * (1 - I2[i]) * (1 - I3[i]) * (1 - I4[i]) *
                log(prob[i][1] + prob[i][2] + prob[i][3] + prob[i][4] +
                    prob[i][5] + prob[i][6] + prob[i][7] + prob[i][8]);
    z[i][9] <- (1 - I1[i]) * I2[i] * I3[i] * I4[i] *
                log(prob[i][1] + prob[i][9]);
    z[i][10] <- (1 - I1[i]) * I2[i] * I3[i] * (1 - I4[i]) *
                  log(prob[i][1] + prob[i][2] + prob[i][9] + prob[i][10]);
    z[i][11] <- (1 - I1[i]) * I2[i] * (1 - I3[i]) * I4[i] *
                  log(prob[i][1] + prob[i][3] + prob[i][9] + prob[i][11]);
    z[i][12] <- (1 - I1[i]) * I2[i] * (1 - I3[i]) * (1 - I4[i]) *
                  log(prob[i][1] + prob[i][2] + prob[i][3] + prob[i][4] +
                      prob[i][9] + prob[i][10] + prob[i][11] + prob[i][12]);
    z[i][13] <- (1 - I1[i]) * (1 - I2[i]) * I3[i] * I4[i] *
                  log(prob[i][1] + prob[i][5] + prob[i][9] + prob[i][13]);
    z[i][14] <- (1 - I1[i]) * (1 - I2[i]) * I3[i] * (1 - I4[i]) *
                  log(prob[i][1] + prob[i][2] + prob[i][5] + prob[i][6] +
                      prob[i][9] + prob[i][10] + prob[i][13] + prob[i][14]);
    z[i][15] <- (1 - I1[i]) * (1 - I2[i]) * (1 - I3[i]) * I4[i] *
                  log(prob[i][1] + prob[i][3] + prob[i][5] + prob[i][7] +
                      prob[i][9] + prob[i][11] + prob[i][13] + prob[i][15]);
    z[i][16] <- (1 - I1[i]) * (1 - I2[i]) * (1 - I3[i]) * (1 - I4[i]) *
                  log(sum(prob[i]));
                  
    ll[i] <- sum(z[i]);
  }
}