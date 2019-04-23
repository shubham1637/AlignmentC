#include <iostream>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <cmath> //sqrt

#include <assert.h>
#define PRECONDITION(condition, message) \
  assert (condition);

struct SimMatrix
{
  std::vector<double> data;
  int n_row; 
  int n_col;
};

/** @brief Will add the outer product of vectors d1 and d2 to the similarity matrix s
 *
 * @param d1 Input vector 1
 * @param d2 Input vector 2
 * @param s Output similarity matrix (outer product will get added to existing matrix)
 *
*/
void addOuterProduct(const std::vector<double>& d1, const std::vector<double>& d2, SimMatrix& s)
{
  PRECONDITION(s.n_row == s.n_col, "Similarity matrix needs to be symmetric");
  PRECONDITION(s.n_row == d1.size(), "Data vector size (vector 1) needs to equal matrix dimension");
  PRECONDITION(s.n_row == d2.size(), "Data vector size (vector 2) needs to equal matrix dimension");

  int n = s.n_row;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      s.data[i*n + j] += d1[i] * d2[j];
    }
  }
} 

void addOuterProductEuclidian(const std::vector<double>& d1, const std::vector<double>& d2, SimMatrix& s)
{
  PRECONDITION(s.n_row == s.n_col, "Similarity matrix needs to be symmetric");
  PRECONDITION(s.n_row == d1.size(), "Data vector size (vector 1) needs to equal matrix dimension");
  PRECONDITION(s.n_row == d2.size(), "Data vector size (vector 2) needs to equal matrix dimension");

  int n = s.n_row;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      s.data[i*n + j] += (d1[i] - d2[j]) * (d1[i] - d2[j]);
    }
  }
} 

void meanNorm(std::vector<double>& v)
{
  double mean = std::accumulate( v.begin(), v.end(), 0.0)/v.size(); 
  for (auto& vv : v) vv /= mean;
}

double meanVecOfVec(const std::vector<std::vector<double>>& vec)
{
  double mean = 0.0; 
  for (const auto& v : vec) mean += std::accumulate( v.begin(), v.end(), 0.0)/v.size(); 
  return mean / vec.size(); 
}


// functor for getting sum of previous result and square of current element
template<typename T>
struct square
{
    T operator()(const T& Left, const T& Right) const
    {   
        return (Left + Right*Right);
    }
};


/** @brief Compute the L2 norm of all numbers in a set of vectors of vectors
 *
 * @note Treats the input as one long vector of which it takes the L2 norm 
 *
*/
double L2NormVecOfVec(const std::vector<std::vector<double>>& vec)
{
  double mean = 0.0; 
  for (const auto& v : vec) mean += std::accumulate( v.begin(), v.end(), 0.0, square<double>());
  return std::sqrt(mean);
}

/** @brief Will compute the outer product between all fragment ions and add each one to the similarity matrix
 *
 * @param d1 List of input chromatograms from run 1
 * @param d2 List of input chromatograms from run 2
 * @param s Output similarity matrix (outer product will get added to existing matrix)
 *
*/
void OuterProdMeanNormAllFunc(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2, SimMatrix& s)
{
  PRECONDITION(!d1.empty(), "Vectors cannot be empty");
  PRECONDITION(d1.size() == d2.size(), "Number of fragments needs to be equal");
  PRECONDITION(s.n_row == s.n_col, "Similarity matrix needs to be symmetric");
  PRECONDITION(s.n_row == d1[0].size(), "Data vector size (vector 1) needs to equal matrix dimension");
  PRECONDITION(s.n_row == d2[0].size(), "Data vector size (vector 2) needs to equal matrix dimension");

  // calculate overall mean for normalization
  // TODO: why divide again by number of samples?
  double mean_d1 = meanVecOfVec(d1);
  double mean_d2 = meanVecOfVec(d2);

  int n_frag = d1.size();
  for (int i = 0; i < n_frag; i++)
  {
    // Copy vector and divide by global mean (alternatively normalize each vector individually)
    std::vector<double> tmp1(d1[i]), tmp2(d2[i]);
    for (auto& vv : tmp1) vv /= mean_d1;
    for (auto& vv : tmp2) vv /= mean_d2;
    // meanNorm(tmp1);
    // meanNorm(tmp2);
    addOuterProduct(tmp1, tmp2, s);
  }


// TODO: is MeanNormA a scalar or a vector?

// OuterProdMeanNormAll6Func <- function(data, pep, runA, runB){
//     num_of_frag <- length(data[[runA]][[pep]])
//     num_of_samplesA <- length(data[[runA]][[pep]][[1]][,1])
//     num_of_samplesB <- length(data[[runB]][[pep]][[1]][,1])
//     MeanNormA <- sapply(data[[runA]][[pep]], function(x) sum(x[,2])/num_of_samplesA)
//     MeanNormA <- mean(MeanNormA)
//     MeanNormB <- sapply(data[[runB]][[pep]], function(x) sum(x[,2])/num_of_samplesB)
//     MeanNormB <- mean(MeanNormB)
//     outerProdList <- list()
//     for (i in 1:num_of_frag){
//         NormIntensityA <- data[[runA]][[pep]][[i]][,2]/MeanNormA
//         NormIntensityB <- data[[runB]][[pep]][[i]][,2]/MeanNormB
//         outerProdList[[i]] <- outer(NormIntensityA, NormIntensityB)
//     }
// return(outerProdList) }
}


/** @brief Will compute the outer product between all fragment ions and add each one to the similarity matrix
 *
 * @param d1 List of input chromatograms from run 1
 * @param d2 List of input chromatograms from run 2
 * @param s Output similarity matrix (outer product will get added to existing matrix)
 *
*/
void OuterProdL2NormAllFunc(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2, SimMatrix& s)
{
  PRECONDITION(!d1.empty(), "Vectors cannot be empty");
  PRECONDITION(d1.size() == d2.size(), "Number of fragments needs to be equal");
  PRECONDITION(s.n_row == s.n_col, "Similarity matrix needs to be symmetric");
  PRECONDITION(s.n_row == d1[0].size(), "Data vector size (vector 1) needs to equal matrix dimension");
  PRECONDITION(s.n_row == d2[0].size(), "Data vector size (vector 2) needs to equal matrix dimension");

  // calculate overall mean for normalization
  // TODO: why divide again by number of samples?
  double mean_d1 = L2NormVecOfVec(d1);
  double mean_d2 = L2NormVecOfVec(d2);

  int n_frag = d1.size();
  for (int i = 0; i < n_frag; i++)
  {
    // Copy vector and divide by global L2 norm
    std::vector<double> tmp1(d1[i]), tmp2(d2[i]);
    for (auto& vv : tmp1) vv /= mean_d1;
    for (auto& vv : tmp2) vv /= mean_d2;
    addOuterProduct(tmp1, tmp2, s);
  }


// TODO: is L2NormA a scalar or a vector?

// OuterProdL2NormAllFunc <- function(data, pep, runA, runB){
//     num_of_frag <- length(data[[runA]][[pep]])
//     L2NormA <- sapply(data[[runA]][[pep]], function(x) x[,2])
//     L2NormA <- sqrt(rowSums(L2NormA^2))
//     L2NormB <- sapply(data[[runB]][[pep]], function(x) x[,2])
//     L2NormB <- sqrt(rowSums(L2NormB^2))
//     outerProdList <- list()
//     for (i in 1:num_of_frag){
//         NormIntensityA <- data[[runA]][[pep]][[i]][,2]/L2NormA
//         NormIntensityA[is.nan(NormIntensityA)] <-0
//         NormIntensityB <- data[[runB]][[pep]][[i]][,2]/L2NormB
//         NormIntensityB[is.nan(NormIntensityB)] <-0
//         outerProdList[[i]] <- outer(NormIntensityA, NormIntensityB)
//     }
//     return(outerProdList) }
}


void OuterProdEuclFunc(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2, SimMatrix& s)
{
  PRECONDITION(!d1.empty(), "Vectors cannot be empty");
  PRECONDITION(d1.size() == d2.size(), "Number of fragments needs to be equal");
  PRECONDITION(s.n_row == s.n_col, "Similarity matrix needs to be symmetric");
  PRECONDITION(s.n_row == d1[0].size(), "Data vector size (vector 1) needs to equal matrix dimension");
  PRECONDITION(s.n_row == d2[0].size(), "Data vector size (vector 2) needs to equal matrix dimension");

  // calculate overall mean for normalization
  // TODO: why divide again by number of samples?
  double mean_d1 = meanVecOfVec(d1);
  double mean_d2 = meanVecOfVec(d2);

  int n_frag = d1.size();
  for (int i = 0; i < n_frag; i++)
  {
    // Copy vector and divide by global L2 norm
    std::vector<double> tmp1(d1[i]), tmp2(d2[i]);
    for (auto& vv : tmp1) vv /= mean_d1;
    for (auto& vv : tmp2) vv /= mean_d2;
    addOuterProductEuclidian(tmp1, tmp2, s);
  }


// OuterProdEuclFunc <- function(data, pep, runA, runB){
//     num_of_frag <- length(data[[runA]][[pep]])
//     num_of_samplesA <- length(data[[runA]][[pep]][[1]][,1])
//     num_of_samplesB <- length(data[[runB]][[pep]][[1]][,1])
//     MeanNormA <- sapply(data[[runA]][[pep]], function(x) sum(x[,2])/num_of_samplesA)
//     MeanNormA <- mean(MeanNormA)
//     MeanNormB <- sapply(data[[runB]][[pep]], function(x) sum(x[,2])/num_of_samplesB)
//     MeanNormB <- mean(MeanNormB)
//     outerProdList <- list()
//     for (i in 1:num_of_frag){
//         NormIntensityA <- data[[runA]][[pep]][[i]][,2]/MeanNormA
//         NormIntensityB <- data[[runB]][[pep]][[i]][,2]/MeanNormB
//         outerProdList[[i]] <- (outer(NormIntensityA, NormIntensityB, FUN = "-"))**2
//     }
// return(outerProdList)}
}

void getSimilarityMatrix(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2, SimMatrix& s, const std::string type)
{
  PRECONDITION(!d1.empty(), "Vectors cannot be empty");
  PRECONDITION(d1.size() == d2.size(), "Number of fragments needs to be equal");
  PRECONDITION(s.n_row == s.n_col, "Similarity matrix needs to be symmetric");
  PRECONDITION(s.data.size() == s.n_col * s.n_row , "Similarity matrix length needs to be consistent");
  PRECONDITION(s.n_row == d1[0].size(), "Data vector size (vector 1) needs to equal matrix dimension");
  PRECONDITION(s.n_row == d2[0].size(), "Data vector size (vector 2) needs to equal matrix dimension");

  if (type == "dotProduct") OuterProdMeanNormAllFunc(d1, d2, s);
  if (type == "cosineAngle") OuterProdL2NormAllFunc(d1, d2, s);
  if (type == "euclidianDist") OuterProdEuclFunc(d1, d2, s);
  if (type == "dotProductMasked") 
  {
    SimMatrix stmp1;
    stmp1.data.resize(s.data.size(), 0);
    SimMatrix stmp2;
    stmp2.data.resize(s.data.size(), 0);


    OuterProdMeanNormAllFunc(d1, d2, stmp1);
    OuterProdL2NormAllFunc(d1, d2, stmp2);

    // TODO: create MASK

    // OuterProdEuclFunc(d1, d2, s);
    //      s1 <- add(OuterProdNormAll6)
    //      OuterProdL2NormAll <- OuterProdL2NormAllFunc(data, pep, runA, runB)
    //      s2 <- cos(2*acos(pmin(add(OuterProdL2NormAll), 1)))
    //      MASK <- (s1 > quantile(s1, dotProdThresh))
    //      AngleGreat <- (((1*MASK)*s2) + (1-MASK)) > cosAngleThresh
  }


}

