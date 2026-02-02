#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]


#define TESTBROKE 1
// #define TESTCOMP 1

using namespace Rcpp;
using namespace R;
using namespace arma;


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]


arma::vec cpp_fmeasure_distribution(
    arma::mat predictions,
    arma::vec reference,
    int n_labels
) {
    Rcpp::Rcout << "cpp_fmeasure_distribution: START" << std::endl;
    // Bookkeeping
    arma::vec F1_dist = arma::vec(predictions.n_cols, arma::fill::zeros);
    // Constants
    double N = reference.n_elem;
    Rcpp::Rcout << "reference.n_elem: " << N << ", predictions.n_rows: " << predictions.n_rows << ", predictions.n_cols: " << predictions.n_cols << std::endl;
    arma::vec refsums_vec = arma::vec(n_labels, arma::fill::zeros);
    for (std::size_t i = 0; i < static_cast<std::size_t>(N); i++) {
        if (reference(i) < 0 || reference(i) >= n_labels) {
            Rcpp::Rcout << "reference out of bounds at i=" << i << ": " << reference(i) << ", nlabels=" << n_labels << std::endl;
        }
        refsums_vec(reference(i))++;
    }
    // Predeclared but transient variables
    arma::mat confusion_matrix;
    arma::vec TP_vec;
    double sumTP;
    arma::vec predsums_vec;
    arma::vec PPV_vec; // Positive Predictive Value AKA Precision
    arma::vec TPR_vec; // True Positive Rate AKA Recall AKA Sensitivity
    double TN;
    arma::vec TNR_vec; // True Negative Rate AKA Specificity
    // Class-wise metrics before taking the mean
    arma::rowvec F1;
    // cout << "here\n";
    for (std::size_t j = 0; j < static_cast<std::size_t>(predictions.n_cols); j++) {
        Rcpp::Rcout << "Processing column j=" << j << std::endl;
        // Re-initialize and fill confusion matrix
        confusion_matrix = arma::mat(n_labels, n_labels, arma::fill::zeros);
        for (std::size_t i = 0; i < static_cast<std::size_t>(reference.n_elem); i++) {
            int pred_idx = predictions(i, j);
            int ref_idx = reference(i);
            if (pred_idx < 0 || pred_idx >= n_labels || ref_idx < 0 || ref_idx >= n_labels) {
                Rcpp::Rcout << "Index out of bounds: pred_idx=" << pred_idx << ", ref_idx=" << ref_idx
                    << ", n_labels=" << n_labels << ", i=" << i << ", j=" << j << std::endl;
            }
            confusion_matrix(pred_idx, ref_idx)++;
        }
        Rcpp::Rcout << "Filled confusion_matrix for j=" << j << std::endl;

        // Compute transient temporary variables
        TP_vec = arma::diagvec(confusion_matrix);
        sumTP = arma::sum(TP_vec);
        predsums_vec = arma::vectorise(arma::sum(confusion_matrix, 1));
        PPV_vec = TP_vec / predsums_vec;
        TPR_vec = TP_vec / refsums_vec;

        // Compute F1
        F1 = arma::conv_to<arma::rowvec>::from( 2 * (PPV_vec % TPR_vec) / (PPV_vec + TPR_vec));
        F1.elem(arma::find_nonfinite(F1)).zeros();
        F1_dist(j) = arma::mean(F1);
        Rcpp::Rcout << "F1_dist[" << j << "] = " << F1_dist(j) << std::endl;
    }
    Rcpp::Rcout << "cpp_fmeasure_distribution: END" << std::endl;
    // cout << "\nPPV_dist\n";
    // PPV_dist.print(cout);
    return F1_dist;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List cpp_classification_metric_distributions(arma::mat predictions,
                                                  arma::vec reference,
                                                  int n_labels) {
  // Bookkeeping
  arma::vec F1_dist = arma::vec(predictions.n_cols, arma::fill::zeros);
  arma::vec jaccard_dist = arma::vec(predictions.n_cols, arma::fill::zeros);
  arma::vec PPV_dist = arma::vec(predictions.n_cols, arma::fill::zeros);
  arma::vec TPR_dist = arma::vec(predictions.n_cols, arma::fill::zeros);
  arma::vec TNR_dist = arma::vec(predictions.n_cols, arma::fill::zeros);
  arma::vec youden_dist = arma::vec(predictions.n_cols, arma::fill::zeros);
  // Constants
  double N = reference.n_elem;
  arma::vec refsums_vec = arma::vec(n_labels, arma::fill::zeros);
  for (std::size_t i = 0; i < static_cast<std::size_t>(N); i++) {
    refsums_vec(reference(i))++;
  }
  // Predeclared but transient variables
  arma::mat confusion_matrix;
  arma::vec TP_vec;
  double sumTP;
  arma::vec predsums_vec;
  arma::vec PPV_vec; // Positive Predictive Value AKA Precision
  arma::vec TPR_vec; // True Positive Rate AKA Recall AKA Sensitivity
  double TN;
  arma::vec TNR_vec; // True Negative Rate AKA Specificity
  // Class-wise metrics before taking the mean
  arma::rowvec F1;
  arma::vec jaccard;
  arma::vec youden;
  // cout << "here\n";
  for (std::size_t j = 0; j < static_cast<std::size_t>(predictions.n_cols); j++) {
    // Re-initialize and fill confusion matrix
    confusion_matrix = arma::mat(n_labels, n_labels, arma::fill::zeros);
    for (std::size_t i = 0; i < static_cast<std::size_t>(reference.n_elem); i++) {
      confusion_matrix(predictions(i, j),reference(i))++;
    }
    // confusion_matrix.print(cout);
    // cout << "\n\n";
    // Compute transient temporary variables
    TP_vec = arma::diagvec(confusion_matrix);
    sumTP = arma::sum(TP_vec);
    predsums_vec = arma::vectorise(arma::sum(confusion_matrix, 1));
    //predsums_vec = arma::vectorise(arma::sum(confusion_matrix, 0));
    PPV_vec = TP_vec / predsums_vec;
    TPR_vec = TP_vec / refsums_vec;
    // Compute F1
    F1 = arma::conv_to<arma::rowvec>::from( 2 * (PPV_vec % TPR_vec) / (PPV_vec + TPR_vec));
    F1.elem( arma::find_nonfinite(F1) ).zeros();
    F1_dist(j) = arma::mean(F1);
    // Compute Jaccard similarity and Youden's J index
    jaccard = arma::vec(n_labels, arma::fill::zeros);
    TNR_vec = arma::vec(n_labels, arma::fill::zeros);
    youden = arma::vec(n_labels, arma::fill::zeros);
    for (int i = 0; i < n_labels; i++) {
      /*cout << "predsum" << i << ": " << sum(confusion_matrix.col(i)) << "\n";
       cout << "refsums_vec" << i << ": " << sum(confusion_matrix.row(i)) << "\n";
       cout << "confusion_matrix(" << i << "," << i << "): " << confusion_matrix(i,i) << "\n\n";*/
      // Jaccard
      if (confusion_matrix(i,i) == 0) {
        jaccard(i) = 0;
      } else {
        //jaccard(i) = confusion_matrix(i,i)/ (sum(confusion_matrix.col(i)) + sum(confusion_matrix.row(i)) - confusion_matrix(i, i));
        jaccard(i) = confusion_matrix(i,i)/ (predsums_vec(i) + refsums_vec(i) - confusion_matrix(i, i));
      }
      // Youden
      TN = N-predsums_vec(i)-refsums_vec(i)+confusion_matrix(i,i);
      TNR_vec(i) = TN/(TN + arma::sum(predsums_vec(i)) - confusion_matrix(i, i));
      youden(i) = TPR_vec(i) + TNR_vec(i) - 1;
    }
    // Compute and store means for each posterior draw
    PPV_dist(j) = arma::mean(PPV_vec);
    TPR_dist(j) = arma::mean(TPR_vec);
    TNR_dist(j) = arma::mean(TNR_vec);
    jaccard_dist(j) = arma::mean(jaccard);
    youden_dist(j) = arma::mean(youden);
    // Trace prints
    // cout << "\n\nTP_vec\n";
    // TP_vec.print(cout);
    // cout << "\n\nTPR_vec\n";
    // TPR_vec.print(cout);
    // cout << "\n\nTNR_vec\n";
    // TNR_vec.print(cout);
    // cout << "\nrefsums_vec\n";
    // refsums_vec.print(cout);
    // cout << "\npredsums_vec\n";
    // predsums_vec.print(cout);
    // cout << "\nPPV_vec\n";
    // PPV_vec.print(cout);
    // cout << "\nPPV_mean\n" << PPV_dist(j) <<"\n";
    
  }
  // cout << "\nPPV_dist\n";
  // PPV_dist.print(cout);
  return Rcpp::List::create(Rcpp::Named("F1") = F1_dist,
                      Rcpp::Named("Jaccard") = jaccard_dist,
                      Rcpp::Named("Youden") = youden_dist,
                      Rcpp::Named("PPV") = PPV_dist,
                      Rcpp::Named("TPR") = TPR_dist,
                      Rcpp::Named("TNR") = TNR_dist);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec cpp_harmonic_rowmeans(arma::mat variables) {
  arma::vec result = arma::vec(variables.n_rows, arma::fill::zeros);
  int n = variables.n_cols;
  //variables.elem(find(variables<0)).zeros();
  for (std::size_t i = 0; i < static_cast<std::size_t>(variables.n_rows); i++) {
    result(i) = n/arma::sum(1/variables.row(i));
  }
  return result;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List cpp_compare_distributions(arma::vec dist1, arma::vec dist2) {
  double score1 = 0;
  double score2 = 0;
  std::size_t n1 = dist1.n_elem;
  std::size_t n2 = dist2.n_elem;
  double n = static_cast<double>(n1)*static_cast<double>(n2);
  for (std::size_t i = 0; i < n1; i++) {
    for (std::size_t j = 0; j < n2; j++) {
      if (dist1(i) < dist2(j)) {
        score2++;
      } else if (dist1(i) > dist2(j)) {
        score1++;
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("Dist1") = score1/n,
                      Rcpp::Named("Dist2") = score2/n);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec cpp_diff_distributions(arma::vec dist1, arma::vec dist2) {
  arma::vec res = arma::vec(dist1.n_elem*dist2.n_elem, arma::fill::zeros);
  std::size_t k = 0;
  for (std::size_t i = 0; i < dist1.n_elem; i++) {
    for (std::size_t j = 0; j < dist2.n_elem; j++) {
      res(k) = dist2(j) - dist1(i);
      k++;
    }
  }
  return res;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec cpp_rmse_dist(arma::mat predictions,
                        arma::mat reference,
                        int ndraws) {
  // Bookkeeping
  arma::vec rmse_dist = arma::vec(predictions.n_rows, arma::fill::zeros);
  int j = -1;
  for (std::size_t i = 0; i < static_cast<std::size_t>(predictions.n_rows); i++) {
    if (i % ndraws == 0) {
      j++;
    }
    // cout << j << "\n";
    rmse_dist(i) = std::sqrt(arma::mean(arma::square(predictions.row(i)-reference.row(j))));
  }
  return rmse_dist;
}

/*** R

# cpp_compare_distributions(c(0,1,2,3,4), c(5,4,3,2,1))
# plot(density(cpp_diff_distributions(c(0,1,2,3,4), c(5,4,3,2,1))))

testmat <- matrix(c(c(3, 4, 1), c(3,3,3), c(2, 3, 2), c(3,4,2)), byrow = TRUE, ncol = 3)
print(cpp_rmse_dist(testmat, matrix(c(c(3,4,2), c(2, 3, 3)), byrow = TRUE, ncol = 3), 2))

# df2 <- data.frame(
#   IDS = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
#   CESD = c(1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
# 
# set.seed(12345)
# Xbin <- sapply(seq(5), function (x) sample(c(0,1), 10, replace = TRUE, prob = c(0.1, 0.9)))
# Xord <- sapply(seq(5), function (x) sample(0:3, 10, replace = TRUE, prob = c(0.1, 0.3, 0.2, 0.4)))
# # #xbin <- sample(c(0,1), 100000, replace = TRUE)
# ybin <- sample(c(0,1), 10, replace = TRUE)
# yord <- sample(0:3, 10, replace = TRUE)

# #print(cpp_fscore(xbin, ybin, 2))
#print(cpp_mean_fscore_distribution(Xbin, ybin, 2))
#print(cpp_mean_fscore_distribution(Xord, yord, 4))


#print(cpp_jaccard_similarity_distribution(Xbin, ybin, 2))
#print(cpp_jaccard_similarity_distribution(Xord, yord, 4))

# fixedres <- cpp_classification_metric_distributions(matrix(IDS, ncol = 1), CESD, 2)
# binres <- cpp_classification_metric_distributions(Xbin, ybin, 2)
# ordres <- cpp_classification_metric_distributions(as.matrix(Xord), yord, 4)




# print(binres)

# print(fixedres$Jaccard) # should be 0.093...


# library(caret)

# IDS_fac <- as.factor(df2$IDS)
# CESD_fac <- factor(df2$CESD, levels = levels(IDS_fac))
# cm <- confusionMatrix(IDS_fac, CESD_fac, positive = NULL)
# print(cm$table)
# cm$byClass["Recall"]
# fixedres$TPR

# yord_fac <- as.factor(yord)
# for (i in 1:5) {
#   xord <- Xord[, i]
#   xord_fac = factor(xord, levels = levels(yord_fac))
#   cm = confusionMatrix(xord_fac, yord_fac)
#   # print(cm$table)
#   # Test F1
#   # F1 = cm$byClass[, "F1"]
#   # F1[which(!is.finite(F1))] = 0
#   # print(mean(F1))
#   # print(mean(cm$byClass[, "F1"]))
#   # Test Youden J
#   # print(mean(cm$byClass[, "Sensitivity"] + cm$byClass[, "Specificity"] - 1))
#   # Test Jaccard similarity
#   # jaccard = sapply(1:4, function(i) {
#   #   return(cm$table[i,i]/(sum(cm$table[i, ]) + sum(cm$table[, i]) - cm$table[i, i]))
#   # })
#   # print(mean(jaccard))
# }
# print(ordres$Jaccard)
*/
