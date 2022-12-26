#pragma once

#include "crn.hpp"

#include <string>
#include <ginac/ginac.h>

class MonotoneDependenciesCalculator{
public:
  /// constructor
  MonotoneDependenciesCalculator(const CRN& crn, const size_t mode):crn(crn), mode(mode){}

  /// calculates monoticity dependencies on the crn
  void run();

  /// logs all output
  /**
   * logs all output into the given dir; the output files will be characterized by the given model name
  */
  void log(const std::string dirname) const;

  /**
   * the output comprises of graph and data files with the naming convention as follows:
   * graph files: <modelname>_graph_<number_leading_zeros zeros><incrementing number of output, starting from zero>.dot
   * data files: <modelname>_output_<number_leading_zeros zeros><incrementing number of output, starting from zero>.txt
  */
  size_t number_leading_zeros = 5;
private:
  void calculate_next_index_vector(std::vector<int>& index_vectors, const std::vector<int>& max_index_vectors, size_t number_species);

  /* TODO: calculate consistency not for lines, but for whole matrices */
  void calculate_consistency(GiNaC::matrix& jacobian_transformed, std::vector<std::vector<int>>& transformed_sign_matrix,
    size_t number_species, bool& is_consistent, int i);

  void calculate_graph_consistency(size_t number_species, size_t number_graph_nodes, std::vector<std::vector<int>>& path_sign_matrix,
    const std::vector<std::vector<int>>& graph, bool& is_consistent);

  void calculate_path_signs(size_t number_species, size_t number_graph_nodes, std::vector<std::vector<int>>& path_sign_matrix,
    const std::vector<std::vector<int>>& graph, bool& is_consistent, std::vector<size_t>& path, int sign);

  /// saves output to log vector
  void save_consistent_transformation(const GiNaC::matrix time_derivatives, const GiNaC::matrix& jacobian_c, const GiNaC::matrix& jacobian_rr,
    const GiNaC::matrix& transformation_matrix,
    const GiNaC::matrix& jacobian_c_transformed, const GiNaC::matrix& jacobian_rr_transformed,
    const std::vector<std::vector<int>>& transformed_sign_matrix_c, const std::vector<std::vector<int>>& transformed_sign_matrix_rr,
    const std::vector<std::vector<int>>& graph, const std::vector<std::vector<int>>& path_sign_matrix);

  const CRN crn;
  const size_t mode;

  std::vector<std::string> log_vector;
  std::vector<std::string> graph_vector;
};

#define LEADINGZEROS 5
