#include "monotonedependenciescalculator.hpp"
#include "math.hpp"

#include <fstream>
#include <iostream>
#include <sstream>

void MonotoneDependenciesCalculator::run(){
  /* generate jacobians, possible vectors for the transformation matrix, ... */
  GiNaC::matrix tmp = Math::jacobian(crn.reaction_list, crn.species_list);
  GiNaC::matrix jacobian_c = crn.stoichiometric.mul(tmp);
  tmp = Math::jacobian(crn.reaction_list, crn.reaction_rate_list);
  GiNaC::matrix jacobian_rr = crn.stoichiometric.mul(tmp);
  GiNaC::matrix transformation_matrix(crn.number_species, crn.number_species);
  std::vector<std::vector<int>> transformation_matrix_ints(crn.number_species, std::vector<int>(crn.number_species, 0));
  std::vector<std::vector<double>> transformation_matrix_doubles(crn.number_species, std::vector<double>(crn.number_species, 0));
  std::vector<std::vector<int>> possible_vectors = Math::binary_matrix(crn.number_species);
  std::vector<int> index_vectors;
  std::vector<int> max_index_vectors;

  /* calculate time derivatives for output */
  GiNaC::matrix tmp2 (crn.number_reactions, 1);
  for(size_t i = 0; i < crn.number_reactions; i++){
    tmp2(i,0) = crn.reaction_list[i];
  }
  GiNaC::matrix time_derivatives = crn.stoichiometric.mul(tmp2);
  size_t number_graph_nodes;
    number_graph_nodes = crn.number_species + crn.number_reaction_rates;
  /*
   * graph matrix; the value of graph[i][j] may be one of {-1,0,1} and describes the sign that an edge between the two nodes is marked with;
   * '0' means no edge.
   * the graph matrix doesn't have to be quadratic as there cannot be paths connecting reaction rates to each other
  */
  std::vector<std::vector<int>> graph (crn.number_species, std::vector<int>(number_graph_nodes, 0));

  for(size_t i = 0; i < crn.number_species; i++){
    index_vectors.push_back(pow(2,i));
    max_index_vectors.push_back(pow(2,crn.number_species) - crn.number_species + i);
  }
  uint64_t counter = 0;
  uint64_t counter2 = 0;
  uint64_t counter3 = 0;
  uint64_t counter4 = 0;

  /* calculate transformation matrices and check consistency */
  while(index_vectors[0] <= max_index_vectors[0]){
    counter++;
    if(!(counter % 10000)){std::cout << "number transformations: " << counter << std::endl;}

    // calculate new transformation matrix from index vectors
    for(size_t i = 0; i < crn.number_species; i++){
      for(size_t j = 0; j < crn.number_species; j++){
        transformation_matrix_doubles[i][j] = possible_vectors[index_vectors[i]-1][j];
        transformation_matrix_ints[i][j] = possible_vectors[index_vectors[i]-1][j];
      }
    }

    if (Math::calculate_determinant(transformation_matrix_doubles, crn.number_species) != 0){
      counter2++;
      if(!(counter2 % 10000)){std::cout << "number dets != 0 : " << counter2 << std::endl;}

      /*
       * the value of (i,j) may be one of {-1,0,1} (0 means unchecked) and describes the sign that all checked PATHS between the two nodes are marked with
       * if a newly found path has a different sign, the graph is not consistent
      */
      std::vector<std::vector<int>> path_sign_matrix (crn.number_species, std::vector<int>(number_graph_nodes, 0));

      for(size_t i = 0; i < crn.number_species; i++){
        for(size_t j = 0; j < crn.number_species; j++){
          transformation_matrix(i,j) = transformation_matrix_ints[i][j];
        }
      }

      GiNaC::matrix jacobian_c_transformed = transformation_matrix.mul(jacobian_c);//.mul(transformation_matrix.inverse());
      GiNaC::matrix jacobian_rr_transformed;
        jacobian_rr_transformed = transformation_matrix.mul(jacobian_rr);

      // build sign matrices and evaluate their consistency
      std::vector<std::vector<int>> transformed_sign_matrix_c (crn.number_species, std::vector<int>(crn.number_species, 0));
      std::vector<std::vector<int>> transformed_sign_matrix_rr (crn.number_species, std::vector<int>(crn.number_reaction_rates, 0));

      bool is_consistent = 1;
      for(size_t i = 0; i < crn.number_species; i++){
        calculate_consistency(jacobian_c_transformed, transformed_sign_matrix_c, crn.number_species, is_consistent, i);
          calculate_consistency(jacobian_rr_transformed, transformed_sign_matrix_rr, crn.number_reaction_rates, is_consistent, i);
      }

      // build graph and evaluate graph consistency
      if(is_consistent){
        counter3++;
        if(!(counter3 % 10000)){std::cout << "number consistent, but not yet graph consistent: " << counter3 << std::endl;}

        // write values to graph
        for(size_t i = 0; i < crn.number_species; i++){
          size_t j = 0;
          while(j < crn.number_species){
            if(i != j){
              graph[i][j] = transformed_sign_matrix_c[i][j];
            }
            j++;
          }
          while(j < number_graph_nodes){
            graph[i][j] = transformed_sign_matrix_rr[i][j - crn.number_species];
            j++;
          }
        }

        calculate_graph_consistency(crn.number_species, number_graph_nodes, path_sign_matrix, graph, is_consistent);
      }

      // if transformation sign matrices and graph are consistent, save results
      if(is_consistent){
        counter4++;
        if(!(counter4 % 10000)){std::cout << "number consistent and graph consistent: " << counter4 << std::endl;}

        save_consistent_transformation(time_derivatives, jacobian_c, jacobian_rr, transformation_matrix, jacobian_c_transformed, jacobian_rr_transformed,
          transformed_sign_matrix_c, transformed_sign_matrix_rr, graph, path_sign_matrix);
      }
    }

    // calculate next index vector
    calculate_next_index_vector(index_vectors, max_index_vectors, crn.number_species);
  }

  std::cout << "final number transformations: " << counter << std::endl;
  std::cout << "final number transformations with determinants not 0 : " << counter2 << std::endl;
  std::cout << "final number consistent, but not yet graph consistent transformations: " << counter3 << std::endl;
  std::cout << "final number consistent and graph consistent transformations: " << counter4 << std::endl;
  std::cout << "done." << std::endl;
}


// TODO: write this in a more comprehensible way
void MonotoneDependenciesCalculator::calculate_next_index_vector(std::vector<int>& index_vectors, const std::vector<int>& max_index_vectors,
  size_t number_species){
  index_vectors[number_species - 1]++;
  int i = number_species - 1;
  while(i >= 0){
    if(index_vectors[i] > max_index_vectors[i]){
      i--;
      if(i >= 0){
        index_vectors[i] = index_vectors[i] + 1;
      }
    }else{
      for(size_t j = i + 1; j < number_species; j++){
        index_vectors[j] = index_vectors[i]+j-i;
        if(index_vectors[j] < pow(2,j)){
          index_vectors[j] = pow(2,j);
          for(size_t k = j + 1; k < number_species; k++){
            index_vectors[k] = pow(2,k);
          }
          break;
        }
      }
      break;
    }
  }
}

// TODO: write this in a more comprehensible way
void MonotoneDependenciesCalculator::calculate_consistency(GiNaC::matrix& jacobian_transformed,
  std::vector<std::vector<int>>& transformed_sign_matrix, size_t number_species, bool& is_consistent, int i){
  if(is_consistent){
    for(size_t j = 0; j < number_species; j++){
      if(jacobian_transformed(i,j) == 0){
      }else if(jacobian_transformed(i,j) == 1){
        transformed_sign_matrix[i][j] = 1;
      }else if(jacobian_transformed(i,j) == -1){
        transformed_sign_matrix[i][j] = -1;
      }else{
        std::ostringstream ss;
        ss << jacobian_transformed(i,j);
        std::string str (ss.str());
        // no - found => positive
        if(str.find('-') == std::string::npos){
          transformed_sign_matrix[i][j] = 1;
        // - minus found, but not at first place => inconsistent
        }else if(str.find('-') != 0){
          is_consistent = 0;
          return;
        // - found at first place and + found => inconsistent
        }else if(str.find('-') != std::string::npos && str.find('+') != std::string::npos){
          is_consistent = 0;
          return;
        // - found at first place and + not found => negative
        }else{
          transformed_sign_matrix[i][j] = -1;
        }
      }
    }
  }
  return;
}


// TODO: this uses "brute force" to determine consistency. there should be a more efficient way to do this.
void MonotoneDependenciesCalculator::calculate_graph_consistency(size_t number_species, size_t number_graph_nodes,
  std::vector<std::vector<int>>& path_sign_matrix, const std::vector<std::vector<int>>& graph, bool& is_consistent){
  // for each two nodes, check all possible paths for their signs and thus discover possible inconsistencies
  // start with iterating through all nodes
  for(size_t i = 0; i < number_graph_nodes; i++){
    std::vector<size_t> path;
    path.push_back(i);
    calculate_path_signs(number_species, number_graph_nodes, path_sign_matrix, graph, is_consistent, path, 1);
    // if no inconsistencies are found, continue with next node
    if(!is_consistent){break;}
  }
}

void MonotoneDependenciesCalculator::calculate_path_signs(size_t number_species, size_t number_graph_nodes,
  std::vector<std::vector<int>>& path_sign_matrix, const std::vector<std::vector<int>>& graph, bool& is_consistent,
  std::vector<size_t>& path, int sign){
  size_t path_end = path.size() - 1;
  for(size_t i = 0; i < number_species; i++){ // search for possible successor nodes
    if(i != path[path_end] && graph[i][path[path_end]] != 0){ // if an edge exists between the last and the new node (=i) within the path
      size_t cycle_index = 0; // initialize
      for(size_t j = 0; j < path.size(); j++){ // search for a possible cycle in the path
        if(i == path[j]){ // if the new node already exists in the path
          cycle_index = j; // index of the node which appears now twice
          break; // go out of this loop, another cycle is not possible
        }
      }

      if(cycle_index == 0){ // if there is no cycle in the path
        if(path_sign_matrix[i][path[0]] != 0){
          if(path_sign_matrix[i][path[0]] == (sign * graph[i][path[path_end]])){ // if there is already a path from node 'path[0]' to node 'i'
            path.push_back(i);
            calculate_path_signs(number_species, number_graph_nodes, path_sign_matrix, graph, is_consistent, path, sign * graph[i][path[path_end]]); // if the existing path sign equals the sign of the new path
            path.pop_back(); // add the new node to the, update the sign and call this function again
          }else{
            is_consistent = 0; return; // else there are paths with different signs => close this function
          }
        }else{
          path_sign_matrix[i][path[0]] = sign * graph[i][path[path_end]]; // else set the sign for paths between these nodes
          path.push_back(i);
          calculate_path_signs(number_species, number_graph_nodes, path_sign_matrix, graph, is_consistent, path, sign * graph[i][path[path_end]]); // add the new node to the, update the sign and call this function again
          path.pop_back();
        }
      }else{ // else there is a cycle which sign has to be determined
        int sign_cycle = 1; // start with the sign +1 and update it with each edge
        for(size_t j = cycle_index; j < path_end; j++){ // for each edge in the path and the cycle
          sign_cycle = sign_cycle * graph[path[j + 1]][path[j]]; // update the sign of the cycle
        }
        sign_cycle = sign_cycle * graph[i][path[path_end]]; // add the sign for the edge which closes the cycle
        if(sign_cycle == -1){ // a negative cycle
          is_consistent = 0; return; // implies paths with different signs => close this function
        }
      }
    }
    if(!is_consistent){return;} // check if one of the previous calls of this function found inconsistent signs => close this function

  }
}


void MonotoneDependenciesCalculator::save_consistent_transformation(const GiNaC::matrix time_derivatives, const GiNaC::matrix& jacobian_c,
  const GiNaC::matrix& jacobian_rr,
  const GiNaC::matrix& transformation_matrix, const GiNaC::matrix& jacobian_c_transformed, const GiNaC::matrix& jacobian_rr_transformed,
  const std::vector<std::vector<int>>& transformed_sign_matrix_c, const std::vector<std::vector<int>>& transformed_sign_matrix_rr,
  const std::vector<std::vector<int>>& graph, const std::vector<std::vector<int>>& path_sign_matrix){
  std::stringstream ss;

  // species
  ss << "species:" << std::endl;
  for(size_t i = 0; i < crn.number_species - 1; i++){
    ss << crn.species_list[i] << "  ";
  }
  ss << crn.species_list[crn.number_species - 1] << std::endl << std::endl;

  // reaction rates
  ss << "reaction rates:" << std::endl;
  for(size_t i = 0; i < crn.number_reaction_rates - 1; i++){
    ss << crn.reaction_rate_list[i] << "  ";
  }
  ss << crn.reaction_rate_list[crn.number_reaction_rates - 1] << std::endl << std::endl;

  // reactions
  ss << "reactions:" << std::endl;
  for(size_t i = 0; i < crn.number_reactions; i++){
    ss << "r" << i << " = " << crn.reaction_list[i] << std::endl;
  }
  ss << std::endl;

  // time derivatives
  ss << "time derivatives:" << std::endl;
  for(size_t i = 0; i < crn.number_species; i++){
    ss << "d/dt " << crn.species_list[i] << " = " << time_derivatives(i,0) << std::endl;
  }
  ss << std::endl;

  // Jacobian with respect to the species
  ss << "Jacobian with respect to the species:" << std::endl;
  for(size_t i = 0; i < crn.number_species; i++){
    for(size_t j = 0; j < crn.number_species - 1; j++){
      ss << jacobian_c(i,j) << "\t";
    }
    ss << jacobian_c(i, crn.number_species - 1) << std::endl;;
  }
  ss << std::endl;

  // Jacobian with respect to the reaction rates
    ss << "Jacobian with respect to the reaction rates:" << std::endl;
    for(size_t i = 0; i < crn.number_species; i++){
      for(size_t j = 0; j < crn.number_reaction_rates - 1; j++){
        ss << jacobian_rr(i,j) << "\t";
      }
      ss << jacobian_rr(i, crn.number_reaction_rates - 1) << std::endl;;
    }
    ss << std::endl;

  // Transformation matrix
  ss << "Transformation matrix:" << std::endl;
  for(size_t i = 0; i < crn.number_species; i++){
    for(size_t j = 0; j < crn.number_species - 1; j++){
      ss << transformation_matrix(i,j) << "\t";
    }
    ss << transformation_matrix(i, crn.number_species - 1) << std::endl;;
  }
  ss << std::endl;

  // Inverse of transformation matrix
  ss << "Inverse Transformation matrix:" << std::endl;
  GiNaC::matrix transformation_matrix_inverse = transformation_matrix.inverse();
  for(size_t i = 0; i < crn.number_species; i++){
    for(size_t j = 0; j < crn.number_species - 1; j++){
      ss << transformation_matrix_inverse(i,j) << "\t";
    }
    ss << transformation_matrix_inverse(i, crn.number_species - 1) << std::endl;;
  }
  ss << std::endl;

  // Species in the new coordinates
  GiNaC::matrix new_coordinates (crn.number_species, 1);
  for(size_t i = 0; i < crn.number_species; i++){
    new_coordinates(i,0) = crn.species_list[i];
  }
  new_coordinates = transformation_matrix.mul(new_coordinates);
  ss << "Species in the new coordinates:" << std::endl;
  for(size_t i = 0; i < crn.number_species - 1; i++){
    ss << new_coordinates(i,0) << "  ";
  }
  ss << new_coordinates(crn.number_species - 1,0) << std::endl << std::endl;

  // Time-derivatives in the new coordinates
/*  GiNaC::matrix tmp = crn.stoichiometric.mul(Math::jacobian(crn.reaction_list, crn.species_list));
  GiNaC::matrix new_coordinates_time_derivatives = transformation_matrix.mul(tmp);
  ss << "Time-derivatives in the new coordinates:" << std::endl;
  for(size_t i = 0; i < crn.number_species; i++){
    for(size_t j = 0; j < crn.number_species - 1; j++){
      ss << new_coordinates_time_derivatives(i,j) << "\t";
    }
    ss << new_coordinates_time_derivatives(i, crn.number_species - 1) << std::endl;;
  }
  ss << std::endl;
*/
  // Jacobian with respect to the species in the new coordinates
  ss << "Jacobian with respect to the species in the new coordinates:" << std::endl;
  for(size_t i = 0; i < crn.number_species; i++){
    for(size_t j = 0; j < crn.number_species - 1; j++){
      ss << jacobian_c_transformed(i,j) << "\t";
    }
    ss << jacobian_c_transformed(i, crn.number_species - 1) << std::endl;;
  }
  ss << std::endl;

  // Jacobian with respect to the reaction rates in the new coordinates
    ss << "Jacobian with respect to the reaction rates in the new coordinates:" << std::endl;
    for(size_t i = 0; i < crn.number_species; i++){
      for(size_t j = 0; j < crn.number_reaction_rates - 1; j++){
        ss << jacobian_rr_transformed(i,j) << "\t";
      }
      ss << jacobian_rr_transformed(i, crn.number_reaction_rates - 1) << std::endl;;
    }
    ss << std::endl;

  // Signs of the entries of the Jacobian with respect to the species in the new coordinates
  ss << "Signs of the entries of the Jacobian with respect to the species in the new coordinates:" << std::endl;
  for(size_t i = 0; i < crn.number_species; i++){
    for(size_t j = 0; j < crn.number_species - 1; j++){
      ss << transformed_sign_matrix_c[i][j] << "\t";
    }
    ss << transformed_sign_matrix_c[i][crn.number_species - 1] << std::endl;
  }
  ss << std::endl;


  // Signs of the entries of the Jacobian with respect to the reaction rates in the new coordinates
    ss << "Signs of the entries of the Jacobian with respect to the reaction rates in the new coordinates:" << std::endl;
    for(size_t i = 0; i < crn.number_species; i++){
      for(size_t j = 0; j < crn.number_reaction_rates - 1; j++){
        ss << transformed_sign_matrix_rr[i][j] << "\t";
      }
      ss << transformed_sign_matrix_rr[i][crn.number_reaction_rates - 1] << std::endl;
    }
    ss << std::endl;

  // Graph matrix
  ss << "Graph matrix:" << std::endl;
  for(size_t i = 0; i < crn.number_species; i++){
    ss << "\t" << new_coordinates(i,0);
  }
  for(size_t i = 0; i < crn.number_reaction_rates; i++){
    ss << "\t" << crn.reaction_rate_list[i];
  }
  ss << std::endl;
  for(size_t i = 0; i < crn.number_species; i++){
    ss << new_coordinates(i,0);
    for(size_t j = 0; j < crn.number_species + crn.number_reaction_rates; j++){
      ss << "\t" << graph[i][j];
    }
    ss << std::endl;
  }
  ss << std::endl;

  // Matrix of the signs for all paths between two nodes
  ss << "Matrix of the signs of all paths between two nodes" << std::endl;
  for(size_t i = 0; i < crn.number_species; i++){
    ss << "\t" << new_coordinates(i,0);
  }
  for(size_t i = 0; i < crn.number_reaction_rates; i++){
    ss << "\t" << crn.reaction_rate_list[i];
  }
  ss << std::endl;
  for(size_t i = 0; i < crn.number_species; i++){
    ss << new_coordinates(i,0);
    for(size_t j = 0; j < crn.number_species + crn.number_reaction_rates; j++){
      ss << "\t" << path_sign_matrix[i][j];
    }
    ss << std::endl;
  }
  ss << std::endl;

  // push string to log vector
  log_vector.push_back(ss.str());



  // string for graph file
  std::stringstream graph_ss;

  // graph name
  graph_ss << "digraph " << crn.model_name << "{" << std::endl << std::endl;

  // species nodes
  graph_ss << "node [shape=ellipse]; ";
  for(size_t i = 0; i < crn.number_species; i++){
    graph_ss << "\"" << new_coordinates(i,0) << "\"; ";
  }
  graph_ss << std::endl << std::endl;

  // links between species
  for(size_t i = 0; i < crn.number_species; i++){
    for(size_t j = 0; j < crn.number_species; j++){
      if(graph[i][j] == 1){
        graph_ss << "\"" << new_coordinates(i,0) << "\" -> \"" << new_coordinates(j,0) << "\" [label=\"+\"];" << std::endl;
      }else if(graph[i][j] == -1){
        graph_ss << "\"" << new_coordinates(i,0) << "\" -> \"" << new_coordinates(j,0) << "\" [label=\"-\"];" << std::endl;
      }
    }
  }
  graph_ss << std::endl;

  // reaction rate nodes
  graph_ss << "node [shape=box]; ";
  for(size_t i = 0; i < crn.number_reaction_rates; i++){
    graph_ss << "\"" << crn.reaction_rate_list[i] << "\";";
  }
  graph_ss << std::endl << std::endl;

  // links between reaction rates and species
  for(size_t i = 0; i < crn.number_species; i++){
    for(size_t j = 0; j < crn.number_reaction_rates; j++){
      if(graph[i][i + j] == 1){
        graph_ss << "\"" << crn.reaction_rate_list[j] << "\" -> \"" << new_coordinates(i,0) << "\" [label=\"+\"];" << std::endl;
      }else if(graph[i][i + j] == -1){
        graph_ss << "\"" << crn.reaction_rate_list[j] << "\" -> \"" << new_coordinates(i,0) << "\" [label=\"-\"];" << std::endl;
      }
    }
  }
  graph_ss << std::endl;

  // graph label and font size
  graph_ss << "label = \""<< crn.model_name << "\";" << std::endl << "fontsize=20;" << std::endl << "}" << std::endl;

  // push string to graph vector
  graph_vector.push_back(graph_ss.str());
}


void MonotoneDependenciesCalculator::log(const std::string dirname)const{
  for(size_t i = 0; i < log_vector.size(); i++){
    std::string leading_zeros = "";
    if(i){
      for(size_t j = pow(10, number_leading_zeros); j > i ; j /= 10){
        leading_zeros += "0";
      }
    }else{
      for(size_t j = 0; j < number_leading_zeros; j++){
        leading_zeros += "0";
      }
    }

    // create names for data file and graph file
    std::stringstream ss, graph_ss;
    if(dirname[dirname.size() - 1] == '/'){
      ss << dirname << crn.model_name << "_output_" << leading_zeros << i << ".txt";
      graph_ss << dirname << crn.model_name << "_graph_" << leading_zeros << i << ".doc";
    }else{
      ss << dirname << "/" << crn.model_name << "_output_" << leading_zeros << i << ".txt";
      graph_ss << dirname << "/" <<  crn.model_name << "_graph_" << leading_zeros << i << ".dot";
    }

    // write to files
    std::ofstream outputfile (ss.str());
    outputfile << log_vector[i] << std::endl;
    outputfile.close();

    std::ofstream graphfile (graph_ss.str());
    graphfile << graph_vector[i] << std::endl;
    graphfile.close();
  }
}

