#pragma once

#include <string>
#include <vector>
#include <ginac/ginac.h>

class CRN
{
public:
  /// constructor
  /**
   * takes a filename of a text file containing all necessary information about the CRN and converts the input to processable state.
   * for information on how the input file should look like, see the input_example.txt file and the README.txt file.
  */
  CRN(const std::string& filename, const std::string& model_name);

  /// name of the file containing all input information
  const std::string filename;

  /// name of the crn model
  /**
   * this name will be used for marking the output data
  */
  const std::string model_name;

  /// number of reaction rates
  size_t number_reaction_rates = 0;
  /// number of species
  size_t number_species = 0;
  /// number of reactions
  size_t number_reactions = 0;

  /// list of all reaction rates
  std::vector<GiNaC::symbol> reaction_rate_list;
  /// list of all species
  std::vector<GiNaC::symbol> species_list;
  /// list of all reactions
  std::vector<GiNaC::ex> reaction_list;
  /// stoichiometric matrix of the crn
  GiNaC::matrix stoichiometric;

  /// prints all reaction rates to stdout
  void print_rr()const;
  /// prints all species to stdout
  void print_s()const;
  /// prints all reactions to stdout
  void print_r()const;
  /// prints raction rates, species, reactions and the stoichiometric matrix to stdout
  void print_all()const;

  /// string in input file that signals that the parsing of reaction rates starts in the next line
  std::string reaction_rates_signal = "reaction rates start";
  /// string in input file that signals that the parsing of species starts in the next line
  std::string species_signal = "species start";
  /// string in input file that signals that the parsing of reactions starts in the next line
  std::string reactions_signal = "reactions start";
  /// string in input file that signals that the parsing ends with this line
  std::string end_signal = "end";
  /// string in input file that signals that no reaction constant exists
  std::string no_constant_signal = "-";
private:
  /// converts a string containing a reaction to a processable GiNaC expression
  GiNaC::ex stringToReaction(std::string reaction_string, const std::string forward_reaction_constant,
    const std::string backward_reaction_constant);

  /// calculates a stoichiometrix matrix from the reactions read from input
  void calculateStoichiometric();

  /// replaces all occurences of a string in a string with a string */
  void replaceAll(std::string& str, const std::string& from, const std::string& to)const;

  /**
   * maps are needed for building the reactions and the stoichiometric matrix from the reaction rates and species.
   * for calculations later however, vectors are used due to better performance.
  */
  std::map<std::string,size_t> reaction_rate_map;
  /**
   * maps are needed for building the reactions and the stoichiometric matrix from the reaction rates and species.
   * for calculations later however, vectors are used due to better performance.
  */
  std::map<std::string,size_t> species_map;
  /// vector containing all reactions as a string
  /** this vector is needed for calculating the stoichiometric matrix */
  std::vector<std::string> reaction_string_vector;
};
