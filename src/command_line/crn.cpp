#include "crn.hpp"
#include "math.hpp"

#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <iostream>

/*
 * TODO:
 * - The way the information is extracted from file is too specific; a more abstract approach would improve flexibility
*/
CRN::CRN(const std::string& filename, const std::string& model_name):filename(filename),model_name(model_name){
/*  this->filename = filename;
  this->model_name = model_name;*/
/*  try{
    file.open(filename);
  }catch(std::ios_base::failure& e){
    std::cout << "Opening or reading input file wasn't possible." << std::endl;
  }*/
  try{
    std::ifstream file(filename);
    if(!file){
      throw std::ios::failure("Error opening file!");
    }

    std::string line;
    int marker = 0;

    // read file line by line
    while(getline(file,line)){
      replaceAll(line, "\t", " ");
      std::string::iterator start = line.begin();
      std::string::iterator end = line.end();
      std::string::iterator current, next, reaction_end;
      const char c = ' ';

      // if the patterns are matched, set the marker so that the following lines get parsed accordingly
      if(line.find(end_signal) != std::string::npos){
        marker = 0; continue;
      }else if(line.find(reaction_rates_signal) != std::string::npos){
        marker = 1; continue;
      }else if(line.find(species_signal) != std::string::npos){
        marker = 2; continue;
      }else if(line.find(reactions_signal) != std::string::npos){
        marker = 3; continue;
      }

      // parse reaction rates
      if(marker == 1){
        next = std::find(start,end,c);
        reaction_rate_list.push_back(GiNaC::symbol(std::string(start, next)));
        reaction_rate_map[std::string(start,next)] = number_reaction_rates;
        number_reaction_rates++;

      // parse species
      }else if(marker == 2){
        next = std::find(start,end,c);
        std::string word = std::string(start, next);
        species_list.push_back(GiNaC::symbol(word));
        species_map[word] = number_species;
        number_species++;

      // parse reactions
      }else if(marker == 3){
        std::string reaction_string;
        std::string forward_reaction_constant = "0"; std::string backward_reaction_constant = "0";
        const char quotes = '\"';
        current = std::find(start, end, quotes);
        reaction_end = std::find(current + 1, end, quotes);
        reaction_string = std::string(current + 1, reaction_end);
        reaction_string_vector.push_back(reaction_string);
        current = reaction_end + 1;
        size_t number_words_after_reaction = 0;
        while(next != end){
          next = std::find(current, end, c);
          std::string word = std::string(current, next);
          if(word == ""){
            //do nothing
          }else if(number_words_after_reaction == 0){
            if(word != no_constant_signal){
              forward_reaction_constant = word;
            }
            number_words_after_reaction++;
          }else if (number_words_after_reaction == 1){
            if(word != no_constant_signal){
              backward_reaction_constant = word;
            }
            number_words_after_reaction++;
          }
          current = next + 1;
        }

        reaction_list.push_back(stringToReaction(reaction_string, forward_reaction_constant, backward_reaction_constant));
        number_reactions++;
      }
    }

    // calculate stoichiometric matrix from reactions
    calculateStoichiometric();

    file.close();
  }catch(const std::exception& e){
    std::cout << e.what() << std::endl;
  }

}


void CRN::calculateStoichiometric(){
  const char c = ' ';
  GiNaC::matrix tmp_matrix(number_species, number_reactions);
  for(size_t i = 0; i < number_reactions; i++){
    for(size_t j = 0; j < number_species; j++){
      std::string reaction_string = reaction_string_vector[i];
      std::string::iterator next = reaction_string.begin();
      std::string::iterator current = reaction_string.begin();
      std::string::iterator end = reaction_string.end();
      std::string word, lastWord;
      bool reaction_sign = 0;
      int tmp_sum = 0;
      while(next != end){
        next = std::find(current, end, c);
        word = std::string(current, next);

        // word is reaction sign
        if(word == "=>" || word == "<=>"){
          reaction_sign = 1;

        // if word is species, test for potential coefficients
        }else if(!word.empty() && species_map.find(word) != species_map.end() && j == species_map[word]){
          if(reaction_sign && !lastWord.empty() && lastWord.find_first_not_of("0123456789") == std::string::npos){
            tmp_sum -= std::stoi(lastWord);
          }else if(reaction_sign){
            tmp_sum -= 1;
          }else if(!lastWord.empty() && lastWord.find_first_not_of("0123456789") == std::string::npos){
            tmp_sum += std::stoi(lastWord);
          }else{
            tmp_sum += 1;
          }
        }
        lastWord = word;
        current = next + 1;
      }
      if(tmp_sum < 0){
        tmp_matrix(j,i) = 1;
      }else if(tmp_sum > 0){
        tmp_matrix(j,i) = -1;
      }else{
        tmp_matrix(j,i) = 0;
      }
    }
  }
  stoichiometric = tmp_matrix;
}


GiNaC::ex CRN::stringToReaction(std::string reaction_string, const std::string forward_reaction_constant,
  const std::string backward_reaction_constant){
  std::string::iterator next = reaction_string.begin();
  std::string::iterator current = reaction_string.begin();
  std::string::iterator end = reaction_string.end();
  std::string word;
  const char c = ' ';
  std::vector<GiNaC::ex> expressionBuffer1;
  std::vector<GiNaC::ex> expressionBuffer2;
  bool reaction_sign = 0;
  GiNaC::ex reaction;
  GiNaC::ex forwardReaction, backwardReaction;

  while(next != end){
    next = std::find(current, end, c);
    word = std::string(current, next);

    // word is reaction sign
    if(word == "=>" || word == "<=>"){
      reaction_sign = 1;

    // word is species
    }else if(!word.empty() && species_map.find(word) != species_map.end()){
      GiNaC::ex tmp = species_list[species_map[word]];
      if(reaction_sign){
        expressionBuffer2.push_back(tmp);
      }else{
        expressionBuffer1.push_back(tmp);
      }

    // word is coefficient
    }else if(!word.empty() && word.find_first_not_of("0123456789") == std::string::npos){
      if(reaction_sign){
        expressionBuffer2.push_back(std::stoi(word));
      }else{
        expressionBuffer1.push_back(std::stoi(word));
      }
    }
    current = next + 1;
  }

  if(forward_reaction_constant != "0"){
    forwardReaction = reaction_rate_list[reaction_rate_map[forward_reaction_constant]];
    for(size_t i = 0; i < expressionBuffer1.size() - 1; i++){
      forwardReaction = forwardReaction * expressionBuffer1[i];
    }
    forwardReaction = forwardReaction * expressionBuffer1[expressionBuffer1.size()-1];
  }

  if(backward_reaction_constant != "0"){
    backwardReaction = reaction_rate_list[reaction_rate_map[backward_reaction_constant]];
    for(size_t i = 0; i < expressionBuffer2.size() - 1; i++){
      backwardReaction = backwardReaction * expressionBuffer2[i];
    }
    backwardReaction = backwardReaction * expressionBuffer2[expressionBuffer2.size()-1];
  }

  reaction = forwardReaction - backwardReaction;

  return reaction;
}


void CRN::replaceAll(std::string& str, const std::string& from, const std::string& to)const{
  if(from.empty()){
    return;
  }
  size_t start_pos = 0;
  while((start_pos = str.find(from, start_pos)) != std::string::npos) {
    str.replace(start_pos, from.length(), to);
    start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
  }
}


void CRN::print_rr()const{
  std::cout << "reaction rate list: " << std::endl << "[";
  for (size_t i = 0; i < reaction_rate_list.size() - 1; i++){
    std::cout << reaction_rate_list[i] << ", ";
  }
  std::cout << reaction_rate_list[reaction_rate_list.size() - 1] << "]" << std::endl;
}


void CRN::print_s()const{
  std::cout << "species list: " << std::endl << "[";
  for (size_t i = 0; i < species_list.size() - 1; i++){
    std::cout << species_list[i] << ", ";
  }
  std::cout << species_list[species_list.size() - 1] << "]" << std::endl;
}


void CRN::print_r()const{
  std::cout << "reaction list: " << std::endl << "[";
  for (size_t i = 0; i < reaction_list.size() - 1; i++){
    std::cout << reaction_list[i] << ", ";
  }
  std::cout << reaction_list[reaction_list.size() - 1] << "]" << std::endl;
}


void CRN::print_all()const{
  print_rr();
  print_s();
  print_r();
  std::cout << "stoichiometric matrix: " << std::endl << stoichiometric << std::endl;
}
