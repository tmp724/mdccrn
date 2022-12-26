#include "crn.hpp"
#include "monotonedependenciescalculator.hpp"

#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>

struct stat info;

int main(int argc, char **argv){
  try{
    if(argc < 4 || argc > 5){
      throw "Bad number of arguments. Please see the README file for usage instructions.";
    }
    if(stat(argv[1], &info)){
      throw "Could not access input file.";
    }
    if(stat(argv[2], &info) && (!(info.st_mode & S_IFDIR))){
      throw "Could not open output directory.";
    }
    const CRN crn1 (argv[1], argv[3]);
    size_t mode = 0;
    if(argc == 5){
      mode = atoi(argv[4]);
    }
    MonotoneDependenciesCalculator mdc(crn1, mode);
    mdc.run();
    mdc.log(argv[2]);
  }catch(const std::exception& e){
    std::cerr << e.what() << std::endl;
  }catch(const char* e){
    std::cerr << e << std::endl;
  }
  return 0;
}

