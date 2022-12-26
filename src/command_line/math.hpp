#pragma once

#include <ginac/ginac.h>

class Math
{
public:
  /// calculates the jacobian matrix of given sets of functions and variables
  static GiNaC::matrix jacobian(const std::vector<GiNaC::ex>& functions, const std::vector<GiNaC::symbol>& variables){
    size_t number_functions = functions.size();
    size_t number_variables = variables.size();
    GiNaC::matrix jacobian_matrix(number_functions, number_variables);
    for(size_t i = 0; i < number_functions; i++){
      for(size_t j = 0; j < number_variables; j++){
        jacobian_matrix(i,j) = functions[i].diff(variables[j]);
      }
    }
    return jacobian_matrix;
  }

  /// returns a matrix of 0s and 1s
  /** takes a number and returns a matrix with size 2^number-1 where each column
   * represents a unique binary number within the range of 1 to 2^number
  */
  static std::vector<std::vector<int>> binary_matrix(size_t number_rows){
    size_t number_columns = pow(2,number_rows) - 1;
    std::vector<std::vector<int>> possible_vectors(number_columns, std::vector<int>(number_rows));
    for(size_t i = number_columns; i > 0; i--){
      size_t tmp = i;
      for(size_t j = 0; j < number_rows; j++){
        possible_vectors[i-1][j] = tmp % 2;
        tmp = tmp/2;
      }
    }

    return possible_vectors;
  }

  /// calculates the determinant of a matrix
  /** returns the determinant of a given square matrix of given size.
   * The chosen algorithm is rather efficient, but for now only works for matrices with size 1-10.
  */
  //TODO: default: throw exception when size > 10
  static int calculate_determinant(const std::vector<std::vector<double>>& m, size_t size){
    // for rounding determinant, as numeric calculation with doubles isn't precise
    double epsilon = 0.1;
    double det;
    std::vector<double> trace_vector (size);
    std::vector<std::vector<double>> T_p (size, std::vector<double>(size));
    std::vector<std::vector<double>> T_temp (size, std::vector<double>(size));

    for(size_t i = 0; i < size; i++) {
      trace_vector[i] = 0;                                                // initialize the trace variables with zeroes
      for(size_t j = 0; j < size; j++) {
        T_p[i][j] = 0;
        T_temp[i][j] = m[i][j];                                             // initialize T_temp as copy of T
      }
    }

    for(size_t i = 0; i < size; i++) {
      trace_vector[0] = trace_vector[0] + m[i][i];
    }
    // calculate
    for(size_t k = 1; k < size-1; k++) {
      // calculate T^k=T^(k-1)*T (T_p=T_temp*T)
      for(size_t i = 0; i < size; i++) {
        for(size_t j = 0; j < size; j++) {
          for(size_t o = 0; o < size; o++) {
            T_p[i][j] = T_p[i][j] + T_temp[i][o]*m[o][j];
          }
        }
      }
      // calculate the trace of T^k (trace of T_p)
      for(size_t i = 0; i < size; i++) {
        trace_vector[k] = trace_vector[k] + T_p[i][i];
      }
      // copy T_p into T_temp and reset T_p
      for(size_t i = 0; i < size; i++) {
        for(size_t j = 0; j < size; j++) {
          T_temp[i][j] = T_p[i][j];                                 // copy T_p into T_temp
          T_p[i][j] = 0;                                           // reset the entries of T_p with zeroes
        }
      }
    }
    // calculate only the trace of T^n=T^(n-1)*T (T_p=T_temp*T)
    for(size_t i = 0; i < size; i++){
      for(size_t o = 0; o < size; o++){
        trace_vector[size-1] = trace_vector[size-1] + T_temp[i][o]*m[o][i];
      }
    }

    switch(size){
      case 1 : det = -trace_vector[0]; break;
      case 2 : det = (trace_vector[0]*trace_vector[0] - trace_vector[1])/2; break;
      case 3 : det = -trace_vector[0]*trace_vector[0]*trace_vector[0]/6 + trace_vector[0]*trace_vector[1]/2 - trace_vector[2]/3; break;
      case 4 : det = pow(trace_vector[0],4)/24 - trace_vector[0]*trace_vector[0]*trace_vector[1]/4 + trace_vector[0]*trace_vector[2]/3 + trace_vector[1]*trace_vector[1]/8 - trace_vector[3]/4; break;
      case 5 : det = trace_vector[0]*trace_vector[3]/4 - trace_vector[4]/5 + trace_vector[1]*trace_vector[2]/6 - trace_vector[0]*trace_vector[1]*trace_vector[1]/8 - trace_vector[0]*trace_vector[0]*trace_vector[2]/6 + trace_vector[0]*trace_vector[0]*trace_vector[0]*trace_vector[1]/12 - pow(trace_vector[0],5)/120; break;
      case 6 : det = pow(trace_vector[0],6)/720 - (pow(trace_vector[0],4)*trace_vector[1])/48 + (pow(trace_vector[0],3)*trace_vector[2])/18 + (trace_vector[0]*trace_vector[0]*trace_vector[1]*trace_vector[1])/16 - (trace_vector[3]*trace_vector[0]*trace_vector[0])/8 - (trace_vector[0]*trace_vector[1]*trace_vector[2])/6 + (trace_vector[4]*trace_vector[0])/5 - pow(trace_vector[1],3)/48 + (trace_vector[3]*trace_vector[1])/8 + trace_vector[2]*trace_vector[2]/18 - trace_vector[5]/6; break;
      case 7 : det = (trace_vector[0]*trace_vector[5])/6 - (pow(trace_vector[0],3)*trace_vector[1]*trace_vector[1])/48 - trace_vector[6]/7 + (trace_vector[1]*trace_vector[4])/10 + (trace_vector[2]*trace_vector[3])/12 + (trace_vector[0]*pow(trace_vector[1],3))/48 - (trace_vector[0]*trace_vector[2]*trace_vector[2])/18 - (trace_vector[1]*trace_vector[1]*trace_vector[2])/24 - (trace_vector[0]*trace_vector[0]*trace_vector[4])/10 + (pow(trace_vector[0],3)*trace_vector[3])/24 - (pow(trace_vector[0],4)*trace_vector[2])/72 + (pow(trace_vector[0],5)*trace_vector[1])/240 - pow(trace_vector[0],7)/5040 + (trace_vector[0]*trace_vector[0]*trace_vector[1]*trace_vector[2])/12 - (trace_vector[0]*trace_vector[1]*trace_vector[3])/8; break;
      case 8 : det = pow(trace_vector[0],8)/40320 - (pow(trace_vector[0],6)*trace_vector[1])/1440 + (pow(trace_vector[0],5)*trace_vector[2])/360 + (pow(trace_vector[0],4)*pow(trace_vector[1],2))/192 - (pow(trace_vector[0],4)*trace_vector[3])/96 - (pow(trace_vector[0],3)*trace_vector[1]*trace_vector[2])/36 + (trace_vector[4]*pow(trace_vector[0],3))/30 - (pow(trace_vector[0],2)*pow(trace_vector[1],3))/96 + (pow(trace_vector[0],2)*trace_vector[1]*trace_vector[3])/16 + (pow(trace_vector[0],2)*pow(trace_vector[2],2))/36 - (trace_vector[5]*pow(trace_vector[0],2))/12 + (trace_vector[0]*pow(trace_vector[1],2)*trace_vector[2])/24 - (trace_vector[4]*trace_vector[0]*trace_vector[1])/10 - (trace_vector[0]*trace_vector[2]*trace_vector[3])/12 + (trace_vector[6]*trace_vector[0])/7 + pow(trace_vector[1],4)/384 - (pow(trace_vector[1],2)*trace_vector[3])/32 - (trace_vector[1]*pow(trace_vector[2],2))/36 + (trace_vector[5]*trace_vector[1])/12 + (trace_vector[4]*trace_vector[2])/15 + pow(trace_vector[3],2)/32 - trace_vector[7]/8; break;
      case 9 : det = (pow(trace_vector[0],3)*pow(trace_vector[1],3))/288 - trace_vector[8]/9 - (pow(trace_vector[0],3)*pow(trace_vector[2],2))/108 - (pow(trace_vector[0],5)*pow(trace_vector[1],2))/960 + (trace_vector[0]*trace_vector[7])/8 + (trace_vector[1]*trace_vector[6])/14 + (trace_vector[2]*trace_vector[5])/18 + (trace_vector[3]*trace_vector[4])/20 - (trace_vector[0]*pow(trace_vector[1],4))/384 - (trace_vector[0]*pow(trace_vector[3],2))/32 + (pow(trace_vector[1],3)*trace_vector[2])/144 - (pow(trace_vector[1],2)*trace_vector[4])/40 - (pow(trace_vector[0],2)*trace_vector[6])/14 + (pow(trace_vector[0],3)*trace_vector[5])/36 - (pow(trace_vector[0],4)*trace_vector[4])/120 + (pow(trace_vector[0],5)*trace_vector[3])/480 - (pow(trace_vector[0],6)*trace_vector[2])/2160 + (pow(trace_vector[0],7)*trace_vector[1])/10080 - pow(trace_vector[2],3)/162 - pow(trace_vector[0],9)/362880 + (trace_vector[0]*trace_vector[1]*pow(trace_vector[2],2))/36 + (trace_vector[0]*pow(trace_vector[1],2)*trace_vector[3])/32 + (pow(trace_vector[0],2)*trace_vector[1]*trace_vector[4])/20 + (pow(trace_vector[0],2)*trace_vector[2]*trace_vector[3])/24 - (pow(trace_vector[0],3)*trace_vector[1]*trace_vector[3])/48 + (pow(trace_vector[0],4)*trace_vector[1]*trace_vector[2])/144 - (pow(trace_vector[0],2)*pow(trace_vector[1],2)*trace_vector[2])/48 - (trace_vector[0]*trace_vector[1]*trace_vector[5])/12 - (trace_vector[0]*trace_vector[2]*trace_vector[4])/15 - (trace_vector[1]*trace_vector[2]*trace_vector[3])/24; break;
      case 10 : det = pow(trace_vector[0],10)/3628800 - (pow(trace_vector[0],8)*trace_vector[1])/80640 + (pow(trace_vector[0],7)*trace_vector[2])/15120 + (pow(trace_vector[0],6)*pow(trace_vector[1],2))/5760 - (pow(trace_vector[0],6)*trace_vector[3])/2880 - (pow(trace_vector[0],5)*trace_vector[1]*trace_vector[2])/720 + (pow(trace_vector[0],5)*trace_vector[4])/600 - (pow(trace_vector[0],4)*pow(trace_vector[1],3))/1152 + (pow(trace_vector[0],4)*trace_vector[1]*trace_vector[3])/192 + (pow(trace_vector[0],4)*pow(trace_vector[2],2))/432 - (trace_vector[5]*pow(trace_vector[0],4))/144 + (pow(trace_vector[0],3)*pow(trace_vector[1],2)*trace_vector[2])/144 - (pow(trace_vector[0],3)*trace_vector[1]*trace_vector[4])/60 - (pow(trace_vector[0],3)*trace_vector[2]*trace_vector[3])/72 + (trace_vector[6]*pow(trace_vector[0],3))/42 + (pow(trace_vector[0],2)*pow(trace_vector[1],4))/768 - (pow(trace_vector[0],2)*pow(trace_vector[1],2)*trace_vector[3])/64 - (pow(trace_vector[0],2)*trace_vector[1]*pow(trace_vector[2],2))/72 + (trace_vector[5]*pow(trace_vector[0],2)*trace_vector[1])/24 + (pow(trace_vector[0],2)*trace_vector[2]*trace_vector[4])/30 + (pow(trace_vector[0],2)*pow(trace_vector[3],2))/64 - (trace_vector[7]*pow(trace_vector[0],2))/16 - (trace_vector[0]*pow(trace_vector[1],3)*trace_vector[2])/144 + (trace_vector[0]*pow(trace_vector[1],2)*trace_vector[4])/40 + (trace_vector[0]*trace_vector[1]*trace_vector[2]*trace_vector[3])/24 - (trace_vector[6]*trace_vector[0]*trace_vector[1])/14 + (trace_vector[0]*pow(trace_vector[2],3))/162 - (trace_vector[5]*trace_vector[0]*trace_vector[2])/18 - (trace_vector[0]*trace_vector[3]*trace_vector[4])/20 + (trace_vector[8]*trace_vector[0])/9 - pow(trace_vector[1],5)/3840 + (pow(trace_vector[1],3)*trace_vector[3])/192 + (pow(trace_vector[1],2)*pow(trace_vector[2],2))/144 - (trace_vector[5]*pow(trace_vector[1],2))/48 - (trace_vector[1]*trace_vector[2]*trace_vector[4])/30 - (trace_vector[1]*pow(trace_vector[3],2))/64 + (trace_vector[7]*trace_vector[1])/16 - (pow(trace_vector[2],2)*trace_vector[3])/72 + (trace_vector[6]*trace_vector[2])/21 + (trace_vector[5]*trace_vector[3])/24 + pow(trace_vector[4],2)/50 - trace_vector[9]/10; break;
    }

    if (std::abs(det) < epsilon){
      det = 0;
    }

    return nearbyint(det);
  }
};
