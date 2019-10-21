/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
//  Copyright (c) 2019 Leonardo Consoni <leonardojc@protonmail.com>            //
//                                                                             //
//  This file is part of Simple-iLQR.                                          //
//                                                                             //
//  Simple-iLQR is free software: you can redistribute it                      //
//  and/or modify it under the terms of the GNU Lesser General Public License  //
//  as published by the Free Software Foundation, either version 3 of the      //
//  License, or (at your option) any later version.                            //
//                                                                             //
//  Simple-iLQR is distributed in the hope that it will                        //
//  be useful, but WITHOUT ANY WARRANTY; without even the implied warranty     //
//  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            //
//  GNU Lesser General Public License for more details.                        //
//                                                                             //
//  You should have received a copy of the GNU Lesser General Public License   //
//  along with Simple-iLQR. If not, see <http://www.gnu.org/licenses/>.        //
//                                                                             //
/////////////////////////////////////////////////////////////////////////////////


#include "ilq_regulator.h"

#include "matrix/matrix.h"

#include <stdlib.h>

struct _ILQRegulatorData
{
  Matrix dynamicModel;
  Matrix inputModel;
  Matrix cost2Go;
  Matrix costWeight;
  double costRatio;
  Matrix gain;
  Matrix aux;
  Matrix state;
};

ILQRegulator ILQR_Create( size_t statesNumber, size_t inputsNumber, double costRatio )
{
  ILQRegulator newRegulator = (ILQRegulator) malloc( sizeof(ILQRegulatorData) );
  
  newRegulator->dynamicModel = Mat_CreateSquare( statesNumber, MATRIX_IDENTITY );
  newRegulator->inputModel = Mat_Create( NULL, statesNumber, inputsNumber );
  
  newRegulator->cost2Go = Mat_CreateSquare( statesNumber, MATRIX_IDENTITY );
  newRegulator->gain = Mat_Create( NULL, inputsNumber, statesNumber );
  newRegulator->costWeight = Mat_CreateSquare( statesNumber, MATRIX_IDENTITY );
  newRegulator->costRatio = costRatio;
  
  newRegulator->aux = Mat_CreateSquare( statesNumber, MATRIX_ZERO );
  
  newRegulator->state = Mat_Create( NULL, statesNumber, 1 );
  
  return newRegulator;
}

void ILQR_SetTransitionFactor( ILQRegulator regulator, size_t oldStateIndex, size_t newStateIndex, double ratio )
{
  if( regulator == NULL ) return;
  
  Mat_SetElement( regulator->dynamicModel, newStateIndex, oldStateIndex, ratio );
}

void ILQR_SetInputFactor( ILQRegulator regulator, size_t inputIndex, size_t stateIndex, double ratio )
{
  if( regulator == NULL ) return;
  
  Mat_SetElement( regulator->inputModel, stateIndex, inputIndex, ratio );
}

bool ILQR_CalculateFeedback( ILQRegulator regulator, double* statesList, double* feedbacksList )
{
  if( regulator == NULL ) return false;
  // At * X * B
  Mat_Dot( regulator->dynamicModel, MATRIX_TRANSPOSE, regulator->cost2Go, MATRIX_KEEP, regulator->aux );
  Mat_Dot( regulator->aux, MATRIX_KEEP, regulator->inputModel, MATRIX_KEEP, regulator->aux );
  Mat_Dot( regulator->aux, MATRIX_KEEP, regulator->gain, MATRIX_KEEP, regulator->aux );
  // Q + At * X * A
  Mat_Dot( regulator->dynamicModel, MATRIX_TRANSPOSE, regulator->cost2Go, MATRIX_KEEP, regulator->cost2Go );
  Mat_Dot( regulator->cost2Go, MATRIX_KEEP, regulator->dynamicModel, MATRIX_KEEP, regulator->cost2Go );
  Mat_Sum( regulator->cost2Go, 1.0, regulator->costWeight, 1.0, regulator->cost2Go );
  // Q + At * X * A - At * X * B * K
  Mat_Sum( regulator->cost2Go, 1.0, regulator->aux, -1.0, regulator->cost2Go );
  // K = ( R + Bt * X * B )^-1 * Bt * X * A
  Mat_Dot( regulator->inputModel, MATRIX_TRANSPOSE, regulator->cost2Go, MATRIX_KEEP, regulator->gain );
  Mat_Dot( regulator->gain, MATRIX_KEEP, regulator->inputModel, MATRIX_KEEP, regulator->gain );
  Mat_Sum( regulator->gain, 1.0, regulator->costWeight, regulator->costRatio, regulator->gain );
  if( Mat_Inverse( regulator->gain, regulator->gain ) == NULL ) return false;
  Mat_Dot( regulator->gain, MATRIX_KEEP, regulator->inputModel, MATRIX_TRANSPOSE, regulator->gain );
  Mat_Dot( regulator->gain, MATRIX_KEEP, regulator->cost2Go, MATRIX_KEEP, regulator->gain );
  Mat_Dot( regulator->gain, MATRIX_KEEP, regulator->dynamicModel, MATRIX_KEEP, regulator->gain );
  // u = K * x
  Mat_SetData( regulator->state, statesList );
  Mat_Dot( regulator->gain, MATRIX_KEEP, regulator->state, MATRIX_KEEP, regulator->aux );
  
  Mat_GetData( regulator->aux, feedbacksList );
  
  return true;
}

void ILQR_Delete( ILQRegulator regulator )
{
  if( regulator == NULL ) return;
  
  Mat_Discard( regulator->dynamicModel );
  Mat_Discard( regulator->inputModel );
  Mat_Discard( regulator->cost2Go );
  Mat_Discard( regulator->gain );
  Mat_Discard( regulator->aux );
  Mat_Discard( regulator->state );
  
  return;
}
