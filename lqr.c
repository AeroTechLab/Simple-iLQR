#include "lqr.h"

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
  
  size_t variablesNumber = ( statesNumber > inputsNumber ) ? statesNumber : inputsNumber;
  newRegulator->aux = Mat_CreateSquare( variablesNumber, MATRIX_ZERO );
  
  return newRegulator;
}

void ILQR_SetPredictionFactor( ILQRegulator regulator, size_t outputIndex, size_t inputIndex, double ratio )
{
  if( regulator == NULL ) return;
  
  Mat_SetElement( regulator->dynamicModel, outputIndex, inputIndex, ratio );
}

void ILQR_SetInputFactor( ILQRegulator regulator, size_t outputIndex, size_t inputIndex, double ratio )
{
  if( regulator == NULL ) return;
  
  Mat_SetElement( regulator->inputModel, outputIndex, inputIndex, ratio );
}

bool ILQR_CalculateFeedback( ILQRegulator regulator, double* statesList, double* inputsList )
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
  if( Mat_Inverse( regulator->gain, regulator->gain ) ) return false;
  Mat_Dot( regulator->gain, MATRIX_KEEP, regulator->inputModel, MATRIX_TRANSPOSE, regulator->gain );
  Mat_Dot( regulator->gain, MATRIX_KEEP, regulator->cost2Go, MATRIX_KEEP, regulator->gain );
  Mat_Dot( regulator->gain, MATRIX_KEEP, regulator->dynamicModel, MATRIX_KEEP, regulator->gain );
  // u = K * x
  Mat_Resize( regulator->aux, Mat_GetWidth( regulator->gain ), 1 );
  Mat_SetData( regulator->aux, statesList );
  Mat_Dot( regulator->gain, MATRIX_KEEP, regulator->aux, MATRIX_KEEP, regulator->aux );
  
  Mat_GetData( regulator->aux, inputsList );
  
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
  
  return;
}
