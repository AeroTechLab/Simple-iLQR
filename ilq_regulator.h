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


/// @file ilq_regulator.h
/// @brief Iterative linear-quadratic regulator implementation
///
/// Calculate optimal state feedback with "cheap control" iterative linear quadratic regulator (iLQR)

#ifndef ILQR_H
#define ILQR_H

#include <stddef.h>
#include <stdbool.h>

typedef struct _ILQRegulatorData ILQRegulatorData;      ///< Single iLQ regulator data structure
typedef ILQRegulatorData* ILQRegulator;                 ///< Opaque reference to iLQ regulator data structure

/// @brief Creates and initializes internal matrices of iLQ regulator data structure 
/// @param[in] statesNumber size (in elements) of the plant internal state (error) vector    
/// @param[in] inputsNumber size (in elements) of the plant inputs vector 
/// @param[in] costRatio state error/input weight ratio desired for regulator cost function
/// @return reference to created linear system
ILQRegulator ILQR_Create( size_t statesNumber, size_t inputsNumber, double costRatio );

/// @brief Defines correlation between two state variables in plant model                           
/// @param[in] regulator reference to regulator
/// @param[in] newStateIndex index (in state vector) of output state variable  
/// @param[in] oldStateIndex index (in state vector) of input state variable                                       
/// @param[in] ratio output/input ratio desired on calculation
void ILQR_SetTransitionFactor( ILQRegulator regulator, size_t newStateIndex, size_t oldStateIndex, double ratio );

/// @brief Defines correlation between input and state variables in plant model         
/// @param[in] regulator reference to regulator
/// @param[in] stateIndex index of the correspondent state variable in the internal state vector
/// @param[in] inputIndex index of the input variable in the internal input vector
/// @param[in] ratio output/input ratio desired on calculation
void ILQR_SetInputFactor( ILQRegulator regulator, size_t stateIndex, size_t inputIndex, double ratio );

/// @brief Calculates control input from state feedback with optimal gain                      
/// @param[in] regulator reference to regulator
/// @param[in] statesList array of state (error) values             
/// @param[out] feedbacksList pointer to array where control feedback will be copied
/// @return @return true on successful calculation, false otherwise
bool ILQR_CalculateFeedback( ILQRegulator regulator, double* statesList, double* feedbacksList );

/// @brief Deallocates memory for internal regulator data
/// @param[in] regulator reference to regulator
void ILQR_Delete( ILQRegulator regulator );

#endif // ILQR_H
