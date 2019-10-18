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

#ifndef ILQR_H
#define ILQR_H

#include <stddef.h>
#include <stdbool.h>

typedef struct _ILQRegulatorData ILQRegulatorData;
typedef ILQRegulatorData* ILQRegulator;

ILQRegulator ILQR_Create( size_t statesNumber, size_t inputsNumber, double costRatio );
void ILQR_SetPredictionFactor( ILQRegulator regulator, size_t outputIndex, size_t inputIndex, double ratio );
void ILQR_SetInputFactor( ILQRegulator regulator, size_t outputIndex, size_t inputIndex, double ratio );
bool ILQR_CalculateFeedback( ILQRegulator regulator, double* statesList, double* inputsList );
void ILQR_Delete( ILQRegulator regulator );

#endif // ILQR_H
