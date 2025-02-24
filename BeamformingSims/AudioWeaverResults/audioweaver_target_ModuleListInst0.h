/*******************************************************************************
*
********************************************************************************
*     audioweaver_target_ModuleListInst0.h  Created on 23-Feb-2025 21:18:12
********************************************************************************
*
*     Description:  Module list for the target
*
*     Copyright: (c) 2023 DSP Concepts, Inc. All rights reserved.
*                         3235 Kifer Road
*                         Santa Clara, CA 95054-1527
*
*******************************************************************************/
#include "AWECore.h"

#define TOTALNUMBEROFCLASSES 7

extern const UINT32 awe_modFftClass;
extern const UINT32 awe_modIfftClass;
extern const UINT32 awe_modRepWinOverlapClass;
extern const UINT32 awe_modSbBeamformerV3Class;
extern const UINT32 awe_modSetWirePropertiesClass;
extern const UINT32 awe_modTypeConversionClass;
extern const UINT32 awe_modWindowAliasClass;


#define LISTOFCLASSOBJECTS \
&awe_modFftClass, \
&awe_modIfftClass, \
&awe_modRepWinOverlapClass, \
&awe_modSbBeamformerV3Class, \
&awe_modSetWirePropertiesClass, \
&awe_modTypeConversionClass, \
&awe_modWindowAliasClass
