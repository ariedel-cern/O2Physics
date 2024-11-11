// Copyright 2019-2024 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.femtobaseder
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef PWGCF_FEMTOUNITED_DATAMODEL_FEMTOCASCADESDERIVED_H_
#define PWGCF_FEMTOUNITED_DATAMODEL_FEMTOCASCADESDERIVED_H_

#include "Framework/ASoA.h"
#include "Framework/Expressions.h"
#include "PWGCF/FemtoUnited/Core/DataTypes.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoBaseDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoTracksDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoVzerosDerived.h"

namespace o2::aod
{
namespace femtocascades
{
// columns for Vzero
DECLARE_SOA_COLUMN(CascadeMass, cascadeMass, float);                           //! Mass of Vzero
DECLARE_SOA_COLUMN(CascaseMask, cascadeMask, femtodatatypes::CascadeMaskType); //! Bitmask for Vzero selections

} // namespace femtocascades

// table for basic cascade information
DECLARE_SOA_TABLE_VERSIONED(FUCascades_001, "AOD", "FUCASCADESS", 1,
                            o2::soa::Index<>,
                            femtobase::CollisionId,
                            femtobase::Pt,
                            femtobase::Eta,
                            femtobase::Phi,
                            femtocascades::CascadeMass,
                            femtobase::Theta<femtobase::Eta>,
                            femtobase::Px<femtobase::Pt, femtobase::Eta>,
                            femtobase::Py<femtobase::Pt, femtobase::Eta>,
                            femtobase::Pz<femtobase::Pt, femtobase::Eta>,
                            femtobase::P<femtobase::Pt, femtobase::Eta>);
using FUCascades = FUCascades_001;

} // namespace o2::aod

#endif // PWGCF_FEMTOUNITED_DATAMODEL_FEMTOCASCADESDERIVED_H_
