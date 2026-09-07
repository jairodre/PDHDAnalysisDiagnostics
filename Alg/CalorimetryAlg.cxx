/**
 * @file CalorimetryAlg.cxx
 * @brief Implements bounds-checked extraction of every calorimetry vector.
 *
 * The largest input-vector size defines row count. Missing elements remain
 * invalid rather than receiving physics-looking sentinel values.
 */
#include "protoduneana/PDHDAnalysisDiagnostics/Alg/CalorimetryAlg.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include <algorithm>
#include <cmath>
#include <utility>
namespace pdhd::diagnostics {
CalorimetryResult CalorimetryAlg::Extract(anab::Calorimetry const &c,
                                          CalorimetryVariant variant,
                                          bool calibrated, std::string producer,
                                          std::string dQdxUnits) {
  CalorimetryResult r;
  auto const &dedx = c.dEdx();
  auto const &dqdx = c.dQdx();
  auto const &rr = c.ResidualRange();
  auto const &pitch = c.TrkPitchVec();
  auto const &xyz = c.XYZ();
  auto const &tp = c.TpIndices();
  auto const &ef = c.Efield();
  auto const &phi = c.Phi();
  std::size_t const n =
      std::max({dedx.size(), dqdx.size(), rr.size(), pitch.size(), xyz.size(),
                tp.size(), ef.size(), phi.size()});
  r.points.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    CalorimetryPoint p;
    p.pointIndex = i;
    if (i < dedx.size() && std::isfinite(dedx[i])) {
      p.dEdxValid = true;
      p.dEdxMeVPerCm = dedx[i];
    }
    if (i < dqdx.size() && std::isfinite(dqdx[i])) {
      p.dQdxValid = true;
      p.dQdx = dqdx[i];
    }
    if (i < rr.size() && std::isfinite(rr[i])) {
      p.residualRangeValid = true;
      p.residualRangeCm = rr[i];
    }
    if (i < pitch.size() && std::isfinite(pitch[i]) && pitch[i] > 0.) {
      p.pitchValid = true;
      p.pitchCm = pitch[i];
      r.summary.sumPitchCm += pitch[i];
      r.summary.maximumPitchCm =
          std::max(r.summary.maximumPitchCm, double(pitch[i]));
    } else
      ++r.summary.invalidPitchCount;
    if (i < xyz.size() && std::isfinite(xyz[i].X()) &&
        std::isfinite(xyz[i].Y()) && std::isfinite(xyz[i].Z())) {
      p.positionValid = true;
      p.xCm = xyz[i].X();
      p.yCm = xyz[i].Y();
      p.zCm = xyz[i].Z();
    }
    if (i < tp.size()) {
      p.trajectoryPointIndexValid = true;
      p.trajectoryPointIndex = tp[i];
    }
    if (i < ef.size() && std::isfinite(ef[i])) {
      p.electricFieldValid = true;
      p.electricField = ef[i];
    }
    if (i < phi.size() && std::isfinite(phi[i])) {
      p.phiValid = true;
      p.phiDegrees = phi[i];
    }
    if (p.dEdxValid && p.dQdxValid && p.residualRangeValid && p.pitchValid &&
        p.positionValid)
      ++r.summary.completePointCount;
    r.points.push_back(p);
  }
  auto const plane = c.PlaneID();
  r.summary.cryostat = plane.Cryostat;
  r.summary.tpc = plane.TPC;
  r.summary.plane = plane.Plane;
  r.summary.kineticEnergyMeV = c.KineticEnergy();
  r.summary.rangeCm = c.Range();
  r.summary.pointCount = n;
  r.summary.variant = variant;
  r.summary.calibrated = calibrated;
  r.summary.producer = std::move(producer);
  r.summary.dQdxUnits = std::move(dQdxUnits);
  r.summary.valid = plane.isValid && n > 0;
  return r;
}
} // namespace pdhd::diagnostics
