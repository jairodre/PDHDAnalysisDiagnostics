/**
 * @file GeometryContainmentAlg.cxx
 * @brief Implements inclusive point containment and slab-based segment
 * clipping.
 *
 * A positive margin shrinks every face. Invalid volumes, negative margins, or
 * margins that remove the volume yield invalid results. Boundaries are
 * inclusive.
 */
#include "protoduneana/PDHDAnalysisDiagnostics/Alg/GeometryContainmentAlg.h"
#include <algorithm>
#include <cmath>

namespace pdhd::diagnostics {
bool AxisAlignedVolume::IsValid() const noexcept {
  return std::isfinite(minX) && std::isfinite(maxX) && minX <= maxX &&
         std::isfinite(minY) && std::isfinite(maxY) && minY <= maxY &&
         std::isfinite(minZ) && std::isfinite(maxZ) && minZ <= maxZ;
}
namespace {
AxisAlignedVolume Shrink(AxisAlignedVolume const &v, double m) {
  return {v.minX + m, v.maxX - m, v.minY + m,
          v.maxY - m, v.minZ + m, v.maxZ - m};
}
} // namespace
bool GeometryContainmentAlg::Contains(AxisAlignedVolume const &volume,
                                      Point3D const &point,
                                      double marginCm) noexcept {
  if (!volume.IsValid() || !std::isfinite(marginCm) || marginCm < 0.)
    return false;
  auto const v = Shrink(volume, marginCm);
  if (!v.IsValid())
    return false;
  return std::isfinite(point.x) && point.x >= v.minX && point.x <= v.maxX &&
         std::isfinite(point.y) && point.y >= v.minY && point.y <= v.maxY &&
         std::isfinite(point.z) && point.z >= v.minZ && point.z <= v.maxZ;
}
ContainmentResult GeometryContainmentAlg::Evaluate(
    AxisAlignedVolume const &volume, Point3D const &start, Point3D const &end,
    std::vector<Point3D> const &trajectory, double marginCm) noexcept {
  ContainmentResult r;
  if (!volume.IsValid() || !std::isfinite(marginCm) || marginCm < 0. ||
      !Shrink(volume, marginCm).IsValid())
    return r;
  r.volumeValid = true;
  r.startContained = Contains(volume, start, marginCm);
  r.endContained = Contains(volume, end, marginCm);
  r.bothEndpointsContained = r.startContained && r.endContained;
  r.hasContainedSample =
      std::any_of(trajectory.begin(), trajectory.end(),
                  [&](auto const &p) { return Contains(volume, p, marginCm); });
  r.allSampledPointsContained =
      !trajectory.empty() &&
      std::all_of(trajectory.begin(), trajectory.end(),
                  [&](auto const &p) { return Contains(volume, p, marginCm); });
  return r;
}
std::optional<std::pair<Point3D, Point3D>>
GeometryContainmentAlg::ClipSegment(AxisAlignedVolume const &volume,
                                    Point3D const &start, Point3D const &end,
                                    double marginCm) noexcept {
  if (!volume.IsValid() || !std::isfinite(marginCm) || marginCm < 0.)
    return std::nullopt;
  auto const v = Shrink(volume, marginCm);
  if (!v.IsValid())
    return std::nullopt;
  double tMin = 0., tMax = 1.;
  auto update = [&](double o, double d, double lo, double hi) {
    if (!std::isfinite(o) || !std::isfinite(d))
      return false;
    if (d == 0.)
      return o >= lo && o <= hi;
    double a = (lo - o) / d, b = (hi - o) / d;
    if (a > b)
      std::swap(a, b);
    tMin = std::max(tMin, a);
    tMax = std::min(tMax, b);
    return tMin <= tMax;
  };
  Point3D const d{end.x - start.x, end.y - start.y, end.z - start.z};
  if (!update(start.x, d.x, v.minX, v.maxX) ||
      !update(start.y, d.y, v.minY, v.maxY) ||
      !update(start.z, d.z, v.minZ, v.maxZ))
    return std::nullopt;
  auto at = [&](double t) {
    return Point3D{start.x + t * d.x, start.y + t * d.y, start.z + t * d.z};
  };
  return std::make_pair(at(tMin), at(tMax));
}
} // namespace pdhd::diagnostics
