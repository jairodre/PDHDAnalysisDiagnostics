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
AxisAlignedVolume Shrink(AxisAlignedVolume const &volume, double marginCm) {
  return {volume.minX + marginCm, volume.maxX - marginCm,
          volume.minY + marginCm, volume.maxY - marginCm,
          volume.minZ + marginCm, volume.maxZ - marginCm};
}
} // namespace
bool GeometryContainmentAlg::Contains(AxisAlignedVolume const &volume,
                                      Point3D const &point,
                                      double marginCm) noexcept {
  bool const marginValid = std::isfinite(marginCm) && marginCm >= 0.;
  if (!volume.IsValid() || !marginValid) {
    return false;
  }

  auto const containedVolume = Shrink(volume, marginCm);
  if (!containedVolume.IsValid()) {
    return false;
  }

  return std::isfinite(point.x) && point.x >= containedVolume.minX &&
         point.x <= containedVolume.maxX && std::isfinite(point.y) &&
         point.y >= containedVolume.minY && point.y <= containedVolume.maxY &&
         std::isfinite(point.z) && point.z >= containedVolume.minZ &&
         point.z <= containedVolume.maxZ;
}
ContainmentResult GeometryContainmentAlg::Evaluate(
    AxisAlignedVolume const &volume, Point3D const &start, Point3D const &end,
    std::vector<Point3D> const &trajectory, double marginCm) noexcept {
  ContainmentResult result;
  bool const marginValid = std::isfinite(marginCm) && marginCm >= 0.;
  if (!volume.IsValid() || !marginValid) {
    return result;
  }

  auto const containedVolume = Shrink(volume, marginCm);
  if (!containedVolume.IsValid()) {
    return result;
  }

  result.volumeValid = true;
  result.startContained = Contains(volume, start, marginCm);
  result.endContained = Contains(volume, end, marginCm);
  result.bothEndpointsContained = result.startContained && result.endContained;
  result.hasContainedSample = std::any_of(
      trajectory.begin(), trajectory.end(),
      [&](Point3D const &point) { return Contains(volume, point, marginCm); });
  result.allSampledPointsContained =
      !trajectory.empty() &&
      std::all_of(trajectory.begin(), trajectory.end(),
                  [&](Point3D const &point) {
                    return Contains(volume, point, marginCm);
                  });
  return result;
}
std::optional<std::pair<Point3D, Point3D>>
GeometryContainmentAlg::ClipSegment(AxisAlignedVolume const &volume,
                                    Point3D const &start, Point3D const &end,
                                    double marginCm) noexcept {
  bool const marginValid = std::isfinite(marginCm) && marginCm >= 0.;
  if (!volume.IsValid() || !marginValid) {
    return std::nullopt;
  }

  auto const clippedVolume = Shrink(volume, marginCm);
  if (!clippedVolume.IsValid()) {
    return std::nullopt;
  }

  double tMin = 0., tMax = 1.;
  auto updateClipInterval = [&](double origin, double direction, double minimum,
                                double maximum) {
    if (!std::isfinite(origin) || !std::isfinite(direction)) {
      return false;
    }
    if (direction == 0.) {
      return origin >= minimum && origin <= maximum;
    }

    double entry = (minimum - origin) / direction;
    double exit = (maximum - origin) / direction;
    if (entry > exit) {
      std::swap(entry, exit);
    }
    tMin = std::max(tMin, entry);
    tMax = std::min(tMax, exit);
    return tMin <= tMax;
  };
  Point3D const segmentDirection{end.x - start.x, end.y - start.y,
                                 end.z - start.z};
  if (!updateClipInterval(start.x, segmentDirection.x, clippedVolume.minX,
                          clippedVolume.maxX) ||
      !updateClipInterval(start.y, segmentDirection.y, clippedVolume.minY,
                          clippedVolume.maxY) ||
      !updateClipInterval(start.z, segmentDirection.z, clippedVolume.minZ,
                          clippedVolume.maxZ)) {
    return std::nullopt;
  }

  auto const pointAt = [&](double parameter) {
    return Point3D{start.x + parameter * segmentDirection.x,
                   start.y + parameter * segmentDirection.y,
                   start.z + parameter * segmentDirection.z};
  };
  return std::make_pair(pointAt(tMin), pointAt(tMax));
}
} // namespace pdhd::diagnostics
