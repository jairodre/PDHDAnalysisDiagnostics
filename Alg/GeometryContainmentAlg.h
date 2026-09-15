/**
 * @file GeometryContainmentAlg.h
 * @brief Detector-independent axis-aligned containment and clipping utilities.
 *
 * Coordinates and margins are in cm. The caller obtains physical detector or
 * fiducial bounds from geometry/configuration; this class never embeds PDHD
 * dimensions. Endpoint containment, both-endpoint containment,
 * sampled-trajectory containment, and visible segment clipping remain separate
 * results.
 */
#ifndef PROTODUNEANA_PDHDANALYSISDIAGNOSTICS_GEOMETRYCONTAINMENTALG_H
#define PROTODUNEANA_PDHDANALYSISDIAGNOSTICS_GEOMETRYCONTAINMENTALG_H

#include <optional>
#include <utility>
#include <vector>

namespace pdhd::diagnostics {
struct Point3D {
  // Cartesian detector coordinate in cm.
  double x = 0.;
  double y = 0.;
  double z = 0.;
};
struct AxisAlignedVolume {
  // Inclusive physical or fiducial bounds in cm, supplied by the caller.
  double minX = 0.;
  double maxX = 0.;
  double minY = 0.;
  double maxY = 0.;
  double minZ = 0.;
  double maxZ = 0.;
  bool IsValid() const noexcept;
};
struct ContainmentResult {
  // Endpoint and sampled-trajectory answers are intentionally independent.
  bool volumeValid = false;
  bool startContained = false;
  bool endContained = false;
  bool bothEndpointsContained = false;
  bool allSampledPointsContained = false;
  bool hasContainedSample = false;
};
class GeometryContainmentAlg {
public:
  // Tests a finite point against an inclusive volume after shrinking every
  // face by a finite non-negative margin. Invalid inputs return false.
  static bool Contains(AxisAlignedVolume const &, Point3D const &,
                       double marginCm = 0.) noexcept;
  // Evaluates endpoints and only the caller-provided trajectory samples; it
  // does not infer unsampled crossings between those points.
  static ContainmentResult Evaluate(AxisAlignedVolume const &, Point3D const &,
                                    Point3D const &,
                                    std::vector<Point3D> const &,
                                    double marginCm = 0.) noexcept;
  // Clips the finite segment to the same inclusive, shrunken volume; nullopt
  // denotes invalid input or no intersection.
  static std::optional<std::pair<Point3D, Point3D>>
  ClipSegment(AxisAlignedVolume const &, Point3D const &, Point3D const &,
              double marginCm = 0.) noexcept;
};
} // namespace pdhd::diagnostics
#endif
