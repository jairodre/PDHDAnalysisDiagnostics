/**
 * @file PDHDCalorimetryPointInventory_module.cc
 * @brief Template for one-row-per-calorimetry-point diagnostic output.
 *
 * Inputs: configured track/shower calorimetry associations for SCE and no-SCE
 * reconstruction. Outputs: object/plane/point keys, XYZ in cm, wire/tick, pitch
 * in cm, dQ/dx with native units, dE/dx in MeV/cm, residual range in cm, kinetic
 * energy/range when supplied, local electric field, calibration variant, and
 * validity. Preserve original vector ordering and record size mismatches. Include
 * point-quality flags so extreme pitch or invalid coordinates can be excluded
 * reproducibly without discarding the parent track.
 *
 * Status: template-only contract; intentionally unbuilt, unscheduled, and
 * without an art plugin registration.
 */

// TODO: Implement through CalorimetryAlg after unit and validity conventions.
