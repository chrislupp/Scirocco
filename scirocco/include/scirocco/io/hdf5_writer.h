/*
    Scirocco: C++ Airfoil Analysis Library

    Copyright 2024-2026 Christopher A. Lupp

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.

    This work has been cleared for public release, distribution unlimited, case
    number: AFRL-2026-0405.

    The views expressed are those of the authors and do not reflect the
    official guidance or position of the United States Government, the
    Department of Defense or of the United States Air Force.

    Statement from DoD: The Appearance of external hyperlinks does not
    constitute endorsement by the United States Department of Defense (DoD) of
    the linked websites, of the information, products, or services contained
    therein. The DoD does not exercise any editorial, security, or other
    control over the information you may find at these locations.
*/
#ifndef SCIROCCO_IO_HDF5_WRITER_H
#define SCIROCCO_IO_HDF5_WRITER_H

#include <memory>
#include <string>

namespace scirocco {

class Scirocco;
struct FoilState;
struct Sensitivity;

/**
 * @brief Writes Scirocco solver state and sensitivities to an HDF5 file.
 *
 * Each call to WriteIteration appends one iteration group to the file.
 * The file is finalized on Close() or destruction.
 */
class Hdf5Writer {
   public:
    /**
     * @brief Opens (or creates) an HDF5 file for writing.
     * @param filepath Path to the HDF5 output file.
     */
    explicit Hdf5Writer(const std::string& filepath);
    ~Hdf5Writer();

    Hdf5Writer(const Hdf5Writer&) = delete;
    Hdf5Writer& operator=(const Hdf5Writer&) = delete;

    /**
     * @brief Writes the current solver state (and optional sensitivities) as
     *        one iteration group.
     * @param solver     Solver whose state is extracted and written.
     * @param sensitivity Optional pointer to sensitivity data; omitted if null.
     */
    void WriteIteration(const Scirocco& solver,
                        const Sensitivity* sensitivity = nullptr);

    /**
     * @brief Writes a FoilState (and optional sensitivities) as one iteration
     *        group.
     * @param state       Foil state to write.
     * @param sensitivity Optional pointer to sensitivity data; omitted if null.
     */
    void WriteIteration(const FoilState& state,
                        const Sensitivity* sensitivity = nullptr);

    /**
     * @brief Flushes and closes the HDF5 file.
     */
    void Close();

    /**
     * @brief Returns the number of iterations written so far.
     * @return Iteration count.
     */
    [[nodiscard]] int NumIterations() const;

   private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

}  // namespace scirocco

#endif  // SCIROCCO_IO_HDF5_WRITER_H
