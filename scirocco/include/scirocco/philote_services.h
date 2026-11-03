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

    Philote-Cpp RPC Service Lifetime
*/
#pragma once

#include <explicit.h>
#include <grpcpp/grpcpp.h>

#include <memory>
#include <stdexcept>

namespace scirocco {

/// Own RPC services separately from the discipline to avoid the shared-pointer
/// cycle in Philote-Cpp's ExplicitDiscipline::RegisterServices. Keep this
/// object alive until after the gRPC server has shut down and been destroyed.
class PhiloteServices {
   public:
    PhiloteServices(
        grpc::ServerBuilder& builder,
        const std::shared_ptr<philote::ExplicitDiscipline>& discipline) {
        if (!discipline) throw std::invalid_argument("Discipline is required");
        discipline_server_.LinkPointers(discipline);
        explicit_server_.LinkPointers(discipline);
        discipline->Initialize();
        builder.RegisterService(&discipline_server_);
        builder.RegisterService(&explicit_server_);
    }

    PhiloteServices(const PhiloteServices&) = delete;
    PhiloteServices& operator=(const PhiloteServices&) = delete;

   private:
    philote::DisciplineServer discipline_server_;
    philote::ExplicitServer explicit_server_;
};

}  // namespace scirocco
