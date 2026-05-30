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

    Scirocco GeometryDiscipline Philote Server
*/

#include <grpcpp/grpcpp.h>

#include <iostream>
#include <memory>
#include <string>

#include "scirocco/geometry_discipline.h"

int main(int argc, char* argv[]) {
    std::string address = "localhost:50052";
    if (argc > 1) {
        address = argv[1];
    }

    auto discipline = std::make_shared<scirocco::GeometryDiscipline>();

    grpc::ServerBuilder builder;
    builder.AddListeningPort(address, grpc::InsecureServerCredentials());
    discipline->RegisterServices(builder);

    std::unique_ptr<grpc::Server> server = builder.BuildAndStart();
    std::cout << "Scirocco geometry-discipline server listening on " << address
              << '\n';
    server->Wait();

    return 0;
}
