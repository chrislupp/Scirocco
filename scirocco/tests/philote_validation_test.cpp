#include <grpcpp/grpcpp.h>
#include <gtest/gtest.h>

#include <chrono>
#include <limits>
#include <memory>

#include "scirocco/kulfan_discipline.h"
#include "scirocco/philote_discipline.h"
#include "scirocco/philote_services.h"

namespace {
google::protobuf::Struct NacaOptions() {
    google::protobuf::Struct options;
    (*options.mutable_fields())["naca"].set_string_value("0012");
    (*options.mutable_fields())["npanel"].set_number_value(60);
    return options;
}

class RpcServer {
   public:
    explicit RpcServer(
        const std::shared_ptr<philote::ExplicitDiscipline>& discipline) {
        grpc::ServerBuilder builder;
        int port = 0;
        builder.AddListeningPort("127.0.0.1:0",
                                 grpc::InsecureServerCredentials(), &port);
        services =
            std::make_unique<scirocco::PhiloteServices>(builder, discipline);
        server = builder.BuildAndStart();
        if (!server) throw std::runtime_error("Test server failed to start");
        client.ConnectChannel(
            grpc::CreateChannel("127.0.0.1:" + std::to_string(port),
                                grpc::InsecureChannelCredentials()));
        client.SetRPCTimeout(std::chrono::seconds(10));
    }
    ~RpcServer() { server->Shutdown(); }
    void Options(const google::protobuf::Struct& options) {
        philote::DisciplineOptions request;
        *request.mutable_options() = options;
        client.SendOptions(request);
    }
    philote::Variables Setup() {
        client.Setup();
        client.GetVariableDefinitions();
        philote::Variables inputs;
        for (const auto& meta : client.GetVariableMetaAll())
            if (meta.type() == philote::kInput)
                inputs.emplace(meta.name(), philote::Variable(meta));
        return inputs;
    }
    std::unique_ptr<scirocco::PhiloteServices> services;
    philote::ExplicitClient client;
    std::unique_ptr<grpc::Server> server;
};
}  // namespace

TEST(PhiloteValidation, RejectsInvalidOptionTypesAndCounts) {
    scirocco::SiroccoDiscipline analysis;
    scirocco::KulfanDiscipline geometry;
    for (double invalid :
         {0.0, -1.0, 4.5, 4097.0, 1e30, std::numeric_limits<double>::infinity(),
          std::numeric_limits<double>::quiet_NaN()}) {
        auto options = NacaOptions();
        (*options.mutable_fields())["npanel"].set_number_value(invalid);
        EXPECT_THROW(analysis.SetOptions(options), std::invalid_argument);
        EXPECT_THROW(geometry.SetOptions(options), std::invalid_argument);
    }
    auto options = NacaOptions();
    (*options.mutable_fields())["npanel"].set_string_value("60");
    EXPECT_THROW(analysis.SetOptions(options), std::invalid_argument);
    EXPECT_THROW(geometry.SetOptions(options), std::invalid_argument);
    options = NacaOptions();
    (*options.mutable_fields())["naca"].set_number_value(12);
    EXPECT_THROW(analysis.SetOptions(options), std::invalid_argument);
    options = NacaOptions();
    (*options.mutable_fields())["n_cst_upper"].set_number_value(3.5);
    EXPECT_THROW(geometry.SetOptions(options), std::invalid_argument);
}

TEST(PhiloteRpc, MalformedCoordinatesAndNonfiniteInputsAllowRecovery) {
    RpcServer rpc(std::make_shared<scirocco::SiroccoDiscipline>());
    google::protobuf::Struct options;
    auto& coords =
        *(*options.mutable_fields())["coordinates"].mutable_struct_value();
    // Missing fields, wrong types, then mismatched lengths must return RPC
    // errors.
    EXPECT_THROW(rpc.Options(options), philote::PhiloteServerError);
    (*coords.mutable_fields())["x"].set_number_value(1);
    EXPECT_THROW(rpc.Options(options), philote::PhiloteServerError);
    auto* x = (*coords.mutable_fields())["x"].mutable_list_value();
    for (int i = 0; i < 4; ++i) x->add_values()->set_number_value(i);
    (*coords.mutable_fields())["z"]
        .mutable_list_value()
        ->add_values()
        ->set_number_value(0);
    EXPECT_THROW(rpc.Options(options), philote::PhiloteServerError);
    auto* z = (*coords.mutable_fields())["z"].mutable_list_value();
    for (int i = 0; i < 3; ++i) z->add_values()->set_number_value(0);
    z->mutable_values(1)->set_string_value("bad");
    EXPECT_THROW(rpc.Options(options), philote::PhiloteServerError);
    z->mutable_values(1)->set_number_value(
        std::numeric_limits<double>::infinity());
    EXPECT_THROW(rpc.Options(options), philote::PhiloteServerError);

    ASSERT_NO_THROW(rpc.Options(NacaOptions()));
    auto inputs = rpc.Setup();
    scirocco::Scirocco foil("0012", 60);
    for (int i = 0; i < 61; ++i) {
        inputs.at("x")(i) = foil.GetFoilCoords()(0, i);
        inputs.at("x")(61 + i) = foil.GetFoilCoords()(1, i);
    }
    inputs.at("alpha")(0) = std::numeric_limits<double>::quiet_NaN();
    inputs.at("Re")(0) = 0;
    inputs.at("Ma")(0) = 0;
    EXPECT_THROW(rpc.client.ComputeFunction(inputs),
                 philote::PhiloteServerError);
    inputs.at("alpha")(0) = 4;
    auto outputs = rpc.client.ComputeFunction(inputs);
    EXPECT_GT(outputs.at("cl")(0), 0.4);
    EXPECT_EQ(outputs.at("cd")(0), 0);
}

TEST(PhiloteRpc, KulfanRecoversFromInvalidPanelOption) {
    RpcServer rpc(std::make_shared<scirocco::KulfanDiscipline>());
    google::protobuf::Struct options;
    (*options.mutable_fields())["npanel"].set_number_value(4.5);
    EXPECT_THROW(rpc.Options(options), philote::PhiloteServerError);
    (*options.mutable_fields())["npanel"].set_number_value(60);
    ASSERT_NO_THROW(rpc.Options(options));
    auto inputs = rpc.Setup();
    for (auto& [name, variable] : inputs)
        for (size_t i = 0; i < variable.Size(); ++i) variable(i) = 0.1;
    auto outputs = rpc.client.ComputeFunction(inputs);
    EXPECT_EQ(outputs.at("x").Size(), 122);
    for (size_t i = 0; i < outputs.at("x").Size(); ++i)
        EXPECT_TRUE(std::isfinite(outputs.at("x")(i)));
}

TEST(PhiloteRpc, ShutdownReleasesBothDisciplineTypes) {
    for (bool analysis : {true, false}) {
        std::shared_ptr<philote::ExplicitDiscipline> discipline;
        if (analysis)
            discipline = std::make_shared<scirocco::SiroccoDiscipline>();
        else
            discipline = std::make_shared<scirocco::KulfanDiscipline>();
        std::weak_ptr<philote::ExplicitDiscipline> observer = discipline;
        {
            RpcServer rpc(discipline);
            discipline.reset();
            ASSERT_FALSE(observer.expired());
            rpc.Options(NacaOptions());
            EXPECT_FALSE(rpc.Setup().empty());
        }
        EXPECT_TRUE(observer.expired());
    }
}

TEST(PhiloteRpc, UnstartedServicesReleaseDiscipline) {
    std::weak_ptr<philote::ExplicitDiscipline> observer;
    {
        grpc::ServerBuilder builder;
        auto discipline = std::make_shared<scirocco::SiroccoDiscipline>();
        observer = discipline;
        scirocco::PhiloteServices services(builder, discipline);
        discipline.reset();
        EXPECT_FALSE(observer.expired());
        // Simulate leaving before a server is started (e.g. bind failure).
    }
    EXPECT_TRUE(observer.expired());
}
