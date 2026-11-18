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
/* Tests for fauXfoil, the XFOIL-style command-line front end. Expected strings
 * are copied from XFOIL 6.99 sessions so the layout stays identical. */

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "format.h"
#include "session.h"

namespace {

namespace fs = std::filesystem;

/// Run a session on a script inside a scratch directory and return stdout.
class FauxfoilTest : public ::testing::Test {
   protected:
    void SetUp() override {
        dir_ =
            fs::temp_directory_path() /
            ("fauxfoil_test_" + std::string(::testing::UnitTest::GetInstance()
                                                ->current_test_info()
                                                ->name()));
        fs::remove_all(dir_);
        fs::create_directories(dir_);
        previous_ = fs::current_path();
        fs::current_path(dir_);
    }
    void TearDown() override {
        fs::current_path(previous_);
        fs::remove_all(dir_);
    }

    static std::string Run(const std::string& script,
                           const std::vector<std::string>& args = {}) {
        std::istringstream in(script);
        std::ostringstream out;
        fauxfoil::Session session(in, out, "0.2.0");
        EXPECT_EQ(session.Run(args), 0);
        return out.str();
    }

    static std::vector<std::string> Lines(const std::string& path) {
        std::ifstream stream(path);
        std::vector<std::string> lines;
        std::string line;
        while (std::getline(stream, line)) lines.push_back(line);
        return lines;
    }

    static bool Has(const std::string& text, const std::string& part) {
        return text.find(part) != std::string::npos;
    }

    fs::path dir_;
    fs::path previous_;
};

TEST(FauxfoilFormat, MatchesGfortranEditDescriptors) {
    EXPECT_EQ(fauxfoil::E(0.1925, 10, 4), "0.1925E+00");
    EXPECT_EQ(fauxfoil::E(-0.7031, 10, 4), "-.7031E+00");
    EXPECT_EQ(fauxfoil::E(1.0e6, 12, 4), "  0.1000E+07");
    EXPECT_EQ(fauxfoil::E(0.0, 10, 4), "0.0000E+00");
    EXPECT_EQ(fauxfoil::G(1.0, 15, 7), "   1.000000    ");
    EXPECT_EQ(fauxfoil::G(0.9916796, 15, 7), "  0.9916796    ");
    EXPECT_EQ(fauxfoil::G(0.00126, 15, 7), "  0.1260000E-02");
    EXPECT_EQ(fauxfoil::G(0.0, 12, 4), "   0.000    ");
    EXPECT_EQ(fauxfoil::F(-0.05202, 9, 4), "  -0.0520");
    EXPECT_EQ(fauxfoil::F(1234.5, 5, 1), "*****");
    EXPECT_EQ(fauxfoil::I(30, 12), "          30");
    EXPECT_EQ(fauxfoil::A("pol.txt", 10), "pol.txt   ");
}

TEST_F(FauxfoilTest, TopLevelLooksLikeXfoil) {
    const std::string out = Run("xyz\nquit\n");
    EXPECT_TRUE(Has(out, "  Caveat computor\n"));
    EXPECT_TRUE(Has(out,
                    "\n   QUIT    Exit program\n\n  .OPER    Direct "
                    "operating point(s)\n"));
    EXPECT_TRUE(Has(out,
                    "   XYCM rr Change CM reference location, currently "
                    " 0.25000 0.00000\n"));
    EXPECT_TRUE(Has(out,
                    "   PANE    Set current-airfoil panel nodes ( 160 ) "
                    "based on curvature\n"));
    EXPECT_TRUE(Has(out, "   U       Unzoom  | \n"));
    EXPECT_TRUE(Has(out, "\n fauXFOIL   c>  "));
    EXPECT_TRUE(
        Has(out, " XYZ  command not recognized.  Type a \"?\" for list\n"));
}

TEST_F(FauxfoilTest, ScriptsThatDisableGraphicsStillRun) {
    /* "PLOP / G / <cr>" opens nearly every XFOIL batch script */
    const std::string out = Run("plop\ng\n\nnaca 0012\nquit\n");
    EXPECT_TRUE(Has(out, "  G raphics-enable flag:        T\n"));
    EXPECT_TRUE(Has(out, "      Option, Value   (or <Return>)    c>  "));
    EXPECT_TRUE(Has(out, " Buffer airfoil set using"));
}

TEST_F(FauxfoilTest, NacaReportsGeometryAndPaneling) {
    const std::string out = Run("naca 0012\nquit\n");
    /* XFOIL 6.99 gives 0.120032 from its own 245-point buffer */
    EXPECT_TRUE(Has(out, " Max thickness =     0.1200"));
    EXPECT_TRUE(Has(out, " Max camber    =    -0.000000") ||
                Has(out, " Max camber    =     0.000000"));
    EXPECT_TRUE(Has(out, "\n Blunt trailing edge.  Gap =  0.00252\n"));
    EXPECT_TRUE(Has(out,
                    "\n Paneling parameters used...\n"
                    "   Number of panel nodes       160\n"));
}

TEST_F(FauxfoilTest, SaveWritesXfoilOrderAndLoadReadsItBack) {
    const std::string out =
        Run("naca 0012\nsave foil.dat\nload foil.dat\nquit\n");
    const std::vector<std::string> lines = Lines("foil.dat");
    ASSERT_EQ(lines.size(), 161u);
    EXPECT_EQ(lines[0], "NACA 0012");
    /* Counterclockwise from the upper trailing edge, in (1X,2G15.7) */
    EXPECT_EQ(lines[1], "    1.000000      0.1260000E-02");
    EXPECT_EQ(lines[160], "    1.000000     -0.1260000E-02");
    EXPECT_TRUE(Has(out, " Labeled airfoil file.  Name:  NACA 0012"));
    EXPECT_TRUE(Has(out,
                    " Number of input coordinate points: 160\n"
                    " Counterclockwise ordering\n"));
    EXPECT_TRUE(Has(out,
                    "  LE  x,y  =  -0.00000   0.00000  |   Chord =   "
                    "1.00000\n") ||
                Has(out,
                    "  LE  x,y  =   0.00000   0.00000  |   Chord =   "
                    "1.00000\n"));
    EXPECT_TRUE(Has(out,
                    " Current airfoil nodes set from buffer airfoil "
                    "nodes ( 160 )\n"));
}

TEST_F(FauxfoilTest, SaveAsksBeforeOverwriting) {
    const std::string out =
        Run("naca 0012\nsave foil.dat\nsave foil.dat\nn\nquit\n");
    EXPECT_TRUE(Has(out,
                    "\n Output file exists.  Overwrite?  Y\n"
                    " Current airfoil not saved.\n"));
}

TEST_F(FauxfoilTest, MissingFileReportsLoadNotCompleted) {
    const std::string out = Run("load nofile.dat\nquit\n");
    EXPECT_TRUE(Has(out, " File OPEN error.  Nonexistent file:  nofile.dat"));
    EXPECT_TRUE(Has(out, " *** LOAD NOT COMPLETED ***\n"));
}

TEST_F(FauxfoilTest, InviscidLiftTarget) {
    const std::string out = Run("naca 2412\noper\ncl 0.5\n\nquit\n");
    EXPECT_TRUE(Has(out, "\n.OPERi   c>  "));
    EXPECT_TRUE(Has(out, " Calculating unit vorticity distributions ...\n"));
    EXPECT_TRUE(Has(out, "      CL =  0.5000\n"));
}

TEST_F(FauxfoilTest, ViscousPolarMatchesXfoilLayout) {
    const std::string out = Run(
        "naca 0012\noper\nvisc 1e6\npacc\npolar.txt\n\naseq 0 2 2\n\nquit\n");
    EXPECT_TRUE(Has(out, "\n M  =    0.0000\n Re =  0.1000E+07\n"));
    EXPECT_TRUE(Has(out, "\n.OPERva   c>  "));
    EXPECT_TRUE(Has(out, " Solving BL system ...\n"));
    EXPECT_TRUE(Has(out, " Side 1  free  transition at x/c = "));
    /* Iterations stream as the solver runs: one block per Newton step */
    EXPECT_TRUE(Has(out, "\n   1   rms: "));
    EXPECT_TRUE(Has(out, "\n   2   rms: "));
    EXPECT_TRUE(Has(out, "       a =  0.000      CL = "));
    EXPECT_TRUE(Has(out,
                    " Point added to stored polar  1\n"
                    " Point written to save file  polar.txt"));
    EXPECT_FALSE(Has(out, "Convergence failed"));

    const std::vector<std::string> lines = Lines("polar.txt");
    ASSERT_EQ(lines.size(), 14u);
    EXPECT_EQ(lines[1], "       fauXfoil      Version 0.2.0");
    EXPECT_EQ(lines[3].size(), 71u);
    EXPECT_EQ(lines[5],
              " 1 1 Reynolds number fixed          Mach number fixed         ");
    EXPECT_EQ(lines[7], " xtrf =   1.000 (top)        1.000 (bottom)  ");
    EXPECT_EQ(lines[8],
              " Mach =   0.000     Re =     1.000 e 6     Ncrit =   9.000  "
              "9.000");
    EXPECT_EQ(lines[10],
              "   alpha    CL        CD       CDp       CM     Top_Xtr  "
              "Bot_Xtr  Top_Itr  Bot_Itr");
    for (std::size_t k = 10; k < lines.size(); ++k)
        EXPECT_EQ(lines[k].size(), 82u) << lines[k];

    /* Symmetric section at zero incidence */
    std::istringstream row(lines[12]);
    double alpha = 1.0;
    double cl = 1.0;
    double cd = 0.0;
    row >> alpha >> cl >> cd;
    EXPECT_EQ(alpha, 0.0);
    EXPECT_NEAR(cl, 0.0, 1e-4);
    EXPECT_GT(cd, 0.004);
    EXPECT_LT(cd, 0.007);
}

TEST_F(FauxfoilTest, ParametersLockedWhileAccumulating) {
    const std::string out =
        Run("naca 0012\noper\npacc\n\n\nmach 0.3\n\nquit\n");
    EXPECT_TRUE(Has(out,
                    "\n * Polar is being accumulated.\n"
                    " * Cannot change its parameters in midstream.\n"));
}

}  // namespace
